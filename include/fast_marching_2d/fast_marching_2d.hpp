#pragma once

// fast_marching_2d: the Fast Marching Method (Sethian 1996) for the
// unit-speed Eikonal equation |grad T| = 1 on a regular 2D grid --
// computes the arrival-time (geometric distance) field from a set of
// seed points, each optionally offset by its own initial arrival time.
//
// Not a general speed-field / time-of-arrival solver: propagation speed
// is fixed at 1 everywhere, so T is exactly the distance from the
// nearest seed (plus that seed's arrival_time offset).

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <vector>

#include "common_geometry/grid.hpp"
#include "common_geometry/types.hpp"

namespace ns_fm2d {

// A source point the arrival-time field propagates from. Need not
// coincide with a grid node -- ComputeArrivalTime initializes its
// containing cell's corners via direct Euclidean distance to `position`,
// not by snapping to the nearest node, so accuracy near the seed isn't
// limited to the grid's own resolution.
struct Seed {
  ns_cg::Vec2d position;
  double arrival_time = 0.0;
};

namespace detail {

enum class NodeState { Far, Trial, Known };

struct HeapEntry {
  double arrival_time;
  std::size_t node_index;
};

// Min-heap on arrival_time: std::priority_queue is a max-heap by
// default, so this comparator is inverted.
struct HeapEntryGreater {
  bool operator()(const HeapEntry& a, const HeapEntry& b) const {
    return a.arrival_time > b.arrival_time;
  }
};

inline std::size_t NodeIndex(std::size_t i, std::size_t j, std::size_t nx) {
  return j * nx + i;
}

// Resolves a world-space point to its containing cell's lower-left node
// index, clamped so points on or outside the grid boundary still
// resolve to the last valid cell.
struct CellLocal {
  std::size_t i, j;
};

inline CellLocal LocateCell(const ns_cg::Vec2d& p, const ns_cg::Vec2d& origin,
                             double dx, double dy, std::size_t nx, std::size_t ny) {
  const double fx = (p.x() - origin.x()) / dx;
  const double fy = (p.y() - origin.y()) / dy;
  const long max_i = static_cast<long>(nx) - 2;
  const long max_j = static_cast<long>(ny) - 2;
  const long i = std::clamp<long>(static_cast<long>(std::floor(fx)), 0, std::max(max_i, 0L));
  const long j = std::clamp<long>(static_cast<long>(std::floor(fy)), 0, std::max(max_j, 0L));
  return CellLocal{static_cast<std::size_t>(i), static_cast<std::size_t>(j)};
}

// Solves the 2D Eikonal update at a node from its known neighbor values
// along x (`tx`, if any) and y (`ty`, if any), with grid spacing
// (dx, dy). At least one of tx/ty must be finite.
inline double EikonalUpdate(double tx, double ty, double dx, double dy) {
  const bool has_x = std::isfinite(tx);
  const bool has_y = std::isfinite(ty);
  if (has_x && !has_y) return tx + dx;
  if (has_y && !has_x) return ty + dy;

  // Both known: solve a*T^2 + b*T + c = 0 for
  // (T-tx)^2/dx^2 + (T-ty)^2/dy^2 = 1.
  const double inv_dx2 = 1.0 / (dx * dx);
  const double inv_dy2 = 1.0 / (dy * dy);
  const double a = inv_dx2 + inv_dy2;
  const double b = -2.0 * (tx * inv_dx2 + ty * inv_dy2);
  const double c = tx * tx * inv_dx2 + ty * ty * inv_dy2 - 1.0;
  const double discriminant = b * b - 4.0 * a * c;

  const double one_d_fallback = std::min(tx + dx, ty + dy);
  if (discriminant < 0.0) return one_d_fallback;

  const double t = (-b + std::sqrt(discriminant)) / (2.0 * a);
  // Causality/upwind check: the 2D solution must be no smaller than
  // either contributing neighbor value, else it isn't physically valid
  // and the 1D update (from whichever neighbor is closer) is used
  // instead.
  if (t < std::max(tx, ty)) return one_d_fallback;
  return t;
}

}  // namespace detail

// Computes the arrival-time field T over an nx x ny grid (spacing dx x
// dy, starting at origin) via the Fast Marching Method under unit
// propagation speed (|grad T| = 1): T at each node is the Euclidean
// distance to the nearest seed, plus that seed's arrival_time.
//
// If `seeds` is empty, or nx/ny < 2, every node is left at +infinity
// (nothing to propagate from).
inline ns_cg::Grid2d<double> ComputeArrivalTime(std::size_t nx, std::size_t ny, double dx,
                                                 double dy, const ns_cg::Vec2d& origin,
                                                 const std::vector<Seed>& seeds) {
  using detail::EikonalUpdate;
  using detail::HeapEntry;
  using detail::HeapEntryGreater;
  using detail::LocateCell;
  using detail::NodeIndex;
  using detail::NodeState;

  constexpr double kInf = std::numeric_limits<double>::infinity();
  const std::size_t num_nodes = nx * ny;

  ns_cg::Grid2d<double> field(nx, ny, dx, dy, origin);
  for (std::size_t j = 0; j < ny; ++j) {
    for (std::size_t i = 0; i < nx; ++i) field.at(i, j) = kInf;
  }
  if (nx < 2 || ny < 2 || seeds.empty()) return field;

  std::vector<NodeState> state(num_nodes, NodeState::Far);
  std::vector<double> distance(num_nodes, kInf);

  std::priority_queue<HeapEntry, std::vector<HeapEntry>, HeapEntryGreater> heap;

  const auto push_trial = [&](std::size_t node, double candidate) {
    if (candidate < distance[node]) {
      distance[node] = candidate;
      if (state[node] == NodeState::Far) state[node] = NodeState::Trial;
      heap.push({candidate, node});
    }
  };

  // --- Initialization: each seed's containing cell's 4 corners get a
  // sub-grid-accurate initial distance (direct Euclidean distance to
  // the seed itself, not the nearest node). ---
  for (const Seed& seed : seeds) {
    const detail::CellLocal cell = LocateCell(seed.position, origin, dx, dy, nx, ny);
    const std::size_t corner_i[4] = {cell.i, cell.i + 1, cell.i + 1, cell.i};
    const std::size_t corner_j[4] = {cell.j, cell.j, cell.j + 1, cell.j + 1};
    for (int c = 0; c < 4; ++c) {
      const std::size_t node = NodeIndex(corner_i[c], corner_j[c], nx);
      const ns_cg::Vec2d corner_pos = field.position(corner_i[c], corner_j[c]);
      const double candidate = seed.arrival_time + (corner_pos - seed.position).norm();
      push_trial(node, candidate);
    }
  }

  // --- Fast Marching main loop. ---
  const int di[4] = {-1, 1, 0, 0};
  const int dj[4] = {0, 0, -1, 1};

  while (!heap.empty()) {
    const HeapEntry entry = heap.top();
    heap.pop();
    const std::size_t node = entry.node_index;
    if (state[node] == NodeState::Known) continue;       // stale entry
    if (entry.arrival_time > distance[node]) continue;    // stale entry

    state[node] = NodeState::Known;

    const std::size_t i = node % nx;
    const std::size_t j = node / nx;

    for (int n = 0; n < 4; ++n) {
      const long ni = static_cast<long>(i) + di[n];
      const long nj = static_cast<long>(j) + dj[n];
      if (ni < 0 || nj < 0 || ni >= static_cast<long>(nx) || nj >= static_cast<long>(ny)) continue;
      const std::size_t neighbor =
          NodeIndex(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj), nx);
      if (state[neighbor] == NodeState::Known) continue;

      const auto known_min = [&](long a_i, long a_j, long b_i, long b_j) {
        double best = kInf;
        if (a_i >= 0 && a_j >= 0 && a_i < static_cast<long>(nx) && a_j < static_cast<long>(ny)) {
          const std::size_t idx =
              NodeIndex(static_cast<std::size_t>(a_i), static_cast<std::size_t>(a_j), nx);
          if (state[idx] == NodeState::Known) best = std::min(best, distance[idx]);
        }
        if (b_i >= 0 && b_j >= 0 && b_i < static_cast<long>(nx) && b_j < static_cast<long>(ny)) {
          const std::size_t idx =
              NodeIndex(static_cast<std::size_t>(b_i), static_cast<std::size_t>(b_j), nx);
          if (state[idx] == NodeState::Known) best = std::min(best, distance[idx]);
        }
        return best;
      };

      const double tx = known_min(static_cast<long>(ni) - 1, nj, static_cast<long>(ni) + 1, nj);
      const double ty = known_min(ni, static_cast<long>(nj) - 1, ni, static_cast<long>(nj) + 1);
      if (!std::isfinite(tx) && !std::isfinite(ty)) continue;  // shouldn't happen

      const double candidate = EikonalUpdate(tx, ty, dx, dy);
      push_trial(neighbor, candidate);
    }
  }

  for (std::size_t j = 0; j < ny; ++j) {
    for (std::size_t i = 0; i < nx; ++i) field.at(i, j) = distance[NodeIndex(i, j, nx)];
  }
  return field;
}

}  // namespace ns_fm2d
