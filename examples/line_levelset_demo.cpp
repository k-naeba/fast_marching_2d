// Builds a level-set (signed distance) field over the whole grid from a
// single straight-line interface: samples the line densely as Fast
// Marching seeds to get the *unsigned* distance (fast_marching_2d has
// no notion of sign), then assigns sign via a simple half-plane test
// against the line's own normal. Since the exact signed distance to an
// infinite line is just a linear function of position, this doubles as
// a correctness check on the simplest possible interface before moving
// on to curved ones.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

#include "fast_marching_2d/fast_marching_2d.hpp"

using ns_cg::Grid2d;
using ns_cg::Vec2d;
using ns_fm2d::ComputeArrivalTime;
using ns_fm2d::Seed;

namespace {

// An (infinite) line through `point` along unit direction `dir`, with
// `normal` (dir rotated +90 degrees) defining the positive side.
struct Line {
  Vec2d point;
  Vec2d dir;
  Vec2d normal;
};

Line MakeLine(const Vec2d& a, const Vec2d& b) {
  const Vec2d dir = (b - a).normalized();
  return Line{a, dir, Vec2d(-dir.y(), dir.x())};
}

double SignedDistanceToLine(const Line& line, const Vec2d& p) {
  return line.normal.dot(p - line.point);
}

// Dense seed points along the line, spanning well past the grid's own
// bounding box on both ends (projected onto the line, plus a margin of
// one grid extent) so every node's nearest point on the *infinite* line
// -- not just the visible segment between the two points that defined
// it -- is covered by a sampled seed.
std::vector<Seed> SampleLineSeeds(const Line& line, const Grid2d<double>& grid, std::size_t nx,
                                   std::size_t ny, double spacing) {
  const Vec2d origin = grid.position(0, 0);
  const Vec2d far_corner = grid.position(nx - 1, ny - 1);
  const Vec2d corners[4] = {origin, Vec2d(far_corner.x(), origin.y()), far_corner,
                             Vec2d(origin.x(), far_corner.y())};

  double t_min = std::numeric_limits<double>::infinity();
  double t_max = -std::numeric_limits<double>::infinity();
  for (const Vec2d& corner : corners) {
    const double t = line.dir.dot(corner - line.point);
    t_min = std::min(t_min, t);
    t_max = std::max(t_max, t);
  }
  const double margin = (far_corner - origin).norm();
  t_min -= margin;
  t_max += margin;

  std::vector<Seed> seeds;
  const int steps = static_cast<int>(std::ceil((t_max - t_min) / spacing));
  seeds.reserve(steps + 1);
  for (int s = 0; s <= steps; ++s) {
    const double t = t_min + (t_max - t_min) * static_cast<double>(s) / steps;
    seeds.push_back(Seed{line.point + t * line.dir, 0.0});
  }
  return seeds;
}

}  // namespace

int main() {
  const Vec2d origin(0.0, 0.0);
  const std::size_t n = 101;
  const double extent = 2.0;
  const double spacing = extent / static_cast<double>(n - 1);

  // A simple straight line crossing the domain diagonally.
  const Line line = MakeLine(Vec2d(0.0, 0.4), Vec2d(2.0, 1.6));

  const Grid2d<double> grid_template(n, n, spacing, spacing, origin);
  // Sample seeds at half the grid spacing for sub-grid accurate
  // initialization (see fast_marching_2d.hpp's own note on Seed).
  const std::vector<Seed> seeds = SampleLineSeeds(line, grid_template, n, n, spacing * 0.5);

  const Grid2d<double> unsigned_distance = ComputeArrivalTime(n, n, spacing, spacing, origin, seeds);

  // Combine into a signed level-set field, comparing against the exact
  // analytic signed distance to the line at every node.
  Grid2d<double> levelset(n, n, spacing, spacing, origin);
  double max_error = 0.0;
  for (std::size_t j = 0; j < n; ++j) {
    for (std::size_t i = 0; i < n; ++i) {
      const Vec2d p = unsigned_distance.position(i, j);
      const double exact = SignedDistanceToLine(line, p);
      const double sign = exact < 0.0 ? -1.0 : 1.0;
      const double value = sign * unsigned_distance.at(i, j);
      levelset.at(i, j) = value;
      max_error = std::max(max_error, std::abs(value - exact));
    }
  }

  std::printf("Level-set field from a straight-line interface, %zux%zu grid\n", n, n);
  std::printf("line: (%.2f, %.2f) -> (%.2f, %.2f)\n", line.point.x(), line.point.y(),
              (line.point + line.dir).x(), (line.point + line.dir).y());
  std::printf("max |levelset - exact signed distance to line| = %.6f\n\n", max_error);

  std::printf("probe along x=1.00 (crosses the line near y=%.2f):\n",
              line.point.y() + line.dir.y() / line.dir.x() * (1.0 - line.point.x()));
  std::printf("%8s %12s %12s\n", "y", "levelset", "exact");
  for (double y : {0.0, 0.3, 0.6, 0.9, 1.0, 1.2, 1.5, 1.8, 2.0}) {
    const Vec2d p(1.0, y);
    const auto i = static_cast<std::size_t>(std::lround((p.x() - origin.x()) / spacing));
    const auto j = static_cast<std::size_t>(std::lround((p.y() - origin.y()) / spacing));
    std::printf("%8.2f %12.6f %12.6f\n", y, levelset.at(i, j), SignedDistanceToLine(line, p));
  }
  return 0;
}
