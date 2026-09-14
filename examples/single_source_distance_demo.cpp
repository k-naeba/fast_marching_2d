// Minimal usage example: computes the arrival-time field from a single
// seed point at the grid's center and prints how its error against the
// analytic Euclidean distance shrinks as grid resolution increases.

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "fast_marching_2d/fast_marching_2d.hpp"

using ns_cg::Grid2d;
using ns_cg::Vec2d;
using ns_fm2d::ComputeArrivalTime;
using ns_fm2d::Seed;

namespace {

double MaxErrorAgainstAnalyticDistance(const Grid2d<double>& field, std::size_t n,
                                        const Vec2d& center) {
  double max_error = 0.0;
  for (std::size_t j = 0; j < n; ++j) {
    for (std::size_t i = 0; i < n; ++i) {
      const Vec2d p = field.position(i, j);
      const double expected = (p - center).norm();
      max_error = std::max(max_error, std::abs(field.at(i, j) - expected));
    }
  }
  return max_error;
}

}  // namespace

int main() {
  const Vec2d origin(0.0, 0.0);
  const Vec2d center(1.0, 1.0);
  const double extent = 2.0;

  std::printf("Single-source arrival-time field, seed at (%.2f, %.2f)\n", center.x(), center.y());
  std::printf("%8s %12s\n", "grid n", "max |T - dist|");
  for (std::size_t n : {21, 41, 81, 161}) {
    const double spacing = extent / static_cast<double>(n - 1);
    const Grid2d<double> field = ComputeArrivalTime(n, n, spacing, spacing, origin, {Seed{center, 0.0}});
    const double error = MaxErrorAgainstAnalyticDistance(field, n, center);
    std::printf("%8zu %12.6f\n", n, error);
  }
  return 0;
}
