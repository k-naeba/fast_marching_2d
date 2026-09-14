#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

#include "fast_marching_2d/fast_marching_2d.hpp"

namespace ns_fm2d {
namespace {

using ns_cg::Grid2d;
using ns_cg::Vec2d;

constexpr double kInf = std::numeric_limits<double>::infinity();

// Max absolute error between the computed field and an analytic
// arrival-time function, over every node of `field`.
double MaxErrorAgainstAnalytic(const Grid2d<double>& field, std::size_t nx, std::size_t ny,
                                const std::function<double(const Vec2d&)>& analytic) {
  double max_error = 0.0;
  for (std::size_t j = 0; j < ny; ++j) {
    for (std::size_t i = 0; i < nx; ++i) {
      const double expected = analytic(field.position(i, j));
      max_error = std::max(max_error, std::abs(field.at(i, j) - expected));
    }
  }
  return max_error;
}

TEST(FastMarchingTest, SingleSeedErrorShrinksWithResolution) {
  const Vec2d center(1.0, 1.0);
  const Vec2d origin(0.0, 0.0);
  const double extent = 2.0;
  const auto analytic = [&](const Vec2d& p) { return (p - center).norm(); };

  const auto error_at_resolution = [&](std::size_t n) {
    const double spacing = extent / static_cast<double>(n - 1);
    const Grid2d<double> field =
        ComputeArrivalTime(n, n, spacing, spacing, origin, {Seed{center, 0.0}});
    return MaxErrorAgainstAnalytic(field, n, n, analytic);
  };

  const double error_coarse = error_at_resolution(21);
  const double error_medium = error_at_resolution(41);
  const double error_fine = error_at_resolution(81);

  EXPECT_LT(error_medium, error_coarse + 1e-9);
  EXPECT_LT(error_fine, error_medium + 1e-9);
}

TEST(FastMarchingTest, ArrivalTimeOffsetShiftsWholeFieldByThatConstant) {
  const Vec2d center(1.0, 1.0);
  const Vec2d origin(0.0, 0.0);
  const std::size_t n = 31;
  const double spacing = 2.0 / static_cast<double>(n - 1);
  const double offset = 5.0;

  const Grid2d<double> base = ComputeArrivalTime(n, n, spacing, spacing, origin, {Seed{center, 0.0}});
  const Grid2d<double> shifted =
      ComputeArrivalTime(n, n, spacing, spacing, origin, {Seed{center, offset}});

  for (std::size_t j = 0; j < n; ++j) {
    for (std::size_t i = 0; i < n; ++i) {
      EXPECT_NEAR(shifted.at(i, j), base.at(i, j) + offset, 1e-9) << "at (" << i << "," << j << ")";
    }
  }
}

TEST(FastMarchingTest, MultiSeedErrorShrinksWithResolutionAgainstMinOfBoth) {
  const Vec2d seed_a(0.3, 1.0);
  const Vec2d seed_b(1.7, 1.0);
  const double offset_a = 0.0, offset_b = 0.2;
  const Vec2d origin(0.0, 0.0);
  const double extent = 2.0;
  const auto analytic = [&](const Vec2d& p) {
    return std::min(offset_a + (p - seed_a).norm(), offset_b + (p - seed_b).norm());
  };

  const auto error_at_resolution = [&](std::size_t n) {
    const double spacing = extent / static_cast<double>(n - 1);
    const Grid2d<double> field = ComputeArrivalTime(
        n, n, spacing, spacing, origin, {Seed{seed_a, offset_a}, Seed{seed_b, offset_b}});
    return MaxErrorAgainstAnalytic(field, n, n, analytic);
  };

  const double error_coarse = error_at_resolution(21);
  const double error_fine = error_at_resolution(81);
  EXPECT_LT(error_fine, error_coarse + 1e-9);
}

TEST(FastMarchingTest, SubGridSeedInitializationIsExactAtContainingCellCorners) {
  // Seed deliberately off a grid node (dx=dy=0.1, seed not a multiple of
  // 0.1): the 4 corners of its containing cell should hold the exact
  // Euclidean distance to the seed, not a value degraded by snapping the
  // seed to its nearest node first.
  const Vec2d origin(0.0, 0.0);
  const double spacing = 0.1;
  const std::size_t n = 21;
  const Vec2d seed_pos(1.05, 1.03);

  const Grid2d<double> field = ComputeArrivalTime(n, n, spacing, spacing, origin, {Seed{seed_pos, 0.0}});

  const std::size_t ci = 10, cj = 10;  // floor(1.05/0.1)=10, floor(1.03/0.1)=10
  const std::size_t corner_i[4] = {ci, ci + 1, ci + 1, ci};
  const std::size_t corner_j[4] = {cj, cj, cj + 1, cj + 1};
  for (int c = 0; c < 4; ++c) {
    const Vec2d corner_pos = field.position(corner_i[c], corner_j[c]);
    const double expected = (corner_pos - seed_pos).norm();
    EXPECT_NEAR(field.at(corner_i[c], corner_j[c]), expected, 1e-9) << "corner " << c;
  }
}

TEST(FastMarchingTest, EmptySeedsLeavesEveryNodeAtInfinity) {
  const Grid2d<double> field = ComputeArrivalTime(5, 5, 1.0, 1.0, Vec2d(0, 0), {});
  for (std::size_t j = 0; j < 5; ++j) {
    for (std::size_t i = 0; i < 5; ++i) {
      EXPECT_EQ(field.at(i, j), kInf) << "at (" << i << "," << j << ")";
    }
  }
}

TEST(FastMarchingTest, TooSmallGridDoesNotCrash) {
  EXPECT_NO_THROW({
    const Grid2d<double> field_1x5 =
        ComputeArrivalTime(1, 5, 1.0, 1.0, Vec2d(0, 0), {Seed{Vec2d(0, 0), 0.0}});
    const Grid2d<double> field_5x1 =
        ComputeArrivalTime(5, 1, 1.0, 1.0, Vec2d(0, 0), {Seed{Vec2d(0, 0), 0.0}});
    const Grid2d<double> field_0x0 = ComputeArrivalTime(0, 0, 1.0, 1.0, Vec2d(0, 0), {});
  });
}

TEST(FastMarchingTest, SeedOutsideGridDomainDoesNotCrashAndStillPropagates) {
  const std::size_t n = 11;
  const Grid2d<double> field =
      ComputeArrivalTime(n, n, 1.0, 1.0, Vec2d(0, 0), {Seed{Vec2d(-50.0, -50.0), 0.0}});
  for (std::size_t j = 0; j < n; ++j) {
    for (std::size_t i = 0; i < n; ++i) {
      EXPECT_TRUE(std::isfinite(field.at(i, j))) << "at (" << i << "," << j << ")";
    }
  }
}

}  // namespace
}  // namespace ns_fm2d
