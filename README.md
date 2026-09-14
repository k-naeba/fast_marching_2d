# fast_marching_2d

The Fast Marching Method (Sethian, 1996) for the unit-speed Eikonal
equation `|grad T| = 1` on a regular 2D grid: computes the arrival-time
(geometric distance) field from a set of seed points, each optionally
offset by its own initial arrival time.

Propagation speed is fixed at 1 everywhere -- this is not a general
speed-field / time-of-arrival solver. `T` at each grid node is exactly
the Euclidean distance to the nearest seed, plus that seed's
`arrival_time`.

Built on [`common_geometry`](../common_geometry) (`Vec2d`, `Grid2d`) as
its only dependency.

## Requirements

- CMake >= 3.20
- A C++17 compiler
- [Eigen3](https://eigen.tuxfamily.org/) (e.g. `brew install eigen` on macOS)
- A sibling checkout of [`common_geometry`](../common_geometry) at
  `../common_geometry` relative to this repository

## Building and testing

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
cd build && ctest --output-on-failure
```

Note: because `common_geometry` is pulled in via `add_subdirectory`,
`ctest` also runs its own test suite alongside this project's.

## Usage

```cpp
#include "fast_marching_2d/fast_marching_2d.hpp"

using namespace ns_fm2d;

// nx x ny grid, spacing dx x dy, starting at origin.
const Grid2d<double> field = ComputeArrivalTime(
    nx, ny, dx, dy, origin,
    {Seed{Vec2d(1.0, 1.0), /*arrival_time=*/0.0}});
```

Seed points need not coincide with a grid node: each seed's containing
cell is initialized via direct Euclidean distance to the seed itself
(not the nearest node), so accuracy near a seed isn't limited by the
grid's own resolution.

See
[`examples/single_source_distance_demo.cpp`](examples/single_source_distance_demo.cpp)
for a complete, runnable example.

## What's here

- `fast_marching_2d.hpp`: the entire implementation (`Seed`,
  `ComputeArrivalTime`, and the internal heap-based Fast Marching
  machinery in `detail`).
- `tests/`: GoogleTest, including convergence tests against the
  analytic Euclidean-distance solution (single- and multi-seed), an
  `arrival_time` offset invariance check, and sub-grid seed placement
  accuracy.
- `examples/single_source_distance_demo.cpp`: minimal usage example.
