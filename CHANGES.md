# QuadriFlow — Windows Build & Performance Improvements

This document summarises the changes made to the original QuadriFlow codebase
to fix the Windows build, enable OpenMP parallelism, add progress output, fix
hangs on large/irregular meshes, and introduce a fast subdivision-refit
workflow.

---

## Critical bug fixes

### `/O2` optimization level (MSVC)
**Impact: ~13.5× speedup across the entire pipeline.**

`CMakeLists.txt` set `CMAKE_CXX_FLAGS_RELEASE` to `-O3`, which is a GCC/Clang
flag. MSVC silently ignores it (warning D9002) and CMake fell back to
`<Optimization>MinSpace</Optimization>` (`/O1`) in the generated `.vcxproj`.
All previous builds on Windows were compiled at size-optimization level with
no speed optimization whatsoever.

Fix: detect MSVC and use `/O2 /Ob2 /Oi /Ot /DNDEBUG`.

Benchmarked on `BruceLee_Decimated90.obj` (218 MB) at `-f 25000 -subdiv 2`:

| Stage | Before (`/O1`) | After (`/O2`) | Speedup |
|---|---|---|---|
| Initialize | 16.6 s | 2.7 s | 6× |
| Orientation field | 17.5 s | 0.6 s | 29× |
| Position field | 230.9 s | 2.2 s | **105×** |
| Index map | 208.2 s | 30.4 s | 7× |
| **Total** | **585 s** | **43 s** | **~13.5×** |

### Windows build compatibility (`cmake_minimum_required`)
The original `cmake_minimum_required(VERSION 3.1)` in both `CMakeLists.txt`
and `3rd/lemon-1.3.1/CMakeLists.txt` is incompatible with CMake ≥ 3.27.
Bumped to `VERSION 3.5` and added `cmake_policy(SET CMP0091 NEW)` to activate
the `CMAKE_MSVC_RUNTIME_LIBRARY` variable, which fixes the MD/MT runtime
mismatch that caused linker errors between lemon and the main binary.

### `TravelField` infinite loop
`TravelField` (used by `-adaptive`'s `EstimateSlope`) had an unbounded
`while (len > 0)` loop. On large or irregular meshes, near-degenerate
triangles (sliver faces, tiny quads) caused `max_len` to round to ~0 in
floating point, so `len -= max_len` never converged — infinite loop on a
single face, causing the tool to appear hung.

Fixes in `src/field-math.hpp`:
- Loop capped at 100 triangle crossings per call (more than sufficient for
  any realistic step size).
- Added a `1e-10` epsilon guard on `max_len` — exits the loop gracefully
  when the ray is grazing an edge rather than looping forever.
- Replaced `exit(0)` on degenerate geometry with a `break` (graceful skip).

---

## OpenMP parallelism

### BUILD_OPENMP — Windows / MSVC support
The original `CMakeLists.txt` only added `-fopenmp` (GCC) when
`BUILD_OPENMP=ON`. MSVC requires `/openmp` and `find_package(OpenMP)`.
Fixed with proper compiler detection so `-DBUILD_OPENMP=ON` works on all
platforms.

### Newly parallelised stages
All loops below now have `#ifdef WITH_OMP / #pragma omp parallel for` guards:

| Function | What was parallelised |
|---|---|
| `optimize_orientations` | Gauss-Seidel phases per level |
| `optimize_positions` | Gauss-Seidel phases per level × iteration |
| `EstimateSlope` | Per-face orientation averaging and slope traversal |
| `ComputeInverseAffine` | Per-face inverse affine computation |
| `optimize_scale` | S initialisation, face-loop triplet assembly (thread-local), RHS build, scale propagation |
| `optimize_positions_dynamic` | `FindNearest` per-vertex, `O_compact` position update |
| `SubdivideAndRefit` | Per-vertex projection onto original mesh |

### `optimize_scale` matrix assembly
Replaced the intermediate `std::map<int,double>` (O(log n) per lookup) with
direct thread-local `Eigen::Triplet` vectors. `setFromTriplets` handles
summing of duplicate entries, making the face loop embarrassingly parallel.
Also fixed the triplet reservation from `F.cols() * 6` to `F.cols() * 24`
(each edge contributes 8 triplets × 3 edges per face).

---

## Progress output

### Stage percentage markers
`src/main.cpp` and `src/optimizer.cpp` now print `N% done` at every
meaningful point in the pipeline:

| Range | Stage |
|---|---|
| 0–4% | Load + Initialize |
| 4–12% | Orientation field (per hierarchy level) |
| 12–91% | Position field (per level × iteration) |
| 91% | Index map start |
| 92–99% | Index map sub-steps (integer constraints, max flow, flip fixing, position refinement, quad extraction) |
| 97–99% | `optimize_positions_dynamic` iterations |
| 99–100% | Write output file |

### Total elapsed time
Prints `Total time: X.X seconds` after the 100% marker.

---

## New feature: `-subdiv N` — Subdivide and Refit

**Usage:** `-subdiv N` (default `2`, pass `0` to disable)

After QuadriFlow produces a low-resolution quad mesh, applies `N` rounds of
midpoint subdivision (each round multiplies face count by 4) and projects
every vertex back onto the original triangle mesh surface.

This is the recommended workflow for producing high-resolution output:

```sh
# Produces ~373k faces in ~43s instead of running directly at 400k
# (which would take hours due to super-linear solver scaling)
quadriflow -i mesh.obj -o out.obj -f 25000 -subdiv 2

# 3 levels → ~1.5M faces
quadriflow -i mesh.obj -o out.obj -f 25000 -subdiv 3

# Disable subdivision
quadriflow -i mesh.obj -o out.obj -f 25000 -subdiv 0
```

**Face count after subdivision:**

| Base faces | `-subdiv 1` | `-subdiv 2` | `-subdiv 3` |
|---|---|---|---|
| 5,000 | ~20k | ~80k | ~320k |
| 10,000 | ~40k | ~160k | ~640k |
| 25,000 | ~100k | ~400k | ~1.6M |

**Implementation** (`src/subdivide-refine.hpp`):
- `ClosestPointOnTriangle` — exact nearest point via Eberly's method.
- `MeshProjector` — 3D grid spatial index sized to ~8 triangles/cell,
  auto-capped at 8M cells (~64 MB). O(1) average query; early shell-exit
  stops the search once a hit is confirmed.
- `SubdivideQuadMesh` — one midpoint-subdivision pass with shared edge
  midpoints via `unordered_map`.
- `SubdivideAndRefit` — calls subdivision `N` times, builds index once,
  projects all vertices in parallel with OMP.

---

## Linear solver: ConjugateGradient

Replaced `Eigen::SimplicialLLT` (serial Cholesky factorization, O(n^1.5)) with
`Eigen::ConjugateGradient` in `optimize_positions_dynamic` and
`optimize_positions_fixed`.

Key properties:
- **No factorization** — CG is iterative, avoiding the O(n^1.5) setup cost.
- **Warm-start** — each outer iteration initializes CG from the previous
  solution (vertex positions in local frame), reducing iteration counts in
  later passes.
- **DiagonalPreconditioner** (Eigen default) — safe on these mesh Laplacian
  systems; the IncompleteCholesky preconditioner was tried by the original
  authors but crashed.

Observed CG convergence (BruceLee, `-f 25000`):

| Outer iter | CG iters | Error |
|---|---|---|
| 1 | 64 | 9.5e-5 |
| 5 | 71 | 9.4e-5 |
| 10 | 57 | 9.5e-5 |

All well within the 300-iteration cap at 1e-4 tolerance. At small face counts
(25k) the matrices are compact enough that Cholesky and CG take similar wall
time. The benefit grows at very large face counts (100k+) where Cholesky
factorization dominates.

CG iteration counts are printed when built with `-DBUILD_LOG=ON`.

---

## Profiling: index map bottleneck analysis

Per-step timing is available in `ComputeIndexMap` and
`optimize_integer_constraints` when built with `-DBUILD_LOG=ON`
(`cmake ... -DBUILD_LOG=ON`). The profiling macros compile to no-ops in
normal Release builds.

### Findings (BruceLee_Decimated90.obj, 218 MB)

**25k faces** (index map = 31 s total):

| Sub-step | Time | Share |
|---|---|---|
| ComputeMaxFlow (Boykov-Kolmogorov) | 24.3 s | **78%** |
| optimize_positions_dynamic (CG ×10) | 1.8 s | 6% |
| FixFlipHierarchy | 1.0 s | 3% |
| AdvancedExtractQuad | 1.0 s | 3% |
| BuildIntegerConstraints | 1.0 s | 3% |

**200k faces** (index map = 222 s total):

| Sub-step | Time | Share |
|---|---|---|
| ComputeMaxFlow | 154 s | 69% |
| optimize_positions_dynamic (CG ×10) | 23 s | 10% |
| AdvancedExtractQuad | 7.0 s | 3% |
| FixFlipHierarchy | 7.5 s | 3% |

### Root cause

The Boykov-Kolmogorov max flow algorithm scales super-linearly with problem
size: 24 s at 25k grows to 154 s at 200k (≈ O(n^1.1) empirically, but with
a large constant). It is inherently serial and has no parallel implementation.

### Hierarchical max flow — attempted, not merged

A hierarchical approach (`level=-1` in `DownsampleEdgeGraph`) was tested:
supply at the finest level dropped from 840 to 2 at 25k, cutting flow time
from 24 s to 4 s. At 100k the coarse-level flow pushes some edge diff values
beyond |diff|=1 — a constraint that `subdivide_edgeDiff` enforces strictly.
All three corrective strategies tried (clamped propagation, clamped capacity,
selective reset) either failed to reduce supply enough or created larger
imbalances at finer levels. A correct hierarchical max flow for this specific
problem requires proper restriction and prolongation operators that preserve
the |diff|≤1 invariant; this remains future work.

### Practical recommendation

For high face counts, use the `-subdiv N` workflow instead of solving directly
at the target resolution. Running at 25k + 2 subdivisions produces ~400k faces
in ~43 s vs hours for a direct 400k solve.

---

## Windows build instructions

See [`README.md`](README.md#windows-build-visual-studio--msvc) for the
complete step-by-step guide including vcpkg dependency installation.

### Quick reference

```powershell
# Dependencies (one-time)
C:\vcpkg\vcpkg.exe install eigen3:x64-windows boost-filesystem:x64-windows `
    boost-system:x64-windows boost-program-options:x64-windows `
    boost-graph:x64-windows

# Configure
cmake -S . -B build `
    "-DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake" `
    "-DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDLL" `
    -DBUILD_OPENMP=ON

# Build
cmake --build build --config Release --parallel

# Add vcpkg DLLs to PATH (one-time, user scope)
[Environment]::SetEnvironmentVariable(
    'PATH',
    [Environment]::GetEnvironmentVariable('PATH','User') + ';C:\vcpkg\installed\x64-windows\bin',
    'User'
)
```
