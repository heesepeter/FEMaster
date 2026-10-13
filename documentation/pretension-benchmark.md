# Issue #2: mesh-preparation baseline (before issue #4)

Measured on 2026-10-13, branch `feature/bolt-pretension`, physics/topology
revision `f1213472`. No issue #2 optimizations have been applied.

The opt-in executable in `benchmarks/pretension_prepare.cpp` constructs a
structured C3D8 bar, three elements long and n by n elements across, with one
merged planar cut at z=1. A fresh model is built and compiled for every run.
Only `prepare_pretension_sections()` is timed using steady_clock; construction,
compile, material assignment, field/output setup and equilibrium solvers are
excluded. No solver is invoked. Timing includes the preparation diagnostic log.
Repeat 0 is warmup; repeats 1–5 are reported. Wall-clock measurements are machine
and background-load dependent, not performance assertions.

| n | Input nodes | Elements | Cut faces | Pairs / generated nodes | Final nodes | Median ms | Min–max ms |
|---|---:|---:|---:|---:|---:|---:|---:|
|4|100|48|16|25|125|0.123|0.120–0.202|
|8|324|192|64|81|405|1.386|1.305–2.227|
|16|1156|768|256|289|1445|14.471|13.139–17.330|
|32|4356|3072|1024|1089|5445|191.938|185.516–195.540|
|64|16900|12288|4096|4225|21125|3625.087|3446.392–3636.374|

Raw data: `benchmarks/pretension_baseline_f1213472.csv`.

Build: Darwin 25.6.0 arm64, Apple clang 21.0.0 (clang-2100.1.1.101),
C++17, CMake Release (`-O3 -DNDEBUG`), double precision ON, fast Eigen compile
ON (`EIGEN_UNROLLING_LIMIT=50`), Eigen no-debug ON, OpenMP ON, Accelerate ON,
MKL/CUDA/cuDSS OFF, static runtime OFF, tests ON. `OMP_NUM_THREADS=1`.
Dependency revisions in `build/upstream-master-openmp`:

| Dependency | Revision |
|---|---|
|Eigen|`bc3b39870ecb690a623a3f49149a358b95c5781d`|
|Spectra|`6841bcbacaa0f0a8446210314e682057a084be4e`|
|argparse|`3eda91b2e1ce7d569f84ba295507c4cd8fd96910`|
|LZ4|`ebb370ca83af193212df4dcbadcc5d87bc0de2f0`|
|GoogleTest|`063de7e9578f82b369302001269680b4b1553359`|

The Eigen headers report 3.5.0; the configured dependency tag is 5.0.1.
Use the recorded revision, rather than the displayed version alone, for exact
reproduction with cached dependency checkouts.

Reproduce on the baseline revision with the benchmark source/CMake target from
step 1 of issue #4 (the target itself does not change preparation):

```sh
cmake -S . -B build/upstream-master-openmp -DCMAKE_BUILD_TYPE=Release \
  -DFEMASTER_BUILD_TESTS=ON -DFEMASTER_BUILD_BENCHMARKS=ON \
  -DFEMASTER_DOUBLE_PRECISION=ON -DFEMASTER_FAST_EIGEN_COMPILE=ON \
  -DFEMASTER_ENABLE_OPENMP=ON -DFEMASTER_ENABLE_ACCELERATE=ON \
  -DFEMASTER_ENABLE_MKL=OFF -DFEMASTER_ENABLE_CUDA=OFF
cmake --build build/upstream-master-openmp --target FEMaster_pretension_benchmark -j 4
OMP_NUM_THREADS=1 build/upstream-master-openmp/bin/FEMaster_pretension_benchmark
```

For a baseline checkout use `f1213472` and copy only the benchmark source and its
opt-in CMake target. Later issue #4 runs include one additional control node per
section; they must not be labeled as the original baseline. Pair matching and
all topology algorithms remain unchanged by issue #4.
