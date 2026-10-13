# Issue #2 — mesh preparation optimization (2026-10-13)

Branch: `feature/bolt-pretension`; direct reference: `914c3a11` (after #4).
Issues #3 (`0fbcca7f`), #1 (`f1213472`) and #4 (`49888970`, `914c3a11`)
were read together with current GitHub issues #1–#4. No remote writes are part
of this work. The original pre-#4 reference and raw CSV remain in
[pretension-benchmark.md](pretension-benchmark.md).

## Implementation and invariants

- `ModelData::append_nodes` resizes each compiled position field once per section,
  reserves the result mapping once, and scans existing semantic IDs once. Cut
  copies remain ordered by original node ID; the generated native control follows
  them, exactly as before. Explicit Abaqus controls are already compiled and are
  not regenerated. All-node sets and reference positions include every new node.
- One dense node-to-element incidence table per section serves explicit-control
  attachment checks, flood fill, retained-side validation and rewiring. Only
  incident elements are visited for each cut node. The table is rebuilt for the
  next section because the preceding cut changes connectivity.
- Geometry/normal validation is a separate step before allocation or rewiring.
  Region separability and affected solid types are also checked before mutation.
  Every owned boundary face on side B follows its owner's new connectivity.
  Validation is per section; a later invalid section does not roll back earlier
  successfully prepared sections. This is not a global transaction guarantee.
- Nearest-node matching and its deterministic tie order remain unchanged.
  No spatial search or additional dependency is introduced. The scalar control
  formulation, input syntax, load-step history, exact-constraint decision,
  output writers and physics code are unchanged.

## Validation

Before editing, all 37 existing Pretension CTest checks passed. After the change,
all 39 pass. New checks cover an oblique 16×16 cut (289 pairs), all six boundary
faces of each solid, high/non-dense semantic result IDs, reference positions,
control-node placement and repeated preparation. A nonplanar cut is rejected
without adding nodes, changing element connectivity or allocating its control.
The existing exterior-cut rejection now also checks unchanged node storage.

Existing tests continue to cover two independent rotated sections, native and
Abaqus inputs, explicit/generated scalar controls and DOF numbering,
force/displacement/LOCK with later external load, NULLSPACE and ELIMINATION,
invalid/mixed/redundant constraints, and RES/FRD/FEMR values and metadata.
No Lagrange support is added; issue #1's restriction remains intact.

The isolated feature suite is used because the previously documented full-suite
build failures in obsolete structural-element tests are outside this change;
see [issue #4 validation](pretension-issue4-validation.md).

## Measurement protocol

All final timings run sequentially, without a concurrent build or another
benchmark. Each size/revision uses a fresh model: repeat 0 is warmup, repeats
1–5 are the recorded samples. Medians and min–max ranges are descriptive,
not test thresholds. Models are identical to the original baseline: C3D8 bars,
three elements long, n×n elements across, a merged planar cut at z=1.
There is one additional scalar control node after #4.

Preparation is timed only around `prepare_pretension_sections()` with
`steady_clock` (includes its diagnostic logging). Construction, compile,
material assignment and equilibrium solving are excluded. The direct baseline
is rebuilt from `914c3a11`, using its original `model_data.cpp` and
`model_pretension.cpp` in an otherwise identical core archive; those are the only
changed core translation units, and the header change is a method declaration,
not a layout change. Executables are linked with the same benchmark/main objects
and dependency archives. An old benchmark executable initially found in the
build directory still produced pre-#4 node counts; its measurements were discarded.

Solver timings come from the existing `solving constraint system` timer in the
CLI, separately from preparation, stiffness assembly, constraint transformation,
recovery and file output. Its resolution is integer milliseconds, so 0 means
less than 1 ms. The same five meshes use E=210000, nu=0, area=1, length=3,
F=1000, fixed end faces, CPU DIRECT and ELIMINATION. Each solve checks
q=FL/(EA)=0.0142857142857 and total cut force=1000. Every RES field row and
semantic ID is also compared between revisions (absolute tolerance 1e-5,
relative tolerance 1e-6); roundoff near zero is allowed.

`pretension_matching.cpp` isolates the unchanged nearest-node kernel with the
same Eigen 3D distance, ordered sets, map insertion and erase operations. It
includes matching-container setup but excludes surface gathering, topology and
position-field allocation. This is a separate kernel experiment, not a precise
in-situ subtraction from total preparation time.

Build: Darwin 25.6.0 arm64, Apple clang 21.0.0, C++17, Release `-O3 -DNDEBUG`,
double precision ON, fast Eigen compile ON (`EIGEN_UNROLLING_LIMIT=50`), Eigen
no-debug ON, OpenMP/Accelerate ON, MKL/CUDA/cuDSS/static runtime OFF,
tests/benchmarks ON; `OMP_NUM_THREADS=1`, CLI `--ncpus` default 1.
The cached dependency revisions match the pre-#4 benchmark document.

## Reproduction

Use identical settings from [the original protocol](pretension-benchmark.md)
for both the `914c3a11` checkout and the optimized checkout. Build targets
`FEMaster`, `FEMaster_pretension_test`, `FEMaster_pretension_benchmark` and
(on the optimized checkout) `FEMaster_pretension_matching`.

```sh
ctest --test-dir build/upstream-master-openmp -R Pretension --output-on-failure
OMP_NUM_THREADS=1 build/upstream-master-openmp/bin/FEMaster_pretension_benchmark
OMP_NUM_THREADS=1 build/upstream-master-openmp/bin/FEMaster_pretension_matching
python3 benchmarks/pretension_solver.py /tmp/pretension-models
python3 benchmarks/pretension_compare.py /path/to/baseline/FEMaster \
  build/upstream-master-openmp/bin/FEMaster /tmp/pretension-models \
  benchmarks/pretension_solver_issue2.csv
```

## Preparation results

One warmup plus five timed samples per row; all times in milliseconds.
Input/final nodes are dense compiled nodes, including the generated control.

| n | Input nodes | Elements | Cut faces | Pairs | Final nodes | After #4 median (min–max) | Optimized median (min–max) | Speedup | Matching median |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|4|100|48|16|25|126|0.163 (0.158–0.188)|0.082 (0.074–0.224)|1.98×|0.008|
|8|324|192|64|81|406|1.307 (1.243–1.351)|0.330 (0.312–0.399)|3.97×|0.033|
|16|1156|768|256|289|1446|12.976 (12.879–13.801)|1.765 (1.454–1.964)|7.35×|0.230|
|32|4356|3072|1024|1089|5446|182.117 (180.886–188.055)|7.588 (6.763–9.231)|24.00×|2.892|
|64|16900|12288|4096|4225|21126|3448.840 (3315.809–3538.210)|82.433 (79.728–84.637)|41.84×|22.298|

Raw samples: `benchmarks/pretension_baseline_914c3a11.csv`,
`benchmarks/pretension_optimized_issue2.csv`, `benchmarks/pretension_matching.csv`.
The old pre-#4 median for n=64 was 3625.087 ms with 21125 final nodes;
it is historical context, not the denominator for the 41.84× direct speedup.
The largest isolated matching median is 22.298 ms versus 82.433 ms for full
preparation. At these measured sizes the simpler unchanged kernel is retained;
it is not the dominant cost, and a spatial index would add matching/tolerance
complexity. Reconsider for larger cuts or a preparation budget below this range.

## Solver and physical results

The equilibrium timer is independent of preparation; there is no claimed solver
speedup. Medians (min–max), in integer milliseconds:

| n | After #4 | Optimized |
|---|---:|---:|
|4|0 (0–0)|0 (0–0)|
|8|0 (0–0)|0 (0–0)|
|16|2 (2–2)|2 (2–2)|
|32|9 (8–9)|10 (9–10)|
|64|46 (45–47)|50 (45–51)|

Raw warmup and repeated measurements, model sizes, q and force are in
`benchmarks/pretension_solver_issue2.csv`. All 60 solves meet the analytical
q/force checks. The field comparator preserves equal NaN markers in
inactive components; they are not treated as physical numerical values.
Every RES row/ID is compared for sample 1 of all five sizes. Maximum absolute
field differences are 1.113e-12, 6.528e-12, 1.102e-11, 1e-7 and 4.464e-11,
respectively. These are within the recorded tolerances. The final 39/39 CTest
run takes 0.39 s. FRD/FEMR values and scalar-control metadata remain covered by
the feature regressions.

## Remaining limits

Matching is still quadratic; incidence storage is linear in total element
connectivity and is rebuilt per section. Batching is per section rather than
across the entire model to preserve ID order and the lifecycle of independent
cuts and explicit controls. Only existing matching planar solid faces are
supported; no automatic cuts through element interiors, nonlinear pretension,
FOLLOWER=YES or external Abaqus-solver comparison is added.
