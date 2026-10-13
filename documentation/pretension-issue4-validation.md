# Issue #4 validation — 2026-10-13

Branch: `feature/bolt-pretension`. Starting revision: `f1213472`.
The working tree was clean. Issue #3 was present as `0fbcca7f`; its removed
fields remain absent. Issue #1 was present as `f1213472`; its exact-constraint
decision is retained. No push, PR, issue comment or issue-state mutation was made.

## Before implementation

All 20 existing Pretension CTest checks passed. The solver-free preparation
baseline was measured before changes to the implementation, with five cut sizes
and six fresh models per size (one warmup, five samples). See
[benchmark protocol and results](pretension-benchmark.md) and the checked-in CSV.
No mesh-preparation optimization from issue #2 is included.

## Step 1: scalar control formulation

Local commit: `49888970`.

All 22 feature checks passed after introducing the scalar control DOF, direct
scalar force, affine pair compatibility, successful-solve gap capture and
native LOAD/LOCK mapping. Explicit and generated nodes were checked; only DOF 1
is active. Existing analytical force/displacement/LOCK/external-load, mixed and
contradictory constraints, supported methods and RES/FRD/FEMR checks stayed green.

## Step 2: Abaqus input and history

All 37 Pretension CTest checks pass, including:

- Native/Abaqus force and displacement input compared for all physical nodes,
  interface gaps and forces, across three steps, with NULLSPACE and ELIMINATION.
- Analytical uniaxial reference: E=210000, nu=0, A=4, L=3;
  F=1000 gives q=0.0035714285714285714. FIXED holds q under later external load.
- Carried force without FIXED, carried force plus FIXED, CLOAD MOD/NEW,
  BOUNDARY NEW releasing FIXED, force and boundary amplitudes, and replacement
  through scalar-node/singleton-NSET aliases.
- Two independent bars/sections with different controls and a rotated cut,
  Assembly with Instance-qualified references, singleton NSET, automatic average
  normal and explicit/reversed normals.
- U1, PTGAP, PTFORC, RF1 and TF accounting, binary field values against RES,
  FRD output and explicit output selection. Control-node metadata uses semantic
  result IDs.
- Rejection of undefined/structural/reused control nodes, ambiguous NSETs,
  transformations, DOFs 2–6 (including zero-valued requests), invalid normals,
  FIXED before a successful solve, inconsistent prescriptions, unsupported
  procedure/section/follower/scope variants and scalar-load balancing.

Reproduce the feature checks with:

```sh
cmake --build build/upstream-master-openmp --target FEMaster FEMaster_pretension_test -j 4
ctest --test-dir build/upstream-master-openmp -R Pretension --output-on-failure
```

Both `095_bolt_pretension.inp` and `096_bolt_pretension_abaqus.inp` also ran through
the CLI successfully using their original Poisson ratio 0.3. Both report
q=0.00329623 across all three steps. These examples use the same physical mesh
and loading; the Abaqus deck carries its original CLOAD through FIXED.

## Broader test limitation

An attempted build of `FEMaster_test` fails before the full suite can run.
Existing `tests/test_elements.cpp` calls at lines 91, 165, 221, 266, 298 and 313,
and `tests/test_elements_solid.cpp` calls starting at line 175 use obsolete
stiffness/internal-force function signatures (too few/too many arguments).
Those tests and their structural-element interfaces are unchanged relative to
`f1213472`. They were not repaired as part of this issue. The application,
complete core library and isolated feature test target build successfully.

## Remaining scope

See [compatibility details](pretension.md): linear static matching planar solid
cuts only; NULLSPACE/ELIMINATION only; explicit existing Abaqus node required;
root/Assembly section cards with compiled Instance references; no PART/INSTANCE
section-card expansion, ELEMENT pretension, FOLLOWER=YES, nonlinear procedures,
general structural FIXED, inertia relief or load rebalancing. Initial uncontrolled
Abaqus sections deliberately retain the native zero gap. No external Abaqus solver
comparison was performed; normal/load conventions were checked against the
referenced documentation and analytical/native models. Interface supports retain
the documented force-recovery ambiguity.
