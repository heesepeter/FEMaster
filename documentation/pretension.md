# Bolt pretension with matching faces

Pretension is currently supported in `LINEARSTATIC` load cases. The interface
must be planar and have matching nodes on both faces. Solid elements may share
the interface nodes initially; FEMaster then duplicates those nodes and rewires
one connected side before creating user fields, constraints and result files.
This separates an existing mesh face; it does not cut through element interiors.

A complete three-step example is in
[`examples/09_constraints/095_bolt_pretension.inp`](../examples/09_constraints/095_bolt_pretension.inp).

## Define the section

For a merged mesh, define one element face set on the cut:

```text
*SURFACE, SFSET=CUT_FACE
2, 2, 1
*PRETENSIONSECTION, NAME=PT1, SURFACE_A=CUT_FACE
```

The surface record above uses the native `surface-id, element-id, face-id`
syntax. Choose faces that completely separate the two sides of the solid mesh.
An exterior face or an incomplete cut is rejected.

For an already separated mesh, specify both coincident face sets:

```text
*PRETENSIONSECTION, NAME=PT1, SURFACE_A=CUT_A, SURFACE_B=CUT_B
```

The sets must have equal node counts and matching coordinates. Overlapping
pretension sections are unsupported. The axial direction follows the normal of
the first face in `SURFACE_A`; use consistent face orientation throughout a
section. Positive `VALUE` applies opposite forces along that direction, or
prescribes a positive relative displacement `(u_B-u_A) dot axis`.

## Load and lock

Inside the first load case, apply the total section force:

```text
*PRETENSION, SECTION=PT1, ACTION=LOAD, CONTROL=FORCE, VALUE=1000
```

Alternatively, prescribe the relative displacement:

```text
*PRETENSION, SECTION=PT1, ACTION=LOAD, CONTROL=DISPLACEMENT, VALUE=0.001
```

The faces remain tied in both transverse directions. Force control leaves one
common axial gap free and distributes the total force equally over the node
pairs. Displacement control prescribes that common gap. A section without a
load action remains tied with zero gap.

In the next load case, lock the solved gap:

```text
*PRETENSION, SECTION=PT1, ACTION=LOCK
```

Force control must complete a solve before it can be locked. The locked gap
persists into later load cases; an additional `LOCK` is harmless. External loads
then change the transmitted section force while preserving the gap. This is a
sequence of linear equilibrium solves; nonlinear bolt/contact behavior is not
included in this version.

## Constraint methods and limitations

Use `*CONSTRAINTMETHOD, TYPE=NULLSPACE` (the default) or `TYPE=ELIMINATION`.
Both enforce the prescribed or locked gap through the existing affine reduction
`u = u_p + T q`, to numerical tolerance. Supports and manual equations may be
combined with pretension. Consistent redundant equations are accepted;
inconsistent equations are rejected before the equilibrium solve and before a
new solved gap is captured. This also applies when newly activated supports
contradict a previously locked gap. Constraint feasibility alone does not
guarantee structural stability.

`LAGRANGE` is currently rejected whenever a pretension section is present,
including an unloaded section. Select `NULLSPACE` or `ELIMINATION` explicitly
if an existing input deck requests Lagrange multipliers. Models without
pretension retain their existing constraint behavior and regularization.

### Decision for the current upstream implementation (issue #1)

The Lagrange implementation scales the constraint rows and assembles
`[K, C_s^T; C_s, -r I]`, with
`r = 1e-10 * max(abs(diag(K)))` (unit fallback for an unusable stiffness scale).
The multiplier diagonal stabilizes redundant constraint rows, but changes the
constraint equation to `C_s u - r lambda = d_s`: prescribed displacements are
therefore approximate. Its build report sets `rank_known=false` and
`feasible=true` without checking either rank or consistency.

Setting `r=0` gives exact multipliers for a suitable independent, consistent
constraint system. It is not reliable for all currently accepted inputs:
duplicate equations make the multiplier system singular, and contradictory
equations have no exact solution. The regression suite demonstrates the
singular duplicate-row case. The previous pretension-specific global switch
also changed support and MPC treatment throughout the model.

`NULLSPACE` already uses sparse QR rank detection and checks the affine residual
against the original equations. `ELIMINATION` substitutes dependent DOFs,
recognizes redundant rows and detects contradictory rows. The pretension path
now also rejects a non-feasible transformer report before solving, rather than
only reporting the residual and continuing.

Selective removal of regularization from pretension rows was considered and
rejected for now: it would still approximate other constraints and would not
resolve dependencies or contradictions spanning supports, MPCs and pretension.
Supporting exact Lagrange later requires rank/consistency handling of the
combined constraint system and validation of the saddle-point solver backends.
No global constraint algorithm is changed in this feature.

## Results

Linear-static default output includes the following nodal fields when a
pretension section exists:

- `PTGAP`: measured axial gap, written on both nodes of each interface pair.
- `PTFORC`: force transmitted by the interface, including constraint reactions
  after locking. Sum the vectors on one face and project onto the section axis
  to obtain the signed axial resultant. Other nodes contain zero.

Keep external supports away from the cut: a support reaction at an interface
node cannot be separated from the cut reaction by this nodal result recovery.

Both fields follow the existing output selection and are available in RES, FRD
and FEMR. To request just displacement and pretension results in a load case:

```text
*NODE FILE
U, PTFORC, PTGAP
```

Generated nodes appear in the result mesh with unique identifiers. User node
sets keep their existing membership; only the global all-nodes region gains
the generated nodes. Boundary surfaces owned by rewired elements follow their
updated node connectivity.

## Verification

The regression tests cover force and displacement control, locking, a subsequent
external load, separate matching faces, invalid cuts, invalid lock ordering,
explicit output selection and binary output. Both supported methods are tested
with force and displacement control followed by LOCK and external loading,
manual MPCs, scaled duplicate equations and equations redundant with supports.
Conflicting support/MPC prescriptions, an equation conflicting with a prescribed
gap and supports conflicting with a locked gap must fail. Lagrange rejection
and the unchanged default multiplier regularization without pretension are
checked separately. A uniaxial bar with Poisson ratio zero provides the analytical reference `gap = force * length / (E * area)`.

With `FEMASTER_BUILD_TESTS=ON`, the feature checks can be built independently:

```sh
cmake --build build/upstream-master-openmp --target FEMaster FEMaster_pretension_test
ctest --test-dir build/upstream-master-openmp -R Pretension --output-on-failure
```

Replace the build directory with the one used by your local configuration.
