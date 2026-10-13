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

The faces remain tied in both transverse directions. Each section owns a control node with only DOF 1 active. Its scalar displacement
q satisfies `(u_B-u_A) dot axis = q` for every pair. Force control applies the
entire force to q, without additional pair loads; displacement control prescribes
q. No artificial mass or stiffness is added. Native sections create a unique
internal node unless `NODE=number` or a singleton NSET names an existing isolated
node. Control-node coordinates do not affect the physics. Element-connected
nodes, shared control nodes and local transformations are rejected. A section without a
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

## Abaqus-compatible input

A complete companion deck is
[`096_bolt_pretension_abaqus.inp`](../examples/09_constraints/096_bolt_pretension_abaqus.inp).
Define an isolated node explicitly, then associate it with a merged cut surface:

```text
*NODE
900001, 0, 0, 0
*SURFACE, TYPE=ELEMENT, NAME=CUT
2, S1
*PRE-TENSION SECTION, SURFACE=CUT, NODE=900001, FOLLOWER=NO
```

`NODE` accepts a scalar compiled reference or an NSET with exactly one node.
Unlike Abaqus's optional node declaration, FEMaster currently requires the node
to exist. Undefined references are rejected, rather than creating a node in an
already compiled Part/Instance map. Native input may still generate its own
internal control node by omitting NODE.

The Abaqus section card is supported at root and Assembly scope. Existing
Part-defined surfaces and nodes can be referenced there as `INSTANCE.CUT`,
`INSTANCE.900001` or `INSTANCE.CONTROL_SET`; assembly NSET aliases also work.
Section cards inside PART or INSTANCE are rejected: automatic per-instance
section creation and local-normal rotation are not implemented. Define one
assembly section and one isolated node for each intended cut. Explicit normal
components are global assembly components.

An optional single data line gives a finite nonzero normal, e.g. `0, 0, -1`.
Without it, the normal is the area-weighted average of the outward face normals;
faces must have consistent orientation. This follows the Abaqus convention:
a positive force acts along that normal on the underlying elements (the rewired
B side of the merged cut), with the opposite force on side A. Positive U1 is
`(u_B-u_A) dot normal`, a tightening adjustment. The normals of native merged
sections and corresponding Abaqus sections therefore give the same sign.
The tests check outward -Z, rotated -X and explicitly reversed +Z directions.
These are analytical and native-equivalence checks; no Abaqus executable was
used for a cross-solver comparison.

Apply the scalar force in a linear static step:

```text
*STEP, NAME=TIGHTEN
*STATIC
1, 1
*CLOAD
900001, 1, 1000
*END STEP
*STEP, NAME=HOLD
*STATIC
1, 1
*BOUNDARY, FIXED
900001, 1, 1
*END STEP
```

Alternatively replace CLOAD with `*BOUNDARY` and
`900001, 1, 1, 0.001` to prescribe U1. FIXED is supported only for Abaqus
pretension-node DOF 1, after a successfully solved prior adjustment, without a
magnitude or amplitude. It snapshots U1 at the start of that step. Its resulting
ordinary boundary condition persists until replaced or removed. Repeated FIXED
snapshots the latest successful U1; an unexecuted or failed solve cannot supply
it. General structural BOUNDARY,FIXED is outside this implementation.

CLOAD and BOUNDARY use their existing independent history families:

- `OP=MOD` replaces the addressed source/DOF and retains other definitions.
  Scalar node references and singleton NSET aliases of one control node share
  a normalized history target, so changing between them replaces the same slot.
- `CLOAD,OP=NEW` removes previous concentrated-load targets, including a
  pretension force, without removing a held boundary.
- `BOUNDARY,OP=NEW` removes previous boundaries, including FIXED. Restate all
  structural supports that must remain; the control DOF becomes free again.
- FIXED does not delete a carried force. With force F and a held U1, the
  scalar reaction is `TF-F`; the force is applied exactly once.
- Without FIXED or a displacement boundary, the carried force remains a force
  prescription. External loading then changes U1 while TF remains prescribed.
- Existing relative STEP TIME amplitudes apply at the linear step end. BOUNDARY
  magnitudes are scaled when materialized, and CLOAD uses its existing amplitude
  evaluation. TOTAL TIME and nonlinear evolution are not supported here.

Before any force/boundary has controlled an Abaqus section, FEMaster preserves
native zero-gap behavior. **This differs from Abaqus's initially free control
DOF.** After the first control definition, removing it releases q (even if all
loads are removed); the section does not silently revert to its initial tie.

Native PRETENSION actions and Abaqus node conditions cannot be mixed on one
section. A force and a displacement on an Abaqus node are valid together;
contradictory displacement equations are rejected by the supported exact
constraint methods. DOFs 2–6, local nodal transformations, CLOAD orientations,
reused or element-connected control nodes, zero/invalid normals, ELEMENT-based
beam/truss sections and FOLLOWER=YES are rejected. Scalar control-node loads
also require inertia relief and load rebalancing to be off.

Reference conventions: [Abaqus 2025 section keyword](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-pre-tensionsection.htm)
and [assembly loads](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEPRCRefMap/simaprc-c-pretension.htm).

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

- `PTGAP`: measured axial gap, written on both nodes of each interface pair and
  on the control node (also available as U1).
- `PTFORC`: force transmitted by the interface, including constraint reactions
  after locking. Sum the vectors on one face and project onto the section axis
  to obtain the signed axial resultant. The control node stores this resultant
  in component 1; `PretensionSection::total_force` exposes it to C++ callers.
  Other nodes contain zero. RF alone is not the total force at the control node.

Keep external supports away from the cut: a support reaction at an interface
node cannot be separated from the cut reaction by this nodal result recovery.

`TF` is a scalar nodal field containing the total signed axial section force
at each control node, zero elsewhere. It satisfies `TF = RF1 + applied CLOAD1`;
U1 gives q. RF at the control node includes the generalized holding reaction,
while physical interface PTFORC retains the original nodal traction recovery.

All three fields follow the existing output selection and are available in RES, FRD
and FEMR. To request just displacement and pretension results in a load case:

```text
*NODE FILE
U, RF, TF, PTFORC, PTGAP
```

The RES section header identifies the scalar control node and DOF 1.
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

Additional checks cover the single active control DOF, explicit and generated
nodes, native/Abaqus equivalence for both exact methods, scalar U/RF/TF,
RES/FRD/FEMR output, two independent sections (one rotated), Assembly/Instance
references, singleton NSETs, explicit and reversed normals, unloaded behavior,
force carryover, CLOAD MOD/NEW, BOUNDARY NEW, FIXED with carried force and
amplitudes. Invalid node references, shared control nodes, structural connectivity,
transformations, DOFs and conflicting adjustments are rejected.

With `FEMASTER_BUILD_TESTS=ON`, the feature checks can be built independently:

```sh
cmake --build build/upstream-master-openmp --target FEMaster FEMaster_pretension_test
ctest --test-dir build/upstream-master-openmp -R Pretension --output-on-failure
```

Replace the build directory with the one used by your local configuration.

The solver-free baseline for issue #2 is recorded in [pretension-benchmark.md](pretension-benchmark.md).
