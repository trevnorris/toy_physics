# Build directive — S11c-b slab/face maps on a one-direction background (Part 1) + round-background scope note (Part 2)

**Status:** **v4** (2026-09-30) · Codex-revised from the orchestrator-written v1–v3 · baseline commit
`c1e96e76` · ⛔ not governing until review-cleared.

You are the builder. Your job is: **build → run → report → stop.**
- ⛔ Do not launch, call or spawn any other AI, agent, reviewer or review process.
- ⛔ Do not commit.
- Write only the files named under *Deliverables*.

## 0 · Declared freezes and execution safeguards

**Put these freezes in the first lines of the report.**

1. `v_bulk_normal_0`, the bulk drain, “appears in no derived operator”
   (`research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md:90–91`). Part 1 has no drain flow.
2. S11c-b “performs no curved-bulk response solve” (spec `:95–97`). Part 1 leaves the bulk at its supplied
   face-trace operands.

**Safeguards.**

1. Run every engine, constructor and comparator job under `scripts/s11c_guarded_run.py`, as `AGENTS.md:8–12`
   requires. Use the default whole-job 2 GiB cap, zero swap, the process/CPU controls, one job at a time, no
   overlapping CAS and no unguarded fallback. Follow the standing no-deadline policy (`AGENTS.md:24–33`): do not
   add a wall-clock, CPU-time, native-alarm, inactivity or progress-stall deadline. A containment refusal or
   contained failure is a result; record the stage and available resource receipt, then stop that job.
2. The suspended-job record says not to launch another job alongside PID 4097233
   (`research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md:5`). Run
   `ps -p 4097233` without signalling it. If it exists in any state, read sources, write the Part 2 scope note
   without its probe, stop before any engine/constructor/comparator launch, and report this stop.
3. Run at most one Wolfram kernel at a time. Preserve its durable log and progress receipts.

## 1 · Supplied physics — unfalsifiable within this build

The supplied object is the S11c-b variable-coefficient slab operator and face response governed by
`research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`, which amends `S11c_a_SHARED_PHYSICS.md` and
`S11b_SHARED_PHYSICS.md`.

The engines define the physical slab coordinates as

```text
X_slab = (u_1, u_2, u_3, theta, e_W, zeta_c).
```

`u_1,u_2,u_3,theta,e_W` are defined in
`research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:143–169,214–233`; `zeta_c` is defined at
`:569–577`. These are the engine coordinates; do not substitute a paraphrased `ζ_±` list. The two independent
face degrees of freedom and their relation to `δW,e_W,ζ_c` are supplied at spec `:59–71`.

The supplied face equations are (spec `:145–154`):

```text
J_s = Λ_A(ω) 𝒜_s + Λ_V(ω) V_s,
Λ_I(ω) = Λ_I⁰/(1-iωτ_I),  I ∈ {A,V,X},
𝒜_s = μ_s - δp_s/ρ_m,          μ_s = μ_θ/ρ_br⁰,
t_s = -(δp_s + Λ_X(ω)𝒜_s) n̂_s,
n̂_s·v_bulk,s = V_s + J_s/ρ_m,
∂_tΣ + ∇_x·(Σv) = -(J_+ + J_-),
Σ = ρ_4D W,                     v = ∂_t u,
δ_vΣ_mat = 0,
δ_vθ + δ_ve_W + ∇_x·δ_vu = 0   (uniform linearisation).
```

Obtain the equations of motion by the supplied method at spec `:152–154`: balance laws, the binding
virtual-displacement rule, variational derivatives with held-fixed fields named, and prescribed external
virtual work—not by placing an irreversible response kernel in an ordinary action.

Keep the accepted energy object intact: 40 records = 10 uniform + 15 `W_BG` first-jet + 15 `MU_R_BG`
first-jet records (step record `research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md:35–37`).
The SymPy engine constructs each accepted record with its coefficient slot at
`scripts/S11c_b_brane_operator_sympy_audit.py:1760–1837` and emits the retained basis at `:4119–4141`.
Before construction, have each engine emit an ordered manifest of those 40 engine records—record label, exact
`COEFFICIENT` slot, its free-symbol list and defining file:line—and use that manifest as the coefficient domain.
Keep every manifested coefficient expression and free symbol live. Also keep live both bookkeepers `eta_bg` and
`sigma_W`, both density representatives, both
anchorings `LAB_HELD` and `MATERIAL_ADVECTED`, and every background jet admitted by R1 below. Retain first order
in wave amplitude and first order in each background bookkeeper. No numerical model point is part of this build.

This supplied physics is not an acceptance test. The build does not test or revise it.

## 2 · Part 1 object and domains

### 2.1 Index convention

Directions 1, 2 and 3 are engine symbol-name labels. SymPy defines `u_1,u_2,u_3` at
`scripts/S11c_b_brane_operator_sympy_audit.py:214–219` and the profile first-jet labels at `:182–199`.
Wolfram defines `spatialCoordinates={xOne,xTwo,xThree}` at
`research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl:239–242`.
SymPy’s `DIRECTIONS=range(3)` at `:57` is a separate 0-based loop index. The report gives each engine’s mapping,
with file:line.

### 2.2 Background and perturbation domains

**R1 is a symmetry/domain restriction, applied before construction.** Every background datum—scalar, polar,
axial, tensor, support, face, bulk-background, domain and boundary datum—is invariant under the full `O(2)` of
rotations and reflections about direction 1, and every background profile depends on `x_1` alone. For profile
jets only, set a jet to zero when any spatial-derivative index has symbol-name label 2 or 3; retain every allowed
direction-1 jet through the engine’s carried order. Do not apply component-name zeroing to general vectors or
tensors. Implement the group-fixed domain itself in each engine before the energy, background-jet maps, face
substrate or variations are constructed. The report states the implementation with file:line.

**Class P.** Every perturbation and every trial/test field is independent of `x_3`. Its `x_2,t` dependence is
`exp(i(k_2 x_2-omega t))`, with `k_2` and `omega` symbolic and live; its `x_1` dependence remains general.

### 2.3 Formal physical-input domain

Use this directive order; every name below comes from the engine coordinate definitions. The slab portion is

```text
u_1, u_2, u_3, theta, e_W, zeta_c
```

and the face-trace portion is

```text
delta_p_plus
d_w_delta_p_plus
delta_v_bulk_plus_1
delta_v_bulk_plus_2
delta_v_bulk_plus_3
delta_v_bulk_plus_4
d_w_delta_v_bulk_plus_1
d_w_delta_v_bulk_plus_2
d_w_delta_v_bulk_plus_3
d_w_delta_v_bulk_plus_4
delta_p_minus
d_w_delta_p_minus
delta_v_bulk_minus_1
delta_v_bulk_minus_2
delta_v_bulk_minus_3
delta_v_bulk_minus_4
d_w_delta_v_bulk_minus_1
d_w_delta_v_bulk_minus_2
d_w_delta_v_bulk_minus_3
d_w_delta_v_bulk_minus_4
```

The trace values and normal jets are defined in
`research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py:95–100,150–165`; S11c-b binds the same names at
`scripts/S11c_b_brane_operator_sympy_audit.py:484–524`. They enter the physical shifted traces at S11c-a
`:578–592,629–646`. These 20 trace coordinates define a **formal operand domain**, not 20 independent physical
bulk degrees of freedom.

### 2.4 Supplied potential-flow pullback

Also emit the pullback of the formal map to the supplied rest-frame bulk equations (spec `:95–97`). Write
`Psi_s=∂_w Phi_s` for the trace of
`φ=Phi(x_1,w) exp(i(k_2x_2-omega t))`. For each face, use the supplied equations

```text
delta_p_s = i rho_m omega Phi_s,
delta_v_bulk_s = (∂_1 Phi_s, i k_2 Phi_s, 0, Psi_s),
∂_w delta_p_s = i rho_m omega Psi_s,
∂_w delta_v_bulk_s =
  (∂_1 Psi_s, i k_2 Psi_s, 0,
   (k_2^2-omega^2/c_s0^2) Phi_s - ∂_1^2 Phi_s).
```

Label the first object `FORMAL_FACE_OPERAND_MAP` and the pulled-back object `POTENTIAL_FLOW_TRACE_MAP`. These are
object names, not result descriptions. The equations in this subsection are supplied physics and are
unfalsifiable within the build.

### 2.5 Physical-output codomain

Flatten component-valued objects in the order written below while retaining the parent engine label. The engine
rows are, in order:

```text
U_BODY_BALANCE[1]
U_BODY_BALANCE[2]
U_BODY_BALANCE[3]
THETA_BALANCE
E_W_BALANCE
ADVECTIVE_MASS_OPERAND
FACE_FLUX_BOUNDARY_OPERANDS
FACE_GENERALIZED_FORCE_ROWS.U[1]
FACE_GENERALIZED_FORCE_ROWS.U[2]
FACE_GENERALIZED_FORCE_ROWS.U[3]
FACE_GENERALIZED_FORCE_ROWS.E_W
FACE_GENERALIZED_FORCE_ROWS.THETA_FACE_FLUX
FACE_GENERALIZED_FORCE_ROWS.CENTER_FACE_GENERALIZED_ROW
FACE_GENERALIZED_FORCE_ROWS.SOURCE_OPERANDS
MU_THETA_FACE_BINDING
```

These labels and their container structure are defined by the SymPy engine at
`scripts/S11c_b_brane_operator_sympy_audit.py:2967–3021,3095–3127`; the generalized-row member names originate at
`:2135–2239`. Preserve each engine’s `LOCAL`, `DIVERGENCE_FLUX`, `SOURCE_OPERAND`, `EVOLUTION_TERM_ORIGINS` and
`EXPANDED` sublabels rather than inventing a replacement row list.

For each ordered pair `(face,dof)`

```text
(plus, DELTA_W), (plus, ZETA_C), (minus, DELTA_W), (minus, ZETA_C)
```

append these engine-defined face objects, in order, flattening vector components 1,2,3,4:

```text
normal_exact[1..4]
measure_exact
face_velocity_exact[1..4]
pressure_trace_exact
bulk_velocity_trace_exact[1..4]
density_trace_exact
current_trace_exact[1..4]
rho4_bg_exact
rhobr_bg_exact
face_normal_raw[tuple_position_1][component_label_1..4]
face_normal_raw[tuple_position_2][component_label_1..4]
face_normal_raw[tuple_position_3][component_label_1..4]
face_measure_raw[tuple_position_1..3]
face_velocity_raw
relative_flux_raw
true_area_flux_raw
pressure_trace_raw
mu_specific_raw
affinity_raw
traction_raw[1..4]
closure_raw
kinematic_raw.OPERAND_A
kinematic_raw.OPERAND_B
kinematic_raw.RESIDUAL
```

The `FaceSource` physical members are defined at
`research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py:600–626`; the named raw objects are
defined at `:878–944`. Do not replace these objects with prose aliases.

Define and print

```text
FORMAL_FACE_OPERAND_MAP := D_(X_slab,X_trace) (physical-output codomain) |_(R1,P),
POTENTIAL_FLOW_TRACE_MAP := FORMAL_FACE_OPERAND_MAP pulled back by §2.4.
```

Print every entry of both maps, including entries computed as zero. Emission is determined only by the ordered
domain/codomain labels, never by an entry’s value.

### 2.6 Separate virtual/test map

The test domain is, in engine order,

```text
delta_v_u_1, delta_v_u_2, delta_v_u_3, delta_v_e_W, delta_v_zeta_c.
```

Those names are defined in S11c-a at
`scripts/S11c_a_interface_geometry_sympy_audit.py:95–100,133–148`; the thickness/centre routing is at `:566–575`.
For the same four `(face,dof)` pairs used in §2.5, the codomain is the four components of the engine member
`virtual_displacement`, defined at `:780–804,840–861`. Define and print its derivative with respect to this test
domain as `VIRTUAL_KINEMATIC_MAP`. It is a separate test-field object, not a row of either physical-input map.
Print every labelled entry.

Construct the stored-energy/face-law weak object once for all physical trial fields and all test fields; then
extract the three named maps. Do not construct a single field’s equation in isolation or split the trial fields
before the weak object exists. The pairing domain is class `P` with compact support in the in-plane interior, as
supplied in spec §3c (`S11c_b_SHARED_PHYSICS.md:312–346`).

## 3 · Supplied FORM controls

Apply R1 and class `P` to the baseline first. Then append each control term below to the stored-energy density and
rerun the complete §2 construction. The control-only datum `e_hat` is introduced after R1 and is not a background
datum in R1’s domain.

```text
K1: a_K1 e_hat_i (∂_k u_i)(∂_k theta),
    e_hat = (sin(beta), 0, cos(beta));

K2: a_K2 theta epsilon_ijk g_i ∂_j u_k,
    g_i = ∂_(y_i) W_bg.
```

Indices `i,j,k` run over engine labels 1,2,3; `epsilon_ijk` is the Levi-Civita symbol. Keep `a_K1`, `a_K2` and
`beta` live, with control coefficients assigned the dimensions required for stored-energy density. For `K1` and
`K2`, print all three complete maps under object tags that name the control and map, followed by the computed
entrywise difference from the baseline. If a complete control payload is byte-identical to its baseline payload,
report that fact after comparison. The controls enter only at stored energy.

## 4 · Engines and comparator

**SymPy import whitelist.** From `scripts/S11c_b_exports.py` or the existing SymPy engine, import only symbol
definitions and the accepted 40-record stored-energy basis. Do not import any operator row, face row, coupling
kernel or other derived payload. R1 must be represented in the imported symbols/basis before construction.

**Blind Wolfram engine.** Import nothing. Re-derive the supplied physics and named objects from the governing
specs; do not transcribe the Python implementation.

**Raw-first comparator.** Join by map name, parent output label, component label and input label. For baseline,
`K1` and `K2`, emit first the raw `operand_PY`, raw `operand_WL` and raw `residual` for every entry, with no
convention map applied. The step record requires the comparator to surface, not normalize, the kinetic whole-row
sign and face-generalized-force whole-row sign (`S11c_b_variable_coefficient_operator.md:112–114`). A separately
tagged convention-mapped diagnostic may follow; state every map and its source file:line. It never replaces the
raw comparison.

If an engine or comparator cannot complete under §0, report the contained result and stop that job. Do not
substitute another computation.

## 5 · The three clauses and structural rule

> **1. The script may PRINT computed objects. It may NOT state conclusions.** An `emit`/`Print` payload is a CAS
> object—an expression, a solved root or a boolean from a symbolic test—not prose describing a result.
>
> **2. PRINT the residual; do NOT assert it.** Compute → emit → then assert.
>
> **3. Interpretation belongs to the STEP RECORD.** The scripts do not editorialise.

> **The only place physical symbols may be combined by hand is the supplied stored energy (including §3), the
> supplied face laws, the supplied potential-flow pullback, and the ansatz/domain. Every other expression involving
> them is reached by computation. Every control re-enters at stored energy, never at a row or result.**

Tag names name objects—field, coordinate, row, map or control—never a value, sign or shape of a result. The scripts,
comparator and report contain no prediction about any entry and no acceptance criterion involving a value.

## 6 · Part 2 — scope note only (⛔ do not build)

The later object is the complete linear operator on a background invariant under `O(3)` in the three in-plane
coordinates about a point, with the bulk response closed. The drain is a live background field: a normal face
component plus a radial in-plane component and no azimuthal component. Its requested map is between the toroidal
trial/test domain and the complete remaining field domain at each angular order for which the toroidal domain
exists.

Write a scope note that:

1. identifies, with file:line, the committed governing source or missing source for the slab energy, face laws and
   responses, curved-bulk closure (including S11c-c1), and drain flow;
2. surfaces how the drain would have to enter the slab, face and bulk equations given spec `:90–91`; do not invent
   that missing specification;
3. states whether the existing three-direction background jets can represent the requested `O(3)` domain at the
   retained order, with file:line;
4. records a small guarded memory/time probe only if §0 permits it, and otherwise records why it was not run;
5. states what cannot be done under the guard.

Then stop. A drain equation or scientific method change requires an orchestrator specification before a build.

## 7 · Deliverables

- `research/pde_ledger_v3/scripts/S11c_d_zinvariant_operator_blocks_sympy.py`, with literal stdout at
  `research/pde_ledger_v3/scripts/out/S11c_d_zinvariant_operator_blocks_sympy.out`
- `research/pde_ledger_v3/mathematica/S11c_d_zinvariant_operator_blocks.wl`, with literal stdout at
  `research/pde_ledger_v3/mathematica/out/S11c_d_zinvariant_operator_blocks.out`
- `research/pde_ledger_v3/scripts/S11c_d_zinvariant_operator_blocks_comparator.py`, with literal stdout at
  `research/pde_ledger_v3/scripts/out/S11c_d_zinvariant_operator_blocks_comparator.out`
- `research/pde_ledger_v3/_measurements/S11c_d_zinvariant_operator_blocks_report.md`, containing:
  - the two declared freezes first;
  - every command, guard invocation, exit code, wall time and peak-memory receipt;
  - output paths;
  - the 40-record manifest;
  - every ordered domain/codomain label and its engine definition file:line;
  - the R1 and class-P implementation in each engine, with file:line;
  - any comparator convention maps and their source file:line;
  - any byte-identical control payloads and contained failures;
  - no interpretation of map entries.
- `research/pde_ledger_v3/_measurements/S11c_d_round_background_scope.md` — the Part 2 scope note.

Stop after writing these.
