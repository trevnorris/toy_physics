# Build directive — S11c-b slab/face maps on a one-direction background (Part 1) + round-background scope note (Part 2)

**Status:** **v5** (2026-10-01) · Codex-revised (v4–v5) from the orchestrator-written v1–v3 · v4 baseline commit
`5693e861` · ⛔ not governing until review-cleared.

You are the builder. Your job is: **build → run → report → stop.**
- ⛔ Do not launch, call or spawn any other AI, agent, reviewer or review process.
- ⛔ Do not commit.
- Write only the files named under *Deliverables*.

## 0 · Declared freezes and execution safeguards

**Put these freezes in the first lines of the report.**

1. `v_bulk_normal_0`, the bulk drain, “appears in no derived operator”
   (`research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md:90–91`). Part 1 has no drain flow.
2. S11c-b “performs no curved-bulk response solve” (spec `:95–97`). Part 1 leaves the bulk at its supplied
   face-trace and projection operands.

**Safeguards.** These are requirements of this directive.

1. Run every engine, constructor and comparator job—including every Wolfram invocation—through
   `scripts/s11c_guarded_run.py`: default whole-job 2 GiB cap, zero swap, process/CPU controls, one job at a time,
   no overlapping CAS and no unguarded fallback. `AGENTS.md:8–12` independently requires the guard for S11c
   Python constructors, validators and export jobs. Wolfram containment under this guard is not established;
   make the first Wolfram action a minimal guarded startup/containment probe, record its receipt, and stop the
   Wolfram lane if it cannot run inside the guard. Follow the no-deadline policy (`AGENTS.md:27–36`): add no
   wall-clock, CPU-time, native-alarm, inactivity or progress-stall deadline. A containment refusal or contained
   failure is a reportable result; record the stage and available resource receipt, then stop that job.
2. The suspended-job record forbids another job alongside PID 4097233
   (`research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md:5`), and the
   later record says that PID was resumed (`…S11c_d_numerical_radiating_equal_speed_feasibility.md:47`). Run
   `ps -p 4097233` without signalling it. If it exists in any state, write the Part 2 scope note without a probe,
   stop before every engine/constructor/comparator launch, and report this stop.
3. Run at most one Wolfram kernel at a time. Preserve its durable log and progress receipts.

## 1 · Supplied physics — unfalsifiable within this build

The supplied object is the S11c-b variable-coefficient slab operator and face response governed by
`research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`, which amends `S11c_a_SHARED_PHYSICS.md` and
`S11b_SHARED_PHYSICS.md`.

The supplied face equations are (spec `:145–154`):

```text
J_s = Λ_A(omega) 𝒜_s + Λ_V(omega) V_s,
Λ_I(omega) = Λ_I^0/(1-i omega tau_I),  I in {A,V,X},
𝒜_s = μ_s - delta_p_s/rho_m,           μ_s = μ_theta/rho_br^0,
t_s = -(delta_p_s + Λ_X(omega)𝒜_s) n_hat_s,
n_hat_s·v_bulk,s = V_s + J_s/rho_m,
partial_t Sigma + div_x(Sigma v) = -(J_+ + J_-),
Sigma = rho_4D W,                        v = partial_t u,
delta_v Sigma_mat = 0,
delta_v theta + delta_v e_W + div_x(delta_v u) = 0  (uniform linearisation).
```

Obtain the equations of motion by the supplied method at spec `:152–154`: balance laws, the binding
virtual-displacement rule, variational derivatives with held-fixed fields named, and prescribed external
virtual work—not by placing an irreversible response kernel in an ordinary action.

The engine, not this directive, enumerates the field coordinates. SymPy declares the base and wave-jet
coordinates at `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:143–169,214–294`, the inherited
face-centre coordinates at `:568–577`, the bulk trace/interior/generic coordinates at `:484–550`, and the
virtual/test coordinates at `:552–599`. Use those declarations and the constructed object's actual free symbols;
do not substitute an authored field list.

Keep the accepted energy object intact: 40 records = 10 uniform + 15 `W_BG` first-jet + 15 `MU_R_BG`
first-jet records (step record `research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md:35–37`).
The SymPy engine constructs each retained record and its coefficient slot at
`scripts/S11c_b_brane_operator_sympy_audit.py:1760–1837` and emits the retained basis at `:4119–4141`.
Before construction, each engine must emit an ordered manifest of the 40 engine records: exact record path,
coefficient slot, free-symbol list and defining file:line. That manifest is the coefficient domain. Keep every
manifested expression and free symbol live, plus both bookkeepers `eta_bg` and `sigma_W`, both density
representatives, both anchorings and every background jet admitted by R1. Retain first order in wave amplitude
and first order in each independent background bookkeeper. No numerical model point is part of this build.

This section supplies physics; the build does not test or revise it.

## 2 · Part 1 object and domains

### 2.1 Cases and index convention

Build every pair in the exact engine case product `BRANCHES × DENSITY_REPS`, whose members are declared at
`scripts/S11c_b_brane_operator_sympy_audit.py:54–55`.

Directions 1, 2 and 3 are engine symbol-name labels. SymPy declares their displacement and profile-jet symbols at
`:215–233,182–199`; Wolfram declares `spatialCoordinates={xOne,xTwo,xThree}` at
`research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl:239–242`. SymPy's
`DIRECTIONS=range(3)` at `:57` is a separate 0-based loop index. The report gives each engine's mapping with
file:line.

### 2.2 Background, trial and test domains

**R1 is a symmetry/domain restriction, applied before construction.** Every background datum—scalar, polar,
axial, tensor, support, face, bulk-background, domain and boundary datum—is fixed by the full `O(2)` of rotations
and reflections about direction 1, and every background profile depends on `x_1` alone. For profile jets only,
set a jet to zero when any spatial-derivative index has symbol-name label 2 or 3; retain every allowed direction-1
jet through the engine's carried order. Do not apply component-name zeroing to general vectors or tensors.
Implement the group-fixed domain in each engine before the energy, background-jet maps, face substrate or
variations are constructed. The report states how, with file:line.

**Class P trial domain.** Physical perturbations are independent of `x_3`, have `x_2,t` dependence
`exp(i(k_2 x_2-omega t))`, and retain general `x_1` dependence. `k_2` and `omega` remain symbolic and live.

**Class P test domain and pairing.** Test fields carry the conjugate factor
`exp(-i(k_2 x_2-omega t))`. Trial and test envelopes have compact support in `x_1` only. Define the weak pairing
per one `x_2` period and per unit `x_3` length (with the `k_2→0` member understood by the continuous Fourier
normalisation). This is the class-P restriction of the supplied weak pairing at spec `:312–346`; it does not
claim compact support in directions in which the Fourier mode is extended.

### 2.3 Code-defined codomain

For each case, reconstruct the complete `S11CB_SLAB_OPERATOR` and `S11CB_MU_THETA_OPERATOR` emitted by the engine
at `:4156–4183`, from the whitelisted inputs in §4. The SymPy operator's code-defined top-level paths are exactly:

```text
U_BODY_BALANCE
THETA_BALANCE
E_W_BALANCE
ADVECTIVE_MASS_OPERAND
FACE_FLUX_BOUNDARY_OPERANDS
FACE_GENERALIZED_FORCE_ROWS
MU_THETA_FACE_BINDING
```

The first four are created at engine `:2416–2433`; the face and binding paths are attached at `:3103–3127`.
`FACE_FLUX_BOUNDARY_OPERANDS` is the entire code-defined substrate bundle, not an authored face-quantity list.
Its exact keys, enumerated by the engine at `:370–390`, are:

```text
background_density_map
face_normal
conormal_deriv
face_measure_shape_deriv
face_velocity
relative_flux
kinematic_balance
traction
virtual_work_shape_deriv
face_shift
projection_shape_deriv
projection_term_origins
projection_static_operand
projection_dynamic_operand
projection_residual
virtual_constraint
evolution_mass_balance
evolution_term_origins
closure_shape_deriv
```

The bundle is filtered and attached by the engine at `:2024–2075,2834–2897`. The face-source physical members
are declared in S11c-a at
`research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py:600–626`; its geometry, trace and law
objects are constructed at `:818–955,1014–1088`. These code containers—not a prose alias list—define the output.

Define `CODE_DEFINED_CODOMAIN` as the ordered scalar-expression leaves of those two reconstructed root objects,
after R1 and class P, retaining each exact tuple/dictionary path as its output label. Emit a manifest containing
every root/path, component position, declaring constructor and file:line. Print every leaf, including a leaf
computed as zero. Output membership is determined only by the code path.

### 2.4 Free-coordinate closure and formal map

For each case define, as an object,

```text
FREE_COORDINATE_CLOSURE :=
  {q in union_free_symbols(CODE_DEFINED_CODOMAIN) : DECLARED_SYMBOLS[q].class == COORDINATE}.
```

Use the engine registry and declaration sites, not a name allowlist. Emit every member with its exact symbol,
engine class/description/dimension, declaring file:line and all codomain paths in which it occurs. Type it by its
declaration source as one of:

```text
SLAB_TRIAL
BULK_FACE_TRACE
BULK_INTERIOR_OR_GENERIC
VIRTUAL_TEST
INDEPENDENT_OR_SPECTRAL
```

The declaration spans are engine `:143–169,214–294,568–577` (slab), `:484–524` (face trace), `:526–550` (bulk
interior/generic), `:552–599` (virtual/test), and S11c-a engine `:76–100` (independent coordinates and `omega`).
The new ansatz symbol `k_2` is `INDEPENDENT_OR_SPECTRAL`. A coordinate in the closure is either a map column or is
listed in `HELD_COORDINATE_MANIFEST` with a reason. Independent coordinates and spectral parameters are held;
every perturbation field/value/jet and every virtual/test field/value/jet is a column. Report the set identity
between the closure and the disjoint union of map columns and held coordinates. Never silently hold a coordinate.

Define and print

```text
FORMAL_WEAK_COORDINATE_MAP :=
  D_(all non-held members of FREE_COORDINATE_CLOSURE) CODE_DEFINED_CODOMAIN |_(R1,P).
```

This full Jacobian deliberately includes typed virtual/test columns wherever the complete weak payload carries
them. It is a formal coordinate map, not a claim that every formal bulk coordinate is an independent physical
bulk degree of freedom.

### 2.5 Supplied rest-frame potential-flow pullback

The equations in this subsection are **supplied to the build and unfalsifiable within it**. They are consequences
of the supplied rest-frame acoustics and mass law (`S11b_SHARED_PHYSICS.md:162–170,183–193`; S11c-b spec
`:95–97`) under class P. For face `s`, define the flat reference height

```text
w_s := s W_0/2,
Phi_s := Phi(x_1,w_s),
Psi_s := partial_w Phi(x_1,w_s),
phi := Phi(x_1,w) exp(i(k_2 x_2-omega t)).
```

The engine's perturbation trace coordinates are based at this flat reference face, exactly as S11c-a engine
`:150–180,576–585` states. Do not evaluate `Phi_s` or `Psi_s` at `s W_bg/2`; the shifted-trace construction moves
the flat-reference data to the background/physical face.

Supply, for each face:

```text
delta_p_s = i rho_m omega Phi_s,
delta_v_bulk_s = (partial_1 Phi_s, i k_2 Phi_s, 0, Psi_s),
partial_w delta_p_s = i rho_m omega Psi_s,
partial_w delta_v_bulk_s =
  (partial_1 Psi_s, i k_2 Psi_s, 0,
   (k_2^2-omega^2/c_s0^2) Phi_s - partial_1^2 Phi_s),

delta_rho_s = i rho_m omega Phi_s/c_s0^2,
partial_w delta_rho_s = i rho_m omega Psi_s/c_s0^2,
delta_j_s = rho_m delta_v_bulk_s,
partial_w delta_j_s = rho_m partial_w delta_v_bulk_s.
```

For the projection's bulk-interior density-time coordinate, supply

```text
partial_t delta_rho(x_1,w) = rho_m omega^2 Phi(x_1,w)/c_s0^2.
```

The engine's `delta_j_bulk_1..4` values are deliberately unlabelled by face while their normal jets are
face-labelled (engine `:519–531`; S11c-a `:174–181`). Preserve that ownership: in a face-labelled output path,
evaluate the single bulk-current field at that path's `w_s`; in a projection path, retain it as the bulk field
`rho_m grad_4 phi`. Do not mint independent plus/minus value coordinates. Emit the path-to-evaluation assignment.

`trace_grad_f_1..4` and `d_w_trace_grad_f_1..4` are engine-declared generic conormal operands (`:533–543`;
S11c-a `:946–954`). No governing source identifies their generic `f` with `phi`; keep them formal in the pulled-back
payload and list them in `HELD_GENERIC_TRACE_MANIFEST`.

Define `POTENTIAL_FLOW_TRACE_MAP` by applying these supplied substitutions, with the path ownership above, to
`FORMAL_WEAK_COORDINATE_MAP`. Emit every entry and both held manifests. The equations themselves are not an
acceptance test.

### 2.6 Separate virtual/test kinematic map

For the same code-defined face cases, take the scalar leaves of `FaceSource.virtual_displacement`, declared at
S11c-a engine `:600–626` and constructed at `:780–861`. Its domain is the subset of the computed
`FREE_COORDINATE_CLOSURE` typed `VIRTUAL_TEST` that those leaves actually contain. Define and print

```text
VIRTUAL_KINEMATIC_MAP :=
  D_(its computed VIRTUAL_TEST domain) (FaceSource.virtual_displacement leaves) |_(R1,P).
```

Emit its domain/codomain manifest and every labelled entry. It is a separate kinematic object even though the
complete weak-coordinate map also types every test coordinate it encounters.

### 2.7 Entry representation

Represent both engines' maps as differential-operator entries after
`partial_2→i k_2`, `partial_3→0`, and `partial_t→-i omega`, while retaining `x_1` jets explicitly. One record is
keyed by

```text
(control, map, branch, density representative,
 exact output root/path/component, input field family, x_1-jet order).
```

The engine field/jet manifest supplies the input-family and jet correspondence. Do not compare pretty-printed
whole rows or infer a column from a result's value. Emit every record, including computed zeros.

## 3 · Supplied FORM controls

Apply R1 and class P to the baseline first. The following are **supplied controls, unfalsifiable within this
build**; their output entries are not specified here.

Append K1 and K2 separately to the stored-energy density and rerun the complete §2 construction:

```text
K1: a_K1 e_hat_i (partial_k u_i)(partial_k theta),
    e_hat = (sin(beta), 0, cos(beta));

K2: a_K2 theta epsilon_ijk g_i partial_j u_k,
    g_i = partial_(y_i) W_bg.
```

The control-only datum `e_hat` is introduced after R1 and is not a background datum in R1's domain. Indices
`i,j,k` run over engine labels 1,2,3; `epsilon_ijk` is the Levi-Civita symbol. Keep `a_K1`, `a_K2` and `beta` live,
with coefficient dimensions required for stored-energy density.

For G3, introduce after R1 a fixed, untransformed direction-3 background first-jet datum:

```text
G3: g_i^(G3) = g_i + a_G3 delta_(i3),
    with a_G3 delta_(i3) held unchanged when the trial/test fields are reflected.
```

Keep `a_G3` live with the same dimension as `g_i`. Unlike K1/K2, G3 is a controlled change of the background
domain. Feed it through every construction that consumes the thickness-background first jet: stored energy,
face geometry, shifted traces, projection, weak rows and kinematics. Do not insert it at a finished row or map.

For K1, K2 and G3, print all three complete maps under tags naming the control and map, then print the computed
entrywise difference from the baseline. If two complete payloads are byte-identical, report that comparison after
emission. No control result is an acceptance value.

## 4 · Engines and comparator

**SymPy import whitelist.** From `scripts/S11c_b_exports.py` or the existing SymPy engine, import only symbol
definitions, registry metadata and the accepted 40-record stored-energy basis. Do not import any operator row,
face row, coupling kernel or other derived payload. Represent R1 in the imported symbols/basis before construction.

**Blind Wolfram engine.** Import nothing. Re-derive the supplied physics and code-defined object schema from the
governing specs and cited engine declarations; do not transcribe a Python result. Apply §0's guarded startup probe
before its first construction.

**Raw-first comparator.** Join the jet-indexed records by the §2.7 key. For baseline, K1, K2 and G3, emit first
the raw `operand_PY`, raw `operand_WL` and raw `residual` for every record, with no coefficient, sign or row
convention applied.

Only after the raw payload, emit separately tagged mapped diagnostics:

1. The P2b coefficient-normalisation convention is the committed map at
   `research/pde_ledger_v3/directives/S11c_b_p2b_gamma_bridge_directive.md:10–17,29–57`:
   `gamma_WL = W_0 gamma_PY` for the `W_BG` family and `gamma_WL = mu_R gamma_PY` for the `MU_R_BG` family, paired
   by the unfolded invariant object rather than emission position. Emit the applied expression-valued mapping and
   its record pairing. P2b's production build remains recorded as deferred work (step record `:54–56,106–115`);
   do not invent a different bridge.
2. The kinetic and face-generalized-force whole-row signs are separate mapped diagnostics. The step record requires
   the raw comparator to surface them, not normalize them (`:112–114`). State every applied map and source
   file:line. A mapped diagnostic never replaces raw output.

If an engine or comparator cannot complete under §0, report the contained result and stop that job. Do not
substitute another computation.

## 5 · The three clauses and structural rule

> **1. The script may PRINT computed objects. It may NOT state conclusions.** An `emit`/`Print` payload is a CAS
> object—an expression, a solved root or a boolean from a symbolic test—not prose describing a result.
>
> **2. PRINT the residual; do NOT assert it.** Compute → emit → then assert.
>
> **3. Interpretation belongs to the STEP RECORD.** The scripts do not editorialise.

> **The only places physical symbols may be combined by hand are the supplied stored energy, supplied face laws,
> supplied potential-flow pullback, ansatz/domain, and the supplied controls in §3. K1/K2 enter at stored energy;
> G3 enters at the background domain. Every other expression involving physical symbols is reached by computation.
> No control is inserted at a finished row, block, map or result.**

Tag names name objects—field, coordinate, row, map or control—never a value, sign or shape of a result. Entry values
are not acceptance criteria. Printing every entry, including zeros, is the object; no entry is selected for
emission by inspecting its value.

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
   retained order, with file:line, including the engine's carried mixed-jet convention;
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
  - the code-defined codomain and free-coordinate-closure manifests, including every typed or held coordinate;
  - the R1 and class-P implementation in each engine, with file:line;
  - the flat-reference-face and unlabelled-current path assignments;
  - every comparator mapping and its source file:line;
  - any byte-identical control payloads and contained failures;
  - no interpretation of map entries.
- `research/pde_ledger_v3/_measurements/S11c_d_round_background_scope.md` — the Part 2 scope note.

Stop after writing these.
