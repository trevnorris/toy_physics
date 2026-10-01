# Build directive — S11c-b slab/face operator blocks on a one-direction background (Part 1) + round-background scope note (Part 2)

You are the builder. Your job is: **build → run → report → stop.**
- ⛔ Do not launch, call or spawn any other AI, agent, reviewer or review process.
- ⛔ Do not commit.
- Write only the files named under *Deliverables*.

## 0 · Execution safeguards

1. **Guarded runs only.** Run every engine, constructor and comparator job under `scripts/s11c_guarded_run.py`,
   as `AGENTS.md:8–12` requires: a whole-job 2 GiB memory cap, one job at a time, no overlapping CAS, and no
   unguarded fallback. A contained failure is a **result**, whether the cap is hit or containment is
   unavailable. Record the measured peak and the stage reached, then stop that job. ⛔ Never relaunch
   unguarded, and ⛔ never raise the cap.
2. **The suspended job.** `research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md:5`
   says no other job may run alongside the suspended S11c-d process (PID 4097233). Run `ps -p 4097233`. If the
   process exists in any state:
   - do the source reading;
   - write the Part 2 scope note, without its probe;
   - stop before running any engine, and report that you stopped for this reason.
3. **Wolfram.** Run at most one kernel at a time, and print observable progress.

## 1 · Governing physics — supplied, and unfalsifiable within this build

**The operator.** It is the S11c-b variable-coefficient slab operator with its face-response coupling, exactly
as governed by `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`, which amends
`S11c_a_SHARED_PHYSICS.md` and `S11b_SHARED_PHYSICS.md`:
- the slab degrees of freedom `u` (three in-plane components, no `w`-component), `θ`, and the two independent
  face variables `ζ_+`, `ζ_-` (spec `:59–71`);
- the face laws and face responses `Λ_A`, `Λ_V`, `Λ_X` (spec `:145–148`);
- **all 40 accepted stored-energy basis terms**, each with its own symbolic coefficient (record `:35–37`);
- both bookkeepers `η` and `σ_W`; both density representatives; both anchorings, `LAB_HELD` and
  `MATERIAL_ADVECTED`;
- the retained order: first in wave amplitude `ε`, first in each of `η` and `σ_W`.

**Method.** Equations of motion are obtained by the spec's method (spec `:151–154`): "balance laws, the binding
virtual-displacement rule, variational derivatives with held-fixed fields named, and prescribed external
virtual work — **not** by putting an irreversible response kernel in an ordinary action."

**Sources.** The step record is `research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md`. The
engines are `scripts/S11c_b_brane_operator_sympy_audit.py` and
`mathematica/S11c_b_brane_operator_mathematica_audit.wl`. The exports are `scripts/S11c_b_exports.py`
(`af560257`).

**Status.** All of the above is supplied:
- per-engine verified;
- cross-engine residual **deferred**;
- the two whole-row sign conventions are cross-engine-**unvalidated** (record `:15–31`).

This build does not test any of it.

**Declared freezes.** State both in the first lines of the report:
1. `v_bulk_normal_0`, the bulk drain, "appears in no derived operator" (spec `:90–91`). Part 1 has **no drain
   flow**.
2. S11c-b performs no bulk response solve. The bulk enters **only** as face operands.

## 2 · Part 1 — the object

**Index convention.** Directions 1, 2 and 3 are the engines' **symbol-name labels**:
- in SymPy, the direction suffix in symbol names such as `u_1` and `w1_profile_d1`
  (`scripts/S11c_b_brane_operator_sympy_audit.py:190–195`, `:420`);
- in Wolfram, `xOne`, `xTwo`, `xThree`.

⚠ The SymPy loop variable `DIRECTIONS = range(3)` (`:57`) is 0-based and is **not** this label. State the mapping
you used in each engine, with file:line.

**Background restriction (R1), applied before construction.**
1. Every background symbol whose name carries direction label 2 or 3 in **any** position is zero. That includes
   mixed jets such as one carrying labels 1 and 2, and every higher jet through the highest order the engines
   carry.
2. Every other background datum depends only on direction 1: density representatives, support bundle, face,
   bulk-background, domain and boundary data.
3. Every background **vector** datum has zero components along directions 2 and 3. This includes the hold force
   `f_hold⁰`, the hold tractions `t_hold,s⁰` and the boundary loads (spec `:195–196`).
4. Every background tensor datum is invariant under rotations about direction 1.

State how R1 is imposed in each engine, with file:line. ⛔ Do not build the unrestricted operator and restrict it
afterwards.

**Perturbation class `P`.** Every perturbation field, in whichever coordinates each engine uses, is independent
of direction 3. Dependence on direction 2 is the form `exp(i(k₂ x₂ − ω t))`, with `k₂` and `ω` symbolic and
live. Dependence on direction 1 is general.

**Object.** One rectangular linear map: the first derivative, at the background and on class `P`, of an **ordered
output vector** with respect to an **ordered input vector**. Name every input and output, with the file:line of its
definition, and print the full matrix.

- **Inputs, in order.**
  - The slab fields `u_1`, `u_2`, `u_3`, `θ`, `ζ_+`, `ζ_-`.
  - For each face `s = ±`, the bulk trace inputs: `δp_s` and the four components of `v_bulk,s`. Without a bulk
    solve these are independent inputs, not outputs.
- **Outputs, in order.**
  1. Every operator row as the engines define it: after the constraint fold (pin B), with `μ_θ` kept as its named
     operand. This includes:
     - the three `U` body-balance rows;
     - the `θ` and `e_W` balance rows;
     - the face generalized-force rows, including `CENTER_FACE_GENERALIZED_ROW` for the independent `ζ_c`
       (`scripts/S11c_b_brane_operator_sympy_audit.py:2237`);
     - the `μ_θ` face binding.
  2. For each face `s = ±`, every face quantity the face laws use (S11c-a `:343–354`, `:365–366`):
     - the components of `n̂_s`;
     - the components of `v_face,s`;
     - the components of `δ_v x_s`, as functions of the virtual fields;
     - `V_s`, `a_s`, `J_s`, `𝒜_s`, `μ_s`;
     - the components of `t_s`.

**Form.** You may use the strong form, or the weak form with trial and test fields in `P` (spec `§3c`). State which.
If you use the weak form, the test-field list mirrors the slab inputs.

⭐ **Construct the operator for all fields at once** from the stored energy and the supplied face laws, by the
method in §1. Then read the blocks off it. ⛔ Do not construct any single field's equation on its own, and ⛔ do
not split the fields into groups before the operator is constructed.

Print **every entry, including any that evaluates to 0, as a computed object**. Emission must not depend on any
entry's value: which entries appear is decided only by the field and row labels.

**Model point.** Keep every physical symbol live: the 40 basis coefficients, `μ_R`, the density representatives,
`c_s0`, `ρ_m`, `B_ρ`, `C`, `k_W`, `κ_W`, `μ_W`, `Λ_A⁰`, `Λ_V⁰`, `Λ_X⁰`, `τ_A`, `τ_V`, `τ_X`, `k₂`, `ω`, every
background jet that survives R1, and the bookkeepers `ε`, `η`, `σ_W`. ⛔ No numeric substitution in constructing
the operator. If an entry is too large to print in closed form:
1. print its exact carrier structure, plus exact rational evaluations at declared witness points;
2. declare every witness point, and state that witnesses are not model values;
3. choose one witness per component of any union-of-loci the carriers define.

## 3 · Part 1 — controls

These are FORM controls. Each one is a term **added to the stored-energy density** with a live coefficient. Then
the whole Part 1 construction is re-run. The term itself is the only hand-written addition. Summation runs over
the in-plane indices `i, j, k ∈ {1, 2, 3}`; `ε_ijk` is the Levi-Civita symbol; `ê ≡ (sin β, 0, cos β)` in
components (1, 2, 3), with `β` live; `g_i ≡ ∂_{y_i} W_bg`.

```text
K1 :   a_K1 · ê_i (∂_k u_i)(∂_k θ)
K2 :   a_K2 · θ · ε_ijk g_i ∂_j u_k
```

`a_K1` and `a_K2` are live symbols carrying whatever units make each term a stored-energy density. For each
control, print the complete block matrix exactly as for the baseline, plus the entry-by-entry difference from
the baseline. If a control's output is byte-identical to the baseline, report that explicitly.

## 4 · Engines and comparator

**SymPy.** From `scripts/S11c_b_exports.py` or the S11c-b SymPy engine, it may import **only**:
- symbol definitions;
- the accepted 40-term stored-energy basis.

⛔ It may not import any exported operator row, coupling kernel, face row, or other derived payload. R1 enters the
energy, the background-jet maps and the face substrate **before** any variation.

**Wolfram.** It is written blind: it **imports nothing**, re-derives from the governing specs, and ⛔ is never a
transcription of the `.py`.

**Comparator.** Join the two engines' printed matrices by input and output label. For the baseline and for each
control, emit **first**, for every entry, the raw `operand_PY`, the raw `operand_WL` and the raw `residual`, with no
convention applied. The record requires the comparator to surface the two whole-row sign conventions, not
normalize them (record `:112–114`): the kinetic-term sign, and the face generalized-force convention. A
convention-mapped diagnostic may follow, under distinct tag names, with each map stated alongside the file:line it
comes from. ⛔ The mapped diagnostic never replaces the raw comparison.

If an engine or the comparator cannot complete under §0, report the measurement and stop that job. ⛔ Do not
replace it with something else.

## 5 · The three clauses, and the structural rule

> **1. The script may PRINT computed objects. It may NOT state conclusions.** An `emit`/`Print` payload must be
> a CAS object — an expression, a solved root, a boolean from a symbolic test. ⛔ Never prose describing a
> result.
> **2. PRINT the residual; do NOT assert it.** Compute → emit → *then* assert.
> **3. Interpretation belongs to the STEP RECORD.** ⛔ The script does not editorialise.

> **The ONLY place the physical symbols may be combined by hand is in constructing the STORED ENERGY (including
> the K1/K2 terms of §3), the SUPPLIED FACE LAWS, and the ANSATZ. Every other expression involving them must be
> REACHED BY COMPUTATION. Every control re-enters the chain at the stored energy, ⛔ never at a result.**

Tag names name the object (field, row, block, control), ⛔ never a value, sign or shape of a result.

## 6 · Part 2 — scope note only (⛔ do not build)

**The later object.** The complete linear operator on a background invariant under all rotations and reflections
of the three in-plane coordinates about a point (`O(3)`). The bulk drain flow is carried as a **live background
field**: a normal component at the faces plus a radial in-plane component, with no azimuthal component. The
blocks of interest are those between the toroidal displacement fields and every other field, at each angular
order. Toroidal fields are tangent to the spheres `r = const` and divergence-free. The bulk response is closed
rather than left as face operands.

**Write a scope note** answering:
1. For each ingredient, which committed spec governs it (file:line) and which is **missing**. The ingredients are:
   - the slab energy;
   - the face laws and face responses;
   - the bulk closure (including the S11c-c1 curved-bulk closure);
   - the drain flow.

   In particular: how the drain flow enters the slab/face/bulk equations, given spec `:90–91`.
2. Whether the existing three-direction background jets can represent an `O(3)`-symmetric profile at the retained
   order, and what changes if not (file:line).
3. A measured memory/time probe on a small case, under §0, only if §0 permits it.
4. What cannot be done under the §0 guard.

Then **stop**.

## 7 · Deliverables

- `research/pde_ledger_v3/scripts/S11c_d_zinvariant_operator_blocks_sympy.py`, with literal stdout at
  `research/pde_ledger_v3/scripts/out/S11c_d_zinvariant_operator_blocks_sympy.out`
- `research/pde_ledger_v3/mathematica/S11c_d_zinvariant_operator_blocks.wl`, with literal stdout at
  `research/pde_ledger_v3/mathematica/out/S11c_d_zinvariant_operator_blocks.out`
- `research/pde_ledger_v3/scripts/S11c_d_zinvariant_operator_blocks_comparator.py`, with literal stdout at
  `research/pde_ledger_v3/scripts/out/S11c_d_zinvariant_operator_blocks_comparator.out`
- `research/pde_ledger_v3/_measurements/S11c_d_zinvariant_operator_blocks_report.md`. It contains:
  - the declared freezes, first;
  - each command run, with its guard invocation, exit code, wall time and measured peak memory;
  - the output paths;
  - the index convention, the R1 implementation, the field list and row definitions, with file:line;
  - any convention maps used by the mapped diagnostic;
  - any byte-identical controls and any contained failures.

  ⛔ It contains no interpretation of the entries.
- `research/pde_ledger_v3/_measurements/S11c_d_round_background_scope.md` — the Part 2 scope note.

Stop after writing these.
