# Build directive — S11c-b operator blocks on z-independent perturbations (Part 1) + round-background scope note (Part 2)

You are the builder. Your job is: **build → run → report → stop.**
- ⛔ Do not launch, call or spawn any other AI, agent, reviewer or review process.
- ⛔ Do not commit.
- Write only the files named under *Deliverables*.

## 0 · Before any engine run

`research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md:5` says no other
job may be launched alongside the suspended S11c-d central-balance process (PID 4097233). Run `ps -p 4097233`.
If that process exists in any state:
1. do the source reading for Parts 1 and 2;
2. write the Part 2 scope note, without its memory probe;
3. stop before running any engine. Report that you stopped for this reason.

## 1 · Governing physics — supplied, and unfalsifiable within this build

**The operator.** It is the S11c-b variable-coefficient slab operator with its face-response coupling, exactly
as governed by `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`, which amends
`S11c_a_SHARED_PHYSICS.md` and `S11b_SHARED_PHYSICS.md`. Specifically:
- the slab degrees of freedom are `u` (three in-plane components, no `w`-component), `θ`, and the two independent
  face variables `ζ_+`, `ζ_-` (spec `:59–71`);
- the face responses are `Λ_A𝒜_s`, `Λ_V V_s` and `Λ_X𝒜_s`;
- the background ansatz is `W_bg(y)` and `μ_R,bg(y)`, with the density representatives `ρ_4D,bg⁰`/`ρ_br,bg⁰`;
- both anchorings are included: `LAB_HELD` and `MATERIAL_ADVECTED`;
- the retained order is first in wave amplitude `ε` and first in each of `η` and `σ_W` (spec `§2a`).

The step record is `research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md`. The engines are:
- `scripts/S11c_b_brane_operator_sympy_audit.py`;
- `mathematica/S11c_b_brane_operator_mathematica_audit.wl`.

The exports are `scripts/S11c_b_exports.py` (`af560257`).

**Status of what you are handed.** All of the above is supplied:
- per-engine verified;
- cross-engine residual **deferred**;
- the two whole-row sign conventions are cross-engine-**unvalidated** (record `:15–31`).

This build does not test any of it.

**Declared freeze.** `v_bulk_normal_0`, the bulk drain, "appears in no derived operator" (spec `:90–91`). Part 1
therefore runs with **no drain flow**. State this in the first lines of the report.

**Coordinates.** `y` is the in-plane coordinate along which the background varies. `x` and `z` are the two
in-plane coordinates along the interface. `w` is the slab-normal coordinate.

## 2 · Part 1 — the object

**Perturbation class `P`.** Every perturbation field is independent of `z`. Dependence on `(x, y, t)` is general,
or the form `exp(i(k_x x − ω t))` × functions of `y`, with `k_x` and `ω` symbolic and live.

**Object.** The complete linear operator of the S11c-b system acting on class `P`, as a matrix over the full field
list the S11c-b engines carry. Name each field and each row, with the file:line of its definition. Use the rows as
the S11c-b engines define them: after the constraint fold (pin B), with `μ_θ` kept as its named operand, and the
face generalized-force rows including the face responses. Keep the bulk face loads (`δp_±`) as symbolic
operands, as S11c-b does. You may present it either as:
- the strong-form operator matrix, or
- the weak-form bilinear form restricted to trial and test fields in `P` (spec `§3c`).

State which, and print **every block**.

⭐ **Construct the operator for all fields at once** from the governing action and face laws, then read the
blocks off it. ⛔ Do not construct any single field's equation on its own, and ⛔ do not split the fields into
groups before the operator is constructed.

Print **every entry, including any that evaluates to 0, as a computed object**. Emission must not depend on any
entry's value: which entries appear is decided only by the field and row labels.

**Model point.** Keep every physical parameter symbolic: `μ_R`, the density representatives, `c_s0`, `ρ_m`,
`B_ρ`, `C`, `k_W`, `κ_W`, `μ_W`, the `Λ`s and their relaxation times, `k_x`, `ω`, and the background jets. ⛔ No
numeric substitution in constructing the operator. If an entry is too large to print in closed form:
1. print its exact carrier structure, plus exact rational evaluations at declared witness points;
2. declare every witness point, and state that witnesses are not model values;
3. choose one witness per component of any union-of-loci the carriers define.

## 3 · Part 1 — controls

These are FORM controls. Each one re-enters at the action, ⛔ never at a result. For each control, add the term
to the stored energy with a **live coefficient**, state the exact term you added, and re-run the whole Part 1
construction.

- **K1 — anisotropy axis.** Add one term that couples field gradients to a fixed in-plane unit vector
  `ê = (0, sin β, cos β)` in `(x, y, z)` components, with `β` live.
- **K2 — parity-odd term.** Add one term containing a single Levi-Civita contraction, built from the fields and
  the background gradient.

For each control, print the complete block matrix exactly as for the baseline, plus the entry-by-entry difference
from the baseline. If a control's output is byte-identical to the baseline, report that explicitly.

## 4 · Engines

**SymPy.** It may import `scripts/S11c_b_exports.py`.

**Wolfram.** It is written blind: it **imports nothing**, re-derives from the governing specs, and ⛔ is never a
transcription of the `.py`.
- Run at most one kernel at a time; the licence has two seats.
- Print observable progress.
- If the kernel's resident memory passes 20 GB, stop it and report the measured peak and the stage reached. ⛔ Do
  not substitute a reduced variant without saying so.

If one engine cannot complete on this machine, report the measurement and stop that engine. ⛔ Do not replace it
with something else.

## 5 · The three clauses, and the structural rule (verbatim)

> **1. The script may PRINT computed objects. It may NOT state conclusions.** An `emit`/`Print` payload must be
> a CAS object — an expression, a solved root, a boolean from a symbolic test. ⛔ Never prose describing a
> result.
> **2. PRINT the residual; do NOT assert it.** Compute → emit → *then* assert.
> **3. Interpretation belongs to the STEP RECORD.** ⛔ The script does not editorialise.

> **The ONLY place the physical symbols may be combined by hand is in CONSTRUCTING THE ACTION and the ANSATZ.
> Every other expression involving them must be REACHED BY COMPUTATION. Every control re-enters the chain at
> the ACTION, ⛔ never at a result.**

Tag names name the object (field, row, block, control), ⛔ never a value, sign or shape of a result.

## 6 · Part 2 — scope note only (⛔ do not build)

**The later object.** The complete linear operator on a background invariant under all rotations and reflections
of the three in-plane coordinates about a point (`O(3)`). The bulk drain flow is carried as a **live background
field**: a normal component at the faces plus a radial in-plane component, with no azimuthal component. The blocks
of interest are those between the toroidal displacement fields and every other field, at each angular order.
Toroidal fields are tangent to the spheres `r = const` and divergence-free.

**Write a scope note** answering:
1. For each ingredient (slab energy, face laws, face responses, bulk closure, drain flow), which committed spec
   governs it, with file:line, and which is **missing**. In particular: how the drain flow enters the
   slab/face/bulk equations, given spec `:90–91`.
2. Whether the existing engines' basis and rows are written for a general background-gradient direction, or only
   for variation along `y`. Give file:line.
3. A measured memory/time probe on a small case, only if §0 permits it.
4. What cannot be done on this 30 GB machine.

Then **stop**.

## 7 · Deliverables

- `research/pde_ledger_v3/scripts/S11c_d_zinvariant_operator_blocks_sympy.py`, with literal stdout at
  `research/pde_ledger_v3/scripts/out/S11c_d_zinvariant_operator_blocks_sympy.out`
- `research/pde_ledger_v3/mathematica/S11c_d_zinvariant_operator_blocks.wl`, with literal stdout at
  `research/pde_ledger_v3/mathematica/out/S11c_d_zinvariant_operator_blocks.out`
- `research/pde_ledger_v3/_measurements/S11c_d_zinvariant_operator_blocks_report.md`. It contains:
  - the declared freeze, first;
  - each command run, with exit code, wall time and peak memory;
  - the output paths;
  - the field list and row definitions, with file:line;
  - the exact K1/K2 terms added;
  - any byte-identical controls.

  ⛔ It contains no interpretation of the entries.
- `research/pde_ledger_v3/_measurements/S11c_d_round_background_scope.md` — the Part 2 scope note.

Stop after writing these.
