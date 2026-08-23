# Independent physics review (2nd pass) — the S11c-a spec (background & geometry)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`

⚠ **This spec was revised once** after a prior review. ⛔ **Do not assume the revision is correct.** A fold
can over-correct, double-count, introduce an internal inconsistency, or fail to actually close what it
claims to. Verify the physics is now right **and** hunt any defect the fold introduced — derive your own
view and let it disagree with the artifact.

## What it is, and why it is load-bearing
S11c is the non-uniform transverse coupling at a brane–bulk interface — *is light's confinement
unconditional?* The prior step **S11b** (uniform) is CLOSED and proved the uniform transverse coupling is
**identically zero**; S11c computes the gradient-driven channel. S11c is a five-sub-step family; **S11c-a**
is the first: it **un-freezes the in-plane background** (`W₀→W₀(x)`, and the twist modulus `μ_R→μ_R(x)`) and
computes the **first-order geometry** — the tilted-face kinematics and the shape-derivative of every
interface law S11b wrote for a flat wall. This spec is orchestrator-written and is the physics authority two
blind engines (SymPy + Wolfram) will read; an error here corrupts both. ⛔ This is a **physics review**, not
a copy-edit.

## Derive your own view FIRST, then read the spec
Before opening the spec, form your own answer to: *what must a correct "background & geometry" sub-step set
up, given S11b's flat-wall interface laws and the folded S11c decision list?* Read:
- `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` — the closed uniform spec S11c-a builds on
  (esp. §1 geometry/fields/conventions `:63-160`, §3 relative flux `:183-207`, §4 face closure/affinity
  `:208-254`, §5 energy `:255-300`, §6 balance laws + virtual constraint + sign conventions + exit semantics
  `:301-380,:512-546`, B0 bulk face response + the physical mass balance `:593-646`, B6 transverse
  `:819-831`).
- `research/pde_ledger_v3/directives/S11c_decisions.md` — the folded decision list; S11c-a must honor
  **N3** (the FULL tilted-face shape-derivative, ⛔ not only normal face motion), **N4** (Eulerian/material
  is a representation S11b already fixed; the physics is profile *anchoring*), **N11b** (the frozen-wall
  reconciliation), **N12** (two-parameter ε,η power counting + background admissibility), **N14** (F9 name
  reservations), **N15** (the variable-coefficient operator + new invariants are S11c-b's), **N6** (the η→0
  limit is a known-vacuous regression, not the main control).
- `research/pde_ledger_v3/V3_STEP_PLAN.md:1179` — the coupling is `∝ k·a`, `∇μ_R≠0` mixes; **form** now with
  a generic profile, **magnitude** needs the interior (`R1`).
- `research/pde_ledger_v3/steps/S11b_interface_coupling_law.md` — the closed uniform result.
- `research/pde_ledger_v3/directives/_measurements/S11c_a_SHARED_PHYSICS.md` — the rule-2 twin; spot-check
  that §1's claims about the imported S11b setup resolve.
There is **no do-not-read list** — read anything that helps you judge the physics.

## Required method — DOCUMENT branch
Read the sources first, form your own view, then read the spec. Quote **both sides** (source and spec line)
for every finding. ⚠ Reading order is a method instruction, not a blindness control.

## What to check — hunt each, and hunt fold-induced defects
1. **Faithful supply of S11b (§1).** Are the imported objects quoted correctly — relative flux, closure,
   affinity `𝒜=μ_s−δp/ρ_m`, the **material** virtual constraint `δ_vθ+δ_ve_W+∇·δ_vu=0` (and that
   `δ_vΣ_E=0` is forbidden), the **physical** evolution mass balance `∂_tΣ+∇·(Σv)=−(J₊+J₋)`, the flat
   virtual work `δ_v𝒲_bulk`, the sign conventions, exit semantics? Any misquote corrupts both engines.
2. **Face geometry & orientation (§3).** Is the outward-normal orientation now unambiguous for **both**
   faces (⛔ not inward for the lower face), for **both** anchoring branches? Is the lower-face parity right?
3. **Anchoring (§2c) as two named branches.** Are lab-held and material-advected each supplied as a
   computable face map, so the `O(εη)` difference is a computed object rather than an engine choice?
4. **The virtual vs physical variations (T-g, T-h).** Does T-g linearize `δ_vΣ_mat=0` with the **virtual**
   `δ_v u` (⛔ not physical `u`), dimensionally consistently (⛔ no `Σ`-dimension term added to the
   dimensionless constraint), and does it capture the `e_W` redefinition that survives even when `∇Σ_E⁰=0`?
   Is the physical evolution balance (T-h) genuinely separate and present?
5. **Completeness of the shape-derivatives (T-a…T-i).** Is any interface law whose shape-derivative changes
   at `O(η)`/`O(εη)` still **missing** — the conormal/normal-derivative operator (T-a′), the tilted
   kinematic mass balance (T-c′), the in-plane virtual work (T-d), the assembled-closure + affinity-
   normalization term (T-i), the dynamic-window `∂_tΩ` (T-f)? Conversely, does any item **pre-state** its
   result form (violating rule 3)?
6. **μ_R(x) (2a).** Is it now actually introduced (varies, at order η, with anchoring), not left optional?
7. **Admissibility (2d).** Is the obligation now a computable, can-fail object (a named background state +
   declared self-supporting/externally-held support), rather than an identically-zero residual — and is the
   S11c-a↔S11c-b seam (premise here, failing check with S11c-b's operator) sound?
8. **Power counting (2e).** Are contrast `η` and slope `|∇W_bg|` correctly distinguished, and is the
   `(ε,η,slope)` order binding on every emitted object?
9. **Scope & controls.** Does S11c-a defer S11c-b…e correctly (exporting geometric boundary operators, ⛔
   not solving DtN), without presupposing a global dispersion relation? Are C-1 (form, off-diagonal drop
   across 3 directions) and C-2 (independence) genuine and able to fail? Is the emit-the-residual (⛔ not
   assert-zero) discipline in §4 correct (rule 5)? Is the no-VERDICT / exit-0-on-physics-disagreement rule
   intact?

## Physics filter
Report a finding only if it catches a way the physics or the sub-step's set-up could be wrong — not
stylistic preference, not "wrong on a different input." ⛔ A leg returning "nothing survives" is weak
evidence: state what you checked, what you derived, and where you looked hardest.

## Output
1. Findings, most-serious first — each with the source quote, the spec line, and the concrete way the
   physics/set-up goes wrong if uncorrected.
2. Is any interface-law shape-derivative still **missing**? Any **new** defect the fold introduced?
3. Anything **under-specified** (prose where an equation is needed) or **over-specified** (a result form
   pre-stated)?
4. Does the scope correctly defer S11c-b…e, and is the S11c-a↔S11c-c geometric-boundary seam right?
