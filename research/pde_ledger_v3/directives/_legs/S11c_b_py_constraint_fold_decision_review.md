# Decision-list review — PY constraint-fold directive (spec-pin B follow-on)

You are one of two independent decision legs (the other is the sibling engine). You are reviewing an
ORCHESTRATOR-written builder decision list BEFORE any engine is changed. Your job is to find every way the
directive is wrong, incomplete, over-reaching, or leaks an expected value — because the builder trusts this list and
nothing downstream re-checks it. A leg that finds nothing is weak evidence; derive and quote, do not rubber-stamp.

## The artifact under review

`research/pde_ledger_v3/directives/S11c_b_py_constraint_fold_directive.md`

## What is SETTLED — do NOT re-litigate (it has its own two decision legs + a Codex verification)

The S11c-b slab U-momentum row is PINNED as object **(B)**: the CONSTRAINT-REDUCED in-plane momentum balance
carrying the reaction `−∇μ_θ`, NOT the raw held-fixed `δU/δu`. `θ`'s virtual displacement is not independent. WL is
already correct; PY is the engine that must change. You are NOT asked whether (B) is right — that is closed. You are
asked whether the directive correctly, completely, and cleanly implements (B) in the PY engine. Pin records (read
for context, do not re-argue): `directives/S11c_b_jet_depth_spec_pin_decision.md` (FOLDED VERDICT section) and
`directives/_measurements/S11c_b_strong_row_jet_depth_reconciliation.md`.

## What you are handed (read the sources and form your own view of the code + spec)

- The directive above.
- Spec: `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` — read §1c (`:113-148`, esp. the SUPPLIED
  constraint `δ_vθ + δ_ve_W + ∇_x·δ_vu = 0` at `:143`, the sourced mass balance at `:142`, and `μ_θ ≡ δU/δθ|held-fixed`
  at `:128,131`), §3a (`:242-270`), §3b (`:272-288`), and §5 controls (`:387-460`).
- Governing method: `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md:337,341-342,426` (the virtual
  constraint, "Do NOT vary U with θ held fixed", "the same multiplier supplies the in-plane restoring force and the
  thickness term", and the `−∇(δU/δθ)` convention check).
- The PY engine: `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py` — read
  `operator_from_density` (`:1968-2062`), `strong_rows` (`:2072-2086`), `committed_strong_rows` (`:2123-2175`, the
  FROZEN reference #88 uses), `build_operator` (`:2211-2337`), `STRONG_ROW_JET_DEPTH`/`COUPLING_JET_DEPTH` (`:59-60`),
  the coupling extraction near `:2569-2594` (`min(background_depth, STRONG_ROW_JET_DEPTH)`), and the mass-balance
  ingredients `evolution_mass_balance` (`:2814,2821`) and `ADVECTIVE_MASS_OPERAND` (`:2043,2060`).
- #88 (must not be contradicted): `directives/S11c_b_88_blast_radius_build_directive.md`,
  `directives/_measurements/S11c_b_88_blast_radius_result.md`.

## Decide, WITH the code and spec quoted for every claim

1. **Correct + complete + no over-reach.** Is the change set (1)-(4) in the directive the correct and COMPLETE
   realization of pin (B)? Specifically: do the U-row AND the thickness row receive the reaction from the SAME
   multiplier `μ_θ` (S11b:341-342)? Is the `θ`-row correctly re-identified as the sourced mass-evolution equation
   (§1c:142), assembled from the ingredients PY already has, rather than a new construction or a deletion? Does `μ_θ`
   correctly stay a SEPARATE constitutive operand (held-fixed, §1c:128), not double-counted by the vector reaction
   `∇μ_θ`? Is raising `STRONG_ROW_JET_DEPTH` 2→3 correct AND necessary (the raw EL is depth-invariant; the reaction
   is the order-3 object)? Name any RIPPLE the directive fails to call out — e.g. the coupling cascade seeing a
   depth-3 operator via `min(background_depth, STRONG_ROW_JET_DEPTH)`; the tower-depth control; the frozen
   `committed_strong_rows` reference used by #88; the §5c uniform-limit operand and `S11B_TRANSVERSE_DISPERSION`; the
   MATERIAL-pullback route; the `energy_origins` per-term loop. Flag over-reach (anything changed that pin (B) does
   not require).
2. **§3b amendment accuracy.** Is the amendment paragraph in change (1) accurate and consistent with §1c? In
   particular: does it correctly keep `μ_θ` held-fixed as the CONSTITUTIVE derivative (§1c:131 "θ may not be
   eliminated through a constraint before this derivative") while the U-ROW carries the reaction `∇μ_θ` — i.e. are
   these two genuinely different objects, or does the paragraph conflate them or contradict §1c? Quote both.
3. **Rule-5 cleanliness.** Does ANY clause — in the directive body or the acceptance — leak an expected value: a
   target row content, a jet-atom count, a specific background-jet atom, or a "match the Wolfram engine" exit
   condition that a builder could iterate toward? The convention SIGN `−∇μ_θ` is SUPPLIED S11b:426 physics stated in
   the spec paragraph and the acceptance does not gate on it — judge whether that is a legitimate supplied premise or
   a leak, and whether the directive's "derive, do not hard-code" instruction plus the raw-EL diagnostic actually
   prevent a hand-typed reaction (build-skill corollary 1: a hand-typed CAS object is still hand-typed). Name any
   leaking clause verbatim.
4. **Consistency.** Does anything in the directive contradict #88's standing result (its energy-basis-completion
   disturbance measurement stands; its full-row/KINETIC adjudication is redone AFTER this fold), the S11b coupling
   law, or the §5c uniform-limit regression against `S11B_INPLANE_EOM`?

## Output

Your verdict per question (CORRECT / DEFECT-with-fix / AMBIGUOUS), each with the governing spec line or code line
quoted. If you would change the directive, give the exact replacement text. A prose claim without a quoted
spec/code citation — or, for a spot-check of the settled facts, a runnable script and its literal stdout — is
discarded. Do not edit any file; this is a review.
