# Question-vet — what is the correct SUFFICIENT test to close per-engine N6 (non-circular, tractable)?

You are helping me settle a **contested framing** before I commit to an instrument. ⛔ Document/reasoning only; ⛔ do
NOT modify the tree; ⛔ run no CAS (you MAY read + reproduce the committed tally parser to check numbers). Working dir
`/var/projects/toy_physics`; paths under `research/pde_ledger_v3/`. This is the E1 question-vet: I own the question and
adjudication; any instrument will be astra-authored + G1-reviewed. ⛔ Do NOT hand me the physics answer (whether N6
holds); tell me whether a sufficient test is NEEDED and, if so, what it is — **non-circular and tractable**.

## What is ESTABLISHED (per-engine SymPy, reconcile instrument `scripts/S11c_c2_N6_reconcile_sympy.py`, both build legs CLEAR)
Tally in `_measurements/S11c_c2_N6_reconcile_adjudication.md` (4 cases; reproduce via
`_measurements/S11c_c2_N6_reconcile_disposition_tally.py` on `/tmp/S11c_c2_N6_reconcile_sympy.*.out`, or regenerate a
case with `python3 scripts/S11c_c2_N6_reconcile_sympy.py --anchoring LAB_HELD --density RHOBR_CONSTANT`):
- **Geometric carrier reconciles:** `CARRIER_BRIDGE_RESIDUAL = C_E − C_M` no-nonzero all 4 cases; a build-leg-verified
  LIVE control (material-normal FORM knife moves it). ⇒ the mechanical/geometric part is representation-invariant.
- **Residual localizes to the constitutive source channel:** `R_N6 = B(C_M, ΔS)`, `SPLIT_CHECK=0`,
  `CARRIER_CHANNEL=CROSS=0`. `ΔS = es − ms` (the diagnostic `source_terms` circuits: Eulerian μ_E + Eulerian V vs
  material μ_M + material V; `μ_M` via `b.material_pullback`, route-2 spec §2/§5).

## The CONTESTED open question (the two adjudication-review legs SPLIT on it)
Whether the nonzero source discrepancy (`SOURCE_BRIDGE_RESIDUAL = es − ms`, certified nonzero in 3 of 4 cases; `R_N6 =
B(C_M,ΔS)`) is the **SANCTIONED covariant transformation** under the material↔Eulerian field redefinition, or a **real
non-covariance**.
- **Codex-sol (over-clear finding):** `R_N6 = B(C_M,ΔS)` is an algebraic tautology given `ΔC=0` — it holds for ANY
  source discrepancy, so localization ≠ proof of sanctioned-ness. Needs an INDEPENDENT test: (a) quotient reduction of
  `R_N6` by the ideal `θ_M−θ_E−a_ρ=0`, `e_M−e_E−h_α=0` (+ prolongations) with a vanishing remainder, OR (b) an
  independently-derived PREDICTED source transformation matched to the actual `ΔS`.
- **Grok (sound finding):** `ΔS` is built from the predeclared maps and vanishes iff they are trivial (the `a_ρ+h_α`
  truth table); a wrong coefficient INSIDE the parent `material_pullback` would share that pattern and is **a
  parent-map question, not an N6 remainder** (the parent map is already reviewed at c1/S11c-a).

## The sub-questions I need settled (this is the crux — a circular or intractable test is worthless)
1. **Is a sufficient test even needed to close per-engine N6, or is Grok right that it is out of N6's scope?** State
   precisely what the per-engine N6 CLAIM is: is it "the two routes agree GIVEN the (parent-verified) builders" — in
   which case carrier-reconciliation + source-localization-to-the-parent-pullback-channel IS the complete per-engine
   statement — or is it "representation invariance holds" which requires validating the source transformation is
   covariant, INDEPENDENTLY of `material_pullback`?
2. **Circularity.** The material route already IMPLEMENTS `θ→θ+a_ρ` via `material_pullback`. So `ΔS = es − ms` is BY
   CONSTRUCTION `b(θ+a_ρ) − b(θ)`. Does option (b) — "independently derive the predicted transformation and match
   `ΔS`" — reduce to re-implementing `material_pullback` (⇒ circular/vacuous)? If NOT, what makes the independent
   derivation genuinely independent (a different route to the covariant transformation of `b` that does not import the
   builder under test)?
3. **Well-posedness of the quotient reduction.** In the current construction the material and Eulerian fields are the
   SAME sampled variables (shared PIT samples), not independent symbols, so `θ_M−θ_E−a_ρ` is already identically
   imposed. Is option (a) well-posed here, or does it require rebuilding both routes with θ_M, θ_E as INDEPENDENT
   variables and imposing the relation as a constraint — a redesign? Which?
4. **Tractability.** If a sufficient test exists and is non-circular, is it tractable (compact carrier + PIT, like the
   reconcile), or does it hit the F/G-style full-symbolic wall? If it risks the wall, say so and give the tractable
   form or say none exists.
5. **What, concretely, would DISTINGUISH sanctioned from a real non-covariance** — i.e. an input under which a genuine
   representation-invariance failure gives a DIFFERENT result than the sanctioned transformation (an able-to-fail
   discriminator). If no such discriminator exists at the retained order, that itself is the answer (N6 is
   necessary-only here).

## Output
For each of 1-5: your position + reasoning, grounded in the cited objects. End with a **one-paragraph VERDICT**: is a
sufficient test needed; if yes, its exact non-circular tractable form (the object it computes + the able-to-fail
discriminator); if no, why carrier-reconciliation + source-localization is the complete per-engine N6 statement.
Brief, evidence-first.
