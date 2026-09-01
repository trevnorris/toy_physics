# Decision-list review — S11c-b #89 WL §3a-repair directive (BEFORE any builder runs)

You are one of two independent decision legs reviewing an **orchestrator-written build directive** before a
Codex builder is launched against it. This is the one artifact the builder trusts — everything downstream is
checked twice, the directive is checked here. Your job is to find defects that would cost a build round or,
worse, produce a wrong-but-passing build. Report them; do not fix them.

## The artifact under review
`research/pde_ledger_v3/directives/S11c_b_89_wl_3a_repair_directive.md`

## What it is trying to do
Repair the Wolfram engine's §3a energy-basis enumeration. The engine hand-codes a fixed family of new
"spurion" invariants (background-first-jet field-bilinears); the spec requires **every** symmetry-allowed
such invariant, and the hand-coded family is a strict subset. The directive asks the builder to replace the
hand-picked family with a systematic complete enumeration and let the engine's existing incremental-`MatrixRank`
selection compute the independent count — **without being told the completed count** (it is withheld; the
only public target is the current incomplete-family regression value `26`).

## What you are handed (read all of it; a directive review checks the directive against its sources)
- The directive above.
- The spec it must satisfy: `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` (esp. §1d ~155–168,
  §3a ~242–270). The directive quotes it — verify the quotes are faithful and the obligation is correctly
  transcribed, not narrowed or widened.
- The engine it repairs: `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
  (esp. the four hand-coded lists L390–436, the rank machinery L590–618, `applyProfile`/`profileDerivativeRules`
  L297–364, `profileRulesRetainingGeneratedJets` L339–361, the consumers L501–576, the emits L1236–1275).
- The prior finding it rests on: `directives/_measurements/S11c_b_87_wl_subspace_result.md` (WL's current
  family is a strict subspace of the correct family).
- The orchestrator's spec-confirm reasoning: `directives/_measurements/S11c_b_89_wl_spec_confirm.md`.
- Context: the sibling PY repair directive `directives/S11c_b_89_sympy_3a_repair_directive.md` and its
  clearance `_measurements/S11c_b_89_sympy_buildleg_clearance.md` (for the parallel structure ONLY; the two
  engines are built independently — the WL engine must not transliterate the PY one).
  ⛔ The completed count is withheld from the builder on purpose. If you compute it yourself while checking,
  do **not** put the number in a form the builder would read; report only whether the directive's *logic* is
  sound. Your report is read by the orchestrator, not pasted into the build.

## The failure modes to hunt (each has cost this program a build round)
1. **Spec fidelity / level-above miss.** Does the directive faithfully state §3a's obligation? Is the
   obligation itself correct, or does the directive (or the spec) encode a subtler error one level up? In
   particular: §3a says "independence is judged as field bilinears." **Does that mean raw linear
   independence (what WL's `MatrixRank` does) or independence modulo total in-plane divergences (a quotient,
   which WL does NOT take)?** If those two readings give different counts for the *complete* family, the
   directive's premise (WL reaches the correct answer by plain rank once the family is complete) is wrong,
   and the two engines will disagree for a reason the directive mis-attributes. Settle which reading §3a
   compels, and whether raw and quotient independence coincide for the complete family (they do only if the
   quotient's nullity is zero with the first jet live). If you can, compute it in an abstract field-bilinear
   space and **save the script + its literal stdout to a named absolute path**; a prose argument that they
   coincide is worth nothing here.
2. **Freezing a varying field (rule 17).** WL's rank reduction (`independentRepresentativeIndices` → the
   standard `applyProfile`) zeroes every 2nd/3rd background jet (`profileDerivativeRules` `higherRules`,
   L304–307). The directive's §5 requires the completed rank to be computed BOTH with that frozen path and
   with the higher-jet-retaining path (`applyBackgroundProfileWithGeneratedJets`) and the two ranks emitted.
   Is that control correctly specified and able to fail? Is there any completed-family invariant whose
   spurion is itself differentiated (a Hessian atom) that the frozen path would silently zero — making the
   first-jet `rankVariables` space inadequate? If so the basis is under-counted by a *frozen* mechanism and
   the directive must say so louder.
3. **Naming the object vs the recipe / leaking the answer.** Does the directive name the object (the complete
   symmetry-allowed family) without handing the builder the specific missing forms or the completed count?
   Grep it yourself for the withheld count and any phrasing that leaks the magnitude or a specific new form.
   Conversely, is the enumeration rule specified *precisely enough* to be complete — or is it vague enough
   that a builder could produce an under- or over-complete family and still exit 0? Under-specification has
   cost this ledger more than contamination.
4. **Designed-away disagreement (rule 1/6).** The directive forbids transliterating the PY enumeration. Is
   that prohibition backed by *absence* (the builder is given the physics, not the PY list) or only by a
   sentence? Could the builder still converge the WL family onto the PY one in a way that makes the
   cross-engine agreement vacuous?
5. **Acceptance that passes with the defect in place.** The only public acceptance is "restrict to the
   original 8 forms ⇒ `26`." Could a builder satisfy that (and every emitted control) while leaving the
   completion wrong or absent? Is the form ablation (restrict-to-8 must MOVE the count) decisive, or could a
   coefficient-only change fake it? Is any control tautological (residual zero for any input) or emitted
   conditionally on its value?
6. **Propagation / scope.** Completing the four lists is claimed to propagate to the count, the operator
   input, and the dimension/coefficient tables via `basisRepresentativeIndices`. Verify against the code that
   nothing keyed to the new-invariant labels is left stale (dimensions L1750, coefficient-dimension map
   L1831, the omissions emit L1275). Is any consumer of the family NOT reached by completing the four lists?
   Is the operator/coupling genuinely free of a separate frozen-jet defect (the directive claims WL's
   symbolic `gradient` makes it so — check that claim against the code)?

## Method
- This is a directive/document review: read the sources, form your own view of what §3a compels and what the
  WL engine does, then judge the directive against that. Quote both sides (directive line ↔ spec/code line)
  for every finding.
- Where you dispute a *physics* claim (raw vs quotient independence; first-jet adequacy; completeness of the
  enumeration rule), **demand computation of yourself**: write a small script in an abstract field-bilinear
  space, run it, and save the script + literal stdout to a named absolute path you report. A prose
  re-derivation is discarded.
- ⛔ Do not edit the directive or the engine. ⛔ Do not run the WL engine or any heavy Mathematica build; if
  you spawn a CAS kernel for your own abstract check, wrap it in `timeout 600` and run one at a time.
- Report a numbered list: `severity — directive:line — problem — why it costs a round — suggested correction`.
  Severity ∈ {BLOCKER (build must not launch), SHOULD-FIX, NIT}. Then a one-line verdict:
  `VERDICT: N findings (B blockers)` or `VERDICT: CLEAR`.
- A leg that finds nothing is weak evidence. If you clear it, name the two or three things you actively tried
  to break and how, so the orchestrator can weigh the clearance.
