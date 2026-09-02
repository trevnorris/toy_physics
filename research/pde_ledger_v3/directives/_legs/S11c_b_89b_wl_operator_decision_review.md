# Decision-list review — S11c-b #89b WL operator-unfreeze directive (BEFORE any builder runs)

You are one of two independent decision legs reviewing an **orchestrator-written build directive** before a
Codex builder is launched against it. This is the one artifact the builder trusts — everything downstream is
checked twice, the directive is checked here. Your job is to find defects that would cost a build round or,
worse, produce a wrong-but-passing build. Report them; do not fix them.

## The artifact under review
`research/pde_ledger_v3/directives/S11c_b_89b_wl_operator_directive.md`

## What it is trying to do
Repair the Wolfram engine's §3b **slab-operator** construction. The engine builds the operator by reducing the
energy with a background-profile rule that **zeroes the 2nd/3rd spatial jets** of the background coefficients
`W_bg`/`μ_R,bg` (`applyProfile` → `profileDerivativeRules` `higherRules`) **before** taking the Euler–Lagrange
variation. The spec (§3b/§3a) requires the coefficients and their jet tower to stay LIVE through the
differentiation, retaining every higher spatial jet at its originating background-amplitude grade
(`η,σ_W ≤ 1`) — the retained-grade background curvature `∂²W_bg`/`∂²μ_R,bg`. The directive asks the builder to
take the EL variation on the LIVE energy and reduce afterward with the jet-retaining path
(`applyBackgroundProfileWithGeneratedJets`), so the emitted operator (and the §3c kernel, μ_θ operator, and
term-origins derived from it) carries the curvature — **without being told the corrected operator's rank or
term counts** (those are withheld; the only public target is that re-freezing reproduces the current frozen
operator). This brings the WL operator to parity with the PY side (PY #89 already un-froze its operator).

## What you are handed (read all of it; a directive review checks the directive against its sources)
- The directive above.
- The spec it must satisfy: `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` — esp. **§3b ~272–288**
  (the divergence-form operator; "do not freeze a coefficient at its constant binding before differentiation"),
  **§3a ~242–256** (first-background-jet order bounds *amplitude* factors, not spatial-derivative order; "a
  second spatial derivative of `W_bg` is still first order in background amplitude … not dropped"), **§1d
  ~163–168**, **§3c ~290–323** (the kernel is a weak block of the §3b operator), **§3d ~325–353** (the
  separate admissibility ε⁰ route). Verify the directive's quotes are faithful and the obligation is
  transcribed, not narrowed or widened.
- The engine it repairs: `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` —
  locate by NAME (line numbers below are post-#89a and approximate): `evaluatedModel` (~L1130), `applyProfile`
  / `profileDerivativeRules` / `higherRules` (~L297–364), `profileRulesRetainingGeneratedJets` /
  `applyBackgroundProfileWithGeneratedJets` (~L339–361), `constrainedRows` / `variationalSource` (~L1022 /
  ~L275), `rawModel` (~L1049), `faceSources` / `evolutionSource` / `faceGeneralizedRows` (~L937 / ~L907 /
  ~L1039), `extractCoupling` (~L1376), the operator emits (~L1798–1838), `operatorFreezeRankDiagnostic`
  (~L1219–1244), and the §5.E dimension residual (~L2269–2290).
- The orchestrator's spec-confirm reasoning: `directives/_measurements/S11c_b_89b_wl_operator_spec_confirm.md`.
- The #89a clearance (the basis is done, the operator was diagnosed-not-repaired):
  `directives/_measurements/S11c_b_89_wl_basis_buildleg_clearance.md`.
- The measured downstream consequence and the sibling parity:
  `directives/_measurements/S11c_b_88_blast_radius_result.md` and the PY directive
  `directives/S11c_b_89_sympy_3a_repair_directive.md`. ⚠ **`S11c_b_88_blast_radius_result.md` contains the
  corrected operator's PY-measured rank-gain / term counts — those are the WITHHELD acceptance values.** If you
  reference them while checking, do **not** recommend putting any of them in the directive, and do **not** put
  them in a form a builder would read. Your report is read by the orchestrator, not pasted into the build. The
  two engines are built independently — the WL engine must NOT transliterate the PY one (rule 1).

## The failure modes to hunt (each has cost this program a build round)
1. **The jet DEPTH the operator actually needs (rule 17, the central physics question).** The directive
   requires the operator to carry the 2nd (and higher) spatial jets of the background at retained grade. Is
   that the *correct* depth? Independently determine which spatial-jet orders the divergence-form EL variation
   of the §3a energy actually generates for the operator, and at what amplitude grade each survives. Then
   check the two profile rules against that: `higherRules` (~L304–307) zeroes `2 ≤ i+j+k ≤ 3` (the freeze);
   `profileRulesRetainingGeneratedJets` (~L339–357) retains `i+j+k ≥ 1` up to order 3 per direction, scaled by
   a single `σ_W`. **Does the live path retain every jet the operator generates, or could the operator
   generate a jet order beyond what that rule covers (a residual freeze the directive does not catch)? Does the
   single-`σ_W` retained-grade scaling correctly bound the amplitude, or does it keep a term that should be
   truncated / drop one that should survive?** If you can, compute the generated jet content in a minimal
   symbolic operator and **save the script + its literal stdout to a named absolute path**; a prose argument is
   discarded.
2. **Does the fix actually un-freeze — or only relocate a freeze?** The directive's fix is "EL on the live
   energy, then reduce with the jet-retaining path." Verify against the code that this ordering is realisable
   for the WHOLE operator — not just the internal bulk rows but the μ_θ operand, the face/flux rows
   (`faceSources`/`faceGeneralizedRows`), and the evolution rows (`evolutionSource`), each of which the current
   `evaluatedModel` `applyProfile`s early (~L1141–1152). Is there any surface where the live energy is not
   available before the variation (e.g. a face term already reduced inside `faceSources`), so the builder would
   be forced to reduce early and re-freeze it? If so the directive under-specifies and the fix would be
   partial.
3. **Naming the object vs the recipe / leaking the answer.** Does the directive name the object (the live-jet
   operator) without handing the builder the corrected rank/term counts? Grep it yourself for any corrected
   operator magnitude. Conversely — is the construction specified *precisely enough* to be complete, or vague
   enough that a builder produces a wrong-but-passing operator and still exits 0? And the reverse over-reach:
   does §3 dictate a *recipe* (e.g. forcing use of `rawModel`) where it should name the object and let the
   engine build it (rule 3)? Under-specification and over-specification both cost rounds.
4. **Blast-radius completeness (scope).** Does the fix reach EVERY emitted operator surface — `SLAB_OPERATOR`,
   `SLAB_OPERATOR_TERM_ORIGINS`, `MU_THETA_OPERATOR`, `COUPLING_KERNEL`, `COUPLING_KERNEL_TERM_ORIGINS`, both
   anchoring branches, both density representatives (~L1798–1838, ~L1951–1987)? Is any consumer left reading a
   frozen operator? Is the term-origin provenance re-derived from the corrected rows or carried stale? Is the
   §3d admissibility ε⁰ route correctly left ALONE (the directive says verify-consistent, not double-fix — is
   that the right call, given it is a separate already-live route)? Are the `formOnly`/corrupted/parametric
   and §5c uniform-limit branches kept working?
5. **Acceptance/controls that pass with the defect in place.** The public regression is "re-freeze reproduces
   the current frozen operator." Could a builder satisfy the mandatory FORM ablation (control 1: re-freeze must
   MOVE the operator) with a change that is not actually the un-freeze — e.g. a coefficient rescale, or a
   cosmetic term that moves the output without carrying the curvature? Is control 2 ("retained-grade curvature
   present": the corrected operator depends on the 2nd-jet atoms the frozen one lacks) **decisive**, or could
   it be faked / is it tautological? Is control 3 (one-sided corruption of one branch's reduction) a genuine
   independence test? Is any control emitted conditionally on its payload's value, or a residual that is zero
   for any input (rule 2 corollary 3)?
6. **The §5.E fold.** The directive replaces a tautological dimension residual (both operands = the same
   factor-metadata sum) with "dimensional analysis of the emitted invariant EXPRESSION vs the factor-metadata
   sum," or deletion. Is the proposed replacement genuinely a *second, independent* route (does reducing the
   actual invariant expression's atom-dimensions re-use the same lookup, making it tautological again)? Would
   the "delete + note" option hide a real dimension error? Confirm the current residual is in fact structural
   zero (read ~L2269–2290).
7. **Withholding.** Is the withheld class correct — basis count 40 public (already emitted/cleared), corrected
   operator rank/term counts withheld? The freeze diagnostic already emits an abstract frozen-vs-live rank;
   does keeping it (control) create a tuning target for the *production* operator, or is it safely a different
   (abstract, minimal) object? Does the directive anywhere state a rank/termcount as a pass condition
   (rule 5)?

## Method
- This is a directive/document review: read the sources, form your own view of what §3b/§3a compel and what
  the WL operator construction does, then judge the directive against that. Quote both sides (directive line ↔
  spec/code line) for every finding.
- Where you dispute a *physics* claim (the jet depth of #1; whether the live path retains what the operator
  generates; whether a control is decisive), **demand computation of yourself**: write a small script (an
  abstract divergence-form operator, or an extraction of the relevant WL functions), run it, and save the
  script + literal stdout to a named absolute path you report. A prose re-derivation is discarded.
- ⛔ Do not edit the directive or the engine. ⛔ Do not run the full WL engine or any heavy Mathematica build;
  if you spawn a CAS kernel for your own abstract check, wrap it in `timeout 600` and run one kernel at a time.
- Report a numbered list: `severity — directive:line — problem — why it costs a round — suggested correction`.
  Severity ∈ {BLOCKER (build must not launch), SHOULD-FIX, NIT}. Then a one-line verdict:
  `VERDICT: N findings (B blockers)` or `VERDICT: CLEAR`.
- A leg that finds nothing is weak evidence. If you clear it, name the two or three things you actively tried
  to break and how, so the orchestrator can weigh the clearance.
