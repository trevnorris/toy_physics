# Pass-2 Batch V.1 (stages 164–175) — summary

**Part V.1 — Microscopic logs, drifts, bundle inversion.** 12/12 verified, all
`material_change: false` (NO downstream staling), 0 stop-cold/blocked, all Codex iter-1 exit 0
(no iter-2), 0 Codex deviations. **This batch CLOSES the 105–175 first-pass orchestrator-direct
TRANSLITERATION WATCH — the watch ended at 175.** **NO checkpoints in range.** **NO status-only
stages** (every stage is dual-engine `.py`+`.wl` — confirmed; a missing engine would have been a
finding). **NO EM-projected stages.** Value reconciliation: **0 misaligned batch-wide**
(164=19, 165=7, 166=11, 167=20, 168=11, 169=10, 170=11, 171=18, 172=12, 173=12, 174=9, 175=11
deliverable values checked) → **ZERO `paper_misalignment` → ZERO paper/notes edits.**

## The user "re-author all 6" decision
The orchestrator ground-truth `.wl`-vs-`.py` read split the 12 into **6 confirmed genuinely
independent (no re-author needed)** and **6 confirmed `.wl` ports.** For close-of-watch strict
dual-engine hygiene the **USER DECIDED "re-author all 6"** of the ports. Codex designed + wrote
the new routes (the `.py` reference engines were left untouched except 165's F2 numeric block);
the verifier confirmed all 6 are genuinely independent routes (not relabels). All
`material_change: false`.

## Disposition
- **6 confirmed genuinely independent (clean, no re-author):** 164, 166, 170, 172, 174, 175.
- **6 confirmed ports → USER-DECIDED re-author all 6:** 165, 167, 168, 169, 171, 173.

## The 6 re-authors (USER-DECIDED, Codex-written, verify-confirmed independent)
- **165** — F1: `eqR`/`eqG` now ASSEMBLED from parent-coupling log-drift scalings
  (`dKs`/`dKq`/`dLam`/`dGs`/`dGq` from primitives → `chanR`/`chanG` emerge via `FullSimplify`, +
  matches-channel guards). F2 `insufficient_verification`: added 4 numeric prefactor checks in
  BOTH engines (`Tm_pref≈1.2715890393387603`, `v_pref≈1.1428896163056477`,
  `ratio_pref≈0.8987885086678338`, `prod_pref≈1.4532859092683434`), `rstar` from canonical
  `√(4107−100π²)/(10π)`. `.py` + `.wl` changed.
- **167** — F1: additive independent numeric-closure route (`numericClosureResiduals` recomputes
  the full primitive drift chain for 5 tuples incl. mixed `{2,-3,5,-7}`, asserts all
  invariants/channels vanish). Token algebraic independence (pure-substitution stage,
  user-accepted) but additive + route-independent. `.wl` only.
- **168** — F1: `epsPerp` now DERIVED via `Coefficient[Expand[-deltaPerpSlip], epsT/epsv/epsL]` +
  boxed-weight `expectZero` guard; decimal `rNum`→canonical `√(4107−100π²)/(10π)`. `.wl` only.
- **169** — F1: decimal radius `rNum`→canonical `√(4107−100π²)/(10π)`; the 3 paper-comparison
  TARGET literals (`0.758035078944663`, `1.00314310113848`, `1.88373219118005`) correctly
  preserved as paper-side values; sphere-avg + matrix-invariant sections already independent.
  `.wl` only.
- **171** — F1: additive `bCombSeries` (BdG Series-route in a fresh `eps2`) + `kScalar`/`gScalar`
  weak-axisymmetric scalar routes (the Z/N bundles ALREADY had independent Series routes from a
  prior pass). `.wl` only.
- **173** — F1: re-authored to direct analytic `D[ratio,eps]/.eps->0` (distinct mechanism from the
  `.py`'s `series().diff()`); SymPy-mirroring `d0A`/`d2A`/`d4A`/`n0A` names removed; even-preserving
  collapse via `Numerator`/`Together`/`Coefficient` (not `Solve`/`First`); the byte-identical
  carry-forward `Print` block replaced with a native one-liner; all 6 `expectZero` targets
  preserved. `.wl` only.

## The 6 confirmed genuinely independent (no re-author)
- **164** — independent series route via perturbation + `Series`/`Coefficient`.
- **166** — matrix-inverse route (batch-7 remediation).
- **170** — `D[,eps]` vs series (BATCH-1 `bda2107` rework; reliability re-run deterministic 14s
  < 600s).
- **172** — implicit-differentiation route.
- **174** — perturbed-expression `D[,eps]` route.
- **175** — `dlogSeries` route for the load-bearing `Sigma_N` (batch-8, user-accepted — the
  residual-port was NOT unilaterally reversed).

## Numbering / arbiter grep
- **Numbering CLEAN:** all 12 cards clean of the +17 `\stagefield{Purpose}` self-label class
  (181–192 absent; correct self-numbers).
- **Arbiter grep on committed outputs CLEAN:** no stale self-epoch band (147–158 = NNN−17), no
  `168π²`/`100π²` class. Canonical Family-1 radius `√(4107−100π²)/(10π)` used correctly
  (165/168/169 de-transcribe the decimal `1.77799353547498` to this closed form).

## Independent-Mirror Set
ADD **165, 168, 171, 173** (genuine re-authored independent routes). **167** = additive
numeric-closure (token but route-independent — note the limited-achievable-independence
rationale). **169** = de-transcription only (sphere-avg/matrix routes already independent). **175**
stays (batch-8 `dlogSeries`). 164/166/170/172/174 confirmed independent (no re-author).

## Infrastructure
- **Orchestrator independent exec re-run all 12** (reliability gate) exit 0, all FAIL=0.
- **Committed outputs:** 7 mma re-authored (165, 167, 168, 169, 171, 173) + 165 sympy refreshed;
  **170 BOTH engines NORMALIZED** (refresh stripped a stray `# exit_code: 0` trailer — the
  IV.2-108 class; deliverables byte-identical, NOT a math change); 164/166/172/174/175
  byte-identical (deterministic).
- **Seat policy held:** 6 `.wl`-touching Codex sessions ran in 3 waves of 2 (165∥167, 168∥169,
  171∥173) under the flock; orchestrator `exec-*` sequential AFTER all Codex done (no seat overlap
  with the orchestrator exec-mathematica). Pass-1 `MANIFEST.yaml` untouched (isolation held).
- 6 trackers synced (PAPER_CLEANUP **P5-16** = ZERO paper/notes edits).

## Close of the transliteration watch
With V.1 closed, the 105–175 first-pass orchestrator-direct transliteration watch ENDS at 175. The
recurring lesson held to the last batch: audit agents UNDER-call transliteration; the orchestrator
ground-truth `.wl`-vs-`.py` read is the backstop; re-author-vs-accept is USER-LEVEL (here the user
escalated to "re-author all 6" of the confirmed ports). All 12 V.1 dual-engine stages are now
genuinely independent.
