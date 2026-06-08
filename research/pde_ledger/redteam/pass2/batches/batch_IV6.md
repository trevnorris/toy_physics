# Pass-2 Batch IV.6 (stages 151–163) — summary

**Part IV.6 — Correction, coevolution, traction, off-family.** 13/13 verified, all
`material_change: false` (NO downstream staling), 0 stop-cold/blocked, all Codex iter-1 exit 0
(no iter-2), 0 Codex deviations. **CHECKPOINT 163 cleared the higher bar.** Status-only by design
(no engines): **153** (consolidation/carry-forward card; values trace to upstream Stage 152→147 —
same class as 103/113/120/124). **NO EM-projected stages.** Value reconciliation: **0 misaligned
batch-wide; NO `paper_misalignment` → ZERO paper/notes edits.**

## Disposition
- **9 dual-engine script-side clean:** 151, 152, 154, 155, 156, 157, 159, 160, 163 (`.wl`
  independence orchestrator-confirmed; value-reconciliation 0-misaligned).
- **1 status-only:** 153 (no scripts — by design; carry-forward from Stage 152, upstream dual-engine).
- **3 script-side findings, all resolved:** 158 (transliteration re-author — ORCHESTRATOR OVERTURN),
  161 (transliteration re-author — audit-flagged), 162 (stale-numbering comment).

## 158 — F1 `mathematica_transliteration` → FULL re-author (USER-AUTHORIZED; ORCHESTRATOR OVERTURN)
**The standout.** The audit agent called 158 **CLEAN**, declining to file under the framing "no
independent route exists for a pure-Taylor-coefficient stage." The orchestrator's ground-truth
`.wl`-vs-`.py` read OVERTURNED that — for the two load-bearing checks the `.wl` ran the IDENTICAL
Taylor-expansion choreography on the IDENTICAL shifted closed forms as the `.py`: "linear delta R law"
`rLin = Normal[Series[rFun /. g->gStar+dg, {dg,0,1}]]` mirrored `R_lin = sp.series(R.subs(g,g_star+dg),
dg,0,2)`; "linear Delta_Q law" `chiLin = Normal[Series[chi,{eps,0,1}]]` mirrored `chi_lin =
sp.series(chi,eps,0,2)`; the gain/slope checks re-typed the same `expand`-and-drop-cross-terms. This is
the IV.5-139 lesson exactly (audit agents UNDER-call transliteration; the "no independent route"
defense is rejected when one IS feasible). The orchestrator surfaced the user-level re-author-vs-accept
call (surfaced, NOT reversed unilaterally; same class as 139/100/117). **User AUTHORIZED re-author.**
Codex re-authored the `.wl` to derive every linear law via analytic differentiation at the base point
(`rSlope = D[rFun,g] /. g->gStar`; `D[mQBase,...]`/`D[piBase,...]` partials; `D[chi,eps]/.eps->0`;
`dRFromDg` from `D[(g-r)^2/(1+r^2),g]/.g->gStar`) — genuinely independent of Series-on-shifted-form,
and can-fail (each compared to a separately hand-typed closed form). All 6 `expectZero` identities and
targets UNCHANGED; all 8 printed numeric coefficients byte-identical in value; **committed `.wl` output
BYTE-IDENTICAL to HEAD** (method-only). `.py` UNCHANGED. **158 ADDED to the Independent-Mirror Set.**
`material_change: false`.

## 161 — F1 `mathematica_transliteration` → route-only re-author (USER-AUTHORIZED; audit-flagged)
The audit agent correctly FILED this one: the `.wl` was a near-line-by-line port — same `eps->0`
derivative trick (`dBW = D[bWPert,eps]/.eps->0`), same `epsKExact`/`epsGExact` closed forms + the same
`12 lW^2 -> Pi^2 a^2 (1+rc)` branch substitution, same `upsilonPi`/`deltaQ` construction, the same
re-typed Stage-159 transport literals (`0.832…`/`1.162…`/`6.429…`), and even a **verbatim-ported inline
comment** (wl:72 ↔ py:81). Only `PolynomialRemainder` vs `.subs(poly,0)` differed in 2/9 checks.
**User AUTHORIZED route-only re-author.** Codex re-authored: `dB_W` via `Normal[Series[…,{eps,0,1}]]`
+ `Coefficient` (vs the `.py`'s `diff`); the D/N even/odd defects via **logarithmic derivatives**
`D[Log[kappaRatio],…]` and `D[Log[gamma0Sym/(1+rc)],…]` (no `PolynomialRemainder` survives — the `dln`
target emerges directly); the Stage-160 prefactor made an explicit named factor; the ported comment
deleted. All 9 `expectZero` identities + targets UNCHANGED; prefactor `(1+r_F1^2)/9 =
0.46236233468786880105…` and the `Delta_Q in (dThat,dS)` coeffs (5.352238871696225 / 10.70447774339245
/ −1.16275838754222) preserved. Committed `.wl` output **form-only** change (the `d eps_kappa`
intermediate now prints the reduced form `(-2 da)/a + (2 dLW)/lW − drc/(1+rc)`; the
`d eps_kappa identity = 0` check passes either way — same value). `.py` UNCHANGED. **161 ADDED to the
Independent-Mirror Set.** `material_change: false`.

## 162 — F1 stale-numbering comment (NOT user-level; Codex-applied, 0-seat)
`scripts/…stage162…_sympy_audit.py:39` `# Exact parent family formulas from Stages 99 and 102.` →
`# Exact parent family formulas from Stage 119 (with gamma0 from Stages 115-116).` Sibling of the
+17-era garble the Step-0 fix corrected in the notes (`Stages 99 and 170`→`Stage 119`); this `.py`
carried it as `99 and 102`. Stage 119 owns the parent-compensation family + balance law +
`L_W/a=(π/2)√((1+𝔯²)/3)`; `γ₀=(1+𝔯²)/9` traces to Stages 115–116. Non-load-bearing comment; the math
and all three `expect_zero` checks are unaffected; committed output byte-identical. `.wl` untouched.
`material_change: false`.

## Orchestrator backstop (transliteration watch — IV.6 = tail of the 105–175 band)
Ground-truth `.wl`-vs-`.py` read on ALL 12 dual-engine stages. **158 was the standout** (CLEAN→findings
overturn; the 139/100 lesson reconfirmed — audit agents UNDER-call transliteration, the orchestrator
read is the backstop, re-author-vs-accept is USER-LEVEL). 161 audit-flagged correctly. 154 borderline
but ACCEPTED (engines vary mechanism: `.py` exact `expand` vs `.wl` `Series` on the shifted-R; `.subs`-
drop vs `Series` on the linearization — meaningfully less port-like than 158; same disposition as 150).
155/156 numerical independence CONFIRMED via committed-output last-digit divergence (155: 30-digit
mpmath vs machine doubles; 156: independent Picard/bisection runs). 151/152/157/159/160/163 confirmed
genuinely independent: **151** sampled recursive-IBP integrator at 5 rational `Pi_*` vs symbolic
`Integrate` with symbolic `piStar` (the documented hung-Mathematica methodology — re-confirmed exec
deterministic <600s, math reformulated not cap-raised); **152** `.wl` re-derives `Pi_*`/`g_*`/`A_T`/`B_T`
via `FindRoot`/`D[` what the `.py` hardcodes as `mp.mpf` literals (MORE independent); **157** native 2×2
`Det[[1,−9σ],[5,−72σ]]=−27σ` non-degeneracy + `Solve[…,Reals]` (StatusOpen even-map defer-to-158 left
as designed, NOT re-opened); **159** Jacobian-`D[` vs `.py` `series`+`Poly` filter; **160** hand-written
total differential vs `.py` automated `.series()`; **163** (checkpoint) `.wl` adds the IFT slope
`−F_r/F_g` + a full `Series` perturbation route the `.py` lacks → strictly stronger; its load-bearing
`4√(1+r*²)` re-derived from `∂F/∂g`, not trusted.

## Infrastructure
- **6 exec exit 0** (reliability gate): 158 sympy byte-identical (`.py` unchanged) + 158 mma
  byte-identical (method-only re-author) + 161 sympy byte-identical (`.py` unchanged) + 161 mma form-only
  + 162 both byte-identical (comment-only). `$RT exec-*` writes `exec_logs/` only → orchestrator
  sed-refreshed every committed `.txt`.
- **Arbiter grep** on committed outputs CLEAN of stale self-epoch (NNN−17 = 134–146) banners; no `168π²`/
  `168%` class anywhere in IV.6; canonical Family-1 radius `√(4107−100π²)/(10π)` correct at stage 157.
- **Numbering:** IV.6 cards CLEAN of the +17 `\stagefield{Purpose}` self-label class (168–180 absent;
  151/162/163 `\stagefield{Purpose}` carry the canonical self-label). Step-0 garble fix (`Stages 99 and
  170`→`Stage 119`) re-confirmed holding in notes 162/163; the sibling `.py` comment garble (`99 and
  102`) fixed this batch (stage 162).
- **Seat policy held:** 158∥161 = 2 `.wl`-touching Codex sessions concurrent at the cap; 162 = 0-seat
  `.py`-only concurrent; orchestrator `exec-*` sequential after all Codex done (no overlap with any
  active `.wl` session).
- 6 trackers synced (PAPER_CLEANUP **P5-15** = ZERO paper/notes edits). Pass-1 `MANIFEST.yaml` untouched
  (isolation held).
