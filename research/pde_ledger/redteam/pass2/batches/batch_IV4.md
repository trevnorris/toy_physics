# Pass-2 Batch IV.4 (stages 127–138) — summary

**Part IV.4 — Penetration, mouth boundary, fixedpoint.** 12/12 verified, all
`material_change: false` (NO downstream staling), 0 stop-cold/blocked, all Codex iter-1 exit 0
(no iter-2). **NO checkpoints in range. NO EM-projected stages.** Status-only by design: **128, 132,
136** (consolidation cards, no scripts of their own; all carry-forwards upstream-verified in both
engines). Value reconciliation: **1 misaligned (127, resolved); all other emitted deliverable values
reconciled.**

## Disposition
- **6 dual-engine script-side clean:** 127, 130, 131, 133, 137, 138 (`.wl` independence orchestrator-
  confirmed; value-reconciliation 0-misaligned on the script side).
- **3 status-only:** 128, 132, 136 (no scripts — by design; carry-forward values consistent;
  132 "180–182"→129–131 & 136 "184–186"→133–135 first-pass attribution corrections re-confirmed holding).
- **3 script-side findings, all resolved:** 129 (re-author + de-insufficient), 134 (re-author), 135 (de-taut).
- **1 notes value typo, resolved (correct-to-script):** 127 x*_exp `…160`→`…161`.

## 127 — F1 `paper_misalignment`/`value_mismatch` → notes digit fix (USER-RESOLVED correct-to-script)
Both engines independently compute `x*_exp = 0.6627654026231614025…`; the stage-127 notes box gave
`0.662765402623160` (15th sig-fig `0`, should be `1`). The companion slab depth in the same box already
matches the engines to 15 dp, fixing the intended precision. **User resolved in favor of the SCRIPT**
(scripts already correct). Codex (session NEW, 0-seat notes-only) edited
`notes/stages/moving_throat_pde_stage127_penetration_families.md:78` → `…161`; Claude reviewed (only the
digit changed, no `.py`/`.wl`/`paper/` touched). Published cards UNAFFECTED. `material_change: false`.

## 129 — F1 `mathematica_transliteration` + F2 `insufficient_verification` → FULL re-author (USER-AUTHORIZED)
Orchestrator ground-truth `.wl`-vs-`.py` read confirmed the audit agent's call: the `.wl` POSTULATED the same
hardcoded boundary-layer profile `sigma = piM·Exp[−piM·z/lM]/(lM(1−Exp[−piM]))` as the `.py` and ran the
identical normalization/zero-flux/residual checks with the identical `v1→piM·thetaSigma/lM` substitution — a
line-by-line port (pass-1 called 129 CLEAN). F2: the chemical-potential→Onsager-current link, the stage's whole
physical motivation, was never asserted (`mu` was dead code in the `.py`, absent from the `.wl`).
**User AUTHORIZED re-author** (re-author-vs-accept = USER-LEVEL; surfaced, not reversed unilaterally; same class
as IV.3-117 / IV.2-107/110/114). Codex (session 019..., iter-1):
- **F1** → the `.wl` now independently DERIVES the profile: `DSolve[thetaSigma·sigmaFn'[z]+v1·sigmaFn[z]==0,
  sigmaFn, z]` (the stationary zero-flux ODE) → normalize via `Solve[Integrate[sigmaGen,{z,0,lM}]==1]` →
  `sigmaDerived`, then `expectZero["derived profile matches boxed sigma_Pi", sigmaDerived − sigma]`. The boxed
  `sigma` survives only as the comparison target + for the downstream (already-present) checks.
- **F2** → added the Onsager-current identity `expectZero["Onsager current from mu identity",
  (−mobility·sigma·D[mu,z]) − jSigma]` (uses pre-substitution `jSigma`, holds for all v1) in the `.wl`, and the
  matching `J_from_mu = −M·sigma·diff(mu,z)`; assert `J_from_mu − J == 0` in the `.py` — resurrecting the
  formerly-dead `mu`.
Committed outputs ADDITIVE (new derivation/Onsager PASS lines); deliverables (boxed σ_Π, normalization=1,
J_σ=0, ODE residual=0, Π_m=V₁L/Θ) UNCHANGED. `material_change: false`. **129 ADDED to the Independent-Mirror Set.**

## 134 — F1 `mathematica_transliteration` → FULL re-author (USER-AUTHORIZED; ORCHESTRATOR-rewritten directive)
The `.wl` kernel `sKernel[p,k]` was character-for-character the SymPy closed form `S(Pi,kappa)`, and `S_shell`/
`S_q` were extracted from that POSTULATED kernel by the same operations in both engines (pass-1 batch-5 removed an
X−X gain-line assert but never made the `.wl` independent; `MATHEMATICA_MIRROR_POLICY` recorded "134 introduced
NO new mirror"). **⚠️ ORCHESTRATOR CATCH:** the audit agent's drafted directive only added a cosmetic
`Limit[kap→Pi/2]`-vs-substitution cross-check on the SAME postulated kernel — that would NOT make the `.wl`
independent. The orchestrator REWROTE the directive to require a genuine derivation (requirement+acceptance only;
Codex designed the route per [[feedback-claude-reviews-codex-codes]]). **User AUTHORIZED re-author.** Codex
(session 019..., iter-1): the `.wl` now solves the scalar mixed Dirichlet/Neumann BVP via
`DSolveValue[{−uFun''[x]+k²uFun[x]==G·sigma, uFun[0]==0, uFun'[1]==0}, uFun[x], x]` (mirroring stage-133's
independent route), reads `derivedKernel = u'(0)/G`, and sources BOTH `S_shell` (κ→0) and `S_q` (κ=Pi/2) from it;
the boxed paper closed form survives ONLY as the RHS of `expectZero["derived mouth kernel equals boxed paper
closed form", derivedKernel − paperKernel]`. Added BVP ODE-residual + both-BC checks. `.py` UNCHANGED. Committed
output additive; deliverables (S_shell=1, S_q numerics, fixed-point law) UNCHANGED. `material_change: false`.
**134 ADDED to the Independent-Mirror Set.**

## 135 — F1 `tautological_check` → de-taut (routine)
The `.wl:78` `expectApprox["closure residual", residual, 0, 10^-14]` was X−X: `residual = piStar −
sigmaStar·(4−sStar)` with `sigmaStar` solved from exactly `piStar == sigmaVar·(4−sStar)`, so the residual ≡ 0 by
construction (the SymPy side had already demoted it to a bare print). Codex DELETED the assertion (no
replacement) and KEPT the residual PRINT — restoring SymPy/Mathematica parity. The kernel value remains
exercised by the numeric anchors (`.wl:74-77`) and range checks. 135's `.wl` confirmed STILL independent
(carry-forward kernel + independent outlet-reduction identity `reducedLaw − expectedLaw`); NOT a mirror change.
`material_change: false`.

## Independence (orchestrator ground-truth `.wl`-vs-`.py` read on ALL 9 dual-engine stages)
**Backstop CONFIRMED every audit-agent call — NO under-calls this batch (contrast IV.1-100).**
- 127 **independent** — `.wl` derives both bias factors from the source integrals
  (`Integrate[(1/x)Cos[Pi z/2],{z,0,x}]`, `Integrate[(Exp[−z/x]/(x(1−Exp[−1/x])))Cos[Pi z/2],{z,0,1}]`) and
  checks them against the closed forms; `.py` postulates them.
- 130 **independent** — `.wl` derives gPi via `Integrate[sigma·f]` (not postulated) and proves global
  monotonicity via native `Reduce[Sin[…]>0]` (different procedure than SymPy's closed-form match).
- 131 **independent** — cleared-denominator polynomial residual `40πp(2pe^p+Π)−20π·gMinus(4p²+Π²)(e^p−1)` with a
  bracketed `FindRoot`, distinct from SymPy's `nsolve` on the rational. Sign verified `(e^p−1)` (NOT the `(1−e^p)`
  trap the heads-up warned of).
- 133 **independent** — `DSolveValue` BVP solve from scratch; `.py` uses a hand-ansatz.
- 137 **independent** — matrix-`Inverse[M_core]` Schur reconstruction of ρ_c/σ_c (mirrors verified stage-114),
  `Series` vs SymPy's `Limit` for the static lane, working non-vacuity guard.
- 138 **independent** — independent algebra; `.wl` adds an `R_q exact formula` check the `.py` lacks.
- 129 & 134 = the two genuine transliterations (the `.wl` *postulated* the same closed form + verified, vs.
  130/137 which *re-derive* via Integrate/matrix-inverse) → re-authored (above). **0 sanctioned mirrors remain in
  IV.4.**

## INFRA
- **6 orchestrator exec exit 0** (reliability/determinism gate: 129 both engines, 134 both, 135 both).
- `$RT exec-*` writes `exec_logs/` only → orchestrator refreshed every committed `.txt` from the authoritative
  run: **129 sympy** stripped a stale pre-pass-2 `# SymPy Audit Output / # Status: PASS` wrapper header (adopts
  the bare pass-2 format); **129 mma** gained the DSolve-derivation + Onsager PASS lines; **135 mma** dropped the
  stale `PASS: closure residual` line (Codex removed the assert from the `.wl` but had not refreshed the `.txt`);
  134 sympy/mma + 135 sympy byte-identical to Codex's run.
- **Arbiter grep on the refreshed committed outputs + changed scripts: CLEAN of stale self-epoch (−17 band
  110–121) self-banners.** Only hits = legit upstream CROSS-refs (137's `notes stage097`/`stage114` Schur
  provenance; 121 r_F1 owner cross-ref in out-of-batch 139) — content-correct, not self-labels.
- **Seat policy held:** 129 + 134 = 2 `.wl`-touching Codex sessions (concurrent, at the 2-seat cap, flock-safe),
  then 135 solo; 127 notes-edit = a separate 0-seat session (no `math -script`, ran concurrent); orchestrator
  exec sequential, after all Codex done (no overlap).
- Pass-1 `MANIFEST.yaml` untouched (isolation held). 6 trackers synced (PAPER_CLEANUP **P5-13**).

## Deferred (NOT this loop)
- **NONE new for IV.4.** The IV.4 cards are CLEAN of the +17 `\stagefield{Purpose}` self-label class that IV.3
  surfaced (grep for "Stage 144–155" found nothing). The 132/136 downstream-attribution corrections (Cluster B)
  were re-confirmed holding (zero residual "180–182"/"184–186").

## Independence outcome
2 newly-independent `.wl` (129 + 134 re-authored from transliterations to DSolve/DSolveValue BVP routes);
**0 sanctioned mirrors**; 129 & 134 added to the Independent-Mirror Set. All 9 dual-engine IV.4 stages now
genuinely independent. NO checkpoint constant changed (no checkpoints in range). Pass-1 `MANIFEST.yaml` untouched.
