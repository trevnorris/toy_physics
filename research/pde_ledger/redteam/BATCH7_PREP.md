# Batch 7 prep digest (recon complete — NOT yet executed)

> Written 2026-05-29 during "prepare for batch 7". Reconnaissance only; the fix loop
> has NOT started. Resume here after the user's explicit GO. Per-stage recon was done
> by 5 clean read-only agents against each `codex_reviews/stage_NNN.md` (authoritative),
> the current `directives/stage_NNN.md` (the ORIGINAL tainted prescriptions — all 2026-05-28,
> `applied:false`, and NONE faithfully encode their review), and the live scripts.

## Scope decision (for the user gate)
- Remaining FINDINGS stages: **148, 150, 157, 166, 175** (5 total — the last of the 29).
- Handoff says "pick next ~4 by ascending #": **batch 7 = {148, 150, 157, 166}**, then **175 as the final batch 8**.
- Alternative the user may prefer: fold 175 in → one final push of all 5. 175 needs a Mathematica-feasibility consult (see below), so it pairs naturally as its own small unit. **DEFAULT = {148,150,157,166}; await user.**

## Per-stage synthesis

### 148 — representative_positive_families  (review: 2 findings | HEAVIEST; has a LIVE bug)
- **R1 (paper_misalignment / insufficient assertion):** the `ξ_*` closed-form bridge `(1−λ_{Π,0}) − ξ_* == 0`. SymPy only `print`s, weakened to `1e-15` (residual 7.8e-17); review wants raw `<1e-25`. Also a **same-source** concern: both sides flow from `gMinus` (not strict X−X, but not independent).
- **R2 (insufficient_verification) — REAL LIVE DIVERGENCE BUG:** the reworked Mathematica `dSigmaOfDeltas`/`dTOfDeltas` route (wl:43-48,53,63,73) does NOT reproduce SymPy's first-order corrections: SymPy `dTU=0.5087563022150839`, `dTD=-0.1169438021518107`, `dTLam` slope `-0.625700104366894`; Mathematica gives `dTU=0.4976…`, `dTD=-0.1144…`, `dTLam=-0.6121…`. Root cause (agent's hypothesis): Mathematica treats `dPi`,`dS` as independent increments but SymPy's `AT` bakes in the `Sp_star = dS/dΠ` cross-term. **Nothing currently asserts cross-engine agreement, so both pass despite diverging.** Needs the correct first-order chain-rule decision (consult).
- **168π² → 100π² carry-forward (ACTIVE):** directive lines **27, 43, 59, 68, 138** carry stale `168`; live SymPy:87 + Mathematica:84 correctly use `100`. Confirmed from first principles: `rF1 = √(12·(37/20)²/π² − 1)`, `12·(37/20)²=4107/100` ⇒ `rF1²=(4107−100π²)/(100π²)`. `100` is forced; `168` is the Codex/stale-notes false-positive. (Do NOT touch wl:36 `MaxIterations->100`.)
- Directive findings_count=3 → really 2 (F2 "subsumed by F1" stub + already-moot hardcoded-float + already-fixed banner). All directive "before" quotes are stale (pre-restructure code). Notes anchor: `notes/stages/...stage148....md` (correct, uses 100π²); `notes/.../review/stage_148_review.md` is MISFILED (reviews Stage 029) — ignore.
- **Consult Qs:** (1) genuine-independence form for the ξ_* bridge (external Stage-126/228 numeric anchor vs same-source re-derive); (2) **the correct dΣ/dT first-order chain rule (the live bug) — free `dS` vs `Sp_star·dΠ`**; (3) can SymPy honestly hit `<1e-25` (raise nsolve precision / go symbolic) or is it nsolve-float-capped near 1e-17.

### 150 — full_profile_residual  (review: 1 finding | LIGHTEST; near-mechanical)
- **R1 (insufficient_verification, display-only):** transcripts (`output/...txt:5`) print the fully-substituted rational instead of the compact `S_q(Pi)=Aq*k−Cq*Pi` (SymPy) / `aq*k−cq*p` (Mathematica) the directive wanted, because `Aq,Cq` are concrete exprs that auto-expand under pprint/fmt.
- The directive's F1 (source de-tautology `Sq=Aq*k−Cq*Pi`) is **already applied + committed (3e2b5c0)**; directive anchors drifted +2 (assert now py:47/wl:46). So the rewrite must PIVOT from the stale source-edit to the transcript-display fix.
- **No `sp.integrate`/`Integrate` anywhere** — "residual" is pointwise derivative/curvature at x=0. No hang/precision risk.
- **Consult Q (one-liner):** prefer display approach (b) — build `Sq` from free symbols `Aq_s,Cq_s`, print, then `.subs` concrete defs for the load-bearing asserts — over approach (a) hardcoded display string (needs a fabricated-display guard). Confirm (b) doesn't disturb `T_q'(0)−S_q==0`.
- Out-of-scope flag: `notes/.../review/stage_150_review.md` mislabeled (body is Stage 031).

### 157 — core_mouth_coevolution_status  (review: 3 findings | HARDEST conceptually; NOT status-only)
- MANIFEST `is_status_only_candidate:false` — carries real load-bearing checks.
- **R1 (insufficient_verification):** SymPy `sol_deltaC` re-solve (py:112-114) duplicates the py:107-110 solve of the SAME homogeneous `{dE2=0,dE4=0}` pair → `deltaC≡0` by self-solved root (tautological). The old `-16 σ_*·dR` projector route was itself tautological (multiplies a known-zero `dR`).
- **R2 (transliteration):** Mathematica `solDeltaC` (wl:102-105) is the SAME hard-coded 2×2 numerator system `{dC−9σ dK, 5dC−72σ dK}` — a mirror; wrong literals would pass in both engines. The `9/72/5/27/243` coefficients have no in-stage provenance (borderline fabricated-literal).
- **R3 (symbol_assumption_error):** wl:93 `$Assumptions` is `Reals`-only; needs `0<sigmaStar<1` (claim invalid at σ=0 for dKappa; `(1−σ)` denominators singular at σ=1). **The directive does NOT address R3** (directive's F3 is an informational mtime item) — F3↔R3 is a full mismatch.
- Banner carry-forward: SymPy:63 says "Stages 138-139" but .wl:50 + docstring + notes say **155-156** → fix SymPy banner to 155-156, refresh both .txt.
- **Consult Qs (real independence judgment):** (1) Is an honest independent `deltaC(family-motion)` route even constructible in-stage, or is the tangent→defect map deferred to Stage 158 (→ then the right move may be RELABEL + defer, like batch-5 134's outlet-to-135)? (2) Can the canonical-even `9/72/5` coefficients be reconstructed from an upstream Galerkin source, or are they not available here? (3) Does adding `0<σ<1` make Solve/FullSimplify return a `ConditionalExpression` that breaks `expectZero` (strip per the idiom memo)?
- Notes anchor real; `notes/.../review/stage_157_review.md` is orphaned (Stage 038) — ignore.

### 166 — bundle_inversion_four_drifts  (review: 1 finding | CLEAN, well-understood)
- **R1 (tautological_check, blocking):** wl:76 `expectZero["matrix round-trip", Total[(Mmat.solVec − v)^2]]` with `solVec=Inverse[Mmat].v` ⇒ `Mmat·Inverse[Mmat]·v − v ≡ 0` for ANY invertible `Mmat` — vacuous; does NOT confirm `Mmat` transcription. Classic "independent route reuses the same primitive" pitfall (= V.1 171 lesson). The V.1 scalarization (`Total[…^2]`) fixed list→scalar but NOT the tautology.
- **Fix (review-prescribed):** replace with a forward-transcription check — `Mmat . {drho,da,dcs,dZ}` vs the HAND-TYPED branch laws `{2*drho, drho+2*da, dZ+2*dcs−2*da, 5*(dcs−da)}` (notes §1), scalarized `Total[(…)^2]`. RHS must be hand-copied, NOT built from any matrix/Solve op.
- Directive findings_count=2 → both F1/F2 already-applied & PASS per review; only the regression R1 remains. SymPy has no matrix route — untouched.
- **Consult Qs:** confirm (1) the hand-typed `fwdLaws` vector matches notes §1 boxes; (2) `Mmat` rows (row3 `{0,-2,2,1}`=`dKq=dZ+2dcs−2da`, row4 `{0,-5,5,0}`=`dP=5(dcs−da)`) transcribe eq1..eq4 correctly; (3) replace vs supplement the round-trip (recommend replace).

### 175 — wall_normalized_load_shape  (review: 1 finding | needs Mathematica-feasibility consult)
- **R1 (transliteration):** Mathematica `Sigma_N` `dlog` block (wl:26-29 helper + wl:95-98) is a line-by-line port of the SymPy `dlog(diff(log()))` route → no independent second-engine protection. F3-step3 (an independent `dlogSeries` route) was BLOCKED in V.1 and waived as a policy mirror; **the review does NOT accept the waiver for this differential block.**
- Directive findings_count=3 → F1 (a simplify-commutes near-tautology) and F2 already resolved + PASS in the live scripts; **do NOT re-prescribe F1/F2** (re-touching F1 risks reintroducing the simplify-commutes trap). Only R1 survives.
- **Crux consult Qs:** (1) is the policy mirror genuinely SANCTIONED here (does the algebraic `N0−Λ²` route at wl:61 + Common-shape + Xi_load checks already give independent coverage?) OR does the differential slope need a real second route? (2) **Mathematica feasibility:** will `dlogSeries[e]:=Coefficient[Normal[Series[Log[e],{eps,0,1}]],eps]` on the exp-laden `(p/δ)/.subsHat/.subsEps` argument robustly land `===0` under FullSimplify? — the applier blocked it precisely because this can't be decided without RUNNING Mathematica. So 175's resolution likely needs a Codex/exec-backed feasibility check, with a mandatory escape clause (accept sanctioned mirror w/ written justification, NOT a forced route that silently changes a passing result).
- Anti-guards on rewrite: assert series-route-DIRECT vs SHAPE-target (`dlog[(λ²/k)/.subsEps]`), NOT series-vs-D[Log] on the same argument (that's a differentiation-method near-tautology). Keep `−kappa` symbolic.

## Execution plan (AFTER user GO — do not start yet)
1. **One Claude+Codex read-only consult** covering the cross-stage math questions, esp. **148-R2 (the live chain-rule divergence)**, **157 independence/defer**, **175 dlogSeries feasibility** (148-R1, 150 display, 166 forward-laws can ride the same consult). Record clean `redteam/codex_reviews/_consult_batch7.md`; delete bloated raw. Escalate to user ONLY if a fix changes the CONCEPTUAL nature (e.g. if 157's honest resolution is a relabel+defer that drops a claimed in-stage check).
2. **Clean audit agent REWRITES each directive** to encode its `codex_review` (drop moot/already-applied items; fix 148's 168→100 incl. claim manifest M1 + self-test; add 157-R3 `0<σ<1`; add 157 SymPy banner 138-139→155-156; set findings_count to the review's count). Orchestrator-review each drafted directive (batch 4/5 caught still-tautological & sign errors this way).
3. Per stage: `set-status NNN fixing` → `codex-invoke NNN directives/stage_NNN.md` (Codex applies+RUNS+iterates ≤600s) → orchestrator `exec-sympy NNN` then `exec-mathematica NNN` SEQUENTIAL + refresh `output/*.txt` → `capture-diff` → `render-verify-prompt NNN > tmp_prompts/verify_prompt_NNN.md` → clean verify agent writes `verifications/stage_NNN.md` → `set-status NNN verified`.
4. Seat-safe pairing (2 Mathematica seats): run 2 Codex in parallel under the flock, gate orchestrator `exec-mathematica` to AFTER each Codex pair; verify agents (0 seats) ∥ next Codex pair. Suggested pairs: **{166,150}** (clean/light) then **{148,157}** (heavy — 148 live bug, 157 independence). 175 in its own unit if not folded in.
5. Sync the 6 trackers, update handoff + memory, commit (only when user asks).

## NOT to do
- Do NOT touch `paper/**.tex` or `notes/stages/*.md` in the fix loop (scripts-only; published-paper edits are Codex-applied + Claude-reviewed in the deferred paper pass).
- Do NOT rewrite directives or consult Codex before the user's GO (this file is recon only).
- The 3 misfiled `notes/.../review/stage_{148,150,157}_review.md` (wrong-stage bodies) are out of scope — flag for a separate cleanup, do not fix in the fix loop.
