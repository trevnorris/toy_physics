# Stage Verification Coverage Baseline

This document is the stage-coverage control sheet for the PDE ledger archive.

Use it together with `STAGE_PROVENANCE_INDEX.md`:

- `STAGE_PROVENANCE_INDEX.md` holds the exact raw note and audit artifact paths.
- `STAGE_VERIFICATION_COVERAGE.md` summarizes the current verification surface,
  exposes the main gaps, and gives us a stable baseline for audit planning.

Snapshot date: `2026-06-03` (**batch VIII.1 close {243–253} — `Part VIII.1 — Relaxed branch, dynamic event chain, cold survival`, the FINAL forward first-pass batch under the v2 paper-grounded auditor WITH the dual-engine rule.** Forward first-pass for stages **243–253** (11 stages). All 11 reached `verified`; **`material_change: false` on all 11**; 0 stop-cold, 0 ultimately-blocked; 10 of 11 codex-invoke exited 0 on iteration 1, **one iter-2 (stage 248)**. **Cumulative coverage rises from 242/253 (95.7%) to 253 of 253 stages red-team verified (100%) — THE FIRST END-TO-END RED-TEAM PASS IS COMPLETE; the entire range 001–253 is now paper-aligned at v2 depth.** **THREE checkpoints this batch — 243, 248, 253 — ALL cleared the higher checkpoint bar (no rubber-stamp): all THREE EXISTING `.wl` were line-by-line transliterations of their `.py` and were RE-AUTHORED to genuinely independent routes** (243 = IBP closure / native `LinearSolve` / `TrigExpand` / `Series`-at-∞ for the leakage-work lane and non-rigid solve, hardcoded `expected*` residuals DELETED; 248 = the §II `Solve`-mirror replaced by a native SATISFACTION route — compiler closed forms verified to satisfy their defining energy equalities via substitution + `FullSimplify` + non-vacuity guard + positive-branch guard; 253 = native `D[Log[V[r]]]` + 5 `Solve` energy/force balances + regrouped threshold / Korringa / screening blocks; no checkpoint constant CHANGED — route-only rewrites land on the SAME values). **Dual-engine: 8 NEW independent-route `.wl` authored (244, 245, 246, 247, 249, 250, 251, 252) + 243, 248, 253's existing transliteration `.wl` RE-AUTHORED → 11 independent Mathematica audits, 0 sanctioned mirrors in VIII.1.** So the Coverage-Totals "Mathematica audits present" raw count rises by 8 to 242 (243/248/253 were already counted). **Dominant defect theme (continuing from VII.2) — the "variable-independence self-test trap"** (support-blindness / independence "verified" by differentiating w.r.t. an absent variable; vacuous), hit at 244-F1 and 245-F1; fixed via free-symbol containment + positive control, or a live-channel positive control. **PLUS tautological/round-trip checks** — 246 (σ_min-vs-own-piecewise), 247 (inverse round-trip + back-substituted λ_L), 249 (subtraction-by-linearity over blank placeholders), 251 (quintic X−X), 252 (gamma_safe_eq X−X), 253 (r_turn_phys expr−expr); 250's defect was a global claim tested only at a single sample point (fixed → global strict monotonicity via `Resolve[ForAll]`). **Stage 248 iter-2 (notable):** iter-1 prescribed a `Reduce`/`ToRules` route which Codex correctly BLOCKED (a Wolfram-version dead end — `ToRules` rejects the domain-predicate conjunction); the ORCHESTRATOR REFRAMED to the native SATISFACTION route above (a Claude math-coverage resolution, NON-conceptual, no paper edit), Codex applied iter-2, independent verify confirmed. **PAPER_MISALIGNMENT (5 typos across 4 stages, ALL notes-only, USER-RESOLVED 2026-06-03; the published paper cards were UNAFFECTED — they carry abstract forms):** NUMERICAL/coefficient TYPOS corrected to the verified SymPy scripts, each CROSS-ENGINE or INTERNALLY CORROBORATED — 247 (notes:406 Δ `210.17750000→142.17750000`; notes' own 9·16−1.35²=142.1775 + adjacent D0=3.76481862), 253 (notes:274 benchmark `187.23361317→119.23361317`; notes' own 65.45193926/0.5489386551=119.2336 + cross-engine `.wl`; AND notes:419 a_int `10.95423247→10.95423248` =4·K_turn), 244 (notes:366 `196√2→128√2`; script E0=16 structure), 248 (notes:506 `×168%→×100%`; notes' own 23.3128% + script ×100 — recurs the stale "168" at 148/232). **NOTE:** the audit agents initially MISATTRIBUTED 247/253 to the PUBLISHED cards but the orchestrator confirmed the cards are CLEAN (247 card=93 lines, 253 card=140 lines) — all 5 are NOTES-ONLY. **In-file stale-label `.wl` fixes that rode the fix loop (single-file, NOT a batch renumber):** 243 `.wl` `STAGE 226→243`, 248 `.wl` `STAGE 231→248`, 253 `.wl` banner→253. **NO NEW NUMBERING RESIDUAL FROM VIII.1** — all 11 notes/stages H1 titles are already canonical (Stage 243…253). The deferred PROJECT-WIDE stem-keyed renumber reconciliation (already logged: VII.1 219/221/222/223-title/227/228/229 + `.wl` banners 221/227; VII.2 notes-title drift 232/234/235/236) STILL STANDS for the post-253 dedicated cleanup pass; NEVER offset-sweep. See `PAPER_CLEANUP_TRACKER.md` P4-55, `MATHEMATICA_MIRROR_POLICY.md`. Previous VII.2 snapshot below.)

Snapshot date prior: `2026-06-03` (**batch VII.2 close {231–242} — third batch of Part VII under the v2 paper-grounded auditor WITH the dual-engine rule, theme "Rigid-mouth orbit-lock / branch-dressing / twin-support".** Forward first-pass for stages **231–242** (12 stages). All 12 reached `verified`; **`material_change: false` on all 12**; 0 stop-cold, 0 ultimately-blocked; 11 of 12 codex-invoke exited 0 on iteration 1, **one iter-2 (stage 238)**. **Cumulative coverage rises from 230/253 (90.9%) to 242 of 253 stages red-team verified (95.7%)** — the entire range 001–242 is now paper-aligned at v2 depth; stages 243–253 remain `pending`. **TWO checkpoints this batch — 239 and 242 — BOTH cleared the higher checkpoint bar (no rubber-stamp): both EXISTING `.wl` were line-by-line transliterations of their `.py` and were RE-AUTHORED to genuinely independent routes** (239 = forward Jacobian of the boxed dependent vector + native `PseudoInverse` left-inverse + `Reduce`/`Equivalent` orbit-lock, vs the `.py`'s backward hardcoded `SrmDep`; 242 = `Resolve[ForAll,Reals]` strict-inequality certificate + `D[]` on real closed forms vs the abstract-ζ device + `logDrift` total-log-differential vs `Exp[t·d]`, so the load-bearing twin-window strict inclusion `C_mix < Pi_tr < 2 C_mix` is now tested STRICTLY on both engines). **Dual-engine: 10 NEW independent-route `.wl` authored (231, 232, 233, 234, 235, 236, 237, 238, 240, 241) + 239 and 242's existing transliteration `.wl` RE-AUTHORED → 12 independent Mathematica audits, 0 sanctioned mirrors in VII.2.** So the Coverage-Totals "Mathematica audits present" raw count rises by 10 to 234 (239 and 242 were already counted). **Dominant defect theme — the "variable-independence self-test trap":** support-blindness / independence "verified" by differentiating an expression w.r.t. variables it never contains (vacuous), hit at 237-F2, 238-F1, 240-F1; fixed by negative-control + leak-detector + structural exclusion (237/238) or by extract-the-weight-from-the-variable-bearing-object (240, weights pulled from `Y_support` which carries `Omega_Q` via its pole). **Stage 238 iter-2 (notable):** F1 (support-blindness) and F4 (its `.wl`) were correctly BLOCKED on iter-1 — the notes define `M_tr=M_mix[1+ζ(1−ε)/(1−ζε)]` (§4) but give NO pre-reduction observable form where `M_tr` cancels, and inventing one was forbidden, so Codex refused to fabricate; the ORCHESTRATOR REFRAMED F1/F4 (Claude math-coverage resolution, NON-conceptual, no paper edit) to a faithful non-vacuous form — (a) negative control `∂_ζ M_tr≠0`; (b) leak detector `∂_ζ(Rtr·M_tr/M_mix)≠0`; (c) exclusion of {ζ,M_mix,M_tr} from the reduced observables — Codex applied iter-2, independent verify confirmed non-vacuity (cross-check: 237-F2 had been fixed the same way — live-channel negative control + `Not[FreeQ[#,Derivative]]` guard). **PAPER_MISALIGNMENT (3 stages, ALL notes-only, USER-RESOLVED 2026-06-03; the published paper cards/appendices were UNAFFECTED — they carry abstract forms):** NUMERICAL/coefficient TYPOS corrected to the verified SymPy scripts, each CROSS-ENGINE CORROBORATED by the new `.wl` independently computing the corrected value — 231 (line-98 `dF/dξ` numerator `240·δ²ξ→189·δ²ξ` and `189·ξ³→121·ξ³`; SymPy `Factor` + `.wl` M1 both give 189/121), 232 (lines 153/157 figure-of-merit prefactor `168→100`, the only value reproducing the notes' own decimals Ξ_χ≈5.5548e5/Ξ_J≈1.2664e5; M4 uses c=100 — mirrors the recurring stale "168" typo at 148 IV.5), 241 (line-577 `ϱ_WΛ` upper bound `193/369→125/369`, matching the notes' own printed decimal 0.338753; M7 computes `ϱ_WΛ|_{β=2/11}=125/369`). **In-file stale-label fixes that rode the fix loop (single-file, NOT a batch renumber):** 233 `.wl` comments (Stage 188/223/224→239/240/241), 239 `.wl` banner `STAGE 222→239` + `Stage221→Stage238` suffixes, 242 `.wl` banner `STAGE 225→242`. **RESIDUAL notes-title renumber drift DEFERRED (NOT fixed this batch — post-253 stem-keyed pass per the project numbering-drift policy; NEVER offset-sweep):** notes self-titles 232 "Stage 249", 234 "Stage 251", 235 "Stage 251/252/253", 236 "Stage 253"; canonical (paper-card/script/MANIFEST) = 232/234/235/236 is ground truth. See `PAPER_CLEANUP_TRACKER.md` P4-54 (continuation of P4-53), `MATHEMATICA_MIRROR_POLICY.md`. Previous VII.1 snapshot below.)

Snapshot date prior: `2026-06-02` (**batch VII.1 close {219–230} — second batch of Part VII under the v2 paper-grounded auditor WITH the dual-engine rule, theme "Mixed-bundle / resonance / branch-packet".** Forward first-pass for stages **219–230** (12 stages). All 12 reached `verified`; **`material_change: false` on all 12**; 0 stop-cold, 0 blocked, all 12 codex-invoke exited 0 on iteration 1 (no iter-2 reworks). **Cumulative coverage rises from 218/253 (86.2%) to 230 of 253 stages red-team verified (90.9%)** — the entire range 001–230 is now paper-aligned at v2 depth; stages 231–253 remain `pending`. **Checkpoint 221 cleared the higher checkpoint bar (no rubber-stamp): its EXISTING transliteration `.wl` was RE-AUTHORED to an independent route (native `D[QPi/DeltaPi,portPi]` derivative, `Residue`, `ComplexExpand`, uncollapsed Breit–Wigner); F1 tautological survival round-trips de-tautologized; F2 deliverable #9 (linear survival window) print-only/tautological → genuinely covered in BOTH engines; F4 notes renumber.** One sanctioned Codex deviation (221-F3): used the native `D[QPi/DeltaPi,portPi]` instead of the directive's leading-minus form, reconciling the Stage-220 identity `∂_Π D_Π = −N` — verified correct. **Dual-engine: 11 NEW independent-route `.wl` authored (219, 220, 222, 223, 224, 225, 226, 227, 228, 229, 230) + 221's existing `.wl` re-authored → 12 independent Mathematica audits, 0 sanctioned mirrors in VII.1.** Every new `.wl` was confirmed independent by a clean verify agent (native primitives via a different decomposition than the SymPy `.py`; e.g. 219 structural family extraction via Collect/CoefficientList, 220 Laurent-support CoefficientRules, 226 Orthogonalize projector not QR, 228 NSolve+Series+implicit-function slopes, 229/230 Resolve[ForAll] universal-quantifier proofs). So the Coverage-Totals "Mathematica audits present" raw count rises by 11 to 224 (221 was already counted). **Script-side findings fixed besides the missing `.wl`:** 220 F2 insufficient → symbolic P_abs perfect-square assert; 221 F1/F2 (above); 223 F2 tautological circular compat-surface → `sp.solve(Eq(K_norm,K_pole),P0_target)`; 224 F2 hardcoded budgets-vs-themselves → defining-relation checks tied to the ceiling; 225 F2 tautological `0==0` one-pole → tied to the one-pole constraint `D4=−3·D0·u2²` (+ a negative-control `expectNonZero`); 230 F2 tautological onset round-trip → `sp.solve(Eq(onset,R_star),delta)`. **PAPER_MISALIGNMENT (7 stages, ALL notes-only, USER-RESOLVED 2026-06-02; the published paper cards/appendices were UNAFFECTED — they carry abstract forms):** NUMERICAL TYPOS corrected to the verified SymPy scripts, each CROSS-ENGINE CORROBORATED by the new `.wl` independently computing the corrected value — 222 (λ_W=0.2 upper-wall R_Q `213.483858657863`→`145.483858657863`), 223 (λ_W=0.2 wall R_Q `206.814136942081`/`205.502546600713`→`138.814136942081`/`137.502546600713`), 227 (i=h rigidity det factor `251+215π²`→`200+147π²`), 228 (δ_1 coeff `247π²/(98π²−25)`→`196π²/(98π²−25)` AND reduced-det `247(251+215π²)`→`196(200+147π²)`), 229 (crossover-cubic leading coeff `189ξ³`→`121ξ³`) — the slips were a systematic +68/+51 additive family; plus RENUMBER (notes labels→canonical) 221 (238→221, 237→220) and 225 (240/241/242→223/224/225 + a filename token) as formal findings, and 220/224/226/230 informational drift renumbered by Codex. **Orchestrator note:** 230's codex log was captured EMPTY (a logging anomaly, NOT a stall — exit 0 recorded, all `## Applied` blocks present, all artifacts built, independently re-confirmed + verified); the `.wl` naming bug from VI.1 did NOT recur (all 12 directive `.wl` targets carried the `_mathematica_audit.wl` suffix — an explicit pre-invoke guard worked). See `PAPER_CLEANUP_TRACKER.md` P4-53, `MATHEMATICA_MIRROR_POLICY.md`. Previous VI.1 snapshot below.)

Snapshot date prior: `2026-06-02` (**batch VI.1 close {201–218} — first batch of Part VI under the v2 paper-grounded auditor WITH the dual-engine rule, theme "Explicit realization, scalar slice, ray ranking".** Forward first-pass for stages **201–218** (18 stages). All 18 reached `verified`; **`material_change: false` on all 18**; 0 stop-cold, 0 blocked, 0 needs_rework left open. **Cumulative coverage rises from 200/253 (79.1%) to 218 of 253 stages red-team verified (86.2%)** — the entire range 001–218 is now paper-aligned at v2 depth; stages 219–253 remain `pending`. **Checkpoints 203 and 218 BOTH cleared the higher checkpoint bar (no rubber-stamp).** **Dual-engine: every stage has an independent second engine; 0 sanctioned mirrors.** 16 stages had NO `.wl` and got a NEW independent-route `.wl` (201, 202, 204–217), so the Coverage-Totals "Mathematica audits present" raw count rises by 16 to 213 (218 and 203 were already counted); checkpoint 218's pre-existing transliteration `.wl` was RE-AUTHORED to an independent route (M4 witness counts even differ across engines, 256/192/65+63 in the `.wl` vs 192/192/64+64 in the `.py`), and 203's `.wl` got a strengthened independent composition check. **Load-bearing corrected constant: the per-envelope lifted Bézout bound `162 = 3⁴·2` (stages 217/218)** — arithmetically forced; the PUBLISHED card (`stage_217.tex`), the Part-VI appendix (`stage_appendix_part06.tex`), and the 217/218 notes had carried WRONG typos (179 and 230) which were corrected to 162 (the SCRIPT was always correct; no script value moved). **One iter-2 timeout rework** (202's first `.wl` symbolic-`Solve` of a transcendental log equation TIMED OUT at the 600s cap → reformulated to a `LinearSolve` of the log-linearized monomial-match system) and one iter-2 forward-ref follow-up (212 R2 renumber missed 247→213). Four orchestrator catches (three `.wl` target paths dropped the `_mathematica_audit` suffix — 201/210 renamed, 212/216 directives fixed pre-build; plus the 202 timeout and 212 forward-ref). Paper/notes edits were Codex-applied + Claude-reviewed this batch (217/218 `230→162`, 217 card/appendix `179→162`, 212 budget typo + stale stage-label renumbers, 214 `218→150`); see `PAPER_CLEANUP_TRACKER.md` P4-52. Previous retro-sweep snapshot below.)

Snapshot date prior: `2026-06-01` (**retro-sweep {121, 122, 123} — DUAL-ENGINE GAP CLOSED for already-verified non-status-only stages.** These 3 stages were already `verified` in batch IV.3 but were SymPy-ONLY. Under the dual-engine rule (a `.wl` is REQUIRED wherever Mathematica CAN independently verify — test = "is it possible," not "is it necessary"), each was retrofitted with a NEW independent-route `.wl` (Codex designed + wrote; Claude reviewed audit + verify). **3/3 verified, `material_change: false` on all 3** — the SymPy `.py` reference engines were UNCHANGED and paper/notes UNTOUCHED; 0 transliterations accepted, 0 iteration-2 reworks, 0 blocked, 0 stop-cold. **121** (geometric_r_selection) 6 PASS — `r_geom` via `Solve` of the Stage-99 length law + positive-branch `Select`, `r_F1` surd via `FullSimplify`/`RootReduce`; **122** (mouth_source_compensation_test) 9 PASS — `g_±` DERIVED by `Solve`-ing the compensation quadratic + branch select, the traction-ratio de-tautology preserved cross-engine with `cStage` a FREE POSITIVE SYMBOL collapsing only when g_nat→1 is substituted last; **123** (parent_normalized_branch_values) 6 PASS — `Xi_v`/`Xi_T` LAWS from `Reduce`-based inversions, the un-squared λ NEGATIVE sign preserved so `Xi_v(F1) ≈ −1.01675633282526` (negative, not the +1.0168 mis-pass). **Cumulative coverage UNCHANGED at 200/253 verified (79.1%)** — a retro-sweep adds 2nd engines to already-`verified` stages, it does NOT add new verified stages; 121/122/123 were the ONLY 3 already-`verified` non-status-only stages still missing a `.wl`, so that dual-engine gap is now CLOSED (the 11 status-only single-engine stages compute nothing a `.wl` could check and legitimately remain single-engine). No checkpoints in range, no new constants, no EM-projection change. 3 new `.wl` files now exist on disk for stages 121/122/123 (the Coverage-Totals "Mathematica audits present" baseline figure was stale at 165 — not maintained since before V.2 — and has now been corrected to the raw count 197 in Coverage-Totals, matching the SymPy raw-count convention). See `MATHEMATICA_MIRROR_POLICY.md` retro-sweep entry. Previous V.3 snapshot below.)

Snapshot date prior: `2026-06-01` (**batch V.3 close — second batch of the resumed first pass; THE DUAL-ENGINE RULE CORRECTION batch.** First-pass paper-grounded audit under v2 prompt for stages **188–200** (`Part V.3 — Branch observables, isotropic target, home stretch`; checkpoint 200). All 13 reached `verified`; **7 clean** (190, 192, 194, 196, 197, 198, 199), **6 dirty** (188, 189, 191, 193, 195, 200). 13 findings, all resolved, 0 blocked, 0 stop-cold, **`material_change: false` on all 13** (new-second-engine / de-tautologized / demoted-to-definition / banner-relabel only — no derived value, constant, identity target, or paper number moved, so no `upstream_stale` propagation). **200 of 253 stages now red-team verified (79.1%)** (the MANIFEST count: 187 in range 001–187 + 13 in range 188–200; 201–253 all `pending`); range 001–200 paper-aligned at v2 depth. **The defining event: a standing-drift correction to the dual-engine rule (user-clarified 2026-06-01)** — a Mathematica `.wl` is REQUIRED wherever Mathematica CAN independently verify (test = "is it possible," not "is it necessary"). Before V.3, stages 188–199 were SymPy-ONLY; consequence: **12 NEW independent-route `.wl` (188–199) were created by Codex, ALL genuine independent routes (0 transliterations accepted)**, and the checkpoint 200's pre-existing `.wl` (caught as a transliteration) was DE-TRANSLITERATED. So all 13 V.3 stages are now genuine dual-engine. **0 sanctioned mirrors.** Labor split: Claude reviews (audit + verify); Codex writes ALL script code. **One iteration-2 rework** (189 Section II, orchestrator-review catch + no-rubber-stamp — the re-pointed selected-branch identity was still a back-definition ≡ 0 for any input; demoted to a printed definition, with the genuine direct-slope bridge `δln T_A²=ε·λ_A·Ξ₁` carrying an INDEPENDENT `Ξ₁`; consult `redteam/codex_reviews/_consult_V3.md`, session 019e843e, CONCUR full). Two stale stage-number banner relabels (189 "172"→"189"; 191 "174"→"191"; string-only, settled canonical convention, NOT escalated). See `MATHEMATICA_MIRROR_POLICY.md` V.3 entry for the full dual-engine retrofit. Previous V.2 snapshot below.)

Snapshot date prior: `2026-05-30` (**batch V.2 close — first batch of the RESUMED FIRST PASS** after the IV.x/V.1 orchestrator-direct integrity remediation (batches 1–8) closed; run under the RESTORED Codex-as-fix-applier contract. First-pass paper-grounded audit under v2 prompt for stages **176–187** (`Part V.2 — Load shape, transfer shape, coherent slippage`). All 12 reached `verified`; 2 clean (179, 180), 10 dirty. **185 is the only checkpoint** (verified at the higher bar). ~24 findings, all resolved, 0 blocked, 0 stop-cold, **0 paper_misalignment**, **`material_change: false` on all 12** (additive/strengthening checks only — no derived value, constant, identity target, or paper number moved, so no `upstream_stale` propagation). **199 of 253 stages now red-team verified (78.7%)** (was 187 after V.1); range 001–187 paper-aligned at v2 depth. Two dominant patterns mirror V.1: **`mathematica_transliteration` (7 stages — 177/178/181/182/184/186/187, plus 176-F2 comment-only), ALL given genuine independent routes, ZERO accepted as sanctioned mirrors** (the converse of V.1's 2); and the same stale `STAGE N−17` script-banner cluster (176→159 … 187→170), fixed script-side in the fix loop (the 2 clean stages 179/180 carry residual cosmetic banner drift, non-blocking). **Two iteration-2 reworks** (verify-wave catches, no rubber-stamp): 181 F1 (vacuous ζ-free Mathematica check → routed through ζ-bearing `rTargetLoaded` + `FreeQ` guard) and the checkpoint **185 F1** (coefficient not load-bearing → independent slippage-law `Theta₁`/`Xi₁` reconstruction; 2nd consult `019e77e6`). Consult: `redteam/codex_reviews/_consult_V2.md` (session 019e77af + iter-2 019e77e6). Previous batch-8 snapshot below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 8 — FINAL** — re-verification of ONE stage whose first-pass fix had been applied orchestrator-direct (Codex bypassed) and has now been Codex-reconciled + re-verified. Stage **175** (V.1-range `wall_normalized_load_shape`) REMEDIATION-`verified` (2026-05-29); both engines exit 0/0; this status **supersedes the earlier orchestrator-direct "verified"** it carried from its first-pass V.1 batch. No new stages added to the verified count (175 was already counted from V.1); this is a verification-hardening / integrity re-pass. **175 is the LAST of the 29 findings stages, so with this close ALL 29 findings stages are now remediated/`verified`** — the only remaining work is the planned full end-to-end second pass. **`material_change: false`** — the edit only ADDS a corroborating independent check; no derived value, constant, identity target, or paper number moved. **NO checkpoints in the batch-8 range** (V.1 has no checkpoints; nearest checkpoints are 096 and 105). Finding closed: **175 R1 mathematica_transliteration** — the Mathematica `Sigma_N` differential block was a line-by-line transliteration of the SymPy block (both extracted the first-order log-slope via the SAME `dlog = D[Log[.],eps]/.eps->0` primitive), so the differential SLOPE identity was singly-routed across engines; in V.1 this was ACCEPTED as a policy mirror (old "175 F3-step3" `dlogSeries` waiver), but batch 8 REJECTED that waiver and added a genuine INDEPENDENT Mathematica Series-route check via a native `dlogSeries[expr_] := Coefficient[Normal[Series[Log[expr], {eps, 0, 1}]], eps]` extractor and exactly ONE new `expectZero["Sigma_N - dln(Lambda^2/K) [series route]", 2*dlogSeries[exprPoverDeltaPhys] - kappa - dlogSeries[(lambda^2/k)/.subsEps]]` (series route compares series-route DIRECT vs the SHAPE target — NOT `dlogSeries` vs `dlog` on the same arg, `-kappa` kept symbolic); the `[series route]` check now sits ALONGSIDE the existing algebraic `N0 - Lambda^2`, common-shape, and `Xi_load` checks. Structurally distinct from SymPy's `sp.diff(sp.log(...))` route; it landed `=== 0`, so the mandatory escape clause (accept the sanctioned mirror with written justification) was AVAILABLE but NOT triggered — independence ACHIEVED, not waived (the converse of the V.1 F3-step3 disposition). The existing `dlog` line was LEFT UNTOUCHED as corroboration, and the SymPy `.py` was LEFT UNTOUCHED as the reference engine: **SymPy 13 PASS (UNCHANGED); Mathematica 13→14 PASS (+1 = the new `[series route]` line)**. Direction settled by ONE Claude+Codex read-only consult (4/4 CONCUR; none conceptual; none escalated to the user), recorded at `redteam/codex_reviews/_consult_batch8.md`; Codex applied iter=1, deviation: none. NO new paper-cleanup items (no `.tex`/paper number touched); a misfiled `notes/stages/review/stage_175_review.md` (H1 "Stage 175" but body points at pre-renumber `stage022_grouped_p2_normalization_bridge` files) was flagged for separate orchestrator/notes repair, same class as the batch-7 trio (logged in `PAPER_CLEANUP_TRACKER.md`). See `redteam/verifications/stage_175.md`. Previous batch-7 snapshot below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 7** — re-verification of four stages whose first-pass fixes had been applied orchestrator-direct (Codex bypassed) and have now been Codex-reconciled + re-verified. Stages **148** (IV.5-range representative_positive_families), **150** (IV.5-range full_profile_residual), **157** (IV.6-range core_mouth_coevolution_status), **166** (V.1-range bundle_inversion_four_drifts) all REMEDIATION-`verified` (2026-05-29); both engines exit 0/0 on every stage; this status **supersedes the earlier orchestrator-direct "verified"** these four carried from their first-pass IV.5/IV.6/V.1 batches. No new stages added to the verified count (all four were already counted from their first-pass batches); this is a verification-hardening / integrity re-pass. **`material_change: false` on 150/157/166; `material_change: true` on 148 — with ZERO downstream propagation** (the defective Mathematica `dT` was corrected to match the already-correct paper/SymPy `dT`, which did NOT move, so no downstream stage is stale; SymPy-derived numerics were always correct). **NO checkpoints in the batch-7 range** (`is_checkpoint: false` / `is_status_only_candidate: false` confirmed in MANIFEST for 148/150/157/166; IV.5/IV.6/V.1 have no checkpoints; nearest checkpoints are 096 and 105). Findings closed per stage: **148** (2 findings) — **F2 insufficient_verification (the LIVE bug)**: the Mathematica first-order traction-shift route (`dSigmaOfDeltas`/`dTOfDeltas`) routed `dG` only through `dPi=−dG/gPrimeStar` and SILENTLY DROPPED the S-follows-Π chain term that SymPy's `AT` carries, so Mathematica computed a WRONG `dTU=0.4976…` while SymPy was correct at `0.5087…` and NOTHING asserted cross-engine agreement (both passed while disagreeing); fixed by DELETING that block and deriving `aT`/`bT` via Mathematica's OWN `D[]` autodiff of `Tm[p]=Sqrt[(9/20)(p/(1−sFormula/4))]` along the S(Π) curve (regenerating the dropped chain term), with BOTH engines now anchoring `A_T`/`B_T` to the published paper literals `A_T=−4.27263956256927`, `B_T=0.134875005736706` (`paper/appendices/stage_appendix_part04.tex:846,848`) and `dTU`/`dTD` agreeing across engines to ~16 digits (no baked SymPy literal in the `.wl`). **F1 paper_misalignment/insufficient_assertion**: the same-source ξ_* bridge check raised to an EXACT symbolic zero on the rF1-forced `100π²` radical (`12·(37/20)²=4107/100`); the stale `168π²` was purged as a directive-doc typo — the scripts were always correct. **150** (1 finding, insufficient_verification, DISPLAY-only) — both transcripts now print the compact slope `S_q(Π)=Aq·k−Cq·Π` (coefficient symbols preserved, built from FREE placeholders then `.subs`/`/.` to the concrete definitions feeding the load-bearing `T_q'(0)−S_q==0` assert), not the expanded rational; the source slope was already correct/committed. **157** (3 findings) — F1/F2 de-tautologized the canonical-even check in BOTH engines: SymPy's duplicate re-solve and Mathematica's mirrored literal 9/72/5 system were replaced by a parallel determinant non-degeneracy assertion `det([[1,−9σ],[5,−72σ]])=−27σ≠0` (genuine fail mode on a mistyped coefficient / rank loss); F3 added the physical branch assumption `0<σ<1` (Mathematica). The SymPy docstring item 6 (was overclaiming "tangent motion kills delta C") was corrected to match the published card's EXISTING deferral of the deviation-to-normalization map to Stage 158, and a carry-forward banner (SymPy "138-139"→"155-156") was fixed. **ESCALATION NOTE:** the consult flagged a CONDITIONAL conceptual-escalate; the orchestrator RESOLVED it AGAINST escalation because the paper card (`stage_157.tex`) is tagged Open and already defers the map — so this is a how-it's-checked/labeling fix, not a conceptual change; **paper card NOT edited**. **166** (1 finding, tautological_check) — replaced the vacuous Mathematica matrix "round-trip" `Total[(Mmat·Inverse[Mmat]·v−v)^2]` (an X−X self-cancellation true for any invertible `Mmat`) with a genuine forward-transcription check `Total[(Mmat.{drho,da,dcs,dZ}−fwdLaws)^2]==0`, `fwdLaws={2 drho, drho+2 da, dZ+2 dcs−2 da, 5(dcs−da)}` HAND-TYPED from the notes §1 boxed laws (not built from `Mmat`/`Solve`/`Inverse`) — a wrong `Mmat` coefficient now fails. Directions settled by ONE Claude+Codex read-only consult (5/6 unconditional CONCUR + 1 conditional conceptual-escalate on 157 RESOLVED against escalation; none disputed; none escalated to the user), recorded at `redteam/codex_reviews/_consult_batch7.md`; NO new paper-cleanup items (the paper `A_T`/`B_T` are correct, scripts fixed to match), and THREE misfiled notes-review artifacts (`stage_{148,150,157}_review.md` bodies point at pre-renumber `stage029`/`stage031`/`stage038` source files) flagged for separate orchestrator/notes repair, logged in `PAPER_CLEANUP_TRACKER.md` P4-48. See `redteam/verifications/stage_{148,150,157,166}.md`. Previous batch-6 snapshot below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 6** — re-verification of four IV.5-range stages whose first-pass fixes had been applied orchestrator-direct (Codex bypassed) and have now been Codex-reconciled + re-verified. Stages **143** (exp-remainder positivity), **144** (threshold-Pi_* independence), **146** (affine-law residual), **147** (rigidity-kernel projection) all REMEDIATION-`verified` (2026-05-29); both engines exit 0/0 on every stage; this status **supersedes the earlier orchestrator-direct "verified"** these four carried from their first-pass IV.5 batch. No new stages added to the verified count (all four were already counted from IV.5); this is a verification-hardening / integrity re-pass. **`material_change: false` on every stage**. **NO checkpoints in the batch-6 range** (IV.5 has no checkpoints; nearest checkpoints are 096 and 105). Findings closed per stage: **143** (2 findings, insufficient_verification) — replaced a cubic-Taylor-coefficient-only "positivity" gate (which passed for a wrong remainder like `Pi^3/6 − Pi^4`) with a GENUINE GLOBAL positivity proof of `exp(Pi)−1−Pi−Pi^2/2 > 0` on `Pi>0` (SymPy Taylor-remainder chain `R(0)=R'(0)=R''(0)=0, R'''(Pi)=exp(Pi)>0`; Mathematica primary `Reduce[...>0,piM,Reals]`→True with the same Taylor backing). **144** (1 finding, transliteration) — replaced a line-by-line SymPy-mirror `Pi_*` block with an INDEPENDENT cleared-denominator bracketed `FindRoot` on `gThresholdResidual[p]=2p(2p e^p+pi)−gMinus(4p^2+pi^2)(e^p−1)` (sign kept as `(e^p−1)`, NOT the §131-F3 `(1−e^p)` trap) + residual-near-zero witness, anchored to the OWNING stage-131 `Pi_*` value `1.50882951349315558300555075595`; MIRROR-policy `needs_user_resolution` cleared (consult Q2: a load-bearing transcendental root is NOT a sanctioned mirror → independent route, precedent §131-F3/§142). **146** (3 findings) — R1 SymPy `nsolve(prec=50)` + tightened affine-residual assert from `1e-15` to raw `<1e-25` (~1e-51); R2 Mathematica removed `Chop[...,10^-6]` and asserts raw `<10^-25` (~1e-58/exact 0) via the directive-authorized SYMBOLIC ENDPOINT-INTEGRAL fallback (sanctioned, genuinely hits `<1e-25`); R3 banner `FINITE-CORRECTION EXPANSION`→`FIRST-ORDER EXPANSION` (`.wl` script fix matching the paper-card title). **CONSULT caveat (146):** the affine residual collapses to `(1−eps)(gPi(Pi_*)−gMinus)` and `gMinus` also feeds `Pi_star`, so this check is NOT an independent guard against `gMinus` (documented scope limit — tests intercept-vs-direct-integral + kernels/source). **147** (5 findings, the heavy one) — R1/R5a CAS-autodiff of `T_m(Pi)=sqrt((9/20)Pi/(1−Sformula/4))` replacing a hand-typed `A_T` chain rule; R2/R3 projection identity `∫W_*(σ−Σ_*)=A_T(ḡ_σ−g_*)+B_T(S̄_σ−S_*)` [σ=2x] + full-symbolic x-independence check + (consult-Q5-ADDED) source-centering assertion `∫Σ_*W_*=0` (the projection identity alone is blind to the centering constants; this orthogonality check tests them); R4/R5b source-moment quadrature integrals vs the closed forms replacing `g_*`/`S_*` resubstitution; R5 independent Wolfram `D[]`/`NIntegrate` (not a SymPy port); R6 (consult Q6) the `|A_T|/B_T` ratio `31.6785` was MISLABELED paper-quoted → relabeled as the script's own computed cross-check (the genuine paper literals `A_T=−4.27263956256927`, `B_T=0.134875005736706` are sourced from `paper/appendices/stage_appendix_part04.tex:846,848`). Directions settled by ONE Claude+Codex read-only consult (4 CONCUR / 2 DISPUTE-resolved, none conceptual, none escalated), recorded at `redteam/codex_reviews/_consult_batch6.md`; NO new paper-cleanup items, logged in `PAPER_CLEANUP_TRACKER.md` P4-47. See `redteam/verifications/stage_{143,144,146,147}.md`. Previous batch-5 snapshot below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 5** — re-verification of four Family-1 mouth-gain/branch stages whose first-pass fixes had been applied orchestrator-direct (Codex bypassed) and have now been Codex-reconciled + re-verified. Stages **134, 137** (IV.4-range canonical_mouth_gain / schur_reduction) and **139, 142** (IV.5-range family1_actual_mouth_gains / mouth_gain_susceptibility) all REMEDIATION-`verified` (2026-05-29); both engines exit 0/0 on every stage; this status **supersedes the earlier orchestrator-direct "verified"** these four carried from their first-pass IV.4/IV.5 batches. No new stages added to the verified count (all four were already counted from IV.4/IV.5); this is a verification-hardening / integrity re-pass. **`material_change: false` on every stage**. **NO checkpoints in the batch-5 range** (`is_checkpoint: false` confirmed in MANIFEST for 134/137/139/142; IV.4/IV.5 have no checkpoints). New check types added per stage: **134** (canonical_mouth_gain) R1 tautological_check — REMOVED the X−X "canonical gain line" SymPy assert (`intercept = Pi_*` vs the literal that defined `Pi_*`), replaced with a provenance comment deferring outlet/gain-line consistency to Stage 135 (the paper-card downgrade itself was done in IV.4); kept the corroborated shell-limit + 3 S_q spot-checks + S_q(Π_*) checks; `.wl` gain-line block relabeled "printed only; not asserted". **137** (schur_reduction) R1/R2/R3 — added an INDEPENDENT **matrix-Schur reconstruction** of `rho_c, sigma_c` from the physical core matrix `M_core=[[K_s,λ],[λ,−K_q·D]]`, `v=(g_s,g_q)`, `δ=vᵀM_core⁻¹v` (SymPy `apart`+`limit`; Mathematica `Inverse[mCore]`+`Normal[Series]`, mirroring the verified OWNER stage 114); replaced the X−X static-limit check with a matrix-route-vs-reduced-envelope identity, and the `S_q=0` outlet check (which deleted M_q) with a nonzero-S_q mixed-channel check `M_q·S_q == −L·σ_c_schur·S_q/Θ` + sign non-vacuity guard — **9 PASS / 0 FAIL**. **139** (family1_actual_mouth_gains) R1/R2/R3/R4 — R1 replaced tautological outlet asserts with an independent `S_q(Π_*)` recompute from the Stage-134 closed-form kernel; R2 relabeled `R_q^comp=1/4` as definitional-consistency + added a **branch-discrimination value anchor** `g_-^F1=0.758035078944662826919680890414` (discriminates the lower from the upper branch), made the Mathematica `gMinus` the DIRECT closed form (sanctioned mirror); R3 (paper_misalignment) resolved no-script-fix (susceptibility closure is a Stage-140 deliverable; existing P4-42); R4 (provenance) resolved scripts-already-correct. **142** (mouth_gain_susceptibility) R1/R2 — R1 relabeled the self-solved `R_q(Π_*)=1/4` as solver-consistency + added a non-tautological **cross-stage independent-Π_* anchor** evaluating `R_q` at STAGE 131's independently-derived `Π_*=1.50882951349315558300555075595` (cleared-denominator FindRoot, batch-4-verified — NOT 142's own nsolve); R2 replaced a `Normal[Series]` self-comparison with an INDEPENDENT **projection integral** `∫₀¹ σ_Π(z)·cos(πz/2) dz` (σ_Π = Stage-129 source law) closing to the hardcoded `gPi` with symbolic-zero residual; 5 external-decimal anchors retained; tolerance 1e-15 kept — **11 PASS / 0 FAIL**. Directions settled by ONE Claude+Codex read-only consult (7/7 CONCUR, none conceptual, none escalated), recorded at `redteam/codex_reviews/_consult_batch5.md`; NO new paper-cleanup items, logged in `PAPER_CLEANUP_TRACKER.md` P4-46. See `redteam/verifications/stage_{134,137,139,142}.md`. Previous batch-4 snapshot below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 4** — re-verification of four stages whose first-pass fixes had been applied orchestrator-direct (Codex bypassed) and have now been Codex-reconciled + re-verified. Stages **125, 126** (IV.3-range) and **130, 131** (IV.4-range) all REMEDIATION-`verified` (2026-05-29); both engines exit 0/0 on every stage; this status **supersedes the earlier orchestrator-direct "verified"** these four carried from their first-pass IV.3/IV.4 batches. No new stages added to the verified count (all four were already counted from IV.3/IV.4); this is a verification-hardening / integrity re-pass. **`material_change: false` on every stage**. Fixed via a 2-parallel Codex run enabled by a user-approved `redteam.sh` MANIFEST flock. **NO checkpoints in the batch-4 range** (IV.3/IV.4 have no checkpoints). Findings closed + PASS-line counts per stage: **125** (positive_source_theorem) F1 insufficient_verification — strengthened the peaked-source check from a weak smallness test `abs(g_a(100))<0.05` (which accepted small NEGATIVE moments) to GENUINE paper-bound range checks `0 ≤ g ≤ 1` (the `≥0` line, NOT wrapped in abs, fails on a sign error in the closed-form moment), retaining `<1/20` as a subordinate trend witness; Mathematica mirrored with two range bounds at the a=100 proxy (deviation: extra precision + `Re[]`/`Chop` on the incomplete-gamma closed form, benign) — **SymPy 3 new check lines / Mathematica 15 PASS**. **126** (positive_source_families) F1 insufficient_verification — replaced endpoint/corner-only positivity SAMPLING with a GLOBAL box-positivity proof (exact affine-in-ξ check `d²σ_ξ/dξ²==0` PLUS both ξ-endpoint slices nonneg over `[0,L]` ⟹ nonneg on the (z,ξ) box; Mathematica uses `Resolve[ForAll[...]]` with the directive-authorized `lM→1` scale-covariant quantifier; an interior-curved perturbation fails the affineness check, non-tautological) — **SymPy +3 lines / Mathematica 7 PASS**. **130** (mouth_bias_map) F1 insufficient_verification — replaced the finite 6-point monotonicity sweep with the notes-§2 FKG/Chebyshev symmetrized-covariance GLOBAL certificate proving `dg_Π/dΠ>0 ∀Π>0` and Π_* uniqueness (density normalization, symmetrized double-integral identity, integrand sign `f'(z)=−(π/2L)sin(πz/2L)<0` on (0,L), consistency `dg/dΠ=−Cov/L`, uniqueness bracket `2/π<g_-<1`; the directive-sanctioned optional SymPy `reduce_inequalities` probe was removed as non-load-bearing, the load-bearing sign arg intact; the double integral simplified under the 600s cap, no fallback) + F2 (boxed `g_Π=2Π(2Πe^Π+π)/((4Π²+π²)(e^Π−1))` closed form) pre-existing reconciled-PASSED, unchanged — **Mathematica 11 PASS / SymPy print+assert guards**. **131** (parent_mouth_threshold) F1 tautological_check — the parent-threshold-identity-at-Π_* check was X−X (purely DEFINITIONAL); per consult option (ii) the Anchor-3 block was DROPPED entirely in both engines + F2 insufficient_verification — replaced weak `g(2Π_*)≠g_-` with 3 branch-discrimination checks (lower-branch membership, singular-branch exclusion `g_nat−g(Π_*)≈Δg_-=0.241964921055337`, upper-branch exclusion `|g(Π_*)−g_+^{F1}|>1` with `g_+^{F1}≈2.79795199200529`) + F3 transliteration — **NOT accepted as a policy-mirror**, added an INDEPENDENT Mathematica route for Π_* via the cleared-denominator residual solved by bracketing `FindRoot[{piM,1.4,1.6}]` (per the IV.5 139/143/144 transcendental-numerical-root precedent); reconciled-PASSED Anchor1 (Π_*=1.50882951349316), Anchor2 (g'(Π_*)=0.0714453558083195), F4 (g_-^{F1} closed form vs 0.758035078944663), banner STAGE 131 unchanged — **SymPy 6 / Mathematica 6 PASS**. The resolution directions came from a Claude+Codex consult (Codex base `019e7594`, all CONCUR, none conceptual), recorded at `redteam/codex_reviews/_consult_batch4.md`; paper-card cross-reference / informational-only status logged in `PAPER_CLEANUP_TRACKER.md` P4-45 (NO new paper-cleanup action items — note specifically that stage 130's paper-card claim "monotonicity, and unique Family-1 point Π_*" is now PROPERLY verified globally). See `redteam/verifications/stage_{125,126,130,131}.md`. Previous batch-3 snapshot below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 3** — re-verification of four IV.3-range stages whose first-pass IV.3 fixes had been **tainted-applied** (edits committed in the IV.3 pass `b4e02d8` while Codex was bypassed) and have now been Codex-reconciled + re-verified. Stages **117, 118, 119, 122** all REMEDIATION-`verified` (2026-05-29); both engines exit 0/0 on every scripted stage; this status **supersedes the earlier tainted "verified"** these four carried from their first-pass IV.3 batch. No new stages added to the verified count (all four were already counted from IV.3); this is a verification-hardening / integrity re-pass. **`material_change: false` on every stage**. **NO checkpoints in the batch-3 range** (IV.3 has no checkpoints). Findings closed + PASS-line counts per stage: **119** (parent_balance) F1 tautological_check de-tautologized (tube-length law `rc → rhat²` link ties L_W to the family `rhat`) + F2 insufficient_verification (`T_m (±)` branch match vs notes §5 closed form; `.wl` uses a `stripCE` ConditionalExpression-stripper) — **8 checks/PASS per engine**, genuine independence confirmed (not X−X). **118** (parent_core) F1 paper_misalignment (λ SIGN) → resolved **direction (a) = MINUS** (consistent with the script's own independently-derived section IV `sq_coeff = −q_* ϱ_s v0 A_q`, the notes' three boxed minus signs, and downstream stage 123's UN-squared λ; NOT conceptual — paper card states no λ sign, so no paper-cleanup item) + F2 insufficient_verification added three asserts (`K_q closed form`, `g_s closed form`, `lambda from bilinear`) — **12 sympy / 14 mathematica PASS**. **CAVEAT (118):** verify agent flagged `K_q closed form` as effectively X−X (target restates the line-82 definition) — **non-blocking** because K_q's value is independently anchored upstream (stage-118 sympy line 50); `g_s closed form` genuine; `lambda from bilinear` modest (sign/factorization consistency). **117** (outlet_core_status, a status-consolidation card) F1 mathematica_transliteration = **accepted policy-mirror** (pure rational series-coefficient classification of which outlet-deformation classes preserve Ŷ₂; no independent-route rewrite forced on a consolidation card) + F2 tautological_check de-tautologized by importing the stage-116 forward tube-length closed form `L_W = pi a sqrt((1+r_c)/3)/2` to build `kappa0 = 4 L_W²/(pi² a²)` (the target-inverting dead `Lw_required` was DELETED) + F3 tautological_check = comment-only provenance fix (γ₀=(1+r_c)/9 corrected from false "Stage 119"/"derived upstream" to a postulated pure-scale ANSATZ; κ₀'s "derived upstream → stage 116" retained; `gamma_c = 1/9` assertion byte-unchanged) + F4 tautological_check reconciles already-wired capstone classification flags (one nontrivial survivor `compensated Robin-mixed core realization`) — **12 PASS per engine**. **CAVEAT (117 F2 honest strength):** F2 is now a falsifiable consistency check against the in-script forward L_W law (a wrong coefficient breaks the O(z⁶) core-collapse residual), NOT an in-stage re-derivation of the D/N BVP eigenvalue (that lives at stage 116). **122** (mouth_source_compensation_test, **SymPy-ONLY, no `.wl`**) F1 tautological_check de-tautologized — bare reciprocal traction ratios (`T_ratio = 1/g±`, never-fail `g·(1/g)−1`) replaced by ratios DERIVED from the stage-119 proportionality law `g = C/T_m` (C a branch-independent constant that cancels symbolically) + the `g_nat = 1` natural-branch ansatz, compared against the notes §5 closed form `1/g±` — **7 `= 0` lines**, confirms the natural-defect uses `100π²` (NOT `168π²`). **CAVEAT (122 moderate strength):** a consistency check of `g_nat=1` + C-branch-independence; would fail if `g_nat≠1` or the normalization were branch-dependent. The 117/118 resolution directions came from a Claude+Codex consult (Codex base `019e74f7`, all CONCUR), recorded at `redteam/codex_reviews/_consult_batch3.md`; paper-card cross-reference / informational-only items logged in `PAPER_CLEANUP_TRACKER.md` P4-44 (NO new paper-cleanup action items). See `redteam/verifications/stage_{117,118,119,122}.md`. Previous batch-2 snapshot below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 2** — re-verification of four more stages whose first-pass fixes had been applied orchestrator-direct while Codex was bypassed. Stages **105, 106, 109, 112** all REMEDIATION-`verified` (Codex applied + clean Claude verify agents, 2026-05-29); both engines exit 0/0 on every stage; this status **supersedes the earlier tainted orchestrator-direct "verified"** these four carried from their first-pass IV.2 batch. No new stages added to the verified count (all four were already counted from IV.2); this is a verification-hardening / integrity re-pass. **`material_change: false` on every stage** (no derived constant moved). Findings closed per stage (from verification frontmatter): **105** (checkpoint) F1 mathematica_transliteration resolved — `.wl` chi_Q derivation rewritten along an INDEPENDENT residue/`Reduce`-witness path (single-ratio retarded module, `Im[]` projection, operator-product polynomial inversion; non-transliterated names) + F2 paper_misalignment (stale "Stage 88" labels → canonical "Stage 105"), 9 Mathematica PASS; **106** F1 paper_misalignment (doc-only citation fix: item (ii)→Stage 102, item (iii) DtN fingerprint→Stage 104, chi_Q=1→Stage 105) + F2 tautological_check resolved (by-construction K0/K2/K4 relation → INDEPENDENT hardcoded `*_target` literal consistency `K0_target·K4_target = 4·K2_target²`) + F3 mathematica_transliteration resolved (`.wl` re-authored on independent `Yret` one-pole omega-series path deriving chi_Q=1 from ω⁵ match) + F4 insufficient_verification (Delta_Q first-order sensitivity slope = −Gamma5_target on both engines), 11 Mathematica PASS; **109** F1/F2 tautological_check resolved ("preservation substitution" de-tautologized on BOTH engines: substitutes the INDEPENDENT closed form `a_5 = -5b/9 - a0/27`; `.wl` also re-derived via separate num/den linearization + denominator series-inversion, `coeff=(chiSeries-1)/eps` intermediate removed) + F3 paper_misalignment (NO script change), 4 Mathematica PASS; **112** F1 paper_misalignment (stale "Stage 95" → "Stage 112") + F2 mathematica_transliteration resolved (independent Stage-92 linearized (b,a0,a5)=(0,3σ,−σγ) cross-check deriving gamma_W=1/9, a route the `.py` lacks) + F3 symbol_assumption_error resolved (`sigma_W != 0` qualifier folded into BOTH engines: preservation `chi_B=1 ⟺ sigma_W(1−9 gamma_W)=0`, nontrivial branch sigma≠0 forces gamma_W=1/9 via `Reduce[...&&sigma!=0]`, degenerate sigma=0 gives chi_B=1 for any gamma), 14 Mathematica PASS. **Paper_misalignment resolutions** (Claude+Codex consult `019e748e`, recorded at `redteam/codex_reviews/_consult_batch2.md`) all non-conceptual: 105 F2 & 112 F1 = canonical-internal-stage-number labels (per IV.4/IV.5 convention); 109 F3 & 106 F1 = paper-card cross-references (no script change), logged in `PAPER_CLEANUP_TRACKER.md` P4-43. See `redteam/verifications/stage_{105,106,109,112}.md`. Previous batch-1 snapshot below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 1** — re-verification of the four stages whose first-pass fixes had been applied orchestrator-direct while Codex was bypassed. Stages **108, 116, 151, 170** all `verified` (116 & 151 in prior commit `e1cdfec`; 170 & 108 this session); both engines exit 0/0 on every stage. Findings closed per stage (from verification frontmatter): **116** (batch IV.3) verdict `verified`, sympy/math 0/0, F1 mathematica_transliteration + F2 insufficient_verification resolved (independent eigenvalue solve of `Cos[u]==0` on `(0,π)`; tautological renorm block → labeled prints); **151** (batch IV.6) verdict `verified`, sympy/math 0/0, F1 insufficient_verification resolved — **methodology: Mathematica = full all-`Pi_star` symbolic authority; SymPy = an EXACT 5-point cross-check at rational `Pi_star ∈ {1/2,1,3/2,2,5/3}` symbolic in `r1,r2,A_T,B_T,gprime`** (supersedes the IV.6-snapshot mpmath description; `material_change: false`); **170** (batch V.1) verdict `verified`, sympy/math 0/0, F1 tautological_check resolved (Section-5 lanes rerouted through the derived Section-2 maps `dkappa_from_du2`/`dgamma_from_dP0`); **108** (batch IV.2) verdict `verified`, sympy/math 0/0, F2 mathematica_transliteration + F3/F4 tautological_check resolved (independent non-`Series` `chiGenAlt` route; Class A anchored to fingerprint `1+z²/9+4z⁴/81+iz⁵/27`; Class C anchored to `-Σ0/27`). **108 F1 = paper-card cleanup, no script change** (direction (a), Claude+Codex concur — see `redteam/codex_reviews/_consult_108_f1.md`; deferred to `PAPER_CLEANUP_TRACKER.md` P4-42). 108 carries `material_change: true` (verification surface strengthened; no derived value changed; no downstream propagation); 116/151/170 `material_change: false`. Deferred paper/doc cleanups (106/139/148 + 112 precision caveat) logged in `PAPER_CLEANUP_TRACKER.md` P4-42. No new stages added to the verified count (all four were already counted from their first-pass batches); this is a verification-hardening / integrity re-pass. Previous V.1 snapshot below.)

Snapshot date prior: `2026-05-28` (batch V.1 close — first-pass paper-grounded audit under v2 prompt for stages 164-175, no checkpoints, no status-only units. 22 findings closed across all 12 (dirty) stages, 0 blocked. 1 material_change (stage 170 — additive new Sec. 5 weak-axisymmetric (1,1/2,-1) signature coverage; no carried/derived result changed, zero downstream propagation). **Cluster A only**: all 12 stages carried a stale `STAGE N-17` script banner (164->147 ... 175->158); 11 fixed inside their per-finding directives, 4 residuals (164.py, 165.py+.wl, 174.py) mass-fixed in place. **No Cluster B** (no body-text citation re-attribution) and **no Cluster C** (no paper-card Checks downgrade). One **non-cluster paper_misalignment** (170 F1, weak-axisymmetric (1,1/2,-1) signature unverified) resolved by user direction (a) = add the missing check to both engines — first non-cluster paper_misalignment in the v2 run. `mathematica_transliteration` on 9/12 (highest v2 share); 7 independent routes added (164/166/170/171/172/173/174), 2 accepted as policy mirrors (169 F2, 175 F3-step3). Three orchestrator catches: (1) 166 matrix round-trip vector residual -> `Total[(...)^2]` scalarization; (2) 175 F1 cross-check reduced to a simplify-commutes identity -> kept only the load-bearing `Sigma_N-dln(Lambda^2/K)`; (3) 171 F1 directive route was tautological (`zCombFormula` rebuilt from the same `D[...]` as `zCombExact`) -> reworked to collected-literal target + independent `Series` route, caught by the first verification (`needs_rework`) and re-verified clean. Thirteenth consecutive batch clear of stop-cold.)

Previous snapshot: `2026-05-27` (batch IV.5 close — first-pass paper-grounded audit under v2 prompt for stages 139-150, no checkpoints, zero material_change. Same renumbering pattern as IV.4 — 11 `.wl` banners + 6 `.py` banners off by -17 (Cluster A) and 4 notes H1 lines off by +102; plus 22 body-text forward-stage citations off by -51 or -102 re-attributed to current upstream (Cluster B). Stage 144 paper card Checks downgraded to carry-forward citations of 135/140 via Cluster C, mirroring IV.4's stage 134. Three status-only stages clean (141, 145, 149). Eleventh consecutive zero-redirection batch. Four orchestrator catches in the rework loop including: stage 148 directive's prescribed closed form was wrong (`168π²` typo from stage 148 notes; correct `100π²` confirmed at stage 126 upstream), stage 139 pitfall #13 recurrence on the `Pi_*)` comment substring, plus three SymPy tolerance/precision/Integrate-fallback adjustments.)

## Scope

- Canonical stage range: `001--253`
- Canonical stage source model:
  - Parts `I--VII`: compact stage cards
  - Part `VIII`: full derivation-stage files

## Coverage Totals

| Metric | Count |
|---|---:|
| Total stages in archive | 253 |
| SymPy audits present | 240 |
| Mathematica audits present | 242 |
| Numerical stress artifacts present | 15 |
| Stage-specific review notes present | 178 |
| Stages with no SymPy audit | 13 |
| Stages with no executable audit | 11 |
| Stages with Mathematica but no SymPy | 2 |

Mathematica counts above are presence counts only. They are **not**
independence counts. See `MATHEMATICA_MIRROR_POLICY.md` for the current rule
that distinguishes secondary replay coverage from genuinely independent
Mathematica mirrors.

The SymPy runner summary reports `241` passing files because it includes the
repo-level `moving_throat_pde_master_sympy_audit.py` in addition to the `240`
stage-level SymPy audits counted above.

### Red-team verification status

Distinct from "audit file present" (counted above) is "red-team verified" —
the stage has passed the `audit → directive → codex → verifier` pipeline run
out of `redteam/`, with both engines independently checking load-bearing
claims and a clean-context verifier agent confirming the directive's intent
was honored. See `redteam/BATCHES.md` for the live batch table.

As of `2026-06-03`: **242** of 253 stages red-team verified (95.7%) — the authoritative MANIFEST count (230 verified in range 001–230 + 12 in range 231–242; stages 243–253 all `pending`). With VII.2 closed (first audit pass for stages 231–242 under v2, the **third batch of Part VII**, theme "Rigid-mouth orbit-lock / branch-dressing / twin-support"), the entire range 001–242 is now paper-aligned at v2 depth. **239 and 242 are the two checkpoints** in VII.2 (BOTH cleared the higher bar — each EXISTING transliteration `.wl` RE-AUTHORED to an independent route: 239 forward-Jacobian + native `PseudoInverse` left-inverse + `Reduce`/`Equivalent` orbit-lock, 242 `Resolve[ForAll,Reals]` strict-inequality certificate + `D[]` on real closed forms + `logDrift`, the twin-window strict inclusion `C_mix < Pi_tr < 2 C_mix` now tested strictly on both engines). All 12 verified; **`material_change: false` on all 12**; 0 stop-cold, 0 ultimately-blocked; 11 of 12 codex-invoke exited 0 on iteration 1, one iter-2 (stage 238, the support-blindness F1/F4 reframe). **Dual-engine: 0 sanctioned mirrors** — 10 NEW independent-route `.wl` (231–238, 240, 241) + 239/242's `.wl` re-authored = 12 independent Mathematica audits. **3 notes-only paper_misalignment numerical/coefficient typos** (231 `240·δ²ξ→189·δ²ξ`/`189·ξ³→121·ξ³`; 232 figure-of-merit prefactor `168→100`; 241 `ϱ_WΛ` upper bound `193/369→125/369` — the published cards/appendices were UNAFFECTED, carrying abstract forms) corrected to the verified SymPy scripts and each cross-engine corroborated by the new `.wl`. Dominant defect theme: the "variable-independence self-test trap" (237-F2, 238-F1, 240-F1). Non-checkpoint constants now independently corroborated by the new `.wl`: 240 `ρ_α=4/3, ζ_req=1/3, Pi_tr=(4/3)C_mix`; 241 `ϱ` windows `1/3, 125/369, 2/3, 250/441`; 232 `Ξ` prefactor=100. RESIDUAL notes-title renumber drift (232/234/235/236) DEFERRED to the post-253 stem-keyed pass. See the VII.2 snapshot above, `MATHEMATICA_MIRROR_POLICY.md`, and `PAPER_CLEANUP_TRACKER.md` P4-54.

As of `2026-06-02`: **230** of 253 stages red-team verified (90.9%) — the authoritative MANIFEST count (218 verified in range 001–218 + 12 in range 219–230; stages 231–253 all `pending`). With VII.1 closed (first audit pass for stages 219–230 under v2, the **second batch of Part VII**, theme "Mixed-bundle / resonance / branch-packet"), the entire range 001–230 is now paper-aligned at v2 depth. **221 is the only checkpoint** in VII.1 (cleared the higher bar — its existing transliteration `.wl` RE-AUTHORED to an independent route, F1 survival round-trips de-tautologized, deliverable #9 linear-survival-window now genuinely covered in BOTH engines). All 12 verified; **`material_change: false` on all 12**; 0 stop-cold, 0 blocked; all 12 codex-invoke exited 0 on iteration 1 (no iter-2 reworks). **Dual-engine: 0 sanctioned mirrors** — 11 NEW independent-route `.wl` (219, 220, 222–230) + 221's `.wl` re-authored = 12 independent Mathematica audits. **7 notes-only paper_misalignment numerical typos** (222/223/227/228/229 — a systematic +68/+51 additive family; the published cards/appendices were UNAFFECTED, carrying abstract forms) corrected to the verified SymPy scripts and each cross-engine corroborated by the new `.wl`, plus notes-label renumbers (221, 225). The VII.1 load-bearing constants now dual-engine-anchored: the corrected R_Q figures (145.483858657863; 138.814136942081/137.502546600713), the i=h rigidity determinant factor (200+147π²), the δ_1 coefficient (196π²/(98π²−25)), the crossover-cubic leading coeff (121ξ³), and the 230 thresholds R_*≈1.229255438463336 / δ_*≈0.723111617875019. 230's codex log was captured EMPTY (a logging anomaly, NOT a stall — exit 0, all `## Applied` blocks present, independently re-confirmed). See the VII.1 snapshot above, `MATHEMATICA_MIRROR_POLICY.md`, and `PAPER_CLEANUP_TRACKER.md` P4-53.

As of `2026-06-02`: **218** of 253 stages red-team verified (86.2%) — the authoritative MANIFEST count (200 verified in range 001–200 + 18 in range 201–218; stages 219–253 all `pending`). With VI.1 closed (first audit pass for stages 201–218 under v2, the **first batch of Part VI**, theme "Explicit realization, scalar slice, ray ranking"), the entire range 001–218 is now paper-aligned at v2 depth. **203 and 218 are the two checkpoints** in VI.1 (BOTH cleared the higher checkpoint bar, no rubber-stamp). All 18 verified; **`material_change: false` on all 18**; 0 stop-cold, 0 blocked, 0 needs_rework left open. **Dual-engine: 0 sanctioned mirrors** — 16 stages got a NEW independent-route `.wl` (201, 202, 204–217), 218's pre-existing transliteration `.wl` was RE-AUTHORED to an independent route (M4 witness counts differ across engines), and 203's `.wl` got a strengthened independent composition check. The load-bearing corrected constant is the per-envelope lifted Bézout bound `162 = 3⁴·2` (217/218), arithmetically forced — the published card/appendix/notes carried wrong typos (179, 230) corrected to 162 (script always correct, no script value moved). One iter-2 timeout rework (202 `.wl` → `LinearSolve` of the log-linearized system after a 600s-cap timeout) and one iter-2 forward-ref follow-up (212). See the VI.1 snapshot above, `MATHEMATICA_MIRROR_POLICY.md`, and `PAPER_CLEANUP_TRACKER.md` P4-52.

As of `2026-06-01`: **200** of 253 stages red-team verified (79.1%) — the authoritative MANIFEST count (187 verified in range 001–187 + 13 in range 188–200; stages 201–253 all `pending`). With V.3 closed (first audit pass for stages 188-200 under v2, the **second batch of the resumed first pass**), the entire range 001-200 is now paper-aligned at v2 depth. **200 is the only checkpoint** in V.3 (verified at the higher bar; F1 transliteration de-transliterated, F2 HIGH tautological de-tautologized). 7 clean (190, 192, 194, 196, 197, 198, 199), 6 dirty (188, 189, 191, 193, 195, 200); 13 findings closed, 0 blocked, 0 stop-cold, **`material_change: false` on all 13**. **The defining event was the dual-engine rule correction**: `.wl` now REQUIRED wherever Mathematica CAN independently verify (test = "is it possible"). Stages 188–199 were previously SymPy-ONLY; **12 new independent-route `.wl` were created (all genuine, ZERO sanctioned mirrors)** and 200's `.wl` was de-transliterated. One iteration-2 rework (189 Section II — back-definition demoted to a printed definition + independent direct-slope bridge). Two stale-banner relabels (189, 191; string-only). See the V.3 snapshot above and `MATHEMATICA_MIRROR_POLICY.md`.

As of `2026-05-30`: **199** of 253 stages red-team verified (78.7%). With V.2 closed (first audit pass for stages 176-187 under v2, the **first batch of the resumed first pass** after the IV.x/V.1 integrity remediation closed, run under the restored Codex-as-fix-applier contract), the entire range 001-187 is now paper-aligned at v2 depth. **185 is the only checkpoint** in V.2 (verified at the higher bar; F1 required an iteration-2 rework making the coefficient load-bearing via an independent slippage-law reconstruction). 2 clean (179, 180), 10 dirty; ~24 findings closed, 0 blocked, 0 stop-cold, **0 paper_misalignment**, **`material_change: false` on all 12**. `mathematica_transliteration` hit 7/12 — **all 7 given genuine independent routes, ZERO accepted as policy mirrors** (the converse of V.1's 2). One other iteration-2 rework (181 F1, vacuous ζ-free Mathematica check). No Cluster B or Cluster C this batch; the stale `STAGE N−17` banner cluster was fixed script-side (clean stages 179/180 carry residual cosmetic banner drift, non-blocking).

As of `2026-05-28`: **187** of 253 stages red-team verified (73.9%). With V.1 closed (first audit pass for stages 164-175 under v2), the entire range 001-175 is now paper-aligned at v2 depth. No checkpoints, no status-only units in V.1. **One material_change** (stage 170 — additive Sec. 5 weak-axisymmetric (1,1/2,-1) signature coverage added per user direction (a); no carried/derived result changed, zero downstream propagation). **First non-cluster paper_misalignment in the v2 run** (170 F1), resolved by adding the missing script check rather than a carry-forward citation — content-disambiguation confirmed the deliverable is stage 170's own (notes Sec. 5), with 171/173 verifying the lane pattern only on other quantities. `mathematica_transliteration` hit 9/12 (highest v2 share); 7 independent routes added, 2 policy mirrors (169, 175). No Cluster B or Cluster C this batch.

As of `2026-05-28`: **175** of 253 stages red-team verified (69.2%). With IV.6 closed (first audit pass for these 13 stages under v2), the entire range 001–163 is now paper-aligned at v2 depth. No checkpoints in IV.6; three stages clean (153 status-only, 159 and 162 substantive-clean). **One material_change** in IV.6 (stage 151 SymPy rewrite from symbolic-integration to mpmath numerical because `sp.integrate` of `Pi_star * exp(-Pi_star*x) * cos(pi*x/2)` over [0,1] with free `Pi_star` hung >30 min CPU); downstream impact zero since the script verifies algebraic identities not numeric carry-forwards. **First forward (downstream) carry-forward in v2 batches** via Cluster C on stage 158 (items 2-3 cite `\ref{stage:159}` / `\ref{stage:162}` / `\ref{stage:163}`). The −51 (IV.4 refs) and −102 (IV.5 refs) body-citation offsets carried over from IV.4/IV.5; new in IV.6 is the **−85 offset** for IV.6 self-references (239-248 → 154-163), disambiguated from −102 by content cross-check (e.g., stage 158's "Stage 241, exact co-evolving compensated Family-1 point" matches stage 156 renormalized canonical branch, not stage 139 family1_actual_mouth_gains).

As of `2026-05-27`: **162** of 253 stages red-team verified. With IV.5 closed (first audit pass for these 12 stages under v2), the entire range 001–150 is now paper-aligned at v2 depth. No checkpoints in IV.5; three stages clean (141, 145, 149 — all status-only). **Zero material_change** in IV.5 — every fix was renumbering, notes re-attribution, paper-card downgrade, or script-substance addition; no derived numerical constants moved. Pitfall #13 (Mathematica `_*)` substring in comment closing prematurely) **re-confirmed** at stage 139 — workaround applied (ASCII-safe names `piStar`, `sQStar`). The two systematic stage-renumbering shifts and the body-text forward-citation pattern that surfaced in IV.4 repeated in IV.5 with the same offsets (banners −17, notes H1 +102, body citations −51 or −102).

| Batch | Range | Stages | Verified | Date |
|---|---|---:|---:|---|
| I.1 | `001--012` | 12 | 12 | 2026-05-21 (v1) / 2026-05-25 (v2 paper-grounded) |
| I.2 | `013--023` | 11 | 11 | 2026-05-21 (v1) / 2026-05-25 (v2 paper-grounded) |
| II.1 | `024--036` | 13 | 13 | 2026-05-22 (v1) / 2026-05-26 (v2 paper-grounded) |
| III.1 | `037--048` | 12 | 12 | 2026-05-22 (v1) / 2026-05-26 (v2 paper-grounded) |
| III.2 | `049--060` | 12 | 12 | 2026-05-22 (v1) / 2026-05-26 (v2 paper-grounded) |
| III.3 | `061--072` | 12 | 12 | 2026-05-22 (v1) / 2026-05-26 (v2 paper-grounded) |
| III.4 | `073--084` | 12 | 12 | 2026-05-25 (v1) / 2026-05-27 (v2 paper-grounded) |
| III.5 | `085--090` | 6 | 6 | 2026-05-27 (v2 paper-grounded, first-pass) |
| IV.1 | `091--102` | 12 | 12 | 2026-05-27 (v2 paper-grounded — first pass) |
| IV.2 | `103--114` | 12 | 12 | 2026-05-27 (v2 paper-grounded — first pass) |
| IV.3 | `115--126` | 12 | 12 | 2026-05-27 (v2 paper-grounded — first pass) |
| IV.4 | `127--138` | 12 | 12 | 2026-05-27 (v2 paper-grounded — first pass) |
| IV.5 | `139--150` | 12 | 12 | 2026-05-27 (v2 paper-grounded — first pass) |
| IV.6 | `151--163` | 13 | 13 | 2026-05-28 (v2 paper-grounded — first pass) |
| V.1 | `164--175` | 12 | 12 | 2026-05-28 (v2 paper-grounded — first pass) |
| V.2 | `176--187` | 12 | 12 | 2026-05-30 (v2 paper-grounded — first pass, resumed; checkpoint 185) |
| V.3 | `188--200` | 13 | 13 | 2026-06-01 (v2 paper-grounded — first pass, resumed; checkpoint 200; dual-engine rule correction) |
| VI.1 | `201--218` | 18 | 18 | 2026-06-02 (v2 paper-grounded — first pass; checkpoints 203, 218; explicit realization, scalar slice, ray ranking) |
| VII.1 | `219--230` | 12 | 12 | 2026-06-02 (v2 paper-grounded — first pass; checkpoint 221; mixed-bundle / resonance / branch-packet) |
| VII.2 | `231--242` | 12 | 12 | 2026-06-03 (v2 paper-grounded — first pass; checkpoints 239, 242; rigid-mouth orbit-lock / branch-dressing / twin-support) |
| VII.3 onward | `243--253` | 11 | 0 | pending |

Cumulative findings closed: ~451 (~219 v1 + 10 v2 from I.1 + 10 v2 from I.2 + 18 v2 from II.1 + 13 v2 from III.1 + 16 v2 from III.2 + 13 v2 from III.3 + 14 v2 from III.4 + 15 from III.5 + 27 from IV.1 + 16 from IV.2 + 27 from IV.3 + 22 from IV.4 (21 resolved + 1 blocked_legitimate on 137 F3 matrix-Schur) + 31 from IV.5 (30 resolved + 1 blocked_legitimate on 144 F4 transliteration policy) + **19 from IV.6 (19 resolved + 0 blocked, 1 material_change on 151 mpmath rewrite)** + **22 from V.1 (22 resolved + 0 blocked, 1 material_change on 170 — additive Sec.5 (1,1/2,-1) signature coverage)**, plus 1 blocked-legitimate from IV.1). Of the 27 IV.3 closes, 7 were `paper_misalignment` (Cluster A 4 notes-side numerical typos at 121/122/123/126 — `168π²→100π²` for three stages and `228→160` for stage 123; Cluster B λ sign flip at 118 — internal script inconsistency; Cluster C integral inequality coverage gap at 125; **Cluster D 117 transliteration + 2 tautological resolved via cite-upstream-and-downgrade**), 3 `mathematica_transliteration` (115, 117 F1, 125 F2 — independent-derivation insertions), 6 `tautological_check` (116 F1, 116 F2, 117 F2, 117 F3, 117 F4, 119 F1), 5 `insufficient_verification` (116 F4, 119 F2, 121 F2, 122 F2, 126 F3), 2 `hardcoded_result` (116 F3, included in F3 count), 1 `script_missing_paper_claim` (121 F3 Ω_W), 1 `paper_missing_script_claim` (126 F2 positivity), 1 `stale_output` banner (126 F4). Plus a 10-site banner-relabel sweep across all IV.3 scripted stages (115, 116, 117, 121, 122, 123, 125, 126 — orchestrator-direct, matching IV.2's pattern).
`tautological_check` dominant overall, `mathematica_transliteration` second.
`hardcoded_result` rose sharply in III.4 to 12 because the Family-1 numerology
cluster 075-084 packs many literal constants; III.5 quieted again (1 hardcoded
in 089 F4). v2 added `paper_misalignment` as the 10th category — **31** items
total across the eight v2 batches (7 in I.1, 3 in I.2, 3 in II.1, 3 in III.1,
2 in III.2, 4 in III.3 — 2 substantive + 2 banner relabels, 7 in III.4 — 4
substantive + 3 audit-flagged banner relabels, plus an 8-stage
orchestrator-direct banner-relabel sweep when the global-renumber leftover
turned out to be pervasive across III.4; **2 in III.5 — both substantive (087 F1
status/checkpoint consolidation, 089 F1 Pe_req=0 chain closure), plus a 12-script
orchestrator-direct banner-relabel sweep**);
zero user redirections in II.1, III.1, III.2, III.3, III.4, III.5, IV.1, IV.2, IV.3, IV.4, IV.5, IV.6, **V.1** (13 consecutive
batches — Codex was bypassed from III.5 onward per the III.4 availability lesson; orchestrator-direct
math-authority worked cleanly because the audit + grep evidence was conclusive. V.1 surfaced the first non-cluster
paper_misalignment (170 F1); the user took the recommended direction (a) — add the missing check — so no conclusion
was redirected and the streak holds).
v2 surfaces `insufficient_verification` prominently — 8 in II.1, 5 in III.1, 8 in III.2,
4 in III.3, 1 in III.4, 1 in III.5, 7 in IV.1, 2 in IV.2, 5 in IV.3, 4 in IV.4 (130 F1, 130 F2, 131 F3, 137 F1), 11 in IV.5 (139 F1, 140 F2, 142 F2, 143 F1, 144 F2, 144 F3, 146 F2, 147 F3, 148 F1) = **56** cumulative.
Stage 060 (v1 `material_change: true`) returned **clean (0 findings)** under III.2 v2.
**Stage 068 (v1 `material_change: true`) returned clean at v2**.
Stages now carrying `material_change: true`: 001, 004
(I.1 v2); 013, 014, 015, 018 (I.2 v2); 045 (III.1 v2 — structural-only, F_tr
export value unchanged); 060 (v1, clean at v2); 068 (v1, clean at v2). II.1,
III.1, III.2, III.3, III.4, III.5 v2 each added **zero** value-changing material_change. **IV.1 added one structural material_change at stage 100** (closure derivation strengthened from tautological cross-check to substantive `mhat_0^2 Gamma_5 = Gamma_5_target` imposition; no derived value changed; downstream stages > 100 not marked `upstream_stale`). **IV.2 added one structural material_change at stage 108** (Cluster B β-parameterized preservation submanifold added; the β=1 reduction already verified previously is unchanged, only the verification surface widened; downstream stages > 108 not marked `upstream_stale`). **IV.3 added one substantive material_change at stage 118** (λ sign flip from + to − in section V's `lam_uniform = qstar*v0*I_sq` closure; internal script inconsistency where section IV's bilinear derivation and the notes' boxed form both had `−` but section V dropped the minus during integration. Resolved by aligning section V with section IV. Downstream Schur reductions use `K_s K_q + λ²` and `(K_s g_q − λ g_s)²` — both sign-invariant under squaring — so no numerical downstream propagation; upstream_stale NOT flagged. Future-batch caveat: revisit if a downstream stage uses λ in an unsquared cross-term).
III.3 v2 introduced one orchestrator hot-fix on stage 064 Mathematica
(`Integrate[]` with symbolic functions does not factor constants — verify
integrands first; pitfall #9 candidate). III.4 v2 introduced one orchestrator
fix on stage 082 SymPy (`sp.nsolve` is unstable for `y tan y = 37` near
`pi/2` and jumps to far-away roots — use `mpmath.findroot(..., solver="bisect")`
instead; pitfall #10 candidate). **III.5 introduced two orchestrator hot-fixes:
(a) on stage 088 SymPy, `Y_rho.subs(omega**2/Omega_Q**2, u)` failed silently
because `sp.simplify` reshapes the denominator into `(Omega_Q**2 - omega**2)`
form and the combined ratio is no longer a syntactic subexpression — fix:
substitute `omega**2 -> u * Omega_Q**2` then `sp.simplify`. (b) On stage 088
Mathematica, a comment containing the substring `stage085_*)` was prematurely
closed by the embedded `*)`, causing `Syntax::sntx` and silently skipping the
F1 assertion and regime trichotomy while still reaching `Exit[0]` (rc=0
masking a partial run) — fix: reword to avoid `*)` substrings in comment
text. New pitfall #11 candidate.**
Pitfall #8 was promoted from candidate to
documented in `codex.md` "Common cross-engine pitfalls" item #1 before
III.3 launched. Pitfalls #6, #7 remain candidates; #9 (Mathematica
`Integrate[]` constant factoring), #10 (SymPy `nsolve` near
singularities), and #11 (Mathematica `*)` substring inside comment body
closes prematurely; verifier must check that all expected PASS lines appear,
not just `rc=0`) added in III.5. **IV.1 added no new pitfall candidates** —
all orchestrator-direct edits applied first-attempt clean (Cluster A docstring
carry-forwards + Cluster B 100 closure derivation + Cluster C 23-site banner
sweep), with verifier PASS-line counting confirming all expected substantive
checks across the 12 stages. **IV.2 added no new pitfall candidates** — all
orchestrator-direct edits applied first-attempt clean (Cluster A 106/109
docstring carry-forwards + Cluster B 108 β-locus extension + Cluster C 24-site
banner sweep), plus a Mathematica parse-bug correction (`chiArg /. beta -> 1 - 1`
parsed as `beta -> 0` at 108 F2; fix: `(chiArg /. beta -> 1) - 1`) which
re-confirms pitfall #11 PASS-line discipline (the buggy line passed by accident).
**IV.3 added pitfall #12 candidate** (Mathematica `Solve[expr == 0, frakG]`
fails with "frakG is not a valid variable" when `frakG` is bound to its definition
`= gQ*Sqrt[kS]/(gS*Sqrt[kQ])` — always introduce a fresh symbol for Solve's
target variable, then substitute back). Also re-confirmed: Mathematica `Minimize[]`
on `cos` over an interval frequently returns the input unevaluated; use boundary-value
checks instead under monotone-decreasing assumptions. Directive-correction during
edit: (a) stage 116 directive's `gamma0_from_D = -I·coeff(z,5)` produced wrong sign;
orchestrator corrected to `+I·coeff(z,5)` (analytic check confirms `+I` is right);
(b) stage 115 directive's `parentFamilyResidual - balanceEq*(kS*kQ)/(gS^2)`
multiplicative factor was wrong because `balanceEq` carries denominator
`kS*(kQ*kS + lam^2)`; correct factor is `(kS*kQ + lam^2)/(gS²*kQ)`. Both
corrections caught by re-running the scripts; sub-agents would otherwise have
silently produced FAIL residuals. See per-batch summaries in
`redteam/batches/batch_<ID>.md`.

### Linear projected-EM update

Snapshot addendum: `2026-05-11`

Stages `004--021` are now canonical linear stages for the projection-first
Maxwell integration, parent-action packet, and retained reduced one-port normal
form.  They are no longer counted as Stage `004` substages.  The old compact
Stage `004` reduced Maxwell/mixed calculation is retained as Stage `021`.

Stages `004--020` have file-for-file SymPy migrations from the derivation-only
`notes/em_projected` scripts through
`step_18_parent_throat_action_weak_axisym_packet`.  The `step_19_*`
branch-export packet and `step_20+` computational/runtime diagnostics remain
excluded from paper coverage.  File-for-file Mathematica mirrors for
Stages `004--020` have not been independently derived yet; Stage `021` retains
the legacy reduced Mathematica audit.

## Coverage By Part

| Part | Stage Range | Total | SymPy | Mathematica | Numerical | Review |
|---|---|---:|---:|---:|---:|---:|
| I | `001--023` | 23 | 23 | 6 | 2 | 6 |
| II | `024--036` | 13 | 13 | 13 | 0 | 13 |
| III | `037--090` | 54 | 53 | 54 | 2 | 54 |
| IV | `091--163` | 73 | 61 | 59 | 8 | 69 |
| V | `164--200` | 37 | 37 | 25 | 1 | 25 |
| VI | `201--218` | 18 | 18 | 2 | 0 | 2 |
| VII | `219--242` | 24 | 24 | 3 | 0 | 3 |
| VIII | `243--253` | 11 | 11 | 3 | 2 | 3 |

## Coverage Classes

| Coverage class | Count | Stage ranges |
|---|---:|---|
| SymPy + Mathematica | 163 | See `STAGE_PROVENANCE_INDEX.md` for the file-by-file list. |
| SymPy only | 77 | Includes projected-EM Stages `004--020` and later stages whose Mathematica mirrors are not yet present. |
| Mathematica only | 2 | `084`, `093` |
| No executable audit | 11 | `103`, `113`, `120`, `124`, `128`, `132`, `136`, `141`, `145`, `149`, `153` |

## Constant Provenance Rule

Coverage counts are not trust grades.

Likewise, `Mathematica present` is not the same thing as `independent second
CAS derivation`. Repo-wide counts track coverage breadth; independence is now a
separate policy classification in `MATHEMATICA_MIRROR_POLICY.md`.

For this archive, an audit should be treated as insufficient until every
nontrivial constant used in it is classified as one of:

- `derived in audit`
- `carried forward with source anchor`
- `probe-only numeric value labeled`

Any unexplained literal should be treated as a verification defect, not a style
issue.

## Immediate Gaps

### 1. No executable audit

These stages currently have neither SymPy nor Mathematica nor numerical stress:

`103`, `113`, `120`, `124`, `128`, `132`, `136`, `141`, `145`, `149`, `153`

These are the first places where the archive has mathematical content without an
executable backstop.

### 2. Mathematica without SymPy

These stages have a Mathematica artifact but no SymPy mirror:

`084`, `093`

These need reconciliation before we can make strong claims about dual-CAS
coverage.

### 3. SymPy-only frontier

The current SymPy-only region is:

`121--123`, `188--199`, `201--202`, `204--217`, `219--220`, `222--238`, `240--241`, `244--247`, `249--252`

Operationally, the main contiguous Mathematica mirror still ends at Stage `187`,
with isolated later checkpoints now hardened at Stages `200`, `203`, `218`,
`221`, `239`, `242`, `243`, `248`, and `253`.

### 4. Review-coverage gap

Stage-specific review notes are missing for:

`121--124`, `188--199`, `201--202`, `204--217`, `219--220`, `222--238`, `240--241`, `244--247`, `249--252`

That means the late-stage archive is not only missing Mathematica coverage, but
also lacks the earlier review-note pattern that exists for most of Stages
`001--187`.

### 5. Numerical stress remains sparse

Numerical stress coverage exists only for the following stage families:

`003/021`, `045--048`, `058/170`, `106`, `125`, `142--144`,
`146--148`, `150--152`, `154`, `155--156`, `157`, `185--187`,
`248`, `253`

This is a narrow slice relative to the symbolic verification surface and should
be treated as targeted spot-checking, not broad numerical validation.

## How To Use This Baseline

1. Use `STAGE_PROVENANCE_INDEX.md` when you need the exact note or script path
   for a stage.
2. Use this baseline when deciding which verification gaps are structural and
   which are only missing a second engine or review note.
3. Use the part-level counts to prioritize the next audit wave.
4. Use the coverage classes to define the future-paper citation support set.
5. Use `CITATION_SUPPORT_SET.md` when deciding which gaps are most important for
   downstream citation hardening.
6. Use `CHECKPOINT_CONSTANT_PROVENANCE.md` for the growing no-magic-numbers log
   on the checkpoint subset.
7. Use `CHECKPOINT_TRUST_AUDIT.md` for the current stage-level trust baseline
   on that checkpoint subset.

## Recommended Next Verification Sequence

1. Use `CHECKPOINT_TRUST_AUDIT.md` as the current checkpoint trust baseline.
2. Reconcile the remaining repo-wide Mathematica-only outliers:
   `084`, `093`, `103`, `113` (last two are IV.2 status-only).
3. Use `MATHEMATICA_MIRROR_POLICY.md` when deciding whether a mirror gap is an
   execution-coverage gap or an independence gap.
4. Then widen the audit wave to the remaining repo-wide gaps.
5. Backfill executable audits for the remaining no-executable stages:
   `128`, `132`, `136`, `141`, `145`, `149`, `153`. (103, 113, 120, 124
   cleared as status-only consolidation cards under IV.2/IV.3.)
6. Use the new `248` / `253` harnesses as the template if numerical-stress
   coverage is widened beyond the current spot-check set.

## Pass 2 — Batch I.1 (2026-06-04)

Pass-2 Batch I.1 (001-012): all 12 re-verified at v2 depth + exhaustive
value-reconciliation augmentation; 12/12 verified. Coverage EXTENDED at stage
009: added a generic-kernel first-moment (μ₁) check in BOTH engines (Gamma
kernel w(u)=u·e^{-u}, μ₁=2; O(ℓ) leading term asserted = q0+ℓ·μ₁·q1) — the
card's general μ₁ formula is now exercised, not just the exponential special
case (μ₁=1). Stage 003 .wl `lRed` consolidated to one parenthesized assignment
(robustness; math no-op, output byte-identical).

## Pass 2 — Batch I.2 (2026-06-04)

Pass-2 Batch I.2 (013-023, `Part I.2 — Maxwell bridge, parent throat action,
reduced one-port`): all 11 re-verified at v2 depth + exhaustive
value-reconciliation augmentation; 11/11 verified, `material_change: false` on
all 11, 0 stop-cold, 0 blocked. Every emitted deliverable value reconciles —
86 values checked batch-wide (013=6, 014=5, 015=4, 016=9, 017=6, 018=2, 019=5,
020=11, 021=11, 022=14, 023=13), 0 MISMATCH, 0 MISSING-DELIVERABLE. All 11 are
genuine dual-engine (both independent engines from pass 1; re-audit confirmed
independence). 10 clean (013-020, 022, 023); one finding at 021 — a low-severity
`stale_output`: the committed Mathematica output `.txt` was a stale
pre-numbering-reconciliation capture (banner "STAGE 004") and the live `.wl:195`
carried a leftover bare `Print["Stage 004 Mathematica audit passed."]` label
that the numbering-recon scan (keyed on the `banner[...]` call + docstrings)
missed; fix = one label-only `.wl:195` edit `Stage 004`→`Stage 021` + output
refresh (math no-op, no assertion/result changed). No coverage count moved (all
11 were already verified+dual-engine from pass 1).

## Pass 2 — Batch II.1 (2026-06-04)

Pass-2 Batch II.1 (024-036, `Part II.1 — Overlap isotropy through continuum
kernel`, 13 stages): all 13 re-verified at v2 depth + exhaustive
value-reconciliation augmentation; 13/13 verified, `material_change: false` on
all 13, 0 stop-cold, 0 blocked. Every emitted deliverable value reconciles —
177 values checked batch-wide, 1 MISALIGNED (024's `4π/122` sixth-moment typo,
resolved below), 0 MISSING-DELIVERABLE. All 13 are genuine dual-engine (both
independent engines from pass 1; e.g. 024's `.wl` integrates sphere moments
directly vs SymPy's pairing-sum); no new `.wl`, 0 sanctioned mirrors. Checkpoints
024 and 036 are in range (both re-verified; no checkpoint constant changed —
κ_*=√5/(7√π) and M^(20) unchanged). **7 clean** (025-031); **6 with
findings→resolved→verified** (024, 032-036). **024 — TWO findings:** (F1)
`paper_misalignment` value_mismatch RESOLVED via Claude+Codex (non-conceptual,
NOTES-ONLY, published `.tex` card UNAFFECTED): notes:213 sixth-moment prefactor
`4π/122`→`4π/105` (textbook 1/(3·5·7)=1/105 in d=3; both engines use 105; the
notes' own κ_*=√5/(7√π) requires 105); Codex applied the single notes hunk, no
script change; (F2) `stale_output` — transcripts refreshed (Mathematica banner
stale `STAGE 007`→`STAGE 024`). **032/033/034/035/036 — each a single
`stale_output`** (NOT `paper_misalignment`): committed transcripts predated the
June-3 numbering commit, leaving residual stale self-labels (SymPy docstring +
closing `All Stage NN checks passed.` print, + a `.wl` Print on 032 + a `.wl`
sub-stage comment on 033); LABEL-ONLY fixes `Stage {15,16,17,18,19}`→`Stage
{32,33,34,35,36}` (matching each file's own canonical banner) + both-engine
re-run; same class as the I.2 stage-021 `.wl:195` fix; zero math change. No
coverage count moved (all 13 were already verified+dual-engine from pass 1).

**⚠️ DISCOVERY (deferred cleanup — see also `PAPER_CLEANUP_TRACKER.md`):** the
"mechanically exhausted (001-253)" numbering reconciliation did NOT clear the
II.1 script/output band. Stale labels REMAIN (left UNTOUCHED — content-dependent,
partly ambiguous, NEVER offset-sweep, reserved for a careful gated per-reference
pass): clean-stage committed OUTPUTS are themselves stale with unflagged stale
self-banners (028 out `STAGE 011`; 030 `STAGE 13 AUDIT COMPLETE`; 031 `STAGE 14
AUDIT COMPLETE`; 025 `Stage 8 … passed`; 026/027 `STAGE-8/9` cross-refs + self
pass-lines); ambiguous self-vs-cross multi-epoch refs (024 `.wl:293 FINAL
STAGE-007 LEDGER` + `.py "Stage-6"`; cross-refs 033 `Stage-15`=032, 036
`Stage-18`=035). The MATH is verified clean regardless. Recommendation: a
dedicated careful numbering pass over the script/output band, separate from the
red-team loop. Reference memory `numbering-drift-root-cause`.

## Pass 2 — Batch III.1 (2026-06-05)

Pass-2 Batch III.1 (037-048, `Part III.1 — Continuum kernel, generalized branch,
rank-2`, 12 stages): all 12 pass-2 re-verified at v2 depth + exhaustive
value-reconciliation augmentation; 12/12 verified, `material_change: false` on
all 12, 0 stop-cold, 0 blocked, 0 Codex deviations, all iter-1 exit 0. Every
emitted deliverable value reconciles — 143 values checked batch-wide, 0
MISALIGNED, 0 MISSING-DELIVERABLE. All 12 are genuine dual-engine (both
independent engines from pass 1; e.g. 039's `.wl` derives δ_split via
`a1Direct/a0Expected−1`, 043's `.wl` uses `Det`/`Limit`/`Series` + numeric sign
points); no new `.wl`, 0 sanctioned mirrors. NO checkpoints and NO EM-projected
stages in range. **2 real script-side math findings→fixed→verified** (039, 043);
**9 label-only self-label fixes→verified** (037,038,040,041,042,044,046,047,048);
**1 deferred-clean→verified** (045). **039 — TWO findings:** (F1)
`tautological_check` — the only `R_U` direction-factor check was identically true
by construction (κ0 cancelled, `R_U` never referenced); replaced (both engines)
with the falsifiable `z1/z0 − (κ1/κ0)·R_U == 0` (`R_U` defined independently from
ρ0,δU). (F2) `insufficient_verification` — the surviving exact product law
`R_target·M_mix = 8Λ(1−ε_W,split)/π²` was only `print`ed; added an `expect_zero`
in both engines. (+ F3 label-only docstring/subbanner.) **043 — ONE finding:**
(F2) `insufficient_verification` — the `M_supp` baseline check re-substituted the
SAME literal `B=8/π²` into both sides of a `B`-symbolic identity (could not
fail, never derived the number); replaced (both engines) with a check that
DERIVES `B=κ0²=(9/11)·σ`, `σ=88/(9π²) ⇒ 8/π²`. (+ F1 label-only docstring/
subbanner.) Both new-check sets PASS in both engines, engines agree, no
regression; `material_change=false` (assertions added/strengthened over
already-correct values, no verified RESULT value moved). The 9 label-only fixes
canonicalized unambiguous SymPy self-labels (docstring filenames/headers,
`.py` sub-stage `NN.k` indices, closing pass-lines, + 3 inline-comment self-refs
the arbiter grep caught — 040:136 `(section 23.2)`→`40.2`, 041:107 `section
24.1`→`41.1`) — every `.py` diff strip-the-number identical to HEAD; the `.wl`
source was already canonical on all 12 (stale `.txt` outputs cured by the
orchestrator re-run alone); same class as the I.2 stage-021 / II.1 032-036 fixes.
045 deferred-clean (source self-labels already canonical + outputs fresh; only
finding = 4 `Stage-27` CROSS-refs to upstream 044, deferred). No coverage count
moved (all 12 already verified+dual-engine from pass 1).

Numbering cross-refs (to OTHER stages) and one ambiguous self-vs-cross ref
(047:121 `the exact Stage-30 support-loading coefficient`) left UNTOUCHED and
deferred to `redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` (PENDING —
content-keyed, never offset-sweep). Reference memory `numbering-drift-root-cause`.

## Pass 2 — Batch III.2 (2026-06-05)

Pass-2 Batch III.2 (049-060, `Part III.2 — Tracking, zeta thresholds, asymmetry,
boost`, 12 stages): all 12 pass-2 re-verified at v2 depth + exhaustive
value-reconciliation augmentation; 12/12 verified, `material_change: false` on
all 12, 0 stop-cold, 0 blocked, 0 Codex deviations, all iter-1 exit 0. Every
emitted deliverable value reconciles — 126 values checked batch-wide, 0
MISALIGNED, 0 MISSING-DELIVERABLE (per stage: 049=10, 050=11, 051=9, 052=12,
053=11, 054=14, 055=9, 056=12, 057=7, 058=11, 059=10, 060=10). All 12 are
genuine dual-engine (both independent engines from pass 1, re-confirmed); no
new `.wl`, 0 sanctioned mirrors. **One checkpoint in range — 051 — re-verified
clean (no certified constant moved).** NO EM-projected stages in range.
**1 real script-side math finding→fixed→verified** (060); **11 label-only
self-label fixes→verified** (049-059, plus 060's F1 self-banner — every one of
the 12 was `verdict:findings`). **060 — ONE finding:** (F2) `tautological_check`
— the `wl:140` `Xi_micro` baseline assertion subtracted `xiMicro` from a verbatim
copy of its own `wl:132` definition (identically 0, could not fail); replaced
with `expectZero["Xi_micro chi-route equals D/M-route", xiMicroFromChi -
xiMicroFromDM]`, a genuine cross-check between two INDEPENDENTLY constructed
routes (susceptibility `chiSigma→1/theta` at `wl:133` vs Einstein/diffusion
`dSigma→mSigma*theta` at `wl:136`, where `mSigma` cancels). The value
`Xi_micro = Λφ²L²/(Θ T_X)` is UNCHANGED, deliverable still covered, 060.wl remains
an independent engine (not a transliteration); `material_change=false`. The 11
label-only fixes canonicalized unambiguous self-labels (NUMBER only, FORMAT
preserved — 2-digit docstrings/closing-prints kept 2-digit, 3-digit
filename-docstrings on 051/052 kept 3-digit, correct 2-digit `STAGE NN` banners
LEFT UNPADDED; 049 also fixed `.wl:93` closing `Stage 32`→`Stage 49`); all
`.py`/`.wl` diffs strip-the-number identical to HEAD except 060.wl's one math
line; same class as the I.2 stage-021 / II.1 032-036 / III.1 9-stage fixes. No
coverage count moved (all 12 already verified+dual-engine from pass 1).

Numbering cross-refs (to OTHER stages) and compound dual-epoch refs (050 `py:61`
`Stage 32`→stage049; 051 `py:20-21` `Stage 050/034` + `py:126`/`wl:87`
`Stage 047/030`; 055 `py:73` `Stage-35`→stage052; 056 `py:7` `Stage-36`→stage053;
059 `py:6` `Stage-41`→stage058 + `py:9`/`75` `Stage-39`→stage056; 060 `py:159`
`Stage-39`→stage056) left UNTOUCHED and deferred to
`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` (PENDING — content-keyed, never
offset-sweep). Reference memory `numbering-drift-root-cause`.

## Pass 2 — Batch III.3 (2026-06-05)

Pass-2 Batch III.3 (061-072, `Part III.3 — Microclosure, gain thresholds,
equilibrium, walls`, 12 stages): all 12 pass-2 re-verified at v2 depth (dual-engine)
+ exhaustive value-reconciliation augmentation; 12/12 verified, `material_change:
false` on all 12, 0 stop-cold, 0 blocked, 0 Codex deviations, all iter-1 exit 0.
Every emitted deliverable value reconciles — 116 values checked batch-wide, 0
MISALIGNED, 0 MISSING-DELIVERABLE (per stage: 061=12, 062=9, 063=8, 064=7, 065=11,
066=9, 067=8, 068=7, 069=9, 070=10, 071=13, 072=13). All 12 are genuine dual-engine
(both independent engines from pass 1, re-confirmed; 072.wl transliteration-leaning
but accepted per policy — a pure substitution+limit-identity stage); no new `.wl`,
0 sanctioned mirrors. **One checkpoint in range — 069 (`final_reduced_verdict`) —
re-verified clean at the higher bar (no certified constant moved; pins no numeric
constant — `Cres2Prim`/`Pres_gap` carried as FREE symbols, never re-asserted).**
NO EM-projected stages in range. **1 real script-side math finding→fixed→verified**
(070); **8 label-only self-label fixes→verified** (061,062,063,064,065,066,068,071);
**3 refresh-only→verified** (067,069,072). **070 — F1 `tautological_check`:** the
`I_1/J_1 = 4πa²ℓ` anchoring check self-cancelled in both engines (the `I_f/H_w`
factor cancels, independent of the sech moment value); made the moments load-bearing
— SymPy asserts `I_f=2/3`, Mathematica asserts `I_f≈2/3`/`I_g≈14/15` against the
already-computed NIntegrate values (`tol=10^-10`) — plus a corrected `wl:86` print
annotation `8/15`→`14/15`. (⚠️ The orchestrator false-positive guard caught the
audit's proposed `I_g=8/15` was WRONG — correct = 14/15 ≈ 0.9333, confirmed by the
30-digit NIntegrate; the script would have FAILED on `8/15`.) 070.wl remains an
INDEPENDENT engine (verify-confirmed); the symbolic `kappa`/`W_wall`/`Xi`
deliverables are UNCHANGED; `material_change=false`. The 8 label-only fixes
canonicalized unambiguous self-labels (NUMBER only, FORMAT preserved — 2-digit
prose docstrings/closing-prints kept 2-digit, 3-digit filename-docstrings on
068/070/071 kept 3-digit, correct 2-digit `STAGE NN` banners LEFT UNPADDED); all
`.py`/`.wl` diffs strip-the-number identical to HEAD except 070's three added math
lines; same class as the I.2 stage-021 / II.1 032-036 / III.1 / III.2 fixes. No
coverage count moved (all 12 already verified+dual-engine from pass 1). INFRA NOTE
(III.2 lesson re-confirmed): `exec-*` refreshes `exec_logs/` not the committed
`output/*.txt` — all 12 committed transcripts predated June-3 content additions;
the orchestrator re-ran all 24 engines (exit 0) and sed-refreshed every output, the
arbiter grep confirming no stale self-banner remains.

Numbering cross-refs (to OTHER stages), variable names, and ambiguous self-vs-cross
refs in `.py`/`.wl` source (063 `Stage-44`; 064 `Stage-45/46`; 065 `Stage-44`; 066
`Stage-48`; 069 `Stage-49/51` family + a 3-digit `Stage 049`; 070 `Stage-47/48` +
the `J1_stage48`/`J1Stage48` variable names) left UNTOUCHED and deferred to
`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` (PENDING — content-keyed, never
offset-sweep). Reference memory `numbering-drift-root-cause`.

## Pass 2 — Batch III.4 (2026-06-05)

Pass-2 Batch III.4 (073-084, `Part III.4 — Family-1 geometry, thresholds,
quadrupole`, 12 stages): all 12 pass-2 re-verified at v2 depth (dual-engine where
applicable) + exhaustive value-reconciliation augmentation; 12/12 verified,
`material_change: false` on all 12, 0 stop-cold, 0 blocked, 0 Codex deviations, all
iter-1 exit 0. Every emitted deliverable value reconciles — **119 values checked
batch-wide, 0 MISALIGNED, 0 MISSING-DELIVERABLE** (per stage: 073=5, 074=6, 075=11,
076=8, 077=9, 078=8, 079=6, 080=5, 081=11, 082=15, 083=21, 084=14). 11 stages are
genuine dual-engine; **084 is single-engine by design** (Mathematica-only
write-up/consolidation skeleton — card says "SymPy audit: none yet",
`is_status_only_candidate: true`; the absent SymPy is not a finding). No new `.wl`,
0 sanctioned mirrors. **NO checkpoints in range; NO EM-projected stages in range; NO
genuine `paper_misalignment` anywhere → ZERO paper/notes edits.**

**4 real script-side findings → fixed → verified:**
- **075 — F1 `tautological_check` (both engines):** `Theta_fail := Upsilon_fail/alpha_r²`
  made the round-trip `Upsilon_fail − alpha_r²·Theta_fail ≡ 0` (could never fail; the
  SymPy comment falsely claimed "not the trivial identity"). Removed both round-trips;
  added SymPy `expect_close` anchors of all 8 deliverables (Delta_0/Delta_inf,
  Upsilon/Xi/Theta fail+suff) against fixed EXTERNAL literals; the reduction lock
  `alpha_r² == 100` retained; the `.wl`'s pre-existing independent `expectApprox`
  battery stays.
- **083 — F1 `tautological_check` + F2 `insufficient_verification` (both engines):**
  `Delta := numer/denom` made `denom·Delta − numer ≡ 0` (could not catch a Cosh↔Sinh
  paste error), the `.wl` `A_F1 independent` check compared `(Pi/2)² ≡ Pi²/4`
  (identical operands), and the SymPy script only PRINTED the deliverable Pe/zeta
  windows. Replaced the SymPy Delta residuals with external-literal anchors, deleted
  the tautological `.wl` `A_F1` check, corrected the misleading "Independent BVP
  derivation" comments, and added the 4 Pe + 5 zeta window numeric asserts (F2).
- **081 — F1 `insufficient_verification` (SymPy):** the inversion `Q` was pinned only
  at `Q(0)`, `Q(1)` (two points don't pin a rational function). Added
  `expect_zero("Q-closedform", Q − (1 + zeta − 2·eps_blk·zeta)/(1 − eps_blk·zeta))`,
  matching the full identity the Mathematica engine already asserts.
- **077 — F1 `symbol_assumption_error` (SymPy):** the symmetric full-line integration
  variable `xi` (cut point `xi_* ≈ −0.3856 < 0`) was declared `positive=True`. Split
  the declaration so `xi` is `real=True` only (`alpha_r`/`lambda_mu` stay positive).
  Dormant latent trap; explicit `(-oo,oo)` bounds set the domain, so `I_f = 1/3` and
  all numeric values are byte-identical — no result moved.

**1 self-label fix → verified:** 074 — SymPy module docstring filename `stage57` →
`stage074` (matches the actual 3-digit filename; banner already canonical) + output
refresh. **1 refresh-only → verified:** 076 — `stale_output` (committed banners
`STAGE 59`/`STAGE 059` → `STAGE 076`; source canonical, no edit). **1 deferred-clean
→ verified:** 080 — the audit flagged five comment/docstring labels, but ALL are
CROSS-references (`Stage-61/62` → stages 078/079), not stage-080 self-labels (its own
self-labels are canonical); orchestrator deferred per Reading-2 (the auditor
over-flagged — the same class was correctly marked CLEAN on sibling 078). **5 clean →
verified:** 073, 078, 079, 082, 084.

All `.py`/`.wl` source diffs are strip-the-number identical to HEAD except the genuine
math additions on 075/077/081/083 (de-taut removals + anchor/assert insertions); no
deliverable result value moved on any stage (confirmed by the committed-output diff:
only banner refreshes, removed tautological-check lines, and new passing-anchor lines).
No coverage count moved (all 12 already verified+dual-engine-or-single-by-design from
pass 1). INFRA NOTE (III.2/III.3 lesson re-confirmed): `exec-*` refreshes `exec_logs/`
not committed `output/*.txt`; the orchestrator re-ran all 23 engines (exit 0) and
sed-refreshed every output, the arbiter grep confirming no stale self-epoch (NNN−17)
banner remains. (073's SymPy engine is the `..._sympy_audit_refresh.py` variant; its
committed output is `..._sympy_audit_refresh.txt` — present, git-tracked, fresh; the
MANIFEST `sympy_output` glob keys on the non-refresh name so it shows `None`, harmless.)

Numbering cross-refs (to OTHER stages) in `.py`/`.wl` source — 078 docstring
`Stage-60`/`Stage-58`; 080 docstring+comments `Stage-61/62` (×5); 081 comments
`Stage-35`→052 / `Stage 63`→080 — left UNTOUCHED and deferred to
`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` (PENDING — content-keyed, never
offset-sweep). Reference memory `numbering-drift-root-cause`.

## Pass 2 — Batch III.5 (2026-06-05)

Pass-2 Batch III.5 (085-090, `Part III.5 — Quadrupole cancellation, loading ratio,
verdict`, 6 stages): all 6 pass-2 re-verified at v2 depth (both engines present on all
6) + exhaustive value-reconciliation augmentation; 6/6 verified, `material_change: false`
on all 6, 0 stop-cold, 0 blocked, 0 Codex deviations beyond one sanctioned helper-add,
all iter-1 exit 0. Every emitted deliverable value reconciles — **59 values checked
batch-wide, 0 MISALIGNED, 0 MISSING-DELIVERABLE** (per stage: 085=13, 086=13, 087=5,
088=8, 089=14, 090=6). All 6 are genuine dual-engine. No new `.wl`, 0 sanctioned mirrors.
**TWO checkpoints in range (089 & 090) — BOTH cleared the higher bar; NO EM-projected
stages in range; NO genuine `paper_misalignment` anywhere → ZERO paper/notes edits; NO
new postulated constant.**

**3 real script-side findings → fixed → verified:**
- **087 — F1 `insufficient_verification` (both engines):** the "upstream stage-086
  cross-check" asserts compared each window literal (`rho_suff/rho_fail/rho_max`)
  against the SAME literal re-typed in the same file (py:58 vs py:73; wl:57 vs wl:61) —
  a hollow self-comparison (both sides move together, can't detect a mistyped literal).
  Replaced the three self-comparisons with can-fail structural relations
  (`rho_suff < rho_fail`, `rho_fail < rho_max`, `0 < rho_max − rho_fail < 1e-6`; gap =
  9.671731e-8) + reworded the overclaiming comments; re-anchored the `.wl` `zeta_*`
  numeric checks (wl:73-75) to `rho_X − 1` so they genuinely test the `epsBlk→0`
  substitution of `zetaReq`. iter-2 reworded the same overclaim surviving in the SymPy
  top docstring (non-printing; transcript byte-identical). The three literal VALUES
  unchanged; `.wl` verify-confirmed STILL INDEPENDENT (not a transliteration).
- **088 — F1 `stale_output`/stale self-label (`.py` docstring only, 0-seat):** the
  pre-renumber OWN-number self-label py:3 `...stage71_...` + py:5 `Stage 71.` (088−17 =
  071 EM-extension drift) → `stage088`/`Stage 088`. Banner/filename/card/`.wl` banner
  all already canonical; `stage085` cross-refs (py:112, wl) correct → untouched.
  Committed outputs byte-identical (docstrings not printed). (Both first-pass
  fragilities re-confirmed ALREADY FIXED, no new finding: the `omega**2 → u` subs lands
  on atomic `omega**2`; the `stage 085` upstream ref carries no `*)` to close the
  comment early — assertion count = PASS-line count, no silent partial run.)
- **089 — F1 + F2 `tautological_check` (CHECKPOINT, higher bar; de-tautologized BOTH,
  not deleted):** F1 — `Q` baked with `eps_blk = 0` made `expect_zero(Q − (1 + zeta))`
  a pure X−X self-cancel; fix made `eps_blk` symbolic, introduced general
  `Q_gen = (1+(1−2 eps_blk) zeta)/(1 − eps_blk zeta)`, and asserted the `eps→0`
  reduction on the GENERAL form (a structural transcription error in `Q` now fails),
  plus a parallel `.wl` reduction assert for engine symmetry. F2 — the boxed `Pe_req=0`
  was verified `0==0`; fix removed the self-check, replaced with a CAN-FAIL positivity
  assertion on the zero-bias success margin `zeta_F1(0) − zeta_min` (= A_F1 − 1/3 ≈
  0.6667, the quantity that FORCES `Pe_req=0`), then printed `Pe_req=0` as the
  consequence. ⭐ ORCHESTRATOR FALSE-POSITIVE/SAFETY GUARD: the audit's drafted
  `sp.Piecewise((0,cond),(sp.nan,True)) → expect_zero` gate was FRAGILE (a failed
  precondition passes SILENTLY since `abs(complex(sp.nan)) > tol` is False) → the
  orchestrator rewrote the directive to the explicit-raise margin form; Codex applied
  it correctly. **SANCTIONED Codex deviation:** the `.wl` had no `expectZero` helper
  (only `expectTrue`/`expectApprox`) → Codex added a minimal one
  (`FullSimplify[Together[Expand[expr]]]`, pass iff `=== 0`); verify-confirmed correct,
  `.wl` STILL INDEPENDENT. After the fixes 089 carries no remaining named tautology,
  both engines substantive, paper alignment exact (boxed `Pe_req = 0` Output retained);
  the Mathematica side independently re-derives the upstream `Pe` via `FindRoot` — a
  robust route, NOT the latent `nsolve`-near-`tan` **pitfall #10** (explicitly avoided).

**3 clean → verified:** 085, 086, 090 (090 = checkpoint, clean, cleared the higher bar
via a substantive own-assertion `zeta_req = rho_alpha − 1` + branch-ordering
inequalities).

All `.py`/`.wl` source diffs are strip-the-number identical to HEAD except the genuine
math additions on 087/089 (de-taut removals + can-fail anchor/assert insertions + the
sanctioned 089 `.wl` helper-add); no deliverable result value moved on any stage
(confirmed by the committed-output diff: only removed tautological-check lines and new
passing-anchor lines; 088's outputs byte-identical). No coverage count moved (all 6
already verified+dual-engine from pass 1). INFRA NOTE (III.2/III.3/III.4 lesson
re-confirmed): `exec-*` refreshes `exec_logs/` not committed `output/*.txt`; the
orchestrator re-ran all 12 engines (exit 0) and sed-refreshed every output, the arbiter
grep confirming no stale self-epoch (NNN−17, 068-073 band) banner remains on any of the
6 (only canonical `STAGE 0{85..90}` banners + the canonical self paper-eq ref
`app-stage089-Pe-zero`).

Numbering cross-refs (to OTHER stages) in `.py`/`.wl`/notes source — 086 `py:37`
`Stages 63-64`→080/081; 089 comments cross-ref `Stage-62`/082/075/074; 090 notes
provenance-tracker `Stage 73` (out of audit scope) — left UNTOUCHED and deferred to
`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` (PENDING — content-keyed, never
offset-sweep). (088's `stage085` refs are correct — no defer.) Reference memory
`numbering-drift-root-cause`.

## Pass 2 — Batch IV.1 (2026-06-05)

Pass-2 Batch IV.1 (091-102, `Part IV.1 — Grouped P_2 geometry, decoupling, contamination`,
12 stages): all 12 pass-2 re-verified at v2 depth + exhaustive value-reconciliation
augmentation; 12/12 verified, `material_change: false` on all 12, 0 stop-cold, 0 blocked,
all iter-1 exit 0. Every emitted deliverable value reconciles — **99 values checked
batch-wide, 0 MISALIGNED, 0 MISSING-DELIVERABLE** (per stage: 091=8, 092=13, 093=5, 094=10,
095=6, 096=6, 097=11, 098=8, 099=7, 100=12, 101=6, 102=5). **ONE checkpoint in range (096)
— cleared the higher bar; NO EM-projected stages in range; NO genuine `paper_misalignment`
anywhere → ZERO paper/notes edits; NO new postulated constant.** **093 is status-only /
Mathematica-only by design** (`is_status_only_candidate`; the absent SymPy is not a finding).
No coverage count moved (all 12 already verified+dual-engine from pass 1).

**4 real script-side findings → fixed → verified:**
- **096 — F1 `insufficient_verification` (CHECKPOINT, higher bar):** SECTION II evaluated
  the carried Stage-092 obstruction formula `c_pole=(1+eps_4)/(4(1+eps_2)^2)` ONLY at the
  degenerate static point `eps_2=eps_4=0`, so its eps-structure was never exercised. Fix
  (both engines): `eps_2`/`eps_4` free symbols + general formula + two CAN-FAIL off-static
  probes (`eps_4=1,eps_2=0 → 1/2`; `eps_2=1,eps_4=0 → 1/16`, both ≠ 1/4) + the static limit
  taken FROM the general formula. SECTION I (l=0⊥l=2 orthogonality, Laplace eigenvalue
  ℓ(ℓ+1)=6) was already substantive, unchanged. All deliverables (1/4, 3/4, 4/3, 1/3, Yhat)
  preserved. Cleared the checkpoint higher bar.
- **093 — F1 `insufficient_verification` (status-only, Mathematica-only):** the SAME
  obstruction-formula de-taut — a symbolic-eps anchor block added to the `.wl`; original four
  deliverable checks byte-identical; stays single-engine by design.
- **094 — F1 `tautological_check` (both engines):** `K_g2`/`K_g4` were bare literal 0 (the
  comment falsely claimed "established by orthogonality"); now the ACCUMULATED (proven-zero)
  l=0↔l=2 overlap moments with can-fail asserts + `c_pole`/`c_geom` pinned individually.
- **100 — F1 `mathematica_transliteration` ⭐ ORCHESTRATOR OVERRIDE (clean→findings):** the
  audit agent called 100 clean on `sp.im(...)` vs `.../I`; the orchestrator's ground-truth
  read OVERTURNED that (the `.wl` was a line-by-line port of the `.py` — same `Series[]`/
  `series()` black box on the same rational `Y`, same coefficient extraction). Per the
  dual-engine rule (user-level call) the **user AUTHORIZED the independent-route rewrite
  2026-06-05** (pass-1 had declined it → `blocked_legitimate`). Codex re-derived `K2`/`K4`/
  `Gamma5` via an INDEPENDENT analytic geometric-series expansion (`1+u+u²` through ω⁵,
  hand-collected order coefficients, NO `Series[yRet,…]` on the full rational); verifier
  confirmed genuine independence, all 4 deliverables still pass via can-fail checks, `chiQ`
  stays a free symbol. **Committed Mathematica output BYTE-IDENTICAL to HEAD** (method
  changed, no emitted value did); `material_change: false`; reference SymPy `.py` untouched.
  **100's `.wl` is now an INDEPENDENT engine** (see `MATHEMATICA_MIRROR_POLICY.md`).

**8 clean → verified:** 091, 092, 095, 097, 098, 099, 101, 102.

The recurring de-taut theme this batch — the obstruction-formula-evaluated-only-at-the-
degenerate-point trap — hit 093/094/096 (sibling 092 was correctly CLEAN, already exercising
the eps-dependence via the first-order expansion). All `.py`/`.wl` source diffs are
strip-the-number identical to HEAD except the genuine math additions on 093/094/096 (de-taut
removals + free-symbol general formula + can-fail off-static probes) and 100's `.wl`
de-transliteration; no deliverable result value moved on any stage. INFRA: `exec-*` refreshes
`exec_logs/` not committed `output/*.txt`; the orchestrator re-ran the 094/096/100 SymPy +
093/094/096/100 Mathematica engines (all exit 0; 094 SymPy and BOTH 100 outputs
byte-identical post-refresh), arbiter grep on all 12 committed outputs = CLEAN.

Numbering CROSS-refs (to OTHER stages) deferred to `redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`
(PENDING — content-keyed, never offset-sweep): 095 notes legacy "Stage 74/75/77/78" narrative;
096 docstring "Stage 092" vs notes "Stage-75" (source cross-ref, not a 096 self-label).
Reference memory `numbering-drift-root-cause`.

## Pass 2 — Batch IV.2 (2026-06-06)

Pass-2 Batch IV.2 (103-114, `Part IV.2 — Outgoing DtN, deformation, robustness, robin`,
12 stages): all 12 pass-2 re-verified at v2 depth + exhaustive value-reconciliation
augmentation; 12/12 verified, `material_change: false` on all 12 (NO downstream staling), 0
stop-cold, 0 blocked, 0 Codex deviations, all iter-1 exit 0. Every emitted deliverable value
reconciles — **86 values checked batch-wide, 0 MISALIGNED, 0 MISSING-DELIVERABLE** (per stage:
104=11, 105=8, 106=11, 107=9, 108=14, 109=4, 110=6, 111=9, 112=7, 114=7; plus status-only 103
[symbolic carry-forward] and 113 [8 notes-stated] reconciled). **TWO checkpoints in range (105,
112) — BOTH cleared the higher bar; NO EM-projected stages in range; NO genuine
`paper_misalignment` anywhere → ZERO paper/notes edits; NO new postulated constant.** **103 and
113 are status-only by design** (no scripts). No coverage count moved (all 10 dual-engine stages
already verified+dual-engine from pass 1; 103/113 stay status-only single-engine).

**4 real script-side findings (5 findings) → fixed → verified:**
- **105 — F1 `insufficient_verification` + F2 `mathematica_transliteration` (CHECKPOINT, higher
  bar):** F2 was an OVER-CALL — the `.wl` retarded half already used `Im[]`-projection +
  `Reduce`/`ToRules` + a factored form and the deformed branch was independent — but BOTH engines
  shared the F1 tautology: they matched the canonical ω⁵ coefficient to a HARDCODED fingerprint
  target `a^5/(27 c_s^5)` (= `sigma_can/4` retyped), never evaluating the actual outgoing l=2
  Hankel DtN fingerprint. DE-TAUTOLOGIZED (both engines): each now DERIVES
  `Lambda_2^out = z d/dz ln h_2^(1)(z)` from the spherical Hankel function via VISIBLY DISTINCT
  constructions — SymPy from `j_2 + i*y_2` (spherical Bessel) closed form, Mathematica from native
  `SphericalHankelH1[2,z]` — asserts the can-fail series `-3 + z²/3 + z⁴/9 + i z⁵/9`, reads the
  DERIVED `Y_2^out` imag z⁵ coefficient (`= 1/27`), and forces `chi_Q=1` by matching the retarded
  ω⁵ coefficient to that DERIVED target (no typed literal on the RHS). `chi_Q=1` UNCHANGED; F2
  dissolves under the distinct per-engine constructions. Cleared the checkpoint higher bar.
- **107 — F1 `mathematica_transliteration` (pass-1 MISSED):** the `.wl` was a full line-by-line
  port (`Series[l0/lambdaDef,{z,0,5}]` on the same normalized ratio). USER-AUTHORIZED re-author:
  the `.wl` now reaches the deformed branch / `chi_Q` via an order-by-order undetermined-coefficient
  `Solve` of `Lambda_def*Y = L0` (`CoefficientList` + `Solve`, unique-solution guard); all `Series`
  removed; verify-confirmed independent; MMA byte-identical.
- **110 — F1 `mathematica_transliteration` (pass-1 MISSED):** full port (`Series[(-3+rho)/lambdaR]`).
  USER-AUTHORIZED re-author to an undetermined-coefficient linear `Solve` of `lambdaR*Y_R = (-3+rho)`
  + a separate rho-jet solve for the linearization; all `Series` removed; verify-confirmed
  independent; MMA byte-identical.
- **114 — F1 `mathematica_transliteration` (pass-1 MISSED):** full port (`Inverse[m]`/`apart` on the
  same matrix). USER-AUTHORIZED re-author to the Schur complement `delta_Lambda(D)` via explicit
  scalar `Solve`-elimination of the 2×2 core system (`Solve` for s, back-sub, `Solve` for q, form
  `g_s*s+g_q*q`); `Inverse[m]` removed, matrix object gone; verify-confirmed independent; MMA
  byte-identical.

**Checkpoint 112's pre-existing `.wl` confirmed INDEPENDENT with NO re-author** (native
`Reduce[presCond==0 && sigma!=0, gamma]` reconstructing the Stage-92 `(b,a0,a5)` to force
`gamma_W=1/9`, a route the `.py` lacks); cleared the checkpoint as-is. 104/106/108/109/111 `.wl`
confirmed independent, unchanged.

**8 clean → verified:** 103 [status-only], 104, 106, 108, 109, 111, 112 [checkpoint], 113
[status-only].

**Dual-engine / independence outcome:** 4 newly-independent `.wl` this batch (105 re-grounded +
107/110/114 re-authored), 0 sanctioned mirrors — added to the Independent-Mirror Set (see
`MATHEMATICA_MIRROR_POLICY.md`). All `.py`/`.wl` source diffs are strip-the-number identical to
HEAD except the genuine de-taut/re-author additions on 105/107/110/114; no deliverable result
value moved on any stage. INFRA: 20 exec runs exit 0 (10 dual-engine × 2 engines; 103/113 have no
scripts); `exec-*` refreshes `exec_logs/` not committed `output/*.txt` → the orchestrator
sed-refreshed every committed `.txt` (107/110/114 + the 5 clean dual-engine MMA outputs
byte-identical; 105 additive [fingerprint lines]; **108 output NORMALIZED** — refresh stripped a
stray `# exit_code: 0` trailer a prior commit had baked in, deliverables identical, NOT a math
change); arbiter grep on all 10 committed outputs = CLEAN (no stale self-epoch NNN−17 = 086–097
banner; the only stage-label hit is 106's correct upstream cross-ref "Carry-in chi_Q = 1 from
stage 105"). Seat policy held: the 4 `.wl`-touching Codex sessions ran in 2 waves of 2 (105∥114,
then 107∥110) under the flock; orchestrator exec sequential after all Codex done.

**Side note (OUT OF SCOPE, awareness only):** `notes/stages/review/stage_103_review.md` is a
stale first-pass review artifact whose body reviews **Stage 035, not 103** (flagged by the 103
audit agent; not a script-audit item) — logged as a residual in `PAPER_CLEANUP_TRACKER.md` (P5-11)
to clean later. Numbering CROSS-refs: none new beyond the correct 106→105 carry-in. Reference
memory `numbering-drift-root-cause`.

## Pass 2 — Batch IV.3 (2026-06-06)

Pass-2 Batch IV.3 (115-126, `Part IV.3 — Core balance, DtN mixed, outlet, positive source`,
12 stages): all 12 pass-2 re-verified at v2 depth + exhaustive value-reconciliation augmentation;
12/12 verified, `material_change: false` on all 12 (NO downstream staling), 0 stop-cold, 0 blocked,
0 Codex deviations, all iter-1 exit 0. **108 deliverable values checked batch-wide, 1 MISALIGNED
(resolved — the r_F1 appendix surd below), 0 MISSING-DELIVERABLE** (per stage: 115=9, 116=9, 117=16,
118=15, 119=6, 120=6, 121=7, 122=12, 123=6, 124=7, 125=6, 126=9). **NO checkpoints in range; NO
EM-projected stages in range; no new postulated constant.** **120 & 124 are legitimately status-only
by design** (notes-only, both engines null, like 103/113) — no coverage count moved (the 10
dual-engine stages stay verified+dual-engine; 120/124 stay status-only). No downstream staling.

**Published-paper value typo RESOLVED (Codex-applied out-of-band, Claude-reviewed):** the Family-1
radius surd `r_F1` was published with a `117π²` typo (`√(4107−117π²)`→`√(4107−100π²)` in
`paper/appendices/stage_appendix_part04.tex:562` + `paper/parts/part04_geometry_retarded_mouth.tex:576`);
`100π²` is arithmetically forced and matches every script/note; flagged by 122 F1 / 123 F1; the
SCRIPTS were always correct → 121/122/123 are script-clean. Pass-2 caught a published-paper arithmetic
typo the first pass MISSED.

**2 script-side findings → fixed → verified (both `material_change: false`):**
- **117 — F1 `mathematica_transliteration` → FULL re-author:** the orchestrator ground-truth
  `.wl`-vs-`.py` read found the `.wl` was a full line-by-line `Series[]`-based transliteration across
  ALL SIX sections (the audit agent scoped the directive to §5; the orchestrator broadened it).
  USER-AUTHORIZED re-author to an independent route — §1–§4 via an undetermined-coefficient
  `jetFromBalance` solve (`den·Y−num==0` order-by-order, NO `Series[ratio]` survives), §5 via an
  explicit 2×2 `coreMatrix` `Solve`-elimination + Schur complement with a `coreSchurResidual===0`
  guard (re-typed `sigmaC` literal GONE); only the permitted final §5 residual `Series` remains.
  **117 GAINED a re-authored independent dual-engine `.wl`** (committed output byte-identical to HEAD,
  method-only); added to the Independent-Mirror Set (see `MATHEMATICA_MIRROR_POLICY.md`).
- **118 — F1 `tautological_check` → de-taut:** the "K_q closed form" check was X−X; `K_q`/`kQ` now
  tied to the independently-computed gradient integral `chiGrad = ∫(χ')²dz` (asserted `=π²/(4 L_W²)`)
  so the closed-form check is load-bearing; printed `K_q` value unchanged; both engines exit 0; `.wl`
  stays independent.

**121/122/123 retro-sweep `.wl` re-confirmed genuinely INDEPENDENT** (Reduce/Solve-with-branch-guard,
all using the correct `100π²`; 123's negative `Xi_v(F1)≈−1.01675633282526` preserved). **8 script-side
clean → verified:** 115, 116, 119, 121, 122, 123, 125, 126 (all dual-engine, independence confirmed).

**Dual-engine / independence outcome:** 1 newly-independent `.wl` this batch (117 re-authored), 0
sanctioned mirrors. All `.py`/`.wl` source diffs are strip-the-number identical to HEAD except the
genuine re-author/de-taut additions on 117/118; no deliverable result value moved on any stage — so
NO downstream staling. INFRA: 20 exec runs exit 0 (10 dual-engine × 2 engines; 120/124 have no
scripts); `exec-*` refreshes `exec_logs/` not committed `output/*.txt` → the orchestrator
sed-refreshed every committed `.txt` (117 + 118 + the 7 clean dual-engine MMA outputs byte-identical;
**116 NORMALIZED** — refresh stripped a stray `# exit_code: 0` trailer a prior commit had baked into
BOTH committed `.txt`, deliverables identical, NOT a math change; the historically-flaky 116 ran
clean/deterministic). Arbiter grep on all 10 committed outputs = CLEAN of stale self-epoch (−17 band
098–109) self-banners (only hit = 116's `gamma0_bare (upstream-carried input, Stage 98)` γ₀-provenance
CROSS-ref → DEFERRED, content-keyed). Seat policy held: 117 + 118 = 2 `.wl`-touching Codex sessions
(concurrent, at the 2-seat cap, flock-safe); the r_F1 paper fix = a separate out-of-band 0-seat Codex
session.

**Deferred (numbering pass; content-keyed, NEVER offset-sweep):** **⚠️ DISCOVERY** — `\stagefield{Purpose}`
card self-labels drift +17 (117 "Stage 134", 119 "Stage 136", 120 "Stage 137", 124 "Stage 141"), a
class the numbering reconciliation MISSED (its scan keyed on `\section`/`\label`, NOT
`\stagefield{Purpose}`); 122/123 cards "Mathematica audit: none yet" status understatement; 121 `.wl`
`stage99TubeLength`; 116 source "Stage 98" γ₀ cross-ref. **Ansatz catalog:** γ₀=(1+r_c)/9 re-confirmed
a POSTULATED pure-scale ANSATZ (not derived). Reference memory `numbering-drift-root-cause`.
