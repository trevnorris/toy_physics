# Checkpoint Constant Provenance Audit

This document records constant-provenance findings for the checkpoint stages in
`CITATION_SUPPORT_SET.md`.

The goal is narrow: make sure the checkpoint audits do not hide unexplained
literals behind apparently passing CAS scripts.

> **Numbering — SCRIPT/OUTPUT band: COMPLETE 2026-06-10.** The dedicated
> content-keyed, label-only stage-number reconciliation of the script/output
> layer (`.py`/`.wl`/`.txt`) is finished — 3 bands committed (`b25cb57` 001–090,
> `7a94cbc` 091–180, `857e255` 181–253). Digit-only label/banner/string edits,
> **zero verification-substance change** (no constant introduced or moved;
> strip-digits proof per band; both engines re-run exit 0). Plan:
> `redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`.

Snapshot date: `2026-06-03` (**batch VIII.1 close {243–253} — `Part VIII.1 — Relaxed branch, dynamic event chain, cold survival`, the FINAL forward first-pass batch under the v2 paper-grounded auditor WITH the dual-engine rule.** First-pass paper-grounded audit on stages 243–253. **Checkpoints `243`, `248`, and `253` ALL returned `verified` at the higher bar with `material_change: false`** (and on all 11 stages). **ALL THREE are ROUTE/COVERAGE changes, NOT constant changes** — each pre-existing line-by-line transliteration `.wl` was RE-AUTHORED to a genuinely independent route, landing on the SAME values. **243** (relaxed-branch lift / leakage-work lane / non-rigid solve / recovery / short-range firewall): IBP closure / native `LinearSolve` / `TrigExpand` / `Series`-at-∞, with the hardcoded `expected*` residuals DELETED. **248** (energy conservation / threshold speeds / Coulomb reference / near-top action): the §II `Solve`-mirror replaced by a native SATISFACTION route (compiler closed forms verified to satisfy their defining energy equalities via substitution + `FullSimplify` + non-vacuity guard + positive-branch guard); this was the batch's one iter-2 (iter-1's `Reduce`/`ToRules` route was a Wolfram-version dead end Codex correctly BLOCKED, orchestrator-reframed — a Claude math-coverage resolution, NON-conceptual, no paper edit). **253** (lattice-turnover / calibration recovery / stiffness map / temperature ceiling / screening ratios): native `D[Log[V[r]]]` + 5 `Solve` energy/force balances + regrouped threshold / Korringa / screening blocks. **No checkpoint constant was introduced or moved at 243, 248, or 253;** the Stage 243/248/253 audit notes below are updated to record the VIII.1 re-authoring (243 retains its `-Sqrt[2]/4` / `Sqrt[2]Sqrt[Pi]/8` / `chi_lam/k_V` / quadratic-vertex forms; 248 retains its symbolic `1/E` / `Pi` near-top-action forms + benchmark-only labels; 253 retains its `Υ_lat` calibration symbol + benchmark-only readbacks). **5 NOTES-ONLY paper_misalignment numerical typos were USER-RESOLVED this batch (correct-to-script; published cards UNAFFECTED; each cross-engine or internally corroborated):** 247 notes:406 Δ `210.17750000→142.17750000` (notes' own 9·16−1.35²=142.1775 + adjacent D0=3.76481862); 253 notes:274 benchmark `187.23361317→119.23361317` (notes' own 65.45193926/0.5489386551=119.2336 + cross-engine `.wl`) AND notes:419 a_int `10.95423247→10.95423248` (=4·K_turn); 244 notes:366 `196√2→128√2` (script E0=16 structure); 248 notes:506 `×168%→×100%` (notes' own 23.3128% + script ×100; recurs the stale "168" at 148/232). The audit agents initially MISATTRIBUTED 247/253 to the published cards but the orchestrator verified the cards are clean — all 5 are NOTES-ONLY. No derived or carried checkpoint constant MOVED, so **cumulative checkpoint-constant provenance is otherwise unchanged from the VII.2 close.** **MILESTONE — with VIII.1 closed, the first end-to-end red-team pass is COMPLETE (253/253 = 100%).** Previous VII.2 entry retained below.)

Snapshot date prior: `2026-06-03` (**batch VII.2 close {231–242} — third batch of Part VII under the v2 paper-grounded auditor WITH the dual-engine rule, theme "Rigid-mouth orbit-lock / branch-dressing / twin-support".** First-pass paper-grounded audit on stages 231–242. **Checkpoints `239` and `242` BOTH returned `verified` at the higher bar with `material_change: false`** (and on all 12 stages). **BOTH 239 and 242 are ROUTE/COVERAGE changes, NOT constant changes** — each pre-existing line-by-line transliteration `.wl` was RE-AUTHORED to a genuinely independent route. **239** (rigid-mouth physical normal form / Cartesian orbit-lock): forward Jacobian of the boxed dependent vector + native `PseudoInverse` left-inverse + `Reduce`/`Equivalent` orbit-lock (vs the `.py`'s backward hardcoded `SrmDep`). **242** (twin-support strict inclusion): `Resolve[ForAll,Reals]` strict-inequality certificate + `D[]` on real closed forms (vs the abstract-ζ device) + `logDrift` total-log-differential (vs `Exp[t·d]`) — so the load-bearing twin-window strict inclusion `C_mix < Pi_tr < 2 C_mix` (the ρ_α-style ratio 4/3 region: `Pi_tr = (4/3) C_mix`, already carried at the 242 entry below) is now tested STRICTLY on BOTH engines. No checkpoint constant was introduced or moved at 239 or 242; the Stage 239 and Stage 242 audit notes below are updated to record the VII.2 re-authoring (242 retains its carried `2/11`, `4/3`, `8/π²`, `2/3`, beta-window, and probe-only labels; 239 retains its symbolic `I_2` / `(0,-V,U-V)` / `(0,1,1)` / `U=0,V=0`). **The VII.2 load-bearing constants worth logging for the no-magic-numbers record (all NON-checkpoint stages, tracked here for completeness, each NOW DUAL-ENGINE-ANCHORED — independently recomputed by the new `.wl` and corroborating the verified SymPy script; the published cards/appendices were UNAFFECTED, carrying abstract forms, so no published value moved):** **240** `ρ_α = 4/3`, `ζ_req = 1/3`, `Pi_tr = (4/3) C_mix`; **241** the `ϱ` windows `1/3, 125/369, 2/3, 250/441`; **232** the figure-of-merit prefactor `100` (`Ξ` block). The 3 notes-only numerical/coefficient typos corrected this batch (231 `240·δ²ξ→189·δ²ξ` / `189·ξ³→121·ξ³`; 232 figure-of-merit prefactor `168→100`; 241 `ϱ_WΛ` upper bound `193/369→125/369`) were corrected to the already-correct SymPy scripts (Codex-applied, Claude-reviewed), each cross-engine corroborated by the new `.wl` — no derived or carried constant MOVED, so **cumulative checkpoint-constant provenance is otherwise unchanged from the VII.1 close (221 retained route-only; 218/203 retained, 162 retained at 217/218).** Previous VII.1 entry retained below.)

Snapshot date prior: `2026-06-02` (**batch VII.1 close {219–230} — second batch of Part VII under the v2 paper-grounded auditor WITH the dual-engine rule, theme "Mixed-bundle / resonance / branch-packet".** First-pass paper-grounded audit on stages 219–230. **Checkpoint `221` returned `verified` at the higher bar with `material_change: false`** (and on all 12 stages). **221 is a ROUTE/COVERAGE change, NOT a constant change** — its pre-existing line-by-line transliteration `.wl` was RE-AUTHORED to an independent route (native `D[QPi/DeltaPi,portPi]` derivative, `Residue`, `ComplexExpand`, uncollapsed Breit–Wigner), F1 survival round-trips de-tautologized, deliverable #9 (linear survival window) brought into genuine dual-engine coverage; one sanctioned Codex deviation (F3) used the native derivative form, reconciling the Stage-220 identity `∂_Π D_Π = −N` — verified correct. No checkpoint constant was introduced or moved at 221; the 221 entry below remains symbolic with only probe-only numeric labels (the Stage 221 audit note below is updated to record the VII.1 re-authoring). **The VII.1 load-bearing constants worth logging for the no-magic-numbers record (all NON-checkpoint stages, tracked here for completeness, each NOW DUAL-ENGINE-ANCHORED — independently recomputed by the new `.wl` and corroborating the verified SymPy script; the published cards/appendices were UNAFFECTED, carrying abstract forms, so no published value moved):** the corrected upper-/wall `R_Q` figures `145.483858657863` (222, λ_W=0.2), `138.814136942081` / `137.502546600713` (223, λ_W=0.2); the `i=h` rigidity determinant factor `200+147π²` (227); the δ_1 coefficient `196π²/(98π²−25)` and reduced-det `196(200+147π²)` (228); the crossover-cubic leading coeff `121ξ³` (229); and the 230 thresholds `R_*≈1.229255438463336` / `δ_*≈0.723111617875019`. The seven notes-only numerical typos these corrected (`213.483858657863`/`206.814136942081`/`205.502546600713`/`251+215π²`/`247π²/(98π²−25)`/`247(251+215π²)`/`189ξ³`) were a systematic +68/+51 additive family corrected to the already-correct SymPy scripts (Codex-applied, Claude-reviewed) — no derived or carried constant MOVED, so **cumulative checkpoint-constant provenance is otherwise unchanged from the VI.1 close (218/203 retained, 162 retained at 217/218).** Previous VI.1 entry retained below.)

Snapshot date prior: `2026-06-02` (**batch VI.1 close {201–218} — first batch of Part VI under the v2 paper-grounded auditor WITH the dual-engine rule, theme "Explicit realization, scalar slice, ray ranking".** First-pass paper-grounded audit on stages 201–218. **Checkpoints `203` and `218` BOTH returned `verified` at the higher bar with `material_change: false`** (and on all 18 stages). **The load-bearing checkpoint constant recorded this batch is the per-envelope lifted Bézout bound `162 = 3⁴·2` (stages 217/218).** It is **arithmetically forced** — the only value consistent with the downstream budgets `2 × 162 = 324`, `1140 + 324 = 1464`, the fallback `2 × 750 = 1500 → 2640`, and the projected-chart per-envelope bound `750 = 5·5·5·6`. **The SCRIPT was always correct (162); the WRONG typos `179` and `230` lived only in the PUBLISHED card (`stage_217.tex`), the Part-VI appendix (`stage_appendix_part06.tex`, `eq:app-part06-five-bezout` + the part06 stage-table row), and the 217/218 notes — all corrected to 162 this batch (Codex-applied, Claude-reviewed). No derived or carried constant MOVED** — this is a paper/notes typo correction aligning the prose to the already-correct script, so **cumulative checkpoint-constant provenance is otherwise unchanged from the V.3 close (200 retained, 105 retained at the higher-bar standard).** Provenance now anchored at the 203/218 checkpoints: `162 = 3⁴·2` (lifted per-envelope Bézout bound), `324 = 2×162` (per-envelope support-5 budget), `1464 = 1140 + 324` (total support-≤5 ledger), `2640 = ... + 2×750` (projected/fallback total), `750 = 5·5·5·6` (projected-chart per-envelope bound); `1140` is the carried Stage-215 support-`≤4` global ledger budget. The Stage 218 provenance entry below is updated (the stale `179` corrected to `162`). The other 16 stages (201, 202, 204–217) are NOT checkpoints. Previous retro-sweep entry retained below.)

Snapshot date prior: `2026-06-01` (**retro-sweep {121, 122, 123} — dual-engine retrofit of already-`verified` SymPy-only stages.** Each gained a NEW independent-route `.wl` (Codex wrote, Claude reviewed); all 3 re-verified, `material_change: false`. **The retro-sweep introduced NO new postulated or derived constant/literal** — the new `.wl` REUSE the existing source-anchored values (121's `r_geom`/`r_F1`/`r_c` from the Stage-99 length law; 122's `g_±` solved from the compensation quadratic with `cStage` a free positive symbol; 123's `Xi_v`/`Xi_T` from `Reduce`-based inversions, `Xi_v(F1) ≈ −1.01675633282526`). **NO checkpoints in the retro-sweep range {121/122/123}** — IV.3 has no checkpoints — so there is no checkpoint-constant provenance to log this sweep. **Cumulative checkpoint-constant provenance is unchanged (105 retained at the higher-bar standard); no checkpoint constant moved.** Previous V.3 entry retained below.)

Snapshot date prior: `2026-06-01` (**batch V.3 close — second batch of the resumed first pass; THE DUAL-ENGINE RULE CORRECTION batch.** First-pass paper-grounded audit on stages 188–200. **Checkpoint stage `200` (reference-free home-stretch theorem) is the only checkpoint in this batch** and returned `verified` at the higher bar with **`material_change: false`**. **V.3 introduced NO new postulated or derived constant/literal in the checkpoint range — 200 is ROUTE-ONLY.** The two 200 findings are both route/structure changes: F1 DE-TRANSLITERATED the `.wl` Section I to a `ratioSubs` helper-monomial-quotient route feeding the `Mderived` Jacobian, and F2 de-tautologized Section III's `Log[a^b]` collapse to the full `ctrMonomial[...]/CtrTarget`. Neither introduces a literal; both REUSE the existing source-anchored Stage 200 constants (`chi_Q=1`; the packet-length/pairing coefficient `2`; the `3`/`5`/`9` carried from `chi_Q = 3(S β^5 + 9 Σ_5)/(3 S − Σ_0)`) already logged in the Stage 200 entry below. **Cumulative checkpoint-constant provenance is unchanged from the IV.2 close (105 retained at the higher-bar standard); no checkpoint constant moved.** The other 12 stages (188–199) are NOT checkpoints (the batch-wide event was the creation of 12 new independent-route `.wl`, all genuine, 0 sanctioned mirrors — a coverage/independence change, not a constant change). Consult `redteam/codex_reviews/_consult_V3.md` (session 019e843e, covering 189 Section II only — not a checkpoint). Previous V.2 entry retained below.)

Snapshot date prior: `2026-05-30` (**batch V.2 close — first batch of the RESUMED FIRST PASS** after the IV.x/V.1 orchestrator-direct integrity remediation (batches 1–8) closed; run under the RESTORED Codex-as-fix-applier contract. First-pass paper-grounded audit on stages 176–187. **Checkpoint stage `185` (microscopic-monomial compiler) is the only checkpoint in this batch** and returned `verified` at the higher bar with **`material_change: false`**. **V.2 introduced NO new postulated or derived constant/literal in the checkpoint range** — the new/strengthened 185 checks REUSE existing symbols (`chi0s`, `deltaUs`, `E_star`, `F_star`, the carried Stage-183 exponents `1+deltaU_*`/`1+chi0_*`, and the source-anchored Stage-183 `2`/`11`/`9`/`4` coefficients already logged in the Stage 185 entry below). The F1 iteration-2 rework reconstructs `Theta₁`/`Xi₁` from the slippage drifts (`chi1_indep=chi0s·Σ_chi`, `deltaU1_indep=deltaUs·Σ_delta`, …) so that the transfer coefficients `C_tr,*`/`A_tr,*` — which are EXISTING paper formulas — become LOAD-BEARING (a wrong typed coefficient now fails via `(C_typed−C_true)·Σ_tr`, Σ_tr ≠ 0), NOT new values; F2's det check `det M_*^(τ,κη,μ)=1+χ0*` likewise introduces no literal. **Cumulative checkpoint-constant provenance is unchanged from the IV.2 close (105 retained at the higher-bar standard); no checkpoint constant moved.** The other 11 stages (176–184, 186, 187) are NOT checkpoints. Consult `redteam/codex_reviews/_consult_V2.md` (session 019e77af + iter-2 019e77e6). Previous batch-8 entry retained below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 8 — FINAL** — re-verification of ONE stage 175 (V.1-range wall_normalized_load_shape), the **LAST of the 29 findings stages**, whose first-pass fix had been applied orchestrator-direct while Codex was bypassed; REMEDIATION-`verified` on 2026-05-29 with **`material_change: false`** (the fix only ADDS a corroborating independent Mathematica Series-route check — a native `dlogSeries[expr_] := Coefficient[Normal[Series[Log[expr], {eps, 0, 1}]], eps]` extractor and one new `[series route]` `expectZero` for the `Sigma_N` differential log-slope). **Stage 175 / batch 8 introduced NO new postulated or derived constant** — the `dlogSeries` fix adds no literal/constant (the existing identity targets and `-kappa` are unchanged and kept symbolic). **NO checkpoints in the batch-8 range {175}** — `is_checkpoint: false` confirmed for 175; V.1 has no checkpoints — so there is no checkpoint-constant provenance to log this batch. **Cumulative checkpoint-constant provenance is unchanged from the IV.2 close (105 retained at the higher-bar standard);** no checkpoint constant moved. With this close all 29 findings stages are now remediated/`verified`. Direction settled by ONE Claude+Codex read-only consult (4/4 CONCUR, none conceptual), recorded at `redteam/codex_reviews/_consult_batch8.md`. Previous batch-7 entry retained below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 7** — re-verification of four stages 148 (IV.5-range representative_positive_families), 150 (IV.5-range full_profile_residual), 157 (IV.6-range core_mouth_coevolution_status), 166 (V.1-range bundle_inversion_four_drifts) whose first-pass fixes had been applied orchestrator-direct while Codex was bypassed; all REMEDIATION-`verified` on 2026-05-29 with **`material_change: false` on 150/157/166 and `material_change: true` on 148 (the defective Mathematica `dT` was corrected to match the already-correct paper/SymPy `dT`, which did not move — ZERO downstream propagation, no stage stale)**. **NO checkpoints in the batch-7 range {148/150/157/166}** — `is_checkpoint: false` confirmed in `redteam/MANIFEST.yaml` for all four; IV.5/IV.6/V.1 have no checkpoints — so there is no checkpoint-constant provenance to log this batch. **Cumulative checkpoint-constant provenance is unchanged from the IV.2 close (105 retained at the higher-bar standard);** no checkpoint constant moved. Provenance item worth recording (non-checkpoint, but tracked here for the no-magic-numbers log): **148** — the rigidity-kernel traction-shift coefficients `A_T = −4.27263956256927` and `B_T = 0.134875005736706` are GENUINE PAPER LITERALS sourced from `paper/appendices/stage_appendix_part04.tex:846,848` (the projected-S-deviation model `A_T(ḡ−g_*)+B_T(S̄−S_*)`, part04:839-840); both engines now anchor `A_T`/`B_T` to these external paper values (no baked SymPy literal in the `.wl`), and the Mathematica side derives its `aT`/`bT` independently via `D[]` autodiff of `Tm[p]=Sqrt[(9/20)(p/(1−sFormula/4))]` along the S(Π) curve before checking against the literals — the same paper literals as the batch-6 stage-147 entry below, now also anchored at 148. (The buggy Mathematica route that DROPPED the S-follows-Π chain term — computing a WRONG `dTU=0.4976…` while SymPy/paper were correct at `0.5087…` — was corrected; the correct `A_T`/`B_T`/`dT` were unchanged.) Directions settled by ONE Claude+Codex read-only consult (5/6 unconditional CONCUR + 1 conditional conceptual-escalate on 157 RESOLVED against escalation), recorded at `redteam/codex_reviews/_consult_batch7.md`. Previous batch-6 entry retained below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 6** — re-verification of four IV.5-range stages 143/144/146/147 (exp-remainder positivity / threshold-Pi_* independence / affine-law residual / rigidity-kernel projection) whose first-pass fixes had been applied orchestrator-direct while Codex was bypassed; all REMEDIATION-`verified` on 2026-05-29 with **`material_change: false` on every stage**. **NO checkpoints in the batch-6 range {143/144/146/147}** — checkpoints in this neighborhood are 096 and 105 only; IV.5 has no checkpoints — so there is no checkpoint-constant provenance to log this batch. **Cumulative checkpoint-constant provenance is unchanged from the IV.2 close (105 retained at the higher-bar standard);** no derived or carried constant moved on any of 143/144/146/147 (the fixes were verification-hardening only). Provenance items worth recording (non-checkpoint, but tracked here for the no-magic-numbers log): (1) **147** — the rigidity-kernel coefficients `A_T = −4.27263956256927` and `B_T = 0.134875005736706` are GENUINE PAPER LITERALS sourced from `paper/appendices/stage_appendix_part04.tex:846,848`; (2) **144** — the threshold root `Pi_*` is anchored to the OWNING stage-131 value `1.50882951349315558300555075595` (not 144's own/truncated value), the correct non-self-anchor; (3) **FLAG (147)** — the `|A_T|/B_T` ratio `31.6785` is a SCRIPT-DERIVED cross-check value, NOT a paper constant; it had been MISLABELED "paper-quoted" and is now correctly relabeled as the script's own computed ratio (paper quotes only `A_T`/`B_T`). Directions settled by ONE Claude+Codex read-only consult (4 CONCUR / 2 DISPUTE-resolved, none conceptual), recorded at `redteam/codex_reviews/_consult_batch6.md`. Previous batch-5 entry retained below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 5** — re-verification of four Family-1 mouth-gain/branch stages 134/137 (IV.4-range canonical_mouth_gain / schur_reduction) and 139/142 (IV.5-range family1_actual_mouth_gains / mouth_gain_susceptibility) whose first-pass fixes had been applied orchestrator-direct while Codex was bypassed; all REMEDIATION-`verified` on 2026-05-29 with **`material_change: false` on every stage**. **NO checkpoints in the batch-5 range {134/137/139/142}** — `is_checkpoint: false` confirmed in `redteam/MANIFEST.yaml` for all four; IV.4 and IV.5 have no checkpoints — so there is no checkpoint-constant provenance to log this batch. **Cumulative checkpoint-constant provenance is unchanged from the IV.2 close (105 retained at the higher-bar standard);** no derived or carried constant moved on any of 134/137/139/142 (the fixes were verification-hardening only — a removed X−X gain-line assert at 134, an independent matrix-Schur reconstruction at 137, a branch-discrimination value anchor + independent `S_q(Π_*)` recompute at 139, and a cross-stage independent-Π_* anchor + independent projection integral at 142; the carried readback `Π_*=1.50882951349316` and the canonical `g_-^F1≈0.758035078944663` are unchanged). Directions settled by ONE Claude+Codex read-only consult (7/7 CONCUR, none conceptual), recorded at `redteam/codex_reviews/_consult_batch5.md`. Previous batch-4 entry retained below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 4** — re-verification of stages 125/126 (IV.3-range positive_source_theorem / positive_source_families) and 130/131 (IV.4-range mouth_bias_map / parent_mouth_threshold) whose first-pass fixes had been applied orchestrator-direct while Codex was bypassed; all REMEDIATION-`verified` on 2026-05-29 with **`material_change: false` on every stage**. **NO checkpoints in the batch-4 range {125/126/130/131}** — IV.3 and IV.4 have no checkpoints — so there is no checkpoint-constant provenance to log this batch. **Cumulative checkpoint-constant provenance is unchanged from the IV.2 close (105 retained at the higher-bar standard);** no derived or carried constant moved on any of 125/126/130/131 (the fixes were verification-hardening only — strengthened positivity/bound checks at 125/126, a global FKG/Chebyshev monotonicity certificate at 130, and dropped-tautology + independent-root-route hardening at 131; the boxed `g_Π` closed form and the readback `Π_*=1.50882951349316` are unchanged). The resolution directions came from a Claude+Codex consult (Codex base `019e7594`, all CONCUR, none conceptual), recorded at `redteam/codex_reviews/_consult_batch4.md`. Previous batch-3 entry retained below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 3** — re-verification of IV.3-range stages 117/118/119/122 whose first-pass IV.3 fixes had been tainted-applied (committed in the IV.3 pass `b4e02d8`) while Codex was bypassed; all REMEDIATION-`verified` on 2026-05-29 with **`material_change: false` on every stage**. **NO checkpoints in the batch-3 range {117–122}** — IV.3 has no checkpoints, and the only nearby checkpoint (105) was batch 2 — so there is no checkpoint-constant provenance to log this batch. **Cumulative checkpoint-constant provenance is unchanged from the IV.2 close (105 retained at the higher-bar standard);** no derived or carried constant moved on any of 117/118/119/122 (e.g. 118's λ-sign resolution to MINUS aligns section V with the script's own already-correct section IV and the notes' boxed forms — no value changed). The 117/118 resolution directions came from a Claude+Codex consult (Codex base `019e74f7`, all CONCUR), recorded at `redteam/codex_reviews/_consult_batch3.md`. Previous batch-2 entry retained below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 2** — re-verification of stages 105/106/109/112 whose first-pass IV.2 fixes had been applied orchestrator-direct while Codex was bypassed. **Checkpoint stage `105` is the only checkpoint in this batch.** 105 was REMEDIATION-`verified` (2026-05-29) with **`material_change: false`** — the `.wl` chi_Q derivation was rewritten along an independent residue/`Reduce`-witness route, but no derived or carried constant changed, so **105's constant provenance is unchanged from the IV.2 close** (`sigma_Q^can = 4 a^5/(27 c_s^5)` still pole-scale-derived in-script; `chi_Q = 1` still derived non-tautologically — now via the residue/`Reduce` witness on the Mathematica side; deformed-branch coefficients still imported with explicit provenance; zero unexplained literals). The IV.2 105 provenance entry below remains accurate. The other three batch-2 stages (106, 109, 112) are NOT checkpoints, so no further checkpoint-constant provenance to log; cumulative unchanged. Previous V.1 entry retained below.)

Snapshot date prior: `2026-05-28` (batch V.1 close — first-pass paper-grounded audit on stages 164-175. **No checkpoints in V.1 range**, so no checkpoint-constant provenance to log; cumulative unchanged. Previous IV.6 entry retained below.)

Snapshot date prior: `2026-05-28` (batch IV.6 close — first-pass paper-grounded
audit on stages 151-163. **No checkpoints in IV.6 range.** Cumulative
checkpoint-constant provenance unchanged from IV.2 close (105 retained at
higher-bar standard). Previous IV.5 entry retained below.

Snapshot date prior: `2026-05-27` (batch IV.5 close — first-pass paper-grounded
audit on stages 139-150. **No checkpoints in IV.5 range.** Cumulative
checkpoint-constant provenance unchanged from IV.2 close (105 retained at
higher-bar standard). Previous IV.4 entry retained below.

Snapshot date prior: `2026-05-27` (batch IV.4 close — first-pass paper-grounded
audit on stages 127-138. **No checkpoints in IV.4 range.** Cumulative
checkpoint-constant provenance unchanged from IV.2 close (105 retained at
higher-bar standard). Previous IV.3 entry retained below.

batch IV.3 close — first-pass paper-grounded
audit on stages 115-126. **No checkpoints in IV.3 range.** Cumulative
checkpoint-constant provenance unchanged from IV.2 close (105 retained at
higher-bar standard). Previous IV.2 entry below.

batch IV.2 close — first-pass paper-grounded
audit on stages 103-114. **Checkpoint stage `105` (chi_Q fix from outgoing
DtN) verified after first-pass cycle at the higher-bar standard.**
Constant provenance assessment for 105: the load-bearing constant
`sigma_Q^can = 4 a^5/(27 c_s^5)` is *derived* in-script from the pole-scale
identity `sigma_Q^can = (9/8)/Omega_Q^5` with `Omega_Q = 3 c_s/(2 a)`,
checked via `expect_zero("sigma_Q^can - 4 a^5/(27 c_s^5)", ...)` (rather than
literal-asserted). The canonical odd-coefficient identification
`chi_Q = 1` is derived non-tautologically via two independent paths: SymPy
solves `sp.Eq(Yret.coeff(omega, 5)/I, a^5/(27 c_s^5))` directly for
`chi_Q`; Mathematica uses `Reduce[c5/I == a^5/(27 c_s^5), chiQ, Reals]` on
an `Apart`-decomposed retarded-module form. The deformed-branch
coefficients `(1, 1/9, 4/81, xi_Q/27)` for `Y_def` are derived in SymPy
via `sp.series(-3/Lambda_def)` and in Mathematica via a structurally
distinct polynomial inversion `Solve[Lambda_def · y_ansatz = -3]`, with
neither engine importing the other's RHS. **Zero unexplained literals in
checkpoint 105's scripts**; every constant is either pole-scale-derived,
imported with explicit provenance from stages 074/088 (now 104) on the
SymPy docstring, or pinned via a paper-quoted carry-in (`Lambda_out`
fingerprint coefficients from stage 104).

Snapshot date prior: `2026-05-27` (batch IV.1 close — first-pass paper-grounded
audit on stages 091-102. **Checkpoint stage `096` (geometry-lane check verdict)
verified after first-pass cycle.** Constant provenance assessment for 096:
the four cardinal constants (`c_pole = 1/4`, `c_geom = 3/4`, `rho_alpha = 4/3`,
`zeta_req = 1/3`) are *derived* in-script from the static-limit hypotheses
`eps_2 = eps_4 = 0` (themselves established by stage 094's orthogonality
output `K_{g,2} = K_{g,4} = 0`) via the `c_pole = (1 + eps_4)/(4(1 + eps_2)^2)`
obstruction formula carried from stage 092. The `Yhat_Q^cons(omega) = 3/4 +
(1/4)/(1 - omega^2/Omega_Q^2)` partial-fraction form is built from these
derived constants, not literal-asserted. The 15-mode orthogonality block at
the top derives `K_{g,2} = K_{g,4} = 0` from explicit angular integrals
(5 Y2A labels × 3 checks each: overlap, Laplace eigenvalue, gradient cross).
The `l = 2` Laplace eigenvalue `6 = ell(ell+1)` is documented in the
constant-provenance docstring. **Zero unexplained literals in checkpoint
096's scripts**; every load-bearing constant is either derived in-script
from orthogonality / obstruction-formula inputs or carry-forward with explicit
source anchor.

Snapshot date prior: `2026-05-27` (batch III.5 close — first-pass paper-grounded
audit on stages 085-090. Checkpoint stages `089` and `090` (both in III.5)
verified after first-pass cycle. **Constant provenance assessment for III.5
checkpoints**: 089's `Pe_suff_chi = 96.5285247264386` and `Pe_fail_chi =
11220.5441626259` are now anchored by an explicit provenance comment
naming `scripts/output/moving_throat_pde_stage082_*_sympy_audit.txt` as the
upstream source (SymPy side; per pitfall #10 SymPy nsolve was not used).
On the Mathematica side, the same Pe values are *rederived* via
`FindRoot[zetaF1[pe] == zetaTarget, {pe, …}]` from notes-quoted
`rho_target - 1`, giving a second-engine independent path. 090's
`c_contact = 3/4` and `c_pole = 1/4` are paper-quoted minimal-isotropic
module coefficients; both engines derive `rho_alpha = 1/c_contact` and
`zeta_req = c_pole/c_contact` from them (no longer hardcoded on the
Mathematica side). Both checkpoints add a `Pe_req = 0` carry-forward
proxy with source-anchored comments. Every load-bearing constant in
III.5 checkpoint scripts is now either derived in-script or
carry-forward with explicit source anchor.

Checkpoint stage `069` falls in III.3's range; returned
**clean (0 findings)** under v2 with no hardcoded constants — every
constant in 069's scripts is either derived (`Cres2` from definition,
`PresGap` via `Solve` in Mathematica) or carried forward with source
anchor (`Cres2 ≈ 0.99441883...` and `Pres ≈ 1.005612487...` carried
from stage 067's sech-Gaussian benchmark with 60-dps mpmath provenance,
`Pe_req` carried from stages 048/049/052). The III.3 v2 sweep surfaced
0 new hardcoded_result findings across the batch. The 2 substantive
paper_misalignment items at stage 062 added: (a) closed-form definition
`C_sp_sq := Osp²/(Nss·Npp)` (no new constant) and (b) Cauchy
parameterization `Osp = cos(θ)·√(Nss·Npp)` (introduces only the
declarative symbol `θ`); the σφ coupling sign flip changed only the sign
of an existing expression. The two banner-relabel items (067, 072) are
pure string changes. III.2 v2 update text remains accurate for that range;
I.1, I.2, II.1, III.1 v2 close text remains accurate for those ranges;
III.4 update remains accurate for batch III.4.)

## Audit Rule

Every constant used in a checkpoint audit should fall into one of three buckets:

- `derived in audit`
- `carried forward with source anchor`
- `probe-only numeric value labeled`

Anything outside those buckets is treated as an audit defect.

## Completed Entries

### Stage 001

- Canonical stage: `paper/stages/stage_001.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `1/(2 sqrt(pi))`
  derived as the normalized real monopole harmonic `Y_00`
- `2 sqrt(pi)`
  derived from the mouth-average extraction rule
  `delta a = q_00 / (2 sqrt(pi))`
- `6`
  derived as the specialization `ell(ell+1)` at `ell = 2`

Audit note:

- The earlier surface-measure ambiguity is now explicit rather than hidden.
  The Stage 001 audits check both the weighted wall-action form and the
  densitized convention actually used by the ledger.
- The current canonical Stage 001 card now states that the ledger adopts the
  densitized one-dimensional convention, rather than leaving that choice
  implicit.
- The added Mathematica mirror rechecks the harmonic normalization, chain-rule
  sign, and weighted-vs-densitized wall-action split in a second CAS.

### Stage 002

- Canonical stage: `paper/stages/stage_002.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage002_breathing_reduction_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `1/(2 sqrt(pi))`
  carried from the normalized `Y_00` convention fixed in Stage 001
- `2 sqrt(pi)`
  carried from the Stage 001 mouth-average bridge
- `4 pi`
  derived in the Stage 002 audit as `(2 sqrt(pi))^2`
- `6`
  derived as the specialization `ell(ell+1)` at `ell = 2`

Audit note:

- The audits check the conservative `(a, L)` matrix reduction and the grouped
  real `P_2` degeneracy without introducing any free numerical coefficients.
- The current Stage 002 note and canonical stage card now state explicitly that
  the Stage 001 surface weight has already been absorbed into the effective
  axial coefficients and profiles before the reduced overlaps are written.
- The added Mathematica mirror rederives the `Y_00` bridge, the `4 pi`
  reduction factor, the conservative two-mode matrix system, and the `ell=2`
  restoring shift independently of the SymPy implementation.

### Stage 003

- Canonical stage: `paper/stages/stage_003.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage003_bdg_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl`
- Numerical stress:
  `scripts/numerical/stage003_021_foundational_stress.py`
  and
  `mathematica/numerical/stage003_021_foundational_stress.wl`
- Verdict: `clean for current symbolic + stress scope`

Constants reviewed:

- none in the theorem path
  the Schur-complement kernels, low-frequency moments, pole formulas, and
  grouped trace / anomaly identities are all checked symbolically in both CAS
  layers
- JSON sample values in
  `scripts/numerical/stage003_021_foundational_samples.json`
  probe-only numeric values labeled

Audit note:

- The current Stage 003 checkpoint is not resting on hidden literals.
- The symbolic theorem path is exact in both CAS layers.
- The shared numerical-stress harness now resolves its repo-local sample JSON
  after the `research/pde_ledger/` move, so the stress layer is runnable rather
  than merely listed in the stage card.
- The sample JSON is used only for perturbative-validity and scaling probes; it
  does not supply any constant used to derive the Stage 003 formulas.
- Red-team batch I.1 (2026-05-21) identified and patched a Wolfram Language
  multi-line continuation defect in the original `.wl`: the `lRed = ...`
  assignment spanning several lines only captured the kinetic terms, missing
  the potential and coupling additions. Downstream results were unaffected
  because the dispersion derivations flowed through `mMat`, `kMat`, `cMat`,
  and `oMat` rather than from `lRed`. The corrected `.wl` now adds the missing
  terms via parenthesised `lRed = lRed + (...)` and verifies all four EL
  residuals against the SymPy convention.

### Stage 022

- Canonical stage: `paper/stages/stage_022.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `27`, `9`, `81`
  derived in the compact outgoing `l=2` fingerprint
  `a^5/(27 c_s^5)`, `a^2/(9 c_s^2)`, and `4 a^4/(81 c_s^4)`
- `54`, `6`, `8`, `15`, `5`
  derived by solving the invariant normalization product against
  `2 G / (5 c^5)` and then collapsing the `K_2`, `K_4` target formulas
- there are no probe-only decimals in the theorem path

Audit note:

- The current SymPy audit explicitly performs the Stage-021 dictionary
  back-substitution round-trip for `N0`, `N2`, and `N4`; the earlier review
  note claiming that gap was stale.
- The Mathematica mirror independently replays the same round-trip and the
  normalization-product solve.
- The checkpoint therefore no longer depends on unverified dictionary
  substitutions or unexplained literals.

### Stage 023

- Canonical stage: `paper/stages/stage_023.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `5`, `20`, `4`
  derived in the grouped `(1,2,2)` metric as the exact `Ggrp` norms of
  `ebar`, `ea`, and `eb`
- `1/5`, `2/5`, `1/10`, `1/2`
  derived from the corresponding rank-one projectors and grouped inverse-map
  formulas under the same weighted metric
- `2`, `3`
  derived in the prefactor-cancellation solve for the exact `N2` and `N4`
  constant-prefactor targets
- `54`, `27`, `5`
  carried forward with source anchor from the Stage-022 normalization bridge;
  Stage `023` now rebuilds `Gamma5_port` through the Stage-021 exact outgoing
  `l=2` branch before using the invariant product
- there are no probe-only decimals in the theorem path

Audit note:

- The earlier “assembled reduced coefficients only” caveat is now stale.
  Both CAS layers explicitly reconstruct representative one-port `Z_n` and
  `N_n` formulas from the underlying `(\Delta, S, Q, H, P)` lane data before
  assembling the grouped packet.
- The Stage-022 outgoing odd coefficient is no longer a bare literal at the
  bundle level. Both CAS layers now rebuild `Gamma5_port` through the same
  Stage-021 exact outgoing route used in Stage `022`.
- The remainder of the checkpoint then stays symbolic: grouped decomposition of
  `D_A0`, `D_A2`, `D_A4`, `N_A0`, isotropic reduction, prefactor constraints,
  anisotropy transport, and monotonicity derivatives.
- The current checkpoint therefore does not depend on hidden literals or on an
  unverified black-box coefficient handoff.

### Stage 024

- Canonical stage: `paper/stages/stage_024.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `sqrt(15/(8 pi))`
  derived as the normalization factor for the real STF `\ell = 2` harmonics
- `4 pi / 15`, `4 pi / 105`
  derived unit-sphere fourth and sixth moments used in the exact overlap
  contractions (sixth-moment prefactor corrected from the stale `4 pi / 122`
  notes typo to the verified `4 pi / 105` in pass-2 batch II.1; see the II.1 entry
  below and PAPER_CLEANUP P5-04)
- `sqrt(5)/(7 sqrt(pi))`
  derived from the exact `Y_20` triple-overlap matrix
- `1`, `1/2`, `-1`
  derived grouped axisymmetric splitting signature from that overlap matrix
- `1/4`, `3/4`
  derived grouped defect weights implied by the `20/21/22` signature
- there are no probe-only decimals in the theorem path

Audit note:

- The earlier notation caveat is stale. The current note explicitly says that
  `H_r = G_{U,r}^2 + G_{W,r}^2` is just the Stage-6 combined gauge/mixed
  strength written with a new letter to avoid collision with Newton's `G`.
- The earlier “tautological Section II” caveat is also stale. Both CAS layers
  now include unequal-lane witness checks showing that grouped defects become
  nonzero off the isotropic locus.
- The stage still treats the radial/axial overlap layer as carried reduced data
  from Stage 6, but that scope limit is explicit in the note and canonical
  stage card. The checkpoint claim itself is the angular closure and splitting
  law, and that claim is now clean.

### Stage 036

- Canonical stage: `paper/stages/stage_036.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `9`, `11`, `8`
  carried forward with source anchor from the Stage-18 selected-branch
  D/N-normal-form coefficients and the `8 / (pi^2 A)` dimensionless scaling
- `18`, `81`
  derived algebraically inside the derivative and endpoint reductions from the
  same closed form `G(xi,delta)`
- `2/9`
  derived in the near-onset expansion coefficient `-2 xi^2 / (9 delta)`
- there are no probe-only decimals in the theorem path

Audit note:

- The earlier boundary-assumption caveat is stale. The current SymPy audit uses
  `xi >= 0` via `nonnegative=True`, and the Mathematica mirror assumes
  `0 <= xi < 1` explicitly.
- Both CAS layers now check the same live theorem packet: exact `G`, exact
  loading split, manifestly positive derivative form, onset value, upper
  endpoint, and near-onset series.
- The checkpoint therefore no longer depends on a hidden `xi > 0` restriction
  that would narrow the onset-boundary claim.

### Stage 089

- Canonical stage: `paper/stages/stage_089.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `4/3`, `1/3`
  carried exact minimal-isotropic loading and support demands from the Stage-71
  precursor packet
- `3.46622291347846`, `3.46752913273870`, `3.46752922945601`
  carried Family-1 loading-ratio window markers from the upstream support-window
  packet
- `2.46752922945601`
  carried Family-1 hard support ceiling from the Stage 63/64 ceiling packet
- `1.00005192880220`
  carried zero-bias Family-1 baseline from the Stage-62 transport map
- there are no probe-only decimals in the theorem path

Audit note:

- This checkpoint is a closed arithmetic theorem conditional on the upstream
  minimal-isotropic module. It does not need to rederive the earlier support
  windows in order to verify the explicit zero-bias Family-1 verdict.
- Both CAS layers keep the carried thresholds explicit and check the exact
  ordering and margin claims without hidden literals.
- The resulting theorem claim is exactly the one stated in the stage card:
  the explicit Family-1 branch succeeds already at `Pe_req = 0` for the
  minimal isotropic passive/outgoing demand.

### Stage 090

- Canonical stage: `paper/stages/stage_090.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage090_updated_reduced_status_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage090_updated_reduced_status_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `3/4`, `1/4`
  carried from the minimal isotropic conservative module fixed upstream in the
  grouped-`P_2` packet
- `4/3`
  derived in the audit as `1 / (3/4)`
- `1/3`
  derived in the audit as `(1/4) / (3/4)`
- `3.46622291347846`
  carried forward as `rho_suff^(chi)` from the Stage 69 support-window audit
- `2.46752922945601`
  carried forward as `zeta_max^(F1)` from the Stage 63/64 support ceiling
- `1.00005192880220`
  carried forward as `A_F1` from the Stage 62 transport map

Audit note:

- This script is intentionally a status-consistency audit. It does not smuggle
  those decimals in as unexplained literals; it declares them as carried
  threshold data and uses them only for the branch-ordering checks that Stage 73
  is supposed to summarize.
- Because the checkpoint claim is itself only that reduced status boundary, and
  the carried inputs are explicit and source-anchored, this narrow audit surface
  is sufficient for citation support.

### Stage 096

- Canonical stage: `paper/stages/stage_096.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.py`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `6`
  derived in the audit as the Laplace-sphere eigenvalue `ell(ell+1)` at
  `ell = 2`
- `1/4`
  derived from the obstruction formula's general free-symbol form, then taken at
  the static limit `eps_2 = eps_4 = 0` (red-team pass-2 batch IV.1, 2026-06-05,
  de-tautologized SECTION II — the formula is now exercised with `eps_2`/`eps_4`
  free + two can-fail off-static probes `1/2` and `1/16`, not only at the static
  point; the value `1/4` did not move)
- `3/4`
  derived as `1 - c_pole`
- `4/3`
  derived as `1 / c_geom`
- `1/3`
  derived as `c_pole / c_geom`

Audit note:

- The script rechecks the isotropic `l=0 <-> l=2` decoupling directly instead
  of treating `eps_2 = eps_4 = 0` as a free status assertion.
- Red-team pass-2 batch IV.1 (2026-06-05) re-verified at the higher checkpoint
  bar with no tier shift: SECTION II's obstruction formula
  `c_pole=(1+eps_4)/(4(1+eps_2)^2)` had been evaluated only at the degenerate
  static point, so its eps-structure was never exercised; the fix adds the
  free-symbol general formula + two can-fail off-static probes + the static limit
  taken FROM the general formula. SECTION I (orthogonality / `ell(ell+1)=6`) was
  already substantive and unchanged. `material_change: false`, no constant moved,
  no new pinned constant.

### Stage 163

- Canonical stage: `paper/stages/stage_163.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage163_off_family_normal_coordinate_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- no free constants in the symbolic theorem path
  the parent compensation defect, `delta_perp` transport, microscopic
  parent-variable formula, outlet-defect transport, and tangent/normal split
  are all checked symbolically in both CAS layers
- `1.77799353547498`, `0.758035078944663`, `4.651033550168876`,
  `0.6703621156734617`
  carried Family-1 readbacks for `r_*`, `g_*`, `Sigma0_can`, and `S_can`,
  used only in the final numerical coefficient banner
- the derived readback coefficients
  `4 sqrt(1+r_*^2)`, `-1/sqrt(1+r_*^2)`, `1/(4 sqrt(1+r_*^2))`,
  `Sigma0_can S_can / sqrt(1+r_*^2)`, and `16 / sqrt(1+r_*^2)`
  are explanatory numeric readbacks, not theorem inputs

Audit note:

- The current symbolic theorem path does not depend on the Family-1 readback
  packet. Those values are only used to print the final numerical transport
  coefficients for the canonical point.
- Because the readbacks are explanatory rather than proof-critical, they are a
  provenance item to track, not a remaining trust defect.
- The checkpoint therefore no longer needs a separate hardening gate before it
  can count as citation support.

### Stage 185

- Canonical stage: `paper/stages/stage_185.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl`
- Numerical stress:
  `scripts/numerical/stage185_187_orbit_stress.py`
  and
  `mathematica/numerical/stage185_187_orbit_stress.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `1 + deltaU_*`, `1 + chi0_*`
  carried symbolic exponents from the Stage 183 tracking compiler
- `E_*`, `F_*`
  carried symbolic coefficients from the Stage 183 nontracking compiler; not
  numerically specialized in either CAS layer
- `2`, `11`, `9`, `4`
  source-anchored coefficients appearing in the carried Stage 183 formulas for
  `E_*` and `F_*`, not unexplained literals
- there are no free decimal literals in the symbolic theorem path

Audit note:

- The Stage 185 audits now reconstruct the primitive microscopic ratios
  `(gamma, c_{etaU}, T_U, K_U, K_eta^{(eff)}, K_W^{(eff)}, lambda_W, mu_W)`
  before rebuilding `chi_0`, `delta_U`, `epsilon_W`,
  `Z_W/Omega_W^2`, and `epsilon_eta`.
- The tracking, nontracking, and dressing monomial laws are then checked both
  from those primitive-ratio compilers and from the direct monomial ratios,
  removing the earlier tautology concern.
- The existing `185--187` numerical stress layer remains secondary; the main
  theorem path is now symbolic in both CAS layers.
- Red-team pass-2 batch V.2 (2026-06-09) re-verified at the higher checkpoint
  bar with no constant moved: the `.wl` was RE-AUTHORED so the load-bearing
  monomial-exponent compilation is now DERIVED via a `monomialExponentVector`
  (substitutes each primitive var → `Exp[logVar]`, takes `Log`/`PowerExpand`,
  reads exponents via `Coefficient`) instead of hand-coded identically to the
  `.py`, and the checkpoint quantity `det ∂(Σ_tr,Σ_nt,Σ_eta)/∂(τ₁,κ_η,μ₁) =
  1+χ₀*` was re-derived (`firstRatioDrift` kept, now acceptable since the
  exponent provenance is independent). `material_change: false`, NO checkpoint
  constant changed, NO new pinned constant introduced — the carried `2`/`11`/`9`/
  `4` Stage-183 coefficients and the symbolic exponents/coefficients are
  preserved; both engines now genuinely independent.

### Stage 200

- Canonical stage: `paper/stages/stage_200.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage200_reference_free_home_stretch_theorem_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `1`
  structural zero-defect target carried by the Packet-A finish line
  `chi_Q = 1`
- `2`
  structural packet-length and pairing coefficient carried by the exact
  two-point orbit transport / cocycle formulas
- `3`, `5`, `9`
  carried from the Stage 200 definition
  `chi_Q = 3 (S beta^5 + 9 Sigma_5) / (3 S - Sigma_0)` and therefore source
  anchored rather than inserted ad hoc

Audit note:

- The Stage 200 audits are symbolic and contain no free decimal literals.
- The added Mathematica mirror checks the same reference-free packet identities,
  mismatch conversions, and linearized compiler in a second CAS.

### Stage 203

- Canonical stage: `paper/stages/stage_203.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `1`
  structural target value carried by the scalar closure condition
  `widehat chi_Q = 1`
- `2`
  structural coefficient carried by the Stage 192 quotient matrix and the
  graph-tangent / repair formulas
- `(1 + deltaU_*) / (1 + chi0_*)`
  carried symbolically from the Stage 192 same-free-quintuple graph formulas
- `E_*`, `F_*`
  carried symbolic coefficients from the monomial compiler; not numerically
  specialized in either audit

Audit note:

- The Stage 203 audits are symbolic and contain no free decimal literals.
- The added Mathematica mirror checks the graph-kernel theorem, graph-error
  packet compiler, inverse compiler, and repair vector in a second CAS.

### Stage 218

- Canonical stage: `paper/stages/stage_218.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `5`
  derived in the audit as the number of primitive free axes
  `{\lambda, c, \gamma, U, W}`
- `30`
  derived in the audit as the total count of nonempty proper support strata
  `Sum[Binomial[5,k], {k,1,4}] = 2^5 - 2`
- `1140`
  carried forward as the Stage 215 support-`<=4` global ledger budget
- `162`
  carried forward as the Stage 217 lifted support-5 per-envelope Bezout bound
  `162 = 3^4 * 2`; arithmetically forced (red-team batch VI.1, 2026-06-02,
  corrected the WRONG `179`/`230` typos that had appeared in the published card,
  the Part-VI appendix, and the 217/218 notes — the SCRIPT always carried `162`,
  so no derived constant moved)
- `750`
  carried forward as the Stage 217 projected-chart fallback per-envelope bound
  `750 = 5 * 5 * 5 * 6`
- `324`, `1500`, `1464`, `2640`
  derived in the audit from the carried Stage 215 and Stage 217 budget packets
  (`324 = 2 * 162`, `1464 = 1140 + 324`, fallback `1500 = 2 * 750`, projected
  total `2640`)

Audit note:

- The Stage 218 audits are symbolic / combinatorial and contain no free decimal
  literals.
- The added Mathematica mirror (RE-AUTHORED from a transliteration in red-team
  batch VI.1, 2026-06-02, to a genuinely independent route — native
  `Subsets`/`SubsetQ`/`Boole`/`Tally` set-combinatorics, `Reduce`/`Resolve`/`ForAll`
  for the splice, and independently-generated regime witnesses whose counts DIFFER
  across engines, 256/192/65+63 vs 192/192/64+64) checks the
  boundary-identification counts, splice theorems, and carried budget arithmetic
  in a second CAS.

### Stage 221

- Canonical stage: `paper/stages/stage_221.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `1/2`
  derived in the audit as the exact maximum of the conservative shape factor
  `r / (1 + r^2)` at `r = 1`
- `eta / (1 + eta^2)`
  derived in the audit as the exact low-loss envelope factor on the domain
  `r >= 1 / eta`
- `1 / (2 Q_* eta)`
  derived in the audit by substituting `gamma_* = omega_* / (2 Q_*)` into the
  low-loss detuning threshold
- `7`, `2`, `3`, `1/4`, `5`, `11`, `40`, `80`, `3/5`
  labeled probe-only numerical values used only for the constructive sign/scale
  sanity slice, not for the symbolic proof path

Audit note:

- The theorem path is symbolic in both CAS layers; there are no free decimal
  literals.
- The Mathematica mirror checks the simple-pole normal form, the carried
  Stage 220 derivative identity, the wall-like specialization, the exact
  line-shape tradeoff laws, and the linear survival-window formulas. Red-team
  batch VII.1 (2026-06-02) RE-AUTHORED this `.wl` from a line-by-line
  transliteration to a genuinely independent route (native
  `D[QPi/DeltaPi,portPi]` derivative, `Residue`, `ComplexExpand`, an uncollapsed
  Breit–Wigner form), de-tautologized the F1 survival round-trips, and brought
  deliverable #9 (the linear survival window) into genuine dual-engine coverage
  (previously print-only/tautological); one sanctioned Codex deviation (F3) used
  the native `D[QPi/DeltaPi,portPi]` instead of the directive's leading-minus
  form, reconciling the Stage-220 identity `∂_Π D_Π = −N` — verified correct. No
  constant introduced or moved (`material_change: false`).

### Stage 239

- Canonical stage: `paper/stages/stage_239.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `I_2`
  derived in the audit as the exact diagonal rigid-mouth packet compiler in the
  physical `(U,V)` chart
- `(0, -V, U - V)`
  derived in the audit from the Stage 236 dependent-plane carrier after the
  substitutions `q_nt = U`, `q_eta = V`
- `(0, 1, 1)`
  carried forward as the Stage 236 equal-drift dressing ray and used explicitly
  as a sourced packet direction rather than as a hidden literal
- `U = 0`, `V = 0`
  derived in the audit as the exact Cartesian orbit-lock conditions in the
  physical logarithmic chart

Audit note:

- The Stage 239 audits are symbolic and contain no probe numerics or free
  decimal literals.
- The added Mathematica mirror checks the diagonal physical chart, the exact
  physical-to-microscopic correction compiler, the support-blindness statements,
  and the Cartesian orbit-lock theorem in a second CAS.
- Red-team batch VII.2 (2026-06-03) re-verified at the higher checkpoint bar with
  no constant change: the pre-existing line-by-line transliteration `.wl` was
  RE-AUTHORED to a genuinely independent route — a forward Jacobian of the boxed
  dependent vector + a native `PseudoInverse` left-inverse + a
  `Reduce`/`Equivalent` orbit-lock (vs the `.py`'s backward hardcoded `SrmDep`).
  `material_change: false`; no derived or carried constant moved (route/coverage
  change only).

### Stage 242

- Canonical stage: `paper/stages/stage_242.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `2 / 11`
  carried forward from the coherent local D/N placement map used in the Stage
  242 physical support variable
- `4 / 3`
  carried forward from the Stage 240 selected support demand
  `Pi_tr = (4 / 3) C_mix`
- `8 / Pi^2`
  carried forward from the Stage 240 mixed-support carrier
  `C_mix = 8 Lambda (1 - epsilon) / Pi^2`
- `2 / 3`
  derived in the audit by substituting the selected support demand into
  `varrho_phys = Pi^2 Pi_tr / (16 Lambda)`
- `1 / (2 + beta^2)` and `beta / (1 + beta + beta^2)`
  derived in the audit by rewriting the Stage 241 thresholds in the realized
  support variable `epsilon`
- `2 deltaU / ((1 + deltaU) (11 + 9 deltaU))`
  derived in the audit by differentiating the carried support variable
  `epsilon = epsilon_W (1 - (2 / 11) deltaU / (1 + deltaU))`
- `3 / 2`, `2 / 3`, `13 / 17`, `1 / 3`, `1 / 5`, `7 / 11`, `1`
  labeled probe-only numerical values used only in the rational sanity sample,
  not in the symbolic proof path

Audit note:

- The Stage 242 theorem path is symbolic in both CAS layers; the explicit
  rational sample point is probe-only.
- The added Mathematica mirror checks the realized support coordinate, the
  threshold rewrites, the support-blind orbit packet, the infinitesimal
  observable compilers, and the exact support/orbit split in a second CAS.
- Red-team batch VII.2 (2026-06-03) re-verified at the higher checkpoint bar with
  no constant change: the pre-existing line-by-line transliteration `.wl` was
  RE-AUTHORED to a genuinely independent route — a `Resolve[ForAll,Reals]`
  strict-inequality certificate + `D[]` on the real closed forms (vs the
  abstract-ζ device) + a `logDrift` total-log-differential (vs `Exp[t·d]`); the
  load-bearing twin-window strict inclusion `C_mix < Pi_tr < 2 C_mix` (with
  `Pi_tr = (4/3) C_mix`) is now tested STRICTLY on both engines.
  `material_change: false`; no derived or carried constant moved (route/coverage
  change only).

### Stage 243

- Canonical stage: `paper/stages/stage_243.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `-Sqrt[2] / 4` and `Sqrt[2] Sqrt[Pi] / 8`
  derived in the audit by evaluating the exact Gaussian leakage and work
  integrals
- `1 / 2`
  used only as the standard quadratic normalization in the non-rigid free
  energy and carried transparently in both CAS scripts
- `chi_lam / k_V`
  derived in the audit as the exact non-rigid ratio `V / U`
- `-a / (4 b)` and `1 - b - a^2 / (8 b)`
  derived in the audit as the interior stationary point and vertex value of the
  compensated source quadratic
- `r^-6`, `Exp[-2 kappa r] / r^4`, `Exp[-4 kappa r] / r^2`
  carried forward from the one-port short-range same-charge kernel verdict and
  checked explicitly for vanishing `r * V(r)` tails

Audit note:

- The Stage 243 audits are symbolic in both CAS layers and contain no probe
  numerics.
- The added Mathematica mirror checks the exact Gaussian leakage/work channel,
  the linear `(U,V)` solve, the recovery slice, and the short-range limit
  firewall in a second CAS.
- Red-team batch VIII.1 (2026-06-03) re-verified at the higher checkpoint bar with
  no constant change: the pre-existing line-by-line transliteration `.wl` was
  RE-AUTHORED to a genuinely independent route — IBP closure / native `LinearSolve`
  / `TrigExpand` / `Series`-at-∞ for the leakage-work lane and the non-rigid solve,
  with the hardcoded `expected*` residuals DELETED. `material_change: false`; no
  derived or carried constant moved (route/coverage change only).
- Pass-2 batch VIII.1 (2026-06-09) re-verified at the higher bar: the pass-1 `.wl`
  re-author was confirmed SUFFICIENT by the orchestrator ground-truth `.wl`-vs-`.py`
  read (the VI.1-218 / VII.1-221 / VII.2-239&242 outcome, NOT V.3-200) — the `.wl`
  derives U/V via BOTH `Solve` AND `LinearSolve` (asserts agreement) + a `Series`
  asymptotic route, vs the `.py`'s `sp.integrate`/`sp.solve`/`sp.limit`. No checkpoint
  constant changed, moved, or re-pinned; `material_change: false`.

### Stage 248

- Canonical stage: `paper/stages/stage_248.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`
- Numerical stress:
  `scripts/numerical/stage248_event_chain_stress.py`
- Numerical stress (Mathematica):
  `mathematica/numerical/stage248_event_chain_stress.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `1 / 2`
  carried transparently as the standard kinetic-energy normalization in the
  reduced one-dimensional event chain
- `1 / E`
  derived in the audit as the exact pure-Coulomb outer turning point
- `Pi`
  derived in the audit through the exact near-top parabolic action integral
- `2.5413906350657705`, `3.272783388968954`, `2.1447620194324593`, `0.4`,
  `0.673752615`, `0.546377065`, `1.23312756`, `23.3128`
  derived benchmark outputs checked against the declared Session-II benchmark
  specialization
- `5.0`, `0.18`, `2.5`, `0.19999794`, `3.42933112`, `0.23944389`,
  `0.39096144`, `0.19039548`, `0.19744614`, `0.30222297`, `0.34437471`,
  `0.42826825`, `2.59221845`, `0.28091705`
  labeled benchmark-only numeric inputs copied from the declared Session-II
  readback, not used as part of the symbolic theorem path
- In `stage248_event_chain_samples.json`, the probe families are parameterized
  by barrier gap, subbarrier fraction, cross-window fraction, `I_new /
  I_Coul`, `\Xi_{\rm turn}`, and `|V'(r_{\rm turn})|`.
  Those are probe-only stress inputs.
- The stress harness then derives `V_{\rm peak}`, `E_{\rm sub}`,
  `v_{\rm crit,new}`, `v_{0,\rm sub}`, `v_{\rm cross}`,
  `r_{\rm turn,Coul}`, `\lambda_{\rm th}`, and the transmission ratio from the
  exact Stage 248 formulas instead of hard-coding those outputs.
- The near-top stress block uses probe-only `(\Delta E, K_{\rm peak}, m_s,
  \hbar_{\rm eff})` tuples and compares direct quadrature against the exact
  parabolic top-action formula.

Audit note:

- The Stage 248 theorem path is symbolic in both CAS layers up through energy
  conservation, threshold-speed compilation, Coulomb reference formulas, and
  near-top action.
- The Session-II numbers are intentionally confined to a benchmark-only
  specialization layer; they are tracked as declared readback values rather
  than hidden derivation inputs.
- The new shared numerical-stress layer now checks three admissible event-chain
  probe families plus three near-top probe families in both CAS layers.
- The Coulomb WKB action is numerically integrated and compared against the
  exact closed form in both CAS layers; the near-top action is also compared
  against direct quadrature in both CAS layers.
- Red-team batch VIII.1 (2026-06-03) re-verified at the higher checkpoint bar with
  no constant change: the pre-existing line-by-line transliteration `.wl` was
  RE-AUTHORED to a genuinely independent route — the §II `Solve`-mirror replaced by
  a native SATISFACTION route (the compiler closed forms are verified to satisfy
  their defining energy equalities via substitution + `FullSimplify`, plus a
  non-vacuity guard and a positive-branch guard); this was the batch's one iter-2
  (iter-1's `Reduce`/`ToRules` route was a Wolfram-version dead end Codex correctly
  BLOCKED, orchestrator-reframed to the satisfaction route — a Claude math-coverage
  resolution, non-conceptual, no paper edit). A NOTES-ONLY benchmark typo (notes:506
  `×168%→×100%`, recurring the stale "168") was corrected to the already-correct
  script. `material_change: false`; no derived or carried constant moved.
- Pass-2 batch VIII.1 (2026-06-09) re-verified at the higher bar: the pass-1 `.wl`
  re-author (the SATISFACTION route — closed forms verified to satisfy the defining
  energy equality + non-vacuity guard `FreeQ[deltaNew,v0]`→fail, vs the `.py`'s
  `sp.solve`) was confirmed SUFFICIENT by the orchestrator ground-truth read. No
  checkpoint constant changed, moved, or re-pinned; `material_change: false`.

### Stage 253

- Canonical stage: `paper/stages/stage_253.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.wl`
- Numerical stress:
  `scripts/numerical/stage253_material_threshold_stress.py`
- Numerical stress (Mathematica):
  `mathematica/numerical/stage253_material_threshold_stress.wl`
- The lattice-turnover compiler, legacy recovery slice, harmonic
  trigger/stiffness compiler, Korringa ceiling, and screening ratios are all
  derived symbolically in both CAS layers.
- `\Upsilon_{\rm lat}` remains an explicit calibration parameter; it is not
  solved for or hidden inside a hard-coded literal.
- `4.79562976` is a carried forward legacy Session-V reduced lattice-rate input
  used only in the declared legacy-slice recovery.
- `0.5489386551062235`, `6.94311167`, `0.75`, and `1.0` are benchmark-only
  Stage 252 slice inputs for `s_c`, `s_0`, `f_{\rm lat}`, and `\mu_\eta`.
- `0.39096144`, `0.42826825`, and `2.73855812` are benchmark-only
  turning-point / force-match inputs carried from the declared Session-V
  benchmark slice.
- `65.45193925961132`, `13.64824695299483`, `119.23361317476524`,
  `8.736185210116078`, `0.9128891530016525`, `2.1908464937104797`,
  `10.95423248`, `5.837462857946154`, and `0.5489386551062235` are
  benchmark-derived outputs computed from those declared inputs.
- The Stage 253 theorem path is symbolic in both CAS layers. All explicit
  decimals are benchmark-only readbacks, not hidden theorem inputs.
- In `stage253_material_threshold_samples.json`, the declared
  `Pi_ep/Pi_chi/Pi_k/Pi_t` targets are probe-only screening margins. The
  stress harnesses derive `\lambda_{\rm ep}\omega_D`, `r_{\rm turn}`,
  `k_{\rm eff}`, and `T` from the exact Stage 253 thresholds rather than
  hard-coding those candidate values directly.
- `\Upsilon_{\rm lat}` in the stress layer is either the raw microscopic slice
  `1`, the exact legacy replay `\gamma_{\rm lat,safe}^{\rm eq} /
  \gamma_{\rm lattice}^{\rm red}`, or an explicit custom calibration probe.
- Red-team batch VIII.1 (2026-06-03) re-verified at the higher checkpoint bar with
  no constant change: the pre-existing line-by-line transliteration `.wl` was
  RE-AUTHORED to a genuinely independent route — native `D[Log[V[r]]]` + 5 `Solve`
  energy/force balances + regrouped threshold / Korringa / screening blocks (vs the
  `.py`'s back-substituted `r_turn_phys` expr−expr round-trip). Two NOTES-ONLY
  benchmark typos were corrected to the already-correct script/card (notes:274
  `187.23361317→119.23361317` — corroborated by the notes' own
  65.45193926/0.5489386551=119.2336 and the cross-engine `.wl`; notes:419 a_int
  `10.95423247→10.95423248` = 4·K_turn). The audit agents initially MISattributed
  these to the published card but the orchestrator confirmed the card (140 lines)
  is clean — the typos were NOTES-ONLY. `material_change: false`; no derived or
  carried constant moved. **Stale-value catch (corrected this close):** the
  benchmark-derived-outputs list above carried a STALE `136.23361317476524`
  (pre-existing since the VII.2 close `5860b3a`) for the micro-product threshold;
  the current verified script and `.wl` both compute `119.23361317476524` (= the
  same value to all 14 fractional digits, only the integer part was corrupted),
  and there is NO `136.2336` anywhere in the stage-253 script/`.wl`/output. The
  list was corrected `136.23361317476524 → 119.23361317476524` to match the
  verified engines (and to be self-consistent with the notes:274 correction in
  this same entry). It is an end-of-ladder readback that feeds no downstream
  algebra, so the correction has no propagation.
- Pass-2 batch VIII.1 (2026-06-09) re-verified at the higher bar — 253 is THE FINAL
  STAGE, closing the full second pass 253/253: the pass-1 `.wl` re-author (the
  lattice-balance `Solve` + `chi_lambda` via `D[Log[V[r]],r]=2/r`, vs the `.py`'s
  division/hardcoded `2/r`) was confirmed SUFFICIENT by the orchestrator ground-truth
  read. The benchmark stays the engines' `119.23361317476524` (stale
  `136.23361317476524`/`187.23361317` confirmed fully purged). No checkpoint constant
  changed, moved, or re-pinned; `material_change: false`.

## Open Follow-up

- Symbolic parity is now closed across the checkpoint support set.
- The stage-level trust baseline now lives in `CHECKPOINT_TRUST_AUDIT.md`.
- Dedicated numerical hardening is now in place for checkpoint stages `248`
  and `253`.
- The next checkpoint pass should decide whether to widen numerical stress
  beyond the current four stress-backed checkpoints (`003`, `185`, `248`,
  `253`) and how to expose that verification surface compactly in the paper.

## Pass 2 — Batch I.1 (2026-06-04)

Pass-2 Batch I.1: no carried paper constant changed. (Stage-009 μ₁=2 is the
Gamma-kernel verification first moment, not a carried paper constant.) No new
provenance items.

## Pass 2 — Batch I.2 (2026-06-04)

Pass-2 Batch I.2: no carried paper constant changed; no new provenance items.
Checkpoints 022 and 023 are in range and were re-verified clean with the
value-reconciliation augmentation (all emitted deliverable values reconcile —
0 misaligned); the 022/023 constant-provenance entries above stand unchanged.
The one finding (021, a NON-checkpoint stage) was a label-only `.wl:195` edit
(`Stage 004`→`Stage 021`) + output refresh — no constant or literal introduced
or moved.

## Pass 2 — Batch II.1 (2026-06-04)

Pass-2 Batch II.1: no carried checkpoint constant changed; no new provenance
items. Checkpoints 024 and 036 are in range and were re-verified with the
value-reconciliation augmentation (177 deliverable values checked batch-wide,
1 misaligned — 024's NOTES-ONLY sixth-moment prefactor typo, below — 0
MISSING-DELIVERABLE); the 024/036 constant-provenance entries above stand
unchanged (κ_*=√5/(7√π) and the M^(20) overlap matrix unchanged). The one
value correction was 024's `paper_misalignment`: the notes' (line 213)
sixth-moment prefactor `4π/122`→`4π/105` — a NOTES-ONLY value typo (textbook
unit-sphere sixth moment 1/(3·5·7)=1/105 in d=3; BOTH engines already use 105,
and the notes' own κ_*=√5/(7√π) requires 105), the published `.tex` card was
UNAFFECTED, and no derived or carried checkpoint constant MOVED (the 024
"Constants reviewed" list above has been corrected from the stale `4 pi / 122`
to the verified `4 pi / 105`; the in-script value was always `4 pi / 105`). The
032-036 fixes were script-label-only `stale_output` (no constant). NO checkpoint
constant introduced or moved.

## Pass 2 — Batch III.1 (2026-06-05)

Pass-2 Batch III.1 (037-048): no checkpoints in range and no new checkpoint
constants. (143 deliverable values checked batch-wide, 0 misaligned, 0
MISSING-DELIVERABLE; `material_change: false` on all 12.) One PROVENANCE
STRENGTHENING worth recording though no value moved: stage 043's `M_supp`
baseline `B=κ0²=8/π²` is now DERIVED in-script from the stage's frozen overlap
constant — `B=κ0²=(9/11)·σ` with `σ=88/(9π²) ⇒ 8/π²`, asserted via
`expect_zero("baseline B = 8/pi^2 from frozen sigma", …)` in both engines —
replacing the prior check that re-substituted the SAME literal `8/π²` into both
sides of an already-`B`-symbolic identity (so it could never fail and never
derived the number). The value `8/π²` is unchanged and matches notes:149+151
exactly; the check now actually exercises it (fails if κ0² were mis-stated). No
other provenance items; no carried checkpoint constant changed.

## Pass 2 — Batch III.2 (2026-06-05)

Pass-2 Batch III.2 (049-060): reviewed — NO checkpoint constant changed. ONE
checkpoint in range (051), re-verified clean with its constants unchanged; no new
checkpoint constant introduced. (126 deliverable values checked batch-wide, 0
misaligned, 0 MISSING-DELIVERABLE; `material_change: false` on all 12.) One
PROVENANCE STRENGTHENING worth recording though no value moved: stage 060's
`Xi_micro` baseline (non-checkpoint) is now cross-checked via two INDEPENDENT
routes — the susceptibility route `chiSigma→1/theta` (`wl:133`) vs the
Einstein/diffusion route `dSigma→mSigma*theta` (`wl:136`, `mSigma` cancels),
asserted via `expectZero["Xi_micro chi-route equals D/M-route", xiMicroFromChi -
xiMicroFromDM]` — replacing the prior `wl:140` check that subtracted `xiMicro`
from a verbatim copy of its own `wl:132` definition (identically 0, could never
fail). The value `Xi_micro = Λφ²L²/(Θ T_X)` is unchanged; the check now actually
exercises it. No other provenance items; no carried checkpoint constant changed.

## Pass 2 — Batch III.3 (2026-06-05)

Pass-2 Batch III.3 (061-072): reviewed — NO checkpoint constant changed. ONE
checkpoint in range (069, `final_reduced_verdict`), re-verified clean at the
higher bar with NO constant pinned: `C_res² = 0.994418…` / `P_res = 1.005612…`
are upstream Stages 067/068 deliverables carried as the FREE symbols `Cres2Prim`
(Mathematica) / `Pres_gap` (SymPy) — never re-asserted against themselves — and
the substantive verdict (three-zone classification, matched-window width, exact
side-band widths, `P_res-1=(1-C²)/C²`, strict interior-point ordering) is a set
of genuine symbolic identities verified by two independent engine routes. No new
checkpoint constant introduced. (116 deliverable values checked batch-wide, 0
misaligned, 0 MISSING-DELIVERABLE; `material_change: false` on all 12.) One
PROVENANCE STRENGTHENING worth recording though no value moved: stage 070's
internal sech-profile moments (non-checkpoint anchor) are now ASSERTED rather than
print-only — `I_f = 2/3` (SymPy + Mathematica) and `I_g = 14/15` (Mathematica,
against the 30-digit NIntegrate value ≈0.9333) — replacing the prior `I_1/J_1 =
4πa²ℓ` ratio check whose `I_f/H_w` factor cancelled (independent of the moment
value, could never fail), plus a corrected `wl:86` print annotation `8/15`→`14/15`
(the orchestrator false-positive guard caught that the audit's proposed `8/15`
was WRONG; correct = 14/15). The moments are an internal anchor (not a deliverable,
not consumed downstream); the symbolic `kappa`/`W_wall`/`Xi` deliverables are
unchanged. No other provenance items; no carried checkpoint constant changed.

## Pass 2 — Batch III.4 (2026-06-05)

Pass-2 Batch III.4 (073-084): reviewed — NO checkpoint constant changed, NO
checkpoint in range. III.4 is the Family-1 geometry / thresholds / quadrupole band
(no `is_checkpoint` stage). 119 deliverable values checked batch-wide, 0 misaligned,
0 MISSING-DELIVERABLE; `material_change: false` on all 12; no new postulated constant
introduced. The script-side math fixes (075, 083) DE-TAUTOLOGIZED checks and added
numeric anchors against fixed EXTERNAL paper/notes literals — they do not pin any NEW
constant; the anchored values were already present and correct. Worth recording (no
value moved): the Family-1 healing-window endpoints `Delta_0 = 1.73302079021525e-4`
and `Delta_inf = 2.01447565540522e-2` are now ANCHORED on the SymPy side in BOTH
stage 075 (threshold-window route) and stage 083 (direct-operator route) — two
independent derivations of the same window, mutually cross-corroborating the literals
(previously the SymPy side of each only PRINTED them; the Mathematica side of 075
already pinned them). All carried Family-1 constants (`kappa_F1 = 12321/5`,
`eta_F1 = 37`, `alpha_r² = 100`) are unchanged and still locked. No carried checkpoint
constant changed.

## Pass 2 — Batch III.5 (2026-06-05)

Pass-2 Batch III.5 (085-090): reviewed — NO checkpoint constant changed, and the TWO
checkpoints in range (**089** and **090**) pin **NO new numeric constant**. 59
deliverable values checked batch-wide, 0 misaligned, 0 MISSING-DELIVERABLE;
`material_change: false` on all 6; no new postulated constant introduced. **089** —
`A_F1` (=1.000051928802195328659334) and the `rho_suff/rho_fail/rho_max` thresholds are
DERIVED/carried with provenance, NOT pinned: the Mathematica side independently re-derives
the upstream `Pe` via `FindRoot` (a robust route — NOT the latent `nsolve`-near-`tan`
pitfall #10), and the boxed `Pe_req = 0` is FORCED by the positive zero-bias success
margin `zeta_F1(0) − zeta_min` (= A_F1 − 1/3 ≈ 0.6667185954688619953260008), not asserted
against itself (the de-tautologized F2 replaced the old `0==0` self-check with this
can-fail positivity assert). The de-tautologized F1 made `eps_blk` symbolic and asserts the
`eps→0` reduction of the general `Q_gen` form, so `Q` is now structurally guarded rather
than baked — no value moved (downstream `rho_*` byte-identical). **090** carries the
Family-1 threshold decimals as labeled CARRY-FORWARDS (source-anchored), matching 089's
in-script derivation; its own substantive assertion `zeta_req = rho_alpha − 1` plus the
branch-ordering inequalities introduce no new constant. No other provenance items; no
carried checkpoint constant changed.

## Pass 2 — Batch IV.1 (2026-06-05)

Pass-2 Batch IV.1 (091-102): reviewed — NO checkpoint constant changed, and the ONE
checkpoint in range (**096**, `geometry_lane_check_verdict`) pins **NO new numeric
constant**. 99 deliverable values checked batch-wide, 0 misaligned; `material_change: false`
on all 12; no new postulated constant introduced. **096** cleared the higher bar after
de-tautologizing its SECTION II: the obstruction formula `c_pole=(1+eps_4)/(4(1+eps_2)^2)`
(carried from Stage 092) had been evaluated only at the degenerate static point
`eps_2=eps_4=0`, collapsing to literal `1/4`, so its eps-structure was never exercised. The
fix makes `eps_2`/`eps_4` free symbols + adds the general formula + two CAN-FAIL off-static
probes (`eps_4=1,eps_2=0 → 1/2`; `eps_2=1,eps_4=0 → 1/16`) + takes the static limit FROM the
general formula. **The obstruction-formula deliverable values are PRESERVED, NOT moved or
re-pinned** — `c_pole = 1/4` (now derived from the general formula at the static limit, not
hardcoded at the static point), `c_geom = 3/4`, `rho_alpha = 4/3`, `zeta_req = 1/3`, and the
`Yhat_Q^cons` closed form remain the same source-anchored constants logged in the Stage 096
entry above; they are now exercised via off-static probes rather than only the static point.
SECTION I (l=0⊥l=2 orthogonality, Laplace eigenvalue `6 = ell(ell+1)` at `ell=2`) was
already substantive and unchanged. **Zero unexplained literals remain in checkpoint 096's
scripts.** (For completeness: the non-checkpoint provenance-adjacent fixes this batch
introduced no new constant either — 093 [status-only] and 094 got the same obstruction-formula
de-taut [094's `K_g2`/`K_g4` now DERIVED from the proven-zero l=0↔l=2 overlap moments rather
than bare literal 0]; 100's `.wl` de-transliteration [orchestrator override, user-authorized
independent geometric-series route] moved no emitted value — committed Mathematica output
byte-identical to HEAD, `chiQ` stays a free symbol.) No carried checkpoint constant changed.

## Pass 2 — Batch IV.2 (2026-06-06)

Pass-2 Batch IV.2 (103-114): reviewed — NO checkpoint constant changed, and neither of the TWO
checkpoints in range pins a new numeric constant. 86 deliverable values checked batch-wide, 0
misaligned; `material_change: false` on all 12; no new postulated constant introduced. **105**
(`chi_Q_from_outgoing_DtN`) owns **`chi_Q = 1`**; **112** (`hybrid_robin_branch`) owns
**`gamma_W = 1/9`** (chi_B=1 ⟺ sigma_W(1−9·gamma_W)=0).

**105** cleared the higher bar after DE-TAUTOLOGIZING its retarded-half DtN match. Previously
BOTH engines matched the canonical ω⁵ coefficient to a HARDCODED fingerprint target
`a^5/(27 c_s^5)` (= `sigma_can/4` retyped), never evaluating the actual outgoing l=2 Hankel DtN
fingerprint. Now both engines DERIVE the fingerprint from the spherical Hankel function via
visibly distinct constructions (SymPy `j_2 + i*y_2` closed form, Mathematica native
`SphericalHankelH1[2,z]`), assert the can-fail series `-3 + z²/3 + z⁴/9 + i z⁵/9`, read the
DERIVED imag z⁵ coefficient `1/27`, and force `chi_Q=1` by matching the retarded ω⁵ coefficient
to that DERIVED target — **`chi_Q = 1` is now DERIVED from the Hankel fingerprint, NOT typed on
the RHS; the value is UNCHANGED.** The in-script `sigma_Q^can = 4 a^5/(27 c_s^5)` remains
pole-scale-derived; zero unexplained literals remain in checkpoint 105's scripts. (Provenance:
the 2026-05-29 integrity-remediation batch-2 re-authored 105's `.wl` to a residue/`Reduce`
path, but the hardcoded-RHS verification-strength tautology survived in BOTH engines until
pass-2 caught it.)

**112** cleared the higher bar with NO re-author — its pre-existing `.wl` is genuinely
independent (native `Reduce[presCond==0 && sigma!=0, gamma]` reconstructing the Stage-92
`(b,a0,a5)` from branch-B coefficients to force `gamma_W = 1/9`, a route the `.py` lacks); the
constant `gamma_W = 1/9` is unchanged. (For completeness: the non-checkpoint fixes this batch
introduced no new constant — 107/110/114's `.wl` were RE-AUTHORED from `Series`/`Inverse[m]`
ports to undetermined-coefficient/Schur `Solve` routes, USER-AUTHORIZED, all MMA outputs
byte-identical, no emitted value moved.) No carried checkpoint constant changed.

## Pass 2 — Batch IV.3 (2026-06-06)

Pass-2 Batch IV.3 (115-126): reviewed — **NO checkpoints in range**, NO checkpoint constant moved.
108 deliverable values checked batch-wide, 1 misaligned (resolved — see below); `material_change:
false` on all 12; no new postulated constant introduced. **Cumulative checkpoint-constant provenance
is unchanged from the IV.2 close (105 retained at the higher-bar standard); no checkpoint constant
moved.** Provenance worth recording (NON-checkpoint, tracked here for the no-magic-numbers log): the
canonical Family-1 geometric radius `r_F1 = √(4107−100π²)/(10π) ≈ 1.77799353547498` — the PUBLISHED
PAPER (`paper/appendices/stage_appendix_part04.tex:562` + `paper/parts/part04_geometry_retarded_mouth.tex:576`)
was CORRECTED to this script-derived form this batch (it had carried a `117π²` typo). `100π²` is
arithmetically forced (`r_geom=√((12/π²)(L/a)²−1)`, `L/a=37/20` ⟹ `12·(37/20)²=4107/100`), is
corroborated by the paper's OWN adjacent numeric `≈1.77799353547498` (the `117π²` form gives ≈1.7295)
and downstream `g_-^{F1}≈0.758035`, and matches EVERY script and note (stages 121/122/123/126/127/142/143/148).
**The SCRIPTS were always correct (`100π²`); the WRONG `117π²` typo lived only in the published paper —
a paper/notes typo correction aligning the prose to the already-correct script (Codex-applied
out-of-band, Claude-reviewed); no derived or carried constant MOVED.** Flagged by 122 F1 / 123 F1 as
`paper_misalignment` (value_mismatch); pass-2 caught a published-paper arithmetic typo the first pass
MISSED. **Re-confirmed (ANSATZ catalog, NON-checkpoint):** `γ₀ = (1+r_c)/9` is a POSTULATED pure-scale
ANSATZ, NOT derived (117 §5 / 116 "Bare outgoing normalization" note) — already the ansatz catalog's
first entry. (For completeness: the two script-side fixes introduced no new constant — 117 F1
`mathematica_transliteration` FULL re-author to an independent route [§1–§4 undetermined-coefficient
solve + §5 2×2 `coreMatrix` Solve-elimination; committed `.wl` byte-identical, no emitted value moved],
118 F1 `tautological_check` de-taut [`K_q`/`kQ` tied to the independently-computed gradient integral
`chiGrad = ∫(χ')²dz`, asserted `=π²/(4 L_W²)`; printed `K_q` value unchanged]. 120 & 124 are status-only
by design.)

## Pass 2 — Batch IV.4 (2026-06-06)

Pass-2 Batch IV.4 (127-138): reviewed — **NO checkpoints in range**, NO checkpoint constant moved, no
new checkpoint constant introduced. 1 value misaligned (resolved — 127's notes digit below), all other
emitted deliverable values reconciled; `material_change: false` on all 12. **Cumulative
checkpoint-constant provenance is unchanged from the IV.2 close (105 retained at the higher-bar
standard); no checkpoint constant moved.** The pass-2 re-confirmed carry-forwards
`Π_*=1.50882951349316`, `S_q(Π_*)≈0.658075937605428`, and `Σ_m*≈0.451485277739090` are owned UPSTREAM
(at 130/131/133/135) and are **NOT checkpoint constants** — they are carried, not pinned, in the IV.4
range. The ONE value-mismatch is NOTES-ONLY and NON-checkpoint: 127 x*_exp `0.662765402623160`→
`0.662765402623161` (15th sig-fig transcription typo; both engines compute …161; user-resolved
correct-to-script; published cards UNAFFECTED — the SCRIPTS were always correct). (For completeness:
the three script-side fixes introduced no new constant — 129's `.wl` RE-AUTHORED from a transliterated
profile port to an independent `DSolve` zero-flux-ODE route + Onsager-current assertion, 134's `.wl`
RE-AUTHORED from a postulated-kernel port to an independent `DSolveValue` mixed-D/N BVP route, 135's
X−X closure residual removed; all committed outputs additive/byte-identical, no emitted value moved.)

## Pass 2 — Batch IV.5 (2026-06-08)

Pass-2 Batch IV.5 (139-150): reviewed — **NO checkpoints in range**, NO checkpoint constant moved, no
new checkpoint constant introduced. 0 values misaligned batch-wide (no `paper_misalignment` anywhere →
ZERO paper/notes edits); `material_change: false` on all 12; no new postulated constant. **Cumulative
checkpoint-constant provenance is unchanged from the IV.2 close (105 retained at the higher-bar
standard); no checkpoint constant moved.** The pass-2 re-confirmed carry-forwards
`Π_*=1.50882951349316`, `S_q(Π_*)≈0.658075937605429` are owned UPSTREAM (at 131/133/134/135) and are
**NOT checkpoint constants** — they are carried, not pinned, in the IV.5 range. No stale `168π²`/`168%`
anywhere (only interior mantissa digits); the canonical Family-1 radius `√(4107−100π²)/(10π)` is used
throughout. (For completeness: the two script-side findings introduced no new constant — 139's `.wl`
RE-AUTHORED from a NUMERICAL mirror [imported `Π_*`/`S_q` literals + a hardcoded Stage-134 kernel closed
form, "sanctioned mirror" comment] to an independent route [`FindRoot` Π_* on `g(Π)=g_minus`, correct
sign `(e^p−1)`; `S_q(Π_*)` via independent source-moment quadrature `∫Σ·K_q` cross-checked against a
symbolic Integrate route; the two literals survive ONLY as `expectApprox` cross-check targets], `.py`
UNCHANGED, deliverables preserved at 1e-12 + gained precision; 146's ε-slope checks de-taut'd in BOTH
engines, no emitted value moved.)
No carried checkpoint constant changed.

## Pass 2 — Batch IV.6 (2026-06-08)

Pass-2 Batch IV.6 (151-163): reviewed — NO checkpoint constant changed, and the ONE checkpoint in range
(**163**, `off_family_normal_coordinate`) pins **NO new numeric constant**. 0 values misaligned
batch-wide (no `paper_misalignment` anywhere → ZERO paper/notes edits); `material_change: false` on all
13; no new postulated constant introduced. **163** cleared the higher bar with BOTH engines confirmed
genuinely independent (the `.wl` adds an implicit-function-theorem slope `−F_r/F_g` + a full `Series`
perturbation route the `.py` lacks → strictly stronger). Its load-bearing constant `4√(1+r*²)` (printed
`8.15966765224253`) is **re-derived from `∂F/∂g` — NOT trusted as a literal, NOT moved or re-pinned**;
the value is UNCHANGED. **Zero unexplained literals remain in checkpoint 163's scripts.** **Cumulative
checkpoint-constant provenance is unchanged from the IV.2 close (105 retained at the higher-bar
standard); no checkpoint constant moved.** (For completeness: the three script-side findings introduced
no new constant — 158's `.wl` RE-AUTHORED [⭐ ORCHESTRATOR OVERTURN of a CLEAN verdict, USER-AUTHORIZED]
from a `Series`-on-shifted-form mirror to independent analytic base-point differentiation [`D[rFun,g]/.
g→gStar`, `D[chi,eps]/.eps→0`, `dRFromDg` from `D[(g−r)²/(1+r²),g]/.g→gStar`]; all `expectZero`
identities/targets unchanged, all 8 printed coefficients byte-identical, committed `.wl` output
byte-identical to HEAD, `.py` unchanged; 161's `.wl` RE-AUTHORED [USER-AUTHORIZED route-only] from a
near-line-by-line port to a Series+Coefficient/logarithmic-derivative route [`D[Log[kappaRatio],…]`], the
prefactor `(1+r_F1²)/9=0.46236233468786880105…` and the `Δ_Q` coeffs `5.352238871696225`/
`10.70447774339245`/`−1.16275838754222` preserved, form-only output change, `.py` unchanged; 162's `.py`
stale-numbering COMMENT fixed [`Stages 99 and 102`→`Stage 119 (with gamma0 from Stages 115-116)`;
non-load-bearing, output byte-identical]. 153 is status-only by design. **0 sanctioned mirrors remain in
IV.6.**)
No carried checkpoint constant changed.

## Pass 2 — Batch V.1 (2026-06-08)

Pass-2 Batch V.1 (164-175): reviewed — **NO checkpoints in range; no checkpoint constant pinned, moved,
or re-pinned this batch; no new postulated constant introduced.** 0 values misaligned batch-wide (no
`paper_misalignment` anywhere → ZERO paper/notes edits); `material_change: false` on all 12. **This batch
CLOSES the 105–175 first-pass orchestrator-direct TRANSLITERATION WATCH at 175.** (For completeness, no
new constant entered the provenance ledger: the orchestrator ground-truth `.wl`-vs-`.py` read split the
12 into 6 confirmed genuinely independent [164/166/170/172/174/175] and 6 confirmed `.wl` ports, and the
**USER DECIDED to "re-author all 6"** ports [165/167/168/169/171/173] for strict close-of-watch
dual-engine hygiene — Codex wrote the routes, `.py` untouched except 165's F2 numeric block, verify
confirmed all 6 genuinely independent. The route changes introduced NO new pinned numeric constant:
165's F2 prefactor checks [`Tm_pref≈1.2715890393387603`, `v_pref≈1.1428896163056477`,
`ratio_pref≈0.8987885086678338`, `prod_pref≈1.4532859092683434`] are re-derived cross-check targets, not
new carried constants; the 165/168/169 decimal→canonical-radius substitution [`1.77799353547498`→
`√(4107−100π²)/(10π)`] is a SCRIPT-side de-transcription of an EXISTING value, not a new or moved constant
[the canonical Family-1 radius is owned upstream at 121]; 169's 3 paper-comparison TARGET literals
[`0.758035078944663`/`1.00314310113848`/`1.88373219118005`] are preserved unchanged as paper-side values.
No checkpoint among the 12 stages, so cumulative checkpoint-constant provenance is unchanged from the IV.2
close [105 retained at the higher-bar standard]. Arbiter grep CLEAN — no `168π²`/`100π²` class. **0
sanctioned mirrors remain in V.1.**)
No carried checkpoint constant changed.

## Pass 2 — Batch V.2 (2026-06-09)

Pass-2 Batch V.2 (176-187): reviewed — ONE checkpoint in range (**185**, `microscopic_monomials`)
re-verified at the higher bar; **NO checkpoint constant changed, moved, or re-pinned; no new
postulated constant introduced.** 0 values misaligned batch-wide (~122 deliverable values; no
`paper_misalignment` anywhere → ZERO paper/notes edits); `material_change: false` on all 12; no new
postulated constant. **185** cleared the higher bar after its `.wl` was RE-AUTHORED so the load-bearing
monomial-exponent compilation is now DERIVED — a `monomialExponentVector` substituting each primitive
var → `Exp[logVar]`, taking `Log`/`PowerExpand`, and reading exponents via `Coefficient` — instead of
hand-coded identically to the `.py`; the checkpoint quantity `det
∂(Σ_tr,Σ_nt,Σ_eta)/∂(τ₁,κ_η,μ₁)=1+χ₀*` was **re-derived, NOT trusted as a literal or moved or
re-pinned** (`firstRatioDrift`, a `D[]`-slope, was kept and is acceptable now that the exponent
provenance is independent). The carried Stage-183 coefficients `2`/`11`/`9`/`4` and the symbolic
exponents/coefficients (`1+deltaU_*`, `1+chi0_*`, `E_*`, `F_*`) are PRESERVED. **Zero unexplained
literals remain in checkpoint 185's scripts.** **Cumulative checkpoint-constant provenance is unchanged
from the IV.2 close (105 retained at the higher-bar standard); no checkpoint constant moved.** (For
completeness: the 11 non-checkpoint findings stages introduced no new constant — the orchestrator
ground-truth `.wl`-vs-`.py` read found 8 ports [179/180/181/183/184/185/186/187] and 4 already
independent [176/177/178/182], and the **USER DECIDED to "re-author all 8" ports**; Codex wrote the
routes, verify confirmed all 8 genuinely independent [184 routes `R_target`'s complement law through a
separately-named `(1-ε_η)` identity object; 186 derives the eta-scaling from the physical
`ε_η=c_η U²/(K_U K_η)` monomial] — no emitted value moved. Arbiter grep CLEAN — no `168π²`/`100π²`
class. **0 sanctioned mirrors remain in V.2.**)
No carried checkpoint constant changed.

## Pass 2 — Batch V.3 (2026-06-09)

Pass-2 Batch V.3 (188-200): reviewed — ONE checkpoint in range (**200**,
`reference_free_home_stretch_theorem`) re-verified at the higher bar; **NO checkpoint constant changed,
moved, or re-pinned; no new postulated constant introduced — 200 pins NO new numeric constant (its
deliverables are SYMBOLIC).** 0 values misaligned batch-wide (171 deliverable values; no `paper_misalignment`
anywhere → ZERO paper/notes edits); `material_change: false` on all 13. **200's `.wl` was RE-AUTHORED but
the change is METHOD-ONLY** — the first-pass de-transliteration was found INSUFFICIENT (BOTH engines still
computed the load-bearing compiler matrix M_* by the SAME autodiff-Jacobian of the SAME log-ratios, SymPy
`q_pair.jacobian(Dvec)` vs Mathematica `Table[D[qPair,Dvec]]`; §III posited the same orbit closed forms; §V
used the same `Series`), so the USER-AUTHORIZED re-author DERIVES M_* from primitive monomial exponent-weight
vectors (§I), SOLVES the orbit via `Coefficient`→`LinearSolve` of the log-linear residual system (§III), and
linearizes Packet-A via a base-point derivative `D[chiPerturbed,eps]/.eps->0` (§V) — all landing on the SAME
values. 200's deliverables are SYMBOLIC and contain no free decimal literals: the carried Stage-192 M_*
(`Mexpected`) matrix, the §III mismatch chart law `q=((1+chi0_*)ln m_T, ln m_mu−ln m_K−F_*ln m_T, −ln m_K)`,
the §IV cocycle, and the §V Packet-A linear coeff `eps(5 eps_beta+dSigma0/(3S)+9 dSigma5/S)` are all
preserved; the source-anchored Stage-200 constants (`chi_Q=1`; the packet-length/pairing coefficient `2`;
the `3`/`5`/`9` carried from `chi_Q = 3(S β^5 + 9 Σ_5)/(3 S − Σ_0)`) already logged in the Stage 200 entry
above are REUSED, not moved. **Zero unexplained literals remain in checkpoint 200's scripts.** **Cumulative
checkpoint-constant provenance is unchanged from the IV.2 close (105 retained at the higher-bar standard); no
checkpoint constant moved.** (For completeness: the OTHER 12 non-checkpoint stages [188-199] introduced no
new constant — confirmed genuinely INDEPENDENT, 0 ports, by the orchestrator read + a dedicated calibration
agent against the 200 standard; only 200 had the Jacobian-vs-Jacobian sameness. Arbiter grep CLEAN — no
stale self-epoch 171–183 banner, no `168π²`/`100π²` class. **0 sanctioned mirrors remain in V.3.**)
No carried checkpoint constant changed.

## Pass 2 — Batch VI.1 (2026-06-09)

Pass-2 Batch VI.1 (201-218): reviewed — TWO checkpoints in range (**203**,
`free_quintuple_scalar_closure_slice_and_crossing_theorem`; **218**,
`full_support_cardinality_5_completion_and_local_mixed_ray_search_closure`) re-verified at the higher bar; **NO
checkpoint constant changed, moved, or re-pinned; no new postulated constant introduced.** 0 values misaligned
batch-wide (222 deliverable values; no `paper_misalignment` anywhere → ZERO substantive paper/notes edits);
`material_change: false` on all 18. **203** cleared the higher bar with its crossing theorem re-derived in-script —
the unique root `τ=1/2` of `32^(2τ-1)-1` (sign change `−31/32 → 31` across the root), confirmed genuinely independent
(log-additive + `Reduce` + target-monomial-invariance on the `.wl` vs power-multiplicative + `solveset` on the `.py`);
the carried structural constants (`widehat chi_Q = 1`; the Stage-192 quotient coefficient `2`; the symbolic
`(1+deltaU_*)/(1+chi0_*)`, `E_*`, `F_*`) are PRESERVED, not moved. **218** cleared the higher bar with its budget
ledger intact: **162 = 3⁴·2** (the lifted per-envelope Bézout bound), `324 = 2·162`, `750 = 5·5·5·6`,
`1500 = 2·750`, `1140` (carried Stage-215 support-≤4 global ledger), `1464 = 1140+324`, `2640` (projected/fallback
total) — every value internally consistent, appendix-exact, and CROSS-CHECKED at 217 (both engines derive 162; zero
surviving 179/230 across all 7 artifacts). The pass-2 orchestrator ground-truth read confirmed 218's pass-1 re-author
SUFFICIENT (M1 `Subsets`/`ContainsAll`/`Boole`/`Tally`/`2^5-2` vs `itertools`; M2-M3 `Resolve[ForAll,…,Reals]` real-QE
vs `simplify_logic`; M4 generated vs hand-listed witnesses; M5 benign shared budget arithmetic). **Zero unexplained
literals remain in checkpoints 203 and 218's scripts.** The Stage 218 provenance entry above already records
`162 = 3⁴·2` and the budget chain; nothing moved this batch. **Cumulative checkpoint-constant provenance is unchanged
from the V.3 close (200/203/218 retained at the higher-bar standard, 162 retained at 217/218); no checkpoint constant
moved.** (For completeness: the 16 non-checkpoint stages [201, 202, 204–217] introduced no new constant — confirmed
genuinely INDEPENDENT except 211, whose `.wl` was a confirmed PORT → USER-AUTHORIZED re-author to a `Resultant`-derived
eliminant route [`material_change: false`, no value moved], a CALIBRATION agent across the 207–216 ray-ranking family
finding NO other ports. Arbiter grep CLEAN — no `168π²`/`100π²` class. **0 sanctioned mirrors remain in VI.1.**)
No carried checkpoint constant changed.

## Pass 2 — Batch VII.1 (2026-06-09)

Pass-2 Batch VII.1 (219-230): reviewed — ONE checkpoint in range (**221**,
`resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window`) re-verified at the higher
bar; **NO checkpoint constant changed, moved, or re-pinned; no new postulated constant introduced.** 0 values
misaligned batch-wide (no `paper_misalignment` anywhere → ZERO substantive paper/notes edits); `material_change: false`
on all 12. **221** cleared the higher bar with deliverables that are SYMBOLIC, not pinned numerics — the
Breit–Wigner normal form, the perfect-square numerator `N = (A·G_W+R·G_U)²/Δ_Π²`, and the dispersive/no-free-lunch
ratio `|Re χ|/|Im χ| = |δ|/γ_*` — all PRESERVED, not moved. The pass-2 orchestrator ground-truth `.wl`-vs-`.py` read
confirmed 221's pass-1 re-author SUFFICIENT (the VI.1-218 outcome, UNLIKE V.3-200): §II DERIVES the Stage-220
derivative identity `dD_Pi/dPi=-N` via native `D[QPi/DeltaPi,portPi]` (an OUTPUT, vs the `.py` POSITING that form as an
INPUT — opposite information flow), plus `.wl`-only `Residue[]` extraction + generic `ComplexExpand` line-shape checks.
**Zero unexplained literals remain in checkpoint 221's scripts.** The benchmark numerics re-confirmed present in the
VII.1 batch — 230's thresholds `R_*≈1.229255438463336` / `δ_*≈0.723111617875019`, the 222/223 R_Q figures, the `i=h`
rigidity det factor `200+147π²`, the δ_1 coeff `196π²/(98π²−25)`, the crossover-cubic leading coeff `121ξ³` (all logged
at the VII.1 first-pass entry above) — HOLD, cross-engine-corroborated; none moved this batch. **Cumulative
checkpoint-constant provenance is unchanged from the VI.1 close (200/203/218/221 retained at the higher-bar standard);
no checkpoint constant moved.** (For completeness: the 11 non-checkpoint stages [219, 220, 222–230] introduced no new
constant — confirmed genuinely INDEPENDENT, 0 ports; **⭐ THE FIRST PASS-2 BATCH NEEDING ZERO SCRIPT CORRECTIONS**
[VII.1's first pass authored the 11 `.wl` FRESH as independent routes, so there was no porting residue], the all-clean
result EARNED via the orchestrator read of 220/221/226 + a calibrated independence-only sweep on the rest. Arbiter grep
CLEAN. **0 sanctioned mirrors remain in VII.1.**)
No carried checkpoint constant changed.

## Pass 2 — Batch VII.2 (2026-06-09)

Pass-2 Batch VII.2 (231-242): reviewed — TWO checkpoints in range (**239**, rigid-mouth physical normal form /
Cartesian orbit-lock; **242**, actual twin-support placement / coherent orbit-lock) re-verified at the higher bar;
**NO checkpoint constant changed, moved, or re-pinned; no new postulated constant introduced.** 0 values misaligned
batch-wide (no NEW `paper_misalignment` → ZERO substantive paper/notes edits); `material_change: false` on all 12.
Both checkpoints cleared the higher bar with deliverables that are SYMBOLIC, not pinned numerics — 239 the
orbit-lock map, 242 the twin-window ratio — both PRESERVED, not moved. The pass-2 orchestrator ground-truth
`.wl`-vs-`.py` read confirmed BOTH pass-1 re-authors SUFFICIENT (the VI.1-218 / VII.1-221 outcome, UNLIKE V.3-200):
**239** DERIVES the compiler via a `D[]`-Jacobian + native `PseudoInverse` left-inverse + `Reduce`/`Equivalent`
orbit-lock (OUTPUTS) vs the `.py` POSITING the literal `SrmDep` + `sp.solve` (opposite information flow); **242**
tests the load-bearing twin-window strict inclusion `C_mix < Pi_tr < 2·C_mix` STRICTLY on BOTH engines — the `.wl`
via `Resolve[ForAll,Reals]` QE + a `FullSimplify`-derived `4/3` vs the `.py`'s `nsimplify` + scalar compare.
**Zero unexplained literals remain in checkpoints 239 and 242's scripts.** The benchmark numerics re-confirmed
present in the VII.2 batch — the 3 first-pass notes-only typos HOLD, cross-engine-corroborated: 231 `dF/dξ` coeffs
`240→189` & `189→121`; 232 figure-of-merit prefactor `168→100` (ZERO surviving 168 — the recurring stale-168 family
is purged here); 241 `ϱ_WΛ` upper bound `193/369→125/369` (plus the 240 windows `ρ_α=4/3, ζ_req=1/3,
Pi_tr=(4/3)C_mix` and 241 `ϱ` windows `1/3, 125/369, 2/3, 250/441` independently corroborated) — none moved this
batch. **Cumulative checkpoint-constant provenance is unchanged from the VII.1 close (200/203/218/221/239/242
retained at the higher-bar standard); no checkpoint constant moved.** (For completeness: the 10 non-checkpoint
stages [231–238, 240, 241] introduced no new constant — confirmed genuinely INDEPENDENT, 0 ports; **THE SECOND
CONSECUTIVE PASS-2 BATCH NEEDING ZERO SCRIPT CORRECTIONS** [VII.2's first pass authored the 10 `.wl` FRESH as
independent routes + re-authored the 2 checkpoints, so there was no porting residue], the all-clean result EARNED via
the orchestrator read of 239/242 + a calibrated independence-only sweep on the rest. Arbiter grep CLEAN. **0
sanctioned mirrors remain in VII.2.**)
No carried checkpoint constant changed.

## Pass 2 — Batch VIII.1 (2026-06-09)

Pass-2 Batch VIII.1 (243-253): reviewed — THREE checkpoints in range (**243**, relaxed-branch lift / leakage-work
lane / non-rigid solve / short-range firewall; **248**, dynamic event chain / energy conservation / threshold speeds
/ near-top action; **253**, THE FINAL STAGE — lattice-turnover / calibration recovery / temperature ceiling / cold
survival) re-verified at the higher bar; **NO checkpoint constant changed, moved, or re-pinned; no new postulated
constant introduced.** 0 values misaligned batch-wide (no NEW `paper_misalignment` → ZERO substantive paper/notes
edits); `material_change: false` on all 11. All three checkpoints cleared the higher bar with deliverables that are
SYMBOLIC, not pinned numerics — 243 the leakage/non-rigid forms, 248 the threshold + near-top action, 253 the
calibration symbol + benchmark-only readbacks — all PRESERVED, not moved. The pass-2 orchestrator ground-truth
`.wl`-vs-`.py` read confirmed ALL THREE pass-1 re-authors SUFFICIENT (the VI.1-218 / VII.1-221 / VII.2-239&242
outcome, UNLIKE V.3-200): **243** derives U/V via BOTH `Solve` AND `LinearSolve` (asserts agreement) + an IBP-closure
cross-check + a `Series` asymptotic route, vs the `.py`'s `sp.integrate`/`sp.solve`/`sp.limit`; **248** POSITS the
threshold closed form and verifies it SATISFIES the defining energy equality with a non-vacuity guard
(`FreeQ[deltaNew,v0]`→fail), vs the `.py`'s `sp.solve` — a non-vacuous satisfaction route; **253** `Solve`s the
physical lattice balance + derives `chi_lambda` via `D[Log[V[r]],r]=2/r`, vs the `.py`'s division/hardcoded `2/r`.
**Zero unexplained literals remain in checkpoints 243, 248, and 253's scripts.** The benchmark numerics re-confirmed
present in the VIII.1 batch — the 5 first-pass notes-only typos HOLD, cross-engine or internally corroborated: 244
`196√2→128√2`; 247 Δ `210.17750000→142.17750000`; 248 `×168%→×100%` (ZERO surviving 168 — the recurring stale-168
family 148/232/248 is purged here); 253 benchmark `187.23361317→119.23361317` (the engines' value is
`119.23361317476524`, the stale `136.23361317476524`/`187.23361317` fully purged) + a_int
`10.95423247→10.95423248` (=4·K_turn) — none moved this batch; the 247/253 pass-1 misattributions to published cards
were re-overruled (cards clean/abstract). **Cumulative checkpoint-constant provenance is unchanged from the VII.2
close (200/203/218/221/239/242/243/248/253 retained at the higher-bar standard); no checkpoint constant moved.** (For
completeness: the 8 non-checkpoint stages [244/245/246/247/249/250/251/252] introduced no new constant — confirmed
genuinely INDEPENDENT, 0 ports; **THE THIRD CONSECUTIVE PASS-2 BATCH NEEDING ZERO SCRIPT CORRECTIONS** [VIII.1's
first pass authored the 8 `.wl` FRESH as independent routes + re-authored the 3 checkpoints, so there was no porting
residue], the all-clean result EARNED via the orchestrator read of 243/248/253 + a calibrated independence-only sweep
on the rest. Arbiter grep CLEAN. **0 sanctioned mirrors remain in VIII.1.** ⭐ **VIII.1 COMPLETES THE FULL SECOND
PASS — cumulative pass-2 253/253 (100%).**)
No carried checkpoint constant changed.