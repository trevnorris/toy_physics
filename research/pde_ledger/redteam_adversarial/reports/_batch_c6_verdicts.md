# Phase C — Batch 6 verdicts (targeted free_choice sweep, stages 200–239)

**Headline: 5/5 units NO FATAL FLAW (24/24 candidates NO). No YES, no defense, no adjudication, no user gate.**

Batch 6 of the PC6 targeted free_choice sweep — the **200–239 "finite-throat illustrative-sample-slice / branch-calibration / figure-of-merit" band**. 24 free_choice candidates (~42 param-checks across two compound bundles) → 5 disjoint stage-cluster units, per-member sub-verdicts. The fit-vs-derive test for free_choice = **posit-dressed-as-derived overclaim** (the stage_192 failure mode): is a CHOSEN value labeled derived/forced/exact/non-tunable/canonical-as-unique in any paper card or per-stage notes line? `Benchmarks: []` is correct for free_choice (no external-match claim). Roster: `_phase_c_batch6_roster.yaml`. Dispatch prompts: `tmp_prompts/_dispatch_b6_u0[1-5]_*.md`.

## Per-unit verdict table

| unit | stages | members | verdict | crux |
|---|---|---|---|---|
| u01_linearization_path | 200/202/203 | 4 | **NO** | eps_beta symbolic linearization direction; C_tr_target a "declared" target-orbit value ("uniquely_fixed" scoped to δ_U, not the target); beta_path / KW_t script-only illustrative IVT-path coordinates. |
| u02_barrier_benchmark_Vknown ⚠SENSITIVE | 222/223 | 3 | **NO** | V_known = "illustrative reduced barrier benchmark"; card never names it; DeltaV_req symbolically inherits illustrative status. Audit-correction (free_choice, not external) holds. |
| u03_primitive_sample_slice ⚠SENSITIVE | 222/225/228 | 5 | **NO** | 20+ sample-slice values all "illustrative … on one admissible sample slice"; card:15 "the displayed slice is a numerical placement"; 225/228 "Frozen input carried forward". |
| u04_fom_normalization_232 ⚠SENSITIVE | 232 | 5 | **NO** | L_over_a=20 = "carried thin-wall lock ℓ/a=1/20" (origin 073 "reference-branch numerical freeze, not a new theorem") — NOT 37/20; N_Q=1 "canonical outgoing-normalization condition"; λ_μ=1 "for the benchmark"; Ξ=100 bare definitional literal. |
| u05_branch_calibration_samples | 224/233/234/236/237 | 7 | **NO** | lambda_20 imported lane-weight convention; Delta_rm / Delta_norm=0 explicit conditional calibration ("If…"/"On a calibrated branch with…"); Xi_1 "prescribed"/"any chosen" free target; T_U "at fixed T_U" held-fixed axis; Rratio_base=7/6 script-only sample; `hat` = scanner artifact. |

**Member tally: 24 NO / 0 YES.** Statuses: all 24 `audit_pending → audited`. Manifest after batch: **189 audited + 2 verdict_logged** (A1 stage_192 + A8 γ₀), 700 provenance_built, 7 scanned.

## Decisive finding of the band

The entire 200–239 band is honestly-disclosed free_choice with **zero posit-as-derived overclaim**. The band's content is dominated by:
- **Illustrative numerical placements** on a posited finite-throat admissible sample branch (the primitive tuple λ_B/R/U/W, M, Ω_U/W, ϖ, K, a, c_s) — explicitly flagged illustrative at notes222:1/7/434 and card stage_222.tex:15 ("the displayed slice is a numerical placement"), carried forward as "Frozen input" at 225/228. Every "exact" word in these stages attaches to the **quartic pole polynomial / residue-linewidth cancellation / survival gate "exact inside that reduced closure"**, never to the slice numbers.
- **Normalization conventions** — N_Q=1 ("the canonical outgoing-normalization condition", notes232:314/316), λ_μ=1 ("for the benchmark", notes232:149).
- **Script-only path/sample coordinates** — beta_path=2^(2τ−1), KW_t=KW0·e^{t·κ_W} (Stage-203 IVT/crossing demo), Rratio_base=7/6 (Stage-237 identity-verification sample).
- **Conditional calibration / held-fixed axes** — Delta_norm=0 ("On a calibrated branch with Δ_norm=0", notes233:203-205), T_U "at fixed T_U" (notes236:22).
- **Prescribed free targets** — Xi_1 "for a prescribed first-order defect value Ξ₁=ξ" realized by a "canonical least-deformation family that realizes any chosen Ξ₁" (notes234:328/372); "canonical" labels the minimum-norm family construction, not a forced value.

## Orchestrator spot-check record (SENSITIVE units, distrust-all-clean backstop)

- **u02 V_known** (the audit-correction site). Ratified: stage_222.tex:13/15 ("illustrative admissible slice"; "the displayed slice is a numerical placement"; "exact in the primitive closure" scoped to the pole polynomial + survival gate, card never names V_known); notes222:350 "Carry forward the same illustrative reduced barrier benchmark", :436 "The benchmark ΔV_req(1) used here is illustrative. The true local barrier requirement must be pulled back from the actual same-charge branch"; notes223:329 carries it forward illustrative. **No card/notes line dresses V_known as derived/forced/external.**
- **u03 sample slice** (high param count). Ratified: notes222:1/7 "illustrative evaluations on one admissible sample slice of the same exact formulas", :7/434 "The sample couplings are illustrative. The real moving-throat PDE must determine the branch-level overlap data λ_B,λ_U,λ_W,λ_R on the actual branch", :6 "with a=c_s=1"; stage_228.tex:9 / notes228:60/64 "Frozen input carried forward" / "Keep the same explicit finite-throat one-port branch used in Stages 223–227"; stage_225.tex:9 / notes225:60/64 likewise. **The illustrative/placement disclosure covers every sample-slice param.**
- **u04 stage-232 literals**. Ratified: notes232:42 defines Λ_ℓ := L/ℓ (not L/a); stage073:42/56/70 establish ℓ/a=1/20 is a "**reference-branch numerical freeze, not a new theorem**" and Λ_ℓ = (L/a)/(ℓ/a) = (37/20)/(1/20) = 37 — so the scanner's "L_over_a=20" is the inverse thin-wall ratio a/ℓ=20, **NOT the batch-1 aspect ratio 37/20** (not re-litigated); N_Q=1 "the canonical outgoing-normalization condition" (notes232:314/316); λ_μ=1 "for the benchmark" (notes232:149); Ξ_χ = 100 Θ_w^(χ) Λ_ℓ² a bare definitional literal (notes232:153). **No literal dressed as derived/forced.**

## OPTIONAL close-out items (all non-fatal; appended to STEP6_CLOSEOUT_PUNCHLIST.md)

1. **Scanner artifacts (candidate-list cleanup, audit-internal, no paper fix):** `fit_stage_224_hat` (the `hat` param is a LaTeX fragment of `\hat m_0`, not a real parameter — mhat_0 is an undetermined lower-bounded branch factor at 224); `fit_stage_200_stale_output` (empty bundle, zero constraints). Same class as batch-2's V_0f token.
2. **Stage-232 scanner dedup / mis-name (D-class, audit-internal):** the two `L_over_a` slugs (`fit_stage232_family1_geometry_normalization` + `fit_stage232_lambda_ell_geometry_literals`) are duplicate scanner candidates for the same literal `20`, and the param name `L_over_a` is a stale_provenance_anchor mis-label of the inverse thin-wall ratio a/ℓ=20. No content fix.
3. **Provenance-fill (Tier-2/3 class, non-overclaim gaps):** `Xi_prefactor_100` (the literal 100 has no in-stage derivation; FoM structure traces to Stage 066); `V_known` float genealogy below 222 undocumented (disclosed illustrative, so not an overclaim); `beta_path` base-2 / `KW_t` exp-ramp / `Rratio_base`=7/6 live only in scripts, not notes (script-only illustrative).
4. **Completeness-critic targets — overclaim-flavored SLUGS that are constraint_kind=internal_consistency (OUT OF SCOPE for the free_choice sweep; to be covered by the internal_consistency completeness-critic spot-check, NOT silently dropped):** 223 `lambda02_rq_dual_engine_overclaim` (R_Q…; slug literally says "overclaim"), 225 `xk_xm_forced_zero` (x_K), 224 `grouped_signature_exact` (b_P0_slope), 232 `known_5pn_data_injection` (rho_alpha_max — confirmed INTERNAL EM-projected chain, not a publication, per the memory audit correction), 200 `mismatch_chart_rederived_not_posited` (m_star), 236 `equal_drift_ray_forced` (k_eta), 228 `k_compat_exact_stiffness` (k_compat), 230 `rigid_slope_carry` (s_minus_den), 224 `kill_test_budgets_slice_anchored` / `stage223_carryover_ceilings` (P_crit…). All confirmed internal_consistency by direct constraint_kind read this batch.

## Non-fatal observation (not a finding)

The u03 audit noted a value-divergence: the Stage-222 illustrative slice (K=3.0, Ω_W=1.4, ϖ=2.0) differs from the carried branch at 223+ (K=K_compat≈24.47 compatibility-solved, Ω_W=7/5). The 223 K_compat is `internal_consistency` (out of scope here) and the slice is honestly bounded as illustrative; this is a stage-boundary role difference, not an overclaim.

## Coverage log (no silent caps)

- **In scope and audited:** all 24 free_choice candidates with a free_choice parameter in band 200–239 (per the band-filtered enumerator; 13 stages: 200, 202, 203, 222, 223, 224, 225, 228, 232, 233, 234, 236, 237).
- **Out of scope by design (logged, not audited this batch):** the internal_consistency candidates in-band (→ completeness-critic, item 4 above); the published_target candidates in-band — 223 `p0_target_*` (P0_target), 228 `rq_requirement_threshold` (rq_req) — already audited in batch 1's `published_target ∪ HIGH` cluster.
- **Remaining free_choice after batch 6:** band 240–259 ≈ 34 candidates → batch 7 (the last free_choice batch), then the internal_consistency completeness-critic spot-check.
