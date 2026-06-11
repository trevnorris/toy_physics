# Step 4 — Phase B chain-curation inputs (from the data-review agent + post-apply state)

Durable record of the seeded-chain corrections and the 813-low-confidence triage, to feed the Phase B genealogy pass (task 4C). Source: clean-agent data review of the regenerated proposals (2026-06-11) + orchestrator post-apply inspection. The PRIMARY families and the dedup map are sound and APPLIED (915 canonicals, 7 aliases, 0 unmapped, overlap 0). Only the seeded D3 chain OVERLAYS need agent curation; final chain membership is determined during genealogy from `notes/`.

## Seeded-chain corrections

### chain_quad_54_5 (42 members) — PRUNE false positives
Concept matching over-pulled. Remove (not genuine 54/5 quadrupole-coefficient members):
- `fit_stage_177_p_0` — P1/P0 ratio (value "38"), not the 54/5 target.
- `fit_stage032_mhat_minus_completely_fixed` — minus-branch source map, different object from mhat_0.
- `fit_stage038_mhat_p0_fitted_scales` — atlas graph directive node, not a derivation.
- `fit_stage225_K_compat_P0_target_pin` — a compatibility *sample* value (~0.00207), not the coefficient.
Borderline (genealogy agent decides from notes): `fit_stage022_mhat0_unity_branch` (parametrizes P0_target but doesn't carry 54/5), `fit_stage100_factorization_forced_script` (enforces mhat₀²·χ_Q·N_Q=1, not the 54/5 target).
NOTE: the chain's `value_divergence` buckets `'1'`/`'38'` are ARTIFACTS of these false positives, not real 54/5 variants. The genuine variant to record is `54` (denominator-stripped retype, e.g. stage019) vs `54/5`.

### chain_aspect_37_20 (44 members) — one false positive + mislabeled divergence
- Remove `fit_stage072_canonical_controls_exact_thresholds` — generic branch thresholds; stage072 precedes the 37/20 fixing at stage073.
- Its `value_divergence` finding conflates FIVE distinct parameters (g_star 0.758, r_F1 1.778, epsilon_r 1/20, L_over_a 37/20, Sigma0_can 4.651) as if one family's variants — it is a multi-parameter overlay. The only genuinely actionable within-parameter divergences are: r_F1 {1.778 ↔ Sqrt[4107−100π²]/(10π)} and L_over_a {37/20 ↔ 20}.

### Genuine coverage misses (ADD via genealogy agent, both confirmed genuine)
- `chain_barrier_222_224` ← `fit_stage224_grouped_signature_exact` (primary target `b_P0_slope`, absent from the chain concept list; key + stage say it belongs — confirm b_P0_slope isn't a genuinely separate barrier sub-parameter).
- `chain_calibration_245_253` ← `fit_stage248_xi_turn_lambda_th_carried_hardcodes` (Ξ_turn=0.34437471; key literally says "lambda_th_carried_hardcodes", stage 248 is a calibration carrier).

### Alias-driven FALSE coverage alarms (IGNORE — canonical + value already in chain)
Post-apply, the coverage check flags these aliased members as "missing" but their canonicals remain in the chain, so coverage/value is intact. Do NOT act on them:
- chain_aspect_37_20: `fit_stage073_wall_fraction_epsilon_r`, `fit_stage121_aspect_ratio_37_over_20`
- chain_barrier_222_224: `fit_stage222_vknown_barrier_benchmark`
- chain_chi_Q_norm: `fit_stage105_chiQ_eq_1_by_matching_card`, `fit_stage105_chi_q_canonical_fix`
- chain_sigma0_transport: `fit_stage168_sigma0_num_digit_variant` (the …876 carrier; value preserved via its canonical — the …867/…876 divergence still holds)

## 813 low-confidence triage (identity is SOUND; ~13/15 correct, low only for missing value)
- **TIER 1 — review before/within Phase B (~60–90):** low-conf candidates inside a seeded chain or a multi-member primary family with a value_divergence finding (the chain_quad_54_5 prunes + the 2 coverage adds + `multiple_plausible_targets` records sitting in chains, e.g. stage072). These shape per-family genealogy.
- **TIER 2 — light value-fill (~232 `multiple_plausible_targets` non-singletons):** identity is a guess; ~1-in-7 may be wrong (observed miss: `fit_stage252_three_to_one_split` → should be split_surface, not f_lat).
- **TIER 3 — defer (~500+ low-conf singletons):** identity from a meaningful candidate_key matching a real parameter; per-candidate genealogy unaffected by the missing value. Fill values opportunistically.

## ~23 genuine duplicates left unmerged (safe under-merge, self-healing)
Confidence-asymmetry pairs (a HIGH-conf valued record beside a LOW-conf same-parameter sibling with no extracted value) were not merged by the strict rule. Re-cluster these within-family during the value-fill pass. Clearest HIGH-HIGH cross-file pair worth merging: stage108 chi_Q=1 (`fit_stage108_chiQ_robust_claim` + `fit_stage_108_chi_q`, both HIGH value 1, different .md reports, line distance >3).

## Tooling nits to fold into the Step-4 close cleanup (LOW, non-blocking)
1. Dead code: `records_value_compatible`, `records_alias_adjacent`, `connected_alias_components`, `records_share_parameter` are orphaned after the fallback removal — delete to prevent accidental re-wiring (this is the exact unsound-merge logic).
2. Seeded-chain coverage check should resolve aliases→canonicals before flagging missing (it produced 6 false alarms post-apply, listed above).
3. `frac_54_5` extraction branch is digit-pattern over-greedy; `mention-fallback` RHS scan is the loosest extraction path — both low false-positive risk, worth tightening or asserting against a verified expression.
