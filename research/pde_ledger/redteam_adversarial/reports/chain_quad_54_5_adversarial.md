# Phase C adversarial report — family `chain_quad_54_5` (54/5 GR-quadrupole external fit)

Batch: batch_c1_external_fit_and_high. Unit: 20 candidates (P0_target / Gamma_5 / N_Q angles, stages 019/022/025/038/106/189/193/195/197/223). Clean adversarial agent, per-member sub-verdicts. Benchmark: Peters 1964 (DOI 10.1103/PhysRev.136.B1224) / Maggiore Vol.1 Ch.3.

## Family verdict: NO fatal flaw
The 54/5 coefficient IS a fit — back-solved from the external GR value 2G/5c^5 via `sp.solve(mhat^2*P0*a^5/(27 c_s^5) == 2G/5c^5, P0)` (stage022 script; card 022:127) — but a fit honestly disclosed as such is not a fatal flaw under the layer-2 test. Across all 20 members it is consistently labeled "target" / "branch test" / StatusOpen, never relabeled derived/forced/non-tunable.

## Per-member verdicts (all NO)
| candidate_id | flaw? | reason |
|---|---|---|
| fit_stage019_P0_target_upstream_pin | NO | "P_{0,target}"; StatusExact only on the algebraic rewrite, value labeled "target" |
| fit_stage019_p0target_54_5_raw_retype | NO | hardcoded 54/5 = published_target retype, labeled "target" not "derived" |
| fit_stage_019_p0_target | NO | import; no card claim of internal derivation |
| fit_stage022_gamma_GR_quadrupole_target | NO | origin; card 022:121 "GR quadrupole target", check#4 "branch test, not a proof" |
| fit_stage022_p0_quadrupole_normalization_target | NO | origin back-solve; disclosed as branch test (022:143) |
| fit_stage_022_p_0 | NO | duplicate origin angle; same honest disclosure |
| fit_stage025_one_mode_normalization_target | NO | card 025:5 StatusOpen "whether branch realizes formula at required value" |
| fit_stage025_target_coefficient_54_5 | NO | 54/5=(2/5)*27 confirmed; "carried forward", StatusOpen |
| fit_stage_025_p0_compact | NO | P0_compact internal; RHS is the GR target, StatusOpen |
| fit_stage038_mhat_p0_fitted_scales | NO | object is R_target not P0_target (mislabel) but still a disclosed 54/5 import |
| fit_stage106_coeffs_fixed_to_targets | NO | Gamma5_target=2G/5c^5 hardcoded; "GR-Burke-Thorne target", StatusOpen |
| fit_stage_106_gamma_5 | NO | clean published_target |
| fit_stage189_normalization_target_2g_5c5 | NO | back-solve to 2G/5c^5; "the product the PDE still has to hit" |
| fit_stage_189_p_0 | NO | duplicate |
| fit_stage193_p0_target_surface_frozen_without_derivation_edge | NO | StatusExactClosure on the surface/collapse; value defined in notes, "target" wording |
| fit_stage195_p0_target_canonical_value | NO | notes define P_0^target := 54/5, Gamma5^target := 2G/5c^5 explicitly as targets |
| fit_stage197_p0_target_quadrupole_normalization | NO | StatusOpen; "exact" applies to Delta_branch=0<=>chi_Q=1 algebra, not the value |
| fit_stage223_p0_target_calibration_coefficient | NO | card 223:9 "the universal target P_{0,target}", kept symbolic |
| fit_stage223_p0_target_import_unprovenanced | NO | provenanced upstream at 022; "universal target" flagged imported |
| fit_stage223_universal_p0_target | NO | "universal" = external GR universality, imported not derived |

## Disclosure assessment (layer-2 crux): ADEQUATE across all 20
Decisive card-face evidence (read directly, not from provenance synthesis): 022:121 "The GR quadrupole target"; 022:143 check#4 "a branch test, not a proof of branch realization"; 025:5 "Whether the actual branch realizes the formula at the required value is StatusOpen"; 193/195/197 carry StatusExactClosure (+StatusOpen at 197) where "exact" attaches to the algebraic collapse/identity, never to the value; 223:9 "the universal target P_{0,target}" as an explicit import.
**Weakest point (NOT an overclaim):** 193/195/197 name P0_target's GR origin only in the notes, not on the card face (derivation-ledger line 13 "states the exact isotropic target surface"); "exact" there is true of the surface algebra and StatusOpen is present, so no card applies derived/forced/exact to the back-solved number.

## Proof of look
Benchmarks fam_0187_p0_target / fam_0093_gamma_5 / fam_0174_n_q (bench_2e5805a359c4): coefficient G/(5c^5), 2/5 radiation-reaction (Peters 1964 / Maggiore Ch.3). Coefficient check: 54/5 = (2/5)*27, the 27 = structural fingerprint a^5/(27 c_s^5), the 2/5 matches GR exactly digit-for-digit (022:127 => m̂0^2 P0 = (2G/5c^5)(27 c_s^5/a^5) = 54 G c_s^5/(5 a^5 c^5)). Fit-vs-derive: 1 normalization knob, fixed by `solve(...==gamma_GR)` => FIT; 20 members = ONE imported number reused, no overdetermination, one knob per radiative phenomenon. Red-flags: rational-vs-rational (no transcendental mismatch); no over-digit match; only validation is symbolic residual==0 (tests algebra not physics) but cards carry StatusOpen, not over-credited.
