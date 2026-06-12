You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage247_benchmark_pinned_to_paper_figures__phase_c_adversarial.md

This is a FAMILY audit (parameter-value family `chain_calibration_245_253`) — cover EVERY member with its OWN sub-verdict.


## Focus
Calibration-end external imports — KEEP MEMBER LINES SEPARATE (these are DISTINCT questions): (a) lambda_L=0.26971918 at stage 247 is back-solved so V_eff(r_soft) reproduces the recorded Session-I softening point V_eff^sess=1.74701126 — notes honest ('audit point, not a theorem') but the SCRIPT wording ('independently derived... pinned to paper figures', lines ~239/248) is misleading since the forward check is tautological; (b) CODATA 1836.15267343 proton/electron mass ratio at stage 250 imported as proton-proxy m_s=mu_eta. Central Q per member: is each honestly disclosed as an imported benchmark, or overclaimed as derived? Check CODATA against benchmarks.yaml fam_0457/fam_0458.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage247_benchmark_pinned_to_paper_figures (param=lambda_L): redteam_adversarial/provenance/fit_stage247_benchmark_pinned_to_paper_figures__lambda_l.yaml
- fit_stage247_lambda_L_backsolve (param=lambda_L): redteam_adversarial/provenance/fit_stage247_lambda_l_backsolve__lambda_l.yaml
- fit_stage247_lambda_L_closure (param=lambda_L): redteam_adversarial/provenance/fit_stage247_lambda_l_closure__lambda_l.yaml
- fit_stage247_lambda_l_fixed_by_recorded_point (param=lambda_L): redteam_adversarial/provenance/fit_stage247_lambda_l_fixed_by_recorded_point__lambda_l.yaml
- fit_stage250_benchmark_masses_and_window (param=chi_peak): redteam_adversarial/provenance/fit_stage250_benchmark_masses_and_window__chi_peak.yaml
- fit_stage250_benchmark_masses_and_window (param=chi_raw): redteam_adversarial/provenance/fit_stage250_benchmark_masses_and_window__chi_raw.yaml
- fit_stage250_benchmark_masses_and_window (param=E_max): redteam_adversarial/provenance/fit_stage250_benchmark_masses_and_window__e_max.yaml
- fit_stage250_benchmark_masses_and_window (param=g_UV): redteam_adversarial/provenance/fit_stage250_benchmark_masses_and_window__g_uv.yaml
- fit_stage250_benchmark_masses_and_window (param=m_s): redteam_adversarial/provenance/fit_stage250_benchmark_masses_and_window__m_s.yaml
- fit_stage250_benchmark_masses_and_window (param=mu_eta): redteam_adversarial/provenance/fit_stage250_benchmark_masses_and_window__mu_eta.yaml
- fit_stage250_benchmark_masses_and_window (param=t_transit_max): redteam_adversarial/provenance/fit_stage250_benchmark_masses_and_window__t_transit_max.yaml
- fit_stage250_benchmark_masses_and_window (param=t_transit_min): redteam_adversarial/provenance/fit_stage250_benchmark_masses_and_window__t_transit_min.yaml
- fit_stage_247_lambda_l_soft (param=lambda_L): redteam_adversarial/provenance/fit_stage_247_lambda_l_soft__lambda_l.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open what you need to ground each finding.

## Benchmarks (use sourced entries, NEVER memory)
If this family claims an external match: GR mass-quadrupole 2G/5c^5 -> benchmarks.yaml family_ids fam_0187_p0_target / fam_0093_gamma_5 / fam_0174_n_q; CODATA proton/electron mass ratio 1836.15267343 -> fam_0457_m_s / fam_0458_m_s_ratio_1836. Read research/pde_ledger/redteam_adversarial/benchmarks.yaml. The rendered prompt already embeds the relevant entry if applicable.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id | YES/NO (fatal flaw?) | one-clause reason` (one row per member).
2. Family verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific external value / broken step / overclaiming card+line.
3. Disclosure assessment (layer-2 crux): honestly labeled target/benchmark/open, or overclaimed derived/forced/exact? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which benchmark entries (if any external-match claim) checked digit-by-digit, and which fit-vs-derive tests you ran.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts — read and reason.
