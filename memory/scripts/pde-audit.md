---
title: PDE audit script catalog
type: script_catalog
status: current
sources:
  - research/pde_audit/CLAIM_CHECK_MATRIX.md
  - research/pde_audit/run_all.sh
  - research/pde_audit/scripts/run_all_audits.sh
  - research/pde_audit/mathematica/run_all_audits.sh
  - research/pde_audit/scripts/stage_v2_01_parent_wall_action_sympy_audit.py
  - research/pde_audit/scripts/stage_v2_12_stf_angular_source_map_sympy_audit.py
  - research/pde_audit/scripts/stage_v2_16_branch_freeze_no_refit_sympy_audit.py
  - research/pde_audit/scripts/stage_v2_23_minimal_open_throat_branch_solver.py
  - research/pde_audit/scripts/stage_v2_23_mesh_convergence_audit.py
  - research/pde_audit/mathematica/stage_v2_23_formula_mathematica_audit.wl
last_updated: 2026-09-03
---

## `research/pde_audit/run_all.sh`

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

- Language: shell
- Source unit: `pde-audit`
- Role: runner
- Computes/builds: top-level fixture, Python-audit, simulation, root-JSON-cleanliness, and Mathematica run records; then coordinates referee-summary construction.
- Inputs: the audit directory and its fixture verifier, Python runner, simulation runner, Mathematica runner, and summary writer.
- Emits: `output/fixture_manifest.txt`, `output/python_audits.txt`, `output/simulation.txt`, `output/root_json_clean.txt`, `output/mathematica_audits.txt`, and the summary writer’s `output/_summary.txt`; each captured check includes metadata, stdout, and `EXIT_CODE`.
- Guards/controls: propagates child failures, rejects root-level generated JSON files, and incorporates summary-writer status. This covers harness execution, not physical correctness.
- Invocation: `not recorded`
- Interpretation source: `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `# PDE V2 Claim-Check Matrix`
- Related topics: `memory/sources/pde-audit.md`

Sources:

- `research/pde_audit/run_all.sh` — nearby identifier `run_top_check`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Reproducibility Artifacts`

## `research/pde_audit/scripts/run_all_audits.sh`

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

- Language: shell
- Source unit: `pde-audit`
- Role: runner
- Computes/builds: per-script Python execution records, aggregate pass/fail counts, and a machine-readable suite summary.
- Inputs: optional substring filter; `stage_v2_*.py`; fixture and artifact directories. It routes manifest/output arguments separately to V2-21, V2-22A, V2-22B, V2-22C, V2-23 mesh convergence, and V2-24, and runs an additional invalid-hard-cap V2-22B scenario.
- Emits: `scripts/output/*.txt`, `scripts/output/_summary.txt`, and artifacts beneath `scripts/output/artifacts/`; it invokes the summary writer for `scripts/output/_summary.json`.
- Guards/controls: records exit codes, counts failures, exercises the explicit invalid-hard-cap scenario, and fails if an audit or summary generation fails. Coverage is runner-level.
- Invocation: `not recorded`
- Interpretation source: `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`
- Related topics: `memory/sources/pde-audit.md`

Sources:

- `research/pde_audit/scripts/run_all_audits.sh` — nearby identifier `run_capture`
- `research/pde_audit/scripts/run_all_audits.sh` — nearby identifier `args_for_script`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`

## `research/pde_audit/mathematica/run_all_audits.sh`

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

- Language: shell
- Source unit: `pde-audit`
- Role: runner
- Computes/builds: environment-probe and selected Mathematica-mirror execution records, aggregate counts, and a Mathematica suite summary.
- Inputs: optional substring filter, environment probe, `stage_v2_*_mathematica_audit.wl`, fixture manifest, and summary writer.
- Emits: `mathematica/output/_environment.txt`, per-script `mathematica/output/*.txt`, `mathematica/output/_summary.txt`, and the invoked summary writer’s `_summary.json`.
- Guards/controls: captures elapsed time and exit status, filters mirrors by basename, and fails on environment, mirror, or summary-generation failure. The matrix identifies these mirrors as secondary execution coverage rather than independent derivations.
- Invocation: `not recorded`
- Interpretation source: `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `# PDE V2 Claim-Check Matrix`
- Related topics: `memory/sources/pde-audit.md`

Sources:

- `research/pde_audit/mathematica/run_all_audits.sh` — nearby identifier `math -script`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `# PDE V2 Claim-Check Matrix`

## `research/pde_audit/scripts/stage_v2_01_parent_wall_action_sympy_audit.py`

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

- Language: Python
- Source unit: `pde-audit`
- Role: validator
- Computes/builds: symbolic Euler derivatives for confinement-only and quadratic wall actions, modal wall operators, boundary momenta, source-free Hamiltonian density, positivity conditions, and the lowest axisymmetric two-mode reduction.
- Inputs: symbolic fields and coefficients declared inside `main`; no external fixture is visible in the supplied excerpt.
- Emits: `PASS:` stdout labels and failure details through assertions; no tracked output path is visible.
- Guards/controls: compares symbolic expressions for equality and checks whether the confinement-only source contains wall-field derivatives. The script explicitly does not solve the moving-throat PDE.
- Invocation: `not recorded`
- Interpretation source: `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`
- Related topics: `memory/sources/pde-audit.md`

Sources:

- `research/pde_audit/scripts/stage_v2_01_parent_wall_action_sympy_audit.py` — function `main`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`

## `research/pde_audit/scripts/stage_v2_12_stf_angular_source_map_sympy_audit.py`

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

- Language: Python
- Source unit: `pde-audit`
- Role: validator
- Computes/builds: a real STF \(l=2\) basis, tensor and angular Gram matrices, source-map and quadrupole-reconstruction residuals, source norms, grouped-metric comparisons, and the \(\Pi\)-to-\(q\) convention map.
- Inputs: symbolic STF tensors, angular variables, and source coefficients defined in the script; no external fixture is visible in the supplied excerpt.
- Emits: matrices, residuals, coefficients, named check results, and summary-verdict stdout; no tracked output path is visible.
- Guards/controls: checks trace-freedom, basis normalization, reconstruction, grouped norms, and convention-map rank. Its own description excludes radial, axial, and outgoing-normalization questions.
- Invocation: `not recorded`
- Interpretation source: `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`
- Related topics: `memory/sources/pde-audit.md`

Sources:

- `research/pde_audit/scripts/stage_v2_12_stf_angular_source_map_sympy_audit.py` — function `tr_inner`
- `research/pde_audit/scripts/stage_v2_12_stf_angular_source_map_sympy_audit.py` — nearby identifier `The STF angular source-map has no angular normalization defect`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`

## `research/pde_audit/scripts/stage_v2_16_branch_freeze_no_refit_sympy_audit.py`

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

- Language: Python
- Source unit: `pde-audit`
- Role: validator
- Computes/builds: Boolean dependency incidence and transitive closure for a freeze-protocol DAG, grouped-P2 target residual expressions, interface identities, and a post-hoc fitting rank diagnostic.
- Inputs: symbolic protocol nodes, target variables, residual expressions, and internally constructed freeze data visible around `main`; no external fixture path is supplied by the excerpt.
- Emits: check records and serialized diagnostic material from `main`; no tracked result path is visible in the supplied excerpt.
- Guards/controls: checks topological ordering and absence of target/residual feedback into branch definitions, then distinguishes frozen comparison from post-hoc parameter fitting. It does not decide whether a PDE branch is physically realized.
- Invocation: `not recorded`
- Interpretation source: `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`
- Related topics: `memory/sources/pde-audit.md`

Sources:

- `research/pde_audit/scripts/stage_v2_16_branch_freeze_no_refit_sympy_audit.py` — function `bool_incidence`
- `research/pde_audit/scripts/stage_v2_16_branch_freeze_no_refit_sympy_audit.py` — function `boolean_transitive_closure`
- `research/pde_audit/scripts/stage_v2_16_branch_freeze_no_refit_sympy_audit.py` — function `main`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`

## `research/pde_audit/scripts/stage_v2_23_minimal_open_throat_branch_solver.py`

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

- Language: Python
- Source unit: `pde-audit`
- Role: builder
- Computes/builds: a target-blind reduced open-throat shape and mode branch, overlap-derived coefficients, response/prefactor quantities, target residuals, a branch-freeze hash, gate results, and symbolic formula checks.
- Inputs: `--grid-points` (default `181`) plus constants, geometry, boundary protocol, couplings, and coefficient family frozen inside `build_branch`.
- Emits: the supplied excerpt shows construction of an observable packet and creation of the artifact directory; the exact completed write operation and full output-path set are not visible.
- Guards/controls: requires at least 25 grid points and evaluates open, stability, outgoing-transfer, and target tolerances; symbolic formula checks run separately. The matrix limits this to a reduced linear-FEM prototype.
- Invocation: `not recorded`
- Interpretation source: `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`
- Related topics: `memory/sources/pde-audit.md`

Sources:

- `research/pde_audit/scripts/stage_v2_23_minimal_open_throat_branch_solver.py` — function `solve_open_shape_profile`
- `research/pde_audit/scripts/stage_v2_23_minimal_open_throat_branch_solver.py` — function `build_branch`
- `research/pde_audit/scripts/stage_v2_23_minimal_open_throat_branch_solver.py` — function `evaluate_gates`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`

## `research/pde_audit/scripts/stage_v2_23_mesh_convergence_audit.py`

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

- Language: Python
- Source unit: `pde-audit`
- Role: comparator
- Computes/builds: repeated reduced-solver summaries, coefficient and residual values, pairwise relative deltas, observed convergence orders, gate aggregates, and maximum mode/stationarity residual.
- Inputs: `--grids` (default `91,121,181,241`), `--finest-pair-rel-tol`, `--residual-tol`, optional `--out-json`, and the V2-23 solver module.
- Emits: an optional JSON report; the Python runner supplies `scripts/output/artifacts/stage_v2_23_mesh_convergence_report.json`.
- Guards/controls: requires at least three strictly increasing grids of at least 25 points; checks finite values, solver gates, residual tolerance, freeze identity, and finest-grid-pair convergence over the named coefficient set.
- Invocation: `not recorded`
- Interpretation source: `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`
- Related topics: `memory/sources/pde-audit.md`

Sources:

- `research/pde_audit/scripts/stage_v2_23_mesh_convergence_audit.py` — function `summarize_run`
- `research/pde_audit/scripts/stage_v2_23_mesh_convergence_audit.py` — function `convergence_orders`
- `research/pde_audit/scripts/run_all_audits.sh` — function `args_for_script`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`

## `research/pde_audit/mathematica/stage_v2_23_formula_mathematica_audit.wl`

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

- Language: Wolfram Language
- Source unit: `pde-audit`
- Role: comparator
- Computes/builds: series coefficients \(u_2,u_4,P_0,P_2,P_4\), one-pole and constant-prefactor residuals, and normalization-equivalence expressions; it also reads the V2-23 tolerance and observable packets.
- Inputs: shared Mathematica audit helpers plus `stage_v2_23_tolerance_report.json` and `stage_v2_23_observable_packet.json`.
- Emits: named check output and `All Stage V2-23 Mathematica checks passed.` to stdout; the Mathematica runner captures per-script output.
- Guards/controls: compares formula identities, checks packet schemas, requires open and stability gates, requires an honest target-packet failure, and requires a nonempty dominant-residual label. This mirror shares the project’s formulas and packet exports and is not an independent derivation.
- Invocation: `not recorded`
- Interpretation source: `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`
- Related topics: `memory/sources/pde-audit.md`

Sources:

- `research/pde_audit/mathematica/stage_v2_23_formula_mathematica_audit.wl` — `subbanner["Formula audit"]`
- `research/pde_audit/mathematica/stage_v2_23_formula_mathematica_audit.wl` — `subbanner["Solver residual packet"]`
- `research/pde_audit/mathematica/run_all_audits.sh` — nearby identifier `math -script`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`
