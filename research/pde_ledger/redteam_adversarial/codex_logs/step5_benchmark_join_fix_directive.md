# Step 5 pre-fan-out fix — Phase C benchmark→candidate join is broken (Codex applies, Claude reviews)

Context: the Phase C pilot render exposed a defect. Every rendered Phase C adversarial prompt has `## Benchmarks` = `[]`, even for the GR-quadrupole and CODATA candidates that DO have sourced benchmark entries. The frozen directive requires external matches be checked against the sourced `benchmarks.yaml` entry ("do not adjudicate a match from model memory"), so this empty section would force all 26 external-fit agents to judge from memory — defeating the benchmark-sourcing infrastructure. Must be fixed before the batch-1 fan-out.

## Root cause (verify, then fix)

`benchmarks.yaml` entries are keyed by **`family_id`** (e.g. `fam_0187_p0_target`, `fam_0093_gamma_5`, `fam_0174_n_q`, `fam_0457_m_s`, `fam_0458_m_s_ratio_1836`) with `candidate_id: null`. But `load_benchmarks_for_candidate(env, candidate_id)` in `.claude/skills/adversarial-audit/lib/core.py` (~line 3963) joins strictly on `e.get("candidate_id") == candidate_id`. Since real benchmark entries carry no `candidate_id`, nothing ever matches → `[]` for every candidate.

The family map `research/pde_ledger/redteam_adversarial/provenance/_family_map.yaml` carries two candidate→family maps, both **list-valued** (a candidate can belong to several families): `primary_candidate_family_map` and `candidate_family_map`. For `fit_stage022_p0_quadrupole_normalization_target`, `primary_candidate_family_map` = `['fam_0187_p0_target']` — exactly the benchmark family. So the join must go through the family map.

## The change (one file: `.claude/skills/adversarial-audit/lib/core.py`)

Fix `load_benchmarks_for_candidate` so it resolves the candidate's families from `_family_map.yaml` (union of `primary_candidate_family_map[candidate_id]` and `candidate_family_map[candidate_id]`, each of which may be a string or a list) and returns every benchmark entry whose `family_id` is in that family set.

Requirements / acceptance criteria (pin behavior, you choose the implementation):
1. A candidate whose family set intersects a benchmark entry's `family_id` gets that entry. Concretely: `fit_stage022_p0_quadrupole_normalization_target` → the `fam_0187_p0_target` (Peters 1964 / Maggiore quadrupole) entry; `fit_stage022_gamma_GR_quadrupole_target` → `fam_0093_gamma_5`; `fit_stage023_nq_target_gr_quadrupole_match` → `fam_0174_n_q`; `fit_stage250_goldilocks_inputs` → `fam_0458_m_s_ratio_1836`; `fit_stage250_benchmark_masses_and_window` → `fam_0457_m_s`. Exactly **26** of the 65 batch-1 candidates should resolve to ≥1 benchmark entry (the GR-quadrupole P0/γ5/N_Q chain + the two CODATA candidates).
2. **Preserve the legacy `candidate_id` match** (dry-run benchmark placeholders created by `update_benchmarks_for_candidate` are candidate-keyed and must still attach). An entry that matches by EITHER family OR candidate_id is returned; dedup so an entry is never returned twice.
3. A candidate with no benchmark family (e.g. `fit_stage089_pe_window_upstream_misattribution`, `fit_stage163_Rstar_quarter`) still gets `[]` — no false attach.
4. Robust if `_family_map.yaml` is absent or a map key is missing: fall back to the current candidate_id-only behavior, no crash.
5. Do not change benchmark INGEST, the family map, the directive, or any other render section. No new files. Loading the family map once per call is acceptable (renders are interactive, low-volume); do not add global mutable caches.

## Verification contract (read-only function calls — do NOT advance any candidate status)

Do NOT re-run `phase-c-render` on real candidates for verification (it mutates `provenance_built → audit_pending`). Verify by importing the function and calling it directly:

1. `python3 -m py_compile .claude/skills/adversarial-audit/lib/core.py` → exits 0.
2. A small read-only python harness that imports `core`, builds the `Env`, and calls `load_benchmarks_for_candidate(env, cid)` for: `fit_stage022_p0_quadrupole_normalization_target` (expect 1 entry, family `fam_0187_p0_target`, source_citation mentions Peters 1964), `fit_stage250_goldilocks_inputs` (expect family `fam_0458_m_s_ratio_1836`, CODATA 1836.15267343), and `fit_stage089_pe_window_upstream_misattribution` (expect `[]`). Print each result's `id`+`family_id`. The harness reads only; it must not call any WRITE_COMMANDS and must not load+save the manifest.
3. Print the count of batch-1 candidates (the 65-candidate `published_target ∪ high` union; you can reuse the union logic by scanning the bundles) for which `load_benchmarks_for_candidate` returns non-empty → expect **26**.
4. `git status` shows only `.claude/skills/adversarial-audit/lib/core.py` dirty (re-rendered prompts under `tmp_prompts/` are transient and may be ignored).

Do NOT commit. Report: the diff (function before/after), the py_compile result, the harness output for the 3 probe candidates + the count=26, and the clean-tree confirmation.
