# PDE-Ledger Corpus Exclusion Validation

Date: 2026-09-01
Target commit inspected: `60e31504b8425840cf6b5e814a0c518d68ab2db6`

## Policy decision

The research memory does not catalog any source under:

- `research/pde_ledger/`
- `research/pde_ledger_v2/`
- `research/pde_ledger_v3/`

The similarly named but separate roots `research/pde/` and
`research/pde_audit/`, plus relevant Stage-1 and other software results,
remain eligible.

## Deterministic enforcement

The scanner rejects hard-excluded paths used as exact source-unit members,
selector roots, derived direct sources, supporting lineages, or active Atlas
migration reanchors. Broad selectors also filter these roots before matching.
Prefix matching is segment-safe, so a path such as
`research/pde_ledger_v30/example.md` is not excluded by name similarity.

The selected configuration resolves to 25 source units and 20 derived tasks.
No configured unit, direct source, supporting lineage, or active migration
requirement uses a PDE-ledger path.

## Retired served content

The previous S10 capsule, v3 status topic, and v3 script catalog were removed.
The S10 unit/page/statement records were also removed from served state.
Ignored sealed transactions may retain historical packets, but they are not
served memory and are not committed.

## Checks

The following checks passed after the policy change:

```text
python3 memory/tools/memory.py status
python3 -m unittest \
  memory.tools.tests.test_memory.MemoryToolTests.test_hard_excluded_roots_reject_unit_members_and_selectors \
  memory.tools.tests.test_memory.MemoryToolTests.test_hard_excluded_roots_reject_direct_sources_and_supporting_lineages \
  memory.tools.tests.test_memory.MemoryToolTests.test_hard_excluded_roots_reject_migration_original_sources \
  memory.tools.tests.test_memory.MemoryToolTests.test_hard_excluded_roots_filter_broad_selectors_and_are_segment_safe
git diff --check
```

The focused test run completed four tests successfully. `memory status`
reported no PDE-ledger source unit or derived task.

## Benchmark migration

Benchmark IDs `B03`, `B04`, `B06`, and `B07` were retired without reuse.
Replacement cases `B11` and `B12` cover material/nonlinear closure and a
scoped Stage-1 source revision using only selected sources. Existing Python
and Wolfram coverage remains in `B09` and `B10`.
