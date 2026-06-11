# Step 4E close cleanup — remove orphaned unsound-merge dead code (Codex applies, Claude reviews)

Context: the Step-4 strict-dedup revision (session 019eb73a) dropped the value-compatibility fallback but left its helper functions orphaned in `core.py`. These functions ARE the exact unsound-merge logic that produced the bad dedup proposal; leaving them invites accidental re-wiring. Remove them.

## The change (one file: `.claude/skills/adversarial-audit/lib/core.py`)

CORRECTION (your first pass correctly halted): these 4 functions are NOT individually orphaned — they form a self-contained DEAD CLUSTER that only references ITSELF:
`connected_alias_components` (the cluster's only external entry point) → `records_alias_adjacent` → {`records_share_parameter`, `records_value_compatible`}.

Delete the WHOLE cluster (all 4 functions). The correctness condition is: the cluster's only external entry point, **`connected_alias_components`, must be referenced nowhere except its own `def`** (your grep already showed exactly that — only line 2677). Re-confirm that one fact, then delete all 4. After deletion, the only remaining mention of any of the 4 names is the descriptive status string at ~core.py:2900 (`"removed_fallbacks": "records_value_compatible/..."`) — that is a string literal, not a call. Reword that string to NOT name a deleted function (e.g. `"removed_fallbacks": "the value-compatibility / no-conflicting-value fallback and parameter-family union merging are disabled"`), preserving its meaning.

If `connected_alias_components` turns out to be referenced anywhere other than its `def` (i.e. something OUTSIDE the cluster calls into it), do NOT delete — report and stop.

Do NOT change any other behavior. Do NOT touch dedup/family/target-resolve/phase-b logic, the seeded-chain coverage check, or the extractor. Do NOT touch `paper/`/`notes/`/`scripts/`/`graph/`; FROZEN directive unchanged. Do NOT commit; leave only the `core.py` deletion dirty.

## Verification contract
1. Show the grep transcript proving each of the 4 functions is referenced ONLY at its own `def` (then deleted).
2. `python3 -m py_compile .claude/skills/adversarial-audit/lib/core.py` → exits 0.
3. Smoke-test that the read-only generators still work unchanged: run `phase-b-build-all --dry` and confirm it still reports `target_set_size: 915`, `target_error_count: 0` (read-only, mutates nothing). Run `target-resolve --out /tmp/_tl_check.yaml` (or a throwaway path) and confirm it still produces a sidecar with the same high/low split, then purge that throwaway.
4. Confirm `git status` shows only `.claude/skills/adversarial-audit/lib/core.py` dirty.

Report: the 4 grep transcripts + deletions (line ranges removed), py_compile result, the `--dry` counts, and the clean-tree confirmation.
