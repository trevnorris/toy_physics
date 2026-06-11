# Step 4C fix directive — align ingest path resolver with the bounded bundle writer (Codex applies, Claude reviews)

Context: your bulk-build change (session `019eb79c`) added `bounded_phase_b_stem` (core.py:252) + `phase_b_bundle_path` (core.py:277) so over-length bundle filenames are shortened with a sha1 suffix — this is REQUIRED (4 stage-252/253 candidates have 316-char stems that exceed the ext4 255-byte filename limit). The WRITER uses `phase_b_bundle_path` (bounded). But the INGEST-side resolver `provenance_path_for_parameter` (core.py:3787, unchanged) still suffix-matches `__{slug_param}.yaml` and falls back to the UNBOUNDED `{id}__{slug_param}.yaml`. For the 4 bounded candidates the recorded bundle ends in `..._<digest>.yaml`, so the suffix match fails and the fallback points at a nonexistent path → `phase-b-ingest` raises `provenance bundle does not exist` for all 4.

Confirmed: exactly 4 affected bundles, both targets (`..._sympy_audit` and `..._mathematica_audit`) at `fit_stage_252_*` and `fit_stage_253_*`. The other 911 candidates use stems ≤220 chars where `bounded_phase_b_stem` returns the legacy `{slug(id)}__{slug(param)}` verbatim, so they are unaffected.

## The fix (surgical, one function)

Make `provenance_path_for_parameter(env, entry, parameter_name)` resolve to the SAME path the writer produced:
1. Compute the canonical expected path `expected = phase_b_bundle_path(env, entry["id"], parameter_name)` (the bounded helper — identical to the legacy path for normal-length candidates).
2. Prefer a recorded path in `entry["paths"]["provenance"]` that equals `env.rel(expected)` (or whose basename equals `expected.name`); return it if present.
3. Otherwise return `expected`.

Drop the fragile `endswith("__{slug_param}.yaml")` suffix heuristic and the hand-rolled unbounded fallback string — both are wrong for bounded names. The new logic must remain byte-identical in OUTPUT path for the 911 normal candidates (since `phase_b_bundle_path` returns the legacy stem there).

Do NOT change `bounded_phase_b_stem`, `phase_b_bundle_path`, the writer, or any other behavior. Do NOT touch `paper/`/`notes/`/`scripts/`/`graph/`; FROZEN directive unchanged. Do NOT commit; leave only the `core.py` fix dirty.

## Verification contract (iterate until clean)
1. Apply the fix.
2. On a THROWAWAY copy of the artifact tree (or `--limit`/`--ids-file` real build that you then fully revert so the real tree is left with ONLY the `core.py` fix dirty): build the bundles for the 4 affected candidates (their ids contain `...stage253_physical_calibration_and_material_threshold_companion...sympy_audit` / `...mathematica_audit` at `fit_stage_252_*` and `fit_stage_253_*`), then assert that for each affected (candidate, parameter), `provenance_path_for_parameter(...)` returns a path that EXISTS and equals the path the writer wrote. Also assert it still returns the correct legacy path for 2–3 normal candidates (e.g. `fit_stage001_canonical_invariant_pair_determined` / `A` style). Paste the transcript.
3. Confirm `git status` shows only `.claude/skills/adversarial-audit/lib/core.py` dirty afterward (no leftover bundles / no manifest mutation).

Report: the new `provenance_path_for_parameter` body, the resolver-equivalence transcript (4 bounded resolve-to-existing + normal-candidate legacy-path preserved), and the clean-tree confirmation.
