Micro-fix iteration (same Step 2 build). The fix review verified F1-F8 all FIXED, overall CLEAN. Three small advisories remain; fix the first, and the other two only if trivial. Same constraints as before; do not commit.

N1 (required): `core.py:1113` — `source_evidence_for_candidate` still contains hardcoded `"Stage 87"` / `"Stage 104"` literals on the CAMPAIGN path (it runs for every candidate, adding chi_Q-episode evidence lines to unrelated candidates). Generalize or remove: evidence selection must derive from the candidate's own parameters/citations, never from episode-specific literals. Grep core.py afterward to confirm zero episode-specific stage/parameter literals remain on any campaign-reachable path (dry-run fixture code behind dry_run flags is fine).

N2 (optional, if trivial): `artifact_has_dry_run_marker` content-heuristic sweep in `purge-dry-run all` could over-delete a historical non-dry-run file that quotes `dry_run: true` early in its text (e.g., a future codex log). Cheap hardening: exempt `codex_logs/` from the marker sweep, or require markers in YAML front-matter position rather than anywhere in the first 20000 chars.

N3 (optional, if trivial): `graph_context_for_candidate` (`core.py:1137`) falls back to `anchor_stages[0]` for every citation; multi-stage candidates get mis-staged fallback variants. If citations can carry their stage, thread it through.

Verify: py_compile + re-run `dry-run --stages 003 104 105` clean (all original bullets), confirm the chi_q candidate's provenance evidence is unchanged in substance (its own citations still found via general logic). Report per-item with file:line.
