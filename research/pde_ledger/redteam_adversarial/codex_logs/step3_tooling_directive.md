# Step 3 tooling directive — Phase A campaign blockers (Codex applies, Claude reviews)

Context: Step 3 (the real Phase A fit-insertion-point scan over stages 001–253) has started. Two runtime gaps in `.claude/skills/adversarial-audit/` block it. Same binding constraints as the Step 2 build: allowed write set = `.claude/skills/adversarial-audit/**` and `research/pde_ledger/redteam_adversarial/**` (plus additive-only edits to the `adversarial:` section of `research/pde_ledger/.redteam-config.yaml` if a config key is genuinely needed); `docs/adversarial_audit_directive.md` is FROZEN — do not touch; zero edits to `.claude/skills/redteam-audit/`; every manifest write through the flock; `timeout 600` everywhere, exit 124 = failure to reformulate, never raise the cap; no OS temp files for anything an agent reads — use `redteam_adversarial/tmp_prompts/`; YAML/markdown for agent-facing files. Do NOT commit — leave the tree dirty for orchestrator review.

## Item A (HIGH, blocker) — full-band `phase-a-scan` cannot finish inside the 600s cap

Evidence: `phase-a-scan --stages 001…253` was killed at exactly 600s (exit 124) having written nothing (fragments dir empty, manifest unchanged). Cause: `graph_scan` (core.py) calls `query_graph_source` per (stage × source-role × path-variant), each a `subprocess.run` of `graph/query_graph.py` that reloads the full atlas YAML — over 1000 cold subprocess invocations for the full band.

Requirement: a SINGLE invocation `phase-a-scan --stages 001 … 253` must complete comfortably under the 600s cap (target: low single-digit minutes). The design choice is yours — load the atlas graph once per scan invocation and answer all source queries in-process, batch the queries into one subprocess call, cache by source path — whatever you choose:
- Fragment/gap SEMANTICS must not change: same fragment fields, same `graph_gap` rule (genuine gaps still logged per deployment doc §7), same per-modality fragment files and union behavior.
- No silent band-limiting, no sampling, no skipping stages to get under the cap.
- Also add `--prefix` support to `phase-a-scan` (default `phase_a`, mirroring `render-phase-a-prompts`) so banded invocations cannot clobber each other's fragment files — defense in depth even after the perf fix.

## Item B (HIGH, blocker) — no ingestion path for agent-emitted Phase A fragments

Evidence: the four blind modality agent prompts instruct agents to emit YAML fragments, and SKILL.md says "Agent fragment incorporation for Step 3 must preserve the same fragment schema and blindness rule" — but `run_phase_a` unions only the fragments its own mechanical scanners produce. There is no CLI path that merges an agent fragment file into `fit_insertion_points.yaml`/MANIFEST, and hand-editing the manifest violates the skill's own lock invariant.

Requirement: a new flock-gated command `phase-a-ingest <fragment.yaml> [<fragment.yaml> ...]` that:
1. Accepts agent-emitted fragment YAML files conforming to the modality-prompt schema (`modality:` + `candidates:` with `candidate_key`, `anchor_stage`, `parameter_names`, `citation:{path,line,excerpt}`, `reason`). Normalize what agents plausibly omit (fill `citation.stage` from `anchor_stage`; `role` optional). Malformed files / malformed candidates fail INFORMATIVELY (say which file, which candidate, what's missing) — never silent-drop.
2. Applies NO `candidate_parameters()` / denylist / fit-signal filtering to agent fragments. This is load-bearing: the review punch list established that bare-greek-symbol fit claims are invisible to the mechanical scanner — the agent modalities ARE the coverage for those. Accept agent-provided `parameter_names` verbatim. If a candidate arrives with empty `parameter_names`, do not silently drop or silently placeholder it: reject that candidate with an informative message listing it (the orchestrator re-asks the agent).
3. Merges with the same key semantics as `merge_fragments`: an agent fragment whose `candidate_key` matches an existing campaign candidate merges INTO it (append modality attribution, citations, fragments — dedup as merge_fragments does); a new key creates a new manifest entry with status `scanned` (id `fit_` + key). Never regress a candidate already past `scanned` — append attribution only and leave its status alone.
4. Accepted modality labels: the four standard ones PLUS `completeness_critic` (the critic's finds are ingested as a fifth scan through this same path).
5. Re-renders `fit_insertion_points.yaml` and `BATCHES.md` under the lock, exactly as `phase-a-scan` does.
6. Companion: a `render-critic [--prefix NAME]` command (or an explicit flag on `phase-a-ingest` — your choice) that re-renders the completeness-critic prompt from the CURRENT manifest union, listing every fragment file present under `phase_a_fragments/` (mechanical + agent), because the critic must run AFTER agent ingestion and `phase-a-scan`'s critic render predates it.
7. SKILL.md Phase A runbook updated to document the real Step-3 flow end to end: mechanical scan → blind agent fragments → `phase-a-ingest` → `render-critic` → critic agent → ingest critic finds via the same command.

## Verification contract (iterate until everything exits 0)

Order matters so the tree ends in the real campaign state:
1. Build both items.
2. Demonstrate ingest + render-critic + the `--prefix` collision guard using DRY-RUN-marked throwaway fixtures only (e.g. a hand-marked test fragment with `dry_run: true` semantics or a scratch prefix), then `purge-dry-run` / remove them — the manifest and `phase_a_fragments/` must end this step containing ZERO test artifacts. Show the purge in the transcript.
3. As the final acceptance step, run the real `phase-a-scan --stages 001 … 253` once, timed; report wall-clock, candidate count, per-modality fragment counts, and graph_gap count. Leave those real artifacts in the tree.
4. Do not run any agents, do not adjudicate anything, do not edit anything under `paper/`, `notes/`, `scripts/` — Phase A judgment is the orchestrator's side.

Report: per-item summary with file:line of every change, design choices made, the timing evidence for Item A, the ingest/critic/purge transcript evidence for Item B, anything not satisfied and why.
