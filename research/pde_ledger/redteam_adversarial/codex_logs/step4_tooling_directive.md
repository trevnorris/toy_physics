# Step 4 tooling directive — Phase B plumbing (Codex applies, Claude reviews)

Context: Step 3 (Phase A) is complete and committed (`0e5d177`): `MANIFEST.yaml` holds **922 candidates, all `scanned`**, with `fit_insertion_points.yaml` rendered. Step 4 is Phase B — the parameter-value provenance ledger + sourced benchmarks. The Step-3-close Claude+Codex consult (record: `BATCHING_DECISIONS.md`) agreed five decisions; this directive builds the tooling those decisions require, which does not yet exist in the skill. The existing `phase-b-build` is per-candidate, emits a *mechanical* evidence bundle with empty `origin_claims`/`constraints`/`provenance_findings` and has **no path to ingest the agent synthesis back**; `update_benchmarks_for_candidate` only emits a dry-run placeholder; there is no dedup/alias or family-map tooling and `render_batches` clobbers the consult decisions on every write.

This is a TOOLING build only. Do NOT apply the real dedup map, do NOT run audit/provenance/benchmark agents, do NOT adjudicate, do NOT edit anything under `paper/`, `notes/`, `scripts/`, `graph/`. Phase B judgment (synthesis, benchmark sourcing, the dedup-map review) is the orchestrator/agent side and happens in later sub-steps.

## Binding constraints (identical to the Step 2/Step 3 builds)

- Allowed write set: `.claude/skills/adversarial-audit/**` and `research/pde_ledger/redteam_adversarial/**`. Additive-only edits to the `adversarial:` section of `research/pde_ledger/.redteam-config.yaml` ONLY if a config key is genuinely needed.
- `docs/adversarial_audit_directive.md` is FROZEN — do not touch. Zero edits to `.claude/skills/redteam-audit/`.
- Every MANIFEST / benchmarks / provenance write that mutates campaign state goes through the existing flock (`_manifest_locked` / `require_manifest_lock`). Read-only proposal generators (Items A, C) must NOT take the lock and must NOT mutate the manifest.
- `timeout 600` everywhere; exit 124 = reformulate the approach, never raise the cap.
- No OS temp files for anything an agent will read — render under `redteam_adversarial/tmp_prompts/`. YAML/markdown for every agent-facing file (no JSON for LLM I/O).
- Do NOT commit. Leave the tree dirty for orchestrator review.

## Design intent (from `docs/adversarial_audit_deployment.md` — build to this, do not re-derive)

Phase B builds **parameter-value provenance only** (taxonomy layer 3 of §3.3: where a numerical value / previously-symbolic parameter was *introduced and what constrained it* — NOT file lineage, NOT stage attribution alone). Per candidate/parameter the genealogy is: (1) where each parameter was first introduced; (2) at what stage/line it was assigned a value; (3) **what constrained the value — internal consistency, a published target, or a free choice** (this is the fit-vs-derive discriminator Phase C consumes); (4) which downstream stages depend on the claim. Two binding rules, bought with real failures: **no attribution without opening the per-stage `notes/` source**, and **a self-contradictory origin claim is itself a finding** (the χ_Q=1 episode). Benchmarks: one entry per claimed external match, with the literature value, its source citation, and how obtained (web lookup / textbook / CODATA) — "matches GR" is **never** adjudicated from model memory.

---

## Item A (blocker) — `dedup-propose`: read-only alias-map proposal (consult decision D1)

A new subcommand `dedup-propose [--out PATH]` (default `redteam_adversarial/provenance/_dedup_proposal.yaml`).

- Reads MANIFEST **read-only** — takes NO lock, mutates NOTHING.
- Groups the `scanned` candidates into alias clusters. Alias rule: **same `anchor_stage` AND same normalized parameter family AND same-or-adjacent cited line** (define and document "adjacent": e.g. identical `(path,line)`, or same `path` within a small line window — your choice, justify it). The same physical insertion point routinely appears as 2–4 candidates (mechanical key + per-modality agent keys); this collapses those.
- **No deletions, ever.** The proposal only *names* a canonical per group and lists its aliases.
- Canonical-selection rule must be deterministic and documented (e.g. most modality attributions, tie-broken by id sort).
- A separate `ambiguous` section for groups that share stage+parameter but differ in value or insertion site — these stay SEPARATE (not aliased), listed for the reviewer.
- The proposal YAML must be reviewable by a clean agent without re-deriving: each group carries `canonical_id`, `aliases: [ids]`, the shared `{anchor_stage, parameter_family, citation}`, and per-member `{candidate_key, parameter_names, citation, modalities}` evidence. Include summary counts (groups, total aliases, resulting canonical count, ambiguous groups).
- Codex measured 408 same-stage/same-parameter duplicate groups; the post-dedup canonical estimate is ~400–550. Report the actual numbers your rule yields.

## Item B (blocker) — `apply-alias-map <map.yaml>`: flock-gated apply (D1)

A new subcommand that consumes a (reviewed) alias map and applies it under the manifest lock.

- For each alias: set `duplicate_of: <canonical_id>` on the alias entry; merge its modality attribution + citations + fragments **onto the canonical** with the same dedup semantics as `merge_fragments`. Aliases **keep** their manifest entries (no deletion) and are excluded from later Phase B/Phase C work selection (they inherit the canonical's verdict at Phase C close).
- Never regress the canonical's status; never silently drop a listed id. Unknown ids / already-aliased ids fail INFORMATIVELY.
- Idempotent: re-applying the same map is a safe no-op.
- Re-render `fit_insertion_points.yaml` and `BATCHES.md` under the lock.

## Item C — `family-build`: parameter-value family map (D2)

A new subcommand `family-build [--out PATH]` (default `redteam_adversarial/provenance/_family_map.yaml`).

- Operates over the **canonical** candidate set (aliases excluded if `duplicate_of` present; if run pre-apply, operate over all `scanned` and note that it should be re-run post-apply — document this).
- Derives families mechanically: group canonicals by **normalized `parameter_names` + shared literal values across stages** (e.g. the 54/5 target is one family covering every stage it touches). Deterministic, documented rule.
- Output `_family_map.yaml`: `families: [{family_id, parameter_names, representative_values, member_candidate_ids, stages_touched}]`. Every canonical maps to **≥1 family or an explicit singleton family** (no canonical left unmapped). Include summary counts.
- Read-only on the manifest (writes only the map file); reviewable by a clean agent.

## Item D (blocker) — `phase-b-ingest <synthesis.yaml> [...]`: flock-gated synthesis ingest

Mirror `phase-a-ingest`'s validation/merge/lock discipline for Phase B synthesis. This is the missing return path: agents read the `phase-b-build` prompt + `notes/`, emit a synthesis YAML, and this command merges it into the candidate's `provenance/*.yaml`.

- Accept agent-emitted synthesis YAML keyed by `candidate_id` (and `parameter_name` to target one provenance file, since `phase-b-build` emits one file per parameter). Schema per parameter:
  - `origin_claims: [{parameter, introduced_at_stage, introduced_at_line, citation:{path,line,excerpt}}]`
  - `constraints: [{parameter, constraint_kind: internal_consistency | published_target | free_choice, evidence_citation:{path,line,excerpt}, note}]`
  - `downstream_dependents: [stage, ...]`
  - `provenance_findings: [{type, severity, summary, citations:[...]}]` (self-contradiction findings live here)
- **Binding validation (load-bearing):** every `origin_claims` and `constraints` entry MUST carry a citation whose `path` is under `notes/` (the "no attribution without the notes" rule). Reject — informatively, all-or-nothing per file, naming the file/claim/what's missing — any claim lacking a `notes/` citation, and any malformed file/claim. Never silent-drop.
- Merge the synthesis into the existing `provenance/<candidate_id>__<param>.yaml` (built by `phase-b-build`) under the lock; do not destroy the mechanical `source_evidence`/`graph_context` already there. Set a `synthesis_status: complete` field (distinct from the mechanical-only bundle, which is `synthesis_status: pending`).
- Status: a real candidate reaches `provenance_built` only once synthesis is ingested. Reconcile with the existing `phase-b-build`, which currently transitions straight to `provenance_built`: change `phase-b-build` so a real (non-dry-run) candidate stops at a pre-synthesis state (mechanical bundle written, `synthesis_status: pending`) and `phase-b-ingest` completes the transition to `provenance_built`. **Preserve the existing dry-run behavior** (`dry-run --stages 003 104 105` must still reach `provenance_built` and the documented dry-run expectations) — dry-run fixtures may keep the current straight-through path.
- Re-render `fit_insertion_points.yaml` / `BATCHES.md` under the lock as the other write commands do.

## Item E (blocker) — `benchmark-ingest <bench.yaml> [...]` + fix `update_benchmarks_for_candidate`

- New flock-gated `benchmark-ingest`: accept sourced benchmark YAML, entries keyed by `family_id` or `candidate_id`. Per entry: `{claim, value, source_type: web_lookup | textbook | CODATA, source_citation (REQUIRED — URL / DOI / CODATA identifier / textbook ref+page), obtained_note, model_value (optional)}`.
  - **Reject (informatively) any entry lacking a real `source_citation`** — this enforces the "never adjudicate a match from model memory" rule. Reject malformed entries. Never silent-drop.
  - Merge into `benchmarks.yaml` under the lock (dedup by entry id; do not clobber unrelated entries).
- Fix `update_benchmarks_for_candidate`: it currently removes all of a candidate's entries and only re-adds for the dry-run chi_Q placeholder, so a `phase-b-build` after a real `benchmark-ingest` would silently wipe real sourced entries. Change it to **preserve real (non-dry-run) ingested entries**; it may continue to manage only the dry-run placeholder under dry-run. Document the new contract.

## Item F — extend `render_batches` (`BATCHES.md`)

`BATCHES.md` is clobbered on every manifest write, losing the D-decisions (which is why `BATCHING_DECISIONS.md` exists as the durable record). Extend the renderer so a regenerated `BATCHES.md` reflects **real campaign state**, not an empty stub: render the dedup canonical/alias counts and the family groupings (read `_family_map.yaml` if present), plus a header pointer to `BATCHING_DECISIONS.md` as the authoritative consult record. If a per-candidate `batch:` field is present, render the batch assignment; otherwise render "batch assignment pending (Step 5)". Keep it modest — this is faithful reflection, not a new planner.

## SKILL.md runbook

Update the Phase B section of `.claude/skills/adversarial-audit/SKILL.md` to document the real Step-4 flow end to end: `dedup-propose` → (orchestrator review) → `apply-alias-map` → `family-build` → per-family `phase-b-build` → provenance agent synthesis → `phase-b-ingest` → benchmark agents → `benchmark-ingest`. Note the binding rules (notes/-citation required; self-contradiction is a finding; benchmarks never from model memory).

## Verification contract (iterate until everything exits 0)

Order matters so the tree ends in the right state:

1. Build all items.
2. Run the **read-only** generators on the REAL manifest (they mutate nothing): `dedup-propose` and `family-build`. Leave `_dedup_proposal.yaml` and `_family_map.yaml` in the tree — these are the proposals the orchestrator/agents will review next sub-step. Report: dedup group count, total aliases, resulting canonical count, ambiguous-group count; family count, singleton count, any unmapped canonical (must be zero).
3. Demonstrate `apply-alias-map`, `phase-b-ingest`, and `benchmark-ingest` using **DRY-RUN-marked throwaway fixtures only** (e.g. `dry-run --stages 104 105` to create fixture candidates, `phase-b-build` them, then ingest a hand-written dry-run-marked synthesis + benchmark fixture; for `apply-alias-map` use a scratch fixture map over dry-run candidates). Then `purge-dry-run` / remove them. The manifest, `provenance/` (except the two real read-only proposals from step 2), and `benchmarks.yaml` must end this step containing **ZERO test/fixture artifacts**. Show the purge in the transcript.
4. Confirm `dry-run --stages 003 104 105` still satisfies every documented dry-run expectation (003 vacuous; 104/105 real candidate; advances scanned→provenance_built→audit_pending; provenance slice generated; Phase C prompt rendered; no binding verdicts) — i.e. the `phase-b-build` state-machine change did not break dry-run. Purge afterward.
5. Do NOT apply the real alias map; do NOT run any agents; do NOT adjudicate; do NOT edit `paper/`/`notes/`/`scripts/`/`graph/`.

Report: per-item summary with file:line of every change; design choices (adjacency rule, canonical selection, family grouping rule, the `phase-b-build`/`phase-b-ingest` state split); the real dedup/family proposal counts from step 2; the fixture + purge transcript evidence for step 3; the dry-run re-confirmation from step 4; anything not satisfied and why.
