---
name: adversarial-audit
description: Layer-2 adversarial fit-vs-derive audit for the PDE ledger. Builds fit-insertion-point candidates, parameter-value provenance, adversarial prompts, Codex defense prompts, and adjudication records without modifying the existing redteam-audit skill.
allowed-tools: Bash, Read, Edit, Write, Agent
user_invocable: true
---

# Adversarial Audit

This skill drives Step 2+ of the PDE ledger layer-2 audit. It is parallel to, and independent from, `.claude/skills/redteam-audit/`. Do not edit the redteam-audit skill while using this one.

The per-pass adversarial directive is frozen at `docs/adversarial_audit_directive.md`. This skill consumes it by path or exact inclusion only.

## State Machine

The tracked unit is a fit-insertion-point candidate, not a stage.

```
pending / scanned / provenance_built / audit_pending / audited /
defense_pending / defended / adjudicating / verdict_logged / blocked
```

## Commands

Set:

```bash
AA=/var/projects/toy_physics/.claude/skills/adversarial-audit/lib/adversarial.sh
```

Run from either `/var/projects/toy_physics` or `/var/projects/toy_physics/research/pde_ledger`; the CLI locates `research/pde_ledger/.redteam-config.yaml`.

- `init` creates `redteam_adversarial/`, `MANIFEST.yaml`, `fit_insertion_points.yaml`, `benchmarks.yaml`, `BATCHES.md`, and required subdirectories.
- `render-phase-a-prompts --stages NNN ... [--prefix NAME]` renders the four blind Phase A modality prompts only, under `tmp_prompts/`.
- `phase-a-scan --stages NNN ...` runs the broad mechanical Phase A modalities, writes YAML fragments, unions candidates, and renders all four modality prompts plus the completeness-critic prompt.
- `phase-b-build <candidate_id>` builds a mechanical parameter-value evidence bundle, renders the Phase B provenance-builder prompt, and updates manifest status to `provenance_built`.
- `phase-c-render <candidate_id>` renders the adversarial prompt under `redteam_adversarial/tmp_prompts/` and updates status to `audit_pending`.
- `set-status <candidate_id> <status> [--note TEXT]` advances a candidate through the state machine under the manifest lock; `blocked` is reachable from any state.
- `codex-defense <candidate_id> <parameter> [iter]` renders a defense prompt and invokes Codex, resuming the session keyed by parameter.
- `graph <show|neighborhood|paths|source> ...` wraps `graph/query_graph.py` using the configured graph path.
- `status` or `summary` prints manifest counts and dry-run state.
- `candidate-info <candidate_id>` dumps one manifest entry.
- `dry-run --stages 003 104 105` exercises Phase A -> Phase B -> Phase C prompt rendering without binding verdicts.
- `purge-dry-run <id|all>` removes dry-run manifest entries and dry-run artifacts in one step.

## Runbooks

### Phase A Scan

For Step 3 real campaign work, do not hand-select a family of stages or parameters. The real sweep is the full configured stage range, currently 001 through 253. Render the blind prompts for the campaign band first:

```bash
$AA render-phase-a-prompts --stages $(seq -f "%03g" 1 253) --prefix campaign_001_253
```

Launch four clean agents, one per rendered prompt: numeric-literal, claim-label, graph, and existing-provenance. Each agent remains blind to the other modality outputs and emits only a YAML fragment. The mechanical `phase-a-scan` path is broad and not seeded to any parameter family; it provides the union baseline and records per-modality attribution. Agent fragment incorporation for Step 3 must preserve the same fragment schema and blindness rule.

For a bounded scan or dry-run fixture:

1. Run:
   ```bash
   $AA phase-a-scan --stages 003 104 105
   ```
2. Inspect `redteam_adversarial/phase_a_fragments/*.yaml` and `redteam_adversarial/tmp_prompts/*_phase_a_*.md`. Each modality is blind to the others.
3. Inspect `redteam_adversarial/fit_insertion_points.yaml`. The union preserves `modality_attribution` and `file_line_citations`.
4. Read the rendered completeness-critic prompt under `redteam_adversarial/tmp_prompts/` and launch a clean agent only for the real campaign sweep.

### Phase B Build

1. Choose a scanned candidate:
   ```bash
   $AA candidate-info <candidate_id>
   ```
2. Run:
   ```bash
   $AA phase-b-build <candidate_id>
   ```
3. Verify the generated `provenance/*.yaml` is parameter-value provenance only. It is a mechanical evidence bundle, not interpretive synthesis.
4. Give the rendered `tmp_prompts/*__phase_b__*.md` prompt to a clean provenance agent. Every attribution must cite the per-stage `notes/` source. If origin claims contradict each other, keep the contradiction as a finding. If graph lookup misses a source, keep `graph_gap`.

### Phase C Batch

1. Ensure candidate status is `provenance_built`.
2. Render:
   ```bash
   $AA phase-c-render <candidate_id>
   ```
3. Give the rendered prompt to a clean adversarial agent. The prompt includes the frozen directive, primary sources, provenance slice, benchmark entries, and graph context.
4. If the adversarial verdict is NO, archive the report and run `$AA set-status <candidate_id> audited`.
5. If the adversarial verdict is YES, run `$AA set-status <candidate_id> defense_pending` and run Codex defense per parameter.
6. Pause for `/compact` between Phase C batches. This is context hygiene, not a user gate.

### Codex Defense

Defense sessions are keyed per parameter. The same parameter recurring across stages resumes the same Codex session.

```bash
$AA codex-defense <candidate_id> chi_Q 1
```

The wrapper captures `codex_session_id:` from Codex output and stores it at:

```yaml
codex_session:
  by_parameter:
    chi_Q:
      session_id:
      log_paths: []
```

The defense log must contain per-finding `DEFEND`, `CONCEDE`, or `PARTIAL` outcomes. A zero-exit Codex wrapper without parsed outcomes is recorded as `defense_pending` and fails informatively.

### Adjudication

Use `prompts/adjudication_consult.md` for the Claude+Codex consult record under `redteam_adversarial/verdicts/`. The signed-off outcomes are `FIND_STANDS`, `FIND_FAILS`, or `PARTIAL`. User escalation happens only for `FIND_STANDS` or conceptual changes, per `docs/adversarial_audit_execution_plan.md`.

### Dry Run

Run exactly:

```bash
$AA dry-run --stages 003 104 105
```

Expected behavior:

- Stage 003 yields zero candidates.
- Stages 104/105 yield at least one real candidate with stage and file-line anchors.
- The candidate advances through `scanned`, `provenance_built`, and `audit_pending`.
- A provenance slice is generated under `redteam_adversarial/provenance/`.
- A Phase C prompt is rendered under `redteam_adversarial/tmp_prompts/`.
- No binding verdict fields are populated.

Dry-run artifacts carry `dry_run: true` and `non_binding: true` where applicable, and are removable with `purge-dry-run`.

## Files Read

- `research/pde_ledger/.redteam-config.yaml`
- `docs/adversarial_audit_directive.md`
- `docs/adversarial_audit_deployment.md`
- `docs/adversarial_audit_execution_plan.md`
- `graph/fluid_universe_derivation_atlas_graph.yaml`
- `graph/query_graph.py`
- `research/pde_ledger/redteam/pass2/RECONCILIATION_AUGMENTATION.md`
- `research/pde_ledger/redteam/pass2/reports/stage_*.md`
- `research/pde_ledger/notes/CHECKPOINT_CONSTANT_PROVENANCE.md`
- Per-candidate paper, notes, script, and output sources.

## Files Written

All writes are under `research/pde_ledger/redteam_adversarial/` except the additive `adversarial:` config section in `research/pde_ledger/.redteam-config.yaml`.

- `MANIFEST.yaml`
- `fit_insertion_points.yaml`
- `benchmarks.yaml`
- `phase_a_fragments/*.yaml`
- `provenance/*.yaml`
- `reports/*.md`
- `defenses/*.md`
- `verdicts/*.md`
- `BATCHES.md`
- `tmp_prompts/*.md`
- `codex_logs/*.txt`

## Invariants

- Every manifest write goes through `flock` serialization in `lib/adversarial.sh`.
- Every script execution path uses `timeout 600`; exit 124 is failure and requires reformulation, never a higher cap.
- If any Wolfram `.wl` execution path is added later, keep at most two concurrent `math -script` processes.
- Runtime-rendered prompts are written under `redteam_adversarial/tmp_prompts/`; operating-system temp prompt files are forbidden.
- Agent-readable files are YAML or markdown. JSON is reserved for machine-to-machine subprocess output only.
- Phase A modalities are mutually blind. The union step records which modality found each candidate.
- Benchmarks come from `benchmarks.yaml`; Phase C must not adjudicate external matches from model memory.
- Defense Codex sessions are resumed per parameter, not per stage and not per candidate.
- Clean adversarial agent per candidate. Keep critique separate from construction.
- Auto mode follows `docs/adversarial_audit_execution_plan.md`: user gate only on `FIND_STANDS` or conceptual change.
