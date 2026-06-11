# Adversarial Audit Skill — Build Directive (Step 2)

**Status:** ACTIVE directive for the Step 2 build (execution plan `docs/adversarial_audit_execution_plan.md`, Step 2).
**Roles (binding):** Codex designs and writes ALL runtime code (scripts, SKILL.md, prompts, templates, config). Claude reviews against the acceptance criteria below. This directive states requirements and acceptance criteria only — internal design (file layout inside `lib/`, language choices, function structure, mechanical-script vs agent-prompt per scan modality) is Codex's call.

## Binding references (read before designing)

- `docs/adversarial_audit_deployment.md` — design doc. §3.1 (three phases, scan modalities, benchmark sourcing), §4 (pipeline shape), §5 (skill design sketch — the authoritative component list), §3.2 (graph backbone), §3.3 (provenance taxonomy).
- `docs/adversarial_audit_execution_plan.md` — operating contract §1–7 and Step 2 exit gate.
- `docs/adversarial_audit_directive.md` — the per-pass adversarial directive. **FROZEN. Consumed by reference/inclusion only. Never edited, never paraphrased in prompts.**
- `.claude/skills/redteam-audit/` — structural template (SKILL.md / `lib/` / `prompts/` / `templates/`, manifest mechanics, flock serialization, Codex invocation + session resume). Reuse its patterns; do NOT modify it.

## Path convention

All `redteam_adversarial/...` paths below are relative to `research/pde_ledger/`. The skill lives at `.claude/skills/adversarial-audit/` (project-root relative).

## Requirements

**R1 — Skill structure.** `.claude/skills/adversarial-audit/` mirroring the redteam-audit shape: `SKILL.md`, `lib/`, `prompts/`, `templates/`. SKILL.md documents commands, runbooks (at minimum: phase-A scan, phase-B build, phase-C batch, status/summary, dry-run), the state machine (R3), files read/written, and invariants — parallel in form to `redteam-audit/SKILL.md`.

**R2 — Config.** Either `.adversarial-audit-config.yaml` at project root or an `adversarial:` section added to the existing `.redteam-config.yaml` — Codex's choice. Whatever the choice, the redteam-audit skill's behavior must be unchanged. Config must carry at minimum: project name; adversarial artifact root (`redteam_adversarial`); paths to the frozen directive, the atlas graph YAML, `query_graph.py`, and the Phase A seeds (`research/pde_ledger/redteam/pass2/RECONCILIATION_AUGMENTATION.md` + its per-stage outputs, `research/pde_ledger/notes/CHECKPOINT_CONSTANT_PROVENANCE.md`); Codex wrapper path + sandbox mode; limits (`max_iterations`, parallelism caps — cap VALUES are decided later by the Step 3 consult, so they must be config-editable, not hardcoded).

**R3 — Manifest.** `redteam_adversarial/MANIFEST.yaml`. The tracked unit is the **fit-insertion-point candidate**, not the stage. Per-candidate state machine exactly as deployment doc §5:

```
pending / scanned / provenance_built / audit_pending / audited /
defense_pending / defended / adjudicating / verdict_logged / blocked
```

Per-candidate entry must include: candidate id; anchor stage(s); file:line citation(s); parameter name(s); batch id (nullable until the Step 3 batching consult); status + status timestamps; `codex_session` (see R6 — keyed per-parameter); paths to its report / defense / verdict artifacts as they are produced; verdict fields (adversarial YES/NO; adjudication FIND_STANDS / FIND_FAILS / PARTIAL). Every manifest write goes through flock serialization (the `_manifest_locked` pattern from `redteam-audit/lib/redteam.sh`) so parallel invocations are race-safe.

**R4 — Report tree.** Under `redteam_adversarial/`: `fit_insertion_points.yaml` (Phase A union output), `benchmarks.yaml`, `provenance/` (one YAML per parameter), `reports/`, `defenses/`, `verdicts/`, `BATCHES.md`, `tmp_prompts/`, `codex_logs/` (already exists). Every file an agent reads or writes is YAML or markdown; JSON only for pure machine-to-machine data.

**R5 — Graph wrappers.** CLI wrapper(s) giving Phase B/C tooling access to `graph/query_graph.py` (at minimum `neighborhood`, `paths`, `show`, `source`), graph path taken from config.

**R6 — Codex invocation.** Reuse the redteam-audit invocation pattern (wrapper invocation, `codex_session_id:` capture from output, `--resume` on subsequent iterations, per-iteration logs under `codex_logs/`), with one deliberate difference: defense sessions are keyed and resumed **per-parameter**, not per-stage or per-candidate — the same parameter recurring across stages reuses one session (deployment doc §7, resolved).

**R7 — Prompts.** Templates for:

- (a) **Phase A scan modalities** — the four from deployment doc §3.1 (numeric-literal, claim-label, graph, existing-provenance cross-check) plus the completeness critic. Each modality runs independently and must be blind to the other modalities' outputs; each emits a YAML candidate fragment; a union/merge step dedupes into `fit_insertion_points.yaml` while preserving per-modality attribution (which modality(ies) found each candidate). Whether a given modality is a mechanical script or an agent prompt is Codex's design call, except the completeness critic, which is agent-based by nature.
- (b) **Phase B provenance builder** — must embed both binding rules from deployment doc §2.6.2/§3.1: *no attribution without opening the per-stage `notes/` source*, and *a self-contradictory origin claim is itself a finding, never silently resolved*. Must also encode the `graph_gap` logging rule (§7) and the three-kind provenance taxonomy (the ledger is parameter-value provenance ONLY).
- (c) **Phase C adversarial wrapper** — assembles, per candidate: the frozen directive (by inclusion or path, unmodified), the candidate's primary source paths, its provenance slice, its `benchmarks.yaml` entry/entries, and graph context. External-match checks point at the sourced benchmark entry — never at model memory.
- (d) **Codex defense prompt** — reads the adversarial report + same sources + provenance; produces DEFEND / CONCEDE / PARTIAL per finding.
- (e) **Adjudication consult record template** — for the structured Claude+Codex consult records under `verdicts/` (same shape as the numbering-escalation consult records).

All runtime-rendered prompts are written under `redteam_adversarial/tmp_prompts/` — **never `/tmp`** (sub-agents cannot read `/tmp`; this supersedes any `/tmp` convention seen in redteam-audit's SKILL.md).

**R8 — Standing rules, encoded.** SKILL.md and any runners must encode: `timeout 600` hard cap on any script execution (exit 124 = failure → reformulate, never raise the cap); ≤2 concurrent `math -script` if any `.wl` execution path exists; clean agent per adversarial pass; the auto-mode contract pointer (user gate ONLY on FIND_STANDS or conceptual change); `/compact` between Phase C batches.

**R9 — Dry-run mode.** A mode exercising the A→B→C plumbing end to end on an explicit small stage list **without producing binding verdicts**. All dry-run artifacts must be unambiguously marked non-binding (mechanism is Codex's choice: e.g., a `dry_run: true` field, an isolated subtree, or both) and trivially removable without touching real campaign state.

## Dry-run specification (part of this build's exit gate)

Run the dry-run against stages **003, 104, 105**:

- 003 is pure scaffolding — expected to yield zero candidates (exercises the "structurally vacuous → no candidate" path honestly; do not force a candidate).
- 104/105 is the known χ_Q = 1 episode (fingerprint at 104, fix at 105; see commit `b052471`) — expected to yield at least one real candidate, walk it through provenance-slice build and an adversarial-prompt render, and advance its manifest entry through the state machine.

The dry-run must demonstrate: ≥1 candidate detected with stage + file:line anchor; a provenance slice generated for it; a Phase C prompt rendered for it under `tmp_prompts/`; manifest states transitioning correctly; zero binding verdicts recorded. Iterate until the dry-run exits 0 end to end.

## Out of scope (binding)

- No edits to: `docs/adversarial_audit_directive.md` (frozen), the deployment doc, the execution plan, this directive, anything under `.claude/skills/redteam-audit/`, `research/pde_ledger/paper/`, `research/pde_ledger/notes/`, stage scripts, `graph/*.yaml`, or `redteam/` (both passes' artifacts are closed records).
- No real Phase A sweep across the 253 stages; no binding verdicts of any kind; no entries in `benchmarks.yaml` beyond dry-run-marked placeholders.
- No hand-authored candidate/provenance/benchmark content presented as generated output.

## Acceptance criteria (Claude review checks each)

1. Skill tree exists at `.claude/skills/adversarial-audit/` with all four components; SKILL.md states the R3 state machine verbatim.
2. Config loads; the redteam-audit skill is byte-identical (or, if the shared-config option was chosen, `.redteam-config.yaml` changed only additively and `redteam-audit/lib/redteam.sh summary` still runs clean against pass-2).
3. Manifest init produces valid YAML with the R3 per-candidate schema; **every** writer path goes through flock serialization.
4. All four scan-modality templates + critic template exist; modalities are mutually blind; union step preserves per-modality attribution.
5. Phase B prompt embeds both binding rules and the `graph_gap` rule; states the parameter-value-only taxonomy.
6. Phase C wrapper includes the frozen directive unmodified and all four context slices; benchmark adjudication references `benchmarks.yaml`, never model recall.
7. Defense invocation resumes Codex sessions per-parameter; session ids captured into the manifest.
8. Runtime prompts render under `redteam_adversarial/tmp_prompts/`; no `/tmp` writes anywhere in the skill.
9. Dry-run against 003/104/105 meets every bullet of the dry-run specification, including the zero-candidate result on 003.
10. `timeout 600` on every script-execution path; no manifest write bypasses the lock.
11. `git status` shows changes ONLY in: `.claude/skills/adversarial-audit/`, `research/pde_ledger/redteam_adversarial/`, and (only if the shared-config option was chosen) `.redteam-config.yaml`.
12. Dry-run artifacts are unambiguously marked non-binding and removable in one step.

## Reporting contract (Codex)

When done, report: (1) every file created/changed with a one-line purpose; (2) the design choices made where this directive left them open (config placement, modality mechanization, dry-run marking mechanism) and why; (3) dry-run transcript evidence for each dry-run-specification bullet; (4) any acceptance criterion not satisfied and why. Do not commit — leave the tree dirty for review.
