---
name: redteam-audit
description: Adversarial red-team audit of a research paper's SymPy and Mathematica scripts. Scripts-only — does NOT touch paper.tex or notes/. Audits one batch at a time (math errors cascade) and halts on any stop-cold or blocking finding, invokes Codex (resuming a per-unit session across iterations) to apply fixes, and verifies. Reusable across papers via .redteam-config.yaml.
allowed-tools: Bash, Read, Edit, Write, Agent
user_invocable: true
---

# Red-Team Script Audit

Adversarial, chunked verification pass for a research paper's math, with the audit unit being **a (SymPy + Mathematica) script pair**. The red-team checks whether each pair is mathematically valid, non-tautological, and that the two engines independently agree.

Doc alignment (paper.tex, notes/) is out of scope and handled manually after script verification lands.

Operates from `.redteam-config.yaml` in the current working directory, so the same skill drives any paper that follows the same per-unit script layout.

## Core principles

1. **Scripts only.** The red-team reads and modifies only `.py` and `.wl` scripts (and their saved `.txt` outputs). It does NOT read paper.tex or notes/. Doc alignment is a manual step done by the user after the scripts are verified.
2. **Sequential batches.** Audit one batch and finish it before starting the next: **math errors cascade — finishing batch 1 may invalidate batch 2's premise.** ⛔ **Never roll forward past a stop-cold verdict or an unresolved blocking finding.** ⚠ The standing *per-chunk user gate* was cut 2026-07-29/30 (`feedback_physics_not_ceremony`): halt for the user at a decision, a blocking finding or a no-go — not at every batch boundary.
3. **Two-AI loop.** Claude (orchestrator + audit/verify agents) and Codex (fix application) communicate via structured markdown directives in `redteam/directives/`. No chat handoff required.
4. **Per-unit codex session.** Each audit unit's fix loop uses a single codex session id, resumed across iterations so codex retains the context of what it's already attempted. The session is stored in the manifest and cleared automatically when upstream changes invalidate the unit's context.
5. **Clean-context audit agents.** Each unit gets a fresh sub-agent for audit and a fresh sub-agent for verify. No cross-unit contamination on the auditor side. (Codex is the opposite — it benefits from resumed context.)
6. **Both engines required.** A unit is `verified` only when SymPy AND Mathematica each independently check the load-bearing claims. Missing scripts are findings, not silent gaps.
7. **Stop-cold verdicts.** `UNFIXABLE` and `CRITICAL_DOWNSTREAM` halt the loop immediately. The framework never pushes past either.
8. **YAML/markdown only.** No JSON in any file an LLM reads or writes.

## Finding categories (9, scripts-only)

The auditor classifies every finding into one of these:

1. `tautological_check` — assertion guaranteed by construction (`x = expr; assert x == expr`)
2. `hardcoded_result` — literal numeric / symbolic answer used with no derivation
3. `symbol_assumption_error` — wrong domain, missing assumption, invalid `simplify` step
4. `missing_branch` — only one branch tested when claim spans a range
5. `engine_disagreement` — sympy and mathematica produce mismatched results
6. `mathematica_transliteration` — `.wl` is a line-by-line port of `.py` rather than independent re-derivation
7. `missing_verification_script` — `missing_sympy`, `missing_mathematica`, or `script_doesnt_cover_claim`
8. `insufficient_verification` — assertions exist but don't actually exercise the claim
9. `stale_output` — saved output `.txt` mtime older than the script's mtime

Plus stop-cold verdicts: `UNFIXABLE`, `CRITICAL_DOWNSTREAM`.

## State machine (per unit)

```
pending
  ↓ audit
audited ─────────────── (no findings) ────────→ verified
  ↓ findings exist
directive_ready
  ↓ codex (fresh session OR --resume existing session)
codex_applied
  ↓ run scripts + verifier agent
verified | needs_rework | blocked_unfixable | blocked_critical_downstream
  ↓ if needs_rework AND iter < max: back to codex (--resume same session, delta directive)
upstream_stale ── (re-trigger audit on this unit; clear codex session — old context is stale)
```

## Commands

All commands except `bootstrap`, `detect`, and `help` require `.redteam-config.yaml` in the current working directory (or any parent).

Setup & inspection: `bootstrap`, `detect`, `init`, `status`, `scan`, `ls`, `batch-info`, `stage-info`, `paths`, `next-batch`, `blocked`, `render-batches`.

Mutation primitives: `set-status`, `set-iter`, `mark-stale-downstream`.

Execution primitives: `exec-sympy`, `exec-mathematica`, `capture-diff`, `codex-invoke`, `codex-reset`.

Prompt rendering: `render-audit-prompt`, `render-verify-prompt`.

Run `redteam.sh help` for the full reference. Run `redteam.sh bootstrap` to set up a new project.

## Runbooks — what the orchestrator does when invoked

When the user invokes this skill as `/redteam-audit <subcommand> <args>`, follow the matching runbook. The orchestrator is Claude; the runbooks below describe Claude's exact actions step by step.

Common preamble for **every** runbook:
- `cd` into the project root (the directory containing `.redteam-config.yaml`).
- Define `RT=/var/projects/toy_physics/.claude/skills/redteam-audit/lib/redteam.sh` for brevity.
- All status checks and mutations go through `$RT` subcommands — never edit `MANIFEST.yaml` directly.

### `audit <unit>` or `audit batch:<ID>`

**Single unit:**

1. Check current status: `$RT stage-info <unit>`. Audit is valid from `pending`, `audited`, or `upstream_stale`. Refuse otherwise with an explanation.
2. Render the audit prompt:
   ```bash
   $RT render-audit-prompt <unit> > /tmp/audit_prompt_<unit>.md
   ```
3. Spawn ONE Agent call (subagent_type: `general-purpose`):
   - description: `"Audit unit <unit>"`
   - prompt: contents of `/tmp/audit_prompt_<unit>.md` (read it and inline the contents — the Agent does not have access to the temp file unless explicitly told)
4. After the Agent returns, read `redteam/reports/stage_<unit>.md` to get its YAML frontmatter.
5. Transition based on verdict:
   - `verdict: clean` → `$RT set-status <unit> verified`, then set audit metadata via stdin-to-set-stage-fields:
     ```bash
     printf 'last_audit_date: "%s"\nlast_audit_findings: 0\n' "$(date -Iseconds)" \
       | python3 .../manifest.py .redteam-config.yaml set-stage-fields <unit>
     ```
   - `verdict: findings` and `stop_cold: null` → `$RT set-status <unit> directive_ready`; record `last_audit_findings: N`.
   - `stop_cold: UNFIXABLE` → `$RT set-status <unit> blocked_unfixable`; **halt the entire pipeline**; report to user with link to the report.
   - `stop_cold: CRITICAL_DOWNSTREAM` → `$RT set-status <unit> blocked_critical_downstream`; **halt**; report.

**Batch:**

1. `$RT batch-info <ID>` to get the unit list.
2. Use `$RT status --batch <ID>` to find which units are pending/auditable.
3. Read `limits.parallel_audit_max` from `.redteam-config.yaml`.
4. Spawn up to `parallel_audit_max` Agent calls **in a single message** (parallel). Each Agent gets one unit's rendered prompt.
5. After all return, ingest each report (steps 4–5 from the single-unit runbook).
6. If ANY stop-cold verdict appeared, halt and report.
7. Otherwise summarize: counts of `verified` (clean) / `directive_ready` / blocked. Recommend the next subcommand (`fix batch:<ID>` if any are `directive_ready`).

### `fix <unit>` or `fix batch:<ID>`

**Single unit (must be in `directive_ready` or `needs_rework`):**

1. Check status. Refuse if not `directive_ready` or `needs_rework`.
2. Get current `iteration_count` from `$RT stage-info <unit>`. Next iter = current + 1.
3. Determine the directive path:
   - On `directive_ready`: `redteam/directives/stage_<unit>.md` (the original).
   - On `needs_rework`: the last verification report (`redteam/verifications/stage_<unit>.md`) contains a `## Delta directive` section. Extract that section and write it to a temp delta file, then pass that.
4. `$RT set-status <unit> fixing`.
5. Invoke codex:
   ```bash
   $RT codex-invoke <unit> <directive-or-delta-path> <next-iter>
   ```
   - First iteration with no `codex_session`: codex.md preamble is prepended; new session is created and stored.
   - Later iterations: `codex-chat --resume <session>` is used (preamble not re-sent; codex remembers it).
6. After codex returns, read the directive file. Confirm it now has `applied: true` in front-matter, or count the `## Applied:` / `## Blocked:` blocks.
7. `$RT set-status <unit> codex_applied`.

**Batch:**

1. List units in `directive_ready` or `needs_rework`, sorted by stage number ascending.
2. For each, run the single-unit fix runbook **sequentially**. Later fixes may depend on earlier ones — do NOT parallelize.

### `verify <unit>` or `verify batch:<ID>`

**Single unit (must be in `codex_applied`):**

1. Check status. Refuse if not `codex_applied`.
2. `$RT set-status <unit> verifying`.
3. Execute the scripts (orchestrator runs them, not the agent):
   ```bash
   $RT exec-sympy <unit>          # captures redteam/exec_logs/stage_<unit>_sympy.log
   $RT exec-mathematica <unit>    # captures redteam/exec_logs/stage_<unit>_mathematica.log
   ```
4. Capture diff:
   ```bash
   $RT capture-diff <unit>        # writes redteam/exec_logs/stage_<unit>_diff.patch
   ```
5. Render the verifier prompt:
   ```bash
   $RT render-verify-prompt <unit> > /tmp/verify_prompt_<unit>.md
   ```
6. Spawn ONE verifier Agent (clean context):
   - description: `"Verify unit <unit>"`
   - prompt: contents of the rendered file.
7. After Agent returns, read `redteam/verifications/stage_<unit>.md` frontmatter.
8. Transition based on `verdict`:
   - `verified` → `$RT set-status <unit> verified`; record `last_verify_date`.
     If `material_change: true`:
     - `$RT mark-stale-downstream <unit>` — this BOTH demotes every non-pending, non-blocked downstream unit's status to `upstream_stale` AND clears each demoted unit's `codex_session` (its old session is about a now-stale premise). Pending units stay pending; blocked units stay blocked.
   - `needs_rework`:
     - If `iteration_count < max_iterations`: leave status `needs_rework`. The `run` subcommand will loop back to fix. If invoked stand-alone, report to user and stop.
     - If `iteration_count >= max_iterations`: `$RT set-status <unit> blocked_unfixable`; halt; report.
   - `blocked_unfixable` → `$RT set-status <unit> blocked_unfixable`; halt; report.

**Batch:**

1. List units in `codex_applied`.
2. Run the exec + capture-diff phase **sequentially** per unit (cheap; each takes seconds to minutes).
3. Spawn verifier Agents in parallel (up to `parallel_verify_max`) — one per unit.
4. Ingest each verification report; apply transitions per unit.
5. If ANY `blocked_unfixable` or unmet max-iter, halt; report.

### `run <unit>` or `run batch:<ID>`

End-to-end audit → fix → verify with iteration up to `limits.max_iterations`:

1. Run `audit` per the runbook above.
2. If any stop-cold, halt.
3. For all units now in `directive_ready`: run `fix`. (Sequential.)
4. For all units now in `codex_applied`: run `verify`. (Parallel where allowed.)
5. For all units now in `needs_rework` AND `iteration_count < max_iterations`: loop back to step 3 (fix → verify).
6. After the loop, **report at the batch boundary** (and halt there only on a stop-cold or unresolved blocking finding). Print:
   - per-unit verdict table
   - link to `redteam/batches/batch_<ID>.md` (if you wrote one)
   - the next-batch suggestion
   - if anything is blocked: "Review the diffs and verifier reports, then invoke `/redteam-audit run batch:<NEXT>` to proceed."

**Never auto-advance past a stop-cold verdict or an unresolved blocking finding.** ⚠ A clean batch does not need user approval to proceed (per-chunk user gate cut 2026-07-29/30); report and continue.

### Iteration cap & halt conditions

The orchestrator halts and reports when ANY of these happen:
- Auditor verdict is `UNFIXABLE` or `CRITICAL_DOWNSTREAM`
- Verifier verdict is `blocked_unfixable`
- A unit hits `needs_rework` at `iteration_count >= max_iterations` (set status to `blocked_unfixable`)
- A batch completes with any unit `blocked_*` — halt and report before the next batch (⚠ a batch in which every unit is `verified` does **not** halt; per-chunk user gate cut 2026-07-29/30)

### Output expectations to the user

After every runbook completes (or halts):
- Print a concise summary: units processed, verdicts, anything blocked
- Link to the relevant report/verification files by path
- State the next recommended subcommand (or "review and confirm")

Avoid burying important verdicts. Stop-cold and blocked units are the top of the report.

## Files this skill reads/writes

Reads from project root (per config):
- `.redteam-config.yaml`
- `scripts/*.py`, `mathematica/*.wl`, `scripts/output/*.txt`, `mathematica/output/*.txt`

Writes under `redteam/`:
- `MANIFEST.yaml` — single source of truth (schema v2)
- `BATCHES.md` — human-readable status overview (auto-regenerated)
- `reports/stage_<NNN>.md` — auditor findings
- `directives/stage_<NNN>.md` — patch contract for Codex (with appended `## Applied`/`## Blocked` blocks after codex runs)
- `verifications/stage_<NNN>.md` — post-codex verifier report
- `batches/batch_<NN>.md` — batch-level summary
- `exec_logs/stage_<NNN>_{sympy,mathematica}.log` — captured script execution output
- `exec_logs/stage_<NNN>_diff.patch` — git diff captured for the verifier
- `codex_logs/<NNN>_iter<N>.txt` — captured codex session output (full transcript per iteration)

## Manifest schema (v2)

Paths are stored **relative to the project root** (the directory containing `.redteam-config.yaml`) so the manifest is portable across checkout locations. CLI helpers resolve to absolute when invoking scripts or building agent prompts. Pre-v2 manifests stored absolute paths; legacy absolute paths still work but new entries are relative.

Each unit entry looks like:

```yaml
'058':
  stage: 58
  stage_str: '058'
  batch_id: III.2
  status: pending
  iteration_count: 0
  upstream_stale: false
  is_checkpoint: false
  is_status_only_candidate: false
  last_audit_date: null
  last_audit_findings: 0
  last_verify_date: null
  # Populated by codex-invoke after first fix:
  # codex_session: <uuid>
  # codex_log_paths: [redteam/codex_logs/058_iter1.txt, ...]
  # last_codex_run: <iso8601>
  # last_codex_exit: 0
  files:
    sympy:              {path: scripts/..._sympy_audit.py, exists: true, mtime: ...}
    mathematica:        {path: mathematica/..._mathematica_audit.wl, exists: true, mtime: ...}
    sympy_output:       {path: scripts/output/..._sympy_audit.txt, exists: true, mtime: ...}
    mathematica_output: {path: mathematica/output/..._mathematica_audit.txt, exists: true, mtime: ...}
```

v1 had `files.tex`, `files.notes`, and absolute paths; v2 drops the doc roles and uses relative paths. `init` produces v2 directly.

## Invariants the CLI enforces

- Manifest writes are atomic (write-temp-then-rename).
- When an upstream unit's verify is `material_change: true`, every non-pending, non-blocked downstream unit is demoted to `status: upstream_stale` (not just flagged) and its `codex_session` is cleared.
- A unit cannot enter `fix` unless in `directive_ready` or `needs_rework`.
- `fix` runs sequentially per unit within a batch.
- `next-batch` returns the earliest batch with any audit-eligible unit (`pending` or `upstream_stale`); cascade-affected batches surface before later ones.
- A batch cannot advance to next-batch state until every unit in it is `verified` or `blocked_*`.

## Codex invocation

Uses `~/.claude/hooks/codex-chat/codex-chat`. The wrapper supports `--resume <session-id>` and auto-extracts the session id from codex's output (`codex_session_id: <uuid>` appended at the end). The skill stores this id in `MANIFEST.yaml` at `stages.<NNN>.codex_session` and passes it back on the next iteration. The directive is piped to codex via stdin.

For a new session, `codex-invoke` prepends the `prompts/codex.md` system prompt to the directive content. For a resumed session, only the directive is sent (codex already has the system prompt in its session context).

Never sets a timeout on the Bash call (codex runs can take minutes).

Session lifecycle:
- First fix iteration on a unit → new session (preamble + directive)
- Iterations 2..N within the same audit pass → `--resume` (directive only)
- Unit went `verified` → `upstream_stale` → `codex_session` cleared automatically
- `--resume` returns "session not found" → fall back to new session, log a warning

## Bootstrapping a new paper

When the user asks to bootstrap a new research paper, drive the workflow rather than running `bootstrap` blindly:

1. `cd` to the project root and run `$RT detect .` to inspect what the scanner found.
2. Show the detected script roles, stage range, and any candidate parts/chapters dirs.
3. Use AskUserQuestion to confirm:
   - Project name
   - Whether to use detected parts/chapters dirs to inform batch boundaries (vs uniform N-stage batches)
   - Default batch size
   - Whether any non-default runners are needed (e.g., `sage` instead of `python3`)
4. Run `$RT bootstrap . --batch-size N [--prefer-role role=glob ...]`
5. Read the generated config and use Edit to add the user-supplied checkpoint stages, status-only candidates, and batch labels.
6. Run `$RT init` to seed the manifest.
7. Show the resulting `BATCHES.md` and the first batch's `status --batch <ID>`.

## Scope notes

- `pde_ledger` is the reference project. 253 stages across 21 batches. Configured and ready.
- `pde_audit` is stage-decomposed; could be onboarded (has notes but no paper.tex — the red-team is scripts-only, so this doesn't matter).
- Future `pde_*` papers will follow the same stage-decomposed pattern.
- **Older 4d_* papers are intentionally out of scope.** They use monolithic paper.tex + topic-keyed scripts and the user decided manual paper-comparison after script verification is the simpler workflow there. The red-team will not be retrofitted for them.

## Memory-respected rules

- Audit agents must NOT run `python3 -c` commentary scripts (per user feedback). They read and reason. Real script execution happens only in the orchestrator's `verify` phase, against the actual `sympy_audit.py` / `mathematica_audit.wl` files.
- Each unit audit and verify gets a clean agent (per user feedback). Codex is the exception — it keeps a resumable session per unit so it remembers what it just tried.
- Scripts are checked specifically for hardcoded values, tautological assertions, and symbol assumption errors (per user feedback).
- Batches are audited sequentially because math errors cascade; never roll forward past a stop-cold verdict (per user feedback). ⚠ The per-chunk user gate between clean batches was cut 2026-07-29/30 (`feedback_physics_not_ceremony`).
- All LLM-readable state is YAML or markdown; no JSON (per user feedback).
- Script-correctness is the primary check; paper/notes alignment is out of scope here — done manually after the red-team verifies the math (per user feedback).

## Implementation files

- `lib/redteam.sh` — CLI entry point. All subcommands dispatched here.
- `lib/_yaml.py` — YAML read/write helper (PyYAML-based; called from bash).
- `lib/scan.py` — filesystem scanner for one unit's script artifacts.
- `lib/manifest.py` — fast single-pass manifest queries and mutations.
- `lib/detect.py` — project-tree detector for bootstrap.
- `lib/bootstrap.py` — config generator from detection report.
- `lib/render_prompt.py` — placeholder substitution for auditor/verifier prompts.
- `prompts/{auditor,codex,verifier}.md` — locked agent prompts (scripts-only).
- `templates/{config.yaml,report.md,directive.md,verification.md,batch_summary.md}` — markdown / yaml skeletons.
