# Red-team pipeline — design and thoroughness

A description of the adversarial verification pipeline applied to the
SymPy + Mathematica scripts that back the moving-throat / PDE-ledger
research line.

Status as of 2026-05-27: 96 / 253 stages verified at v2 depth (~38%)
through 7 batches (I.1, I.2, II.1, III.1, III.2, III.3, batch III.4 in
progress). Cumulative findings closed: ~299.

## 1. Goal

For every research stage, prove or disprove that the stage's two
verification scripts — one SymPy (`scripts/..._sympy_audit.py`), one
Mathematica (`mathematica/..._mathematica_audit.wl`) — **actually
verify** the claim the paper states the stage proves. "Verify" here
means three things at once:

1. The two engines independently derive the claim from the physical
   premises, not from each other.
2. The bottom-line assertions are non-tautological and exercise the
   stage's load-bearing identities.
3. The identities, constants, and parameter ranges in the scripts
   match the identities, constants, and parameter ranges in the paper
   card and source notes.

A script that exits 0 but checks a different identity than the paper
requires is a defect. A script that checks the right identity but does
it tautologically is a defect. A script whose constants disagree with
the paper card's constants is a defect — direction of resolution
(update script or update paper) is the user's call, not the system's.

## 2. The unit of work

An **audit unit** is a single research stage and its (SymPy +
Mathematica) script pair. The unit is the audit boundary, the codex
session boundary, and the verifier's grain. Cross-stage reasoning is
explicitly disallowed at audit time — each audit agent reads only the
unit's paper card, source notes, part appendix, and the two scripts.
This isolation is enforced by:

- The `auditor.md` reading list is closed and explicit (paper card,
  notes for this stage, part appendix, the two scripts, two saved
  outputs). The trackers in `notes/` are out of scope.
- Each audit and each verify get a fresh clean-context sub-agent. No
  cross-stage agent re-use.
- Codex's session is per-unit and is cleared when upstream changes
  invalidate the unit's context.

## 3. Roles

Three actors collaborate, each restricted to a competence boundary:

- **Claude (orchestrator + audit + verify):** routes work, parses
  YAML front-matter, spawns sub-agents, decides batch transitions,
  applies orchestrator hot-fixes when a Codex iteration produces a
  flawed but locally-correctable result.
- **Codex (math authority):** receives a structured markdown
  directive, applies it to the scripts, iterates until `python3` and
  `math -script` exit 0. Codex *never* edits `paper/`, `notes/`, or
  other prose documents.
- **User (conceptual authority):** receives `paper_misalignment`
  questions through a structured resolutions doc, picks a direction
  per finding (paper-side fix vs script-side fix vs both-from-third-
  source review). The pipeline halts until the user has chosen.

The split is the load-bearing piece: math correctness questions go to
Codex because Codex can re-derive; convention/direction questions go
to the user because no engine can decide between two consistent
choices.

## 4. The 10 finding categories

Every finding is classified into exactly one of the following. The
CLI parses category names verbatim from the auditor's YAML output.

1. `tautological_check` — assertion guaranteed by construction. E.g.,
   `x = expr; assert x == expr`. Cannot fail no matter what the
   physics is.
2. `hardcoded_result` — literal numeric or pre-baked symbolic value
   used as the answer with no derivation and no provenance comment.
3. `symbol_assumption_error` — wrong domain (real vs complex,
   positive vs unrestricted), or `simplify` under aggressive
   assumptions that hides a branch.
4. `missing_branch` — claim spans a parameter range, only one branch
   is exercised.
5. `engine_disagreement` — SymPy and Mathematica produce mismatched
   results on the same claim.
6. `mathematica_transliteration` — `.wl` is a line-by-line port of
   `.py` rather than an independent re-derivation. Violates the
   second-engine policy (both engines must derive the result
   independently from the physical premises, not echo each other's
   algebra).
7. `missing_verification_script` — one engine absent or substantially
   absent. Subtypes: `missing_sympy`, `missing_mathematica`,
   `script_doesnt_cover_claim`.
8. `insufficient_verification` — assertions exist but only check a
   limiting case, a single ratio, or two expressions both derived
   from the same source.
9. `stale_output` — saved `.txt` output mtime older than the script's
   mtime.
10. `paper_misalignment` — script's verified claim doesn't match the
    paper's stated claim, even if the script is internally
    consistent. Five subtypes:
    - `target_mismatch` — script's load-bearing assertion verifies a
      different identity than the paper's `\stagefield{Output}`
      requires.
    - `value_mismatch` — script uses a numeric constant whose value
      differs from what the paper or notes state.
    - `script_missing_paper_claim` — paper says the stage proves both
      X and Y; script only tests X.
    - `paper_missing_script_claim` — script tests Y, but the paper
      and notes don't mention Y.
    - `notes_contradicts_script` — source notes (which informed the
      .tex) state an intent the script does not honor.

Two **stop-cold verdicts** halt the entire pipeline:

- `UNFIXABLE` — math cannot be reconciled within the unit's scope
  even if Codex had unlimited edits. Example: the script's premise is
  internally inconsistent.
- `CRITICAL_DOWNSTREAM` — a finding's correction would mathematically
  propagate to later stages. Halts before invoking Codex; the user
  reviews the cascade.

`paper_misalignment` is **never** auto-`CRITICAL_DOWNSTREAM` — it
requires user resolution first (which sets the direction), then the
downstream-impact analysis happens.

## 5. The two-pass design (v1 then v2)

Each stage is audited twice:

- **v1 (scripts-only):** auditor reads only the scripts and saved
  outputs. Catches `tautological_check`, `hardcoded_result`,
  `mathematica_transliteration`, `insufficient_verification`, and the
  other purely-script categories.
- **v2 (paper-grounded):** auditor reads paper card + notes +
  appendix **first**, builds a mental model of what the stage claims,
  *then* reads the scripts. Adds the `paper_misalignment` category
  on top of the v1 set.

The v2 pass routinely catches defects that v1 missed:

- A "fix" that v1 marked as discharged but didn't actually exercise
  the underlying physics (most common pattern: a "round-trip" check
  that algebraically reduces to `1 = 1`).
- Constants in the script that disagree with the paper card by
  factors or by tens of percent, where neither side asserts the value
  numerically so v1 saw two consistent integers and moved on.
- Stale banner / docstring strings left over from a global stage
  renumber (the script's math is correct but the script's self-label
  is wrong).

Empirically, v2 finds new substantive defects on roughly 1-2 of every
12 stages re-audited, plus banner-relabel items.

## 6. The fix loop

When findings exist, the auditor writes a structured directive
(`redteam/directives/stage_<NNN>.md`) with one block per finding.
For non-`paper_misalignment` items the directive is a patch contract:

- `Target:` file:line
- `Issue:` one paragraph
- `Required change:` step-by-step edit instructions
- `Verification:` how the next verifier confirms the fix landed

For `paper_misalignment` items the directive's body is a
`## Resolve before fix_loop` block — a question routed to the user.
The orchestrator refuses to invoke Codex on the unit until the user
has chosen a direction.

The orchestrator then:

1. **User gate.** For every `paper_misalignment`, build a YAML
   resolutions doc, get a Codex recommendation (math authority), and
   present the user a single multiple-choice question per finding.
   The user picks (a), (b), (c), or "other".
2. **Codex apply.** For each unit with non-`paper_misalignment`
   findings, invoke `$RT codex-invoke <NNN> <directive> <iter>`.
   Codex resumes the unit's session and applies the directive.
3. **Auto-iteration.** Codex iterates the scripts until both `python3`
   and `math -script` exit 0. The system trusts Codex to keep
   pumping; the verifier reviews substance afterward, not whether the
   scripts run.
4. **Orchestrator hot-fix.** If Codex's applied edit is locally
   broken but the fix is a canonical idiom from the pitfall library
   (`codex.md`), the orchestrator may apply the edit directly rather
   than spending another Codex iteration on a known pattern.

## 7. The verifier

After Codex applies, a fresh sub-agent (the verifier) reads:

- The auditor's report and directive.
- The applied diff (`redteam/exec_logs/stage_<NNN>_diff.patch`).
- Codex's apply log.
- The current scripts and saved outputs.

The verifier's job is **destination verification, not prescription
verification**. It checks whether what landed in the file actually
exercises the paper's claim. It does not check whether Codex followed
the directive's literal instructions. This split matters: a few times
per batch Codex notices a math error in the directive and applies a
corrected fix; the verifier confirms the math, not the directive
adherence.

The verifier produces a verdict (`verified`, `needs_rework`,
`blocked_unfixable`, `blocked_critical_downstream`). On `needs_rework`
the auditor writes a delta directive and Codex resumes its session.

## 8. State machine (per unit)

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
  ↓ if needs_rework AND iter < max: back to codex (--resume, delta directive)
upstream_stale ── (re-trigger audit on this unit; clear codex session)
```

The `upstream_stale` transition fires when a later batch's findings
mutate a result an earlier stage relied on. The system never silently
rolls a downstream verification past an upstream change.

## 9. Batch boundaries

Stages are grouped into batches of 12 (with some variation at part
boundaries). Each batch is processed in five phases:

1. **Audit wave** — 12 parallel auditor sub-agents, one per stage.
2. **Triage** — orchestrator partitions findings into
   `paper_misalignment` (user gate) and script-side (Codex fix).
3. **User gate** — resolutions doc with Codex recommendations, single
   batched question to the user.
4. **Fix loop** — sequential Codex invocations per dirty stage.
   Orchestrator hot-fixes when warranted.
5. **Verify wave** — parallel verifier sub-agents on every dirty
   stage. Batch closes when all dirty stages return `verified`.

**Sequential audit chunks** is a hard rule: the system never rolls
forward to the next batch automatically. The user must explicitly
approve each batch transition. The reason is mathematical: a fix in
batch N may invalidate a premise batch N+1 relied on, and that
invalidation can be invisible at audit time. Halting between batches
forces a human review of the cumulative state.

## 10. Codex pitfall library

`codex.md` is a growing list of cross-engine pitfalls that the
auditor and Codex both consult. Each pitfall has a name, a symptom,
and a canonical defense. Currently 8 documented, 1 candidate:

- **#1** Trivial `simplify` under aggressive assumptions. Defense:
  audit `$Assumptions` and `assume(...)` lists; never widen
  assumptions to make a `simplify` pass.
- **#2** Symmetry / parity traps. Defense: identify integrand parity
  before asserting an integral vanishes.
- **#3** Variable-independence trap. Defense: for every
  `sp.diff(EXPR, VAR)`, list the symbols `EXPR` actually depends on
  and confirm `VAR` is among them. Otherwise the derivative is
  identically zero and the `assert_nonzero` will spuriously fail.
- **#4** Mathematica `ConditionalExpression` wrapper. Defense:
  preemptively strip `ConditionalExpression[e_, _] :> e` from
  `Solve`/`Reduce` outputs adjacent to downstream substitutions.
- **#5** Mathematica `Limit` is non-deterministic at poles. Defense:
  use `1/pi1 == 0` for pole checks rather than `=!= Infinity`.
- **#6** SymPy `simplify` may not zero-test transcendental
  combinations; use `trigsimp` and `expand_trig` first.
- **#7** Numerical-tolerance asserts should compare residual to a
  named threshold, not to a hard-coded literal.
- **#8** Heavy BVP `dsolve`/`DSolve` calls. Defense: replace with
  symbolic-integral or numerical-sweep equivalents; check for
  pre-existing `Integrate` to mirror.
- **#9 (candidate)** Mathematica `Integrate[]` does not factor
  constant multipliers when the integrand contains unspecified
  symbolic functions. Defense: verify integrand equality first, then
  assign the closed-form integral once.

Each pitfall is promoted to the canonical list only after at least
one production occurrence forces the fix.

## 11. Prose trackers

Six markdown trackers under `notes/` carry process state across
batches:

- `MATHEMATICA_MIRROR_POLICY.md` — per-stage Independent-Mirror Set
  entries; documents which `.wl` files have been certified as
  independent re-derivations rather than transliterations.
- `CHECKPOINT_TRUST_AUDIT.md` — per-checkpoint snapshot tier;
  tracks which stage-checkpoint outputs downstream batches are
  permitted to trust without re-deriving.
- `CHECKPOINT_CONSTANT_PROVENANCE.md` — every literal constant carried
  forward from a checkpoint, with its derivation provenance.
- `PAPER_CLEANUP_TRACKER.md` — per-batch summary row, plus a
  cumulative change-log of paper-side edits required by
  `paper_misalignment` resolutions.
- `EM_PROJECTED_INTEGRATION_TRACKER.md` — alignment of the projected-
  EM core (stages 004-021) against the rest of the ledger; updated
  every batch even if the batch is out of range, so the linear core
  has a continuous trust record.
- `STAGE_VERIFICATION_COVERAGE.md` — the single source-of-truth
  coverage table: per-batch row of stages verified, findings closed,
  dates, and material-change flag.

Every batch close updates all six trackers, even if a tracker has no
new substance to report. The absence of an update is itself a defect.

## 12. Constraints that prevent specific failure modes

Some failure modes are common enough that the system encodes them as
hard constraints rather than relying on agent judgment:

- **No JSON for LLM I/O.** Every file an LLM reads or writes is YAML
  or markdown. JSON is reserved for pure machine-to-machine data
  (e.g., the `MANIFEST.yaml` cache is YAML; raw codex apply logs are
  text). Reason: JSON's punctuation overhead inflates context for no
  benefit and trips agents on unicode-quote rendering.
- **No fake scripts.** Agents may not run `python3 -c '...'` or
  `math -script -e '...'` commentary scripts during an audit. They
  read and reason from the file contents. Reason: a CAS commentary
  call can hide a real error behind a passing one-liner.
- **Mathematica single-seat.** Only one `math -script` invocation is
  allowed across the whole system at a time. Reason: the licensed
  Mathematica kernel is single-seat; concurrent runs lock-out one
  another.
- **No parallel `exec-sympy` / `exec-mathematica`.** These commands
  write to `MANIFEST.yaml`; parallel calls race and corrupt it.
  Refresh outputs via direct `python3` / `math -script` instead.
- **Audit prompts under project root.** Sub-agents cannot read
  `/tmp/`; audit and verify prompts are rendered to
  `redteam/tmp_prompts/` under the project root.
- **Codex iterates until clean.** Codex runs its applied scripts
  until they exit 0; the verifier reviews substance afterward, not
  whether the scripts run.

## 13. What the pipeline catches

Categorized by frequency in the 7 batches processed to date:

- **High frequency** (every batch):
  - Tautological checks that look substantive (a "round-trip" that's
    actually `expr == expr`).
  - Hardcoded numeric literals copied byte-for-byte from the other
    engine's saved output.
  - Mathematica scripts that are line-by-line ports of the SymPy
    script.
  - Stale `.txt` output transcripts older than the script.

- **Medium frequency** (most batches):
  - `paper_misalignment` from a stale paper-side prose line vs an
    updated script (or vice versa).
  - Stale banner / docstring strings from the global stage renumber.
  - `simplify` calls under assumptions that mask a branch.
  - Numerical-tolerance asserts compared against a hardcoded literal
    rather than a named threshold.

- **Low frequency** (one or two batches each):
  - Material algebra errors that v1 missed and v2 caught by reading
    the paper first (e.g., a sign error in a coefficient).
  - CAS-engine pitfalls promoted to the pitfall library (currently
    8 documented, 1 candidate from batch III.3).
  - Orchestrator hot-fixes for canonical idioms Codex didn't apply
    preemptively.

- **Stop-cold (so far observed)**:
  - Zero `UNFIXABLE` verdicts in 96 stages.
  - Zero `CRITICAL_DOWNSTREAM` verdicts in 96 stages.

## 14. What the pipeline doesn't catch

Calibration matters: the system is thorough at the per-stage scope
and at the cross-engine scope, but there are failure modes it does
not address by construction:

- **Cross-stage premise drift.** If stage N's premise was true in
  isolation but stage N+5's combination of premises silently relies
  on a different reading of stage N's premise, the per-stage audit
  cannot see it. Defense: the upstream-stale cascade plus the
  human-gated batch transitions.
- **Paper-side prose-only errors.** If the paper card's prose
  description is wrong but the boxed equations and stated constants
  are right, the auditor reads the prose and may flag a non-
  existent issue. Defense: the auditor is required to quote both
  sides of any `paper_misalignment` finding; the user can dismiss a
  false-positive.
- **Externally-cited constants.** If a stage cites an external
  numerical constant (e.g., a Riemann zeta value, a fine-structure
  constant), the audit can confirm the script uses the same digits
  the paper states but cannot confirm the digits themselves. Defense:
  source-anchor comments are required on every literal constant.
- **Physics-validity errors at the modeling layer.** The system
  verifies that the scripts prove what the paper says they prove,
  not that what the paper says is itself a correct physics claim.
  That's the layer the user reviews directly.

## 15. Coverage and reliability

As of 2026-05-27:

- **Stages verified at v2 depth:** 96 / 253 (~38%).
- **Batches closed at v2:** 7 (I.1, I.2, II.1, III.1, III.2, III.3;
  III.4 in progress).
- **Cumulative findings closed:** ~299 across 96 stages (avg ~3.1
  findings/stage).
- **Material-change events:** 2 across all batches (stages 060 and
  068 at v1). Both confirmed clean at v2 — the v1 derivation-route
  rewrite did not change any printed symbolic content or numeric
  value.
- **Stop-cold events:** 0.
- **User-redirection rate:** 4 consecutive batches with zero
  user-redirection (every Codex math recommendation accepted by the
  user without modification).
- **Orchestrator hot-fixes:** 2 across all batches (stage 081
  `ConditionalExpression` strip retrofit; stage 064 `Integrate`
  factoring workaround, the pitfall #9 candidate).

The user-redirection rate is the system's most useful health metric.
A high rate would suggest Codex is hallucinating math and the user is
catching it; a zero rate sustained across batches suggests the math-
authority delegation is correctly placed.

## 16. Files and code

- Skill definition: `.claude/skills/redteam-audit/SKILL.md`
- Auditor prompt: `.claude/skills/redteam-audit/prompts/auditor.md`
- Verifier prompt: `.claude/skills/redteam-audit/prompts/verifier.md`
- Codex pitfall library: `.claude/skills/redteam-audit/prompts/codex.md`
- CLI entry point: `.claude/skills/redteam-audit/lib/redteam.sh`
- Manifest cache: `research/pde_ledger/redteam/MANIFEST.yaml`
- Per-stage reports: `research/pde_ledger/redteam/reports/stage_<NNN>.md`
- Per-stage directives: `research/pde_ledger/redteam/directives/stage_<NNN>.md`
- Per-batch summaries: `research/pde_ledger/redteam/batches/batch_<ID>.md`
- Per-stage codex logs: `research/pde_ledger/redteam/codex_logs/<NNN>_iter<N>.txt`
- Six prose trackers: `research/pde_ledger/notes/*.md`
