---
name: agent-roles
description: Keep each engine in the multi-agent pipeline strictly to its assigned task. Defines the per-role BEHAVIORAL contract (a leg's no-spawn / no-tree-write is held by the packet Bounds line + provenance + INCONCLUSIVE, ⛔ not a per-call sandbox or tool-strip the launch does not set), packet-absence (no cross-role instructions in a packet), the one-shot review-leg calling convention, PROGRESS-AWARE monitoring keyed on DELIVERABLE growth (⛔ never a top-level wall-clock kill — legitimate runs take hours), and the rule that the orchestrator never authors or self-clears an instrument. ⛔ This skill is for the ORCHESTRATOR only — never hand it (or the build / review-legs skills) to a builder or a review leg.
user_invocable: true
---

# Agent roles — keep each engine in its lane

⚠⚠ **Measured over-reach (2026-09-07), which this exists to prevent:**
1. The **builder** (astra), handed a directive that mentioned *"reviewed by fresh-Claude + Grok,"* spawned its OWN
   grok + claude review subprocesses and ran review-until-clear on its own output — invalid self-review, ~1h45m
   wasted after the real deliverables were done.
2. A **review leg** (grok) built its own `run_all` / `watch` ablation-orchestration scripts and **deadlocked twice**,
   instead of running a bounded review and returning a verdict.
3. The **orchestrator** (Claude) reacted by AUTHORING load-bearing tooling — a kill-watchdog, monitor loops — and
   **self-cleared** them: the single-engine bias the pipeline exists to remove.

⭐ **Enforce roles by CAPABILITY and PACKET CONTENTS first, behavioral prose second.** A prohibition alone is not a
control (`CLAUDE.md` S1). Below, "capability" = something actually bounded at launch (sandbox, `-p`, which tools are
handed); "behavioral" = what the packet asks. ⛔ Do not call a behavioral ask an enforced bound — where a bound is
only behavioral, say so, and rely on the **provenance gate** (below) + **progress-aware monitoring** to catch a
violation.

## The three roles

| Role | Behavioral contract | Capability bound that DOES exist |
|---|---|---|
| **Builder** (codex/astra) | Write its named deliverable(s); run each once; print/report. ⛔ Never review/spawn/iterate/commit. | Launched non-interactively; but under `--sandbox danger-full-access` (needed for Mathematica) it CAN still exec `grok`/`claude` and read `.claude/skills/` — so its "no self-review" is **behavioral**, backstopped below. |
| **Review leg** (codex / grok / fresh-Claude) | Read the ONE target artifact **+ its enumerated input manifest**; run the NAMED checks (writing bounded /tmp check-scripts is its *method*); write ONE report **+ named evidence files under a write fence**. ⛔ Never build `run_all`/`watch`/supervisor orchestration, spawn agents, modify the tree, repair, or iterate. | **Behavioral for every leg** — held by the packet's Bounds line + the provenance gate + INCONCLUSIVE-on-violation. A leg **must** write /tmp check-scripts + stdout (its method), so ⛔ `-s read-only` is wrong (it blocks those writes); the Agent tool has no allowlist arg and grok runs `bypassPermissions`. ⚠ A *real* fs fence would be `--sandbox workspace-write` with cwd in a /tmp scratch (repo outside the writable root, /tmp writable) — a hardening not yet wired into the launch; until it is, ⛔ do not claim a hard bound. |
| **Orchestrator** (Claude) | Write packets/directives/specs/**step records**; launch; **mechanical fact-lookup**; adjudicate; advance state. ⛔ Never author a CAS instrument, a build-deliverable, or an ops **decision** script; never review or clear its own output. | Self-imposed; the discipline in this skill is the enforcement. |

**Provenance gate (the real backstop for a behavioral bound):** every artifact records its author; **its author's
own "review" of it is disregarded entirely** — it is not evidence and cannot clear anything; only an *independent*
leg's verdict advances state. So even if a builder spawns its own review, that output is discarded on sight, and the
wasted work is caught as **no deliverable progress** (monitoring, below).

## Packet absence — necessary, not sufficient (S1)
- **Builder packet** = the artifact spec + supplied equations + "run once, report." ⛔ **No** "reviewed by…", other
  engine, "until clear," or authorship-pairing sentence — that one line caused failure 1.
- **Gate order:** first `rg -ni 'review|grok|claude|\bleg\b|fresh (claude|agent)|review-until-clear|subagent|/build|/review-legs'`
  over the *authored body* and strip every hit; **then** append the fixed end-of-task epilogue (below) — the epilogue
  is a constant, not authored content, so it is exempt from the strip. (⛔ Do NOT put `comparator` in the strip regex —
  a T7-comparator build legitimately names it.)
- **End-of-task epilogue (verbatim):** *"YOUR TASK ENDS AT build→run-once→report. Do not review, launch or spawn any
  agent or second engine, iterate, read the /build or /review-legs skills and act on them, or commit."*
- ⚠ **Honest limit:** the gate removes the *trigger*; under `danger-full-access` it does **not** remove the builder's
  ability to spawn or to read the skills off disk. That residual is closed **not** by an auto-kill watchdog (failure
  3) but by the **provenance gate** (its self-review is discarded) **+ progress-aware monitoring** (out-of-role
  children show up as no-deliverable-progress and are stopped as a judgment).
- ⛔ **Never hand the `build`, `review-legs`, or `agent-roles` skills to a builder or a leg** — for a full-access
  builder the on-disk tree is effectively in-packet, so keep these read-relevant only to the orchestrator's own runs.

## Calling a review leg so it stays a review (one-shot)
"One-shot" ≠ the CLI name — `codex exec`, `grok`, and a fresh `Agent` are all agentic. It means:
1. **Be honest about what the launch bounds:** the "no spawn / no tree-write" is currently **behavioral for every
   leg** — held by the packet's Bounds line + the provenance gate + INCONCLUSIVE-on-violation. ⛔ `-s read-only`
   is NOT the bound (the leg must write its /tmp check-scripts + stdout); the Agent tool has no allowlist arg; grok
   runs `bypassPermissions`. ⛔ Do not claim a hard bound the launch does not set.
2. **Closed packet:** exact artifact path, the input manifest, the NAMED checks (or "run these committed scripts"),
   a single `REPORT_PATH` (+ a `/tmp` evidence fence), and the verdict schema (`VERDICT: SOUND|NOT-SOUND|INCONCLUSIVE`,
   `FINDINGS[]`, `EVIDENCE[]`, does-it-change-what-is-computed-or-claimed).
3. **Write fence** (behavioral — see the launch note above): only `REPORT_PATH` + a `/tmp` scratch should be
   written; the working tree stays untouched (ablate copies).
4. **Bounds block — the last thing in every review packet** (the canonical wording is the `## Bounds` block in the
   `review-legs` prompt template): it must require, at minimum, *"write your report and exit; do not spawn agents,
   build orchestration/watch/supervisor scripts, modify the working tree, repair, or iterate."*
5. A leg that violates scope or makes no progress is **INCONCLUSIVE** — the orchestrator relaunches a fresh one-shot;
   the leg never manages its own retries. Review-until-clear is the **orchestrator's** loop, never the leg's.

## Progress-aware monitoring — ⛔⛔ never a TOP-LEVEL wall-clock kill
⚠ Legitimate agent/build runs are long — CAS steps take 20+ min; a codex build ran ~3 hours while it iterated. ⛔ Do
**not** wall-clock-kill the top-level agent/build. Watch it **progress-aware** instead:
- **Progress = growth of the NAMED DELIVERABLE** (the target file/`.out`/report bytes) — ⛔ **not** CPU and ⛔ **not**
  raw log/heartbeat growth (a runaway self-review or a livelocked `FullSimplify` burns CPU and grows a log while
  producing nothing — treating that as progress reopens failure 1).
- **Check cadence:** a soft checkpoint every ~10 min. **No new deliverable bytes for ~3 consecutive checks** (≈30
  min) ⇒ **investigate**, do not kill.
- **Investigate = look, then judge:** read the partial output + the process tree yourself, including a dumb
  `pgrep`/`/proc/<pid>/cmdline` scan for **out-of-role children** (e.g. `grok`/`claude` under a builder's session =
  role violation). Then *you* decide — keep waiting, accept the partial, relaunch, or **stop specific out-of-role
  children as a deliberate judgment**. ⛔ Never an automatic timeout-kill; the monitor's only job is to **summon you
  to look** (the `Monitor` tool emitting one event on no-deliverable-progress is the mechanism).
- ⭐ **Per-kernel budgets are SEPARATE and stay:** a child CAS kernel wrapped in `timeout 600` (a *predeclared
  bounded sub-computation*, per `review-legs`) is NOT a top-level wall-clock kill — keep it. The prohibition is
  against killing the **agent/build** on elapsed time, not against bounding one child kernel.

## The orchestrator never authors OR self-clears an instrument
- ✅ **No review needed — an individual OBSERVATION command:** a single `test -s` / file-mtime read / `pgrep`+`%cpu` /
  `ls` / a lone `timeout` on one child kernel. Its output is not validation evidence and it decides nothing.
- ⛔ **Instrument — Codex-written + G1-reviewed BEFORE use:** ANY script that *composes* observations into a
  **decision** — a loop/branch that concludes done / hung / failed / passed, gates a launch, collects evidence, or
  adjudicates a physics/process claim. A watchdog or a stall/health monitor **is** an instrument (E1: the author
  fixes the framing then trusts the output). The dumb-primitive exemption covers the *primitives*, ⛔ never their
  composition into a status verdict.
- ⛔⛔ **Forbidden:** the orchestrator authors a decision-script **and** uses it as a gate in the same session.
- Default when something hangs: **bound + investigate the callee**, ⛔ not grow the caller with a supervisor.

## Launching a long detached job (dumb primitives, ⛔ no kill loop)
`run_in_background` is reaped on this box, so a long build/leg is launched `setsid`-detached with a `DONE`-marker —
a thin one-liner, **not** a supervisor:
```bash
setsid bash -c 'cd /var/projects/toy_physics; <cmd> > "$LOG" 2>&1 < /dev/null; echo "EXIT=$?" >> "$LOG"; touch "$DONE"' >/dev/null 2>&1 < /dev/null &
```
Then watch `$DONE` + deliverable growth with the progress-aware `Monitor` above. ⛔ No embedded watchdog, no
grok/claude kill loop, no status-deciding composition in the launcher.

## Minimum standing rule (above everything)
**The process that WRITES an artifact is never in the process tree that REVIEWS it; a review process spawns nothing
and writes nothing but its report + named evidence; and the orchestrator never both authors and trusts a
decision-script.** Enforce with capability bounds where they exist, packet-absence, the provenance gate, and
progress-aware investigation — ⛔ never with an in-session tool the orchestrator both writes and clears.

Related: `[[feedback_builder_directive_no_orchestrator_process]]`, `[[feedback_thin_orchestrator_definition]]`,
`[[feedback_orchestrator_never_authors_cas_script]]`, `[[feedback_review_agents]]`, `.claude/skills/build/SKILL.md`,
`.claude/skills/review-legs/SKILL.md`.
