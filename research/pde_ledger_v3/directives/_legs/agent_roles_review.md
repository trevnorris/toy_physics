# Review — new role-discipline process docs (orchestrator-written)

Review these process/discipline docs for soundness before they are committed. You have repo read access. Answer
concisely; ground each point in a cited line. This is NOT physics — it is pipeline role discipline.

## Artifacts (read all)
- NEW skill: `/var/projects/toy_physics/.claude/skills/agent-roles/SKILL.md`
- Revised: `/var/projects/toy_physics/.claude/skills/build/SKILL.md` (the §isolation section + the step-2 gate + the
  step-3 launch note were changed; a `astra_launch.sh` kill-watchdog was REMOVED)
- CLAUDE.md: the new "Role discipline (agent-roles skill)" paragraph near the operational-runbooks note (end of file)

## Background (the failures these must prevent)
1. A **builder** (codex) handed a directive mentioning "reviewed by fresh-Claude+Grok" spawned its own grok+claude
   review and ran review-until-clear on its own output (invalid self-review, ~1h45m wasted).
2. A **review leg** (grok) built its own run_all/watch orchestration and deadlocked twice instead of a bounded
   review + verdict.
3. The **orchestrator** (Claude) reacted by authoring a kill-watchdog + monitors and self-cleared them.
Also: Claude Code one-shot (`-p`) is NOT per-token, so this is role discipline + reliability, not cost.

## Questions
1. **Do these rules actually stop all three failures?** Point to the specific rule that blocks each. Any failure
   still reachable?
2. **Internally consistent + consistent with existing CLAUDE.md controls** (S1 blindness-by-absence; G1 whatever
   writes does not review; E1 orchestrator never authors the instrument; G3 launch-on-sight; G4 orchestrator owns
   the review loop)? Flag any contradiction or redundant/weaker restatement.
3. **The progress-aware monitoring rule** (⛔ never wall-clock-kill; check "was additional work done?" → progress =
   leave, no-progress = investigate not kill): is it sound and implementable purely with dumb primitives
   (`test -s`, file mtime, `pgrep`+%cpu, `ls`)? Does it leave any hang un-handled?
4. **Was dropping the watchdog correct** — is packet-absence (the step-2 gate + end-of-task clause) genuinely
   sufficient to stop a builder's self-review, or does removing the watchdog reopen failure 1?
5. **Any remaining over-reach, gap, ambiguity, or place a builder/leg could still be told to do the wrong thing?**
   Is the `agent-roles` tool-allowlist table complete and correct?

## Output
Per-question findings with cited lines, then a final line: **SOUND** (ready to commit) or **NOT-SOUND** (exact
file + line + fix for each blocking issue).
