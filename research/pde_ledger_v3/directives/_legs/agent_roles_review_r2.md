# Review ROUND 2 — role-discipline docs (confirm the round-1 fixes)

Round 1 came back NOT-SOUND. Every finding was folded. Confirm each fix is correct and no new defect was introduced.
Repo read access; answer concisely; cite lines. Not physics — pipeline role discipline.

## Files
- `.claude/skills/agent-roles/SKILL.md` (rewritten)
- `.claude/skills/build/SKILL.md` (§isolation, step-2 gate, step-3 launcher, step-4)
- `.claude/skills/review-legs/SKILL.md` (new "## Bounds" terminal line in the prompt template, ~line 71)
- `CLAUDE.md` (the "Role discipline (agent-roles skill)" paragraph near the end)

## Confirm each round-1 fix
1. **Self-referential gate/clause:** is it now resolved by gating the *authored body* first, then appending the
   fixed end-of-task clause as an exempt constant (agent-roles "Gate order"; build step-2)? Can a conforming packet
   now pass its own gate?
2. **Overclaim:** do the docs now state HONESTLY that packet-absence removes only the *trigger* (not spawn capability
   under danger-full-access), with the backstop = provenance-gate (self-review discarded) + progress-aware monitoring
   (agent-roles "Honest limit"; build §isolation)? No remaining "packet-absence is sufficient" claim?
3. **Progress = deliverable growth:** does the monitoring rule now key on NAMED-DELIVERABLE growth (⛔ not CPU/log),
   define a cadence (~10min ×3), name the summon mechanism (Monitor), and add the out-of-role-children check
   (agent-roles "Progress-aware monitoring")? Is a runaway self-review now caught (no deliverable growth)?
4. **Wall-clock-kill scoping:** is "never wall-clock-kill" now scoped to the TOP-LEVEL agent/build, with per-kernel
   `timeout 600` budgets explicitly preserved (agent-roles; CLAUDE.md)? No conflict with review-legs' kernel budget?
5. **setsid+DONE launcher restored:** does build step-3 now give a thin `setsid`+DONE one-liner (dumb primitives, no
   kill loop), and does step-4 no longer claim "the harness re-invokes you"?
6. **Dumb-primitive exemption narrowed:** is it now clear that a single OBSERVATION is exempt but any COMPOSITION that
   DECIDES status is an instrument (agent-roles; CLAUDE.md) — so it cannot authorize a watchdog?
7. **review-legs terminal line:** does the prompt template now end with the mandatory "write report and exit; no
   spawn/orchestration/iterate" bounds line?
8. **Honesty/scope:** is the role table now honestly labeled behavioral-contract + capability-bound-that-exists (not
   a false "enforced allowlist")? Is agent-roles marked orchestrator-only (never handed to a builder/leg)? Is
   `comparator` removed from the strip regexes? Is "orchestrator writes step records + mechanical fact-lookup"
   consistent with E1 (not an over-broad "never author build-code")?

## New defects?
Any new contradiction, gap, or place an agent could still be told to leave its lane?

## Output
Per-item RESOLVED/still-open with cited lines, then a final line: **SOUND** (ready to commit) or **NOT-SOUND** (exact
file+line+fix).
