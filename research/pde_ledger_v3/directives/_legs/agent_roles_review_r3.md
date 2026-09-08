# Review ROUND 3 (confirm) — role-discipline docs

Round 2 confirmed the DESIGN sound; the four remaining items were wording/consistency and are now folded. Confirm
each, and flag any NEW defect. Repo read access; concise; cite lines.

## Files
`.claude/skills/agent-roles/SKILL.md`, `.claude/skills/build/SKILL.md`, `.claude/skills/review-legs/SKILL.md`, `CLAUDE.md`

## Confirm the four round-2 fixes
1. **build step-3 self-contradiction** — is it now ONE coherent launch method: a thin `setsid`+`DONE` one-liner with
   a generic `<CMD>` (no hardcoded model/effort/sandbox; `xhigh` default; `danger-full-access` only for Mathematica),
   and NO leftover "`run_in_background: true` / no `&` / harness notifies on exit" contradiction? (build step-3.)
2. **review-leg capability overclaim** — does the capability column now say honestly that only the **codex** leg has a
   real `-s read-only` fs bound, while **grok / fresh-Claude-Agent** legs are **behavioral** (no per-call tool-strip;
   grok uses bypassPermissions), held by the Bounds line + provenance + INCONCLUSIVE — with no false "stripped tools"
   claim (check the frontmatter too)? (agent-roles table row + frontmatter.)
3. **CLAUDE.md over-statement** — is "Enforce by CAPABILITY and PACKET ABSENCE" replaced by "capability bounds where
   the launch actually sets them, packet absence, the provenance gate, and progress-aware investigation"?
4. **terminal-line verbatim mismatch** — does agent-roles now require a Bounds *block* "at minimum …" pointing to the
   `review-legs` `## Bounds` template as canonical, rather than demanding a byte-identical line?

## New defects
Any new contradiction, hardcode, honesty gap, or place an agent could still be told to leave its lane? (Non-blocking
already-noted: review-legs' own parallel-launch section still uses run_in_background for CLI legs — say if that now
rises to blocking given step-3's setsid, or stays non-blocking.)

## Output
Per-item RESOLVED/open with cited lines, then **SOUND** (ready to commit) or **NOT-SOUND** (file+line+fix).
