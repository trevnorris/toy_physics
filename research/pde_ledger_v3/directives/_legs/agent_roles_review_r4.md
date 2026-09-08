# Review ROUND 4 (final confirm) — role-discipline docs

Rounds 2–3 confirmed the design SOUND; the last open item (a "stripped tools" capability overclaim that survived in
the calling-convention section, plus the review-legs codex command missing `-s read-only`) is now folded. This is a
final check. Repo read access; concise; cite lines.

## Files
`.claude/skills/agent-roles/SKILL.md`, `.claude/skills/review-legs/SKILL.md` (and `.claude/skills/build/SKILL.md`, `CLAUDE.md` unchanged since round 3)

## Confirm
1. **No "stripped tools" overclaim remains** anywhere — the calling-convention section (agent-roles ~:52-59) now
   matches the honest table row (:29): codex leg = real `-s read-only` fs bound; a Claude `-p` launch = real
   non-interactive bound; Agent has no per-call tool-strip and grok is `bypassPermissions`, so those are
   behavioral (Bounds + provenance + INCONCLUSIVE). The write-fence is labeled real-for-codex / behavioral-otherwise.
2. **review-legs codex command** now includes `-s read-only` (the one real fs bound).
3. **Any NEW contradiction or overclaim** introduced by these edits? Any place an agent could still be told it has a
   bound the launch does not set, or be told to leave its lane?

(Already agreed non-blocking, do not re-raise as blocking: review-legs' in-session CLI parallel-launch still uses
`run_in_background` — that is the correct scope for a short in-session leg, distinct from the hours-scale detached
builder.)

## Output
Per-item RESOLVED/open + any new defect, then **SOUND** (ready to commit) or **NOT-SOUND** (file+line+fix).
