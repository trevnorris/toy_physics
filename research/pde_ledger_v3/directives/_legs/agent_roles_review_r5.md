# Review ROUND 5 (final confirm) — role-discipline docs

Round 4 found that `-s read-only` on a codex leg blocks its required /tmp evidence writes. Resolution: the
sandbox-bound claim was REMOVED — a review leg's no-spawn/no-tree-write is now uniformly **behavioral** (packet
Bounds line + provenance gate + INCONCLUSIVE), with a note that a real fence would be `workspace-write` scoped to a
/tmp cwd (not yet wired). Confirm this subtractive fix is clean. Repo read access; concise; cite lines.

## Confirm
1. No `-s read-only` remains in any launch command, and no prose still calls a leg "read-only" as if enforced
   (agent-roles table :29 + calling-convention :55 + write-fence :62 + frontmatter; review-legs codex command ~:104
   and the prose ~:101; build; CLAUDE.md).
2. The write-fence / no-spawn is stated **behaviorally and consistently** for every leg (Bounds + provenance +
   INCONCLUSIVE), with the honest "a real fence would be workspace-write scoped to /tmp, not yet wired" caveat — no
   residual claim of a hard bound the launch does not set.
3. Any NEW contradiction introduced by this removal?

## Output
Per-item RESOLVED/open, then **SOUND** (ready to commit) or **NOT-SOUND** (file+line+fix).
