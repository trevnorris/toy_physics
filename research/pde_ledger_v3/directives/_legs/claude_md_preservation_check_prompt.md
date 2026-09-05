# Preservation check — CLAUDE.md restructure (did the rewrite drop or weaken any control?)

`CLAUDE.md` is the governing rule list for this repository. It was just **restructured** (17 rules regrouped
into M/E/G/S sections + a quick-reference + an artifact→review-discipline table + an evidence ledger). The
restructure must change **presentation only, never authority**: ⛔ no control dropped, weakened, or narrowed;
every measured-evidence war-story preserved; every one of the 17 original rules mapped. Your job is to verify
that — or find exactly what was lost or altered. ⛔ Do not rubber-stamp.

## Files
- **OLD (source of truth for content)**: `/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/CLAUDE_old.md` (101 lines, the pre-restructure rules).
- **NEW (the artifact under review)**: `/var/projects/toy_physics/CLAUDE.md` (the restructured version).
- **DESIGN + preservation map** (the intended mapping, for reference — but verify against OLD, not just this):
  `/var/projects/toy_physics/CLAUDE_streamline_proposal_2026-09-05.md`.

## What to verify (be exhaustive and adversarial)
1. **Every normative clause in OLD survives in NEW** — possibly relocated/reworded, but with the same
   obligation and the same strength. Go rule by rule (OLD rules 1–17 + preamble + the "Repository
   infrastructure" section). For each, point to where its content lives in NEW, or flag it as
   DROPPED/WEAKENED/NARROWED with the exact OLD line and what's missing.
2. **Every measured war-story / dated failure in OLD is preserved verbatim** (the 2026-08-12 four-designs,
   the 2026-08-09 six-defects, the 26=26 coincidence, etc.) — check the NEW "Evidence ledger". Flag any
   lost, truncated, or altered evidence.
3. **The "Repository infrastructure" section is verbatim** (the annex/GIN policy, both `out/*.out` globs, the
   `datalad`/`git push` commands, "never `git add -f`", "never annex an `*_exports.py`", the `gin` credential,
   the `project-datalad-gin-out-storage` link). Flag any change.
4. **No NEW error introduced by the rewrite** — a mis-stated rule, a wrong disambiguation, or an accidental
   weakening in the artifact→review-discipline table or the quick reference. In particular check that:
   - the **A1 disambiguation** is correct: rule-7's "one pass, fold and go" applies to the **decision list
     only**; a physics **spec / script / record** is **review-until-clear**. Is that stated correctly and
     consistently across the quick-ref, the table, G2, and G4?
   - the artifact table's reviewer column matches OLD rule 7's authorship pairings (Orchestrator → Codex +
     Grok; Codex → fresh Claude + Grok) and does not invent a third pairing or reduce the leg count.
   - the "commit the reviewed baseline before repair" and "no commit before both legs" clauses coexist
     without contradiction (OLD rule 9 + the build-skill baseline rule).
5. **Nothing was silently added that changes authority** — a new obligation not traceable to OLD or to the
   skills is itself a finding (the restructure is presentation-only).

## Method
This is a document fidelity check: OLD is the content source of truth, NEW is the artifact. Quote both sides
for every finding (OLD line → NEW location, or "absent"). Rank findings by severity: DROPPED-CONTROL (a
normative obligation lost) > WEAKENED (obligation present but softer/narrower) > MISSTATED (a new error) >
NIT. No SymPy needed.

## Output
A per-rule preservation table (OLD rule → NEW location, verdict). Then the ranked findings. End with an
explicit verdict: **PRESERVATION-COMPLETE (nothing dropped, weakened, or misstated — safe to commit)** or
**BLOCK** with the exact list of losses/alterations to fix. If it is genuinely faithful, say so plainly.
