# Compact-prep doc audit — S11c-a close

You are auditing whether this session's DOC + MEMORY updates accurately reflect the real repository state,
before a context compaction. Report ONLY discrepancies: a doc claim that contradicts git / the final run /
the engines, a stale or misleading statement, an inaccurate sha/path/number, an overclaim, or a
citation to a file that does not say what the doc attributes to it. If a doc is accurate, say so briefly.
Do NOT fix anything — just report. Run from repo root `/var/projects/toy_physics`.

## Background
S11c-a (the first S11c sub-step, background & geometry) was just closed. The two-engine method surfaced
FIVE single-engine defects, all fixed and committed. An earlier "the two engines AGREE; every residual
representational" conclusion was retracted (it rested on a blanket-collapse classifier that hid two
findings). The step record was rewritten; its two review legs (Codex + Grok) caught SIX overclaims in the
rewrite, which were then folded. Your job is to confirm the folded record + the STATUS/memory updates are
honest and correct.

## Ground truth to check
- `git log --oneline -8` — HEAD should be `3b552426` (step-record doc-prep), atop `cccb4f9e` (WL
  density-grounding, 5th fix). The five fix shas: c36beac4, 6fae82b8, 49b5c525, 8c1a5ed1, cccb4f9e —
  confirm each exists and does what the record says (`git show --stat`).
- The final comparator run `~/.s11_build/comparator_final_cccb4f9e.out`: `RUN_ACCOUNTING` should show
  `families=39 families_with_unpaired=11 parse_failed=0`. Grep the 11 `^ACCOUNTING` lines with
  `py_only`/`wl_only` > 0 — confirm they are ALL control/bookkeeping families (ADMISSIBILITY_PREMISE,
  REP_INVARIANCE_*, CONTROL_INDEPENDENCE_*, CONTROL_FORM_*, DIMENSIONS), NONE a physics family. Confirm
  `FACE_SHIFT {join=160, py_only=0, wl_only=0, axis_set_mismatch=0}` and
  `UNIFORM_LIMIT_RESIDUAL {join=332, …=0}`.

## Docs/memories to audit (absolute paths)
1. `research/pde_ledger_v3/steps/S11c_a_interface_shape_derivatives.md` — the rewritten step record. Verify
   the SIX previously-caught defects are actually fixed and none re-introduced:
   (a) PROJECTION_* and FACE_SHIFT are NOT in a "clean zero (both engines identical)" list — they are
       under representational-identity classes (FACE_SHIFT density verified 0; current = applied-at-face ↔
       bare-symbol; PROJECTION under the δρ identity). Confirm the "Clean zero" list holds only literally-
       zero families.
   (b) The supplied/bookkeeping families (BACKGROUND_STATE, FACE_MAP_*, BACKGROUND_DENSITY_MAP,
       ADMISSIBILITY) are classified as structural/naming, not derived physics.
   (c) UNIFORM_LIMIT is described as fully joined (332), NOT unpaired; only REP_INVARIANCE/
       CONTROL_INDEPENDENCE carry the PY-missing-DENSITY unpaired.
   (d) The central "no physics disagreement survives" claim is attributed to the comparator raw + per-family
       reconciliations, NOT solely to `S11c_a_handcoded_comparison.py`; that script is cited at `bb2a050a`
       (not cccb4f9e) and described as FLAGging the bridge-reducible families by design.
   (e) CONORMAL verdict-A is cited to `directives/_measurements/S11c_a_comparator_reemit_plan.md:82`.
   (f) The control battery says CONTROL_FORM bites 16/18 (DEAD only for FACE_SHIFT + PROJECTION_STATIC),
       CONTROL_INDEPENDENCE bites all 6, HOLD holds — matching
       `directives/_measurements/S11c_a_bite_liveness.py` output and the CORRECTION section appended to
       `directives/_measurements/S11c_a_control_battery_result.md`.
   Flag any surviving overclaim, any claim not backed by a cited measurement, or any physics misstatement.
   In particular: does the record anywhere still assert a clean unqualified "engines agree / 39 clean"?
2. `STATUS.md` top (`## CURRENT FRONT`) — verify it now says S11c-a CLOSED (`3b552426`), lists the five
   fixes, states the honest cross-comparison (physics families clean; 11 control/bookkeeping keying owed),
   NEXT = S11c-b, and does not still say "pending this session" or repeat the retracted overclaim.
3. `/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_a_state.md` and the
   MEMORY.md index line for it — do they match the corrected/closed state (5 fixes, honest, owed items,
   NEXT=S11c-b)? Any stale "engines agree / 3 fixes / step record pending" wording?

## Specifically flag
- Any surviving unqualified "engines agree / every residual representational / 39 clean".
- Any claim that a control/bookkeeping keying asymmetry is a physics disagreement, or vice-versa.
- Any sha/path/number git or the run does not support.
- Any measurement citation whose file does not contain what the doc says.

Return a short list of discrepancies (file + quoted text + the ground truth it contradicts), or "in the
clear" if none.
