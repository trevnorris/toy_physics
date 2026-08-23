# Resumability check for a /compact hand-off (S11b closed → S11c is the next front)

The unified step **S11b** ("the linear brane–bulk interface coupling law") was just CLOSED, and the deferred
non-uniform step C has been **renamed `S11c`**. Verify that a FRESH context (after /compact) can correctly
resume S11c from the committed docs, with no contradiction or stale pointer that would mislead it. This is a
consistency/resumability audit — ⛔ not a physics review.

## Read (form your own view of "where are we and what's next")
- `STATUS.md` — the top "CURRENT FRONT — S11b" block and its NEXT.
- `research/pde_ledger_v3/steps/S11c_SCOPE.md` — the new S11c scope doc (the "next plan").
- `research/pde_ledger_v3/V3_STEP_PLAN.md` — the `#s11b-split` section (search "S11b IS SPLIT").
- `research/pde_ledger_v3/steps/S11b_RUN_CHECKLIST.md` — the closed run checklist.
- `research/pde_ledger_v3/steps/S11b_HANDOFF.md` — the historical C-requirements source (the scope doc
  consolidates it).
- Memory (out-of-repo; read with `--sandbox danger-full-access`):
  `/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11b_interface_law_result.md`
  and `.../memory/MEMORY.md` (the S11b hook line).

## Check
1. **Consistency.** Do STATUS, V3_STEP_PLAN, the scope doc, the checklist, and the memory all AGREE that
   S11b is CLOSED (`565b3fe8`, all 13 steps) and that **S11c** (not "S11b-C", not "part of S11b") is the
   next step? Flag any live contradiction — e.g. a doc still asserting "S11b is not closed until C", "C is
   not a separate step / no renumbering", or "the immediate NEXT is the blind Wolfram engine". (⚠ A clearly
   marked *superseded/history* note is fine; a live assertion is not.)
2. **Resumability.** Is the resume path unambiguous: start at `S11c_SCOPE.md` → author the **S11c decision
   list** (orchestrator-written, 2 legs, rule 7 TRIGGER) → spec → engines (SymPy + blind WL) → comparator →
   step record → card? Does the scope doc resolve every doc it points at?
3. **Completeness.** Does `S11c_SCOPE.md` capture the C requirements faithfully against `S11b_HANDOFF.md` —
   the five requirements (tilted faces; Eulerian-vs-material density; plane-waves-not-eigenmodes / ⛔ no
   global dispersion relation; uniform-limit control known-vacuous; falsification-first bench-top-optics
   bound), the carry-ins (`O(v₀|q_n|/ω)`, the frozen-wall reconciliation), the G14 scope boundary (nonlinear
   radiation audit is NOT S11c), and the "split finer than S11b's 11 revisions" lesson? Anything dropped or
   misstated?
4. **Stale pointers.** Any forward-looking/live doc still pointing at S11b work as NEXT, or naming the next
   step "S11b-C" where it should be `S11c`? (⛔ Historical/archival mentions of "C" or "S11b-C" in closed
   records are fine — flag only LIVE forward pointers.)

## Output
Report ONLY genuine resumability/consistency problems (a fresh context would be misled or blocked), ⛔ not
style. For each: the file + line, the problem, and the fix. If clean, say so explicitly per check. End with a
one-line bottom line: can a fresh context correctly resume S11c as written, or the specific fixes needed.
