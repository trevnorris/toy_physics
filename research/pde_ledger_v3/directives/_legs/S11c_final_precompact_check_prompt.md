# FINAL pre-/compact verification — S11c prep + the staleness fixes

We are about to `/compact`. This is the last check before losing context. Confirm a FRESH post-/compact
context can open **S11c** (the non-uniform transverse coupling; S11b is CLOSED) with NO contradiction, stale
live pointer, or gap that would block or mislead it. ⛔ This is a consistency/resumability audit, ⛔ not a
physics review. ⚠ A prior audit flagged six stale live surfaces; they were then fixed — your job is to
confirm the fixes landed, introduced nothing new, and that the whole hand-off is now clean.

## What changed (review these diffs first)
```
git -C /var/projects/toy_physics show fd706dae   # prep S11c: new steps/S11c_SCOPE.md + rename forward front
git -C /var/projects/toy_physics show 38f90b1a   # clear the six stale live surfaces the prior check flagged
```
(Also relevant, already committed: 565b3fe8 closed S11b; the card/step-record trail before it.)

## Then read the CURRENT state of the live hand-off surfaces
- `STATUS.md` — the top "CURRENT FRONT — S11c" block AND the deeper "[HISTORICAL SNAPSHOT 2026-07-31 …]"
  section (confirm the latter is now clearly marked historical, not live).
- `research/pde_ledger_v3/steps/S11c_SCOPE.md` — the S11c plan (source of record for the next front).
- `research/pde_ledger_v3/steps/S11b_RUN_CHECKLIST.md`, `research/pde_ledger_v3/steps/S11b_HANDOFF.md`
  (its top banner), `research/pde_ledger_v3/LAUNCH_PROMPT.md` (its top banner + the old "NEXT: S11b" line),
  `research/pde_ledger_v3/V3_STEP_PLAN.md` (the `#s11b-split` section + line ~850 MacCullagh rec).
- Memory (out-of-repo; read with `--sandbox danger-full-access`):
  `/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11b_interface_law_result.md`
  and `.../memory/MEMORY.md` (the S11b hook line).

## Confirm, explicitly (say PASS/FAIL for each)
1. **The six previously-flagged surfaces are RESOLVED** — no live text still says "S11b is being rebuilt",
   "the immediate next is the blind Wolfram engine", "NEXT: S11b", "SymPy X-1 repair (the current NEXT)", or
   names the next step "S11b-C" as a live pointer. (⚠ Superseded/historical/DEAD-marked mentions are fine.)
2. **No NEW contradiction introduced** by the fixes — do STATUS, the scope doc, the checklist, the handoff
   banner, LAUNCH_PROMPT, V3_STEP_PLAN, and the memory all now AGREE that S11b is CLOSED (`565b3fe8`) and
   S11c is the next step, starting at `steps/S11c_SCOPE.md`? Any hash/step-count/naming disagreement?
3. **Resume path is unambiguous** — a fresh context that reads STATUS's top block → `S11c_SCOPE.md` lands on
   "author the S11c decision list (2 legs, rule 7 TRIGGER), then spec → SymPy+blind-WL engines → comparator →
   step record → card". Every doc the scope references resolves.
4. **Nothing material dropped** — the S11c scope still faithfully carries the five requirements, the carry-ins
   (`O(v₀|q_n|/ω)`, frozen-wall reconciliation), the G14 boundary, and the split-finer lesson.

## Output
PASS/FAIL per item, with file+line and an exact quote for any residual problem and its fix. If fully clean,
say so explicitly. End with a one-line bottom line: are we CLEAR to /compact, or the specific fixes needed
first.
