# REVIEW ROUND 2 — confirmation pass. Read-only, one pass.

Round 1 reviewed the v3 ledger plan set at `2db79549`. Both reviewers returned blocking findings.
They have been folded. **HEAD is now `18eaf2ad`.** Your job: check whether the fixes are real, and
whether folding them introduced new problems.

Repository root `/var/projects/toy_physics`, branch `ledger-v2-rebuild`.

**Read, in order:**
1. `research/pde_ledger_v3/V3_STEP_PLAN.md`   ← the main object
2. `research/pde_ledger_v3/CHARTER.md`
3. `research/pde_ledger_v3/DEFECT_REGISTER.md`
4. `research/pde_ledger_v3/SESSION_REASONING.md`

`git diff 2db79549..HEAD -- research/pde_ledger_v3/` shows exactly what changed.

## The round-1 blockers that were folded — verify each is ACTUALLY fixed

1. **S12 used a drain law the user ruled out.** Now the dynamical `Γ_B` order-conversion law; the
   frozen-wall `ρ`-sink is marked RULED OUT. Check against `stage045_nonvariational_block_prep.md`.
2. **S16 overclaimed the worldtube reduction.** Now a *conditional response-side monopole
   approximation*, with supplied inputs listed. Check against `4d_1pn_full.tex` / `4d_2_5pn.tex`.
3. **The mass–radius claim was false.** `SESSION_REASONING` §5 and register `D1` now say no *unique,
   target-blind, empirically successful* law survives, several conditional slopes remain. Check
   against `notes/lepton_mass_notes.md` — **all 3239 lines matter here, not the first 200.**
4. **No pre-S1 registry surgery.** Added as S0.5 (`c_γ`/`λγ` leave the medium contract; recompute
   acceptance) plus a substrate-action step before S2.
5. **Cone lock had no step.** Added as S20a.
6. **S22's debt table was incomplete.** Expanded.
7. **S15 preceded the theorem licensing it.** Order corrected.
8. **S19 must audit, not inherit, the orphan `tadpole.md` note.**

## Operating constraints
- **READ ONLY.** Do not modify or commit. **One pass, no clarifying questions** (you will hang).
- You may run Python/SymPy; write scratch outside the repo.
- ⭐ Open the cited file rather than trusting the document's summary of it.

## Output

```
# Round 2 — <your name/model>

## VERDICT
<CLEAN — SAFE TO EXECUTE / STILL BLOCKING / REGRESSED> + 2-3 sentences

## PER-BLOCKER  (1-8 above)
<one line each: FIXED / PARTIAL / NOT FIXED / NEW PROBLEM, with file:line>

## NEW BLOCKING
<anything the folding introduced, or that round 1 missed. Empty is legitimate.>

## THE FIRST THREE STEPS
S0.5, S1, S2 will be executed side by side with the user next session. What concretely goes wrong
in those three? Be specific.

## MATH FLAGS
<table: claim | file:line | your result | agree/disagree>
```

## Standards
A matching number is not evidence. Dimensional agreement is not physical agreement — ask whether both
sides are indexed by the same thing. Falsification is welcome and first-class. Apparatus above physics
has killed two efforts here. Absence of a denial is not evidence.

⛔ Do not be agreeable. A finding that stops execution is the most valuable output.
