# S11 requirements pass (pass 2): directive for Codex

**Author:** Claude (orchestrator), 2026-10-05. **Builder:** Codex. `AGENTS.md` applies.

**Deliverable:** `SUBSTRATE_REQUIREMENTS.md` gains an entry for every object that a kept S11, S11b-A, S11b-B
or S11c result rests on, made by the register's own schema and method.

**Why now:**
- **Pass 2 never ran.** `SUBSTRATE_REQUIREMENTS.md:299–300` schedules it.
- **This closes S11.** It is step 1 of the light-sector close-out.
- **Review:** this directive has had its two review legs. The entries you write get their own two legs at
  the STOP.

## Read

- **Schema and method:** `SUBSTRATE_REQUIREMENTS.md`:
  - "Entry schema", including its status values;
  - the pass method under "Pass 2";
  - the "second route".

  Follow them as written. Where this directive and the register disagree, the register wins, and you report
  the disagreement.
- **Records:**
  - `steps/S11_stray_longitudinal.md`;
  - `steps/S11b_interface_coupling_law.md`, `steps/S11bA_interface_response.md` and
    `steps/S11bB_interface_assembly.md`, all three in full;
  - `steps/S11c_PARTIAL_CLOSEOUT.md` and the S11c a, b, c1, c2 and d records;
  - for rule 6 only: `steps/S9_light_requires_shear.md` and `V3_STEP_PLAN.md:208–250` (S2 through O-02).

## Rules

1. **Rest-on test.** A record-derived entry needs a kept result that rests on the object
   (`SUBSTRATE_REQUIREMENTS.md:291–293`: a requirement is not everything a later step would find useful).
   - **S11c's kept results** are the items its closeout marks CONDITIONAL. Items it marks UNRESOLVED or OPEN
     are under "Claims future work must not inherit". They source no entry. List each in the STOP report,
     with its mark and the owner the closeout names.
   - **Numeric inputs** that the S11c calculations used and no kept result rests on get no entry. List them
     once, at the end of the pass-2 section, in a table titled "Inputs a future nonuniform calculation must
     define". Give each the owner the closeout names, or "no owner named". Cite where each is recorded; do
     not copy values.
2. **Name each object as its record names it.** Where a record says a coefficient is one basis's
   representative, name the basis-invariant object (for example `S11b_interface_coupling_law.md:76–87`).
3. **Merge, don't duplicate.** If an existing entry already covers the object, add the new source to it.
4. **Sideways targets are fine.** Examples: S11 → S11b, S11c → S12, S11c → Q2/S22.
5. **Status** takes only the schema's values. A qualification such as CONDITIONAL goes in the entry's text.
6. **The bulk equation-of-state exponent `n`.** Check whether any kept S9, S11, S11b or S11c result rests on
   `n`.
   - **If one does:** write the entry, with that result as source and the plan's S2 as target.
   - **If none does:** write no entry. Add a short note to the pass-2 paragraph that says:
     - `n_eos = 5` is inherited. Its original light-side justification took light's phase speed to be the
       bulk sound speed: `research/1pn_optics/paper/1pn_optics.tex:720–737` and `:975–978`. The paper's
       competing constructions, and its reasons for preferring `n = 5` among them, are at `:2165–2343`.
     - Whether, and where, `n` enters light bending and delay under v3's light picture is open. Give the owner
       the records name, or write "no owner named".
     - O-02 (`V3_STEP_PLAN.md:249–250`) stays open. This pass does not decide it.
7. **Second route.** Run it only on a prior result that a record in the reading list identifies as reproduced
   by the sector, with the objection taken from that record or from a source the record cites.
   - A framing or exploratory document may supply a candidate only. That includes the closeout's "Exploratory
     notes" and `docs/s11_maccullagh_differentiation.md`. List candidates in the STOP report, never as
     entries.
   - For a step where no record identifies a reproduced prior result, report "no identified prior-art
     route".
8. **Update the file's status line and its pass-2 paragraph** to say what was done. Change the STATUS.md row
   that points to `SUBSTRATE_REQUIREMENTS` only if its wording becomes wrong.
9. **Scope:**
   - Edit only `SUBSTRATE_REQUIREMENTS.md`, and STATUS.md where rule 8 requires it.
   - No calculations, no new physics, and no reviews by you.
   - Where records disagree about what is required, quote both sides in the STOP report. Do not resolve it.
   - One commit. Do not push.

## STOP report (at most one page)

- Entry counts before and after.
- New entries, by target step.
- Existing entries that gained a source.
- Size of the "Inputs a future nonuniform calculation must define" table.
- The S11c UNRESOLVED and OPEN items, with their owners.
- Rule 6: the entry, or the note.
- Rule 7, per step: entries, listed candidates, or "no identified prior-art route".
- Disagreements between records, both sides quoted.
- The commit ID.

Claude then runs the two review legs.
