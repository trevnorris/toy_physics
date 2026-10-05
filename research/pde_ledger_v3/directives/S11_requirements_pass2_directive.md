# S11 requirements pass (pass 2): directive for Codex

**Author:** Claude (orchestrator), 2026-10-05. **Builder:** Codex. `AGENTS.md` applies.

**Deliverable:** `SUBSTRATE_REQUIREMENTS.md` gains an entry for every substrate or sibling-step object that a
kept S11, S11b or S11c result rests on. Then every brane, bulk and interface quantity the light sector uses
has either an owner step or an explicit "no owner named".

**Why now:**
- **Pass 2 was skipped.** `SUBSTRATE_REQUIREMENTS.md:299–300` schedules pass 2 over S11 and S11b, and it never
  ran. Quantities first introduced in S11b and S11c have no entry: bulk mass density, the face law, thickness
  dynamics, `μ_S`, the energy-basis coefficients, and the drain velocity.
- **This closes S11.** It is step 1 of the light-sector close-out.
- **Classification:** this directive states no physics. It points at the register's own schema and method. The
  physics is in the entries, which get two review legs (a fresh Claude agent and Grok) at the STOP.

## Read

- **Schema and method:** `SUBSTRATE_REQUIREMENTS.md`:
  - "Entry schema";
  - the pass method under "Pass 2";
  - the "second route" (prior-art failure modes).

  Follow them as written.
- **Records:**
  - `steps/S11_stray_longitudinal.md`;
  - `steps/S11b_interface_coupling_law.md`, plus `steps/S11bA_interface_response.md` and
    `steps/S11bB_interface_assembly.md` where S11b points to them;
  - `steps/S11c_PARTIAL_CLOSEOUT.md` and the S11c a, b, c1, c2 and d records.
- **Starting checklist:** `_scratch/light_inventory/brane_bulk_inventory_2026-10-05.md`. ⛔ It is not a source.
  Every entry cites the step record, or the spec or report that the record cites.

## Rules

1. **Rest-on test.** An entry needs a kept result that rests on the object. The register says a requirement
   is not everything a later step would find useful.
   - S11c closed PARTIAL. Its kept results are the items its closeout marks CONDITIONAL, UNRESOLVED or OPEN.
     An entry whose only source is an UNRESOLVED result states that status.
   - S11c development values that no kept result rests on get **no entry**. List them once, at the end of the
     pass-2 section, in a table titled "Inputs a future nonuniform calculation must define". Give each the
     owner its records name, or write "no owner named". Do not copy their numeric values; cite where they are
     recorded.
2. **Merge, don't duplicate.** If an existing entry already covers the object, add the new source to it.
3. **Sideways targets are fine.** Examples: S11 → S11b, S11c → S12, S11c → Q2/S22.
4. **Status** is OPEN unless a record shows the object delivered.
5. **Add one specific entry** for the bulk equation-of-state exponent:
   - **Source:** S9, which replaced "light travels at the bulk sound speed" with a brane shear wave at
     `c_γ² = μ_R/ρ_br`.
   - **Targets:** the proposed light-bending requirement step, and the S1/S1.5 bulk block.
   - **Requirement:** re-establish the light-side justification of `n_eos = 5` under v3's light picture.
   - **Evidence for the original justification:** `research/1pn_optics/paper/1pn_optics.tex:731–736, 977–978`.
   - **Evidence for the "provisional" label:** `V3_STEP_PLAN.md:246–247`.
6. **Second route.** Run it for S11 (prior art on an elastic aether's longitudinal mode) and for S11b/S11c
   (prior art on interface coupling and mode conversion). Use only prior art that the records or docs already
   cite. Name the source of each new obligation.
7. **Update the file's status line and its pass-2 paragraph** to say what was done. Change the STATUS.md row
   that points to `SUBSTRATE_REQUIREMENTS` only if its wording becomes wrong.
8. **Scope:**
   - Edit only `SUBSTRATE_REQUIREMENTS.md`, and STATUS.md where rule 7 requires it.
   - No calculations, no new physics, and no reviews by you.
   - One commit. Do not push.

## STOP report (at most one page)

- Entry counts before and after.
- New entries, by target step.
- Existing entries that gained a source.
- Size of the "Inputs a future nonuniform calculation must define" table.
- Any place where records disagree about what is required. List these; do not resolve them.
- The commit ID.

Claude then runs the two review legs.
