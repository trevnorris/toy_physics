# O2: input-contract authoring directive (orchestrator)

**Author:** Claude (orchestrator), 2026-10-06. **Status:** folded once after one Codex + Grok pass (`CLAUDE.md` G2;
both NEEDS REVISION: the author may no longer adopt a restricted form, the scope is row 2's, and the premise
labels are the decision list's).
The contract the author writes is then reviewed until clear by a fresh Claude agent and Grok.

## Task
Write `research/pde_ledger_v3/directives/O2_input_contract.md`: sub-step 2 of the cleared inventory
`directives/O2_steady_brane_balance_scoping.md` (`59855a38`; its §7 row 2 is the deliverable), applying the
premise decision list `directives/O2_premise_decision_list.md` (`77d2c39a`).

## The object
One input contract for O2, on v9's setting (inventory §1). It covers exactly the inputs the inventory's §7 row 2
names, plus the decision list's two assignments to this contract:
- the brane material and its branch;
- the conservative momentum and normal-stress operands, with the inertia and normal material response;
- the material-reference evolution, which is the relaxation response under premise 1;
- the geometry, including the embedding and O4;
- the density, stiffness and projection inputs;
- the bulk-traction premises;
- the energy balance premise 1 asks for.

The drive, the exchange momentum and the full face/support balance belong to sub-step 3's spec.

For each input, give:
- **Status:** recorded, supplied, postulated, adopted premise (user, 2026-10-06), or OPEN.
- **Form:** a supplied or recorded relation as an equation, with its source lines and verified domain. An OPEN
  input is a named operand. Its arguments are not restricted beyond what a record supplies.
- **Label:** premises 1–4 carry the decision list's label (lines 13–14), with its date.
- **Owner:** where the inventory or plan assigns one.

## What binds the author
- The decision list's premises, its OPEN list and its intended claim. Apply them; do not re-decide them. This
  includes O4 and the grades as that list records them.
- The underlying records, not the inventory's summary of them. Cite source lines from the records.
- **Do not invent a law, and add no premise.** A form not supplied by a record stays a general named operand. If
  a further restriction seems needed, report it as a question for the user. Do not adopt it.
- Missing S1.5, S8 or material content is returned to its owner or kept as an operand. This is not a substrate
  reconstruction or an S11c repair.

## Constraints
- ⛔ No derivation and no CAS.
- ⛔ No expected value, sign, cancellation or outcome (`CLAUDE.md` M2).
- Keep every varying quantity live (M3). A relation that holds only with a quantity frozen is reported, not
  adopted.
- Prior art is not a premise (M3).
- Physics only.

Commit nothing. Stop after writing the file. The final message is a short summary.
