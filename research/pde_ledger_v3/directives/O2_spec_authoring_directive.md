# O2: live-balance spec authoring directive (orchestrator)

**Author:** Claude (orchestrator), 2026-10-06. **Status:** folded once after one Codex + Grok pass (`CLAUDE.md` G2;
Codex CLEAR, Grok NEEDS REVISION: inputs keep the contract's statuses, and the spec chooses no sheet/slab region).
The spec the author writes is then reviewed until clear by two non-author legs.

## Task
Write `research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md`: sub-step 3 of the inventory
`directives/O2_steady_brane_balance_scoping.md` (`59855a38`; its §7 row 3 is the deliverable). Both engines will
read this spec in sub-step 5.

Its inputs are:
- the accepted input contract `directives/O2_input_contract.md` (`217a92e9`);
- the premise decision list `directives/O2_premise_decision_list.md` (`77d2c39a`).

## The object
The live steady momentum and support balance `ℬ_hold^live`, in-plane and normal, on v9's setting (inventory §1),
with `V` and `j_n` live. It is conditional on premises 1–4, untruncated, with every OPEN operand named. The
engines construct it; the spec says what they construct and from what.

## What the spec supplies
- **Every input, as an equation or a named operand,** with the status and domain the contract gives it, and its
  source in the contract. Premises 1–4 stay adopted substrate inputs; OPEN operands stay named general unknowns;
  a recorded relation keeps its recorded domain. Label an input supplied (the build cannot test it) only where
  the contract does. The spec must stand alone, because the Wolfram engine is blind and reads nothing else.
- **The sheet/slab map.** How the OPEN `𝒥_map` (O6) enters. The spec chooses neither a sharp-sheet nor a
  finite-slab reduction, since the decision list leaves that map a general unknown.
- **How each input enters.** This is the obligation the contract's §1 leaves to this spec, so that nothing is
  counted twice or dropped. The drive, the face/support loads and O3 enter as the decision list's premises 2–3
  declare.
- **The energy pairing,** as the decision list's premise 1 requires.
- **The model point:** the regime, every recorded freeze or restriction used, and its transfer limits
  (inventory row 3).
- **The interfaces** to S12, S14a/S14, S16 and S21 (inventory row 3).
- **What each engine prints.** The balance's terms per component, each traceable to the input it came from,
  with its domain qualifications. Scripts print computed objects and state no conclusions.
- **What is deferred to the build,** stated as such.

## What binds the author
- The contract and the decision list. Apply them; do not re-decide them. Add no premise; a further restriction
  that seems needed is reported as a question for the user.
- ⛔ The spec does not write the assembled balance. A balance typed into the spec would be a hand-typed payload,
  not a computed one (`CLAUDE.md` E2).
- ⛔ No expected value, sign, cancellation or outcome, and no acceptance criterion referencing one (M2).
- Keep every varying quantity live (M3). A relation that holds only with a quantity frozen is reported, not
  adopted.
- Prior art is not a premise (M3).
- Physics only. Implementation belongs to the build directive (sub-step 4).

Commit nothing. Stop after writing the file. The final message is a short summary.
