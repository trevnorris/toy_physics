# O2: the steady brane balance under drain flow, scoping directive (orchestrator)

**Author:** Claude (orchestrator), 2026-10-06. **Status:** folded once after one Codex + Grok pass (`CLAUDE.md`
G2; both NEEDS REVISION: an outcome sentence removed, and the object pointed at v9's O2).

## Why
S9b spec v9, preserved as `05b1a5d5` (`research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md` at that commit,
with `directives/S9b_linked_brane_sources.md`), lists the live steady momentum and support balance as an OPEN
operand, O2 (`ℬ_hold^live`, at `:287–292` of that file). The user decided on 2026-10-06 to derive this law
before continuing. This directive asks for an inventory of what that derivation needs.

## Task: an inventory, not a derivation
Write `research/pde_ledger_v3/directives/O2_steady_brane_balance_scoping.md`, covering:

1. **The object.** O2 exactly as v9 names it (`:287–292`), on v9's setting (`:113–114`: far field, steady, `V`
   live). Keep the driving force `F_drive`, its coupling to the independent `GM`, and the exchange momentum
   `Π_n` (O3) unsupplied unless a record supplies them. Classify each of O1 and O3–O7 separately: is it an
   operand of O2, an input to it, or independent of it? Name the object, not a method for obtaining it.
2. **What the records supply.** List every verified ingredient (mass balance, momentum balance, normal force
   balance, stress or constitutive law, interface traction, bulk response, exchange momentum), each with its
   source lines and the domain it was verified on (uniform, perturbative, held, static or live).
   - Sources: `research/pde_ledger_v3/` (steps, directives, `V3_STEP_PLAN.md`, `SUBSTRATE_REQUIREMENTS.md`),
     S9b v9 and `directives/S9b_linked_brane_sources.md`.
   - The v2 tree, read-only, where the v3 plan cites it.
3. **What is missing.** Each missing ingredient, and whether it could follow from supplied ingredients or needs
   a new premise.
4. **Where the plan files it.** Find the plan step that owns this law, if any (`V3_STEP_PLAN.md`), and say
   whether this work belongs there.
5. **Prior art, as an oracle only** (`CLAUDE.md` M3). Name published results that a derived law could be checked
   against, with citations. Mark any source you did not open as unverified. ⛔ No prior-art result is a premise.
6. **Decisions for the user.** Each premise only the user can supply, stated as a question with the options.
7. **Proposed sub-steps.** Ordered, each with its deliverable and its review gate.

## Constraints
- ⛔ No derivation and no CAS. ⛔ No expected value, sign, cancellation or outcome (M2).
- Keep every varying quantity live (M3). A relation that holds only with something frozen is reported, not
  adopted.
- Physics only.

Commit nothing. Stop after writing the file. The final message is a short summary.
