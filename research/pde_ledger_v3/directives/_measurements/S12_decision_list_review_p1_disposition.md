# S12 decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_decision_list_review_p1_lookups.md`, beside this file
(generator `_scratch/s12/gen/dl_p1_lookups.sh`). "H-n" cites a section of that file.

## The pass

- **Artifact.** `directives/S12_decision_list.md` as first written, with its lookups `S12_setup_lookups.md`. Both
  were uncommitted at review and are preserved unchanged in the commit that adds this file. Both hashes matched the
  baseline after the pass, and the repository was unchanged (H1).
- **Prompt and legs.** `_scratch/s12/dl_review_prompt_p1.md`, sha256 `509c33d8…` (H1). Both legs received it
  identically.
  - **Codex:** `gpt-6.1-sol`, effort xhigh. Report `/var/projects/_review_scratch/s12-dl-p1-codex/REVIEW.txt`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/s12-dl-p1-grok/grok_review.txt`.
- **Verdicts.** Codex: FINDINGS, 4. Grok: FINDINGS, 3. One finding coincides, giving six distinct findings (H1).
- **What holds.** Both legs report that the two balances are quoted faithfully, that their dimensions agree, and that
  the projection, the mass link, the flow and the puncture interior stay with their owners.

## Findings

**K1. The U1 encoding carries no U1 content.** (Codex 1, High.) **CONFIRMED.**
- **The list.** The order a point "would" have is a live target `F(w − ζ, W)`, and the actual order is `F` plus a
  live deviation (H2).
- **The problem.** With the deviation unrestricted, `F` plus a deviation is any order field at all, including order
  outside the slab. So the encoding imposes nothing from U1. Read the other way, as relaxation toward one target set
  by geometry alone, it excludes interior order that depends on history or strain, which U1 does not exclude.
- **Also needed.** U1's content is that order vanishes outside the slab's thickness. That has content only if the
  thickness is the extent within which shear can be held at all, set by the slab, and not simply wherever order
  happens to be nonzero.
- **Repair.** State U1 as an admissibility condition on `χ_B`. Keep the interior order live, history and strain
  included. Drop the target.

**K2. A vanishing deviation does not give a source limit.** (Codex 2, Medium.) **CONFIRMED.**
- **The list.** The case where order follows the slab instantly "is reached as a limit", and the deviation-free
  object "follows as the limit case" (H3).
- **The problem.** The balance contains derivatives of the order. A deviation that tends to zero does not control
  its derivatives, so no unique source limit follows from order convergence alone.
- **Repair.** Where order is a given function of position relative to the slab, label the case a conditional
  specialization, which assumes ordering is fast. Do not call it an earned limit.

**K3. The conversion-partner inventory drops OPEN channels and merges power obligations.** (Codex 3 and Grok 1,
Medium.) **CONFIRMED.**
- **The list.** It gives O2 premise 3's carried velocity and files the energy of ordering under `R-O2-01` and
  `R-S12-01` (H4).
- **The sources.** Premise 3 "closes the transported momentum only". The additional momentum partners and their
  reaction system stay OPEN with S12 (`O2_premise_decision_list.md:32–35`). The carried total energy and the
  conversion-energy partners belong to S12 (`O2_SHARED_PHYSICS.md:280`). The bulk-direction carry is not `j_n`
  times one velocity (`:243`). `R-S12-01` is separate from S12's source partners and from `R-O2-01` (register
  `:521`). `R-O2-01` assigns its whole obligation to no one (`:568`) (H4).
- **Repair.** Point at these channels as OPEN. Keep the two power obligations separate, and assign neither to S12 as
  a whole.

**K4. The OPEN list merges O6's native identifications into the S14a bridge.** (Codex 4, Medium.) **CONFIRMED.**
- **The sources.** O2's O6 row hands S12 "dynamical conversion and distinct native source/boundary data", with
  "S14a's projected order-loss/far-field bridge" kept distinct (`O2_steady_brane_balance.md:399`). O2's S12 interface
  lists "O6 identifications" (`O2_SHARED_PHYSICS.md:605`). S14a's list is the projected source, `J`, the `J → Q` map
  and the return law (`V3_STEP_PLAN.md:636–637`) (H5).
- **Repair.** Keep the native O6 identifications as a distinct OPEN S12 handoff, inventoried and not computed here.
  S14a keeps its bridge.

**K5. The boundary inventory replaced the IR/return conditions with a galactic outer edge.** (Grok 2, Medium.)
**CONFIRMED.**
- **The list.** "The mouth, the collar, the outer edge and the `𝔅` branch selection", the outer edge being a body
  inside a galaxy's flow inside an expanding brane. The same list says galaxy-scale behaviour is not a requirement
  here (H6).
- **The sources.** The plan names "the mouth, collar and IR conditions, and the `𝔅` branch selection"
  (`V3_STEP_PLAN.md:612`). Stage 045 gives two-sided IR return conditions (`:48`). The two-throat handoff's IR/return
  scheme covers one- versus two-sided domains, the return or absorber law and far-field truncation
  (`docs/two_throat_simulation_handoff_spec.md:340`) (H6).
- **Repair.** Use the plan's inventory, with the IR/return scheme OPEN. U3 is a statement about the IR condition on
  the bulk. The galactic setting is not imposed.

**K6. The S9b/O2 domain is not a resting brane.** (Grok 3, Low.) **CONFIRMED.**
- **The sources.** The mass is at rest and its drain flows. The profiles are Eulerian steady, and lab time is in the
  brane's far-field rest frame (`O2_steady_brane_balance.md:74–75`; `S9b_light_bending.md:84`) (H7).
- **Repair.** Describe the domain as recorded, and name its restrictions (one isolated body, steady profiles, the
  far-field rest frame) as restrictions.

## Outcome

Six findings, all confirmed. The list is folded once (G2). K1 and K2 change what the build computes, so the revised
list gets its own review before any builder starts (spec row). The reviewed list is preserved as is.
