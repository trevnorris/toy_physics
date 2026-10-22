# S16 decision list, review pass 4: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S16_decision_list_review_p4_lookups.md`, beside this file
(generator `_scratch/s14a/gen/s16_dl_p4_lookups.sh`). "J-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout.

## The pass

- **Artifact.** Revision 4 of `directives/S16_decision_list.md` (sha256 `9ca7480d…`), with its setup lookups and the
  pass 3 dispositions. The list and lookups were uncommitted on `c9da00e8`. All three matched the baseline after the
  pass. The working tree was unchanged apart from the lookups file (J1).
- **Role.** The spec content of a pre-builder decision list, reviewed until clear (spec row). This review is scoped
  to the pass 3 repair, as the user chose: "Repair, then scoped review (Recommended)".
- **Prompt and legs.** `_scratch/s16/dl/dl_review_prompt_p4.md`, sha256 `294016d4…` (J1). Both legs received it
  identically, in parallel:
  - Codex (`gpt-6.1-sol`, xhigh, workspace-write, working directory its scratch):
    `/var/projects/_review_scratch/s16-dl-p4-codex/REVIEW.md`, exit 0;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s16-dl-p4-grok/REVIEW.md`, exit 0.
- **Verdicts.** Codex: CLEAR. Grok: CLEAR (J2).

## The repairs, as the legs checked them

- **F1 holds.** Both legs found the frame named and each piece one object (J3).
  - Codex: the received momentum in the equation's coordinates is `A`. The reading in coordinates moving with `X_c`
    is `A + T`, where `T` is the printed piece from the motion of `X_c`. Item 3 excludes both pieces, and the face
    exchange stays undivided.
  - Grok: two routes give the same at-rest flux and the same two pieces, with residual `0`. Neither piece is
    identically zero. With the radii fixed, the moving-frame reading equals `dP/dt`, residual `0`.
- **F2 holds.** The debt sentence names the plan's six supplied inputs and keeps the printed count, moments and
  fluxes apart: `DEBT_SENTENCE_FILES_PRINTED_AS_DEBTS False` (J3).
- **F3 holds.** The requirement states both comparisons: `COMPACT_HAS_SEPARATION True`,
  `COMPACT_HAS_FIELD_SCALE True` (J3).
- **What the repair touches.** Codex: reversing the six repair replacements reproduces revision 3's whole body.
  Grok: four hunks, with the model point, the OPEN list and the intended claim untouched. Grok also found
  `NO_MASS_STILL True`, `FREEZES_NONE True` and `ITEM3_NOT_EQUATED_TO_P True` (J3).
- **The user's decisions.** Both legs found the five quotations matching their filed sources, and Grok found the
  decision section unchanged from revision 3 (J3).

## Disposition: cleared

Nothing outstanding changes what the build computes or what the record may claim. Revision 4 of the list is cleared
as governing for the S16 build directive, at the commit that adds this file.

## The list's review history

| Pass | Scope | Codex | Grok | Preserved / repaired |
|---|---|---|---|---|
| 1 | G2 decision-list pass, revision 1 | 5 findings | 4 findings | D1–D6; preserved `ffaa2824`; user decisions on passive mass and mass; G2 fold `ef82fb3b` |
| 2 | spec content, revision 2 | 2 findings | 1 finding | E1–E3, dispositions `6ad2d816` |
| 3 | spec content, full, revision 3 | CLEAR | 3 findings | F1–F3; preserved `c9da00e8` |
| 4 | spec content, scoped, revision 4 | CLEAR | CLEAR | cleared |

The first Codex launch of pass 1 stopped before any witness ran. It is kept aside and not counted (pass 1
dispositions, N2).

## Next

A build directive for Codex, derived from this list. The directive carries the builder's side only. The record's
side, which covers the inventories, the requirements to file and the prior-art pointers, stays out of it.
