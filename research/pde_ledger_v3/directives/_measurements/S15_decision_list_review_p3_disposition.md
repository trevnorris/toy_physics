# S15 decision list, review pass 3: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_decision_list_review_p3_lookups.md`, beside this file
(generator `_scratch/s15/gen/s15_dl_p3_lookups.sh`). "R-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout; the legs' scripts were not re-run.

## The pass

- **Artifact.** Revision 3 of `directives/S15_decision_list.md` (sha256 `ccabfa67…`), committed at `e393430a`, with
  its setup lookups and the pass 2 disposition and lookups. All four matched the baseline after the pass. The working
  tree was unchanged apart from the lookups file (R0).
- **Role.** The spec content of a pre-builder decision list, reviewed until clear (spec row). This review is scoped
  to the pass 2 repair, as the user chose: "Scoped pass 3 (Recommended)".
- **Prompt and legs.** `_scratch/s15/dl/dl_review_prompt_p3.md`, sha256 `b144f964…` (R0). Both legs received it
  identically, in parallel:
  - Codex (`gpt-6.1-sol`, xhigh): `/var/projects/_review_scratch/s15-dl-p3-codex/REVIEW.md`, exit 0;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-dl-p3-grok/REVIEW.md`, exit 0.
- **Verdicts.** Codex: CLEAR. Grok: CLEAR (R1).

## The repairs, as the legs checked them

- **D6 holds.** The filing line cites the step as "S15" with the plan's lines; "attract" has 0 matches (R2). Grok's
  text audit finds "inverse", "square" and "pull" only in the record section. Outside it, "expect" and "outcome" occur
  only in the list's rule that no name or acceptance criterion carries a value or an outcome (`:126–127`) (R2).
- **D7 holds.** Both legs built each item-1 piece from the supplied equations and listed the live quantities it
  contains (R3):
  - Grok: the pressure pieces contain `p_br`, and the first-order expansion contains one derivative of it, so a
    condition on `p_br` is one on `c_comp`. The viscous piece contains `η` and `ζ`. The transport pieces contain the
    motion of `X_2` and the radius and its rate. Each is named in the repaired sentences, or covered by "every live
    quantity the piece contains".
  - Codex: the same membership, by piece, with the sensitivities of the pressure and viscous pieces to `p_br`, `η` and
    `ζ` printed.
- **Nothing else changed.** The sections from "Supplied equations" through "What the build computes" are identical
  between revisions 2 and 3 (R4).

## Disposition: cleared

Nothing outstanding changes what the build computes or what the record may claim. Revision 3 of the list is cleared
as governing for the S15 build directive, at the commit that adds this file.

## The list's review history

| Pass | Scope | Codex | Grok | Preserved / repaired |
|---|---|---|---|---|
| 1 | G2 decision-list pass, revision 1 | 3 findings | 2 findings | D1–D5; preserved `d6e4f74a`; G2 fold `cd853e16` |
| 2 | spec content, full, revision 2 | 1 finding | 1 finding | D6–D7, dispositions `a674659c`; repair `e393430a` |
| 3 | spec content, scoped, revision 3 | CLEAR | CLEAR | cleared |

The user's three decisions of 2026-10-24 are quoted verbatim in the pass 1 disposition ("The user's decisions").

## Next

A build directive for Codex, derived from this list. The directive carries the builder's side only. The record's
side, which covers the user's decisions, the inventories, the requirements to file and the prior-art pointers, stays
out of it.
