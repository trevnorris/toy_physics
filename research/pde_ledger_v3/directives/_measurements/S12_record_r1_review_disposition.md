# S12 record, review round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_record_r1_review_lookups.md`, beside this file (generator
`_scratch/s12/gen/rec_r1_lookups.sh`). "K-n" cites a section of that file.

## The round

- **Artifact.** The same author's repair of the round-0 record (`8608bc71`), under round 0's dispositions
  (`S12_record_r0_review_disposition.md`). It is preserved unchanged in the commit that adds this file. The six
  files matched the round-1 baseline after the round, and the working tree was unchanged (K1).
- **Prompt and legs.** `_scratch/s12/rec/s12_record_review_prompt_r1.md`, sha256 `a42166c7…` (K1). It asked for a
  full review of the revised record. Neither leg was handed round 0's findings. Both received the prompt
  identically, in parallel:
  - a fresh Claude agent (Opus), report `/var/projects/_review_scratch/s12-record-r1-claude/REVIEW.md`;
  - Grok (`grok-4.7`, effort high), report `/var/projects/_review_scratch/s12-record-r1-grok/REVIEW.md`.
- **Verdicts.** Claude: NEEDS REVISION, 3 findings. Grok: CLEAR (K2).
- **Both legs** found the following sound, each from its own enumeration and a byte-identical rerun of the generator:
  - the comparison account;
  - the round-0 C2 and G1 repairs;
  - the domain objects and the premises' standings;
  - the status lines.

## Findings

All three concern owners. The first two come from the round-0 repair, which removed S12 as the deliverer more widely
than the sources support. Round 0's C1 "What must be true" invited this. It said no obligation S12 does not deliver
may show S12 as its deliverer, and it did not keep the owners that the sources and the user give.

**F1. The A13 owner chain no longer matches its sources.** (Claude F1.) **CONFIRMED.**
- **The record.** It now reads "S1 gate, propagated through S5; O2's branch row also names S12 as a consumer". At
  `8608bc71` it read "S1 gate, propagated through S5/S12" (K3).
- **The sources** (K3):
  - O2's shared physics: "S1's A13 gate, propagated through S5/S12" and "Owners: S1, S5 and S12". The file never
    uses the word "consumer".
  - The O2 input contract: "propagated at S5 and S12".
  - The O2 record: "propagated S5/S12".
  - The plan: "the conversion balance must use the branch chosen at S1".
  - The step's decision list retains A13 as OPEN.
- **What must be true.** The A13 row gives the owners the sources give, quoted. It states that the printed objects
  are on the real-fraction branch until S1's gate decides.

**F2. Owners the sources and the user named were changed without being reported as disagreements.** (Claude F2.)
**CONFIRMED.**
- **The S9b oracle record**, accepted, gives "The drain's source and boundary data belong to S12 (R-S12-02)" and "is
  again S12's object (R-S12-02)" (K4). The record reports these handoffs with no live owner named and lists no
  oracle-record disagreement.
- **The register's `R-S22-01`.**
  - At `560c6c51`, its target named S12, for the throat's order-conversion source and its boundary data. It added:
    "These owners were chosen by the user on 2026-10-19" (K4).
  - The repair removed S12 from that target. The register's text says the assignment is "retained as historical
    provenance" and that who should deliver is "a question for the user" (K4).
  - The record never mentions the user's assignment (count 0, K4). Its row gives "no live owner named … formerly
    assigned to S12".
- **What must be true.**
  - No owner the sources or the user named is replaced without the user's decision.
  - Every place the record's owners differ from a source's is reported as a disagreement, both sides quoted.

**F3. `R-S21-01`'s by-target row and its entry disagree on what is owed.** (Claude F3.) **CONFIRMED.**
- **The row.** The by-target row now adds "comparison of force sectors' flow requirements" (K5).
- **The entry.** Its requirement field does not carry it (count 0, K5). The entry says the framing provenance "adds
  no claim of cross-sector compatibility".
- **What must be true.** The row and the entry say the same thing about what `R-S21-01` owes, with its source.

## The leg disagreement

Grok found the owner statements sound. The source and register lines in K3–K4 show the replaced owners and the
unreported user assignment, so the findings stand.

## Outcome

Three findings, all confirmed. None changes a computed object or a cross-engine claim. Together they turn on one
question only the user can answer: once S12's record states the requirements and leaves the forms open, who
delivers the native conversion forms, the boundary data and the throat's conversion source? Earlier the sources and
the user named S12. This is the record's second review round with findings, so the user decides how to finish, and
the ownership question goes to the user with it. The reviewed version is preserved as-is in this commit.
