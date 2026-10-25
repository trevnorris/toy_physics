# S15 comparator decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_comparator_decision_list_review_p1_lookups.md`, beside this
file (generator `_scratch/s15/gen/s15_cmpdl_p1_lookups.sh`). "M-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout; the legs' scripts were not re-run.

## The pass

- **Artifact.** `directives/S15_comparator_decision_list.md` and its lookups
  `directives/_measurements/S15_comparator_setup_lookups.md` (generator `_scratch/s15/gen/s15_cmp_setup_lookups.sh`),
  uncommitted on `3192123b`. Both matched the baseline after the pass, and the working tree was unchanged apart from
  this file and the lookups file (M0).
- **Gate.** G2: a pre-builder decision list, orchestrator-written, so its legs are Codex + Grok. One pass, folded once.
- **Prompt and legs.** `_scratch/s15/cmp/cmp_dl_review_prompt_p1.md`, sha256 `d4f1b5dc…` (M0). Both legs received it
  identically, in parallel:
  - Codex (`gpt-6.1-sol`, effort xhigh): `/var/projects/_review_scratch/s15-cmpdl-p1-codex/REVIEW.md`, exit 0;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-cmpdl-p1-grok/REVIEW.md`, exit 0.
- **Verdicts.** Codex: CLEAR. Grok: CLEAR (M1). Neither leg reports a finding.

## What the legs checked

- **Codex** wrote its own source-side requirements before reading the list (`PRE_ARTIFACT_REQUIREMENTS.md`), then
  cross-walked the list against the sources (`13_contract_crosswalk.py`). It confirmed that the engines, the S15
  decision list and the S16 comparator are byte-identical to the commits the list names (M2). It notes that comparing
  prose conditions as text is textual evidence only, not semantic equivalence, and that this is what the list states.
- **Grok** checked each of the seven questions against the streams:
  - every one of the seven label sets is present in both streams for every grouped piece, component and surface, and
    for the region's momentum and rate groups (`group_coverage.py`: `absent 0`);
  - the streams share 401 tag names, none of them a group tag, so a join on shared names would compare pieces and miss
    every group; the list's binding rules forbid that join (M2);
  - the 208 `dimension_argument/name` payloads are distinct, while the 27 `partition/unlabelled_coefficient` payloads
    are identical; the list identifies records by content or by their place as the source defines it, never by stream
    position alone (M2);
  - it mapped each comparator failure the prompt named to the control in item 7 that would catch it.

## Examined by Grok and not filed

- One group's rate is printed as `0` by one engine and as a Leibniz sum of integrals of `0` by the other. Item 5
  prints the residual as computed, with evaluation beside it, and item 1 reads both notations as the objects they are.
- The prose conditions state the same order in different words. Item 4 compares them as text, and the list does not
  let the record read a text residual as a failed domain or a failed order.
- `balance_parts` lists the inner surface first in both engines. Item 2 forbids pairing by position alone, and the
  inner/outer repoint control catches a comparator that follows list order.

## Decision

**Cleared at this version.** No finding, so the fold G2 allows has nothing to fold. The list governs the comparator
build as written. Next: the build brief, its leak probe and contract check, then the Codex build. The comparator is
then reviewed by fresh Claude + Grok until clear.
