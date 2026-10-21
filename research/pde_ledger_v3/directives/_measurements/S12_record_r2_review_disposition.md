# S12 record, scoped review round 2: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_record_r2_review_lookups.md`, beside this file (generator
`_scratch/s12/gen/rec_r2_lookups.sh`). "K-n" cites a section of that file.

## The round

- **Artifact.** The same author's repair of the round-1 record (`4c37473c`), under round 1's dispositions and the
  user's decisions ("Keep S12"; "Fix, scoped review"; `S12_record_r1_review_disposition.md`, committed `cbdb3aa5`).
  It is preserved unchanged in the commit that adds this file. The six files matched the round-2 baseline after the
  round, and the working tree was unchanged (K1).
- **Prompt and legs.** `_scratch/s12/rec/s12_record_review_prompt_r2.md`, sha256 `dc2a40a1…` (K1). The round was
  scoped to:
  - the repair;
  - every owner statement against its source;
  - nothing else having moved since `4c37473c`;
  - the generator rerun.

  Both legs received it identically, in parallel. A fresh Claude agent (Opus) wrote
  `/var/projects/_review_scratch/s12-record-r2-claude/REVIEW.md`, and Grok (`grok-4.7`, effort high) wrote
  `/var/projects/_review_scratch/s12-record-r2-grok/REVIEW.md`.
- **Verdicts.** Claude: NEEDS REVISION, 2 findings. Grok: CLEAR (K2).
- **What both legs found sound:**
  - F1–F3 are repaired under "Keep S12". Every owner the sources and the user gave is restored, with differences
    reported as disagreements.
  - Nothing else moved since `4c37473c`. The round-0 repairs still hold.
  - The generator rerun is byte-identical.

## Findings

**H1. The register hands the S12 requirements on without their branch condition.** (Claude 1.) **CONFIRMED.**
- **The record** states the condition: "The printed objects are on the real-fraction branch until S1's gate decides"
  (K3, §6 l.287).
- **The sources** (K3):
  - The step's decision list: "The committed law is written for a real ordered fraction; a phase-carrying branch
    would add content not computed here".
  - The plan: "A13 (the conversion balance must use the branch chosen at S1)".
- **The register.**
  - The S12 pass says the requirements "constrain that future delivery". Under "Keep S12", that delivery comes
    later, possibly after S1's gate decides.
  - Neither the pass nor `R-S12-02` contains "real-fraction", "phase-carrying" or "branch". "A13" appears only once
    in the pass, as the owner disagreement (K3).
- **Origin.** The omission predates this round. Under the user's decision it now governs a future delivery.
- **What must be true.** Wherever the register hands the S12 requirements on (the S12 pass and `R-S12-02`), it
  states that they rest on the real-fraction branch and that A13 stays OPEN, as the sources give it.

**H2. An edit to the accepted light-closing pass is not reported.** (Claude 2.) **CONFIRMED.**
- **The edit.** At `e070c7f3` the pass's bullet read "the first two belong to S12 (`R-S12-02`) and to O2's carried
  handoffs". It now splits the material-destination question between S12 (`R-S12-02`) and `R-S21-01`, as the oracle
  record l.166 does (K4).
- **How it is described.**
  - The S12 pass describes its light-closing edits only as "The light-closeout status lines now record its
    acceptance".
  - The record says "The earlier light-closeout acceptance edits at `e070c7f3` are retained".
  - The round-0 record's "No other scientific record is upgraded" has been absent since round 1 (K4).
- **Origin.** The edit came in round 1 and was kept in round 2. The new text matches the oracle record. The defect is
  that the change to an accepted pass is not reported.
- **What must be true.** The record and the S12 pass report every edit made to an accepted pass. For this one, they
  quote both the accepted text and the new text, give the source that settles it, and cite a retrieval.

## The leg disagreement

Grok found the register sound. K3 and K4 show the missing branch condition and the unreported edit, so the findings
stand.

## Outcome

Two findings, both confirmed, both text in the register's S12 hand-on and in the record's account of its own edits.
Neither changes a computed object, a cross-engine claim or an owner. The reviewed version is preserved as-is in this
commit. The user decides how to finish.
