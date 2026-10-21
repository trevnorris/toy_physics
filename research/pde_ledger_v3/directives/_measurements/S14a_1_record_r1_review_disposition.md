# S14a part 1 record, review round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_1_record_r1_review_lookups.md`, beside this file
(generator `_scratch/s14a/gen/rec_r1_lookups.sh`). "L-n" cites a section of that file. Every computed result quoted
here is a leg's saved stdout.

## The round

- **Artifact.** The record after repair 1 of round 0's findings T1–T4 (`S14a_1_record_r0_review_disposition.md`; the
  reviewed round-0 version is preserved at `25cb3a86`). The repair went to the same author in a resumed session,
  which exited 0 (331,514 tokens for the session, L1). It changed the record, the generator, the measurements file
  and the register; the S12 record, the build's run record and `STATUS.md` are as at `25cb3a86` (L1). Reviewed as the
  seven files of the round-1 baseline, which all matched after the round. The working tree was unchanged apart from
  the lookups file itself (L1).
- **Prompt and legs.** `_scratch/s14a/rec/s14a1_record_review_prompt_r1.md`, sha256 `98897157…` (L1): item 0 checks
  T1–T4 and that nothing else moved since `25cb3a86`; items 1–10 review the whole record again. A fresh Claude agent
  (Opus) and Grok (`grok-4.7`, effort high), in parallel, with the identical prompt (G1). Reports:
  `/var/projects/_review_scratch/s14a1-record-r1-claude/REVIEW.md` and
  `/var/projects/_review_scratch/s14a1-record-r1-grok/REVIEW.md`.
- **Verdicts.** Claude: CLEAR, no findings. Grok: CLEAR, no findings (L1).

## The repair, checked

- **T1.** Each statement of the seven unresolved rows' sample check now says that only the 18 quadratic-form samples
  of 54 can move when a derivative is taken at another point; for b173, 18 of its 36 quadratic samples lie on the
  unresolved left side. The register's S14a pass says the same. The reused sentence "The affine and cubic forms are
  unchanged" is now read as describing the comparator's repair, not the samples (L2).
- **T2.** "The bound `0≤f≤1` enters no computed object" (record l.59, L2).
- **T3.** "Commit nothing" occurs 0 times in the record (L2).
- **T4.** "material-loss" occurs 0 times. The headline reads "the two channels by which material leaves the slab"
  (record l.10, L2).
- **Nothing else moved.** Both legs found that `git diff 25cb3a86` touches only what T1–T4 and the record's account of
  edits need. Both legs reran the generator into scratch and matched the filed measurements file byte for byte (L2).

## Notes below the filter, carried open

None changes what may be claimed or what is handed on. Both are of the kind S12's J2 was.
- **The l.59 citation.** The sentence "The no-concentrated-input restriction does not remove concentrated source
  parts" cites "(M17, T2)". Its source is the S12 decision list's item 1, retrieved in M01 (Claude note 2, L3).
- **A quote retrieved second-hand.** Record §6 quotes the earlier population directive through the accepted S12
  record and pass; M13 does not retrieve it from the directive itself (Claude note 3, L3).

The Claude leg's other note, on where `R-S12-02` names the S14a objects, finds the entry and its by-target row in
agreement (L3). Its rerun left an ignored bytecode cache under `scripts/__pycache__/`; it is not tracked (L3).

## Stale status text, carried to the next record list

As with S12, the accepted files describe this record and pass as "written for source-first review": the record
header (l.3), the register's status (l.8) and S14a part 1 pass heading (l.1483), and `STATUS.md` (l.11–12) (L4). The
next step's record list updates these lines to this acceptance, with the two notes above.

## Outcome: accepted

Both legs are clear in round 1. T1–T4 hold, and nothing outstanding changes what is computed or what may be claimed
from S14a part 1, or what the record, the register or the status lines hand on (`CLAUDE.md` G4). The record, its
measurements file and generator, the S14a part 1 register pass and the status and provenance edits are accepted as
reviewed in round 1, byte for byte as in the round-1 baseline. Two full review rounds preceded acceptance, with one
author.
