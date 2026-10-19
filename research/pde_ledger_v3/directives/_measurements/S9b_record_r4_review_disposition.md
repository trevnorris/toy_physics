# S9b record, review round 4 (scoped): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S9b_record_r4_review_lookups.md` beside this file (generator
`_scratch/s9b_build/gen/s9b_record_r4_lookups.sh`).

## The round
- **Artifacts.** The record, its measurements file and generator, and the S9b register pass, as repaired in round 4
  by the same author. The brief, `_scratch/s9b_build/record/resume_record_r4.md`, pointed at the round-2
  disposition's lines for the user decision instead of restating it. The repair answers the round-3 disposition's
  "What must be true". The reviewed previous version is `d877e857`. Four files in the baseline
  `_scratch/s9b_build/record/s9b_record_review_baseline_r4.sha256`, all OK at review.
- **Scope.** This is the fifth review of the record. After four full reviews, the user chose a scoped review for it,
  from the options scoped, full or stop. The legs were asked to check three things: the two required changes; that
  every new statement cites a retrieval and the generator reproduces the filed measurements; and that every hunk of
  the diff against `d877e857` belongs to those changes.
- **Prompt and legs.** The prompt is `_scratch/s9b_build/record/s9b_record_review_prompt_r4.md` (sha256 `ecd7e2e5…`),
  new for the scoped round. The legs were a fresh Claude agent (opus) and Grok (`grok-4.7`, effort high), launched
  identically, with scratch under `/var/projects/_review_scratch/s9b-record-r4-{claude,grok}/`. Reports:
  `_scratch/s9b_build/record/s9b_record_review_r4_claude.md` and `…_r4_grok.txt`. **Both returned CLEAR.**
- **Leg note.** The Claude leg reports that the harness saved one long diff output under `~/.claude/projects/`, and
  that it did not read it. It re-saved the diffs in its scratch directory. It regenerated the measurements file in an
  isolated `git clone --shared` in its scratch directory, not in the working tree.

## Findings
None from either leg.

**Observation, not counted (Claude).** The §6 A-class cell (561) lists the 31 conditional F_ rows that carry the
unresolved brane-zero value samples. It then says "this limit is not assigned to every general or restricted A row".
The leg's retrieval finds no general or R_ A row with an unresolved brane-zero value sample (0 of 44 and 0 of 30).
**No repair.** The sentence is literally true, and the cell's exhaustive list already says which rows carry the
limit, so nothing claimable changes. "Every" could read as "any" in a later revision.

**Sound, both legs, each with its own scripts and stdout.**
- **Item 1.** The record (771–778) and the register (1076–1083) cite the round-2 disposition at `4760e6e2`, lines
  19–21, which M21 retrieves verbatim. Both state the check's handling as that source states it. Their oracle
  sections are identical, and every other oracle stays "no owner named".
- **Item 2.** Both legs' own retrievals over the comparator output, written before reading the revised cells, give
  123 unresolved `brane_zero_at_point` value samples on 41 rows:
  - 105 on the 31 conditional F_ A rows and the 4 conditional F_ B_RESIDUAL rows (residual not formed);
  - 18 on the 6 branch-existence and path-traversal rows (predicate samples).

  The three §6 cells name exactly these rows, and M22 prints the same split. That split also accounts for the
  comparator r4 acceptance's count of 105, which the record quotes: that acceptance counted only the samples with a
  residual.
- **Regeneration.** Both legs reproduce the filed measurements file byte for byte (sha256 `158fd122…`).
- **Nothing else changed.** Every hunk of the diff against `d877e857` belongs to item 1, item 2, or the
  regenerated output of the changed retrievals (M10's census pin, M21, M22).

**Clear.** No finding is outstanding that changes what the record claims or hands on. This is the record's fifth
review: four full, then one scoped by the user's choice. The record, its measurements and generator, and the S9b
register pass are accepted at this version.

## Handed on from acceptance
- The river/acoustic-flow check of the flow-only condition. It runs after this acceptance, as the user decided: a
  question agreed with Codex, a Codex-written check, review under G1, then a report. No oracle result is stated
  anywhere until then.
- Then the light register pass, then S12.
