# S9b record, review round 3: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S9b_record_r3_review_lookups.md` beside this file (generator
`_scratch/s9b_build/gen/s9b_record_r3_lookups.sh`).

## The round
- **Artifacts.** The record, its measurements file and generator, and the S9b register pass, as repaired in round 3
  by the same author (brief `_scratch/s9b_build/record/resume_record_r3.md`, against the round-2 disposition's "What
  must be true"). The reviewed previous version is `4760e6e2`. Four files in the baseline
  `_scratch/s9b_build/record/s9b_record_review_baseline_r3.sha256`, all OK at review.
- **Prompt and legs.** The prompt is `_scratch/s9b_build/record/s9b_record_review_prompt_r3.md` (sha256 `42dc76e2…`):
  round 2's prompt with the baseline and the previous version advanced. The legs were a fresh Claude agent (opus) and
  Grok (`grok-4.7`, effort high), launched identically, with scratch under
  `/var/projects/_review_scratch/s9b-record-r3-{claude,grok}/`. Reports:
  `_scratch/s9b_build/record/s9b_record_review_r3_claude.md` and `…_r3_grok.txt`. Both returned **NEEDS REVISION**.
- **Review scope.** This is the fourth review of the record. The user chose a full review for it, from the options
  full, scoped or stop.
- **Leg note.** The Claude leg reports that the harness saved one long `git diff` of the register under
  `~/.claude/projects/`, and that it read that file and nothing else under that path. The file held that diff only.

## Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Claude 1, Grok 1 | **The river/acoustic-flow check is handed on with "no owner named", and the user decision is stated without a source.** The record (766–772) and the register pass (1073–1079) say the 2026-10-19 decision names no owner. The written decision is the round-2 disposition at `4760e6e2`, lines 19–21: the user chose "Queue it after record", and "Codex writes the check and it is reviewed under G1". Neither the record nor the register cites that disposition, and neither retrieves it. | **Confirmed, both legs.** The lookups show the disposition's lines and zero citations of it in the three files. **The cause is my round-3 brief.** It gave the decision as a dated fact, with no source and no owner (brief lines 3–5), and asked for "its owner or 'no owner named'". The author recorded what it was handed. So this round repairs my previous round's prose, not the author's work. The other oracles stay "no owner named", and both legs agree that part is right. |
| 2 | Claude 2 | **The §6 row-class table puts the brane-zero gap on the wrong rows.** The cell for the A observables "and conditional F_FLOW/F_LIVE copies" (559) says unresolved brane-zero samples remain. The cell for B_RESIDUAL "and R_*/conditional F_* copies" (561) gives "Formed value samples, with printed differences retained", with no brane-zero limit. The leg's retrieval script places all 105 unresolved brane-zero value samples in the conditional F_ rows: 93 in F_ A rows and 12 in F_ B_RESIDUAL rows. The general and R_ rows formed all of theirs. | **Confirmed.** The lookups count the comparator's "operand evaluations unavailable for subtraction" reason in each F_ B_RESIDUAL row (3 each) and in F_FLOW_A_DEFLECTION_G020 (3). They find none in the general and R_ rows checked. So the cell at 561 omits a limit its own class carries. The cell at 559 states a limit true only of the class's F_ copies. The aggregate count of 105 is already quoted at 628–631, so nothing is lost from the totals. What is wrong is which rows the cells say carry the gap. The cell at 561 has been unchanged since round 0 (`36f59ff1`), and the four earlier full reviews passed it. |

**Raised by the Claude leg and not counted (I agree).** Register line 61 says "A postulate with a named retirement
condition generates a requirement", and the supplied local wave equation has no entry of its own. The handoff is not
lost: §7 (682) names S11c as owner, `R-S8-01`'s S9b scope says it "stays the unfinished S11c handoff", and the
round-0 legs found the dependences owned. Whether line 61 applies here stays an open question for the register; it
changes no owner and nothing handed on from S9b.

**Sound, both legs, each with its own scripts and stdout.**
- Both legs reran the generator, and it reproduces the filed measurements file byte for byte (sha256 `b1c4393a…`).
- Both legs' comparator enumerations match §6, including the six line kinds, the stratum inventory (19,534 entries,
  two empty), the coverage labels and Appendices B and C.
- Round 2's repairs hold:
  - The comparator's line-kind table and stratum inventory are reported with labels and sizes.
  - The open-oracle column is consistent between the record and the register, except the owner in finding 1.
- The single-mechanism restrictions, the item-5 count, the branch bindings, the forward case, the premises, D5's
  obligation, the dependences and the register census are as the previous rounds found them.

**Not clear.** Both findings change what the record hands on: the first changes an owner, the second where a stated
limit applies. Finding 1 repairs my brief, so this round repaired my previous round's wording rather than converging
on the artifact. Finding 2 has been in the record since round 0. The reviewed round-3 record is preserved before any
repair overwrites it.

## What must be true after the repair
1. The record and the register cite the round-2 disposition at `4760e6e2` (lines 19–21) for the user decision. Each
   carries a retrieval of those lines, and each states the check's handling as that source states it. Every other
   oracle stays as it is.
2. Each row-class cell that states the brane-zero limit names the rows that carry it, as the comparator output prints
   them.
3. Nothing else changes. The measurements file is regenerated by the final generator.
