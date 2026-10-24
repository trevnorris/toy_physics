# S16 record decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S16_record_decision_list_review_p1_lookups.md`, beside this
file (generator `_scratch/s16/gen/recdl_p1_lookups.sh`). "N-n" cites a section of that file.

## The pass

- **Artifact.** `directives/S16_record_decision_list.md`, orchestrator-written, with its setup lookups
  `S16_record_setup_lookups.md`. Both matched the baseline after the pass. The working tree was unchanged (N1).
- **Role.** A pre-builder decision list for a physics-bearing record (G2: one Codex + Grok pass, fold once). After
  Codex writes the record, fresh Claude + Grok review it until clear.
- **Prompt and legs.** `_scratch/s16/rec/recdl_review_prompt_p1.md` (N1). Both legs received it identically and ran in
  parallel:
  - Codex (`gpt-6.1-sol`, xhigh) wrote `/var/projects/_review_scratch/s16-recdl-p1-codex/REVIEW.md`;
  - Grok (`grok-4.7`, high) wrote `/var/projects/_review_scratch/s16-recdl-p1-grok/REVIEW.md`.

  Both exited 0 (N1).
- **Verdicts.** Codex: CLEAR. Grok: CLEAR (N2).

## Findings

None.

## Notes

- **Codex's reading order.** Codex reports that its first command printed the whole list before it had finished
  reading the sources. It then completed the source reading and saved a source checklist before evaluating the list
  (N3). Reading order is a method instruction, not a control (`.claude/skills/review-legs/SKILL.md`), so the verdict
  stands with that limitation recorded.
- **What Grok checked beyond the list's own wording.**
  - The comparator output's row classes are each named by item 2, item 3 or the STOP report. They include the
    unjoined reasons and the bindings below their parsed count.
  - The stale S14a part 2 statements that wrap across lines, which a single-line search splits, are covered by
    item 10's "every statement" rule; L10 does not claim to list them all.
  - Item 9's merge rule leaves to the register what `R-S8-06` says about the inertial response.

## Outcome

Both legs are clear, so there is nothing to fold. The list is committed as reviewed and cleared. The record author
starts from that commit.
