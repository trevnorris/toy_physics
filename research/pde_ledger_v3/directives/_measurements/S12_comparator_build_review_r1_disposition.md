# S12 comparator build, review round 1 (after repair 1): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_comparator_build_review_r1_lookups.md`, beside this file
(generator `_scratch/s12/gen/cmp_build_r1_lookups.sh`). "E-n" cites a section of that file.

## The round

- **Why.** Round 0 confirmed four findings, F1–F4, with notes N1–N2 (`S12_comparator_build_review_r0_disposition.md`).
  The same builder repaired them. Round 0 is preserved at `f1d075b8`.
- **Artifact.** The comparator after repair 1, uncommitted at review on top of `f1d075b8` (E1):
  - `scripts/S12_cross_engine_comparator.py` and `scripts/test_S12_cross_engine_comparator.py`;
  - `directives/_measurements/S12_comparator_run_record.md` and `S12_comparator_builder_report.md`;
  - the output in `_scratch/s12/cmp_work/` (`r1_*`), which is scratch and is not committed.

  The four repository files are committed unchanged with this file. All 39 baseline hashes matched after both legs,
  and the working tree was unchanged apart from this file and the lookups (E1).
- **Prompt and legs.** `_scratch/s12/cmp/review/s12_cmp_build_review_prompt_r1.md`, sha256 `52ba2ab2…` (E1). Item 0
  asked whether each round 0 "What must be true" holds, and the method asked each leg to undo each repair alone. Both
  legs received the prompt identically and ran in parallel, each under a 25 GiB pooled cap.
  - **Claude:** a fresh opus agent. Report `/var/projects/_review_scratch/s12-cmp-r1-claude/REVIEW.md`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/s12-cmp-r1-grok/grok_stdout.txt`.
- **Verdicts.** Claude: NEEDS REVISION, 2 findings. Grok: no findings (E2).

## The round 0 repairs

Both legs found that F1–F4 and N1–N2 hold on the real streams and in the tests. Undoing each repair alone on a copy
made at least one test fail and moved the production output, in each leg. Neither leg found a repair that narrowed,
reconciled or collapsed anything outside its target.

## Findings

**G1. A label that is a list of strings is read as raw text, so a notation difference is printed as a text
difference.** (Claude 1.) **CONFIRMED.**
- **The list.** Item 4: "A label (a premise line, a scope line, a comparison note) is compared as text." The round 0
  F1 "What must be true": "Where the two engines print the same object in different notations, the reader builds the
  same object."
- **The gap.**
  - A SymPy label is read as its raw line text (l.346–347). A Wolfram label is parsed, and where it is a list, its raw
    text is used (l.1147–1148). The two are then compared character by character (l.1093–1098).
  - Both engines print some labels as a list of strings, in their own notations: `dimension_basis` and `completion`
    (E3).
  - The comparison then measures the notation, starting at the first character. Where the strings themselves differ,
    the printed first difference is the bracket, not the place where the strings differ.
  - Grok recorded two such rows as unequal as text, without examining why (E3).
- **What must be true.**
  - A label that is a list of strings is read as the same object from both notations, so a notation difference alone
    is not printed as a text difference.
  - A difference in any string's content is still printed, at that string.

**G2. The changed-label control does not go through the production label path.** (Claude 2.) **CONFIRMED.**
- **The list.** Item 6: "Each fixture is serialized in its engine's own format and goes through the same reader,
  extraction, join and comparison path as the production run", and "each control fails when the comparison step it
  targets is disabled". "A changed label" is a required mutation.
- **The gap.**
  - The control uses a tag that is not a label tag, with a quoted SymPy text that is not SymPy's own label format
    (test l.236–239). The production label branch (l.346) runs only for label tags.
  - Claude's ablation A7 blanked the text that branch reads. Production changed on every label row, and all 36 tests
    passed (E4).
- **What must be true.** A changed-label control uses each engine's own label serialization. It goes through the
  production label path, for a label tag and for a label that is a list of strings. It fails when that path is
  disabled.

## Notes

- **N4.** In the Claude leg, a guarded probe overlapped a guarded ablation job. Both were admitted under the pooled
  budget, neither was killed, and no kernel ran (E5). It has no effect on the evidence.

## On the leg without findings

Grok checked the repairs and the value, distribution and dimension rows in depth. It recorded G1's two rows as unequal
as text without examining why (E3), and it did not ablate the label path (G2). Both findings rest on the comparator's
source lines, the two transcripts' label lines, and Claude's ablation and its stdout.

## Outcome

Two findings, confirmed. Both are confined to label rows. Both legs found every value, distribution, residual and
dimension row they checked sound. The reviewed build is preserved unchanged in the commit that adds this file. It is
not accepted. This was the comparator's second review round, so a further repair and review would be the third. How to
finish is raised with the user.
