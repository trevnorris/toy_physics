# S12 comparator build, review round 2 (scoped, after repair 2): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_comparator_build_review_r2_lookups.md`, beside this file
(generator `_scratch/s12/gen/cmp_build_r2_lookups.sh`). "H-n" cites a section of that file.

## The round

- **Why.** Round 1 confirmed two findings, G1 and G2, both in label rows (`S12_comparator_build_review_r1_disposition.md`).
  As the comparator's third review round, how to finish was raised with the user, who chose: "Fix labels, scoped
  review". The same builder repaired G1 and G2. Round 1 is preserved at `e4dae246`.
- **Artifact.** The comparator after repair 2, uncommitted at review on top of `e4dae246` (H1):
  - `scripts/S12_cross_engine_comparator.py` and `scripts/test_S12_cross_engine_comparator.py`;
  - `directives/_measurements/S12_comparator_run_record.md` and `S12_comparator_builder_report.md`;
  - the output in `_scratch/s12/cmp_work/` (`r2_*`), which is scratch and is not committed.

  The four repository files are committed unchanged with this file. All 57 baseline hashes matched after both legs,
  and the working tree was unchanged apart from this file and the lookups (H1).
- **Prompt and legs.** `_scratch/s12/cmp/review/s12_cmp_build_review_prompt_r2.md`, sha256 `ccf82b6a…` (H1). It asked
  about the repair, the label path, whether anything else moved since round 1, the label controls and
  reproducibility. Both legs received it identically and ran in parallel, each under a 25 GiB pooled cap.
  - **Claude:** a fresh opus agent. Report `/var/projects/_review_scratch/s12-cmp-r2-claude/REVIEW.md`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/s12-cmp-r2-grok/grok_stdout.txt`.
- **Verdicts.** Claude: NEEDS REVISION, 1 finding. Grok: CLEAR, no findings (H2).

## What both legs found sound

- **G1 and G2 hold.** Each leg read every shared label row in both engines' notations with its own script, and
  matched the comparator on all of them.
  - A list of strings is now read string by string from both notations.
  - A notation difference alone gives no text difference, and a content difference is printed at its string.
  - The label controls use each engine's own serialization and go through the production label path. They fail when
    that path is blanked on either side, or when G1 or G2 is undone (H5).
- **Nothing else moved.** Each leg compared round 1's output with round 2's, record by record. Outside the four
  list-label rows, only the configuration record (one added key) and the resource record changed (H4). Blanking the
  label path moved no record outside the label rows (H5).
- **Reproducibility.** The run record matches the guard records, and each leg's guarded rerun reproduced the output
  apart from its resource record.

## Findings

**J1. The repaired `engine` row pairs components that the engines build as different objects.** (Claude 1.)
**CONFIRMED.**
- **The list.** Item 1: "the comparator establishes which components of the two payloads denote the same object,
  citing the line in each engine that builds them".
- **The gap.**
  - The repair gave the `engine` label two paired parts: a "runtime version" and an "engine-specific version/platform
    text" (l.315–317).
  - The cited engine lines build `(sys.version, s.__version__)` in SymPy and `{$Version, $SystemID}` in Wolfram
    (H3). So the second pair joins SymPy's library version with Wolfram's platform identifier, and the comparator
    prints a text difference between them.
  - In round 1 this row was one text comparison. The pairing is new in repair 2, in the code that repair changed.
  - Grok's parse shows both of the row's parts as unequal, without examining why (H3).
- **Effect.** The row identifies which programs ran. Its two engines' parts are not counterparts, so no comparison of
  them carries information. No value, distribution, residual or dimension row is affected, and no physics claim
  changes. What it affects is whether item 1 holds for every row.
- **What must be true.** A component pair is printed only where the cited lines build the same object. Otherwise the
  components stay printed, unpaired.

## Outcome

One finding, confirmed. It is confined to the `engine` identification row. The reviewed build is preserved unchanged
in the commit that adds this file. It is not accepted. This was the comparator's third review round, and the finding
is in code that the last repair changed. How to finish is raised with the user.
