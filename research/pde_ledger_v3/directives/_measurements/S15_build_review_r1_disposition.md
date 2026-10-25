# S15 build, review round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_build_review_r1_lookups.md`, beside this file (generator
`_scratch/s15/gen/s15_build_r1_lookups.sh`). "M-n" cites a section of that file. Every computed result quoted here is
a leg's saved stdout; the legs' scripts were not re-run.

## The round

- **Artifact.** Repair 1 of the S15 build, made in the builder's session against H1
  (`S15_build_review_r0_disposition.md`), uncommitted on the reviewed baseline `efe5ce73`. The five files are
  `scripts/S15_two_drain_region.py`, `mathematica/S15_two_drain_region.wl`, their transcripts and
  `directives/_measurements/S15_build_run_record.md`. All five matched the r1 baseline after the round. The working
  tree was unchanged apart from this file and the lookups file (M0).
- **The runs.** The repair's guarded runs, `_scratch/s11c/s15/sympy_h1_01` and `wolfram_h1_01`, exited 0 with limits
  verified and empty stderr (M2).
- **Prompt and legs.** `_scratch/s15/build_review/s15_build_review_prompt_r1.md`, sha256 `b8bb2e03…` (M0). It is the
  round-0 prompt with one added item 0, the repair check. Both legs received it identically, in parallel:
  - a fresh Claude agent (Opus): `/var/projects/_review_scratch/s15-build-r1-claude/REVIEW.md`;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-build-r1-grok/REVIEW.md`, exit 0.
- **Verdicts.** Claude: CLEAR. Grok: CLEAR (M1). Neither leg reports a finding.

## H1: repaired

H1's "What must be true" (`S15_build_review_r0_disposition.md:58–61`) holds in both engines.
- **What each engine computes** (M3). Each still forms only `rate − rhs` (`.py:294–297`, `.wl:169–171`). Neither
  reduces it with the supplied local balance. Each now prints, beside each component's residual, a label saying no
  route in that engine reduces it (`.py:298`, `.wl:172`; transcript `.py.out:1294/1349/1404`, `.wl.out:1254/1308/1362`).
  `balance_construction` now lists the printed objects only (`.py:311`, `.wl:177`).
- **The round-0 sentence is gone.** Its phrase counts are 1 in each round-0 transcript and 0 in each repaired one
  (M3).
- **Nothing else changed.** Against the round-0 transcripts (sha256 `85d99aed…`, `8774462d…`), each engine removes one
  line, the old `balance_construction`, and adds four: the three reduction labels and the new `balance_construction`.
  Every other line is byte-identical (M4). No object of items 1 or 2 changed.
- **The repair check fails for the intended reason** (M5). Each leg's check fails on the repair undone, where the
  round-0 claim stays byte-identical while a FORM ablation opens the closure on the engine's own objects:
  - Claude: dropping `X_2`'s motion from `v_S` moves the closure from `-4.809e-13` to `3.799e-01` (component 1). The
    check gives `H1_CHECK FAILS` on both engines with the repair undone, and `passes` on both repaired engines, and
    under the supplied balance corrupted.
  - Grok: `CHECK: FAIL` on both round-0 transcripts (`claim_lines: 1`, `reduction_lines: 0`) and `CHECK: PASS` on
    both repaired ones (`claim_lines: 0`, `reduction_lines: 3`). Undoing the repair on a Wolfram copy and corrupting
    the stress's sign restores the round-0 claim byte for byte while the residual moves (`UNDONE_CLAIM_EQUAL_R0
    True`, `UNDONE_RESIDUAL_EQUAL_R0 False`).
- **The legs agree.** Both find that the remaining lines stay true under one-sided corruption of the supplied
  balance, while the balance operands, residual and equation move.

## What the legs checked and found sound

Both legs derived the contract's objects first, from the supplied balance, with every field live and both centres and
both radii moving.
- **Item 1.** The Claude leg evaluated every printed object on its live instance, 315 objects plus 3 residuals (SymPy)
  and 318 (Wolfram), against its own derivation. Its largest object difference is `4.1e-13` in each engine. Its
  larger residual differences, such as `-6.809e-09` in the filed tail, come from its six-digit print of its own route
  (M6). Grok's derivation
  confirms the normals, `v_S`, the four rest and two transport pieces, the groups, the expansion and the frame.
- **Items 2 to 8.** The rate, the parts and `∫f` are printed apart, `f` undivided. Dimensions are computed. The domain
  and pressure conditions are printed. The supplied forms are as stated. There is no force law and no verdict
  token, and nothing is emitted conditionally on a value. The Wolfram engine reads nothing and differs in
  construction. Both legs' guarded reruns reproduce both transcripts byte for byte: `4839b124…` and `811bea8c…` (M2,
  M6).
- **Item 9.** Every FORM ablation moved the output: the inner normal flipped (133 changed tags per engine, Grok),
  `X_2`'s motion dropped (116 per engine, Claude), two-label terms misfiled (186 group objects, SymPy, Claude), and a
  linear `p_br` (67, Wolfram, Claude).

## Notes (not findings)

- **For the record.** `region/component_i/balance_equation` (`Eq(rate, rhs)` / `lhs == rhs`) is an equation object,
  printed unreduced beside the reduction label. The record does not read it, or the residual, as the balance holding.
  Whether the parts and `∫f` make up the rate is shown only by the legs' outside computations (M5), not by the
  engines.
- **`c_comp`.** The round-0 note stands: the free symbol enters no object, and the record states `c_comp` only as
  `p_br`'s derivative.
- **The Claude leg's stopped run.** The Claude leg's own Wolfram evaluator (`claude_evalwl_full01`, unit
  `s11c-guard-293a6f9936b7`) looped endlessly on an empty `$ScriptCommandLine`. Its prompt forbids stopping guarded
  runs, so it reported the unit and exited. The orchestrator checked that the unit ran the leg's script, not an
  artifact, and stopped it. The guard wrote its outcome (exit 1, limits verified). The orchestrator deleted that run's
  2.49 GB stdout, after a listing showed the one file. The artifact baseline and status were unchanged (M7). The leg
  redid the job as `claude_evalwl_full02` and finished.

## Decision

**Accepted.** Nothing outstanding changes what either engine computes or what may be claimed from its transcripts (G4).
The S15 build is accepted at the repaired version. Next: the S15 cross-engine comparator, which starts with its own
decision list.
