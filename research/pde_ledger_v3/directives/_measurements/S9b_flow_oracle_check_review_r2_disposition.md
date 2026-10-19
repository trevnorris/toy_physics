# S9b flow-only oracle check, review round 2 (scoped): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S9b_flow_oracle_check_review_r2_lookups.md` beside this file
(generator `_scratch/s9b_oracle/gen/check_r2_lookups.sh`).

## The round
- **Scope.** The user chose to finish round 1's finding with the same author and a two-leg review scoped to the
  repair. While the repair was in progress, the user made three more choices: bounded zero tests, then sampled
  numeric tests, and finally exact zero tests with the builder's Hamilton–Lisle FORM control dropped. The briefs are
  `_scratch/s9b_oracle/build/resume_oracle_check_{3..7}.md`; the last one, `resume_oracle_check_7.md`, governs the
  reviewed version.
- **Artifacts.** The check after its second repair, by the same Codex author: eleven files in the baseline
  `_scratch/s9b_oracle/check_review/s9b_oracle_check_review_baseline_r2.sha256`, all OK at review and at acceptance.
  Compared with round 1, the value checker `S9b_flow_oracle_check_values.wl` and its transcript
  `scripts/out/S9b_flow_oracle_check_values.out` are new. The Hamilton–Lisle control transcript
  `scripts/out/S9b_flow_oracle_check_hl_spatial_form.out` is removed. The transcripts come from guard runs `run_032`
  (baseline), `run_033` (Visser control) and `run_034` (value checker). The reviewed previous version is `2caec139`.
- **Prompt and legs.** The prompt is `_scratch/s9b_oracle/check_review/s9b_oracle_check_review_prompt_r2.md`
  (sha256 `d25f4f4d…`). It states the repair as properties Z1–Z4 and asks for the legs' own census, sampling, value
  comparison against `2caec139`, and ablations of the zero test. Both legs got the identical prompt: first a fresh
  Claude agent (opus), then Grok (`grok-4.7`, effort high). Their scratch is under
  `/var/projects/_review_scratch/s9b-oracle-check-r2-{claude,grok}/`, hard-linked at
  `_scratch/s9b_oracle/check_review/r2_{claude,grok}_evidence/`. Reports: `…/check_review/r2_claude.md` and
  `…/r2_grok.txt`. **Claude: CLEAR. Grok: CLEAR.** Neither leg reports a finding.

## Findings

None from either leg.

**Sound, both legs, each with its own scripts and stdout.**
- **Z1.** Each of the 2806 residuals has its printed assumptions just before it and its zero test just after, with
  matching tags and domain ids.
  - Both legs' censuses give the same four classes:
    - 253 residuals print `0`;
    - 108 print in a nonzero form with a `True` test;
    - 40 have a `False` test;
    - 2405 have a test that reduces to an equation.
  - Neither leg finds a test that failed or did not evaluate.
  - Claude samples all 2806 residuals at points that satisfy their printed domains. Its sampler finds no anomaly
    and a witness for every domain.
  - Grok re-runs 31 of the tests under the printed assumptions, and every one matches. Grok also shows that `False`
    under those assumptions simplifies to `False`.
  - Claude's SymPy check covers 44 residuals across the four classes.
- **Z2.** Operands, assumptions and residual values are byte-identical to `2caec139`: Claude compares the 2806
  triples field by field, and Grok aligns 9097 lines with no mismatch. With the zero-test lines dropped, the only
  other difference is the `PROVENANCE` line. Both legs reproduce the baseline from the final code. The zero tests
  are identical to those of `run_019`.
- **The zero test measures.** Both legs ablate it three ways: an operand shift, removal of a needed assumption, and
  a FORM change. Claude's FORM change is to the Hamilton–Lisle (2) flow profile; Grok's is to the square-root
  velocity factor in an exchange operand. Each ablation moves the affected tests: a test goes from `True` to `False`,
  or from `True` to a reduced equation.
- **Z3.** The value checker reports deltas of values. Claude rewrites one control residual as a textually different
  but equal zero, and every delta prints `0`. Grok's independent recomputation of the 18 deltas matches the printed
  deltas and zero tests. The Visser control transcript comes from the final code.
- **Z4.** The driver's `--control` choices no longer include the Hamilton–Lisle control. The Wolfram script no
  longer rewrites the Hamilton–Lisle line element. Its transcript is removed. The controls record cites `run_023` and
  `run_025`, and both legs confirm what it says about them: exit 1, the 50 GiB cap, `No more memory available.`, and
  the recorded peaks. No transcript stands in for that control, and no bound is claimed to contain a zero test.
- **Nothing else changed.** Every hunk in the diff against `2caec139` belongs to Z1–Z4 or to the regenerated
  output.

**Notes from the Claude leg that it says change no claim.** I agree with each:
- **The value checker compares only `EXCHANGE_BALANCE_OPERATOR_RESIDUAL`.** The Visser control also moves the
  zero tests of 8 `EXCHANGE_ORACLE_OPERATOR_SOURCE_RESIDUAL` lines. The controls record says nothing of those lines,
  so nothing is misreported; they are only not counted.
- **The optical comparison has no FORM control in the accepted version.** This follows from the user's choice, and
  no record claims one. The optical residual values are byte-identical to `2caec139`. Round 1 reviewed those values,
  including Grok's FORM change to Hamilton–Lisle (1), which moved 64 of the 256 optical residual lines. Their zero tests use the same function as
  every other residual, and the ablations above show that function responding.
- **The 2405 conditional tests rest on measurement.** That they are correct is shown by the legs' sampling and
  re-runs, not guaranteed by the simplifier.

**Clear.** Nothing outstanding changes what the check computes or what may be claimed from its output. The check is
accepted at this version.
