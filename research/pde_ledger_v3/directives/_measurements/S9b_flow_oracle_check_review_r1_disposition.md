# S9b flow-only oracle check, review round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S9b_flow_oracle_check_review_r1_lookups.md` beside this file
(generator `_scratch/s9b_oracle/gen/check_r1_lookups.sh`).

## The round
- **Artifacts.** The check after its first repair, by the same Codex author: ten files in the baseline
  `_scratch/s9b_oracle/check_review/s9b_oracle_check_review_baseline_r1.sha256`, all OK at review. These are the
  seven round-0 files, the new controls record `directives/_measurements/S9b_flow_oracle_check_controls.md`, and two
  control transcripts `scripts/out/S9b_flow_oracle_check_{hl_spatial_form,visser_form}.out`. The reviewed previous
  version is `0ea65625`. The repair brief, `_scratch/s9b_oracle/build/resume_oracle_check_2.md`, stated round 0's
  four findings as properties R1–R4.
- **Prompt and legs.** The prompt is `_scratch/s9b_oracle/check_review/s9b_oracle_check_review_prompt_r1.md`:
  round 0's, with an added "The repair" section listing R1–R4. The legs were a fresh Claude agent (opus) and Grok
  (`grok-4.7`, effort high), run one after the other. Scratch is under
  `/var/projects/_review_scratch/s9b-oracle-check-r1-{claude,grok}/`, hard-linked at
  `_scratch/s9b_oracle/check_review/r1_{claude,grok}_evidence/`. Reports: `…/check_review/r1_claude.md` (report and
  addendum) and `…/r1_grok.txt`. **Claude: NEEDS REVISION (one finding). Grok: CLEAR.**
- **Prompt change during the round.** The Claude leg's own value-diff job deadlocked. Its unit
  (`s11c-guard-83eae5d41f12`) sat at the guard's 32-task cap with every thread idle. I stopped the unit, and the
  guard wrote its outcome with exit 1. The check limits its own Wolfram threads (its Wolfram script, lines 9–10); the
  leg's script did not. This explains the round-0 stall too, whose unit showed 31 of 32 tasks. I added an operational
  note saying this to the prompt before Grok ran. The version the Claude leg received is kept as
  `…_prompt_r1_v0.md` (sha256 `03f0a92b…`); the final version is sha256 `c1119659…`; the diff is the one paragraph.
  I sent the Claude leg the same note and asked it to finish the two checks the deadlock had blocked. Its addendum
  reruns the final code and runs its five ablations; its verdict and finding are unchanged.

## Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Claude 1 | **The residual simplifier does not reduce residuals under the domain it prints, so residuals that vanish there print in a nonzero form.** `domainSimplify[x_, ass_] := Refine[Cancel[Together[x]], ass]` (Wolfram script line 24) is applied to every emitted residual (lines 41–45). It does not reduce radicals under the printed inequalities, and it does not use the case's equality conditions. | **Confirmed.** Eight of the 18 `EXCHANGE_BALANCE_OPERATOR_RESIDUAL` lines print a nonzero radical form, and ten print `0`. One example is line 3960, under the `RESIDUAL_ASSUMPTIONS` printed at line 3959. The legs agree on the value. The Claude leg's `scripts/residual_sweep.wl` finds every one of the eight equal to zero under its printed assumptions. Grok's own `recompute/recompute_residuals.wl` prints `OPERATOR_RESIDUAL_TOGETHER: 0` for the example, and Grok reads it as "an identity". Grok did not count this as a finding. I do, because the printed residual is the check's measurement. As printed, a reader would see a leftover `ρ_br′` term between S9b's balance and Visser (5). Seeing that it vanishes takes algebra the reviewed check does not do. The Claude leg's A3 ablation prints the same zero residual as `0` in one run and in unreduced form in another, so the printed form depends on how the expression is written, not on its value. The sweep also finds 100 `FLOW_SIGNED_RESIDUAL` lines, on the zero-locus cases, that vanish on their printed case domain but print in a nonzero form (for example, line 639). The controls record counts FORM-control sensitivity by text: "Textually differing balance residuals: 8 of 18". The record says itself that this is a text count. This changes what may be claimed. Round 0's R1 asked that every residual be simplified under the domain its operands carry; the repair printed that domain but did not reduce the residuals under it. |

**Sound, both legs, each with its own scripts and stdout.**
- **R1.** Each residual prints its operands and its assumptions: the engine's declaration and the case domain. A
  domain change in a copy of a transcript reaches the residuals (Claude A4).
- **R2.** Visser (5) is evaluated on S9b's printed velocity with `ρ_br` live, and no oracle density is substituted.
  The builder's `visser_form` control moves the comparison-2 values (Claude's value-level diff), and Grok's diff
  shows the operator and the balance residuals moving.
- **R3.** Binding labels name their source objects. A change to the PY speed-identification moves the binding and
  its dependent residuals (Claude A5).
- **R4.** The quoted (2) and (4) profiles are compared as quoted. The orientation from (1) is a separately labelled
  operand with its own residual. A change to (2) moves only the Hamilton–Lisle (2) rows (Claude A2; Grok's oracle
  profile ablation). A FORM change to (1) moves the drift (Grok).
- **The rest.** Operands are read from the `503736ff` transcripts. A missing transcript fails before any residual is
  printed. Quotes and equation numbers match the manuscripts. The integration constant, `ρ_br` and `M(r)` stay live.
  Each leg's recomputations match the printed lines. One-sided corruptions and FORM ablations move only the rows that
  depend on what was changed. No verdict token is emitted. A guarded rerun reproduces the transcript; Grok's matches
  byte for byte, and Claude's differs only in the copy's `PROVENANCE` paths.

**Not clear.** One finding changes what may be claimed. It is in the material the first repair changed, and it is an
incomplete repair of round 0's finding 1, not a new defect. The reviewed version is preserved before any further
repair.
