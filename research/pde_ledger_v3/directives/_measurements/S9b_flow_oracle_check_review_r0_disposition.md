# S9b flow-only oracle check, review round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S9b_flow_oracle_check_review_r0_lookups.md` beside this file
(generator `_scratch/s9b_oracle/gen/check_r0_lookups.sh`).

## The round
- **Artifacts.** The check as first built, seven files, Codex-written to the decision list at `cd61815d`: the script
  `scripts/S9b_flow_oracle_check.py`, its Wolfram script and two data files under
  `scripts/S9b_flow_oracle_check_data/`, the transcript `scripts/out/S9b_flow_oracle_check.out`, and the run record
  and identification notes `directives/_measurements/S9b_flow_oracle_check_{run,identifications}.md`. Baseline
  `_scratch/s9b_oracle/check_review/s9b_oracle_check_review_baseline_r0.sha256`, all OK at review.
- **Prompt and legs.** The prompt is `_scratch/s9b_oracle/check_review/s9b_oracle_check_review_prompt_r0.md`
  (sha256 `09d28745…`). The legs were a fresh Claude agent (opus) and Grok (`grok-4.7`, effort high), given the same
  prompt and run one after the other, because both ablate the Wolfram script. Scratch is under
  `/var/projects/_review_scratch/s9b-oracle-check-r0-{claude,grok}/`, hard-linked at
  `_scratch/s9b_oracle/check_review/r0_{claude,grok}_evidence/`. Reports: `_scratch/s9b_oracle/check_review/r0_claude.md`
  (with the leg's addendum) and `…/r0_grok.txt`. **Both returned NEEDS REVISION.**

## Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Claude 1 | **The signed flow residuals are simplified without the engines' declared domain.** Every residual goes through `residual[a,b,ass]` with `ass = r > 0` (Wolfram script lines 17, 111, 148–159). Both engines declare `c0 > 0`, and each case attaches its own algebraic domain. The exchange block applies the case domain (lines 179–183); the flow block does not. | **Confirmed.** The PY transcript declares `c_0` positive (345 occurrences), and WL's `PROFILE_ANSATZ` domain contains `c0 > 0`. None of the 416 `FLOW_SIGNED_RESIDUAL` lines prints `0`. On the retained branch, line 114 prints a squared residual `0`, while line 123 prints a signed residual containing `c0*Sqrt[GM/c0^2]`. The leg's census (`derive/d6_domain_census.wl` and its stdout) applies each engine's declared domain and each case's attached domain, and finds printed residuals that change under them. The instrument simplifies its own residuals; a reader of the transcript should not have to do that algebra. This changes what may be claimed. |
| 2 | Claude 2 | **Comparison 2's oracle side is the right-hand side `0` of Visser (5); the operator of (5) never enters.** Wolfram script line 165: `oracleExchange=Last[First[visContinuity]]`, used at lines 169 and 186. | **Confirmed.** Output line 559 prints `EXCHANGE_ORACLE_OPERAND: 0`, and line 565 says "only zero source is compared". The leg's FORM ablation of (5)'s operator (`densityField` → `densityField^2` in line 37) moves only `PROVENANCE` and `ORACLE_TRANSCRIBED_OBJECTS`, so no comparison output depends on (5) beyond its right-hand side. The decision list names the equation (5) as the oracle operand (lines 84–90). Its residual wording, "S9b's exchange less the oracle's", allowed this reading, which is the list's part in it. As built, comparison 2 compares S9b's exchange with a literal `0`. It cannot support a statement relating S9b's balance to Visser's continuity, including about §7's density. The leg's own separate operator computation is a reviewer's computation and is not used here; the check must compute the comparison itself. |
| 3 | Claude 3 | **The PY `s9b_A` binding comes from the oracle's normalized local metric, not from S9b.** Wolfram script line 99 binds `py`s9b$A->baseHL[[2,2]]`, with `baseHL=localMetric[lineMetric[hlLine],1]` (line 68). | **Confirmed, label.** Output line 1238 prints `py`s9b$A -> 1`. The map entry `PY_dispersion_chi_A_U` names the spatial metric as A's binding. The leg's `derive/d8_localmetric_identity.py` shows that `localMetric`'s `[[2,2]]` entry equals its factor argument for a general metric. So the bound value is fixed by the normalization, not by S9b's spatial metric or by the oracle's `g_rr`. The leg reports that no printed residual changes because of this. The defect is the stated source of the binding, which changes what may be claimed about the identification a residual rests on. |
| 4 | Grok 1 | **Hamilton–Lisle's signed operand is (1)'s coordinate drift, attributed to (2).** Signed operands are `drift = −g_tr/g_rr` of (1) after the clock pullback (Wolfram script lines 64, 139, 145). They are printed as `{"HL2", …}` (output lines 48, 78), and the flow map prints `"orientation" -> "HL2"`. | **Confirmed.** In the manuscript, (1) is `−dt_ff² + (dr + v_ff dt_ff)² + …`. (2) and (4) give `v_ff` as a non-negative root, and (4) is introduced as an "infall velocity". The orientation of the signed operand therefore comes from (1), but the output attributes it to (2). The decision list names (2) and (4) as comparison 1's oracle operands. The identification map prints HL's velocity-and-frame entry as `UNMATCHED`. Visser's orientations are printed both ways; Hamilton–Lisle gets one. This changes what may be claimed about which oracle object the signed comparison uses, and it prints one orientation where the frame identification is open (list item 3). |

**Observations, not counted.**
- **A stalled ablation.** The Claude leg's FORM ablation of the spatial part of Hamilton–Lisle (1), run in its own copy,
  stopped printing after transcript line 107. Its kernel then sat idle: 4 s of CPU in 49 minutes, with stdout
  unchanged since 09:12:59. I stopped the unit (`s11c-guard-c84ad906def9`) at 10:03. The guard wrote its outcome with
  exit 1 and released the reservation. The leg does not use that run as a control. The driver reads the kernel's
  output line by line (`S9b_flow_oracle_check.py` lines 293–300), so this was not a full pipe. The cause is not
  established, and it is passed to the builder as information. The production run is unaffected.
- **Rerun from a clean checkout.** The check reads its inputs from the gitignored `_scratch/s9b_oracle_check/`
  (script line 15). Both legs reran path-retargeted copies and reproduced the transcript apart from `PROVENANCE`.
  This changes neither what is computed nor what may be claimed.
- **Leg procedure.** The Claude leg ran one 0.5 s SymPy job under the guard while its F4 run was live, against the
  prompt's one-job-at-a-time rule. Only one kernel ran at a time. No finding depends on it.

**Sound, both legs, each with its own scripts and stdout.**
- S9b operands are read from the `503736ff` transcripts. A missing tag or transcript raises. The WL inactive integral
  is activated with its printed integrand and limits.
- The integration constant and `ρ_br` stay live. No oracle density is substituted. Hamilton–Lisle continuity is
  `UNMATCHED`.
- The oracle quotes match the manuscripts in equation number, section and text, and the Wolfram transcription adds no
  term.
- Each leg's recomputed residuals match the printed lines.
- One-sided corruptions of each transcript and of an oracle equation move only their own blocks. FORM ablations of a
  PY profile, a WL profile, an oracle profile and WL's metric move the printed output.
- Operands, map and residuals are printed before the guards. No verdict token is emitted. The run record's command,
  log, cap and transcript sha256 match the guard record.

**Not clear.** Four findings change what is computed or what may be claimed. The reviewed version is preserved before
the repair. The repair goes to the same author (first repair round) as statements of what must be true, and the
repaired check gets its own review.
