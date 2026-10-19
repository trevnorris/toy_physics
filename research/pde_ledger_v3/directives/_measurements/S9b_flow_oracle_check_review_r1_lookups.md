# Measurements: round-1 review of the S9b flow-only oracle check (generated 2026-10-19 13:06)

Generator: `_scratch/s9b_oracle/gen/check_r1_lookups.sh` (sha256sum/grep/sed/cut/wc/cat only). The legs' scratch
directories are hard-linked at `_scratch/s9b_oracle/check_review/r1_{claude,grok}_evidence/`.

## Baseline, prompts and verdicts
```
$ (cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s9b_oracle/check_review/s9b_oracle_check_review_baseline_r1.sha256)
research/pde_ledger_v3/scripts/S9b_flow_oracle_check.py: OK
research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl: OK
research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/identification_map.json: OK
research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/oracle_quotes.json: OK
research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out: OK
research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check_hl_spatial_form.out: OK
research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check_visser_form.out: OK
research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_run.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_identifications.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_controls.md: OK
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_oracle/check_review/s9b_oracle_check_review_prompt_r1_v0.md /var/projects/toy_physics/_scratch/s9b_oracle/check_review/s9b_oracle_check_review_prompt_r1.md | cut -c1-64
03f0a92b968a6d34d8e863ccf0665f423aa97225b4c2a1122291d0aee423eff5
c111965929a3b37921b32e7a1e5551848fb823eda3b6e5d61a24753b50754ad6
```

```
$ diff /var/projects/toy_physics/_scratch/s9b_oracle/check_review/s9b_oracle_check_review_prompt_r1_v0.md /var/projects/toy_physics/_scratch/s9b_oracle/check_review/s9b_oracle_check_review_prompt_r1.md
108c108,110
<   Run Wolfram scripts as `math -script <file>`. Run one kernel at a time: before each Wolfram run, check
---
>   Run Wolfram scripts as `math -script <file>`. The guard caps each job at 32 tasks, and a Wolfram kernel that starts
>   more threads than that deadlocks with every thread idle. In any Wolfram script you write or run, set the kernel's
>   thread limits first, as the check does (its Wolfram script, lines 9–10). Run one kernel at a time: before each Wolfram run, check
```

```
$ grep -o 'Verdict: *[A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r1_claude.md
Verdict: NEEDS REVISION
Verdict: NEEDS REVISION
```

```
$ grep -o 'Verdict: *[A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r1_grok.txt | head -1
Verdict: CLEAR
```

```
$ grep -n -A2 '^## Findings' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r1_grok.txt
5:## Findings
6-
7-None.
```

## Claude 1: the simplifier and residuals printed in unreduced form
```
$ sed -n '24p;41,45p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl
domainSimplify[x_, ass_] := Refine[Cancel[Together[x]], ass];
report[tag_,a_,b_,dom_] := Module[{ans},
 progress[{tag,dom["id"]}];
 emit["RESIDUAL_OPERANDS",{a,b}];
 emit["RESIDUAL_ASSUMPTIONS",dom];
 ans=domainSimplify[a-b,dom["CAS"]];emit[tag,ans];ans];
```

```
$ sed -n '3955,3960p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
EXCHANGE_ORACLE_OPERATOR_REQUIRED_SOURCE: -((Sqrt[py`GM]*py`rho$br$profile[py`s9b$r])/(Sqrt[2]*py`s9b$r^(3/2))) + (2*Sqrt[2]*Sqrt[py`GM/py`s9b$r]*py`rho$br$profile[py`s9b$r])/py`s9b$r + Sqrt[2]*Sqrt[py`GM/py`s9b$r]*Derivative[1][py`rho$br$profile][py`s9b$r]
EXCHANGE_S9B_EVALUATED_OPERAND: (Sqrt[py`GM]*(3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r]))/(Sqrt[2]*py`s9b$r^(3/2))
WORKING_ON: {"EXCHANGE_BALANCE_OPERATOR_RESIDUAL", "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1"}
RESIDUAL_OPERANDS: {(Sqrt[py`GM]*(3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r]))/(Sqrt[2]*py`s9b$r^(3/2)), -((Sqrt[py`GM]*py`rho$br$profile[py`s9b$r])/(Sqrt[2]*py`s9b$r^(3/2))) + (2*Sqrt[2]*Sqrt[py`GM/py`s9b$r]*py`rho$br$profile[py`s9b$r])/py`s9b$r + Sqrt[2]*Sqrt[py`GM/py`s9b$r]*Derivative[1][py`rho$br$profile][py`s9b$r]}
RESIDUAL_ASSUMPTIONS: <|"id" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "CAS" -> Element[py`s9b$delta[py`s9b$r], Reals] && Element[py`s9b$V[py`s9b$r], Reals] && Element[py`s9b$xi[py`s9b$r], Reals] && py`rho$br$profile[py`s9b$r] > 0 && py`rho$br$profile[py`s9b$r] >= 0 && py`rho$br$profile[py`s9b$r] != 0 && Element[py`rho$br$profile[py`s9b$r], Reals] && Element[py`s9b$f[py`s9b$r], Reals
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: (Sqrt[2]*Sqrt[py`GM]*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] - Sqrt[2]*Sqrt[py`GM/py`s9b$r]*py`s9b$r^(3/2)*Derivative[1][py`rho$br$profile][py`s9b$r])/py`s9b$r^(3/2)
```

```
$ sed -n '3960p;3980p;7972p;7992p;11447p;11467p;13836p;13856p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: (Sqrt[2]*Sqrt[py`GM]*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] - Sqrt[2]*Sqrt[py`GM/py`s9b$r]*py`s9b$r^(3/2)*Derivative[1][py`rho$br$profile][py`s9b$r])/py`s9b$r^(3/2)
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: (-(Sqrt[2]*Sqrt[py`GM]*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r]) + Sqrt[2]*Sqrt[py`GM/py`s9b$r]*py`s9b$r^(3/2)*Derivative[1][py`rho$br$profile][py`s9b$r])/py`s9b$r^(3/2)
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: (Sqrt[2]*Sqrt[py`GM]*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] - Sqrt[2]*Sqrt[py`GM/py`s9b$r]*py`s9b$r^(3/2)*Derivative[1][py`rho$br$profile][py`s9b$r])/py`s9b$r^(3/2)
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: (-(Sqrt[2]*Sqrt[py`GM]*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r]) + Sqrt[2]*Sqrt[py`GM/py`s9b$r]*py`s9b$r^(3/2)*Derivative[1][py`rho$br$profile][py`s9b$r])/py`s9b$r^(3/2)
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: (Sqrt[2]*Sqrt[GM]*r*Derivative[1][rhoBr][r] - Sqrt[2]*Sqrt[GM/r]*r^(3/2)*Derivative[1][rhoBr][r])/r^(3/2)
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: (-(Sqrt[2]*Sqrt[GM]*r*Derivative[1][rhoBr][r]) + Sqrt[2]*Sqrt[GM/r]*r^(3/2)*Derivative[1][rhoBr][r])/r^(3/2)
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: (Sqrt[2]*Sqrt[GM]*r*Derivative[1][rhoBr][r] - Sqrt[2]*Sqrt[GM/r]*r^(3/2)*Derivative[1][rhoBr][r])/r^(3/2)
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: (-(Sqrt[2]*Sqrt[GM]*r*Derivative[1][rhoBr][r]) + Sqrt[2]*Sqrt[GM/r]*r^(3/2)*Derivative[1][rhoBr][r])/r^(3/2)
```

```
$ grep -c '^EXCHANGE_BALANCE_OPERATOR_RESIDUAL' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
18
```

```
$ grep -c '^EXCHANGE_BALANCE_OPERATOR_RESIDUAL: 0$' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
10
```

```
$ sed -n '635,639p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
FLOW_SQUARED_RESIDUAL: 0
WORKING_ON: {"FLOW_SIGNED_RESIDUAL", "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/3"}
RESIDUAL_OPERANDS: {0, Sqrt[2]*py`c$0*Sqrt[py`GM/(py`c$0^2*py`s9b$velocity$zero$radius)]}
RESIDUAL_ASSUMPTIONS: <|"id" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/3", "CAS" -> Element[py`s9b$delta[py`s9b$r], Reals] && Element[py`s9b$V[py`s9b$r], Reals] && Element[py`s9b$xi[py`s9b$r], Reals] && py`rho$br$profile[py`s9b$r] > 0 && py`rho$br$profile[py`s9b$r] >= 0 && py`rho$br$profile[py`s9b$r] != 0 && Element[py`rho$br$profile[py`s9b$r], Reals] && Element[py`s9b$f[py`s9b$r], Reals
FLOW_SIGNED_RESIDUAL: -(Sqrt[2]*Sqrt[py`GM/py`s9b$velocity$zero$radius])
```

```
$ grep -E '^SWEEP_COUNTS|^SWEEP_DISGUISED_BY_TAG' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r1_claude_evidence/scripts/residual_sweep_baseline.stdout
SWEEP_COUNTS: residuals=2806 printed_literal_zero=253 printed_nonzero_but_Simplify_zero=108 provenance_nonzero=0 timeouts=0
SWEEP_DISGUISED_BY_TAG: <|"FLOW_SIGNED_RESIDUAL" -> 100, "EXCHANGE_BALANCE_OPERATOR_RESIDUAL" -> 8|>
```

```
$ grep -E 'baseline=|control =|simplified difference=' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r1_claude_evidence/scripts/control_value_diff_visser_form.stdout | head -6
COUNTS baseline=2806 control=2806
  baseline=(Sqrt[2]*Sqrt[py`GM]*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] - Sqrt[2]*Sqrt[py`GM/py`s9b$r]*py`s9b$r^(3/2)*Derivative[1][py`rho$br$profile][py`s9b$r])/py`s9b$r^(3/2)
  control =Sqrt[2]*Sqrt[py`GM/py`s9b$r]*Derivative[1][py`rho$br$profile][py`s9b$r]
  simplified difference=-(Sqrt[2]*Sqrt[py`GM/py`s9b$r]*Derivative[1][py`rho$br$profile][py`s9b$r])
  baseline=(3*Sqrt[2]*Sqrt[py`GM]*py`rho$br$profile[py`s9b$r] + 2*Sqrt[2]*Sqrt[py`GM/py`s9b$r]*py`s9b$r^(3/2)*Derivative[1][py`rho$br$profile][py`s9b$r])/(2*py`s9b$r^(3/2))
  control =(3*Sqrt[py`GM]*py`rho$br$profile[py`s9b$r])/(Sqrt[2]*py`s9b$r^(3/2))
```

```
$ grep -B1 -A1 'EXCHANGE_BALANCE_OPERATOR_RESIDUAL' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r1_claude_evidence/ablations/kernel_a3_aligned.stdout | head -8
REPLACED_PAIRS 1351
MOVED_RESIDUAL ('EXCHANGE_BALANCE_OPERATOR_RESIDUAL', 'WL-or-no-py-symbol') 2
MOVED_RESIDUAL ('EXCHANGE_ORACLE_OPERATOR_SOURCE_RESIDUAL', 'WL-or-no-py-symbol') 2
--
PAIR
  committed: EXCHANGE_BALANCE_OPERATOR_RESIDUAL: (Sqrt[2]*Sqrt[GM]*r*Derivative[1][rhoBr][r] - Sqrt[2]*Sqrt[GM/r]*r^(3/2)*Derivative[1][rhoBr][r])/r^(3/2)
  other    : EXCHANGE_BALANCE_OPERATOR_RESIDUAL: 0
PAIR
```

```
$ sed -n '20,27p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_controls.md
CONTINUITY_ORACLE_EVALUATED_EQUATION: (2*py`rho$br$profile[py`s9b$r]*py`s9b$V[py`s9b$r])/py`s9b$r + py`s9b$V[py`s9b$r]*Derivative[1][py`rho$br$profile][py`s9b$r] + py`rho$br$profile[py`s9b$r]*Derivative[1][py`s9b$V][py`s9b$r] == 0
CONTINUITY_ORACLE_EVALUATED_EQUATION: py`rho$br$profile[py`s9b$r]*((2*py`s9b$V[py`s9b$r])/py`s9b$r + Derivative[1][py`s9b$V][py`s9b$r]) == 0
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: (Sqrt[2]*Sqrt[py`GM]*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] - Sqrt[2]*Sqrt[py`GM/py`s9b$r]*py`s9b$r^(3/2)*Derivative[1][py`rho$br$profile][py`s9b$r])/py`s9b$r^(3/2)
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: Sqrt[2]*Sqrt[py`GM/py`s9b$r]*Derivative[1][py`rho$br$profile][py`s9b$r]
```

Textually differing balance residuals: 8 of 18. Text inequality is the FORM sensitivity observation, not a claim that all unevaluated radical forms are canonical. Original source-exchange residuals and density zero sets remain printed.

```

## Grok's reading of the same residual
```
$ grep -n 'OPERATOR_RESIDUAL_TOGETHER' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r1_grok_evidence/recompute/recompute_residuals.math.stdout
16:OPERATOR_RESIDUAL_TOGETHER: 0
```

```
$ grep -n 'unreduced' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r1_grok.txt
25:`derive/12_audit_transcript.stdout`. The unreduced `ρ'` difference is an identity: an independent `Together` of the two printed writings is 0 (`recompute/recompute_residuals.math.stdout`, `OPERATOR_RESIDUAL_TOGETHER: 0`).
```

