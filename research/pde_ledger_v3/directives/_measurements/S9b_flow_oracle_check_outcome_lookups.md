# Measurements: outcome of the accepted S9b flow-only oracle check (generated 2026-10-19 16:16)

Generator: `_scratch/s9b_oracle/gen/check_outcome_lookups.sh`. Transcript `scripts/out/S9b_flow_oracle_check.out`
at commit `8d91f448`. Literal retrieval and counting only; the interpretation is in the report to the user, not here.

```
$ sha256sum /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | cut -c1-64
e66ff63512e4f044bab2a60eb327665db865179531fe6077622de8fab553a25e
```

## S9b's flow-only velocities (both engines, every case the check reads)
```
$ grep '^FLOW_S9B_CASE: <|"Velocity"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sed -E 's/, "ImpliedExchange".*//' | sort | uniq -c
    416 FLOW_S9B_CASE: <|"Velocity" -> 0
    104 FLOW_S9B_CASE: <|"Velocity" -> (Sqrt[2]*Sqrt[py`GM])/Sqrt[py`s9b$r]
    104 FLOW_S9B_CASE: <|"Velocity" -> -((Sqrt[2]*Sqrt[py`GM])/Sqrt[py`s9b$r])
```

```
$ grep '^FLOW_S9B_CASE: <|"Domain"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sed -E 's/"Component".*//' | sort | uniq -c | cut -c1-200
    104 FLOW_S9B_CASE: <|"Domain" -> velocityIntegrationConstant == 0 && -(r*Inactive[Integrate][(-4*GM)/tailRadius^3, {tailRadius, r, Infinity}]) == 0, 
    208 FLOW_S9B_CASE: <|"Domain" -> velocityIntegrationConstant == 0 && -(r*Inactive[Integrate][(-4*GM)/tailRadius^3, {tailRadius, r, Infinity}]) > 0, 
```

## Comparison 1: the oracle profiles as quoted
```
$ grep -m1 '^ORACLE_QUOTED_HL2' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
ORACLE_QUOTED_HL2: <|"citation" -> "Hamilton-Lisle, same local copy, p.1, section I, equation (2)", "quoted" -> "beta = (2GM/r)^(1/2)", "conditions" -> "Spherical Schwarzschild mass M; beta in light-speed units. G M here is an oracle parameter, not S9b's orbital GM."|>
```

```
$ grep -m1 '^ORACLE_QUOTED_HL4' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
ORACLE_QUOTED_HL4: <|"citation" -> "Hamilton-Lisle, same local copy, p.3, section II.A, equation (4)", "quoted" -> "beta(r) = [2M(r)/r]^(1/2)", "conditions" -> "Spherical mass-function metric (3); units c = G = 1. M(r) is live, not assumed constant."|>
```

```
$ grep -m1 '^ORACLE_QUOTED_V7background' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
ORACLE_QUOTED_V7background: <|"citation" -> "Visser, same local copy, p.15, section 7, unnumbered background equations immediately preceding (53)", "quoted" -> "c = 1; v = sqrt(2GM/r); rho |v| proportional to 1/r^2; rho proportional to r^(-3/2); c^2 = dp/drho; p = p_infinity + c^2 rho; p-p_infinity proportional to c^2 r^(-3/2)", "conditions" -> "c is positio
```

```
$ grep -m1 '^ORACLE_QUOTED_V54' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
ORACLE_QUOTED_V54: <|"citation" -> "Visser, same local copy, p.16, section 8, equation (54)", "quoted" -> "v = c r_0^2/r^2", "conditions" -> "Incompressible, spherically symmetric fluid. rho is position independent; the barotropic assumption makes pressure and c position independent. Continuity gives v proportional to 1/r^2. For the bubble realization use th
```

## Comparison 1: zero tests of the squared-profile residuals, per domain (distinct values and counts)
```
$ grep -F '"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sed -E 's/.*"value" -> (.*)\|>$/\1/' | sort | uniq -c
      4 2*py`GM*py`s9b$r^3 == clockRate^2*radius0O^4*soundO^2
      4 2*py`GM*py`s9b$r^3 == radius0O^4*soundO^2
      4 py`c$0^2*radius0O^4 == 2*py`GM*py`s9b$r^3
      3 py`GM == clockRate^2*gmH
      4 py`GM == clockRate^2*gmV
      3 py`GM == clockRate^2*massH[py`s9b$r]
      3 py`GM == gmH
      4 py`GM == gmV
      3 py`GM == massH[py`s9b$r]
      3 py`GM == py`c$0^2*gmH
      4 py`GM == py`c$0^2*gmV
      3 py`GM == py`c$0^2*massH[py`s9b$r]
     10 True
```

```
$ grep -F '"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sed -E 's/.*"value" -> (.*)\|>$/\1/' | sort | uniq -c
      4 2*GM*r^3 == c0^2*radius0O^4
      4 2*GM*r^3 == clockRate^2*radius0O^4*soundO^2
      4 2*GM*r^3 == radius0O^4*soundO^2
      3 GM == c0^2*gmH
      4 GM == c0^2*gmV
      3 GM == c0^2*massH[r]
      3 GM == clockRate^2*gmH
      4 GM == clockRate^2*gmV
      3 GM == clockRate^2*massH[r]
      3 GM == gmH
      4 GM == gmV
      3 GM == massH[r]
     10 True
```

```
$ grep -F '"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/1"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sed -E 's/.*"value" -> (.*)\|>$/\1/' | sort | uniq -c
      4 2*py`GM*py`s9b$r^3 == clockRate^2*radius0O^4*soundO^2
      4 2*py`GM*py`s9b$r^3 == radius0O^4*soundO^2
      4 py`c$0^2*radius0O^4 == 2*py`GM*py`s9b$r^3
      3 py`GM == clockRate^2*gmH
      4 py`GM == clockRate^2*gmV
      3 py`GM == clockRate^2*massH[py`s9b$r]
      3 py`GM == gmH
      4 py`GM == gmV
      3 py`GM == massH[py`s9b$r]
      3 py`GM == py`c$0^2*gmH
      4 py`GM == py`c$0^2*gmV
      3 py`GM == py`c$0^2*massH[py`s9b$r]
     10 True
```

```
$ grep -F '"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/general"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sed -E 's/.*"value" -> (.*)\|>$/\1/' | sort | uniq -c
      3 2*py`GM + py`s9b$r^2*py`s9b$velocity$integration$constant == 2*clockRate^2*gmH
      4 2*py`GM + py`s9b$r^2*py`s9b$velocity$integration$constant == 2*clockRate^2*gmV
      3 2*py`GM + py`s9b$r^2*py`s9b$velocity$integration$constant == 2*clockRate^2*massH[py`s9b$r]
      3 2*py`GM + py`s9b$r^2*py`s9b$velocity$integration$constant == 2*gmH
      4 2*py`GM + py`s9b$r^2*py`s9b$velocity$integration$constant == 2*gmV
      3 2*py`GM + py`s9b$r^2*py`s9b$velocity$integration$constant == 2*massH[py`s9b$r]
      3 2*py`GM + py`s9b$r^2*py`s9b$velocity$integration$constant == 2*py`c$0^2*gmH
      4 2*py`GM + py`s9b$r^2*py`s9b$velocity$integration$constant == 2*py`c$0^2*gmV
      3 2*py`GM + py`s9b$r^2*py`s9b$velocity$integration$constant == 2*py`c$0^2*massH[py`s9b$r]
      4 2*py`GM*py`s9b$r^3 + py`s9b$r^5*py`s9b$velocity$integration$constant == clockRate^2*radius0O^4*soundO^2
      4 2*py`GM*py`s9b$r^3 + py`s9b$r^5*py`s9b$velocity$integration$constant == radius0O^4*soundO^2
      4 py`c$0^2*radius0O^4 == 2*py`GM*py`s9b$r^3 + py`s9b$r^5*py`s9b$velocity$integration$constant
     10 py`s9b$velocity$integration$constant == 0
```

```
$ grep -F '"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/general"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sed -E 's/.*"value" -> (.*)\|>$/\1/' | sort | uniq -c
      3 2*GM + r^2*velocityIntegrationConstant == 2*c0^2*gmH
      4 2*GM + r^2*velocityIntegrationConstant == 2*c0^2*gmV
      3 2*GM + r^2*velocityIntegrationConstant == 2*c0^2*massH[r]
      3 2*GM + r^2*velocityIntegrationConstant == 2*clockRate^2*gmH
      4 2*GM + r^2*velocityIntegrationConstant == 2*clockRate^2*gmV
      3 2*GM + r^2*velocityIntegrationConstant == 2*clockRate^2*massH[r]
      3 2*GM + r^2*velocityIntegrationConstant == 2*gmH
      4 2*GM + r^2*velocityIntegrationConstant == 2*gmV
      3 2*GM + r^2*velocityIntegrationConstant == 2*massH[r]
      4 2*GM*r^3 + r^5*velocityIntegrationConstant == c0^2*radius0O^4
      4 2*GM*r^3 + r^5*velocityIntegrationConstant == clockRate^2*radius0O^4*soundO^2
      4 2*GM*r^3 + r^5*velocityIntegrationConstant == radius0O^4*soundO^2
     10 velocityIntegrationConstant == 0
```

## Comparison 1: PY deflection case 1 against Hamilton-Lisle (2), operands, residual and zero test
```
$ sed -n '99,108p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -E '^(RESIDUAL_OPERANDS|FLOW_SQUARED_RESIDUAL|FLOW_SIGNED_RESIDUAL|RESIDUAL_ZERO_TEST)'
FLOW_SQUARED_RESIDUAL: (2*(py`GM - gmH))/py`s9b$r
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> py`GM == gmH|>
RESIDUAL_OPERANDS: {-(Sqrt[2]*Sqrt[py`GM/py`s9b$r]), Sqrt[2]*Sqrt[gmH/py`s9b$r]}
FLOW_SIGNED_RESIDUAL: -(Sqrt[2]*Sqrt[py`GM/py`s9b$r]) - (Sqrt[2]*Sqrt[gmH])/Sqrt[py`s9b$r]
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SIGNED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> Sqrt[py`GM] + Sqrt[gmH] == 0|>
```

```
$ sed -n '167,168p;188,197p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -E '^(FLOW_ORACLE_OPERAND|RESIDUAL_OPERANDS|FLOW_SIGNED_RESIDUAL|RESIDUAL_ZERO_TEST)'
FLOW_ORACLE_OPERAND: <|"orientation" -> <|"label" -> "coordinate_drift", "equations" -> {"HL1", "HL2"}|>, "value" -> -(Sqrt[2]*Sqrt[gmH/py`s9b$r])|>
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> py`GM == gmH|>
RESIDUAL_OPERANDS: {-(Sqrt[2]*Sqrt[py`GM/py`s9b$r]), -(Sqrt[2]*Sqrt[gmH/py`s9b$r])}
FLOW_SIGNED_RESIDUAL: -(Sqrt[2]*Sqrt[py`GM/py`s9b$r]) + (Sqrt[2]*Sqrt[gmH])/Sqrt[py`s9b$r]
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SIGNED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> py`GM == gmH|>
```

## Comparison 1: PY deflection case 1 against Visser section 8, first block
```
$ sed -n '3670,3681p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -E '^(RESIDUAL_OPERANDS|RESIDUAL_ZERO_TEST)'
RESIDUAL_OPERANDS: {(2*py`GM)/py`s9b$r, (radius0O^4*soundO^2)/py`s9b$r^4}
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> 2*py`GM*py`s9b$r^3 == radius0O^4*soundO^2|>
RESIDUAL_OPERANDS: {-(Sqrt[2]*Sqrt[py`GM/py`s9b$r]), (radius0O^2*soundO)/py`s9b$r^2}
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SIGNED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> Sqrt[2]*Sqrt[py`GM*py`s9b$r^3] + radius0O^2*soundO == 0|>
```

## Comparison 2: continuity and exchange
```
$ grep -m1 '^ORACLE_QUOTED_V5' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
ORACLE_QUOTED_V5: <|"citation" -> "Visser, same local copy, p.4, section 2.1, equation (5)", "quoted" -> "partial_t rho + nabla dot (rho v) = 0", "conditions" -> "Background continuity; steady exterior specialization of the stated equation. No off-sheet source appears."|>
```

```
$ grep -m1 '^CONTINUITY_ORACLE_EVALUATED_EQUATION' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
CONTINUITY_ORACLE_EVALUATED_EQUATION: (2*py`rho$br$profile[py`s9b$r]*py`s9b$V[py`s9b$r])/py`s9b$r + py`s9b$V[py`s9b$r]*Derivative[1][py`rho$br$profile][py`s9b$r] + py`rho$br$profile[py`s9b$r]*Derivative[1][py`s9b$V][py`s9b$r] == 0
```

```
$ grep -m1 '^SOURCE_EXCHANGE_SIGN' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
SOURCE_EXCHANGE_SIGN: -1
```

```
$ sed -n '4740,4760p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -E '^(EXCHANGE_S9B_EVALUATED_OPERAND|EXCHANGE_BALANCE_OPERATOR_RESIDUAL|EXCHANGE_RESIDUAL|RESIDUAL_ZERO_TEST|DENSITY_ODE_SOLUTIONS_ON_OPEN_COMPONENT|HL_CONTINUITY_COMPARISON)'
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: (Sqrt[2]*Sqrt[py`GM]*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] - Sqrt[2]*Sqrt[py`GM/py`s9b$r]*py`s9b$r^(3/2)*Derivative[1][py`rho$br$profile][py`s9b$r])/py`s9b$r^(3/2)
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "EXCHANGE_BALANCE_OPERATOR_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> True|>
EXCHANGE_RESIDUAL: (Sqrt[py`GM]*(3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r]))/(Sqrt[2]*py`s9b$r^(3/2))
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "EXCHANGE_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> 3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] == 0|>
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "EXCHANGE_ORACLE_OPERATOR_SOURCE_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> 3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] == 0|>
DENSITY_ODE_SOLUTIONS_ON_OPEN_COMPONENT: {{py`rho$br$profile[py`s9b$r] -> C[1]/py`s9b$r^(3/2)}}
HL_CONTINUITY_COMPARISON: <|"oracle" -> "UNMATCHED: Hamilton-Lisle states no material continuity/exchange operand.", "S9b" -> (Sqrt[py`GM]*(3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r]))/(Sqrt[2]*py`s9b$r^(3/2)), "residual" -> "UNMATCHED"|>
EXCHANGE_S9B_EVALUATED_OPERAND: -((Sqrt[py`GM]*(3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r]))/(Sqrt[2]*py`s9b$r^(3/2)))
```

```
$ grep '^RESIDUAL_ZERO_TEST: <|"residual_tag" -> "EXCHANGE_BALANCE_OPERATOR_RESIDUAL"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -o '"value" -> [A-Za-z]*' | sort | uniq -c
     18 "value" -> True
```

```
$ grep '^RESIDUAL_ZERO_TEST: <|"residual_tag" -> "EXCHANGE_RESIDUAL"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sed -E 's/.*"domain" -> "([^"]*)", "value" -> (.*)\|>$/\1 | \2/'
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1 | 3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] == 0
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/2 | 3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] == 0
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/3 | True
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/4 | True
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/5 | True
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/6 | True
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/1 | 3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] == 0
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/2 | 3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] == 0
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/3 | True
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/4 | True
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/5 | True
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/6 | True
WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1 | 3*rhoBr[r] + 2*r*Derivative[1][rhoBr][r] == 0
WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/2 | 3*rhoBr[r] + 2*r*Derivative[1][rhoBr][r] == 0
WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/3 | True
WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/1 | 3*rhoBr[r] + 2*r*Derivative[1][rhoBr][r] == 0
WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/2 | 3*rhoBr[r] + 2*r*Derivative[1][rhoBr][r] == 0
WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/3 | True
```

```
$ grep '^DENSITY_ODE_SOLUTIONS_ON_OPEN_COMPONENT' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sort | uniq -c
      4 DENSITY_ODE_SOLUTIONS_ON_OPEN_COMPONENT: {{py`rho$br$profile[py`s9b$r] -> C[1]/py`s9b$r^(3/2)}}
      4 DENSITY_ODE_SOLUTIONS_ON_OPEN_COMPONENT: {{rhoBr[r] -> C[1]/r^(3/2)}}
```

```
$ grep -m1 '^WL_FLOW_SPEED_RELATIONS_FROM_ATTACHED_DOMAIN' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
WL_FLOW_SPEED_RELATIONS_FROM_ATTACHED_DOMAIN: {c0 == Sqrt[muPerp[r]/rhoBr[r]]}
```

## Comparison 3: optical characteristics
```
$ grep '^OPTICAL_INPUT_STATUS' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sort | uniq -c
      2 OPTICAL_INPUT_STATUS: "SUPPLIED leading-eikonal input; not an S9b result"
```

```
$ grep '^WORKING_ON: {"OPTICAL_CHARACTERISTIC_RESIDUAL"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sort | uniq -c
    192 WORKING_ON: {"OPTICAL_CHARACTERISTIC_RESIDUAL", "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/general"}
     64 WORKING_ON: {"OPTICAL_CHARACTERISTIC_RESIDUAL", "WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/general"}
```

```
$ grep -c '^OPTICAL_CHARACTERISTIC_RESIDUAL: 0$' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
3
```

```
$ grep -n '^OPTICAL_CHARACTERISTIC_RESIDUAL: 0$' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
10866:OPTICAL_CHARACTERISTIC_RESIDUAL: 0
10890:OPTICAL_CHARACTERISTIC_RESIDUAL: 0
17075:OPTICAL_CHARACTERISTIC_RESIDUAL: 0
```

```
$ sed -n '10862,10867p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -E '^(OPTICAL_IDENTIFICATION_MAP|RESIDUAL_OPERANDS|OPTICAL_CHARACTERISTIC_RESIDUAL|RESIDUAL_ZERO_TEST)'
RESIDUAL_OPERANDS: {-(py`c$0^2*(otherQ^2 + radialQ^2 + tangentQ^2)) + (frequencyQ - radialQ*py`s9b$V[py`s9b$r])^2, Sqrt[soundO^(-2)]*(frequencyQ^2/Sqrt[soundO^(-2)] - (py`c$0^2*otherQ^2)/Sqrt[soundO^(-2)] - (py`c$0^2*radialQ^2)/Sqrt[soundO^(-2)] - (py`c$0^2*tangentQ^2)/Sqrt[soundO^(-2)] - (2*frequencyQ*radialQ*py`s9b$V[py`s9b$r])/Sqrt[soundO^(-2)] + (radialQ
OPTICAL_CHARACTERISTIC_RESIDUAL: 0
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "OPTICAL_CHARACTERISTIC_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/general", "value" -> True|>
```

```
$ sed -n '9810,9820p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -E '^(OPTICAL_IDENTIFICATION_MAP|RESIDUAL_OPERANDS|RESIDUAL_ZERO_TEST)'
RESIDUAL_OPERANDS: {(frequencyQ - radialQ*py`s9b$U)^2 - py`c$0^2*(otherQ^2 + radialQ^2 + tangentQ^2), (2*py`c$0^2*gmH*radialQ^2 + frequencyQ^2*py`s9b$r - py`c$0^2*otherQ^2*py`s9b$r - py`c$0^2*radialQ^2*py`s9b$r + 2*Sqrt[2]*py`c$0*frequencyQ*radialQ*Sqrt[gmH/py`s9b$r]*py`s9b$r - py`c$0^2*py`s9b$r*tangentQ^2)/py`s9b$r}
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "OPTICAL_CHARACTERISTIC_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/general", "value" -> radialQ*(2*py`c$0^2*gmH*radialQ + 2*Sqrt[2]*py`c$0*frequencyQ*Sqrt[gmH*py`s9b$r] + py`s9b$r*py`s9b$U*(2*frequencyQ - radialQ*py`s9b$U)) == 0|>
OPTICAL_IDENTIFICATION_MAP: <|"oracle" -> "HL2", "clock" -> {"conditional_speed_match", py`c$0}, "binding_operands" -> {<|"rule" -> py`s9b$chi -> py`c$0^2, "object" -> "PY_LOCAL_S9B_SPEED_IDENTIFICATION / equality LHS after flow restriction (conditional chi binding)"|>, <|"rule" -> py`s9b$A -> 1, "object" -> "Hamilton-Lisle (1) / coordinate spatial rr coeffi
RESIDUAL_OPERANDS: {-(py`c$0^2*(otherQ^2 + radialQ^2 + tangentQ^2)) + (frequencyQ - radialQ*py`s9b$V[py`s9b$r])^2, (2*py`c$0^2*gmH*radialQ^2 + frequencyQ^2*py`s9b$r - py`c$0^2*otherQ^2*py`s9b$r - py`c$0^2*radialQ^2*py`s9b$r + 2*Sqrt[2]*py`c$0*frequencyQ*radialQ*Sqrt[gmH/py`s9b$r]*py`s9b$r - py`c$0^2*py`s9b$r*tangentQ^2)/py`s9b$r}
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "OPTICAL_CHARACTERISTIC_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/general", "value" -> radialQ*(2*py`c$0*(py`c$0*gmH*radialQ + Sqrt[2]*frequencyQ*Sqrt[gmH*py`s9b$r]) + py`s9b$r*py`s9b$V[py`s9b$r]*(2*frequencyQ - radialQ*py`s9b$V[py`s9b$r])) == 0|>
OPTICAL_IDENTIFICATION_MAP: <|"oracle" -> "HL2", "clock" -> {"conditional_speed_match", py`c$0}, "binding_operands" -> {<|"rule" -> py`s9b$chi -> py`c$0^2, "object" -> "PY_LOCAL_S9B_SPEED_IDENTIFICATION / equality LHS after flow restriction (conditional chi binding)"|>, <|"rule" -> py`s9b$A -> 1, "object" -> "Unit radial coefficient imposed by normalization 
```

## The identification map and order statement
```
$ grep -m1 '^SUPPLIED_ORDER' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
SUPPLIED_ORDER: <|"S9b" -> "Exterior flow-only restriction; Part B retained first optical order: velocity through its squared contribution, with supplied half-order velocity counting. General integration constant remains live in a separately printed, counting-excluded family.", "Optics" -> "Supplied leading-eikonal dispersion, not an S9b result; coefficients
```

```
$ grep -o '"[a-zA-Z0-9_]*" -> "UNMATCHED' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sort | uniq -c
    224 "frame" -> "UNMATCHED
    432 "mass" -> "UNMATCHED
     18 "oracle" -> "UNMATCHED
    208 "orientation_identification" -> "UNMATCHED
     32 "residual" -> "UNMATCHED
     14 "S9b" -> "UNMATCHED
     18 "status" -> "UNMATCHED
```

