# Grounding for S11c-b carrier ablation-harness directive — mechanical fact-lookups (E1)

> Revised 2026-09-07: the directive's "Build & review discipline" section was rewritten to end the builder's task
> at build→run-once→report (no review-launch — it had caused astra to self-review). No site claim changed; the
> fact-lookups below still ground every cited engine site.

Every site the directive cites, confirmed by grep line-location / verbatim read in the live engines
(unchanged since `dc4d4977`). Mechanical fact-lookups only (existence + line + verbatim retrieval;
no derived/algebraic predicate). Each command is shown inline below with its literal stdout — re-run any to reproduce.
Engines: `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`, `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`.

## SymPy sites
```
$ grep -n 'def build_operator\|^if __name__\|def named_tuple_row\|def casify\|def substrate_substitutions\|def filtered_substrate\|def face_generalized_force_rows' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
603:def casify(value: object) -> object:
1995:def substrate_substitutions(
2024:def filtered_substrate(
2135:def face_generalized_force_rows(
2584:def named_tuple_row(rows: sp.Tuple, name: str) -> object:
2723:def build_operator(
5585:if __name__ == "__main__":
```

```
$ grep -n 'closure_residuals\|closure_residual_sum\|mass_balance_source - closure_residual_sum' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py | head -6
2815:    closure_residuals = tuple(
2829:    closure_residual_sum = sp.Add(*closure_residuals)
2834:    mass_balance_source = sp.expand(mass_balance_source - closure_residual_sum)
2839:                sp.expand(value - closure_residual_sum)
```

```
$ grep -n 'Lambda_A_0\|Lambda_V_0' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py | head -6
408:    ("Lambda_A_0", dim_div(DIM_FLUX, DIM_AFFINITY)),
409:    ("Lambda_V_0", dim_div(DIM_FLUX, DIM_VELOCITY)),
```

```
$ grep -n 'operator\["U_BODY_BALANCE"\] =\|operator\["E_W_BALANCE"\] =\|operator\["THETA_BALANCE"\] =\|+ face_u\[a\]\|+ face_e' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py | head
2918:    operator["U_BODY_BALANCE"] = sp.Tuple(
2954:    operator["E_W_BALANCE"] = sp.Tuple(
2998:    operator["U_BODY_BALANCE"] = sp.Tuple(
3005:                        + face_u[a]
3021:                        + face_u[a]
3028:    operator["E_W_BALANCE"] = sp.Tuple(
3031:            sp.expand(named_tuple_row(reduced_e_balance, "LOCAL") + face_e),
3039:            sp.expand(named_tuple_row(reduced_e_balance, "EXPANDED") + face_e),
3042:    operator["THETA_BALANCE"] = sp.Tuple(
```

```
$ grep -n 'delta_p_plus\|delta_p_minus\|d_w_delta_p\|face_name in' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py | head
484:for face_name in ("plus", "minus"):
489:        f"d_w_delta_p_{face_name}",
```

```
$ grep -n 'e_kinetic =\|- e_kinetic' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py | head
2889:    e_kinetic = epsilon * mu_W * W_bg**2 * e_tt
2970:            sp.expand(reduced_e["EXPANDED"] - e_kinetic),
```

```
$ grep -n 'result = casify(operator)\|emit_primary("SLAB_OPERATOR"\|OPERATOR_PRIMARY_CASES\[case\] = retained_grade' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py | head
3261:    result = casify(operator), origins, mu_theta_value
4121:            OPERATOR_PRIMARY_CASES[case] = retained_grade(operator)
4142:    emit_primary("SLAB_OPERATOR", OPERATOR_PRIMARY_CASES, "slab_operator")
```

## WL sites
```
$ grep -n 'evaluatedModel\[route_String\|frozenEvaluatedModel\[route_String\|faceSources\[route_String\|Main variable-coefficient' research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl
1017:faceSources[route_String, branch_String, sign_Integer,
1230:frozenEvaluatedModel[route_String, branch_String, density_String,
1317:evaluatedModel[route_String, branch_String, density_String,
1878:(* Main variable-coefficient objects.                                    *)
```

```
$ grep -n 'flux = lambdaAResponse affinity\|tractionPressure = pressureField\|virtualWork = ' research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl | head
1080:  flux = lambdaAResponse affinity + lambdaVResponse normalVelocity;
1081:  tractionPressure = pressureField[sign] + lambdaXResponse affinity;
1082:  virtualWork = -graphMeasure /. wave -> 0;
1083:  virtualWork = virtualWork tractionPressure virtualNormalDisplacement;
2862:    flux = lambdaAResponse affinity + lambdaVResponse velocity;
```

```
$ grep -n 'pressureField\[1\] :=\|pressureField\[-1\] :=\|projectedFaceFlux\[faceAssoc\|MASS_EVOLUTION_ROW ->' research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl | head
1014:pressureField[1] := pressureUpper[xOne, xTwo, xThree, time];
1015:pressureField[-1] := pressureLower[xOne, xTwo, xThree, time];
1094:projectedFaceFlux[faceAssociation_Association] := Total[Values[
```

```
$ grep -n 'kineticEwLive =\|kineticEwLive +\|THICKNESS_ROW ->' research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl | head
1343:  kineticEwLive = muW WZero^2 D[eWField, {time, 2}];
1347:    "THICKNESS_ROW" -> (kineticEwLive + rowsLive["EW_INTERNAL"] +
```

```
$ grep -n 'processed = evaluatedModel\|extractCouplingData\[processed\]\|kernelOriginsFromOrigins\|S11CB_PRIMARIES_ONLY' research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl | head
868:If[!StringQ[Environment["S11CB_PRIMARIES_ONLY"]],
1826:kernelOriginsFromOrigins[origins_Association] := Association[
1920:If[!StringQ[Environment["S11CB_PRIMARIES_ONLY"]],
2196:  processed = evaluatedModel["EULERIAN", branch, density];
2199:  kernelData = extractCouplingData[processed];
2203:  AssociateTo[mainKernelOrigins, key -> kernelOriginsFromOrigins[
2206:  If[StringQ[Environment["S11CB_PRIMARIES_ONLY"]],
2313:If[!StringQ[Environment["S11CB_PRIMARIES_ONLY"]],
2426:If[StringQ[Environment["S11CB_PRIMARIES_ONLY"]], Quit[0]];
2557:    frozenKernelOrigins = kernelOriginsFromOrigins[frozenModel["ORIGINS"]];
```

