# Measurements behind S11b_comparator_build_directive.md (rule 2) — FORMAT facts, not values
# (This twin is the orchestrator's record; it is NOT handed to the builder, which stays blind.)

## Engines emit some objects in DIFFERENT structures (delta #3) — PY positional Tuple vs WL Association
```
$ grep -h "^PY_S11B_ZPERM_SLICE_MAP:" <pyout> | head -c 140 ; grep -h "^WL_..." <wlout>
PY_S11B_ZPERM_SLICE_MAP: Tuple(Mul(Integer(-1), Symbol('Lambda_A_0', real=True), Pow(Symbol('rho_m', positive=True), Integer(-1))), Equality
WL_S11B_ZPERM_SLICE_MAP: <|"RAW_PRESSURE_COEFFICIENT" -> lambdaA0/(rhoM*(1 - I*omega*tauCommon)), "MAPPING" -> {{lambdaPressureCoefficient -
```

## No dim(W_0) tag exists; DIM_ tags are in different shells (delta #5)
```
(no DIM_W_0 tag in either engine)
PY_S11B_DIM_THICKNESS_RESPONSE: ImmutableDenseMatrix([[Integer(2)], [Integer(2)], [Integer(-1)]])
WL_S11B_DIM_THICKNESS_RESPONSE: <|"VECTOR_L_T_M" -> {3, 2, -1}, "ROUTE" -> HoldForm[thicknessDisplacementDimension - gen
```

## parse_mathematica RAISES on <|...|> (delta #4)
```
$ python3 -c "from sympy.parsing.mathematica import parse_mathematica as p; p(chr(60)+chr(124)+...)"
RAISES: SyntaxError
```

## Native booleans present in both streams (delta #2) — counts (structural presence, not a physics value)
```
PY True/False lines: 108 ; WL: 26
S10 residual scores False==0 as agreement (the bug delta #2 fixes):
False==0 -> True ; True==1 -> True
```
