(* pathA_29 v3 brane<->bulk return gate, Mathematica engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
Off[Limit::alimv];

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_29_brane_bulk_return.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
sympyJson = FileNameJoin[{scratchDir, "pathA_29_brane_bulk_return_sympy.json"}];
jsonOut = FileNameJoin[{scratchDir, "pathA_29_brane_bulk_return_mathematica.json"}];

If[! FileExistsQ[sympyJson], fail["missing SymPy JSON for engine agreement: " <> sympyJson]];

$Assumptions = omega > 0 && cS > 0 && d > 0 && a > 0 && r > 0 &&
  rhoB > 0 && epsilon0 > 0 && epsilon1 > 0 && M0 > 0 &&
  kWarp > 0 && W > 0 && Element[{D1, GammaUniform, w, m}, Reals] &&
  m >= 0;

nonzeroQ[expr_] := ! TrueQ[FullSimplify[expr == 0]];

omegaOrder[expr_, max_] := Module[{series, n, coeff},
  series = Expand[expr];
  For[n = 0, n <= max, n++,
    coeff = FullSimplify[Coefficient[series, omega, n]];
    If[nonzeroQ[coeff], Return[n]]
  ];
  fail["no nonzero omega coefficient through " <> ToString[max] <> " in " <> ToString[series, InputForm]]
];

radialExponent[flow_] := Module[{p},
  p = FullSimplify[-Limit[r D[Log[flow], r], r -> Infinity]];
  If[! NumericQ[p], fail["could not extract radial exponent from " <> ToString[flow, InputForm]]];
  p
];

radialResidual[expr_, kappaSquared_] := FullSimplify[D[expr, {r, 2}] + 2 D[expr, r]/r + kappaSquared expr];

classifyGate[plist_] := Module[{target, eq},
  target = 2;
  eq = TrueQ[FullSimplify[# == target]] & /@ plist;
  Which[
    And @@ eq, "RETURN_RESIDUAL_PREDICTION",
    ! Or @@ eq, "RETURN_NOGO",
    True, "BC_DEPENDENT"
  ]
];

z = FullSimplify[omega d/cS];
outgoingBasis = Exp[I omega w/cS];
returningBasis = Exp[-I omega w/cS];
forwardPhase = FullSimplify[(outgoingBasis /. w -> d)/(outgoingBasis /. w -> 0)];
returnPhase = FullSimplify[(returningBasis /. w -> 0)/(returningBasis /. w -> d)];
transportPhase = FullSimplify[forwardPhase returnPhase];
tau = FullSimplify[D[Log[transportPhase], omega]/I];
If[! TrueQ[FullSimplify[Exp[I omega tau] == transportPhase]],
  fail["round-trip transit phase was not solved from the Helmholtz basis"]
];
alpha0 = FullSimplify[1/(1 + epsilon0)];
alpha1 = FullSimplify[1/(1 + epsilon1)];
neighborFraction0 = FullSimplify[epsilon0/(1 + epsilon0)];
neighborFraction1 = FullSimplify[epsilon1/(1 + epsilon1)];

T0Full = FullSimplify[alpha0 transportPhase];
T1Full = FullSimplify[alpha1 transportPhase];
T0Series = Expand[Normal@Series[T0Full, {omega, 0, 4}]];
T1Series = Expand[Normal@Series[T1Full, {omega, 0, 2}]];
T0DC = FullSimplify[Limit[T0Full, omega -> 0]];
T1DC = FullSimplify[Limit[T1Full, omega -> 0]];
If[! TrueQ[FullSimplify[T0DC == alpha0 && T1DC == alpha1]],
  fail["DC transfer values did not follow from continuity fractions"]
];

sigma0 = omegaOrder[T0Series - T0DC, 4];
sigma1 = omegaOrder[T1Series - T1DC, 2];
nu0 = omegaOrder[1 - T0Series, 4];
nu1 = omegaOrder[1 - T1Series, 2];
If[nu0 =!= 0 || nu1 =!= 0, fail["finite DC sink did not leave omega^0 deviation-from-one"]];

rawKernel0 = I a omega/cS;
rawKernel1 = I a^3 omega^3/(2 cS^3);
rawKernel2 = I a^5 omega^5/(27 cS^5);
R0 = FullSimplify[-M0 T0Series];
R1 = FullSimplify[-D1 T1Series];
A0Res = FullSimplify[rawKernel0 M0 (1 - T0Series)];
A1Res = FullSimplify[rawKernel1 D1 (1 - T1Series)];
pRes0 = 1 + nu0;
pRes1 = 3 + nu1;

JReturn0 = FullSimplify[alpha0 M0];
JNeighbor0 = FullSimplify[neighborFraction0 M0];
steadyBalance0 = FullSimplify[M0 - JReturn0 - JNeighbor0];
JReturn1 = FullSimplify[alpha1 D1];
JNeighbor1 = FullSimplify[neighborFraction1 D1];
steadyBalance1 = FullSimplify[D1 - JReturn1 - JNeighbor1];
If[steadyBalance0 =!= 0 || steadyBalance1 =!= 0, fail["steady circulation balance failed"]];

ZThroat = -M0;
ZReturn = FullSimplify[M0 T0DC];
ZReplenishmentLocalized = 0;
ZBoundaryDof = 0;
Z = FullSimplify[ZThroat + ZReturn + ZReplenishmentLocalized + ZBoundaryDof];
ZLocalFormula = FullSimplify[-M0 (1 - T0DC)];
If[! TrueQ[FullSimplify[Z == ZLocalFormula]], fail["local-channel Z reduction failed"]];
ZNegativeCertificate = FullSimplify[-Z (1 + epsilon0)/(M0 epsilon0)];
If[ZNegativeCertificate =!= 1, fail["Z sign certificate failed"]];
strictT0Limit = FullSimplify[Limit[T0Full, epsilon0 -> 0, Direction -> "FromAbove"]];
strictT0Series = Expand[Normal@Series[strictT0Limit, {omega, 0, 4}]];
strictNu0 = omegaOrder[1 - strictT0Series, 4];
strictPRes0 = 1 + strictNu0;
strictPRes1 = 3 + strictNu0;
strictZLocalLimit = FullSimplify[Limit[Z, epsilon0 -> 0, Direction -> "FromAbove"]];
If[strictZLocalLimit =!= 0, fail["perfect-return local Z limit failed"]];

f0Abs = FullSimplify[1/Sqrt[d]];
f1Abs = FullSimplify[Sqrt[2/d] Cos[Pi w/d]];
m0Abs = 0;
m1Abs = FullSimplify[Pi/d];
norm0Abs = FullSimplify[Integrate[f0Abs^2, {w, 0, d}]];
norm1Abs = FullSimplify[Integrate[f1Abs^2, {w, 0, d}]];
absResiduals = {
  FullSimplify[D[f0Abs, {w, 2}] + m0Abs^2 f0Abs],
  FullSimplify[D[f0Abs, w] /. w -> 0],
  FullSimplify[D[f0Abs, w] /. w -> d],
  FullSimplify[D[f1Abs, {w, 2}] + m1Abs^2 f1Abs],
  FullSimplify[D[f1Abs, w] /. w -> 0],
  FullSimplify[D[f1Abs, w] /. w -> d]
};
If[norm0Abs =!= 1 || norm1Abs =!= 1 || Or @@ (# =!= 0 & /@ absResiduals),
  fail["destructuring/absorbing transverse solve failed"]
];

f0Bloch = FullSimplify[1/Sqrt[d]];
f1Bloch = FullSimplify[Sqrt[2/d] Cos[2 Pi w/d]];
m0Bloch = 0;
m1Bloch = FullSimplify[2 Pi/d];
norm0Bloch = FullSimplify[Integrate[f0Bloch^2, {w, 0, d}]];
norm1Bloch = FullSimplify[Integrate[f1Bloch^2, {w, 0, d}]];
blochResiduals = {
  FullSimplify[D[f0Bloch, {w, 2}] + m0Bloch^2 f0Bloch],
  FullSimplify[(f0Bloch /. w -> d) - (f0Bloch /. w -> 0)],
  FullSimplify[(D[f0Bloch, w] /. w -> d) - (D[f0Bloch, w] /. w -> 0)],
  FullSimplify[D[f1Bloch, {w, 2}] + m1Bloch^2 f1Bloch],
  FullSimplify[(f1Bloch /. w -> d) - (f1Bloch /. w -> 0)],
  FullSimplify[(D[f1Bloch, w] /. w -> d) - (D[f1Bloch, w] /. w -> 0)]
};
If[norm0Bloch =!= 1 || norm1Bloch =!= 1 || Or @@ (# =!= 0 & /@ blochResiduals),
  fail["Bloch transverse solve failed"]
];

dynamicGeneral = FullSimplify[
  DSolveValue[gDynamic''[r] + 2 gDynamic'[r]/r + (omega/cS)^2 gDynamic[r] == 0, gDynamic[r], r]
];
greenZeroDynamic = FullSimplify[dynamicGeneral /. {C[1] -> 0, C[2] -> I omega/(2 Pi d cS)}];
If[radialResidual[greenZeroDynamic, (omega/cS)^2] =!= 0, fail["dynamic 3D radial solve residual failed"]];
greenDynamicLimit = FullSimplify[Limit[greenZeroDynamic, omega -> 0]];
flowDynamicLimit = FullSimplify[-D[greenDynamicLimit, r]];
pDynamic = radialExponent[flowDynamicLimit];

staticAbsGeneral = FullSimplify[
  DSolveValue[gAbsStatic''[r] + 2 gAbsStatic'[r]/r == 0, gAbsStatic[r], r]
];
greenZeroStatic = FullSimplify[staticAbsGeneral /. {C[1] -> -1/(4 Pi d), C[2] -> 0}];
If[radialResidual[greenZeroStatic, 0] =!= 0, fail["static destructuring 3D radial solve residual failed"]];
flowZeroStatic = FullSimplify[-D[greenZeroStatic, r]];
pOmegaZeroSolve = radialExponent[flowZeroStatic];
If[! TrueQ[FullSimplify[pDynamic == pOmegaZeroSolve]], fail["static-dynamic exponents disagree"]];

staticBlochGeneral = FullSimplify[
  DSolveValue[gBlochStatic''[r] + 2 gBlochStatic'[r]/r == 0, gBlochStatic[r], r]
];
greenBlochZeroStatic = FullSimplify[staticBlochGeneral /. {C[1] -> -1/(4 Pi d), C[2] -> 0}];
If[radialResidual[greenBlochZeroStatic, 0] =!= 0, fail["static Bloch 3D radial solve residual failed"]];
flowBlochZeroStatic = FullSimplify[-D[greenBlochZeroStatic, r]];

massiveAbsGeneral = FullSimplify[
  DSolveValue[gAbsMassive''[r] + 2 gAbsMassive'[r]/r - m1Abs^2 gAbsMassive[r] == 0, gAbsMassive[r], r]
];
massiveAbsTail = FullSimplify[massiveAbsGeneral /. {C[1] -> 1/(4 Pi), C[2] -> 0}];
If[radialResidual[massiveAbsTail, -m1Abs^2] =!= 0, fail["destructuring Yukawa radial solve residual failed"]];

massiveBlochGeneral = FullSimplify[
  DSolveValue[gBlochMassive''[r] + 2 gBlochMassive'[r]/r - m1Bloch^2 gBlochMassive[r] == 0, gBlochMassive[r], r]
];
massiveBlochTail = FullSimplify[massiveBlochGeneral /. {C[1] -> 1/(4 Pi), C[2] -> 0}];
If[radialResidual[massiveBlochTail, -m1Bloch^2] =!= 0, fail["Bloch Yukawa radial solve residual failed"]];

candidateResidual = radialResidual[greenZeroStatic, 0];
perturbedResidual = radialResidual[greenZeroStatic/r^4, 0];
If[candidateResidual =!= 0, fail["candidate radial solution failed the 3D radial operator"]];
If[! nonzeroQ[perturbedResidual], fail["counterfactual wrong falloff passed the 3D radial operator"]];

pAbs = radialExponent[flowZeroStatic];
pBloch = radialExponent[flowBlochZeroStatic];
pEq2Indicator = Boole[TrueQ[FullSimplify[pAbs == 2]]];
If[pEq2Indicator =!= 1 || ! TrueQ[FullSimplify[pAbs == pBloch]],
  fail["localizing DC-sink p derivation failed"]
];

greenAbsStatic = FullSimplify[Z greenZeroStatic + CAbs1 massiveAbsTail];
greenBlochStatic = FullSimplify[Z greenBlochZeroStatic + CBloch1 massiveBlochTail];

antiMeasure = Exp[2 kWarp w];
antiZeroNormCutoff = FullSimplify[Integrate[antiMeasure, {w, 0, W}]];
antiZeroNormLimit = FullSimplify[Limit[antiZeroNormCutoff, W -> Infinity]];
If[antiZeroNormLimit =!= Infinity, fail["anti-localizing zero mode unexpectedly normalizable"]];
continuumGreenStatic = FullSimplify[Integrate[Exp[-m r], {m, 0, Infinity}]/(4 Pi r)];
continuumFlow = FullSimplify[-D[continuumGreenStatic, r]];
pDelocalizing = radialExponent[continuumFlow];
delocalizingVerdict = classifyGate[{pDelocalizing}];
If[delocalizingVerdict =!= "RETURN_NOGO", fail["no-go reachable control failed"]];

headline = classifyGate[{pAbs, pBloch}];
If[headline =!= "RETURN_RESIDUAL_PREDICTION", fail["unexpected v3 headline: " <> headline]];

(* Dimensional homogeneity with dimensions ordered {M,L,T}. *)
dimA = {0, 1, 0};
dimD = {0, 1, 0};
dimR = {0, 1, 0};
dimCS = {0, 1, -1};
dimOmega = {0, 0, -1};
dimM0 = {0, 0, -1};
dimD1 = {0, 1, -1};
dimRho = {0, -3, 0};
dimMass = {1, 0, 0};
dimK = {1, 14, -2};
If[dimOmega + dimD - dimCS =!= {0, 0, 0}, fail["omega*d/c_s dimension check failed"]];
If[dimA + dimOmega - dimCS + dimM0 =!= dimM0, fail["A0 dimension check failed"]];
If[3 dimA + 3 dimOmega - 3 dimCS + dimD1 =!= dimD1, fail["A1 dimension check failed"]];
If[dimM0 - dimRho - 2 dimR =!= {0, 1, -1}, fail["radial flow dimension check failed"]];
If[dimK + 4 dimRho - dimMass =!= {0, 2, -2}, fail["c_s^2 dimension check failed"]];
KEos = 1/500;
rhoFrozenMEq1 = FullSimplify[(1/(5 KEos))^(1/4)];
csSquaredFrozen = FullSimplify[5 KEos rhoFrozenMEq1^4];
If[rhoFrozenMEq1 =!= Sqrt[10] || csSquaredFrozen =!= 1, fail["frozen c_s^2 slice check failed"]];

sympyResults = Import[sympyJson, "RawJSON"];
sympyExprs = sympyResults["engine_agreement"]["mathematica_exprs"];
sympyDigest = sympyResults["engine_agreement"]["sympy_expression_digest"];

assertEngine[name_, actual_] := Module[{expectedText, expectedExpr},
  expectedText = sympyExprs[name];
  If[! StringQ[expectedText], fail["missing SymPy expression for " <> name]];
  expectedExpr = ToExpression[expectedText, InputForm];
  If[! TrueQ[FullSimplify[actual == expectedExpr]],
    fail["engine disagreement " <> name <> ": Mathematica got " <> ToString[actual, InputForm] <>
      ", SymPy exported " <> expectedText]
  ];
];

engineKeys = {
  "T0_at_DC",
  "T1_at_DC",
  "p",
  "p_eq_2",
  "destructuring_p",
  "bloch_p",
  "dynamic_limit_exponent",
  "static_solve_exponent",
  "zero_mode.r_dependence",
  "green_function",
  "destructuring_green_function",
  "bloch_green_function",
  "delocalizing_p",
  "A0_residual",
  "A1_residual",
  "Z"
};

actuals = <|
  "T0_at_DC" -> T0DC,
  "T1_at_DC" -> T1DC,
  "p" -> pAbs,
  "p_eq_2" -> pEq2Indicator,
  "destructuring_p" -> pAbs,
  "bloch_p" -> pBloch,
  "dynamic_limit_exponent" -> pDynamic,
  "static_solve_exponent" -> pOmegaZeroSolve,
  "zero_mode.r_dependence" -> greenZeroStatic,
  "green_function" -> greenZeroDynamic,
  "destructuring_green_function" -> greenAbsStatic,
  "bloch_green_function" -> greenBlochStatic,
  "delocalizing_p" -> pDelocalizing,
  "A0_residual" -> A0Res,
  "A1_residual" -> A1Res,
  "Z" -> Z
|>;

Scan[assertEngine[#, actuals[#]] &, engineKeys];

required = {"p", "p_eq_2", "dynamic_limit_exponent", "static_solve_exponent", "zero_mode.r_dependence", "green_function"};
If[! SubsetQ[engineKeys, required], fail["engine checked quantities missing directive-required entries"]];

out = <|
  "schema" -> "pathA_29_brane_bulk_return_mathematica/v3",
  "status" -> "PASS",
  "headline" -> headline,
  "checked_quantities" -> engineKeys,
  "sympy_expression_digest" -> sympyDigest,
  "computed" -> <|
    "p" -> ToString[pAbs, InputForm],
    "p_eq_2" -> TrueQ[pEq2Indicator == 1],
    "dynamic_limit_exponent" -> ToString[pDynamic, InputForm],
    "static_solve_exponent" -> ToString[pOmegaZeroSolve, InputForm],
    "delocalizing_p" -> ToString[pDelocalizing, InputForm],
    "no_go_reachable" -> delocalizingVerdict
  |>
|>;

Export[jsonOut, out, "RawJSON"];
Print["PASS pathA_29_brane_bulk_return_mathematica"];
Print[ExportString[<|"json" -> jsonOut, "headline" -> headline|>, "RawJSON"]];
