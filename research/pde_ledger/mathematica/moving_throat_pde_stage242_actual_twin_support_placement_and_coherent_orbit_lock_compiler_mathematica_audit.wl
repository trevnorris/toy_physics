ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

normalizeExpr[expr_] := Module[{res},
  res = If[
    MatrixQ[expr],
    Map[FullSimplify[#, Assumptions -> $Assumptions] &, expr, {2}],
    If[
      VectorQ[expr],
      Map[FullSimplify[#, Assumptions -> $Assumptions] &, expr],
      FullSimplify[expr, Assumptions -> $Assumptions]
    ]
  ];
  res
];

allZeroQ[expr_] := If[
  ListQ[expr],
  And @@ Flatten[Map[TrueQ[# == 0] &, expr, {ArrayDepth[expr]}]],
  TrueQ[expr == 0]
];

prettyArray[arr_] := If[VectorQ[arr], MatrixForm[{arr}], MatrixForm[arr]];

expectZero[name_String, expr_] := Module[{res},
  res = normalizeExpr[expr];
  If[ListQ[res],
    Print[name, " ="];
    Print[prettyArray[res]];
    If[allZeroQ[res], pass[name], fail[name, res]],
    Print[name, " = ", fmt[res]];
    If[allZeroQ[res], pass[name], fail[name, res]]
  ];
];

banner["STAGE 242 — ACTUAL TWIN-SUPPORT PLACEMENT AND COHERENT ORBIT-LOCK COMPILER"];

expectTrue[name_String, statement_] := Module[{res},
  res = FullSimplify[statement, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

Clear[
  chi0, deltaU, ZW, epsW, epsEta, Lambda, Lambda0, zeta, beta,
  chi0Ref, deltaURef, ZWRef, epsWRef, epsEtaRef, LambdaRef,
  Bstar, Cstar, dchi0, ddeltaU, dZW, depsW, depsEta, dLambda,
  lambdaWin, epsilonWin, thetaVar, xiVar, rVar
];

$Assumptions =
  Element[
    {
      chi0, deltaU, ZW, epsW, epsEta, Lambda, Lambda0, zeta, beta,
      chi0Ref, deltaURef, ZWRef, epsWRef, epsEtaRef, LambdaRef,
      Bstar, Cstar, dchi0, ddeltaU, dZW, depsW, depsEta, dLambda,
      thetaVar, xiVar, rVar
    },
    Reals
  ] &&
  chi0 > 0 && deltaU > 0 && ZW > 0 && 0 < epsW < 1 &&
  0 < epsEta < 1 && Lambda > 0 && Lambda0 > 0 && beta > 0 &&
  chi0Ref > 0 && deltaURef > 0 && ZWRef > 0 && 0 < epsWRef < 1 &&
  0 < epsEtaRef < 1 && LambdaRef > 0 && Bstar != 0 && Cstar != 0;

logDrift[expr_, variables_List, driftSymbols_List] := Module[{terms},
  terms = MapThread[#1 D[Log[expr], #1] #2 &, {variables, driftSymbols}];
  FullSimplify[Together[Plus @@ terms], Assumptions -> $Assumptions]
];

subbanner["M1-M2. Placement coordinate and sigma consistency"];

epsilonSelected = FullSimplify[
  epsW (1 - (2/11) deltaU/(1 + deltaU)),
  Assumptions -> $Assumptions
];
mixedCapacity = FullSimplify[
  8 Lambda (1 - epsilonSelected)/Pi^2,
  Assumptions -> $Assumptions
];
traceLoad = FullSimplify[(4/3) mixedCapacity, Assumptions -> $Assumptions];
rhoSelected = FullSimplify[
  Pi^2 traceLoad/(16 Lambda),
  Assumptions -> $Assumptions
];
sigmaFromRho = FullSimplify[
  4/(3 rhoSelected) - 2,
  Assumptions -> $Assumptions
];
sigmaFromEpsilon = FullSimplify[
  2 epsilonSelected/(1 - epsilonSelected),
  Assumptions -> $Assumptions
];

expectZero[
  "M1 selected coordinate rho",
  rhoSelected - (2/3) (1 - epsilonSelected)
];
expectZero[
  "M2 sigma from placement epsilon",
  sigmaFromRho - sigmaFromEpsilon
];
expectZero[
  "M2 independent sigma-rho consistency",
  sigmaFromEpsilon - (4/(3 rhoSelected) - 2)
];

subbanner["M3. Threshold rewrites and strict selected window"];

rhoWallLambda = FullSimplify[
  2 (1 + beta^2)/(3 (2 + beta^2)),
  Assumptions -> $Assumptions
];
rhoUnitLambda = FullSimplify[
  2 (1 + beta^2)/(3 (1 + beta + beta^2)),
  Assumptions -> $Assumptions
];
epsilonFromRho[r_] := FullSimplify[1 - (3/2) r, Assumptions -> $Assumptions];

expectZero[
  "M3 epsilon_WLambda threshold rewrite",
  epsilonFromRho[rhoWallLambda] - 1/(2 + beta^2)
];
expectZero[
  "M3 epsilon_ULambda threshold rewrite",
  epsilonFromRho[rhoUnitLambda] - beta/(1 + beta + beta^2)
];

demandRatio = FullSimplify[traceLoad/mixedCapacity, Assumptions -> $Assumptions];
expectZero["M3 selected demand ratio", demandRatio - 4/3];
windowCertificate = Resolve[
  ForAll[
    {lambdaWin, epsilonWin},
    Implies[
      lambdaWin > 0 && 0 < epsilonWin < 1,
      1 < ((4/3) (8 lambdaWin (1 - epsilonWin)/Pi^2))/
        (8 lambdaWin (1 - epsilonWin)/Pi^2) < 2
    ]
  ],
  Reals
];
expectTrue["M3 strict twin-window inclusion by Resolve", windowCertificate];

subbanner["M4-M5. Reduced bridge and support-blind closed packet"];

rescaledWallCharge = FullSimplify[ZW Lambda0/Lambda, Assumptions -> $Assumptions];
trailObservable = FullSimplify[
  (1 + chi0/(1 + deltaU))/(1 + chi0),
  Assumptions -> $Assumptions
];
targetObservable = FullSimplify[
  Lambda (1 - epsEta) (1 - epsilonSelected)^2/(ZW (1 + chi0)^2),
  Assumptions -> $Assumptions
];
targetWithRescaledCharge = FullSimplify[
  Lambda0 (1 - epsEta) (1 - epsilonSelected)^2/
    (rescaledWallCharge (1 + chi0)^2),
  Assumptions -> $Assumptions
];

expectZero["M4 reduced-state bridge for target observable", targetObservable - targetWithRescaledCharge];
expectZero["M5 zeta derivative of epsilon", D[epsilonSelected, zeta]];
expectZero["M5 zeta derivative of throat ratio", D[trailObservable, zeta]];
expectZero["M5 zeta derivative of target ratio", D[targetObservable, zeta]];

epsilonRef = FullSimplify[
  epsWRef (1 - (2/11) deltaURef/(1 + deltaURef)),
  Assumptions -> $Assumptions
];
trailRef = FullSimplify[
  (1 + chi0Ref/(1 + deltaURef))/(1 + chi0Ref),
  Assumptions -> $Assumptions
];
targetRef = FullSimplify[
  LambdaRef (1 - epsEtaRef) (1 - epsilonRef)^2/
    (ZWRef (1 + chi0Ref)^2),
  Assumptions -> $Assumptions
];
closedObservablePacket = {
  -Cstar Log[trailObservable/trailRef],
  Bstar Log[trailObservable/trailRef] +
    Log[(1 - epsEta)/(1 - epsEtaRef)] -
    Log[targetObservable/targetRef],
  Log[epsEta/epsEtaRef]
};
expectZero[
  "M5 closed-form q-packet zeta derivatives",
  D[closedObservablePacket, zeta]
];

subbanner["M6. Total logarithmic differentials and orbit compiler"];

epsilonLogDrift = logDrift[
  epsilonSelected,
  {epsW, deltaU},
  {depsW, ddeltaU}
];
epsilonLogFormula = FullSimplify[
  depsW - (2 deltaU/((1 + deltaU) (11 + 9 deltaU))) ddeltaU,
  Assumptions -> $Assumptions
];
expectZero["M6 dln epsilon by total log differential", epsilonLogDrift - epsilonLogFormula];

trailLogDrift = logDrift[
  trailObservable,
  {chi0, deltaU},
  {dchi0, ddeltaU}
];
trailLogFormula = FullSimplify[
  -(
    chi0 deltaU/((1 + chi0) (1 + deltaU) (1 + chi0 + deltaU))
  ) ((1 + deltaU) dchi0 + (1 + chi0) ddeltaU),
  Assumptions -> $Assumptions
];
expectZero["M6 dln R_tr by total log differential", trailLogDrift - trailLogFormula];

targetLogDrift = logDrift[
  targetObservable,
  {Lambda, ZW, epsEta, chi0, epsW, deltaU},
  {dLambda, dZW, depsEta, dchi0, depsW, ddeltaU}
];
targetLogFormula = FullSimplify[
  dLambda - dZW - 2 chi0/(1 + chi0) dchi0 -
    epsEta/(1 - epsEta) depsEta -
    2 epsilonSelected/(1 - epsilonSelected) epsilonLogDrift,
  Assumptions -> $Assumptions
];
expectZero["M6 dln R_target by total log differential", targetLogDrift - targetLogFormula];

trailSigma = (1 + deltaU) dchi0 + (1 + chi0) ddeltaU;
trailCoefficient = FullSimplify[
  chi0 deltaU/((1 + chi0) (1 + deltaU) (1 + chi0 + deltaU)),
  Assumptions -> $Assumptions
];
thetaFromPacket = FullSimplify[-trailCoefficient trailSigma, Assumptions -> $Assumptions];
cEta = FullSimplify[epsEta/(1 - epsEta), Assumptions -> $Assumptions];
xiFromDefinition = FullSimplify[-targetLogDrift - cEta depsEta, Assumptions -> $Assumptions];
xiFormula = FullSimplify[
  -dLambda + dZW + 2 chi0/(1 + chi0) dchi0 +
    2 epsilonSelected/(1 - epsilonSelected) epsilonLogDrift,
  Assumptions -> $Assumptions
];
rFromDefinition = FullSimplify[-xiFromDefinition - cEta depsEta, Assumptions -> $Assumptions];

expectZero["M6 Theta_1 packet form equals dln R_tr", thetaFromPacket - trailLogDrift];
expectZero["M6 Xi_1 definition", xiFromDefinition - xiFormula];
expectZero["M6 R_1 definition", rFromDefinition - targetLogDrift];

directLogPacket = {trailLogDrift, targetLogDrift, depsEta};
orbitCompiler = {
  {1, 0, 0},
  {0, -1, -cEta},
  {0, 1, 0}
};
orbitLogPacket = {thetaFromPacket, xiFromDefinition, rFromDefinition};
expectZero[
  "M6 orbit compiler maps direct packet",
  orbitLogPacket - orbitCompiler . directLogPacket
];
expectZero["M6 orbit compiler determinant", Det[orbitCompiler] - cEta];
expectZero[
  "M6 LinearSolve recovers direct packet",
  LinearSolve[orbitCompiler, orbitLogPacket] - directLogPacket
];
formalRecovered = FullSimplify[
  LinearSolve[orbitCompiler, {thetaVar, xiVar, rVar}],
  Assumptions -> $Assumptions
];
expectZero[
  "M6 formal inverse map",
  formalRecovered - {
    thetaVar,
    rVar,
    -((1 - epsEta) (xiVar + rVar))/epsEta
  }
];

subbanner["M7. Two-packet split"];

mixedPacket = FullSimplify[
  8 ZW (1 + chi0)^2/(Pi^2 (1 - epsEta) (1 - epsilonSelected)),
  Assumptions -> $Assumptions
];
supportFactor = FullSimplify[
  1 + zeta (1 - epsilonSelected)/(1 - zeta epsilonSelected),
  Assumptions -> $Assumptions
];
throatPacket = FullSimplify[mixedPacket supportFactor, Assumptions -> $Assumptions];

expectZero["M7 target times mixed packet equals C_mix", targetObservable mixedPacket - mixedCapacity];
expectZero[
  "M7 support-packet zeta derivative",
  D[throatPacket, zeta] -
    mixedPacket (1 - epsilonSelected)/(1 - zeta epsilonSelected)^2
];

subbanner["Probe-only rational sample point"];

sampleRules = {
  chi0 -> 3/2,
  deltaU -> 2/3,
  ZW -> 13/17,
  epsW -> 1/3,
  epsEta -> 1/5,
  Lambda -> 7/11,
  Lambda0 -> 5/7,
  zeta -> 1
};

Print["epsilon = ", fmt[FullSimplify[epsilonSelected /. sampleRules]]];
Print["rhoSelected = ", fmt[FullSimplify[rhoSelected /. sampleRules]]];
Print["sigmaFromRho = ", fmt[FullSimplify[sigmaFromRho /. sampleRules]]];
Print["trailObservable = ", fmt[FullSimplify[trailObservable /. sampleRules]]];
Print["targetObservable = ", fmt[FullSimplify[targetObservable /. sampleRules]]];
Print["mixedPacket = ", fmt[FullSimplify[mixedPacket /. sampleRules]]];
Print["throatPacket = ", fmt[FullSimplify[throatPacket /. sampleRules]]];

Print[""];
Print["All Stage 242 symbolic checks passed."];
