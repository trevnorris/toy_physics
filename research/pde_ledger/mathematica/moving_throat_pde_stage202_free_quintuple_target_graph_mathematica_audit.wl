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

scalarReduce[expr_] := Module[{res},
  res = TimeConstrained[
    Simplify[Expand[PowerExpand[expr]], Assumptions -> $Assumptions],
    20,
    $Failed
  ];
  If[res === $Failed, fail["timed reduction", expr], res]
];

normalize[expr_] := If[ListQ[expr], normalize /@ expr, scalarReduce[expr]];

allZeroQ[expr_] := If[
  ListQ[expr],
  And @@ (allZeroQ /@ expr),
  TrueQ[expr == 0]
];

pretty[expr_] := If[ListQ[expr], MatrixForm[{expr}], fmt[expr]];

expectZero[name_String, expr_] := Module[{res},
  res = normalize[expr];
  If[ListQ[res],
    Print[name, " ="];
    Print[pretty[res]],
    Print[name, " = ", pretty[res]]
  ];
  If[allZeroQ[res], pass[name], fail[name, res]];
];

expectTrue[name_String, expr_] := Module[{res},
  res = TimeConstrained[
    Simplify[expr, Assumptions -> $Assumptions],
    20,
    $Failed
  ];
  If[res === $Failed, fail["timed truth reduction", expr]];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 202 - FREE-QUINTUPLE TARGET GRAPH MATHEMATICA AUDIT"];

Clear[
  lamW, cetaU, gamma, KU, Keta, KW, muW, TU, L, sigma,
  chi0s, deltaUs, Estar, Fstar,
  Ctrtgt, Cnttgt, epsEtatgt, chiQ
];

$Assumptions =
  Element[
    {
      lamW, cetaU, gamma, KU, Keta, KW, muW, TU, L, sigma,
      chi0s, deltaUs, Estar, Fstar,
      Ctrtgt, Cnttgt, epsEtatgt, chiQ
    },
    Reals
  ] &&
  lamW > 0 && cetaU > 0 && gamma > 0 && KU > 0 &&
  Keta > 0 && KW > 0 && muW > 0 && TU > 0 &&
  L > 0 && sigma > 0 &&
  chi0s > 0 && deltaUs > 0 && Estar > 0 && Fstar > 0 &&
  Ctrtgt > 0 && Cnttgt > 0 && epsEtatgt > 0;

logSplit[lt_] := 2 Log[Pi] + lt - 2 Log[L] - Log[KU];
logTrackingBase = Log[gamma] + Log[cetaU] - Log[KU];
logDriveBase = 2 Log[gamma] + 2 Log[lamW] + Log[sigma] - Log[KU] - Log[KW];

logCtr[lt_] := (1 + deltaUs) logTrackingBase + (1 + chi0s) logSplit[lt];
logCnt[lt_, lk_, lmu_] := (
  2 Log[lamW] + lmu - lk - 2 Log[KW]
  + Estar logDriveBase
  - Fstar logSplit[lt]
);
logEta[lk_] := 2 Log[cetaU] - Log[KU] - lk;

splitU = Pi^2 TU/(L^2 KU);
trackingBase = gamma cetaU/KU;
driveBase = gamma^2 lamW^2 sigma/(KU KW);

Ctr = trackingBase^(1 + deltaUs) splitU^(1 + chi0s);
Cnt = (lamW^2 muW/(Keta KW^2)) driveBase^Estar splitU^(-Fstar);
epsEta = cetaU^2/(KU Keta);

subbanner["I. Linear log-system target graph"];

depMatrix = {
  {1 + chi0s, 0, 0},
  {0, 1, 0},
  {-Fstar, -1, 1}
};

depRhs = {
  Log[Ctrtgt]
    - (1 + deltaUs) logTrackingBase
    - (1 + chi0s) (2 Log[Pi] - 2 Log[L] - Log[KU]),
  2 Log[cetaU] - Log[KU] - Log[epsEtatgt],
  Log[Cnttgt]
    - 2 Log[lamW] + 2 Log[KW]
    - Estar logDriveBase
    + Fstar (2 Log[Pi] - 2 Log[L] - Log[KU])
};

{logTGraph, logKetaGraph, logMuGraph} = normalize[LinearSolve[depMatrix, depRhs]];

logDeltaUGraph = normalize[logSplit[logTGraph]];
deltaUGraph = Exp[logDeltaUGraph];
TGraph = Exp[logTGraph];
KetaGraph = Exp[logKetaGraph];
muGraph = Exp[logMuGraph];

Print["log(deltaU_graph) = ", fmt[logDeltaUGraph]];
Print["log(T_graph) = ", fmt[logTGraph]];
Print["log(Keta_graph) = ", fmt[logKetaGraph]];
Print["log(mu_graph) = ", fmt[logMuGraph]];

expectZero[
  "graph log-system residual",
  depMatrix . {logTGraph, logKetaGraph, logMuGraph} - depRhs
];

subbanner["II. M1-M4 direct graph reconstruction checks"];

expectZero[
  "M1 tracking monomial log residual",
  logCtr[logTGraph] - Log[Ctrtgt]
];

expectZero[
  "M2 eta dressing log residual",
  logEta[logKetaGraph] - Log[epsEtatgt]
];

expectZero[
  "M3 nontracking monomial log residual",
  logCnt[logTGraph, logKetaGraph, logMuGraph] - Log[Cnttgt]
];

muGraphFlat = scalarReduce[logMuGraph];
expectTrue[
  "M4 mu_graph is free of L and Pi",
  FreeQ[muGraphFlat, L] && FreeQ[muGraphFlat, Pi]
];
expectZero["M4 d log(mu_graph)/dL", D[logMuGraph, L]];

subbanner["III. M5 graph-error identities"];

qtr = logCtr[Log[TU]] - Log[Ctrtgt];
qnt = logCnt[Log[TU], Log[Keta], Log[muW]] - Log[Cnttgt];
qeta = logEta[Log[Keta]] - Log[epsEtatgt];

ET = Log[TU] - logTGraph;
EK = Log[Keta] - logKetaGraph;
Emu = Log[muW] - logMuGraph;

expectZero[
  "M5 graph-error quotient packet",
  {
    ET - qtr/(1 + chi0s),
    EK + qeta,
    Emu - (qnt - qeta + Fstar qtr/(1 + chi0s))
  }
];

subbanner["IV. M6 canonical same-free-quintuple projection equals graph"];

logXProjCan = {
  Log[lamW],
  Log[cetaU],
  Log[gamma],
  Log[KU],
  Log[Keta] + qeta,
  Log[KW],
  Log[muW] - qnt + qeta - Fstar qtr/(1 + chi0s),
  Log[TU] - qtr/(1 + chi0s)
};

logXGraph = {
  Log[lamW],
  Log[cetaU],
  Log[gamma],
  Log[KU],
  logKetaGraph,
  Log[KW],
  logMuGraph,
  logTGraph
};

projectionLogResiduals = MapThread[#1 - #2 &, {logXProjCan, logXGraph}];
projectionEquations = Thread[normalize[projectionLogResiduals] == ConstantArray[0, 8]];

expectZero["M6 projection component log residuals", projectionLogResiduals];
expectTrue["M6 threaded component equations", And @@ projectionEquations];

subbanner["V. M7 repair vector and reduced-family packet"];

repairFromQuotients = {
  0,
  0,
  0,
  0,
  qeta,
  0,
  -qnt + qeta - Fstar qtr/(1 + chi0s),
  -qtr/(1 + chi0s)
};

repairFromErrors = {0, 0, 0, 0, -EK, 0, -Emu, -ET};
repairResiduals = MapThread[#1 - #2 &, {repairFromQuotients, repairFromErrors}];

expectZero["M7 repair-vector rewrite residuals", repairResiduals];

reducedFamilyPacket = {chiQ - 1, ET, EK, Emu};
graphLogRules = {
  Log[TU] -> logTGraph,
  Log[Keta] -> logKetaGraph,
  Log[muW] -> logMuGraph,
  chiQ -> 1
};

expectZero[
  "M7 reduced-family packet on the target graph",
  reducedFamilyPacket /. graphLogRules
];

banner["STAGE 202 MATHEMATICA AUDIT PASSED"];
Exit[0];
