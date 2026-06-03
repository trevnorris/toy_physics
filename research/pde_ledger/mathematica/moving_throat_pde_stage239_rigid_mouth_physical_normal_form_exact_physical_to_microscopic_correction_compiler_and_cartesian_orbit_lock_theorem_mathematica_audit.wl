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

expectTrue[name_String, test_] := Module[{res},
  res = FullSimplify[test, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 239 — RIGID-MOUTH PHYSICAL NORMAL FORM"];

Clear[
  U, V, T2, T2ref, epsEta, epsEtaRef, Rtarget, RtargetRef, Lambda0,
  chi0, deltaU, ZW, OmegaW2, eps, Rtr, RtrRef, zeta, Mmix,
  h, dlnT2, dlnEps, DeltaT, DeltaKeta, DeltaMu,
  T2sbFn, epsEtaSbFn
];

$Assumptions = (
  Element[{U, V, zeta, Mmix, h, dlnT2, dlnEps}, Reals] &&
  Element[
    {
      T2, T2ref, epsEta, epsEtaRef, Rtarget, RtargetRef, Lambda0,
      chi0, deltaU, ZW, OmegaW2, eps, Rtr, RtrRef,
      DeltaT, DeltaKeta, DeltaMu
    },
    Reals
  ] &&
  T2 > 0 && T2ref > 0 &&
  0 < epsEta < 1 && 0 < epsEtaRef < 1 &&
  Rtarget > 0 && RtargetRef > 0 && Lambda0 > 0 &&
  chi0 > 0 && deltaU > 0 && ZW > 0 && OmegaW2 > 0 &&
  eps > 0 && eps != 1 && Rtr > 0 && RtrRef > 0
);

physicalCoordinates = {U, V};
uLog = Log[T2/T2ref];
vLog = Log[epsEta/epsEtaRef];
Mphys = IdentityMatrix[2];

subbanner["I. Chart, target ratio, and physical projectors"];

expectZero["M1 transfer chart exponent", Exp[uLog] - T2/T2ref];
expectZero["M1 eta chart exponent", Exp[vLog] - epsEta/epsEtaRef];
expectZero["M1 diagonal packet map", Mphys . physicalCoordinates - {U, V}];

targetFromBranch = Lambda0 (1 - epsEta)/T2;
targetRefFromBranch = Lambda0 (1 - epsEtaRef)/T2ref;
ratioFromPremise = FullSimplify[
  targetFromBranch/targetRefFromBranch,
  Assumptions -> $Assumptions
];
ratioUV = ((1 - epsEtaRef Exp[V])/(1 - epsEtaRef)) Exp[-U];
ratioPremiseInChart = FullSimplify[
  ratioFromPremise /. {T2 -> T2ref Exp[U], epsEta -> epsEtaRef Exp[V]},
  Assumptions -> $Assumptions
];

expectZero["M2 selected branch product identity", targetFromBranch T2 - Lambda0 (1 - epsEta)];
expectZero["M2 target quotient in UV chart", ratioPremiseInChart - ratioUV];

projectorT = DiagonalMatrix[{1, 0}];
projectorEta = DiagonalMatrix[{0, 1}];
ratioAlong[coords_] := FullSimplify[
  ratioUV /. Thread[physicalCoordinates -> coords],
  Assumptions -> $Assumptions
];
ratioTransfer = ratioAlong[projectorT . physicalCoordinates];
ratioDressing = ratioAlong[projectorEta . physicalCoordinates];

expectZero["M3 transfer projector squares", projectorT . projectorT - projectorT];
expectZero["M3 eta projector squares", projectorEta . projectorEta - projectorEta];
expectZero["M3 projector cross term TE", projectorT . projectorEta];
expectZero["M3 projector cross term ET", projectorEta . projectorT];
expectZero["M3 projectors complete chart", projectorT + projectorEta - IdentityMatrix[2]];
expectZero["M3 finite-leg product", ratioUV - ratioTransfer ratioDressing];

subbanner["II. Dependent correction and native left inverse"];

dependentDelta = {
  0,
  -physicalCoordinates[[2]],
  physicalCoordinates[[1]] - physicalCoordinates[[2]]
};
compilerJacobian = Table[
  D[dependentDelta[[row]], physicalCoordinates[[col]]],
  {row, Length[dependentDelta]},
  {col, Length[physicalCoordinates]}
];

expectZero["M4 Jacobian rebuilds dependent vector", compilerJacobian . physicalCoordinates - dependentDelta];
expectZero["M4 compiler passes through physical identity", compilerJacobian . Mphys - compilerJacobian];

deltaSymbols = {DeltaT, DeltaKeta, DeltaMu};
leftNative = PseudoInverse[compilerJacobian];

expectZero["M5 native left inverse matrix", leftNative . compilerJacobian - IdentityMatrix[2]];
expectZero["M5 native inverse coordinate formulas", leftNative . deltaSymbols - {DeltaMu - DeltaKeta, -DeltaKeta}];
expectZero["M5 dependent vector returns UV", leftNative . dependentDelta - physicalCoordinates];

subbanner["III. Stage 238 branch carried into the chart"];

T2Stage238 = ZW (1 + chi0)^2/(OmegaW2 (1 - eps)^2);
RtargetStage238 = Lambda0 (1 - epsEta)/T2Stage238;
qNtRigidStage238 = FullSimplify[
  PowerExpand[
    Log[(1 - epsEta)/(1 - epsEtaRef)] -
    Log[RtargetStage238/targetRefFromBranch]
  ],
  Assumptions -> $Assumptions
];
UStage238 = Log[T2Stage238/T2ref];
VStage238 = Log[epsEta/epsEtaRef];
xPhysStage238 = {UStage238, VStage238};
yDepStage238 = compilerJacobian . {qNtRigidStage238, VStage238};

expectZero["M1 carried q_nt exponential", Exp[qNtRigidStage238] - T2Stage238/T2ref];
expectZero["M1 carried q_nt equals U", PowerExpand[qNtRigidStage238 - UStage238]];
expectZero["M1 carried q_eta equals V", PowerExpand[VStage238 - vLog]];
expectZero["M2 carried target identity", RtargetStage238 T2Stage238 - Lambda0 (1 - epsEta)];
expectZero["M4 carried dependent map", PowerExpand[yDepStage238 - compilerJacobian . xPhysStage238]];

subbanner["IV. Axis images and correction packets"];

transferAxisImage = compilerJacobian . (projectorT . physicalCoordinates);
etaAxisImage = compilerJacobian . (projectorEta . physicalCoordinates);

expectZero["M6 transfer axis image", transferAxisImage - {0, 0, U}];
expectZero["M6 eta axis image", etaAxisImage + V {0, 1, 1}];
expectZero["M6 axis sum reconstructs dependent vector", dependentDelta - transferAxisImage - etaAxisImage];

deltaStatic = -transferAxisImage;
deltaEtaRest = -etaAxisImage;
deltaOrbit = -dependentDelta;

expectZero["M7 static correction packet", deltaStatic - {0, 0, -U}];
expectZero["M7 eta-rest correction packet", deltaEtaRest - {0, V, V}];
expectZero["M7 orbit correction packet", deltaOrbit - {0, V, V - U}];
expectZero["M7 correction packet split", deltaOrbit - deltaStatic - deltaEtaRest];
expectZero["M7 orbit correction cancels defect", dependentDelta + deltaOrbit];

subbanner["V. Support-blind propagation"];

supportVars = {zeta, Mmix};
T2Support = T2sbFn[zeta, Mmix];
epsSupport = epsEtaSbFn[zeta, Mmix];
uSupport = Log[T2Support/T2ref];
vSupport = Log[epsSupport/epsEtaRef];
supportRules = {
  Derivative[1, 0][T2sbFn][zeta, Mmix] -> 0,
  Derivative[0, 1][T2sbFn][zeta, Mmix] -> 0,
  Derivative[1, 0][epsEtaSbFn][zeta, Mmix] -> 0,
  Derivative[0, 1][epsEtaSbFn][zeta, Mmix] -> 0
};
supportResiduals[vec_] := FullSimplify[
  Flatten[Outer[D, vec, supportVars]] /. supportRules,
  Assumptions -> $Assumptions
];

dependentSupport = compilerJacobian . {uSupport, vSupport};
staticSupport = -(compilerJacobian . {uSupport, 0});
orbitSupport = -dependentSupport;

expectZero["M8 Stage238 T2 zeta derivative", D[T2Stage238, zeta]];
expectZero["M8 Stage238 T2 Mmix derivative", D[T2Stage238, Mmix]];
expectZero["M8 chain rule for U support variables", supportResiduals[{uSupport}]];
expectZero["M8 chain rule for V support variables", supportResiduals[{vSupport}]];
expectZero["M8 dependent correction support gradient", supportResiduals[dependentSupport]];
expectZero["M8 static correction support gradient", supportResiduals[staticSupport]];
expectZero["M8 orbit correction support gradient", supportResiduals[orbitSupport]];

subbanner["VI. Cartesian orbit lock and infinitesimal packet"];

orbitLockReduced = Reduce[
  dependentDelta[[2]] == 0 && dependentDelta[[3]] == 0,
  physicalCoordinates,
  Reals
];
Print["M9 Reduce orbit-lock set = ", fmt[orbitLockReduced]];
expectTrue["M9 lock set is chart origin", Equivalent[orbitLockReduced, U == 0 && V == 0]];

expectZero["M9 origin removes dependent vector", dependentDelta /. {U -> 0, V -> 0}];
expectZero["M9 transfer equality fixes U", uLog /. T2 -> T2ref];
expectZero["M9 eta equality fixes V", vLog /. epsEta -> epsEtaRef];
expectZero["M9 target quotient at origin", (ratioUV /. {U -> 0, V -> 0}) - 1];

perturbedDependent = dependentDelta /. Thread[
  physicalCoordinates -> {
    Log[(T2 Exp[h dlnT2])/T2],
    Log[(epsEta Exp[h dlnEps])/epsEta]
  }
];
firstOrderDependent = FullSimplify[
  PowerExpand[D[perturbedDependent, h] /. h -> 0],
  Assumptions -> $Assumptions
];

expectZero["M9 first-order dependent vector", firstOrderDependent - {0, -dlnEps, dlnT2 - dlnEps}];

Print[""];
Print["All Stage 239 symbolic checks passed."];
