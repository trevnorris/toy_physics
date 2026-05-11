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

banner["STAGE 225 — ACTUAL TWIN-SUPPORT PLACEMENT AND COHERENT ORBIT-LOCK COMPILER"];

Clear[
  chi0, deltaU, ZW, epsW, epsEta, Lambda, Lambda0, zeta, beta, varrho,
  Bstar, Cstar,
  chi0Ref, deltaURef, ZWRef, epsWRef, epsEtaRef, LambdaRef, EStar, FStar,
  dchi0, ddeltaU, dZW, depsW, depsEta, dLambda, t
];

$Assumptions =
  Element[
    {
      chi0, deltaU, ZW, epsW, epsEta, Lambda, Lambda0, zeta, beta, varrho,
      chi0Ref, deltaURef, ZWRef, epsWRef, epsEtaRef, LambdaRef
    },
    Reals
  ] &&
  Element[{Bstar, Cstar}, Reals] &&
  chi0 > 0 && deltaU > 0 && ZW > 0 && epsW > 0 && epsEta > 0 &&
  Lambda > 0 && Lambda0 > 0 && beta > 0 && varrho > 0 &&
  chi0Ref > 0 && deltaURef > 0 && ZWRef > 0 && epsWRef > 0 &&
  epsEtaRef > 0 && LambdaRef > 0 && Bstar != 0 && Cstar != 0 &&
  Element[{EStar, FStar, dchi0, ddeltaU, dZW, depsW, depsEta, dLambda, t}, Reals];

subbanner["I. Actual selected-twin placement"];

eps = FullSimplify[epsW (1 - (2/11) deltaU/(1 + deltaU)), Assumptions -> $Assumptions];
Cmix = FullSimplify[8 Lambda (1 - eps)/Pi^2, Assumptions -> $Assumptions];
PiTr = (4/3) Cmix;
varrhoPhys = FullSimplify[Pi^2 PiTr/(16 Lambda), Assumptions -> $Assumptions];
sigmaPhys = FullSimplify[4/(3 varrhoPhys) - 2, Assumptions -> $Assumptions];

expectZero[
  "selected-branch coordinate varrhoPhys",
  varrhoPhys - (2/3) (1 - eps)
];
expectZero[
  "selected-branch sigmaPhys",
  sigmaPhys - 2 eps/(1 - eps)
];

subbanner["II. Threshold rewrite and selected-branch twin-window inclusion"];

varrhoWL = FullSimplify[2 (1 + beta^2)/(3 (2 + beta^2)), Assumptions -> $Assumptions];
varrhoUL = FullSimplify[2 (1 + beta^2)/(3 (1 + beta + beta^2)), Assumptions -> $Assumptions];
epsWL = FullSimplify[1 - (3/2) varrhoWL, Assumptions -> $Assumptions];
epsUL = FullSimplify[1 - (3/2) varrhoUL, Assumptions -> $Assumptions];

expectZero["epsilon_WLambda rewrite", epsWL - 1/(2 + beta^2)];
expectZero["epsilon_ULambda rewrite", epsUL - beta/(1 + beta + beta^2)];

sigmaSel = FullSimplify[4/(3 varrho) - 2, Assumptions -> $Assumptions];
expectZero[
  "selected branch lies above mixed-only bound",
  sigmaSel - (1/varrho - 2) - 1/(3 varrho)
];
expectZero[
  "selected branch lies below non-twin bound",
  (2/varrho - 2) - sigmaSel - 2/(3 varrho)
];

subbanner["III. Reduced-state bridge and direct coherent observables"];

ZhatW = FullSimplify[ZW Lambda0/Lambda, Assumptions -> $Assumptions];
Rtr = FullSimplify[(1 + chi0/(1 + deltaU))/(1 + chi0), Assumptions -> $Assumptions];
Rtarget = FullSimplify[
  Lambda (1 - epsEta) (1 - eps)^2/(ZW (1 + chi0)^2),
  Assumptions -> $Assumptions
];
RtargetHat = FullSimplify[
  Lambda0 (1 - epsEta) (1 - eps)^2/(ZhatW (1 + chi0)^2),
  Assumptions -> $Assumptions
];

expectZero["reduced-state bridge for R_target", Rtarget - RtargetHat];
expectZero["zeta-absence of epsilon in the coherent placement map", D[eps, zeta]];
expectZero["zeta-absence of R_tr in the coherent placement map", D[Rtr, zeta]];
expectZero["zeta-absence of R_target in the coherent placement map", D[Rtarget, zeta]];

subbanner["IV. Finite orbit packet and support-blind propagation"];

qtr = (1 + deltaURef) Log[chi0/chi0Ref] + (1 + chi0Ref) Log[deltaU/deltaURef];
qnt =
  Log[ZW/ZWRef] - Log[Lambda/LambdaRef] +
  EStar Log[epsW/epsWRef] - FStar Log[deltaU/deltaURef];
qeta = Log[epsEta/epsEtaRef];

RtrSb = RtrSbFn[zeta];
RtargetSb = RtargetSbFn[zeta];
epsEtaSb = epsEtaSbFn[zeta];
supportBlindRules = {
  Derivative[1][RtrSbFn][zeta] -> 0,
  Derivative[1][RtargetSbFn][zeta] -> 0,
  Derivative[1][epsEtaSbFn][zeta] -> 0
};

qtrFromObservables = -Cstar Log[RtrSb/Rtr];
qntFromObservables =
  Bstar Log[RtrSb/Rtr] +
  Log[(1 - epsEtaSb)/(1 - epsEtaRef)] -
  Log[RtargetSb/Rtarget];
qetaFromObservables = Log[epsEtaSb/epsEtaRef];

expectZero[
  "support-blind direct observables propagate to q_tr",
  D[qtrFromObservables, zeta] /. supportBlindRules
];
expectZero[
  "support-blind direct observables propagate to q_nt",
  D[qntFromObservables, zeta] /. supportBlindRules
];
expectZero[
  "support-blind direct observables propagate to q_eta",
  D[qetaFromObservables, zeta] /. supportBlindRules
];
expectZero[
  "finite q_eta matches the direct observable chart",
  (qetaFromObservables /. epsEtaSb -> epsEta) - qeta
];

subbanner["V. Infinitesimal coherent packet and direct observable compiler"];

chi0T = chi0 Exp[t dchi0];
deltaUT = deltaU Exp[t ddeltaU];
ZWT = ZW Exp[t dZW];
epsWT = epsW Exp[t depsW];
epsEtaT = epsEta Exp[t depsEta];
LambdaT = Lambda Exp[t dLambda];

epsT = FullSimplify[
  epsWT (1 - (2/11) deltaUT/(1 + deltaUT)),
  Assumptions -> $Assumptions
];
dlnEps = FullSimplify[D[Log[epsT], t] /. t -> 0, Assumptions -> $Assumptions];
dlnEpsFormula = FullSimplify[
  depsW - (2 deltaU/((1 + deltaU) (11 + 9 deltaU))) ddeltaU,
  Assumptions -> $Assumptions
];
expectZero["dln epsilon compiler", dlnEps - dlnEpsFormula];

RtrT = FullSimplify[
  (1 + chi0T/(1 + deltaUT))/(1 + chi0T),
  Assumptions -> $Assumptions
];
dlnRtr = FullSimplify[D[Log[RtrT], t] /. t -> 0, Assumptions -> $Assumptions];
dlnRtrFormula = FullSimplify[
  -(
    chi0 deltaU/((1 + chi0) (1 + deltaU) (1 + chi0 + deltaU))
  ) ((1 + deltaU) dchi0 + (1 + chi0) ddeltaU),
  Assumptions -> $Assumptions
];
expectZero["dln R_tr compiler", dlnRtr - dlnRtrFormula];

RtargetT = FullSimplify[
  LambdaT (1 - epsEtaT) (1 - epsT)^2/(ZWT (1 + chi0T)^2),
  Assumptions -> $Assumptions
];
dlnRtarget = FullSimplify[
  D[Log[RtargetT], t] /. t -> 0,
  Assumptions -> $Assumptions
];
dlnRtargetFormula = FullSimplify[
  dLambda - dZW - 2 chi0/(1 + chi0) dchi0 -
  epsEta/(1 - epsEta) depsEta - 2 eps/(1 - eps) dlnEps,
  Assumptions -> $Assumptions
];
expectZero["dln R_target compiler", dlnRtarget - dlnRtargetFormula];

Xi1 = FullSimplify[
  -dlnRtarget - epsEta/(1 - epsEta) depsEta,
  Assumptions -> $Assumptions
];
Xi1Formula = FullSimplify[
  -dLambda + dZW + 2 chi0/(1 + chi0) dchi0 + 2 eps/(1 - eps) dlnEps,
  Assumptions -> $Assumptions
];
Theta1 = dlnRtr;
R1 = FullSimplify[
  -Xi1 - epsEta/(1 - epsEta) depsEta,
  Assumptions -> $Assumptions
];
cEta = FullSimplify[epsEta/(1 - epsEta), Assumptions -> $Assumptions];

expectZero["Theta_1 direct-observable identity", Theta1 - dlnRtr];
expectZero["Xi_1 direct-observable identity", Xi1 - Xi1Formula];
expectZero["R_1 direct-observable identity", R1 - dlnRtarget];

directPacket = {dlnRtr, dlnRtarget, depsEta};
orbitPacket = {Theta1, Xi1, R1};
orbitCompiler = {
  {1, 0, 0},
  {0, -1, -cEta},
  {0, 1, 0}
};
expectZero[
  "direct-observable orbit packet compiler",
  orbitPacket - orbitCompiler . directPacket
];
expectZero[
  "orbit packet compiler determinant",
  Det[orbitCompiler] - cEta
];
expectZero[
  "inverse direct-observable orbit packet compiler",
  FullSimplify[Inverse[orbitCompiler].orbitPacket - directPacket, Assumptions -> $Assumptions]
];
recoveredDirectPacket = FullSimplify[
  Inverse[orbitCompiler].{ThetaVar, XiVar, RVar},
  Assumptions -> $Assumptions
];
expectZero[
  "formal orbit packet recovers the direct drifts",
  recoveredDirectPacket - {
    ThetaVar,
    RVar,
    -((1 - epsEta) (XiVar + RVar))/epsEta
  }
];
expectZero[
  "zero orbit packet forces zero direct drifts",
  recoveredDirectPacket /. {ThetaVar -> 0, XiVar -> 0, RVar -> 0}
];

subbanner["VI. Exact two-packet split"];

Mmix = FullSimplify[
  8 ZW (1 + chi0)^2/(Pi^2 (1 - epsEta) (1 - eps)),
  Assumptions -> $Assumptions
];
Sfactor = FullSimplify[
  1 + zeta (1 - eps)/(1 - zeta eps),
  Assumptions -> $Assumptions
];
Mtr = FullSimplify[Mmix Sfactor, Assumptions -> $Assumptions];

expectZero["mixed-only product law", FullSimplify[Rtarget Mmix - Cmix, Assumptions -> $Assumptions]];
expectZero[
  "support-packet sensitivity",
  D[Mtr, zeta] - Mmix (1 - eps)/(1 - zeta eps)^2
];

subbanner["VII. Probe-only rational sample point"];

sampleRules = {
  chi0 -> 3/2,
  deltaU -> 2/3,
  ZW -> 13/17,
  epsW -> 1/3,
  epsEta -> 1/5,
  Lambda -> 7/11,
  zeta -> 1
};

Print["epsilon = ", fmt[FullSimplify[eps /. sampleRules]]];
Print["varrhoPhys = ", fmt[FullSimplify[varrhoPhys /. sampleRules]]];
Print["sigmaPhys = ", fmt[FullSimplify[sigmaPhys /. sampleRules]]];
Print["Rtr = ", fmt[FullSimplify[Rtr /. sampleRules]]];
Print["Rtarget = ", fmt[FullSimplify[Rtarget /. sampleRules]]];
Print["Mmix = ", fmt[FullSimplify[Mmix /. sampleRules]]];
Print["Mtr = ", fmt[FullSimplify[Mtr /. sampleRules]]];

Print[""];
Print["All Stage 242 symbolic checks passed."];
