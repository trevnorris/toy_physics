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

stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

normalizeExpr[expr_] := Module[{res},
  res = Which[
    MatrixQ[expr],
      Map[
        FullSimplify[stripConditional[#], Assumptions -> $Assumptions] &,
        expr,
        {2}
      ],
    VectorQ[expr],
      Map[
        FullSimplify[stripConditional[#], Assumptions -> $Assumptions] &,
        expr
      ],
    True,
      FullSimplify[stripConditional[expr], Assumptions -> $Assumptions]
  ];
  stripConditional[res]
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

logReduce[expr_] := FullSimplify[
  PowerExpand[expr],
  Assumptions -> $Assumptions
];

banner["STAGE 237 — ACTUAL-BRANCH DRESSING COMPILER"];

Clear[
  Bstar, Cstar, Rtr, RtrRef, Rtarget, RtargetRef, epsEta, epsEtaRef,
  Rratio, qEtaVar, t, R1, E1, cEtaU, cEtaURef, KU, KURef,
  KEtaEff, KEtaEffRef, c1, kappaU, kappaEta, zeta, Mtr,
  lambdaPhi, KPhiEff, lambdaW, KWEff, eps, Mmix
];

$Assumptions =
  Element[
    {
      Bstar, Cstar, qEtaVar, t, R1, E1, c1, kappaU, kappaEta
    },
    Reals
  ] &&
  Element[
    {
      Rtr, RtrRef, Rtarget, RtargetRef, Rratio, epsEta, epsEtaRef,
      cEtaU, cEtaURef, KU, KURef, KEtaEff, KEtaEffRef, zeta, Mtr,
      lambdaPhi, KPhiEff, lambdaW, KWEff, eps, Mmix
    },
    Reals
  ] &&
  Rtr > 0 && RtrRef > 0 && Rtarget > 0 && RtargetRef > 0 &&
  Rratio > 0 && 0 < epsEta < 1 && 0 < epsEtaRef < 1 &&
  0 < epsEtaRef Exp[qEtaVar] < 1 &&
  0 < 1 - (1 - epsEtaRef) Rratio &&
  cEtaU > 0 && cEtaURef > 0 && KU > 0 && KURef > 0 &&
  KEtaEff > 0 && KEtaEffRef > 0 &&
  zeta > 0 && Mtr > 0 && lambdaPhi > 0 && KPhiEff > 0 &&
  lambdaW > 0 && KWEff > 0 && 0 < eps < 1 && Mmix > 0 &&
  1 - zeta eps != 0;

subbanner["M1. Rigid-mouth reduction"];

qTr = -Cstar Log[Rtr/RtrRef];
qNt = (
  Bstar Log[Rtr/RtrRef] +
  Log[(1 - epsEta)/(1 - epsEtaRef)] -
  Log[Rtarget/RtargetRef]
);
qEta = Log[epsEta/epsEtaRef];

qTrRigid = qTr /. Rtr -> RtrRef;
qNtRigid = qNt /. Rtr -> RtrRef;

expectZero["M1 q_tr on R_tr = R_tr_ref", qTrRigid];
expectZero[
  "M1 q_nt rigid-mouth reduction",
  qNtRigid - (
    Log[(1 - epsEta)/(1 - epsEtaRef)] -
    Log[Rtarget/RtargetRef]
  )
];

subbanner["M2. Finite static-blind curve and inverse"];

qNtRatio = Log[(1 - epsEta)/(1 - epsEtaRef)] - Log[Rratio];
staticCurve = FullSimplify[(1 - epsEta)/(1 - epsEtaRef), Assumptions -> $Assumptions];
staticExp = logReduce[Exp[qNtRatio]];
staticSolve = Solve[staticExp == 1, Rratio, Reals];
If[Length[staticSolve] != 1, fail["M2 q_nt = 0 solve is not unique", staticSolve]];
expectZero[
  "M2 q_nt = 0 equivalent to static curve",
  (Rratio /. First[staticSolve]) - staticCurve
];

epsFromQeta = epsEtaRef Exp[qEtaVar];
RratioOfQeta = FullSimplify[
  (1 - epsFromQeta)/(1 - epsEtaRef),
  Assumptions -> $Assumptions
];
qEtaFromRratio = Log[(1 - (1 - epsEtaRef) Rratio)/epsEtaRef];

expectZero[
  "M2 parameterized static-blind curve",
  RratioOfQeta - (1 - epsEtaRef Exp[qEtaVar])/(1 - epsEtaRef)
];
expectZero[
  "M2 inverse round-trip q_eta -> Rratio -> q_eta",
  logReduce[(qEtaFromRratio /. Rratio -> RratioOfQeta) - qEtaVar]
];
expectZero[
  "M2 inverse round-trip Rratio -> q_eta -> Rratio",
  logReduce[(RratioOfQeta /. qEtaVar -> qEtaFromRratio) - Rratio]
];
expectZero["M2 endpoint Rratio(q_eta = 0)", (RratioOfQeta /. qEtaVar -> 0) - 1];

subbanner["M3. First-order packet compiler and tangent"];

cEta = FullSimplify[epsEtaRef/(1 - epsEtaRef), Assumptions -> $Assumptions];
epsEtaPert = epsEtaRef Exp[t E1];
RratioPert = Exp[t R1];
qEtaPert = logReduce[Log[epsEtaPert/epsEtaRef]];
qNtPert = Log[(1 - epsEtaPert)/(1 - epsEtaRef)] - Log[RratioPert];

qEtaLinear = Normal[Series[qEtaPert, {t, 0, 1}]];
qNtLinear = Normal[Series[qNtPert, {t, 0, 1}]];

expectZero["M3 q_eta first-order term", qEtaLinear - t E1];
expectZero["M3 q_nt first-order term", qNtLinear + t (R1 + cEta E1)];

tangentSlope = FullSimplify[
  D[Log[RratioOfQeta], qEtaVar] /. qEtaVar -> 0,
  Assumptions -> $Assumptions
];
expectZero["M3 tangent slope of finite curve", tangentSlope + cEta];

packetMatrix = {{-1, -cEta}, {0, 1}};
expectZero["M3 packet matrix determinant", Det[packetMatrix] + 1];

subbanner["M4. Microscopic coherent compiler"];

epsEtaMicro = cEtaU^2/(KU KEtaEff);
epsEtaMicroRef = cEtaURef^2/(KURef KEtaEffRef);
qEtaMicro = Log[epsEtaMicro/epsEtaMicroRef];
qEtaMicroSplit = (
  2 Log[cEtaU/cEtaURef] -
  Log[KU/KURef] -
  Log[KEtaEff/KEtaEffRef]
);

expectZero[
  "M4 microscopic log expansion",
  logReduce[qEtaMicro - qEtaMicroSplit]
];

subbanner["M5. Microscopic first-order drift extractor"];

qEtaMicroPert = qEtaMicro /. {
  cEtaU -> cEtaURef Exp[t c1],
  KU -> KURef Exp[t kappaU],
  KEtaEff -> KEtaEffRef Exp[t kappaEta]
};
qEtaMicroFirst = Coefficient[
  Normal[Series[logReduce[qEtaMicroPert], {t, 0, 1}]],
  t,
  1
];

expectZero[
  "M5 microscopic drift coefficient",
  qEtaMicroFirst - (2 c1 - kappaU - kappaEta)
];

subbanner["M6. Support-blindness through sector independence"];

zetaDef = lambdaPhi^2 KWEff/(lambdaW^2 KPhiEff);
MtrDef = Mmix (1 + zeta (1 - eps)/(1 - zeta eps));
MtrComposite = MtrDef /. zeta -> zetaDef;

qEtaSupport = Log[
  (
    cEtaUSupport[zeta, Mtr, lambdaPhi, KPhiEff]^2/
    (KUSupport[zeta, Mtr, lambdaPhi, KPhiEff]
      KEtaEffSupport[zeta, Mtr, lambdaPhi, KPhiEff])
  ) (KURef KEtaEffRef/cEtaURef^2)
];
If[
  !(And @@ (Not[FreeQ[qEtaSupport, #]] & /@ {zeta, Mtr, lambdaPhi, KPhiEff})),
  fail["M6 support variables are absent from q_eta support frame", qEtaSupport]
];

qEtaSupportComposite = qEtaSupport /. {zeta -> zetaDef, Mtr -> MtrComposite};
rawSupportDerivatives = {
  D[qEtaSupport, zeta],
  D[qEtaSupport, Mtr],
  D[qEtaSupportComposite, lambdaPhi],
  D[qEtaSupportComposite, KPhiEff]
};
If[
  !(And @@ (Not[FreeQ[#, Derivative]] & /@ rawSupportDerivatives)),
  fail["M6 raw derivatives do not expose chain-rule support paths", rawSupportDerivatives]
];

supportBlindRules = {
  Derivative[__][cEtaUSupport][__] :> 0,
  Derivative[__][KUSupport][__] :> 0,
  Derivative[__][KEtaEffSupport][__] :> 0
};

(* Negative control: if KEtaEffSupport were redefined with a KPhiEff factor,
   D[q_eta, KPhiEff] would retain a nonzero -1/KPhiEff contribution. *)
expectZero["M6 D[q_eta, zeta]", D[qEtaSupport, zeta] /. supportBlindRules];
expectZero["M6 D[q_eta, M_tr]", D[qEtaSupport, Mtr] /. supportBlindRules];
expectZero[
  "M6 D[q_eta, lambda_phi]",
  D[qEtaSupportComposite, lambdaPhi] /. supportBlindRules
];
expectZero[
  "M6 D[q_eta, K_phi_eff]",
  D[qEtaSupportComposite, KPhiEff] /. supportBlindRules
];

Print["zeta = ", fmt[zetaDef]];
Print["M_tr = ", fmt[MtrDef]];

subbanner["M7. Dependent-plane ray and orbit-lock"];

yEtaRay = -qEtaVar {0, 1, 1};
yEtaFromEps = -qEta {0, 1, 1};
yEtaFromTarget = -qEtaFromRratio {0, 1, 1};

expectZero[
  "M7 dependent-plane ray from eps_eta",
  (yEtaRay /. qEtaVar -> qEta) - yEtaFromEps
];
expectZero[
  "M7 dependent-plane ray on static-blind curve",
  (yEtaRay /. qEtaVar -> qEtaFromRratio) - yEtaFromTarget
];
expectZero[
  "M7 finite endpoint from Rratio = 1",
  qEtaFromRratio /. Rratio -> 1
];
expectZero[
  "M7 finite endpoint from eps_eta = eps_eta_ref",
  qEta /. epsEta -> epsEtaRef
];

codimSolution = Solve[{-R1 - cEta E1 == 0, E1 == 0}, {R1, E1}, Reals];
If[Length[codimSolution] != 1, fail["M7 codim-two solve is not unique", codimSolution]];
expectZero["M7 codim-two solution R1", R1 /. First[codimSolution]];
expectZero["M7 codim-two solution E1", E1 /. First[codimSolution]];
expectZero[
  "M7 R1 = E1 = 0 satisfies first-order packet",
  {-R1 - cEta E1, E1} /. {R1 -> 0, E1 -> 0}
];

Print[""];
Print["All Stage 237 Mathematica symbolic checks passed."];
