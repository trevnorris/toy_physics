ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

clean[expr_] := Module[{res},
  res = FullSimplify[Together[expr], Assumptions -> $Assumptions];
  res = res /. ConditionalExpression[e_, _] :> e;
  FullSimplify[res, Assumptions -> $Assumptions]
];

expectZero[expr_, label_] := Module[{res},
  res = clean[expr];
  Print[label, " residual = ", fmt[res]];
  If[
    TrueQ[FullSimplify[res == 0, Assumptions -> $Assumptions]],
    Print["[ok] ", label],
    Print["FAIL: ", label];
    Exit[1]
  ];
];

expectNonZero[expr_, label_] := Module[{res},
  res = clean[expr];
  Print[label, " witness = ", fmt[res]];
  If[
    TrueQ[FullSimplify[res == 0, Assumptions -> $Assumptions]],
    Print["FAIL: ", label];
    Exit[1],
    Print["[ok] ", label]
  ];
];

expectSupportFree[expr_, label_] := Module[{free},
  free = FreeQ[expr, zeta] && FreeQ[expr, Mmix];
  Print[label, " support-free = ", free];
  If[TrueQ[free], Print["[ok] ", label], Print["FAIL: ", label]; Exit[1]];
];

Clear[
  chi0, deltaU, ZW, OmegaW2, eps, epsEta, Lambda0, Bstar, Cstar,
  RtrRef, T2ref, epsEtaRef, RtargetRef, zeta, Mmix, h, dlnchi,
  dlndelta, dlnZW, dlnOm2, dlneps, dlnepsEta
];

$Assumptions = (
  Element[
    {
      chi0, deltaU, ZW, OmegaW2, eps, epsEta, Lambda0, Bstar, Cstar,
      RtrRef, T2ref, epsEtaRef, RtargetRef, zeta, Mmix, h, dlnchi,
      dlndelta, dlnZW, dlnOm2, dlneps, dlnepsEta
    },
    Reals
  ] &&
  chi0 > 0 && deltaU > 0 && ZW > 0 && OmegaW2 > 0 &&
  0 < eps < 1 && 0 < epsEta < 1 && 0 < epsEtaRef < 1 &&
  Lambda0 > 0 && RtrRef > 0 && T2ref > 0 && RtargetRef > 0 &&
  Cstar != 0 && Mmix != 0 && 1 - zeta eps != 0
);

Print["Stage 238 Mathematica audit"];

Rtr = (1 + chi0/(1 + deltaU))/(1 + chi0);
T2 = ZW (1 + chi0)^2/(OmegaW2 (1 - eps)^2);
Rtarget = Lambda0 OmegaW2 (1 - epsEta) (1 - eps)^2/(ZW (1 + chi0)^2);

expectZero[
  Rtarget T2 - Lambda0 (1 - epsEta),
  "M1 coherent transfer-shape identity"
];

qTr = -Cstar Log[Rtr/RtrRef];
qNt = (
  Bstar Log[Rtr/RtrRef]
  + Log[(1 - epsEta)/(1 - epsEtaRef)]
  - Log[Rtarget/RtargetRef]
);
qEta = Log[epsEta/epsEtaRef];

Mtr = Mmix (1 + zeta (1 - eps)/(1 - zeta eps));
expectZero[
  D[Mtr, zeta] - Mmix (1 - eps)/(1 - zeta eps)^2,
  "M2 support channel derivative in zeta"
];
expectNonZero[D[Mtr, zeta], "M2 negative control zeta live"];
expectNonZero[D[Mtr, Mmix], "M2 negative control Mmix live"];

RtrLeak = Rtr Mtr/Mmix;
expectNonZero[D[RtrLeak, zeta], "M2 leak detector catches Rtr support contamination"];

qNtFactored = Log[T2/T2ref] - (Bstar/Cstar) qTr;

expectSupportFree[Rtr, "M2 structural exclusion R_tr"];
expectSupportFree[T2, "M2 structural exclusion T^2"];
expectSupportFree[epsEta, "M2 structural exclusion eps_eta"];
expectSupportFree[qTr, "M2 structural exclusion q_tr"];
expectSupportFree[qNtFactored, "M2 structural exclusion q_nt"];
expectSupportFree[qEta, "M2 structural exclusion q_eta"];

expectZero[D[Rtr, zeta], "M2 support-blind R_tr zeta derivative"];
expectZero[D[T2, zeta], "M2 support-blind T^2 zeta derivative"];
expectZero[D[epsEta, zeta], "M2 support-blind eps_eta zeta derivative"];
expectZero[D[qTr, zeta], "M2 support-blind q_tr zeta derivative"];
expectZero[D[qNtFactored, zeta], "M2 support-blind q_nt zeta derivative"];
expectZero[D[qEta, zeta], "M2 support-blind q_eta zeta derivative"];
expectZero[D[Rtr, Mmix], "M2 support-blind R_tr Mmix derivative"];
expectZero[D[T2, Mmix], "M2 support-blind T^2 Mmix derivative"];
expectZero[D[epsEta, Mmix], "M2 support-blind eps_eta Mmix derivative"];
expectZero[D[qTr, Mmix], "M2 support-blind q_tr Mmix derivative"];
expectZero[D[qNtFactored, Mmix], "M2 support-blind q_nt Mmix derivative"];
expectZero[D[qEta, Mmix], "M2 support-blind q_eta Mmix derivative"];

qFactor = qNt + (Bstar/Cstar) qTr;
qFactorResidual = PowerExpand[
  qFactor - Log[T2/T2ref] /. RtargetRef -> Lambda0 (1 - epsEtaRef)/T2ref
];
expectZero[qFactorResidual, "M3 finite packet log factorization"];

qNtRigid = PowerExpand[
  Log[(1 - epsEta)/(1 - epsEtaRef)]
  - Log[Rtarget/(Lambda0 (1 - epsEtaRef)/T2ref)]
];
expectZero[
  qNtRigid - Log[T2/T2ref],
  "M3 rigid-mouth q_nt transfer-shape reduction"
];
expectZero[
  qEta - Log[epsEta/epsEtaRef],
  "M3 rigid-mouth q_eta dressing reduction"
];

logRtr = Log[1 + chi0/(1 + deltaU)] - Log[1 + chi0];
logT2 = Log[ZW] + 2 Log[1 + chi0] - Log[OmegaW2] - 2 Log[1 - eps];
logRtarget = (
  Log[Lambda0] + Log[OmegaW2] + Log[1 - epsEta] + 2 Log[1 - eps]
  - Log[ZW] - 2 Log[1 + chi0]
);

perturbCore = {
  chi0 -> chi0 Exp[h dlnchi],
  deltaU -> deltaU Exp[h dlndelta],
  ZW -> ZW Exp[h dlnZW],
  OmegaW2 -> OmegaW2 Exp[h dlnOm2],
  eps -> eps Exp[h dlneps],
  epsEta -> epsEta Exp[h dlnepsEta]
};

dlnRtrCalc = clean[D[logRtr /. perturbCore, h] /. h -> 0];
dlnT2Calc = clean[D[logT2 /. perturbCore, h] /. h -> 0];
dlnRtargetCalc = clean[D[logRtarget /. perturbCore, h] /. h -> 0];

trackingCondition = (1 + deltaU) dlnchi + (1 + chi0) dlndelta;
prefactor = chi0 deltaU/((1 + chi0) (1 + deltaU) (1 + chi0 + deltaU));
dlnRtrExpected = -prefactor trackingCondition;
dlnT2Expected = (
  dlnZW - dlnOm2 + 2 chi0/(1 + chi0) dlnchi + 2 eps/(1 - eps) dlneps
);
dlnRtargetExpected = (
  dlnOm2 - dlnZW - 2 chi0/(1 + chi0) dlnchi
  - 2 eps/(1 - eps) dlneps - epsEta/(1 - epsEta) dlnepsEta
);

expectZero[dlnRtrCalc - dlnRtrExpected, "M4 first-order delta ln R_tr compiler"];
expectZero[dlnT2Calc - dlnT2Expected, "M4 first-order delta ln T^2 compiler"];
expectZero[
  dlnRtargetCalc - dlnRtargetExpected,
  "M4 first-order delta ln R_target compiler"
];

qTrFirst = -Cstar dlnRtrCalc;
qNtFirst = Bstar dlnRtrCalc - dlnRtargetCalc - epsEta/(1 - epsEta) dlnepsEta;
qEtaFirst = clean[D[Log[epsEta Exp[h dlnepsEta]], h] /. h -> 0];
qNtFirstRigid = -dlnRtargetCalc - epsEta/(1 - epsEta) dlnepsEta;

expectZero[
  qNtFirst + (Bstar/Cstar) qTrFirst - dlnT2Calc,
  "M5 first-order packet relation"
];
expectZero[
  qNtFirstRigid - dlnT2Calc,
  "M5 rigid first-order q_nt transfer-shape reduction"
];
expectZero[qEtaFirst - dlnepsEta, "M5 rigid first-order q_eta dressing compiler"];

expectZero[
  dlnRtrCalc + prefactor trackingCondition,
  "M6 tracking gate factorization from derived drift"
];
expectZero[
  dlnRtrCalc /. dlndelta -> -(1 + deltaU)/(1 + chi0) dlnchi,
  "M6 tracking gate locus kills derived drift"
];

expectZero[
  PowerExpand[Exp[qNtRigid] - T2/T2ref],
  "M7 finite transfer-shape gate"
];
expectZero[
  PowerExpand[Exp[qEta] - epsEta/epsEtaRef],
  "M7 finite dressing gate"
];

Print["All Stage 238 Mathematica checks passed."];
