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

stripConditional[expr_] := expr /. ConditionalExpression[value_, _] :> value;

normalizeExpr[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

expectZero[name_String, expr_] := Module[{res},
  res = normalizeExpr[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = stripConditional[FullSimplify[cond, Assumptions -> $Assumptions]];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === True], pass[name], fail[name, res]];
];

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = N[Abs[value - target], 20];
  Print[name, " diff = ", fmt[diff], " (tol = ", fmt[tol], ")"];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 250 - CROSSING VS COLLAPSE GOLDILOCKS WINDOW COMPILER"];

Clear[En, Vpeak, V0, v0, x, ms, lambdaEff, gUV, chiPeak, muEta, alpha];

paramsPositive =
  ms > 0 && lambdaEff > 0 && gUV > 0 && chiPeak > 0 && muEta > 0 && alpha > 0;

$Assumptions =
  Element[{En, Vpeak, V0, v0, x, ms, lambdaEff, gUV, chiPeak, muEta, alpha}, Reals] &&
  paramsPositive && x > 0;

subbanner["M1. Crossing-time form and global monotonicity"];

tcrossGap = FullSimplify[lambdaEff Sqrt[ms/(2 x)], Assumptions -> $Assumptions];
tcross = FullSimplify[tcrossGap /. x -> En - Vpeak, Assumptions -> $Assumptions && En > Vpeak];
dtcrossDE = FullSimplify[D[tcross, En], Assumptions -> $Assumptions && En > Vpeak];

expectZero["M1 t_cross form", tcross - lambdaEff Sqrt[ms/(2 (En - Vpeak))]];
Print["M1 dt_cross/dE = ", fmt[dtcrossDE]];

m1Monotone = Resolve[
  ForAll[{En, Vpeak, ms, lambdaEff},
    Implies[
      ms > 0 && lambdaEff > 0 && En > Vpeak,
      D[lambdaEff Sqrt[ms/(2 (En - Vpeak))], En] < 0
    ]
  ],
  Reals
];
m1GapMonotone = Resolve[
  ForAll[{x, ms, lambdaEff},
    Implies[
      ms > 0 && lambdaEff > 0 && x > 0,
      D[lambdaEff Sqrt[ms/(2 x)], x] < 0
    ]
  ],
  Reals
];
expectTrue["M1 global dt_cross/dE < 0", m1Monotone];
expectTrue["M1 global dt_cross/dE < 0 for x>0", m1GapMonotone];

subbanner["M2. Collapse-time compiler"];

tcollapse = FullSimplify[Sqrt[muEta/(gUV chiPeak)], Assumptions -> $Assumptions];
gammaColl = FullSimplify[1/tcollapse, Assumptions -> $Assumptions];

Print["M2 t_collapse = ", fmt[tcollapse]];
expectZero["M2 Gamma_coll identity", gammaColl - Sqrt[gUV chiPeak/muEta]];

subbanner["M3. Survival ratio and edge from Solve/Reduce"];

Sratio = FullSimplify[tcross/tcollapse, Assumptions -> $Assumptions && En > Vpeak];
Eedge = FullSimplify[Vpeak + ms gUV chiPeak lambdaEff^2/(2 muEta), Assumptions -> $Assumptions];
edgeSolve = Solve[Sratio^2 == 1, En, Reals];
edgeFromSolve = FullSimplify[stripConditional[En /. First[edgeSolve]], Assumptions -> $Assumptions];
edgeReduce = Reduce[
  Sratio^2 == 1 && En > Vpeak && ms > 0 && lambdaEff > 0 &&
    gUV > 0 && chiPeak > 0 && muEta > 0,
  En,
  Reals
];
m3Unique = Resolve[
  ForAll[{En, Vpeak, ms, lambdaEff, gUV, chiPeak, muEta},
    Implies[
      ms > 0 && lambdaEff > 0 && gUV > 0 && chiPeak > 0 && muEta > 0 &&
        En > Vpeak && Sratio^2 == 1,
      En == Eedge
    ]
  ],
  Reals
];

Print["M3 S(E) = ", fmt[Sratio]];
Print["M3 Reduce edge = ", fmt[edgeReduce]];
expectZero["M3 E_edge from Solve", edgeFromSolve - Eedge];
expectTrue["M3 unique positive-domain edge", m3Unique];

subbanner["M4. Global one-sided Goldilocks window"];

windowEquality = (Sratio < 1) == (En > Eedge);
windowEquivalence = Equivalent[Sratio < 1, En > Eedge];
m4Window = Resolve[
  ForAll[{En, Vpeak, ms, lambdaEff, gUV, chiPeak, muEta},
    Implies[
      ms > 0 && lambdaEff > 0 && gUV > 0 && chiPeak > 0 && muEta > 0 &&
        En > Vpeak,
      windowEquivalence
    ]
  ],
  Reals
];
expectTrue["M4 S(E)<1 iff E>E_edge globally", m4Window];

subbanner["M5. Heavy-throat cancellation"];

EedgeHeavy = FullSimplify[Eedge /. muEta -> alpha ms, Assumptions -> $Assumptions];

expectZero["M5 dE_edge/dm_s after mu_eta=alpha m_s", D[EedgeHeavy, ms]];
expectZero[
  "M5 heavy-throat edge value",
  EedgeHeavy - (Vpeak + gUV chiPeak lambdaEff^2/(2 alpha))
];

subbanner["M6. Speed-space identity"];

Elaunch = FullSimplify[1/2 ms v0^2 + V0, Assumptions -> $Assumptions];
vcritSq = FullSimplify[2 (Vpeak - V0)/ms, Assumptions -> $Assumptions];
speedResidual = FullSimplify[
  2 (Eedge - V0)/ms - (vcritSq + lambdaEff^2 gUV chiPeak/muEta),
  Assumptions -> $Assumptions
];

Print["M6 E_launch = ", fmt[Elaunch]];
expectZero["M6 v_safe,min^2 identity", speedResidual];

subbanner["M7. Session-III benchmark numerics"];

VpeakNum = 3.42933112;
V0Num = 0.19999794;
lambdaEffNum = 0.42826825;
gUVNum = 0.95;
chiPeakNum = 21.73204372;
msNum = 1836.15267343;
muEtaNum = 1836.15267343;
lambdaRawNum = 1;
chiRawNum = 50.74399964;

benchRules = {
  Vpeak -> VpeakNum,
  V0 -> V0Num,
  lambdaEff -> lambdaEffNum,
  gUV -> gUVNum,
  chiPeak -> chiPeakNum,
  ms -> msNum,
  muEta -> muEtaNum
};
rawRules = {
  Vpeak -> VpeakNum,
  V0 -> V0Num,
  lambdaEff -> lambdaRawNum,
  gUV -> gUVNum,
  chiPeak -> chiRawNum,
  ms -> msNum,
  muEta -> muEtaNum
};

tcollapseNum = N[tcollapse /. benchRules, 20];
EedgeNum = N[Eedge /. benchRules, 20];
vcritNum = N[Sqrt[2 (VpeakNum - V0Num)/msNum], 20];
vsafeNum = N[Sqrt[2 (EedgeNum - V0Num)/msNum], 20];
ratioNum = N[vsafeNum/vcritNum, 20];
SedgeNum = N[Sratio /. Join[benchRules, {En -> EedgeNum}], 20];
tcollapseRawNum = N[tcollapse /. rawRules, 20];
EedgeRawNum = N[Eedge /. rawRules, 20];

expectApprox["M7 t_collapse", tcollapseNum, 9.43066476, 5*10^-9];
expectApprox["M7 E_safe,min", EedgeNum, 5.32265943, 5*10^-8];
expectApprox["M7 v_safe,min", vsafeNum, 0.07469791, 5*10^-9];
expectApprox["M7 v_safe/v_crit", ratioNum, 1.25948037, 5*10^-9];
expectApprox["M7 t_collapse_raw", tcollapseRawNum, 6.17163516, 5*10^-9];
expectApprox["M7 E_safe,min_raw", EedgeRawNum, 27.53273095, 5*10^-8];
expectApprox["M7 S(E_edge)", SedgeNum, 1, 5*10^-12];
expectTrue["M7 raw edge above relaxed edge", EedgeRawNum > EedgeNum];

Print[""];
Print["All Stage 250 Mathematica checks passed."];
