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

stripConditional[expr_] := (
  expr
  /. ConditionalExpression[e_, _] :> e
  /. Piecewise[{{e_, _}}, Indeterminate] :> e
);

cleanResidual[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

cleanLogResidual[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[PowerExpand[res], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanResidual[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectLogZero[name_String, expr_] := Module[{res},
  res = cleanLogResidual[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  res = stripConditional[res];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

linCoeff[expr_] := FullSimplify[
  Coefficient[Normal[Series[expr, {eps, 0, 1}]], eps, 1],
  Assumptions -> $Assumptions
];

$Assumptions = (
  cStar > 0 && bStar > 0 &&
  rtrRef > 0 && rtargetRef > 0 &&
  epsEtaRef > 0 && epsEtaStar > 0 &&
  rtr > 0 && rtarget > 0 && epsEta > 0 &&
  Element[{qTr, qNt, qEta, r1, R1, E1, eps, xi}, Reals] &&
  Element[vareps, Reals] && vareps != 0
);

banner["Stage 234 Mathematica audit: direct branch-observable static gate"];

cEta = epsEtaStar/(1 - epsEtaStar);

qTrChart[rad_] := -cStar*Log[rad/rtrRef];
qNtChart[rad_, dressing_, target_, refDressing_] := (
  bStar*Log[rad/rtrRef]
  + Log[(1 - dressing)/(1 - refDressing)]
  - Log[target/rtargetRef]
);
qEtaChart[dressing_, refDressing_] := Log[dressing/refDressing];

subbanner["M1. Exact finite quotient chart and inverse roundtrip"];

rtrInv = rtrRef*Exp[-qTr/cStar];
epsEtaInv = epsEtaRef*Exp[qEta];
rtargetInv = rtargetRef*Exp[-qNt - (bStar/cStar)*qTr]*(1 - epsEtaInv)/(1 - epsEtaRef);

expectLogZero["M1 q_tr roundtrip", qTrChart[rtrInv] - qTr];
expectLogZero["M1 q_nt roundtrip", qNtChart[rtrInv, epsEtaInv, rtargetInv, epsEtaRef] - qNt];
expectLogZero["M1 q_eta roundtrip", qEtaChart[epsEtaInv, epsEtaRef] - qEta];

Print["q_tr chart = ", fmt[qTrChart[rtr]]];
Print["q_nt chart = ", fmt[qNtChart[rtr, epsEta, rtarget, epsEtaRef]]];
Print["q_eta chart = ", fmt[qEtaChart[epsEta, epsEtaRef]]];

subbanner["M2. First-order linearization by Series coefficients"];

rtrPert = rtrRef*Exp[eps*r1];
rtargetPert = rtargetRef*Exp[eps*R1];
epsEtaPert = epsEtaStar*Exp[eps*E1];

qTrSeriesExpr = qTrChart[rtrPert];
qNtSeriesExpr = qNtChart[rtrPert, epsEtaPert, rtargetPert, epsEtaStar];
qEtaSeriesExpr = qEtaChart[epsEtaPert, epsEtaStar];

qTr1 = linCoeff[qTrSeriesExpr];
qNt1 = linCoeff[qNtSeriesExpr];
qEta1 = linCoeff[qEtaSeriesExpr];

expectZero["M2 q_tr^(1) + C_star r1", qTr1 + cStar*r1];
expectZero["M2 q_nt^(1) - expected", qNt1 - (bStar*r1 - cEta*E1 - R1)];
expectZero["M2 q_eta^(1) - E1", qEta1 - E1];

Print["q_tr^(1) = ", fmt[qTr1]];
Print["q_nt^(1) = ", fmt[FullSimplify[qNt1, Assumptions -> $Assumptions]]];
Print["q_eta^(1) = ", fmt[qEta1]];

subbanner["M3. Triangular compiler and inverse drift map"];

theta1 = FullSimplify[-qTr1/cStar, Assumptions -> $Assumptions];
xi1FromCoefficients = FullSimplify[qNt1 + (bStar/cStar)*qTr1, Assumptions -> $Assumptions];
rCal1 = FullSimplify[-xi1FromCoefficients - cEta*qEta1, Assumptions -> $Assumptions];
e1Inverse = FullSimplify[-(1 - epsEtaStar)/epsEtaStar*(rCal1 + xi1FromCoefficients), Assumptions -> $Assumptions];

expectZero["M3 Theta1 - r1", theta1 - r1];
expectZero["M3 Xi1 + R1 + c_eta E1", xi1FromCoefficients + R1 + cEta*E1];
expectZero["M3 Rcal1 - R1", rCal1 - R1];
expectZero["M3 inverse E1 map", e1Inverse - E1];

Print["Theta1 = ", fmt[theta1]];
Print["Xi1 = ", fmt[xi1FromCoefficients]];
Print["Rcal1 = ", fmt[rCal1]];

subbanner["M4. Exact R_tr cancellation"];

xi1Raw = qNt1 + (bStar/cStar)*qTr1;
expectZero["M4 D[Xi1, r1]", D[xi1Raw, r1]];
expectTrue["M4 simplified Xi1 contains no r1", FreeQ[FullSimplify[xi1Raw, Assumptions -> $Assumptions], r1]];

subbanner["M5. Rigid-mouth reduction"];

rigidRule = {r1 -> 0};
qNtRigidFinite = qNtChart[rtrRef, epsEta, rtarget, epsEtaRef];
qNtRigidFiniteExpected = Log[(1 - epsEta)/(1 - epsEtaRef)] - Log[rtarget/rtargetRef];

expectZero["M5 Theta1 rigid", theta1 /. rigidRule];
expectZero["M5 Xi1 rigid equals q_nt^(1) rigid", (xi1FromCoefficients /. rigidRule) - (qNt1 /. rigidRule)];
expectLogZero["M5 finite rigid-mouth q_nt", qNtRigidFinite - qNtRigidFiniteExpected];

Print["Rigid-mouth finite q_nt = ", fmt[qNtRigidFinite]];

subbanner["M6. Two-observable strip form and ceiling ordering"];

robust = 367930328492646/1000000000000000; (* 0.367930328492646 *)
nonempty = 737619063660757/1000000000000000; (* 0.737619063660757 *)
xi1TwoObservable = -R1 - cEta*E1;

robustLower = -cEta*E1 - robust/Abs[vareps];
robustUpper = -cEta*E1 + robust/Abs[vareps];
nonemptyLower = -cEta*E1 - nonempty/Abs[vareps];
nonemptyUpper = -cEta*E1 + nonempty/Abs[vareps];

expectZero["M6 robust lower edge", (xi1TwoObservable /. R1 -> robustLower) - robust/Abs[vareps]];
expectZero["M6 robust upper edge", (xi1TwoObservable /. R1 -> robustUpper) + robust/Abs[vareps]];
expectZero["M6 nonempty lower edge", (xi1TwoObservable /. R1 -> nonemptyLower) - nonempty/Abs[vareps]];
expectZero["M6 nonempty upper edge", (xi1TwoObservable /. R1 -> nonemptyUpper) + nonempty/Abs[vareps]];
expectTrue["M6 0 < robust < nonempty", 0 < robust < nonempty];

Print["Robust strip = ", fmt[{robustLower, robustUpper}]];
Print["Nonempty strip = ", fmt[{nonemptyLower, nonemptyUpper}]];

subbanner["M7. Canonical direct-branch families"];

expectZero["M7 pure-target Xi1 - xi", (xi1TwoObservable /. {R1 -> -xi, E1 -> 0}) - xi];
expectZero["M7 pure-dressing Xi1 - xi", (xi1TwoObservable /. {R1 -> 0, E1 -> -xi/cEta}) - xi];

balancedMinimum = FullSimplify[
  Minimize[{R1^2 + E1^2, -R1 - cEta*E1 == xi}, {R1, E1}, Reals],
  Assumptions -> $Assumptions
];
balancedRules = stripConditional[balancedMinimum[[2]]];
balancedR1 = stripConditional[R1 /. balancedRules];
balancedE1 = stripConditional[E1 /. balancedRules];
expectedBalancedR1 = -xi/(1 + cEta^2);
expectedBalancedE1 = -cEta*xi/(1 + cEta^2);

expectZero["M7 balanced R1 minimizer", balancedR1 - expectedBalancedR1];
expectZero["M7 balanced E1 minimizer", balancedE1 - expectedBalancedE1];
expectZero["M7 balanced Xi1 - xi", (xi1TwoObservable /. {R1 -> balancedR1, E1 -> balancedE1}) - xi];

Print["Balanced Minimize result = ", fmt[balancedMinimum]];
Print["Balanced expected family = ", fmt[{expectedBalancedR1, expectedBalancedE1}]];

Print[""];
Print["All Stage 234 Mathematica checks passed."];
