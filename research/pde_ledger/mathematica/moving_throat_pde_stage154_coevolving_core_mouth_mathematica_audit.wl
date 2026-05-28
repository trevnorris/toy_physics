ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 154 — EXACT CO-EVOLVING CORE-MOUTH MAP"];

Clear[g, r, dg];
$Assumptions = Element[{g, r, dg}, Reals];

rFun = (g - r)^2/(1 + r^2);
Print["R(g) = ", fmt[FullSimplify[rFun]]];

gStar = r - Sqrt[1 + r^2]/2;
expectZero["R(g_star) - 1/4", (rFun /. g -> gStar) - 1/4];

rShiftSeries = Normal[Series[rFun /. g -> gStar + dg, {dg, 0, 2}]];
rShiftExpected = 1/4 - dg/Sqrt[1 + r^2] + dg^2/(1 + r^2);
expectZero["exact shifted R formula", rShiftSeries - rShiftExpected];

banner["Linearized slope identity"];

Clear[sigma0, dSigma0, sStar, dS, rStar, dR];
$Assumptions = Element[{sigma0, dSigma0, sStar, dS, rStar, dR}, Reals];

piExpr = (sigma0 + dSigma0) * (1 - (rStar + dR) * (sStar + dS));
Clear[epsLin];
piExprEps = piExpr /. {dSigma0 -> epsLin*dSigma0, dR -> epsLin*dR, dS -> epsLin*dS};
piLin = Expand[Normal[Series[piExprEps, {epsLin, 0, 1}]] /. epsLin -> 1];
pi0 = sigma0*(1 - rStar*sStar);
dPi = Expand[piLin - pi0];
dPiExpected = Expand[(1 - rStar*sStar)*dSigma0 - sigma0*(rStar*dS + sStar*dR)];
expectZero["dPi identity", dPi - dPiExpected];

Print[""];
Print["Carry-forward formulas:"];
Print["  R(g) = (g-r)^2/(1+r^2)"];
Print["  g = g_*  <=>  R = 1/4 on the lower branch"];
Print["  delta R = -delta g/sqrt(1+r^2) + delta g^2/(1+r^2)"];
Print["  delta Pi = (1-R_*S_*) delta Sigma0 - Sigma0_* (R_* delta S + S_* delta R)"];

Print[""];
Print["Stage 154 Mathematica audit passed."];

Exit[0];
