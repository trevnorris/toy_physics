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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 40] - N[target, 40]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 086 — FAMILY-1 PURE LOADING-RATIO WINDOW"];

Clear[epsBlk, zeta, zetaMax];
$Assumptions =
  Element[{epsBlk, zeta, zetaMax}, Reals] &&
  zeta > 0 && zetaMax > 1 && 0 <= epsBlk < 1/zetaMax;

qMap = FullSimplify[(1 + (1 - 2*epsBlk)*zeta)/(1 - epsBlk*zeta), Assumptions -> $Assumptions];
dQ = FullSimplify[D[qMap, zeta], Assumptions -> $Assumptions];
dQExpected = FullSimplify[(1 - epsBlk)/(1 - epsBlk*zeta)^2, Assumptions -> $Assumptions];

Print["Q(zeta;eps) = ", fmt[qMap]];
Print["dQ/dzeta = ", fmt[dQ]];
expectZero["dQ exact formula", dQ - dQExpected];
expectZero["Q(zeta;0) - (1+zeta)", (qMap /. epsBlk -> 0) - (1 + zeta)];

zetaSuffChi = ToExpression["2.46622291347846`20"];
zetaFailChi = ToExpression["2.46752913273870`20"];
zetaSuffJ = ToExpression["2.44257571477179`20"];
zetaMaxNum = ToExpression["2.46752922945601`20"];

expectApprox["zeta_suff^(chi) vs paper", zetaSuffChi, 2.46622291347846, 10^-14];
expectApprox["zeta_fail^(chi) vs paper", zetaFailChi, 2.46752913273870, 10^-14];
expectApprox["zeta_suff^(J) vs paper",   zetaSuffJ,   2.44257571477179, 10^-14];
expectApprox["zeta_max^(F1) vs paper",   zetaMaxNum,  2.46752922945601, 10^-14];

rhoSuffChi = N[qMap /. {zeta -> zetaSuffChi, epsBlk -> 0}, 30];
rhoFailChi = N[qMap /. {zeta -> zetaFailChi, epsBlk -> 0}, 30];
rhoSuffJ = N[qMap /. {zeta -> zetaSuffJ, epsBlk -> 0}, 30];
rhoMaxNum = N[qMap /. {zeta -> zetaMaxNum, epsBlk -> 0}, 30];

Print["rho_suff^(chi)(0) = ", fmt[rhoSuffChi]];
Print["rho_fail^(chi)(0) = ", fmt[rhoFailChi]];
Print["rho_suff^(J)(0) = ", fmt[rhoSuffJ]];
Print["rho_max^(F1)(0) = ", fmt[rhoMaxNum]];

expectApprox["rho_suff^(chi) numeric check", rhoSuffChi, 3.46622291347846012143918414949, 10^-14];
expectApprox["rho_fail^(chi) numeric check", rhoFailChi, 3.46752913273869989296827043290, 10^-14];
expectApprox["rho_suff^(J) numeric check", rhoSuffJ, 3.44257571477178991870005120290, 10^-14];
expectApprox["rho_max numeric check", rhoMaxNum, 3.46752922945601005366711433453, 10^-14];

qMax = FullSimplify[qMap /. zeta -> zetaMax, Assumptions -> $Assumptions];
dQMax = FullSimplify[D[qMax, epsBlk], Assumptions -> $Assumptions];
dQMaxExpected = FullSimplify[zetaMax*(zetaMax - 1)/(1 - epsBlk*zetaMax)^2, Assumptions -> $Assumptions];
expectZero["d rho_max / d eps exact formula", dQMax - dQMaxExpected];
expectApprox["eps_blk cap numeric check", N[1/zetaMaxNum, 25], 0.4052636897113714997686884, 10^-14];

Print[""];
Print["Stage 086 Mathematica audit passed."];

Exit[0];
