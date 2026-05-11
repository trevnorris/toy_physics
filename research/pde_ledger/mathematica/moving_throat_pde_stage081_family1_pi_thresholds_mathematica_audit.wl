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
  diff = Abs[N[value, 50] - N[target, 50]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS"];

Clear[zeta, piTr, cMix, epsBlk];
$Assumptions =
  Element[{zeta, piTr, cMix, epsBlk}, Reals] &&
  zeta >= 0 && piTr > 0 && cMix > 0 && epsBlk >= 0 && epsBlk < 1;

zetaExpr = FullSimplify[(piTr - cMix)/(cMix - epsBlk*(2*cMix - piTr)), Assumptions -> $Assumptions];
piOfZeta = FullSimplify[piTr /. First[Solve[zeta == zetaExpr, piTr]], Assumptions -> $Assumptions];
qq = FullSimplify[(1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta), Assumptions -> $Assumptions];
dqq = FullSimplify[D[qq, zeta], Assumptions -> $Assumptions];

Print["Pi_of_zeta = ", fmt[piOfZeta]];
Print["Q(zeta;eps_blk) = ", fmt[qq]];
expectZero["Q(0)-1", (qq /. zeta -> 0) - 1];
expectZero["Q(1)-2", (qq /. zeta -> 1) - 2];
expectZero["dQ/dzeta exact formula", dqq - (1 - epsBlk)/(1 - epsBlk*zeta)^2];

zetaSuffChi1 = ToExpression["2.46622291347846`20"];
zetaFailChi1 = ToExpression["2.46752913273870`20"];
zetaSuffJ1 = ToExpression["2.44257571477179`20"];
zetaFailJ1 = ToExpression["2.46752736855058`20"];
zetaMaxF1 = ToExpression["2.46752922945601`20"];

piSuffChiOverC = FullSimplify[qq /. zeta -> zetaSuffChi1];
piFailChiOverC = FullSimplify[qq /. zeta -> zetaFailChi1];
piSuffJOverC = FullSimplify[qq /. zeta -> zetaSuffJ1];
piFailJOverC = FullSimplify[qq /. zeta -> zetaFailJ1];
piMaxOverC = FullSimplify[qq /. zeta -> zetaMaxF1];

Print["Pi_suff^(chi)/C_mix = ", fmt[piSuffChiOverC]];
Print["Pi_fail^(chi)/C_mix = ", fmt[piFailChiOverC]];
Print["Pi_suff^(J)/C_mix = ", fmt[piSuffJOverC]];
Print["Pi_fail^(J)/C_mix = ", fmt[piFailJOverC]];
Print["Pi_max^(F1)/C_mix = ", fmt[piMaxOverC]];

expectApprox["Pi_suff^(chi)/C_mix at eps=0", N[piSuffChiOverC /. epsBlk -> 0, 40], ToExpression["3.4662229134784601214`25"], 10^-14];
expectApprox["Pi_fail^(chi)/C_mix at eps=0", N[piFailChiOverC /. epsBlk -> 0, 40], ToExpression["3.4675291327386998930`25"], 10^-14];
expectApprox["Pi_suff^(J)/C_mix at eps=0", N[piSuffJOverC /. epsBlk -> 0, 40], ToExpression["3.4425757147717899187`25"], 10^-14];
expectApprox["Pi_fail^(J)/C_mix at eps=0", N[piFailJOverC /. epsBlk -> 0, 40], ToExpression["3.4675273685505798582`25"], 10^-14];
expectApprox["Pi_max^(F1)/C_mix at eps=0", N[piMaxOverC /. epsBlk -> 0, 40], ToExpression["3.4675292294560100537`25"], 10^-14];

epsCeiling = N[1/zetaMaxF1, 40];
Print["Blocking ceiling eps_blk < ", fmt[epsCeiling]];
expectApprox["blocking ceiling numeric check", epsCeiling, ToExpression["0.40526368971137149977`25"], 10^-14];

Print[""];
Print["Stage 081 Mathematica audit passed."];

Exit[0];
