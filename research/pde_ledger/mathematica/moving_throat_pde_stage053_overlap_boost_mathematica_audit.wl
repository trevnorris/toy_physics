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

banner["STAGE 036 — EXACT OVERLAP-BOOST WINDOW"];

Clear[s, ell, alpha];
$Assumptions = Element[{s, ell, alpha}, Reals] && ell > 0 && alpha > 0;

chi0 = Sqrt[2/ell] Sin[Pi s/(2 ell)];
iW = FullSimplify[Integrate[chi0, {s, 0, ell}], Assumptions -> $Assumptions];
chi0Max = FullSimplify[Sqrt[2/ell], Assumptions -> $Assumptions];
omegaMax = FullSimplify[ell chi0Max/iW, Assumptions -> $Assumptions];
aIMax = FullSimplify[omegaMax^2, Assumptions -> $Assumptions];

Print["chi0(s) = ", fmt[chi0]];
Print["I_W = ", fmt[iW]];
Print["max chi0 = ", fmt[chi0Max]];
Print["Omega_max = ", fmt[omegaMax]];
Print["A_I,max = ", fmt[aIMax]];
expectZero["Omega_max - Pi/2", omegaMax - Pi/2];
expectZero["A_I,max - Pi^2/4", aIMax - Pi^2/4];

sigmaAlpha = FullSimplify[alpha Exp[alpha s/ell]/(Exp[alpha] - 1), Assumptions -> $Assumptions];
sigmaTotal = FullSimplify[Integrate[sigmaAlpha, {s, 0, ell}], Assumptions -> $Assumptions];
iAlpha = FullSimplify[Integrate[sigmaAlpha chi0, {s, 0, ell}], Assumptions -> $Assumptions];
omegaAlpha = FullSimplify[iAlpha/iW, Assumptions -> $Assumptions];
omegaAlphaSimple = FullSimplify[
  Pi alpha (2 alpha Exp[alpha] + Pi)/((4 alpha^2 + Pi^2) (Exp[alpha] - 1)),
  Assumptions -> $Assumptions
];

Print["sigma_alpha(s) = ", fmt[sigmaAlpha]];
Print["Integral sigma_alpha ds = ", fmt[sigmaTotal]];
Print["I_alpha = ", fmt[iAlpha]];
Print["Omega_alpha = ", fmt[omegaAlpha]];
expectZero["same total source strength", sigmaTotal - ell];
expectZero["Omega_alpha closed form", omegaAlpha - omegaAlphaSimple];

omega0 = FullSimplify[Limit[omegaAlpha, alpha -> 0], Assumptions -> ell > 0];
omegaInf = FullSimplify[Limit[omegaAlpha, alpha -> Infinity], Assumptions -> ell > 0];
seriesSmall = FullSimplify[Normal[Series[omegaAlpha, {alpha, 0, 2}]], Assumptions -> ell > 0];
linearCoeff = FullSimplify[2/Pi - 1/2];

Print["Omega_alpha(alpha->0) = ", fmt[omega0]];
Print["Omega_alpha(alpha->+infty) = ", fmt[omegaInf]];
Print["small-alpha series = ", fmt[seriesSmall]];
Print["linear coefficient = ", fmt[linearCoeff]];
expectZero["uniform branch value", omega0 - 1];
expectZero["antinode concentration limit", omegaInf - Pi/2];
expectZero["linear coefficient - (4-Pi)/(2Pi)", linearCoeff - (4 - Pi)/(2 Pi)];

Clear[zetaReq];
$Assumptions = Element[{zetaReq}, Reals] && zetaReq > 0;
criterion = FullSimplify[aIMax - zetaReq, Assumptions -> $Assumptions];
Print["A_I,max - zeta_req = ", fmt[criterion]];

Print[""];
Print["Stage 053 Mathematica audit passed."];

Exit[0];
