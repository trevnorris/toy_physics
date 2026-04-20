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

banner["STAGE 120 — EXPLICIT CORE-TO-MOUTH GAIN MAP"];

Clear[lM, thetaSigma, kS, kQ, lam, gS, gQ];
$Assumptions =
  Element[{lM, thetaSigma, kS, kQ, lam, gS, gQ}, Reals] &&
  lM > 0 && thetaSigma > 0 && kS > 0 && kQ > 0 && gS > 0 && gQ > 0;

rhoC = FullSimplify[gS^2/kS, Assumptions -> $Assumptions];
sigmaC = FullSimplify[(kS*gQ - lam*gS)^2/(kS*(kS*kQ + lam^2)), Assumptions -> $Assumptions];

mS = FullSimplify[lM*rhoC/thetaSigma, Assumptions -> $Assumptions];
mQ = FullSimplify[-lM*sigmaC/thetaSigma, Assumptions -> $Assumptions];

Print["rho_c = ", fmt[rhoC]];
Print["sigma_c = ", fmt[sigmaC]];
Print["M_s = ", fmt[mS]];
Print["M_q = ", fmt[mQ]];

sigmaCRc = FullSimplify[(kS*gQ - lam*gS)^2/(kS^2*kQ*(1 + lam^2/(kS*kQ))), Assumptions -> $Assumptions];
Print["sigma_c (r_c form) = ", fmt[sigmaCRc]];
expectZero["sigma_c equivalence with r_c form", sigmaC - sigmaCRc];

Print[""];
Print["Final explicit gain map verified."];

Exit[0];
