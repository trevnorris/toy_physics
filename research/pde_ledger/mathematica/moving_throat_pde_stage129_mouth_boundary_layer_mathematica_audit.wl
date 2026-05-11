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

banner["STAGE 112 — EXPLICIT GNLS + LOCALIZED-MAXWELL MOUTH BOUNDARY LAYER"];

Clear[z, lM, piM, thetaSigma, v1, mobility, sigmaStar];
$Assumptions =
  Element[{z, lM, piM, thetaSigma, v1, mobility, sigmaStar}, Reals] &&
  lM > 0 && piM > 0 && thetaSigma > 0 && mobility > 0 && sigmaStar > 0;

sigma = piM*Exp[-piM*z/lM]/(lM*(1 - Exp[-piM]));
jSigma = -mobility*(thetaSigma*D[sigma, z] + v1*sigma);

Print["sigma_Pi(z) = ", fmt[sigma]];
Print["Normalization = ", fmt[FullSimplify[Integrate[sigma, {z, 0, lM}], Assumptions -> $Assumptions]]];
expectZero["profile normalization", Integrate[sigma, {z, 0, lM}] - 1];

jSub = FullSimplify[jSigma /. v1 -> piM*thetaSigma/lM, Assumptions -> $Assumptions];
Print["Zero-flux current J_sigma = ", fmt[jSub]];
expectZero["zero-flux current", jSub];

residual = FullSimplify[(thetaSigma*D[sigma, z] + v1*sigma) /. v1 -> piM*thetaSigma/lM, Assumptions -> $Assumptions];
Print["Stationary zero-flux ODE residual = ", fmt[residual]];
expectZero["boundary-layer ODE residual", residual];

Print[""];
Print["Derived electrochemical bias:"];
Print["Pi_m = V1*L/Theta_sigma"];

Print[""];
Print["Stage 129 Mathematica audit passed."];

Exit[0];
