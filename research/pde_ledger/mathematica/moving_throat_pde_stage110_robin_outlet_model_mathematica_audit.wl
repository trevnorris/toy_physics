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

banner["STAGE 110 — EXPLICIT ISOTROPIC ROBIN OUTLET MODEL"];

Clear[z, rho];
$Assumptions = Element[{z, rho}, Reals] && rho != 3;

lambdaOut = -3 + z^2/3 + z^4/9 + I*z^5/9;
lambdaR = lambdaOut + rho;
yR = FullSimplify[(-3 + rho)/lambdaR, Assumptions -> $Assumptions];
yRSeries = Expand[Normal[Series[yR, {z, 0, 5}]]];

c2 = FullSimplify[Coefficient[yRSeries, z, 2], Assumptions -> $Assumptions];
c4 = FullSimplify[Coefficient[yRSeries, z, 4], Assumptions -> $Assumptions];
c5 = FullSimplify[Coefficient[yRSeries, z, 5]/I, Assumptions -> $Assumptions];
chiR = FullSimplify[c5/(1/27), Assumptions -> $Assumptions];
chiRLinear = Expand[Normal[Series[chiR, {rho, 0, 2}]]];

Print["Y_R(z) = ", fmt[yRSeries]];
Print["c2 = ", fmt[c2]];
Print["c4 = ", fmt[c4]];
Print["c5 = ", fmt[c5]];
Print["chi_Q^R = ", fmt[chiR]];
Print["chi_Q^R linearized = ", fmt[chiRLinear]];

expectZero["c2 - 1/(9 - 3 rho)", c2 - 1/(9 - 3*rho)];
expectZero["c4 - (4 - rho)/(9 (3 - rho)^2)", c4 - (4 - rho)/(9*(3 - rho)^2)];
expectZero["c5 - 1/(27 - 9 rho)", c5 - 1/(27 - 9*rho)];
expectZero["chi_Q^R - 3/(3 - rho)", chiR - 3/(3 - rho)];
expectZero["chi_Q^R linearized - (1 + rho/3 + rho^2/9)", chiRLinear - (1 + rho/3 + rho^2/9)];

Print[""];
Print["Stage 110 Mathematica audit passed."];

Exit[0];
