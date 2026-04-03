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

banner["STAGE 094 — MIXED SIDE-CHANNEL POLE"];

Clear[z, sigma, kappa, gamma];
$Assumptions =
  Element[{z, sigma, kappa, gamma}, Reals] &&
  sigma != -3;

lambdaOut = -3 + z^2/3 + z^4/9 + I*z^5/9;
lambdaMix = Expand[Normal[Series[lambdaOut - sigma/(1 - kappa*z^2 - I*gamma*z^5), {z, 0, 5}]]];

l0 = FullSimplify[Coefficient[lambdaMix, z, 0], Assumptions -> $Assumptions];
l2 = FullSimplify[Coefficient[lambdaMix, z, 2], Assumptions -> $Assumptions];
l4 = FullSimplify[Coefficient[lambdaMix, z, 4], Assumptions -> $Assumptions];
l5 = FullSimplify[Coefficient[lambdaMix, z, 5]/I, Assumptions -> $Assumptions];

Print["Lambda_mix(z) = ", fmt[lambdaMix]];
Print["L0 = ", fmt[l0]];
Print["L2 = ", fmt[l2]];
Print["L4 = ", fmt[l4]];
Print["L5 = ", fmt[l5]];

kappaMatch = FullSimplify[kappa /. First[Solve[-l2/l0 == 1/9, kappa, Reals]], Assumptions -> $Assumptions];
sigmaMatch = FullSimplify[
  sigma /. First[Solve[((l2^2/l0^2 - l4/l0) /. kappa -> kappaMatch) == 4/81, sigma, Reals]],
  Assumptions -> $Assumptions
];
chiMix = FullSimplify[(-l5/l0)/(1/27), Assumptions -> $Assumptions];
chiMixLinear = Expand[Normal[Series[chiMix, {sigma, 0, 1}]]];

Print["kappa from z^2 matching = ", fmt[kappaMatch]];
Print["sigma from z^4 matching = ", fmt[sigmaMatch]];
Print["chi_Q^mix = ", fmt[chiMix]];
Print["chi_Q^mix linearized = ", fmt[chiMixLinear]];

expectZero["kappa_match + 1/9", kappaMatch + 1/9];
expectZero["sigma_match", sigmaMatch];
expectZero["chi_Q^mix - 3 (1 - 9 sigma gamma)/(3 + sigma)", chiMix - 3*(1 - 9*sigma*gamma)/(3 + sigma)];
expectZero["chi_Q^mix linearized - (1 - sigma (1/3 + 9 gamma))", chiMixLinear - (1 - sigma*(1/3 + 9*gamma))];

Print[""];
Print["Stage 094 Mathematica audit passed."];

Exit[0];
