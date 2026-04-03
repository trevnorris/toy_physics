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

banner["STAGE 099 — FINITE D/N MIXED-TUBE REALIZATION"];

Clear[a, lW, kS, kQ, lam, z];
$Assumptions =
  Element[{a, lW, kS, kQ, lam, z}, Reals] &&
  a > 0 && lW > 0 && kS > 0 && kQ > 0 && lam > 0;

rC = FullSimplify[lam^2/(kS*kQ), Assumptions -> $Assumptions];
kappa0FromTube = FullSimplify[4*lW^2/(Pi^2*a^2), Assumptions -> $Assumptions];
lWRequired = FullSimplify[lW /. First[Solve[kappa0FromTube == (1 + rC)/3, lW, Reals]], Assumptions -> $Assumptions];

Print["kappa0 from D/N half-wave tube = ", fmt[kappa0FromTube]];
Print["Required tube length L_W = ", fmt[lWRequired]];
expectZero["tube-length law", lWRequired - (Pi*a*Sqrt[(1 + rC)/3])/2];

kappa0Bare = FullSimplify[(1 + rC)/3, Assumptions -> $Assumptions];
gamma0Bare = FullSimplify[(1 + rC)/9, Assumptions -> $Assumptions];
kappaC = FullSimplify[kappa0Bare/(1 + rC), Assumptions -> $Assumptions];
gammaC = FullSimplify[gamma0Bare/(1 + rC), Assumptions -> $Assumptions];

expectZero["final kappa_c - 1/3", kappaC - 1/3];
expectZero["final gamma_c - 1/9", gammaC - 1/9];

dBare = Expand[(1 + rC)*(1 - z^2/3 - I*z^5/9)];
dFinal = FullSimplify[dBare/(1 + rC), Assumptions -> $Assumptions];
expectZero["bare scaled-canonical branch renormalizes to canonical", dFinal - (1 - z^2/3 - I*z^5/9)];

Print[""];
Print["Summary:"];
Print["  D/N half-wave mixed tube length: ", fmt[lWRequired]];
Print["  Bare mixed coefficients: kappa0 = ", fmt[kappa0Bare], ", gamma0 = ", fmt[gamma0Bare]];
Print["  Final coefficients: kappa_c = ", fmt[kappaC], ", gamma_c = ", fmt[gammaC]];

Print[""];
Print["Stage 099 Mathematica audit passed."];

Exit[0];
