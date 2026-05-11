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

banner["STAGE 097 — CONCRETE TWO-CHANNEL CORE OUTLET MODEL"];

Clear[kS, kQ, lam, gS, gQ, kappa0, gamma0, z, dSym];
$Assumptions =
  Element[{kS, kQ, lam, gS, gQ, kappa0, gamma0, z, dSym}, Reals] &&
  kS > 0 && kQ > 0 && kappa0 > 0 && gamma0 > 0;

m = {{kS, lam}, {lam, -kQ*dSym}};
c = {gS, gQ};

deltaD = FullSimplify[Apart[c.Inverse[m].c, dSym], Assumptions -> $Assumptions];
Print["delta_Lambda(D) = ", fmt[deltaD]];

rhoC = FullSimplify[gS^2/kS, Assumptions -> $Assumptions];
rC = FullSimplify[lam^2/(kS*kQ), Assumptions -> $Assumptions];
sigmaTilde = FullSimplify[(kS*gQ - lam*gS)^2/(kS^2*kQ), Assumptions -> $Assumptions];
sigmaC = FullSimplify[sigmaTilde/(1 + rC), Assumptions -> $Assumptions];
kappaC = FullSimplify[kappa0/(1 + rC), Assumptions -> $Assumptions];
gammaC = FullSimplify[gamma0/(1 + rC), Assumptions -> $Assumptions];

targetD = FullSimplify[rhoC - sigmaTilde/(dSym + rC), Assumptions -> $Assumptions];
expectZero["Schur form identity", deltaD - targetD];

dBare = 1 - kappa0*z^2 - I*gamma0*z^5;
deltaZ = FullSimplify[deltaD /. dSym -> dBare, Assumptions -> $Assumptions];
targetZ = FullSimplify[rhoC - sigmaC/(1 - kappaC*z^2 - I*gammaC*z^5), Assumptions -> $Assumptions];
expectZero["low-frequency normalized outlet identity", deltaZ - targetZ];

Print[""];
Print["Exact core-level identifications:"];
Print["rho_c   = ", fmt[rhoC]];
Print["sigma_c = ", fmt[sigmaC]];
Print["kappa_c = ", fmt[kappaC]];
Print["gamma_c = ", fmt[gammaC]];

Print[""];
Print["Stage 114 Mathematica audit passed."];

Exit[0];
