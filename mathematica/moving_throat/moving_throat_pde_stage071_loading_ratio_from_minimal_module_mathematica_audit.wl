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

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 071 — LOADING-RATIO EXTRACTION FROM THE MINIMAL ISOTROPIC MODULE"];

Clear[rhoAlpha, omega, omegaQ, alphaReq, alphaMix, c0, c1, cMix];
$Assumptions =
  Element[{rhoAlpha, omega, omegaQ, alphaReq, alphaMix, c0, c1, cMix}, Reals] &&
  rhoAlpha > 0 && omegaQ > 0 && alphaReq > 0 && alphaMix > 0 && c0 > 0 && c1 > 0 && cMix > 0;

yLoading = FullSimplify[alphaMix/alphaReq + ((alphaReq - alphaMix)/alphaReq)/(1 - omega^2/omegaQ^2), Assumptions -> $Assumptions];
yRho = FullSimplify[1/rhoAlpha + ((rhoAlpha - 1)/rhoAlpha)/(1 - omega^2/omegaQ^2), Assumptions -> $Assumptions];

Print["Y_loading(omega) = ", fmt[yLoading]];
Print["Y_rho(omega) = ", fmt[yRho]];
expectZero["loading form - rho form", (yLoading /. alphaReq -> rhoAlpha*alphaMix) - yRho];

c0FromRho = FullSimplify[1/rhoAlpha, Assumptions -> $Assumptions];
c1FromRho = FullSimplify[(rhoAlpha - 1)/rhoAlpha, Assumptions -> $Assumptions];
expectZero["contact-plus-pole reconstruction", yRho - (c0FromRho + c1FromRho/(1 - omega^2/omegaQ^2))];
expectZero["c0 + c1 - 1", c0FromRho + c1FromRho - 1];

rhoFromC0 = 1/c0;
rhoFromC1 = 1/(1 - c1);
zetaFromC = c1/c0;
expectZero["rho(c0(rho)) - rho", (rhoFromC0 /. c0 -> c0FromRho) - rhoAlpha];
expectZero["rho(c1(rho)) - rho", (rhoFromC1 /. c1 -> c1FromRho) - rhoAlpha];
expectZero["zeta(c(rho)) - (rho-1)", (zetaFromC /. {c0 -> c0FromRho, c1 -> c1FromRho}) - (rhoAlpha - 1)];

rhoMin = FullSimplify[rhoFromC0 /. c0 -> 3/4];
zetaMin = FullSimplify[zetaFromC /. {c0 -> 3/4, c1 -> 1/4}];
piMin = FullSimplify[(4/3)*cMix];

Print["rho_alpha(minimal isotropic module) = ", fmt[rhoMin]];
Print["zeta_req(minimal isotropic module) = ", fmt[zetaMin]];
expectZero["rho_min - 4/3", rhoMin - 4/3];
expectZero["zeta_min - 1/3", zetaMin - 1/3];
expectZero["Pi_tr/C_mix - 4/3", piMin/cMix - 4/3];
expectTrue["C_mix < Pi_tr < 2 C_mix", cMix < piMin < 2*cMix];

Print[""];
Print["Stage 071 Mathematica audit passed."];

Exit[0];
