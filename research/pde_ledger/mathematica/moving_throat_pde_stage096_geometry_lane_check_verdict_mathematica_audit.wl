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
  res = FullSimplify[Together[Expand[expr]]];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 096 — GEOMETRY-LANE CHECK VERDICT"];

lapS2[expr_, th_, ph_] := FullSimplify[
  (1/Sin[th]) D[Sin[th] D[expr, th], th] + D[expr, {ph, 2}]/Sin[th]^2
];

dOmegaIntegral[expr_, th_, ph_] := FullSimplify[
  Integrate[expr Sin[th], {ph, 0, 2 Pi}, {th, 0, Pi}]
];

Clear[theta, phi, omega, omegaQ];
$Assumptions = Element[{theta, phi, omega, omegaQ}, Reals] && 0 < theta < Pi && omegaQ > 0;

y00 = 1/(2 Sqrt[Pi]);
y20 = Sqrt[5/(16 Pi)] (3 Cos[theta]^2 - 1);
y21c = Sqrt[15/(4 Pi)] Sin[theta] Cos[theta] Cos[phi];
y21s = Sqrt[15/(4 Pi)] Sin[theta] Cos[theta] Sin[phi];
y22c = Sqrt[15/(16 Pi)] Sin[theta]^2 Cos[2 phi];
y22s = Sqrt[15/(16 Pi)] Sin[theta]^2 Sin[2 phi];

banner["SECTION I — l=0 / l=2 ORTHOGONALITY AND LAPLACE EIGENVALUE"];
Do[
  name = pair[[1]];
  y2m = pair[[2]];
  expectZero["<Y00|Y" <> name <> ">", dOmegaIntegral[y00*y2m, theta, phi]];
  expectZero["(-Delta)Y" <> name <> " - 6 Y" <> name, -lapS2[y2m, theta, phi] - 6 y2m];
  expectZero["<Y00|(-Delta)Y" <> name <> ">", dOmegaIntegral[y00*(-lapS2[y2m, theta, phi]), theta, phi]];
  ,
  {pair, {{"20", y20}, {"21c", y21c}, {"21s", y21s}, {"22c", y22c}, {"22s", y22s}}}
];

eps2 = 0;
eps4 = 0;
cPole = FullSimplify[(1 + eps4)/(4*(1 + eps2)^2)];
cGeom = FullSimplify[1 - cPole];
rhoAlpha = FullSimplify[1/cGeom];
zetaReq = FullSimplify[cPole/cGeom];
yhatCons = FullSimplify[cGeom + cPole/(1 - omega^2/omegaQ^2)];
yhatExpected = FullSimplify[3/4 + (1/4)/(1 - omega^2/omegaQ^2)];

banner["SECTION II — ISOTROPIC CONSERVATIVE CARRIER"];

Print["c_pole = ", fmt[cPole]];
Print["c_geom = ", fmt[cGeom]];
Print["Yhat_Q^cons(omega) = ", fmt[yhatCons]];
Print["rho_alpha = ", fmt[rhoAlpha]];
Print["zeta_req = ", fmt[zetaReq]];

Print["eps2 = ", fmt[eps2]];
Print["eps4 = ", fmt[eps4]];
expectZero["c_pole - 1/4", cPole - 1/4];
expectZero["c_geom - 3/4", cGeom - 3/4];
expectZero["Yhat_Q^cons - [3/4 + (1/4)/(1 - omega^2/Omega_Q^2)]", yhatCons - yhatExpected];
expectZero["rho_alpha - 4/3", rhoAlpha - 4/3];
expectZero["zeta_req - 1/3", zetaReq - 1/3];

Print[""];
Print["Stage 096 Mathematica audit passed."];

Print[""];
Print["HYPOTHESIS CARRIED"];
Print["These results are conditional on the Part III minimal isotropic module"];
Print["and the grouped real P_2 carrier. The card is a derivation ledger entry,"];
Print["not an unconditional actual-branch theorem."];

Exit[0];
