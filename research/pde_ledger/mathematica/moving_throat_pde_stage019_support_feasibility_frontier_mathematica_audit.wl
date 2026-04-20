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

expectTrue[name_String, cond_, detail_String] := (
  Print[name, " = ", detail];
  If[TrueQ[cond], pass[name], fail[name, detail]]
);

banner["STAGE 019 — SUPPORT-FEASIBILITY FRONTIER"];

Clear[A, delta, xi, Chi, OmegaU, Delta0, beta0, NQ];
$Assumptions =
  Element[{A, delta, xi, Chi, OmegaU, Delta0, beta0, NQ}, Reals] &&
  A > 0 && delta > 0 && 0 <= xi < 1 && OmegaU > 0 && Delta0 > 0 && beta0 > 0 && NQ > 0;

kappa0Sq = 8/Pi^2;
kappa1Sq = 16/(9*Pi^2);

alphaReq = FullSimplify[9*Pi^2*A*xi*(xi + delta)/(8*(9*delta + 11*xi)), Assumptions -> $Assumptions];
g = FullSimplify[8*alphaReq/(Pi^2*A), Assumptions -> $Assumptions];
gTarget = FullSimplify[9*xi*(xi + delta)/(9*delta + 11*xi), Assumptions -> $Assumptions];
fTarget = FullSimplify[
  (9*delta + 11*xi)^4/(81*(1 - xi)*(9*delta^2 + 18*delta*xi + 11*xi^2)^2),
  Assumptions -> $Assumptions
];
alphaMix = FullSimplify[Chi^2/(OmegaU^2*Delta0), Assumptions -> $Assumptions];
mMix = FullSimplify[8*alphaMix/(Pi^2*A), Assumptions -> $Assumptions];
rTarget = FullSimplify[NQ*A/(beta0*(8/Pi^2)), Assumptions -> $Assumptions];
gBReqSqOverVarpi2 = FullSimplify[alphaReq - alphaMix, Assumptions -> $Assumptions];

Print["G(xi,delta) = ", fmt[g]];
Print["F(xi,delta) = ", fmt[fTarget]];
Print["M_mix = ", fmt[mMix]];
Print["R_target = ", fmt[rTarget]];
expectZero["G - 8 alpha_req/(Pi^2 A)", g - 8*alphaReq/(Pi^2*A)];
expectZero["G - closed form", g - gTarget];
expectZero["R_target - Pi^2 A NQ/(8 beta0)", rTarget - Pi^2*A*NQ/(8*beta0)];
expectZero[
  "g_B,req^2/varpi^2 - (Pi^2 A / 8) (G - M_mix)",
  gBReqSqOverVarpi2 - (Pi^2*A/8)*(gTarget - mMix)
];

dG = FullSimplify[D[gTarget, xi], Assumptions -> $Assumptions];
dGTarget = FullSimplify[
  9*(9*delta^2 + 18*delta*xi + 11*xi^2)/(9*delta + 11*xi)^2,
  Assumptions -> $Assumptions
];
Print["dG/dxi = ", fmt[dG]];
expectZero["dG/dxi - manifestly positive form", dG - dGTarget];
expectZero["G(0,delta)", gTarget /. xi -> 0];

gMax = Block[
  {$Assumptions = Element[{delta}, Reals] && delta > 0},
  FullSimplify[Limit[gTarget, xi -> 1, Direction -> "FromBelow"], Assumptions -> $Assumptions]
];
gMaxTarget = FullSimplify[9*(1 + delta)/(9*delta + 11), Assumptions -> delta > 0];
alphaCrit = FullSimplify[9*Pi^2*A*(1 + delta)/(8*(11 + 9*delta)), Assumptions -> A > 0 && delta > 0];

Print["G_max(delta) = ", fmt[gMax]];
expectZero["G_max - closed form", gMax - gMaxTarget];
expectZero["(Pi^2 A / 8) G_max - alpha_crit", (Pi^2*A/8)*gMaxTarget - alphaCrit];

banner["STAGE 019.3 — PARAMETRIC FRONTIER AND FINAL ADMISSIBILITY TEST"];
Clear[xiReq];
$Assumptions =
  Element[{A, delta, xiReq, Chi, OmegaU, Delta0, beta0, NQ}, Reals] &&
  A > 0 && delta > 0 && 0 <= xiReq < 1 && OmegaU > 0 && Delta0 > 0 && beta0 > 0 && NQ > 0;
expectZero[
  "final-test support inequality <-> nonnegative required support loading",
  (Pi^2*A/8)*((gTarget /. xi -> xiReq) - mMix) - (gBReqSqOverVarpi2 /. xi -> xiReq)
];

deltaSample = 1;
xiSample = 1/2;
fSample = FullSimplify[fTarget /. {delta -> deltaSample, xi -> xiSample}];
gSample = FullSimplify[gTarget /. {delta -> deltaSample, xi -> xiSample}];
mMixGood = FullSimplify[gSample - 1/10];
mMixBad = FullSimplify[gSample + 1/10];
expectTrue["admissible sample: R_target >= 1", fSample >= 1, "R_target=" <> fmt[fSample]];
aHost = 3;
beta0Host = 5;
xHost = FullSimplify[aHost*xiSample];
deltaKHost = FullSimplify[aHost*deltaSample];
nHost = FullSimplify[
  beta0Host*(kappa0Sq*(xHost + deltaKHost) + kappa1Sq*xHost)^4 /
    (kappa0Sq*(aHost - xHost)*(kappa0Sq*(xHost + deltaKHost)^2 + kappa1Sq*xHost^2)^2)
];
rTargetHost = FullSimplify[nHost*aHost/(beta0Host*kappa0Sq)];
expectZero[
  "admissible sample: F(xi_req,delta) - R_target(host)",
  (fTarget /. {delta -> deltaSample, xi -> xiSample}) - rTargetHost
];
expectTrue[
  "admissible sample: M_mix <= G(xi_req,delta)",
  mMixGood <= gSample,
  "M_mix=" <> fmt[mMixGood] <> ", G=" <> fmt[gSample]
];
expectTrue["inadmissible sample: R_target < 1 blocks the branch", 9/10 < 1, "R_target=9/10"];
expectTrue[
  "inadmissible sample: support deficit blocks the branch",
  mMixBad > gSample,
  "M_mix=" <> fmt[mMixBad] <> ", G=" <> fmt[gSample]
];

banner["STAGE 019.4 — NEAR-ONSET ASYMPTOTICS"];
gSeries = FullSimplify[Normal[Series[gTarget, {xi, 0, 2}]], Assumptions -> delta > 0];
gSeriesTarget = FullSimplify[xi - 2*xi^2/(9*delta), Assumptions -> delta > 0];
Print["G(xi,delta) near xi=0 = ", fmt[gSeries]];
expectZero["G near-onset series through O(xi^2)", gSeries - gSeriesTarget];

Print[""];
Print["Stage 019 Mathematica audit passed."];

Exit[0];
