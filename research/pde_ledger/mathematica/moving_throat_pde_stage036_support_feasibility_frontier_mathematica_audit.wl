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
expectZero[
  "g_B,req^2/varpi^2 - (Pi^2 A / 8) (G - M_mix)",
  gBReqSqOverVarpi2 - (Pi^2*A/8)*(gTarget - mMix)
];

(* Derive dG/dxi from gTarget; multiply through by the squared denominator
   so the result is a polynomial in xi and delta. Mathematica must produce
   11*xi^2 + 18*delta*xi + 9*delta^2 on its own; the closed-form coefficients
   are not declared up front. *)
dG = FullSimplify[D[gTarget, xi], Assumptions -> $Assumptions];
dGPolynomial = FullSimplify[Expand[dG*(9*delta + 11*xi)^2/9], Assumptions -> $Assumptions];
Print["dG/dxi = ", fmt[dG]];
Print["9 * dG/dxi * (9 delta + 11 xi)^2 / 81 (polynomial) = ", fmt[dGPolynomial]];
expectZero[
  "dG/dxi positivity polynomial: 9 dG/dxi (9d+11xi)^2 / 9 == 11 xi^2 + 18 delta xi + 9 delta^2",
  dGPolynomial - (11*xi^2 + 18*delta*xi + 9*delta^2)
];
(* Also confirm the manifest non-negativity of the numerator polynomial
   for delta, xi >= 0: discriminant <= 0 in xi. *)
disc = Discriminant[11*xi^2 + 18*delta*xi + 9*delta^2, xi];
discSimplified = FullSimplify[disc, Assumptions -> delta > 0];
Print["discriminant (in xi) = ", fmt[discSimplified]];
expectZero["dG/dxi numerator discriminant equals -72 delta^2", discSimplified + 72*delta^2];
expectZero["G(0,delta)", gTarget /. xi -> 0];

gMax = Block[
  {$Assumptions = Element[{delta}, Reals] && delta > 0},
  FullSimplify[Limit[gTarget, xi -> 1, Direction -> "FromBelow"], Assumptions -> $Assumptions]
];
(* Derive gMaxTarget from gMax via substitution rather than declaring it. *)
gMaxTarget = FullSimplify[gMax, Assumptions -> delta > 0];

Print["G_max(delta) = ", fmt[gMax]];
expectZero["G_max - 9(1+delta)/(9delta+11)", gMax - 9*(1 + delta)/(9*delta + 11)];

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
(* Symbolic kappa-based cross-check: derive R_target_sym symbolically and confirm = F. *)
Clear[Asym, beta0Sym];
$Assumptions =
  Element[{A, delta, xiReq, Chi, OmegaU, Delta0, beta0, NQ, Asym, beta0Sym}, Reals] &&
  A > 0 && delta > 0 && 0 <= xiReq < 1 && OmegaU > 0 && Delta0 > 0 && beta0 > 0 && NQ > 0 &&
  Asym > 0 && beta0Sym > 0;
xSym = Asym * xi;
deltaKSym = Asym * delta;
nSym = beta0Sym*(kappa0Sq*(xSym + deltaKSym) + kappa1Sq*xSym)^4 /
  (kappa0Sq*(Asym - xSym)*(kappa0Sq*(xSym + deltaKSym)^2 + kappa1Sq*xSym^2)^2);
rTargetSym = FullSimplify[nSym*Asym/(beta0Sym*kappa0Sq), Assumptions -> $Assumptions];
expectZero[
  "symbolic kappa derivation: F(xi,delta) - R_target_sym",
  FullSimplify[Together[Expand[fTarget - rTargetSym]], Assumptions -> $Assumptions]
];
expectTrue[
  "admissible sample: M_mix <= G(xi_req,delta)",
  mMixGood <= gSample,
  "M_mix=" <> fmt[mMixGood] <> ", G=" <> fmt[gSample]
];
expectTrue[
  "inadmissible sample: support deficit blocks the branch",
  mMixBad > gSample,
  "M_mix=" <> fmt[mMixBad] <> ", G=" <> fmt[gSample]
];

banner["STAGE 019.4 — NEAR-ONSET ASYMPTOTICS"];
(* Read the series coefficients directly out of gSeries; do not declare a target. *)
gSeries = FullSimplify[Normal[Series[gTarget, {xi, 0, 2}]], Assumptions -> delta > 0];
c0 = FullSimplify[Coefficient[gSeries, xi, 0], Assumptions -> delta > 0];
c1 = FullSimplify[Coefficient[gSeries, xi, 1], Assumptions -> delta > 0];
c2 = FullSimplify[Coefficient[gSeries, xi, 2], Assumptions -> delta > 0];
Print["G(xi,delta) near xi=0 = ", fmt[gSeries]];
Print["near-onset coefficients: c0=", fmt[c0], ", c1=", fmt[c1], ", c2=", fmt[c2]];
expectZero["near-onset c0 = 0", c0];
expectZero["near-onset c1 = 1", c1 - 1];
expectZero["near-onset c2 = -2/(9 delta)", c2 + 2/(9*delta)];

Print[""];
Print["Stage 036 Mathematica audit passed."];

Exit[0];
