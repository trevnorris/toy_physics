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

banner["STAGE 035 — DIMENSIONLESS NORMALIZATION LOCUS"];

Clear[A, beta0, NQ, xi, delta, Chi, OmegaU, Delta0];
$Assumptions =
  Element[{A, beta0, NQ, xi, delta, Chi, OmegaU, Delta0}, Reals] &&
  A > 0 && beta0 > 0 && NQ > 0 && delta > 0 && 0 <= xi < 1 &&
  OmegaU > 0 && Delta0 > 0;

kappa0Sq = 8/Pi^2;
kappa1Sq = 16/(9*Pi^2);
eta = FullSimplify[kappa1Sq/kappa0Sq, Assumptions -> $Assumptions];
x = FullSimplify[A*xi, Assumptions -> $Assumptions];
deltaKSub = FullSimplify[A*delta, Assumptions -> $Assumptions];

alphaX = FullSimplify[
  x*(x + deltaKSub)/(kappa0Sq*(x + deltaKSub) + kappa1Sq*x),
  Assumptions -> $Assumptions
];
nX = FullSimplify[
  beta0*(kappa0Sq*(x + deltaKSub) + kappa1Sq*x)^4/
    (kappa0Sq*(A - x)*(kappa0Sq*(x + deltaKSub)^2 + kappa1Sq*x^2)^2),
  Assumptions -> $Assumptions
];

f = FullSimplify[nX/(beta0*kappa0Sq/A), Assumptions -> $Assumptions];
fTarget = FullSimplify[
  (9*delta + 11*xi)^4/(81*(1 - xi)*(9*delta^2 + 18*delta*xi + 11*xi^2)^2),
  Assumptions -> $Assumptions
];
rTarget = FullSimplify[NQ*A/(beta0*kappa0Sq), Assumptions -> $Assumptions];
rTargetClosed = FullSimplify[Pi^2*A*NQ/(8*beta0), Assumptions -> $Assumptions];

Print["eta = ", fmt[eta]];
Print["F(xi,delta) = ", fmt[f]];
Print["R_target = ", fmt[rTarget]];
expectZero["F - closed D/N form", f - fTarget];
expectZero["R_target - Pi^2 A NQ/(8 beta0)", rTarget - rTargetClosed];

dF = FullSimplify[D[f, xi], Assumptions -> $Assumptions];
dFTarget = FullSimplify[
  (9*delta + 11*xi)^3*(81*delta^3 + 72*delta^2 + 189*delta^2*xi + 297*delta*xi^2 + 121*xi^3)/
    (81*(1 - xi)^2*(9*delta^2 + 18*delta*xi + 11*xi^2)^3),
  Assumptions -> $Assumptions
];
Print["dF/dxi = ", fmt[dF]];
expectZero["dF/dxi - manifestly positive form", dF - dFTarget];
expectZero["F(0,delta) - 1", (fTarget /. xi -> 0) - 1];

softConst = Block[
  {$Assumptions = Element[{delta}, Reals] && delta > 0},
  FullSimplify[Limit[(1 - xi)*f, xi -> 1, Direction -> "FromBelow"], Assumptions -> $Assumptions]
];
softConstTarget = FullSimplify[
  (9*delta + 11)^4/(81*(9*delta^2 + 18*delta + 11)^2),
  Assumptions -> delta > 0
];
Print["softening constant C(delta) = ", fmt[softConst]];
expectZero["softening constant - closed form", softConst - softConstTarget];

epsSoft = Symbol["epsSoft"];
softScaledSeries = FullSimplify[
  Normal[Series[(epsSoft*f /. xi -> 1 - epsSoft), {epsSoft, 0, 0}]],
  Assumptions -> delta > 0 && epsSoft > 0
];
expectZero["near-softening asymptotic coefficient", softScaledSeries - softConstTarget];

alphaReq = FullSimplify[alphaX, Assumptions -> $Assumptions];
alphaReqTarget = FullSimplify[
  9*Pi^2*A*xi*(xi + delta)/(8*(9*delta + 11*xi)),
  Assumptions -> $Assumptions
];
alphaCrit = Block[
  {$Assumptions = Element[{A, delta}, Reals] && A > 0 && delta > 0},
  FullSimplify[Limit[alphaReq, xi -> 1, Direction -> "FromBelow"], Assumptions -> $Assumptions]
];
alphaCritTarget = FullSimplify[9*Pi^2*A*(1 + delta)/(8*(11 + 9*delta)), Assumptions -> A > 0 && delta > 0];

Print["alpha_req(xi,delta) = ", fmt[alphaReq]];
Print["alpha_crit = ", fmt[alphaCrit]];
expectZero["alpha_req - closed D/N form", alphaReq - alphaReqTarget];
expectZero["alpha_crit - closed form", alphaCrit - alphaCritTarget];

alphaMix = FullSimplify[Chi^2/(OmegaU^2*Delta0), Assumptions -> $Assumptions];
gBReqSqOverVarpi2 = FullSimplify[alphaReqTarget - alphaMix, Assumptions -> $Assumptions];
Print["g_B,req^2 / varpi^2 = ", fmt[gBReqSqOverVarpi2]];

fSeries = FullSimplify[Normal[Series[f, {xi, 0, 2}]], Assumptions -> delta > 0];
fSeriesTarget = FullSimplify[
  1 + (1 + 8/(9*delta))*xi + (1 + 8/(9*delta) - 28/(27*delta^2))*xi^2,
  Assumptions -> delta > 0
];
alphaSeries = FullSimplify[Normal[Series[alphaReq, {xi, 0, 2}]], Assumptions -> A > 0 && delta > 0];
alphaSeriesTarget = FullSimplify[Pi^2*A*xi/8 - Pi^2*A*xi^2/(36*delta), Assumptions -> A > 0 && delta > 0];

Print["F(xi,delta) near xi=0 = ", fmt[fSeries]];
Print["alpha_req near xi=0 = ", fmt[alphaSeries]];
expectZero["near-onset F series through O(xi^2)", fSeries - fSeriesTarget];
expectZero["alpha_req near-onset series through O(xi^2)", alphaSeries - alphaSeriesTarget];

Print[""];
Print["Stage 035 Mathematica audit passed."];

Exit[0];
