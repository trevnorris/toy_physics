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

banner["STAGE 147 — FIRST-ORDER RIGIDITY KERNEL"];

Clear[p, x, eps, gBar, sBar];
$Assumptions = Element[{p, x, eps, gBar, sBar}, Reals] && p > 0 && 0 <= x <= 1;

kap = Pi/2;
c = Cos[Pi*x/2];
kq = Cosh[kap*(1 - x)]/Cosh[kap];

gFormula = 2*p*(2*p*Exp[p] + Pi)/((4*p^2 + Pi^2)*(Exp[p] - 1));
sFormula = p*(kap*Tanh[kap] + p*(Exp[-p]*Sech[kap] - 1))/((1 - Exp[-p])*(kap^2 - p^2));

rF1 = Sqrt[(12*(37/20)^2)/Pi^2 - 1];
gMinus = FullSimplify[rF1 - Sqrt[1 + rF1^2]/2, Assumptions -> True];
pStar = p /. FindRoot[gFormula == gMinus, {p, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
gStar = N[gFormula /. p -> pStar, 40];
sStar = N[sFormula /. p -> pStar, 40];
gPrimeStar = N[D[gFormula, p] /. p -> pStar, 40];
sPrimeStar = N[D[sFormula, p] /. p -> pStar, 40];

sigmaStar = N[pStar/(1 - sStar/4), 40];
tStar = N[Sqrt[9*sigmaStar/20], 40];

aT = N[
  -(9/(40*tStar))*(1/(gPrimeStar*(1 - sStar/4)) + pStar*sPrimeStar/(4*gPrimeStar*(1 - sStar/4)^2)),
  30
];
bT = N[(9/(40*tStar))*pStar/(4*(1 - sStar/4)^2), 30];

banner["Stage 147 audit: first-order rigidity coefficients"];
Print["Pi_* = ", fmt[N[pStar, 30]]];
Print["Sigma_* = ", fmt[sigmaStar]];
Print["T_* = ", fmt[tStar]];
Print["A_T = ", fmt[aT]];
Print["B_T = ", fmt[bT]];
Print["|A_T|/B_T = ", fmt[N[Abs[aT]/bT, 20]]];

(* --- Audit assertions: numerical anchor against paper-quoted literals --- *)
aTPaper = -4.27263956256927`30;
bTPaper = 0.134875005736706`30;
ratioCrosscheck = 31.6785`20;
expectZero["A_T vs paper -4.27263956256927",
  If[Abs[aT - aTPaper] < 10^-12, 0, aT - aTPaper]];
expectZero["B_T vs paper 0.134875005736706",
  If[Abs[bT - bTPaper] < 10^-12, 0, bT - bTPaper]];
expectZero["|A_T|/B_T computed ratio cross-check (not a paper literal) vs 31.6785",
  If[Abs[Abs[aT]/bT - ratioCrosscheck] < 10^-3, 0, Abs[aT]/bT - ratioCrosscheck]];

(* --- Audit assertion: A_T from symbolic differentiation of T_m(p) (independent) --- *)
(* Build T_m as a symbolic function of p from sFormula and let D[] derive dT_m/dp. *)
(* Retuning identity: A_T = -(dT_m/dp)/(dg/dp) at pStar. Independent of the hand-  *)
(* written closed-form aT (wl:49-52). *)
tmOfP = Sqrt[(9/20)*(p/(1 - sFormula/4))];
dTmDp = D[tmOfP, p];
dgDp = D[gFormula, p];
aTAutodiff = N[-(dTmDp /. p -> pStar)/(dgDp /. p -> pStar), 30];
expectZero["A_T closed form vs autodiff of T_m(p)",
  If[Abs[aTAutodiff - aT] < 10^-20, 0, aTAutodiff - aT]];

dT = FullSimplify[eps*(aT*(gBar - gMinus) + bT*(sBar - sStar))];
Print["delta T_m = ", fmt[dT]];

wCenter = FullSimplify[aT*(c - gMinus) + bT*(kq - sStar)];
Print["Centered rigidity kernel W_*(x) = ", fmt[wCenter]];

(* --- Audit assertion: full symbolic x-independence of the centering offset --- *)
(* W_*(x) - (aT c(x) + bT Kq(x)) must be CONSTANT in x (the centering shift). Check the *)
(* symbolic derivative in x is identically zero -- holds for all x, not one sample. *)
wStar = aT*(c - gMinus) + bT*(kq - sStar);
expectZero["W_* centering offset is x-independent",
  Chop[FullSimplify[D[wStar - (aT*c + bT*kq), x]], 10^-25]];

(* --- Audit assertion: rigidity-kernel projection identity (independent quadrature) --- *)
(* Mirror of the SymPy R2 check but via NIntegrate (independent engine + numerical    *)
(* primitive). smallsigma(x) = 2 x is normalized and non-canonical. LHS integrates the *)
(* kernel against (smallSigma - sigmaStar); RHS is the algebraic two-moment formula.   *)
sigmaStarX = pStar*Exp[-pStar*x]/(1 - Exp[-pStar]);
smallSigmaX = 2*x;
normS = NIntegrate[smallSigmaX, {x, 0, 1}, WorkingPrecision -> 40];
normSigma = NIntegrate[sigmaStarX, {x, 0, 1}, WorkingPrecision -> 40];
expectZero["deformation normalized", Chop[(normS - 1) + (normSigma - 1), 10^-30]];
lhsProj = NIntegrate[wStar*(smallSigmaX - sigmaStarX), {x, 0, 1}, WorkingPrecision -> 40];
gBarS = NIntegrate[smallSigmaX*c, {x, 0, 1}, WorkingPrecision -> 40];
sBarS = NIntegrate[smallSigmaX*kq, {x, 0, 1}, WorkingPrecision -> 40];
rhsMoment = aT*(gBarS - gStar) + bT*(sBarS - sStar);
expectZero["kernel projection reproduces two-moment traction shift",
  If[Abs[lhsProj - rhsMoment] < 10^-20, 0, lhsProj - rhsMoment]];

(* --- Audit assertion: source-centering of the rigidity kernel (independent) --- *)
(* CONSULT Q5 (batch 6): the projection identity is blind to the centering constants *)
(* (they vanish against (smallSigma - sigmaStar)); the D[...,x] check only proves the *)
(* offset is constant in x. The kernel's centering condition is orthogonality to the  *)
(* canonical source: integral Sigma_* W_* == 0. Dropping the constants leaves          *)
(* A_T g_* + B_T S_* != 0, so this assertion DOES test them.                           *)
centerResid = NIntegrate[
  Evaluate[SetPrecision[sigmaStarX*wStar, 60]],
  {x, 0, 1},
  WorkingPrecision -> 60,
  AccuracyGoal -> 30,
  PrecisionGoal -> 30,
  MaxRecursion -> 30
];
expectZero["rigidity kernel W_* is source-centered (integral Sigma_* W_* = 0)",
  If[Abs[centerResid] < 10^-20, 0, centerResid]];

(* --- Audit assertion: g_*, S_* equal their source-moment integrals (NIntegrate) --- *)
(* Appendix eq:app-part04-gbar-Sbar: g_* = integral Sigma_*(x) c(x) dx, S_* = integral *)
(* Sigma_*(x) Kq(x) dx, Sigma_*(x) = pStar e^(-pStar x)/(1-e^(-pStar)). Compare to the  *)
(* closed forms gFormula, sFormula at pStar via an independent numerical integrator.    *)
sigmaStarXm = pStar*Exp[-pStar*x]/(1 - Exp[-pStar]);
gStarMoment = NIntegrate[sigmaStarXm*c, {x, 0, 1}, WorkingPrecision -> 40];
sStarMoment = NIntegrate[sigmaStarXm*kq, {x, 0, 1}, WorkingPrecision -> 40];
expectZero["g_* equals source-moment integral",
  If[Abs[gStarMoment - gStar] < 10^-20, 0, gStarMoment - gStar]];
expectZero["S_* equals source-moment integral",
  If[Abs[sStarMoment - sStar] < 10^-20, 0, sStarMoment - sStar]];

Print[""];
Print["Stage 147 Mathematica audit passed."];

Exit[0];
