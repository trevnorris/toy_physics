ClearAll["Global`*"];
$HistoryLength = 0;
$MaxExtraPrecision = 1000;

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  detail -> ", fmt[detail]]];
  Exit[1];
);

n80[s_String] := Module[{parts},
  parts = StringSplit[s, {"e", "E"}];
  If[Length[parts] == 2,
    ToExpression[parts[[1]] <> "`80"] 10^ToExpression[parts[[2]]],
    ToExpression[s <> "`80"]
  ]
];

expectClose[name_String, actual_, expected_, tol_] := Module[{diff},
  diff = N[Abs[actual - expected], 50];
  Print[name, " residual = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

expectSmall[name_String, actual_, tol_] := Module[{res},
  res = N[Abs[actual], 50];
  Print[name, " residual = ", fmt[res]];
  If[TrueQ[res <= tol], pass[name], fail[name, res]];
];

expectPositive[name_String, actual_] := Module[{res},
  res = N[actual, 50];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res > 0], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := (
  Print[name, " = ", fmt[cond]];
  If[TrueQ[cond], pass[name], fail[name, cond]];
);

Print["=== Stage 232 Mathematica audit: known 5PN data injection and branch verdict ==="];

(* M1: refreshed geometry. *)
x01 = N[BesselJZero[0, 1], 100];
LambdaEll = N[20 Sqrt[2] Pi/x01, 100];
chiS = N[LambdaEll/2, 100];
eta = LambdaEll;
kappa = N[(9/5) LambdaEll^2, 100];

Print["M1 x01 = ", fmt[N[x01, 30]]];
Print["M1 Lambda_ell = ", fmt[N[LambdaEll, 30]]];
Print["M1 chi_s = ", fmt[N[chiS, 30]]];
Print["M1 eta = ", fmt[N[eta, 30]]];
Print["M1 kappa = ", fmt[N[kappa, 30]]];
expectClose["M1 Lambda_ell decimal", LambdaEll, n80["36.94973154240256"], 5 10^-12];
expectClose["M1 kappa decimal", kappa, n80["2457.508789900114"], 5 10^-12];

(* M2: Robin support ceiling. *)
yRoot = y /. FindRoot[
  y Tan[y] == eta,
  {y, n80["1.5"]},
  WorkingPrecision -> 90,
  AccuracyGoal -> 55,
  PrecisionGoal -> 55
];
AK = N[(kappa + Pi^2/4)/(kappa + yRoot^2), 90];
zetaMax = N[AK Pi^2/4, 90];

Print["M2 y root = ", fmt[N[yRoot, 30]]];
Print["M2 A_K = ", fmt[N[AK, 30]]];
Print["M2 zeta_max = ", fmt[N[zetaMax, 30]]];
expectTrue["M2 lowest Robin root lies in (0, Pi/2)", 0 < yRoot < Pi/2];
expectSmall["M2 y Tan[y] - eta", yRoot Tan[yRoot] - eta, 5 10^-55];
expectClose["M2 y decimal", yRoot, n80["1.5294278190457656"], 5 10^-13];
expectClose["M2 A_K decimal", AK, n80["1.0000521380385143"], 5 10^-13];
expectClose["M2 zeta_max decimal", zetaMax, n80["2.4675297457259358"], 5 10^-13];

(* M3: support-drop kernel and endpoint limits by native integration. *)
alpha = N[Sqrt[kappa], 100];
denom = N[alpha Sinh[alpha] + eta Cosh[alpha], 100];
kernel[x_] := (
  Cosh[alpha x] + (eta/alpha) Sinh[alpha x] - Cosh[alpha (1 - x)]
)/denom;
sourceStable[x_, p_] := p Exp[p (x - 1)]/(1 - Exp[-p]);

Clear[pVar, xVar];
deltaIntegralExpr = Integrate[
  kernel[xVar] sourceStable[xVar, pVar],
  {xVar, 0, 1},
  GenerateConditions -> False,
  Assumptions -> pVar > 0 && pVar != alpha
];
If[!FreeQ[deltaIntegralExpr, Integrate],
  fail["M3 native Integrate produced an unevaluated integral", deltaIntegralExpr]
];

deltaByIntegral[p_?NumericQ] := N[deltaIntegralExpr /. pVar -> SetPrecision[p, 80], 70];

delta0Formula = N[eta (Cosh[alpha] - 1)/(alpha^2 denom), 90];
deltaInfFormula = N[(Cosh[alpha] + (eta/alpha) Sinh[alpha] - 1)/denom, 90];
delta0FromLimit = N[Limit[deltaIntegralExpr, pVar -> 0, Direction -> "FromAbove"], 80];
deltaInfFromLimit = N[Limit[deltaIntegralExpr, pVar -> Infinity], 80];
delta0FromUniformSource = N[Integrate[kernel[xVar], {xVar, 0, 1}, GenerateConditions -> False], 80];

Print["M3 Delta integral expression head = ", fmt[Head[deltaIntegralExpr]]];
Print["M3 Delta_0 = ", fmt[N[delta0Formula, 30]]];
Print["M3 Delta_inf = ", fmt[N[deltaInfFormula, 30]]];
expectClose["M3 Delta_0 formula vs Pe->0 integral limit", delta0Formula, delta0FromLimit, 5 10^-35];
expectClose["M3 Delta_0 formula vs uniform-source integral", delta0Formula, delta0FromUniformSource, 5 10^-35];
expectClose["M3 Delta_inf formula vs Pe->infinity integral limit", deltaInfFormula, deltaInfFromLimit, 5 10^-35];
expectClose["M3 Delta_0 decimal", delta0Formula, n80["1.7377393923469950e-4"], 5 10^-16];
expectClose["M3 Delta_inf decimal", deltaInfFormula, n80["2.0172162594593645e-2"], 5 10^-16];

(* M4: 100-based figures of merit. *)
ThetaWChi = n80["4.06863235008162"];
ThetaWJ = n80["0.927552032539308"];
cPrefactor = 100;
XiChi = N[cPrefactor ThetaWChi LambdaEll^2, 90];
XiJ = N[cPrefactor ThetaWJ LambdaEll^2, 90];

Print["M4 Xi_chi = ", fmt[N[XiChi, 30]]];
Print["M4 Xi_J = ", fmt[N[XiJ, 30]]];
expectClose["M4 Xi_chi decimal", XiChi, n80["5.5548332017764099e5"], 5 10^-8];
expectClose["M4 Xi_J decimal", XiJ, n80["1.2663707072528143e5"], 5 10^-8];

(* M5: fixed-point roots from the integral-derived Delta(Pe). *)
solveBranch[xi_, guess_, label_String] := Module[{lo, hi, root, residual},
  lo = N[xi delta0Formula, 80];
  hi = N[xi deltaInfFormula, 80];
  root = pe /. FindRoot[
    pe == xi deltaByIntegral[pe],
    {pe, SetPrecision[guess, 80]},
    WorkingPrecision -> 80,
    AccuracyGoal -> 40,
    PrecisionGoal -> 40,
    MaxIterations -> 100
  ];
  residual = root - xi deltaByIntegral[root];
  Print[label, " bracket lo = ", fmt[N[lo, 25]], " root = ", fmt[N[root, 25]], " hi = ", fmt[N[hi, 25]]];
  expectTrue[label <> " root lies in support bracket", lo <= root <= hi];
  expectSmall[label <> " Pe - Xi Delta(Pe)", residual, 1 10^-35];
  {root, lo, hi}
];

{PeChi, loChi, hiChi} = solveBranch[XiChi, n80["11155"], "M5 chi"];
{PeJ, loJ, hiJ} = solveBranch[XiJ, n80["2505"], "M5 J"];

expectClose["M5 Pe_*^(chi) decimal", PeChi, n80["11155.7265863205869"], 5 10^-10];
expectClose["M5 Pe_*^(J) decimal", PeJ, n80["2504.9703142859238"], 5 10^-10];

(* M6: physical support ratios. *)
omegaPe[p_] := Pi p (2 p + Pi Exp[-p])/((4 p^2 + Pi^2) (1 - Exp[-p]));
zetaPhysChi = N[AK omegaPe[PeChi]^2, 90];
zetaPhysJ = N[AK omegaPe[PeJ]^2, 90];
rhoAlphaMaxChi = N[1 + zetaPhysChi, 90];
rhoAlphaMaxJ = N[1 + zetaPhysJ, 90];

Print["M6 zeta_phys^(chi) = ", fmt[N[zetaPhysChi, 30]]];
Print["M6 rho_alpha,max^(chi) = ", fmt[N[rhoAlphaMaxChi, 30]]];
Print["M6 zeta_phys^(J) = ", fmt[N[zetaPhysJ, 30]]];
Print["M6 rho_alpha,max^(J) = ", fmt[N[rhoAlphaMaxJ, 30]]];
expectClose["M6 zeta_phys^(chi) decimal", zetaPhysChi, n80["2.4675296478814376"], 5 10^-13];
expectClose["M6 zeta_phys^(J) decimal", zetaPhysJ, n80["2.467527805167508"], 5 10^-13];

(* M7: injected support margins and ratios. *)
zetaReq = 1/3;
rhoAlphaReq = 4/3;
marginZetaChi = N[zetaPhysChi - zetaReq, 90];
marginZetaJ = N[zetaPhysJ - zetaReq, 90];
marginRhoChi = N[rhoAlphaMaxChi - rhoAlphaReq, 90];
marginRhoJ = N[rhoAlphaMaxJ - rhoAlphaReq, 90];
gapToCeilingChi = N[zetaMax - zetaPhysChi, 90];
gapToCeilingJ = N[zetaMax - zetaPhysJ, 90];
ratioZetaChi = N[zetaPhysChi/zetaReq, 90];
ratioZetaJ = N[zetaPhysJ/zetaReq, 90];
ratioRhoChi = N[rhoAlphaMaxChi/rhoAlphaReq, 90];
ratioRhoJ = N[rhoAlphaMaxJ/rhoAlphaReq, 90];

Print["M7 margin_zeta^(chi) = ", fmt[N[marginZetaChi, 30]]];
Print["M7 margin_rho^(chi) = ", fmt[N[marginRhoChi, 30]]];
Print["M7 margin_zeta^(J) = ", fmt[N[marginZetaJ, 30]]];
Print["M7 margin_rho^(J) = ", fmt[N[marginRhoJ, 30]]];
Print["M7 ceiling gap chi = ", fmt[N[gapToCeilingChi, 30]]];
Print["M7 ceiling gap J = ", fmt[N[gapToCeilingJ, 30]]];
expectPositive["M7 zeta margin chi", marginZetaChi];
expectPositive["M7 zeta margin J", marginZetaJ];
expectPositive["M7 rho margin chi", marginRhoChi];
expectPositive["M7 rho margin J", marginRhoJ];
expectPositive["M7 ceiling gap chi", gapToCeilingChi];
expectPositive["M7 ceiling gap J", gapToCeilingJ];
expectClose["M7 margin_zeta^(chi) decimal", marginZetaChi, n80["2.1341963145481043"], 5 10^-13];
expectClose["M7 margin_zeta^(J) decimal", marginZetaJ, n80["2.1341944718341751"], 5 10^-13];
expectClose["M7 margin_rho^(chi) decimal", marginRhoChi, n80["2.1341963145481043"], 5 10^-13];
expectClose["M7 margin_rho^(J) decimal", marginRhoJ, n80["2.1341944718341751"], 5 10^-13];
expectClose["M7 ceiling gap chi decimal", gapToCeilingChi, n80["9.784449817674381e-8"], 5 10^-13];
expectClose["M7 ceiling gap J decimal", gapToCeilingJ, n80["1.9405584273504838e-6"], 5 10^-13];
expectClose["M7 zeta ratio chi decimal", ratioZetaChi, n80["7.402588943644313"], 5 10^-13];
expectClose["M7 zeta ratio J decimal", ratioZetaJ, n80["7.402583415502525"], 5 10^-13];
expectClose["M7 rho ratio chi decimal", ratioRhoChi, n80["2.600647235911078"], 5 10^-13];
expectClose["M7 rho ratio J decimal", ratioRhoJ, n80["2.600645853875631"], 5 10^-13];

Print["All Stage 232 Mathematica checks passed."];
Exit[0];
