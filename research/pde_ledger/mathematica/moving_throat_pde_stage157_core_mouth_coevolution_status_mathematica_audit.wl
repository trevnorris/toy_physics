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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 40] - N[target, 40]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

loadConstants[] := Module[{root, path, raw},
  root = DirectoryName[$InputFileName, 2];
  path = FileNameJoin[{root, "scripts", "numerical", "stage155_156_fixedpoint_samples.json"}];
  raw = Import[path, "RawJSON"];
  raw["constants"]
];

banner["STAGE 157 — CORE-MOUTH COEVOLUTION STATUS"];

constants = loadConstants[];

banner["1. Exact Family-1 compensation identities"];
Clear[r, g];
$Assumptions = Element[{r, g}, Reals] && r > 0;
rFun = FullSimplify[(g - r)^2/(1 + r^2)];
gStar = FullSimplify[r - Sqrt[1 + r^2]/2];
expectZero["R(g_*) - 1/4", (rFun /. g -> gStar) - 1/4];

banner["2. Carry-forward numerical basepoint from Stages 155-156"];
rF1Exact = Sqrt[4107 - 100 Pi^2]/(10 Pi);
gStarExact = FullSimplify[rF1Exact - Sqrt[1 + rF1Exact^2]/2];
rF1 = SetPrecision[constants["rF1"], 30];
gStarNum = SetPrecision[constants["g_star"], 30];
sigma0Star = SetPrecision[constants["Sigma0_star"], 30];
tHatStar = SetPrecision[constants["T_hat_star"], 30];
piStar = SetPrecision[constants["Pi_star"], 30];
sigma0Can = SetPrecision[constants["Sigma0_can_expected"], 30];
sCan = SetPrecision[constants["S_can_expected"], 30];
piCan = SetPrecision[constants["Pi_can_expected"], 30];
tHatCan = SetPrecision[constants["T_hat_can_expected"], 30];

expectApprox["r_F1 radical check", N[rF1Exact, 25], rF1, 10^-12];
expectApprox["g_* lower-branch check", N[gStarExact, 25], gStarNum, 10^-12];
expectApprox["self-matched traction law", sigma0Can, (20/9) tHatCan^2, 10^-11];

rCan = N[((gStarNum - rF1)^2)/(1 + rF1^2), 30];
expectApprox["R_can = 1/4", rCan, 1/4, 10^-12];
expectApprox["Pi_can identity", piCan, sigma0Can (1 - sCan/4), 10^-12];

Print["Sigma0_star = ", fmt[sigma0Star]];
Print["T_hat_star  = ", fmt[tHatStar]];
Print["Pi_star     = ", fmt[piStar]];
Print["Sigma0_can  = ", fmt[sigma0Can]];
Print["T_hat_can   = ", fmt[tHatCan]];
Print["S_can       = ", fmt[sCan]];
Print["Pi_can      = ", fmt[piCan]];

If[!(sigma0Can > sigma0Star && tHatCan > tHatStar && piCan > piStar),
  fail["renormalized point ordering", {sigma0Can, sigma0Star, tHatCan, tHatStar, piCan, piStar}]
];

banner["3. Tangent-on-family and even-preservation handoff"];
Clear[dg, dr];
$Assumptions = Element[{r, dg, dr}, Reals] && r > 0;
s = Sqrt[1 + r^2];
gminus = FullSimplify[r - s/2];
gp = FullSimplify[D[gminus, r]];
dR = FullSimplify[(D[rFun, g] /. g -> gminus) dg + (D[rFun, r] /. g -> gminus) dr];
expectZero["tangent motion keeps delta R = 0", dR /. dg -> gp dr];

Clear[sigmaStar, deltaC, dKappa];
$Assumptions = Element[{sigmaStar, deltaC, dKappa}, Reals];
dE2 = (deltaC - 9 sigmaStar dKappa)/(27 (1 - sigmaStar));
dE4 = (5 deltaC - 72 sigmaStar dKappa)/(243 (1 - sigmaStar));
evenPreservation = Solve[{dE2 == 0, dE4 == 0}, {deltaC, dKappa}, Reals];
Print["canonical-even preservation solutions = ", fmt[evenPreservation]];
If[evenPreservation =!= {{deltaC -> 0, dKappa -> 0}},
  fail["canonical-even preservation", evenPreservation]
];

deltaCFromTangent = FullSimplify[-16 sigmaStar (dR /. dg -> gp dr)];
expectZero["tangent motion kills delta C", deltaCFromTangent];

banner["4. Stage 158 expansion point"];
Clear[sigma0, dSigma0, sStar, dS];
$Assumptions = Element[{sigma0, dSigma0, sStar, dS}, Reals];
piExpr = (sigma0 + dSigma0) (1 - (1/4) (sStar + dS));
piLin = Expand[piExpr] /. dSigma0 dS -> 0;
pi0 = sigma0 (1 - sStar/4);
dPiExpected = Expand[(1 - sStar/4) dSigma0 - sigma0 dS/4];
expectZero["Stage 158 tangent expansion packet", (piLin - pi0) - dPiExpected];
Print["dPi_tan at the renormalized canonical point = ", fmt[dPiExpected /. {sigma0 -> sigma0Can, sStar -> sCan}]];

Print[""];
Print["Open branch note:"];
Print["  This capstone closes the reduced co-evolving classification and basepoint data."];
Print["  It does not assert that the full moving-throat PDE has already realized the"];
Print["  balance-selected branch; that remains the open microscopic question carried forward."];

Print[""];
Print["Stage 157 Mathematica audit passed."];

Exit[0];
