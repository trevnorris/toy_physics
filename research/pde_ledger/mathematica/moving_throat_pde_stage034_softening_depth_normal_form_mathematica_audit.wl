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

banner["STAGE 017 — SOFTENING-DEPTH NORMAL FORM"];

Clear[kappa0Sq, kappa1Sq, A, DeltaK, beta0, lambda, x];
$Assumptions =
  Element[{kappa0Sq, kappa1Sq, A, DeltaK, beta0, lambda, x}, Reals] &&
  kappa0Sq > 0 && kappa1Sq > 0 && A > 0 && DeltaK > 0 && beta0 > 0 &&
  0 < lambda < A && 0 <= x < A;

alphaLambda = FullSimplify[
  1/(kappa0Sq/(A - lambda) + kappa1Sq/(A + DeltaK - lambda)),
  Assumptions -> $Assumptions
];
s1 = FullSimplify[
  kappa0Sq/(A - lambda) + kappa1Sq/(A + DeltaK - lambda),
  Assumptions -> $Assumptions
];
s1p = FullSimplify[D[s1, lambda], Assumptions -> $Assumptions];
sLambda = FullSimplify[s1^2/s1p, Assumptions -> $Assumptions];
nLambda = FullSimplify[beta0*sLambda^2/(kappa0Sq*lambda), Assumptions -> $Assumptions];

alphaX = FullSimplify[
  x*(x + DeltaK)/(kappa0Sq*(x + DeltaK) + kappa1Sq*x),
  Assumptions -> $Assumptions
];
sX = FullSimplify[
  (kappa0Sq*(x + DeltaK) + kappa1Sq*x)^2/(kappa0Sq*(x + DeltaK)^2 + kappa1Sq*x^2),
  Assumptions -> $Assumptions
];
nX = FullSimplify[
  beta0*(kappa0Sq*(x + DeltaK) + kappa1Sq*x)^4/
    (kappa0Sq*(A - x)*(kappa0Sq*(x + DeltaK)^2 + kappa1Sq*x^2)^2),
  Assumptions -> $Assumptions
];

Print["alpha(lambda) = ", fmt[alphaLambda]];
Print["alpha(x) = ", fmt[alphaX]];
Print["s_-(x) = ", fmt[sX]];
Print["N_-(x) = ", fmt[nX]];

expectZero["alpha(lambda=A-x) - alpha(x)", (alphaLambda /. lambda -> A - x) - alphaX];
expectZero["s(lambda=A-x) - s(x)", (sLambda /. lambda -> A - x) - sX];
expectZero["N(lambda=A-x) - N(x)", (nLambda /. lambda -> A - x) - nX];

dAlphaDx = FullSimplify[D[alphaX, x], Assumptions -> $Assumptions];
dAlphaTarget = FullSimplify[
  (kappa0Sq*(x + DeltaK)^2 + kappa1Sq*x^2)/(kappa0Sq*(x + DeltaK) + kappa1Sq*x)^2,
  Assumptions -> $Assumptions
];
Print["d alpha / dx = ", fmt[dAlphaDx]];
expectZero["dalpha/dx - manifestly positive form", dAlphaDx - dAlphaTarget];
expectZero["s_x * d alpha / dx - 1", sX*dAlphaDx - 1];

Clear[gBsqOverVarpi2, Chi, OmegaU, Delta0];
$Assumptions =
  Element[{kappa0Sq, kappa1Sq, A, DeltaK, beta0, x, gBsqOverVarpi2, Chi, OmegaU, Delta0}, Reals] &&
  kappa0Sq > 0 && kappa1Sq > 0 && A > 0 && DeltaK > 0 && beta0 > 0 &&
  0 <= x < A && OmegaU > 0 && Delta0 > 0;

alphaMix = FullSimplify[Chi^2/(OmegaU^2*Delta0), Assumptions -> $Assumptions];
gBSolution = Solve[alphaMix + gBsqOverVarpi2 == alphaX, gBsqOverVarpi2, Reals];
If[Length[gBSolution] != 1, fail["support-loading solve count", Length[gBSolution]]];
gBReqSqOverVarpi2 = FullSimplify[gBsqOverVarpi2 /. First[gBSolution], Assumptions -> $Assumptions];

Print["alpha_mix = ", fmt[alphaMix]];
Print["g_B,req^2 / varpi^2 = ", fmt[gBReqSqOverVarpi2]];
expectZero[
  "solved g_B,req^2/varpi^2 - (alpha(x) - alpha_mix)",
  gBReqSqOverVarpi2 - (alphaX - alphaMix)
];
expectZero[
  "alpha_mix + g_B,req^2/varpi^2 - alpha(x)",
  alphaMix + gBReqSqOverVarpi2 - alphaX
];

Print[""];
Print["Stage 034 Mathematica audit passed."];

Exit[0];
