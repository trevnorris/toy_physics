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

banner["FIRST-ORDER SELF-CONSISTENT SOURCE CORRECTION"];

Clear[x, piStar, r1, r2, gPrime, aT, bT, epsilon];
$Assumptions = piStar > 0 && Element[{r1, r2, gPrime, aT, bT}, Reals] && gPrime != 0;

Phi[x_] := piStar*x + epsilon*(r1*x + r2*x^2);
unnorm[x_] := Exp[-Phi[x]];
Z = Integrate[unnorm[x], {x, 0, 1}, Assumptions -> piStar > 0];
SigmaFull[x_] := unnorm[x]/Z;
SigmaSeries = Normal[Series[SigmaFull[x], {epsilon, 0, 1}]];
SigmaStar = Coefficient[SigmaSeries, epsilon, 0];
deltaSigma = Coefficient[SigmaSeries, epsilon, 1];

cKernel[x_] := Cos[Pi*x/2];
kKernel[x_] := Cosh[Pi*(1 - x)/2]/Cosh[Pi/2];
RResidual[x_] := r1*x + r2*x^2;

mean[f_] := Integrate[SigmaStar*f, {x, 0, 1}, Assumptions -> piStar > 0];

rBar = FullSimplify[mean[RResidual[x]]];
cBar = FullSimplify[mean[cKernel[x]]];
kBar = FullSimplify[mean[kKernel[x]]];
cRBar = FullSimplify[mean[cKernel[x]*RResidual[x]]];
kRBar = FullSimplify[mean[kKernel[x]*RResidual[x]]];
covCR = FullSimplify[cRBar - cBar*rBar];
covKR = FullSimplify[kRBar - kBar*rBar];

expectZero["<deltaSigma>_*  (centering, from Series)", Integrate[deltaSigma, {x, 0, 1}, Assumptions -> piStar > 0]];
expectZero["deltaSigma + SigmaStar*(R - <R>)", deltaSigma + SigmaStar*(RResidual[x] - rBar)];

deltaGInt = Integrate[cKernel[x]*deltaSigma, {x, 0, 1}, Assumptions -> piStar > 0];
deltaSInt = Integrate[kKernel[x]*deltaSigma, {x, 0, 1}, Assumptions -> piStar > 0];
expectZero["deltaGInt + Cov(c,R)", deltaGInt + covCR];
expectZero["deltaSInt + Cov(K,R)", deltaSInt + covKR];

deltaPi = -deltaGInt/gPrime;
deltaT = aT*deltaGInt + bT*deltaSInt;
expectZero["deltaPi - Cov(c,R)/gPrime", deltaPi - covCR/gPrime];
expectZero["deltaT + aT*Cov(c,R) + bT*Cov(K,R)", deltaT + aT*covCR + bT*covKR];

Print[""];
Print["Theorem:"];
Print["  Once the full mouth residual R_*(x) is known, the selected first-order"];
Print["  source correction is completely determined by Cov_*(c,R_*) and Cov_*(K_q,R_*)."];

Print[""];
Print["Stage 151 Mathematica audit passed."];

Exit[0];
