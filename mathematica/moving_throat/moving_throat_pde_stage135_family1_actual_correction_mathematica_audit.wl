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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 50] - N[target, 50]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["ACTUAL FAMILY-1 MOUTH CORRECTION"];

Clear[p, t];

kap = Pi/2;
gFormula[p_] := 2*p*(2*p*Exp[p] + Pi)/((4*p^2 + Pi^2)*(Exp[p] - 1));
sFormula[p_] := p*(kap*Tanh[kap] + p*(Exp[-p]*Sech[kap] - 1))/((1 - Exp[-p])*(kap^2 - p^2));
rF1 = Sqrt[(12*(37/20)^2)/Pi^2 - 1];
gMinus = rF1 - Sqrt[1 + rF1^2]/2;

pStar = p /. FindRoot[gFormula[p] == gMinus, {p, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
sStar = N[sFormula[pStar], 50];
gStar = N[gFormula[pStar], 50];
gPrimeStar = N[D[gFormula[p], p] /. p -> pStar, 50];
sPrimeStar = N[D[sFormula[p], p] /. p -> pStar, 50];
sigmaMStar = N[pStar/(4 - sStar), 50];
tStar = N[Sqrt[9*(pStar/(1 - sStar/4))/20], 50];
aT = N[
  -(9/(40*tStar))*(1/(gPrimeStar*(1 - sStar/4)) + pStar*sPrimeStar/(4*gPrimeStar*(1 - sStar/4)^2)),
  40
];
bT = N[(9/(40*tStar))*pStar/(4*(1 - sStar/4)^2), 40];

ts[x_] := (1 - Exp[-pStar*x])/(pStar*(1 - Exp[-pStar])) - x*Exp[-pStar]/(1 - Exp[-pStar]);
tq[x_] := Module[{cQ, aQ},
  cQ = pStar/((1 - Exp[-pStar])*(kap^2 - pStar^2));
  aQ = cQ*(kap*Sinh[kap] + pStar*Exp[-pStar])/(kap*Cosh[kap]);
  aQ*Sinh[kap*x] - cQ*Cosh[kap*x] + cQ*Exp[-pStar*x]
];
sQKernel[x_] := Cosh[(Pi/2)*(1 - x)]/Cosh[Pi/2];
cKernel[x_] := Cos[Pi*x/2];
sigmaStar[x_] := pStar*Exp[-pStar*x]/(1 - Exp[-pStar]);
rStar[x_] := sigmaMStar*(4*ts[x] - tq[x]) - pStar*x;

Print["R'(0) forward diff at h=1e-6 : ", fmt[N[(rStar[10^-6] - rStar[0])/10^-6, 30]]];
Print["R'(0) forward diff at h=1e-5 : ", fmt[N[(rStar[10^-5] - rStar[0])/10^-5, 30]]];
Print["R(0) = ", fmt[N[rStar[0], 30]]];
Print["R(1) = ", fmt[N[rStar[1], 30]]];

nint[expr_] := Quiet[
  NIntegrate[
    Evaluate[expr],
    {t, 0, 1},
    WorkingPrecision -> 60,
    AccuracyGoal -> 20,
    PrecisionGoal -> 20,
    Method -> {"GlobalAdaptive", "MaxErrorIncreases" -> 100}
  ],
  NIntegrate::precw
];

avgR = nint[sigmaStar[t]*rStar[t]];
avgC = nint[sigmaStar[t]*cKernel[t]];
avgS = nint[sigmaStar[t]*sQKernel[t]];

covCR = nint[sigmaStar[t]*(cKernel[t] - avgC)*(rStar[t] - avgR)];
covSR = nint[sigmaStar[t]*(sQKernel[t] - avgS)*(rStar[t] - avgR)];

deltaG = -covCR;
deltaS = -covSR;
deltaPi = covCR/gPrimeStar;
deltaT = aT*deltaG + bT*deltaS;

Print["Cov_*(c,R_*)      = ", fmt[N[covCR, 30]]];
Print["Cov_*(K_q,R_*)    = ", fmt[N[covSR, 30]]];
Print["delta g_act       = ", fmt[N[deltaG, 30]]];
Print["delta S_act       = ", fmt[N[deltaS, 30]]];
Print["delta Pi_act      = ", fmt[N[deltaPi, 30]]];
Print["delta Tm_act      = ", fmt[N[deltaT, 30]]];
Print["Pi_corr           = ", fmt[N[pStar + deltaPi, 30]]];
Print["T_corr            = ", fmt[N[tStar + deltaT, 30]]];

sigma1Unnorm[x_] := Exp[-pStar*x - rStar[x]];
z1 = nint[sigma1Unnorm[t]];
g1 = nint[(sigma1Unnorm[t]/z1)*cKernel[t]];
s1 = nint[(sigma1Unnorm[t]/z1)*sQKernel[t]];
deltaPi1 = -(g1 - gStar)/gPrimeStar;
deltaT1 = aT*(g1 - gStar) + bT*(s1 - sStar);

Print["g_1               = ", fmt[N[g1, 30]]];
Print["S_1               = ", fmt[N[s1, 30]]];
Print["Pi_1              = ", fmt[N[pStar + deltaPi1, 30]]];
Print["T_1               = ", fmt[N[tStar + deltaT1, 30]]];

gU = N[2/Pi, 50];
sU = N[2*Tanh[Pi/2]/Pi, 50];
gD = N[Pi/4, 50];
sD = N[(1 + Sinh[Pi/2])/(2*Cosh[Pi/2]), 50];
deltaPiU = -(gU - gStar)/gPrimeStar;
deltaPiD = -(gD - gStar)/gPrimeStar;
deltaTU = aT*(gU - gStar) + bT*(sU - sStar);
deltaTD = aT*(gD - gStar) + bT*(sD - sStar);
lambdaEffPi = (deltaPi - deltaPiU)/(deltaPiD - deltaPiU);
lambdaEffT = (deltaT - deltaTU)/(deltaTD - deltaTU);

Print["lambda_eff^(Pi)   = ", fmt[N[lambdaEffPi, 30]]];
Print["lambda_eff^(T)    = ", fmt[N[lambdaEffT, 30]]];

linearizedDeltaG = nint[sigmaStar[t]*(1 - (rStar[t] - avgR))*cKernel[t]] - gStar;
expectApprox["delta g linearized covariance consistency", linearizedDeltaG, deltaG, 10^-9];
expectApprox["delta Pi_act scale", deltaPi, 0.907084414842908, 10^-6];
expectApprox["delta Tm_act scale", deltaT, 0.271653979462338, 10^-6];
expectApprox["lambda_eff^(Pi) scale", lambdaEffPi, 0.380487632771110, 10^-6];
expectApprox["lambda_eff^(T) scale", lambdaEffT, 0.378939241176339, 10^-6];

Print[""];
Print["Conclusion:"];
Print["  The full compensated mouth potential broadens the source relative to the tangent exponential,"];
Print["  shifts the canonical Family-1 point upward, and is well approximated by a moderate positive"];
Print["  broadening toward the uniform family (lambda_eff ≈ 0.38)."];

Print[""];
Print["Stage 135 Mathematica audit passed."];

Exit[0];
