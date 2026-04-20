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

expectMatrixZero[name_String, expr_?MatrixQ] := Module[{res, zero},
  res = Map[FullSimplify[Together[Expand[#]], Assumptions -> $Assumptions] &, expr, {2}];
  zero = ConstantArray[0, Dimensions[res]];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === zero], pass[name], fail[name, res]];
];

banner["STAGE 012 — DYNAMIC LOADING"];

Clear[
  omega, K0, M0, M1, DeltaKax, varpi, OmegaU, OmegaW,
  lambdaB, lambdaU, lambdaW, lambdaR, portPi, GammaPort
];
$Assumptions =
  Element[
    {omega, K0, M0, M1, DeltaKax, varpi, OmegaU, OmegaW, lambdaB, lambdaU, lambdaW, lambdaR, portPi, GammaPort},
    Reals
  ] &&
  K0 > 0 && M0 > 0 && M1 > 0 && DeltaKax > 0 && varpi > 0 &&
  OmegaU > 0 && OmegaW > 0 && GammaPort > 0;

K1 = K0 + DeltaKax;
kappa0 = 2*Sqrt[2]/Pi;
kappa1 = -4/(3*Pi);
kappa0Sq = FullSimplify[kappa0^2, Assumptions -> $Assumptions];
kappa1Sq = FullSimplify[kappa1^2, Assumptions -> $Assumptions];
sigma = FullSimplify[kappa0Sq + kappa1Sq, Assumptions -> $Assumptions];
xiConst = FullSimplify[kappa0Sq - kappa1Sq, Assumptions -> $Assumptions];
eta = FullSimplify[kappa0*kappa1, Assumptions -> $Assumptions];
v = {kappa0, kappa1};
i2 = IdentityMatrix[2];

dbare = {{K0 - M0*omega^2, 0}, {0, K1 - M1*omega^2}};
aphi = varpi^2 - omega^2;
aU = OmegaU^2 - omega^2;
aW = OmegaW^2 - omega^2 - portPi;
deltaUW = FullSimplify[aU*aW - lambdaR^2*sigma, Assumptions -> $Assumptions];
xiShift = FullSimplify[lambdaU^2/aU, Assumptions -> $Assumptions];
alpha = FullSimplify[
  lambdaB^2/aphi + (aU*lambdaW + lambdaR*lambdaU)^2/(aU*deltaUW),
  Assumptions -> $Assumptions
];

Print["D_bare(omega) = ", fmt[dbare]];
Print["Xi(omega) = ", fmt[xiShift]];
Print["alpha(omega) = ", fmt[alpha]];
Print["Delta_UW(omega) = ", fmt[deltaUW]];

mint = {
  {aU, 0, -lambdaR*kappa0, 0},
  {0, aU, -lambdaR*kappa1, 0},
  {-lambdaR*kappa0, -lambdaR*kappa1, aW, 0},
  {0, 0, 0, aphi}
};
cMat = {
  {lambdaU, 0},
  {0, lambdaU},
  {lambdaW*kappa0, lambdaW*kappa1},
  {lambdaB*kappa0, lambdaB*kappa1}
};

sigmaWall = FullSimplify[Transpose[cMat].LinearSolve[mint, cMat], Assumptions -> $Assumptions];
sigmaExpected = FullSimplify[xiShift*i2 + alpha*Outer[Times, v, v], Assumptions -> $Assumptions];
expectMatrixZero["Sigma - (Xi I + alpha vv^T)", sigmaWall - sigmaExpected];
expectZero["sigma - 88/(9 Pi^2)", sigma - 88/(9*Pi^2)];
expectZero["xi - 56/(9 Pi^2)", xiConst - 56/(9*Pi^2)];
expectZero["eta + 8 Sqrt[2]/(3 Pi^2)", eta + 8*Sqrt[2]/(3*Pi^2)];

xi0 = FullSimplify[xiShift /. {omega -> 0, portPi -> 0}, Assumptions -> $Assumptions];
delta0 = FullSimplify[deltaUW /. {omega -> 0, portPi -> 0}, Assumptions -> $Assumptions];
alpha0 = FullSimplify[alpha /. {omega -> 0, portPi -> 0}, Assumptions -> $Assumptions];
k0t = FullSimplify[K0 - xi0, Assumptions -> $Assumptions];
k1t = FullSimplify[K1 - xi0, Assumptions -> $Assumptions];

Print["Xi_0 = ", fmt[xi0]];
Print["Delta_0 = ", fmt[delta0]];
Print["alpha_0 = ", fmt[alpha0]];
expectZero["DeltaK_tilde - DeltaK_ax", (k1t - k0t) - DeltaKax];

keff0 = {
  {k0t - alpha0*kappa0Sq, -alpha0*kappa0*kappa1},
  {-alpha0*kappa0*kappa1, k1t - alpha0*kappa1Sq}
};
theta = Symbol["theta"];
q = {Cos[theta], Sin[theta]};
energy = FullSimplify[(q.keff0.q)/2, Assumptions -> $Assumptions];
dEnergy = FullSimplify[TrigExpand[D[energy, theta]], Assumptions -> $Assumptions];
tan2Theta = FullSimplify[2*alpha0*eta/(DeltaKax + alpha0*xiConst), Assumptions -> $Assumptions];
stationarity = FullSimplify[
  (DeltaKax + alpha0*xiConst)*Sin[2*theta] - 2*alpha0*eta*Cos[2*theta],
  Assumptions -> $Assumptions
];

expectZero["dE/dtheta - stationarity/2", dEnergy - stationarity/2];
expectZero[
  "-tan(2 theta_-) - manifestly positive form",
  -tan2Theta - 2*alpha0*(-eta)/(DeltaKax + alpha0*xiConst)
];
Print["tan(2 theta_-) = ", fmt[tan2Theta]];

al = Symbol["alphaLoad"];
detTemplate = FullSimplify[
  Det[{{k0t - al*kappa0Sq, -al*kappa0*kappa1}, {-al*kappa0*kappa1, k1t - al*kappa1Sq}}],
  Assumptions -> $Assumptions
];
alphaCrit = FullSimplify[k0t*k1t/(k1t*kappa0Sq + k0t*kappa1Sq), Assumptions -> $Assumptions];
Print["alpha_crit = ", fmt[alphaCrit]];
expectZero["det(alpha_crit)", detTemplate /. al -> alphaCrit];

alphaCons = FullSimplify[alpha /. portPi -> 0, Assumptions -> $Assumptions];
beta = FullSimplify[D[alpha, portPi] /. portPi -> 0, Assumptions -> $Assumptions];
deltaCons = FullSimplify[aU*(OmegaW^2 - omega^2) - lambdaR^2*sigma, Assumptions -> $Assumptions];
betaClean = FullSimplify[(aU*lambdaW + lambdaR*lambdaU)^2/deltaCons^2, Assumptions -> $Assumptions];

Print["alpha_cons(omega) = ", fmt[alphaCons]];
Print["beta(omega) = ", fmt[beta]];
expectZero["beta - clean transfer factor", beta - betaClean];

eps = Symbol["eps"];
alphaPiSeries = FullSimplify[Normal[Series[alpha /. portPi -> eps, {eps, 0, 1}]], Assumptions -> $Assumptions];
expectZero["alpha - (alpha_cons + beta portPi) at O(portPi)", (alphaPiSeries /. eps -> portPi) - (alphaCons + betaClean*portPi)];

alphaOut = FullSimplify[alpha /. portPi -> I*GammaPort*omega^5, Assumptions -> $Assumptions];
alphaOutSeries = FullSimplify[Normal[Series[alphaOut, {omega, 0, 5}]], Assumptions -> $Assumptions];
beta5 = FullSimplify[(betaClean /. omega -> 0)*GammaPort, Assumptions -> $Assumptions];
beta5Target = FullSimplify[GammaPort*(OmegaU^2*lambdaW + lambdaR*lambdaU)^2/delta0^2, Assumptions -> $Assumptions];

Print["alpha_out(omega) through O(omega^5) = ", fmt[alphaOutSeries]];
Print["beta_5 = ", fmt[beta5]];
expectZero["beta_5 - GammaPort (OmegaU^2 lambdaW + lambdaR lambdaU)^2/Delta0^2", beta5 - beta5Target];
expectZero["extracted beta_5 - expected beta_5", Coefficient[alphaOutSeries, omega, 5]/I - beta5];

kappaThetaSq = FullSimplify[q.Outer[Times, v, v].q, Assumptions -> $Assumptions];
Print["kappa(theta)^2 = ", fmt[kappaThetaSq]];

discTemplate = FullSimplify[(DeltaKax + al*xiConst)^2 + 4*al^2*eta^2, Assumptions -> $Assumptions];
trTemplate = FullSimplify[k0t + k1t - al*sigma, Assumptions -> $Assumptions];
lambdaMinusTemplate = FullSimplify[(trTemplate - Sqrt[discTemplate])/2, Assumptions -> $Assumptions];
kappaSelSq = FullSimplify[-D[lambdaMinusTemplate, al], Assumptions -> $Assumptions];

Print["kappa_sel^2 = ", fmt[kappaSelSq]];
expectZero["weak-loading kappa_sel^2 - kappa0^2", (kappaSelSq /. al -> 0) - kappa0Sq];
expectZero["strong-loading kappa_sel^2 - sigma", FullSimplify[Limit[kappaSelSq, al -> Infinity], Assumptions -> DeltaKax > 0] - sigma];

oddProjection = FullSimplify[-I*beta5*kappaSelSq*omega^5, Assumptions -> $Assumptions];
Print["delta D_-^(odd)(omega) template = ", fmt[oddProjection]];

Print[""];
Print["Stage 012 Mathematica audit passed."];

Exit[0];
