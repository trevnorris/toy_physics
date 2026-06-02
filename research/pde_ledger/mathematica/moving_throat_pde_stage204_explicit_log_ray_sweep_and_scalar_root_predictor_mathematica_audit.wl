ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

normalizeScalar[expr_] := Module[{res},
  res = FullSimplify[PowerExpand[Together[Expand[expr]]], Assumptions -> $Assumptions];
  res = res /. ConditionalExpression[e_, _] :> e;
  FullSimplify[res, Assumptions -> $Assumptions]
];

normalizeExpr[expr_] := If[
  ListQ[expr],
  Map[normalizeScalar, expr, {ArrayDepth[expr]}],
  normalizeScalar[expr]
];

allZeroQ[expr_] := And @@ Flatten[Map[TrueQ[# === 0] &, {expr}, {-1}]];

pretty[expr_] := If[
  MatrixQ[expr],
  MatrixForm[expr],
  If[VectorQ[expr], MatrixForm[{expr}], fmt[expr]]
];

expectZero[name_String, expr_] := Module[{res},
  res = normalizeExpr[expr];
  If[ListQ[res],
    Print[name, " ="];
    Print[pretty[res]];
    If[allZeroQ[res], pass[name], fail[name, res]],
    Print[name, " = ", fmt[res]];
    If[allZeroQ[res], pass[name], fail[name, res]]
  ];
];

banner["STAGE 204 -- EXPLICIT LOG-RAY SWEEP AND SCALAR ROOT PREDICTOR"];

Clear[
  chi0Star, deltaUStar, eStar, fStar, tau, sLambda, sC, sGamma, sU, sW,
  lambda0, c0, gamma0, kU0, kW0, ctrTarget, cntTarget, epsEtaTarget,
  L, sigma, eps, epsStep, L0, Y0, Y1
];

$Assumptions = (
  Element[{tau, sLambda, sC, sGamma, sU, sW, eps, epsStep}, Reals] &&
  chi0Star > 0 && deltaUStar > 0 && eStar > 0 && fStar > 0 &&
  lambda0 > 0 && c0 > 0 && gamma0 > 0 && kU0 > 0 && kW0 > 0 &&
  ctrTarget > 0 && cntTarget > 0 && epsEtaTarget > 0 &&
  L > 0 && sigma > 0 && L0 > 0 && Y0 > 0 && Y1 > 0
);

logRate[expr_] := FullSimplify[
  D[PowerExpand[Log[expr]], tau],
  Assumptions -> $Assumptions
];

lambdaW[z_] := lambda0 Exp[sLambda z];
cEtaU[z_] := c0 Exp[sC z];
gammaRay[z_] := gamma0 Exp[sGamma z];
kURay[z_] := kU0 Exp[sU z];
kWRay[z_] := kW0 Exp[sW z];

subbanner["M1. Constant free log-slopes"];
expectZero["M1 d log lambda_W/dtau - s_lambda", logRate[lambdaW[tau]] - sLambda];
expectZero["M1 d log c_etaU/dtau - s_c", logRate[cEtaU[tau]] - sC];
expectZero["M1 d log gamma/dtau - s_gamma", logRate[gammaRay[tau]] - sGamma];
expectZero["M1 d log K_U/dtau - s_U", logRate[kURay[tau]] - sU];
expectZero["M1 d log K_W/dtau - s_W", logRate[kWRay[tau]] - sW];

subbanner["M2. Dependent delta exponent from graph log derivative"];
aStar = FullSimplify[(1 + deltaUStar)/(1 + chi0Star), Assumptions -> $Assumptions];
deltaGraph[z_] := (
  ctrTarget/((gammaRay[z] cEtaU[z]/kURay[z])^(1 + deltaUStar))
)^(1/(1 + chi0Star));

sigmaDelta = FullSimplify[-aStar (sGamma + sC - sU), Assumptions -> $Assumptions];
sigmaDeltaRecovered = logRate[deltaGraph[tau]];
expectZero["M2 d sigma_delta/dtau", D[sigmaDeltaRecovered, tau]];
expectZero["M2 sigma_delta recovered - formula", sigmaDeltaRecovered - sigmaDelta];

subbanner["M3. Dependent T, K_eta, and mu exponents from graph log derivatives"];
tGraph[z_] := L^2 kURay[z] deltaGraph[z]/Pi^2;
kEtaGraph[z_] := cEtaU[z]^2/(kURay[z] epsEtaTarget);
muGraph[z_] := (
  cntTarget cEtaU[z]^2 kWRay[z]^2/(epsEtaTarget kURay[z] lambdaW[z]^2) *
  ((gammaRay[z]^2 lambdaW[z]^2 sigma)/(kURay[z] kWRay[z]))^(-eStar) *
  deltaGraph[z]^fStar
);

sigmaT = FullSimplify[sU + sigmaDelta, Assumptions -> $Assumptions];
sigmaKeta = FullSimplify[2 sC - sU, Assumptions -> $Assumptions];
sigmaMu = FullSimplify[
  2 sC - sU + 2 sW - 2 sLambda -
  eStar (2 sGamma + 2 sLambda - sU - sW) + fStar sigmaDelta,
  Assumptions -> $Assumptions
];

sigmaTRecovered = logRate[tGraph[tau]];
sigmaKetaRecovered = logRate[kEtaGraph[tau]];
sigmaMuRecovered = logRate[muGraph[tau]];
expectZero["M3 d sigma_T/dtau", D[sigmaTRecovered, tau]];
expectZero["M3 sigma_T recovered - formula", sigmaTRecovered - sigmaT];
expectZero["M3 d sigma_Keta/dtau", D[sigmaKetaRecovered, tau]];
expectZero["M3 sigma_Keta recovered - formula", sigmaKetaRecovered - sigmaKeta];
expectZero["M3 d sigma_mu/dtau", D[sigmaMuRecovered, tau]];
expectZero["M3 sigma_mu recovered - formula", sigmaMuRecovered - sigmaMu];

subbanner["M4. Finite target monomial invariance on the lifted ray"];
deltaFromTiming[kU_, tU_] := Pi^2 tU/(L^2 kU);
ctrMonomial[lambda_, c_, gamma_, kU_, kEta_, kW_, mu_, tU_] := (
  (gamma c/kU)^(1 + deltaUStar) deltaFromTiming[kU, tU]^(1 + chi0Star)
);
cntMonomial[lambda_, c_, gamma_, kU_, kEta_, kW_, mu_, tU_] := (
  lambda^2 mu/(kEta kW^2) *
  ((gamma^2 lambda^2 sigma)/(kU kW))^eStar *
  deltaFromTiming[kU, tU]^(-fStar)
);
epsEtaMonomial[lambda_, c_, gamma_, kU_, kEta_, kW_, mu_, tU_] := c^2/(kU kEta);

xGraphRay = {
  lambdaW[tau],
  cEtaU[tau],
  gammaRay[tau],
  kURay[tau],
  kEtaGraph[0] Exp[sigmaKeta tau],
  kWRay[tau],
  muGraph[0] Exp[sigmaMu tau],
  tGraph[0] Exp[sigmaT tau]
};

ctrRay = ctrMonomial @@ xGraphRay;
cntRay = cntMonomial @@ xGraphRay;
epsEtaRay = epsEtaMonomial @@ xGraphRay;
expectZero["M4 Ctr(tau) - Ctr_target", ctrRay - ctrTarget];
expectZero["M4 Cnt(tau) - Cnt_target", cntRay - cntTarget];
expectZero["M4 eps_eta(tau) - epsEta_target", epsEtaRay - epsEtaTarget];

subbanner["M5. Quotient-map kernel"];
mStar = {
  {0, 1 + deltaUStar, 1 + deltaUStar, -(2 + chi0Star + deltaUStar), 0, 0, 0, 1 + chi0Star},
  {2 (1 + eStar), 0, 2 eStar, fStar - eStar, -1, -(2 + eStar), 1, -fStar},
  {0, 2, 0, -1, -1, 0, 0, 0}
};
dxRay = {sLambda, sC, sGamma, sU, sigmaKeta, sW, sigmaMu, sigmaT};
expectZero["M5 Mstar.dxRay", mStar.dxRay];

subbanner["M6. Primitive free-direction table"];
recoveredExponents = {sigmaDeltaRecovered, sigmaTRecovered, sigmaKetaRecovered, sigmaMuRecovered};
primitiveTable = {
  {"M6 e_lambda", {sLambda -> 1, sC -> 0, sGamma -> 0, sU -> 0, sW -> 0}, {0, 0, 0, -2 - 2 eStar}},
  {"M6 e_c", {sLambda -> 0, sC -> 1, sGamma -> 0, sU -> 0, sW -> 0}, {-aStar, -aStar, 2, 2 - fStar aStar}},
  {"M6 e_gamma", {sLambda -> 0, sC -> 0, sGamma -> 1, sU -> 0, sW -> 0}, {-aStar, -aStar, 0, -2 eStar - fStar aStar}},
  {"M6 e_U", {sLambda -> 0, sC -> 0, sGamma -> 0, sU -> 1, sW -> 0}, {aStar, 1 + aStar, -1, -1 + eStar + fStar aStar}},
  {"M6 e_W", {sLambda -> 0, sC -> 0, sGamma -> 0, sU -> 0, sW -> 1}, {0, 0, 0, 2 + eStar}}
};
Do[
  expectZero[entry[[1]] <> " exponents", (recoveredExponents /. entry[[2]]) - entry[[3]]],
  {entry, primitiveTable}
];

subbanner["M7. First-order root predictors and log-ray completeness"];
phi0 = 1 + eps;
phi1 = L0 (1 + eps);
tauAff = (1 - phi0)/phi1;
tauLog = -Log[phi0]/L0;
predictorDifference = tauLog - tauAff;
Print["M7 series tau_log - tau_aff = ", fmt[Normal[Series[predictorDifference, {eps, 0, 4}]]]];
expectZero["M7 predictor constant coefficient", SeriesCoefficient[predictorDifference, {eps, 0, 0}]];
expectZero["M7 predictor linear coefficient", SeriesCoefficient[predictorDifference, {eps, 0, 1}]];
expectZero["M7 predictor quadratic coefficient", SeriesCoefficient[predictorDifference, {eps, 0, 2}] + 1/(2 L0)];
expectZero["M7 predictor cubic coefficient", SeriesCoefficient[predictorDifference, {eps, 0, 3}] - 2/(3 L0)];

localRayFirstOrder = Normal[Series[Y0 Exp[(Y1/Y0) epsStep], {epsStep, 0, 1}]];
expectZero["M7 first-order log-ray completeness", localRayFirstOrder - (Y0 + Y1 epsStep)];

banner["STAGE 204 MATHEMATICA AUDIT PASSED"];
Exit[0];
