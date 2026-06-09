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
  res = FullSimplify[Expand[expr], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 185 — DIRECT MICROSCOPIC MONOMIALS"];

Clear[
  chi0s, deltaUs, epsWs, epss, lam1, c1, gam1, kU, keta, kW, mu1, tau1,
  epsVar, lamScale, gammaRef, cetaURef, tuRef, kuRef, ketaRef, kweffRef,
  lamWRef, muWRef, gammaVar, cetaUVar, tuVar, kuVar, ketaVar, kweffVar,
  lamWVar, muWVar, logGamma, logCetaU, logTU, logKU, logKeta, logKWeff,
  logLamW, logMuW
];

$Assumptions =
  Element[
    {
      chi0s, deltaUs, epsWs, epss, lam1, c1, gam1, kU, keta, kW, mu1, tau1,
      epsVar, lamScale, gammaRef, cetaURef, tuRef, kuRef, ketaRef, kweffRef,
      lamWRef, muWRef
    },
    Reals
  ] && chi0s > 0 && deltaUs > 0 && epsWs > 0 && 0 < epss < 1 &&
  gammaRef > 0 && cetaURef > 0 && tuRef > 0 && kuRef > 0 &&
  ketaRef > 0 && kweffRef > 0 && lamWRef > 0 && muWRef > 0;

firstRatioDrift[ratio_] := FullSimplify[(D[ratio, epsVar] /. epsVar -> 0)/lamScale, Assumptions -> $Assumptions];

sigmaChi = gam1 + c1 - kU;
sigmaDelta = tau1 - kU;
sigmaZ = 2*lam1 + mu1 - keta - 2*kW;
sigmaEps = 2*gam1 + 2*lam1 - kU - kW;
sigmaEta = 2*c1 - kU - keta;

sigmaTr = (1 + chi0s)*sigmaDelta + (1 + deltaUs)*sigmaChi;
eStar = 2*epsWs/(1 - epss)*(11 + 9*deltaUs)/(11*(1 + deltaUs));
fStar = 2*chi0s/(1 + deltaUs) + 4*epsWs*deltaUs/(11*(1 - epss)*(1 + deltaUs)^2);
sigmaNt = sigmaZ + eStar*sigmaEps - fStar*sigmaDelta;

Print["Sigma_tr = ", fmt[sigmaTr]];
Print["Sigma_nt = ", fmt[sigmaNt]];
Print["Sigma_eta = ", fmt[sigmaEta]];

gammaP = gammaRef*Exp[epsVar*lamScale*gam1];
cetaUP = cetaURef*Exp[epsVar*lamScale*c1];
tuP = tuRef*Exp[epsVar*lamScale*tau1];
kuP = kuRef*Exp[epsVar*lamScale*kU];
ketaP = ketaRef*Exp[epsVar*lamScale*keta];
kweffP = kweffRef*Exp[epsVar*lamScale*kW];
lamWP = lamWRef*Exp[epsVar*lamScale*lam1];
muWP = muWRef*Exp[epsVar*lamScale*mu1];

gammaRatio = FullSimplify[gammaP/gammaRef, Assumptions -> $Assumptions];
cetaURatio = FullSimplify[cetaUP/cetaURef, Assumptions -> $Assumptions];
tuRatio = FullSimplify[tuP/tuRef, Assumptions -> $Assumptions];
kuRatio = FullSimplify[kuP/kuRef, Assumptions -> $Assumptions];
ketaRatio = FullSimplify[ketaP/ketaRef, Assumptions -> $Assumptions];
kweffRatio = FullSimplify[kweffP/kweffRef, Assumptions -> $Assumptions];
lamWRatio = FullSimplify[lamWP/lamWRef, Assumptions -> $Assumptions];
muWRatio = FullSimplify[muWP/muWRef, Assumptions -> $Assumptions];

chi00 = gammaRef*cetaURef/kuRef;
chi0P = gammaP*cetaUP/kuP;
chiRatio = FullSimplify[chi0P/chi00, Assumptions -> $Assumptions];

deltaU0 = tuRef/kuRef;
deltaUP = tuP/kuP;
deltaURatio = FullSimplify[deltaUP/deltaU0, Assumptions -> $Assumptions];

epsW0 = gammaRef^2*lamWRef^2/(kuRef*kweffRef);
epsWP = gammaP^2*lamWP^2/(kuP*kweffP);
epsWRatio = FullSimplify[epsWP/epsW0, Assumptions -> $Assumptions];

zratio0 = lamWRef^2*muWRef/(ketaRef*kweffRef^2);
zratioP = lamWP^2*muWP/(ketaP*kweffP^2);
zratio = FullSimplify[zratioP/zratio0, Assumptions -> $Assumptions];

epseta0 = FullSimplify[cetaURef^2/(kuRef*ketaRef), Assumptions -> $Assumptions];
epsetaP = cetaUP^2/(kuP*ketaP);
epsetaRatio = FullSimplify[epsetaP/epseta0, Assumptions -> $Assumptions];

sigmaChiDirect = firstRatioDrift[chiRatio];
sigmaDeltaDirect = firstRatioDrift[deltaURatio];
sigmaEpsDirect = firstRatioDrift[epsWRatio];
sigmaZDirect = firstRatioDrift[zratio];
sigmaEtaDirect = firstRatioDrift[epsetaRatio];
sigmaGammaDirect = firstRatioDrift[gammaRatio];
sigmaCetaUDirect = firstRatioDrift[cetaURatio];
sigmaTUDirect = firstRatioDrift[tuRatio];
sigmaKUDirect = firstRatioDrift[kuRatio];
sigmaKetaDirect = firstRatioDrift[ketaRatio];
sigmaKWeffDirect = firstRatioDrift[kweffRatio];
sigmaLamWDirect = firstRatioDrift[lamWRatio];
sigmaMuWDirect = firstRatioDrift[muWRatio];

primitiveVars = {gammaVar, cetaUVar, tuVar, kuVar, ketaVar, kweffVar, lamWVar, muWVar};
primitiveRatios = {gammaRatio, cetaURatio, tuRatio, kuRatio, ketaRatio, kweffRatio, lamWRatio, muWRatio};
primitiveDrifts = {
  sigmaGammaDirect, sigmaCetaUDirect, sigmaTUDirect, sigmaKUDirect,
  sigmaKetaDirect, sigmaKWeffDirect, sigmaLamWDirect, sigmaMuWDirect
};
logVars = {logGamma, logCetaU, logTU, logKU, logKeta, logKWeff, logLamW, logMuW};

monomialExponentVector[monomial_] := Module[{logForm},
  logForm = Expand[
    PowerExpand[
      Log[monomial /. Thread[primitiveVars -> (Exp /@ logVars)]]
    ]
  ];
  FullSimplify[Coefficient[logForm, #], Assumptions -> $Assumptions] & /@ logVars
];

ratioFromExponentVector[vec_] := FullSimplify[
  Times @@ MapThread[#1^#2 &, {primitiveRatios, vec}],
  Assumptions -> $Assumptions
];

driftFromExponentVector[vec_] := FullSimplify[
  vec.primitiveDrifts,
  Assumptions -> $Assumptions
];

chiMonomial = gammaVar*cetaUVar/kuVar;
deltaUMonomial = tuVar/kuVar;
epsWMonomial = gammaVar^2*lamWVar^2/(kuVar*kweffVar);
zMonomial = lamWVar^2*muWVar/(ketaVar*kweffVar^2);
epsEtaMonomial = cetaUVar^2/(kuVar*ketaVar);

banner["Primitive microscopic ratios"];
expectZero["d ln gamma - gamma1", sigmaGammaDirect - gam1];
expectZero["d ln c_etaU - c1", sigmaCetaUDirect - c1];
expectZero["d ln T_U - tau1", sigmaTUDirect - tau1];
expectZero["d ln K_U - kU", sigmaKUDirect - kU];
expectZero["d ln K_eta - keta", sigmaKetaDirect - keta];
expectZero["d ln K_W^(eff) - kW", sigmaKWeffDirect - kW];
expectZero["d ln lambda_W - lambda1", sigmaLamWDirect - lam1];
expectZero["d ln mu_W - mu1", sigmaMuWDirect - mu1];
expectZero["chi0 ratio from primitive ratios", chiRatio - gammaRatio*cetaURatio/kuRatio];
expectZero["delta_U ratio from primitive ratios", deltaURatio - tuRatio/kuRatio];
expectZero[
  "epsilon_W ratio from primitive ratios",
  epsWRatio - gammaRatio^2*lamWRatio^2/(kuRatio*kweffRatio)
];
expectZero[
  "Z_W/Omega_W^2 ratio from primitive ratios",
  zratio - lamWRatio^2*muWRatio/(ketaRatio*kweffRatio^2)
];
expectZero[
  "epsilon_eta ratio from primitive ratios",
  epsetaRatio - cetaURatio^2/(kuRatio*ketaRatio)
];

banner["Microscopic ratio drifts"];
expectZero["d ln chi0 - Sigma_chi", sigmaChiDirect - sigmaChi];
expectZero["d ln delta_U - Sigma_delta", sigmaDeltaDirect - sigmaDelta];
expectZero["d ln epsilon_W - Sigma_eps", sigmaEpsDirect - sigmaEps];
expectZero["d ln (Z_W/Omega_W^2) - Sigma_Z", sigmaZDirect - sigmaZ];
expectZero["d ln epsilon_eta - Sigma_eta", sigmaEtaDirect - sigmaEta];

banner["Tracking monomial"];
ctrRatio = FullSimplify[chiRatio^(1 + deltaUs)*deltaURatio^(1 + chi0s), Assumptions -> $Assumptions];
ctrMonomial = chiMonomial^(1 + deltaUs)*deltaUMonomial^(1 + chi0s);
ctrExponentVector = monomialExponentVector[ctrMonomial];
ctrRatioPrimitive = ratioFromExponentVector[ctrExponentVector];
sigmaTrDirect = firstRatioDrift[ctrRatio];
sigmaTrCompiled = driftFromExponentVector[ctrExponentVector];
expectZero["C_tr,* ratio from primitive coordinates", ctrRatio - ctrRatioPrimitive];
expectZero["d ln C_tr,* (primitive compiler) - Sigma_tr", sigmaTrCompiled - sigmaTr];
expectZero["d ln C_tr,* - Sigma_tr", sigmaTrDirect - sigmaTr];

banner["Nontracking monomial"];
cntRatio = FullSimplify[zratio*epsWRatio^eStar*deltaURatio^(-fStar), Assumptions -> $Assumptions];
cntMonomial = zMonomial*epsWMonomial^eStar*deltaUMonomial^(-fStar);
cntExponentVector = monomialExponentVector[cntMonomial];
cntRatioPrimitive = ratioFromExponentVector[cntExponentVector];
sigmaNtDirect = firstRatioDrift[cntRatio];
sigmaNtCompiled = driftFromExponentVector[cntExponentVector];
expectZero["C_nt,* ratio from primitive coordinates", cntRatio - cntRatioPrimitive];
expectZero["d ln C_nt,* (primitive compiler) - Sigma_nt", sigmaNtCompiled - sigmaNt];
expectZero["d ln C_nt,* - Sigma_nt", sigmaNtDirect - sigmaNt];

banner["Dressing monomial"];
epsEtaExponentVector = monomialExponentVector[epsEtaMonomial];
epsetaRatioPrimitive = ratioFromExponentVector[epsEtaExponentVector];
sigmaEtaCompiled = driftFromExponentVector[epsEtaExponentVector];
expectZero[
  "epsilon_eta ratio from primitive coordinates",
  epsetaRatio - epsetaRatioPrimitive
];
expectZero["d ln epsilon_eta (primitive compiler) - Sigma_eta", sigmaEtaCompiled - sigmaEta];
expectZero["d ln epsilon_eta - Sigma_eta", sigmaEtaDirect - sigmaEta];

banner["Observable triangular law in microscopic monomials"];
cTrStar = chi0s*deltaUs/((1 + chi0s)*(1 + deltaUs)*(1 + chi0s + deltaUs));
aTrStar = 2*chi0s/((1 + chi0s)*(1 + deltaUs));
chi1Indep = FullSimplify[chi0s*sigmaChi, Assumptions -> $Assumptions];
deltaU1Indep = FullSimplify[deltaUs*sigmaDelta, Assumptions -> $Assumptions];
theta1 = FullSimplify[
  -(chi0s*(1 + chi0s)*deltaU1Indep + deltaUs*(1 + deltaUs)*chi1Indep)/
    ((1 + chi0s)*(1 + deltaUs)*(1 + chi0s + deltaUs)),
  Assumptions -> $Assumptions
];
xi1 = FullSimplify[
  sigmaZ + 2*chi0s/(1 + chi0s)*sigmaChi + eStar*sigmaEps -
    4*epsWs*deltaUs/(11*(1 - epss)*(1 + deltaUs)^2)*sigmaDelta,
  Assumptions -> $Assumptions
];
rcombo = FullSimplify[-epseta0/(1 - epseta0)*sigmaEta, Assumptions -> $Assumptions];
rcomboRatio = FullSimplify[(1 - epseta0*epsetaRatio)/(1 - epseta0), Assumptions -> $Assumptions];
rcomboDirect = firstRatioDrift[rcomboRatio];

expectZero["Theta_1 independent slippage law", theta1 - (-cTrStar*sigmaTr)];
expectZero["Xi_1 independent slippage law", xi1 - (aTrStar*sigmaTr + sigmaNt)];
expectZero["Theta_1 monomial law", theta1 - (-cTrStar*sigmaTrCompiled)];
expectZero["Xi_1 monomial law", xi1 - (aTrStar*sigmaTrCompiled + sigmaNtCompiled)];
expectZero["R_1 + Xi_1 complement law", rcomboDirect - rcombo];
Print["Theta1 = ", fmt[theta1]];
Print["Xi1 = ", fmt[xi1]];
Print["R1 + Xi1 = ", fmt[rcombo]];

banner["Exact zero-defect compatibility solve"];
ledgerSigmaTr = sigmaTrCompiled;
ledgerSigmaNt = sigmaNtCompiled;
ledgerSigmaEta = sigmaEtaCompiled;
mStarMinor = {
  {D[ledgerSigmaTr, tau1], D[ledgerSigmaTr, keta], D[ledgerSigmaTr, mu1]},
  {D[ledgerSigmaNt, tau1], D[ledgerSigmaNt, keta], D[ledgerSigmaNt, mu1]},
  {D[ledgerSigmaEta, tau1], D[ledgerSigmaEta, keta], D[ledgerSigmaEta, mu1]}
};
expectZero["det M_*^(tau,keta,mu) - (1+chi0s)", Det[mStarMinor] - (1 + chi0s)];
tauSol = FullSimplify[tau1 /. First[Solve[ledgerSigmaTr == 0, tau1, Reals]], Assumptions -> $Assumptions];
ketaSol = FullSimplify[keta /. First[Solve[ledgerSigmaEta == 0, keta, Reals]], Assumptions -> $Assumptions];
muSol = FullSimplify[mu1 /. First[Solve[ledgerSigmaNt == 0, mu1, Reals]], Assumptions -> $Assumptions];
muSolFull = FullSimplify[muSol /. {tau1 -> tauSol, keta -> ketaSol}, Assumptions -> $Assumptions];

Print["tau1 = ", fmt[tauSol]];
Print["kappa_eta = ", fmt[ketaSol]];
Print["mu1 = ", fmt[muSol]];
Print["mu1 on full zero-defect branch = ", fmt[muSolFull]];

expectZero["tracking substitution", ledgerSigmaTr /. tau1 -> tauSol];
expectZero["dressing substitution", ledgerSigmaEta /. keta -> ketaSol];
expectZero[
  "nontracking substitution",
  FullSimplify[
    ledgerSigmaNt /. {tau1 -> tauSol, keta -> ketaSol, mu1 -> muSolFull},
    Assumptions -> $Assumptions
  ]
];

Print[""];
Print["Carry-forward formulas:"];
Print["  C_tr,*  = chi0^(1+deltaU*) deltaU^(1+chi0*)"];
Print["  C_nt,*  = (Z_W/Omega_W^2) eps_W^(E_*) delta_U^(-F_*)"];
Print["  epsilon_eta = c_{etaU}^2 / (K_U K_eta^(eff))"];
Print["  Zero defect iff these three microscopic monomials are invariant at first grouped order."];

Exit[0];
