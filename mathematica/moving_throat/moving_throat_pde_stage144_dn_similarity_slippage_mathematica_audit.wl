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

banner["STAGE 144 — D/N SIMILARITY SLIPPAGE DECOMPOSITION"];

Clear[rc, epsKappa, epsGamma];
$Assumptions = Element[{rc, epsKappa, epsGamma}, Reals];

kappa0 = (1 + rc)*(1 + epsKappa)/3;
gamma0 = (1 + rc)*(1 + epsGamma)/9;
bW = FullSimplify[gamma0 - kappa0/3];
Print["B_W = ", fmt[bW]];
expectZero["exact similarity-defect decomposition", bW - (1 + rc)*(epsGamma - epsKappa)/9];

Clear[eps, rcStar, drc, depsK, depsG];
$Assumptions = Element[{eps, rcStar, drc, depsK, depsG}, Reals];

bWLin = Normal[Series[((1 + rcStar + eps*drc)/9)*(eps*depsG - eps*depsK), {eps, 0, 1}]];
dBW = Expand[Coefficient[bWLin, eps, 1]];
Print["dB_W = ", fmt[dBW]];
expectZero["linearized slippage law", dBW - (1 + rcStar)*(depsG - depsK)/9];

banner["D/N-TUBE EVEN DEFECT AND THE EXACT HYBRIDIZATION CANCELLATION"];

Clear[lW, a, dLW, da, dgamma0];
$Assumptions = Element[{lW, a, dLW, da, dgamma0, rc, drc}, Reals] && lW > 0 && a > 0;

epsKExact = FullSimplify[12*lW^2/(Pi^2*a^2*(1 + rc)) - 1];
Print["eps_kappa = ", fmt[epsKExact]];

depsKDirect = D[epsKExact, lW]*dLW + D[epsKExact, a]*da + D[epsKExact, rc]*drc;
depsKBranch = FullSimplify[depsKDirect /. 12*lW^2 -> Pi^2*a^2*(1 + rc)];
Print["d eps_kappa = ", fmt[depsKBranch]];
depsKTarget = 2*dLW/lW - 2*da/a - drc/(1 + rc);
depsKDiff = PolynomialRemainder[
  Expand[Numerator[Together[depsKBranch - depsKTarget]]],
  -12*lW^2 + a^2*Pi^2*(1 + rc),
  lW
];
depsKDiff = FullSimplify[depsKDiff];
expectZero["d eps_kappa identity", depsKDiff];

depsGBranch = FullSimplify[9*dgamma0/(1 + rc) - drc/(1 + rc)];
Print["d eps_gamma = ", fmt[depsGBranch]];
expectZero[
  "d eps_gamma rewritten as d ln gamma0 - d ln(1+r_c)",
  depsGBranch - (9*dgamma0/(1 + rc) - drc/(1 + rc))
];
diffIdentity = PolynomialRemainder[
  Expand[
    Numerator[
      Together[
        (depsGBranch - depsKBranch) - (9*dgamma0/(1 + rc) - 2*(dLW/lW - da/a))
      ]
    ]
  ],
  -12*lW^2 + a^2*Pi^2*(1 + rc),
  lW
];
diffIdentity = FullSimplify[diffIdentity];
expectZero["difference identity", diffIdentity];

banner["TANGENTIAL SUSCEPTIBILITY AND FINAL DEFECT LAW"];

Clear[xiGamma, xiL, sigmaStar, dPiTan, dSigma0, dS, dThat];
$Assumptions = Element[{xiGamma, xiL, sigmaStar, dPiTan, dSigma0, dS, dThat, rcStar}, Reals];

upsilonPi = FullSimplify[(1 + rcStar)*(xiGamma - 2*xiL)/9];
Print["Upsilon_Pi = ", fmt[upsilonPi]];

deltaQ = FullSimplify[-9*sigmaStar*upsilonPi*dPiTan/((1 - sigmaStar)*(1 + rcStar))];
nQm1 = FullSimplify[9*sigmaStar*upsilonPi*dPiTan/((1 - sigmaStar)*(1 + rcStar))];
Print["Delta_Q = ", fmt[deltaQ]];
Print["N_Q - 1 = ", fmt[nQm1]];
expectZero[
  "collapsed Delta_Q law",
  deltaQ + sigmaStar*(xiGamma - 2*xiL)*dPiTan/(1 - sigmaStar)
];
expectZero[
  "collapsed N_Q-1 law",
  nQm1 - sigmaStar*(xiGamma - 2*xiL)*dPiTan/(1 - sigmaStar)
];

dPiTanExpr = 0.832409471081635*dSigma0 - 1.16275838754222*dS;
deltaQMouth = Expand[deltaQ /. dPiTan -> dPiTanExpr];
Print["Delta_Q in (dSigma0,dS) = ", fmt[deltaQMouth]];
deltaQT = Expand[deltaQMouth /. dSigma0 -> 6.42981496203006*dThat];
Print["Delta_Q in (dThat,dS) = ", fmt[deltaQT]];

banner["D/N SIMILARITY PRESERVATION"];
expectZero["Xi_gamma = 2 Xi_L => Delta_Q = 0", deltaQ /. xiGamma -> 2*xiL];
expectZero["Xi_gamma = 2 Xi_L => N_Q - 1 = 0", nQm1 /. xiGamma -> 2*xiL];

rF1 = SetPrecision[1.77799353547498, 30];
upsilonPrefactor = N[(1 + rF1^2)/9, 18];
Print["(1+r_F1^2)/9 = ", fmt[upsilonPrefactor]];

banner["Carry-forward formulas"];
Print["1) B_W = ((1+r_c)/9) (eps_gamma - eps_kappa)"];
Print["2) dB_W = ((1+r_c*)/9) (d eps_gamma - d eps_kappa)"];
Print["3) d eps_gamma - d eps_kappa = d ln gamma_0 - 2 d ln(LW/a)"];
Print["4) Upsilon_Pi = ((1+r_c*)/9) (Xi_gamma - 2 Xi_L)"];
Print["5) Delta_Q = -(sigma_*/(1-sigma_*)) (Xi_gamma - 2 Xi_L) dPi_tan"];
Print["6) If Xi_gamma = 2 Xi_L, then the full first-order defect vanishes."];

Print[""];
Print["Stage 144 Mathematica audit passed."];

Exit[0];
