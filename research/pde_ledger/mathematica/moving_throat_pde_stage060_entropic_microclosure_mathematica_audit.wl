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

banner["STAGE 043 — ENTROPIC MICROCLOSURE"];

Clear[s, ell, theta, sigmaStar, lambdaPhi, tX, kX, kM, mSigma, dSigma, chiSigma, kappa, eta, deltaDrop, phi0, pe];
$Assumptions =
  Element[{s, ell, theta, sigmaStar, lambdaPhi, tX, kX, kM, mSigma, dSigma, chiSigma, kappa, eta, deltaDrop, phi0, pe}, Reals] &&
  ell > 0 && theta > 0 && sigmaStar > 0 && lambdaPhi > 0 && tX > 0 &&
  kX > 0 && kM > 0 && mSigma > 0 && dSigma > 0 && chiSigma > 0 &&
  kappa > 0 && eta > 0 && pe > 0;

sigma = Symbol["sigma"];
phi = Symbol["phi"];
fSigma = theta*sigma*(Log[sigma/sigmaStar] - 1) - lambdaPhi*sigma*phi;
muExpr = FullSimplify[D[fSigma, sigma], Assumptions -> $Assumptions];

Print["f_sigma = ", fmt[fSigma]];
Print["mu_sigma^(chem) = ", fmt[muExpr]];
expectZero["mu - [Theta log(sigma/sigma_*) - Lambda phi]", muExpr - (theta*Log[sigma/sigmaStar] - lambdaPhi*phi)];

sigmaField = sigma[s];
phiField = phi[s];
jExpr = Expand[-mSigma*sigmaField*D[(theta*Log[sigmaField/sigmaStar] - lambdaPhi*phiField), s]];
Print["J expanded = ", fmt[jExpr]];
expectZero[
  "Onsager current decomposition",
  jExpr - (-mSigma*theta*D[sigmaField, s] + mSigma*lambdaPhi*sigmaField*D[phiField, s])
];
expectZero[
  "current with D_sigma substitution",
  FullSimplify[jExpr /. (mSigma*theta) -> dSigma, Assumptions -> $Assumptions] -
    (-dSigma*D[sigmaField, s] + mSigma*lambdaPhi*sigmaField*D[phiField, s])
];

a = FullSimplify[lambdaPhi*deltaDrop/(theta*ell), Assumptions -> $Assumptions];
cNorm = Symbol["Cnorm"];
sigmaTrial = cNorm*Exp[a*s];
phiAff = phi0 + deltaDrop*s/ell;
jAff = FullSimplify[
  jExpr /. {
    phiField -> phiAff,
    Derivative[1][phi][s] -> deltaDrop/ell
  },
  Assumptions -> $Assumptions
];
expectZero[
  "J=0 solved by exponential family",
  FullSimplify[
    jAff /. {
      sigmaField -> sigmaTrial,
      Derivative[1][sigma][s] -> D[sigmaTrial, s]
    },
    Assumptions -> $Assumptions
  ]
];

sigmaPe = FullSimplify[pe*Exp[pe*s]/(Exp[pe] - 1), Assumptions -> pe > 0];
expectZero["normalized Sigma_Pe family", Integrate[sigmaPe, {s, 0, 1}] - 1];
expectZero["Pe identification", a*ell - lambdaPhi*deltaDrop/theta];

deltaSupport = Symbol["DeltaSupport"];
xiMicro = FullSimplify[(lambdaPhi/ theta)*(lambdaPhi*ell^2*deltaSupport/tX)/deltaSupport, Assumptions -> $Assumptions];
Print["Xi_micro = ", fmt[xiMicro]];
expectZero["Xi_micro - Lambda^2 L^2/(Theta T_X)", xiMicro - lambdaPhi^2*ell^2/(theta*tX)];
expectZero["Xi_micro susceptibility form", FullSimplify[xiMicro /. theta -> 1/chiSigma, Assumptions -> $Assumptions] - chiSigma*lambdaPhi^2*ell^2/tX];
expectZero["Xi_micro phenomenological form", FullSimplify[xiMicro /. theta -> dSigma/mSigma, Assumptions -> $Assumptions] - mSigma*lambdaPhi^2*ell^2/(dSigma*tX)];

muFun = mu[s];
jFun = j[s];
localIdentity = FullSimplify[-muFun*D[jFun, s] - (-D[muFun*jFun, s] + D[muFun, s]*jFun), Assumptions -> $Assumptions];
expectZero["integration-by-parts identity", localIdentity];
muS = Symbol["muS"];
jOnsager = -mSigma*sigma*muS;
expectZero["Onsager dissipation density", FullSimplify[muS*jOnsager + jOnsager^2/(mSigma*sigma), Assumptions -> $Assumptions]];

phiHead = Unique["phiFun"];
sigmaHead = Unique["sigmaFun"];
phiFun = phiHead[s];
sigmaFun = sigmaHead[s];
fFull = theta*sigmaFun*(Log[sigmaFun/sigmaStar] - 1) - lambdaPhi*sigmaFun*phiFun + 1/2*tX*D[phiFun, s]^2 + 1/2*kX*phiFun^2;
elPhi = FullSimplify[D[fFull, phiFun] - D[D[fFull, D[phiFun, s]], s], Assumptions -> $Assumptions];
Print["Euler-Lagrange support bulk equation = ", fmt[elPhi]];
expectZero["support bulk equation form", elPhi - (-lambdaPhi*sigmaFun + kX*phiFun - tX*D[phiFun, {s, 2}])];

Print[""];
Print["Stage 060 Mathematica audit passed."];

Exit[0];
