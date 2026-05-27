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

banner["STAGE 045 — PARENT ACTION GAIN"];

Clear[rho, capitalK, m, nPoly];
$Assumptions =
  Element[{rho, capitalK, m, nPoly}, Reals] &&
  Element[nPoly, Integers] && rho > 0 && capitalK > 0 && m > 0 && nPoly > 1;

(* Preserve the stage normalization U = K rho^n/(n - 1), with c_s^2 = d(K rho^n)/d rho / m. *)
uGeneral = capitalK*rho^nPoly/(nPoly - 1);
hGeneral = FullSimplify[D[uGeneral, rho], Assumptions -> $Assumptions];
hPrimeGeneral = FullSimplify[D[hGeneral, rho], Assumptions -> $Assumptions];
csSqGeneral = FullSimplify[(1/m)*D[capitalK*rho^nPoly, rho], Assumptions -> $Assumptions];

expectZero["h'(rho) = m c_s^2 / rho (general polytrope)", hPrimeGeneral - m*csSqGeneral/rho];

subsN5 = {nPoly -> 5};
Print["Specializing to n=5:"];
Print["  U(rho) = ", fmt[FullSimplify[uGeneral /. subsN5, Assumptions -> $Assumptions]]];
Print["  h(rho) = ", fmt[FullSimplify[hGeneral /. subsN5, Assumptions -> $Assumptions]]];
Print["  h'(rho) = ", fmt[FullSimplify[hPrimeGeneral /. subsN5, Assumptions -> $Assumptions]]];
Print["  c_s^2(rho) = ", fmt[FullSimplify[csSqGeneral /. subsN5, Assumptions -> $Assumptions]]];
expectZero["n=5 specialization of h' = m c_s^2/rho", (hPrimeGeneral - m*csSqGeneral/rho) /. subsN5];

csSqWrong = FullSimplify[(1/m)*D[capitalK*rho^(nPoly + 1), rho], Assumptions -> $Assumptions];
residualWrong = FullSimplify[hPrimeGeneral - m*csSqWrong/rho, Assumptions -> $Assumptions];
residualWrongN5 = FullSimplify[residualWrong /. subsN5, Assumptions -> $Assumptions];
If[FullSimplify[residualWrongN5, Assumptions -> $Assumptions] === 0,
  fail["Inconsistency check failed to detect wrong exponent", residualWrongN5],
  Print["inconsistency probe nonzero: ", fmt[residualWrongN5]]
];

Clear[rhoStar, csStarSq, nSS, nPP, oSP, gPhi, kX, tX, ell, kappa, sigma, phi];
$Assumptions =
  Element[{rhoStar, csStarSq, nSS, nPP, oSP, gPhi, kX, tX, ell, kappa, sigma, phi, m}, Reals] &&
  rhoStar > 0 && csStarSq > 0 && nSS > 0 && nPP > 0 && oSP > 0 &&
  gPhi > 0 && kX > 0 && tX > 0 && ell > 0 && kappa > 0 && m > 0;

thetaSigma = (m*csStarSq/rhoStar)*nSS;
lambdaPhi = gPhi*oSP;

(* Independent susceptibility-route derivation (notes section 4):
   chi_sigma^(eff) = 1/Theta_sigma, then G_micro = chi_sigma^(eff) * Lambda_phi^2 / K_X *)
chiSigmaEff = 1/thetaSigma;
gainViaSusceptibility = FullSimplify[chiSigmaEff * lambdaPhi^2 / kX, Assumptions -> $Assumptions];

(* Action-coefficient route. *)
sParent = (1/2)*thetaSigma*sigma^2 - lambdaPhi*sigma*phi + (1/2)*kX*phi^2;
sigmaStar = sigma /. First[Solve[D[sParent, sigma] == 0, sigma]];
sigmaStar = sigmaStar /. ConditionalExpression[e_, _] :> e;
sEff = Expand[sParent /. sigma -> sigmaStar];
phiQuadraticCoeff = FullSimplify[Coefficient[sEff, phi, 2], Assumptions -> $Assumptions];
gainFromAction = FullSimplify[(kX - 2*phiQuadraticCoeff)/kX, Assumptions -> $Assumptions];
gainFromSeries = FullSimplify[(kX - 2*SeriesCoefficient[Series[sEff, {phi, 0, 2}], 2])/kX, Assumptions -> $Assumptions];
gClosed = FullSimplify[rhoStar*gPhi^2*oSP^2/(m*csStarSq*kX*nSS), Assumptions -> $Assumptions];

Print["Theta_sigma = ", fmt[thetaSigma]];
Print["Lambda_phi = ", fmt[lambdaPhi]];
Print["chi_sigma^(eff) = ", fmt[chiSigmaEff]];
Print["G_micro via susceptibility route = ", fmt[gainViaSusceptibility]];
Print["sigmaStar = ", fmt[sigmaStar]];
Print["S_eff = ", fmt[sEff]];
Print["G_micro from action = ", fmt[gainFromAction]];
expectZero["G_micro via susceptibility route vs closed form", gainViaSusceptibility - gClosed];
expectZero["G_micro: action route equals susceptibility route", gainFromAction - gainViaSusceptibility];
expectZero["Mathematica two-route consistency", gainFromAction - gainFromSeries];
expectZero["gMicro from parent action vs closed form", gainFromAction - gClosed];

(* Second equality of boxed eq:app-stage062-Gmicro *)
cSpSq = oSP^2 / (nSS * nPP);
gMicroFactored = (rhoStar * gPhi^2 * nPP / (m * csStarSq * kX)) * cSpSq;
expectZero["Second equality of boxed G_micro: closed vs factored form", gClosed - gMicroFactored];

(* Cauchy-Schwarz parameterization: O_{sigma phi} = cos(theta) sqrt(N_ss N_pp) *)
theta = Symbol["theta"];
cSpSqCos = cSpSq /. oSP -> Cos[theta] * Sqrt[nSS * nPP];
expectZero["C_{sigma phi}^2 via Cauchy parameterization equals cos^2(theta)",
           FullSimplify[cSpSqCos - Cos[theta]^2, Assumptions -> $Assumptions]];
Print["C_sp_sq Cauchy parameterization yields Cos[theta]^2 in [0, 1] (Cauchy-Schwarz bound)."];

Print["Coherence factor (definition):  C_(sigma phi)^2 := O_sp^2 / (N_ss N_pp)"];

xiMicro = FullSimplify[kappa*gainFromAction, Assumptions -> $Assumptions];
xiTarget = FullSimplify[rhoStar*gPhi^2*oSP^2*ell^2/(m*csStarSq*tX*nSS), Assumptions -> $Assumptions];
Print["Xi_micro = ", fmt[xiMicro]];
kappaRules = Solve[xiMicro == xiTarget, kappa];
kappaSolved = FullSimplify[(kappa /. First[kappaRules]) /. ConditionalExpression[e_, _] :> e, Assumptions -> $Assumptions];
Print["kappa solved from xiMicro = xiTarget: ", fmt[kappaSolved]];
If[FullSimplify[kappaSolved == kX*ell^2/tX, Assumptions -> $Assumptions] === True,
  pass["kappa solution equals kX ell^2/tX"],
  fail["Unexpected kappa solution", kappaRules]
];

Print[""];
Print["Stage 062 Mathematica audit passed."];

Exit[0];
