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

Clear[rho, capitalK, m];
$Assumptions = Element[{rho, capitalK, m}, Reals] && rho > 0 && capitalK > 0 && m > 0;

uRho = capitalK*rho^5/4;
hRho = FullSimplify[D[uRho, rho], Assumptions -> $Assumptions];
hPrime = FullSimplify[D[hRho, rho], Assumptions -> $Assumptions];
csSq = FullSimplify[(1/m)*D[capitalK*rho^5, rho], Assumptions -> $Assumptions];

Print["U(rho) = ", fmt[uRho]];
Print["h(rho) = ", fmt[hRho]];
Print["h'(rho) = ", fmt[hPrime]];
Print["c_s^2(rho) = ", fmt[csSq]];
expectZero["h'(rho) - m c_s^2 / rho", hPrime - m*csSq/rho];

Clear[rhoStar, csStarSq, nSS, nPP, oSP, gPhi, kX, tX, ell, kappa, c2];
$Assumptions =
  Element[{rhoStar, csStarSq, nSS, nPP, oSP, gPhi, kX, tX, ell, kappa, c2, m}, Reals] &&
  rhoStar > 0 && csStarSq > 0 && nSS > 0 && nPP > 0 && oSP > 0 &&
  gPhi > 0 && kX > 0 && tX > 0 && ell > 0 && kappa > 0 && c2 >= 0 && m > 0;

thetaSigma = (m*csStarSq/rhoStar)*nSS;
lambdaPhi = gPhi*oSP;
chiEff = FullSimplify[1/thetaSigma, Assumptions -> $Assumptions];
gMicro = FullSimplify[chiEff*lambdaPhi^2/kX, Assumptions -> $Assumptions];
gExpected = FullSimplify[rhoStar*gPhi^2*oSP^2/(m*csStarSq*kX*nSS), Assumptions -> $Assumptions];

Print["Theta_sigma = ", fmt[thetaSigma]];
Print["Lambda_phi = ", fmt[lambdaPhi]];
Print["chi_sigma^(eff) = ", fmt[chiEff]];
Print["G_micro = ", fmt[gMicro]];
expectZero["G_micro - expected parent formula", gMicro - gExpected];

c2Def = FullSimplify[oSP^2/(nSS*nPP), Assumptions -> $Assumptions];
gCoherence = FullSimplify[rhoStar*gPhi^2*nPP*c2/(m*csStarSq*kX), Assumptions -> $Assumptions];
Print["C_(sigma phi)^2 definition = ", fmt[c2Def]];
expectZero["coherence-factor decomposition", FullSimplify[gExpected /. oSP^2 -> c2*nSS*nPP, Assumptions -> $Assumptions] - gCoherence];

xiMicro = FullSimplify[kappa*gMicro, Assumptions -> $Assumptions];
xiExpected = FullSimplify[rhoStar*gPhi^2*oSP^2*ell^2/(m*csStarSq*tX*nSS), Assumptions -> $Assumptions];
Print["Xi_micro = ", fmt[xiMicro]];
expectZero["Xi_micro - parent projected formula", FullSimplify[xiMicro /. kappa -> kX*ell^2/tX, Assumptions -> $Assumptions] - xiExpected];

Print[""];
Print["Stage 062 Mathematica audit passed."];

Exit[0];
