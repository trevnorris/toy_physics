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

banner["STAGE 161 — OUTGOING-PORT CO-LOADING THEOREM"];

banner["1. Exact weak-axisymmetric slope of N_{A,0}^{(r)}"];
Clear[eps, lam, p0r, d0r, pR, dR];
$Assumptions = Element[{eps, lam, p0r, d0r, pR, dR}, Reals] && p0r > 0 && d0r > 0;

pAr = p0r*(1 + eps*lam*pR);
dAr = d0r*(1 + eps*lam*dR);
n0r = p0r^2/d0r^2;
nAr = Expand[pAr^2/dAr^2];
nuFromSeries = FullSimplify[(Normal[Series[nAr, {eps, 0, 1}]]/n0r - 1)/(eps*lam), Assumptions -> $Assumptions];
expectZero["nu_r - 2(p_r-d_r)", nuFromSeries - 2*(pR - dR)];

Print["N_{A,0}^{(r)} = ", fmt[Normal[Series[nAr, {eps, 0, 1}]]]];
Print["nu_r = ", fmt[FullSimplify[2*(pR - dR), Assumptions -> $Assumptions]]];

banner["2. Weighted collapse Xi_1 = nubar_N - kappa_1"];
Clear[rho1, rho2, rho3, nu1, nu2, nu3, kappa1];
$Assumptions = Element[{rho1, rho2, rho3, nu1, nu2, nu3, kappa1}, Reals];
xi = Expand[rho1*(nu1 - kappa1) + rho2*(nu2 - kappa1) + rho3*(nu3 - kappa1)];
nuBar = Expand[rho1*nu1 + rho2*nu2 + rho3*nu3];
expectZero[
  "Xi_1 - (nubar_N - kappa_1)",
  (xi /. rho3 -> 1 - rho1 - rho2) - ((nuBar /. rho3 -> 1 - rho1 - rho2) - kappa1)
];
Print["Xi_1 = ", fmt[xi]];
Print["nubar_N = ", fmt[nuBar]];

banner["3. Exact formulas for p_r and d_r from actual port data"];
Clear[ou2, ow2, r, gu, gw, oU, oW, rr, gU, gW];
$Assumptions = Element[{eps, lam, ou2, ow2, r, gu, gw, oU, oW, rr, gU, gW}, Reals] &&
  ou2 > 0 && ow2 > 0 && r > 0 && gu > 0 && gw > 0;

p = ou2*gw + r*gu;
delta = ou2*ow2 - r^2;

pA = ou2*gw*(1 + eps*lam*(oU + gW)) + r*gu*(1 + eps*lam*(rr + gU));
dA = ou2*ow2*(1 + eps*lam*(oU + oW)) - r^2*(1 + 2*eps*lam*rr);

pFromSeries = FullSimplify[(Normal[Series[pA/p, {eps, 0, 1}]] - 1)/(eps*lam), Assumptions -> $Assumptions];
dFromSeries = FullSimplify[(Normal[Series[dA/delta, {eps, 0, 1}]] - 1)/(eps*lam), Assumptions -> $Assumptions];

alpha = FullSimplify[ou2*gw/p, Assumptions -> $Assumptions];
beta = FullSimplify[r*gu/p, Assumptions -> $Assumptions];
chi = FullSimplify[ou2*ow2/delta, Assumptions -> $Assumptions];
zeta = FullSimplify[r^2/delta, Assumptions -> $Assumptions];

pExpected = FullSimplify[alpha*(oU + gW) + beta*(rr + gU), Assumptions -> $Assumptions];
dExpected = FullSimplify[chi*(oU + oW) - 2*zeta*rr, Assumptions -> $Assumptions];

expectZero["p_r formula", pFromSeries - pExpected];
expectZero["d_r formula", dFromSeries - dExpected];
expectZero["alpha+beta-1", alpha + beta - 1];
expectZero["chi-zeta-1", chi - zeta - 1];

Print["alpha_r = ", fmt[alpha]];
Print["beta_r  = ", fmt[beta]];
Print["chi_r   = ", fmt[chi]];
Print["zeta_r  = ", fmt[zeta]];
Print["p_r     = ", fmt[pExpected]];
Print["d_r     = ", fmt[dExpected]];

banner["4. Equivalence with the Stage-160 slippage formula"];
iRExpr = FullSimplify[r*gu/(ou2*gw), Assumptions -> $Assumptions];
hRExpr = FullSimplify[r^2/(ou2*ow2), Assumptions -> $Assumptions];

mR = gW - oW - kappa1/2;
iR = rr + gU - oU - gW;
hR = 2*rr - oU - oW;

nuDirect = FullSimplify[2*pExpected - 2*dExpected, Assumptions -> $Assumptions];
nuExpected = FullSimplify[
  kappa1 + 2*mR + 2*iRExpr*iR/(1 + iRExpr) + 2*hRExpr*hR/(1 - hRExpr),
  Assumptions -> $Assumptions
];
expectZero["nu_r - [kappa1 + sigma_r]", nuDirect - nuExpected];

sigmaR = FullSimplify[nuExpected - kappa1, Assumptions -> $Assumptions];
Print["Ir = ", fmt[iRExpr]];
Print["Hr = ", fmt[hRExpr]];
Print["nu_r = ", fmt[nuExpected]];
Print["sigma_r = ", fmt[sigmaR]];
Print["Hence Xi_1 = sum_r rho_r^(N) sigma_r = sum_r rho_r^(N)(nu_r-kappa1)."];

Print[""];
Print["Stage 161 Mathematica audit passed."];

Exit[0];
