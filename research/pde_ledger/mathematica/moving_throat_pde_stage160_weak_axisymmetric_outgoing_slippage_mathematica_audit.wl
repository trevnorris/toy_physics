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

banner["STAGE 160 — WEAK-AXISYMMETRIC OUTGOING-SLIPPAGE COLLAPSE"];

Clear[eps, lam, k, ou2, ow2, r, gu, gw, kappa1, gW1, gU1, r1, oU1, oW1];
$Assumptions = Element[{eps, lam, k, ou2, ow2, r, gu, gw, kappa1, gW1, gU1, r1, oU1, oW1}, Reals] &&
  k > 0 && ou2 > 0 && ow2 > 0 && r > 0 && gu > 0 && gw > 0;

mCal = gw/(Sqrt[k]*ow2);
iCal = r*gu/(ou2*gw);
hCal = r^2/(ou2*ow2);
lambda0 = (ou2*gw + r*gu)/(ou2*ow2 - r^2);

kP = k*Exp[eps*lam*kappa1];
gwP = gw*Exp[eps*lam*gW1];
guP = gu*Exp[eps*lam*gU1];
rP = r*Exp[eps*lam*r1];
ou2P = ou2*Exp[eps*lam*oU1];
ow2P = ow2*Exp[eps*lam*oW1];

mCalP = gwP/(Sqrt[kP]*ow2P);
iCalP = rP*guP/(ou2P*gwP);
hCalP = rP^2/(ou2P*ow2P);
lambdaP = (ou2P*gwP + rP*guP)/(ou2P*ow2P - rP^2);

mR = gW1 - oW1 - kappa1/2;
iR = r1 + gU1 - oU1 - gW1;
hR = 2*r1 - oU1 - oW1;

dlnMExact = FullSimplify[D[Log[mCalP/mCal], eps] /. eps -> 0, Assumptions -> $Assumptions];
dlnIExact = FullSimplify[D[Log[iCalP/iCal], eps] /. eps -> 0, Assumptions -> $Assumptions];
dlnHExact = FullSimplify[D[Log[hCalP/hCal], eps] /. eps -> 0, Assumptions -> $Assumptions];

expectZero["weak-axisymmetric d ln M", dlnMExact - lam*mR];
expectZero["weak-axisymmetric d ln I", dlnIExact - lam*iR];
expectZero["weak-axisymmetric d ln H", dlnHExact - lam*hR];

banner["Portwise outgoing-defect amplitude"];
sigmaExact = FullSimplify[D[Log[(lambdaP^2/kP)/(lambda0^2/k)], eps] /. eps -> 0, Assumptions -> $Assumptions];
sigmaR = FullSimplify[2*mR + 2*iCal*iR/(1 + iCal) + 2*hCal*hR/(1 - hCal), Assumptions -> $Assumptions];
expectZero["Sigma_{A,r} = lambda_A sigma_r", sigmaExact - lam*sigmaR];

banner["Grouped weak-axisymmetric trace/anomaly law"];
Clear[xi1, epsilon];
$Assumptions = Element[{xi1, epsilon}, Reals];
lam20 = 1;
lam21 = 1/2;
lam22 = -1;
xi20 = epsilon*lam20*xi1;
xi21 = epsilon*lam21*xi1;
xi22 = epsilon*lam22*xi1;

xiBar = FullSimplify[(xi20 + 2*xi21 + 2*xi22)/5, Assumptions -> $Assumptions];
aXi = FullSimplify[(2*xi20 - xi21 - xi22)/10, Assumptions -> $Assumptions];
bXi = FullSimplify[(xi21 - xi22)/2, Assumptions -> $Assumptions];

expectZero["grouped trace vanishes", xiBar];
expectZero["a_Xi - eps Xi1/4", aXi - epsilon*xi1/4];
expectZero["b_Xi - 3 eps Xi1/4", bXi - 3*epsilon*xi1/4];
expectZero["axisymmetric defect law b_Xi - 3 a_Xi", bXi - 3*aXi];

banner["Outgoing-prefactor slope identification"];
Clear[p0, d0, n0, d1, n1];
$Assumptions = Element[{eps, lam, p0, d0, n0, d1, n1}, Reals] && p0 > 0 && d0 > 0 && n0 > 0;
p0Exact = n0/d0;
dA = d0*Exp[eps*lam*d1];
nA = n0*Exp[eps*lam*n1];
pA = FullSimplify[nA/dA, Assumptions -> $Assumptions];
pSlopeExact = FullSimplify[D[pA, eps] /. eps -> 0, Assumptions -> $Assumptions]/lam;
xiLoad = FullSimplify[n1 - d1, Assumptions -> $Assumptions];
expectZero["P_A slope = P0 * Xi_load", pSlopeExact - p0Exact*xiLoad];
expectZero["(P1/P0) - Xi_load", pSlopeExact/p0Exact - xiLoad];

Print[""];
Print["Carry-forward formulas:"];
Print["  m_r = gW1 - oW1 - kappa1/2"];
Print["  i_r = r1 + gU1 - oU1 - gW1"];
Print["  h_r = 2 r1 - oU1 - oW1"];
Print["  sigma_r = 2 m_r + 2 I/(1+I) i_r + 2 H/(1-H) h_r"];
Print["  Xi_A = eps lambda_A Xi1  =>  abar=0, a=eps Xi1/4, b=3 eps Xi1/4"];
Print["  On the conservative-shape-preserving branch, Xi1 = P1/P0."];

Print[""];
Print["Stage 160 Mathematica audit passed."];

Exit[0];
