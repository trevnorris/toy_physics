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

banner["STAGE 179 — WALL-NORMALIZED TRANSFER-SHAPE THEOREM"];

Clear[k, ou2, ow2, gw, gu, r];
$Assumptions = Element[{k, ou2, ow2, gw, gu, r}, Reals] &&
  k > 0 && ou2 > 0 && ow2 > 0 && gw > 0 && gu > 0;

p = ou2*gw + r*gu;
delta = ou2*ow2 - r^2;
n0 = FullSimplify[p^2/delta^2, Assumptions -> $Assumptions];

gWh = gw/(Sqrt[k]*ow2);
gUh = gu/(Sqrt[k]*Sqrt[ou2]*Sqrt[ow2]);
rHat = r/(Sqrt[ou2]*Sqrt[ow2]);
t = FullSimplify[(gWh + rHat*gUh)/(1 - rHat^2), Assumptions -> $Assumptions];

Print["P   = ", fmt[p]];
Print["Delta = ", fmt[delta]];
Print["T   = ", fmt[t]];
expectZero["N0/K - T^2", n0/k - t^2];

banner["Weak-axisymmetric slope identity"];
Clear[eps, lam, kappa1, gW, gU, oU, oW, rr];
$Assumptions = Element[{eps, lam, kappa1, gW, gU, oU, oW, rr, k, ou2, ow2, gw, gu, r}, Reals] &&
  k > 0 && ou2 > 0 && ow2 > 0 && gw > 0 && gu > 0;

kA = k*(1 + eps*lam*kappa1);
ou2A = ou2*(1 + eps*lam*oU);
ow2A = ow2*(1 + eps*lam*oW);
gwA = gw*(1 + eps*lam*gW);
guA = gu*(1 + eps*lam*gU);
rA = r*(1 + eps*lam*rr);

pA = Expand[ou2A*gwA + rA*guA];
deltaA = Expand[ou2A*ow2A - rA^2];
n0A = Expand[pA^2/deltaA^2];
nuDirect = FullSimplify[(D[Log[n0A], eps] /. eps -> 0)/lam, Assumptions -> $Assumptions];
Print["nu_direct = ", fmt[nuDirect]];

w = FullSimplify[gW - oW - kappa1/2, Assumptions -> $Assumptions];
u = FullSimplify[gU - oU/2 - oW/2 - kappa1/2, Assumptions -> $Assumptions];
c = FullSimplify[rr - oU/2 - oW/2, Assumptions -> $Assumptions];

alpha = FullSimplify[gWh/(gWh + rHat*gUh), Assumptions -> $Assumptions];
beta = FullSimplify[rHat*gUh/(gWh + rHat*gUh), Assumptions -> $Assumptions];
tau = FullSimplify[alpha*w + beta*(u + c) + 2*rHat^2*c/(1 - rHat^2), Assumptions -> $Assumptions];
nuExpected = FullSimplify[kappa1 + 2*tau, Assumptions -> $Assumptions];

Print["tau = ", fmt[tau]];
Print["nu_expected = ", fmt[nuExpected]];
expectZero["nu_direct - (kappa1 + 2 tau)", nuDirect - nuExpected];

banner["Exact equivalence to the Stage 176/160/161 slippage formulas"];
iHat = FullSimplify[rHat*gUh/gWh, Assumptions -> $Assumptions];
hHat = FullSimplify[rHat^2, Assumptions -> $Assumptions];
m = w;
i = FullSimplify[(u + c) - w, Assumptions -> $Assumptions];
h = FullSimplify[2*c, Assumptions -> $Assumptions];
tauSlippage = FullSimplify[m + iHat*i/(1 + iHat) + hHat*h/(1 - hHat), Assumptions -> $Assumptions];
expectZero["tau - slippage form", tau - tauSlippage];
expectZero["(nu-kappa1) - 2*tau_slippage", (nuExpected - kappa1) - 2*tauSlippage];

banner["Weighted defect identity"];
Clear[rho1, rho2, tau1, tau2, tau3];
$Assumptions = Element[{rho1, rho2, tau1, tau2, tau3, kappa1}, Reals];
rho3 = 1 - rho1 - rho2;
xi = rho1*(kappa1 + 2*tau1) + rho2*(kappa1 + 2*tau2) + rho3*(kappa1 + 2*tau3) - kappa1;
xiExpected = 2*(rho1*tau1 + rho2*tau2 + rho3*tau3);
expectZero["Xi_1 - 2 weighted tau", xi - xiExpected];

Print[""];
Print["Carry-forward formulas:"];
Print["  N0^(r) = K * T_r^2"];
Print["  T_r = (Ghat_W + Rhat Ghat_U)/(1-Rhat^2)"];
Print["  nu_r = kappa1 + 2 tau_r"];
Print["  tau_r = m_r + I_r/(1+I_r) i_r + H_r/(1-H_r) h_r"];
Print["  Xi_1 = 2 sum_r rho_r^(N) tau_r"];

Print[""];
Print["Stage 179 Mathematica audit passed."];

Exit[0];
