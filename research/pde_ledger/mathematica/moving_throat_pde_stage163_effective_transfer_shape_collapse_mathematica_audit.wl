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

banner["STAGE 163 — EFFECTIVE TRANSFER-SHAPE COLLAPSE"];

Clear[t1, t2, tau1, tau2, eps, lam];
$Assumptions = Element[{t1, t2, tau1, tau2, eps, lam}, Reals] && t1 > 0 && t2 > 0;
teff2 = t1^2*Exp[2*eps*lam*tau1] + t2^2*Exp[2*eps*lam*tau2];
xiEff = FullSimplify[(D[Log[teff2], eps] /. eps -> 0)/lam, Assumptions -> $Assumptions];

rho1 = FullSimplify[t1^2/(t1^2 + t2^2), Assumptions -> $Assumptions];
rho2 = FullSimplify[t2^2/(t1^2 + t2^2), Assumptions -> $Assumptions];
xiExpected = FullSimplify[2*(rho1*tau1 + rho2*tau2), Assumptions -> $Assumptions];
expectZero["multi-port effective-shape identity", xiEff - xiExpected];

banner["One-port continuum transfer shape"];
Clear[muW, muEta, kEta, kW, zW, rho, epsW, omegaW2];
$Assumptions = Element[{muW, muEta, kEta, kW, zW, rho, epsW, omegaW2}, Reals] &&
  muW > 0 && muEta > 0 && kEta > 0 && kW > 0 && zW > 0 && omegaW2 > 0;

k0 = kEta/muEta;
beta0 = (muW/muEta)*(kEta/kW)*zW*(1 + rho)^2/(1 - epsW)^2;
t2Direct = FullSimplify[beta0/k0, Assumptions -> $Assumptions];
t2Expected = FullSimplify[(muW/kW)*zW*(1 + rho)^2/(1 - epsW)^2, Assumptions -> $Assumptions];
expectZero["T^2 = beta0/K0 -> muW/KW form", t2Direct - t2Expected];
expectZero[
  "T^2 = ZW(1+rho)^2 / [OmegaW^2 (1-epsW)^2]",
  t2Direct - ((zW*(1 + rho)^2)/(omegaW2*(1 - epsW)^2) /. omegaW2 -> kW/muW)
];

banner["Selected-branch reformulation"];
Clear[gConst, cs, a, cSpeed, epsEta, rTarget];
$Assumptions = Element[{gConst, cs, a, cSpeed, epsEta, rTarget, muW, kW, zW, rho, epsW}, Reals] &&
  gConst > 0 && cs > 0 && a > 0 && cSpeed > 0 && muW > 0 && kW > 0 && zW > 0;

lambdaNorm = FullSimplify[27*Pi^2*gConst*cs^5*kW/(20*a^5*cSpeed^5*muW), Assumptions -> $Assumptions];
rTargetDef = FullSimplify[lambdaNorm*(1 - epsEta)*(1 - epsW)^2/(zW*(1 + rho)^2), Assumptions -> $Assumptions];
t2Selected = FullSimplify[27*Pi^2*gConst*cs^5*(1 - epsEta)/(20*a^5*cSpeed^5*rTarget), Assumptions -> $Assumptions];
expectZero[
  "selected-branch T^2 identity",
  (t2Direct /. rTarget -> rTargetDef) - (t2Selected /. rTarget -> rTargetDef)
];

banner["Weak-axisymmetric slope laws"];
Clear[zetaW, omegaW, rho1s, epsW1, e, eta1, r1];
$Assumptions = Element[{zetaW, omegaW, rho1s, epsW1, e, eta1, r1, lam, zW, rho, epsW, muW, kW}, Reals] &&
  zW > 0 && muW > 0 && kW > 0;

t2Pert =
  zW*Exp[e*lam*zetaW]*(1 + rho + e*lam*rho1s)^2/
  ((kW/muW)*Exp[e*lam*omegaW]*(1 - epsW - e*lam*epsW1)^2);
xiDirect = FullSimplify[(D[Log[t2Pert], e] /. e -> 0)/lam, Assumptions -> $Assumptions];
xiDirectExpected = FullSimplify[zetaW - omegaW + 2*rho1s/(1 + rho) + 2*epsW1/(1 - epsW), Assumptions -> $Assumptions];
expectZero["direct slope law", xiDirect - xiDirectExpected];

t2SelPert = 27*Pi^2*gConst*cs^5*(1 - epsEta - e*lam*eta1)/(20*a^5*cSpeed^5*rTarget*Exp[e*lam*r1]);
xiSel = FullSimplify[(D[Log[t2SelPert], e] /. e -> 0)/lam, Assumptions -> $Assumptions];
xiSelExpected = FullSimplify[-eta1/(1 - epsEta) - r1, Assumptions -> $Assumptions];
expectZero["selected-branch slope law", xiSel - xiSelExpected];

Print[""];
Print["Carry-forward formulas:"];
Print["  T_eff^2 = sum_r T_r^2 = N_0/K"];
Print["  Xi_1    = delta ln(T_eff^2)/(eps lambda_A)"];
Print["  On the one-port continuum branch,"];
Print["    T^2 = Z_W (1+rho)^2 / [Omega_W^2 (1-eps_W)^2]"];
Print["        = (27 pi^2 G c_s^5 / (20 a^5 c^5)) * (1-eps_eta)/R_target"];
Print["  Hence"];
Print["    Xi_1 = zeta_W - omega_W + 2 rho_1/(1+rho) + 2 epsW_1/(1-eps_W)"];
Print["        = - eta_1/(1-eps_eta) - R_1"];

Print[""];
Print["Stage 163 Mathematica audit passed."];

Exit[0];
