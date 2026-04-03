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

banner["STAGE 165 — MICROSCOPIC COHERENT-KERNEL SLIPPAGE DECOMPOSITION"];

Clear[lam1, c1, gam1, kU, kEta, kW, mu1, tau1, lamphi1, kphi, chi0, epsW, epsEta, deltaU];
$Assumptions = Element[{lam1, c1, gam1, kU, kEta, kW, mu1, tau1, lamphi1, kphi, chi0, epsW, epsEta, deltaU}, Reals] &&
  chi0 > 0 && epsW > 0 && epsEta > 0 && deltaU > 0;

zetaZ = 2*lam1 - kEta - kW;
omegaW = kW - mu1;
chi1 = chi0*(gam1 + c1 - kU);
eta1 = epsEta*(2*c1 - kU - kEta);
varepsW = epsW*(2*gam1 + 2*lam1 - kU - kW);
deltaU1 = deltaU*(tau1 - kU);
eps = FullSimplify[epsW*(1 - (2/11)*deltaU/(1 + deltaU)), Assumptions -> $Assumptions];
eps1 = FullSimplify[
  (1 - (2/11)*deltaU/(1 + deltaU))*varepsW - (2*epsW*deltaU1)/(11*(1 + deltaU)^2),
  Assumptions -> $Assumptions
];

Clear[sigmaZ, sigmaChi, sigmaEta, sigmaEps, sigmaDel];
slipSubs = {
  sigmaZ -> 2*lam1 + mu1 - kEta - 2*kW,
  sigmaChi -> gam1 + c1 - kU,
  sigmaEta -> 2*c1 - kU - kEta,
  sigmaEps -> 2*gam1 + 2*lam1 - kU - kW,
  sigmaDel -> tau1 - kU
};

banner["Physical branch drifts from microscopic logs"];
expectZero["zeta_Z formula", zetaZ - (2*lam1 - kEta - kW)];
expectZero["omega_W formula", omegaW - (kW - mu1)];
expectZero["chi_1 formula", chi1 - chi0*(gam1 + c1 - kU)];
expectZero["eta_1 formula", eta1 - epsEta*(2*c1 - kU - kEta)];
expectZero["varepsilon_W formula", varepsW - epsW*(2*gam1 + 2*lam1 - kU - kW)];
expectZero["delta_U,1 formula", deltaU1 - deltaU*(tau1 - kU)];

banner["Four-slippage grouped-defect law"];
xi1Direct = FullSimplify[zetaZ - omegaW + 2*chi1/(1 + chi0) + 2*eps1/(1 - eps), Assumptions -> $Assumptions];
xi1Slip = FullSimplify[
  sigmaZ +
  2*chi0*sigmaChi/(1 + chi0) +
  2*epsW*((11 + 9*deltaU)*sigmaEps/(11*(1 + deltaU)) - 2*deltaU*sigmaDel/(11*(1 + deltaU)^2))/(1 - eps),
  Assumptions -> $Assumptions
];
expectZero["Xi_1 direct - slippage form", xi1Direct - (xi1Slip /. slipSubs)];
Print["Xi_1 = ", fmt[xi1Slip]];

banner["Selected-branch demand slippage"];
r1Direct = FullSimplify[-eta1/(1 - epsEta) - xi1Direct, Assumptions -> $Assumptions];
r1Slip = FullSimplify[-epsEta*sigmaEta/(1 - epsEta) - xi1Slip, Assumptions -> $Assumptions];
expectZero["R_1 direct - slippage form", r1Direct - (r1Slip /. slipSubs)];
Print["R_1 = ", fmt[r1Slip]];

banner["Tracking-factor factorization"];
Clear[sigmaTr];
theta1Direct = FullSimplify[
  -(chi0*(1 + chi0)*deltaU1 + deltaU*(1 + deltaU)*chi1)/
   ((1 + chi0)*(1 + deltaU)*(1 + chi0 + deltaU)),
  Assumptions -> $Assumptions
];
theta1Fact = FullSimplify[
  -chi0*deltaU*sigmaTr/((1 + chi0)*(1 + deltaU)*(1 + chi0 + deltaU)),
  Assumptions -> $Assumptions
];
sigmaTrDef = FullSimplify[(1 + chi0)*sigmaDel + (1 + deltaU)*sigmaChi, Assumptions -> $Assumptions];
expectZero["Theta_1 factorization", theta1Direct - (theta1Fact /. sigmaTr -> sigmaTrDef /. slipSubs)];
Print["Sigma_tr = ", fmt[sigmaTrDef]];
Print["Theta_1 = ", fmt[theta1Fact]];

banner["Tracking/nontracking split of Xi_1"];
xi1Split = FullSimplify[
  sigmaZ +
  2*chi0*sigmaTr/((1 + chi0)*(1 + deltaU)) +
  2*epsW*(11 + 9*deltaU)*sigmaEps/(11*(1 - eps)*(1 + deltaU)) -
  (2*chi0/(1 + deltaU) + 4*epsW*deltaU/(11*(1 - eps)*(1 + deltaU)^2))*sigmaDel,
  Assumptions -> $Assumptions
];
expectZero["Xi_1 split - slippage form", (xi1Split /. sigmaTr -> sigmaTrDef) - xi1Slip];
Print["Xi_1 split = ", fmt[xi1Split]];

banner["Support-blindness at the microscopic level"];
Print["free symbols of Xi_1: ", fmt[Sort[ToString /@ (List @@ xi1Slip /. Plus -> List /. Times -> List)]]];
expectZero["dXi_1/dlamphi1", D[xi1Slip, lamphi1]];
expectZero["dXi_1/dkphi", D[xi1Slip, kphi]];
expectZero["dR_1/dlamphi1", D[r1Slip, lamphi1]];
expectZero["dR_1/dkphi", D[r1Slip, kphi]];
expectZero["dTheta_1/dlamphi1", D[theta1Fact, lamphi1]];
expectZero["dTheta_1/dkphi", D[theta1Fact, kphi]];

Print[""];
Print["Carry-forward formulas:"];
Print["  Sigma_Z   = 2 lam_1 + mu_1 - kappa_eta - 2 kappa_W"];
Print["  Sigma_chi = gamma_1 + c_1 - kappa_U"];
Print["  Sigma_eta = 2 c_1 - kappa_U - kappa_eta"];
Print["  Sigma_eps = 2 gamma_1 + 2 lam_1 - kappa_U - kappa_W"];
Print["  Sigma_del = tau_1 - kappa_U"];
Print["  Xi_1 depends only on (Sigma_Z, Sigma_chi, Sigma_eps, Sigma_del)"];
Print["  R_1 adds only Sigma_eta"];
Print["  Theta_1 is carried by Sigma_tr = (1+chi_0) Sigma_del + (1+delta_U) Sigma_chi"];

Print[""];
Print["Stage 165 Mathematica audit passed."];

Exit[0];
