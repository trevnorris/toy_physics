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

banner["STAGE 182 — MICROSCOPIC COHERENT-KERNEL SLIPPAGE DECOMPOSITION"];

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
expectZero["Sigma_chi = chi_1/chi_0", (sigmaChi /. slipSubs) - chi1/chi0];
expectZero["Sigma_eta = eta_1/eps_eta", (sigmaEta /. slipSubs) - eta1/epsEta];
expectZero["Sigma_del = delta_U,1/delta_U", (sigmaDel /. slipSubs) - deltaU1/deltaU];
expectZero["Sigma_eps = varepsilon_W/eps_W", (sigmaEps /. slipSubs) - varepsW/epsW];
expectZero["Sigma_Z = zeta_Z - omega_W", (sigmaZ /. slipSubs) - (zetaZ - omegaW)];

banner["Four-slippage grouped-defect law"];
xi1Direct = FullSimplify[zetaZ - omegaW + 2*chi1/(1 + chi0) + 2*eps1/(1 - eps), Assumptions -> $Assumptions];
sigmaEqns = {
  sigmaZ == 2*lam1 + mu1 - kEta - 2*kW,
  sigmaChi == gam1 + c1 - kU,
  sigmaEta == 2*c1 - kU - kEta,
  sigmaEps == 2*gam1 + 2*lam1 - kU - kW,
  sigmaDel == tau1 - kU
};
sigmaSolve = First[Solve[sigmaEqns, {mu1, gam1, kEta, kW, tau1}]];
logSyms = {lam1, c1, gam1, kU, kEta, kW, mu1, tau1};
xi1DirectInSigmas = Collect[
  FullSimplify[xi1Direct /. sigmaSolve, Assumptions -> $Assumptions],
  {sigmaZ, sigmaChi, sigmaEta, sigmaEps, sigmaDel},
  FullSimplify[#, Assumptions -> $Assumptions]&
];
xi1DirectInSigmas = xi1DirectInSigmas /. ConditionalExpression[e_, _] :> e;
xi1DirectInSigmas = FullSimplify[xi1DirectInSigmas, Assumptions -> $Assumptions];
If[FreeQ[xi1DirectInSigmas, Alternatives @@ logSyms],
  pass["Xi_1 microscopic gauges eliminated"],
  fail["Xi_1 microscopic gauges eliminated", xi1DirectInSigmas]
];
expectZero["Xi_1 coeff sigma_Z", Coefficient[Expand[xi1DirectInSigmas], sigmaZ] - 1];
expectZero["Xi_1 coeff sigma_chi", Coefficient[Expand[xi1DirectInSigmas], sigmaChi] - 2*chi0/(1 + chi0)];
expectZero["Xi_1 coeff sigma_eta", Coefficient[Expand[xi1DirectInSigmas], sigmaEta]];
expectZero["Xi_1 coeff sigma_eps", Coefficient[Expand[xi1DirectInSigmas], sigmaEps] - 2*epsW*(11 + 9*deltaU)/(11*(1 - eps)*(1 + deltaU))];
expectZero["Xi_1 coeff sigma_del", Coefficient[Expand[xi1DirectInSigmas], sigmaDel] + 4*epsW*deltaU/(11*(1 - eps)*(1 + deltaU)^2)];
expectZero["Xi_1 constant term", xi1DirectInSigmas /. {sigmaZ -> 0, sigmaChi -> 0, sigmaEta -> 0, sigmaEps -> 0, sigmaDel -> 0}];
xi1Slip = xi1DirectInSigmas;
expectZero["Xi_1 direct - derived slippage form", xi1Direct - (xi1Slip /. slipSubs)];
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
xi1Split = Collect[
  FullSimplify[
    xi1Slip /. sigmaChi -> (sigmaTr - (1 + chi0)*sigmaDel)/(1 + deltaU),
    Assumptions -> $Assumptions
  ],
  {sigmaZ, sigmaTr, sigmaEps, sigmaDel},
  FullSimplify[#, Assumptions -> $Assumptions]&
];
xi1Split = xi1Split /. ConditionalExpression[e_, _] :> e;
xi1Split = FullSimplify[xi1Split, Assumptions -> $Assumptions];
expectZero["Xi_1 split coeff sigma_Z", Coefficient[Expand[xi1Split], sigmaZ] - 1];
expectZero["Xi_1 split coeff sigma_tr", Coefficient[Expand[xi1Split], sigmaTr] - 2*chi0/((1 + chi0)*(1 + deltaU))];
expectZero["Xi_1 split coeff sigma_eps", Coefficient[Expand[xi1Split], sigmaEps] - 2*epsW*(11 + 9*deltaU)/(11*(1 - eps)*(1 + deltaU))];
expectZero["Xi_1 split coeff sigma_del", Coefficient[Expand[xi1Split], sigmaDel] + 2*chi0/(1 + deltaU) + 4*epsW*deltaU/(11*(1 - eps)*(1 + deltaU)^2)];
expectZero["Xi_1 split - slippage form", (xi1Split /. sigmaTr -> sigmaTrDef) - xi1Slip];
Print["Xi_1 split = ", fmt[xi1Split]];

banner["Support-blindness at the microscopic level"];
Print["free symbols of Xi_1: ", fmt[Sort[ToString /@ (List @@ xi1Slip /. Plus -> List /. Times -> List)]]];
(* Support-blindness: no support-lane drift enters the microscopic-log defect
   construction. The zeta-cancellation mechanism lives upstream in Stage 249. *)
Module[{forms = {{"Xi_1 direct", xi1Direct}, {"R_1 direct", r1Direct}, {"Theta_1 direct", theta1Direct}}},
  Do[
    With[{label = forms[[i, 1]], form = forms[[i, 2]]},
      If[FreeQ[form, lamphi1] && FreeQ[form, kphi],
        pass[label <> " support-blind"],
        fail[label <> " support-blind", form]
      ]
    ],
    {i, Length[forms]}
  ]
];

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
Print["Stage 182 Mathematica audit passed."];

Exit[0];
