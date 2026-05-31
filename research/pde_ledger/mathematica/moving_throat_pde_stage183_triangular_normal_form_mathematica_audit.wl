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

expectNonzero[name_String, expr_] := Module[{res},
  res = FullSimplify[expr, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], fail[name, res], pass[name]];
];

banner["STAGE 166 — TRIANGULAR NORMAL FORM OF THE COHERENT DEFECT"];

Clear[chi0, epsW, epsEta, deltaU, sigmaZ, sigmaChi, sigmaEta, sigmaEps, sigmaDel, sigmaTr, sigmaNT];
$Assumptions = Element[{chi0, epsW, epsEta, deltaU, sigmaZ, sigmaChi, sigmaEta, sigmaEps, sigmaDel, sigmaTr, sigmaNT}, Reals] &&
  chi0 > 0 && epsW > 0 && epsEta > 0 && deltaU > 0;

eps = FullSimplify[epsW*(1 - (2/11)*deltaU/(1 + deltaU)), Assumptions -> $Assumptions];
theta1 = FullSimplify[
  -chi0*deltaU*sigmaTr/((1 + chi0)*(1 + deltaU)*(1 + chi0 + deltaU)),
  Assumptions -> $Assumptions
];
xi1 = FullSimplify[
  sigmaZ +
  2*chi0*sigmaTr/((1 + chi0)*(1 + deltaU)) +
  2*epsW*(11 + 9*deltaU)*sigmaEps/(11*(1 - eps)*(1 + deltaU)) -
  (2*chi0/(1 + deltaU) + 4*epsW*deltaU/(11*(1 - eps)*(1 + deltaU)^2))*sigmaDel,
  Assumptions -> $Assumptions
];
r1 = FullSimplify[-epsEta*sigmaEta/(1 - epsEta) - xi1, Assumptions -> $Assumptions];

banner["Branch-adapted nontracking slippage"];
sigmaNTDef = FullSimplify[
  sigmaZ +
  2*epsW*(11 + 9*deltaU)*sigmaEps/(11*(1 - eps)*(1 + deltaU)) -
  (2*chi0/(1 + deltaU) + 4*epsW*deltaU/(11*(1 - eps)*(1 + deltaU)^2))*sigmaDel,
  Assumptions -> $Assumptions
];
Print["Sigma_nt = ", fmt[sigmaNTDef]];

aTr = FullSimplify[2*chi0/((1 + chi0)*(1 + deltaU)), Assumptions -> $Assumptions];
cTr = FullSimplify[chi0*deltaU/((1 + chi0)*(1 + deltaU)*(1 + chi0 + deltaU)), Assumptions -> $Assumptions];
expectZero["Xi_1 - (A_tr Sigma_tr + Sigma_nt)", xi1 - (aTr*sigmaTr + sigmaNTDef)];
expectZero["R_1 + Xi_1 + eps_eta/(1-eps_eta) Sigma_eta", r1 + xi1 + epsEta*sigmaEta/(1 - epsEta)];

banner["Triangular observable ledger"];
Print["Theta_1 = ", fmt[theta1]];
Print["Xi_1 = ", fmt[FullSimplify[aTr*sigmaTr + sigmaNT, Assumptions -> $Assumptions]]];
Print["R_1 + Xi_1 = ", fmt[FullSimplify[-epsEta*sigmaEta/(1 - epsEta), Assumptions -> $Assumptions]]];

banner["Exact inverse reconstruction"];
sigmaTrInv = FullSimplify[-((1 + chi0)*(1 + deltaU)*(1 + chi0 + deltaU))*theta1/(chi0*deltaU), Assumptions -> $Assumptions];
expectZero["Sigma_tr inverse", sigmaTrInv - sigmaTr];

ratio = FullSimplify[aTr/cTr, Assumptions -> $Assumptions];
Print["A_tr/C_tr = ", fmt[ratio]];
expectZero["A_tr/C_tr - 2(1+chi0+deltaU)/deltaU", ratio - 2*(1 + chi0 + deltaU)/deltaU];

sigmaNTInv = FullSimplify[xi1 + ratio*theta1, Assumptions -> $Assumptions];
expectZero["Sigma_nt inverse", sigmaNTInv - sigmaNTDef];

sigmaEtaInv = FullSimplify[-(1 - epsEta)*(r1 + xi1)/epsEta, Assumptions -> $Assumptions];
expectZero["Sigma_eta inverse", sigmaEtaInv - sigmaEta];

banner["Triple-rigidity theorem"];
(* Rigidity holds iff the triangular map is invertible on the branch, i.e. iff
   each diagonal prefactor is nonzero there. We test that non-trivial content;
   the inverse round-trips above confirm full invertibility. *)
dressingPref = FullSimplify[epsEta/(1 - epsEta), Assumptions -> $Assumptions];
expectNonzero["C_tr (Theta_1 <- Sigma_tr prefactor) nonzero on branch", cTr];
expectNonzero["A_tr (Xi_1 <- Sigma_tr feed-through) nonzero on branch", aTr];
expectNonzero["eps_eta/(1-eps_eta) (R_1+Xi_1 <- Sigma_eta prefactor) nonzero on branch", dressingPref];

Print[""];
Print["Carry-forward formulas:"];
Print["  Sigma_tr = (1+chi_0) Sigma_del + (1+delta_U) Sigma_chi"];
Print["  Sigma_nt = Sigma_Z"];
Print["             + [2 eps_W/(1-eps)] * [(11+9 delta_U)/(11(1+delta_U))] Sigma_eps"];
Print["             - [2 chi_0/(1+delta_U) + 4 eps_W delta_U/(11(1-eps)(1+delta_U)^2)] Sigma_del"];
Print["  Theta_1   = -C_tr Sigma_tr"];
Print["  Xi_1      = A_tr Sigma_tr + Sigma_nt"];
Print["  R_1+Xi_1  = -(eps_eta/(1-eps_eta)) Sigma_eta"];
Print["  Sigma_tr  = -((1+chi_0)(1+delta_U)(1+chi_0+delta_U)/(chi_0 delta_U)) Theta_1"];
Print["  Sigma_nt  = Xi_1 + 2(1+chi_0+delta_U)/delta_U * Theta_1"];
Print["  Sigma_eta = -((1-eps_eta)/eps_eta) (R_1 + Xi_1)"];

Print[""];
Print["Stage 183 Mathematica audit passed."];

Exit[0];
