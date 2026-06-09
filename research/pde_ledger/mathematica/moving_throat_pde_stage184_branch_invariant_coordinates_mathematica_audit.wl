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

banner["STAGE 184 — EXACT BRANCH-INVARIANT COORDINATES"];

Clear[
  chi0, deltaU, epsEta, theta1, xi1, sigmaEta, r1,
  rTr, tShape, rTarget, lam0
];
$Assumptions = (
  Element[{chi0, deltaU, epsEta, theta1, xi1, sigmaEta, r1, rTr, tShape, rTarget, lam0}, Reals]
  && chi0 > 0 && deltaU > 0 && epsEta > 0
  && rTr > 0 && tShape > 0 && rTarget > 0 && lam0 > 0
);

bStar = FullSimplify[2*(1 + chi0 + deltaU)/deltaU, Assumptions -> $Assumptions];
cStar = FullSimplify[(1 + chi0)*(1 + deltaU)*(1 + chi0 + deltaU)/(chi0*deltaU), Assumptions -> $Assumptions];
Print["B_* = ", fmt[bStar]];
Print["C_* = ", fmt[cStar]];

sigmaTr = FullSimplify[-cStar*theta1, Assumptions -> $Assumptions];
sigmaNT = FullSimplify[xi1 + bStar*theta1, Assumptions -> $Assumptions];

branchVars = {rTr, tShape, epsEta, rTarget};
branchVelocities = {theta1*rTr, xi1*tShape, sigmaEta*epsEta, r1*rTarget};

firstVariation[expr_] := FullSimplify[
  (D[expr, #] & /@ branchVars).branchVelocities,
  Assumptions -> $Assumptions
];

logDrift[expr_] := FullSimplify[
  firstVariation[expr]/expr,
  Assumptions -> $Assumptions
];

banner["Exact branch identities"];
selectedBranchClosedForm = 1 - epsEta;
productDriftResidual = FullSimplify[
  logDrift[rTarget*tShape] - logDrift[selectedBranchClosedForm],
  Assumptions -> $Assumptions
];
Print["product-drift residual before branch law = ", fmt[productDriftResidual]];
targetDriftLaw = r1 -> FullSimplify[-xi1 - epsEta*sigmaEta/(1 - epsEta), Assumptions -> $Assumptions];
expectZero[
  "R_target T^2 drift - delta ln(1 - eps_eta)",
  productDriftResidual /. targetDriftLaw
];

banner["Tracking invariant"];
tTr = rTr^(-cStar);
dlnTtr = logDrift[tTr];
Print["delta ln T_* = ", fmt[dlnTtr]];
expectZero["delta ln T_* - Sigma_tr", dlnTtr - sigmaTr];

banner["Corrected nontracking composite"];
nTr = tShape*rTr^bStar;
dlnNtr = logDrift[nTr];
Print["delta ln N_* = ", fmt[dlnNtr]];
expectZero["delta ln N_* - Sigma_nt", dlnNtr - sigmaNT];

banner["Dressing coordinate and selected-branch complement"];
dlnEpsEta = logDrift[epsEta];
Print["delta ln eps_eta = ", fmt[dlnEpsEta]];
expectZero["delta ln eps_eta - Sigma_eta", dlnEpsEta - sigmaEta];

dlnEcomp = logDrift[selectedBranchClosedForm];
Print["delta ln[(R_target T^2)/Lambda0] = ", fmt[dlnEcomp]];
expectZero["selected-branch complement identity", dlnEcomp + epsEta*sigmaEta/(1 - epsEta)];

banner["Composite zero-defect theorem"];
expectZero["Sigma_tr as branch-invariant log drift", sigmaTr - dlnTtr];
expectZero["Sigma_nt as branch-invariant log drift", sigmaNT - dlnNtr];
expectZero["Sigma_eta as branch-invariant log drift", sigmaEta - dlnEpsEta];

zeroMap1 = FullSimplify[dlnTtr /. theta1 -> 0, Assumptions -> $Assumptions];
zeroMap2 = FullSimplify[dlnNtr /. {theta1 -> 0, xi1 -> 0}, Assumptions -> $Assumptions];
zeroMap3 = FullSimplify[dlnEpsEta /. sigmaEta -> 0, Assumptions -> $Assumptions];
expectZero["delta ln T_* | Theta1=0", zeroMap1];
expectZero["delta ln N_* | Theta1=Xi1=0", zeroMap2];
expectZero["delta ln eps_eta | Sigma_eta=0", zeroMap3];

Print[""];
Print["Carry-forward formulas:"];
Print["  T_*  = R_tr^(-C_*)"];
Print["  N_*  = T^2 R_tr^(B_*)"];
Print["  D_*  = eps_eta"];
Print["  delta ln T_* = Sigma_tr"];
Print["  delta ln N_* = Sigma_nt"];
Print["  delta ln D_* = Sigma_eta"];
Print["  (R_target T^2)/Lambda0 = 1 - eps_eta"];
Print["  delta ln[(R_target T^2)/Lambda0] = -(eps_eta/(1-eps_eta)) Sigma_eta"];
Print["  Full zero-defect theorem: invariance of (R_tr, N_*, eps_eta) at first grouped order."];

Print[""];
Print["Stage 184 Mathematica audit passed."];

Exit[0];
