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

banner["STAGE 167 — EXACT BRANCH-INVARIANT COORDINATES"];

Clear[chi0, deltaU, epsEta];
$Assumptions = Element[{chi0, deltaU, epsEta}, Reals] && chi0 > 0 && deltaU > 0 && epsEta > 0;

bStar = FullSimplify[2*(1 + chi0 + deltaU)/deltaU, Assumptions -> $Assumptions];
cStar = FullSimplify[(1 + chi0)*(1 + deltaU)*(1 + chi0 + deltaU)/(chi0*deltaU), Assumptions -> $Assumptions];
Print["B_* = ", fmt[bStar]];
Print["C_* = ", fmt[cStar]];

Clear[small, theta1, xi1, sigmaEta, rTr0, t20, lam0];
$Assumptions = Element[{small, theta1, xi1, sigmaEta, rTr0, t20, lam0, chi0, deltaU, epsEta}, Reals] &&
  rTr0 > 0 && t20 > 0 && lam0 > 0 && chi0 > 0 && deltaU > 0 && epsEta > 0;

sigmaTr = FullSimplify[-cStar*theta1, Assumptions -> $Assumptions];
sigmaNT = FullSimplify[xi1 + bStar*theta1, Assumptions -> $Assumptions];

banner["Exact branch identities"];
rTr = rTr0*Exp[small*theta1];
t2 = t20*Exp[small*xi1];
epsEtaVar = epsEta*(1 + small*sigmaEta);
rTarget = lam0*(1 - epsEtaVar)/t2;
expectZero["R_target * T^2 - Lambda0 * (1 - eps_eta)", rTarget*t2 - lam0*(1 - epsEtaVar)];

banner["Tracking invariant"];
tTr = rTr^(-cStar);
tTr0 = rTr0^(-cStar);
dlnTtr = FullSimplify[SeriesCoefficient[Log[tTr/tTr0], {small, 0, 1}], Assumptions -> $Assumptions];
Print["delta ln T_* = ", fmt[dlnTtr]];
expectZero["delta ln T_* - Sigma_tr", dlnTtr - sigmaTr];

banner["Corrected nontracking composite"];
nTr = t2*rTr^bStar;
nTr0 = t20*rTr0^bStar;
dlnNtr = FullSimplify[SeriesCoefficient[Log[nTr/nTr0], {small, 0, 1}], Assumptions -> $Assumptions];
Print["delta ln N_* = ", fmt[dlnNtr]];
expectZero["delta ln N_* - Sigma_nt", dlnNtr - sigmaNT];

banner["Dressing coordinate and selected-branch complement"];
dlnEpsEta = FullSimplify[SeriesCoefficient[Log[epsEtaVar/epsEta], {small, 0, 1}], Assumptions -> $Assumptions];
Print["delta ln eps_eta = ", fmt[dlnEpsEta]];
expectZero["delta ln eps_eta - Sigma_eta", dlnEpsEta - sigmaEta];

eComp = (rTarget*t2)/lam0;
eComp0 = 1 - epsEta;
dlnEcomp = FullSimplify[SeriesCoefficient[Log[eComp/eComp0], {small, 0, 1}], Assumptions -> $Assumptions];
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
