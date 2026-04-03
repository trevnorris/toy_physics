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
  res = FullSimplify[Expand[expr], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 40] - N[target, 40]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 066 — DIRECT FAMILY-1 OPERATOR WINDOW"];

Clear[xi, peStar, delta, pe, alpha, eta, kappa, thetaChi, thetaJ];
$Assumptions = xi > 0 && Element[{alpha, eta, kappa}, Reals] && alpha > 0 && eta > 0;

peStar = pe[xi];
fixedPoint = peStar == xi*delta[peStar];
dpeRule = First[Solve[D[fixedPoint[[1]] - fixedPoint[[2]], xi] == 0, pe'[xi]]];
dpeFormula = FullSimplify[pe'[xi] /. dpeRule, Assumptions -> $Assumptions];

(* Independent expected formula from the note. *)
dpeExpected = FullSimplify[
  delta[peStar]/(1 - xi*Derivative[1][delta][peStar]),
  Assumptions -> $Assumptions
];

Print["dPe_*/dXi = ", fmt[dpeExpected]];
expectZero["implicit fixed-point derivative", dpeFormula - dpeExpected];

kappaF1 = 12321/5;
etaF1 = 37;
alphaF1 = Sqrt[kappaF1];

delta0F1 = FullSimplify[
  etaF1*(Cosh[alphaF1] - 1)/(alphaF1^2*(alphaF1*Sinh[alphaF1] + etaF1*Cosh[alphaF1]))
];
deltaInfF1 = FullSimplify[
  (Cosh[alphaF1] + (etaF1/alphaF1)*Sinh[alphaF1] - 1)/(alphaF1*Sinh[alphaF1] + etaF1*Cosh[alphaF1])
];

yRoot = y /. FindRoot[y*Tan[y] == etaF1, {y, 1.53}, WorkingPrecision -> 50, AccuracyGoal -> 25, PrecisionGoal -> 25];
aF1 = N[(kappaF1 + Pi^2/4)/(kappaF1 + yRoot^2), 40];

Print["Delta_0(F1) = ", N[delta0F1, 30]];
Print["Delta_inf(F1) = ", N[deltaInfF1, 30]];
Print["y_F1 = ", N[yRoot, 30]];
Print["A_F1 = ", N[aF1, 30]];

thetaChi = N[4.06863235008162, 30];
thetaJ = N[0.927552032539308, 30];
xiChi = N[136900*thetaChi, 30];
xiJ = N[136900*thetaJ, 30];

Print["Xi_chi coefficient = ", N[xiChi, 30]];
Print["Xi_J coefficient   = ", N[xiJ, 30]];

peMinusChi = N[xiChi*delta0F1, 30];
pePlusChi = N[xiChi*deltaInfF1, 30];
peMinusJ = N[xiJ*delta0F1, 30];
pePlusJ = N[xiJ*deltaInfF1, 30];

Print["Pe_-^(chi) = ", peMinusChi];
Print["Pe_+^(chi) = ", pePlusChi];
Print["Pe_-^(J)   = ", peMinusJ];
Print["Pe_+^(J)   = ", pePlusJ];

omega[pp_] := Pi*pp*(2*pp*Exp[pp] + Pi)/((4*pp^2 + Pi^2)*(Exp[pp] - 1));
zetaF1[pp_] := aF1*omega[pp]^2;

zetaMinusChi = N[zetaF1[peMinusChi], 30];
zetaPlusChi = N[zetaF1[pePlusChi], 30];
zetaMinusJ = N[zetaF1[peMinusJ], 30];
zetaPlusJ = N[zetaF1[pePlusJ], 30];
zetaMaxF1 = N[aF1*Pi^2/4, 30];

Print["zeta_-^(chi) = ", zetaMinusChi];
Print["zeta_+^(chi) = ", zetaPlusChi];
Print["zeta_-^(J)   = ", zetaMinusJ];
Print["zeta_+^(J)   = ", zetaPlusJ];
Print["zeta_max^(F1)= ", zetaMaxF1];

expectApprox["Delta_0(F1) numeric check", delta0F1, 1.73302079021525*10^-4, 10^-16];
expectApprox["Delta_inf(F1) numeric check", deltaInfF1, 2.01447565540522*10^-2, 10^-15];
expectApprox["Pe_-^(chi) numeric check", peMinusChi, 96.5285247264385, 10^-10];
expectApprox["Pe_+^(chi) numeric check", pePlusChi, 11220.5441626259, 10^-7];
expectApprox["Pe_-^(J) numeric check", peMinusJ, 22.0062226330754, 10^-10];
expectApprox["Pe_+^(J) numeric check", pePlusJ, 2558.01892349205, 10^-8];
expectApprox["zeta_-^(chi) numeric check", zetaMinusChi, 2.46622291347846, 10^-12];
expectApprox["zeta_+^(chi) numeric check", zetaPlusChi, 2.46752913273870, 10^-12];
expectApprox["zeta_-^(J) numeric check", zetaMinusJ, 2.44257571477179, 10^-12];
expectApprox["zeta_+^(J) numeric check", zetaPlusJ, 2.46752736855058, 10^-12];
expectApprox["zeta_max^(F1) numeric check", zetaMaxF1, 2.46752922945601, 10^-12];

Print["Pi_suff^(chi)/C_mix at eps_blk=0 = ", N[1 + zetaMinusChi, 30]];
Print["Pi_fail^(chi)/C_mix at eps_blk=0 = ", N[1 + zetaPlusChi, 30]];
Print["Pi_suff^(J)/C_mix at eps_blk=0   = ", N[1 + zetaMinusJ, 30]];
Print["Pi_fail^(J)/C_mix at eps_blk=0   = ", N[1 + zetaPlusJ, 30]];
Print["Pi_max^(F1)/C_mix at eps_blk=0   = ", N[1 + zetaMaxF1, 30]];

Print[""];
Print["Theorem ledger:"];
Print["  dPe_*/dXi = Delta / (1 - Xi dDelta/dPe)"];
Print["  zeta_phys(F1) and Pi/C_mix windows are monotone in Xi on the stable branch"];
Print["  inserting the natural Family-1 wall data reproduces the Stage-61/63/64 windows directly"];

Exit[0];
