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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 40] - N[target, 40]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 067 — FULL REDUCED MOVING-THROAT PDE WRITE-UP SKELETON"];

Clear[piTr, cMix, epsBlk, zeta, pe, kappa, eta, y, thetaW, upsilonW];
$Assumptions =
  Element[{piTr, cMix, epsBlk, zeta, pe, kappa, eta, y, thetaW, upsilonW}, Reals] &&
  cMix > 0 && pe > 0 && kappa > 0 && eta > 0 && thetaW > 0 && upsilonW > 0;

omegaPe = FullSimplify[
  Pi*pe*(2*pe*Exp[pe] + Pi)/((4*pe^2 + Pi^2)*(Exp[pe] - 1)),
  Assumptions -> $Assumptions
];
zetaPhys = FullSimplify[omegaPe^2*(kappa + Pi^2/4)/(kappa + y^2), Assumptions -> $Assumptions];
zetaReq = FullSimplify[(piTr - cMix)/(cMix - epsBlk*(2*cMix - piTr)), Assumptions -> $Assumptions];
qMap = FullSimplify[(1 + (1 - 2*epsBlk)*zeta)/(1 - epsBlk*zeta), Assumptions -> $Assumptions];
rQuad = FullSimplify[zetaReq - zetaPhys, Assumptions -> $Assumptions];

Print["zeta_phys(Pe,eta;kappa) = ", fmt[zetaPhys]];
Print["zeta_req(Pi_tr,C_mix,eps_blk) = ", fmt[zetaReq]];
Print["Q(zeta;eps_blk) = ", fmt[qMap]];
Print["R_quad = ", fmt[rQuad]];

expectZero["inverse demand map", (zetaReq /. piTr -> cMix*qMap) - zeta];

lambdaEll = 37;
kappaF1 = 12321/5;
etaF1 = 37;
xiF1FromUpsilon = lambdaEll^2*upsilonW;
xiF1FromTheta = 100*lambdaEll^2*thetaW;

Print["kappa_F1 = ", fmt[kappaF1]];
Print["eta_F1 = ", fmt[etaF1]];
Print["Xi_F1 from Upsilon_w = ", fmt[xiF1FromUpsilon]];
Print["Xi_F1 from Theta_w = ", fmt[xiF1FromTheta]];
expectZero["Xi_F1(Upsilon|Upsilon_w->100 Theta_w) - Xi_F1(Theta)", (xiF1FromUpsilon /. upsilonW -> 100*thetaW) - xiF1FromTheta];

zetaMinusChi = ToExpression["2.46622291347846`20"];
zetaPlusChi = ToExpression["2.46752913273870`20"];
zetaMinusJ = ToExpression["2.44257571477179`20"];
zetaPlusJ = ToExpression["2.46752736855058`20"];
zetaMaxF1 = ToExpression["2.46752922945601`20"];
yF1 = y /. FindRoot[y*Tan[y] - 37, {y, 153/100}, WorkingPrecision -> 40];
zetaPhysF1Limit = Limit[zetaPhys /. {kappa -> kappaF1, eta -> etaF1, y -> yF1}, pe -> Infinity];
zetaPhysF1Numeric = N[zetaPhysF1Limit, 20];
Print["zeta_phys(Pe->oo, kappa_F1, eta_F1, y_F1) = ", fmt[zetaPhysF1Numeric]];
expectApprox["zeta_phys at Family-1 (Pe->oo limit) matches upstream zeta_max^(F1)", zetaPhysF1Numeric, zetaMaxF1, 10^-10];

Print["zeta_-^(chi) = ", fmt[zetaMinusChi]];
Print["zeta_+^(chi) = ", fmt[zetaPlusChi]];
Print["zeta_-^(J) = ", fmt[zetaMinusJ]];
Print["zeta_+^(J) = ", fmt[zetaPlusJ]];
Print["zeta_max^(F1) = ", fmt[zetaMaxF1]];

expectZero["chi-window ordering positive (zeta_+^chi > zeta_-^chi)", If[TrueQ[zetaPlusChi > zetaMinusChi], 0, 1]];
expectZero["hard-ceiling gap positive (zeta_max^F1 > zeta_+^chi)", If[TrueQ[zetaMaxF1 > zetaPlusChi], 0, 1]];
expectZero["J-window ordering positive (zeta_+^J > zeta_-^J)", If[TrueQ[zetaPlusJ > zetaMinusJ], 0, 1]];
expectZero["fail-side J below chi (zeta_+^J <= zeta_+^chi)", If[TrueQ[zetaPlusJ <= zetaPlusChi], 0, 1]];

Print[""];
Print["Stage 084 Mathematica audit passed."];

Exit[0];
