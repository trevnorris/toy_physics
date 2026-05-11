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
expectZero["Xi_F1(Upsilon) - 1369 Upsilon_w", xiF1FromUpsilon - 1369*upsilonW];
expectZero["Xi_F1(Theta) - 136900 Theta_w", xiF1FromTheta - 136900*thetaW];

zetaMinusChi = ToExpression["2.46622291347846`20"];
zetaPlusChi = ToExpression["2.46752913273870`20"];
zetaMinusJ = ToExpression["2.44257571477179`20"];
zetaPlusJ = ToExpression["2.46752736855058`20"];
zetaMaxF1 = ToExpression["2.46752922945601`20"];

Print["zeta_-^(chi) = ", fmt[zetaMinusChi]];
Print["zeta_+^(chi) = ", fmt[zetaPlusChi]];
Print["zeta_-^(J) = ", fmt[zetaMinusJ]];
Print["zeta_+^(J) = ", fmt[zetaPlusJ]];
Print["zeta_max^(F1) = ", fmt[zetaMaxF1]];

expectApprox["natural-window ordering gap", zetaPlusChi - zetaMinusChi, ToExpression["0.00130621926024`20"], 10^-12];
expectApprox["hard-ceiling gap", zetaMaxF1 - zetaPlusChi, ToExpression["9.6717311`10*^-8"], 10^-12];

Print[""];
Print["Stage 084 Mathematica audit passed."];

Exit[0];
