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

banner["STAGE 028 — COHERENT LOCAL TRACKING"];

Clear[lamW, lamPhi, gamma, muEta, muU, muW, muPhi, KU, gU, chi0, deltaU];
$Assumptions =
  Element[{lamW, lamPhi, gamma, muEta, muU, muW, muPhi, KU, gU, chi0, deltaU}, Reals] &&
  lamW > 0 && lamPhi > 0 && gamma > 0 && muEta > 0 && muU > 0 && muW > 0 &&
  muPhi > 0 && KU > 0 && chi0 > 0 && deltaU > 0;

gW = lamW/Sqrt[muEta*muW];
gR = gamma*lamW/Sqrt[muU*muW];
gB = lamPhi/Sqrt[muEta*muPhi];
gS = gamma*lamPhi/Sqrt[muU*muPhi];

rho0 = FullSimplify[gR*gU/(KU*gW), Assumptions -> $Assumptions];
sigma0 = FullSimplify[gU*gS/(KU*gB), Assumptions -> $Assumptions];
rTr = FullSimplify[(1 + chi0/(1 + deltaU))/(1 + chi0), Assumptions -> $Assumptions];

expectZero["g_B g_R - g_W g_S", gB*gR - gW*gS];
Print["rho_0 = ", fmt[rho0]];
Print["sigma_0 = ", fmt[sigma0]];
expectZero["rho_0 - sigma_0", rho0 - sigma0];

Print["R_tr = ", fmt[rTr]];
expectZero[
  "1 - R_tr - chi0 deltaU/((1+chi0)(1+deltaU))",
  (1 - rTr) - chi0*deltaU/((1 + chi0)*(1 + deltaU))
];
expectZero[
  "R_tr - 1/(1+deltaU) - deltaU/((1+chi0)(1+deltaU))",
  (rTr - 1/(1 + deltaU)) - deltaU/((1 + chi0)*(1 + deltaU))
];

Clear[ZW, ZPhi, epsEta, epsWSplit, epsPhiSplit];
$Assumptions =
  $Assumptions &&
  Element[{ZW, ZPhi, epsEta, epsWSplit, epsPhiSplit}, Reals] &&
  ZW > 0 && ZPhi > 0;

mMix = FullSimplify[8*ZW*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - epsWSplit)), Assumptions -> $Assumptions];
mSupp = FullSimplify[8*ZPhi*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - epsPhiSplit)), Assumptions -> $Assumptions];
mTr = FullSimplify[mMix + mSupp, Assumptions -> $Assumptions];
mTrExpected = FullSimplify[
  8*(1 + chi0)^2/(Pi^2*(1 - epsEta))*(ZW/(1 - epsWSplit) + ZPhi/(1 - epsPhiSplit)),
  Assumptions -> $Assumptions
];

Print["M_mix = ", fmt[mMix]];
Print["M_supp = ", fmt[mSupp]];
Print["M_tr = ", fmt[mTr]];
expectZero["M_tr - expected", mTr - mTrExpected];

Clear[xi, delta, lambda0, rU, rPhi, mMixSym, mSuppSym, mTrSym, rTarget];
$Assumptions =
  Element[{xi, delta, lambda0, rU, rPhi, mMixSym, mSuppSym, mTrSym, rTarget}, Reals] &&
  xi > 0 && delta > 0 && lambda0 > 0;

branchEq = FullSimplify[
  mSuppSym - (xi*(delta + xi) - mMixSym*(delta + (1 + lambda0*rU^2)*xi))/
    (delta + (1 + lambda0*rPhi^2)*xi - mMixSym*lambda0*(rU - rPhi)^2),
  Assumptions -> $Assumptions
];
branchTrack = Together[FullSimplify[branchEq /. rPhi -> rU, Assumptions -> $Assumptions]];
numTrack = Expand[Numerator[branchTrack]];
denTrack = Factor[Denominator[branchTrack]];
collapsedNum = Expand[xi^2 + (delta - mTrSym*(1 + lambda0*rU^2))*xi - delta*mTrSym];

Print["tracking numerator = ", fmt[numTrack]];
Print["tracking denominator = ", fmt[denTrack]];
expectZero["tracking quadratic collapse", numTrack + (collapsedNum /. mTrSym -> mMixSym + mSuppSym)];

mTrReq = FullSimplify[xi*(delta + xi)/(delta + (1 + lambda0*rU^2)*xi), Assumptions -> $Assumptions];
Print["M_tr required on tracking branch = ", fmt[mTrReq]];
expectZero["G_tr generic formula", mTrReq - xi*(delta + xi)/(delta + (1 + lambda0*rU^2)*xi)];

lambda0DN = 2/9;
gTrDN = FullSimplify[mTrReq /. lambda0 -> lambda0DN, Assumptions -> $Assumptions];
gTrExpected = FullSimplify[9*xi*(delta + xi)/(9*delta + (9 + 2*rU^2)*xi), Assumptions -> $Assumptions];
expectZero["G_tr D/N specialization", gTrDN - gTrExpected];

fTrack = FullSimplify[
  (delta + (1 + lambda0*rU^2)*xi)^2*(delta + (1 + lambda0*rU)*xi)^2/
    ((1 - xi)*((delta + xi)^2 + lambda0*rU^2*xi^2)^2),
  Assumptions -> $Assumptions
];
fTrExpected = FullSimplify[
  (9*delta + (9 + 2*rU^2)*xi)^2*(9*delta + (9 + 2*rU)*xi)^2/
    (81*(1 - xi)*(9*delta^2 + 18*delta*xi + (9 + 2*rU^2)*xi^2)^2),
  Assumptions -> $Assumptions
];
expectZero["F_tr normalization law", FullSimplify[fTrack /. lambda0 -> lambda0DN, Assumptions -> $Assumptions] - fTrExpected];
Print["coherent normalization residual = ", fmt[FullSimplify[rTarget - fTrExpected, Assumptions -> $Assumptions]]];

Print[""];
Print["Stage 028 Mathematica audit passed."];

Exit[0];
