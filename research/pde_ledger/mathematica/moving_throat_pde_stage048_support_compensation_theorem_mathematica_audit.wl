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

banner["STAGE 031 — SUPPORT COMPENSATION THEOREM"];

Clear[xi, delta, r, eps, zeta, mMix, sReq, sCrit];
$Assumptions =
  Element[{xi, delta, r, eps, zeta, mMix, sReq, sCrit}, Reals] &&
  0 < xi < 1 && delta > 0 && r > 0 && 0 < eps < 1 && zeta > 0 &&
  mMix > 0 && sCrit > sReq && sReq > 1;

gTr = FullSimplify[9*xi*(delta + xi)/(9*delta + (9 + 2*r^2)*xi), Assumptions -> $Assumptions];
fTr = FullSimplify[
  (9*delta + (9 + 2*r^2)*xi)^2*(9*delta + (9 + 2*r)*xi)^2/
    (81*(1 - xi)*(9*delta^2 + 18*delta*xi + (9 + 2*r^2)*xi^2)^2),
  Assumptions -> $Assumptions
];
mCrit = Block[
  {$Assumptions = Element[{delta, r}, Reals] && delta > 0 && r > 0},
  FullSimplify[Limit[gTr, xi -> 1, Direction -> "FromBelow"], Assumptions -> $Assumptions]
];
dGdXi = Factor[D[gTr, xi]];

Print["G_tr = ", fmt[gTr]];
Print["M_crit = ", fmt[mCrit]];
Print["dG_tr/dxi = ", fmt[dGdXi]];
expectZero[
  "dG_tr/dxi formula",
  dGdXi - 9*(2*r^2*xi^2 + 9*delta^2 + 18*delta*xi + 9*xi^2)/(2*r^2*xi + 9*delta + 9*xi)^2
];
expectZero["F_tr(xi=0)-1", FullSimplify[fTr /. xi -> 0, Assumptions -> delta > 0 && r > 0] - 1];

softCoeff = Block[
  {$Assumptions = Element[{delta, r}, Reals] && delta > 0 && r > 0},
  FullSimplify[Limit[(1 - xi)*fTr, xi -> 1, Direction -> "FromBelow"], Assumptions -> $Assumptions]
];
softCoeffExpected = FullSimplify[
  (9*delta + 9 + 2*r^2)^2*(9*delta + 9 + 2*r)^2/(81*(9*delta^2 + 18*delta + 9 + 2*r^2)^2),
  Assumptions -> delta > 0 && r > 0
];
Print["softening coefficient for F_tr = ", fmt[softCoeff]];
expectZero["(1-xi) F_tr softening coefficient", softCoeff - softCoeffExpected];

mGap = Factor[FullSimplify[mCrit - gTr, Assumptions -> $Assumptions]];
Print["M_crit - G_tr = ", fmt[mGap]];
expectZero[
  "M_crit - G_tr formula",
  mGap - 9*(1 - xi)*(2*r^2*xi + 9*delta^2 + 9*delta*xi + 9*delta + 9*xi)/
    ((2*r^2 + 9*delta + 9)*(2*r^2*xi + 9*delta + 9*xi))
];

sEnhance = FullSimplify[1 + zeta*(1 - eps)/(1 - zeta*eps), Assumptions -> $Assumptions];
Print["S(zeta;eps) = ", fmt[sEnhance]];
expectZero["S(zeta=0)-1", FullSimplify[sEnhance /. zeta -> 0, Assumptions -> 0 < eps < 1] - 1];
expectZero["dS/dzeta - (1-eps)/(1-zeta eps)^2", D[sEnhance, zeta] - (1 - eps)/(1 - zeta*eps)^2];

poleCoeff = Block[
  {$Assumptions = Element[{eps}, Reals] && 0 < eps < 1},
  FullSimplify[Limit[(1/eps - zeta)*sEnhance, zeta -> 1/eps, Direction -> "FromBelow"], Assumptions -> $Assumptions]
];
Print["pole coefficient for S = ", fmt[poleCoeff]];
expectZero["pole coefficient formula", poleCoeff - (1 - eps)/eps^2];

zetaSolutions = Solve[sEnhance == sReq, zeta, Reals];
If[Length[zetaSolutions] != 1, fail["inverse-map solve count", Length[zetaSolutions]]];
zetaReq = FullSimplify[zeta /. First[zetaSolutions], Assumptions -> $Assumptions];
zetaCrit = FullSimplify[(sCrit - 1)/(1 + eps*(sCrit - 2)), Assumptions -> $Assumptions];

Print["zeta_req = ", fmt[zetaReq]];
Print["zeta_crit = ", fmt[zetaCrit]];
expectZero["inverse map S(zeta_req)-S_req", FullSimplify[sEnhance /. zeta -> zetaReq, Assumptions -> $Assumptions] - sReq];
expectZero["inverse map S(zeta_crit)-S_crit", FullSimplify[sEnhance /. zeta -> zetaCrit, Assumptions -> $Assumptions] - sCrit];

marginPole = Factor[FullSimplify[1/eps - zetaReq, Assumptions -> $Assumptions]];
marginBranch = Factor[FullSimplify[zetaCrit - zetaReq, Assumptions -> $Assumptions]];
Print["1/eps - zeta_req = ", fmt[marginPole]];
Print["zeta_crit - zeta_req = ", fmt[marginBranch]];
expectZero["pole margin formula", marginPole - (1 - eps)/(eps*(1 + eps*(sReq - 2)))];
expectZero[
  "branch margin formula",
  marginBranch - (sCrit - sReq)*(1 - eps)/((1 + eps*(sCrit - 2))*(1 + eps*(sReq - 2)))
];

dXiDzeta = FullSimplify[mMix*D[sEnhance, zeta]/D[gTr, xi], Assumptions -> $Assumptions];
Print["dxi_phys/dzeta = ", fmt[Factor[dXiDzeta]]];
expectZero[
  "dxi_phys/dzeta formula",
  dXiDzeta - mMix*(1 - eps)*(2*r^2*xi + 9*delta + 9*xi)^2/
    ((1 - zeta*eps)^2*9*(2*r^2*xi^2 + 9*delta^2 + 18*delta*xi + 9*xi^2))
];

Print[""];
Print["Stage 048 Mathematica audit passed."];

Exit[0];
