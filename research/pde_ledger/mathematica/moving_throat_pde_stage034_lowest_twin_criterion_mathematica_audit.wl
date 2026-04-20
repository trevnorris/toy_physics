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

banner["STAGE 034 — LOWEST TWIN CRITERION"];

Clear[xi, delta, r, lambda, eps, epsEta, chi0, mMix, zW, pBranch];
$Assumptions =
  Element[{xi, delta, r, lambda, eps, epsEta, chi0, mMix, zW, pBranch}, Reals] &&
  0 < xi < 1 && delta > 0 && r > 0 && lambda > 0 && 0 < eps < 1 &&
  0 < epsEta < 1 && chi0 > -1 && mMix > 0 && zW > 0 && pBranch > 0;

gTr = FullSimplify[9 xi (xi + delta)/(9 delta + (9 + 2 r^2) xi), Assumptions -> $Assumptions];
fTr = FullSimplify[
  (9 delta + (9 + 2 r^2) xi)^2 (9 delta + (9 + 2 r) xi)^2/
    (81 (1 - xi) (9 delta^2 + 18 delta xi + (9 + 2 r^2) xi^2)^2),
  Assumptions -> $Assumptions
];
piTr = FullSimplify[fTr gTr, Assumptions -> $Assumptions];
piExpected = FullSimplify[
  xi (xi + delta) (9 delta + (9 + 2 r) xi)^2 (9 delta + (9 + 2 r^2) xi)/
    (9 (1 - xi) (9 delta^2 + 18 delta xi + (9 + 2 r^2) xi^2)^2),
  Assumptions -> $Assumptions
];

Print["G_tr = ", fmt[gTr]];
Print["F_tr = ", fmt[fTr]];
Print["Pi_tr = ", fmt[piTr]];
expectZero["Pi_tr - expected closed form", piTr - piExpected];

pi0 = FullSimplify[Limit[piTr, xi -> 0, Direction -> "FromAbove"], Assumptions -> delta > 0 && r > 0];
pi1 = Limit[piTr, xi -> 1, Direction -> "FromBelow"];
Print["Pi_tr(xi->0+) = ", fmt[pi0]];
Print["Pi_tr(xi->1-) = ", fmt[pi1]];
If[pi0 =!= 0, fail["Pi_tr(xi->0+)"], None];
If[pi1 =!= Infinity, fail["Pi_tr(xi->1-)"], None];

cMix = FullSimplify[8 lambda (1 - eps)/Pi^2, Assumptions -> $Assumptions];
sReq = FullSimplify[pBranch/cMix, Assumptions -> $Assumptions];
zetaReq = FullSimplify[(sReq - 1)/(1 + eps (sReq - 2)), Assumptions -> $Assumptions];
zetaReqBranch = FullSimplify[zetaReq /. pBranch -> piTr, Assumptions -> $Assumptions];

Print["C_mix = ", fmt[cMix]];
Print["S_req = ", fmt[sReq]];
Print["zeta_req(Pi,C_mix,eps) = ", fmt[zetaReq]];
Print["Physical-branch zeta_req = ", fmt[zetaReqBranch]];
expectZero["zeta_req at Pi = C_mix", zetaReq /. pBranch -> cMix];
expectZero["zeta_req at Pi = 2 C_mix minus 1", (zetaReq /. pBranch -> 2 cMix) - 1];
expectZero["zeta_req - 1 at Pi = 2 C_mix", ((zetaReq - 1) /. pBranch -> 2 cMix)];

lambdaTwinReq = FullSimplify[Pi^2 piTr/(16 (1 - eps)), Assumptions -> $Assumptions];
mMixTwinReq = FullSimplify[gTr/2, Assumptions -> $Assumptions];
zWTwinReq = FullSimplify[Pi^2 (1 - epsEta) (1 - eps) gTr/(16 (1 + chi0)^2), Assumptions -> $Assumptions];
mMixFromZW = FullSimplify[8 zW (1 + chi0)^2/(Pi^2 (1 - epsEta) (1 - eps)), Assumptions -> $Assumptions];

Print["Lambda_twin,req = ", fmt[lambdaTwinReq]];
Print["M_mix^(twin,req) = ", fmt[mMixTwinReq]];
Print["Z_W^(twin,req) = ", fmt[zWTwinReq]];
expectZero["M_mix(Z_W^(twin,req)) - G_tr/2", (mMixFromZW /. zW -> zWTwinReq) - gTr/2];

xi2x = FullSimplify[
  (2 mMix (9 + 2 r^2) - 9 delta + Sqrt[(2 mMix (9 + 2 r^2) - 9 delta)^2 + 648 mMix delta])/18,
  Assumptions -> $Assumptions
];
Print["xi_(2x) = ", fmt[xi2x]];
expectZero["G_tr(xi_(2x)) - 2 M_mix", (gTr /. xi -> xi2x) - 2 mMix];

Print[""];
Print["Stage 034 Mathematica audit passed."];

Exit[0];
