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

banner["STAGE 047 — EQUILIBRIUM ALIGNMENT"];

Clear[gPhi, kX, npp, i1, i2, hw];
$Assumptions =
  Element[{gPhi, kX, npp, i1, i2, hw}, Reals] &&
  gPhi > 0 && kX > 0 && npp > 0 && i1 > 0 && i2 > 0 && hw > 0;

osp = FullSimplify[gPhi*i1, Assumptions -> $Assumptions];
nss = FullSimplify[gPhi^2*i2, Assumptions -> $Assumptions];
c2 = FullSimplify[osp^2/(npp*nss), Assumptions -> $Assumptions];
gEq = FullSimplify[gPhi^2*i1/kX, Assumptions -> $Assumptions];

Print["O_(sigma phi) = ", fmt[osp]];
Print["N_(sigma sigma) = ", fmt[nss]];
Print["C_(sigma phi)^2 = ", fmt[c2]];
Print["G_eq = ", fmt[gEq]];

constSubs = {i1 -> npp/hw, i2 -> npp/hw^2};
c2Const = FullSimplify[c2 /. constSubs, Assumptions -> $Assumptions];
gEqConst = FullSimplify[gEq /. constSubs, Assumptions -> $Assumptions];

Print["C^2 | H=const = ", fmt[c2Const]];
Print["G_eq | H=const = ", fmt[gEqConst]];
expectZero["matched-layer coherence", c2Const - 1];
expectZero["matched-layer gain vs Stage-45 best-alignment formula", gEqConst - gPhi^2*npp/(kX*hw)];

Clear[w1, w2, h1, h2];
$Assumptions =
  Element[{w1, w2, h1, h2}, Reals] && w1 > 0 && w2 > 0 && h1 > 0 && h2 > 0;

nppDisc = FullSimplify[w1 + w2, Assumptions -> $Assumptions];
i1Disc = FullSimplify[w1/h1 + w2/h2, Assumptions -> $Assumptions];
i2Disc = FullSimplify[w1/h1^2 + w2/h2^2, Assumptions -> $Assumptions];
gapDisc = FullSimplify[nppDisc*i2Disc - i1Disc^2, Assumptions -> $Assumptions];
gapExpected = FullSimplify[w1*w2*(h1 - h2)^2/(h1^2*h2^2), Assumptions -> $Assumptions];

Print["N_pp I2 - I1^2 = ", fmt[gapDisc]];
Print["expected gap = ", fmt[gapExpected]];
expectZero["two-point Cauchy gap identity", gapDisc - gapExpected];

Clear[theta, lambda, phi, sigma];
$Assumptions =
  Element[{theta, lambda, phi, sigma, kX, hw, gPhi, npp}, Reals] &&
  theta > 0 && lambda > 0 && kX > 0 && hw > 0 && gPhi > 0 && npp > 0;

f = 1/2*theta*sigma^2 - lambda*sigma*phi + 1/2*kX*phi^2;
sigmaStat = FullSimplify[sigma /. First[Solve[D[f, sigma] == 0, sigma]], Assumptions -> $Assumptions];
fEff = FullSimplify[f /. sigma -> sigmaStat, Assumptions -> $Assumptions];

Print["sigma_stat = ", fmt[sigmaStat]];
Print["F_eff = ", fmt[fEff]];
expectZero["effective support softening", fEff - 1/2*(kX - lambda^2/theta)*phi^2];

thetaEq = FullSimplify[hw*(nss /. constSubs), Assumptions -> gPhi > 0 && npp > 0 && hw > 0];
lambdaEq = FullSimplify[gPhi*(osp /. constSubs), Assumptions -> gPhi > 0 && npp > 0 && hw > 0];
softEq = FullSimplify[lambdaEq^2/thetaEq, Assumptions -> gPhi > 0 && npp > 0 && hw > 0];

Print["(Lambda^2/Theta)_eq = ", fmt[softEq]];
expectZero["equilibrium softening equals g_phi^2 I1", softEq - gPhi^2*(npp/hw)];

Print[""];
Print["Stage 064 Mathematica audit passed."];

Exit[0];
