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

(* --- Local static linear-response closure --- *)
Clear[yLoc, chiPhiLoc, hLoc, sigmaLoc];
Block[{},
  fLoc = (1/2) hLoc[yLoc] sigmaLoc^2 - gPhi chiPhiLoc[yLoc] sigmaLoc;
  closureSolutions = Solve[D[fLoc, sigmaLoc] == 0, sigmaLoc];
  If[Length[closureSolutions] =!= 1, fail["linear-response minimiser must be unique"]];
  chiSigmaClosure = sigmaLoc /. First[closureSolutions];
  Print["closure chi_sigma = ", fmt[chiSigmaClosure]];
  expectZero["closure law chi_sigma = g_phi chi_phi/H", chiSigmaClosure - gPhi*chiPhiLoc[yLoc]/hLoc[yLoc]];
];

(* --- Integral-level overlap invariants (concrete Gaussian profile) --- *)
Clear[yInt, lInt];
chiPhiG = Exp[-yInt^2/(2 lInt^2)];
hG = hw;
nppIntCheck = Integrate[chiPhiG^2, {yInt, -Infinity, Infinity}, Assumptions -> lInt > 0];
i1IntCheck  = Integrate[chiPhiG^2/hG, {yInt, -Infinity, Infinity}, Assumptions -> lInt > 0 && hw > 0];
i2IntCheck  = Integrate[chiPhiG^2/hG^2, {yInt, -Infinity, Infinity}, Assumptions -> lInt > 0 && hw > 0];
chiSigmaG = gPhi chiPhiG/hG;
ospIntCheck = Integrate[chiSigmaG*chiPhiG, {yInt, -Infinity, Infinity}, Assumptions -> lInt > 0 && hw > 0];
nssIntCheck = Integrate[chiSigmaG^2, {yInt, -Infinity, Infinity}, Assumptions -> lInt > 0 && hw > 0];
Print["Npp_int = ", fmt[nppIntCheck], "  I1_int = ", fmt[i1IntCheck], "  I2_int = ", fmt[i2IntCheck]];
expectZero["overlap O = g_phi * I1 (integral form)", ospIntCheck - gPhi*i1IntCheck];
expectZero["overlap N_ss = g_phi^2 * I2 (integral form)", nssIntCheck - gPhi^2*i2IntCheck];

osp = FullSimplify[gPhi*i1, Assumptions -> $Assumptions];
nss = FullSimplify[gPhi^2*i2, Assumptions -> $Assumptions];
c2 = FullSimplify[osp^2/(npp*nss), Assumptions -> $Assumptions];
gEq = FullSimplify[gPhi^2*i1/kX, Assumptions -> $Assumptions];

Print["O_(sigma phi) = ", fmt[osp]];
Print["N_(sigma sigma) = ", fmt[nss]];
Print["C_(sigma phi)^2 = ", fmt[c2]];
Print["G_eq = ", fmt[gEq]];

banner["MATCHED-LAYER INTEGRAL REDUCTION (concrete Gaussian profile)"];

Clear[y, L];
$Assumptions = Element[{y, L, hw}, Reals] && L > 0 && hw > 0;

chiPhiY = Exp[-y^2/(2 L^2)];
nppInt = Integrate[chiPhiY^2, {y, -Infinity, Infinity}, Assumptions -> L > 0];
i1Int  = Integrate[chiPhiY^2/hw, {y, -Infinity, Infinity}, Assumptions -> L > 0 && hw > 0];
i2Int  = Integrate[chiPhiY^2/hw^2, {y, -Infinity, Infinity}, Assumptions -> L > 0 && hw > 0];

Print["Npp_int = ", fmt[nppInt]];
Print["I1_int  = ", fmt[i1Int]];
Print["I2_int  = ", fmt[i2Int]];
expectZero["matched-layer I1 reduction", i1Int - nppInt/hw];
expectZero["matched-layer I2 reduction", i2Int - nppInt/hw^2];

$Assumptions =
  Element[{gPhi, kX, npp, i1, i2, hw}, Reals] &&
  gPhi > 0 && kX > 0 && npp > 0 && i1 > 0 && i2 > 0 && hw > 0;

c2Const = FullSimplify[(gPhi*i1Int)^2/(nppInt*gPhi^2*i2Int), Assumptions -> lInt > 0 && hw > 0 && gPhi > 0 && npp > 0];
gEqConst = FullSimplify[gPhi^2*i1Int/kX, Assumptions -> lInt > 0 && hw > 0 && gPhi > 0 && kX > 0];

Print["C^2 | H=const = ", fmt[c2Const]];
Print["G_eq | H=const = ", fmt[gEqConst]];
expectZero["matched-layer coherence", c2Const - 1];
expectZero["matched-layer gain vs Stage-45 best-alignment formula", gEqConst - gPhi^2*nppInt/(kX*hw)];

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

thetaAbs = FullSimplify[hw*nss, Assumptions -> gPhi > 0 && kX > 0 && npp > 0 && hw > 0 && i1 > 0 && i2 > 0];
lambdaAbs = FullSimplify[gPhi*osp, Assumptions -> gPhi > 0 && i1 > 0];
softAbs = FullSimplify[lambdaAbs^2/thetaAbs, Assumptions -> gPhi > 0 && kX > 0 && npp > 0 && hw > 0 && i1 > 0 && i2 > 0];
Print["Lambda^2/Theta (closure form) = ", fmt[softAbs]];
softMatched = FullSimplify[(softAbs /. i2 -> i1^2/npp) /. npp -> i1*hw, Assumptions -> gPhi > 0 && hw > 0 && i1 > 0];
Print["Lambda^2/Theta (matched layer) = ", fmt[softMatched]];
expectZero["equilibrium softening equals g_phi^2 I1", softMatched - gPhi^2*i1];

Print[""];
Print["Stage 064 Mathematica audit passed."];

Exit[0];
