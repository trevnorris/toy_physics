ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
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

expectMatrixZero[name_String, expr_] := Module[{res, target},
  res = FullSimplify[Map[Together[Expand[#]] &, expr, {2}], Assumptions -> $Assumptions];
  target = ConstantArray[0, Dimensions[res]];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === target], pass[name], fail[name, res]];
];

banner["STAGE 003 — MINIMAL BdG-WALL COUPLING"];

subbanner["I. Axisymmetric wall + scalar BdG reduction"];

Clear[t, omega, qaFun, qLFun, xaFun, xbFun, maa, maL, mLL, kaa, kaL, kLL, wa, wb, c1a, c1b, c2a, c2b];
$Assumptions =
  Element[{t, omega, maa, maL, mLL, kaa, kaL, kLL, wa, wb, c1a, c1b, c2a, c2b}, Reals] &&
  wa > 0 && wb > 0;

qa = qaFun[t];
qL = qLFun[t];
xa = xaFun[t];
xb = xbFun[t];

lRed =
  1/2 maa D[qa, t]^2 + maL D[qa, t] D[qL, t] + 1/2 mLL D[qL, t]^2
  - 1/2 kaa qa^2 - kaL qa qL - 1/2 kLL qL^2
  + 1/2 D[xa, t]^2 - 1/2 wa^2 xa^2
  + 1/2 D[xb, t]^2 - 1/2 wb^2 xb^2
  + c1a qa xa + c1b qa xb + c2a qL xa + c2b qL xb;

	Clear[qa0, qL0, xa0, xb0, vqa, vqL, vxa, vxb, aqa, aqL, axa, axb];
	lRed = lRed + (
	  -1/2 kaa qa^2 - kaL qa qL - 1/2 kLL qL^2
	  + 1/2 D[xa, t]^2 - 1/2 wa^2 xa^2
	  + 1/2 D[xb, t]^2 - 1/2 wb^2 xb^2
	  + c1a qa xa + c1b qa xb + c2a qL xa + c2b qL xb
	);
	lAlg = lRed /. {qaFun[t] -> qa0, qLFun[t] -> qL0, xaFun[t] -> xa0, xbFun[t] -> xb0,
	                Derivative[1][qaFun][t] -> vqa, Derivative[1][qLFun][t] -> vqL,
	                Derivative[1][xaFun][t] -> vxa, Derivative[1][xbFun][t] -> vxb};
	timeD[expr_] := (
	  D[expr, qa0] vqa + D[expr, qL0] vqL + D[expr, xa0] vxa + D[expr, xb0] vxb
	  + D[expr, vqa] aqa + D[expr, vqL] aqL + D[expr, vxa] axa + D[expr, vxb] axb
	);
	backEL = {qa0 -> qaFun[t], qL0 -> qLFun[t], xa0 -> xaFun[t], xb0 -> xbFun[t],
	          vqa -> Derivative[1][qaFun][t], vqL -> Derivative[1][qLFun][t],
	          vxa -> Derivative[1][xaFun][t], vxb -> Derivative[1][xbFun][t],
	          aqa -> Derivative[2][qaFun][t], aqL -> Derivative[2][qLFun][t],
	          axa -> Derivative[2][xaFun][t], axb -> Derivative[2][xbFun][t]};
	elQa = FullSimplify[(timeD[D[lAlg, vqa]] - D[lAlg, qa0]) /. backEL, Assumptions -> $Assumptions];
	elQL = FullSimplify[(timeD[D[lAlg, vqL]] - D[lAlg, qL0]) /. backEL, Assumptions -> $Assumptions];
	elXa = FullSimplify[(timeD[D[lAlg, vxa]] - D[lAlg, xa0]) /. backEL, Assumptions -> $Assumptions];
	elXb = FullSimplify[(timeD[D[lAlg, vxb]] - D[lAlg, xb0]) /. backEL, Assumptions -> $Assumptions];
	
	expectZero["qa equation",
	  elQa - (maa qaFun''[t] + maL qLFun''[t] + kaa qaFun[t] + kaL qLFun[t]
	          - c1a xaFun[t] - c1b xbFun[t])];
	expectZero["qL equation",
	  elQL - (maL qaFun''[t] + mLL qLFun''[t] + kaL qaFun[t] + kLL qLFun[t]
	          - c2a xaFun[t] - c2b xbFun[t])];
	expectZero["xa equation",
	  elXa - (xaFun''[t] + wa^2 xaFun[t] - c1a qaFun[t] - c2a qLFun[t])];
	expectZero["xb equation",
	  elXb - (xbFun''[t] + wb^2 xbFun[t] - c1b qaFun[t] - c2b qLFun[t])];
	
	mMat = {{maa, maL}, {maL, mLL}};
	kMat = {{kaa, kaL}, {kaL, kLL}};
	cMat = {{c1a, c1b}, {c2a, c2b}};
	oMat = DiagonalMatrix[{wa^2, wb^2}];
	dEff = FullSimplify[kMat - omega^2 mMat - cMat . LinearSolve[oMat - omega^2 IdentityMatrix[2], Transpose[cMat]], Assumptions -> $Assumptions];
	
	(* derive D_eff by eliminating Xa, Xb from the EL equations *)
	Clear[Qa, QL, Xa, Xb];
	ansatz = {qaFun[t] -> Qa Exp[-I omega t], qLFun[t] -> QL Exp[-I omega t],
	          xaFun[t] -> Xa Exp[-I omega t], xbFun[t] -> Xb Exp[-I omega t],
	          qaFun''[t] -> -omega^2 Qa Exp[-I omega t], qLFun''[t] -> -omega^2 QL Exp[-I omega t],
	          xaFun''[t] -> -omega^2 Xa Exp[-I omega t], xbFun''[t] -> -omega^2 Xb Exp[-I omega t]};
	elQaF = FullSimplify[(elQa /. ansatz)/Exp[-I omega t], Assumptions -> $Assumptions];
	elQLF = FullSimplify[(elQL /. ansatz)/Exp[-I omega t], Assumptions -> $Assumptions];
	elXaF = FullSimplify[(elXa /. ansatz)/Exp[-I omega t], Assumptions -> $Assumptions];
	elXbF = FullSimplify[(elXb /. ansatz)/Exp[-I omega t], Assumptions -> $Assumptions];
	xsol = Solve[{elXaF == 0, elXbF == 0}, {Xa, Xb}][[1]];
	elQaRed = FullSimplify[elQaF /. xsol, Assumptions -> $Assumptions];
	elQLRed = FullSimplify[elQLF /. xsol, Assumptions -> $Assumptions];
	dDerived = {{Coefficient[elQaRed, Qa], Coefficient[elQaRed, QL]},
	            {Coefficient[elQLRed, Qa], Coefficient[elQLRed, QL]}};
	expectMatrixZero["derived D0_eff vs Deff", dDerived - dEff];
	manual = {
	  {
	    kaa - maa omega^2 - c1a^2/(wa^2 - omega^2) - c1b^2/(wb^2 - omega^2),
    kaL - maL omega^2 - c1a c2a/(wa^2 - omega^2) - c1b c2b/(wb^2 - omega^2)
  },
  {
    kaL - maL omega^2 - c1a c2a/(wa^2 - omega^2) - c1b c2b/(wb^2 - omega^2),
    kLL - mLL omega^2 - c2a^2/(wa^2 - omega^2) - c2b^2/(wb^2 - omega^2)
  }
};

Print["D0_eff(omega) = ", fmt[dEff]];
expectMatrixZero["D0_eff - manual form", dEff - manual];

dEffSeries = Map[Normal[Series[#, {omega, 0, 4}]] &, dEff, {2}] // Expand;
kEff = {
  {kaa - c1a^2/wa^2 - c1b^2/wb^2, kaL - c1a c2a/wa^2 - c1b c2b/wb^2},
  {kaL - c1a c2a/wa^2 - c1b c2b/wb^2, kLL - c2a^2/wa^2 - c2b^2/wb^2}
};
mEff = {
  {maa + c1a^2/wa^4 + c1b^2/wb^4, maL + c1a c2a/wa^4 + c1b c2b/wb^4},
  {maL + c1a c2a/wa^4 + c1b c2b/wb^4, mLL + c2a^2/wa^4 + c2b^2/wb^4}
};
nEff = {
  {c1a^2/wa^6 + c1b^2/wb^6, c1a c2a/wa^6 + c1b c2b/wb^6},
  {c1a c2a/wa^6 + c1b c2b/wb^6, c2a^2/wa^6 + c2b^2/wb^6}
};
expectMatrixZero["series match", dEffSeries - (kEff - omega^2 mEff - omega^4 nEff)];

subbanner["II. One wall mode + one stable BdG mode"];

Clear[m, k, varpi2, g, eps, om2, w2, delta];
	$Assumptions =
	  Element[{m, k, varpi2, g, eps, om2, w2, delta}, Reals] &&
	  m > 0 && k > 0 && varpi2 > 0 && eps > 0 && om2 > 0 && delta > 0;
	
	dispersion = Expand[(k - m w2) (varpi2 - w2) - g^2];
	
	(* derive dispersion from a single-mode Lagrangian *)
	Clear[qFun, xFun, Q, X];
	lOne = m/2 D[qFun[t], t]^2 - k/2 qFun[t]^2 + 1/2 D[xFun[t], t]^2 - varpi2/2 xFun[t]^2 + g qFun[t] xFun[t];
	elQ = D[D[lOne, D[qFun[t], t]], t] - D[lOne, qFun[t]];
	elX = D[D[lOne, D[xFun[t], t]], t] - D[lOne, xFun[t]];
	ansatz1 = {qFun[t] -> Q Exp[-I omega t], xFun[t] -> X Exp[-I omega t],
	           qFun''[t] -> -omega^2 Q Exp[-I omega t], xFun''[t] -> -omega^2 X Exp[-I omega t]};
	elQF = FullSimplify[(elQ /. ansatz1)/Exp[-I omega t], Assumptions -> $Assumptions];
	elXF = FullSimplify[(elX /. ansatz1)/Exp[-I omega t], Assumptions -> $Assumptions];
	xSol1 = Solve[elXF == 0, X][[1]];
	elQRed = FullSimplify[elQF /. xSol1, Assumptions -> $Assumptions];
	dispersionDerived = FullSimplify[Together[elQRed (varpi2 - omega^2) / Q], Assumptions -> $Assumptions];
	expectZero["derived dispersion vs (k - m w2)(varpi2 - w2) - g^2",
	           (dispersionDerived /. omega^2 -> w2) - dispersion];
	minusRoot = FullSimplify[(om2 + varpi2 - Sqrt[(om2 - varpi2)^2 + 4 g^2/m])/2, Assumptions -> $Assumptions];
	plusRoot = FullSimplify[(om2 + varpi2 + Sqrt[(om2 - varpi2)^2 + 4 g^2/m])/2, Assumptions -> $Assumptions];

Print["omega_-^2 = ", fmt[minusRoot]];
Print["omega_+^2 = ", fmt[plusRoot]];
expectZero["minus-root solves dispersion", dispersion /. {k -> m om2, w2 -> minusRoot}];
expectZero["plus-root solves dispersion", dispersion /. {k -> m om2, w2 -> plusRoot}];
expectZero["sum of roots", FullSimplify[minusRoot + plusRoot - (om2 + varpi2), Assumptions -> $Assumptions]];
expectZero["product of roots", FullSimplify[minusRoot plusRoot - (om2 varpi2 - g^2/m), Assumptions -> $Assumptions]];

wallLike = FullSimplify[(minusRoot /. {varpi2 -> om2 + delta, g -> eps g}), Assumptions -> $Assumptions];
matterLike = FullSimplify[(plusRoot /. {varpi2 -> om2 + delta, g -> eps g}), Assumptions -> $Assumptions];
wallSeries = Expand[Normal[Series[wallLike, {eps, 0, 2}]]];
matterSeries = Expand[Normal[Series[matterLike, {eps, 0, 2}]]];

Print["wall-like pole through O(g^2) = ", fmt[wallSeries]];
Print["matter-like pole through O(g^2) = ", fmt[matterSeries]];
expectZero["wall-like shift", wallSeries - (om2 - eps^2 g^2/(m delta))];
expectZero["matter-like shift", matterSeries - (om2 + delta + eps^2 g^2/(m delta))];

subbanner["III. Grouped real P2 + BdG couplings"];

Clear[k20, k21, k22, m20, m21, m22, g20, g21, g22, w20, w21, w22, k2, m2, g2, wq];
$Assumptions =
  Element[{omega, k20, k21, k22, m20, m21, m22, g20, g21, g22, w20, w21, w22, k2, m2, g2, wq}, Reals] &&
  w20 > 0 && w21 > 0 && w22 > 0 && wq > 0;

d20 = FullSimplify[k20 - m20 omega^2 - g20^2/(w20^2 - omega^2), Assumptions -> $Assumptions];
d21 = FullSimplify[k21 - m21 omega^2 - g21^2/(w21^2 - omega^2), Assumptions -> $Assumptions];
d22 = FullSimplify[k22 - m22 omega^2 - g22^2/(w22^2 - omega^2), Assumptions -> $Assumptions];
dP2 = DiagonalMatrix[{d20, d21, d22}];
dP2s = Map[Normal[Series[#, {omega, 0, 4}]] &, dP2, {2}] // Expand;
d2coeffRaw = Map[Coefficient[#, omega, 2] &, dP2s, {2}];
d2coeffMat = DiagonalMatrix[{d2coeffRaw[[1, 1]], 2 d2coeffRaw[[2, 2]], 2 d2coeffRaw[[3, 3]]}];
T0 = (1/Sqrt[5]) DiagonalMatrix[{1, Sqrt[2], Sqrt[2]}];
Ta = (1/Sqrt[10]) DiagonalMatrix[{2, -1, -1}];
Tb = (1/Sqrt[2]) DiagonalMatrix[{0, 1, -1}];

d20s = dP2s[[1, 1]];
d21s = dP2s[[2, 2]];
d22s = dP2s[[3, 3]];

(* projections onto representation-theoretic basis *)
d2Bar = FullSimplify[Tr[d2coeffMat]/5, Assumptions -> $Assumptions];
a2 = FullSimplify[(2 d2coeffRaw[[1,1]] - d2coeffRaw[[2,2]] - d2coeffRaw[[3,3]])/10,
                  Assumptions -> $Assumptions];
b2 = FullSimplify[(d2coeffRaw[[2,2]] - d2coeffRaw[[3,3]])/2,
                  Assumptions -> $Assumptions];

Print["D20 = ", fmt[d20s]];
Print["D21 = ", fmt[d21s]];
Print["D22 = ", fmt[d22s]];
Print["d2bar = ", fmt[d2Bar]];
Print["a2 = ", fmt[a2]];
Print["b2 = ", fmt[b2]];

isoSubs = {k20 -> k2, k21 -> k2, k22 -> k2, m20 -> m2, m21 -> m2, m22 -> m2, g20 -> g2, g21 -> g2, g22 -> g2, w20 -> wq, w21 -> wq, w22 -> wq};
expectZero["isotropic a2", a2 /. isoSubs];
expectZero["isotropic b2", b2 /. isoSubs];
expectZero["isotropic D20-D21", (d20 - d21) /. isoSubs];
expectZero["isotropic D21-D22", (d21 - d22) /. isoSubs];

subbanner["IV. Harmonic selection rule"];

Clear[th, ph];
$Assumptions = Element[{th, ph}, Reals];

dOmega = Sin[th];
y00 = 1/(2 Sqrt[Pi]);
y20 = Sqrt[5] (3 Cos[th]^2 - 1)/(4 Sqrt[Pi]);
y21c = Sqrt[15] Sin[th] Cos[th] Cos[ph]/(2 Sqrt[Pi]);
y22c = Sqrt[15] Sin[th]^2 Cos[2 ph]/(4 Sqrt[Pi]);

i0020 = FullSimplify[Integrate[Integrate[y00 y20 dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
i0021c = FullSimplify[Integrate[Integrate[y00 y21c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
i0022c = FullSimplify[Integrate[Integrate[y00 y22c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
i2021c = FullSimplify[Integrate[Integrate[y20 y21c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
i2022c = FullSimplify[Integrate[Integrate[y20 y22c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
i21c22c = FullSimplify[Integrate[Integrate[y21c y22c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
norm00 = FullSimplify[Integrate[Integrate[y00 y00 dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
norm20 = FullSimplify[Integrate[Integrate[y20 y20 dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
norm21c = FullSimplify[Integrate[Integrate[y21c y21c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
norm22c = FullSimplify[Integrate[Integrate[y22c y22c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];

yList = {y00, y20, y21c, y22c};
overlap = FullSimplify[
  Table[Integrate[Integrate[yList[[i]] yList[[j]] dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}],
        {i, 1, 4}, {j, 1, 4}],
  Assumptions -> $Assumptions];
expectMatrixZero["spherical harmonic overlap matrix - identity",
                 overlap - IdentityMatrix[4]];

Print[""];
Print["Stage 003 Mathematica audit passed."];

Exit[0];
