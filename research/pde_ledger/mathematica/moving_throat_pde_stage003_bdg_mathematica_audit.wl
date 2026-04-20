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

staticL = lRed /. {D[qa, t] -> 0, D[qL, t] -> 0, D[xa, t] -> 0, D[xb, t] -> 0};
staticTmp = staticL /. {qa -> qa0, qL -> qL0, xa -> xa0, xb -> xb0};
staticBack = {qa0 -> qa, qL0 -> qL, xa0 -> xa, xb0 -> xb};
qad = D[qa, t];
qLd = D[qL, t];
xad = D[xa, t];
xbd = D[xb, t];
lVel = Expand[lRed /. {qad -> vqa, qLd -> vqL, xad -> vxa, xbd -> vxb}];

expectZero["qa kinetic coefficient", Coefficient[lVel, vqa^2] - maa/2];
expectZero["qL kinetic coefficient", Coefficient[lVel, vqL^2] - mLL/2];
expectZero["qa-qL mixed kinetic coefficient", Coefficient[lVel, vqa vqL] - maL];

mMat = {{maa, maL}, {maL, mLL}};
kMat = {{kaa, kaL}, {kaL, kLL}};
cMat = {{c1a, c1b}, {c2a, c2b}};
oMat = DiagonalMatrix[{wa^2, wb^2}];
dEff = FullSimplify[kMat - omega^2 mMat - cMat . LinearSolve[oMat - omega^2 IdentityMatrix[2], Transpose[cMat]], Assumptions -> $Assumptions];
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

d20s = Expand[Normal[Series[d20, {omega, 0, 4}]]];
d21s = Expand[Normal[Series[d21, {omega, 0, 4}]]];
d22s = Expand[Normal[Series[d22, {omega, 0, 4}]]];

d220 = FullSimplify[Coefficient[d20s, omega, 2], Assumptions -> $Assumptions];
d221 = FullSimplify[Coefficient[d21s, omega, 2], Assumptions -> $Assumptions];
d222 = FullSimplify[Coefficient[d22s, omega, 2], Assumptions -> $Assumptions];
d2Bar = FullSimplify[(d220 + 2 d221 + 2 d222)/5, Assumptions -> $Assumptions];
a2 = FullSimplify[(2 d220 - d221 - d222)/10, Assumptions -> $Assumptions];
b2 = FullSimplify[(d221 - d222)/2, Assumptions -> $Assumptions];

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
i2021c = FullSimplify[Integrate[Integrate[y20 y21c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
i2022c = FullSimplify[Integrate[Integrate[y20 y22c dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];
norm20 = FullSimplify[Integrate[Integrate[y20 y20 dOmega, {ph, 0, 2 Pi}], {th, 0, Pi}], Assumptions -> $Assumptions];

expectZero["Y00-Y20 cross integral", i0020];
expectZero["Y20-Y21c cross integral", i2021c];
expectZero["Y20-Y22c cross integral", i2022c];
expectZero["Y20 norm - 1", norm20 - 1];

Print[""];
Print["Stage 003 Mathematica audit passed."];

Exit[0];
