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

banner["STAGE 025 — RANK-2 SELECTED-MODE NORMALIZATION"];

Clear[a0, delta, xi, m, q, r, t, rU, eps];
$Assumptions =
  Element[{a0, delta, xi, m, q, r, t, rU, eps}, Reals] &&
  a0 > 0 && delta > 0 && xi > 0 && rU > 0 && eps > 0;

lambda0 = 2/9;

subbanner["1. Exact selected-mode eigenvector ratio after inserting Stage-24 support loading"];

nReq = FullSimplify[
  (xi (delta + xi) - m (delta + (1 + q^2) xi))/(delta + (1 + r^2) xi - m (q - r)^2),
  Assumptions -> $Assumptions
];

ratioRow1 = FullSimplify[(xi - m - nReq)/(m q + nReq r), Assumptions -> $Assumptions];
ratioRow2 = FullSimplify[(m q + nReq r)/(delta + xi - m q^2 - nReq r^2), Assumptions -> $Assumptions];
ratioExpected = FullSimplify[(m (q - r) + r xi)/(delta + xi - m q (q - r)), Assumptions -> $Assumptions];

Print["e1/e0 from row 1 = ", fmt[ratioRow1]];
Print["e1/e0 from row 2 = ", fmt[ratioRow2]];
Print["expected e1/e0 = ", fmt[ratioExpected]];
expectZero["row1 - expected", ratioRow1 - ratioExpected];
expectZero["row2 - expected", ratioRow2 - ratioExpected];

subbanner["2. Exact overlap formulas and generalized rank-2 normalization"];

dQR = FullSimplify[(delta + xi - m q (q - r))^2 + (m (q - r) + r xi)^2, Assumptions -> $Assumptions];
zOverlap = FullSimplify[(1 + q ratioExpected)^2/(1 + ratioExpected^2), Assumptions -> $Assumptions];
sOverlap = FullSimplify[(1 + t ratioExpected)^2/(1 + ratioExpected^2), Assumptions -> $Assumptions];

zExpected = FullSimplify[(delta + (1 + q r) xi)^2/dQR, Assumptions -> $Assumptions];
sExpected = FullSimplify[
  (delta + (1 + r t) xi - m (q - r) (q - t))^2/dQR,
  Assumptions -> $Assumptions
];
fGeneral = FullSimplify[zOverlap sOverlap/(1 - xi), Assumptions -> $Assumptions];
fExpected = FullSimplify[
  (delta + (1 + q r) xi)^2 (delta + (1 + r t) xi - m (q - r) (q - t))^2/
    ((1 - xi) dQR^2),
  Assumptions -> $Assumptions
];

Print["(z.e_-)^2 / z0^2 = ", fmt[zOverlap]];
Print["(s.e_-)^2 / s0^2 = ", fmt[sOverlap]];
Print["F_(q,r,t) = ", fmt[fGeneral]];
expectZero["Z_overlap - expected", zOverlap - zExpected];
expectZero["S_overlap - expected", sOverlap - sExpected];
expectZero["F_general - expected", fGeneral - fExpected];

subbanner["3. Tracking-support collapse back to Stage 23"];

fTrack = FullSimplify[fExpected /. r -> q, Assumptions -> $Assumptions];
fStage23 = FullSimplify[
  (delta + (1 + q^2) xi)^2 (delta + (1 + q t) xi)^2/
    ((1 - xi) ((delta + xi)^2 + q^2 xi^2)^2),
  Assumptions -> $Assumptions
];

Print["F_track = ", fmt[fTrack]];
expectZero["tracking collapse", fTrack - fStage23];

subbanner["4. Source-tied split-U specialization"];

qTiedQR = lambda0 rU;
a1 = FullSimplify[delta + (1 + qTiedQR) xi, Assumptions -> $Assumptions];
b1 = FullSimplify[delta + (1 + lambda0) xi - m lambda0 (rU - 1)^2, Assumptions -> $Assumptions];
dSrc = FullSimplify[
  (delta + xi - m lambda0 rU (rU - 1))^2 + lambda0 (xi + m (rU - 1))^2,
  Assumptions -> $Assumptions
];
fSrc = FullSimplify[a1^2 b1^2/((1 - xi) dSrc^2), Assumptions -> $Assumptions];
fSrcDirect = FullSimplify[
  fExpected /. {q -> Sqrt[lambda0] rU, r -> Sqrt[lambda0], t -> Sqrt[lambda0]},
  Assumptions -> $Assumptions
];

Print["F_src = ", fmt[fSrc]];
expectZero["source-tied specialization", fSrcDirect - fSrc];

subbanner["5. Exact flat-U recovery"];

fFlat = FullSimplify[
  (delta + (1 + lambda0) xi)^4/((1 - xi) ((delta + xi)^2 + lambda0 xi^2)^2),
  Assumptions -> $Assumptions
];
expectZero["F_src(R_U=1) - F_flat", (fSrc /. rU -> 1) - fFlat];

subbanner["6. First-order source-tied deformation about R_U = 1"];

nSrc = FullSimplify[
  (xi (delta + xi) - m (delta + (1 + lambda0 rU^2) xi))/
    (delta + (1 + lambda0) xi - m lambda0 (rU - 1)^2),
  Assumptions -> $Assumptions
];
gFlat = FullSimplify[xi (delta + xi)/(delta + (1 + lambda0) xi), Assumptions -> $Assumptions];
hNSrc = FullSimplify[D[nSrc, rU] /. rU -> 1, Assumptions -> $Assumptions];
hNExpected = FullSimplify[-2 lambda0 m xi/(delta + (1 + lambda0) xi), Assumptions -> $Assumptions];

fRatio = FullSimplify[fSrc/fFlat, Assumptions -> $Assumptions];
hFSrc = FullSimplify[D[fRatio, rU] /. rU -> 1, Assumptions -> $Assumptions];
hFExpected = FullSimplify[
  2 lambda0 (xi ((delta + xi)^2 + lambda0 xi^2) + 2 m delta (delta + (1 + lambda0) xi))/
    ((delta + (1 + lambda0) xi) ((delta + xi)^2 + lambda0 xi^2)),
  Assumptions -> $Assumptions
];

nLinear = FullSimplify[(gFlat - m) + eps hNExpected, Assumptions -> $Assumptions];
fLinear = FullSimplify[1 + eps hFExpected, Assumptions -> $Assumptions];

Print["H_n^(src) = ", fmt[hNSrc]];
Print["H_F^(src) = ", fmt[hFSrc]];
expectZero["H_n^(src) - expected", hNSrc - hNExpected];
expectZero["H_F^(src) - expected", hFSrc - hFExpected];
expectZero[
  "n_src - linear expansion",
  Expand[Normal[Series[nSrc /. rU -> 1 + eps, {eps, 0, 1}]] - nLinear]
];
expectZero[
  "F_src/F_flat - linear expansion",
  Expand[Normal[Series[fRatio /. rU -> 1 + eps, {eps, 0, 1}]] - fLinear]
];

Print[""];
Print["Stage 042 Mathematica audit passed."];

Exit[0];
