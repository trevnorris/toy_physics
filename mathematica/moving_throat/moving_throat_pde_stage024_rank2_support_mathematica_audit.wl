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

banner["STAGE 024 — RANK-2 SUPPORT COMPLETION"];

Clear[a0, delta, xi, m, n, q, r, rU];
$Assumptions =
  Element[{a0, delta, xi, m, n, q, r, rU}, Reals] &&
  a0 > 0 && delta > 0 && xi > 0 && rU > 0;

lambda0 = 2/9;

subbanner["1. Exact rank-2 determinant and support-loading theorem"];

mMat = {
  {1 - m - n, -(m q + n r)},
  {-(m q + n r), 1 + delta - m q^2 - n r^2}
};

lam = 1 - xi;
detExpr = Expand[Det[mMat - lam IdentityMatrix[2]]];
detExpected = Expand[
  xi (delta + xi) - m (delta + (1 + q^2) xi) - n (delta + (1 + r^2) xi) + m n (q - r)^2
];

Print["det(M - (1-xi)I) = ", fmt[detExpr]];
expectZero["determinant decomposition", detExpr - detExpected];

nReq = FullSimplify[n /. First[Solve[detExpected == 0, n]], Assumptions -> $Assumptions];
nExpected = FullSimplify[
  (xi (delta + xi) - m (delta + (1 + q^2) xi))/(delta + (1 + r^2) xi - m (q - r)^2),
  Assumptions -> $Assumptions
];

Print["n_req = ", fmt[nReq]];
expectZero["n_req - expected", nReq - nExpected];

subbanner["2. Exact monotonicity with respect to mixed baseline loading"];

dndm = FullSimplify[D[nExpected, m], Assumptions -> $Assumptions];
monotoneExpected = FullSimplify[
  -(delta + (1 + q r) xi)^2/(delta + (1 + r^2) xi - m (q - r)^2)^2,
  Assumptions -> $Assumptions
];

Print["d n_req / d m = ", fmt[dndm]];
expectZero["dn/dm - expected", dndm - monotoneExpected];

subbanner["3. Tracking theorem: support follows the mixed direction"];

nTrack = FullSimplify[nExpected /. r -> q, Assumptions -> $Assumptions];
gQ = FullSimplify[xi (delta + xi)/(delta + (1 + q^2) xi), Assumptions -> $Assumptions];

Print["n_req(r=q) = ", fmt[nTrack]];
expectZero["tracking collapse", nTrack - (gQ - m)];

subbanner["4. Source-tied theorem: support remains aligned with the original source vector"];

q2Src = lambda0 rU^2;
r2Src = lambda0;
qrSrc = lambda0 rU;
qrDiffSq = lambda0 (rU - 1)^2;

nSrc = FullSimplify[
  (xi (delta + xi) - m (delta + (1 + q2Src) xi))/(delta + (1 + r2Src) xi - m qrDiffSq),
  Assumptions -> $Assumptions
];
nSrcExpected = FullSimplify[
  (xi (delta + xi) - m (delta + (1 + lambda0 rU^2) xi))/
    (delta + (1 + lambda0) xi - m lambda0 (rU - 1)^2),
  Assumptions -> $Assumptions
];

regThreshold = FullSimplify[(delta + (1 + lambda0) xi)/(lambda0 (rU - 1)^2), Assumptions -> $Assumptions];
numZeroThreshold = FullSimplify[xi (delta + xi)/(delta + (1 + lambda0 rU^2) xi), Assumptions -> $Assumptions];
dndmSrc = FullSimplify[D[nSrcExpected, m], Assumptions -> $Assumptions];
dndmSrcExpected = FullSimplify[
  -(delta + (1 + lambda0 rU) xi)^2/(delta + (1 + lambda0) xi - m lambda0 (rU - 1)^2)^2,
  Assumptions -> $Assumptions
];

Print["n_req^(src) = ", fmt[nSrc]];
expectZero["source-tied formula", nSrc - nSrcExpected];
Print["mixed-loading pole threshold = ", fmt[regThreshold]];
Print["mixed-loading ceiling from n_req >= 0 = ", fmt[numZeroThreshold]];
Print["d n_req^(src) / d m = ", fmt[dndmSrc]];
expectZero["source-tied dn/dm", dndmSrc - dndmSrcExpected];

Print[""];
Print["Stage 024 Mathematica audit passed."];

Exit[0];
