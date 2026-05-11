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

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, cond]];
];

banner["STAGE 081 — FAMILY-1 SUPPORT IS AUTOMATIC"];

Clear[epsBlk, zMax];
$Assumptions =
  Element[{epsBlk, zMax}, Reals] &&
  zMax > 1 && 0 <= epsBlk < 1/zMax;

zetaReq = FullSimplify[(4/3 - 1)/(1 - epsBlk*(2 - 4/3)), Assumptions -> $Assumptions];
dZ = FullSimplify[D[zetaReq, epsBlk], Assumptions -> $Assumptions];
zetaEdge = FullSimplify[zetaReq /. epsBlk -> 1/zMax, Assumptions -> $Assumptions];
gap = FullSimplify[zMax - zetaEdge, Assumptions -> $Assumptions];
gapExpected = FullSimplify[3*zMax*(zMax - 1)/(3*zMax - 2), Assumptions -> $Assumptions];

Print["zeta_req(eps) = ", fmt[zetaReq]];
Print["d zeta_req / d eps = ", fmt[dZ]];
Print["zeta_edge = ", fmt[zetaEdge]];
Print["zmax - zeta_edge = ", fmt[Factor[gap]]];

expectZero["zeta_req - 1/(3-2 eps)", zetaReq - 1/(3 - 2*epsBlk)];
expectZero["d zeta_req / d eps exact formula", dZ - 2/(3 - 2*epsBlk)^2];
expectZero["gap factorization", gap - gapExpected];
expectTrue["automatic-support gap is positive for zmax > 1", gap > 0];

zMaxF1 = SetPrecision[2.46752922945601, 20];
zetaEdgeF1 = N[zetaEdge /. zMax -> zMaxF1, 30];
gapF1 = N[zMaxF1 - zetaEdgeF1, 30];

Print["Family-1 zeta_edge = ", fmt[zetaEdgeF1]];
Print["Family-1 margin = ", fmt[gapF1]];
expectApprox["Family-1 zeta_edge numeric check", zetaEdgeF1, 0.456730991107963169017835980412, 10^-15];
expectApprox["Family-1 margin numeric check", gapF1, 2.01079823834804688464927835412, 10^-15];

Print[""];
Print["Stage 098 Mathematica audit passed."];

Exit[0];
