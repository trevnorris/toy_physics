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

banner["STAGE 033 — PHYSICAL ZETA VS ZETA_REQ"];

Clear[sReq, eps, zeta, n, x];
$Assumptions =
  Element[{sReq, eps, zeta, x}, Reals] && sReq > 0 && 0 < eps < 1 && x > 0 &&
  Element[n, Integers] && n >= 1;

zetaReq = FullSimplify[(sReq - 1)/(1 + eps (sReq - 2)), Assumptions -> $Assumptions];
sEnhance = FullSimplify[1 + zeta (1 - eps)/(1 - eps zeta), Assumptions -> $Assumptions];

Print["zeta_req = ", fmt[zetaReq]];
Print["S(zeta;eps) = ", fmt[sEnhance]];
expectZero["S(1;eps) - 2", (sEnhance /. zeta -> 1) - 2];
criterion = FullSimplify[zetaReq - 1, Assumptions -> $Assumptions];
Print["zeta_req - 1 = ", fmt[criterion]];
expectZero[
  "zeta_req - 1 - (1-eps)(S_req-2)/(1+eps(S_req-2))",
  criterion - (1 - eps) (sReq - 2)/(1 + eps (sReq - 2))
];

zetaN = FullSimplify[1/((2 n + 1)^2 (1 + x n (n + 1))), Assumptions -> $Assumptions];
xEq = FullSimplify[(1/(((2 n + 1)^2) zetaReq) - 1)/(n (n + 1)), Assumptions -> $Assumptions];

Print["zeta_n^(twin) = ", fmt[zetaN]];
Print["x_max(n;zeta_req) = ", fmt[xEq]];
expectZero["zeta_n^(twin)(x_max) - zeta_req", (zetaN /. x -> xEq) - zetaReq];
expectZero[
  "x_max - [1/((2n+1)^2 zeta_req)-1]/[n(n+1)]",
  xEq - (1/(((2 n + 1)^2) zetaReq) - 1)/(n (n + 1))
];

supp = FullSimplify[(2 n + 1)^2 zetaN, Assumptions -> $Assumptions];
Print["(2n+1)^2 zeta_n^(twin) = ", fmt[supp]];
expectZero["suppression factor - 1/(1+x n(n+1))", supp - 1/(1 + x n (n + 1))];

sN = FullSimplify[sEnhance /. zeta -> zetaN, Assumptions -> $Assumptions];
sNMax = FullSimplify[1 + (1 - eps)/((2 n + 1)^2 - eps), Assumptions -> $Assumptions];
Print["S_n^(twin) = ", fmt[sN]];
Print["S_n^(max) = ", fmt[sNMax]];
expectZero["S_n^(twin)(x=0) - S_n^(max)", (sN /. x -> 0) - sNMax];
Print["S_1^(max) = ", fmt[sNMax /. n -> 1]];
Print["S_2^(max) = ", fmt[sNMax /. n -> 2]];

Print[""];
Print["Stage 050 Mathematica audit passed."];

Exit[0];
