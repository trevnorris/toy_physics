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
  (* Strip ConditionalExpression wrapper: under $Assumptions, a result of
     the form ConditionalExpression[0, cond] is identically zero on the
     declared domain.  Solve[]/Reduce[] often introduce these wrappers when
     auxiliary inequalities are nontrivial. *)
  res = res /. ConditionalExpression[e_, _] :> e;
  res = FullSimplify[res, Assumptions -> $Assumptions];
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
xEqSolution = x /. First[Solve[zetaN == zetaReq, x]];
xEqSolution = xEqSolution /. ConditionalExpression[e_, _] :> e;
xEq = FullSimplify[xEqSolution, Assumptions -> $Assumptions];
xEqClosedForm = (1/(((2 n + 1)^2) zetaReq) - 1)/(n (n + 1));

Print["zeta_n^(twin) = ", fmt[zetaN]];
dZetaNdx = FullSimplify[D[zetaN, x], Assumptions -> $Assumptions];
dZetaNdxTarget = -n (n + 1) / ((2 n + 1)^2 (1 + n (n + 1) x)^2);
expectZero[
  "d zeta_n^(twin) / dx + n(n+1)/[(2n+1)^2 (1 + x n(n+1))^2]",
  dZetaNdx - dZetaNdxTarget
];
Print["x_max(n;zeta_req) = ", fmt[xEq]];
expectZero["zeta_n^(twin)(x_max) - zeta_req", (zetaN /. x -> xEq) - zetaReq];
expectZero[
  "x_max (from Solve) - [1/((2n+1)^2 zeta_req)-1]/[n(n+1)]",
  xEq - xEqClosedForm
];

admissibilityNum = FullSimplify[
  -xEq n (n + 1) - ((2 n + 1)^2 zetaReq - 1)/((2 n + 1)^2 zetaReq),
  Assumptions -> $Assumptions
];
Print["admissibility numerator residual = ", fmt[admissibilityNum]];
expectZero[
  "x_max non-negativity reduces to (2n+1)^2 zeta_req <= 1",
  admissibilityNum
];

sN = FullSimplify[sEnhance /. zeta -> zetaN, Assumptions -> $Assumptions];
sNMax = FullSimplify[1 + (1 - eps)/((2 n + 1)^2 - eps), Assumptions -> $Assumptions];
Print["S_n^(twin) = ", fmt[sN]];
Print["S_n^(max) = ", fmt[sNMax]];
expectZero["S_n^(twin)(x=0) - S_n^(max)", (sN /. x -> 0) - sNMax];
ceilingDiff = FullSimplify[sNMax - sN, Assumptions -> $Assumptions];
Print["S_n^(max) - S_n^(twin) = ", fmt[ceilingDiff]];
ceilingDiffTarget =
  ((1 - eps) (2 n + 1)^2 n (n + 1) x) /
  (((2 n + 1)^2 - eps) ((2 n + 1)^2 (1 + n (n + 1) x) - eps));
expectZero[
  "S_n^(max) - S_n^(twin) factored form",
  ceilingDiff - ceilingDiffTarget
];
Print["S_1^(max) = ", fmt[sNMax /. n -> 1]];
Print["S_2^(max) = ", fmt[sNMax /. n -> 2]];

Print[""];
Print["Stage 050 Mathematica audit passed."];

Exit[0];
