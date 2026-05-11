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

normalizeExpr[expr_] := Module[{res},
  res = If[
    MatrixQ[expr],
    Map[FullSimplify[#, Assumptions -> $Assumptions] &, expr, {2}],
    If[
      VectorQ[expr],
      Map[FullSimplify[#, Assumptions -> $Assumptions] &, expr],
      FullSimplify[expr, Assumptions -> $Assumptions]
    ]
  ];
  res
];

allZeroQ[expr_] := If[
  ListQ[expr],
  And @@ Flatten[Map[TrueQ[# == 0] &, expr, {ArrayDepth[expr]}]],
  TrueQ[expr == 0]
];

prettyArray[arr_] := If[VectorQ[arr], MatrixForm[{arr}], MatrixForm[arr]];

expectZero[name_String, expr_] := Module[{res},
  res = normalizeExpr[expr];
  If[ListQ[res],
    Print[name, " ="];
    Print[prettyArray[res]];
    If[allZeroQ[res], pass[name], fail[name, res]],
    Print[name, " = ", fmt[res]];
    If[allZeroQ[res], pass[name], fail[name, res]]
  ];
];

banner["STAGE 226 — RELAXED-CONSTRAINT BRANCH DECLARATION AND SHORT-RANGE OPEN-SYSTEM COMPILER"];

Clear[
  w, ellW, j0, E0, U, V, kU, kV, chiLam, fU, z, a, b, y, x, kappa,
  C6, C4, C2, D6, D4, D2
];

$Assumptions =
  Element[{w, ellW, j0, E0, U, V, z, a, b, y, C6, C4, C2, D6, D4, D2}, Reals] &&
  Element[{kU, kV, chiLam, fU, x, kappa}, Reals] &&
  kU > 0 && kV > 0 && chiLam > 0 && fU > 0 && x > 0 && kappa > 0;

subbanner["I. Exact open-system leakage/work lane"];

W = Exp[-w^2]/Sqrt[Pi];
jw = ellW j0 w Exp[-w^2];
Ew = E0 w Exp[-w^2];

boundary = Block[
  {$Assumptions = Element[{ellW, j0}, Reals]},
  FullSimplify[Limit[W jw, w -> Infinity] - Limit[W jw, w -> -Infinity]]
];
Sleak = FullSimplify[
  Integrate[D[W, w] jw, {w, -Infinity, Infinity}],
  Assumptions -> $Assumptions
];
Wwork = FullSimplify[
  Integrate[jw Ew, {w, -Infinity, Infinity}],
  Assumptions -> $Assumptions
];

expectedSleak = -Sqrt[2] ellW j0/4;
expectedWwork = Sqrt[2] Sqrt[Pi] E0 ellW j0/8;

expectZero["boundary term", boundary];
expectZero["exact leakage scalar", Sleak - expectedSleak];
expectZero["exact work scalar", Wwork - expectedWwork];
expectZero["S_leak vanishes on ell_w = 0", Sleak /. ellW -> 0];
expectZero["W_work vanishes on ell_w = 0", Wwork /. ellW -> 0];

subbanner["II. Exact non-rigid U/V response and drain"];

FUV = (1/2) kU U^2 + (1/2) kV V^2 - chiLam U V - fU U;
stationarity = {D[FUV, U], D[FUV, V]};
uvSol = First[Solve[stationarity == {0, 0}, {U, V}]];
Usol = FullSimplify[U /. uvSol, Assumptions -> $Assumptions];
Vsol = FullSimplify[V /. uvSol, Assumptions -> $Assumptions];
Uexpected = FullSimplify[kV fU/(kU kV - chiLam^2), Assumptions -> $Assumptions];
Vexpected = FullSimplify[chiLam fU/(kU kV - chiLam^2), Assumptions -> $Assumptions];
HUV = {
  {D[FUV, U, U], D[FUV, U, V]},
  {D[FUV, V, U], D[FUV, V, V]}
};
detH = FullSimplify[Det[HUV], Assumptions -> $Assumptions];
ratioVU = FullSimplify[Vsol/Usol, Assumptions -> $Assumptions];
drainUV = FullSimplify[chiLam Usol Vsol, Assumptions -> $Assumptions];
drainExpected = FullSimplify[
  chiLam^2 kV fU^2/(kU kV - chiLam^2)^2,
  Assumptions -> $Assumptions
];

expectZero["U solution", Usol - Uexpected];
expectZero["V solution", Vsol - Vexpected];
expectZero["det Hessian", detH - (kU kV - chiLam^2)];
expectZero["V/U ratio", ratioVU - chiLam/kV];
expectZero["positive drain compiler", drainUV - drainExpected];
expectZero["U vanishes on f_U = 0", Usol /. fU -> 0];
expectZero["V vanishes on f_U = 0", Vsol /. fU -> 0];
expectZero["V vanishes on chi = 0", Vsol /. chiLam -> 0];
expectZero["U reduces to f_U / k_U on chi = 0", (Usol /. chiLam -> 0) - fU/kU];

subbanner["III. Compensated sign-changing source profile"];

varsigma = 1 + a Cos[Pi z] + b Cos[2 Pi z];
varsigmaMean = FullSimplify[
  Integrate[varsigma, {z, 0, 1}],
  Assumptions -> $Assumptions
];
varsigmaY = Expand[
  varsigma /. {
    Cos[Pi z] -> y,
    Cos[2 Pi z] -> 2 y^2 - 1
  }
];
yStar = FullSimplify[-a/(4 b), Assumptions -> $Assumptions];
varsigmaVertex = FullSimplify[varsigmaY /. y -> yStar, Assumptions -> $Assumptions];
boundaryLeft = FullSimplify[varsigma /. z -> 0, Assumptions -> $Assumptions];
boundaryRight = FullSimplify[varsigma /. z -> 1, Assumptions -> $Assumptions];

expectZero["unit-mean source normalization", varsigmaMean - 1];
expectZero["trivial source slice", (varsigma /. {a -> 0, b -> 0}) - 1];
expectZero[
  "quadratic source rewrite",
  varsigmaY - (1 - b + a y + 2 b y^2)
];
expectZero["interior stationary point", yStar + a/(4 b)];
expectZero["vertex value", varsigmaVertex - (1 - b - a^2/(8 b))];
expectZero["left boundary value", boundaryLeft - (1 + a + b)];
expectZero["right boundary value", boundaryRight - (1 - a + b)];

subbanner["IV. Codimension-three recovery slice"];

recoveryRules = {ellW -> 0, fU -> 0, a -> 0, b -> 0};
Urec = FullSimplify[Usol /. recoveryRules, Assumptions -> $Assumptions];
Vrec = FullSimplify[Vsol /. recoveryRules, Assumptions -> $Assumptions];
Drec = FullSimplify[drainUV /. recoveryRules, Assumptions -> $Assumptions];
Srec = FullSimplify[Sleak /. recoveryRules, Assumptions -> $Assumptions];
Wrec = FullSimplify[Wwork /. recoveryRules, Assumptions -> $Assumptions];
varsigmaRec = FullSimplify[varsigma /. recoveryRules, Assumptions -> $Assumptions];

expectZero["recovered S_leak", Srec];
expectZero["recovered W_work", Wrec];
expectZero["recovered U", Urec];
expectZero["recovered V", Vrec];
expectZero["recovered drain", Drec];
expectZero["recovered varsigma", varsigmaRec - 1];

subbanner["V. Short-range kernel invariant"];

SQ = x^-3;
SY = Exp[-2 kappa x]/x;
QQ = FullSimplify[SQ^2, Assumptions -> $Assumptions];
QY = FullSimplify[SQ SY, Assumptions -> $Assumptions];
YY = FullSimplify[SY^2, Assumptions -> $Assumptions];

deltaVStat = -(1/2) (C6 QQ + 2 C4 QY + C2 YY);
VdynReal = -(1/2) (D6 QQ + 2 D4 QY + D2 YY);
limitQQ = FullSimplify[
  Limit[x QQ, x -> Infinity, Assumptions -> kappa > 0],
  Assumptions -> $Assumptions
];
limitQY = FullSimplify[
  Limit[x QY, x -> Infinity, Assumptions -> kappa > 0],
  Assumptions -> $Assumptions
];
limitYY = FullSimplify[
  Limit[x YY, x -> Infinity, Assumptions -> kappa > 0],
  Assumptions -> $Assumptions
];
limitStat = FullSimplify[
  Limit[x deltaVStat, x -> Infinity, Assumptions -> kappa > 0],
  Assumptions -> $Assumptions
];
limitDyn = FullSimplify[
  Limit[x VdynReal, x -> Infinity, Assumptions -> kappa > 0],
  Assumptions -> $Assumptions
];

expectZero["QQ source product", QQ - x^-6];
expectZero["QY source product", QY - Exp[-2 kappa x]/x^4];
expectZero["YY source product", YY - Exp[-4 kappa x]/x^2];
expectZero["lim x QQ", limitQQ];
expectZero["lim x QY", limitQY];
expectZero["lim x YY", limitYY];
expectZero["lim x deltaV_stat", limitStat];
expectZero["lim x Re V_dyn", limitDyn];

Print[""];
Print["All Stage 243 symbolic checks passed."];
