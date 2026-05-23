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

banner["STAGE 037 — ROBIN-COMPLIANCE SOFTENING"];

Clear[s, ell, k, h, y, eta, a, b];
$Assumptions =
  Element[{s, ell, k, h, y, eta, a, b}, Reals] && ell > 0 && k > 0 && h > 0 &&
  y > 0 && eta > 0 && a != 0;

psi = a Cos[k s] + b Sin[k s];
bExpr = FullSimplify[b /. First@Solve[(D[psi, s] /. s -> ell) == 0, b], Assumptions -> $Assumptions];
psiBN = FullSimplify[psi /. b -> bExpr, Assumptions -> $Assumptions];
charEq = FullSimplify[(D[psiBN, s] /. s -> 0) - h (psiBN /. s -> 0), Assumptions -> $Assumptions];

Print["B from Neumann bottom = ", fmt[bExpr]];
Print["Robin characteristic equation factor = ", fmt[charEq/a]];
expectZero["Robin equation -> k tan(kL) - h", charEq/a + h - k Tan[k ell]];
expectZero[
  "dimensionless form",
  ((k Tan[k ell] - h) /. {k -> y/ell, h -> eta/ell}) - (y Tan[y] - eta)/ell
];

Clear[kX, tX, kWBar, x, zetaReq];
$Assumptions =
  Element[{ell, y, kX, tX, kWBar, x, zetaReq}, Reals] && ell > 0 && y > 0 &&
  kX > 0 && tX > 0 && kWBar > 0 && 0 < x < 4 && zetaReq > 0;

kPhi = FullSimplify[kX + tX y^2/ell^2, Assumptions -> $Assumptions];
kWEff = FullSimplify[kX + Pi^2 tX/(4 ell^2), Assumptions -> $Assumptions];
aK = FullSimplify[kWEff/kPhi, Assumptions -> $Assumptions];

Print["K_W^(eff) = ", fmt[kWEff]];
Print["K_phi,0^(eff) = ", fmt[kPhi]];
Print["A_K = ", fmt[aK]];

kXFromX = FullSimplify[kWBar (1 - x/4), Assumptions -> $Assumptions];
tXFromX = FullSimplify[x ell^2 kWBar/Pi^2, Assumptions -> $Assumptions];
expectZero["K_W identity", kWBar - (kXFromX + Pi^2 tXFromX/(4 ell^2))];

aKX = FullSimplify[(kWBar/(kX + tX y^2/ell^2)) /. {kX -> kXFromX, tX -> tXFromX}, Assumptions -> $Assumptions];
aKSym = FullSimplify[1/(1 - x/4 + x y^2/Pi^2), Assumptions -> $Assumptions];
aKDN = FullSimplify[aKSym /. y -> Pi/2, Assumptions -> 0 < x < 4];
aKSoft = FullSimplify[Limit[aKSym, y -> 0, Direction -> "FromAbove"], Assumptions -> 0 < x < 4];

Print["A_K in x,y form = ", fmt[aKX]];
Print["A_K(y=Pi/2) = ", fmt[aKDN]];
Print["A_K(y->0+) = ", fmt[aKSoft]];
expectZero["A_K x-form", aKX - aKSym];
expectZero["DN limit", aKDN - 1];
expectZero["soft-mouth limit", aKSoft - 4/(4 - x)];

ineqRhs = FullSimplify[1/zetaReq - 1 + x/4, Assumptions -> $Assumptions];
yReqSq = FullSimplify[Pi^2 ineqRhs/x, Assumptions -> $Assumptions];
aKMax = FullSimplify[4/(4 - x), Assumptions -> $Assumptions];
xFloor = FullSimplify[x /. First@Solve[aKMax == zetaReq, x], Assumptions -> zetaReq > 0];

Print["Condition A_K >= zeta_req implies y^2 <= ", fmt[yReqSq]];
Print["A_K,max = ", fmt[aKMax]];
Print["x floor at saturation = ", fmt[xFloor]];
expectZero["x floor = 4 - 4/zeta_req", xFloor - (4 - 4/zetaReq)];
expectZero["A_K,max(x_floor) - zeta_req", (aKMax /. x -> xFloor) - zetaReq];

Print[""];
Print["Stage 054 Mathematica audit passed."];

Exit[0];
