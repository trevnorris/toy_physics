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

banner["STAGE 040 — PHYSICAL (Pe, kappa, eta) PLACEMENT MAP"];

Clear[pe, kappa, y, zetaReq, x];
$Assumptions =
  Element[{pe, kappa, y, zetaReq, x}, Reals] && pe > 0 && kappa > 0 && y > 0 && zetaReq > 0 && x > 0;

omegaPe = FullSimplify[
  Pi pe (2 pe Exp[pe] + Pi)/((4 pe^2 + Pi^2) (Exp[pe] - 1)),
  Assumptions -> $Assumptions
];
aKX = FullSimplify[1/(1 - x/4 + x y^2/Pi^2), Assumptions -> $Assumptions];
xSub = FullSimplify[Pi^2/(kappa + Pi^2/4), Assumptions -> $Assumptions];
aKKappa = FullSimplify[aKX /. x -> xSub, Assumptions -> $Assumptions];
zetaPhys = FullSimplify[omegaPe^2 aKKappa, Assumptions -> $Assumptions];

Print["x(kappa) = ", fmt[xSub]];
Print["A_K(eta;kappa) = ", fmt[aKKappa]];
Print["zeta_0^(Pe+R) = ", fmt[zetaPhys]];
expectZero["A_K - (kappa+Pi^2/4)/(kappa+y^2)", aKKappa - (kappa + Pi^2/4)/(kappa + y^2)];
expectZero[
  "zeta_phys - Omega_Pe^2*(kappa+Pi^2/4)/(kappa+y^2)",
  zetaPhys - omegaPe^2 (kappa + Pi^2/4)/(kappa + y^2)
];

dKappa = FullSimplify[D[zetaPhys, kappa], Assumptions -> $Assumptions];
dY = FullSimplify[D[zetaPhys, y], Assumptions -> $Assumptions];
Print["partial_kappa zeta = ", fmt[dKappa]];
Print["partial_y zeta = ", fmt[dY]];
expectZero[
  "partial_kappa identity",
  dKappa - omegaPe^2 (y^2 - Pi^2/4)/(kappa + y^2)^2
];
expectZero[
  "partial_y identity",
  dY + 2 omegaPe^2 y (kappa + Pi^2/4)/(kappa + y^2)^2
];

zetaMax = FullSimplify[(Pi^2/4) (kappa + Pi^2/4)/kappa, Assumptions -> $Assumptions];
kappaMax = FullSimplify[Pi^4/(4 (4 zetaReq - Pi^2)), Assumptions -> zetaReq > Pi^2/4];

Print["zeta_max(kappa) = ", fmt[zetaMax]];
Print["kappa_max(zeta_req) = ", fmt[kappaMax]];
expectZero[
  "zeta_max - limit(Pe->inf, y->0)",
  zetaMax - FullSimplify[Limit[Limit[zetaPhys, pe -> Infinity], y -> 0, Direction -> "FromAbove"], Assumptions -> kappa > 0]
];
expectZero["kappa_max identity", kappaMax - Pi^4/(4 (4 zetaReq - Pi^2))];
expectZero["zeta_max(kappa_max) - zeta_req", (zetaMax /. kappa -> kappaMax) - zetaReq];

omegaReqSq = FullSimplify[zetaReq (kappa + y^2)/(kappa + Pi^2/4), Assumptions -> $Assumptions];
yReqSq = FullSimplify[(omegaPe^2/zetaReq) (kappa + Pi^2/4) - kappa, Assumptions -> $Assumptions];
kappaReq = FullSimplify[(omegaPe^2 Pi^2/4 - zetaReq y^2)/(zetaReq - omegaPe^2), Assumptions -> $Assumptions];

Print["Omega_req^2 = ", fmt[omegaReqSq]];
Print["y_req^2 = ", fmt[yReqSq]];
Print["kappa_req = ", fmt[kappaReq]];
expectZero[
  "kappa_req defining equation",
  zetaReq - FullSimplify[(omegaPe^2 (kappa + Pi^2/4)/(kappa + y^2)) /. kappa -> kappaReq, Assumptions -> $Assumptions]
];
expectZero[
  "kappa_req identity",
  kappaReq - (omegaPe^2 Pi^2/4 - zetaReq y^2)/(zetaReq - omegaPe^2)
];
expectZero[
  "y_req identity",
  yReqSq - ((omegaPe^2/zetaReq) (kappa + Pi^2/4) - kappa)
];

Print[""];
Print["Stage 057 Mathematica audit passed."];

Exit[0];
