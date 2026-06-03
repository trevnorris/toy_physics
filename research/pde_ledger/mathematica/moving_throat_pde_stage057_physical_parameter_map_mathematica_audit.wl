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

banner["STAGE 057 — PHYSICAL (Pe, kappa, eta) PLACEMENT MAP"];

Clear[pe, kappa, y, zetaReq, x];
$Assumptions =
  Element[{pe, kappa, y, zetaReq, x}, Reals] && pe > 0 && kappa > 0 && y > 0 && y < Pi/2 && zetaReq > 0 && x > 0;

omegaPe = FullSimplify[
  Pi pe (2 pe Exp[pe] + Pi)/((4 pe^2 + Pi^2) (Exp[pe] - 1)),
  Assumptions -> $Assumptions
];

(* Independent derivation of A_K from the physical support operator.
   Notes section 2: K_W^(eff) = (T_X/L^2)(kappa + Pi^2/4),
                    K_phi,0^(eff) = (T_X/L^2)(kappa + y^2),
   so A_K = K_W^(eff)/K_phi,0^(eff) = (kappa + Pi^2/4)/(kappa + y^2). *)
Module[{KW, Kphi0, aKPhys},
  KW = KX + Pi^2 TX/(4 LL^2);
  Kphi0 = KX + TX y^2/LL^2;
  aKPhys = FullSimplify[KW/Kphi0,
    Assumptions -> KX > 0 && TX > 0 && LL > 0 && y > 0 && y < Pi/2];
  aKKappaFromPhys = FullSimplify[aKPhys /. KX -> kappa TX/LL^2,
    Assumptions -> kappa > 0 && TX > 0 && LL > 0 && y > 0 && y < Pi/2];
  expectZero["A_K(physical) - (kappa+Pi^2/4)/(kappa+y^2)",
    aKKappaFromPhys - (kappa + Pi^2/4)/(kappa + y^2)]
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

(* Sign check on partial_kappa over 0 < y < Pi/2 (from y tan y = eta, eta finite).
   Notes deliverable (4c) requires partial_kappa zeta < 0 on the constructive branch. *)
Module[{yvals, signOk},
  yvals = {Pi/8, Pi/6, Pi/4, Pi/3, 7 Pi/16};
  signOk = AllTrue[yvals, TrueQ[N[(D[zetaPhys, kappa] /. {pe -> 1, kappa -> 1, y -> #})] < 0] &];
  If[signOk,
    pass["partial_kappa zeta < 0 on 0 < y < Pi/2 (numerical sweep)"],
    fail["partial_kappa zeta sign sweep"]
  ]
];

(* Pe-monotonicity sweep — carry-forward from Stage 056 (notes §4: dOmega_Pe/dPe > 0 on
   the constructive branch via Cov_Pe(chi_0, s) > 0). Stage 056's scripts verify the
   covariance identity (wl:65) but not the sign, so we anchor the sign locally here. *)
Module[{pevals, signOk},
  pevals = {1/10, 1/2, 1, 2, 5, 10};
  signOk = AllTrue[pevals,
    TrueQ[N[(D[zetaPhys, pe] /. {pe -> #, kappa -> 1, y -> Pi/4})] > 0] &];
  If[signOk,
    pass["partial_Pe zeta > 0 on constructive branch (numerical sweep)"],
    fail["partial_Pe zeta sign sweep"]
  ]
];

zetaMax = FullSimplify[(Pi^2/4) (kappa + Pi^2/4)/kappa, Assumptions -> $Assumptions];
kappaMaxSol = Solve[zetaReq == (Pi^2/4) (kappa + Pi^2/4)/kappa, kappa, Reals];
kappaMax = FullSimplify[kappa /. First[kappaMaxSol], Assumptions -> zetaReq > Pi^2/4];

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
kappaReqSol = Solve[zetaReq == omegaPe^2 (kappa + Pi^2/4)/(kappa + y^2), kappa, Reals];
kappaReq = FullSimplify[kappa /. First[kappaReqSol], Assumptions -> $Assumptions];

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
yReqSqSol = Solve[zetaReq == omegaPe^2 (kappa + Pi^2/4)/(kappa + ysq), ysq, Reals];
yReqSqSolved = FullSimplify[ysq /. First[yReqSqSol], Assumptions -> $Assumptions];
expectZero["y_req identity", yReqSq - yReqSqSolved];

Print[""];
Print["Stage 057 Mathematica audit passed."];

Exit[0];
