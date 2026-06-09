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

banner["STAGE 187 - EXACT ORBIT-QUOTIENT CLOSURE"];

Clear[
  dl, dc, dg, du, deta, dw, dm, dt,
  chiStar, deltaStar, eStar, fStar,
  lambdaW, cEtaU, gamma, kU, kEtaEff, kWEff, muW, tU,
  ell, sigma
];

deltaVector = {dl, dc, dg, du, deta, dw, dm, dt};
primitiveSymbols = {lambdaW, cEtaU, gamma, kU, kEtaEff, kWEff, muW, tU};
primitiveNames = {"lambda", "c", "gamma", "KU", "Keta", "KW", "mu", "T"};

$Assumptions = (
  Element[Join[deltaVector, {chiStar, deltaStar, eStar, fStar}], Reals] &&
  chiStar > 0 && deltaStar > 0 && eStar > 0 && fStar > 0 &&
  Element[Join[primitiveSymbols, {ell, sigma}], Reals] &&
  And @@ Thread[Join[primitiveSymbols, {ell, sigma}] > 0]
);

baseState = AssociationThread[primitiveNames -> primitiveSymbols];
targetState = AssociationThread[
  primitiveNames -> MapThread[#1 Exp[#2] &, {primitiveSymbols, deltaVector}]
];

cTrInvariant[state_Association] := (
  ((state["gamma"] state["c"])/state["KU"])^(1 + deltaStar) *
  ((Pi^2 state["T"])/(ell^2 state["KU"]))^(1 + chiStar)
);

cNtInvariant[state_Association] := (
  (state["lambda"]^2 state["mu"]/(state["Keta"] state["KW"]^2)) *
  ((state["gamma"]^2 state["lambda"]^2 sigma)/(state["KU"] state["KW"]))^eStar *
  ((Pi^2 state["T"])/(ell^2 state["KU"]))^(-fStar)
);

epsEtaInvariant[state_Association] := (
  state["c"]^2/(state["KU"] state["Keta"])
);

finiteLogRatio[invariant_] := Module[{ratio, logRatio},
  ratio = FullSimplify[invariant[targetState]/invariant[baseState], Assumptions -> $Assumptions];
  logRatio = PowerExpand[Log[ratio]];
  FullSimplify[Together[Expand[logRatio]], Assumptions -> $Assumptions]
];

rowLabels = {"C_tr", "C_nt", "epsilon_eta"};
derivedRows = finiteLogRatio /@ {cTrInvariant, cNtInvariant, epsEtaInvariant};

coefficientRow[expr_] := Coefficient[Expand[expr], #] & /@ deltaVector;
mStar = coefficientRow /@ derivedRows;
matrixRows = Expand[mStar.deltaVector];

Print["Exact finite log-ratio equations derived from physical monomials:"];
Do[
  Print["row_", rowLabels[[i]], " = ", fmt[derivedRows[[i]]]],
  {i, Length[rowLabels]}
];

Do[
  expectZero[
    "D1 " <> rowLabels[[i]] <> " log ratio - M_* Delta x",
    derivedRows[[i]] - matrixRows[[i]]
  ],
  {i, Length[rowLabels]}
];

minor = mStar[[All, {5, 7, 8}]];
Print["det selected minor (Delta_eta, Delta_mu, Delta_T) = ", fmt[Det[minor]]];
expectZero["selected minor determinant", Det[minor] - (1 + chiStar)];

solveLinearFor[expr_, var_] := FullSimplify[
  -((expr /. var -> 0)/Coefficient[expr, var]),
  Assumptions -> $Assumptions
];

detaSolved = solveLinearFor[derivedRows[[3]], deta];
dtSolved = solveLinearFor[derivedRows[[1]], dt];
dmSolved = solveLinearFor[derivedRows[[2]] /. {deta -> detaSolved, dt -> dtSolved}, dm];
finiteFiberRules = {deta -> detaSolved, dt -> dtSolved, dm -> dmSolved};

Print[""];
Print["Exact finite fibre solution:"];
Print["Delta_eta = ", fmt[detaSolved]];
Print["Delta_T = ", fmt[dtSolved]];
Print["Delta_mu = ", fmt[dmSolved]];

detaExpected = 2*dc - du;
dtExpected = du - (1 + deltaStar)*(dg + dc - du)/(1 + chiStar);
dmExpected = FullSimplify[
  2*dc - du + 2*dw - 2*dl
  - eStar*(2*dg + 2*dl - du - dw)
  - fStar*((1 + deltaStar)/(1 + chiStar))*(dg + dc - du),
  Assumptions -> $Assumptions
];

expectZero["Delta_eta finite law", detaSolved - detaExpected];
expectZero["Delta_T finite law", dtSolved - dtExpected];
expectZero["Delta_mu finite law", dmSolved - dmExpected];

Do[
  expectZero[
    "D4 " <> rowLabels[[i]] <> " row after finite fibre solve",
    derivedRows[[i]] /. finiteFiberRules
  ],
  {i, Length[rowLabels]}
];

banner["Finite orbit interpretation"];
Print["The three physical monomial equalities derive the rank-3 matrix condition"];
Print["M_* Delta x = 0, with Delta x an exact finite log-ratio vector."];
Print["The selected minor is positive, so the finite similarity fibre is unique."];

banner["Carry-forward formulas"];
Print["  Delta_eta = 2 Delta_c - Delta_U"];
Print["  Delta_T   = Delta_U - ((1+deltaU_*)/(1+chi_*))(Delta_gamma + Delta_c - Delta_U)"];
Print["  Delta_mu  = 2 Delta_c - Delta_U + 2 Delta_W - 2 Delta_lambda"];
Print["              - E_*(2 Delta_gamma + 2 Delta_lambda - Delta_U - Delta_W)"];
Print["              - F_*((1+deltaU_*)/(1+chi_*))(Delta_gamma + Delta_c - Delta_U)"];
Print[""];
Print["Conclusion:"];
Print["  Equal values of (C_tr, C_nt, epsilon_eta) are exactly equivalent to lying on"];
Print["  the same finite similarity orbit G_*. The weak-axisymmetric defect therefore"];
Print["  lives in the exact three-dimensional orbit quotient."];

Exit[0];
