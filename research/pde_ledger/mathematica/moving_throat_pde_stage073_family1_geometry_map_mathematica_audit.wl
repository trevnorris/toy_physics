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

banner["STAGE 073 — FAMILY-1 GEOMETRY MAP"];

Clear[lSym, aSym, ellSym];
Module[{lambdaStarSym, ellOverASym, lambdaEllSym},
  lambdaStarSym = lSym/aSym;
  ellOverASym = ellSym/aSym;
  lambdaEllSym = FullSimplify[lambdaStarSym/ellOverASym,
    Assumptions -> lSym > 0 && aSym > 0 && ellSym > 0];
  expectZero["Lambda_ell - L/ell (symbolic)", lambdaEllSym - lSym/ellSym];
];

epsilonR = 1/20;
lambdaStar = 37/20;
ellOverA = epsilonR;
lambdaEll = FullSimplify[lambdaStar/ellOverA];

Print["epsilon_r = ", fmt[epsilonR]];
Print["ell/a = ", fmt[ellOverA]];
Print["L/a = ", fmt[lambdaStar]];
Print["Lambda_ell = L/ell = ", fmt[lambdaEll]];
expectZero["Lambda_ell - 37", lambdaEll - 37];

Clear[km, tx, len, ell];
$Assumptions =
  Element[{km, tx, len, ell}, Reals] &&
  km > 0 && tx > 0 && len > 0 && ell > 0;

etaSym = km*len/tx;
eta = FullSimplify[etaSym /. km -> tx/ell, Assumptions -> $Assumptions];
Print["eta under K_m = T_X/ell -> ", fmt[eta]];
expectZero["eta - L/ell", eta - len/ell];
expectZero["eta(reference) - 37", (eta /. (len/ell) -> 37) - 37];

Print[""];
Print["Stage 073 Mathematica audit passed."];

Exit[0];
