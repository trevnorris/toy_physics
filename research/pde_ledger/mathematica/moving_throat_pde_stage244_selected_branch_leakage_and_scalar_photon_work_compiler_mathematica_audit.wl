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

fmt[expr_] := ToString[InputForm[expr]];
pass[name_String] := Print["PASS: ", name];

fail[name_String, detail_] := (
  Print["FAIL: ", name];
  Print["  residual -> ", fmt[detail]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  res = res /. ConditionalExpression[value_, _] :> value;
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, statement_] := Module[{res},
  res = FullSimplify[statement, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 244 - SELECTED-BRANCH LEAKAGE AND SCALAR-PHOTON WORK COMPILER"];

Clear[
  w, lam, E0, muw, rho0, q, Lam, eps, varrho, etaleak,
  Rtr, Rtarget, epseta
];

$Assumptions =
  Element[{w, eps, varrho, etaleak, Rtr, Rtarget, epseta}, Reals] &&
  Element[{lam, E0, muw, rho0, q, Lam}, Reals] &&
  lam > 0 && E0 > 0 && muw > 0 && rho0 > 0 && q > 0 &&
  Lam > 0 && 0 < eps < 1;

subbanner["M1-M4. Gaussian leakage and one-mode work"];

projector = Exp[-w^2/lam^2]/(lam Sqrt[Pi]);
oddLane = 2 w Exp[-w^2/lam^2]/(Sqrt[Pi] lam^3);
fieldW = -E0 oddLane;
currentW = muw rho0 fieldW;
chargeCurrentW = q currentW;

boundaryTerm = FullSimplify[
  Limit[
    projector currentW,
    w -> Infinity,
    Assumptions -> lam > 0 && E0 > 0 && muw > 0 && rho0 > 0
  ] -
    Limit[
      projector currentW,
      w -> -Infinity,
      Assumptions -> lam > 0 && E0 > 0 && muw > 0 && rho0 > 0
    ],
  Assumptions -> $Assumptions
];
leakageIntegral = FullSimplify[
  Integrate[
    D[projector, w] currentW,
    {w, -Infinity, Infinity},
    Assumptions -> $Assumptions,
    GenerateConditions -> False
  ],
  Assumptions -> $Assumptions
];
bulkIntegral = FullSimplify[
  Integrate[
    chargeCurrentW fieldW,
    {w, -Infinity, Infinity},
    Assumptions -> $Assumptions,
    GenerateConditions -> False
  ],
  Assumptions -> $Assumptions
];

leakageClosed = Sqrt[2] E0 muw rho0/(2 Sqrt[Pi] lam^3);
bulkClosed = Sqrt[2] E0^2 muw q rho0/(2 Sqrt[Pi] lam^3);
sessionClosed = 2 E0^2 muw q rho0/lam^2;

expectZero["M1 boundary term", boundaryTerm];
expectZero["M2 leakage compiler", leakageIntegral - leakageClosed];
expectZero["M3 bulk work", bulkIntegral - bulkClosed];
expectZero["M4 work-leak relation", bulkIntegral - q E0 leakageIntegral];

subbanner["M5. Session scalar and quadratic law"];

sessionFromBulk = FullSimplify[
  2 Sqrt[2 Pi] lam bulkIntegral,
  Assumptions -> $Assumptions
];
sessionFromLeak = FullSimplify[
  4 Pi q lam^4 leakageIntegral^2/(muw rho0),
  Assumptions -> $Assumptions
];

expectZero[
  "M5 session scalar",
  sessionFromBulk - 2 E0^2 muw q rho0/lam^2
];
expectZero[
  "M5 quadratic law",
  sessionFromBulk - sessionFromLeak
];

subbanner["M6. Selected-support pullback"];

mixCapacity = 8 Lam (1 - eps)/Pi^2;
traceDemand = FullSimplify[(4/3) mixCapacity, Assumptions -> $Assumptions];
traceDemandRho = FullSimplify[
  traceDemand /. eps -> 1 - (3/2) varrho,
  Assumptions -> $Assumptions
];
pullField = etaleak traceDemand;

SpullEps = FullSimplify[
  leakageClosed /. E0 -> pullField,
  Assumptions -> $Assumptions
];
WbulkpullEps = FullSimplify[
  bulkClosed /. E0 -> pullField,
  Assumptions -> $Assumptions
];
WsesspullEps = FullSimplify[
  sessionClosed /. E0 -> pullField,
  Assumptions -> $Assumptions
];

Spull = FullSimplify[
  SpullEps /. eps -> 1 - (3/2) varrho,
  Assumptions -> $Assumptions
];
Wbulkpull = FullSimplify[
  WbulkpullEps /. eps -> 1 - (3/2) varrho,
  Assumptions -> $Assumptions
];
Wsesspull = FullSimplify[
  WsesspullEps /. eps -> 1 - (3/2) varrho,
  Assumptions -> $Assumptions
];

expectZero[
  "M6 Pi_tr",
  traceDemand - 32 Lam (1 - eps)/(3 Pi^2)
];
expectZero[
  "M6 Pi_tr varrho",
  traceDemandRho - 16 Lam varrho/Pi^2
];
expectZero[
  "M6 S pullback epsilon",
  SpullEps -
    16 Sqrt[2] etaleak muw rho0 Lam (1 - eps)/(3 Pi^(5/2) lam^3)
];
expectZero[
  "M6 Wbulk pullback epsilon",
  WbulkpullEps -
    512 Sqrt[2] etaleak^2 muw q rho0 Lam^2 (1 - eps)^2/
      (9 Pi^(9/2) lam^3)
];
expectZero[
  "M6 Wsess pullback epsilon",
  WsesspullEps -
    2048 etaleak^2 muw q rho0 Lam^2 (1 - eps)^2/
      (9 Pi^4 lam^2)
];

subbanner["M7. Support/orbit split"];

expectTrue[
  "Spull orbit-free",
  FreeQ[Spull, Rtr] && FreeQ[Spull, Rtarget] && FreeQ[Spull, epseta]
];
expectTrue[
  "Spull support-dependent",
  Not[FreeQ[Spull, Lam]] && Not[FreeQ[Spull, varrho]] &&
    Not[FreeQ[Spull, etaleak]]
];
expectTrue[
  "Wbulkpull orbit-free",
  FreeQ[Wbulkpull, Rtr] && FreeQ[Wbulkpull, Rtarget] &&
    FreeQ[Wbulkpull, epseta]
];
expectTrue[
  "Wbulkpull support-dependent",
  Not[FreeQ[Wbulkpull, Lam]] && Not[FreeQ[Wbulkpull, varrho]] &&
    Not[FreeQ[Wbulkpull, etaleak]]
];
expectTrue[
  "Wsesspull orbit-free",
  FreeQ[Wsesspull, Rtr] && FreeQ[Wsesspull, Rtarget] &&
    FreeQ[Wsesspull, epseta]
];
expectTrue[
  "Wsesspull support-dependent",
  Not[FreeQ[Wsesspull, Lam]] && Not[FreeQ[Wsesspull, varrho]] &&
    Not[FreeQ[Wsesspull, etaleak]]
];

subbanner["M8. Transport-orientation parity"];

expectZero[
  "M8 Spull odd in etaleak",
  (Spull /. etaleak -> -etaleak) + Spull
];
expectZero[
  "M8 Wbulkpull even in etaleak",
  (Wbulkpull /. etaleak -> -etaleak) - Wbulkpull
];
expectZero[
  "M8 Wsesspull even in etaleak",
  (Wsesspull /. etaleak -> -etaleak) - Wsesspull
];

subbanner["M9. Recovery slice"];

expectZero["M9 Spull eta=0", Spull /. etaleak -> 0];
expectZero["M9 Wbulkpull eta=0", Wbulkpull /. etaleak -> 0];
expectZero["M9 Wsesspull eta=0", Wsesspull /. etaleak -> 0];

Print[""];
Print["All Stage 244 Mathematica checks passed."];
