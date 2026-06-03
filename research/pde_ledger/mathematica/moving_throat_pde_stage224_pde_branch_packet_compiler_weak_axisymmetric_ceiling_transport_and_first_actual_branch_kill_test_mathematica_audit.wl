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
fmt[expr_] := StringReplace[ToString[InputForm[expr]], "Global`" -> ""];
stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  detail -> ", fmt[detail]]];
  Exit[1];
);

cleanExprUnder[expr_, asum_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Cancel[Together[res]], Assumptions -> asum];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> asum]
];

expectZeroUnder[name_String, expr_, asum_] := Module[{res},
  res = cleanExprUnder[expr, asum];
  Print[name, " residual = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectZero[name_String, expr_] := expectZeroUnder[name, expr, $Assumptions];

expectTrueUnder[name_String, cond_, asum_] := Module[{res},
  res = stripConditional[FullSimplify[cond, Assumptions -> asum]];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

exactDecimal[s_String] := Module[
  {txt, sign, parts, whole, frac, scale, value},
  txt = s;
  sign = 1;
  If[StringStartsQ[txt, "-"],
    sign = -1;
    txt = StringDrop[txt, 1]
  ];
  parts = StringSplit[txt, "."];
  whole = ToExpression[parts[[1]]];
  If[Length[parts] == 1,
    value = whole,
    frac = parts[[2]];
    scale = 10^StringLength[frac];
    value = (whole scale + ToExpression[frac])/scale
  ];
  sign Rationalize[value, 0]
];

Clear[
  deltaNorm, tQuad, mhat0, pcrit,
  p20Lane, p21Lane, p22Lane,
  pbarUnknown, aUnknown, bUnknown,
  aSeed, bSeed, eps, xi1
];

$Assumptions = (
  Element[
    {
      deltaNorm, tQuad, mhat0, pcrit,
      p20Lane, p21Lane, p22Lane,
      pbarUnknown, aUnknown, bUnknown,
      aSeed, bSeed, eps, xi1
    },
    Reals
  ]
  && deltaNorm > 0 && tQuad > 0 && mhat0 > 0 && pcrit > 0
);

banner["STAGE 224 -- BRANCH PACKET COMPILER MATHEMATICA AUDIT"];

laneMatrix = {{1, 4, 0}, {1, -1, 1}, {1, -1, -1}};
pbarIso = (deltaNorm + tQuad)/mhat0^2;

subbanner["M1. Grouped inverse map from solved forward system"];

inverseRule = First[
  Solve[
    Thread[{p20Lane, p21Lane, p22Lane} == laneMatrix . {pbarUnknown, aUnknown, bUnknown}],
    {pbarUnknown, aUnknown, bUnknown}
  ]
];
packetLaneRules = {
  p20Lane -> pbarIso + 4 aSeed,
  p21Lane -> pbarIso - aSeed + bSeed,
  p22Lane -> pbarIso - aSeed - bSeed
};

Print["solved inverse rule = ", fmt[inverseRule]];
expectZero["M1 solved Pbar recovers input", (pbarUnknown /. inverseRule /. packetLaneRules) - pbarIso];
expectZero["M1 solved aP recovers input", (aUnknown /. inverseRule /. packetLaneRules) - aSeed];
expectZero["M1 solved bP recovers input", (bUnknown /. inverseRule /. packetLaneRules) - bSeed];

subbanner["M2. Isotropic ceiling rearrangement"];

expectZero[
  "M2 ceiling residual",
  (pbarIso - pcrit) mhat0^2 - (deltaNorm - (mhat0^2 pcrit - tQuad))
];

subbanner["M3. Weak-axisymmetric compiler via LinearSolve"];

lambdaList = {1, 1/2, -1};
waLanes = FullSimplify[pbarIso (1 + eps xi1 #) & /@ lambdaList, Assumptions -> $Assumptions];
waRecovered = FullSimplify[LinearSolve[laneMatrix, waLanes], Assumptions -> $Assumptions];
aPwa = waRecovered[[2]];
bPwa = waRecovered[[3]];

Print["WA lanes = ", fmt[waLanes]];
Print["LinearSolve recovered {Pbar, aP, bP} = ", fmt[waRecovered]];
expectZero["M3 aP weak-axisymmetric value", aPwa - eps pbarIso xi1/4];
expectZero["M3 bP weak-axisymmetric value", bPwa - 3 eps pbarIso xi1/4];
expectZero["M3 bP = 3 aP", bPwa - 3 aPwa];

subbanner["M4. Weak-axisymmetric lanes re-expand"];

waReexpanded = FullSimplify[laneMatrix . {pbarIso, aPwa, bPwa}, Assumptions -> $Assumptions];
expectZero["M4 P20 re-expansion", waReexpanded[[1]] - waLanes[[1]]];
expectZero["M4 P21 re-expansion", waReexpanded[[2]] - waLanes[[2]]];
expectZero["M4 P22 re-expansion", waReexpanded[[3]] - waLanes[[3]]];

subbanner["M5. Robust max-lane collapse by sign cases"];

positiveCase = $Assumptions && eps xi1 > 0;
negativeCase = $Assumptions && eps xi1 < 0;
absXiPositive = Refine[Abs[eps xi1], positiveCase];
absXiNegative = Refine[Abs[eps xi1], negativeCase];
absAPositive = Refine[Abs[aPwa], positiveCase];
absANegative = Refine[Abs[aPwa], negativeCase];

Print["Refined Abs[eps xi1] for eps xi1 > 0: ", fmt[absXiPositive]];
Print["Refined Abs[eps xi1] for eps xi1 < 0: ", fmt[absXiNegative]];

expectTrueUnder[
  "M5 positive sign P20 dominates",
  waLanes[[1]] >= waLanes[[2]] && waLanes[[1]] >= waLanes[[3]],
  positiveCase
];
expectZeroUnder[
  "M5 positive sign max lane equals Pbar(1+Abs[eps xi1])",
  waLanes[[1]] - pbarIso (1 + absXiPositive),
  positiveCase
];
expectTrueUnder[
  "M5 negative sign P22 dominates",
  waLanes[[3]] >= waLanes[[1]] && waLanes[[3]] >= waLanes[[2]],
  negativeCase
];
expectZeroUnder[
  "M5 negative sign max lane equals Pbar(1+Abs[eps xi1])",
  waLanes[[3]] - pbarIso (1 + absXiNegative),
  negativeCase
];
expectZeroUnder[
  "M5 positive sign robust aP form",
  pbarIso (1 + absXiPositive) - (pbarIso + 4 absAPositive),
  positiveCase
];
expectZeroUnder[
  "M5 negative sign robust aP form",
  pbarIso (1 + absXiNegative) - (pbarIso + 4 absANegative),
  negativeCase
];

subbanner["M6. Stage 223 carried headroom budget relations"];

(* Ceilings and compatibility point carried from Stage 223 dynamic-survival-window audit. *)
barP0compat = exactDecimal["0.002069792318062885"];
ceilingData = {
  {"both_10", exactDecimal["0.0028313316855593175"]},
  {"one_10", exactDecimal["0.0035965105896846573"]},
  {"both_30", exactDecimal["0.00817339430971383"]},
  {"one_30", exactDecimal["0.0116633929790174"]}
};

Do[
  Module[{label, pcritValue, epsXiBudget, aBudget},
    {label, pcritValue} = entry;
    epsXiBudget = pcritValue/barP0compat - 1;
    aBudget = (pcritValue - barP0compat)/4;
    Print[
      label,
      ": |eps xi1| <= ", fmt[N[epsXiBudget, 24]],
      " ; |aP| <= ", fmt[N[aBudget, 24]]
    ];
    expectZero[
      "M6 " <> label <> " epsXi defining relation",
      (epsXiBudget + 1) barP0compat - pcritValue
    ];
    expectZero[
      "M6 " <> label <> " aP defining relation",
      barP0compat + 4 aBudget - pcritValue
    ];
  ],
  {entry, ceilingData}
];

Print[""];
Print["All Stage 224 Mathematica checks passed."];
