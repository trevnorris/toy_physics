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
stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

cleanExpr[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Cancel[Together[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

zeroQ[expr_] := And @@ (TrueQ[# === 0] & /@ Flatten[{expr}]);

expectZero[name_String, expr_] := Module[{res},
  res = cleanExpr[expr];
  Print[name, " = ", fmt[res]];
  If[zeroQ[res], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

totalDegree[poly_, vars_List] := Module[{rules, powers},
  rules = CoefficientRules[Expand[poly], vars];
  If[Length[rules] == 0, Return[-Infinity]];
  powers = rules[[All, 1]];
  Max[Total /@ powers]
];

banner["STAGE 217 -- FIVE-COORDINATE SIMPLEX OPTIMIZER MATHEMATICA AUDIT"];

ClearAll[
  H0, r, s, t, u, y, kL, kc, kg, kU, kW, k, nu, nuD, nuO
];

coefNames = CharacterRange["A", "O"];
coefSymbols = AssociationThread[
  coefNames -> Quiet[(ToExpression["Global`" <> #] & /@ coefNames), General::shdw]
];
{coefA, coefB, coefC, coefD, coefE, coefF, coefG, coefH,
  coefI, coefJ, coefK, coefL, coefM, coefN, coefO} = coefSymbols /@ coefNames;
coefList = Values[coefSymbols];

ratioVars = {r, s, t, u};
liftVars = {r, s, t, u, y};
slopeVars = {kc, kg, kU, kW};

S = 1 + r^2 + s^2 + t^2 + u^2;
Klin = kL + r kc + s kg + t kU + u kW;
delta = (
  coefA + coefB r + coefC s + coefD t + coefE u
    + coefF r^2 + coefG r s + coefH r t + coefI r u
    + coefJ s^2 + coefK s t + coefL s u
    + coefM t^2 + coefN t u + coefO u^2
);
sqrtDelta = Sqrt[delta];
phi = (Klin + sqrtDelta)/Sqrt[S];

$Assumptions = (
  Element[
    Join[{H0, r, s, t, u, y, kL, kc, kg, kU, kW, k, nu, nuD, nuO}, coefList],
    Reals
  ]
  && H0 > 0
  && r > 0 && s > 0 && t > 0 && u > 0 && y > 0
  && kL > 0 && kc > 0 && kg > 0 && kU > 0 && kW > 0 && k > 0
  && delta > 0
);

metricNumerators = MapThread[
  FullSimplify[S #2 - #1 Klin, Assumptions -> $Assumptions] &,
  {ratioVars, slopeVars}
];
envelopeNumerators = FullSimplify[
    S D[delta, #] - D[S, #] delta,
    Assumptions -> $Assumptions
  ] & /@ ratioVars;

{Mr, Ms, Mt, Mu} = metricNumerators;
{Lr, Ls, Lt, Lu} = envelopeNumerators;

subbanner["M1. Exact stationary-numerator identities"];

stationaryNumerators = 2 sqrtDelta metricNumerators + envelopeNumerators;
scaledDerivativeNumerators = FullSimplify[
    Cancel[Together[2 sqrtDelta S^(3/2) D[phi, #]]],
    Assumptions -> $Assumptions
  ] & /@ ratioVars;

Print["M1 cleared derivative numerators = ", fmt[scaledDerivativeNumerators]];
Scan[
  Function[{row},
    expectZero[
      "M1 stationary numerator identity " <> row[[1]],
      row[[2]] - row[[3]]
    ]
  ],
  Transpose[{{"r", "s", "t", "u"}, scaledDerivativeNumerators, stationaryNumerators}]
];

subbanner["M2. Five codimension-one quadruple faces"];

axes = {"lambda", "c", "gamma", "U", "W"};
faces = Table[Delete[axes, i], {i, Length[axes]}];

Print["M2 primitive axes = ", fmt[axes]];
Print["M2 generated faces = ", fmt[faces]];
expectTrue["M2 face count equals five", Length[faces] == 5];
expectTrue["M2 every face has cardinality four", And @@ Thread[Length /@ faces == 4]];

subbanner["M3. Lifted degree pattern and Bezout product"];

liftedPolys = Expand[2 # y + #2] & @@@ Transpose[{metricNumerators, envelopeNumerators}];
deltaLift = y^2 - delta;
liftedDegrees = totalDegree[#, liftVars] & /@ Append[liftedPolys, deltaLift];
liftedProduct = Times @@ liftedDegrees;

Print["M3 lifted degrees {F_r,F_s,F_t,F_u,F_delta} = ", fmt[liftedDegrees]];
Print["M3 lifted Bezout product = ", liftedProduct];
Print["M3 literal 3^4*2 evaluates to ", 3^4*2];
expectTrue["M3 lifted degree tuple is {3,3,3,3,2}", liftedDegrees === {3, 3, 3, 3, 2}];
expectTrue["M3 lifted product equals 3^4*2", liftedProduct == 3^4*2];

subbanner["M4. Projected square-root-free degree pattern and bound"];

Crs = Expand[Ms Lr - Mr Ls];
Crt = Expand[Mt Lr - Mr Lt];
Cru = Expand[Mu Lr - Mr Lu];
Sr = Expand[Lr^2 - 4 Mr^2 delta];
projectedDegrees = totalDegree[#, ratioVars] & /@ {Crs, Crt, Cru, Sr};
projectedProduct = Times @@ projectedDegrees;

Print["M4 projected degrees {C_rs,C_rt,C_ru,S_r} = ", fmt[projectedDegrees]];
Print["M4 projected one-chart Bezout product = ", projectedProduct];
Print["M4 literal 5*5*5*6 evaluates to ", 5*5*5*6];
expectTrue["M4 projected degree tuple is {5,5,5,6}", projectedDegrees === {5, 5, 5, 6}];
expectTrue["M4 projected product equals 5*5*5*6", projectedProduct == 5*5*5*6];

subbanner["M5. Diagonal-isotropic gradient-optimal reduction"];

diagRules = {
  coefA -> kL^2 - 2 H0 nu,
  coefB -> 2 kL kc,
  coefC -> 2 kL kg,
  coefD -> 2 kL kU,
  coefE -> 2 kL kW,
  coefF -> kc^2 - 2 H0 nu,
  coefG -> 2 kc kg,
  coefH -> 2 kc kU,
  coefI -> 2 kc kW,
  coefJ -> kg^2 - 2 H0 nu,
  coefK -> 2 kg kU,
  coefL -> 2 kg kW,
  coefM -> kU^2 - 2 H0 nu,
  coefN -> 2 kU kW,
  coefO -> kW^2 - 2 H0 nu
};
gradientRatioRules = {r -> kc/kL, s -> kg/kL, t -> kU/kL, u -> kW/kL};
diagMetric = FullSimplify[metricNumerators /. diagRules, Assumptions -> $Assumptions];
diagEnvelope = FullSimplify[envelopeNumerators /. diagRules, Assumptions -> $Assumptions];
diagKlin = Klin;

Scan[
  Function[{row},
    expectZero[
      "M5 L_" <> row[[1]] <> "(diag) - 2 Klin M_" <> row[[1]],
      row[[3]] - 2 diagKlin row[[2]]
    ]
  ],
  Transpose[{{"r", "s", "t", "u"}, diagMetric, diagEnvelope}]
];

Scan[
  Function[{row},
    expectZero[
      "M5 M_" <> row[[1]] <> " at gradient ratios",
      row[[2]] /. gradientRatioRules
    ];
    expectZero[
      "M5 L_" <> row[[1]] <> " at gradient ratios",
      row[[3]] /. gradientRatioRules
    ];
  ],
  Transpose[{{"r", "s", "t", "u"}, diagMetric, diagEnvelope}]
];

subbanner["M6. Fivefold-symmetric equal-mix reduction"];

symRules = {
  kL -> k,
  kc -> k,
  kg -> k,
  kU -> k,
  kW -> k,
  coefA -> k^2 - 2 H0 nuD,
  coefB -> 2 k^2 - 4 H0 nuO,
  coefC -> 2 k^2 - 4 H0 nuO,
  coefD -> 2 k^2 - 4 H0 nuO,
  coefE -> 2 k^2 - 4 H0 nuO,
  coefF -> k^2 - 2 H0 nuD,
  coefG -> 2 k^2 - 4 H0 nuO,
  coefH -> 2 k^2 - 4 H0 nuO,
  coefI -> 2 k^2 - 4 H0 nuO,
  coefJ -> k^2 - 2 H0 nuD,
  coefK -> 2 k^2 - 4 H0 nuO,
  coefL -> 2 k^2 - 4 H0 nuO,
  coefM -> k^2 - 2 H0 nuD,
  coefN -> 2 k^2 - 4 H0 nuO,
  coefO -> k^2 - 2 H0 nuD
};
equalMixRules = {r -> 1, s -> 1, t -> 1, u -> 1};
symMetricAtEqualMix = FullSimplify[
  metricNumerators /. symRules /. equalMixRules,
  Assumptions -> $Assumptions
];
symEnvelopeAtEqualMix = FullSimplify[
  envelopeNumerators /. symRules /. equalMixRules,
  Assumptions -> $Assumptions
];

Scan[
  Function[{row},
    expectZero["M6 M_" <> row[[1]] <> " at equal mix", row[[2]]];
    expectZero["M6 L_" <> row[[1]] <> " at equal mix", row[[3]]];
  ],
  Transpose[{{"r", "s", "t", "u"}, symMetricAtEqualMix, symEnvelopeAtEqualMix}]
];

Print[""];
Print["All Stage 217 Mathematica audit checks passed."];
Exit[0];
