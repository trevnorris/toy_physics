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

cleanScalar[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[res], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

cleanTensor[expr_] := If[
  ListQ[expr],
  Map[cleanScalar, expr, {ArrayDepth[expr]}],
  cleanScalar[expr]
];

zeroTensorQ[expr_] := And @@ (TrueQ[# === 0] & /@ Flatten[{expr}]);

expectZero[name_String, expr_] := Module[{res},
  res = cleanTensor[expr];
  If[ListQ[res],
    Print[name, " = ", fmt[res]],
    Print[name, " = ", fmt[res]]
  ];
  If[zeroTensorQ[res], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 210 -- THREE-COORDINATE MIXED-SIMPLEX MATHEMATICA AUDIT"];

Clear[
  ai, aj, ak, r, s, nu, eps, x, y, z, lag,
  ki, kj, kk, H0, uii, uij, uik, ujj, ujk, ukk
];

$Assumptions = (
  Element[
    {
      ai, aj, ak, r, s, nu, x, y, z, lag,
      ki, kj, kk, H0, uii, uij, uik, ujj, ujk, ukk
    },
    Reals
  ] &&
  ki > 0 && kj > 0 && kk > 0 &&
  ai >= 0 && aj >= 0 && ak >= 0 &&
  H0 > 0 &&
  r >= 0 && s >= 0 && nu >= 0
);

aVec = {ai, aj, ak};
kVec = {ki, kj, kk};
kSimplex = aVec . kVec;

normalizeRaw[raw_List] := FullSimplify[raw/Sqrt[raw . raw], Assumptions -> $Assumptions];
slopeOf[vec_List] := FullSimplify[vec . kVec, Assumptions -> $Assumptions];

edgeIJ = normalizeRaw[{1, r, 0}];
edgeIK = normalizeRaw[{1, 0, s}];
edgeJK = normalizeRaw[{0, 1, nu}];

subbanner["M1. Edge normalization"];

expectZero["M1 edge ij normalization", edgeIJ . edgeIJ - 1];
expectZero["M1 edge ik normalization", edgeIK . edgeIK - 1];
expectZero["M1 edge jk normalization", edgeJK . edgeJK - 1];

subbanner["M2. Boundary slope reductions from limiting edge patches"];

slopeIJLimit = FullSimplify[
  Limit[slopeOf[normalizeRaw[{1, r, eps}]], eps -> 0, Direction -> "FromAbove"],
  Assumptions -> $Assumptions
];
slopeIKLimit = FullSimplify[
  Limit[slopeOf[normalizeRaw[{1, eps, s}]], eps -> 0, Direction -> "FromAbove"],
  Assumptions -> $Assumptions
];
slopeJKLimit = FullSimplify[
  Limit[slopeOf[normalizeRaw[{eps, 1, nu}]], eps -> 0, Direction -> "FromAbove"],
  Assumptions -> $Assumptions
];

expectZero["M2 ij slope boundary residual", slopeIJLimit - (ki + r kj)/Sqrt[1 + r^2]];
expectZero["M2 ik slope boundary residual", slopeIKLimit - (ki + s kk)/Sqrt[1 + s^2]];
expectZero["M2 jk slope boundary residual", slopeJKLimit - (kj + nu kk)/Sqrt[1 + nu^2]];

subbanner["M3. Gradient-optimal point from Lagrange stationarity"];

kNorm = Sqrt[kVec . kVec];
lagrangian = (
  x ki + y kj + z kk
  - lag (x^2 + y^2 + z^2 - 1)
);
stationarityEquations = {
  D[lagrangian, x] == 0,
  D[lagrangian, y] == 0,
  D[lagrangian, z] == 0,
  x^2 + y^2 + z^2 == 1
};
stationaritySolutions = Solve[stationarityEquations, {x, y, z, lag}, Reals];
positiveStationary = SelectFirst[
  stationaritySolutions,
  TrueQ[
    FullSimplify[
      (x /. #) > 0 && (y /. #) > 0 && (z /. #) > 0,
      Assumptions -> $Assumptions
    ]
  ] &
];
If[MissingQ[positiveStationary], fail["M3 positive Lagrange branch", stationaritySolutions]];

stationaryVec = FullSimplify[{x, y, z} /. positiveStationary, Assumptions -> $Assumptions];
gradVec = kVec/kNorm;
stationarityResiduals = {
  ki - 2 lag x,
  kj - 2 lag y,
  kk - 2 lag z,
  x^2 + y^2 + z^2 - 1
} /. positiveStationary;
cauchyResidual = (
  (kVec . kVec) (aVec . aVec) - (aVec . kVec)^2
  - ((ki aj - kj ai)^2 + (ki ak - kk ai)^2 + (kj ak - kk aj)^2)
);

Print["M3 positive Lagrange branch = ", fmt[stationaryVec]];
expectZero["M3 Lagrange stationarity residuals", stationarityResiduals];
expectZero["M3 stationary branch minus gradient vector", stationaryVec - gradVec];
expectZero["M3 gradient normalization", gradVec . gradVec - 1];
expectZero["M3 Cauchy identity residual", cauchyResidual];

subbanner["M4. Maximum slope value and interior ratios"];

expectZero["M4 maximum slope residual", stationaryVec . kVec - kNorm];
expectZero["M4 ratio aj/ai residual", stationaryVec[[2]]/stationaryVec[[1]] - kj/ki];
expectZero["M4 ratio ak/ai residual", stationaryVec[[3]]/stationaryVec[[1]] - kk/ki];

subbanner["M5. Strict interior gain via Pythagorean decompositions"];

expectZero["M5 Kijk^2 - Kij^2 - kk^2", (ki^2 + kj^2 + kk^2) - (ki^2 + kj^2) - kk^2];
expectZero["M5 Kijk^2 - Kik^2 - kj^2", (ki^2 + kj^2 + kk^2) - (ki^2 + kk^2) - kj^2];
expectZero["M5 Kijk^2 - Kjk^2 - ki^2", (ki^2 + kj^2 + kk^2) - (kj^2 + kk^2) - ki^2];

subbanner["M6. Cross-leverage identity, Cauchy slack, and screen values"];

wSigma = 2 (ai aj + ai ak + aj ak);
sumA = Total[aVec];
normA2 = aVec . aVec;
cauchySlack = (
  (ai - aj)^2 + (ai - ak)^2 + (aj - ak)^2
);
equalMix = {1, 1, 1}/Sqrt[3];
pairMix = {1, 1, 0}/Sqrt[2];
uniqueEqualMixCondition = FullSimplify[
  Reduce[
    ai >= 0 && aj >= 0 && ak >= 0 &&
      normA2 == 1 && wSigma == 2,
    {ai, aj, ak},
    Reals
  ],
  Assumptions -> $Assumptions
];

expectZero["M6 wSigma sum-square identity", wSigma - (sumA^2 - normA2)];
expectZero["M6 Cauchy slack identity", 3 normA2 - sumA^2 - cauchySlack];
expectZero["M6 slack bound identity on unit sphere", cauchySlack - (2 normA2 - wSigma)];
expectZero["M6 equal-mix screen value", (wSigma /. Thread[aVec -> equalMix]) - 2];
expectZero["M6 pairwise equal-edge screen value", (wSigma /. Thread[aVec -> pairMix]) - 1];
Print["M6 wSigma == 2 unit-simplex condition = ", fmt[uniqueEqualMixCondition]];
expectTrue[
  "M6 equal-mix uniqueness condition",
  Equivalent[
    uniqueEqualMixCondition,
    ai == 1/Sqrt[3] && aj == 1/Sqrt[3] && ak == 1/Sqrt[3]
  ]
];

subbanner["M7. Curvature law and discriminant numerator by coefficient extraction"];

hMat = {
  {uii, uij, uik},
  {uij, ujj, ujk},
  {uik, ujk, ukk}
};
kappaOf[vec_List] := FullSimplify[vec . hMat . vec, Assumptions -> $Assumptions];

denRS = 1 + r^2 + s^2;
ratioVec = {1, r, s}/Sqrt[denRS];
slopeRS = slopeOf[ratioVec];
kappaRS = kappaOf[ratioVec];

coefA = ki^2 - 2 H0 uii;
coefB = 2 ki kj - 4 H0 uij;
coefC = 2 ki kk - 4 H0 uik;
coefD = kj^2 - 2 H0 ujj;
coefE = 2 kj kk - 4 H0 ujk;
coefF = kk^2 - 2 H0 ukk;
deltaSharp = coefA + coefB r + coefC s + coefD r^2 + coefE r s + coefF s^2;

deltaCleared = Expand[
  FullSimplify[
    denRS (slopeRS^2 - 2 H0 kappaRS),
    Assumptions -> $Assumptions
  ]
];
deltaSeries = Expand[Normal[Series[deltaCleared, {r, 0, 2}, {s, 0, 2}]]];
deltaCoefficientList = CoefficientList[deltaSeries, {r, s}];
coeffRS[i_Integer, j_Integer] := Coefficient[Coefficient[deltaSeries, r, i], s, j];

Print["M7 CoefficientList in {r,s} = ", fmt[deltaCoefficientList]];
expectZero["M7 coefficient A residual", coeffRS[0, 0] - coefA];
expectZero["M7 coefficient B residual", coeffRS[1, 0] - coefB];
expectZero["M7 coefficient C residual", coeffRS[0, 1] - coefC];
expectZero["M7 coefficient D residual", coeffRS[2, 0] - coefD];
expectZero["M7 coefficient E residual", coeffRS[1, 1] - coefE];
expectZero["M7 coefficient F residual", coeffRS[0, 2] - coefF];

subbanner["M8. Curvature edge reductions and diagonal-neutral case"];

kappaGeneral = kappaOf[aVec];
kappaIJLimit = FullSimplify[
  Limit[kappaOf[normalizeRaw[{1, r, eps}]], eps -> 0, Direction -> "FromAbove"],
  Assumptions -> $Assumptions
];
kappaIKLimit = FullSimplify[
  Limit[kappaOf[normalizeRaw[{1, eps, s}]], eps -> 0, Direction -> "FromAbove"],
  Assumptions -> $Assumptions
];
kappaJKLimit = FullSimplify[
  Limit[kappaOf[normalizeRaw[{eps, 1, nu}]], eps -> 0, Direction -> "FromAbove"],
  Assumptions -> $Assumptions
];

expectZero["M8 ij curvature edge residual", kappaIJLimit - (uii + 2 uij r + ujj r^2)/(1 + r^2)];
expectZero["M8 ik curvature edge residual", kappaIKLimit - (uii + 2 uik s + ukk s^2)/(1 + s^2)];
expectZero["M8 jk curvature edge residual", kappaJKLimit - (ujj + 2 ujk nu + ukk nu^2)/(1 + nu^2)];
expectZero[
  "M8 diagonal-neutral curvature residual",
  (kappaGeneral /. {uij -> 0, uik -> 0, ujk -> 0}) -
    (uii ai^2 + ujj aj^2 + ukk ak^2)
];

subbanner["M9. Fixed-simplex root map and boundary reductions"];

tauRatioDirect = FullSimplify[
  2 H0/(slopeRS + Sqrt[slopeRS^2 - 2 H0 kappaRS]),
  Assumptions -> $Assumptions
];
tauRatioExpected = (
  2 H0 Sqrt[denRS]/
    (ki + r kj + s kk + Sqrt[deltaSharp])
);

tauIJExpected = (
  2 H0 Sqrt[1 + r^2]/
    (ki + r kj + Sqrt[
      ki^2 - 2 H0 uii +
        (2 ki kj - 4 H0 uij) r +
        (kj^2 - 2 H0 ujj) r^2
    ])
);
tauIKExpected = (
  2 H0 Sqrt[1 + s^2]/
    (ki + s kk + Sqrt[
      ki^2 - 2 H0 uii +
        (2 ki kk - 4 H0 uik) s +
        (kk^2 - 2 H0 ukk) s^2
    ])
);
tauJKDirect = FullSimplify[
  2 H0/(slopeJKLimit + Sqrt[slopeJKLimit^2 - 2 H0 kappaJKLimit]),
  Assumptions -> $Assumptions
];
tauJKExpected = (
  2 H0 Sqrt[1 + nu^2]/
    (kj + nu kk + Sqrt[
      kj^2 - 2 H0 ujj +
        (2 kj kk - 4 H0 ujk) nu +
        (kk^2 - 2 H0 ukk) nu^2
    ])
);

expectZero["M9 ratio-coordinate tau residual", tauRatioDirect - tauRatioExpected];
expectZero["M9 ij boundary tau residual", (tauRatioExpected /. s -> 0) - tauIJExpected];
expectZero["M9 ik boundary tau residual", (tauRatioExpected /. r -> 0) - tauIKExpected];
expectZero["M9 jk boundary tau residual", tauJKDirect - tauJKExpected];

banner["STAGE 210 MATHEMATICA AUDIT PASSED"];
Exit[0];
