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

stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

normalizeExpr[expr_, assumptions_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[res], Assumptions -> assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> assumptions]
];

zeroQ[expr_] := If[ListQ[expr], And @@ (zeroQ /@ Flatten[expr]), TrueQ[expr == 0]];

expectZeroUnder[name_String, expr_, assumptions_] := Module[{res},
  res = normalizeExpr[expr, assumptions];
  Print[name, " = ", fmt[res]];
  If[zeroQ[res], pass[name], fail[name, res]];
];

expectZero[name_String, expr_] := expectZeroUnder[name, expr, $Assumptions];

expectTrue[name_String, cond_] := (
  Print[name, " = ", cond];
  If[TrueQ[cond], pass[name], fail[name, cond]];
);

coeff3[poly_, powers_List] :=
  Coefficient[
    Coefficient[
      Coefficient[poly, r, powers[[1]]],
      s,
      powers[[2]]
    ],
    t,
    powers[[3]]
  ];

banner["STAGE 213 -- FOUR-COORDINATE MIXED SIMPLEX AND SUPPORT-CARDINALITY-4 GATE"];

ClearAll[
  ai, aj, ak, al, r, s, t, u, v, ki, kj, kk, kl, H0,
  uii, uij, uik, uil, ujj, ujk, ujl, ukk, ukl, ull,
  x1, x2, x3, x4, mu
];

$Assumptions = (
  Element[
    {ai, aj, ak, al, r, s, t, u, v, ki, kj, kk, kl, H0,
      uii, uij, uik, uil, ujj, ujk, ujl, ukk, ukl, ull},
    Reals
  ]
  && ai >= 0 && aj >= 0 && ak >= 0 && al >= 0
  && r >= 0 && s >= 0 && t >= 0 && u >= 0 && v >= 0
  && ki > 0 && kj > 0 && kk > 0 && kl > 0 && H0 > 0
);

kVec = {ki, kj, kk, kl};
aVec = {ai, aj, ak, al};
kSimplex[vec_List] := kVec.vec;

Hblock = {
  {uii, uij, uik, uil},
  {uij, ujj, ujk, ujl},
  {uik, ujk, ukk, ukl},
  {uil, ujl, ukl, ull}
};
kappa[vec_List] := vec.Hblock.vec;

subbanner["M1. Combinatorial ledger"];

axes = {"lambda", "c", "gamma", "U", "W"};
triples = Subsets[axes, {3}];
quadruples = Subsets[axes, {4}];
tripleIncidences = Association @ Table[
  tri -> Count[quadruples, quad_ /; SubsetQ[quad, tri]],
  {tri, triples}
];
axisIncidences = Association @ Table[
  axis -> Count[quadruples, quad_ /; MemberQ[quad, axis]],
  {axis, axes}
];

Print["M1 primitive axes = ", axes];
Print["M1 primitive quadruples = ", quadruples];
Print["M1 triple incidences = ", tripleIncidences];
Print["M1 axis incidences = ", axisIncidences];
expectZero["M1 Binomial[5,4] quadruple count", Length[quadruples] - Binomial[5, 4]];
expectTrue["M1 every triple is contained in exactly two quadruples", And @@ Thread[Values[tripleIncidences] == 2]];
expectTrue["M1 every axis appears in exactly four quadruples", And @@ Thread[Values[axisIncidences] == 4]];

subbanner["M2. Face slope reductions"];

faceData = {
  <|
    "name" -> "ijk",
    "vector" -> {1, r, s, 0}/Sqrt[1 + r^2 + s^2],
    "slope" -> (ki + r kj + s kk)/Sqrt[1 + r^2 + s^2]
  |>,
  <|
    "name" -> "ijl",
    "vector" -> {1, r, 0, t}/Sqrt[1 + r^2 + t^2],
    "slope" -> (ki + r kj + t kl)/Sqrt[1 + r^2 + t^2]
  |>,
  <|
    "name" -> "ikl",
    "vector" -> {1, 0, s, t}/Sqrt[1 + s^2 + t^2],
    "slope" -> (ki + s kk + t kl)/Sqrt[1 + s^2 + t^2]
  |>,
  <|
    "name" -> "jkl",
    "vector" -> {0, 1, u, v}/Sqrt[1 + u^2 + v^2],
    "slope" -> (kj + u kk + v kl)/Sqrt[1 + u^2 + v^2]
  |>
};

Scan[
  Function[face,
    expectZero["M2 face " <> face["name"] <> " unit normalization", face["vector"].face["vector"] - 1];
    expectZero["M2 face " <> face["name"] <> " slope reduction", kSimplex[face["vector"]] - face["slope"]];
  ],
  faceData
];

subbanner["M3. Gradient-optimal ray by Lagrange equations"];

lagVars = {x1, x2, x3, x4};
lagObjective = kVec.lagVars;
lagNorm = lagVars.lagVars;
lagGradient = D[lagObjective - mu (lagNorm - 1), #] & /@ lagVars;
lagSystem = Join[
  Thread[lagGradient == 0],
  {lagNorm == 1}
];
lagSolutions = Solve[lagSystem, Append[lagVars, mu], Reals];
positiveBranch = SelectFirst[
  lagSolutions,
  TrueQ[FullSimplify[(mu /. #) > 0, Assumptions -> $Assumptions]] &
];
If[MissingQ[positiveBranch], fail["M3 positive Lagrange branch was not found", lagSolutions]];

kNormSq = kVec.kVec;
kNorm = Sqrt[kNormSq];
lagPoint = lagVars /. positiveBranch;
lagValue = lagObjective /. positiveBranch;
activeSupportValues = Association @ Table[
  support -> Sqrt[Total[kVec[[support]]^2]],
  {support, Rest[Subsets[Range[4]]]}
];

Print["M3 positive Lagrange point = ", fmt[lagPoint]];
Print["M3 active-support Lagrange values = ", fmt[activeSupportValues]];
expectZero["M3 positive Lagrange normalization", lagPoint.lagPoint - 1];
expectZero["M3 maximum value from positive branch minus ||k||", lagValue - kNorm];
expectZero["M3 ratio a_j/a_i - k_j/k_i", lagPoint[[2]]/lagPoint[[1]] - kj/ki];
expectZero["M3 ratio a_k/a_i - k_k/k_i", lagPoint[[3]]/lagPoint[[1]] - kk/ki];
expectZero["M3 ratio a_l/a_i - k_l/k_i", lagPoint[[4]]/lagPoint[[1]] - kl/ki];
Scan[
  Function[support,
    expectTrue[
      "M3 full-support value dominates active support " <> fmt[support],
      FullSimplify[kNormSq - activeSupportValues[support]^2 >= 0, Assumptions -> $Assumptions]
    ];
  ],
  Keys[activeSupportValues]
];

subbanner["M4. Synergy gaps"];

faceNormSq = <|
  "ijk" -> ki^2 + kj^2 + kk^2,
  "ijl" -> ki^2 + kj^2 + kl^2,
  "ikl" -> ki^2 + kk^2 + kl^2,
  "jkl" -> kj^2 + kk^2 + kl^2
|>;
expectZero["M4 ||k||^2 - ||k_ijk||^2 - k_l^2", kNormSq - faceNormSq["ijk"] - kl^2];
expectZero["M4 ||k||^2 - ||k_ijl||^2 - k_k^2", kNormSq - faceNormSq["ijl"] - kk^2];
expectZero["M4 ||k||^2 - ||k_ikl||^2 - k_j^2", kNormSq - faceNormSq["ikl"] - kj^2];
expectZero["M4 ||k||^2 - ||k_jkl||^2 - k_i^2", kNormSq - faceNormSq["jkl"] - ki^2];

subbanner["M5. Cross-leverage identity and bound"];

wSigma[vec_List] := 2 Total[Times @@@ Subsets[vec, {2}]];
sumA = Total[aVec];
normASq = aVec.aVec;
cauchySlack = Total[(Subtract @@ #)^2 & /@ Subsets[aVec, {2}]];
maxLeverage = Maximize[
  {wSigma[lagVars], lagVars.lagVars == 1 && And @@ Thread[lagVars >= 0]},
  lagVars,
  Reals
];

Print["M5 Maximize result = ", fmt[maxLeverage]];
expectZero["M5 wSigma identity", wSigma[aVec] - (sumA^2 - normASq)];
expectZero["M5 Cauchy slack identity", 4 normASq - sumA^2 - cauchySlack];
expectZero["M5 constrained maximum value - 3", First[maxLeverage] - 3];
expectZero["M5 maximizer equals four-way equal mix", (lagVars /. Last[maxLeverage]) - {1/2, 1/2, 1/2, 1/2}];
expectZero["M5 wSigma triple-equal point - 2", wSigma[{1, 1, 1, 0}/Sqrt[3]] - 2];
expectZero["M5 wSigma pair-equal edge - 1", wSigma[{1, 1, 0, 0}/Sqrt[2]] - 1];

subbanner["M6. Curvature quadratic form and diagonal-neutral reduction"];

kappaGeneral = Expand[kappa[aVec]];
kappaDiagonal = kappaGeneral /. {uij -> 0, uik -> 0, uil -> 0, ujk -> 0, ujl -> 0, ukl -> 0};
expectZero["M6 diagonal-neutral curvature reduction", kappaDiagonal - (uii ai^2 + ujj aj^2 + ukk ak^2 + ull al^2)];

subbanner["M7. Discriminant numerator by coefficient extraction"];

normRst = 1 + r^2 + s^2 + t^2;
zRst = {1, r, s, t};
aRst = zRst/Sqrt[normRst];
kRst = kSimplex[aRst];
kappaRst = kappa[aRst];
discNumerator = normalizeExpr[normRst (kRst^2 - 2 H0 kappaRst), $Assumptions];
monomialPowers = <|
  "A" -> {0, 0, 0},
  "B" -> {1, 0, 0},
  "C" -> {0, 1, 0},
  "D" -> {0, 0, 1},
  "E" -> {2, 0, 0},
  "F" -> {1, 1, 0},
  "G" -> {1, 0, 1},
  "Hhat" -> {0, 2, 0},
  "I" -> {0, 1, 1},
  "J" -> {0, 0, 2}
|>;
monomialFromPowers[powers_List] := r^powers[[1]] s^powers[[2]] t^powers[[3]];
extractedCoeffs = Association @ KeyValueMap[
  #1 -> normalizeExpr[coeff3[discNumerator, #2], $Assumptions] &,
  monomialPowers
];
quadraticKernel = Outer[Times, kVec, kVec] - 2 H0 Hblock;
matrixCoeffs = <|
  "A" -> quadraticKernel[[1, 1]],
  "B" -> 2 quadraticKernel[[1, 2]],
  "C" -> 2 quadraticKernel[[1, 3]],
  "D" -> 2 quadraticKernel[[1, 4]],
  "E" -> quadraticKernel[[2, 2]],
  "F" -> 2 quadraticKernel[[2, 3]],
  "G" -> 2 quadraticKernel[[2, 4]],
  "Hhat" -> quadraticKernel[[3, 3]],
  "I" -> 2 quadraticKernel[[3, 4]],
  "J" -> quadraticKernel[[4, 4]]
|>;
deltaSharp = Expand[
  Total[
    KeyValueMap[
      extractedCoeffs[#1] monomialFromPowers[#2] &,
      monomialPowers
    ]
  ]
];

Print["M7 extracted coefficient set = ", fmt[extractedCoeffs]];
Scan[
  Function[name,
    expectZero["M7 coefficient " <> name <> " matches boxed set", extractedCoeffs[name] - matrixCoeffs[name]];
  ],
  Keys[monomialPowers]
];
expectZero["M7 reconstructed extracted Delta sharp equals numerator", discNumerator - deltaSharp];

subbanner["M8. Certified tau bracket and ratio form"];

linearRst = kVec.zRst;
tauPatch = 2 H0/(kRst + Sqrt[kRst^2 - 2 H0 kappaRst]);
tauRatio = 2 H0 Sqrt[normRst]/(linearRst + Sqrt[deltaSharp]);
expectZeroUnder[
  "M8 tau ratio-form residual",
  tauPatch - tauRatio,
  $Assumptions && deltaSharp > 0
];

subbanner["M9. Face collapses of the bracket"];

faceDeltaIJK = Expand[
  extractedCoeffs["A"] + extractedCoeffs["B"] r + extractedCoeffs["C"] s
    + extractedCoeffs["E"] r^2 + extractedCoeffs["F"] r s + extractedCoeffs["Hhat"] s^2
];
faceDeltaIJL = Expand[
  extractedCoeffs["A"] + extractedCoeffs["B"] r + extractedCoeffs["D"] t
    + extractedCoeffs["E"] r^2 + extractedCoeffs["G"] r t + extractedCoeffs["J"] t^2
];
faceDeltaIKL = Expand[
  extractedCoeffs["A"] + extractedCoeffs["C"] s + extractedCoeffs["D"] t
    + extractedCoeffs["Hhat"] s^2 + extractedCoeffs["I"] s t + extractedCoeffs["J"] t^2
];

expectZero["M9 t=0 Delta sharp face polynomial", (deltaSharp /. t -> 0) - faceDeltaIJK];
expectZero["M9 s=0 Delta sharp face polynomial", (deltaSharp /. s -> 0) - faceDeltaIJL];
expectZero["M9 r=0 Delta sharp face polynomial", (deltaSharp /. r -> 0) - faceDeltaIKL];

expectZeroUnder[
  "M9 tau t=0 ijk face collapse",
  (tauRatio /. t -> 0)
    - 2 H0 Sqrt[1 + r^2 + s^2]/(ki + r kj + s kk + Sqrt[faceDeltaIJK]),
  $Assumptions && faceDeltaIJK > 0
];
expectZeroUnder[
  "M9 tau s=0 ijl face collapse",
  (tauRatio /. s -> 0)
    - 2 H0 Sqrt[1 + r^2 + t^2]/(ki + r kj + t kl + Sqrt[faceDeltaIJL]),
  $Assumptions && faceDeltaIJL > 0
];
expectZeroUnder[
  "M9 tau r=0 ikl face collapse",
  (tauRatio /. r -> 0)
    - 2 H0 Sqrt[1 + s^2 + t^2]/(ki + s kk + t kl + Sqrt[faceDeltaIKL]),
  $Assumptions && faceDeltaIKL > 0
];

normUV = 1 + u^2 + v^2;
aUV = {0, 1, u, v}/Sqrt[normUV];
kUV = kSimplex[aUV];
kappaUV = kappa[aUV];
faceDeltaJKL = Expand[
  extractedCoeffs["E"] + extractedCoeffs["F"] u + extractedCoeffs["G"] v
    + extractedCoeffs["Hhat"] u^2 + extractedCoeffs["I"] u v + extractedCoeffs["J"] v^2
];
discJKL = normalizeExpr[normUV (kUV^2 - 2 H0 kappaUV), $Assumptions];
tauJKL = 2 H0/(kUV + Sqrt[kUV^2 - 2 H0 kappaUV]);

expectZero["M9 jkl Delta sharp face polynomial", discJKL - faceDeltaJKL];
expectZeroUnder[
  "M9 direct jkl tau face collapse",
  tauJKL - 2 H0 Sqrt[normUV]/(kj + u kk + v kl + Sqrt[faceDeltaJKL]),
  $Assumptions && faceDeltaJKL > 0
];

Print["All Stage 213 Mathematica audit checks passed."];
