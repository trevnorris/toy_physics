(* Stage 209 Mathematica audit: pairwise ratio optimizer and mixed-ray winner theorem. *)
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
  Print[name, " = ", fmt[res]];
  If[zeroTensorQ[res], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 209 -- PAIRWISE RATIO OPTIMIZER AND MIXED-RAY WINNER THEOREM"];

Clear[ki, kj, H0, u, v, w, r, kappaStar, z];

$Assumptions = (
  Element[{u, v, w, kappaStar, z}, Reals] &&
  ki > 0 && kj > 0 && H0 > 0 && r > 0
);

den[x_] := 1 + x^2;
raySlope[x_] := (ki + x kj)/Sqrt[den[x]];
rayCurvature[x_] := (u + 2 v x + w x^2)/den[x];
discNumerator[x_] := Expand[
  Numerator[Together[den[x] (raySlope[x]^2 - 2 H0 rayCurvature[x])]]
];

discPoly = discNumerator[r];
discCoefficients = CoefficientList[discPoly, r];
abcByCoefficientExtraction = {
  discCoefficients[[1]],
  discCoefficients[[2]],
  discCoefficients[[3]]
};
abcClosedForm = {
  ki^2 - 2 H0 u,
  2 ki kj - 4 H0 v,
  kj^2 - 2 H0 w
};

discFromCoefficients[x_] := (
  abcByCoefficientExtraction[[1]]
  + abcByCoefficientExtraction[[2]] x
  + abcByCoefficientExtraction[[3]] x^2
);

closureRoot[x_] := 2 H0/(
  raySlope[x] + Sqrt[raySlope[x]^2 - 2 H0 rayCurvature[x]]
);
coefficientRoot[x_] := 2 H0 Sqrt[den[x]]/(
  ki + x kj + Sqrt[discFromCoefficients[x]]
);
denominatorFunctional[x_] := (
  ki + x kj + Sqrt[discFromCoefficients[x]]
)/Sqrt[den[x]];

subbanner["M1. Algebraic form of the certified root"];
expectZero[
  "M1 closure root minus coefficient-discriminant form",
  closureRoot[r] - coefficientRoot[r]
];

subbanner["M2. Discriminant numerator reduction"];
expectZero[
  "M2 discriminant coefficient extraction residuals",
  abcByCoefficientExtraction - abcClosedForm
];
expectZero[
  "M2 discriminant numerator polynomial residual",
  den[r] (raySlope[r]^2 - 2 H0 rayCurvature[r]) - discFromCoefficients[r]
];

subbanner["M3. Denominator-functional identity"];
expectZero[
  "M3 closure root minus 2H0/Phi",
  closureRoot[r] - 2 H0/denominatorFunctional[r]
];

subbanner["M4. Stationary numerator and derivative law"];
phiDerivative = D[denominatorFunctional[r], r];
polyPartFromDifferentiatedDiscriminant = FullSimplify[
  D[discFromCoefficients[r], r] den[r] - 2 r discFromCoefficients[r],
  Assumptions -> $Assumptions
];
stationaryNumeratorFromDerivative = FullSimplify[
  2 den[r]^(3/2) Sqrt[discFromCoefficients[r]] phiDerivative,
  Assumptions -> $Assumptions
];
manifestStationaryNumerator = (
  2 (kj - ki r) Sqrt[discFromCoefficients[r]]
  + abcByCoefficientExtraction[[2]]
  + 2 (abcByCoefficientExtraction[[3]] - abcByCoefficientExtraction[[1]]) r
  - abcByCoefficientExtraction[[2]] r^2
);

expectTrue[
  "M4 Phi derivative is not identically zero before reductions",
  D[denominatorFunctional[r], r] =!= 0
];
expectZero[
  "M4 scaled derivative numerator minus manifest N",
  stationaryNumeratorFromDerivative - manifestStationaryNumerator
];
expectZero[
  "M4 Phi derivative law residual",
  phiDerivative -
    manifestStationaryNumerator/(
      2 den[r]^(3/2) Sqrt[discFromCoefficients[r]]
    )
];

subbanner["M5. Quartic elimination degree and factorization"];
quarticByResultant = Factor[
  Resultant[
    polyPartFromDifferentiatedDiscriminant + 2 (kj - ki r) z,
    z^2 - discFromCoefficients[r],
    z
  ]
];

Print["M5 resultant quartic = ", fmt[quarticByResultant]];
expectTrue[
  "M5 resultant quartic degree is 4",
  Exponent[Expand[quarticByResultant], r] == 4
];
expectZero[
  "M5 resultant factorization identity",
  quarticByResultant -
    (polyPartFromDifferentiatedDiscriminant -
      2 (kj - ki r) Sqrt[discFromCoefficients[r]]) *
    (polyPartFromDifferentiatedDiscriminant +
      2 (kj - ki r) Sqrt[discFromCoefficients[r]])
];
expectZero[
  "M5 plus radical factor equals derivative numerator",
  stationaryNumeratorFromDerivative -
    (polyPartFromDifferentiatedDiscriminant +
      2 (kj - ki r) Sqrt[discFromCoefficients[r]])
];

subbanner["M6. Diagonal-neutral reduction"];
diagonalNeutralRules = {u -> kappaStar, v -> 0, w -> kappaStar};
gradientRatioFromSlopeSolve = r /. Solve[D[raySlope[r], r] == 0, r, Reals];
expectTrue[
  "M6 raw tau derivative is not identically zero before diagonal-neutral reduction",
  D[coefficientRoot[r], r] =!= 0
];
expectTrue[
  "M6 gradient ratio recovered from slope stationarity",
  MemberQ[
    FullSimplify[gradientRatioFromSlopeSolve, Assumptions -> $Assumptions],
    kj/ki
  ]
];
expectZero[
  "M6 diagonal-neutral curvature is constant",
  (rayCurvature[r] /. diagonalNeutralRules) - kappaStar
];
expectZero[
  "M6 gradient-optimal stationarity",
  D[coefficientRoot[r] /. diagonalNeutralRules, r] /. r -> kj/ki
];

subbanner["M7. Pair-symmetry reduction"];
pairSymmetryRules = {kj -> ki, w -> u};
pairSymmetricTau[x_] := coefficientRoot[x] /. pairSymmetryRules;
expectZero[
  "M7 reciprocal ray invariance",
  pairSymmetricTau[r] - pairSymmetricTau[1/r]
];
expectZero[
  "M7 differentiated reciprocal identity at equal mix",
  D[pairSymmetricTau[r] - pairSymmetricTau[1/r], r] /. r -> 1
];
expectZero[
  "M7 equal-mix stationarity",
  D[pairSymmetricTau[r], r] /. r -> 1
];

banner["STAGE 209 MATHEMATICA AUDIT PASSED"];
Exit[0];
