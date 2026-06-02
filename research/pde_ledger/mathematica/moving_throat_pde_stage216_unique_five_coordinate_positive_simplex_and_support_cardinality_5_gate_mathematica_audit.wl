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

cleanExpr[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

zeroQ[expr_] := If[ListQ[expr], And @@ (zeroQ /@ Flatten[expr]), TrueQ[expr === 0]];

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

banner["STAGE 216 -- UNIQUE FIVE-COORDINATE POSITIVE SIMPLEX AND SUPPORT-CARDINALITY-5 GATE"];

ClearAll[
  kLam, kC, kGam, kU, kW, a1, a2, a3, a4, a5,
  mu, H0, k, kappa, tau
];

kVec = {kLam, kC, kGam, kU, kW};
aVec = {a1, a2, a3, a4, a5};
axisNames = {"lambda", "c", "gamma", "U", "W"};

$Assumptions = (
  Element[Join[kVec, aVec, {mu, H0, k, kappa, tau}], Reals]
  && And @@ Thread[kVec > 0]
  && H0 > 0 && k > 0 && kappa > 0
  && k^2 > 2 H0 kappa
);

S = kVec.kVec;

subbanner["M1-M2. Lagrange-derived five-coordinate optimum"];

lagrangian = kVec.aVec - mu (aVec.aVec - 1);
stationarity = Thread[(D[lagrangian, #] & /@ aVec) == 0];
lagrangeRules = Solve[Join[stationarity, {aVec.aVec == 1}], Append[aVec, mu], Reals];
positiveRule = SelectFirst[
  lagrangeRules,
  TrueQ[FullSimplify[(mu /. #) > 0, Assumptions -> $Assumptions]] &
];
If[MissingQ[positiveRule], fail["M1-M2 positive Lagrange branch was not found", lagrangeRules]];

lagrangePoint = aVec /. positiveRule;
lagrangeValue = kVec.aVec /. positiveRule;
stationaryValues = FullSimplify[kVec.aVec /. lagrangeRules, Assumptions -> $Assumptions];

Print["positive Lagrange point = ", fmt[lagrangePoint]];
Print["stationary slope values = ", fmt[stationaryValues]];
expectZero["M1 gradient-optimal unit norm", lagrangePoint.lagrangePoint - 1];
expectZero["M2 positive Lagrange value minus Sqrt[S]", lagrangeValue - Sqrt[S]];
expectTrue[
  "M2 positive branch dominates the other stationary branch",
  And @@ (FullSimplify[lagrangeValue >= #, Assumptions -> $Assumptions] & /@ stationaryValues)
];

subbanner["M3. Per-face first-order gaps"];

faceRows = Table[
  {axisNames[[idx]], kVec[[idx]], Total[Delete[kVec^2, idx]]},
  {idx, 1, Length[kVec]}
];

Scan[
  Function[row,
    expectZero["M3 face gap for " <> row[[1]], S - row[[3]] - row[[2]]^2];
    expectTrue["M3 strict positive improvement for " <> row[[1]], row[[2]]^2 > 0];
  ],
  faceRows
];

subbanner["M4. Cross-leverage and Cauchy slack as quadratic forms"];

ones = ConstantArray[1, Length[aVec]];
offDiagonalForm = ConstantArray[1, {Length[aVec], Length[aVec]}] - IdentityMatrix[Length[aVec]];
laplacianForm = Length[aVec] IdentityMatrix[Length[aVec]] - Outer[Times, ones, ones];
wSigma = aVec.offDiagonalForm.aVec;
pairDifferenceEnergy = Total[(Subtract @@ #)^2 & /@ Subsets[aVec, {2}]];

expectZero[
  "M4 off-diagonal quadratic identity",
  wSigma - ((ones.aVec)^2 - aVec.aVec)
];
expectZero[
  "M4 complete-graph Laplacian slack identity",
  aVec.laplacianForm.aVec - pairDifferenceEnergy
];

subbanner["M5. Equal-mix barycenter leverage"];

barycenter = Normalize[ones];
leverageEigenvalues = Eigenvalues[offDiagonalForm];
Print["off-diagonal leverage eigenvalues = ", fmt[leverageEigenvalues]];
expectZero["M5 barycenter unit norm", barycenter.barycenter - 1];
expectZero["M5 wSigma at barycenter minus 4", barycenter.offDiagonalForm.barycenter - 4];
expectZero["M5 top quadratic-form eigenvalue minus 4", Max[leverageEigenvalues] - 4];

subbanner["M6. Quadratic-solve certified bracket"];

quadraticResidual = 1/2 kappa tau^2 - k tau + H0;
rootRules = Solve[quadraticResidual == 0, tau, Reals];
roots = FullSimplify[stripConditional[tau /. rootRules], Assumptions -> $Assumptions];
statedTau = 2 H0/(k + Sqrt[k^2 - 2 H0 kappa]);
solvedTau = SelectFirst[
  roots,
  TrueQ[FullSimplify[# == statedTau, Assumptions -> $Assumptions]] &
];
If[MissingQ[solvedTau], fail["M6 stated bracket root was not returned by Solve", roots]];
otherTau = SelectFirst[
  roots,
  !TrueQ[FullSimplify[# == solvedTau, Assumptions -> $Assumptions]] &
];
If[MissingQ[otherTau], fail["M6 second quadratic root was not returned by Solve", roots]];

Print["quadratic roots from Solve = ", fmt[roots]];
expectZero["M6 solved smaller root minus stated bracket form", solvedTau - statedTau];
expectZero["M6 stated bracket quadratic residual", quadraticResidual /. tau -> statedTau];
expectTrue["M6 stated bracket root is positive", statedTau > 0];
expectTrue["M6 stated bracket root is the smaller positive root", statedTau < otherTau];

Print[""];
Print["All Stage 216 Mathematica audit checks passed."];
Exit[0];
