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

stripCE[expr_] := expr /. ConditionalExpression[e_, _] :> e;

cleanScalar[expr_] := Module[{res},
  res = stripCE[expr];
  res = FullSimplify[Together[res], Assumptions -> $Assumptions];
  res = stripCE[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

cleanTensor[expr_] := If[
  ListQ[expr],
  Map[cleanScalar, expr, {ArrayDepth[expr]}],
  cleanScalar[expr]
];

zeroTensorQ[expr_] := If[
  ListQ[expr],
  And @@ (TrueQ[# === 0] & /@ Flatten[expr]),
  TrueQ[expr === 0]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanTensor[expr];
  Print[name, " = ", fmt[res]];
  If[zeroTensorQ[res], pass[name], fail[name, res]];
];

expectTrue[name_String, claim_] := Module[{res},
  res = FullSimplify[stripCE[claim], Assumptions -> $Assumptions];
  res = stripCE[res];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 208 -- MATHEMATICA PAIRWISE MIXED-RAY AUDIT"];

Clear[
  ki, kj, r, hii, hij, hjj,
  hiiLo, hijLo, hjjLo, hiiHi, hijHi, hjjHi, H0
];

$Assumptions =
  Element[
    {
      ki, kj, r, hii, hij, hjj,
      hiiLo, hijLo, hjjLo, hiiHi, hijHi, hjjHi, H0
    },
    Reals
  ] &&
  ki > 0 && kj > 0 && r >= 0 && H0 > 0;

sRay[x_] := {1, x}/Sqrt[1 + x^2];
gradient = {-ki, -kj};
hessianPair = {{hii, hij}, {hij, hjj}};

orientedSlope = cleanScalar[gradient . sRay[r]];
positiveSlope = cleanScalar[-orientedSlope];

subbanner["M1-M3. Pairwise slope law and gradient optimizer"];

Print["s(r) = ", fmt[sRay[r]]];
Print["K_ij(r) = ", fmt[orientedSlope]];
Print["k_ij(r) = ", fmt[positiveSlope]];

expectZero[
  "M1 oriented slope law",
  orientedSlope + (ki + r*kj)/Sqrt[1 + r^2]
];
expectZero[
  "M1 positive slope magnitude",
  positiveSlope - (ki + r*kj)/Sqrt[1 + r^2]
];

slopeDerivative = cleanScalar[D[positiveSlope, r]];
expectZero[
  "M2 slope derivative law",
  slopeDerivative - (kj - ki*r)/(1 + r^2)^(3/2)
];

stationarySolutions = Solve[
  {D[positiveSlope, r] == 0, r >= 0},
  r,
  Reals
];
rGrad = cleanScalar[r /. First[stationarySolutions]];
stationaryRegion = Reduce[
  D[positiveSlope, r] == 0 && r >= 0,
  r,
  Reals
];
stationaryRegionOnDomain = FullSimplify[
  stationaryRegion,
  Assumptions -> $Assumptions
];

Print["Solve stationarity roots = ", fmt[stationarySolutions]];
Print["Reduce stationarity region = ", fmt[stationaryRegion]];
Print["domain-reduced stationarity region = ", fmt[stationaryRegionOnDomain]];

expectZero["M3 unique stationary root count", Length[stationarySolutions] - 1];
expectTrue[
  "M3 stationarity region equals r == kj/ki",
  Equivalent[stationaryRegionOnDomain, r == kj/ki]
];
expectZero["M3 gradient-optimal ratio", rGrad - kj/ki];
expectZero[
  "M3 gradient-optimal slope square",
  positiveSlope /. r -> rGrad // (#^2 - ki^2 - kj^2 &)
];

subbanner["M4-M5. Mixed curvature weights and cross-weight extremum"];

mixedCurvature = cleanScalar[sRay[r] . hessianPair . sRay[r]];
curvatureTarget = (hii + 2*r*hij + r^2*hjj)/(1 + r^2);
curvatureNumerator = Expand[Together[mixedCurvature*(1 + r^2)]];
linearHijNumerator = Normal[Series[curvatureNumerator, {hij, 0, 1}]];

weightI = cleanScalar[Coefficient[curvatureNumerator, hii]/(1 + r^2)];
weightX = cleanScalar[Coefficient[linearHijNumerator, hij]/(1 + r^2)];
weightJ = cleanScalar[Coefficient[curvatureNumerator, hjj]/(1 + r^2)];

Print["H1(r) = ", fmt[mixedCurvature]];
Print["coefficient weights {w_i, w_x, w_j} = ", fmt[{weightI, weightX, weightJ}]];

expectZero["M4 mixed curvature decomposition", mixedCurvature - curvatureTarget];
expectZero[
  "M4 diagonal neutrality",
  (mixedCurvature /. hij -> 0) - (hii + r^2*hjj)/(1 + r^2)
];
expectZero["M4 coefficient-recovered cross weight", weightX - 2*r/(1 + r^2)];

crossWeightDerivative = cleanScalar[D[weightX, r]];
expectZero[
  "M5 cross-weight derivative",
  crossWeightDerivative - 2*(1 - r^2)/(1 + r^2)^2
];
expectZero["M5 cross-weight derivative at r=1", crossWeightDerivative /. r -> 1];
expectZero["M5 cross-weight value at r=1", (weightX /. r -> 1) - 1];

subbanner["M6-M7. Canonical screen rays and coincidence condition"];

gradientDirection = cleanTensor[sRay[rGrad]];
unitGradientDirection = {ki, kj}/Sqrt[ki^2 + kj^2];
gradientSlope = cleanScalar[gradient . unitGradientDirection];
gradientCurvature = cleanScalar[unitGradientDirection . hessianPair . unitGradientDirection];

equalDirection = cleanTensor[sRay[1]];
equalSlope = cleanScalar[gradient . equalDirection];
equalCurvature = cleanScalar[equalDirection . hessianPair . equalDirection];

Print["s(r_grad) = ", fmt[gradientDirection]];
Print["s(1) = ", fmt[equalDirection]];

expectZero[
  "M6 gradient-optimal direction",
  gradientDirection - unitGradientDirection
];
expectZero[
  "M6 gradient-optimal slope",
  gradientSlope + Sqrt[ki^2 + kj^2]
];
expectZero[
  "M6 gradient-optimal curvature",
  gradientCurvature -
    (ki^2*hii + 2*ki*kj*hij + kj^2*hjj)/(ki^2 + kj^2)
];
expectZero["M6 equal-mix direction", equalDirection - {1, 1}/Sqrt[2]];
expectZero["M6 equal-mix slope", equalSlope + (ki + kj)/Sqrt[2]];
expectZero[
  "M6 equal-mix curvature",
  equalCurvature - (hii + 2*hij + hjj)/2
];

expectZero["M7 ratio coincidence difference", rGrad - 1 - (kj - ki)/ki];
expectZero[
  "M7 coincident directions when ki == kj",
  (gradientDirection - equalDirection) /. kj -> ki
];

subbanner["M8. Entrywise curvature envelopes"];

loHessian = {{hiiLo, hijLo}, {hijLo, hjjLo}};
hiHessian = {{hiiHi, hijHi}, {hijHi, hjjHi}};

kappaLoFromRayleigh = cleanScalar[sRay[r] . loHessian . sRay[r]];
kappaHiFromRayleigh = cleanScalar[sRay[r] . hiHessian . sRay[r]];
kappaLo = (hiiLo + 2*r*hijLo + r^2*hjjLo)/(1 + r^2);
kappaHi = (hiiHi + 2*r*hijHi + r^2*hjjHi)/(1 + r^2);

expectZero["M8 lower envelope Rayleigh form", kappaLoFromRayleigh - kappaLo];
expectZero["M8 upper envelope Rayleigh form", kappaHiFromRayleigh - kappaHi];
expectZero[
  "M8 lower envelope weighted form",
  kappaLo - (weightI*hiiLo + weightX*hijLo + weightJ*hjjLo)
];
expectZero[
  "M8 upper envelope weighted form",
  kappaHi - (weightI*hiiHi + weightX*hijHi + weightJ*hjjHi)
];

subbanner["M9. Certified bracket root closure"];

rootMap[h_, k_, c_] := 2*h/(k + Sqrt[k^2 - 2*c*h]);
tauLo = rootMap[H0, positiveSlope, kappaLo];
tauHi = rootMap[H0, positiveSlope, kappaHi];

closureResidual[tau_, kappa_] := H0 - positiveSlope*tau + (1/2)*kappa*tau^2;

Print["tau_lo(r) = ", fmt[tauLo]];
Print["tau_hi(r) = ", fmt[tauHi]];

expectZero[
  "M9 lower bracket closure quadratic",
  closureResidual[tauLo, kappaLo]
];
expectZero[
  "M9 upper bracket closure quadratic",
  closureResidual[tauHi, kappaHi]
];

banner["STAGE 208 MATHEMATICA AUDIT COMPLETED SUCCESSFULLY"];
