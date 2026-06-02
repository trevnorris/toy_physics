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

cleanWith[expr_, assumptions_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Cancel[Together[res]], Assumptions -> assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> assumptions]
];

zeroQ[expr_] := And @@ (TrueQ[# === 0] & /@ Flatten[{expr}]);

expectZeroUnder[name_String, expr_, assumptions_] := Module[{res},
  res = cleanWith[expr, assumptions];
  Print[name, " = ", fmt[res]];
  If[zeroQ[res], pass[name], fail[name, res]];
];

expectZero[name_String, expr_] := expectZeroUnder[name, expr, $Assumptions];

expectTrueUnder[name_String, cond_, assumptions_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := expectTrueUnder[name, cond, $Assumptions];

totalDegree[poly_, vars_List] := Module[{mons},
  mons = DeleteCases[MonomialList[Expand[poly], vars], 0];
  If[Length[mons] == 0,
    -Infinity,
    Max[Table[Total[Table[Exponent[mon, var], {var, vars}]], {mon, mons}]]
  ]
];

banner["STAGE 214 -- FOUR-COORDINATE SIMPLEX OPTIMIZER MATHEMATICA AUDIT"];

ClearAll[
  r, s, t, y, ki, kj, kk, kl, H0, u, k, ud, ux,
  aa, bb, cc, dd, ee, ff, gg, hh, ii, jj
];

ratioVars = {r, s, t};
liftVars = {r, s, t, y};
coefVars = {aa, bb, cc, dd, ee, ff, gg, hh, ii, jj};
kVec = {ki, kj, kk, kl};
ratioVec = {1, r, s, t};
den = ratioVec.ratioVec;
linearK = kVec.ratioVec;
quadEnvelope = (
  aa + bb r + cc s + dd t + ee r^2 + ff r s + gg r t
    + hh s^2 + ii s t + jj t^2
);
phi = (linearK + Sqrt[quadEnvelope])/Sqrt[den];
tau = 2 H0/phi;

$Assumptions = (
  Element[Join[ratioVars, {y, ki, kj, kk, kl, H0, u, k, ud, ux}, coefVars], Reals]
  && ki > 0 && kj > 0 && kk > 0 && kl > 0
  && H0 > 0 && k > 0
  && den > 0 && quadEnvelope > 0
);

subbanner["M1. Exact derivative laws and stationary numerators"];

metricNumerators = Table[
  FullSimplify[den D[linearK, var] - D[den, var] linearK/2, Assumptions -> $Assumptions],
  {var, ratioVars}
];
envelopeNumerators = Table[
  FullSimplify[den D[quadEnvelope, var] - D[den, var] quadEnvelope, Assumptions -> $Assumptions],
  {var, ratioVars}
];
stationaryNumerators = 2 metricNumerators Sqrt[quadEnvelope] + envelopeNumerators;
derivativeDen = 2 den^(3/2) Sqrt[quadEnvelope];
directNumerators = Table[
  FullSimplify[derivativeDen D[phi, var], Assumptions -> $Assumptions],
  {var, ratioVars}
];

Print["M1 direct derivative numerator list = ", fmt[directNumerators]];
expectZero["M1 direct numerator r minus N_r", directNumerators[[1]] - stationaryNumerators[[1]]];
expectZero["M1 direct numerator s minus N_s", directNumerators[[2]] - stationaryNumerators[[2]]];
expectZero["M1 direct numerator t minus N_t", directNumerators[[3]] - stationaryNumerators[[3]]];
expectZero["M1 derivative law r residual", D[phi, r] - stationaryNumerators[[1]]/derivativeDen];
expectZero["M1 derivative law s residual", D[phi, s] - stationaryNumerators[[2]]/derivativeDen];
expectZero["M1 derivative law t residual", D[phi, t] - stationaryNumerators[[3]]/derivativeDen];

subbanner["M2. Lifted polynomial degree ledger"];

liftedPolys = Expand[2 # y + #2] & @@@ Transpose[{metricNumerators, envelopeNumerators}];
deltaLift = y^2 - quadEnvelope;
liftedDegrees = totalDegree[#, liftVars] & /@ Append[liftedPolys, deltaLift];

Print["M2 lifted degrees {F_r,F_s,F_t,F_Delta} = ", fmt[liftedDegrees]];
expectZero["M2 degree F_r minus 3", liftedDegrees[[1]] - 3];
expectZero["M2 degree F_s minus 3", liftedDegrees[[2]] - 3];
expectZero["M2 degree F_t minus 3", liftedDegrees[[3]] - 3];
expectZero["M2 degree F_Delta minus 2", liftedDegrees[[4]] - 2];

subbanner["M3. Preferred lifted Bezout bound"];

liftedBezout = Times @@ liftedDegrees;
Print["M3 lifted degree product = ", liftedBezout];
expectZero["M3 lifted Bezout product minus 3*3*3*2", liftedBezout - 3*3*3*2];

subbanner["M4. Projected square-root-free eliminant degrees"];

crossEliminants = {
  -Cancel[Resultant[liftedPolys[[1]], liftedPolys[[2]], y]/2],
  -Cancel[Resultant[liftedPolys[[1]], liftedPolys[[3]], y]/2],
  -Cancel[Resultant[liftedPolys[[2]], liftedPolys[[3]], y]/2]
};
crossDefinitions = {
  metricNumerators[[2]] envelopeNumerators[[1]]
    - metricNumerators[[1]] envelopeNumerators[[2]],
  metricNumerators[[3]] envelopeNumerators[[1]]
    - metricNumerators[[1]] envelopeNumerators[[3]],
  metricNumerators[[3]] envelopeNumerators[[2]]
    - metricNumerators[[2]] envelopeNumerators[[3]]
};
squareEliminants = Resultant[#, deltaLift, y] & /@ liftedPolys;
squareDefinitions = envelopeNumerators^2 - 4 metricNumerators^2 quadEnvelope;

expectZero["M4 C_rs resultant minus definition", crossEliminants[[1]] - crossDefinitions[[1]]];
expectZero["M4 C_rt resultant minus definition", crossEliminants[[2]] - crossDefinitions[[2]]];
expectZero["M4 C_st resultant minus definition", crossEliminants[[3]] - crossDefinitions[[3]]];
expectZero["M4 S_r resultant minus definition", squareEliminants[[1]] - squareDefinitions[[1]]];
expectZero["M4 S_s resultant minus definition", squareEliminants[[2]] - squareDefinitions[[2]]];
expectZero["M4 S_t resultant minus definition", squareEliminants[[3]] - squareDefinitions[[3]]];

crossDegrees = totalDegree[#, ratioVars] & /@ crossEliminants;
squareDegrees = totalDegree[#, ratioVars] & /@ squareEliminants;
Print["M4 cross degrees {C_rs,C_rt,C_st} = ", fmt[crossDegrees]];
Print["M4 square degrees {S_r,S_s,S_t} = ", fmt[squareDegrees]];
expectZero["M4 C_rs degree minus 5", crossDegrees[[1]] - 5];
expectZero["M4 C_rt degree minus 5", crossDegrees[[2]] - 5];
expectZero["M4 C_st degree minus 5", crossDegrees[[3]] - 5];
expectZero["M4 S_r degree minus 6", squareDegrees[[1]] - 6];
expectZero["M4 S_s degree minus 6", squareDegrees[[2]] - 6];
expectZero["M4 S_t degree minus 6", squareDegrees[[3]] - 6];

subbanner["M5. Projected one-chart Bezout product"];

projectedProduct = crossDegrees[[1]] crossDegrees[[2]] squareDegrees[[1]];
Print["M5 projected one-chart degree product = ", projectedProduct];
expectZero["M5 projected product minus 5*5*6", projectedProduct - 5*5*6];

subbanner["M6. Diagonal-isotropic reduction and gradient ray"];

isoRules = {
  aa -> ki^2 - 2 H0 u,
  bb -> 2 ki kj,
  cc -> 2 ki kk,
  dd -> 2 ki kl,
  ee -> kj^2 - 2 H0 u,
  ff -> 2 kj kk,
  gg -> 2 kj kl,
  hh -> kk^2 - 2 H0 u,
  ii -> 2 kk kl,
  jj -> kl^2 - 2 H0 u
};
isoDelta = quadEnvelope /. isoRules;
isoDeltaExpected = linearK^2 - 2 H0 u den;
kRST = linearK/Sqrt[den];
tauIso = tau /. isoRules;
tauExpected = 2 H0/(kRST + Sqrt[kRST^2 - 2 H0 u]);
isoAssumptions = (
  $Assumptions && isoDeltaExpected > 0 && kRST^2 - 2 H0 u > 0
);

expectZero["M6 diagonal-isotropic Delta residual", isoDelta - isoDeltaExpected];
expectZeroUnder["M6 diagonal-isotropic tau residual", tauIso - tauExpected, isoAssumptions];

kNorm = Sqrt[kVec.kVec];
gradRay = kVec/kNorm;
Print["M6 gradient ray = ", fmt[gradRay]];
expectZero["M6 gradient ray unit norm", gradRay.gradRay - 1];
expectZero["M6 gradient ray slope residual", gradRay.kVec - kNorm];

subbanner["M7. Full-symmetry equal-mix stationary point"];

symRules = {
  ki -> k,
  kj -> k,
  kk -> k,
  kl -> k,
  aa -> k^2 - 2 H0 ud,
  bb -> 2 k^2 - 4 H0 ux,
  cc -> 2 k^2 - 4 H0 ux,
  dd -> 2 k^2 - 4 H0 ux,
  ee -> k^2 - 2 H0 ud,
  ff -> 2 k^2 - 4 H0 ux,
  gg -> 2 k^2 - 4 H0 ux,
  hh -> k^2 - 2 H0 ud,
  ii -> 2 k^2 - 4 H0 ux,
  jj -> k^2 - 2 H0 ud
};
equalMixNumerators = stationaryNumerators /. symRules /. {r -> 1, s -> 1, t -> 1};
equalMixRay = {1, 1, 1, 1}/2;

expectZero["M7 symmetric N_r(1,1,1)", equalMixNumerators[[1]]];
expectZero["M7 symmetric N_s(1,1,1)", equalMixNumerators[[2]]];
expectZero["M7 symmetric N_t(1,1,1)", equalMixNumerators[[3]]];
expectZero["M7 equal-mix unit norm", equalMixRay.equalMixRay - 1];

Print[""];
Print["All Stage 214 Mathematica audit checks passed."];
Exit[0];
