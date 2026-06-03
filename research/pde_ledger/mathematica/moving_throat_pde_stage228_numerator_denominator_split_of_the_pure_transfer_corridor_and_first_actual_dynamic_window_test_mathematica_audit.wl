ClearAll["Global`*"];
$HistoryLength = 0;
$MaxExtraPrecision = 1000;

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  detail -> ", fmt[detail]]];
  Exit[1];
);

stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

cleanScalar[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanScalar[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

expectClose[name_String, actual_, expected_, tol_] := Module[{diff},
  diff = N[Abs[actual - expected], 40];
  Print[name, " actual = ", fmt[N[actual, 18]], " expected = ", fmt[N[expected, 18]], " diff = ", fmt[diff]];
  If[TrueQ[diff < tol], pass[name], fail[name, diff]];
];

expectSmall[name_String, actual_, tol_] := Module[{mag},
  mag = N[Abs[actual], 40];
  Print[name, " magnitude = ", fmt[mag]];
  If[TrueQ[mag < tol], pass[name], fail[name, mag]];
];

expectVectorClose[name_String, actual_List, expected_List, tol_] := Module[{diffs, maxDiff},
  diffs = N[Abs[actual - expected], 40];
  maxDiff = Max[diffs];
  Print[name, " actual = ", fmt[N[actual, 14]], " expected = ", fmt[expected], " maxdiff = ", fmt[maxDiff]];
  If[TrueQ[maxDiff < tol], pass[name], fail[name, diffs]];
];

$Assumptions = (
  Element[{eps, y, omega, xK, xM, xLB, xV, xLU, xLW, xLR, xOU, xOW}, Reals]
  && kappa > 0 && K > 0 && M > 0 && lamB > 0 && varpi > 0
  && lamU > 0 && lamW > 0 && lamR > 0 && OmU > 0 && OmW > 0
  && a > 0 && cs > 0
);

Print["=== Stage 228 Mathematica audit: numerator/denominator split and dynamic-window test ==="];

mixedVars = {xLU, xLW, xLR, xOU, xOW};

capC = kappa lamB;
gU = lamU;
gW = kappa lamW;
rCouple = kappa lamR;

deltaGap = OmU^2 OmW^2 - rCouple^2;
s2 = OmU^2 + OmW^2;
hMix = gU^2 + gW^2;
qMix = gU^2 OmW^2 + 2 gU gW rCouple + gW^2 OmU^2;
pMix = OmU^2 gW + rCouple gU;

b0 = capC^2/varpi^2;
b2 = capC^2/varpi^4;
b4 = capC^2/varpi^6;
z0 = qMix/deltaGap;
z2 = (qMix s2 - hMix deltaGap)/deltaGap^2;
z4 = (qMix (s2^2 - deltaGap) - s2 hMix deltaGap)/deltaGap^3;
n0Expr = pMix^2/deltaGap^2;
p0Expr = n0Expr/(K - b0 - z0);

growthRules = {
  K -> K Exp[eps xK],
  M -> M Exp[eps xM],
  lamB -> lamB Exp[eps xLB],
  varpi -> varpi Exp[eps xV],
  lamU -> lamU Exp[eps xLU],
  lamW -> lamW Exp[eps xLW],
  lamR -> lamR Exp[eps xLR],
  OmU -> OmU Exp[eps xOU],
  OmW -> OmW Exp[eps xOW]
};

linCoeff[expr_] := FullSimplify[
  Coefficient[Normal[Series[expr /. growthRules, {eps, 0, 1}]], eps, 1],
  Assumptions -> $Assumptions
];

d01 = linCoeff[K - b0 - z0];
d21 = linCoeff[-(M + b2 + z2)];
d41 = linCoeff[-(b4 + z4)];
n01 = linCoeff[n0Expr];
p01 = linCoeff[pMix];
delta01 = linCoeff[deltaGap];

sampleRulesNoK = {
  kappa -> 2 Sqrt[2]/Pi,
  lamB -> 1/2,
  lamU -> 3/10,
  lamW -> 2/5,
  lamR -> 1/4,
  OmU -> 1,
  OmW -> 7/5,
  varpi -> 2,
  M -> 1,
  a -> 1,
  cs -> 1
};

b0Sample = FullSimplify[b0 /. sampleRulesNoK];
b2Sample = FullSimplify[b2 /. sampleRulesNoK];
b4Sample = FullSimplify[b4 /. sampleRulesNoK];
z0Sample = FullSimplify[z0 /. sampleRulesNoK];
z2Sample = FullSimplify[z2 /. sampleRulesNoK];
z4Sample = FullSimplify[z4 /. sampleRulesNoK];
kCompat = FullSimplify[
  b0Sample + z0Sample + 3 (1 + b2Sample + z2Sample)^2/(b4Sample + z4Sample)
];
sampleRulesFull = Join[sampleRulesNoK, {K -> kCompat}];
zeroNonmixed = {xK -> 0, xM -> 0, xLB -> 0, xV -> 0};

restrictMixed[expr_] := Expand[FullSimplify[(expr /. sampleRulesFull) /. zeroNonmixed]];
rowOf[expr_] := FullSimplify[Coefficient[Expand[expr], #] & /@ mixedVars];

d01Mixed = restrictMixed[d01];
d21Mixed = restrictMixed[d21];
d41Mixed = restrictMixed[d41];
n01Mixed = restrictMixed[n01];
p01Mixed = restrictMixed[p01];
delta01Mixed = restrictMixed[delta01];

pBase = FullSimplify[pMix /. sampleRulesFull];
deltaBase = FullSimplify[deltaGap /. sampleRulesFull];
n0Base = FullSimplify[n0Expr /. sampleRulesFull];

xi1 = Expand[FullSimplify[n01Mixed/n0Base]];
pi1 = Expand[FullSimplify[p01Mixed/pBase]];
delta1 = Expand[FullSimplify[delta01Mixed/deltaBase]];

Print["k_compat = ", fmt[kCompat]];
Print["Xi_1 = ", fmt[xi1]];
Print["pi_1 = ", fmt[pi1]];
Print["delta_1 = ", fmt[delta1]];

expectZero["M1 split identity", xi1 - 2 (pi1 - delta1)];

piCoeff = rowOf[pi1];
deltaCoeff = rowOf[delta1];
xiCoeff = rowOf[xi1];
expectedPi = {3/19, 16/19, 3/19, 32/19, 0};
expectedDelta = {0, 0, 50/(25 - 98 Pi^2), 196 Pi^2/(98 Pi^2 - 25), 196 Pi^2/(98 Pi^2 - 25)};

Print["pi_1 row = ", fmt[piCoeff]];
Print["delta_1 row = ", fmt[deltaCoeff]];
expectZero["M2 exact pi_1 row", Total[(piCoeff - expectedPi)^2]];
expectZero["M3 exact delta_1 row", Total[(deltaCoeff - expectedDelta)^2]];

transferMatrix = rowOf /@ {d01Mixed, d21Mixed, d41Mixed};
numMatrix = Join[transferMatrix, {piCoeff}];
denMatrix = Join[transferMatrix, {deltaCoeff}];
bothMatrix = Join[transferMatrix, {piCoeff, deltaCoeff}];

transferRank = MatrixRank[transferMatrix];
numRank = MatrixRank[numMatrix];
denRank = MatrixRank[denMatrix];
bothRank = MatrixRank[bothMatrix];
transferNull = NullSpace[transferMatrix];
numNull = NullSpace[numMatrix];
denNull = NullSpace[denMatrix];
bothNull = NullSpace[bothMatrix];

Print["transfer matrix = ", fmt[transferMatrix]];
expectZero["M4 pure-transfer rank", transferRank - 3];
expectZero["M4 pure-transfer nullity", Length[transferNull] - 2];
expectZero["M4 numerator-rigid rank", numRank - 4];
expectZero["M4 numerator-rigid nullity", Length[numNull] - 1];
expectZero["M4 denominator-rigid rank", denRank - 4];
expectZero["M4 denominator-rigid nullity", Length[denNull] - 1];
expectZero["M4 both-rigid rank", bothRank - 5];
expectZero["M4 both-rigid nullity", Length[bothNull]];

detReducedRaw = Factor[FullSimplify[Det[{piCoeff . Transpose[transferNull], deltaCoeff . Transpose[transferNull]}]]];
orientedTransferNull = transferNull;
If[TrueQ[N[detReducedRaw, 50] < 0], orientedTransferNull[[1]] = -orientedTransferNull[[1]]];
transferBasisColumns = Transpose[orientedTransferNull];
detReduced = Factor[FullSimplify[Det[{piCoeff . transferBasisColumns, deltaCoeff . transferBasisColumns}]]];
expectedDet = Factor[
  196 (200 + 147 Pi^2) (80000 + 343225 Pi^2 + 43218 Pi^4)/
    (475 (8670000 + 14894275 Pi^2 + 2117682 Pi^4))
];
Print["det[(pi_1, delta_1)|_pure transfer] = ", fmt[detReduced]];
expectZero["M5 reduced determinant", detReduced - expectedDet];
expectTrue["M5 reduced determinant nonzero", detReduced != 0];

unitPositiveXi[vec_] := Module[{unit, xiVal},
  unit = FullSimplify[vec/Sqrt[vec . vec]];
  xiVal = N[xiCoeff . unit, 50];
  If[TrueQ[xiVal < 0], -unit, unit]
];

vNum = unitPositiveXi[First[numNull]];
vDen = unitPositiveXi[First[denNull]];

expectedVNum = {-0.55551149, 0.31814576, -0.65766801, -0.04533730, -0.39447126};
expectedVDen = {-0.26583993, 0.18448137, 0.94454459, 0.04984499, -0.02543112};
expectVectorClose["M6 v_num", N[vNum, 30], expectedVNum, 10^-8];
expectVectorClose["M6 v_den", N[vDen, 30], expectedVDen, 10^-8];

deltaNum = N[deltaCoeff . vNum, 50];
xiNum = N[xiCoeff . vNum, 50];
piDen = N[piCoeff . vDen, 50];
xiDen = N[xiCoeff . vDen, 50];

expectZero["M6 pi_1(v_num)", piCoeff . vNum];
expectClose["M6 delta_1(v_num)", deltaNum, -0.86805617, 10^-8];
expectClose["M6 Xi_1(v_num)", xiNum, 1.73611235, 10^-8];
expectZero["M6 delta_1(v_den)", deltaCoeff . vDen];
expectClose["M6 pi_1(v_den)", piDen, 0.34646608, 10^-8];
expectClose["M6 Xi_1(v_den)", xiDen, 0.69293215, 10^-8];
expectClose["M6 Xi_1(v_num) + 2 delta_1(v_num)", xiNum + 2 deltaNum, 0, 10^-10];
expectClose["M6 Xi_1(v_den) - 2 pi_1(v_den)", xiDen - 2 piDen, 0, 10^-10];

fYExpr = Expand[
  (((K - M y) (varpi^2 - y) - capC^2) ((OmU^2 - y) (OmW^2 - y) - rCouple^2)
    - (varpi^2 - y) (gU^2 (OmW^2 - y) + 2 gU gW rCouple + gW^2 (OmU^2 - y)))
];
nYExpr = ((OmU^2 - y) gW + rCouple gU)^2/
  (((OmU^2 - y) (OmW^2 - y) - rCouple^2)^2);
rqYExpr = FullSimplify[27 cs^5/(a^5 y^(5/2) nYExpr)];

positiveRootsY[rules_] := Module[{poly, sol, vals},
  poly = N[FullSimplify[fYExpr /. rules], 60];
  sol = y /. NSolve[poly == 0, y, WorkingPrecision -> 60];
  vals = Select[sol, TrueQ[Abs[Im[N[#]]] < 10^-40 && Re[N[#]] > 0] &];
  Sort[Re /@ vals]
];

baseRootsY = positiveRootsY[sampleRulesFull];
wallY = Take[baseRootsY, -2];
wallOmega = Sqrt /@ wallY;
rqBase = N[(rqYExpr /. sampleRulesFull) /. y -> #, 50] & /@ wallY;
p0Base = N[p0Expr /. sampleRulesFull, 50];

Print["positive omega roots = ", fmt[N[Sqrt /@ baseRootsY, 18]]];
Print["wall omega roots = ", fmt[N[wallOmega, 18]]];
Print["wall R_Q values = ", fmt[N[rqBase, 18]]];
Print["P_0 = ", fmt[N[p0Base, 18]]];
expectClose["M7 omega_-", wallOmega[[1]], 1.99753568, 2 10^-8];
expectClose["M7 omega_+", wallOmega[[2]], 4.94905432, 2 10^-8];
expectClose["M7 R_Q,-", rqBase[[1]], 30.19990756, 2 10^-8];
expectClose["M7 R_Q,+", rqBase[[2]], 36.17118648, 2 10^-8];
expectClose["M7 P_0", p0Base, 0.00206979232, 2 10^-11];

rulesAlong[dir_List] := {
  kappa -> 2 Sqrt[2]/Pi,
  lamB -> 1/2,
  lamU -> (3/10) Exp[eps dir[[1]]],
  lamW -> (2/5) Exp[eps dir[[2]]],
  lamR -> (1/4) Exp[eps dir[[3]]],
  OmU -> Exp[eps dir[[4]]],
  OmW -> (7/5) Exp[eps dir[[5]]],
  varpi -> 2,
  M -> 1,
  K -> kCompat,
  a -> 1,
  cs -> 1
};

branchSlopes[dir_List] := Module[
  {rules, fBranch, rqBranch, p0Branch, fEps, fY, rqEps, rqY, dy, omegaSlopes, rqSlopes, p0Slope},
  rules = rulesAlong[dir];
  fBranch = fYExpr /. rules;
  rqBranch = rqYExpr /. rules;
  p0Branch = p0Expr /. rules;
  fEps = D[fBranch, eps];
  fY = D[fBranch, y];
  rqEps = D[Log[rqBranch], eps];
  rqY = D[Log[rqBranch], y];
  dy = N[-((fEps /. {eps -> 0, y -> #})/(fY /. {eps -> 0, y -> #})), 60] & /@ wallY;
  omegaSlopes = N[dy/(2 wallY), 50];
  rqSlopes = N[
    Table[
      (rqEps + rqY dy[[i]]) /. {eps -> 0, y -> wallY[[i]]},
      {i, Length[wallY]}
    ],
    50
  ];
  p0Slope = N[D[Log[p0Branch], eps] /. eps -> 0, 50];
  {omegaSlopes, rqSlopes, p0Slope}
];

numData = branchSlopes[vNum];
denData = branchSlopes[vDen];
numOmegaSlopes = numData[[1]];
numRqSlopes = numData[[2]];
numP0Slope = numData[[3]];
denOmegaSlopes = denData[[1]];
denRqSlopes = denData[[2]];
denP0Slope = denData[[3]];

Print["numerator-rigid omega slopes = ", fmt[N[numOmegaSlopes, 18]]];
Print["numerator-rigid R_Q slopes = ", fmt[N[numRqSlopes, 18]]];
Print["numerator-rigid dln P_0 = ", fmt[N[numP0Slope, 18]]];
Print["denominator-rigid omega slopes = ", fmt[N[denOmegaSlopes, 18]]];
Print["denominator-rigid R_Q slopes = ", fmt[N[denRqSlopes, 18]]];
Print["denominator-rigid dln P_0 = ", fmt[N[denP0Slope, 18]]];

expectClose["M8 numerator dln P_0", numP0Slope, 1.73611235, 10^-7];
expectClose["M8 numerator dln P_0 equals Xi_1", numP0Slope, xiNum, 5 10^-7];
expectClose["M8 numerator dln R_Q,+", numRqSlopes[[2]], -0.52346582, 10^-6];
expectClose["M8 numerator dln R_Q,-", numRqSlopes[[1]], 0.71358484, 10^-6];
expectSmall["M8 numerator omega_- slope", numOmegaSlopes[[1]], 5 10^-5];
expectSmall["M8 numerator omega_+ slope", numOmegaSlopes[[2]], 5 10^-5];

expectClose["M8 denominator dln P_0", denP0Slope, 0.69293215, 10^-7];
expectClose["M8 denominator dln P_0 equals Xi_1", denP0Slope, xiDen, 5 10^-7];
expectClose["M8 denominator dln R_Q,+", denRqSlopes[[2]], -0.35245541, 10^-6];
expectClose["M8 denominator dln R_Q,-", denRqSlopes[[1]], -0.23169484, 10^-6];
expectSmall["M8 denominator omega_- slope", denOmegaSlopes[[1]], 5 10^-5];
expectSmall["M8 denominator omega_+ slope", denOmegaSlopes[[2]], 5 10^-5];

finiteDynamicCeilings[rqVals_List, slopes_List, threshold_] := Module[
  {finite = {}, improving = False, tau},
  Do[
    If[TrueQ[N[slopes[[i]]] < 0],
      tau = N[Log[rqVals[[i]]/threshold]/(-slopes[[i]]), 40];
      finite = Append[finite, tau],
      improving = True
    ],
    {i, Length[slopes]}
  ];
  If[Length[finite] == 0,
    {Infinity, Infinity},
    {Min[finite], If[improving, Infinity, Max[finite]]}
  ]
];

threshold = 21.854566296358396;
budgetBoth = 0.367930328492646;
budgetNonempty = 0.737619063660757;

numDynamic = finiteDynamicCeilings[rqBase, numRqSlopes, threshold];
denDynamic = finiteDynamicCeilings[rqBase, denRqSlopes, threshold];
numStatic = {budgetBoth/xiNum, budgetNonempty/xiNum};
denStatic = {budgetBoth/xiDen, budgetNonempty/xiDen};

Print["numerator-rigid dynamic ceilings = ", fmt[N[numDynamic, 18]]];
Print["denominator-rigid dynamic ceilings = ", fmt[N[denDynamic, 18]]];
Print["numerator-rigid static ceilings = ", fmt[N[numStatic, 18]]];
Print["denominator-rigid static ceilings = ", fmt[N[denStatic, 18]]];

expectClose["M9 numerator both dynamic", numDynamic[[1]], 0.96253269, 10^-6];
expectTrue["M9 numerator nonempty dynamic is Infinity", numDynamic[[2]] === Infinity];
expectClose["M9 denominator both dynamic", denDynamic[[1]], 1.39592653, 10^-6];
expectClose["M9 denominator nonempty dynamic", denDynamic[[2]], 1.42955095, 10^-6];
expectClose["M9 numerator both static", numStatic[[1]], 0.21192772, 5 10^-8];
expectClose["M9 numerator nonempty static", numStatic[[2]], 0.42486828, 5 10^-8];
expectClose["M9 denominator both static", denStatic[[1]], 0.53097598, 5 10^-8];
expectClose["M9 denominator nonempty static", denStatic[[2]], 1.06448959, 5 10^-8];
expectTrue["M9 numerator both dynamic exceeds transported static", numDynamic[[1]] > numStatic[[1]]];
expectTrue["M9 denominator both dynamic exceeds transported static", denDynamic[[1]] > denStatic[[1]]];
expectTrue["M9 denominator nonempty dynamic exceeds transported static", denDynamic[[2]] > denStatic[[2]]];

Print[""];
Print["Stage 228 Mathematica audit completed successfully."];
