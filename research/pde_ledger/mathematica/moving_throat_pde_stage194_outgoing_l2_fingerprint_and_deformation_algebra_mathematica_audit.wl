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
stripCE[expr_] := expr /. ConditionalExpression[e_, _] :> e;

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[stripCE[expr]]], Assumptions -> $Assumptions];
  res = stripCE[res];
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

positiveRoot[name_String, roots_List] := Module[{clean, selected},
  clean = DeleteDuplicates[FullSimplify[stripCE /@ roots, Assumptions -> $Assumptions]];
  selected = Select[clean, TrueQ[FullSimplify[# > 0, Assumptions -> $Assumptions]] &];
  If[Length[selected] =!= 1, fail[name, clean]];
  First[selected]
];

banner["STAGE 194 — OUTGOING l=2 FINGERPRINT AND DEFORMATION ALGEBRA"];

Clear[x, freq, radius, soundSpeed, chiQ, poleScale, sigmaCan];
Clear[scaleS, betaStretch, sigma0, sigma2, sigma4, sigma5];
Clear[bigG, lightSpeed, baseK0];

$Assumptions =
  Element[
    {x, freq, radius, soundSpeed, chiQ, poleScale, sigmaCan,
     scaleS, betaStretch, sigma0, sigma2, sigma4, sigma5,
     bigG, lightSpeed, baseK0},
    Reals
  ] &&
  radius > 0 && soundSpeed > 0 && poleScale > 0 && sigmaCan > 0 &&
  bigG > 0 && lightSpeed > 0 &&
  scaleS != 0 && betaStretch != 0 && 3*scaleS - sigma0 != 0;

(* I. Exact outgoing spherical l=2 DtN fingerprint. *)
subbanner["I. Exact outgoing spherical l=2 DtN fingerprint"];

hankel2 = FunctionExpand[SphericalHankelH1[2, x]];
lambdaExact = FullSimplify[x*D[hankel2, x]/hankel2, Assumptions -> $Assumptions];
yNormalized = FullSimplify[-3/lambdaExact, Assumptions -> $Assumptions];

lambdaSeries = Expand[Normal[Series[lambdaExact, {x, 0, 7}]]];
ySeries = Expand[Normal[Series[yNormalized, {x, 0, 7}]]];

lambdaTarget = (
  -3
  + x^2/3
  + x^4/9
  + I*x^5/9
  - 2*x^6/27
  - I*x^7/27
);
yTarget = (
  1
  + x^2/9
  + 4*x^4/81
  + I*x^5/27
  - 11*x^6/729
  - I*x^7/243
);

Print["h_2^(1)(x) = ", fmt[hankel2]];
Print["Lambda_2^out(x) series = ", fmt[lambdaSeries]];
Print["Yhat_2^out(x) series = ", fmt[ySeries]];

expectZero["Lambda_out series - SymPy target", lambdaSeries - lambdaTarget];
expectZero["Yhat_out series - SymPy target", ySeries - yTarget];
expectZero["static DtN slot + 3", Coefficient[lambdaSeries, x, 0] + 3];

(* II. Exact matching to the retarded grouped-P2 one-pole module. *)
subbanner["II. Retarded one-pole matching"];

retardedKernel = (
  3/4
  + 1/(4*(1 - freq^2/poleScale^2 - I*chiQ*sigmaCan*freq^5))
);
retardedSeriesRaw = Expand[Normal[Series[retardedKernel, {freq, 0, 5}]]];
outgoingOmegaSeries = Expand[
  Normal[Series[yNormalized /. x -> radius*freq/soundSpeed, {freq, 0, 5}]]
];

poleCandidates = poleScale /. Solve[
  Coefficient[retardedSeriesRaw, freq, 2] ==
    Coefficient[outgoingOmegaSeries, freq, 2],
  poleScale,
  Reals
];
poleFromEven = positiveRoot["unique positive pole from omega^2 coefficient", poleCandidates];

sigmaCandidates = sigmaCan /. Solve[
  (Coefficient[retardedSeriesRaw, freq, 5]/I /. chiQ -> 1) ==
    Coefficient[outgoingOmegaSeries, freq, 5]/I,
  sigmaCan,
  Reals
];
sigmaFromOdd = positiveRoot["unique positive sigma from outgoing omega^5 coefficient", sigmaCandidates];

retardedSeries = FullSimplify[
  retardedSeriesRaw /. {poleScale -> poleFromEven, sigmaCan -> sigmaFromOdd},
  Assumptions -> $Assumptions
];

retardedTarget = (
  1
  + radius^2*freq^2/(9*soundSpeed^2)
  + 4*radius^4*freq^4/(81*soundSpeed^4)
  + I*chiQ*radius^5*freq^5/(27*soundSpeed^5)
);

Print["positive pole from even matching = ", fmt[poleFromEven]];
Print["positive sigma from odd normalization = ", fmt[sigmaFromOdd]];
Print["retarded one-pole series = ", fmt[retardedSeries]];

expectZero["Omega_Q - 3 c_s/(2 a)", poleFromEven - 3*soundSpeed/(2*radius)];
expectZero["sigma_Q^can - 4 a^5/(27 c_s^5)", sigmaFromOdd - 4*radius^5/(27*soundSpeed^5)];
expectZero["sigma_Q^can - 9/(8 Omega_Q^5)", sigmaFromOdd - 9/(8*poleFromEven^5)];
expectZero["Yret series - SymPy low-frequency target", retardedSeries - retardedTarget];
expectZero[
  "retarded branch with chi_Q=1 - outgoing DtN branch",
  (retardedSeries /. chiQ -> 1) - outgoingOmegaSeries
];

oddCoeffRet = Coefficient[retardedSeries, freq, 5];
oddCoeffOut = Coefficient[outgoingOmegaSeries, freq, 5];
chiCondition = FullSimplify[
  Reduce[ComplexExpand[oddCoeffRet/I == oddCoeffOut/I], chiQ, Reals],
  Assumptions -> $Assumptions
];
Print["chi_Q condition from odd coefficient = ", fmt[chiCondition]];
expectZero[
  "odd-coefficient matching fixes chi_Q - 1",
  (oddCoeffRet - oddCoeffOut)*27*soundSpeed^5/(I*radius^5) - (chiQ - 1)
];

(* III. Exact isotropic DtN deformation algebra. *)
subbanner["III. Isotropic DtN deformation algebra"];

lambdaWindow5 = Expand[Normal[Series[lambdaExact, {x, 0, 5}]]];
deformedOperator = Expand[
  scaleS*(lambdaWindow5 /. x -> betaStretch*x)
  + sigma0
  + sigma2*x^2
  + sigma4*x^4
  + I*sigma5*x^5
];

slot0 = FullSimplify[Coefficient[deformedOperator, x, 0], Assumptions -> $Assumptions];
slot2 = FullSimplify[Coefficient[deformedOperator, x, 2], Assumptions -> $Assumptions];
slot4 = FullSimplify[Coefficient[deformedOperator, x, 4], Assumptions -> $Assumptions];
slot5 = FullSimplify[Coefficient[deformedOperator, x, 5]/I, Assumptions -> $Assumptions];

normalizedDeformed = Expand[Normal[Series[slot0/deformedOperator, {x, 0, 5}]]];
compilerFormula = (
  1
  - slot2*x^2/slot0
  + (slot2^2/slot0^2 - slot4/slot0)*x^4
  - I*slot5*x^5/slot0
);

Print["deformed DtN operator through x^5 = ", fmt[deformedOperator]];
Print["L0 = ", fmt[slot0]];
Print["L2 = ", fmt[slot2]];
Print["L4 = ", fmt[slot4]];
Print["L5 = ", fmt[slot5]];
Print["Yhat_2^def(x) series = ", fmt[normalizedDeformed]];

expectZero["L0 formula", slot0 - (-3*scaleS + sigma0)];
expectZero["L2 formula", slot2 - (scaleS*betaStretch^2/3 + sigma2)];
expectZero["L4 formula", slot4 - (scaleS*betaStretch^4/9 + sigma4)];
expectZero["L5 formula", slot5 - (scaleS*betaStretch^5/9 + sigma5)];
expectZero["Ydef series - exact compiler", normalizedDeformed - compilerFormula];

evenSolutions = Solve[
  {
    Coefficient[normalizedDeformed, x, 2] == 1/9,
    Coefficient[normalizedDeformed, x, 4] == 4/81
  },
  {sigma2, sigma4},
  Reals
];
If[Length[evenSolutions] =!= 1, fail["unique canonical-even deformation solution", evenSolutions]];
evenRule = First[evenSolutions];

sigma2Target = -(3*scaleS*betaStretch^2 - 3*scaleS + sigma0)/9;
sigma4Target = -(3*scaleS*betaStretch^4 - 3*scaleS + sigma0)/27;

Print["Sigma_2 from canonical-even matching = ", fmt[FullSimplify[sigma2 /. evenRule, Assumptions -> $Assumptions]]];
Print["Sigma_4 from canonical-even matching = ", fmt[FullSimplify[sigma4 /. evenRule, Assumptions -> $Assumptions]]];

expectZero["Sigma_2 formula", (sigma2 /. evenRule) - sigma2Target];
expectZero["Sigma_4 formula", (sigma4 /. evenRule) - sigma4Target];

chiFromDeformation = FullSimplify[
  Coefficient[normalizedDeformed, x, 5]/(I/27),
  Assumptions -> $Assumptions
];
chiEven = FullSimplify[chiFromDeformation /. evenRule, Assumptions -> $Assumptions];
chiTarget = 3*(scaleS*betaStretch^5 + 9*sigma5)/(3*scaleS - sigma0);
chiMinusOneTarget = (3*scaleS*(betaStretch^5 - 1) + sigma0 + 27*sigma5)/(3*scaleS - sigma0);

Print["chi_Q from canonical-even deformed branch = ", fmt[chiEven]];
expectZero["chi_Q deformation law", chiEven - chiTarget];
expectZero["chi_Q - 1 deformation law", (chiEven - 1) - chiMinusOneTarget];

(* IV. Carry-forward corollary: canonical invariant tuple. *)
subbanner["IV. Canonical invariant tuple"];

k0Candidates = baseK0 /. Solve[
  9*baseK0/(32*poleFromEven^5) == 2*bigG/(5*lightSpeed^5),
  baseK0,
  Reals
];
k0Canonical = positiveRoot["unique positive Kbar_0 from Gamma_5 relation", k0Candidates];
k2Canonical = FullSimplify[k0Canonical/(4*poleFromEven^2), Assumptions -> $Assumptions];
k4Canonical = FullSimplify[k0Canonical/(4*poleFromEven^4), Assumptions -> $Assumptions];
gamma5Canonical = FullSimplify[9*k0Canonical/(32*poleFromEven^5), Assumptions -> $Assumptions];

Print["Kbar_0 = ", fmt[k0Canonical]];
Print["Kbar_2 = ", fmt[k2Canonical]];
Print["Kbar_4 = ", fmt[k4Canonical]];
Print["Gammabar_5 = ", fmt[gamma5Canonical]];

expectZero["Kbar_0 - 54 G c_s^5/(5 a^5 c^5)", k0Canonical - 54*bigG*soundSpeed^5/(5*radius^5*lightSpeed^5)];
expectZero["Kbar_2 - 6 G c_s^3/(5 a^3 c^5)", k2Canonical - 6*bigG*soundSpeed^3/(5*radius^3*lightSpeed^5)];
expectZero["Kbar_4 - 8 G c_s/(15 a c^5)", k4Canonical - 8*bigG*soundSpeed/(15*radius*lightSpeed^5)];
expectZero["Gammabar_5 - 2 G/(5 c^5)", gamma5Canonical - 2*bigG/(5*lightSpeed^5)];

banner["STAGE 194 MATHEMATICA LEDGER"];
Print["1. The native SphericalHankelH1 branch gives the exact outgoing l=2 DtN series."];
Print["2. Coefficient matching of the retarded one-pole module fixes chi_Q = 1."];
Print["3. The canonical-even isotropic DtN deformation leaves only beta, Sigma_0, and Sigma_5"];
Print["   in the outgoing-normalization scalar."];
Print["4. The canonical point-particle branch returns the carried invariant tuple."];

Exit[0];
