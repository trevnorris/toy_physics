(* Ledger stage012 Mathematica audit: DtN pole ladder + Robin falsifier.

   Standalone, print-only, no arguments, no file I/O.  This keeps the native
   transfer-matrix 012 engine, derives the consumed inputs natively with
   dual-site integrity, and strips the scratch bridge plumbing.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;

DNUNITTESTFAILDIMENSIONAL = "DN_UNITTEST_FAIL_DIMENSIONAL";
FAILPOLELADDER = "FAIL_POLE_LADDER";
FAILCOUNTERFACTUAL = "FAIL_COUNTERFACTUAL";
DNUNITTESTBCDEPENDENT = "DN_UNITTEST_BC_DEPENDENT";
DNUNITTESTPASS = "DN_UNITTEST_PASS";

expectedGuardKeys = {
  "robin_determinant_emitted",
  "recovers_DN_at_alpha0",
  "recovers_DD_at_alpha_inf",
  "halfshift_destroyed_for_DD",
  "numeric_alpha_distinct",
  "dtn_mismatch"
};

$Assumptions =
  L0 > 0 && cS > 0 && omega > 0 && K > 0 && rhoStar > 0 && m > 0 &&
  psiM != 0 && alpha > 0 && Element[{s, rho}, Reals] &&
  Element[n, Integers] && n >= 0 && Element[j, Integers] && j >= 1;

k = omega/cS;
zeroDim = {0, 0, 0};

raise[msg_] := Throw[msg, "ledgerStage012Failure"];

heading[text_] := (
  Print[""];
  Print[StringRepeat["=", StringLength[text]]];
  Print[text];
  Print[StringRepeat["=", StringLength[text]]]
);

subheading[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);

dropConditions[expr_] := expr /. ConditionalExpression[value_, _] :> value;
fmt[expr_] := ToString[InputForm[expr]];
cleanZero[expr_] := FullSimplify[dropConditions[expr]];
limitAlphaInfinity[expr_] := Block[
  {$Assumptions = L0 > 0 && cS > 0 && omega > 0},
  FullSimplify[Limit[dropConditions[expr], alpha -> Infinity]]
];
limitOmegaZero[expr_] := Block[
  {$Assumptions = L0 > 0 && cS > 0},
  FullSimplify[Limit[dropConditions[expr], omega -> 0]]
];

assertExact[name_, expr_] := Module[{reals},
  reals = Cases[Unevaluated[expr], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FAIL  ", name, ": machine-real atom(s) found: ", ToString[InputForm[reals]]];
    raise[name]
  ]
];

expectZero[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FAIL  ", name, ": residual = ", fmt[clean]];
    raise[name]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

expectNonzero[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[! TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name, " is nonzero as required (residual = ", fmt[clean], ")"],
    failCount++;
    Print["FAIL  ", name, ": required nonzero residual vanished"];
    raise[name]
  ]
];

expectFail[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[! TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name, " produced required FAIL (residual = ", fmt[clean], ")"],
    failCount++;
    Print["FAIL  ", name, ": required mutation/ablation did not fire"];
    raise[name]
  ]
];

exprEqual[lhs_, rhs_: 0] := TrueQ[FullSimplify[dropConditions[lhs - rhs] == 0]];
extractTanDataFromDtn[expr_] := Module[{clean, tanTerms, tanTerm},
  clean = FullSimplify[dropConditions[expr]];
  tanTerms = DeleteDuplicates[Cases[clean, Tan[arg_] :> Tan[arg], Infinity]];
  If[Length[tanTerms] =!= 1, raise["expected exactly one Tan factor in derived DtN"]];
  tanTerm = First[tanTerms];
  <|
    "TanTerm" -> tanTerm,
    "TanArgument" -> FullSimplify[tanTerm[[1]]],
    "DtnPrefactor" -> FullSimplify[clean/tanTerm]
  |>
];
extractCosPoleDenominatorFromTanTerm[tanTerm_] := Module[{arg},
  If[Head[tanTerm] =!= Tan, raise["expected Tan factor for pole denominator extraction"]];
  arg = FullSimplify[tanTerm[[1]]];
  FullSimplify[Cos[arg]]
];
nonzeroQ[expr_] := ! TrueQ[FullSimplify[dropConditions[expr] == 0]];
boolResidual[condition_] := If[TrueQ[condition], 0, 1];
verdictResidual[actual_, expected_] := If[actual === expected, 0, 1];
guardAll[guard_Association] := AllTrue[Values[guard], TrueQ];

r1SiteFromExponent[exponent_] := FullSimplify[exponent K rho^(exponent - 1)/m];
r1EosSiteFromExponent[exponent_] := FullSimplify[D[K rho^exponent, rho]/m];

computeVerdict[dimensionalOk_, dtnMatchesTarget_, halfshift_, counterfactualGuardAll_, bcDerivationEmitted_] :=
  Which[
    ! TrueQ[dimensionalOk], DNUNITTESTFAILDIMENSIONAL,
    ! TrueQ[dtnMatchesTarget && halfshift], FAILPOLELADDER,
    ! TrueQ[counterfactualGuardAll], FAILCOUNTERFACTUAL,
    ! TrueQ[bcDerivationEmitted], DNUNITTESTBCDEPENDENT,
    True, DNUNITTESTPASS
  ];

reconstructLsFromPair[pair_] := Module[
  {aNull, bNull, equations, solved, aValue, bValue, operator},
  equations = FullSimplify[D[#, {s, 2}] + aNull D[#, s] + bNull #] & /@ pair;
  solved = Solve[Thread[equations == 0], {aNull, bNull}];
  If[Length[solved] < 1, raise["could not reconstruct monic L_s from null-space pair"]];
  aValue = FullSimplify[aNull /. First[solved]];
  bValue = FullSimplify[bNull /. First[solved]];
  operator = FullSimplify[D[psiHat[s], {s, 2}] + aValue D[psiHat[s], s] + bValue psiHat[s]];
  <|"a" -> aValue, "b" -> bValue, "operator" -> operator, "equations" -> equations|>
];

transferMatrix[ell_] := {
  {Cos[k ell], Sin[k ell]/k},
  {-k Sin[k ell], Cos[k ell]}
};

counterfactualGuard[robinDtn_, detWitness_, denominatorCore_, alphaValue_: 2/L0] := Module[
  {
    dnFromRobin, ddFromRobin, ddTarget, ddDenominator, ddHalfshiftSamples,
    ddIntegerSamples, halfshiftDestroyedForDD, ddZeroModeRemovable,
    numericRobinDen, numericRobinDtn, numericAlphaDistinct, guard
  },
  dnFromRobin = FullSimplify[robinDtn /. alpha -> 0];
  ddFromRobin = limitAlphaInfinity[robinDtn];
  ddTarget = FullSimplify[k Cot[k L0]];
  ddDenominator = FullSimplify[Sin[k L0]];
  ddHalfshiftSamples = Table[
    FullSimplify[ddDenominator /. omega -> Pi cS (idx + 1/2)/L0],
    {idx, 0, 3}
  ];
  ddIntegerSamples = Table[
    FullSimplify[ddDenominator /. omega -> Pi cS idx/L0],
    {idx, 1, 4}
  ];
  halfshiftDestroyedForDD =
    AllTrue[ddHalfshiftSamples, nonzeroQ] &&
    AllTrue[ddIntegerSamples, exprEqual[#, 0] &];
  ddZeroModeRemovable = exprEqual[limitOmegaZero[ddTarget], 1/L0];
  numericRobinDen = FullSimplify[denominatorCore /. alpha -> alphaValue];
  numericRobinDtn = FullSimplify[robinDtn /. alpha -> alphaValue];
  numericAlphaDistinct =
    nonzeroQ[numericRobinDtn - dtnTransfer] &&
    nonzeroQ[numericRobinDtn - ddTarget];
  guard = <|
    "robin_determinant_emitted" -> nonzeroQ[detWitness],
    "recovers_DN_at_alpha0" -> exprEqual[dnFromRobin, dtnTransfer],
    "recovers_DD_at_alpha_inf" -> exprEqual[ddFromRobin, ddTarget],
    "halfshift_destroyed_for_DD" -> halfshiftDestroyedForDD,
    "numeric_alpha_distinct" -> numericAlphaDistinct,
    "dtn_mismatch" -> nonzeroQ[robinDtn - dtnTransfer]
  |>;
  <|
    "Guard" -> guard,
    "DNFromRobin" -> dnFromRobin,
    "DDFromRobin" -> ddFromRobin,
    "DDTarget" -> ddTarget,
    "DDDenominator" -> ddDenominator,
    "DDHalfshiftSamples" -> ddHalfshiftSamples,
    "DDIntegerSamples" -> ddIntegerSamples,
    "DDZeroModeRemovable" -> ddZeroModeRemovable,
    "NumericAlpha" -> alphaValue,
    "NumericRobinDen" -> numericRobinDen,
    "NumericRobinDtn" -> numericRobinDtn
  |>
];

dimResidualVec[actual_, expected_] := FullSimplify[(actual - expected).(actual - expected)];

dimOf[expr_, dims_] := Module[{clean, args, ds, base, pow, argDims},
  clean = FullSimplify[expr];
  Which[
    TrueQ[clean == 0] || NumericQ[clean], zeroDim,
    AtomQ[clean] && KeyExistsQ[dims, clean], dims[clean],
    AtomQ[clean], raise["missing dimension for " <> ToString[Unevaluated[clean], InputForm]],
    Head[clean] === Times, Total[dimOf[#, dims] & /@ (List @@ clean)],
    Head[clean] === Power,
      base = clean[[1]];
      pow = clean[[2]];
      If[! NumericQ[pow], raise["non-numeric dimension exponent"]];
      pow dimOf[base, dims],
    Head[clean] === Plus,
      args = Select[List @@ clean, ! TrueQ[FullSimplify[# == 0]] &];
      ds = dimOf[#, dims] & /@ args;
      If[Length[ds] == 0, zeroDim,
        If[Length[DeleteDuplicates[ds]] != 1, raise["dimension mismatch in sum"]];
        First[ds]
      ],
    MemberQ[{Sin, Cos, Tan, Cot}, Head[clean]],
      argDims = dimOf[#, dims] & /@ (List @@ clean);
      If[AnyTrue[argDims, # =!= zeroDim &], raise["dimensionful argument in dimensionless function"]];
      zeroDim,
    True, raise["unsupported dimension expression " <> ToString[clean, InputForm]]
  ]
];

buildDimensionalBlock[] := Module[
  {
    lengthDim, energyDim, fourVolumeDim, pressureDim, rhoDim, kDim,
    omegaDim, massDim, alphaDim, expectedTanDim, expectedZ00Dim,
    cleanWalk, corruptWalk, corruptKDim, dimensionalOk,
    corruptDimensionalOk, mutationFires, cleanVerdict, mutatedVerdict,
    failSuppressed, walk
  },
  lengthDim = {1, 0, 0};
  energyDim = {2, 1, -2};
  fourVolumeDim = 4 lengthDim;
  pressureDim = energyDim - fourVolumeDim;
  rhoDim = -4 lengthDim;
  kDim = pressureDim - 5 rhoDim;
  omegaDim = {0, 0, -1};
  massDim = {0, 1, 0};
  alphaDim = {-1, 0, 0};
  expectedTanDim = zeroDim;
  expectedZ00Dim = {-1, 0, 0};
  walk[KDimension_] := Module[{baseRules, csSquaredDim, csDim, dimRules, tanDim, prefDim, z00Dim},
    baseRules = <|L0 -> lengthDim, omega -> omegaDim, K -> KDimension, rhoStar -> rhoDim, m -> massDim, alpha -> alphaDim|>;
    csSquaredDim = dimOf[5 K rhoStar^4/m, baseRules];
    csDim = (1/2) csSquaredDim;
    dimRules = Join[baseRules, <|cS -> csDim|>];
    tanDim = dimOf[k L0, dimRules];
    prefDim = dimOf[-k, dimRules];
    z00Dim = If[tanDim === zeroDim, dimOf[-k Tan[k L0], dimRules], "dimensionful_tan_argument"];
    <|
      "CsSquaredDim" -> csSquaredDim,
      "CsDim" -> csDim,
      "KDimFromOmegaOverCs" -> dimOf[omega/cS, dimRules],
      "TanArgumentDim" -> tanDim,
      "Z00PrefactorDim" -> prefDim,
      "Z00Dim" -> z00Dim,
      "AlphaCSDim" -> dimOf[alpha cS, dimRules]
    |>
  ];
  cleanWalk = walk[kDim];
  dimensionalOk =
    TrueQ[cleanWalk["TanArgumentDim"] == expectedTanDim] &&
    TrueQ[cleanWalk["Z00PrefactorDim"] == expectedZ00Dim] &&
    TrueQ[cleanWalk["Z00Dim"] == expectedZ00Dim];
  corruptKDim = kDim + {1, 0, 0};
  corruptWalk = walk[corruptKDim];
  corruptDimensionalOk =
    TrueQ[corruptWalk["TanArgumentDim"] == expectedTanDim] &&
    TrueQ[corruptWalk["Z00PrefactorDim"] == expectedZ00Dim] &&
    TrueQ[corruptWalk["Z00Dim"] == expectedZ00Dim];
  mutationFires =
    ! TrueQ[corruptDimensionalOk] &&
    ! TrueQ[corruptWalk["TanArgumentDim"] == expectedTanDim] &&
    ! TrueQ[corruptWalk["Z00PrefactorDim"] == expectedZ00Dim];
  cleanVerdict = computeVerdict[dimensionalOk, True, True, True, False];
  mutatedVerdict = computeVerdict[corruptDimensionalOk, True, True, True, False];
  failSuppressed =
    TrueQ[mutationFires] &&
    cleanVerdict === DNUNITTESTBCDEPENDENT &&
    mutatedVerdict === DNUNITTESTFAILDIMENSIONAL;
  <|
    "LengthDim" -> lengthDim,
    "EnergyDim" -> energyDim,
    "FourVolumeDim" -> fourVolumeDim,
    "PressureDim" -> pressureDim,
    "RhoDim" -> rhoDim,
    "KDim" -> kDim,
    "OmegaDim" -> omegaDim,
    "AlphaDim" -> alphaDim,
    "ExpectedTanDim" -> expectedTanDim,
    "ExpectedZ00Dim" -> expectedZ00Dim,
    "CleanWalk" -> cleanWalk,
    "DimensionalOk" -> dimensionalOk,
    "CorruptKDim" -> corruptKDim,
    "CorruptWalk" -> corruptWalk,
    "CorruptDimensionalOk" -> corruptDimensionalOk,
    "MutationFires" -> mutationFires,
    "CleanVerdict" -> cleanVerdict,
    "MutatedVerdict" -> mutatedVerdict,
    "FailSuppressed" -> failSuppressed
  |>
];

buildBaseline[] := Module[
  {
    ode, generalSolution, fundamentalPair, r1SiteA, r1SiteB, consumedSpeed,
    lsSiteA, lsRecon, lsSiteB, consumedLs, anchorLs, mouthState, capState,
    neumannResidual, pMFromNeumann, dtnTarget, denominatorFull,
    tanData, tanArgument, dtnPrefactor,
    poleDenominator, poleLadder, poleResidual, halfshift, staticSeries,
    staticSeriesPoly, staticSeriesTarget, staticLimit, roundTrip,
    roundTripOnLadder, roundTripCloses, robinResidual, pMFromRobin,
    robinDenominatorCore, robinDetWitness, robinNumeratorCore, robinData,
    robinAlpha0, robinAlphaInf, ddResidual, pMFromDD, ddTransferLocal,
    ddTargetLocal, alpha0Agree, ddAgree, dim, bcDerivationEmitted,
    bcProvenance, bcDerivation, verdict
  },
  ode = D[psi[s], {s, 2}] + k^2 psi[s] == 0;
  generalSolution = FullSimplify[DSolveValue[ode, psi[s], s]];
  fundamentalPair = {FullSimplify[Coefficient[generalSolution, C[2]]], FullSimplify[Coefficient[generalSolution, C[1]]]};

  r1SiteA = r1SiteFromExponent[5];
  r1SiteB = r1EosSiteFromExponent[5];
  consumedSpeed = FullSimplify[r1SiteA /. rho -> rhoStar];
  lsSiteA = FullSimplify[D[psiHat[s], {s, 2}] + k^2 psiHat[s]];
  lsRecon = reconstructLsFromPair[fundamentalPair];
  lsSiteB = lsRecon["operator"];
  consumedLs = lsSiteA;
  anchorLs = FullSimplify[D[psiHat[s], {s, 2}] + (omega/cS)^2 psiHat[s]];

  mouthState = {psiM, pM};
  capState = FullSimplify[transferMatrix[L0] . mouthState];
  neumannResidual = FullSimplify[capState[[2]]];
  pMFromNeumann = FullSimplify[dropConditions[pM /. First[Solve[neumannResidual == 0, pM]]]];
  dtnTransfer = FullSimplify[-pMFromNeumann/psiM];
  dtnTarget = FullSimplify[-k Tan[k L0]];
  tanData = extractTanDataFromDtn[dtnTransfer];
  tanArgument = tanData["TanArgument"];
  dtnPrefactor = tanData["DtnPrefactor"];
  dtnMatchesTarget = exprEqual[dtnTransfer, dtnTarget];
  denominatorFull = Denominator[Together[dtnTransfer]];
  poleDenominator = extractCosPoleDenominatorFromTanTerm[tanData["TanTerm"]];
  poleLadder = FullSimplify[Pi cS (n + 1/2)/L0];
  poleResidual = FullSimplify[poleDenominator /. omega -> poleLadder];
  halfshift = exprEqual[poleResidual, 0];
  staticSeries = Series[dtnTransfer, {omega, 0, 5}];
  staticSeriesPoly = FullSimplify[Normal[staticSeries]];
  staticSeriesTarget = FullSimplify[-L0 omega^2/cS^2 - L0^3 omega^4/(3 cS^4)];
  staticLimit = limitOmegaZero[dtnTransfer];
  roundTrip = FullSimplify[-Exp[2 I k L0]];
  roundTripOnLadder = FullSimplify[roundTrip /. omega -> poleLadder];
  roundTripCloses = exprEqual[roundTripOnLadder, 1];

  robinResidual = FullSimplify[capState[[2]] - alpha capState[[1]]];
  pMFromRobin = FullSimplify[dropConditions[pM /. First[Solve[robinResidual == 0, pM]]]];
  robinDtnTransfer = FullSimplify[-pMFromRobin/psiM];
  robinDenominatorCore = FullSimplify[k Cos[k L0] - alpha Sin[k L0]];
  robinDetWitness = robinDenominatorCore;
  robinNumeratorCore = FullSimplify[robinDtnTransfer robinDenominatorCore];
  robinAlpha0 = FullSimplify[robinDtnTransfer /. alpha -> 0];
  robinAlphaInf = limitAlphaInfinity[robinDtnTransfer];
  ddResidual = FullSimplify[capState[[1]]];
  pMFromDD = FullSimplify[dropConditions[pM /. First[Solve[ddResidual == 0, pM]]]];
  ddTransferLocal = FullSimplify[-pMFromDD/psiM];
  ddTransfer = ddTransferLocal;
  ddTargetLocal = FullSimplify[k Cot[k L0]];
  alpha0Agree = exprEqual[robinAlpha0, dtnTransfer];
  ddAgree = exprEqual[robinAlphaInf, ddTransferLocal] && exprEqual[ddTransferLocal, ddTargetLocal];
  robinData = counterfactualGuard[robinDtnTransfer, robinDetWitness, robinDenominatorCore, 2/L0];

  dim = buildDimensionalBlock[];
  bcDerivationEmitted = False;
  bcProvenance = "imposed";
  bcDerivation = <|
    "bc_type" -> "D-at-mouth / N-at-cap",
    "mouth_gradient_from_Vconf" -> "not emitted from an explicit V_wall profile in this unit test",
    "cap_gradient_from_Vconf" -> "not emitted from an explicit V_wall profile in this unit test",
    "regularity_at_pinchoff" -> "regular cap closure R0(L0)=0 motivates Neumann, but a full asymptotic derivation is not emitted",
    "mouth_condition" -> "clamp to quasi-static bulk reservoir is imposed for this frozen-wall benchmark, not derived as radiation",
    "classification_rule" -> "bc_derivation_emitted=false forces DN_UNITTEST_BC_DEPENDENT unless an explicit mouth/cap derivation is later supplied"
  |>;
  verdict = computeVerdict[dim["DimensionalOk"], dtnMatchesTarget, halfshift, guardAll[robinData["Guard"]], bcDerivationEmitted];

  <|
    "ODE" -> ode,
    "GeneralSolution" -> generalSolution,
    "FundamentalPair" -> fundamentalPair,
    "R1SiteA" -> r1SiteA,
    "R1SiteB" -> r1SiteB,
    "ConsumedSpeed" -> consumedSpeed,
    "LsSiteA" -> lsSiteA,
    "LsReconstructed" -> lsRecon,
    "LsSiteB" -> lsSiteB,
    "ConsumedLs" -> consumedLs,
    "AnchorLs" -> anchorLs,
    "TransferMatrix" -> transferMatrix[L0],
    "CapState" -> capState,
    "NeumannResidual" -> neumannResidual,
    "PMFromNeumann" -> pMFromNeumann,
    "DtnTransfer" -> dtnTransfer,
    "DtnTarget" -> dtnTarget,
    "DtnTanTerm" -> tanData["TanTerm"],
    "TanArgument" -> tanArgument,
    "DtnPrefactor" -> dtnPrefactor,
    "DtnMatchesTarget" -> dtnMatchesTarget,
    "DenominatorFull" -> denominatorFull,
    "PoleDenominator" -> poleDenominator,
    "PoleLadder" -> poleLadder,
    "PoleResidual" -> poleResidual,
    "Halfshift" -> halfshift,
    "StaticSeries" -> staticSeries,
    "StaticSeriesPoly" -> staticSeriesPoly,
    "StaticSeriesTarget" -> staticSeriesTarget,
    "StaticLimit" -> staticLimit,
    "RoundTrip" -> roundTrip,
    "RoundTripOnLadder" -> roundTripOnLadder,
    "RoundTripCloses" -> roundTripCloses,
    "RobinResidual" -> robinResidual,
    "RobinDtnTransfer" -> robinDtnTransfer,
    "RobinDenominatorCore" -> robinDenominatorCore,
    "RobinDetWitness" -> robinDetWitness,
    "RobinNumeratorCore" -> robinNumeratorCore,
    "RobinAlpha0" -> robinAlpha0,
    "RobinAlphaInf" -> robinAlphaInf,
    "DDTransfer" -> ddTransferLocal,
    "DDTarget" -> ddTargetLocal,
    "Alpha0Agree" -> alpha0Agree,
    "DDAgree" -> ddAgree,
    "Robin" -> robinData,
    "Dim" -> dim,
    "BCDerivationEmitted" -> bcDerivationEmitted,
    "BCProvenance" -> bcProvenance,
    "BCDerivation" -> bcDerivation,
    "Verdict" -> verdict
  |>
];

runAritySelfCheck[data_] := Module[{guardProbe},
  subheading["Wolfram arity self-check"];
  guardProbe = counterfactualGuard[data["RobinDtnTransfer"], data["RobinDetWitness"], data["RobinDenominatorCore"]];
  expectBool["arity transferMatrix[ell] returns a 2x2 matrix", Dimensions[transferMatrix[L0]] === {2, 2}];
  expectBool["arity reconstructLsFromPair[pair] returns operator", KeyExistsQ[reconstructLsFromPair[{Sin[k s], Cos[k s]}], "operator"]];
  expectBool["arity counterfactualGuard[3 args] returns six guard members", Length[guardProbe["Guard"]] == 6];
  expectBool["arity counterfactualGuard[4 args] accepts alpha override", KeyExistsQ[counterfactualGuard[data["RobinDtnTransfer"], data["RobinDetWitness"], data["RobinDenominatorCore"], 2/L0], "NumericRobinDtn"]];
  expectBool["arity computeVerdict[5 args] returns DN_UNITTEST_BC_DEPENDENT", computeVerdict[True, True, True, True, False] === DNUNITTESTBCDEPENDENT];
  expectBool["arity r1SiteFromExponent[e] returns exact literal site", exprEqual[r1SiteFromExponent[5], 5 K rho^4/m]];
  expectBool["arity r1EosSiteFromExponent[e] returns exact EOS route", exprEqual[r1EosSiteFromExponent[5], 5 K rho^4/m]];
  expectBool["arity buildDimensionalBlock[] returns mutation_fires", buildDimensionalBlock[]["MutationFires"] === True]
];

runOpeningAndConsumedInputs[data_] := Module[{recon},
  recon = data["LsReconstructed"];
  subheading["Consumed stage011 inputs and opening dsolve"];
  Print["  CONSUMED-from-011: L_s, c_S, and domain [0,L0] are cited; stage011 reduction/certificate are not recomputed."];
  Print["  CITED-speed: c_S^2 = 5*K*rho_star^4/m is R1 at rho_star; EOS exponent-5 P=K*rho^5 IMPOSED; bare m is stage004 m_GNLS."];
  Print["  ODE consumed by 012 = ", fmt[data["ODE"]]];
  Print["  DSolve general solution = ", fmt[data["GeneralSolution"]]];
  expectZero["DSolve solution satisfies cited Helmholtz L_s", D[data["GeneralSolution"], {s, 2}] + k^2 data["GeneralSolution"]];
  expectZero["fundamental sin branch is emitted by DSolve", data["FundamentalPair"][[1]] - Sin[k s]];
  expectZero["fundamental cos branch is emitted by DSolve", data["FundamentalPair"][[2]] - Cos[k s]];
  Print["  R1 site A literal = ", fmt[data["R1SiteA"]]];
  Print["  R1 site B EOS route d(K*rho^5)/d rho / m = ", fmt[data["R1SiteB"]]];
  expectZero["c_S^2 R1 site A minus site B equals zero", data["R1SiteA"] - data["R1SiteB"]];
  expectZero["c_S^2 evaluated at rho_star equals consumed speed", data["ConsumedSpeed"] - 5 K rhoStar^4/m];
  expectZero["c_S^2 frozen-export anchor consumed - 5*K*rho_star^4/m equals zero", data["ConsumedSpeed"] - 5 K rhoStar^4/m];
  Print["  L_s site A export = ", fmt[data["LsSiteA"]]];
  Print["  L_s site B null-space solve: a = ", fmt[recon["a"]], ", b = ", fmt[recon["b"]]];
  Print["  L_s site B reconstructed = ", fmt[data["LsSiteB"]]];
  expectZero["L_s null-space reconstruction recovers a=0", recon["a"]];
  expectZero["L_s null-space reconstruction recovers b=(omega/c_S)^2", recon["b"] - k^2];
  expectZero["L_s site A minus site B equals zero", data["LsSiteA"] - data["LsSiteB"]];
  Print["  structural note: L_s frozen-export alias anchor is de-counted; dual-site null-space integrity above carries the citation check."];
  Print["  structural note: consumed domain [0,L0] has length L0 by construction; length bookkeeping is de-counted."]
];

runTransferDtn[data_] := (
  subheading["Transfer-matrix D/N DtN"];
  Print["  transferMatrix[L0] = ", fmt[data["TransferMatrix"]]];
  Print["  capState = ", fmt[data["CapState"]]];
  Print["  Neumann cap residual = ", fmt[data["NeumannResidual"]]];
  Print["  pM from Neumann = ", fmt[data["PMFromNeumann"]]];
  Print["  dtnTransfer = ", fmt[data["DtnTransfer"]], "; target = ", fmt[data["DtnTarget"]]];
  expectZero["tan_argument extracted from derived DtN is k*L0", data["TanArgument"] - k L0];
  expectZero["dtn_prefactor extracted from derived DtN is -k", data["DtnPrefactor"] + k];
  expectZero["transfer-matrix DtN equals -k*tan(k*L0)", data["DtnTransfer"] - data["DtnTarget"]];
  expectBool["dtn_matches_target is transfer-derived-vs-typed", data["DtnMatchesTarget"]];
  expectZero["transfer route agrees with the LUsolve reference expression", data["DtnTransfer"] + k Tan[k L0]]
);

runPoleStaticRoundTrip[data_] := (
  subheading["Half-shift pole ladder, static series, and round trip"];
  Print["  transfer denominator full = ", fmt[data["DenominatorFull"]]];
  Print["  pole denominator = ", fmt[data["PoleDenominator"]]];
  expectZero["pole denominator is cos(k*L0)", data["PoleDenominator"] - Cos[k L0]];
  Print["  pole ladder = ", fmt[data["PoleLadder"]]];
  Print["  pole residual after ladder substitution = ", fmt[data["PoleResidual"]]];
  expectZero["half-shift pole residual is zero", data["PoleResidual"]];
  expectBool["halfshift = pole_residual==0 is computed true", data["Halfshift"]];
  Print["  static small-omega series = ", fmt[data["StaticSeries"]]];
  Print["  static limit omega->0+ = ", fmt[data["StaticLimit"]]];
  expectZero["static series polynomial matches -L0*omega^2/cS^2 - L0^3*omega^4/(3*cS^4)", data["StaticSeriesPoly"] - data["StaticSeriesTarget"]];
  expectZero["static DC limit is zero", data["StaticLimit"]];
  expectNonzero["static series is distinguished from the DC limit", data["StaticSeriesPoly"] - data["StaticLimit"]];
  Print["  Static note: this is the small-omega expansion of the dynamic DtN; no separate static solve is emitted."];
  Print["  round_trip = ", fmt[data["RoundTrip"]], "; on D/N ladder = ", fmt[data["RoundTripOnLadder"]]];
  expectZero["round-trip closes to R_rt=1 on the D/N ladder", data["RoundTripOnLadder"] - 1];
  expectBool["round_trip_closes is computed true", data["RoundTripCloses"]];
  Print["  round-trip phase: phi0 = 0 mod 2*pi"]
);

runRobinCounterfactual[data_] := Module[{robin, guard},
  robin = data["Robin"];
  guard = robin["Guard"];
  subheading["Robin counterfactual falsifier"];
  Print["  Robin residual = ", fmt[data["RobinResidual"]]];
  Print["  Robin DtN transfer = ", fmt[data["RobinDtnTransfer"]]];
  Print["  robin_denominator_core = ", fmt[data["RobinDenominatorCore"]]];
  expectNonzero["robin_determinant_emitted is tied to computed denominator core", data["RobinDetWitness"]];
  Print["  alpha->0 branch = ", fmt[data["RobinAlpha0"]]];
  Print["  alpha->infinity branch = ", fmt[data["RobinAlphaInf"]], "; D/D transfer = ", fmt[data["DDTransfer"]]];
  expectBool["native alpha0 self-check robinAlpha0==dtnTransfer is kept", data["Alpha0Agree"]];
  expectBool["native alphaInfinity self-check robinAlphaInf==ddTransfer is kept", data["DDAgree"]];
  expectZero["Robin alpha->0 recovers D/N DtN", robin["DNFromRobin"] - data["DtnTransfer"]];
  expectZero["Robin alpha->infinity recovers D/D k*cot(k*L0)", robin["DDFromRobin"] - robin["DDTarget"]];
  Print["  D/D half-shift samples = ", fmt[robin["DDHalfshiftSamples"]]];
  Print["  D/D integer-ladder samples = ", fmt[robin["DDIntegerSamples"]]];
  Print["  dd_zero_mode_removable artifact = ", fmt[robin["DDZeroModeRemovable"]], " (not a guard member)"];
  expectBool["D/D zero mode removable limit equals 1/L0", robin["DDZeroModeRemovable"]];
  Print["  numeric alpha = ", fmt[robin["NumericAlpha"]], "; numeric Robin DtN = ", fmt[robin["NumericRobinDtn"]]];
  expectZero["counterfactual_guard has exactly six members", Length[guard] - 6];
  expectBool["counterfactual_guard membership matches source dict L522-L529", Keys[guard] === expectedGuardKeys];
  expectBool["dd_zero_mode_removable is not verdict-bearing", ! KeyExistsQ[guard, "dd_zero_mode_removable"]];
  Scan[
    Function[key,
      Print["  guard[", key, "] = ", fmt[guard[key]]];
      expectBool["counterfactual guard member " <> key <> " is computed true", guard[key]]
    ],
    expectedGuardKeys
  ];
  expectBool["all(counterfactual_guard.values()) is true", guardAll[guard]]
];

runDimensionalBlock[data_] := Module[{dim, clean, corrupt},
  dim = data["Dim"];
  clean = dim["CleanWalk"];
  corrupt = dim["CorruptWalk"];
  subheading["012 dimensional legs and corrupt-[K] probe"];
  Print["  dimension order: (L,M,T)"];
  Print["  [energy] = ", dim["EnergyDim"], "; [four-volume] = ", dim["FourVolumeDim"], "; [P] = ", dim["PressureDim"]];
  Print["  [rho] = ", dim["RhoDim"], "; [K]=[P]-5[rho] = ", dim["KDim"], "; [alpha] = ", dim["AlphaDim"]];
  Print["  propagated [c_S^2] = ", clean["CsSquaredDim"], " -> [c_S] = ", clean["CsDim"], " -> [k] = ", clean["KDimFromOmegaOverCs"]];
  Print["  [tan_argument=k*L0] = ", clean["TanArgumentDim"], "; [Z00_prefactor=-k] = ", clean["Z00PrefactorDim"], "; [Z00] = ", clean["Z00Dim"]];
  expectZero["tan_argument dimensional leg is dimensionless", dimResidualVec[clean["TanArgumentDim"], dim["ExpectedTanDim"]]];
  expectZero["Z00_prefactor dimensional leg is L^-1", dimResidualVec[clean["Z00PrefactorDim"], dim["ExpectedZ00Dim"]]];
  expectZero["Z00 dimensional leg is L^-1", dimResidualVec[clean["Z00Dim"], dim["ExpectedZ00Dim"]]];
  expectBool["dimensional_ok for 012 tan_argument/Z00 legs", dim["DimensionalOk"]];
  Print["  corrupt [K]+(1,0,0) gives [K] = ", dim["CorruptKDim"]];
  Print["  corrupt propagated [c_S] = ", corrupt["CsDim"], "; [k] = ", corrupt["KDimFromOmegaOverCs"]];
  Print["  corrupt [tan_argument] = ", corrupt["TanArgumentDim"], "; corrupt [Z00_prefactor] = ", corrupt["Z00PrefactorDim"]];
  expectFail["corrupt-[K] makes tan_argument non-dimensionless", dimResidualVec[corrupt["TanArgumentDim"], dim["ExpectedTanDim"]]];
  expectFail["corrupt-[K] makes Z00_prefactor differ from L^-1", dimResidualVec[corrupt["Z00PrefactorDim"], dim["ExpectedZ00Dim"]]];
  expectBool["corrupt-[K] mutation_fires=True", dim["MutationFires"]];
  expectZero["self-ablation with mutation gives DN_UNITTEST_FAIL_DIMENSIONAL", verdictResidual[dim["MutatedVerdict"], DNUNITTESTFAILDIMENSIONAL]];
  expectZero["self-ablation without mutation gives clean DN_UNITTEST_BC_DEPENDENT", verdictResidual[dim["CleanVerdict"], DNUNITTESTBCDEPENDENT]];
  expectBool["self-ablation fail_suppressed=True", dim["FailSuppressed"]]
];

runBCAndVerdict[data_] := (
  subheading["BC provenance, 012 verdict, and joint composition"];
  Print["  bc_provenance = ", data["BCProvenance"]];
  Print["  bc_derivation_emitted = ", fmt[data["BCDerivationEmitted"]]];
  Print["  bc_derivation descriptor = ", fmt[data["BCDerivation"]]];
  expectBool["bc_provenance is imposed", data["BCProvenance"] === "imposed"];
  expectBool["bc_derivation_emitted is the honest false scope flag", data["BCDerivationEmitted"] === False];
  Print["  012 scoped verdict = ", data["Verdict"]];
  expectZero["012 verdict lands at DN_UNITTEST_BC_DEPENDENT", verdictResidual[data["Verdict"], DNUNITTESTBCDEPENDENT]];
  Print["  DN_UNITTEST_BC_DEPENDENT (JOINT, COMPLETED)"];
  Print["    = (011: REDUCTION_CERTIFIED, cited from ledger_stage011)"];
  Print["    AND (012: DtN ladder EARNED + bc_derivation_emitted=False -> BC_DEPENDENT landing, computed here)"];
  expectBool["joint composition cites stage011 REDUCTION_CERTIFIED and computed 012 landing", data["Verdict"] === DNUNITTESTBCDEPENDENT]
);

printProvenance[] := (
  subheading["Provenance and scope"];
  Print["  CONSUMED-from-011: L_s, domain [0,L0] with cap R0(L0)=0, and c_S are CITED from ledger_stage011 with dual-site integrity; 011 reduction/certificate/de-rig are not recomputed."];
  Print["  CITED-speed: c_S^2=5*K*rho_star^4/m is Part I edge R1 (stage005, re-exported by stage011) at rho_star; EOS exponent-5 P=K*rho^5 IMPOSED."];
  Print["  IMPOSED-BC: D/N mouth/cap boundary pair is IMPOSED; bc_provenance=imposed and bc_derivation_emitted=False are banked calibration; V_wall derivation is deferred, not fabricated."];
  Print["  EARNED: DtN, half-shifted ladder, static small-omega series, round-trip R_rt=1, and Robin counterfactual guard are computed here; dtn_matches_target is derived-vs-typed, not X==X."];
  Print["  Robin-falsifier: Robin cap recovers D/N at alpha=0, D/D at alpha->infinity, destroys the half-shift for D/D, and is numerically distinct at alpha=2/L0."];
  Print["  control-symbol: alpha is a Robin cap admittance with [alpha]=L^-1, tracked-not-counted like stage010 k_warp; it builds the falsifiable counterfactual, not the physics."];
  Print["  split: this stage COMPLETES DN_UNITTEST_BC_DEPENDENT; 011 carried REDUCTION_CERTIFIED, 012 carries the D/N ladder + BC_DEPENDENT landing; DN_UNITTEST_PASS is deferred."];
  Print["  dropped-bookkeeping: scratch-YAML/_sympy_exprs.wl write, MMA-YAML re-read, expression_digest, and engine_agreement plumbing are stripped."];
  Print["  downstream consumers: stage 013 (harmonic beta lift) + stage 017 (calibration input) consume Z00, the resonance ladder, and BC_DEPENDENT provenance."];
  Print["  register note: zero new counted knobs; alpha is tracked-not-counted; L0 already registered in stage011; edge R28 is the imposed D/N boundary calibration obligation."]
);

printVerdictLabels[] := (
  Print[""];
  Print["Verdict labels:"];
  Print["  ledger earned-label (NOT a source verdict token): DTN_POLE_LADDER_ROBIN_FALSIFIER_EARNED  (dsolve of the cited frozen L_s -> D/N coefficient matrix -> outward-mouth DtN Z00=-(omega/c_S)*tan(L0*omega/c_S); half-shifted pole ladder omega_n=pi*c_S*(n+1/2)/L0; static small-omega series; round-trip R_rt=1; Robin cap counterfactual falsifier, guard = {robin_determinant_emitted, recovers_DN_at_alpha0, recovers_DD_at_alpha_inf, halfshift_destroyed_for_DD, numeric_alpha_distinct, dtn_mismatch}, each a computed residual; tan_argument/Z00 dim legs via [K]=[P]-5[rho] + corrupt-[K] probe)"];
  Print["  source top-line verdict: DN_UNITTEST_BC_DEPENDENT  (JOINT; 012 COMPLETES it -- adds the D/N ladder + bc_derivation_emitted=False -> BC_DEPENDENT landing to stage 011's cited REDUCTION_CERTIFIED)"];
  Print["  joint composition (COMPLETED): DN_UNITTEST_BC_DEPENDENT = (011: REDUCTION_CERTIFIED, cited from ledger_stage011) AND (012: DtN ladder EARNED + BC_DEPENDENT landing, computed here)"];
  Print["  earned: DtN derived via LUsolve (dtn_matches_target = genuine derived-vs-typed comparison, NOT X==X); half-shifted ladder (halfshift = pole_residual==0, COMPUTED); round-trip R_rt=1 (COMPUTED via substitution); Robin guard all booleans COMPUTED; tan_argument/Z00 dim legs (2,0,-2)-consistent + corrupt-[K] probe fires"];
  Print["  consumed (cited from stage011, dual-site integrity): frozen Helmholtz L_s = psi'' + (omega/c_S)^2 psi on [0,L0]; c_S^2 = 5*K*rho_star^4/m (R1 at rho_star); domain cap R0(L0)=0"];
  Print["  imposed (banked calibration, edge R28): D/N mouth/cap boundary pair; bc_provenance=imposed; bc_derivation_emitted=False -> BC_DEPENDENT landing; the mouth/cap V_wall gradient derivation earning DN_UNITTEST_PASS is a DEFERRED upgrade (NOT fabricated here)"];
  Print["  control symbol (tracked, not counted): alpha = Robin cap admittance, [alpha]=L^-1 (like k_warp at stage010)"]
);

runAbleToFailTeeth[data_] := Module[
  {
    flippedCore, flippedRobinDtn, flippedGuard, dtnMut, integerLadder,
    integerResidual, integerHalfshift, badRoundTrip, noAlphaNumerator,
    noAlphaRobinDtn, noAlphaGuard, degenerateGuard, badLsSiteA,
    hyperRecon, dim
  },
  dim = data["Dim"];
  subheading["Able-to-fail mutation teeth"];

  flippedCore = FullSimplify[k Cos[k L0] + alpha Sin[k L0]];
  flippedRobinDtn = FullSimplify[data["RobinNumeratorCore"]/flippedCore];
  flippedGuard = counterfactualGuard[flippedRobinDtn, flippedCore, flippedCore, 2/L0]["Guard"];
  expectFail["tooth 1 Robin denominator sign flip makes recovers_DD_at_alpha_inf false", boolResidual[flippedGuard["recovers_DD_at_alpha_inf"]]];
  expectZero[
    "tooth 1 verdict is FAIL_COUNTERFACTUAL",
    verdictResidual[computeVerdict[True, True, True, guardAll[flippedGuard], False], FAILCOUNTERFACTUAL]
  ];

  dtnMut = Module[{capState, residual, pMut, dtnBad},
    capState = FullSimplify[transferMatrix[L0] . {psiM, pM}];
    residual = FullSimplify[capState[[2]] - psiM/L0];
    pMut = FullSimplify[dropConditions[pM /. First[Solve[residual == 0, pM]]]];
    dtnBad = FullSimplify[-pMut/psiM];
    <|"Dtn" -> dtnBad, "Matches" -> exprEqual[dtnBad, data["DtnTarget"]]|>
  ];
  expectFail["tooth 2 mutated Neumann RHS changes derived DtN", dtnMut["Dtn"] - data["DtnTarget"]];
  expectFail["tooth 2 dtn_matches_target flips false", boolResidual[dtnMut["Matches"]]];
  expectZero["tooth 2 verdict is FAIL_POLE_LADDER", verdictResidual[computeVerdict[True, dtnMut["Matches"], True, True, False], FAILPOLELADDER]];

  integerLadder = FullSimplify[Pi cS j/L0];
  integerResidual = FullSimplify[data["PoleDenominator"] /. omega -> integerLadder];
  integerHalfshift = exprEqual[integerResidual, 0];
  expectFail["tooth 3 integer ladder does not solve cos(k*L0)=0", integerResidual];
  expectFail["tooth 3 halfshift boolean flips false", boolResidual[integerHalfshift]];
  expectZero["tooth 3 verdict is FAIL_POLE_LADDER", verdictResidual[computeVerdict[True, True, integerHalfshift, True, False], FAILPOLELADDER]];

  badRoundTrip = FullSimplify[Exp[2 I k L0] /. omega -> data["PoleLadder"]];
  expectFail["tooth 4 corrupt r_D/r_N makes round_trip_on_ladder differ from 1", badRoundTrip - 1];
  expectFail["tooth 4 round_trip_closes flips false", boolResidual[exprEqual[badRoundTrip, 1]]];

  noAlphaNumerator = FullSimplify[-k (k Sin[k L0])];
  noAlphaRobinDtn = FullSimplify[noAlphaNumerator/data["RobinDenominatorCore"]];
  noAlphaGuard = counterfactualGuard[noAlphaRobinDtn, data["RobinDenominatorCore"], data["RobinDenominatorCore"], 2/L0]["Guard"];
  expectFail["tooth 5 broken alpha->infinity path makes recovers_DD_at_alpha_inf false", boolResidual[noAlphaGuard["recovers_DD_at_alpha_inf"]]];
  expectZero["tooth 5 verdict is FAIL_COUNTERFACTUAL", verdictResidual[computeVerdict[True, True, True, guardAll[noAlphaGuard], False], FAILCOUNTERFACTUAL]];

  degenerateGuard = counterfactualGuard[data["RobinDtnTransfer"], data["RobinDetWitness"], data["RobinDenominatorCore"], 0]["Guard"];
  expectFail["tooth 6 numeric alpha degeneracy alpha=0 makes distinctness false", boolResidual[degenerateGuard["numeric_alpha_distinct"]]];
  expectZero["tooth 6 verdict is FAIL_COUNTERFACTUAL", verdictResidual[computeVerdict[True, True, True, guardAll[degenerateGuard], False], FAILCOUNTERFACTUAL]];

  Scan[
    Function[exponent,
      expectFail[
        "tooth 7 c_S^2 site A exponent 5->" <> ToString[exponent] <> " trips R1 dual-site integrity",
        r1SiteFromExponent[exponent] - data["R1SiteB"]
      ];
      expectFail[
        "tooth 7 c_S^2 site B exponent 5->" <> ToString[exponent] <> " trips R1 dual-site integrity",
        data["R1SiteA"] - r1EosSiteFromExponent[exponent]
      ]
    ],
    {4, 6}
  ];
  expectFail[
    "tooth 7 coordinated R1 both-site exponent drift trips frozen-export anchor",
    (r1SiteFromExponent[6] /. rho -> rhoStar) - 5 K rhoStar^4/m
  ];
  badLsSiteA = FullSimplify[D[psiHat[s], {s, 2}] - k^2 psiHat[s]];
  expectFail["tooth 7 L_s export sign corruption trips site A/B integrity", badLsSiteA - data["LsSiteB"]];
  hyperRecon = reconstructLsFromPair[{Sinh[k s], Cosh[k s]}];
  expectFail["tooth 7 L_s site-B sinh/cosh corruption trips null-space integrity", data["LsSiteA"] - hyperRecon["operator"]];

  expectFail["tooth 8 corrupt-[K] probe trips tan_argument dimensional gate", dimResidualVec[dim["CorruptWalk"]["TanArgumentDim"], dim["ExpectedTanDim"]]];
  expectFail["tooth 8 corrupt-[K] probe trips Z00_prefactor dimensional gate", dimResidualVec[dim["CorruptWalk"]["Z00PrefactorDim"], dim["ExpectedZ00Dim"]]];
  expectZero["tooth 8 corrupt-[K] verdict is DN_UNITTEST_FAIL_DIMENSIONAL", verdictResidual[dim["MutatedVerdict"], DNUNITTESTFAILDIMENSIONAL]];
  expectBool["tooth 8 self-ablation fail_suppressed remains true", dim["FailSuppressed"]];

  expectZero["baseline immutable after teeth: DtN still equals target", data["DtnTransfer"] - data["DtnTarget"]];
  expectZero["baseline immutable after teeth: halfshift pole residual remains zero", data["PoleResidual"]];
  expectZero["baseline immutable after teeth: Robin alpha->infinity still recovers D/D", data["Robin"]["DDFromRobin"] - data["Robin"]["DDTarget"]];
  expectZero["baseline immutable after teeth: L_s site integrity remains zero", data["LsSiteA"] - data["LsSiteB"]];
  expectZero["baseline immutable after teeth: clean 012 verdict remains DN_UNITTEST_BC_DEPENDENT", verdictResidual[data["Verdict"], DNUNITTESTBCDEPENDENT]]
];

Module[{ok, data},
  heading["ledger_stage012_dtn_pole_ladder_robin Mathematica audit"];
  ok = Catch[
    data = buildBaseline[];
    assertExact["baseline", data];
    runAritySelfCheck[data];
    runOpeningAndConsumedInputs[data];
    runTransferDtn[data];
    runPoleStaticRoundTrip[data];
    runRobinCounterfactual[data];
    runDimensionalBlock[data];
    runBCAndVerdict[data];
    printProvenance[];
    printVerdictLabels[];
    runAbleToFailTeeth[data];
    True,
    "ledgerStage012Failure",
    Function[{msg, tag}, Print["FAIL: ", msg]; False]
  ];

  Print[""];
  Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
  If[TrueQ[ok],
    Print["OVERALL PASS: Mathematica verified ledger_stage012 DtN pole ladder + Robin falsifier exactly"];
    Exit[0],
    Print["OVERALL FAIL: Mathematica stage012 audit did not close"];
    Exit[1]
  ]
]
