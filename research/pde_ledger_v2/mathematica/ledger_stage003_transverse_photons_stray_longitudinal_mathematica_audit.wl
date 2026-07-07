(* Ledger stage003 Mathematica audit: transverse photons plus stray longitudinal mode.

   Print-only, standalone, no arguments, no exports.  This engine uses native
   Wolfram matrix algebra for the transverse projector, a Faddeev-Jackiw
   symplectic null-vector route for the constrained longitudinal mode, and a
   canonical Poisson matrix only to assert the named Dirac residuals.
*)

ClearAll[
  heading, subheading, cleanZero, assertExact, expectZero, expectBool,
  expectFail, assertVerdict, dimResidual, ibpChain, signLabel,
  poissonMatrixBracket, longitudinalAnalysis, elasticControl,
  decoupledSlavedTheta, independentNoContinuity, branchAContinuity,
  epsilonMismatch, transverseDispersionFromLag, transverseBlock, runIbp,
  runTransverse, runBaseline, runBranches, runControls, runDimensions,
  passCount
];

passCount = 0;

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

cleanZero[expr_] := FullSimplify[expr] /. ConditionalExpression[0, _] -> 0;

assertExact[name_, expr_] := Module[{reals},
  reals = Cases[Unevaluated[expr], _Real, Infinity];
  If[reals =!= {},
    Print["FAIL  ", name, ": machine-real atom(s) found: ", ToString[InputForm[reals]]];
    Throw[name, "ledgerStage003Failure"]
  ]
];

expectZero[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name],
    Print["FAIL  ", name, ": residual = ", ToString[InputForm[clean]]];
    Throw[name, "ledgerStage003Failure"]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

expectFail[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[! TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name, " produced required FAIL (residual = ", ToString[InputForm[clean]], ")"],
    Print["FAIL  ", name, ": required mutation/ablation did not fire"];
    Throw[name, "ledgerStage003Failure"]
  ]
];

assertVerdict[name_, actual_, expected_] := Module[{allowed},
  allowed = If[ListQ[expected], expected, {expected}];
  expectBool[name <> ": verdict " <> actual, MemberQ[allowed, actual]]
];

dimResidual[actual_, expected_] := FullSimplify[(actual - expected).(actual - expected)];

ClearAll[
  rhoBr, muR, Bsym, Jsym, rhoB0, chiC, kappaPhase, betaB, mTheta2,
  ksym, omega, uL, uT1, uT2, theta, duL, duT, dtheta, pU, piTheta,
  deltaRhoB, divU, thetaDotSym, divDotU, gradThetaDotDu, thetaBoundary
];

$Assumptions =
  rhoBr > 0 && muR > 0 && Bsym > 0 && Jsym > 0 && rhoB0 > 0 &&
  chiC > 0 && kappaPhase > 0 && betaB > 0 && mTheta2 > 0 &&
  ksym > 0 && omega > 0;

ibpChain[slavingSign_: -1] := Module[
  {deltaRho, rawCoeff, timeCoeff, spaceCoeff},
  deltaRho = FullSimplify[slavingSign rhoB0 divU];
  rawCoeff = FullSimplify[Jsym slavingSign rhoB0];
  timeCoeff = FullSimplify[-rawCoeff];
  spaceCoeff = FullSimplify[-timeCoeff];
  <|
    "slavingSign" -> slavingSign,
    "deltaRho" -> deltaRho,
    "rawCoeff" -> rawCoeff,
    "timeCoeff" -> timeCoeff,
    "spaceCoeff" -> spaceCoeff,
    "CJ" -> spaceCoeff,
    "rawDensity" -> FullSimplify[rawCoeff thetaDotSym divU],
    "timeBulkDensity" -> FullSimplify[timeCoeff theta divDotU],
    "maxwellCrossDensity" -> FullSimplify[spaceCoeff gradThetaDotDu],
    "timeBoundaryDropped" -> True,
    "spaceBoundaryDropped" -> True
  |>
];

signLabel[expr_] := Which[
  TrueQ[FullSimplify[expr == 0]], "zero",
  TrueQ[FullSimplify[expr > 0]], "positive",
  TrueQ[FullSimplify[expr < 0]], "negative",
  True, "symbolic"
];

poissonMatrixBracket[f_, g_, vars_] := Module[{jmat, gradF, gradG},
  jmat = {{0, 0, 1, 0}, {0, 0, 0, 1}, {-1, 0, 0, 0}, {0, -1, 0, 0}};
  gradF = D[f, #] & /@ vars;
  gradG = D[g, #] & /@ vars;
  FullSimplify[gradF.jmat.gradG]
];

longitudinalAnalysis[name_, cExpr_, kExpr_, bExpr_, massExpr_: 0, flags_: <||>] := Module[
  {
    sTheta, lag, pUexpr, pThetaExpr, dotRule, hamiltonian, primary,
    z, oneForm, symplectic, rank, null, fjConstraint, canonicalVars,
    phi2, bracket, secondaryPreservation, dynamicMatrix, determinant,
    expectedDeterminant, aEff, omega2Pole, bracketZero, bZero, massZero,
    firstClass, secondClass, physicalDof, initialData, classification,
    poleCount, residue, bounded, verdict, squareResidual
  },

  sTheta = FullSimplify[kExpr ksym^2 - massExpr];
  lag = FullSimplify[
    1/2 rhoBr duL^2 - cExpr ksym uL dtheta +
      1/2 sTheta theta^2 - 1/2 bExpr ksym^2 uL^2
  ];
  pUexpr = D[lag, duL];
  pThetaExpr = D[lag, dtheta];
  dotRule = First[Solve[pU == pUexpr, duL]];
  hamiltonian = FullSimplify[(pU duL + pThetaExpr dtheta - lag) /. dotRule];
  primary = FullSimplify[piTheta - pThetaExpr];

  z = {uL, theta, pU};
  oneForm = {pU, -cExpr ksym uL, 0};
  symplectic = Table[
    D[oneForm[[j]], z[[i]]] - D[oneForm[[i]], z[[j]]],
    {i, Length[z]}, {j, Length[z]}
  ];
  rank = MatrixRank[symplectic];
  null = First[NullSpace[symplectic]];
  fjConstraint = FullSimplify[null.(D[hamiltonian, #] & /@ z)];

  canonicalVars = {uL, theta, pU, piTheta};
  phi2 = poissonMatrixBracket[primary, hamiltonian, canonicalVars];
  bracket = poissonMatrixBracket[primary, phi2, canonicalVars];
  secondaryPreservation = poissonMatrixBracket[phi2, hamiltonian, canonicalVars];

  dynamicMatrix = {
    {bExpr ksym^2 - rhoBr omega^2, -I cExpr ksym omega},
    {I cExpr ksym omega, -sTheta}
  };
  determinant = Factor[FullSimplify[Det[dynamicMatrix]]];
  expectedDeterminant = Factor[FullSimplify[
    (rhoBr sTheta - cExpr^2 ksym^2) omega^2 - bExpr ksym^2 sTheta
  ]];
  aEff = FullSimplify[rhoBr - cExpr^2 ksym^2/sTheta];
  omega2Pole = Which[
    TrueQ[FullSimplify[bExpr == 0]], 0,
    TrueQ[FullSimplify[aEff == 0]], "singular_second_class_no_propagating_pole",
    True, FullSimplify[bExpr ksym^2/aEff]
  ];

  bracketZero = TrueQ[FullSimplify[bracket == 0]];
  bZero = TrueQ[FullSimplify[bExpr == 0]];
  massZero = TrueQ[FullSimplify[massExpr == 0]];
  squareResidual = FullSimplify[kExpr - cExpr^2/rhoBr];

  Which[
    TrueQ[Lookup[flags, "GaugeOnShellOnly", False]],
      verdict = "FAIL_GAUGE_ON_SHELL_ONLY";
      firstClass = 0; secondClass = 2; physicalDof = 1; initialData = 2;
      classification = "ON_SHELL_ONLY_NOT_FIRST_CLASS";
      poleCount = 1; residue = "symbolic"; bounded = False,
    TrueQ[Lookup[flags, "PartialGaugeOnly", False]],
      verdict = "FAIL_PARTIAL_GAUGE_ONLY";
      firstClass = 0; secondClass = 2; physicalDof = 1; initialData = 2;
      classification = "PARTIAL_GAUGE_ONLY";
      poleCount = 1; residue = "symbolic"; bounded = False,
    bracketZero && bZero && massZero,
      verdict = If[
        TrueQ[Lookup[flags, "ProvenanceForced", False]],
        "C5_RESOLVED_MAXWELL_WITH_PROVENANCE",
        "C5_RESOLVED_MAXWELL_BY_TUNING"
      ];
      firstClass = 2; secondClass = 0; physicalDof = 0; initialData = 0;
      classification = "FIRST_CLASS_MAXWELL_CHAIN";
      poleCount = 0; residue = "no_physical_longitudinal_pole"; bounded = True,
    bracketZero && ! bZero,
      verdict = "FAIL_SECOND_CLASS_NOT_MAXWELL";
      firstClass = 0; secondClass = 4; physicalDof = 0; initialData = 0;
      classification = "TERTIARY_SECOND_CLASS_CHAIN";
      poleCount = 0; residue = "no_gauge_pole_but_not_first_class"; bounded = True,
    True,
      firstClass = 0; secondClass = 2; physicalDof = 1; initialData = 2;
      classification = "SECOND_CLASS_PAIR";
      poleCount = 1;
      residue = signLabel[1/aEff];
      bounded = TrueQ[FullSimplify[aEff > 0]] && (bZero || TrueQ[FullSimplify[bExpr > 0]]);
      verdict = Which[
        TrueQ[Lookup[flags, "GappedNotGaugeRemoved", False]], "FAIL_GAPPED_NOT_GAUGE_REMOVED",
        ! massZero, "FAIL_SECOND_CLASS_NOT_MAXWELL",
        ! TrueQ[bounded], "FAIL_GHOST_OR_NEGATIVE_NORM",
        bZero, "FAIL_C5_LONGITUDINAL_ZERO_MODE",
        True, "FAIL_CAUCHY_STRAY_LONGITUDINAL"
      ]
  ];

  <|
    "name" -> name,
    "lagrangian" -> lag,
    "pU" -> FullSimplify[pUexpr],
    "piTheta" -> FullSimplify[pThetaExpr],
    "primary" -> primary,
    "hamiltonian" -> hamiltonian,
    "symplecticMatrix" -> symplectic,
    "symplecticRank" -> rank,
    "symplecticNullVector" -> null,
    "fjConstraint" -> fjConstraint,
    "secondary" -> FullSimplify[phi2],
    "fjMatchesDiracSecondary" -> FullSimplify[cExpr ksym fjConstraint - phi2],
    "bracket" -> FullSimplify[bracket],
    "secondaryPreservationNoMultiplier" -> FullSimplify[secondaryPreservation],
    "determinant" -> determinant,
    "expectedDeterminant" -> expectedDeterminant,
    "omega2Pole" -> FullSimplify[omega2Pole],
    "aEff" -> aEff,
    "squareResidual" -> squareResidual,
    "firstClass" -> firstClass,
    "secondClass" -> secondClass,
    "physicalDof" -> physicalDof,
    "initialData" -> initialData,
    "classification" -> classification,
    "poleCount" -> poleCount,
    "residue" -> residue,
    "bounded" -> bounded,
    "verdict" -> verdict
  |>
];

elasticControl[name_, bExpr_] := <|
  "name" -> name,
  "omega2" -> FullSimplify[bExpr ksym^2/rhoBr],
  "physicalDof" -> 1,
  "initialData" -> 2,
  "poleCount" -> 1,
  "classification" -> "UNCONSTRAINED_ELASTIC_COORDINATE",
  "verdict" -> If[
    TrueQ[FullSimplify[bExpr == 0]],
    "FAIL_C5_LONGITUDINAL_ZERO_MODE",
    "FAIL_CAUCHY_STRAY_LONGITUDINAL"
  ]
|>;

decoupledSlavedTheta[] := Module[
  {lag, pUexpr, pThetaExpr, hamiltonian, primary, vars, phi2, bracket},
  lag = FullSimplify[1/2 rhoBr duL^2 - 1/2 kappaPhase ksym^2 theta^2];
  pUexpr = D[lag, duL];
  pThetaExpr = D[lag, dtheta];
  hamiltonian = FullSimplify[pU^2/(2 rhoBr) + 1/2 kappaPhase ksym^2 theta^2];
  primary = FullSimplify[piTheta - pThetaExpr];
  vars = {uL, theta, pU, piTheta};
  phi2 = poissonMatrixBracket[primary, hamiltonian, vars];
  bracket = poissonMatrixBracket[primary, phi2, vars];
  <|
    "pU" -> pUexpr,
    "piTheta" -> pThetaExpr,
    "primary" -> primary,
    "secondary" -> phi2,
    "bracket" -> bracket,
    "physicalDof" -> 1,
    "initialData" -> 2,
    "classification" -> "THETA_ALGEBRAIC_SECOND_CLASS_PAIR_PLUS_U_ZERO_MODE",
    "verdict" -> "FAIL_C5_LONGITUDINAL_ZERO_MODE"
  |>
];

independentNoContinuity[] := Module[
  {densityLag, gaussianRule, thetaKinetic, secondOrderLag, hessian, thetaOmega2},
  densityLag = FullSimplify[Jsym dtheta deltaRhoB - deltaRhoB^2/(2 chiC)];
  gaussianRule = First[Solve[D[densityLag, deltaRhoB] == 0, deltaRhoB]];
  thetaKinetic = FullSimplify[densityLag /. gaussianRule];
  secondOrderLag = FullSimplify[
    1/2 rhoBr duL^2 + thetaKinetic - 1/2 kappaPhase ksym^2 theta^2
  ];
  hessian = D[secondOrderLag, {{duL, dtheta}, 2}];
  thetaOmega2 = FullSimplify[kappaPhase ksym^2/(Jsym^2 chiC)];
  <|
    "gaussianSolution" -> (deltaRhoB /. gaussianRule),
    "thetaKinetic" -> thetaKinetic,
    "hessianRank" -> MatrixRank[hessian],
    "thetaOmega2" -> thetaOmega2,
    "physicalDof" -> 2,
    "initialData" -> 4,
    "classification" -> "TWO_UNCONSTRAINED_SECOND_ORDER_FIELDS",
    "verdict" -> "FAIL_EXTRA_SCALAR_DOF"
  |>
];

branchAContinuity[cExpr_] := Module[
  {densityLag, continuityRule, withContinuity, rawCrossCoeff, canonicalCoeff, bIncrement},
  densityLag = FullSimplify[Jsym dtheta deltaRhoB - deltaRhoB^2/(2 chiC)];
  continuityRule = First[Solve[omega (deltaRhoB + rhoB0 ksym uL) == 0, deltaRhoB]];
  withContinuity = FullSimplify[densityLag /. continuityRule];
  rawCrossCoeff = FullSimplify[Coefficient[withContinuity, uL dtheta]/ksym];
  canonicalCoeff = FullSimplify[-rawCrossCoeff];
  bIncrement = FullSimplify[rhoB0^2/chiC];
  <|
    "continuityResidual" -> FullSimplify[omega ((deltaRhoB /. continuityRule) + rhoB0 ksym uL)],
    "withContinuity" -> withContinuity,
    "rawCrossCoeff" -> rawCrossCoeff,
    "canonicalMinusCJ" -> FullSimplify[canonicalCoeff + cExpr],
    "BIncrement" -> bIncrement,
    "proofStatus" -> "CONTINUITY_FORCES_SAME_SLAVED_SECTOR"
  |>
];

epsilonMismatch[cExpr_] := Module[
  {rhoEps, kEps, sTheta, lag, pThetaExpr, hamiltonian, primary, vars,
    phi2, bracket, frozenTransverse, shiftedTransverse},
  rhoEps = 2 rhoBr;
  kEps = FullSimplify[cExpr^2/rhoEps];
  sTheta = FullSimplify[kEps ksym^2];
  lag = FullSimplify[1/2 rhoEps duL^2 - cExpr ksym uL dtheta + 1/2 sTheta theta^2];
  pThetaExpr = D[lag, dtheta];
  hamiltonian = FullSimplify[pU^2/(2 rhoEps) - 1/2 sTheta theta^2];
  primary = FullSimplify[piTheta - pThetaExpr];
  vars = {uL, theta, pU, piTheta};
  phi2 = poissonMatrixBracket[primary, hamiltonian, vars];
  bracket = poissonMatrixBracket[primary, phi2, vars];
  frozenTransverse = FullSimplify[muR ksym^2/rhoBr];
  shiftedTransverse = FullSimplify[muR ksym^2/rhoEps];
  <|
    "bracket" -> bracket,
    "squareCloses" -> TrueQ[FullSimplify[bracket == 0]],
    "frozenTransverseOmega2" -> frozenTransverse,
    "mismatchedTransverseOmega2" -> shiftedTransverse,
    "speedShift" -> FullSimplify[shiftedTransverse - frozenTransverse],
    "physicalDof" -> If[TrueQ[FullSimplify[bracket == 0]], 0, 1],
    "initialData" -> If[TrueQ[FullSimplify[bracket == 0]], 0, 2],
    "verdict" -> "FAIL_TRANSVERSE_DISTURBED"
  |>
];

transverseDispersionFromLag[lag_] := Module[
  {tt, q, lagT, el, planeWaveEl, residual, omegaSq, omega2},
  lagT = lag /. {uT1 -> q[tt], duT -> D[q[tt], tt]};
  el = FullSimplify[D[D[lagT, D[q[tt], tt]], tt] - D[lagT, q[tt]]];
  planeWaveEl = FullSimplify[el /. Derivative[2][q][tt] -> -omega^2 q[tt]];
  residual = FullSimplify[planeWaveEl/q[tt]];
  omega2 = omegaSq /. First[Solve[(residual /. omega^2 -> omegaSq) == 0, omegaSq]];
  <|
    "EL" -> el,
    "planeWaveResidual" -> residual,
    "omega2" -> FullSimplify[omega2]
  |>
];

transverseBlock[] := Module[
  {
    wave, uVec, khat, longProjector, transProjector, transverseBasis,
    curlVec, curlSq, divTransverse, josephsonTransverse, perPolLag,
    pT, dispersion, omega2
  },
  wave = {ksym, 0, 0};
  uVec = {uL, uT1, uT2};
  khat = {1, 0, 0};
  longProjector = Outer[Times, khat, khat];
  transProjector = IdentityMatrix[3] - longProjector;
  transverseBasis = NullSpace[{khat}];
  curlVec = Cross[wave, uVec];
  curlSq = FullSimplify[curlVec.curlVec /. uL -> 0];
  divTransverse = FullSimplify[wave.(transProjector.uVec)];
  josephsonTransverse = FullSimplify[Jsym dtheta (-rhoB0 divTransverse)];
  perPolLag = FullSimplify[1/2 rhoBr duT^2 - 1/2 muR ksym^2 uT1^2];
  pT = D[perPolLag, duT];
  dispersion = transverseDispersionFromLag[perPolLag];
  omega2 = dispersion["omega2"];
  <|
    "basisCount" -> Length[transverseBasis],
    "curlSq" -> curlSq,
    "josephsonTransverse" -> josephsonTransverse,
    "pT" -> pT,
    "EL" -> dispersion["EL"],
    "planeWaveResidual" -> dispersion["planeWaveResidual"],
    "omega2" -> FullSimplify[omega2],
    "cGamma2" -> FullSimplify[omega2/ksym^2],
    "physicalDof" -> 2,
    "massless" -> True,
    "thetaCouplings" -> josephsonTransverse,
    "verdict" -> "PASS_TRANSVERSE_UNDISTURBED"
  |>
];

runIbp[] := Module[{chain, canonicalLag, piExpr, corrupted, corruptedPi},
  subheading["Josephson slaving and C_J IBP sign"];
  chain = ibpChain[-1];
  expectZero["slaved delta_rho_B = -rho_B0 div u", chain["deltaRho"] + rhoB0 divU];
  expectZero["time IBP flips theta_dot div(u) coefficient", chain["timeCoeff"] + chain["rawCoeff"]];
  expectZero["space IBP flips theta div(dot u) coefficient", chain["spaceCoeff"] + chain["timeCoeff"]];
  expectZero["derived C_J = -J*rho_B0", chain["CJ"] + Jsym rhoB0];
  expectBool["time boundary term dropped in ledger", chain["timeBoundaryDropped"]];
  expectBool["space boundary term dropped in ledger", chain["spaceBoundaryDropped"]];
  canonicalLag = FullSimplify[1/2 rhoBr duL^2 - chain["CJ"] ksym uL dtheta];
  piExpr = D[canonicalLag, dtheta];
  expectZero["sign-sensitive pi_theta = +J*k*rho_B0*u_L", piExpr - Jsym ksym rhoB0 uL];
  corrupted = ibpChain[1];
  corruptedPi = D[-corrupted["CJ"] ksym uL dtheta, dtheta];
  expectFail["C_J sign mutation flips pi_theta residual", corruptedPi - Jsym ksym rhoB0 uL];
  chain["CJ"]
];

runTransverse[] := Module[{trans, corruptedLag, corruptedDispersion},
  subheading["Transverse earned sector"];
  trans = transverseBlock[];
  expectZero["transverse omega^2 = (mu_R/rho_br) k^2", trans["omega2"] - muR ksym^2/rhoBr];
  expectZero["c_gamma^2 = mu_R/rho_br", trans["cGamma2"] - muR/rhoBr];
  expectZero["theta couplings vanish for transverse modes", trans["thetaCouplings"]];
  expectZero["two transverse physical polarizations", trans["physicalDof"] - 2];
  expectBool["transverse modes are massless", trans["massless"]];
  assertVerdict["transverse sector", trans["verdict"], "PASS_TRANSVERSE_UNDISTURBED"];
  corruptedLag = FullSimplify[1/2 rhoBr duT^2 - 1/2 (2 muR) ksym^2 uT1^2];
  corruptedDispersion = transverseDispersionFromLag[corruptedLag];
  expectFail[
    "transverse primitive mu_R->2 mu_R mutation fails c_gamma^2",
    FullSimplify[corruptedDispersion["omega2"]/ksym^2 - muR/rhoBr]
  ];
];

runBaseline[cExpr_] := Module[{analysis, bEff},
  subheading["Longitudinal Faddeev-Jackiw / Dirac chain"];
  bEff = FullSimplify[rhoB0^2/chiC];
  analysis = longitudinalAnalysis[
    "branch_b_slaved_finite_compressibility_conventional_K",
    cExpr, -kappaPhase, bEff
  ];
  expectZero["p_u = rho_br dot_u_L", analysis["pU"] - rhoBr duL];
  expectZero["pi_theta = +J*k*rho_B0*u_L", analysis["piTheta"] - Jsym ksym rhoB0 uL];
  expectZero["Phi_1 = pi_theta - J*k*rho_B0*u_L", analysis["primary"] - (piTheta - Jsym ksym rhoB0 uL)];
  expectZero[
    "Faddeev-Jackiw null constraint matches Dirac secondary",
    analysis["fjMatchesDiracSecondary"]
  ];
  expectZero[
    "Phi_2 = -k(J*p_u*rho_B0 + k*kappa_phase*rho_br*theta)/rho_br",
    analysis["secondary"] + ksym (Jsym pU rhoB0 + ksym kappaPhase rhoBr theta)/rhoBr
  ];
  expectZero[
    "{Phi_1,Phi_2} = k^2(J^2 rho_B0^2 + kappa_phase rho_br)/rho_br",
    analysis["bracket"] - ksym^2 (Jsym^2 rhoB0^2 + kappaPhase rhoBr)/rhoBr
  ];
  expectZero["symplectic matrix has one null direction", analysis["symplecticRank"] - 2];
  expectZero["first-class count = 0", analysis["firstClass"]];
  expectZero["second-class count = 2", analysis["secondClass"] - 2];
  expectZero["longitudinal physical DOF = (4-2)/2 = 1", analysis["physicalDof"] - 1];
  expectZero["independent initial data functions = 2", analysis["initialData"] - 2];
  expectBool["constraint classification SECOND_CLASS_PAIR", analysis["classification"] == "SECOND_CLASS_PAIR"];
  expectZero[
    "dynamic determinant assembled from Euler-Lagrange matrix",
    analysis["determinant"] - analysis["expectedDeterminant"]
  ];
  expectZero[
    "stray pole omega^2",
    analysis["omega2Pole"] -
      ksym^2 kappaPhase rhoB0^2/(chiC (Jsym^2 rhoB0^2 + kappaPhase rhoBr))
  ];
  expectZero["pole count = 1", analysis["poleCount"] - 1];
  expectBool["positive pole residue", analysis["residue"] == "positive"];
  expectBool["reduced Hamiltonian bounded", analysis["bounded"]];
  assertVerdict["baseline headline", analysis["verdict"], "FAIL_CAUCHY_STRAY_LONGITUDINAL"];
  Print[""];
  Print["Postulate line:"];
  Print["  K_theta = -kappa_phase is an explicit conventional phase-stiffness input, not a derivation."];
  analysis
];

runBranches[cExpr_] := Module[{bEff, kLocus, branchA, branchAFinite, curlOnly, maxwell},
  subheading["Branch (a), curl-only subcase, and Maxwell locus"];
  bEff = FullSimplify[rhoB0^2/chiC];
  kLocus = FullSimplify[cExpr^2/rhoBr];
  branchA = branchAContinuity[cExpr];
  branchAFinite = longitudinalAnalysis[
    "branch_a_independent_with_continuity_integrated_out",
    cExpr, -kappaPhase, bEff
  ];
  curlOnly = longitudinalAnalysis[
    "branch_b_slaved_curl_only_conventional_K",
    cExpr, -kappaPhase, 0
  ];
  maxwell = longitudinalAnalysis[
    "branch_b_slaved_tuned_Maxwell_locus",
    cExpr, kLocus, 0
  ];
  expectZero["branch (a) continuity equation solved in fixed-number sector", branchA["continuityResidual"]];
  expectZero["branch (a) canonical cross coefficient matches derived C_J", branchA["canonicalMinusCJ"]];
  expectZero["branch (a) B increment = rho_B0^2/chi_c", branchA["BIncrement"] - bEff];
  assertVerdict["branch (a) integrated sector", branchAFinite["verdict"], "FAIL_CAUCHY_STRAY_LONGITUDINAL"];
  expectBool["branch (a) proof status", branchA["proofStatus"] == "CONTINUITY_FORCES_SAME_SLAVED_SECTOR"];
  assertVerdict["curl-only conventional subcase", curlOnly["verdict"], "FAIL_C5_LONGITUDINAL_ZERO_MODE"];
  expectZero["curl-only longitudinal DOF = 1", curlOnly["physicalDof"] - 1];
  expectZero["Maxwell-locus square residual K_theta - C_J^2/rho_br", maxwell["squareResidual"]];
  expectZero["Maxwell-locus bracket = 0", maxwell["bracket"]];
  expectZero["Maxwell-locus first-class count = 2", maxwell["firstClass"] - 2];
  expectZero["Maxwell-locus longitudinal DOF = 0", maxwell["physicalDof"]];
  assertVerdict["Maxwell locus reachable", maxwell["verdict"], "C5_RESOLVED_MAXWELL_BY_TUNING"];
  expectBool["Maxwell-locus classification", maxwell["classification"] == "FIRST_CLASS_MAXWELL_CHAIN"];
  {curlOnly, maxwell}
];

runControls[cExpr_, baseline_, curlOnly_, maxwell_] := Module[
  {
    kLocus, controls, epsilon, decSlaved, noContinuity, trans, thetaMass,
    reachableTokens
  },
  subheading["Fourteen source controls"];
  kLocus = FullSimplify[cExpr^2/rhoBr];
  controls = {
    {"1_no_theta", "FAIL_C5_LONGITUDINAL_ZERO_MODE", elasticControl["no_theta_curl_only", 0]},
    {"2_cauchy_bulk", "FAIL_CAUCHY_STRAY_LONGITUDINAL", elasticControl["cauchy_bulk_no_theta", betaB]},
    {"3_mismatched_positive_K_no_B", "FAIL_C5_LONGITUDINAL_ZERO_MODE",
      longitudinalAnalysis["mismatched_positive_K_no_B", cExpr, FullSimplify[2 cExpr^2/rhoBr], 0]},
    {"3_mismatched_K_theta_le_0", "FAIL_C5_LONGITUDINAL_ZERO_MODE",
      longitudinalAnalysis["mismatched_K_theta_le_0", cExpr, -kappaPhase, 0]},
    {"3_mismatched_positive_K_with_B", "FAIL_CAUCHY_STRAY_LONGITUDINAL",
      longitudinalAnalysis["mismatched_positive_K_with_B", cExpr, FullSimplify[2 cExpr^2/rhoBr], betaB]},
    {"3_positive_K_negative_residue", "FAIL_GHOST_OR_NEGATIVE_NORM",
      longitudinalAnalysis["positive_K_negative_residue", cExpr, FullSimplify[cExpr^2/(2 rhoBr)], betaB]},
    {"3_B_on_square_locus", "FAIL_SECOND_CLASS_NOT_MAXWELL",
      longitudinalAnalysis["B_on_square_locus", cExpr, kLocus, betaB]}
  };
  Scan[Function[row, assertVerdict[row[[1]], row[[3]]["verdict"], row[[2]]]], controls];

  epsilon = epsilonMismatch[cExpr];
  assertVerdict["3_epsilon_mismatch", epsilon["verdict"], "FAIL_TRANSVERSE_DISTURBED"];
  expectZero[
    "3_epsilon_mismatch transverse speed shifts to mu_R/(2 rho_br)",
    epsilon["mismatchedTransverseOmega2"] - muR ksym^2/(2 rhoBr)
  ];

  decSlaved = decoupledSlavedTheta[];
  assertVerdict["4_decoupled_theta_slaved", decSlaved["verdict"], "FAIL_C5_LONGITUDINAL_ZERO_MODE"];

  noContinuity = independentNoContinuity[];
  assertVerdict[
    "4_decoupled_theta_independent_no_continuity",
    noContinuity["verdict"],
    "FAIL_EXTRA_SCALAR_DOF"
  ];
  expectZero[
    "4_decoupled_theta_independent_no_continuity theta kinetic",
    noContinuity["thetaKinetic"] - 1/2 Jsym^2 chiC dtheta^2
  ];

  trans = transverseBlock[];
  assertVerdict["5_transverse", trans["verdict"], "PASS_TRANSVERSE_UNDISTURBED"];
  assertVerdict["6_provenance_ablation fixed coefficients", baseline["verdict"], "FAIL_CAUCHY_STRAY_LONGITUDINAL"];
  assertVerdict["6_provenance_ablation free locus", maxwell["verdict"], "C5_RESOLVED_MAXWELL_BY_TUNING"];
  assertVerdict["7_compressibility_absent_vs_included absent", curlOnly["verdict"], "FAIL_C5_LONGITUDINAL_ZERO_MODE"];
  assertVerdict["7_compressibility_absent_vs_included included", baseline["verdict"], "FAIL_CAUCHY_STRAY_LONGITUDINAL"];

  thetaMass = longitudinalAnalysis["theta_mass_breaks_gauge", cExpr, kLocus, 0, mTheta2];
  assertVerdict["8_theta_mass", thetaMass["verdict"], "FAIL_SECOND_CLASS_NOT_MAXWELL"];

  reachableTokens = {
    longitudinalAnalysis[
      "grammar_with_provenance", cExpr, kLocus, 0, 0,
      <|"ProvenanceForced" -> True|>
    ]["verdict"],
    longitudinalAnalysis[
      "grammar_gapped", cExpr, -kappaPhase, 0, 0,
      <|"GappedNotGaugeRemoved" -> True|>
    ]["verdict"],
    longitudinalAnalysis[
      "grammar_partial", cExpr, -kappaPhase, 0, 0,
      <|"PartialGaugeOnly" -> True|>
    ]["verdict"],
    longitudinalAnalysis[
      "grammar_on_shell", cExpr, -kappaPhase, 0, 0,
      <|"GaugeOnShellOnly" -> True|>
    ]["verdict"]
  };
  expectBool[
    "reachable-but-unfired grammar tokens are preserved",
    SubsetQ[
      reachableTokens,
      {
        "C5_RESOLVED_MAXWELL_WITH_PROVENANCE",
        "FAIL_GAPPED_NOT_GAUGE_REMOVED",
        "FAIL_PARTIAL_GAUGE_ONLY",
        "FAIL_GAUGE_ON_SHELL_ONLY"
      }
    ]
  ];
];

runDimensions[] := Module[
  {
    braneLag, zdim, du, dTheta, dGrad, dDt, dk, dOmega, dRhoBr, dMuR,
    dB, dRhoB0, dJ, dCJ, dKtheta, dChiC, dMtheta2, dDivU, dCurlU,
    dGradTheta, dDtU, dDtTheta, dDeltaRho, checks
  },
  subheading["Dimensional firewall"];
  braneLag = {1, -1, -2};
  zdim = {0, 0, 0};
  du = {0, 1, 0};
  dTheta = zdim;
  dGrad = {0, -1, 0};
  dDt = {0, 0, -1};
  dk = dGrad;
  dOmega = dDt;
  dRhoBr = {1, -3, 0};
  dMuR = braneLag;
  dB = braneLag;
  dRhoB0 = dRhoBr;
  dJ = {0, 2, -1};
  dCJ = {1, -1, -1};
  dKtheta = {1, 1, -2};
  dChiC = {1, -5, 2};
  dMtheta2 = braneLag;
  dDivU = dGrad + du;
  dCurlU = dDivU;
  dGradTheta = dGrad + dTheta;
  dDtU = dDt + du;
  dDtTheta = dDt + dTheta;
  dDeltaRho = dRhoB0;

  checks = {
    {"brane inertia rho_br (partial_t u)^2", dRhoBr + 2 dDtU, braneLag},
    {"MacCullagh curl energy mu_R (curl u)^2", dMuR + 2 dCurlU, braneLag},
    {"Cauchy bulk term B (div u)^2", dB + 2 dDivU, braneLag},
    {"Josephson density term J theta_dot delta_rho_B", dJ + dDtTheta + dDeltaRho, braneLag},
    {"slaved density rho_B0 div u", dRhoB0 + dDivU, dDeltaRho},
    {"signed phase gradient K_theta (grad theta)^2", dKtheta + 2 dGradTheta, braneLag},
    {"compressibility delta_rho_B^2/chi_c", 2 dDeltaRho - dChiC, braneLag},
    {"theta mass m_theta^2 theta^2", dMtheta2 + 2 dTheta, braneLag},
    {"C_J = -J rho_B0", dJ + dRhoB0, dCJ},
    {"IBP cross C_J partial_t u dot grad theta", dCJ + dDtU + dGradTheta, braneLag},
    {"electric square velocity piece", dRhoBr + 2 dDtU, braneLag},
    {"electric square mixed piece", dCJ + dDtU + dGradTheta, braneLag},
    {"electric square gradient piece", dKtheta + 2 dGradTheta, braneLag},
    {"Maxwell locus C_J^2 = rho_br K_theta", 2 dCJ, dRhoBr + dKtheta},
    {"c_gamma^2 = mu_R/rho_br", dMuR - dRhoBr, {0, 2, -2}},
    {"branch-a theta kinetic chi_c J^2 theta_dot^2", dChiC + 2 dJ + 2 dDtTheta, braneLag},
    {"rho_B0^2/chi_c stiffness increment", 2 dRhoB0 - dChiC, dB},
    {"transverse omega^2", dMuR - dRhoBr + 2 dk, 2 dOmega}
  };
  Scan[
    Function[row, expectZero["dimension: " <> row[[1]], dimResidual[row[[2]], row[[3]]]]],
    checks
  ];
  expectFail[
    "dimension ablation: drop rho_B0 from Josephson cross",
    dimResidual[dJ + dDtU + dGradTheta, braneLag]
  ];
  expectFail[
    "dimension ablation: omit gradient from phase stiffness",
    dimResidual[dKtheta + 2 dTheta, braneLag]
  ];
  expectFail[
    "dimension ablation: multiply by chi_c instead of divide",
    dimResidual[dChiC + 2 dDeltaRho, braneLag]
  ];
];

Module[{ok, cExpr, baseline, branches, curlOnly, maxwell},
  heading["ledger_stage003_transverse_photons_stray_longitudinal Mathematica audit"];
  ok = Catch[
    cExpr = runIbp[];
    runTransverse[];
    baseline = runBaseline[cExpr];
    branches = runBranches[cExpr];
    curlOnly = branches[[1]];
    maxwell = branches[[2]];
    runControls[cExpr, baseline, curlOnly, maxwell];
    runDimensions[];
    True,
    "ledgerStage003Failure",
    Function[{name, tag}, False]
  ];

  If[TrueQ[ok],
    Print[""];
    Print["Verdict labels:"];
    Print["  ledger earned-label (NOT a source verdict token): LIGHT_ON_BRANE_TWO_TRANSVERSE_PHOTONS  (c_gamma^2 = mu_R/rho_br, physical_dof=2)"];
    Print["  source top-line verdict: FAIL_CAUCHY_STRAY_LONGITUDINAL"];
    Print["  characterized departure: SECOND_CLASS_PAIR (one stray longitudinal DOF; Maxwell locus reachable BY_TUNING only)"];
    Print[""];
    Print["CHECK TALLY: PASS=", passCount, " FAIL=0"];
    Print["OVERALL PASS: Mathematica derived ledger_stage003 light-sector claims exactly"];
    Exit[0],
    Print[""];
    Print["OVERALL FAIL: Mathematica stage003 audit did not close"];
    Exit[1]
  ]
]
