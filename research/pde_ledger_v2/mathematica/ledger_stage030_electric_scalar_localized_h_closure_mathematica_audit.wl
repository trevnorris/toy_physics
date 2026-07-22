(* Ledger stage030 Mathematica audit: electric scalar + localized-H closure.

   Standalone, print-only, no arguments, no imports or exports.  This is an
   independent Wolfram route: Integrate derives N0, DSolve derives ker(A),
   Limit derives the continuum threshold, and the coupled kernel is checked
   with PositiveDefiniteMatrixQ, CharacteristicPolynomial, and Eigenvalues.

   The optional internal ablation switch is the environment variable
   LEDGER_STAGE030_MUTATION.  Accepted values are the predicate names in the
   ablations association below.
*)

ClearAll[
  raise, heading, subheading, assertExact, cleanZero, expectZero, expectBool,
  primitive, dimResidual, matrixZeroQ, positiveExactQ, passCount, failCount,
  mutationEnvironment, activeMutation, ablations
];

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE030_MUTATION";

raise[name_] := Throw[name, "ledgerStage030Failure"];

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

assertExact[name_, expression_] := Module[{reals},
  reals = Cases[Unevaluated[expression], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FIRST_FAILURE=", name];
    Print["FAIL  ", name, ": machine-real atom(s) found: ",
      ToString[InputForm[reals]]];
    raise[name]
  ]
];

cleanZero[expression_] := FullSimplify[expression] /.
  ConditionalExpression[0, _] -> 0;

expectZero[name_, residual_, evidence_: None] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FIRST_FAILURE=", name];
    Print["FAIL  ", name, ": residual = ", ToString[InputForm[clean]]];
    If[evidence =!= None,
      Print["      evidence = ", ToString[InputForm[evidence]]]
    ];
    raise[name]
  ]
];

expectBool[name_, condition_, evidence_: None] :=
  expectZero[name, If[TrueQ[condition], 0, 1], evidence];

ablations = <|
  "PASS_TRANSVERSE_FACTORIZATION" -> <|
    "Primitive" -> "factor_A_coefficient", "Value" -> 3,
    "Description" -> "A superpotential coefficient: 2 -> 3"|>,
  "PASS_PARENT_TRANSVERSE_ZERO_EIGENVALUE" -> <|
    "Primitive" -> "f0_power", "Value" -> 1,
    "Description" -> "f0 sech exponent: 2 -> 1"|>,
  "PASS_UNIQUE_NORMALIZABLE_KERNEL" -> <|
    "Primitive" -> "kernel_derivative_coefficient", "Value" -> 2,
    "Description" -> "coefficient of d/dw in A f=0: 1 -> 2"|>,
  "PASS_POSITIVE_CONTINUUM_GAP" -> <|
    "Primitive" -> "V_asymptote_coefficient", "Value" -> 0,
    "Description" -> "V_H asymptote coefficient: 4 -> 0"|>,
  "PASS_ZERO_MODE_NORM" -> <|
    "Primitive" -> "norm_inner_weight", "Value" -> 3,
    "Description" -> "parent-H norm weight: 2 -> 3"|>,
  "PASS_NORMALIZED_ZERO_MODE_PROJECTION" -> <|
    "Primitive" -> "projection_inner_weight", "Value" -> 3,
    "Description" -> "projection weight in h=P0 H: 2 -> 3"|>,
  "PASS_PHYSICAL_H_NORM" -> <|
    "Primitive" -> "physical_M4_star", "Value" -> -3/160,
    "Description" -> "M4*: 3/160 -> -3/160"|>,
  "PASS_REDUCED_H_INERTIA" -> <|
    "Primitive" -> "inertia_M4_star", "Value" -> 3/320,
    "Description" -> "M4*: 3/160 -> 3/320 in M_h=N0 M4"|>,
  "PASS_PARENT_STIFFNESS_DIMENSIONAL_CONSISTENCY" -> <|
    "Primitive" -> "parent_K4_dimension", "Value" -> {1, -2, 1},
    "Description" -> "[K4]: M L^2 T^-2 -> M L T^-2"|>,
  "PASS_REDUCED_H_STIFFNESS" -> <|
    "Primitive" -> "declared_Kh_star", "Value" -> 2,
    "Description" -> "declared K_h*: 1 -> 2"|>,
  "PASS_REDUCED_SPEED_PRESERVATION" -> <|
    "Primitive" -> "speed_Kh_N0_power", "Value" -> 2,
    "Description" -> "K_h reduction: N0 K4 -> N0^2 K4"|>,
  "PASS_STABILITY" -> <|
    "Primitive" -> "stability_Chu_star", "Value" -> 2,
    "Description" -> "C_hu*: 1/2 -> 2"|>,
  "PASS_POSITIVE_GENERALIZED_WAVE_SPEEDS" -> <|
    "Primitive" -> "wave_Aeff_star", "Value" -> 2,
    "Description" ->
      "A_eff*: 1 -> 2 while the Sylvester margin stays positive"|>,
  "PASS_REDUCED_H_MASSLESSNESS" -> <|
    "Primitive" -> "reduced_h_mass_star", "Value" -> 1,
    "Description" -> "reduced k^0 h^2 coefficient*: 0 -> 1"|>,
  "PASS_CONSERVATIVE_HESSIAN_SYMMETRY" -> <|
    "Primitive" -> "u_equation_cross_star", "Value" -> 3/4,
    "Description" -> "u_L Euler cross coefficient*: 1/2 -> 3/4"|>,
  "PASS_DIMENSIONAL_HOMOGENEITY" -> <|
    "Primitive" -> "firewall_Chu_dimension", "Value" -> {0, -1, 1},
    "Description" -> "[C_hu] in mix: M T^-2 -> M T^-1"|>
|>;

activeMutation = Environment[mutationEnvironment];
If[! StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];

primitive[predicate_String, name_String, cardValue_] := If[
  activeMutation === predicate &&
    Lookup[ablations[predicate], "Primitive"] === name,
  Lookup[ablations[predicate], "Value"],
  cardValue
];

dimResidual[actual_List, expected_List] :=
  FullSimplify[(actual - expected).(actual - expected)];
matrixZeroQ[matrix_List] := And @@ (TrueQ[# === 0] & /@ Flatten[matrix]);
positiveExactQ[expression_] := TrueQ[FullSimplify[expression > 0]];

Module[
  {
    ok,
    dimZero = {0, 0, 0}, dimL = {1, 0, 0}, dimT = {0, 1, 0},
    dimMass = {0, 0, 1}, dimEnergy, dimVelocity, dimH, dimh, dimU,
    dimN0, dimM4, dimK4, dimAEff, dimBEff, dimMH, dimKH, dimCHU,
    dimCE, bulkDensity, braneDensity, dimD,
    bEffStar, kHStar, stabilityChuStar, stiffnessStar, dStar,
    stiffnessEigenvalues, stabilityCondition, aEffStar, mHStar,
    inertiaStar, waveStiffnessStar, dynamicMatrix, directDeterminant,
    characteristicPolynomial, waveRoots, expectedRoots, zStar,
    unchangedMargin, waveCondition,
    w, ell, f0, canonicalN0, normWeight, normN0, projectionWeight,
    p0f0, projectedH, idempotenceResidual, h0, ellStar, n0Star,
    baseM4Star, physicalM4Star, physicalNormStar, inertiaM4Star,
    reducedMassStar, cEStar, k4DerivedStar, parentK4Dimension,
    speedKhN0Power, parentSpeedSquaredStar, reducedSpeedSquaredStar,
    khProjectedStar, khDeclaredStar,
    parentPotential, factorCoefficient, factorSuperpotential, testField,
    aOnTest, adagAOnTest, operatorOnTest, factorResidual, f0Power,
    zeroModeTest, zeroModeResidual, kernelDerivativeCoefficient,
    kernelEquation, kernelRules, kernelShape, kernelConstants,
    normalizedKernelShape, kernelRatio, kernelResidual, kernelNorm,
    asymptoteCoefficient, asymptoticPotential, gapPlus, gapMinus, gap,
    omega, k, reducedHMassStar, qScalar, qAtOrigin, gradU, gradH,
    conservativeDensity, actionHessian, uCrossStar, eulerOperator,
    firewallCHU, dimensionTerms, dimensionTargets, dimensionResiduals
  },

  dimEnergy = 2 dimL - 2 dimT + dimMass;
  dimVelocity = dimL - dimT;
  dimH = -dimL;
  dimh = dimZero;
  dimU = dimL;
  dimN0 = -dimL;
  dimM4 = dimMass;
  dimK4 = dimEnergy;
  dimAEff = {-3, 0, 1};
  dimBEff = {-1, -2, 1};
  dimMH = {-1, 0, 1};
  dimKH = {1, -2, 1};
  dimCHU = {0, -2, 1};
  dimCE = dimVelocity;
  bulkDensity = dimEnergy - 4 dimL;
  braneDensity = dimEnergy - 3 dimL;
  dimD = {0, -4, 2};

  heading["ledger_stage030_electric_scalar_localized_h_closure Mathematica audit"];

  ok = Catch[
    If[activeMutation =!= "" && ! KeyExistsQ[ablations, activeMutation],
      failCount++;
      Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
      Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
      raise["UNKNOWN_MUTATION"]
    ];

    If[activeMutation =!= "",
      Print["ACTIVE_MUTATION=", activeMutation];
      Print["MUTATED_CARD_PRIMITIVE=",
        Lookup[ablations[activeMutation], "Description"]]
    ];

    (* Native matrix route comes first, unlike the SymPy spectral-first route. *)
    subheading[
      "Coupled (u_L,h) kernel via Wolfram positive-definiteness and spectrum"];
    bEffStar = 2;
    kHStar = 1;
    stabilityChuStar = primitive[
      "PASS_STABILITY", "stability_Chu_star", 1/2];
    stiffnessStar = {
      {bEffStar, stabilityChuStar}, {stabilityChuStar, kHStar}};
    dStar = FullSimplify[Det[stiffnessStar]];
    stiffnessEigenvalues = FullSimplify[Eigenvalues[stiffnessStar]];
    stabilityCondition =
      TrueQ[dStar === 7/4] && positiveExactQ[dStar] &&
      TrueQ[PositiveDefiniteMatrixQ[stiffnessStar]] &&
      And @@ (positiveExactQ /@ stiffnessEigenvalues) &&
      TrueQ[dimResidual[dimBEff + dimKH, dimD] === 0] &&
      TrueQ[dimResidual[2 dimCHU, dimD] === 0];
    expectBool[
      "PASS_STABILITY",
      stabilityCondition,
      <|"D*=B_eff* K_h*-C_hu*^2" -> dStar,
        "physical [D]_[L,T,M]" -> dimD,
        "PositiveDefiniteMatrixQ[K*]" ->
          PositiveDefiniteMatrixQ[stiffnessStar],
        "Eigenvalues[K*]" -> stiffnessEigenvalues|>
    ];

    aEffStar = primitive[
      "PASS_POSITIVE_GENERALIZED_WAVE_SPEEDS", "wave_Aeff_star", 1];
    mHStar = 1;
    inertiaStar = DiagonalMatrix[{aEffStar, mHStar}];
    waveStiffnessStar = {{2, 1/2}, {1/2, 1}};
    directDeterminant = Collect[
      Expand[Det[waveStiffnessStar - zStar inertiaStar]], zStar];
    dynamicMatrix = LinearSolve[inertiaStar, waveStiffnessStar];
    characteristicPolynomial = Collect[
      CharacteristicPolynomial[dynamicMatrix, zStar], zStar];
    waveRoots = Sort[FullSimplify[Eigenvalues[dynamicMatrix]]];
    expectedRoots = Sort[{(3 - Sqrt[2])/2, (3 + Sqrt[2])/2}];
    unchangedMargin = Det[waveStiffnessStar];
    waveCondition =
      TrueQ[inertiaStar === IdentityMatrix[2]] &&
      TrueQ[directDeterminant === zStar^2 - 3 zStar + 7/4] &&
      TrueQ[characteristicPolynomial === zStar^2 - 3 zStar + 7/4] &&
      And @@ MapThread[
        TrueQ[FullSimplify[#1 == #2]] &,
        {waveRoots, expectedRoots}] &&
      And @@ (positiveExactQ /@ waveRoots) &&
      TrueQ[PositiveDefiniteMatrixQ[inertiaStar]] &&
      TrueQ[unchangedMargin === 7/4] && positiveExactQ[unchangedMargin];
    expectBool[
      "PASS_POSITIVE_GENERALIZED_WAVE_SPEEDS",
      waveCondition,
      <|"M*" -> inertiaStar,
        "det(K*-z* M*)" -> directDeterminant,
        "CharacteristicPolynomial[M*^-1 K*]" -> characteristicPolynomial,
        "Eigenvalues[M*^-1 K*]=c_pm*^2" -> waveRoots,
        "unchanged D*" -> unchangedMargin|>
    ];
    Print[
      "      native kernel route: PositiveDefiniteMatrixQ / CharacteristicPolynomial / Eigenvalues"];

    (* N0 and all reductions are integrated before the transverse ODE route. *)
    subheading["Integrated zero-mode projection and reduction relations"];
    f0 = Sech[w/ell]^2/ell;
    canonicalN0 = FullSimplify[
      Integrate[2 f0^2, {w, -Infinity, Infinity},
        Assumptions -> ell > 0, GenerateConditions -> False],
      Assumptions -> ell > 0];

    normWeight = primitive["PASS_ZERO_MODE_NORM", "norm_inner_weight", 2];
    normN0 = FullSimplify[
      Integrate[normWeight f0^2, {w, -Infinity, Infinity},
        Assumptions -> ell > 0, GenerateConditions -> False],
      Assumptions -> ell > 0];
    expectBool[
      "PASS_ZERO_MODE_NORM",
      TrueQ[FullSimplify[normN0 == 8/(3 ell), ell > 0]] &&
        TrueQ[FullSimplify[normN0 > 0, ell > 0]],
      <|"N0 via Integrate" -> normN0, "target" -> 8/(3 ell)|>
    ];

    projectionWeight = primitive[
      "PASS_NORMALIZED_ZERO_MODE_PROJECTION",
      "projection_inner_weight", 2];
    p0f0 = FullSimplify[
      Integrate[projectionWeight f0 f0, {w, -Infinity, Infinity},
        Assumptions -> ell > 0, GenerateConditions -> False]/canonicalN0,
      Assumptions -> ell > 0];
    projectedH = FullSimplify[
      Integrate[projectionWeight f0 (h0 f0), {w, -Infinity, Infinity},
        Assumptions -> ell > 0, GenerateConditions -> False]/canonicalN0,
      Assumptions -> ell > 0];
    idempotenceResidual = FullSimplify[p0f0^2 - p0f0];
    expectBool[
      "PASS_NORMALIZED_ZERO_MODE_PROJECTION",
      TrueQ[p0f0 === 1] && TrueQ[projectedH === h0] &&
        TrueQ[idempotenceResidual === 0],
      <|"h=P0 H for H=h0 f0" -> projectedH,
        "P0 f0" -> p0f0, "P0^2 f0-P0 f0" -> idempotenceResidual|>
    ];

    ellStar = 1/20;
    n0Star = FullSimplify[canonicalN0 /. ell -> ellStar];
    baseM4Star = 3/160;
    physicalM4Star = primitive[
      "PASS_PHYSICAL_H_NORM", "physical_M4_star", baseM4Star];
    physicalNormStar = FullSimplify[physicalM4Star n0Star];
    expectBool[
      "PASS_PHYSICAL_H_NORM",
      positiveExactQ[physicalNormStar] &&
        TrueQ[dimResidual[dimM4 + dimN0, dimMH] === 0],
      <|"M4* N0*" -> physicalNormStar,
        "[M4 N0]_[L,T,M]" -> dimM4 + dimN0|>
    ];

    inertiaM4Star = primitive[
      "PASS_REDUCED_H_INERTIA", "inertia_M4_star", baseM4Star];
    reducedMassStar = FullSimplify[n0Star inertiaM4Star];
    expectBool[
      "PASS_REDUCED_H_INERTIA",
      TrueQ[reducedMassStar === 1] && positiveExactQ[reducedMassStar] &&
        TrueQ[dimResidual[dimN0 + dimM4, dimMH] === 0],
      <|"N0*" -> n0Star, "M4*" -> inertiaM4Star,
        "N0* M4*=M_h*" -> reducedMassStar|>
    ];

    cEStar = 1;
    k4DerivedStar = FullSimplify[baseM4Star cEStar^2];
    parentK4Dimension = primitive[
      "PASS_PARENT_STIFFNESS_DIMENSIONAL_CONSISTENCY",
      "parent_K4_dimension", dimK4];
    expectBool[
      "PASS_PARENT_STIFFNESS_DIMENSIONAL_CONSISTENCY",
      TrueQ[dimResidual[
        dimM4 + 2 dimCE, parentK4Dimension] === 0],
      <|"[M4 c_E^2]_[L,T,M]" -> dimM4 + 2 dimCE,
        "declared [K4]_[L,T,M]" -> parentK4Dimension,
        "value dependency" -> "K4=M4 c_E^2 is definitional"|>
    ];

    speedKhN0Power = primitive[
      "PASS_REDUCED_SPEED_PRESERVATION", "speed_Kh_N0_power", 1];
    khProjectedStar = FullSimplify[n0Star^speedKhN0Power k4DerivedStar];
    parentSpeedSquaredStar = FullSimplify[k4DerivedStar/baseM4Star];
    reducedSpeedSquaredStar = FullSimplify[
      khProjectedStar/reducedMassStar];
    expectBool[
      "PASS_REDUCED_SPEED_PRESERVATION",
      TrueQ[FullSimplify[parentSpeedSquaredStar - cEStar^2] === 0] &&
        TrueQ[FullSimplify[
          reducedSpeedSquaredStar - parentSpeedSquaredStar] === 0] &&
        TrueQ[dimResidual[dimMH + 2 dimCE, dimKH] === 0],
      <|"c_E*^2" -> cEStar^2,
        "K4*/M4*" -> parentSpeedSquaredStar,
        "K_h*/M_h*" -> reducedSpeedSquaredStar,
        "M_h*=N0* M4*" -> reducedMassStar,
        "K_h*=N0*^p K4*" -> khProjectedStar,
        "p" -> speedKhN0Power|>
    ];

    khDeclaredStar = primitive[
      "PASS_REDUCED_H_STIFFNESS", "declared_Kh_star", 1];
    expectBool[
      "PASS_REDUCED_H_STIFFNESS",
      TrueQ[khProjectedStar === khDeclaredStar] &&
        positiveExactQ[khDeclaredStar] &&
        TrueQ[dimResidual[dimN0 + dimK4, dimKH] === 0],
      <|"N0* K4*" -> khProjectedStar, "K_h*" -> khDeclaredStar|>
    ];

    Print[
      "      Integrate route: N0*=160/3, M4*=K4*=3/160, M_h*=K_h*=c_E*=1"];

    (* Native differential-equation and asymptotic route. *)
    subheading["Localized parent H via operator action, DSolve, and Limit"];
    parentPotential = (4 - 6 Sech[w/ell]^2)/ell^2;
    factorCoefficient = primitive[
      "PASS_TRANSVERSE_FACTORIZATION", "factor_A_coefficient", 2];
    factorSuperpotential = factorCoefficient Tanh[w/ell]/ell;
    aOnTest = D[testField[w], w] +
      factorSuperpotential testField[w];
    adagAOnTest = -D[aOnTest, w] + factorSuperpotential aOnTest;
    operatorOnTest = -D[testField[w], {w, 2}] +
      parentPotential testField[w];
    factorResidual = FullSimplify[
      adagAOnTest - operatorOnTest,
      Assumptions -> ell > 0 && Element[w, Reals]];
    expectZero[
      "PASS_TRANSVERSE_FACTORIZATION",
      factorResidual,
      <|"O_perp-A^dagger A on testField" -> factorResidual,
        "nonnegative form" -> "<A psi,A psi> >= 0"|>
    ];
    Print[
      "      O_perp=A^dagger A>=0 with A=d/dw+2 tanh(w/ell)/ell"];

    f0Power = primitive[
      "PASS_PARENT_TRANSVERSE_ZERO_EIGENVALUE", "f0_power", 2];
    zeroModeTest = Sech[w/ell]^f0Power/ell;
    zeroModeResidual = FullSimplify[
      -D[zeroModeTest, {w, 2}] + parentPotential zeroModeTest,
      Assumptions -> ell > 0 && Element[w, Reals]];
    expectZero[
      "PASS_PARENT_TRANSVERSE_ZERO_EIGENVALUE",
      zeroModeResidual,
      <|"f0_power" -> f0Power, "O_perp f0" -> zeroModeResidual|>
    ];

    kernelDerivativeCoefficient = primitive[
      "PASS_UNIQUE_NORMALIZABLE_KERNEL",
      "kernel_derivative_coefficient", 1];
    kernelEquation =
      kernelDerivativeCoefficient Derivative[1][testField][w] +
        2 Tanh[w/ell] testField[w]/ell == 0;
    kernelRules = DSolve[kernelEquation, testField, w];
    kernelShape = testField[w] /. First[kernelRules];
    kernelConstants = DeleteDuplicates[Cases[kernelRules, C[_], Infinity]];
    normalizedKernelShape = kernelShape /. C[1] -> 1;
    kernelRatio = FullSimplify[
      kernelShape/f0,
      Assumptions -> ell > 0 && Element[w, Reals]];
    kernelResidual = FullSimplify[
      kernelDerivativeCoefficient D[kernelShape, w] +
        2 Tanh[w/ell] kernelShape/ell,
      Assumptions -> ell > 0 && Element[w, Reals]];
    kernelNorm = FullSimplify[
      Integrate[normalizedKernelShape^2, {w, -Infinity, Infinity},
        Assumptions -> ell > 0, GenerateConditions -> False],
      Assumptions -> ell > 0];
    expectBool[
      "PASS_UNIQUE_NORMALIZABLE_KERNEL",
      TrueQ[Length[kernelRules] === 1] &&
        TrueQ[Length[kernelConstants] === 1] &&
        TrueQ[kernelResidual === 0] && FreeQ[kernelRatio, w] &&
        TrueQ[FullSimplify[kernelNorm > 0, ell > 0]],
      <|"DSolve rules" -> kernelRules,
        "integration constants" -> kernelConstants,
        "DSolve basis/f0" -> kernelRatio,
        "unit-basis norm" -> kernelNorm|>
    ];
    Print[
      "      DSolve gives one integration constant and ker(A)=span{f0}"];

    asymptoteCoefficient = primitive[
      "PASS_POSITIVE_CONTINUUM_GAP", "V_asymptote_coefficient", 4];
    asymptoticPotential =
      (asymptoteCoefficient - 6 Sech[w/ell]^2)/ell^2;
    gapPlus = FullSimplify[
      Limit[asymptoticPotential, w -> Infinity], Assumptions -> ell > 0];
    gapMinus = FullSimplify[
      Limit[asymptoticPotential, w -> -Infinity], Assumptions -> ell > 0];
    gap = 4/ell^2;
    expectBool[
      "PASS_POSITIVE_CONTINUUM_GAP",
      TrueQ[FullSimplify[gapPlus == gap, ell > 0]] &&
        TrueQ[FullSimplify[gapMinus == gap, ell > 0]] &&
        TrueQ[FullSimplify[gap > 0, ell > 0]],
      <|"Limit[V_H,w->+Infinity]" -> gapPlus,
        "Limit[V_H,w->-Infinity]" -> gapMinus, "gap" -> gap|>
    ];
    Print[
      "      unique normalizable zero mode lies below continuum threshold 4/ell^2>0"];

    subheading["Reduced masslessness and action-Hessian symmetry"];
    reducedHMassStar = primitive[
      "PASS_REDUCED_H_MASSLESSNESS", "reduced_h_mass_star", 0];
    qScalar = {
      {omega^2 - 2 k^2, -(1/2) k^2},
      {-(1/2) k^2, omega^2 - k^2 - reducedHMassStar}
    };
    qAtOrigin = qScalar /. {omega -> 0, k -> 0};
    expectBool[
      "PASS_REDUCED_H_MASSLESSNESS",
      TrueQ[reducedHMassStar === 0] && matrixZeroQ[qAtOrigin],
      <|"k^0 h^2 coefficient*" -> reducedHMassStar,
        "Q_s(0,0)" -> qAtOrigin|>
    ];

    conservativeDensity = gradU^2 + gradH^2/2 + gradU gradH/2;
    actionHessian = FullSimplify[
      D[conservativeDensity, {{gradU, gradH}, 2}]];
    uCrossStar = primitive[
      "PASS_CONSERVATIVE_HESSIAN_SYMMETRY",
      "u_equation_cross_star", 1/2];
    eulerOperator = {{2, uCrossStar}, {1/2, 1}};
    expectBool[
      "PASS_CONSERVATIVE_HESSIAN_SYMMETRY",
      TrueQ[actionHessian === Transpose[actionHessian]] &&
        TrueQ[eulerOperator === Transpose[eulerOperator]] &&
        TrueQ[actionHessian === eulerOperator],
      <|"action Hessian" -> actionHessian,
        "displayed Euler operator" -> eulerOperator|>
    ];

    subheading[
      "Units-restored electric-scalar dimensional firewall [L,T,M]"];
    firewallCHU = primitive[
      "PASS_DIMENSIONAL_HOMOGENEITY",
      "firewall_Chu_dimension", dimCHU];
    dimensionTerms = <|
      "H_kin" -> dimM4 + 2 (dimH - dimT),
      "H_grad" -> dimK4 + 2 (dimH - dimL),
      "H_pot" -> dimK4 + dimH + (dimH - 2 dimL),
      "u_kin" -> dimAEff + 2 (dimU - dimT),
      "u_stiff" -> dimBEff + 2 (dimU - dimL),
      "h_kin" -> dimMH + 2 (dimh - dimT),
      "h_stiff" -> dimKH + 2 (dimh - dimL),
      "mix" -> firewallCHU + (dimU - dimL) + (dimh - dimL)
    |>;
    dimensionTargets = <|
      "H_kin" -> bulkDensity, "H_grad" -> bulkDensity,
      "H_pot" -> bulkDensity, "u_kin" -> braneDensity,
      "u_stiff" -> braneDensity, "h_kin" -> braneDensity,
      "h_stiff" -> braneDensity, "mix" -> braneDensity
    |>;
    dimensionResiduals = AssociationMap[
      dimResidual[dimensionTerms[#], dimensionTargets[#]] &,
      Keys[dimensionTerms]];
    expectZero[
      "PASS_DIMENSIONAL_HOMOGENEITY",
      Total[Values[dimensionResiduals]],
      AssociationMap[
        <|"actual" -> dimensionTerms[#],
          "target" -> dimensionTargets[#]|> &,
        Keys[dimensionTerms]]
    ];
    Print[
      "      checked terms={H_kin,H_grad,H_pot,u_kin,u_stiff,h_kin,h_stiff,mix}"];

    subheading["Scope and verdict"];
    Print[
      "SCOPE: static algebraic prerequisites of the charge electric-scalar sector only."];
    Print["  DEFERRED: mouth SOURCE mechanism -> stage031."];
    Print[
      "  EXCLUDED: assembled one-body/two-body/far-field claims -> Part VII."];
    Print[
      "  OUT-OF-SCOPE / UNDISCHARGED CROSS-SECTOR G0 PREDICATE (Part VII de-dup obligation):"];
    Print[
      "    planar-wall factorization; bulk U/Fourier/phase predicates; drain kernel/controllers."];
    Print["VERDICT_TOKEN: ELECTRIC_SCALAR_CLOSURE_STATIC"];

    If[activeMutation =!= "",
      failCount++;
      Print["FIRST_FAILURE=MUTATION_DID_NOT_FIRE"];
      Print["FAIL  MUTATION_DID_NOT_FIRE: ", activeMutation];
      raise["MUTATION_DID_NOT_FIRE"]
    ];

    True,
    "ledgerStage030Failure",
    Function[{message, tag}, False]
  ];

  Print[""];
  Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
  If[TrueQ[ok],
    Print[
      "OVERALL PASS: Mathematica verified stage030 electric-scalar localized-H static closure"];
    Exit[0],
    Print["OVERALL FAIL: Mathematica stage030 audit did not close"];
    Exit[1]
  ]
]
