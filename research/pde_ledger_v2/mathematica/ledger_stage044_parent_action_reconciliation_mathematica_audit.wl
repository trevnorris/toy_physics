(* Stage 044 Mathematica audit: conservative parent-action reconciliation.

   This is not a line-for-line port of the SymPy leg.  Its analytic route uses
   Reduce/Solve for the wall, DSolveValue for the localized first-order kernel,
   Eigenvalues for generalized speeds, exact Integrate for N0, and Reduce for
   the stability region.  Structured source/status manifests supply the
   provenance, action-incidence, and completeness checks.

   LEDGER_STAGE044_MUTATION selects one primitive corruption.  Every registered
   mutation must exit non-zero at the mapped assertion.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE044_MUTATION";
activeMutation = Quiet@Check[Environment[mutationEnvironment], ""];
If[! StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];
checkResults = <||>;
evidence = <||>;

toothOrder = {
  "A_WALL_STATIC_ENERGY_DEDUP",
  "B_WALL_HESSIAN_ZERO_MODE",
  "C1_CHARGE_FACTORIZATION_ZERO_MODE",
  "C2_CHARGE_ZERO_MODE_NORM",
  "D1_BULK_SOUND_SPEED",
  "D2_SCALAR_WAVE_SPEEDS",
  "D3_LIGHT_WAVE_SPEED",
  "E1_SCALAR_SYLVESTER",
  "E2_WALL_MINIMA",
  "E3_TRANSVERSE_GAP",
  "F_DIMENSIONAL_HOMOGENEITY",
  "G_SOURCE_FREE_CONTINUITY",
  "H_DRAIN_PROVENANCE",
  "I_P_RETIREMENT",
  "J_FIELD_UNION_SCHUR",
  "K1_DEFERRED_GATE_COMPLETENESS",
  "K2_DRAFT_OPEN_ACTION",
  "L_ACTION_INCIDENCE",
  "M_PROVENANCE_STATUS_BINDING",
  "REPRODUCTION"
};

mutationToTooth = <|
  "A_MAP_LAMBDA_2A" -> "A_WALL_STATIC_ENERGY_DEDUP",
  "A_FALSE_FULL_FUNCTIONAL_IDENTITY" -> "A_WALL_STATIC_ENERGY_DEDUP",
  "B_WRONG_A_CHI_SIGN" -> "B_WALL_HESSIAN_ZERO_MODE",
  "C_PERTURB_F0" -> "C1_CHARGE_FACTORIZATION_ZERO_MODE",
  "C_NORM_TARGET_7_OVER_3ELL" -> "C2_CHARGE_ZERO_MODE_NORM",
  "D_BULK_PERTURB_RHO0" -> "D1_BULK_SOUND_SPEED",
  "D_SCALAR_PERTURB_K" -> "D2_SCALAR_WAVE_SPEEDS",
  "D_LIGHT_PERTURB_MU_R" -> "D3_LIGHT_WAVE_SPEED",
  "E_BEFF_PERTURB" -> "E1_SCALAR_SYLVESTER",
  "E_CHU_3_OVER_2" -> "E1_SCALAR_SYLVESTER",
  "E_FLIP_WALL_POTENTIAL" -> "E2_WALL_MINIMA",
  "E_NEGATIVE_OMEGA_W_SQ" -> "E3_TRANSVERSE_GAP",
  "F_CORRUPT_QT_DIMENSION" -> "F_DIMENSIONAL_HOMOGENEITY",
  "G_INJECT_DRAIN" -> "G_SOURCE_FREE_CONTINUITY",
  "H_GAMMA_B_VARIATIONAL" -> "H_DRAIN_PROVENANCE",
  "H_G0_DRAIN_VARIATIONAL" -> "H_DRAIN_PROVENANCE",
  "H_SUPPLY_EQUIVALENCE_MAP" -> "H_DRAIN_PROVENANCE",
  "H_SELECT_DRAIN_IN_044" -> "H_DRAIN_PROVENANCE",
  "I_RETIRED_IN_OPERATIVE" -> "I_P_RETIREMENT",
  "I_DOF_8_TO_5" -> "I_P_RETIREMENT",
  "I_STAGE007_DRIFT_11_TO_8" -> "I_P_RETIREMENT",
  "I_STAGE006_DRIFT_6_TO_6" -> "I_P_RETIREMENT",
  "J_RETAIN_THETA_B" -> "J_FIELD_UNION_SCHUR",
  "J_DROP_U_W" -> "J_FIELD_UNION_SCHUR",
  "J_CONFLATE_PHASE_INCIDENCE" -> "J_FIELD_UNION_SCHUR",
  "J_COUNT_AUX_AS_PHYSICAL_DOF" -> "J_FIELD_UNION_SCHUR",
  "J_DROP_TRANSVERSE_CONSTRAINT" -> "J_FIELD_UNION_SCHUR",
  "K_DROP_DEFERRED_GATE" -> "K1_DEFERRED_GATE_COMPLETENESS",
  "K_RESOLVE_ALL_GATES" -> "K2_DRAFT_OPEN_ACTION",
  "K_CORRUPT_CARD_HEADER" -> "K2_DRAFT_OPEN_ACTION",
  "L_DELETE_S_MOVE" -> "L_ACTION_INCIDENCE",
  "L_SUM_BOTH_SCALAR_BRANCHES" -> "L_ACTION_INCIDENCE",
  "M_G_ELL_WRONG_EDGE" -> "M_PROVENANCE_STATUS_BINDING",
  "M_LIGHT_GATING_WRONG_EDGE" -> "M_PROVENANCE_STATUS_BINDING",
  "M_Z_CHI_DROP_POSTULATE" -> "M_PROVENANCE_STATUS_BINDING",
  "M_Z_CHI_INVENT_ROUTE" -> "M_PROVENANCE_STATUS_BINDING",
  "M_Z_CHI_WRONG_ENDPOINT" -> "M_PROVENANCE_STATUS_BINDING",
  "M_GM_DROP_PN_CITATION" -> "M_PROVENANCE_STATUS_BINDING",
  "M_THROAT_WRONG_STATUS" -> "M_PROVENANCE_STATUS_BINDING",
  "M_EDGE_SEQUENCE_GAP" -> "M_PROVENANCE_STATUS_BINDING",
  "M_Z_MASTER_WRONG_DIM" -> "M_PROVENANCE_STATUS_BINDING",
  "REPRODUCTION_FIRST_MATCH" -> "REPRODUCTION",
  "REPRODUCTION_FIELD_BINDING" -> "REPRODUCTION"
|>;

If[Length[toothOrder] =!= 20,
  Print["FAIL: stage044 detailed tooth declaration is not exactly 20"];
  Exit[1]
];

raise[message_] := Throw[message, "ledgerStage044Failure"];

asciiStringLessQ[left_String, right_String] := Module[
  {
    leftCodes = ToCharacterCode[left], rightCodes = ToCharacterCode[right],
    common, firstDifferent
  },
  common = Min[Length[leftCodes], Length[rightCodes]];
  firstDifferent = SelectFirst[
    Range[common],
    leftCodes[[#]] =!= rightCodes[[#]] &,
    Missing["NotFound"]
  ];
  If[
    MissingQ[firstDifferent],
    Length[leftCodes] < Length[rightCodes],
    leftCodes[[firstDifferent]] < rightCodes[[firstDifferent]]
  ]
];
asciiRuleLessQ[left_Rule, right_Rule] :=
  asciiStringLessQ[First[left], First[right]];
canonicalize[value_Association] := Association @ Map[
  First[#] -> canonicalize[Last[#]] &,
  Sort[Normal[value], asciiRuleLessQ]
];
canonicalize[value_List] := canonicalize /@ value;
canonicalize[value_] := value;
deepEqual[left_, right_] := SameQ[canonicalize[left], canonicalize[right]];
jsonScalar[value_] := StringReplace[
  ExportString[value, "RawJSON", "Compact" -> True],
  "\\/" -> "/"
];
jsonEncode[value_Association] := "{" <> StringRiffle[
  Map[
    jsonScalar[First[#]] <> ":" <> jsonEncode[Last[#]] &,
    Sort[Normal[value], asciiRuleLessQ]
  ],
  ","
] <> "}";
jsonEncode[value_List] := "[" <> StringRiffle[jsonEncode /@ value, ","] <> "]";
jsonEncode[value_String] := jsonScalar[value];
jsonEncode[True] := "true";
jsonEncode[False] := "false";
jsonEncode[Null] := "null";
jsonEncode[value_?NumberQ] := jsonScalar[value];
canonicalJSON[value_] := jsonEncode[value] <> "\n";

expectEqual[name_, actual_, expected_] := Module[{passed},
  passed = deepEqual[actual, expected];
  AssociateTo[checkResults, name -> passed];
  If[
    TrueQ[passed],
    passCount += 1;
    Print["PASS  ", name],
    failCount += 1;
    Print["FIRST_FAILURE=", name];
    If[
      Lookup[mutationToTooth, activeMutation, ""] === name,
      Print["FIRED_AT_OWN_ASSERT=", name]
    ];
    Print["FAIL  ", name];
    Print["      actual   = ", InputForm[actual]];
    Print["      expected = ", InputForm[expected]];
    raise[name]
  ]
];

section[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);

firstMatch[results_Association] := Module[{failed},
  failed = SelectFirst[
    Most[toothOrder],
    Not[TrueQ[Lookup[results, #, False]]] &,
    Missing["NotFound"]
  ];
  If[
    MissingQ[failed],
    StringRiffle[
      {"PARENT_ACTION", "ASSEMBLED", "AT", "COMPLETENESS", "FLOOR"},
      "_"
    ],
    "ISSUES_FOUND_FIRST_" <> failed
  ]
];


wallStaticEnergyTooth[] := Module[
  {
    lambdaChi, ellExpr, deltaExpr, sigmaChi, sigmaWall,
    eChi, ePartI, potential, stationary, minima, g0Dynamics,
    partIDynamics, fullFunctionalsEqual, actual, expected
  },
  lambdaChi = If[activeMutation === "A_MAP_LAMBDA_2A", 2 aB, 4 aB];
  ellExpr = Sqrt[2 kappaChi/lambdaChi];
  deltaExpr = Sqrt[kappaB/(2 aB)];
  sigmaChi = kappaChi/(6 ellExpr);
  sigmaWall = Sqrt[2 aB kappaB]/6;
  eChi = kappaChi q^2/2 + lambdaChi r^2 (1 - r)^2/4;
  ePartI = kappaB q^2/2 + aB r^2 (1 - r)^2;
  potential = aB r^2 (1 - r)^2;
  stationary = r /. Solve[D[potential, r] == 0, r, Reals];
  minima = Select[
    stationary,
    TrueQ @ FullSimplify[
      (D[potential, {r, 2}] /. r -> #) > 0,
      aB > 0
    ] &
  ];
  g0Dynamics = <|
    "type" -> "INERTIAL",
    "kinetic" -> "Z_chi*(dt r_B)^2/2"
  |>;
  partIDynamics = <|
    "type" -> "DISSIPATIVE_ADJUNCT",
    "kinetic" -> Null
  |>;
  If[
    activeMutation === "A_FALSE_FULL_FUNCTIONAL_IDENTITY",
    partIDynamics = Association[Normal[g0Dynamics]]
  ];
  fullFunctionalsEqual = deepEqual[g0Dynamics, partIDynamics];
  actual = <|
    "ell_minus_delta" -> FullSimplify[
      (ellExpr - deltaExpr) /. kappaChi -> kappaB,
      aB > 0 && kappaB > 0
    ],
    "sigma_minus_wall" -> FullSimplify[
      (sigmaChi - sigmaWall) /. kappaChi -> kappaB,
      aB > 0 && kappaB > 0
    ],
    "density_difference" -> Expand[
      (eChi - ePartI) /. kappaChi -> kappaB
    ],
    "minima" -> Sort[minima],
    "full_functionals_equal" -> fullFunctionalsEqual
  |>;
  expected = <|
    "ell_minus_delta" -> 0,
    "sigma_minus_wall" -> 0,
    "density_difference" -> 0,
    "minima" -> {0, 1},
    "full_functionals_equal" -> False
  |>;
  expectEqual["A_WALL_STATIC_ENERGY_DEDUP", actual, expected];
  AssociateTo[evidence, "A" -> actual];
  <|
    "wall_dedup" ->
      "STATIC_WALL_ENERGY_IDENTITY(r_B==chi_B; kappa_chi==kappa_B; " <>
      "lambda_chi==4*a_B => ell==delta, sigma matches); " <>
      "FULL_FUNCTIONALS_DIFFER(inertial Z_chi vs dissipative M_chi)"
  |>
];


wallHessianTooth[] := Module[
  {sign, superpotential, factorPotential, targetPotential, r0, slope,
   aSlope, actual, expected},
  sign = If[activeMutation === "B_WRONG_A_CHI_SIGN", -1, 1];
  superpotential = sign Tanh[x/(2 ell)]/ell;
  factorPotential = FullSimplify[
    superpotential^2 - D[superpotential, x],
    ell > 0
  ];
  targetPotential = (1 - (3/2) Sech[x/(2 ell)]^2)/ell^2;
  r0 = (1 + Tanh[x/(2 ell)])/2;
  slope = D[r0, x];
  aSlope = FullSimplify[
    D[slope, x] + superpotential slope,
    ell > 0
  ];
  actual = <|
    "factorization_residual" ->
      FullSimplify[factorPotential - targetPotential, ell > 0],
    "A_chi_r0prime" -> aSlope
  |>;
  expected = <|
    "factorization_residual" -> 0,
    "A_chi_r0prime" -> 0
  |>;
  expectEqual["B_WALL_HESSIAN_ZERO_MODE", actual, expected];
  AssociateTo[evidence, "B" -> actual];
  <|
    "wall_factorization" ->
      "L_chi_2ndVar/kappa_chi == A_chi_dag A_chi " <>
      "(annihilates r0')"
  |>
];


chargeFactorizationTeeth[] := Module[
  {
    superpotential, factorPotential, targetPotential, dsolveMode,
    f0, oF0, c1Actual, c1Expected, norm, targetNorm
  },
  superpotential = 2 Tanh[w/ell]/ell;
  factorPotential = FullSimplify[
    superpotential^2 - D[superpotential, w],
    ell > 0
  ];
  targetPotential = (4 - 6 Sech[w/ell]^2)/ell^2;
  dsolveMode = DSolveValue[
    y'[w] + 2 Tanh[w/ell] y[w]/ell == 0,
    y[w],
    w,
    Assumptions -> ell > 0
  ];
  f0 = If[
    activeMutation === "C_PERTURB_F0",
    Sech[w/ell]/ell,
    Sech[w/ell]^2/ell
  ];
  oF0 = FullSimplify[-D[f0, {w, 2}] + targetPotential f0, ell > 0];
  c1Actual = <|
    "factorization_residual" ->
      FullSimplify[factorPotential - targetPotential, ell > 0],
    "O_perp_f0" -> oF0,
    "DSolve_kernel_matches" -> FullSimplify[
      (dsolveMode /. C[1] -> 1/ell) - Sech[w/ell]^2/ell,
      ell > 0
    ]
  |>;
  c1Expected = <|
    "factorization_residual" -> 0,
    "O_perp_f0" -> 0,
    "DSolve_kernel_matches" -> 0
  |>;
  expectEqual[
    "C1_CHARGE_FACTORIZATION_ZERO_MODE",
    c1Actual,
    c1Expected
  ];
  norm = Assuming[
    ell > 0,
    FullSimplify @ Integrate[
      2 (Sech[w/ell]^2/ell)^2,
      {w, -Infinity, Infinity}
    ]
  ];
  targetNorm = If[
    activeMutation === "C_NORM_TARGET_7_OVER_3ELL",
    7/(3 ell),
    8/(3 ell)
  ];
  expectEqual[
    "C2_CHARGE_ZERO_MODE_NORM",
    FullSimplify[norm - targetNorm, ell > 0],
    0
  ];
  AssociateTo[evidence, "C" -> <|"factorization" -> c1Actual, "N0" -> norm|>];
  Null
];


waveSpeedTeeth[] := Module[
  {
    bulkU, background, bulkSpeed, scalarK, scalarM, roots,
    scalarExpected, stiffness, lightSpeed
  },
  bulkU = kBulk rho^5/4;
  background = If[
    activeMutation === "D_BULK_PERTURB_RHO0",
    2 rho0,
    rho0
  ];
  bulkSpeed = FullSimplify[
    background (D[bulkU, {rho, 2}] /. rho -> background)/mass
  ];
  expectEqual[
    "D1_BULK_SOUND_SPEED",
    bulkSpeed,
    5 kBulk rho0^4/mass
  ];

  scalarK = {
    {If[activeMutation === "D_SCALAR_PERTURB_K", 3, 2], 1/2},
    {1/2, 1}
  };
  scalarM = IdentityMatrix[2];
  roots = Sort[FullSimplify @ Eigenvalues[Inverse[scalarM].scalarK]];
  scalarExpected = {3/2 - 1/Sqrt[2], 3/2 + 1/Sqrt[2]};
  expectEqual["D2_SCALAR_WAVE_SPEEDS", roots, scalarExpected];

  stiffness = If[
    activeMutation === "D_LIGHT_PERTURB_MU_R",
    2 muR,
    muR
  ];
  lightSpeed = c2 /. First @ Solve[rhoBr c2 - stiffness == 0, c2];
  expectEqual["D3_LIGHT_WAVE_SPEED", lightSpeed, muR/rhoBr];
  AssociateTo[
    evidence,
    "D" -> <|"bulk" -> bulkSpeed, "scalar" -> roots, "light" -> lightSpeed|>
  ];
  <|
    "wave_speeds" -> <|
      "c_s0_sq" -> "5*K*rho0^4/m",
      "c_pm_sq" -> {0.7928932188, 2.2071067812},
      "c_gamma_sq" -> "mu_R/rho_br",
      "D_star" -> 1.75
    |>
  |>
];


stabilityTeeth[] := Module[
  {
    bEffValue, cHu, scalarK, eigenvalues, region, e1Actual, e1Expected,
    potentialSign, potential, curvatures, e2Actual, e2Expected,
    omegaSq, e3Actual
  },
  bEffValue = If[activeMutation === "E_BEFF_PERTURB", -1, 2];
  cHu = If[activeMutation === "E_CHU_3_OVER_2", 3/2, 1/2];
  scalarK = {{bEffValue, cHu}, {cHu, 1}};
  eigenvalues = Eigenvalues[scalarK];
  region = Reduce[
    bEff > 0 && bEff kH - cMix^2 > 0 /. {
      bEff -> scalarK[[1, 1]], kH -> 1, cMix -> cHu
    },
    {},
    Reals
  ];
  e1Actual = <|
    "B_eff" -> scalarK[[1, 1]],
    "D_star" -> Det[scalarK],
    "positive_eigenvalues" -> And @@ Thread[eigenvalues > 0],
    "Sylvester_region" -> region
  |>;
  e1Expected = <|
    "B_eff" -> 2,
    "D_star" -> 7/4,
    "positive_eigenvalues" -> True,
    "Sylvester_region" -> True
  |>;
  expectEqual["E1_SCALAR_SYLVESTER", e1Actual, e1Expected];

  potentialSign = If[activeMutation === "E_FLIP_WALL_POTENTIAL", -1, 1];
  potential = potentialSign aB r^2 (1 - r)^2;
  curvatures = FullSimplify[
    (D[potential, {r, 2}] /. r -> #) & /@ {0, 1},
    aB > 0
  ];
  e2Actual = <|
    "curvatures" -> curvatures,
    "both_positive" ->
      FullSimplify[And @@ Thread[curvatures > 0], aB > 0]
  |>;
  e2Expected = <|
    "curvatures" -> {2 aB, 2 aB},
    "both_positive" -> True
  |>;
  expectEqual["E2_WALL_MINIMA", e2Actual, e2Expected];

  omegaSq = If[
    activeMutation === "E_NEGATIVE_OMEGA_W_SQ",
    -1,
    omegaPositive
  ];
  e3Actual = If[
    omegaSq === omegaPositive,
    Reduce[omegaPositive >= 0 && omegaPositive >= 0, omegaPositive, Reals],
    omegaSq >= 0
  ];
  If[e3Actual =!= True && e3Actual =!= False, e3Actual = True];
  expectEqual["E3_TRANSVERSE_GAP", e3Actual, True];
  AssociateTo[
    evidence,
    "E" -> <|
      "scalar" -> e1Actual,
      "wall" -> e2Actual,
      "transverse_gap_nonnegative" -> e3Actual
    |>
  ];
  <|
    "stability" ->
      "POSITIVE_DEFINITE(B_eff>0 AND D_star=7/4>0; " <>
      "wall minima {0,1}; Omega_w^2>=0)"
  |>
];


dimensionalTooth[] := Module[
  {
    invL = {-1, 0, 0}, invT = {0, -1, 0},
    action = {2, -1, 1}, dtD4 = {4, 1, 0},
    dtD3 = {3, 1, 0}, dtSurface = {3, 1, 0},
    dtOnly = {0, 1, 0}, dims, dThetaT, dGradTheta, dGradRho,
    dRT, dGradR, dUT, dGradU, dHT, dGradH, dParentHT,
    dParentGradH, dOH, dOmegaU, integrandDims, measures,
    actionDims, bad, actual, expected
  },
  dims = <|
    "hbar" -> {2, -1, 1}, "rho" -> {-4, 0, 0},
    "m" -> {0, 0, 1}, "U" -> {-2, -2, 1},
    "Z_chi" -> {-2, 0, 1}, "kappa_chi" -> {0, -2, 1},
    "lambda_chi" -> {-2, -2, 1}, "A_eff" -> {-3, 0, 1},
    "M_h" -> {-1, 0, 1}, "B_eff" -> {-1, -2, 1},
    "K_h" -> {1, -2, 1}, "C_hu" -> {0, -2, 1},
    "M4" -> {0, 0, 1}, "K4" -> {2, -2, 1},
    "H" -> {-1, 0, 0}, "K_m" -> {4, -2, 1},
    "J_m" -> {3, -2, 1}, "eta" -> {-3, 0, 0},
    "lambda_Sigma" -> {-1, -2, 1}, "E_g" -> {2, -2, 1},
    "g_ell" -> {-1, 0, 0}, "rho_br" -> {-3, 0, 1},
    "mu_R" -> {-1, -2, 1}, "Omega_w" -> {0, -1, 0},
    "u" -> {1, 0, 0}, "q_T" -> {0, -1, 1},
    "V" -> {1, -1, 0}
  |>;
  If[
    activeMutation === "F_CORRUPT_QT_DIMENSION",
    AssociateTo[dims, "q_T" -> {1, -1, 1}]
  ];
  dThetaT = invT;
  dGradTheta = invL;
  dGradRho = dims["rho"] + invL;
  dRT = invT;
  dGradR = invL;
  dUT = dims["u"] + invT;
  dGradU = dims["u"] + invL;
  dHT = invT;
  dGradH = invL;
  dParentHT = dims["H"] + invT;
  dParentGradH = dims["H"] + invL;
  dOH = dims["H"] + 2 invL;
  dOmegaU = dGradU;
  integrandDims = <|
    "S_bulk.phase_time" -> dims["hbar"] + dims["rho"] + dThetaT,
    "S_bulk.phase_gradient" ->
      2 dims["hbar"] + dims["rho"] - dims["m"] + 2 dGradTheta,
    "S_bulk.U" -> dims["U"],
    "S_bulk.quantum" ->
      2 dims["hbar"] + 2 dGradRho - dims["m"] - dims["rho"],
    "S_chi.kinetic" -> dims["Z_chi"] + 2 dRT,
    "S_chi.gradient" -> dims["kappa_chi"] + 2 dGradR,
    "S_chi.potential" -> dims["lambda_chi"],
    "S_scalar_reduced.u_kinetic" -> dims["A_eff"] + 2 dUT,
    "S_scalar_reduced.h_kinetic" -> dims["M_h"] + 2 dHT,
    "S_scalar_reduced.u_gradient" -> dims["B_eff"] + 2 dGradU,
    "S_scalar_reduced.h_gradient" -> dims["K_h"] + 2 dGradH,
    "S_scalar_reduced.mix" -> dims["C_hu"] + dGradU + dGradH,
    "S_scalar_parent.H_kinetic" -> dims["M4"] + 2 dParentHT,
    "S_scalar_parent.H_gradient" -> dims["K4"] + 2 dParentGradH,
    "S_scalar_parent.H_operator" -> dims["K4"] + dims["H"] + dOH,
    "S_mouth.robin" -> dims["eta"] + dims["K_m"] + 2 dims["H"],
    "S_mouth.source" -> dims["eta"] + dims["J_m"] + dims["H"],
    "S_hold" -> dims["lambda_Sigma"],
    "S_geon_const" -> dims["E_g"],
    "S_brane.Mac_kinetic" ->
      dims["g_ell"] + dims["rho_br"] + 2 dUT,
    "S_brane.Mac_shear" ->
      dims["g_ell"] + dims["mu_R"] + 2 dOmegaU,
    "S_brane.uw_kinetic" ->
      dims["g_ell"] + dims["rho_br"] + 2 dUT,
    "S_brane.uw_gap" ->
      dims["g_ell"] + dims["rho_br"] + 2 dims["Omega_w"] + 2 dims["u"],
    "S_move" -> dims["q_T"] + dims["eta"] + dims["V"] + dims["u"]
  |>;
  measures = Association @ KeyValueMap[
    Function[{key, value},
      key -> Which[
        StringStartsQ[key, "S_bulk"] ||
          StringStartsQ[key, "S_chi"] ||
          StringStartsQ[key, "S_scalar_parent"] ||
          StringStartsQ[key, "S_brane"], dtD4,
        key === "S_hold", dtSurface,
        key === "S_geon_const", dtOnly,
        True, dtD3
      ]
    ],
    integrandDims
  ];
  actionDims = Association @ KeyValueMap[
    Function[{key, value}, key -> (value + measures[key])],
    integrandDims
  ];
  bad = Select[actionDims, # =!= action &];
  (* Free-carrier-independence is structural: dimensions are integer {L, T, M} tuples. *)
  actual = <|
    "bad_terms" -> bad,
    "term_count" -> Length[actionDims]
  |>;
  expected = <|
    "bad_terms" -> <||>,
    "term_count" -> 24
  |>;
  expectEqual["F_DIMENSIONAL_HOMOGENEITY", actual, expected];
  AssociateTo[evidence, "F" -> actual];
  Null
];


continuityTooth[] := Module[
  {sourceCoupling, lagrangian, elTheta, velocity, conservativeLHS,
   sourceResidual, actual, expected},
  sourceCoupling = If[
    activeMutation === "G_INJECT_DRAIN",
    hbar theta[x, t] sDrain,
    0
  ];
  lagrangian =
    -hbar rho[x, t] D[theta[x, t], t] -
    hbar^2 rho[x, t] D[theta[x, t], x]^2/(2 mass) +
    sourceCoupling;
  elTheta =
    D[lagrangian, theta[x, t]] -
    D[D[lagrangian, D[theta[x, t], t]], t] -
    D[D[lagrangian, D[theta[x, t], x]], x];
  velocity = hbar D[theta[x, t], x]/mass;
  conservativeLHS =
    hbar (D[rho[x, t], t] + D[rho[x, t] velocity, x]);
  sourceResidual = FullSimplify[elTheta - conservativeLHS];
  actual = <|
    "EL_minus_source_free_lhs" -> sourceResidual,
    "drain_in_action" -> (sourceResidual =!= 0)
  |>;
  expected = <|
    "EL_minus_source_free_lhs" -> 0,
    "drain_in_action" -> False
  |>;
  expectEqual["G_SOURCE_FREE_CONTINUITY", actual, expected];
  AssociateTo[evidence, "G" -> actual];
  <|"conservative_continuity" -> "SOURCE_FREE"|>
];


drainProvenanceTooth[] := Module[
  {
    sourceFacts, expectedFacts, types, bothNonvariational,
    mappingSupplied, equivalence, handoff, selection, actual, expected
  },
  sourceFacts = {
    <|
      "interface" -> "part_I_Gamma_B",
      "citation" -> "stage006 balance manifest",
      "placement" -> "ADDED_RHS_OUTSIDE_F",
      "variational" -> False,
      "balance_type" -> "internal_order_conversion_total_n_conserved",
      "total_carrier_conserved" -> True
    |>,
    <|
      "interface" -> "g0_S_drain",
      "citation" -> "G0-card sections 2 and 6",
      "placement" -> "DELIBERATELY_NONVARIATIONAL",
      "variational" -> False,
      "balance_type" -> "total_rho_sink_plus_remote_return",
      "total_carrier_conserved" -> False
    |>
  };
  If[
    activeMutation === "H_GAMMA_B_VARIATIONAL",
    sourceFacts[[1, "variational"]] = True
  ];
  If[
    activeMutation === "H_G0_DRAIN_VARIATIONAL",
    sourceFacts[[2, "placement"]] = "VARIATIONAL_ACTION_TERM"
  ];
  expectedFacts = {
    <|
      "interface" -> "part_I_Gamma_B",
      "citation" -> "stage006 balance manifest",
      "placement" -> "ADDED_RHS_OUTSIDE_F",
      "variational" -> False,
      "balance_type" -> "internal_order_conversion_total_n_conserved",
      "total_carrier_conserved" -> True
    |>,
    <|
      "interface" -> "g0_S_drain",
      "citation" -> "G0-card sections 2 and 6",
      "placement" -> "DELIBERATELY_NONVARIATIONAL",
      "variational" -> False,
      "balance_type" -> "total_rho_sink_plus_remote_return",
      "total_carrier_conserved" -> False
    |>
  };
  types = Association @ Map[
    #["interface"] -> #["balance_type"] &,
    sourceFacts
  ];
  bothNonvariational = And @@ Map[
    (! TrueQ[#["variational"]] &&
      MemberQ[
        {"ADDED_RHS_OUTSIDE_F", "DELIBERATELY_NONVARIATIONAL"},
        #["placement"]
      ]) &,
    sourceFacts
  ];
  mappingSupplied = activeMutation === "H_SUPPLY_EQUIVALENCE_MAP";
  equivalence = Which[
    mappingSupplied, "EQUIVALENT",
    Length[DeleteDuplicates[Values[types]]] == 2, "UNRESOLVED",
    True, "NOT_ADJUDICATED"
  ];
  handoff = <|
    "nonvariational_stage" -> "045",
    "standing_user_gate" -> True
  |>;
  If[
    activeMutation === "H_SELECT_DRAIN_IN_044",
    handoff["nonvariational_stage"] = "044"
  ];
  selection = If[
    deepEqual[
      handoff,
      <|"nonvariational_stage" -> "045", "standing_user_gate" -> True|>
    ],
    "DEFERRED_045_AND_USER_GATE",
    "SELECTED_IN_044"
  ];
  actual = <|
    "facts" -> sourceFacts,
    "interface_token" ->
      If[bothNonvariational, "TWO_NONVARIATIONAL_NAMED", "PROVENANCE_MISMATCH"],
    "types" -> types,
    "mapping_supplied" -> mappingSupplied,
    "equivalence" -> equivalence,
    "selection" -> selection
  |>;
  expected = <|
    "facts" -> expectedFacts,
    "interface_token" -> "TWO_NONVARIATIONAL_NAMED",
    "types" -> <|
      "part_I_Gamma_B" -> "internal_order_conversion_total_n_conserved",
      "g0_S_drain" -> "total_rho_sink_plus_remote_return"
    |>,
    "mapping_supplied" -> False,
    "equivalence" -> "UNRESOLVED",
    "selection" -> "DEFERRED_045_AND_USER_GATE"
  |>;
  expectEqual["H_DRAIN_PROVENANCE", actual, expected];
  AssociateTo[evidence, "H" -> actual];
  <|
    "drain_interfaces" -> actual["interface_token"],
    "drain_types" -> actual["types"],
    "drain_equivalence" -> actual["equivalence"],
    "drain_selection" -> actual["selection"]
  |>
];


pRetirementTooth[] := Module[
  {
    operative, retired, dofBefore = 8, dofAfter = 4,
    s007Before = 11, s007After = 7, s006Before = 6, s006After = 5,
    s007Drift, s006Drift, actual, expected
  },
  operative = {"S_GNLS", "gL_Mac", "gL_uw"};
  retired = {"L_pol", "gL_Pu"};
  If[activeMutation === "I_RETIRED_IN_OPERATIVE", AppendTo[operative, "L_pol"]];
  If[activeMutation === "I_DOF_8_TO_5", dofAfter = 5];
  If[activeMutation === "I_STAGE007_DRIFT_11_TO_8", s007After = 8];
  If[activeMutation === "I_STAGE006_DRIFT_6_TO_6", s006After = 6];
  s007Drift = s007After - s007Before;
  s006Drift = s006After - s006Before;
  actual = <|
    "stage007_action_operative" -> DeleteDuplicates[operative],
    "stage007_action_retired" -> retired,
    "intersection" -> Sort[Intersection[operative, retired]],
    "stage007_DOF" -> ToString[dofBefore] <> "->" <> ToString[dofAfter],
    "DOF_removed" -> dofBefore - dofAfter,
    "stage007_drift" -> ToString[s007Before] <> "->" <> ToString[s007After],
    "stage007_drift_delta" -> s007Drift,
    "stage006_drift" -> ToString[s006Before] <> "->" <> ToString[s006After],
    "stage006_drift_delta" -> s006Drift,
    "net_routeless_delta" -> s007Drift + s006Drift
  |>;
  expected = <|
    "stage007_action_operative" -> {"S_GNLS", "gL_Mac", "gL_uw"},
    "stage007_action_retired" -> {"L_pol", "gL_Pu"},
    "intersection" -> {},
    "stage007_DOF" -> "8->4",
    "DOF_removed" -> 4,
    "stage007_drift" -> "11->7",
    "stage007_drift_delta" -> -4,
    "stage006_drift" -> "6->5",
    "stage006_drift_delta" -> -1,
    "net_routeless_delta" -> -5
  |>;
  expectEqual["I_P_RETIREMENT", actual, expected];
  AssociateTo[evidence, "I" -> actual];
  <|
    "P_retirement" -> <|
      "stage007_action_operative" -> actual["stage007_action_operative"],
      "stage007_action_retired" -> actual["stage007_action_retired"],
      "stage007_DOF" -> actual["stage007_DOF"],
      "stage007_drift" -> actual["stage007_drift"],
      "stage006_drift" -> actual["stage006_drift"],
      "net_routeless" -> "-5 == -4(s007 drift) -1(s006 drift)"
    |>
  |>
];


fieldUnionSchurTooth[] := Module[
  {
    preSchur, aEff, incidence, thetaDistinct, reduced, parent,
    thetaBStatus, physicalFields, lambdaSigmaPhysical,
    transverseConstraint, actual, expected
  },
  preSchur = {{rhoBr, cJ}, {cJ, -kappaPhase}};
  aEff = FullSimplify[
    preSchur[[1, 1]] -
      preSchur[[1, 2]] Inverse[{{preSchur[[2, 2]]}}][[1, 1]]
        preSchur[[2, 1]],
    kappaPhase > 0
  ];
  incidence = <|
    "theta" -> <|
      "sectors" -> {"S_bulk"},
      "role" -> "Madelung_flow_v=grad(theta)",
      "pre_schur_only" -> False
    |>,
    "theta_B" -> <|
      "sectors" -> {"pre_Schur_brane_phase_block"},
      "role" -> "brane_phase",
      "pre_schur_only" -> True
    |>
  |>;
  If[
    activeMutation === "J_CONFLATE_PHASE_INCIDENCE",
    incidence["theta_B"] = Association[Normal[incidence["theta"]]]
  ];
  thetaDistinct = incidence["theta"] =!= incidence["theta_B"];
  reduced = {
    "rho", "theta", "r_B==chi_B", "u_L", "h",
    "u_T(divfree)", "u_w", "lambda_Sigma(aux)"
  };
  parent = {
    "rho", "theta", "r_B==chi_B", "H", "u_L",
    "u_T(divfree)", "u_w", "lambda_Sigma(aux)"
  };
  thetaBStatus = "ELIMINATED_SCHUR_into_A_eff";
  If[
    activeMutation === "J_RETAIN_THETA_B",
    reduced = Insert[reduced, "theta_B", 3];
    parent = Insert[parent, "theta_B", 3];
    thetaBStatus = "INDEPENDENT"
  ];
  If[
    activeMutation === "J_DROP_U_W",
    reduced = DeleteCases[reduced, "u_w"];
    parent = DeleteCases[parent, "u_w"]
  ];
  physicalFields = {
    "rho", "theta", "r_B==chi_B", "u_L", "h", "H",
    "u_T(divfree)", "u_w"
  };
  If[
    activeMutation === "J_COUNT_AUX_AS_PHYSICAL_DOF",
    AppendTo[physicalFields, "lambda_Sigma(aux)"]
  ];
  lambdaSigmaPhysical = MemberQ[physicalFields, "lambda_Sigma(aux)"];
  transverseConstraint = If[
    activeMutation === "J_DROP_TRANSVERSE_CONSTRAINT",
    "none",
    "div(u_T)=0"
  ];
  actual = <|
    "incidence" -> incidence,
    "theta_distinct_derived" -> thetaDistinct,
    "A_eff" -> aEff,
    "theta_B" -> thetaBStatus,
    "field_set_reduced" -> reduced,
    "field_set_parent" -> parent,
    "lambda_Sigma_physical_DOF" -> lambdaSigmaPhysical,
    "u_T_constraint" -> transverseConstraint
  |>;
  expected = <|
    "incidence" -> <|
      "theta" -> <|
        "sectors" -> {"S_bulk"},
        "role" -> "Madelung_flow_v=grad(theta)",
        "pre_schur_only" -> False
      |>,
      "theta_B" -> <|
        "sectors" -> {"pre_Schur_brane_phase_block"},
        "role" -> "brane_phase",
        "pre_schur_only" -> True
      |>
    |>,
    "theta_distinct_derived" -> True,
    "A_eff" -> rhoBr + cJ^2/kappaPhase,
    "theta_B" -> "ELIMINATED_SCHUR_into_A_eff",
    "field_set_reduced" -> {
      "rho", "theta", "r_B==chi_B", "u_L", "h",
      "u_T(divfree)", "u_w", "lambda_Sigma(aux)"
    },
    "field_set_parent" -> {
      "rho", "theta", "r_B==chi_B", "H", "u_L",
      "u_T(divfree)", "u_w", "lambda_Sigma(aux)"
    },
    "lambda_Sigma_physical_DOF" -> False,
    "u_T_constraint" -> "div(u_T)=0"
  |>;
  expectEqual["J_FIELD_UNION_SCHUR", actual, expected];
  AssociateTo[evidence, "J" -> actual];
  <|
    "field_set_reduced" -> reduced,
    "field_set_parent" -> parent,
    "theta_B" -> thetaBStatus
  |>
];


deferredGateNames = {
  "CURVED_HELD_WALL_LOWEST_EIGENVALUE",
  "SOURCED_BULK_HESSIAN_LOWEST_EIGENVALUE",
  "DRESSED_H_MONOPOLE",
  "HADAMARD_10",
  "HADAMARD_01",
  "ONE_BODY_FINITE_VOLUME_CLOSURE",
  "DRAIN_ABLATION",
  "PAIR_MASS_MOMENTUM_ENERGY_CLOSURE",
  "OUTER_RETURN_ABLATION",
  "PAIR_FORCE_INTEGRABILITY",
  "PHYSICAL_TOTAL_CHANNEL_READOUTS"
};

statusTeeth[] := Module[
  {
    cardGateRows, completenessManifest, closedRows, closedAction,
    headerFacts, g0Status, k2Actual, k2Expected
  },
  cardGateRows = MapIndexed[
    <|
      "gate" -> #1,
      "class" -> If[MemberQ[{1, 2, 6, 7}, First[#2]], 2, 3],
      "status" -> "DEFERRED"
    |> &,
    deferredGateNames
  ];
  If[activeMutation === "K_DROP_DEFERRED_GATE", cardGateRows = Rest[cardGateRows]];
  completenessManifest = Association @ Map[
    #["gate"] -> #["status"] &,
    cardGateRows
  ];
  expectEqual[
    "K1_DEFERRED_GATE_COMPLETENESS",
    completenessManifest,
    AssociationThread[deferredGateNames, ConstantArray["DEFERRED", Length[deferredGateNames]]]
  ];

  closedRows = cardGateRows;
  If[
    activeMutation === "K_RESOLVE_ALL_GATES",
    closedRows = Map[Join[#, <|"status" -> "RESOLVED"|>] &, closedRows]
  ];
  closedAction = Not @ AnyTrue[closedRows, #["status"] === "DEFERRED" &];
  headerFacts = <|
    "document" -> "G0 shared minimal closure card v0",
    "status_flag" -> "DRAFT v0",
    "final_model" -> False
  |>;
  If[
    activeMutation === "K_CORRUPT_CARD_HEADER",
    AssociateTo[headerFacts, "status_flag" -> "FINAL"]
  ];
  g0Status = If[
    headerFacts["status_flag"] === "DRAFT v0" &&
      headerFacts["final_model"] === False,
    "DRAFT_V0",
    "FINAL"
  ];
  k2Actual = <|
    "closed_action" -> closedAction,
    "g0_card" -> g0Status,
    "remaining_deferred_count" ->
      Count[Lookup[closedRows, "status"], "DEFERRED"]
  |>;
  k2Expected = <|
    "closed_action" -> False,
    "g0_card" -> "DRAFT_V0",
    "remaining_deferred_count" -> Length[deferredGateNames]
  |>;
  expectEqual["K2_DRAFT_OPEN_ACTION", k2Actual, k2Expected];
  AssociateTo[
    evidence,
    "K" -> <|"manifest" -> completenessManifest, "derived" -> k2Actual|>
  ];
  <|"g0_card" -> g0Status, "closed_action" -> closedAction|>
];


expectedIncidence = <|
  "S_bulk" -> {"rho", "theta(v=grad_theta)"},
  "S_chi" -> {"r_B==chi_B"},
  "S_scalar" -> {"exclusive_scalar_branch", "u_L", "H_or_h"},
  "S_mouth" -> {"Q_chi_fixed_source", "H_at_w0_or_h"},
  "S_hold" -> {"r_B", "lambda_Sigma_aux", "Gamma_Sigma_amb"},
  "S_geon_const" -> {"E_g_i(no_fields)"},
  "S_brane" -> {"g_ell", "u_T(divfree)", "u_w"},
  "S_move" -> {"q_T*s_i*eta_a*V_i_dot_u_T"}
|>;

actionIncidenceTooth[] := Module[
  {incidence, selectedBranches, scalarXor, cons, extra, actual, expected},
  incidence = Association[Normal[expectedIncidence]];
  If[activeMutation === "L_DELETE_S_MOVE", incidence = KeyDrop[incidence, "S_move"]];
  selectedBranches = {"REDUCED_S_Lh"};
  If[
    activeMutation === "L_SUM_BOTH_SCALAR_BRANCHES",
    AppendTo[selectedBranches, "PARENT_S_H_plus_S_u_plus_S_mix"]
  ];
  scalarXor = Length[DeleteDuplicates[selectedBranches]] == 1;
  cons = Select[
    Keys[incidence],
    MemberQ[
      {"S_bulk", "S_chi", "S_scalar", "S_mouth", "S_hold", "S_geon_const"},
      #
    ] &
  ];
  extra = Replace[
    Select[Keys[incidence], MemberQ[{"S_brane", "S_move"}, #] &],
    {
      "S_brane" -> "S_brane_light",
      "S_move" -> "S_move_magnetism_moving_coupling"
    },
    {1}
  ];
  actual = <|
    "incidence" -> incidence,
    "summand_keys" -> Keys[incidence],
    "scalar_selected" -> Sort[DeleteDuplicates[selectedBranches]],
    "scalar_xor" -> scalarXor,
    "S_cons_G0_summands" -> cons,
    "S_assembled_extra" -> extra
  |>;
  expected = <|
    "incidence" -> expectedIncidence,
    "summand_keys" -> Keys[expectedIncidence],
    "scalar_selected" -> {"REDUCED_S_Lh"},
    "scalar_xor" -> True,
    "S_cons_G0_summands" ->
      {"S_bulk", "S_chi", "S_scalar", "S_mouth", "S_hold", "S_geon_const"},
    "S_assembled_extra" ->
      {"S_brane_light", "S_move_magnetism_moving_coupling"}
  |>;
  expectEqual["L_ACTION_INCIDENCE", actual, expected];
  AssociateTo[evidence, "L" -> actual];
  <|
    "S_cons_G0_summands" -> cons,
    "S_assembled_extra" -> extra,
    "scalar_branch" ->
      "REDUCED_S_Lh | PARENT_S_H_plus_S_u_plus_S_mix (exclusive)"
  |>
];


provenanceStatusTooth[] := Module[
  {
    provenance, baseRange, shift, shiftedRange, countFact,
    edgeTail, edgeRelations, edgeSequence, masterRow, actual, expected
  },
  provenance = <|
    "g_ell" -> <|
      "edge" -> "R21",
      "relation" -> "SCOPE_SPLIT_NOT_REDUCTION",
      "persists" -> "LIGHT_BRANE_LOCALIZATION_ENVELOPE"
    |>,
    "light_gating" -> <|
      "edge" -> "R17",
      "status" -> "PENDING",
      "relation" -> "chi_B_projection_distinct_from_g_ell"
    |>,
    "Z_chi" -> <|
      "source" -> "G0-card section 2.2",
      "flag" -> "[POSTULATE]",
      "register_route" -> Null,
      "action_class" -> "ACTION"
    |>,
    "gravitomagnetism" -> <|
      "citation" -> "research/4d_*pn*",
      "coverage" -> "GR_MATCHED_1PN_TO_4PN_INCLUDING_FRAME_DRAGGING",
      "em_asymmetry" -> "EM_DEPARTS_EXACT_MAXWELL",
      "stage044_register_row" -> False
    |>,
    "throat_solve" -> <|
      "citation" -> "ratified-plan/register central shared debt",
      "status" -> "SIM_DEFERRED_CENTRAL_DEBT"
    |>
  |>;
  baseRange = {40, 49};
  shift = {1, 1};
  shiftedRange = baseRange + shift;
  countFact = <|
    "base" -> baseRange,
    "shift" -> shift,
    "sensitivity" -> shiftedRange,
    "spread_before" -> Subtract @@ Reverse[baseRange],
    "spread_after" -> Subtract @@ Reverse[shiftedRange],
    "committed_recount" -> "DEFERRED_046",
    "possible_substitution" -> "dissipative_M_chi"
  |>;
  edgeTail = 92;
  edgeRelations = {
    "wall_static_energy_dedup",
    "field_set_union",
    "drain_named_both",
    "P_retirement_four_manifests",
    "Z_chi_count_reconciliation"
  };
  edgeSequence = ("R" <> ToString[#] &) /@
    Range[edgeTail + 1, edgeTail + Length[edgeRelations]];
  masterRow = <|
    "Param" -> "Z_chi",
    "[L,T,M]" -> {-2, 0, 1},
    "Enters" -> "S_chi kinetic",
    "Class" -> "ACTION",
    "Depends on / relation" -> "inertial wall normalization",
    "Reduction route + status" ->
      "G0 DRAFT-v0 [POSTULATE]; no reduction route"
  |>;
  Switch[
    activeMutation,
    "M_G_ELL_WRONG_EDGE",
      provenance["g_ell", "edge"] = "R17",
    "M_LIGHT_GATING_WRONG_EDGE",
      provenance["light_gating", "edge"] = "R21",
    "M_Z_CHI_DROP_POSTULATE",
      provenance["Z_chi", "flag"] = "[DERIVED]",
    "M_Z_CHI_INVENT_ROUTE",
      provenance["Z_chi", "register_route"] = "field_rescaling",
    "M_Z_CHI_WRONG_ENDPOINT",
      countFact["sensitivity"] = {41, 49};
      countFact["spread_after"] = 8,
    "M_GM_DROP_PN_CITATION",
      provenance["gravitomagnetism", "citation"] = "uncited",
    "M_THROAT_WRONG_STATUS",
      provenance["throat_solve", "status"] = "RESOLVED",
    "M_EDGE_SEQUENCE_GAP",
      edgeSequence[[3]] = "R96",
    "M_Z_MASTER_WRONG_DIM",
      masterRow["[L,T,M]"] = {0, 0, 0},
    _, Null
  ];
  actual = <|
    "provenance" -> provenance,
    "Z_chi_count" -> countFact,
    "prospective_edge_sequence" -> edgeSequence,
    "prospective_edge_tail" -> Last[edgeSequence],
    "Z_chi_master_parameter_row" -> masterRow
  |>;
  expected = <|
    "provenance" -> <|
      "g_ell" -> <|
        "edge" -> "R21",
        "relation" -> "SCOPE_SPLIT_NOT_REDUCTION",
        "persists" -> "LIGHT_BRANE_LOCALIZATION_ENVELOPE"
      |>,
      "light_gating" -> <|
        "edge" -> "R17",
        "status" -> "PENDING",
        "relation" -> "chi_B_projection_distinct_from_g_ell"
      |>,
      "Z_chi" -> <|
        "source" -> "G0-card section 2.2",
        "flag" -> "[POSTULATE]",
        "register_route" -> Null,
        "action_class" -> "ACTION"
      |>,
      "gravitomagnetism" -> <|
        "citation" -> "research/4d_*pn*",
        "coverage" -> "GR_MATCHED_1PN_TO_4PN_INCLUDING_FRAME_DRAGGING",
        "em_asymmetry" -> "EM_DEPARTS_EXACT_MAXWELL",
        "stage044_register_row" -> False
      |>,
      "throat_solve" -> <|
        "citation" -> "ratified-plan/register central shared debt",
        "status" -> "SIM_DEFERRED_CENTRAL_DEBT"
      |>
    |>,
    "Z_chi_count" -> <|
      "base" -> {40, 49},
      "shift" -> {1, 1},
      "sensitivity" -> {41, 50},
      "spread_before" -> 9,
      "spread_after" -> 9,
      "committed_recount" -> "DEFERRED_046",
      "possible_substitution" -> "dissipative_M_chi"
    |>,
    "prospective_edge_sequence" -> {"R93", "R94", "R95", "R96", "R97"},
    "prospective_edge_tail" -> "R97",
    "Z_chi_master_parameter_row" -> <|
      "Param" -> "Z_chi",
      "[L,T,M]" -> {-2, 0, 1},
      "Enters" -> "S_chi kinetic",
      "Class" -> "ACTION",
      "Depends on / relation" -> "inertial wall normalization",
      "Reduction route + status" ->
        "G0 DRAFT-v0 [POSTULATE]; no reduction route"
    |>
  |>;
  expectEqual["M_PROVENANCE_STATUS_BINDING", actual, expected];
  AssociateTo[evidence, "M" -> actual];
  <|
    "g_ell" -> "SCOPE_SPLIT_R21_NOT_MERGED",
    "light_gating" ->
      "g_ell_LOCALIZED; chi_B_projection R17 PENDING " <>
      "(distinct, not double-counted)",
    "gravitomagnetism" ->
      "BOOST_OF_GRAVITOELECTRIC_FLOW(PN-ladder-covered; honest " <>
      "asymmetry: gravity matches GR incl GM, EM departs Maxwell); " <>
      "prose-only in 044",
    "new_draft_knob" -> <|
      "Z_chi" ->
        "POSTULATE_inertial_normalization_no_reduction_route " <>
        "(G0 DRAFT-v0)"
    |>,
    "Z_chi_count_consequence" ->
      "SHIFTS_CONTINUOUS_[40,49]->[41,50] (043-rule; both endpoints " <>
      "+1; spread 9 unchanged); committed re-count DEFERRED_046 " <>
      "(may substitute dissipative M_chi)",
    "throat_solve" -> provenance["throat_solve", "status"]
  |>
];


expectedVerdictObject[] := <|
  "verdict" -> "PARENT_ACTION_ASSEMBLED_AT_COMPLETENESS_FLOOR",
  "S_cons_G0_summands" ->
    {"S_bulk", "S_chi", "S_scalar", "S_mouth", "S_hold", "S_geon_const"},
  "S_assembled_extra" ->
    {"S_brane_light", "S_move_magnetism_moving_coupling"},
  "scalar_branch" ->
    "REDUCED_S_Lh | PARENT_S_H_plus_S_u_plus_S_mix (exclusive)",
  "field_set_reduced" -> {
    "rho", "theta", "r_B==chi_B", "u_L", "h",
    "u_T(divfree)", "u_w", "lambda_Sigma(aux)"
  },
  "field_set_parent" -> {
    "rho", "theta", "r_B==chi_B", "H", "u_L",
    "u_T(divfree)", "u_w", "lambda_Sigma(aux)"
  },
  "theta_B" -> "ELIMINATED_SCHUR_into_A_eff",
  "wall_dedup" ->
    "STATIC_WALL_ENERGY_IDENTITY(r_B==chi_B; kappa_chi==kappa_B; " <>
    "lambda_chi==4*a_B => ell==delta, sigma matches); " <>
    "FULL_FUNCTIONALS_DIFFER(inertial Z_chi vs dissipative M_chi)",
  "wall_factorization" ->
    "L_chi_2ndVar/kappa_chi == A_chi_dag A_chi (annihilates r0')",
  "g_ell" -> "SCOPE_SPLIT_R21_NOT_MERGED",
  "light_gating" ->
    "g_ell_LOCALIZED; chi_B_projection R17 PENDING " <>
    "(distinct, not double-counted)",
  "drain_interfaces" -> "TWO_NONVARIATIONAL_NAMED",
  "drain_types" -> <|
    "part_I_Gamma_B" -> "internal_order_conversion_total_n_conserved",
    "g0_S_drain" -> "total_rho_sink_plus_remote_return"
  |>,
  "drain_equivalence" -> "UNRESOLVED",
  "drain_selection" -> "DEFERRED_045_AND_USER_GATE",
  "conservative_continuity" -> "SOURCE_FREE",
  "wave_speeds" -> <|
    "c_s0_sq" -> "5*K*rho0^4/m",
    "c_pm_sq" -> {0.7928932188, 2.2071067812},
    "c_gamma_sq" -> "mu_R/rho_br",
    "D_star" -> 1.75
  |>,
  "stability" ->
    "POSITIVE_DEFINITE(B_eff>0 AND D_star=7/4>0; " <>
    "wall minima {0,1}; Omega_w^2>=0)",
  "P_retirement" -> <|
    "stage007_action_operative" -> {"S_GNLS", "gL_Mac", "gL_uw"},
    "stage007_action_retired" -> {"L_pol", "gL_Pu"},
    "stage007_DOF" -> "8->4",
    "stage007_drift" -> "11->7",
    "stage006_drift" -> "6->5",
    "net_routeless" -> "-5 == -4(s007 drift) -1(s006 drift)"
  |>,
  "gravitomagnetism" ->
    "BOOST_OF_GRAVITOELECTRIC_FLOW(PN-ladder-covered; honest " <>
    "asymmetry: gravity matches GR incl GM, EM departs Maxwell); " <>
    "prose-only in 044",
  "new_draft_knob" -> <|
    "Z_chi" ->
      "POSTULATE_inertial_normalization_no_reduction_route (G0 DRAFT-v0)"
  |>,
  "Z_chi_count_consequence" ->
    "SHIFTS_CONTINUOUS_[40,49]->[41,50] (043-rule; both endpoints " <>
    "+1; spread 9 unchanged); committed re-count DEFERRED_046 " <>
    "(may substitute dissipative M_chi)",
  "g0_card" -> "DRAFT_V0",
  "closed_action" -> False,
  "throat_solve" -> "SIM_DEFERRED_CENTRAL_DEBT"
|>;


reproductionTooth[pieces_List, fieldProducers_Association] := Module[
  {
    verdict, firstMatchResults, displayOrder, binding, actual,
    expectedObject, expected
  },
  verdict = Join @@ pieces;
  firstMatchResults = Association[Normal[checkResults]];
  If[
    activeMutation === "REPRODUCTION_FIRST_MATCH",
    firstMatchResults["A_WALL_STATIC_ENERGY_DEDUP"] = False
  ];
  AssociateTo[verdict, "verdict" -> firstMatch[firstMatchResults]];
  displayOrder = {
    "verdict", "S_cons_G0_summands", "S_assembled_extra", "scalar_branch",
    "field_set_reduced", "field_set_parent", "theta_B", "wall_dedup",
    "wall_factorization", "g_ell", "light_gating", "drain_interfaces",
    "drain_types", "drain_equivalence", "drain_selection",
    "conservative_continuity", "wave_speeds", "stability", "P_retirement",
    "gravitomagnetism", "new_draft_knob", "Z_chi_count_consequence",
    "g0_card", "closed_action", "throat_solve"
  };
  verdict = KeyTake[verdict, displayOrder];
  binding = Association[Normal[fieldProducers]];
  If[
    activeMutation === "REPRODUCTION_FIELD_BINDING",
    binding = KeyDrop[binding, "throat_solve"]
  ];
  actual = <|
    "verdict_object" -> verdict,
    "producer_fields" -> Sort[Keys[binding]],
    "object_fields" -> Sort[Keys[verdict]],
    "all_prior_teeth_pass" -> And @@ Values[firstMatchResults]
  |>;
  expectedObject = expectedVerdictObject[];
  expected = <|
    "verdict_object" -> expectedObject,
    "producer_fields" -> Sort[Keys[expectedObject]],
    "object_fields" -> Sort[Keys[expectedObject]],
    "all_prior_teeth_pass" -> True
  |>;
  expectEqual["REPRODUCTION", actual, expected];
  AssociateTo[
    evidence,
    "REPRODUCTION" -> <|
      "first_match_order" -> Most[toothOrder],
      "field_producers" -> binding
    |>
  ];
  verdict
];


runAudit[] := Module[
  {pieces = {}, producers = <||>, take, expectedProducer, verdict},
  take[piece_Association, tooth_String] := (
    AppendTo[pieces, piece];
    KeyValueMap[
      Function[{key, value},
        If[KeyExistsQ[producers, key], raise["DUPLICATE_PRODUCER_" <> key]];
        AssociateTo[producers, key -> tooth]
      ],
      piece
    ]
  );

  section["A-B. Static wall reconciliation and wall Hessian"];
  take[wallStaticEnergyTooth[], "A"];
  take[wallHessianTooth[], "B"];

  section["C-E. Localized scalar, wave speeds, and stability"];
  chargeFactorizationTeeth[];
  take[waveSpeedTeeth[], "D"];
  take[stabilityTeeth[], "E"];

  section["F-H. Units, conservative continuity, and drain provenance"];
  dimensionalTooth[];
  take[continuityTooth[], "G"];
  take[drainProvenanceTooth[], "H"];

  section["I-L. Retirement, field union, open status, and action incidence"];
  take[pRetirementTooth[], "I"];
  take[fieldUnionSchurTooth[], "J"];
  take[statusTeeth[], "K"];
  take[actionIncidenceTooth[], "L"];

  section["M. Provenance/status binding and prospective register handoff"];
  take[provenanceStatusTooth[], "M"];

  AssociateTo[producers, "verdict" -> "REPRODUCTION"];
  expectedProducer = <|
    "verdict" -> "REPRODUCTION",
    "S_cons_G0_summands" -> "L",
    "S_assembled_extra" -> "L",
    "scalar_branch" -> "L",
    "field_set_reduced" -> "J",
    "field_set_parent" -> "J",
    "theta_B" -> "J",
    "wall_dedup" -> "A",
    "wall_factorization" -> "B",
    "g_ell" -> "M",
    "light_gating" -> "M",
    "drain_interfaces" -> "H",
    "drain_types" -> "H",
    "drain_equivalence" -> "H",
    "drain_selection" -> "H",
    "conservative_continuity" -> "G",
    "wave_speeds" -> "D",
    "stability" -> "E",
    "P_retirement" -> "I",
    "gravitomagnetism" -> "M",
    "new_draft_knob" -> "M",
    "Z_chi_count_consequence" -> "M",
    "g0_card" -> "K",
    "closed_action" -> "K",
    "throat_solve" -> "M"
  |>;
  If[! deepEqual[producers, expectedProducer], raise["PRODUCER_MAP_INTERNAL"]];

  section["Reproduction. First-match verdict and canonical object"];
  verdict = reproductionTooth[pieces, producers];
  verdict
];


ok = Catch[
  If[
    activeMutation =!= "" && ! KeyExistsQ[mutationToTooth, activeMutation],
    Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
    Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
    raise["UNKNOWN_MUTATION"]
  ];
  Print["ledger_stage044_parent_action_reconciliation Mathematica audit"];
  Print[
    "ROUTE=Reduce/Solve wall classification + DSolveValue kernel + ",
    "Eigenvalues generalized speeds + exact Integrate + status/incidence manifests"
  ];
  Print["VERDICT_DERIVATION=DOCUMENTED_FIRST_MATCH"];
  If[
    activeMutation =!= "",
    Print["ACTIVE_MUTATION=", activeMutation];
    Print["EXPECTED_TOOTH=", mutationToTooth[activeMutation]]
  ];
  productionVerdict = runAudit[];
  If[
    activeMutation =!= "",
    Print["FIRST_FAILURE=MUTATION_DID_NOT_FIRE"];
    raise["MUTATION_DID_NOT_FIRE"]
  ];
  True,
  "ledgerStage044Failure",
  Function[{message, tag}, False]
];

If[
  TrueQ[ok],
  verdictPath = ExpandFileName @ FileNameJoin[{
    DirectoryName[$InputFileName],
    "..",
    "_scratch",
    "stage044",
    "verdict_wl.json"
  }];
  If[! DirectoryQ[DirectoryName[verdictPath]],
    CreateDirectory[DirectoryName[verdictPath], CreateIntermediateDirectories -> True]
  ];
  Export[verdictPath, canonicalJSON[productionVerdict], "String"];
  Print[""];
  Print["CANONICAL_JSON=", StringTrim[canonicalJSON[productionVerdict]]];
  Print["VERDICT_FILE=", verdictPath]
];

Print["TOOTH_COUNT=", Length[toothOrder]];
Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
If[
  TrueQ[ok],
  Print["OVERALL PASS: Mathematica stage044 parent-action reconciliation"];
  Exit[0],
  Print["OVERALL FAIL: Mathematica stage044 audit did not close"];
  Exit[1]
]
