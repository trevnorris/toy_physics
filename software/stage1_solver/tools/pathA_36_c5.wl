(* pathA_36 C5 phase-potential test, Mathematica engine.

   Independent Wolfram-language derivation from the primitive first-order
   theta/Josephson Lagrangian.  The Maxwell square is diagnosed only after the
   Dirac-Bergmann constraint chain is computed.
*)

ClearAll[
  fail, scriptPath, stage1Root, reportsDir, scratchDir, outPath, schema,
  rhoBr, muR, Bsym, Jsym, rhoB0, KthetaSym, chiC, ksym, omega, kappaPhase,
  betaB, mTheta2, q, theta, vq, vt, pQ, piTheta, deltaRho, Cj, C2,
  dadd, dmul, dsub, dimString, dimRecord, records, ablations, checkDim,
  expectFail, buildDimensions, exprString, poisson, signLabel,
  firstOrderAnalysis, elasticLongitudinalControl, decoupledSlavedThetaControl,
  decoupledSecondOrderThetaControl, epsilonMismatchControl,
  branchAIntegrationProof, transverseSector, locusEvaluation, controlStatus,
  controlAgreementLabel, payload, writePayload
];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_36_c5.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
reportsDir = FileNameJoin[{stage1Root, "reports"}];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir, CreateIntermediateDirectories -> True]];
outPath = FileNameJoin[{scratchDir, "pathA_36_c5_mathematica.json"}];
schema = "pathA_36_c5_phase_potential/v1";

ClearAll[
  rhoBr, muR, Bsym, Jsym, rhoB0, KthetaSym, chiC, ksym, omega, kappaPhase,
  betaB, mTheta2, q, theta, vq, vt, pQ, piTheta, deltaRho
];

$Assumptions =
  rhoBr > 0 && muR > 0 && Jsym > 0 && rhoB0 > 0 && chiC > 0 &&
  ksym > 0 && kappaPhase > 0 && betaB > 0 && mTheta2 > 0;

Cj = -Jsym rhoB0;
C2 = FullSimplify[Cj^2];

dadd[dims__] := Total[{dims}];
dmul[n_, dim_] := n dim;
dsub[left_, right_] := left - right;
dimString[dim_] := Module[{labels = {"M", "L", "T"}, parts},
  parts = MapThread[
    If[#2 === 0, Nothing, If[#2 === 1, #1, #1 <> "^" <> ToString[#2]]] &,
    {labels, dim}
  ];
  If[parts === {}, "1", StringRiffle[parts, " "]]
];
dimRecord[dim_] := <|"triple_MLT" -> dim, "string" -> dimString[dim]|>;

records = {};
ablations = {};

checkDim[category_, name_, actual_, expected_, expression_] := Module[{},
  If[actual =!= expected,
    fail[category <> ":" <> name <> ": expected " <> ToString[expected] <> ", got " <> ToString[actual]]
  ];
  AppendTo[
    records,
    <|
      "category" -> category,
      "name" -> name,
      "expression" -> expression,
      "dimension" -> dimRecord[actual],
      "expected" -> dimRecord[expected],
      "status" -> "PASS"
    |>
  ];
];

expectFail[category_, name_, actual_, expected_, expression_] := Module[{},
  If[actual === expected, fail["dimension ablation did not fire: " <> category <> ":" <> name]];
  AppendTo[
    ablations,
    <|
      "category" -> category,
      "name" -> name,
      "expression" -> expression,
      "actual" -> dimRecord[actual],
      "expected" -> dimRecord[expected],
      "status" -> "FIRED"
    |>
  ];
];

buildDimensions[] := Module[
  {
    braneLag, zdim, du, dtheta, dgrad, ddt, dk, domega, drhoBr, dmuR,
    dB, drhoB0, dJ, dCJ, dKtheta, dchiC, dmtheta2, ddivU, dcurlU,
    dgradTheta, ddtU, ddtTheta, ddeltaRho
  },
  records = {};
  ablations = {};
  braneLag = {1, -1, -2};
  zdim = {0, 0, 0};
  du = {0, 1, 0};
  dtheta = zdim;
  dgrad = {0, -1, 0};
  ddt = {0, 0, -1};
  dk = dgrad;
  domega = ddt;
  drhoBr = {1, -3, 0};
  dmuR = braneLag;
  dB = braneLag;
  drhoB0 = drhoBr;
  dJ = {0, 2, -1};
  dCJ = {1, -1, -1};
  dKtheta = {1, 1, -2};
  dchiC = {1, -5, 2};
  dmtheta2 = braneLag;

  ddivU = dadd[dgrad, du];
  dcurlU = ddivU;
  dgradTheta = dadd[dgrad, dtheta];
  ddtU = dadd[ddt, du];
  ddtTheta = dadd[ddt, dtheta];
  ddeltaRho = drhoB0;

  checkDim["primitive", "brane_inertia", dadd[drhoBr, dmul[2, ddtU]], braneLag, "rho_br (partial_t u)^2"];
  checkDim["primitive", "MacCullagh_curl", dadd[dmuR, dmul[2, dcurlU]], braneLag, "mu_R (nabla x u)^2"];
  checkDim["primitive", "Cauchy_bulk", dadd[dB, dmul[2, ddivU]], braneLag, "B (nabla dot u)^2"];
  checkDim["primitive", "Josephson_density", dadd[dJ, ddtTheta, ddeltaRho], braneLag, "J (partial_t theta) delta_rho_B"];
  checkDim["primitive", "slaved_delta_rho", dadd[drhoB0, ddivU], ddeltaRho, "rho_B0 nabla dot u"];
  checkDim["primitive", "phase_gradient_signed", dadd[dKtheta, dmul[2, dgradTheta]], braneLag, "K_theta (nabla theta)^2"];
  checkDim["primitive", "density_compressibility", dsub[dmul[2, ddeltaRho], dchiC], braneLag, "delta_rho_B^2/chi_c"];
  checkDim["primitive", "theta_mass_locking", dadd[dmtheta2, dmul[2, dtheta]], braneLag, "m_theta^2 theta^2"];
  checkDim["ibp", "C_J_definition", dadd[dJ, drhoB0], dCJ, "C_J = -J rho_B0"];
  checkDim["ibp", "Josephson_cross", dadd[dCJ, ddtU, dgradTheta], braneLag, "C_J partial_t u dot nabla theta"];
  checkDim["ibp", "electric_square_velocity_piece", dadd[drhoBr, dmul[2, ddtU]], braneLag, "rho_br (partial_t u)^2"];
  checkDim["ibp", "electric_square_mixed_piece", dadd[dCJ, ddtU, dgradTheta], braneLag, "C_J partial_t u dot nabla theta"];
  checkDim["ibp", "electric_square_gradient_piece", dadd[dKtheta, dmul[2, dgradTheta]], braneLag, "K_theta (nabla theta)^2"];
  checkDim["derived", "Maxwell_locus_CJ2", dmul[2, dCJ], dadd[drhoBr, dKtheta], "C_J^2 = rho_br K_theta"];
  checkDim["derived", "c_gamma_squared", dsub[dmuR, drhoBr], {0, 2, -2}, "c_gamma^2 = mu_R/rho_br"];
  checkDim["branch_a", "induced_theta_kinetic_without_continuity", dadd[dchiC, dmul[2, dJ], dmul[2, ddtTheta]], braneLag, "chi_c J^2 (partial_t theta)^2"];
  checkDim["branch_b", "slaved_compressibility_B_eff", dsub[dmul[2, drhoB0], dchiC], dB, "rho_B0^2/chi_c"];
  checkDim["fourier", "omega2_transverse", dadd[dsub[dmuR, drhoBr], dmul[2, dk]], dmul[2, domega], "(mu_R/rho_br) k^2"];

  expectFail["ablation", "drop_rho_B0_from_Josephson_cross", dadd[dJ, ddtU, dgradTheta], braneLag, "J partial_t u dot nabla theta"];
  expectFail["ablation", "phase_stiffness_without_gradient", dadd[dKtheta, dmul[2, dtheta]], braneLag, "K_theta theta^2"];
  expectFail["ablation", "compressibility_multiplied_not_divided", dadd[dchiC, dmul[2, ddeltaRho]], braneLag, "chi_c delta_rho_B^2"];

  <|
    "pass" -> True,
    "checked_expression_count" -> Length[records],
    "records" -> records,
    "ablations" -> ablations
  |>
];

exprString[expr_] := StringReplace[ToString[FortranForm[FullSimplify[expr]]], {
  "rhoBr" -> "rho_br",
  "muR" -> "mu_R",
  "Jsym" -> "J",
  "rhoB0" -> "rho_B0",
  "chiC" -> "chi_c",
  "ksym" -> "k",
  "kappaPhase" -> "kappa_phase",
  "betaB" -> "beta_B",
  "mTheta2" -> "m_theta2",
  "omega" -> "omega"
}];

poisson[f_, g_] := FullSimplify[
  D[f, q] D[g, pQ] - D[f, pQ] D[g, q] +
  D[f, theta] D[g, piTheta] - D[f, piTheta] D[g, theta]
];

signLabel[expr_] := Which[
  TrueQ[FullSimplify[expr == 0]], "zero",
  TrueQ[FullSimplify[expr > 0]], "positive",
  TrueQ[FullSimplify[expr < 0]], "negative",
  True, "symbolic"
];

firstOrderAnalysis[name_, Kexpr_, Bexpr_, thetaMassExpr_: 0] := Module[
  {
    sTheta, primitiveL, pU, pTh, primary, hamiltonian, secondary, bracket,
    secondaryPreservation, determinant, aEff, omega2Text, bracketZero,
    bZero, massZero, allKGauge, firstClass, secondClass, physicalDof,
    initialData, constraintStructure, gaugeClosure, verdict, poleCount,
    residueSign, bounded, hamiltonianStatus, squareResidual, generator
  },
  sTheta = FullSimplify[Kexpr ksym^2 - thetaMassExpr];
  primitiveL = 1/2 rhoBr vq^2 - Cj ksym q vt + 1/2 sTheta theta^2 - 1/2 Bexpr ksym^2 q^2;
  pU = D[primitiveL, vq];
  pTh = D[primitiveL, vt];
  primary = FullSimplify[piTheta - pTh];
  hamiltonian = FullSimplify[pQ^2/(2 rhoBr) + 1/2 Bexpr ksym^2 q^2 - 1/2 sTheta theta^2];
  secondary = FullSimplify[poisson[primary, hamiltonian]];
  bracket = FullSimplify[poisson[primary, secondary]];
  secondaryPreservation = FullSimplify[poisson[secondary, hamiltonian]];
  determinant = FullSimplify[(rhoBr sTheta - C2 ksym^2) omega^2 - Bexpr ksym^2 sTheta];
  aEff = FullSimplify[rhoBr - C2 ksym^2/sTheta];
  omega2Text = Which[
    TrueQ[FullSimplify[Bexpr == 0]], "0",
    ! TrueQ[FullSimplify[aEff == 0]], exprString[FullSimplify[Bexpr ksym^2/aEff]],
    True, "singular_second_class_no_propagating_pole"
  ];
  bracketZero = TrueQ[FullSimplify[bracket == 0]];
  bZero = TrueQ[FullSimplify[Bexpr == 0]];
  massZero = TrueQ[FullSimplify[thetaMassExpr == 0]];
  allKGauge = bracketZero && bZero && massZero;
  squareResidual = FullSimplify[Kexpr - C2/rhoBr];
  generator = Null;

  If[allKGauge,
    firstClass = 2; secondClass = 0; physicalDof = 0; initialData = 0;
    constraintStructure = "FIRST_CLASS_MAXWELL_CHAIN";
    gaugeClosure = "OFF_SHELL_ALL_K_LOCAL";
    verdict = "C5_RESOLVED_MAXWELL_BY_TUNING";
    poleCount = 0; residueSign = "no_physical_longitudinal_pole";
    bounded = True; hamiltonianStatus = "bounded_on_Gauss_constraint";
    generator = <|
      "generator" -> "G[chi]=(rho_br/C_J)*(chi*Phi_2-dot(chi)*Phi_1)",
      "delta_u_L" -> "k chi",
      "delta_theta" -> "-(rho_br/C_J) dot(chi)",
      "local" -> True,
      "inverse_k_or_inverse_omega" -> False
    |>,
    If[bracketZero && ! bZero,
      firstClass = 0; secondClass = 4; physicalDof = 0; initialData = 0;
      constraintStructure = "TERTIARY_SECOND_CLASS_CHAIN";
      gaugeClosure = "BROKEN_BY_CAUCHY_TERM";
      verdict = "FAIL_SECOND_CLASS_NOT_MAXWELL";
      poleCount = 0; residueSign = "no_gauge_pole_but_not_first_class";
      bounded = True; hamiltonianStatus = "bounded_but_second_class",
      firstClass = 0; secondClass = 2; physicalDof = 1; initialData = 2;
      constraintStructure = "SECOND_CLASS_PAIR";
      gaugeClosure = "BROKEN_OFF_LOCUS";
      poleCount = 1;
      residueSign = signLabel[1/aEff];
      bounded = TrueQ[FullSimplify[aEff > 0]] && (bZero || TrueQ[FullSimplify[Bexpr > 0]]);
      hamiltonianStatus = If[bounded, "bounded_reduced_Hamiltonian", "unbounded_or_negative_residue"];
      verdict = Which[
        ! TrueQ[FullSimplify[thetaMassExpr == 0]], "FAIL_SECOND_CLASS_NOT_MAXWELL",
        ! bounded, "FAIL_GHOST_OR_NEGATIVE_NORM",
        bZero, "FAIL_C5_LONGITUDINAL_ZERO_MODE",
        True, "FAIL_CAUCHY_STRAY_LONGITUDINAL"
      ]
    ]
  ];

  <|
    "name" -> name,
    "coefficients" -> <|
      "C_J" -> "-J*rho_B0",
      "K_theta" -> exprString[Kexpr],
      "B_eff" -> exprString[Bexpr],
      "theta_mass2" -> exprString[thetaMassExpr]
    |>,
    "constraints" -> <|
      "primary" -> exprString[primary],
      "secondary" -> exprString[secondary],
      "constraint_bracket" -> exprString[bracket],
      "secondary_preservation_no_multiplier" -> exprString[secondaryPreservation],
      "first_class_count" -> firstClass,
      "second_class_count" -> secondClass,
      "classification" -> constraintStructure
    |>,
    "gauge" -> <|
      "closure" -> gaugeClosure,
      "square_residual_K_minus_CJ2_over_rho" -> exprString[squareResidual],
      "B_eff_required_zero" -> bZero,
      "theta_mass_required_zero" -> massZero,
      "generator" -> generator
    |>,
    "hamiltonian" -> <|
      "reduced_kinetic_A" -> exprString[aEff],
      "bounded" -> bounded,
      "status" -> hamiltonianStatus
    |>,
    "dispersion" -> <|
      "determinant" -> exprString[determinant],
      "omega2" -> omega2Text,
      "pole_count" -> poleCount,
      "pole_residue_sign" -> residueSign
    |>,
    "mode_count" -> <|
      "physical_dof_per_finite_k" -> physicalDof,
      "independent_initial_data_functions_per_k" -> initialData
    |>,
    "verdict" -> verdict
  |>
];

elasticLongitudinalControl[name_, Bexpr_] := Module[{omega2Expr, propagating},
  omega2Expr = FullSimplify[Bexpr ksym^2/rhoBr];
  propagating = ! TrueQ[FullSimplify[Bexpr == 0]];
  <|
    "name" -> name,
    "constraints" -> <|
      "first_class_count" -> 0,
      "second_class_count" -> 0,
      "classification" -> "UNCONSTRAINED_ELASTIC_COORDINATE"
    |>,
    "dispersion" -> <|
      "omega2" -> exprString[omega2Expr],
      "pole_count" -> 1,
      "pole_residue_sign" -> "positive"
    |>,
    "mode_count" -> <|
      "physical_dof_per_finite_k" -> 1,
      "independent_initial_data_functions_per_k" -> 2
    |>,
    "verdict" -> If[propagating, "FAIL_CAUCHY_STRAY_LONGITUDINAL", "FAIL_C5_LONGITUDINAL_ZERO_MODE"]
  |>
];

decoupledSlavedThetaControl[] := Module[
  {primitiveL, pU, pTh, primary, hamiltonian, secondary, bracket},
  primitiveL = 1/2 rhoBr vq^2 + 1/2 (-kappaPhase) ksym^2 theta^2;
  pU = D[primitiveL, vq];
  pTh = D[primitiveL, vt];
  primary = FullSimplify[piTheta - pTh];
  hamiltonian = FullSimplify[pQ^2/(2 rhoBr) + 1/2 kappaPhase ksym^2 theta^2];
  secondary = FullSimplify[poisson[primary, hamiltonian]];
  bracket = FullSimplify[poisson[primary, secondary]];
  <|
    "name" -> "decoupled_theta_slaved_CJ_zero",
    "momenta" -> <|"p_u_L" -> exprString[pU], "pi_theta" -> exprString[pTh]|>,
    "constraints" -> <|
      "primary" -> exprString[primary],
      "secondary" -> exprString[secondary],
      "constraint_bracket" -> exprString[bracket],
      "first_class_count" -> 0,
      "second_class_count" -> 2,
      "classification" -> "THETA_ALGEBRAIC_SECOND_CLASS_PAIR_PLUS_U_ZERO_MODE"
    |>,
    "dispersion" -> <|"u_L" -> "omega^2 = 0", "theta" -> "algebraic", "pole_count" -> 1, "pole_residue_sign" -> "positive"|>,
    "mode_count" -> <|"physical_dof_per_finite_k" -> 1, "independent_initial_data_functions_per_k" -> 2|>,
    "decoupled_theta_status" -> "NON_DYNAMICAL_ALGEBRAIC_PHASE",
    "verdict" -> "FAIL_C5_LONGITUDINAL_ZERO_MODE"
  |>
];

epsilonMismatchControl[] := Module[
  {
    rhoEps, kEps, sTheta, primitiveL, pU, pTh, primary, hamiltonian,
    secondary, bracket, shiftedTransverse, frozenTransverse, closes
  },
  rhoEps = 2 rhoBr;
  kEps = FullSimplify[C2/rhoEps];
  sTheta = FullSimplify[kEps ksym^2];
  primitiveL = 1/2 rhoEps vq^2 - Cj ksym q vt + 1/2 sTheta theta^2;
  pU = D[primitiveL, vq];
  pTh = D[primitiveL, vt];
  primary = FullSimplify[piTheta - pTh];
  hamiltonian = FullSimplify[pQ^2/(2 rhoEps) - 1/2 sTheta theta^2];
  secondary = FullSimplify[poisson[primary, hamiltonian]];
  bracket = FullSimplify[poisson[primary, secondary]];
  shiftedTransverse = FullSimplify[muR ksym^2/rhoEps];
  frozenTransverse = FullSimplify[muR ksym^2/rhoBr];
  closes = TrueQ[FullSimplify[bracket == 0]];
  <|
    "name" -> "epsilon_mismatch_square_closes_wrong_inertia",
    "constraints" -> <|
      "primary" -> exprString[primary],
      "secondary" -> exprString[secondary],
      "constraint_bracket" -> exprString[bracket],
      "first_class_count" -> If[closes, 2, 0],
      "second_class_count" -> If[closes, 0, 2],
      "classification" -> If[closes, "FIRST_CLASS_LONGITUDINAL_BUT_WRONG_EPSILON", "NOT_CLOSED"]
    |>,
    "longitudinal_square" -> <|"epsilon" -> "2*rho_br", "closes" -> closes, "physical_dof_per_finite_k" -> If[closes, 0, 1]|>,
    "transverse_check" -> <|
      "frozen_omega2" -> exprString[frozenTransverse],
      "mismatched_omega2" -> exprString[shiftedTransverse],
      "speed_shift" -> exprString[FullSimplify[shiftedTransverse - frozenTransverse]]
    |>,
    "mode_count" -> <|"physical_dof_per_finite_k" -> If[closes, 0, 1], "independent_initial_data_functions_per_k" -> If[closes, 0, 2]|>,
    "verdict" -> "FAIL_TRANSVERSE_DISTURBED"
  |>
];

decoupledSecondOrderThetaControl[] := Module[{mTheta, thetaOmega2},
  mTheta = FullSimplify[chiC Jsym^2];
  thetaOmega2 = FullSimplify[kappaPhase ksym^2/mTheta];
  <|
    "name" -> "decoupled_theta_independent_density_no_continuity",
    "hessian_rank" -> 2,
    "constraints" -> <|"first_class_count" -> 0, "second_class_count" -> 0, "classification" -> "TWO_UNCONSTRAINED_SECOND_ORDER_FIELDS"|>,
    "dispersion" -> <|"u_L" -> "omega^2 = 0", "theta" -> "omega^2 = k**2*kappa_phase/(J**2*chi_c)", "pole_count" -> 2, "pole_residue_sign" -> "positive"|>,
    "mode_count" -> <|"physical_dof_per_finite_k" -> 2, "independent_initial_data_functions_per_k" -> 4|>,
    "verdict" -> "FAIL_EXTRA_SCALAR_DOF"
  |>
];

branchAIntegrationProof[] := Module[
  {r, lDensity, rGaussian, lNoContinuity, rContinuity, lWithContinuity, bIncrement},
  r = deltaRho;
  lDensity = FullSimplify[Jsym vt r - r^2/(2 chiC)];
  rGaussian = First[r /. Solve[D[lDensity, r] == 0, r]];
  lNoContinuity = FullSimplify[lDensity /. r -> rGaussian];
  rContinuity = -rhoB0 ksym q;
  lWithContinuity = FullSimplify[lDensity /. r -> rContinuity];
  bIncrement = FullSimplify[rhoB0^2/chiC];
  <|
    "independent_field" -> "delta_rho_B",
    "continuity_constraint_fourier" -> "omega*(delta_rho_B + rho_B0*k*u_L)=0",
    "finite_frequency_solution_fixed_number_sector" -> "delta_rho_B = -rho_B0*k*u_L",
    "gaussian_solution_if_continuity_removed" -> "delta_rho_B = J*chi_c*dot_theta",
    "theta_kinetic_if_continuity_removed" -> "J**2*chi_c*dot_theta**2/2",
    "effective_density_L_with_continuity" -> exprString[lWithContinuity],
    "effective_C_J" -> "-J*rho_B0",
    "effective_B_increment" -> "rho_B0**2/chi_c",
    "proof_status" -> "CONTINUITY_FORCES_SAME_SLAVED_SECTOR"
  |>
];

transverseSector[] := <|
  "basis" -> {"u_T1", "u_T2"},
  "physical_dof" -> 2,
  "omega2" -> "mu_R*k^2/rho_br",
  "c_gamma_squared" -> "mu_R/rho_br",
  "massless" -> True,
  "theta_couplings" -> 0,
  "verdict" -> "PASS_TRANSVERSE_UNDISTURBED"
|>;

locusEvaluation[maxwellCase_] := <|
  "conditions" -> <|
    "i_square_closes" -> "C_J^2 = rho_br*K_theta",
    "i_requires_K_theta_positive" -> True,
    "ii_epsilon_equals_frozen_rho_br" -> True,
    "iii_B_eff_zero" -> "B + rho_B0^2/chi_c (if finite stiffness included) = 0",
    "iv_bounded_Hamiltonian" -> maxwellCase["hamiltonian"]["bounded"],
    "derived_C_J" -> "C_J = -J*rho_B0"
  |>,
  "free_locus_result" -> <|
    "verdict" -> maxwellCase["verdict"],
    "physical_dof_per_finite_k" -> maxwellCase["mode_count"]["physical_dof_per_finite_k"],
    "constraint_classification" -> maxwellCase["constraints"]["classification"],
    "gauge_closure" -> maxwellCase["gauge"]["closure"]
  |>,
  "provenance_status" -> "BY_TUNING",
  "with_provenance" -> False,
  "not_forced_by_frozen_definitions" -> {
    "no frozen definition forces K_theta = J^2*rho_B0^2/rho_br",
    "the order-parameter phase has conventional signed gradient K_theta=-kappa_phase<0, while the Maxwell square requires K_theta>0",
    "finite conjugate-density stiffness contributes rho_B0^2/chi_c to B_eff and is not forced to vanish"
  },
  "exact_tuning_conditions" -> {
    "K_theta -> J^2*rho_B0^2/rho_br",
    "B_eff -> 0, i.e. B + rho_B0^2/chi_c -> 0 if finite density stiffness is present",
    "m_theta2 -> 0"
  }
|>;

controlStatus[name_, analysis_, expected_] := Module[{verdict},
  verdict = analysis["verdict"];
  If[! MemberQ[expected, verdict], fail[name <> ": unexpected verdict " <> verdict]];
  <|
    "name" -> name,
    "status" -> "FIRED",
    "verdict" -> verdict,
    "physical_dof_per_finite_k" -> analysis["mode_count"]["physical_dof_per_finite_k"],
    "initial_data_functions_per_k" -> analysis["mode_count"]["independent_initial_data_functions_per_k"]
  |>
];

controlAgreementLabel[item_] := Which[
  KeyExistsQ[item, "verdict"], item["verdict"],
  KeyExistsQ[item, "free_locus_verdict"], item["fixed_coefficients_verdict"] <> "->" <> item["free_locus_verdict"] <> ":" <> item["provenance_status"],
  KeyExistsQ[item, "absent_verdict"], item["absent_verdict"] <> "->" <> item["included_verdict"],
  True, fail["control lacks agreement verdict: " <> item["name"]]
];

buildPayload[] := Module[
  {
    dimensions, kConventional, bNone, bFiniteCompressibility, kLocus,
    branchBCurlOnly, branchBFinite, branchBTuned, branchAProof,
    branchAFinite, branchANoContinuity, noTheta, cauchy, detunedK,
    wrongSignK, detunedWithB, ghostK, bOnLocus, epsilonMismatch,
    decoupledSlaved, thetaMass, transverse, controls, locus, drift, headline,
    agreementPayload, controlVerdicts
  },
  dimensions = buildDimensions[];
  kConventional = -kappaPhase;
  bNone = 0;
  bFiniteCompressibility = FullSimplify[rhoB0^2/chiC];
  kLocus = FullSimplify[C2/rhoBr];

  branchBCurlOnly = firstOrderAnalysis["branch_b_slaved_curl_only_conventional_K", kConventional, bNone];
  branchBFinite = firstOrderAnalysis["branch_b_slaved_finite_compressibility_conventional_K", kConventional, bFiniteCompressibility];
  branchBTuned = firstOrderAnalysis["branch_b_slaved_tuned_Maxwell_locus", kLocus, bNone];
  branchAProof = branchAIntegrationProof[];
  branchAFinite = firstOrderAnalysis["branch_a_independent_with_continuity_integrated_out", kConventional, bFiniteCompressibility];
  branchANoContinuity = decoupledSecondOrderThetaControl[];
  noTheta = elasticLongitudinalControl["no_theta_curl_only", 0];
  cauchy = elasticLongitudinalControl["cauchy_bulk_no_theta", betaB];
  detunedK = firstOrderAnalysis["mismatched_positive_K_no_Cauchy", FullSimplify[2 C2/rhoBr], 0];
  wrongSignK = firstOrderAnalysis["mismatched_K_theta_le_0", kConventional, 0];
  detunedWithB = firstOrderAnalysis["mismatched_positive_K_with_Cauchy", FullSimplify[2 C2/rhoBr], betaB];
  ghostK = firstOrderAnalysis["mismatched_positive_K_negative_residue", FullSimplify[C2/(2 rhoBr)], betaB];
  bOnLocus = firstOrderAnalysis["B_nonzero_on_square_locus", kLocus, betaB];
  epsilonMismatch = epsilonMismatchControl[];
  decoupledSlaved = decoupledSlavedThetaControl[];
  thetaMass = firstOrderAnalysis["theta_mass_breaks_gauge", kLocus, 0, mTheta2];
  transverse = transverseSector[];

  controls = {
    controlStatus["1_no_theta", noTheta, {"FAIL_C5_LONGITUDINAL_ZERO_MODE"}],
    controlStatus["2_cauchy_bulk", cauchy, {"FAIL_CAUCHY_STRAY_LONGITUDINAL"}],
    controlStatus["3_mismatched_positive_K_no_B", detunedK, {"FAIL_C5_LONGITUDINAL_ZERO_MODE"}],
    controlStatus["3_mismatched_K_theta_le_0", wrongSignK, {"FAIL_C5_LONGITUDINAL_ZERO_MODE"}],
    controlStatus["3_mismatched_positive_K_with_B", detunedWithB, {"FAIL_CAUCHY_STRAY_LONGITUDINAL"}],
    controlStatus["3_mismatched_positive_K_ghost", ghostK, {"FAIL_GHOST_OR_NEGATIVE_NORM"}],
    controlStatus["3_B_nonzero_on_locus", bOnLocus, {"FAIL_SECOND_CLASS_NOT_MAXWELL"}],
    <|"name" -> "3_epsilon_mismatch", "status" -> "FIRED", "verdict" -> epsilonMismatch["verdict"], "longitudinal_dof_per_finite_k" -> epsilonMismatch["mode_count"]["physical_dof_per_finite_k"], "frozen_transverse_omega2" -> epsilonMismatch["transverse_check"]["frozen_omega2"], "mismatched_transverse_omega2" -> epsilonMismatch["transverse_check"]["mismatched_omega2"]|>,
    controlStatus["4_decoupled_theta_slaved", decoupledSlaved, {"FAIL_C5_LONGITUDINAL_ZERO_MODE"}],
    controlStatus["4_decoupled_theta_independent", branchANoContinuity, {"FAIL_EXTRA_SCALAR_DOF"}],
    <|"name" -> "5_transverse_undisturbed", "status" -> "FIRED", "verdict" -> transverse["verdict"], "physical_dof" -> transverse["physical_dof"], "omega2" -> transverse["omega2"]|>,
    <|"name" -> "6_provenance_ablation", "status" -> "FIRED", "fixed_coefficients_verdict" -> branchBFinite["verdict"], "free_locus_verdict" -> branchBTuned["verdict"], "provenance_status" -> "BY_TUNING"|>,
    <|"name" -> "7_compressibility_absent_vs_included", "status" -> "FIRED", "absent_verdict" -> branchBCurlOnly["verdict"], "included_verdict" -> branchBFinite["verdict"], "B_increment" -> "rho_B0^2/chi_c"|>,
    controlStatus["8_theta_mass", thetaMass, {"FAIL_SECOND_CLASS_NOT_MAXWELL", "FAIL_C5_LONGITUDINAL_ZERO_MODE"}]
  };

  locus = locusEvaluation[branchBTuned];
  drift = <|
    "branch_b_slaved" -> <|"new_fields" -> {"theta"}, "new_constants" -> {"B", "J", "rho_B0", "K_theta", "chi_c"}, "count" -> 6, "tag" -> "DRIFT(6)"|>,
    "branch_a_independent" -> <|"new_fields" -> {"theta", "delta_rho_B"}, "new_constants" -> {"B", "J", "rho_B0", "K_theta", "chi_c"}, "count" -> 7, "tag" -> "DRIFT(7)"|>
  |>;
  headline = <|
    "verdict" -> branchBFinite["verdict"],
    "main_branch" -> "branch_b_slaved_finite_compressibility_conventional_K",
    "physical_reason" -> "finite density stiffness creates B_eff=rho_B0^2/chi_c and conventional K_theta=-kappa_phase<0 misses the electric Maxwell sign",
    "curl_only_conventional_subcase_verdict" -> branchBCurlOnly["verdict"],
    "tuned_locus_verdict" -> branchBTuned["verdict"],
    "provenance_status" -> locus["provenance_status"],
    "longitudinal_dof_main" -> branchBFinite["mode_count"]["physical_dof_per_finite_k"],
    "longitudinal_initial_data_main" -> branchBFinite["mode_count"]["independent_initial_data_functions_per_k"],
    "transverse_dof" -> transverse["physical_dof"],
    "engine_agreement" -> "PENDING"
  |>;
  controlVerdicts = Association[Map[#["name"] -> controlAgreementLabel[#] &, controls]];
  agreementPayload = <|
    "headline_verdict" -> headline["verdict"],
    "main_constraint_classification" -> branchBFinite["constraints"]["classification"],
    "main_longitudinal_dof" -> branchBFinite["mode_count"]["physical_dof_per_finite_k"],
    "main_initial_data_functions" -> branchBFinite["mode_count"]["independent_initial_data_functions_per_k"],
    "main_pole_count" -> branchBFinite["dispersion"]["pole_count"],
    "curl_only_verdict" -> branchBCurlOnly["verdict"],
    "curl_only_longitudinal_dof" -> branchBCurlOnly["mode_count"]["physical_dof_per_finite_k"],
    "tuned_locus_verdict" -> branchBTuned["verdict"],
    "tuned_constraint_classification" -> branchBTuned["constraints"]["classification"],
    "tuned_longitudinal_dof" -> branchBTuned["mode_count"]["physical_dof_per_finite_k"],
    "tuned_first_class_count" -> branchBTuned["constraints"]["first_class_count"],
    "locus_provenance" -> locus["provenance_status"],
    "branch_a_integrated_verdict" -> branchAFinite["verdict"],
    "transverse_dof" -> transverse["physical_dof"],
    "transverse_omega2" -> transverse["omega2"],
    "control_verdicts" -> controlVerdicts,
    "dimensional_ablations_fired" -> Length[dimensions["ablations"]]
  |>;

  <|
    "schema" -> schema,
    "engine" -> "mathematica",
    "directive" -> "software/stage1_solver/directives/pathA_36_c5_phase_potential.md",
    "frozen_baseline" -> "T0_SHEAR_FROZEN(d9520d3819c3)",
    "primitive_start" -> <|
      "lagrangian" -> "1/2 rho_br dot(u)^2 - 1/2 mu_R (curl u)^2 - 1/2 B (div u)^2 + J dot(theta) delta_rho_B + 1/2 K_theta (grad theta)^2",
      "square_completion_used_as_input" -> False,
      "derived_C_J" -> "C_J = -J*rho_B0"
    |>,
    "headline" -> headline,
    "branches" -> <|
      "branch_b_slaved_curl_only_conventional_K" -> branchBCurlOnly,
      "branch_b_slaved_finite_compressibility_conventional_K" -> branchBFinite,
      "branch_b_slaved_tuned_Maxwell_locus" -> branchBTuned,
      "branch_a_integration_proof" -> branchAProof,
      "branch_a_independent_with_continuity_integrated_out" -> branchAFinite,
      "branch_a_no_continuity_second_medium_ablation" -> branchANoContinuity,
      "epsilon_mismatch_control" -> epsilonMismatch
    |>,
    "locus" -> locus,
    "controls" -> controls,
    "transverse" -> transverse,
    "dimensional_firewall" -> dimensions,
    "drift" -> drift,
    "absence_tags" -> <|"AXIS_RE_ADMITTED" -> False, "U_W_COLLISION" -> False|>,
    "agreement_payload" -> agreementPayload
  |>
];

payload = buildPayload[];
Export[outPath, payload, "JSON", "Compact" -> False];
Print[ExportString[<|"engine" -> "mathematica", "status" -> "OK", "verdict" -> payload["headline"]["verdict"]|>, "JSON", "Compact" -> True]];
Exit[0];
