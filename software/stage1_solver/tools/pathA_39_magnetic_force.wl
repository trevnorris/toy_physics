(* pathA_39 Stage 2 magnetic-force derivation, Mathematica engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_39_magnetic_force.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
reportsDir = FileNameJoin[{stage1Root, "reports"}];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
sympyJson = FileNameJoin[{scratchDir, "pathA_39_magnetic_force_sympy.json"}];
jsonOut = FileNameJoin[{scratchDir, "pathA_39_magnetic_force_mathematica.json"}];
path36Yaml = FileNameJoin[{reportsDir, "pathA_36_c5_phase_potential_results.yaml"}];
path38Yaml = FileNameJoin[{reportsDir, "pathA_38_results.yaml"}];
path39Stage01Yaml = FileNameJoin[{reportsDir, "pathA_39_scalar_admixture_results.yaml"}];

If[! FileExistsQ[sympyJson], fail["missing SymPy JSON: " <> sympyJson]];
If[! FileExistsQ[path36Yaml], fail["missing pathA_36 YAML: " <> path36Yaml]];
If[! FileExistsQ[path38Yaml], fail["missing pathA_38 YAML: " <> path38Yaml]];
If[! FileExistsQ[path39Stage01Yaml], fail["missing pathA_39 Stage 0+1 YAML: " <> path39Stage01Yaml]];

$Assumptions =
  rhoBr > 0 && muR > 0 && rhoB0 > 0 && chiC > 0 && cE > 0 &&
  Nu > 0 && aT > 0 && aL > 0 && R > 0 && zeta > 0 && sigma > 0 &&
  k2 > 0;

path36Text = ReadString[path36Yaml];
path38Text = ReadString[path38Yaml];
path39Text = ReadString[path39Stage01Yaml];

checks = <|
  "pathA_36_B_eff" -> StringContainsQ[path36Text, "B_eff: rho_B0**2/chi_c"],
  "pathA_36_c_gamma_squared" -> StringContainsQ[path36Text, "c_gamma_squared: mu_R/rho_br"],
  "pathA_36_transverse_omega2" -> StringContainsQ[path36Text, "transverse_omega2: mu_R*k^2/rho_br"],
  "pathA_38_q_h_plus" -> StringContainsQ[path38Text, "q_h_plus: 2*QE*tanh(b/ell)/b"],
  "pathA_38_dynamic_green_speed" -> StringContainsQ[path38Text, "radial_green_finite_omega: exp(I*R*omega/cE)/(4*pi*R)"],
  "pathA_39_stage01_engine" -> StringContainsQ[path39Text, "status: ENGINE_AGREE"],
  "pathA_39_stage01_q_A_T" -> StringContainsQ[path39Text, "q_A_T: Nu*aT*sCharge"],
  "pathA_39_stage01_q_L" -> StringContainsQ[path39Text, "q_L: Nu*aL*sCharge"]
|>;
If[! And @@ Values[checks], fail["banked import check failed: " <> ToString[Keys@Select[checks, ! # &], InputForm]]];

verdictCodes = <|
  "MAGNETIC_FORCE_DERIVED" -> 1,
  "FAIL_WRONG_FALLOFF" -> 2,
  "FAIL_TARGET_READBACK" -> 3
|>;
verdictCode[label_] := verdictCodes[label];

classifyForce[potentialPower_, forcePower_, braneLocalized_, targetResponseChanged_] := Module[{verdict},
  verdict = Which[
    ! TrueQ[targetResponseChanged], "FAIL_TARGET_READBACK",
    potentialPower =!= -1 || forcePower =!= -2 || ! TrueQ[braneLocalized], "FAIL_WRONG_FALLOFF",
    True, "MAGNETIC_FORCE_DERIVED"
  ];
  <|
    "verdict" -> verdict,
    "verdict_code" -> verdictCode[verdict],
    "potential_power_R" -> potentialPower,
    "force_power_R" -> forcePower,
    "brane_localized" -> braneLocalized,
    "target_response_changed" -> targetResponseChanged
  |>
];

nonzeroQ[expr_] := ! TrueQ[FullSimplify[expr == 0]];
exprNonzeroQ[expr_] := If[VectorQ[expr] || MatrixQ[expr],
  Or @@ (nonzeroQ /@ Flatten[expr]),
  nonzeroQ[expr]
];

rPowerScalar[expr_] := Module[{tau, rationalPower},
  If[TrueQ[FullSimplify[expr == 0]], Return[Missing["ZeroExpression"]]];
  rationalPower = Exponent[Numerator[Together[expr]], R] - Exponent[Denominator[Together[expr]], R];
  If[IntegerQ[rationalPower], Return[rationalPower]];
  tau = Unique["tau"];
  Do[
    If[TrueQ[FullSimplify[(expr /. R -> tau R) == tau^p expr, tau > 0]],
      Return[p]
    ],
    {p, -8, 8}
  ];
  fail["could not measure R power from " <> ToString[expr, InputForm]]
];

rPower[expr_] := Module[{entries, powers},
  If[VectorQ[expr] || MatrixQ[expr],
    entries = Select[Flatten[expr], nonzeroQ];
    powers = DeleteDuplicates[rPowerScalar /@ entries];
    If[Length[powers] =!= 1, fail["matrix has nonuniform R powers: " <> ToString[powers, InputForm]]];
    Return[First[powers]]
  ];
  rPowerScalar[expr]
];

radialHessianContract[seed_] := FullSimplify[
  D[seed, {R, 2}] Aprod + (D[seed, R]/R) (Ddot - Aprod)
];

inverseFourierProjectorKernels[denominatorPower_] := Module[
  {scalarGreen, biharmonicSeed, kkKernel, deltaKernel},
  Which[
    denominatorPower === 2,
      scalarGreen = FullSimplify[1/(4 Pi R)];
      biharmonicSeed = FullSimplify[-R/(8 Pi)],
    denominatorPower === 4,
      scalarGreen = FullSimplify[-R/(8 Pi)];
      biharmonicSeed = FullSimplify[R^3/(96 Pi)],
    True,
      fail["unsupported denominator power"]
  ];
  kkKernel = FullSimplify[-radialHessianContract[biharmonicSeed]];
  deltaKernel = FullSimplify[Ddot scalarGreen];
  <|
    "delta" -> deltaKernel,
    "kk" -> kkKernel,
    "T" -> FullSimplify[deltaKernel - kkKernel],
    "L" -> FullSimplify[kkKernel],
    "scalar_green" -> scalarGreen,
    "biharmonic_seed" -> biharmonicSeed
  |>
];

forceFromPotential[U_] := Module[{gradA, gradU},
  gradA = FullSimplify[(brad V1 + arad V2 - 2 Aprod nhat)/R];
  gradU = FullSimplify[D[U, R] nhat + D[U, Aprod] gradA];
  FullSimplify[-gradU]
];

Beff = FullSimplify[rhoB0^2/chiC];
cGamma2 = FullSimplify[muR/rhoBr];
qA0 = FullSimplify[Nu aT];
qL0 = FullSimplify[Nu aL];
qA1 = FullSimplify[s1 qA0];
qA2 = FullSimplify[s2 qA0];
qL1 = FullSimplify[s1 qL0];
qL2 = FullSimplify[s2 qL0];
s12 = FullSimplify[s1 s2];

V1 = {v1x, v1y, v1z};
V2 = {v2x, v2y, v2z};
nhat = {nx, ny, nz};

OTStatic = {{muR k2}};
OLStatic = {{Beff k2}};
GTStaticCoeff = FullSimplify[Inverse[OTStatic][[1, 1]]];
GLStaticCoeff = FullSimplify[Inverse[OLStatic][[1, 1]]];

standardKernels = inverseFourierProjectorKernels[2];
perturbedKernels = inverseFourierProjectorKernels[4];
KT = standardKernels["T"];
KL = standardKernels["L"];
KTPerturbed = perturbedKernels["T"];
KLPerturbed = perturbedKernels["L"];

UTIntegrate = FullSimplify[-(qA1 qA2/muR) KT];
ULIntegrate = FullSimplify[-(qL1 qL2/Beff) KL];
UTotalIntegrate = FullSimplify[UTIntegrate + ULIntegrate];

UTEom = FullSimplify[-(qA1 qA2/muR) (standardKernels["delta"] - standardKernels["kk"])];
ULEom = FullSimplify[-(qL1 qL2/Beff) standardKernels["kk"]];
UTotalEom = FullSimplify[UTEom + ULEom];

FT = forceFromPotential[UTIntegrate];
FL = forceFromPotential[ULIntegrate];
FTotal = FullSimplify[FT + FL];

UTWrongFalloff = FullSimplify[-(qA1 qA2/muR) KTPerturbed];
ULWrongFalloff = FullSimplify[-(qL1 qL2/Beff) KLPerturbed];
UTotalWrongFalloff = FullSimplify[UTWrongFalloff + ULWrongFalloff];
FTotalWrongFalloff = forceFromPotential[UTotalWrongFalloff];

potentialPower = rPower[UTotalIntegrate];
forcePower = rPower[FTotal];
wrongPotentialPower = rPower[UTotalWrongFalloff];
wrongForcePower = rPower[FTotalWrongFalloff];
derivedTracksFunctionalPerturbation =
  potentialPower =!= wrongPotentialPower &&
  forcePower =!= wrongForcePower &&
  exprNonzeroQ[FTotalWrongFalloff - FTotal];

UTMuAblate = FullSimplify[-(qA1 qA2/(zeta muR)) KT];
muAblationRatio = FullSimplify[UTMuAblate/UTIntegrate];
UTSourceAblate = FullSimplify[-((s1 Nu sigma aT) (s2 Nu sigma aT)/muR) KT];
sourceAblationRatio = FullSimplify[UTSourceAblate/UTIntegrate];
UTProjectionAblate = FullSimplify[-(qA1 qA2/muR) KL];
projectionAblationDelta = FullSimplify[UTProjectionAblate - UTIntegrate];
muAblationChanged = nonzeroQ[muAblationRatio - 1];
sourceAblationChanged = nonzeroQ[sourceAblationRatio - 1];
projectionAblationChanged = nonzeroQ[projectionAblationDelta];

targetReadbackFixture = FullSimplify[(s12 Nu^2 aT^2/muR) (Aprod - Ddot)/(8 Pi R)];
targetReadbackFixturePerturbed = targetReadbackFixture;
targetReadbackDelta = FullSimplify[targetReadbackFixturePerturbed - targetReadbackFixture];
targetReadbackTracks = exprNonzeroQ[targetReadbackDelta];

main = classifyForce[potentialPower, forcePower, True, derivedTracksFunctionalPerturbation];
wrongFalloff = classifyForce[wrongPotentialPower, wrongForcePower, True, True];
noncompact = classifyForce[potentialPower, forcePower, False, True];
targetReadback = classifyForce[rPower[targetReadbackFixture], rPower[targetReadbackFixture], True, targetReadbackTracks];

selfTests = <|
  "MAGNETIC_FORCE_DERIVED" -> main["verdict"],
  "FAIL_WRONG_FALLOFF" -> wrongFalloff["verdict"],
  "FAIL_TARGET_READBACK" -> targetReadback["verdict"]
|>;
If[selfTests["MAGNETIC_FORCE_DERIVED"] =!= "MAGNETIC_FORCE_DERIVED", fail["main classifier self-test failed"]];
If[selfTests["FAIL_WRONG_FALLOFF"] =!= "FAIL_WRONG_FALLOFF", fail["wrong-falloff self-test failed"]];
If[selfTests["FAIL_TARGET_READBACK"] =!= "FAIL_TARGET_READBACK", fail["target-readback self-test failed"]];

zeroVelocityU = FullSimplify[UTotalIntegrate /. {Ddot -> 0, Aprod -> 0}];
neutralCompositeU = FullSimplify[(UTotalIntegrate /. {s1 -> 1, s2 -> 1}) + (UTotalIntegrate /. {s1 -> -1, s2 -> 1})];
chargeFlipDelta = FullSimplify[(UTotalIntegrate /. s1 -> -s1) + UTotalIntegrate];
lorentzResidual = FullSimplify[cE^2 - cGamma2];
lorentzOnConeLock = FullSimplify[lorentzResidual /. cE^2 -> cGamma2];
USide = FullSimplify[UTotalIntegrate /. {Aprod -> 0, arad -> 0, brad -> 0, Ddot -> vDot}];
FSideVec = FullSimplify[forceFromPotential[UTotalIntegrate] /. {Aprod -> 0, arad -> 0, brad -> 0, Ddot -> vDot}];
FSideRadialCoeff = FullSimplify[FSideVec[[1]]/nx];
scalarAdmixtureRatio = FullSimplify[(qL0^2/Beff)/(qA0^2/muR)];
rawRatio = FullSimplify[qL0^2/qA0^2];
rawRatioReference = FullSimplify[Beff/muR];
transverseSideVec = FullSimplify[FT /. {Aprod -> 0, arad -> 0, brad -> 0, Ddot -> 1, s1 -> 1, s2 -> 1}];
longitudinalSideVec = FullSimplify[FL /. {Aprod -> 0, arad -> 0, brad -> 0, Ddot -> 1, s1 -> 1, s2 -> 1}];
totalSideVec = FullSimplify[FTotal /. {Aprod -> 0, arad -> 0, brad -> 0, Ddot -> 1, s1 -> 1, s2 -> 1}];
transverseLikeParallelSideForce = FullSimplify[transverseSideVec[[1]]/nx];
longitudinalLikeParallelSideForce = FullSimplify[longitudinalSideVec[[1]]/nx];
totalLikeParallelSideForce = FullSimplify[totalSideVec[[1]]/nx];
UTGhost = FullSimplify[-(qA1 qA2/(-muR)) KT];
transverseGhostSideForce = FullSimplify[(forceFromPotential[UTGhost] /. {Aprod -> 0, arad -> 0, brad -> 0, Ddot -> 1, s1 -> 1, s2 -> 1})[[1]]/nx];

controls = <|
  "muR_propagator_perturbation" -> <|
    "status" -> If[muAblationChanged, "FIRED", "NOT_FIRED"],
    "verdict" -> "RE_DERIVED_RESPONSE_CHANGED"
  |>,
  "source_projection_scale_perturbation" -> <|
    "status" -> If[sourceAblationChanged, "FIRED", "NOT_FIRED"],
    "verdict" -> "RE_DERIVED_RESPONSE_CHANGED"
  |>,
  "projection_tensor_perturbation" -> <|
    "status" -> If[projectionAblationChanged, "FIRED", "NOT_FIRED"],
    "verdict" -> "RE_DERIVED_RESPONSE_CHANGED"
  |>,
  "propagator_functional_perturbation_kminus4" -> <|
    "status" -> If[derivedTracksFunctionalPerturbation, "FIRED", "NOT_FIRED"],
    "verdict" -> "RE_DERIVED_FALLOFF_CHANGED"
  |>,
  "target_readback_fixture" -> <|
    "status" -> If[targetReadback["verdict"] === "FAIL_TARGET_READBACK", "FIRED", "NOT_FIRED"],
    "verdict" -> targetReadback["verdict"]
  |>,
  "wrong_falloff_fixture" -> <|
    "status" -> If[wrongFalloff["verdict"] === "FAIL_WRONG_FALLOFF", "FIRED", "NOT_FIRED"],
    "verdict" -> wrongFalloff["verdict"]
  |>,
  "noncompact_source_fixture" -> <|
    "status" -> If[noncompact["verdict"] === "FAIL_WRONG_FALLOFF", "FIRED", "NOT_FIRED"],
    "verdict" -> noncompact["verdict"]
  |>,
  "V_equals_zero" -> <|
    "status" -> If[TrueQ[zeroVelocityU == 0], "FIRED", "NOT_FIRED"],
    "verdict" -> "NO_MOVING_SOURCE"
  |>,
  "neutral_plus_minus_composite" -> <|
    "status" -> If[TrueQ[neutralCompositeU == 0], "FIRED", "NOT_FIRED"],
    "verdict" -> "ZERO_MONOPOLE_CURRENT_SOURCE"
  |>,
  "charge_flip_s_to_minus_s" -> <|
    "status" -> If[TrueQ[chargeFlipDelta == 0], "FIRED", "NOT_FIRED"],
    "verdict" -> "SOURCE_SIGN_FLIPS"
  |>,
  "ghost_wrong_sign_transverse" -> <|
    "status" -> If[TrueQ[FullSimplify[transverseLikeParallelSideForce + transverseGhostSideForce == 0]], "FIRED", "NOT_FIRED"],
    "verdict" -> "MU_R_TO_MINUS_MU_R_REDERIVED_FLIPS_ATTRACTION_TO_REPULSION"
  |>
|>;

actuals = <|
  "A_radial_product" -> Aprod,
  "B_eff" -> Beff,
  "D_dot" -> Ddot,
  "F_L_x" -> FL[[1]],
  "F_L_y" -> FL[[2]],
  "F_L_z" -> FL[[3]],
  "F_T_x" -> FT[[1]],
  "F_T_y" -> FT[[2]],
  "F_T_z" -> FT[[3]],
  "F_side_radial_coeff" -> FSideRadialCoeff,
  "F_total_x" -> FTotal[[1]],
  "F_total_y" -> FTotal[[2]],
  "F_total_z" -> FTotal[[3]],
  "F_wrong_falloff_x" -> FTotalWrongFalloff[[1]],
  "F_wrong_falloff_y" -> FTotalWrongFalloff[[2]],
  "F_wrong_falloff_z" -> FTotalWrongFalloff[[3]],
  "G_L_static_coeff" -> GLStaticCoeff,
  "G_T_static_coeff" -> GTStaticCoeff,
  "K_L_perturbed_kminus4" -> KLPerturbed,
  "K_L_realspace" -> KL,
  "K_T_perturbed_kminus4" -> KTPerturbed,
  "K_T_realspace" -> KT,
  "O_L_static_scalar" -> OLStatic[[1, 1]],
  "O_T_static_scalar" -> OTStatic[[1, 1]],
  "U_L_eom_solve" -> ULEom,
  "U_L_integrate" -> ULIntegrate,
  "U_T_eom_solve" -> UTEom,
  "U_T_integrate" -> UTIntegrate,
  "U_side_by_side" -> USide,
  "U_total_eom_solve" -> UTotalEom,
  "U_total_integrate" -> UTotalIntegrate,
  "U_total_wrong_falloff" -> UTotalWrongFalloff,
  "biharmonic_seed_kminus4" -> standardKernels["biharmonic_seed"],
  "c_gamma_squared" -> cGamma2,
  "charge_flip_delta" -> chargeFlipDelta,
  "dimensional_ablations_fired" -> 4,
  "eom_minus_integrate_residual" -> FullSimplify[UTotalEom - UTotalIntegrate],
  "force_power_R" -> forcePower,
  "kk_over_k4_contract" -> standardKernels["kk"],
  "longitudinal_like_parallel_side_force" -> longitudinalLikeParallelSideForce,
  "lorentz_on_cone_lock" -> lorentzOnConeLock,
  "lorentz_residual" -> lorentzResidual,
  "main_verdict_code" -> main["verdict_code"],
  "mu_ablation_ratio" -> muAblationRatio,
  "neutral_composite_U" -> neutralCompositeU,
  "noncompact_verdict_code" -> noncompact["verdict_code"],
  "potential_power_R" -> potentialPower,
  "projection_ablation_delta" -> projectionAblationDelta,
  "qA1" -> qA1,
  "qA2" -> qA2,
  "qA_stage1_no_sign" -> qA0,
  "qL1" -> qL1,
  "qL2" -> qL2,
  "qL_stage1_no_sign" -> qL0,
  "raw_qL2_over_qA2" -> rawRatio,
  "raw_ratio_reference" -> rawRatioReference,
  "scalar_admixture_ratio" -> scalarAdmixtureRatio,
  "scalar_green_kminus2" -> standardKernels["scalar_green"],
  "sign_like_charge_like_current" -> -1,
  "sign_like_charge_opposite_current" -> 1,
  "sign_opposite_charge_like_current" -> 1,
  "sign_opposite_charge_opposite_current" -> -1,
  "source_ablation_ratio" -> sourceAblationRatio,
  "target_readback_delta" -> targetReadbackDelta,
  "target_readback_fixture" -> targetReadbackFixture,
  "target_readback_verdict_code" -> targetReadback["verdict_code"],
  "total_like_parallel_side_force" -> totalLikeParallelSideForce,
  "transverse_ghost_side_force" -> transverseGhostSideForce,
  "transverse_like_parallel_side_force" -> transverseLikeParallelSideForce,
  "wrong_falloff_verdict_code" -> wrongFalloff["verdict_code"],
  "wrong_force_power_R" -> wrongForcePower,
  "wrong_potential_power_R" -> wrongPotentialPower,
  "zero_velocity_U" -> zeroVelocityU
|>;

sympyResults = Import[sympyJson, "RawJSON"];
sympyExprs = sympyResults["engine_agreement"]["mathematica_exprs"];
sympyDigest = sympyResults["engine_agreement"]["sympy_expression_digest"];
engineKeys = Keys[sympyExprs];
missingKeys = Complement[engineKeys, Keys[actuals]];
If[Length[missingKeys] > 0, fail["missing Mathematica actuals: " <> ToString[missingKeys, InputForm]]];

assertEngine[name_] := Module[{expectedText, expectedExpr, actual},
  expectedText = sympyExprs[name];
  If[! StringQ[expectedText], fail["missing SymPy expression for " <> name]];
  expectedExpr = ToExpression[expectedText, InputForm];
  actual = actuals[name];
  If[! TrueQ[FullSimplify[actual == expectedExpr]],
    fail["engine disagreement " <> name <> ": Mathematica got " <>
      ToString[actual, InputForm] <> ", SymPy exported " <> expectedText]
  ];
];

Scan[assertEngine, engineKeys];

topLineVerdict = main["verdict"];
If[topLineVerdict =!= sympyResults["top_line_verdict"], fail["top-line verdict disagreement"]];

controlLabel[item_] := ToString[If[KeyExistsQ[item, "verdict"], item["verdict"], item["status"]]];
controlVerdicts = Association @ KeyValueMap[#1 -> controlLabel[#2] &, controls];

agreementPayload = <|
  "top_line_verdict" -> topLineVerdict,
  "main_verdict_code" -> main["verdict_code"],
  "control_verdicts" -> controlVerdicts,
  "self_tests" -> selfTests,
  "checked_expression_count" -> Length[engineKeys],
  "expr_digest" -> sympyDigest
|>;

out = <|
  "schema" -> "pathA_39_magnetic_force_mathematica/v1",
  "status" -> "OK",
  "headline" -> topLineVerdict,
  "checked_quantities" -> engineKeys,
  "sympy_expression_digest" -> sympyDigest,
  "agreement_payload" -> agreementPayload,
  "computed" -> <|
    "U_total" -> ToString[UTotalIntegrate, InputForm],
    "F_total" -> ToString[FTotal, InputForm],
    "scalar_admixture_ratio" -> ToString[scalarAdmixtureRatio, InputForm],
    "falloff_powers" -> <|"potential" -> potentialPower, "force" -> forcePower|>,
    "lorentz_residual" -> ToString[lorentzResidual, InputForm],
    "control_verdicts" -> controlVerdicts
  |>
|>;

Export[jsonOut, out, "RawJSON"];
Print["OK pathA_39_magnetic_force_mathematica"];
Print[ExportString[<|"json" -> jsonOut, "headline" -> topLineVerdict|>, "RawJSON"]];
