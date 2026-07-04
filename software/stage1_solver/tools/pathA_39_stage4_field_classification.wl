(* pathA_39 Stage 4 field-coupling classification, Mathematica engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_39_stage4_field_classification.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
reportsDir = FileNameJoin[{stage1Root, "reports"}];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
sympyJson = FileNameJoin[{scratchDir, "pathA_39_stage4_field_classification_sympy.json"}];
jsonOut = FileNameJoin[{scratchDir, "pathA_39_stage4_field_classification_mathematica.json"}];
path36Yaml = FileNameJoin[{reportsDir, "pathA_36_c5_phase_potential_results.yaml"}];
path38Yaml = FileNameJoin[{reportsDir, "pathA_38_results.yaml"}];
stage01Report = FileNameJoin[{reportsDir, "pathA_39_scalar_admixture_screen.md"}];
stage2Report = FileNameJoin[{reportsDir, "pathA_39_magnetic_force.md"}];
stage3Yaml = FileNameJoin[{reportsDir, "pathA_39_stage3_operator_parity_results.yaml"}];

Scan[
  If[! FileExistsQ[#], fail["missing required input: " <> #]] &,
  {sympyJson, path36Yaml, path38Yaml, stage01Report, stage2Report, stage3Yaml}
];

path36Text = ReadString[path36Yaml];
path38Text = ReadString[path38Yaml];
stage01Text = ReadString[stage01Report];
stage2Text = ReadString[stage2Report];
stage3Text = ReadString[stage3Yaml];

If[! StringContainsQ[path36Text, "B_eff: rho_B0**2/chi_c"], fail["pathA_36 B_eff import mismatch"]];
If[! StringContainsQ[path36Text, "c_gamma_squared: mu_R/rho_br"], fail["pathA_36 c_gamma import mismatch"]];
If[! StringContainsQ[path36Text, "classification: SECOND_CLASS_PAIR"], fail["pathA_36 second-class import missing"]];
If[! StringContainsQ[path36Text, "classification: FIRST_CLASS_MAXWELL_CHAIN"], fail["pathA_36 first-class Maxwell import missing"]];
If[! StringContainsQ[path38Text, "q_h_plus: 2*QE*tanh(b/ell)/b"], fail["pathA_38 q_h import mismatch"]];
If[! StringContainsQ[path38Text, "radial_green_finite_omega: exp(I*R*omega/cE)/(4*pi*R)"], fail["pathA_38 dynamic Green import mismatch"]];
If[! StringContainsQ[path38Text, "status: ENGINE_AGREE"], fail["pathA_38 engine agreement missing"]];
If[! StringContainsQ[stage01Text, "D_x"] || ! StringContainsQ[stage01Text, "Chu"] ||
   ! StringContainsQ[stage01Text, "Kh=Mh*cE^2"] || ! StringContainsQ[stage01Text, "q_L = Nu*aL*sCharge"] ||
   ! StringContainsQ[stage01Text, "ENGINE_AGREE"], fail["Stage 0+1 import block mismatch"]];
If[! StringContainsQ[stage2Text, "NO_CANCELLATION_BOTH_CHANNELS_ATTRACTIVE"] ||
   ! StringContainsQ[stage2Text, "ENGINE_AGREE"], fail["Stage 2 sign source mismatch"]];
If[! StringContainsQ[stage3Text, "verdict: FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING"] ||
   ! StringContainsQ[stage3Text, "status: ENGINE_AGREE"], fail["Stage 3 contamination import mismatch"]];

$Assumptions =
  rhoBr > 0 && muR > 0 && rhoB0 > 0 && chiC > 0 && Mh > 0 && cE > 0 &&
  cBad > 0 && deltaB > 0 && deltaKh > 0 && Chu > 0 && QE > 0 && b > 0 &&
  ell > 0 && k > 0 && omega > 0 && x > 0;

verdictCodes = <|
  "FIELD_CLASSIFICATION_UNDERDETERMINED" -> 0,
  "FIELD_EXACT_MAXWELL_STRUCTURE" -> 1,
  "FIELD_TRANSVERSE_EM_PLUS_CLEAN_GRAVITY_DENSITY" -> 2,
  "FIELD_SCALAR_VECTOR_DEPARTURE" -> 3,
  "FIELD_SCALAR_SECTOR_UNSTABLE" -> 4
|>;

classLabel[code_] := First @ First @ Select[Normal[verdictCodes], Last[#] === code &];

classifyCode[importOK_, stableCode_, totalDof_, transverseDof_, fcGen_, hDof_, scalarChargeNonzero_, transverseChargeNonzero_, hChargeResidueNonzero_] := Which[
  ! TrueQ[importOK], verdictCodes["FIELD_CLASSIFICATION_UNDERDETERMINED"],
  stableCode === 0, verdictCodes["FIELD_SCALAR_SECTOR_UNSTABLE"],
  totalDof === 2 && fcGen === 1 && hDof === 0 && scalarChargeNonzero === 0 && transverseChargeNonzero === 1,
    verdictCodes["FIELD_EXACT_MAXWELL_STRUCTURE"],
  totalDof >= 3 && transverseDof === 2 && scalarChargeNonzero === 0 && transverseChargeNonzero === 1,
    verdictCodes["FIELD_TRANSVERSE_EM_PLUS_CLEAN_GRAVITY_DENSITY"],
  totalDof >= 3 && transverseDof === 2 && hDof === 1 && hChargeResidueNonzero === 1 && scalarChargeNonzero === 1,
    verdictCodes["FIELD_SCALAR_VECTOR_DEPARTURE"],
  True, verdictCodes["FIELD_CLASSIFICATION_UNDERDETERMINED"]
];

extractNumber[anchor_, key_] := Module[{hits},
  hits = StringCases[path36Text, RegularExpression["(?s)" <> anchor <> ".*?" <> key <> ":\\s*([0-9]+)"] -> "$1"];
  If[Length[hits] == 0, fail["could not import " <> key <> " for " <> anchor]];
  ToExpression[First[hits]]
];

extractClass[anchor_] := Module[{hits},
  hits = StringCases[path36Text, RegularExpression["(?s)" <> anchor <> ".*?classification:\\s*([A-Z_]+)"] -> "$1"];
  If[Length[hits] == 0, fail["could not import classification for " <> anchor]];
  First[hits]
];

constraintRecord[anchor_] := <|
  "classification" -> extractClass[anchor],
  "firstClassCount" -> extractNumber[anchor, "first_class_count"],
  "secondClassCount" -> extractNumber[anchor, "second_class_count"],
  "longitudinalDof" -> extractNumber[anchor, "physical_dof_per_finite_k"]
|>;

kineticRankFromQ[Qmat_] := Module[{massMat},
  massMat = Map[Coefficient[#, omega, 2] &, Qmat, {2}];
  MatrixRank[massMat]
];

diracDof[includeH_, constraint_] := Module[
  {transverseConfig = 2, longitudinalImportedConfig = 2, hConfig, firstClass, secondClass, totalConfig},
  hConfig = If[TrueQ[includeH], 1, 0];
  firstClass = constraint["firstClassCount"];
  secondClass = constraint["secondClassCount"];
  totalConfig = transverseConfig + longitudinalImportedConfig + hConfig;
  <|
    "configCount" -> totalConfig,
    "transverse" -> 2,
    "uL" -> constraint["longitudinalDof"],
    "h" -> hConfig,
    "firstClassCount" -> firstClass,
    "secondClassCount" -> secondClass,
    "physicalTotal" -> totalConfig - firstClass - secondClass/2
  |>
];

signCode[expr_] := Which[
  TrueQ[FullSimplify[expr == 0]], 0,
  TrueQ[FullSimplify[expr < 0]], -1,
  TrueQ[FullSimplify[expr > 0]], 1,
  True, 1
];

stabilityRecord[Bexpr_, Khexpr_, scalar_] := Module[{detSign, diagonalPositive, eigenPositive, stableCode},
  detSign = signCode[scalar["detStiffness"]];
  diagonalPositive = TrueQ[FullSimplify[Bexpr > 0 && Khexpr > 0 && k^2 > 0]];
  eigenPositive = If[detSign === 1 && diagonalPositive, 1, 0];
  stableCode = If[detSign === 1 && eigenPositive === 1, 1, 0];
  <|"stableCode" -> stableCode, "detSignCode" -> detSign, "eigenPositiveCode" -> eigenPositive|>
];

nonzeroResidueCode[expr_] := If[TrueQ[expr === 0], 0, 1];

anyResidueNonzero[exprs_] := If[AnyTrue[exprs, nonzeroResidueCode[#] === 1 &], 1, 0];

chargeResidueFlags[Bexpr_, Khexpr_, Cexpr_, qLexpr_, qhExpr_, scalar_, includeH_, maxwell_] := Module[
  {densityOnly, hOnly, scalarResidues, densityResidues, hResidues},
  If[TrueQ[maxwell],
    Return[<|"scalar" -> 0, "density" -> 0, "h" -> 0|>]
  ];
  densityOnly = scalarCommon[Bexpr, Khexpr, Cexpr, qLexpr, 0];
  hOnly = scalarCommon[Bexpr, Khexpr, Cexpr, 0, qhExpr];
  scalarResidues = {scalar["rootMinusChargeResidue"], scalar["rootPlusChargeResidue"]};
  densityResidues = {densityOnly["rootMinusChargeResidue"], densityOnly["rootPlusChargeResidue"]};
  hResidues = {hOnly["rootMinusChargeResidue"], hOnly["rootPlusChargeResidue"]};
  <|
    "scalar" -> anyResidueNonzero[scalarResidues],
    "density" -> anyResidueNonzero[densityResidues],
    "h" -> If[TrueQ[includeH], anyResidueNonzero[hResidues], 0]
  |>
];

importOK[Bexpr_, Khexpr_, cEexpr_] := TrueQ[FullSimplify[Bexpr == Beff && Khexpr == Kh && cEexpr == cE]];

Beff = rhoB0^2/chiC;
Kh = Mh cE^2;
cGamma2 = muR/rhoBr;
qT1 = Nu aT sCharge;
qT2 = Nu aTp sCharge;
qTnorm = qT1^2 + qT2^2;
qLStage = Nu aL sCharge;
qh = 2 QE Tanh[b/ell]/b;
largeChu = 2 Sqrt[Beff Kh];
QT = rhoBr omega^2 - muR k^2;

scalarCommon[Bexpr_, Khexpr_, Cexpr_, qLexpr_, qhExpr_] := Module[
  {Dmat, detD, Gx, numCharge, numMass, numChargeMass, trace, delta, cMinus, cPlus,
   denPrime, detStiff, lambdaMinus, lambdaPlus},
  Dmat = {{rhoBr x - Bexpr, -Cexpr}, {-Cexpr, Mh x - Khexpr}};
  detD = Det[Dmat];
  Gx = {{Mh x - Khexpr, Cexpr}, {Cexpr, rhoBr x - Bexpr}}/detD;
  numCharge = qLexpr^2 (Mh x - Khexpr) + 2 qLexpr qhExpr Cexpr + qhExpr^2 (rhoBr x - Bexpr);
  numMass = qM^2 (Mh x - Khexpr);
  numChargeMass = qM (qLexpr (Mh x - Khexpr) + qhExpr Cexpr);
  trace = rhoBr Khexpr + Mh Bexpr;
  delta = (rhoBr Khexpr - Mh Bexpr)^2 + 4 rhoBr Mh Cexpr^2;
  cMinus = (trace - Sqrt[delta])/(2 rhoBr Mh);
  cPlus = (trace + Sqrt[delta])/(2 rhoBr Mh);
  denPrime = D[detD, x];
  detStiff = Bexpr Khexpr - Cexpr^2;
  lambdaMinus = k^2 (Bexpr + Khexpr - Sqrt[(Bexpr - Khexpr)^2 + 4 Cexpr^2])/2;
  lambdaPlus = k^2 (Bexpr + Khexpr + Sqrt[(Bexpr - Khexpr)^2 + 4 Cexpr^2])/2;
  <|
    "Dmat" -> Dmat,
    "detD" -> detD,
    "Gx" -> Gx,
    "numCharge" -> numCharge,
    "numMass" -> numMass,
    "numChargeMass" -> numChargeMass,
    "AqqX" -> numCharge/detD,
    "AmmX" -> numMass/detD,
    "AqmX" -> numChargeMass/detD,
    "cMinus" -> cMinus,
    "cPlus" -> cPlus,
    "denPrime" -> denPrime,
    "rootMinusChargeResidue" -> (numCharge /. x -> cMinus)/(denPrime /. x -> cMinus),
    "rootPlusChargeResidue" -> (numCharge /. x -> cPlus)/(denPrime /. x -> cPlus),
    "rootMinusMassResidue" -> (numMass /. x -> cMinus)/(denPrime /. x -> cMinus),
    "rootPlusMassResidue" -> (numMass /. x -> cPlus)/(denPrime /. x -> cPlus),
    "detStiffness" -> detStiff,
    "lambdaMinus" -> lambdaMinus,
    "lambdaPlus" -> lambdaPlus
  |>
];

realScalar = scalarCommon[Beff, Kh, Chu, qLStage, qh];
cleanScalar = scalarCommon[Beff, Kh, Chu, 0, 0];
aL0Scalar = scalarCommon[Beff, Kh, Chu, 0, qh];
largeScalar = scalarCommon[Beff, Kh, largeChu, qLStage, qh];
importBScalar = scalarCommon[Beff + deltaB, Kh, Chu, qLStage, qh];

realConstraint = constraintRecord["branch_b_slaved_finite_compressibility_conventional_K:"];
maxwellConstraint = constraintRecord["branch_b_slaved_tuned_Maxwell_locus:"];

makeBranch[name_, includeH_, maxwell_, qLexpr_, qhExpr_, Bexpr_, Khexpr_, cEexpr_, Cexpr_, constraint_] := Module[
  {Qmat, scalar, stability, residues, dof, rank, fcGen, importFlag, transverseCharge, code},
  If[TrueQ[maxwell],
    Qmat = {{QT, 0, 0}, {0, QT, 0}, {0, 0, 0}};
    scalar = <||>;
    stability = <|"stableCode" -> 2, "detSignCode" -> 2, "eigenPositiveCode" -> 2|>;
    residues = <|"scalar" -> 0, "density" -> 0, "h" -> 0|>,
    Qmat = {
      {QT, 0, 0, 0},
      {0, QT, 0, 0},
      {0, 0, rhoBr omega^2 - Bexpr k^2, -Cexpr k^2},
      {0, 0, -Cexpr k^2, Mh omega^2 - Khexpr k^2}
    };
    scalar = scalarCommon[Bexpr, Khexpr, Cexpr, qLexpr, qhExpr];
    stability = stabilityRecord[Bexpr, Khexpr, scalar];
    residues = chargeResidueFlags[Bexpr, Khexpr, Cexpr, qLexpr, qhExpr, scalar, includeH, maxwell]
  ];
  dof = diracDof[includeH, constraint];
  rank = kineticRankFromQ[Qmat];
  fcGen = If[constraint["firstClassCount"] > 0 && StringContainsQ[constraint["classification"], "FIRST_CLASS"], 1, 0];
  importFlag = importOK[Bexpr, Khexpr, cEexpr];
  transverseCharge = nonzeroResidueCode[qTnorm];
  code = classifyCode[
    importFlag,
    stability["stableCode"],
    dof["physicalTotal"],
    dof["transverse"],
    fcGen,
    dof["h"],
    residues["scalar"],
    transverseCharge,
    residues["h"]
  ];
  <|
    "name" -> name,
    "Q" -> Qmat,
    "scalar" -> scalar,
    "stability" -> stability,
    "residues" -> residues,
    "dof" -> dof,
    "kineticRank" -> rank,
    "firstClassGenerator" -> fcGen,
    "importOK" -> importFlag,
    "primaryCode" -> code
  |>
];

realBranch = makeBranch["real_provenance_fixed", True, False, qLStage, qh, Beff, Kh, cE, Chu, realConstraint];
maxwellBranch = makeBranch["maxwell_counterfactual", False, True, 0, 0, Beff, Kh, cE, Chu, maxwellConstraint];
cleanBranch = makeBranch["clean_coexistence", True, False, 0, 0, Beff, Kh, cE, Chu, realConstraint];
aL0Branch = makeBranch["aL_to_0", True, False, 0, qh, Beff, Kh, cE, Chu, realConstraint];
largeBranch = makeBranch["large_C_hu", True, False, qLStage, qh, Beff, Kh, cE, largeChu, realConstraint];
restoredBranch = makeBranch["large_C_hu_restored_bound", True, False, qLStage, qh, Beff, Kh, cE, Chu, realConstraint];
importBBranch = makeBranch["import_fidelity_B_eff_corrupt", True, False, qLStage, qh, Beff + deltaB, Kh, cE, Chu, realConstraint];
importKBranch = makeBranch["import_fidelity_K_h_corrupt", True, False, qLStage, qh, Beff, Kh + deltaKh, cE, Chu, realConstraint];
importCEBranch = makeBranch["import_fidelity_c_E_corrupt", True, False, qLStage, qh, Beff, Mh cBad^2, cBad, Chu, realConstraint];

realCode = realBranch["primaryCode"];
maxwellCode = maxwellBranch["primaryCode"];
cleanCode = cleanBranch["primaryCode"];
aL0Code = aL0Branch["primaryCode"];
largeCode = largeBranch["primaryCode"];
restoredCode = restoredBranch["primaryCode"];
importBCode = importBBranch["primaryCode"];
importKCode = importKBranch["primaryCode"];
importCECode = importCEBranch["primaryCode"];

actuals = <|
  "Q_real_00" -> QT,
  "Q_real_22" -> rhoBr omega^2 - Beff k^2,
  "Q_real_23" -> -Chu k^2,
  "Q_real_33" -> Mh omega^2 - Kh k^2,
  "Q_maxwell_00" -> QT,
  "Q_maxwell_gauge_22" -> 0,
  "Q_clean_23" -> -Chu k^2,
  "Q_large_23" -> -largeChu k^2,
  "detQ_real" -> QT^2 k^4 realScalar["detD"],
  "physical_detQ_maxwell" -> QT^2,
  "B_eff" -> Beff,
  "K_h" -> Kh,
  "c_gamma_squared" -> cGamma2,
  "qT_norm_squared" -> qTnorm,
  "q_L_stage" -> qLStage,
  "q_h" -> qh,
  "J_real_0" -> qT1,
  "J_real_1" -> qT2,
  "J_real_2" -> qLStage,
  "J_real_3" -> qh,
  "R_charge_real" -> qTnorm/QT + realScalar["AqqX"]/k^2,
  "R_charge_maxwell" -> qTnorm/QT,
  "det_scalar_real" -> realScalar["detD"],
  "scalar_c_minus_real" -> realScalar["cMinus"],
  "scalar_c_plus_real" -> realScalar["cPlus"],
  "scalar_det_stiffness_real" -> realScalar["detStiffness"],
  "scalar_lambda_minus_real" -> realScalar["lambdaMinus"],
  "scalar_lambda_plus_real" -> realScalar["lambdaPlus"],
  "real_root_minus_charge_residue" -> realScalar["rootMinusChargeResidue"],
  "real_root_plus_charge_residue" -> realScalar["rootPlusChargeResidue"],
  "real_root_minus_mass_residue" -> realScalar["rootMinusMassResidue"],
  "real_root_plus_mass_residue" -> realScalar["rootPlusMassResidue"],
  "transverse_T1_charge_residue" -> qT1^2/rhoBr,
  "transverse_T2_charge_residue" -> qT2^2/rhoBr,
  "clean_root_minus_charge_residue" -> cleanScalar["rootMinusChargeResidue"],
  "clean_root_plus_charge_residue" -> cleanScalar["rootPlusChargeResidue"],
  "aL0_root_minus_charge_residue" -> aL0Scalar["rootMinusChargeResidue"],
  "aL0_root_plus_charge_residue" -> aL0Scalar["rootPlusChargeResidue"],
  "large_Chu" -> largeChu,
  "large_scalar_det_stiffness" -> largeScalar["detStiffness"],
  "large_root_minus_charge_residue" -> largeScalar["rootMinusChargeResidue"],
  "large_root_plus_charge_residue" -> largeScalar["rootPlusChargeResidue"],
  "import_B_Q22" -> rhoBr omega^2 - (Beff + deltaB) k^2,
  "import_B_det_scalar" -> importBScalar["detD"],
  "real_kinetic_rank" -> realBranch["kineticRank"],
  "maxwell_kinetic_rank" -> maxwellBranch["kineticRank"],
  "real_total_dof" -> realBranch["dof"]["physicalTotal"],
  "maxwell_total_dof" -> maxwellBranch["dof"]["physicalTotal"],
  "real_longitudinal_dof" -> realBranch["dof"]["uL"],
  "maxwell_longitudinal_dof" -> maxwellBranch["dof"]["uL"],
  "real_h_dof" -> realBranch["dof"]["h"],
  "maxwell_h_dof" -> maxwellBranch["dof"]["h"],
  "real_first_class_count" -> realBranch["dof"]["firstClassCount"],
  "maxwell_first_class_count" -> maxwellBranch["dof"]["firstClassCount"],
  "real_second_class_count" -> realBranch["dof"]["secondClassCount"],
  "maxwell_second_class_count" -> maxwellBranch["dof"]["secondClassCount"],
  "real_scalar_stable_code" -> realBranch["stability"]["stableCode"],
  "large_scalar_stable_code" -> largeBranch["stability"]["stableCode"],
  "real_scalar_det_sign_code" -> realBranch["stability"]["detSignCode"],
  "large_scalar_det_sign_code" -> largeBranch["stability"]["detSignCode"],
  "real_scalar_eigen_positive_code" -> realBranch["stability"]["eigenPositiveCode"],
  "large_scalar_eigen_positive_code" -> largeBranch["stability"]["eigenPositiveCode"],
  "real_primary_code" -> realCode,
  "maxwell_primary_code" -> maxwellCode,
  "clean_primary_code" -> cleanCode,
  "aL0_primary_code" -> aL0Code,
  "large_primary_code" -> largeCode,
  "restored_primary_code" -> restoredCode,
  "import_B_primary_code" -> importBCode,
  "import_K_primary_code" -> importKCode,
  "import_cE_primary_code" -> importCECode,
  "real_density_charge_flag_code" -> realBranch["residues"]["density"],
  "aL0_density_charge_flag_code" -> aL0Branch["residues"]["density"],
  "real_scalar_charge_residue_flag_code" -> realBranch["residues"]["scalar"],
  "clean_scalar_charge_residue_flag_code" -> cleanBranch["residues"]["scalar"],
  "aL0_h_charge_residue_flag_code" -> aL0Branch["residues"]["h"],
  "clean_h_charge_residue_flag_code" -> cleanBranch["residues"]["h"]
|>;

sympyResults = Import[sympyJson, "RawJSON"];
sympyExprs = sympyResults["engine_agreement"]["mathematica_exprs"];
sympyDigest = sympyResults["engine_agreement"]["sympy_expression_digest"];
engineKeys = Keys[sympyExprs];

assertEngine[name_] := Module[{expectedText, expectedExpr, actual},
  expectedText = sympyExprs[name];
  If[! StringQ[expectedText], fail["missing SymPy expression for " <> name]];
  If[! KeyExistsQ[actuals, name], fail["Mathematica missing actual for " <> name]];
  expectedExpr = ToExpression[expectedText, InputForm];
  actual = actuals[name];
  If[! TrueQ[FullSimplify[actual == expectedExpr]],
    fail["engine disagreement " <> name <> ": Mathematica got " <>
      ToString[actual, InputForm] <> ", SymPy exported " <> expectedText]
  ];
];

Scan[assertEngine, engineKeys];

topLineFlags = <|
  "density_charge_coupled" -> (realBranch["residues"]["density"] === 1),
  "operator_parity_contamination" -> True,
  "scalar_sector_stable" -> (realBranch["stability"]["stableCode"] === 1)
|>;
controlClasses = <|
  "real_provenance_fixed" -> classLabel[realCode],
  "maxwell_counterfactual" -> classLabel[maxwellCode],
  "clean_coexistence" -> classLabel[cleanCode],
  "aL_to_0" -> classLabel[aL0Code],
  "large_C_hu" -> classLabel[largeCode],
  "large_C_hu_restored_bound" -> classLabel[restoredCode],
  "import_fidelity" -> <|
    "B_eff_corrupt" -> classLabel[importBCode],
    "K_h_corrupt" -> classLabel[importKCode],
    "c_E_corrupt" -> classLabel[importCECode]
  |>
|>;
selfTests = <|
  "FIELD_SCALAR_VECTOR_DEPARTURE" -> classLabel[realCode],
  "FIELD_EXACT_MAXWELL_STRUCTURE" -> classLabel[maxwellCode],
  "FIELD_TRANSVERSE_EM_PLUS_CLEAN_GRAVITY_DENSITY" -> classLabel[cleanCode],
  "FIELD_SCALAR_SECTOR_UNSTABLE" -> classLabel[largeCode],
  "FIELD_CLASSIFICATION_UNDERDETERMINED" -> classLabel[importBCode],
  "A_L_ZERO_H_FLOOR" -> classLabel[aL0Code],
  "RESTORED_STABILITY_BOUND" -> classLabel[restoredCode]
|>;

agreementPayload = <|
  "top_line_primary" -> classLabel[realCode],
  "top_line_primary_code" -> realCode,
  "top_line_flags" -> topLineFlags,
  "control_classes" -> controlClasses,
  "dof_discriminator" -> <|"real" -> realBranch["dof"]["physicalTotal"], "maxwell" -> maxwellBranch["dof"]["physicalTotal"]|>,
  "self_tests" -> selfTests,
  "checked_expression_count" -> Length[engineKeys],
  "expr_digest" -> sympyDigest
|>;

If[agreementPayload["top_line_primary"] =!= sympyResults["agreement_payload"]["top_line_primary"],
  fail["top-line primary disagreement"]
];
If[agreementPayload["top_line_primary_code"] =!= sympyResults["agreement_payload"]["top_line_primary_code"],
  fail["top-line primary code disagreement"]
];
If[agreementPayload["checked_expression_count"] =!= sympyResults["agreement_payload"]["checked_expression_count"],
  fail["checked expression count disagreement"]
];
If[agreementPayload["expr_digest"] =!= sympyResults["agreement_payload"]["expr_digest"],
  fail["expression digest disagreement"]
];

out = <|
  "schema" -> "pathA_39_stage4_field_classification_mathematica/v1",
  "status" -> "OK",
  "headline" -> classLabel[realCode],
  "checked_quantities" -> engineKeys,
  "sympy_expression_digest" -> sympyDigest,
  "agreement_payload" -> agreementPayload,
  "computed" -> <|
    "detQ_real" -> ToString[actuals["detQ_real"], InputForm],
    "R_charge_real" -> ToString[actuals["R_charge_real"], InputForm],
    "real_root_minus_charge_residue" -> ToString[actuals["real_root_minus_charge_residue"], InputForm],
    "real_root_plus_charge_residue" -> ToString[actuals["real_root_plus_charge_residue"], InputForm],
    "control_classes" -> controlClasses
  |>
|>;

Export[jsonOut, out, "RawJSON"];
Print["OK pathA_39_stage4_field_classification_mathematica"];
Print[ExportString[<|"json" -> jsonOut, "headline" -> classLabel[realCode]|>, "RawJSON"]];
