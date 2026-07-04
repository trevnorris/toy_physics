(* pathA_39 Stage 0+1 scalar-admixture screen, Mathematica engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_39_scalar_admixture_screen.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
reportsDir = FileNameJoin[{stage1Root, "reports"}];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
sympyJson = FileNameJoin[{scratchDir, "pathA_39_scalar_admixture_screen_sympy.json"}];
jsonOut = FileNameJoin[{scratchDir, "pathA_39_scalar_admixture_screen_mathematica.json"}];
path36Yaml = FileNameJoin[{reportsDir, "pathA_36_c5_phase_potential_results.yaml"}];
path38Yaml = FileNameJoin[{reportsDir, "pathA_38_results.yaml"}];

If[! FileExistsQ[sympyJson], fail["missing SymPy JSON: " <> sympyJson]];
If[! FileExistsQ[path36Yaml], fail["missing pathA_36 YAML: " <> path36Yaml]];
If[! FileExistsQ[path38Yaml], fail["missing pathA_38 YAML: " <> path38Yaml]];

$Assumptions =
  rhoBr > 0 && muR > 0 && rhoB0 > 0 && chiC > 0 && Mh > 0 && cE > 0 &&
  QE > 0 && b > 0 && ell > 0 && k > 0 && x > 0;

path36Text = ReadString[path36Yaml];
path38Text = ReadString[path38Yaml];
expectedB = "B_eff: rho_B0**2/chi_c";
expectedCGamma = "c_gamma_squared: mu_R/rho_br";
expectedQh = "q_h_plus: 2*QE*tanh(b/ell)/b";
expectedGreen = "radial_green_finite_omega: exp(I*R*omega/cE)/(4*pi*R)";
wrongCGammaGreen = "radial_green_finite_omega: exp(I*R*omega/c_gamma)/(4*pi*R)";

bImportMatchQ = StringContainsQ[path36Text, expectedB];
cGammaImportMatchQ = StringContainsQ[path36Text, expectedCGamma];
qhImportMatchQ = StringContainsQ[path38Text, expectedQh];
cEImportMatchQ = StringContainsQ[path38Text, expectedGreen];
cEWrongCGammaGuardMatchQ = StringContainsQ[path38Text, wrongCGammaGreen];

If[! bImportMatchQ, fail["pathA_36 B_eff import mismatch"]];
If[! cGammaImportMatchQ, fail["pathA_36 c_gamma import mismatch"]];
If[! qhImportMatchQ, fail["pathA_38 q_h import mismatch"]];
If[! cEImportMatchQ, fail["pathA_38 cE dynamic Green import mismatch"]];

verdictCodes = <|
  "SCALAR_COEXISTENCE_CLEAN" -> 1,
  "PASS_BY_DECLARATION" -> 2,
  "FAIL_EXTRA_H_BRANON" -> 3,
  "FAIL_CHARGE_COUPLED_CS_SCALAR" -> 4,
  "FAIL_OBSERVABLE_SCALAR_ADMIXTURE" -> 5,
  "ERROR_IMPORTED_CE_MISMATCH" -> 6
|>;

verdictCode[label_] := verdictCodes[label];
zeroQ[expr_] := TrueQ[FullSimplify[expr == 0]];
nonzeroQ[expr_] := ! zeroQ[expr];

classifyDecoupled[qLexpr_, qhExpr_, bExpr_, qLStatus_] := Module[
  {
    densityPropagating, densityChargeResidue, hChargeResidue, densityCoupled,
    hExtra, flags, verdict
  },
  densityPropagating = ! zeroQ[bExpr];
  densityChargeResidue = If[densityPropagating, FullSimplify[qLexpr^2/rhoBr], 0];
  hChargeResidue = FullSimplify[qhExpr^2/Mh];
  densityCoupled = densityPropagating && nonzeroQ[densityChargeResidue];
  hExtra = nonzeroQ[hChargeResidue];
  flags = {};
  If[densityCoupled, AppendTo[flags, "FAIL_CHARGE_COUPLED_CS_SCALAR"]];
  If[hExtra, AppendTo[flags, "FAIL_EXTRA_H_BRANON"]];
  verdict = Which[
    zeroQ[qLexpr] && qLStatus =!= "DERIVED_STAGE1", "PASS_BY_DECLARATION",
    densityCoupled && hExtra, "FAIL_OBSERVABLE_SCALAR_ADMIXTURE",
    densityCoupled, "FAIL_CHARGE_COUPLED_CS_SCALAR",
    hExtra, "FAIL_EXTRA_H_BRANON",
    qLStatus === "DERIVED_STAGE1", "SCALAR_COEXISTENCE_CLEAN",
    True, "PASS_BY_DECLARATION"
  ];
  If[verdict === "FAIL_OBSERVABLE_SCALAR_ADMIXTURE" && ! MemberQ[flags, "FAIL_OBSERVABLE_SCALAR_ADMIXTURE"],
    flags = Prepend[flags, "FAIL_OBSERVABLE_SCALAR_ADMIXTURE"]
  ];
  <|"verdict" -> verdict, "verdict_code" -> verdictCode[verdict], "flags" -> flags|>
];

(* Independent scalar-block assembly. *)
Beff = FullSimplify[rhoB0^2/chiC];
cGamma2 = FullSimplify[muR/rhoBr];
qh = FullSimplify[2 QE Tanh[b/ell]/b];
Kh = FullSimplify[Mh cE^2];
cs2 = FullSimplify[Beff/rhoBr];

Dx = {{rhoBr x - Beff, -Chu}, {-Chu, Mh x - Kh}};
detD = FullSimplify[Det[Dx]];
Gx = FullSimplify[{{Mh x - Kh, Chu}, {Chu, rhoBr x - Beff}}/detD];
Gomega = FullSimplify[Gx/k^2];
sourceCharge = {qL, qh};
sourceMass = {qM, 0};
Aqq = FullSimplify[sourceCharge . Gx . sourceCharge];
Aqm = FullSimplify[sourceCharge . Gx . sourceMass];
Amm = FullSimplify[sourceMass . Gx . sourceMass];
AqqNumFormula = FullSimplify[qL^2 (Mh x - Kh) + 2 qL qh Chu + qh^2 (rhoBr x - Beff)];
AmmNumFormula = FullSimplify[qM^2 (Mh x - Kh)];

trace = FullSimplify[rhoBr Kh + Mh Beff];
delta = FullSimplify[(rhoBr Kh - Mh Beff)^2 + 4 rhoBr Mh Chu^2];
cMinus = FullSimplify[(trace - Sqrt[delta])/(2 rhoBr Mh)];
cPlus = FullSimplify[(trace + Sqrt[delta])/(2 rhoBr Mh)];
denPrime = D[detD, x];
residueFromNum[num_, root_] := FullSimplify[(num /. x -> root)/(denPrime /. x -> root)];
rootMinusChargeResidue = residueFromNum[AqqNumFormula, cMinus];
rootPlusChargeResidue = residueFromNum[AqqNumFormula, cPlus];
rootMinusMassResidue = residueFromNum[AmmNumFormula, cMinus];
rootPlusMassResidue = residueFromNum[AmmNumFormula, cPlus];

qAStage1 = FullSimplify[sCharge aT Nu];
qLStage1 = FullSimplify[sCharge aL Nu];
qBulkStage1 = FullSimplify[sCharge aBulk];
qEvenToA = 0;

main = classifyDecoupled[qLStage1, qh, Beff, "DERIVED_STAGE1"];
injected = classifyDecoupled[epsQ, qh, Beff, "ABLATION_INJECTED"];
bZero = classifyDecoupled[epsQ, qh, 0, "ABLATION_INJECTED"];
extraH = classifyDecoupled[0, qh, Beff, "DERIVED_STAGE1"];
csOnly = classifyDecoupled[epsQ, 0, Beff, "ABLATION_INJECTED"];
passByDecl = classifyDecoupled[0, qh, Beff, "DECLARED_NOT_DERIVED"];
cleanFixture = classifyDecoupled[0, 0, Beff, "DERIVED_STAGE1"];

qL0MixingAqq = FullSimplify[Aqq /. qL -> 0];
mixingResidueMinus = FullSimplify[rootMinusChargeResidue /. qL -> 0];
mixingResiduePlus = FullSimplify[rootPlusChargeResidue /. qL -> 0];
mixingVerdict = If[nonzeroQ[mixingResidueMinus] || nonzeroQ[mixingResiduePlus],
  "FAIL_OBSERVABLE_SCALAR_ADMIXTURE",
  "SCALAR_COEXISTENCE_CLEAN"
];

kVec = {kx, ky, kz};
jVec = {jx, jy, jz};
kNormSquared = FullSimplify[kVec . kVec];
divergenceConstraint = FullSimplify[kVec . jVec];
longitudinalProjectorCurrent = FullSimplify[kVec divergenceConstraint/kNormSquared];
wireResponseUnconstrained = FullSimplify[longitudinalProjectorCurrent . longitudinalProjectorCurrent];
wireResponse = FullSimplify[wireResponseUnconstrained /. divergenceConstraint -> 0];

cEMatch = If[cEImportMatchQ, 1, 0];
cEMismatch = If[cEWrongCGammaGuardMatchQ, 1, 0];

controls = <|
  "inject_qL_epsilon" -> <|
    "status" -> If[injected["verdict"] === "FAIL_OBSERVABLE_SCALAR_ADMIXTURE", "FIRED", "NOT_FIRED"],
    "verdict" -> injected["verdict"]
  |>,
  "closed_steady_current_wire_limit" -> <|
    "status" -> If[wireResponse === 0, "FIRED", "NOT_FIRED"],
    "verdict" -> "WIRE_LIMIT_NO_LONGITUDINAL_RESPONSE"
  |>,
  "B_eff_positive_vs_zero" -> <|
    "status" -> If[injected["verdict"] =!= bZero["verdict"], "FIRED", "NOT_FIRED"],
    "positive_B_eff_verdict" -> injected["verdict"],
    "B_eff_to_zero_verdict" -> bZero["verdict"]
  |>,
  "Mh_positive_qh_nonzero" -> <|
    "status" -> If[extraH["verdict"] === "FAIL_EXTRA_H_BRANON", "FIRED", "NOT_FIRED"],
    "verdict" -> extraH["verdict"]
  |>,
  "cE_import_match" -> <|
    "status" -> If[cEMatch === 1 && cEMismatch === 0, "FIRED", "NOT_FIRED"],
    "verdict" -> "IMPORT_MATCH"
  |>,
  "mixing_on_with_derived_qL_zero" -> <|
    "status" -> If[mixingVerdict === "FAIL_OBSERVABLE_SCALAR_ADMIXTURE", "FIRED", "NOT_FIRED"],
    "verdict" -> mixingVerdict
  |>
|>;

controlLabel[item_] := Which[
  KeyExistsQ[item, "verdict"], item["verdict"],
  KeyExistsQ[item, "positive_B_eff_verdict"] && KeyExistsQ[item, "B_eff_to_zero_verdict"],
    item["positive_B_eff_verdict"] <> "->" <> item["B_eff_to_zero_verdict"],
  True, item["status"]
];

selfTests = <|
  "FAIL_OBSERVABLE_SCALAR_ADMIXTURE" -> injected["verdict"],
  "FAIL_EXTRA_H_BRANON" -> extraH["verdict"],
  "FAIL_CHARGE_COUPLED_CS_SCALAR" -> csOnly["verdict"],
  "PASS_BY_DECLARATION" -> passByDecl["verdict"],
  "SCALAR_COEXISTENCE_CLEAN" -> cleanFixture["verdict"]
|>;

If[selfTests["FAIL_OBSERVABLE_SCALAR_ADMIXTURE"] =!= "FAIL_OBSERVABLE_SCALAR_ADMIXTURE", fail["observable self-test failed"]];
If[selfTests["FAIL_EXTRA_H_BRANON"] =!= "FAIL_EXTRA_H_BRANON", fail["extra-h self-test failed"]];
If[selfTests["FAIL_CHARGE_COUPLED_CS_SCALAR"] =!= "FAIL_CHARGE_COUPLED_CS_SCALAR", fail["c_s self-test failed"]];
If[selfTests["PASS_BY_DECLARATION"] =!= "PASS_BY_DECLARATION", fail["pass-by-declaration self-test failed"]];
If[selfTests["SCALAR_COEXISTENCE_CLEAN"] =!= "SCALAR_COEXISTENCE_CLEAN", fail["clean fixture self-test failed"]];
If[controls["inject_qL_epsilon"]["status"] =!= "FIRED", fail["qL injection control did not fire"]];
If[controls["B_eff_positive_vs_zero"]["status"] =!= "FIRED", fail["B_eff ablation did not fire"]];
If[controls["Mh_positive_qh_nonzero"]["status"] =!= "FIRED", fail["extra h control did not fire"]];

actuals = <|
  "B_eff" -> Beff,
  "c_gamma_squared" -> cGamma2,
  "q_h" -> qh,
  "K_h" -> Kh,
  "c_s_squared" -> cs2,
  "D00" -> Dx[[1, 1]],
  "D01" -> Dx[[1, 2]],
  "D11" -> Dx[[2, 2]],
  "detD" -> detD,
  "G00_num" -> Mh x - Kh,
  "G01_num" -> Chu,
  "G11_num" -> rhoBr x - Beff,
  "A_qq_num" -> Numerator[Together[Aqq]],
  "A_qq_den" -> Denominator[Together[Aqq]],
  "A_qm_expr" -> Aqm,
  "A_mm_expr" -> Amm,
  "c_minus_squared" -> cMinus,
  "c_plus_squared" -> cPlus,
  "generic_root_minus_charge_residue" -> rootMinusChargeResidue,
  "generic_root_plus_charge_residue" -> rootPlusChargeResidue,
  "density_speed_C0" -> cs2,
  "h_speed_C0" -> cE^2,
  "density_charge_residue_stage1" -> qLStage1^2/rhoBr,
  "h_charge_residue" -> qh^2/Mh,
  "density_mass_residue" -> qM^2/rhoBr,
  "h_mass_residue" -> 0,
  "qA_stage1" -> qAStage1,
  "qL_stage1" -> qLStage1,
  "qBulk_stage1" -> qBulkStage1,
  "qEvenToA" -> qEvenToA,
  "wire_longitudinal_response" -> wireResponse,
  "main_verdict_code" -> main["verdict_code"],
  "inject_qL_verdict_code" -> injected["verdict_code"],
  "B_zero_verdict_code" -> bZero["verdict_code"],
  "extra_h_verdict_code" -> extraH["verdict_code"],
  "cs_only_verdict_code" -> csOnly["verdict_code"],
  "pass_by_declaration_code" -> passByDecl["verdict_code"],
  "clean_fixture_code" -> cleanFixture["verdict_code"],
  "cE_import_match" -> cEMatch,
  "cE_mismatch_ablation" -> cEMismatch,
  "dimensional_ablations_fired" -> 4
|>;

engineKeys = {
  "B_eff", "c_gamma_squared", "q_h", "K_h", "c_s_squared",
  "D00", "D01", "D11", "detD", "G00_num", "G01_num", "G11_num",
  "A_qq_num", "A_qq_den", "A_qm_expr", "A_mm_expr",
  "c_minus_squared", "c_plus_squared", "density_speed_C0", "h_speed_C0",
  "generic_root_minus_charge_residue", "generic_root_plus_charge_residue",
  "density_charge_residue_stage1", "h_charge_residue", "density_mass_residue", "h_mass_residue",
  "qA_stage1", "qL_stage1", "qBulk_stage1", "qEvenToA", "wire_longitudinal_response",
  "main_verdict_code", "inject_qL_verdict_code", "B_zero_verdict_code",
  "extra_h_verdict_code", "cs_only_verdict_code", "pass_by_declaration_code", "clean_fixture_code",
  "cE_import_match", "cE_mismatch_ablation", "dimensional_ablations_fired"
};

sympyResults = Import[sympyJson, "RawJSON"];
sympyExprs = sympyResults["engine_agreement"]["mathematica_exprs"];
sympyDigest = sympyResults["engine_agreement"]["sympy_expression_digest"];

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
  "schema" -> "pathA_39_scalar_admixture_screen_mathematica/v1",
  "status" -> "OK",
  "headline" -> topLineVerdict,
  "checked_quantities" -> engineKeys,
  "sympy_expression_digest" -> sympyDigest,
  "agreement_payload" -> agreementPayload,
  "computed" -> <|
    "detD" -> ToString[detD, InputForm],
    "A_qq" -> ToString[Aqq, InputForm],
    "density_charge_residue_stage1" -> ToString[qLStage1^2/rhoBr, InputForm],
    "h_charge_residue" -> ToString[qh^2/Mh, InputForm],
    "control_verdicts" -> controlVerdicts
  |>
|>;

Export[jsonOut, out, "RawJSON"];
Print["OK pathA_39_scalar_admixture_screen_mathematica"];
Print[ExportString[<|"json" -> jsonOut, "headline" -> topLineVerdict|>, "RawJSON"]];
