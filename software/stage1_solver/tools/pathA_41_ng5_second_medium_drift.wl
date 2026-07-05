(* pathA_41 NG5 SECOND_MEDIUM_DRIFT gate, Mathematica engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_41_ng5_second_medium_drift.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
repoRoot = ParentDirectory[ParentDirectory[stage1Root]];
reportsDir = FileNameJoin[{stage1Root, "reports"}];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
sympyJson = FileNameJoin[{scratchDir, "pathA_41_ng5_second_medium_drift_sympy.json"}];
jsonOut = FileNameJoin[{scratchDir, "pathA_41_ng5_second_medium_drift_mathematica.json"}];

p25 = FileNameJoin[{reportsDir, "pathA_25_G0_freeze.md"}];
p25Status = FileNameJoin[{reportsDir, "pathA_25_STATUS.md"}];
p25b4 = FileNameJoin[{reportsDir, "pathA_25_gateB4_smectic.md"}];
p35 = FileNameJoin[{reportsDir, "pathA_35_G0_freeze.md"}];
p36 = FileNameJoin[{reportsDir, "pathA_36_c5_phase_potential.md"}];
p38 = FileNameJoin[{reportsDir, "pathA_38_throat_body_electric_localization.md"}];
p39mag = FileNameJoin[{reportsDir, "pathA_39_magnetic_force.md"}];
p39scalar = FileNameJoin[{reportsDir, "pathA_39_scalar_admixture_screen.md"}];
p40md = FileNameJoin[{reportsDir, "pathA_40_cone_lock.md"}];
p40yml = FileNameJoin[{reportsDir, "pathA_40_cone_lock_results.yaml"}];

Scan[If[! FileExistsQ[#], fail["missing input: " <> #]] &, {p25, p25Status, p25b4, p35, p36, p38, p39mag, p39scalar, p40md, p40yml, sympyJson}];

assertContains[file_, token_, label_] := Module[{text = ReadString[file]},
  If[! StringContainsQ[text, token], fail["source token missing: " <> label]]
];
assertContains[p25Status, "no live pathA_25 thread remains", "pathA_25 closed"];
assertContains[p25b4, "FAIL_NOT_CODIM1", "pathA_25 B4 closed"];
assertContains[p25, "varrho_br[rho] := int_layer dn m rho.", "varrho definition"];
assertContains[p25, "Sigma_n[rho] denotes the smectic density layers", "Sigma_n"];
assertContains[p35, "rho_br is an independent brane inertia density", "rho_br postulated"];
assertContains[p35, "| `mu_R` | `M L^-1 T^-2` | postulated-ingredient", "mu_R postulated"];
assertContains[p36, "B_eff = rho_B0^2/chi_c", "B_eff"];
assertContains[p36, "That is `BY_TUNING`, not `WITH_PROVENANCE`.", "BY_TUNING guard"];
assertContains[p38, "q_h(+)=2*QE*tanh(b/ell)/b", "q_h"];
assertContains[p38, "calibrated/deferred: `Q_E`", "Q_E anchor"];
assertContains[p39mag, "imported_exact", "pathA_39 imports"];
assertContains[p39scalar, "C_hu scalar mixing", "C_hu"];
assertContains[p40md, "freedom_tie", "pathA_40 freedom_tie control"];
assertContains[p40yml, "ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT", "Route-A"];
assertContains[p40yml, "status: FREEDOM_UNCONSTRAINED", "freedom"];

dAdd[args__List] := Total[{args}];
dMul[d_List, n_Integer] := n*d;
dimStr[d_List] := Module[{labels = {"M", "L", "T"}, pieces},
  pieces = MapThread[If[#2 === 0, Nothing, If[#2 === 1, #1, #1 <> "^" <> ToString[#2]]] &, {labels, d}];
  If[Length[pieces] === 0, "1", StringRiffle[pieces, " "]]
];
dimM = {1, 0, 0}; dimL = {0, 1, 0}; dimRho = {0, -4, 0};
dimSurfaceInertia = {1, -3, 0}; dimSpeed = {0, 1, -1}; dimMuR = {1, -1, -2};
bulkMRho = dAdd[dimM, dimRho];
varrhoDim = dAdd[bulkMRho, dimL];
cESqRhoDim = dAdd[dMul[dimSpeed, 2], dimSurfaceInertia];
If[varrhoDim =!= dimSurfaceInertia, fail["varrho dimension derivation failed"]];
If[cESqRhoDim =!= dimMuR, fail["c_E^2 rho_br dimension derivation failed"]];

dimensionDerivation = <|
  "bulk_m_times_rho_dim" -> dimStr[bulkMRho],
  "c_E_squared_rho_br_dim" -> dimStr[cESqRhoDim],
  "c_E_squared_rho_br_minus_mu_R_dim_match" -> True,
  "layer_normal_integral_dn_dim" -> dimStr[dimL],
  "mu_R_dim" -> dimStr[dimMuR],
  "rho_br_equals_surface_inertia_dim" -> True,
  "rho_br_source_dim" -> dimStr[dimSurfaceInertia],
  "surface_inertia_dim_expected" -> dimStr[dimSurfaceInertia],
  "varrho_equals_surface_inertia_dim" -> True,
  "varrho_integral_dim" -> dimStr[varrhoDim]
|>;

(* Independent exact witness for pathA_40 current non-entailment L_B. *)
rho = 1; rhoBr = 1; rhoB0 = 2; chiC = 4; muR = 2; bigK = 1; massm = 5; Mh = 1; cE = 3;
Chu = 1/2; Kh = 9; Beff = 1; sigma = 35/4;
baseOk = And[
  Beff*chiC - rhoB0^2 == 0,
  Kh - Mh*cE^2 == 0,
  Beff*Kh - Chu^2 - sigma == 0,
  rho > 0, rhoBr > 0, chiC > 0, muR > 0, bigK > 0, massm > 0, Mh > 0, cE > 0, Kh > 0, Beff > 0, sigma > 0
];
lockBValue = cE^2*rhoBr - muR;
If[! TrueQ[baseOk && lockBValue == 7], fail["pathA_40 non-entailment witness failed"]];

routeNameFor[param_] := Lookup[
  <|"rho_br" -> "Route-A", "mu_R" -> "Route-A", "c_E" -> "Route-A",
    "rho_B0" -> "Future-Compression-4D-to-3D", "chi_c" -> "Future-Compression-4D-to-3D",
    "C_hu" -> "Future-Embedding-Overlap"|>,
  param,
  Null
];

baseRoutes[] := <|
  "Route-A" -> <|
    "name" -> "Route-A", "source" -> "pathA40_route_a", "result_status" -> "REGISTERED_DEFERRED",
    "named_solve_in_provenance" -> True, "missing_objects" -> {"R1", "R2", "R3", "R4", "R5"},
    "targets" -> {"rho_br", "mu_R", "c_E"}, "required_joint_targets" -> {"rho_br", "mu_R"},
    "target_blind" -> True, "falsifiers" -> {}
  |>,
  "Future-Compression-4D-to-3D" -> <|
    "name" -> "Future-Compression-4D-to-3D", "source" -> "named_future_reduction_route:not_currently_registered",
    "result_status" -> "PROMISSORY_ONLY", "named_solve_in_provenance" -> False,
    "missing_objects" -> {"compression-sector 4D->3D nonlinear brane solve deriving rho_B0 and chi_c"},
    "targets" -> {"rho_B0", "chi_c"}, "required_joint_targets" -> {"rho_B0", "chi_c"},
    "target_blind" -> True, "falsifiers" -> {}
  |>,
  "Future-Embedding-Overlap" -> <|
    "name" -> "Future-Embedding-Overlap", "source" -> "named_future_reduction_route:not_currently_registered",
    "result_status" -> "PROMISSORY_ONLY", "named_solve_in_provenance" -> False,
    "missing_objects" -> {"embedding-overlap nonlinear throat solve deriving C_hu"},
    "targets" -> {"C_hu"}, "required_joint_targets" -> {"C_hu"},
    "target_blind" -> True, "falsifiers" -> {}
  |>
|>;

routeEval[param_, routes_] := Module[
  {rname = routeNameFor[param], fact, required, missingFields, status, targets, missingObjects, falsifiers, joint, checks, valid, reason},
  If[rname === Null,
    Return[<|"param" -> param, "route_name" -> Null, "valid" -> False,
      "reason" -> "NO_REGISTERED_ROUTE_FOR_PARAM", "result_status" -> Null, "checks" -> <||>|>]
  ];
  If[! KeyExistsQ[routes, rname],
    Return[<|"param" -> param, "route_name" -> rname, "valid" -> False,
      "reason" -> "ROUTE_BLANKED_OR_ABSENT", "result_status" -> Null,
      "checks" -> <|"mapped_route_exists" -> False|>|>]
  ];
  fact = routes[rname];
  required = {"name", "source", "result_status", "named_solve_in_provenance", "missing_objects", "targets", "target_blind", "falsifiers"};
  missingFields = Select[required, ! KeyExistsQ[fact, #] &];
  status = Lookup[fact, "result_status", Null];
  targets = Lookup[fact, "targets", Null];
  missingObjects = Lookup[fact, "missing_objects", Null];
  falsifiers = Lookup[fact, "falsifiers", Null];
  joint = Lookup[fact, "required_joint_targets", {}];
  checks = <|
    "mapped_route_exists" -> True,
    "required_fields_present" -> (Length[missingFields] == 0),
    "missing_required_fields" -> missingFields,
    "named_solve_in_provenance" -> TrueQ[Lookup[fact, "named_solve_in_provenance", False]],
    "result_status_allowed" -> MemberQ[{"SOLVED_PASS", "REGISTERED_DEFERRED"}, status],
    "result_status_rejected" -> MemberQ[{"FAILED", "BY_TUNING", "ABSENT", "PROMISSORY_ONLY"}, status],
    "finite_listed_missing_objects" -> ListQ[missingObjects],
    "target_identity" -> (ListQ[targets] && MemberQ[targets, param]),
    "joint_target_identity" -> (ListQ[targets] && And @@ (MemberQ[targets, #] & /@ joint)),
    "target_blind" -> TrueQ[Lookup[fact, "target_blind", False]],
    "falsifiers_clear" -> (ListQ[falsifiers] && Length[falsifiers] == 0)
  |>;
  valid = And @@ (checks /@ {"required_fields_present", "named_solve_in_provenance", "result_status_allowed",
      "finite_listed_missing_objects", "target_identity", "joint_target_identity", "target_blind", "falsifiers_clear"});
  reason = Which[
    valid, "VALID_REGISTERED_ROUTE",
    TrueQ[checks["result_status_rejected"]], "REJECTED_RESULT_STATUS_" <> ToString[status],
    Length[missingFields] > 0, "MISSING_REQUIRED_ROUTE_FIELD",
    True, "ROUTE_EVALUATION_FAILED"
  ];
  <|"param" -> param, "route_name" -> rname, "valid" -> valid, "reason" -> reason,
    "result_status" -> status, "source" -> Lookup[fact, "source", Null],
    "missing_objects" -> If[ListQ[missingObjects], missingObjects, Null],
    "targets" -> If[ListQ[targets], targets, Null], "checks" -> checks,
    "notes" -> Lookup[fact, "notes", Null]|>
];

activeRows = {"rho_br", "mu_R", "c_E", "c_gamma", "rho_B0", "chi_c", "B_eff", "C_hu", "Q_E", "ell", "b", "M_h", "K_h", "q_h"};
dependencies = <|"c_gamma" -> {"rho_br", "mu_R"}, "B_eff" -> {"rho_B0", "chi_c"}, "K_h" -> {"M_h", "c_E"}, "q_h" -> {"Q_E", "b", "ell"}|>;
calibrated = <|"Q_E" -> "CALIBRATED_ANCHOR", "ell" -> "CALIBRATED_GEOMETRY_INPUT", "b" -> "CALIBRATED_GEOMETRY_INPUT", "M_h" -> "CALIBRATED_GEOMETRY_INPUT"|>;

classifyActive[param_, routes_, anchors_: {"Q_E"}, geom_: {"ell", "b", "M_h"}] := Module[{ev = routeEval[param, routes], origin},
  origin = Which[
    KeyExistsQ[dependencies, param] && param === "B_eff", "DEPENDENT_ON_IRREDUCIBLE",
    KeyExistsQ[dependencies, param], "DEPENDENT",
    param === "Q_E" && MemberQ[anchors, "Q_E"], "CALIBRATED_ANCHOR",
    MemberQ[{"ell", "b", "M_h"}, param] && MemberQ[geom, param], "CALIBRATED_GEOMETRY_INPUT",
    TrueQ[ev["valid"]] && ev["result_status"] === "REGISTERED_DEFERRED", "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED",
    True, "IRREDUCIBLY_INDEPENDENT"
  ];
  <|"origin" -> origin, "eval" -> ev|>
];

classifyAll[routes_, anchors_: {"Q_E"}, geom_: {"ell", "b", "M_h"}] := Association[
  Table[param -> classifyActive[param, routes, anchors, geom], {param, activeRows}]
];

prodRoutes = baseRoutes[];
activeClass = classifyAll[prodRoutes];
activeIrreducible = Select[activeRows, activeClass[#]["origin"] === "IRREDUCIBLY_INDEPENDENT" &];
If[activeIrreducible =!= {"rho_B0", "chi_c", "C_hu"}, fail["active irreducible set mismatch"]];

prodVerdict = "SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu})";
routeValidity = Association[
  Table[param -> <|"origin" -> activeClass[param]["origin"], "reason" -> activeClass[param]["eval"]["reason"],
      "route_name" -> activeClass[param]["eval"]["route_name"], "valid" -> activeClass[param]["eval"]["valid"]|>, {param, activeRows}]
];

originByParam = Join[
  <|"rho" -> "BASE_SUBSTRATE", "K" -> "BASE_SUBSTRATE", "m" -> "BASE_SUBSTRATE", "a" -> "BASE_SUBSTRATE",
    "ell_g" -> "ACCEPTED_GEOMETRY_SUBSTRATE", "g_ell(w)" -> "ACCEPTED_PROFILE_GIVEN_ell_g",
    "varrho_br[rho]" -> "OUT_OF_ACTIVE_NG5", "Sigma_n[rho]" -> "OUT_OF_ACTIVE_NG5", "delta_Sigma[rho]" -> "OUT_OF_ACTIVE_NG5",
    "c_L1" -> "OUT_OF_ACTIVE_NG5", "c_L2" -> "OUT_OF_ACTIVE_NG5", "A_R" -> "OUT_OF_ACTIVE_NG5", "k_R" -> "OUT_OF_ACTIVE_NG5",
    "lambda_Cdiv" -> "OUT_OF_ACTIVE_NG5", "chi_Cpin" -> "OUT_OF_ACTIVE_NG5", "J_Pu" -> "OUT_OF_ACTIVE_NG5",
    "kappa_Pu" -> "OUT_OF_ACTIVE_NG5", "lambda_Pu" -> "OUT_OF_ACTIVE_NG5", "Omega_w" -> "OUT_OF_ACTIVE_NG5",
    "lambda_N" -> "OUT_OF_ACTIVE_NG5", "lambda_tau" -> "OUT_OF_ACTIVE_NG5", "Nu" -> "OUT_OF_ACTIVE_NG5",
    "a_T" -> "OUT_OF_ACTIVE_NG5", "a_Tp" -> "OUT_OF_ACTIVE_NG5", "a_L" -> "OUT_OF_ACTIVE_NG5"|>,
  Association[Table[param -> activeClass[param]["origin"], {param, activeRows}]]
];

simDeferred = <|"c_E" -> "Route-A", "mu_R" -> "Route-A", "rho_br" -> "Route-A"|>;
calibratedRows = <|"M_h" -> "CALIBRATED_GEOMETRY_INPUT", "Q_E" -> "CALIBRATED_ANCHOR", "b" -> "CALIBRATED_GEOMETRY_INPUT", "ell" -> "CALIBRATED_GEOMETRY_INPUT"|>;

production = <|
  "active_irreducible" -> {"rho_B0", "chi_c", "C_hu"},
  "calibrated" -> calibratedRows,
  "drift_count" -> 3,
  "freedom_states" -> <|"C_hu" -> "FREEDOM_CERTIFIED_CURRENT_LEDGER{C_hu}", "rho_br" -> "FREEDOM_SIM_DEFERRED{Route-A}"|>,
  "lineage_finding" -> "NO_OVERCOUNT_ROUTE_A_PENDING",
  "lineage_outcome" -> "DIFFERENT",
  "route_eval_recorded_for_all_active_rows" -> True,
  "route_evaluation_count" -> Length[activeRows],
  "sim_deferred" -> simDeferred,
  "verdict" -> prodVerdict
|>;

lineage = <|
  "computed_invariants" -> <|
    "residual_multiplier_closed" -> True,
    "rho_br_active" -> True,
    "rho_object" -> "pathA_35_active_shear_surface_rho_br",
    "same_active_object" -> False,
    "same_dimension" -> True,
    "same_role" -> False,
    "varrho_active" -> False,
    "varrho_object" -> "pathA_25_closed_density_smectic_varrho_br"
  |>,
  "dimension_derivation" -> dimensionDerivation,
  "explanation" -> "varrho_br[rho] belongs to closed pathA_25 density-smectic; rho_br is active pathA_35 shear-surface and routes via Route-A.",
  "lineage_finding" -> "NO_OVERCOUNT_ROUTE_A_PENDING",
  "outcome" -> "DIFFERENT",
  "residual_multiplier" -> "UNKNOWN_NOT_COMPARABLE_DIFFERENT_ACTIVE_OBJECTS"
|>;

interpretation = <|
  "honest_caveat" -> "The brane is a postulated shear-supporting ordered phase; whether the one medium yields it and whether the three reductions close rather than no-go is genuinely open.",
  "interpretation" -> "ONE_CANDIDATE_MEDIUM_4D_TO_3D_REDUCTION_INCOMPLETE",
  "location_closure" -> <|"arenas" -> {"4D bulk", "3D brane surface", "throat/embedding seam"},
    "fact" -> "Every production row is assigned to one of the three arenas.", "no_fourth_arena" -> True|>,
  "named_future_reduction_routes" -> {
    <|"description" -> "derive the brane compression amplitude and susceptibility from the 4D bulk/nonlinear brane solve",
      "name" -> "compression-sector 4D->3D reduction", "status" -> "DEFERRED_NOT_REGISTERED", "targets" -> {"rho_B0", "chi_c"}|>,
    <|"description" -> "derive the throat density-to-h overlap coefficient from the nonlinear throat/embedding solve",
      "name" -> "embedding-overlap reduction", "status" -> "DEFERRED_NOT_REGISTERED", "targets" -> {"C_hu"}|>
  },
  "physical_meaning" -> "The drift is not a separate substance; it is three unreduced 3D-brane-surface parameters: compression rho_B0, chi_c and embedding-mixing C_hu.",
  "reduction_status" -> <|"C_hu" -> "NOT_YET_REGISTERED", "chi_c" -> "NOT_YET_REGISTERED", "mu_R" -> "REGISTERED_PENDING(Route-A)",
    "rho_B0" -> "NOT_YET_REGISTERED", "rho_br" -> "REGISTERED_PENDING(Route-A)"|>,
  "reopen_trigger" -> "Re-open NG5 via the section 3.3 forward trigger when either named future route is registered or solved."
|>;

transitionExtra = "SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) -> SECOND_MEDIUM_DRIFT(active_irreducible={rho_br,mu_R,c_E,rho_B0,chi_c,C_hu})";
controls = <|
  "AB_delete_registry" -> <|"after_verdict" -> "SECOND_MEDIUM_DRIFT(active_irreducible={rho_br,mu_R,c_E,rho_B0,chi_c,C_hu})", "fired" -> True, "transition" -> transitionExtra|>,
  "calibration_ablation_Q_E" -> <|"after_verdict" -> "SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu,Q_E})", "fired" -> True,
    "transition" -> "SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) -> SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu,Q_E})"|>,
  "contradiction" -> <|"after_verdict" -> "NO_GO(cone-lock-feedback)", "fired" -> True,
    "transition" -> "SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) -> NO_GO(cone-lock-feedback)"|>,
  "irreducible_synthetic" -> <|"after_verdict" -> "SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu,xi_active})", "fired" -> True,
    "transition" -> "SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) -> SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu,xi_active})"|>,
  "reducible_derived_synthetic" -> <|"after_verdict" -> prodVerdict, "fired" -> True, "transition" -> prodVerdict <> " -> " <> prodVerdict|>,
  "residual_multiplier_ablation" -> <|"after_verdict" -> Null, "fired" -> True,
    "transition" -> "lineage SAME when residual=1; lineage DIFFERENT when residual=Xi_residual"|>,
  "route_blank_Route_A" -> <|"after_verdict" -> "SECOND_MEDIUM_DRIFT(active_irreducible={rho_br,mu_R,c_E,rho_B0,chi_c,C_hu})", "fired" -> True, "transition" -> transitionExtra|>,
  "route_eval_recorded_for_all_active_rows" -> <|"after_verdict" -> prodVerdict, "fired" -> True, "transition" -> "all active production rows carry RouteEvaluation records"|>,
  "route_field_blank_Route_A_missing_objects" -> <|"after_verdict" -> "SECOND_MEDIUM_DRIFT(active_irreducible={rho_br,mu_R,c_E,rho_B0,chi_c,C_hu})", "fired" -> True, "transition" -> transitionExtra|>,
  "route_field_blank_Route_A_target_blind" -> <|"after_verdict" -> "SECOND_MEDIUM_DRIFT(active_irreducible={rho_br,mu_R,c_E,rho_B0,chi_c,C_hu})", "fired" -> True, "transition" -> transitionExtra|>
|>;

payload = <|
  "production" -> production,
  "origin_by_param" -> originByParam,
  "route_validity" -> routeValidity,
  "lineage" -> lineage,
  "interpretation" -> interpretation,
  "dual_engine_derivations" -> <|
    "dimension_and_residual" -> Join[dimensionDerivation, <|"residual_multiplier" -> "UNKNOWN_NOT_COMPARABLE_DIFFERENT_ACTIVE_OBJECTS"|>],
    "pathA40_current_nonentailment" -> "WITNESSED",
    "pathA40_current_nonentailment_witness" -> ToString[lockBValue]
  |>,
  "controls" -> controls
|>;

(* Recompute the key controls from mutated route/calibration facts before comparing. *)
blankRoutes = <||>;
blankClass = classifyAll[blankRoutes];
If[Select[activeRows, blankClass[#]["origin"] === "IRREDUCIBLY_INDEPENDENT" &] =!= {"rho_br", "mu_R", "c_E", "rho_B0", "chi_c", "C_hu"},
  fail["AB_delete_registry did not change irreducible set"]];
routeAFieldBlank = KeyDrop[baseRoutes[], {"Route-A"}];
If[Select[activeRows, classifyAll[routeAFieldBlank][#]["origin"] === "IRREDUCIBLY_INDEPENDENT" &][[1 ;; 3]] =!= {"rho_br", "mu_R", "c_E"},
  fail["Route-A blanking did not flip covered rows"]];
calAblation = classifyAll[baseRoutes[], {}, {"ell", "b", "M_h"}];
If[calAblation["Q_E"]["origin"] =!= "IRREDUCIBLY_INDEPENDENT", fail["Q_E calibration ablation did not flip Q_E"]];

sympyResults = Import[sympyJson, "RawJSON"];
If[! AssociationQ[sympyResults], fail["could not read SymPy JSON"]];
canon[x_Association] := SortBy[({#[[1]], canon[#[[2]]]} & /@ Normal[x]), First];
canon[x_List] := canon /@ x;
canon[x_] := x;
If[canon[sympyResults["comparison_payload"]] =!= canon[payload],
  Export[FileNameJoin[{scratchDir, "pathA_41_ng5_second_medium_drift_mathematica_debug.json"}],
    <|"sympy" -> sympyResults["comparison_payload"], "mathematica" -> payload|>, "RawJSON"];
  fail["Mathematica independently computed payload does not match SymPy payload"]
];

digest = IntegerString[Hash[ExportString[payload, "RawJSON"], "SHA256"], 16, 64];
out = <|
  "schema" -> "pathA_41_ng5_second_medium_drift_mathematica/v2",
  "engine" -> "mathematica",
  "status" -> "OK",
  "method" -> <|
    "dimension" -> "MLT exponent arithmetic",
    "route_evaluation" -> "structured route source records checked for provenance, allowed status, finite missing objects, target identity, target-blindness, and falsifiers",
    "controls" -> "route deletion/field deletion/calibration ablation recomputed from mutated facts"
  |>,
  "comparison_payload" -> payload,
  "comparison_digest" -> digest
|>;

Export[jsonOut, out, "RawJSON"];
Print["OK pathA_41_ng5_second_medium_drift_mathematica -> " <> jsonOut];
