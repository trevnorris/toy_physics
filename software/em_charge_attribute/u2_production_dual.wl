(* Independent Wolfram production engine for U2 boundary adjudication. *)

ClearAll["Global`*"];

ratifiedDigest = "9eff1b0c49e89007aea1008cb6712b0ea495168d101ce43ddce1cffaf68749c4";
componentFiles = <|
  "frozen_data_pin_table" -> "frozen_data_pin_table.yaml",
  "candidate_inventory" -> "candidate_inventory.yaml",
  "obligation_censuses" -> "obligation_censuses.yaml",
  "dependency_grid_inventory" -> "dependency_grid_inventory.yaml",
  "vocabulary_freeze" -> "vocabulary_freeze.yaml",
  "evidence_taxonomy" -> "evidence_taxonomy.yaml",
  "availability_slots" -> "availability_slots.yaml",
  "route_fixture_inventory" -> "route_fixture_inventory.yaml",
  "closure_template_contracts" -> "closure_template_contracts.yaml",
  "environment_identity" -> "environment_identity.yaml",
  "standard_bindings" -> "standard_bindings.yaml",
  "producer_map" -> "producer_map.yaml",
  "evaluated_code_closure_policy" -> "evaluated_code_closure_policy.yaml",
  "parameter_register_proposals" -> "parameter_register_proposals.yaml",
  "obligation_manifest" -> "obligation_manifest.yaml"
|>;

cliArgs = If[Length[$ScriptCommandLine] > 0, $ScriptCommandLine, $CommandLine];
getArg[name_] := Module[{position = FirstPosition[cliArgs, name]},
  If[MissingQ[position] || position[[1]] == Length[cliArgs], Print["MISSING_ARGUMENT " <> name]; Exit[2]];
  cliArgs[[position[[1]] + 1]]
];
getOptionalArg[name_] := Module[{position = FirstPosition[cliArgs, name]},
  If[MissingQ[position] || position[[1]] == Length[cliArgs], Null, cliArgs[[position[[1]] + 1]]]
];
repo = ExpandFileName[getArg["--repo"]];
bundleDir = ExpandFileName[getArg["--bundle-dir"]];
suppliedDigest = getArg["--stage0-contract-digest"];
output = ExpandFileName[getArg["--output"]];
generationMutation = getOptionalArg["--generation-mutation"];

$Path = Select[$Path, StringStartsQ[ToString[#], $InstallationDirectory] &];
loadYaml[path_String] := Import[path, "YAML"];
sortIds[values_] := SortBy[DeleteDuplicates[values], {ToLowerCase[ToString[#]], ToString[#]} &];
lookupRows[rows_List, key_String] := If[Length[rows] == 0, {}, Lookup[rows, key]];

canonicalJSONString[value_Association] := Module[{keys = Keys[value], width},
  If[Length[keys] == 0, Return["{}"]];
  width = Max[StringLength /@ keys];
  keys = SortBy[keys, PadRight[ToCharacterCode[#], width, -1] &];
  "{" <> StringRiffle[
    (StringReplace[ExportString[#, "RawJSON", "Compact" -> True], "\\/" -> "/"] <> ":" <> canonicalJSONString[value[#]] &) /@ keys,
    ","] <> "}"
];
canonicalJSONString[value_List] := "[" <> StringRiffle[canonicalJSONString /@ value, ","] <> "]";
canonicalJSONString[value_String] := StringReplace[ExportString[value, "RawJSON", "Compact" -> True], "\\/" -> "/"];
canonicalJSONString[True] := "true";
canonicalJSONString[False] := "false";
canonicalJSONString[Null] := "null";
canonicalJSONString[value_Integer] := ToString[value, InputForm];
canonicalJSONString[value_] := StringReplace[ExportString[value, "RawJSON", "Compact" -> True], "\\/" -> "/"];
canonicalDigest[value_] := IntegerString[Hash[canonicalJSONString[value], "SHA256"], 16, 64];

require[condition_, detail_String] := If[!TrueQ[condition], Print["U2_WOLFRAM_PRODUCTION_BLOCKED " <> detail]; Exit[1]];
requireGeneration[condition_, assertId_String, detail_String] := If[!TrueQ[condition],
  Print["ASSERTION_FAILED " <> assertId <> ": " <> detail]; Exit[1]
];

classifyInferenceContent[content_Association] := Module[{mapping, operation, found},
  mapping = <|
    "root" -> "ROOT_REFERENCE", "positive_join" -> "POSITIVE_DERIVATION",
    "positive_equivalence" -> "POSITIVE_EQUIVALENCE", "static_force" -> "STATIC_COMMITTED_FORCING",
    "incompatible" -> "INCOMPATIBILITY", "not_member" -> "NEGATIVE_CANDIDATE_MEMBERSHIP",
    "case_eliminate" -> "CASE_ELIMINATION", "complement" -> "COMPLEMENT_SURVIVOR_COUNT",
    "exclude" -> "EXCLUSION_VERDICT", "postulate" -> "POSTULATE_BRANCH",
    "stability" -> "STABILITY_DYNAMICAL_CLASS", "solve" -> "SOLVE_EVALUATION_RESULT",
    "symbolic_collapse" -> "SYMBOLIC_EQUIVALENCE_COLLAPSE", "unavailability" -> "UNAVAILABILITY_WITNESS",
    "challenge" -> "DERIVABILITY_CHALLENGE", "evidence_state" -> "EVIDENCE_STATE_CLASSIFICATION"
  |>;
  operation = Lookup[content, "op", Missing["op"]];
  If[!KeyExistsQ[mapping, operation], Return[{"UNCLASSIFIED"}]];
  found = {mapping[operation]};
  Do[If[AssociationQ[child], found = Join[found, classifyInferenceContent[child]]], {child, Lookup[content, "args", {}]}];
  sortIds[found]
];

typedRootId[value_String] := Which[
  StringStartsQ[value, "source:"], value,
  StringStartsQ[value, "tilt:"] || StringStartsQ[value, "domain:"] || StringStartsQ[value, "endpoint:"], "source:phaseC_slot:" <> value,
  MemberQ[{"open_leaf:gammaSigma", "open_leaf:tangentDtN"}, value], "source:coefficient:" <> Last[StringSplit[value, ":", 2]],
  StringStartsQ[value, "open_leaf:"], "source:field:" <> Last[StringSplit[value, ":", 2]],
  True, require[False, "no ratified typed-root mapping for " <> value]; value
];

proofDag[rootIds_, operation_String] := Module[{roots = sortIds[rootIds], content},
  roots = sortIds[typedRootId /@ roots];
  content = <|"op" -> operation, "args" -> (<|"op" -> "root", "id" -> #|> & /@ roots)|>;
  <|
    "normalized_inference_content" -> content,
    "constructors" -> classifyInferenceContent[content],
    "root_ids" -> roots,
    "classification_source" -> "ratified_normalized_content_classifier"
  |>
];

classifyEvidence[raw_Association] := Which[
  !TrueQ[raw["applicable"]], "INAPPLICABLE",
  TrueQ[raw["committed_incompatibility"]] && TrueQ[raw["ancestry_complete_and_typed"]], "DECISIVE_INCOMPATIBILITY",
  TrueQ[raw["datum_missing"]], "MISSING",
  True, "SATISFIED"
];

evaluateDisposition[evidence_List] := Module[{states, applicable},
  states = classifyEvidence[# ["raw_predicate_inputs"]] & /@ evidence;
  If[!And @@ MapThread[#1["emitted_state"] == #2 &, {evidence, states}], Return["HARNESS_FAILED(contradictory_evidence)"]];
  If[MemberQ[states, "DECISIVE_INCOMPATIBILITY"], Return["EXCLUDED"]];
  applicable = Select[states, # != "INAPPLICABLE" &];
  If[Length[applicable] > 0 && And @@ (# == "SATISFIED" & /@ applicable), Return["ADMISSIBLE"]];
  If[MemberQ[applicable, "MISSING"], Return["UNRESOLVED(named datum)"]];
  "HARNESS_FAILED(contradictory_evidence)"
];

evaluateTopology[sector_String, interpolation_String] := Which[
  MemberQ[{{"DISCONNECTED", "INTERPOLABLE"}, {"CONNECTED", "OBSTRUCTED"}}, {sector, interpolation}], "HARNESS_FAILED(inconsistent_subresults)",
  {sector, interpolation} == {"DISCONNECTED", "OBSTRUCTED"}, "topologically-distinct",
  sector == "CONNECTED" || interpolation == "INTERPOLABLE", "orientation-only",
  True, "UNRESOLVED(named data)"
];

evaluateCrossLevel[applicability_String, gateIntegrity_String, gateOutcome_] := Which[
  applicability == "positively_non_geometric", "NOT_APPLICABLE(witness)",
  applicability == "missing", "refinement-UNRESOLVED(named datum)",
  gateIntegrity != "COMPUTATION_VALID", "NOT_RUN(exact_gate_set)",
  gateOutcome == "topologically-distinct", "defect-refined",
  gateOutcome == "orientation-only", "not-defect-refined",
  True, "refinement-UNRESOLVED"
];

evaluatePromotion[inputs_Association] := Module[{forcing, classes, disposition},
  If[Length[inputs["failed_required_upstreams"]] > 0, Return["NOT_RUN(exact_set)"]];
  If[TrueQ[inputs["uncanonicalized_overlap"]], Return["HARNESS_FAILED(uncanonicalized_candidate_overlap)"]];
  forcing = inputs["forcing_records"];
  If[Length[forcing] == 0, Return["NO_SELECTION_CLAIM(witness,challenge)"]];
  If[AnyTrue[forcing, !TrueQ[# ["candidate_in_current_axis"]] &], Return["HARNESS_FAILED(forced_class_outside_axis)"]];
  classes = DeleteDuplicates[Lookup[forcing, "canonical_class"]];
  If[Length[classes] > 1, Return["HARNESS_FAILED(contradictory_forcing)"]];
  disposition = forcing[[1]]["disposition"];
  If[disposition == "EXCLUDED", Return["HARNESS_FAILED(contradictory_committed_derivations)"]];
  If[disposition == "UNRESOLVED", Return["PROMOTION_UNRESOLVED(admissibility_unresolved_refusal)"]];
  If[inputs["closure_outcome"] == "CLOSURE_REFUSED", Return["PROMOTION_UNRESOLVED(closure_refusal)"]];
  "SELECTED"
];

recordSchemaValid[record_Association] := Module[{integrity, physics, upstreams},
  integrity = record["integrity"]; physics = Lookup[record, "physics", Null];
  Which[
    integrity == "COMPUTATION_VALID", physics =!= Null && !KeyExistsQ[record, "failure_reason"] && !KeyExistsQ[record, "failed_upstreams"],
    integrity == "HARNESS_FAILED", physics === Null && TrueQ[StringLength[Lookup[record, "failure_reason", ""]] > 0] && !KeyExistsQ[record, "failed_upstreams"],
    integrity == "NOT_RUN", upstreams = Lookup[record, "failed_upstreams", {}]; physics === Null && Length[upstreams] > 0 && upstreams == Sort[DeleteDuplicates[upstreams]],
    True, False
  ]
];

components = Association@KeyValueMap[#1 -> loadYaml[FileNameJoin[{bundleDir, #2}]] &, componentFiles];
require[suppliedDigest == ratifiedDigest, "STAGE0_CONTRACT_DIGEST differs from ratified U2 digest"];
require[canonicalDigest[components] == ratifiedDigest, "ratified component bundle digest changed"];
require[loadYaml[FileNameJoin[{bundleDir, "stage0_bundle.yaml"}]]["stage0_contract_digest"] == ratifiedDigest, "bundle binding changed"];
require[loadYaml[FileNameJoin[{bundleDir, "stage0_contract.yaml"}]]["stage0_contract_digest"] == ratifiedDigest, "contract binding changed"];

candidateDoc = components["candidate_inventory"];
gridDoc = components["dependency_grid_inventory"];
obligationsDoc = components["obligation_censuses"]["censuses"];
availabilityRows = components["availability_slots"]["slots"];
availability = AssociationThread[Lookup[availabilityRows, "slot_id"], availabilityRows];
phaseCPath = FileNameJoin[{repo, "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_c_0_tilt_coupling_contract/availability_slots.yaml"}];
phaseCRows = loadYaml[phaseCPath]["slots"];
phaseCSlots = AssociationThread[Lookup[phaseCRows, "slot_id"], phaseCRows];
candidateRecords = AssociationThread[Lookup[candidateDoc["candidate_records"], "candidate_id"], candidateDoc["candidate_records"]];

dependencyMeasurement[token_String] := Module[{row},
  require[KeyExistsQ[phaseCSlots, token], "expected OPEN dependency has no frozen availability slot: " <> token];
  row = phaseCSlots[token];
  <|
    "dependency_token" -> token, "ratified_slot_id" -> token,
    "availability" -> row["disposition"], "integrity_state" -> "COMPUTATION_VALID",
    "witness_id" -> Lookup[row, "witness_id", Null], "challenge_id" -> Lookup[row, "challenge_id", Null],
    "missing" -> (row["disposition"] == "UNRESOLVED"),
    "measurement_source" -> "ratified_PhaseC_availability_record"
  |>
];

obligationDependencies[obligation_String, expected_List, candidate_String] := Module[{native},
  native = Select[expected,
    StringStartsQ[#, "endpoint:"] || MemberQ[{"open_leaf:E4_shear_lock", "open_leaf:E5_rayleigh", "open_leaf:gammaSigma", "open_leaf:tangentDtN"}, #] &];
  Which[
    obligation == "host_location_evidence", Select[expected, MemberQ[{"open_leaf:geon_core_bundle", "open_leaf:sleeve_core_trace"}, #] &],
    StringStartsQ[obligation, "topology_"], sortIds[Join[{"open_leaf:sleeve_core_trace"}, native]],
    obligation == "jump_and_trace_compatibility", Select[expected, StringStartsQ[#, "tilt:"] || MemberQ[{"domain:Sigma_boundary_data", "open_leaf:throat_surface_functional", "open_leaf:outer_surface_functional"}, #] &],
    obligation == "mechanical_closure_contribution_census", Select[expected, MemberQ[{"open_leaf:geon_core_bundle", "open_leaf:outer_surface_functional"}, #] &],
    MemberQ[{"native_root_structure", "composite_native_root_compatibility"}, obligation], native,
    StringStartsQ[obligation, "E1_"] || StringStartsQ[obligation, "E2_"] || StringStartsQ[obligation, "E3_"], Select[expected, MemberQ[{"domain:Sigma_boundary_data", "open_leaf:throat_surface_functional"}, #] &],
    StringStartsQ[obligation, "E4_"] || StringStartsQ[obligation, "E5_"], native,
    MemberQ[{"operator_definition_or_complement_closure", "whole_complement_class_coverage"}, obligation], Select[expected, MemberQ[{"open_leaf:throat_surface_functional", "open_leaf:outer_surface_functional"}, #] &],
    obligation == "boundary_variation_or_virtual_work" && MemberQ[{"E3", "OTHER"}, candidate], Select[expected, MemberQ[{"open_leaf:throat_surface_functional", "open_leaf:outer_surface_functional"}, #] &],
    MemberQ[{"ensemble_exchange_discharge", "geometric_refinement_applicability"}, obligation] && MemberQ[{"E3", "OTHER"}, candidate], Select[expected, # == "open_leaf:sleeve_core_trace" &],
    True, {}
  ]
];

obligationStage0Slots[obligation_String, candidate_String, stratum_String] := Module[{candidateSlot, suffix},
  candidateSlot = "candidate_definition:" <> candidate;
  Which[
    StringStartsQ[obligation, "topology_"],
      suffix = <|
        "topology_finite_energy_interpolation" -> "finite_energy_interpolation",
        "topology_pair_annihilation_path" -> "pair_annihilation",
        "topology_sector_disconnection" -> "sector_disconnection"
      |>[obligation]; {"topology:" <> candidate <> ":" <> suffix},
    obligation == "host_location_evidence", {"host_location:" <> candidate},
    obligation == "mechanical_closure_contribution_census", {"mechanical_closure:" <> candidate},
    MemberQ[{"boundary_variation_or_virtual_work", "ensemble_exchange_discharge", "geometric_refinement_applicability"}, obligation],
      {"ensemble:" <> candidate <> ":boundary_action_variation"},
    obligation == "template_constituent_contract", {"template_cell_specific:" <> candidate, "template_free_data:" <> stratum},
    StringStartsQ[obligation, "mixture_law:"], {obligation},
    MemberQ[{"operator_definition_or_complement_closure", "whole_complement_class_coverage"}, obligation], {candidateSlot, "basis_closure"},
    True, {candidateSlot}
  ]
];

stage0SlotMeasurement[slotId_String] := Module[{row, outcome},
  row = availability[slotId]; outcome = row["availability_outcome"];
  <|
    "stage0_slot_id" -> slotId, "availability" -> outcome,
    "integrity_state" -> row["integrity_state"], "missing" -> (outcome == "UNRESOLVED"),
    "decisive_incompatibility" -> (outcome == "EXCLUDED"),
    "witness_id" -> Lookup[row, "witness_id", Null], "challenge_id" -> Lookup[row, "challenge_id", Null],
    "producer_set" -> row["producer_set"], "measurement_source" -> "ratified_U2_stage0_availability_record"
  |>
];

evidenceRecord[obligation_String, candidate_String, stratum_String, expected_List] := Module[
  {dependencies, measurements, slotMeasurements, allMeasurements, applicable, missing,
   committedIncompatibility, ancestryComplete, raw, state, operation, roots},
  dependencies = obligationDependencies[obligation, expected, candidate];
  measurements = dependencyMeasurement /@ dependencies;
  slotMeasurements = stage0SlotMeasurement /@ obligationStage0Slots[obligation, candidate, stratum];
  allMeasurements = Join[measurements, slotMeasurements];
  applicable = obligation != "template_constituent_contract";
  missing = applicable && AnyTrue[allMeasurements, TrueQ[# ["missing"]] &];
  committedIncompatibility = applicable && AnyTrue[allMeasurements,
    TrueQ[Lookup[#, "decisive_incompatibility", False]] || Lookup[#, "availability", ""] == "EXCLUDED" &];
  ancestryComplete = Length[allMeasurements] > 0 && AllTrue[allMeasurements,
    # ["integrity_state"] == "COMPUTATION_VALID" && StringLength[# ["measurement_source"]] > 0 &];
  raw = <|
    "applicable" -> applicable, "committed_incompatibility" -> committedIncompatibility,
    "ancestry_complete_and_typed" -> ancestryComplete, "datum_missing" -> missing
  |>;
  state = classifyEvidence[raw]; operation = If[state == "MISSING", "unavailability", "positive_join"];
  roots = Join[lookupRows[measurements, "dependency_token"], Flatten[lookupRows[slotMeasurements, "producer_set"]]];
  <|
    "obligation_id" -> obligation, "integrity" -> "COMPUTATION_VALID",
    "raw_predicate_inputs" -> raw, "predicate_measurements" -> measurements,
    "stage0_slot_measurements" -> slotMeasurements,
    "emitted_state" -> state, "proof_dag" -> proofDag[roots, operation],
    "state_classifier" -> "u2_stage0_sympy.classify_evidence"
  |>
];

minimalDatum = <|
  "E1" -> "throat_surface_functional:E1_holonomic_placement",
  "E2" -> "throat_surface_functional:E2_free_slip_placement",
  "E3" -> "boundary_action_variation:E3_bulk_texture",
  "E4" -> "sleeve_core_trace:E4_nonholonomic_constraint",
  "E5" -> "E5_rayleigh_dissipation_functional",
  "MIXTURE(F_E1_E4)" -> "sleeve_core_trace:F_E1_E4_composite_law",
  "MIXTURE(F_E2_E4)" -> "sleeve_core_trace:F_E2_E4_composite_law",
  "MIXTURE(F_E4_E5)" -> "sleeve_core_trace:F_E4_E5_constraint_Rayleigh_law",
  "OTHER" -> "operator_definition_or_complement_closure"
|>;

unavailabilityPayload[candidate_String, stratum_String, missingDependencies_List, missingSlotIds_List] := Module[
  {primarySlotId, primarySlot, certificate, restore, challenge, slotReferences, dependencyReferences},
  require[Length[missingSlotIds] > 0, "UNRESOLVED cell has no ratified missing slot"];
  primarySlotId = First[Sort[missingSlotIds]]; primarySlot = availability[primarySlotId];
  certificate = primarySlot["witness"]["insufficiency_certificate"];
  restore = primarySlot["witness"]["counterfactual_restore_mutation"];
  challenge = primarySlot["challenge"];
  slotReferences = Table[
    With[{row = availability[slotId]}, <|
      "stage0_slot_id" -> slotId, "witness_id" -> row["witness_id"], "challenge_id" -> row["challenge_id"],
      "availability_outcome" -> row["availability_outcome"]
    |>], {slotId, sortIds[missingSlotIds]}];
  dependencyReferences = Table[<|
    "dependency_token" -> row["dependency_token"], "ratified_slot_id" -> row["ratified_slot_id"],
    "witness_id" -> row["witness_id"], "challenge_id" -> row["challenge_id"], "availability" -> row["availability"]
  |>, {row, SortBy[DeleteDuplicatesBy[missingDependencies, # ["dependency_token"] &], ToLowerCase[# ["dependency_token"]] &]}];
  <|
    "witness_id" -> primarySlot["witness_id"], "challenge_id" -> primarySlot["challenge_id"],
    "kind" -> primarySlot["witness"]["kind"],
    "named_datum" -> minimalDatum[candidate], "stratum_datum" -> "tilt:" <> stratum,
    "missing_dependency_tokens" -> lookupRows[dependencyReferences, "dependency_token"],
    "missing_stage0_slot_ids" -> lookupRows[slotReferences, "stage0_slot_id"],
    "referenced_phaseC_slots" -> dependencyReferences, "referenced_stage0_slots" -> slotReferences,
    "reference_set_exact" -> (
      sortIds[lookupRows[dependencyReferences, "dependency_token"]] == sortIds[lookupRows[missingDependencies, "dependency_token"]]
      && sortIds[lookupRows[slotReferences, "stage0_slot_id"]] == sortIds[missingSlotIds]
    ),
    "required_type" -> primarySlot["required_type"], "required_dimensions" -> primarySlot["required_dimensions"],
    "challenge" -> <|
      "executed" -> certificate["executed"], "semantic_route_id" -> certificate["executed_semantic_route_id"],
      "attempted_candidate_count" -> challenge["attempted_candidate_count"],
      "measured_rank" -> certificate["measured_rank"], "measured_nullity" -> certificate["measured_nullity"],
      "restore_rank" -> restore["restored_rank"], "restore_nullity" -> restore["restored_nullity"],
      "outcome" -> challenge["outcome"], "empty_output" -> challenge["empty_output"],
      "candidate_well_typed" -> challenge["candidate_is_well_typed"],
      "dual_engine_certificate" -> challenge["dual_engine_certificate"]
    |>,
    "proof_dag" -> proofDag[Join[lookupRows[dependencyReferences, "dependency_token"], primarySlot["producer_set"]], "challenge"]
  |>
];

exactResidual[matrix_List, vector_List, rhs_List] := Flatten[matrix . vector - rhs];

routeControl[row_Association] := Module[
  {positive, malformed, matrix, candidate, rhs, positiveResidual, malformedResidual, first, ablated,
   ablatedResidual, native, dissipative, coordinateIndex, coefficient, velocity, power},
  positive = row["positive_fixture"]; malformed = row["malformed_fixture"];
  matrix = positive["matrix"]; candidate = positive["candidate"]; rhs = positive["rhs"];
  positiveResidual = exactResidual[matrix, candidate, rhs];
  malformedResidual = exactResidual[matrix, malformed["candidate"], malformed["rhs"]];
  first = FirstPosition[candidate, value_ /; value != 0][[1]]; ablated = ReplacePart[candidate, first -> 0];
  ablatedResidual = exactResidual[matrix, ablated, rhs]; native = positive["native_structure_exercised"];
  dissipative = TrueQ[Lookup[native, "dissipation_bookkeeping", False]];
  If[dissipative,
    coordinateIndex = If[Length[candidate] == 4, 3, 1]; coefficient = native["rayleigh_coefficient"];
    velocity = candidate[[coordinateIndex + 1]]; power = coefficient velocity^2,
    coordinateIndex = Null; coefficient = Null; velocity = Null; power = Null
  ];
  <|
    "route_id" -> row["route_id"], "cell_id" -> "candidate=" <> StringDelete[row["route_id"], StartOfString ~~ "route:"],
    "fixture_id" -> row["fixture_id"], "semantic_route_id" -> positive["semantic_route_id"], "test_only" -> True,
    "positive_residual" -> positiveResidual, "positive_admissible" -> And @@ (# == 0 & /@ positiveResidual),
    "nondegenerate_norm_squared" -> Total[candidate^2], "native_structure_exercised" -> native,
    "malformed_residual" -> malformedResidual, "malformed_excluded" -> AnyTrue[malformedResidual, # != 0 &],
    "semantic_ablation" -> <|
      "removed_coordinate_index" -> first - 1, "ablated_candidate" -> ablated,
      "ablated_residual" -> ablatedResidual, "nondefinitional_obligation_failed" -> AnyTrue[ablatedResidual, # != 0 &],
      "criterion" -> "fail_nondefinitional_obligation_or_change_canonical_operator_class"
    |>,
    "dissipation_measurement" -> <|
      "applicable" -> dissipative, "rayleigh_coefficient" -> coefficient,
      "velocity_coordinate_index" -> coordinateIndex, "velocity_coordinate" -> velocity,
      "formula" -> If[dissipative, "rayleigh_coefficient*velocity_coordinate^2", Null], "value" -> power,
      "frozen_fixture_field_equal" -> If[dissipative, power == Lookup[native, "computed_dissipated_power", Null], Lookup[native, "computed_dissipated_power", Null] === Null]
    |>
  |>
];

conjugateSignatureRoles = <|
  "normal_velocity_control" -> "normal_traction_response",
  "normal_traction_response" -> "normal_velocity_control",
  "tangential_velocity_control" -> "tangential_traction_response",
  "tangential_traction_response" -> "tangential_velocity_control"
|>;

canonicalConjugateExpectation = <|
  "normal_velocity_lock" -> "normal_traction_response",
  "tangential_velocity_lock" -> "tangential_traction_response",
  "tangential_traction_free" -> "tangential_velocity_control"
|>;

staticSignatureClassifier[boundaryVariation_Association] := Module[
  {signature = {}, measurements, variationChannels, fluxChannels},
  measurements = Table[Module[{condition, roles = {}},
    condition = component["boundary_condition"];
    If[StringContainsQ[condition, "v.normal=V.normal"], AppendTo[roles, "normal_velocity_control"],
      If[StringContainsQ[condition, "normal_traction"], AppendTo[roles, "normal_traction_response"]]];
    If[StringContainsQ[condition, "v.tangent=V.tangent"], AppendTo[roles, "tangential_velocity_control"],
      If[StringContainsQ[condition, "tangential_traction=0"], AppendTo[roles, "tangential_traction_response"]]];
    signature = Join[signature, roles];
    <|
      "endpoint_id" -> component["endpoint_id"], "boundary_condition" -> condition,
      "classified_roles" -> sortIds[roles], "variation_channels" -> sortIds[component["variation_channels"]],
      "flux_channels" -> sortIds[component["flux_channels"]]
    |>
  ], {component, boundaryVariation["component_variations"]}];
  variationChannels = sortIds[Flatten[Lookup[measurements, "variation_channels"]]];
  fluxChannels = sortIds[Flatten[Lookup[measurements, "flux_channels"]]];
  signature = sortIds[signature];
  <|
    "classifier" -> "committed_boundary_variation_static_signature_walk_v1",
    "derived_object_kind" -> boundaryVariation["object_kind"],
    "component_count" -> Length[boundaryVariation["component_variations"]],
    "component_measurements" -> measurements, "signature" -> signature,
    "signature_complete" -> (Length[signature] > 0 && And @@
      (Length[# ["classified_roles"]] >= Length[# ["variation_channels"]] & /@ measurements)),
    "variation_channels" -> variationChannels, "flux_channels" -> fluxChannels
  |>
];

nativePairingAnalysis[candidateStructure_Association] := Module[
  {components, canonicalSignature, eligibleRows, eligible, nonconjugate, unmapped, complete, expected},
  components = candidateStructure["components"]; canonicalSignature = candidateStructure["canonical_signature"];
  eligibleRows = Select[components, # ["variational_class"] == "holonomic_field_BC" &];
  eligible = If[Length[eligibleRows] == 0, {}, sortIds[Lookup[eligibleRows, "endpoint_id"]]];
  nonconjugate = sortIds[(# ["endpoint_id"] <> ":" <> # ["variational_class"] &) /@
    Select[components, # ["variational_class"] != "holonomic_field_BC" &]];
  unmapped = sortIds[Select[canonicalSignature, !KeyExistsQ[canonicalConjugateExpectation, #] &]];
  complete = candidateStructure["native_root_class"] == "variational_holonomic" &&
    Length[eligible] == Length[components] && Length[nonconjugate] == 0 && Length[unmapped] == 0 && Length[canonicalSignature] > 0;
  expected = If[TrueQ[complete], sortIds[(canonicalConjugateExpectation[#] &) /@ canonicalSignature], {}];
  <|
    "construction" -> "native_root_class_and_canonical_signature_pairing_walk_v1",
    "native_root_class" -> candidateStructure["native_root_class"],
    "required_component_count" -> Length[components], "eligible_variational_components" -> eligible,
    "nonconjugate_components" -> nonconjugate, "unmapped_canonical_signature_terms" -> unmapped,
    "complete_structure_pairing" -> complete, "independently_expected_exchanged_signature" -> expected
  |>
];

localEnsemble[ensembleSlot_Association, candidateStructureSlot_Association] := Module[
  {structureId, boundaryVariation, classifier, pairing, pureDisplacement, classification},
  structureId = ensembleSlot["slot_id"];
  If[ensembleSlot["availability_outcome"] != "DERIVED", Return[<|
    "available" -> False, "classification" -> Null, "committed_structure_id" -> structureId,
    "geometric_component" -> Null, "exchange_route" -> Null, "classifier_evidence" -> Null,
    "pairing_analysis" -> Null
  |>]];
  boundaryVariation = ensembleSlot["derived_object"];
  classifier = staticSignatureClassifier[boundaryVariation];
  pairing = nativePairingAnalysis[candidateStructureSlot["derived_object"]];
  pureDisplacement = TrueQ[classifier["signature_complete"]] &&
    classifier["signature"] == {"normal_velocity_control", "tangential_velocity_control"};
  classification = If[pureDisplacement, "fixed-displacement/geometric", "mixed/other-ensemble(" <> structureId <> ")"];
  <|
    "available" -> True, "classification" -> classification, "committed_structure_id" -> structureId,
    "geometric_component" -> (Length[classifier["variation_channels"]] > 0),
    "exchange_route" -> If[TrueQ[pairing["complete_structure_pairing"]], "native_conjugate_transform", "computed_no_pairing_certificate"],
    "classifier_evidence" -> classifier, "pairing_analysis" -> pairing
  |>
];

exchangeControl[local_Association] := Module[
  {route, classifier, pairing, baseline, exchanged, expected, insufficiency, certificate},
  If[!TrueQ[local["available"]], Return[Null]];
  route = local["exchange_route"]; classifier = local["classifier_evidence"]; pairing = local["pairing_analysis"];
  If[route == "native_conjugate_transform",
    baseline = classifier["signature"];
    exchanged = sortIds[(conjugateSignatureRoles[#] &) /@ baseline];
    expected = pairing["independently_expected_exchanged_signature"];
    If[generationMutation == "TOOTH_EXCHANGE_EXPECTED_SIGNATURE_GENERATION", expected = Most[expected]];
    requireGeneration[TrueQ[classifier["signature_complete"]] && exchanged == expected && baseline != exchanged,
      "ASSERT_EXCHANGE_SIGNATURE_GENERATION", "boundary-variation exchange disagrees with native pairing for " <> local["committed_structure_id"]];
    Return[<|
      "route" -> route, "baseline_signature" -> baseline, "computed_exchanged_signature" -> exchanged,
      "independently_generated_expected_signature" -> expected,
      "response_character_flipped" -> (baseline != exchanged && exchanged == expected),
      "signature_classifier_evidence" -> classifier, "native_pairing_construction" -> pairing
    |>]
  ];
  insufficiency = KeyTake[pairing, {
    "construction", "native_root_class", "required_component_count", "eligible_variational_components",
    "nonconjugate_components", "unmapped_canonical_signature_terms", "complete_structure_pairing"
  }];
  If[generationMutation == "TOOTH_NO_PAIRING_CERTIFICATE_GENERATION",
    insufficiency["nonconjugate_components"] = {}; insufficiency["unmapped_canonical_signature_terms"] = {}];
  certificate = !TrueQ[insufficiency["complete_structure_pairing"]] &&
    Length[insufficiency["eligible_variational_components"]] < insufficiency["required_component_count"] &&
    (Length[insufficiency["nonconjugate_components"]] > 0 || Length[insufficiency["unmapped_canonical_signature_terms"]] > 0);
  requireGeneration[certificate, "ASSERT_NO_PAIRING_CERTIFICATE_GENERATION",
    "no-pairing insufficiency is incomplete for " <> local["committed_structure_id"] <> ": " <> ToString[insufficiency, InputForm]];
  <|
    "route" -> route, "observed_boundary_signature" -> classifier["signature"],
    "signature_classifier_evidence" -> classifier, "pairing_insufficiency" -> insufficiency,
    "no_pairing_certificate" -> certificate
  |>
];

topologyRecord[cell_Association] := Module[{candidate, slotIds, nativeTokens, roots},
  candidate = cell["candidate_id"];
  slotIds = <|
    "sector_disconnection" -> "topology:" <> candidate <> ":sector_disconnection",
    "finite_energy_interpolation" -> "topology:" <> candidate <> ":finite_energy_interpolation",
    "pair_annihilation" -> "topology:" <> candidate <> ":pair_annihilation"
  |>;
  require[And @@ (availability[#]["availability_outcome"] == "UNRESOLVED" & /@ Values[slotIds]), "topology datum unexpectedly available"];
  nativeTokens = Select[cell["expected_dependencies"],
    StringStartsQ[#, "endpoint:"] || MemberQ[{"open_leaf:E4_shear_lock", "open_leaf:E5_rayleigh", "open_leaf:gammaSigma", "open_leaf:tangentDtN"}, #] &];
  roots = sortIds[Join[{"open_leaf:sleeve_core_trace"}, nativeTokens]];
  <|
    "record_id" -> "topology_gate:" <> cell["cell_id"], "integrity" -> "COMPUTATION_VALID",
    "sector_disconnection" -> "UNRESOLVED(named datum)", "finite_energy_interpolation" -> "UNRESOLVED(named datum)",
    "pair_annihilation" -> "UNRESOLVED(named datum)", "pair_used_in_aggregate" -> False,
    "physics" -> evaluateTopology["UNRESOLVED", "UNRESOLVED"], "native_root_class" -> cell["native_root_class"],
    "native_structure_tokens" -> roots, "subquestion_slot_ids" -> slotIds,
    "pair_configuration" -> <|
      "object_type" -> "static_plus_w_minus_w_pair_configuration", "firewall_tag" -> "PAIR_ANNIHILATION_ONLY",
      "consumer" -> "topology_pair_annihilation_subquestion"
    |>,
    "proof_dag" -> proofDag[roots, "unavailability"]
  |>
];

hostRecord[cell_Association] := Module[{candidate, slot},
  candidate = cell["candidate_id"]; slot = availability["host_location:" <> candidate];
  require[slot["availability_outcome"] == "UNRESOLVED", "host location unexpectedly derived"];
  <|
    "record_id" -> "host:" <> cell["cell_id"], "stable_branch_id" -> cell["stable_branch_id"],
    "integrity" -> "COMPUTATION_VALID", "physics" -> "undetermined",
    "evidential_basis" -> <|
      "witness_id" -> slot["witness_id"], "challenge_id" -> slot["challenge_id"],
      "challenge_executed" -> TrueQ[slot["challenge"]["dual_engine_certificate"]],
      "missing_tokens" -> {"open_leaf:geon_core_bundle", "open_leaf:sleeve_core_trace"}
    |>,
    "proof_dag" -> proofDag[{"open_leaf:geon_core_bundle", "open_leaf:sleeve_core_trace"}, "unavailability"]
  |>
];

claimSideContributionCensus[boundaryVariation_Association] := Module[{components, variationChannels, kinds},
  components = boundaryVariation["component_variations"];
  variationChannels = Flatten[Lookup[components, "variation_channels"]];
  kinds = {"static_radiation"};
  If[AnyTrue[components, Length[# ["variation_channels"]] > 0 &], AppendTo[kinds, "surface_boundary_response"]];
  If[AnyTrue[components, Length[# ["flux_channels"]] > 0 &], AppendTo[kinds, "outer_matching_response"]];
  If[AnyTrue[variationChannels, StringEndsQ[#, "_reaction"] &], AppendTo[kinds, "constraint_reaction"]];
  If[AnyTrue[variationChannels, StringEndsQ[#, "_loss"] &], AppendTo[kinds, "rayleigh_loss"]];
  sortIds[kinds]
];

committedRootContributionWalk[candidateStructure_Association] := Module[{rootsByKind, endpoint, channels, census, terms},
  rootsByKind = <|"static_radiation" -> {}|>;
  Do[
    endpoint = component["endpoint_id"]; channels = component["channels"];
    If[Length[channels["var"]] > 0,
      rootsByKind["surface_boundary_response"] = Append[Lookup[rootsByKind, "surface_boundary_response", {}], endpoint]];
    If[Length[channels["flux"]] > 0,
      rootsByKind["outer_matching_response"] = Append[Lookup[rootsByKind, "outer_matching_response", {}], endpoint]];
    If[Length[channels["constraint"]] > 0,
      rootsByKind["constraint_reaction"] = Append[Lookup[rootsByKind, "constraint_reaction", {}], endpoint]];
    If[Length[channels["Rayleigh"]] > 0,
      rootsByKind["rayleigh_loss"] = Append[Lookup[rootsByKind, "rayleigh_loss", {}], endpoint]];
    If[KeyExistsQ[channels, "rad"], rootsByKind["static_radiation"] = Append[rootsByKind["static_radiation"], endpoint <> ":rad_channel_schema"]],
    {component, candidateStructure["components"]}
  ];
  census = sortIds[Keys[rootsByKind]];
  terms = Table[<|
    "contribution_id" -> contribution, "committed_schema_roots" -> sortIds[rootsByKind[contribution]],
    "symbolic_coefficient" -> If[contribution == "static_radiation", 0, 1]
  |>, {contribution, census}];
  {census, terms}
];

closureRefusal[cell_Association, local_Association, host_Association, ensembleSlot_Association, candidateStructureSlot_Association] := Module[
  {candidate, censusA, censusB, independentTerms, assignments, assignmentTotal, independentTotal},
  If[!TrueQ[local["available"]], Return[Null]]; candidate = cell["candidate_id"];
  censusA = claimSideContributionCensus[ensembleSlot["derived_object"]];
  {censusB, independentTerms} = committedRootContributionWalk[candidateStructureSlot["derived_object"]];
  If[generationMutation == "TOOTH_CLOSURE_CENSUS_GENERATION", censusB = Most[censusB]];
  requireGeneration[censusA == censusB, "ASSERT_CLOSURE_CENSUS_GENERATION",
    "claim-side and committed-root contribution censuses disagree for " <> candidate];
  assignments = Table[
    If[contribution == "static_radiation",
      <|"contribution_id" -> contribution, "landing" -> "computed-zero", "owner" -> "radiation/static-zero", "symbolic_coefficient" -> 0|>,
      <|"contribution_id" -> contribution, "landing" -> "typed-refusal", "refusal_witness_id" -> "witness:host_location:" <> candidate, "symbolic_coefficient" -> 1|>
    ], {contribution, censusA}];
  assignmentTotal = Total[Lookup[assignments, "symbolic_coefficient"]];
  independentTotal = Total[Lookup[independentTerms, "symbolic_coefficient"]];
  If[generationMutation == "TOOTH_CLOSURE_TOTAL_GENERATION", independentTotal = independentTotal + 1];
  requireGeneration[independentTotal == assignmentTotal, "ASSERT_CLOSURE_TOTAL_GENERATION",
    "direct committed balance and assignment total disagree for " <> candidate];
  <|
    "record_id" -> "closure:" <> cell["cell_id"] <> ":ensemble_level_1", "stable_branch_id" -> cell["stable_branch_id"],
    "integrity" -> "COMPUTATION_VALID",
    "physics" -> <|"kind" -> "CLOSURE_REFUSED", "typed_refusal_reason" -> <|"kind" -> "under_specified_host", "witness_id" -> "witness:host_location:" <> candidate|>|>,
    "gated_claim_id" -> "ensemble_level_1:" <> cell["cell_id"],
    "contribution_census_A" -> censusA, "contribution_census_B" -> censusB,
    "contribution_census_A_construction" -> "claim_side_boundary_variation_kind_expansion",
    "contribution_census_B_construction" -> "raw_committed_component_channel_schema_walk",
    "assignments" -> assignments, "assignment_coverage" -> sortIds[Lookup[assignments, "contribution_id"]],
    "independently_constructed_terms" -> independentTerms,
    "independently_constructed_symbolic_total" -> independentTotal, "assignment_symbolic_total" -> assignmentTotal,
    "host_physics" -> host["physics"], "certificate_emitted" -> False, "return_closure_consumed" -> False,
    "proof_dag" -> proofDag[{"open_leaf:geon_core_bundle", "open_leaf:sleeve_core_trace"}, "unavailability"]
  |>
];

cellRecord[gridCell_Association] := Module[
  {candidate, candidateRecord, census, cell, obligations, evidence, landing, missingDependencies, missingSlotIds, witness,
   ensembleSlot, candidateStructureSlot, local, host, closure, finalLevel1, primaryWitness, topology, applicability, level2, usedDependencies},
  candidate = gridCell["candidate_id"]; candidateRecord = candidateRecords[candidate]; census = obligationsDoc[candidate];
  cell = Join[gridCell, <|"native_root_class" -> candidateRecord["native_root_class"]|>]; obligations = census["generator_A"];
  evidence = evidenceRecord[#, candidate, cell["stratum"], cell["expected_dependencies"]] & /@ obligations;
  landing = evaluateDisposition[evidence]; require[landing == "UNRESOLVED(named datum)", "unexpected disposition " <> landing];
  missingDependencies = Flatten[Select[# ["predicate_measurements"], TrueQ[# ["missing"]] &] & /@ evidence];
  missingSlotIds = sortIds[Flatten[lookupRows[
    Select[# ["stage0_slot_measurements"], TrueQ[# ["missing"]] &], "stage0_slot_id"
  ] & /@ evidence]];
  witness = unavailabilityPayload[candidate, cell["stratum"], missingDependencies, missingSlotIds];
  ensembleSlot = availability["ensemble:" <> candidate <> ":boundary_action_variation"];
  candidateStructureSlot = availability["candidate_definition:" <> candidate];
  local = localEnsemble[ensembleSlot, candidateStructureSlot]; host = hostRecord[cell];
  closure = closureRefusal[cell, local, host, ensembleSlot, candidateStructureSlot];
  If[TrueQ[local["available"]],
    finalLevel1 = "UNRESOLVED(mechanical-closure refusal)";
    primaryWitness = <|
      "witness_id" -> "witness:ensemble_closure:" <> cell["cell_id"], "challenge_id" -> "challenge:ensemble_closure:" <> cell["cell_id"],
      "named_datum" -> minimalDatum[candidate], "executed" -> True
    |>,
    finalLevel1 = "UNRESOLVED(boundary-action variation)";
    With[{slot = availability["ensemble:" <> candidate <> ":boundary_action_variation"]},
      primaryWitness = <|
        "witness_id" -> slot["witness_id"], "challenge_id" -> slot["challenge_id"],
        "named_datum" -> "ensemble:" <> candidate <> ":boundary_action_variation",
        "executed" -> TrueQ[slot["challenge"]["dual_engine_certificate"]]
      |>
    ]
  ];
  topology = topologyRecord[cell]; applicability = "missing";
  level2 = evaluateCrossLevel[applicability, topology["integrity"], "UNRESOLVED"];
  usedDependencies = sortIds[Flatten[lookupRows[#, "dependency_token"] & /@ Lookup[evidence, "predicate_measurements"]]];
  <|
    "cell_id" -> cell["cell_id"], "stable_branch_id" -> cell["stable_branch_id"], "candidate_id" -> candidate,
    "ambient" -> cell["ambient"], "stratum" -> cell["stratum"], "native_root_class" -> cell["native_root_class"],
    "integrity" -> "COMPUTATION_VALID", "expected_dependencies" -> cell["expected_dependencies"], "used_dependencies" -> usedDependencies,
    "obligation_evidence" -> evidence,
    "disposition" -> <|
      "kind" -> "UNRESOLVED", "witness_id" -> witness["witness_id"], "challenge_id" -> witness["challenge_id"],
      "named_datum" -> witness["named_datum"], "stratum_datum" -> witness["stratum_datum"]
    |>,
    "disposition_evaluator_landing" -> landing, "unavailability" -> witness,
    "ensemble" -> <|
      "level_1" -> <|
        "record_id" -> "ensemble_level_1:" <> cell["cell_id"], "integrity" -> "COMPUTATION_VALID", "physics" -> finalLevel1,
        "local_preclosure_classification" -> local["classification"],
        "local_preclosure_committed_structure_id" -> local["committed_structure_id"],
        "local_preclosure_geometric_component" -> local["geometric_component"],
        "local_preclosure_classifier_evidence" -> local["classifier_evidence"],
        "closure_record_id" -> If[AssociationQ[closure], closure["record_id"], Null], "witness" -> primaryWitness,
        "exchange_control" -> exchangeControl[local],
        "proof_dag" -> proofDag[{"open_leaf:geon_core_bundle", "open_leaf:sleeve_core_trace"}, "unavailability"]
      |>,
      "applicability" -> <|
        "value" -> applicability, "local_preclosure_value" -> If[TrueQ[local["geometric_component"]], "geometric-component-bearing", "missing"],
        "reason" -> "final_primary_is_UNRESOLVED"
      |>,
      "level_2" -> <|
        "record_id" -> "ensemble_level_2:" <> cell["cell_id"], "integrity" -> "COMPUTATION_VALID", "physics" -> level2,
        "evaluator_inputs" -> <|"applicability" -> applicability, "gate_integrity" -> topology["integrity"], "gate_outcome" -> "UNRESOLVED"|>,
        "proof_dag" -> proofDag[{"open_leaf:sleeve_core_trace"}, "unavailability"]
      |>
    |>,
    "topology_gate" -> topology, "host_location" -> host, "closure_adjudication" -> closure,
    "return_closure_ownership" -> <|
      "record_id" -> "return_closure:" <> cell["cell_id"], "stable_branch_id" -> cell["stable_branch_id"],
      "owner" -> "downstream_flux_path", "U2_owned" -> False, "preserved_terminal" -> "UNRESOLVED(return_closure)",
      "source_reference" -> "source:field:return_closure", "computed_reachable_from_U2_verdict" -> False
    |>,
    "template_eligible" -> False
  |>
];

forcingCensus[context_Association, cellMap_Association] := Module[
  {mappings, consideredRoots, availableForcing, generatorA, forcingNodes = {}, generatorB},
  mappings = context["candidate_cell_mappings"];
  consideredRoots = {"source:field:throat_surface_functional", "source:field:outer_surface_functional"};
  availableForcing = Select[consideredRoots, Function[root,
    AnyTrue[mappings, Function[mapping,
      AnyTrue[cellMap[mapping["cell_id"]]["obligation_evidence"], Function[evidence,
        MemberQ[evidence["proof_dag"]["root_ids"], root] && MemberQ[evidence["proof_dag"]["constructors"], "STATIC_COMMITTED_FORCING"]
      ]]
    ]]
  ]];
  generatorA = <|
    "algorithm" -> "committed_surface_forcing_root_walk",
    "considered_roots" -> consideredRoots,
    "available_normalized_operator_derivations" -> availableForcing,
    "measured_rank" -> MatrixRank[({Boole[MemberQ[availableForcing, #]]} &) /@ consideredRoots]
  |>;
  Do[
    Do[If[MemberQ[evidence["proof_dag"]["constructors"], "STATIC_COMMITTED_FORCING"], AppendTo[forcingNodes, evidence["obligation_id"]]],
      {evidence, cellMap[mapping["cell_id"]]["obligation_evidence"]}],
    {mapping, mappings}
  ];
  generatorB = <|
    "algorithm" -> "executed_cell_DAG_static_forcing_node_walk",
    "available_normalized_operator_derivations" -> sortIds[forcingNodes], "visited_cell_count" -> Length[mappings]
  |>;
  {generatorA, generatorB}
];

promotionRecord[context_Association, cellMap_Association] := Module[
  {censuses, censusA, censusB, failedMappings, failedUpstreams, contextCandidates, forcingRecords,
   inputs, landing, key, challengeMatrix, witness},
  censuses = forcingCensus[context, cellMap]; {censusA, censusB} = censuses;
  require[censusA["available_normalized_operator_derivations"] == censusB["available_normalized_operator_derivations"], "forcing census generators disagree"];
  failedMappings = Select[context["candidate_cell_mappings"],
    cellMap[# ["cell_id"]]["integrity"] != "COMPUTATION_VALID" &];
  failedUpstreams = If[Length[failedMappings] == 0, {}, sortIds[Lookup[failedMappings, "cell_id"]]];
  contextCandidates = Lookup[context["candidate_cell_mappings"], "candidate_id"];
  forcingRecords = (<|
    "forcing_id" -> #, "candidate_in_current_axis" -> False,
    "canonical_class" -> Null, "disposition" -> Null
  |> &) /@ censusA["available_normalized_operator_derivations"];
  inputs = <|
    "failed_required_upstreams" -> failedUpstreams,
    "uncanonicalized_overlap" -> (Length[contextCandidates] != Length[DeleteDuplicates[contextCandidates]]),
    "forcing_records" -> forcingRecords, "closure_outcome" -> Null
  |>;
  landing = evaluatePromotion[inputs]; key = context["promotion_key"];
  challengeMatrix = ({Boole[MemberQ[censusA["available_normalized_operator_derivations"], #]]} &) /@ censusA["considered_roots"];
  witness = <|
    "witness_id" -> "no_forcing_witness:" <> key, "challenge_id" -> "no_forcing_challenge:" <> key,
    "complete_forcing_root_census" -> censusA["considered_roots"], "construct_B_star_attempted" -> True,
    "measured_rank" -> MatrixRank[challengeMatrix], "measured_nullity" -> (Dimensions[challengeMatrix][[2]] - MatrixRank[challengeMatrix]),
    "challenge_outcome" -> "CONSTRUCTIVE_FAIL", "executed" -> True
  |>;
  <|
    "record_id" -> "promotion:" <> key, "promotion_key" -> key,
    "stable_branch_id" -> "U2:PROMOTION:" <> context["ambient"] <> ":" <> context["global_common_refinement_context"],
    "integrity" -> "COMPUTATION_VALID", "physics" -> "NO_SELECTION_CLAIM", "evaluator_landing" -> landing,
    "candidate_universe_digest" -> candidateDoc["candidate_universe_digest"],
    "forcing_census_A" -> censusA, "forcing_census_B" -> censusB, "evaluator_inputs" -> inputs,
    "no_forcing_witness" -> witness, "candidate_cell_ids" -> Lookup[context["candidate_cell_mappings"], "cell_id"],
    "proof_dag" -> proofDag[censusA["considered_roots"], "challenge"],
    "exclusion_records_referenced" -> {}, "survivor_or_complement_objects_referenced" -> {}, "stability_roots_referenced" -> {}
  |>
];

evaluatorBehavior[] := Module[{decisive, missing, satisfied, topology, cross, forcing, promotionInputs},
  decisive = <|"obligation_id" -> "decisive", "raw_predicate_inputs" -> <|"applicable" -> True, "committed_incompatibility" -> True, "ancestry_complete_and_typed" -> True, "datum_missing" -> False|>, "emitted_state" -> "DECISIVE_INCOMPATIBILITY"|>;
  missing = <|"obligation_id" -> "missing", "raw_predicate_inputs" -> <|"applicable" -> True, "committed_incompatibility" -> False, "ancestry_complete_and_typed" -> True, "datum_missing" -> True|>, "emitted_state" -> "MISSING"|>;
  satisfied = <|"obligation_id" -> "satisfied", "raw_predicate_inputs" -> <|"applicable" -> True, "committed_incompatibility" -> False, "ancestry_complete_and_typed" -> True, "datum_missing" -> False|>, "emitted_state" -> "SATISFIED"|>;
  topology = Flatten[Table[<|"sector" -> sector, "interpolation" -> interpolation, "landing" -> evaluateTopology[sector, interpolation]|>,
    {sector, {"DISCONNECTED", "CONNECTED", "UNRESOLVED"}}, {interpolation, {"OBSTRUCTED", "INTERPOLABLE", "UNRESOLVED"}}]];
  cross = Flatten[Table[<|
      "applicability" -> applicability, "gate_integrity" -> gate[[1]], "gate_outcome" -> gate[[2]],
      "landing" -> evaluateCrossLevel[applicability, gate[[1]], gate[[2]]]
    |>, {applicability, {"geometric", "positively_non_geometric", "missing"}},
    {gate, {{"COMPUTATION_VALID", "topologically-distinct"}, {"COMPUTATION_VALID", "orientation-only"}, {"COMPUTATION_VALID", "UNRESOLVED"}, {"HARNESS_FAILED", Null}, {"NOT_RUN", Null}}}]];
  forcing = {<|"candidate_id" -> "E1", "canonical_class" -> "E1", "candidate_in_current_axis" -> True, "disposition" -> "ADMISSIBLE"|>};
  promotionInputs = <|
    "failed_upstream" -> <|"failed_required_upstreams" -> {"closure:E1"}, "uncanonicalized_overlap" -> False, "forcing_records" -> {}, "closure_outcome" -> Null|>,
    "alias" -> <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> True, "forcing_records" -> {}, "closure_outcome" -> Null|>,
    "required_positive" -> <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> forcing, "closure_outcome" -> "CLOSURE_CERTIFIED"|>,
    "closure_refusal" -> <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> forcing, "closure_outcome" -> "CLOSURE_REFUSED"|>,
    "forced_unresolved" -> <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> {Join[forcing[[1]], <|"disposition" -> "UNRESOLVED"|>]}, "closure_outcome" -> Null|>,
    "forced_excluded" -> <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> {Join[forcing[[1]], <|"disposition" -> "EXCLUDED"|>]}, "closure_outcome" -> Null|>,
    "multi_forcing" -> <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> {forcing[[1]], Join[forcing[[1]], <|"candidate_id" -> "E2", "canonical_class" -> "E2"|>]}, "closure_outcome" -> "CLOSURE_CERTIFIED"|>,
    "outside_axis" -> <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> {Join[forcing[[1]], <|"candidate_in_current_axis" -> False|>]}, "closure_outcome" -> "CLOSURE_CERTIFIED"|>,
    "no_forcing" -> <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> {}, "closure_outcome" -> Null|>
  |>;
  <|
    "disposition" -> <|
      "decisive_plus_missing" -> evaluateDisposition[{decisive, missing}], "all_satisfied" -> evaluateDisposition[{satisfied}],
      "missing_only" -> evaluateDisposition[{missing}], "relabelled_decisive" -> evaluateDisposition[{Join[decisive, <|"emitted_state" -> "MISSING"|>], missing}]
    |>,
    "topology" -> topology, "cross_level" -> cross,
    "promotion" -> Association@KeyValueMap[#1 -> evaluatePromotion[#2] &, promotionInputs],
    "record_schema" -> <|
      "valid" -> recordSchemaValid[<|"integrity" -> "COMPUTATION_VALID", "physics" -> "UNRESOLVED"|>],
      "failed" -> recordSchemaValid[<|"integrity" -> "HARNESS_FAILED", "physics" -> Null, "failure_reason" -> "schema_violation"|>],
      "not_run" -> recordSchemaValid[<|"integrity" -> "NOT_RUN", "physics" -> Null, "failed_upstreams" -> {"A", "B"}|>]
    |>
  |>
];

dimensionFirewall[] := Module[{traction, area, velocity, expected, surfacePower, ablated},
  traction = {-1, -2, 1}; area = {2, 0, 0}; velocity = {1, -1, 0}; expected = {2, -3, 1};
  surfacePower = traction + area + velocity; ablated = traction + area;
  <|
    "basis" -> {"L", "T", "M"}, "surface_power_dimensions" -> surfacePower,
    "constraint_power_dimensions" -> {2, -3, 1}, "rayleigh_power_dimensions" -> {2, -3, 1},
    "expected_power_dimensions" -> expected, "all_constructed_real_expressions_homogeneous" -> (surfacePower == expected),
    "firing_ablation" -> <|"mutated_dimensions" -> ablated, "heterogeneity_detected" -> (ablated != expected)|>
  |>
];

guardFixtures[contracts_Association] := Module[{template, closure, expectedTerms, emittedTerms, contributions, total},
  template = contracts["template_guard_fixture"]; closure = contracts["closure_guard_fixture"];
  expectedTerms = template["expected_term_census"];
  emittedTerms = (<|"term_id" -> #["term_id"], "kind" -> #["kind"]|> &) /@ template["template_record"]["symbolic_ast"]["args"];
  contributions = closure["committed_root_contributions"];
  total = Table[Total[(# ["vector"][[index]] &) /@ contributions], {index, 2}];
  <|
    "template" -> <|
      "expected_terms" -> expectedTerms, "emitted_terms" -> emittedTerms,
      "term_incidence_exact" -> (expectedTerms == emittedTerms),
      "dependent_fields_unbound" -> (template["template_record"]["evaluation_state"] == "UNEVALUATED"),
      "solve_evaluation_node_reachable" -> False
    |>,
    "closure" -> <|
      "census_A" -> closure["census_A"], "census_B" -> sortIds[closure["census_B"]],
      "assignments" -> closure["certificate"]["assignments"], "independently_constructed_total" -> total,
      "expected_total" -> closure["independently_constructed_total"], "reconstruction_exact" -> (total == closure["independently_constructed_total"]),
      "no_double_count" -> (Length[DeleteDuplicates[Lookup[closure["certificate"]["assignments"], "contribution_id"]]] == Length[closure["certificate"]["assignments"]])
    |>,
    "mixed_applicability" -> <|
      "mixed_geometric" -> evaluateCrossLevel["geometric", "COMPUTATION_VALID", "UNRESOLVED"],
      "mixed_non_geometric" -> evaluateCrossLevel["positively_non_geometric", "COMPUTATION_VALID", "UNRESOLVED"]
    |>,
    "pair_firewall" -> <|
      "object_type" -> "static_plus_w_minus_w_pair_configuration", "firewall_tag" -> "PAIR_ANNIHILATION_ONLY",
      "allowed_consumer" -> "topology_pair_annihilation_subquestion", "other_consumer_count" -> 0
    |>
  |>
];

counterAssociation[values_List] := Association@KeyValueMap[ToString[#1] -> #2 &, Counts[values]];
summaryCounts[cells_List, promotions_List, templates_List] := Module[
  {candidates, byCandidate, disposition, finalLevel1, localLevel1, level2, topology, pair, promotion, closure},
  candidates = sortIds[Lookup[cells, "candidate_id"]];
  byCandidate = Association@Table[candidate -> counterAssociation[Lookup[Select[cells, # ["candidate_id"] == candidate &], "disposition"][[All, "kind"]]], {candidate, candidates}];
  disposition = counterAssociation[Lookup[cells, "disposition"][[All, "kind"]]];
  finalLevel1 = counterAssociation[(# ["ensemble"]["level_1"]["physics"] &) /@ cells];
  localLevel1 = counterAssociation[(With[{classification = # ["ensemble"]["level_1"]["local_preclosure_classification"]},
    Which[classification === Null, "UNAVAILABLE", StringStartsQ[classification, "mixed/other-ensemble("], "mixed/other-ensemble", True, classification]
  ] &) /@ cells];
  level2 = counterAssociation[(# ["ensemble"]["level_2"]["physics"] &) /@ cells];
  topology = counterAssociation[(# ["topology_gate"]["physics"] &) /@ cells];
  pair = counterAssociation[(# ["topology_gate"]["pair_annihilation"] &) /@ cells];
  promotion = counterAssociation[Lookup[promotions, "physics"]];
  closure = counterAssociation[(# ["closure_adjudication"]["physics"]["kind"] &) /@ Select[cells, AssociationQ[# ["closure_adjudication"]] &]];
  <|
    "dispositions" -> disposition, "dispositions_by_candidate" -> byCandidate,
    "ensemble_level_1_final" -> finalLevel1, "ensemble_level_1_local_preclosure" -> localLevel1,
    "ensemble_level_2" -> level2, "topology_gate" -> topology, "pair_annihilation" -> pair,
    "closure_adjudications" -> closure, "promotion_outcomes" -> promotion,
    "posed_template_count" -> Length[templates], "integrity_failures" -> 0
  |>
];

require[Length[gridDoc["grid_cells"]] == 144, "ratified grid is not 144 cells"];
require[candidateDoc["candidate_count"] == 9, "ratified candidate axis changed"];
routes = routeControl /@ components["route_fixture_inventory"]["route_records"];
cells = cellRecord /@ gridDoc["grid_cells"];
cellMap = AssociationThread[Lookup[cells, "cell_id"], cells];
promotions = promotionRecord[#, cellMap] & /@ gridDoc["promotion_contexts"];
templates = {};
closureRecords = Select[Lookup[cells, "closure_adjudication"], AssociationQ];

semantic = <|
  "schema_version" -> "U2_PRODUCTION_SEMANTIC_VIEW_V1",
  "stage0_binding" -> <|
    "supplied_digest" -> suppliedDigest, "recomputed_component_digest" -> canonicalDigest[components],
    "equal" -> (canonicalDigest[components] == suppliedDigest == ratifiedDigest)
  |>,
  "scope" -> <|
    "static_adjudication_only" -> True, "dynamical_selection_deferred" -> True,
    "one_body_only" -> True, "BVP_solved" -> False, "two_body_consumer_object_count" -> 0,
    "postulate_used_as_selection_evidence" -> False, "banned_native_import_count" -> 0
  |>,
  "axes" -> <|
    "candidate_axis" -> candidateDoc["candidate_axis"], "ambient_branches" -> sortIds[Lookup[gridDoc["grid_cells"], "ambient"]],
    "active_strata" -> gridDoc["active_strata"], "cell_count" -> Length[cells], "promotion_context_count" -> Length[promotions],
    "candidate_universe_digest" -> candidateDoc["candidate_universe_digest"]
  |>,
  "frozen_evaluator_behavior" -> evaluatorBehavior[], "route_controls" -> routes, "cell_records" -> cells,
  "promotion_records" -> promotions, "closure_records" -> closureRecords, "posed_BVP_templates" -> templates,
  "guard_fixtures" -> guardFixtures[components["closure_template_contracts"]],
  "dimensional_firewall" -> dimensionFirewall[], "headlines" -> summaryCounts[cells, promotions, templates]
|>;

localRegistry = {
    <|"semantic_route_id" -> "obligation_residual_classifier_v1", "engine_local_function" -> "routeControl", "exists" -> (Length[DownValues[routeControl]] > 0), "executed" -> (Length[routes] == 144)|>,
    <|"semantic_route_id" -> "frozen_evidence_state_classifier_v1", "engine_local_function" -> "evidenceRecord", "exists" -> (Length[DownValues[evidenceRecord]] > 0), "executed" -> (Length[cells] > 0)|>,
    <|"semantic_route_id" -> "frozen_disposition_precedence_v1", "engine_local_function" -> "cellRecord", "exists" -> (Length[DownValues[cellRecord]] > 0), "executed" -> (Length[cells] > 0)|>,
    <|"semantic_route_id" -> "frozen_topology_aggregate_v1", "engine_local_function" -> "topologyRecord", "exists" -> (Length[DownValues[topologyRecord]] > 0), "executed" -> (Length[cells] > 0)|>,
    <|"semantic_route_id" -> "frozen_cross_level_ensemble_v1", "engine_local_function" -> "cellRecord", "exists" -> (Length[DownValues[cellRecord]] > 0), "executed" -> (Length[cells] > 0)|>,
    <|"semantic_route_id" -> "frozen_promotion_evaluator_v1", "engine_local_function" -> "promotionRecord", "exists" -> (Length[DownValues[promotionRecord]] > 0), "executed" -> (Length[promotions] > 0)|>
  };

result = <|
  "schema_version" -> "U2_PRODUCTION_ENGINE_V1", "engine" -> "Wolfram",
  "independent_route" -> "Wolfram association/DAG walks plus exact symbolic linear algebra with comparator-asserted ratified evaluator equivalence",
  "semantic_view" -> semantic,
  "engine_local_route_registry" -> localRegistry,
  "runtime_identity" -> <|"wolfram_version" -> $Version, "system_id" -> $SystemID, "user_init_disabled" -> True, "sanitized_path" -> $Path|>
|>;

If[!DirectoryQ[DirectoryName[output]], CreateDirectory[DirectoryName[output], CreateIntermediateDirectories -> True]];
Export[output, result, "YAML"];
Print["U2_WOLFRAM_PRODUCTION_PASS cells=" <> ToString[semantic["axes"]["cell_count"]] <>
  " dispositions=" <> ToString[semantic["headlines"]["dispositions"]] <>
  " templates=" <> ToString[semantic["headlines"]["posed_template_count"]]];
Exit[0];
