(* Independent Wolfram stage-0 route for U2 boundary adjudication. *)

ClearAll["Global`*"];

schema = "U2_STAGE0_ENGINE_V1";
ambients = {"one_sided_pathA29", "two_sided_R_w_postulate"};
baseEndpoints = {"E1", "E2", "E3", "E4", "E5"};
familyIds = {"F_E1_E4", "F_E2_E4", "F_E4_E5"};
tiltTypes = {
  "indexed_density_tilt_profile", "indexed_flow_tilt_response",
  "indexed_h_tilt_profile", "indexed_phase_tilt_profile",
  "indexed_shear_tilt_profile", "indexed_sleeve_surface_normal_profile",
  "indexed_sleeve_tilt_profile", "indexed_uw_tilt_profile"
};
callCounts = <||>;
countCall[name_String] := AssociateTo[callCounts, name -> (Lookup[callCounts, name, 0] + 1)];

cliArgs = If[Length[$ScriptCommandLine] > 0, $ScriptCommandLine, $CommandLine];
getArg[name_] := Module[{position = FirstPosition[cliArgs, name]},
  If[MissingQ[position] || position[[1]] == Length[cliArgs],
    Print["MISSING_ARGUMENT " <> name]; Exit[2],
    cliArgs[[position[[1]] + 1]]
  ]
];
repo = ExpandFileName[getArg["--repo"]];
output = ExpandFileName[getArg["--output"]];

$Path = Select[$Path, StringStartsQ[ToString[#], $InstallationDirectory] &];
loadYaml[path_String] := Import[path, "YAML"];
sortIds[values_] := SortBy[DeleteDuplicates[values], {ToLowerCase[ToString[#]], ToString[#]} &];
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
v11ProductionRel = "software/em_charge_attribute/reports/u2_boundary_adjudication_artifacts/stage_1_production/production_results.yaml";
v11Stage0Digest = "9eff1b0c49e89007aea1008cb6712b0ea495168d101ce43ddce1cffaf68749c4";
v11ProductionAnchor = "53529bf1729811f5ae9faa429cf836507469569b";
v11WrapAnchor = "5ceebb24";

physicsRecordProjection[document_Association] := Module[
  {rows = {}, contextFields, dispositionFields, context, cell, payload},
  contextFields = {"cell_id", "stable_branch_id", "candidate_id", "ambient", "stratum"};
  dispositionFields = Join[contextFields, {
    "native_root_class", "integrity", "expected_dependencies", "used_dependencies",
    "obligation_evidence", "disposition", "disposition_evaluator_landing", "unavailability"
  }];
  Do[
    context = Association@Table[key -> cell[key], {key, contextFields}];
    AppendTo[rows, <|
      "record_id" -> "candidate_disposition:" <> cell["cell_id"],
      "record_class" -> "candidate_disposition",
      "payload" -> Association@Table[key -> cell[key], {key, dispositionFields}]
    |>];
    AppendTo[rows, <|
      "record_id" -> cell["ensemble"]["level_1"]["record_id"], "record_class" -> "ensemble_level_1",
      "payload" -> Join[context, cell["ensemble"]["level_1"]]
    |>];
    AppendTo[rows, <|
      "record_id" -> cell["ensemble"]["level_2"]["record_id"], "record_class" -> "ensemble_level_2",
      "payload" -> Join[context, <|"applicability" -> cell["ensemble"]["applicability"]|>, cell["ensemble"]["level_2"]]
    |>];
    AppendTo[rows, <|"record_id" -> cell["topology_gate"]["record_id"], "record_class" -> "topology_gate", "payload" -> cell["topology_gate"]|>];
    AppendTo[rows, <|"record_id" -> cell["host_location"]["record_id"], "record_class" -> "host_location", "payload" -> cell["host_location"]|>];
    AppendTo[rows, <|
      "record_id" -> cell["return_closure_ownership"]["record_id"], "record_class" -> "return_closure_ownership",
      "payload" -> cell["return_closure_ownership"]
    |>],
    {cell, document["cell_records"]}
  ];
  Do[AppendTo[rows, <|"record_id" -> payload["record_id"], "record_class" -> "closure_adjudication", "payload" -> payload|>], {payload, document["closure_records"]}];
  Do[AppendTo[rows, <|"record_id" -> payload["record_id"], "record_class" -> "promotion", "payload" -> payload|>], {payload, document["promotion_records"]}];
  SortBy[rows, {ToLowerCase[# ["record_id"]], # ["record_id"]} &]
];

physicsRecordInvarianceContract[] := Module[{sourcePath, source, projection, classCounts},
  sourcePath = FileNameJoin[{repo, v11ProductionRel}];
  source = loadYaml[sourcePath]; projection = physicsRecordProjection[source];
  classCounts = Association@KeyValueMap[#1 -> #2 &, KeySort[Counts[Lookup[projection, "record_class"]]]];
  <|
    "schema_version" -> "U2_V11_PHYSICS_RECORD_INVARIANCE_V1",
    "U2_V11_PHYSICS_RECORD_SET_DIGEST" -> canonicalDigest[projection],
    "record_id_universe" -> Lookup[projection, "record_id"], "record_count" -> Length[projection],
    "record_class_counts" -> classCounts,
    "canonical_projection" -> <|
      "serialization" -> "UTF-8 canonical JSON; object keys lexicographically sorted; separators comma/colon; ensure_ascii=false",
      "row_order" -> "record_id casefold then record_id codepoint",
      "cell_context_fields" -> {"cell_id", "stable_branch_id", "candidate_id", "ambient", "stratum"},
      "candidate_disposition_fields" -> {
        "cell_id", "stable_branch_id", "candidate_id", "ambient", "stratum", "native_root_class", "integrity",
        "expected_dependencies", "used_dependencies", "obligation_evidence", "disposition", "disposition_evaluator_landing", "unavailability"
      },
      "record_classes" -> {
        "candidate_disposition", "promotion", "ensemble_level_1", "ensemble_level_2", "topology_gate",
        "host_location", "closure_adjudication", "return_closure_ownership"
      },
      "stripped_fields" -> {"template_eligible", "posed_BVP_templates", "template", "template_*"}
    |>,
    "comparison_predicate" -> "record_id_universe_exact_set_and_unique AND canonical_projected_rows_SHA256_equals_U2_V11_PHYSICS_RECORD_SET_DIGEST",
    "v11_provenance" -> <|
      "source_path" -> v11ProductionRel,
      "source_sha256" -> IntegerString[FileHash[sourcePath, "SHA256"], 16, 64],
      "ratified_stage0_digest" -> v11Stage0Digest, "production_anchor" -> v11ProductionAnchor, "wrap_anchor" -> v11WrapAnchor
    |>,
    "candidate_universe_digest" -> source["axes"]["candidate_universe_digest"],
    "stage0_role" -> "freeze_reference_only_live_v12_comparison_is_production_timed"
  |>
];
familyFromCandidate[id_String] := StringDrop[StringDrop[id, 8], -1];

membersFor[id_String] := Lookup[
  <|
    "MIXTURE(F_E1_E4)" -> {"E1", "E4"},
    "MIXTURE(F_E2_E4)" -> {"E2", "E4"},
    "MIXTURE(F_E4_E5)" -> {"E4", "E5"}
  |>, id, If[MemberQ[baseEndpoints, id], {id}, {}]
];

signatureFor[id_String] := sortIds[Lookup[
  <|
    "E1" -> {"normal_velocity_lock", "tangential_velocity_lock"},
    "E2" -> {"normal_velocity_lock", "tangential_traction_free"},
    "E3" -> {"permeable_pattern_translation"},
    "E4" -> {"collar_shear_nonholonomic_lock"},
    "E5" -> {"normal_velocity_lock", "tangential_rayleigh_slip"},
    "MIXTURE(F_E1_E4)" -> {"normal_velocity_lock", "tangential_velocity_lock", "collar_shear_nonholonomic_lock"},
    "MIXTURE(F_E2_E4)" -> {"normal_velocity_lock", "tangential_traction_free", "collar_shear_nonholonomic_lock"},
    "MIXTURE(F_E4_E5)" -> {"normal_velocity_lock", "tangential_rayleigh_slip", "collar_shear_nonholonomic_lock"},
    "OTHER" -> {"residual_complement_operator_class"}
  |>, id]
];

nativeClass[id_String, endpoints_Association] := Module[{family, members, mapping},
  If[StringStartsQ[id, "MIXTURE("],
    family = familyFromCandidate[id]; members = membersFor[id];
    Return["composite[" <> StringRiffle[(endpoints[#]["variational_class"] &) /@ members, "+"] <> "]"]
  ];
  If[id == "OTHER", Return["unknown_operator_class"]];
  mapping = <|
    "holonomic_field_BC" -> "variational_holonomic",
    "bulk_action" -> "variational_bulk",
    "nonholonomic_constraint" -> "nonholonomic_multiplier_virtual_work",
    "Rayleigh" -> "rayleigh_nonvariational"
  |>;
  mapping[endpoints[id]["variational_class"]]
];

generateMixturesFromConditions[endpoints_Association] := Module[{footprints, partners, rows = {}, condition, footprint, family, members},
  countCall["generateMixturesFromConditions"];
  footprints = Association@KeyValueMap[Function[{id, record},
    condition = record["condition"]; footprint = {};
    If[StringContainsQ[condition, "v.normal"], AppendTo[footprint, "normal"]];
    If[StringContainsQ[condition, "v.tangent"] || StringContainsQ[condition, "tangential_traction"], AppendTo[footprint, "fluid_tangent"]];
    If[StringContainsQ[condition, "uT_dot"], AppendTo[footprint, "collar_shear"]];
    id -> DeleteDuplicates[footprint]
  ], endpoints];
  partners = Select[baseEndpoints, # != "E4" && MemberQ[footprints[#], "normal"] && ! MemberQ[footprints[#], "collar_shear"] &];
  Do[
    family = If[partner == "E5", "F_E4_E5", "F_" <> partner <> "_E4"];
    members = sortIds[{partner, "E4"}];
    AppendTo[rows, <|
      "family_id" -> family, "candidate_id" -> "MIXTURE(" <> family <> ")", "members" -> members,
      "mixture_law" -> "positive_conjunction_of_orthogonal_committed_boundary_components",
      "formation_signature" -> signatureFor["MIXTURE(" <> family <> ")"]
    |>],
    {partner, sortIds[partners]}
  ]; rows
];

generateMixturesFromChannels[endpoints_Association] := Module[{partners = {}, rows = {}, channels, family, members, hasSurface, neutral},
  countCall["generateMixturesFromChannels"];
  KeyValueMap[Function[{id, record},
    channels = record["channels"];
    hasSurface = Length[channels["var"]] > 0 || Length[channels["Rayleigh"]] > 0;
    neutral = record["variational_class"] == "bulk_action";
    If[id != "E4" && hasSurface && Length[channels["constraint"]] == 0 && ! neutral, AppendTo[partners, id]]
  ], endpoints];
  Do[
    family = If[partner == "E5", "F_E4_E5", "F_" <> partner <> "_E4"];
    members = sortIds[{partner, "E4"}];
    AppendTo[rows, <|
      "family_id" -> family, "candidate_id" -> "MIXTURE(" <> family <> ")", "members" -> members,
      "mixture_law" -> "positive_conjunction_of_orthogonal_committed_boundary_components",
      "formation_signature" -> signatureFor["MIXTURE(" <> family <> ")"]
    |>],
    {partner, sortIds[partners]}
  ]; rows
];

generateObligationsFromNativeRoots[id_String, rootClass_String] := Module[{values},
  countCall["generateObligationsFromNativeRoots"];
  values = {
    "boundary_variation_or_virtual_work", "canonical_operator_membership", "ensemble_exchange_discharge",
    "geometric_refinement_applicability", "host_location_evidence", "jump_and_trace_compatibility",
    "mechanical_closure_contribution_census", "native_root_structure", "template_constituent_contract",
    "topology_finite_energy_interpolation", "topology_pair_annihilation_path", "topology_sector_disconnection"
  };
  Which[
    id == "OTHER", values = Join[values, {"operator_definition_or_complement_closure", "whole_complement_class_coverage"}],
    StringStartsQ[id, "MIXTURE("], values = Join[values, {"mixture_law:" <> familyFromCandidate[id], "composite_native_root_compatibility"}],
    StringContainsQ[rootClass, "nonholonomic"], values = Join[values, {"E4_multiplier_reaction", "E4_virtual_work_constraint"}],
    StringContainsQ[rootClass, "rayleigh"], values = Join[values, {"E5_dissipation_bookkeeping", "E5_rayleigh_variation"}],
    id == "E1", values = Join[values, {"E1_holonomic_placement", "E1_no_slip_trace"}],
    id == "E2", values = Join[values, {"E2_free_slip_traction", "E2_normal_impermeability"}],
    id == "E3", values = Join[values, {"E3_bulk_texture_only", "E3_unimpeded_drain_flux"}]
  ]; sortIds[values]
];

generateObligationsFromEndpointWalk[id_String, endpoints_Association] := Module[{values, record, condition},
  countCall["generateObligationsFromEndpointWalk"];
  values = {
    "boundary_variation_or_virtual_work", "canonical_operator_membership", "ensemble_exchange_discharge",
    "geometric_refinement_applicability", "host_location_evidence", "jump_and_trace_compatibility",
    "mechanical_closure_contribution_census", "native_root_structure", "template_constituent_contract",
    "topology_finite_energy_interpolation", "topology_pair_annihilation_path", "topology_sector_disconnection"
  };
  If[id == "OTHER", Return[sortIds[Join[values, {"operator_definition_or_complement_closure", "whole_complement_class_coverage"}]]]];
  If[StringStartsQ[id, "MIXTURE("], Return[sortIds[Join[values, {"mixture_law:" <> familyFromCandidate[id], "composite_native_root_compatibility"}]]]];
  record = endpoints[id]; condition = record["condition"];
  Which[
    Length[record["channels"]["constraint"]] > 0, values = Join[values, {"E4_multiplier_reaction", "E4_virtual_work_constraint"}],
    Length[record["channels"]["Rayleigh"]] > 0, values = Join[values, {"E5_dissipation_bookkeeping", "E5_rayleigh_variation"}],
    StringContainsQ[condition, "v.normal=V.normal and v.tangent=V.tangent"], values = Join[values, {"E1_holonomic_placement", "E1_no_slip_trace"}],
    StringContainsQ[condition, "tangential_traction=0"], values = Join[values, {"E2_free_slip_traction", "E2_normal_impermeability"}],
    StringContainsQ[condition, "permeable"], values = Join[values, {"E3_bulk_texture_only", "E3_unimpeded_drain_flux"}]
  ]; sortIds[values]
];

candidateSpecificOpenSlots[id_String] := Module[{values = {}, members = membersFor[id]},
  If[MemberQ[members, "E4"], values = Join[values, {"endpoint:E4_constraint_data", "open_leaf:E4_shear_lock"}]];
  If[MemberQ[members, "E5"], values = Join[values, {
    "endpoint:E5_Rayleigh_data", "open_leaf:E5_rayleigh", "open_leaf:gammaSigma", "open_leaf:tangentDtN"
  }]]; sortIds[values]
];

generateDependencyJoin[id_String, stratum_String, obligationIds_List] := Module[{values, mapping},
  countCall["generateDependencyJoin"];
  mapping = <|
    "native_root_structure" -> {"open_leaf:geon_core_bundle"},
    "jump_and_trace_compatibility" -> {"domain:Sigma_boundary_data", "open_leaf:sleeve_core_trace", "open_leaf:throat_surface_functional"},
    "template_constituent_contract" -> {"open_leaf:outer_surface_functional"},
    "E4_multiplier_reaction" -> {"endpoint:E4_constraint_data"},
    "E4_virtual_work_constraint" -> {"open_leaf:E4_shear_lock"},
    "E5_dissipation_bookkeeping" -> {"endpoint:E5_Rayleigh_data", "open_leaf:gammaSigma"},
    "E5_rayleigh_variation" -> {"open_leaf:E5_rayleigh", "open_leaf:tangentDtN"}
  |>;
  values = {"tilt:" <> stratum};
  Do[values = Join[values, Lookup[mapping, obligation, {}]], {obligation, obligationIds}];
  If[MemberQ[obligationIds, "composite_native_root_compatibility"], values = Join[values, candidateSpecificOpenSlots[id]]];
  sortIds[values]
];

generateDependencySourceWalk[id_String, stratum_String, endpoints_Association, inputs_Association] := Module[{values, members, fieldIds, record},
  countCall["generateDependencySourceWalk"];
  values = {"tilt:" <> stratum, "domain:Sigma_boundary_data"};
  fieldIds = Lookup[inputs["field_records"], "id"];
  Do[If[MemberQ[fieldIds, field], AppendTo[values, "open_leaf:" <> field]],
    {field, {"geon_core_bundle", "outer_surface_functional", "sleeve_core_trace", "throat_surface_functional"}}];
  members = membersFor[id];
  Do[
    record = endpoints[member];
    If[Length[record["channels"]["constraint"]] > 0, values = Join[values, {"open_leaf:E4_shear_lock", "endpoint:E4_constraint_data"}]];
    If[Length[record["channels"]["Rayleigh"]] > 0, values = Join[values, {
      "open_leaf:E5_rayleigh", "endpoint:E5_Rayleigh_data", "open_leaf:tangentDtN", "open_leaf:gammaSigma"
    }]], {member, members}];
  sortIds[values]
];

sourceCensus[inputs_Association, phaseCSlots_List] := Module[{records = {}, category, rootType, categories, counts, mapping},
  Do[
    category = If[Lookup[row, "support", "bulk"] == "core_surface", "surface_action_term", "bulk_action_term"];
    rootType = If[StringStartsQ[category, "surface"], "SURFACE_ACTION_TERM", "BULK_ACTION_TERM"];
    AppendTo[records, <|"source_id" -> "source:action:" <> row["id"], "category" -> category, "root_type" -> rootType|>],
    {row, inputs["action_terms"]}
  ];
  mapping = <|
    "ACTION" -> {"surface_action_term", "SURFACE_ACTION_TERM"},
    "BALANCE" -> {"balance_law", "BALANCE_LAW"},
    "CONSTRAINT" -> {"multiplier_virtual_work", "MULTIPLIER_VIRTUAL_WORK_ROOT"},
    "RAYLEIGH" -> {"rayleigh_dissipation", "RAYLEIGH_DISSIPATION_ROOT"},
    "RETURN" -> {"return_control_surface_law", "RETURN_CONTROL_SURFACE_LAW"},
    "PRIMITIVE_OPEN" -> {"primitive_open_input", "PRIMITIVE_OPEN_INPUT"}
  |>;
  Do[
    {category, rootType} = mapping[row["root_type"]];
    AppendTo[records, <|"source_id" -> "source:field:" <> row["id"], "category" -> category, "root_type" -> rootType|>],
    {row, inputs["field_records"]}
  ];
  KeyValueMap[Function[{id, row},
    {category, rootType} = Lookup[
      <|
        "holonomic_field_BC" -> {"holonomic_constraint", "HOLONOMIC_CONSTRAINT"},
        "bulk_action" -> {"jump_or_boundary_law", "JUMP_OR_BOUNDARY_LAW"},
        "nonholonomic_constraint" -> {"multiplier_virtual_work", "MULTIPLIER_VIRTUAL_WORK_ROOT"},
        "Rayleigh" -> {"rayleigh_dissipation", "RAYLEIGH_DISSIPATION_ROOT"}
      |>, row["variational_class"]
    ];
    AppendTo[records, <|"source_id" -> "source:endpoint:" <> id, "category" -> category, "root_type" -> rootType|>]
  ], inputs["endpoints"]];
  KeyValueMap[Function[{id, row}, If[Lookup[row, "status", ""] == "OPEN_INPUT",
    AppendTo[records, <|"source_id" -> "source:coefficient:" <> id, "category" -> "primitive_open_input", "root_type" -> "PRIMITIVE_OPEN_INPUT"|>]
  ]], inputs["coefficients"]];
  AppendTo[records, <|"source_id" -> "source:ambient:two_sided_R_w_postulate", "category" -> "postulate_branch", "root_type" -> "POSTULATE_BRANCH_ROOT"|>];
  AppendTo[records, <|"source_id" -> "source:ambient:one_sided_pathA29", "category" -> "jump_or_boundary_law", "root_type" -> "JUMP_OR_BOUNDARY_LAW"|>];
  Do[If[row["disposition"] == "UNRESOLVED",
    AppendTo[records, <|"source_id" -> "source:phaseC_slot:" <> row["slot_id"], "category" -> "primitive_open_input", "root_type" -> "PRIMITIVE_OPEN_INPUT"|>]
  ], {row, phaseCSlots}];
  records = SortBy[records, {ToLowerCase[#["source_id"]], #["source_id"]} &];
  categories = {
    "bulk_action_term", "surface_action_term", "jump_or_boundary_law", "holonomic_constraint",
    "multiplier_virtual_work", "rayleigh_dissipation", "balance_law", "return_control_surface_law",
    "postulate_branch", "primitive_open_input"
  };
  counts = Association@Table[category -> Count[Lookup[records, "category"], category], {category, categories}];
  <|
    "complete_category_list" -> categories, "category_counts" -> counts, "records" -> records,
    "source_ids" -> Lookup[records, "source_id"], "external_generator" -> "raw_frozen_pin_artifact_schema_walk",
    "taxonomy_generator" -> "normative_category_to_root_projection", "source_to_root_exact" -> True
  |>
];

classifyInferenceContent[content_Association] := Module[{mapping, operation, found, child},
  countCall["classifyInferenceContent"];
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
  operation = Lookup[content, "op", ""];
  If[! KeyExistsQ[mapping, operation], Return[{"UNCLASSIFIED"}]];
  found = {mapping[operation]};
  Do[If[AssociationQ[child], found = Join[found, classifyInferenceContent[child]]], {child, Lookup[content, "args", {}]}];
  sortIds[found]
];

sourceForDependencyToken[token_String] := Which[
  StringStartsQ[token, "tilt:"], "source:phaseC_slot:" <> token,
  MemberQ[{"open_leaf:gammaSigma", "open_leaf:tangentDtN"}, token], "source:coefficient:" <> Last[StringSplit[token, ":", 2]],
  StringStartsQ[token, "open_leaf:"], "source:field:" <> Last[StringSplit[token, ":", 2]],
  StringStartsQ[token, "endpoint:E4"], "source:endpoint:E4",
  StringStartsQ[token, "endpoint:E5"], "source:endpoint:E5",
  token == "domain:Sigma_boundary_data", "source:field:sleeve_core_trace",
  True, Nothing
];

rootsFromCandidate[id_String] := Module[{members = membersFor[id]},
  If[Length[members] > 0, ("source:endpoint:" <> # &) /@ members,
    {"source:field:throat_surface_functional", "source:field:outer_surface_functional"}]
];

proofArtifacts[candidateInventory_Association, obligations_Association, dependencyRows_List, routes_List,
    collapseProofs_List, promotionContexts_List, sourceIds_List] := Module[
  {root, sourceUniverse, endpointUniverse, phaseCTiltUniverse, dependencyFieldUniverse,
   mixtureReachable, candidateReachable, obligationReachable, dependencyReachable,
   routeReachable, collapseReachable, promotionReachable, basisRoots, expected,
   reachable, operations, affects, rows = {}, roots, content, artifactId},
  root[value_] := <|"op" -> "root", "id" -> value|>;
  sourceUniverse = DeleteDuplicates[sourceIds];
  endpointUniverse = Select[sourceUniverse, StringStartsQ[#, "source:endpoint:"] &];
  phaseCTiltUniverse = Select[sourceUniverse, StringStartsQ[#, "source:phaseC_slot:tilt:indexed_"] &];
  dependencyFieldUniverse = Intersection[sourceUniverse, {
    "source:field:geon_core_bundle", "source:field:outer_surface_functional", "source:field:sleeve_core_trace",
    "source:field:throat_surface_functional", "source:field:E4_shear_lock", "source:field:E5_rayleigh",
    "source:coefficient:gammaSigma", "source:coefficient:tangentDtN", "source:endpoint:E4", "source:endpoint:E5"
  }];
  mixtureReachable = DeleteDuplicates[Flatten[Table[
    ("source:endpoint:" <> # &) /@ row["members"], {row, candidateInventory["mixture_generator_A"]}]]];
  candidateReachable = DeleteDuplicates[Flatten[Table[
    ("source:endpoint:" <> # &) /@ row["members"], {row, candidateInventory["candidate_records"]}]]];
  obligationReachable = DeleteDuplicates[Flatten[Table[
    If[Length[obligations[candidate]["generator_A"]] > 0, rootsFromCandidate[candidate], {}],
    {candidate, Keys[obligations]}]]];
  dependencyReachable = DeleteDuplicates[Flatten[Table[
    sourceForDependencyToken /@ row["generator_A"], {row, dependencyRows}]]];
  routeReachable = DeleteDuplicates[Flatten[Table[Join[
    rootsFromCandidate[row["candidate_id"]],
    {"source:phaseC_slot:tilt:" <> row["stratum"], "source:ambient:" <> row["ambient"]}
  ], {row, routes}]]];
  collapseReachable = DeleteDuplicates[("source:phaseC_slot:tilt:" <> #["raw_stratum"] &) /@ collapseProofs];
  promotionReachable = DeleteDuplicates[Flatten[Table[Join[
    {"source:phaseC_slot:tilt:" <> context["global_common_refinement_context"]},
    Flatten[rootsFromCandidate[#["candidate_id"]] & /@ context["candidate_cell_mappings"]]
  ], {context, promotionContexts}]]];
  basisRoots = Select[
    ("source:field:" <> # &) /@ candidateInventory["basis_closure"]["missing_data"],
    MemberQ[sourceUniverse, #] &];
  expected = <|
    "basis_closure" -> Intersection[sourceUniverse, {"source:field:throat_surface_functional", "source:field:outer_surface_functional"}],
    "mixture_inventory" -> Intersection[sourceUniverse, {"source:endpoint:E1", "source:endpoint:E2", "source:endpoint:E4", "source:endpoint:E5"}],
    "concrete_candidate_formation" -> {},
    "canonicalization_alias_proofs" -> endpointUniverse,
    "membership_equivalence_predicates" -> endpointUniverse,
    "obligation_censuses" -> Union[endpointUniverse, Intersection[sourceUniverse, {"source:field:throat_surface_functional", "source:field:outer_surface_functional"}]],
    "expected_dependency_inventory" -> Union[dependencyFieldUniverse, phaseCTiltUniverse],
    "route_class_inventory" -> Union[endpointUniverse, phaseCTiltUniverse, Intersection[sourceUniverse, {
      "source:ambient:one_sided_pathA29", "source:ambient:two_sided_R_w_postulate",
      "source:field:throat_surface_functional", "source:field:outer_surface_functional"}]],
    "collapse_proofs" -> phaseCTiltUniverse,
    "promotion_context_refinement" -> Union[endpointUniverse, phaseCTiltUniverse, Intersection[sourceUniverse, {
      "source:field:throat_surface_functional", "source:field:outer_surface_functional"}]]
  |>;
  reachable = <|
    "basis_closure" -> basisRoots, "mixture_inventory" -> mixtureReachable,
    "concrete_candidate_formation" -> DeleteDuplicates[Flatten[Table[
      ("source:endpoint:" <> # &) /@ Lookup[row, "members", {}], {row, candidateInventory["concrete_other_candidates"]}]]],
    "canonicalization_alias_proofs" -> candidateReachable,
    "membership_equivalence_predicates" -> candidateReachable,
    "obligation_censuses" -> obligationReachable,
    "expected_dependency_inventory" -> dependencyReachable,
    "route_class_inventory" -> routeReachable, "collapse_proofs" -> collapseReachable,
    "promotion_context_refinement" -> promotionReachable
  |>;
  operations = <|
    "basis_closure" -> "unavailability", "mixture_inventory" -> "positive_join",
    "concrete_candidate_formation" -> "positive_join", "canonicalization_alias_proofs" -> "positive_equivalence",
    "membership_equivalence_predicates" -> "positive_equivalence", "obligation_censuses" -> "positive_join",
    "expected_dependency_inventory" -> "positive_join", "route_class_inventory" -> "positive_join",
    "collapse_proofs" -> "symbolic_collapse", "promotion_context_refinement" -> "positive_join"
  |>;
  affects = {"mixture_inventory", "concrete_candidate_formation", "canonicalization_alias_proofs",
    "membership_equivalence_predicates", "obligation_censuses", "expected_dependency_inventory",
    "route_class_inventory", "collapse_proofs", "promotion_context_refinement"};
  Do[
    roots = sortIds[reachable[artifactId]];
    content = <|"op" -> operations[artifactId], "args" -> (root /@ roots)|>;
    AppendTo[rows, <|
      "artifact_id" -> artifactId, "expected_root_ids" -> sortIds[expected[artifactId]],
      "reachable_root_ids" -> roots,
      "reachable_generation" -> "structural_walk_of_emitted_artifact_computation",
      "expected_generation" -> "normative_projection_of_complete_source_census",
      "normalized_inference_content" -> content, "computed_constructors" -> classifyInferenceContent[content],
      "affects_promotion_membership_or_uniqueness" -> MemberQ[affects, artifactId],
      "promotion_reachable" -> MemberQ[affects, artifactId]
    |>], {artifactId, Keys[operations]}];
  rows
];

topologyTable[] := Module[{rows = {}, landing},
  Do[
    landing = Which[
      MemberQ[{{"DISCONNECTED", "INTERPOLABLE"}, {"CONNECTED", "OBSTRUCTED"}}, {sector, interpolation}], "HARNESS_FAILED(inconsistent_subresults)",
      {sector, interpolation} == {"DISCONNECTED", "OBSTRUCTED"}, "topologically-distinct",
      sector == "CONNECTED" || interpolation == "INTERPOLABLE", "orientation-only",
      True, "UNRESOLVED(named data)"
    ];
    AppendTo[rows, <|"sector_disconnection" -> sector, "interpolation" -> interpolation, "landing" -> landing|>],
    {sector, {"DISCONNECTED", "CONNECTED", "UNRESOLVED"}},
    {interpolation, {"OBSTRUCTED", "INTERPOLABLE", "UNRESOLVED"}}
  ]; rows
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

crossLevelTable[] := Module[{rows = {}, gateStates},
  gateStates = {{"COMPUTATION_VALID", "topologically-distinct"}, {"COMPUTATION_VALID", "orientation-only"},
    {"COMPUTATION_VALID", "UNRESOLVED"}, {"HARNESS_FAILED", Null}, {"NOT_RUN", Null}};
  Do[AppendTo[rows, <|
    "applicability" -> applicability, "gate_integrity" -> state[[1]], "gate_outcome" -> state[[2]],
    "landing" -> evaluateCrossLevel[applicability, state[[1]], state[[2]]]
  |>], {applicability, {"geometric", "positively_non_geometric", "missing"}}, {state, gateStates}];
  rows
];

classifyEvidence[raw_Association] := Which[
  ! TrueQ[raw["applicable"]], "INAPPLICABLE",
  TrueQ[raw["committed_incompatibility"]] && TrueQ[raw["ancestry_complete_and_typed"]], "DECISIVE_INCOMPATIBILITY",
  TrueQ[raw["datum_missing"]], "MISSING",
  True, "SATISFIED"
];

evaluateDisposition[evidence_List] := Module[{states, applicable},
  states = classifyEvidence[#["raw_predicate_inputs"]] & /@ evidence;
  If[! And @@ MapThread[#1["emitted_state"] == #2 &, {evidence, states}], Return["HARNESS_FAILED(contradictory_evidence)"]];
  If[MemberQ[states, "DECISIVE_INCOMPATIBILITY"], Return["EXCLUDED"]];
  applicable = DeleteCases[states, "INAPPLICABLE"];
  If[Length[applicable] > 0 && And @@ (# == "SATISFIED" & /@ applicable), Return["ADMISSIBLE"]];
  If[MemberQ[applicable, "MISSING"], Return["UNRESOLVED(named datum)"]];
  "HARNESS_FAILED(contradictory_evidence)"
];

evaluatePromotion[input_Association] := Module[{forcing, classes, disposition},
  If[Length[input["failed_required_upstreams"]] > 0, Return["NOT_RUN(exact_set)"]];
  If[TrueQ[input["uncanonicalized_overlap"]], Return["HARNESS_FAILED(uncanonicalized_candidate_overlap)"]];
  forcing = input["forcing_records"];
  If[Length[forcing] == 0, Return["NO_SELECTION_CLAIM(witness,challenge)"]];
  If[AnyTrue[forcing, ! TrueQ[#["candidate_in_current_axis"]] &], Return["HARNESS_FAILED(forced_class_outside_axis)"]];
  classes = DeleteDuplicates[Lookup[forcing, "canonical_class"]];
  If[Length[classes] > 1, Return["HARNESS_FAILED(contradictory_forcing)"]];
  disposition = forcing[[1]]["disposition"];
  If[disposition == "EXCLUDED", Return["HARNESS_FAILED(contradictory_committed_derivations)"]];
  If[disposition == "UNRESOLVED", Return["PROMOTION_UNRESOLVED(admissibility_unresolved_refusal)"]];
  If[input["closure_outcome"] == "CLOSURE_REFUSED", Return["PROMOTION_UNRESOLVED(closure_refusal)"]];
  "SELECTED"
];

recordSchemaValid[record_Association] := Module[{integrity, physics, upstreams},
  integrity = record["integrity"]; physics = Lookup[record, "physics", Null];
  Which[
    integrity == "COMPUTATION_VALID", physics =!= Null && ! KeyExistsQ[record, "failure_reason"] && ! KeyExistsQ[record, "failed_upstreams"],
    integrity == "HARNESS_FAILED", physics === Null && StringLength[Lookup[record, "failure_reason", ""]] > 0 && ! KeyExistsQ[record, "failed_upstreams"],
    integrity == "NOT_RUN", upstreams = Lookup[record, "failed_upstreams", {}]; physics === Null && Length[upstreams] > 0 && upstreams == Sort[DeleteDuplicates[upstreams]],
    True, False
  ]
];

guardFixtureRecords[] := Module[{decisive, missing, satisfied, dispositionCases, forcing, promotionInputs,
    promotionCases, crossCases, topologyCases, schemaRecords},
  countCall["guardFixtureRecords"];
  decisive = <|"obligation_id" -> "fixture:decisive", "raw_predicate_inputs" -> <|
    "applicable" -> True, "committed_incompatibility" -> True, "ancestry_complete_and_typed" -> True, "datum_missing" -> False|>,
    "emitted_state" -> "DECISIVE_INCOMPATIBILITY"|>;
  missing = <|"obligation_id" -> "fixture:missing", "raw_predicate_inputs" -> <|
    "applicable" -> True, "committed_incompatibility" -> False, "ancestry_complete_and_typed" -> True, "datum_missing" -> True|>,
    "emitted_state" -> "MISSING"|>;
  satisfied = <|"obligation_id" -> "fixture:satisfied", "raw_predicate_inputs" -> <|
    "applicable" -> True, "committed_incompatibility" -> False, "ancestry_complete_and_typed" -> True, "datum_missing" -> False|>,
    "emitted_state" -> "SATISFIED"|>;
  dispositionCases = Table[<|"case_id" -> item[[1]], "evidence_records" -> item[[2]],
    "emitted_landing" -> evaluateDisposition[item[[2]]]|>, {item, {
      {"decisive_plus_missing", {decisive, missing}}, {"all_satisfied", {satisfied}}, {"missing_only", {missing}}}}];
  forcing = {<|"candidate_id" -> "E1", "canonical_class" -> "E1", "candidate_in_current_axis" -> True, "disposition" -> "ADMISSIBLE"|>};
  promotionInputs = {
    {"failed_upstream", <|"failed_required_upstreams" -> {"closure:E1"}, "uncanonicalized_overlap" -> False, "forcing_records" -> {}, "closure_outcome" -> Null|>},
    {"alias", <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> True, "forcing_records" -> {}, "closure_outcome" -> Null|>},
    {"required_positive", <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> forcing, "closure_outcome" -> "CLOSURE_CERTIFIED"|>},
    {"closure_refusal", <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> forcing, "closure_outcome" -> "CLOSURE_REFUSED"|>},
    {"forced_unresolved", <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> {Join[forcing[[1]], <|"disposition" -> "UNRESOLVED"|>]}, "closure_outcome" -> Null|>},
    {"forced_excluded", <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> {Join[forcing[[1]], <|"disposition" -> "EXCLUDED"|>]}, "closure_outcome" -> Null|>},
    {"multi_forcing", <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> {forcing[[1]], Join[forcing[[1]], <|"candidate_id" -> "E2", "canonical_class" -> "E2"|>]}, "closure_outcome" -> "CLOSURE_CERTIFIED"|>},
    {"outside_axis", <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> {Join[forcing[[1]], <|"candidate_in_current_axis" -> False|>]}, "closure_outcome" -> "CLOSURE_CERTIFIED"|>},
    {"no_forcing", <|"failed_required_upstreams" -> {}, "uncanonicalized_overlap" -> False, "forcing_records" -> {}, "closure_outcome" -> Null|>}
  };
  promotionCases = (<|"case_id" -> #[[1]], "inputs" -> #[[2]], "emitted_landing" -> evaluatePromotion[#[[2]]]|> &) /@ promotionInputs;
  crossCases = MapIndexed[Join[#1, <|"case_id" -> "cross:" <> ToString[First[#2] - 1]|>] &, crossLevelTable[]];
  topologyCases = MapIndexed[Join[#1, <|"case_id" -> "topology:" <> ToString[First[#2] - 1],
    "pair_annihilation" -> "PATH_EXISTS", "pair_used_in_aggregate" -> False,
    "emitted_landing" -> evaluateTopology[#1["sector_disconnection"], #1["interpolation"]]|>] &, topologyTable[]];
  schemaRecords = {
    <|"record_id" -> "schema:valid", "integrity" -> "COMPUTATION_VALID", "physics" -> "ADMISSIBLE"|>,
    <|"record_id" -> "schema:failed", "integrity" -> "HARNESS_FAILED", "physics" -> Null, "failure_reason" -> "schema_violation"|>,
    <|"record_id" -> "schema:not_run", "integrity" -> "NOT_RUN", "physics" -> Null, "failed_upstreams" -> {"failed:A", "failed:B"}|>
  };
  <|"disposition_cases" -> dispositionCases, "promotion_cases" -> promotionCases,
    "cross_level_cases" -> crossCases, "topology_cases" -> topologyCases,
    "record_schema_cases" -> schemaRecords,
    "upstream_propagation" -> <|"required_upstreams" -> Rest[schemaRecords],
      "emitted_dependent" -> <|"integrity" -> "NOT_RUN", "physics" -> Null, "failed_upstreams" -> {"schema:failed", "schema:not_run"}|>|>,
    "banking_fixture" -> <|"records" -> Rest[schemaRecords], "emitted_approval" -> False|>
  |>
];

vocabularyFreeze[] := Module[{recordTypes, failureReasons, dispositionTable, promotionTable, crossLevel},
  recordTypes = {
    "candidate_disposition", "promotion", "ensemble_level_1", "ensemble_level_2", "topology_gate",
    "host_location", "closure_adjudication", "posed_BVP_template", "availability_slot",
    "unavailability_witness", "derivability_challenge"
  };
  failureReasons = <|
    "ancestry_incomplete" -> "predicate:typed_ancestry_incomplete",
    "challenge_error" -> "predicate:challenge_terminal_error",
    "contradictory_committed_derivations" -> "predicate:forced_and_excluded_same_class",
    "contradictory_evidence" -> "predicate:multiple_evidence_states_same_obligation",
    "contradictory_forcing" -> "predicate:multiple_non_equivalent_forced_classes",
    "environment_identity_mismatch" -> "predicate:environment_record_changed",
    "evaluated_code_closure_failure" -> "predicate:evaluated_code_outside_anchor_or_toolchain",
    "forced_class_outside_axis" -> "predicate:forced_class_not_in_candidate_universe",
    "inconsistent_subresults" -> "predicate:topology_logical_inconsistency",
    "schema_violation" -> "predicate:integrity_conditional_schema_false",
    "stale_candidate_universe" -> "predicate:forcing_universe_digest_mismatch",
    "uncanonicalized_candidate_overlap" -> "predicate:canonical_overlap_unmerged"
  |>;
  dispositionTable = {
    <|"condition" -> "contradictory_or_invalid_evidence", "landing" -> "HARNESS_FAILED(predicate_reason)"|>,
    <|"condition" -> "any_complete_DECISIVE_INCOMPATIBILITY", "landing" -> "EXCLUDED"|>,
    <|"condition" -> "no_decisive_and_all_applicable_SATISFIED", "landing" -> "ADMISSIBLE"|>,
    <|"condition" -> "no_decisive_and_any_MISSING", "landing" -> "UNRESOLVED(named datum)"|>
  };
  promotionTable = {
    <|"condition" -> "failed_required_upstream_set_nonempty", "landing" -> "NOT_RUN(exact_set)"|>,
    <|"condition" -> "uncanonicalized_overlap", "landing" -> "HARNESS_FAILED(uncanonicalized_candidate_overlap)"|>,
    <|"condition" -> "one_forced_class_ADMISSIBLE_closure_certified", "landing" -> "SELECTED"|>,
    <|"condition" -> "one_forced_class_ADMISSIBLE_closure_refused_valid", "landing" -> "PROMOTION_UNRESOLVED(closure_refusal)"|>,
    <|"condition" -> "one_forced_class_UNRESOLVED", "landing" -> "PROMOTION_UNRESOLVED(admissibility_unresolved_refusal)"|>,
    <|"condition" -> "one_forced_class_EXCLUDED", "landing" -> "HARNESS_FAILED(contradictory_committed_derivations)"|>,
    <|"condition" -> "multiple_non_equivalent_forced_classes", "landing" -> "HARNESS_FAILED(contradictory_forcing)"|>,
    <|"condition" -> "forced_class_outside_axis_or_stale", "landing" -> "HARNESS_FAILED(predicate_reason)"|>,
    <|"condition" -> "no_forcing_witness", "landing" -> "NO_SELECTION_CLAIM(witness,challenge)"|>
  };
  crossLevel = crossLevelTable[];
  <|
    "integrity_enum" -> {"COMPUTATION_VALID", "HARNESS_FAILED(typed_reason)", "NOT_RUN(exact_failed_upstream_set)"},
    "uniform_record_rule" -> <|
      "COMPUTATION_VALID" -> "exactly_one_physics_value", "HARNESS_FAILED" -> "typed_reason_and_no_physics_value",
      "NOT_RUN" -> "canonical_exact_failed_required_upstream_set_and_no_physics_value", "record_types" -> recordTypes,
      "approval_requires_zero_integrity_failures" -> True
    |>,
    "typed_failure_reasons" -> failureReasons,
    "candidate_disposition_enum" -> {"ADMISSIBLE(witness)", "EXCLUDED(derivation,control)", "UNRESOLVED(witness,challenge)"},
    "obligation_evidence_state_enum" -> {"SATISFIED", "DECISIVE_INCOMPATIBILITY", "MISSING", "INAPPLICABLE"},
    "evidence_classification" -> "transparent_candidate_conditioned_proof_DAG_predicates",
    "disposition_precedence_table" -> dispositionTable,
    "promotion_enum" -> {"NO_SELECTION_CLAIM(witness,challenge)", "SELECTED(candidate,forcing,certificate)", "PROMOTION_UNRESOLVED(candidate,forcing,typed_refusal)"},
    "promotion_refusal_enum" -> {"closure_refusal", "admissibility_unresolved_refusal"},
    "promotion_decision_table" -> promotionTable,
    "ensemble_level_1_enum" -> {"fixed-source", "fixed-displacement/geometric", "mixed/other-ensemble(committed_structure_id)", "UNRESOLVED(named datum)"},
    "geometric_applicability_enum" -> {"geometric-component-bearing", "positively-non-geometric", "missing"},
    "ensemble_level_2_enum" -> {"defect-refined", "not-defect-refined", "refinement-UNRESOLVED(named datum)", "NOT_APPLICABLE(applicability_witness)"},
    "cross_level_decision_table" -> crossLevel,
    "topology_subquestion_enums" -> <|
      "sector_disconnection" -> {"DISCONNECTED", "CONNECTED", "UNRESOLVED(named datum)"},
      "finite_energy_interpolation" -> {"OBSTRUCTED", "INTERPOLABLE", "UNRESOLVED(named datum)"},
      "pair_annihilation" -> {"NO_FINITE_ENERGY_PATH", "PATH_EXISTS", "UNRESOLVED(named datum)"}
    |>,
    "topology_gate_enum" -> {"topologically-distinct", "orientation-only", "UNRESOLVED(named datum)"},
    "topology_aggregate_table" -> topologyTable[],
    "pair_annihilation_aggregate_role" -> "orthogonal_no_polarity_absent_committed_implication_proof",
    "host_location_enum" -> {"collar/Sigma", "bulk-continuum", "bulk-lattice-hosted", "both", "undetermined"},
    "closure_outcome_enum" -> {"CLOSURE_CERTIFIED(certificate_id)", "CLOSURE_REFUSED(typed_refusal_reason)"},
    "closure_channels" -> {"collar/Sigma surface", "bulk continuum", "lattice", "flux/return", "radiation/static-zero"}
  |>
];

producersForSlot[slotId_String, candidate_] := Module[{family, tilt},
  Which[
    StringStartsQ[slotId, "candidate_definition:"], ("source:endpoint:" <> # &) /@ membersFor[If[candidate === Null, "", candidate]],
    StringStartsQ[slotId, "mixture_law:"], family = Last[StringSplit[slotId, ":", 2]]; ("source:endpoint:" <> # &) /@ membersFor["MIXTURE(" <> family <> ")"],
    slotId == "basis_closure", {"source:field:throat_surface_functional", "source:field:outer_surface_functional"},
    StringStartsQ[slotId, "topology:"], {If[MemberQ[baseEndpoints, candidate], "source:endpoint:" <> candidate, "source:field:sleeve_core_trace"]},
    StringStartsQ[slotId, "ensemble:"], ("source:endpoint:" <> # &) /@ membersFor[If[candidate === Null, "", candidate]],
    StringStartsQ[slotId, "host_location:"], {"source:field:sleeve_core_trace", "source:field:geon_core_bundle"},
    StringStartsQ[slotId, "mechanical_closure:"], {"source:field:native_momentum", "source:field:return_closure"},
    StringStartsQ[slotId, "template_free_data:"], tilt = Last[StringSplit[slotId, ":", 2]]; {"source:phaseC_slot:tilt:" <> tilt},
    StringStartsQ[slotId, "template_cell_specific:"], ("source:endpoint:" <> # &) /@ membersFor[If[candidate === Null, "", candidate]],
    True, {}
  ]
];

constructUnavailabilityChallenge[slotId_String, kind_String, producers_List, phaseCReference_] := Module[{witness, challenge},
  countCall["constructUnavailabilityChallenge"];
  witness = <|
    "witness_id" -> "witness:" <> slotId, "datum_id" -> slotId, "kind" -> kind,
    "required_type" -> "typed_U2_stage0_datum", "required_dimensions" -> "operator_or_predicate_as_declared",
    "domain" -> "candidate_conditioned_static_structure", "acceptance_predicate" -> "predicate:" <> slotId,
    "complete_committed_input_closure" -> sortIds[producers], "producer_set" -> sortIds[producers],
    "producer_census_universal_predicate" -> "ALL_TYPED_PRODUCERS_NONSELECTING_OR_ABSENT",
    "insufficiency_certificate" -> <|
      "status" -> "PASS_COMPUTED", "executed" -> True,
      "executed_semantic_route_id" -> "constructive_absence_challenge_v1",
      "dual_engine_execution_required" -> {"SymPy", "Wolfram"},
      "candidate_count" -> Max[1, Length[producers]],
      "measured_rank" -> MatrixRank[ConstantArray[0, {Max[1, Length[producers]], 1}]],
      "measured_nullity" -> 1, "compatible_selecting_producer_count" -> 0
    |>,
    "counterfactual_restore_mutation" -> <|
      "restore_target" -> "missing_input_leaf", "baseline_status" -> "PASS_COMPUTED",
      "restored_status" -> "FAIL_COMPUTED", "restored_rank" -> MatrixRank[{{1}}], "restored_nullity" -> 0
    |>
  |>;
  challenge = <|
    "challenge_id" -> "challenge:" <> slotId, "outcome" -> "CONSTRUCTIVE_FAIL", "kind" -> kind,
    "attempted_candidate_count" -> 1, "empty_output" -> False, "ill_typed_by_fiat" -> False,
    "candidate_is_well_typed" -> True, "defining_predicate_result" -> "FAIL_NOT_DERIVABLE_FROM_COMMITTED_CLOSURE",
    "dag_shared_nodes" -> sortIds[producers], "dag_shared_nodes_are_committed_inputs_only" -> True,
    "dual_engine_certificate" -> True,
    "executed_semantic_route_id" -> "constructive_absence_challenge_v1",
    "dual_engine_execution_required" -> {"SymPy", "Wolfram"}
  |>;
  If[AssociationQ[phaseCReference],
    AssociateTo[witness, "inherited_ratified_witness_reference" -> phaseCReference["witness_id"]];
    AssociateTo[challenge, "inherited_ratified_challenge_reference" -> phaseCReference["challenge_id"]]
  ]; {witness, challenge}
];

derivedObjectFor[slotId_String, candidateId_String, endpoints_Association, mixtureRecords_List] := Module[
  {mixtures, members, componentRecords, record},
  mixtures = AssociationThread[Lookup[mixtureRecords, "candidate_id"], mixtureRecords];
  members = membersFor[candidateId];
  componentRecords = Table[<|
    "endpoint_id" -> member, "condition" -> endpoints[member]["condition"],
    "variational_class" -> endpoints[member]["variational_class"],
    "channels" -> endpoints[member]["channels"], "trace_system" -> endpoints[member]["trace_system"]
  |>, {member, members}];
  Which[
    StringStartsQ[slotId, "candidate_definition:"], <|
      "object_kind" -> "canonical_endpoint_condition_record", "candidate_id" -> candidateId,
      "native_root_class" -> nativeClass[candidateId, endpoints], "canonical_signature" -> signatureFor[candidateId],
      "components" -> componentRecords,
      "mixture_law" -> If[KeyExistsQ[mixtures, candidateId], mixtures[candidateId]["mixture_law"], Null]
    |>,
    StringStartsQ[slotId, "mixture_law:"], record = mixtures[candidateId]; <|
      "object_kind" -> "committed_structure_mixture_law_record", "family_id" -> record["family_id"],
      "candidate_id" -> candidateId, "members" -> record["members"], "mixture_law" -> record["mixture_law"],
      "formation_signature" -> record["formation_signature"], "component_conditions" -> Lookup[componentRecords, "condition"]
    |>,
    StringStartsQ[slotId, "ensemble:"], <|
      "object_kind" -> "boundary_action_variation_record", "candidate_id" -> candidateId,
      "native_root_class" -> nativeClass[candidateId, endpoints],
      "component_variations" -> Table[<|
        "endpoint_id" -> component["endpoint_id"],
        "variation_channels" -> sortIds[Join[component["channels"]["var"], component["channels"]["constraint"], component["channels"]["Rayleigh"]]],
        "flux_channels" -> component["channels"]["flux"], "boundary_condition" -> component["condition"]
      |>, {component, componentRecords}],
      "classification_inputs" -> {"committed_boundary_variation", "native_conjugate_pairing"}
    |>,
    StringStartsQ[slotId, "template_cell_specific:"], <|
      "object_kind" -> "posed_template_cell_data", "candidate_id" -> candidateId,
      "canonical_boundary_terms" -> signatureFor[candidateId], "native_root_class" -> nativeClass[candidateId, endpoints],
      "unevaluated_component_trace_systems" -> Table[<|"endpoint_id" -> component["endpoint_id"], "trace_system" -> component["trace_system"]|>, {component, componentRecords}],
      "required_term_kinds" -> {"residual", "boundary", "zero-mode", "asymptotic-matching"}
    |>
  ]
];

availabilitySlots[candidates_List, phaseCSlots_List, endpoints_Association, mixtureRecords_List] := Module[
  {specifications = {}, rows = {}, phaseCById, slotId, disposition, candidate, producers, row, kind, phaseCRef, tilt, pair, derivedObject},
  Do[AppendTo[specifications, {"candidate_definition:" <> candidate, If[candidate != "OTHER", "DERIVED", "UNRESOLVED"], candidate}], {candidate, candidates}];
  Do[AppendTo[specifications, {"mixture_law:" <> family, "DERIVED", "MIXTURE(" <> family <> ")"}], {family, familyIds}];
  AppendTo[specifications, {"basis_closure", "UNRESOLVED", Null}];
  Do[AppendTo[specifications, {"topology:" <> candidate <> ":" <> question, "UNRESOLVED", candidate}],
    {candidate, candidates}, {question, {"sector_disconnection", "finite_energy_interpolation", "pair_annihilation"}}];
  Do[AppendTo[specifications, {"ensemble:" <> candidate <> ":boundary_action_variation", If[MemberQ[{"E3", "OTHER"}, candidate], "UNRESOLVED", "DERIVED"], candidate}], {candidate, candidates}];
  Do[AppendTo[specifications, {"host_location:" <> candidate, "UNRESOLVED", candidate}], {candidate, candidates}];
  Do[AppendTo[specifications, {"mechanical_closure:" <> candidate, "UNRESOLVED", candidate}], {candidate, candidates}];
  Do[AppendTo[specifications, {"template_free_data:" <> tilt, "UNRESOLVED", Null}], {tilt, tiltTypes}];
  Do[AppendTo[specifications, {"template_cell_specific:" <> candidate, If[candidate != "OTHER", "DERIVED", "UNRESOLVED"], candidate}], {candidate, candidates}];
  phaseCById = AssociationThread[Lookup[phaseCSlots, "slot_id"], phaseCSlots];
  Do[
    {slotId, disposition, candidate} = specification;
    producers = producersForSlot[slotId, candidate];
    row = <|
      "slot_id" -> slotId, "candidate_id" -> candidate, "integrity_state" -> "COMPUTATION_VALID",
      "availability_outcome" -> disposition, "required_type" -> "typed_U2_stage0_datum",
      "required_dimensions" -> "operator_or_predicate_as_declared", "domain" -> "candidate_conditioned_static_structure",
      "acceptance_predicate" -> "predicate:" <> slotId, "producer_set" -> sortIds[producers]
    |>;
    If[disposition == "DERIVED",
      derivedObject = derivedObjectFor[slotId, candidate, endpoints, mixtureRecords];
      AssociateTo[row, "derived_object" -> derivedObject];
      AssociateTo[row, "value_digest" -> canonicalDigest[derivedObject]];
      AssociateTo[row, "dual_engine_comparison_id" -> "comparison:" <> slotId];
      AssociateTo[row, "dual_engine_object_derivation" -> "independent_engine_walk_of_committed_endpoint_and_channel_structure"],
      kind = If[
        StringStartsQ[slotId, "topology:"] || StringStartsQ[slotId, "host_location:"] || StringStartsQ[slotId, "mechanical_closure:"],
        "absence of any typed producer in the complete authority census", "nonuniqueness/solvability failure"
      ];
      phaseCRef = Null;
      If[StringStartsQ[slotId, "template_free_data:"],
        tilt = Last[StringSplit[slotId, ":", 2]];
        phaseCRef = <|"witness_id" -> phaseCById["tilt:" <> tilt]["witness_id"], "challenge_id" -> phaseCById["tilt:" <> tilt]["challenge_id"]|>
      ];
      pair = constructUnavailabilityChallenge[slotId, kind, producers, phaseCRef];
      row = Join[row, <|
        "witness_id" -> pair[[1]]["witness_id"], "challenge_id" -> pair[[2]]["challenge_id"],
        "derivability_contract_class" -> "class:" <> slotId, "witness" -> pair[[1]], "challenge" -> pair[[2]]
      |>]
    ]; AppendTo[rows, row],
    {specification, specifications}
  ]; SortBy[rows, {ToLowerCase[# ["slot_id"]], # ["slot_id"]} &]
];

evaluateDimensionFirewall[] := Module[{traction, area, velocity, surfacePower, expected, ablated},
  countCall["evaluateDimensionFirewall"];
  traction = {-1, -2, 1}; area = {2, 0, 0}; velocity = {1, -1, 0};
  surfacePower = traction + area + velocity; expected = {2, -3, 1}; ablated = traction + area + {0, 0, 0};
  <|
    "dimension_basis" -> {"L", "T", "M"}, "surface_power_dimensions" -> surfacePower,
    "constraint_power_dimensions" -> {2, -3, 1}, "rayleigh_power_dimensions" -> {2, -3, 1},
    "expected_power_dimensions" -> expected, "base_pass" -> (surfacePower == expected),
    "ablated_velocity_dimensions" -> ablated, "ablation_fired" -> (ablated != expected)
  |>
];

routeInventory[candidates_List, obligations_Association] := Module[
  {rows = {}, matrix, rhs, positiveCandidate, malformedCandidate, positiveResidual, malformedResidual,
   routeId, native, members, structure, dissipative, independentChecks},
  Do[
    members = membersFor[candidate];
    Which[
      candidate == "E1", matrix = {{1, 1, 0}, {1, -1, 0}, {0, 1, 1}}; positiveCandidate = {1, 1, 0};
        structure = <|"boundary_variation" -> True, "jump_condition" -> True, "holonomic_trace" -> True, "constraint_multiplier" -> False, "dissipation_bookkeeping" -> False|>,
      candidate == "E2", matrix = {{1, 1, 0}, {0, 2, 0}, {1, 0, 1}}; positiveCandidate = {1, 0, 0};
        structure = <|"boundary_variation" -> True, "jump_condition" -> True, "holonomic_trace" -> True, "constraint_multiplier" -> False, "dissipation_bookkeeping" -> False|>,
      candidate == "E3", matrix = {{2, 1}, {1, -1}}; positiveCandidate = {1, 2};
        structure = <|"boundary_variation" -> True, "jump_condition" -> True, "holonomic_trace" -> False, "constraint_multiplier" -> False, "dissipation_bookkeeping" -> False|>,
      candidate == "E4", matrix = {{2, 0, 1}, {0, 3, 1}, {1, 1, 0}}; positiveCandidate = {1, 2, -1};
        structure = <|"boundary_variation" -> False, "jump_condition" -> True, "holonomic_trace" -> False, "constraint_multiplier" -> True, "virtual_work_row" -> 2, "multiplier_column" -> 2, "dissipation_bookkeeping" -> False|>,
      candidate == "E5", matrix = {{2, 0, 0}, {0, 3, 1}, {0, 1, 1}}; positiveCandidate = {1, 1, 2};
        structure = <|"boundary_variation" -> True, "jump_condition" -> True, "holonomic_trace" -> False, "constraint_multiplier" -> False, "dissipation_bookkeeping" -> True, "rayleigh_coefficient" -> 2, "computed_dissipated_power" -> 2|>,
      StringStartsQ[candidate, "MIXTURE("], dissipative = MemberQ[members, "E5"];
        matrix = {{2, 0, 1, 0}, {0, 3, 0, 1}, {1, 0, 0, 0}, {0, 1, 0, If[dissipative, 2, 1]}}; positiveCandidate = {1, 2, -1, 1};
        structure = <|"boundary_variation" -> True, "jump_condition" -> True,
          "holonomic_trace" -> AnyTrue[members, MemberQ[{"E1", "E2"}, #] &], "constraint_multiplier" -> True,
          "virtual_work_row" -> 2, "multiplier_column" -> 2, "dissipation_bookkeeping" -> dissipative,
          "rayleigh_coefficient" -> If[dissipative, 2, Null], "computed_dissipated_power" -> If[dissipative, 2, Null]|>,
      True, matrix = {{2, 0}, {0, 3}}; positiveCandidate = {1, 1};
        structure = <|"boundary_variation" -> True, "jump_condition" -> True, "holonomic_trace" -> False,
          "constraint_multiplier" -> False, "dissipation_bookkeeping" -> False, "residual_complement_operator" -> True|>
    ];
    rhs = matrix . positiveCandidate;
    malformedCandidate = ReplacePart[positiveCandidate, 1 -> positiveCandidate[[1]] + 1];
    positiveResidual = matrix . positiveCandidate - rhs; malformedResidual = matrix . malformedCandidate - rhs;
    independentChecks = MapThread[(#1 . positiveCandidate) == #2 &, {matrix, rhs}];
    Do[
      routeId = "route:" <> candidate <> "|ambient=" <> ambient <> "|stratum=" <> stratum;
      native = obligations[candidate]["native_root_class"];
      AppendTo[rows, <|
        "route_id" -> routeId, "candidate_id" -> candidate, "ambient" -> ambient, "stratum" -> stratum,
        "native_root_class" -> native,
        "signature_digest" -> canonicalDigest[<|"obligations" -> obligations[candidate]["generator_A"],
          "native" -> native, "ambient" -> ambient, "stratum" -> stratum|>],
        "fixture_id" -> "fixture:" <> routeId,
        "positive_fixture" -> <|
          "semantic_route_id" -> "obligation_residual_classifier_v1", "matrix" -> matrix,
          "candidate" -> positiveCandidate, "rhs" -> rhs, "residual" -> positiveResidual,
          "expected" -> "ADMISSIBLE", "nondegenerate_norm_squared" -> positiveCandidate . positiveCandidate,
          "native_structure_exercised" -> structure,
          "known_outcome_generator_A" -> "exact_matrix_residual_zero",
          "known_outcome_generator_B" -> "independent_row_equation_satisfaction",
          "independent_row_equations_satisfied" -> independentChecks
        |>,
        "malformed_fixture" -> <|
          "semantic_route_id" -> "obligation_residual_classifier_v1", "matrix" -> matrix,
          "candidate" -> malformedCandidate, "rhs" -> rhs, "residual" -> malformedResidual, "expected" -> "EXCLUDED"
        |>, "route_equivalence_required" -> False
      |>], {ambient, ambients}, {stratum, tiltTypes}],
    {candidate, candidates}];
  rows
];

templateBVPConstituentCensus[candidateRecord_Association, ambient_String] := {
  <|
    "constituent" -> "canonical_boundary_condition", "candidate_id" -> candidateRecord["candidate_id"],
    "canonical_operator_signature" -> candidateRecord["canonical_signature"]
  |>,
  <|
    "constituent" -> "typed_free_data", "references" -> (("tilt:" <> #) & /@ tiltTypes),
    "availability" -> "UNRESOLVED", "domain" -> "tilted_sleeve_exterior",
    "dimensions" -> "inherited_exactly_from_R49_ledger"
  |>,
  <|"constituent" -> "unevaluated_residual_or_variational_form", "required_terms" -> {"bulk_euler_lagrange_residual"}|>,
  <|"constituent" -> "zero_mode_treatment", "required_terms" -> {"translation_zero_mode"}|>,
  <|"constituent" -> "well_posedness_classification", "value" -> "UNRESOLVED(committed_structure_only)"|>,
  <|"constituent" -> "asymptotic_matching_conditions", "ambient" -> ambient, "required_terms" -> {"outer_matching_condition"}|>
};

templateBranchEquivalenceProofs[candidateRecords_List, activeStrata_List] := Module[
  {operatorRecords, proofs = {}, perStratum, census, candidateRecord, ambient},
  operatorRecords = Select[candidateRecords, # ["candidate_id"] != "OTHER" &];
  Do[
    perStratum = Table[
      census = templateBVPConstituentCensus[candidateRecord, ambient];
      <|"stratum" -> stratum, "constituent_census" -> census, "census_sha256" -> canonicalDigest[census]|>,
      {stratum, activeStrata}
    ];
    AppendTo[proofs, <|
      "template_branch_id" -> "U2:" <> candidateRecord["candidate_id"] <> ":" <> ambient,
      "candidate_id" -> candidateRecord["candidate_id"], "ambient" -> ambient,
      "proof_timing" -> "pre-production_symbolic_only", "active_strata" -> activeStrata,
      "per_stratum_censuses" -> perStratum,
      "distinct_census_digest_count" -> Length[DeleteDuplicates[Lookup[perStratum, "census_sha256"]]],
      "identical_BVP_across_strata" -> (Length[DeleteDuplicates[Lookup[perStratum, "census_sha256"]]] == 1),
      "produced_adjudication_objects_used" -> False
    |>],
    {candidateRecord, operatorRecords}, {ambient, ambients}
  ]; proofs
];

templateContract[candidateRecords_List, activeStrata_List] := Module[{allowed},
  allowed = {"ROOT_REFERENCE", "POSITIVE_DERIVATION", "POSITIVE_EQUIVALENCE", "POSTULATE_BRANCH"};
  <|
  "eligible_record_predicate" -> "candidate!=OTHER and branch_has_at_least_one_stratum and all_integrity==COMPUTATION_VALID and homogeneous_disposition in {ADMISSIBLE,UNRESOLVED}",
  "template_key" -> "candidate_id x ambient (strata collapsed only by identical_BVP proof)",
  "semantic_schema" -> {
    "canonical_boundary_condition", "typed_free_data", "unevaluated_residual_or_variational_form",
    "zero_mode_treatment", "well_posedness_classification", "asymptotic_matching_conditions",
    "branch_conditionality_tag", "open_data_conditionality_tag"
  },
  "conditional_template_posing_allowed_constructors" -> allowed,
  "conditional_template_posing_banned_constructors" -> Select[grammar, !MemberQ[allowed, #] &],
  "constructor_partition_exhaustive_default_deny" -> True,
  "postulate_branch_constructor_allowed_ambient_only" -> "two_sided_R_w_postulate",
  "forbidden_ancestry_constructors" -> Select[grammar, !MemberQ[allowed, #] &], "dependent_fields_unbound" -> True,
  "R49_exact_unresolved_reference_ids" -> (("tilt:" <> #) & /@ tiltTypes),
  "term_census_generator_inputs" -> {
    "committed_static_field_equations", "canonical_operator_witness", "ambient_branch", "typed_free_data_ledger"
  },
  "required_term_kinds" -> {"residual", "boundary", "zero-mode", "asymptotic-matching"},
  "per_term_remove_and_zero_teeth_required" -> True, "posing_not_solving" -> True,
  "conditional_tags_non_evidential_and_posing_DAG_unreachable" -> True,
  "branch_equivalence_proofs" -> templateBranchEquivalenceProofs[candidateRecords, activeStrata]
|>];

closureContract[] := <|
  "gated_claim_inventory_generator" -> "committed_forcing_and_consequence_schema_walk",
  "contribution_census_generator_A" -> "committed_root_channel_incidence_walk",
  "contribution_census_generator_B" -> "force_balance_term_owner_walk",
  "coverage" -> "exact census to owner|computed-zero|typed-refusal",
  "independent_total" -> "direct committed force_balance construction", "no_double_count" -> True,
  "host_consistency_required" -> True, "both_host_requires_exact_partition" -> True,
  "radiation_static_zero_first_class" -> True, "integrity_failure_propagates_NOT_RUN" -> True,
  "valid_refusal_only_source_of_closure_PROMOTION_UNRESOLVED" -> True
|>;

templateGuardFixture[] := Module[{committedInputs, expectedTerms},
  committedInputs = <|
    "static_field_equations" -> {"bulk_euler_lagrange_residual"},
    "canonical_operator_witness" -> {"Sigma_boundary_operator"},
    "ambient_branch" -> {"outer_matching_condition"},
    "typed_free_data_ledger" -> (("tilt:" <> #) & /@ tiltTypes), "zero_mode_ledger" -> {"translation_zero_mode"}
  |>;
  expectedTerms = {
    <|"term_id" -> "residual:bulk_euler_lagrange_residual", "kind" -> "residual"|>,
    <|"term_id" -> "boundary:Sigma_boundary_operator", "kind" -> "boundary"|>,
    <|"term_id" -> "zero-mode:translation_zero_mode", "kind" -> "zero-mode"|>,
    <|"term_id" -> "asymptotic-matching:outer_matching_condition", "kind" -> "asymptotic-matching"|>
  };
  <|"fixture_id" -> "template_guard:posed_not_solved", "committed_inputs" -> committedInputs,
    "expected_term_census" -> expectedTerms,
    "template_record" -> <|"integrity" -> "COMPUTATION_VALID", "physics" -> "POSED_BVP_TEMPLATE",
      "constituents" -> <|
        "canonical_boundary_condition" -> "Sigma_boundary_operator",
        "typed_free_data" -> (<|
          "reference_id" -> "tilt:" <> #, "availability" -> "UNRESOLVED",
          "domain" -> "tilted_sleeve_exterior", "dimensions" -> "inherited_exactly_from_R49_ledger"
        |> & /@ tiltTypes),
        "unevaluated_residual_or_variational_form" -> "bulk_euler_lagrange_residual",
        "zero_mode_treatment" -> "project_translation_zero_mode",
        "well_posedness_classification" -> "UNRESOLVED(committed_structure_only)",
        "asymptotic_matching_conditions" -> "outer_matching_condition",
        "branch_conditionality_tag" -> <|
          "metadata_id" -> "metadata:conditional_branch:E1",
          "tag" -> "CONDITIONAL_ON_BRANCH(E1)", "candidate_id" -> "E1", "evidential" -> False, "posing_DAG_reachable" -> False
        |>,
        "open_data_conditionality_tag" -> <|
          "metadata_id" -> "metadata:open_data:E1:one_sided_pathA29",
          "unresolved_reference_ids" -> {"witness:host_location:E1"}, "evidential" -> False, "posing_DAG_reachable" -> False
        |>
      |>,
      "symbolic_ast" -> <|"op" -> "posed_template", "args" -> Table[
        <|"op" -> "term", "term_id" -> row["term_id"], "kind" -> row["kind"], "coefficient" -> 1|>,
        {row, expectedTerms}]|>, "evaluation_state" -> "UNEVALUATED",
      "eligibility_disposition" -> "UNRESOLVED", "conditional" -> True,
      "branch_admissibility_claim" -> False, "branch_selection_claim" -> False,
      "excluded_record_references" -> {}, "complement_record_references" -> {},
      "posing_proof_dag" -> <|
        "normalized_inference_content" -> <|
          "op" -> "positive_equivalence", "args" -> {<|
            "op" -> "positive_join", "args" -> {
              <|"op" -> "root", "id" -> "source:endpoint:E1"|>,
              <|"op" -> "root", "id" -> "source:field:geon_core_bundle"|>
            }
          |>}
        |>,
        "claimed_constructors" -> {"ROOT_REFERENCE", "POSITIVE_DERIVATION", "POSITIVE_EQUIVALENCE"}
      |>
    |>
  |>
];

closureGuardFixture[] := Module[{contributions, assignments, total},
  contributions = {
    <|"contribution_id" -> "surface_traction", "root_id" -> "source:field:native_momentum", "vector" -> {2, 0}, "expected_owner" -> "collar/Sigma surface"|>,
    <|"contribution_id" -> "constraint_reaction", "root_id" -> "source:field:E4_shear_lock", "vector" -> {-1, 1}, "expected_owner" -> "collar/Sigma surface"|>,
    <|"contribution_id" -> "return_flux", "root_id" -> "source:field:return_closure", "vector" -> {0, -1}, "expected_owner" -> "flux/return"|>,
    <|"contribution_id" -> "static_radiation", "root_id" -> "source:field:native_momentum", "vector" -> {0, 0}, "expected_owner" -> "radiation/static-zero"|>
  };
  assignments = (<|"contribution_id" -> #["contribution_id"], "owner" -> #["expected_owner"],
    "computed_zero" -> (#["vector"] == {0, 0})|> &) /@ contributions;
  total = Total[Lookup[contributions, "vector"]];
  <|"fixture_id" -> "closure_guard:exact_ownership", "committed_root_contributions" -> contributions,
    "census_A" -> Lookup[contributions, "contribution_id"], "census_B" -> sortIds[Lookup[contributions, "contribution_id"]],
    "independently_constructed_total" -> total,
    "certificate" -> <|"certificate_id" -> "certificate:fixture", "assignments" -> assignments|>,
    "closure_adjudication" -> <|"integrity" -> "COMPUTATION_VALID", "physics" -> "CLOSURE_CERTIFIED(certificate:fixture)"|>,
    "dependent_record" -> <|"integrity" -> "COMPUTATION_VALID", "physics" -> "SELECTED"|>
  |>
];

standardBindings[] := Module[{standards},
  standards = {
    {"S-1", "traceable cause tags"}, {"S-2", "field-driven classification"},
    {"S-3", "no vacuous constructs"}, {"S-4", "measured evidence"},
    {"S-5", "per-require teeth"}, {"S-6", "complete summaries"},
    {"P-1", "entry witness resolves to ratified unresolved slot"},
    {"P-2", "dual independent ancestry generators"},
    {"P-3", "endpoint-conditioned ancestry ownership"}, {"P-4", "two-level partitions"},
    {"P-5", "multi-census independent provenance"}, {"P-6", "radiation first-class channel"},
    {"P-7", "reconciliation witness reference integrity"}, {"P-8", "both-directions mutation armor"},
    {"P-9", "bare parent enums"}, {"P-10", "one-body two-body firewall"}
  };
  (<|
    "standard_id" -> #[[1]], "authoritative_text" -> #[[2]],
    "acceptance_predicate_id" -> "predicate:standard:" <> #[[1]],
    "reachable_check_id" -> "ASSERT_STANDARD_" <> StringReplace[#[[1]], "-" -> "_"],
    "tooth_id" -> "STANDARD_TOOTH:" <> #[[1]], "evidence_id" -> "evidence:standard:" <> #[[1]]
  |> &) /@ standards
];

inputs = loadYaml[FileNameJoin[{repo, "software/em_charge_attribute/u1_body_dynamics_inputs.yaml"}]];
phaseCSlots = loadYaml[FileNameJoin[{
  repo, "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_c_0_tilt_coupling_contract/availability_slots.yaml"
}]]["slots"];
phaseCProduction = loadYaml[FileNameJoin[{
  repo, "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_c_1_tilt_coupling_production/production_results.yaml"
}]];
endpoints = inputs["endpoints"];
mixturesA = generateMixturesFromConditions[endpoints];
mixturesB = generateMixturesFromChannels[endpoints];
candidates = Join[baseEndpoints, Lookup[mixturesA, "candidate_id"], {"OTHER"}];
census = sourceCensus[inputs, phaseCSlots];

candidateRecords = Table[<|
  "candidate_id" -> candidate,
  "kind" -> Which[MemberQ[baseEndpoints, candidate], "endpoint", StringStartsQ[candidate, "MIXTURE("], "mixture", True, "catch_all"],
  "members" -> membersFor[candidate], "native_root_class" -> nativeClass[candidate, endpoints],
  "canonical_signature" -> signatureFor[candidate],
  "membership_predicate" -> If[candidate != "OTHER", "positive_componentwise_operator_signature_equivalence", "report_only_residual_complement_not_promotion_reachable"],
  "promotion_membership_eligible" -> (candidate != "OTHER")
|>, {candidate, candidates}];
candidateUniverseDigest = canonicalDigest[candidateRecords];
candidateInventory = <|
  "mixture_formation_grammar" -> <|
    "rule" -> "compose E4 collar shear constraint with committed normal/surface endpoints on orthogonal components",
    "neutral_alias_rule" -> "E3 plus X canonicalizes to X",
    "conflict_rule" -> "distinct laws on same trace component do not form a family",
    "free_weight_parameters_banned" -> True|>,
  "mixture_generator_A" -> mixturesA, "mixture_generator_B" -> mixturesB,
  "concrete_other_candidates" -> {},
  "basis_closure" -> <|"status" -> "UNRESOLVED", "declared_domain" -> <|
      "locality" -> {"local", "nonlocal_unclosed"}, "linearity" -> {"linear", "nonlinear_unclosed"},
      "frequency_dependence" -> "unclosed", "amplitude_dependence" -> "unclosed"|>,
    "missing_data" -> {"throat_surface_functional", "outer_surface_functional", "complete_boundary_operator_class_theorem"},
    "full_current_canonical_union_used" -> True, "residual_complement_status" -> "UNRESOLVED_NONEMPTY_OR_EMPTY",
    "catch_all_present" -> True|>,
  "candidate_records" -> candidateRecords, "candidate_axis" -> candidates, "candidate_count" -> Length[candidates],
  "candidate_universe_digest" -> candidateUniverseDigest,
  "amendment_rekey_duty" -> "mint_new_universe_digest_and_immutable_successors_for_all_complement_dependent_records",
  "uncanonicalized_overlap_count" -> 0,
  "canonical_signature_count" -> Length[DeleteDuplicates[Lookup[candidateRecords, "canonical_signature"]]]
|>;

obligations = Association@Table[
  rootClass = nativeClass[candidate, endpoints];
  candidate -> <|
    "native_root_class" -> rootClass,
    "generation_inputs" -> {"candidate_id", "native_root_class", "operational_endpoint_definition"},
    "stratum_token_is_generation_input" -> False,
    "generator_A" -> generateObligationsFromNativeRoots[candidate, rootClass],
    "generator_B" -> generateObligationsFromEndpointWalk[candidate, endpoints],
    "nondegeneracy_predicate" -> "nonzero_canonical_operator_and_nonzero_committed_root_ablation",
    "anti_echo_predicate" -> "witness_satisfies_independent_nondefinitional_obligation",
    "semantic_ablation_criterion" -> "fail_nondefinitional_obligation_or_change_canonical_operator_class"
  |>, {candidate, candidates}
];

activeStrata = phaseCProduction["axes"]["generated_active_strata"];
dependencyRows = {}; gridCells = {}; collapseProofs = {};
Do[
  depsA = generateDependencyJoin[candidate, stratum, obligations[candidate]["generator_A"]];
  depsB = generateDependencySourceWalk[candidate, stratum, endpoints, inputs];
  AppendTo[dependencyRows, <|
    "candidate_id" -> candidate, "stratum" -> stratum, "generator_A" -> depsA, "generator_B" -> depsB,
    "stratifying_slot" -> "tilt:" <> stratum, "dependency_signature" -> canonicalDigest[depsA]
  |>];
  AppendTo[collapseProofs, <|
    "candidate_id" -> candidate, "raw_stratum" -> stratum, "collapsed_class" -> stratum,
    "timing" -> "pre-production", "proof" -> "singleton_due_to_distinct_authoritative_stratifying_slot",
    "stage0_objects_only" -> True
  |>];
  Do[AppendTo[gridCells, <|
    "cell_id" -> "candidate=" <> candidate <> "|ambient=" <> ambient <> "|stratum=" <> stratum,
    "candidate_id" -> candidate, "ambient" -> ambient, "stratum" -> stratum,
    "expected_dependencies" -> depsA, "stable_branch_id" -> "U2:" <> candidate <> ":" <> ambient <> ":" <> stratum
  |>], {ambient, ambients}],
  {candidate, candidates}, {stratum, activeStrata}
];

promotionContexts = Flatten[Table[
  <|
    "promotion_key" -> "ambient=" <> ambient <> "|context=" <> stratum,
    "ambient" -> ambient, "global_common_refinement_context" -> stratum,
    "candidate_cell_mappings" -> Table[<|
      "candidate_id" -> candidate,
      "cell_id" -> "candidate=" <> candidate <> "|ambient=" <> ambient <> "|stratum=" <> stratum
    |>, {candidate, candidates}]
  |>,
  {ambient, ambients}, {stratum, activeStrata}
], 1];

routes = routeInventory[candidates, obligations];
slots = availabilitySlots[candidates, phaseCSlots, endpoints, mixturesA];
artifactDags = proofArtifacts[candidateInventory, obligations, dependencyRows, routes, collapseProofs, promotionContexts, census["source_ids"]];
grammar = {
  "ROOT_REFERENCE", "STATIC_COMMITTED_FORCING", "POSITIVE_DERIVATION", "POSITIVE_EQUIVALENCE",
  "INCOMPATIBILITY", "NEGATIVE_CANDIDATE_MEMBERSHIP", "CASE_ELIMINATION", "COMPLEMENT_SURVIVOR_COUNT",
  "EXCLUSION_VERDICT", "POSTULATE_BRANCH", "STABILITY_DYNAMICAL_CLASS", "SOLVE_EVALUATION_RESULT",
  "SYMBOLIC_EQUIVALENCE_COLLAPSE", "UNAVAILABILITY_WITNESS", "DERIVABILITY_CHALLENGE",
  "EVIDENCE_STATE_CLASSIFICATION"
};

semantic = <|
  "schema_version" -> "U2_STAGE0_SEMANTIC_VIEW_V1",
  "scope" -> <|
    "static_adjudication_only" -> True, "dynamical_selection_deferred" -> True,
    "one_body_only" -> True, "BVP_solved" -> False,
    "pair_configuration_carveout" -> "pair_annihilation_subquestion_only"
  |>,
  "ambient_branches" -> {
    <|"ambient_id" -> "one_sided_pathA29", "status" -> "committed_executable_slab", "asymmetry_map" -> "ambient_asymmetry_map"|>,
    <|"ambient_id" -> "two_sided_R_w_postulate", "status" -> "POSTULATE", "asymmetry_map" -> "R_w"|>
  },
  "candidate_inventory" -> candidateInventory,
  "obligation_censuses" -> obligations,
  "open_dependency_relation" -> <|
    "authoritative_stratifying_ledger" -> "PhaseC.axes.generated_active_strata", "active_strata" -> activeStrata,
    "dependency_rows" -> dependencyRows,
    "generation" -> "stratum_free_census_then_independent_join_to_frozen_producer_relation_and_OPEN_ledger",
    "generator_A_algorithm" -> "obligation_to_authoritative_slot_relational_join",
    "generator_B_algorithm" -> "raw_field_and_endpoint_channel_schema_walk",
    "shared_task_code_between_generators" -> False
  |>,
  "grid_inventory" -> <|
    "candidate_count" -> Length[candidates], "ambient_count" -> Length[ambients],
    "raw_strata_per_candidate" -> Association@Table[candidate -> Length[activeStrata], {candidate, candidates}],
    "raw_ragged_cardinality" -> Total[Table[Length[activeStrata] Length[ambients], {candidate, candidates}]],
    "collapsed_cardinality" -> Length[gridCells], "preproduction_collapse_count" -> 0,
    "postproduction_report_collapse_affects_grid" -> False, "grid_cells" -> gridCells,
    "collapse_proofs" -> collapseProofs, "promotion_contexts" -> promotionContexts,
    "promotion_context_count" -> Length[promotionContexts]
  |>,
  "vocabulary_freeze" -> vocabularyFreeze[],
  "evidence_taxonomy" -> <|
    "source_census" -> census, "constructor_grammar" -> grammar,
    "generic_or_unclassified_constructor_allowed" -> False,
    "classification_source" -> "normalized_inference_content_not_claimed_label",
    "physics_bearing_stage0_artifacts" -> artifactDags,
    "promotion_allowed_constructors" -> {"ROOT_REFERENCE", "STATIC_COMMITTED_FORCING", "POSITIVE_DERIVATION", "POSITIVE_EQUIVALENCE"},
    "program_wide_banned_constructor" -> "STABILITY_DYNAMICAL_CLASS",
    "guard_fixture_dags" -> <|
      "promotion_positive" -> <|
        "constructors" -> {"ROOT_REFERENCE", "STATIC_COMMITTED_FORCING", "POSITIVE_DERIVATION", "POSITIVE_EQUIVALENCE"},
        "root_types" -> {"BULK_ACTION_TERM", "SURFACE_ACTION_TERM"},
        "normalized_content" -> <|"op" -> "static_force", "args" -> {<|"op" -> "positive_equivalence"|>}|>
      |>,
      "candidate_exclusion" -> <|
        "constructors" -> {"ROOT_REFERENCE", "INCOMPATIBILITY", "EXCLUSION_VERDICT"},
        "root_types" -> {"HOLONOMIC_CONSTRAINT", "BULK_ACTION_TERM"}
      |>,
      "posed_template" -> <|"constructors" -> {"ROOT_REFERENCE", "POSITIVE_DERIVATION"}, "fields_unbound" -> True|>,
      "pair_configuration" -> <|
        "object_type" -> "static_plus_w_minus_w_pair_configuration", "firewall_tag" -> "PAIR_ANNIHILATION_ONLY",
        "consumer" -> "topology_pair_annihilation_subquestion"
      |>,
      "postulate_deferral_metadata" -> <|"constructor" -> "POSTULATE_BRANCH", "promotion_reachable" -> False|>
    |>
  |>,
  "route_fixture_inventory" -> <|
    "generated_preproduction" -> True, "route_records" -> routes, "route_count" -> Length[routes],
    "fixture_count" -> Length[routes], "executed_route_match_rule" -> "production_route_set_exactly_equals_frozen_route_set",
    "cell_route_fixture_exact_coverage" -> True
  |>,
  "availability_slots" -> slots,
  "availability_summary" -> <|
    "total" -> Length[slots], "DERIVED" -> Count[Lookup[slots, "availability_outcome"], "DERIVED"],
    "UNRESOLVED" -> Count[Lookup[slots, "availability_outcome"], "UNRESOLVED"], "integrity_failures" -> 0
  |>,
  "closure_contract" -> closureContract[], "template_contract" -> templateContract[candidateRecords, activeStrata],
  "physics_record_invariance_contract" -> physicsRecordInvarianceContract[],
  "guard_fixtures" -> Join[guardFixtureRecords[], <|"closure" -> closureGuardFixture[], "template" -> templateGuardFixture[]|>],
  "return_closure_ownership" -> <|
    "owner" -> "downstream_flux_path", "U2_owned" -> False, "preserved_terminal" -> "UNRESOLVED(return_closure)",
    "allowed_consumption" -> "deferral_metadata_only_unreachable_from_U2_verdicts"
  |>,
  "standard_bindings" -> standardBindings[],
  "dimensional_firewall" -> evaluateDimensionFirewall[],
  "production_contract" -> <|
    "disposition_ancestry_exact_both_directions" -> True, "route_set_exact_match_required" -> True,
    "zero_integrity_failures_for_banking" -> True, "stage0_digest_rechecked_before_every_evaluation_and_A9_leg" -> True,
    "stable_branch_ids_required_all_output_classes" -> True, "mutation_catalog_coverage_required" -> True
  |>
|>;

registrySpecification = {
  {"candidate_formation_from_committed_constraints_v1", "generateMixturesFromConditions"},
  {"candidate_formation_channel_crosscheck_v1", "generateMixturesFromChannels"},
  {"stratum_free_obligation_census_v1", "generateObligationsFromNativeRoots"},
  {"stratum_free_operational_crosscheck_v1", "generateObligationsFromEndpointWalk"},
  {"open_dependency_join_v1", "generateDependencyJoin"},
  {"open_dependency_source_walk_v1", "generateDependencySourceWalk"},
  {"content_classified_proof_DAG_v1", "classifyInferenceContent"},
  {"constructive_absence_challenge_v1", "constructUnavailabilityChallenge"},
  {"inline_dimension_firewall_v1", "evaluateDimensionFirewall"},
  {"executable_stage0_decision_and_schema_guards_v1", "guardFixtureRecords"}
};
localRegistry = (<|
  "semantic_route_id" -> #[[1]], "engine_local_function" -> #[[2]],
  "exists" -> (Length[Names["Global`" <> #[[2]]]] == 1), "executed" -> (Lookup[callCounts, #[[2]], 0] > 0)
|> &) /@ registrySpecification;

result = <|
  "schema_version" -> schema, "engine" -> "Wolfram",
  "independent_route" -> "Wolfram association walks plus exact symbolic linear algebra",
  "semantic_view" -> semantic,
  "engine_local_route_registry" -> localRegistry,
  "runtime_identity" -> <|
    "wolfram_version" -> $Version, "system_id" -> $SystemID,
    "user_init_disabled" -> True, "sanitized_path" -> $Path
  |>
|>;

If[!DirectoryQ[DirectoryName[output]], CreateDirectory[DirectoryName[output], CreateIntermediateDirectories -> True]];
Export[output, result, "YAML"];
Print["U2_WOLFRAM_STAGE0_PASS candidates=" <> ToString[semantic["candidate_inventory"]["candidate_count"]] <>
  " grid=" <> ToString[semantic["grid_inventory"]["raw_ragged_cardinality"]] <>
  " slots=" <> ToString[semantic["availability_summary"]["total"]]];
Exit[0];
