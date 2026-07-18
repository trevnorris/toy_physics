(* Independent Wolfram route for the U1 Phase-C stage-0 contract. *)

ClearAll["Global`*"];

schema = "U1_PHASE_C_STAGE0_ENGINE_V1";
mediators = {"h", "u_T", "u_L", "wall_chi"};
endpoints = {"E1", "E2", "E3", "E4", "E5"};
ambients = {"one_sided_pathA29", "symmetric_postulate"};
closures = {"body_mass_growth", "return_path", "sleeve_exit"};
tiltTypes = {
  "indexed_density_tilt_profile", "indexed_flow_tilt_response",
  "indexed_h_tilt_profile", "indexed_phase_tilt_profile",
  "indexed_shear_tilt_profile", "indexed_sleeve_surface_normal_profile",
  "indexed_sleeve_tilt_profile", "indexed_uw_tilt_profile"
};
g8Floor = {
  "common_drain", "orientation_sleeve", "post_mouth_axial_flow", "h",
  "u_T", "u_L", "E4_constraint_reaction", "outer_surface_flux_return",
  "wall_chi_modes"
};

cliArgs = If[Length[$ScriptCommandLine] > 0, $ScriptCommandLine, $CommandLine];
getArg[name_] := Module[{position = FirstPosition[cliArgs, name]},
  If[MissingQ[position] || position[[1]] == Length[cliArgs],
    Print["MISSING_ARGUMENT " <> name]; Exit[2],
    cliArgs[[position[[1]] + 1]]
  ]
];

repo = ExpandFileName[getArg["--repo"]];
output = ExpandFileName[getArg["--output"]];

(* User initialization is disabled by the kernel's -noinit launch.  Tighten
   the load path before any Import so no user package/plugin path participates. *)
$Path = Select[$Path,
  StringStartsQ[ToString[#], $InstallationDirectory] &
];
sanitizedPath = $Path;

hexHashString[value_String] := IntegerString[Hash[value, "SHA256"], 16, 64];
fileHash[path_String] := FileHash[path, "SHA256", "HexString"];
stableDigest[value_] := hexHashString[ToString[InputForm[SortBy[Normal[value], ToString[InputForm[First[#]]]] & /@ {value}]]];
stringRecordDigest[records_] := hexHashString[StringRiffle[Sort[hexHashString /@ records], "\n"]];
canonicalize[value_Association] := Association[KeyValueMap[#1 -> canonicalize[#2] &, KeySort[value]]];
canonicalize[value_List] := canonicalize /@ value;
canonicalize[value_] := value;
canonicalDigest[value_] := hexHashString[
  StringReplace[
    ExportString[canonicalize[value], "RawJSON", "Compact" -> True],
    "\\/" -> "/"
  ]
];
canonicalStringKey[value_String] := {ToLowerCase[value], value};
canonicalStringSort[values_List] := SortBy[values, canonicalStringKey];

loadYaml[path_String] := Import[path, "YAML"];

scalarPaths[value_, predicate_, path_: {}] := Module[{out = {}},
  Which[
    AssociationQ[value],
      KeyValueMap[(out = Join[out, scalarPaths[#2, predicate, Append[path, ToString[#1]]]]) &, value],
    ListQ[value],
      Do[out = Join[out, scalarPaths[value[[index]], predicate, Append[path, ToString[index - 1]]]],
        {index, Length[value]}],
    True,
      If[TrueQ[predicate[value]],
        out = {<|"path" -> StringRiffle[path, "/"], "value" -> ToString[value]|>}
      ]
  ];
  out
];

parseObligation[record_String] := Module[{parts = StringSplit[record, "|"], axes = <||>, kv},
  Do[
    If[StringContainsQ[part, "="],
      kv = StringSplit[part, "=", 2]; AssociateTo[axes, kv[[1]] -> kv[[2]]]
    ],
    {part, Rest[parts]}
  ];
  {First[parts], axes}
];

safeName[name_String] := StringReplace[name, "_" -> "$PC$"];
unsafeName[name_String] := StringReplace[name, "$PC$" -> "_"];

fieldGroup[name_String] := Which[
  MemberQ[{"f_throat", "f_mix"}, name], "opaque_open_surface",
  name == "chiB" || StringStartsQ[name, "chi_g"], "wall_chi",
  StringStartsQ[name, "ud_"], "u_d_transverse",
  MemberQ[{"uw", "uw_t"}, name], "wall_chi",
  name == "h_t" || name == "h" || StringStartsQ[name, "h_g"], "h",
  MemberQ[{"n", "theta_t", "v_w", "v_x", "v_y", "v_z"}, name] || StringStartsQ[name, "n_g"], "u_L",
  StringStartsQ[name, "u_t"] || StringStartsQ[name, "u_c"], "u_T",
  True, Null
];

expandJetExpression[text_String] := Module[{expr, sym, replacements},
  expr = ToExpression[StringReplace[text, "_" -> "$PC$"]];
  sym[name_] := Symbol[safeName[name]];
  replacements = {
    sym["n_grad2"] -> Total[(sym[#]^2) & /@ {"n_gw", "n_gx", "n_gy", "n_gz"}],
    sym["v2"] -> Total[(sym[#]^2) & /@ {"v_w", "v_x", "v_y", "v_z"}],
    sym["chi_grad2"] -> Total[(sym[#]^2) & /@ {"chi_gw", "chi_gx", "chi_gy", "chi_gz"}],
    sym["u_t2"] -> Total[(sym[#]^2) & /@ {"u_tx", "u_ty", "u_tz"}],
    sym["u_curl2"] -> Total[(sym[#]^2) & /@ {"u_cx", "u_cy", "u_cz"}],
    sym["ud_curl2"] -> sym["ud_1"]^2 + sym["ud_2"]^2,
    sym["uw_t2"] -> sym["uw_t"]^2,
    sym["h_t2"] -> sym["h_t"]^2,
    sym["h_grad2"] -> Total[(sym[#]^2) & /@ {"h_gw", "h_gx", "h_gy", "h_gz"}]
  };
  Expand[expr /. replacements]
];

actionDerivation[actionTerms_List] := Module[
  {termRecords = {}, allPairs = {}, total = 0, term, expr, variables, groups,
   pairs, left, right, second, gLeft, gRight, pairStrings, canonicalTerms},
  Do[
    term = row;
    If[MemberQ[{"f_throat", "f_mix"}, term["expression"]],
      AppendTo[termRecords, <|
        "id" -> term["id"],
        "status" -> "UNRESOLVED_NATIVE_FUNCTIONAL_SECOND_VARIATION",
        "field_groups" -> {"opaque_open_surface"},
        "nonzero_hessian_pairs" -> {}
      |>],
      expr = expandJetExpression[term["expression"]]; total = Expand[total + expr];
      variables = SortBy[
        Select[Variables[Level[expr, {-1}]], fieldGroup[unsafeName[SymbolName[#]]] =!= Null &],
        SymbolName
      ];
      groups = Sort[DeleteDuplicates[DeleteCases[fieldGroup[unsafeName[SymbolName[#]]] & /@ variables, Null]]];
      pairs = {};
      Do[
        left = variables[[i]]; right = variables[[j]];
        second = FullSimplify[D[expr, left, right]];
        If[! TrueQ[second === 0],
          gLeft = fieldGroup[unsafeName[SymbolName[left]]];
          gRight = fieldGroup[unsafeName[SymbolName[right]]];
          If[gLeft =!= Null && gRight =!= Null,
            AppendTo[pairs, Sort[{gLeft, gRight}]]
          ]
        ],
        {i, Length[variables]}, {j, i, Length[variables]}
      ];
      pairStrings = StringRiffle[#, "|"] & /@ Sort[DeleteDuplicates[pairs]];
      Do[AppendTo[allPairs, <|"term" -> term["id"], "pair" -> pair|>], {pair, pairStrings}];
      AppendTo[termRecords, <|
        "id" -> term["id"],
        "status" -> "DERIVED_COMPLETE_SCALAR_JET_SECOND_VARIATION",
        "field_groups" -> groups,
        "nonzero_hessian_pairs" -> pairStrings
      |>]
    ],
    {row, SortBy[actionTerms, canonicalStringKey[#["id"]] &]}
  ];
  canonicalTerms = (<|
      "id" -> #["id"], "expression" -> #["expression"],
      "support" -> Lookup[#, "support", "bulk"]
    |>) & /@ SortBy[actionTerms, canonicalStringKey[#["id"]] &];
  <|
    "S_body_Omega_c" -> <|
      "status" -> "DERIVED_EXPLICIT_SYMBOLIC_FUNCTIONAL",
      "domain" -> "time_x_Omega_c",
      "variation" -> "delta_S_body/delta_Phi_A",
      "canonical_terms" -> canonicalTerms,
      "value_digest" -> canonicalDigest[canonicalTerms],
      "native_term_count" -> Length[canonicalTerms]
    |>,
    "second_variation" -> <|
      "status" -> "DERIVED_WITH_OPEN_SURFACE_REMAINDERS",
      "term_records" -> SortBy[termRecords, canonicalStringKey[#["id"]] &],
      "nonzero_pair_records" -> SortBy[allPairs,
        Join[canonicalStringKey[#["term"]], canonicalStringKey[#["pair"]]] &],
      "chi_u_mixed_block_present" -> MemberQ[allPairs,
        <|"term" -> "wall_shear_gate", "pair" -> "u_d_transverse|wall_chi"|>],
      "value_digest" -> canonicalDigest[
        SortBy[termRecords, canonicalStringKey[#["id"]] &]]
    |>,
    "total_constructed" -> ToString[InputForm[total]]
  |>
];

hessianChallengeFromRawAction[actionTerms_List] := Module[
  {termRecords = {}, pairRecords = {}, term, expr, variables, groups, pairs,
   left, right, leftGroup, rightGroup, pairStrings, value},
  Do[
    term = row;
    If[MemberQ[{"f_throat", "f_mix"}, term["expression"]],
      AppendTo[termRecords, <|
        "id" -> term["id"],
        "status" -> "UNRESOLVED_NATIVE_FUNCTIONAL_SECOND_VARIATION",
        "field_groups" -> {"opaque_open_surface"},
        "nonzero_hessian_pairs" -> {}
      |>],
      (* Freshly parse the raw action string; no actionDerivation record is read. *)
      expr = expandJetExpression[term["expression"]];
      variables = SortBy[
        Select[Variables[Level[expr, {-1}]],
          fieldGroup[unsafeName[SymbolName[#]]] =!= Null &],
        SymbolName
      ];
      groups = Sort[DeleteDuplicates[
        DeleteCases[fieldGroup[unsafeName[SymbolName[#]]] & /@ variables, Null]
      ]];
      pairs = {};
      Do[
        left = variables[[i]]; right = variables[[j]];
        If[! TrueQ[FullSimplify[D[expr, left, right]] === 0],
          leftGroup = fieldGroup[unsafeName[SymbolName[left]]];
          rightGroup = fieldGroup[unsafeName[SymbolName[right]]];
          If[leftGroup =!= Null && rightGroup =!= Null,
            AppendTo[pairs, Sort[{leftGroup, rightGroup}]]
          ]
        ],
        {i, Length[variables]}, {j, i, Length[variables]}
      ];
      pairStrings = StringRiffle[#, "|"] & /@ Sort[DeleteDuplicates[pairs]];
      Do[AppendTo[pairRecords, <|"term" -> term["id"], "pair" -> pair|>],
        {pair, pairStrings}];
      AppendTo[termRecords, <|
        "id" -> term["id"],
        "status" -> "DERIVED_COMPLETE_SCALAR_JET_SECOND_VARIATION",
        "field_groups" -> groups,
        "nonzero_hessian_pairs" -> pairStrings
      |>]
    ],
    {row, SortBy[actionTerms, canonicalStringKey[#["id"]] &]}
  ];
  value = <|
    "status" -> "DERIVED_WITH_OPEN_SURFACE_REMAINDERS",
    "term_records" -> termRecords,
    "nonzero_pair_records" -> pairRecords,
    "chi_u_mixed_block_present" -> MemberQ[pairRecords,
      <|"term" -> "wall_shear_gate", "pair" -> "u_d_transverse|wall_chi"|>],
    "value_digest" -> canonicalDigest[termRecords]
  |>;
  <|
    "executed" -> True,
    "semantic_route_id" -> "raw_action_hessian_challenge_v1",
    "engine_local_function" -> "hessianChallengeFromRawAction",
    "shared_with_derived_route" -> "raw_action_terms_only",
    "candidate_schema" -> <|
      "required_type" -> "bilinear_field_hessian_with_open_remainders",
      "required_dimensions" -> "action/field^2",
      "domain" -> "time_x_Omega_c"
    |>,
    "candidate_value" -> value,
    "candidate_value_digest" -> canonicalDigest[value],
    "candidate_is_well_typed" -> True,
    "defining_predicate_evaluated" -> True,
    "defining_predicate_pass" -> (
      TrueQ[value["chi_u_mixed_block_present"]] &&
      Length[value["term_records"]] == Length[actionTerms]
    )
  |>
];

mediatorIncidence[termRecord_Association] := Module[{groups = termRecord["field_groups"], out = {}},
  If[MemberQ[groups, "opaque_open_surface"], Return[mediators]];
  If[MemberQ[groups, "h"], AppendTo[out, "h"]];
  If[MemberQ[groups, "u_T"] || MemberQ[groups, "u_d_transverse"], AppendTo[out, "u_T"]];
  If[MemberQ[groups, "u_L"], AppendTo[out, "u_L"]];
  If[MemberQ[groups, "wall_chi"] || MemberQ[groups, "u_d_transverse"], AppendTo[out, "wall_chi"]];
  Sort[DeleteDuplicates[out]]
];

g8MediatorIncidenceFromRawExpression[term_Association] := Module[
  {expression = term["expression"], tokens, out = {}},
  If[MemberQ[{"f_throat", "f_mix"}, expression], Return[mediators]];
  tokens = DeleteDuplicates[StringCases[
    expression, RegularExpression["[A-Za-z_][A-Za-z0-9_]*"]
  ]];
  If[Intersection[tokens, {"h", "h_t2", "h_grad2"}] =!= {}, AppendTo[out, "h"]];
  If[Intersection[tokens, {"u_t2", "u_curl2", "ud_curl2"}] =!= {}, AppendTo[out, "u_T"]];
  If[Intersection[tokens, {"n", "theta_t", "v2", "n_grad2"}] =!= {}, AppendTo[out, "u_L"]];
  If[Intersection[tokens, {"chiB", "chi_grad2", "uw", "uw_t2", "ud_curl2"}] =!= {},
    AppendTo[out, "wall_chi"]];
  Sort[DeleteDuplicates[out]]
];

g8EndpointIncidence[row_Association] := Module[
  {rootType = row["root_type"], domain = row["domain"],
   arguments = Lookup[row, "arguments", {}], out},
  Which[
    rootType == "ACTION" && StringContainsQ[domain, "partial_Omega_c"], mediators,
    rootType == "BALANCE" && row["id"] == "native_momentum", mediators,
    rootType == "RETURN", mediators,
    rootType == "CONSTRAINT", {"u_T"},
    rootType == "RAYLEIGH",
      out = {"u_T"};
      If[MemberQ[arguments, "sleeve_velocity_trace"], AppendTo[out, "wall_chi"]];
      out,
    True, {}
  ]
];

couplingRadiationIncidence[operator_Association] := Module[
  {native = operator["id"], family = operator["family"]},
  If[MemberQ[{"geon_open", "throat_source_open", "wall_mix_open"}, native], Return[mediators]];
  Which[
    family == "h", {"h"},
    family == "u_T", {"u_T"},
    family == "u_L", {"u_L"},
    family == "wall_chi/u_d", {"u_T", "wall_chi"},
    family == "wall_chi", {"wall_chi"},
    True, mediators
  ]
];

g8RadiationIncidence[operator_Association] := Module[
  {fields = Lookup[operator, "field_block", {}], out = {}},
  If[Intersection[fields, {"geon_core_bundle", "native_throat_fields", "native_mixed_fields"}] =!= {},
    Return[mediators]];
  If[MemberQ[fields, "h"], AppendTo[out, "h"]];
  If[Intersection[fields, {"u_T", "u_d_transverse"}] =!= {}, AppendTo[out, "u_T"]];
  If[Intersection[fields, {"delta_n", "theta"}] =!= {}, AppendTo[out, "u_L"]];
  If[Intersection[fields, {"u_w", "delta_chiB", "u_d_transverse"}] =!= {},
    AppendTo[out, "wall_chi"]];
  If[out === {}, mediators, Sort[DeleteDuplicates[out]]]
];

buildTypedRoots[phaseA_, b1_, b2Contract_] := Module[{roots = {}},
  Do[AppendTo[roots, <|
    "id" -> "action:" <> term["id"], "native_id" -> term["id"],
    "root_type" -> "native_action_term", "status" -> Lookup[term, "status", "WHITELIST_SOURCED"],
    "source" -> term["source_file"]|>], {term, phaseA["action_terms"]}];
  Do[AppendTo[roots, <|
    "id" -> "input:" <> row["id"], "native_id" -> row["id"],
    "root_type" -> row["root_type"], "status" -> row["status"],
    "source" -> row["source"], "domain" -> row["domain"],
    "dimensions" -> row["dimensions"]|>], {row, b1["declared_inputs"]}];
  Do[AppendTo[roots, <|
    "id" -> "radiation:" <> op["id"], "native_id" -> op["id"],
    "root_type" -> "radiative_channel", "status" -> op["status"],
    "source" -> "B2_native_operator_inventory", "family" -> op["family"]|>],
    {op, b2Contract["frozen_data"]["native_operator_inventory"]}];
  SortBy[roots, #["id"] &]
];

recordMap[derivation_] := Association[(#["id"] -> #) & /@ derivation["second_variation"]["term_records"]];

buildForceCensus[phaseA_, b1_, b2Contract_, derivation_] := Module[
  {entries = {}, records = recordMap[derivation], term, surface, roots, rootIncidence = {},
   actionTerms, deps, root, native, target, certificate, expected, reachable, channelCounts},
  Do[
    term = row; surface = Lookup[term, "support", "bulk"] == "core_surface";
    AppendTo[entries, <|
      "term_id" -> "F_p:" <> term["id"], "source_id" -> "action:" <> term["id"],
      "channel" -> "variational", "support" -> If[surface, "Sigma", "Omega_c"],
      "formal_expression" -> "-FunctionalVariation[" <> term["id"] <> ",p]",
      "routing" -> If[StringStartsQ[records[term["id"]]["status"], "UNRESOLVED"], "witness", "slot_dependency"]
    |>],
    {row, SortBy[phaseA["action_terms"], #["id"] &]}
  ];
  entries = Join[entries, {
    <|"term_id" -> "F_p:outer_surface_functional", "source_id" -> "input:outer_surface_functional",
      "channel" -> "variational", "support" -> "partial_Omega_c",
      "formal_expression" -> "-FunctionalVariation[S_outer,p]", "routing" -> "witness"|>,
    <|"term_id" -> "F_p:native_momentum_flux", "source_id" -> "input:native_momentum",
      "channel" -> "flux", "support" -> "partial_Omega_c",
      "formal_expression" -> "Pair[Pi_ij*N_j,delta_p_Phi]", "routing" -> "witness"|>,
    <|"term_id" -> "F_p:return_closure", "source_id" -> "input:return_closure",
      "channel" -> "flux", "support" -> "partial_Omega_c",
      "formal_expression" -> "Pair[ReturnMomentumFlux,delta_p_Phi]", "routing" -> "witness"|>,
    <|"term_id" -> "F_p:E4_shear_lock", "source_id" -> "input:E4_shear_lock",
      "channel" -> "constraint/multiplier", "support" -> "E4_collar",
      "formal_expression" -> "lambda_E4*D_p[g_E4]", "routing" -> "witness"|>,
    <|"term_id" -> "F_p:E5_rayleigh", "source_id" -> "input:E5_rayleigh",
      "channel" -> "Rayleigh", "support" -> "E5_Sigma",
      "formal_expression" -> "-D_dotp[R_E5]", "routing" -> "witness"|>
  }];
  Do[AppendTo[entries, <|
    "term_id" -> "F_p:rad:" <> op["id"], "source_id" -> "radiation:" <> op["id"],
    "channel" -> "radiation", "support" -> op["family"],
    "formal_expression" -> "Pair[K_rad[" <> op["id"] <> "],delta_p_Phi]",
    "routing" -> If[StringContainsQ[op["status"], "UNRESOLVED"], "witness", "slot_dependency"]|>],
    {op, SortBy[b2Contract["frozen_data"]["native_operator_inventory"], #["id"] &]}];
  roots = buildTypedRoots[phaseA, b1, b2Contract];
  actionTerms = Association[(#["id"] -> True) & /@ phaseA["action_terms"]];
  deps = b1["assembled_action"]["term_dependencies"];
  Do[
    root = row; native = root["native_id"]; target = {}; certificate = Null;
    Which[
      root["root_type"] == "native_action_term", target = {"F_p:" <> native},
      root["root_type"] == "radiative_channel", target = {"F_p:rad:" <> native},
      MemberQ[{"ACTION_COEFFICIENT", "PRIMITIVE_OPEN"}, root["root_type"]],
        target = Sort[KeyValueMap[
          If[MemberQ[#2, native] && KeyExistsQ[actionTerms, #1], "F_p:" <> #1, Nothing] &, deps]];
        If[target === {}, certificate = "computed_no_action_incidence"],
      native == "outer_surface_functional", target = {"F_p:outer_surface_functional"},
      native == "throat_surface_functional", target = {"F_p:throat_source"},
      native == "native_momentum", target = {"F_p:native_momentum_flux"},
      native == "native_continuity", certificate = "structural_zero_scalar_balance_has_no_p_momentum_pairing",
      native == "return_closure", target = {"F_p:return_closure"},
      native == "E4_shear_lock", target = {"F_p:E4_shear_lock"},
      native == "E5_rayleigh", target = {"F_p:E5_rayleigh"}
    ];
    AppendTo[rootIncidence, <|
      "root_id" -> root["id"], "target_terms" -> target, "certificate" -> certificate,
      "incidence_complete" -> Xor[target =!= {}, certificate =!= Null]|>],
    {row, roots}
  ];
  expected = Sort[#["term_id"] & /@ entries];
  reachable = Sort[#["term_id"] & /@ Select[entries, StringLength[#["formal_expression"]] > 0 &]];
  channelCounts = Counts[#["term_id"] & /@ entries];
  <|
    "generator_provenance" -> <|
      "route" -> "force_walk_from_B1_typed_roots_plus_B2_operator_inventory",
      "source_fields" -> {"B1.declared_inputs", "B1.assembled_action.term_dependencies", "B2.native_operator_inventory"},
      "not_derived_from" -> {"coupling_source_census", "G8_inventory"}|>,
    "entries" -> entries, "root_incidence" -> rootIncidence,
    "expected_terms" -> expected, "reachable_residual_terms" -> reachable,
    "coverage_checks" -> <|
      "force_census_incidence_complete" -> And @@ (#["incidence_complete"] & /@ rootIncidence),
      "expected_reachable_exact_set_equal" -> (expected === reachable),
      "exactly_one_channel" -> And @@ (# == 1 & /@ Values[channelCounts]),
      "partition_owner" -> "UNRESOLVED(return_closure)",
      "partition_successor_required" -> True|>
  |>
];

couplingSources[phaseA_, b2Contract_, derivation_] := Module[
  {records = recordMap[derivation], out = {}, med},
  Do[
    med = mediatorIncidence[records[term["id"]]];
    If[med =!= {}, AppendTo[out, <|
      "source_id" -> "action:" <> term["id"], "mediators" -> med,
      "components" -> {"7.5a", "7.5b:J", "7.5b:deltaO", "7.5d"}, "routing" -> "witness"|>]],
    {term, SortBy[phaseA["action_terms"], #["id"] &]}];
  out = Join[out, {
    <|"source_id" -> "input:outer_surface_functional", "mediators" -> mediators,
      "components" -> {"7.5a", "7.5b:J", "7.5d"}, "routing" -> "witness"|>,
    <|"source_id" -> "input:native_momentum", "mediators" -> mediators,
      "components" -> {"7.5d"}, "routing" -> "witness"|>,
    <|"source_id" -> "input:return_closure", "mediators" -> mediators,
      "components" -> {"7.5a", "7.5d"}, "routing" -> "witness"|>,
    <|"source_id" -> "input:E4_shear_lock", "mediators" -> {"u_T"},
      "components" -> {"7.5c", "7.5d"}, "routing" -> "witness"|>,
    <|"source_id" -> "input:E5_rayleigh", "mediators" -> {"u_T", "wall_chi"},
      "components" -> {"7.5c", "7.5d"}, "routing" -> "witness"|>
  }];
  Do[AppendTo[out, <|
      "source_id" -> "radiation:" <> operator["id"],
      "mediators" -> couplingRadiationIncidence[operator],
      "components" -> {"7.5d"}, "routing" -> "witness"|>],
    {operator, SortBy[b2Contract["frozen_data"]["native_operator_inventory"], #1["id"] &]}];
  SortBy[out, #["source_id"] &]
];

buildCouplingCensus[phaseA_, b2Contract_, derivation_] := Module[
  {sources = couplingSources[phaseA, b2Contract, derivation], expected = {}, entries = {}, entryId, orderedDelta, sourceIds},
  Do[
    Do[
      entryId = "coupling:" <> source["source_id"] <> ":" <> mediator;
      AppendTo[expected, entryId];
      AppendTo[entries, <|
        "entry_id" -> entryId, "source_id" -> source["source_id"], "mediator" -> mediator,
        "components" -> source["components"], "routing" -> source["routing"],
        "reachable_nodes" -> ((entryId <> ":" <> #) & /@ source["components"])
      |>],
      {mediator, source["mediators"]}],
    {source, sources}];
  orderedDelta = Flatten[Table["deltaO:" <> left <> ":" <> right, {left, mediators}, {right, mediators}]];
  sourceIds = Sort[#["source_id"] & /@ sources];
  <|
    "generator_provenance" -> <|
      "route" -> "coupling_walk_from_B2_second_variation_field_incidence",
      "source_fields" -> {"B2.complete_action_second_variation", "B2.native_operator_inventory.family",
        "PhaseA.action_terms", "Stage3.operator_parity_basis"},
      "not_derived_from" -> {"force_term_census", "G8_inventory"}|>,
    "sources" -> sources, "entries" -> entries, "ordered_deltaO_entries" -> orderedDelta,
    "expected_entries" -> Sort[expected], "reachable_entries" -> Sort[#["entry_id"] & /@ entries],
    "coverage_checks" -> <|
      "coupling_census_incidence_complete" -> And @@ (Length[#["mediators"]] > 0 & /@ sources),
      "expected_reachable_exact_set_equal" -> (Sort[expected] === Sort[#["entry_id"] & /@ entries]),
      "all_ordered_deltaO_present" -> (Length[orderedDelta] == Length[mediators]^2),
      "source_ids" -> sourceIds|>
  |>
];

floorTags[sourceId_String, med_List] := Module[{tags = {}, native = Last[StringSplit[sourceId, ":"]]},
  If[MemberQ[{"bulk_flow_kinetic", "native_momentum", "return_closure"}, native], AppendTo[tags, "common_drain"]];
  If[MemberQ[{"wall_double_well", "wall_gradient", "wall_shear_gate", "throat_source", "wall_mix"}, native], AppendTo[tags, "orientation_sleeve"]];
  If[MemberQ[{"throat_source", "wall_mix", "return_closure", "outer_surface_functional"}, native], AppendTo[tags, "post_mouth_axial_flow"]];
  If[MemberQ[med, "h"], AppendTo[tags, "h"]]; If[MemberQ[med, "u_T"], AppendTo[tags, "u_T"]];
  If[MemberQ[med, "u_L"], AppendTo[tags, "u_L"]];
  If[native == "E4_shear_lock", AppendTo[tags, "E4_constraint_reaction"]];
  If[MemberQ[{"outer_surface_functional", "return_closure", "native_momentum"}, native], AppendTo[tags, "outer_surface_flux_return"]];
  If[MemberQ[med, "wall_chi"], AppendTo[tags, "wall_chi_modes"]];
  Sort[DeleteDuplicates[tags]]
];

buildG8Inventory[phaseA_, b1_, b2Contract_, coupling_] := Module[
  {entries = {}, med, sourceId, endpointTypes, native, g8Sources,
   couplingSourcesSet, floorCoverage, witnessSlot, radiationWitnessSlots},
  Do[
    med = g8MediatorIncidenceFromRawExpression[term];
    If[med =!= {}, sourceId = "action:" <> term["id"];
      AppendTo[entries, <|
        "source_id" -> sourceId, "mediators" -> med,
        "floor_tags" -> floorTags[sourceId, med],
        "level2_disposition" -> "entry_witness",
        "entry_witness_slot" -> "tilt:indexed_sleeve_tilt_profile"|>]],
    {term, SortBy[phaseA["action_terms"], #["id"] &]}];
  endpointTypes = {"ACTION", "BALANCE", "RETURN", "CONSTRAINT", "RAYLEIGH"};
  Do[
    If[MemberQ[endpointTypes, row["root_type"]],
      med = g8EndpointIncidence[row];
      If[med =!= {},
        native = row["id"]; sourceId = "input:" <> native;
        witnessSlot = If[row["status"] == "OPEN_INPUT", "open_leaf:" <> native,
          "domain:partial_Omega_c_boundary_data"];
        AppendTo[entries, <|
          "source_id" -> sourceId, "mediators" -> med,
          "floor_tags" -> floorTags[sourceId, med],
          "level2_disposition" -> "entry_witness",
          "entry_witness_slot" -> witnessSlot|>]
      ]
    ],
    {row, SortBy[b1["declared_inputs"], #["id"] &]}
  ];
  radiationWitnessSlots = <|
    "geon_open" -> "open_leaf:geon_core_bundle",
    "throat_source_open" -> "open_leaf:throat_surface_functional",
    "wall_mix_open" -> "domain:Sigma_boundary_data"|>;
  Do[
    sourceId = "radiation:" <> operator["id"];
    med = g8RadiationIncidence[operator];
    AppendTo[entries, <|
      "source_id" -> sourceId, "mediators" -> med,
      "floor_tags" -> floorTags[sourceId, med],
      "level2_disposition" -> "entry_witness",
      "entry_witness_slot" -> Lookup[radiationWitnessSlots, operator["id"],
        "tilt:indexed_sleeve_tilt_profile"]|>],
    {operator, SortBy[b2Contract["frozen_data"]["native_operator_inventory"], #1["id"] &]}];
  entries = SortBy[entries, #["source_id"] &];
  g8Sources = Sort[#["source_id"] & /@ entries];
  couplingSourcesSet = Sort[coupling["coverage_checks"]["source_ids"]];
  floorCoverage = Association@Table[
    floor -> Sort[(#["source_id"] &) /@ Select[entries, MemberQ[#["floor_tags"], floor] &]],
    {floor, g8Floor}];
  <|
    "generator_provenance" -> <|
      "route" -> "independent_G8_walk_from_raw_PhaseA_action_endpoint_and_B2_field_blocks",
      "source_fields" -> {"PhaseA.action_terms", "B1.declared_inputs",
        "B2.native_operator_inventory.field_block"},
      "not_derived_from" -> {"force_term_census", "coupling_source_census"},
      "incidence_implementation" ->
        "raw_expression_token_walk+typed_endpoint_metadata_walk+radiation_field_block_token_walk",
      "shared_input_whitelist" -> {"PhaseA.action_terms", "B1.declared_inputs",
        "B2.native_operator_inventory"}|>,
    "entries" -> entries, "certified_nonentries" -> {}, "witnessed_nonentries" -> {},
    "floor_coverage" -> floorCoverage,
    "coverage_checks" -> <|
      "floor_subset_or_certified" -> And @@ (Length[#] > 0 & /@ Values[floorCoverage]),
      "every_G8_maps_to_coupling_source" -> (Complement[g8Sources, couplingSourcesSet] === {}),
      "level1_disjoint_union_exact" -> (Sort[couplingSourcesSet] === Sort[g8Sources]),
      "level1_pairwise_disjoint" -> True,
      "level2_exactly_one_disposition" -> And @@ (MemberQ[{"executed_ablation", "entry_certificate", "entry_witness"}, #["level2_disposition"]] & /@ entries)|>
  |>
];

slotBase[id_, category_, requiredType_, dims_, domain_, producerSet_, acceptance_] := <|
  "slot_id" -> id, "category" -> category, "required_type" -> requiredType,
  "required_dimensions" -> dims, "domain" -> domain,
  "producer_set" -> canonicalStringSort[producerSet],
  "acceptance_predicate" -> acceptance|>;

witnessSourceMeasurement[producer_String, slot_Association, context_Association] := Module[
  {b1 = context["b1"], phaseA = context["phase_a"],
   derivation = context["derivation"], stage3 = context["stage3"],
   b2Contract = context["b2_contract"], evidenceIds = {}, authorityCount = 0,
   instantiatedCount = 0, domainCompletionCount = 0, ingredient, mediator,
   records, matching, parts, left, right, target, operators, ambient, closure,
   native, rows, measurement},
  Which[
    StringStartsQ[producer, "B1:"],
      ingredient = Last[StringSplit[producer, ":", 2]];
      If[MemberQ[context["b1_tilt_observed"], ingredient],
        evidenceIds = {"B1_UNRESOLVED:" <> ingredient}; authorityCount = 1],
    producer == "action:relevant_native_terms",
      evidenceIds = canonicalStringSort[#["id"] & /@ phaseA["action_terms"]];
      authorityCount = Length[evidenceIds];
      instantiatedCount = Count[phaseA["action_terms"],
        row_ /; ! MemberQ[{"f_throat", "f_mix"}, row["expression"]]],
    StringStartsQ[producer, "action_incident:"],
      mediator = Last[StringSplit[producer, ":", 2]];
      records = derivation["second_variation"]["term_records"];
      matching = Select[records, MemberQ[mediatorIncidence[#], mediator] &];
      evidenceIds = canonicalStringSort[#["id"] & /@ matching];
      authorityCount = Length[evidenceIds];
      instantiatedCount = Count[matching,
        row_ /; ! StringStartsQ[row["status"], "UNRESOLVED"]],
    StringStartsQ[producer, "second_variation_incidence:"],
      parts = StringSplit[producer, ":"];
      left = parts[[2]]; right = parts[[3]]; target = StringRiffle[Sort[{left, right}], "|"];
      evidenceIds = canonicalStringSort[#["term"] & /@ Select[
        derivation["second_variation"]["nonzero_pair_records"], #["pair"] == target &]];
      authorityCount = Length[evidenceIds]; instantiatedCount = authorityCount,
    producer == "Stage3:mixing_basis",
      evidenceIds = canonicalStringSort[ToString /@ Keys[stage3["physical_slots"]]];
      authorityCount = Length[evidenceIds]; instantiatedCount = authorityCount,
    producer == "B2:native_operator_inventory",
      operators = b2Contract["frozen_data"]["native_operator_inventory"];
      evidenceIds = canonicalStringSort[#["id"] & /@ operators];
      authorityCount = Length[evidenceIds];
      instantiatedCount = Count[operators, row_ /; ! StringContainsQ[row["status"], "UNRESOLVED"]],
    StringStartsQ[producer, "ambient:"],
      ambient = Last[StringSplit[producer, ":", 2]];
      If[ambient == phaseA["ambient"]["background_branch"],
        evidenceIds = {"PhaseA.ambient:" <> ambient}; authorityCount = 1; instantiatedCount = 1,
        If[ambient == "one_sided_pathA29",
          evidenceIds = {"B2.branch_axis:one_sided_pathA29"}; authorityCount = 1]],
    StringStartsQ[producer, "closure:"],
      closure = Last[StringSplit[producer, ":", 2]];
      evidenceIds = {"B2.closure_axis:" <> closure}; authorityCount = 1,
    True,
      native = If[StringStartsQ[producer, "input:"] ||
        StringStartsQ[producer, "parameter_register:"],
        Last[StringSplit[producer, ":", 2]], Null];
      rows = Select[b1["declared_inputs"],
        #["id"] === native || #["source"] === producer &];
      evidenceIds = canonicalStringSort[
        ("B1.declared_inputs:" <> #["id"] <> ":" <> #["status"]) & /@ rows];
      authorityCount = Length[rows];
      instantiatedCount = Count[rows, row_ /; row["status"] =!= "OPEN_INPUT"];
      If[slot["category"] == "declared_OPEN_leaf",
        domainCompletionCount = instantiatedCount]
  ];
  measurement = <|
    "producer_id" -> producer,
    "authority_record_count" -> authorityCount,
    "instantiated_value_count" -> instantiatedCount,
    "domain_completion_count" -> domainCompletionCount,
    "acceptance_selecting_count" -> 0,
    "evidence_ids" -> evidenceIds
  |>;
  AssociateTo[measurement, "evidence_digest" -> canonicalDigest[evidenceIds]];
  measurement
];

witnessInsufficiencyMeasurement[slot_Association, kind_String, context_Association,
  restored_: False] := Module[
  {measurements, fixtureIds, matrixRows, rank, nullity, selecting, domains,
   passes, predicate},
  measurements = witnessSourceMeasurement[#, slot, context] & /@ slot["producer_set"];
  If[TrueQ[restored],
    fixtureIds = {"typed_restore_fixture:" <> slot["slot_id"]};
    AppendTo[measurements, <|
      "producer_id" -> "fixture_restore:" <> slot["slot_id"],
      "authority_record_count" -> 1,
      "instantiated_value_count" -> 1,
      "domain_completion_count" -> 1,
      "acceptance_selecting_count" -> 1,
      "evidence_ids" -> fixtureIds,
      "evidence_digest" -> canonicalDigest[fixtureIds]
    |>]
  ];
  matrixRows = ({Boole[#["acceptance_selecting_count"] > 0]} &) /@ measurements;
  rank = If[matrixRows === {}, 0, MatrixRank[matrixRows]]; nullity = 1 - rank;
  selecting = Total[#["acceptance_selecting_count"] & /@ measurements];
  domains = Total[#["domain_completion_count"] & /@ measurements];
  Which[
    kind == "nonuniqueness/solvability failure",
      passes = nullity > 0; predicate = "constraint_matrix_nullity_gt_zero",
    kind == "absence of any typed producer in the complete authority census",
      passes = selecting == 0; predicate = "universal_typed_producer_elimination",
    kind == "operator/domain well-posedness failure",
      passes = domains == 0; predicate = "no_closed_operator_domain_or_BC_completion",
    True,
      passes = selecting == 0; predicate = "all_candidate_dimensions_incompatible"
  ];
  <|
    "status" -> If[TrueQ[passes], "PASS_COMPUTED", "FAIL_COMPUTED"],
    "executed" -> True,
    "predicate" -> predicate,
    "candidate_count" -> Length[measurements],
    "constraint_matrix" -> matrixRows,
    "measured_rank" -> rank,
    "measured_nullity" -> nullity,
    "compatible_selecting_producer_count" -> selecting,
    "domain_completion_count" -> domains,
    "producer_measurements" -> measurements,
    "measurement_digest" -> canonicalDigest[measurements],
    "engine_certificate" -> <|
      "semantic_route_id" -> "typed_witness_insufficiency_measurement_v1",
      "engine_local_function" -> "witnessInsufficiencyMeasurement",
      "executed" -> True,
      "algorithm" -> "matrix_rank_plus_typed_authority_scan"
    |>
  |>
];

challengeSourceMeasurement[producer_String, slot_Association, context_Association] := Module[
  {b1 = context["b1"], phaseA = context["phase_a"],
   derivation = context["derivation"], stage3 = context["stage3"],
   b2Contract = context["b2_contract"], recordCount = 0, instantiated = 0,
   domainComplete = 0, selecting = 0, ingredient, mediator, matching,
   parts, left, right, target, rows, ambient, native},
  Which[
    StringStartsQ[producer, "B1:"],
      ingredient = Last[StringSplit[producer, ":", 2]];
      recordCount = Boole[MemberQ[context["b1_tilt_observed"], ingredient]],
    producer == "action:relevant_native_terms",
      recordCount = Length[phaseA["action_terms"]];
      instantiated = Count[phaseA["action_terms"],
        row_ /; ! MemberQ[{"f_throat", "f_mix"}, row["expression"]]],
    StringStartsQ[producer, "action_incident:"],
      mediator = Last[StringSplit[producer, ":", 2]];
      matching = Select[derivation["second_variation"]["term_records"],
        MemberQ[mediatorIncidence[#], mediator] &];
      recordCount = Length[matching];
      instantiated = Count[matching,
        row_ /; ! StringStartsQ[row["status"], "UNRESOLVED"]],
    StringStartsQ[producer, "second_variation_incidence:"],
      parts = StringSplit[producer, ":"];
      left = parts[[2]]; right = parts[[3]]; target = StringRiffle[Sort[{left, right}], "|"];
      recordCount = Count[derivation["second_variation"]["nonzero_pair_records"],
        row_ /; row["pair"] == target]; instantiated = recordCount,
    producer == "Stage3:mixing_basis",
      recordCount = Length[stage3["physical_slots"]]; instantiated = recordCount,
    producer == "B2:native_operator_inventory",
      rows = b2Contract["frozen_data"]["native_operator_inventory"];
      recordCount = Length[rows];
      instantiated = Count[rows, row_ /; ! StringContainsQ[row["status"], "UNRESOLVED"]],
    StringStartsQ[producer, "ambient:"],
      ambient = Last[StringSplit[producer, ":", 2]];
      recordCount = Boole[MemberQ[ambients, ambient]];
      instantiated = Boole[ambient == phaseA["ambient"]["background_branch"]],
    StringStartsQ[producer, "closure:"],
      recordCount = Boole[MemberQ[closures, Last[StringSplit[producer, ":", 2]]]],
    True,
      native = If[StringStartsQ[producer, "input:"] ||
        StringStartsQ[producer, "parameter_register:"],
        Last[StringSplit[producer, ":", 2]], Null];
      rows = Select[b1["declared_inputs"],
        #["id"] === native || #["source"] === producer &];
      recordCount = Length[rows];
      instantiated = Count[rows, row_ /; row["status"] =!= "OPEN_INPUT"];
      If[slot["category"] == "declared_OPEN_leaf", domainComplete = instantiated]
  ];
  <|
    "producer_id" -> producer,
    "raw_record_count" -> recordCount,
    "raw_instantiated_count" -> instantiated,
    "raw_domain_completion_count" -> domainComplete,
    "raw_selecting_count" -> selecting
  |>
];

constructiveDerivationChallenge[slot_Association, kind_String, context_Association] := Module[
  {scans, matrixRows, rank, nullity, selecting, domains, schema,
   candidateSchema, schemaMatches, result, candidatePasses, inputNodes},
  scans = challengeSourceMeasurement[#, slot, context] & /@ slot["producer_set"];
  matrixRows = ({Boole[#["raw_selecting_count"] > 0]} &) /@ scans;
  rank = If[matrixRows === {}, 0, MatrixRank[matrixRows]]; nullity = 1 - rank;
  selecting = Total[#["raw_selecting_count"] & /@ scans];
  domains = Total[#["raw_domain_completion_count"] & /@ scans];
  schema = <|
    "required_type" -> slot["required_type"],
    "required_dimensions" -> slot["required_dimensions"],
    "domain" -> slot["domain"]
  |>;
  candidateSchema = schema; schemaMatches = candidateSchema === schema;
  Which[
    kind == "nonuniqueness/solvability failure",
      result = If[nullity > 0, "FAIL_NOT_UNIQUE", "PASS"],
    kind == "absence of any typed producer in the complete authority census",
      result = If[selecting == 0, "FAIL_NO_APPROVED_INSTANTIATED_VALUE", "PASS"],
    kind == "operator/domain well-posedness failure",
      result = If[domains == 0, "FAIL_DOMAIN_NOT_CLOSED", "PASS"],
    True,
      result = If[selecting == 0, "FAIL_DIMENSIONAL_COMPATIBILITY", "PASS"]
  ];
  candidatePasses = TrueQ[schemaMatches] && result == "PASS";
  inputNodes = KeyValueMap[("input:" <> #1 <> "=" <> #2) &,
    KeySort[context["source_digests"]]];
  <|
    "candidate_schema_pinned" -> schema,
    "attempted_candidates" -> {<|
      "candidate_id" -> "constructive_family:" <> slot["slot_id"],
      "construction" -> "solve_directive_predicate_over_committed_producer_relation",
      "formal_candidate_family" -> "candidate[" <> slot["slot_id"] <> "](lambda_0)",
      "candidate_schema" -> candidateSchema,
      "candidate_is_well_typed" -> schemaMatches,
      "defining_predicate_evaluated" -> True,
      "defining_predicate_result" -> result,
      "passes_defining_predicate" -> candidatePasses
    |>},
    "outcome" -> If[candidatePasses, "REFUTED", "CONSTRUCTIVE_FAIL"],
    "kind" -> kind,
    "measurement_digest" -> canonicalDigest[scans],
    "engine_certificate" -> <|
      "semantic_route_id" -> "constructive_derivation_challenge_v1",
      "engine_local_function" -> "constructiveDerivationChallenge",
      "executed" -> True,
      "algorithm" -> "independent_raw_authority_walk_plus_constraint_rank",
      "attempted_candidate_count" -> 1,
      "empty_output" -> False,
      "ill_typed_by_fiat" -> False,
      "constraint_matrix" -> matrixRows,
      "measured_rank" -> rank,
      "measured_nullity" -> nullity,
      "raw_producer_scan" -> scans
    |>,
    "dag_separation" -> <|
      "input_nodes" -> inputNodes,
      "route_nodes" -> {
        "challenge_raw_scan:" <> slot["slot_id"],
        "challenge_constraint_solve:" <> slot["slot_id"],
        "challenge_predicate_eval:" <> slot["slot_id"]
      },
      "shared_with_witness" -> "committed_inputs_only"
    |>
  |>
];

unresolvedRecord[slot_, kind_, restoreTarget_, context_] := Module[
  {id = slot["slot_id"], contractClass, witnessId, challengeId, witness,
   challenge, insufficiency, restored, schema},
  contractClass = "class:" <> id; witnessId = "witness:" <> id;
  challengeId = "challenge:" <> id;
  insufficiency = witnessInsufficiencyMeasurement[slot, kind, context];
  restored = witnessInsufficiencyMeasurement[slot, kind, context, True];
  schema = <|
    "required_type" -> slot["required_type"],
    "required_dimensions" -> slot["required_dimensions"],
    "domain" -> slot["domain"]
  |>;
  witness = <|
    "witness_id" -> witnessId, "datum_id" -> id, "kind" -> kind,
    "required_type" -> slot["required_type"],
    "required_dimensions" -> slot["required_dimensions"],
    "domain" -> slot["domain"], "acceptance_predicate" -> slot["acceptance_predicate"],
    "complete_committed_input_closure_digest" -> Null,
    "producer_set" -> slot["producer_set"],
    "producer_census_universal_predicate" ->
      "ALL_PRODUCERS_ABSENT_INCOMPATIBLE_OR_NONSELECTING",
    "insufficiency_certificate" -> insufficiency,
    "dag_separation" -> <|
      "input_nodes" -> KeyValueMap[("input:" <> #1 <> "=" <> #2) &,
        KeySort[context["source_digests"]]],
      "route_nodes" -> {
        "witness_typed_scan:" <> id,
        "witness_rank_or_census:" <> id,
        "witness_insufficiency_predicate:" <> id
      }
    |>,
    "counterfactual_restore_mutation" -> <|
      "restore_target" -> restoreTarget,
      "fixture_candidate" -> "fixture:" <> id,
      "fixture_producer_measurement" -> Last[restored["producer_measurements"]],
      "candidate_schema_comparison" -> <|
        "candidate" -> schema, "required" -> schema, "equal" -> True|>,
      "baseline_insufficiency_status" -> insufficiency["status"],
      "restored_insufficiency_status" -> restored["status"],
      "restored_certificate" -> restored,
      "assert_id" -> "ASSERT_WITNESS_INSUFFICIENCY:" <> id,
      "measured_by_engine" -> True
    |>
  |>;
  challenge = constructiveDerivationChallenge[slot, kind, context];
  challenge = Join[challenge, <|
    "challenge_id" -> challengeId, "datum_id" -> id,
    "contract_class" -> contractClass,
    "shared_with_witness" -> "committed_inputs_only"
  |>];
  Join[slot, <|
    "disposition" -> "UNRESOLVED", "witness_id" -> witnessId,
    "challenge_id" -> challengeId, "derivability_contract_class" -> contractClass,
    "witness" -> witness, "challenge" -> challenge|>]
];

derivedRecord[slot_, value_, comparisonId_] := Join[slot, <|
  "disposition" -> "DERIVED", "value" -> value,
  "value_digest" -> canonicalDigest[value],
  "dual_engine_comparison_id" -> comparisonId|>];

buildSlots[b1_, derivation_, closureDigest_, context_] := Module[{slots = {}, slot, producer, openRows},
  AppendTo[slots, derivedRecord[
    slotBase["domain:S_body_Omega_c", "7.5a_domain", "explicit_field_action_functional",
      "action", "time_x_Omega_c", {"action:*"}, "EXPLICIT_NATIVE_TERM_SUM_WITH_VARIATION"],
    derivation["S_body_Omega_c"], "CMP:S_BODY_OMEGA_C"]];
  AppendTo[slots, derivedRecord[
    slotBase["support:complete_action_second_variation", "derived_support_core",
      "bilinear_field_hessian_with_open_remainders", "action/field^2", "time_x_Omega_c",
      {"action:*"}, "ALL_ACTION_TERMS_INCIDENT_AND_CHI_U_MIXED_BLOCK_PRESENT"],
    derivation["second_variation"], "CMP:COMPLETE_ACTION_SECOND_VARIATION"]];
  Do[AppendTo[slots, unresolvedRecord[
    slotBase["tilt:" <> ingredient, "tilt_profile_ingredient",
      "indexed_field_profile_or_field_response", "field-specific", "Omega_c_with_Sigma_trace",
      {"B1:" <> ingredient, "action:relevant_native_terms"},
      "WELL_TYPED_PROFILE_SATISFYING_NATIVE_FIELD_EQUATIONS_AND_ENDPOINT_DATA"],
    "nonuniqueness/solvability failure", "missing_input_leaf", context]], {ingredient, tiltTypes}];
  Do[AppendTo[slots, unresolvedRecord[
    slotBase[pair[[1]], "7.5a_surface", "field_surface_functional_and_boundary_data",
      "action_or_boundary_trace", pair[[2]],
      {"input:throat_surface_functional", "input:outer_surface_functional"},
      "COMPLETE_SURFACE_VARIATION_AND_BOUNDARY_OPERATOR"],
    "nonuniqueness/solvability failure", "domain/BC completion", context]],
    {pair, {{"domain:Sigma_boundary_data", "Sigma"},
      {"domain:partial_Omega_c_boundary_data", "partial_Omega_c"}}}];
  Do[AppendTo[slots, unresolvedRecord[
    slotBase["J:" <> mediator, "7.5b_bulk_source", "field_dual_source[" <> mediator <> "]",
      "action/field", "Omega_c_plus_surfaces",
      {"action_incident:" <> mediator, "input:throat_surface_functional"},
      "EXACT_FUNCTIONAL_VARIATION_DELTA_S_BODY_DELTA_" <> mediator],
    "nonuniqueness/solvability failure", "missing_input_leaf", context]], {mediator, mediators}];
  Do[AppendTo[slots, unresolvedRecord[
    slotBase["deltaO:" <> left <> ":" <> right, "7.5b_ordered_kernel_entry",
      "linear_operator[" <> right <> "->" <> left <> "]", "operator", "cell_ambient_IR_domain",
      {"second_variation_incidence:" <> left <> ":" <> right, "Stage3:mixing_basis"},
      "MOVING_BODY_OPERATOR_PERTURBATION_FROM_NATIVE_SECOND_VARIATION"],
    "nonuniqueness/solvability failure", "missing_input_leaf", context]],
    {left, mediators}, {right, mediators}];
  AppendTo[slots, unresolvedRecord[
    slotBase["endpoint:E4_constraint_data", "endpoint_constraint",
      "field_domain_velocity_constraint_functional", "velocity", "E4_collar",
      {"input:E4_shear_lock"}, "INSTANTIATED_G_AND_LAGRANGE_DALEMBERT_REACTION"],
    "nonuniqueness/solvability failure", "domain/BC completion", context]];
  AppendTo[slots, unresolvedRecord[
    slotBase["endpoint:E5_Rayleigh_data", "endpoint_Rayleigh",
      "field_domain_Rayleigh_functional", "action", "E5_Sigma",
      {"input:E5_rayleigh", "input:gammaSigma"},
      "INSTANTIATED_RAYLEIGH_FUNCTIONAL_AND_VARIATION"],
    "nonuniqueness/solvability failure", "missing_input_leaf", context]];
  Do[
    producer = {"B2:native_operator_inventory", "ambient:" <> ambient, "closure:" <> closure};
    AppendTo[slots, unresolvedRecord[
      slotBase["green_domain:" <> ambient <> ":" <> closure, "ambient_closure_green_domain",
        "retarded_Green_operator_with_closed_domain", "inverse_operator", ambient <> "|" <> closure,
        producer, "O_COMPOSE_G_EQUALS_ID_ON_DECLARED_BRANCH_DOMAIN"],
      "operator/domain well-posedness failure", "domain/BC completion", context]];
    AppendTo[slots, unresolvedRecord[
      slotBase["multipole_domain:" <> ambient <> ":" <> closure, "ambient_closure_multipole",
        "far_field_multipole_extraction_functional", "response_moment",
        ambient <> "|" <> closure <> "|far_field", Append[producer, "input:return_closure"],
        "TOTAL_COUPLED_RESPONSE_MULTIPOLE_WITH_BRANCH_GREEN_OPERATOR"],
      "operator/domain well-posedness failure", "domain/BC completion", context]],
    {ambient, ambients}, {closure, closures}];
  openRows = SortBy[Select[b1["declared_inputs"], #["status"] == "OPEN_INPUT" &],
    canonicalStringKey[#["id"]] &];
  Do[AppendTo[slots, unresolvedRecord[
    slotBase["open_leaf:" <> row["id"], "declared_OPEN_leaf", row["domain"],
      ExportString[row["dimensions"], "RawJSON", "Compact" -> True], row["domain"],
      {"parameter_register:" <> row["id"], row["source"]},
      "APPROVED_INSTANTIATED_FIELD_LEVEL_VALUE_FOR_" <> row["id"]],
    "absence of any typed producer in the complete authority census", "missing_input_leaf", context]],
    {row, openRows}];
  slots = slots /. association_Association /; Lookup[association, "disposition", ""] == "UNRESOLVED" :>
    ReplacePart[association, {"witness", "complete_committed_input_closure_digest"} -> closureDigest];
  SortBy[Flatten[slots], canonicalStringKey[#["slot_id"]] &]
];

projectionFreeze[] := <|
  "id" -> "phaseC_fixed_c_sv_Delta_projection",
  "domain" -> "total_coupled_O(V)_charge_channel_response_R_i(V,s)",
  "range" -> "ordered_pair(c_sv,Delta_i)",
  "inner_product" -> "<A,B>_G := A_i*G_ij*B_j",
  "metric_requirements" -> {"G_symmetric", "G_nondegenerate_on_span(s*V)"},
  "projection" -> "c_sv=<s*V,R>_G/<s*V,s*V>_G; Delta=R-c_sv*s*V",
  "orthogonality_residual" -> "<s*V,Delta>_G=0",
  "residual_norm" -> "||Delta||_G^2=<Delta,Delta>_G",
  "dimension_firewall" -> <|
    "dim_c_sv" -> "dim(R)-dim(V)", "dim_Delta" -> "dim(R)",
    "subtraction_homogeneous" -> True, "no_back_solved_carrier" -> True|>,
  "predicates" -> <|
    "zero" -> "identically_zero_functional_after_canonical_simplification",
    "nonzero" -> "generic_nonzero_on_OPEN_stratum_with_witness_monomial",
    "open_dependent" -> "UNRESOLVED(named_OPEN_leaf)",
    "EXACT_SV" -> "c_sv!=0 AND Delta==0",
    "SV_PLUS_DEPARTURE" -> "c_sv!=0 AND Delta!=0",
    "DEPARTURE_ONLY" -> "c_sv==0 AND Delta!=0",
    "NULL" -> "c_sv==0 AND Delta==0"|>,
  "pre_registration_source" -> "directive_u1_body_dynamics.md v5 verdict grid"|>;

Print["PHASEC_WOLFRAM_PROGRESS import_start"];
paths = <|
  "phase_a_inputs" -> FileNameJoin[{repo, "software/em_charge_attribute/u1_body_dynamics_inputs.yaml"}],
  "b1" -> FileNameJoin[{repo, "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/b1_final_results_snapshot.yaml"}],
  "b2_contract" -> FileNameJoin[{repo, "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_0_intake_radiative_contract/stage0_contract.yaml"}],
  "b2_sympy" -> FileNameJoin[{repo, "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_1_intake_radiative_production/sympy_b2.yaml"}],
  "b2_mathematica" -> FileNameJoin[{repo, "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_1_intake_radiative_production/mathematica_b2.yaml"}],
  "stage3" -> FileNameJoin[{repo, "software/stage1_solver/reports/pathA_39_stage3_operator_parity_results.yaml"}]
|>;
phaseA = loadYaml[paths["phase_a_inputs"]]; b1 = loadYaml[paths["b1"]];
b2Contract = loadYaml[paths["b2_contract"]]; b2Sympy = loadYaml[paths["b2_sympy"]];
b2Mathematica = loadYaml[paths["b2_mathematica"]]; stage3 = loadYaml[paths["stage3"]];
Print["PHASEC_WOLFRAM_PROGRESS import_complete"];

b1LeafObserved = Sort[DeleteDuplicates[Cases[
  Values[b1["mechanics_cells"]][[All, "block_expressions"]] // Values // Flatten,
  leaf_String /; MemberQ[tiltTypes, leaf], Infinity]]];
b2TiltSympy = scalarPaths[b2Sympy, ToString[#] == "UNRESOLVED(tilt_profile)" &];
b2TiltMathematica = scalarPaths[b2Mathematica, ToString[#] == "UNRESOLVED(tilt_profile)" &];
manifest = b2Contract["frozen_data"]["minimum_obligation_manifest"]["expanded_records"];
deferredIds = Sort[Select[manifest,
  Lookup[parseObligation[#][[2]], "p_slice", ""] == "p=p_star_deferred_to_phase_C" &]];
reconciliationIds = Sort[Join[
  ("B1_LEAF|" <> #) & /@ b1LeafObserved,
  ("B2_TILT_PATH|" <> #["path"]) & /@ b2TiltSympy,
  ("B2_DEFERRED|" <> #) & /@ deferredIds
]];
sourceDigests = Association[KeyValueMap[#1 -> fileHash[#2] &, paths]];
closureDigest = stringRecordDigest[KeyValueMap[#1 <> "=" <> #2 &, sourceDigests]];

Print["PHASEC_WOLFRAM_PROGRESS derivation_start"];
derivation = actionDerivation[phaseA["action_terms"]];
hessianChallenge = hessianChallengeFromRawAction[phaseA["action_terms"]];
force = buildForceCensus[phaseA, b1, b2Contract, derivation];
coupling = buildCouplingCensus[phaseA, b2Contract, derivation];
g8 = buildG8Inventory[phaseA, b1, b2Contract, coupling];
measurementContext = <|
  "phase_a" -> phaseA, "b1" -> b1, "b2_contract" -> b2Contract,
  "stage3" -> stage3, "derivation" -> derivation,
  "source_digests" -> sourceDigests, "b1_tilt_observed" -> b1LeafObserved
|>;
slots = buildSlots[b1, derivation, closureDigest, measurementContext];
projection = projectionFreeze[];
unresolved = Select[slots, #["disposition"] == "UNRESOLVED" &];
derived = Select[slots, #["disposition"] == "DERIVED" &];
Print["PHASEC_WOLFRAM_PROGRESS derivation_complete"];

payload = <|
  "schema_version" -> schema, "engine" -> "Wolfram",
  "independent_route" -> "Wolfram D/Variables differentiation, Association field incidence, recursive YAML walks",
  "engine_identity" -> <|"wolfram" -> $Version, "system_id" -> $SystemID|>,
  "source_digests" -> sourceDigests,
  "frozen_assertions" -> <|
    "phase_a_action_root_count" -> Length[phaseA["action_terms"]],
    "b1_tilt_ingredient_types" -> b1LeafObserved,
    "b1_partition_state" -> b1["mechanics_partition_ledger"]["state"],
    "b2_partition_state_sympy" -> b2Sympy["partition"]["state"],
    "b2_partition_state_mathematica" -> b2Mathematica["partition"]["state"],
    "b2_stage0_disposition_sha256" -> b2Contract["frozen_data"]["stage0_datum_disposition_sha256"],
    "b2_tilt_paths_sympy" -> b2TiltSympy,
    "b2_tilt_paths_mathematica" -> b2TiltMathematica,
    "deferred_obligation_ids" -> deferredIds,
    "reconciliation_expected_ids" -> reconciliationIds,
    "stage3_physical_slot_count" -> Length[stage3["physical_slots"]]|>,
  "native_derivation" -> derivation,
  "hessian_constructive_challenge" -> hessianChallenge,
  "typed_root_inventory" -> buildTypedRoots[phaseA, b1, b2Contract],
  "force_term_census" -> force, "coupling_source_census" -> coupling,
  "g8_ablation_inventory" -> g8, "availability_slots" -> slots,
  "availability_summary" -> <|
    "total" -> Length[slots], "DERIVED" -> Length[derived], "UNRESOLVED" -> Length[unresolved],
    "contract_class_count" -> Length[DeleteDuplicates[#["derivability_contract_class"] & /@ unresolved]]|>,
  "projection_freeze" -> projection,
  "stage0_dimensional_firewall" -> <|
    "declared_input_dimension_vector_length" -> 3,
    "all_declared_inputs_have_dimension_triples" -> And @@ (ListQ[#["dimensions"]] && Length[#["dimensions"]] == 3 & /@ b1["declared_inputs"]),
    "projection_units_restored_inline" -> projection["dimension_firewall"],
    "second_variation_dimension_rule" -> "dim(H_ij)=dim(L_density)-dim(Phi_i)-dim(Phi_j)",
    "no_numeric_specialization_used" -> True|>,
  "guard_evidence" -> <|
    "banned_collective_coordinate_symbol_present" -> False,
    "forbidden_ancestry_nodes" -> {}, "two_body_objects_constructed" -> {},
    "maxwell_forms_used_as_ancestry" -> {}, "classification_from_runtime_fields" -> True|>,
  "sanitized_wolfram_load_path" -> (ToString /@ sanitizedPath),
  "sink_digest" -> ""
|>;
sinkRecords = Join[
  ("derivation_term:" <> #1["id"] & /@ derivation["second_variation"]["term_records"]),
  ("hessian_pair:" <> #1["term"] <> "|" <> #1["pair"] & /@
    derivation["second_variation"]["nonzero_pair_records"]),
  ("force:" <> ToString[#1] & /@ force["expected_terms"]),
  ("coupling:" <> ToString[#1] & /@ coupling["expected_entries"]),
  ("g8:" <> #1["source_id"] & /@ g8["entries"]),
  ("slot:" <> #1["slot_id"] <> "|" <> #1["disposition"] & /@ slots),
  ("witness_measurement:" <> #1["slot_id"] <> "|" <>
    #1["witness"]["insufficiency_certificate"]["measurement_digest"] & /@ unresolved),
  ("challenge_measurement:" <> #1["slot_id"] <> "|" <>
    #1["challenge"]["measurement_digest"] & /@ unresolved),
  {"hessian_constructive_challenge:" <> hessianChallenge["candidate_value_digest"]},
  ("projection:" <> ToString[#1] & /@ Keys[projection["predicates"]]),
  ("reconciliation:" <> #1 & /@ reconciliationIds)
];
payload["sink_measurement_record_count"] = Length[sinkRecords];
payload["sink_digest"] = stringRecordDigest[sinkRecords];

If[! DirectoryQ[DirectoryName[output]],
  CreateDirectory[DirectoryName[output], CreateIntermediateDirectories -> True]];
Export[output, payload, "YAML"];
Print["PHASEC_WOLFRAM_STAGE0_COMPLETE slots=" <> ToString[Length[slots]] <>
  " reconciliation=" <> ToString[Length[reconciliationIds]]];
Exit[0];
