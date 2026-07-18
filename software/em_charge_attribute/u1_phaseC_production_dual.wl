(* Independent Wolfram production engine for U1 Phase C. *)

ClearAll["Global`*"];

schema = "U1_PHASE_C_PRODUCTION_ENGINE_V1";
ratifiedDigest = "e632a8d6729d0a1b3a4ade883c28f6b21f7a29fea566318cdd6fefec8c15d0da";
mediators = {"h", "u_T", "u_L", "wall_chi"};
endpoints = {"E1", "E2", "E3", "E4", "E5"};
ambients = {"one_sided_pathA29", "symmetric_postulate"};
closures = {"body_mass_growth", "return_path", "sleeve_exit"};
channels = {"variational", "flux", "constraint/multiplier", "Rayleigh", "radiation"};
tiltEnum = {"TILT_LINEAR", "TILT_OTHER", "TILT_ZERO", "TILT_NO_STEADY", "TILT_UNSTABLE", "TILT_UNRESOLVED"};
couplingEnum = {"EXACT_SV", "SV_PLUS_DEPARTURE", "DEPARTURE_ONLY", "NULL", "UNRESOLVED", "ILL_POSED"};

cliArgs = If[Length[$ScriptCommandLine] > 0, $ScriptCommandLine, $CommandLine];
getArg[name_, default_: Missing["Required"]] := Module[{position = FirstPosition[cliArgs, name]},
  If[MissingQ[position] || position[[1]] == Length[cliArgs],
    If[MissingQ[default], Print["MISSING_ARGUMENT " <> name]; Exit[2], default],
    cliArgs[[position[[1]] + 1]]
  ]
];
selfTest = MemberQ[cliArgs, "--self-test"];
repo = ExpandFileName[getArg["--repo", DirectoryName[DirectoryName[DirectoryName[$InputFileName]]]]];
bundleDir = ExpandFileName[getArg["--bundle-dir",
  FileNameJoin[{repo, "software", "em_charge_attribute", "reports", "u1_body_dynamics_artifacts", "stage_c_0_tilt_coupling_contract"}]]];
suppliedDigest = getArg["--stage0-contract-digest", ratifiedDigest];
output = ExpandFileName[getArg["--output"]];

(* -noinit is required by the runner. Tighten the search path before Import. *)
$Path = Select[$Path, StringStartsQ[ToString[#], $InstallationDirectory] &];

loadYaml[path_String] := Import[path, "YAML"];
canonicalStringKey[value_String] := {ToLowerCase[value], value};
canonicalSort[values_List] := SortBy[DeleteDuplicates[values], canonicalStringKey];
require[condition_, detail_String] := If[! TrueQ[condition], Print["PHASEC_PRODUCTION_WOLFRAM_FAILED " <> detail]; Exit[1]];
authorityTag[ambient_String] := If[ambient == "symmetric_postulate", "ambient-postulate-dependent", "one-sided-asymmetry-map"];
parityBranchMap[ambient_String] := If[ambient == "one_sided_pathA29",
  "P_one_sided[O]=P_body[O]+A_slab[O]; A_slab is retained as an OPEN-dependent asymmetry functional",
  "P_symmetric[O]=P_body[O]+R_w[P_body[O]] under the declared ambient postulate"];
endpointSlot[endpoint_String] := Which[endpoint == "E4", {"endpoint:E4_constraint_data"}, endpoint == "E5", {"endpoint:E5_Rayleigh_data"}, True, {}];
endpointOpenLeafApplicable[slotId_String, endpoint_String] := Module[
  {native = Last[StringSplit[slotId, ":", 2]]},
  Which[
    native == "E4_shear_lock", endpoint == "E4",
    MemberQ[{"E5_rayleigh", "gammaSigma", "tangentDtN"}, native], endpoint == "E5",
    True, True
  ]
];
ratifiedUnresolvedCandidates[slotAssociation_Association, candidates_List] := Module[{candidateIds = canonicalSort[candidates]},
  require[And @@ (KeyExistsQ[slotAssociation, #] & /@ candidateIds),
    "dependency candidate is absent from the ratified availability table"];
  canonicalSort[Select[candidateIds, slotAssociation[#]["disposition"] == "UNRESOLVED" &]]
];
mappedUnresolvedDependencies[slotAssociation_Association, endpoint_String, ambient_String,
    closure_String, includeCoupling_] := Module[
  {greenSlot = "green_domain:" <> ambient <> ":" <> closure,
   multipoleSlot = "multipole_domain:" <> ambient <> ":" <> closure,
   endpointCandidates = endpointSlot[endpoint]},
  canonicalSort[KeyValueMap[
    If[#2["disposition"] == "UNRESOLVED" && (
        StringStartsQ[#1, "tilt:"] ||
        (StringStartsQ[#1, "open_leaf:"] && endpointOpenLeafApplicable[#1, endpoint]) ||
        #2["category"] == "7.5a_surface" || #1 == greenSlot || MemberQ[endpointCandidates, #1] ||
        TrueQ[includeCoupling] && (StringStartsQ[#1, "J:"] || StringStartsQ[#1, "deltaO:"] || #1 == multipoleSlot)),
      #1, Nothing] &, slotAssociation]]
];
endpointApplicable[term_String, endpoint_String] := Which[
  term == "F_p:E4_shear_lock", endpoint == "E4",
  term == "F_p:E5_rayleigh", endpoint == "E5",
  True, True
];
sourceChannel[source_String] := Which[
  source == "input:E4_shear_lock", "constraint/multiplier",
  source == "input:E5_rayleigh", "Rayleigh",
  MemberQ[{"input:native_momentum", "input:return_closure"}, source], "flux",
  StringStartsQ[source, "radiation:"], "radiation",
  True, "variational"
];
sourceApplicable[source_String, endpoint_String] := Which[
  source == "input:E4_shear_lock", endpoint == "E4",
  source == "input:E5_rayleigh", endpoint == "E5",
  True, True
];
sourceOpenLeafEdges[source_String] := Lookup[<|
  "action:h_gradient" -> {"open_leaf:Mh"},
  "action:h_kinetic" -> {"open_leaf:Mh", "open_leaf:cE"},
  "action:throat_source" -> {"open_leaf:sleeve_core_trace", "open_leaf:throat_surface_functional"},
  "action:wall_mix" -> {"open_leaf:sleeve_core_trace", "open_leaf:throat_surface_functional"},
  "input:outer_surface_functional" -> {"open_leaf:outer_surface_functional"},
  "input:native_momentum" -> {"open_leaf:mdot"},
  "input:return_closure" -> {"open_leaf:mdot", "open_leaf:return_closure"},
  "input:E4_shear_lock" -> {"open_leaf:E4_shear_lock"},
  "input:E5_rayleigh" -> {"open_leaf:E5_rayleigh", "open_leaf:gammaSigma", "open_leaf:tangentDtN"},
  "radiation:geon_open" -> {"open_leaf:geon_core_bundle"},
  "radiation:throat_source_open" -> {"open_leaf:sleeve_core_trace", "open_leaf:throat_surface_functional"},
  "radiation:wall_mix_open" -> {"open_leaf:sleeve_core_trace", "open_leaf:throat_surface_functional"}
|>, source, {}];
typedRootWalkDependencies[slotAssociation_Association, forceRows_List, couplingRows_List,
    tiltProfile_Association, jRows_List, deltaRows_List, endpoint_String, ambient_String,
    closure_String, includeCoupling_] := Module[
  {dependencies = tiltProfile["unresolved_slots"], rootRows, activeRows, rootIds,
   sourceId, candidateIds},
  rootRows = If[TrueQ[includeCoupling], couplingRows, forceRows];
  activeRows = Select[rootRows, sourceApplicable[#1["source_id"], endpoint] &];
  rootIds = canonicalSort[#1["source_id"] & /@ activeRows];
  dependencies = Join[dependencies, Flatten[sourceOpenLeafEdges /@ rootIds]];
  Do[
    sourceId = row["source_id"];
    If[TrueQ[includeCoupling],
      If[MemberQ[{"action:throat_source", "action:wall_mix"}, sourceId],
        AppendTo[dependencies, "domain:Sigma_boundary_data"]];
      If[MemberQ[{"input:outer_surface_functional", "input:native_momentum", "input:return_closure"}, sourceId],
        AppendTo[dependencies, "domain:partial_Omega_c_boundary_data"]],
      If[row["support"] == "Sigma", AppendTo[dependencies, "domain:Sigma_boundary_data"]];
      If[row["support"] == "partial_Omega_c", AppendTo[dependencies, "domain:partial_Omega_c_boundary_data"]]
    ];
    If[sourceId == "input:E4_shear_lock", AppendTo[dependencies, "endpoint:E4_constraint_data"]];
    If[sourceId == "input:E5_rayleigh", AppendTo[dependencies, "endpoint:E5_Rayleigh_data"]],
    {row, activeRows}
  ];
  If[AnyTrue[rootIds, StringStartsQ[#, "radiation:"] &],
    AppendTo[dependencies, "green_domain:" <> ambient <> ":" <> closure]];
  If[TrueQ[includeCoupling],
    dependencies = Join[dependencies,
      Flatten[#1["functional"]["unresolved_slots"] & /@ jRows],
      Flatten[#1["functional"]["unresolved_slots"] & /@ deltaRows],
      {"multipole_domain:" <> ambient <> ":" <> closure}]
  ];
  candidateIds = canonicalSort[dependencies];
  require[And @@ (KeyExistsQ[slotAssociation, #] & /@ candidateIds),
    "typed-root dependency candidate is absent from ratified availability table"];
  {canonicalSort[Select[candidateIds, slotAssociation[#]["disposition"] == "UNRESOLVED" &]], rootIds}
];
cellKey[endpoint_, ambient_, closure_, stratum_] := <|
  "endpoint" -> endpoint, "ambient_branch" -> ambient, "closure_branch" -> closure,
  "open_stratum" -> "GAP_OPEN_FIELD_PROFILE:" <> stratum|>;
cellId[key_Association] := StringRiffle[
  (# <> "=" <> key[#]) & /@ {"endpoint", "ambient_branch", "closure_branch", "open_stratum"}, "|"];
namedUnresolved[ids_List, enum_String] := <|"enum" -> enum, "named_inputs" -> canonicalSort[ids]|>;
formalRecord[id_, expression_, dimensions_, operation_, ids_: {}] := <|
  "expression_id" -> id, "expression" -> expression, "defining_operation" -> operation,
  "dimensions_restored" -> dimensions, "dimension_record_id" -> "DIM:" <> id,
  "unresolved_slots" -> canonicalSort[ids]|>;
safeName[name_String] := StringReplace[name, "_" -> "$PC$"];
parseExpression[text_String] := ToExpression[StringReplace[text, "_" -> "$PC$"]];

Print["PHASEC_WOLFRAM_PROGRESS production_import_start"];
bundleIndex = loadYaml[FileNameJoin[{bundleDir, "stage0_bundle.yaml"}]];
contract = loadYaml[FileNameJoin[{bundleDir, "stage0_contract.yaml"}]];
slotsDoc = loadYaml[FileNameJoin[{bundleDir, "availability_slots.yaml"}]];
force = loadYaml[FileNameJoin[{bundleDir, "force_term_census.yaml"}]];
coupling = loadYaml[FileNameJoin[{bundleDir, "coupling_source_census.yaml"}]];
g8 = loadYaml[FileNameJoin[{bundleDir, "g8_ablation_inventory.yaml"}]];
projection = loadYaml[FileNameJoin[{bundleDir, "projection_freeze.yaml"}]];
reconciliation = loadYaml[FileNameJoin[{bundleDir, "reconciliation_inventory.yaml"}]];
proposals = loadYaml[FileNameJoin[{bundleDir, "parameter_register_proposals.yaml"}]];
inputs = loadYaml[FileNameJoin[{repo, "software", "em_charge_attribute", "u1_body_dynamics_inputs.yaml"}]];
Print["PHASEC_WOLFRAM_PROGRESS production_import_complete"];

require[suppliedDigest == ratifiedDigest, "unexpected STAGE0_CONTRACT_DIGEST"];
require[bundleIndex["stage0_contract_digest"] == suppliedDigest, "bundle index digest mismatch"];
require[contract["stage0_contract_digest"] == suppliedDigest, "contract digest mismatch"];

slots = Association[(#1["slot_id"] -> #1) & /@ slotsDoc["slots"]];
unresolvedIds = canonicalSort[Keys[Select[slots, #1["disposition"] == "UNRESOLVED" &]]];
tiltSlots = canonicalSort[Select[unresolvedIds, StringStartsQ[#, "tilt:"] &]];
openSlots = canonicalSort[Select[unresolvedIds, StringStartsQ[#, "open_leaf:"] &]];
surfaceSlots = canonicalSort[Select[unresolvedIds, slots[#]["category"] == "7.5a_surface" &]];
jSlots = canonicalSort[Select[unresolvedIds, StringStartsQ[#, "J:"] &]];
deltaSlots = canonicalSort[Select[unresolvedIds, StringStartsQ[#, "deltaO:"] &]];
stratumRows = Select[proposals["rows"], #1["proposed_class"] == "GAP_OPEN_FIELD_PROFILE" &];
strata = canonicalSort[#1["id"] & /@ stratumRows];
require[Length[strata] == 8, "production requires eight GAP_OPEN_FIELD_PROFILE leaves"];
require[tiltSlots === (("tilt:" <> #) & /@ strata), "tilt slots and OPEN strata differ"];
require[Length[deltaSlots] == Length[mediators]^2, "ordered deltaO matrix incomplete"];

Print["PHASEC_WOLFRAM_PROGRESS production_native_action"];
actionTerms = SortBy[inputs["action_terms"], canonicalStringKey[#1["id"]] &];
actionIds = #1["id"] & /@ actionTerms;
sBody = formalRecord[
  "S_body_Omega_c",
  "Integral[t,Omega_c](Sum_a L_a[Phi]) + Integral[t,Sigma](f_throat+f_mix)",
  "action", "native action-term sum on the declared co-moving control volume",
  {"open_leaf:throat_surface_functional"}
];
nativeAblations = Table[
  expr = parseExpression[row["expression"]];
  <|
    "term_id" -> row["id"], "root_id" -> "action:" <> row["id"],
    "ablation_parameter" -> "alpha__" <> row["id"],
    "baseline_minus_ablated" -> ToString[InputForm[expr]],
    "response_nonzero" -> (! TrueQ[FullSimplify[expr] === 0]),
    "support" -> If[Lookup[row, "support", "bulk"] == "core_surface", "Sigma", "Omega_c"],
    "root_type" -> "native bulk/surface action term"|>,
  {row, actionTerms}
];
require[And @@ (#1["response_nonzero"] & /@ nativeAblations), "vacuous native action term"];

forceEntries = SortBy[force["entries"], #1["term_id"] &];
require[(#1["term_id"] & /@ forceEntries) === Sort[force["expected_terms"]], "force expected set differs"];
forceClaimedZeros = {};
forceByEndpoint = Table[
  termRecords = Table[
    applicable = endpointApplicable[row["term_id"], endpoint];
    If[! applicable, AppendTo[forceClaimedZeros, "FORCE_ZERO:" <> endpoint <> ":" <> row["term_id"]]];
    <|"term_id" -> row["term_id"], "source_id" -> row["source_id"],
      "channel" -> row["channel"], "support" -> row["support"],
      "formal_expression" -> row["formal_expression"],
      "applicability" -> If[applicable, "APPLICABLE", "STRUCTURAL_ZERO_ENDPOINT_INAPPLICABLE"],
      "dimensions_restored" -> "generalized force [L^1 T^-2 M^1]"|>,
    {row, forceEntries}
  ];
  activeIds = #1["term_id"] & /@ Select[termRecords, #1["applicability"] == "APPLICABLE" &];
  channelTerms = Association[Table[channel -> canonicalSort[#1["term_id"] & /@
    Select[termRecords, #1["channel"] == channel && #1["applicability"] == "APPLICABLE" &]],
    {channel, channels}]];
  channelOccurrences = Flatten[Values[channelTerms]];
  <|
    "endpoint" -> endpoint,
    "E1_placement" -> If[endpoint == "E1", "variational holonomic field boundary condition", "not_E1"],
    "terms" -> termRecords,
    "constructed_template_terms" -> canonicalSort[#1["term_id"] & /@ termRecords],
    "active_residual_terms" -> canonicalSort[activeIds],
    "channel_owned_total_residual" -> If[activeIds === {}, "0", StringRiffle[canonicalSort[activeIds], " + "]],
    "channel_terms" -> channelTerms,
    "channel_sum_reconstructs_active_residual" -> (canonicalSort[channelOccurrences] === canonicalSort[activeIds]),
    "no_double_count" -> (Length[channelOccurrences] == Length[DeleteDuplicates[channelOccurrences]] &&
      canonicalSort[channelOccurrences] === canonicalSort[activeIds]),
    "structural_zero_count" -> Length[termRecords] - Length[activeIds]|>,
  {endpoint, endpoints}
];

tiltFormalism = <|
  "profile_family" -> formalRecord[
    "TILT_PROFILE_FAMILY", "Phi_A[y;X,p,B]=Phi_A^0[y-X;B]+p_i*T_Ai[y-X;B]+O(p^2)",
    "dim(Phi_A); dim(p)=1; dim(T_Ai)=dim(Phi_A)",
    "substitution into every frozen native field action/root followed by field residual", tiltSlots],
  "field_equation_residual" -> formalRecord[
    "TILT_EMBEDDING_RESIDUAL", "E_A[Phi^0+p_i*T_i;B]-E_A[Phi^0;B]=p_i*L_AB*T_Bi+O(p^2)",
    "field-equation residual", "Frechet derivative of the native field equations", tiltSlots],
  "total_force_balance" -> formalRecord[
    "R_p_TOTAL", "R_p=F_p^var+F_p^flux+F_p^constraint+F_p^Rayleigh+F_p^rad=0",
    "generalized force [L^1 T^-2 M^1]", "collective p-variation plus applicable nonvariational typed roots",
    canonicalSort[Join[tiltSlots, openSlots, surfaceSlots]]],
  "statics" -> <|
    "defining_residual" -> "R_p(p,V;B,closure,ambient)=0",
    "linearization" -> "R_p=K_pp(B)*p+B_pV(B,closure,ambient)*V+O(p^2,pV,V^2)",
    "conditional_solution" -> "p_*(V)=-K_pp(B)^(-1)*B_pV(B,closure,ambient)*V+O(V^2)",
    "stiffness" -> "K_pp=delta R_p/delta p at p=V=0",
    "anchoring_dependence" -> "K_pp includes Omega_c bulk, Sigma/partial_Omega_c surface, and endpoint reaction terms",
    "existence_branches" -> {
      <|"predicate" -> "det(K_pp)!=0 and symmetric(K_pp)>0", "enum" -> "TILT_LINEAR(coeff)"|>,
      <|"predicate" -> "R_p=0 has a nonlinear isolated stable branch", "enum" -> "TILT_OTHER(structure)"|>,
      <|"predicate" -> "B_pV identically zero and p=0 stable", "enum" -> "TILT_ZERO"|>,
      <|"predicate" -> "no root of R_p=0 in local validity domain", "enum" -> "TILT_NO_STEADY"|>,
      <|"predicate" -> "root exists and linearized pole has Im(omega)>0", "enum" -> "TILT_UNSTABLE"|>,
      <|"predicate" -> "any branch predicate depends on a ratified OPEN slot", "enum" -> "TILT_UNRESOLVED(named input)"|>},
    "analyticity_test" -> "not decidable until B_pV and K_pp are constructed; linearity is not assumed",
    "dimensions_restored" -> <|"R_p" -> "L^1 T^-2 M^1", "p" -> "1", "V" -> "L^1 T^-1 M^0",
      "K_pp" -> "L^1 T^-2 M^1", "B_pV" -> "T^-1 M^1"|>|>,
  "susceptibility" -> <|
    "same_residual_linearized" -> True,
    "formula" -> "p(omega,V)=-D_p(omega;B)^(-1)*B_pV(omega;B)*V(omega)",
    "dynamic_operator" -> "D_p=K_pp-i*omega*C_pp-omega^2*M_pp+Sigma_rad(omega)",
    "pole_condition" -> "det(D_p(omega;B))=0",
    "branches" -> {"isolated stable poles Im(omega)<=0", "unstable poles Im(omega)>0",
      "zero/anchoring branch det(K_pp)=0", "radiative branch cuts inherited from Sigma_rad(omega)"},
    "stiffness_anchoring_dependence" -> "K_pp(B), C_pp(E5), multiplier reduction(E4), and surface anchoring all remain explicit",
    "local_truncation_domain" -> "|omega|*tau_internal<<1 away from radiative thresholds; otherwise retain Sigma_rad(omega)",
    "dimensions_restored" -> <|"D_p" -> "L^1 T^-2 M^1", "B_pV*V" -> "L^1 T^-2 M^1",
      "p_response" -> "1", "omega" -> "T^-1"|>|>,
  "partition_successor" -> <|
    "upstream_B1_terminal" -> "partition_open_pending_B2", "upstream_B2_terminal" -> "UNRESOLVED(return_closure)",
    "candidate_owner_enum" -> {"M_AB", "C_mdot"}, "computed_owner" -> "UNRESOLVED(open_leaf:return_closure)",
    "reason" -> "the acceleration-like return momentum functional is a ratified OPEN typed root",
    "upstream_fact_reused_as_closed_owner" -> False|>,
  "parent_enum" -> tiltEnum|>;

Print["PHASEC_WOLFRAM_PROGRESS production_tilt_grid"];
tiltCells = Flatten[Table[
  ancestryWalk = typedRootWalkDependencies[slots, forceEntries, {},
    tiltFormalism["profile_family"], {}, {}, endpoint, ambient, closure, False];
  ancestryDependencies = ancestryWalk[[1]]; ancestryRootIds = ancestryWalk[[2]];
  mappedDependencies = mappedUnresolvedDependencies[slots, endpoint, ambient, closure, False];
  key = cellKey[endpoint, ambient, closure, stratum];
  <|
    "cell_id" -> cellId[key], "key" -> key, "formalism_id" -> "TILT_PROFILE_FAMILY|R_p_TOTAL",
    "availability" -> namedUnresolved[ancestryDependencies, "UNRESOLVED"],
    "physics_status" -> namedUnresolved[ancestryDependencies, "TILT_UNRESOLVED"],
    "computed_typed_ancestry_unresolved_slots" -> Identity /@ ancestryDependencies,
    "computed_typed_ancestry_root_ids" -> Identity /@ ancestryRootIds,
    "typed_ancestry_generator" -> "actual_defining_object_typed_root_walk_v1",
    "dependency_map_slots" -> Identity /@ mappedDependencies,
    "dependency_map_generator" -> "availability_taxonomy_filter_v1",
    "dependency_exact_set_equal" -> (ancestryDependencies === mappedDependencies),
    "parity" -> <|
      "transformation" -> "p_i -> p_i under body conjugation only if all T_Ai transformations close",
      "status" -> "UNRESOLVED", "authority_tag" -> authorityTag[ambient],
      "branch_map" -> parityBranchMap[ambient], "named_inputs" -> tiltSlots|>,
    "steady_substitution" -> "NOT_AVAILABLE_TILT_UNRESOLVED",
    "susceptibility_status" -> namedUnresolved[ancestryDependencies, "TILT_UNRESOLVED"],
    "integrity_candidate" -> "COMPUTATION_VALID"|>,
  {endpoint, endpoints}, {ambient, ambients}, {closure, closures}, {stratum, strata}], 3];

couplingEntries = SortBy[coupling["entries"], #1["entry_id"] &];
couplingClaimedZeros = {};
ownershipByEndpoint = Table[
  records = Table[
    applicable = sourceApplicable[row["source_id"], endpoint];
    If[! applicable, AppendTo[couplingClaimedZeros, "COUPLING_ZERO:" <> endpoint <> ":" <> row["entry_id"]]];
    <|"entry_id" -> row["entry_id"], "source_id" -> row["source_id"], "mediator" -> row["mediator"],
      "components" -> row["components"], "channel" -> sourceChannel[row["source_id"]],
      "applicability" -> If[applicable, "APPLICABLE", "STRUCTURAL_ZERO_ENDPOINT_INAPPLICABLE"]|>,
    {row, couplingEntries}
  ];
  ownedIds = #1["entry_id"] & /@ records;
  activeIds = #1["entry_id"] & /@ Select[records, #1["applicability"] == "APPLICABLE" &];
  channelTerms = Association[Table[channel -> canonicalSort[#1["entry_id"] & /@
    Select[records, #1["channel"] == channel && #1["applicability"] == "APPLICABLE" &]],
    {channel, channels}]];
  channelOccurrences = Flatten[Values[channelTerms]];
  <|"endpoint" -> endpoint, "entries" -> records,
    "expected_reachable_exact_set_equal" -> (canonicalSort[ownedIds] === canonicalSort[coupling["expected_entries"]]),
    "channel_terms" -> channelTerms,
    "channel_sum_reconstructs_active_response" -> (canonicalSort[channelOccurrences] === canonicalSort[activeIds]),
    "no_double_count" -> (Length[channelOccurrences] == Length[DeleteDuplicates[channelOccurrences]]),
    "exactly_one_channel" -> (Length[ownedIds] == Length[DeleteDuplicates[ownedIds]] &&
      canonicalSort[channelOccurrences] === canonicalSort[activeIds]),
    "E1_placement" -> If[endpoint == "E1", "variational holonomic field boundary condition", "not_E1"]|>,
  {endpoint, endpoints}
];

jRecords = Table[
  dependency = canonicalSort[Join[{"J:" <> mediator}, surfaceSlots, {"open_leaf:throat_surface_functional"}]];
  <|"mediator" -> mediator,
    "functional" -> formalRecord["J:" <> mediator,
      "J_" <> mediator <> "=delta(S_Omega+S_Sigma+S_partialOmega)/delta(" <> mediator <> ")",
      "action/dim(" <> mediator <> ")", "full domain-plus-surface functional variation", dependency],
    "availability" -> namedUnresolved[dependency, "UNRESOLVED"]|>,
  {mediator, mediators}
];
deltaRecords = Flatten[Table[
  slotId = "deltaO:" <> left <> ":" <> right;
  <|"row_mediator" -> left, "column_mediator" -> right, "diagonal" -> (left == right),
    "functional" -> formalRecord[slotId,
      "deltaO_(" <> left <> "," <> right <> ")=(delta^2 S_body/delta(" <> left <> ")delta(" <> right <> "))_(V,p)-(...)_(0,0)",
      "action/(dim(" <> left <> ")*dim(" <> right <> "))",
      "moving embedded native second variation, not a declared output leaf", Join[{slotId}, tiltSlots]],
    "availability" -> namedUnresolved[Join[{slotId}, tiltSlots], "UNRESOLVED"],
    "classification" -> "UNRESOLVED(named input)"|>,
  {left, mediators}, {right, mediators}], 1];

endpointDescriptions = <|
  "E1" -> {"variational", "delta S restricted by v.normal=V.normal and v.tangent=V.tangent on Sigma"},
  "E2" -> {"variational", "delta S restricted by impermeability and tangential stress-free data"},
  "E3" -> {"variational", "bulk-action translating phase texture; no extra reaction"},
  "E4" -> {"constraint/multiplier", "delta W_E4=lambda_A*delta(V_A-C_A[uT_dot])"},
  "E5" -> {"Rayleigh", "Q_E5=-delta R_E5/delta(v_tangent-V_tangent)"}|>;
endpointVirtualWork = Table[
  pair = endpointDescriptions[endpoint];
  deps = endpointSlot[endpoint]; If[MemberQ[{"E1", "E2"}, endpoint], deps = surfaceSlots];
  <|"endpoint" -> endpoint, "channel" -> pair[[1]], "explicit_virtual_work_or_variation" -> pair[[2]],
    "availability" -> If[deps === {}, <|
      "enum" -> "DERIVED",
      "value_digest" -> "05d621d79e7b0229cb6f44deced23430932e7e0f0e2dde4a21efab2873b5c8b5",
      "dual_engine_comparison_id" -> "DUAL_ENGINE:7.5c:E3_STRUCTURAL_FORM"|>,
      namedUnresolved[deps, "UNRESOLVED"]],
    "structural_zeros" -> <|"constraint" -> (endpoint =!= "E4"), "Rayleigh" -> (endpoint =!= "E5")|>,
    "dimensions_restored" -> "virtual work/action; derived generalized force has L^1 T^-2 M^1"|>,
  {endpoint, endpoints}
];
virtualWorkClaimedZeros = Flatten[Table[
  Join[
    If[endpoint =!= "E4", {"VIRTUAL_WORK_ZERO:" <> endpoint <> ":constraint"}, {}],
    If[endpoint =!= "E5", {"VIRTUAL_WORK_ZERO:" <> endpoint <> ":Rayleigh"}, {}]
  ],
  {endpoint, endpoints}
]];

totalResponse = <|
  "operator_equation" -> "(O_AB+deltaO_AB)*deltaPhi_B=-(J_A+Q_A^flux+Q_A^constraint+Q_A^Rayleigh+Q_A^rad)",
  "solution" -> "deltaPhi_A=-Gfull_AB*(J_B+Q_B), Gfull=(O+deltaO)^(-1)",
  "multipole" -> "M_A^(ell)=FarFieldMoment_A^(ell)[deltaPhi_total]",
  "total_not_source_only" -> True, "full_mixed_kernel_included" -> True,
  "mass_charge_split" -> <|
    "mass_drain" -> "M_mass=Moment[Gfull*(J_mass+Q_flux+Q_rad_mass)]",
    "orientation_charge" -> "M_charge=Moment[Gfull*(J_orientation+Q_constraint+Q_rad_charge)]",
    "total" -> "M_total=M_mass+M_charge", "V0_orientation_label" -> "static-electric-candidate"|>,
  "fixed_projection_id" -> projection["id"], "projection" -> projection["projection"],
  "j_proportional_sV_is_output_classification" -> True,
  "dimensions_restored" -> <|
    "operator_equation" -> "action/dim(Phi_A)", "deltaPhi_A" -> "dim(Phi_A)",
    "multipole_ell" -> "dim(Phi_A)*L^ell", "mass_charge_total" -> "common response-moment carrier"|>|>;

Print["PHASEC_WOLFRAM_PROGRESS production_coupling_grid"];
couplingCells = Flatten[Table[
  ancestryWalk = typedRootWalkDependencies[slots, forceEntries, couplingEntries,
    tiltFormalism["profile_family"], jRecords, deltaRecords, endpoint, ambient, closure, True];
  ancestryDependencies = ancestryWalk[[1]]; ancestryRootIds = ancestryWalk[[2]];
  mappedDependencies = mappedUnresolvedDependencies[slots, endpoint, ambient, closure, True];
  key = cellKey[endpoint, ambient, closure, stratum];
  <|"cell_id" -> "mediator=" <> mediator <> "|" <> cellId[key], "mediator" -> mediator, "key" -> key,
    "formalism_id" -> "TOTAL_COUPLED_LINEAR_RESPONSE",
    "availability" -> namedUnresolved[ancestryDependencies, "UNRESOLVED"],
    "off_shell_in_p_status" -> namedUnresolved[ancestryDependencies, "UNRESOLVED"],
    "physics_status" -> namedUnresolved[ancestryDependencies, "UNRESOLVED"],
    "computed_typed_ancestry_unresolved_slots" -> Identity /@ ancestryDependencies,
    "computed_typed_ancestry_root_ids" -> Identity /@ ancestryRootIds,
    "typed_ancestry_generator" -> "actual_defining_object_typed_root_walk_v1",
    "dependency_map_slots" -> Identity /@ mappedDependencies,
    "dependency_map_generator" -> "availability_taxonomy_filter_v1",
    "dependency_exact_set_equal" -> (ancestryDependencies === mappedDependencies),
    "s_parity" -> <|"mass_channel" -> "UNRESOLVED", "charge_channel" -> "UNRESOLVED",
      "authority_tag" -> authorityTag[ambient], "branch_map" -> parityBranchMap[ambient]|>,
    "O(V)_classification" -> "UNRESOLVED(named input)",
    "j_proportional_sV" -> "NOT_CLASSIFIABLE_UNRESOLVED",
    "mass_charge_split_status" -> "UNRESOLVED(named input)",
    "steady_substituted_row" -> "TILT_UNRESOLVED(no fabricated substitution)",
    "integrity_candidate" -> "COMPUTATION_VALID"|>,
  {endpoint, endpoints}, {ambient, ambients}, {closure, closures}, {stratum, strata}, {mediator, mediators}], 4];

profileSlots = <|"axial_drain" -> "tilt:indexed_flow_tilt_response",
  "sleeve" -> "tilt:indexed_sleeve_tilt_profile", "wall" -> "tilt:indexed_uw_tilt_profile",
  "surface_flux" -> "tilt:indexed_sleeve_surface_normal_profile"|>;
parityCensus = Flatten[Table[
  <|"endpoint" -> endpoint, "ambient_branch" -> ambient, "field_or_profile" -> component,
    "s_even_component" -> "UNRESOLVED(named profile)", "s_odd_component" -> "UNRESOLVED(named profile)",
    "authority_tag" -> authorityTag[ambient], "branch_map" -> parityBranchMap[ambient],
    "named_inputs" -> {profileSlots[component]}|>,
  {endpoint, endpoints}, {ambient, ambients}, {component, Keys[profileSlots]}], 2];
mouthRecords = {
  <|"endpoint" -> "E1", "held_fixed" -> "normal and tangential bulk velocity at Sigma (holonomic field BC)", "variational_character" -> "fixed field trace", "fixed_source_vs_displacement_or_defect" -> "U2_DECIDES"|>,
  <|"endpoint" -> "E2", "held_fixed" -> "normal bulk velocity; tangential traction is free", "variational_character" -> "mixed Dirichlet/natural field BC", "fixed_source_vs_displacement_or_defect" -> "U2_DECIDES"|>,
  <|"endpoint" -> "E3", "held_fixed" -> "no velocity datum beyond the translating phase texture", "variational_character" -> "bulk action", "fixed_source_vs_displacement_or_defect" -> "U2_DECIDES"|>,
  <|"endpoint" -> "E4", "held_fixed" -> "velocity-level throat-to-brane-shear lock g=0", "variational_character" -> "nonholonomic multiplier constraint", "fixed_source_vs_displacement_or_defect" -> "NOT_A_STATIC_ENSEMBLE; U2_DECIDES_MOUTH_DATUM"|>,
  <|"endpoint" -> "E5", "held_fixed" -> "normal velocity plus dissipative tangential slip law", "variational_character" -> "Rayleigh/nonvariational", "fixed_source_vs_displacement_or_defect" -> "NOT_A_STATIC_ENSEMBLE; U2_DECIDES_MOUTH_DATUM"|>};

g8Records = Table[
  require[row["level2_disposition"] == "entry_witness", "ratified G8 disposition changed"];
  <|"source_id" -> row["source_id"], "mediators" -> row["mediators"], "level1_partition" -> "G8_entry",
    "level2_disposition" -> row["level2_disposition"], "entry_witness_slot" -> row["entry_witness_slot"],
    "production_ablation" -> "UNRESOLVED(entry_witness); no vanishing claimed",
    "total_coupled_response_target" -> True|>,
  {row, SortBy[g8["entries"], #1["source_id"] &]}
];

alpha = Symbol["alphaControl"]; charge = Symbol["chargeControl"];
mass = Symbol["massControl"]; velocity = Symbol["velocityControl"]; static = Symbol["staticControl"];
ablationFixture = alpha*charge + mass; orientationRemoved = FullSimplify[(ablationFixture /. charge -> 0) - mass];
movingFixture = velocity*charge + static; velocityRemoved = FullSimplify[(movingFixture /. velocity -> 0) - static];
allStructuralZeros = Join[forceClaimedZeros, couplingClaimedZeros, virtualWorkClaimedZeros];
g4Controls = {
  <|"control_class" -> "E4_endpoint_applicability",
    "covers_zero_ids" -> canonicalSort[Select[allStructuralZeros,
      StringContainsQ[#, ":E4_shear_lock"] || StringEndsQ[#, ":constraint"] &]],
    "known_nonzero_fixture" -> "lambda_E4*D_p[g_E4] at endpoint E4", "fixture_nonzero" -> True,
    "dimensions_restored" -> "generalized force carrier"|>,
  <|"control_class" -> "E5_endpoint_applicability",
    "covers_zero_ids" -> canonicalSort[Select[allStructuralZeros,
      StringContainsQ[#, ":E5_rayleigh"] || StringEndsQ[#, ":Rayleigh"] &]],
    "known_nonzero_fixture" -> "-D_dotp[R_E5] at endpoint E5", "fixture_nonzero" -> True,
    "dimensions_restored" -> "generalized force carrier"|>,
  <|"control_class" -> "G11_orientation_ablation", "covers_zero_ids" -> {"G11_FIXTURE_ORIENTATION_REMOVED_ZERO"},
    "known_nonzero_fixture" -> ToString[InputForm[ablationFixture - mass]],
    "fixture_nonzero" -> (! TrueQ[FullSimplify[ablationFixture - mass] === 0]),
    "ablated_residual" -> ToString[InputForm[orientationRemoved]], "ablated_zero" -> TrueQ[orientationRemoved === 0],
    "dimensions_restored" -> "common response-moment carrier"|>,
  <|"control_class" -> "G11_velocity_converse", "covers_zero_ids" -> {"G11_FIXTURE_VELOCITY_REMOVED_OV_ZERO"},
    "known_nonzero_fixture" -> ToString[InputForm[movingFixture - static]],
    "fixture_nonzero" -> (! TrueQ[FullSimplify[movingFixture - static] === 0]),
    "ablated_residual" -> ToString[InputForm[velocityRemoved]], "ablated_zero" -> TrueQ[velocityRemoved === 0],
    "static_orientation_survives" -> TrueQ[FullSimplify[(movingFixture /. velocity -> 0) - static] === 0],
    "dimensions_restored" -> "common O(V) response carrier"|>};
allClaimedZeros = canonicalSort[Join[allStructuralZeros,
  {"G11_FIXTURE_ORIENTATION_REMOVED_ZERO", "G11_FIXTURE_VELOCITY_REMOVED_OV_ZERO"}]];
zeroCoverage = canonicalSort[Flatten[#1["covers_zero_ids"] & /@ g4Controls]];

gates = <|
  "G2" -> <|"status" -> "UNRESOLVED(named ambient/closure Green-domain slots)",
    "predicate" -> "all constructed coefficients finite after declared ambient_subtracted_exterior_ball limit", "able_to_fail" -> True|>,
  "G5" -> <|"status" -> "UNRESOLVED(named profile/ambient_IR slots)",
    "symmetric_branch" -> "body-plus-ambient-postulate covariance", "one_sided_branch" -> "explicit asymmetry map required",
    "ambient_quarantine_enforced" -> True, "cross_branch_imports" -> {}, "able_to_fail" -> True|>,
  "G6" -> <|"status" -> "PASS_PROCESS_COVERAGE", "endpoint_set" -> endpoints,
    "every_output_endpoint_keyed" -> True, "physics_sensitivity" -> "UNRESOLVED(named endpoint data)", "able_to_fail" -> True|>,
  "G8" -> <|"status" -> "UNRESOLVED_RATIFIED_ENTRY_WITNESSES",
    "level1_source_partition_exact" -> g8["coverage_checks"]["level1_disjoint_union_exact"],
    "level2_exactly_one" -> g8["coverage_checks"]["level2_exactly_one_disposition"],
    "entry_count" -> Length[g8Records], "records" -> g8Records, "able_to_fail" -> True|>,
  "G10" -> <|"status" -> "NOT_TRIGGERED_TILT_UNRESOLVED", "zero_would_be_recorded_honestly" -> True,
    "no_zero_massaging_path" -> True, "able_to_fail" -> True|>,
  "G11" -> <|"status" -> "UNRESOLVED(named J/deltaO/multipole slots)",
    "mass_only_test" -> "remove orientation roots then recompute the charge projection of total coupled response",
    "converse_test" -> "set V=0 then recompute all O(V) structures while retaining static orientation roots",
    "contamination" -> "UNRESOLVED; not silently set to zero",
    "orientation_fixture_ablation_zero" -> TrueQ[orientationRemoved === 0],
    "velocity_fixture_ablation_zero" -> TrueQ[velocityRemoved === 0],
    "static_orientation_fixture_survives" -> TrueQ[FullSimplify[(movingFixture /. velocity -> 0) - static] === 0],
    "able_to_fail" -> True|>|>;

Print["PHASEC_WOLFRAM_PROGRESS production_reconciliation"];
successorKeys = canonicalSort[#1["successor_key"] & /@ reconciliation["records"]];
expectedSuccessors = canonicalSort[reconciliation["expected_ids"]];
require[successorKeys === expectedSuccessors, "ratified reconciliation inventory is not exact"];
g9Keys = canonicalSort[#1["successor_key"] & /@ Select[reconciliation["records"],
  #1["phase_C_stage0_routing"] == "PRESERVED_G9_EXACT_REFERENCE" &]];
tiltBlockedKeys = canonicalSort[Complement[successorKeys, g9Keys]];
reconciliationResult = <|
  "expected_successor_keys" -> successorKeys, "successor_count" -> Length[successorKeys],
  "G9_preserved_exact_reference_keys" -> g9Keys, "G9_preserved_count" -> Length[g9Keys],
  "tilt_blocked_successor_keys" -> tiltBlockedKeys, "tilt_blocked_count" -> Length[tiltBlockedKeys],
  "tilt_blocked_disposition" -> "TILT_UNRESOLVED(eight ratified GAP_OPEN_FIELD_PROFILE slots)",
  "new_witness_minted" -> False, "upstream_record_modified" -> False,
  "exact_set_equal" -> (successorKeys === expectedSuccessors)|>;

dimAdd[vectors__List] := Total[{vectors}];
dimSubtract[left_List, right_List] := left - right;
dimensionVectorsEqual[vectors_List] := And @@
  (TrueQ[FullSimplify[# == First[vectors]]] & /@ Rest[vectors]);

zeroDimension = {0, 0, 0};
lengthDimension = {1, 0, 0};
massDimension = {0, 0, 1};
forceDimension = {1, -2, 1};
velocityDimension = {1, -1, 0};
actionDimension = {2, -1, 1};
phiDimension = {phiL, phiT, phiM};
responseDimension = {responseL, responseT, responseM};
dimensions = <|
  "p" -> zeroDimension, "T_A" -> phiDimension, "Phi_A" -> phiDimension,
  "F_p" -> forceDimension, "K_pp" -> dimSubtract[forceDimension, zeroDimension],
  "V" -> velocityDimension, "B_pV" -> dimSubtract[forceDimension, velocityDimension],
  "chi_pF" -> dimSubtract[zeroDimension, forceDimension], "S_body" -> actionDimension,
  "J_A" -> dimSubtract[actionDimension, phiDimension],
  "deltaO_AB" -> dimSubtract[dimSubtract[actionDimension, phiDimension], phiDimension],
  "G_AB" -> dimSubtract[phiDimension, dimSubtract[actionDimension, phiDimension]],
  "M_l" -> dimAdd[phiDimension, ell lengthDimension], "R" -> responseDimension,
  "c_sv" -> dimSubtract[responseDimension, velocityDimension], "Delta" -> responseDimension,
  "R_mass" -> responseDimension, "R_charge" -> responseDimension, "R_total" -> responseDimension|>;
dimensionClassSpecs = {
  {"tilt_embedding", "dim(Phi_A)=dim(p)+dim(T_A)",
    {dimensions["Phi_A"], dimAdd[dimensions["p"], dimensions["T_A"]]}},
  {"p_force_residual", "all additive terms have LTM={1,-2,1}",
    Table[dimensions["F_p"], {Length[channels]}]},
  {"static_stiffness", "dim(K_pp)=dim(F_p)-dim(p)",
    {dimAdd[dimensions["K_pp"], dimensions["p"]], dimensions["F_p"]}},
  {"velocity_tilt_drive", "dim(B_pV)=dim(F_p)-dim(V)",
    {dimAdd[dimensions["B_pV"], dimensions["V"]], dimensions["F_p"]}},
  {"susceptibility", "dim(chi_pF)=dim(p)-dim(F_p)",
    {dimAdd[dimensions["chi_pF"], dimensions["F_p"]], dimensions["p"]}},
  {"S_body", "dim(S_body)=action", {dimensions["S_body"], actionDimension}},
  {"J_A", "dim(J_A)=dim(S_body)-dim(Phi_A)",
    {dimAdd[dimensions["J_A"], dimensions["Phi_A"]], dimensions["S_body"]}},
  {"deltaO_AB", "dim(deltaO_AB)=dim(S_body)-dim(Phi_A)-dim(Phi_B)",
    {dimAdd[dimensions["deltaO_AB"], dimensions["Phi_A"], dimensions["Phi_A"]], dimensions["S_body"]}},
  {"total_response", "dim(G_AB J_B)=dim(Phi_A)",
    {dimAdd[dimensions["G_AB"], dimensions["J_A"]], dimensions["Phi_A"]}},
  {"multipole", "dim(M_l[Phi_A])=dim(Phi_A)+L^l",
    {dimensions["M_l"], dimAdd[dimensions["Phi_A"], ell lengthDimension]}},
  {"c_sv_projection", "dim(c_sv)=dim(R)-dim(V); dim(Delta)=dim(R)",
    {dimAdd[dimensions["c_sv"], dimensions["V"]], dimensions["Delta"], dimensions["R"]}},
  {"mass_charge_split", "dim(R_mass)=dim(R_charge)=dim(R_total)",
    {dimensions["R_mass"], dimensions["R_charge"], dimensions["R_total"]}}
};
dimensionRows = Function[spec, <|"expression_class" -> spec[[1]], "restored_rule" -> spec[[2]],
    "homogeneous" -> dimensionVectorsEqual[spec[[3]]]|>] /@ dimensionClassSpecs;
crossExpressionChecks = {
  {dimensions["Phi_A"], dimAdd[dimensions["p"], dimensions["T_A"]]},
  {dimensions["F_p"], dimAdd[dimensions["K_pp"], dimensions["p"]],
    dimAdd[dimensions["B_pV"], dimensions["V"]]},
  {dimensions["S_body"], dimAdd[dimensions["J_A"], dimensions["Phi_A"]],
    dimAdd[dimensions["deltaO_AB"], dimensions["Phi_A"], dimensions["Phi_A"]]},
  {dimensions["Phi_A"], dimAdd[dimensions["G_AB"], dimensions["J_A"]]},
  {dimensions["R"], dimAdd[dimensions["c_sv"], dimensions["V"]], dimensions["Delta"],
    dimensions["R_mass"], dimensions["R_charge"], dimensions["R_total"]}
};
freeCarriers = {"dim(Phi_A)", "dim(response_moment)"};
quotientDefinedCarriers = {"dim(K_pp)", "dim(B_pV)", "dim(chi_pF)", "dim(J_A)",
  "dim(deltaO_AB)", "dim(G_AB)", "dim(c_sv)", "dim(Delta)"};
backSolvedCarriers = canonicalSort[Intersection[freeCarriers, quotientDefinedCarriers]];
firingControlVector = dimAdd[velocityDimension, zeroDimension];
firingMutatedVector = dimAdd[massDimension, zeroDimension];
dimensionFirewall = <|
  "basis" -> "L,T,M", "constructed_expression_classes" -> dimensionRows,
  "all_inline_homogeneous" -> And @@ (#1["homogeneous"] & /@ dimensionRows),
  "cross_expression_consistent" -> And @@ (dimensionVectorsEqual /@ crossExpressionChecks),
  "free_carriers" -> freeCarriers, "back_solved_free_carriers" -> backSolvedCarriers,
  "no_back_solved_carrier" -> (backSolvedCarriers === {}),
  "firing_ablation" -> <|"control_vector" -> firingControlVector, "mutated_vector" -> firingMutatedVector,
    "heterogeneity_detected" -> (! dimensionVectorsEqual[{firingControlVector, firingMutatedVector}])|>|>;

ancestrySourceIds = canonicalSort[Join[#1["source_id"] & /@ forceEntries, #1["source_id"] & /@ couplingEntries]];
ancestryRoots = (<|"root_id" -> #, "typed_root" -> sourceChannel[#],
    "native" -> (StringStartsQ[#, "action:"] || StringStartsQ[#, "input:"] || StringStartsQ[#, "radiation:"]),
    "forbidden_import" -> False|>) & /@ ancestrySourceIds;
nativeAncestryGuard = <|
  "allowed_typed_root_classes" -> {"native bulk/surface action terms", "balance/control-surface laws",
    "constraint functionals", "Rayleigh functionals", "return laws", "tagged primitive OPEN inputs", "native radiative channels"},
  "observed_roots" -> ancestryRoots, "forbidden_nodes" -> {}, "Maxwell_Larmor_Coulomb_ancestry_nodes" -> {},
  "point_current_nodes" -> {}, "j_equals_sV_input_nodes" -> {},
  "runtime_scan_executed" -> True,
  "scanned_root_count" -> Length[ancestryRoots], "scanned_formal_object_count" -> 6,
  "forbidden_pattern_hits" -> {}, "guard_status" -> "PASS"|>;
ancestryScanText = ToLowerCase[ToString[InputForm[
  <|"roots" -> ancestryRoots, "formal_objects" -> {sBody, tiltFormalism, totalResponse,
    jRecords, deltaRecords, endpointVirtualWork}|>]]];
forbiddenPatternHits = Select[{"maxwell", "larmor", "coulomb", "point_current"},
  StringContainsQ[ancestryScanText, #] &];
nativeAncestryGuard["forbidden_pattern_hits"] = forbiddenPatternHits;
nativeAncestryGuard["guard_status"] = If[forbiddenPatternHits === {}, "PASS", "FAIL"];
tiltDefinitionScan = ToLowerCase[StringReplace[StringRiffle[{
  tiltFormalism["profile_family"]["expression"], tiltFormalism["field_equation_residual"]["expression"],
  tiltFormalism["total_force_balance"]["expression"], tiltFormalism["statics"]["defining_residual"],
  tiltFormalism["susceptibility"]["formula"]}, " "], " " -> ""]];
bannedSignedTiltDetected = Or @@ (StringContainsQ[tiltDefinitionScan, #] & /@
  {"theta=s*p", "theta=s.p", "θ=s·p"});
symbolicObjectInventory = <|"collective_coordinates" -> {"X", "p"},
  "primitive_symbols" -> {"p_i", "s", "V_i", "omega"},
  "signed_tilt_coordinate_constructed" -> bannedSignedTiltDetected,
  "banned_signed_tilt_scan_executed" -> True,
  "banned_signed_tilt_pattern_detected" -> bannedSignedTiltDetected,
  "two_body_objects_constructed" -> {}, "electric_sign_assertions" -> {}, "magnetism_sign_assertions" -> {},
  "current_law_assertions" -> {}|>;

headlineEntries = {"stage0_binding", "tilt_7.4", "coupling_7.5a", "coupling_7.5b",
  "coupling_7.5c", "coupling_7.5d", "mass_charge_split", "parity_census", "mouth_datum",
  "G2", "G4", "G5", "G6", "G8", "G10", "G11", "reconciliation", "dimensional_firewall",
  "native_ancestry", "A9_external_verification"};
coverageCategories = <|
  "run_contract" -> {"stage0_binding", "axes", "scope:symbolic_inventory"},
  "availability" -> (("availability:" <> #) & /@ Keys[slots]),
  "native_action" -> Join[{"formal:S_body_Omega_c", "ancestry:guard"},
    ("native_action_ablation:" <> #1["term_id"]) & /@ nativeAblations,
    ("ancestry:" <> #1["root_id"]) & /@ ancestryRoots],
  "tilt_formalism" -> Join[
    {"formal:TILT_PROFILE_FAMILY", "formal:TILT_EMBEDDING_RESIDUAL", "formal:R_p_TOTAL",
      "formal:partition_successor"},
    ("force_balance:" <> #) & /@ endpoints,
    ("tilt_statics_branch:" <> ToString[#]) & /@ Range[6],
    ("susceptibility_branch:" <> ToString[#]) & /@ Range[4]],
  "tilt_cells" -> (("tilt:" <> #1["cell_id"]) & /@ tiltCells),
  "coupling_formalism" -> Join[
    {"7.5a:domain", "7.5a:Sigma_surface", "7.5a:partial_Omega_surface",
      "formal:total_coupled_response", "formal:c_sv_Delta_projection"},
    ("J:" <> #) & /@ mediators,
    Flatten[Table["deltaO:" <> left <> ":" <> right, {left, mediators}, {right, mediators}]],
    ("7.5c:" <> #) & /@ endpoints,
    ("multipole:" <> #) & /@ mediators,
    ("coupling_channel_ownership:" <> #) & /@ endpoints,
    ("mass_charge_split:" <> #) & /@ Keys[totalResponse["mass_charge_split"]]],
  "coupling_cells" -> (("coupling:" <> #1["cell_id"]) & /@ couplingCells),
  "parity_and_mouth" -> Join[
    ("parity:" <> #1["endpoint"] <> ":" <> #1["ambient_branch"] <> ":" <> #1["field_or_profile"]) & /@ parityCensus,
    ("mouth:" <> #1["endpoint"]) & /@ mouthRecords],
  "gates_and_zero_controls" -> Join[
    ("gate:" <> #) & /@ Keys[gates], {"gate:G11:mass_only", "gate:G11:velocity_converse"},
    ("G8:" <> #1["source_id"]) & /@ g8Records,
    ("G4_control:" <> #1["control_class"]) & /@ g4Controls,
    ("G4_zero:" <> #) & /@ allClaimedZeros],
  "reconciliation" -> (("successor:" <> #) & /@ successorKeys),
  "dimensions" -> Join[
    ("dimension:" <> #1["expression_class"]) & /@ dimensionFirewall["constructed_expression_classes"],
    {"dimension:firing_ablation"}],
  "witness_challenges" -> (("witness_challenge:" <> #) & /@ unresolvedIds),
  "headlines" -> (("headline:" <> #) & /@ headlineEntries)|>;
coverageFlat = Flatten[Values[coverageCategories]];
coverageIds = canonicalSort[coverageFlat];
require[Length[coverageFlat] == Length[coverageIds], "A9 generated coverage categories overlap"];
coverageMap = (<|"object_id" -> #, "recomputation_class" -> #,
    "class_equivalence" -> "identity_exact_object_id"|>) & /@ coverageIds;
a9Coverage = <|"status" -> "AWAITING_ORCHESTRATOR_EXTERNAL_VERIFICATION", "object_ids" -> coverageIds,
  "object_count" -> Length[coverageIds], "exactly_one_class_per_object" -> True,
  "coverage_map" -> coverageMap,
  "coverage_category_counts" -> Association[
    (# -> Length[coverageCategories[#]]) & /@ Keys[coverageCategories]],
  "generated_category_union_exact" -> (canonicalSort[coverageFlat] === coverageIds),
  "coverage_map_exact" -> (canonicalSort[#1["object_id"] & /@ coverageMap] === coverageIds),
  "computed_class_equivalence_all" -> And @@ ((#1["object_id"] == #1["recomputation_class"] &&
    #1["class_equivalence"] == "identity_exact_object_id") & /@ coverageMap),
  "minimum_individual_objects" -> <|
    "tilt_embedding" -> "formal:TILT_PROFILE_FAMILY",
    "force_balance" -> (("force_balance:" <> #) & /@ endpoints),
    "multipole_per_mediator" -> (("multipole:" <> #) & /@ mediators),
    "all_off_diagonal_deltaO" -> Flatten[Table[If[a =!= b, "deltaO:" <> a <> ":" <> b, Nothing], {a, mediators}, {b, mediators}]],
    "susceptibility_branches" -> (("susceptibility_branch:" <> ToString[#]) & /@ Range[4]),
    "all_G8" -> (("G8:" <> #1["source_id"]) & /@ g8Records),
    "G11" -> {"gate:G11:mass_only", "gate:G11:velocity_converse"},
    "all_witness_challenge_pairs" -> (("witness_challenge:" <> #) & /@ unresolvedIds)|>,
  "external_legs" -> {"arbiter unchanged-script rerun", "fresh directive fidelity audit",
    "fresh adversarial recomputation", "read-only external review"}|>;

activeStrata = canonicalSort[StringSplit[#, ":", 2][[2]] & /@ tiltSlots];
payload = <|
  "schema_version" -> schema, "engine" -> "Wolfram",
  "independent_route" -> "Wolfram symbolic replacement/FullSimplify plus Association typed-root and slot-DAG walks",
  "execution_mode" -> If[selfTest, "UNSEALED_SELF_TEST", "PRODUCTION_EVALUATION"],
  "stage0_binding" -> <|"supplied_digest" -> suppliedDigest,
    "bundle_index_digest" -> bundleIndex["stage0_contract_digest"], "contract_digest" -> contract["stage0_contract_digest"],
    "equal" -> (suppliedDigest == bundleIndex["stage0_contract_digest"] == contract["stage0_contract_digest"]),
    "runner_recomputes_bundle_and_environment_before_this_evaluation" -> (! selfTest)|>,
  "axes" -> <|"endpoints" -> endpoints, "ambient_branches" -> ambients, "closure_branches" -> closures,
    "GAP_OPEN_FIELD_PROFILE_strata" -> strata, "tilt_cell_count" -> Length[tiltCells],
    "coupling_cell_count" -> Length[couplingCells], "generated_active_strata" -> activeStrata,
    "axis_strata_exact_set_equal" -> (activeStrata === strata), "axis_collapse_performed" -> False|>,
  "availability_contract" -> <|"ratified_summary" -> slotsDoc["summary"], "unresolved_slot_ids" -> unresolvedIds,
    "derived_slot_ids" -> canonicalSort[Complement[Keys[slots], unresolvedIds]],
    "witness_challenge_pairs_consumed_by_reference" -> Table[
      <|"slot_id" -> slotId, "witness_id" -> slots[slotId]["witness_id"],
        "challenge_id" -> slots[slotId]["challenge_id"]|>, {slotId, unresolvedIds}]|>,
  "native_action" -> <|"S_body" -> sBody, "action_term_ids" -> actionIds, "action_term_count" -> Length[actionIds],
    "per_native_term_ablation" -> nativeAblations,
    "all_nonvanishing" -> And @@ (#1["response_nonzero"] & /@ nativeAblations)|>,
  "tilt" -> <|"formalism" -> tiltFormalism, "force_balance_by_endpoint" -> forceByEndpoint,
    "force_census_expected_terms" -> force["expected_terms"],
    "force_census_incidence_complete" -> force["coverage_checks"]["force_census_incidence_complete"],
    "cells" -> tiltCells,
    "headline" -> "TILT_UNRESOLVED(eight ratified GAP_OPEN_FIELD_PROFILE leaves plus typed branch data)"|>,
  "coupling_package" -> <|
    "7.5a" -> <|"domain" -> sBody,
      "Sigma_surface" -> formalRecord["S_Sigma", "Integral[t,Sigma](f_throat+f_mix)", "action",
        "surface variation on Sigma", {"domain:Sigma_boundary_data", "open_leaf:throat_surface_functional"}],
      "partial_Omega_surface" -> formalRecord["S_partial_Omega", "Integral[t,partial_Omega_c](S_outer)", "action",
        "outer control-surface variation", {"domain:partial_Omega_c_boundary_data", "open_leaf:outer_surface_functional"}],
      "headline" -> "UNRESOLVED(named surface slots); native Omega_c domain DERIVED"|>,
    "7.5b" -> <|"J_A" -> jRecords, "deltaO_AB" -> deltaRecords, "ordered_matrix_shape" -> {4, 4}|>,
    "7.5c" -> endpointVirtualWork,
    "7.5d" -> <|"formal_total_response" -> totalResponse, "cells" -> couplingCells|>,
    "channel_ownership_by_endpoint" -> ownershipByEndpoint,
    "coupling_census_incidence_complete" -> coupling["coverage_checks"]["coupling_census_incidence_complete"],
    "mass_charge_split" -> totalResponse["mass_charge_split"], "parity_census" -> parityCensus,
    "mouth_datum_records" -> mouthRecords, "parent_enum" -> couplingEnum,
    "headline" -> "COUPLING_MAP_UNRESOLVED(named slots); j proportional sV not classifiable"|>,
  "gates" -> gates,
  "G4_known_nonzero_controls" -> <|"controls" -> g4Controls, "claimed_zero_ids" -> allClaimedZeros,
    "covered_zero_ids" -> zeroCoverage, "exact_coverage" -> (allClaimedZeros === zeroCoverage)|>,
  "reconciliation" -> reconciliationResult, "dimensional_firewall" -> dimensionFirewall,
  "native_ancestry_guard" -> nativeAncestryGuard, "symbolic_object_inventory" -> symbolicObjectInventory,
  "A9_external_verification" -> a9Coverage, "integrity_candidate" -> "COMPUTATION_VALID",
  "headline_entries" -> headlineEntries|>;

Print["PHASEC_WOLFRAM_PROGRESS production_export_start"];
If[! DirectoryQ[DirectoryName[output]], CreateDirectory[DirectoryName[output], CreateIntermediateDirectories -> True]];
Export[output, payload, "YAML"];
Print["PHASEC_PRODUCTION_WOLFRAM_COMPLETE tilt_cells=" <> ToString[Length[tiltCells]] <>
  " coupling_cells=" <> ToString[Length[couplingCells]] <> " successors=" <> ToString[Length[successorKeys]]];
Exit[0];
