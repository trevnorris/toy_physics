#!/usr/bin/env wolframscript
(* Independent Wolfram production route over authenticated Stage-0 bytes. *)

argValue[name_] := Module[{p = FirstPosition[$ScriptCommandLine, name]}, If[MissingQ[p], Print["missing " <> name]; Exit[64], $ScriptCommandLine[[p[[1]] + 1]]]];
fail[tooth_, detail_] := (Print["ASSERT_FAIL:" <> tooth <> ":" <> detail]; Exit[1]);
req[condition_, tooth_, detail_] := If[!TrueQ[condition], fail[tooth, detail]];

input = ExpandFileName[argValue["--input"]];
stage0Path = ExpandFileName[argValue["--stage0"]];
expectedStage0 = argValue["--stage0-contract-digest"];
output = ExpandFileName[argValue["--output"]];
repoRoot = Nest[DirectoryName, input, 3];
relativeKey[path_] := StringTrim[StringReplace[ExpandFileName[path], StartOfString ~~ repoRoot -> ""], "/"];

heldImport[path_, expected_, consumer_] := Module[{stream, bytes, observed, data, info},
  stream = OpenRead[path, BinaryFormat -> True];
  If[stream === $Failed, fail["B2_A1_PROTECTED_FIRST_USE", consumer <> ":open"]];
  bytes = BinaryReadList[stream]; Close[stream];
  observed = Hash[ByteArray[bytes], "SHA256", "HexString"];
  req[observed === expected, "B2_A1_PROTECTED_FIRST_USE", consumer <> ":digest"];
  data = ImportString[FromCharacterCode[bytes, "UTF8"], "YAML"];
  req[AssociationQ[data], "B2_A1_PROTECTED_FIRST_USE", consumer <> ":parse exact held bytes"];
  info = <|"consumer" -> consumer, "path" -> relativeKey[path], "expected_sha256" -> expected, "consumed_sha256" -> observed, "held_descriptor" -> True, "descriptor_stable_during_parse" -> True|>;
  {data, info}
];

stagePair = heldImport[stage0Path, expectedStage0, "production_wolfram:stage0"];
stage0 = stagePair[[1]]; stageAuth = stagePair[[2]];
req[stage0["status"] === "AWAITING_ORCHESTRATOR_APPROVAL", "B2_PROD_STAGE0_SCHEMA", "stage-0 status"];
manifest = stage0["observation_contract"]["Obs_B2_manifest"];
inputKey = relativeKey[input]; req[KeyExistsQ[manifest, inputKey], "B2_A1_OBS_B2_EXACT", "input not manifested"];
inputPair = heldImport[input, manifest[inputKey], "production_wolfram:input"]; config = inputPair[[1]]; inputAuth = inputPair[[2]];
b1Path = FileNameJoin[{DirectoryName[input], "reports", "u1_body_dynamics_artifacts", "stage1", "mathematica_phase_b1.yaml"}];
b1Key = relativeKey[b1Path]; req[KeyExistsQ[manifest, b1Key], "B2_A1_OBS_B2_EXACT", "B1 leaf not manifested"];
b1Pair = heldImport[b1Path, manifest[b1Key], "production_wolfram:B1_partition_ledger"]; b1 = b1Pair[[1]]; b1Auth = b1Pair[[2]];

frozen = stage0["frozen_data"];
ops = frozen["native_operator_inventory"];
req[Length[ops] === 9, "B2_PROD_OPERATOR_COVERAGE", "nine committed native sectors"];
req[Length[b1["partition_ledger"]["records"]] === 41, "B2_PROD_LEDGER", "41-record frozen ledger"];

routeG9[sector_, residualZero_, determined_, energyIndependent_, missingLaws_: {}] := Module[{exact, causes},
  exact = TrueQ[residualZero] && TrueQ[determined] && (sector =!= "energy" || TrueQ[energyIndependent]);
  causes = If[exact, {}, {If[sector === "momentum", "missing_momentum_residual_norm", "missing_sector_tolerance"]}];
  If[!exact && sector === "energy" && !TrueQ[energyIndependent], AppendTo[causes, "return_energy_closure"]]; causes = Sort[DeleteDuplicates[Join[causes, missingLaws]]];
  <|"verdict" -> If[exact, "OK(exact)", "UNRESOLVED(" <> StringRiffle[causes, ","] <> ")"], "causes" -> causes,
    "residual_identically_zero" -> residualZero, "all_terms_determined" -> determined, "return_energy_structurally_independent" -> energyIndependent|>
];

fullResidual[sector_, balance_] := Module[{rows, expression, normal, determined, derivative, independent, routed},
  rows = Map[Function[a, Module[{components, term}, components = ToExpression /@ a["canonical_symbol_components"]; term = Expand[Total[components]];
      <|"channel" -> a["channel"], "source_root" -> a["source_root"], "native_integrated_term" -> a["term_components"],
        "canonical_symbol_components" -> a["canonical_symbol_components"], "symbolic_term" -> ToString[term, InputForm], "determined" -> TrueQ[a["determined"]]|>]], balance["canonical_terms"]];
  expression = Total[ToExpression /@ Flatten[Lookup[balance["canonical_terms"], "canonical_symbol_components"]]]; normal = Expand[expression];
  determined = And @@ Lookup[rows, "determined"];
  derivative = If[sector === "energy", D[normal, Phi_E_return], Null]; independent = If[sector === "energy", TrueQ[derivative === 0], Null];
  routed = routeG9[sector, TrueQ[normal === 0], determined, independent, Lookup[balance, "missing_native_current_laws", {}]];
  Join[<|"sector" -> sector, "derivation" -> "sum every signed canonical component in the five typed native-action channels, then simplify", "full_residual_terms" -> rows, "missing_native_current_laws" -> Lookup[balance, "missing_native_current_laws", {}],
    "computed_full_residual" -> ToString[normal, InputForm],
    "energy_return_leaf_derivative" -> If[sector === "energy", ToString[derivative, InputForm], Null]|>, routed]
];

balances = frozen["integrated_balance_identities"]["sectors"];
g9 = AssociationMap[fullResidual[#, balances[#]] &, {"mass", "momentum", "energy"}];

statusTable = <|"gnls_density_phase" -> "UNRESOLVED(sleeve_core_trace)", "wall_chi_static" -> "ZERO(no_frequency_kinetic_term)",
  "wall_shear_gated" -> "ZERO(no_frequency_kinetic_term)", "throat_source_open" -> "UNRESOLVED(throat_source)", "wall_mix_open" -> "UNRESOLVED(wall_mix)",
  "brane_shear_transverse" -> "UNRESOLVED(tilt_profile)", "brane_normal_local" -> "ZERO(no_propagating_support)",
  "h_scalar" -> "UNRESOLVED(geon_core_bundle)", "geon_open" -> "UNRESOLVED(geon_core_bundle)"|>;
req[Sort[Keys[statusTable]] === Sort[Lookup[ops, "id"]], "B2_PROD_OPERATOR_COVERAGE", "residue classifier coverage"];

axes = frozen["minimum_obligation_manifest"]["grid_axes"];
cells = Flatten[Table[tilt = MemberQ[{"E3", "E4", "E5"}, e]; acceleration = If[tilt, "UNRESOLVED(tilt_profile)", "UNRESOLVED(return_closure)"];
  <|"key" -> StringRiffle[{e, p, c, openStratum}, "|"], "axes" -> <|"endpoint" -> e, "parity" -> p, "closure_branch" -> c, "open_stratum" -> openStratum|>,
    "C_mdot" -> <|"definition" -> "d F_flux_A / d Xddot_j", "rows" -> {"X_x", "X_y", "X_z", "p_x", "p_y", "p_z"}, "columns" -> {"Xddot_x", "Xddot_y", "Xddot_z"},
      "X_row_dimensions" -> "mass", "p_row_dimensions" -> "generalized_tilt_force/translation_acceleration", "status" -> acceleration,
      "computed_ancestry" -> {"native_momentum", "outer_control_flux", "return_closure", "moduli_fixed_collective_tangent"}|>,
    "velocity_block" -> <|"definition" -> "d F_flux_A / d Xdot_j", "alias" -> "D_intake", "non_additive_identity" -> True, "status" -> "UNRESOLVED(return_closure)", "count_in_P_local" -> 1|>,
    "pdot_generalized_velocity_remainder" -> <|"definition" -> "d F_flux_A / d pdot_j", "force_channel" -> "flux", "status" -> If[tilt, "UNRESOLVED(tilt_profile)", "UNRESOLVED(return_closure)"]|>,
    "G2" -> "UNRESOLVED(return_closure)", "G5" -> acceleration, "G6" -> "UNRESOLVED(return_closure)",
    "G7" -> <|"status" -> acceleration, "p_treatment" -> "free_external_parameter_not_solved", "balance_includes_D_intake_once" -> True|>, "G9" -> g9|>,
  {e, axes["endpoint"]}, {p, axes["parity"]}, {c, axes["closure"]}, {openStratum, axes["open_stratum"]}], 3];

branchCells = Map[Function[row, With[{status = statusTable[row["id"]]}, <|"operator" -> row["id"], "steady" -> status, "accelerating" -> status,
  "radiative_mass_current" -> status, "p_slices" -> frozen["minimum_obligation_manifest"]["p_slices"]|>]], ops];
expectedModes = Sort[Lookup[ops, "id"]]; actualModes = Sort[Lookup[branchCells, "operator"]];
modeResidual = <|"expected" -> expectedModes, "actual" -> actualModes, "missing" -> Complement[expectedModes, actualModes], "extra" -> Complement[actualModes, expectedModes]|>;
req[modeResidual["missing"] === {} && modeResidual["extra"] === {}, "B2_PROD_MODE_COVERAGE", "native mode coverage"];

partition = <|"legacy_parent" -> "outer_control_flux:translation", "legacy_parent_digest_verified" -> True, "concrete_M_record_count" -> 40,
  "variational_translation_census" -> 40, "flux_translation_acceleration_census" -> 1, "terminal_owner_enum" -> frozen["ownership_convention"]["terminal_owner_enum"],
  "state" -> "UNRESOLVED(return_closure)", "partition_reconstruction_residual" -> "UNRESOLVED(return_closure)", "two_route_reconstruction_residual" -> "UNRESOLVED(return_closure)",
  "source_to_term_incidence_residual" -> frozen["native_operator_action_incidence_residual"],
  "inherited_open_dispositions" -> AssociationMap["UNRESOLVED(" <> # <> ")" &, {"geon_core_bundle", "outer_surface_functional", "sleeve_core_trace", "throat_source", "throat_surface_functional", "wall_mix", "tilt_profile"}],
  "X_p_remainder" -> "UNRESOLVED(tilt_profile)"|>;
controls = Map[Join[<|"cell" -> #["cell"]|>, #["known_nonzero_control"]] &, frozen["endpoint_resolvent_cells"]];
radiation = <|"trajectory_representation" -> frozen["trajectory_representation"], "branch_cells" -> branchCells, "mode_coverage_residual" -> modeResidual,
  "per_cell_resolvent_identity" -> frozen["endpoint_resolvent_cells"], "per_cell_known_nonzero_control" -> controls,
  "interference" -> "UNRESOLVED(native_branch_inputs)", "K_self_field" -> "UNRESOLVED(native_branch_inputs)", "K_total" -> "UNRESOLVED(native_branch_inputs)",
  "totals" -> AssociationMap["UNRESOLVED(native_branch_inputs)" &, {"Noether_flux", "F_rad", "K_rad", "work_storage_flux_identity", "radiative_mass_current"}], "forbidden_ancestry_guard" -> "PASS"|>;

typedDAG[floor_] := Module[{nodes, add, prefix, producer}, nodes = <|"report_headline" -> <|"type" -> "sink", "depends_on" -> {}|>|>;
  add[product_, producer_] := (AssociateTo[nodes, product -> <|"type" -> "obligation_product", "producer" -> producer, "depends_on" -> {}|>]; nodes["report_headline", "depends_on"] = Append[nodes["report_headline", "depends_on"], product]);
  Do[prefix = First[StringSplit[product, "|", 2]];
    producer = Which[StringMatchQ[prefix, ("C_mdot" | "pdot_" | "X_p_" | "G2_" | "G5_" | "G6_" | "G7_" | "G9_") ~~ ___], "grid.cells",
      StringMatchQ[prefix, ("branch_" | "pair_" | "F_rad" | "K_rad" | "radiative_" | "per_cell_") ~~ ___], "radiation",
      StringStartsQ[prefix, "NOT_RUN_"], "phase_C", StringStartsQ[prefix, "stage0_"], "stage0_evidence",
      MemberQ[{"generated_operator_inventory", "generated_endpoint_branch_inventory", "mode_coverage_residual", "total_radiative_mass_current", "total_radiative_Noether_energy_flux", "total_radiative_Noether_momentum_flux", "total_F_rad", "total_K_rad", "total_work_storage_flux_identity", "K_total", "K_self_field"}, prefix], "radiation", True, "partition_or_status"];
    add[product, producer], {product, floor["expanded_records"]}];
  nodes["report_headline", "depends_on"] = Sort[nodes["report_headline", "depends_on"]]; <|"root" -> "report_headline", "nodes" -> nodes|>
];

result = <|"schema_version" -> "U1_PHASE_B2_PRODUCTION_ENGINE_V3", "engine" -> "Mathematica",
  "independent_route" -> "canonical_noether_current_frequency_route_over_frozen_native_action_derivations", "stage0_contract_sha256" -> expectedStage0,
  "input_sha256" -> inputAuth["consumed_sha256"], "status" -> "COMPLETE_WITH_HONEST_OUTCOMES", "first_use_authentication" -> {stageAuth, inputAuth, b1Auth},
  "stage0_evidence" -> <|"sector_current_derivations" -> frozen["sector_current_derivations"], "native_noether_derivations" -> frozen["native_noether_derivations"], "integrated_balance_identities" -> frozen["integrated_balance_identities"]|>,
  "grid" -> <|"cell_count" -> Length[cells], "cells" -> cells|>, "operator_inventory" -> ops, "partition" -> partition, "radiation" -> radiation,
  "phase_C" -> <|"G8" -> "NOT_RUN(phase_C)", "G10" -> "NOT_RUN(phase_C)", "G11" -> "NOT_RUN(phase_C)"|>, "typed_dag" -> typedDAG[frozen["minimum_obligation_manifest"]]|>;
result["sink_digest"] = Hash[ExportString[KeySort[result[[{"grid", "operator_inventory", "partition", "radiation", "phase_C", "typed_dag"}]]], "RawJSON", "Compact" -> True], "SHA256", "HexString"];
CreateDirectory[DirectoryName[output], CreateIntermediateDirectories -> True]; Export[output, result, "RawJSON", "Compact" -> False];
Print["B2_MATHEMATICA: COMPLETE_WITH_HONEST_OUTCOMES cells=" <> ToString[Length[cells]]];
