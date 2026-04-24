Get[FileNameJoin[{DirectoryName[$InputFileName], "pde_v2_mathematica_common.wl"}]];

banner["Stage V2-22B solver-handoff validator Mathematica audit"];

schemaPath = fixturePath["stage_v2_22b_solver_output_schema.json"];
validPath = fixturePath["stage_v2_22b_sample_solver_output_valid.json"];
invalidPath = fixturePath["stage_v2_22b_sample_solver_output_invalid_hardcap.json"];

checkTrue["schema fixture hash", fileSHA256[schemaPath] == "d3714c84920027a1432b6b510c434fca801020e99cda04fcc0e8b6c2dd3849ca"];
checkTrue["valid fixture hash", fileSHA256[validPath] == "c15682194fafb3d129027780b67ca4ba62c16bd03bfb7e54283044799a353fad"];
checkTrue["invalid fixture hash", fileSHA256[invalidPath] == "2c194f20f56f1cbfddcf3190079dc81025f966a5f99c20d82d32ffa3bacba121"];

schema = jsonImport[schemaPath];
valid = jsonImport[validPath];
invalid = jsonImport[invalidPath];

subbanner["Schema and valid packet"];
checkTrue["schema declares required fields", Length[schema["required_top_level_keys"]] >= 8];
checkTrue["schema requires solver metadata", MemberQ[schema["required_top_level_keys"], "solver_metadata"]];
checkTrue["valid schema", valid["schema"] == "stage_v2_22b_solver_handoff/v1"];
checkTrue["valid pre-target freeze", TrueQ[valid["freeze"]["pre_target_freeze"]]];
checkTrue["valid target blind", TrueQ[valid["freeze"]["target_blind"]]];
checkTrue["valid open boundary", valid["geometry"]["boundary_class"] == "open_impedance" && valid["geometry"]["R_exit"] > 0];
grid = valid["grid"]["points"];
checkTrue["valid monotone grid endpoints", First[grid] == 0 && Last[grid] > First[grid]];
checkTrue["valid profile lengths match grid", And @@ ((Length[#] == Length[grid]) & /@ Values[valid["profiles"]])];
metadata = valid["solver_metadata"];
checkTrue["valid metadata mesh length", metadata["mesh_points"] == Length[grid]];
checkTrue["valid metadata provenance strings", And @@ (StringLength[metadata[#]] > 0 & /@ {"exporter", "coefficient_family", "source_commit"})];
checkTrue["valid metadata residual finite", NumberQ[metadata["nonlinear_residual_norm"]] && metadata["nonlinear_residual_norm"] >= 0];

port = valid["mixed_ports"][[1]];
delta = port["Omega_U"]^2 port["Omega_W"]^2 - port["lambda_R"]^2;
checkTrue["valid mixed block Delta positive", delta > 0, delta];

subbanner["Invalid hard-cap negative control"];
checkTrue["invalid schema", invalid["schema"] == "stage_v2_22b_solver_handoff/v1"];
checkTrue["invalid hard cap is present", invalid["geometry"]["boundary_class"] == "hard_cap"];
checkTrue["invalid exit radius rejected", invalid["geometry"]["R_exit"] == 0];

Print["All Stage V2-22B Mathematica fixture checks passed."];
