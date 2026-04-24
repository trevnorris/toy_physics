Get[FileNameJoin[{DirectoryName[$InputFileName], "pde_v2_mathematica_common.wl"}]];

banner["Stage V2-22A profile-to-coefficient adapter Mathematica audit"];

path = fixturePath["stage_v2_22a_profile_input_manifest.json"];
manifest = jsonImport[path];
checkTrue["fixture schema", manifest["schema"] == "stage_v2_22a_profile_adapter/v1"];
checkTrue["fixture hash", fileSHA256[path] == "ecbd9f008bde0621cad62361ce98708041c4c98ebfd8b4ecdd8b04b579d0aacd"];
checkTrue["three profile branches", Length[manifest["branches"]] == 3];

names = manifest["branches"][[All, "name"]];
checkTrue["analytic fixture present", MemberQ[names, "analytic_isotropic_DN_profile_demo"]];
checkTrue["sampled fixture present", MemberQ[names, "sampled_isotropic_DN_profile_demo"]];
checkTrue["weak-axisymmetric fixture present", MemberQ[names, "weak_axisymmetric_profile_slope_demo"]];

sampled = SelectFirst[manifest["branches"], #["name"] == "sampled_isotropic_DN_profile_demo" &];
profiles = sampled["profiles"];
sampleLengths = Length /@ (profiles[#]["samples"] & /@ Select[Keys[profiles], profiles[#]["kind"] == "sampled" &]);
checkTrue["sampled profiles have common length", Length[Union[sampleLengths]] == 1, sampleLengths];
checkTrue["sampled profiles have dense grid", First[sampleLengths] >= 801, First[sampleLengths]];
checkTrue["sampled boundary class open", sampled["geometry"]["boundary_class"] == "open_impedance"];
checkTrue["sampled exit positive", sampled["geometry"]["R_exit"] > 0];

weak = SelectFirst[manifest["branches"], #["name"] == "weak_axisymmetric_profile_slope_demo" &];
checkTrue["weak branch has weak-axisymmetric data", KeyExistsQ[weak, "weak_axisymmetric"]];

Print["All Stage V2-22A Mathematica fixture checks passed."];
