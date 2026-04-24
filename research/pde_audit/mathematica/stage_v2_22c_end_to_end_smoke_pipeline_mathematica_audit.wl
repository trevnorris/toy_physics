Get[FileNameJoin[{DirectoryName[$InputFileName], "pde_v2_mathematica_common.wl"}]];

banner["Stage V2-22C end-to-end smoke-pipeline Mathematica audit"];

validPath = fixturePath["stage_v2_22c_valid_solver_packet.json"];
invalidPath = fixturePath["stage_v2_22c_invalid_hardcap_solver_packet.json"];
valid = jsonImport[validPath];
invalid = jsonImport[invalidPath];

subbanner["Fixture controls"];
checkTrue["valid packet hash", fileSHA256[validPath] == "c15682194fafb3d129027780b67ca4ba62c16bd03bfb7e54283044799a353fad"];
checkTrue["invalid packet hash", fileSHA256[invalidPath] == "2c194f20f56f1cbfddcf3190079dc81025f966a5f99c20d82d32ffa3bacba121"];
checkTrue["valid open branch", valid["geometry"]["boundary_class"] == "open_impedance" && valid["geometry"]["R_exit"] > 0];
checkTrue["invalid hard-cap branch", invalid["geometry"]["boundary_class"] == "hard_cap" && invalid["geometry"]["R_exit"] == 0];
checkTrue["valid target blind", TrueQ[valid["freeze"]["target_blind"]] && TrueQ[valid["freeze"]["pre_target_freeze"]]];
checkTrue["valid solver metadata present", KeyExistsQ[valid, "solver_metadata"] && valid["solver_metadata"]["mesh_points"] == Length[valid["grid"]["points"]]];

subbanner["Orchestration formula audit"];
Clear[x20, x21, x22, eps, x0, x1, D0, D2, D4, N0, N2, N4, G, cs, a, c, mhat, Sport];
$Assumptions = D0 != 0 && Element[{x20, x21, x22, eps, x0, x1, D0, D2, D4, N0, N2, N4}, Reals] &&
  G > 0 && cs > 0 && a > 0 && c > 0 && mhat > 0 && Sport > 0;
xb = weightedBar[x20, x21, x22]; ax = groupA[x20, x21, x22]; bx = groupB[x20, x21, x22];
checkEqual["group inverse 20", xb + 4 ax, x20];
checkEqual["group inverse 21", xb - ax + bx, x21];
checkEqual["group inverse 22", xb - ax - bx, x22];
y20 = x0 + eps x1; y21 = x0 + eps x1/2; y22 = x0 - eps x1;
checkZero["axisymmetric b equals 3a", groupB[y20, y21, y22] - 3 groupA[y20, y21, y22]];
N2const = 2 D2 N0/D0;
N4const = N0 (D2^2 + 2 D0 D4)/D0^2;
checkZero["constant P2", prefactorP2[D0, D2, N0, N2const]];
checkZero["constant P4", prefactorP4[D0, D2, D4, N0, N2const, N4const]];
Ptarget = 54 G cs^5/(5 a^5 c^5);
gammaEff = mhat^2 Sport Ptarget a^5/(27 cs^5);
gammaGR = Sport mhat^2 2 G/(5 c^5);
checkZero["target gamma equivalence", gammaEff - gammaGR];

Print["All Stage V2-22C Mathematica checks passed."];
