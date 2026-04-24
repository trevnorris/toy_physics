Get[FileNameJoin[{DirectoryName[$InputFileName], "pde_v2_mathematica_common.wl"}]];

banner["Stage V2-21 branch-extraction fixture Mathematica audit"];

path = fixturePath["stage_v2_21_sample_branch_manifest.json"];
manifest = jsonImport[path];
checkTrue["fixture schema", manifest["schema"] == "stage_v2_21_branch_extraction_fixture/v1"];
checkTrue["fixture branch count", Length[manifest["branches"]] == 2];
checkTrue["fixture hash", fileSHA256[path] == "5aff3aac38f31841c1d089c5383c2da809667cba2ac36ad39948fc95cc9ee7c2"];

laneCoeff[branch_, lane_String] := branch["lanes"][lane]["direct_coefficients"];
directResiduals[coeff_, const_] := Module[
  {K, M, B0, B2, B4, Z0, Z2, Z4, N0, N2, N4, D0, D2, D4, A, C, P0, P2, P4, Ptarget, Rnorm},
  K = coeff["K"]; M = coeff["M"]; B0 = coeff["B0"]; B2 = coeff["B2"]; B4 = coeff["B4"];
  Z0 = coeff["Z0"]; Z2 = coeff["Z2"]; Z4 = coeff["Z4"]; N0 = coeff["N0"]; N2 = coeff["N2"]; N4 = coeff["N4"];
  D0 = K - B0 - Z0; D2 = -(M + B2 + Z2); D4 = -(B4 + Z4);
  A = M + B2 + Z2; C = B4 + Z4;
  P0 = prefactorP0[D0, N0]; P2 = prefactorP2[D0, D2, N0, N2]; P4 = prefactorP4[D0, D2, D4, N0, N2, N4];
  Ptarget = 54 const["G"] const["c_s"]^5/(5 const["a"]^5 const["c"]^5);
  Rnorm = const["mhat0"]^2 const["S_port"] P0 - Ptarget;
  <|"D0" -> D0, "R_pole" -> D0 C - 3 A^2, "R_norm" -> Rnorm, "P2" -> P2, "P4" -> P4|>
];

cal = manifest["branches"][[1]];
const = cal["target"]["constants"];
res20 = directResiduals[laneCoeff[cal, "20"], const];
subbanner["Calibrated direct-coefficient branch"];
checkTrue["calibrated open boundary", cal["geometry"]["boundary_class"] == "open_impedance" && cal["geometry"]["R_exit"] > 0];
checkTrue["calibrated D0 positive", res20["D0"] > 0];
checkTrue["calibrated one-pole", nearZeroQ[res20["R_pole"]]];
checkTrue["calibrated normalization", nearZeroQ[res20["R_norm"]]];
checkTrue["calibrated P2", nearZeroQ[res20["P2"]]];
checkTrue["calibrated P4", nearZeroQ[res20["P4"]]];

demo = manifest["branches"][[2]];
subbanner["Primitive open-throat demo branch"];
checkTrue["demo open boundary", demo["geometry"]["boundary_class"] == "open_impedance" && demo["geometry"]["R_exit"] > 0];
checkTrue["demo is not direct-coefficient target fixture", !KeyExistsQ[demo["lanes"]["20"], "direct_coefficients"]];
checkTrue["demo has BdG modes", Length[demo["lanes"]["20"]["bdg_modes"]] == 2];
checkTrue["demo has mixed port", Length[demo["lanes"]["20"]["mixed_ports"]] == 1];

Print["All Stage V2-21 Mathematica fixture checks passed."];
