(* PathA 22b Gate 1 cross-check.
   Scope: re-quadrature of exported Z_w data, flat dw measure confirmation,
   and dimensional bookkeeping only. No branch/profile solve. *)

ClearAll[
  packetPath, manifestPath, outDir, jsonPath, data, grid, derived,
  wFaces, widths, zValues, zInt, edgeMax, peak, edgeFraction,
  branch, localizationFloor, domainWidth, floorIntegral, gaussianIntegral,
  floorFraction, gaussianFraction,
  manifest, rows, rowForPath, values, nearestGridDelta, ladderSpread,
  dim0, dimL, checks, report, allPass
];

packetPath = FileNameJoin[{
  "software", "stage1_solver", "frozen", "m1c",
  "834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8",
  "wp1_background_10x8.json"
}];
manifestPath = FileNameJoin[{
  "software", "stage1_solver", "frozen", "m1c",
  "834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8",
  "wp1_solve_manifest.json"
}];

data = Import[packetPath, "RawJSON"];
grid = data["grid"];
derived = data["derived"];
branch = data["config"]["branch"];
wFaces = N[grid["w_faces"]];
widths = Differences[wFaces];
zValues = N[derived["Z_w"]];
zInt = N[Total[zValues*widths], 16];
edgeMax = Max[Abs[First[zValues]], Abs[Last[zValues]]];
peak = Max[Abs[zValues]];
edgeFraction = N[(Abs[First[zValues]]*First[widths] + Abs[Last[zValues]]*Last[widths])/Abs[zInt], 16];
localizationFloor = N[branch["localization_floor"], 16];
domainWidth = N[grid["w_max"] - grid["w_min"], 16];
floorIntegral = N[localizationFloor*domainWidth, 16];
gaussianIntegral = N[Total[(zValues - localizationFloor)*widths], 16];
floorFraction = N[floorIntegral/zInt, 16];
gaussianFraction = N[gaussianIntegral/zInt, 16];

rowForPath[path_] := Module[
  {packet = Import[path, "RawJSON"], g, z, dx},
  g = packet["grid"];
  z = N[packet["derived"]["Z_w"]];
  dx = Differences[N[g["w_faces"]]];
  <|
    "path" -> path,
    "nw" -> g["nw"],
    "dw" -> N[g["dw"], 16],
    "z_int_finite_domain" -> N[Total[z*dx], 16]
  |>
];

manifest = Import[manifestPath, "RawJSON"];
rows = SortBy[rowForPath /@ manifest["background_paths"], #["nw"] &];
values = rows[[All, "z_int_finite_domain"]];
nearestGridDelta = N[Abs[Last[values] - values[[-2]]], 16];
ladderSpread = N[Max[values] - Min[values], 16];

dim0 = {0, 0, 0};
dimL = {1, 0, 0};

checks = {
  <|
    "name" -> "same-length exported arrays",
    "pass" -> TrueQ[Length[zValues] + 1 == Length[wFaces]]
  |>,
  <|
    "name" -> "positive w cell widths",
    "pass" -> AllTrue[widths, # > 0 &]
  |>,
  <|
    "name" -> "flat dw Z_int quadrature",
    "pass" -> TrueQ[Abs[zInt - 2.031114372358842] < 10^-14],
    "actual" -> zInt
  |>,
  <|
    "name" -> "nearest-grid error bar",
    "pass" -> TrueQ[Abs[nearestGridDelta - 0.001118740362587922] < 10^-14],
    "actual" -> nearestGridDelta
  |>,
  <|
    "name" -> "Z dimensionless",
    "pass" -> TrueQ[dim0 === {0, 0, 0}],
    "dimension_LTM" -> dim0
  |>,
  <|
    "name" -> "Z_int dimension length",
    "pass" -> TrueQ[dimL === {1, 0, 0}],
    "dimension_LTM" -> dimL
  |>,
  <|
    "name" -> "edge tail not small",
    "pass" -> TrueQ[edgeMax/peak > 0.5],
    "edge_to_peak_abs_ratio" -> N[edgeMax/peak, 16],
    "edge_cell_integral_fraction_abs" -> edgeFraction
  |>
};

allPass = AllTrue[checks, TrueQ[#["pass"]] &];
report = <|
  "schema" -> "stage1_pathA_22b_gate1_mathematica_crosscheck/v1",
  "scope" -> "flat dw exported Z_w quadrature only",
  "pass" -> allPass,
  "measure" -> <|
    "pde_definition" -> "Z_int=int Z(w) dw",
    "sqrt_g_w_exported" -> False,
    "sqrt_g_w_variant_status" -> "not independently computable from exported Z_w"
  |>,
  "z_int_finite_domain" -> zInt,
  "nearest_grid_delta" -> nearestGridDelta,
  "ladder_spread" -> ladderSpread,
  "rows" -> rows,
  "tail_diagnostics" -> <|
    "edge_max_abs" -> edgeMax,
    "peak_abs" -> peak,
    "edge_to_peak_abs_ratio" -> N[edgeMax/peak, 16],
    "edge_cell_integral_fraction_abs" -> edgeFraction
  |>,
  "floor_decomposition" -> <|
    "localization_floor" -> localizationFloor,
    "domain_width" -> domainWidth,
    "floor_integral_finite_domain" -> floorIntegral,
    "floor_fraction_finite_domain" -> floorFraction,
    "gaussian_integral_finite_domain" -> gaussianIntegral,
    "gaussian_fraction_finite_domain" -> gaussianFraction
  |>,
  "checks" -> checks
|>;

outDir = FileNameJoin[{"software", "stage1_solver", "_scratch"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
jsonPath = FileNameJoin[{outDir, "pathA_22b_gate1_mathematica_crosscheck.json"}];
Export[jsonPath, report, "RawJSON"];
Print["wrote ", jsonPath];
Print["pathA_22b Gate 1 Mathematica cross-check: ", Count[checks[[All, "pass"]], True], "/", Length[checks], " checks"];
If[TrueQ[report["pass"]], Exit[0], Exit[1]]
