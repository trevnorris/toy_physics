#!/usr/bin/env wolframscript

(* MT15-05 / Spike-3: direct-derived V2-22B packet producer.

   This script does not derive new physics. It imports the M1b derived
   BdG/wall sectors and the Spike-2 derived Maxwell transfer, then emits a
   V2-22B handoff packet that uses the direct-coefficient V2-21 lane.

   The packet must not contain target residual outputs. A separate intentionally
   invalid packet with an injected R_norm key is emitted only as a validator
   teeth test.
*)

$HistoryLength = 0;

ClearAll["Global`*"];

scriptPath =
  If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
    $InputFileName,
    FileNameJoin[{Directory[], "software/stage1_solver/mathematica/mt15_05_spike3_clean_rnorm.wls"}]
  ];
scriptDir = DirectoryName[scriptPath];
repoRoot = AbsoluteFileName[FileNameJoin[{scriptDir, "..", "..", ".."}]];
runDir = FileNameJoin[{scriptDir, "runs", "mt15_05_spike3_clean_rnorm"}];
If[! DirectoryQ[runDir], CreateDirectory[runDir, CreateIntermediateDirectories -> True]];

m1bRunDir = FileNameJoin[{scriptDir, "runs", "mt15_02_bdg_wall_derivation"}];
spike2RunDir = FileNameJoin[{scriptDir, "runs", "mt15_04_spike2_transfer_n0"}];

m1bPacketPath = FileNameJoin[{m1bRunDir, "mt15_02_v2_22b_bdg_wall_packet.json"}];
m1bDiagnosticsPath = FileNameJoin[{m1bRunDir, "mt15_02_derivation_diagnostics.json"}];
spike2DiagnosticsPath = FileNameJoin[{spike2RunDir, "mt15_04_spike2_transfer_n0_diagnostics.json"}];

packetPath = FileNameJoin[{runDir, "mt15_05_v2_22b_direct_derived_packet.json"}];
badPacketPath = FileNameJoin[{runDir, "mt15_05_v2_22b_direct_derived_forbidden_R_norm_packet.json"}];
diagnosticsPath = FileNameJoin[{runDir, "mt15_05_direct_derived_diagnostics.json"}];
reportPath = FileNameJoin[{scriptDir, "mt15_05_spike3_clean_rnorm_report.md"}];

ClearAll[
  toJsonNumber, sourceCommitString, sha256File, requireFile, importJson,
  formatNumber, markdownTable, boolQ
];

toJsonNumber[x_] := Module[{n = Quiet[Check[N[x], $Failed]]},
  If[n === $Failed || ! NumericQ[n],
    Print["ERROR: nonnumeric value reached JSON export: ", InputForm[x]];
    Exit[4]
  ];
  N[Re[n], 17]
];

sourceCommitString[] := Module[{proc, out},
  proc = Quiet[RunProcess[{"git", "-C", repoRoot, "rev-parse", "--short", "HEAD"}, "StandardOutput"]];
  out = If[StringQ[proc], StringTrim[proc], ""];
  If[StringLength[out] > 0, out, "unavailable"]
];

sha256File[path_String] := Module[{text, hash},
  If[! FileExistsQ[path], Return["missing"]];
  text = Quiet[Check[Import[path, "Text"], $Failed]];
  If[text === $Failed || ! StringQ[text], Return["unavailable"]];
  hash = Hash[text, "SHA256"];
  IntegerString[hash, 16, 64]
];

requireFile[path_String, label_String] := If[! FileExistsQ[path],
  Print["ERROR: required ", label, " file is missing: ", path];
  Exit[2]
];

importJson[path_String] := Module[{data = Quiet[Check[Import[path, "RawJSON"], $Failed]]},
  If[data === $Failed || ! AssociationQ[data],
    Print["ERROR: failed to import JSON association: ", path];
    Exit[3]
  ];
  data
];

formatNumber[x_] := ToString[FortranForm[N[x, 17]]];

markdownTable[headers_List, rows_List] := StringRiffle[
  Join[
    {StringJoin["| ", StringRiffle[headers, " | "], " |"]},
    {StringJoin["| ", StringRiffle[ConstantArray["---", Length[headers]], " | "], " |"]},
    (StringJoin["| ", StringRiffle[ToString /@ #, " | "], " |"] & /@ rows)
  ],
  "\n"
];

boolQ[x_] := TrueQ[x];

Scan[requireFile @@ # &, {
  {m1bPacketPath, "M1b V2-22B packet"},
  {m1bDiagnosticsPath, "M1b diagnostics"},
  {spike2DiagnosticsPath, "Spike-2 transfer diagnostics"}
}];

m1bPacket = importJson[m1bPacketPath];
m1bDiagnostics = importJson[m1bDiagnosticsPath];
spike2Diagnostics = importJson[spike2DiagnosticsPath];

If[m1bPacket["schema"] =!= "stage_v2_22b_solver_handoff/v1",
  Print["ERROR: unexpected M1b packet schema: ", m1bPacket["schema"]];
  Exit[5]
];

wallK = toJsonNumber[m1bPacket["wall"]["K"]];
wallM = toJsonNumber[m1bPacket["wall"]["M"]];
diagK = toJsonNumber[m1bDiagnostics["wall_weak_form"]["active_wall_K_scalar"]];
diagM = toJsonNumber[m1bDiagnostics["wall_weak_form"]["active_wall_M_scalar"]];
If[Abs[wallK - diagK] > 10^-10 || Abs[wallM - diagM] > 10^-10,
  Print["ERROR: M1b packet wall K/M does not match diagnostics."];
  Exit[6]
];

b0 = toJsonNumber[m1bDiagnostics["bdg_moments"]["B0"]];
b2 = toJsonNumber[m1bDiagnostics["bdg_moments"]["B2"]];
b4 = toJsonNumber[m1bDiagnostics["bdg_moments"]["B4"]];

z0 = toJsonNumber[spike2Diagnostics["coefficients"]["Z0"]];
z2 = toJsonNumber[spike2Diagnostics["coefficients"]["Z2"]];
z4 = toJsonNumber[spike2Diagnostics["coefficients"]["Z4"]];
n0 = toJsonNumber[spike2Diagnostics["coefficients"]["N0"]];
n2 = toJsonNumber[spike2Diagnostics["coefficients"]["N2"]];
n4 = toJsonNumber[spike2Diagnostics["coefficients"]["N4"]];
gammaPort = toJsonNumber[spike2Diagnostics["open_dtn"]["Gamma_port"]];

passChecks = spike2Diagnostics["validation"]["pass_checks"];
pureGauge = spike2Diagnostics["validation"]["pure_gauge"];
basisInvariance = spike2Diagnostics["validation"]["basis_invariance"];
v209Regression = spike2Diagnostics["validation"]["v2_09_regression"];
hardCapContrast = spike2Diagnostics["validation"]["hard_cap_contrast"];
meshRows = spike2Diagnostics["mesh_convergence"]["rows"];
maxGreenResidual = Max[toJsonNumber /@ (Lookup[#, "maxGreenResidual"] & /@ meshRows)];

operatorGaugeResidualMetrics = Join[
  passChecks,
  <|
    "max_green_residual" -> toJsonNumber[maxGreenResidual],
    "pure_gauge_physical_field_norm" -> toJsonNumber[pureGauge["physical_field_norm"]],
    "pure_gauge_outgoing_flux_abs" -> toJsonNumber[Abs[pureGauge["outgoing_flux"]]],
    "basis_max_relative_Z_delta" -> toJsonNumber[basisInvariance["max_relative_Z_delta"]],
    "basis_max_relative_N_delta" -> toJsonNumber[basisInvariance["max_relative_N_delta"]],
    "v2_09_relative_Z_error" -> toJsonNumber[v209Regression["relative_Z_error"]],
    "v2_09_relative_N0_error" -> toJsonNumber[v209Regression["relative_N0_error"]],
    "hard_cap_to_open_boundary_norm_ratio" -> toJsonNumber[hardCapContrast["hard_cap_to_open_boundary_norm_ratio"]]
  |>
];

sourceHashes = <|
  "m1b_packet_sha256" -> sha256File[m1bPacketPath],
  "m1b_diagnostics_sha256" -> sha256File[m1bDiagnosticsPath],
  "spike2_diagnostics_sha256" -> sha256File[spike2DiagnosticsPath],
  "m1b_script_sha256" -> sha256File[FileNameJoin[{scriptDir, "mt15_02_bdg_wall_derivation.wls"}]],
  "spike2_script_sha256" -> sha256File[FileNameJoin[{scriptDir, "mt15_04_spike2_transfer_n0.wls"}]],
  "spike3_script_sha256" -> sha256File[scriptPath]
|>;

sourcePaths = <|
  "m1b_packet" -> m1bPacketPath,
  "m1b_diagnostics" -> m1bDiagnosticsPath,
  "spike2_diagnostics" -> spike2DiagnosticsPath
|>;

derivedBdgWall = <|
  "status" -> "derived_bdg_wall_m1b",
  "coefficients" -> <|
    "K" -> wallK,
    "M" -> wallM,
    "B0" -> b0,
    "B2" -> b2,
    "B4" -> b4
  |>,
  "machinery_label" -> m1bDiagnostics["integrity_notes"]["machinery_not_physics"],
  "source_paths" -> sourcePaths,
  "source_hashes" -> KeyTake[sourceHashes, {"m1b_packet_sha256", "m1b_diagnostics_sha256", "m1b_script_sha256"}]
|>;

derivedMaxwellTransfer = <|
  "status" -> "derived_green_function_transfer",
  "gauge_convention" -> "localized_H_equals_Z_G_equals_0_quotient_with_m1b_zero_spatial_background",
  "flux_normalization" -> <|
    "Gamma_port" -> gammaPort,
    "formula" -> spike2Diagnostics["open_dtn"]["Gamma_port_formula"],
    "legacy_v2_09_reduced_formula" -> spike2Diagnostics["open_dtn"]["legacy_v2_09_reduced_formula"],
    "legacy_factor_gap" -> toJsonNumber[spike2Diagnostics["open_dtn"]["legacy_factor_gap"]],
    "derived_not_fitted_to_chi_Q" -> boolQ[spike2Diagnostics["open_dtn"]["derived_not_fitted_to_chi_Q"]]
  |>,
  "coefficients" -> <|
    "Z0" -> z0,
    "Z2" -> z2,
    "Z4" -> z4,
    "N0" -> n0,
    "N2" -> n2,
    "N4" -> n4
  |>,
  "operator_gauge_residual_metrics" -> operatorGaugeResidualMetrics,
  "source_paths" -> sourcePaths,
  "source_hashes" -> KeyTake[sourceHashes, {"spike2_diagnostics_sha256", "spike2_script_sha256"}],
  "machinery_label" -> spike2Diagnostics["machinery_not_physics"]
|>;

directPacket = Join[
  KeyDrop[m1bPacket, {"mixed_ports"}],
  <|
    "branch_id" -> "mt15_spike3_clean_direct_derived_coefficients",
    "freeze" -> Join[
      m1bPacket["freeze"],
      <|
        "parent_action" -> "GNLS M1b derived BdG/wall sectors plus Spike-2 derived Maxwell Green-function transfer",
        "maxwell_transfer_status" -> "derived_green_function_transfer",
        "mixed_ports_status" -> "absent_for_direct_derived_N0_path"
      |>
    ],
    "wall" -> Join[
      m1bPacket["wall"],
      <|"derivation_status" -> "derived_weak_form_galerkin_m1b"|>
    ],
    "derived_bdg_wall_coefficients" -> derivedBdgWall,
    "derived_maxwell_transfer" -> derivedMaxwellTransfer,
    "direct_derived_extension" -> <|
      "status" -> "v2_22b_direct_derived_coefficient_path",
      "no_posited_ports_in_N0_path" -> True,
      "v2_21_lane" -> "direct_coefficients",
      "target_blind_constants" -> True
    |>,
    "solver_metadata" -> Join[
      KeyDrop[m1bPacket["solver_metadata"], {"mixed_ports_status"}],
      <|
        "exporter" -> "mt15_05_spike3_clean_rnorm.wls",
        "coefficient_family" -> "m1b_bdg_wall_plus_spike2_direct_derived_maxwell_transfer_v1",
        "source_commit" -> sourceCommitString[],
        "export_timestamp" -> DateString[Now, {"ISODate", "T", "Time", " ", "TimeZoneNameShort"}],
        "bdg_status" -> "derived_from_m1b",
        "wall_status" -> "derived_from_m1b",
        "maxwell_transfer_status" -> "derived_green_function_transfer_spike2",
        "mixed_ports_status" -> "absent_for_direct_derived_N0_path",
        "notes" -> "CLEAN structure packet: K/M/Bn are imported from M1b derived diagnostics and Zn/Nn from Spike-2 derived Green-function transfer. No posited mixed_ports key is present, and N0 is not re-posited. This remains a MACHINERY packet on the M1b smoke background, not a physical result until M1c."
      |>
    ]
  |>
];

badPacket = Join[
  directPacket,
  <|
    "derived_maxwell_transfer" -> Join[
      directPacket["derived_maxwell_transfer"],
      <|"R_norm" -> 0.0|>
    ]
  |>
];

directCoefficients = <|
  "K" -> wallK, "M" -> wallM,
  "B0" -> b0, "B2" -> b2, "B4" -> b4,
  "Z0" -> z0, "Z2" -> z2, "Z4" -> z4,
  "N0" -> n0, "N2" -> n2, "N4" -> n4
|>;

d0 = wallK - b0 - z0;
d2 = -(wallM + b2 + z2);
d4 = -(b4 + z4);
p0 = n0/d0;
p2 = (d0 n2 - 2 d2 n0)/d0^2;
p4 = (d0^2 n4 - 2 d0 (d2 n2 + d4 n0) + 3 d2^2 n0)/d0^3;
rPole = d0 (b4 + z4) - 3 (wallM + b2 + z2)^2;
constants = m1bPacket["constants"];
rNormMachinery = constants["mhat0"]^2 constants["S_port"] p0 -
  54 constants["G"] constants["c_s"]^5/(5 constants["a"]^5 constants["c"]^5);

diagnostics = <|
  "schema" -> "mt15_05_spike3_direct_derived_diagnostics/v1",
  "packet_path" -> packetPath,
  "bad_packet_path" -> badPacketPath,
  "source_paths" -> sourcePaths,
  "source_hashes" -> sourceHashes,
  "direct_coefficients" -> directCoefficients,
  "derived_bdg_wall_coefficients" -> derivedBdgWall,
  "derived_maxwell_transfer" -> derivedMaxwellTransfer,
  "machinery_label" -> "CLEAN structure: no posited ports in the N0 path. The number is MACHINERY on the M1b smoke background and is not physical until M1c.",
  "formula_preview_not_exported_in_packet" -> <|
    "D0" -> toJsonNumber[d0],
    "P2" -> toJsonNumber[p2],
    "P4" -> toJsonNumber[p4],
    "R_pole" -> toJsonNumber[rPole],
    "R_norm" -> toJsonNumber[rNormMachinery]
  |>
|>;

Export[packetPath, directPacket, "JSON"];
Export[badPacketPath, badPacket, "JSON"];
Export[diagnosticsPath, diagnostics, "JSON"];

coefficientRows = {
  {"K", formatNumber[wallK], "M1b mt15_02 wall weak form"},
  {"M", formatNumber[wallM], "M1b mt15_02 wall weak form"},
  {"B0", formatNumber[b0], "M1b mt15_02 BdG moments"},
  {"B2", formatNumber[b2], "M1b mt15_02 BdG moments"},
  {"B4", formatNumber[b4], "M1b mt15_02 BdG moments"},
  {"Z0", formatNumber[z0], "Spike-2 mt15_04 Green transfer"},
  {"Z2", formatNumber[z2], "Spike-2 mt15_04 Green transfer"},
  {"Z4", formatNumber[z4], "Spike-2 mt15_04 Green transfer"},
  {"N0", formatNumber[n0], "Spike-2 mt15_04 retarded transfer"},
  {"N2", formatNumber[n2], "Spike-2 mt15_04 retarded transfer"},
  {"N4", formatNumber[n4], "Spike-2 mt15_04 retarded transfer"}
};

report = StringRiffle[
  {
    "# MT15-05 Spike-3 Clean R_norm Direct-Derived Packet",
    "",
    "Generated from `software/stage1_solver/mathematica/mt15_05_spike3_clean_rnorm.wls`.",
    "",
    "**LOUD LABEL:** CLEAN structure packet: no posited `mixed_ports` in the `N0` path. The computed residual is a MACHINERY number on the M1b smoke background, not a physical falsification result until M1c.",
    "",
    "## Inputs",
    "",
    markdownTable[{"term", "value", "provenance"}, coefficientRows],
    "",
    "## Packet",
    "",
    "- Direct-derived packet: `" <> packetPath <> "`",
    "- Forbidden-field teeth packet: `" <> badPacketPath <> "`",
    "- Diagnostics: `" <> diagnosticsPath <> "`",
    "- `mixed_ports` key present in direct packet: `" <> ToString[KeyExistsQ[directPacket, "mixed_ports"]] <> "`",
    "- `derived_maxwell_transfer.status`: `" <> directPacket["derived_maxwell_transfer"]["status"] <> "`",
    "- `derived_bdg_wall_coefficients.status`: `" <> directPacket["derived_bdg_wall_coefficients"]["status"] <> "`",
    "- `Gamma_port`: `" <> formatNumber[gammaPort] <> "`",
    "",
    "## Direct Formula Preview",
    "",
    "These values are not exported in the V2-22B packet; the Python V2 chain remains the external residual judge.",
    "",
    "- `D0 = K - B0 - Z0 = " <> formatNumber[d0] <> "`",
    "- `R_norm = " <> formatNumber[rNormMachinery] <> "`",
    "- `R_pole = " <> formatNumber[rPole] <> "`",
    "- `P2 = " <> formatNumber[p2] <> "`",
    "- `P4 = " <> formatNumber[p4] <> "`",
    "",
    "## Term/Provenance Map",
    "",
    "- M1b (`mt15_02_bdg_wall_derivation.wls`): `K`, `M`, `B0`, `B2`, `B4`, wall/profile/BdG provenance.",
    "- Spike-2 (`mt15_04_spike2_transfer_n0.wls`): `Z0`, `Z2`, `Z4`, `N0`, `N2`, `N4`, `Gamma_port`, Green/gauge validation checks.",
    "- V2-22B extension: additive direct-derived branch. Legacy `mixed_ports` packets still follow the old path; this packet omits `mixed_ports` and generates V2-21 `direct_coefficients`.",
    "",
    "## V2 Chain",
    "",
    "Run the V2-22C chain on the direct-derived packet to produce the authoritative observable packet. This report records the packet-side machinery preview only until that run is completed."
  },
  "\n"
];

Export[reportPath, report, "Text"];

Print["MT15-05 Spike-3 direct-derived packet export"];
Print["packet_path: ", packetPath];
Print["bad_packet_path: ", badPacketPath];
Print["diagnostics_path: ", diagnosticsPath];
Print["report_path: ", reportPath];
Print["mixed_ports_present_in_direct_packet: ", KeyExistsQ[directPacket, "mixed_ports"]];
Print["D0: ", formatNumber[d0]];
Print["R_norm_preview: ", formatNumber[rNormMachinery]];
Print["R_pole_preview: ", formatNumber[rPole]];
Print["P2_preview: ", formatNumber[p2]];
Print["P4_preview: ", formatNumber[p4]];
