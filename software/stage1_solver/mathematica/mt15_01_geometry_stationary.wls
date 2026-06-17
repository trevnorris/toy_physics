#!/usr/bin/env wolframscript

(* M1a stationary open-throat geometry/profile producer.
   This script only exports a V2-22B JSON handoff packet.  The Python V2
   validator/extractor remains an external judge and is not imported here. *)

$HistoryLength = 0;

ClearAll["Global`*"];

scriptPath =
  If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
    $InputFileName,
    FileNameJoin[{Directory[], "software/stage1_solver/mathematica/mt15_01_geometry_stationary.wls"}]
  ];
scriptDir = DirectoryName[scriptPath];
repoRoot = AbsoluteFileName[FileNameJoin[{scriptDir, "..", "..", ".."}]];
runDir = FileNameJoin[{scriptDir, "runs", "mt15_01_geometry_stationary"}];
If[! DirectoryQ[runDir], CreateDirectory[runDir, CreateIntermediateDirectories -> True]];

packetPath = FileNameJoin[{runDir, "mt15_01_v2_22b_placeholder_packet.json"}];
badPacketPath = FileNameJoin[{runDir, "mt15_01_v2_22b_forbidden_R_norm_packet.json"}];
diagnosticsPath = FileNameJoin[{runDir, "mt15_01_geometry_diagnostics.json"}];

ClearAll[trapz, weightedNorm, weightedInner, normalizeProfile, smootherStep, toJsonNumber,
  sourceCommitString];

trapz[xs_List, ys_List] /; Length[xs] == Length[ys] && Length[xs] >= 2 :=
  Total[
    Table[
      0.5 (ys[[i]] + ys[[i + 1]]) (xs[[i + 1]] - xs[[i]]),
      {i, 1, Length[xs] - 1}
    ]
  ];

weightedInner[xs_List, weights_List, f_List, g_List] :=
  trapz[xs, weights f g];

weightedNorm[xs_List, weights_List, f_List] :=
  weightedInner[xs, weights, f, f];

normalizeProfile[xs_List, weights_List, f_List] := Module[{norm = weightedNorm[xs, weights, f]},
  If[! NumericQ[norm] || norm <= 0,
    Print["ERROR: nonpositive norm while normalizing profile: ", norm];
    Exit[2]
  ];
  N[f/Sqrt[norm], 17]
];

smootherStep[x_] := 3 x^2 - 2 x^3;
toJsonNumber[x_] := N[x, 17];

sourceCommitString[] := Module[
  {gitDir, headPath, head, refName, refPath, packedPath, packedLines, match, sha},
  gitDir = FileNameJoin[{repoRoot, ".git"}];
  headPath = FileNameJoin[{gitDir, "HEAD"}];
  If[! FileExistsQ[headPath], Return["unavailable"]];

  head = Quiet[Check[StringTrim[Import[headPath, "Text"]], ""]];
  If[! StringQ[head] || StringLength[head] == 0, Return["unavailable"]];

  sha = If[StringStartsQ[head, "ref:"],
    refName = StringTrim[StringDrop[head, 4]];
    refPath = FileNameJoin[Join[{gitDir}, FileNameSplit[refName]]];
    If[FileExistsQ[refPath],
      Quiet[Check[StringTrim[Import[refPath, "Text"]], ""]],
      packedPath = FileNameJoin[{gitDir, "packed-refs"}];
      If[! FileExistsQ[packedPath],
        "",
        packedLines = StringSplit[Quiet[Check[Import[packedPath, "Text"], ""]], "\n"];
        match = Select[packedLines, StringMatchQ[#, Except["#"] ~~ ___ ~~ " " ~~ refName] &];
        If[Length[match] > 0, First[StringSplit[First[match]]], ""]
      ]
    ],
    head
  ];

  If[StringQ[sha] && StringLength[sha] >= 7, StringTake[sha, 7], "unavailable"]
];

(* Placeholder, target-blind dimensional constants and geometry knobs. *)
L = 2.0;
a = 1.0;
RExit = 1.65;
meshPoints = 129;

GConst = 1.0;
csConst = 1.0;
cConst = 2.5;
mhat0Const = 0.85;
SPortConst = 1.2;
thetaTailConst = 1.0;

points = N[Table[L i/(meshPoints - 1), {i, 0, meshPoints - 1}], 17];
xGrid = points/L;
radiusValues = N[a + (RExit - a) (smootherStep /@ xGrid), 17];
weightValues = N[radiusValues^2, 17];

rawWall = N[Sin[Pi points/L], 17];
wallChiEta = normalizeProfile[points, weightValues, rawWall];

rawBdg = N[Sin[Pi points/(2 L)], 17];
bdgPhi = normalizeProfile[points, weightValues, rawBdg];

rawU = N[Cos[Pi points/(4 L)], 17];
uProfile = normalizeProfile[points, weightValues, rawU];

rawW = N[Sin[3 Pi points/(4 L)] + 0.15 Sin[Pi points/L], 17];
wProfile = normalizeProfile[points, weightValues, rawW];

wallNorm = weightedNorm[points, weightValues, wallChiEta];
bdgNorm = weightedNorm[points, weightValues, bdgPhi];
uNorm = weightedNorm[points, weightValues, uProfile];
wNorm = weightedNorm[points, weightValues, wProfile];
Iuw = weightedInner[points, weightValues, uProfile, wProfile];

lambdaR = 0.18;
omegaU = 3.1;
omegaW = 4.2;
deltaEff = omegaU^2 omegaW^2 - (lambdaR Iuw)^2;

geometryDiagnostics = <|
  "L" -> toJsonNumber[L],
  "mesh_points" -> meshPoints,
  "grid_first" -> toJsonNumber[First[points]],
  "grid_last" -> toJsonNumber[Last[points]],
  "grid_strictly_increasing" -> TrueQ[And @@ Thread[Differences[points] > 0]],
  "R0_first_equals_a" -> TrueQ[Abs[First[radiusValues] - a] < 10^-14],
  "R0_last_equals_R_exit" -> TrueQ[Abs[Last[radiusValues] - RExit] < 10^-14],
  "R_mouth" -> toJsonNumber[a],
  "R_exit" -> toJsonNumber[RExit],
  "R_exit_positive" -> TrueQ[RExit > 0],
  "R0_monotone_nondecreasing" -> TrueQ[And @@ Thread[Differences[radiusValues] >= -10^-14]],
  "R0_strictly_increasing_on_grid" -> TrueQ[And @@ Thread[Differences[radiusValues] > 0]],
  "weight_min" -> toJsonNumber[Min[weightValues]],
  "weight_max" -> toJsonNumber[Max[weightValues]],
  "wall_norm" -> toJsonNumber[wallNorm],
  "bdg_norm" -> toJsonNumber[bdgNorm],
  "u_norm" -> toJsonNumber[uNorm],
  "w_norm" -> toJsonNumber[wNorm],
  "I_u_w" -> toJsonNumber[Iuw],
  "lambda_R" -> toJsonNumber[lambdaR],
  "Omega_U" -> toJsonNumber[omegaU],
  "Omega_W" -> toJsonNumber[omegaW],
  "Delta_eff" -> toJsonNumber[deltaEff]
|>;

packet = <|
  "schema" -> "stage_v2_22b_solver_handoff/v1",
  "branch_id" -> "mt15_m1a_pre_freeze_placeholder_open_throat_machinery",
  "freeze" -> <|
    "pre_target_freeze" -> True,
    "target_blind" -> True,
    "parent_action" -> "GNLS + localized Maxwell + effective wall closure placeholder",
    "parent_action_status" -> "effective_wall_closure",
    "gauge_convention" -> "localized_H_equals_Z",
    "boundary_protocol" -> "open_impedance_AC_reflecting_DC_leaking",
    "source_map_convention" -> "canonical real STF grouped P2 basis",
    "support_profile_family" -> "stationary smooth open-throat placeholder profiles"
  |>,
  "geometry" -> <|
    "L" -> toJsonNumber[L],
    "R_mouth" -> toJsonNumber[a],
    "R_exit" -> toJsonNumber[RExit],
    "boundary_class" -> "open_impedance",
    "Y_L_limit" -> 0.0,
    "exit_model" -> "impedance_mismatch_open_exit"
  |>,
  "constants" -> <|
    "G" -> toJsonNumber[GConst],
    "c_s" -> toJsonNumber[csConst],
    "c" -> toJsonNumber[cConst],
    "a" -> toJsonNumber[a],
    "mhat0" -> toJsonNumber[mhat0Const],
    "S_port" -> toJsonNumber[SPortConst],
    "theta_tail" -> toJsonNumber[thetaTailConst]
  |>,
  "grid" -> <|
    "coordinate" -> "w",
    "points" -> (toJsonNumber /@ points)
  |>,
  "wall" -> <|
    "K" -> 2.25,
    "M" -> 0.72
  |>,
  "profiles" -> <|
    "radius_R0" -> (toJsonNumber /@ radiusValues),
    "weight" -> (toJsonNumber /@ weightValues),
    "wall_chi_eta" -> (toJsonNumber /@ wallChiEta)
  |>,
  "bdg_modes" -> {
    <|
      "name" -> "m1a_placeholder_bdg_halfwave",
      "lambda_B" -> 0.34,
      "varpi" -> 2.6,
      "profile_values" -> (toJsonNumber /@ bdgPhi)
    |>
  },
  "mixed_ports" -> {
    <|
      "name" -> "m1a_placeholder_open_port",
      "lambda_U" -> 0.24,
      "Omega_U" -> toJsonNumber[omegaU],
      "u_values" -> (toJsonNumber /@ uProfile),
      "lambda_W" -> 0.31,
      "Omega_W" -> toJsonNumber[omegaW],
      "w_values" -> (toJsonNumber /@ wProfile),
      "lambda_R" -> toJsonNumber[lambdaR]
    |>
  },
  "solver_metadata" -> <|
    "exporter" -> "mt15_01_geometry_stationary.wls",
    "mesh_points" -> meshPoints,
    "coefficient_family" -> "v2_22b_stationary_open_throat_placeholder_machinery_v1",
    "source_commit" -> sourceCommitString[],
    "export_timestamp" -> DateString[Now, {"ISODate", "T", "Time", " ", "TimeZoneNameShort"}],
    "stationary_profile_residual_norm" -> 0.0,
    "notes" -> "Pre-freeze target-blind machinery export. Geometry, modal, and port values are placeholders for the M1a handoff only and are to be replaced by coupled eigensolve data in M1b."
  |>,
  "normalization_tolerances" -> <|
    "profile_norm_tol" -> 0.005,
    "Delta_tol" -> 1.0*^-12
  |>
|>;

badPacket = Join[packet, <|"R_norm" -> 0.0|>];

Export[packetPath, packet, "JSON"];
Export[badPacketPath, badPacket, "JSON"];
Export[diagnosticsPath, geometryDiagnostics, "JSON"];

Print["MT15 M1a stationary open-throat packet export"];
Print["packet_path: ", packetPath];
Print["bad_packet_path: ", badPacketPath];
Print["diagnostics_path: ", diagnosticsPath];
Print["mesh_points: ", meshPoints];
Print["R0(0): ", ToString[NumberForm[First[radiusValues], {16, 12}]], "  a: ", ToString[NumberForm[a, {16, 12}]]];
Print["R0(L): ", ToString[NumberForm[Last[radiusValues], {16, 12}]], "  R_exit: ", ToString[NumberForm[RExit, {16, 12}]]];
Print["grid_strictly_increasing: ", geometryDiagnostics["grid_strictly_increasing"]];
Print["R0_monotone_nondecreasing: ", geometryDiagnostics["R0_monotone_nondecreasing"]];
Print["wall_norm: ", ToString[NumberForm[wallNorm, {16, 12}]]];
Print["bdg_norm: ", ToString[NumberForm[bdgNorm, {16, 12}]]];
Print["I_u_w: ", ToString[NumberForm[Iuw, {16, 12}]]];
Print["Delta_eff: ", ToString[NumberForm[deltaEff, {16, 12}]]];
