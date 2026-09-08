(* Manifest — copied from the orchestrator-owned knife-list:
   K_A: faceSources, flux assignment; drop lambdaAResponse affinity.
   K_T: faceSources, virtualWork product assignment; replace the whole product
        with 0.
   K_W: pressureField[-1]; replace pressureLower[xOne,xTwo,xThree,time]
        with pressureUpper[xOne,xTwo,xThree,time].
   EXTRACTOR_ORDER: apply the pressure-zero rules before differentiation.
   DEAD_PATH: evaluatedModel, THICKNESS_ROW construction; delete kineticEwLive.
   RESCALE_K_A: faceSources, flux assignment; multiply lambdaAResponse affinity
                by 2, retaining lambdaVResponse normalVelocity.
   RESCALE_K_T: faceSources, virtualWork product assignment; multiply by 2.
   RESCALE_K_W: pressureField[-1]; multiply pressureLower[...] by 2.

   1. The harness may **PRINT** computed objects (carriers, diffs). It may ⛔ NOT state conclusions — no `PASS`, no
      verdict, no "bites"/"fails". Interpretation is the orchestrator's, from the printed triples.
   2. **Print operand and residual, then guard.** Emit `baseline`, `corrupted`, and their `diff`; a residual asserted
      zero/nonzero carries no information.
   3. Interpretation belongs to the review / step record, ⛔ not the script.

   Invoke through ../scripts/S11c_b_carrier_ablation_harness_wl.py. Each worker
   uses one fresh kernel; the driver wraps every invocation in timeout 600.
   Only the source prefix above the production main-object marker is loaded.
*)

Begin["S11cbCarrierHarness`"];
$HistoryLength = 0;

emit[name_String, object_] := WriteString[First[$Output],
  name <> ": " <> ToString[object, InputForm, PageWidth -> Infinity] <> "\n"];

stop[reason_String, operands_:<||>] := (
  emit["HARNESS_GUARD", <|"REASON" -> reason, "OPERANDS" -> operands|>];
  Exit[2]
);

sha256[path_String] := IntegerString[FileHash[path, "SHA256"], 16, 64];
objectDigest[object_] := IntegerString[Hash[object, "SHA256"], 16, 64];

marker = "(* Main variable-coefficient objects.                                    *)";
slotNames = {"delta_p_plus", "d_w_delta_p_plus", "delta_p_minus",
  "d_w_delta_p_minus"};
rowNames = {"U_MOMENTUM_ROWS[1]", "U_MOMENTUM_ROWS[2]",
  "U_MOMENTUM_ROWS[3]", "MASS_EVOLUTION_ROW", "THICKNESS_ROW"};
notApplicable = Missing["NotApplicable", "NoInPlaneWPressureSlot"];
constructionModes = {"K_A", "K_T", "K_W", "DEAD_PATH",
  "RESCALE_K_A", "RESCALE_K_T", "RESCALE_K_W"};

manifest = <|
  "K_A" -> <|"FUNCTION" -> "faceSources",
    "FORM" -> "drop lambdaAResponse affinity from the flux assignment"|>,
  "K_T" -> <|"FUNCTION" -> "faceSources",
    "FORM" -> "set the whole virtualWork product to 0"|>,
  "K_W" -> <|"FUNCTION" -> "pressureField[-1]",
    "FORM" -> "redefine the lower face from pressureUpper[...]"|>|>;

(* Exact strings constrain both the function and the single assignment.
   A source mismatch terminates before any engine definitions are loaded. *)
fluxBefore = "  flux = lambdaAResponse affinity + lambdaVResponse normalVelocity;";
workBefore = "  virtualWork = virtualWork tractionPressure virtualNormalDisplacement;";
pressureBefore =
  "pressureField[-1] := pressureLower[xOne, xTwo, xThree, time];";
kineticBefore =
  "    \"THICKNESS_ROW\" -> (kineticEwLive + rowsLive[\"EW_INTERNAL\"] +\n" <>
  "      faceRowsLive[\"EW_FACE\"]),";

patches = <|
  "K_A" -> {"faceSources", fluxBefore,
    "  flux = lambdaVResponse normalVelocity;"},
  "K_T" -> {"faceSources", workBefore, "  virtualWork = 0;"},
  "K_W" -> {"pressureField[-1]", pressureBefore,
    "pressureField[-1] := pressureUpper[xOne, xTwo, xThree, time];"},
  "DEAD_PATH" -> {"evaluatedModel", kineticBefore,
    "    \"THICKNESS_ROW\" -> (rowsLive[\"EW_INTERNAL\"] +\n" <>
      "      faceRowsLive[\"EW_FACE\"]),"},
  "RESCALE_K_A" -> {"faceSources", fluxBefore,
    "  flux = 2 lambdaAResponse affinity + lambdaVResponse normalVelocity;"},
  "RESCALE_K_T" -> {"faceSources", workBefore,
    "  virtualWork = 2 virtualWork tractionPressure virtualNormalDisplacement;"},
  "RESCALE_K_W" -> {"pressureField[-1]", pressureBefore,
    "pressureField[-1] := 2 pressureLower[xOne, xTwo, xThree, time];"}|>;

functionBounds = <|
  "faceSources" -> {
    "faceSources[route_String, branch_String, sign_Integer,",
    "projectedFaceFlux[faceAssociation_Association]"},
  "evaluatedModel" -> {
    "evaluatedModel[route_String, branch_String, density_String,",
    "activateSpatialDivergences[expression_Association]"},
  "pressureField[-1]" -> {pressureBefore, "faceSources[route_String"}|>;

validateSites[source_String] := Module[{cut, bounds, start, finish, region},
  If[StringCount[source, marker] =!= 1,
    stop["definitions marker count", <|"COUNT" -> StringCount[source, marker]|>]];
  cut = First[First[StringPosition[source, marker]]];
  Do[
    bounds = functionBounds[spec[[1]]];
    start = StringPosition[source, bounds[[1]]];
    finish = StringPosition[source, bounds[[2]]];
    If[Length[start] =!= 1 || Length[finish] =!= 1,
      stop["named function bounds", <|"FUNCTION" -> spec[[1]],
        "START_COUNT" -> Length[start], "END_COUNT" -> Length[finish]|>]];
    If[!TrueQ[start[[1, 1]] < finish[[1, 1]] < cut],
      stop["named function position", <|"FUNCTION" -> spec[[1]]|>]];
    region = StringTake[source, {start[[1, 1]], finish[[1, 1]] - 1}];
    If[StringCount[source, spec[[2]]] =!= 1 ||
        StringCount[region, spec[[2]]] =!= 1,
      stop["construction assignment count", <|"FUNCTION" -> spec[[1]],
        "SITE" -> spec[[2]], "SOURCE_COUNT" -> StringCount[source, spec[[2]]],
        "FUNCTION_COUNT" -> StringCount[region, spec[[2]]]|>]],
    {spec, Values[patches]}];
  cut
];

nativeAtoms[] := {
  Global`pressureUpper[Global`xOne, Global`xTwo, Global`xThree, Global`time],
  Global`pressureLower[Global`xOne, Global`xTwo, Global`xThree, Global`time]};

scalarRows[operator_Association] := Join[operator["U_MOMENTUM_ROWS"],
  {operator["MASS_EVOLUTION_ROW"], operator["THICKNESS_ROW"]}];

extractCarrier[operator_Association, zeroFirst_:False] := Module[
  {atoms = nativeAtoms[], zeroRules, rows, coefficients},
  zeroRules = Thread[atoms -> 0];
  rows = scalarRows[operator];
  coefficients = Table[
    If[TrueQ[zeroFirst], D[rows[[r]] /. zeroRules, atoms[[p]]],
      D[rows[[r]], atoms[[p]]] /. zeroRules],
    {r, Length[rows]}, {p, Length[atoms]}];
  AssociationThread[rowNames, Map[Function[pair,
    AssociationThread[slotNames,
      {pair[[1]], notApplicable, pair[[2]], notApplicable}]], coefficients]]
];

carrierDifference[corrupted_Association, baseline_Association] :=
  AssociationMap[Function[row, AssociationMap[Function[slot,
    If[MemberQ[{"d_w_delta_p_plus", "d_w_delta_p_minus"}, slot],
      notApplicable, Expand[corrupted[row][slot] - baseline[row][slot]]]],
    slotNames]], rowNames];

carrierShapeQ[carrier_] := AssociationQ[carrier] &&
  Keys[carrier] === rowNames && AllTrue[Values[carrier],
    AssociationQ[#] && Keys[#] === slotNames &];

emitTriple[name_String, baseline_, corrupted_] := Module[{diff},
  diff = carrierDifference[corrupted, baseline];
  emit[name, <|"baseline" -> baseline, "corrupted" -> corrupted,
    "diff" -> diff|>];
  If[!And @@ (carrierShapeQ /@ {baseline, corrupted, diff}),
    stop["carrier shape", <|"TRIPLE" -> name|>]];
  If[!FreeQ[{baseline, corrupted, diff}, $Aborted | $Failed | Indeterminate],
    stop["carrier evaluation", <|"TRIPLE" -> name|>]]
];

worker[mode_String, enginePath_String, directory_String] := Module[
  {source, copiedPath, inputPath, cut, prefixPath, spec, operator, carrier,
    zeroFirstCarrier, record, exportPath},
  source = Import[enginePath, "Text"];
  If[!StringQ[source], stop["engine source read"]];
  cut = validateSites[source];
  inputPath = enginePath;
  If[mode =!= "BASELINE",
    copiedPath = FileNameJoin[{directory, mode <> "_engine.wl"}];
    If[mode === "UNABLATED_COPY", CopyFile[enginePath, copiedPath],
      If[!KeyExistsQ[patches, mode], stop["worker mode", <|"MODE" -> mode|>]];
      spec = patches[mode];
      source = StringReplace[source, spec[[2]] -> spec[[3]]];
      Export[copiedPath, source, "Text"];
      emit["CONSTRUCTION_PATCH", <|"MODE" -> mode,
        "FUNCTION" -> spec[[1]], "baseline_source" -> spec[[2]],
        "corrupted_source" -> spec[[3]]|>]];
    inputPath = copiedPath;
    source = Import[inputPath, "Text"];
    cut = First[First[StringPosition[source, marker]]]
  ];
  prefixPath = FileNameJoin[{directory, mode <> "_definitions.wl"}];
  Export[prefixPath, StringTake[source, cut - 1], "Text"];
  emit["DEFINITIONS_INPUT", <|"MODE" -> mode,
    "ENGINE_SOURCE_SHA256" -> sha256[inputPath],
    "DEFINITIONS_SHA256" -> sha256[prefixPath],
    "ENTRYPOINT" -> "evaluatedModel[\"EULERIAN\", \"MATERIAL_ADVECTED\", \"RHO4_CONSTANT\"][\"OPERATOR\"]"|>];
  Block[{$Context = "Global`", $ContextPath = {"System`", "Global`"}},
    Get[prefixPath]];
  operator = Global`evaluatedModel["EULERIAN", "MATERIAL_ADVECTED",
    "RHO4_CONSTANT"]["OPERATOR"];
  If[!AssociationQ[operator], stop["operator evaluation"]];
  carrier = extractCarrier[operator];
  record = <|"MODE" -> mode, "OPERATOR_SHA256" -> objectDigest[operator],
    "CARRIER" -> carrier|>;
  If[mode === "BASELINE",
    zeroFirstCarrier = extractCarrier[operator, True];
    AssociateTo[record, "ZERO_FIRST_CARRIER" -> zeroFirstCarrier]];
  emit["COMPUTED_OBJECT", record];
  If[!carrierShapeQ[carrier], stop["worker carrier shape"]];
  exportPath = FileNameJoin[{directory, mode <> ".wxf"}];
  Export[exportPath, record, "WXF"];
  If[!FileExistsQ[exportPath], stop["worker export"]]
];

report[directory_String] := Module[{baseline, copy, current, digestComparison},
  baseline = Import[FileNameJoin[{directory, "BASELINE.wxf"}], "WXF"];
  copy = Import[FileNameJoin[{directory, "UNABLATED_COPY.wxf"}], "WXF"];
  emit["MANIFEST", manifest];
  emit["PRESSURE_SLOTS", AssociationThread[slotNames,
    {nativeAtoms[][[1]], notApplicable, nativeAtoms[][[2]], notApplicable}]];
  digestComparison = <|"baseline" -> baseline["OPERATOR_SHA256"],
    "corrupted" -> copy["OPERATOR_SHA256"],
    "same_object_digest" -> SameQ[baseline["OPERATOR_SHA256"],
      copy["OPERATOR_SHA256"]],
    "same_carrier" -> SameQ[baseline["CARRIER"], copy["CARRIER"]]|>;
  emit["CANONICAL_COPY_OBJECT_COMPARISON", digestComparison];
  emitTriple["CANONICAL_COPY_CARRIER", baseline["CARRIER"], copy["CARRIER"]];
  If[!TrueQ[digestComparison["same_object_digest"]] ||
      !TrueQ[digestComparison["same_carrier"]], stop["canonical copy drift"]];
  Do[
    current = Import[FileNameJoin[{directory, mode <> ".wxf"}], "WXF"];
    emitTriple[mode, baseline["CARRIER"], current["CARRIER"]],
    {mode, Take[constructionModes, 3]}];
  emitTriple["EXTRACTOR_ORDER", baseline["CARRIER"],
    baseline["ZERO_FIRST_CARRIER"]];
  Do[
    current = Import[FileNameJoin[{directory, mode <> ".wxf"}], "WXF"];
    emitTriple[mode, baseline["CARRIER"], current["CARRIER"]],
    {mode, Drop[constructionModes, 3]}]
];

mode = Environment["S11CB_CARRIER_MODE"];
directory = Environment["S11CB_CARRIER_WORKDIR"];
enginePath = Environment["S11CB_CARRIER_ENGINE"];
If[!StringQ[mode] || !StringQ[directory] || !StringQ[enginePath],
  stop["driver environment"]];
If[mode === "REPORT", report[directory], worker[mode, enginePath, directory]];
End[];
Exit[0];
