(* WL sibling provenance probe (rule 16): where does the order-3 background jet
   in WL's U_MOMENTUM_ROWS come from — the BULK slot (would be a real disagreement
   with PY's order-2-complete strong rows) or the FACE/KINETIC slot (an object-boundary
   artifact the row_residual comparator already excludes)?

   PRINTS the engine's own backgroundJetOrder over each provenance slot; states nothing.
   Reduced scale (basisRepresentativeIndices = {16}), one case: EULERIAN/LAB_HELD/RHO4_CONSTANT.
   Loader = the proven StringTake-between-markers pattern from the #89b re-review harness. *)

ClearAll["Global`*"];
$HistoryLength = 0;
$Messages = {OutputStream["stderr", 2]};

enginePath = "/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl";
engineSource = Import[enginePath, "Text"];
markerStart[marker_String] := First[First[StringPosition[engineSource, marker]]];

definitionStart = markerStart["sharedObject[name_String]"];
rankSetupStart = markerStart["baseEnergyDataByBranch ="];
postRankStart = markerStart["epsilonPseudoscalarCandidates["];
mainExecutionStart = markerStart["(* Main variable-coefficient objects."];
commonDefinitions = StringTake[engineSource, {definitionStart, rankSetupStart - 1}] <>
  "basisRepresentativeIndices = {16};\n" <>
  StringTake[engineSource, {postRankStart, mainExecutionStart - 1}];
ToExpression[commonDefinitions, InputForm];

printLiteral[label_String, value_] := WriteString[First[$Output],
  label <> ": " <> ToString[value, InputForm, PageWidth -> Infinity] <> "\n"];

processed = evaluatedModel["EULERIAN", "LAB_HELD", "RHO4_CONSTANT"];
origins = processed["ORIGINS"];                       (* activated-then-reduced origins *)
fullU = processed["OPERATOR"]["U_MOMENTUM_ROWS"];     (* = kineticU + U_INTERNAL + U_FACE, reduced *)

slotU[slot_] := Lookup[Lookup[origins, slot, <||>], "U_MOMENTUM_ROWS", {0, 0, 0}];

report[label_String, expr_] := Module[{atoms, maxOrder, order3},
  atoms = backgroundJetAtomsIn[expr];
  maxOrder = If[atoms === {}, 0, Max[backgroundJetOrder /@ atoms]];
  order3 = Select[atoms, backgroundJetOrder[#] >= 3 &];
  printLiteral[label <> "_MAX_JET_ORDER", maxOrder];
  printLiteral[label <> "_N_ATOMS", Length[atoms]];
  printLiteral[label <> "_ORDER3PLUS_ATOMS", order3];
];

printLiteral["ORIGINS_SLOT_KEYS", Keys[origins]];
report["FULL_U_MOMENTUM_ROWS", fullU];       (* cross-check: should show the order-3 the build-review saw *)
report["ORIGIN_KINETIC_U", slotU["KINETIC"]];
report["ORIGIN_BULK_ENERGY_U", slotU["BULK_ENERGY"]];  (* <- maps to PY strong rows; hypothesis: <=2 *)
report["ORIGIN_FACE_U", slotU["FACE"]];                (* <- hypothesis: carries the order-3 *)

(* provenance reconstruction check: do the slots sum to the full U row? *)
reconstruction = fullU - (slotU["KINETIC"] + slotU["BULK_ENERGY"] + slotU["FACE"]);
printLiteral["FULL_MINUS_SLOTSUM_U (expect 0 vector)", Simplify[reconstruction]];

WriteString[First[$Output], "DONE\n"];
Quit[];
