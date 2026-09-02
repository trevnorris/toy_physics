(* Follow-up: activate the HELD constraint Div on U_INTERNAL and compare
   to unconstrained EL and to -grad(mu_theta).  Same reduced record 16.
   Prints operands; states no conclusion. *)

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

reduceExpr[expr_] := finalBackgroundReduction[expr, "BASE"];

report[label_String, expr_] := Module[{atoms, maxOrder, order3, liveDepth, nHeldDiv},
  liveDepth = generatedBackgroundDerivativeDepth[expr];
  nHeldDiv = Count[expr, HoldPattern[Inactive[Div][___]], {0, Infinity}];
  atoms = backgroundJetAtomsIn[reduceExpr[expr]];
  maxOrder = If[atoms === {}, 0, Max[backgroundJetOrder /@ atoms]];
  order3 = Select[atoms, backgroundJetOrder[#] >= 3 &];
  printLiteral[label <> "_HELD_DIV_COUNT", nHeldDiv];
  printLiteral[label <> "_LIVE_DERIVATIVE_DEPTH", liveDepth];
  printLiteral[label <> "_REDUCED_MAX_JET_ORDER", maxOrder];
  printLiteral[label <> "_REDUCED_N_ATOMS", Length[atoms]];
  printLiteral[label <> "_REDUCED_ORDER3PLUS_ATOMS", order3];
];

energyData = constructEnergyData["LAB_HELD", widthBase, modulusBase];
energyLive = energyData["RECORDS"][[16, "ENERGY_TERM"]];
printLiteral["RECORD_16_LABEL", energyData["RECORDS"][[16, "LABEL"]]];

constraintLive = virtualConstraintSource["EULERIAN", "LAB_HELD", "RHO4_CONSTANT",
  widthBase, False];
rowsLive = constrainedRowsWithLiveEnergyEL[energyLive, constraintLive];

uUnconstAct = activateSpatialDivergences[
  variationalSource[energyLive, #] & /@ uVector];
muAct = activateSpatialDivergences[variationalSource[energyLive, thetaField]];
gradMu = D[muAct, #] & /@ spatialCoordinates;

printLiteral["CONSTRAINT_LINEAR_IN_VIRTUALS",
  Expand[constraintLive]];

report["U_INTERNAL_BEFORE_SECOND_ACTIVATE", rowsLive["U_INTERNAL"]];
uInternalAct = activateSpatialDivergences[rowsLive["U_INTERNAL"]];
report["U_INTERNAL_AFTER_SECOND_ACTIVATE", uInternalAct];

report["U_EL_UNCONSTRAINED_ACTIVATED", uUnconstAct];
report["GRAD_MU_THETA", gradMu];

(* constrained_activated - unconstrained *)
diff = uInternalAct - uUnconstAct;
report["ACTIVATED_UINTERNAL_MINUS_UNCONSTRAINED", diff];

(* constrained_activated - unconstrained + grad mu  (if U_INTERNAL absorbed -grad mu) *)
diffPlus = uInternalAct - uUnconstAct + gradMu;
report["UINTERNAL_MINUS_UEL_PLUS_GRADMU", diffPlus];

(* How many Inactive[Div] were sitting in U_INTERNAL before second activate? *)
printLiteral["U_INTERNAL_HELD_DIV_VECTORS",
  Cases[rowsLive["U_INTERNAL"], HoldPattern[Inactive[Div][v_, _]] :> v,
    {0, Infinity}]];

WriteString[First[$Output], "DONE\n"];
Quit[];
