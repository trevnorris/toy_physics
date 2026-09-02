(* Compare the raw energy U Euler derivative with WL's constrained U_INTERNAL. *)
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
report[label_String, expr_] := Module[{atoms, orders, order3},
  atoms = backgroundJetAtomsIn[expr];
  orders = backgroundJetOrder /@ atoms;
  order3 = Select[atoms, backgroundJetOrder[#] >= 3 &];
  printLiteral[label <> "_MAX_JET_ORDER", If[orders === {}, 0, Max[orders]]];
  printLiteral[label <> "_ORDER3PLUS_ATOMS", order3];
];

energyData = constructEnergyData["LAB_HELD", widthBase, modulusBase];
selectedRecord = energyData["RECORDS"][[16]];
energyLive = selectedRecord["ENERGY_TERM"];
constraintLive = virtualConstraintSource[
  "EULERIAN", "LAB_HELD", "RHO4_CONSTANT", widthBase, False];

rawULive = activateSpatialDivergences[
  variationalSource[energyLive, #]] & /@ uVector;
muThetaLive = activateSpatialDivergences[
  variationalSource[energyLive, thetaField]];
constrainedLive = constrainedRowsWithLiveEnergyEL[
  energyLive, constraintLive]["U_INTERNAL"];

energy = finalBackgroundReduction[energyLive];
uFlux = finalBackgroundReduction[Flatten[Table[
  D[energyLive, D[uVector[[a]], spatialCoordinates[[i]]]],
  {a, directions}, {i, directions}]]];
rawU = finalBackgroundReduction[rawULive];
muTheta = finalBackgroundReduction[muThetaLive];
constrainedU = finalBackgroundReduction[
  activateSpatialDivergences[constrainedLive]];
correction = Simplify[constrainedU - rawU];

printLiteral["SELECTED_RECORD_LABEL", selectedRecord["LABEL"]];
report["ENERGY_DENSITY", energy];
report["DU_COEFFICIENT", uFlux];
report["RAW_U_ENERGY_EL", rawU];
report["THETA_ENERGY_EL", muTheta];
report["WL_U_INTERNAL_CONSTRAINED", constrainedU];
report["CONSTRAINED_MINUS_RAW_U", correction];
printLiteral["RAW_EQUALS_U_INTERNAL", Simplify[rawU - constrainedU]];
WriteString[First[$Output], "DONE\n"];
Quit[];
