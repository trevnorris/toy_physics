(* WL probe: split the order-3 U-row into
     (i)  the energy density
     (ii) the unconstrained energy EL  variationalSource[energy, u] activated
     (iii) constrained U_INTERNAL
     (iv)  mu_theta and grad(mu_theta)
   Reduced scale: basisRepresentativeIndices = {16} (same as the
   provenance probe).  Prints operands; states no conclusion. *)

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

report[label_String, expr_] := Module[{atoms, maxOrder, order3, liveDepth},
  liveDepth = generatedBackgroundDerivativeDepth[expr];
  atoms = backgroundJetAtomsIn[reduceExpr[expr]];
  maxOrder = If[atoms === {}, 0, Max[backgroundJetOrder /@ atoms]];
  order3 = Select[atoms, backgroundJetOrder[#] >= 3 &];
  printLiteral[label <> "_LIVE_DERIVATIVE_DEPTH", liveDepth];
  printLiteral[label <> "_REDUCED_MAX_JET_ORDER", maxOrder];
  printLiteral[label <> "_REDUCED_N_ATOMS", Length[atoms]];
  printLiteral[label <> "_REDUCED_ORDER3PLUS_ATOMS", order3];
];

energyData = constructEnergyData["LAB_HELD", widthBase, modulusBase];
printLiteral["RECORD_16_LABEL", energyData["RECORDS"][[16, "LABEL"]]];
printLiteral["RECORD_16_ORIGIN", energyData["RECORDS"][[16, "ORIGIN"]]];

energyLive = energyData["RECORDS"][[16, "ENERGY_TERM"]];
report["ENERGY_DENSITY", energyLive];

(* Unconstrained first-order energy EL for u, activated once. *)
uUnconstHeld = variationalSource[energyLive, #] & /@ uVector;
uUnconstAct = activateSpatialDivergences[uUnconstHeld];
report["U_EL_UNCONSTRAINED_HELD", uUnconstHeld];
report["U_EL_UNCONSTRAINED_ACTIVATED", uUnconstAct];

(* theta EL *)
muHeld = variationalSource[energyLive, thetaField];
muAct = activateSpatialDivergences[muHeld];
report["MU_THETA_HELD", muHeld];
report["MU_THETA_ACTIVATED", muAct];

gradMu = D[muAct, #] & /@ spatialCoordinates;
report["GRAD_MU_THETA", gradMu];

(* Constrained U_INTERNAL via the engine constructor (same path as evaluatedModel). *)
constraintLive = virtualConstraintSource["EULERIAN", "LAB_HELD", "RHO4_CONSTANT",
  widthBase, False];
rowsLive = constrainedRowsWithLiveEnergyEL[energyLive, constraintLive];
report["U_INTERNAL_CONSTRAINED", rowsLive["U_INTERNAL"]];
report["MU_THETA_FROM_CONSTRAINED_ROWS", rowsLive["MU_THETA"]];

(* Difference: constrained U minus unconstrained U. *)
diffCU = rowsLive["U_INTERNAL"] - uUnconstAct;
report["CONSTRAINED_MINUS_UNCONSTRAINED_U", diffCU];

(* Does constrained - unconstrained equal -grad(mu) up to lower-order terms?
   Residual against -grad mu. *)
diffVsGradMu = rowsLive["U_INTERNAL"] - uUnconstAct + gradMu;
report["UINTERNAL_MINUS_UEL_PLUS_GRADMU", diffVsGradMu];

(* Second activation of already-activated unconstrained EL (idempotence). *)
uUnconstTwice = activateSpatialDivergences[uUnconstAct];
report["U_EL_UNCONSTRAINED_DOUBLE_ACTIVATED", uUnconstTwice];
printLiteral["DOUBLE_ACTIVATE_MINUS_SINGLE (expect 0)",
  Simplify[uUnconstTwice - uUnconstAct]];

WriteString[First[$Output], "DONE\n"];
Quit[];
