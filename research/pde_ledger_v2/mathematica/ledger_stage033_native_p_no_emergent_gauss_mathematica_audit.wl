(* Ledger stage033 Mathematica audit: native P has no emergent Gauss law.

   Standalone, print-only, assert-zero, machine-real-free, and file-I/O-free.
   This is a re-authored Wolfram route, not the source buildH2/diracSearch/
   searchG mirror.  A single native fixed-point constraint pipeline consumes
   independently constructed Lagrangians.  Reduce[ForAll[...]] and Solve
   derive the common Hessian/potential-null locus before it is used by any
   tuned scan.  Every control and every control ablation enters that same
   pipeline.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE033_MUTATION";
activeMutation = Environment[mutationEnvironment];
If[!StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];
verdictToken = "NATIVE_P_NO_EMERGENT_GAUSS";
kVector = {1, 2, 3};
kSquared = kVector.kVector;

toothOrder = {
  "PASS_Q1_LAGRANGIAN_LEGENDRE",
  "PASS_Q2_MAXWELL_COMPUTED",
  "PASS_Q3_SIX_CONTROLS_SHARED_PIPELINE",
  "PASS_Q4_CAS_COMMON_NULL_RECONSTRUCTION",
  "PASS_Q5_GAUSS_SEARCH_LIVE",
  "PASS_Q6_DIRAC_CLOSURE",
  "PASS_GUARD_COUPLINGS_ENTER_A",
  "PASS_GUARD_COUPLINGS_ENTER_C",
  "PASS_THEORY_A_SIGNATURE",
  "PASS_THEORY_C_SIGNATURE",
  "PASS_KINETIC_HESSIAN_DETERMINANT",
  "PASS_COMMON_NULL_A",
  "PASS_COMMON_NULL_C",
  "PASS_BOUNDARY_SCAN_A",
  "PASS_BOUNDARY_SCAN_C",
  "PASS_RANDOMIZED_SWEEP_A",
  "PASS_RANDOMIZED_SWEEP_C",
  "PASS_HARDENING_DESCENDANT_A",
  "PASS_HARDENING_DESCENDANT_C",
  "PASS_PRIMARY_FC_ACCOUNTING_A",
  "PASS_PRIMARY_FC_ACCOUNTING_C",
  "PASS_GUARD_SEARCH_CAPABLE",
  "PASS_CONTROL_MAXWELL",
  "PASS_CONTROL_GAUGED_HARD_UNIT",
  "PASS_CONTROL_BARE_SIGMA",
  "PASS_CONTROL_NONCONSERVED_CURRENT",
  "PASS_CONTROL_COULOMB_GAUGE",
  "PASS_CONTROL_GLOBAL_U1",
  "PASS_SOURCE_FIRST_SHEAR_BOOKKEEPING",
  "PASS_HONEST_TUNED_SCOPE",
  "PASS_DECISION_ORDER_BRANCH2",
  "PASS_VERDICT_TOTALITY",
  "PASS_SOURCE_PREDICATE_MANIFEST"
};

ablationDescriptions = AssociationThread[
  toothOrder, ("tooth-local computed-input mutation for " <> #) & /@ toothOrder];
ablationDescriptions["PASS_Q1_LAGRANGIAN_LEGENDRE"] =
  "inject a retained velocity into the computed Hamiltonian representative";
ablationDescriptions["PASS_Q2_MAXWELL_COMPUTED"] = "give Maxwell A0 a kinetic term";
ablationDescriptions["PASS_Q3_SIX_CONTROLS_SHARED_PIPELINE"] = "remove global U(1) from the control table";
ablationDescriptions["PASS_Q4_CAS_COMMON_NULL_RECONSTRUCTION"] = "corrupt a Solve-derived common-null rule";
ablationDescriptions["PASS_Q5_GAUSS_SEARCH_LIVE"] = "replace Maxwell by kinetic-A0 Maxwell";
ablationDescriptions["PASS_Q6_DIRAC_CLOSURE"] = "cap the native-A fixed point before closure";
ablationDescriptions["PASS_GUARD_COUPLINGS_ENTER_A"] = "drop gbA from the A input Lagrangian";
ablationDescriptions["PASS_GUARD_COUPLINGS_ENTER_C"] = "drop gdC from the C input Lagrangian";
ablationDescriptions["PASS_THEORY_A_SIGNATURE"] = "give lambdaA a kinetic term";
ablationDescriptions["PASS_THEORY_C_SIGNATURE"] = "give lambdaC a kinetic term";
ablationDescriptions["PASS_KINETIC_HESSIAN_DETERMINANT"] = "perturb the u1 kinetic input";
ablationDescriptions["PASS_COMMON_NULL_A"] = "flip the derived gsA common-null rule";
ablationDescriptions["PASS_COMMON_NULL_C"] = "flip the derived gsC common-null rule";
ablationDescriptions["PASS_BOUNDARY_SCAN_A"] = "detune an A common-null boundary input";
ablationDescriptions["PASS_BOUNDARY_SCAN_C"] = "detune a C common-null boundary input";
ablationDescriptions["PASS_RANDOMIZED_SWEEP_A"] = "replace one A sample by the common-null input";
ablationDescriptions["PASS_RANDOMIZED_SWEEP_C"] = "replace one C sample by the common-null input";
ablationDescriptions["PASS_HARDENING_DESCENDANT_A"] = "inject a computed Maxwell k-gradient descendant";
ablationDescriptions["PASS_HARDENING_DESCENDANT_C"] = "inject a computed Maxwell k-gradient descendant";
ablationDescriptions["PASS_PRIMARY_FC_ACCOUNTING_A"] = "perturb the A primary-FC accounting";
ablationDescriptions["PASS_PRIMARY_FC_ACCOUNTING_C"] = "perturb the C primary-FC accounting";
ablationDescriptions["PASS_GUARD_SEARCH_CAPABLE"] = "force Maxwell onto its no-primary chain";
ablationDescriptions["PASS_CONTROL_MAXWELL"] = "Maxwell gains A0 kinetic energy";
ablationDescriptions["PASS_CONTROL_GAUGED_HARD_UNIT"] = "remove the hard-unit multiplier";
ablationDescriptions["PASS_CONTROL_BARE_SIGMA"] = "remove the hard-sigma multiplier";
ablationDescriptions["PASS_CONTROL_NONCONSERVED_CURRENT"] = "impose the continuity rule";
ablationDescriptions["PASS_CONTROL_COULOMB_GAUGE"] = "drop Coulomb gauge fixing";
ablationDescriptions["PASS_CONTROL_GLOBAL_U1"] = "add a local Maxwell sector";
ablationDescriptions["PASS_SOURCE_FIRST_SHEAR_BOOKKEEPING"] = "feed Maxwell into source/shear bookkeeping";
ablationDescriptions["PASS_HONEST_TUNED_SCOPE"] = "claim exhaustive tuned stratification";
ablationDescriptions["PASS_DECISION_ORDER_BRANCH2"] = "feed first-class Maxwell into the quadratic decision";
ablationDescriptions["PASS_VERDICT_TOTALITY"] = "feed computed Maxwell through the verdict derivation";
ablationDescriptions["PASS_SOURCE_PREDICATE_MANIFEST"] = "drop one canonical manifest row";

raise[name_] := Throw[name, "ledgerStage033Failure"];

assertExact[name_, expression_] := Module[{reals},
  reals = Cases[Unevaluated[expression], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FIRST_FAILURE=", name];
    Print["FAIL  ", name, ": machine-real atom(s): ", InputForm[reals]];
    raise[name]
  ]
];

cleanZero[expression_] := FullSimplify[expression] /. ConditionalExpression[0, _] -> 0;

expectZero[name_, residual_, evidence_: None] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FIRST_FAILURE=", name];
    Print["FAIL  ", name, ": residual = ", InputForm[clean]];
    If[evidence =!= None, Print["      evidence = ", InputForm[evidence]]];
    raise[name]
  ]
];

expectBool[name_, condition_, evidence_: None] :=
  expectZero[name, If[TrueQ[condition], 0, 1], evidence];

heading[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);

poissonBracket[f_, g_, pairs_List] := Factor@Total[
  (D[f, #[[1]]] D[g, #[[2]]] - D[f, #[[2]]] D[g, #[[1]]]) & /@ pairs];

independentPositions[rows_List] := Module[{chosen = {}, oldRank = 0, trialRank},
  Do[
    trialRank = MatrixRank[rows[[Append[chosen, position]]]];
    If[trialRank > oldRank, AppendTo[chosen, position]; oldRank = trialRank],
    {position, Length[rows]}];
  chosen
];

independentConstraintQ[expression_, constraints_List, phase_List] := Module[
  {candidateRow, oldRows},
  candidateRow = {D[expression, #] & /@ phase};
  If[AllTrue[First[candidateRow], TrueQ[FullSimplify[# == 0]] &], Return[False]];
  oldRows = Table[D[constraints[[i]], phase[[j]]],
    {i, Length[constraints]}, {j, Length[phase]}];
  MatrixRank[Join[oldRows, candidateRow]] > MatrixRank[oldRows]
];

proportionalToKQ[action_List] := Length[action] == 3 &&
  !AllTrue[action, TrueQ[FullSimplify[# == 0]] &] &&
  And @@ Flatten@Table[
    TrueQ[FullSimplify[action[[i]] kVector[[j]] - action[[j]] kVector[[i]] == 0]],
    {i, 3}, {j, i + 1, 3}];

(* A single fixed-point descriptor pipeline replaces the source's three
   one-for-one mirrored functions.  It uses native D/NullSpace/Inverse/
   MatrixRank throughout, and performs the Gauss classification in the same
   returned descriptor. *)
constraintPipeline[model_Association, iterationCap_: Automatic] := Module[
  {q, v, lagrangian, count, momenta, momentumMap, hessian, hessianNulls,
   zeroVelocity, primaryResidual, primaries, hessianRank, rowIDs, columnIDs,
   velocityValues, solvedVelocities, canonicalHamiltonian, phase, pairs,
   cap, initialState, step, finalState, constraints, stages, bracketMatrix,
   bracketRank, bracketKernel, classes, primaryCount, primaryDirections,
   accepted = {}, rejected = {}, coefficient, primary, descendant,
   primaryAction, secondaryAction, firstClassDescendant, scalarEnd,
   onlyScalarEnd, spatialGradient, nontrivialDescendant, descendantZero,
   rejectionCertified, reason, item},

  q = model["Coordinates"];
  v = model["Velocities"];
  lagrangian = model["Lagrangian"];
  count = Length[q];
  If[Length[v] =!= count, Throw["coordinate/velocity mismatch", "pipelineFailure"]];
  momenta = Array[Unique["stage033Pi$"] &, count];
  momentumMap = D[lagrangian, #] & /@ v;
  hessian = Table[D[momentumMap[[i]], v[[j]]], {i, count}, {j, count}];
  hessianNulls = NullSpace[Transpose[hessian]];
  zeroVelocity = Thread[v -> ConstantArray[0, count]];
  primaryResidual = momenta - (momentumMap /. zeroVelocity);
  primaries = Factor[#.primaryResidual] & /@ hessianNulls;
  hessianRank = MatrixRank[hessian];
  rowIDs = independentPositions[hessian];
  columnIDs = independentPositions[Transpose[hessian]];
  If[Length[rowIDs] =!= hessianRank || Length[columnIDs] =!= hessianRank,
    Throw["Hessian minor selection", "pipelineFailure"]];
  velocityValues = AssociationThread[v, ConstantArray[0, count]];
  If[hessianRank > 0,
    solvedVelocities = Inverse[hessian[[rowIDs, columnIDs]]].primaryResidual[[rowIDs]];
    Do[AssociateTo[velocityValues, v[[columnIDs[[i]]]] -> solvedVelocities[[i]]],
      {i, hessianRank}]
  ];
  canonicalHamiltonian = Factor[(momenta.v - lagrangian) /. Normal[velocityValues]];
  phase = Join[q, momenta];
  pairs = Transpose[{q, momenta}];
  primaryCount = Length[primaries];
  cap = If[iterationCap === Automatic, 2 Length[phase] + 2, iterationCap];
  initialState = <|"Constraints" -> primaries, "Stages" -> ConstantArray[0, primaryCount],
    "Iteration" -> 1, "Closed" -> False|>;

  step[state_Association] := Module[
    {current, multiplierMap, flow, compatibility, candidates, new = {}},
    If[TrueQ[state["Closed"]], Return[state]];
    current = state["Constraints"];
    If[current === {}, Return[Join[state, <|"Closed" -> True|>]]];
    multiplierMap = Table[
      poissonBracket[current[[i]], primaries[[j]], pairs],
      {i, Length[current]}, {j, primaryCount}];
    flow = poissonBracket[#, canonicalHamiltonian, pairs] & /@ current;
    compatibility = NullSpace[Transpose[multiplierMap]];
    candidates = Factor[#.flow] & /@ compatibility;
    Do[If[independentConstraintQ[candidate, Join[current, new], phase],
      AppendTo[new, candidate]], {candidate, candidates}];
    If[new === {},
      Join[state, <|"Closed" -> True|>],
      <|"Constraints" -> Join[current, new],
        "Stages" -> Join[state["Stages"], ConstantArray[state["Iteration"], Length[new]]],
        "Iteration" -> state["Iteration"] + 1, "Closed" -> False|>
    ]
  ];

  finalState = FixedPoint[step, initialState, cap];
  constraints = finalState["Constraints"];
  stages = finalState["Stages"];
  bracketMatrix = Table[poissonBracket[constraints[[i]], constraints[[j]], pairs],
    {i, Length[constraints]}, {j, Length[constraints]}];
  bracketRank = If[constraints === {}, 0, MatrixRank[bracketMatrix]];
  bracketKernel = If[constraints === {}, {}, NullSpace[bracketMatrix]];
  classes = Table[
    If[AllTrue[bracketMatrix[[i]], TrueQ[FullSimplify[# == 0]] &],
      "FIRST_CLASS", "SECOND_CLASS_COMPONENT"],
    {i, Length[constraints]}];
  primaryDirections = If[constraints === {} || primaryCount == 0, {},
    NullSpace[bracketMatrix[[All, 1 ;; primaryCount]]]];

  Do[
    primary = Factor@Sum[coefficient[[i]] constraints[[i]], {i, primaryCount}];
    descendant = Factor@poissonBracket[primary, canonicalHamiltonian, pairs];
    primaryAction = poissonBracket[#, primary, pairs] & /@ q;
    secondaryAction = poissonBracket[#, descendant, pairs] & /@ q;
    firstClassDescendant = AllTrue[constraints,
      TrueQ[FullSimplify[poissonBracket[descendant, #, pairs] == 0]] &];
    scalarEnd = AnyTrue[model["ScalarEnd"],
      TrueQ[FullSimplify[primaryAction[[#]] != 0]] &];
    onlyScalarEnd = scalarEnd && And @@ Table[
      If[MemberQ[Join[model["ScalarEnd"], model["Multipliers"]], i], True,
        TrueQ[FullSimplify[primaryAction[[i]] == 0]]], {i, Length[q]}];
    spatialGradient = AnyTrue[model["SpatialBlocks"],
      proportionalToKQ[secondaryAction[[#]]] &];
    nontrivialDescendant = AnyTrue[secondaryAction,
      TrueQ[FullSimplify[# != 0]] &];
    descendantZero = TrueQ[FullSimplify[descendant == 0]];
    rejectionCertified = TrueQ[descendantZero || !spatialGradient];
    reason = Which[descendantZero, "DESCENDANT_ZERO", !spatialGradient,
      "SECONDARY_ACTION_NON_GRADIENT", True, "OTHER_GAUSS_CRITERION"];
    item = <|"Primary" -> primary, "Descendant" -> descendant,
      "PrimaryAction" -> primaryAction, "SecondaryAction" -> secondaryAction,
      "SpatialGradient" -> spatialGradient, "DescendantZero" -> descendantZero,
      "SecondaryNonGradient" -> !spatialGradient,
      "RejectionCertified" -> rejectionCertified, "Reason" -> reason|>;
    If[TrueQ[firstClassDescendant && onlyScalarEnd && spatialGradient && nontrivialDescendant],
      AppendTo[accepted, item], AppendTo[rejected, item]],
    {coefficient, primaryDirections}];

  <|"Model" -> model, "Momenta" -> momenta, "MomentumMap" -> momentumMap,
    "Hessian" -> hessian, "HessianRank" -> hessianRank,
    "HessianNullspace" -> hessianNulls, "Primaries" -> primaries,
    "Hamiltonian" -> canonicalHamiltonian, "Pairs" -> pairs,
    "Constraints" -> constraints, "Stages" -> stages,
    "PrimaryCount" -> primaryCount, "Closed" -> finalState["Closed"],
    "BracketMatrix" -> bracketMatrix, "BracketRank" -> bracketRank,
    "FirstClass" -> Length[bracketKernel], "SecondClass" -> bracketRank,
    "Classes" -> classes, "PrimaryFirstClass" -> Length[primaryDirections],
    "GaussCandidates" -> Length[accepted],
    "AdditionalG" -> (Length[accepted] > 0), "Candidates" -> accepted,
    "Rejected" -> rejected|>
];

nativeLagrangian[theory_String, substitutions_: {}] := Module[
  {pFields, uFields, scalar, lambda, sigma, radial, auxiliary, coordinates,
   velocities, pVelocity, uVelocity, gt, gs, gd, gb, rhoU, alphaP, alphaU,
   alphaB, transverseMetric, kineticDensity, potentialDensity, lagrangian,
   couplings},
  pFields = Array[Symbol["p" <> theory <> ToString[#]] &, 3];
  uFields = Array[Symbol["u" <> theory <> ToString[#]] &, 3];
  scalar = Symbol["s" <> theory];
  lambda = Symbol["lambda" <> theory];
  sigma = Symbol["sigma" <> theory];
  radial = Symbol["b" <> theory];
  auxiliary = Symbol["xi" <> theory];
  coordinates = Join[pFields, uFields, {scalar, lambda, sigma, radial},
    If[theory === "C", {auxiliary}, {}]];
  velocities = (Symbol["vel" <> SymbolName[#]] &) /@ coordinates;
  pVelocity = velocities[[1 ;; 3]];
  uVelocity = velocities[[4 ;; 6]];
  gt = Symbol["gt" <> theory];
  gs = Symbol["gs" <> theory];
  gd = Symbol["gd" <> theory];
  gb = Symbol["gb" <> theory];
  rhoU = Symbol["rhoU" <> theory];
  alphaP = Symbol["alphaP" <> theory];
  alphaU = Symbol["alphaU" <> theory];
  alphaB = Symbol["alphaB" <> theory];
  transverseMetric = kSquared IdentityMatrix[3] - Outer[Times, kVector, kVector];
  kineticDensity = (pVelocity.pVelocity + velocities[[7]]^2 +
      rhoU uVelocity.uVelocity)/2 + gt pVelocity.uVelocity +
    If[theory === "A", velocities[[10]]^2/2, 0];
  potentialDensity = alphaP pFields.pFields/2 +
    alphaU uFields.transverseMetric.uFields/2 +
    gs kSquared pFields.uFields + gd (kVector.pFields)^2/2 +
    alphaB radial^2/2 + If[theory === "A", gb radial kVector.pFields, 0];
  lagrangian = kineticDensity - potentialDensity + lambda scalar +
    sigma kVector.uFields + If[theory === "C", auxiliary radial, 0];
  couplings = Join[{gt, gs, gd}, If[theory === "A", {gb}, {}]];
  If[substitutions =!= {},
    lagrangian = Factor[lagrangian /. substitutions];
    couplings = Select[couplings, !MemberQ[First /@ substitutions, #] &]
  ];
  <|"Key" -> "native_" <> theory, "Theory" -> theory,
    "Coordinates" -> coordinates, "Velocities" -> velocities,
    "Lagrangian" -> lagrangian, "Couplings" -> couplings,
    "ScalarEnd" -> {}, "SpatialBlocks" -> {{1, 2, 3}, {4, 5, 6}},
    "Multipliers" -> Join[{8, 9}, If[theory === "C", {11}, {}]],
    "Symbols" -> <|"gt" -> gt, "gs" -> gs, "gd" -> gd, "gb" -> gb,
      "rhoU" -> rhoU, "alphaP" -> alphaP, "alphaU" -> alphaU,
      "alphaB" -> alphaB|>|>
];

maxwellLagrangian[key_String : "maxwell", kineticA0_: False,
  imposeConservation_: False] := Module[
  {coordinates, velocities, spatial, spatialVelocity, electric, magnetic,
   lagrangian, currentData},
  coordinates = Array[Symbol[StringReplace[key, "_" -> ""] <> "A" <> ToString[# - 1]] &, 4];
  velocities = (Symbol["vel" <> SymbolName[#]] &) /@ coordinates;
  spatial = coordinates[[2 ;; 4]];
  spatialVelocity = velocities[[2 ;; 4]];
  electric = spatialVelocity - kVector coordinates[[1]];
  magnetic = kSquared IdentityMatrix[3] - Outer[Times, kVector, kVector];
  lagrangian = electric.electric/2 - spatial.magnetic.spatial/2 +
    If[kineticA0, velocities[[1]]^2/2, 0];
  currentData = None;
  If[key === "nonconserved_current",
    lagrangian += ncRho coordinates[[1]] + {ncJ1, ncJ2, ncJ3}.spatial;
    currentData = <|"Rho" -> ncRho, "RhoDot" -> ncRhoDot,
      "ImposeConservation" -> imposeConservation|>
  ];
  <|"Key" -> key, "Coordinates" -> coordinates, "Velocities" -> velocities,
    "Lagrangian" -> lagrangian, "Couplings" -> {}, "ScalarEnd" -> {1},
    "SpatialBlocks" -> {{2, 3, 4}}, "Multipliers" -> {},
    "CurrentData" -> currentData|>
];

gaugedHardLagrangian[removeHard_: False] := Module[
  {coordinates, velocities, spatial, spatialVelocity, phi, phiVelocity,
   rotatedPhi, magnetic, lagrangian},
  coordinates = {ghA0, ghA1, ghA2, ghA3, ghPhi1, ghPhi2};
  If[!removeHard, coordinates = Append[coordinates, ghLambda]];
  velocities = (Symbol["vel" <> SymbolName[#]] &) /@ coordinates;
  spatial = coordinates[[2 ;; 4]];
  spatialVelocity = velocities[[2 ;; 4]];
  phi = coordinates[[5 ;; 6]];
  phiVelocity = velocities[[5 ;; 6]];
  rotatedPhi = {-phi[[2]], phi[[1]]};
  magnetic = kSquared IdentityMatrix[3] - Outer[Times, kVector, kVector];
  lagrangian = (spatialVelocity - kVector coordinates[[1]]).
      (spatialVelocity - kVector coordinates[[1]])/2 -
    spatial.magnetic.spatial/2 +
    (phiVelocity - coordinates[[1]] rotatedPhi).
      (phiVelocity - coordinates[[1]] rotatedPhi)/2 - phi.phi/2 +
    If[removeHard, 0, ghLambda (phi.phi - 1)/2];
  <|"Key" -> "gauged_hard_unit", "Coordinates" -> coordinates,
    "Velocities" -> velocities, "Lagrangian" -> lagrangian,
    "Couplings" -> {}, "ScalarEnd" -> {1}, "SpatialBlocks" -> {{2, 3, 4}},
    "Multipliers" -> If[removeHard, {}, {7}]|>
];

bareSigmaLagrangian[removeHard_: False] := Module[
  {phi, coordinates, velocities, lagrangian},
  phi = {bsPhi1, bsPhi2, bsPhi3, bsPhi4};
  coordinates = Join[phi, If[removeHard, {}, {bsLambda}]];
  velocities = (Symbol["vel" <> SymbolName[#]] &) /@ coordinates;
  lagrangian = velocities[[1 ;; 4]].velocities[[1 ;; 4]]/2 - phi.phi/2 +
    If[removeHard, 0, bsLambda (phi.phi - 1)/2];
  <|"Key" -> "bare_sigma", "Coordinates" -> coordinates,
    "Velocities" -> velocities, "Lagrangian" -> lagrangian,
    "Couplings" -> {}, "ScalarEnd" -> {}, "SpatialBlocks" -> {},
    "Multipliers" -> If[removeHard, {}, {5}]|>
];

coulombGaugeLagrangian[] := Module[
  {coordinates, velocities, spatial, spatialVelocity, magnetic, lagrangian},
  coordinates = {gfA0, gfA1, gfA2, gfA3, gfEta, gfZeta};
  velocities = (Symbol["vel" <> SymbolName[#]] &) /@ coordinates;
  spatial = coordinates[[2 ;; 4]];
  spatialVelocity = velocities[[2 ;; 4]];
  magnetic = kSquared IdentityMatrix[3] - Outer[Times, kVector, kVector];
  lagrangian = (spatialVelocity - kVector coordinates[[1]]).
      (spatialVelocity - kVector coordinates[[1]])/2 -
    spatial.magnetic.spatial/2 + gfEta gfA0 + gfZeta kVector.spatial;
  <|"Key" -> "gauge_fixed_maxwell", "Coordinates" -> coordinates,
    "Velocities" -> velocities, "Lagrangian" -> lagrangian,
    "Couplings" -> {}, "ScalarEnd" -> {1}, "SpatialBlocks" -> {{2, 3, 4}},
    "Multipliers" -> {5, 6}|>
];

globalU1Lagrangian[addMaxwell_: False] := If[addMaxwell,
  maxwellLagrangian["global_only_ablation"],
  <|"Key" -> "global_only", "Coordinates" -> {globalX, globalY},
    "Velocities" -> {velGlobalX, velGlobalY},
    "Lagrangian" -> (velGlobalX^2 + velGlobalY^2)/2 -
      (globalX^2 + globalY^2)/2,
    "Couplings" -> {}, "ScalarEnd" -> {}, "SpatialBlocks" -> {},
    "Multipliers" -> {}|>
];

controlInput[key_String, ablate_: False] := Switch[key,
  "maxwell", maxwellLagrangian["maxwell", ablate],
  "gauged_hard_unit", gaugedHardLagrangian[ablate],
  "bare_sigma", bareSigmaLagrangian[ablate],
  "nonconserved_current", maxwellLagrangian["nonconserved_current", False, ablate],
  "gauge_fixed_maxwell", If[ablate, maxwellLagrangian["gauge_fixed_ablation"], coulombGaugeLagrangian[]],
  "global_only", globalU1Lagrangian[ablate]
];

controlResult[key_String, ablate_: False] := Module[
  {model, pipeline, currentData, primaryDirections, primary, gauss,
   preservation = 0, continuitySolutions, inconsistent = False, token},
  model = controlInput[key, ablate];
  pipeline = constraintPipeline[model];
  currentData = Lookup[model, "CurrentData", None];
  If[AssociationQ[currentData],
    primaryDirections = NullSpace[
      pipeline["BracketMatrix"][[All, 1 ;; pipeline["PrimaryCount"]]]];
    If[primaryDirections === {}, Throw["current primary absent", "pipelineFailure"]];
    primary = Factor@Sum[
      primaryDirections[[1, i]] pipeline["Constraints"][[i]],
      {i, pipeline["PrimaryCount"]}];
    gauss = Factor@poissonBracket[primary, pipeline["Hamiltonian"], pipeline["Pairs"]];
    preservation = Factor[
      D[gauss, currentData["Rho"]] currentData["RhoDot"] +
      poissonBracket[gauss, pipeline["Hamiltonian"], pipeline["Pairs"]]];
    If[TrueQ[currentData["ImposeConservation"]],
      continuitySolutions = Solve[preservation == 0, currentData["RhoDot"]];
      If[continuitySolutions === {}, Throw["continuity solve", "pipelineFailure"]];
      preservation = Factor[preservation /. First[continuitySolutions]]
    ];
    inconsistent = !TrueQ[FullSimplify[preservation == 0]]
  ];
  token = Switch[key,
    "maxwell", If[pipeline["GaussCandidates"] > 0 && pipeline["SecondClass"] == 0,
      "FIRST_CLASS_GAUSS", "NO_PRIMARY_GAUSS_CHAIN"],
    "gauged_hard_unit", If[pipeline["GaussCandidates"] > 0 &&
      pipeline["FirstClass"] > 0 && pipeline["SecondClass"] > 0, "MIXED", "NOT_MIXED"],
    "bare_sigma", If[pipeline["FirstClass"] == 0 && pipeline["SecondClass"] > 0 &&
      pipeline["GaussCandidates"] == 0, "SECOND_CLASS_RADIAL_NO_GAUSS", "BAD_SIGMA_CLASS"],
    "nonconserved_current", If[pipeline["GaussCandidates"] > 0 && inconsistent,
      "INCONSISTENT_PRESERVATION", "CURRENT_CONSISTENT"],
    "gauge_fixed_maxwell", If[pipeline["FirstClass"] == 0 && pipeline["SecondClass"] > 0 &&
      pipeline["GaussCandidates"] == 0, "SECOND_CLASS_NO_LOCAL_GAUGE", "LOCAL_GAUGE_REMAINS"],
    "global_only", If[Length[pipeline["HessianNullspace"]] == 0 &&
      pipeline["FirstClass"] == 0 && pipeline["GaussCandidates"] == 0,
      "GLOBAL_CHARGE_NO_LOCAL_GAUSS", "LOCAL_GAUGE_PRESENT"]
  ];
  <|"Classification" -> token, "FC" -> pipeline["FirstClass"],
    "SC" -> pipeline["SecondClass"], "G" -> pipeline["GaussCandidates"],
    "Nullity" -> Length[pipeline["HessianNullspace"]],
    "Preservation" -> preservation, "Pipeline" -> pipeline|>
];

couplingLocations[pipeline_Association] := Association@Table[
  coupling -> DeleteCases[{
    If[!FreeQ[pipeline["MomentumMap"], coupling], "momentum_map", Nothing],
    If[!FreeQ[pipeline["Hessian"], coupling], "hessian", Nothing],
    If[!FreeQ[pipeline["Constraints"], coupling], "constraints", Nothing],
    If[!FreeQ[pipeline["BracketMatrix"], coupling], "bracket_matrix", Nothing]
  }, Nothing],
  {coupling, pipeline["Model"]["Couplings"]}];

commonNullCertificate[theory_String] := Module[
  {model, symbols, qPhysical, potential, potentialMatrix, x, y,
   transverseVector, results = <||>, nullVector, residual,
   coefficientEquations, solveRules, quantifiedReduction, solveLogic},
  model = nativeLagrangian[theory];
  symbols = model["Symbols"];
  qPhysical = model["Coordinates"][[1 ;; 6]];
  potential = -model["Lagrangian"] /.
    Thread[model["Velocities"] -> ConstantArray[0, Length[model["Velocities"]]]];
  potentialMatrix = Table[D[potential, qPhysical[[i]], qPhysical[[j]]],
    {i, 6}, {j, 6}];
  transverseVector = {x, y, -(x + 2 y)/3};
  Do[
    nullVector = Join[-sign transverseVector, transverseVector];
    residual = Factor /@ (potentialMatrix.nullVector);
    coefficientEquations = DeleteDuplicates@DeleteCases[
      FullSimplify /@ Flatten@Table[{
        Coefficient[residual[[i]], x], Coefficient[residual[[i]], y],
        residual[[i]] /. {x -> 0, y -> 0}}, {i, Length[residual]}], 0];
    solveRules = Solve[Thread[coefficientEquations == 0],
      {symbols["alphaP"], symbols["gs"]}, Reals];
    quantifiedReduction = Reduce[
      ForAll[{x, y}, And @@ Thread[residual == 0]],
      {symbols["alphaP"], symbols["gs"]}, Reals];
    solveLogic = Or @@ Table[
      And @@ Thread[{symbols["alphaP"], symbols["gs"]} ==
        ({symbols["alphaP"], symbols["gs"]} /. rule)],
      {rule, solveRules}];
    AssociateTo[results, sign -> <|"Residual" -> residual,
      "SolveRules" -> solveRules, "ReduceResult" -> quantifiedReduction,
      "SolveLogic" -> solveLogic, "TransverseParameters" -> {x, y}|>],
    {sign, {-1, 1}}];
  results
];

scanNativePoint[theory_String, rules_List] := Module[{pipeline},
  pipeline = constraintPipeline[nativeLagrangian[theory, rules]];
  <|"FC" -> pipeline["FirstClass"], "G" -> pipeline["GaussCandidates"],
    "Additional" -> pipeline["AdditionalG"],
    "Nullity" -> Length[pipeline["HessianNullspace"]],
    "Pipeline" -> pipeline, "Rules" -> rules|>
];

nativeFamilyAnalysis[theory_String] := Module[
  {model, symbols, openPipeline, physicalHessian, determinant,
   componentDeterminant, degeneracyReduction, common, boundary = {},
   hardening = {}, sign, commonRules, rankDropRule, noncommonRules,
   point, seed, sweep, alphaUValue, alphaPValue, sampledRules,
   couplingMap, row, label, rules},
  Print["PROGRESS MATHEMATICA THEORY-", theory, " OPEN PIPELINE"];
  model = nativeLagrangian[theory];
  symbols = model["Symbols"];
  openPipeline = constraintPipeline[model];
  physicalHessian = openPipeline["Hessian"][[1 ;; 6, 1 ;; 6]];
  determinant = Factor[Det[physicalHessian]];
  componentDeterminant = Factor[Det[physicalHessian[[{1, 4}, {1, 4}]]]];
  degeneracyReduction = Reduce[
    determinant == 0 && symbols["rhoU"] > 0, symbols["gt"], Reals];
  Print["PROGRESS MATHEMATICA THEORY-", theory, " Solve/Reduce COMMON NULL"];
  common = commonNullCertificate[theory];
  Do[
    commonRules = Join[{
      symbols["gt"] -> sign, symbols["rhoU"] -> 1},
      First[common[sign]["SolveRules"]]];
    rankDropRule = First@Solve[
      symbols["alphaP"] + 14 symbols["alphaU"] - 28 sign symbols["gs"] == 0,
      symbols["gs"]];
    noncommonRules = {
      symbols["gt"] -> sign, symbols["rhoU"] -> 1,
      symbols["alphaU"] -> 1, symbols["alphaP"] -> 1,
      rankDropRule[[1]] /. {symbols["alphaU"] -> 1, symbols["alphaP"] -> 1},
      symbols["gd"] -> 0};
    If[theory === "A", noncommonRules = Append[noncommonRules, symbols["gb"] -> 0]];
    Do[
      label = row[[1]];
      rules = row[[2]];
      Print["PROGRESS MATHEMATICA THEORY-", theory, " BOUNDARY sign=", sign,
        " lane=", label];
      point = scanNativePoint[theory, rules];
      point = Join[point, <|"Sign" -> sign, "Label" -> label|>];
      AppendTo[boundary, point];
      If[label === "common",
        AppendTo[hardening, <|"Sign" -> sign, "FC" -> point["FC"],
          "PrimaryFC" -> point["Pipeline"]["PrimaryFirstClass"],
          "Rejected" -> point["Pipeline"]["Rejected"],
          "Candidates" -> point["Pipeline"]["Candidates"]|>]
      ],
      {row, {
        {"generic", {symbols["gt"] -> sign, symbols["rhoU"] -> 1}},
        {"rankdrop", noncommonRules},
        {"common", commonRules}
      }}],
    {sign, {-1, 1}}];

  seed = If[theory === "A", 260713, 260715];
  Print["PROGRESS MATHEMATICA THEORY-", theory, " BlockRandom SWEEP seed=", seed];
  sweep = BlockRandom[
    SeedRandom[seed];
    Flatten[Table[
      Table[
        alphaUValue = RandomInteger[{1, 5}];
        alphaPValue = RandomInteger[{1, 40}];
        While[alphaPValue == 14 alphaUValue,
          alphaPValue = RandomInteger[{1, 40}]];
        rankDropRule = First@Solve[
          symbols["alphaP"] + 14 symbols["alphaU"] - 28 sign symbols["gs"] == 0,
          symbols["gs"]];
        sampledRules = {
          symbols["gt"] -> sign, symbols["rhoU"] -> 1,
          symbols["alphaU"] -> alphaUValue,
          symbols["alphaP"] -> alphaPValue,
          rankDropRule[[1]] /. {symbols["alphaU"] -> alphaUValue,
            symbols["alphaP"] -> alphaPValue}, symbols["gd"] -> 0};
        If[theory === "A", sampledRules = Append[sampledRules, symbols["gb"] -> 0]];
        point = scanNativePoint[theory, sampledRules];
        Join[point, <|"Sign" -> sign, "Sample" -> sample|>],
        {sample, 3}],
      {sign, {-1, 1}}], 1]
  ];
  couplingMap = couplingLocations[openPipeline];
  <|"Theory" -> theory, "Model" -> model, "Symbols" -> symbols,
    "Open" -> openPipeline, "Determinant" -> determinant,
    "ComponentDeterminant" -> componentDeterminant,
    "DegeneracyReduce" -> degeneracyReduction, "Common" -> common,
    "Boundary" -> boundary, "Sweep" -> sweep, "Seed" -> seed,
    "Hardening" -> hardening, "CouplingLocations" -> couplingMap|>
];

pipelineSignature[pipeline_Association] := {
  pipeline["HessianRank"], Length[pipeline["HessianNullspace"]],
  Length[pipeline["Constraints"]], pipeline["Stages"], pipeline["Classes"],
  pipeline["BracketRank"], pipeline["FirstClass"], pipeline["SecondClass"],
  pipeline["GaussCandidates"], pipeline["AdditionalG"],
  pipeline["FirstClass"] == Length[pipeline["Constraints"]] - pipeline["BracketRank"]
};

expectedFamilySignature["A"] = {8, 2, 8, {0, 0, 1, 1, 2, 2, 3, 3},
  ConstantArray["SECOND_CLASS_COMPONENT", 8], 8, 0, 8, 0, False, True};
expectedFamilySignature["C"] = {7, 4, 12,
  {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3},
  ConstantArray["SECOND_CLASS_COMPONENT", 12], 12, 0, 12, 0, False, True};

verdictFromComputed[openSearches_List, tunedSearches_List] := Which[
  AnyTrue[openSearches, TrueQ[#] &], "FIRST_CLASS_GENERIC_EM_CANDIDATE",
  AnyTrue[tunedSearches, TrueQ[#] &], "FIRST_CLASS_TUNED_INVERSE_DESIGN",
  True, verdictToken
];

sourceToothIDs = {
  "Q1_LAGRANGIAN_LEGENDRE", "Q2_MAXWELL_COMPUTED", "Q3_SIX_CONTROLS",
  "Q4_WOLFRAM_INDEPENDENCE", "Q5_G_SEARCH", "Q6_DIRAC_CLOSURE",
  "GUARD_COUPLINGS_ENTER_A", "GUARD_COUPLINGS_ENTER_C", "GUARD_SEARCH_CAPABLE",
  "HARDENING_DESCENDANT_A", "HARDENING_DESCENDANT_C",
  "THEORY_A_HESSIAN_RANK_NULLITY", "THEORY_A_CONSTRAINT_COUNT",
  "THEORY_A_STAGES", "THEORY_A_ALL_SECOND_CLASS", "THEORY_A_PB_FC_SC",
  "THEORY_A_G_SEARCH", "THEORY_C_HESSIAN_RANK_NULLITY",
  "THEORY_C_CONSTRAINT_COUNT", "THEORY_C_STAGES", "THEORY_C_ALL_SECOND_CLASS",
  "THEORY_C_PB_FC_SC", "THEORY_C_G_SEARCH", "KINETIC_HESSIAN_DETERMINANT",
  "CONTROL_MAXWELL", "CONTROL_GAUGED_HARD_UNIT", "CONTROL_BARE_SIGMA",
  "CONTROL_NONCONSERVED_CURRENT", "CONTROL_COULOMB_GAUGE", "CONTROL_GLOBAL_U1",
  "ABLATION_MAXWELL", "ABLATION_GAUGED_HARD_UNIT", "ABLATION_BARE_SIGMA",
  "ABLATION_NONCONSERVED_CURRENT", "ABLATION_COULOMB_GAUGE", "ABLATION_GLOBAL_U1",
  "BOUNDARY_SCAN_A", "BOUNDARY_SCAN_C", "RANDOM_SWEEP_A", "RANDOM_SWEEP_C",
  "SEEDS_AND_CARDINALITY", "COMMON_NULL_A", "COMMON_NULL_C",
  "SOURCE_FIRST_ORDERING", "SHEAR_DUPLICATE", "DECISION_ORDER_BRANCH2",
  "HONEST_TUNED_SCOPE", "VERDICT_TOTALITY", "ARGPARSE_OUT_DIR_HARNESS",
  "JSON_ARTIFACT_WRITING", "CROSS_ENGINE_FILE_COMPARATOR"
};

sourceManifest = {
  {"Q1_LAGRANGIAN_LEGENDRE", "PRESERVED", "STAGE033_SHARED_PIPELINE"},
  {"Q2_MAXWELL_COMPUTED", "PRESERVED", "STAGE033_CONTROL_MAXWELL"},
  {"Q3_SIX_CONTROLS", "REPLACED_BY_STRONGER", "STAGE033_NATIVE_ABLATIONS"},
  {"Q4_WOLFRAM_INDEPENDENCE", "REPLACED_BY_STRONGER", "STAGE033_REAUTHORED_WOLFRAM"},
  {"Q5_G_SEARCH", "PRESERVED", "STAGE033_GAUSS_SEARCH"},
  {"Q6_DIRAC_CLOSURE", "PRESERVED", "STAGE033_DIRAC_CLOSURE"},
  {"GUARD_COUPLINGS_ENTER_A", "PRESERVED", "STAGE033_COUPLING_GUARD_A"},
  {"GUARD_COUPLINGS_ENTER_C", "PRESERVED", "STAGE033_COUPLING_GUARD_C"},
  {"GUARD_SEARCH_CAPABLE", "PRESERVED", "STAGE033_SEARCH_CAPABILITY"},
  {"HARDENING_DESCENDANT_A", "PRESERVED", "STAGE033_DESCENDANT_A"},
  {"HARDENING_DESCENDANT_C", "PRESERVED", "STAGE033_DESCENDANT_C"},
  {"THEORY_A_HESSIAN_RANK_NULLITY", "PRESERVED", "STAGE033_SIGNATURE_A"},
  {"THEORY_A_CONSTRAINT_COUNT", "PRESERVED", "STAGE033_SIGNATURE_A"},
  {"THEORY_A_STAGES", "PRESERVED", "STAGE033_SIGNATURE_A"},
  {"THEORY_A_ALL_SECOND_CLASS", "PRESERVED", "STAGE033_SIGNATURE_A"},
  {"THEORY_A_PB_FC_SC", "REPLACED_BY_STRONGER", "STAGE033_RANK_NULLITY_A"},
  {"THEORY_A_G_SEARCH", "PRESERVED", "STAGE033_GAUSS_ABSENCE_A"},
  {"THEORY_C_HESSIAN_RANK_NULLITY", "PRESERVED", "STAGE033_SIGNATURE_C"},
  {"THEORY_C_CONSTRAINT_COUNT", "PRESERVED", "STAGE033_SIGNATURE_C"},
  {"THEORY_C_STAGES", "PRESERVED", "STAGE033_SIGNATURE_C"},
  {"THEORY_C_ALL_SECOND_CLASS", "PRESERVED", "STAGE033_SIGNATURE_C"},
  {"THEORY_C_PB_FC_SC", "REPLACED_BY_STRONGER", "STAGE033_RANK_NULLITY_C"},
  {"THEORY_C_G_SEARCH", "PRESERVED", "STAGE033_GAUSS_ABSENCE_C"},
  {"KINETIC_HESSIAN_DETERMINANT", "REPLACED_BY_STRONGER", "STAGE033_FULL_AND_COMPONENT_DET"},
  {"CONTROL_MAXWELL", "PRESERVED", "STAGE033_CONTROL_MAXWELL"},
  {"CONTROL_GAUGED_HARD_UNIT", "PRESERVED", "STAGE033_CONTROL_GAUGED"},
  {"CONTROL_BARE_SIGMA", "PRESERVED", "STAGE033_CONTROL_SIGMA"},
  {"CONTROL_NONCONSERVED_CURRENT", "PRESERVED", "STAGE033_CONTROL_CURRENT"},
  {"CONTROL_COULOMB_GAUGE", "PRESERVED", "STAGE033_CONTROL_COULOMB"},
  {"CONTROL_GLOBAL_U1", "PRESERVED", "STAGE033_CONTROL_GLOBAL"},
  {"ABLATION_MAXWELL", "REPLACED_BY_STRONGER", "STAGE033_ENV_OWN_ASSERT"},
  {"ABLATION_GAUGED_HARD_UNIT", "REPLACED_BY_STRONGER", "STAGE033_ENV_OWN_ASSERT"},
  {"ABLATION_BARE_SIGMA", "REPLACED_BY_STRONGER", "STAGE033_ENV_OWN_ASSERT"},
  {"ABLATION_NONCONSERVED_CURRENT", "REPLACED_BY_STRONGER", "STAGE033_ENV_OWN_ASSERT"},
  {"ABLATION_COULOMB_GAUGE", "REPLACED_BY_STRONGER", "STAGE033_ENV_OWN_ASSERT"},
  {"ABLATION_GLOBAL_U1", "REPLACED_BY_STRONGER", "STAGE033_ENV_OWN_ASSERT"},
  {"BOUNDARY_SCAN_A", "PRESERVED", "STAGE033_BOUNDARY_A"},
  {"BOUNDARY_SCAN_C", "PRESERVED", "STAGE033_BOUNDARY_C"},
  {"RANDOM_SWEEP_A", "PRESERVED", "STAGE033_SWEEP_A"},
  {"RANDOM_SWEEP_C", "PRESERVED", "STAGE033_SWEEP_C"},
  {"SEEDS_AND_CARDINALITY", "PRESERVED", "STAGE033_FIXED_SEEDS"},
  {"COMMON_NULL_A", "REPLACED_BY_STRONGER", "STAGE033_SYMBOLIC_SOLVE_A"},
  {"COMMON_NULL_C", "REPLACED_BY_STRONGER", "STAGE033_SYMBOLIC_SOLVE_C"},
  {"SOURCE_FIRST_ORDERING", "PRESERVED", "STAGE033_SOURCE_FIRST"},
  {"SHEAR_DUPLICATE", "PRESERVED", "STAGE033_NO_G_BOOKKEEPING"},
  {"DECISION_ORDER_BRANCH2", "PRESERVED", "STAGE033_QUADRATIC_BRANCH2"},
  {"HONEST_TUNED_SCOPE", "REPLACED_BY_STRONGER", "STAGE033_ARGUED_SCANNED_SCOPE"},
  {"VERDICT_TOTALITY", "REPLACED_BY_STRONGER", "STAGE033_COMPUTED_VERDICT"},
  {"ARGPARSE_OUT_DIR_HARNESS", "SCOPED_OUT", "STAGE033_PRINT_ONLY_CONTRACT"},
  {"JSON_ARTIFACT_WRITING", "SCOPED_OUT", "STAGE033_ZERO_FILE_IO"},
  {"CROSS_ENGINE_FILE_COMPARATOR", "SCOPED_OUT", "STAGE033_INDEPENDENT_TOKENS"}
};

expectedManifestCounts = KeySort@<|
  "PRESERVED" -> 33, "REPLACED_BY_STRONGER" -> 15, "SCOPED_OUT" -> 3|>;
committedManifestDigest = "6b191e77fefe24c9000445f01e4e2c6154ab1bb9b15bb40a6d1515dc463a7e9d";

lexicographicCodeLess[left_List, right_List] := Module[{limit, index},
  limit = Min[Length[left], Length[right]];
  index = SelectFirst[Range[limit], left[[#]] =!= right[[#]] &, Missing["NotFound"]];
  If[MissingQ[index], Length[left] < Length[right], left[[index]] < right[[index]]]
];

canonicalManifestText[manifest_] := Module[{rows},
  rows = StringRiffle[#, "|"] & /@ manifest;
  StringRiffle[Sort[rows,
    lexicographicCodeLess[ToCharacterCode[#1], ToCharacterCode[#2]] &], "\n"]
];

manifestSHA256[manifest_] :=
  IntegerString[Hash[canonicalManifestText[manifest], "SHA256"], 16, 64];

ok = Catch[
  If[activeMutation =!= "" && !MemberQ[toothOrder, activeMutation],
    Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
    Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
    raise["UNKNOWN_MUTATION"]
  ];

  Print["ledger_stage033_native_p_no_emergent_gauss Mathematica audit"];
  Print["PIPELINE=independent Lagrangian -> native fixed-point constraint descriptor -> Gauss classification"];
  Print["COMMON_NULL_ROUTE=Reduce[ForAll[...]] + Solve coefficient system (no inserted tuning)"];
  If[activeMutation =!= "",
    Print["ACTIVE_MUTATION=", activeMutation];
    Print["MUTATED_PRIMITIVE=", ablationDescriptions[activeMutation]]
  ];

  heading["Computed input-Lagrangian pipeline and six controls"];
  controlKeys = {"maxwell", "gauged_hard_unit", "bare_sigma",
    "nonconserved_current", "gauge_fixed_maxwell", "global_only"};
  baselineControls = Association@Table[key -> controlResult[key, False], {key, controlKeys}];
  ablatedControls = Association@Table[key -> controlResult[key, True], {key, controlKeys}];
  maxwellLive = baselineControls["maxwell"];

  theoryA = nativeFamilyAnalysis["A"];
  Print["PROGRESS MATHEMATICA THEORY-A OPEN+BOUNDARY+SWEEP COMPLETE"];
  theoryC = nativeFamilyAnalysis["C"];
  Print["PROGRESS MATHEMATICA THEORY-C OPEN+BOUNDARY+SWEEP COMPLETE"];

  legendreTest = theoryA["Open"]["Hamiltonian"];
  If[activeMutation === "PASS_Q1_LAGRANGIAN_LEGENDRE",
    legendreTest += First[theoryA["Model"]["Velocities"]]];
  expectBool["PASS_Q1_LAGRANGIAN_LEGENDRE",
    Length[theoryA["Open"]["Primaries"]] == Length[theoryA["Open"]["HessianNullspace"]] &&
    Length[theoryC["Open"]["Primaries"]] == Length[theoryC["Open"]["HessianNullspace"]] &&
    FreeQ[legendreTest, Alternatives @@ theoryA["Model"]["Velocities"]]];

  q2Result = controlResult["maxwell", activeMutation === "PASS_Q2_MAXWELL_COMPUTED"];
  q2Tuple = {q2Result["Classification"], q2Result["FC"], q2Result["SC"],
    q2Result["G"], q2Result["Nullity"]};
  expectBool["PASS_Q2_MAXWELL_COMPUTED",
    q2Tuple === {"FIRST_CLASS_GAUSS", 2, 0, 1, 1} &&
    q2Result["Pipeline"]["BracketMatrix"] === ConstantArray[0, {2, 2}], q2Tuple];

  controlTable = baselineControls;
  If[activeMutation === "PASS_Q3_SIX_CONTROLS_SHARED_PIPELINE",
    controlTable = KeyDrop[controlTable, "global_only"]];
  expectBool["PASS_Q3_SIX_CONTROLS_SHARED_PIPELINE",
    Keys[controlTable] === controlKeys &&
    AllTrue[Values[controlTable], AssociationQ[#["Pipeline"]] &], Keys[controlTable]];

  q4OK = True;
  Do[
    item = family["Common"][sign];
    rulesTest = item["SolveRules"];
    If[activeMutation === "PASS_Q4_CAS_COMMON_NULL_RECONSTRUCTION" &&
      family["Theory"] === "A" && sign == 1,
      rulesTest = {Thread[
        {family["Symbols"]["alphaP"], family["Symbols"]["gs"]} ->
        {family["Symbols"]["alphaP"] /. First[rulesTest],
          -(family["Symbols"]["gs"] /. First[rulesTest])}]}];
    q4OK = q4OK && Length[rulesTest] == 1 &&
      And @@ (TrueQ[FullSimplify[# == 0]] & /@
        (item["Residual"] /. First[rulesTest])) &&
      item["ReduceResult"] =!= False &&
      TrueQ[FullSimplify[item["ReduceResult"] /. First[item["SolveRules"]],
        Element[family["Symbols"]["alphaU"], Reals]]],
    {family, {theoryA, theoryC}}, {sign, {-1, 1}}];
  expectBool["PASS_Q4_CAS_COMMON_NULL_RECONSTRUCTION", q4OK];

  q5Result = controlResult["maxwell", activeMutation === "PASS_Q5_GAUSS_SEARCH_LIVE"];
  expectBool["PASS_Q5_GAUSS_SEARCH_LIVE",
    q5Result["Pipeline"]["GaussCandidates"] > 0 &&
    TrueQ[q5Result["Pipeline"]["AdditionalG"]], q5Result["Classification"]];

  q6Closed = TrueQ[theoryA["Open"]["Closed"]] && TrueQ[theoryC["Open"]["Closed"]];
  If[activeMutation === "PASS_Q6_DIRAC_CLOSURE",
    q6Closed = TrueQ[constraintPipeline[theoryA["Model"], 1]["Closed"]]];
  expectBool["PASS_Q6_DIRAC_CLOSURE", q6Closed];

  heading["Native THEORY-A/C symbolic open stratum"];
  Do[
    family = row[[1]];
    couplingTooth = row[[2]];
    locationsTest = family["CouplingLocations"];
    If[activeMutation === couplingTooth,
      targetCoupling = If[family["Theory"] === "A",
        family["Symbols"]["gb"], family["Symbols"]["gd"]];
      badModel = Join[family["Model"],
        <|"Lagrangian" -> (family["Model"]["Lagrangian"] /. targetCoupling -> 0)|>];
      locationsTest = couplingLocations[constraintPipeline[badModel]]
    ];
    expectedCouplings = If[family["Theory"] === "A", {gtA, gsA, gdA, gbA},
      {gtC, gsC, gdC}];
    expectBool[couplingTooth,
      Sort[Keys[locationsTest]] === Sort[expectedCouplings] &&
      AllTrue[Values[locationsTest], Length[#] > 0 &], locationsTest],
    {row, {
      {theoryA, "PASS_GUARD_COUPLINGS_ENTER_A"},
      {theoryC, "PASS_GUARD_COUPLINGS_ENTER_C"}
    }}];

  Do[
    family = row[[1]];
    signatureTooth = row[[2]];
    signatureTest = pipelineSignature[family["Open"]];
    If[activeMutation === signatureTooth,
      badModel = Join[family["Model"], <|"Lagrangian" ->
        family["Model"]["Lagrangian"] + family["Model"]["Velocities"][[8]]^2/2|>];
      signatureTest = pipelineSignature[constraintPipeline[badModel]]
    ];
    expectBool[signatureTooth,
      signatureTest === expectedFamilySignature[family["Theory"]], signatureTest],
    {row, {
      {theoryA, "PASS_THEORY_A_SIGNATURE"},
      {theoryC, "PASS_THEORY_C_SIGNATURE"}
    }}];

  determinantA = {theoryA["Determinant"], theoryA["ComponentDeterminant"]};
  determinantC = {theoryC["Determinant"], theoryC["ComponentDeterminant"]};
  If[activeMutation === "PASS_KINETIC_HESSIAN_DETERMINANT",
    badModel = Join[theoryA["Model"], <|"Lagrangian" ->
      theoryA["Model"]["Lagrangian"] + theoryA["Model"]["Velocities"][[4]]^2/2|>];
    badPipeline = constraintPipeline[badModel];
    determinantA = {
      Factor@Det[badPipeline["Hessian"][[1 ;; 6, 1 ;; 6]]],
      Factor@Det[badPipeline["Hessian"][[{1, 4}, {1, 4}]]]
    }
  ];
  expectBool["PASS_KINETIC_HESSIAN_DETERMINANT",
    TrueQ[FullSimplify[determinantA[[1]] ==
      (theoryA["Symbols"]["rhoU"] - theoryA["Symbols"]["gt"]^2)^3]] &&
    TrueQ[FullSimplify[determinantA[[2]] ==
      theoryA["Symbols"]["rhoU"] - theoryA["Symbols"]["gt"]^2]] &&
    TrueQ[FullSimplify[determinantC[[1]] ==
      (theoryC["Symbols"]["rhoU"] - theoryC["Symbols"]["gt"]^2)^3]] &&
    TrueQ[FullSimplify[determinantC[[2]] ==
      theoryC["Symbols"]["rhoU"] - theoryC["Symbols"]["gt"]^2]] &&
    theoryA["DegeneracyReduce"] =!= False && theoryC["DegeneracyReduce"] =!= False &&
    TrueQ[Resolve[ForAll[{gtA, rhoUA}, Implies[rhoUA > 0,
      Equivalent[determinantA[[1]] == 0, gtA^2 == rhoUA]]], Reals]] &&
    TrueQ[Resolve[ForAll[{gtC, rhoUC}, Implies[rhoUC > 0,
      Equivalent[determinantC[[1]] == 0, gtC^2 == rhoUC]]], Reals]],
    <|"A" -> determinantA, "C" -> determinantC|>];

  heading["Solve/Reduce common-null locus and tuned scans"];
  Do[
    family = row[[1]];
    commonTooth = row[[2]];
    commonOK = True;
    Do[
      item = family["Common"][sign];
      rulesTest = item["SolveRules"];
      If[activeMutation === commonTooth && sign == 1,
        rulesTest = {Thread[
          {family["Symbols"]["alphaP"], family["Symbols"]["gs"]} ->
          {family["Symbols"]["alphaP"] /. First[rulesTest],
            -(family["Symbols"]["gs"] /. First[rulesTest])}]}];
      commonOK = commonOK && Length[rulesTest] == 1 &&
        TrueQ[FullSimplify[(family["Symbols"]["alphaP"] /. First[rulesTest]) ==
          14 family["Symbols"]["alphaU"]]] &&
        TrueQ[FullSimplify[(family["Symbols"]["gs"] /. First[rulesTest]) ==
          sign family["Symbols"]["alphaU"]]] &&
        And @@ (TrueQ[FullSimplify[# == 0]] & /@
          (item["Residual"] /. First[rulesTest])),
      {sign, {-1, 1}}];
    expectBool[commonTooth, commonOK],
    {row, {
      {theoryA, "PASS_COMMON_NULL_A"}, {theoryC, "PASS_COMMON_NULL_C"}
    }}];

  Do[
    family = row[[1]];
    boundaryTooth = row[[2]];
    pointsTest = family["Boundary"];
    If[activeMutation === boundaryTooth,
      badRules = {family["Symbols"]["gt"] -> 1, family["Symbols"]["rhoU"] -> 1,
        family["Symbols"]["alphaP"] -> 14, family["Symbols"]["alphaU"] -> 1,
        family["Symbols"]["gs"] -> 3/2, family["Symbols"]["gd"] -> 0};
      If[family["Theory"] === "A", badRules = Append[badRules, family["Symbols"]["gb"] -> 0]];
      pointsTest = ReplacePart[pointsTest, 3 -> scanNativePoint[family["Theory"], badRules]]
    ];
    observed = ({#["FC"], #["G"], #["Additional"], #["Nullity"]} &) /@ pointsTest;
    nullityTarget = If[family["Theory"] === "A", 5, 7];
    expected = Join @@ ConstantArray[{{0, 0, False, nullityTarget},
      {0, 0, False, nullityTarget}, {2, 0, False, nullityTarget}}, 2];
    expectBool[boundaryTooth, observed === expected, observed],
    {row, {
      {theoryA, "PASS_BOUNDARY_SCAN_A"}, {theoryC, "PASS_BOUNDARY_SCAN_C"}
    }}];

  Do[
    family = row[[1]];
    sweepTooth = row[[2]];
    sweepTest = family["Sweep"];
    If[activeMutation === sweepTooth,
      sweepTest = ReplacePart[sweepTest, 1 -> family["Boundary"][[3]]]];
    observed = ({#["FC"], #["G"], #["Additional"], #["Nullity"]} &) /@ sweepTest;
    signCounts = Counts[(#["Sign"] &) /@ family["Sweep"]];
    nullityTarget = If[family["Theory"] === "A", 5, 7];
    expectBool[sweepTooth,
      Length[observed] == 6 && signCounts === <|-1 -> 3, 1 -> 3|> &&
      observed === ConstantArray[{0, 0, False, nullityTarget}, 6] &&
      family["Seed"] == If[family["Theory"] === "A", 260713, 260715],
      <|"Seed" -> family["Seed"], "Signs" -> signCounts, "Observed" -> observed|>],
    {row, {
      {theoryA, "PASS_RANDOMIZED_SWEEP_A"}, {theoryC, "PASS_RANDOMIZED_SWEEP_C"}
    }}];

  maxwellCandidate = First[maxwellLive["Pipeline"]["Candidates"]];
  Do[
    family = row[[1]];
    descendantTooth = row[[2]];
    accountingTooth = row[[3]];
    recordsTest = Map[Join[#, <|"Rejected" -> #["Rejected"]|>] &, family["Hardening"]];
    If[activeMutation === descendantTooth,
      recordsTest = ReplacePart[recordsTest, {1, "Rejected"} ->
        Append[recordsTest[[1, "Rejected"]], maxwellCandidate]];
      recordsTest = ReplacePart[recordsTest, {1, "PrimaryFC"} ->
        (recordsTest[[1, "PrimaryFC"]] + 1)]
    ];
    If[activeMutation === accountingTooth,
      recordsTest = ReplacePart[recordsTest, {1, "PrimaryFC"} ->
        (recordsTest[[1, "PrimaryFC"]] + 1)]
    ];
    descendantOK = Length[recordsTest] == 2 &&
      Total[Length[#["Rejected"]] & /@ recordsTest] == 4 &&
      AllTrue[Flatten[(#["Rejected"] &) /@ recordsTest],
        TrueQ[#["DescendantZero"]] && TrueQ[#["SecondaryNonGradient"]] &&
        TrueQ[#["RejectionCertified"]] && #["Reason"] === "DESCENDANT_ZERO" &];
    accountingOK = AllTrue[recordsTest,
      #["PrimaryFC"] == Length[#["Rejected"]] + Length[#["Candidates"]] &];
    expectBool[descendantTooth, descendantOK,
      <|"Strata" -> Length[recordsTest],
        "Directions" -> Total[Length[#["Rejected"]] & /@ recordsTest]|>];
    expectBool[accountingTooth, accountingOK],
    {row, {
      {theoryA, "PASS_HARDENING_DESCENDANT_A", "PASS_PRIMARY_FC_ACCOUNTING_A"},
      {theoryC, "PASS_HARDENING_DESCENDANT_C", "PASS_PRIMARY_FC_ACCOUNTING_C"}
    }}];

  heading["Six same-pipeline controls and native own-assert ablations"];
  searchProbe = If[activeMutation === "PASS_GUARD_SEARCH_CAPABLE",
    controlResult["maxwell", True], baselineControls["maxwell"]];
  expectBool["PASS_GUARD_SEARCH_CAPABLE",
    searchProbe["G"] > 0 && baselineControls["gauged_hard_unit"]["G"] > 0,
    <|"Maxwell" -> searchProbe["G"],
      "Gauged" -> baselineControls["gauged_hard_unit"]["G"]|>];

  controlExpected = <|
    "maxwell" -> {"FIRST_CLASS_GAUSS", 2, 0, 1, 1},
    "gauged_hard_unit" -> {"MIXED", 2, 4, 1, 2},
    "bare_sigma" -> {"SECOND_CLASS_RADIAL_NO_GAUSS", 0, 4, 0, 1},
    "nonconserved_current" -> {"INCONSISTENT_PRESERVATION", 2, 0, 1, 1},
    "gauge_fixed_maxwell" -> {"SECOND_CLASS_NO_LOCAL_GAUGE", 0, 8, 0, 3},
    "global_only" -> {"GLOBAL_CHARGE_NO_LOCAL_GAUSS", 0, 0, 0, 0}|>;
  controlPredicates = <|
    "maxwell" -> "PASS_CONTROL_MAXWELL",
    "gauged_hard_unit" -> "PASS_CONTROL_GAUGED_HARD_UNIT",
    "bare_sigma" -> "PASS_CONTROL_BARE_SIGMA",
    "nonconserved_current" -> "PASS_CONTROL_NONCONSERVED_CURRENT",
    "gauge_fixed_maxwell" -> "PASS_CONTROL_COULOMB_GAUGE",
    "global_only" -> "PASS_CONTROL_GLOBAL_U1"|>;
  Do[
    tooth = controlPredicates[key];
    testedControl = If[activeMutation === tooth, ablatedControls[key], baselineControls[key]];
    observed = {testedControl["Classification"], testedControl["FC"], testedControl["SC"],
      testedControl["G"], testedControl["Nullity"]};
    ablatedObserved = {ablatedControls[key]["Classification"], ablatedControls[key]["FC"],
      ablatedControls[key]["SC"], ablatedControls[key]["G"], ablatedControls[key]["Nullity"]};
    extraControl = If[key === "nonconserved_current",
      TrueQ[FullSimplify[testedControl["Preservation"] ==
        -ncJ1 - 2 ncJ2 - 3 ncJ3 + ncRhoDot]], True];
    expectBool[tooth, observed === controlExpected[key] &&
      ablatedObserved =!= controlExpected[key] && extraControl,
      <|"Observed" -> observed, "Ablated" -> ablatedObserved|>];
    Print["      ", key, ": ", First[controlExpected[key]], " FC=", observed[[2]],
      " SC=", observed[[3]], " G=", observed[[4]], " NULL=", observed[[5]]];
    Print["      ABLATION ", key, ": FIRED_AT_OWN_ASSERT"],
    {key, controlKeys}];

  heading["Bookkeeping, honest scope, branch decision, and verdict"];
  anyNativeG = AnyTrue[Join[
      {theoryA["Open"]["AdditionalG"], theoryC["Open"]["AdditionalG"]},
      (#["Additional"] &) /@ Join[theoryA["Boundary"], theoryA["Sweep"],
        theoryC["Boundary"], theoryC["Sweep"]]], TrueQ];
  If[activeMutation === "PASS_SOURCE_FIRST_SHEAR_BOOKKEEPING",
    anyNativeG = TrueQ[maxwellLive["Pipeline"]["AdditionalG"]]];
  sourceFirst = <|"SearchedSourceFreeFirst" -> True,
    "jAAdded" -> anyNativeG, "Sourced" -> anyNativeG|>;
  shearDuplicate = If[anyNativeG, "REQUIRES_G", "NOT_APPLICABLE_NO_G"];
  expectBool["PASS_SOURCE_FIRST_SHEAR_BOOKKEEPING",
    sourceFirst === <|"SearchedSourceFreeFirst" -> True,
      "jAAdded" -> False, "Sourced" -> False|> &&
    shearDuplicate === "NOT_APPLICABLE_NO_G",
    <|"SourceFirst" -> sourceFirst, "Shear" -> shearDuplicate|>];

  scope = <|"OpenStratum" -> "FULLY_SYMBOLIC_ALL_RETAINED_COUPLINGS",
    "TunedStratum" -> "ARGUED_PLUS_FIXED_SEED_SCANNED",
    "BoundaryPerFamily" -> 6, "RandomTotal" -> 12,
    "ExhaustiveTuned" -> False, "MissedMeasureZero" -> "TUNED_INVERSE_DESIGN",
    "GenericNoGoDecisive" -> True|>;
  If[activeMutation === "PASS_HONEST_TUNED_SCOPE",
    scope["ExhaustiveTuned"] = True];
  expectBool["PASS_HONEST_TUNED_SCOPE",
    scope === <|"OpenStratum" -> "FULLY_SYMBOLIC_ALL_RETAINED_COUPLINGS",
      "TunedStratum" -> "ARGUED_PLUS_FIXED_SEED_SCANNED",
      "BoundaryPerFamily" -> 6, "RandomTotal" -> 12,
      "ExhaustiveTuned" -> False, "MissedMeasureZero" -> "TUNED_INVERSE_DESIGN",
      "GenericNoGoDecisive" -> True|>, scope];

  branchFC = If[activeMutation === "PASS_DECISION_ORDER_BRANCH2",
    maxwellLive["FC"], theoryA["Open"]["FirstClass"] + theoryC["Open"]["FirstClass"]];
  branch = If[branchFC == 0, "BRANCH_2_QUADRATIC_ABSENCE",
    "BRANCH_1_LINEARIZED_GAUGE_PRESENT"];
  expectBool["PASS_DECISION_ORDER_BRANCH2",
    branch === "BRANCH_2_QUADRATIC_ABSENCE", branch];

  openFlags = {theoryA["Open"]["AdditionalG"], theoryC["Open"]["AdditionalG"]};
  If[activeMutation === "PASS_VERDICT_TOTALITY",
    openFlags = {maxwellLive["Pipeline"]["AdditionalG"], theoryC["Open"]["AdditionalG"]}];
  tunedFlags = (#["Additional"] &) /@ Join[theoryA["Boundary"], theoryA["Sweep"],
    theoryC["Boundary"], theoryC["Sweep"]];
  verdictA = verdictFromComputed[openFlags, tunedFlags];
  If[activeMutation === "PASS_VERDICT_TOTALITY",
    expectBool["PASS_VERDICT_TOTALITY", verdictA === verdictToken,
      <|"Rederived" -> verdictA,
        "RequiredMutationToken" -> "FIRST_CLASS_GENERIC_EM_CANDIDATE"|>],
    verdictC = verdictFromComputed[{theoryC["Open"]["AdditionalG"]},
      (#["Additional"] &) /@ Join[theoryC["Boundary"], theoryC["Sweep"]]];
    expectBool["PASS_VERDICT_TOTALITY",
      verdictA === verdictC === verdictToken, <|"A" -> verdictA, "C" -> verdictC|>]
  ];

  heading["Canonical source-to-stage predicate manifest"];
  manifestTest = If[activeMutation === "PASS_SOURCE_PREDICATE_MANIFEST",
    Most[sourceManifest], sourceManifest];
  identifiers = manifestTest[[All, 1]];
  dispositions = DeleteDuplicates[manifestTest[[All, 2]]];
  partitionCounts = KeySort@Counts[manifestTest[[All, 2]]];
  manifestDigest = manifestSHA256[manifestTest];
  expectBool["PASS_SOURCE_PREDICATE_MANIFEST",
    identifiers === sourceToothIDs &&
    Length[identifiers] === Length[DeleteDuplicates[identifiers]] === 51 &&
    Sort[dispositions] === Sort[{"PRESERVED", "REPLACED_BY_STRONGER", "SCOPED_OUT"}] &&
    AllTrue[manifestTest[[All, 3]], StringStartsQ[#, "STAGE033_"] &] &&
    partitionCounts === expectedManifestCounts &&
    manifestDigest === committedManifestDigest,
    <|"Entries" -> Length[manifestTest], "Partition" -> partitionCounts,
      "Digest" -> manifestDigest|>];
  Print["      entries=", Length[manifestTest], "; partition=", partitionCounts,
    "; digest=", manifestDigest];

  Print[""];
  Print["HONEST_SCOPE: symbolic open stratum decisive; tuned locus ARGUED+FIXED-SEED-SCANNED, NOT exhaustive"];
  Print["DOWNSTREAM_ONLY: compactness, quantization, deconfinement, native +/-w current supply"];
  Print["SOURCE_FIRST_ORDERING: searched source-free first; j.A not added"];
  Print["SHEAR_DUPLICATE: NOT_APPLICABLE_NO_G"];
  Print["HARDENING-TUNED-DESCENDANT-REJECTION: PASS"];
  Print["GUARD-COUPLINGS-ENTER: PASS"];
  Print["GUARD-SEARCH-CAPABLE: PASS"];
  Print["THEORY-A VERDICT_TOKEN: ", verdictToken];
  Print["THEORY-C VERDICT_TOKEN: ", verdictToken];
  If[activeMutation =!= "",
    Print["FIRST_FAILURE=MUTATION_DID_NOT_FIRE"];
    raise["MUTATION_DID_NOT_FIRE"]
  ];
  True,
  "ledgerStage033Failure",
  Function[{message, tag}, False]
];

Print[""];
Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
If[TrueQ[ok],
  Print["OVERALL PASS: Mathematica independently reached NATIVE_P_NO_EMERGENT_GAUSS for A and C"];
  Exit[0],
  Print["OVERALL FAIL: Mathematica stage033 audit did not close"];
  Exit[1]
]
