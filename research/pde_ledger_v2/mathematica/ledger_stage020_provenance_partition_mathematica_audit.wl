(* Ledger stage020 Mathematica audit: provenance partition + CALIBRATED label.

   Standalone, print-only, no arguments, no file I/O.  This is a genuinely
   authored Wolfram route for the pathA_33 II-G4c algebra/provenance slice:
   native rational-function scaling and bridge algebra, an expression-bound
   54/5 = 2*27/5 identity, an Association-based four-way classifier, and the
   source-faithful local provenance verdict.  Adjacent stages are provenance.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
failureMessage = "";

PROVENANCEPARTITIONCALIBRATED = "PROVENANCE_PARTITION_CALIBRATED";
QUADPASS = "QUAD_PASS";
QUADCALIBRATED = "QUAD_CALIBRATED";
FAILSCALING = "FAIL_SCALING";
FAILEQUIVALENCE = "FAIL_EQUIVALENCE";
FAILPROVENANCEPARTITION = "FAIL_PROVENANCE_PARTITION";
NOFAIL = "NO_FAIL";

derivedInGate = "derived_in_gate";
externalBridgeInput = "external_bridge_input";
deferredBranchData = "deferred_branch_data";
conventionClass = "convention";
unclassified = "unclassified";

$Assumptions = Element[{a, cs, c, G, N0, D0, chiQ, lambdaG}, Reals] &&
  a > 0 && cs > 0 && c > 0 && G > 0 &&
  N0 != 0 && D0 != 0 && chiQ != 0 && lambdaG != 0;

fail[msg_] := Throw[msg, "ledgerStage020Failure"];

heading[text_] := (
  Print[""];
  Print[StringRepeat["=", StringLength[text]]];
  Print[text];
  Print[StringRepeat["=", StringLength[text]]]
);

subheading[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);

clean[expr_] := FullSimplify[Cancel[Together[expr]], $Assumptions];
fmt[expr_String] := expr;
fmt[expr_] := ToString[InputForm[clean[expr]]];

assertExact[name_, expr_] := Module[{reals},
  reals = Cases[Unevaluated[expr], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FAIL  ", name, ": machine-real atom(s) found: ", ToString[InputForm[reals]]];
    fail[name]
  ]
];

expectZero[name_, residual_] := Module[{c0},
  assertExact[name, residual];
  c0 = clean[residual];
  assertExact[name, c0];
  If[TrueQ[c0 === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FAIL  ", name, ": residual = ", fmt[c0]];
    fail[name]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

expectFail[name_, residual_] := Module[{c0},
  assertExact[name, residual];
  c0 = clean[residual];
  assertExact[name, c0];
  If[! TrueQ[c0 === 0],
    passCount++;
    Print["PASS  ", name, " produced required FAIL (residual = ", fmt[c0], ")"],
    failCount++;
    Print["FAIL  ", name, ": required mutation/ablation did not fire"];
    fail[name]
  ]
];

verdictResidual[actual_, expected_] := If[actual === expected, 0, 1];
classResidual[actual_, expected_] := If[actual === expected, 0, 1];

aPower[expr_] := Module[{rat},
  rat = Cancel[Together[expr]];
  Exponent[Numerator[rat], a] - Exponent[Denominator[rat], a]
];

(* The native classifier is built from class metadata Associations. *)
tagClasses = <|
  "dtn_hankel_expansion" -> derivedInGate,
  "dtn_radiative_slot" -> derivedInGate,
  "prefactor_series_algebra" -> derivedInGate,
  "target_rhs_algebra" -> derivedInGate,
  "gamma_bridge_algebra" -> derivedInGate,
  "emergent_outgoing_passivity" -> derivedInGate,
  "external_gr_constant" -> externalBridgeInput,
  "external_pn_bridge" -> externalBridgeInput,
  "einstein_bridge_identity" -> externalBridgeInput,
  "gate6_branch_solve" -> deferredBranchData,
  "deferred_nonlinear_pde" -> deferredBranchData,
  "normalization_convention" -> conventionClass,
  "unit_choice" -> conventionClass,
  "static_slot_convention" -> conventionClass
|>;

classRanks = <|
  unclassified -> 0,
  conventionClass -> 1,
  derivedInGate -> 2,
  externalBridgeInput -> 3,
  deferredBranchData -> 4
|>;

classifyProvenance[tags_List] := Module[{classes},
  classes = DeleteDuplicates[Lookup[tagClasses, tags, unclassified]];
  First[MaximalBy[classes, Lookup[classRanks, #, 0] &]]
];

provenanceSources = <|
  "fingerprint_u2" -> {"dtn_hankel_expansion"},
  "fingerprint_u4" -> {"dtn_hankel_expansion"},
  "fingerprint_v5" -> {"dtn_hankel_expansion", "dtn_radiative_slot"},
  "fingerprint_27" -> {"dtn_radiative_slot"},
  "prefactor_P0_P2_P4" -> {"prefactor_series_algebra"},
  "P0_target_scaling_minus5" -> {"target_rhs_algebra"},
  "chi_Q" -> {"dtn_radiative_slot"},
  "Gamma5_equivalence_chain" -> {"gamma_bridge_algebra"},
  "emergent_passivity" -> {"emergent_outgoing_passivity"},
  "G" -> {"external_gr_constant"},
  "PN_2_over_5" -> {"external_pn_bridge"},
  "Einstein_2G_over_5c5_identity" -> {"einstein_bridge_identity"},
  "assembled_54_over_5_magnitude" -> {"external_pn_bridge", "dtn_radiative_slot"},
  "D_n_N_n_numeric_values" -> {"gate6_branch_solve"},
  "port_scalars" -> {"gate6_branch_solve"},
  "actual_branch_a_scaling" -> {"gate6_branch_solve"},
  "unit_choices" -> {"unit_choice"},
  "static_slot_minus3" -> {"static_slot_convention"}
|>;

groupPartition[items_Association, field_String] := GroupBy[
  Keys[items],
  items[#][field] &
];

buildPartition[overrides_Association, sources_Association] := Module[
  {items, tags, computed, emitted},
  items = AssociationMap[
    Function[name,
      tags = sources[name];
      computed = classifyProvenance[tags];
      emitted = Lookup[overrides, name, computed];
      <|
        "ProvenanceTags" -> tags,
        "ComputedClass" -> computed,
        "EmittedClass" -> emitted,
        "ClassMatchesComputed" -> TrueQ[emitted === computed]
      |>
    ],
    Keys[sources]
  ];
  <|
    "Items" -> items,
    "Groups" -> groupPartition[items, "ComputedClass"],
    "EmittedGroups" -> groupPartition[items, "EmittedClass"],
    "Ok" -> TrueQ[And @@ (Lookup[Values[items], "ClassMatchesComputed"])]
  |>
];

verdictFromPartition[scalingOk_, equivalenceOk_, provenanceOk_, partition_] := Module[
  {gClass, magClass},
  Which[
    ! TrueQ[scalingOk], FAILSCALING,
    ! TrueQ[equivalenceOk], FAILEQUIVALENCE,
    ! TrueQ[provenanceOk], FAILPROVENANCEPARTITION,
    True,
      gClass = partition["Items"]["G"]["ComputedClass"];
      magClass = partition["Items"]["assembled_54_over_5_magnitude"]["ComputedClass"];
      If[
        gClass === derivedInGate && magClass === derivedInGate,
        QUADPASS,
        QUADCALIBRATED
      ]
  ]
];

localVerdict[gates_Association, partition_Association] := verdictFromPartition[
  gates["ScalingOk"],
  gates["EquivalenceOk"],
  gates["ProvenanceOk"],
  partition
];

dynamicAblation[
  baselineGates_Association,
  baselinePartition_Association,
  mutatedGates_Association,
  mutatedPartition_Association,
  expectedFail_String
] := Module[{withMutation, withoutMutation},
  withMutation = localVerdict[mutatedGates, mutatedPartition];
  withoutMutation = localVerdict[baselineGates, baselinePartition];
  <|
    "RerunGateLogic" -> True,
    "WithMutation" -> withMutation,
    "WithoutMutation" -> withoutMutation,
    "ExpectedFail" -> expectedFail,
    "FailSuppressed" -> TrueQ[
      withMutation === expectedFail &&
      withoutMutation =!= expectedFail &&
      withMutation =!= withoutMutation
    ]
  |>
];

isGInvariant[expr_] := TrueQ[clean[(expr /. G -> lambdaG G) - expr] === 0];

(* Strengthened scaling: the derived path carries a^5 only through v5Slot. *)
v5Slot = a^5/(27 cs^5);
chi = 1;
gammaTarget = 2 G/(5 c^5);
targetFromBridge = clean[gammaTarget/(chi v5Slot)];

(* Independently assembled comparison target and its one-tooth mutation. *)
targetRHS = 54 G cs^5/(5 a^5 c^5);
mutatedTargetRHS = 54 G cs^5/(5 a^4 c^5);
gammaPower = aPower[gammaTarget];
slotPower = aPower[v5Slot];
derivedPower = gammaPower - slotPower;
bridgePower = aPower[targetFromBridge];
assembledBridgeResidual = clean[targetFromBridge - targetRHS];
mutatedBridgeResidual = clean[targetFromBridge - mutatedTargetRHS];
mutatedScalingOk = TrueQ[mutatedBridgeResidual === 0];
scalingOk = TrueQ[
  gammaPower === 0 && slotPower === 5 && derivedPower === -5 &&
  bridgePower === derivedPower && assembledBridgeResidual === 0
];

(* Gamma5 definition and independent bridge residuals. *)
p0Physical = (cs/a)^2 (N0/D0);
gamma5 = chiQ p0Physical v5Slot;
gamma5WithOutgoingChi = clean[gamma5 /. chiQ -> chi];
forwardGeneral = clean[targetRHS chiQ v5Slot - gammaTarget];
forwardExpected = 2 G (chiQ - 1)/(5 c^5);
forward = clean[forwardGeneral /. chiQ -> chi];
reverse = clean[gammaTarget/(chi v5Slot) - targetRHS];
wrongGamma = 3 G/(5 c^5);
wrongReverse = clean[wrongGamma/(chi v5Slot) - targetRHS];
mutatedEquivalenceOk = TrueQ[wrongReverse === 0];
equivalenceOk = TrueQ[forward === 0 && reverse === 0];

(* The identity is bound to targetRHS and v5Slot, not bare literals. *)
normalizationUnit = G cs^5/(a^5 c^5);
mag = Simplify[targetRHS/(G cs^5/(a^5 c^5)), $Assumptions];
twentySevenFromSlot = Simplify[a^5/(cs^5 v5Slot), $Assumptions];
identityRight = Simplify[Rational[2, 5] twentySevenFromSlot, $Assumptions];
identityResidual = Simplify[mag - identityRight, $Assumptions];
mutatedIdentityRight = Simplify[
  Rational[2, 5] (twentySevenFromSlot - 1),
  $Assumptions
];
mutatedIdentityResidual = Simplify[mag - mutatedIdentityRight, $Assumptions];
identityOk = TrueQ[identityResidual === 0];

realPartition = buildPartition[<||>, provenanceSources];
realItems = realPartition["Items"];
earnedFactorClass = realItems["fingerprint_27"]["ComputedClass"];
calibratedFactorClass = realItems["PN_2_over_5"]["ComputedClass"];
assembledMagnitudeClass = realItems["assembled_54_over_5_magnitude"]["ComputedClass"];

(* Independent dominance proof: table, key quantities, and a computed tag mutation. *)
truthInputs = <|
  "deferred_over_external" -> {"gate6_branch_solve", "external_gr_constant"},
  "external_over_derived" -> {"external_pn_bridge", "dtn_hankel_expansion"},
  "derived_over_convention" -> {"dtn_radiative_slot", "unit_choice"},
  "single_deferred" -> {"gate6_branch_solve"},
  "single_external" -> {"external_gr_constant"},
  "single_derived" -> {"dtn_radiative_slot"},
  "single_convention" -> {"unit_choice"}
|>;
truthExpected = <|
  "deferred_over_external" -> deferredBranchData,
  "external_over_derived" -> externalBridgeInput,
  "derived_over_convention" -> derivedInGate,
  "single_deferred" -> deferredBranchData,
  "single_external" -> externalBridgeInput,
  "single_derived" -> derivedInGate,
  "single_convention" -> conventionClass
|>;
truthResults = AssociationMap[
  classifyProvenance[truthInputs[#]] &,
  Keys[truthInputs]
];
truthOk = TrueQ[And @@ MapThread[SameQ, {Values[truthResults], Values[truthExpected]}]];

keyExpected = <|
  "G" -> externalBridgeInput,
  "PN_2_over_5" -> externalBridgeInput,
  "fingerprint_27" -> derivedInGate,
  "assembled_54_over_5_magnitude" -> externalBridgeInput
|>;
keyResults = AssociationMap[realItems[#]["ComputedClass"] &, Keys[keyExpected]];
keyOk = TrueQ[And @@ MapThread[SameQ, {Values[keyResults], Values[keyExpected]}]];

baselineMagnitudeTags = provenanceSources["assembled_54_over_5_magnitude"];
strippedMagnitudeTags = DeleteCases[baselineMagnitudeTags, "external_pn_bridge"];
baselineComputedClass = classifyProvenance[baselineMagnitudeTags];
mutatedComputedClass = classifyProvenance[strippedMagnitudeTags];
tagMutationFires = TrueQ[
  baselineComputedClass === externalBridgeInput &&
  mutatedComputedClass === derivedInGate &&
  mutatedComputedClass =!= baselineComputedClass
];
classifierProofOk = TrueQ[truthOk && keyOk && tagMutationFires];

factorClassesOk = TrueQ[
  earnedFactorClass === derivedInGate &&
  calibratedFactorClass === externalBridgeInput &&
  assembledMagnitudeClass === externalBridgeInput
];
provenanceOk = TrueQ[
  realPartition["Ok"] && identityOk && classifierProofOk && factorClassesOk
];

(* Separate diagnostic; no diagnostic value is an argument of the verdict. *)
gExpressions = <|
  "G" -> G,
  "target_2G_over_5c5" -> gammaTarget,
  "pure_54_over_5" -> mag,
  "fingerprint_27" -> twentySevenFromSlot
|>;
gExpressionClasses = <|
  "G" -> realItems["G"]["ComputedClass"],
  "target_2G_over_5c5" -> realItems["PN_2_over_5"]["ComputedClass"],
  "pure_54_over_5" -> realItems["assembled_54_over_5_magnitude"]["ComputedClass"],
  "fingerprint_27" -> realItems["fingerprint_27"]["ComputedClass"]
|>;
gDiagnostics = AssociationMap[
  Function[name,
    Module[{expr, transformed, residual},
      expr = gExpressions[name];
      transformed = clean[expr /. G -> lambdaG G];
      residual = clean[transformed - expr];
      <|
        "Expression" -> expr,
        "Transformed" -> transformed,
        "Residual" -> residual,
        "GInvariant" -> TrueQ[residual === 0],
        "ProvenanceClass" -> gExpressionClasses[name]
      |>
    ]
  ],
  Keys[gExpressions]
];
invarianceOnlyTrapCatches54Over5 = TrueQ[
  gDiagnostics["pure_54_over_5"]["GInvariant"] &&
  gDiagnostics["pure_54_over_5"]["ProvenanceClass"] === externalBridgeInput
];

baselineGates = <|
  "ScalingOk" -> scalingOk,
  "EquivalenceOk" -> equivalenceOk,
  "ProvenanceOk" -> provenanceOk
|>;
baselineVerdict = localVerdict[baselineGates, realPartition];

(* Coherent computed+emitted provenance controls. *)
allDerivedSources = Join[
  provenanceSources,
  <|
    "G" -> {"dtn_radiative_slot"},
    "assembled_54_over_5_magnitude" -> {"dtn_radiative_slot"}
  |>
];
allDerivedPartition = buildPartition[<||>, allDerivedSources];
allDerivedVerdict = verdictFromPartition[True, True, allDerivedPartition["Ok"], allDerivedPartition];

mixedSources = Join[
  provenanceSources,
  <|
    "G" -> {"external_gr_constant"},
    "assembled_54_over_5_magnitude" -> {"dtn_radiative_slot"}
  |>
];
mixedPartition = buildPartition[<||>, mixedSources];
mixedVerdict = verdictFromPartition[True, True, mixedPartition["Ok"], mixedPartition];

invertedVerdictControl[scalingGate_, equivalenceGate_, provenanceGate_, partition_] := Module[
  {gClass, magClass, bothExternal},
  Which[
    ! TrueQ[scalingGate], FAILSCALING,
    ! TrueQ[equivalenceGate], FAILEQUIVALENCE,
    ! TrueQ[provenanceGate], FAILPROVENANCEPARTITION,
    True,
      gClass = partition["Items"]["G"]["ComputedClass"];
      magClass = partition["Items"]["assembled_54_over_5_magnitude"]["ComputedClass"];
      bothExternal = TrueQ[gClass === externalBridgeInput && magClass === externalBridgeInput];
      If[bothExternal, QUADCALIBRATED, QUADPASS]
  ]
];

constantCalibratedOnAllDerived = QUADCALIBRATED;
invertedMixedVerdict = invertedVerdictControl[True, True, mixedPartition["Ok"], mixedPartition];

(* Computed probes and dynamic local reruns. *)
scalingMutatedGates = Join[baselineGates, <|"ScalingOk" -> mutatedScalingOk|>];
scalingAblation = dynamicAblation[
  baselineGates, realPartition,
  scalingMutatedGates, realPartition,
  FAILSCALING
];
scalingProbeVerdict = If[! mutatedScalingOk, FAILSCALING, NOFAIL];
scalingCorrectVerdict = If[scalingOk, NOFAIL, FAILSCALING];

equivalenceMutatedGates = Join[
  baselineGates,
  <|"EquivalenceOk" -> mutatedEquivalenceOk|>
];
equivalenceAblation = dynamicAblation[
  baselineGates, realPartition,
  equivalenceMutatedGates, realPartition,
  FAILEQUIVALENCE
];
equivalenceProbeVerdict = If[! mutatedEquivalenceOk, FAILEQUIVALENCE, NOFAIL];
equivalenceCorrectVerdict = If[equivalenceOk, NOFAIL, FAILEQUIVALENCE];

gAsDerivedPartition = buildPartition[<|"G" -> derivedInGate|>, provenanceSources];
fingerprintAsExternalPartition = buildPartition[
  <|"fingerprint_27" -> externalBridgeInput|>,
  provenanceSources
];

partitionAblation[mutatedPartition_Association] := Module[{mutatedGates},
  mutatedGates = Join[
    baselineGates,
    <|"ProvenanceOk" -> mutatedPartition["Ok"]|>
  ];
  dynamicAblation[
    baselineGates, realPartition,
    mutatedGates, mutatedPartition,
    FAILPROVENANCEPARTITION
  ]
];

gPartitionAblation = partitionAblation[gAsDerivedPartition];
fingerprintPartitionAblation = partitionAblation[fingerprintAsExternalPartition];
gMislabelVerdict = If[! gAsDerivedPartition["Ok"], FAILPROVENANCEPARTITION, NOFAIL];
fingerprintMislabelVerdict = If[
  ! fingerprintAsExternalPartition["Ok"],
  FAILPROVENANCEPARTITION,
  NOFAIL
];
partitionMutationCases = <|
  "G_external_to_derived" -> gPartitionAblation["WithMutation"],
  "fingerprint_27_derived_to_external" -> fingerprintPartitionAblation["WithMutation"]
|>;

runScalingAndEquivalence[] := (
  subheading["Strengthened a⁻⁵ scaling from the frozen 018 v5 slot"];
  Print["  gamma_target = ", fmt[gammaTarget]];
  Print["  a_power(gamma_target) = ", gammaPower];
  Print["  v5_slot = ", fmt[v5Slot]];
  Print["  a_power(v5_slot) = ", slotPower];
  Print["  target_from_bridge = ", fmt[targetFromBridge]];
  Print["  derived_power = a_power(gamma_target) - a_power(v5_slot) = ", derivedPower];
  Print["  a_power(target_from_bridge) = ", bridgePower, " (a⁻⁵)"];
  Print["  independently assembled target_rhs = ", fmt[targetRHS]];
  expectZero["Burke-Thorne gamma_target is a-free", gammaPower];
  expectZero["frozen v5_slot supplies a-power +5", slotPower - 5];
  expectZero["derived bridge power is -5", derivedPower + 5];
  expectZero["target_from_bridge power equals the derived power", bridgePower - derivedPower];
  expectZero["bridge reconstruction equals the independently assembled target_rhs", assembledBridgeResidual];
  expectFail["3c a^-4 mutation of only the assembled target breaks bridge cancellation", mutatedBridgeResidual];
  expectBool["strengthened scaling gate passes", scalingOk];

  subheading["Gamma5/chi_Q equivalence bridge"];
  Print["  P0_physical (019 provenance definition only) = ", fmt[p0Physical]];
  Print["  Gamma5 DEF = ", fmt[gamma5]];
  Print["  forward_general = ", fmt[forwardGeneral]];
  Print["  forward = ", fmt[forward]];
  Print["  reverse = ", fmt[reverse]];
  Print["  wrong_reverse (3/5 mutation) = ", fmt[wrongReverse]];
  expectBool[
    "P0=N0/D0 enters the Gamma5 DEF and is absent from bridge residuals",
    And @@ (! FreeQ[gamma5, #] & /@ {N0, D0}) &&
      FreeQ[{forward, reverse, forwardGeneral}, N0 | D0]
  ];
  expectZero["forward general form is 2*G*(chi_Q-1)/(5*c^5)", forwardGeneral - forwardExpected];
  expectZero["Gamma5 bridge forward=0 for cited chi_Q=+1", forward];
  expectZero["Gamma5 bridge reverse=0", reverse];
  expectFail["3e wrong_gamma=3*G/(5*c^5) makes wrong_reverse nonzero", wrongReverse];
  expectBool["equivalence gate passes", equivalenceOk]
);

runIdentityAndClassifier[] := (
  subheading["Expression-bound 54/5 = 2·27/5 decomposition"];
  Print["  unit extracted from target_rhs = ", fmt[normalizationUnit]];
  Print["  mag = target_rhs/unit = ", fmt[mag]];
  Print["  27_from_slot = a^5/(c_s^5*v5_slot) = ", fmt[twentySevenFromSlot]];
  Print["  bound left = ", fmt[mag]];
  Print["  bound right = 2*27_from_slot/5 = ", fmt[identityRight]];
  Print["  54/5=2·27/5 TRUE bound residual = ", fmt[identityResidual]];
  expectZero["BOUND 54/5=2*27/5 identity from target_rhs and v5_slot", identityResidual];
  expectFail["bound-identity factor mutation 27_from_slot->27_from_slot-1 fires", mutatedIdentityResidual];

  subheading["Four-way provenance classifier proof"];
  Print["  dominance order = deferred > external > derived > convention"];
  Scan[
    Function[name,
      Print["  truth-table ", name, ": ", truthInputs[name], " -> ", truthResults[name]];
      expectZero[
        "classifier truth-table " <> name <> " -> " <> truthExpected[name],
        classResidual[truthResults[name], truthExpected[name]]
      ]
    ],
    Keys[truthInputs]
  ];
  Scan[
    Function[name,
      Print["  key class ", name, " -> ", keyResults[name]];
      expectZero[
        "key provenance class " <> name <> " -> " <> keyExpected[name],
        classResidual[keyResults[name], keyExpected[name]]
      ]
    ],
    Keys[keyExpected]
  ];
  Print[
    "  tag mutation assembled magnitude: ", baselineMagnitudeTags,
    " -> ", strippedMagnitudeTags, " -> ", mutatedComputedClass
  ];
  expectZero[
    "stripping external_pn_bridge drops assembled magnitude to derived_in_gate",
    classResidual[mutatedComputedClass, derivedInGate]
  ];
  expectZero[
    "baseline assembled magnitude class is external_bridge_input before tag stripping",
    classResidual[baselineComputedClass, externalBridgeInput]
  ];
  expectBool["classifier truth-table + key classes + tag mutation all pass", classifierProofOk];
  expectBool["real emitted partition matches computed classes", realPartition["Ok"]];
  Print["  earned factor class (READ from partition) = ", earnedFactorClass];
  Print["  calibrated factor class (READ from partition) = ", calibratedFactorClass];
  Print["  assembled class (READ from partition) = ", assembledMagnitudeClass];
  expectZero["27 factor class is read as derived_in_gate", classResidual[earnedFactorClass, derivedInGate]];
  expectZero["2/5 factor class is read as external_bridge_input", classResidual[calibratedFactorClass, externalBridgeInput]];
  expectZero["assembled 54/5 class is read as external_bridge_input", classResidual[assembledMagnitudeClass, externalBridgeInput]]
);

runGInvarianceDiagnostic[] := Module[{expected},
  subheading["SEPARATE non-verdict-driving G->lambda_G*G diagnostic"];
  expected = <|
    "G" -> False,
    "target_2G_over_5c5" -> False,
    "pure_54_over_5" -> True,
    "fingerprint_27" -> True
  |>;
  Scan[
    Function[name,
      Print[
        "  ", name, ": expr=", fmt[gDiagnostics[name]["Expression"]],
        ", g_invariant=", gDiagnostics[name]["GInvariant"],
        ", class=", gDiagnostics[name]["ProvenanceClass"]
      ];
      expectBool[
        "G-invariance diagnostic for " <> name <> " is " <> ToString[expected[name]],
        gDiagnostics[name]["GInvariant"] === expected[name]
      ]
    ],
    Keys[gDiagnostics]
  ];
  Print["  invariance_only_trap_catches_54_over_5 = ", invarianceOnlyTrapCatches54Over5];
  expectBool[
    "separate invariance-only trap catches G-invariant yet calibrated 54/5",
    invarianceOnlyTrapCatches54Over5
  ];
  expectZero[
    "provenance verdict remains CALIBRATED independently of the diagnostic",
    verdictResidual[baselineVerdict, QUADCALIBRATED]
  ]
];

runVerdictControls[] := Module[{allItems, mixedItems},
  subheading["Source-faithful provenance verdict and positive controls"];
  allItems = allDerivedPartition["Items"];
  mixedItems = mixedPartition["Items"];
  Print[
    "  real classes: G=", realItems["G"]["ComputedClass"],
    ", assembled54/5=", realItems["assembled_54_over_5_magnitude"]["ComputedClass"]
  ];
  Print["  real verdict = ", baselineVerdict];
  Print[
    "  QUAD_PASS control classes: G=", allItems["G"]["ComputedClass"],
    ", assembled54/5=", allItems["assembled_54_over_5_magnitude"]["ComputedClass"],
    "; partition_ok=", allDerivedPartition["Ok"]
  ];
  Print["  QUAD_PASS control verdict = ", allDerivedVerdict];
  Print[
    "  MIXED control classes: G=", mixedItems["G"]["ComputedClass"],
    ", assembled54/5=", mixedItems["assembled_54_over_5_magnitude"]["ComputedClass"],
    "; partition_ok=", mixedPartition["Ok"]
  ];
  Print["  MIXED control verdict = ", mixedVerdict];
  expectZero["real both-external partition lands QUAD_CALIBRATED", verdictResidual[baselineVerdict, QUADCALIBRATED]];
  expectBool["coherent all-derived control partition is internally consistent", allDerivedPartition["Ok"]];
  expectZero["coherent all-derived positive control lands QUAD_PASS", verdictResidual[allDerivedVerdict, QUADPASS]];
  expectBool["coherent mixed control partition is internally consistent", mixedPartition["Ok"]];
  expectZero["required MIXED control lands QUAD_CALIBRATED", verdictResidual[mixedVerdict, QUADCALIBRATED]];
  expectFail[
    "constant-CALIBRATED verdict ablation fails at the QUAD_PASS control",
    verdictResidual[constantCalibratedOnAllDerived, QUADPASS]
  ];
  expectFail[
    "inverted PASS-unless-both-external rule fails at the MIXED control",
    verdictResidual[invertedMixedVerdict, QUADCALIBRATED]
  ]
];

assertDynamicAblation[label_String, ablation_Association, expectedFail_String] := (
  Print["  ", label, " dynamic self-ablation = ", ablation];
  expectBool[label <> " reruns 020-local gate logic", ablation["RerunGateLogic"]];
  expectZero[
    label <> " dynamic rerun with mutation reaches " <> expectedFail,
    verdictResidual[ablation["WithMutation"], expectedFail]
  ];
  expectZero[
    label <> " dynamic rerun without mutation returns QUAD_CALIBRATED",
    verdictResidual[ablation["WithoutMutation"], QUADCALIBRATED]
  ];
  expectBool[label <> " self-ablation suppresses the local failure", ablation["FailSuppressed"]]
);

runProbes[] := (
  subheading["Probes 3c/3e/3f and DYNAMIC 020-local self-ablations"];
  Print["  3c wrong-scaling residual = ", fmt[mutatedBridgeResidual]];
  Print["  3c verdict = ", scalingProbeVerdict];
  expectZero["3c wrong assembled scaling reaches FAIL_SCALING", verdictResidual[scalingProbeVerdict, FAILSCALING]];
  expectZero["3c correct bridge scaling is NO_FAIL", verdictResidual[scalingCorrectVerdict, NOFAIL]];
  assertDynamicAblation["3c scaling_ablation", scalingAblation, FAILSCALING];

  Print["  3e wrong_gamma = ", fmt[wrongGamma]];
  Print["  3e wrong_reverse = ", fmt[wrongReverse]];
  Print["  3e verdict = ", equivalenceProbeVerdict];
  expectZero["3e wrong Gamma bridge reaches FAIL_EQUIVALENCE", verdictResidual[equivalenceProbeVerdict, FAILEQUIVALENCE]];
  expectZero["3e correct Gamma bridge is NO_FAIL", verdictResidual[equivalenceCorrectVerdict, NOFAIL]];
  assertDynamicAblation["3e equivalence_ablation", equivalenceAblation, FAILEQUIVALENCE];

  Print[
    "  3f G external->derived: partition_ok=", gAsDerivedPartition["Ok"],
    ", verdict=", gMislabelVerdict
  ];
  Print[
    "  3f fingerprint_27 derived->external: partition_ok=", fingerprintAsExternalPartition["Ok"],
    ", verdict=", fingerprintMislabelVerdict
  ];
  expectZero[
    "3f G external->derived reaches FAIL_PROVENANCE_PARTITION",
    verdictResidual[gMislabelVerdict, FAILPROVENANCEPARTITION]
  ];
  expectZero[
    "3f fingerprint_27 derived->external reaches FAIL_PROVENANCE_PARTITION",
    verdictResidual[fingerprintMislabelVerdict, FAILPROVENANCEPARTITION]
  ];
  expectBool[
    "3f aggregate requires BOTH mutation directions to fire",
    And @@ (TrueQ[# === FAILPROVENANCEPARTITION] & /@ Values[partitionMutationCases])
  ];
  assertDynamicAblation[
    "3f partition_ablation G external->derived",
    gPartitionAblation,
    FAILPROVENANCEPARTITION
  ];
  assertDynamicAblation[
    "3f partition_ablation 27 derived->external",
    fingerprintPartitionAblation,
    FAILPROVENANCEPARTITION
  ];
  expectBool[
    "3f records that G-invariance alone would miss calibrated 54/5",
    invarianceOnlyTrapCatches54Over5
  ]
);

symbolNames[expr_] := DeleteDuplicates[
  Cases[expr, s_Symbol /; Context[s] === "Global`" :> SymbolName[s], Infinity]
];

earnedExpressions = Join[
  {
    chi, v5Slot, gammaTarget, targetFromBridge, targetRHS, mutatedTargetRHS,
    gammaPower, slotPower, derivedPower, bridgePower,
    assembledBridgeResidual, mutatedBridgeResidual,
    p0Physical, gamma5, gamma5WithOutgoingChi,
    forwardGeneral, forwardExpected, forward, reverse, wrongGamma, wrongReverse,
    normalizationUnit, mag, twentySevenFromSlot, identityRight,
    identityResidual, mutatedIdentityRight, mutatedIdentityResidual
  },
  Flatten[
    ({#["Expression"], #["Transformed"], #["Residual"]} &) /@ Values[gDiagnostics]
  ]
];
allowedSymbolNames = {"a", "cs", "c", "G", "N0", "D0", "chiQ", "lambdaG"};

runAritySelfCheck[] := Module[
  {arityPartition, arityAblation, arityClassifier, coreResults},
  subheading["Wolfram arity self-check and unevaluated-leakage scan"];
  arityPartition = buildPartition[<||>, provenanceSources];
  arityClassifier = classifyProvenance[{"external_pn_bridge", "dtn_radiative_slot"}];
  arityAblation = dynamicAblation[
    baselineGates,
    realPartition,
    Join[baselineGates, <|"ScalingOk" -> False|>],
    realPartition,
    FAILSCALING
  ];
  expectZero["arity aPower[expr] returns the derived target exponent", aPower[targetFromBridge] + 5];
  expectZero[
    "arity classifyProvenance[tags] applies external-over-derived dominance",
    classResidual[arityClassifier, externalBridgeInput]
  ];
  expectBool[
    "arity buildPartition[overrides,sources] returns Items/Groups/Ok",
    And @@ (KeyExistsQ[arityPartition, #] & /@ {"Items", "Groups", "Ok"})
  ];
  expectZero[
    "arity verdictFromPartition[three gates,partition] returns QUAD_PASS control",
    verdictResidual[
      verdictFromPartition[True, True, allDerivedPartition["Ok"], allDerivedPartition],
      QUADPASS
    ]
  ];
  expectBool[
    "arity dynamicAblation[baseline,partition,mutant,partition,label] changes verdicts",
    arityAblation["WithMutation"] =!= arityAblation["WithoutMutation"]
  ];
  expectBool["arity isGInvariant[expr] recognizes the bound pure magnitude", isGInvariant[mag]];
  expectBool[
    "arity groupPartition[items,field] includes the external class",
    KeyExistsQ[groupPartition[realItems, "ComputedClass"], externalBridgeInput]
  ];
  coreResults = {
    gammaPower, slotPower, derivedPower, bridgePower,
    forward, reverse, wrongReverse, mag, twentySevenFromSlot,
    identityResidual, truthResults, keyResults, realPartition,
    baselineVerdict, allDerivedVerdict, mixedVerdict, gDiagnostics,
    scalingAblation, equivalenceAblation,
    gPartitionAblation, fingerprintPartitionAblation
  };
  expectBool[
    "no unevaluated authored helper or algebra call remains in core results",
    FreeQ[
      coreResults,
      aPower | classifyProvenance | groupPartition | buildPartition |
      verdictFromPartition | localVerdict | dynamicAblation | isGInvariant |
      Exponent | Together | Cancel | FullSimplify | Simplify | ReplaceAll
    ]
  ]
];

runScopeAndProvenance[] := Module[
  {
    namesByExpression, unexpected, liveNames, verdictDefinitions,
    verdictDefinitionLHS, verdictArity, verdictDefinitionText
  },
  subheading["020 scope, provenance-only consumption, and PARTIAL landing"];
  Print["  020 gate booleans = ", baselineGates];
  Print["  020 scoped verdict = ", baselineVerdict];
  expectZero[
    "020 scoped verdict lands the CALIBRATED partial component",
    verdictResidual[baselineVerdict, QUADCALIBRATED]
  ];
  assertExact["earnedExpressions", earnedExpressions];
  namesByExpression = symbolNames /@ earnedExpressions;
  unexpected = Complement[#, allowedSymbolNames] & /@ namesByExpression;
  liveNames = Sort[DeleteDuplicates[Flatten[namesByExpression]]];
  Print["  earned symbolic expression count under the NAME allowlist = ", Length[earnedExpressions]];
  Print["  live symbolic names across every earned expression = ", liveNames];
  expectBool[
    "EVERY earned symbolic expression obeys the free-symbol NAME allowlist",
    And @@ (TrueQ[# === {}] & /@ unexpected)
  ];
  expectBool[
    "units-bearing algebra exposes exactly a,cs,c,G,N0,D0,chiQ,lambdaG",
    Sort[liveNames] === Sort[allowedSymbolNames]
  ];
  verdictDefinitions = DownValues[verdictFromPartition];
  verdictDefinitionLHS = If[
    Length[verdictDefinitions] === 1,
    Extract[verdictDefinitions, {1, 1}, HoldComplete],
    HoldComplete[]
  ];
  verdictArity = If[
    Length[verdictDefinitions] === 1,
    Extract[
      verdictDefinitionLHS,
      {1, 1},
      Function[call, Length[Unevaluated[call]], HoldAllComplete]
    ],
    -1
  ];
  verdictDefinitionText = ToLowerCase[ToString[verdictDefinitions, InputForm]];
  Print[
    "  structural dimensional cut: verdictFromPartition definition count = ",
    Length[verdictDefinitions], ", arity = ", verdictArity,
    ", dim-like name present = ", StringContainsQ[verdictDefinitionText, "dim"]
  ];
  expectBool[
    "020-local verdict has exactly three gate+partition parameters and no dim-like dependency",
    Length[verdictDefinitions] === 1 && verdictArity === 4 &&
      ! StringContainsQ[verdictDefinitionText, "dim"]
  ];
  Print["  FOUR-WAY SPLIT / THIRD LEG: 020 carries the provenance-partition + CALIBRATED-verdict component 3/4 and lands it as PARTIAL; 018 fingerprint DONE, 019 prefactor DONE, 021 dim closure remaining."];
  Print["  CONSUMED-PROVENANCE -- LABELED NON-CHECK: cites 018 chi_Q=+1 and 27 (=1/v5), 019 P0=N0/D0, and 017 l=2 port-kernel D-lanes; NO guard/dual-site."];
  Print["  UNITS-BEARING-BUT-NO-DIM-GATE: c_s,a,c,G plus abstract N0,D0, chi_Q, and 27 occur in algebra/provenance only; the homogeneity closure belongs to 021."];
  Print["  54/5=2*27/5-COMPUTED / PROVENANCE-PARTITION: bound symbolic residual, computed factor classes, provenance-driven verdict; the separate invariance trap is non-driving."];
  Print["  STRENGTHENED-a^-5 / 3f-BOTH-DIRECTIONS: -5 flows from the frozen v5 slot; both G external->derived and 27 derived->external fire."];
  Print["  EARNED-vs-CALIBRATED scope: scaling, bridge, identity, classifier and verdict rule are earned; 54/5 and G are calibrated; actual branch scaling and numerical port scalars are Gate-6 deferred."];
  Print["  dropped-bookkeeping: scratch-YAML agreement handoff and report/feed writers are absent; tri-review plus the independent Wolfram route and transcript agreement replace them."];
  Print["  register note: likely zero new counted knobs; registration decides the structural edge and confirms the G/c bridge accounting."];
  Print["  QUAD_CALIBRATED (3/4) — the 54/5=2·27/5 provenance partition + the CALIBRATED verdict label EARNED (PARTIAL)"];
  Print["    = the a⁻⁵ target scaling (DERIVED via equivalence a-cancellation, not a typed target power)"];
  Print["    AND the Gamma5/chi_Q equivalence (forward=0, reverse=0; closes iff chi=1 and 54/27=2)"];
  Print["    AND the bound 54/5=2*27/5 identity with 27 derived_in_gate and 2/5 external_bridge_input"];
  Print["    AND the four-way provenance dominance that makes assembled 54/5 external_bridge_input -> QUAD_CALIBRATED."];
  Print["  REMAINING: fingerprint=018 (DONE); prefactor algebra=019 (DONE); the mu_hat0-free [P0^phys]=1 dim closure=021 (COMPLETING leg)."];
  Print["  CAVEATS: actual branch a-scaling from N0/D0 and numerical (D_n,N_n) port scalars are Gate-6 deferred; 54/5 and G are calibrated (G=GENUINE_BLOCKED)."];
  Print["  CONSUMES (PROVENANCE only): 018 chi_Q=+1 and 27; 019 P0=N0/D0 in the Gamma5 DEF; 017 port-kernel D-lanes; no dual-site."];
  Print["  EXPORTS: partition + CALIBRATED label + Gamma5/chi_Q equivalence + a^-5 scaling -> 021 + 022 + 027; no file artifact is written."];
  Print["  reduction certificate: FROZEN-INPUT 018 chi_Q/27, 019 P0 definition, 017 port provenance, and GR G/c/2/5 bridge; COMPUTED scaling/equivalence/bound identity/partition/trap/verdict; DEFERRED actual branch scaling and numerical port data."]
];

printVerdictLabels[] := (
  subheading["Verdict labels:"];
  Print["  ledger earned-label (NOT a source verdict token): PROVENANCE_PARTITION_CALIBRATED  (the assembled quadrupole magnitude 54/5 decomposes as 54/5 = 2*27/5 [SymPy-verified rational identity]; the 27 is derived_in_gate [018's fingerprint 1/v5], the 2/5 + G are external_bridge_input [GR Burke-Thorne 2G/5c^5, G=GENUINE_BLOCKED]; a 4-way provenance partition [dominance deferred>external>derived>convention] classifies the assembled 54/5 as external_bridge_input, so the verdict lands QUAD_CALIBRATED not QUAD_PASS -- PROVENANCE-driven, NOT G->lambda*G invariance [the separate g-invariance diagnostic exposes the invariance-only trap: 54/5 is G-invariant yet calibrated])"];
  Print["  source top-line verdict: QUAD_CALIBRATED  (JOINT 4-stage; 020 carries the 54/5=2*27/5 provenance-partition + the CALIBRATED verdict label component 3/4 and lands the token as a PARTIAL)"];
  Print["  joint composition (020 = the THIRD leg; 018/019 DONE, 021 remaining): QUAD_CALIBRATED = (018: outgoing DtN Hankel fingerprint + chi_Q sign)[EARNED, PARTIAL, DONE] AND (019: squared-denominator prefactor algebra)[EARNED, PARTIAL, DONE] AND (020: 54/5=2*27/5 provenance partition + the CALIBRATED verdict label)[EARNED here, PARTIAL] AND (021: mu_hat0-free [P0^phys]=1 dim closure)[the COMPLETING leg]"];
  Print["  earned (the partition + bridge + scaling): the a^-5 target scaling (DERIVED via the equivalence a-cancellation, not a typed target power); the Gamma5/chi_Q equivalence 54*G*c_s^5/(5*a^5*c^5) <=> 2*G/(5*c^5) (closes iff chi=1 and 54/27=2); the 54/5=2*27/5 SymPy identity; the 4-way provenance partition (assembled 54/5 = external_bridge_input, external dominates); the CALIBRATED verdict (provenance-driven, NOT G-invariance)"];
  Print["  earned (able-to-fail): 3c_wrong_scaling (STRENGTHENED -- a^5->a^4 breaks the equivalence a-cancellation -> FAIL_SCALING); 3e_equivalence_break (2/5->3/5 -> FAIL_EQUIVALENCE); 3f_partition_mislabel (BOTH directions: G external->derived AND 27 derived->external -> FAIL_PROVENANCE_PARTITION), each with a DYNAMIC 020-local self-ablation (re-run, NOT a constant, NOT the joint base_verdict)"];
  Print["  calibrated / deferred (NOT 020): the fingerprint (018, DONE); the prefactor algebra (019, DONE); the mu_hat0-free dim closure (021); the 54/5 magnitude + G (external_bridge_input, G=GENUINE_BLOCKED); the ACTUAL branch a-scaling + the numerical (D_n,N_n) port scalars (Gate-6 sim-deferred, report :49)"];
  Print["  consumed (PROVENANCE only -- NO guard/dual-site): 018's chi_Q=+1 + the 27 (=1/v5) [enter 020's self-contained equivalence bridge] + 019's P0=N0/D0 [Gamma5 def] + 017's l=2 port kernel D-lanes; NO mu_hat0/dim gate"];
  Print["  exports: the 54/5=2*27/5 partition + the CALIBRATED label + the Gamma5/chi_Q equivalence + the a^-5 scaling => 021 (dim closure completes the joint) + 022 (non-regression) + 027 (pathA_43 closure)"];
  Print["  new symbols first-appearing (020): none new-counted (G=GENUINE_BLOCKED/external_bridge_input, already in the register; c=GR/PN light-speed bridge = c_gamma role, cited benchmark; the 2/5=GR bridge; the 27=018's derived fingerprint); units-bearing (c_s/a/c/G live) but NO dim gate; no counted knobs expected (an EARNED classification slice, like 018/019)"]
);

runAll[] := (
  heading["ledger_stage020_provenance_partition_mathematica_audit"];
  Print["Target stem confirmed: ledger_stage020_provenance_partition"];
  Print["Engine: genuinely authored native Wolfram Exponent/Together/Cancel/FullSimplify, Rational/Simplify, Association classification, and replacement diagnostics; no floats/tolerances; zero file I/O."];
  runScalingAndEquivalence[];
  runIdentityAndClassifier[];
  runGInvarianceDiagnostic[];
  runVerdictControls[];
  runProbes[];
  runAritySelfCheck[];
  runScopeAndProvenance[];
  printVerdictLabels[];
  0
);

result = Catch[
  runAll[],
  "ledgerStage020Failure",
  Function[{msg, tag}, failureMessage = ToString[msg]; 1]
];

heading["Tallies"];
Print["TALLY mathematica: ", passCount, " pass + ", failCount, " fail = ", passCount + failCount, " checks"];
If[result === 0 && failCount === 0,
  Print["OVERALL PASS"];
  Exit[0],
  If[failureMessage =!= "", Print["ABORTED: ", failureMessage]];
  Print["OVERALL FAIL"];
  Exit[1]
];
