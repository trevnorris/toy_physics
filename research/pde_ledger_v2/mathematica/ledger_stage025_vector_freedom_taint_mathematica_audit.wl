(* Ledger stage025 Wolfram audit: density-port vector-freedom taint proof.

   Standalone, print-only, no arguments, and no file I/O.  This is a native
   re-author: the decisive taint is obtained by directed graph reachability
   from the expression's live symbol vertices to provenance-tag vertices.
   The redundant ablation witness uses partial derivatives, not substitution.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
failureMessage = "";

localVerdict = "DENSITY_PORT_VECTOR_FREE";
jointPartial = "DENSITY_PORT_HOSTED (2/4, VECTOR-FREEDOM — N0_den is computationally vector-free; 024 = derivation, 026 = continuity lineage, 027 = port checks + closure)";
baselineTags = {
  "continuity_interface", "pathA_29_bulk", "pathA_32_wall"
};
taggedCarrierTag = "pathA_34_dimensionless_free_carrier";

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

clean[expr_] := Factor[Cancel[Together[expr]]];

externalSymbolNames = <|
  "cS" -> "c_s", "rhoEff" -> "rho_eff", "xiQ" -> "Xi_Q",
  "etaQ" -> "eta_q", "etaPhi" -> "eta_phi",
  "varpiQ2" -> "varpi_q2", "varpiPhi2" -> "varpi_Phi2",
  "lambdaC" -> "lambda_c", "phi2" -> "Phi2",
  "OmegaU" -> "Omega_U", "OmegaW" -> "Omega_W",
  "Rmix" -> "R_mix", "gU" -> "g_U", "gW" -> "g_W",
  "Aw" -> "A_w", "Fmuw" -> "F_muw",
  "omegaWall" -> "omega_wall", "omegaRho" -> "omega_rho",
  "rMix" -> "r_mix", "gRho" -> "g_rho", "gQold" -> "g_qold",
  "sigmaHidden" -> "sigma_hidden", "freeCarrier" -> "free_carrier",
  "missingRider" -> "missing_rider", "foreignSubject" -> "foreign_subject"
|>;

externalName[symbol_Symbol] := Module[{raw = SymbolName[symbol]},
  If[KeyExistsQ[externalSymbolNames, raw], externalSymbolNames[raw], raw]
];

fmt[expr_String] := expr;
fmt[expr_] := StringReplace[
  ToString[InputForm[clean[expr]]],
  Normal[externalSymbolNames]
];

fmtNames[items_List] := "{" <> StringRiffle[Sort[items], ", "] <> "}";

fail[] := Throw[failureMessage, "ledgerStage025Failure"];

expectZero[name_, residual_] := Module[{c},
  c = clean[residual];
  If[TrueQ[c === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    failureMessage = name <> ": residual = " <> fmt[c];
    Print["FAIL  ", failureMessage];
    fail[]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

(* Stage024 export symbols, coordinate metadata, and stage025 rig fixtures. *)
vectorSymbols = {Aw, Fmuw, Jw, U, W, OmegaU, OmegaW, Rmix, gU, gW};
relabelSymbols = {omegaWall, omegaRho, rMix, gRho, gQold};

hostContractNames = Sort[{
  "a", "c_s", "rho_eff", "I25", "Xi_Q", "eta_q", "eta_phi",
  "varpi_q2", "varpi_Phi2", "lambda_c"
}];
densityHostUniverseNames = Sort[Join[hostContractNames, {"q2", "Phi2"}]];
vectorHostNames = Sort[externalName /@ vectorSymbols];

allAuthoredSymbols = {
  a, cS, rhoEff, I25, xiQ, etaQ, etaPhi, varpiQ2, varpiPhi2,
  lambdaC, q2, phi2, Aw, Fmuw, Jw, U, W, OmegaU, OmegaW, Rmix,
  gU, gW, omegaWall, omegaRho, rMix, gRho, gQold, sigmaHidden,
  freeCarrier, missingRider, foreignSubject
};
externalSymbolByName = AssociationThread[externalName /@ allAuthoredSymbols, allAuthoredSymbols];

(* Reshaped pathA_43 provenance ledger; I25 is attached from momentValid. *)
baseSourceTags = <|
  "a" -> {"pathA_29_bulk", "pathA_32_wall"},
  "c_s" -> {"pathA_29_bulk"},
  "rho_eff" -> {"pathA_29_bulk"},
  "Xi_Q" -> {"continuity_interface"},
  "eta_q" -> {"continuity_interface"},
  "eta_phi" -> {"continuity_interface"},
  "varpi_q2" -> {"pathA_32_wall"},
  "varpi_Phi2" -> {"pathA_29_bulk"},
  "lambda_c" -> {
    "continuity_interface", "pathA_29_bulk", "pathA_32_wall"
  },
  "free_carrier" -> {},
  "sigma_hidden" -> {"vector_port"}
|>;
Scan[
  Function[symbol, AssociateTo[baseSourceTags, externalName[symbol] -> {"vector_port"}]],
  Join[vectorSymbols, relabelSymbols]
];

citedN0Den[] :=
  I25^2 xiQ^2 cS^4 rhoEff (etaPhi varpiQ2 + etaQ lambdaC)^2/
    (a^7 (lambdaC^2 - varpiPhi2 varpiQ2)^2);

sourceTagLedger[
  momentName_String,
  momentValid_,
  freeCarrierTagged_: False
] := Module[{ledger = Association[Normal[baseSourceTags]]},
  If[TrueQ[freeCarrierTagged],
    AssociateTo[ledger, "free_carrier" -> {taggedCarrierTag}]
  ];
  AssociateTo[
    ledger,
    momentName -> If[
      TrueQ[momentValid],
      {"continuity_interface", "pathA_32_wall"},
      {}
    ]
  ];
  ledger
];

(* Native provenance representation: explicit symbol -> tag directed graph. *)
makeProvenanceGraph[ledger_Association] := Module[
  {symbolVertices, tagVertices, edges},
  symbolVertices = ("sym:" <> #) & /@ Keys[ledger];
  tagVertices = ("tag:" <> #) & /@ DeleteDuplicates[Flatten[Values[ledger]]];
  edges = Flatten[KeyValueMap[
    Function[{name, tags},
      (DirectedEdge["sym:" <> name, "tag:" <> #]) & /@ tags
    ],
    ledger
  ]];
  Graph[Join[symbolVertices, tagVertices], edges]
];

makeProvenance[ledger_Association] := <|
  "Ledger" -> ledger,
  "Graph" -> makeProvenanceGraph[ledger]
|>;

globalSymbolNames[expr_] := Sort[DeleteDuplicates[
  externalName /@ Cases[
    HoldComplete[expr],
    symbol_Symbol /; Context[symbol] === "Global`",
    Infinity
  ]
]];

reachableTags[graph_Graph, name_String] := Module[{vertices},
  vertices = VertexOutComponent[graph, "sym:" <> name];
  Sort[StringDrop[#, 4] & /@ Select[vertices, StringStartsQ[#, "tag:"] &]]
];

(* The directed reachability result itself is verdict-bearing. *)
analyzeProvenance[expr_, provenance_Association] := Module[
  {ledger, graph, live, known, missing, present, empty, reached, taint,
   vectorHosts},
  ledger = provenance["Ledger"];
  graph = provenance["Graph"];
  live = globalSymbolNames[expr];
  known = Keys[ledger];
  missing = Sort[Complement[live, known]];
  present = Complement[live, missing];
  empty = Sort[Select[present, reachableTags[graph, #] === {} &]];
  reached = DeleteDuplicates[Flatten[
    VertexOutComponent[graph, "sym:" <> #] & /@ present
  ]];
  taint = Sort[StringDrop[#, 4] & /@
    Select[reached, StringStartsQ[#, "tag:"] &]];
  vectorHosts = Sort[Intersection[live, vectorHostNames]];
  <|
    "LiveSymbols" -> live,
    "TaintSet" -> taint,
    "MissingSourceSymbols" -> missing,
    "EmptySourceSymbols" -> empty,
    "VectorHostSymbols" -> vectorHosts
  |>
];

sourceMapCompleteQ[analysis_Association] :=
  analysis["MissingSourceSymbols"] === {} &&
  analysis["EmptySourceSymbols"] === {};

vectorFreeQ[analysis_Association] :=
  sourceMapCompleteQ[analysis] &&
  analysis["VectorHostSymbols"] === {} &&
  FreeQ[analysis["TaintSet"], "vector_port"];

(* Independent redundant witness: derivative independence, never rig routing. *)
derivativeAblationDelta[expr_, provenance_Association] := Module[
  {ledger, graph, taggedNames, candidateNames, liveNames, liveVectorNames,
   symbols, derivatives},
  ledger = provenance["Ledger"];
  graph = provenance["Graph"];
  taggedNames = Select[
    Keys[ledger],
    MemberQ[reachableTags[graph, #], "vector_port"] &
  ];
  candidateNames = DeleteDuplicates[Join[vectorHostNames, taggedNames]];
  liveNames = globalSymbolNames[expr];
  liveVectorNames = Intersection[liveNames, candidateNames];
  If[liveVectorNames === {}, Return[0]];
  symbols = externalSymbolByName[#] & /@ liveVectorNames;
  derivatives = clean[D[expr, #]] & /@ symbols;
  If[And @@ (TrueQ[# === 0] & /@ derivatives), 0, derivatives]
];

rigAssert[name_String, condition_] := If[
  ! TrueQ[condition],
  Throw[name, "stage025RigAssertion"]
];

probeAssertion[held_HoldComplete, expectedName_String] := Module[{fired},
  fired = Catch[
    ReleaseHold[held];
    "NO_ASSERT_FIRED",
    "stage025RigAssertion",
    Function[{message, tag}, message]
  ];
  {TrueQ[fired === expectedName], fired}
];

assertionPasses[held_HoldComplete] := TrueQ[Catch[
  ReleaseHold[held];
  True,
  "stage025RigAssertion",
  Function[{message, tag}, False]
]];

exerciseRig[
  label_String,
  assertionName_String,
  rigHeld_HoldComplete,
  neutralHeld_HoldComplete
] := Module[{probe, caught, firedName, neutralPass, outcome},
  probe = probeAssertion[rigHeld, assertionName];
  caught = probe[[1]];
  firedName = probe[[2]];
  neutralPass = assertionPasses[neutralHeld];
  expectBool[
    "META " <> label <> " routed assertion fires and neutering stops it",
    caught && neutralPass
  ];
  outcome = "CAUGHT at " <> ToString[firedName] <> "; neutralized=PASS";
  Print["RIG ", label, ": ", outcome];
  outcome
];

definitionArity[function_Symbol] := Module[{definitions, lhs},
  definitions = DownValues[function];
  If[Length[definitions] =!= 1, Return[-1]];
  lhs = Extract[definitions, {1, 1}, HoldComplete];
  Extract[
    lhs,
    {1, 1},
    Function[call, Length[Unevaluated[call]], HoldAllComplete]
  ]
];

heldCallArity[held_HoldComplete] := Extract[
  held,
  {1},
  Function[call, Length[Unevaluated[call]], HoldAllComplete]
];

arityScan[held_HoldComplete] :=
  definitionArity[analyzeProvenance] === heldCallArity[held];

leakageFreeQ[objects_] := FreeQ[
  objects,
  analyzeProvenance | sourceTagLedger | makeProvenanceGraph |
  makeProvenance | derivativeAblationDelta | authoredLeak | Solve |
  FullSimplify
];

runAudit[] := Module[
  {momentValid, n0Den, ledger, provenance, analysis, baselineAncestryOk,
   ablationDelta, outcomes, relabelP, relabelDelta, relabelExpr,
   relabelAnalysis, relabelNeuteredLedger, relabelNeutralAnalysis,
   hiddenExpr, hiddenAnalysis, hiddenNeuteredLedger, hiddenNeutralAnalysis,
   injectedAnalysis, riderExpr, riderAnalysis, taggedLedger,
   taggedProvenance, taggedAnalysis, missingExpr, missingAnalysis,
   repairedLedger, repairedAnalysis, expectedFourTags, fAssertName,
   fProbe, oldPort, oldPortAnalysis, properCall, badCall, actualTranscript,
   leakedTranscript, corruptSubject, localOK},

  momentValid = True; (* Typed forward reference; stage026 earns it. *)
  n0Den = clean[citedN0Den[]];
  ledger = sourceTagLedger["I25", momentValid, False];
  provenance = makeProvenance[ledger];
  analysis = analyzeProvenance[n0Den, provenance];

  subheading["I. cited stage024 subject-integrity contract"];
  expectBool[
    "I cited N0_den live-symbol contract equals stage024's exact 10-symbol export",
    analysis["LiveSymbols"] === hostContractNames
  ];
  Print[
    "DIAGNOSTIC (not counted): density-host metadata structural relationship: ",
    "DENSITY_HOST_UNIVERSE minus HOST_CONTRACT = ",
    fmtNames[Complement[densityHostUniverseNames, hostContractNames]],
    " (10-vs-12)"
  ];

  subheading["P1/P2. computed provenance-taint proof"];
  baselineAncestryOk = analysis["TaintSet"] === Sort[baselineTags];
  expectBool[
    "P1 baseline_ancestry_ok has exactly the three density tags",
    baselineAncestryOk
  ];
  expectBool[
    "P2 source_map_complete has no missing or empty provenance",
    sourceMapCompleteQ[analysis]
  ];
  expectBool[
    "P2 vector_host_symbols is empty",
    analysis["VectorHostSymbols"] === {}
  ];
  expectBool[
    "P2 computed taint excludes vector_port",
    FreeQ[analysis["TaintSet"], "vector_port"]
  ];
  expectBool[
    "P2 vector_free combines completeness, hosts, and computed taint",
    vectorFreeQ[analysis]
  ];
  Print[
    "DIAGNOSTIC (not counted): moment_valid=True is a typed forward reference ",
    "cited from stage026; the LOCAL verdict is conditional on it"
  ];

  (* Redundant witness only: deliberately not an expect* gate. *)
  ablationDelta = derivativeAblationDelta[n0Den, provenance];

  subheading["A-I. routed rig assertions and coupling meta-tests"];
  outcomes = <||>;

  relabelP = omegaWall^2 gRho + rMix gQold;
  relabelDelta = omegaWall^2 omegaRho^2 - rMix^2;
  relabelExpr = clean[relabelP^2/relabelDelta^2];
  relabelAnalysis = analyzeProvenance[relabelExpr, provenance];
  relabelNeuteredLedger = Association[Normal[ledger]];
  Scan[
    Function[symbol,
      AssociateTo[
        relabelNeuteredLedger,
        externalName[symbol] -> {"continuity_interface"}
      ]
    ],
    relabelSymbols
  ];
  relabelNeutralAnalysis = analyzeProvenance[
    relabelExpr, makeProvenance[relabelNeuteredLedger]
  ];
  outcomes["A"] = exerciseRig[
    "A relabel_rig",
    "A relabel_rig computed-taint assert: vector_port absent",
    HoldComplete[rigAssert[
      "A relabel_rig computed-taint assert: vector_port absent",
      FreeQ[relabelAnalysis["TaintSet"], "vector_port"]
    ]],
    HoldComplete[rigAssert[
      "A relabel_rig computed-taint assert: vector_port absent",
      FreeQ[relabelNeutralAnalysis["TaintSet"], "vector_port"]
    ]]
  ];

  hiddenExpr = clean[n0Den sigmaHidden];
  hiddenAnalysis = analyzeProvenance[hiddenExpr, provenance];
  hiddenNeuteredLedger = Association[Normal[ledger]];
  AssociateTo[
    hiddenNeuteredLedger,
    "sigma_hidden" -> {"continuity_interface"}
  ];
  hiddenNeutralAnalysis = analyzeProvenance[
    hiddenExpr, makeProvenance[hiddenNeuteredLedger]
  ];
  outcomes["B"] = exerciseRig[
    "B hidden_vector",
    "B hidden_vector computed-taint assert: vector_port absent",
    HoldComplete[rigAssert[
      "B hidden_vector computed-taint assert: vector_port absent",
      FreeQ[hiddenAnalysis["TaintSet"], "vector_port"]
    ]],
    HoldComplete[rigAssert[
      "B hidden_vector computed-taint assert: vector_port absent",
      FreeQ[hiddenNeutralAnalysis["TaintSet"], "vector_port"]
    ]]
  ];

  injectedAnalysis = analyzeProvenance[
    clean[n0Den OmegaU/OmegaW], provenance
  ];
  outcomes["C"] = exerciseRig[
    "C vector_injection",
    "C vector_injection vector-host assert: vector_host_symbols empty",
    HoldComplete[rigAssert[
      "C vector_injection vector-host assert: vector_host_symbols empty",
      injectedAnalysis["VectorHostSymbols"] === {}
    ]],
    HoldComplete[rigAssert[
      "C vector_injection vector-host assert: vector_host_symbols empty",
      analysis["VectorHostSymbols"] === {}
    ]]
  ];

  riderExpr = clean[n0Den freeCarrier];
  riderAnalysis = analyzeProvenance[riderExpr, provenance];
  taggedLedger = sourceTagLedger["I25", momentValid, True];
  taggedProvenance = makeProvenance[taggedLedger];
  taggedAnalysis = analyzeProvenance[riderExpr, taggedProvenance];
  outcomes["D"] = exerciseRig[
    "D provenance_less_rider",
    "D provenance_less_rider source_map_complete assert",
    HoldComplete[rigAssert[
      "D provenance_less_rider source_map_complete assert",
      sourceMapCompleteQ[riderAnalysis]
    ]],
    HoldComplete[rigAssert[
      "D provenance_less_rider source_map_complete assert",
      sourceMapCompleteQ[taggedAnalysis]
    ]]
  ];

  missingExpr = clean[n0Den missingRider];
  missingAnalysis = analyzeProvenance[missingExpr, provenance];
  repairedLedger = Association[Normal[ledger]];
  AssociateTo[
    repairedLedger,
    "missing_rider" -> {taggedCarrierTag}
  ];
  repairedAnalysis = analyzeProvenance[
    missingExpr, makeProvenance[repairedLedger]
  ];
  outcomes["E"] = exerciseRig[
    "E missing_symbol",
    "E missing_symbol source_map_complete assert",
    HoldComplete[rigAssert[
      "E missing_symbol source_map_complete assert",
      sourceMapCompleteQ[missingAnalysis]
    ]],
    HoldComplete[rigAssert[
      "E missing_symbol source_map_complete assert",
      sourceMapCompleteQ[repairedAnalysis]
    ]]
  ];

  expectedFourTags = Sort[Append[baselineTags, taggedCarrierTag]];
  fAssertName = "F tagged_carrier source_map_complete flip assert";
  fProbe = probeAssertion[
    HoldComplete[rigAssert[fAssertName, sourceMapCompleteQ[riderAnalysis]]],
    fAssertName
  ];
  expectBool[
    "F tagged_carrier passes P2 with four tags and flips when tag is neutered",
    vectorFreeQ[taggedAnalysis] &&
    taggedAnalysis["TaintSet"] === expectedFourTags &&
    taggedAnalysis["TaintSet"] =!= Sort[baselineTags] &&
    TrueQ[fProbe[[1]]]
  ];
  outcomes["F"] =
    "PASS P2; taint=" <> fmtNames[expectedFourTags] <>
    "; tag-neutered=FLIP at " <> ToString[fProbe[[2]]];
  Print["RIG F tagged_carrier: ", outcomes["F"]];

  oldPort = clean[OmegaU^2 gW + Rmix gU];
  oldPortAnalysis = analyzeProvenance[oldPort, provenance];
  outcomes["G"] = exerciseRig[
    "G raw_vector_port",
    "G raw_vector_port vector-host assert: vector_host_symbols empty",
    HoldComplete[rigAssert[
      "G raw_vector_port vector-host assert: vector_host_symbols empty",
      oldPortAnalysis["VectorHostSymbols"] === {}
    ]],
    HoldComplete[rigAssert[
      "G raw_vector_port vector-host assert: vector_host_symbols empty",
      analysis["VectorHostSymbols"] === {}
    ]]
  ];

  properCall = HoldComplete[analyzeProvenance[n0Den, provenance]];
  badCall = HoldComplete[analyzeProvenance[n0Den]];
  outcomes["HArity"] = exerciseRig[
    "H' arity_scanner",
    "H' definition/call arity scanner assert",
    HoldComplete[rigAssert[
      "H' definition/call arity scanner assert",
      arityScan[badCall]
    ]],
    HoldComplete[rigAssert[
      "H' definition/call arity scanner assert",
      arityScan[properCall]
    ]]
  ];

  actualTranscript = {
    n0Den, analysis, ablationDelta, Sort[analysis["TaintSet"]]
  };
  leakedTranscript = Append[actualTranscript, authoredLeak[n0Den]];
  outcomes["HLeak"] = exerciseRig[
    "H' leakage_scanner",
    "H' unevaluated-leakage transcript scanner assert",
    HoldComplete[rigAssert[
      "H' unevaluated-leakage transcript scanner assert",
      leakageFreeQ[leakedTranscript]
    ]],
    HoldComplete[rigAssert[
      "H' unevaluated-leakage transcript scanner assert",
      leakageFreeQ[actualTranscript]
    ]]
  ];

  corruptSubject = clean[n0Den /. rhoEff -> foreignSubject];
  outcomes["I"] = exerciseRig[
    "I subject_integrity",
    "I subject_integrity exact host-contract assert",
    HoldComplete[rigAssert[
      "I subject_integrity exact host-contract assert",
      globalSymbolNames[corruptSubject] === hostContractNames
    ]],
    HoldComplete[rigAssert[
      "I subject_integrity exact host-contract assert",
      globalSymbolNames[n0Den] === hostContractNames
    ]]
  ];

  expectBool[
    "H' actual transcript has no unevaluated authored-helper leakage",
    leakageFreeQ[actualTranscript]
  ];
  localOK = TrueQ[momentValid && baselineAncestryOk && vectorFreeQ[analysis]];
  expectBool["LOCAL DENSITY_PORT_VECTOR_FREE = P1 and P2", localOK];

  <|
    "N0Den" -> n0Den,
    "Analysis" -> analysis,
    "AblationDelta" -> ablationDelta,
    "AblationWitnessOK" -> TrueQ[ablationDelta === 0],
    "Outcomes" -> outcomes,
    "MomentValid" -> momentValid,
    "LocalOK" -> localOK
  |>
];

emit[data_Association] := Module[{analysis = data["Analysis"]},
  subheading["Computed density-port vector-freedom transcript"];
  Print["consumes: stage024 N0_den (cited canonical factored export)"];
  Print["N0_den (canonical factored): ", fmt[data["N0Den"]]];
  Print["N0_den live symbols: ", fmtNames[analysis["LiveSymbols"]]];
  Print["host contract (10): ", fmtNames[hostContractNames]];
  Print["allowable density-host universe (12): ", fmtNames[densityHostUniverseNames]];
  Print["computed taint set: ", fmtNames[analysis["TaintSet"]]];
  Print["missing_source_symbols: ", fmtNames[analysis["MissingSourceSymbols"]]];
  Print["empty_source_symbols: ", fmtNames[analysis["EmptySourceSymbols"]]];
  Print["vector_host_symbols: ", fmtNames[analysis["VectorHostSymbols"]]];
  Print["source_map_complete: ", sourceMapCompleteQ[analysis]];
  Print["baseline_ancestry_ok (P1): ", analysis["TaintSet"] === Sort[baselineTags]];
  Print["vector_free (P2): ", vectorFreeQ[analysis]];
  Print["ablation delta (redundant witness only): ", fmt[data["AblationDelta"]]];
  Print["ablation witness invariant: ", data["AblationWitnessOK"], " (de-counted)"];
  Print["moment_valid: True (typed forward reference; LOCAL verdict is conditional until stage026 discharges it)"];

  subheading["Verdict labels"];
  Print["LOCAL_AUDIT_VERDICT: ", localVerdict];
  Print["JOINT_LANDING_LABEL (PARTIAL): ", jointPartial]
];

runAll[] := Module[{data},
  heading["ledger_stage025_vector_freedom_taint_mathematica_audit"];
  Print["Target stem confirmed: ledger_stage025_vector_freedom_taint"];
  Print["Engine: Wolfram directed provenance-graph reachability; zero file I/O."];
  data = runAudit[];
  emit[data];
  0
];

result = Catch[
  runAll[],
  "ledgerStage025Failure",
  Function[{message, tag}, 1]
];

heading["Tallies"];
Print[
  "TALLY mathematica: ", passCount, " pass + ", failCount,
  " fail = ", passCount + failCount, " checks"
];
If[result === 0 && failCount === 0,
  Print["OVERALL PASS"];
  Exit[0],
  If[failureMessage =!= "", Print["ABORTED: ", failureMessage]];
  Print["OVERALL FAIL"];
  Exit[1]
];
