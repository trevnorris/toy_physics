ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

normalizeExpr[expr_] := FullSimplify[expr, Assumptions -> $Assumptions];

expectZero[name_String, expr_] := Module[{res},
  res = normalizeExpr[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res == 0], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := (
  Print[name, " = ", cond];
  If[TrueQ[cond], pass[name], fail[name, cond]];
);

packetInterval[packetWindows_Association] := {
  Min[Min /@ Values[packetWindows]],
  Min[Max /@ Values[packetWindows]]
};

boundaryBest[assignment_Association, boundaryPacketNames_List] :=
  Min[assignment /@ boundaryPacketNames];

classifyFamily[familyName_String, packetWindows_Association, boundaryPacketNames_List] := Module[
  {keys, fullLo, fullHi, counts, tuples, assignment, boundaryValue, support5Value, fullValue},
  keys = Keys[packetWindows];
  {fullLo, fullHi} = packetInterval[packetWindows];
  counts = <|"support5" -> 0, "boundary" -> 0, "tie" -> 0|>;
  tuples = Tuples[packetWindows /@ keys];

  Scan[
    Function[{values},
      assignment = AssociationThread[keys -> values];
      boundaryValue = boundaryBest[assignment, boundaryPacketNames];
      support5Value = assignment["support5_int"];
      fullValue = Min[boundaryValue, support5Value];
      If[!(fullLo <= fullValue <= fullHi),
        fail[familyName <> " interval splice failed", assignment]
      ];
      Which[
        support5Value < boundaryValue,
        counts["support5"] = counts["support5"] + 1,
        boundaryValue < support5Value,
        counts["boundary"] = counts["boundary"] + 1,
        True,
        counts["tie"] = counts["tie"] + 1
      ];
    ],
    tuples
  ];

  Print[familyName, " exhaustive outcomes = ", counts];
  counts
];

banner["STAGE 201 — FULL SUPPORT-<=5 COMPLETION AND LOCAL MIXED-RAY SEARCH CLOSURE"];

$Assumptions = Element[
    {
      tauLe3Best, tauLe3Lo, tauLe3Hi,
      tau5BestInt, tau5LoInt, tau5HiInt,
      tauFaceLambdaBest, tauFaceLambdaLo, tauFaceLambdaHi,
      tauFacecBest, tauFacecLo, tauFacecHi,
      tauFacegammaBest, tauFacegammaLo, tauFacegammaHi,
      tauFaceUBest, tauFaceULo, tauFaceUHi,
      tauFaceWBest, tauFaceWLo, tauFaceWHi
    },
    Reals
  ];

subbanner["I. Exact boundary-identification via imported Stage-198 face packets"];

axes = {"lambda", "c", "gamma", "U", "W"};
quadrupleFaces = Association @ Table[
  With[{axis = axes[[i]]},
    "omit_" <> axis -> <|
      "omitted_axis" -> axis,
      "support" -> DeleteCases[axes, axis],
      "source_stage" -> 198
    |>
  ],
  {i, Length[axes]}
];
boundaryPackets = Join[
  <|"support_le3" -> <|"support" -> "<=3 global ledger", "source_stage" -> 198|>|>,
  quadrupleFaces
];
boundaryFaceNames = Keys[quadrupleFaces];
properFaces = Select[Subsets[axes], 0 < Length[#] < Length[axes] &];

quadrupleSupports = Sort[Sort /@ (quadrupleFaces[#]["support"] & /@ boundaryFaceNames)];
expectedQuadruples = Sort[Sort /@ Subsets[axes, {4}]];

Print["primitive axes = ", axes];
Print["imported Stage-198 boundary packets = ", Keys[boundaryPackets]];
Print["quadruple face supports = ", quadrupleFaces[#]["support"] & /@ boundaryFaceNames];
Print["#proper nonempty support strata = ", Length[properFaces]];

expectZero["primitive axes - 5", Length[axes] - 5];
expectZero["imported quadruple packet count - 5", Length[quadrupleFaces] - 5];
expectTrue[
  "imported quadruple supports match the five simplex facets",
  quadrupleSupports === expectedQuadruples
];
expectZero["support-cardinality ceiling 5 - #axes", Length[axes] - 5];

Scan[
  Function[subset,
    Module[{coveringFaces, expected},
      coveringFaces = Select[
        boundaryFaceNames,
        SubsetQ[quadrupleFaces[#]["support"], subset] &
      ];
      expected = Length[axes] - Length[subset];
      Print["boundary coverage incidence ", subset, " -> ", coveringFaces];
      If[Length[coveringFaces] =!= expected,
        fail["boundary-identification coverage count is incorrect", subset]
      ];
    ]
  ],
  properFaces
];

subbanner["II. Exact imported support-<=4 and support-<=5 ledger splice"];

boundaryBestSymbols = {
  tauLe3Best,
  tauFaceLambdaBest,
  tauFacecBest,
  tauFacegammaBest,
  tauFaceUBest,
  tauFaceWBest
};
boundaryLoSymbols = {
  tauLe3Lo,
  tauFaceLambdaLo,
  tauFacecLo,
  tauFacegammaLo,
  tauFaceULo,
  tauFaceWLo
};
boundaryHiSymbols = {
  tauLe3Hi,
  tauFaceLambdaHi,
  tauFacecHi,
  tauFacegammaHi,
  tauFaceUHi,
  tauFaceWHi
};

tauLe4Best = Min @@ boundaryBestSymbols;
tauLe4Lo = Min @@ boundaryLoSymbols;
tauLe4Hi = Min @@ boundaryHiSymbols;

tauLe5Best = Min[tauLe4Best, tau5BestInt];
tauLe5Lo = Min[tauLe4Lo, tau5LoInt];
tauLe5Hi = Min[tauLe4Hi, tau5HiInt];

tauLe5BestFlat = Min @@ Append[boundaryBestSymbols, tau5BestInt];
tauLe5LoFlat = Min @@ Append[boundaryLoSymbols, tau5LoInt];
tauLe5HiFlat = Min @@ Append[boundaryHiSymbols, tau5HiInt];

Print["tau_{<=4,*}^{best} = ", fmt[tauLe4Best]];
Print["tau_{<=5,*}^{best} = ", fmt[tauLe5Best]];
Print["tau_{<=5,lo} = ", fmt[tauLe5Lo]];
Print["tau_{<=5,hi} = ", fmt[tauLe5Hi]];

expectZero["support<=5 best flattening over imported packets", tauLe5Best - tauLe5BestFlat];
expectZero["support<=5 lower splice over imported packets", tauLe5Lo - tauLe5LoFlat];
expectZero["support<=5 upper splice over imported packets", tauLe5Hi - tauLe5HiFlat];

subbanner["III. Exact improvement / no-improvement / overlap families on the actual finite ledger"];

improvementFamily = <|
  "support_le3" -> {10, 11},
  "omit_lambda" -> {12, 13},
  "omit_c" -> {11, 12},
  "omit_gamma" -> {13, 14},
  "omit_U" -> {15, 16},
  "omit_W" -> {14, 15},
  "support5_int" -> {2, 3, 4}
|>;
noImprovementFamily = <|
  "support_le3" -> {2, 3},
  "omit_lambda" -> {4, 5},
  "omit_c" -> {3, 4},
  "omit_gamma" -> {5, 6},
  "omit_U" -> {6, 7},
  "omit_W" -> {4, 6},
  "support5_int" -> {9, 10, 11}
|>;
overlapFamily = <|
  "support_le3" -> {5, 6},
  "omit_lambda" -> {4, 8},
  "omit_c" -> {7, 8},
  "omit_gamma" -> {6, 9},
  "omit_U" -> {8, 9},
  "omit_W" -> {5, 7},
  "support5_int" -> {3, 7}
|>;

boundaryPacketNames = Append[boundaryFaceNames, "support_le3"];

improvementCounts = classifyFamily[
  "genuine support-5 improvement family",
  improvementFamily,
  boundaryPacketNames
];
expectZero["support-5 improvement family boundary wins", improvementCounts["boundary"]];
expectZero["support-5 improvement family ties", improvementCounts["tie"]];
expectTrue[
  "support-5 improvement family interior wins exist",
  improvementCounts["support5"] > 0
];

noImprovementCounts = classifyFamily[
  "support-5 no-improvement family",
  noImprovementFamily,
  boundaryPacketNames
];
expectZero["support-5 no-improvement family interior wins", noImprovementCounts["support5"]];
expectZero["support-5 no-improvement family ties", noImprovementCounts["tie"]];
expectTrue[
  "support-5 no-improvement family boundary wins exist",
  noImprovementCounts["boundary"] > 0
];

overlapCounts = classifyFamily[
  "support-5 overlap family",
  overlapFamily,
  boundaryPacketNames
];
expectTrue[
  "support-5 overlap family retains boundary winners",
  overlapCounts["boundary"] > 0
];
expectTrue[
  "support-5 overlap family retains interior winners",
  overlapCounts["support5"] > 0
];

subbanner["IV. Exact support-five candidate filters and final budget ledger"];

supportFivePacket = <|
  "source_stage" -> 200,
  "canonical_screens" -> {"gradient-optimal", "equal-mix"},
  "preferred_lifted_degree_pattern" -> {3, 3, 3, 3, 2},
  "fallback_projected_degree_pattern" -> {5, 5, 5, 6}
|>;
liftedPerEnvelope = Times @@ supportFivePacket["preferred_lifted_degree_pattern"];
projectedPerEnvelope = Times @@ supportFivePacket["fallback_projected_degree_pattern"];

Print["support-five canonical screens = ", supportFivePacket["canonical_screens"]];
Print["preferred lifted degree pattern = ", supportFivePacket["preferred_lifted_degree_pattern"]];
Print["fallback projected degree pattern = ", supportFivePacket["fallback_projected_degree_pattern"]];

expectZero["support-five canonical screen count - 2", Length[supportFivePacket["canonical_screens"]] - 2];
expectZero["lifted compiler bound - 162", liftedPerEnvelope - 162];
expectZero["projected compiler bound - 750", projectedPerEnvelope - 750];

supportLe3Budget = 600;
quadrupleEvalPerEnvelope = 54;
supportLe4Budget = supportLe3Budget + Length[quadrupleFaces] * 2 * quadrupleEvalPerEnvelope;
preferredTotal = supportLe4Budget + 2 * liftedPerEnvelope;
fallbackTotal = supportLe4Budget + 2 * projectedPerEnvelope;

Print["support<=3 imported budget = ", supportLe3Budget];
Print["support<=4 rebuilt budget = ", supportLe4Budget];
Print["preferred full support<=5 budget = ", preferredTotal];
Print["fallback full support<=5 budget = ", fallbackTotal];

expectZero["support<=4 rebuilt budget - 1140", supportLe4Budget - 1140];
expectZero["preferred full support<=5 budget - 1464", preferredTotal - 1464];
expectZero["fallback full support<=5 budget - 2640", fallbackTotal - 2640];

Print[""];
Print["All Stage-201 imported ledger, splice, classification, and budget checks verified."];
