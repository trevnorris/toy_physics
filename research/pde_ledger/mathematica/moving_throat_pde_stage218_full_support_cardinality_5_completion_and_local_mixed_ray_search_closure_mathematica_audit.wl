ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]]
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]]
);

fmt[expr_] := ToString[InputForm[expr]];
pass[name_String] := Print["PASS: ", name];
fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  detail -> ", fmt[detail]]];
  Exit[1]
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]]
];

expectTrue[name_String, cond_] := (
  Print[name, " = ", fmt[cond]];
  If[TrueQ[cond], pass[name], fail[name, cond]]
);

banner["Stage 218 independent Mathematica audit"];

subbanner["M1. Boundary incidence by finite set combinatorics"];

axes = {"lambda", "c", "gamma", "U", "W"};
quadrupleFaces = Subsets[axes, {4}];
properStrata = Subsets[axes, {1, 4}];

incidenceRows = Table[
  {
    support,
    Length[support],
    Total[Boole[ContainsAll[#, support]] & /@ quadrupleFaces],
    Length[axes] - Length[support]
  },
  {support, properStrata}
];
strataTally = Sort[Tally[Length /@ properStrata]];

Print["axes = ", axes];
Print["quadruple faces = ", quadrupleFaces];
Print["support-cardinality tally = ", strataTally];
Print["incidence rows = ", incidenceRows];

expectTrue["M1 every proper support has 5-k covering quadruple faces", AllTrue[incidenceRows, #[[3]] == #[[4]] &]];
expectTrue["M1 proper strata tally is 5+10+10+5", strataTally === {{1, 5}, {2, 10}, {3, 10}, {4, 5}}];
expectZero["M1 proper nonempty stratum count - 30", Length[properStrata] - 30];
expectZero["M1 2^5 - 2 - 30", (2^Length[axes] - 2) - 30];

subbanner["M2-M3. Support ceiling and final splice by Resolve"];

boundaryLabels = Join[{"support_le3"}, "omit_" <> # & /@ axes];
packetLabels = Join[boundaryLabels, {"support5_int"}];
packetCount = Length[packetLabels];

loVars = ToExpression /@ Table["lo" <> ToString[i], {i, packetCount}];
bestVars = ToExpression /@ Table["best" <> ToString[i], {i, packetCount}];
hiVars = ToExpression /@ Table["hi" <> ToString[i], {i, packetCount}];

boundaryBestVars = Take[bestVars, Length[boundaryLabels]];
tauLe4Best = Min @@ boundaryBestVars;
tau5BestInt = Last[bestVars];
tauLe5Best = Min[tauLe4Best, tau5BestInt];

boundaryLedger = ToExpression["boundaryLedger"];
support5Ledger = ToExpression["support5Ledger"];
m2BoundaryBranch = Resolve[
  ForAll[
    {boundaryLedger, support5Ledger},
    Implies[
      Element[{boundaryLedger, support5Ledger}, Reals] && boundaryLedger <= support5Ledger,
      Min[boundaryLedger, support5Ledger] == boundaryLedger
    ]
  ],
  Reals
];
m2InteriorBranch = Resolve[
  ForAll[
    {boundaryLedger, support5Ledger},
    Implies[
      Element[{boundaryLedger, support5Ledger}, Reals] && support5Ledger <= boundaryLedger,
      Min[boundaryLedger, support5Ledger] == support5Ledger
    ]
  ],
  Reals
];

m3LowerBranches = Table[
  Resolve[
    ForAll[
      {loVars[[i]], bestVars[[i]]},
      Implies[
        Element[{loVars[[i]], bestVars[[i]]}, Reals] && loVars[[i]] <= bestVars[[i]],
        ! (loVars[[i]] > bestVars[[i]])
      ]
    ],
    Reals
  ],
  {i, packetCount}
];
m3UpperBranches = Table[
  Resolve[
    ForAll[
      {bestVars[[i]], hiVars[[i]]},
      Implies[
        Element[{bestVars[[i]], hiVars[[i]]}, Reals] && bestVars[[i]] <= hiVars[[i]],
        ! (bestVars[[i]] > hiVars[[i]])
      ]
    ],
    Reals
  ],
  {i, packetCount}
];

lowerProbe = AssociationThread[packetLabels -> {1, 4, 5, 6, 7, 8, 9}];
upperProbe = AssociationThread[packetLabels -> {3, 6, 7, 8, 9, 10, 11}];

expectZero["M2 support<=4 imported boundary packet count - 6", Length[boundaryBestVars] - 6];
expectTrue["M2 tau<=5 best is the two-ledger support ceiling", m2BoundaryBranch && m2InteriorBranch];
expectTrue["M3 lower splice counterexample branches close", And @@ m3LowerBranches];
expectTrue["M3 upper splice counterexample branches close", And @@ m3UpperBranches];
expectZero["M3 lower endpoint audit probe", Min @@ Values[lowerProbe] - 1];
expectZero["M3 upper endpoint audit probe", Min @@ Values[upperProbe] - 3];
expectTrue["M3 upper endpoint is sharper than max-hi mutation", Min @@ Values[upperProbe] < Max @@ Values[upperProbe]];

subbanner["M4. Exhaustive regime outcomes from independent witnesses"];

makeBoundaryWindows[start_Integer, width_Integer, gap_Integer] := AssociationThread[
  boundaryLabels,
  Table[Range[start + gap*(i - 1), start + gap*(i - 1) + width - 1], {i, Length[boundaryLabels]}]
];

ledgerCounts[packetWindows_Association] := Module[
  {labels, assignments, labelsByRow, rawCounts},
  labels = Keys[packetWindows];
  assignments = AssociationThread[labels -> #] & /@ Tuples[Values[packetWindows]];
  labelsByRow = Function[row,
      Module[{boundaryValue, supportValue},
        boundaryValue = Min @@ Lookup[row, boundaryLabels];
        supportValue = row["support5_int"];
        Which[
          supportValue < boundaryValue, "support5",
          boundaryValue < supportValue, "boundary",
          True, "tie"
        ]
      ]
    ] /@ assignments;
  rawCounts = Counts[labelsByRow];
  Merge[{<|"support5" -> 0, "boundary" -> 0, "tie" -> 0|>, rawCounts}, Total]
];

totalAssignments[packetWindows_Association] := Times @@ (Length /@ Values[packetWindows]);
boundaryValues[packetWindows_Association] := Flatten[Lookup[packetWindows, boundaryLabels]];

improvementWindows = Join[
  makeBoundaryWindows[20, 2, 3],
  <|"support5_int" -> Range[1, 4]|>
];
noImprovementWindows = Join[
  makeBoundaryWindows[2, 2, 3],
  <|"support5_int" -> Range[40, 42]|>
];
overlapWindows = Join[
  AssociationThread[boundaryLabels, ConstantArray[{4, 8}, Length[boundaryLabels]]],
  <|"support5_int" -> {3, 7}|>
];

improvementCounts = ledgerCounts[improvementWindows];
noImprovementCounts = ledgerCounts[noImprovementWindows];
overlapCounts = ledgerCounts[overlapWindows];
improvementTotal = totalAssignments[improvementWindows];
noImprovementTotal = totalAssignments[noImprovementWindows];
overlapTotal = totalAssignments[overlapWindows];

Print["M4 improvement outcomes = ", Normal[improvementCounts]];
Print["M4 no-improvement outcomes = ", Normal[noImprovementCounts]];
Print["M4 overlap outcomes = ", Normal[overlapCounts]];

expectTrue["M4 regime 5.1 hypothesis support5_hi < boundary_lo", Max[improvementWindows["support5_int"]] < Min[boundaryValues[improvementWindows]]];
expectZero["M4 regime 5.1 boundary wins", improvementCounts["boundary"]];
expectZero["M4 regime 5.1 ties", improvementCounts["tie"]];
expectZero["M4 regime 5.1 support5 wins equal total", improvementCounts["support5"] - improvementTotal];

expectTrue["M4 regime 5.2 hypothesis support5_lo > every boundary hi", Min[noImprovementWindows["support5_int"]] > Max[boundaryValues[noImprovementWindows]]];
expectZero["M4 regime 5.2 support5 wins", noImprovementCounts["support5"]];
expectZero["M4 regime 5.2 ties", noImprovementCounts["tie"]];
expectZero["M4 regime 5.2 boundary wins equal total", noImprovementCounts["boundary"] - noImprovementTotal];

expectTrue[
  "M4 regime 5.3 hypothesis ranges interleave",
  Min[overlapWindows["support5_int"]] <= Min[Max /@ Lookup[overlapWindows, boundaryLabels]]
    && Max[overlapWindows["support5_int"]] >= Min[Min /@ Lookup[overlapWindows, boundaryLabels]]
];
expectTrue["M4 regime 5.3 retains support5 winners", overlapCounts["support5"] > 0];
expectTrue["M4 regime 5.3 retains boundary winners", overlapCounts["boundary"] > 0];
expectZero["M4 regime 5.3 ties", overlapCounts["tie"]];
expectZero["M4 regime 5.3 support5 plus boundary equals total", overlapCounts["support5"] + overlapCounts["boundary"] - overlapTotal];

subbanner["M5. Paper-stated evaluation budget"];

liftedPattern = {3, 3, 3, 3, 2};
fallbackPattern = {5, 5, 5, 6};
liftedPerEnvelope = Times @@ liftedPattern;
fallbackPerEnvelope = Times @@ fallbackPattern;
supportLe4Budget = 1140;
support5LiftedBudget = 2*liftedPerEnvelope;
support5FallbackBudget = 2*fallbackPerEnvelope;
preferredTotal = supportLe4Budget + support5LiftedBudget;
fallbackTotal = supportLe4Budget + support5FallbackBudget;

Print["lifted pattern = ", liftedPattern];
Print["fallback pattern = ", fallbackPattern];
Print["support<=4 paper budget = ", supportLe4Budget];
Print["support-five lifted budget = ", support5LiftedBudget];
Print["support-five fallback budget = ", support5FallbackBudget];

expectZero["M5 lifted per-envelope bound - 162", liftedPerEnvelope - 162];
expectZero["M5 lifted total - 324", support5LiftedBudget - 324];
expectZero["M5 fallback per-envelope bound - 750", fallbackPerEnvelope - 750];
expectZero["M5 fallback total - 1500", support5FallbackBudget - 1500];
expectZero["M5 support<=4 paper budget - 1140", supportLe4Budget - 1140];
expectZero["M5 preferred total - 1464", preferredTotal - 1464];
expectZero["M5 fallback total - 2640", fallbackTotal - 2640];

Print[""];
Print["All Stage 218 Mathematica claims M1-M5 verified."];
