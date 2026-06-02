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

stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

cleanScalar[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanScalar[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

$Assumptions = True;

banner["STAGE 215 -- FULL PRIMITIVE-QUADRUPLE RANKING AND THE UP-TO-FOUR-COORDINATE SIEVE"];

subbanner["M1. Combinatorial ledger"];

axes = {lambda, c, gamma, U, W};
triples = Subsets[axes, {3}];
quadruples = Subsets[axes, {4}];

quadFaceCounts = Table[
  Count[triples, tri_ /; SubsetQ[quad, tri]],
  {quad, quadruples}
];
tripleIncidenceCounts = Table[
  Count[quadruples, quad_ /; SubsetQ[quad, tri]],
  {tri, triples}
];
axisIncidenceCounts = Table[
  Total[Boole[MemberQ[#, axis] & /@ quadruples]],
  {axis, axes}
];

Print["M1 axes = ", fmt[axes]];
Print["M1 triples = ", fmt[triples]];
Print["M1 quadruples = ", fmt[quadruples]];
Print["M1 quadruple face counts = ", fmt[quadFaceCounts]];
Print["M1 triple incidence counts = ", fmt[tripleIncidenceCounts]];
Print["M1 axis incidence counts = ", fmt[axisIncidenceCounts]];

expectZero["M1 #triples - Binomial[5,3]", Length[triples] - Binomial[5, 3]];
expectZero["M1 #quadruples - Binomial[5,4]", Length[quadruples] - Binomial[5, 4]];
expectTrue["M1 every quadruple has four triple faces", AllTrue[quadFaceCounts, # == 4 &]];
expectTrue["M1 every triple lies in two quadruples", AllTrue[tripleIncidenceCounts, # == 2 &]];
expectTrue["M1 every axis lies in four quadruples", AllTrue[axisIncidenceCounts, # == 4 &]];

subbanner["M2. Min flattening"];

m2LoVars = {iotaLo, aLo, bLo, cLo, dLo};
m2HiVars = {iotaHi, aHi, bHi, cHi, dHi};
loNested = Min[iotaLo, Min[aLo, bLo, cLo, dLo]];
loFlat = Min[iotaLo, aLo, bLo, cLo, dLo];
hiNested = Min[iotaHi, Min[aHi, bHi, cHi, dHi]];
hiFlat = Min[iotaHi, aHi, bHi, cHi, dHi];

Print["M2 lo nested Min = ", fmt[loNested]];
Print["M2 lo flat Min = ", fmt[loFlat]];
Print["M2 hi nested Min = ", fmt[hiNested]];
Print["M2 hi flat Min = ", fmt[hiFlat]];

expectTrue[
  "M2 lo envelope nested Min collapse",
  FullSimplify[loNested == loFlat, Assumptions -> Element[m2LoVars, Reals]]
];
expectTrue[
  "M2 hi envelope nested Min collapse",
  FullSimplify[hiNested == hiFlat, Assumptions -> Element[m2HiVars, Reals]]
];

subbanner["M3. Boundary-splice / full-simplex interval"];

m3Theorem = Resolve[
  ForAll[
    {betaLo, betaHi, iotaLo, iotaHi, bBest, iBest},
    Implies[
      betaLo <= betaHi && iotaLo <= iotaHi &&
        betaLo <= bBest <= betaHi && iotaLo <= iBest <= iotaHi,
      Min[betaLo, iotaLo] <= Min[bBest, iBest] <= Min[betaHi, iotaHi]
    ]
  ],
  Reals
];
expectTrue["M3 continuous quantified splice theorem", m3Theorem];

subbanner["M4. Local classification orderings"];

m4Interior = Resolve[
  ForAll[
    {betaLo, betaHi, iotaLo, iotaHi, bBest, iBest},
    Implies[
      betaLo <= betaHi && iotaLo <= iotaHi && iotaHi < betaLo &&
        betaLo <= bBest <= betaHi && iotaLo <= iBest <= iotaHi,
      iBest < bBest
    ]
  ],
  Reals
];
m4Boundary = Resolve[
  ForAll[
    {betaLo, betaHi, iotaLo, iotaHi, bBest, iBest},
    Implies[
      betaLo <= betaHi && iotaLo <= iotaHi && iotaLo > betaHi &&
        betaLo <= bBest <= betaHi && iotaLo <= iBest <= iotaHi,
      bBest < iBest
    ]
  ],
  Reals
];
expectTrue["M4 interior-certified order", m4Interior];
expectTrue["M4 boundary-certified order", m4Boundary];

subbanner["M5. Primitive-quadruple ranking and unique winner"];

m5Pairwise = Resolve[
  ForAll[
    {L1, U1, L2, U2, x, y},
    Implies[
      L1 <= U1 && L2 <= U2 && U1 < L2 &&
        L1 <= x <= U1 && L2 <= y <= U2,
      x < y
    ]
  ],
  Reals
];
expectTrue["M5 pairwise certified interval ranking", m5Pairwise];

lVars = {L1, L2, L3, L4, L5};
uVars = {U1, U2, U3, U4, U5};
xVars = {x1, x2, x3, x4, x5};
allIntervalVars = Join[lVars, uVars, xVars];

uniqueWinnerTheorem[star_Integer] := Module[
  {orderedIntervals, chosenValues, minOtherLower, winnerPremise, winnerConclusion},
  orderedIntervals = And @@ Thread[lVars <= uVars];
  chosenValues = And @@ Table[lVars[[p]] <= xVars[[p]] <= uVars[[p]], {p, 1, 5}];
  minOtherLower = Min @@ Delete[lVars, star];
  winnerPremise = uVars[[star]] < minOtherLower;
  winnerConclusion = And @@ Table[
    If[p == star, True, xVars[[star]] < xVars[[p]]],
    {p, 1, 5}
  ];
  With[
    {vars = allIntervalVars, premise = orderedIntervals && chosenValues && winnerPremise,
      conclusion = winnerConclusion},
    Resolve[ForAll[vars, Implies[premise, conclusion]], Reals]
  ]
];

uniqueWinnerResults = Table[uniqueWinnerTheorem[star], {star, 1, 5}];
Print["M5 unique-winner quantified results = ", fmt[uniqueWinnerResults]];
expectTrue[
  "M5 unique certified winner over five quadruple intervals",
  AllTrue[uniqueWinnerResults, TrueQ]
];

subbanner["M6. Global support<=4 splice and improvement tests"];

m6Splice = Resolve[
  ForAll[
    {tau3Lo, tau3Hi, tau4Lo, tau4Hi, sBest, qBest},
    Implies[
      tau3Lo <= tau3Hi && tau4Lo <= tau4Hi &&
        tau3Lo <= sBest <= tau3Hi && tau4Lo <= qBest <= tau4Hi,
      Min[tau3Lo, tau4Lo] <= Min[sBest, qBest] <= Min[tau3Hi, tau4Hi]
    ]
  ],
  Reals
];
m6Improvement = Resolve[
  ForAll[
    {tau3Lo, tau3Hi, tau4Lo, tau4Hi, sBest, qBest},
    Implies[
      tau3Lo <= tau3Hi && tau4Lo <= tau4Hi && tau4Hi < tau3Lo &&
        tau3Lo <= sBest <= tau3Hi && tau4Lo <= qBest <= tau4Hi,
      qBest < sBest
    ]
  ],
  Reals
];
m6NoImprovement = Resolve[
  ForAll[
    {tau3Lo, tau3Hi, tau4Lo, tau4Hi, sBest, qBest},
    Implies[
      tau3Lo <= tau3Hi && tau4Lo <= tau4Hi && tau4Lo > tau3Hi &&
        tau3Lo <= sBest <= tau3Hi && tau4Lo <= qBest <= tau4Hi,
      sBest < qBest
    ]
  ],
  Reals
];
expectTrue["M6 support<=4 splice theorem", m6Splice];
expectTrue["M6 certified four-coordinate improvement", m6Improvement];
expectTrue["M6 certified no-improvement", m6NoImprovement];

subbanner["M7. Finite budget reconstructed from factors"];

degreePattern = {3, 3, 3, 2};
supportLe3Factors = {{10, 12}, {10, 48}};
quadEvalPerEnvelope = Product[degreePattern[[idx]], {idx, 1, Length[degreePattern]}];
supportLe3Budget = Total[Times @@@ supportLe3Factors];
quadEvalTotal = Binomial[5, 4]*2*quadEvalPerEnvelope;
fullEvalTotal = supportLe3Budget + quadEvalTotal;

Print["M7 degree pattern = ", fmt[degreePattern]];
Print["M7 support<=3 factor products = ", fmt[supportLe3Factors]];
Print["M7 per-envelope candidate bound = ", fmt[quadEvalPerEnvelope]];
Print["M7 support<=3 budget = ", fmt[supportLe3Budget]];
Print["M7 interior quadruple budget = ", fmt[quadEvalTotal]];
Print["M7 full support<=4 budget = ", fmt[fullEvalTotal]];

expectZero["M7 per-envelope 3*3*3*2 factor product", quadEvalPerEnvelope - Times @@ degreePattern];
expectZero["M7 support<=3 10*12+10*48 factor sum", supportLe3Budget - Dot[{10, 10}, {12, 48}]];
expectZero["M7 interior quadruple budget - 540", quadEvalTotal - 540];
expectZero["M7 full support<=4 budget - 1140", fullEvalTotal - 1140];

Print[""];
Print["Stage 215 Mathematica audit passed."];

Exit[0];
