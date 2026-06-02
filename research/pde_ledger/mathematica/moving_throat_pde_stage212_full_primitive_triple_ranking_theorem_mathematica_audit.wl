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
  If[!MissingQ[detail], Print["  result -> ", fmt[detail]]];
  Exit[1];
);

expectExact[name_String, result_, expected_] := (
  Print[name, " = ", fmt[result]];
  If[TrueQ[result === expected], pass[name], fail[name, result]];
);

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === True], pass[name], fail[name, res]];
];

banner["STAGE 212 -- FULL PRIMITIVE-TRIPLE RANKING AND THE UP-TO-THREE-COORDINATE SIEVE"];

subbanner["M1. Combinatorial ledger"];

axes = {lambdaAxis, cAxis, gammaAxis, UAxis, WAxis};
pairs = Subsets[axes, {2}];
triples = Subsets[axes, {3}];

pairIncidenceCounts = Table[
  Count[triples, t_ /; SubsetQ[t, pair]],
  {pair, pairs}
];
axisIncidenceCounts = Table[
  Count[triples, t_ /; MemberQ[t, axis]],
  {axis, axes}
];

Print["M1 axes = ", fmt[axes]];
Print["M1 primitive pairs = ", fmt[pairs]];
Print["M1 primitive triples = ", fmt[triples]];
Print["M1 pair incidence counts = ", fmt[pairIncidenceCounts]];
Print["M1 axis incidence counts = ", fmt[axisIncidenceCounts]];

expectExact["M1 #pairs", Length[pairs], Binomial[5, 2]];
expectExact["M1 #triples", Length[triples], Binomial[5, 3]];
expectTrue["M1 every primitive pair lies in three triples", AllTrue[pairIncidenceCounts, # == 3 &]];
expectTrue["M1 every primitive axis lies in six triples", AllTrue[axisIncidenceCounts, # == 6 &]];

subbanner["M2. Boundary-splice nested Min identity"];

m2LoIdentity = Min[iotaLo, Min[tauIJLo, tauIKLo, tauJKLo]] ===
  Min[iotaLo, tauIJLo, tauIKLo, tauJKLo];
m2HiIdentity = Min[iotaHi, Min[tauIJHi, tauIKHi, tauJKHi]] ===
  Min[iotaHi, tauIJHi, tauIKHi, tauJKHi];

Print["M2 lo nested Min identity = ", fmt[m2LoIdentity]];
If[m2LoIdentity =!= True, fail["M2 lo nested Min identity", m2LoIdentity]];
pass["M2 lo nested Min identity"];

Print["M2 hi nested Min identity = ", fmt[m2HiIdentity]];
If[m2HiIdentity =!= True, fail["M2 hi nested Min identity", m2HiIdentity]];
pass["M2 hi nested Min identity"];

subbanner["M3. Symbolic interval and order theorems"];

m3a = Resolve[
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
Print["M3a local full-simplex sandwich = ", fmt[m3a]];
If[m3a =!= True, fail["M3a local full-simplex sandwich", m3a]];
pass["M3a local full-simplex sandwich"];

m3bInterior = Resolve[
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
Print["M3b interior-certified order = ", fmt[m3bInterior]];
If[m3bInterior =!= True, fail["M3b interior-certified order", m3bInterior]];
pass["M3b interior-certified order"];

m3bBoundary = Resolve[
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
Print["M3b boundary-certified order = ", fmt[m3bBoundary]];
If[m3bBoundary =!= True, fail["M3b boundary-certified order", m3bBoundary]];
pass["M3b boundary-certified order"];

m3c = Resolve[
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
Print["M3c primitive-triple ranking = ", fmt[m3c]];
If[m3c =!= True, fail["M3c primitive-triple ranking", m3c]];
pass["M3c primitive-triple ranking"];

m3dSandwich = Resolve[
  ForAll[
    {pLo, pHi, tLo, tHi, pBest, tBest},
    Implies[
      pLo <= pHi && tLo <= tHi &&
        pLo <= pBest <= pHi && tLo <= tBest <= tHi,
      Min[pLo, tLo] <= Min[pBest, tBest] <= Min[pHi, tHi]
    ]
  ],
  Reals
];
Print["M3d support<=3 splice sandwich = ", fmt[m3dSandwich]];
If[m3dSandwich =!= True, fail["M3d support<=3 splice sandwich", m3dSandwich]];
pass["M3d support<=3 splice sandwich"];

m3dImprovement = Resolve[
  ForAll[
    {pLo, pHi, tLo, tHi, pBest, tBest},
    Implies[
      pLo <= pHi && tLo <= tHi && tHi < pLo &&
        pLo <= pBest <= pHi && tLo <= tBest <= tHi,
      tBest < pBest
    ]
  ],
  Reals
];
Print["M3d three-coordinate improvement = ", fmt[m3dImprovement]];
If[m3dImprovement =!= True, fail["M3d three-coordinate improvement", m3dImprovement]];
pass["M3d three-coordinate improvement"];

m3dNoImprovement = Resolve[
  ForAll[
    {pLo, pHi, tLo, tHi, pBest, tBest},
    Implies[
      pLo <= pHi && tLo <= tHi && tLo > pHi &&
        pLo <= pBest <= pHi && tLo <= tBest <= tHi,
      pBest < tBest
    ]
  ],
  Reals
];
Print["M3d no-improvement order = ", fmt[m3dNoImprovement]];
If[m3dNoImprovement =!= True, fail["M3d no-improvement order", m3dNoImprovement]];
pass["M3d no-improvement order"];

subbanner["M4. Finite evaluation budget"];

pairwiseTotal = 10*12;
tripleInteriorTotal = 10*48;
fullTotal = pairwiseTotal + tripleInteriorTotal;

Print["M4 pairwise total = ", fmt[pairwiseTotal]];
Print["M4 triple interior total = ", fmt[tripleInteriorTotal]];
Print["M4 full support<=3 total = ", fmt[fullTotal]];

expectExact["M4 10*12", pairwiseTotal, 120];
expectExact["M4 10*48", tripleInteriorTotal, 480];
expectExact["M4 120+480", fullTotal, 600];

Print[""];
Print["Stage 212 Mathematica audit passed."];

Exit[0];
