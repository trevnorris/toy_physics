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

stripCE[expr_] := expr /. ConditionalExpression[e_, _] :> e;

cleanScalar[expr_] := FullSimplify[
  Together[Expand[stripCE[expr]]],
  Assumptions -> $Assumptions
];

cleanLogScalar[expr_] := FullSimplify[
  Together[Expand[PowerExpand[stripCE[expr]]]],
  Assumptions -> $Assumptions
];

cleanTensor[expr_] := If[
  ListQ[expr],
  Map[cleanScalar, expr, {Length[Dimensions[expr]]}],
  cleanScalar[expr]
];

cleanLogTensor[expr_] := If[
  ListQ[expr],
  Map[cleanLogScalar, expr, {Length[Dimensions[expr]]}],
  cleanLogScalar[expr]
];

zeroTensorQ[expr_] := If[
  ListQ[expr],
  And @@ (TrueQ[# === 0] & /@ Flatten[expr]),
  TrueQ[expr === 0]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanTensor[expr];
  Print[name, " = ", fmt[res]];
  If[zeroTensorQ[res], pass[name], fail[name, res]];
];

expectZeroLog[name_String, expr_] := Module[{res},
  res = cleanLogTensor[expr];
  Print[name, " = ", fmt[res]];
  If[zeroTensorQ[res], pass[name], fail[name, res]];
];

expectTrue[name_String, claim_] := Module[{res},
  res = FullSimplify[stripCE[claim], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 191 -- MINIMAL PDE DATA PACKET MATHEMATICA AUDIT"];

Clear[
  w,
  d0, d2, d4, n0, n2, n4,
  x20, x21, x22, xMean, xAnom, xSplit,
  d200, d202, d204, n200, n202, n204,
  d210, d212, d214, n210, n212, n214,
  d220, d222, d224, n220, n222, n224,
  dIso0, dIso2, dIso4, nIso0, nIso2, nIso4,
  grav, cSound, scaleA, cLight, mhat,
  chi0Star, fStar, mT, mK, mMu, rTr, rNt, rEta,
  qTr, qNt, qEta, ellT, ellK, ellMu
];

$Assumptions = (
  Element[
    {
      w,
      d0, d2, d4, n0, n2, n4,
      x20, x21, x22, xMean, xAnom, xSplit,
      d200, d202, d204, n200, n202, n204,
      d210, d212, d214, n210, n212, n214,
      d220, d222, d224, n220, n222, n224,
      dIso0, dIso2, dIso4, nIso0, nIso2, nIso4,
      grav, cSound, scaleA, cLight, mhat,
      chi0Star, fStar, mT, mK, mMu, rTr, rNt, rEta,
      qTr, qNt, qEta, ellT, ellK, ellMu
    },
    Reals
  ] &&
  d0 != 0 && d200 != 0 && d210 != 0 && d220 != 0 && dIso0 != 0 &&
  grav > 0 && cSound > 0 && scaleA > 0 && cLight > 0 && mhat > 0 &&
  chi0Star > 0 && fStar > 0 &&
  mT > 0 && mK > 0 && mMu > 0 && rTr > 0 && rNt > 0 && rEta > 0
);

responsePair[denData_] := Module[{shape},
  shape = denData[[1]]/(denData[[1]] + denData[[2]]*w^2 + denData[[3]]*w^4);
  cleanTensor[
    {
      Coefficient[Normal[Series[shape, {w, 0, 4}]], w, 2],
      Coefficient[Normal[Series[shape, {w, 0, 4}]], w, 4]
    }
  ]
];

prefactorTriple[laneData_] := Module[{den, num, pref},
  den = laneData[[1]] + laneData[[2]]*w^2 + laneData[[3]]*w^4;
  num = laneData[[4]] + laneData[[5]]*w^2 + laneData[[6]]*w^4;
  pref = laneData[[1]]*num/den^2;
  cleanTensor[
    {
      Coefficient[Normal[Series[pref, {w, 0, 4}]], w, 0],
      Coefficient[Normal[Series[pref, {w, 0, 4}]], w, 2],
      Coefficient[Normal[Series[pref, {w, 0, 4}]], w, 4]
    }
  ]
];

closedResponsePair[denData_] := cleanTensor[
  {
    -denData[[2]]/denData[[1]],
    (denData[[2]]^2 - denData[[1]]*denData[[3]])/denData[[1]]^2
  }
];

closedPrefactorTriple[laneData_] := cleanTensor[
  {
    laneData[[4]]/laneData[[1]],
    (laneData[[1]]*laneData[[5]] - 2*laneData[[2]]*laneData[[4]])/laneData[[1]]^2,
    (
      laneData[[1]]^2*laneData[[6]]
      - 2*laneData[[1]]*(laneData[[2]]*laneData[[5]] + laneData[[3]]*laneData[[4]])
      + 3*laneData[[2]]^2*laneData[[4]]
    )/laneData[[1]]^3
  }
];

traceVec = {1, 1, 1};
anomVec = {4, -1, -1};
splitVec = {0, 1, -1};
groupMetric = DiagonalMatrix[{1, 2, 2}];

coordinatesFromBasis[triple_] := Module[{sol},
  sol = First[
    Solve[
      Thread[triple == xMean*traceVec + xAnom*anomVec + xSplit*splitVec],
      {xMean, xAnom, xSplit},
      Reals
    ]
  ];
  cleanTensor[{xMean, xAnom, xSplit} /. sol]
];

sympyGroupCoords[triple_] := cleanTensor[
  {
    (triple[[1]] + 2*triple[[2]] + 2*triple[[3]])/5,
    (2*triple[[1]] - triple[[2]] - triple[[3]])/10,
    (triple[[2]] - triple[[3]])/2
  }
];

weightedProjector[v_] := cleanTensor[
  Outer[Times, v, v.groupMetric]/(v.groupMetric.v)
];

subbanner["I. Native Taylor coefficient compilers"];

nativeResponse = responsePair[{d0, d2, d4}];
targetResponse = closedResponsePair[{d0, d2, d4}];
nativePrefactor = prefactorTriple[{d0, d2, d4, n0, n2, n4}];
targetPrefactor = closedPrefactorTriple[{d0, d2, d4, n0, n2, n4}];
singleResponse = d0/(d0 + d2*w^2 + d4*w^4);
singlePrefactor = d0*(n0 + n2*w^2 + n4*w^4)/(d0 + d2*w^2 + d4*w^4)^2;
responseDerivativeRoute = cleanTensor[
  {
    (D[singleResponse, {w, 2}]/2) /. w -> 0,
    (D[singleResponse, {w, 4}]/24) /. w -> 0
  }
];
prefactorDerivativeRoute = cleanTensor[
  {
    singlePrefactor /. w -> 0,
    (D[singlePrefactor, {w, 2}]/2) /. w -> 0,
    (D[singlePrefactor, {w, 4}]/24) /. w -> 0
  }
];
oneLanePoleDefect = cleanTensor[nativeResponse[[2]] - 4*nativeResponse[[1]]^2];

Print["response coefficients from Series/Coefficient = ", fmt[nativeResponse]];
Print["prefactor coefficients from Series/Coefficient = ", fmt[nativePrefactor]];
expectZero["response Series/Coefficient - derivative route", nativeResponse - responseDerivativeRoute];
expectZero["prefactor Series/Coefficient - derivative route", nativePrefactor - prefactorDerivativeRoute];
expectZero["response coefficients - SymPy compiler", nativeResponse - targetResponse];
expectZero["prefactor coefficients - SymPy compiler", nativePrefactor - targetPrefactor];
expectZero[
  "one-lane pole defect + (3 d2^2 + d0 d4)/d0^2",
  oneLanePoleDefect + (3*d2^2 + d0*d4)/d0^2
];

subbanner["II. Weighted grouped trace/anomaly basis and projectors"];

groupNative = coordinatesFromBasis[{x20, x21, x22}];
groupTarget = sympyGroupCoords[{x20, x21, x22}];
reconstructedTriple = cleanTensor[
  groupNative[[1]]*traceVec + groupNative[[2]]*anomVec + groupNative[[3]]*splitVec
];

projTrace = weightedProjector[traceVec];
projAnom = weightedProjector[anomVec];
projSplit = weightedProjector[splitVec];

projTraceTarget = {{1, 2, 2}, {1, 2, 2}, {1, 2, 2}}/5;
projAnomTarget = {{16, -8, -8}, {-4, 2, 2}, {-4, 2, 2}}/20;
projSplitTarget = {{0, 0, 0}, {0, 2, -2}, {0, -2, 2}}/4;

Print["basis coordinates from Solve = ", fmt[groupNative]];
expectZero["group coordinates - SymPy formulas", groupNative - groupTarget];
expectZero["basis reconstruction - original grouped triple", reconstructedTriple - {x20, x21, x22}];
expectZero["trace basis G-orthogonal to anomaly basis", traceVec.groupMetric.anomVec];
expectZero["trace basis G-orthogonal to split basis", traceVec.groupMetric.splitVec];
expectZero["anomaly basis G-orthogonal to split basis", anomVec.groupMetric.splitVec];
expectZero["P_trace - SymPy projector", projTrace - projTraceTarget];
expectZero["P_anomaly - SymPy projector", projAnom - projAnomTarget];
expectZero["P_split - SymPy projector", projSplit - projSplitTarget];
expectZero["P_trace + P_anomaly + P_split - I", projTrace + projAnom + projSplit - IdentityMatrix[3]];
expectZero["P_trace^2 - P_trace", projTrace.projTrace - projTrace];
expectZero["P_anomaly^2 - P_anomaly", projAnom.projAnom - projAnom];
expectZero["P_split^2 - P_split", projSplit.projSplit - projSplit];
expectZero["P_trace P_anomaly", projTrace.projAnom];
expectZero["P_trace P_split", projTrace.projSplit];
expectZero["P_anomaly P_split", projAnom.projSplit];

subbanner["III. Packet A to Delta_branch"];

packetALanes = {
  {d200, d202, d204, n200, n202, n204},
  {d210, d212, d214, n210, n212, n214},
  {d220, d222, d224, n220, n222, n224}
};
denLanes = packetALanes[[All, {1, 2, 3}]];

responseNativeByLane = responsePair /@ denLanes;
prefactorNativeByLane = prefactorTriple /@ packetALanes;
responseTargetByLane = closedResponsePair /@ denLanes;
prefactorTargetByLane = closedPrefactorTriple /@ packetALanes;

nu2ByLane = responseNativeByLane[[All, 1]];
nu4ByLane = responseNativeByLane[[All, 2]];
p0ByLane = prefactorNativeByLane[[All, 1]];

nu2Coords = coordinatesFromBasis[nu2ByLane];
nu4Coords = coordinatesFromBasis[nu4ByLane];
p0Coords = coordinatesFromBasis[p0ByLane];

targetNu2Coords = sympyGroupCoords[responseTargetByLane[[All, 1]]];
targetNu4Coords = sympyGroupCoords[responseTargetByLane[[All, 2]]];
targetP0Coords = sympyGroupCoords[prefactorTargetByLane[[All, 1]]];

deltaBranchNative = cleanTensor[
  {
    nu2Coords[[2]],
    nu2Coords[[3]],
    nu4Coords[[2]],
    nu4Coords[[3]],
    p0Coords[[2]],
    p0Coords[[3]],
    nu4Coords[[1]] - 4*nu2Coords[[1]]^2,
    mhat^2*p0Coords[[1]] - 54*grav*cSound^5/(5*scaleA^5*cLight^5)
  }
];

deltaBranchTarget = cleanTensor[
  {
    targetNu2Coords[[2]],
    targetNu2Coords[[3]],
    targetNu4Coords[[2]],
    targetNu4Coords[[3]],
    targetP0Coords[[2]],
    targetP0Coords[[3]],
    targetNu4Coords[[1]] - 4*targetNu2Coords[[1]]^2,
    mhat^2*targetP0Coords[[1]] - 54*grav*cSound^5/(5*scaleA^5*cLight^5)
  }
];

isoRules = {
  d200 -> dIso0, d202 -> dIso2, d204 -> dIso4, n200 -> nIso0, n202 -> nIso2, n204 -> nIso4,
  d210 -> dIso0, d212 -> dIso2, d214 -> dIso4, n210 -> nIso0, n212 -> nIso2, n214 -> nIso4,
  d220 -> dIso0, d222 -> dIso2, d224 -> dIso4, n220 -> nIso0, n222 -> nIso2, n224 -> nIso4
};

isoResponse = responsePair[{dIso0, dIso2, dIso4}];
isoPrefactor = prefactorTriple[{dIso0, dIso2, dIso4, nIso0, nIso2, nIso4}];
p0Target = 54*grav*cSound^5/(5*scaleA^5*cLight^5*mhat^2);
onePoleNormRules = {dIso4 -> -3*dIso2^2/dIso0, nIso0 -> dIso0*p0Target};

Print["Delta_branch from native coefficient and basis compilers = ", fmt[deltaBranchNative]];
expectZero["native lane response coefficients - SymPy lane targets", responseNativeByLane - responseTargetByLane];
expectZero["native lane prefactor coefficients - SymPy lane targets", prefactorNativeByLane - prefactorTargetByLane];
expectZero["Delta_branch native route - SymPy packet", deltaBranchNative - deltaBranchTarget];
expectZero["anisotropy coordinates vanish on isotropic Packet A", deltaBranchNative[[1 ;; 6]] /. isoRules];
expectZero["mean nu2 on isotropic Packet A - one-lane nu2", (nu2Coords[[1]] /. isoRules) - isoResponse[[1]]];
expectZero["mean nu4 on isotropic Packet A - one-lane nu4", (nu4Coords[[1]] /. isoRules) - isoResponse[[2]]];
expectZero["mean P0 on isotropic Packet A - one-lane P0", (p0Coords[[1]] /. isoRules) - isoPrefactor[[1]]];
expectZero[
  "Delta_branch on isotropic one-pole normalized Packet A",
  deltaBranchNative /. isoRules /. onePoleNormRules
];

subbanner["IV. Packet B interconversion and Delta_orbit"];

qVars = {qTr, qNt, qEta};
mVars = {mT, mK, mMu};
rVars = {rTr, rNt, rEta};

rFromMTarget = {mT^(1 + chi0Star), mMu/(mK*mT^fStar), 1/mK};
qFromMNative = cleanLogTensor[Log /@ rFromMTarget];
qFromMTarget = {
  (1 + chi0Star)*Log[mT],
  Log[mMu] - Log[mK] - fStar*Log[mT],
  -Log[mK]
};

logToQ = {
  (1 + chi0Star)*ellT,
  ellMu - ellK - fStar*ellT,
  -ellK
};
logSolution = First[Solve[Thread[logToQ == qVars], {ellT, ellK, ellMu}, Reals]];
mFromQNative = cleanLogTensor[Exp /@ ({ellT, ellK, ellMu} /. logSolution)];
mFromQTarget = {
  Exp[qTr/(1 + chi0Star)],
  Exp[-qEta],
  Exp[qNt - qEta + fStar*qTr/(1 + chi0Star)]
};

rFromQNative = Exp /@ qVars;
qFromRNative = Log /@ rVars;
deltaOrbit = qVars;

Print["q from residual mismatches through invariant ratios = ", fmt[qFromMNative]];
Print["m from quotient coordinates by solved log equations = ", fmt[mFromQNative]];
expectZeroLog["q_from_m - SymPy quotient formulas", qFromMNative - qFromMTarget];
expectZeroLog["m_from_q - SymPy mismatch formulas", mFromQNative - mFromQTarget];
expectZeroLog["R_from_m after m_from_q - R_from_q", (rFromMTarget /. Thread[mVars -> mFromQNative]) - rFromQNative];
expectZeroLog["q_from_m after m_from_q - q", (qFromMTarget /. Thread[mVars -> mFromQNative]) - qVars];
expectZeroLog["q_from_R after R_from_q - q", (qFromRNative /. Thread[rVars -> rFromQNative]) - qVars];
expectZeroLog["m_from_q at orbit lock", (mFromQNative /. Thread[qVars -> {0, 0, 0}]) - {1, 1, 1}];
expectZeroLog["R_from_q at orbit lock", (rFromQNative /. Thread[qVars -> {0, 0, 0}]) - {1, 1, 1}];
expectZeroLog["q_from_m at mismatch lock", qFromMTarget /. Thread[mVars -> {1, 1, 1}]];
expectZero["Delta_orbit solve zero-set - quotient lock", qVars /. First[Solve[Thread[deltaOrbit == {0, 0, 0}], qVars, Reals]]];

subbanner["V. Two-packet home-stretch ledger"];

branchSlots = Array[br, 8];
orbitSlots = Array[or, 3];
formalClosureSplit = LogicalExpand[
  (And @@ Thread[Join[branchSlots, orbitSlots] == 0]) \[Equivalent]
    ((And @@ Thread[branchSlots == 0]) && (And @@ Thread[orbitSlots == 0]))
];
reducedClosurePacket = Join[deltaBranchNative, deltaOrbit];

expectZero["reduced closure packet branch head - Delta_branch", reducedClosurePacket[[1 ;; 8]] - deltaBranchNative];
expectZero["reduced closure packet orbit tail - Delta_orbit", reducedClosurePacket[[9 ;; 11]] - deltaOrbit];
expectTrue["formal closure zero-set iff both packet zero-sets vanish", formalClosureSplit];

banner["STAGE 191 MATHEMATICA LEDGER"];
Print["1. Taylor coefficients are recovered with SeriesCoefficient and match the SymPy compiler formulas."];
Print["2. Grouped coordinates are recovered by solving the weighted basis decomposition."];
Print["3. The weighted projectors are built from the diag(1,2,2) inner product and satisfy the projector algebra."];
Print["4. Packet A compiles to Delta_branch and vanishes on the isotropic one-pole normalized branch."];
Print["5. Packet B interconverts through solved log coordinates and Delta_orbit is the quotient-coordinate zero test."];
Print["6. The final reduced closure packet is exactly the branch packet followed by the orbit packet."];

Exit[0];
