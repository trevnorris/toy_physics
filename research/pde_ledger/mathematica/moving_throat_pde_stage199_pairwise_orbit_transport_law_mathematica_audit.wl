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
  Together[Expand[PowerExpand[stripCE[expr]]]],
  Assumptions -> $Assumptions
];

cleanTensor[expr_] := If[
  ListQ[expr],
  Map[cleanScalar, expr, {Length[Dimensions[expr]]}],
  cleanScalar[expr]
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

expectTrue[name_String, claim_] := Module[{res},
  res = FullSimplify[stripCE[claim], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 199 -- MATHEMATICA PAIRWISE ORBIT-TRANSPORT AUDIT"];

Clear[
  chi, delU, eStar, fStar,
  rLam, rC, rGam, rU, rK, rW, rMu, rT,
  rLam21, rC21, rGam21, rU21, rK21, rW21, rMu21, rT21,
  rLam32, rC32, rGam32, rU32, rK32, rW32, rMu32, rT32,
  baseT, baseK, baseMu,
  ellLam, ellC, ellGam, ellU, ellK, ellW, ellMu, ellT,
  tauLog, kappaLog, muLog
];

$Assumptions =
  Element[
    {
      chi, delU, eStar, fStar,
      rLam, rC, rGam, rU, rK, rW, rMu, rT,
      rLam21, rC21, rGam21, rU21, rK21, rW21, rMu21, rT21,
      rLam32, rC32, rGam32, rU32, rK32, rW32, rMu32, rT32,
      baseT, baseK, baseMu,
      ellLam, ellC, ellGam, ellU, ellK, ellW, ellMu, ellT,
      tauLog, kappaLog, muLog
    },
    Reals
  ] &&
  chi > 0 && delU > 0 && eStar > 0 && fStar > 0 &&
  rLam > 0 && rC > 0 && rGam > 0 && rU > 0 &&
  rK > 0 && rW > 0 && rMu > 0 && rT > 0 &&
  rLam21 > 0 && rC21 > 0 && rGam21 > 0 && rU21 > 0 &&
  rK21 > 0 && rW21 > 0 && rMu21 > 0 && rT21 > 0 &&
  rLam32 > 0 && rC32 > 0 && rGam32 > 0 && rU32 > 0 &&
  rK32 > 0 && rW32 > 0 && rMu32 > 0 && rT32 > 0 &&
  baseT > 0 && baseK > 0 && baseMu > 0;

ratioVars = {rLam, rC, rGam, rU, rK, rW, rMu, rT};
logVars = {ellLam, ellC, ellGam, ellU, ellK, ellW, ellMu, ellT};
deltaLogs = Log /@ ratioVars;
logRules = Thread[logVars -> deltaLogs];

(* Ordered microscopic basis: lambda, c, gamma, U, K_eta, W, mu, T. *)
chiCoreWeights = {0, 1, 1, -1, 0, 0, 0, 0};
thermalWeights = {0, 0, 0, -1, 0, 0, 0, 1};
nontrackPrefactorWeights = {2, 0, 0, 0, -1, -2, 1, 0};
epsWWeights = {2, 0, 2, -1, 0, -1, 0, 0};
epsEtaWeights = {0, 2, 0, -1, -1, 0, 0, 0};

compilerRows = {
  (1 + delU)*chiCoreWeights + (1 + chi)*thermalWeights,
  nontrackPrefactorWeights + eStar*epsWWeights - fStar*thermalWeights,
  epsEtaWeights
};
compilerTarget = {
  {0, 1 + delU, 1 + delU, -(2 + chi + delU), 0, 0, 0, 1 + chi},
  {2*(1 + eStar), 0, 2*eStar, fStar - eStar, -1, -(2 + eStar), 1, -fStar},
  {0, 2, 0, -1, -1, 0, 0, 0}
};

logForm[expr_] := cleanScalar[Log[expr]];
rowLogRatios = cleanTensor[compilerRows . deltaLogs];

subbanner["I. Pairwise monomial compiler"];

expectZero["derived M_* rows - SymPy compiler rows", compilerRows - compilerTarget];

ctrClosed = (rGam*rC/rU)^(1 + delU)*(rT/rU)^(1 + chi);
cntClosed =
  (rLam^2*rMu/(rK*rW^2))*
  (rGam^2*rLam^2/(rU*rW))^eStar*
  (rT/rU)^(-fStar);
epsClosed = rC^2/(rU*rK);

Print["log C_tr ratio from exponent ledger = ", fmt[rowLogRatios[[1]]]];
Print["log C_nt ratio from exponent ledger = ", fmt[rowLogRatios[[2]]]];
Print["log eps_eta ratio from exponent ledger = ", fmt[rowLogRatios[[3]]]];

expectZero["ln C_tr ratio - ln SymPy closed form", rowLogRatios[[1]] - logForm[ctrClosed]];
expectZero["ln C_nt ratio - ln SymPy closed form", rowLogRatios[[2]] - logForm[cntClosed]];
expectZero["ln eps_eta ratio - ln SymPy closed form", rowLogRatios[[3]] - logForm[epsClosed]];

subbanner["II. Transport factors from the native same-orbit solve"];

depSlots = {8, 5, 7};      (* T, K_eta, mu *)
freeSlots = {1, 2, 3, 4, 6};
depLogVars = logVars[[depSlots]];
freeLogVars = logVars[[freeSlots]];

sameOrbitSolve = First[
  Solve[Thread[compilerRows . logVars == ConstantArray[0, 3]], depLogVars, Reals]
];

phiLogTemplate = cleanTensor[depLogVars /. sameOrbitSolve];
phiLogFor[ratios_List] := cleanTensor[
  phiLogTemplate /. Thread[freeLogVars -> (Log /@ ratios[[freeSlots]])]
];

{phiTLog, phiKLog, phiMuLog} = phiLogFor[ratioVars];
alpha = (1 + delU)/(1 + chi);

phiTClosed = rU*(rU/(rGam*rC))^alpha;
phiKClosed = rC^2/rU;
phiMuClosed =
  phiKClosed*rW^2/rLam^2*
  (rGam^2*rLam^2/(rU*rW))^(-eStar)*
  (phiTClosed/rU)^fStar;
phiMuExpanded =
  rC^(2 - fStar*alpha)*
  rGam^(-2*eStar - fStar*alpha)*
  rLam^(-2*(1 + eStar))*
  rU^(-1 + eStar + fStar*alpha)*
  rW^(2 + eStar);

Print["ln Phi_T from Solve = ", fmt[phiTLog]];
Print["ln Phi_K from Solve = ", fmt[phiKLog]];
Print["ln Phi_mu from Solve = ", fmt[phiMuLog]];

expectZero["ln Phi_T solve - ln SymPy transport", phiTLog - logForm[phiTClosed]];
expectZero["ln Phi_K solve - ln SymPy transport", phiKLog - logForm[phiKClosed]];
expectZero["ln Phi_mu solve - ln SymPy factorized transport", phiMuLog - logForm[phiMuClosed]];
expectZero["ln Phi_mu solve - ln expanded monomial form", phiMuLog - logForm[phiMuExpanded]];

sameOrbitLogs = cleanTensor[deltaLogs /. {
  Log[rT] -> phiTLog,
  Log[rK] -> phiKLog,
  Log[rMu] -> phiMuLog
}];
expectZero["same-orbit compiler packet", compilerRows . sameOrbitLogs];

subbanner["III. Reference-independent mismatch collapse"];

mTLog = cleanScalar[Log[rT] - phiTLog];
mKLog = cleanScalar[Log[rK] - phiKLog];
mMuLog = cleanScalar[Log[rMu] - phiMuLog];

Print["ln m_T = ", fmt[mTLog]];
Print["ln m_K = ", fmt[mKLog]];
Print["ln m_mu = ", fmt[mMuLog]];

expectZero["ln C_tr ratio - (1+chi) ln m_T", rowLogRatios[[1]] - (1 + chi)*mTLog];
expectZero["ln eps_eta ratio + ln m_K", rowLogRatios[[3]] + mKLog];
expectZero[
  "ln C_nt ratio - [ln m_mu - ln m_K - F ln m_T]",
  rowLogRatios[[2]] - (mMuLog - mKLog - fStar*mTLog)
];

subbanner["IV. q-chart and finite projector split"];

qPair = rowLogRatios;
qTr = qPair[[1]];
qNt = qPair[[2]];
qEta = qPair[[3]];

expectZero["q_tr - (1+chi) ln m_T", qTr - (1 + chi)*mTLog];
expectZero["q_eta + ln m_K", qEta + mKLog];
expectZero["q_nt - [ln m_mu - ln m_K - F ln m_T]", qNt - (mMuLog - mKLog - fStar*mTLog)];

freeSelector = UnitVector[8, #] & /@ freeSlots;
constrainedSystem = Join[compilerRows, freeSelector];
sectionFromSolve = cleanTensor[
  Transpose[
    Table[
      LinearSolve[
        constrainedSystem,
        Join[UnitVector[3, col], ConstantArray[0, Length[freeSlots]]]
      ],
      {col, 3}
    ]
  ]
];
quotientProjector = cleanTensor[sectionFromSolve . compilerRows];
orbitProjector = cleanTensor[IdentityMatrix[8] - quotientProjector];

qPart = cleanTensor[quotientProjector . deltaLogs];
oPart = cleanTensor[orbitProjector . deltaLogs];
expectedQPart = {0, 0, 0, 0, mKLog, 0, mMuLog, mTLog};
expectedOPart = {Log[rLam], Log[rC], Log[rGam], Log[rU], phiKLog, Log[rW], phiMuLog, phiTLog};

Print["Q_quot Delta x = ", fmt[qPart]];
Print["O_orb Delta x = ", fmt[oPart]];

expectZero["Q_quot Delta x - SymPy mismatch support", qPart - expectedQPart];
expectZero["O_orb Delta x - SymPy transported orbit part", oPart - expectedOPart];
expectZero["Q + O - Delta x", qPart + oPart - deltaLogs];
expectZero["M_* O_orb Delta x", compilerRows . oPart];
expectZero["M_* Q_quot Delta x - q", compilerRows . qPart - qPair];
expectZero["Q^2 - Q", quotientProjector . quotientProjector - quotientProjector];
expectZero["O^2 - O", orbitProjector . orbitProjector - orbitProjector];
expectZero["Q O", quotientProjector . orbitProjector];
expectZero["O Q", orbitProjector . quotientProjector];

subbanner["V. Pairwise restoration map"];

restoredTLog = cleanScalar[Log[baseT] + Log[rT] - qTr/(1 + chi)];
restoredKLog = cleanScalar[Log[baseK] + Log[rK] + qEta];
restoredMuLog = cleanScalar[Log[baseMu] + Log[rMu] - qNt + qEta - fStar*qTr/(1 + chi)];

expectZero["ln T_restore - ln(Phi_T T_1)", restoredTLog - (Log[baseT] + phiTLog)];
expectZero["ln K_restore - ln(Phi_K K_1)", restoredKLog - (Log[baseK] + phiKLog)];
expectZero["ln mu_restore - ln(Phi_mu mu_1)", restoredMuLog - (Log[baseMu] + phiMuLog)];

subbanner["VI. Cocycle and composition laws"];

r21 = {rLam21, rC21, rGam21, rU21, rK21, rW21, rMu21, rT21};
r32 = {rLam32, rC32, rGam32, rU32, rK32, rW32, rMu32, rT32};
r31 = MapThread[Times, {r32, r21}];

deltaFrom[ratios_List] := Log /@ ratios;
qFrom[ratios_List] := cleanTensor[compilerRows . deltaFrom[ratios]];
mLogsFrom[ratios_List] := Module[{phiLogs = phiLogFor[ratios]},
  cleanTensor[
    {
      Log[ratios[[8]]] - phiLogs[[1]],
      Log[ratios[[5]]] - phiLogs[[2]],
      Log[ratios[[7]]] - phiLogs[[3]]
    }
  ]
];

phi21 = phiLogFor[r21];
phi32 = phiLogFor[r32];
phi31 = phiLogFor[r31];
m21 = mLogsFrom[r21];
m32 = mLogsFrom[r32];
m31 = mLogsFrom[r31];

expectZero["ln Phi^31 - ln Phi^32 - ln Phi^21", phi31 - phi32 - phi21];
expectZero["ln m^31 - ln m^32 - ln m^21", m31 - m32 - m21];
expectZero["Delta x^31 - Delta x^32 - Delta x^21", deltaFrom[r31] - deltaFrom[r32] - deltaFrom[r21]];
expectZero["q^31 - q^32 - q^21", qFrom[r31] - qFrom[r32] - qFrom[r21]];

subbanner["VII. Two-point orbit lock"];

zeroQLogSolve = First[
  Solve[Thread[compilerRows . logVars == ConstantArray[0, 3]], depLogVars, Reals]
];
zeroProjectorSolve = First[
  Solve[Thread[quotientProjector . logVars == ConstantArray[0, 8]], depLogVars, Reals]
];

expectZero[
  "Q Delta x == 0 and M Delta x == 0 solve laws agree",
  (depLogVars /. zeroProjectorSolve) - (depLogVars /. zeroQLogSolve)
];
expectZero["Q Delta x vanishes on same-orbit law", quotientProjector . (logVars /. zeroQLogSolve)];
expectZero["M Delta x vanishes on Q-kernel law", compilerRows . (logVars /. zeroProjectorSolve)];

mismatchToQ = {
  {1 + chi, 0, 0},
  {-fStar, -1, 1},
  {0, -1, 0}
};
Print["det(mismatch-to-q chart) = ", fmt[Factor[Det[mismatchToQ]]]];
expectZero["det mismatch-to-q chart - (1+chi)", Factor[Det[mismatchToQ]] - (1 + chi)];
lockMismatchSolve = First[
  Solve[
    Thread[mismatchToQ . {tauLog, kappaLog, muLog} == ConstantArray[0, 3]],
    {tauLog, kappaLog, muLog},
    Reals
  ]
];
expectZero["q == 0 forces zero logarithmic mismatch", {tauLog, kappaLog, muLog} /. lockMismatchSolve];

subbanner["VIII. Fixed-free-coordinate reduction"];

freeEqualRules = {rLam -> 1, rC -> 1, rGam -> 1, rU -> 1, rW -> 1};

expectZero["ln Phi_T at equal free coordinates", phiTLog /. freeEqualRules];
expectZero["ln Phi_K at equal free coordinates", phiKLog /. freeEqualRules];
expectZero["ln Phi_mu at equal free coordinates", phiMuLog /. freeEqualRules];
expectZero["ln m_T at equal free coordinates - ln r_T", (mTLog /. freeEqualRules) - Log[rT]];
expectZero["ln m_K at equal free coordinates - ln r_K", (mKLog /. freeEqualRules) - Log[rK]];
expectZero["ln m_mu at equal free coordinates - ln r_mu", (mMuLog /. freeEqualRules) - Log[rMu]];

banner["STAGE 199 MATHEMATICA LEDGER"];
Print["1. Primitive monomial exponent vectors reproduce the SymPy pairwise ratios."];
Print["2. A native same-orbit Solve derives Phi_T, Phi_K, and Phi_mu."];
Print["3. Dividing raw dependent ratios by those transport laws gives the mismatch triple."];
Print["4. LinearSolve constructs the Stage 192 quotient section and finite pairwise projectors."];
Print["5. The restoration, cocycle, orbit-lock, and fixed-base reduction laws all match the SymPy conclusions."];

Exit[0];
