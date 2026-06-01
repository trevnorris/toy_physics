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

samePositive[name_String, lhs_, rhs_] := expectZeroLog[name, Log[lhs/rhs]];

banner["STAGE 198 -- EXACT FINITE ORBIT LAW MATHEMATICA AUDIT"];

Clear[
  lambdaW, cEtaU, gammaC, kU, kEta, kW, muW, tU, ell, sigmaP,
  chiS, deltaS, eS, fS, cTrS, cNtS, epsS,
  zK, zT, zM, sT, sK, sM, tau, kappa, rho,
  dLam, dC, dG, dU, dK, dW, dMu, dT, qTr, qNt, qEta
];

$Assumptions = (
  Element[
    {
      lambdaW, cEtaU, gammaC, kU, kEta, kW, muW, tU, ell, sigmaP,
      chiS, deltaS, eS, fS, cTrS, cNtS, epsS, sT, sK, sM
    },
    Reals
  ] &&
  lambdaW > 0 && cEtaU > 0 && gammaC > 0 && kU > 0 &&
  kEta > 0 && kW > 0 && muW > 0 && tU > 0 &&
  ell > 0 && sigmaP > 0 && chiS > 0 && deltaS > 0 &&
  eS > 0 && fS > 0 && cTrS > 0 && cNtS > 0 &&
  epsS > 0 && sT > 0 && sK > 0 && sM > 0 &&
  Element[
    {zK, zT, zM, tau, kappa, rho, dLam, dC, dG, dU, dK, dW, dMu, dT,
     qTr, qNt, qEta},
    Reals
  ]
);

cTrLaw = (
  (gammaC*cEtaU/kU)^(1 + deltaS) *
  (Pi^2*tU/(ell^2*kU))^(1 + chiS)
);

cNtLaw = (
  (lambdaW^2*muW/(kEta*kW^2)) *
  (gammaC^2*lambdaW^2*sigmaP/(kU*kW))^eS *
  (Pi^2*tU/(ell^2*kU))^(-fS)
);

epsEtaLaw = cEtaU^2/(kU*kEta);

subbanner["I. Orbit solve from log-linear invariant equations"];

depLogVars = {zK, zT, zM};
depLogRules = {kEta -> Exp[zK], tU -> Exp[zT], muW -> Exp[zM]};
invariantLogResiduals = cleanLogTensor[
  Log[({epsEtaLaw, cTrLaw, cNtLaw}/{epsS, cTrS, cNtS}) /. depLogRules]
];
depCoefficientMatrix = Table[
  Coefficient[invariantLogResiduals[[row]], depLogVars[[col]]],
  {row, 3}, {col, 3}
];
depOffset = cleanLogTensor[
  invariantLogResiduals /. Thread[depLogVars -> ConstantArray[0, 3]]
];
logOrbitSolution = cleanLogTensor[
  LinearSolve[depCoefficientMatrix, -depOffset]
];
orbitFromLogs = PowerExpand[Exp /@ logOrbitSolution];

Print["log invariant residuals = ", fmt[invariantLogResiduals]];
Print["dependent log coefficient matrix = ", fmt[depCoefficientMatrix]];
Print["log orbit solution = ", fmt[logOrbitSolution]];

expectZero["det(dependent log matrix) + (1+chiS)", Det[depCoefficientMatrix] + (1 + chiS)];

kEtaSympy = cEtaU^2/(kU*epsS);
tUSympy = (
  (ell^2*kU/Pi^2) *
  (cTrS/(gammaC*cEtaU/kU)^(1 + deltaS))^(1/(1 + chiS))
);
muWSympy = (
  cNtS*kEtaSympy*kW^2/lambdaW^2 *
  (gammaC^2*lambdaW^2*sigmaP/(kU*kW))^(-eS) *
  (Pi^2*tUSympy/(ell^2*kU))^fS
);

orbitRules = Thread[{kEta, tU, muW} -> orbitFromLogs];

samePositive["K_eta orbit agrees with SymPy target", orbitFromLogs[[1]], kEtaSympy];
samePositive["T_U orbit agrees with SymPy target", orbitFromLogs[[2]], tUSympy];
samePositive["mu_W orbit agrees with SymPy target", orbitFromLogs[[3]], muWSympy];
samePositive["epsilon_eta at log-solved orbit", epsEtaLaw /. orbitRules, epsS];
samePositive["C_tr at log-solved orbit", cTrLaw /. orbitRules, cTrS];
samePositive["C_nt at log-solved orbit", cNtLaw /. orbitRules, cNtS];

subbanner["II. Multiplicative mismatch ratios"];

scaleRules = {tU -> sT*tUSympy, kEta -> sK*kEtaSympy, muW -> sM*muWSympy};

samePositive[
  "C_tr(actual)/C_tr_star agrees with sT^(1+chiS)",
  (cTrLaw /. scaleRules)/cTrS,
  sT^(1 + chiS)
];
samePositive[
  "epsilon_eta(actual)/epsilon_eta_star agrees with 1/sK",
  (epsEtaLaw /. scaleRules)/epsS,
  1/sK
];
samePositive[
  "C_nt(actual)/C_nt_star agrees with sM/(sK sT^fS)",
  (cNtLaw /. scaleRules)/cNtS,
  sM/(sK*sT^fS)
];

subbanner["III. Log chart and Stage 192 drift rows from finite ratios"];

logActualRules = {
  tU -> Exp[tau]*tUSympy,
  kEta -> Exp[kappa]*kEtaSympy,
  muW -> Exp[rho]*muWSympy
};
packetFromRatios = cleanLogTensor[
  Log[({cTrLaw/cTrS, cNtLaw/cNtS, epsEtaLaw/epsS} /. logActualRules)]
];
packetExpected = {(1 + chiS)*tau, rho - kappa - fS*tau, -kappa};

expectZeroLog["Packet-B log chart agrees with SymPy q chart", packetFromRatios - packetExpected];

driftVars = {dLam, dC, dG, dU, dK, dW, dMu, dT};
driftRules = {
  lambdaW -> Exp[dLam]*lambdaW,
  cEtaU -> Exp[dC]*cEtaU,
  gammaC -> Exp[dG]*gammaC,
  kU -> Exp[dU]*kU,
  kEta -> Exp[dK]*kEta,
  kW -> Exp[dW]*kW,
  muW -> Exp[dMu]*muW,
  tU -> Exp[dT]*tU
};
finiteLogPacket = cleanLogTensor[
  Log[(({cTrLaw, cNtLaw, epsEtaLaw} /. driftRules)/{cTrLaw, cNtLaw, epsEtaLaw})]
];
compilerByJacobian = Table[
  cleanLogScalar[D[finiteLogPacket[[row]], driftVars[[col]]]],
  {row, 3}, {col, 8}
];
compilerExpected = {
  {0, 1 + deltaS, 1 + deltaS, -(2 + chiS + deltaS), 0, 0, 0, 1 + chiS},
  {2*(1 + eS), 0, 2*eS, fS - eS, -1, -(2 + eS), 1, -fS},
  {0, 2, 0, -1, -1, 0, 0, 0}
};
pureDependentMismatch = {0, 0, 0, 0, kappa, 0, rho, tau};

Print["finite log packet = ", fmt[finiteLogPacket]];
Print["Jacobian compiler = ", fmt[compilerByJacobian]];
expectZero["finite-ratio Jacobian agrees with SymPy M_*", compilerByJacobian - compilerExpected];
expectZero[
  "M_* pure dependent mismatch agrees with q chart",
  compilerByJacobian . pureDependentMismatch - packetExpected
];

subbanner["IV. Restoration as a dependent-column linear solve"];

dependentColumns = compilerByJacobian[[All, {5, 7, 8}]];
restorationSolve = cleanTensor[LinearSolve[dependentColumns, -packetExpected]];
restorationExpected = {
  packetExpected[[3]],
  -packetExpected[[2]] + packetExpected[[3]] - fS*packetExpected[[1]]/(1 + chiS),
  -packetExpected[[1]]/(1 + chiS)
};

Print["dependent correction exponents {K_eta, mu_W, T_U} = ", fmt[restorationSolve]];
expectZero[
  "restoration exponents agree with SymPy restoration map",
  restorationSolve - restorationExpected
];
expectZero[
  "restoration kills the Packet-B coordinates",
  dependentColumns . restorationSolve + packetExpected
];

samePositive[
  "T_U restoration returns orbit",
  Exp[tau + restorationSolve[[3]]]*tUSympy,
  tUSympy
];
samePositive[
  "K_eta restoration returns orbit",
  Exp[kappa + restorationSolve[[1]]]*kEtaSympy,
  kEtaSympy
];
samePositive[
  "mu_W restoration returns orbit",
  Exp[rho + restorationSolve[[2]]]*muWSympy,
  muWSympy
];

subbanner["V. Inverse chart and finite orbit-lock criterion"];

inverseChart = First[Solve[
  Thread[packetExpected == {qTr, qNt, qEta}],
  {tau, kappa, rho},
  Reals
]];
inverseExpected = {qTr/(1 + chiS), -qEta, qNt - qEta + fS*qTr/(1 + chiS)};

Print["inverse chart = ", fmt[inverseChart]];
expectZero[
  "inverse chart agrees with SymPy target",
  ({tau, kappa, rho} /. inverseChart) - inverseExpected
];
expectZero[
  "direct chart after inverse returns q",
  (packetExpected /. Thread[{tau, kappa, rho} -> ({tau, kappa, rho} /. inverseChart)]) -
    {qTr, qNt, qEta}
];

mismatchFromQ = Exp /@ ({tau, kappa, rho} /. inverseChart);
expectZero["m_T at q=0 - 1", (mismatchFromQ[[1]] /. qTr -> 0) - 1];
expectZero["m_K at q=0 - 1", (mismatchFromQ[[2]] /. qEta -> 0) - 1];
expectZero[
  "m_mu at q=0 - 1",
  (mismatchFromQ[[3]] /. {qTr -> 0, qNt -> 0, qEta -> 0}) - 1
];

zeroPacket = And @@ Thread[packetExpected == ConstantArray[0, 3]];
zeroLogs = tau == 0 && kappa == 0 && rho == 0;
expectTrue[
  "finite orbit-lock equivalence in log coordinates",
  Equivalent[zeroPacket, zeroLogs]
];

banner["STAGE 198 MATHEMATICA LEDGER"];
Print["1. The dependent triple is solved by a log-linear invariant system with determinant -(1+chiS)."];
Print["2. The solved orbit agrees with the SymPy formulas for K_eta, T_U, and mu_W."];
Print["3. Positive mismatch scalings give the exact C_tr, epsilon_eta, and C_nt ratios."];
Print["4. The Packet-B chart and the Stage 192 drift compiler are recovered from finite log-ratio Jacobians."];
Print["5. A dependent-column linear solve gives the restoration exponents and returns the orbit point."];
Print["6. Solve confirms the inverse chart and the finite orbit-lock condition q=0 iff all log mismatches vanish."];

Exit[0];
