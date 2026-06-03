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

cleanResidual[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

zeroResidualQ[res_] := If[
  ListQ[res],
  TrueQ[And @@ (TrueQ[# === 0] & /@ Flatten[res])],
  TrueQ[res === 0]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanResidual[expr];
  Print[name, " = ", fmt[res]];
  If[zeroResidualQ[res], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  res = stripConditional[res];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

$Assumptions = (
  Element[
    {epsEtaStar, chi0Star, FStar, qNt, qEta, qTr, R1, E1,
     dT, dK, dMu, uNt, uEta},
    Reals
  ] &&
  0 < epsEtaStar < 1
);

banner[
  "Stage 236 Mathematica audit: rigid-mouth dependent-plane projectors and restoration gap"
];

cEta = epsEtaStar/(1 - epsEtaStar);
obsVec = {R1, E1};
packetCoords = {qNt, qEta};
deltaCoords = {dT, dK, dMu};
deltaVec = deltaCoords;

rowFromLinearAction[expr_] := (
  Coefficient[FullSimplify[expr, Assumptions -> $Assumptions], #] & /@
    deltaCoords
);
matrixFromLinearAction[action_List] := rowFromLinearAction /@ action;

subbanner["M1. Rigid-mouth packet map"];

mRm = {{-1, -cEta}, {0, 1}};
qRmFromObs = FullSimplify[mRm . obsVec, Assumptions -> $Assumptions];

expectZero[
  "M1 packet action",
  qRmFromObs - {-R1 - cEta E1, E1}
];

Print["M_rm . {R1, E1} = ", fmt[qRmFromObs]];

subbanner["M2. Dependent section on q_tr = 0"];

generalSection = {
  qTr/(1 + chi0Star),
  -qEta,
  FStar qTr/(1 + chi0Star) + qNt - qEta
};
yRm = FullSimplify[generalSection /. qTr -> 0, Assumptions -> $Assumptions];
sRmDep = {{0, 0}, {0, -1}, {1, -1}};

expectZero[
  "M2 rigid-mouth section",
  yRm - {0, -qEta, qNt - qEta}
];
expectZero[
  "M2 section matrix action",
  yRm - sRmDep . packetCoords
];

Print["y_rm = ", fmt[yRm]];

subbanner["M3. Direct-observable-to-dependent compiler by basis actions"];

compilerColumns = FullSimplify[
  sRmDep . (mRm . #) & /@ IdentityMatrix[2],
  Assumptions -> $Assumptions
];
cRmDep = Transpose[compilerColumns];
cExpected = {{0, 0}, {0, -1}, {-1, -1/(1 - epsEtaStar)}};
yFromObs = FullSimplify[cRmDep . obsVec, Assumptions -> $Assumptions];

expectZero[
  "M3 basis-action compiler entries",
  cRmDep - cExpected
];
expectZero[
  "M3 product and basis-action compiler agree on generic input",
  sRmDep . (mRm . obsVec) - cRmDep . obsVec
];
expectZero[
  "M3 compiler action",
  yFromObs - {0, -E1, -R1 - E1/(1 - epsEtaStar)}
];

Print["C_rm_dep = ", fmt[cRmDep]];
Print["C_rm_dep . {R1, E1} = ", fmt[yFromObs]];

subbanner["M4. Recovery map from the Delta_T = 0 plane"];

planeBlock = sRmDep[[{2, 3}, All]];
recoveredByLinearSolve = FullSimplify[
  LinearSolve[planeBlock, {dK, dMu}],
  Assumptions -> $Assumptions
];
recoveredBySolve = FullSimplify[
  {uNt, uEta} /. First[
    Solve[sRmDep . {uNt, uEta} == {0, dK, dMu}, {uNt, uEta}]
  ],
  Assumptions -> $Assumptions
];
lRecovered = matrixFromLinearAction[recoveredByLinearSolve];
lExpected = {{0, -1, 1}, {0, -1, 0}};

expectZero[
  "M4 LinearSolve and Solve recovery agree",
  recoveredByLinearSolve - recoveredBySolve
];
expectZero[
  "M4 recovered left-inverse matrix",
  lRecovered - lExpected
];
expectZero[
  "M4 L_rm_dep . S_rm_dep",
  lRecovered . sRmDep - IdentityMatrix[2]
];
expectZero[
  "M4 recovery on y_rm",
  lRecovered . yRm - packetCoords
];

Print["Recovered {q_nt, q_eta} from {0, Delta_Keta, Delta_mu} = ",
  fmt[recoveredByLinearSolve]];
Print["L_rm_dep = ", fmt[lRecovered]];

subbanner["M5. Dependent-plane projectors from recovered coordinates"];

selectNt = DiagonalMatrix[{1, 0}];
selectEta = DiagonalMatrix[{0, 1}];
recoveredGeneric = lRecovered . deltaVec;
pNtAction = FullSimplify[
  sRmDep . (selectNt . recoveredGeneric),
  Assumptions -> $Assumptions
];
pEtaAction = FullSimplify[
  sRmDep . (selectEta . recoveredGeneric),
  Assumptions -> $Assumptions
];
pNtDep = matrixFromLinearAction[pNtAction];
pEtaDep = matrixFromLinearAction[pEtaAction];
pNtExpected = {{0, 0, 0}, {0, 0, 0}, {0, -1, 1}};
pEtaExpected = {{0, 0, 0}, {0, 1, 0}, {0, 1, 0}};

expectZero[
  "M5 P_nt entries from recovered action",
  pNtDep - pNtExpected
];
expectZero[
  "M5 P_eta entries from recovered action",
  pEtaDep - pEtaExpected
];
expectZero[
  "M5 P_nt idempotent on generic vector",
  pNtDep . (pNtDep . deltaVec) - pNtDep . deltaVec
];
expectZero[
  "M5 P_eta idempotent on generic vector",
  pEtaDep . (pEtaDep . deltaVec) - pEtaDep . deltaVec
];
expectZero[
  "M5 P_nt after P_eta annihilates generic vector",
  pNtDep . (pEtaDep . deltaVec)
];
expectZero[
  "M5 P_eta after P_nt annihilates generic vector",
  pEtaDep . (pNtDep . deltaVec)
];
expectZero[
  "M5 P_nt eigenvalues",
  Sort[Eigenvalues[pNtDep]] - {0, 0, 1}
];
expectZero[
  "M5 P_eta eigenvalues",
  Sort[Eigenvalues[pEtaDep]] - {0, 0, 1}
];
expectZero["M5 P_nt rank", MatrixRank[pNtDep] - 1];
expectZero["M5 P_eta rank", MatrixRank[pEtaDep] - 1];

Print["P_nt_dep = ", fmt[pNtDep]];
Print["P_eta_dep = ", fmt[pEtaDep]];

subbanner["M6. Plane completeness"];

planeIdentity = DiagonalMatrix[{0, 1, 1}];

expectZero[
  "M6 projector sum",
  pNtDep + pEtaDep - planeIdentity
];
expectZero[
  "M6 projector sum acts as identity on Delta_T = 0",
  (pNtDep + pEtaDep) . {0, dK, dMu} - {0, dK, dMu}
];

subbanner["M7. Packet decomposition"];

yNt = FullSimplify[pNtDep . yRm, Assumptions -> $Assumptions];
yEta = FullSimplify[pEtaDep . yRm, Assumptions -> $Assumptions];

expectZero["M7 y_nt", yNt - {0, 0, qNt}];
expectZero["M7 y_eta", yEta + qEta {0, 1, 1}];
expectZero["M7 y_rm decomposition", yRm - yNt - yEta];

Print["y_nt = ", fmt[yNt]];
Print["y_eta = ", fmt[yEta]];

subbanner["M8. Static strip and equal-drift ray"];

deltaKetaRm = yRm[[2]];
deltaMuRm = yRm[[3]];

expectZero[
  "M8 Delta_mu - Delta_Keta residual",
  (deltaMuRm - deltaKetaRm) - qNt
];
expectTrue[
  "M8 Delta_mu equals Delta_Keta iff q_nt is zero",
  Equivalent[deltaMuRm == deltaKetaRm, qNt == 0]
];

stripRuleByE = {R1 -> -cEta E1};
stripRuleByR = {E1 -> -R1/cEta};
stripVectorByE = FullSimplify[yFromObs /. stripRuleByE, Assumptions -> $Assumptions];
stripVectorByR = FullSimplify[yFromObs /. stripRuleByR, Assumptions -> $Assumptions];

expectZero[
  "M8 equal-drift strip vector by E1",
  stripVectorByE - {0, -E1, -E1}
];
expectZero[
  "M8 equal-drift ray by E1",
  stripVectorByE - ((R1/cEta) {0, 1, 1} /. stripRuleByE)
];
expectZero[
  "M8 equal-drift ray by R1",
  stripVectorByR - (R1/cEta) {0, 1, 1}
];

normEta = FullSimplify[yEta . yEta, Assumptions -> $Assumptions];
stripNormByE = FullSimplify[stripVectorByE . stripVectorByE, Assumptions -> $Assumptions];
stripNormByR = FullSimplify[stripVectorByR . stripVectorByR, Assumptions -> $Assumptions];

expectZero["M8 eta norm", normEta - 2 qEta^2];
expectZero["M8 strip norm by E1", stripNormByE - 2 E1^2];
expectZero[
  "M8 strip norm by R1",
  stripNormByR - 2 R1^2/cEta^2
];

Print["Equal-drift strip vector = ", fmt[stripVectorByE]];
Print["||y_eta||^2 = ", fmt[normEta]];

subbanner["M9. Microscopic correction compilers"];

deltaYStatic = FullSimplify[-yNt, Assumptions -> $Assumptions];
deltaYOrbit = FullSimplify[-yRm, Assumptions -> $Assumptions];

expectZero[
  "M9 static correction",
  deltaYStatic - {0, 0, -qNt}
];
expectZero[
  "M9 static correction leaves eta packet",
  yRm + deltaYStatic - yEta
];
expectZero[
  "M9 orbit correction decomposition",
  deltaYOrbit - (deltaYStatic + qEta {0, 1, 1})
];
expectZero[
  "M9 orbit correction closes y_rm",
  yRm + deltaYOrbit
];

Print["Delta_y_static = ", fmt[deltaYStatic]];
Print["Delta_y_orbit = ", fmt[deltaYOrbit]];

Print[""];
Print["All Stage 236 Mathematica checks passed."];
