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

cleanScalar[expr_, ass_: True] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[res], Assumptions -> ass];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> ass]
];

expectZero[name_String, expr_, ass_: True] := Module[{res},
  res = cleanScalar[expr, ass];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_, ass_: True] := Module[{res},
  res = FullSimplify[Refine[cond, ass], Assumptions -> ass];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

expectFalse[name_String, cond_, ass_: True] := Module[{res},
  res = FullSimplify[Refine[cond, ass], Assumptions -> ass];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === False], pass[name], fail[name, res]];
];

banner["STAGE 206 -- CERTIFIED RAY RANKING AND LOCAL BRACKETING"];

Clear[
  H0, k, c, cL, cU, c0, cmid, eta, a, Klower, tau,
  tauLoA, tauHiA, tauStarA, tauLoB, tauHiB, tauStarB
];

baseAssumptions = (
  Element[{H0, k, c, cL, cU, c0, cmid, eta, a, Klower}, Reals] &&
  H0 > 0 && k > 0
);
K0 = -k;
disc[cc_] := K0^2 - 2 cc H0;
branchAssumptions[cc_] := baseAssumptions && disc[cc] > 0;
q[cc_, x_] := H0 + K0 x + 1/2 cc x^2;

subbanner["M1. Solve-selected oriented root map"];

solveRoots = stripConditional[tau /. Solve[q[c, tau] == 0, tau]];
linearLimit[root_] := Assuming[H0 > 0 && k > 0, Limit[root, c -> 0]];
selectedRoot = SelectFirst[
  solveRoots,
  TrueQ[FullSimplify[linearLimit[#] == H0/k, Assumptions -> H0 > 0 && k > 0]] &
];
If[MissingQ[selectedRoot], fail["M1 physical branch selection", solveRoots]];

rootMap = -2 H0/(K0 + Sign[K0] Sqrt[disc[c]]);
Print["M1 solve roots = ", fmt[solveRoots]];
Print["M1 selected root = ", fmt[selectedRoot]];
Print["M1 oriented root map = ", fmt[rootMap]];
expectZero["M1 closed root solves quadratic", q[c, rootMap], branchAssumptions[c]];
expectZero["M1 closed root - Solve-selected branch", rootMap - selectedRoot, branchAssumptions[c]];

subbanner["M2. Zero-curvature limit"];

expectZero[
  "M2 zero-curvature limit",
  Limit[rootMap, c -> 0] - H0/k,
  H0 > 0 && k > 0
];

subbanner["M3. Curvature monotonicity"];

dRootDc = D[rootMap, c];
expectZero[
  "M3 curvature derivative identity",
  dRootDc - rootMap^2/(2 Sqrt[disc[c]]),
  branchAssumptions[c]
];
expectTrue[
  "M3 curvature derivative sign",
  dRootDc > 0,
  branchAssumptions[c]
];

subbanner["M4. Descent sign at the bracket endpoint"];

expectZero[
  "M4 endpoint slope identity",
  (K0 + c rootMap) + Sqrt[disc[c]],
  branchAssumptions[c]
];
expectTrue[
  "M4 strict endpoint descent on simple-root branch",
  K0 + c rootMap < 0,
  branchAssumptions[c]
];

subbanner["M5. Bracket endpoints"];

tauLo = rootMap /. c -> cL;
tauHi = rootMap /. c -> cU;
expectZero["M5 tau_lo solves lower quadratic", q[cL, tauLo], branchAssumptions[cL]];
expectZero["M5 tau_hi solves upper quadratic", q[cU, tauHi], branchAssumptions[cU]];
expectZero[
  "M5 degenerate-envelope collapse",
  (tauLo /. cL -> c0) - (tauHi /. cU -> c0),
  baseAssumptions && disc[c0] > 0
];

subbanner["M6. Small-envelope width law"];

tauAt[cc_] := rootMap /. c -> cc;
width = tauAt[cmid + eta/2] - tauAt[cmid - eta/2];
widthSeries = Normal[Series[width, {eta, 0, 3}]];
leadingCoeff = Coefficient[widthSeries, eta, 1];
eta2Coeff = Coefficient[widthSeries, eta, 2];
tauMid = tauAt[cmid];
midAssumptions = baseAssumptions && disc[cmid] > 0;
expectZero[
  "M6 leading width coefficient",
  leadingCoeff - tauMid^2/(2 Sqrt[disc[cmid]]),
  midAssumptions
];
expectZero["M6 eta^2 coefficient cancellation", eta2Coeff, midAssumptions];
tau0 = H0/k;
expectZero[
  "M6 zero-mean width leading term",
  (leadingCoeff /. cmid -> 0) eta - tau0^2 eta/(2 k),
  H0 > 0 && k > 0 && Element[eta, Reals]
];

subbanner["M7. Turning-ray roots"];

turningAssumptions = H0 > 0 && a > 0 && Element[{H0, a}, Reals];
turningEquation = (H0 + 1/2 Klower tau^2 == 0) /. Klower -> -a;
turningRoots = stripConditional[tau /. Solve[turningEquation, tau]];
turningRoot = SelectFirst[
  turningRoots,
  TrueQ[Refine[# > 0, turningAssumptions]] &
];
If[MissingQ[turningRoot], fail["M7 positive turning-root selection", turningRoots]];
tauTp = Sqrt[2 H0/a];
Print["M7 turning roots = ", fmt[turningRoots]];
expectZero["M7 turning root solves quadratic", H0 - 1/2 a tauTp^2, turningAssumptions];
expectZero["M7 Solve-selected turning root", turningRoot - tauTp, turningAssumptions];
expectZero[
  "M7 turning derivative identity",
  D[tauTp, a] + tauTp/(2 a),
  turningAssumptions
];

subbanner["F3a. Pairwise ray-ordering theorem"];

pairwiseVars = {tauLoA, tauHiA, tauStarA, tauLoB, tauHiB, tauStarB};
pairwiseHypotheses = (
  tauLoA <= tauStarA <= tauHiA &&
  tauLoB <= tauStarB <= tauHiB &&
  tauHiA < tauLoB
);
pairwiseConclusion = tauStarA < tauStarB;
pairwiseTheorem = With[
  {vars = pairwiseVars, premise = pairwiseHypotheses, conclusion = pairwiseConclusion},
  Resolve[ForAll[vars, Implies[premise, conclusion]], Reals]
];
relaxedTheorem = With[
  {
    vars = pairwiseVars,
    premise = tauLoA <= tauStarA <= tauHiA && tauLoB <= tauStarB <= tauHiB,
    conclusion = pairwiseConclusion
  },
  Resolve[ForAll[vars, Implies[premise, conclusion]], Reals]
];
expectTrue["F3a pairwise certified interval ordering", pairwiseTheorem];
expectFalse["F3a ordering without separation is not tautological", relaxedTheorem];

subbanner["F3b. Local search-sieve admissibility test"];

monotoneAdmissible[h_, slope_, cLo_, cHi_, validT_, upperTau_] := Refine[
  h > 0 && slope < 0 && cLo <= cHi &&
  slope^2 - 2 cLo h >= 0 && slope^2 - 2 cHi h >= 0 &&
  upperTau <= validT
];

turningAdmissible[h_, slope_, cLo_, cHi_, validT_, upperTauTp_] := Refine[
  h > 0 && slope == 0 && cLo <= cHi && cHi < 0 && upperTauTp <= validT
];

localAdmissible[h_, slope_, cLo_, cHi_, validT_, upperTau_, upperTauTp_] := Refine[
  monotoneAdmissible[h, slope, cLo, cHi, validT, upperTau] ||
  turningAdmissible[h, slope, cLo, cHi, validT, upperTauTp]
];

monotoneCaseTau = FullSimplify[rootMap /. {H0 -> 1, k -> 3, c -> 1}];
turningCaseTau = Sqrt[2*2/1];
expectTrue[
  "F3b monotone admissible bracket",
  localAdmissible[1, -3, 0, 1, 1, monotoneCaseTau, 0]
];
expectTrue[
  "F3b turning admissible bracket",
  localAdmissible[2, 0, -3, -1, 3, 0, turningCaseTau]
];
expectFalse[
  "F3b single-clause tau_hi violation",
  localAdmissible[1, -3, 0, 1, 1, 2, 0]
];

banner["STAGE 206 MATHEMATICA AUDIT PASSED"];
Exit[0];
