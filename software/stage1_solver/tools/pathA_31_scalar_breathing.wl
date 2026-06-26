(* pathA_31 Gate 2 scalar-breathing reduction, Mathematica engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
boolText[x_] := If[TrueQ[x], "true", "false"];
quoteText[x_] := "'" <> StringReplace[If[StringQ[x], x, ToString[x, InputForm]], "'" -> "''"] <> "'";
numText[x_] := StringReplace[ToString[N[x, 17], InputForm], "*^" -> "e"];

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_31_scalar_breathing.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
sympyExprFile = FileNameJoin[{scratchDir, "pathA_31_sympy_exprs.wl"}];
yamlOut = FileNameJoin[{scratchDir, "pathA_31_mathematica_results.yaml"}];

If[! FileExistsQ[sympyExprFile], fail["missing SymPy expression export: " <> sympyExprFile]];
Get[sympyExprFile];

$Assumptions = L0 > 0 && beta > 0 && muEta > 0 && Tw > 0 && rAL > 0 &&
  rhoStar > 0 && Vp0 > 0 && ellC > 0 && kappa > 0 && chi > 0 &&
  sigmaA > 0 && sigmaL > 0 && Element[w, Reals];

kEta = FullSimplify[Tw beta^2];
src = FullSimplify[rhoStar Vp0/ellC];

alphaA = FullSimplify[
  DSolveValue[{-D[y[w], {w, 2}] + beta^2 y[w] == 0, y[0] == 1, y[L0] == 0}, y[w], w]
];
alphaL = FullSimplify[
  DSolveValue[{-D[y[w], {w, 2}] + beta^2 y[w] == 0, y[0] == 0, y[L0] == rAL}, y[w], w]
];

massAA = FullSimplify[4 Pi Integrate[muEta alphaA^2, {w, 0, L0}]];
massAL = FullSimplify[4 Pi Integrate[muEta alphaA alphaL, {w, 0, L0}]];
massLL = FullSimplify[4 Pi Integrate[muEta alphaL^2, {w, 0, L0}]];

stiffAA = FullSimplify[4 Pi Integrate[Tw D[alphaA, w]^2 + kEta alphaA^2, {w, 0, L0}]];
stiffAL = FullSimplify[4 Pi Integrate[Tw D[alphaA, w] D[alphaL, w] + kEta alphaA alphaL, {w, 0, L0}]];
stiffLL = FullSimplify[4 Pi Integrate[Tw D[alphaL, w]^2 + kEta alphaL^2, {w, 0, L0}]];
massDet = FullSimplify[Det[{{massAA, massAL}, {massAL, massLL}}]];
stiffDet = FullSimplify[Det[{{stiffAA, stiffAL}, {stiffAL, stiffLL}}]];

forceDistA = FullSimplify[4 Pi Integrate[src alphaA, {w, 0, L0}]];
forceDistL = FullSimplify[4 Pi Integrate[src alphaL, {w, 0, L0}]];
vConf = (Vp0/ellC) (r - (R0[w] + qa alphaA + qL alphaL));
forceLegacyA = FullSimplify[-4 Pi Integrate[rhoStar (D[vConf, qa] /. {qa -> 0, qL -> 0}), {w, 0, L0}]];
forceLegacyL = FullSimplify[-4 Pi Integrate[rhoStar (D[vConf, qL] /. {qa -> 0, qL -> 0}), {w, 0, L0}]];

wrongZeroAlpha = 0;
wrongZeroMAA = FullSimplify[4 Pi Integrate[muEta wrongZeroAlpha^2, {w, 0, L0}]];
wrongZeroMAL = FullSimplify[4 Pi Integrate[muEta wrongZeroAlpha alphaL, {w, 0, L0}]];
wrongZeroMDet = FullSimplify[Det[{{wrongZeroMAA, wrongZeroMAL}, {wrongZeroMAL, massLL}}]];
wrongConstMAA = FullSimplify[4 Pi Integrate[muEta, {w, 0, L0}]];
wrongConstMAL = FullSimplify[4 Pi Integrate[muEta alphaL, {w, 0, L0}]];
wrongConstMDet = FullSimplify[Det[{{wrongConstMAA, wrongConstMAL}, {wrongConstMAL, massLL}}]];
wrongZeroFA = FullSimplify[4 Pi Integrate[src wrongZeroAlpha, {w, 0, L0}]];
wrongConstFA = FullSimplify[4 Pi Integrate[src, {w, 0, L0}]];

egeom = 1/2 kappa (dL - chi da)^2 + 1/2 sigmaA da^2 + 1/2 sigmaL dL^2;
legacyH = FullSimplify[D[egeom, {{da, dL}, 2}]];

massMatrix = {{massAA, massAL}, {massAL, massLL}};
stiffMatrix = {{stiffAA, stiffAL}, {stiffAL, stiffLL}};
zeroPattern[m_] := Map[TrueQ[FullSimplify[# == 0]] &, m, {2}];
legacyRank = If[TrueQ[FullSimplify[Det[legacyH] != 0]], 2, MatrixRank[legacyH]];
legacyZeroPattern = zeroPattern[legacyH];
legacySymmetric = TrueQ[FullSimplify[legacyH == Transpose[legacyH]]];
legacyOffdiagNegative = TrueQ[FullSimplify[legacyH[[1, 2]] < 0]];
structureB = beta L0;
structureMAAPositiveForm = FullSimplify[2 Pi muEta (Sinh[structureB] Cosh[structureB] - structureB)/(beta Sinh[structureB]^2)];
structureMDetPositiveForm = FullSimplify[4 Pi^2 muEta^2 rAL^2 (Sinh[structureB]^2 - structureB^2)/(beta^2 Sinh[structureB]^2)];
structureMLeadingPositive = TrueQ[FullSimplify[massAA - structureMAAPositiveForm == 0]];
structureMDetPositive = TrueQ[FullSimplify[massDet - structureMDetPositiveForm == 0]];
structureMPosdef = TrueQ[structureMLeadingPositive && structureMDetPositive];
structureKSymmetric = TrueQ[FullSimplify[stiffMatrix == Transpose[stiffMatrix]]];
structureKOffdiagNegative = TrueQ[FullSimplify[stiffAL < 0]];
structureKRank = If[TrueQ[FullSimplify[stiffDet != 0]], 2, MatrixRank[stiffMatrix]];
structureKRankMatchesLegacy = TrueQ[structureKRank == legacyRank];
structureKZeroPatternMatchesLegacy = TrueQ[zeroPattern[stiffMatrix] == legacyZeroPattern];
structureKStructureOk = TrueQ[
  structureKSymmetric && structureKOffdiagNegative && structureKRankMatchesLegacy &&
    structureKZeroPatternMatchesLegacy && legacySymmetric && legacyOffdiagNegative
];
structureMProbeMatrix = {{-massAA, massAL}, {massAL, massLL}};
structureMProbeLeadingPositive = TrueQ[FullSimplify[structureMProbeMatrix[[1, 1]] - structureMAAPositiveForm == 0]];
structureMProbeDetPositive = TrueQ[FullSimplify[Det[structureMProbeMatrix] - structureMDetPositiveForm == 0]];
structureMProbePosdef = TrueQ[structureMProbeLeadingPositive && structureMProbeDetPositive];
structureKProbeMatrix = {{stiffAA, -stiffAL}, {-stiffAL, stiffLL}};
structureKProbeOffdiagNegative = TrueQ[FullSimplify[structureKProbeMatrix[[1, 2]] < 0]];
structureKProbeRank = If[TrueQ[FullSimplify[Det[structureKProbeMatrix] != 0]], 2, MatrixRank[structureKProbeMatrix]];
structureKProbeStructureOk = TrueQ[
  TrueQ[FullSimplify[structureKProbeMatrix == Transpose[structureKProbeMatrix]]] &&
    structureKProbeOffdiagNegative && TrueQ[structureKProbeRank == legacyRank] &&
    TrueQ[zeroPattern[structureKProbeMatrix] == legacyZeroPattern] &&
    legacySymmetric && legacyOffdiagNegative
];

checks = <|
  "alpha_a" -> TrueQ[FullSimplify[alphaA - sympyAlphaA == 0]],
  "alpha_L" -> TrueQ[FullSimplify[alphaL - sympyAlphaL == 0]],
  "M_aa" -> TrueQ[FullSimplify[massAA - sympyMassAA == 0]],
  "M_aL" -> TrueQ[FullSimplify[massAL - sympyMassAL == 0]],
  "M_LL" -> TrueQ[FullSimplify[massLL - sympyMassLL == 0]],
  "K_aa" -> TrueQ[FullSimplify[stiffAA - sympyStiffAA == 0]],
  "K_aL" -> TrueQ[FullSimplify[stiffAL - sympyStiffAL == 0]],
  "K_LL" -> TrueQ[FullSimplify[stiffLL - sympyStiffLL == 0]],
  "M_det" -> TrueQ[FullSimplify[massDet - sympyMassDet == 0]],
  "K_det" -> TrueQ[FullSimplify[stiffDet - sympyStiffDet == 0]],
  "force_dist_a" -> TrueQ[FullSimplify[forceDistA - sympyForceDistA == 0]],
  "force_dist_L" -> TrueQ[FullSimplify[forceDistL - sympyForceDistL == 0]],
  "force_legacy_a" -> TrueQ[FullSimplify[forceLegacyA - sympyForceLegacyA == 0]],
  "force_legacy_L" -> TrueQ[FullSimplify[forceLegacyL - sympyForceLegacyL == 0]],
  "wrong_zero_M_det" -> TrueQ[FullSimplify[wrongZeroMDet - sympyWrongZeroMDet == 0]],
  "wrong_const_M_det" -> TrueQ[FullSimplify[wrongConstMDet - sympyWrongConstMDet == 0]],
  "wrong_zero_F_a" -> TrueQ[FullSimplify[wrongZeroFA - sympyWrongZeroFA == 0]],
  "wrong_const_F_a" -> TrueQ[FullSimplify[wrongConstFA - sympyWrongConstFA == 0]],
  "legacy_H_aa" -> TrueQ[FullSimplify[legacyH[[1, 1]] - sympyLegacyHAA == 0]],
  "legacy_H_aL" -> TrueQ[FullSimplify[legacyH[[1, 2]] - sympyLegacyHAL == 0]],
  "legacy_H_LL" -> TrueQ[FullSimplify[legacyH[[2, 2]] - sympyLegacyHLL == 0]],
  "structure_M_posdef" -> TrueQ[FullSimplify[structureMPosdef == sympyMPosdef]],
  "structure_K_symmetric" -> TrueQ[FullSimplify[structureKSymmetric == sympyKSymmetric]],
  "structure_K_offdiag_negative" -> TrueQ[FullSimplify[structureKOffdiagNegative == sympyKOffdiagNegative]],
  "structure_K_structure_ok" -> TrueQ[FullSimplify[structureKStructureOk == sympyKStructureOk]],
  "structure_K_rank" -> TrueQ[FullSimplify[structureKRank - sympyKRank == 0]],
  "structure_probe_M_posdef" -> TrueQ[FullSimplify[structureMProbePosdef == sympyStructureProbeMPosdef]],
  "structure_probe_K_structure_ok" -> TrueQ[FullSimplify[structureKProbeStructureOk == sympyStructureProbeKStructureOk]]
|>;

If[! And @@ Values[checks],
  Print["Agreement checks: ", checks];
  fail["Mathematica symbolic route disagrees with SymPy artifacts"]
];

ClearAll[profileFunctions, galerkinRow];
profileFunctions[b_?NumericQ, profile_String] := Module[
  {aa, daa, al, dal},
  al = Function[x, Sinh[b x]/Sinh[b]];
  dal = Function[x, b Cosh[b x]/Sinh[b]];
  Switch[profile,
    "baseline",
      aa = Function[x, Sinh[b (1 - x)]/Sinh[b]];
      daa = Function[x, -b Cosh[b (1 - x)]/Sinh[b]],
    "degenerate_zero",
      aa = Function[x, 0.];
      daa = Function[x, 0.],
    "constant_one",
      aa = Function[x, 1.];
      daa = Function[x, 0.],
    _, fail["unknown profile " <> profile]
  ];
  {{aa, al}, {daa, dal}}
];

galerkinRow[b_?NumericQ, nModes_Integer, profile_String] := Module[
  {pair, funcs, ders, kList, mFull, mass, stiff, active, ma, ka, vals, vecs,
   ord, sub, selector, gram, pinv, overlaps, gap, min12, i, j, n, k},
  pair = profileFunctions[b, profile];
  funcs = pair[[1]];
  ders = pair[[2]];
  kList = Table[(n - 1/2) Pi, {n, 1, nModes}];
  mFull = 2 + nModes;
  mass = ConstantArray[0., {mFull, mFull}];
  stiff = ConstantArray[0., {mFull, mFull}];
  Do[
    mass[[i, j]] = NIntegrate[funcs[[i]][x] funcs[[j]][x], {x, 0, 1},
      AccuracyGoal -> 11, PrecisionGoal -> 11, MaxRecursion -> 12];
    stiff[[i, j]] = NIntegrate[ders[[i]][x] ders[[j]][x] + b^2 funcs[[i]][x] funcs[[j]][x], {x, 0, 1},
      AccuracyGoal -> 11, PrecisionGoal -> 11, MaxRecursion -> 12];
    mass[[j, i]] = mass[[i, j]];
    stiff[[j, i]] = stiff[[i, j]],
    {i, 1, 2}, {j, i, 2}
  ];
  Do[
    k = kList[[n]];
    Do[
      mass[[i, 2 + n]] = NIntegrate[funcs[[i]][x] Sin[k x], {x, 0, 1},
        AccuracyGoal -> 11, PrecisionGoal -> 11, MaxRecursion -> 12];
      stiff[[i, 2 + n]] = NIntegrate[ders[[i]][x] k Cos[k x] + b^2 funcs[[i]][x] Sin[k x], {x, 0, 1},
        AccuracyGoal -> 11, PrecisionGoal -> 11, MaxRecursion -> 12];
      mass[[2 + n, i]] = mass[[i, 2 + n]];
      stiff[[2 + n, i]] = stiff[[i, 2 + n]],
      {i, 1, 2}
    ];
    mass[[2 + n, 2 + n]] = 1/2;
    stiff[[2 + n, 2 + n]] = (k^2 + b^2)/2,
    {n, 1, nModes}
  ];
  mass = N[(mass + Transpose[mass])/2, 20];
  stiff = N[(stiff + Transpose[stiff])/2, 20];
  active = Select[Range[mFull], mass[[#, #]] > 10^-13 &];
  ma = mass[[active, active]];
  ka = stiff[[active, active]];
  {vals, vecs} = Eigensystem[{ka, ma}];
  ord = Ordering[vals];
  vals = vals[[ord]];
  vecs = vecs[[ord]];
  sub = Flatten[Position[active, 1 | 2]];
  selector = ConstantArray[0., {Length[active], Length[sub]}];
  Do[selector[[sub[[j]], j]] = 1., {j, 1, Length[sub]}];
  gram = Transpose[selector].ma.selector;
  pinv = PseudoInverse[gram, Tolerance -> 10^-12];
  overlaps = Table[
    Module[{v = vecs[[j]], norm, coeff, pnorm},
      norm = Re[v.ma.v];
      coeff = pinv.Transpose[selector].ma.v;
      pnorm = Re[coeff.gram.coeff];
      Sqrt[Max[0., Min[1., N[pnorm/norm]]]]
    ],
    {j, 1, 2}
  ];
  min12 = Min[vals[[1]], vals[[2]]];
  gap = (vals[[3]] - vals[[2]])/vals[[2]];
  N[{b, overlaps[[1]], overlaps[[2]], min12, gap}, 16]
];

mmaSweep = galerkinRow[#, 16, "baseline"] & /@ sympyBetaSweep[[All, 1]];
mmaSelected = galerkinRow[sympySelectedOverlap[[1]], 16, "baseline"];
mmaWrongZero = (galerkinRow[sympySelectedOverlap[[1]], 16, "degenerate_zero"])[[2 ;; 3]];
mmaWrongConst = (galerkinRow[sympySelectedOverlap[[1]], 16, "constant_one"])[[2 ;; 3]];

maxSweepDelta = Max[Abs[Flatten[mmaSweep - sympyBetaSweep]]];
selectedDelta = Max[Abs[mmaSelected - sympySelectedOverlap]];
wrongZeroDelta = Max[Abs[mmaWrongZero - sympyWrongZeroOverlap]];
wrongConstDelta = Max[Abs[mmaWrongConst - sympyWrongConstOverlap]];
maxNumericDelta = Max[maxSweepDelta, selectedDelta, wrongZeroDelta, wrongConstDelta];
numericChecks = <|
  "beta_sweep" -> TrueQ[maxSweepDelta < 10^-6],
  "selected_overlap" -> TrueQ[selectedDelta < 10^-6],
  "wrong_zero_overlap" -> TrueQ[wrongZeroDelta < 10^-6],
  "wrong_const_overlap" -> TrueQ[wrongConstDelta < 10^-6]
|>;

If[! And @@ Values[numericChecks],
  Print["Numeric deltas: ", {maxSweepDelta, selectedDelta, wrongZeroDelta, wrongConstDelta}];
  fail["Mathematica numeric Galerkin route disagrees with SymPy/SciPy artifacts"]
];

yaml = StringRiffle[{
  "schema: pathA_31_scalar_breathing_mathematica/v3",
  "route: native_DSolveValue_BVP_Integrate_plus_NIntegrate_generalized_Eigensystem",
  "sympy_expression_digest: " <> quoteText[sympyExpressionDigest],
  "alpha_a: " <> quoteText[alphaA],
  "alpha_L: " <> quoteText[alphaL],
  "M_AB:",
  "  aa: " <> quoteText[massAA],
  "  aL: " <> quoteText[massAL],
  "  LL: " <> quoteText[massLL],
  "K_AB:",
  "  aa: " <> quoteText[stiffAA],
  "  aL: " <> quoteText[stiffAL],
  "  LL: " <> quoteText[stiffLL],
  "selected_overlap:",
  "  beta_L0: " <> numText[mmaSelected[[1]]],
  "  o_1: " <> numText[mmaSelected[[2]]],
  "  o_2: " <> numText[mmaSelected[[3]]],
  "  min_omega12_squared: " <> numText[mmaSelected[[4]]],
  "  gap: " <> numText[mmaSelected[[5]]],
  "max_numeric_abs_delta: " <> numText[maxNumericDelta],
  "engine_agreement:",
  "  status: pass",
  "  checks:",
  "    alpha_a: " <> boolText[checks["alpha_a"]],
  "    alpha_L: " <> boolText[checks["alpha_L"]],
  "    M_aa: " <> boolText[checks["M_aa"]],
  "    M_aL: " <> boolText[checks["M_aL"]],
  "    M_LL: " <> boolText[checks["M_LL"]],
  "    K_aa: " <> boolText[checks["K_aa"]],
  "    K_aL: " <> boolText[checks["K_aL"]],
  "    K_LL: " <> boolText[checks["K_LL"]],
  "    M_det: " <> boolText[checks["M_det"]],
  "    K_det: " <> boolText[checks["K_det"]],
  "    force_dist_a: " <> boolText[checks["force_dist_a"]],
  "    force_dist_L: " <> boolText[checks["force_dist_L"]],
  "    force_legacy_a: " <> boolText[checks["force_legacy_a"]],
  "    force_legacy_L: " <> boolText[checks["force_legacy_L"]],
  "    wrong_zero_M_det: " <> boolText[checks["wrong_zero_M_det"]],
  "    wrong_const_M_det: " <> boolText[checks["wrong_const_M_det"]],
  "    wrong_zero_F_a: " <> boolText[checks["wrong_zero_F_a"]],
  "    wrong_const_F_a: " <> boolText[checks["wrong_const_F_a"]],
  "    legacy_H_aa: " <> boolText[checks["legacy_H_aa"]],
  "    legacy_H_aL: " <> boolText[checks["legacy_H_aL"]],
  "    legacy_H_LL: " <> boolText[checks["legacy_H_LL"]],
  "    structure_M_posdef: " <> boolText[checks["structure_M_posdef"]],
  "    structure_K_symmetric: " <> boolText[checks["structure_K_symmetric"]],
  "    structure_K_offdiag_negative: " <> boolText[checks["structure_K_offdiag_negative"]],
  "    structure_K_structure_ok: " <> boolText[checks["structure_K_structure_ok"]],
  "    structure_K_rank: " <> boolText[checks["structure_K_rank"]],
  "    structure_probe_M_posdef: " <> boolText[checks["structure_probe_M_posdef"]],
  "    structure_probe_K_structure_ok: " <> boolText[checks["structure_probe_K_structure_ok"]],
  "  numeric_checks:",
  "    beta_sweep: " <> boolText[numericChecks["beta_sweep"]],
  "    selected_overlap: " <> boolText[numericChecks["selected_overlap"]],
  "    wrong_zero_overlap: " <> boolText[numericChecks["wrong_zero_overlap"]],
  "    wrong_const_overlap: " <> boolText[numericChecks["wrong_const_overlap"]],
  "structure:",
  "  M_posdef: " <> boolText[structureMPosdef],
  "  K_symmetric: " <> boolText[structureKSymmetric],
  "  K_offdiag_negative: " <> boolText[structureKOffdiagNegative],
  "  K_structure_ok: " <> boolText[structureKStructureOk],
  "  K_rank: " <> ToString[structureKRank],
  "  structure_probe_M_posdef: " <> boolText[structureMProbePosdef],
  "  structure_probe_K_structure_ok: " <> boolText[structureKProbeStructureOk],
  "degenerate_counterfactual:",
  "  derived: true",
  "  wrong_zero_M_det: " <> quoteText[wrongZeroMDet],
  "  wrong_zero_F_a: " <> quoteText[wrongZeroFA]
}, "\n"] <> "\n";

Export[yamlOut, yaml, "Text"];
Print["verdict: mathematica_engine_pass"];
Print["yaml: ", yamlOut];
Exit[0];
