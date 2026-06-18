#!/usr/bin/env wolframscript

(* MT15-06 / M1c-prep cross-engine producer.
   Imports a self-consistent torch WP1 background and re-runs the derived
   BdG/wall plus Spike-2 transfer chain with the gauge fields active.
   Target-blind machinery only: no Gate A, no frozen physical packet. *)

$HistoryLength = 0;
ClearAll["Global`*"];

scriptPath =
  If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
    $InputFileName,
    FileNameJoin[{Directory[], "software/stage1_solver/mathematica/mt15_06_m1c_prep_crossengine.wls"}]
  ];
scriptDir = DirectoryName[scriptPath];
repoRoot = AbsoluteFileName[FileNameJoin[{scriptDir, "..", "..", ".."}]];
runDir = FileNameJoin[{scriptDir, "runs", "mt15_06_m1c_prep_crossengine"}];
If[! DirectoryQ[runDir], CreateDirectory[runDir, CreateIntermediateDirectories -> True]];

defaultBackgroundPath = FileNameJoin[{
  repoRoot, "software", "stage1_solver", "runs", "m1c_prep_background_export",
  "background_bundles", "m1c_prep_wp1_background_latest.json"
}];
backgroundPath =
  If[StringQ[Environment["M1C_BACKGROUND_JSON"]] && StringLength[Environment["M1C_BACKGROUND_JSON"]] > 0,
    Environment["M1C_BACKGROUND_JSON"],
    If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], defaultBackgroundPath]
  ];

physicalFreezePath =
  If[StringQ[Environment["M1C_PHYSICAL_FREEZE_JSON"]] && StringLength[Environment["M1C_PHYSICAL_FREEZE_JSON"]] > 0,
    Environment["M1C_PHYSICAL_FREEZE_JSON"],
    ""
  ];
physicalMode = StringLength[physicalFreezePath] > 0;
physicalOutputDir =
  If[StringQ[Environment["M1C_FROZEN_OUTPUT_DIR"]] && StringLength[Environment["M1C_FROZEN_OUTPUT_DIR"]] > 0,
    Environment["M1C_FROZEN_OUTPUT_DIR"],
    runDir
  ];
If[physicalMode && ! DirectoryQ[physicalOutputDir],
  CreateDirectory[physicalOutputDir, CreateIntermediateDirectories -> True]
];

packetPath = If[physicalMode,
  FileNameJoin[{physicalOutputDir, "m1c_v2_22b_physical_frozen_packet.json"}],
  FileNameJoin[{runDir, "mt15_06_v2_22b_m1c_prep_direct_derived_packet.json"}]
];
diagnosticsPath = If[physicalMode,
  FileNameJoin[{physicalOutputDir, "m1c_physical_derived_diagnostics.json"}],
  FileNameJoin[{runDir, "mt15_06_m1c_prep_crossengine_diagnostics.json"}]
];
reportPath = If[physicalMode,
  FileNameJoin[{physicalOutputDir, "m1c_physical_derived_report.md"}],
  FileNameJoin[{scriptDir, "mt15_06_m1c_prep_crossengine_report.md"}]
];

ClearAll[
  finiteQ, toReal, toJsonNumber, toJsonVector, toJsonMatrix, keySortRecursive,
  contentHash, sourceCommitString, sha256File, formatScientific, formatNumber,
  markdownTable, smootherStep, trapz, weightedInner, weightedNorm,
  normalizeProfile, bracketIndex, interp1, interp2, localInterpolationIndices,
  lagrangeInterpolate, interp1Order, interp2Order, overlapLength,
  radialOverlapVolume, projectScalarCells, pointVector, sourceMassNo4Pi,
  cellIndex, rotateByMaxComponent, diagSparse, firstDerivativeMatrix,
  secondDerivativeMatrix, blockRow, weightedGram, boundaryMatrix,
  buildMaxwellSystem, buildWallMode, linearizedCurrentFormula,
  nonlinearNumberCurrent, buildMatterSource, currentFrechetCheck,
  selfEnergySolve, fitEvenCoefficients, estimateOrderTriplet,
  fieldNormFromComponents, makeGaugeVector, pureGaugeCheck,
  basisInvarianceCheck, v209Regression, hardCapContrast, transferForGrid,
  runBdgWallSector, extractionForBackgroundInterpolationOrder,
  coefficientRelativeDeltaTable
];

finiteQ[x_] := Quiet[Check[NumberQ[N[x]] && TrueQ[Abs[N[x]] < Infinity], False]];
toReal[x_] := Module[{n = Quiet[Check[N[x], $Failed]]},
  If[n === $Failed || ! NumericQ[n],
    Print["ERROR: nonnumeric value reached JSON export: ", InputForm[x]];
    Exit[4]
  ];
  N[Re[n], 17]
];
toJsonNumber[x_] := toReal[x];
toJsonVector[x_List] := toJsonNumber /@ x;
toJsonMatrix[x_List] := Map[toJsonNumber, x, {2}];

keySortRecursive[x_Association] := Association @ KeyValueMap[#1 -> keySortRecursive[#2] &, KeySort[x]];
keySortRecursive[x_List] := keySortRecursive /@ x;
keySortRecursive[x_] := x;
contentHash[x_Association] := IntegerString[Hash[ToString[InputForm[keySortRecursive[KeyDrop[x, {"content_hash"}]]]], "SHA256"], 16, 64];

sourceCommitString[] := Module[{proc, out},
  proc = Quiet[RunProcess[{"git", "-C", repoRoot, "rev-parse", "--short", "HEAD"}, "StandardOutput"]];
  out = If[StringQ[proc], StringTrim[proc], ""];
  If[StringLength[out] > 0, out, "unavailable"]
];

sha256File[path_String] := Module[{text, hash},
  If[! FileExistsQ[path], Return["missing"]];
  text = Quiet[Check[Import[path, "Text"], $Failed]];
  If[text === $Failed || ! StringQ[text], Return["unavailable"]];
  hash = Hash[text, "SHA256"];
  IntegerString[hash, 16, 64]
];

formatScientific[x_] := StringReplace[ToString[FortranForm[N[x, 8]]], " " -> ""];
formatNumber[x_] := ToString[FortranForm[N[x, 17]]];
markdownTable[headers_List, rows_List] := StringRiffle[
  Join[
    {"| " <> StringRiffle[headers, " | "] <> " |"},
    {"| " <> StringRiffle[ConstantArray["---", Length[headers]], " | "] <> " |"},
    ("| " <> StringRiffle[ToString /@ #, " | "] <> " |") & /@ rows
  ],
  "\n"
];

If[! FileExistsQ[backgroundPath],
  Print["ERROR: background JSON not found: ", backgroundPath];
  Exit[2]
];
background = Quiet[Check[Import[backgroundPath, "RawJSON"], $Failed]];
If[background === $Failed || ! AssociationQ[background],
  Print["ERROR: failed to import background JSON: ", backgroundPath];
  Exit[3]
];
If[! MemberQ[{"stage1_m1c_prep_wp1_background/v1", "stage1_m1c_physical_wp1_background/v1"}, background["schema"]],
  Print["ERROR: unexpected background schema: ", background["schema"]];
  Exit[3]
];
physicalFreeze = <||>;
If[physicalMode,
  If[! FileExistsQ[physicalFreezePath],
    Print["ERROR: physical freeze JSON not found: ", physicalFreezePath];
    Exit[3]
  ];
  physicalFreeze = Quiet[Check[Import[physicalFreezePath, "RawJSON"], $Failed]];
  If[physicalFreeze === $Failed || ! AssociationQ[physicalFreeze],
    Print["ERROR: failed to import physical freeze JSON: ", physicalFreezePath];
    Exit[3]
  ];
  If[! TrueQ[background["m1c_physical_frozen_run"]],
    Print["ERROR: physical mode requires a physical frozen WP1 background"];
    Exit[3]
  ];
  If[background["freeze_hash"] =!= background["freeze_sheet_sha256"],
    Print["ERROR: physical background freeze hash fields disagree"];
    Exit[3]
  ];
];

cfg = background["config"]["branch"];
gridBg = background["grid"];
fieldsBg = background["fields"];
derivedBg = background["derived"];

rFacesSrc = N[gridBg["r_faces"], 17];
rCentersSrc = N[gridBg["r_centers"], 17];
wFacesSrc = N[gridBg["w_faces"], 17];
wCentersSrc = N[gridBg["w_centers"], 17];
rhoMatrix = N[derivedBg["rho0"], 17];
psiRMatrix = N[fieldsBg["psi_R0"], 17];
psiIMatrix = N[fieldsBg["psi_I0"], 17];
a0Matrix = N[fieldsBg["A_00"], 17];
arMatrix = N[fieldsBg["A_r0"], 17];
awMatrix = N[fieldsBg["A_w0"], 17];
r0ValuesSrc = N[derivedBg["R0_w"], 17];
zValuesSrc = N[derivedBg["Z_w"], 17];

nrSrc = Length[rCentersSrc];
nwSrc = Length[wCentersSrc];
wMin = toReal[gridBg["w_min"]];
wMax = toReal[gridBg["w_max"]];
wLength = wMax - wMin;
rDomainMax = toReal[gridBg["r_max"]];

hbar = toReal[cfg["hbar"]];
particleMass = toReal[cfg["particle_mass"]];
gaugeCharge = toReal[cfg["gauge_charge"]];
mu0 = toReal[cfg["mu0"]];
xiGauge = toReal[cfg["xi"]];
eosK = toReal[Last[cfg["continuation_K_values"]]];
chemicalPotential = toReal[background["solver"]["chemical_potential"]];
radialWallStrength = toReal[cfg["radial_wall_strength"]];
axialTrapStrength = toReal[cfg["axial_trap_strength"]];

ell = 2;
Lharm = ell (ell + 1);
GConst = If[physicalMode, toReal[physicalFreeze["constants"]["G"]], 1.0];
csConst = If[physicalMode, toReal[physicalFreeze["constants"]["c_s"]], 1.0];
cConst = If[physicalMode, toReal[physicalFreeze["constants"]["c"]], 2.5];
aConst = If[physicalMode, toReal[physicalFreeze["geometry"]["a"]], toReal[First[r0ValuesSrc]]];
mhat0Const = If[physicalMode, toReal[physicalFreeze["constants"]["mhat0"]], 0.85];
SPortConst = If[physicalMode, toReal[physicalFreeze["constants"]["S_port"]], 1.2];
thetaTailConst = If[physicalMode, toReal[physicalFreeze["constants"]["theta_tail"]], 1.0];
phaseRegularizer = 0.05;
rExitConst = If[physicalMode, toReal[physicalFreeze["geometry"]["R_exit"]], toReal[Last[r0ValuesSrc]]];

smootherStep[x_] := 3 x^2 - 2 x^3;
trapz[xs_List, ys_List] /; Length[xs] == Length[ys] && Length[xs] >= 2 :=
  Total[Table[0.5 (ys[[i]] + ys[[i + 1]]) (xs[[i + 1]] - xs[[i]]), {i, 1, Length[xs] - 1}]];
weightedInner[xs_List, weights_List, f_List, g_List] := trapz[xs, weights f g];
weightedNorm[xs_List, weights_List, f_List] := weightedInner[xs, weights, f, f];
normalizeProfile[xs_List, weights_List, f_List, label_String] := Module[{norm = weightedNorm[xs, weights, f]},
  If[! NumericQ[norm] || norm <= 0,
    Print["ERROR: nonpositive norm while normalizing ", label, ": ", norm];
    Exit[5]
  ];
  N[f/Sqrt[norm], 17]
];

bracketIndex[xs_List, x_?NumericQ] := Module[{n = Length[xs], k},
  If[n < 2, Return[1]];
  If[x <= xs[[1]], Return[1]];
  If[x >= xs[[-1]], Return[n - 1]];
  k = FirstCase[Range[n - 1], i_ /; xs[[i]] <= x <= xs[[i + 1]], n - 1];
  k
];

interp1[xs_List, ys_List, x_?NumericQ] := Module[{i, x1, x2, t},
  i = bracketIndex[xs, x];
  x1 = xs[[i]]; x2 = xs[[i + 1]];
  t = If[x2 == x1, 0.0, (x - x1)/(x2 - x1)];
  N[(1.0 - t) ys[[i]] + t ys[[i + 1]], 17]
];

interp2[m_List, r_?NumericQ, w_?NumericQ] := Module[{i, j, r1, r2, w1, w2, tr, tw},
  i = bracketIndex[rCentersSrc, r];
  j = bracketIndex[wCentersSrc, w];
  r1 = rCentersSrc[[i]]; r2 = rCentersSrc[[i + 1]];
  w1 = wCentersSrc[[j]]; w2 = wCentersSrc[[j + 1]];
  tr = If[r2 == r1, 0.0, (r - r1)/(r2 - r1)];
  tw = If[w2 == w1, 0.0, (w - w1)/(w2 - w1)];
  N[
    (1.0 - tr) (1.0 - tw) m[[i, j]] +
      tr (1.0 - tw) m[[i + 1, j]] +
      (1.0 - tr) tw m[[i, j + 1]] +
      tr tw m[[i + 1, j + 1]],
    17
  ]
];

localInterpolationIndices[xs_List, x_?NumericQ, order_Integer] := Module[
  {n = Length[xs], count, i, start},
  count = Min[n, Max[2, order + 1]];
  If[n <= count, Return[Range[n]]];
  i = bracketIndex[xs, x];
  start = i - Floor[(count - 1)/2];
  start = Max[1, Min[start, n - count + 1]];
  Range[start, start + count - 1]
];

lagrangeInterpolate[xs_List, ys_List, x_?NumericQ] := Module[{n = Length[xs]},
  Sum[
    ys[[i]] Product[
      If[j == i, 1.0, (x - xs[[j]])/(xs[[i]] - xs[[j]])],
      {j, 1, n}
    ],
    {i, 1, n}
  ]
];

interp1Order[xs_List, ys_List, x_?NumericQ, order_Integer] := Module[{idx},
  idx = localInterpolationIndices[xs, x, order];
  N[lagrangeInterpolate[xs[[idx]], ys[[idx]], x], 17]
];

interp2Order[m_List, r_?NumericQ, w_?NumericQ, order_Integer] := Module[
  {rIdx, wIdx, rNodes, wNodes, valuesAtW},
  rIdx = localInterpolationIndices[rCentersSrc, r, order];
  wIdx = localInterpolationIndices[wCentersSrc, w, order];
  rNodes = rCentersSrc[[rIdx]];
  wNodes = wCentersSrc[[wIdx]];
  valuesAtW = Table[
    interp1Order[rNodes, m[[rIdx, wIdx[[jw]]]], r, order],
    {jw, 1, Length[wIdx]}
  ];
  N[interp1Order[wNodes, valuesAtW, w, order], 17]
];

r0Geometry[w_] := If[physicalMode,
  Module[{x = Clip[(w - wMin)/wLength, {0.0, 1.0}]},
    aConst + (rExitConst - aConst) (3.0 x^2 - 2.0 x^3)
  ],
  interp1[wCentersSrc, r0ValuesSrc, w]
];
zWeight[w_] := interp1[wCentersSrc, zValuesSrc, w];

overlapLength[a1_, a2_, b1_, b2_] := Max[0.0, Min[a2, b2] - Max[a1, b1]];
radialOverlapVolume[a1_, a2_, b1_, b2_] := Module[{lo = Max[a1, b1], hi = Min[a2, b2]},
  If[hi <= lo, 0.0, (hi^3 - lo^3)/3.0]
];

projectScalarCells[m_List, rFacesT_List, wFacesT_List] := Module[
  {nrT = Length[rFacesT] - 1, nwT = Length[wFacesT] - 1, values, weights, mass, vol, acc, ovR, ovW},
  values = Flatten[
    Table[
      vol = (rFacesT[[i + 1]]^3 - rFacesT[[i]]^3)/3.0 (wFacesT[[j + 1]] - wFacesT[[j]]);
      acc = Sum[
        ovR = radialOverlapVolume[rFacesT[[i]], rFacesT[[i + 1]], rFacesSrc[[is]], rFacesSrc[[is + 1]]];
        ovW = overlapLength[wFacesT[[j]], wFacesT[[j + 1]], wFacesSrc[[js]], wFacesSrc[[js + 1]]];
        m[[is, js]] ovR ovW,
        {is, 1, nrSrc}, {js, 1, nwSrc}
      ];
      If[vol > 0, acc/vol, 0.0],
      {j, 1, nwT}, {i, 1, nrT}
    ]
  ];
  weights = Flatten[
    Table[
      (rFacesT[[i + 1]]^3 - rFacesT[[i]]^3)/3.0 (wFacesT[[j + 1]] - wFacesT[[j]]),
      {j, 1, nwT}, {i, 1, nrT}
    ]
  ];
  mass = Total[values weights];
  <|"values" -> N[values, 17], "weights" -> N[weights, 17], "mass" -> N[mass, 17]|>
];

sourceMassNo4Pi[m_List] := Sum[
  m[[is, js]] ((rFacesSrc[[is + 1]]^3 - rFacesSrc[[is]]^3)/3.0) (wFacesSrc[[js + 1]] - wFacesSrc[[js]]),
  {is, 1, nrSrc}, {js, 1, nwSrc}
];

pointVector[m_List, rVals_List, wVals_List] := N[MapThread[interp2[m, #1, #2] &, {rVals, wVals}], 17];

sourceRhoMass = N[sourceMassNo4Pi[rhoMatrix], 17];
sourceA0Linf = toReal[derivedBg["A_00_linf"]];
sourceArLinf = toReal[derivedBg["A_r0_linf"]];
sourceAwLinf = toReal[derivedBg["A_w0_linf"]];

fieldSelfMatchMax = Max[
  Flatten[
    Table[
      Abs[interp2[psiRMatrix, rCentersSrc[[i]], wCentersSrc[[j]]] - psiRMatrix[[i, j]]],
      {i, 1, nrSrc}, {j, 1, nwSrc}
    ]
  ]
];
r0SelfMatchMax = Max[Abs[MapThread[interp1[wCentersSrc, r0ValuesSrc, #1] - #2 &, {wCentersSrc, r0ValuesSrc}]]];

diagSparse[v_List] := SparseArray[Band[{1, 1}] -> N[v, 17], {Length[v], Length[v]}];
firstDerivativeMatrix[n_Integer, h_?NumericQ] := Module[{rules},
  rules = Join[
    Table[If[i < n, {i, i + 1} -> 1.0/(2.0 h), Nothing], {i, 1, n}],
    Table[If[i > 1, {i, i - 1} -> -1.0/(2.0 h), Nothing], {i, 1, n}]
  ];
  SparseArray[rules, {n, n}]
];
secondDerivativeMatrix[n_Integer, h_?NumericQ] := Module[{rules},
  rules = Join[
    Table[{i, i} -> -2.0/h^2, {i, 1, n}],
    Table[If[i < n, {i, i + 1} -> 1.0/h^2, Nothing], {i, 1, n}],
    Table[If[i > 1, {i, i - 1} -> 1.0/h^2, Nothing], {i, 1, n}]
  ];
  SparseArray[rules, {n, n}]
];
blockRow[blocks_List] := ArrayFlatten[{blocks}];
weightedGram[m_, wdiag_, sign_] := sign ConjugateTranspose[m].wdiag.m;

buildWallMode[] := Module[
  {
    meshPoints = 129, points, xGrid, radiusValues, weightValues, muEta,
    Tw, TOmega, KEta, Vl, YL0, basisCount = 3, basisValues,
    basisDerivatives, wallMMatrix, wallKMatrix, wallEigenvalues,
    wallEigenvectors, wallPairs, wallCoeffRaw, wallShapeRaw,
    wallWeightNorm, wallCoeff, wallChiEta
  },
  points = N[Table[wMin + wLength i/(meshPoints - 1), {i, 0, meshPoints - 1}], 17];
  xGrid = (points - wMin)/wLength;
  radiusValues = N[r0Geometry /@ points, 17];
  weightValues = N[radiusValues^2, 17];
  muEta = If[physicalMode, ConstantArray[1.0, meshPoints], N[1.0 + 0.08 Cos[Pi xGrid]^2, 17]];
  Tw = If[physicalMode, ConstantArray[1.0, meshPoints], N[1.1 + 0.10 Sin[Pi xGrid], 17]];
  TOmega = If[physicalMode, ConstantArray[1.0, meshPoints], N[0.35 + 0.03 Cos[2 Pi xGrid], 17]];
  KEta = If[physicalMode, ConstantArray[1.0, meshPoints], N[0.45 + 0.08 256 xGrid^4 (1 - xGrid)^4, 17]];
  Vl = N[KEta + Lharm TOmega, 17];
  YL0 = 0.0;
  basisValues = Table[N[xGrid^p, 17], {p, 1, basisCount}];
  basisDerivatives = Table[
    N[If[p == 1, ConstantArray[1.0/wLength, meshPoints], (p/wLength) xGrid^(p - 1)], 17],
    {p, 1, basisCount}
  ];
  wallMMatrix = Table[
    trapz[points, muEta basisValues[[i]] basisValues[[j]]],
    {i, 1, basisCount}, {j, 1, basisCount}
  ];
  wallKMatrix = Table[
    trapz[
      points,
      Tw basisDerivatives[[i]] basisDerivatives[[j]] + Vl basisValues[[i]] basisValues[[j]]
    ] + YL0 basisValues[[i, -1]] basisValues[[j, -1]],
    {i, 1, basisCount}, {j, 1, basisCount}
  ];
  {wallEigenvalues, wallEigenvectors} = Eigensystem[{N[wallKMatrix, 30], N[wallMMatrix, 30]}];
  wallPairs = SortBy[Transpose[{wallEigenvalues, wallEigenvectors}], Re[First[#]] &];
  wallCoeffRaw = Re[wallPairs[[1, 2]]];
  wallShapeRaw = Total[Table[wallCoeffRaw[[i]] basisValues[[i]], {i, 1, basisCount}]];
  If[wallShapeRaw[[-1]] < 0, wallCoeffRaw = -wallCoeffRaw; wallShapeRaw = -wallShapeRaw];
  wallWeightNorm = weightedNorm[points, weightValues, wallShapeRaw];
  wallCoeff = N[wallCoeffRaw/Sqrt[wallWeightNorm], 17];
  wallChiEta = Total[Table[wallCoeff[[i]] basisValues[[i]], {i, 1, basisCount}]];
  <|
    "points" -> points,
    "packet_points" -> N[points - wMin, 17],
    "radius" -> radiusValues,
    "weight" -> weightValues,
    "chi" -> N[wallChiEta, 17],
    "interp" -> Interpolation[Transpose[{points, N[wallChiEta, 17]}], InterpolationOrder -> 1],
    "K" -> N[wallCoeff . wallKMatrix . wallCoeff, 17],
    "M" -> N[wallCoeff . wallMMatrix . wallCoeff, 17],
    "MMatrix" -> wallMMatrix,
    "KMatrix" -> wallKMatrix,
    "coeff" -> wallCoeff,
    "norm" -> N[weightedNorm[points, weightValues, wallChiEta], 17],
    "MMinEigenvalue" -> Min[Eigenvalues[N[wallMMatrix, 30]]]
  |>
];

wallMode = buildWallMode[];

cellIndex[i_, j_, nr_Integer] := (j - 1) nr + i;
rotateByMaxComponent[z_List] := Module[{k, phase},
  k = First[Ordering[Abs[z], -1]];
  phase = If[Abs[z[[k]]] > 0, Exp[-I Arg[z[[k]]]], 1.0];
  phase z
];

runBdgWallSector[useA0_: True] := Module[
  {
    nr = 12, nw = 15, nCells, dr, dw, rFaces, wFaces, rCenters,
    wCenters, radialVolumes, radialFaceAreas, wQuadWeights, laplacianMatrix,
    row, vol, coeff, bdgRadii, bdgW, bdgR0, vconf, projectedRho,
    rho0, psi0, psiR0, psiI0, A0Values, ArValues, AwValues, enthalpy,
    enthalpyPrime, enthalpyLinearization, angularMatrix,
    covariantLaplacianMatrix, kineticMatrix, singleParticleMatrix,
    bdgDiagonalMatrix, bdgPairMatrix, bdgMatrix, stationaryResidualNorm,
    bdgEigenvalues, bdgEigenvectors, bdgPairs, positiveBdgPairs,
    modesToExport, cellVolumeWeights, driveCoeff, modeRecords = {},
    bdgPacketModes = {}, couplingList = {}, varpiList = {},
    bdgResidualList = {}, modeIdx, eigenvalue, rotatedVec, u0, v0,
    symplecticNormBefore, scale, modeVec, u, v, symplecticNormAfter,
    uv, uvImagFraction, uvReal, residualAbs, residualRel, densityByW,
    densityOnPacketGrid, rawPhi, phiNorm, adapterOverlap,
    signedPhysicalCoupling, cAlpha, lambdaB, varpi, B0, B2, B4,
    stieltjesGap, spatialGaugeNorm
  },
  nCells = nr nw;
  dr = rDomainMax/nr;
  dw = wLength/nw;
  rFaces = N[Table[dr i, {i, 0, nr}], 17];
  wFaces = N[Table[wMin + dw j, {j, 0, nw}], 17];
  rCenters = N[Table[(i - 0.5) dr, {i, 1, nr}], 17];
  wCenters = N[Table[wMin + (j - 0.5) dw, {j, 1, nw}], 17];
  radialVolumes = N[Table[(rFaces[[i + 1]]^3 - rFaces[[i]]^3)/3.0, {i, 1, nr}], 17];
  radialFaceAreas = N[rFaces^2, 17];
  wQuadWeights = ConstantArray[dw, nw];
  laplacianMatrix = ConstantArray[0.0, {nCells, nCells}];
  Do[
    row = cellIndex[i, j, nr];
    vol = radialVolumes[[i]];
    If[i > 1,
      coeff = radialFaceAreas[[i]]/(dr vol);
      laplacianMatrix[[row, cellIndex[i - 1, j, nr]]] += coeff;
      laplacianMatrix[[row, row]] -= coeff;
    ];
    If[i < nr,
      coeff = radialFaceAreas[[i + 1]]/(dr vol);
      laplacianMatrix[[row, cellIndex[i + 1, j, nr]]] += coeff;
      laplacianMatrix[[row, row]] -= coeff;,
      coeff = radialFaceAreas[[i + 1]]/(dr vol);
      laplacianMatrix[[row, row]] -= 2.0 coeff;
    ];
    If[j == 1,
      laplacianMatrix[[row, cellIndex[i, j + 1, nr]]] += 1.0/dw^2;
      laplacianMatrix[[row, row]] -= 1.0/dw^2;,
      If[j == nw,
        laplacianMatrix[[row, cellIndex[i, j - 1, nr]]] += 1.0/dw^2;
        laplacianMatrix[[row, row]] -= 1.0/dw^2;,
        laplacianMatrix[[row, cellIndex[i, j - 1, nr]]] += 1.0/dw^2;
        laplacianMatrix[[row, cellIndex[i, j + 1, nr]]] += 1.0/dw^2;
        laplacianMatrix[[row, row]] -= 2.0/dw^2;
      ];
    ];
    ,
    {j, 1, nw}, {i, 1, nr}
  ];
  bdgRadii = Flatten[Table[rCenters[[i]], {j, 1, nw}, {i, 1, nr}]];
  bdgW = Flatten[Table[wCenters[[j]], {j, 1, nw}, {i, 1, nr}]];
  bdgR0 = N[r0Geometry /@ bdgW, 17];
  vconf = N[
    radialWallStrength (bdgRadii/bdgR0)^4 +
      0.5 axialTrapStrength ((bdgW - 0.5 (wMin + wMax))/wLength)^2,
    17
  ];
  projectedRho = projectScalarCells[rhoMatrix, rFaces, wFaces];
  rho0 = projectedRho["values"];
  psi0 = N[Sqrt[Map[Max[#, 0.0] &, rho0]], 17];
  psiR0 = pointVector[psiRMatrix, bdgRadii, bdgW];
  psiI0 = pointVector[psiIMatrix, bdgRadii, bdgW];
  A0Values = If[TrueQ[useA0], pointVector[a0Matrix, bdgRadii, bdgW], ConstantArray[0.0, nCells]];
  ArValues = pointVector[arMatrix, bdgRadii, bdgW];
  AwValues = pointVector[awMatrix, bdgRadii, bdgW];
  enthalpy = N[(5 eosK/4.0) rho0^4, 17];
  enthalpyPrime = N[5 eosK rho0^3, 17];
  enthalpyLinearization = N[rho0 enthalpyPrime, 17];
  angularMatrix = DiagonalMatrix[N[Lharm/bdgRadii^2, 17]];
  covariantLaplacianMatrix = laplacianMatrix - angularMatrix;
  kineticMatrix = N[-(hbar^2/(2 particleMass)) covariantLaplacianMatrix, 17];
  singleParticleMatrix = kineticMatrix + DiagonalMatrix[vconf + gaugeCharge A0Values - chemicalPotential + enthalpy];
  bdgDiagonalMatrix = singleParticleMatrix + DiagonalMatrix[enthalpyLinearization];
  bdgPairMatrix = DiagonalMatrix[enthalpyLinearization];
  bdgMatrix = ArrayFlatten[{{bdgDiagonalMatrix, bdgPairMatrix}, {-bdgPairMatrix, -bdgDiagonalMatrix}}];
  stationaryResidualNorm = N[Norm[singleParticleMatrix . psi0]/Max[Norm[psi0], 10^-30], 17];
  {bdgEigenvalues, bdgEigenvectors} = Eigensystem[N[bdgMatrix, 30]];
  bdgPairs = SortBy[Transpose[{bdgEigenvalues, bdgEigenvectors}], Re[First[#]] &];
  positiveBdgPairs = Select[bdgPairs, Re[#[[1]]] > 10^-8 && Abs[Im[#[[1]]]] < 10^-7 &];
  If[Length[positiveBdgPairs] < 2,
    Print["ERROR: fewer than two positive real BdG eigenmodes found: ", Length[positiveBdgPairs]];
    Exit[6]
  ];
  modesToExport = Min[3, Length[positiveBdgPairs]];
  cellVolumeWeights = N[Flatten[Table[radialVolumes[[i]] wQuadWeights[[j]], {j, 1, nw}, {i, 1, nr}]], 17];
  driveCoeff = N[4.0 radialWallStrength bdgRadii^4/bdgR0^5, 17];
  spatialGaugeNorm = N[Sqrt[Total[(ArValues^2 + AwValues^2) cellVolumeWeights]], 17];
  Do[
    eigenvalue = positiveBdgPairs[[modeIdx, 1]];
    rotatedVec = rotateByMaxComponent[positiveBdgPairs[[modeIdx, 2]]];
    u0 = rotatedVec[[1 ;; nCells]];
    v0 = rotatedVec[[nCells + 1 ;; 2 nCells]];
    symplecticNormBefore = Re[Total[cellVolumeWeights (Abs[u0]^2 - Abs[v0]^2)]];
    scale = 1.0/Sqrt[Max[Abs[symplecticNormBefore], 10^-30]];
    modeVec = scale rotatedVec;
    u = modeVec[[1 ;; nCells]];
    v = modeVec[[nCells + 1 ;; 2 nCells]];
    symplecticNormAfter = Re[Total[cellVolumeWeights (Abs[u]^2 - Abs[v]^2)]];
    uv = u + v;
    uvImagFraction = N[Norm[Im[uv]]/Max[Norm[uv], 10^-30], 17];
    uvReal = Re[uv];
    residualAbs = Norm[bdgMatrix . modeVec - eigenvalue modeVec];
    residualRel = N[residualAbs/Max[Norm[bdgMatrix . modeVec], Abs[eigenvalue] Norm[modeVec], 1.0], 17];
    densityByW = Table[
      Total[Table[
        radialVolumes[[i]] driveCoeff[[cellIndex[i, j, nr]]] psi0[[cellIndex[i, j, nr]]] uvReal[[cellIndex[i, j, nr]]],
        {i, 1, nr}
      ]],
      {j, 1, nw}
    ];
    densityOnPacketGrid = (interp1[wCenters, densityByW, #] &) /@ wallMode["points"];
    rawPhi = N[densityOnPacketGrid/wallMode["weight"], 17];
    phiNorm = normalizeProfile[wallMode["points"], wallMode["weight"], rawPhi, "bdg axial profile"];
    adapterOverlap = weightedInner[wallMode["points"], wallMode["weight"], wallMode["chi"], phiNorm];
    signedPhysicalCoupling = -trapz[wallMode["points"], wallMode["chi"] densityOnPacketGrid];
    If[adapterOverlap signedPhysicalCoupling < 0,
      phiNorm = -phiNorm;
      adapterOverlap = -adapterOverlap;
    ];
    cAlpha = Abs[signedPhysicalCoupling];
    lambdaB = cAlpha/Max[Abs[adapterOverlap], 10^-30];
    varpi = N[Re[eigenvalue]/hbar, 17];
    AppendTo[couplingList, cAlpha];
    AppendTo[varpiList, varpi];
    AppendTo[bdgResidualList, residualRel];
    AppendTo[modeRecords, <|
      "name" -> "bdg_mode_" <> ToString[modeIdx],
      "varpi" -> toJsonNumber[varpi],
      "raw_eigenvalue" -> toJsonNumber[Re[eigenvalue]],
      "eigenvalue_imag_abs" -> toJsonNumber[Abs[Im[eigenvalue]]],
      "relative_eigen_residual" -> toJsonNumber[residualRel],
      "symplectic_norm_before_abs_scaling" -> toJsonNumber[symplecticNormBefore],
      "symplectic_norm_after_abs_scaling" -> toJsonNumber[symplecticNormAfter],
      "uv_imag_fraction_after_phase_rotation" -> toJsonNumber[uvImagFraction],
      "wall_drive_overlap_signed" -> toJsonNumber[signedPhysicalCoupling],
      "c_alpha" -> toJsonNumber[cAlpha],
      "adapter_I_eta_phi" -> toJsonNumber[adapterOverlap],
      "lambda_B_convention_factor" -> toJsonNumber[lambdaB],
      "profile_norm_weighted" -> toJsonNumber[weightedNorm[wallMode["points"], wallMode["weight"], phiNorm]]
    |>];
    AppendTo[bdgPacketModes, <|
      "name" -> "m1c_prep_imported_background_bdg_mode_" <> ToString[modeIdx],
      "lambda_B" -> toJsonNumber[lambdaB],
      "varpi" -> toJsonNumber[varpi],
      "profile_values" -> toJsonVector[phiNorm],
      "derivation_status" -> "derived_eigensolve_on_imported_torch_background",
      "lambda_B_status" -> "normalization_convention_factor_c_alpha_equals_lambda_B_times_I_eta_phi"
    |>];
    ,
    {modeIdx, 1, modesToExport}
  ];
  B0 = N[Total[Table[couplingList[[i]]^2/varpiList[[i]]^2, {i, 1, Length[couplingList]}]], 17];
  B2 = N[Total[Table[couplingList[[i]]^2/varpiList[[i]]^4, {i, 1, Length[couplingList]}]], 17];
  B4 = N[Total[Table[couplingList[[i]]^2/varpiList[[i]]^6, {i, 1, Length[couplingList]}]], 17];
  stieltjesGap = N[B0 B4 - B2^2, 17];
  <|
    "grid" -> {nr, nw},
    "A0_active" -> TrueQ[useA0],
    "projected_rho_mass" -> projectedRho["mass"],
    "projected_rho_mass_relative_error" -> Abs[projectedRho["mass"] - sourceRhoMass]/Max[Abs[sourceRhoMass], 10^-30],
    "stationary_residual_norm_on_projected_operator" -> stationaryResidualNorm,
    "spatial_gauge_norm" -> spatialGaugeNorm,
    "positive_real_modes_found" -> Length[positiveBdgPairs],
    "modes_exported" -> modesToExport,
    "max_relative_eigen_residual" -> Max[bdgResidualList],
    "min_exported_varpi" -> Min[varpiList],
    "modes" -> modeRecords,
    "packet_modes" -> bdgPacketModes,
    "B0" -> B0, "B2" -> B2, "B4" -> B4,
    "stieltjes_gap" -> stieltjesGap,
    "stieltjes_pass" -> TrueQ[stieltjesGap >= -10^-12]
  |>
];

Print["MT15-06 M1c-prep cross-engine imported background"];
Print["background_path: ", backgroundPath];
Print["torch background hash: ", background["content_hash"]];
Print["torch residual linf: ", formatScientific[background["residuals"]["coupled_stationary_linf"]]];
Print["source A0/Ar/Aw linf: ", StringRiffle[formatScientific /@ {sourceA0Linf, sourceArLinf, sourceAwLinf}, ", "]];

bdgActive = runBdgWallSector[True];
bdgZeroA0 = runBdgWallSector[False];

boundaryMatrix[grid_Association, fm_Association] := Module[
  {
    nr = grid["nr"], nw = grid["nw"], n = grid["n"],
    wExitScalar, wExitAngular, rOuterScalar, rOuterAngular,
    cell, wExitWeight, rExitWeight
  },
  wExitScalar = ConstantArray[0.0, n];
  wExitAngular = ConstantArray[0.0, n];
  rOuterScalar = ConstantArray[0.0, n];
  rOuterAngular = ConstantArray[0.0, n];
  Do[
    cell = (j - 1) nr + i;
    If[j == nw,
      wExitWeight = grid["dr"] grid["r"][[cell]]^2 grid["z"][[cell]];
      wExitScalar[[cell]] += wExitWeight;
      wExitAngular[[cell]] += Lharm wExitWeight;
    ];
    If[i == nr,
      rExitWeight = grid["dw"] grid["r"][[cell]]^2 grid["z"][[cell]];
      rOuterScalar[[cell]] += rExitWeight;
      rOuterAngular[[cell]] += Lharm rExitWeight;
    ];
    ,
    {j, 1, nw}, {i, 1, nr}
  ];
  ConjugateTranspose[fm["E_r"]] . diagSparse[wExitScalar] . fm["E_r"] +
    ConjugateTranspose[fm["E_E"]] . diagSparse[wExitAngular + rOuterAngular] . fm["E_E"] +
    ConjugateTranspose[fm["E_B"]] . diagSparse[wExitAngular + rOuterAngular] . fm["E_B"] +
    ConjugateTranspose[fm["E_w"]] . diagSparse[rOuterScalar] . fm["E_w"]
];

rMinTransfer = Max[0.08 rDomainMax, rFacesSrc[[2]]];
rOuterTransfer = rDomainMax;
mouthRadiusForPort = If[physicalMode, aConst, r0Geometry[wMin]];
omegaQ = 3.0 csConst/(2.0 mouthRadiusForPort);
gammaPort = 9.0/(8.0 omegaQ^5);
gammaPortLegacyV209 = mouthRadiusForPort^5/(27.0 csConst^5);
gammaPortCompact = 4.0 mouthRadiusForPort^5/(27.0 csConst^5);
gammaPortFactorVsLegacy = gammaPort/gammaPortLegacyV209;
yOutBranchDescription =
  "Frozen target-blind M1c-prep D2 branch: Y_out,l=2(omega)=i Gamma_port omega^5 on physical tangential lanes; Gamma_port=sigma_Q^can=4 a^5/(27 c_s^5).";

gridData[nr_Integer, nw_Integer] := Module[
  {dr, dw, rVals1, wVals1, rVals, wVals, zVals, r0Vals, quadWeights},
  dr = (rOuterTransfer - rMinTransfer)/(nr + 1);
  dw = wLength/(nw + 1);
  rVals1 = N[Table[rMinTransfer + i dr, {i, 1, nr}], 17];
  wVals1 = N[Table[wMin + i dw, {i, 1, nw}], 17];
  rVals = N[Flatten[Table[rVals1[[i]], {j, 1, nw}, {i, 1, nr}]], 17];
  wVals = N[Flatten[Table[wVals1[[j]], {j, 1, nw}, {i, 1, nr}]], 17];
  zVals = N[zWeight /@ wVals, 17];
  r0Vals = N[r0Geometry /@ wVals, 17];
  quadWeights = N[dr dw rVals^2 zVals, 17];
  <|
    "nr" -> nr, "nw" -> nw, "n" -> nr nw,
    "dr" -> dr, "dw" -> dw,
    "r1" -> rVals1, "w1" -> wVals1,
    "r" -> rVals, "w" -> wVals,
    "z" -> zVals, "R0" -> r0Vals,
    "w0" -> quadWeights,
    "h" -> Max[dr, dw]
  |>
];

buildMaxwellSystem[nr_Integer, nw_Integer, omega_?NumericQ] := Module[
  {
    grid, n, total, drMat1, dwMat1, drrMat1, dwwMat1, idR, idW, idN,
    zN, drMat, dwMat, drrMat, dwwMat, rVals, invR, invR2, rDiag,
    r2Diag, invRDiag, invR2Diag, w0Diag, lw0Diag, wLaneVals, wLane,
    wLaneInv, divR, mG, mEr, mEE, mEB, mEw, mCr, mCE, mCB, mBB,
    mBr, mBE, kAction, kFieldPositive, boundaryB, anchorVals,
    anchorMatrix, kCons, kRet, offsets, anchorCell
  },
  grid = gridData[nr, nw];
  n = grid["n"];
  total = 5 n;
  drMat1 = firstDerivativeMatrix[nr, grid["dr"]];
  dwMat1 = firstDerivativeMatrix[nw, grid["dw"]];
  drrMat1 = secondDerivativeMatrix[nr, grid["dr"]];
  dwwMat1 = secondDerivativeMatrix[nw, grid["dw"]];
  idR = SparseArray[Band[{1, 1}] -> 1.0, {nr, nr}];
  idW = SparseArray[Band[{1, 1}] -> 1.0, {nw, nw}];
  idN = SparseArray[Band[{1, 1}] -> 1.0, {n, n}];
  zN = SparseArray[{}, {n, n}];
  drMat = KroneckerProduct[idW, drMat1];
  dwMat = KroneckerProduct[dwMat1, idR];
  drrMat = KroneckerProduct[idW, drrMat1];
  dwwMat = KroneckerProduct[dwwMat1, idR];
  rVals = grid["r"];
  invR = N[1.0/rVals, 17];
  invR2 = N[invR^2, 17];
  rDiag = diagSparse[rVals];
  r2Diag = diagSparse[rVals^2];
  invRDiag = diagSparse[invR];
  invR2Diag = diagSparse[invR2];
  w0Diag = diagSparse[grid["w0"]];
  lw0Diag = Lharm w0Diag;
  wLaneVals = Join[grid["w0"], grid["w0"], Lharm grid["w0"], Lharm grid["w0"], grid["w0"]];
  wLane = diagSparse[wLaneVals];
  wLaneInv = diagSparse[1.0/wLaneVals];
  divR = invR2Diag . drMat . r2Diag;
  mG = blockRow[{I omega idN, divR, -Lharm invRDiag, zN, dwMat}];
  mEr = blockRow[{-drMat, I omega idN, zN, zN, zN}];
  mEE = blockRow[{-invRDiag, zN, I omega idN, zN, zN}];
  mEB = blockRow[{zN, zN, zN, I omega idN, zN}];
  mEw = blockRow[{-dwMat, zN, zN, zN, I omega idN}];
  mCr = blockRow[{zN, -dwMat, zN, zN, drMat}];
  mCE = blockRow[{zN, zN, -dwMat, zN, invRDiag}];
  mCB = blockRow[{zN, zN, zN, -dwMat, zN}];
  mBB = blockRow[{zN, -invRDiag, invRDiag . drMat . rDiag, zN, zN}];
  mBr = blockRow[{zN, zN, zN, -Lharm invRDiag, zN}];
  mBE = blockRow[{zN, zN, zN, -invRDiag . drMat . rDiag, zN}];
  kAction =
    weightedGram[mEr, w0Diag, 1.0] + weightedGram[mEE, lw0Diag, 1.0] +
    weightedGram[mEB, lw0Diag, 1.0] + weightedGram[mEw, w0Diag, 1.0] -
    weightedGram[mCr, w0Diag, 1.0] - weightedGram[mCE, lw0Diag, 1.0] -
    weightedGram[mCB, lw0Diag, 1.0] - weightedGram[mBB, lw0Diag, 1.0] -
    weightedGram[mBr, w0Diag, 1.0] - weightedGram[mBE, lw0Diag, 1.0] -
    (1.0/xiGauge) weightedGram[mG, w0Diag, 1.0];
  kFieldPositive =
    weightedGram[mEr, w0Diag, 1.0] + weightedGram[mEE, lw0Diag, 1.0] +
    weightedGram[mEB, lw0Diag, 1.0] + weightedGram[mEw, w0Diag, 1.0] +
    weightedGram[mCr, w0Diag, 1.0] + weightedGram[mCE, lw0Diag, 1.0] +
    weightedGram[mCB, lw0Diag, 1.0] + weightedGram[mBB, lw0Diag, 1.0] +
    weightedGram[mBr, w0Diag, 1.0] + weightedGram[mBE, lw0Diag, 1.0];
  boundaryB = boundaryMatrix[grid, <|"E_r" -> mEr, "E_E" -> mEE, "E_B" -> mEB, "E_w" -> mEw|>];
  anchorVals = ConstantArray[0.0, total];
  anchorCell = Ceiling[n/2];
  anchorVals[[anchorCell]] = 10^-8 grid["w0"][[anchorCell]];
  anchorMatrix = diagSparse[anchorVals];
  kCons = kAction + anchorMatrix;
  kRet = kCons + I gammaPort omega^5 boundaryB;
  offsets = <|"phi" -> 0, "a_r" -> n, "a_E" -> 2 n, "a_B" -> 3 n, "a_w" -> 4 n|>;
  Join[
    grid,
    <|
      "omega" -> omega, "total" -> total,
      "Dr" -> drMat, "Dw" -> dwMat, "Drr" -> drrMat, "Dww" -> dwwMat,
      "W0" -> w0Diag, "WLane" -> wLane, "WLaneInv" -> wLaneInv,
      "WLaneValues" -> wLaneVals, "KAction" -> kAction, "KCons" -> kCons,
      "KRet" -> kRet, "ACons" -> wLaneInv . kCons, "ARet" -> wLaneInv . kRet,
      "BoundaryB" -> boundaryB, "BoundaryLossStrong" -> wLaneInv . boundaryB,
      "AnchorMatrix" -> anchorMatrix, "ScalarGaugeAnchor" -> <|"cell" -> anchorCell, "strength" -> 10^-8|>,
      "KFieldPositive" -> kFieldPositive,
      "FieldMatrices" -> <|
        "G" -> mG, "E_r" -> mEr, "E_E" -> mEE, "E_B" -> mEB, "E_w" -> mEw,
        "C_r" -> mCr, "C_E" -> mCE, "C_B" -> mCB,
        "B_B" -> mBB, "B_r" -> mBr, "B_E" -> mBE
      |>,
      "Offsets" -> offsets
    |>
  ]
];

linearizedCurrentFormula[sys_Association, psiR0_List, psiI0_List, dpsiR_List, dpsiI_List, ar0_List, aw0_List, dar_List, daw_List] := Module[
  {gradPsiRr, gradPsiIr, gradPsiRw, gradPsiIw, gradDPsiRr, gradDPsiIr, gradDPsiRw, gradDPsiIw, density0, deltaDensity, phaseR, phaseW, gaugeR, gaugeW, jr, jw},
  gradPsiRr = sys["Dr"] . psiR0;
  gradPsiIr = sys["Dr"] . psiI0;
  gradPsiRw = sys["Dw"] . psiR0;
  gradPsiIw = sys["Dw"] . psiI0;
  gradDPsiRr = sys["Dr"] . dpsiR;
  gradDPsiIr = sys["Dr"] . dpsiI;
  gradDPsiRw = sys["Dw"] . dpsiR;
  gradDPsiIw = sys["Dw"] . dpsiI;
  density0 = psiR0^2 + psiI0^2;
  deltaDensity = 2.0 (psiR0 dpsiR + psiI0 dpsiI);
  phaseR = dpsiR gradPsiIr + psiR0 gradDPsiIr - dpsiI gradPsiRr - psiI0 gradDPsiRr;
  phaseW = dpsiR gradPsiIw + psiR0 gradDPsiIw - dpsiI gradPsiRw - psiI0 gradDPsiRw;
  gaugeR = dar density0 + ar0 deltaDensity;
  gaugeW = daw density0 + aw0 deltaDensity;
  jr = (hbar/particleMass) phaseR - (gaugeCharge/particleMass) gaugeR;
  jw = (hbar/particleMass) phaseW - (gaugeCharge/particleMass) gaugeW;
  <|"deltaRho" -> deltaDensity, "jr" -> jr, "jw" -> jw|>
];

nonlinearNumberCurrent[sys_Association, psiR_List, psiI_List, ar_List, aw_List] := Module[{density, jr, jw},
  density = psiR^2 + psiI^2;
  jr = (hbar/particleMass) (psiR (sys["Dr"] . psiI) - psiI (sys["Dr"] . psiR)) - (gaugeCharge/particleMass) ar density;
  jw = (hbar/particleMass) (psiR (sys["Dw"] . psiI) - psiI (sys["Dw"] . psiR)) - (gaugeCharge/particleMass) aw density;
  <|"jr" -> jr, "jw" -> jw|>
];

buildMatterSource[sys_Association, omega_?NumericQ, spatialSign_: 1.0] := Module[
  {
    n = sys["n"], idN, lap, bdgRadii, bdgW, bdgR0, vconf, rFacesT, wFacesT,
    projectedRho, rho0, psiR0, psiI0, psi0, A0Values, ArValues, AwValues,
    enthalpy, enthalpyPrime, enthalpyLin, kinetic, single, lPlus,
    lPhase, etaVals, deltaV, drive, dpsiR, dpsiI, zero, current,
    deltaRho, jr, jw, jE, jVec, ampResidual, phaseResidual,
    sourceNorm, continuityResidual
  },
  idN = SparseArray[Band[{1, 1}] -> 1.0, {n, n}];
  lap = sys["Drr"] + diagSparse[2.0/sys["r"]] . sys["Dr"] + sys["Dww"] - diagSparse[Lharm/sys["r"]^2];
  bdgRadii = sys["r"];
  bdgW = sys["w"];
  bdgR0 = sys["R0"];
  vconf = N[
    radialWallStrength (bdgRadii/bdgR0)^4 +
      0.5 axialTrapStrength ((bdgW - 0.5 (wMin + wMax))/wLength)^2,
    17
  ];
  rFacesT = N[Join[
      {rMinTransfer},
      Table[0.5 (sys["r1"][[i]] + sys["r1"][[i + 1]]), {i, 1, sys["nr"] - 1}],
      {rOuterTransfer}
    ], 17];
  wFacesT = N[Join[
      {wMin},
      Table[0.5 (sys["w1"][[i]] + sys["w1"][[i + 1]]), {i, 1, sys["nw"] - 1}],
      {wMax}
    ], 17];
  projectedRho = projectScalarCells[rhoMatrix, rFacesT, wFacesT];
  rho0 = projectedRho["values"];
  psiR0 = pointVector[psiRMatrix, bdgRadii, bdgW];
  psiI0 = pointVector[psiIMatrix, bdgRadii, bdgW];
  psi0 = N[Sqrt[Map[Max[#, 0.0] &, rho0]], 17];
  A0Values = pointVector[a0Matrix, bdgRadii, bdgW];
  ArValues = spatialSign pointVector[arMatrix, bdgRadii, bdgW];
  AwValues = spatialSign pointVector[awMatrix, bdgRadii, bdgW];
  enthalpy = N[(5 eosK/4.0) rho0^4, 17];
  enthalpyPrime = N[5 eosK rho0^3, 17];
  enthalpyLin = N[rho0 enthalpyPrime, 17];
  kinetic = N[-(hbar^2/(2 particleMass)) lap, 17];
  single = kinetic + diagSparse[vconf + gaugeCharge A0Values - chemicalPotential + enthalpy];
  lPlus = single + diagSparse[2.0 enthalpyLin] + (omega^2) idN;
  lPhase = single + phaseRegularizer idN + (omega^2) idN;
  etaVals = N[wallMode["interp"] /@ bdgW, 17];
  deltaV = N[-4.0 radialWallStrength bdgRadii^4/bdgR0^5 etaVals, 17];
  drive = N[-deltaV psi0, 17];
  dpsiR = LinearSolve[lPlus, drive];
  dpsiI = I omega LinearSolve[lPhase, dpsiR];
  zero = ConstantArray[0.0, n];
  current = linearizedCurrentFormula[sys, psiR0, psiI0, dpsiR, dpsiI, ArValues, AwValues, zero, zero];
  deltaRho = current["deltaRho"];
  jr = current["jr"];
  jw = current["jw"];
  jE = (hbar/particleMass) psi0 dpsiI/sys["r"];
  jVec = Join[gaugeCharge deltaRho, gaugeCharge jr, gaugeCharge jE, zero, gaugeCharge jw];
  ampResidual = Norm[lPlus . dpsiR - drive]/Max[Norm[drive], 10^-30];
  phaseResidual = Norm[lPhase . (dpsiI/(I Max[omega, 10^-30])) - dpsiR]/Max[Norm[dpsiR], 10^-30];
  sourceNorm = Sqrt[Re[Conjugate[jVec].(sys["WLane"].jVec)]];
  continuityResidual = Norm[-I omega deltaRho + sys["Dr"] . jr + sys["Dw"] . jw - Lharm jE/sys["r"]]/
    Max[Norm[deltaRho] + Norm[jr] + Norm[jw] + Norm[jE], 10^-30];
  <|
    "jVec" -> jVec, "deltaRho" -> deltaRho, "jr" -> jr, "jE" -> jE, "jB" -> zero, "jw" -> jw,
    "dpsiR" -> dpsiR, "dpsiI" -> dpsiI, "psi0" -> psi0, "deltaV" -> deltaV,
    "ampResidual" -> N[ampResidual, 17], "phaseResidual" -> N[phaseResidual, 17],
    "sourceNorm" -> N[sourceNorm, 17], "continuityResidual" -> N[continuityResidual, 17],
    "projected_rho_mass" -> projectedRho["mass"], "usedReturnSourceSEtaA" -> False
  |>
];

currentFrechetCheck[] := Module[
  {sys, n, rr, ww, psiR0, psiI0, ar0, aw0, dpsiR, dpsiI, dar, daw, eps = 10^-6, plus, minus, fdR, fdW, lin, errR, errW, normR, normW},
  sys = buildMaxwellSystem[5, 4, 0.0];
  n = sys["n"]; rr = sys["r"]; ww = sys["w"];
  psiR0 = pointVector[psiRMatrix, rr, ww] + 0.03 rr;
  psiI0 = pointVector[psiIMatrix, rr, ww] + 0.02 ww;
  ar0 = pointVector[arMatrix, rr, ww] + 0.04 rr + 0.01 ww;
  aw0 = pointVector[awMatrix, rr, ww] - 0.03 rr + 0.02 ww;
  dpsiR = 0.05 rr^2 - 0.02 ww;
  dpsiI = -0.03 rr + 0.04 ww^2;
  dar = 0.02 rr (1.0 + ww);
  daw = -0.015 (1.0 + rr) ww;
  plus = nonlinearNumberCurrent[sys, psiR0 + eps dpsiR, psiI0 + eps dpsiI, ar0 + eps dar, aw0 + eps daw];
  minus = nonlinearNumberCurrent[sys, psiR0 - eps dpsiR, psiI0 - eps dpsiI, ar0 - eps dar, aw0 - eps daw];
  fdR = (plus["jr"] - minus["jr"])/(2.0 eps);
  fdW = (plus["jw"] - minus["jw"])/(2.0 eps);
  lin = linearizedCurrentFormula[sys, psiR0, psiI0, dpsiR, dpsiI, ar0, aw0, dar, daw];
  errR = Norm[lin["jr"] - fdR]; errW = Norm[lin["jw"] - fdW];
  normR = Max[Norm[fdR], 10^-30]; normW = Max[Norm[fdW], 10^-30];
  <|"grid" -> {sys["nr"], sys["nw"]}, "relative_error_jr" -> N[errR/normR, 17], "relative_error_jw" -> N[errW/normW, 17], "max_relative_error" -> N[Max[errR/normR, errW/normW], 17], "matches_step8c_frechet_formula" -> TrueQ[Max[errR/normR, errW/normW] < 10^-7]|>
];

selfEnergySolve[sys_Association, jVec_List, mode_String] := Module[{op, solution, residual, sigma, boundaryQuadratic, flux},
  op = If[mode == "retarded", sys["ARet"], sys["ACons"]];
  solution = LinearSolve[op, mu0 jVec];
  residual = Norm[op . solution - mu0 jVec]/Max[Norm[mu0 jVec], 10^-30];
  sigma = Conjugate[jVec].(sys["WLane"].solution);
  boundaryQuadratic = Re[Conjugate[solution].(sys["BoundaryB"].solution)];
  flux = gammaPort sys["omega"]^5 boundaryQuadratic/(2.0 mu0);
  <|"solution" -> solution, "sigma" -> N[sigma, 17], "residual" -> N[residual, 17], "boundaryQuadratic" -> N[boundaryQuadratic, 17], "outgoingFlux" -> N[flux, 17]|>
];

fitEvenCoefficients[xs_List, ys_List] := Module[{mat, coeff},
  mat = N[Transpose[{ConstantArray[1.0, Length[xs]], xs, xs^2}], 30];
  coeff = LeastSquares[mat, N[ys, 30]];
  N[coeff, 17]
];
estimateOrderTriplet[hs_List, vals_List] := Module[{d12, d23},
  d12 = Abs[vals[[1]] - vals[[2]]]; d23 = Abs[vals[[2]] - vals[[3]]];
  If[d12 <= 0 || d23 <= 0, Indeterminate, N[Log[d12/d23]/Log[hs[[1]]/hs[[2]]], 17]]
];

fieldNormFromComponents[sys_Association, vec_List] := Module[{fm = sys["FieldMatrices"], w0, lw0, terms, compVals},
  w0 = sys["W0"]; lw0 = Lharm w0;
  terms = {{fm["E_r"], w0}, {fm["E_E"], lw0}, {fm["E_B"], lw0}, {fm["E_w"], w0}, {fm["C_r"], w0}, {fm["C_E"], lw0}, {fm["C_B"], lw0}, {fm["B_B"], lw0}, {fm["B_r"], w0}, {fm["B_E"], lw0}};
  compVals = Table[
    With[{component = terms[[i, 1]] . vec, weight = terms[[i, 2]]}, Conjugate[component] . (weight . component)],
    {i, 1, Length[terms]}
  ];
  Sqrt[Max[0.0, Re[Total[compVals]]]]
];
makeGaugeVector[sys_Association, chi_List] := Module[{omega = sys["omega"], drMat = sys["Dr"], dwMat = sys["Dw"], rVals = sys["r"], n = sys["n"]},
  Join[I omega chi, drMat . chi, chi/rVals, ConstantArray[0.0, n], dwMat . chi]
];
pureGaugeCheck[] := Module[{sys, chi, gaugeVec, source, solution, fieldNorm, flux, sourceNorm},
  sys = buildMaxwellSystem[6, 5, 0.043];
  chi = N[MapThread[Sin[Pi ((#1 - rMinTransfer)/(rOuterTransfer - rMinTransfer))]^2 Sin[Pi (#2 - wMin)/wLength]^2 &, {sys["r"], sys["w"]}], 17];
  gaugeVec = makeGaugeVector[sys, chi];
  source = sys["ARet"] . gaugeVec;
  solution = LinearSolve[sys["ARet"], source];
  fieldNorm = fieldNormFromComponents[sys, solution];
  flux = gammaPort sys["omega"]^5 Re[Conjugate[solution].(sys["BoundaryB"].solution)]/(2.0 mu0);
  sourceNorm = Sqrt[Re[Conjugate[source].(sys["WLane"].source)]];
  <|"grid" -> {sys["nr"], sys["nw"]}, "omega" -> sys["omega"], "source_norm" -> N[sourceNorm, 17], "physical_field_norm" -> N[fieldNorm, 17], "outgoing_flux" -> N[flux, 17], "physical_transfer_zero" -> TrueQ[fieldNorm < 10^-10 && Abs[flux] < 10^-18]|>
];

basisInvarianceCheck[] := Module[
  {sys, matter, jVec, omega, cons, sigmaBase, wSqrt, wInvSqrt, aTildeCons, sTilde, bTilde, dim, deltas, seed = 240617, qRows, q, sigRot, solRot, nRot, k, uBase, nBase},
  omega = 0.039; sys = buildMaxwellSystem[4, 4, omega]; matter = buildMatterSource[sys, omega]; jVec = matter["jVec"];
  cons = selfEnergySolve[sys, jVec, "conservative"]; sigmaBase = cons["sigma"];
  wSqrt = diagSparse[Sqrt[sys["WLaneValues"]]]; wInvSqrt = diagSparse[1.0/Sqrt[sys["WLaneValues"]]];
  aTildeCons = Normal[wSqrt . sys["ACons"] . wInvSqrt]; bTilde = Normal[wInvSqrt . sys["BoundaryB"] . wInvSqrt]; sTilde = Sqrt[sys["WLaneValues"]] jVec;
  uBase = LinearSolve[aTildeCons, sTilde]; nBase = Re[Conjugate[uBase] . bTilde . uBase];
  dim = Length[sTilde]; SeedRandom[seed];
  deltas = Table[
    qRows = Orthogonalize[RandomReal[NormalDistribution[0, 1], {dim, dim}]];
    q = N[qRows, 17];
    solRot = LinearSolve[q . aTildeCons . Transpose[q], q . sTilde];
    sigRot = Conjugate[q . sTilde] . solRot;
    nRot = Re[Conjugate[solRot] . (q . bTilde . Transpose[q]) . solRot];
    <|"rotation" -> k, "relative_Z_delta" -> N[Abs[sigRot - sigmaBase]/Max[Abs[sigmaBase], 10^-30], 17], "relative_N_delta" -> N[Abs[nRot - nBase]/Max[Abs[nBase], 10^-30], 17]|>,
    {k, 1, 3}
  ];
  <|"grid" -> {sys["nr"], sys["nw"]}, "omega" -> omega, "seed" -> seed, "base_Z" -> Re[sigmaBase], "base_N" -> N[nBase, 17], "rotations" -> deltas, "max_relative_Z_delta" -> Max[#["relative_Z_delta"] & /@ deltas], "max_relative_N_delta" -> Max[#["relative_N_delta"] & /@ deltas], "pass" -> TrueQ[Max[#["relative_Z_delta"] & /@ deltas] < 10^-10 && Max[#["relative_N_delta"] & /@ deltas] < 10^-10]|>
];

v209Regression[] := Module[{omegaU = 3.25, omegaW = 4.35, rMix = 0.08, gU = 0.18, gW = 0.13, aBlock, bBlock, kInt, g, ellVec, inv, delta, qNum, zDirect, zAnalytic, nDirect, nAnalytic},
  aBlock = omegaU^2; bBlock = omegaW^2; kInt = {{aBlock, -rMix}, {-rMix, bBlock}}; g = {gU, gW}; ellVec = {0.0, 1.0};
  inv = Inverse[kInt]; delta = aBlock bBlock - rMix^2; qNum = gU^2 bBlock + 2.0 gU gW rMix + gW^2 aBlock;
  zDirect = g . inv . g; zAnalytic = qNum/delta; nDirect = (ellVec . inv . g)^2; nAnalytic = (aBlock gW + rMix gU)^2/delta^2;
  <|"Omega_U" -> omegaU, "Omega_W" -> omegaW, "R" -> rMix, "g_U" -> gU, "g_W" -> gW, "Delta" -> delta, "Z_direct" -> zDirect, "Z_Q_over_Delta" -> zAnalytic, "N0_direct" -> nDirect, "N0_square_over_Delta_squared" -> nAnalytic, "relative_Z_error" -> Abs[zDirect - zAnalytic]/Max[Abs[zAnalytic], 10^-30], "relative_N0_error" -> Abs[nDirect - nAnalytic]/Max[Abs[nAnalytic], 10^-30], "pass" -> TrueQ[Abs[zDirect - zAnalytic]/Max[Abs[zAnalytic], 10^-30] < 10^-13 && Abs[nDirect - nAnalytic]/Max[Abs[nAnalytic], 10^-30] < 10^-13]|>
];

hardCapContrast[omega_?NumericQ] := Module[{sys, matter, open, hardPenalty = 10^6, hardOp, hardSolution, openNorm, hardNorm, ratio},
  sys = buildMaxwellSystem[7, 7, omega]; matter = buildMatterSource[sys, omega]; open = selfEnergySolve[sys, matter["jVec"], "retarded"];
  hardOp = sys["WLaneInv"] . (sys["KCons"] + hardPenalty sys["BoundaryB"]);
  hardSolution = Quiet[LinearSolve[hardOp, mu0 matter["jVec"]]];
  openNorm = Sqrt[Max[0.0, Re[Conjugate[open["solution"]].(sys["BoundaryB"].open["solution"])]]];
  hardNorm = Sqrt[Max[0.0, Re[Conjugate[hardSolution].(sys["BoundaryB"].hardSolution)]]];
  ratio = hardNorm/Max[openNorm, 10^-30];
  <|"grid" -> {sys["nr"], sys["nw"]}, "omega" -> omega, "open_boundary_norm" -> N[openNorm, 17], "hard_cap_boundary_norm" -> N[hardNorm, 17], "hard_cap_to_open_boundary_norm_ratio" -> N[ratio, 17], "open_not_hard_cap" -> TrueQ[openNorm > 10^-10 && ratio < 10^-4], "outgoing_flux" -> open["outgoingFlux"], "outgoing_flux_positive" -> TrueQ[open["outgoingFlux"] > 0]|>
];

transferForGrid[nr_Integer, nw_Integer, window_?NumericQ, spatialSign_: 1.0] := Module[
  {fractions = {0.55, 0.70, 0.85, 1.00, 1.15}, omegas, rows, omega, sys, matter, cons, ret, xVals, zVals, nVals, zCoeffs, nCoeffs},
  omegas = N[window fractions, 17];
  rows = Table[
    sys = buildMaxwellSystem[nr, nw, omega];
    matter = buildMatterSource[sys, omega, spatialSign];
    cons = selfEnergySolve[sys, matter["jVec"], "conservative"];
    ret = selfEnergySolve[sys, matter["jVec"], "retarded"];
    <|"omega" -> omega, "sigmaCons" -> cons["sigma"], "sigmaRet" -> ret["sigma"], "NofOmega" -> N[-Im[ret["sigma"]]/(gammaPort omega^5), 17], "consResidual" -> cons["residual"], "retResidual" -> ret["residual"], "ampResidual" -> matter["ampResidual"], "phaseResidual" -> matter["phaseResidual"], "continuityResidual" -> matter["continuityResidual"], "outgoingFlux" -> ret["outgoingFlux"], "boundaryQuadratic" -> ret["boundaryQuadratic"], "sourceNorm" -> matter["sourceNorm"]|>,
    {omega, omegas}
  ];
  xVals = omegas^2; zVals = Re[#["sigmaCons"]] & /@ rows; nVals = #["NofOmega"] & /@ rows;
  zCoeffs = fitEvenCoefficients[xVals, zVals]; nCoeffs = fitEvenCoefficients[xVals, nVals];
  <|
    "grid" -> {nr, nw}, "h" -> gridData[nr, nw]["h"], "window" -> window,
    "omegas" -> omegas, "Z0" -> zCoeffs[[1]], "Z2" -> zCoeffs[[2]], "Z4" -> zCoeffs[[3]],
    "N0" -> nCoeffs[[1]], "N2" -> nCoeffs[[2]], "N4" -> nCoeffs[[3]], "rows" -> rows,
    "maxGreenResidual" -> Max[Max[#["consResidual"] & /@ rows], Max[#["retResidual"] & /@ rows]],
    "maxBdgAmpResidual" -> Max[#["ampResidual"] & /@ rows],
    "maxBdgPhaseResidual" -> Max[#["phaseResidual"] & /@ rows],
    "maxContinuityResidual" -> Max[#["continuityResidual"] & /@ rows],
    "outgoingFluxAtMaxOmega" -> Last[rows]["outgoingFlux"],
    "boundaryQuadraticAtMaxOmega" -> Last[rows]["boundaryQuadratic"],
    "sourceNormAtMinOmega" -> First[rows]["sourceNorm"]
  |>
];

extractionForBackgroundInterpolationOrder[order_Integer, transferGrid_List, window_?NumericQ] := Block[
  {r0Geometry, zWeight, pointVector, wallMode, omegaQ, gammaPort, gammaPortLegacyV209,
    gammaPortCompact, gammaPortFactorVsLegacy},
  Module[{bdg, transfer},
    r0Geometry[w_?NumericQ] := If[physicalMode,
      Module[{x = Clip[(w - wMin)/wLength, {0.0, 1.0}]},
        aConst + (rExitConst - aConst) (3.0 x^2 - 2.0 x^3)
      ],
      interp1Order[wCentersSrc, r0ValuesSrc, w, order]
    ];
    zWeight[w_?NumericQ] := interp1Order[wCentersSrc, zValuesSrc, w, order];
    pointVector[m_List, rVals_List, wVals_List] :=
      N[MapThread[interp2Order[m, #1, #2, order] &, {rVals, wVals}], 17];
    wallMode = buildWallMode[];
    omegaQ = 3.0 csConst/(2.0 r0Geometry[wMin]);
    gammaPort = 9.0/(8.0 omegaQ^5);
    gammaPortLegacyV209 = r0Geometry[wMin]^5/(27.0 csConst^5);
    gammaPortCompact = 4.0 r0Geometry[wMin]^5/(27.0 csConst^5);
    gammaPortFactorVsLegacy = gammaPort/gammaPortLegacyV209;
    bdg = runBdgWallSector[True];
    transfer = transferForGrid[transferGrid[[1]], transferGrid[[2]], window];
    <|
      "interpolation_order" -> order,
      "scope" -> "diagnostic-only imported smooth-field/background interpolation pass; primary order-1 extraction is unchanged",
      "transfer_grid" -> transferGrid,
      "coefficients" -> <|
        "K" -> wallMode["K"], "M" -> wallMode["M"],
        "B0" -> bdg["B0"], "B2" -> bdg["B2"], "B4" -> bdg["B4"],
        "Z0" -> transfer["Z0"], "Z2" -> transfer["Z2"], "Z4" -> transfer["Z4"],
        "N0" -> transfer["N0"], "N2" -> transfer["N2"], "N4" -> transfer["N4"]
      |>,
      "diagnostics" -> <|
        "max_bdg_eigen_residual" -> bdg["max_relative_eigen_residual"],
        "max_green_residual" -> transfer["maxGreenResidual"],
        "max_bdg_amp_residual" -> transfer["maxBdgAmpResidual"],
        "max_bdg_phase_residual" -> transfer["maxBdgPhaseResidual"],
        "gamma_port" -> gammaPort
      |>
    |>
  ]
];

coefficientRelativeDeltaTable[primary_Association, comparison_Association, keys_List] := Association @ Table[
  With[{p = primary[key], c = comparison[key], delta = comparison[key] - primary[key]},
    key -> <|
      "order1_primary" -> toJsonNumber[p],
      "order2_diagnostic" -> toJsonNumber[c],
      "absolute_delta_order2_minus_order1" -> toJsonNumber[delta],
      "relative_delta_order2_minus_order1" -> toJsonNumber[delta/Max[Abs[p], 10^-30]],
      "relative_abs_delta" -> toJsonNumber[Abs[delta]/Max[Abs[p], 10^-30]]
    |>
  ],
  {key, keys}
];

Print["BdG max rel eigen residual: ", formatScientific[bdgActive["max_relative_eigen_residual"]]];
Print["BdG B0/B2/B4 active: ", StringRiffle[formatScientific /@ {bdgActive["B0"], bdgActive["B2"], bdgActive["B4"]}, ", "]];

currentCheck = currentFrechetCheck[];
meshWindow = 0.040;
meshGrids = {{7, 7}, {9, 9}, {11, 11}};
meshTransfers = Table[
  Print["transfer mesh grid ", grid[[1]], "x", grid[[2]], " window ", meshWindow];
  transferForGrid[grid[[1]], grid[[2]], meshWindow],
  {grid, meshGrids}
];
primaryTransfer = Last[meshTransfers];
spatialFlipTransfer = transferForGrid[11, 11, meshWindow, -1.0];

meshHs = #["h"] & /@ meshTransfers;
meshZ0s = #["Z0"] & /@ meshTransfers;
meshN0s = #["N0"] & /@ meshTransfers;
meshZOrders = <|"Z0" -> estimateOrderTriplet[meshHs, meshZ0s], "N0" -> estimateOrderTriplet[meshHs, meshN0s]|>;

pureGauge = pureGaugeCheck[];
basisInvariant = basisInvarianceCheck[];
v209 = v209Regression[];
capContrast = hardCapContrast[0.040];

passChecks = <|
  "current_frechet_matches_step8c" -> currentCheck["matches_step8c_frechet_formula"],
  "outgoing_flux_positive" -> capContrast["outgoing_flux_positive"],
  "open_not_hard_cap" -> capContrast["open_not_hard_cap"],
  "pure_gauge_zero_physical_transfer" -> pureGauge["physical_transfer_zero"],
  "basis_invariance" -> basisInvariant["pass"],
  "v2_09_regression" -> v209["pass"],
  "green_residuals_small" -> TrueQ[primaryTransfer["maxGreenResidual"] < 10^-8],
  "bdg_residuals_small" -> TrueQ[primaryTransfer["maxBdgAmpResidual"] < 10^-8 && primaryTransfer["maxBdgPhaseResidual"] < 10^-8 && bdgActive["max_relative_eigen_residual"] < 10^-10],
  "N0_positive" -> TrueQ[primaryTransfer["N0"] > 0]
|>;

A0Changed = <|
  "B0_delta_active_minus_zeroA0" -> toJsonNumber[bdgActive["B0"] - bdgZeroA0["B0"]],
  "B2_delta_active_minus_zeroA0" -> toJsonNumber[bdgActive["B2"] - bdgZeroA0["B2"]],
  "B4_delta_active_minus_zeroA0" -> toJsonNumber[bdgActive["B4"] - bdgZeroA0["B4"]],
  "max_varpi_floor_active" -> toJsonNumber[bdgActive["min_exported_varpi"]],
  "max_varpi_floor_zeroA0" -> toJsonNumber[bdgZeroA0["min_exported_varpi"]]
|>;
spatialGaugeSensitivity = <|
  "Z0_delta_spatial_sign_flip" -> toJsonNumber[spatialFlipTransfer["Z0"] - primaryTransfer["Z0"]],
  "N0_delta_spatial_sign_flip" -> toJsonNumber[spatialFlipTransfer["N0"] - primaryTransfer["N0"]],
  "relative_Z0_delta" -> toJsonNumber[Abs[spatialFlipTransfer["Z0"] - primaryTransfer["Z0"]]/Max[Abs[primaryTransfer["Z0"]], 10^-30]],
  "relative_N0_delta" -> toJsonNumber[Abs[spatialFlipTransfer["N0"] - primaryTransfer["N0"]]/Max[Abs[primaryTransfer["N0"]], 10^-30]]
|>;

Kwall = wallMode["K"]; Mwall = wallMode["M"];
B0 = bdgActive["B0"]; B2 = bdgActive["B2"]; B4 = bdgActive["B4"];
Z0 = primaryTransfer["Z0"]; Z2 = primaryTransfer["Z2"]; Z4 = primaryTransfer["Z4"];
N0 = primaryTransfer["N0"]; N2 = primaryTransfer["N2"]; N4 = primaryTransfer["N4"];
d0 = Kwall - B0 - Z0;
d2 = -(Mwall + B2 + Z2);
d4 = -(B4 + Z4);
p0 = N0/d0;
p2 = (d0 N2 - 2 d2 N0)/d0^2;
p4 = (d0^2 N4 - 2 d0 (d2 N2 + d4 N0) + 3 d2^2 N0)/d0^3;
rPole = d0 (B4 + Z4) - 3 (Mwall + B2 + Z2)^2;
rNormMachinery = mhat0Const^2 SPortConst p0 - 54 GConst csConst^5/(5 aConst^5 cConst^5);

directCoefficients = <|"K" -> Kwall, "M" -> Mwall, "B0" -> B0, "B2" -> B2, "B4" -> B4, "Z0" -> Z0, "Z2" -> Z2, "Z4" -> Z4, "N0" -> N0, "N2" -> N2, "N4" -> N4|>;
interpSensitivityKeys = {"K", "M", "B0", "B2", "B4", "Z0", "Z2", "Z4", "N0", "N2", "N4"};
Print["interp sensitivity diagnostic: order-2 background interpolation on primary grid ", primaryTransfer["grid"]];
order2InterpolationExtraction = extractionForBackgroundInterpolationOrder[2, primaryTransfer["grid"], meshWindow];
interpolationMethodDeltas = coefficientRelativeDeltaTable[
  directCoefficients,
  order2InterpolationExtraction["coefficients"],
  interpSensitivityKeys
];
interpolationMethodCarryForwardNote =
  If[physicalMode,
    "Order2-vs-order1 interpolation-method sensitivity is a diagnostic on the same frozen physical torch background. It is distinct from the WP1 background-resolution term and is carried explicitly into the M1c section J numerical-only budget.",
    "Order2-vs-order1 interpolation-method sensitivity is a diagnostic on the same coarse 10x8 torch background. It is distinct from, and smaller than, the unresolved background-resolution error of that coarse background. The production background-resolution term and full interpolation error belong in the M1c-run section J error budget; this M1c-prep machinery result is not resolution-converged."
  ];
interpolationMethodSensitivity = <|
  "primary_interpolation_order" -> 1,
  "diagnostic_interpolation_order" -> order2InterpolationExtraction["interpolation_order"],
  "scope" -> order2InterpolationExtraction["scope"],
  "same_torch_background_content_hash" -> background["content_hash"],
  "note" -> interpolationMethodCarryForwardNote,
  "order2_coefficients" -> Map[toJsonNumber, order2InterpolationExtraction["coefficients"]],
  "relative_deltas_order2_vs_order1" -> interpolationMethodDeltas,
  "diagnostics" -> Map[toJsonNumber, order2InterpolationExtraction["diagnostics"]],
  "max_relative_abs_delta" -> toJsonNumber[
    Max[("relative_abs_delta" /. Normal[#]) & /@ Values[interpolationMethodDeltas]]
  ]
|>;
operatorMetrics = Join[
  passChecks,
  <|
    "max_green_residual" -> toJsonNumber[primaryTransfer["maxGreenResidual"]],
    "max_bdg_amp_residual" -> toJsonNumber[primaryTransfer["maxBdgAmpResidual"]],
    "max_bdg_phase_residual" -> toJsonNumber[primaryTransfer["maxBdgPhaseResidual"]],
    "max_bdg_eigen_residual" -> toJsonNumber[bdgActive["max_relative_eigen_residual"]],
    "max_continuity_residual" -> toJsonNumber[primaryTransfer["maxContinuityResidual"]],
    "pure_gauge_physical_field_norm" -> toJsonNumber[pureGauge["physical_field_norm"]],
    "pure_gauge_outgoing_flux_abs" -> toJsonNumber[Abs[pureGauge["outgoing_flux"]]],
    "basis_max_relative_Z_delta" -> toJsonNumber[basisInvariant["max_relative_Z_delta"]],
    "basis_max_relative_N_delta" -> toJsonNumber[basisInvariant["max_relative_N_delta"]],
    "v2_09_relative_Z_error" -> toJsonNumber[v209["relative_Z_error"]],
    "v2_09_relative_N0_error" -> toJsonNumber[v209["relative_N0_error"]],
    "hard_cap_to_open_boundary_norm_ratio" -> toJsonNumber[capContrast["hard_cap_to_open_boundary_norm_ratio"]]
  |>
];

sourceHashes = <|
  "torch_background_content_hash" -> background["content_hash"],
  "torch_background_sha256" -> sha256File[backgroundPath],
  "mt15_06_script_sha256" -> sha256File[scriptPath],
  "freeze_sheet_sha256" -> If[physicalMode, sha256File[physicalFreezePath], "not_physical"]
|>;
packetBranchId = If[physicalMode,
  "m1c_decision_07_physical_frozen_direct_derived",
  "mt15_06_m1c_prep_imported_torch_background_direct_derived"
];
packetParentStatus = If[physicalMode,
  "effective_closure",
  "m1c_prep_effective_closure_placeholder_unfrozen"
];
packetMachineryLabel = If[physicalMode,
  "M1c PHYSICAL frozen decision-07 effective-closure branch; target-blind freeze hash written before solve.",
  "M1c-prep machinery on real self-consistent WP1 background but placeholder unfrozen params; NOT physical."
];

derivedBdgWall = <|
  "status" -> If[physicalMode, "derived_bdg_wall_m1c_physical_frozen_imported_torch_background", "derived_bdg_wall_m1c_prep_imported_torch_background"],
  "coefficients" -> <|"K" -> toJsonNumber[Kwall], "M" -> toJsonNumber[Mwall], "B0" -> toJsonNumber[B0], "B2" -> toJsonNumber[B2], "B4" -> toJsonNumber[B4]|>,
  "machinery_label" -> packetMachineryLabel,
  "source_paths" -> <|"torch_background" -> backgroundPath|>,
  "source_hashes" -> sourceHashes
|>;

derivedMaxwellTransfer = <|
  "status" -> "derived_green_function_transfer",
  "gauge_convention" -> "localized_H_equals_Z_imported_torch_background_A0_active_spatial_gauge_measured",
  "flux_normalization" -> <|
    "Gamma_port" -> toJsonNumber[gammaPort],
    "formula" -> "sigma_Q^can=9/(8 Omega_Q^5)=4 a^5/(27 c_s^5)",
    "legacy_v2_09_reduced_formula" -> "a^5/(27 c_s^5)",
    "legacy_factor_gap" -> toJsonNumber[gammaPortFactorVsLegacy],
    "derived_not_fitted_to_chi_Q" -> True
  |>,
  "coefficients" -> <|"Z0" -> toJsonNumber[Z0], "Z2" -> toJsonNumber[Z2], "Z4" -> toJsonNumber[Z4], "N0" -> toJsonNumber[N0], "N2" -> toJsonNumber[N2], "N4" -> toJsonNumber[N4]|>,
  "operator_gauge_residual_metrics" -> operatorMetrics,
  "source_paths" -> <|"torch_background" -> backgroundPath|>,
  "source_hashes" -> sourceHashes,
  "machinery_label" -> packetMachineryLabel
|>;

packet = <|
  "schema" -> "stage_v2_22b_solver_handoff/v1",
  "branch_id" -> packetBranchId,
  "freeze" -> <|
    "pre_target_freeze" -> True,
    "target_blind" -> True,
    "parent_action" -> "GNLS imported self-consistent torch WP1 background plus M1b BdG/wall and Spike-2 direct Maxwell transfer",
    "parent_action_status" -> packetParentStatus,
    "freeze_hash" -> If[physicalMode, background["freeze_hash"], "not_physical"],
    "m1c_physical_frozen_run" -> physicalMode,
    "claim_scope" -> If[physicalMode, "pre_registered_effective_closure_branch", "m1c_prep_machinery_only"],
    "gauge_convention" -> "localized_H_equals_Z_imported_torch_background_A0_active",
    "boundary_protocol" -> "open_impedance_AC_reflecting_DC_leaking",
    "source_map_convention" -> "canonical real STF grouped P2 basis",
    "support_profile_family" -> "derived weak-form wall eigen-shape on imported torch WP1 geometry"
  |>,
  "geometry" -> <|
    "L" -> toJsonNumber[wLength],
    "R_mouth" -> toJsonNumber[If[physicalMode, aConst, r0Geometry[wMin]]],
    "R_exit" -> toJsonNumber[If[physicalMode, physicalFreeze["geometry"]["R_exit"], r0Geometry[wMax]]],
    "boundary_class" -> "open_impedance",
    "Y_L_limit" -> 0.0,
    "exit_model" -> "impedance_mismatch_open_exit"
  |>,
  "constants" -> <|"G" -> GConst, "c_s" -> csConst, "c" -> cConst, "a" -> toJsonNumber[aConst], "mhat0" -> mhat0Const, "S_port" -> SPortConst, "theta_tail" -> thetaTailConst|>,
  "grid" -> <|"coordinate" -> "w", "points" -> toJsonVector[wallMode["packet_points"]]|>,
  "wall" -> <|"K" -> toJsonNumber[Kwall], "M" -> toJsonNumber[Mwall], "derivation_status" -> "derived_weak_form_galerkin_on_imported_background_geometry"|>,
  "profiles" -> <|"radius_R0" -> toJsonVector[wallMode["radius"]], "weight" -> toJsonVector[wallMode["weight"]], "wall_chi_eta" -> toJsonVector[wallMode["chi"]]|>,
  "bdg_modes" -> bdgActive["packet_modes"],
  "derived_bdg_wall_coefficients" -> derivedBdgWall,
  "derived_maxwell_transfer" -> derivedMaxwellTransfer,
  "direct_derived_extension" -> <|"status" -> "v2_22b_direct_derived_coefficient_path", "no_posited_ports_in_N0_path" -> True, "v2_21_lane" -> "direct_coefficients", "target_blind_constants" -> True|>,
  "solver_metadata" -> <|
    "exporter" -> "mt15_06_m1c_prep_crossengine.wls",
    "mesh_points" -> Length[wallMode["packet_points"]],
    "coefficient_family" -> If[physicalMode, "m1c_physical_decision_07_frozen_bdg_wall_plus_spike2_direct_transfer_v1", "m1c_prep_imported_torch_background_bdg_wall_plus_spike2_direct_transfer_v1"],
    "source_commit" -> sourceCommitString[],
    "nonlinear_residual_norm" -> toJsonNumber[background["residuals"]["coupled_stationary_linf"]],
    "linear_residual_norm" -> toJsonNumber[bdgActive["max_relative_eigen_residual"]],
    "continuation_residual_norm" -> toJsonNumber[background["residuals"]["coupled_stationary_linf"]],
    "bdg_status" -> "derived_from_imported_background_finite_dimensional_BdG_eigensolve",
    "wall_status" -> "derived_from_weak_form_integrals",
    "maxwell_transfer_status" -> "derived_green_function_transfer_spike2_on_imported_background",
    "mixed_ports_status" -> "absent_for_direct_derived_N0_path",
    "notes" -> If[physicalMode, "M1c PHYSICAL frozen packet: decision-07 freeze hash was written before solve; no target residual outputs are present in this V2-22B handoff.", "M1c-prep only: real self-consistent torch WP1 background with PLACEHOLDER unfrozen params. This packet proves machinery and contains no target residual outputs."]
  |>,
  "normalization_tolerances" -> <|"profile_norm_tol" -> 0.005, "Delta_tol" -> 1.0*^-12|>
|>;
packet = Join[packet, <|"content_hash" -> contentHash[packet]|>];

diagnostics = <|
  "schema" -> If[physicalMode, "m1c_physical_frozen_derived_diagnostics/v1", "mt15_06_m1c_prep_crossengine_diagnostics/v1"],
  "background_path" -> backgroundPath,
  "packet_path" -> packetPath,
  "freeze_hash" -> If[physicalMode, background["freeze_hash"], "not_physical"],
  "m1c_physical_frozen_run" -> physicalMode,
  "claim_scope" -> If[physicalMode, "pre_registered_effective_closure_branch", "m1c_prep_machinery_only"],
  "machinery_label" -> If[physicalMode, "PHYSICAL frozen decision-07 effective-closure branch; conditional on frozen target-blind posits; no value mutated after residual.", "LOUD LABEL: real self-consistent torch WP1 background but PLACEHOLDER unfrozen params -> MACHINERY, NOT physical; proves cross-engine pipeline only."],
  "torch_background" -> <|
    "content_hash" -> background["content_hash"],
    "grid" -> {gridBg["nr"], gridBg["nw"]},
    "residual_linf" -> toJsonNumber[background["residuals"]["coupled_stationary_linf"]],
    "converged" -> background["solver"]["converged"],
    "A_00_linf" -> toJsonNumber[sourceA0Linf],
    "A_r0_linf" -> toJsonNumber[sourceArLinf],
    "A_w0_linf" -> toJsonNumber[sourceAwLinf]
  |>,
	  "import_fidelity" -> <|
	    "smooth_field_interpolation_order" -> 1,
	    "density_projection" -> "conservative piecewise-constant finite-volume overlap using r^2 dr dw measure without the common 4pi factor",
	    "interpolation_convergence_status" -> "not_claimed_for_m1c_prep",
	    "source_rho_mass_no_4pi" -> toJsonNumber[sourceRhoMass],
	    "bdg_projected_rho_mass_relative_error" -> toJsonNumber[bdgActive["projected_rho_mass_relative_error"]],
	    "source_point_self_match" -> <|
	      "note" -> "Interpolant reproduces its own source grid nodes; this is trivially exact by construction and is NOT an interpolation-convergence study.",
	      "psi_R0_max_abs" -> toJsonNumber[fieldSelfMatchMax],
	      "R0_w_max_abs" -> toJsonNumber[r0SelfMatchMax]
	    |>
	  |>,
	  "interpolation_method_sensitivity" -> interpolationMethodSensitivity,
    "direct_coefficients" -> Map[toJsonNumber, directCoefficients],
    "derived_bdg_wall_coefficients" -> derivedBdgWall,
    "derived_maxwell_transfer" -> derivedMaxwellTransfer,
	  "bdg_wall" -> <|
    "A0_active" -> True,
    "A0_zero_comparison" -> A0Changed,
    "grid" -> bdgActive["grid"],
    "positive_real_modes_found" -> bdgActive["positive_real_modes_found"],
    "modes_exported" -> bdgActive["modes_exported"],
    "max_relative_eigen_residual" -> toJsonNumber[bdgActive["max_relative_eigen_residual"]],
    "stieltjes_gap" -> toJsonNumber[bdgActive["stieltjes_gap"]],
    "stieltjes_pass" -> bdgActive["stieltjes_pass"],
    "B0" -> toJsonNumber[B0], "B2" -> toJsonNumber[B2], "B4" -> toJsonNumber[B4],
    "wall_K" -> toJsonNumber[Kwall], "wall_M" -> toJsonNumber[Mwall], "wall_M_positive" -> TrueQ[wallMode["MMinEigenvalue"] > 0]
  |>,
  "spike2_transfer" -> <|
    "current_frechet" -> currentCheck,
    "primary_grid" -> primaryTransfer["grid"],
    "mesh_convergence_rows" -> (Function[row, <|
        "grid" -> row["grid"],
        "h" -> toJsonNumber[row["h"]],
        "Z0" -> toJsonNumber[row["Z0"]],
        "Z2" -> toJsonNumber[row["Z2"]],
        "Z4" -> toJsonNumber[row["Z4"]],
        "N0" -> toJsonNumber[row["N0"]],
        "N2" -> toJsonNumber[row["N2"]],
        "N4" -> toJsonNumber[row["N4"]],
        "maxGreenResidual" -> toJsonNumber[row["maxGreenResidual"]],
        "maxBdgAmpResidual" -> toJsonNumber[row["maxBdgAmpResidual"]],
        "maxBdgPhaseResidual" -> toJsonNumber[row["maxBdgPhaseResidual"]],
        "maxContinuityResidual" -> toJsonNumber[row["maxContinuityResidual"]]
      |>] /@ meshTransfers),
    "coefficients" -> Map[toJsonNumber, <|"Z0" -> Z0, "Z2" -> Z2, "Z4" -> Z4, "N0" -> N0, "N2" -> N2, "N4" -> N4|>],
    "mesh_orders" -> Map[toJsonNumber, meshZOrders],
    "pure_gauge" -> pureGauge,
    "basis_invariance" -> basisInvariant,
    "v2_09_regression" -> v209,
    "hard_cap_contrast" -> capContrast,
    "spatial_gauge_sign_sensitivity" -> spatialGaugeSensitivity,
    "pass_checks" -> passChecks
  |>,
  "direct_formula_preview_not_exported_in_packet" -> <|"D0" -> toJsonNumber[d0], "P2" -> toJsonNumber[p2], "P4" -> toJsonNumber[p4], "R_pole" -> toJsonNumber[rPole], "R_norm" -> toJsonNumber[rNormMachinery]|>,
  "content_hashes" -> <|"packet_content_hash" -> packet["content_hash"], "background_content_hash" -> background["content_hash"]|>,
  "source_hashes" -> sourceHashes
|>;
diagnostics = Join[diagnostics, <|"content_hash" -> contentHash[diagnostics]|>];

Export[packetPath, packet, "JSON"];
Export[diagnosticsPath, diagnostics, "JSON"];

gateRows = ({#, If[passChecks[#], "PASS", "FAIL"]} & /@ Keys[passChecks]);
coeffRows = {
  {"K", formatScientific[Kwall]}, {"M", formatScientific[Mwall]},
  {"B0", formatScientific[B0]}, {"B2", formatScientific[B2]}, {"B4", formatScientific[B4]},
  {"Z0", formatScientific[Z0]}, {"Z2", formatScientific[Z2]}, {"Z4", formatScientific[Z4]},
  {"N0", formatScientific[N0]}, {"N2", formatScientific[N2]}, {"N4", formatScientific[N4]}
};
interpSensitivityRows = ({
    #,
    formatScientific[interpolationMethodDeltas[#]["order1_primary"]],
    formatScientific[interpolationMethodDeltas[#]["order2_diagnostic"]],
    formatScientific[interpolationMethodDeltas[#]["relative_delta_order2_minus_order1"]],
    formatScientific[interpolationMethodDeltas[#]["relative_abs_delta"]]
  } & /@ interpSensitivityKeys);

report = StringRiffle[
  {
    If[physicalMode, "# M1c Physical Frozen Derived Chain Report", "# MT15-06 M1c-Prep Cross-Engine Background Report"],
    "",
    If[physicalMode,
      "**LOUD LABEL:** PHYSICAL frozen decision-07 effective-closure branch. The freeze hash was written before the WP1 solve; the result is conditional on the frozen target-blind posits and modest CPU resolution.",
      "**LOUD LABEL:** real self-consistent torch WP1 background, but PLACEHOLDER unfrozen params. This is MACHINERY, NOT physical; it proves the cross-engine pipeline only. No Gate A, no frozen committed packet, no physical claim."
    ],
    "",
    "## Torch Background",
    "",
    "- Background JSON: `" <> backgroundPath <> "`",
    "- Content hash: `" <> background["content_hash"] <> "`",
    "- Grid: `" <> ToString[gridBg["nr"]] <> "x" <> ToString[gridBg["nw"]] <> "`",
    "- Coupled stationary residual linf: `" <> formatScientific[background["residuals"]["coupled_stationary_linf"]] <> "`",
    "- Gauge field linf `(A_00,A_r0,A_w0)`: `" <> StringRiffle[formatScientific /@ {sourceA0Linf, sourceArLinf, sourceAwLinf}, ", "] <> "`",
    "",
    "## Import Fidelity",
    "",
    "- Smooth-field interpolation order for the primary extraction: `1`.",
    "- Conservative rho projection: finite-volume overlap in the imported `r^2 dr dw` measure.",
    "- Source-point self-match only: the interpolant reproduces its own source grid nodes, so this is trivially exact by construction and is NOT an interpolation-convergence study.",
    "- Source-point self-match `psi_R0` max abs delta: `" <> formatScientific[fieldSelfMatchMax] <> "`",
    "- Source-point self-match `R0(w)` max abs delta: `" <> formatScientific[r0SelfMatchMax] <> "`",
    "- BdG projected rho mass relative error: `" <> formatScientific[bdgActive["projected_rho_mass_relative_error"]] <> "`",
    "",
    "## Interpolation-Method Sensitivity",
    "",
    "Diagnostic-only re-extraction on the same torch background using order-2 imported smooth-field interpolation; the primary order-1 extraction above remains unchanged.",
    "",
    markdownTable[{"coefficient", "order1 primary", "order2 diagnostic", "relative delta", "relative abs delta"}, interpSensitivityRows],
    "",
    "- Max relative abs delta: `" <> formatScientific[interpolationMethodSensitivity["max_relative_abs_delta"]] <> "`",
    "- Carry-forward note: " <> interpolationMethodCarryForwardNote,
    "",
    "## Gauge Activation",
    "",
    "The imported background activates scalar `A_00` in the BdG `q A0` diagonal. The isotropic WP1 solve returns spatial gauges with linf values shown above; they are numerically zero at this resolution, so the spatial-gauge sign flip is a near-null sensitivity check rather than a large effect.",
    "",
    "- BdG active-minus-zero-A0 deltas `(B0,B2,B4)`: `" <> StringRiffle[formatScientific /@ {A0Changed["B0_delta_active_minus_zeroA0"], A0Changed["B2_delta_active_minus_zeroA0"], A0Changed["B4_delta_active_minus_zeroA0"]}, ", "] <> "`",
    "- Spatial gauge sign-flip relative deltas `(Z0,N0)`: `" <> StringRiffle[formatScientific /@ {spatialGaugeSensitivity["relative_Z0_delta"], spatialGaugeSensitivity["relative_N0_delta"]}, ", "] <> "`",
    "",
    "## Coefficients",
    "",
    markdownTable[{"coefficient", "value"}, coeffRows],
    "",
    "## A-not-zero Revalidation",
    "",
    markdownTable[{"check", "status"}, gateRows],
    "",
    "- BdG max relative eigen residual: `" <> formatScientific[bdgActive["max_relative_eigen_residual"]] <> "`",
    "- Stieltjes `B0*B4-B2^2`: `" <> formatScientific[bdgActive["stieltjes_gap"]] <> "`",
    "- Current Frechet max relative error: `" <> formatScientific[currentCheck["max_relative_error"]] <> "`",
    "- Green residual max: `" <> formatScientific[primaryTransfer["maxGreenResidual"]] <> "`",
    "- Pure gauge physical field norm / flux: `" <> formatScientific[pureGauge["physical_field_norm"]] <> " / " <> formatScientific[pureGauge["outgoing_flux"]] <> "`",
    "- Basis invariance max relative Z/N: `" <> formatScientific[basisInvariant["max_relative_Z_delta"]] <> " / " <> formatScientific[basisInvariant["max_relative_N_delta"]] <> "`",
    "- V2-09 relative Z/N errors: `" <> formatScientific[v209["relative_Z_error"]] <> " / " <> formatScientific[v209["relative_N0_error"]] <> "`",
    "- N0: `" <> formatScientific[N0] <> "`",
    "",
    "## Direct Formula Preview",
    "",
    If[physicalMode,
      "These values are physical-run diagnostics computed after the freeze; they are not exported inside the V2-22B handoff packet because the validator forbids target residual fields there.",
      "These preview values are diagnostics only and are not exported inside the V2-22B packet."
    ],
    "",
    "- `D0 = " <> formatNumber[d0] <> "`",
    "- `R_norm = " <> formatNumber[rNormMachinery] <> "`" <> If[physicalMode, " (PHYSICAL frozen diagnostic; V2 remains authoritative)", " (MACHINERY, NOT physical)"],
    "- `R_pole = " <> formatNumber[rPole] <> "`",
    "- `P2 = " <> formatNumber[p2] <> "`",
    "- `P4 = " <> formatNumber[p4] <> "`",
    "",
    "## Artifacts",
    "",
    "- Direct-derived packet: `" <> packetPath <> "`",
    "- Diagnostics JSON: `" <> diagnosticsPath <> "`",
    "- Packet content hash: `" <> packet["content_hash"] <> "`",
    "- Diagnostics content hash: `" <> diagnostics["content_hash"] <> "`"
  },
  "\n"
];
Export[reportPath, report, "Text"];

Print["packet_path: ", packetPath];
Print["diagnostics_path: ", diagnosticsPath];
Print["report_path: ", reportPath];
Print["packet_content_hash: ", packet["content_hash"]];
Print["diagnostics_content_hash: ", diagnostics["content_hash"]];
Print["pass_checks: ", passChecks];
If[physicalMode,
  Print["R_norm_physical_frozen_preview: ", formatNumber[rNormMachinery]],
  Print["R_norm_machinery_not_physical: ", formatNumber[rNormMachinery]]
];
