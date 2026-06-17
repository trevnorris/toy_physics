#!/usr/bin/env wolframscript

(* MT15-04 / Spike-2: open-exit DtN, forward wall->matter->gauge
   source current, and basis-invariant Z0/N0 extraction from the retarded
   Maxwell self-energy.

   Scope is deliberately limited to Decision-05 Spike-2:
   - no V2-22B packet emission,
   - no R_norm computation,
   - no S_eta^(A) return source,
   - target-blind M1b smoke background only.
*)

$HistoryLength = 0;

ClearAll["Global`*"];

scriptPath =
  If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
    $InputFileName,
    FileNameJoin[{Directory[], "software/stage1_solver/mathematica/mt15_04_spike2_transfer_n0.wls"}]
  ];
scriptDir = DirectoryName[scriptPath];
repoRoot = AbsoluteFileName[FileNameJoin[{scriptDir, "..", "..", ".."}]];
runDir = FileNameJoin[{scriptDir, "runs", "mt15_04_spike2_transfer_n0"}];
If[! DirectoryQ[runDir], CreateDirectory[runDir, CreateIntermediateDirectories -> True]];

diagnosticsPath = FileNameJoin[{runDir, "mt15_04_spike2_transfer_n0_diagnostics.json"}];
reportPath = FileNameJoin[{scriptDir, "mt15_04_spike2_transfer_n0_report.md"}];

ClearAll[
  smootherStep, trapz, weightedInner, weightedNorm, normalizeProfile,
  toJsonNumber, toJsonVector, toJsonMatrix, sourceCommitString,
  diagSparse, firstDerivativeMatrix, secondDerivativeMatrix, gridData,
  blockRow, weightedGram, boundaryMatrix, buildMaxwellSystem,
  buildWallMode, buildMatterSource, linearizedCurrentFormula,
  nonlinearNumberCurrent, currentFrechetCheck, selfEnergySolve,
  fitEvenCoefficients, evenFitConditionNumber, transferForGrid, estimateOrderTriplet,
  fieldNormFromComponents, makeGaugeVector, pureGaugeCheck,
  basisInvarianceCheck, v209Regression, hardCapContrast,
  formatScientific, formatFixed, markdownTable, finiteQ
];

ell = 2;
Lharm = ell (ell + 1);
xiGauge = 1.0;

(* Target-blind M1b geometry and smoke constants. *)
wLength = 2.0;
aMouth = 1.0;
rExit = 1.65;
rMinTransfer = 0.28;
zFloor = 0.08;
zAmp = 0.92;
zCenter = 0.5 wLength;
zSigma = 0.38 wLength;

hbar = 1.0;
particleMass = 1.0;
gaugeCharge = 0.35;
mu0 = 1.0;
eosK = 0.04;
chemicalPotential = 0.05;
radialWallStrength = 0.65;
axialTrapStrength = 0.12;
phaseRegularizer = 0.05;

GConst = 1.0;
csConst = 1.0;
cConst = 2.5;
aConst = aMouth;

(* D4 normalization pinning.  The compact canonical retarded quadrupole
   denominator uses sigma_Q^can=9/(8 Omega_Q^5), Omega_Q=3 c_s/(2a).
   Therefore Gamma_port is 4 a^5/(27 c_s^5).  The older V2-09
   a^5/(27 c_s^5) is the c1=1/4 pole-residue-reduced coefficient, not the
   physical outgoing-port denominator normalization used here. *)
omegaQ = 3.0 csConst/(2.0 aConst);
gammaPort = 9.0/(8.0 omegaQ^5);
gammaPortLegacyV209 = aConst^5/(27.0 csConst^5);
gammaPortCompact = 4.0 aConst^5/(27.0 csConst^5);
gammaPortFactorVsLegacy = gammaPort/gammaPortLegacyV209;

(* Frozen target-blind outgoing DtN/admittance branch.  It is deliberately
   declared before any transfer solve.  The i sign is chosen for an e^{-iwt}
   convention so K_ret = K_cons + i Gamma_port omega^5 B gives
   Im <j,G_ret j> <= 0 and positive outgoing flux. *)
yOutBranchDescription =
  "Frozen target-blind D2 branch: Y_out,l=2(omega)=i Gamma_port omega^5 on physical tangential lanes at w=L and r=R_exit; Gamma_port=sigma_Q^can=4 a^5/(27 c_s^5).";

smootherStep[x_] := 3 x^2 - 2 x^3;
r0Geometry[w_] := aMouth + (rExit - aMouth) smootherStep[w/wLength];
zWeight[w_] := zFloor + zAmp Exp[-((w - zCenter)/zSigma)^2];

trapz[xs_List, ys_List] /; Length[xs] == Length[ys] && Length[xs] >= 2 :=
  Total[
    Table[
      0.5 (ys[[i]] + ys[[i + 1]]) (xs[[i + 1]] - xs[[i]]),
      {i, 1, Length[xs] - 1}
    ]
  ];

weightedInner[xs_List, weights_List, f_List, g_List] := trapz[xs, weights f g];
weightedNorm[xs_List, weights_List, f_List] := weightedInner[xs, weights, f, f];

normalizeProfile[xs_List, weights_List, f_List, label_String] := Module[{norm = weightedNorm[xs, weights, f]},
  If[! NumericQ[norm] || norm <= 0,
    Print["ERROR: nonpositive norm while normalizing ", label, ": ", norm];
    Exit[2]
  ];
  N[f/Sqrt[norm], 17]
];

finiteQ[x_] := Quiet[Check[NumberQ[N[x]] && TrueQ[Abs[N[x]] < Infinity], False]];

toJsonNumber[x_] := Module[{n = Quiet[Check[N[x], $Failed]]},
  If[n === $Failed || ! NumericQ[n],
    Print["ERROR: nonnumeric JSON export value: ", InputForm[x]];
    Exit[4]
  ];
  If[Abs[Im[n]] < 10^-12, N[Re[n], 17], <|"re" -> N[Re[n], 17], "im" -> N[Im[n], 17]|>]
];
toJsonVector[x_List] := toJsonNumber /@ x;
toJsonMatrix[x_List] := Map[toJsonNumber, x, {2}];

sourceCommitString[] := Module[
  {gitDir, headPath, head, refName, refPath, packedPath, packedLines, match, sha},
  gitDir = FileNameJoin[{repoRoot, ".git"}];
  headPath = FileNameJoin[{gitDir, "HEAD"}];
  If[! FileExistsQ[headPath], Return["unavailable"]];

  head = Quiet[Check[StringTrim[Import[headPath, "Text"]], ""]];
  If[! StringQ[head] || StringLength[head] == 0, Return["unavailable"]];

  sha = If[StringStartsQ[head, "ref:"],
    refName = StringTrim[StringDrop[head, 4]];
    refPath = FileNameJoin[Join[{gitDir}, FileNameSplit[refName]]];
    If[FileExistsQ[refPath],
      Quiet[Check[StringTrim[Import[refPath, "Text"]], ""]],
      packedPath = FileNameJoin[{gitDir, "packed-refs"}];
      If[! FileExistsQ[packedPath],
        "",
        packedLines = StringSplit[Quiet[Check[Import[packedPath, "Text"], ""]], "\n"];
        match = Select[packedLines, StringMatchQ[#, Except["#"] ~~ ___ ~~ " " ~~ refName] &];
        If[Length[match] > 0, First[StringSplit[First[match]]], ""]
      ]
    ],
    head
  ];

  If[StringQ[sha] && StringLength[sha] >= 7, StringTake[sha, 7], "unavailable"]
];

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

gridData[nr_Integer, nw_Integer] := Module[
  {dr, dw, rVals1, wVals1, rVals, wVals, zVals, r0Vals, quadWeights},
  dr = (rExit - rMinTransfer)/(nr + 1);
  dw = wLength/(nw + 1);
  rVals1 = N[Table[rMinTransfer + i dr, {i, 1, nr}], 17];
  wVals1 = N[Table[i dw, {i, 1, nw}], 17];
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
    "z" -> zVals,
    "R0" -> r0Vals,
    "w0" -> quadWeights,
    "h" -> Max[dr, dw]
  |>
];

blockRow[blocks_List] := ArrayFlatten[{blocks}];
weightedGram[m_, wdiag_, sign_] := sign ConjugateTranspose[m].wdiag.m;

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
    (* w=L exit: physical tangential electric fields E_r, E_E, E_B. *)
    If[j == nw,
      wExitWeight = grid["dr"] grid["r"][[cell]]^2 grid["z"][[cell]];
      wExitScalar[[cell]] += wExitWeight;
      wExitAngular[[cell]] += Lharm wExitWeight;
    ];
    (* radial outer boundary: physical tangential electric fields E_E, E_B, E_w. *)
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
    weightedGram[mEr, w0Diag, 1.0] +
    weightedGram[mEE, lw0Diag, 1.0] +
    weightedGram[mEB, lw0Diag, 1.0] +
    weightedGram[mEw, w0Diag, 1.0] -
    weightedGram[mCr, w0Diag, 1.0] -
    weightedGram[mCE, lw0Diag, 1.0] -
    weightedGram[mCB, lw0Diag, 1.0] -
    weightedGram[mBB, lw0Diag, 1.0] -
    weightedGram[mBr, w0Diag, 1.0] -
    weightedGram[mBE, lw0Diag, 1.0] -
    (1.0/xiGauge) weightedGram[mG, w0Diag, 1.0];

  kFieldPositive =
    weightedGram[mEr, w0Diag, 1.0] +
    weightedGram[mEE, lw0Diag, 1.0] +
    weightedGram[mEB, lw0Diag, 1.0] +
    weightedGram[mEw, w0Diag, 1.0] +
    weightedGram[mCr, w0Diag, 1.0] +
    weightedGram[mCE, lw0Diag, 1.0] +
    weightedGram[mCB, lw0Diag, 1.0] +
    weightedGram[mBB, lw0Diag, 1.0] +
    weightedGram[mBr, w0Diag, 1.0] +
    weightedGram[mBE, lw0Diag, 1.0];

  boundaryB = boundaryMatrix[
    grid,
    <|"E_r" -> mEr, "E_E" -> mEE, "E_B" -> mEB, "E_w" -> mEw|>
  ];
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
      "omega" -> omega,
      "total" -> total,
      "Dr" -> drMat,
      "Dw" -> dwMat,
      "Drr" -> drrMat,
      "Dww" -> dwwMat,
      "W0" -> w0Diag,
      "WLane" -> wLane,
      "WLaneInv" -> wLaneInv,
      "WLaneValues" -> wLaneVals,
      "KAction" -> kAction,
      "KCons" -> kCons,
      "KRet" -> kRet,
      "ACons" -> wLaneInv . kCons,
      "ARet" -> wLaneInv . kRet,
      "BoundaryB" -> boundaryB,
      "BoundaryLossStrong" -> wLaneInv . boundaryB,
      "AnchorMatrix" -> anchorMatrix,
      "ScalarGaugeAnchor" -> <|"cell" -> anchorCell, "strength" -> 10^-8|>,
      "KFieldPositive" -> kFieldPositive,
      "FieldMatrices" -> <|
        "G" -> mG,
        "E_r" -> mEr, "E_E" -> mEE, "E_B" -> mEB, "E_w" -> mEw,
        "C_r" -> mCr, "C_E" -> mCE, "C_B" -> mCB,
        "B_B" -> mBB, "B_r" -> mBr, "B_E" -> mBE
      |>,
      "Offsets" -> offsets
    |>
  ]
];

buildWallMode[] := Module[
  {
    meshPoints = 129, points, xGrid, radiusValues, weightValues, muEta,
    Tw, TOmega, KEta, Vl, YL0, basisCount = 3, basisValues,
    basisDerivatives, wallMMatrix, wallKMatrix, wallEigenvalues,
    wallEigenvectors, wallPairs, wallCoeffRaw, wallShapeRaw,
    wallWeightNorm, wallCoeff, wallChiEta
  },
  points = N[Table[wLength i/(meshPoints - 1), {i, 0, meshPoints - 1}], 17];
  xGrid = points/wLength;
  radiusValues = N[aMouth + (rExit - aMouth) (smootherStep /@ xGrid), 17];
  weightValues = N[radiusValues^2, 17];
  muEta = N[1.0 + 0.08 Cos[Pi xGrid]^2, 17];
  Tw = N[1.1 + 0.10 Sin[Pi xGrid], 17];
  TOmega = N[0.35 + 0.03 Cos[2 Pi xGrid], 17];
  KEta = N[0.45 + 0.08 256 xGrid^4 (1 - xGrid)^4, 17];
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
    "radius" -> radiusValues,
    "weight" -> weightValues,
    "chi" -> N[wallChiEta, 17],
    "interp" -> Interpolation[Transpose[{points, N[wallChiEta, 17]}], InterpolationOrder -> 1],
    "K" -> N[wallCoeff . wallKMatrix . wallCoeff, 17],
    "M" -> N[wallCoeff . wallMMatrix . wallCoeff, 17],
    "MMatrix" -> wallMMatrix,
    "KMatrix" -> wallKMatrix,
    "coeff" -> wallCoeff,
    "norm" -> N[weightedNorm[points, weightValues, wallChiEta], 17]
  |>
];

wallMode = buildWallMode[];

linearizedCurrentFormula[
  sys_Association, psiR0_List, psiI0_List, dpsiR_List, dpsiI_List,
  ar0_List, aw0_List, dar_List, daw_List
] := Module[
  {
    gradPsiRr, gradPsiIr, gradPsiRw, gradPsiIw, gradDPsiRr, gradDPsiIr,
    gradDPsiRw, gradDPsiIw, density0, deltaDensity, phaseR, phaseW,
    gaugeR, gaugeW, jr, jw
  },
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

nonlinearNumberCurrent[
  sys_Association, psiR_List, psiI_List, ar_List, aw_List
] := Module[
  {density, jr, jw},
  density = psiR^2 + psiI^2;
  jr = (hbar/particleMass) (psiR (sys["Dr"] . psiI) - psiI (sys["Dr"] . psiR)) -
    (gaugeCharge/particleMass) ar density;
  jw = (hbar/particleMass) (psiR (sys["Dw"] . psiI) - psiI (sys["Dw"] . psiR)) -
    (gaugeCharge/particleMass) aw density;
  <|"jr" -> jr, "jw" -> jw|>
];

buildMatterSource[sys_Association, omega_?NumericQ] := Module[
  {
    n = sys["n"], idN, invRDiag, lap, bdgRadii, bdgW, bdgR0, vconf,
    rho0, psi0, psiI0, enthalpy, enthalpyPrime, enthalpyLin, kinetic,
    single, lPlus, lPhase, etaVals, deltaV, drive, dpsiR, dpsiI,
    zero, current, deltaRho, jr, jw, jE, jVec, ampResidual, phaseResidual,
    sourceNorm, continuityResidual
  },
  idN = SparseArray[Band[{1, 1}] -> 1.0, {n, n}];
  invRDiag = diagSparse[1.0/sys["r"]];
  lap = sys["Drr"] + diagSparse[2.0/sys["r"]] . sys["Dr"] + sys["Dww"] -
    diagSparse[Lharm/sys["r"]^2];
  bdgRadii = sys["r"];
  bdgW = sys["w"];
  bdgR0 = sys["R0"];
  vconf = N[
    radialWallStrength (bdgRadii/bdgR0)^4 +
      0.5 axialTrapStrength ((bdgW - 0.5 wLength)/wLength)^2,
    17
  ];
  rho0 = N[
    0.035 + 0.18 Exp[-(bdgRadii/(0.72 bdgR0))^2] (1.0 + 0.06 Cos[Pi bdgW/wLength])^2,
    17
  ];
  psi0 = N[Sqrt[rho0], 17];
  psiI0 = ConstantArray[0.0, n];
  enthalpy = N[(5 eosK/4) rho0^4, 17];
  enthalpyPrime = N[5 eosK rho0^3, 17];
  enthalpyLin = N[rho0 enthalpyPrime, 17];
  kinetic = N[-(hbar^2/(2 particleMass)) lap, 17];
  single = kinetic + diagSparse[vconf - chemicalPotential + enthalpy];
  lPlus = single + diagSparse[2.0 enthalpyLin] + (omega^2) idN;
  lPhase = single + phaseRegularizer idN + (omega^2) idN;
  etaVals = N[wallMode["interp"] /@ bdgW, 17];
  deltaV = N[-4.0 radialWallStrength bdgRadii^4/bdgR0^5 etaVals, 17];
  drive = N[-deltaV psi0, 17];
  dpsiR = LinearSolve[lPlus, drive];
  dpsiI = I omega LinearSolve[lPhase, dpsiR];
  zero = ConstantArray[0.0, n];
  current = linearizedCurrentFormula[sys, psi0, psiI0, dpsiR, dpsiI, zero, zero, zero, zero];
  deltaRho = current["deltaRho"];
  jr = current["jr"];
  jw = current["jw"];
  jE = (hbar/particleMass) psi0 dpsiI/sys["r"];
  jVec = Join[
    gaugeCharge deltaRho,
    gaugeCharge jr,
    gaugeCharge jE,
    zero,
    gaugeCharge jw
  ];
  ampResidual = Norm[lPlus . dpsiR - drive]/Max[Norm[drive], 10^-30];
  phaseResidual = Norm[lPhase . (dpsiI/(I Max[omega, 10^-30])) - dpsiR]/Max[Norm[dpsiR], 10^-30];
  sourceNorm = Sqrt[Re[Conjugate[jVec].(sys["WLane"].jVec)]];
  continuityResidual = Norm[-I omega deltaRho + sys["Dr"] . jr + sys["Dw"] . jw - Lharm jE/sys["r"]]/
    Max[Norm[deltaRho] + Norm[jr] + Norm[jw] + Norm[jE], 10^-30];
  <|
    "jVec" -> jVec,
    "deltaRho" -> deltaRho,
    "jr" -> jr,
    "jE" -> jE,
    "jB" -> zero,
    "jw" -> jw,
    "dpsiR" -> dpsiR,
    "dpsiI" -> dpsiI,
    "psi0" -> psi0,
    "deltaV" -> deltaV,
    "ampResidual" -> N[ampResidual, 17],
    "phaseResidual" -> N[phaseResidual, 17],
    "sourceNorm" -> N[sourceNorm, 17],
    "continuityResidual" -> N[continuityResidual, 17],
    "usedReturnSourceSEtaA" -> False
  |>
];

currentFrechetCheck[] := Module[
  {
    sys, n, rr, ww, psiR0, psiI0, ar0, aw0, dpsiR, dpsiI, dar, daw,
    eps = 10^-6, plus, minus, fdR, fdW, lin, errR, errW, normR, normW
  },
  sys = buildMaxwellSystem[5, 4, 0.0];
  n = sys["n"];
  rr = sys["r"];
  ww = sys["w"];
  psiR0 = 1.0 + 0.07 rr + 0.03 ww;
  psiI0 = 0.11 + 0.02 rr - 0.04 ww;
  ar0 = 0.04 rr + 0.01 ww;
  aw0 = -0.03 rr + 0.02 ww;
  dpsiR = 0.05 rr^2 - 0.02 ww;
  dpsiI = -0.03 rr + 0.04 ww^2;
  dar = 0.02 rr (1.0 + ww);
  daw = -0.015 (1.0 + rr) ww;
  plus = nonlinearNumberCurrent[
    sys, psiR0 + eps dpsiR, psiI0 + eps dpsiI, ar0 + eps dar, aw0 + eps daw
  ];
  minus = nonlinearNumberCurrent[
    sys, psiR0 - eps dpsiR, psiI0 - eps dpsiI, ar0 - eps dar, aw0 - eps daw
  ];
  fdR = (plus["jr"] - minus["jr"])/(2.0 eps);
  fdW = (plus["jw"] - minus["jw"])/(2.0 eps);
  lin = linearizedCurrentFormula[sys, psiR0, psiI0, dpsiR, dpsiI, ar0, aw0, dar, daw];
  errR = Norm[lin["jr"] - fdR];
  errW = Norm[lin["jw"] - fdW];
  normR = Max[Norm[fdR], 10^-30];
  normW = Max[Norm[fdW], 10^-30];
  <|
    "grid" -> {sys["nr"], sys["nw"]},
    "relative_error_jr" -> N[errR/normR, 17],
    "relative_error_jw" -> N[errW/normW, 17],
    "max_relative_error" -> N[Max[errR/normR, errW/normW], 17],
    "matches_step8c_frechet_formula" -> TrueQ[Max[errR/normR, errW/normW] < 10^-7]
  |>
];

selfEnergySolve[sys_Association, jVec_List, mode_String] := Module[
  {op, solution, residual, sigma, boundaryQuadratic, flux},
  op = If[mode == "retarded", sys["ARet"], sys["ACons"]];
  solution = LinearSolve[op, mu0 jVec];
  residual = Norm[op . solution - mu0 jVec]/Max[Norm[mu0 jVec], 10^-30];
  sigma = Conjugate[jVec].(sys["WLane"].solution);
  boundaryQuadratic = Re[Conjugate[solution].(sys["BoundaryB"].solution)];
  flux = gammaPort sys["omega"]^5 boundaryQuadratic/(2.0 mu0);
  <|
    "solution" -> solution,
    "sigma" -> N[sigma, 17],
    "residual" -> N[residual, 17],
    "boundaryQuadratic" -> N[boundaryQuadratic, 17],
    "outgoingFlux" -> N[flux, 17]
  |>
];

fitEvenCoefficients[xs_List, ys_List] := Module[{mat, coeff},
  mat = N[Transpose[{ConstantArray[1.0, Length[xs]], xs, xs^2}], 30];
  coeff = LeastSquares[mat, N[ys, 30]];
  N[coeff, 17]
];

evenFitConditionNumber[window_?NumericQ] := Module[{fractions = {0.55, 0.70, 0.85, 1.00, 1.15}, xVals, mat, sv},
  xVals = N[(window fractions)^2, 30];
  mat = N[Transpose[{ConstantArray[1.0, Length[xVals]], xVals, xVals^2}], 30];
  sv = SingularValueList[mat];
  N[Max[sv]/Min[sv], 17]
];

transferForGrid[nr_Integer, nw_Integer, window_?NumericQ] := Module[
  {
    fractions = {0.55, 0.70, 0.85, 1.00, 1.15}, omegas, rows, omega,
    sys, matter, cons, ret, xVals, zVals, nVals, zCoeffs, nCoeffs,
    maxConsResidual, maxRetResidual, maxAmpResidual, maxPhaseResidual,
    maxContinuityResidual, fluxAtMax, boundaryAtMax, sourceNormAtMin
  },
  omegas = N[window fractions, 17];
  rows = Table[
    sys = buildMaxwellSystem[nr, nw, omega];
    matter = buildMatterSource[sys, omega];
    cons = selfEnergySolve[sys, matter["jVec"], "conservative"];
    ret = selfEnergySolve[sys, matter["jVec"], "retarded"];
    <|
      "omega" -> omega,
      "sigmaCons" -> cons["sigma"],
      "sigmaRet" -> ret["sigma"],
      "NofOmega" -> N[-Im[ret["sigma"]]/(gammaPort omega^5), 17],
      "consResidual" -> cons["residual"],
      "retResidual" -> ret["residual"],
      "ampResidual" -> matter["ampResidual"],
      "phaseResidual" -> matter["phaseResidual"],
      "continuityResidual" -> matter["continuityResidual"],
      "outgoingFlux" -> ret["outgoingFlux"],
      "boundaryQuadratic" -> ret["boundaryQuadratic"],
      "sourceNorm" -> matter["sourceNorm"]
    |>,
    {omega, omegas}
  ];
  xVals = omegas^2;
  zVals = Re[#["sigmaCons"]] & /@ rows;
  nVals = #["NofOmega"] & /@ rows;
  zCoeffs = fitEvenCoefficients[xVals, zVals];
  nCoeffs = fitEvenCoefficients[xVals, nVals];
  maxConsResidual = Max[#["consResidual"] & /@ rows];
  maxRetResidual = Max[#["retResidual"] & /@ rows];
  maxAmpResidual = Max[#["ampResidual"] & /@ rows];
  maxPhaseResidual = Max[#["phaseResidual"] & /@ rows];
  maxContinuityResidual = Max[#["continuityResidual"] & /@ rows];
  fluxAtMax = Last[rows]["outgoingFlux"];
  boundaryAtMax = Last[rows]["boundaryQuadratic"];
  sourceNormAtMin = First[rows]["sourceNorm"];
  <|
    "grid" -> {nr, nw},
    "h" -> gridData[nr, nw]["h"],
    "window" -> window,
    "omegas" -> omegas,
    "Z0" -> zCoeffs[[1]], "Z2" -> zCoeffs[[2]], "Z4" -> zCoeffs[[3]],
    "N0" -> nCoeffs[[1]], "N2" -> nCoeffs[[2]], "N4" -> nCoeffs[[3]],
    "rows" -> rows,
    "maxGreenResidual" -> Max[maxConsResidual, maxRetResidual],
    "maxBdgAmpResidual" -> maxAmpResidual,
    "maxBdgPhaseResidual" -> maxPhaseResidual,
    "maxContinuityResidual" -> maxContinuityResidual,
    "outgoingFluxAtMaxOmega" -> fluxAtMax,
    "boundaryQuadraticAtMaxOmega" -> boundaryAtMax,
    "sourceNormAtMinOmega" -> sourceNormAtMin
  |>
];

estimateOrderTriplet[hs_List, vals_List] := Module[{d12, d23},
  d12 = Abs[vals[[1]] - vals[[2]]];
  d23 = Abs[vals[[2]] - vals[[3]]];
  If[d12 <= 0 || d23 <= 0,
    Indeterminate,
    N[Log[d12/d23]/Log[hs[[1]]/hs[[2]]], 17]
  ]
];

fieldNormFromComponents[sys_Association, vec_List] := Module[
  {fm = sys["FieldMatrices"], w0, lw0, terms, compVals},
  w0 = sys["W0"];
  lw0 = Lharm w0;
  terms = {
    {fm["E_r"], w0}, {fm["E_E"], lw0}, {fm["E_B"], lw0}, {fm["E_w"], w0},
    {fm["C_r"], w0}, {fm["C_E"], lw0}, {fm["C_B"], lw0},
    {fm["B_B"], lw0}, {fm["B_r"], w0}, {fm["B_E"], lw0}
  };
  compVals = Table[
    With[{component = terms[[i, 1]] . vec, weight = terms[[i, 2]]},
      Conjugate[component] . (weight . component)
    ],
    {i, 1, Length[terms]}
  ];
  Sqrt[Max[0.0, Re[Total[compVals]]]]
];

makeGaugeVector[sys_Association, chi_List] := Module[
  {omega = sys["omega"], drMat = sys["Dr"], dwMat = sys["Dw"], rVals = sys["r"], n = sys["n"]},
  Join[
    I omega chi,
    drMat . chi,
    chi/rVals,
    ConstantArray[0.0, n],
    dwMat . chi
  ]
];

pureGaugeCheck[] := Module[
  {sys, chi, gaugeVec, source, solution, fieldNorm, flux, sourceNorm},
  sys = buildMaxwellSystem[6, 5, 0.043];
  chi = N[
    MapThread[
      Sin[Pi ((#1 - rMinTransfer)/(rExit - rMinTransfer))]^2 Sin[Pi #2/wLength]^2 &,
      {sys["r"], sys["w"]}
    ],
    17
  ];
  gaugeVec = makeGaugeVector[sys, chi];
  source = sys["ARet"] . gaugeVec;
  solution = LinearSolve[sys["ARet"], source];
  fieldNorm = fieldNormFromComponents[sys, solution];
  flux = gammaPort sys["omega"]^5 Re[Conjugate[solution].(sys["BoundaryB"].solution)]/(2.0 mu0);
  sourceNorm = Sqrt[Re[Conjugate[source].(sys["WLane"].source)]];
  <|
    "grid" -> {sys["nr"], sys["nw"]},
    "omega" -> sys["omega"],
    "source_norm" -> N[sourceNorm, 17],
    "physical_field_norm" -> N[fieldNorm, 17],
    "outgoing_flux" -> N[flux, 17],
    "physical_transfer_zero" -> TrueQ[fieldNorm < 10^-10 && Abs[flux] < 10^-18]
  |>
];

basisInvarianceCheck[] := Module[
  {
    sys, matter, jVec, omega, cons, ret, sigmaBase, nBase, wSqrt, wInvSqrt,
    aTildeCons, sTilde, bTilde, dim, deltas, seed = 240617, qRows,
    q, sigRot, solRot, nRot, k, uBase
  },
  omega = 0.039;
  sys = buildMaxwellSystem[4, 4, omega];
  matter = buildMatterSource[sys, omega];
  jVec = matter["jVec"];
  cons = selfEnergySolve[sys, jVec, "conservative"];
  sigmaBase = cons["sigma"];
  wSqrt = diagSparse[Sqrt[sys["WLaneValues"]]];
  wInvSqrt = diagSparse[1.0/Sqrt[sys["WLaneValues"]]];
  aTildeCons = Normal[wSqrt . sys["ACons"] . wInvSqrt];
  bTilde = Normal[wInvSqrt . sys["BoundaryB"] . wInvSqrt];
  sTilde = Sqrt[sys["WLaneValues"]] jVec;
  uBase = LinearSolve[aTildeCons, sTilde];
  nBase = Re[Conjugate[uBase] . bTilde . uBase];
  dim = Length[sTilde];
  SeedRandom[seed];
  deltas = Table[
    qRows = Orthogonalize[RandomReal[NormalDistribution[0, 1], {dim, dim}]];
    q = N[qRows, 17];
    solRot = LinearSolve[q . aTildeCons . Transpose[q], q . sTilde];
    sigRot = Conjugate[q . sTilde] . solRot;
    nRot = Re[Conjugate[solRot] . (q . bTilde . Transpose[q]) . solRot];
    <|
      "rotation" -> k,
      "relative_Z_delta" -> N[Abs[sigRot - sigmaBase]/Max[Abs[sigmaBase], 10^-30], 17],
      "relative_N_delta" -> N[Abs[nRot - nBase]/Max[Abs[nBase], 10^-30], 17]
    |>,
    {k, 1, 3}
  ];
  <|
    "grid" -> {sys["nr"], sys["nw"]},
    "omega" -> omega,
    "seed" -> seed,
    "base_Z" -> Re[sigmaBase],
    "base_N" -> N[nBase, 17],
    "rotations" -> deltas,
    "max_relative_Z_delta" -> Max[#["relative_Z_delta"] & /@ deltas],
    "max_relative_N_delta" -> Max[#["relative_N_delta"] & /@ deltas],
    "pass" -> TrueQ[
      Max[#["relative_Z_delta"] & /@ deltas] < 10^-10 &&
      Max[#["relative_N_delta"] & /@ deltas] < 10^-10
    ]
  |>
];

v209Regression[] := Module[
  {
    omegaU = 3.25, omegaW = 4.35, rMix = 0.08, gU = 0.18, gW = 0.13,
    aBlock, bBlock, kInt, g, ellVec, inv, delta, qNum, zDirect,
    zAnalytic, nDirect, nAnalytic
  },
  aBlock = omegaU^2;
  bBlock = omegaW^2;
  kInt = {{aBlock, -rMix}, {-rMix, bBlock}};
  g = {gU, gW};
  ellVec = {0.0, 1.0};
  inv = Inverse[kInt];
  delta = aBlock bBlock - rMix^2;
  qNum = gU^2 bBlock + 2.0 gU gW rMix + gW^2 aBlock;
  zDirect = g . inv . g;
  zAnalytic = qNum/delta;
  nDirect = (ellVec . inv . g)^2;
  nAnalytic = (aBlock gW + rMix gU)^2/delta^2;
  <|
    "Omega_U" -> omegaU,
    "Omega_W" -> omegaW,
    "R" -> rMix,
    "g_U" -> gU,
    "g_W" -> gW,
    "Delta" -> delta,
    "Z_direct" -> zDirect,
    "Z_Q_over_Delta" -> zAnalytic,
    "N0_direct" -> nDirect,
    "N0_square_over_Delta_squared" -> nAnalytic,
    "relative_Z_error" -> Abs[zDirect - zAnalytic]/Max[Abs[zAnalytic], 10^-30],
    "relative_N0_error" -> Abs[nDirect - nAnalytic]/Max[Abs[nAnalytic], 10^-30],
    "pass" -> TrueQ[
      Abs[zDirect - zAnalytic]/Max[Abs[zAnalytic], 10^-30] < 10^-13 &&
      Abs[nDirect - nAnalytic]/Max[Abs[nAnalytic], 10^-30] < 10^-13
    ]
  |>
];

hardCapContrast[omega_?NumericQ] := Module[
  {sys, matter, open, hardPenalty = 10^6, hardOp, hardSolution, openNorm, hardNorm, ratio},
  sys = buildMaxwellSystem[7, 7, omega];
  matter = buildMatterSource[sys, omega];
  open = selfEnergySolve[sys, matter["jVec"], "retarded"];
  hardOp = sys["WLaneInv"] . (sys["KCons"] + hardPenalty sys["BoundaryB"]);
  hardSolution = Quiet[LinearSolve[hardOp, mu0 matter["jVec"]]];
  openNorm = Sqrt[Max[0.0, Re[Conjugate[open["solution"]].(sys["BoundaryB"].open["solution"])]]];
  hardNorm = Sqrt[Max[0.0, Re[Conjugate[hardSolution].(sys["BoundaryB"].hardSolution)]]];
  ratio = hardNorm/Max[openNorm, 10^-30];
  <|
    "grid" -> {sys["nr"], sys["nw"]},
    "omega" -> omega,
    "open_boundary_norm" -> N[openNorm, 17],
    "hard_cap_boundary_norm" -> N[hardNorm, 17],
    "hard_cap_to_open_boundary_norm_ratio" -> N[ratio, 17],
    "open_not_hard_cap" -> TrueQ[openNorm > 10^-10 && ratio < 10^-4],
    "outgoing_flux" -> open["outgoingFlux"],
    "outgoing_flux_positive" -> TrueQ[open["outgoingFlux"] > 0]
  |>
];

formatScientific[x_] := StringReplace[ToString[FortranForm[N[x, 8]]], " " -> ""];
formatFixed[x_] := ToString[NumberForm[N[x, 8], {10, 6}]];

markdownTable[headers_List, rows_List] := StringRiffle[
  Join[
    {"| " <> StringRiffle[headers, " | "] <> " |"},
    {"| " <> StringRiffle[ConstantArray["---", Length[headers]], " | "] <> " |"},
    ("| " <> StringRiffle[ToString /@ #, " | "] <> " |") & /@ rows
  ],
  "\n"
];

Print["MT15-04 Spike-2 open DtN + forward j_eta + Z0/N0 transfer"];
Print["frozen Y_out,l=2: ", yOutBranchDescription];
Print["Gamma_port derived from compact sigma_Q^can: ", formatScientific[gammaPort]];

currentCheck = currentFrechetCheck[];
Print["current Frechet max relative error: ", formatScientific[currentCheck["max_relative_error"]]];

meshWindow = 0.040;
meshGrids = {{7, 7}, {9, 9}, {11, 11}};
meshTransfers = Table[
  Print["transfer mesh grid ", grid[[1]], "x", grid[[2]], " window ", meshWindow];
  transferForGrid[grid[[1]], grid[[2]], meshWindow],
  {grid, meshGrids}
];

z0DiagnosticGrid = {13, 13};
Print["diagnostic-only Z0 mesh grid ", z0DiagnosticGrid[[1]], "x", z0DiagnosticGrid[[2]], " window ", meshWindow];
z0DiagnosticTransfer = transferForGrid[z0DiagnosticGrid[[1]], z0DiagnosticGrid[[2]], meshWindow];

omegaWindows = {0.070, 0.0525, 0.040};
windowTransfers = Table[
  Print["omega-window transfer grid 9x9 window ", window];
  transferForGrid[9, 9, window],
  {window, omegaWindows}
];

primaryTransfer = Last[meshTransfers];

meshHs = #["h"] & /@ meshTransfers;
meshZ0s = #["Z0"] & /@ meshTransfers;
meshN0s = #["N0"] & /@ meshTransfers;
meshZOrders = <|
  "Z0" -> estimateOrderTriplet[meshHs, meshZ0s],
  "N0" -> estimateOrderTriplet[meshHs, meshN0s]
|>;

z0ExtendedTransfers = Append[meshTransfers, z0DiagnosticTransfer];
z0ExtendedHs = #["h"] & /@ z0ExtendedTransfers;
z0ExtendedZ0s = #["Z0"] & /@ z0ExtendedTransfers;
z0ExtendedN0s = #["N0"] & /@ z0ExtendedTransfers;
z0SuccessiveDiffs = N[Abs[Differences[z0ExtendedZ0s]], 17];
n0SuccessiveDiffs = N[Abs[Differences[z0ExtendedN0s]], 17];
z0ExtendedOrder91113 = estimateOrderTriplet[Take[z0ExtendedHs, -3], Take[z0ExtendedZ0s, -3]];
n0ExtendedOrder91113 = estimateOrderTriplet[Take[z0ExtendedHs, -3], Take[z0ExtendedN0s, -3]];
z0ResidualFloor = Max[#["maxGreenResidual"] & /@ z0ExtendedTransfers];
z0FitCondition = evenFitConditionNumber[meshWindow];
z0ConditionedResidualScale = z0FitCondition z0ResidualFloor Max[Abs[z0ExtendedZ0s]];
z0MinDiffToResidualRatio = Min[z0SuccessiveDiffs]/Max[z0ResidualFloor, 10^-300];
z0MinDiffToConditionedResidualRatio = Min[z0SuccessiveDiffs]/Max[z0ConditionedResidualScale, 10^-300];
z0OrderClassification = If[
  z0MinDiffToConditionedResidualRatio < 10^3,
  "roundoff_or_fit_floor_limited",
  "resolved_preasymptotic_conservative_spatial_boundary_discretization"
];
z0OrderNote =
  If[z0OrderClassification == "roundoff_or_fit_floor_limited",
    "The Z0 mesh differences are close to the Green residual floor after the omega-fit conditioning reference, so the sub-2 order should be treated as near-floor on this tiny M1b smoke Z0 rather than as a resolved spatial rate.",
    "A diagnostic-only 13x13 grid leaves the reported 7x7/9x9/11x11 rows unchanged and gives a 9x9/11x11/13x13 Z0 order of " <>
      formatScientific[z0ExtendedOrder91113] <> ". The Z0 successive differences are " <>
      StringRiffle[formatScientific /@ z0SuccessiveDiffs, ", "] <> ", while the largest Green-solve residual is " <>
      formatScientific[z0ResidualFloor] <> " and the conditioned residual-scale reference is " <>
      formatScientific[z0ConditionedResidualScale] <> ". These gaps rule out round-off or omega-fit floor control. The added point sharpens Z0 from the 7x7/9x9/11x11 order but does not make the original sub-2 estimate disappear into noise, so this is a resolved coarse/pre-asymptotic conservative spatial/boundary-discretization effect on the tiny M1b smoke Z0 rather than a persistent proof of sub-2 asymptotic order. Since Z0 is fitted from Re Sigma_cons before the imaginary open-DtN loss is added, the open-port DtN is not the cause; plausible contributors are the near-r_min 1/r^2 coefficients and boundary closures. Carry this forward for M1c because Z0 precision feeds D0/R_norm."
  ];

windowHs = omegaWindows;
windowZ0s = #["Z0"] & /@ windowTransfers;
windowN0s = #["N0"] & /@ windowTransfers;
windowOrders = <|
  "Z0" -> estimateOrderTriplet[windowHs, windowZ0s],
  "N0" -> estimateOrderTriplet[windowHs, windowN0s]
|>;

pureGauge = pureGaugeCheck[];
basisInvariant = basisInvarianceCheck[];
v209 = v209Regression[];
capContrast = hardCapContrast[0.040];

canonicalChiQCheck = <|
  "reference" -> "Sigma_ref^ret(omega)=-i sigma_Q^can omega^5 with N0_ref=1",
  "Gamma_port_used" -> gammaPort,
  "legacy_v2_09_residue_reduced_coefficient" -> gammaPortLegacyV209,
  "factor_vs_legacy" -> gammaPortFactorVsLegacy,
  "chi_Q_reference" -> 1.0,
  "fit_not_used" -> True,
  "note" -> "chi_Q=1 is an independent reference check after deriving Gamma_port from Poynting-positive compact sigma_Q^can; it is not used to set Gamma_port."
|>;

passChecks = <|
  "current_frechet_matches_step8c" -> currentCheck["matches_step8c_frechet_formula"],
  "outgoing_flux_positive" -> capContrast["outgoing_flux_positive"],
  "open_not_hard_cap" -> capContrast["open_not_hard_cap"],
  "pure_gauge_zero_physical_transfer" -> pureGauge["physical_transfer_zero"],
  "basis_invariance" -> basisInvariant["pass"],
  "v2_09_regression" -> v209["pass"],
  "green_residuals_small" -> TrueQ[primaryTransfer["maxGreenResidual"] < 10^-8],
  "bdg_residuals_small" -> TrueQ[primaryTransfer["maxBdgAmpResidual"] < 10^-8 && primaryTransfer["maxBdgPhaseResidual"] < 10^-8],
  "N0_positive" -> TrueQ[primaryTransfer["N0"] > 0]
|>;

identityChecks = <|
  "gamma_port_derived_not_fitted" -> <|
    "status" -> True,
    "note" -> "Construction identity: Gamma_port is computed solely from c_s and a through sigma_Q^can; this records code structure plus the transliteration-fidelity audit, not a runtime physics failure mode."
  |>,
  "gamma_port_uses_compact_sigma_Q_can" -> <|
    "status" -> TrueQ[Abs[gammaPort - gammaPortCompact] < 10^-15],
    "note" -> "Construction identity: deterministic equality Gamma_port=4 a^5/(27 c_s^5)."
  |>,
  "legacy_factor_four_recorded" -> <|
    "status" -> TrueQ[Abs[gammaPortFactorVsLegacy - 4.0] < 10^-12],
    "note" -> "Construction identity: deterministic factor-four ledger against the V2-09 residue-reduced coefficient."
  |>
|>;

termMap = {
  <|"term" -> "open exit/radial DtN", "implementation" -> "K_ret=K_cons+i Gamma_port omega^5 B_tan, with B_tan assembled from E_r/E_E/E_B traces at w=L and E_E/E_B/E_w traces at r=R_exit", "source" -> "Decision 05 D2; scratch log D2 lines 9617-9648"|>,
  <|"term" -> "outgoing flux", "implementation" -> "Phi_out=(Gamma_port omega^5/2mu0) A^dagger B_tan A >= 0, equivalent to D2 Poynting form under frozen DtN branch", "source" -> "scratch log D2 physical flux lines 9635-9644"|>,
  <|"term" -> "wall drive", "implementation" -> "delta V_conf=-4 V_radial r^4/R0^5 chi_eta(w)", "source" -> "M1b mt15_02_bdg_wall_derivation.wls driveCoeff/delta V ledger"|>,
  <|"term" -> "BdG response", "implementation" -> "M1b smoke L_plus/L_phase response on the same VSH transfer grid; no self-consistent M1c claim", "source" -> "mt15_02_bdg_wall_derivation.wls finite-dimensional BdG sector"|>,
  <|"term" -> "phasor current", "implementation" -> "Step-8c Frechet formula with delta A_i=0 for the external j_eta source and A_i0=0 on the M1b smoke background", "source" -> "step8c report lines 31-32; scratch log D3 lines 9652-9687"|>,
  <|"term" -> "S_eta^(A)", "implementation" -> "not used; forward source only", "source" -> "Decision 04 and scratch log D3 lines 9684-9687"|>,
  <|"term" -> "self-energy extraction", "implementation" -> "Sigma_A=<j_eta,A>, Z_n from Re Sigma_cons, N_n from -Im Sigma_ret/(Gamma_port omega^5)", "source" -> "Decision 05 D4; scratch log D4 lines 9690-9722"|>,
  <|"term" -> "Gamma_port", "implementation" -> "Gamma_port=sigma_Q^can=9/(8 Omega_Q^5)=4a^5/(27c_s^5), not fitted to chi_Q", "source" -> "compact lines 2559-2568; notes lines 41400-41415"|>
};

meshRows = Table[
  {
    ToString[t["grid"][[1]]] <> "x" <> ToString[t["grid"][[2]]],
    formatFixed[t["h"]],
    formatScientific[t["Z0"]],
    formatScientific[t["Z2"]],
    formatScientific[t["Z4"]],
    formatScientific[t["N0"]],
    formatScientific[t["N2"]],
    formatScientific[t["N4"]],
    formatScientific[t["maxGreenResidual"]]
  },
  {t, meshTransfers}
];

windowRows = Table[
  {
    formatFixed[t["window"]],
    formatScientific[t["Z0"]],
    formatScientific[t["N0"]],
    formatScientific[t["maxGreenResidual"]],
    formatScientific[t["outgoingFluxAtMaxOmega"]]
  },
  {t, windowTransfers}
];

gateRows = ({#, If[passChecks[#], "PASS", "FAIL"]} & /@ Keys[passChecks]);
identityRows = ({#, If[identityChecks[#]["status"], "TRUE", "FALSE"], identityChecks[#]["note"]} & /@ Keys[identityChecks]);
termRows = ({#["term"], #["implementation"], #["source"]} & /@ termMap);

diagnostics = <|
  "schema" -> "mt15_04_spike2_transfer_n0/v1",
  "scope" -> "Spike-2 only: open DtN, forward j_eta source, Z0/N0 extraction. No V2-22B packet and no R_norm.",
  "source_commit" -> sourceCommitString[],
  "target_blind" -> True,
  "machinery_not_physics" -> "Z0/N0 below are machinery/method-demonstration numbers on the non-self-consistent M1b smoke background; physical values require M1c.",
  "z0_order_note" -> z0OrderNote,
  "open_dtn" -> <|
    "boundary_class" -> "open_impedance",
    "Y_out_l2" -> yOutBranchDescription,
    "mouth" -> "no injected Maxwell drive; scalar gauge anchor only",
    "radial_outer" -> "same frozen outgoing tangential DtN matrix on a_E/a_B/a_w",
    "Gamma_port" -> toJsonNumber[gammaPort],
    "Gamma_port_formula" -> "sigma_Q^can=9/(8 Omega_Q^5)=4 a^5/(27 c_s^5)",
    "legacy_v2_09_reduced_formula" -> "a^5/(27 c_s^5)",
    "legacy_factor_gap" -> toJsonNumber[gammaPortFactorVsLegacy],
    "derived_not_fitted_to_chi_Q" -> True
  |>,
  "forward_source" -> <|
    "chain" -> "eta -> delta V_conf -> M1b BdG response -> step8c Frechet delta J_psi -> Maxwell Green solve",
    "delta_V_conf" -> "delta V_conf=-4 V_radial r^4/R0^5 eta",
    "current_formula" -> "delta rho=2(psi_R0 delta psi_R+psi_I0 delta psi_I); delta j_i=(hbar/m)(...)-(q/m)(delta A_i rho0+A_i0 delta rho)",
    "delta_A_in_external_source" -> "set to zero while forming j_eta; Maxwell back-reaction is the solved field, not an added source",
    "S_eta_A_used" -> False,
    "current_frechet_check" -> currentCheck
  |>,
  "coefficients" -> <|
    "primary_grid" -> primaryTransfer["grid"],
    "primary_window" -> toJsonNumber[primaryTransfer["window"]],
    "Z0" -> toJsonNumber[primaryTransfer["Z0"]],
    "Z2" -> toJsonNumber[primaryTransfer["Z2"]],
    "Z4" -> toJsonNumber[primaryTransfer["Z4"]],
    "N0" -> toJsonNumber[primaryTransfer["N0"]],
    "N2" -> toJsonNumber[primaryTransfer["N2"]],
    "N4" -> toJsonNumber[primaryTransfer["N4"]],
    "rows" -> (Map[toJsonNumber, #, {2}] & /@ ({#["omega"], Re[#["sigmaCons"]], Im[#["sigmaRet"]], #["NofOmega"], #["outgoingFlux"]} & /@ primaryTransfer["rows"]))
  |>,
  "mesh_convergence" -> <|
    "window" -> toJsonNumber[meshWindow],
    "rows" -> (AssociationThread[
        {"grid", "h", "Z0", "Z2", "Z4", "N0", "N2", "N4", "maxGreenResidual"},
        {#["grid"], toJsonNumber[#["h"]], toJsonNumber[#["Z0"]], toJsonNumber[#["Z2"]], toJsonNumber[#["Z4"]], toJsonNumber[#["N0"]], toJsonNumber[#["N2"]], toJsonNumber[#["N4"]], toJsonNumber[#["maxGreenResidual"]]}
      ] & /@ meshTransfers),
    "observed_orders_from_successive_differences" -> Map[toJsonNumber, meshZOrders],
    "z0_order_diagnostic" -> <|
      "classification" -> z0OrderClassification,
      "z0_order_note" -> z0OrderNote,
      "diagnostic_probe_row" -> AssociationThread[
        {"grid", "h", "Z0", "Z2", "Z4", "N0", "N2", "N4", "maxGreenResidual"},
        {z0DiagnosticTransfer["grid"], toJsonNumber[z0DiagnosticTransfer["h"]], toJsonNumber[z0DiagnosticTransfer["Z0"]], toJsonNumber[z0DiagnosticTransfer["Z2"]], toJsonNumber[z0DiagnosticTransfer["Z4"]], toJsonNumber[z0DiagnosticTransfer["N0"]], toJsonNumber[z0DiagnosticTransfer["N2"]], toJsonNumber[z0DiagnosticTransfer["N4"]], toJsonNumber[z0DiagnosticTransfer["maxGreenResidual"]]}
      ],
      "z0_successive_differences_7_9_11_13" -> toJsonVector[z0SuccessiveDiffs],
      "n0_successive_differences_7_9_11_13" -> toJsonVector[n0SuccessiveDiffs],
      "z0_order_9_11_13" -> toJsonNumber[z0ExtendedOrder91113],
      "n0_order_9_11_13" -> toJsonNumber[n0ExtendedOrder91113],
      "max_green_residual_floor_reference" -> toJsonNumber[z0ResidualFloor],
      "omega_fit_condition_number" -> toJsonNumber[z0FitCondition],
      "conditioned_residual_scale_reference" -> toJsonNumber[z0ConditionedResidualScale],
      "min_z0_diff_to_green_residual_ratio" -> toJsonNumber[z0MinDiffToResidualRatio],
      "min_z0_diff_to_conditioned_residual_ratio" -> toJsonNumber[z0MinDiffToConditionedResidualRatio],
      "open_dtn_causality_note" -> "Z0 is fitted from Re Sigma_cons using the conservative operator before the imaginary open-DtN boundary loss is added; the open-port DtN can affect N0/flux but is not the direct source of the Z0 mesh order."
    |>
  |>,
  "omega_window_convergence" -> <|
    "grid" -> {9, 9},
    "windows" -> toJsonVector[omegaWindows],
    "rows" -> (AssociationThread[
        {"window", "Z0", "N0", "maxGreenResidual", "outgoingFluxAtMaxOmega"},
        {toJsonNumber[#["window"]], toJsonNumber[#["Z0"]], toJsonNumber[#["N0"]], toJsonNumber[#["maxGreenResidual"]], toJsonNumber[#["outgoingFluxAtMaxOmega"]]}
      ] & /@ windowTransfers),
    "observed_orders_from_successive_window_differences" -> Map[toJsonNumber, windowOrders]
  |>,
  "normalization" -> <|
    "Gamma_port_derivation" -> "Poynting-positive compact retarded denominator uses sigma_Q^can=9/(8 Omega_Q^5). With Omega_Q=3 c_s/(2a), Gamma_port=4 a^5/(27 c_s^5). V2-09 a^5/(27 c_s^5) is the c1=1/4 residue-reduced wall coefficient.",
    "canonical_chi_Q_reference_check" -> canonicalChiQCheck
  |>,
  "validation" -> <|
    "pure_gauge" -> pureGauge,
    "basis_invariance" -> basisInvariant,
    "v2_09_regression" -> v209,
    "hard_cap_contrast" -> capContrast,
    "pass_checks" -> passChecks,
    "identity_checks" -> identityChecks
  |>,
  "term_map" -> termMap
|>;

Export[diagnosticsPath, diagnostics, "JSON"];

meshTable = markdownTable[
  {"grid", "h", "Z0", "Z2", "Z4", "N0", "N2", "N4", "max residual"},
  meshRows
];
windowTable = markdownTable[
  {"omega window", "Z0", "N0", "max residual", "positive flux"},
  windowRows
];
gateTable = markdownTable[{"check", "status"}, gateRows];
identityTable = markdownTable[{"identity check", "value", "note"}, identityRows];
termTable = markdownTable[{"term", "implementation", "source"}, termRows];

report = StringRiffle[
  {
    "# MT15-04 Spike-2 Transfer N0 Report",
    "",
    "Status: Spike-2 only. No V2-22B packet, no `R_norm`, no V2-22B extension, and no `S_eta^(A)` return source.",
    "",
    "**MACHINERY LABEL:** the coefficients below are method-demonstration numbers on the non-self-consistent M1b smoke background. They are not physical branch values; M1c is required for that.",
    "",
    "## Frozen Open DtN",
    "",
    "`boundary_class=open_impedance`. The frozen outgoing branch is `" <> yOutBranchDescription <> "`.",
    "",
    "Implementation: `K_ret=K_cons+i Gamma_port omega^5 B_tan`, where `B_tan` is assembled from the physical tangential electric-field traces `E_r/E_E/E_B` at `w=L` and `E_E/E_B/E_w` at `r=R_exit`. The mouth has no injected Maxwell drive; only a scalar gauge anchor is added. This is an open impedance/DtN loss term, not a hard cap.",
    "",
    "## Gamma Port",
    "",
    "The compact retarded quadrupole denominator gives `sigma_Q^can=9/(8 Omega_Q^5)` with `Omega_Q=3 c_s/(2a)`, hence `Gamma_port=4 a^5/(27 c_s^5)`. The older V2-09 `a^5/(27 c_s^5)` is the residue-reduced wall coefficient after the minimal isotropic pole weight `c1=1/4`; it is not the physical outgoing-port denominator normalization used for Spike-2 extraction.",
    "",
    "- `Gamma_port`: `" <> formatScientific[gammaPort] <> "`.",
    "- legacy V2-09 reduced coefficient: `" <> formatScientific[gammaPortLegacyV209] <> "`.",
    "- factor gap recorded: `" <> formatFixed[gammaPortFactorVsLegacy] <> "`.",
    "- chi_Q reference check: `1.0`, computed after the derivation on `Sigma_ref=-i sigma_Q^can omega^5`; it was not used as a fit knob.",
    "",
    "The Gamma-port value identities are construction checks, not runtime physics gates. `Gamma_port` being derived rather than fitted is established by code structure (`Gamma_port` depends only on `c_s` and `a`) together with the transliteration-fidelity audit; the runtime identity table only records the deterministic formula/value ledger.",
    "",
    "## Forward Source",
    "",
    "The chain is `eta -> delta V_conf -> delta psi -> delta J_psi`: the wall mode is the M1b weak-form wall eigen-shape, `delta V_conf=-4 V_radial r^4/R0^5 eta`, the M1b smoke BdG amplitude/phase response is solved on the transfer grid, and the Step-8c phasor Frechet current forms `j_eta`. The external source uses `delta A_i=0` and `A_i0=0`; the solved Maxwell field is the response, not an extra source. `S_eta^(A)` is not used.",
    "",
    "Current Frechet check against a central finite difference of the nonlinear gauge-covariant current: max relative error `" <> formatScientific[currentCheck["max_relative_error"]] <> "`.",
    "",
    "## Coefficients",
    "",
    meshTable,
    "",
    "Primary reported values are the finest mesh row (`" <> ToString[primaryTransfer["grid"][[1]]] <> "x" <> ToString[primaryTransfer["grid"][[2]]] <> "`, omega window `" <> formatFixed[primaryTransfer["window"]] <> "`):",
    "",
    "- `Z0,Z2,Z4 = " <> StringRiffle[formatScientific /@ {primaryTransfer["Z0"], primaryTransfer["Z2"], primaryTransfer["Z4"]}, ", "] <> "`.",
    "- `N0,N2,N4 = " <> StringRiffle[formatScientific /@ {primaryTransfer["N0"], primaryTransfer["N2"], primaryTransfer["N4"]}, ", "] <> "`.",
    "",
    "Mesh convergence orders from successive-difference stabilization: `Z0 -> " <> formatScientific[meshZOrders["Z0"]] <> "`, `N0 -> " <> formatScientific[meshZOrders["N0"]] <> "`.",
    "",
    "Z0 order diagnostic: " <> z0OrderNote,
    "",
    "Omega-window refinement on the `9x9` grid:",
    "",
    windowTable,
    "",
    "Omega-window orders from successive-difference stabilization: `Z0 -> " <> formatScientific[windowOrders["Z0"]] <> "`, `N0 -> " <> formatScientific[windowOrders["N0"]] <> "`.",
    "",
    "## Validation",
    "",
    "Genuine pass checks:",
    "",
    gateTable,
    "",
    "Construction identity checks:",
    "",
    identityTable,
    "",
    "- Basis invariance: max relative `Z` delta `" <> formatScientific[basisInvariant["max_relative_Z_delta"]] <> "`, max relative `N` delta `" <> formatScientific[basisInvariant["max_relative_N_delta"]] <> "`.",
    "- V2-09 regression: relative `Z` error `" <> formatScientific[v209["relative_Z_error"]] <> "`, relative `N0` error `" <> formatScientific[v209["relative_N0_error"]] <> "`.",
    "- Outgoing flux: `" <> formatScientific[capContrast["outgoing_flux"]] <> "`; hard-cap/open boundary norm ratio `" <> formatScientific[capContrast["hard_cap_to_open_boundary_norm_ratio"]] <> "`.",
    "- Pure gauge: physical field norm `" <> formatScientific[pureGauge["physical_field_norm"]] <> "`, outgoing flux `" <> formatScientific[pureGauge["outgoing_flux"]] <> "`.",
    "",
    "## Term Map",
    "",
    termTable,
    "",
    "Diagnostics JSON: `" <> diagnosticsPath <> "`."
  },
  "\n"
];

Export[reportPath, report, "Text"];

Print["diagnostics_path: ", diagnosticsPath];
Print["report_path: ", reportPath];
Print["machinery_not_physics: coefficients are M1b smoke machinery, not physical branch values"];
Print["Gamma_port: ", formatScientific[gammaPort], " factor_vs_legacy: ", formatFixed[gammaPortFactorVsLegacy]];
Print["primary Z0/Z2/Z4: ", StringRiffle[formatScientific /@ {primaryTransfer["Z0"], primaryTransfer["Z2"], primaryTransfer["Z4"]}, ", "]];
Print["primary N0/N2/N4: ", StringRiffle[formatScientific /@ {primaryTransfer["N0"], primaryTransfer["N2"], primaryTransfer["N4"]}, ", "]];
Print["mesh orders Z0/N0: ", formatScientific[meshZOrders["Z0"]], " / ", formatScientific[meshZOrders["N0"]]];
Print["diagnostic 13x13 Z0/N0: ", formatScientific[z0DiagnosticTransfer["Z0"]], " / ", formatScientific[z0DiagnosticTransfer["N0"]]];
Print["diagnostic Z0 order 9/11/13: ", formatScientific[z0ExtendedOrder91113], " classification: ", z0OrderClassification];
Print["omega-window orders Z0/N0: ", formatScientific[windowOrders["Z0"]], " / ", formatScientific[windowOrders["N0"]]];
Print["basis invariance max rel Z/N: ", formatScientific[basisInvariant["max_relative_Z_delta"]], " / ", formatScientific[basisInvariant["max_relative_N_delta"]]];
Print["V2-09 rel Z/N errors: ", formatScientific[v209["relative_Z_error"]], " / ", formatScientific[v209["relative_N0_error"]]];
Print["outgoing_flux_positive: ", capContrast["outgoing_flux_positive"], " open_not_hard_cap: ", capContrast["open_not_hard_cap"]];
Print["pure_gauge_zero_physical_transfer: ", pureGauge["physical_transfer_zero"]];
Print["pass_checks: ", passChecks];
Print["identity_checks: ", identityChecks];
