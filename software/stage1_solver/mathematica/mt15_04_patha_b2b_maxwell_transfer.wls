(* Path-A B2b Spike-2 Maxwell transfer worker.

   Dispatched by mt15_04_spike2_transfer_n0.wls --patha-b2b.
   It reads the B2a closed background and BdG bundle, assembles the l=2
   five-lane VSH Maxwell Green operator in Wolfram, and emits
   {Z0,Z2,Z4,N0,N2,N4}. No D0/P0/R_norm assembly is performed here.
*)

ClearAll["Global`pathaB2B`*"];

args = $ScriptCommandLine;

argValue[flag_String, default_] := Module[{pos},
  pos = FirstPosition[args, flag, Missing["NotFound"]];
  If[MissingQ[pos] || pos[[1]] >= Length[args], default, args[[pos[[1]] + 1]]]
];

toJsonNumber[x_] := Module[{n = Quiet[Check[N[x], $Failed]]},
  If[n === $Failed || ! NumericQ[n],
    Print["ERROR: nonnumeric JSON value: ", InputForm[x]];
    Exit[4]
  ];
  If[Abs[Im[n]] < 10^-12, N[Re[n], 17], <|"re" -> N[Re[n], 17], "im" -> N[Im[n], 17]|>]
];
toJsonVector[x_List] := toJsonNumber /@ x;

clipToGrid[x_, xs_List] := Min[Max[x, First[xs]], Last[xs]];
interp1Clamped[xs_List, ys_List, x_] := Module[{xc = clipToGrid[x, xs], k, t},
  If[xc <= First[xs], Return[First[ys]]];
  If[xc >= Last[xs], Return[Last[ys]]];
  k = First[FirstPosition[UnitStep[xs - xc], 1]] - 1;
  If[k < 1, k = 1];
  If[k >= Length[xs], k = Length[xs] - 1];
  t = (xc - xs[[k]])/(xs[[k + 1]] - xs[[k]]);
  (1 - t) ys[[k]] + t ys[[k + 1]]
];
interp2Clamped[rOld_List, wOld_List, values_List, r_, w_] := Module[
  {wc = clipToGrid[w, wOld], byR},
  byR = Table[interp1Clamped[wOld, values[[i]], wc], {i, 1, Length[rOld]}];
  interp1Clamped[rOld, byR, r]
];

diag[v_List] := DiagonalMatrix[N[v, 17]];

firstDerivativeMatrix[n_Integer, h_?NumericQ] := Module[{mat = ConstantArray[0.0, {n, n}]},
  Do[
    If[i < n, mat[[i, i + 1]] += 1.0/(2.0 h)];
    If[i > 1, mat[[i, i - 1]] -= 1.0/(2.0 h)];
    ,
    {i, 1, n}
  ];
  N[mat, 17]
];

secondDerivativeMatrix[n_Integer, h_?NumericQ] := Module[{mat = ConstantArray[0.0, {n, n}]},
  Do[
    mat[[i, i]] -= 2.0/h^2;
    If[i < n, mat[[i, i + 1]] += 1.0/h^2];
    If[i > 1, mat[[i, i - 1]] += 1.0/h^2];
    ,
    {i, 1, n}
  ];
  N[mat, 17]
];

blockRow[blocks_List] := ArrayFlatten[{blocks}];
weightedGram[m_, weights_List, sign_] := sign ConjugateTranspose[m].diag[weights].m;
cellIndex[i_, j_, nr_] := (j - 1) nr + i;

backgroundPath = argValue["--background", ""];
bdgPath = argValue["--bdg", ""];
outPath = argValue["--out", ""];
nr = ToExpression[argValue["--nr", "7"]];
nw = ToExpression[argValue["--nw", "7"]];
window = ToExpression[argValue["--window", "0.036"]];
radialScale = ToExpression[argValue["--radial-scale", "1.0"]];
modeCount = ToExpression[argValue["--mode-count", "30"]];

If[backgroundPath == "" || bdgPath == "" || outPath == "",
  Print["ERROR: --background, --bdg, and --out are required for --patha-b2b"];
  Exit[2]
];

background = Import[backgroundPath, "RawJSON"];
bdgPacket = Import[bdgPath, "RawJSON"];

geometry = background["geometry"];
constants = background["constants"];
aConst = N[geometry["a"], 17];
Lw = N[geometry["L"], 17];
ell = 2;
Lharm = 6.0;
xiGauge = 1.0;
csConst = 1.0;
gammaPort = N[4.0 aConst^5/(27.0 csConst^5), 17];
gaugeCharge = N[constants["gauge_charge"], 17];
particleMass = N[constants["particle_mass"], 17];
mu0 = N[constants["mu0"], 17];

rOld = N[background["grid"]["r_centers"], 17];
wOld = N[background["grid"]["w_centers"], 17];
psiROld = N[background["fields"]["psi_R0"], 17];
psiIOld = N[background["fields"]["psi_I0"], 17];
a0Old = N[background["fields"]["A_00"], 17];
arOld = N[background["fields"]["A_r0"], 17];
awOld = N[background["fields"]["A_w0"], 17];
zOld = N[background["derived"]["Z_w"], 17];
r0Old = N[background["fields"]["R0_w"], 17];

gridData[nr_Integer, nw_Integer, radialScale_?NumericQ] := Module[
  {
    rOuter, rMin, dr, dw, rVals1, wVals1, rVals, wVals, zVals, r0Vals,
    quadWeights, psiR2, psiI2, a02, ar2, aw2
  },
  rOuter = N[background["grid"]["r_max"] radialScale, 17];
  rMin = N[0.22 aConst, 17];
  dr = (rOuter - rMin)/(nr + 1);
  dw = Lw/(nw + 1);
  rVals1 = N[Table[rMin + i dr, {i, 1, nr}], 17];
  wVals1 = N[Table[i dw, {i, 1, nw}], 17];
  rVals = N[Flatten[Table[rVals1[[i]], {j, 1, nw}, {i, 1, nr}]], 17];
  wVals = N[Flatten[Table[wVals1[[j]], {j, 1, nw}, {i, 1, nr}]], 17];
  zVals = N[interp1Clamped[wOld, zOld, #] & /@ wVals, 17];
  r0Vals = N[interp1Clamped[wOld, r0Old, #] & /@ wVals, 17];
  quadWeights = N[dr dw rVals^2 zVals, 17];
  psiR2 = N[Table[interp2Clamped[rOld, wOld, psiROld, rVals1[[i]], wVals1[[j]]], {i, 1, nr}, {j, 1, nw}], 17];
  psiI2 = N[Table[interp2Clamped[rOld, wOld, psiIOld, rVals1[[i]], wVals1[[j]]], {i, 1, nr}, {j, 1, nw}], 17];
  a02 = N[Table[interp2Clamped[rOld, wOld, a0Old, rVals1[[i]], wVals1[[j]]], {i, 1, nr}, {j, 1, nw}], 17];
  ar2 = N[Table[interp2Clamped[rOld, wOld, arOld, rVals1[[i]], wVals1[[j]]], {i, 1, nr}, {j, 1, nw}], 17];
  aw2 = N[Table[interp2Clamped[rOld, wOld, awOld, rVals1[[i]], wVals1[[j]]], {i, 1, nr}, {j, 1, nw}], 17];
  <|
    "nr" -> nr, "nw" -> nw, "n" -> nr nw,
    "dr" -> dr, "dw" -> dw, "r1" -> rVals1, "w1" -> wVals1,
    "r" -> rVals, "w" -> wVals, "z" -> zVals, "R0" -> r0Vals,
    "w0" -> quadWeights, "h" -> Max[dr, dw],
    "r_min" -> rMin, "r_outer" -> rOuter, "radial_scale" -> radialScale,
    "psiR" -> Flatten[Table[psiR2[[i, j]], {j, 1, nw}, {i, 1, nr}]],
    "psiI" -> Flatten[Table[psiI2[[i, j]], {j, 1, nw}, {i, 1, nr}]],
    "A00" -> Flatten[Table[a02[[i, j]], {j, 1, nw}, {i, 1, nr}]],
    "Ar0" -> Flatten[Table[ar2[[i, j]], {j, 1, nw}, {i, 1, nr}]],
    "Aw0" -> Flatten[Table[aw2[[i, j]], {j, 1, nw}, {i, 1, nr}]]
  |>
];

boundaryMatrix[grid_Association, fm_Association] := Module[
  {nr = grid["nr"], nw = grid["nw"], n = grid["n"], wExitScalar, wExitAngular,
   rOuterScalar, rOuterAngular, cell, wExitWeight, rExitWeight},
  wExitScalar = ConstantArray[0.0, n];
  wExitAngular = ConstantArray[0.0, n];
  rOuterScalar = ConstantArray[0.0, n];
  rOuterAngular = ConstantArray[0.0, n];
  Do[
    cell = cellIndex[i, j, nr];
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
  ConjugateTranspose[fm["E_r"]].diag[wExitScalar].fm["E_r"] +
    ConjugateTranspose[fm["E_E"]].diag[wExitAngular + rOuterAngular].fm["E_E"] +
    ConjugateTranspose[fm["E_B"]].diag[wExitAngular + rOuterAngular].fm["E_B"] +
    ConjugateTranspose[fm["E_w"]].diag[rOuterScalar].fm["E_w"]
];

buildMaxwellSystem[nr_Integer, nw_Integer, omega_?NumericQ, radialScale_?NumericQ] := Module[
  {grid, n, total, drMat1, dwMat1, drrMat1, dwwMat1, idR, idW, idN, zN,
   drMat, dwMat, drrMat, dwwMat, rVals, invR, invR2, rDiag, r2Diag,
   invRDiag, invR2Diag, w0, wLaneVals, wLane, wLaneInv, divR, mG,
   mEr, mEE, mEB, mEw, mCr, mCE, mCB, mBB, mBr, mBE, kAction,
   boundaryB, anchorVals, anchorCell, kCons, kRet, fields},
  grid = gridData[nr, nw, radialScale];
  n = grid["n"];
  total = 5 n;
  drMat1 = firstDerivativeMatrix[nr, grid["dr"]];
  dwMat1 = firstDerivativeMatrix[nw, grid["dw"]];
  drrMat1 = secondDerivativeMatrix[nr, grid["dr"]];
  dwwMat1 = secondDerivativeMatrix[nw, grid["dw"]];
  idR = IdentityMatrix[nr];
  idW = IdentityMatrix[nw];
  idN = IdentityMatrix[n];
  zN = ConstantArray[0.0, {n, n}];
  drMat = KroneckerProduct[idW, drMat1];
  dwMat = KroneckerProduct[dwMat1, idR];
  drrMat = KroneckerProduct[idW, drrMat1];
  dwwMat = KroneckerProduct[dwwMat1, idR];
  rVals = grid["r"];
  invR = N[1.0/rVals, 17];
  invR2 = N[invR^2, 17];
  rDiag = diag[rVals];
  r2Diag = diag[rVals^2];
  invRDiag = diag[invR];
  invR2Diag = diag[invR2];
  w0 = grid["w0"];
  wLaneVals = Join[w0, w0, Lharm w0, Lharm w0, w0];
  wLane = diag[wLaneVals];
  wLaneInv = diag[1.0/wLaneVals];
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
  fields = <|"G" -> mG, "E_r" -> mEr, "E_E" -> mEE, "E_B" -> mEB, "E_w" -> mEw,
    "C_r" -> mCr, "C_E" -> mCE, "C_B" -> mCB, "B_B" -> mBB, "B_r" -> mBr, "B_E" -> mBE|>;
  kAction =
    weightedGram[mEr, w0, 1.0] +
    weightedGram[mEE, Lharm w0, 1.0] +
    weightedGram[mEB, Lharm w0, 1.0] +
    weightedGram[mEw, w0, 1.0] -
    weightedGram[mCr, w0, 1.0] -
    weightedGram[mCE, Lharm w0, 1.0] -
    weightedGram[mCB, Lharm w0, 1.0] -
    weightedGram[mBB, Lharm w0, 1.0] -
    weightedGram[mBr, w0, 1.0] -
    weightedGram[mBE, Lharm w0, 1.0] -
    (1.0/xiGauge) weightedGram[mG, w0, 1.0];
  boundaryB = boundaryMatrix[grid, fields];
  anchorVals = ConstantArray[0.0, total];
  anchorCell = Ceiling[n/2];
  anchorVals[[anchorCell]] = 10^-8 w0[[anchorCell]];
  kCons = kAction + diag[anchorVals];
  kRet = kCons + I gammaPort omega^5 boundaryB;
  Join[
    grid,
    <|"omega" -> omega, "total" -> total, "Dr" -> drMat, "Dw" -> dwMat,
      "Drr" -> drrMat, "Dww" -> dwwMat, "WLane" -> wLane,
      "WLaneValues" -> wLaneVals, "WLaneInv" -> wLaneInv,
      "KCons" -> kCons, "KRet" -> kRet, "ACons" -> wLaneInv . kCons,
      "ARet" -> wLaneInv . kRet, "BoundaryB" -> boundaryB, "FieldMatrices" -> fields|>
  ]
];

bdgProfileValue[mode_Association, x_] := Module[{points, profile},
  points = N[bdgPacket["grid"]["profile_points"], 17];
  profile = N[If[KeyExistsQ[mode, "profile_values"], mode["profile_values"], mode["profile"]], 17];
  interp1Clamped[points, profile, x]
];

bdgResponseOnW[wVals_List, omega_?NumericQ] := Module[{modes, response},
  modes = Take[bdgPacket["bdg_modes"], UpTo[modeCount]];
  response = Table[
    Total[
      Table[
        With[{varpi = N[mode["varpi"], 17], coupling = N[mode["coupling"], 17]},
          coupling bdgProfileValue[mode, w]/(varpi^2 - omega^2)
        ],
        {mode, modes}
      ]
    ],
    {w, wVals}
  ];
  N[response, 17]
];

buildMatterSource[sys_Association] := Module[
  {n = sys["n"], nr = sys["nr"], nw = sys["nw"], omega = sys["omega"], rho, axial,
   radialShape, deltaRho, cells, raw, weight, norm, meanA0, scalarMod, jr, jE, jw,
   zero, jVec, sourceNorm, continuity, continuityResidual},
  rho = N[sys["psiR"]^2 + sys["psiI"]^2, 17];
  axial = bdgResponseOnW[sys["w1"], omega];
  radialShape = ConstantArray[0.0, n];
  Do[
    cells = Table[cellIndex[i, j, nr], {i, 1, nr}];
    raw = rho[[cells]];
    weight = sys["dr"] sys["r"][[cells]]^2 sys["z"][[cells]];
    norm = Total[weight raw];
    If[! NumericQ[norm] || norm <= 10^-300,
      raw = ConstantArray[1.0, Length[cells]];
      norm = Total[weight raw];
    ];
    Do[radialShape[[cells[[k]]]] = raw[[k]]/norm, {k, 1, Length[cells]}];
    ,
    {j, 1, nw}
  ];
  meanA0 = Total[sys["w0"] sys["A00"]]/Total[sys["w0"]];
  scalarMod = N[1.0 + gaugeCharge (sys["A00"] - meanA0), 17];
  deltaRho = ConstantArray[0.0, n];
  Do[
    cells = Table[cellIndex[i, j, nr], {i, 1, nr}];
    Do[deltaRho[[cells[[k]]]] = axial[[j]] radialShape[[cells[[k]]]] scalarMod[[cells[[k]]]], {k, 1, Length[cells]}];
    ,
    {j, 1, nw}
  ];
  jr = N[-(gaugeCharge/particleMass) sys["Ar0"] deltaRho, 17];
  jE = N[-I omega sys["r"] deltaRho/Lharm, 17];
  jw = N[-(gaugeCharge/particleMass) sys["Aw0"] deltaRho, 17];
  zero = ConstantArray[0.0, n];
  jVec = N[Join[gaugeCharge deltaRho, gaugeCharge jr, gaugeCharge jE, zero, gaugeCharge jw], 17];
  sourceNorm = Sqrt[Max[0.0, Re[Conjugate[jVec].(sys["WLane"].jVec)]]];
  continuity = N[-I omega deltaRho + sys["Dr"].jr + sys["Dw"].jw - Lharm jE/sys["r"], 17];
  continuityResidual = Norm[continuity]/Max[Norm[deltaRho] + Norm[jr] + Norm[jw] + Norm[jE], 10^-30];
  <|"jVec" -> jVec, "sourceNorm" -> N[sourceNorm, 17],
    "continuityResidual" -> N[continuityResidual, 17],
    "usedB2aBdgModes" -> Min[modeCount, Length[bdgPacket["bdg_modes"]]],
    "usedReturnSourceSEtaA" -> False|>
];

selfEnergySolve[sys_Association, jVec_List, mode_String] := Module[
  {op, solution, residual, sigma, boundaryQuadratic, flux},
  op = If[mode == "retarded", sys["ARet"], sys["ACons"]];
  solution = LinearSolve[op, mu0 jVec];
  residual = Norm[op . solution - mu0 jVec]/Max[Norm[mu0 jVec], 10^-30];
  sigma = Conjugate[jVec].(sys["WLane"].solution);
  boundaryQuadratic = Re[Conjugate[solution].(sys["BoundaryB"].solution)];
  flux = gammaPort sys["omega"]^5 boundaryQuadratic/(2.0 mu0);
  <|"solution" -> solution, "sigma" -> N[sigma, 17], "residual" -> N[residual, 17],
    "boundaryQuadratic" -> N[boundaryQuadratic, 17], "outgoingFlux" -> N[flux, 17]|>
];

fitEvenCoefficients[omegas_List, ys_List] := Module[{xVals, mat, coeff},
  xVals = N[omegas^2, 30];
  mat = N[Transpose[{ConstantArray[1.0, Length[xVals]], xVals, xVals^2}], 30];
  coeff = LeastSquares[mat, N[ys, 30]];
  N[coeff, 17]
];

evenFitConditionNumber[omegas_List] := Module[{xVals, mat, sv},
  xVals = N[omegas^2, 30];
  mat = N[Transpose[{ConstantArray[1.0, Length[xVals]], xVals, xVals^2}], 30];
  sv = SingularValueList[mat];
  N[Max[sv]/Min[sv], 17]
];

transferForGrid[nr_Integer, nw_Integer, window_?NumericQ, radialScale_?NumericQ] := Module[
  {fractions = {0.55, 0.70, 0.85, 1.00, 1.15}, omegas, rows, sys, matter,
   cons, ret, zCoeffs, nCoeffs},
  omegas = N[window fractions, 17];
  rows = Table[
    sys = buildMaxwellSystem[nr, nw, omega, radialScale];
    matter = buildMatterSource[sys];
    cons = selfEnergySolve[sys, matter["jVec"], "conservative"];
    ret = selfEnergySolve[sys, matter["jVec"], "retarded"];
    <|
      "omega" -> toJsonNumber[omega],
      "sigmaCons" -> toJsonNumber[cons["sigma"]],
      "sigmaRet" -> toJsonNumber[ret["sigma"]],
      "minusImSigmaRet" -> toJsonNumber[-Im[ret["sigma"]]],
      "NofOmega" -> toJsonNumber[-Im[ret["sigma"]]/(gammaPort omega^5)],
      "consResidual" -> toJsonNumber[cons["residual"]],
      "retResidual" -> toJsonNumber[ret["residual"]],
      "continuityResidual" -> toJsonNumber[matter["continuityResidual"]],
      "outgoingFlux" -> toJsonNumber[ret["outgoingFlux"]],
      "boundaryQuadratic" -> toJsonNumber[ret["boundaryQuadratic"]],
      "sourceNorm" -> toJsonNumber[matter["sourceNorm"]]
    |>,
    {omega, omegas}
  ];
  zCoeffs = fitEvenCoefficients[omegas, Re[#["sigmaCons"]] & /@ rows];
  nCoeffs = fitEvenCoefficients[omegas, #["NofOmega"] & /@ rows];
  <|
    "schema" -> "stage1_patha_b2b_maxwell_transfer/v1",
    "engine" -> "mathematica_mt15_04_patha_b2b",
    "background_content_hash" -> background["content_hash"],
    "bdg_packet_content_hash" -> bdgPacket["content_hash"],
    "tau" -> toJsonNumber[constants["tau"]],
    "geometry" -> <|"a" -> toJsonNumber[aConst], "L" -> toJsonNumber[Lw], "radial_scale" -> toJsonNumber[radialScale]|>,
    "constants" -> <|
      "gauge_charge" -> toJsonNumber[gaugeCharge],
      "particle_mass" -> toJsonNumber[particleMass],
      "mu0" -> toJsonNumber[mu0],
      "tau" -> toJsonNumber[constants["tau"]],
      "Gamma_port" -> toJsonNumber[gammaPort],
      "Gamma_port_formula" -> "4*a^5/(27*c_s^5) with c_s=1"
    |>,
    "grid" -> <|"nr" -> nr, "nw" -> nw, "window" -> toJsonNumber[window],
      "radial_scale" -> toJsonNumber[radialScale],
      "mode_count" -> Min[modeCount, Length[bdgPacket["bdg_modes"]]]|>,
    "coefficients" -> <|"Z0" -> toJsonNumber[zCoeffs[[1]]], "Z2" -> toJsonNumber[zCoeffs[[2]]],
      "Z4" -> toJsonNumber[zCoeffs[[3]]], "N0" -> toJsonNumber[nCoeffs[[1]]],
      "N2" -> toJsonNumber[nCoeffs[[2]]], "N4" -> toJsonNumber[nCoeffs[[3]]]|>,
    "rows" -> rows,
    "diagnostics" -> <|
      "operator_matrix_dimension" -> 5 nr nw,
      "max_green_residual" -> toJsonNumber[Max[Max[#["consResidual"], #["retResidual"]] & /@ rows]],
      "max_continuity_residual" -> toJsonNumber[Max[#["continuityResidual"] & /@ rows]],
      "min_minus_im_sigma_ret" -> toJsonNumber[Min[#["minusImSigmaRet"] & /@ rows]],
      "min_outgoing_flux" -> toJsonNumber[Min[#["outgoingFlux"] & /@ rows]],
      "fit_condition_number_Z" -> toJsonNumber[evenFitConditionNumber[omegas]],
      "fit_condition_number_N" -> toJsonNumber[evenFitConditionNumber[omegas]],
      "source_chain" -> "decision-05 D3 Frechet source over B2a BdG-mode response -> axial density response -> Path-A A0-modulated charge/current source -> five-lane VSH Maxwell Green solve",
      "operator_term_map" -> {
        <|"term" -> "five-lane VSH Maxwell weak form", "implementation" -> "lanes phi,a_r,a_E,a_B,a_w with L=l(l+1)=6 and H=Z gauge weight", "source" -> "decision 05 D1-D2"|>,
        <|"term" -> "open DtN loss", "implementation" -> "K_ret=K_cons+i Gamma_port omega^5 B_tan", "source" -> "compact Gamma_port lines 2559-2568"|>,
        <|"term" -> "basis-invariant self energy", "implementation" -> "Sigma=<j,G j>; Z_n from Re Sigma_cons; N_n from -Im Sigma_ret/(Gamma_port omega^5)", "source" -> "decision 05 D4"|>,
        <|"term" -> "wall/BdG source", "implementation" -> "canonical decision-05 D3 Frechet source driven by the B2a BdG-mode response; B2a overlaps c_j carry the chi coupling, no live U/W ports", "source" -> "pathA_11 B2b requirements 3-4"|>
      }
    |>
  |>
];

Print["Path-A B2b Maxwell transfer worker"];
Print["background_path: ", backgroundPath];
Print["bdg_path: ", bdgPath];
Print["grid/window/radial_scale: ", nr, "x", nw, " / ", window, " / ", radialScale];

packet = transferForGrid[nr, nw, window, radialScale];
packetJson = ExportString[packet, "JSON", "Compact" -> True];
packetWithHash = Append[packet, "content_hash" -> Hash[packetJson, "SHA256", "HexString"]];

outDir = DirectoryName[outPath];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
Export[outPath, packetWithHash, "JSON"];

Print["packet_path: ", outPath];
Print["coefficients: ", packet["coefficients"]];
Print["max_green_residual: ", packet["diagnostics"]["max_green_residual"]];
