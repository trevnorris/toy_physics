(* Path-A B2a BdG derivation worker, dispatched by mt15_02_bdg_wall_derivation.wls. *)

ClearAll["Global`pathaB2A`*"];

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
  N[Re[n], 17]
];
toJsonVector[x_List] := toJsonNumber /@ x;
relativeChange[new_, old_] := Abs[new - old]/Max[Abs[new], 10^-300];
momentAssoc[count_Integer] := Module[{m = Min[count, Length[couplingList]], b0, b2, b4},
  b0 = N[Total[Table[couplingList[[i]]^2/varpiList[[i]]^2, {i, 1, m}]], 17];
  b2 = N[Total[Table[couplingList[[i]]^2/varpiList[[i]]^4, {i, 1, m}]], 17];
  b4 = N[Total[Table[couplingList[[i]]^2/varpiList[[i]]^6, {i, 1, m}]], 17];
  <|
    "B0" -> toJsonNumber[b0],
    "B2" -> toJsonNumber[b2],
    "B4" -> toJsonNumber[b4],
    "B0_B4_minus_B2_squared" -> toJsonNumber[b0 b4 - b2^2]
  |>
];
momentRelErrors[current_Association, reference_Association] := <|
  "B0" -> toJsonNumber[relativeChange[current["B0"], reference["B0"]]],
  "B2" -> toJsonNumber[relativeChange[current["B2"], reference["B2"]]],
  "B4" -> toJsonNumber[relativeChange[current["B4"], reference["B4"]]]
|>;
maxRelError[assoc_Association] := Max[Values[assoc]];

trapz[xs_List, ys_List] /; Length[xs] == Length[ys] && Length[xs] >= 2 := Total[
  Table[0.5 (ys[[i]] + ys[[i + 1]]) (xs[[i + 1]] - xs[[i]]), {i, 1, Length[xs] - 1}]
];

weightedNorm[xs_List, weights_List, f_List] := trapz[xs, weights f f];

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

rotateByMaxComponent[z_List] := Module[{k, phase},
  k = First[Ordering[Abs[z], -1]];
  phase = If[Abs[z[[k]]] > 0, Exp[-I Arg[z[[k]]]], 1.0];
  phase z
];

cellIndex[i_, j_, nr_] := (j - 1) nr + i;

backgroundPath = argValue["--background", ""];
wallInputPath = argValue["--wall-input", ""];
outPath = argValue["--out", ""];
nr = ToExpression[argValue["--nr", "10"]];
nw = ToExpression[argValue["--nw", "10"]];
modesToExport = ToExpression[argValue["--modes", "30"]];

If[backgroundPath == "" || wallInputPath == "" || outPath == "",
  Print["ERROR: --background, --wall-input, and --out are required for --patha-b2a"];
  Exit[2]
];

background = Import[backgroundPath, "RawJSON"];
wallInput = Import[wallInputPath, "RawJSON"];

constants = background["constants"];
geometry = background["geometry"];
hbar = N[constants["hbar"], 17];
particleMass = N[constants["particle_mass"], 17];
gaugeCharge = N[constants["gauge_charge"], 17];
eosK = N[constants["eos_K"], 17];
tau = N[constants["tau"], 17];
radialWallStrength = N[constants["radial_wall_strength"], 17];
axialTrapStrength = N[constants["axial_trap_strength"], 17];
chemicalPotential = N[background["solver"]["chemical_potential"], 17];
L = N[geometry["L"], 17];
a = N[geometry["a"], 17];
ellFactor = 6.0;

rOld = N[background["grid"]["r_centers"], 17];
wOld = N[background["grid"]["w_centers"], 17];
psiROld = N[background["fields"]["psi_R0"], 17];
psiIOld = N[background["fields"]["psi_I0"], 17];
a0Old = N[background["fields"]["A_00"], 17];
arOld = N[background["fields"]["A_r0"], 17];
awOld = N[background["fields"]["A_w0"], 17];
r0Old = N[background["fields"]["R0_w"], 17];

rMax = N[background["grid"]["r_max"], 17];
rFaces = N[Table[rMax (i - 1)/nr, {i, 1, nr + 1}], 17];
wFaces = N[Table[L (j - 1)/nw, {j, 1, nw + 1}], 17];
dr = rFaces[[2]] - rFaces[[1]];
dw = wFaces[[2]] - wFaces[[1]];
rCenters = N[Table[0.5 (rFaces[[i]] + rFaces[[i + 1]]), {i, 1, nr}], 17];
wCenters = N[Table[0.5 (wFaces[[j]] + wFaces[[j + 1]]), {j, 1, nw}], 17];
radialFaceAreas = N[4.0 Pi rFaces^2, 17];
radialVolumes = N[(4.0 Pi/3.0) Table[rFaces[[i + 1]]^3 - rFaces[[i]]^3, {i, 1, nr}], 17];

psiR = Table[interp2Clamped[rOld, wOld, psiROld, rCenters[[i]], wCenters[[j]]], {i, 1, nr}, {j, 1, nw}];
psiI = Table[interp2Clamped[rOld, wOld, psiIOld, rCenters[[i]], wCenters[[j]]], {i, 1, nr}, {j, 1, nw}];
a0Values2D = Table[interp2Clamped[rOld, wOld, a0Old, rCenters[[i]], wCenters[[j]]], {i, 1, nr}, {j, 1, nw}];
arValues2D = Table[interp2Clamped[rOld, wOld, arOld, rCenters[[i]], wCenters[[j]]], {i, 1, nr}, {j, 1, nw}];
awValues2D = Table[interp2Clamped[rOld, wOld, awOld, rCenters[[i]], wCenters[[j]]], {i, 1, nr}, {j, 1, nw}];
r0Values = Table[interp1Clamped[wOld, r0Old, wCenters[[j]]], {j, 1, nw}];

nCells = nr nw;
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
bdgR0 = Flatten[Table[r0Values[[j]], {j, 1, nw}, {i, 1, nr}]];
psiVec = Flatten[Table[psiR[[i, j]] + I psiI[[i, j]], {j, 1, nw}, {i, 1, nr}]];
rho0 = N[Abs[psiVec]^2, 17];
A0Values = Flatten[Table[a0Values2D[[i, j]], {j, 1, nw}, {i, 1, nr}]];
spatialGaugeNorm = Max[
  Max[Abs[Flatten[arValues2D]]],
  Max[Abs[Flatten[awValues2D]]]
];

Vconf = N[
  radialWallStrength (bdgRadii/bdgR0)^4 +
    0.5 axialTrapStrength ((bdgW - 0.5 L)/L)^2,
  17
];
driveCoeff = N[4.0 radialWallStrength bdgRadii^4/bdgR0^5, 17];
enthalpy = N[(5.0 eosK/4.0) rho0^4, 17];
dhdrho = N[5.0 eosK rho0^3, 17];
enthalpyLinearization = N[rho0 dhdrho, 17];

angularMatrix = DiagonalMatrix[N[ellFactor/bdgRadii^2, 17]];
kineticMatrix = N[-(hbar^2/(2.0 particleMass)) (laplacianMatrix - angularMatrix), 17];
diagBlock = N[
  kineticMatrix + DiagonalMatrix[Vconf + gaugeCharge A0Values - chemicalPotential + enthalpy + enthalpyLinearization],
  17
];
pairBlock = DiagonalMatrix[N[dhdrho psiVec^2, 17]];
bdgMatrix = ArrayFlatten[{
    {diagBlock, pairBlock},
    {-Conjugate[pairBlock], -Conjugate[diagBlock]}
  }];

Print["MT15 Path-A B2a BdG derivation export"];
Print["background_path: ", backgroundPath];
Print["wall_input_path: ", wallInputPath];
Print["building dense BdG eigensystem with dimension: ", 2 nCells];

{bdgEigenvalues, bdgEigenvectors} = Eigensystem[N[bdgMatrix, 30]];
bdgPairs = SortBy[Transpose[{bdgEigenvalues, bdgEigenvectors}], Re[First[#]] &];
positiveBdgPairs = Select[
  bdgPairs,
  Re[#[[1]]] > 10^-8 && Abs[Im[#[[1]]]] < 10^-7 &
];

If[Length[positiveBdgPairs] < modesToExport,
  Print["ERROR: fewer positive real BdG eigenmodes found than requested: ", Length[positiveBdgPairs]];
  Exit[3]
];

cellVolumeWeights = N[Flatten[Table[radialVolumes[[i]] dw, {j, 1, nw}, {i, 1, nr}]], 17];
points = N[wallInput["profile_points"], 17];
wallChiEta = N[wallInput["chi"], 17];
weightValues = N[wallInput["weight"], 17];

modeRecords = {};
couplingList = {};
varpiList = {};
bdgResidualList = {};

Do[
  eigenvalue = positiveBdgPairs[[modeIdx, 1]];
  rotatedVec = rotateByMaxComponent[positiveBdgPairs[[modeIdx, 2]]];
  u0 = rotatedVec[[1 ;; nCells]];
  v0 = rotatedVec[[nCells + 1 ;; 2 nCells]];
  symplecticNormBefore = Re[Total[cellVolumeWeights (Abs[u0]^2 - Abs[v0]^2)]];
  scale = 1.0/Sqrt[Max[Abs[symplecticNormBefore], 10^-300]];
  modeVec = scale rotatedVec;
  u = modeVec[[1 ;; nCells]];
  v = modeVec[[nCells + 1 ;; 2 nCells]];
  symplecticNormAfter = Re[Total[cellVolumeWeights (Abs[u]^2 - Abs[v]^2)]];

  residualAbs = Norm[bdgMatrix . modeVec - eigenvalue modeVec];
  residualRel = N[
    residualAbs/Max[Norm[bdgMatrix . modeVec], Abs[eigenvalue] Norm[modeVec], 1.0],
    17
  ];

  densityResponse = N[Re[Conjugate[psiVec] u + psiVec v], 17];
  driveResponse = N[driveCoeff densityResponse, 17];
  densityByW = Table[
    Total[Table[radialVolumes[[i]] driveResponse[[cellIndex[i, j, nr]]], {i, 1, nr}]],
    {j, 1, nw}
  ];
  densityOnPacketGrid = Table[interp1Clamped[wCenters, densityByW, points[[k]]], {k, 1, Length[points]}];
  signedPhysicalCoupling = -trapz[points, wallChiEta densityOnPacketGrid];
  cAlpha = Abs[signedPhysicalCoupling];
  varpi = N[Re[eigenvalue]/hbar, 17];

  AppendTo[couplingList, cAlpha];
  AppendTo[varpiList, varpi];
  If[modeIdx <= modesToExport,
    rawPhi = N[densityOnPacketGrid/weightValues, 17];
    phiNormValue = weightedNorm[points, weightValues, rawPhi];
    If[! NumericQ[phiNormValue] || phiNormValue <= 0,
      Print["ERROR: nonpositive profile norm for mode ", modeIdx, ": ", phiNormValue];
      Exit[5]
    ];
    phiNorm = N[rawPhi/Sqrt[phiNormValue], 17];
    adapterOverlap = trapz[points, weightValues wallChiEta phiNorm];
    If[adapterOverlap < 0,
      phiNorm = -phiNorm;
      adapterOverlap = -adapterOverlap;
    ];
    If[adapterOverlap <= 10^-300,
      Print["ERROR: near-zero adapter overlap for mode ", modeIdx];
      Exit[6]
    ];
    lambdaB = cAlpha/adapterOverlap;
    AppendTo[bdgResidualList, residualRel];
    AppendTo[modeRecords, <|
      "name" -> "patha_b2a_bdg_mode_" <> ToString[modeIdx],
      "lambda_B" -> toJsonNumber[lambdaB],
      "varpi" -> toJsonNumber[varpi],
      "profile" -> toJsonVector[phiNorm],
      "profile_values" -> toJsonVector[phiNorm],
      "overlap_I_eta_phi" -> toJsonNumber[adapterOverlap],
      "coupling" -> toJsonNumber[cAlpha],
      "wall_drive_overlap_signed" -> toJsonNumber[signedPhysicalCoupling],
      "raw_eigenvalue" -> toJsonNumber[Re[eigenvalue]],
      "eigenvalue_imag_abs" -> toJsonNumber[Abs[Im[eigenvalue]]],
      "relative_eigen_residual" -> toJsonNumber[residualRel],
      "symplectic_norm_before_abs_scaling" -> toJsonNumber[symplecticNormBefore],
      "symplectic_norm_after_abs_scaling" -> toJsonNumber[symplecticNormAfter],
      "profile_norm_weighted" -> toJsonNumber[weightedNorm[points, weightValues, phiNorm]]
    |>];
  ];
  ,
  {modeIdx, 1, Length[positiveBdgPairs]}
];

modalTruncationTolerance = 1.0*^-4;
allPositiveModeCount = Length[couplingList];
modalSweepCounts = Select[{3, 5, 8, 15, 30}, # < allPositiveModeCount &];
If[! MemberQ[modalSweepCounts, allPositiveModeCount], AppendTo[modalSweepCounts, allPositiveModeCount]];
allPositiveMoments = momentAssoc[allPositiveModeCount];
exportedMoments = momentAssoc[modesToExport];
selectedError = momentRelErrors[exportedMoments, allPositiveMoments];
modalSweepRows = {};
previousMoments = Null;
Do[
  currentMoments = momentAssoc[count];
  previousRel = If[previousMoments === Null, Null, momentRelErrors[currentMoments, previousMoments]];
  allRel = momentRelErrors[currentMoments, allPositiveMoments];
  AppendTo[modalSweepRows, <|
    "M" -> count,
    "label" -> If[count == allPositiveModeCount, "all_positive", ToString[count]],
    "B0" -> currentMoments["B0"],
    "B2" -> currentMoments["B2"],
    "B4" -> currentMoments["B4"],
    "B0_B4_minus_B2_squared" -> currentMoments["B0_B4_minus_B2_squared"],
    "rel_change_from_previous" -> previousRel,
    "max_rel_change_from_previous" -> If[previousRel === Null, Null, toJsonNumber[maxRelError[previousRel]]],
    "rel_error_vs_all_positive" -> allRel,
    "max_rel_error_vs_all_positive" -> toJsonNumber[maxRelError[allRel]]
  |>];
  previousMoments = currentMoments;
  ,
  {count, modalSweepCounts}
];
eligibleRows = Select[
  modalSweepRows,
  #["label"] =!= "all_positive" && #["max_rel_error_vs_all_positive"] <= modalTruncationTolerance &
];
eligibleCounts = If[Length[eligibleRows] > 0, eligibleRows[[All, "M"]], {}];
smallestCountMeetingTolerance = If[Length[eligibleCounts] > 0, Min[eligibleCounts], allPositiveModeCount];
modalTruncation = <|
  "sweep_counts" -> modalSweepCounts,
  "sweep" -> modalSweepRows,
  "all_positive_mode_count" -> allPositiveModeCount,
  "all_positive_moments" -> allPositiveMoments,
  "exported_mode_count" -> modesToExport,
  "exported_moments" -> exportedMoments,
  "truncation_error_at_export" -> selectedError,
  "max_truncation_error_at_export" -> toJsonNumber[maxRelError[selectedError]],
  "tolerance" -> toJsonNumber[modalTruncationTolerance],
  "smallest_count_meeting_tolerance" -> smallestCountMeetingTolerance,
  "selection_rule" -> "exported K must have max_n |B_n(K)-B_n(all)|/|B_n(all)| <= tolerance; all-positive is the reference tail check"
|>;
B0 = exportedMoments["B0"];
B2 = exportedMoments["B2"];
B4 = exportedMoments["B4"];
stieltjesGap = exportedMoments["B0_B4_minus_B2_squared"];
maxBdgResidual = Max[bdgResidualList];

packet = <|
  "schema" -> "stage1_patha_b2a_bdg_bundle/v1",
  "engine" -> "mathematica_mt15_02_patha_b2a",
  "background_content_hash" -> background["content_hash"],
  "wall_input_hash" -> wallInput["content_hash"],
  "tau" -> toJsonNumber[tau],
  "geometry" -> <|"a" -> toJsonNumber[a], "L" -> toJsonNumber[L]|>,
  "constants" -> <|
    "eos_K" -> toJsonNumber[eosK],
    "hbar" -> toJsonNumber[hbar],
    "particle_mass" -> toJsonNumber[particleMass],
    "gauge_charge" -> toJsonNumber[gaugeCharge],
    "tau" -> toJsonNumber[tau],
    "radial_wall_strength" -> toJsonNumber[radialWallStrength],
    "axial_trap_strength" -> toJsonNumber[axialTrapStrength]
  |>,
  "grid" -> <|
    "nr" -> nr,
    "nw" -> nw,
    "dr" -> toJsonNumber[dr],
    "dw" -> toJsonNumber[dw],
    "r_centers" -> toJsonVector[rCenters],
    "w_centers" -> toJsonVector[wCenters],
    "profile_points" -> toJsonVector[points]
  |>,
  "wall" -> <|
    "source" -> wallInput["source"],
    "K" -> toJsonNumber[wallInput["K"]],
    "M" -> toJsonNumber[wallInput["M"]],
    "chi" -> toJsonVector[wallChiEta],
    "weight" -> toJsonVector[weightValues],
    "normalization_chi_T_W_chi" -> toJsonNumber[wallInput["normalization_chi_T_W_chi"]],
    "trapz_weight_norm" -> toJsonNumber[wallInput["trapz_weight_norm"]],
    "orientation_integral" -> toJsonNumber[wallInput["orientation_integral"]]
  |>,
  "bdg_modes" -> modeRecords,
  "bdg_moments" -> <|
    "B0" -> toJsonNumber[B0],
    "B2" -> toJsonNumber[B2],
    "B4" -> toJsonNumber[B4],
    "B0_B4_minus_B2_squared" -> toJsonNumber[stieltjesGap]
  |>,
  "diagnostics" -> <|
    "operator_matrix_dimension" -> 2 nCells,
    "positive_real_modes_found" -> Length[positiveBdgPairs],
    "modes_exported" -> modesToExport,
    "max_relative_eigen_residual" -> toJsonNumber[maxBdgResidual],
    "modal_truncation" -> modalTruncation,
    "spatial_gauge_linf" -> toJsonNumber[spatialGaugeNorm],
    "operator_term_map" -> {
      <|"term" -> "kinetic", "implementation" -> "-(hbar^2/(2m))*(laplacian-l(l+1)/r^2), l(l+1)=6", "source" -> "notes/moving_throat_pde_program_compact.md:1406 and existing mt15_02 form"|>,
      <|"term" -> "EOS", "implementation" -> "h=(5K/4)rho^4, dh/drho=5K rho^3", "source" -> "notes/moving_throat_pde_program_compact.md:578"|>,
      <|"term" -> "wall drive", "implementation" -> "delta V_conf=-4 V_radial r^4/R0^5 eta", "source" -> "notes/moving_throat_pde_program_compact.md:1080"|>
    }
  |>
|>;

packetJson = ExportString[packet, "JSON", "Compact" -> True];
packetWithHash = Append[packet, "content_hash" -> Hash[packetJson, "SHA256", "HexString"]];

outDir = DirectoryName[outPath];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
Export[outPath, packetWithHash, "JSON"];

Print["packet_path: ", outPath];
Print["BdG positive modes found: ", Length[positiveBdgPairs], " exported: ", modesToExport];
Print["BdG varpi exported: ", ToString[NumberForm[Take[varpiList, modesToExport], {16, 10}]]];
Print["BdG max relative eigen-residual: ", ToString[ScientificForm[maxBdgResidual, 6]]];
Print["BdG moments B0/B2/B4: ", ToString[NumberForm[{B0, B2, B4}, {16, 10}]]];
Print["BdG modal truncation max rel at export: ", ToString[ScientificForm[modalTruncation["max_truncation_error_at_export"], 6]]];
Print["Stieltjes Cauchy-Schwarz identity B0*B4-B2^2: ", ToString[ScientificForm[stieltjesGap, 6]]];
