#!/usr/bin/env wolframscript

(* M1b BdG + wall derivation producer.
   This is pre-freeze, target-blind machinery.  It derives the BdG modes from
   a finite-dimensional eigensolve and derives the wall weak-form coefficients.
   The mixed Maxwell ports are deliberately posited placeholders. *)

$HistoryLength = 0;

ClearAll["Global`*"];

scriptPath =
  If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
    $InputFileName,
    FileNameJoin[{Directory[], "software/stage1_solver/mathematica/mt15_02_bdg_wall_derivation.wls"}]
  ];
scriptDir = DirectoryName[scriptPath];
repoRoot = AbsoluteFileName[FileNameJoin[{scriptDir, "..", "..", ".."}]];

If[MemberQ[$ScriptCommandLine, "--patha-b2a"],
  Get[FileNameJoin[{scriptDir, "mt15_02_patha_b2a_bdg_wall_derivation.wls"}]];
  Exit[0]
];

runDir = FileNameJoin[{scriptDir, "runs", "mt15_02_bdg_wall_derivation"}];
If[! DirectoryQ[runDir], CreateDirectory[runDir, CreateIntermediateDirectories -> True]];

packetPath = FileNameJoin[{runDir, "mt15_02_v2_22b_bdg_wall_packet.json"}];
diagnosticsPath = FileNameJoin[{runDir, "mt15_02_derivation_diagnostics.json"}];

ClearAll[
  trapz, weightedInner, weightedNorm, normalizeProfile, smootherStep, toJsonNumber,
  toJsonVector, toJsonMatrix, sourceCommitString, cellIndex, rotateByMaxComponent
];

trapz[xs_List, ys_List] /; Length[xs] == Length[ys] && Length[xs] >= 2 :=
  Total[
    Table[
      0.5 (ys[[i]] + ys[[i + 1]]) (xs[[i + 1]] - xs[[i]]),
      {i, 1, Length[xs] - 1}
    ]
  ];

weightedInner[xs_List, weights_List, f_List, g_List] :=
  trapz[xs, weights f g];

weightedNorm[xs_List, weights_List, f_List] :=
  weightedInner[xs, weights, f, f];

normalizeProfile[xs_List, weights_List, f_List, label_String] := Module[{norm = weightedNorm[xs, weights, f]},
  If[! NumericQ[norm] || norm <= 0,
    Print["ERROR: nonpositive norm while normalizing ", label, ": ", norm];
    Exit[2]
  ];
  N[f/Sqrt[norm], 17]
];

smootherStep[x_] := 3 x^2 - 2 x^3;
toJsonNumber[x_] := Module[{n = Quiet[Check[N[x], $Failed]]},
  If[n === $Failed || ! NumericQ[n],
    Print["ERROR: nonnumeric value reached JSON export: ", InputForm[x]];
    Exit[4]
  ];
  N[Re[n], 17]
];
toJsonVector[x_List] := toJsonNumber /@ x;
toJsonMatrix[x_List] := Map[toJsonNumber, x, {2}];

rotateByMaxComponent[z_List] := Module[{k, phase},
  k = First[Ordering[Abs[z], -1]];
  phase = If[Abs[z[[k]]] > 0, Exp[-I Arg[z[[k]]]], 1.0];
  phase z
];

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

(* Target-blind geometry and constitutive placeholders. *)
L = 2.0;
a = 1.0;
RExit = 1.65;
meshPoints = 129;

hbar = 1.0;
particleMass = 1.0;
gaugeCharge = 0.35;
eosK = 0.04;
chemicalPotential = 0.05;
ell = 2;
ellFactor = ell (ell + 1);
radialWallStrength = 0.65;
axialTrapStrength = 0.12;

GConst = 1.0;
csConst = 1.0;
cConst = 2.5;
mhat0Const = 0.85;
SPortConst = 1.2;
thetaTailConst = 1.0;

points = N[Table[L i/(meshPoints - 1), {i, 0, meshPoints - 1}], 17];
xGrid = points/L;
radiusValues = N[a + (RExit - a) (smootherStep /@ xGrid), 17];
weightValues = N[radiusValues^2, 17];

(* Wall weak-form matrices from V2-20 lines 118-136. *)
muEta = N[1.0 + 0.08 Cos[Pi xGrid]^2, 17];
Tw = N[1.1 + 0.10 Sin[Pi xGrid], 17];
TOmega = N[0.35 + 0.03 Cos[2 Pi xGrid], 17];
KEta = N[0.45 + 0.08 256 xGrid^4 (1 - xGrid)^4, 17];
Vl = N[KEta + ellFactor TOmega, 17];
YL0 = 0.0;

basisCount = 3;
basisValues = Table[N[xGrid^p, 17], {p, 1, basisCount}];
basisDerivatives = Table[
  N[If[p == 1, ConstantArray[1.0/L, meshPoints], (p/L) xGrid^(p - 1)], 17],
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
wallChiEta = N[wallShapeRaw/Sqrt[wallWeightNorm], 17];
wallMScalar = N[wallCoeff . wallMMatrix . wallCoeff, 17];
wallKScalar = N[wallCoeff . wallKMatrix . wallCoeff, 17];
wallMMinEigenvalue = Min[Eigenvalues[N[wallMMatrix, 30]]];

(* Finite-dimensional BdG operator.  The spatial gauge background is zero in
   this target-blind smoke branch, so the minimal spatial covariant terms are
   represented but vanish.  A0 is also zero; the q A0 term remains in the
   operator ledger as q*A0Values. *)
nr = 12;
nw = 15;
nCells = nr nw;
rMax = RExit;
dr = rMax/nr;
dw = L/(nw - 1);
rFaces = N[Table[dr i, {i, 0, nr}], 17];
rCenters = N[Table[(i - 0.5) dr, {i, 1, nr}], 17];
wCenters = N[Table[L (j - 1)/(nw - 1), {j, 1, nw}], 17];
wQuadWeights = N[Table[If[j == 1 || j == nw, dw/2, dw], {j, 1, nw}], 17];
radialVolumes = N[Table[(rFaces[[i + 1]]^3 - rFaces[[i]]^3)/3, {i, 1, nr}], 17];
radialFaceAreas = N[rFaces^2, 17];

cellIndex[i_, j_] := (j - 1) nr + i;

laplacianMatrix = ConstantArray[0.0, {nCells, nCells}];
Do[
  row = cellIndex[i, j];
  vol = radialVolumes[[i]];

  If[i > 1,
    coeff = radialFaceAreas[[i]]/(dr vol);
    laplacianMatrix[[row, cellIndex[i - 1, j]]] += coeff;
    laplacianMatrix[[row, row]] -= coeff;
  ];
  If[i < nr,
    coeff = radialFaceAreas[[i + 1]]/(dr vol);
    laplacianMatrix[[row, cellIndex[i + 1, j]]] += coeff;
    laplacianMatrix[[row, row]] -= coeff;,
    coeff = radialFaceAreas[[i + 1]]/(dr vol);
    laplacianMatrix[[row, row]] -= 2.0 coeff;
  ];

  If[j == 1,
    laplacianMatrix[[row, cellIndex[i, j + 1]]] += 1.0/dw^2;
    laplacianMatrix[[row, row]] -= 1.0/dw^2;,
    If[j == nw,
      laplacianMatrix[[row, cellIndex[i, j - 1]]] += 1.0/dw^2;
      laplacianMatrix[[row, row]] -= 1.0/dw^2;,
      laplacianMatrix[[row, cellIndex[i, j - 1]]] += 1.0/dw^2;
      laplacianMatrix[[row, cellIndex[i, j + 1]]] += 1.0/dw^2;
      laplacianMatrix[[row, row]] -= 2.0/dw^2;
    ];
  ];
  ,
  {j, 1, nw}, {i, 1, nr}
];

bdgRadii = Flatten[Table[rCenters[[i]], {j, 1, nw}, {i, 1, nr}]];
bdgW = Flatten[Table[wCenters[[j]], {j, 1, nw}, {i, 1, nr}]];
bdgX = bdgW/L;
bdgR0 = N[a + (RExit - a) (smootherStep /@ bdgX), 17];

Vconf = N[
  radialWallStrength (bdgRadii/bdgR0)^4 +
    0.5 axialTrapStrength ((bdgW - 0.5 L)/L)^2,
  17
];

rho0 = N[
  0.035 + 0.18 Exp[-(bdgRadii/(0.72 bdgR0))^2] (1.0 + 0.06 Cos[Pi bdgW/L])^2,
  17
];
psi0 = N[Sqrt[rho0], 17];
A0Values = ConstantArray[0.0, nCells];
spatialGaugeNorm = 0.0;

enthalpy = N[(5 eosK/4) rho0^4, 17];
enthalpyPrime = N[5 eosK rho0^3, 17];
enthalpyLinearization = N[rho0 enthalpyPrime, 17];

angularMatrix = DiagonalMatrix[N[ellFactor/bdgRadii^2, 17]];
covariantLaplacianMatrix = laplacianMatrix - angularMatrix;
kineticMatrix = N[-(hbar^2/(2 particleMass)) covariantLaplacianMatrix, 17];
singleParticleMatrix = kineticMatrix + DiagonalMatrix[Vconf + gaugeCharge A0Values - chemicalPotential + enthalpy];
bdgDiagonalMatrix = singleParticleMatrix + DiagonalMatrix[enthalpyLinearization];
bdgPairMatrix = DiagonalMatrix[enthalpyLinearization];
bdgMatrix = ArrayFlatten[{
    {bdgDiagonalMatrix, bdgPairMatrix},
    {-bdgPairMatrix, -bdgDiagonalMatrix}
  }];

stationaryResidualNorm = N[
  Norm[singleParticleMatrix . psi0]/Max[Norm[psi0], 10^-30],
  17
];

Print["MT15 M1b BdG + wall derivation export"];
Print["building dense BdG eigensystem with dimension: ", 2 nCells];

{bdgEigenvalues, bdgEigenvectors} = Eigensystem[N[bdgMatrix, 30]];
bdgPairs = SortBy[Transpose[{bdgEigenvalues, bdgEigenvectors}], Re[First[#]] &];
positiveBdgPairs = Select[
  bdgPairs,
  Re[#[[1]]] > 10^-8 && Abs[Im[#[[1]]]] < 10^-7 &
];

If[Length[positiveBdgPairs] < 2,
  Print["ERROR: fewer than two positive real BdG eigenmodes found: ", Length[positiveBdgPairs]];
  Exit[3]
];

modesToExport = Min[3, Length[positiveBdgPairs]];
cellVolumeWeights = N[Flatten[Table[radialVolumes[[i]] wQuadWeights[[j]], {j, 1, nw}, {i, 1, nr}]], 17];
driveCoeff = N[4.0 radialWallStrength bdgRadii^4/bdgR0^5, 17];

modeRecords = {};
bdgPacketModes = {};
couplingList = {};
varpiList = {};
bdgResidualList = {};

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
  residualRel = N[
    residualAbs/Max[Norm[bdgMatrix . modeVec], Abs[eigenvalue] Norm[modeVec], 1.0],
    17
  ];

  densityByW = Table[
    Total[
      Table[
        radialVolumes[[i]] driveCoeff[[cellIndex[i, j]]] psi0[[cellIndex[i, j]]] uvReal[[cellIndex[i, j]]],
        {i, 1, nr}
      ]
    ],
    {j, 1, nw}
  ];
  densityOnPacketGrid = Interpolation[Transpose[{wCenters, densityByW}], InterpolationOrder -> 1] /@ points;
  rawPhi = N[densityOnPacketGrid/weightValues, 17];
  phiNorm = normalizeProfile[points, weightValues, rawPhi, "bdg axial profile"];
  adapterOverlap = weightedInner[points, weightValues, wallChiEta, phiNorm];
  signedPhysicalCoupling = -trapz[points, wallChiEta densityOnPacketGrid];
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
    "profile_norm_weighted" -> toJsonNumber[weightedNorm[points, weightValues, phiNorm]]
  |>];
  AppendTo[bdgPacketModes, <|
    "name" -> "m1b_derived_bdg_mode_" <> ToString[modeIdx],
    "lambda_B" -> toJsonNumber[lambdaB],
    "varpi" -> toJsonNumber[varpi],
    "profile_values" -> toJsonVector[phiNorm],
    "derivation_status" -> "derived_eigensolve",
    "lambda_B_status" -> "normalization_convention_factor_c_alpha_equals_lambda_B_times_I_eta_phi"
  |>];
  ,
  {modeIdx, 1, modesToExport}
];

B0 = N[Total[Table[couplingList[[i]]^2/varpiList[[i]]^2, {i, 1, Length[couplingList]}]], 17];
B2 = N[Total[Table[couplingList[[i]]^2/varpiList[[i]]^4, {i, 1, Length[couplingList]}]], 17];
B4 = N[Total[Table[couplingList[[i]]^2/varpiList[[i]]^6, {i, 1, Length[couplingList]}]], 17];
stieltjesGap = N[B0 B4 - B2^2, 17];
maxBdgResidual = Max[bdgResidualList];
minExportedVarpi = Min[varpiList];

(* Mixed-Maxwell ports: posited placeholders, converted to packet lambda factors
   so the V2 adapter recovers the declared effective R, g_U, and g_W. *)
uProfile = wallChiEta;
wRaw = N[Sin[Pi points/(2 L)] + 0.20 Sin[3 Pi points/(2 L)], 17];
wProfile = normalizeProfile[points, weightValues, wRaw, "mixed W profile"];
If[weightedInner[points, weightValues, wallChiEta, wProfile] < 0, wProfile = -wProfile];

IetaU = weightedInner[points, weightValues, wallChiEta, uProfile];
IetaW = weightedInner[points, weightValues, wallChiEta, wProfile];
Iuw = weightedInner[points, weightValues, uProfile, wProfile];

gUPosited = 0.18;
gWPosited = 0.13;
RPosited = 0.08;
omegaUPosited = 3.25;
omegaWPosited = 4.35;

lambdaU = gUPosited/IetaU;
lambdaW = gWPosited/IetaW;
lambdaR = RPosited/Iuw;
deltaEff = omegaUPosited^2 omegaWPosited^2 - RPosited^2;

machineryNotPhysicsNote = "The exported varpi_alpha, c_alpha, and B0/B2/B4 are eigenmodes/couplings of the operator linearized around a NON-EQUILIBRIUM target-blind background (stationary residual about 243). They are correctly-derived MACHINERY, NOT physical frequencies/couplings; physical values require the self-consistent background at M1c.";
portDominanceNote = "R_norm, R_pole, R_P2, and R_P4 are partially/largely determined by the POSITED mixed ports (g_U=0.18, g_W=0.13, R=0.08, Omega_U=3.25, Omega_W=4.35), so they are NOT clean tests of the derived BdG/wall sectors.";
rNormProximityNote = "The proximity |R_norm| ~= 0.1106 to the GR target constant 0.110592 is a calibration artifact of the near-zero posited N0 (N0 ~= 4.8e-5, making the mhat0^2*S_port*N0/D0 transfer term about 2.8e-5), NOT validation and NOT tuning.";
cAlphaCaveatNote = "c_alpha is a defensible realization of the delta V_conf wall<->BdG coupling (symplectically-normalized mode density-response projected through the drive onto the wall shape): derived engineering realization; canonical Schur-complement equivalence pending M1c.";
stieltjesIdentityNote = "B0*B4-B2^2 >= 0 by Cauchy-Schwarz for any positive (c,varpi); this is a structural sanity identity that catches a coding bug only, NOT a physics validation gate.";

packet = <|
  "schema" -> "stage_v2_22b_solver_handoff/v1",
  "branch_id" -> "mt15_m1b_pre_freeze_derived_bdg_wall_posited_mixed_ports",
  "freeze" -> <|
    "pre_target_freeze" -> True,
    "target_blind" -> True,
    "parent_action" -> "GNLS BdG eigensolve + weak-form wall closure; mixed Maxwell ports carried as posited placeholders",
    "parent_action_status" -> "effective_wall_closure",
    "gauge_convention" -> "localized_H_equals_Z_with_zero_spatial_gauge_background_for_m1b_smoke",
    "boundary_protocol" -> "open_impedance_AC_reflecting_DC_leaking",
    "source_map_convention" -> "canonical real STF grouped P2 basis",
    "support_profile_family" -> "derived weak-form wall eigen-shape on target-blind open-throat geometry"
  |>,
  "geometry" -> <|
    "L" -> toJsonNumber[L],
    "R_mouth" -> toJsonNumber[a],
    "R_exit" -> toJsonNumber[RExit],
    "boundary_class" -> "open_impedance",
    "Y_L_limit" -> toJsonNumber[YL0],
    "exit_model" -> "impedance_mismatch_open_exit"
  |>,
  "constants" -> <|
    "G" -> toJsonNumber[GConst],
    "c_s" -> toJsonNumber[csConst],
    "c" -> toJsonNumber[cConst],
    "a" -> toJsonNumber[a],
    "mhat0" -> toJsonNumber[mhat0Const],
    "S_port" -> toJsonNumber[SPortConst],
    "theta_tail" -> toJsonNumber[thetaTailConst]
  |>,
  "grid" -> <|
    "coordinate" -> "w",
    "points" -> toJsonVector[points]
  |>,
  "wall" -> <|
    "K" -> toJsonNumber[wallKScalar],
    "M" -> toJsonNumber[wallMScalar],
    "derivation_status" -> "derived_weak_form_galerkin"
  |>,
  "profiles" -> <|
    "radius_R0" -> toJsonVector[radiusValues],
    "weight" -> toJsonVector[weightValues],
    "wall_chi_eta" -> toJsonVector[wallChiEta]
  |>,
  "bdg_modes" -> bdgPacketModes,
  "mixed_ports" -> {
    <|
      "name" -> "m1b_posited_mixed_port_open_frontier",
      "lambda_U" -> toJsonNumber[lambdaU],
      "Omega_U" -> toJsonNumber[omegaUPosited],
      "u_values" -> toJsonVector[uProfile],
      "lambda_W" -> toJsonNumber[lambdaW],
      "Omega_W" -> toJsonNumber[omegaWPosited],
      "w_values" -> toJsonVector[wProfile],
      "lambda_R" -> toJsonNumber[lambdaR],
      "mixed_ports_status" -> "posited_not_derived",
      "effective_R_posited" -> toJsonNumber[RPosited],
      "effective_g_U_posited" -> toJsonNumber[gUPosited],
      "effective_g_W_posited" -> toJsonNumber[gWPosited]
    |>
  },
  "solver_metadata" -> <|
    "exporter" -> "mt15_02_bdg_wall_derivation.wls",
    "mesh_points" -> meshPoints,
    "coefficient_family" -> "m1b_target_blind_derived_bdg_wall_posited_mixed_ports_v1",
    "source_commit" -> sourceCommitString[],
    "export_timestamp" -> DateString[Now, {"ISODate", "T", "Time", " ", "TimeZoneNameShort"}],
    "nonlinear_residual_norm" -> toJsonNumber[stationaryResidualNorm],
    "linear_residual_norm" -> toJsonNumber[maxBdgResidual],
    "bdg_status" -> "derived_from_finite_dimensional_BdG_eigensolve",
    "wall_status" -> "derived_from_weak_form_integrals",
    "mixed_ports_status" -> "posited_not_derived",
    "notes" -> StringRiffle[
      {
        machineryNotPhysicsNote,
        "BdG + wall DERIVED; mixed ports POSITED; mixed_ports_status remains posited_not_derived.",
        portDominanceNote,
        rNormProximityNote,
        cAlphaCaveatNote
      },
      " "
    ]
  |>,
  "normalization_tolerances" -> <|
    "profile_norm_tol" -> 0.005,
    "Delta_tol" -> 1.0*^-12
  |>
|>;

diagnostics = <|
  "schema" -> "mt15_m1b_derivation_diagnostics/v1",
  "packet_path" -> packetPath,
  "integrity_notes" -> <|
    "machinery_not_physics" -> machineryNotPhysicsNote,
    "residual_port_dominance" -> portDominanceNote,
    "R_norm_proximity_artifact" -> rNormProximityNote,
    "c_alpha_status" -> cAlphaCaveatNote,
    "stieltjes_identity" -> stieltjesIdentityNote
  |>,
  "background" -> <|
    "status" -> "target_blind_engineering_smoke_not_full_2d_newton_solution",
    "machinery_not_physics_note" -> machineryNotPhysicsNote,
    "geometry" -> "R0(w)=a+(R_exit-a)(3x^2-2x^3)",
    "rho0_formula" -> "rho0=0.035+0.18 Exp[-(r/(0.72 R0(w)))^2] (1+0.06 Cos[Pi w/L])^2",
    "psi0" -> "Sqrt[rho0]",
    "stationary_residual_norm" -> toJsonNumber[stationaryResidualNorm]
  |>,
  "operator_term_map" -> {
    <|"term" -> "kinetic", "implementation" -> "-(hbar^2/(2m))*(laplacian-l(l+1)/r^2)", "source" -> "compact L1406-L1428; p2_tangent.py L240-L245"|>,
    <|"term" -> "confinement", "implementation" -> "V_radial*(r/R0)^4+0.5*V_axial*((w-L/2)/L)^2", "source" -> "coupled_branch.py L260-L270 and confinement_potential_torch"|>,
    <|"term" -> "EOS enthalpy", "implementation" -> "h=(5K/4) rho^4", "source" -> "compact L578-L582; physics.py L46-L47"|>,
    <|"term" -> "EOS linearization", "implementation" -> "delta(h psi)=h delta_psi+(rho0 dh/drho)(delta_psi+delta_psi*)", "source" -> "linearization of compact EOS L578-L582"|>,
    <|"term" -> "chemical potential", "implementation" -> "-mu on the single-particle diagonal", "source" -> "coupled_branch.py L263-L270"|>,
    <|"term" -> "gauge coupling", "implementation" -> "q A0 retained; spatial covariant derivative terms vanish because Ar=Aw=0 in this smoke background", "source" -> "compact L569-L573; coupled_branch.py L253-L258"|>,
    <|"term" -> "wall drive", "implementation" -> "delta V_conf=-4 V_radial r^4/R0^5 eta", "source" -> "compact L1080-L1085, L1424-L1428; p2_tangent.py L165-L183"|>
  },
  "bdg_discretization" -> <|
    "radial_cells" -> nr,
    "axial_points" -> nw,
    "matrix_dimension" -> 2 nCells,
    "positive_real_modes_found" -> Length[positiveBdgPairs],
    "modes_exported" -> modesToExport,
    "spatial_gauge_norm" -> toJsonNumber[spatialGaugeNorm],
    "min_exported_varpi" -> toJsonNumber[minExportedVarpi],
    "max_relative_eigen_residual" -> toJsonNumber[maxBdgResidual]
  |>,
  "bdg_modes" -> modeRecords,
  "bdg_moments" -> <|
    "B0" -> toJsonNumber[B0],
    "B2" -> toJsonNumber[B2],
    "B4" -> toJsonNumber[B4],
    "B0_B4_minus_B2_squared" -> toJsonNumber[stieltjesGap],
    "stieltjes_identity_cauchy_schwarz_guaranteed" -> <|
      "value_key" -> "B0_B4_minus_B2_squared",
      "note" -> stieltjesIdentityNote
    |>
  |>,
  "wall_weak_form" -> <|
    "basis" -> {"w/L", "(w/L)^2", "(w/L)^3"},
    "mu_eta" -> "1+0.08 Cos[Pi w/L]^2",
    "T_w" -> "1.1+0.10 Sin[Pi w/L]",
    "T_Omega" -> "0.35+0.03 Cos[2 Pi w/L]",
    "K_eta" -> "0.45+0.08*256*x^4*(1-x)^4",
    "V_l" -> "K_eta+l(l+1)T_Omega with l=2",
    "Y_L_0" -> toJsonNumber[YL0],
    "M_matrix" -> toJsonMatrix[wallMMatrix],
    "K_matrix" -> toJsonMatrix[wallKMatrix],
    "M_eigenvalues" -> toJsonVector[Eigenvalues[N[wallMMatrix, 30]]],
    "M_positive" -> TrueQ[wallMMinEigenvalue > 0],
    "active_wall_coefficients_weight_normalized" -> toJsonVector[wallCoeff],
    "active_wall_K_scalar" -> toJsonNumber[wallKScalar],
    "active_wall_M_scalar" -> toJsonNumber[wallMScalar],
    "active_wall_weight_norm" -> toJsonNumber[weightedNorm[points, weightValues, wallChiEta]]
  |>,
  "mixed_ports" -> <|
    "status" -> "posited_not_derived",
    "frontier_note" -> "No canonical 1-D mixed Maxwell eigenproblem is used here; g_U/g_W are Path-A-blocked placeholders.",
    "residual_port_dominance_note" -> portDominanceNote,
    "Omega_U" -> toJsonNumber[omegaUPosited],
    "Omega_W" -> toJsonNumber[omegaWPosited],
    "R" -> toJsonNumber[RPosited],
    "g_U" -> toJsonNumber[gUPosited],
    "g_W" -> toJsonNumber[gWPosited],
    "lambda_U_packet_factor" -> toJsonNumber[lambdaU],
    "lambda_W_packet_factor" -> toJsonNumber[lambdaW],
    "lambda_R_packet_factor" -> toJsonNumber[lambdaR],
    "I_eta_u" -> toJsonNumber[IetaU],
    "I_eta_w" -> toJsonNumber[IetaW],
    "I_u_w" -> toJsonNumber[Iuw],
    "Delta_eff" -> toJsonNumber[deltaEff],
    "Delta_eff_gt_1e_minus_12" -> TrueQ[deltaEff > 10^-12]
  |>,
  "packet_checks" -> <|
    "wall_profile_weight_norm" -> toJsonNumber[weightedNorm[points, weightValues, wallChiEta]],
    "u_profile_weight_norm" -> toJsonNumber[weightedNorm[points, weightValues, uProfile]],
    "w_profile_weight_norm" -> toJsonNumber[weightedNorm[points, weightValues, wProfile]],
    "R_exit_positive" -> TrueQ[RExit > 0],
    "boundary_class" -> "open_impedance"
  |>
|>;

Export[packetPath, packet, "JSON"];
Export[diagnosticsPath, diagnostics, "JSON"];

Print["packet_path: ", packetPath];
Print["diagnostics_path: ", diagnosticsPath];
Print["background_status: target-blind engineering smoke, not full 2-D Newton"];
Print["machinery_not_physics: ", machineryNotPhysicsNote];
Print["BdG positive modes found: ", Length[positiveBdgPairs], " exported: ", modesToExport];
Print["BdG varpi exported: ", ToString[NumberForm[varpiList, {16, 10}]]];
Print["BdG max relative eigen-residual: ", ToString[ScientificForm[maxBdgResidual, 6]]];
Print["BdG moments B0/B2/B4: ", ToString[NumberForm[{B0, B2, B4}, {16, 10}]]];
Print["Stieltjes Cauchy-Schwarz identity B0*B4-B2^2: ", ToString[ScientificForm[stieltjesGap, 6]]];
Print["Wall scalar K/M: ", ToString[NumberForm[{wallKScalar, wallMScalar}, {16, 10}]]];
Print["Wall M positive: ", TrueQ[wallMMinEigenvalue > 0]];
Print["Mixed ports status: posited_not_derived"];
Print["Mixed effective Omega_U/Omega_W/R/g_U/g_W: ", ToString[NumberForm[{omegaUPosited, omegaWPosited, RPosited, gUPosited, gWPosited}, {16, 10}]]];
Print["Mixed Delta_eff: ", ToString[NumberForm[deltaEff, {16, 10}]]];
