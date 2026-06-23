(* pathA_23 Stage 3b over-count and curvature-localization audit.

   CONDITIONAL path: the rotational/MacCullagh brane law is the user-postulated
   constitutive package.  This script checks whether that brane stress enters
   the bulk-field Euler variations, and whether flat in-plane light transmits
   tangential traction to an ideal scalar bulk.  Negative controls are included
   for the load-bearing decisions. *)

ClearAll[
  dimString, checkDim, checkBoolean, homogeneous, zeroVecQ, nonzeroVecQ,
  containsQ, L, T, M, Q, energy, braneLag, bulkLag, Bdim, uDim, rho3Dim,
  rho4Dim, velocityDim, stressBulkDim, braneForceDim, bulkForceDim, kDim,
  curvatureDim, kx, ky, kz, k, kvec, k2, pTnum, ux, uy, uz, uvec, Uamp,
  phi, axx, axy, axz, ayx, ayy, ayz, azx, azy, azz, A, gradVars, muR,
  Bl, rho, theta, vx, vy, vz, vw, vpar, udotw, sx, sy, sz, slopes, alphaU,
  Jx, Jy, Jz, J0, Phi, beta, eta, piN, Lpsi, Lbulk, Lbrane, Scpl, Urot,
  rvec, sigmaR, sigmaExpected, sigmaFourier, divSigmaR, divSigmaExpected,
  Ldeclared, LnoBrane, bulkVars, bulkVariationsDeclared,
  bulkVariationsNoBrane, badVelocityCoupling, badRhoBraneCoupling,
  badBulkVelocityVariation, badBulkRhoVariation, currentCoupling,
  currentBulkVelocityVariation, badCurrentVelocityCoupling, badCurrentVar,
  mG, P, sigQwx, sigQwy, sigQwz, sigQ, Tflat, pressureTractionFlat,
  convectiveTractionFlat, quantumTractionFlat, totalTractionFlat,
  hbar, rho0, rho0p, deltaRho, deltaRhoW, rhoBr, sigQwaLinear,
  deltaRhoShear, sigQShearGeneral, transverseRules, compressiveRules,
  vWCompat, tFlatLight, tFlatBadNormalFlow, n4, eTangents, v4, Tideal,
  pressureTiltTraction, tangentFlowRules, tiltedIdealTraction,
  tiltedIdealTractionLinear, fixedChartLinearTraction, Kxx, Kxy, Kxz,
  Kyx, Kyy, Kyz, Kzx, Kzy, Kzz, x1, x2, x3, Kmat, xvec, localSlope,
  Axx, Axy, Axz, Ayx, Ayy, Ayz, Azx, Azy, Azz, Amix, curvatureMixing,
  curvatureOffRules, curvatureOnRules, Lmix, FrefDim, epsAmpDim,
  epsPowerDim, dimChecks, variationalChecks, tractionChecks, flatChecks,
  curvatureChecks, negativeChecks, checks, allPass, sigmaRNotBulk,
  sigmaRSourcesBulk, flatLightNoLeak, curvatureLocalized, primaryToken,
  report, outDir, jsonPath
];

dimString[d_] := Module[
  {labels = {"L", "T", "M", "Q"}, pairs, nonzero},
  pairs = Transpose[{labels, d}];
  nonzero = Select[pairs, #[[2]] =!= 0 &];
  If[
    nonzero === {},
    "1",
    StringRiffle[
      Map[If[#[[2]] === 1, #[[1]], #[[1]] <> "^" <> ToString[InputForm[#[[2]]]]] &, nonzero],
      " "
    ]
  ]
];

checkDim[name_, actual_, expected_, note_: ""] := <|
  "name" -> name,
  "pass" -> TrueQ[actual === expected],
  "expected" -> dimString[expected],
  "actual" -> dimString[actual],
  "note" -> note
|>;

checkBoolean[name_, actual_, expected_: True, note_: ""] := <|
  "name" -> name,
  "pass" -> TrueQ[actual === expected],
  "expected" -> expected,
  "actual" -> actual,
  "note" -> note
|>;

homogeneous[name_, terms_Association, expected_, note_: ""] := <|
  "name" -> name,
  "pass" -> AllTrue[Values[terms], TrueQ[# === expected] &],
  "expected" -> dimString[expected],
  "terms" -> Association @ KeyValueMap[#1 -> dimString[#2] &, terms],
  "note" -> note
|>;

zeroVecQ[expr_] := AllTrue[Flatten[FullSimplify[expr]], TrueQ[# === 0] &];
nonzeroVecQ[expr_] := Not[zeroVecQ[expr]];
containsQ[expr_, sym_] := Not[FreeQ[expr, sym]];

L = {1, 0, 0, 0};
T = {0, 1, 0, 0};
M = {0, 0, 1, 0};
Q = {0, 0, 0, 1};

energy = M + 2 L - 2 T;
braneLag = energy - 3 L;
bulkLag = energy - 4 L;
Bdim = -L;
uDim = L;
rho3Dim = M - 3 L;
rho4Dim = M - 4 L;
velocityDim = L - T;
stressBulkDim = bulkLag;
braneForceDim = braneLag - L;
bulkForceDim = bulkLag - L;
kDim = -L;
curvatureDim = -L;

kvec = {kx, ky, kz};
k2 = kx^2 + ky^2 + kz^2;
pTnum = k2 IdentityMatrix[3] - Outer[Times, kvec, kvec];
uvec = {ux, uy, uz};
vpar = {vx, vy, vz};
slopes = {sx, sy, sz};

A = {
  {axx, axy, axz},
  {ayx, ayy, ayz},
  {azx, azy, azz}
};
gradVars = Flatten[A];

(* Rotational/MacCullagh stress derived from the brane energy. *)
rvec = {ayz - azy, azx - axz, axy - ayx};
Urot = (muR/2) (rvec . rvec);
sigmaR = FullSimplify[Table[D[Urot, A[[a, b]]], {a, 1, 3}, {b, 1, 3}]];
sigmaExpected = muR (A - Transpose[A]);
sigmaFourier = FullSimplify[
  sigmaR /. Thread[gradVars -> Flatten[Outer[Times, I kvec, uvec]]]
];
divSigmaR = FullSimplify[
  Table[Sum[I kvec[[a]] sigmaFourier[[a, b]], {a, 1, 3}], {b, 1, 3}]
];
divSigmaExpected = FullSimplify[muR (kvec (kvec . uvec) - k2 uvec)];

(* Full declared density in finite-thickness form.  piN represents the
   Stage-0 normal coupling functional Pi_n[psi,Sigma].  It may depend on bulk
   fields and geometry, but the brane rotational energy is an independent
   functional of u_a. *)
bulkVars = {rho, theta, vx, vy, vz, vw};
Lbulk = Lpsi[rho, theta, vx, vy, vz, vw];
Lbrane = -Urot;
Scpl = uw piN[rho, theta, vx, vy, vz, vw, sx, sy, sz] +
  alphaU (Jx ux + Jy uy + Jz uz) - J0 Phi;
Ldeclared = Lbulk + Bl (Lbrane + Scpl);
LnoBrane = Lbulk + Bl Scpl;
bulkVariationsDeclared = FullSimplify[D[Ldeclared, #] & /@ bulkVars];
bulkVariationsNoBrane = FullSimplify[D[LnoBrane, #] & /@ bulkVars];

(* Able-to-fail controls: these are not in the declared action.  If either
   were added, the bulk variation would contain the rotational stress. *)
badVelocityCoupling = Bl beta (vpar . divSigmaR);
badRhoBraneCoupling = Bl eta rho Urot;
badBulkVelocityVariation = FullSimplify[
  D[Ldeclared + badVelocityCoupling, #] & /@ {vx, vy, vz}
];
badBulkRhoVariation = FullSimplify[D[Ldeclared + badRhoBraneCoupling, rho]];

currentCoupling = Bl (alphaU (Jx ux + Jy uy + Jz uz) - J0 Phi);
currentBulkVelocityVariation = FullSimplify[D[currentCoupling, #] & /@ {vx, vy, vz, vw}];
badCurrentVelocityCoupling = Bl beta (Jx vx + Jy vy + Jz vz);
badCurrentVar = FullSimplify[D[badCurrentVelocityCoupling, #] & /@ {vx, vy, vz}];

(* Flat-surface tangential traction t_a = n_i T_ia with n = w-hat. *)
sigQ = {sigQwx, sigQwy, sigQwz};
Tflat = mG rho Outer[Times, {vx, vy, vz, vw}, {vx, vy, vz, vw}] +
  P IdentityMatrix[4] +
  {
    {0, 0, 0, sigQwx},
    {0, 0, 0, sigQwy},
    {0, 0, 0, sigQwz},
    {sigQwx, sigQwy, sigQwz, 0}
  };
pressureTractionFlat = Table[(UnitVector[4, 4] . (P IdentityMatrix[4]) . UnitVector[4, a]), {a, 1, 3}];
convectiveTractionFlat = FullSimplify[Table[mG rho vw vpar[[a]], {a, 1, 3}]];
quantumTractionFlat = sigQ;
totalTractionFlat = FullSimplify[Table[Tflat[[4, a]], {a, 1, 3}]];

(* Quantum stress for rho = rho0(w) + deltaRho(x,w).  The linear mixed
   normal/tangential component is proportional to the density perturbation.
   A density-preserving transverse brane shear has deltaRho = -i rho_br k.u=0. *)
sigQwaLinear = FullSimplify[
  (hbar^2/(4 mG)) I kvec (rho0p deltaRho/rho0 - deltaRhoW)
];
deltaRhoShear = -I rhoBr (kvec . uvec);
sigQShearGeneral = FullSimplify[sigQwaLinear /. {deltaRho -> deltaRhoShear, deltaRhoW -> 0}];
transverseRules = {kx -> 0, ky -> 0, kz -> k, ux -> Uamp, uy -> 0, uz -> 0};
compressiveRules = {kx -> 0, ky -> 0, kz -> k, ux -> 0, uy -> 0, uz -> Uamp};

(* Kinematic normal compatibility.  With no normal brane motion,
   v_n = v_w - v_parallel.s = dot u_w gives v_w = dot u_w + v_parallel.s. *)
vWCompat = udotw + vpar . slopes;
tFlatLight = FullSimplify[
  (mG rho vWCompat vpar + sigQShearGeneral) /.
    {udotw -> 0, sx -> 0, sy -> 0, sz -> 0} /. transverseRules
];
tFlatBadNormalFlow = FullSimplify[(mG rho vw vpar + sigQShearGeneral) /. transverseRules];

(* Constant tilt/free-slip control.  Pressure contributes no tangential force,
   and a velocity tangent to the tilted surface has n.v = 0, so the ideal
   convective traction vanishes too. *)
n4 = {-sx, -sy, -sz, 1};
eTangents = Table[UnitVector[4, a] + slopes[[a]] UnitVector[4, 4], {a, 1, 3}];
v4 = {vx, vy, vz, vw};
Tideal = mG rho Outer[Times, v4, v4] + P IdentityMatrix[4];
pressureTiltTraction = FullSimplify[Table[n4 . (P IdentityMatrix[4]) . eTangents[[a]], {a, 1, 3}]];
tangentFlowRules = {vw -> vpar . slopes};
tiltedIdealTraction = FullSimplify[Table[n4 . Tideal . eTangents[[a]], {a, 1, 3}] /. tangentFlowRules];
tiltedIdealTractionLinear = FullSimplify[
  Normal[Series[tiltedIdealTraction /. {sx -> eps sx, sy -> eps sy, sz -> eps sz}, {eps, 0, 1}]] /. eps -> 1
];
fixedChartLinearTraction = FullSimplify[
  (mG rho vw vpar + (P IdentityMatrix[3] - (mG rho Outer[Times, vpar, vpar] + P IdentityMatrix[3])) . slopes) /.
    tangentFlowRules
];

(* Curvature localization in a local tangent frame.  A flat tilted patch has no
   invariant leak; over a finite patch the residual slope is Delta s ~= K x.
   The unresolved defect/background anisotropy channel scales with that local
   slope and therefore vanishes with K. *)
Kmat = {
  {Kxx, Kxy, Kxz},
  {Kyx, Kyy, Kyz},
  {Kzx, Kzy, Kzz}
};
xvec = {x1, x2, x3};
localSlope = Kmat . xvec;
Amix = {
  {Axx, Axy, Axz},
  {Ayx, Ayy, Ayz},
  {Azx, Azy, Azz}
};
curvatureMixing = FullSimplify[Amix . localSlope];
curvatureOffRules = Thread[Flatten[Kmat] -> ConstantArray[0, 9]];
curvatureOnRules = {
  Kxy -> 0, Kxz -> 0, Kyx -> 0, Kyy -> 0, Kyz -> 0,
  Kzx -> 0, Kzy -> 0, Kzz -> 0,
  x1 -> Lmix, x2 -> 0, x3 -> 0
};

FrefDim = braneForceDim;
epsAmpDim = stressBulkDim - FrefDim; (* stress-force ratio, then times K Lmix *)
epsPowerDim = 2 (curvatureDim + L) + stressBulkDim - FrefDim;

dimChecks = {
  homogeneous[
    "declared finite-thickness density pieces",
    <|
      "bulk scalar Lpsi" -> bulkLag,
      "B_l U_R" -> Bdim + braneLag,
      "B_l u_w Pi_n" -> Bdim + uDim + stressBulkDim,
      "B_l J^a alpha_u u_a / J0 Phi" -> Bdim + braneLag
    |>,
    bulkLag
  ],
  checkDim["sigma_R stress", braneLag, braneLag],
  checkDim["D_b sigma_R brane force density", braneForceDim, braneLag - L],
  checkDim["bulk tangential traction T_wa", stressBulkDim, braneForceDim],
  checkDim["curvature K_ab", curvatureDim, -L],
  checkDim["geometric amplitude factor K L_mix", curvatureDim + L, {0, 0, 0, 0}],
  checkDim["epsilon amplitude stress/reference factor before K L_mix", epsAmpDim, {0, 0, 0, 0}],
  checkDim["epsilon power-style curvature-squared factor", epsPowerDim, {0, 0, 0, 0}]
};

variationalChecks = {
  checkBoolean["sigma_R_ab derives from U_R", zeroVecQ[Flatten[sigmaR - sigmaExpected]]],
  checkBoolean["D_b sigma_R_ba Fourier form derives correctly", zeroVecQ[divSigmaR - divSigmaExpected]],
  checkBoolean["declared bulk variations equal variations of S_psi + S_cpl only", zeroVecQ[bulkVariationsDeclared - bulkVariationsNoBrane]],
  checkBoolean["declared bulk variations contain no mu_R", Not[containsQ[bulkVariationsDeclared, muR]]],
  checkBoolean["declared bulk variations contain no D_b sigma_R_ba", Not[containsQ[bulkVariationsDeclared, divSigmaR]]],
  checkBoolean["declared source-current coupling has no direct bulk-velocity variation", zeroVecQ[currentBulkVelocityVariation]]
};

tractionChecks = {
  checkBoolean["flat pressure traction n_i P delta_ia is tangential-zero", zeroVecQ[pressureTractionFlat]],
  checkBoolean["flat ideal tangential traction is convective plus quantum only", zeroVecQ[totalTractionFlat - convectiveTractionFlat - quantumTractionFlat]],
  checkBoolean["quantum mixed traction is nonzero for generic density perturbation", nonzeroVecQ[sigQwaLinear]],
  checkBoolean["quantum mixed traction vanishes for density-preserving transverse shear", zeroVecQ[FullSimplify[sigQShearGeneral /. transverseRules]]],
  checkBoolean["quantum mixed traction would survive for compressive density perturbation", nonzeroVecQ[FullSimplify[sigQShearGeneral /. compressiveRules]]]
};

flatChecks = {
  checkBoolean["kinematic compatibility gives v_w(light)=0 on flat brane", TrueQ[FullSimplify[vWCompat /. {udotw -> 0, sx -> 0, sy -> 0, sz -> 0}] === 0]],
  checkBoolean["flat transverse in-plane light has T_wa=0 including quantum stress", zeroVecQ[tFlatLight]],
  checkBoolean["flat light would leak if an independent v_w were present", nonzeroVecQ[tFlatBadNormalFlow]]
};

curvatureChecks = {
  checkBoolean["constant tilted ideal free-slip surface has zero pressure tangential traction", zeroVecQ[pressureTiltTraction]],
  checkBoolean["constant tilted ideal free-slip surface has zero convective tangential traction for v.n=0", zeroVecQ[tiltedIdealTraction]],
  checkBoolean["fixed-chart linear slope terms cancel for tangent ideal flow", zeroVecQ[fixedChartLinearTraction]],
  checkBoolean["curvature-local mixing vanishes when K_ab=0", zeroVecQ[FullSimplify[curvatureMixing /. curvatureOffRules]]],
  checkBoolean["curvature-local mixing is generically nonzero when K_ab != 0", nonzeroVecQ[FullSimplify[curvatureMixing /. curvatureOnRules]]]
};

negativeChecks = {
  checkBoolean["negative control: v_a D_b sigma_R_ba coupling would source bulk velocity equation", containsQ[badBulkVelocityVariation, muR] && nonzeroVecQ[badBulkVelocityVariation]],
  checkBoolean["negative control: rho U_R coupling would source bulk density equation", containsQ[badBulkRhoVariation, muR] && Not[TrueQ[FullSimplify[badBulkRhoVariation] === 0]]],
  checkBoolean["negative control: forbidden J^a-to-bulk-velocity coupling would be detected", zeroVecQ[badCurrentVar] === False]
};

checks = Join[dimChecks, variationalChecks, tractionChecks, flatChecks, curvatureChecks, negativeChecks];
allPass = AllTrue[checks, TrueQ[#["pass"]] &];

sigmaRNotBulk = zeroVecQ[bulkVariationsDeclared - bulkVariationsNoBrane] &&
  Not[containsQ[bulkVariationsDeclared, muR]];
sigmaRSourcesBulk = Not[sigmaRNotBulk];
flatLightNoLeak = zeroVecQ[tFlatLight];
curvatureLocalized = zeroVecQ[FullSimplify[curvatureMixing /. curvatureOffRules]] &&
  nonzeroVecQ[FullSimplify[curvatureMixing /. curvatureOnRules]];

primaryToken = Which[
  allPass && sigmaRNotBulk && flatLightNoLeak && curvatureLocalized,
    "OVER_COUNT_CONFIRMED_CURVATURE_LOCALIZED",
  allPass && (sigmaRSourcesBulk || Not[flatLightNoLeak]),
    "INTRINSIC_LEAK_CONFIRMED_FATAL",
  allPass,
    "MIXED_PARTIAL",
  True,
    "SCRIPT_CHECK_FAILURE"
];

report = <|
  "stage" -> "pathA_23 Stage 3b over-count and curvature localization",
  "conditional_flag" -> "CONDITIONAL",
  "primary_token" -> primaryToken,
  "sigma_R_bulk_source_token" -> If[sigmaRNotBulk, "SIGMA_R_NOT_A_BULK_SOURCE", "SIGMA_R_SOURCES_BULK"],
  "flat_light_token" -> If[flatLightNoLeak, "LIGHT_FREE_SLIPS_NO_LEAK", "LIGHT_LEAKS_FLAT"],
  "checks_passed" -> Count[checks[[All, "pass"]], True],
  "checks_total" -> Length[checks],
  "all_pass" -> allPass,
  "load_bearing_expressions" -> <|
    "bulk_variations_declared" -> ToString[InputForm[bulkVariationsDeclared]],
    "bulk_variations_without_brane" -> ToString[InputForm[bulkVariationsNoBrane]],
    "div_sigma_R_fourier" -> ToString[InputForm[divSigmaR]],
    "flat_tangential_traction" -> ToString[InputForm[totalTractionFlat]],
    "flat_pressure_traction" -> ToString[InputForm[pressureTractionFlat]],
    "flat_convective_traction" -> ToString[InputForm[convectiveTractionFlat]],
    "quantum_sigma_Q_wa_linear" -> ToString[InputForm[sigQwaLinear]],
    "flat_light_T_wa" -> ToString[InputForm[tFlatLight]],
    "tilted_ideal_traction_tangent_flow" -> ToString[InputForm[tiltedIdealTraction]],
    "curvature_mixing_A_K_x" -> ToString[InputForm[curvatureMixing]],
    "epsilon_force_amplitude_scaling" -> "(||A_mix||/F_ref) |K| L_mix",
    "epsilon_power_scaling_if_energy_fraction" -> "(||A_mix||/F_ref) (|K| L_mix)^2, when epsilon is defined as leaked power/probability rather than force amplitude"
  |>,
  "negative_controls" -> <|
    "bad_velocity_sigma_coupling_variation" -> ToString[InputForm[badBulkVelocityVariation]],
    "bad_rho_brane_coupling_variation" -> ToString[InputForm[badBulkRhoVariation]],
    "bad_current_velocity_variation" -> ToString[InputForm[badCurrentVar]],
    "bad_flat_normal_flow_traction" -> ToString[InputForm[tFlatBadNormalFlow]]
  |>,
  "checks" -> checks
|>;

outDir = FileNameJoin[{"software", "stage1_solver", "_scratch"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
jsonPath = FileNameJoin[{outDir, "pathA_23_stage3b_overcount_and_curvature_mathematica.json"}];
Export[jsonPath, report, "JSON"];

If[allPass, Exit[0], Exit[1]];
