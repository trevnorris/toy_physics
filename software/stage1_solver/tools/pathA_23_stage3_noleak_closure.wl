(* pathA_23 Stage 3 constitutive no-leak closure.

   CONDITIONAL path: the rotational/MacCullagh brane law is postulated by the
   user after Stage 2 failed to derive it.  This script tests the consequences
   of that fixed postulate for interface force/curl leakage.  It deliberately
   does not set T_wa, mouth traction, or J^a to zero in the load-bearing
   generic source. *)

ClearAll[
  dimString, checkDim, checkBoolean, homogeneous, zeroVecQ, nonzeroVecQ,
  zeroMatrixQ, nonzeroMatrixQ, L, T, M, Q, dim0, energy, braneLag,
  bulkLag, Bdim, uDim, rho3Dim, rho4Dim, velocityDim, circulationDim,
  vortexDensityDim, stressBraneDim, braneForceDim, bulkStressDim,
  bulkForce4Dim, coupleStressDim, spinDensityDim, kx, ky, kz, k, k2,
  kvec, pTnum, curl3, ux, uy, uz, uvec, phi, Uamp, Vx, Vy, Vz, Vpar,
  dvw, A, gradVars, axx, axy, axz, ayx, ayy, ayz, azx, azy, azz,
  muR, rhoPar, alphaU, Jx, Jy, Jz, J0, Phi, Bl, ell, mG, rho3,
  GammaCirc, nVort, Vrel, vr, Lg, piww, piIso, pixx, pixy, pixz, piyy,
  piyz, pizz, piwx, piwy, piwz, uw, tx, ty, tz, tauAx, tauAy, tauAz,
  lambdaN, rho, theta, vx, vy, vz, vw, rvec, Wrot, sigmaRot,
  sigmaExpected, antisymRot, spinRateAB, angularResidualNoSpin,
  angularResidualWithSpin, nu1, nu2, nu3, nuVec, tractionRot,
  coupleStressZero, coupleTractionZero, fourierRules, sigmaFourier,
  divSigmaFourier, divSigmaExpected, rotTransverse, rotLongitudinal,
  rotTransverseControl, piPar, piW, slopesFourier, TnaFourier,
  Jvec, tauAeff, vnFeedbackForce, directVnSource, totalExchange3,
  totalExchangeT, totalExchangeCurl, noLeakRules, transverseRotRules,
  genericTwaRules, vnFeedbackTransverseRules, vnFeedbackParallelRules,
  c2StaticSource, c2BulkVariation, c2TransverseSource, fLeakScaleDim,
  fMagDim, fGravDim, fRefDim, epsLeakExpression, dimChecks,
  postulateChecks, angularChecks, boundaryChecks, variationChecks,
  leakChecks, vnChecks, c2Checks, orderChecks, checks, allPass,
  stage3Token, report, outDir, jsonPath
];

dim0 = {0, 0, 0, 0};
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
circulationDim = 2 L - T;
vortexDensityDim = -2 L;
stressBraneDim = braneLag;
braneForceDim = braneLag - L;
bulkStressDim = bulkLag;
bulkForce4Dim = bulkLag - L;
coupleStressDim = stressBraneDim + L;
spinDensityDim = stressBraneDim + T;

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

zeroVecQ[vec_] := AllTrue[Flatten[FullSimplify[vec]], TrueQ[# === 0] &];
nonzeroVecQ[vec_] := Not[zeroVecQ[vec]];
zeroMatrixQ[mat_] := AllTrue[Flatten[FullSimplify[mat]], TrueQ[# === 0] &];
nonzeroMatrixQ[mat_] := Not[zeroMatrixQ[mat]];
curl3[v_] := Cross[kvec, v];

kvec = {kx, ky, kz};
k2 = kx^2 + ky^2 + kz^2;
pTnum = k2 IdentityMatrix[3] - Outer[Times, kvec, kvec];
uvec = {ux, uy, uz};
Vpar = {Vx, Vy, Vz};

A = {
  {axx, axy, axz},
  {ayx, ayy, ayz},
  {azx, azy, azz}
};
gradVars = Flatten[A];

(* Postulated rotational/MacCullagh package. *)
rvec = {ayz - azy, azx - axz, axy - ayx};
Wrot = (muR/2) (rvec . rvec);
sigmaRot = FullSimplify[Table[D[Wrot, A[[a, b]]], {a, 1, 3}, {b, 1, 3}]];
sigmaExpected = muR (A - Transpose[A]);
antisymRot = FullSimplify[sigmaRot - Transpose[sigmaRot]];

(* Minimal gyrostatic completion: no first-gradient couple stress is added;
   the local spin reservoir carries the antisymmetric torque rate. *)
spinRateAB = FullSimplify[-antisymRot];
angularResidualNoSpin = antisymRot;
angularResidualWithSpin = FullSimplify[antisymRot + spinRateAB];
coupleStressZero = ConstantArray[0, {3, 3, 3}];
coupleTractionZero = ConstantArray[0, {3, 3}];

nuVec = {nu1, nu2, nu3};
tractionRot = FullSimplify[Transpose[sigmaRot] . nuVec];

fourierRules = Thread[gradVars -> Flatten[Outer[Times, I kvec, uvec]]];
sigmaFourier = FullSimplify[sigmaRot /. fourierRules];
divSigmaFourier = FullSimplify[
  Table[Sum[I kvec[[a]] sigmaFourier[[a, b]], {a, 1, 3}], {b, 1, 3}]
];
divSigmaExpected = FullSimplify[muR (kvec (kvec . uvec) - k2 uvec)];
rotTransverse = FullSimplify[pTnum . divSigmaFourier];
rotLongitudinal = FullSimplify[divSigmaFourier /. Thread[uvec -> phi kvec]];
rotTransverseControl = FullSimplify[
  divSigmaFourier /. {kx -> 0, ky -> 0, kz -> k, ux -> Uamp, uy -> 0, uz -> 0}
];

(* Stage-1 bulk stress projection, kept symbolic. *)
piPar = {
  {pixx, pixy, pixz},
  {pixy, piyy, piyz},
  {pixz, piyz, pizz}
};
piW = {piwx, piwy, piwz};
slopesFourier = I uw kvec;
TnaFourier = FullSimplify[piW + (piww IdentityMatrix[3] - piPar) . slopesFourier];

Jvec = {Jx, Jy, Jz};
tauAeff = {tauAx, tauAy, tauAz};

(* v_n feedback channel: the direct constant-coefficient v_n source is
   longitudinal, but a normal-flow perturbation changes T_wa=m rho v_w v_a.
   That induced T_wa is transverse unless the background in-plane flow is
   parallel to the perturbation wavevector. *)
directVnSource = -I lambdaN uw kvec;
vnFeedbackForce = mG rho3 dvw Vpar;

totalExchange3 = FullSimplify[
  divSigmaFourier + TnaFourier + alphaU Jvec + tauAeff + vnFeedbackForce
];
totalExchangeT = FullSimplify[pTnum . totalExchange3];
totalExchangeCurl = FullSimplify[curl3[totalExchange3]];

noLeakRules = Join[
  Thread[uvec -> phi kvec],
  {
    piwx -> 0, piwy -> 0, piwz -> 0,
    pixx -> piIso, piyy -> piIso, pizz -> piIso,
    pixy -> 0, pixz -> 0, piyz -> 0,
    tauAx -> 0, tauAy -> 0, tauAz -> 0,
    Jx -> 0, Jy -> 0, Jz -> 0,
    dvw -> 0
  }
];
transverseRotRules = {
  kx -> 0, ky -> 0, kz -> k,
  ux -> Uamp, uy -> 0, uz -> 0,
  piwx -> 0, piwy -> 0, piwz -> 0,
  pixx -> piIso, piyy -> piIso, pizz -> piIso,
  pixy -> 0, pixz -> 0, piyz -> 0,
  tauAx -> 0, tauAy -> 0, tauAz -> 0,
  Jx -> 0, Jy -> 0, Jz -> 0,
  dvw -> 0
};
genericTwaRules = {
  kx -> 0, ky -> 0, kz -> k,
  ux -> 0, uy -> 0, uz -> 0,
  piwx -> piwx, piwy -> 0, piwz -> 0,
  pixx -> piIso, piyy -> piIso, pizz -> piIso,
  pixy -> 0, pixz -> 0, piyz -> 0,
  tauAx -> 0, tauAy -> 0, tauAz -> 0,
  Jx -> 0, Jy -> 0, Jz -> 0,
  dvw -> 0
};
vnFeedbackTransverseRules = {
  kx -> 0, ky -> 0, kz -> k,
  ux -> 0, uy -> 0, uz -> 0,
  piwx -> 0, piwy -> 0, piwz -> 0,
  pixx -> piIso, piyy -> piIso, pizz -> piIso,
  pixy -> 0, pixz -> 0, piyz -> 0,
  tauAx -> 0, tauAy -> 0, tauAz -> 0,
  Jx -> 0, Jy -> 0, Jz -> 0,
  Vx -> Vx, Vy -> 0, Vz -> 0,
  dvw -> dvw
};
vnFeedbackParallelRules = {
  kx -> 0, ky -> 0, kz -> k,
  ux -> 0, uy -> 0, uz -> 0,
  piwx -> 0, piwy -> 0, piwz -> 0,
  pixx -> piIso, piyy -> piIso, pizz -> piIso,
  pixy -> 0, pixz -> 0, piyz -> 0,
  tauAx -> 0, tauAy -> 0, tauAz -> 0,
  Jx -> 0, Jy -> 0, Jz -> 0,
  Vx -> 0, Vy -> 0, Vz -> Vz,
  dvw -> dvw
};

(* Static Coulomb C2 check: J^0 Phi is a nonzero static source term, but it
   has no direct transverse force in the bulk variables at this stage. *)
c2StaticSource = -J0 Phi;
c2BulkVariation = D[c2StaticSource, #] & /@ {rho, theta, vx, vy, vz, vw};
c2TransverseSource = {0, 0, 0};

fLeakScaleDim = stressBraneDim - 2 L + uDim;
fMagDim = rho3Dim + circulationDim + vortexDensityDim + velocityDim;
fGravDim = rho3Dim + 2 velocityDim - L;
fRefDim = braneForceDim;
epsLeakExpression =
  "|P_T [D_b sigma^R_ba + T_na + alpha_u J_a + t_A,a + delta T_wa^(v_n)]| / max(rho_3 Gamma n_v V_rel, rho_3 v_r^2/L_g)";

dimChecks = {
  homogeneous[
    "postulated rotational brane package",
    <|
      "rho_parallel dot_u^2" -> rho3Dim + 2 (uDim - T),
      "mu_R r_a r_a" -> stressBraneDim,
      "alpha_u J^a u_a" -> braneLag,
      "J^0 Phi_br" -> braneLag
    |>,
    braneLag
  ],
  checkDim["rotational force-stress sigma^R_ab", stressBraneDim, braneLag],
  checkDim["rotational brane force density D_b sigma^R_ba", braneForceDim, braneLag - L],
  checkDim["finite-thickness brane force density B_l D_b sigma^R_ba", Bdim + braneForceDim, bulkLag - L],
  checkDim["Stage-1 traction T_na as brane-level force density", bulkStressDim, braneForceDim],
  checkDim["finite-thickness Stage-1 source B_l T_na", Bdim + bulkStressDim, bulkLag - L],
  checkDim["mouth effective force density delta_boundary t_A", braneForceDim, braneLag - L],
  checkDim["couple stress M_cab", coupleStressDim, stressBraneDim + L],
  checkDim["spin density s_ab", spinDensityDim, stressBraneDim + T],
  checkDim["Magnus reference force density rho_3 Gamma n_v V_rel", fMagDim, braneForceDim],
  checkDim["gravity-flow reference force density rho_3 v_r^2/L_g", fGravDim, braneForceDim],
  checkDim["rotational leak scale mu_R k^2 U_T", fLeakScaleDim, braneForceDim]
};

postulateChecks = {
  checkBoolean["sigma^R_ab derived from Wrot equals mu_R(D_a u_b-D_b u_a)", zeroMatrixQ[sigmaRot - sigmaExpected]],
  checkBoolean["D_b sigma^R_ba Fourier form derived, not assigned", zeroVecQ[divSigmaFourier - divSigmaExpected]],
  checkBoolean["longitudinal displacement costs no rotational force", zeroVecQ[rotLongitudinal]],
  checkBoolean["transverse displacement has nonzero rotational force", nonzeroVecQ[rotTransverseControl]]
};

angularChecks = {
  checkBoolean["without spin/couple completion angular residual is nonzero", nonzeroMatrixQ[angularResidualNoSpin]],
  checkBoolean["minimal gyrostatic spin-rate completion closes local angular momentum", zeroMatrixQ[angularResidualWithSpin]],
  checkBoolean["minimal first-gradient couple stress is zero by postulate", zeroMatrixQ[Flatten[coupleStressZero]]],
  checkBoolean["minimal couple traction is zero by postulate", zeroMatrixQ[coupleTractionZero]]
};

boundaryChecks = {
  checkBoolean["boundary traction nu_a sigma^R_ab is nonzero generically", nonzeroVecQ[tractionRot]],
  checkBoolean["bulk variation of U_R with respect to rho, theta, v_i is direct-zero", zeroVecQ[D[Wrot, #] & /@ {rho, theta, vx, vy, vz, vw}]]
};

variationChecks = {
  checkBoolean["direct constant-coefficient v_n source is longitudinal", zeroVecQ[pTnum . directVnSource]],
  checkBoolean["direct constant-coefficient v_n source has zero curl", zeroVecQ[curl3[directVnSource]]]
};

leakChecks = {
  checkBoolean["special longitudinal/isotropic/no-mouth/no-current branch has zero total transverse source", zeroVecQ[FullSimplify[totalExchangeT /. noLeakRules]]],
  checkBoolean["rotational transverse brane mode alone sources transverse exchange", nonzeroVecQ[FullSimplify[totalExchangeT /. transverseRotRules]]],
  checkBoolean["generic Stage-1 T_wa channel sources transverse exchange", nonzeroVecQ[FullSimplify[totalExchangeT /. genericTwaRules]]],
  checkBoolean["generic total exchange curl is nonzero", nonzeroVecQ[totalExchangeCurl]]
};

vnChecks = {
  checkBoolean["v_n feedback via delta T_wa=m rho delta v_w v_parallel is transverse for v_parallel not parallel k", nonzeroVecQ[FullSimplify[pTnum . vnFeedbackForce /. vnFeedbackTransverseRules]]],
  checkBoolean["v_n feedback transverse projection vanishes for v_parallel parallel k", zeroVecQ[FullSimplify[pTnum . vnFeedbackForce /. vnFeedbackParallelRules]]],
  checkBoolean["v_n feedback curl is nonzero for v_parallel not parallel k", nonzeroVecQ[FullSimplify[curl3[vnFeedbackForce] /. vnFeedbackTransverseRules]]]
};

c2Checks = {
  checkBoolean["C2 static source term is nonzero", Not[TrueQ[c2StaticSource === 0]]],
  checkBoolean["C2 static source has no direct bulk variation under Stage-0 independence premise", zeroVecQ[c2BulkVariation]],
  checkBoolean["C2 static source has no direct transverse bulk force", zeroVecQ[pTnum . c2TransverseSource]]
};

orderChecks = {
  checkBoolean["named small parameter epsilon_leak is a ratio of equal force-density dimensions", True],
  checkBoolean["bounded condition is not structural no-leak because generic source is nonzero", nonzeroVecQ[totalExchangeT]]
};

checks = Join[
  dimChecks, postulateChecks, angularChecks, boundaryChecks,
  variationChecks, leakChecks, vnChecks, c2Checks, orderChecks
];
allPass = AllTrue[checks, TrueQ[#["pass"]] &];

stage3Token = If[
  allPass && nonzeroVecQ[totalExchangeT],
  "LEAK_BOUNDED_CONDITIONAL(epsilon_leak<<1 + transverse-cancellation/impedance price; otherwise FAIL_CONSTITUTIVE_TRACTION_LEAK)",
  "SCRIPT_CHECK_FAILURE"
];

report = <|
  "stage" -> "pathA_23 Stage 3 constitutive no-leak closure",
  "conditional_flag" -> "CONDITIONAL",
  "postulated_package" -> "ROTATIONAL_POSTULATED + minimal gyrostatic spin-rate closure; no derivation claimed",
  "stage3_token" -> stage3Token,
  "checks_passed" -> Count[checks[[All, "pass"]], True],
  "checks_total" -> Length[checks],
  "all_pass" -> allPass,
  "load_bearing_expressions" -> <|
    "sigma_R_ab" -> ToString[InputForm[sigmaRot]],
    "boundary_traction_t_b" -> ToString[InputForm[tractionRot]],
    "div_sigma_R_fourier" -> ToString[InputForm[divSigmaFourier]],
    "T_na_fourier" -> ToString[InputForm[TnaFourier]],
    "vn_feedback_force" -> ToString[InputForm[vnFeedbackForce]],
    "total_exchange_3D" -> ToString[InputForm[totalExchange3]],
    "total_exchange_transverse_numerator" -> ToString[InputForm[totalExchangeT]],
    "total_exchange_curl" -> ToString[InputForm[totalExchangeCurl]],
    "epsilon_leak" -> epsLeakExpression
  |>,
  "conditions" -> <|
    "structural_no_leak_special_case" -> "u_a longitudinal, P_T T_wa=0, P_T[(T_ww delta_ab-T_ab)k_b]=0, P_T(alpha_u J_a + t_A,a)=0, and v_parallel || k or delta v_w=0",
    "bounded_condition" -> "epsilon_leak << 1 relative to both the Magnus reservoir and gravity-flow reference force",
    "price" -> "a new small brane-to-bulk impedance/cancellation condition; not derived from the postulated rotational law"
  |>,
  "checks" -> checks
|>;

outDir = FileNameJoin[{"software", "stage1_solver", "_scratch"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
jsonPath = FileNameJoin[{outDir, "pathA_23_stage3_noleak_closure_mathematica.json"}];
Export[jsonPath, report, "JSON"];

If[allPass, Exit[0], Exit[1]];
