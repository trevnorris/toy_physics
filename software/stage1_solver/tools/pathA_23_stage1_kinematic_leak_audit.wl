(* pathA_23 Stage 1 kinematic/coupling leak audit.

   This redo derives the coupling sources from the Stage-0 coupling and the
   scalar-sector stress projection. It does not assume Pi_n is a scalar.
   Scope: pre-constitutive, linear in brane slope unless explicitly labelled
   O(slope^2) / deferred. *)

ClearAll[
  dimString, checkDim, checkBoolean, homogeneous, zeroVecQ, nonzeroVecQ,
  containsSymbolQ, linSlope, linBrane, L, T, M, Q, energy, bulkLag,
  braneLag, Bdim, uDim, uwDim, stressDim, velocityDim, lambdaVnDim,
  phiDim, aSpatialDim, alphaUDim, j0Dim, jaDim, kx, ky, kz, k2, kvec,
  pTnum, curl3, eps, sx, sy, sz, slopes, uw, Bl, rho, theta, mG, Pbulk,
  dpinnDrho, lambdaN, J0, Phi, alphaU, ux, uy, uz, Jx, Jy, Jz,
  tx, ty, tz, lambdaBr, muBr, muR, lambdaC, muC, kappaC, Ac, varpiGap,
  pixx, pixy, pixz, piyy, piyz, pizz, piwx, piwy, piwz, piww,
  piIso, piAniso, vx, vy, vz, vw, sigNN, Pi4, piPar, piW, n4,
  eTangents, piNNExact, piNaExact, piNNLin, piNaLin, expectedPiNNLin,
  expectedPiNaLin, fourierSlopes, tNaFourier, tNaTransverse,
  tNaCurl, isoRules, tNaIso, densityScalarVariation,
  densityEulerSource, normalWorkScalarLongitudinal, v4, vN,
  piNNHydroVelocity, dPiNNdV, dLNormaldVA, dLNormaldVALinearBrane,
  vnScalarCoupling, dVnCouplingdVA, vnVelocityInPlaneSource,
  currentCoupling, currentBulkVariations, currentBraneVariation,
  staticCoulombCoupling, staticCoulombBulkSource, sourceSet,
  candidateCoeffs, dimChecks, projectionDerivationChecks,
  variationChecks, c2Checks, candidateChecks, negativeChecks,
  stressChannelChecks, bendingChecks, checks, allPass, c2Pass,
  scalarSourcesLongitudinal, stressProjectionProvenNoLeak,
  mouthTractionDeferred, stressProjectionDeferred, vnFeedbackDeferred,
  bendingO2Deferred, knownStressLeak, knownMouthLeak, stage1Token,
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
      Map[
        If[#[[2]] === 1, #[[1]], #[[1]] <> "^" <> ToString[InputForm[#[[2]]]]] &,
        nonzero
      ],
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
containsSymbolQ[expr_, sym_] := Not[FreeQ[expr, sym]];

(* Linearize in the embedding slope s_a = D_a u_w. *)
linSlope[expr_] := Expand[
  Normal[Series[expr /. {sx -> eps sx, sy -> eps sy, sz -> eps sz}, {eps, 0, 1}]] /. eps -> 1
];

(* Linearize in the brane displacement and its slope together. *)
linBrane[expr_] := Expand[
  Normal[
    Series[
      expr /. {uw -> eps uw, sx -> eps sx, sy -> eps sy, sz -> eps sz},
      {eps, 0, 1}
    ]
  ] /. eps -> 1
];

L = {1, 0, 0, 0};
T = {0, 1, 0, 0};
M = {0, 0, 1, 0};
Q = {0, 0, 0, 1};

energy = M + 2 L - 2 T;
bulkLag = energy - 4 L;
braneLag = energy - 3 L;
Bdim = -L;

uDim = L;
uwDim = L;
stressDim = energy - 4 L;
velocityDim = L - T;
lambdaVnDim = braneLag - velocityDim;
aSpatialDim = M + L - T - Q;
phiDim = M + 2 L - 2 T - Q;
alphaUDim = aSpatialDim - uDim;
j0Dim = Q - 3 L;
jaDim = Q - 2 L - T;

kvec = {kx, ky, kz};
k2 = kx^2 + ky^2 + kz^2;
pTnum = k2 IdentityMatrix[3] - Outer[Times, kvec, kvec];
curl3[v_] := Cross[kvec, v];

slopes = {sx, sy, sz};
fourierSlopes = I uw kvec;

(* General scalar-sector matter stress representative from the kept L_psi:
   Pi_ij = m_GNLS rho v_i v_j + delta_ij P(rho) + sigma_Q,ij.
   The symbolic projection below keeps the components general, then the
   hydrodynamic velocity variation checks the convective part explicitly. *)
Pi4 = {
  {pixx, pixy, pixz, piwx},
  {pixy, piyy, piyz, piwy},
  {pixz, piyz, pizz, piwz},
  {piwx, piwy, piwz, piww}
};
piPar = Pi4[[1 ;; 3, 1 ;; 3]];
piW = {piwx, piwy, piwz};
n4 = {-sx, -sy, -sz, 1};
eTangents = Table[UnitVector[4, a] + slopes[[a]] UnitVector[4, 4], {a, 1, 3}];

piNNExact = Expand[n4 . Pi4 . n4];
piNaExact = Expand[Table[n4 . Pi4 . eTangents[[a]], {a, 1, 3}]];
piNNLin = linSlope[piNNExact];
piNaLin = linSlope /@ piNaExact;
expectedPiNNLin = piww - 2 slopes . piW;
expectedPiNaLin = piW + (piww IdentityMatrix[3] - piPar) . slopes;

tNaFourier = FullSimplify[piNaLin /. Thread[slopes -> fourierSlopes]];
tNaTransverse = FullSimplify[pTnum . tNaFourier];
tNaCurl = FullSimplify[curl3[tNaFourier]];

isoRules = {
  piwx -> 0, piwy -> 0, piwz -> 0,
  pixx -> piIso, piyy -> piIso, pizz -> piIso,
  pixy -> 0, pixz -> 0, piyz -> 0
};
tNaIso = FullSimplify[tNaFourier /. isoRules];

(* Variation of the coupling density B_l u_w Pi_nn.
   A scalar variation of Pi_nn gives a scalar source; in an Euler equation its
   force is a Fourier gradient. The factor I is kept. *)
densityScalarVariation = Bl uw dpinnDrho;
densityEulerSource = I densityScalarVariation kvec;
normalWorkScalarLongitudinal = I Bl piww uw kvec;

(* Velocity variation from the actual convective stress in Pi_nn. *)
v4 = {vx, vy, vz, vw};
vN = n4 . v4;
piNNHydroVelocity = mG rho vN^2 + Pbulk + sigNN;
dPiNNdV = linSlope /@ Table[D[piNNHydroVelocity, v4[[i]]], {i, 1, 4}];
dLNormaldVA = Bl uw dPiNNdV[[1 ;; 3]];
dLNormaldVALinearBrane = linBrane /@ dLNormaldVA;

(* Optional scalar normal-flow dependence, expressed through v_n. *)
vnScalarCoupling = Bl lambdaN vN;
dVnCouplingdVA = linSlope /@ Table[D[vnScalarCoupling, v4[[i]]], {i, 1, 3}];
vnVelocityInPlaneSource = FullSimplify[dVnCouplingdVA /. Thread[slopes -> fourierSlopes]];

(* Brane current source. At Stage 1 J_br and Phi_br are brane/throat fields,
   not derived functions of rho, theta, or v_i. *)
currentCoupling = Bl (alphaU (Jx ux + Jy uy + Jz uz) - J0 Phi);
currentBulkVariations = D[currentCoupling, #] & /@ {rho, theta, vx, vy, vz, vw};
currentBraneVariation = D[currentCoupling, #] & /@ {ux, uy, uz};
staticCoulombCoupling = -Bl J0 Phi;
staticCoulombBulkSource = {0, 0, 0};

sourceSet = <|
  "density scalar variation from B_l u_w Pi_nn" -> densityEulerSource,
  "normal-work scalar force from B_l u_w Pi_nn" -> normalWorkScalarLongitudinal,
  "optional v_n scalar tangential velocity variation" -> vnVelocityInPlaneSource,
  "direct source-current bulk variation" -> staticCoulombBulkSource
|>;

candidateCoeffs = <|
  "Cauchy/Navier" -> {lambdaBr, muBr},
  "rotational/MacCullagh" -> {muR},
  "Cosserat/micropolar" -> {lambdaC, muC, kappaC, Ac, varpiGap}
|>;

dimChecks = {
  checkDim["stress projection Pi_nn or Pi_na", stressDim, energy - 4 L],
  checkDim["normal work brane density u_w Pi_nn", uwDim + stressDim, braneLag],
  checkDim["finite-thickness normal work bulk density B_l u_w Pi_nn", Bdim + uwDim + stressDim, bulkLag],
  checkDim["finite-thickness traction B_l Pi_na", Bdim + stressDim, energy - 5 L],
  checkDim["velocity-intersection coefficient Lambda_n", lambdaVnDim + velocityDim, braneLag],
  checkDim["bulk representation of Lambda_n v_n", Bdim + lambdaVnDim + velocityDim, bulkLag],
  checkDim["A_a^br = alpha_u u_a", alphaUDim + uDim, aSpatialDim],
  homogeneous[
    "source-current brane density",
    <|"J^a A_a^br" -> jaDim + aSpatialDim, "J^0 Phi_br" -> j0Dim + phiDim|>,
    braneLag
  ]
};

projectionDerivationChecks = {
  checkBoolean[
    "derived Pi_nn projection matches T_ww - 2 s_b T_wb at O(slope)",
    TrueQ[FullSimplify[piNNLin - expectedPiNNLin] === 0],
    True
  ],
  checkBoolean[
    "derived Pi_na projection matches T_wa + (T_ww delta_ab - T_ab) s_b at O(slope)",
    zeroVecQ[piNaLin - expectedPiNaLin],
    True
  ],
  checkBoolean[
    "derived Pi_na transverse projector is generically nonzero",
    nonzeroVecQ[tNaTransverse],
    True,
    "This is why Stage 1 cannot claim no-leak without profile/traction assumptions."
  ],
  checkBoolean[
    "derived Pi_na in-plane curl is generically nonzero",
    nonzeroVecQ[tNaCurl],
    True
  ],
  checkBoolean[
    "isotropic in-plane stress with zero normal-tangent stress gives longitudinal Pi_na",
    zeroVecQ[pTnum . tNaIso],
    True
  ],
  checkBoolean[
    "isotropic in-plane stress with zero normal-tangent stress gives zero curl",
    zeroVecQ[curl3[tNaIso]],
    True
  ]
};

variationChecks = Join[
  Flatten @ KeyValueMap[
    {
      checkBoolean[
        #1 <> " has zero transverse projector",
        zeroVecQ[pTnum . #2],
        True,
        "This check applies only to scalar/normal pieces, not to the derived Pi_na traction channel."
      ],
      checkBoolean[
        #1 <> " has zero in-plane curl",
        zeroVecQ[curl3[#2]],
        True
      ]
    } &,
    sourceSet
  ],
  {
    checkBoolean[
      "B_l u_w Pi_nn convective dL/dv_a is O(u_w slope), hence absent at linear brane order",
      zeroVecQ[dLNormaldVALinearBrane],
      True
    ],
    checkBoolean[
      "source-current variation with respect to rho, theta, and v_i is zero under Stage-1 independence premise",
      AllTrue[currentBulkVariations, TrueQ[FullSimplify[#] === 0] &],
      True
    ],
    checkBoolean[
      "source-current brane variation is nonzero and proportional to J^a",
      TrueQ[currentBraneVariation === Bl alphaU {Jx, Jy, Jz}],
      True
    ]
  }
];

c2Checks = {
  checkBoolean[
    "static charge coupling is not zeroed",
    Not[TrueQ[staticCoulombCoupling === 0]],
    True,
    "The C2 pass uses -B_l J^0 Phi_br as a nonzero static/longitudinal source."
  ],
  checkBoolean[
    "static charge coupling contains J0",
    containsSymbolQ[staticCoulombCoupling, J0],
    True
  ],
  checkBoolean[
    "static Coulomb direct bulk variation has zero transverse source",
    zeroVecQ[pTnum . staticCoulombBulkSource],
    True
  ],
  checkBoolean[
    "static Coulomb direct bulk variation has zero in-plane curl",
    zeroVecQ[curl3[staticCoulombBulkSource]],
    True
  ]
};

candidateChecks = KeyValueMap[
  Function[
    {candidate, coeffs},
    checkBoolean[
      "derived coupling/projection sources independent of " <> candidate <> " constitutive coefficients",
      AllTrue[
        Join[Values[sourceSet], {tNaFourier}],
        Function[src, FreeQ[src, Alternatives @@ coeffs]]
      ],
      True,
      "Stage 1 uses S_cpl plus scalar bulk stress only; constitutive brane traction/couple-stress remains later-stage."
    ]
  ],
  candidateCoeffs
];

negativeChecks = {
  checkBoolean[
    "negative control: arbitrary tangential source is detected by transverse projector",
    nonzeroVecQ[pTnum . {tx, ty, tz}],
    True
  ],
  checkBoolean[
    "negative control: arbitrary tangential source is detected by curl",
    nonzeroVecQ[curl3[{tx, ty, tz}]],
    True
  ],
  checkBoolean[
    "negative control: forbidden direct J^a-to-bulk-velocity coupling is detected by transverse projector",
    nonzeroVecQ[pTnum . (alphaU {Jx, Jy, Jz})],
    True
  ],
  checkBoolean[
    "negative control: forbidden direct J^a-to-bulk-velocity coupling is detected by curl",
    nonzeroVecQ[curl3[alphaU {Jx, Jy, Jz}]],
    True
  ]
};

stressChannelChecks = {
  checkBoolean[
    "reachable FAIL condition: known generic T_wa would trigger transverse-source failure",
    nonzeroVecQ[tNaTransverse /. {
      piwx -> tx, piwy -> ty, piwz -> tz,
      pixx -> 0, piyy -> 0, pizz -> 0, pixy -> 0, pixz -> 0, piyz -> 0
    }],
    True
  ],
  checkBoolean[
    "reachable FAIL condition: anisotropic in-plane stress would trigger transverse-source failure",
    nonzeroVecQ[tNaTransverse /. {
      piwx -> 0, piwy -> 0, piwz -> 0,
      pixx -> piAniso, piyy -> 0, pizz -> 0, pixy -> 0, pixz -> 0, piyz -> 0
    }],
    True
  ]
};

bendingChecks = {
  checkBoolean[
    "u_w slope is retained with Fourier factor I",
    containsSymbolQ[fourierSlopes, uw] && containsSymbolQ[fourierSlopes, I],
    True
  ],
  checkBoolean[
    "normal scalar bending source is longitudinal when Pi_nn is scalar/isotropic",
    zeroVecQ[pTnum . normalWorkScalarLongitudinal],
    True
  ],
  checkBoolean[
    "optional v_n scalar source is longitudinal at constant Lambda_n",
    zeroVecQ[pTnum . vnVelocityInPlaneSource],
    True
  ],
  checkBoolean[
    "optional v_n scalar source has zero curl at constant Lambda_n",
    zeroVecQ[curl3[vnVelocityInPlaneSource]],
    True
  ]
};

checks = Join[
  dimChecks,
  projectionDerivationChecks,
  variationChecks,
  c2Checks,
  candidateChecks,
  negativeChecks,
  stressChannelChecks,
  bendingChecks
];

allPass = AllTrue[checks, TrueQ[#["pass"]] &];
c2Pass = Not[TrueQ[staticCoulombCoupling === 0]] &&
  zeroVecQ[pTnum . staticCoulombBulkSource] &&
  zeroVecQ[curl3[staticCoulombBulkSource]];
scalarSourcesLongitudinal = AllTrue[Values[sourceSet], zeroVecQ[pTnum . #] && zeroVecQ[curl3[#]] &];
stressProjectionProvenNoLeak = zeroVecQ[tNaTransverse] && zeroVecQ[tNaCurl];

mouthTractionDeferred = True;
stressProjectionDeferred = Not[stressProjectionProvenNoLeak];
vnFeedbackDeferred = True;
bendingO2Deferred = True;
knownStressLeak = False;
knownMouthLeak = False;

stage1Token = Which[
  Not[allPass], "SCRIPT_CHECK_FAILED",
  knownStressLeak || knownMouthLeak, "FAIL_COUPLING_SOURCES_TRANSVERSE",
  Not[c2Pass], "FAIL_COUPLING_SOURCES_TRANSVERSE",
  Not[scalarSourcesLongitudinal], "FAIL_COUPLING_LEAKS_INDEPENDENT_OF_CONSTITUTIVE",
  stressProjectionDeferred || mouthTractionDeferred || vnFeedbackDeferred || bendingO2Deferred,
    "LEAK_CONDITIONS_DEFERRED",
  True, "NO_KINEMATIC_LEAK_FOR_ALL_CANDIDATES"
];

report = <|
  "schema" -> "pathA_23_stage1_kinematic_leak_audit_mathematica/v2",
  "scope" -> "Stage-1 pre-constitutive coupling leak audit; stress projection derived; linear in brane slope",
  "pass" -> allPass,
  "stage1_token" -> stage1Token,
  "derived_sources" -> <|
    "Pi_nn_linear" -> ToString[InputForm[piNNLin]],
    "Pi_na_linear" -> ToString[InputForm[piNaLin]],
    "Pi_na_fourier" -> ToString[InputForm[tNaFourier]],
    "Pi_na_transverse_projector" -> ToString[InputForm[tNaTransverse]],
    "Pi_na_curl" -> ToString[InputForm[tNaCurl]],
    "density_scalar_euler_source" -> ToString[InputForm[densityEulerSource]],
    "v_n_tangential_velocity_source" -> ToString[InputForm[vnVelocityInPlaneSource]],
    "source_current_bulk_variations" -> ToString[InputForm[currentBulkVariations]],
    "source_current_brane_variation" -> ToString[InputForm[currentBraneVariation]],
    "static_coulomb_coupling" -> ToString[InputForm[staticCoulombCoupling]]
  |>,
  "no_leak_conditions" -> <|
    "stress_projection" ->
      "P_T^{ab} Pi_wb = 0 and P_T^{ab} Pi_bc k_c = 0, equivalently curl(k, Pi_w + i u_w (Pi_ww I - Pi_parallel) k)=0",
    "mouth" -> "P_T^{ab} t_A,b = 0 and curl(k,t_A)=0 at each mouth",
    "v_n_feedback" -> "normal-flow perturbation must not generate bulk vorticity through the Euler/vorticity evolution"
  |>,
  "deferred_quantities" -> {
    "actual scalar stress components Pi_wa and anisotropic Pi_ab on the selected throat/brane branch",
    "mouth tangential traction t_A^a from S_mouth",
    "O(slope^2) bending/curvature corrections",
    "v_n -> bulk-vorticity feedback through the scalar Euler/vorticity equations"
  },
  "outcomes" -> <|
    "scalar_sources_longitudinal" -> scalarSourcesLongitudinal,
    "stress_projection_proven_no_leak" -> stressProjectionProvenNoLeak,
    "stress_projection_deferred" -> stressProjectionDeferred,
    "mouth_traction_deferred" -> mouthTractionDeferred,
    "c2_sector_decomposition" -> c2Pass,
    "reachable_fail_token_for_known_transverse_stress_or_mouth" -> "FAIL_COUPLING_SOURCES_TRANSVERSE"
  |>,
  "checks" -> checks
|>;

outDir = FileNameJoin[{"software", "stage1_solver", "_scratch"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
jsonPath = FileNameJoin[{outDir, "pathA_23_stage1_kinematic_leak_audit_mathematica.json"}];
Export[jsonPath, report, "RawJSON"];
Print["wrote ", jsonPath];
Print[
  "pathA_23 Stage 1 Mathematica leak audit: ",
  Count[checks[[All, "pass"]], True], "/", Length[checks], " checks; token ",
  stage1Token
];
If[TrueQ[report["pass"]] && stage1Token =!= "SCRIPT_CHECK_FAILED", Exit[0], Exit[1]]
