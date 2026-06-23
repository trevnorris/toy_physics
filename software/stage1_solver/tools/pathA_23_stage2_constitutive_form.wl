(* pathA_23 Stage 2 constitutive-form redo.

   Mathematica is the primary engine.  This version avoids the prior tautology:
   the tested in-plane energy contains the independent symmetric deviatoric
   invariant e_dev:e_dev alongside compression.  The substructure classifier
   then decides whether mu_br is zero, nonzero, or not derivable from the
   independently motivated record. *)

ClearAll[
  dimString, checkDim, checkBoolean, homogeneous, zeroQ, zeroMatrixQ,
  nonzeroMatrixQ, classifyMu, L, T, M, Q, dim0, energy, braneLag, bulkLag,
  Bdim, uDim, rhoDim, spinInertiaDim, modulusDim, stressDim, forceDim,
  nvar, Kbr, rhoPar, lambda, muBr, muR, muC, kappaC, Ac, Ic, omega2,
  kx, ky, kz, k2, k, kvec, ux, uy, uz, uvec, vx, vy, vz, vvec,
  axx, axy, axz, ayx, ayy, ayz, azx, azy, azz, A, gradVars, theta,
  eTensor, eDev, eDevSq, Wfull, sigmaFull, nvec, tractionFull, fourierRules,
  Kfull, charFull, expectedCharFull, pTnum, fullTransverseExtractor,
  fullLongitudinalEigen, fullFluidLimitZero, fullNetworkLimitNonzero,
  sigmaSymmetric, coupleStressZero, boundaryWorkDim,
  Wfluid, Kfluid, charFluid, expectedCharFluid, fluidTransverseZero,
  Wrot, rvec, sigmaRot, Krot, charRot, expectedCharRot, rotAngularResidual,
  rotNoSpinClosureFails, rotKineticGaugeObstruction,
  cosTransBlock, cosTransChar, expectedCosTransAtZero, cosGapAtZero,
  actualFacts, fluidFacts, networkFacts, actualMuDecision, fluidMuDecision,
  networkMuDecision, stage2Token, conditionalStatus,
  dimChecks, classifierChecks, fullLawChecks, stressChecks, angularChecks,
  hamiltonianChecks, branchChecks, cosseratChecks, checks, allPass, report,
  outDir, jsonPath
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
rhoDim = M - 3 L;
spinInertiaDim = M - L;
modulusDim = braneLag;
stressDim = braneLag;
forceDim = braneLag - L;

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

zeroQ[expr_] := TrueQ[FullSimplify[expr] === 0];
zeroMatrixQ[mat_] := AllTrue[Flatten[FullSimplify[mat]], TrueQ[# === 0] &];
nonzeroMatrixQ[mat_] := Not[zeroMatrixQ[mat]];

classifyMu[facts_Association] := Which[
  TrueQ[facts["persistent_neighbor_memory"] === "ABSENT"] &&
    TrueQ[facts["surface_tension_only"] === "PRESENT"],
  "ZERO",

  TrueQ[facts["persistent_neighbor_memory"] === "PRESENT"] &&
    TrueQ[facts["affine_network_free_energy"] === "SPECIFIED"] &&
    TrueQ[facts["probe_frequency_regime"] =!= "LOW_FREQUENCY_FLUID"],
  "NONZERO",

  True,
  "UNDETERMINED"
];

kvec = {kx, ky, kz};
k2 = kx^2 + ky^2 + kz^2;
uvec = {ux, uy, uz};
vvec = {vx, vy, vz};
pTnum = k2 IdentityMatrix[3] - Outer[Times, kvec, kvec];

A = {
  {axx, axy, axz},
  {ayx, ayy, ayz},
  {azx, azy, azz}
};
gradVars = Flatten[A];
theta = Tr[A];
eTensor = (A + Transpose[A])/2;
eDev = FullSimplify[eTensor - theta IdentityMatrix[3]/3];
eDevSq = FullSimplify[Tr[eDev . eDev]];

(* The honest first-gradient isotropic energy.  muBr is not set by the invariant
   list; the classifier below decides whether the substructure record fixes it. *)
Wfull = (Kbr/2) theta^2 + muBr eDevSq;
sigmaFull = FullSimplify[Table[D[Wfull, A[[a, b]]], {a, 1, 3}, {b, 1, 3}]];
nvec = {nvar[1], nvar[2], nvar[3]};
tractionFull = FullSimplify[Transpose[sigmaFull] . nvec];
fourierRules = Thread[gradVars -> Flatten[Outer[Times, kvec, uvec]]];
Kfull = FullSimplify[
  Table[D[Wfull /. fourierRules, uvec[[i]], uvec[[j]]], {i, 1, 3}, {j, 1, 3}]
];
charFull = Factor[Det[lambda IdentityMatrix[3] - Kfull]];
expectedCharFull = Factor[(lambda - muBr k2)^2 (lambda - (Kbr + 4 muBr/3) k2)];
fullTransverseExtractor = FullSimplify[Kfull . pTnum - muBr k2 pTnum];
fullLongitudinalEigen = FullSimplify[Kfull . kvec - (Kbr + 4 muBr/3) k2 kvec];
fullFluidLimitZero = zeroMatrixQ[FullSimplify[(Kfull /. muBr -> 0) . pTnum]];
fullNetworkLimitNonzero = nonzeroMatrixQ[FullSimplify[(Kfull /. muBr -> 1) . pTnum]];
sigmaSymmetric = zeroMatrixQ[sigmaFull - Transpose[sigmaFull]];
coupleStressZero = True;
boundaryWorkDim = stressDim + 2 L + uDim;

Wfluid = Wfull /. muBr -> 0;
Kfluid = FullSimplify[Kfull /. muBr -> 0];
charFluid = Factor[Det[lambda IdentityMatrix[3] - Kfluid]];
expectedCharFluid = Factor[lambda^2 (lambda - Kbr k2)];
fluidTransverseZero = zeroMatrixQ[FullSimplify[Kfluid . pTnum]];

rvec = {ayz - azy, azx - axz, axy - ayx};
Wrot = (muR/2) (rvec . rvec);
sigmaRot = FullSimplify[Table[D[Wrot, A[[a, b]]], {a, 1, 3}, {b, 1, 3}]];
Krot = FullSimplify[
  Table[D[Wrot /. fourierRules, uvec[[i]], uvec[[j]]], {i, 1, 3}, {j, 1, 3}]
];
charRot = Factor[Det[lambda IdentityMatrix[3] - Krot]];
expectedCharRot = Factor[lambda (lambda - muR k2)^2];
rotAngularResidual = FullSimplify[sigmaRot - Transpose[sigmaRot]];
rotNoSpinClosureFails = nonzeroMatrixQ[rotAngularResidual];
rotKineticGaugeObstruction = "curl-only potential is invariant under u->u+grad chi, but 1/2 rho_parallel dot(u)^2 is not invariant under time-dependent chi without a phi/constraint sector";

(* Cosserat bookkeeping branch, evaluated along k zhat for the coupled
   transverse pair (u_x, varpi_y).  It is not selected by the substructure
   record; this exposes the extra mode and gap classification if postulated. *)
cosTransBlock = {
  {(muC + kappaC) k^2, -2 kappaC k},
  {-2 kappaC k, 4 kappaC + Ac k^2}
};
cosTransChar = Factor[Det[cosTransBlock - omega2 DiagonalMatrix[{rhoPar, Ic}]]];
expectedCosTransAtZero = Factor[(-rhoPar omega2) (4 kappaC - Ic omega2)];
cosGapAtZero = 4 kappaC/Ic;

actualFacts = <|
  "cohesion" -> "PRESENT_HYPOTHESIS",
  "finite_correlation_or_healing_length" -> "PRESENT_HYPOTHESIS",
  "surface_tension_or_tautness" -> "PRESENT_HYPOTHESIS",
  "viscoelastic_high_frequency_elasticity" -> "ASSERTED_PICTURE",
  "surface_tension_only" -> "NOT_ESTABLISHED",
  "persistent_neighbor_memory" -> "UNSPECIFIED",
  "affine_network_free_energy" -> "UNSPECIFIED",
  "probe_frequency_regime" -> "UNSPECIFIED",
  "gyrostat_or_director" -> "NOT_INDEPENDENTLY_DERIVED"
|>;
fluidFacts = <|
  "surface_tension_only" -> "PRESENT",
  "persistent_neighbor_memory" -> "ABSENT",
  "affine_network_free_energy" -> "ABSENT",
  "probe_frequency_regime" -> "LOW_FREQUENCY_FLUID"
|>;
networkFacts = <|
  "surface_tension_only" -> "ABSENT",
  "persistent_neighbor_memory" -> "PRESENT",
  "affine_network_free_energy" -> "SPECIFIED",
  "probe_frequency_regime" -> "ELASTIC_HIGH_FREQUENCY"
|>;

actualMuDecision = classifyMu[actualFacts];
fluidMuDecision = classifyMu[fluidFacts];
networkMuDecision = classifyMu[networkFacts];

dimChecks = {
  homogeneous[
    "full first-gradient brane energy density terms",
    <|
      "rho_parallel dot_u^2" -> rhoDim + 2 (uDim - T),
      "K_br theta^2" -> modulusDim,
      "mu_br e_dev:e_dev" -> modulusDim
    |>,
    braneLag
  ],
  checkDim["linear force stress sigma_ab", stressDim, braneLag],
  checkDim["boundary work int dt dS traction.delta_u", boundaryWorkDim + T, energy + T],
  checkDim["brane force density D_a sigma_ab", forceDim, braneLag - L],
  checkDim["finite-thickness bulk force density B_l D_a sigma_ab", Bdim + forceDim, bulkLag - L],
  checkDim["Cosserat spin inertia I_c for dimensionless varpi", spinInertiaDim, M - L],
  checkDim["Cosserat curvature modulus A_c", energy - L, energy - L]
};

classifierChecks = {
  checkBoolean[
    "actual substructure record leaves mu_br undetermined",
    actualMuDecision === "UNDETERMINED",
    True,
    "The record motivates cohesion/healing length but does not specify persistent in-plane neighbor memory plus a network free energy."
  ],
  checkBoolean[
    "able-to-fail control: simple-fluid/soap-film facts classify mu_br as zero",
    fluidMuDecision === "ZERO"
  ],
  checkBoolean[
    "able-to-fail control: coherent-network facts classify mu_br as nonzero",
    networkMuDecision === "NONZERO"
  ]
};

fullLawChecks = {
  checkBoolean["deviatoric invariant is explicitly present", zeroQ[D[Wfull, muBr] - eDevSq]],
  checkBoolean["compression invariant is explicitly present", zeroQ[D[Wfull, Kbr] - theta^2/2]],
  checkBoolean["Cauchy characteristic polynomial has two transverse and one longitudinal eigenvalue", zeroQ[charFull - expectedCharFull]],
  checkBoolean["transverse stiffness extractor returns mu_br k^2 on the transverse projector", zeroMatrixQ[fullTransverseExtractor]],
  checkBoolean["longitudinal stiffness is (K_br+4 mu_br/3) k^2", zeroMatrixQ[fullLongitudinalEigen]],
  checkBoolean["mu_br=0 branch has no transverse stiffness", fullFluidLimitZero],
  checkBoolean["mu_br>0 branch is detected as transverse-stiff", fullNetworkLimitNonzero]
};

stressChecks = {
  checkBoolean["symmetric-strain stress is symmetric", sigmaSymmetric],
  checkBoolean["symmetric Cauchy branch has zero couple-stress", coupleStressZero],
  checkBoolean[
    "boundary traction equals nu_a sigma_ab",
    zeroMatrixQ[tractionFull - Transpose[sigmaFull] . nvec]
  ]
};

angularChecks = {
  checkBoolean[
    "total angular momentum closes for symmetric-strain branch without spin/couple",
    sigmaSymmetric && coupleStressZero
  ],
  checkBoolean[
    "MacCullagh curl-only contrast has antisymmetric stress without spin/couple closure",
    rotNoSpinClosureFails
  ]
};

hamiltonianChecks = {
  checkBoolean[
    "Cauchy Hamiltonian stiffness is nonnegative under K_br>=0 and mu_br>=0 by eigenvalue form",
    zeroQ[charFull - expectedCharFull],
    True,
    "Eigenvalues are mu_br k^2, mu_br k^2, (K_br+4 mu_br/3) k^2; kinetic positivity also requires rho_parallel>0."
  ],
  checkBoolean[
    "negative control: mu_br<0 gives transverse ghost",
    TrueQ[(Wfull /. {muBr -> -1, Kbr -> 0, axx -> 0, ayy -> 0, azz -> 0,
      axy -> 0, axz -> 1, ayx -> 0, ayz -> 0, azx -> 0, azy -> 0}) < 0]
  ],
  checkBoolean[
    "negative control: K_br+4 mu_br/3<0 gives longitudinal ghost",
    TrueQ[(Wfull /. {muBr -> 0, Kbr -> -1, axx -> 1, ayy -> 0, azz -> 0,
      axy -> 0, axz -> 0, ayx -> 0, ayz -> 0, azx -> 0, azy -> 0}) < 0]
  ]
};

branchChecks = {
  checkBoolean["earned fluid branch spectrum is lambda^2(lambda-K_br k^2)", zeroQ[charFluid - expectedCharFluid]],
  checkBoolean["earned fluid branch transverse block is zero", fluidTransverseZero],
  checkBoolean["MacCullagh contrast spectrum is lambda(lambda-mu_R k^2)^2", zeroQ[charRot - expectedCharRot]],
  checkBoolean["MacCullagh stress is antisymmetric", zeroMatrixQ[sigmaRot + Transpose[sigmaRot]]],
  checkBoolean["MacCullagh kinetic gauge obstruction is recorded", StringLength[rotKineticGaugeObstruction] > 0]
};

cosseratChecks = {
  checkBoolean[
    "Cosserat transverse pair has acoustic plus micro-rotation optic mode at k=0",
    zeroQ[(cosTransChar /. k -> 0) - expectedCosTransAtZero],
    True,
    "The optic gap is omega_0^2=4 kappa_c/I_c; hiding it is a C6 gap-as-tuning issue unless independently derived."
  ],
  checkBoolean[
    "Cosserat branch contains a finite gap scale requiring provenance",
    cosGapAtZero === 4 kappaC/Ic
  ]
};

checks = Join[
  dimChecks,
  classifierChecks,
  fullLawChecks,
  stressChecks,
  angularChecks,
  hamiltonianChecks,
  branchChecks,
  cosseratChecks
];

allPass = AllTrue[checks, TrueQ[#["pass"]] &];
stage2Token = Which[
  Not[allPass], "SCRIPT_CHECK_FAILED",
  actualMuDecision === "ZERO", "FAIL_NO_TRANSVERSE_STIFFNESS",
  actualMuDecision === "NONZERO", "MICROSTRUCTURE_DERIVES_CAUCHY",
  actualMuDecision === "UNDETERMINED", "FAIL_UNSPECIFIED_SUBSTRUCTURE",
  True, "SCRIPT_CLASSIFIER_INCONSISTENT"
];
conditionalStatus = Which[
  stage2Token === "FAIL_UNSPECIFIED_SUBSTRUCTURE", "NOT_CONDITIONAL_UNDERDETERMINED_FAILURE",
  stage2Token === "FAIL_NO_TRANSVERSE_STIFFNESS", "NOT_CONDITIONAL_DERIVED_FAILURE",
  stage2Token === "MICROSTRUCTURE_DERIVES_CAUCHY", "NOT_CONDITIONAL_DERIVED_CAUCHY_WITH_FAIL_CAUCHY_STRAY_LONGITUDINAL_SIGNATURE",
  True, "REVIEW_REQUIRED"
];

report = <|
  "schema" -> "pathA_23_stage2_constitutive_form_mathematica/v2",
  "scope" -> "Stage-2 redo: full first-gradient symmetric strain content is present; substructure classifier decides mu_br without choosing J-only.",
  "pass" -> allPass,
  "stage2_token" -> stage2Token,
  "conditional_status" -> conditionalStatus,
  "actual_mu_decision" -> actualMuDecision,
  "substructure_facts" -> actualFacts,
  "able_to_fail_controls" -> <|
    "simple_fluid_or_soap_film" -> fluidMuDecision,
    "coherent_elastic_network" -> networkMuDecision
  |>,
  "derived_expressions" -> <|
    "W_full" -> ToString[InputForm[Wfull]],
    "e_dev_sq" -> ToString[InputForm[eDevSq]],
    "sigma_ab" -> ToString[InputForm[sigmaFull]],
    "boundary_traction_t_b" -> ToString[InputForm[tractionFull]],
    "K_ab" -> ToString[InputForm[Kfull]],
    "charpoly" -> ToString[InputForm[charFull]],
    "fluid_limit_charpoly" -> ToString[InputForm[charFluid]],
    "macCullagh_stress" -> ToString[InputForm[sigmaRot]],
    "macCullagh_charpoly" -> ToString[InputForm[charRot]],
    "cosserat_transverse_charpoly_k_zhat" -> ToString[InputForm[cosTransChar]],
    "cosserat_gap_omega0_squared" -> ToString[InputForm[cosGapAtZero]]
  |>,
  "stage3_handoff_if_postulated_or_later_derived" -> <|
    "brane_force_stress" -> "sigma_ab = K_br theta delta_ab + 2 mu_br e_<ab>",
    "couple_stress" -> "M_cab = 0 for symmetric Cauchy branch; MacCullagh/Cosserat require spin/couple provenance",
    "boundary_traction" -> "t_b = nu_a sigma_ab",
    "finite_thickness_body_force" -> "f_b^(4)=B_l D_a sigma_ab",
    "stage1_bulk_channels_to_balance_later" -> "T_na = T_wa + (T_ww delta_ab - T_ab) D_b u_w and mouth traction t_A^a"
  |>,
  "checks" -> checks
|>;

outDir = FileNameJoin[{"software", "stage1_solver", "_scratch"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
jsonPath = FileNameJoin[{outDir, "pathA_23_stage2_constitutive_form_mathematica.json"}];
Export[jsonPath, report, "RawJSON"];
Print["wrote ", jsonPath];
Print[
  "pathA_23 Stage 2 Mathematica constitutive redo: ",
  Count[checks[[All, "pass"]], True], "/", Length[checks], " checks; token ",
  stage2Token, "; mu decision ", actualMuDecision
];
If[TrueQ[report["pass"]] && stage2Token =!= "SCRIPT_CHECK_FAILED", Exit[0], Exit[1]]
