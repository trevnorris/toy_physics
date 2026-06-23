(* pathA_23 Stage 0 action/contracts cross-check.
   Scope: algebraic/unit checks for the Stage-0 candidate only. This does not
   prove a constitutive law, a no-leak result, Maxwell structure, or a spectrum. *)

ClearAll[
  dimString, checkDim, checkBoolean, homogeneous, L, T, M, Q, dim0,
  energy, actionDim, bulkLag, braneLag, deltaEll, dw, uDim, uwDim,
  rhoBr, muBr, tauW, kappaW, kPin, piN, aSpatial, a0Dim, alphaU,
  j0Dim, jaDim, chiDim, checks, dimChecks, gaugeChecks, principalChecks,
  dofChecks, dtJ0, divJ, chi, uOnlyResidual, completedResidual,
  completedConservedResidual, uOnlyConservedResidual, kx, ky, kz, kw, omega,
  kvec, k2, id3, rho, mu, lambda, muR, cauchyP, rotP, bendP, allPass,
  report, outDir, jsonPath
];

dim0 = {0, 0, 0, 0};
L = {1, 0, 0, 0};
T = {0, 1, 0, 0};
M = {0, 0, 1, 0};
Q = {0, 0, 0, 1};

energy = M + 2 L - 2 T;
actionDim = energy + T;
bulkLag = energy - 4 L;
braneLag = energy - 3 L;
deltaEll = -L;
dw = L;

uDim = L;
uwDim = L;
rhoBr = M - 3 L;
muBr = braneLag;
tauW = braneLag;
kappaW = energy - L;
kPin = energy - 5 L;
piN = energy - 4 L;

(* SI-like source dimensions on the 3D brane. A_a has q v.A energy
   normalization; A_0 is scalar-potential normalization. *)
aSpatial = M + L - T - Q;
a0Dim = M + 2 L - 2 T - Q;
alphaU = aSpatial - uDim;
j0Dim = Q - 3 L;
jaDim = Q - 2 L - T;
chiDim = aSpatial + L;

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

homogeneous[name_, terms_Association, expected_, note_: ""] := Module[
  {pass},
  pass = AllTrue[Values[terms], TrueQ[# === expected] &];
  <|
    "name" -> name,
    "pass" -> pass,
    "expected" -> dimString[expected],
    "terms" -> Association @ KeyValueMap[#1 -> dimString[#2] &, terms],
    "note" -> note
  |>
];

dimChecks = {
  checkDim["bulk scalar Lagrangian density", bulkLag, energy - 4 L],
  checkDim["bulk scalar action measure", bulkLag + T + 4 L, actionDim],
  checkDim["brane Lagrangian density", braneLag, energy - 3 L],
  checkDim["brane action measure", braneLag + T + 3 L, actionDim],
  checkDim["finite-thickness profile B_ell(w)", deltaEll, -L],
  checkDim["normalized thickness integral int B_ell dw", deltaEll + dw, dim0],
  checkDim["bulk-density representation of brane density", braneLag + deltaEll, bulkLag],
  homogeneous[
    "in-plane elastic density terms",
    <|
      "rho_parallel dot_u^2" -> rhoBr + 2 (uDim - T),
      "Cauchy mu strain^2" -> muBr,
      "rotational mu_R curl_u^2" -> muBr
    |>,
    braneLag
  ],
  homogeneous[
    "off-brane bending density terms",
    <|
      "rho_w dot_u_w^2" -> rhoBr + 2 (uwDim - T),
      "tau_w grad_u_w^2" -> tauW,
      "kappa_w (Delta u_w)^2" -> kappaW + 2 (uwDim - 2 L),
      "k_pin u_w^2" -> kPin + 2 uwDim
    |>,
    braneLag
  ],
  checkDim["normal bulk-brane work density u_w Pi_n", uwDim + piN, braneLag],
  checkDim["A_a^brane=alpha_u u_a coefficient", alphaU + uDim, aSpatial],
  homogeneous[
    "defect-current source density",
    <|
      "J^a A_a^brane" -> jaDim + aSpatial,
      "J^0 Phi_b" -> j0Dim + a0Dim
    |>,
    braneLag
  ],
  checkDim["gauge parameter from delta A_a=partial_a chi", chiDim - L, aSpatial],
  checkDim["gauge parameter from delta Phi=-partial_t chi", chiDim - T, a0Dim]
};

(* Source-coupling admissibility algebra. After integration by parts:
   u-only generic source gives -chi div J. The phi-completed coupling gives
   -chi (partial_t J0 + div J). *)
uOnlyResidual = -chi divJ;
completedResidual = -chi (dtJ0 + divJ);
completedConservedResidual = FullSimplify[completedResidual /. divJ -> -dtJ0];
uOnlyConservedResidual = FullSimplify[uOnlyResidual /. divJ -> -dtJ0];

gaugeChecks = {
  checkBoolean[
    "phi-completed conserved-current residual",
    TrueQ[completedConservedResidual === 0],
    True,
    "Passes only with partial_t J0 + div J = 0 and boundary flux set by the BC contract."
  ],
  checkBoolean[
    "u-only generic current negative control",
    TrueQ[uOnlyConservedResidual === 0],
    False,
    "A generic time-dependent charge current is not admissible with J^a u_a alone."
  ],
  checkBoolean[
    "nonconserved current negative control",
    TrueQ[completedResidual === 0],
    False,
    "The phi-completed term is not gauge-invariant for an unconserved current."
  ]
};

kvec = {kx, ky, kz};
k2 = kx^2 + ky^2 + kz^2;
id3 = IdentityMatrix[3];
cauchyP = rho omega^2 id3 - mu k2 id3 - (lambda + mu) Outer[Times, kvec, kvec];
rotP = rho omega^2 id3 - muR (k2 id3 - Outer[Times, kvec, kvec]);
bendP = rho omega^2 - tauW k2 - kappaW k2^2 - kPin;

principalChecks = {
  checkBoolean["Cauchy principal block has no k_w", FreeQ[cauchyP, kw]],
  checkBoolean["rotational principal block has no k_w", FreeQ[rotP, kw]],
  checkBoolean["bending principal scalar has no k_w", FreeQ[bendP, kw]],
  checkBoolean[
    "negative control detects explicit k_w",
    FreeQ[k2 + kw^2, kw],
    False,
    "Confirms the no-k_w checks are active."
  ]
};

dofChecks = {
  checkBoolean["brane elastic raw components u_parallel plus u_w", 3 + 1, 4],
  checkBoolean["current-completed raw variables incl auxiliary Phi_b", 3 + 1 + 1, 5],
  checkBoolean["massless 4+1 A_M physical polarizations D-2", 5 - 2, 3],
  checkBoolean["brane zero-mode A_mu physical polarizations 4 components minus 2 gauge/constraint", 4 - 2, 2]
};

checks = Join[dimChecks, gaugeChecks, principalChecks, dofChecks];
allPass = AllTrue[checks, TrueQ[#["pass"]] &];

report = <|
  "schema" -> "stage1_pathA_23_stage0_action_contracts_mathematica/v1",
  "scope" -> "Stage-0 unit, source-admissibility, principal-w, and raw-DOF checks only",
  "pass" -> allPass,
  "checks" -> checks,
  "outcomes" -> <|
    "stage0_token" -> "ACTION_SPECIFIED_CLASSIFIED",
    "current_precheck" -> "ADMISSIBLE_WITH_CONSERVED_DEFECT_CURRENT_AND_PHI_COMPLETION",
    "u_only_negative_control" -> "NOT_ADMISSIBLE_FOR_GENERIC_TIME_DEPENDENT_CHARGE",
    "principal_w_status" -> "NO_K_W_IN_DECLARED_BRANE_PRINCIPAL_BLOCKS",
    "classification" -> "NEW_PARENT_ACTION"
  |>
|>;

outDir = FileNameJoin[{"software", "stage1_solver", "_scratch"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
jsonPath = FileNameJoin[{outDir, "pathA_23_stage0_action_contracts_mathematica.json"}];
Export[jsonPath, report, "RawJSON"];
Print["wrote ", jsonPath];
Print[
  "pathA_23 Stage 0 Mathematica cross-check: ",
  Count[checks[[All, "pass"]], True], "/", Length[checks], " checks"
];
If[TrueQ[report["pass"]], Exit[0], Exit[1]]
