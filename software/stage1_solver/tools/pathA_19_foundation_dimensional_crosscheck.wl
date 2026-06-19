(* Algebraic/dimensional cross-check for pathA_19 foundation claims.
   Scope: flux dimensions, pin-counting, dictionary/action homogeneity.
   This script intentionally does not verify prose classifications,
   boundary-condition choices, or defect-mass emergence. *)

ClearAll[
  dimString, checkDim, checkTrue, homogeneous, factorToReach,
  dim0, L, T, M, action, energy, force, velocity, rho4, rho3, psi,
  kEos, numberRate, lpsiDensity, electricField, magneticField
];

dim0 = {0, 0, 0};
L = {1, 0, 0};
T = {0, 1, 0};
M = {0, 0, 1};
action = {2, -1, 1};
energy = {2, -2, 1};
force = {1, -2, 1};
velocity = L - T;
rho4 = -4 L;
rho3 = -3 L;
psi = -2 L;
kEos = energy - 4 rho4;
numberRate = -T;
lpsiDensity = energy - 4 L;
electricField = force;
magneticField = action - 2 L;

dimString[d_] := Module[
  {labels = {"L", "T", "M"}, pairs, nonzero},
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

factorToReach[expected_, actual_] := expected - actual;

checkDim[name_, actual_, expected_, note_: ""] := <|
  "name" -> name,
  "pass" -> TrueQ[actual === expected],
  "expected" -> dimString[expected],
  "actual" -> dimString[actual],
  "factor_needed_to_reach_expected" -> dimString[factorToReach[expected, actual]],
  "note" -> note
|>;

checkTrue[name_, actual_, expected_, note_: ""] := <|
  "name" -> name,
  "pass" -> TrueQ[actual === expected],
  "expected" -> expected,
  "actual" -> actual,
  "note" -> note
|>;

homogeneous[name_, terms_Association, note_: ""] := Module[
  {values = Values[terms], expected, pass},
  expected = First[values];
  pass = AllTrue[values, TrueQ[# === expected] &];
  <|
    "name" -> name,
    "pass" -> pass,
    "expected" -> dimString[expected],
    "terms" -> Association @ KeyValueMap[#1 -> dimString[#2] &, terms],
    "note" -> note
  |>
];

Module[
  {
    checks, flags, pinMatrixLTM, pinMatrixLT, pinRelationDim,
    maxwellCoeff, maxwellCoeffLT, actionLT, energyLT, forceLT,
    massLT, kEosLT, lpsiDensityLT, electricFieldLT, magneticFieldLT,
    g4, g3, formal4DTarget, observed3DTarget, lt3DTarget,
    outDir, jsonPath, report, allPass
  },

  maxwellCoeff = lpsiDensity - 2 electricField;
  actionLT = {2, -1, 0};
  energyLT = {2, -2, 0};
  forceLT = {1, -2, 0};
  massLT = dim0;
  kEosLT = energyLT - 4 rho4;
  lpsiDensityLT = energyLT - 4 L;
  electricFieldLT = forceLT;
  magneticFieldLT = actionLT - 2 L;
  maxwellCoeffLT = lpsiDensityLT - 2 electricFieldLT;

  pinMatrixLTM = Transpose[{{1, 0, 0}, {1, -1, 0}, {2, -1, 1}, {0, 0, 1}}];
  pinMatrixLT = Transpose[{{1, 0}, {1, -1}, {2, -1}, {0, 0}}];
  pinRelationDim = L + velocity + M - action;

  g4 = {4, -2, -1};
  g3 = {3, -2, -1};
  formal4DTarget = g4 + 5 velocity - 5 L - 5 velocity;
  observed3DTarget = g3 + 5 velocity - 5 L - 5 velocity;
  lt3DTarget = {3, -2, 0} + 5 velocity - 5 L - 5 velocity;

  checks = {
    checkDim["4D bulk number flux", rho4 + velocity + 3 L, numberRate],
    checkDim["3D brane-reduced number flux", rho3 + velocity + 2 L, numberRate],
    checkDim["4D bulk volumetric flux rho^-1 J", -rho4 + numberRate, 4 L - T],
    checkDim["3D brane volumetric flux rho3^-1 J", -rho3 + numberRate, 3 L - T],
    checkDim["constituent mass flux m_GNLS J", M + numberRate, M - T],
    checkDim["conditional hbar J/c_gamma^2 conversion", action + numberRate - 2 velocity, M],
    homogeneous[
      "bulk continuity",
      <|"partial_t rho" -> rho4 - T, "div_4(rho v)" -> rho4 + velocity - L|>
    ],
    checkDim["sound speed law K rho^4 / m", kEos + 4 rho4 - M, 2 velocity],
    checkDim["enthalpy scale m c_s^2", M + 2 velocity, energy],
    checkDim["healing length sqrt(hbar^2/(m h0))", (2 action - M - energy)/2, L],
    homogeneous[
      "GNLS Lagrangian density",
      <|
        "i hbar psi partial_t psi" -> action - T + 2 psi,
        "hbar^2/m |D psi|^2" -> 2 action - M + 2 (psi - L),
        "V rho" -> energy + rho4,
        "K rho^5" -> kEos + 5 rho4
      |>
    ],
    checkDim["gauge coupling q Ai/hbar", {1, -1, 1} - action, -L],
    homogeneous[
      "Maxwell sector with c factors",
      <|
        "coeff E^2" -> maxwellCoeff + 2 electricField,
        "coeff c^2 B^2" -> maxwellCoeff + 2 velocity + 2 magneticField,
        "A0 J0" -> energy + rho4,
        "Ai Ji" -> {1, -1, 1} + rho4 + velocity
      |>
    ],
    homogeneous[
      "wall action density",
      <|
        "mu eta_dot^2" -> {-1, 0, 1} + 2 velocity,
        "T_w eta_w^2" -> force,
        "K_eta eta^2" -> {-1, -2, 1} + 2 L
      |>
    ],
    checkTrue["pin rank with L,T,M", MatrixRank[pinMatrixLTM], 3],
    checkDim["pin null relation a c_s m_GNLS / hbar", pinRelationDim, dim0],
    checkTrue["pin rank in LT representation", MatrixRank[pinMatrixLT], 2],
    homogeneous[
      "LT representation local GNLS terms",
      <|
        "i hbar psi partial_t psi" -> actionLT - T + 2 psi,
        "hbar^2/m |D psi|^2" -> 2 actionLT - massLT + 2 (psi - L),
        "V rho" -> energyLT + rho4,
        "K rho^5" -> kEosLT + 5 rho4
      |>,
      "Representation only; not an m_GNLS emergence derivation."
    ],
    homogeneous[
      "LT representation local Maxwell terms",
      <|
        "coeff E^2" -> maxwellCoeffLT + 2 electricFieldLT,
        "coeff c^2 B^2" -> maxwellCoeffLT + 2 velocity + 2 magneticFieldLT
      |>
    ],
    homogeneous[
      "LT representation local wall terms",
      <|
        "mu eta_dot^2" -> {-1, 0, 0} + 2 velocity,
        "T_w eta_w^2" -> forceLT,
        "K_eta eta^2" -> {-1, -2, 0} + 2 L
      |>
    ]
  };

  flags = {
    <|
      "name" -> "formal_4D_R_norm_target_not_dimensionless_without_conversion",
      "actual" -> dimString[formal4DTarget],
      "factor_needed_to_reach_dimensionless" -> dimString[-formal4DTarget]
    |>,
    <|
      "name" -> "observed_3D_GR_target_not_dimensionless_without_conversion",
      "actual" -> dimString[observed3DTarget],
      "factor_needed_to_reach_dimensionless" -> dimString[-observed3DTarget]
    |>,
    <|
      "name" -> "LT_R_norm_gate_fails_without_new_conversion_factor",
      "actual" -> dimString[lt3DTarget],
      "factor_needed_to_reach_dimensionless" -> dimString[-lt3DTarget]
    |>
  };

  allPass = AllTrue[checks, TrueQ[#["pass"]] &];
  report = <|
    "schema" -> "stage1_pathA_19_foundation_mathematica_crosscheck/v1",
    "scope" -> "algebraic/dimensional only",
    "pass" -> allPass,
    "checks" -> checks,
    "flags" -> flags
  |>;

  outDir = FileNameJoin[{"software", "stage1_solver", "_scratch"}];
  If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
  jsonPath = FileNameJoin[{outDir, "pathA_19_foundation_mathematica_crosscheck.json"}];
  Export[jsonPath, report, "RawJSON"];

  Print["wrote " <> jsonPath];
  If[allPass, Print["pathA_19 Mathematica algebra cross-check: PASS"], Print["pathA_19 Mathematica algebra cross-check: FAIL"]];
  If[allPass, Exit[0], Exit[1]]
]
