(* Algebraic/dimensional cross-check for pathA_20 velocity-constant claims.
   Scope: [c_s], [c_gamma], wave-dispersion dimensions, flux dimensions,
   role-catalog dimensions, and simple ratio/nozzle algebra.
   This script intentionally does not verify physical/provenance judgments:
   c_gamma/c_s closure, transonic profile existence, no-net-accretion topology,
   hbar provenance, or standing-wave ontology. *)

ClearAll[
  dimString, checkDim, checkTrue, checkExpr, homogeneous, factorToReach,
  dim0, L, T, M, action, energy, velocity, rho4, rho3, psi, kEos,
  numberRate, waveNumber, pressure, rho, lambdaGamma
];

dim0 = {0, 0, 0};
L = {1, 0, 0};
T = {0, 1, 0};
M = {0, 0, 1};
action = {2, -1, 1};
energy = {2, -2, 1};
velocity = L - T;
rho4 = -4 L;
rho3 = -3 L;
psi = -2 L;
kEos = energy - 4 rho4;
numberRate = -T;
waveNumber = -L;
pressure = energy - 4 L;

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
  "expected" -> ToString[InputForm[expected]],
  "actual" -> ToString[InputForm[actual]],
  "note" -> note
|>;

checkExpr[name_, actual_, expected_, note_: ""] := Module[
  {residual = FullSimplify[actual - expected]},
  <|
    "name" -> name,
    "pass" -> TrueQ[residual === 0],
    "expected" -> ToString[InputForm[expected]],
    "actual" -> ToString[InputForm[FullSimplify[actual]]],
    "residual" -> ToString[InputForm[residual]],
    "note" -> note
  |>
];

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
    checks, algebra, cstarOverC0, rhostarOverRho0, idealFluxFactor,
    outDir, jsonPath, report, allDimPass, allAlgPass
  },

  checks = {
    checkDim[
      "sound-speed squared law c_s^2=5*K*rho^4/m_GNLS",
      kEos + 4 rho4 - M,
      2 velocity,
      "Relative to imposed P=K*rho^5."
    ],
    checkDim[
      "sound speed c_s=sqrt(5*K*rho^4/m_GNLS)",
      (kEos + 4 rho4 - M)/2,
      velocity
    ],
    homogeneous[
      "stationary quantum-Bernoulli additive terms",
      <|
        "0.5*m_GNLS*v_b^2" -> M + 2 velocity,
        "h(rho)" -> energy,
        "V_conf" -> energy,
        "Q" -> 2 action - M - 2 L
      |>
    ],
    homogeneous[
      "bulk continuity equation with v_b",
      <|
        "partial_t rho" -> rho4 - T,
        "div_4(rho v_b)" -> rho4 + velocity - L
      |>
    ],
    checkDim[
      "Madelung background velocity v_b=(hbar/m_GNLS)*grad(theta)",
      action - M - L,
      velocity
    ],
    checkDim["photon/gauge-wave speed c_gamma", velocity, velocity],
    homogeneous[
      "massless gauge-wave dispersion omega^2=c_gamma^2*k^2",
      <|
        "omega^2" -> -2 T,
        "c_gamma^2*k^2" -> 2 velocity + 2 waveNumber
      |>
    ],
    homogeneous[
      "trapped-mode wave dispersion omega^2=c_gamma^2*(k_parallel^2+k_perp^2)",
      <|
        "omega^2" -> -2 T,
        "c_gamma^2*k_parallel^2" -> 2 velocity + 2 waveNumber,
        "c_gamma^2*k_perp^2" -> 2 velocity + 2 waveNumber
      |>
    ],
    checkDim[
      "trapped-mode group velocity d omega/d k",
      velocity + waveNumber - waveNumber,
      velocity
    ],
    checkDim["ratio c_gamma/c_s", velocity - velocity, dim0],
    checkDim["tail factor (c/c_s)^3 with c=c_gamma", 3 (velocity - velocity), dim0],
    checkDim[
      "4D-bulk candidate sonic number flux rho_* c_s,* A_3,*",
      rho4 + velocity + 3 L,
      numberRate
    ],
    checkDim[
      "3D-brane candidate sonic number flux rho_3,* c_s,* A_2,*",
      rho3 + velocity + 2 L,
      numberRate
    ],
    checkDim["background pressure P0=K*rho0^5", kEos + 5 rho4, pressure],
    checkDim["pin relation hbar=m_GNLS*c_s0*a", M + velocity + L, action],
    checkDim["healing-length relation hbar=m_GNLS*c_s0*xi_h/sqrt(2)", M + velocity + L, action],
    checkDim["circulation kappa=int v_b dl", velocity + L, action - M],
    checkDim["phase-momentum exchange p=hbar*grad(theta)", action - L, M + velocity],
    checkDim[
      "quantum pressure Q=-hbar^2/(2m)*laplacian(sqrt(rho))/sqrt(rho)",
      2 action - M - 2 L,
      energy
    ],
    checkDim["candidate mass bridge hbar*J/c_gamma^2", action + numberRate - 2 velocity, M],
    checkDim["cycle-rate bridge h*J_nu/c_gamma^2", action + numberRate - 2 velocity, M]
  };

  cstarOverC0 = Sqrt[1/3];
  rhostarOverRho0 = Sqrt[cstarOverC0];
  idealFluxFactor = FullSimplify[cstarOverC0*rhostarOverRho0];

  algebra = {
    checkExpr[
      "state dependence d ln c_s / d ln rho for n=5",
      FullSimplify[rho*D[rho^2, rho]/rho^2],
      2
    ],
    checkExpr[
      "conditional ideal no-Q/no-V sonic c_s,* / c_s0",
      cstarOverC0,
      1/Sqrt[3],
      "Not accepted as the branch verdict."
    ],
    checkExpr[
      "conditional ideal no-Q/no-V sonic rho_* / rho0",
      rhostarOverRho0,
      3^(-1/4),
      "Not accepted without the actual branch profile."
    ],
    checkExpr[
      "conditional ideal no-Q/no-V flux factor Jcrit/(rho0*c_s0*A*)",
      idealFluxFactor,
      3^(-3/4),
      "Conditional Euler-nozzle factor, not the pathA_20 flux_law_verdict."
    ],
    checkExpr[
      "tail factor with lambda_gamma=c_gamma/c_s",
      lambdaGamma^3,
      lambdaGamma^3
    ]
  };

  allDimPass = AllTrue[checks, TrueQ[#["pass"]] &];
  allAlgPass = AllTrue[algebra, TrueQ[#["pass"]] &];
  report = <|
    "schema" -> "stage1_pathA_20_velocity_constants_mathematica_crosscheck/v1",
    "scope" -> "algebraic/dimensional only",
    "pass" -> TrueQ[allDimPass && allAlgPass],
    "checks" -> checks,
    "algebraic_checks" -> algebra,
    "summary" -> <|
      "dimensional_total" -> Length[checks],
      "dimensional_passed" -> Count[checks[[All, "pass"]], True],
      "algebraic_total" -> Length[algebra],
      "algebraic_passed" -> Count[algebra[[All, "pass"]], True]
    |>
  |>;

  outDir = FileNameJoin[{"software", "stage1_solver", "_scratch"}];
  If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
  jsonPath = FileNameJoin[{outDir, "pathA_20_velocity_constants_mathematica_crosscheck.json"}];
  Export[jsonPath, report, "RawJSON"];

  Print["wrote " <> jsonPath];
  If[
    TrueQ[allDimPass && allAlgPass],
    Print["pathA_20 Mathematica algebra/dimensional cross-check: PASS"],
    Print["pathA_20 Mathematica algebra/dimensional cross-check: FAIL"]
  ];
  If[TrueQ[allDimPass && allAlgPass], Exit[0], Exit[1]]
]
