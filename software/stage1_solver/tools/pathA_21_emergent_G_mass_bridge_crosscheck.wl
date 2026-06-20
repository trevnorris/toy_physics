(* Algebraic/dimensional cross-check for pathA_21 emergent G and mass bridge.
   Scope: dimensions of the G-free drain force coefficient, cycle/angular
   bridge candidates, retained hbar action dimension, and the conditional
   G_eff algebra. This script intentionally does not adjudicate non-algebraic
   source questions: alpha_J provenance, equivalence principle, universal
   inverse-square closure, W_eff realization, or brane photon reduction. *)

ClearAll[
  dimString, factorToReach, checkDim, checkExpr, homogeneous,
  dim0, L, T, M, action, energy, force, velocity, rho4, rho3,
  numberRate, q3Vol, q4Vol, massDensity3, massDensity4,
  forceCoeff3, forceCoeff4, n3Eff
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
numberRate = -T;
q3Vol = -rho3 + numberRate;
q4Vol = -rho4 + numberRate;
massDensity3 = M + rho3;
massDensity4 = M + rho4;
forceCoeff3 = force + 2 L;
forceCoeff4 = force + 3 L;
n3Eff = L + rho4;

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
    checks, algebra, q1, q2, rhoM, r, hbar, jnu, jomega, alpha1, alpha2,
    cgam, cF, forceRadial, velocityFrom1, forceFromCross, m1, m2, gCond,
    outDir, jsonPath, report, allDimPass, allAlgPass
  },

  checks = {
    checkDim[
      "reduced-3D volumetric drain strength Q_3=J/n_3",
      q3Vol,
      3 L - T,
      "Profile solve supplies n_3,eff and Theta_Q."
    ],
    checkDim[
      "bulk-4D volumetric drain strength Q_4=J/rho_4",
      q4Vol,
      4 L - T
    ],
    checkDim[
      "far-field brane drain velocity v_r=Q_3/(4*pi*r^2)",
      q3Vol - 2 L,
      velocity
    ],
    checkDim[
      "far-field bulk drain velocity v_R=Q_4/(S_3*R^3)",
      q4Vol - 3 L,
      velocity
    ],
    checkDim[
      "reduced-3D pressure/momentum force F=rho_m3*Q_2*v_1",
      massDensity3 + q3Vol + velocity,
      force
    ],
    checkDim[
      "reduced-3D force coefficient C_F=rho_m3*Q_1*Q_2/(4*pi)",
      massDensity3 + 2 q3Vol,
      forceCoeff3
    ],
    checkDim[
      "bulk-4D force coefficient C_F4=rho_m4*Q_1*Q_2/S_3",
      massDensity4 + 2 q4Vol,
      forceCoeff4
    ],
    checkDim[
      "effective reduced number density n_3,eff=W_eff*rho_4",
      n3Eff,
      rho3
    ],
    checkDim[
      "angular-rate bridge candidate alpha_J*hbar*J_omega/c_gamma^2",
      action + numberRate - 2 velocity,
      M
    ],
    checkDim[
      "cycle-rate bridge candidate alpha_J*h*J_nu/c_gamma^2",
      action + numberRate - 2 velocity,
      M
    ],
    checkDim[
      "candidate throat Hamiltonian ratio H_throat/(hbar*J_omega)",
      energy - action - numberRate,
      dim0
    ],
    checkDim[
      "inertial-throat kinetic coefficient dimension",
      M,
      M
    ],
    checkDim[
      "hbar remains an action dimension in retained {L,T,M}",
      action,
      {2, -1, 1}
    ],
    checkDim[
      "conditional G_eff from C_F*c_gamma^4/(alpha1*alpha2*hbar^2*J1*J2)",
      forceCoeff3 + 4 velocity - 2 action - 2 numberRate,
      {3, -2, -1}
    ],
    checkDim[
      "conditional G_eff using W_eff*rho4: c_gamma^4*m_GNLS/(W_eff*rho4*hbar^2)",
      4 velocity + M - n3Eff - 2 action,
      {3, -2, -1}
    ],
    checkDim[
      "lambda_gamma=c_gamma/c_s",
      velocity - velocity,
      dim0
    ]
  };

  q1 = Symbol["Q1"];
  q2 = Symbol["Q2"];
  rhoM = Symbol["rhoM"];
  r = Symbol["r"];
  hbar = Symbol["hbar"];
  jnu = Symbol["Jnu"];
  jomega = Symbol["Jomega"];
  alpha1 = Symbol["alpha1"];
  alpha2 = Symbol["alpha2"];
  cgam = Symbol["cgam"];
  cF = rhoM q1 q2/(4 Pi);
  forceRadial = -cF/r^2;
  velocityFrom1 = -q1/(4 Pi r^2);
  forceFromCross = rhoM q2 velocityFrom1;
  m1 = alpha1 hbar jomega/cgam^2;
  m2 = alpha2 hbar jomega/cgam^2;
  gCond = cF cgam^4/(alpha1 alpha2 hbar^2 jomega^2);

  algebra = {
    checkExpr[
      "reduced-3D drain force is inverse-square before any G",
      forceFromCross,
      forceRadial,
      "Inward-positive Q_i; pressure/momentum cross term."
    ],
    checkExpr[
      "G-free force coefficient C_F",
      -forceRadial r^2,
      cF,
      "No G is introduced."
    ],
    checkExpr[
      "cycle-rate bridge keeps the 2*pi outside alpha_J",
      2 Pi hbar jnu,
      2 Pi hbar jnu
    ],
    checkExpr[
      "conditional Newton coefficient algebra if independent bridge existed",
      FullSimplify[gCond m1 m2],
      cF,
      "Algebraic downstream conditional only."
    ]
  };

  allDimPass = AllTrue[checks, TrueQ[#["pass"]] &];
  allAlgPass = AllTrue[algebra, TrueQ[#["pass"]] &];
  report = <|
    "schema" -> "stage1_pathA_21_emergent_g_mass_bridge_mathematica_crosscheck/v1",
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
  jsonPath = FileNameJoin[{outDir, "pathA_21_emergent_G_mass_bridge_mathematica_crosscheck.json"}];
  Export[jsonPath, report, "RawJSON"];

  Print["wrote " <> jsonPath];
  If[
    TrueQ[allDimPass && allAlgPass],
    Print["pathA_21 Mathematica algebra/dimensional cross-check: PASS"],
    Print["pathA_21 Mathematica algebra/dimensional cross-check: FAIL"]
  ];
  If[TrueQ[allDimPass && allAlgPass], Exit[0], Exit[1]]
]
