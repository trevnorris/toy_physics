(* Algebraic/dimensional cross-check for pathA_21b force closure and BVP.
   Scope: Gauss-solved drain powers, Pi_cross force coefficient dimensions,
   stationary GNLS term homogeneity, EOS derivative algebra, and projection
   dimensions. This script does not prove source anchoring, profile selection,
   or attractive sign; those are documented as source-chain/residual claims. *)

ClearAll[
  dimString, factorToReach, checkDim, checkExpr, homogeneous,
  dim0, L, T, M, action, energy, force, velocity, rho4, rho3, numberRate,
  q3Vol, q4Vol, massDensity3, massDensity4, stress3, stress4,
  forceCoeff3, forceCoeff4, mu, a0
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
q3Vol = numberRate - rho3;
q4Vol = numberRate - rho4;
massDensity3 = M + rho3;
massDensity4 = M + rho4;
stress3 = force - 2 L;
stress4 = force - 3 L;
forceCoeff3 = force + 2 L;
forceCoeff4 = force + 3 L;
mu = energy;
a0 = energy;

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
    checks, algebra, j, rho, r, radius4, v, kEos, rhoSym, mGnls, deltaV2,
    omega2, omega3, v3Solution, v4Solution, q1, q2, n3, r12,
    v1FromGauss, forceFromSurface, forceExpected, hExpr, dhDrho,
    densityResponse, outDir, jsonPath, report, allDimPass, allAlgPass
  },

  checks = {
    checkDim[
      "reduced-3D Gauss solve velocity J/(N_infty,3*Omega2*r^2)",
      numberRate - rho3 - 2 L,
      velocity
    ],
    checkDim[
      "bulk-4D Gauss solve velocity J/(rho_infty,4*Omega3*R^3)",
      numberRate - rho4 - 3 L,
      velocity
    ],
    checkDim[
      "reduced-3D momentum-flux stress m_GNLS*N_infty,3*v^2",
      massDensity3 + 2 velocity,
      stress3
    ],
    checkDim[
      "bulk-4D momentum-flux stress m_GNLS*rho_infty,4*v^2",
      massDensity4 + 2 velocity,
      stress4
    ],
    checkDim[
      "reduced-3D force coefficient m_GNLS*N_infty,3*Q1*Q2",
      massDensity3 + 2 q3Vol,
      forceCoeff3
    ],
    checkDim[
      "bulk-4D force coefficient m_GNLS*rho_infty,4*Q1*Q2",
      massDensity4 + 2 q4Vol,
      forceCoeff4
    ],
    homogeneous[
      "stationary GNLS eigenvalue equation terms",
      <|
        "kinetic" -> 2 action - M - 2 L,
        "V_conf" -> energy,
        "h(rho0)" -> energy,
        "q*A_00" -> a0,
        "mu" -> mu
      |>
    ],
    homogeneous[
      "stationary GNLS density-weighted Euler/Bernoulli terms",
      <|
        "m_GNLS*v^2" -> M + 2 velocity,
        "h(rho0)" -> energy,
        "Q(rho0)" -> energy,
        "V_conf" -> energy,
        "q*A_00" -> a0
      |>
    ],
    homogeneous[
      "Pi_cross reduced-3D stress terms",
      <|
        "convective" -> massDensity3 + 2 velocity,
        "pressure" -> energy - 4 L + L,
        "quantum" -> rho3 + energy,
        "confinement" -> rho3 + energy
      |>,
      "Reduced pressure includes the transverse-width reduction in the 3D lane."
    ],
    checkDim[
      "N_infty,3 from transverse integration of rho_infty,4",
      L + rho4,
      rho3
    ]
  };

  j = Symbol["J"];
  rho = Symbol["rho"];
  r = Symbol["r"];
  radius4 = Symbol["R"];
  v = Symbol["v"];
  kEos = Symbol["K"];
  rhoSym = Symbol["rhoSym"];
  mGnls = Symbol["mGNLS"];
  deltaV2 = Symbol["deltaV2"];
  omega2 = 4 Pi;
  omega3 = 2 Pi^2;
  v3Solution = First[v /. Solve[rho omega2 r^2 v == -j, v]];
  v4Solution = First[v /. Solve[rho omega3 radius4^3 v == -j, v]];
  q1 = Symbol["Q1"];
  q2 = Symbol["Q2"];
  n3 = Symbol["N3"];
  r12 = Symbol["r12"];
  v1FromGauss = v3Solution /. {j -> q1 n3, rho -> n3, r -> r12};
  forceFromSurface = mGnls n3 q2 v1FromGauss;
  forceExpected = -mGnls n3 q1 q2/(4 Pi r12^2);
  hExpr = (5/4) kEos rhoSym^4;
  dhDrho = D[hExpr, rhoSym];
  densityResponse = FullSimplify[(-mGnls deltaV2/2)/dhDrho];

  algebra = {
    checkExpr[
      "Gauss solve gives reduced-3D drain velocity with r^(-2) power",
      v3Solution,
      -j/(4 Pi rho r^2),
      "Solved from rho*v*Omega2*r^2=-J."
    ],
    checkExpr[
      "Gauss solve gives bulk-4D drain velocity with R^(-3) power",
      v4Solution,
      -j/(2 Pi^2 rho radius4^3),
      "Solved from rho*v*Omega3*R^3=-J."
    ],
    checkExpr[
      "Pi_cross surface impulse gives reduced force coefficient",
      forceFromSurface,
      forceExpected,
      "Uses the Gauss-solved v1 and drain-2 momentum uptake."
    ],
    checkExpr[
      "EOS enthalpy derivative used for pressure-response sign",
      dhDrho,
      5 kEos rhoSym^3,
      "Sign uses stable K>0 and rho>0 in the report, not Mathematica assumptions."
    ],
    checkExpr[
      "Bernoulli density response to increased speed",
      densityResponse,
      -mGnls deltaV2/(10 kEos rhoSym^3),
      "For stable K>0, rho>0, deltaV2>0 this is negative."
    ]
  };

  allDimPass = AllTrue[checks, TrueQ[#["pass"]] &];
  allAlgPass = AllTrue[algebra, TrueQ[#["pass"]] &];
  report = <|
    "schema" -> "stage1_pathA_21b_force_closure_and_profile_bvp_mathematica_crosscheck/v1",
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
  jsonPath = FileNameJoin[{outDir, "pathA_21b_force_closure_and_profile_bvp_mathematica_crosscheck.json"}];
  Export[jsonPath, report, "RawJSON"];

  Print["wrote " <> jsonPath];
  If[
    TrueQ[allDimPass && allAlgPass],
    Print["pathA_21b Mathematica algebra/dimensional cross-check: PASS"],
    Print["pathA_21b Mathematica algebra/dimensional cross-check: FAIL"]
  ];
  If[TrueQ[allDimPass && allAlgPass], Exit[0], Exit[1]]
]
