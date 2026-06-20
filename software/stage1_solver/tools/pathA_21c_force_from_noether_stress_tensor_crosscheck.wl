(* Algebraic/dimensional cross-check for pathA_21c Noether stress force.
   Scope: dimensions of the momentum-balance law, EOS/quantum identities,
   and the explicit angular traction integral that yields the reduced-3D
   Q1 Q2/r^2 structure and matter-stress sign. This script does not prove
   source anchoring, profile selection, V_conf body-force closure, or Maxwell
   cancellation; those remain prose/residual claims in the report. *)

ClearAll[
  dimString, factorToReach, checkDim, checkExpr, homogeneous,
  dim0, L, T, M, action, energy, force, velocity, rho4, rho3, numberRate,
  q3Vol, q4Vol, massDensity3, massDensity4, stress3, stress4,
  forceDensity3, forceDensity4, energyGradient, quantumPrefactor
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
forceDensity3 = force - 3 L;
forceDensity4 = force - 4 L;
energyGradient = energy - L;
quantumPrefactor = 2 action - M - 2 L;

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
    checks, algebra, mGnls, n3, rhoInf4, q1, q2, r12, radius4, v1,
    kEos, rhoSym, vdot, hExpr, pExpr, dhDrho, dPDrho,
    deltaRhoCross, deltaPressureCross, hbar, x, m, rhoFn,
    quantumPotential, sigmaQ, quantumResidual, omega2, omega3, d3, d4,
    convectiveFlux3, pressureFlux3, totalFlux3, forceAlongV13,
    v1FromGauss3, forceAlongRhat3, convectiveFlux4, pressureFlux4,
    totalFlux4, forceAlongV14, v1FromGauss4, forceAlongRhat4,
    outDir, jsonPath, report, allDimPass, allAlgPass
  },

  checks = {
    homogeneous[
      "reduced-3D momentum-balance terms",
      <|
        "partial_t(m*N*v_i)" -> massDensity3 + velocity - T,
        "partial_j Pi_ij" -> stress3 - L,
        "V_conf body force" -> rho3 + energyGradient
      |>
    ],
    homogeneous[
      "bulk-4D momentum-balance terms",
      <|
        "partial_t(m*rho*v_i)" -> massDensity4 + velocity - T,
        "partial_J Pi_iJ" -> stress4 - L,
        "V_conf body force" -> rho4 + energyGradient
      |>
    ],
    homogeneous[
      "reduced-3D Noether stress representatives",
      <|
        "convective" -> massDensity3 + 2 velocity,
        "pressure" -> (energy - 4 L) + L,
        "quantum" -> quantumPrefactor + rho3
      |>,
      "Compact lane pressure is reduced by one transverse length."
    ],
    homogeneous[
      "Euler force-per-volume identity terms in reduced 3D",
      <|
        "m*N*acceleration" -> massDensity3 + velocity - T,
        "N*grad h" -> rho3 + energyGradient,
        "N*grad Q" -> rho3 + energyGradient,
        "N*grad V_conf" -> rho3 + energyGradient,
        "N*q(E+vB)" -> rho3 + energyGradient
      |>
    ],
    checkDim["reduced-3D surface traction integral gives force", stress3 + 2 L, force],
    checkDim["bulk-4D surface traction integral gives force", stress4 + 3 L, force],
    checkDim[
      "reduced-3D Noether force coefficient m*N*Q1*Q2",
      massDensity3 + 2 q3Vol,
      force + 2 L
    ],
    checkDim[
      "bulk-4D Noether force coefficient m*rho4*Q1*Q2",
      massDensity4 + 2 q4Vol,
      force + 3 L
    ],
    checkDim["V_conf body-force volume term dimension", forceDensity3 + 3 L, force]
  };

  mGnls = Symbol["mGNLS"];
  n3 = Symbol["N3"];
  rhoInf4 = Symbol["rhoInf4"];
  q1 = Symbol["Q1"];
  q2 = Symbol["Q2"];
  r12 = Symbol["r12"];
  radius4 = Symbol["R"];
  v1 = Symbol["v1"];
  kEos = Symbol["K"];
  rhoSym = Symbol["rho"];
  vdot = Symbol["vdot"];
  hbar = Symbol["hbar"];
  x = Symbol["x"];
  m = Symbol["m"];
  rhoFn = rho[x];
  omega2 = 4 Pi;
  omega3 = 2 Pi^2;
  d3 = 3;
  d4 = 4;

  hExpr = (5/4) kEos rhoSym^4;
  pExpr = kEos rhoSym^5;
  dhDrho = D[hExpr, rhoSym];
  dPDrho = D[pExpr, rhoSym];
  deltaRhoCross = -mGnls vdot/dhDrho;
  deltaPressureCross = FullSimplify[dPDrho deltaRhoCross];

  quantumPotential = -(hbar^2/(2 m)) D[Sqrt[rhoFn], {x, 2}]/Sqrt[rhoFn];
  sigmaQ = (hbar^2/(4 m)) ((D[rhoFn, x]^2)/rhoFn - D[rhoFn, {x, 2}]);
  quantumResidual = FullSimplify[D[sigmaQ, x] - rhoFn D[quantumPotential, x]];

  convectiveFlux3 = -(1 + 1/d3) mGnls n3 q2 v1;
  pressureFlux3 = (1/d3) mGnls n3 q2 v1;
  totalFlux3 = FullSimplify[convectiveFlux3 + pressureFlux3];
  forceAlongV13 = -totalFlux3;
  v1FromGauss3 = -q1/(omega2 r12^2);
  forceAlongRhat3 = FullSimplify[forceAlongV13 /. v1 -> v1FromGauss3];

  convectiveFlux4 = -(1 + 1/d4) mGnls rhoInf4 q2 v1;
  pressureFlux4 = (1/d4) mGnls rhoInf4 q2 v1;
  totalFlux4 = FullSimplify[convectiveFlux4 + pressureFlux4];
  forceAlongV14 = -totalFlux4;
  v1FromGauss4 = -q1/(omega3 radius4^3);
  forceAlongRhat4 = FullSimplify[forceAlongV14 /. v1 -> v1FromGauss4];

  algebra = {
    checkExpr[
      "EOS identity dP/drho = rho*dh/drho",
      dPDrho,
      rhoSym dhDrho
    ],
    checkExpr[
      "Bernoulli/EOS pressure cross term",
      deltaPressureCross,
      -mGnls rhoSym vdot
    ],
    checkExpr[
      "Madelung quantum stress divergence",
      quantumResidual,
      0
    ],
    checkExpr[
      "reduced-3D convective angular traction factor",
      convectiveFlux3,
      -(4/3) mGnls n3 q2 v1
    ],
    checkExpr[
      "reduced-3D pressure angular traction factor",
      pressureFlux3,
      (1/3) mGnls n3 q2 v1
    ],
    checkExpr[
      "reduced-3D total flux from Noether stress",
      totalFlux3,
      -mGnls n3 q2 v1
    ],
    checkExpr[
      "reduced-3D force structure after Gauss substitution",
      forceAlongRhat3,
      -mGnls n3 q1 q2/(4 Pi r12^2)
    ],
    checkExpr[
      "bulk-4D total flux from Noether stress",
      totalFlux4,
      -mGnls rhoInf4 q2 v1
    ],
    checkExpr[
      "bulk-4D force structure after Gauss substitution",
      forceAlongRhat4,
      -mGnls rhoInf4 q1 q2/(2 Pi^2 radius4^3)
    ]
  };

  allDimPass = AllTrue[checks, TrueQ[#["pass"]] &];
  allAlgPass = AllTrue[algebra, TrueQ[#["pass"]] &];
  report = <|
    "schema" -> "stage1_pathA_21c_force_from_noether_stress_tensor_mathematica_crosscheck/v1",
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
  jsonPath = FileNameJoin[{outDir, "pathA_21c_force_from_noether_stress_tensor_mathematica_crosscheck.json"}];
  Export[jsonPath, report, "RawJSON"];
  Print["wrote ", jsonPath];
  Print[
    "pathA_21c Mathematica cross-check: ",
    report["summary", "dimensional_passed"], "/", report["summary", "dimensional_total"],
    " dimensional; ",
    report["summary", "algebraic_passed"], "/", report["summary", "algebraic_total"],
    " algebraic"
  ];
  If[TrueQ[report["pass"]], Exit[0], Exit[1]]
]
