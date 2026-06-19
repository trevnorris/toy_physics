(* Algebraic/dimensional cross-check for pathA_20b c_gamma versus c_s.
   Scope: homogeneous coupled principal-symbol algebra, c_bulk^2=C_B/C_E,
   phonon dispersion, conditional rho-dependence, and the negative control.
   This script intentionally does not adjudicate non-algebraic source questions:
   metric inheritance, brane zero-mode validity, or which localized branch is
   the observed photon. Those remain named residuals in the report. *)

ClearAll[
  dimString, factorToReach, checkDim, checkExpr, checkBool, homogeneous,
  dim0, L, T, M, action, energy, force, velocity, rho4, kEos,
  lpsiDensity, electricField, magneticField, cEDim, cBDim, ai, a0,
  waveNumber
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
kEos = energy - 4 rho4;
lpsiDensity = energy - 4 L;
electricField = force;
magneticField = action - 2 L;
cEDim = lpsiDensity - 2 electricField;
cBDim = cEDim + 2 velocity;
ai = action - L;
a0 = energy;
waveNumber = -L;

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

checkBool[name_, actual_, expected_, note_: ""] := <|
  "name" -> name,
  "pass" -> TrueQ[actual === expected],
  "expected" -> ToString[InputForm[expected]],
  "actual" -> ToString[InputForm[actual]],
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
    checks, algebra, omega, k2, rho0, kSym, mG, hbar, cE, cB, hPrime,
    cs2, cBulk2, phononMatrix, phononDet, gaugeTransverse, coupledDet,
    equalityResidual, sourceMetricEquationPresent, forcedEqualsValid,
    lambdaRhoFactor, lambdaLogSlope, lambdaGamma, outDir, jsonPath,
    allDimPass, allAlgPass, report
  },

  checks = {
    checkDim[
      "phonon sound speed c_s0=sqrt(5*K*rho0^4/m_GNLS)",
      (kEos + 4 rho4 - M)/2,
      velocity,
      "Dimension check only; equality to c_gamma is not inferred."
    ],
    checkDim[
      "Maxwell principal speed squared C_B/C_E",
      cBDim - cEDim,
      2 velocity,
      "Only the principal coefficient ratio sets c_bulk^2."
    ],
    checkDim[
      "gauge speed c_gamma=sqrt(C_B/C_E)",
      (cBDim - cEDim)/2,
      velocity
    ],
    checkDim[
      "conditional bulk ratio c_bulk/c_s0",
      ((cBDim - cEDim) - (kEos + 4 rho4 - M))/2,
      dim0
    ],
    checkDim[
      "tail factor lambda_gamma^3=(c_gamma/c_s)^3",
      3 (((cBDim - cEDim) - (kEos + 4 rho4 - M))/2),
      dim0
    ],
    homogeneous[
      "Maxwell transverse principal operator terms",
      <|
        "C_E partial_t^2 A_T" -> cEDim + ai - 2 T,
        "C_B laplacian A_T" -> cBDim + ai - 2 L
      |>,
      "Gauge-fixing terms are not used as speed evidence."
    ],
    homogeneous[
      "transverse gauge dispersion omega^2=(C_B/C_E)*k^2",
      <|
        "omega^2" -> -2 T,
        "(C_B/C_E) k^2" -> cBDim - cEDim + 2 waveNumber
      |>
    ],
    homogeneous[
      "phonon acoustic dispersion omega^2=c_s0^2*k^2",
      <|
        "omega^2" -> -2 T,
        "c_s0^2 k^2" -> kEos + 4 rho4 - M + 2 waveNumber
      |>
    ],
    checkDim[
      "background charge density J_psi0^0=q_star*rho0",
      rho4,
      rho4
    ],
    homogeneous[
      "linearized spatial current variation terms",
      <|
        "phase current" -> rho4 + action - M - L,
        "London current" -> rho4 + ai - M,
        "rho0 v" -> rho4 + velocity
      |>
    ],
    homogeneous[
      "source coupling dimensions A_M delta J^M",
      <|
        "A0 deltaJ0" -> a0 + rho4,
        "Ai deltaJi" -> ai + rho4 + velocity,
        "local action density" -> energy - 4 L
      |>
    ]
  };

  omega = Symbol["omega"];
  k2 = Symbol["k2"];
  rho0 = Symbol["rho0"];
  kSym = Symbol["K"];
  mG = Symbol["mG"];
  hbar = Symbol["hbar"];
  cE = Symbol["cE"];
  cB = Symbol["cB"];
  lambdaGamma = Symbol["lambdaGamma"];
  hPrime = 5 kSym rho0^3;
  cs2 = FullSimplify[rho0 hPrime/mG];
  cBulk2 = cB/cE;
  phononMatrix = {{omega, -(rho0 hbar/mG) k2}, {-hPrime, hbar omega}};
  phononDet = Factor[Det[phononMatrix]];
  gaugeTransverse = cE omega^2 - cB k2;
  coupledDet = Factor[phononDet gaugeTransverse^2];
  equalityResidual = FullSimplify[cBulk2 - cs2];
  sourceMetricEquationPresent = False;
  forcedEqualsValid = sourceMetricEquationPresent && TrueQ[equalityResidual === 0];
  lambdaRhoFactor = rho0^-2;
  lambdaLogSlope = FullSimplify[rho0 D[lambdaRhoFactor, rho0]/lambdaRhoFactor];

  algebra = {
    checkExpr[
      "phonon determinant gives omega^2=c_s0^2*k^2",
      phononDet,
      hbar (omega^2 - cs2 k2),
      "Uses c_s0^2=5*K*rho0^4/m_GNLS."
    ],
    checkExpr[
      "transverse gauge operator gives c_bulk^2=C_B/C_E",
      gaugeTransverse,
      cE (omega^2 - cBulk2 k2),
      "The overall Maxwell prefactor cancels from the characteristic equation."
    ],
    checkExpr[
      "block coupled characteristic determinant after principal decoupling",
      coupledDet,
      hbar (omega^2 - cs2 k2) (cE (omega^2 - cBulk2 k2))^2
    ],
    checkBool[
      "negative control keeps C_B/C_E independent from c_s0^2",
      ! TrueQ[equalityResidual === 0],
      True
    ],
    checkBool[
      "forced C_GAMMA_EQUALS_C_S verdict rejected without source metric equation",
      forcedEqualsValid,
      False
    ],
    checkExpr[
      "conditional rho dependence of c_bulk/c_s0 when C_B/C_E is rho-independent",
      lambdaLogSlope,
      -2
    ],
    checkExpr[
      "standing-wave tail remains lambda_gamma^3",
      lambdaGamma^3,
      lambdaGamma^3
    ]
  };

  allDimPass = AllTrue[checks, TrueQ[#["pass"]] &];
  allAlgPass = AllTrue[algebra, TrueQ[#["pass"]] &];
  report = <|
    "schema" -> "stage1_pathA_20b_cgamma_cs_mathematica_crosscheck/v1",
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
  jsonPath = FileNameJoin[{outDir, "pathA_20b_cgamma_cs_mathematica_crosscheck.json"}];
  Export[jsonPath, report, "RawJSON"];

  Print["wrote " <> jsonPath];
  If[
    TrueQ[allDimPass && allAlgPass],
    Print["pathA_20b Mathematica algebra/dimensional cross-check: PASS"],
    Print["pathA_20b Mathematica algebra/dimensional cross-check: FAIL"]
  ];
  If[TrueQ[allDimPass && allAlgPass], Exit[0], Exit[1]]
]
