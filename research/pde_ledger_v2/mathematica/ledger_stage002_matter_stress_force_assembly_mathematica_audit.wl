(* Ledger stage002 Mathematica audit: matter-stress force assembly.
   Print-only, no exports.  This route consumes the stage001 geometry
   primitives as inputs and assembles the force with native vector/matrix
   traction algebra plus Solve-based Gauss substitution. *)

ClearAll[
  heading, subheading, cleanZero, expectZero, squareLength,
  residualForUnits, dimensionalResiduals, surfaceCatalog, cfgValue,
  surfaceTractionLane, eosAndQuantumBlock, expectedAlgebra,
  actualAlgebra, algebraResiduals, allResiduals, mutationRows,
  mutationMustFail, omegaTwo, omegaThree, muThree, muFour,
  shellTwo, shellThree, mGnls, nInf, rhoBulk, qOne, qTwo,
  rSep, rBulk, vOneMag, kEos, rhoSym, vDot, hbarSym,
  mQuantum, x, amp, vDrain
];

heading[text_] := (
  Print[""];
  Print[StringRepeat["=", StringLength[text]]];
  Print[text];
  Print[StringRepeat["=", StringLength[text]]]
);

subheading[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);

cleanZero[expr_] := FullSimplify[expr] /. ConditionalExpression[0, _] -> 0;

expectZero[name_, residual_] := Module[{clean = cleanZero[residual]},
  If[TrueQ[clean === 0],
    Print["PASS  ", name],
    Print["FAIL  ", name, ": residual = ", ToString[InputForm[clean]]];
    Throw[name, "ledgerStage002Failure"]
  ]
];

squareLength[vec_] := FullSimplify[vec.vec];

residualForUnits[terms_Association, expected_] := FullSimplify[
  Total[Map[squareLength[# - expected] &, Values[terms]]]
];

dimensionalResiduals[overrides_: <||>] := Module[
  {
    l, t, mass, action, energy, force, velocity, rho3, rho4,
    numberRate, q3Vol, q4Vol, massDensity3, massDensity4,
    stress3, stress4, forceDensity3, forceDensity4,
    energyGradient, quantumPrefactor, want
  },

  l = {1, 0, 0};
  t = {0, 1, 0};
  mass = {0, 0, 1};
  action = 2 l - t + mass;
  energy = action - t;
  force = energy - l;
  velocity = l - t;
  rho3 = -3 l;
  rho4 = -4 l;
  numberRate = -t;
  q3Vol = numberRate - rho3;
  q4Vol = numberRate - rho4;
  massDensity3 = mass + rho3;
  massDensity4 = mass + rho4;
  stress3 = force - 2 l;
  stress4 = force - 3 l;
  forceDensity3 = force - 3 l;
  forceDensity4 = force - 4 l;
  energyGradient = energy - l;
  quantumPrefactor = 2 action - mass - 2 l;
  want[name_, default_] := Lookup[overrides, name, default];

  <|
    "reduced-3D momentum-balance terms" -> residualForUnits[
      <|
        "partial_t(m*N*v_i)" -> massDensity3 + velocity - t,
        "partial_j Pi_ij" -> stress3 - l,
        "V_conf body force" -> rho3 + energyGradient
      |>,
      want["reduced-3D momentum-balance terms", {-2, -2, 1}]
    ],
    "bulk-4D momentum-balance terms" -> residualForUnits[
      <|
        "partial_t(m*rho*v_i)" -> massDensity4 + velocity - t,
        "partial_J Pi_iJ" -> stress4 - l,
        "V_conf body force" -> rho4 + energyGradient
      |>,
      want["bulk-4D momentum-balance terms", {-3, -2, 1}]
    ],
    "reduced-3D Noether stress representatives" -> residualForUnits[
      <|
        "convective" -> massDensity3 + 2 velocity,
        "pressure" -> (energy - 4 l) + l,
        "quantum" -> quantumPrefactor + rho3
      |>,
      want["reduced-3D Noether stress representatives", {-1, -2, 1}]
    ],
    "Euler force-per-volume identity terms in reduced 3D" -> residualForUnits[
      <|
        "m*N*acceleration" -> massDensity3 + velocity - t,
        "N*grad h" -> rho3 + energyGradient,
        "N*grad Q" -> rho3 + energyGradient,
        "N*grad V_conf" -> rho3 + energyGradient,
        "N*q(E+vB)" -> rho3 + energyGradient
      |>,
      want["Euler force-per-volume identity terms in reduced 3D", {-2, -2, 1}]
    ],
    "reduced-3D surface traction integral gives force" -> residualForUnits[
      <|"stress3 + 2L" -> stress3 + 2 l|>,
      want["reduced-3D surface traction integral gives force", force]
    ],
    "bulk-4D surface traction integral gives force" -> residualForUnits[
      <|"stress4 + 3L" -> stress4 + 3 l|>,
      want["bulk-4D surface traction integral gives force", force]
    ],
    "reduced-3D Noether force coefficient m*N*Q1*Q2" -> residualForUnits[
      <|"massDensity3 + 2*q3Vol" -> massDensity3 + 2 q3Vol|>,
      want["reduced-3D Noether force coefficient m*N*Q1*Q2", force + 2 l]
    ],
    "bulk-4D Noether force coefficient m*rho4*Q1*Q2" -> residualForUnits[
      <|"massDensity4 + 2*q4Vol" -> massDensity4 + 2 q4Vol|>,
      want["bulk-4D Noether force coefficient m*rho4*Q1*Q2", force + 3 l]
    ],
    "V_conf body-force volume term dimension" -> residualForUnits[
      <|"forceDensity3 + 3L" -> forceDensity3 + 3 l|>,
      want["V_conf body-force volume term dimension", force]
    ]
  |>
];

surfaceCatalog[] := <|
  3 -> <|
    "omegaSymbol" -> omegaTwo,
    "omegaValue" -> 4 Pi,
    "momentSymbol" -> muThree,
    "momentValue" -> 1/3,
    "shellRadius" -> shellTwo
  |>,
  4 -> <|
    "omegaSymbol" -> omegaThree,
    "omegaValue" -> 2 Pi^2,
    "momentSymbol" -> muFour,
    "momentValue" -> 1/4,
    "shellRadius" -> shellThree
  |>
|>;

cfgValue[cfg_Association, key_, default_] := Lookup[cfg, key, default];

surfaceTractionLane[ambientDim_, density_, separation_, cfg_: <||>] := Module[
  {
    spec, omega, momentPressure, momentConvective, shell, area,
    e1, v1Vector, momentTensor, drainSign, forceSign,
    drainCoeff, convectiveVector, pressureVector, totalVector,
    gaussRule, forceVector, rules
  },

  spec = surfaceCatalog[][ambientDim];
  omega = spec["omegaSymbol"];
  momentPressure = spec["momentSymbol"];
  momentConvective = If[
    TrueQ[ambientDim === 3],
    cfgValue[cfg, "ConvectiveMoment3", momentPressure],
    momentPressure
  ];
  shell = spec["shellRadius"];
  area = omega shell^(ambientDim - 1);
  e1 = UnitVector[ambientDim, 1];
  v1Vector = vOneMag e1;
  momentTensor[mu_] := area mu IdentityMatrix[ambientDim];

  drainSign = cfgValue[cfg, "DrainSign", -1];
  forceSign = cfgValue[cfg, "ForceFluxSign", -1];
  drainCoeff = drainSign qTwo/(omega shell^(ambientDim - 1));

  convectiveVector = mGnls density drainCoeff (
    area v1Vector + momentTensor[momentConvective].v1Vector
  );
  pressureVector = -mGnls density drainCoeff (
    momentTensor[momentPressure].v1Vector
  );

  rules = {
    omega -> spec["omegaValue"],
    spec["momentSymbol"] -> spec["momentValue"]
  };
  totalVector = FullSimplify[(convectiveVector + pressureVector) /. rules];
  gaussRule = First[Solve[omega separation^(ambientDim - 1) vDrain + qOne == 0, vDrain]];
  forceVector = FullSimplify[
    (forceSign totalVector /. vOneMag -> (vDrain /. gaussRule)) /. rules
  ];

  <|
    "convective" -> FullSimplify[(convectiveVector[[1]] /. rules)],
    "pressure" -> FullSimplify[(pressureVector[[1]] /. rules)],
    "total" -> FullSimplify[totalVector[[1]]],
    "force" -> FullSimplify[forceVector[[1]]]
  |>
];

eosAndQuantumBlock[] := Module[
  {
    enthalpy, pressure, dhDrho, dPDrho, densityShift, pressureShift,
    densityProfile, quantumPotential, quantumStress
  },

  enthalpy = (5/4) kEos rhoSym^4;
  pressure = kEos rhoSym^5;
  dhDrho = D[enthalpy, rhoSym];
  dPDrho = D[pressure, rhoSym];
  densityShift = -mGnls vDot/dhDrho;
  pressureShift = FullSimplify[dPDrho densityShift];

  densityProfile = amp[x]^2;
  quantumPotential = -(hbarSym^2/(2 mQuantum)) D[amp[x], {x, 2}]/amp[x];
  quantumStress = (hbarSym^2/(4 mQuantum)) (
    D[densityProfile, x]^2/densityProfile - D[densityProfile, {x, 2}]
  );

  <|
    "EOS identity dP/drho = rho*dh/drho" -> FullSimplify[dPDrho - rhoSym dhDrho],
    "Bernoulli/EOS pressure cross term" -> pressureShift,
    "Madelung quantum stress divergence" -> FullSimplify[
      D[quantumStress, x] - densityProfile D[quantumPotential, x]
    ]
  |>
];

expectedAlgebra[overrides_: <||>] := Module[
  {surfaces, moment3, omega2, omega3, base},

  surfaces = surfaceCatalog[];
  moment3 = surfaces[3]["momentValue"];
  omega2 = surfaces[3]["omegaValue"];
  omega3 = surfaces[4]["omegaValue"];

  base = <|
    "EOS identity dP/drho = rho*dh/drho" -> 0,
    "Bernoulli/EOS pressure cross term" -> -mGnls rhoSym vDot,
    "Madelung quantum stress divergence" -> 0,
    "reduced-3D convective angular traction factor" ->
      -(1 + moment3) mGnls nInf qTwo vOneMag,
    "reduced-3D pressure angular traction factor" ->
      moment3 mGnls nInf qTwo vOneMag,
    "reduced-3D total flux from Noether stress" ->
      -mGnls nInf qTwo vOneMag,
    "reduced-3D force structure after Gauss substitution" ->
      -mGnls nInf qOne qTwo/(omega2 rSep^2),
    "bulk-4D total flux from Noether stress" ->
      -mGnls rhoBulk qTwo vOneMag,
    "bulk-4D force structure after Gauss substitution" ->
      -mGnls rhoBulk qOne qTwo/(omega3 rBulk^3)
  |>;
  Join[base, overrides]
];

actualAlgebra[cfg_: <||>] := Module[
  {local, reduced, bulk},

  local = eosAndQuantumBlock[];
  reduced = surfaceTractionLane[3, nInf, rSep, cfg];
  bulk = surfaceTractionLane[4, rhoBulk, rBulk, cfg];

  Join[
    local,
    <|
      "reduced-3D convective angular traction factor" -> reduced["convective"],
      "reduced-3D pressure angular traction factor" -> reduced["pressure"],
      "reduced-3D total flux from Noether stress" -> reduced["total"],
      "reduced-3D force structure after Gauss substitution" -> reduced["force"],
      "bulk-4D total flux from Noether stress" -> bulk["total"],
      "bulk-4D force structure after Gauss substitution" -> bulk["force"]
    |>
  ]
];

algebraResiduals[cfg_: <||>, overrides_: <||>] := Module[
  {actual, expected},
  actual = actualAlgebra[cfg];
  expected = expectedAlgebra[overrides];
  Association @ KeyValueMap[#1 -> cleanZero[actual[#1] - #2] &, expected]
];

allResiduals[cfg_: <||>, dimOverrides_: <||>, algebraOverrides_: <||>] := Join[
  dimensionalResiduals[dimOverrides],
  algebraResiduals[cfg, algebraOverrides]
];

mutationRows[] := Module[
  {reducedForce, bulkForce, dimForce, convective, total},

  reducedForce = "reduced-3D force structure after Gauss substitution";
  bulkForce = "bulk-4D force structure after Gauss substitution";
  dimForce = "reduced-3D surface traction integral gives force";
  convective = "reduced-3D convective angular traction factor";
  total = "reduced-3D total flux from Noether stress";

  {
    <|
      "label" -> "flip sign of reduced-3D force law",
      "config" -> <||>,
      "required" -> {reducedForce},
      "algebraOverrides" -> <|
        reducedForce -> mGnls nInf qOne qTwo/(4 Pi rSep^2)
      |>
    |>,
    <|
      "label" -> "change 4*pi -> 2*pi in reduced-3D force law",
      "config" -> <||>,
      "required" -> {reducedForce},
      "algebraOverrides" -> <|
        reducedForce -> -mGnls nInf qOne qTwo/(2 Pi rSep^2)
      |>
    |>,
    <|
      "label" -> "change 2*pi^2 -> 4*pi^2 in bulk-4D force law",
      "config" -> <||>,
      "required" -> {bulkForce},
      "algebraOverrides" -> <|
        bulkForce -> -mGnls rhoBulk qOne qTwo/(4 Pi^2 rBulk^3)
      |>
    |>,
    <|
      "label" -> "corrupt reduced-3D surface-traction dimensional exponent",
      "config" -> <||>,
      "required" -> {dimForce},
      "dimOverrides" -> <|dimForce -> {2, -2, 1}|>
    |>,
    <|
      "label" -> "corrupt consumed convective second moment 1/3 -> 1/2",
      "config" -> <|"ConvectiveMoment3" -> 1/2|>,
      "required" -> {convective, total, reducedForce}
    |>,
    <|
      "label" -> "flip F = -int Pi.n force-definition sign",
      "config" -> <|"ForceFluxSign" -> 1|>,
      "required" -> {reducedForce}
    |>,
    <|
      "label" -> "flip v2 drain-inflow Gauss sign",
      "config" -> <|"DrainSign" -> 1|>,
      "required" -> {reducedForce}
    |>
  }
];

mutationMustFail[row_] := Module[
  {residuals, outcomes, failedAll, missing},

  residuals = allResiduals[
    Lookup[row, "config", <||>],
    Lookup[row, "dimOverrides", <||>],
    Lookup[row, "algebraOverrides", <||>]
  ];

  outcomes = Association @ Table[
    check -> Catch[
      expectZero[check, residuals[check]];
      "SURVIVED",
      "ledgerStage002Failure",
      Function[{name, tag}, "FAILED"]
    ],
    {check, row["required"]}
  ];

  failedAll = AllTrue[Values[outcomes], TrueQ[# === "FAILED"] &];
  If[failedAll,
    Print[
      "PASS  mutation probe: ", row["label"],
      " produced required FAIL (", StringRiffle[row["required"], ", "], ")"
    ];
    True,
    missing = Keys[Select[outcomes, # =!= "FAILED" &]];
    Print[
      "FAIL  mutation probe: ", row["label"],
      " survived for ", StringRiffle[missing, ", "]
    ];
    False
  ]
];

Module[
  {dimKeys, algebraKeys, residuals, baselineOutcome, probeOutcome},

  heading["ledger_stage002_matter_stress_force_assembly Mathematica audit"];
  residuals = allResiduals[];
  dimKeys = Keys[dimensionalResiduals[]];
  algebraKeys = Keys[expectedAlgebra[]];

  baselineOutcome = Catch[
    subheading["Baseline dimensional residuals"];
    Scan[Function[key, expectZero[key, residuals[key]]], dimKeys];

    subheading["Baseline algebraic residuals"];
    Scan[Function[key, expectZero[key, residuals[key]]], algebraKeys];

    Print[""];
    Print["Verdict labels:"];
    Print["  leading matter-stress sign: FORCE_ATTRACTIVE_DERIVED"];
    Print["  full sign: SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE"];
    Print["  acceptance: PASS_WITH_NAMED_RESIDUALS"];
    True,
    "ledgerStage002Failure",
    Function[{name, tag}, False]
  ];

  subheading["Able-to-fail mutation probe"];
  probeOutcome = AllTrue[mutationRows[], mutationMustFail];

  If[TrueQ[baselineOutcome && probeOutcome],
    Print[""];
    Print["OVERALL PASS: Mathematica derived the stage002 matter-stress force assembly exactly"];
    Exit[0],
    Print[""];
    Print["OVERALL FAIL: Mathematica stage002 audit did not close"];
    Exit[1]
  ]
]
