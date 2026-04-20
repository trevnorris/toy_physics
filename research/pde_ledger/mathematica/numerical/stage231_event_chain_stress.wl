(* Numerical stress harness for the Stage 231 event-chain compiler. *)

ClearAll["Global`*"];
$HistoryLength = 0;

rootDir = DirectoryName[$InputFileName];
configPath = FileNameJoin[{
   rootDir, "..", "..", "..", "scripts", "numerical",
   "stage231_event_chain_samples.json"
}];
config = Import[ExpandFileName[configPath], "RawJSON"];

If[config["schema"] =!= "moving_throat_numerical_stage231_v2",
  Print["Unexpected config schema."];
  Exit[1];
];

fmt[x_] := ToString @ NumberForm[N[x, 16], {Infinity, 12}, ExponentFunction -> (Null &)];
nearQ[lhs_, rhs_, tol_] := Abs[N[lhs - rhs, 30]] <= tol (1 + Abs[N[rhs, 30]]);

require[label_, condition_, detail_] := Module[{status},
  status = If[TrueQ[condition], "PASS", "FAIL"];
  Print["[", status, "] ", label, ": ", detail];
  If[!TrueQ[condition], Throw[$Failed]]
];

loweredBarrierPotential[r_, strength_, scale_] := 1/r - strength Exp[-r/scale];
loweredBarrierPrime[r_, strength_, scale_] := -1/r^2 + strength Exp[-r/scale]/scale;
loweredBarrierDoublePrime[r_, strength_, scale_] := 2/r^3 - strength Exp[-r/scale]/scale^2;
xiProfile[r_, xiFloor_, xiAmp_, xiScale_] := xiFloor + xiAmp Exp[-r/xiScale];

scanRoots[f_Function, lower_, upper_, samples_: 4000] := Module[
  {xs, vals, roots = {}, candidate, i},
  xs = N[Subdivide[lower, upper, samples], 50];
  vals = N[f /@ xs, 50];
  Do[
    candidate = Missing["NoRoot"];
    If[Abs[vals[[i]]] < 10^-20,
      candidate = xs[[i]],
      If[vals[[i]] vals[[i + 1]] < 0,
        candidate = Module[{left = xs[[i]], right = xs[[i + 1]], fLeft = vals[[i]], fRight = vals[[i + 1]], mid, fMid},
          Do[
            mid = (left + right)/2;
            fMid = f[mid];
            If[fLeft fMid <= 0,
              right = mid;
              fRight = fMid;,
              left = mid;
              fLeft = fMid;
            ],
            {80}
          ];
          (left + right)/2
        ]
      ]
    ];
    If[
      !MissingQ[candidate] &&
      NumericQ[candidate] &&
      lower < candidate < upper &&
      AllTrue[roots, Abs[candidate - #] > 10^-8 &],
      AppendTo[roots, candidate]
    ];
    ,
    {i, 1, Length[xs] - 1}
  ];
  Sort[roots]
];

coulombActionClosed[m_, hbar_, energy_, rContact_] :=
  Sqrt[2 m]/hbar (
    Pi/(2 Sqrt[energy]) -
    Sqrt[rContact (1 - energy rContact)] -
    ArcSin[Sqrt[energy rContact]]/Sqrt[energy]
  );

coulombActionQuad[m_, hbar_, energy_, rContact_] := Module[{turningPoint},
  turningPoint = 1/energy;
  Quiet[
    NIntegrate[
      Evaluate[N[Sqrt[2 m (1/r - energy)]/hbar, 40]],
      {r, rContact, (rContact + turningPoint)/2, turningPoint},
      WorkingPrecision -> 60,
      AccuracyGoal -> 20,
      PrecisionGoal -> 20,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
    ],
    NIntegrate::precw
  ]
];

loweredActionQuad[m_, hbar_, energy_, strength_, scale_, rMinus_, rPeak_, rPlus_] :=
  Quiet[
    NIntegrate[
      Evaluate[N[Sqrt[2 m (loweredBarrierPotential[r, strength, scale] - energy)]/hbar, 40]],
      {r, rMinus, rPeak, rPlus},
      WorkingPrecision -> 60,
      AccuracyGoal -> 20,
      PrecisionGoal -> 20,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
    ],
    NIntegrate::precw
  ];

classifyBarrierCase[case_Association] := Module[
  {
    mS, hbarEff, r0, rContact, strength, scale, margin,
    xiFloor, xiAmp, xiScale, epsilon = 10^-5,
    result, peakRoots, rPeak, vLaunch, vContact, vPeak, kPeak, energy,
    turningRoots, rMinus, rPlus, vCrit, vContactCoul, vSub, vCross,
    iNew, iCoulClosed, iCoulQuad, transmissionRatio, xiTurn, xiContact,
    xiLaunch, slopeExact, hStep, derivativeFD, lambdaTh, transportWindow,
    transportDelta, upperTurnings, upper2Turnings, lowerTurnings, lower2Turnings,
    transportFD, transportExact
  },
  mS = N[case["m_s"], 40];
  hbarEff = N[case["hbar_eff"], 40];
  r0 = N[case["r0"], 40];
  rContact = N[case["r_contact"], 40];
  strength = N[case["lowering_strength"], 40];
  scale = N[case["lowering_scale"], 40];
  margin = N[case["energy_margin_fraction"], 40];
  xiFloor = N[case["xi_floor"], 40];
  xiAmp = N[case["xi_amp"], 40];
  xiScale = N[case["xi_scale"], 40];

  result = <|"name" -> case["name"], "status" -> "failure", "reason" -> "unclassified"|>;

  peakRoots = scanRoots[
    Function[{radius}, loweredBarrierPrime[radius, strength, scale]],
    rContact + epsilon,
    r0 - epsilon
  ];
  AssociateTo[result, "peak_roots" -> peakRoots];
  If[Length[peakRoots] == 0,
    AssociateTo[result, "reason" -> "no_peak"];
    Return[result];
  ];
  If[Length[peakRoots] > 1,
    AssociateTo[result, "reason" -> "multiple_stationary_points"];
    Return[result];
  ];

  rPeak = First[peakRoots];
  vLaunch = loweredBarrierPotential[r0, strength, scale];
  vContact = loweredBarrierPotential[rContact, strength, scale];
  vPeak = loweredBarrierPotential[rPeak, strength, scale];
  kPeak = -loweredBarrierDoublePrime[rPeak, strength, scale];
  AssociateTo[result, <|
    "r_peak" -> rPeak,
    "V0" -> vLaunch,
    "V_contact" -> vContact,
    "V_peak" -> vPeak,
    "K_peak" -> kPeak
  |>];

  If[!(vLaunch < vContact < vPeak),
    AssociateTo[result, "reason" -> "contact_not_below_peak"];
    Return[result];
  ];
  If[!(0 < margin < 1),
    AssociateTo[result, "reason" -> "invalid_energy_margin"];
    Return[result];
  ];

  energy = vContact + margin (vPeak - vContact);
  AssociateTo[result, "E_sub" -> energy];
  If[energy rContact >= 1,
    AssociateTo[result, "reason" -> "coulomb_inadmissible"];
    Return[result];
  ];

  turningRoots = scanRoots[
    Function[{radius}, loweredBarrierPotential[radius, strength, scale] - energy],
    rContact + epsilon,
    r0 - epsilon
  ];
  AssociateTo[result, "turning_roots" -> turningRoots];
  If[Length[turningRoots] =!= 2,
    AssociateTo[result, "reason" -> "missing_turning_pair"];
    Return[result];
  ];

  {rMinus, rPlus} = turningRoots;
  If[!(rContact < rMinus < rPeak < rPlus < r0),
    AssociateTo[result, "reason" -> "turning_pair_ordering"];
    Return[result];
  ];

  vCrit = Sqrt[2 (vPeak - vLaunch)/mS];
  vContactCoul = Sqrt[2 (1/rContact - 1/r0)/mS];
  vSub = Sqrt[2 (energy - vLaunch)/mS];
  vCross = vCrit + N[case["v_cross_fraction"], 40] (vContactCoul - vCrit);
  iNew = loweredActionQuad[mS, hbarEff, energy, strength, scale, rMinus, rPeak, rPlus];
  iCoulClosed = coulombActionClosed[mS, hbarEff, energy, rContact];
  iCoulQuad = coulombActionQuad[mS, hbarEff, energy, rContact];
  transmissionRatio = Exp[-2 (iNew - iCoulClosed)];
  xiTurn = xiProfile[rPlus, xiFloor, xiAmp, xiScale];
  xiContact = xiProfile[rContact, xiFloor, xiAmp, xiScale];
  xiLaunch = xiProfile[r0, xiFloor, xiAmp, xiScale];
  slopeExact = Abs[loweredBarrierPrime[rPlus, strength, scale]];

  hStep = Min[
    (rPlus - rPeak)/80,
    (r0 - rPlus)/80,
    (rPlus - rMinus)/160,
    rPlus/500
  ];
  derivativeFD = Abs[
    (
      -loweredBarrierPotential[rPlus + 2 hStep, strength, scale] +
      8 loweredBarrierPotential[rPlus + hStep, strength, scale] -
      8 loweredBarrierPotential[rPlus - hStep, strength, scale] +
      loweredBarrierPotential[rPlus - 2 hStep, strength, scale]
    )/(12 hStep)
  ];
  lambdaTh = energy/slopeExact;

  transportWindow = Min[energy - vContact, vPeak - energy];
  transportDelta = 0.05 transportWindow;
  upperTurnings = scanRoots[
    Function[{radius}, loweredBarrierPotential[radius, strength, scale] - (energy + transportDelta)],
    rContact + epsilon,
    r0 - epsilon
  ];
  upper2Turnings = scanRoots[
    Function[{radius}, loweredBarrierPotential[radius, strength, scale] - (energy + 2 transportDelta)],
    rContact + epsilon,
    r0 - epsilon
  ];
  lowerTurnings = scanRoots[
    Function[{radius}, loweredBarrierPotential[radius, strength, scale] - (energy - transportDelta)],
    rContact + epsilon,
    r0 - epsilon
  ];
  lower2Turnings = scanRoots[
    Function[{radius}, loweredBarrierPotential[radius, strength, scale] - (energy - 2 transportDelta)],
    rContact + epsilon,
    r0 - epsilon
  ];
  If[
    Length[upper2Turnings] =!= 2 ||
    Length[upperTurnings] =!= 2 ||
    Length[lowerTurnings] =!= 2 ||
    Length[lower2Turnings] =!= 2,
    AssociateTo[result, "reason" -> "transport_window_breakdown"];
    Return[result];
  ];
  transportFD = (
    -upper2Turnings[[2]] +
    8 upperTurnings[[2]] -
    8 lowerTurnings[[2]] +
    lower2Turnings[[2]]
  )/(12 transportDelta);
  transportExact = 1/loweredBarrierPrime[rPlus, strength, scale];

  AssociateTo[result, <|
    "status" -> "success",
    "reason" -> "success",
    "r_minus" -> rMinus,
    "r_plus" -> rPlus,
    "v_sub" -> vSub,
    "v_crit" -> vCrit,
    "v_cross" -> vCross,
    "v_contact_coul" -> vContactCoul,
    "I_new" -> iNew,
    "I_coul_closed" -> iCoulClosed,
    "I_coul_quad" -> iCoulQuad,
    "transmission_ratio" -> transmissionRatio,
    "Xi_turn" -> xiTurn,
    "Xi_contact" -> xiContact,
    "Xi_launch" -> xiLaunch,
    "lambda_th" -> lambdaTh,
    "dV_turn_exact" -> slopeExact,
    "dV_turn_fd" -> derivativeFD,
    "transport_fd" -> transportFD,
    "transport_exact" -> transportExact
  |>];
  result
];

parabolicActionClosed[m_, hbar_, deltaE_, kPeak_] :=
  Pi deltaE Sqrt[m/kPeak]/hbar;

parabolicActionQuad[m_, hbar_, deltaE_, kPeak_] := Module[{yTurn},
  yTurn = Sqrt[2 deltaE/kPeak];
  Quiet[
    NIntegrate[
      Evaluate[N[Sqrt[2 m (deltaE - kPeak y^2/2)]/hbar, 40]],
      {y, -yTurn, 0, yTurn},
      WorkingPrecision -> 60,
      AccuracyGoal -> 20,
      PrecisionGoal -> 20,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
    ],
    NIntegrate::precw
  ]
];

eventChainBlock[] := Module[
  {
    actionTol = N[config["tolerances"]["action_rel_tol"], 30],
    derivativeTol = N[config["tolerances"]["derivative_rel_tol"], 30],
    transportTol = N[config["tolerances"]["transport_rel_tol"], 30]
  },
  Print[""];
  Print["=== Stage 231: dynamic event-chain stress ==="];

  Scan[
    Function[{case},
      Module[{result, expectedStatus, expectedReason},
        result = classifyBarrierCase[case];
        expectedStatus = case["expected_status"];
        expectedReason = Lookup[case, "expected_reason", Missing["NotAvailable"]];

        Print[""];
        Print[case["name"], ": expected=", expectedStatus];
        Print[
          "  raw barrier data: strength=", fmt[case["lowering_strength"]],
          ", scale=", fmt[case["lowering_scale"]],
          ", r_contact=", fmt[case["r_contact"]],
          ", r0=", fmt[case["r0"]]
        ];
        Print[
          "  classification: status=", result["status"],
          ", reason=", result["reason"]
        ];

        If[expectedStatus === "failure",
          require[
            case["name"] <> " expected failure status",
            result["status"] === "failure",
            "status=" <> result["status"]
          ];
          require[
            case["name"] <> " expected failure reason",
            result["reason"] === expectedReason,
            "reason=" <> result["reason"] <> ", expected=" <> expectedReason
          ],

          require[
            case["name"] <> " expected success status",
            result["status"] === "success",
            "status=" <> result["status"] <> ", reason=" <> result["reason"]
          ];

          Print[
            "  peak/turning data: r_peak=", fmt[result["r_peak"]],
            ", r_minus=", fmt[result["r_minus"]],
            ", r_plus=", fmt[result["r_plus"]]
          ];
          Print[
            "  barrier heights: V0=", fmt[result["V0"]],
            ", V_contact=", fmt[result["V_contact"]],
            ", V_peak=", fmt[result["V_peak"]],
            ", Esub=", fmt[result["E_sub"]]
          ];
          Print[
            "  actions: I_new=", fmt[result["I_new"]],
            ", I_coul_quad=", fmt[result["I_coul_quad"]],
            ", I_coul_closed=", fmt[result["I_coul_closed"]],
            ", T_ratio=", fmt[result["transmission_ratio"]]
          ];
          Print[
            "  diagnostics: Xi_turn=", fmt[result["Xi_turn"]],
            ", lambda_th=", fmt[result["lambda_th"]],
            ", |V'(r_+)|=", fmt[result["dV_turn_exact"]]
          ];

          require[
            case["name"] <> " dynamic corridor ordering",
            result["v_sub"] < result["v_crit"] < result["v_cross"] < result["v_contact_coul"],
            "v_sub=" <> fmt[result["v_sub"]] <>
            ", v_crit=" <> fmt[result["v_crit"]] <>
            ", v_cross=" <> fmt[result["v_cross"]] <>
            ", v_contact=" <> fmt[result["v_contact_coul"]]
          ];
          require[
            case["name"] <> " Coulomb action quadrature",
            nearQ[result["I_coul_quad"], result["I_coul_closed"], actionTol],
            "quad=" <> fmt[result["I_coul_quad"]] <>
            ", closed=" <> fmt[result["I_coul_closed"]]
          ];
          require[
            case["name"] <> " lowered branch beats Coulomb",
            result["I_new"] < result["I_coul_closed"] && result["transmission_ratio"] > 1,
            "I_new=" <> fmt[result["I_new"]] <>
            ", I_coul=" <> fmt[result["I_coul_closed"]] <>
            ", ratio=" <> fmt[result["transmission_ratio"]]
          ];
          require[
            case["name"] <> " Xi_turn derived from monotone profile",
            result["Xi_launch"] < result["Xi_turn"] < result["Xi_contact"],
            "Xi_launch=" <> fmt[result["Xi_launch"]] <>
            ", Xi_turn=" <> fmt[result["Xi_turn"]] <>
            ", Xi_contact=" <> fmt[result["Xi_contact"]]
          ];
          require[
            case["name"] <> " derivative at r_plus",
            nearQ[result["dV_turn_fd"], result["dV_turn_exact"], derivativeTol],
            "fd=" <> fmt[result["dV_turn_fd"]] <>
            ", exact=" <> fmt[result["dV_turn_exact"]]
          ];
          require[
            case["name"] <> " turning transport law",
            nearQ[result["transport_fd"], result["transport_exact"], transportTol],
            "fd=" <> fmt[result["transport_fd"]] <>
            ", exact=" <> fmt[result["transport_exact"]]
          ]
        ];
      ]
    ],
    config["barrier_cases"]
  ];
];

nearTopBlock[] := Module[{tol},
  Print[""];
  Print["=== Stage 231: near-top parabolic action stress ==="];
  tol = N[config["tolerances"]["near_top_rel_tol"], 30];

  Scan[
    Function[{case},
      Module[{mS, hbarEff, deltaE, kPeak, closed, quad, yTurn},
        mS = N[case["m_s"]];
        hbarEff = N[case["hbar_eff"]];
        deltaE = N[case["DeltaE"]];
        kPeak = N[case["Kpeak"]];
        closed = parabolicActionClosed[mS, hbarEff, deltaE, kPeak];
        quad = parabolicActionQuad[mS, hbarEff, deltaE, kPeak];
        yTurn = Sqrt[2 deltaE/kPeak];

        Print[""];
        Print[
          case["name"], ": m_s=", fmt[mS], ", hbar_eff=", fmt[hbarEff],
          ", DeltaE=", fmt[deltaE], ", Kpeak=", fmt[kPeak],
          ", y_turn=", fmt[yTurn]
        ];
        require[
          case["name"] <> " near-top quadrature",
          nearQ[quad, closed, tol],
          "quad=" <> fmt[quad] <> ", closed=" <> fmt[closed]
        ];
      ]
    ],
    config["near_top_cases"]
  ];
];

If[
  Catch[
    Print["Loaded config from ", ExpandFileName[configPath]];
    eventChainBlock[];
    nearTopBlock[];
    Print[""];
    Print["All Stage-231 numerical stress checks passed."];
    "ok"
  ] === "ok",
  Null,
  Print[""];
  Print["Stage-231 numerical stress harness failed."];
  Exit[1]
];
