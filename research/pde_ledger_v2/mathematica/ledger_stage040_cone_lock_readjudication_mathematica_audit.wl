(* Ledger stage040 Mathematica audit: cone-lock re-adjudication.

   Standalone, print-only, assert-zero, and cross-engine-file-I/O-free.
   The production interpretation is CALIBRATED / UNCOMMITTED: neither cone
   lock is an earned equality.

   This engine uses Mathematica-native semialgebraic primitives.  Resolve
   proves Exists/ForAll statements, FindInstance supplies non-entailment,
   RegionDimension[ImplicitRegion[...]] computes codimension, and
   FullSimplify/Solve derive the field overlay.  It does not read SymPy.

   Tooth-local ablation uses LEDGER_STAGE040_MUTATION.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE040_MUTATION";
activeMutation = Quiet@Check[Environment[mutationEnvironment], ""];
If[! StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];

toothOrder = {
  "ROUTE_A_GRADE",
  "ROUTE_B_STATUS",
  "FREEDOM_STATUS",
  "SECTOR_SAT",
  "LOCK_SAT",
  "PROVENANCE_LOCK_A",
  "LOCK_A_WITNESS_VALUE",
  "PROVENANCE_LOCK_B",
  "LOCK_B_WITNESS_VALUE",
  "DELTA_R",
  "PRODUCTION_VERDICT",
  "FIELD_OVERLAY_DET_ON_CONE",
  "FIELD_OVERLAY_POLES",
  "OPEN_110_CARRY",
  "EARNED_VS_CALIBRATED_PARTITION",
  "CONDITIONALITY_FREEDOM",
  "CTRL_WELL_POSED",
  "CTRL_ABSENT",
  "CTRL_PARTIAL_INVENTORY",
  "CTRL_FORCED_LOCK",
  "CTRL_A_ONLY_PARTIAL",
  "CTRL_B_ONLY_PARTIAL",
  "CTRL_OVER_CONSTRAINED",
  "CTRL_FREEDOM_TIE",
  "VERDICT_REDERIVATION",
  "DIMENSION_HOMOGENEITY",
  "DUAL_ENGINE_TERMS",
  "SOURCE_TO_STAGE_MANIFEST"
};

ablationDescriptions = <|
  "ROUTE_A_GRADE" ->
    "mark R1 present in the private production inventory",
  "ROUTE_B_STATUS" ->
    "inject the rejected thin-plate over-import relation",
  "FREEDOM_STATUS" ->
    "inject a freedom-tie relation into the private freedom inventory",
  "SECTOR_SAT" ->
    "inject the positive over-stability contradiction",
  "LOCK_SAT" ->
    "inject the freedom-tie contradiction into E plus both locks",
  "PROVENANCE_LOCK_A" ->
    "inject the target-blind force-A bridge",
  "LOCK_A_WITNESS_VALUE" ->
    "change mu_R in the exact A non-entailment witness",
  "PROVENANCE_LOCK_B" ->
    "inject the target-blind force-B bridge",
  "LOCK_B_WITNESS_VALUE" ->
    "change c_E and its dependent base witness values",
  "DELTA_R" ->
    "drop lock B from the private after region",
  "PRODUCTION_VERDICT" ->
    "feed decide the independently computed force-B case",
  "FIELD_OVERLAY_DET_ON_CONE" ->
    "use q=B_eff/rho_br; determinant residual survives but the shear guard fails",
  "FIELD_OVERLAY_POLES" ->
    "drop the C_hu coupling from the private pole polynomial",
  "OPEN_110_CARRY" ->
    "replace the on-cone residual by a C_hu-independent residual",
  "EARNED_VS_CALIBRATED_PARTITION" ->
    "restore the source bug by inserting L_A into earned equalities",
  "CONDITIONALITY_FREEDOM" ->
    "suppress the tied branch before recomputing its verdict",
  "CTRL_WELL_POSED" ->
    "remove R5 from the synthetic well-posed inventory",
  "CTRL_ABSENT" ->
    "complete R1..R5 in the absent inventory",
  "CTRL_PARTIAL_INVENTORY" ->
    "complete R3..R5 in the partial inventory",
  "CTRL_FORCED_LOCK" ->
    "drop the synthetic bridge",
  "CTRL_A_ONLY_PARTIAL" ->
    "drop the force-A bridge",
  "CTRL_B_ONLY_PARTIAL" ->
    "drop the force-B bridge",
  "CTRL_OVER_CONSTRAINED" ->
    "drop the over-stability contradiction",
  "CTRL_FREEDOM_TIE" ->
    "drop the freedom-tie contradiction",
  "VERDICT_REDERIVATION" ->
    "replace only the forced-lock table row by a force-A computed case",
  "DIMENSION_HOMOGENEITY" ->
    "change [c_E] from L T^-1 to L^2 T^-1",
  "DUAL_ENGINE_TERMS" ->
    "drop the locally computed lock-B witness lane",
  "SOURCE_TO_STAGE_MANIFEST" ->
    "drop one source predicate, mis-scope another, and drop one live registry tooth"
|>;

raise[message_] := Throw[message, "ledgerStage040Failure"];

expectBool[name_, condition_, evidence_: None] := If[
  TrueQ[condition],
  passCount += 1;
  Print["PASS  ", name],
  failCount += 1;
  Print["FIRST_FAILURE=", name];
  If[activeMutation === name, Print["FIRED_AT_OWN_ASSERT=", name]];
  Print["FAIL  ", name, ": residual = 1"];
  If[evidence =!= None, Print["      evidence = ", evidence]];
  raise[name]
];

section[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);


(* ---------------------------------------------------------------------- *)
(* Neutral source facts and the semialgebraic cone-lock system.            *)
(* ---------------------------------------------------------------------- *)

baseVars = {
  rho, rhoBr, rhoB0, chiC, muR, bigK, massM,
  Mh, cE, Chu, Kh, Beff, sigma
};
basePositive = {
  rho, rhoBr, chiC, muR, bigK, massM,
  Mh, cE, Kh, Beff, sigma
};
baseEquations = {
  Beff*chiC - rhoB0^2,
  Kh - Mh*cE^2,
  Beff*Kh - Chu^2 - sigma
};
lockA = massM*muR - 5*bigK*rho^4*rhoBr;
lockB = cE^2*rhoBr - muR;
lockExpression = <|"A" -> lockA, "B" -> lockB|>;

rKinds = <|
  "R1" -> "nonlinear_action_with_gnls_and_u",
  "R2" -> "in_plane_shear_profile_fu",
  "R3" -> "common_normalization_rho_mu",
  "R4" -> "reduction_mu_over_rho_equals_cs",
  "R5" -> "no_residual_geometric_factor"
|>;

baselineKinds[] := {
  "candidate_bridge",
  "h_goldstone_profile_imported",
  "postulated_mu_rho",
  "surface_shear_postulated",
  "missing_closed_fu",
  "route_b_distinct_dof",
  "over_import_guard",
  "freedom_free_parameter:C_hu",
  "freedom_free_parameter:rho_br"
};

caseKinds[case_String] := Switch[
  case,
  "production",
    baselineKinds[],
  "well_posed",
    {
      "candidate_bridge",
      Sequence @@ Values[rKinds],
      "route_b_distinct_dof",
      "over_import_guard",
      "freedom_free_parameter:C_hu",
      "freedom_free_parameter:rho_br"
    },
  "absent",
    {
      "postulated_mu_rho",
      "route_b_distinct_dof",
      "over_import_guard",
      "freedom_free_parameter:C_hu",
      "freedom_free_parameter:rho_br"
    },
  "partial_inventory",
    {
      "candidate_bridge",
      rKinds["R1"],
      rKinds["R2"],
      "postulated_mu_rho",
      "route_b_distinct_dof",
      "over_import_guard",
      "freedom_free_parameter:C_hu",
      "freedom_free_parameter:rho_br"
    },
  "forced_lock",
    Append[baselineKinds[], "forced_lock_fake_relation"],
  "a_only_partial",
    Append[baselineKinds[], "force_A_fake_relation"],
  "b_only_partial",
    Append[baselineKinds[], "force_B_fake_relation"],
  "over_constrained",
    Append[baselineKinds[], "over_stability_relation"],
  "freedom_tie",
    Append[baselineKinds[], "freedom_tie_relation"],
  "a_inconclusive",
    Append[baselineKinds[], "inconclusive_A_witness"],
  _,
    raise["unknown case: " <> case]
];

routeA[kinds_List] := Module[{present, missing, grade},
  present = Association@KeyValueMap[
    #1 -> MemberQ[kinds, #2] &,
    rKinds
  ];
  If[! TrueQ[present["R4"]], present["R5"] = False];
  missing = Keys@Select[present, Not@*TrueQ];
  grade = Which[
    missing === {},
      "ROUTE_A_WELL_POSED",
    MemberQ[kinds, "candidate_bridge"],
      "ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT",
    True,
      "ROUTE_A_ABSENT"
  ];
  <|"grade" -> grade, "missing" -> missing|>
];

routeB[kinds_List] := If[
  MemberQ[kinds, "route_b_distinct_dof"] &&
    ! MemberQ[kinds, "thin_plate_over_import_relation"],
  "ROUTE_B_CLOSED_CHECKED_NEGATIVE",
  "ROUTE_B_GUARD_FAIL"
];

freedomRecord[kinds_List] := Module[{free},
  If[
    MemberQ[kinds, "freedom_tie_relation"],
    Return[<|"status" -> "FREEDOM_TIED", "freeParams" -> {}|>]
  ];
  free = Sort@StringDelete[
    Select[
      kinds,
      StringStartsQ[#, "freedom_free_parameter:"] &
    ],
    "freedom_free_parameter:"
  ];
  <|
    "status" ->
      If[
        free === {},
        "FREEDOM_TIED",
        "FREEDOM_UNCONSTRAINED{" <>
          StringRiffle[free, ","] <> "}"
      ],
    "freeParams" -> free
  |>
];

sourceRelations[kinds_List] := Module[
  {equations = {}, variables = {}, positive = {}},
  If[
    MemberQ[kinds, "forced_lock_fake_relation"],
    equations = Join[
      equations,
      {
        muR - rhoBr*tauForced,
        cE^2 - tauForced,
        massM*tauForced - 5*bigK*rho^4
      }
    ];
    AppendTo[variables, tauForced];
    AppendTo[positive, tauForced]
  ];
  If[
    MemberQ[kinds, "force_A_fake_relation"],
    equations = Join[
      equations,
      {
        muR - rhoBr*tauA,
        massM*tauA - 5*bigK*rho^4
      }
    ];
    AppendTo[variables, tauA];
    AppendTo[positive, tauA]
  ];
  If[
    MemberQ[kinds, "force_B_fake_relation"],
    equations = Join[
      equations,
      {muR - rhoBr*tauB, cE^2 - tauB}
    ];
    AppendTo[variables, tauB];
    AppendTo[positive, tauB]
  ];
  If[
    MemberQ[kinds, "over_stability_relation"],
    AppendTo[equations, Chu^2 - Beff*Kh - etaOver];
    AppendTo[variables, etaOver];
    AppendTo[positive, etaOver]
  ];
  If[
    MemberQ[kinds, "freedom_tie_relation"],
    equations = Join[
      equations,
      {
        Chu^2 - qhSq,
        qhSq*rhoBr - 2*Beff*Mh*muR
      }
    ];
    AppendTo[variables, qhSq];
    AppendTo[positive, qhSq]
  ];
  <|
    "equations" -> equations,
    "variables" -> variables,
    "positive" -> positive
  |>
];

variablesFor[kinds_List] :=
  Join[baseVars, sourceRelations[kinds]["variables"]];

positiveFor[kinds_List] :=
  Join[basePositive, sourceRelations[kinds]["positive"]];

equationsFor[kinds_List, locks_List] := Join[
  baseEquations,
  sourceRelations[kinds]["equations"],
  Lookup[lockExpression, locks]
];

domainFor[kinds_List] :=
  And @@ Thread[positiveFor[kinds] > 0] &&
  rhoB0 != 0;

formulaFor[kinds_List, locks_List] :=
  And @@ Thread[equationsFor[kinds, locks] == 0] &&
  domainFor[kinds];

satRecord[kinds_List, locks_List] :=
  satRecord[kinds, locks] = Module[{variables, result},
    variables = variablesFor[kinds];
    result = TimeConstrained[
      Resolve[
        Exists[
          Evaluate[variables],
          formulaFor[kinds, locks]
        ],
        Reals
      ],
      120,
      "UNKNOWN"
    ];
    <|
      "status" -> Which[
        result === True, "SAT",
        result === False, "UNSAT",
        True, "SAT_UNCERTIFIED"
      ],
      "resolveResult" -> result
    |>
  ];

baseWitness[] := {
  rho -> 1,
  rhoBr -> 1,
  rhoB0 -> 2,
  chiC -> 4,
  muR -> 1,
  bigK -> 1,
  massM -> 5,
  Mh -> 1,
  cE -> 1,
  Chu -> 1/2,
  Kh -> 1,
  Beff -> 1,
  sigma -> 3/4
};

unlockedWitness[] := Join[
  DeleteCases[
    baseWitness[],
    (muR | cE | Kh | sigma) -> _
  ],
  {muR -> 2, cE -> 3, Kh -> 9, sigma -> 35/4}
];

candidateWitness[kinds_List, lockName_String] := Module[{witness},
  If[
    MemberQ[kinds, "forced_lock_fake_relation"] ||
    (lockName === "A" &&
      (MemberQ[kinds, "force_A_fake_relation"] ||
       MemberQ[kinds, "inconclusive_A_witness"])) ||
    (lockName === "B" &&
      MemberQ[kinds, "force_B_fake_relation"]) ||
    MemberQ[kinds, "over_stability_relation"] ||
    MemberQ[kinds, "freedom_tie_relation"],
    Return[Missing["NotAvailable"]]
  ];
  witness = Which[
    MemberQ[kinds, "force_A_fake_relation"],
      Join[
        DeleteCases[baseWitness[], (cE | Kh | sigma) -> _],
        {cE -> 2, Kh -> 4, sigma -> 15/4, tauA -> 1}
      ],
    MemberQ[kinds, "force_B_fake_relation"],
      Join[
        DeleteCases[
          baseWitness[],
          (muR | cE | Kh | sigma) -> _
        ],
        {
          muR -> 4,
          cE -> 2,
          Kh -> 4,
          sigma -> 15/4,
          tauB -> 4
        }
      ],
    True,
      unlockedWitness[]
  ];
  witness
];

witnessValidQ[kinds_List, locks_List, witness_List] := TrueQ[
  FullSimplify[formulaFor[kinds, locks] /. witness]
];

entailRecord[kinds_List, lockName_String] :=
  entailRecord[kinds, lockName] = Module[
    {
      variables,
      baseFormula,
      lock,
      universal,
      instance,
      witness,
      value,
      valid
    },
    variables = variablesFor[kinds];
    baseFormula = formulaFor[kinds, {}];
    lock = lockExpression[lockName];
    universal = TimeConstrained[
      Resolve[
        ForAll[
          Evaluate[variables],
          Implies[baseFormula, lock == 0]
        ],
        Reals
      ],
      120,
      "UNKNOWN"
    ];
    If[
      universal === True,
      Return[
        <|
          "status" -> "ENTAILED",
          "universal" -> universal
        |>
      ]
    ];
    instance = If[
      MemberQ[kinds, "inconclusive_A_witness"] &&
        lockName === "A",
      $Failed,
      TimeConstrained[
        FindInstance[
          baseFormula && lock != 0,
          Evaluate[variables],
          Reals,
          1
        ],
        120,
        $Failed
      ]
    ];
    witness = candidateWitness[kinds, lockName];
    If[
      MatchQ[witness, _Missing],
      Return[
        <|
          "status" -> "INCONCLUSIVE",
          "universal" -> universal,
          "instance" -> instance
        |>
      ]
    ];
    value = FullSimplify[lock /. witness];
    valid = witnessValidQ[kinds, {}, witness];
    If[
      universal === False &&
        ListQ[instance] &&
        Length[instance] > 0 &&
        valid &&
        TrueQ[value != 0],
      <|
        "status" -> "WITNESSED",
        "universal" -> universal,
        "instance" -> instance,
        "witness" -> witness,
        "value" -> value
      |>,
      <|
        "status" -> "INCONCLUSIVE",
        "universal" -> universal,
        "instance" -> instance,
        "witness" -> witness,
        "value" -> value,
        "valid" -> valid
      |>
    ]
  ];

regionDimensionFor[kinds_List, locks_List] :=
  regionDimensionFor[kinds, locks] = Module[
    {variables, region},
    variables = variablesFor[kinds];
    region = ImplicitRegion[
      formulaFor[kinds, locks],
      Evaluate[variables]
    ];
    TimeConstrained[
      RegionDimension[region],
      120,
      "DIMENSION_UNCERTIFIED"
    ]
  ];

dimensionRecord[kinds_List, afterLocks_List : {"A", "B"}] :=
  dimensionRecord[kinds, afterLocks] = Module[
    {beforeSat, afterSat, before, after},
    beforeSat = satRecord[kinds, {}]["status"];
    afterSat = satRecord[kinds, afterLocks]["status"];
    If[
      beforeSat =!= "SAT" || afterSat =!= "SAT",
      Return[<|"status" -> "NOT_RUN"|>]
    ];
    before = regionDimensionFor[kinds, {}];
    after = regionDimensionFor[kinds, afterLocks];
    If[
      IntegerQ[before] && IntegerQ[after],
      <|
        "status" -> "CERTIFIED",
        "dimBefore" -> before,
        "dimAfter" -> after,
        "deltaR" -> before - after
      |>,
      <|
        "status" -> "DIMENSION_UNCERTIFIED",
        "dimBefore" -> before,
        "dimAfter" -> after
      |>
    ]
  ];

decideOutput[
  routeAGrade_String,
  sectorStatus_String,
  lockStatus_String,
  provA_String,
  provB_String,
  dimensionStatus_String,
  deltaR_
] := Module[{riders},
  riders = {};
  Which[
    sectorStatus === "UNSAT",
      <|"verdict" -> "NO_GO(sector-ledger)", "riders" -> riders|>,
    routeAGrade === "ROUTE_A_WELL_POSED",
      <|
        "verdict" -> "HALT_ROUTE_A_WELL_POSED",
        "riders" -> riders
      |>,
    lockStatus === "UNSAT",
      <|"verdict" -> "NO_GO(cone-lock)", "riders" -> riders|>,
    True,
      If[
        lockStatus === "SAT_UNCERTIFIED",
        AppendTo[riders, "SAT_UNCERTIFIED"]
      ];
      If[
        provA === "INCONCLUSIVE",
        AppendTo[riders, "ENTAILMENT_INCONCLUSIVE(L_A)"]
      ];
      If[
        provB === "INCONCLUSIVE",
        AppendTo[riders, "ENTAILMENT_INCONCLUSIVE(L_B)"]
      ];
      Which[
        riders =!= {},
          <|
            "verdict" -> "CONE_LOCK_PROVENANCE_INCONCLUSIVE",
            "riders" -> riders
          |>,
        dimensionStatus =!= "CERTIFIED" ||
          lockStatus =!= "SAT",
          <|
            "verdict" -> "CONE_LOCK_PROVENANCE_INCONCLUSIVE",
            "riders" -> riders
          |>,
        provA === "ENTAILED" &&
          provB === "ENTAILED" &&
          deltaR === 0,
          <|"verdict" -> "CONE_LOCK_DERIVED", "riders" -> riders|>,
        provA === "ENTAILED" &&
          provB === "WITNESSED" &&
          deltaR === 1,
          <|
            "verdict" ->
              "CONE_LOCK_PARTIAL(derived=A, calibrated=B)",
            "riders" -> riders
          |>,
        provB === "ENTAILED" &&
          provA === "WITNESSED" &&
          deltaR === 1,
          <|
            "verdict" ->
              "CONE_LOCK_PARTIAL(derived=B, calibrated=A)",
            "riders" -> riders
          |>,
        provA === "WITNESSED" &&
          provB === "WITNESSED" &&
          deltaR === 2,
          <|
            "verdict" -> "CONE_LOCK_CALIBRATED",
            "riders" -> riders
          |>,
        provA === "WITNESSED" &&
          provB === "WITNESSED" &&
          deltaR === 1,
          <|
            "verdict" ->
              "CONE_LOCK_SHARED_CALIBRATION(delta_r=1, derived=none)",
            "riders" -> riders
          |>,
        True,
          <|
            "verdict" -> "CONE_LOCK_PROVENANCE_INCONCLUSIVE",
            "riders" -> riders
          |>
      ]
  ]
];

runCase[kinds_List] :=
  runCase[kinds] = Module[
    {
      routeAResult,
      routeBResult,
      freedomResult,
      sector,
      locked,
      provA,
      provB,
      dimension,
      decision
    },
    routeAResult = routeA[kinds];
    routeBResult = routeB[kinds];
    freedomResult = freedomRecord[kinds];
    If[
      routeAResult["grade"] === "ROUTE_A_WELL_POSED",
      sector = "NOT_RUN";
      locked = "NOT_RUN";
      provA = "NOT_RUN";
      provB = "NOT_RUN";
      dimension = <|"status" -> "NOT_RUN"|>,
      sector = satRecord[kinds, {}]["status"];
      locked = satRecord[kinds, {"A", "B"}]["status"];
      If[
        sector === "SAT" && locked === "SAT",
        provA = entailRecord[kinds, "A"]["status"];
        provB = entailRecord[kinds, "B"]["status"];
        dimension = dimensionRecord[kinds],
        provA = "NOT_RUN";
        provB = "NOT_RUN";
        dimension = <|"status" -> "NOT_RUN"|>
      ]
    ];
    decision = decideOutput[
      routeAResult["grade"],
      sector,
      locked,
      provA,
      provB,
      dimension["status"],
      Lookup[dimension, "deltaR", Null]
    ];
    <|
      "routeA" -> routeAResult,
      "routeB" -> routeBResult,
      "freedom" -> freedomResult,
      "sectorSat" -> sector,
      "lockSat" -> locked,
      "provA" -> provA,
      "provB" -> provB,
      "dimension" -> dimension,
      "verdict" -> decision["verdict"],
      "riders" -> decision["riders"]
    |>
  ];


(* ---------------------------------------------------------------------- *)
(* FullSimplify/Solve field overlay; dimension and partition helpers.      *)
(* ---------------------------------------------------------------------- *)

fieldOverlay[
  competingCone_: False,
  dropCouplingForPoles_: False,
  openConditionMutation_: False
] := Module[
  {
    qGamma,
    determinantQ,
    detOnCone,
    shearResidual,
    coupling,
    polePolynomial,
    roots,
    deltas,
    truePolynomial,
    coefficients,
    expectedDiscriminant,
    expectedRoots,
    expectedDeltas,
    coincidenceResidual,
    coincidenceCondition
  },
  qGamma = If[TrueQ[competingCone], Beff/rhoBr, muR/rhoBr];
  determinantQ =
    ((rhoBr*qCone - Beff)*(Mh*qCone - Kh) - Chu^2)*kWave^4;
  detOnCone = FullSimplify[
    determinantQ /. {
      qCone -> qGamma,
      Kh -> Mh*muR/rhoBr
    }
  ];
  shearResidual = FullSimplify[rhoBr*qGamma - muR];

  coupling = If[TrueQ[dropCouplingForPoles], 0, Chu^2];
  polePolynomial =
    (rhoBr*qCone - Beff)*(Mh*qCone - Kh) - coupling;
  roots = qCone /. Solve[polePolynomial == 0, qCone];
  deltas = FullSimplify[
    (# - muR/rhoBr) /. muR -> Kh*rhoBr/Mh
  ] & /@ roots;

  truePolynomial =
    (rhoBr*qCone - Beff)*(Mh*qCone - Kh) - Chu^2;
  coefficients = CoefficientList[truePolynomial, qCone];
  expectedDiscriminant = FullSimplify[
    coefficients[[2]]^2 -
      4*coefficients[[3]]*coefficients[[1]]
  ];
  expectedRoots = FullSimplify /@ {
    (-coefficients[[2]] - Sqrt[expectedDiscriminant])/
      (2*coefficients[[3]]),
    (-coefficients[[2]] + Sqrt[expectedDiscriminant])/
      (2*coefficients[[3]])
  };
  expectedDeltas = FullSimplify[
    (# - muR/rhoBr) /. muR -> Kh*rhoBr/Mh
  ] & /@ expectedRoots;

  coincidenceResidual = If[
    TrueQ[openConditionMutation],
    -kWave^4,
    detOnCone
  ];
  coincidenceCondition = FullSimplify[
    Reduce[
      (coincidenceResidual /. kWave -> 1) == 0,
      Chu,
      Reals
    ]
  ];
  <|
    "detOnCone" -> detOnCone,
    "shearResidual" -> shearResidual,
    "deltas" -> deltas,
    "expectedDeltas" -> expectedDeltas,
    "coincidenceCondition" -> coincidenceCondition
  |>
];

dimAdd[left_List, right_List] := MapThread[Plus, {left, right}];
dimScale[power_Integer, value_List] := power*value;

dimensionState[mutateCE_: False] := Module[
  {
    dims,
    speedSquared,
    lockADims,
    lockBDims,
    stabilityDims,
    homogeneous
  },
  dims = <|
    "rho" -> {1, -3, 0},
    "rho_br" -> {1, -3, 0},
    "mu_R" -> {1, -1, -2},
    "c_E" -> {0, If[TrueQ[mutateCE], 2, 1], -1},
    "c_gamma" -> {0, 1, -1},
    "m" -> {1, 0, 0},
    "B_eff" -> {1, -1, -2},
    "M_h" -> {1, -1, 0},
    "K_h" -> {1, 1, -2},
    "C_hu" -> {1, 0, -2},
    "sigma" -> {2, 0, -4}
  |>;
  speedSquared = dimScale[2, dims["c_gamma"]];
  dims["K"] = dimAdd[
    dimAdd[speedSquared, dims["m"]],
    dimScale[-4, dims["rho"]]
  ];
  lockBDims = {
    dimAdd[dimScale[2, dims["c_E"]], dims["rho_br"]],
    dims["mu_R"]
  };
  lockADims = {
    dimAdd[dims["m"], dims["mu_R"]],
    dimAdd[
      dimAdd[dims["K"], dimScale[4, dims["rho"]]],
      dims["rho_br"]
    ]
  };
  stabilityDims = {
    dimAdd[dims["B_eff"], dims["K_h"]],
    dimScale[2, dims["C_hu"]],
    dims["sigma"]
  };
  homogeneous =
    SameQ @@ lockADims &&
    SameQ @@ lockBDims &&
    SameQ @@ stabilityDims;
  <|
    "lockA" -> lockADims,
    "lockB" -> lockBDims,
    "stability" -> stabilityDims,
    "homogeneous" -> homogeneous
  |>
];

baseEarned = {
  "c_s^2=5*K*rho^4/m",
  "c_gamma^2=mu_R/rho_br",
  "B_eff=rho_B0^2/chi_c",
  "K_h=M_h*c_E^2",
  "B_eff*K_h-C_hu^2=sigma"
};

partitionFromStatuses[statuses_Association] := Module[
  {earned, calibrated},
  earned = baseEarned;
  calibrated = {};
  Do[
    If[
      statuses[name] === "ENTAILED",
      earned = Union[earned, {"L_" <> name}]
    ];
    If[
      statuses[name] === "WITNESSED",
      calibrated = Union[calibrated, {"L_" <> name}]
    ],
    {name, {"A", "B"}}
  ];
  <|"earned" -> Sort[earned], "calibrated" -> Sort[calibrated]|>
];

partitionBiconditionalQ[
  statuses_Association,
  partition_Association
] := And @@ Flatten@Table[
  {
    MemberQ[partition["earned"], "L_" <> name] ===
      (statuses[name] === "ENTAILED"),
    MemberQ[partition["calibrated"], "L_" <> name] ===
      (statuses[name] === "WITNESSED")
  },
  {name, {"A", "B"}}
];


(* ---------------------------------------------------------------------- *)
(* Bounded, engine-qualified source-to-stage manifest.                     *)
(* ---------------------------------------------------------------------- *)

sourcePredicateTotal = 22;
sourcePredicateUniverse = {
  "pathA40.py::route_a_inventory",
  "pathA40.py::route_b_check",
  "pathA40.py::freedom_check",
  "pathA40.py::sat_decision(sector)",
  "pathA40.py::sat_decision(locks)",
  "pathA40.py::entailment_status(A)",
  "pathA40.py::entailment_status(B)",
  "pathA40.py::dimension_delta",
  "pathA40.py::decide",
  "pathA40.py+pathA40.wl::field_overlay",
  "pathA40.py::control.well_posed",
  "pathA40.py::control.absent",
  "pathA40.py::control.partial_inventory",
  "pathA40.py::control.forced_lock",
  "pathA40.py::control.a_only_partial",
  "pathA40.py::control.b_only_partial",
  "pathA40.py::control.over_constrained",
  "pathA40.py::control.freedom_tie",
  "pathA40.py::ledger.earned_equalities",
  "pathA40.py+pathA40.wl::filesystem_token_scans",
  "pathA40.py+pathA40.wl::harness_and_file_writing",
  "pathA40.py+pathA40.wl::cross_engine_payload_and_cross_read"
};
preservedSourceIds = Take[sourcePredicateUniverse, 18];
replacedSourceIds = {sourcePredicateUniverse[[19]]};
scopedSourceIds = Take[sourcePredicateUniverse, -3];
stageOnlyIds = {
  "stage040::CONDITIONALITY_FREEDOM",
  "stage040::DIMENSION_HOMOGENEITY",
  "stage040::DUAL_ENGINE_TERMS",
  "stage040::SOURCE_TO_STAGE_MANIFEST"
};

manifestRows[] := Join[
  ({#, "preserved-folded"} & /@ preservedSourceIds),
  ({#, "replaced-by-stronger"} & /@ replacedSourceIds),
  ({#, "scoped-out"} & /@ scopedSourceIds),
  ({#, "newly-added"} & /@ stageOnlyIds)
];

expectedManifestCategory = <|
  "pathA40.py::route_a_inventory" -> "preserved-folded",
  "pathA40.py::route_b_check" -> "preserved-folded",
  "pathA40.py::freedom_check" -> "preserved-folded",
  "pathA40.py::sat_decision(sector)" -> "preserved-folded",
  "pathA40.py::sat_decision(locks)" -> "preserved-folded",
  "pathA40.py::entailment_status(A)" -> "preserved-folded",
  "pathA40.py::entailment_status(B)" -> "preserved-folded",
  "pathA40.py::dimension_delta" -> "preserved-folded",
  "pathA40.py::decide" -> "preserved-folded",
  "pathA40.py+pathA40.wl::field_overlay" -> "preserved-folded",
  "pathA40.py::control.well_posed" -> "preserved-folded",
  "pathA40.py::control.absent" -> "preserved-folded",
  "pathA40.py::control.partial_inventory" -> "preserved-folded",
  "pathA40.py::control.forced_lock" -> "preserved-folded",
  "pathA40.py::control.a_only_partial" -> "preserved-folded",
  "pathA40.py::control.b_only_partial" -> "preserved-folded",
  "pathA40.py::control.over_constrained" -> "preserved-folded",
  "pathA40.py::control.freedom_tie" -> "preserved-folded",
  "pathA40.py::ledger.earned_equalities" ->
    "replaced-by-stronger",
  "pathA40.py+pathA40.wl::filesystem_token_scans" -> "scoped-out",
  "pathA40.py+pathA40.wl::harness_and_file_writing" -> "scoped-out",
  "pathA40.py+pathA40.wl::cross_engine_payload_and_cross_read" ->
    "scoped-out",
  "stage040::CONDITIONALITY_FREEDOM" -> "newly-added",
  "stage040::DIMENSION_HOMOGENEITY" -> "newly-added",
  "stage040::DUAL_ENGINE_TERMS" -> "newly-added",
  "stage040::SOURCE_TO_STAGE_MANIFEST" -> "newly-added"
|>;

manifestState[rows_List, liveRegistry_List] := Module[
  {
    sourceRows,
    identifiers,
    sourceIdentifiers,
    prefixes,
    valid
  },
  sourceRows = Select[
    rows,
    MemberQ[sourcePredicateUniverse, First[#]] &
  ];
  identifiers = First /@ rows;
  sourceIdentifiers = First /@ sourceRows;
  prefixes = Union[
    First@StringSplit[#, "::", 2] & /@ identifiers
  ];
  valid =
    Length[sourceRows] === sourcePredicateTotal &&
    Sort[sourceIdentifiers] === Sort[sourcePredicateUniverse] &&
    DuplicateFreeQ[sourceIdentifiers] &&
    DuplicateFreeQ[identifiers] &&
    Association[Rule @@@ rows] === expectedManifestCategory &&
    prefixes === Sort[
      {"pathA40.py", "pathA40.py+pathA40.wl", "stage040"}
    ] &&
    Length[liveRegistry] === 28;
  <|
    "valid" -> valid,
    "sourceCount" -> Length[sourceRows],
    "registrySize" -> Length[liveRegistry],
    "partition" -> Counts[Last /@ rows],
    "prefixes" -> prefixes
  |>
];


(* ---------------------------------------------------------------------- *)
(* The 28 executable teeth.                                                *)
(* ---------------------------------------------------------------------- *)

runAssertions[] := Module[
  {
    productionKinds,
    production,
    productionProvA,
    productionProvB,
    productionDimension,
    productionField,
    routeAKinds,
    routeAActual,
    routeBKinds,
    routeBActual,
    freedomKinds,
    freedomActual,
    sectorKinds,
    sectorActual,
    lockSatKinds,
    lockSatActual,
    provAKinds,
    provAActual,
    witnessA,
    witnessAValid,
    witnessAValue,
    provBKinds,
    provBActual,
    witnessB,
    witnessBValid,
    witnessBValue,
    witnessRCone,
    witnessRConeLockValue,
    dimensionActual,
    verdictCase,
    verdictPair,
    detField,
    detResidualOK,
    shearGuardOK,
    poleField,
    polesOK,
    openField,
    openOK,
    productionStatuses,
    productionPartition,
    forcedStatuses,
    forcedPartition,
    partitionOK,
    tiedKinds,
    tiedCase,
    conditionality,
    wellPosedKinds,
    wellPosedCase,
    absentKinds,
    absentCase,
    absentTuple,
    partialKinds,
    partialCase,
    partialTuple,
    controlSpecs,
    controlKinds,
    control,
    controlDelta,
    verdictRows,
    computedVerdictRows,
    actualCaseName,
    computedCase,
    dimState,
    localInventory,
    expectedInventory,
    inventoryOK,
    liveManifest,
    liveRegistry,
    manifestActual
  },
  productionKinds = caseKinds["production"];
  production = runCase[productionKinds];
  productionProvA = entailRecord[productionKinds, "A"];
  productionProvB = entailRecord[productionKinds, "B"];
  productionDimension = dimensionRecord[productionKinds];
  productionField = fieldOverlay[];

  section["Computed production inventory and algebra"];
  routeAKinds = If[
    activeMutation === "ROUTE_A_GRADE",
    Append[productionKinds, rKinds["R1"]],
    productionKinds
  ];
  routeAActual = routeA[routeAKinds];
  expectBool[
    "ROUTE_A_GRADE",
    routeAActual === <|
      "grade" ->
        "ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT",
      "missing" -> {"R1", "R2", "R3", "R4", "R5"}
    |>,
    routeAActual
  ];

  routeBKinds = If[
    activeMutation === "ROUTE_B_STATUS",
    Append[productionKinds, "thin_plate_over_import_relation"],
    productionKinds
  ];
  routeBActual = routeB[routeBKinds];
  expectBool[
    "ROUTE_B_STATUS",
    routeBActual === "ROUTE_B_CLOSED_CHECKED_NEGATIVE",
    routeBActual
  ];

  freedomKinds = If[
    activeMutation === "FREEDOM_STATUS",
    Append[productionKinds, "freedom_tie_relation"],
    productionKinds
  ];
  freedomActual = freedomRecord[freedomKinds];
  expectBool[
    "FREEDOM_STATUS",
    freedomActual === <|
      "status" -> "FREEDOM_UNCONSTRAINED{C_hu,rho_br}",
      "freeParams" -> {"C_hu", "rho_br"}
    |>,
    freedomActual
  ];

  sectorKinds = If[
    activeMutation === "SECTOR_SAT",
    Append[productionKinds, "over_stability_relation"],
    productionKinds
  ];
  sectorActual = satRecord[sectorKinds, {}]["status"];
  expectBool[
    "SECTOR_SAT",
    sectorActual === "SAT",
    sectorActual
  ];

  lockSatKinds = If[
    activeMutation === "LOCK_SAT",
    Append[productionKinds, "freedom_tie_relation"],
    productionKinds
  ];
  lockSatActual =
    satRecord[lockSatKinds, {"A", "B"}]["status"];
  expectBool[
    "LOCK_SAT",
    lockSatActual === "SAT",
    lockSatActual
  ];

  provAKinds = If[
    activeMutation === "PROVENANCE_LOCK_A",
    Append[productionKinds, "force_A_fake_relation"],
    productionKinds
  ];
  provAActual = entailRecord[provAKinds, "A"];
  expectBool[
    "PROVENANCE_LOCK_A",
    provAActual["status"] === "WITNESSED" &&
      provAActual["universal"] === False,
    provAActual
  ];

  witnessA = unlockedWitness[];
  If[
    activeMutation === "LOCK_A_WITNESS_VALUE",
    witnessA = Join[
      DeleteCases[witnessA, muR -> _],
      {muR -> 1}
    ]
  ];
  witnessAValid = witnessValidQ[productionKinds, {}, witnessA];
  witnessAValue = FullSimplify[lockA /. witnessA];
  expectBool[
    "LOCK_A_WITNESS_VALUE",
    witnessAValid && witnessAValue === 5,
    <|"value" -> witnessAValue, "valid" -> witnessAValid|>
  ];

  provBKinds = If[
    activeMutation === "PROVENANCE_LOCK_B",
    Append[productionKinds, "force_B_fake_relation"],
    productionKinds
  ];
  provBActual = entailRecord[provBKinds, "B"];
  expectBool[
    "PROVENANCE_LOCK_B",
    provBActual["status"] === "WITNESSED" &&
      provBActual["universal"] === False,
    provBActual
  ];

  witnessB = unlockedWitness[];
  If[
    activeMutation === "LOCK_B_WITNESS_VALUE",
    witnessB = Join[
      DeleteCases[witnessB, (cE | Kh | sigma) -> _],
      {cE -> 2, Kh -> 4, sigma -> 15/4}
    ]
  ];
  witnessBValid = witnessValidQ[productionKinds, {}, witnessB];
  witnessBValue = FullSimplify[lockB /. witnessB];
  witnessRCone = FullSimplify[
    cE^2*rhoBr/muR /. witnessB
  ];
  witnessRConeLockValue = FullSimplify[
    muR*(cE^2*rhoBr/muR - 1) /. witnessB
  ];
  expectBool[
    "LOCK_B_WITNESS_VALUE",
    witnessBValid &&
      witnessBValue === 7 &&
      witnessRCone === 9/2 &&
      witnessRCone =!= 1 &&
      witnessRConeLockValue === 7,
    <|
      "value" -> witnessBValue,
      "rCone" -> witnessRCone,
      "muR*(rCone-1)" -> witnessRConeLockValue,
      "valid" -> witnessBValid
    |>
  ];

  dimensionActual = dimensionRecord[
    productionKinds,
    If[activeMutation === "DELTA_R", {"A"}, {"A", "B"}]
  ];
  expectBool[
    "DELTA_R",
    dimensionActual["status"] === "CERTIFIED" &&
      dimensionActual["dimBefore"] === 10 &&
      dimensionActual["dimAfter"] === 8 &&
      dimensionActual["deltaR"] === 2,
    dimensionActual
  ];

  verdictCase = If[
    activeMutation === "PRODUCTION_VERDICT",
    runCase[caseKinds["b_only_partial"]],
    production
  ];
  verdictPair = {
    verdictCase["verdict"],
    verdictCase["riders"]
  };
  expectBool[
    "PRODUCTION_VERDICT",
    verdictPair === {"CONE_LOCK_CALIBRATED", {}},
    verdictPair
  ];

  section["Computed field overlay and OPEN_110"];
  detField = fieldOverlay[
    activeMutation === "FIELD_OVERLAY_DET_ON_CONE",
    False,
    False
  ];
  detResidualOK = TrueQ[
    FullSimplify[
      detField["detOnCone"] == -Chu^2*kWave^4
    ]
  ];
  shearGuardOK = TrueQ[
    FullSimplify[detField["shearResidual"] == 0]
  ];
  expectBool[
    "FIELD_OVERLAY_DET_ON_CONE",
    detResidualOK && shearGuardOK,
    <|
      "detOnCone" -> detField["detOnCone"],
      "detResidualPassed" -> detResidualOK,
      "directShearResidual" -> detField["shearResidual"]
    |>
  ];

  poleField = fieldOverlay[
    False,
    activeMutation === "FIELD_OVERLAY_POLES",
    False
  ];
  polesOK =
    Length[poleField["deltas"]] === 2 &&
    And @@ (
      Function[
        expected,
        AnyTrue[
          poleField["deltas"],
          TrueQ[
            FullSimplify[
              # == expected,
              Assumptions -> rhoBr > 0 && Mh > 0
            ]
          ] &
        ]
      ] /@ poleField["expectedDeltas"]
    );
  expectBool[
    "FIELD_OVERLAY_POLES",
    polesOK,
    poleField["deltas"]
  ];

  openField = fieldOverlay[
    False,
    False,
    activeMutation === "OPEN_110_CARRY"
  ];
  openOK = TrueQ[
    FullSimplify[
      Equivalent[
        openField["coincidenceCondition"],
        Chu == 0
      ]
    ]
  ];
  expectBool[
    "OPEN_110_CARRY",
    openOK,
    openField["coincidenceCondition"]
  ];

  section[
    "Earned/calibrated partition and freedom conditionality"
  ];
  productionStatuses = <|
    "A" -> productionProvA["status"],
    "B" -> productionProvB["status"]
  |>;
  productionPartition =
    partitionFromStatuses[productionStatuses];
  If[
    activeMutation === "EARNED_VS_CALIBRATED_PARTITION",
    productionPartition["earned"] = Union[
      productionPartition["earned"],
      {"L_A"}
    ]
  ];
  forcedStatuses = <|
    "A" -> entailRecord[caseKinds["forced_lock"], "A"]["status"],
    "B" -> entailRecord[caseKinds["forced_lock"], "B"]["status"]
  |>;
  forcedPartition = partitionFromStatuses[forcedStatuses];
  partitionOK =
    Intersection[
      productionPartition["earned"],
      productionPartition["calibrated"]
    ] === {} &&
    partitionBiconditionalQ[
      productionStatuses,
      productionPartition
    ] &&
    partitionBiconditionalQ[
      forcedStatuses,
      forcedPartition
    ] &&
    productionPartition["calibrated"] === {"L_A", "L_B"} &&
    productionPartition["earned"] === Sort[baseEarned] &&
    forcedPartition["earned"] ===
      Sort[Join[baseEarned, {"L_A", "L_B"}]] &&
    forcedPartition["calibrated"] === {};
  expectBool[
    "EARNED_VS_CALIBRATED_PARTITION",
    partitionOK,
    <|
      "productionStatuses" -> productionStatuses,
      "productionPartition" -> productionPartition,
      "forcedStatuses" -> forcedStatuses,
      "forcedPartition" -> forcedPartition
    |>
  ];

  tiedKinds = If[
    activeMutation === "CONDITIONALITY_FREEDOM",
    productionKinds,
    caseKinds["freedom_tie"]
  ];
  tiedCase = runCase[tiedKinds];
  conditionality =
    production["freedom"]["status"] ===
      "FREEDOM_UNCONSTRAINED{C_hu,rho_br}" &&
    production["verdict"] === "CONE_LOCK_CALIBRATED" &&
    tiedCase["freedom"]["status"] === "FREEDOM_TIED" &&
    tiedCase["verdict"] === "NO_GO(cone-lock)";
  expectBool[
    "CONDITIONALITY_FREEDOM",
    conditionality,
    <|
      "free" -> {
        production["freedom"]["status"],
        production["verdict"]
      },
      "tied" -> {
        tiedCase["freedom"]["status"],
        tiedCase["verdict"]
      }
    |>
  ];

  section["Individually computed falsifiability controls"];
  wellPosedKinds = caseKinds["well_posed"];
  If[
    activeMutation === "CTRL_WELL_POSED",
    wellPosedKinds = DeleteCases[
      wellPosedKinds,
      rKinds["R5"]
    ]
  ];
  wellPosedCase = runCase[wellPosedKinds];
  expectBool[
    "CTRL_WELL_POSED",
    wellPosedCase["verdict"] === "HALT_ROUTE_A_WELL_POSED",
    wellPosedCase["verdict"]
  ];

  absentKinds = caseKinds["absent"];
  If[
    activeMutation === "CTRL_ABSENT",
    absentKinds = Join[absentKinds, Values[rKinds]]
  ];
  absentCase = runCase[absentKinds];
  absentTuple = {
    absentCase["routeA"]["grade"],
    absentCase["routeA"]["missing"],
    absentCase["verdict"],
    Lookup[absentCase["dimension"], "deltaR", Null]
  };
  expectBool[
    "CTRL_ABSENT",
    absentTuple === {
      "ROUTE_A_ABSENT",
      {"R1", "R2", "R3", "R4", "R5"},
      "CONE_LOCK_CALIBRATED",
      2
    },
    absentTuple
  ];

  partialKinds = caseKinds["partial_inventory"];
  If[
    activeMutation === "CTRL_PARTIAL_INVENTORY",
    partialKinds = Join[
      partialKinds,
      {rKinds["R3"], rKinds["R4"], rKinds["R5"]}
    ]
  ];
  partialCase = runCase[partialKinds];
  partialTuple = {
    partialCase["routeA"]["grade"],
    partialCase["routeA"]["missing"],
    partialCase["verdict"],
    Lookup[partialCase["dimension"], "deltaR", Null]
  };
  expectBool[
    "CTRL_PARTIAL_INVENTORY",
    partialTuple === {
      "ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT",
      {"R3", "R4", "R5"},
      "CONE_LOCK_CALIBRATED",
      2
    },
    partialTuple
  ];

  controlSpecs = {
    {
      "CTRL_FORCED_LOCK",
      "forced_lock",
      "CONE_LOCK_DERIVED",
      0
    },
    {
      "CTRL_A_ONLY_PARTIAL",
      "a_only_partial",
      "CONE_LOCK_PARTIAL(derived=A, calibrated=B)",
      1
    },
    {
      "CTRL_B_ONLY_PARTIAL",
      "b_only_partial",
      "CONE_LOCK_PARTIAL(derived=B, calibrated=A)",
      1
    },
    {
      "CTRL_OVER_CONSTRAINED",
      "over_constrained",
      "NO_GO(sector-ledger)",
      Null
    },
    {
      "CTRL_FREEDOM_TIE",
      "freedom_tie",
      "NO_GO(cone-lock)",
      Null
    }
  };
  Do[
    controlKinds = If[
      activeMutation === spec[[1]],
      productionKinds,
      caseKinds[spec[[2]]]
    ];
    control = runCase[controlKinds];
    controlDelta =
      Lookup[control["dimension"], "deltaR", Null];
    expectBool[
      spec[[1]],
      control["verdict"] === spec[[3]] &&
        If[spec[[4]] === Null, True, controlDelta === spec[[4]]],
      <|
        "verdict" -> control["verdict"],
        "deltaR" -> controlDelta,
        "sectorSat" -> control["sectorSat"],
        "lockSat" -> control["lockSat"]
      |>
    ],
    {spec, controlSpecs}
  ];

  section[
    "Verdict table, dimensions, local inventory, and manifest"
  ];
  verdictRows = {
    {
      "production",
      "production",
      {"CONE_LOCK_CALIBRATED", {}}
    },
    {
      "well_posed",
      "well_posed",
      {"HALT_ROUTE_A_WELL_POSED", {}}
    },
    {
      "absent",
      "absent",
      {"CONE_LOCK_CALIBRATED", {}}
    },
    {
      "partial_inventory",
      "partial_inventory",
      {"CONE_LOCK_CALIBRATED", {}}
    },
    {
      "forced_lock",
      "forced_lock",
      {"CONE_LOCK_DERIVED", {}}
    },
    {
      "a_only_partial",
      "a_only_partial",
      {
        "CONE_LOCK_PARTIAL(derived=A, calibrated=B)",
        {}
      }
    },
    {
      "b_only_partial",
      "b_only_partial",
      {
        "CONE_LOCK_PARTIAL(derived=B, calibrated=A)",
        {}
      }
    },
    {
      "over_constrained",
      "over_constrained",
      {"NO_GO(sector-ledger)", {}}
    },
    {
      "freedom_tie",
      "freedom_tie",
      {"NO_GO(cone-lock)", {}}
    },
    {
      "a_inconclusive",
      "a_inconclusive",
      {
        "CONE_LOCK_PROVENANCE_INCONCLUSIVE",
        {"ENTAILMENT_INCONCLUSIVE(L_A)"}
      }
    }
  };
  computedVerdictRows = Table[
    actualCaseName = If[
      activeMutation === "VERDICT_REDERIVATION" &&
        row[[1]] === "forced_lock",
      "a_only_partial",
      row[[2]]
    ];
    computedCase = runCase[caseKinds[actualCaseName]];
    {
      row[[1]],
      {computedCase["verdict"], computedCase["riders"]},
      row[[3]]
    },
    {row, verdictRows}
  ];
  expectBool[
    "VERDICT_REDERIVATION",
    And @@ (#[[2]] === #[[3]] & /@ computedVerdictRows),
    computedVerdictRows
  ];

  dimState = dimensionState[
    activeMutation === "DIMENSION_HOMOGENEITY"
  ];
  expectBool[
    "DIMENSION_HOMOGENEITY",
    dimState["homogeneous"],
    dimState
  ];

  localInventory = <|
    "verdict" -> production["verdict"],
    "riders" -> production["riders"],
    "deltaR" -> productionDimension["deltaR"],
    "dimBefore" -> productionDimension["dimBefore"],
    "dimAfter" -> productionDimension["dimAfter"],
    "provA" -> productionProvA["status"],
    "provB" -> productionProvB["status"],
    "valueA" -> productionProvA["value"],
    "valueB" -> productionProvB["value"],
    "freedom" -> production["freedom"]["status"],
    "routeA" -> {
      production["routeA"]["grade"],
      production["routeA"]["missing"]
    },
    "routeB" -> production["routeB"],
    "detOnCone" -> productionField["detOnCone"],
    "poles" -> productionField["deltas"]
  |>;
  If[
    activeMutation === "DUAL_ENGINE_TERMS",
    KeyDropFrom[localInventory, "valueB"]
  ];
  expectedInventory = <|
    "verdict" -> "CONE_LOCK_CALIBRATED",
    "riders" -> {},
    "deltaR" -> 2,
    "dimBefore" -> 10,
    "dimAfter" -> 8,
    "provA" -> "WITNESSED",
    "provB" -> "WITNESSED",
    "valueA" -> 5,
    "valueB" -> 7,
    "freedom" -> "FREEDOM_UNCONSTRAINED{C_hu,rho_br}",
    "routeA" -> {
      "ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT",
      {"R1", "R2", "R3", "R4", "R5"}
    },
    "routeB" -> "ROUTE_B_CLOSED_CHECKED_NEGATIVE",
    "detOnCone" -> -Chu^2*kWave^4,
    "poles" -> productionField["expectedDeltas"]
  |>;
  inventoryOK =
    Sort[Keys[localInventory]] === Sort[Keys[expectedInventory]] &&
    And @@ Table[
      If[
        key === "poles",
        Length[localInventory[key]] === 2 &&
          And @@ (
            Function[
              expected,
              AnyTrue[
                localInventory[key],
                TrueQ[
                  FullSimplify[
                    # == expected,
                    Assumptions -> rhoBr > 0 && Mh > 0
                  ]
                ] &
              ]
            ] /@ expectedInventory[key]
          ),
        TrueQ[
          FullSimplify[
            localInventory[key] == expectedInventory[key]
          ]
        ]
      ],
      {key, Keys[expectedInventory]}
    ];
  expectBool[
    "DUAL_ENGINE_TERMS",
    inventoryOK,
    localInventory
  ];

  liveManifest = manifestRows[];
  liveRegistry = toothOrder;
  If[
    activeMutation === "SOURCE_TO_STAGE_MANIFEST",
    liveManifest = Select[
      liveManifest,
      First[#] =!= sourcePredicateUniverse[[1]] &
    ];
    liveManifest = Replace[
      liveManifest,
      {
        sourcePredicateUniverse[[2]],
        _
      } :> {
        sourcePredicateUniverse[[2]],
        "scoped-out"
      },
      {1}
    ];
    liveRegistry = Most[toothOrder]
  ];
  manifestActual = manifestState[liveManifest, liveRegistry];
  expectBool[
    "SOURCE_TO_STAGE_MANIFEST",
    manifestActual["valid"],
    manifestActual
  ];

  Print[""];
  Print["VERDICT=", production["verdict"]];
  Print[
    "INTERPRETATION=CALIBRATED/UNCOMMITTED; NO committed cone lock; both locks are available calibration choices, not earned facts"
  ];
  Print[
    "DIMENSION=",
    productionDimension["dimBefore"],
    "->",
    productionDimension["dimAfter"],
    "; delta_r=",
    productionDimension["deltaR"]
  ];
  Print[
    "LOCK_PROVENANCE=L_A:",
    productionProvA["status"],
    "@",
    productionProvA["value"],
    ";L_B:",
    productionProvB["status"],
    "@",
    productionProvB["value"]
  ];
  Print[
    "R_CONE=OPEN_R71; witness=9/2; c_E=c_gamma is NOT settled"
  ];
  Print[
    "OPEN_110=OFF_CONE_under_AB_proportional_to_C_hu_squared_OPEN_110"
  ];
  Print[
    "CONDITIONAL_PROVENANCE=FREEDOM_CERTIFIED_CURRENT_LEDGER{C_hu} [stage041 citation only];FREEDOM_SIM_DEFERRED{rho_br}; no stage041 cross-read"
  ];
  Print[
    "ROUTE_A_DEBT=REGISTERED_DEFERRED nonlinear throat solve missing {R1,R2,R3,R4,R5}"
  ];
  Print["SOURCE_PREDICATE_TOTAL=22"];
  Print["EXECUTABLE_TOOTH_TOTAL=28"];
  production["verdict"]
];


ok = Catch[
  If[
    Length[toothOrder] =!= 28,
    raise["TOOTH_REGISTRY_DECLARATION"]
  ];
  If[
    activeMutation =!= "" &&
      ! MemberQ[toothOrder, activeMutation],
    Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
    Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
    raise["UNKNOWN_MUTATION"]
  ];

  Print[
    "ledger_stage040_cone_lock_readjudication Mathematica audit"
  ];
  Print[
    "ROUTE=Resolve Exists/ForAll + FindInstance + CAD RegionDimension[ImplicitRegion] + FullSimplify/Solve"
  ];
  Print["FILE_IO=none; CROSS_ENGINE_COMPARE=none"];
  If[
    activeMutation =!= "",
    Print["ACTIVE_MUTATION=", activeMutation];
    Print[
      "MUTATED_PRIMITIVE=",
      ablationDescriptions[activeMutation]
    ]
  ];

  finalVerdict = runAssertions[];
  If[
    activeMutation =!= "",
    Print["FIRST_FAILURE=MUTATION_DID_NOT_FIRE"];
    raise["MUTATION_DID_NOT_FIRE"]
  ];
  True,
  "ledgerStage040Failure",
  Function[{message, tag}, False]
];

Print[""];
Print["TOOTH_COUNT=", Length[toothOrder]];
Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
If[
  TrueQ[ok],
  Print[
    "OVERALL PASS: Mathematica independently reached ",
    finalVerdict
  ];
  Exit[0],
  Print[
    "OVERALL FAIL: Mathematica stage040 audit did not close"
  ];
  Exit[1]
]
