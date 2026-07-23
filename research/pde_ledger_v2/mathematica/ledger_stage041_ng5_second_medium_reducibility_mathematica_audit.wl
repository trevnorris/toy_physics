(* Ledger stage041 Mathematica audit: NG5 one-medium reducibility gate.

   Standalone, print-only, assert-zero, and cross-engine-file-I/O-free.

   This engine is deliberately organized differently from both the stage
   SymPy engine and the source Mathematica mirror.  Rows are partitioned by
   computed category sets; GroupBy assembles origin groups, Select builds the
   active subsets, and SubsetQ evaluates joint route targets.  FirstCase
   resolves ordered category/reason tables.  There is no per-row Which
   classifier, no hardcoded production result, and no cross-engine read.

   Tooth-local runtime ablation uses LEDGER_STAGE041_MUTATION.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE041_MUTATION";
activeMutation = Quiet@Check[Environment[mutationEnvironment], ""];
If[! StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];

routeConjuncts = {
  "required_fields_present",
  "named_solve_in_provenance",
  "result_status_allowed",
  "finite_listed_missing_objects",
  "target_identity",
  "joint_target_identity",
  "target_blind",
  "falsifiers_clear"
};

conjunctTeeth = {
  "CTRL_ROUTE_CONJUNCT_REQUIRED_FIELDS_PRESENT",
  "CTRL_ROUTE_CONJUNCT_NAMED_SOLVE_IN_PROVENANCE",
  "CTRL_ROUTE_CONJUNCT_RESULT_STATUS_ALLOWED",
  "CTRL_ROUTE_CONJUNCT_FINITE_LISTED_MISSING_OBJECTS",
  "CTRL_ROUTE_CONJUNCT_TARGET_IDENTITY",
  "CTRL_ROUTE_CONJUNCT_JOINT_TARGET_IDENTITY",
  "CTRL_ROUTE_CONJUNCT_TARGET_BLIND",
  "CTRL_ROUTE_CONJUNCT_FALSIFIERS_CLEAR"
};

toothOrder = {
  "ORIGIN_CLASSIFICATION",
  "FULL_ORIGIN_LEDGER",
  "FULL_ROUTE_EVALUATION_LEDGER",
  "PRODUCTION_VERDICT",
  "DRIFT_COUNT",
  "FREEDOM_STATE_RHO_BR",
  "FREEDOM_STATE_C_HU",
  "LOCATION_CLOSURE",
  "ARENA_LABELS",
  "LINEAGE_FINDING",
  "REDUCTION_STATUS",
  "ANTI_ABSORPTION_DRIFT_SIDE",
  "CTRL_AB_DELETE_REGISTRY",
  "CTRL_ROUTE_BLANK_ROUTE_A",
  "CTRL_ROUTE_FIELD_BLANK_TARGET_BLIND",
  "CTRL_ROUTE_FIELD_BLANK_MISSING_OBJECTS",
  Sequence @@ conjunctTeeth,
  "CTRL_CALIBRATION_ABLATION_Q_E",
  "CTRL_IRREDUCIBLE_SYNTHETIC",
  "CTRL_REDUCIBLE_DERIVED_SYNTHETIC",
  "CTRL_CONTRADICTION",
  "CTRL_ROUTE_EVAL_RECORDED",
  "CTRL_RESIDUAL_MULTIPLIER_ABLATION",
  "CTRL_LOCATION_CLOSURE_OUT_OF_ARENA",
  "DIMENSION_HOMOGENEITY",
  "DUAL_ENGINE_TERMS",
  "VERDICT_REDERIVATION",
  "SOURCE_TO_STAGE_MANIFEST"
};

ablationDescriptions = <|
  "ORIGIN_CLASSIFICATION" ->
    "make Future-Embedding-Overlap a valid deferred route before category partitioning",
  "FULL_ORIGIN_LEDGER" ->
    "redirect B_eff's dependency to a non-irreducible row before category partitioning",
  "FULL_ROUTE_EVALUATION_LEDGER" ->
    "give c_gamma an absent route hint so its computed reason changes while its origin does not",
  "PRODUCTION_VERDICT" ->
    "make C_hu reducible-deferred in the private verdict witness",
  "DRIFT_COUNT" ->
    "inject xi_active into the private category sets feeding the count",
  "FREEDOM_STATE_RHO_BR" ->
    "remove Route-A in the private rho_br freedom-state witness",
  "FREEDOM_STATE_C_HU" ->
    "make Future-Embedding-Overlap fully valid and REGISTERED_DEFERRED",
  "LOCATION_CLOSURE" ->
    "inject a row labelled unassigned",
  "ARENA_LABELS" ->
    "relabel C_hu to the still-admissible 3D brane surface",
  "LINEAGE_FINDING" ->
    "replace production lineage facts by one same active object with a free residual multiplier",
  "REDUCTION_STATUS" ->
    "drop c_E while constructing the reduction-status pending map",
  "ANTI_ABSORPTION_DRIFT_SIDE" ->
    "rename the surviving historical rho_br entry to rho_B0 before retirement and projection",
  "CTRL_AB_DELETE_REGISTRY" ->
    "neutralize the all-route deletion control",
  "CTRL_ROUTE_BLANK_ROUTE_A" ->
    "neutralize the Route-A deletion control",
  "CTRL_ROUTE_FIELD_BLANK_TARGET_BLIND" ->
    "neutralize the target_blind field-deletion control",
  "CTRL_ROUTE_FIELD_BLANK_MISSING_OBJECTS" ->
    "neutralize the missing_objects field-deletion control",
  "CTRL_ROUTE_CONJUNCT_REQUIRED_FIELDS_PRESENT" ->
    "neutralize deletion of the otherwise-unused required source field",
  "CTRL_ROUTE_CONJUNCT_NAMED_SOLVE_IN_PROVENANCE" ->
    "neutralize the named-solve=False isolating control",
  "CTRL_ROUTE_CONJUNCT_RESULT_STATUS_ALLOWED" ->
    "neutralize the rejected-status isolating control",
  "CTRL_ROUTE_CONJUNCT_FINITE_LISTED_MISSING_OBJECTS" ->
    "neutralize the non-list missing_objects isolating control",
  "CTRL_ROUTE_CONJUNCT_TARGET_IDENTITY" ->
    "neutralize dropping only c_E from Route-A targets",
  "CTRL_ROUTE_CONJUNCT_JOINT_TARGET_IDENTITY" ->
    "neutralize the sentinel joint target",
  "CTRL_ROUTE_CONJUNCT_TARGET_BLIND" ->
    "neutralize target_blind=False while retaining the field",
  "CTRL_ROUTE_CONJUNCT_FALSIFIERS_CLEAR" ->
    "neutralize insertion of one Route-A falsifier",
  "CTRL_CALIBRATION_ABLATION_Q_E" ->
    "neutralize removal of the Q_E calibration anchor",
  "CTRL_IRREDUCIBLE_SYNTHETIC" ->
    "suppress xi_active injection",
  "CTRL_REDUCIBLE_DERIVED_SYNTHETIC" ->
    "change p_syn local_status away from SYNTHETIC_EARNED_RELATION",
  "CTRL_CONTRADICTION" ->
    "omit the cone locks so freedom_tie has an exact positive SAT witness",
  "CTRL_ROUTE_EVAL_RECORDED" ->
    "drop c_gamma's RouteEvaluation record from the active-record projection",
  "CTRL_RESIDUAL_MULTIPLIER_ABLATION" ->
    "close the free branch residual multiplier to 1",
  "CTRL_LOCATION_CLOSURE_OUT_OF_ARENA" ->
    "suppress loc_sentinel injection",
  "DIMENSION_HOMOGENEITY" ->
    "change the cited speed dimension from L T^-1 to L^2 T^-1",
  "DUAL_ENGINE_TERMS" ->
    "drop the locally computed lock_B witness from the canonical inventory",
  "VERDICT_REDERIVATION" ->
    "leave one future route's named-solve provenance false in the conditional witness",
  "SOURCE_TO_STAGE_MANIFEST" ->
    "drop one source predicate and mis-scope a second predicate"
|>;

If[Length[toothOrder] =!= 35,
  Print["FAIL: stage041 tooth declaration is not exactly 35"];
  Exit[1]
];

raise[message_] := Throw[message, "ledgerStage041Failure"];

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

allowedStatuses = {"SOLVED_PASS", "REGISTERED_DEFERRED"};
rejectedStatuses = {"FAILED", "BY_TUNING", "ABSENT", "PROMISSORY_ONLY"};
arenas = {"4D bulk", "3D brane surface", "throat/embedding seam"};
expectedTrio = {"rho_B0", "chi_c", "C_hu"};
routeATargets = {"rho_br", "mu_R", "c_E"};
expectedVerdict =
  "SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu})";


(* ---------------------------------------------------------------------- *)
(* Source facts as neutral rows and route records.                         *)
(* ---------------------------------------------------------------------- *)

makeRow[
  param_String,
  localStatus_String,
  options___Rule
] := Join[
  <|
    "param" -> param,
    "localStatus" -> localStatus,
    "active" -> False,
    "outOfActive" -> False,
    "dependencies" -> {},
    "routeHint" -> None,
    "calibrationAnchor" -> None,
    "calibratedGeometry" -> False,
    "location" -> "unassigned"
  |>,
  Association[options]
];

baseRoutes[] := {
  <|
    "name" -> "Route-A",
    "source" ->
      "pathA_40:ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT",
    "result_status" -> "REGISTERED_DEFERRED",
    "named_solve_in_provenance" -> True,
    "missing_objects" -> {"R1", "R2", "R3", "R4", "R5"},
    "targets" -> {"rho_br", "mu_R", "c_E"},
    "required_joint_targets" -> {"rho_br", "mu_R"},
    "target_blind" -> True,
    "falsifiers" -> {}
  |>,
  <|
    "name" -> "Future-Compression-4D-to-3D",
    "source" -> "named_future_reduction_route:not_currently_registered",
    "result_status" -> "PROMISSORY_ONLY",
    "named_solve_in_provenance" -> False,
    "missing_objects" ->
      {"compression-sector 4D-to-3D nonlinear brane solve"},
    "targets" -> {"rho_B0", "chi_c"},
    "required_joint_targets" -> {"rho_B0", "chi_c"},
    "target_blind" -> True,
    "falsifiers" -> {}
  |>,
  <|
    "name" -> "Future-Embedding-Overlap",
    "source" -> "named_future_reduction_route:not_currently_registered",
    "result_status" -> "PROMISSORY_ONLY",
    "named_solve_in_provenance" -> False,
    "missing_objects" ->
      {"embedding-overlap nonlinear throat solve deriving C_hu"},
    "targets" -> {"C_hu"},
    "required_joint_targets" -> {"C_hu"},
    "target_blind" -> True,
    "falsifiers" -> {}
  |>
};

updateRoute[routes_List, routeName_String, transform_] :=
  (If[#["name"] === routeName, transform[#], #] &) /@ routes;

validFutureRoutes[futureStatus_String, routeAStatus_String] := Module[
  {routes = baseRoutes[]},
  routes = updateRoute[
    routes,
    "Future-Compression-4D-to-3D",
    Function[fact, Join[
      fact,
      <|
        "named_solve_in_provenance" -> True,
        "result_status" -> futureStatus
      |>
    ]]
  ];
  routes = updateRoute[
    routes,
    "Future-Embedding-Overlap",
    Function[fact, Join[
      fact,
      <|
        "named_solve_in_provenance" -> True,
        "result_status" -> futureStatus
      |>
    ]]
  ];
  updateRoute[
    routes,
    "Route-A",
    Function[fact, Join[
      fact,
      <|
        "named_solve_in_provenance" -> True,
        "result_status" -> routeAStatus
      |>
    ]]
  ]
];

buildRows[options_: <||>] := Module[
  {
    bEffMutation,
    cGammaRouteMutation,
    cHuArenaMutation,
    addXi,
    addPSyn,
    pSynStatusMutation,
    addLocationSentinel,
    rows,
    drivers
  },
  bEffMutation = TrueQ[Lookup[options, "bEffDependencyMutation", False]];
  cGammaRouteMutation =
    TrueQ[Lookup[options, "cGammaRouteMutation", False]];
  cHuArenaMutation =
    TrueQ[Lookup[options, "cHuArenaMutation", False]];
  addXi = TrueQ[Lookup[options, "addXi", False]];
  addPSyn = TrueQ[Lookup[options, "addPSyn", False]];
  pSynStatusMutation =
    TrueQ[Lookup[options, "pSynStatusMutation", False]];
  addLocationSentinel =
    TrueQ[Lookup[options, "addLocationSentinel", False]];

  rows = {
    makeRow["rho", "BASE_SUBSTRATE", "location" -> "4D bulk"],
    makeRow["K", "BASE_SUBSTRATE", "location" -> "4D bulk"],
    makeRow["m", "BASE_SUBSTRATE", "location" -> "4D bulk"],
    makeRow["a", "BASE_SUBSTRATE", "location" -> "4D bulk"],
    makeRow[
      "ell_g",
      "ACCEPTED_GEOMETRY_SUBSTRATE",
      "location" -> "throat/embedding seam"
    ],
    makeRow[
      "g_ell(w)",
      "ACCEPTED_PROFILE_GIVEN_ell_g",
      "location" -> "throat/embedding seam"
    ],
    makeRow[
      "varrho_br[rho]",
      "PATHA25_CLOSED_DENSITY_SMECTIC_OBJECT",
      "outOfActive" -> True,
      "location" -> "3D brane surface"
    ],
    makeRow[
      "Sigma_n[rho]",
      "PATHA25_CLOSED_DENSITY_SMECTIC_OBJECT",
      "outOfActive" -> True,
      "location" -> "3D brane surface"
    ],
    makeRow[
      "delta_Sigma[rho]",
      "PATHA25_CLOSED_DENSITY_SMECTIC_OBJECT",
      "outOfActive" -> True,
      "location" -> "3D brane surface"
    ],
    makeRow[
      "rho_br",
      "POSTULATED_SHEAR_SURFACE_REDUCTION_TARGET",
      "active" -> True,
      "routeHint" -> "Route-A",
      "location" -> "3D brane surface"
    ],
    makeRow[
      "mu_R",
      "POSTULATED_SHEAR_SURFACE_REDUCTION_TARGET",
      "active" -> True,
      "routeHint" -> "Route-A",
      "location" -> "3D brane surface"
    ],
    makeRow[
      "c_E",
      "ROUTE_A_CONE_LOCK_REDUCTION_TARGET",
      "active" -> True,
      "routeHint" -> "Route-A",
      "location" -> "throat/embedding seam"
    ],
    makeRow[
      "c_gamma",
      "DEPENDENT_SPEED",
      "active" -> True,
      "dependencies" -> {"rho_br", "mu_R"},
      "routeHint" ->
        If[cGammaRouteMutation, "Absent-C-Gamma-Route", None],
      "location" -> "3D brane surface"
    ],
    makeRow[
      "rho_B0",
      "NO_CURRENT_REGISTERED_REDUCTION",
      "active" -> True,
      "routeHint" -> "Future-Compression-4D-to-3D",
      "location" -> "3D brane surface"
    ],
    makeRow[
      "chi_c",
      "NO_CURRENT_REGISTERED_REDUCTION",
      "active" -> True,
      "routeHint" -> "Future-Compression-4D-to-3D",
      "location" -> "3D brane surface"
    ],
    makeRow[
      "B_eff",
      "DEPENDENT_ON_COMPRESSION_INPUTS",
      "active" -> True,
      "dependencies" ->
        If[bEffMutation, {"c_gamma"}, {"rho_B0", "chi_c"}],
      "location" -> "3D brane surface"
    ],
    makeRow[
      "C_hu",
      "NO_CURRENT_REGISTERED_REDUCTION",
      "active" -> True,
      "routeHint" -> "Future-Embedding-Overlap",
      "location" ->
        If[
          cHuArenaMutation,
          "3D brane surface",
          "throat/embedding seam"
        ]
    ],
    makeRow[
      "Q_E",
      "DECLARED_CALIBRATION_ANCHOR",
      "active" -> True,
      "calibrationAnchor" -> "Q_E",
      "location" -> "throat/embedding seam"
    ],
    makeRow[
      "ell",
      "CALIBRATED_GEOMETRY_SOURCE_INPUT",
      "active" -> True,
      "calibratedGeometry" -> True,
      "location" -> "throat/embedding seam"
    ],
    makeRow[
      "b",
      "CALIBRATED_GEOMETRY_SOURCE_INPUT",
      "active" -> True,
      "calibratedGeometry" -> True,
      "location" -> "throat/embedding seam"
    ],
    makeRow[
      "M_h",
      "CALIBRATED_GEOMETRY_SOURCE_INPUT",
      "active" -> True,
      "calibratedGeometry" -> True,
      "location" -> "throat/embedding seam"
    ],
    makeRow[
      "K_h",
      "DEPENDENT_STIFFNESS",
      "active" -> True,
      "dependencies" -> {"M_h", "c_E"},
      "location" -> "throat/embedding seam"
    ],
    makeRow[
      "q_h",
      "DEPENDENT_PROJECTION",
      "active" -> True,
      "dependencies" -> {"Q_E", "b", "ell"},
      "location" -> "throat/embedding seam"
    ]
  };

  drivers = {
    "c_L1", "c_L2", "A_R", "k_R", "lambda_Cdiv", "chi_Cpin",
    "J_Pu", "kappa_Pu"
  };
  rows = Join[
    rows,
    (
      makeRow[
        #,
        "PATHA25_DRIVER_OUT_OF_ACTIVE_SHEAR_SURFACE_NG5",
        "outOfActive" -> True,
        "location" -> "3D brane surface"
      ] &
    ) /@ drivers,
    {
      makeRow[
        "lambda_Pu",
        "FROZEN_UNUSED_OUT_OF_ACTIVE_NG5",
        "outOfActive" -> True,
        "location" -> "3D brane surface"
      ],
      makeRow[
        "Omega_w",
        "FROZEN_UNUSED_OUT_OF_ACTIVE_NG5",
        "outOfActive" -> True,
        "location" -> "3D brane surface"
      ],
      makeRow[
        "lambda_N",
        "WALL_INTERNAL_EXCLUDED_FROM_ACTIVE_NG5",
        "outOfActive" -> True,
        "location" -> "throat/embedding seam"
      ],
      makeRow[
        "lambda_tau",
        "WALL_INTERNAL_EXCLUDED_FROM_ACTIVE_NG5",
        "outOfActive" -> True,
        "location" -> "throat/embedding seam"
      ],
      makeRow[
        "Nu",
        "SOURCE_NORMALIZATION_OUTSIDE_ACTIVE_ORIGIN_COMPETITION",
        "outOfActive" -> True,
        "location" -> "throat/embedding seam"
      ],
      makeRow[
        "a_T",
        "SOURCE_NORMALIZATION_OUTSIDE_ACTIVE_ORIGIN_COMPETITION",
        "outOfActive" -> True,
        "location" -> "throat/embedding seam"
      ],
      makeRow[
        "a_Tp",
        "SOURCE_NORMALIZATION_OUTSIDE_ACTIVE_ORIGIN_COMPETITION",
        "outOfActive" -> True,
        "location" -> "throat/embedding seam"
      ],
      makeRow[
        "a_L",
        "SOURCE_NORMALIZATION_OUTSIDE_ACTIVE_ORIGIN_COMPETITION",
        "outOfActive" -> True,
        "location" -> "throat/embedding seam"
      ]
    }
  ];

  If[
    addXi,
    rows = Append[
      rows,
      makeRow[
        "xi_active",
        "SYNTHETIC_ACTIVE_INPUT",
        "active" -> True,
        "location" -> "3D brane surface"
      ]
    ]
  ];
  If[
    addPSyn,
    rows = Append[
      rows,
      makeRow[
        "p_syn",
        If[
          pSynStatusMutation,
          "SYNTHETIC_ACTIVE_INPUT",
          "SYNTHETIC_EARNED_RELATION"
        ],
        "active" -> True,
        "location" -> "4D bulk"
      ]
    ]
  ];
  If[
    addLocationSentinel,
    rows = Append[
      rows,
      makeRow[
        "loc_sentinel",
        "BASE_SUBSTRATE",
        "location" -> "unassigned"
      ]
    ]
  ];
  rows
];


(* ---------------------------------------------------------------------- *)
(* Route evaluation and category-set classifier.                           *)
(* ---------------------------------------------------------------------- *)

routeEvaluate[row_Association, routes_List] := Module[
  {
    routeName,
    selected,
    fact,
    required,
    missingFields,
    status,
    missingObjects,
    targets,
    requiredJoint,
    falsifiers,
    checks,
    valid,
    reasonTable,
    reason
  },
  routeName = row["routeHint"];
  If[
    routeName === None,
    Return[
      <|
        "route_name" -> None,
        "valid" -> False,
        "reason" -> "NO_REGISTERED_ROUTE_FOR_PARAM",
        "result_status" -> None,
        "checks" -> <||>
      |>
    ]
  ];
  selected = Select[routes, #["name"] === routeName &];
  If[
    selected === {},
    Return[
      <|
        "route_name" -> routeName,
        "valid" -> False,
        "reason" -> "ROUTE_BLANKED_OR_ABSENT",
        "result_status" -> None,
        "checks" -> <|"mapped_route_exists" -> False|>
      |>
    ]
  ];
  fact = First[selected];
  required = {
    "name",
    "source",
    "result_status",
    "named_solve_in_provenance",
    "missing_objects",
    "targets",
    "target_blind",
    "falsifiers"
  };
  missingFields = Select[required, ! KeyExistsQ[fact, #] &];
  status = Lookup[fact, "result_status", None];
  missingObjects = Lookup[fact, "missing_objects", None];
  targets = Lookup[fact, "targets", None];
  requiredJoint = Lookup[fact, "required_joint_targets", {}];
  falsifiers = Lookup[fact, "falsifiers", None];
  checks = <|
    "required_fields_present" -> (missingFields === {}),
    "named_solve_in_provenance" ->
      TrueQ[Lookup[fact, "named_solve_in_provenance", False]],
    "result_status_allowed" -> MemberQ[allowedStatuses, status],
    "finite_listed_missing_objects" -> ListQ[missingObjects],
    "target_identity" ->
      (ListQ[targets] && MemberQ[targets, row["param"]]),
    "joint_target_identity" ->
      (ListQ[targets] && SubsetQ[targets, requiredJoint]),
    "target_blind" -> TrueQ[Lookup[fact, "target_blind", False]],
    "falsifiers_clear" ->
      (ListQ[falsifiers] && falsifiers === {})
  |>;
  valid = And @@ Lookup[checks, routeConjuncts];
  reasonTable = {
    {valid, "VALID_REGISTERED_ROUTE"},
    {
      MemberQ[rejectedStatuses, status],
      "REJECTED_RESULT_STATUS_" <> ToString[status]
    },
    {missingFields =!= {}, "MISSING_REQUIRED_ROUTE_FIELD"},
    {True, "ROUTE_EVALUATION_FAILED"}
  };
  reason = FirstCase[reasonTable, {True, value_} :> value];
  <|
    "route_name" -> routeName,
    "valid" -> valid,
    "reason" -> reason,
    "result_status" -> status,
    "checks" -> checks
  |>
];

paramsOf[records_List] := If[
  records === {},
  {},
  Lookup[records, "param"]
];

classifyByCategorySets[
  rows_List,
  routes_List,
  anchors_List,
  geometry_List
] := Module[
  {
    activeRows,
    activeParams,
    allParams,
    routeEvaluations,
    outSet,
    baseSet,
    geometryBaseSet,
    profileSet,
    dependencyRows,
    dependencySet,
    anchorSet,
    calibratedGeometrySet,
    earnedSet,
    solvedSet,
    deferredSet,
    independentSet,
    dependentOnIrreducibleSet,
    ordinaryDependentSet,
    partitions,
    originAssociation,
    groupedRows,
    activeIrreducible
  },
  activeRows = Select[rows, TrueQ[#["active"]] &];
  activeParams = paramsOf[activeRows];
  allParams = paramsOf[rows];
  routeEvaluations = Association[
    (#["param"] -> routeEvaluate[#, routes] &) /@ activeRows
  ];

  outSet = paramsOf[Select[rows, TrueQ[#["outOfActive"]] &]];
  baseSet = paramsOf[
    Select[rows, #["localStatus"] === "BASE_SUBSTRATE" &]
  ];
  geometryBaseSet = paramsOf[
    Select[
      rows,
      #["localStatus"] === "ACCEPTED_GEOMETRY_SUBSTRATE" &
    ]
  ];
  profileSet = paramsOf[
    Select[
      rows,
      #["localStatus"] === "ACCEPTED_PROFILE_GIVEN_ell_g" &
    ]
  ];
  dependencyRows = Select[activeRows, #["dependencies"] =!= {} &];
  dependencySet = paramsOf[dependencyRows];
  anchorSet = paramsOf[
    Select[
      activeRows,
      #["calibrationAnchor"] =!= None &&
        MemberQ[anchors, #["calibrationAnchor"]] &
    ]
  ];
  calibratedGeometrySet = paramsOf[
    Select[
      activeRows,
      TrueQ[#["calibratedGeometry"]] &&
        MemberQ[geometry, #["param"]] &
    ]
  ];
  earnedSet = paramsOf[
    Select[
      activeRows,
      #["localStatus"] === "SYNTHETIC_EARNED_RELATION" &
    ]
  ];
  solvedSet = Select[
    activeParams,
    TrueQ[routeEvaluations[#]["valid"]] &&
      routeEvaluations[#]["result_status"] === "SOLVED_PASS" &
  ];
  deferredSet = Select[
    activeParams,
    TrueQ[routeEvaluations[#]["valid"]] &&
      routeEvaluations[#]["result_status"] === "REGISTERED_DEFERRED" &
  ];
  independentSet = Complement[
    activeParams,
    Union[
      dependencySet,
      anchorSet,
      calibratedGeometrySet,
      earnedSet,
      solvedSet,
      deferredSet
    ]
  ];
  dependentOnIrreducibleSet = paramsOf[
    Select[
      dependencyRows,
      #["param"] === "B_eff" &&
        ! DisjointQ[#["dependencies"], independentSet] &
    ]
  ];
  ordinaryDependentSet =
    Complement[dependencySet, dependentOnIrreducibleSet];

  partitions = {
    {"OUT_OF_ACTIVE_NG5", outSet},
    {"BASE_SUBSTRATE", baseSet},
    {"ACCEPTED_GEOMETRY_SUBSTRATE", geometryBaseSet},
    {"ACCEPTED_PROFILE_GIVEN_ell_g", profileSet},
    {"DEPENDENT_ON_IRREDUCIBLE", dependentOnIrreducibleSet},
    {"DEPENDENT", ordinaryDependentSet},
    {"CALIBRATED_ANCHOR", anchorSet},
    {"CALIBRATED_GEOMETRY_INPUT", calibratedGeometrySet},
    {"REDUCIBLE_DERIVED", Union[earnedSet, solvedSet]},
    {"REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED", deferredSet},
    {"IRREDUCIBLY_INDEPENDENT", independentSet}
  };
  originAssociation = Association[
    (
      # -> FirstCase[
        partitions,
        {category_, members_} /; MemberQ[members, #] :> category
      ]
    ) & /@ allParams
  ];
  groupedRows = GroupBy[rows, originAssociation[#["param"]] &];
  activeIrreducible = Select[
    activeParams,
    originAssociation[#] === "IRREDUCIBLY_INDEPENDENT" &
  ];
  <|
    "rows" -> rows,
    "activeRows" -> activeRows,
    "routeEvaluations" -> routeEvaluations,
    "origin" -> originAssociation,
    "originGroups" -> groupedRows,
    "activeIrreducible" -> activeIrreducible
  |>
];

runCase[
  rows_: Automatic,
  routes_: Automatic,
  anchors_: {"Q_E"},
  geometry_: {"ell", "b", "M_h"}
] := classifyByCategorySets[
  If[rows === Automatic, buildRows[], rows],
  If[routes === Automatic, baseRoutes[], routes],
  anchors,
  geometry
];

simDeferredMap[classified_Association] := Association[
  (
    # -> classified["routeEvaluations"][#]["route_name"]
  ) & /@ Select[
    Lookup[classified["activeRows"], "param"],
    classified["origin"][#] ===
      "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED" &
  ]
];

calibratedMap[classified_Association] := Association[
  (
    # -> classified["origin"][#]
  ) & /@ Select[
    Lookup[classified["rows"], "param"],
    MemberQ[
      {"CALIBRATED_ANCHOR", "CALIBRATED_GEOMETRY_INPUT"},
      classified["origin"][#]
    ] &
  ]
];

decideVerdict[classified_Association, noGo_: False] := Module[
  {
    irreducible,
    simDeferred,
    calibrated,
    verdict,
    driftCount,
    reportedIrreducible,
    simText,
    calibratedText
  },
  irreducible = classified["activeIrreducible"];
  simDeferred = simDeferredMap[classified];
  calibrated = calibratedMap[classified];
  If[
    TrueQ[noGo],
    verdict = "NO_GO(cone-lock-feedback)";
    driftCount = 0;
    reportedIrreducible = {},
    If[
      irreducible =!= {},
      verdict =
        "SECOND_MEDIUM_DRIFT(active_irreducible={" <>
        StringRiffle[irreducible, ","] <> "})";
      driftCount = Length[irreducible];
      reportedIrreducible = irreducible,
      If[
        Length[simDeferred] > 0,
        simText = StringRiffle[
          KeyValueMap[#1 <> "->" <> #2 &, simDeferred],
          ","
        ];
        calibratedText = StringRiffle[
          KeyValueMap[#1 <> "->" <> #2 &, calibrated],
          ","
        ];
        verdict =
          "ONE_MEDIUM_CONDITIONAL(sim-deferred={" <>
          simText <> "}; calibrated={" <> calibratedText <> "})";
        driftCount = 0;
        reportedIrreducible = {},
        verdict = "ONE_MEDIUM_CONSISTENT";
        driftCount = 0;
        reportedIrreducible = {}
      ]
    ]
  ];
  <|
    "verdict" -> verdict,
    "driftCount" -> driftCount,
    "activeIrreducible" -> reportedIrreducible,
    "simDeferred" -> simDeferred,
    "calibrated" -> calibrated
  |>
];

freedomState[classified_Association, param_String] := Module[
  {evaluation = classified["routeEvaluations"][param], table},
  table = {
    {
      TrueQ[evaluation["valid"]] &&
        evaluation["result_status"] === "REGISTERED_DEFERRED",
      "FREEDOM_SIM_DEFERRED{" <> ToString[evaluation["route_name"]] <> "}"
    },
    {
      TrueQ[evaluation["valid"]] &&
        evaluation["result_status"] === "SOLVED_PASS",
      "FREEDOM_REDUCED{" <> ToString[evaluation["route_name"]] <> "}"
    },
    {True, "FREEDOM_CERTIFIED_CURRENT_LEDGER{" <> param <> "}"}
  };
  FirstCase[table, {True, value_} :> value]
];

locationClosure[rows_List] := Module[{offending},
  offending = Select[rows, ! MemberQ[arenas, #["location"]] &];
  <|
    "noFourthArena" -> (offending === {}),
    "offending" -> If[
      offending === {},
      {},
      Lookup[offending, {"param", "location"}]
    ]
  |>
];

arenaLabels[rows_List] := Association[
  (#["param"] -> #["location"] &) /@
    Select[rows, MemberQ[expectedTrio, #["param"]] &]
];


(* ---------------------------------------------------------------------- *)
(* Independent lineage, freeze projection, dimension, and algebra objects. *)
(* ---------------------------------------------------------------------- *)

productionLineageFacts[] := <|
  "varrhoObject" -> "pathA_25_closed_density_smectic_varrho_br",
  "rhoObject" -> "pathA_35_active_shear_surface_rho_br",
  "varrhoActive" -> False,
  "rhoActive" -> True,
  "varrhoRole" -> "density_smectic_layer_kinetic",
  "rhoRole" -> "shear_surface_kinetic",
  "varrhoDimension" -> {1, -3, 0},
  "rhoDimension" -> {1, -3, 0},
  "residualMultiplier" -> None,
  "routeAPending" -> True
|>;

residualLineageFacts[multiplier_String] := <|
  "varrhoObject" -> "same_active_shear_surface_object",
  "rhoObject" -> "same_active_shear_surface_object",
  "varrhoActive" -> True,
  "rhoActive" -> True,
  "varrhoRole" -> "shear_surface_kinetic",
  "rhoRole" -> "shear_surface_kinetic",
  "varrhoDimension" -> {1, -3, 0},
  "rhoDimension" -> {1, -3, 0},
  "residualMultiplier" -> multiplier,
  "routeAPending" -> True
|>;

lineageAdjudication[facts_Association] := Module[
  {
    sameActiveObject,
    sameRole,
    sameDimension,
    residualClosed,
    same,
    orderedFindings
  },
  sameActiveObject =
    TrueQ[facts["varrhoActive"]] &&
    TrueQ[facts["rhoActive"]] &&
    facts["varrhoObject"] === facts["rhoObject"];
  sameRole = facts["varrhoRole"] === facts["rhoRole"];
  sameDimension =
    facts["varrhoDimension"] === facts["rhoDimension"];
  residualClosed =
    MemberQ[{None, "1"}, facts["residualMultiplier"]];
  same =
    sameActiveObject && sameRole && sameDimension && residualClosed;
  orderedFindings = {
    {same, {"SAME", "OVERCOUNT_OR_SMUGGLE_CONTROL"}},
    {
      sameActiveObject && ! residualClosed,
      {
        "DIFFERENT",
        "DIFFERENT_BRANE_INERTIA_OBJECTS_RESIDUAL_MULTIPLIER"
      }
    },
    {
      ! TrueQ[facts["varrhoActive"]] &&
        TrueQ[facts["rhoActive"]] &&
        facts["varrhoObject"] =!= facts["rhoObject"] &&
        TrueQ[facts["routeAPending"]],
      {"DIFFERENT", "NO_OVERCOUNT_ROUTE_A_PENDING"}
    },
    {True, {"DIFFERENT", "DIFFERENT_BRANE_INERTIA_OBJECTS"}}
  };
  FirstCase[orderedFindings, {True, value_} :> value]
];

lineageFromLedger[
  classified_Association,
  facts_Association
] := lineageAdjudication[
  Join[
    facts,
    <|
      "routeAPending" ->
        (
          classified["origin"]["rho_br"] ===
            "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED"
        )
    |>
  ]
];

reductionStatus[
  classified_Association,
  excluded_List
] := Association[
  KeyValueMap[
    If[
      MemberQ[excluded, #1],
      Nothing,
      #1 -> "REGISTERED_PENDING(" <> #2 <> ")"
    ] &,
    simDeferredMap[classified]
  ]
];

historicalStructuralPostulates = {
  {
    "postulate_1",
    "imposed w-hat axis + w=0 surface (conceded-wall)"
  },
  {
    "postulate_2",
    "u^a same-medium surface collective, tangentially free-slip"
  },
  {
    "postulate_3",
    "T0 P^i reused as the Cosserat micro-rotation reservoir"
  },
  {"postulate_4", "baseline P^i spin-wave status = massless"},
  {
    "postulate_5",
    "w-hat-dependent parity-even P-u operator structural cost"
  },
  {
    "postulate_6",
    "no C5 phi analog / no longitudinal constraint"
  }
};

retiredDriftKeys = {
  "lambda_Pu", "postulate_3", "postulate_4", "postulate_5"
};
expectedOperativeKeys = {
  "rho_br",
  "mu_R",
  "Omega_w",
  "g_ell",
  "postulate_1",
  "postulate_2",
  "postulate_6"
};

freezeHistory[renameSurvivor_: False] := Module[{history},
  history = Join[
    {
      <|"key" -> "rho_br", "name" -> "rho_br"|>,
      <|"key" -> "mu_R", "name" -> "mu_R"|>,
      <|"key" -> "lambda_Pu", "name" -> "lambda_Pu"|>,
      <|"key" -> "Omega_w", "name" -> "Omega_w"|>,
      <|"key" -> "g_ell", "name" -> "g_ell(w)"|>
    },
    (<|"key" -> #[[1]], "name" -> #[[2]]|> &) /@
      historicalStructuralPostulates
  ];
  If[
    TrueQ[renameSurvivor],
    history = (
      If[
        #["key"] === "rho_br",
        Join[#, <|"name" -> "rho_B0"|>],
        #
      ] &
    ) /@ history
  ];
  history
];

antiAbsorptionState[renameSurvivor_: False] := Module[
  {
    historical,
    historicalByKey,
    operativeKeys,
    operative,
    operativeNames
  },
  historical = freezeHistory[renameSurvivor];
  historicalByKey = Association[
    (#["key"] -> # &) /@ historical
  ];
  operativeKeys = Complement[Keys[historicalByKey], retiredDriftKeys];
  operative = Select[
    historical,
    MemberQ[operativeKeys, #["key"]] &
  ];
  operativeNames = DeleteDuplicates[Lookup[operative, "name"]];
  <|
    "operativeKeys" -> operativeKeys,
    "operativeCount" -> Length[operative],
    "operativeNames" -> operativeNames,
    "disjoint" -> DisjointQ[expectedTrio, operativeNames]
  |>
];

currentNonentailmentWitness[] := Module[
  {values, residuals, positive, lockB},
  values = <|
    "rho" -> 1,
    "rho_br" -> 1,
    "rho_B0" -> 2,
    "chi_c" -> 4,
    "mu_R" -> 2,
    "K" -> 1,
    "m" -> 5,
    "M_h" -> 1,
    "c_E" -> 3,
    "C_hu" -> 1/2,
    "K_h" -> 9,
    "B_eff" -> 1,
    "sigma" -> 35/4
  |>;
  residuals = {
    values["B_eff"] values["chi_c"] - values["rho_B0"]^2,
    values["K_h"] - values["M_h"] values["c_E"]^2,
    values["B_eff"] values["K_h"] -
      values["C_hu"]^2 - values["sigma"]
  };
  positive = And @@ Thread[Values[values] > 0];
  lockB = values["c_E"]^2 values["rho_br"] - values["mu_R"];
  <|
    "valid" -> (And @@ Thread[Simplify[residuals] == 0] && positive),
    "lockB" -> lockB
  |>
];

freedomTieStatus[includeLocks_: True] := Module[
  {B, Kh, sigma, residual, expectedNegative, witness, residuals, positive},
  If[
    TrueQ[includeLocks],
    residual = Factor[B Kh - 2 B Kh - sigma];
    expectedNegative = Factor[-(B Kh + sigma)];
    Return[
      <|
        "status" ->
          If[residual === expectedNegative, "UNSAT", "UNRESOLVED"],
        "certificateResidual" -> residual
      |>
    ]
  ];
  witness = <|
    "rho_br" -> 1,
    "mu_R" -> 1,
    "M_h" -> 1,
    "c_E" -> 2,
    "K_h" -> 4,
    "B_eff" -> 1,
    "C_hu" -> Sqrt[2],
    "sigma" -> 2,
    "q_h_sq" -> 2
  |>;
  residuals = {
    witness["K_h"] - witness["M_h"] witness["c_E"]^2,
    witness["B_eff"] witness["K_h"] -
      witness["C_hu"]^2 - witness["sigma"],
    witness["C_hu"]^2 - witness["q_h_sq"],
    witness["q_h_sq"] witness["rho_br"] -
      2 witness["B_eff"] witness["M_h"] witness["mu_R"]
  };
  positive = And @@ Thread[Values[witness] > 0];
  <|
    "status" ->
      If[
        And @@ Thread[Simplify[residuals] == 0] && positive,
        "SAT",
        "UNRESOLVED"
      ],
    "witnessResiduals" -> residuals
  |>
];

dimensionState[corrupt_: False] := Module[
  {
    dimMass,
    dimBulkDensity,
    dimLayer,
    dimRhoBr,
    dimSpeed,
    varrho,
    coneProduct
  },
  dimMass = {1, 0, 0};
  dimBulkDensity = {0, -4, 0};
  dimLayer = {0, 1, 0};
  dimRhoBr = {1, -3, 0};
  dimSpeed = {0, If[TrueQ[corrupt], 2, 1], -1};
  varrho = dimMass + dimBulkDensity + dimLayer;
  coneProduct = 2 dimSpeed + dimRhoBr;
  {varrho, dimRhoBr, coneProduct}
];


(* ---------------------------------------------------------------------- *)
(* Exact production ledgers and control helpers.                           *)
(* ---------------------------------------------------------------------- *)

expectedFullOrigin = <|
  "rho" -> "BASE_SUBSTRATE",
  "K" -> "BASE_SUBSTRATE",
  "m" -> "BASE_SUBSTRATE",
  "a" -> "BASE_SUBSTRATE",
  "ell_g" -> "ACCEPTED_GEOMETRY_SUBSTRATE",
  "g_ell(w)" -> "ACCEPTED_PROFILE_GIVEN_ell_g",
  "varrho_br[rho]" -> "OUT_OF_ACTIVE_NG5",
  "Sigma_n[rho]" -> "OUT_OF_ACTIVE_NG5",
  "delta_Sigma[rho]" -> "OUT_OF_ACTIVE_NG5",
  "rho_br" -> "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED",
  "mu_R" -> "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED",
  "c_E" -> "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED",
  "c_gamma" -> "DEPENDENT",
  "rho_B0" -> "IRREDUCIBLY_INDEPENDENT",
  "chi_c" -> "IRREDUCIBLY_INDEPENDENT",
  "B_eff" -> "DEPENDENT_ON_IRREDUCIBLE",
  "C_hu" -> "IRREDUCIBLY_INDEPENDENT",
  "Q_E" -> "CALIBRATED_ANCHOR",
  "ell" -> "CALIBRATED_GEOMETRY_INPUT",
  "b" -> "CALIBRATED_GEOMETRY_INPUT",
  "M_h" -> "CALIBRATED_GEOMETRY_INPUT",
  "K_h" -> "DEPENDENT",
  "q_h" -> "DEPENDENT",
  "c_L1" -> "OUT_OF_ACTIVE_NG5",
  "c_L2" -> "OUT_OF_ACTIVE_NG5",
  "A_R" -> "OUT_OF_ACTIVE_NG5",
  "k_R" -> "OUT_OF_ACTIVE_NG5",
  "lambda_Cdiv" -> "OUT_OF_ACTIVE_NG5",
  "chi_Cpin" -> "OUT_OF_ACTIVE_NG5",
  "J_Pu" -> "OUT_OF_ACTIVE_NG5",
  "kappa_Pu" -> "OUT_OF_ACTIVE_NG5",
  "lambda_Pu" -> "OUT_OF_ACTIVE_NG5",
  "Omega_w" -> "OUT_OF_ACTIVE_NG5",
  "lambda_N" -> "OUT_OF_ACTIVE_NG5",
  "lambda_tau" -> "OUT_OF_ACTIVE_NG5",
  "Nu" -> "OUT_OF_ACTIVE_NG5",
  "a_T" -> "OUT_OF_ACTIVE_NG5",
  "a_Tp" -> "OUT_OF_ACTIVE_NG5",
  "a_L" -> "OUT_OF_ACTIVE_NG5"
|>;

expectedRouteLedger = <|
  "rho_br" -> <|
    "route_name" -> "Route-A",
    "valid" -> True,
    "reason" -> "VALID_REGISTERED_ROUTE"
  |>,
  "mu_R" -> <|
    "route_name" -> "Route-A",
    "valid" -> True,
    "reason" -> "VALID_REGISTERED_ROUTE"
  |>,
  "c_E" -> <|
    "route_name" -> "Route-A",
    "valid" -> True,
    "reason" -> "VALID_REGISTERED_ROUTE"
  |>,
  "c_gamma" -> <|
    "route_name" -> None,
    "valid" -> False,
    "reason" -> "NO_REGISTERED_ROUTE_FOR_PARAM"
  |>,
  "rho_B0" -> <|
    "route_name" -> "Future-Compression-4D-to-3D",
    "valid" -> False,
    "reason" -> "REJECTED_RESULT_STATUS_PROMISSORY_ONLY"
  |>,
  "chi_c" -> <|
    "route_name" -> "Future-Compression-4D-to-3D",
    "valid" -> False,
    "reason" -> "REJECTED_RESULT_STATUS_PROMISSORY_ONLY"
  |>,
  "B_eff" -> <|
    "route_name" -> None,
    "valid" -> False,
    "reason" -> "NO_REGISTERED_ROUTE_FOR_PARAM"
  |>,
  "C_hu" -> <|
    "route_name" -> "Future-Embedding-Overlap",
    "valid" -> False,
    "reason" -> "REJECTED_RESULT_STATUS_PROMISSORY_ONLY"
  |>,
  "Q_E" -> <|
    "route_name" -> None,
    "valid" -> False,
    "reason" -> "NO_REGISTERED_ROUTE_FOR_PARAM"
  |>,
  "ell" -> <|
    "route_name" -> None,
    "valid" -> False,
    "reason" -> "NO_REGISTERED_ROUTE_FOR_PARAM"
  |>,
  "b" -> <|
    "route_name" -> None,
    "valid" -> False,
    "reason" -> "NO_REGISTERED_ROUTE_FOR_PARAM"
  |>,
  "M_h" -> <|
    "route_name" -> None,
    "valid" -> False,
    "reason" -> "NO_REGISTERED_ROUTE_FOR_PARAM"
  |>,
  "K_h" -> <|
    "route_name" -> None,
    "valid" -> False,
    "reason" -> "NO_REGISTERED_ROUTE_FOR_PARAM"
  |>,
  "q_h" -> <|
    "route_name" -> None,
    "valid" -> False,
    "reason" -> "NO_REGISTERED_ROUTE_FOR_PARAM"
  |>
|>;

routeLedger[classified_Association] := Association[
  (
    Function[param,
      param -> KeyTake[
        classified["routeEvaluations"][param],
        {"route_name", "valid", "reason"}
      ]
    ][#]
  ) & /@ Lookup[classified["activeRows"], "param"]
];

recordedState[classified_Association, dropParam_: None] := Module[
  {active, recorded},
  active = Lookup[classified["activeRows"], "param"];
  recorded = Select[
    active,
    # =!= dropParam &&
      AssociationQ[Lookup[classified["routeEvaluations"], #, None]] &
  ];
  {active, recorded}
];

routeConjunctControl[conjunct_String, neutralize_: False] := Module[
  {
    routes,
    transform,
    classified,
    expectedFlips,
    actualFlips,
    isolated,
    checks,
    falseNamed,
    expectedFalse
  },
  routes = baseRoutes[];
  If[
    ! TrueQ[neutralize],
    transform = FirstCase[
      {
        {
          "required_fields_present",
          Function[fact, KeyDrop[fact, {"source"}]]
        },
        {
          "named_solve_in_provenance",
          Function[fact, Join[
            fact,
            <|"named_solve_in_provenance" -> False|>
          ]]
        },
        {
          "result_status_allowed",
          Function[fact, Join[
            fact,
            <|"result_status" -> "BY_TUNING"|>
          ]]
        },
        {
          "finite_listed_missing_objects",
          Function[fact, Join[
            fact,
            <|"missing_objects" -> "not-a-list"|>
          ]]
        },
        {
          "target_identity",
          Function[fact, Join[
            fact,
            <|"targets" -> {"rho_br", "mu_R"}|>
          ]]
        },
        {
          "joint_target_identity",
          Function[fact, Join[
            fact,
            <|
              "required_joint_targets" ->
                {"rho_br", "mu_R", "joint_sentinel"}
            |>
          ]]
        },
        {
          "target_blind",
          Function[fact, Join[fact, <|"target_blind" -> False|>]]
        },
        {
          "falsifiers_clear",
          Function[fact, Join[
            fact,
            <|"falsifiers" -> {"synthetic-falsifier"}|>
          ]]
        }
      },
      {conjunct, value_} :> value
    ];
    routes = updateRoute[routes, "Route-A", transform]
  ];
  classified = runCase[Automatic, routes];
  expectedFlips =
    If[conjunct === "target_identity", {"c_E"}, routeATargets];
  actualFlips = Select[
    routeATargets,
    classified["origin"][#] === "IRREDUCIBLY_INDEPENDENT" &
  ];
  isolated = And @@ (
    Function[param,
      checks = classified["routeEvaluations"][param]["checks"];
      falseNamed = Select[
        routeConjuncts,
        ! TrueQ[checks[#]] &
      ];
      expectedFalse =
        If[MemberQ[expectedFlips, param], {conjunct}, {}];
      Sort[falseNamed] === Sort[expectedFalse]
    ] /@ routeATargets
  );
  <|
    "fired" ->
      (Sort[actualFlips] === Sort[expectedFlips] && isolated),
    "actualFlips" -> actualFlips,
    "isolated" -> isolated
  |>
];


(* ---------------------------------------------------------------------- *)
(* Bounded source-to-stage manifest.                                       *)
(* ---------------------------------------------------------------------- *)

sourcePredicateTotal = 37;
sourcePredicateUniverse = {
  "pathA41.py::dim_derivations",
  "pathA41.py::registered_route_for",
  "pathA41.py::classify_origin",
  "pathA41.py::classify_rows",
  "pathA41.py::decide_verdict",
  "pathA41.py::location_closure",
  "pathA41.py::lineage_adjudication",
  "pathA41.py::freedom_states",
  "pathA41.py::interpretation.reduction_status",
  "pathA41.py::interpretation.physical_meaning",
  "pathA41.py::production_location_closure_raise",
  "pathA41.py::p40_freedom_tie_no_go",
  "pathA41.py::p40_current_nonentailment",
  "pathA41.py::control.AB_delete_registry",
  "pathA41.py::control.route_blank_Route_A",
  "pathA41.py::control.route_field_blank_Route_A_target_blind",
  "pathA41.py::control.route_field_blank_Route_A_missing_objects",
  "pathA41.py::control.calibration_ablation_Q_E",
  "pathA41.py::control.irreducible_synthetic",
  "pathA41.py::control.reducible_derived_synthetic",
  "pathA41.py::control.location_closure_out_of_arena",
  "pathA41.py::control.contradiction",
  "pathA41.py::control.residual_multiplier_ablation",
  "pathA41.py::control.route_eval_recorded_for_all_active_rows",
  "pathA41.wl::hardcoded_production",
  "pathA41.wl::literal_controls",
  "pathA41.py::argparse_compare_harness",
  "pathA41.py::file_writing_and_main",
  "pathA41.py::build_results_yaml_write_report",
  "pathA41.py::comparison_payload_digest_sha_count",
  "pathA41.py::compare_with_mathematica",
  "pathA41.py::filesystem_source_token_scans",
  "pathA41.wl::filesystem_assertContains_scans",
  "pathA41.py::importlib_pathA40",
  "pathA41.py::extraction_audit_narration",
  "pathA41.py::comparison_map_json_serialization",
  "pathA41.wl::Import_canon_cross_read"
};

preservedSourceIDs = Join[
  Take[sourcePredicateUniverse, 7],
  {
    sourcePredicateUniverse[[11]],
    sourcePredicateUniverse[[12]],
    sourcePredicateUniverse[[13]]
  },
  Take[sourcePredicateUniverse, {14, 24}]
];
replacedSourceIDs = sourcePredicateUniverse[[{8, 9, 10, 25, 26}]];
scopedSourceIDs = Take[sourcePredicateUniverse, {27, 37}];
newStageIDs = Join[
  {
    "stage041::ANTI_ABSORPTION_DRIFT_SIDE",
    "stage041::DUAL_ENGINE_TERMS",
    "stage041::SOURCE_TO_STAGE_MANIFEST"
  },
  ("stage041::" <> # &) /@ conjunctTeeth
];

manifestRows[] := Join[
  ({#, "preserved-folded"} &) /@ preservedSourceIDs,
  ({#, "replaced-by-stronger"} &) /@ replacedSourceIDs,
  ({#, "scoped-out"} &) /@ scopedSourceIDs,
  ({#, "newly-added"} &) /@ newStageIDs
];

expectedManifestCategory = Association[
  (#[[1]] -> #[[2]] &) /@ manifestRows[]
];

manifestGuard[mutate_: False] := Module[
  {
    rows,
    dropped,
    moved,
    identifiers,
    sourceRows,
    sourceIdentifiers,
    stageIdentifiers,
    categories,
    disjoint,
    prefixes
  },
  rows = manifestRows[];
  If[
    TrueQ[mutate],
    dropped = First[scopedSourceIDs];
    rows = Select[rows, #[[1]] =!= dropped &];
    moved = First[preservedSourceIDs];
    rows = (
      If[
        #[[1]] === moved,
        {#[[1]], "replaced-by-stronger"},
        #
      ] &
    ) /@ rows
  ];
  identifiers = rows[[All, 1]];
  sourceRows = Select[rows, MemberQ[sourcePredicateUniverse, #[[1]]] &];
  sourceIdentifiers = sourceRows[[All, 1]];
  stageIdentifiers = Select[identifiers, MemberQ[newStageIDs, #] &];
  categories = Association[
    (
      # -> Cases[rows, {identifier_, #} :> identifier]
    ) & /@ {
      "preserved-folded",
      "replaced-by-stronger",
      "newly-added",
      "scoped-out"
    }
  ];
  disjoint =
    DuplicateFreeQ[Flatten[Values[categories]]];
  prefixes = First[StringSplit[#, "::"]] & /@ identifiers;
  <|
    "ok" -> (
      Length[sourceRows] === sourcePredicateTotal &&
      Sort[sourceIdentifiers] === Sort[sourcePredicateUniverse] &&
      DuplicateFreeQ[sourceIdentifiers] &&
      Sort[stageIdentifiers] === Sort[newStageIDs] &&
      disjoint &&
      Association[(#[[1]] -> #[[2]] &) /@ rows] ===
        expectedManifestCategory &&
      And @@ (MemberQ[{"pathA41.py", "pathA41.wl", "stage041"}, #] & /@
        prefixes)
    ),
    "sourceCount" -> Length[sourceRows],
    "partition" -> Counts[rows[[All, 2]]],
    "disjoint" -> disjoint
  |>
];


(* ---------------------------------------------------------------------- *)
(* Executable assertions.                                                  *)
(* ---------------------------------------------------------------------- *)

runAssertions[] := Module[
  {
    productionRows,
    production,
    productionResult,
    productionIrreducible,
    originCase,
    originRoutes,
    fullOriginCase,
    routeCase,
    verdictCase,
    verdictRoutes,
    verdictActual,
    driftCase,
    driftCount,
    rhoFreedomCase,
    rhoFreedom,
    cHuFreedomCase,
    cHuFreedomRoutes,
    cHuFreedom,
    closureRows,
    closure,
    labelRows,
    labels,
    expectedLabels,
    lineageFacts,
    lineage,
    qEIsolation,
    pending,
    simKeys,
    anti,
    blankAll,
    blankAllIrreducible,
    noRouteARoutes,
    noRouteA,
    fieldTooth,
    fieldName,
    fieldRoutes,
    fieldCase,
    fieldFired,
    conjunctMap,
    conjunctControl,
    qECase,
    qEIrreducible,
    xiRows,
    xiCase,
    pSynRows,
    pSynCase,
    witness,
    tie,
    contradictionVerdict,
    recordState,
    closedLineage,
    freeLineage,
    controlRows,
    productionClosure,
    controlClosure,
    locationControlFired,
    dimensions,
    expectedDimensions,
    antiProduction,
    productionLineage,
    localInventory,
    expectedInventory,
    blankVerdict,
    qEVerdict,
    noGoVerdict,
    conditionalRoutes,
    conditionalVerdict,
    consistentVerdict,
    actualVerdicts,
    expectedConditional,
    expectedVerdicts,
    manifest,
    arenaGroups,
    physicalMeaning
  },
  productionRows = buildRows[];
  production = runCase[productionRows];
  productionResult = decideVerdict[production];
  productionIrreducible = production["activeIrreducible"];

  section["Computed production origin and RouteEvaluation ledgers"];
  originCase = production;
  If[
    activeMutation === "ORIGIN_CLASSIFICATION",
    originRoutes = updateRoute[
      baseRoutes[],
      "Future-Embedding-Overlap",
      Function[fact, Join[
        fact,
        <|
          "named_solve_in_provenance" -> True,
          "result_status" -> "REGISTERED_DEFERRED"
        |>
      ]]
    ];
    originCase = runCase[Automatic, originRoutes]
  ];
  expectBool[
    "ORIGIN_CLASSIFICATION",
    originCase["activeIrreducible"] === expectedTrio,
    <|
      "Computed" -> originCase["activeIrreducible"],
      "Ratified" -> expectedTrio
    |>
  ];

  fullOriginCase =
    If[
      activeMutation === "FULL_ORIGIN_LEDGER",
      runCase[buildRows[<|"bEffDependencyMutation" -> True|>]],
      production
    ];
  expectBool[
    "FULL_ORIGIN_LEDGER",
    fullOriginCase["origin"] === expectedFullOrigin,
    fullOriginCase["origin"]
  ];

  routeCase =
    If[
      activeMutation === "FULL_ROUTE_EVALUATION_LEDGER",
      runCase[buildRows[<|"cGammaRouteMutation" -> True|>]],
      production
    ];
  expectBool[
    "FULL_ROUTE_EVALUATION_LEDGER",
    routeLedger[routeCase] === expectedRouteLedger,
    routeLedger[routeCase]
  ];

  section["Production verdict, count, freedom, locations, and lineage"];
  verdictCase = production;
  If[
    activeMutation === "PRODUCTION_VERDICT",
    verdictRoutes = updateRoute[
      baseRoutes[],
      "Future-Embedding-Overlap",
      Function[fact, Join[
        fact,
        <|
          "named_solve_in_provenance" -> True,
          "result_status" -> "REGISTERED_DEFERRED"
        |>
      ]]
    ];
    verdictCase = runCase[Automatic, verdictRoutes]
  ];
  verdictActual = decideVerdict[verdictCase]["verdict"];
  expectBool[
    "PRODUCTION_VERDICT",
    verdictActual === expectedVerdict,
    verdictActual
  ];

  driftCase =
    If[
      activeMutation === "DRIFT_COUNT",
      runCase[buildRows[<|"addXi" -> True|>]],
      production
    ];
  driftCount = Length[driftCase["activeIrreducible"]];
  expectBool["DRIFT_COUNT", driftCount === 3, driftCount];

  rhoFreedomCase =
    If[
      activeMutation === "FREEDOM_STATE_RHO_BR",
      runCase[Automatic, {}],
      production
    ];
  rhoFreedom = freedomState[rhoFreedomCase, "rho_br"];
  expectBool[
    "FREEDOM_STATE_RHO_BR",
    rhoFreedom === "FREEDOM_SIM_DEFERRED{Route-A}",
    rhoFreedom
  ];

  cHuFreedomCase = production;
  If[
    activeMutation === "FREEDOM_STATE_C_HU",
    cHuFreedomRoutes = updateRoute[
      baseRoutes[],
      "Future-Embedding-Overlap",
      Function[fact, Join[
        fact,
        <|
          "named_solve_in_provenance" -> True,
          "result_status" -> "REGISTERED_DEFERRED"
        |>
      ]]
    ];
    cHuFreedomCase = runCase[Automatic, cHuFreedomRoutes]
  ];
  cHuFreedom = freedomState[cHuFreedomCase, "C_hu"];
  expectBool[
    "FREEDOM_STATE_C_HU",
    cHuFreedom === "FREEDOM_CERTIFIED_CURRENT_LEDGER{C_hu}",
    cHuFreedom
  ];

  closureRows =
    If[
      activeMutation === "LOCATION_CLOSURE",
      buildRows[<|"addLocationSentinel" -> True|>],
      productionRows
    ];
  closure = locationClosure[closureRows];
  expectBool[
    "LOCATION_CLOSURE",
    TrueQ[closure["noFourthArena"]],
    closure
  ];

  labelRows =
    If[
      activeMutation === "ARENA_LABELS",
      buildRows[<|"cHuArenaMutation" -> True|>],
      productionRows
    ];
  labels = arenaLabels[labelRows];
  expectedLabels = <|
    "rho_B0" -> "3D brane surface",
    "chi_c" -> "3D brane surface",
    "C_hu" -> "throat/embedding seam"
  |>;
  expectBool["ARENA_LABELS", labels === expectedLabels, labels];

  lineageFacts =
    If[
      activeMutation === "LINEAGE_FINDING",
      residualLineageFacts["Xi_residual"],
      productionLineageFacts[]
    ];
  lineage = lineageFromLedger[production, lineageFacts];
  qEIsolation = lineageFromLedger[
    runCase[Automatic, Automatic, {}],
    productionLineageFacts[]
  ];
  expectBool[
    "LINEAGE_FINDING",
    lineage === {"DIFFERENT", "NO_OVERCOUNT_ROUTE_A_PENDING"} &&
      qEIsolation === {"DIFFERENT", "NO_OVERCOUNT_ROUTE_A_PENDING"},
    <|"Production" -> lineage, "QEIsolation" -> qEIsolation|>
  ];

  pending = reductionStatus[
    production,
    If[activeMutation === "REDUCTION_STATUS", {"c_E"}, {}]
  ];
  simKeys = Keys[simDeferredMap[production]];
  expectBool[
    "REDUCTION_STATUS",
    Sort[Keys[pending]] === Sort[routeATargets] &&
      Sort[simKeys] === Sort[routeATargets],
    <|"Pending" -> pending, "SimKeys" -> simKeys|>
  ];

  anti = antiAbsorptionState[
    activeMutation === "ANTI_ABSORPTION_DRIFT_SIDE"
  ];
  expectBool[
    "ANTI_ABSORPTION_DRIFT_SIDE",
    Sort[anti["operativeKeys"]] === Sort[expectedOperativeKeys] &&
      anti["operativeCount"] === 7 &&
      TrueQ[anti["disjoint"]],
    anti
  ];

  section["Individually asserted controls"];
  blankAll =
    If[
      activeMutation === "CTRL_AB_DELETE_REGISTRY",
      production,
      runCase[Automatic, {}]
    ];
  blankAllIrreducible = blankAll["activeIrreducible"];
  expectBool[
    "CTRL_AB_DELETE_REGISTRY",
    SubsetQ[blankAllIrreducible, productionIrreducible] &&
      Length[blankAllIrreducible] > Length[productionIrreducible] &&
      blankAllIrreducible === {
        "rho_br", "mu_R", "c_E", "rho_B0", "chi_c", "C_hu"
      },
    blankAllIrreducible
  ];

  noRouteARoutes = baseRoutes[];
  If[
    activeMutation =!= "CTRL_ROUTE_BLANK_ROUTE_A",
    noRouteARoutes =
      Select[noRouteARoutes, #["name"] =!= "Route-A" &]
  ];
  noRouteA = runCase[Automatic, noRouteARoutes];
  expectBool[
    "CTRL_ROUTE_BLANK_ROUTE_A",
    And @@ (
      noRouteA["origin"][#] === "IRREDUCIBLY_INDEPENDENT" &
    ) /@ routeATargets,
    Association[(# -> noRouteA["origin"][#] &) /@ routeATargets]
  ];

  Do[
    fieldTooth = pair[[1]];
    fieldName = pair[[2]];
    fieldRoutes = baseRoutes[];
    If[
      activeMutation =!= fieldTooth,
      fieldRoutes = updateRoute[
        fieldRoutes,
        "Route-A",
        Function[fact, KeyDrop[fact, {fieldName}]]
      ]
    ];
    fieldCase = runCase[Automatic, fieldRoutes];
    fieldFired = And @@ (
      Function[param,
        fieldCase["origin"][param] === "IRREDUCIBLY_INDEPENDENT" &&
        fieldCase["routeEvaluations"][param]["reason"] ===
          "MISSING_REQUIRED_ROUTE_FIELD"
      ] /@ routeATargets
    );
    expectBool[
      fieldTooth,
      fieldFired,
      Association[
        (# -> fieldCase["routeEvaluations"][#] &) /@ routeATargets
      ]
    ],
    {pair, {
      {"CTRL_ROUTE_FIELD_BLANK_TARGET_BLIND", "target_blind"},
      {"CTRL_ROUTE_FIELD_BLANK_MISSING_OBJECTS", "missing_objects"}
    }}
  ];

  conjunctMap = AssociationThread[conjunctTeeth, routeConjuncts];
  Do[
    conjunctControl = routeConjunctControl[
      conjunctMap[tooth],
      activeMutation === tooth
    ];
    expectBool[tooth, TrueQ[conjunctControl["fired"]], conjunctControl],
    {tooth, conjunctTeeth}
  ];

  qECase = runCase[
    Automatic,
    Automatic,
    If[
      activeMutation === "CTRL_CALIBRATION_ABLATION_Q_E",
      {"Q_E"},
      {}
    ]
  ];
  qEIrreducible = qECase["activeIrreducible"];
  expectBool[
    "CTRL_CALIBRATION_ABLATION_Q_E",
    qECase["origin"]["Q_E"] === "IRREDUCIBLY_INDEPENDENT" &&
      qEIrreducible === {"rho_B0", "chi_c", "C_hu", "Q_E"} &&
      Length[qEIrreducible] === 4,
    <|"Origin" -> qECase["origin"]["Q_E"], "Irreducible" -> qEIrreducible|>
  ];

  xiRows = buildRows[
    <|
      "addXi" ->
        (activeMutation =!= "CTRL_IRREDUCIBLE_SYNTHETIC")
    |>
  ];
  xiCase = runCase[xiRows];
  expectBool[
    "CTRL_IRREDUCIBLE_SYNTHETIC",
    MemberQ[xiCase["activeIrreducible"], "xi_active"] &&
      xiCase["routeEvaluations"]["xi_active"]["reason"] ===
        "NO_REGISTERED_ROUTE_FOR_PARAM",
    xiCase["activeIrreducible"]
  ];

  pSynRows = buildRows[
    <|
      "addPSyn" -> True,
      "pSynStatusMutation" ->
        (activeMutation === "CTRL_REDUCIBLE_DERIVED_SYNTHETIC")
    |>
  ];
  pSynCase = runCase[pSynRows];
  expectBool[
    "CTRL_REDUCIBLE_DERIVED_SYNTHETIC",
    pSynCase["origin"]["p_syn"] === "REDUCIBLE_DERIVED" &&
      pSynCase["activeIrreducible"] === expectedTrio,
    <|
      "Origin" -> pSynCase["origin"]["p_syn"],
      "Irreducible" -> pSynCase["activeIrreducible"]
    |>
  ];

  witness = currentNonentailmentWitness[];
  tie = freedomTieStatus[
    activeMutation =!= "CTRL_CONTRADICTION"
  ];
  contradictionVerdict = decideVerdict[
    production,
    tie["status"] === "UNSAT"
  ]["verdict"];
  expectBool[
    "CTRL_CONTRADICTION",
    TrueQ[witness["valid"]] &&
      witness["lockB"] === 7 &&
      tie["status"] === "UNSAT" &&
      contradictionVerdict === "NO_GO(cone-lock-feedback)",
    <|
      "Witness" -> witness,
      "Tie" -> tie,
      "Verdict" -> contradictionVerdict
    |>
  ];

  recordState = recordedState[
    production,
    If[
      activeMutation === "CTRL_ROUTE_EVAL_RECORDED",
      "c_gamma",
      None
    ]
  ];
  expectBool[
    "CTRL_ROUTE_EVAL_RECORDED",
    Sort[recordState[[1]]] === Sort[recordState[[2]]],
    <|"Active" -> recordState[[1]], "Recorded" -> recordState[[2]]|>
  ];

  closedLineage = lineageAdjudication[residualLineageFacts["1"]];
  freeLineage = lineageAdjudication[
    residualLineageFacts[
      If[
        activeMutation === "CTRL_RESIDUAL_MULTIPLIER_ABLATION",
        "1",
        "Xi_residual"
      ]
    ]
  ];
  expectBool[
    "CTRL_RESIDUAL_MULTIPLIER_ABLATION",
    closedLineage === {"SAME", "OVERCOUNT_OR_SMUGGLE_CONTROL"} &&
      freeLineage === {
        "DIFFERENT",
        "DIFFERENT_BRANE_INERTIA_OBJECTS_RESIDUAL_MULTIPLIER"
      },
    <|"Closed" -> closedLineage, "Free" -> freeLineage|>
  ];

  controlRows = buildRows[
    <|
      "addLocationSentinel" ->
        (activeMutation =!=
          "CTRL_LOCATION_CLOSURE_OUT_OF_ARENA")
    |>
  ];
  productionClosure = locationClosure[productionRows];
  controlClosure = locationClosure[controlRows];
  locationControlFired =
    TrueQ[productionClosure["noFourthArena"]] &&
    ! TrueQ[controlClosure["noFourthArena"]] &&
    controlClosure["offending"] =!= {};
  expectBool[
    "CTRL_LOCATION_CLOSURE_OUT_OF_ARENA",
    locationControlFired,
    <|"Production" -> productionClosure, "Control" -> controlClosure|>
  ];

  section["Build-global dimensions, local inventory, and verdict table"];
  dimensions = dimensionState[
    activeMutation === "DIMENSION_HOMOGENEITY"
  ];
  expectedDimensions = {
    {1, -3, 0},
    {1, -3, 0},
    {1, -1, -2}
  };
  expectBool[
    "DIMENSION_HOMOGENEITY",
    dimensions === expectedDimensions,
    dimensions
  ];

  antiProduction = antiAbsorptionState[];
  productionLineage =
    lineageFromLedger[production, productionLineageFacts[]][[2]];
  localInventory = <|
    "activeIrreducible" -> productionIrreducible,
    "driftCount" -> productionResult["driftCount"],
    "freedomRhoBr" -> freedomState[production, "rho_br"],
    "freedomCHu" -> freedomState[production, "C_hu"],
    "noFourthArena" -> productionClosure["noFourthArena"],
    "lineage" -> productionLineage,
    "antiAbsorption" -> antiProduction["disjoint"],
    "dimensionState" -> dimensionState[],
    "lockBWitness" -> witness["lockB"]
  |>;
  If[
    activeMutation === "DUAL_ENGINE_TERMS",
    localInventory = KeyDrop[localInventory, {"lockBWitness"}]
  ];
  expectedInventory = <|
    "activeIrreducible" -> expectedTrio,
    "driftCount" -> 3,
    "freedomRhoBr" -> "FREEDOM_SIM_DEFERRED{Route-A}",
    "freedomCHu" -> "FREEDOM_CERTIFIED_CURRENT_LEDGER{C_hu}",
    "noFourthArena" -> True,
    "lineage" -> "NO_OVERCOUNT_ROUTE_A_PENDING",
    "antiAbsorption" -> True,
    "dimensionState" -> expectedDimensions,
    "lockBWitness" -> 7
  |>;
  expectBool[
    "DUAL_ENGINE_TERMS",
    localInventory === expectedInventory,
    localInventory
  ];

  blankVerdict = decideVerdict[runCase[Automatic, {}]]["verdict"];
  qEVerdict = decideVerdict[
    runCase[Automatic, Automatic, {}]
  ]["verdict"];
  noGoVerdict = decideVerdict[
    production,
    freedomTieStatus[True]["status"] === "UNSAT"
  ]["verdict"];
  conditionalRoutes =
    validFutureRoutes["SOLVED_PASS", "REGISTERED_DEFERRED"];
  If[
    activeMutation === "VERDICT_REDERIVATION",
    conditionalRoutes = updateRoute[
      conditionalRoutes,
      "Future-Embedding-Overlap",
      Function[fact, Join[
        fact,
        <|"named_solve_in_provenance" -> False|>
      ]]
    ]
  ];
  conditionalVerdict = decideVerdict[
    runCase[Automatic, conditionalRoutes]
  ]["verdict"];
  consistentVerdict = decideVerdict[
    runCase[
      Automatic,
      validFutureRoutes["SOLVED_PASS", "SOLVED_PASS"]
    ]
  ]["verdict"];
  actualVerdicts = {
    productionResult["verdict"],
    blankVerdict,
    qEVerdict,
    noGoVerdict,
    conditionalVerdict,
    consistentVerdict
  };
  expectedConditional =
    "ONE_MEDIUM_CONDITIONAL(sim-deferred={" <>
    "rho_br->Route-A,mu_R->Route-A,c_E->Route-A" <>
    "}; calibrated={" <>
    "Q_E->CALIBRATED_ANCHOR," <>
    "ell->CALIBRATED_GEOMETRY_INPUT," <>
    "b->CALIBRATED_GEOMETRY_INPUT," <>
    "M_h->CALIBRATED_GEOMETRY_INPUT})";
  expectedVerdicts = {
    expectedVerdict,
    "SECOND_MEDIUM_DRIFT(active_irreducible={rho_br,mu_R,c_E,rho_B0,chi_c,C_hu})",
    "SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu,Q_E})",
    "NO_GO(cone-lock-feedback)",
    expectedConditional,
    "ONE_MEDIUM_CONSISTENT"
  };
  expectBool[
    "VERDICT_REDERIVATION",
    actualVerdicts === expectedVerdicts,
    <|"Actual" -> actualVerdicts, "Expected" -> expectedVerdicts|>
  ];

  section["Bounded source-to-stage predicate manifest"];
  manifest = manifestGuard[
    activeMutation === "SOURCE_TO_STAGE_MANIFEST"
  ];
  expectBool[
    "SOURCE_TO_STAGE_MANIFEST",
    TrueQ[manifest["ok"]],
    manifest
  ];

  arenaGroups = GroupBy[
    Normal[labels],
    Last -> First
  ];
  physicalMeaning = StringRiffle[
    KeyValueMap[
      #1 <> ": " <> StringRiffle[#2, ","] &,
      arenaGroups
    ],
    "; "
  ];
  Print[""];
  Print["VERDICT=", productionResult["verdict"]];
  Print[
    "ACTIVE_IRREDUCIBLE={",
    StringRiffle[productionIrreducible, ","],
    "}"
  ];
  Print["DRIFT_COUNT=", productionResult["driftCount"]];
  Print[
    "FREEDOM_STATES=",
    freedomState[production, "rho_br"],
    ";",
    freedomState[production, "C_hu"]
  ];
  Print[
    "NO_FOURTH_ARENA=",
    productionClosure["noFourthArena"]
  ];
  Print["LINEAGE=", productionLineage];
  Print["POST_D16_DRIFT_TOKEN=POST_D16_DRIFT(7)"];
  Print["ANTI_ABSORPTION_DRIFT_SIDE=DISJOINT"];
  Print["PHYSICAL_MEANING_FROM_ARENA_MAP=", physicalMeaning];
  Print[
    "Q_E_NEAR_MISS=without calibration anchor, active_irreducible={rho_B0,chi_c,C_hu,Q_E}; drift_count=4"
  ];
  Print[
    "HONEST_CAVEAT=bookkeeping/labelling closure only; dynamical brane-genesis remains open and is handed to the deferred throat/embedding solve"
  ];
  Print[
    "FRAMING=CHARACTERIZED-DEPARTURE; no_fourth_arena is labelling-level"
  ];
  Print["SOURCE_PREDICATE_TOTAL=", sourcePredicateTotal];
  Print["EXECUTABLE_TOOTH_TOTAL=35"];
  productionResult["verdict"]
];

ok = Catch[
  If[
    activeMutation =!= "" && ! MemberQ[toothOrder, activeMutation],
    Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
    Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
    raise["UNKNOWN_MUTATION"]
  ];

  Print[
    "ledger_stage041_ng5_second_medium_reducibility Mathematica audit"
  ];
  Print[
    "ROUTE=GroupBy origin categories + Select set partitions + SubsetQ route targets + FirstCase resolution"
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
  "ledgerStage041Failure",
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
    "OVERALL FAIL: Mathematica stage041 audit did not close"
  ];
  Exit[1]
]
