(* Ledger stage043 Mathematica audit: irreducible count reported as a range.

   Standalone, print-only, assert-zero, and cross-engine-file-I/O-free.
   Register classifications and midway buckets are transcribed below by name.
   Association/GroupBy/Counts/Select/Tally perform the count partition.

   The independence diagnostic is materially distinct from the SymPy route:
   CAD RegionDimension[ImplicitRegion[eqs == 0 && vars > 0]] computes the
   positive-real dimensions directly.  Its honest scope is only M and Wchi.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE043_MUTATION";
activeMutation = Quiet@Check[Environment[mutationEnvironment], ""];
If[! StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];

toothOrder = {
  "BUCKET_PARTITION",
  "LOW_ENDPOINT",
  "HIGH_ENDPOINT",
  "RANGE_SPREAD_IS_C1_PLUS_C2",
  "BASE_CONTINUOUS_27_36",
  "E_CONTINUOUS_IS_13",
  "C1_ATTRIBUTION",
  "C2_ATTRIBUTION",
  "CONVENTION_OPEN_NOT_IMPOSED",
  "REDUCTION_DEBT_COUNTED_ONCE",
  "R49_OUT_OF_SCOPE",
  "DISCRETE_POSTULATE_COUNT",
  "DELTA_R_INDEPENDENCE",
  "DELTA_R_IS_DIAGNOSTIC_NOT_SUBTRACTION",
  "R35_DERIVED_IN_FORM_PLUS_RIDER",
  "KNOB_PROVENANCE_WELL_FORMED",
  "REGISTER_TO_COUNT_MANIFEST",
  "DUAL_ENGINE_TERMS",
  "IRREDUCIBLE_COUNT_RANGE",
  "COUNT_RANGE_REDERIVATION"
};

ablationDescriptions = <|
  "BUCKET_PARTITION" ->
    "double-list REG:a:hbar in the private (b) bucket",
  "LOW_ENDPOINT" ->
    "drop REG:a:hbar from the private low bucket membership",
  "HIGH_ENDPOINT" ->
    "drop REG:C2:Btilde0 from the private strict supplement",
  "RANGE_SPREAD_IS_C1_PLUS_C2" ->
    "drop one C1 member only from the private attribution count",
  "BASE_CONTINUOUS_27_36" ->
    "reclassify EOS-exponent-5 as structural-no-knob",
  "E_CONTINUOUS_IS_13" ->
    "remove rho_B0 from the private Parts-III--VI extension",
  "C1_ATTRIBUTION" ->
    "remove Mtilde from the private C1 contributor set",
  "C2_ATTRIBUTION" ->
    "count the two shorthand labels instead of six support scalars",
  "CONVENTION_OPEN_NOT_IMPOSED" ->
    "collapse the private reported range to its high endpoint",
  "REDUCTION_DEBT_COUNTED_ONCE" ->
    "exclude mu_R from the private continuous assembly",
  "R49_OUT_OF_SCOPE" ->
    "fold one R49 field profile into the private continuous assembly",
  "DISCRETE_POSTULATE_COUNT" ->
    "drop D1 from the private discrete-postulate itemization",
  "DELTA_R_INDEPENDENCE" ->
    "inject K=rho0 and mu_eta=T_w into the private baseline blocks",
  "DELTA_R_IS_DIAGNOSTIC_NOT_SUBTRACTION" ->
    "subtract the tested delta-r diagnostic from the private count",
  "R35_DERIVED_IN_FORM_PLUS_RIDER" ->
    "flip the private R35 label to PENDING-debt",
  "KNOB_PROVENANCE_WELL_FORMED" ->
    "slip derived c_s0 into the private counted set",
  "REGISTER_TO_COUNT_MANIFEST" ->
    "mis-scope hbar as DERIVED against its fixed ratified category",
  "DUAL_ENGINE_TERMS" ->
    "drop the locally computed spread from the canonical inventory",
  "IRREDUCIBLE_COUNT_RANGE" ->
    "corrupt the assembled object's private spread field",
  "COUNT_RANGE_REDERIVATION" ->
    "drop hbar from a private source bucket and re-run the full pipeline"
|>;

If[Length[toothOrder] =!= 20,
  Print["FAIL: stage043 tooth declaration is not exactly 20"];
  Exit[1]
];

raise[message_] := Throw[message, "ledgerStage043Failure"];

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

associationEqual[left_Association, right_Association] :=
  SameQ[KeySort[left], KeySort[right]];


(* ---------------------------------------------------------------------- *)
(* Self-contained register-to-count facts.                                *)
(* ---------------------------------------------------------------------- *)

catContinuous = "continuous-counted-irreducible";
catDebt = "reduction-debt-counted-once";
catDerived = "derived-not-counted";
catBridge = "calibrated-external-bridge-not-counted";
catDeparture = "departure-no-knob";
catConvention = "convention";
catDiscrete = "discrete-structural-postulate";
catOpen = "extension-convention-open";
catR49 = "parallel-track-out-of-scope";
catStructural = "structural-no-knob";

statusToCategory = <|
  "ACTION" -> catContinuous,
  "FREE-UNREDUCED-ROUTELESS" -> catContinuous,
  "CALIB" -> catContinuous,
  "CALIB-ANCHOR" -> catContinuous,
  "FREE-UNREDUCED-DEBT" -> catDebt,
  "DERIVED" -> catDerived,
  "EXTERNAL-BRIDGE" -> catBridge,
  "DEPARTURE" -> catDeparture,
  "CONV" -> catConvention,
  "STRUCTURAL-POSTULATE" -> catDiscrete,
  "DERIVED-IN-FORM-UNEXECUTED" -> catOpen,
  "DOWNSTREAM-DEFERRED-OPEN" -> catOpen,
  "PARALLEL-TRACK" -> catR49,
  "STRUCTURAL" -> catStructural
|>;

expectedCategoryCounts = <|
  catContinuous -> 22,
  catDebt -> 18,
  catDerived -> 40,
  catBridge -> 10,
  catDeparture -> 4,
  catConvention -> 2,
  catDiscrete -> 11,
  catOpen -> 9,
  catR49 -> 8,
  catStructural -> 28
|>;

makeFacts[
  identifiers_List,
  status_String,
  category_String,
  metadata_Association : <||>
] := Map[
  Function[identifier,
    Join[
      <|
        "id" -> identifier,
        "status" -> status,
        "expectedCategory" -> category,
        "part" -> "",
        "bucket" -> "",
        "registerLabel" -> "",
        "rider" -> "",
        "frozenCalibrationInput" -> False,
        "momentIntegralExecuted" -> True,
        "shorthand" -> ""
      |>,
      metadata
    ]
  ],
  identifiers
];

manifestFacts = Join[
  (* Parts I--II low buckets: 4 + 15 + 14 + 1 = 34. *)
  makeFacts[
    {"REG:a:hbar", "REG:a:m_GNLS", "REG:a:K_EOS", "REG:a:rho0"},
    "ACTION", catContinuous,
    <|"part" -> "I-II", "bucket" -> "a"|>
  ],
  makeFacts[
    {
      "REG:b:mu_R", "REG:b:rho_br", "REG:b:mu_R_4D",
      "REG:b:d_slab", "REG:b:mu_eta", "REG:b:T_w",
      "REG:b:beta", "REG:b:Vp0_over_lc", "REG:b:T_Omega",
      "REG:b:beta2_profile", "REG:b:epsilon0", "REG:b:epsilon1",
      "REG:b:K0c", "REG:b:K1_eta_direction",
      "REG:b:K1_TOmega_direction"
    },
    "FREE-UNREDUCED-DEBT", catDebt,
    <|"part" -> "I-II", "bucket" -> "b"|>
  ],
  makeFacts[
    {
      "REG:c:a_B", "REG:c:kappa_B", "REG:c:Omega_w",
      "REG:c:g_l_width", "REG:c:C_E", "REG:c:C_B",
      "REG:c:W_slab"
    },
    "FREE-UNREDUCED-ROUTELESS", catContinuous,
    <|"part" -> "I-II", "bucket" -> "c"|>
  ],
  makeFacts[
    {
      "REG:disc:EOS_exponent_5",
      "REG:disc:chi_B_field",
      "REG:disc:Gamma_B_law",
      "REG:disc:chi_B_gating",
      "REG:disc:G0_postulate_1_axis_surface",
      "REG:disc:G0_postulate_2_free_slip",
      "REG:disc:G0_postulate_6_no_C5_phi"
    },
    "STRUCTURAL-POSTULATE", catDiscrete,
    <|"part" -> "I-II", "bucket" -> "c"|>
  ],
  makeFacts[
    {"REG:force:force_magnitude_norm"},
    "CALIB", catContinuous,
    <|"part" -> "I-II", "bucket" -> "force-mag"|>
  ],

  (* Parts III--VI continuous extension: 5 + 7 + 1 + 0 = 13. *)
  makeFacts[
    {
      "REG:E:III:rho_B0", "REG:E:III:chi_c", "REG:E:III:B",
      "REG:E:III:J", "REG:E:III:K_theta_over_kappa_phase"
    },
    "FREE-UNREDUCED-ROUTELESS", catContinuous,
    <|"part" -> "III"|>
  ],
  makeFacts[
    {
      "REG:E:IV:C_hu", "REG:E:IV:M4", "REG:E:IV:ell_over_a",
      "REG:E:IV:K_m", "REG:E:IV:J_m"
    },
    "ACTION", catContinuous,
    <|"part" -> "IV"|>
  ],
  makeFacts[
    {"REG:E:IV:c_E", "REG:E:IV:Q_E"},
    "FREE-UNREDUCED-DEBT", catDebt,
    <|"part" -> "IV"|>
  ],
  makeFacts[
    {"REG:E:V:q_T"},
    "FREE-UNREDUCED-DEBT", catDebt,
    <|"part" -> "V"|>
  ],

  (* Frozen D1--D4, separate from the continuous variety. *)
  makeFacts[
    {
      "REG:disc:D1:H_existence_and_Poschl_Teller_law",
      "REG:disc:D2:s_i_Z2_topology",
      "REG:disc:D3:R63_mouth_BC_class",
      "REG:disc:D4:R68_tau_d_time_arrow"
    },
    "STRUCTURAL-POSTULATE", catDiscrete,
    <|"part" -> "III-VI"|>
  ],

  (* Complete not-counted bulk of the 152-entry bounded manifest. *)
  makeFacts[
    {
      "REG:derived:c_s0", "REG:derived:xi_h", "REG:derived:h0",
      "REG:derived:lambda_gamma", "REG:derived:K_eta",
      "REG:derived:delta_wall", "REG:derived:sigma_wall",
      "REG:derived:Z_return", "REG:derived:C_J", "REG:derived:B_eff",
      "REG:derived:M_h", "REG:derived:K_h", "REG:derived:Q_chi",
      "REG:derived:xi_no_sqrt2", "REG:derived:R_mouth_cancelled",
      "REG:derived:Z0_ret_alias", "REG:derived:Z1_ret_alias",
      "REG:derived:K1_sum", "REG:derived:Mtilde_definition",
      "REG:derived:Ktilde_definition", "REG:derived:Ttilde_definition",
      "REG:derived:rho_eff_projection", "REG:derived:chi_Q_outgoing",
      "REG:derived:lambda_m_SO3", "REG:derived:u2_DtN",
      "REG:derived:u4_DtN", "REG:derived:v5_DtN",
      "REG:derived:Kbar4_moment", "REG:derived:p_localization",
      "REG:derived:force_sign", "REG:derived:K4_from_M4_cE",
      "REG:derived:CJ_IBP_sign", "REG:derived:charge_value_R57",
      "REG:derived:qT_product", "REG:derived:mouth_class_fillers",
      "REG:derived:wall_kink_profile", "REG:derived:round_trip_Rrt",
      "REG:derived:DtN_pole_ladder", "REG:derived:cone_ratio",
      "REG:derived:field_overlay_determinant"
    },
    "DERIVED", catDerived
  ],
  makeFacts[
    {
      "REG:bridge:G_Newton", "REG:bridge:PN_54_over_5",
      "REG:bridge:PN_2_over_5", "REG:bridge:benchmark_force",
      "REG:bridge:benchmark_orbit", "REG:bridge:benchmark_flux",
      "REG:bridge:external_mass_anchor",
      "REG:bridge:external_charge_anchor",
      "REG:bridge:external_time_anchor",
      "REG:bridge:external_length_anchor"
    },
    "EXTERNAL-BRIDGE", catBridge
  ],
  makeFacts[
    {
      "REG:departure:R66_native_P_no_emergent_Gauss",
      "REG:departure:R73_bT_T_even",
      "REG:departure:R81_scalar_sim_gated",
      "REG:departure:light_stray_longitudinal"
    },
    "DEPARTURE", catDeparture
  ],
  makeFacts[
    {"REG:conv:a_pin", "REG:conv:lambda_gamma_ratio"},
    "CONV", catConvention
  ],
  makeFacts[
    {
      "REG:R49:profile_E0", "REG:R49:profile_E1",
      "REG:R49:profile_E2", "REG:R49:profile_E3",
      "REG:R49:profile_B0", "REG:R49:profile_B1",
      "REG:R49:profile_B2", "REG:R49:profile_B3"
    },
    "PARALLEL-TRACK", catR49
  ],
  makeFacts[
    {
      "REG:struct:R23_return_targets", "REG:struct:R25_slab_survival",
      "REG:struct:R26_frozen_validity", "REG:struct:R27_xi_lc_firewall",
      "REG:struct:R28_BC_dependent", "REG:struct:R31_truncation_window",
      "REG:struct:R42_scalar_reduction", "REG:struct:R50_H_parent",
      "REG:struct:R53_kernel_stability", "REG:struct:R63_mouth_data",
      "REG:struct:R65_charge_blocker", "REG:struct:R67_drain_coupling",
      "REG:struct:R69_return_bookkeeping", "REG:struct:R71_open_rcone",
      "REG:struct:R74_NG5_trio", "REG:struct:R75_shared_throat_route",
      "REG:struct:R76_QE_sensitivity", "REG:struct:R77_lineage_seal",
      "REG:struct:R78_cone_adjudication", "REG:struct:R79_OPEN_110",
      "REG:struct:R80_freedom_conditionality",
      "REG:struct:R82_falsifier_set", "REG:struct:R83_guard_A",
      "REG:struct:R84_VI_VII_handoff", "REG:struct:z_g_opaque",
      "REG:struct:z_b_opaque", "REG:struct:control_inventory",
      "REG:struct:census_bookkeeping"
    },
    "STRUCTURAL", catStructural
  ]
];

c1Ids = {"REG:C1:Mtilde", "REG:C1:Ktilde", "REG:C1:Ttilde_Omega"};
c2Specs = {
  {"REG:C2:Btilde0", "Btilde"},
  {"REG:C2:Btilde2", "Btilde"},
  {"REG:C2:Btilde4", "Btilde"},
  {"REG:C2:Ztilde0", "Ztilde"},
  {"REG:C2:Ztilde2", "Ztilde"},
  {"REG:C2:Ztilde4", "Ztilde"}
};

manifestFacts = Join[
  manifestFacts,
  makeFacts[
    c1Ids,
    "DERIVED-IN-FORM-UNEXECUTED", catOpen,
    <|
      "part" -> "I-II",
      "bucket" -> "C1",
      "registerLabel" -> "DERIVED-in-form",
      "rider" -> "C1_contributor",
      "frozenCalibrationInput" -> True,
      "momentIntegralExecuted" -> False,
      "shorthand" -> "R35"
    |>
  ],
  Map[
    Function[spec,
      First@makeFacts[
        {spec[[1]]},
        "DOWNSTREAM-DEFERRED-OPEN", catOpen,
        <|
          "part" -> "I-II",
          "bucket" -> "C2",
          "registerLabel" -> "downstream-deferred",
          "rider" -> "C2_contributor",
          "shorthand" -> spec[[2]]
        |>
      ]
    ],
    c2Specs
  ]
];

categoryOf[row_Association] :=
  Lookup[statusToCategory, row["status"], "UNCLASSIFIED"];

replaceFact[rows_List, identifier_String, changes_Association] :=
  Map[
    If[#["id"] === identifier, Join[#, changes], #] &,
    rows
  ];

rowIndex[rows_List] := AssociationThread[Lookup[rows, "id"], rows];

idsOf[rows_List] := If[rows === {}, {}, Lookup[rows, "id"]];

manifestState[rows_List] := Module[
  {
    ids, computedMap, expectedMap, counts, unique, engineQualified,
    mapOK, countsOK
  },
  ids = Lookup[rows, "id"];
  computedMap = AssociationThread[ids, categoryOf /@ rows];
  expectedMap = AssociationThread[ids, Lookup[rows, "expectedCategory"]];
  counts = Counts[Values[computedMap]];
  unique = DuplicateFreeQ[ids] && Length[ids] == 152;
  engineQualified = And @@ Map[
    MatchQ[StringSplit[#, ":"], {"REG", __String}] &,
    ids
  ];
  mapOK = associationEqual[computedMap, expectedMap];
  countsOK = associationEqual[counts, expectedCategoryCounts];
  <|
    "unique" -> unique,
    "engineQualified" -> engineQualified,
    "perRowMap" -> mapOK,
    "counts" -> counts,
    "countsMatch" -> countsOK,
    "valid" -> TrueQ[unique && engineQualified && mapOK && countsOK]
  |>
];

bucketSets[rows_List] := Module[{names, grouped},
  names = {"a", "b", "c", "force-mag", "C1", "C2"};
  grouped = GroupBy[
    Select[rows, MemberQ[names, #["bucket"]] &],
    #["bucket"] &,
    Lookup[#, "id"] &
  ];
  AssociationMap[DeleteDuplicates@Lookup[grouped, #, {}] &, names]
];

bucketPartitionState[buckets_Association] := Module[
  {expected, counts, flattened, tally, rawLow, rawHigh, disjoint},
  expected = <|
    "a" -> 4, "b" -> 15, "c" -> 14,
    "force-mag" -> 1, "C1" -> 3, "C2" -> 6
  |>;
  counts = AssociationMap[Length[buckets[#]] &, Keys[expected]];
  flattened = Flatten[Lookup[buckets, Keys[expected]]];
  tally = Tally[flattened];
  disjoint = AllTrue[tally, Last[#] == 1 &];
  rawLow = Total[Lookup[counts, {"a", "b", "c", "force-mag"}]];
  rawHigh = Total[Values[counts]];
  <|
    "counts" -> counts,
    "rawLow" -> rawLow,
    "rawHigh" -> rawHigh,
    "disjoint" -> disjoint,
    "valid" -> TrueQ[
      associationEqual[counts, expected] &&
      disjoint &&
      rawLow == 34 &&
      rawHigh == 43
    ]
  |>
];

deriveCountState[rows_List, forceStrict_: False] := Module[
  {
    buckets, bucketState, categories, baseDiscrete,
    parts, extensionByPart, extensionIds,
    lowRawIds, highRawIds, effectiveLowRawIds,
    lowBaseIds, highBaseIds, lowComponents, highComponents,
    eCounts
  },
  buckets = bucketSets[rows];
  bucketState = bucketPartitionState[buckets];
  categories = AssociationThread[Lookup[rows, "id"], categoryOf /@ rows];
  baseDiscrete = idsOf[
    Select[
      rows,
      #["part"] === "I-II" && categoryOf[#] === catDiscrete &
    ]
  ];
  parts = {"III", "IV", "V", "VI"};
  extensionByPart = AssociationMap[
    Function[part,
      idsOf[
        Select[
          rows,
          #["part"] === part &&
            MemberQ[{catContinuous, catDebt}, categoryOf[#]] &
        ]
      ]
    ],
    parts
  ];
  extensionIds = DeleteDuplicates@Flatten[Values[extensionByPart]];
  lowRawIds = Union @@ Lookup[buckets, {"a", "b", "c", "force-mag"}];
  highRawIds = Union[lowRawIds, buckets["C1"], buckets["C2"]];
  effectiveLowRawIds = If[TrueQ[forceStrict], highRawIds, lowRawIds];
  lowBaseIds = Complement[effectiveLowRawIds, baseDiscrete];
  highBaseIds = Complement[highRawIds, baseDiscrete];
  lowComponents = Join[Sort[lowBaseIds], Sort[extensionIds]];
  highComponents = Join[Sort[highBaseIds], Sort[extensionIds]];
  eCounts = AssociationMap[Length[extensionByPart[#]] &, parts];
  <|
    "buckets" -> buckets,
    "bucketState" -> bucketState,
    "rawLow" -> Length[effectiveLowRawIds],
    "rawHigh" -> Length[highRawIds],
    "baseDiscreteIds" -> baseDiscrete,
    "baseContinuous" -> {Length[lowBaseIds], Length[highBaseIds]},
    "extensionByPart" -> extensionByPart,
    "extensionIds" -> extensionIds,
    "EItemization" -> eCounts,
    "EContinuous" -> Total[Values[eCounts]],
    "lowComponents" -> lowComponents,
    "highComponents" -> highComponents,
    "lowC" -> Length[lowComponents],
    "highC" -> Length[highComponents],
    "spread" -> Length[highComponents] - Length[lowComponents],
    "categories" -> categories
  |>
];


(* ---------------------------------------------------------------------- *)
(* Mathematica-native CAD RegionDimension diagnostic.                     *)
(* ---------------------------------------------------------------------- *)

mediumVars = {
  hbarM, massM, bigK, rho0M, cs0M,
  xiHM, h0M, scaleAM, cGammaM, lambdaGammaM
};
mediumBaseline = {
  massM*cs0M^2 - 5*bigK*rho0M^4,
  scaleAM*massM*cs0M - hbarM,
  xiHM^2*massM^2*cs0M^2 - 2*hbarM^2,
  4*h0M - massM*cs0M^2,
  lambdaGammaM*cs0M - cGammaM
};

wallVars = {
  muEtaM, tWM, betaM, kEtaM, vp0LcM,
  tOmegaM, aBM, kappaBM, deltaM, sigmaWallM
};
wallBaseline = {
  kEtaM - tWM*betaM^2,
  2*aBM*deltaM^2 - kappaBM,
  36*sigmaWallM^2 - 2*aBM*kappaBM
};

blocks = <|
  "M" -> <|
    "vars" -> mediumVars,
    "baseline" -> mediumBaseline,
    "forced" -> bigK - rho0M
  |>,
  "Wchi" -> <|
    "vars" -> wallVars,
    "baseline" -> wallBaseline,
    "forced" -> muEtaM - tWM
  |>
|>;

regionDimensionCAD[vars_List, equations_List] :=
  regionDimensionCAD[vars, equations] = Module[
    {nonzero, formula, region, result},
    nonzero = DeleteCases[Expand /@ equations, 0];
    formula = And @@ Join[
      Thread[nonzero == 0],
      Thread[vars > 0]
    ];
    region = ImplicitRegion[formula, Evaluate[vars]];
    result = TimeConstrained[
      RegionDimension[region],
      120,
      "DIMENSION_UNCERTIFIED"
    ];
    If[! IntegerQ[result],
      Print["DIMENSION_UNCERTIFIED: equations=", equations];
      raise["CAD_REGION_DIMENSION"]
    ];
    result
  ];

dimensionRecordCAD[block_Association, equations_List] := Module[
  {before, after},
  before = regionDimensionCAD[block["vars"], {}];
  after = regionDimensionCAD[block["vars"], equations];
  <|
    "dim_before" -> before,
    "dim_after" -> after,
    "delta_r" -> before - after
  |>
];

expectedDeltaBaseline = <|
  "M" -> <|"dim_before" -> 10, "dim_after" -> 5, "delta_r" -> 5|>,
  "Wchi" -> <|"dim_before" -> 10, "dim_after" -> 7, "delta_r" -> 3|>
|>;
expectedDeltaForced = <|
  "M" -> <|"dim_before" -> 10, "dim_after" -> 4, "delta_r" -> 6|>,
  "Wchi" -> <|"dim_before" -> 10, "dim_after" -> 6, "delta_r" -> 4|>
|>;

deltaDiagnostic[forceBaselineDependency_: False] := Module[
  {baseline, forced, guards, token},
  baseline = AssociationMap[
    Function[name,
      With[{block = blocks[name]},
        dimensionRecordCAD[
          block,
          If[
            TrueQ[forceBaselineDependency],
            Append[block["baseline"], block["forced"]],
            block["baseline"]
          ]
        ]
      ]
    ],
    Keys[blocks]
  ];
  forced = AssociationMap[
    Function[name,
      With[{block = blocks[name]},
        dimensionRecordCAD[
          block,
          Append[block["baseline"], block["forced"]]
        ]
      ]
    ],
    Keys[blocks]
  ];
  guards = TrueQ[
    associationEqual[baseline, expectedDeltaBaseline] &&
    associationEqual[forced, expectedDeltaForced] &&
    forced["M"]["delta_r"] > baseline["M"]["delta_r"] &&
    forced["Wchi"]["delta_r"] > baseline["Wchi"]["delta_r"]
  ];
  token = If[
    guards,
    "CONFIRMED_IN_TESTED_M_AND_WCHI_BLOCKS",
    "DEVIATION_DETECTED"
  ];
  <|
    "baseline" -> baseline,
    "forced" -> forced,
    "guards" -> guards,
    "token" -> token
  |>
];


(* ---------------------------------------------------------------------- *)
(* Headline assembly.                                                      *)
(* ---------------------------------------------------------------------- *)

deriveQESensitivity[rows_List, state_Association] := Module[
  {qeIds, trioIds, declaredSplit, undeclaredSplit, before, after, neutral},
  qeIds = Lookup[
    Select[
      rows,
      #["id"] === "REG:E:IV:Q_E" && categoryOf[#] === catDebt &
    ],
    "id"
  ];
  trioIds = {
    "REG:E:III:rho_B0", "REG:E:III:chi_c", "REG:E:IV:C_hu"
  };
  declaredSplit = Length@Select[rows, MemberQ[trioIds, #["id"]] &];
  undeclaredSplit = declaredSplit + Length[qeIds];
  before = Length[state["lowComponents"]];
  after = Length@Union[state["lowComponents"], qeIds];
  neutral = before == after;
  <|
    "declaredSplit" -> declaredSplit,
    "undeclaredSplit" -> undeclaredSplit,
    "totalNeutral" -> neutral,
    "token" -> If[
      declaredSplit == 3 && undeclaredSplit == 4 && neutral,
      "TOTAL_NEUTRAL_ROUTELESS_SHIFT_3_to_4",
      "QE_SENSITIVITY_DEVIATION"
    ]
  |>
];

deriveR49Scope[rows_List, state_Association] := Module[
  {profiles, overlap, valid},
  profiles = Lookup[Select[rows, categoryOf[#] === catR49 &], "id"];
  overlap = Intersection[profiles, state["highComponents"]];
  valid = Length[profiles] == 8 && overlap === {};
  <|
    "profiles" -> profiles,
    "overlap" -> overlap,
    "valid" -> valid,
    "token" -> If[
      valid,
      "OUT_OF_BUILT_V2_SCOPE",
      "R49_SCOPE_VIOLATION"
    ]
  |>
];

deriveConvention[state_Association] := Module[{openMembers},
  openMembers = Union[state["buckets"]["C1"], state["buckets"]["C2"]];
  If[Length[openMembers] > 0 && state["lowC"] < state["highC"], "OPEN", "IMPOSED"]
];

assembleRange[
  rows_List,
  state_Association,
  diagnostic_Association
] := Module[{discreteCount, qe, r49},
  discreteCount = Length@Select[rows, categoryOf[#] === catDiscrete &];
  qe = deriveQESensitivity[rows, state];
  r49 = deriveR49Scope[rows, state];
  <|
    "continuous_codimension" -> {state["lowC"], state["highC"]},
    "base_continuous" -> state["baseContinuous"],
    "E_continuous" -> state["EContinuous"],
    "discrete_postulate_count" -> discreteCount,
    "convention" -> deriveConvention[state],
    "spread" -> state["spread"],
    "C1_cardinality" -> Length[state["buckets"]["C1"]],
    "C2_cardinality" -> Length[state["buckets"]["C2"]],
    "delta_r_independence" -> diagnostic["token"],
    "Q_E_undeclared" -> qe["token"],
    "R49_parallel_track" -> r49["token"]
  |>
];

expectedRangeObject = <|
  "continuous_codimension" -> {40, 49},
  "base_continuous" -> {27, 36},
  "E_continuous" -> 13,
  "discrete_postulate_count" -> 11,
  "convention" -> "OPEN",
  "spread" -> 9,
  "C1_cardinality" -> 3,
  "C2_cardinality" -> 6,
  "delta_r_independence" -> "CONFIRMED_IN_TESTED_M_AND_WCHI_BLOCKS",
  "Q_E_undeclared" -> "TOTAL_NEUTRAL_ROUTELESS_SHIFT_3_to_4",
  "R49_parallel_track" -> "OUT_OF_BUILT_V2_SCOPE"
|>;

rangeObjectValid[object_Association] := Module[{low, high},
  {low, high} = object["continuous_codimension"];
  TrueQ[
    associationEqual[object, expectedRangeObject] &&
    low < high &&
    high - low == 9 &&
    object["spread"] ==
      object["C1_cardinality"] + object["C2_cardinality"]
  ]
];


runAssertions[] := Module[
  {
    baseRows, state, diagnostic, production,
    privateBuckets, bucketActual,
    lowRows, lowActual, highRows, highActual,
    spreadRows, spreadC1Count, spreadC2Count,
    baseRowsPrivate, baseActual, eRows, eActual,
    c1Rows, c1Members, c2Rows, c2Members,
    conventionObject, conventionLow, conventionHigh,
    debtIds, debtComponents, debtTally, debtOccurrences,
    excludedControl, doubledControl,
    r49, r49Components, r49Overlap,
    discreteRows, discreteState, discreteIds,
    discreteBaseIds, discreteExtensionIds, continuousUnion,
    deltaActual, countWithDiagnostic,
    r35Rows, r35Labels, r35Riders,
    index, countedForProvenance, allowedStatuses, badProvenance,
    manifestRows, manifestActual,
    qe, inventory, expectedInventory,
    rangePrivate, rederiveRows, rederivedState,
    rederivedObject, strictState, positiveForceStrict
  },

  baseRows = manifestFacts;
  state = deriveCountState[baseRows];
  diagnostic = deltaDiagnostic[];
  production = assembleRange[baseRows, state, diagnostic];

  section["A. Parts-I--II buckets and continuous endpoints"];

  privateBuckets = Association @@ Normal[bucketSets[baseRows]];
  If[activeMutation === "BUCKET_PARTITION",
    privateBuckets["b"] = Append[privateBuckets["b"], "REG:a:hbar"]
  ];
  bucketActual = bucketPartitionState[privateBuckets];
  expectBool["BUCKET_PARTITION", bucketActual["valid"], bucketActual];

  lowRows = baseRows;
  If[activeMutation === "LOW_ENDPOINT",
    lowRows = replaceFact[lowRows, "REG:a:hbar", <|"bucket" -> ""|>]
  ];
  lowActual = deriveCountState[lowRows];
  expectBool[
    "LOW_ENDPOINT",
    lowActual["lowC"] == 40,
    <|
      "lowC" -> lowActual["lowC"],
      "rawLow" -> lowActual["rawLow"],
      "baseDiscrete" -> lowActual["baseDiscreteIds"],
      "baseContinuous" -> lowActual["baseContinuous"],
      "EItemization" -> lowActual["EItemization"],
      "extensionByPart" -> lowActual["extensionByPart"]
    |>
  ];

  highRows = baseRows;
  If[activeMutation === "HIGH_ENDPOINT",
    highRows = replaceFact[
      highRows, "REG:C2:Btilde0", <|"bucket" -> ""|>
    ]
  ];
  highActual = deriveCountState[highRows];
  expectBool[
    "HIGH_ENDPOINT",
    highActual["highC"] == 49,
    <|"highC" -> highActual["highC"], "rawHigh" -> highActual["rawHigh"]|>
  ];

  spreadRows = baseRows;
  If[activeMutation === "RANGE_SPREAD_IS_C1_PLUS_C2",
    spreadRows = replaceFact[
      spreadRows, "REG:C1:Mtilde", <|"bucket" -> ""|>
    ]
  ];
  spreadC1Count = Length[bucketSets[spreadRows]["C1"]];
  spreadC2Count = Length[bucketSets[baseRows]["C2"]];
  expectBool[
    "RANGE_SPREAD_IS_C1_PLUS_C2",
    state["spread"] == 9 &&
      state["spread"] == spreadC1Count + spreadC2Count,
    <|
      "spread" -> state["spread"],
      "C1" -> spreadC1Count,
      "C2" -> spreadC2Count
    |>
  ];

  baseRowsPrivate = baseRows;
  If[activeMutation === "BASE_CONTINUOUS_27_36",
    baseRowsPrivate = replaceFact[
      baseRowsPrivate,
      "REG:disc:EOS_exponent_5",
      <|"status" -> "STRUCTURAL"|>
    ]
  ];
  baseActual = deriveCountState[baseRowsPrivate];
  expectBool[
    "BASE_CONTINUOUS_27_36",
    baseActual["baseContinuous"] === {27, 36} &&
      Length[baseActual["baseDiscreteIds"]] == 7 &&
      baseActual["rawLow"] == 34 &&
      baseActual["rawHigh"] == 43,
    <|
      "base" -> baseActual["baseContinuous"],
      "pullout" -> Length[baseActual["baseDiscreteIds"]],
      "raw" -> {baseActual["rawLow"], baseActual["rawHigh"]}
    |>
  ];

  eRows = baseRows;
  If[activeMutation === "E_CONTINUOUS_IS_13",
    eRows = replaceFact[eRows, "REG:E:III:rho_B0", <|"part" -> ""|>]
  ];
  eActual = deriveCountState[eRows];
  expectBool[
    "E_CONTINUOUS_IS_13",
    associationEqual[
      eActual["EItemization"],
      <|"III" -> 5, "IV" -> 7, "V" -> 1, "VI" -> 0|>
    ] &&
      eActual["EContinuous"] == 13,
    <|
      "itemization" -> eActual["EItemization"],
      "EContinuous" -> eActual["EContinuous"]
    |>
  ];

  section["B. OPEN convention contributors"];

  c1Rows = baseRows;
  If[activeMutation === "C1_ATTRIBUTION",
    c1Rows = replaceFact[c1Rows, "REG:C1:Mtilde", <|"bucket" -> ""|>]
  ];
  c1Members = Lookup[
    Select[
      c1Rows,
      #["bucket"] === "C1" &&
        TrueQ[#["frozenCalibrationInput"]] &&
        ! TrueQ[#["momentIntegralExecuted"]] &
    ],
    "id"
  ];
  expectBool[
    "C1_ATTRIBUTION",
    Sort[c1Members] === Sort[c1Ids] && Length[c1Members] == 3,
    c1Members
  ];

  c2Rows = Select[baseRows, #["bucket"] === "C2" &];
  c2Members = If[
    activeMutation === "C2_ATTRIBUTION",
    DeleteDuplicates[Lookup[c2Rows, "shorthand"]],
    Lookup[c2Rows, "id"]
  ];
  expectBool[
    "C2_ATTRIBUTION",
    Sort[c2Members] === Sort[c2Specs[[All, 1]]] &&
      Length[c2Members] == 6,
    c2Members
  ];

  conventionObject = Association @@ Normal[production];
  If[activeMutation === "CONVENTION_OPEN_NOT_IMPOSED",
    conventionObject["continuous_codimension"] =
      {state["highC"], state["highC"]}
  ];
  {conventionLow, conventionHigh} =
    conventionObject["continuous_codimension"];
  expectBool[
    "CONVENTION_OPEN_NOT_IMPOSED",
    conventionLow < conventionHigh &&
      Lookup[conventionObject, "convention", "MISSING"] === "OPEN",
    conventionObject
  ];

  section["C. Debt, discrete, and out-of-scope accounting"];

  debtIds = {
    "REG:b:mu_R", "REG:b:rho_br",
    "REG:E:IV:c_E", "REG:E:IV:Q_E", "REG:E:V:q_T"
  };
  debtComponents = state["lowComponents"];
  If[activeMutation === "REDUCTION_DEBT_COUNTED_ONCE",
    debtComponents = DeleteCases[debtComponents, "REG:b:mu_R"]
  ];
  debtTally = Counts[debtComponents];
  debtOccurrences = AssociationMap[Lookup[debtTally, #, 0] &, debtIds];
  excludedControl = Select[
    state["lowComponents"], ! MemberQ[debtIds, #] &
  ];
  doubledControl = Join[
    state["lowComponents"], {"REG:b:mu_R", "REG:b:rho_br"}
  ];
  expectBool[
    "REDUCTION_DEBT_COUNTED_ONCE",
    Length[debtComponents] == 40 &&
      And @@ Thread[Values[debtOccurrences] == 1] &&
      Length[excludedControl] == 35 &&
      Length[doubledControl] == 42,
    <|
      "productionLow" -> Length[debtComponents],
      "occurrences" -> debtOccurrences,
      "excludeControl" -> Length[excludedControl],
      "doubleControl" -> Length[doubledControl]
    |>
  ];

  r49 = deriveR49Scope[baseRows, state];
  r49Components = state["lowComponents"];
  If[activeMutation === "R49_OUT_OF_SCOPE",
    r49Components = Append[r49Components, First@Sort[r49["profiles"]]]
  ];
  r49Overlap = Intersection[r49Components, r49["profiles"]];
  expectBool[
    "R49_OUT_OF_SCOPE",
    Length[r49["profiles"]] == 8 &&
      Length[r49Components] == 40 &&
      r49Overlap === {} &&
      r49["token"] === "OUT_OF_BUILT_V2_SCOPE",
    <|
      "profileCount" -> Length[r49["profiles"]],
      "continuousCount" -> Length[r49Components],
      "overlap" -> r49Overlap
    |>
  ];

  discreteRows = baseRows;
  If[activeMutation === "DISCRETE_POSTULATE_COUNT",
    discreteRows = replaceFact[
      discreteRows,
      "REG:disc:D1:H_existence_and_Poschl_Teller_law",
      <|"status" -> "STRUCTURAL"|>
    ]
  ];
  discreteState = deriveCountState[discreteRows];
  discreteIds = Lookup[
    Select[discreteRows, categoryOf[#] === catDiscrete &],
    "id"
  ];
  discreteBaseIds = Lookup[
    Select[
      discreteRows,
      categoryOf[#] === catDiscrete && #["part"] === "I-II" &
    ],
    "id"
  ];
  discreteExtensionIds = Complement[discreteIds, discreteBaseIds];
  continuousUnion = Union[discreteState["highComponents"]];
  expectBool[
    "DISCRETE_POSTULATE_COUNT",
    Length[discreteIds] == 11 &&
      Length[discreteBaseIds] == 7 &&
      Length[discreteExtensionIds] == 4 &&
      Intersection[discreteIds, continuousUnion] === {},
    <|
      "total" -> Length[discreteIds],
      "base" -> Length[discreteBaseIds],
      "DIII_VI" -> Length[discreteExtensionIds],
      "continuousOverlap" -> Intersection[discreteIds, continuousUnion]
    |>
  ];

  section["D. Scoped delta-r independence diagnostic"];

  deltaActual = deltaDiagnostic[
    activeMutation === "DELTA_R_INDEPENDENCE"
  ];
  expectBool[
    "DELTA_R_INDEPENDENCE",
    deltaActual["guards"] &&
      deltaActual["token"] ===
        "CONFIRMED_IN_TESTED_M_AND_WCHI_BLOCKS",
    deltaActual
  ];

  countWithDiagnostic = state["lowC"];
  If[activeMutation === "DELTA_R_IS_DIAGNOSTIC_NOT_SUBTRACTION",
    countWithDiagnostic -= Total[
      Lookup[Values[diagnostic["baseline"]], "delta_r"]
    ]
  ];
  expectBool[
    "DELTA_R_IS_DIAGNOSTIC_NOT_SUBTRACTION",
    countWithDiagnostic == 40 &&
      countWithDiagnostic == Length[state["lowComponents"]] &&
      diagnostic["token"] ===
        "CONFIRMED_IN_TESTED_M_AND_WCHI_BLOCKS",
    <|
      "withDiagnostic" -> countWithDiagnostic,
      "partitionOnly" -> Length[state["lowComponents"]]
    |>
  ];

  section["E. R35, provenance, manifest, and local engine inventory"];

  r35Rows = baseRows;
  If[activeMutation === "R35_DERIVED_IN_FORM_PLUS_RIDER",
    r35Rows = replaceFact[
      r35Rows,
      "REG:C1:Mtilde",
      <|"registerLabel" -> "PENDING-debt"|>
    ]
  ];
  r35Labels = DeleteDuplicates@Lookup[
    Select[r35Rows, #["bucket"] === "C1" &],
    "registerLabel"
  ];
  r35Riders = DeleteDuplicates@Lookup[
    Select[r35Rows, #["bucket"] === "C1" &],
    "rider"
  ];
  expectBool[
    "R35_DERIVED_IN_FORM_PLUS_RIDER",
    r35Labels === {"DERIVED-in-form"} &&
      r35Riders === {"C1_contributor"},
    <|"labels" -> r35Labels, "riders" -> r35Riders|>
  ];

  index = rowIndex[baseRows];
  countedForProvenance = Union[state["lowComponents"]];
  If[activeMutation === "KNOB_PROVENANCE_WELL_FORMED",
    countedForProvenance = Append[
      countedForProvenance, "REG:derived:c_s0"
    ]
  ];
  allowedStatuses = {
    "ACTION", "FREE-UNREDUCED-ROUTELESS",
    "FREE-UNREDUCED-DEBT", "CALIB", "CALIB-ANCHOR"
  };
  badProvenance = Association@Cases[
    countedForProvenance,
    identifier_ /; ! MemberQ[
      allowedStatuses, index[identifier]["status"]
    ] :> identifier -> index[identifier]["status"]
  ];
  expectBool[
    "KNOB_PROVENANCE_WELL_FORMED",
    Length[badProvenance] == 0 &&
      Length[countedForProvenance] == 40,
    badProvenance
  ];

  manifestRows = baseRows;
  If[activeMutation === "REGISTER_TO_COUNT_MANIFEST",
    manifestRows = replaceFact[
      manifestRows, "REG:a:hbar", <|"status" -> "DERIVED"|>
    ]
  ];
  manifestActual = manifestState[manifestRows];
  expectBool[
    "REGISTER_TO_COUNT_MANIFEST",
    manifestActual["valid"],
    manifestActual
  ];

  qe = deriveQESensitivity[baseRows, state];
  inventory = <|
    "bucket_counts" -> state["bucketState"]["counts"],
    "raw_parts_I_II" -> {
      state["bucketState"]["rawLow"],
      state["bucketState"]["rawHigh"]
    },
    "base_continuous" -> state["baseContinuous"],
    "E_itemization" -> state["EItemization"],
    "E_continuous" -> state["EContinuous"],
    "continuous_codimension" -> {state["lowC"], state["highC"]},
    "spread" -> state["spread"],
    "C1_cardinality" -> Length[state["buckets"]["C1"]],
    "C2_cardinality" -> Length[state["buckets"]["C2"]],
    "discrete_postulate_count" ->
      production["discrete_postulate_count"],
    "reduction_debt_counted_once" -> Sort[debtIds],
    "Q_E_shift" -> {qe["declaredSplit"], qe["undeclaredSplit"]},
    "delta_r_result" -> diagnostic["token"],
    "IRREDUCIBLE_COUNT_RANGE" -> production
  |>;
  expectedInventory = <|
    "bucket_counts" -> <|
      "a" -> 4, "b" -> 15, "c" -> 14,
      "force-mag" -> 1, "C1" -> 3, "C2" -> 6
    |>,
    "raw_parts_I_II" -> {34, 43},
    "base_continuous" -> {27, 36},
    "E_itemization" -> <|"III" -> 5, "IV" -> 7, "V" -> 1, "VI" -> 0|>,
    "E_continuous" -> 13,
    "continuous_codimension" -> {40, 49},
    "spread" -> 9,
    "C1_cardinality" -> 3,
    "C2_cardinality" -> 6,
    "discrete_postulate_count" -> 11,
    "reduction_debt_counted_once" -> Sort[debtIds],
    "Q_E_shift" -> {3, 4},
    "delta_r_result" -> "CONFIRMED_IN_TESTED_M_AND_WCHI_BLOCKS",
    "IRREDUCIBLE_COUNT_RANGE" -> expectedRangeObject
  |>;
  If[activeMutation === "DUAL_ENGINE_TERMS",
    inventory = KeyDrop[inventory, "spread"]
  ];
  expectBool[
    "DUAL_ENGINE_TERMS",
    associationEqual[inventory, expectedInventory],
    inventory
  ];

  section["F. Assembled range and full rederivation"];

  rangePrivate = Association @@ Normal[production];
  If[activeMutation === "IRREDUCIBLE_COUNT_RANGE",
    rangePrivate["spread"] = rangePrivate["spread"] + 1
  ];
  expectBool[
    "IRREDUCIBLE_COUNT_RANGE",
    rangeObjectValid[rangePrivate],
    rangePrivate
  ];

  rederiveRows = baseRows;
  If[activeMutation === "COUNT_RANGE_REDERIVATION",
    rederiveRows = replaceFact[
      rederiveRows, "REG:a:hbar", <|"bucket" -> ""|>
    ]
  ];
  rederivedState = deriveCountState[rederiveRows];
  rederivedObject = assembleRange[
    rederiveRows, rederivedState, diagnostic
  ];
  strictState = deriveCountState[rederiveRows, True];
  positiveForceStrict =
    strictState["lowC"] == 49 &&
    strictState["lowC"] == rederivedState["highC"];
  expectBool[
    "COUNT_RANGE_REDERIVATION",
    associationEqual[rederivedObject, expectedRangeObject] &&
      positiveForceStrict,
    <|
      "rederived" -> rederivedObject,
      "forceStrictLow" -> strictState["lowC"],
      "namedHigh" -> rederivedState["highC"]
    |>
  ];

  production
];


ok = Catch[
  If[
    activeMutation =!= "" && ! MemberQ[toothOrder, activeMutation],
    Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
    Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
    raise["UNKNOWN_MUTATION"]
  ];

  Print["ledger_stage043_irreducible_count_range Mathematica audit"];
  Print[
    "ROUTE=Association/GroupBy/Counts/Select/Tally partition + ",
    "CAD RegionDimension[ImplicitRegion]"
  ];
  Print["FILE_IO=none; CROSS_ENGINE_COMPARE=none"];
  Print["DIMENSIONAL_HOMOGENEITY=N/A_INTEGER_COUNT_STAGE"];
  If[
    activeMutation =!= "",
    Print["ACTIVE_MUTATION=", activeMutation];
    Print[
      "MUTATED_PRIMITIVE=",
      ablationDescriptions[activeMutation]
    ]
  ];

  productionRange = runAssertions[];
  If[
    activeMutation =!= "",
    Print["FIRST_FAILURE=MUTATION_DID_NOT_FIRE"];
    raise["MUTATION_DID_NOT_FIRE"]
  ];
  True,
  "ledgerStage043Failure",
  Function[{message, tag}, False]
];

If[TrueQ[ok],
  Print[""];
  Print[
    "IRREDUCIBLE_COUNT_RANGE=",
    ExportString[productionRange, "RawJSON", "Compact" -> True]
  ];
  Print[
    "RAW_PARTS_I_II=[34,43]; BASE_CONTINUOUS=[27,36]; ",
    "E_CONTINUOUS=13"
  ];
  Print["CONTINUOUS_CODIMENSION=[40,49]; SPREAD=9; C1=3; C2=6"];
  Print[
    "DISCRETE_POSTULATE_COUNT=11 (=7+4); ",
    "NOT_IN_CONTINUOUS_VARIETY"
  ];
  Print["REDUCTION_DEBT_COUNTED_ONCE={mu_R,rho_br,c_E,q_T,Q_E}"];
  Print[
    "R49_PARALLEL_TRACK=OUT_OF_BUILT_V2_SCOPE; PROFILE_DEBT=8_to_9"
  ];
  Print["Q_E_UNDECLARED=TOTAL_NEUTRAL_ROUTELESS_SHIFT_3_to_4"];
  Print["K_THETA_KAPPA_PHASE_AS_TWO_SENSITIVITY=[41,50]"];
  Print[
    "DELTA_R_SCOPE=CONFIRMED_IN_TESTED_M_AND_WCHI_BLOCKS; ",
    "M:(10,5,5)->(10,4,6); Wchi:(10,7,3)->(10,6,4)"
  ];
];

Print["TOOTH_COUNT=", Length[toothOrder]];
Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
If[
  TrueQ[ok],
  Print[
    "OVERALL PASS: Mathematica independently assembled the count range"
  ];
  Exit[0],
  Print["OVERALL FAIL: Mathematica stage043 audit did not close"];
  Exit[1]
]
