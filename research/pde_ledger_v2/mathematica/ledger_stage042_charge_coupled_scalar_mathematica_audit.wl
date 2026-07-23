(* Ledger stage042 Mathematica audit: charge-coupled scalar consistency map.

   Standalone, print-only, assert-zero, and cross-engine-file-I/O-free.

   This engine is deliberately organized around Mathematica-native primitives
   that differ materially from both the stage SymPy engine and the thin source
   Wolfram mirror.  LinearSolve supplies the three sourced quadratic forms
   without Inverse.  A logarithmic cE derivative followed by Limit extracts
   radiation falloff.  Ordered condition tables are resolved by FirstCase and
   the verdict is assembled from set-membership predicates.  Guard A uses
   Position to discover numeric/string leaves and Cases to recover Association
   key paths.  Controls and deletion rows are recomputed from Associations; no
   Python file or cross-engine payload is read.

   Tooth-local runtime ablation uses LEDGER_STAGE042_MUTATION.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE042_MUTATION";
activeMutation = Quiet@Check[Environment[mutationEnvironment], ""];
If[! StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];

toothOrder = {
  "STATIC_DETERMINANT",
  "A_QQ_RESIDUAL",
  "A_QM_RESIDUAL",
  "A_MM_RESIDUAL",
  "DECOUPLED_A_QM_NO_H_MASS",
  "MIXED_QL0_EP_RISK_TERM",
  "DECOUPLED_A_MM_DENSITY",
  "QM_FLIP_A_QM_ODD",
  "A_MM_EVEN_SIGNED_FLIP",
  "H_MASS_PROJECTION",
  "RADIATION_BARE_RATIO",
  "RADIATION_BARE_EXPONENT",
  "RADIATION_PINNED_EXPONENT",
  "RADIATION_WRONG_NORM_DISCRIMINATES",
  "STABILITY_BOUND",
  "CHANNEL_H_EP",
  "CHANNEL_RADIATION",
  "CHANNEL_UNIVERSALITY",
  "CHANNEL_UL_EP",
  "CHANNEL_PREFERRED_FRAME",
  "SUBTAGS",
  "EP_ADJUDICATION",
  "PRODUCTION_VERDICT",
  "GUARD_A_PRODUCTION",
  "GUARD_A_NEGATIVE_FIXTURES",
  "GUARD_A_DIRECT_INJECTION",
  "GUARD_A_FORGED_FLAG",
  "GUARD_A_BYPASS_REGRESSIONS",
  "GUARD_A_SCOPE_DENYLIST",
  "CTRL_B",
  "CTRL_C",
  "CTRL_D",
  "CTRL_E",
  "CTRL_F",
  "CTRL_G",
  "CTRL_H",
  "CTRL_I",
  "CTRL_J",
  "CTRL_K_WITHOUT_PINNED_KH",
  "CTRL_K",
  "DELETION_H_EP_PILLAR",
  "DELETION_GATED_INDIVIDUAL_STABLE",
  "DELETION_COLLECTIVE_NATURALLY_HIDDEN",
  "REACHABLE_FALSIFIERS",
  "LORENTZ_NECESSARY_NOT_SUFFICIENT",
  "DIMENSION_HOMOGENEITY",
  "DUAL_ENGINE_TERMS",
  "VERDICT_REDERIVATION",
  "SOURCE_TO_STAGE_MANIFEST"
};

ablationDescriptions = <|
  "STATIC_DETERMINANT" ->
    "add one to the typed determinant target",
  "A_QQ_RESIDUAL" ->
    "add one to the typed A_qq closed form",
  "A_QM_RESIDUAL" ->
    "add one to the typed A_qm closed form",
  "A_MM_RESIDUAL" ->
    "add one to the typed A_mm closed form",
  "DECOUPLED_A_QM_NO_H_MASS" ->
    "add one to the decoupled A_qm target",
  "MIXED_QL0_EP_RISK_TERM" ->
    "add one to the mixed q_L=0 target",
  "DECOUPLED_A_MM_DENSITY" ->
    "add one to the decoupled A_mm target",
  "QM_FLIP_A_QM_ODD" ->
    "perturb the q_M sign-flip substitution",
  "A_MM_EVEN_SIGNED_FLIP" ->
    "perturb the q_M flip in both A_mm checks",
  "H_MASS_PROJECTION" ->
    "corrupt the fixture h-mass projection target",
  "RADIATION_BARE_RATIO" ->
    "add one to the typed bare flux-ratio target",
  "RADIATION_BARE_EXPONENT" ->
    "insert an extra c_E in the bare ratio",
  "RADIATION_PINNED_EXPONENT" ->
    "use M_h=K_h/c_E instead of K_h/c_E^2",
  "RADIATION_WRONG_NORM_DISCRIMINATES" ->
    "replace the wrong-normalization branch by the bare branch",
  "STABILITY_BOUND" ->
    "corrupt the typed strict-stability slack",
  "CHANNEL_H_EP" ->
    "break the private static-Coulomb input",
  "CHANNEL_RADIATION" ->
    "select the corrupt-speed branch",
  "CHANNEL_UNIVERSALITY" ->
    "make b/ell species indexed",
  "CHANNEL_UL_EP" ->
    "escalate the private u_L channel to NO_GO",
  "CHANNEL_PREFERRED_FRAME" ->
    "require large-c_E suppression",
  "SUBTAGS" ->
    "remove nonzero C_hu from the private subtag ledger",
  "EP_ADJUDICATION" ->
    "report unqualified decoupled-floor EP safety",
  "PRODUCTION_VERDICT" ->
    "break static-Coulomb match in the private production-verdict witness",
  "GUARD_A_PRODUCTION" ->
    "inject a numeric P_h/P_EM into the production tree",
  "GUARD_A_NEGATIVE_FIXTURES" ->
    "authorize all five fixtures through a neutralized shared guard",
  "GUARD_A_DIRECT_INJECTION" ->
    "authorize the directly injected guarded numbers",
  "GUARD_A_FORGED_FLAG" ->
    "let a forged local earned_parent_action flag authorize emission",
  "GUARD_A_BYPASS_REGRESSIONS" ->
    "disable list/tuple descent in the numeric-leaf walk",
  "GUARD_A_SCOPE_DENYLIST" ->
    "add diagnostic_gain to the scoped denylist",
  "CTRL_B" ->
    "neutralize the h bare-mass-residue fixture",
  "CTRL_C" ->
    "neutralize C_hu in the mixed-risk fixture",
  "CTRL_D" ->
    "neutralize the corrupt-speed selector",
  "CTRL_E" ->
    "neutralize the large-c_E preferred-frame selector",
  "CTRL_F" ->
    "neutralize the species-indexed b/ell selector",
  "CTRL_G" ->
    "neutralize the strict-stability violation",
  "CTRL_H" ->
    "restore the static-Coulomb match",
  "CTRL_I" ->
    "neutralize the q_M sign flip",
  "CTRL_J" ->
    "neutralize the wrong-normalization selector",
  "CTRL_K_WITHOUT_PINNED_KH" ->
    "treat c_E=c_gamma alone as sufficient to escalate radiation",
  "CTRL_K" ->
    "drop the pinned-K_h fact",
  "DELETION_H_EP_PILLAR" ->
    "leave h_EP unstubbed",
  "DELETION_GATED_INDIVIDUAL_STABLE" ->
    "replace one individual earned stub by NO_GO",
  "DELETION_COLLECTIVE_NATURALLY_HIDDEN" ->
    "omit universality from the collective earned stub",
  "REACHABLE_FALSIFIERS" ->
    "drop control B from the computed active falsifier map",
  "LORENTZ_NECESSARY_NOT_SUFFICIENT" ->
    "treat c_E=c_gamma alone as an earned radiation state",
  "DIMENSION_HOMOGENEITY" ->
    "change [c_E] from L T^-1 to L^2 T^-1",
  "DUAL_ENGINE_TERMS" ->
    "drop the locally computed determinant inventory term",
  "VERDICT_REDERIVATION" ->
    "neutralize the computed control-B witness while retaining its named verdict",
  "SOURCE_TO_STAGE_MANIFEST" ->
    "mis-scope one preserved source predicate as replaced-by-stronger"
|>;

If[
  Length[toothOrder] =!= 49 ||
    Sort[Keys[ablationDescriptions]] =!= Sort[toothOrder],
  Print["FAIL: stage042 tooth/ablation declaration mismatch"];
  Exit[1]
];

raise[message_] := Throw[message, "ledgerStage042Failure"];

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

verdictMapped = "SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED";
verdictNoGo = "NO_GO_CONSISTENCY";
hEpSafe = "EARNED_SAFE_MASS_CHANNEL_ON_DECOUPLED_FLOOR";
hEpFifth = "FIFTH_FORCE_TRIGGERED";
radiationSim = "SIM_GATED";
radiationTension = "FALSIFIABLE_TENSION";
radiationAudit = "AUDIT_ONLY_NOT_EARNED";
universalitySim = "SIM_GATED_REQUIRED_UNIVERSALITY";
universalityTension = "FALSIFIABLE_TENSION";
ulSim = "SIM_GATED";
pfSim = "SIM_GATED";
pfTension = "PREFERRED_FRAME_TENSION";
earnedSafe = "EARNED_SAFE";


(* ---------------------------------------------------------------------- *)
(* Native symbolic core: LinearSolve, not Inverse.                         *)
(* ---------------------------------------------------------------------- *)

stiffness = {{Beff, Chu}, {Chu, Kh}};
determinant = Factor[Det[stiffness]];
chargeSource = {qL, qh};
massSource = {qM, mh};
chargeResponse = FullSimplify[LinearSolve[stiffness, chargeSource]];
massResponse = FullSimplify[LinearSolve[stiffness, massSource]];
Aqq = Factor[chargeSource.chargeResponse];
Aqm = Factor[chargeSource.massResponse];
Amm = Factor[massSource.massResponse];

typedAqq =
  (Beff*qh^2 - 2*Chu*qL*qh + Kh*qL^2)/determinant;
typedAqm =
  (Beff*mh*qh - Chu*mh*qL - Chu*qM*qh + Kh*qL*qM)/
    determinant;
typedAmm =
  (Beff*mh^2 - 2*Chu*mh*qM + Kh*qM^2)/determinant;

acceleration = omega^2*d;
scalarGradient = qh*acceleration/(Mh*cE^2*r);
scalarPower = Factor[Mh*cE*scalarGradient^2*r^2];
emGradient = QE*acceleration/(cg^2*r);
emPower = Factor[cg*emGradient^2*r^2];
ratioBare = Factor[scalarPower/emPower];
ratioPinned =
  Factor[(ratioBare /. Mh -> Kh/cE^2)*Kh*kappa/qh^2];
corruptGradient = qh*acceleration/(Mh*cE*r);
corruptPower = Factor[Mh*cE*corruptGradient^2*r^2];
ratioCorrupt = Factor[corruptPower/emPower];
ratioWrong =
  Factor[(ratioBare /. Mh -> Kh/cE^2)*Kh*kappa/qh^2];

positiveParameters =
  And @@ Thread[{cg, Mh, Kh, QE, qh, kappa} > 0];
positiveAssumptions = cE > 0 && positiveParameters;

speedFalloff[expr_] := FullSimplify[
  Limit[
    FullSimplify[
      -cE*D[Log[Together[expr]], cE],
      Assumptions -> positiveAssumptions
    ],
    cE -> Infinity,
    Assumptions -> positiveParameters
  ],
  Assumptions -> positiveAssumptions
];

bareExponent = speedFalloff[ratioBare];
pinnedExponent = speedFalloff[ratioPinned];
corruptExponent = speedFalloff[ratioCorrupt];
wrongExponent = speedFalloff[ratioWrong];
mixedGenerator = Factor[Aqm /. {qL -> 0, mh -> 0}];
mixedGeneratorNonzero = ! TrueQ[PossibleZeroQ[mixedGenerator]];


(* ---------------------------------------------------------------------- *)
(* Ordered rule-table channel adjudication and set-membership verdict.     *)
(* ---------------------------------------------------------------------- *)

productionFacts[] := <|
  "static_coulomb_match" -> True,
  "h_bare_mass_residue_zero" -> True,
  "h_bare_mass_residue_fixture" -> False,
  "pinned_kh_fact" -> False,
  "commit_cE_equals_cgamma" -> False,
  "radiation_derivation_ok" -> True,
  "radiation_corrupt_speed" -> False,
  "radiation_wrong_normalization" -> False,
  "q_ratio_global_earned" -> False,
  "species_indexed_b_over_ell" -> False,
  "uL_no_go" -> False,
  "uL_escalated" -> False,
  "suppression_requires_large_cE" -> False,
  "C_hu_nonzero" -> True,
  "qM_nonzero" -> True,
  "qM_sign" -> 1,
  "stability_violation" -> False,
  "production_laundering_attempt" -> False
|>;

makeFacts[overrides_: <||>] := Join[productionFacts[], overrides];

resolveOrdered[table_List] :=
  FirstCase[table, {True, state_} :> state, "UNRESOLVED"];

deriveChannels[facts_Association, commitAloneSufficient_: False] :=
 Module[
  {
    hRules, radiationRules, universalityRules, ulRules, pfRules,
    hState, radiationState, universalityState, ulState, pfState,
    mixedRisk
  },
  hRules = {
    {! TrueQ[facts["static_coulomb_match"]], "NO_GO"},
    {
      TrueQ[facts["h_bare_mass_residue_fixture"]] ||
        ! TrueQ[facts["h_bare_mass_residue_zero"]],
      hEpFifth
    },
    {True, hEpSafe}
  };
  radiationRules = {
    {TrueQ[facts["radiation_corrupt_speed"]], radiationAudit},
    {
      TrueQ[facts["pinned_kh_fact"]] &&
        TrueQ[facts["commit_cE_equals_cgamma"]],
      radiationTension
    },
    {
      TrueQ[commitAloneSufficient] &&
        TrueQ[facts["commit_cE_equals_cgamma"]],
      earnedSafe
    },
    {TrueQ[facts["radiation_derivation_ok"]], radiationSim},
    {True, radiationAudit}
  };
  universalityRules = {
    {
      TrueQ[facts["species_indexed_b_over_ell"]],
      universalityTension
    },
    {TrueQ[facts["q_ratio_global_earned"]], earnedSafe},
    {True, universalitySim}
  };
  ulRules = {
    {TrueQ[facts["uL_no_go"]], "NO_GO"},
    {TrueQ[facts["uL_escalated"]], universalityTension},
    {True, ulSim}
  };
  pfRules = {
    {TrueQ[facts["suppression_requires_large_cE"]], pfTension},
    {True, pfSim}
  };
  hState = resolveOrdered[hRules];
  radiationState = resolveOrdered[radiationRules];
  universalityState = resolveOrdered[universalityRules];
  ulState = resolveOrdered[ulRules];
  pfState = resolveOrdered[pfRules];
  mixedRisk =
    TrueQ[facts["C_hu_nonzero"]] &&
    TrueQ[facts["qM_nonzero"]] &&
    mixedGeneratorNonzero;
  <|
    "h_EP" -> hState,
    "radiation" -> radiationState,
    "universality" -> universalityState,
    "u_L_EP" -> ulState,
    "preferred_frame" -> pfState,
    "CHERENKOV_DEFERRED" -> (radiationState === radiationSim),
    "MIXED_SCALAR_EP_RISK" -> mixedRisk,
    "LAUNDERING_REFUSED" ->
      TrueQ[facts["production_laundering_attempt"]],
    "STABILITY_VIOLATED" -> TrueQ[facts["stability_violation"]]
  |>
 ];

channelNames = {
  "h_EP", "radiation", "universality", "u_L_EP",
  "preferred_frame"
};

channelStateAssociation[channels_Association] :=
  AssociationMap[channels[#] &, channelNames];

verdictFromChannels[channels_Association] := Module[
  {
    states, noGoQ, tensionTests, tensionChannels, gatedStates,
    gatedOrCostQ, allEarnedQ, verdictRules
  },
  states = Values[channelStateAssociation[channels]];
  noGoQ =
    MemberQ[states, "NO_GO"] ||
    TrueQ[channels["LAUNDERING_REFUSED"]] ||
    TrueQ[channels["STABILITY_VIOLATED"]];
  tensionTests = <|
    "h_EP" -> (channels["h_EP"] === hEpFifth),
    "radiation" -> (channels["radiation"] === radiationTension),
    "universality" ->
      (channels["universality"] === universalityTension),
    "u_L_EP" -> (channels["u_L_EP"] === universalityTension)
  |>;
  tensionChannels = Select[
    {"h_EP", "radiation", "universality", "u_L_EP"},
    TrueQ[tensionTests[#]] &
  ];
  gatedStates = {
    channels["radiation"],
    channels["universality"],
    channels["u_L_EP"],
    channels["preferred_frame"]
  };
  gatedOrCostQ =
    ! DisjointQ[
      gatedStates,
      {radiationSim, radiationAudit, universalitySim, ulSim,
       pfSim, pfTension}
    ];
  allEarnedQ =
    SubsetQ[{earnedSafe}, DeleteDuplicates[gatedStates]] &&
    SubsetQ[DeleteDuplicates[gatedStates], {earnedSafe}];
  verdictRules = {
    {noGoQ, verdictNoGo},
    {
      Length[tensionChannels] > 0,
      "FALSIFIABLE_TENSION(channel=" <>
        StringRiffle[tensionChannels, ","] <> ")"
    },
    {
      channels["h_EP"] === hEpSafe && gatedOrCostQ,
      verdictMapped
    },
    {
      channels["h_EP"] === hEpSafe && allEarnedQ,
      "NATURALLY_HIDDEN"
    },
    {True, verdictNoGo}
  };
  resolveOrdered[verdictRules]
];

epAdjudication[
  channels_Association,
  unqualified_: False
] := Module[{fullStatus},
  fullStatus = resolveOrdered[{
    {
      channels["h_EP"] === hEpSafe &&
        channels["universality"] === earnedSafe,
      earnedSafe
    },
    {True, "NOT_EARNED"}
  }];
  <|
    "mass_channel" -> channels["h_EP"],
    "full_decoupled_floor_EP_safety" -> fullStatus,
    "unqualified_decoupled_floor_EP_safe_reported" ->
      TrueQ[unqualified]
  |>
];


(* ---------------------------------------------------------------------- *)
(* Guard A: Position numeric/string leaves plus Cases key-path recovery.   *)
(* ---------------------------------------------------------------------- *)

normalizeToken[value_] := ToLowerCase@StringReplace[
  ToString[value, InputForm],
  {
    "/" -> "_over_",
    "-" -> "_",
    " " -> "_",
    "." -> "_",
    "\"" -> ""
  }
];

guardedFields = {
  "m_h", "c_e", "k_h", "p_h_over_p_em", "ph_over_pem",
  "power_ratio", "flux_ratio", "ep_magnitude",
  "fifth_force_magnitude", "residue_floor"
};

numericLiteralPattern =
  RegularExpression[
    "(^|[^A-Za-z0-9_])[-+]?(([0-9]+(\\.[0-9]*)?)|(\\.[0-9]+))([eE][-+]?[0-9]+)?([^A-Za-z0-9_]|$)"
  ];

positionKeyTokens[position_List] :=
  Cases[position, Key[key_] :> normalizeToken[key], Infinity];

guardedPositionQ[position_List, extraFields_List] := Module[
  {tokens, deny},
  tokens = positionKeyTokens[position];
  deny = Join[guardedFields, extraFields];
  AnyTrue[
    tokens,
    Function[token,
      MemberQ[deny, token] ||
        AnyTrue[
          {"m_h", "c_e", "k_h"},
          Function[prefix,
            token === prefix ||
              StringStartsQ[token, prefix <> "_"]
          ]
        ]
    ]
  ]
];

findGuardedNumericPositions[
  data_,
  descendSequences_: True,
  scanStrings_: True,
  extraFields_: {}
] := Module[
  {numericPositions, stringPositions, candidates},
  numericPositions =
    Position[data, leaf_ /; NumericQ[leaf], {0, Infinity},
      Heads -> False];
  If[
    ! TrueQ[descendSequences],
    numericPositions =
      Select[numericPositions, FreeQ[#, _Integer] &]
  ];
  stringPositions = If[
    TrueQ[scanStrings],
    Position[
      data,
      leaf_String /; StringContainsQ[leaf, numericLiteralPattern],
      {0, Infinity},
      Heads -> False
    ],
    {}
  ];
  candidates = DeleteDuplicates@Join[numericPositions, stringPositions];
  Select[candidates, guardedPositionQ[#, extraFields] &]
];

absentSubstrate[] := <|
  "facts" -> <|"h_time_kinetic_parent_action" -> "ABSENT"|>,
  "provenance_objects" ->
    <|"h_time_kinetic_parent_action" -> None|>
|>;

earnedSubstrate[] := <|
  "facts" ->
    <|"h_time_kinetic_parent_action" -> "DERIVED_PARENT_ACTION"|>,
  "provenance_objects" -> <|
    "h_time_kinetic_parent_action" -> <|
      "status" -> "EARNED",
      "source" -> "substrate:h_time_kinetic_parent_action"
    |>
  |>
|>;

parentActionEarned[substrate_Association] := Module[
  {fact, provenance},
  fact = Lookup[
    substrate["facts"],
    "h_time_kinetic_parent_action",
    Missing["NotAvailable"]
  ];
  provenance = Lookup[
    substrate["provenance_objects"],
    "h_time_kinetic_parent_action",
    None
  ];
  ! MemberQ[{None, "ABSENT", "MISSING", "NOT_EARNED"}, fact] &&
    AssociationQ[provenance] &&
    Lookup[provenance, "status", ""] === "EARNED" &&
    StringQ[Lookup[provenance, "source", None]] &&
    StringLength[Lookup[provenance, "source", ""]] > 0
];

guardScan[
  payload_Association,
  substrate_Association,
  allowForgedLocal_: False,
  descendSequences_: True,
  scanStrings_: True,
  extraFields_: {}
] := Module[
  {forged, authorized, paths, result},
  forged = TrueQ[Lookup[payload, "earned_parent_action", False]];
  authorized =
    parentActionEarned[substrate] ||
    (TrueQ[allowForgedLocal] && forged);
  paths = If[
    authorized,
    {},
    findGuardedNumericPositions[
      payload,
      descendSequences,
      scanStrings,
      extraFields
    ]
  ];
  result = resolveOrdered[{
    {authorized, "MAGNITUDE_EMISSION_ALLOWED"},
    {Length[paths] > 0, "LAUNDERING_REFUSED"},
    {True, "NO_GUARDED_MAGNITUDE_EMITTED"}
  }];
  <|"result" -> result, "paths" -> paths, "authorized" -> authorized|>
];

negativeFixtures = {
  {
    "M_h_from_N0",
    <|"assembled_results_fixture" -> <|"M_h" -> 1|>|>
  },
  {
    "M_h_from_K_parallel",
    <|"assembled_results_fixture" -> <|"M_h" -> 2|>|>
  },
  {
    "c_E_from_c_gamma_Green",
    <|"assembled_results_fixture" -> <|"c_E" -> 3|>|>
  },
  {
    "K_h_from_N0_cgamma2",
    <|"assembled_results_fixture" -> <|"K_h" -> 4|>|>
  },
  {
    "residue_floor_emission",
    <|"assembled_results_fixture" -> <|"residue_floor" -> 5|>|>
  }
};

guardNegativeVector[neutralized_: False] := Module[{substrate},
  substrate = If[TrueQ[neutralized], earnedSubstrate[], absentSubstrate[]];
  (guardScan[#[[2]], substrate]["result"] &) /@ negativeFixtures
];


(* ---------------------------------------------------------------------- *)
(* Data-driven controls and deletion sensitivity.                          *)
(* ---------------------------------------------------------------------- *)

controlResults[mutate_: ""] := Module[
  {
    ledgers, channels, result, sign, aqmBase, signedBase, ammBase,
    aqmCase, signedCase, ammCase, k0Channels
  },
  ledgers = <|
    "B" -> makeFacts[<|
      "h_bare_mass_residue_zero" -> (mutate === "CTRL_B"),
      "h_bare_mass_residue_fixture" -> (mutate =!= "CTRL_B"),
      "C_hu_nonzero" -> False
    |>],
    "C" -> makeFacts[<|
      "C_hu_nonzero" -> (mutate =!= "CTRL_C"),
      "qM_nonzero" -> True
    |>],
    "D" -> makeFacts[<|
      "radiation_corrupt_speed" -> (mutate =!= "CTRL_D")
    |>],
    "E" -> makeFacts[<|
      "suppression_requires_large_cE" -> (mutate =!= "CTRL_E")
    |>],
    "F" -> makeFacts[<|
      "species_indexed_b_over_ell" -> (mutate =!= "CTRL_F")
    |>],
    "G" -> makeFacts[<|
      "stability_violation" -> (mutate =!= "CTRL_G")
    |>],
    "H" -> makeFacts[<|
      "static_coulomb_match" -> (mutate === "CTRL_H")
    |>],
    "I" -> makeFacts[<|
      "qM_sign" -> If[mutate === "CTRL_I", 1, -1]
    |>],
    "J" -> makeFacts[<|
      "radiation_wrong_normalization" -> (mutate =!= "CTRL_J")
    |>],
    "K_without_pinned_Kh" -> makeFacts[<|
      "commit_cE_equals_cgamma" -> True
    |>],
    "K" -> makeFacts[<|
      "pinned_kh_fact" -> (mutate =!= "CTRL_K"),
      "commit_cE_equals_cgamma" -> True
    |>]
  |>;
  channels = AssociationMap[
    deriveChannels[ledgers[#]] &,
    Keys[ledgers]
  ];
  k0Channels = deriveChannels[
    ledgers["K_without_pinned_Kh"],
    mutate === "CTRL_K_WITHOUT_PINNED_KH"
  ];
  channels["K_without_pinned_Kh"] = k0Channels;
  result = AssociationMap[
    <|
      "channels" -> channels[#],
      "verdict" -> verdictFromChannels[channels[#]]
    |> &,
    Keys[channels]
  ];
  result["D"] = Join[
    result["D"],
    <|
      "baseline_exp" -> bareExponent,
      "selected_exp" ->
        If[mutate === "CTRL_D", bareExponent, corruptExponent]
    |>
  ];
  sign = If[mutate === "CTRL_I", 1, -1];
  aqmBase = Aqm /. mh -> 0;
  signedBase = Factor[(Amm /. mh -> 0)/qM];
  ammBase = Amm /. mh -> 0;
  aqmCase = aqmBase /. qM -> sign*qM;
  signedCase = signedBase /. qM -> sign*qM;
  ammCase = ammBase /. qM -> sign*qM;
  result["I"] = Join[
    result["I"],
    <|
      "aqm_flipped" -> TrueQ[FullSimplify[aqmCase + aqmBase == 0]],
      "signed_flipped" ->
        TrueQ[FullSimplify[signedCase + signedBase == 0]],
      "amm_even" -> TrueQ[FullSimplify[ammCase - ammBase == 0]]
    |>
  ];
  result["J"] = Join[
    result["J"],
    <|
      "baseline_exp" -> bareExponent,
      "selected_exp" ->
        If[mutate === "CTRL_J", bareExponent, wrongExponent]
    |>
  ];
  result["_baseline"] = <|
    "verdict" ->
      verdictFromChannels[deriveChannels[productionFacts[]]]
  |>;
  result
];

stubChannels[channels_Association, updates_Association] :=
  Join[channels, updates];

deletionOutcomes[production_Association, mutate_: ""] := Module[
  {
    baseline, hStub, hVerdict, individual, names, candidate,
    candidateVerdict, collectiveUpdates, collective,
    collectiveVerdict
  },
  baseline = verdictFromChannels[production];
  hStub = If[
    mutate === "DELETION_H_EP_PILLAR",
    production,
    stubChannels[production, <|"h_EP" -> earnedSafe|>]
  ];
  hVerdict = verdictFromChannels[hStub];
  names = {"radiation", "universality", "u_L_EP", "preferred_frame"};
  individual = Association@MapIndexed[
    Function[{name, index},
      candidate = stubChannels[
        production,
        <|
          name -> If[
            mutate === "DELETION_GATED_INDIVIDUAL_STABLE" &&
              First[index] === 1,
            "NO_GO",
            earnedSafe
          ]
        |>
      ];
      candidateVerdict = verdictFromChannels[candidate];
      name -> {candidateVerdict, candidateVerdict =!= baseline}
    ],
    names
  ];
  collectiveUpdates = <|
    "radiation" -> earnedSafe,
    "universality" -> earnedSafe,
    "u_L_EP" -> earnedSafe,
    "preferred_frame" -> earnedSafe
  |>;
  If[
    mutate === "DELETION_COLLECTIVE_NATURALLY_HIDDEN",
    collectiveUpdates = KeyDrop[collectiveUpdates, "universality"]
  ];
  collective = stubChannels[production, collectiveUpdates];
  collectiveVerdict = verdictFromChannels[collective];
  <|
    "h_EP" -> {hVerdict, hVerdict =!= baseline},
    "individual" -> individual,
    "collective" ->
      {collectiveVerdict, collectiveVerdict =!= baseline},
    "collective_channels" -> collective
  |>
];

reachableFalsifierMap[controlMap_Association] := Module[
  {baseline, names},
  baseline = controlMap["_baseline", "verdict"];
  names = {
    "B", "C", "D", "E", "F", "G", "H", "I", "J",
    "K_without_pinned_Kh", "K"
  };
  Association@Cases[
    names,
    name_ /; controlMap[name, "verdict"] =!= baseline :>
      (name -> controlMap[name, "verdict"])
  ]
];


(* ---------------------------------------------------------------------- *)
(* MLT dimension object.                                                   *)
(* ---------------------------------------------------------------------- *)

dimensionGuard[mutate_: False] := Module[
  {
    zeroDim, stiffnessDim, speedDim, lengthDim, frequencyDim,
    chargeDim, dims, detDim, detTerms, aqqTerms, aqmTerms,
    ammTerms, accelDim, scalarGradientDim, scalarPowerDim,
    emGradientDim, emPowerDim, ratioDim, quadraticTarget,
    physicalPower, homogeneous
  },
  zeroDim = {0, 0, 0};
  stiffnessDim = {1, 0, 0};
  speedDim = {0, 1, -1};
  lengthDim = {0, 1, 0};
  frequencyDim = {0, 0, -1};
  chargeDim = {1/2, 3/2, -1};
  dims = <|
    "B" -> stiffnessDim,
    "C" -> stiffnessDim,
    "K" -> stiffnessDim,
    "qL" -> chargeDim,
    "qh" -> chargeDim,
    "qM" -> chargeDim,
    "mh" -> chargeDim,
    "cE" -> {0, If[TrueQ[mutate], 2, 1], -1},
    "cg" -> speedDim,
    "Mh" -> zeroDim,
    "QE" -> chargeDim,
    "omega" -> frequencyDim,
    "d" -> lengthDim,
    "r" -> lengthDim
  |>;
  detTerms = {
    dims["B"] + dims["K"],
    2*dims["C"]
  };
  detDim = First[detTerms];
  aqqTerms = {
    dims["B"] + 2*dims["qh"],
    dims["C"] + dims["qL"] + dims["qh"],
    dims["K"] + 2*dims["qL"]
  };
  aqmTerms = {
    dims["B"] + dims["mh"] + dims["qh"],
    dims["C"] + dims["mh"] + dims["qL"],
    dims["C"] + dims["qM"] + dims["qh"],
    dims["K"] + dims["qL"] + dims["qM"]
  };
  ammTerms = {
    dims["B"] + 2*dims["mh"],
    dims["C"] + dims["mh"] + dims["qM"],
    dims["K"] + 2*dims["qM"]
  };
  accelDim = 2*dims["omega"] + dims["d"];
  scalarGradientDim =
    dims["qh"] + accelDim - dims["Mh"] -
      2*dims["cE"] - dims["r"];
  scalarPowerDim =
    dims["Mh"] + dims["cE"] + 2*scalarGradientDim +
      2*dims["r"];
  emGradientDim =
    dims["QE"] + accelDim - 2*dims["cg"] - dims["r"];
  emPowerDim =
    dims["cg"] + 2*emGradientDim + 2*dims["r"];
  ratioDim =
    2*dims["qh"] - dims["Mh"] - 2*dims["QE"] +
      3*dims["cg"] - 3*dims["cE"];
  quadraticTarget = 2*chargeDim - stiffnessDim;
  physicalPower = {1, 2, -3};
  homogeneous =
    SameQ @@ detTerms &&
    First[detTerms] === 2*stiffnessDim &&
    SameQ @@ aqqTerms &&
    SameQ @@ aqmTerms &&
    SameQ @@ ammTerms &&
    First[aqqTerms] - detDim === quadraticTarget &&
    First[aqmTerms] - detDim === quadraticTarget &&
    First[ammTerms] - detDim === quadraticTarget &&
    scalarPowerDim === physicalPower &&
    emPowerDim === physicalPower &&
    ratioDim === zeroDim;
  <|
    "homogeneous" -> homogeneous,
    "determinant_terms" -> detTerms,
    "Aqq_terms" -> aqqTerms,
    "Aqm_terms" -> aqmTerms,
    "Amm_terms" -> ammTerms,
    "scalar_power" -> scalarPowerDim,
    "em_power" -> emPowerDim,
    "ratio" -> ratioDim
  |>
];


(* ---------------------------------------------------------------------- *)
(* Bounded source-to-stage manifest.                                       *)
(* ---------------------------------------------------------------------- *)

sourceCoreIds = {
  "pathA42.py::core.A_qq_residual_zero",
  "pathA42.py::core.A_qm_residual_zero",
  "pathA42.py::core.A_mm_residual_zero",
  "pathA42.py::core.A_qm_decoupled_h_mass_zero",
  "pathA42.py::core.A_qm_mixed_qL0_nonzero_term",
  "pathA42.py::core.A_mm_decoupled_no_h_mass",
  "pathA42.py::core.mass_h_projection_zero",
  "pathA42.py::core.fixture_h_projection_nonzero_symbol",
  "pathA42.py::core.qM_flip_A_qm_residual_zero",
  "pathA42.py::core.A_mm_even_in_qM",
  "pathA42.py::core.A_mm_signed_projection_flips",
  "pathA42.py::core.stability_bound_strict",
  "pathA42.py::core.stability_violation_condition",
  "pathA42.py::core.radiation_bare_flux_ratio_matches",
  "pathA42.py::core.radiation_wrong_normalization_recomputed",
  "pathA42.py::core.radiation_bare_exponent",
  "pathA42.py::core.radiation_pinned_Kh_exponent",
  "pathA42.py::core.radiation_corrupt_speed_exponent",
  "pathA42.py::core.radiation_wrong_normalization_exponent"
};

sourceChannelIds = (
  "pathA42.py::channel." <> # &
) /@ {"h_EP", "radiation", "universality", "u_L_EP", "preferred_frame"};

sourceAdjudicationIds = {
  "pathA42.py::verdict_from_channels",
  "pathA42.py::ep_adjudication",
  "pathA42.py::guard_A_serialization_predicate"
};

sourceControlIds = (
  "pathA42.py::control." <> # &
) /@ {
  "A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
  "K_without_pinned_Kh", "K"
};

sourceDeletionIds = (
  "pathA42.py::deletion." <> # &
) /@ {
  "h_EP", "radiation", "universality", "u_L_EP",
  "preferred_frame", "gated_channels_collective"
};

replacedSourceIds = {
  "pathA42.wl::thin_mirror_coverage_gap"
};

scopedSourceIds = {
  "pathA42.py::argparse_compare_harness",
  "pathA42.py::file_writing_and_report_persistence",
  "pathA42.py::main_and_run_math_script",
  "pathA42.py::comparison_payload_digest_sha_leaf_compare",
  "pathA42.py::substrate_report_reads_and_assertContains_scans",
  "pathA42.py::substrate_and_report_narration_objects",
  "pathA42.wl::Export_comparison_payload"
};

preservedSourceIds = Join[
  sourceCoreIds,
  sourceChannelIds,
  sourceAdjudicationIds,
  sourceControlIds,
  sourceDeletionIds
];

sourcePredicateUniverse = Join[
  preservedSourceIds,
  replacedSourceIds,
  scopedSourceIds
];
sourcePredicateTotal = 53;

newStageIds = {
  "stage042.dual::REACHABLE_FALSIFIERS",
  "stage042.dual::LORENTZ_NECESSARY_NOT_SUFFICIENT",
  "stage042.dual::GUARD_A_SCOPE_DENYLIST",
  "stage042.dual::DIMENSION_HOMOGENEITY",
  "stage042.dual::DUAL_ENGINE_TERMS",
  "stage042.dual::SOURCE_TO_STAGE_MANIFEST"
};

If[
  Length[sourcePredicateUniverse] =!= sourcePredicateTotal,
  Print["FAIL: stage042 source predicate universe count is not 53"];
  Exit[1]
];

ratifiedCategory[identifier_String] := FirstCase[
  {
    {MemberQ[preservedSourceIds, identifier], "preserved-folded"},
    {MemberQ[replacedSourceIds, identifier], "replaced-by-stronger"},
    {MemberQ[scopedSourceIds, identifier], "scoped-out"},
    {MemberQ[newStageIds, identifier], "newly-added"}
  },
  {True, category_} :> category,
  "UNCLASSIFIED"
];

ratifiedManifestMap = AssociationMap[
  ratifiedCategory,
  Join[sourcePredicateUniverse, newStageIds]
];

computedManifestRows[mutate_: False] := Module[
  {moved},
  moved = First[preservedSourceIds];
  (
    # -> If[
      TrueQ[mutate] && # === moved,
      "replaced-by-stronger",
      ratifiedCategory[#]
    ] &
  ) /@ Join[sourcePredicateUniverse, newStageIds]
];

manifestGuards[mutate_: False] := Module[
  {
    rows, identifiers, sourceIdentifiers, sourceCounts, categories,
    categoryNames, disjoint, coverage, categoryMap
  },
  rows = computedManifestRows[mutate];
  identifiers = First /@ rows;
  sourceIdentifiers =
    Select[identifiers, MemberQ[sourcePredicateUniverse, #] &];
  sourceCounts = Counts[sourceIdentifiers];
  categories = AssociationMap[
    Function[category,
      Cases[rows, (identifier_ -> category) :> identifier]
    ],
    {
      "preserved-folded", "replaced-by-stronger", "newly-added",
      "scoped-out"
    }
  ];
  categoryNames = Keys[categories];
  disjoint = And @@ Flatten@Table[
    Intersection[
      categories[categoryNames[[left]]],
      categories[categoryNames[[right]]]
    ] === {},
    {left, 1, Length[categoryNames]},
    {right, left + 1, Length[categoryNames]}
  ];
  coverage = <|
    "ok" -> (
      Length[sourceIdentifiers] === sourcePredicateTotal &&
      Sort[sourceIdentifiers] === Sort[sourcePredicateUniverse] &&
      AllTrue[Values[sourceCounts], # === 1 &] &&
      Sort[identifiers] === Sort[Keys[ratifiedManifestMap]] &&
      disjoint
    ),
    "source_count" -> Length[sourceIdentifiers],
    "disjoint" -> disjoint,
    "partition" -> Counts[Last /@ rows]
  |>;
  categoryMap = Association[rows];
  {
    coverage,
    <|
      "ok" -> (categoryMap === ratifiedManifestMap),
      "computed" -> categoryMap,
      "ratified" -> ratifiedManifestMap
    |>
  }
];


(* ---------------------------------------------------------------------- *)
(* Local computed-term inventory and first-match witness table.            *)
(* ---------------------------------------------------------------------- *)

canonicalInventory[
  productionChannels_Association,
  productionVerdict_String,
  negativeGuardVector_List,
  directGuard_String,
  forgedGuard_List,
  bypassGuard_List,
  falsifiers_Association,
  deletion_Association,
  mutate_: False
] := Module[{terms},
  terms = {
    "verdict" -> productionVerdict,
    "channels" -> Normal[channelStateAssociation[productionChannels]],
    "det" -> Factor[determinant],
    "Aqq" -> Factor[Aqq],
    "Aqm" -> Factor[Aqm],
    "Amm" -> Factor[Amm],
    "exponents" ->
      {bareExponent, pinnedExponent, corruptExponent, wrongExponent},
    "guard_negative" -> negativeGuardVector,
    "guard_direct" -> directGuard,
    "guard_forged" -> forgedGuard,
    "guard_bypass" -> bypassGuard,
    "reachable" -> Sort[Normal[falsifiers]],
    "deletion_h" -> deletion["h_EP"],
    "deletion_individual" -> Normal[deletion["individual"]],
    "deletion_collective" -> deletion["collective"]
  };
  If[TrueQ[mutate], terms = DeleteCases[terms, "det" -> _]];
  terms
];

expectedInventory[] := {
  "verdict" -> verdictMapped,
  "channels" -> {
    "h_EP" -> hEpSafe,
    "radiation" -> radiationSim,
    "universality" -> universalitySim,
    "u_L_EP" -> ulSim,
    "preferred_frame" -> pfSim
  },
  "det" -> (Beff*Kh - Chu^2),
  "Aqq" -> Factor[typedAqq],
  "Aqm" -> Factor[typedAqm],
  "Amm" -> Factor[typedAmm],
  "exponents" -> {3, 1, 1, 1},
  "guard_negative" -> ConstantArray["LAUNDERING_REFUSED", 5],
  "guard_direct" -> "LAUNDERING_REFUSED",
  "guard_forged" ->
    {"LAUNDERING_REFUSED", "MAGNITUDE_EMISSION_ALLOWED"},
  "guard_bypass" ->
    {"LAUNDERING_REFUSED", "LAUNDERING_REFUSED"},
  "reachable" -> Sort[Normal@<|
    "B" -> "FALSIFIABLE_TENSION(channel=h_EP)",
    "F" -> "FALSIFIABLE_TENSION(channel=universality)",
    "G" -> verdictNoGo,
    "H" -> verdictNoGo,
    "K" -> "FALSIFIABLE_TENSION(channel=radiation)"
  |>],
  "deletion_h" -> {verdictNoGo, True},
  "deletion_individual" -> {
    "radiation" -> {verdictMapped, False},
    "universality" -> {verdictMapped, False},
    "u_L_EP" -> {verdictMapped, False},
    "preferred_frame" -> {verdictMapped, False}
  },
  "deletion_collective" -> {"NATURALLY_HIDDEN", True}
};

symbolicInventoryValueQ[value_] :=
  ! FreeQ[value, Beff | Kh | Chu | qL | qh | qM | mh];

inventoryValueEqual[left_, right_] := If[
  symbolicInventoryValueQ[left] || symbolicInventoryValueQ[right],
  TrueQ[FullSimplify[left - right == 0]],
  SameQ[left, right]
];

inventoriesEqual[actual_List, expected_List] :=
  Length[actual] === Length[expected] &&
  (First /@ actual) === (First /@ expected) &&
  And @@ MapThread[
    inventoryValueEqual[Last[#1], Last[#2]] &,
    {actual, expected}
  ];

verdictWitnesses[mutate_: False] := Module[
  {
    production, bFacts, b, f, k, g, h, bg, bkf,
    allEarned, natural
  },
  production =
    verdictFromChannels[deriveChannels[productionFacts[]]];
  bFacts = makeFacts[If[
    TrueQ[mutate],
    <|
      "h_bare_mass_residue_zero" -> True,
      "h_bare_mass_residue_fixture" -> False
    |>,
    <|
      "h_bare_mass_residue_zero" -> False,
      "h_bare_mass_residue_fixture" -> True
    |>
  ]];
  b = verdictFromChannels[deriveChannels[bFacts]];
  f = verdictFromChannels[deriveChannels[makeFacts[
    <|"species_indexed_b_over_ell" -> True|>
  ]]];
  k = verdictFromChannels[deriveChannels[makeFacts[<|
    "pinned_kh_fact" -> True,
    "commit_cE_equals_cgamma" -> True
  |>]]];
  g = verdictFromChannels[deriveChannels[makeFacts[
    <|"stability_violation" -> True|>
  ]]];
  h = verdictFromChannels[deriveChannels[makeFacts[
    <|"static_coulomb_match" -> False|>
  ]]];
  bg = verdictFromChannels[deriveChannels[makeFacts[<|
    "h_bare_mass_residue_zero" -> False,
    "h_bare_mass_residue_fixture" -> True,
    "stability_violation" -> True
  |>]]];
  bkf = verdictFromChannels[deriveChannels[makeFacts[<|
    "h_bare_mass_residue_zero" -> False,
    "h_bare_mass_residue_fixture" -> True,
    "pinned_kh_fact" -> True,
    "commit_cE_equals_cgamma" -> True,
    "species_indexed_b_over_ell" -> True
  |>]]];
  allEarned = stubChannels[
    deriveChannels[productionFacts[]],
    <|
      "radiation" -> earnedSafe,
      "universality" -> earnedSafe,
      "u_L_EP" -> earnedSafe,
      "preferred_frame" -> earnedSafe
    |>
  ];
  natural = verdictFromChannels[allEarned];
  {production, b, f, k, g, h, bg, bkf, natural}
];


runAssertions[] := Module[
  {
    productionChannels, productionVerdict, typedDet, target,
    aqmDecoupled, aqmMixed, ammDecoupled, aqmMass, oddFlip,
    oddResidual, ammMass, signedAmm, evenFlip, evenResidual,
    signedResidual, massProjection, fixtureProjection,
    fixtureTarget, targetBare, bareRatioCase, bareExpCase,
    pinnedRatioCase, pinnedExpCase, wrongRatioCase, wrongExpCase,
    discriminator, stabilityTarget, stableSample, violatedSample,
    stabilityActual, hCase, radiationCase, universalityCase,
    ulCase, pfCase, subtagCase, subtagActual, epCase, epActual,
    verdictFacts, verdictCase, productionTree, productionGuard,
    negativeVector, directSubstrate, directGuard, forgedPayload,
    forgedNegative, forgedPositive, forgedActual, tupleGuard,
    stringGuard, bypassActual, scopeExtra, scopeGuard,
    normalControls, controlMutate, controlCase, bActual, cActual,
    dActual, eActual, fActual, gActual, hActual, iActual, jActual,
    k0Actual, kActual, normalDeletion, deletionMutate,
    deletionCase, expectedIndividual, falsifiers,
    expectedFalsifiers, mappedControls, reachableActual,
    reachableExpected, commitOnly, lorentzActual, dims,
    forgedInventory, bypassInventory, inventoryActual,
    inventoryExpected, witnessActual, witnessExpected,
    manifestPair, coverageGuard, categoryGuard
  },
  productionChannels = deriveChannels[productionFacts[]];
  productionVerdict = verdictFromChannels[productionChannels];

  section["Static exchange and source projections"];
  typedDet = Beff*Kh - Chu^2 +
    If[activeMutation === "STATIC_DETERMINANT", 1, 0];
  expectBool[
    "STATIC_DETERMINANT",
    TrueQ[FullSimplify[determinant - typedDet == 0]],
    Factor[determinant - typedDet]
  ];

  target = typedAqq +
    If[activeMutation === "A_QQ_RESIDUAL", 1, 0];
  expectBool[
    "A_QQ_RESIDUAL",
    TrueQ[FullSimplify[Aqq - target == 0]],
    Factor[Aqq - target]
  ];
  target = typedAqm +
    If[activeMutation === "A_QM_RESIDUAL", 1, 0];
  expectBool[
    "A_QM_RESIDUAL",
    TrueQ[FullSimplify[Aqm - target == 0]],
    Factor[Aqm - target]
  ];
  target = typedAmm +
    If[activeMutation === "A_MM_RESIDUAL", 1, 0];
  expectBool[
    "A_MM_RESIDUAL",
    TrueQ[FullSimplify[Amm - target == 0]],
    Factor[Amm - target]
  ];

  aqmDecoupled = Factor[Aqm /. {Chu -> 0, mh -> 0}];
  target = qL*qM/Beff +
    If[activeMutation === "DECOUPLED_A_QM_NO_H_MASS", 1, 0];
  expectBool[
    "DECOUPLED_A_QM_NO_H_MASS",
    TrueQ[FullSimplify[aqmDecoupled - target == 0]],
    Factor[aqmDecoupled - target]
  ];

  aqmMixed = Factor[Aqm /. {qL -> 0, mh -> 0}];
  target = -Chu*qh*qM/determinant +
    If[activeMutation === "MIXED_QL0_EP_RISK_TERM", 1, 0];
  expectBool[
    "MIXED_QL0_EP_RISK_TERM",
    TrueQ[FullSimplify[aqmMixed - target == 0]],
    Factor[aqmMixed - target]
  ];

  ammDecoupled = Factor[Amm /. {Chu -> 0, mh -> 0}];
  target = qM^2/Beff +
    If[activeMutation === "DECOUPLED_A_MM_DENSITY", 1, 0];
  expectBool[
    "DECOUPLED_A_MM_DENSITY",
    TrueQ[FullSimplify[ammDecoupled - target == 0]],
    Factor[ammDecoupled - target]
  ];

  aqmMass = Aqm /. mh -> 0;
  oddFlip = If[
    activeMutation === "QM_FLIP_A_QM_ODD",
    -qM + 1,
    -qM
  ];
  oddResidual = Factor[(aqmMass /. qM -> oddFlip) + aqmMass];
  expectBool[
    "QM_FLIP_A_QM_ODD",
    TrueQ[FullSimplify[oddResidual == 0]],
    oddResidual
  ];

  ammMass = Amm /. mh -> 0;
  signedAmm = Factor[ammMass/qM];
  evenFlip = If[
    activeMutation === "A_MM_EVEN_SIGNED_FLIP",
    -qM + 1,
    -qM
  ];
  evenResidual = Factor[(ammMass /. qM -> evenFlip) - ammMass];
  signedResidual =
    Factor[(signedAmm /. qM -> evenFlip) + signedAmm];
  expectBool[
    "A_MM_EVEN_SIGNED_FLIP",
    TrueQ[FullSimplify[evenResidual == 0]] &&
      TrueQ[FullSimplify[signedResidual == 0]],
    <|"even" -> evenResidual, "signed" -> signedResidual|>
  ];

  massProjection = {0, 1}.{qM, 0};
  fixtureProjection = {0, 1}.{qM, mh};
  fixtureTarget = mh +
    If[activeMutation === "H_MASS_PROJECTION", 1, 0];
  expectBool[
    "H_MASS_PROJECTION",
    TrueQ[FullSimplify[massProjection == 0]] &&
      TrueQ[FullSimplify[fixtureProjection - fixtureTarget == 0]],
    <|
      "production" -> massProjection,
      "fixture_residual" -> fixtureProjection - fixtureTarget
    |>
  ];

  section["Radiation speed scaling and strict stability"];
  targetBare =
    qh^2/(Mh*QE^2)*(cg/cE)^3 +
    If[activeMutation === "RADIATION_BARE_RATIO", 1, 0];
  expectBool[
    "RADIATION_BARE_RATIO",
    TrueQ[FullSimplify[ratioBare - targetBare == 0]],
    Factor[ratioBare - targetBare]
  ];

  bareRatioCase = ratioBare*
    If[activeMutation === "RADIATION_BARE_EXPONENT", cE, 1];
  bareExpCase = speedFalloff[bareRatioCase];
  expectBool[
    "RADIATION_BARE_EXPONENT",
    bareExpCase === 3,
    bareExpCase
  ];

  pinnedRatioCase = If[
    activeMutation === "RADIATION_PINNED_EXPONENT",
    Factor[(ratioBare /. Mh -> Kh/cE)*Kh/qh^2],
    ratioPinned
  ];
  pinnedExpCase = speedFalloff[pinnedRatioCase];
  expectBool[
    "RADIATION_PINNED_EXPONENT",
    pinnedExpCase === 1,
    pinnedExpCase
  ];

  wrongRatioCase = If[
    activeMutation === "RADIATION_WRONG_NORM_DISCRIMINATES",
    ratioBare,
    ratioWrong
  ];
  wrongExpCase = speedFalloff[wrongRatioCase];
  discriminator = {
    wrongExpCase,
    bareExponent,
    wrongExpCase =!= bareExponent
  };
  expectBool[
    "RADIATION_WRONG_NORM_DISCRIMINATES",
    discriminator === {1, 3, True},
    discriminator
  ];

  stabilityTarget = Beff*Kh - Chu^2 +
    If[activeMutation === "STABILITY_BOUND", 1, 0];
  stableSample =
    determinant /. {Beff -> 2, Kh -> 3, Chu -> 1};
  violatedSample =
    determinant /. {Beff -> 2, Kh -> 2, Chu -> 2};
  stabilityActual = {
    TrueQ[FullSimplify[determinant - stabilityTarget == 0]],
    TrueQ[stableSample > 0],
    TrueQ[violatedSample <= 0]
  };
  expectBool[
    "STABILITY_BOUND",
    stabilityActual === {True, True, True},
    stabilityActual
  ];

  section["Five computed channels, subtags, and production verdict"];
  hCase = deriveChannels[makeFacts[<|
    "static_coulomb_match" ->
      (activeMutation =!= "CHANNEL_H_EP")
  |>]];
  expectBool[
    "CHANNEL_H_EP",
    hCase["h_EP"] === hEpSafe,
    hCase["h_EP"]
  ];

  radiationCase = deriveChannels[makeFacts[<|
    "radiation_corrupt_speed" ->
      (activeMutation === "CHANNEL_RADIATION")
  |>]];
  expectBool[
    "CHANNEL_RADIATION",
    radiationCase["radiation"] === radiationSim,
    radiationCase["radiation"]
  ];

  universalityCase = deriveChannels[makeFacts[<|
    "species_indexed_b_over_ell" ->
      (activeMutation === "CHANNEL_UNIVERSALITY")
  |>]];
  expectBool[
    "CHANNEL_UNIVERSALITY",
    universalityCase["universality"] === universalitySim,
    universalityCase["universality"]
  ];

  ulCase = deriveChannels[makeFacts[<|
    "uL_no_go" -> (activeMutation === "CHANNEL_UL_EP")
  |>]];
  expectBool[
    "CHANNEL_UL_EP",
    ulCase["u_L_EP"] === ulSim,
    ulCase["u_L_EP"]
  ];

  pfCase = deriveChannels[makeFacts[<|
    "suppression_requires_large_cE" ->
      (activeMutation === "CHANNEL_PREFERRED_FRAME")
  |>]];
  expectBool[
    "CHANNEL_PREFERRED_FRAME",
    pfCase["preferred_frame"] === pfSim,
    pfCase["preferred_frame"]
  ];

  subtagCase = deriveChannels[makeFacts[<|
    "C_hu_nonzero" -> (activeMutation =!= "SUBTAGS")
  |>]];
  subtagActual = {
    subtagCase["CHERENKOV_DEFERRED"],
    subtagCase["MIXED_SCALAR_EP_RISK"]
  };
  expectBool[
    "SUBTAGS",
    subtagActual === {True, True},
    subtagActual
  ];

  epCase = epAdjudication[
    productionChannels,
    activeMutation === "EP_ADJUDICATION"
  ];
  epActual = {
    epCase["mass_channel"],
    epCase["full_decoupled_floor_EP_safety"],
    epCase["unqualified_decoupled_floor_EP_safe_reported"]
  };
  expectBool[
    "EP_ADJUDICATION",
    epActual === {hEpSafe, "NOT_EARNED", False},
    epActual
  ];

  verdictFacts = makeFacts[If[
    activeMutation === "PRODUCTION_VERDICT",
    <|"static_coulomb_match" -> False|>,
    <||>
  ]];
  verdictCase = verdictFromChannels[deriveChannels[verdictFacts]];
  expectBool[
    "PRODUCTION_VERDICT",
    verdictCase === verdictMapped,
    verdictCase
  ];

  section["Guard A scoped laundering-refusal denylist"];
  productionTree = <|
    "radiation" -> <|
      "magnitude" ->
        "SIM_GATED_BY_GUARD_A_NO_NUMERIC_POWER_RATIO_EMITTED",
      "attempted_guarded_magnitude_emission" -> False
    |>
  |>;
  If[
    activeMutation === "GUARD_A_PRODUCTION",
    productionTree["radiation", "P_h/P_EM"] = 1
  ];
  productionGuard = guardScan[productionTree, absentSubstrate[]];
  expectBool[
    "GUARD_A_PRODUCTION",
    productionGuard["result"] ===
      "NO_GUARDED_MAGNITUDE_EMITTED",
    productionGuard
  ];

  negativeVector = guardNegativeVector[
    activeMutation === "GUARD_A_NEGATIVE_FIXTURES"
  ];
  expectBool[
    "GUARD_A_NEGATIVE_FIXTURES",
    negativeVector === ConstantArray["LAUNDERING_REFUSED", 5],
    negativeVector
  ];

  directSubstrate = If[
    activeMutation === "GUARD_A_DIRECT_INJECTION",
    earnedSubstrate[],
    absentSubstrate[]
  ];
  directGuard = guardScan[
    <|"P_h/P_EM" -> 1, "M_h" -> 1|>,
    directSubstrate
  ]["result"];
  expectBool[
    "GUARD_A_DIRECT_INJECTION",
    directGuard === "LAUNDERING_REFUSED",
    directGuard
  ];

  forgedPayload = <|
    "earned_parent_action" -> True,
    "P_h_over_P_EM" -> 1
  |>;
  forgedNegative = guardScan[
    forgedPayload,
    absentSubstrate[],
    activeMutation === "GUARD_A_FORGED_FLAG"
  ]["result"];
  forgedPositive = guardScan[
    <|"P_h_over_P_EM" -> 1|>,
    earnedSubstrate[]
  ]["result"];
  forgedActual = {forgedNegative, forgedPositive};
  expectBool[
    "GUARD_A_FORGED_FLAG",
    forgedActual ===
      {"LAUNDERING_REFUSED", "MAGNITUDE_EMISSION_ALLOWED"},
    forgedActual
  ];

  tupleGuard = guardScan[
    <|"assembled" -> <|"power_ratio" -> {0.137}|>|>,
    absentSubstrate[],
    False,
    activeMutation =!= "GUARD_A_BYPASS_REGRESSIONS"
  ]["result"];
  stringGuard = guardScan[
    <|"assembled" -> <|"power_ratio" -> "P_h/P_EM = 0.137"|>|>,
    absentSubstrate[]
  ]["result"];
  bypassActual = {tupleGuard, stringGuard};
  expectBool[
    "GUARD_A_BYPASS_REGRESSIONS",
    bypassActual ===
      {"LAUNDERING_REFUSED", "LAUNDERING_REFUSED"},
    bypassActual
  ];

  scopeExtra = If[
    activeMutation === "GUARD_A_SCOPE_DENYLIST",
    {"diagnostic_gain"},
    {}
  ];
  scopeGuard = guardScan[
    <|"diagnostic_gain" -> 7|>,
    absentSubstrate[],
    False,
    True,
    True,
    scopeExtra
  ]["result"];
  expectBool[
    "GUARD_A_SCOPE_DENYLIST",
    scopeGuard === "NO_GUARDED_MAGNITUDE_EMITTED",
    scopeGuard
  ];

  section["Twelve source controls (A is Guard A; B through K recomputed)"];
  normalControls = controlResults[];
  controlMutate = If[
    MemberQ[
      {
        "CTRL_B", "CTRL_C", "CTRL_D", "CTRL_E", "CTRL_F",
        "CTRL_G", "CTRL_H", "CTRL_I", "CTRL_J",
        "CTRL_K_WITHOUT_PINNED_KH", "CTRL_K"
      },
      activeMutation
    ],
    activeMutation,
    ""
  ];
  controlCase = controlResults[controlMutate];

  bActual = {
    controlCase["B", "channels", "h_EP"],
    controlCase["B", "verdict"]
  };
  expectBool[
    "CTRL_B",
    bActual === {
      hEpFifth,
      "FALSIFIABLE_TENSION(channel=h_EP)"
    },
    bActual
  ];

  cActual = {
    controlCase["C", "channels", "MIXED_SCALAR_EP_RISK"],
    controlCase["C", "verdict"]
  };
  expectBool[
    "CTRL_C",
    cActual === {True, verdictMapped},
    cActual
  ];

  dActual = {
    controlCase["D", "baseline_exp"],
    controlCase["D", "selected_exp"],
    controlCase["D", "channels", "radiation"],
    controlCase["D", "verdict"]
  };
  expectBool[
    "CTRL_D",
    dActual === {3, 1, radiationAudit, verdictMapped},
    dActual
  ];

  eActual = {
    controlCase["E", "channels", "preferred_frame"],
    controlCase["E", "verdict"]
  };
  expectBool[
    "CTRL_E",
    eActual === {pfTension, verdictMapped},
    eActual
  ];

  fActual = {
    controlCase["F", "channels", "universality"],
    controlCase["F", "verdict"]
  };
  expectBool[
    "CTRL_F",
    fActual === {
      universalityTension,
      "FALSIFIABLE_TENSION(channel=universality)"
    },
    fActual
  ];

  gActual = {
    controlCase["G", "channels", "STABILITY_VIOLATED"],
    controlCase["G", "verdict"]
  };
  expectBool[
    "CTRL_G",
    gActual === {True, verdictNoGo},
    gActual
  ];

  hActual = {
    controlCase["H", "channels", "h_EP"],
    controlCase["H", "verdict"]
  };
  expectBool[
    "CTRL_H",
    hActual === {"NO_GO", verdictNoGo},
    hActual
  ];

  iActual = {
    controlCase["I", "aqm_flipped"],
    controlCase["I", "signed_flipped"],
    controlCase["I", "amm_even"],
    controlCase["I", "verdict"]
  };
  expectBool[
    "CTRL_I",
    iActual === {True, True, True, verdictMapped},
    iActual
  ];

  jActual = {
    controlCase["J", "baseline_exp"],
    controlCase["J", "selected_exp"],
    controlCase["J", "verdict"]
  };
  expectBool[
    "CTRL_J",
    jActual === {3, 1, verdictMapped},
    jActual
  ];

  k0Actual = {
    controlCase["K_without_pinned_Kh", "channels", "radiation"],
    controlCase["K_without_pinned_Kh", "verdict"]
  };
  expectBool[
    "CTRL_K_WITHOUT_PINNED_KH",
    k0Actual === {radiationSim, verdictMapped},
    k0Actual
  ];

  kActual = {
    controlCase["K", "channels", "radiation"],
    controlCase["K", "verdict"]
  };
  expectBool[
    "CTRL_K",
    kActual === {
      radiationTension,
      "FALSIFIABLE_TENSION(channel=radiation)"
    },
    kActual
  ];

  section["Deletion sensitivity and reachable falsifiers"];
  normalDeletion = deletionOutcomes[productionChannels];
  deletionMutate = If[
    MemberQ[
      {
        "DELETION_H_EP_PILLAR",
        "DELETION_GATED_INDIVIDUAL_STABLE",
        "DELETION_COLLECTIVE_NATURALLY_HIDDEN"
      },
      activeMutation
    ],
    activeMutation,
    ""
  ];
  deletionCase = deletionOutcomes[productionChannels, deletionMutate];
  expectBool[
    "DELETION_H_EP_PILLAR",
    deletionCase["h_EP"] === {verdictNoGo, True},
    deletionCase["h_EP"]
  ];

  expectedIndividual = <|
    "radiation" -> {verdictMapped, False},
    "universality" -> {verdictMapped, False},
    "u_L_EP" -> {verdictMapped, False},
    "preferred_frame" -> {verdictMapped, False}
  |>;
  expectBool[
    "DELETION_GATED_INDIVIDUAL_STABLE",
    deletionCase["individual"] === expectedIndividual,
    deletionCase["individual"]
  ];
  expectBool[
    "DELETION_COLLECTIVE_NATURALLY_HIDDEN",
    deletionCase["collective"] === {"NATURALLY_HIDDEN", True},
    deletionCase["collective"]
  ];

  falsifiers = reachableFalsifierMap[normalControls];
  If[
    activeMutation === "REACHABLE_FALSIFIERS",
    falsifiers = KeyDrop[falsifiers, "B"]
  ];
  expectedFalsifiers = <|
    "B" -> "FALSIFIABLE_TENSION(channel=h_EP)",
    "F" -> "FALSIFIABLE_TENSION(channel=universality)",
    "G" -> verdictNoGo,
    "H" -> verdictNoGo,
    "K" -> "FALSIFIABLE_TENSION(channel=radiation)"
  |>;
  mappedControls = AssociationMap[
    normalControls[#, "verdict"] &,
    {"C", "D", "E", "I", "J", "K_without_pinned_Kh"}
  ];
  reachableActual = {falsifiers, mappedControls};
  reachableExpected = {
    expectedFalsifiers,
    AssociationMap[
      verdictMapped &,
      {"C", "D", "E", "I", "J", "K_without_pinned_Kh"}
    ]
  };
  expectBool[
    "REACHABLE_FALSIFIERS",
    reachableActual === reachableExpected,
    reachableActual
  ];

  commitOnly = deriveChannels[
    makeFacts[<|"commit_cE_equals_cgamma" -> True|>],
    activeMutation === "LORENTZ_NECESSARY_NOT_SUFFICIENT"
  ];
  lorentzActual = {
    commitOnly["radiation"],
    commitOnly["preferred_frame"],
    verdictFromChannels[commitOnly]
  };
  expectBool[
    "LORENTZ_NECESSARY_NOT_SUFFICIENT",
    lorentzActual === {radiationSim, pfSim, verdictMapped},
    lorentzActual
  ];

  section["Build-global guards and first-match witness table"];
  dims = dimensionGuard[
    activeMutation === "DIMENSION_HOMOGENEITY"
  ];
  expectBool[
    "DIMENSION_HOMOGENEITY",
    TrueQ[dims["homogeneous"]],
    dims
  ];

  forgedInventory = {
    guardScan[forgedPayload, absentSubstrate[]]["result"],
    guardScan[
      <|"P_h_over_P_EM" -> 1|>,
      earnedSubstrate[]
    ]["result"]
  };
  bypassInventory = {
    guardScan[
      <|"assembled" -> <|"power_ratio" -> {0.137}|>|>,
      absentSubstrate[]
    ]["result"],
    guardScan[
      <|
        "assembled" ->
          <|"power_ratio" -> "P_h/P_EM = 0.137"|>
      |>,
      absentSubstrate[]
    ]["result"]
  };
  inventoryActual = canonicalInventory[
    productionChannels,
    productionVerdict,
    guardNegativeVector[],
    guardScan[
      <|"P_h/P_EM" -> 1, "M_h" -> 1|>,
      absentSubstrate[]
    ]["result"],
    forgedInventory,
    bypassInventory,
    reachableFalsifierMap[normalControls],
    normalDeletion,
    activeMutation === "DUAL_ENGINE_TERMS"
  ];
  inventoryExpected = expectedInventory[];
  expectBool[
    "DUAL_ENGINE_TERMS",
    inventoriesEqual[inventoryActual, inventoryExpected],
    <|
      "actual_names" -> (First /@ inventoryActual),
      "expected_names" -> (First /@ inventoryExpected)
    |>
  ];

  witnessActual = verdictWitnesses[
    activeMutation === "VERDICT_REDERIVATION"
  ];
  witnessExpected = {
    verdictMapped,
    "FALSIFIABLE_TENSION(channel=h_EP)",
    "FALSIFIABLE_TENSION(channel=universality)",
    "FALSIFIABLE_TENSION(channel=radiation)",
    verdictNoGo,
    verdictNoGo,
    verdictNoGo,
    "FALSIFIABLE_TENSION(channel=h_EP,radiation,universality)",
    "NATURALLY_HIDDEN"
  };
  expectBool[
    "VERDICT_REDERIVATION",
    witnessActual === witnessExpected,
    <|"actual" -> witnessActual, "ratified" -> witnessExpected|>
  ];

  section["Bounded source-to-stage predicate manifest"];
  manifestPair = manifestGuards[
    activeMutation === "SOURCE_TO_STAGE_MANIFEST"
  ];
  coverageGuard = manifestPair[[1]];
  categoryGuard = manifestPair[[2]];
  expectBool[
    "SOURCE_TO_STAGE_MANIFEST",
    TrueQ[coverageGuard["ok"]],
    Join[<|"assertion" -> "coverage/count"|>, coverageGuard]
  ];
  expectBool[
    "SOURCE_TO_STAGE_MANIFEST",
    TrueQ[categoryGuard["ok"]],
    <|"assertion" -> "per-predicate category map"|>
  ];

  Print[""];
  Print["VERDICT=", productionVerdict];
  KeyValueMap[
    Print["CHANNEL_", #1, "=", #2] &,
    channelStateAssociation[productionChannels]
  ];
  Print[
    "EP_ADJUDICATION=mass-channel-only; full-decoupled-floor-EP=NOT_EARNED; unqualified-safe-report=False"
  ];
  Print[
    "RADIATION_EXPONENTS=bare:", bareExponent,
    ";pinned:", pinnedExponent,
    ";corrupt:", corruptExponent,
    ";wrong_norm:", wrongExponent
  ];
  Print[
    "RADIATION_MAGNITUDE=SIM_GATED_BY_GUARD_A_NO_NUMERIC_POWER_RATIO_EMITTED"
  ];
  Print[
    "REACHABLE_FALSIFIERS=",
    StringRiffle[
      KeyValueMap[#1 <> "->" <> #2 &, reachableFalsifierMap[normalControls]],
      ","
    ]
  ];
  Print[
    "NATURALLY_HIDDEN=DEFINED_BUT_UNREACHABLE_IN_PRODUCTION;collective-four-channel-earned-stub-only"
  ];
  Print[
    "LORENTZ=c_E=c_gamma_NECESSARY_NOT_SUFFICIENT;real-radiating-extra-scalar;magnitude-remains-SIM_GATED-without-pinned-K_h"
  ];
  Print[
    "GUARD_A_SCOPE=DENYLIST_ONLY:{M_h,c_E,K_h,P_h/P_EM,EP-magnitude,residue-floor};unrelated-numeric-fields-pass"
  ];
  Print[
    "GATE_L_CONNECTION=EARNED_CONNECTION;same-embedding-direction-family-as-light"
  ];
  Print[
    "HARD_WALL=deferred-throat-embedding-solve-pins-{M_h,c_E,K_h};on-pinning-check-stellar-cooling-and-BBN/CMB"
  ];
  Print[
    "FRAMING=FIRST-CLASS_CHARACTERIZED_DEPARTURE;MAPPED+SIM-GATED_NOT_RESOLVED;several-reachable-falsifiers"
  ];
  Print["SOURCE_PREDICATE_TOTAL=", sourcePredicateTotal];
  Print["EXECUTABLE_TOOTH_TOTAL=49"];
  productionVerdict
];


ok = Catch[
  If[
    activeMutation =!= "" && ! MemberQ[toothOrder, activeMutation],
    Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
    Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
    raise["UNKNOWN_MUTATION"]
  ];

  Print["ledger_stage042_charge_coupled_scalar Mathematica audit"];
  Print[
    "ROUTE=LinearSolve exchange + logarithmic-cE Limit falloff + ordered FirstCase channel tables + set-membership verdict + Position/Cases Guard-A walk"
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
  "ledgerStage042Failure",
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
    "OVERALL FAIL: Mathematica stage042 audit did not close"
  ];
  Exit[1]
]
