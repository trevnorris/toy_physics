(* Ledger stage032 Mathematica audit: BC ensembles and sealed R1 landing.

   Standalone, print-only, assert-zero, machine-real-free, and file-I/O-free.
   This route is independent of the SymPy engine: native Solve constructs all
   four ensemble members, native Reduce certifies the MIXED case split, and a
   native typed ladder/oracle plus Hash recomputes the exhaustive digest.
   Stage031's response, self-response, shell, 1/R^2 form, and named core-holder
   fact are consumed symbolically and are not re-derived.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE032_MUTATION";
activeMutation = Environment[mutationEnvironment];
If[!StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];

verdict = "R1_REQUIRED(bc_selection)";
magnitudeCoblocker = "R1_REQUIRED(magnitude)";
committedDigest = "7627417ace0f99a17187a90efa2523ca9b68df8b64f7960d38be2dc6f553ac49";
committedManifestDigest = "e2cfd11b41cdd3d95111f25215b16e917b1ce0a9623929619338a1a83fa81014";
expectedTotal = 23040;
expectedCounts = KeySort@<|
  "MECHANISM_FALSIFIED(wrong_range)" -> 1008,
  "MECHANISM_FALSIFIED(wrong_sign)" -> 504,
  "NO_GO(sector)" -> 11520,
  "R1_REQUIRED(bc_selection)" -> 1152,
  "R1_REQUIRED(magnitude)" -> 252,
  "R1_REQUIRED(mixed_bc_parameters)" -> 1152,
  "R1_REQUIRED(sign_and_magnitude)" -> 2856,
  "R1_REQUIRED(subleading)" -> 504,
  "R1_REQUIRED(variant_selection)" -> 3840,
  "SIGN_EARNED" -> 252|>;
expectedManifestCounts = KeySort@<|
  "PRESERVED" -> 22,
  "REPLACED_BY_STRONGER" -> 12,
  "SCOPED_OUT" -> 10|>;

toothOrder = {
  "PASS_STAGE031_CONSUMED_INPUTS",
  "PASS_STAGE030_TRANSITIVE_D_WITNESS",
  "PASS_SCOPE_FIREWALL",
  "PASS_NONZERO_HA_REQUIRES_CORE_HOLDER",
  "PASS_A_V_FORMULA",
  "PASS_A_J_FORMULA",
  "PASS_A_M_FORMULA",
  "PASS_A_MIXED_FORMULA",
  "PASS_WEAK_SIGNS_GENERAL",
  "PASS_STRICT_SIGNS_NONDEGENERATE_WITNESS",
  "PASS_DEGENERATE_Z_G_ZERO",
  "PASS_A_M_INDEFINITE",
  "PASS_MIXED_ENDPOINTS",
  "PASS_MIXED_INTERIOR_ZERO",
  "PASS_MIXED_REGIME_Q_GT_G",
  "PASS_MIXED_REGIME_Q_EQ_G",
  "PASS_MIXED_REGIME_Q_LT_G",
  "PASS_MIXED_FULL_DOMAIN_SPANS",
  "PASS_DECIDED_CONDITIONAL_TYPING",
  "PASS_AMEND_REPLACE",
  "PASS_AMEND_ADD",
  "PASS_ZERO_LEDGER_13",
  "PASS_S_HOLD_SCOPE",
  "PASS_INTERNAL_INCONSISTENCY_NONE",
  "PASS_BC_ACTUAL_CLASSIFIER",
  "PASS_BC_FREE_CONTROL",
  "PASS_BC_VALUE_CONTROL",
  "PASS_BC_MONOPOLE_CONTROL",
  "PASS_BC_MIXED_CONTROL",
  "PASS_ADMISSIBLE_CLASSES",
  "PASS_REPLACE_TOTALS",
  "PASS_ADD_TOTALS",
  "PASS_VARIANT_UNRESOLVED",
  "PASS_NO_DOUBLE_COUNT",
  "PASS_OUTCOME_NOT_INVARIANT",
  "PASS_MAGNITUDE_FREE_FACTOR",
  "PASS_QMAG_DENSITY_HOOK",
  "PASS_QMAG_MONOPOLE_HOOK",
  "PASS_QMAG_MODULUS_HOOK",
  "PASS_QMAG_CLOSE_RANGE_HOOK",
  "PASS_MAGNITUDE_COBLOCKER",
  "PASS_RANGE_SIGN_FLIP",
  "PASS_RANGE_ZERO_TOUCH",
  "PASS_RANGE_SUBDOMINANT",
  "PASS_UNITS_A",
  "PASS_UNITS_U",
  "PASS_UNITS_F",
  "PASS_TYPED_NEUTRAL_FACTS",
  "PASS_VERDICT_TOTALITY",
  "PASS_VERDICT_PRECEDENCE",
  "PASS_TRUTH_TABLE_TOTAL",
  "PASS_TRUTH_TABLE_COUNTS",
  "PASS_TRUTH_TABLE_DIGEST",
  "PASS_LANDING_OWNERSHIP",
  "PASS_PRODUCTION_LANDING",
  "PASS_TARGET_BLINDNESS",
  "PASS_SOURCE_PREDICATE_MANIFEST"
};
mutationOrder = Append[toothOrder, "MANIFEST_MISDISPOSITION"];
ablationDescriptions = AssociationThread[
  mutationOrder, ("tooth-local upstream mutation for " <> #) & /@ mutationOrder];
ablationDescriptions["MANIFEST_MISDISPOSITION"] =
  "flip AMEND_REPLACE from PRESERVED to SCOPED_OUT while retaining all 44 identifiers";

raise[name_] := Throw[name, "ledgerStage032Failure"];

assertExact[name_, expression_] := Module[{reals},
  reals = Cases[Unevaluated[expression], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FIRST_FAILURE=", name];
    Print["FAIL  ", name, ": machine-real atom(s): ", InputForm[reals]];
    raise[name]
  ]
];

cleanZero[expression_] := FullSimplify[expression] /. ConditionalExpression[0, _] -> 0;

expectZero[name_, residual_, evidence_: None] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FIRST_FAILURE=", name];
    Print["FAIL  ", name, ": residual = ", InputForm[clean]];
    If[evidence =!= None, Print["      evidence = ", InputForm[evidence]]];
    raise[name]
  ]
];

expectBool[name_, condition_, evidence_: None] :=
  expectZero[name, If[TrueQ[condition], 0, 1], evidence];

mut[predicate_, canonical_, altered_] := If[activeMutation === predicate, altered, canonical];

heading[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);

canon[expression_] := ToString[Factor[Cancel[FullSimplify[expression]]], InputForm];

(* ---------------------------------------------------------------------- *)
(* Consumed symbolic inputs and the ownership firewall.                    *)
(* ---------------------------------------------------------------------- *)

consumedStage031 = {
  "m", "m_gg=B_eff*z_g^2/D", "det(m)=z_g^2/D", "S_gg", "A=m_gg*C",
  "NONZERO_HA_REQUIRES_CORE_HOLDER", "s1*s2/(4*pi*R^2)"
};
transitiveStage030 = {"D", "D*=7/4"};
ownedStage032 = {
  "A_V", "A_J", "A_M", "A_MIXED", "Q_AMEND_consistency", "Q_BC_typed_facts",
  "Q_COMBINE", "Q_MAG", "section4_landing"
};
forbiddenRederivations = {
  "response_matrix_derivation", "one_over_R2_derivation", "mouth_mechanism_derivation"
};

(* ---------------------------------------------------------------------- *)
(* Independent native ensemble construction.                              *)
(* ---------------------------------------------------------------------- *)

vRules = First@Solve[{
    sgg y1 + eps mgg y2 == s1 phi,
    eps mgg y1 + sgg y2 == s2 phi
  }, {y1, y2}];
vFunctional = -({s1 phi, s2 phi}.{y1, y2} /. vRules)/2;
derivedAV = Factor[Coefficient[Normal@Series[vFunctional, {eps, 0, 1}], eps]/(s1 s2)];

mRules = First@Solve[{
    yM == q + g,
    aM == mgg yM^2 - 2 g mgg yM
  }, {yM, aM}];
derivedAM = Factor[aM /. mRules];

jRules = First@Solve[{
    yJ == j + g,
    aJ == mgg yJ^2 - 2 yJ mgg yJ
  }, {yJ, aJ}];
derivedAJ = Factor[aJ /. jRules];

mixedRules = First@Solve[{
    yMix == q + g,
    aMix == mgg yMix^2 - 2 (g + lambda q) mgg yMix
  }, {yMix, aMix}];
derivedAMixed = Factor[aMix /. mixedRules];

av = mgg phi^2/sgg^2;
aj = -mgg (j + g)^2;
am = mgg (q^2 - g^2);
amixed = mgg ((1 - 2 lambda) q^2 - 2 lambda q g - g^2);
coefficients = <|"V" -> av, "J" -> aj, "M" -> am, "MIXED" -> amixed|>;
mixedRoot = Factor[First[lambda /. Solve[(amixed/mgg) == 0, lambda, Reals]]];

signExact[expression_] := Which[
  TrueQ[FullSimplify[expression == 0]], 0,
  TrueQ[FullSimplify[expression > 0]], 1,
  TrueQ[FullSimplify[expression < 0]], -1,
  True, Indeterminate
];

strictSignContext[zgValue_, mggValue_, sggValue_, phiValue_, jValue_, gValue_] :=
  TrueQ[zgValue =!= 0 && mggValue > 0 && sggValue > 0 && phiValue =!= 0 && jValue + gValue =!= 0];

(* Native Reduce certificates used below; these are intentionally independent
   of the SymPy substitution/sign route. *)
weakVRegion = Reduce[mgg >= 0 && sgg > 0 && Element[phi, Reals] && av < 0,
  {mgg, sgg, phi}, Reals];
weakJRegion = Reduce[mgg >= 0 && Element[{j, g}, Reals] && aj > 0,
  {mgg, j, g}, Reals];

gtExpression = amixed /. {q -> g0 + d, g -> g0, mgg -> m0};
gtPositiveRegion = Reduce[
  d > 0 && g0 > 0 && m0 > 0 && 0 <= lambda <= 1 && gtExpression > 0,
  lambda, Reals];
gtZeroRegion = Reduce[
  d > 0 && g0 > 0 && m0 > 0 && 0 <= lambda <= 1 && gtExpression == 0,
  lambda, Reals];
gtNegativeRegion = Reduce[
  d > 0 && g0 > 0 && m0 > 0 && 0 <= lambda <= 1 && gtExpression < 0,
  lambda, Reals];

eqExpression = amixed /. {q -> g0, g -> g0, mgg -> m0};
eqZeroRegion = Reduce[
  g0 > 0 && m0 > 0 && 0 <= lambda <= 1 && eqExpression == 0,
  lambda, Reals];
eqNegativeRegion = Reduce[
  g0 > 0 && m0 > 0 && 0 < lambda <= 1 && eqExpression < 0,
  lambda, Reals];

ltExpression = amixed /. {q -> g0, g -> g0 + d, mgg -> m0};
ltNonnegativeRegion = Reduce[
  d > 0 && g0 > 0 && m0 > 0 && 0 <= lambda <= 1 && ltExpression >= 0,
  lambda, Reals];

(* ---------------------------------------------------------------------- *)
(* Source Q-block data.                                                    *)
(* ---------------------------------------------------------------------- *)

zeroLedger = {
  "bulk r_BH,r_B^2H^2,Hrho,Hdelta-rho,Hdt-theta,Hgrad-theta outside Omega_mouth",
  "dynamic J_m modulation and neighbor-source response",
  "r_B-u_L,r_B-divu,r_B-u_T,H-u_T,u_L-u_T and scalar-transverse mixing",
  "cross kinetic and one-time-derivative Berry terms",
  "u_L^2,h^2,higher gradients and cubic/quartic terms",
  "independent primitive B(divu)^2",
  "independent theta_B and phase-drain terms",
  "wall/collar/storage/dissipation terms",
  "dynamic drain and return-kernel responses",
  "direct drain sources and direct h contribution to e_c",
  "field-dependent geon derivatives",
  "viscosity/drag/E2/E3 terms",
  "E4/E5/E1 and mixture terms"
};

replaceAmendment = <|
  "Source" -> "core_holder_retypes_existing_h_source_BC",
  "New" -> {}, "SHold" -> "r_B-1/2 only", "Zeros" -> zeroLedger|>;
addAmendment = <|
  "Source" -> "existing_h_source_BC_unchanged",
  "New" -> {"core_embedding_h_holding_row:R_w_odd"},
  "SHold" -> "r_B-1/2 only", "Zeros" -> zeroLedger|>;

amendmentSectors[rep_, add_] := Module[{sectors = {}},
  If[rep["Source"] =!= "core_holder_retypes_existing_h_source_BC" || rep["New"] =!= {},
    AppendTo[sectors, "replace_ledger"]];
  If[add["Source"] =!= "existing_h_source_BC_unchanged" || Length[add["New"]] =!= 1,
    AppendTo[sectors, "add_ledger"]];
  If[rep["SHold"] =!= "r_B-1/2 only" || add["SHold"] =!= "r_B-1/2 only",
    AppendTo[sectors, "S_hold"]];
  If[rep["Zeros"] =!= zeroLedger || add["Zeros"] =!= zeroLedger,
    AppendTo[sectors, "zero_ledger"]];
  sectors
];

bcDomain = {"DIRICHLET_VALUE", "FIXED_MONOPOLE", "FIXED_SOURCE", "MIXED",
  "UNDETERMINED_ANALYTICALLY"};
outcomeDomain = {"POSITIVE_R2", "NEGATIVE_R2", "NULL", "POSITIVE_WRONG_RANGE",
  "NEGATIVE_WRONG_RANGE", "outcome_not_invariant"};
variantDomain = {"replace", "add", "both", "variant_unresolved"};
magnitudeDomain = {"magnitude_no_free_factor", "magnitude_free_factor"};
tierDomain = {"tier_A_gaps_closed", "tier_A_conditional"};

bcFactQ[BCStatus[kind_, missing_List]] := MemberQ[bcDomain, kind] && VectorQ[missing, StringQ];
bcFactQ[_] := False;
outcomeFactQ[OutcomeFact[value_]] := MemberQ[outcomeDomain, value];
outcomeFactQ[_] := False;
variantFactQ[VariantFact[value_]] := MemberQ[variantDomain, value];
variantFactQ[_] := False;
magnitudeFactQ[MagnitudeFact[value_]] := MemberQ[magnitudeDomain, value];
magnitudeFactQ[_] := False;
tierFactQ[TierFact[value_]] := MemberQ[tierDomain, value];
tierFactQ[_] := False;

factText[BCStatus[kind_, _List]] := kind;
factText[OutcomeFact[value_]] := value;
factText[VariantFact[value_]] := value;
factText[MagnitudeFact[value_]] := value;
factText[TierFact[value_]] := value;

makeEvidence[value_, conormal_, mixed_, variation_, stationary_, barrier_, missing_] := <|
  "Value" -> value, "Conormal" -> conormal, "Mixed" -> mixed,
  "Variation" -> variation, "Stationary" -> stationary,
  "Barrier" -> barrier, "Missing" -> missing|>;

classifyBC[evidence_] := Which[
  TrueQ[evidence["Value"]] && !TrueQ[evidence["Variation"]], BCStatus["DIRICHLET_VALUE", {}],
  TrueQ[evidence["Conormal"]] && TrueQ[evidence["Variation"]], BCStatus["FIXED_MONOPOLE", {}],
  TrueQ[evidence["Mixed"]] && TrueQ[evidence["Variation"]], BCStatus["MIXED", {}],
  evidence["Missing"] === {} && TrueQ[evidence["Variation"]] &&
    !TrueQ[evidence["Stationary"]] && !TrueQ[evidence["Barrier"]], BCStatus["FIXED_SOURCE", {}],
  True, BCStatus["UNDETERMINED_ANALYTICALLY", evidence["Missing"]]
];

actualEvidence = makeEvidence[False, False, False, True, False, False,
  {"missing parent-throat/boundary closure"}];
freeControl = makeEvidence[False, False, False, True, False, False, {}];
valueControl = makeEvidence[True, False, False, False, True, True, {}];
monopoleControl = makeEvidence[False, True, False, True, False, True, {}];
mixedControl = makeEvidence[False, False, True, True, False, True, {}];
admissibleClasses = {"DIRICHLET_VALUE", "FIXED_MONOPOLE", "FIXED_SOURCE", "MIXED"};

replaceTotals = <|
  "DIRICHLET_VALUE" -> mgg phi^2/sgg^2,
  "FIXED_MONOPOLE" -> mgg q^2,
  "FIXED_SOURCE" -> -mgg j^2,
  "MIXED" -> mgg (1 - 2 lambda) q^2|>;
addTotals = <|
  "DIRICHLET_VALUE" -> av,
  "FIXED_MONOPOLE" -> am,
  "FIXED_SOURCE" -> aj,
  "MIXED" -> amixed|>;
combineFacts = <|
  "DIRICHLET_VALUE" -> <|"replace" -> OutcomeFact["POSITIVE_R2"],
    "add" -> OutcomeFact["POSITIVE_R2"]|>,
  "FIXED_MONOPOLE" -> <|"replace" -> OutcomeFact["POSITIVE_R2"],
    "add" -> OutcomeFact["outcome_not_invariant"]|>,
  "FIXED_SOURCE" -> <|"replace" -> OutcomeFact["NEGATIVE_R2"],
    "add" -> OutcomeFact["NEGATIVE_R2"]|>,
  "MIXED" -> <|"replace" -> OutcomeFact["outcome_not_invariant"],
    "add" -> OutcomeFact["outcome_not_invariant"]|>|>;
variantRealization = VariantFact["variant_unresolved"];
magnitudeFact = MagnitudeFact["magnitude_free_factor"];
tierFact = TierFact["tier_A_conditional"];
hooks = <|
  "density" -> "NO(no_local_prediction)",
  "radial_monopole" -> "UNDETERMINED(core source/conormal not fixed)",
  "modulus" -> "NO(continuous core modulus)",
  "close_range" -> "UNDETERMINED(out_of_scope_R_comparable_to_r_e)"|>;

sampleOutcome[values_List] := If[Length[DeleteDuplicates[Sign[values]]] == 1,
  "CONSTANT_OUTCOME", "outcome_not_invariant"];

(* ---------------------------------------------------------------------- *)
(* Independent typed landing and ladder.                                  *)
(* ---------------------------------------------------------------------- *)

typedLandingFactsQ[LandingFacts[qbc_, outcomes_, variant_, magnitude_, tier_, sectors_]] :=
  bcFactQ[qbc] && AssociationQ[outcomes] && Keys[outcomes] === admissibleClasses &&
  And @@ Table[
    AssociationQ[outcomes[cls]] && Keys[outcomes[cls]] === {"replace", "add"} &&
      outcomeFactQ[outcomes[cls]["replace"]] && outcomeFactQ[outcomes[cls]["add"]],
    {cls, admissibleClasses}] && variantFactQ[variant] && magnitudeFactQ[magnitude] &&
  tierFactQ[tier] && VectorQ[sectors, StringQ];
typedLandingFactsQ[_] := False;

neutralTokens[LandingFacts[qbc_, outcomes_, variant_, magnitude_, tier_, sectors_]] :=
  Join[{qbc}, Flatten@Table[outcomes[cls][which], {cls, admissibleClasses},
    {which, {"replace", "add"}}], {variant, magnitude, tier}, sectors];

productionTuple[facts_LandingFacts] := Module[
  {qbc, outcomes, variant, magnitude, tier, sectors, allOutcomes, overall,
   summarize, realized, realizedOutcomes, mixedOutcomes},
  If[!typedLandingFactsQ[facts], Return[$Failed]];
  {qbc, outcomes, variant, magnitude, tier, sectors} = List @@ facts;
  allOutcomes = Flatten@Table[outcomes[cls][which], {cls, admissibleClasses},
    {which, {"replace", "add"}}];
  overall = If[Length[DeleteDuplicates[allOutcomes]] == 1, First[allOutcomes],
    OutcomeFact["outcome_not_invariant"]];
  summarize[which_] := Module[{values = outcomes[#][which] & /@ admissibleClasses},
    If[Length[DeleteDuplicates[values]] == 1, First[values], overall]];
  realized = If[MemberQ[{VariantFact["replace"], VariantFact["add"]}, variant],
    {factText[variant]}, {"replace", "add"}];
  realizedOutcomes = Flatten@Table[outcomes[cls][which], {cls, admissibleClasses},
    {which, realized}];
  mixedOutcomes = outcomes["MIXED"][#] & /@ realized;
  <|
    "qbc" -> factText[qbc],
    "replace_outcome" -> factText[summarize["replace"]],
    "add_outcome" -> factText[summarize["add"]],
    "variant" -> factText[variant],
    "magnitude" -> factText[magnitude],
    "tier" -> factText[tier],
    "internal" -> (sectors =!= {}),
    "all_classes_agree" -> (Length[DeleteDuplicates[realizedOutcomes]] == 1),
    "mixed_range_invariant" -> AllTrue[mixedOutcomes,
      # =!= OutcomeFact["outcome_not_invariant"] &],
    "overall_outcome" -> factText[overall]|>
];

adjudicate[case_Association, illegalPrecedence_: False] := Module[
  {qbc, ro, ao, variant, magnitude, tier, internal, classesAgree, mixedInvariant,
   comparison, variantResolved, classResolved, magnitudeDetermined, unconditional},
  {qbc, ro, ao, variant, magnitude, tier, internal, classesAgree, mixedInvariant} =
    Lookup[case, {"qbc", "replace_outcome", "add_outcome", "variant", "magnitude",
      "tier", "internal", "all_classes_agree", "mixed_range_invariant"}];
  If[TrueQ[internal], Return["NO_GO(sector)"]];
  Which[
    variant === "replace", comparison = ro; variantResolved = True,
    variant === "add", comparison = ao; variantResolved = True,
    ro === ao, comparison = ro; variantResolved = True,
    True, comparison = None; variantResolved = False
  ];
  classResolved = Which[
    qbc === "UNDETERMINED_ANALYTICALLY", TrueQ[classesAgree],
    qbc === "MIXED", tier === "tier_A_gaps_closed" && TrueQ[mixedInvariant],
    True, tier === "tier_A_gaps_closed"
  ];
  magnitudeDetermined = comparison =!= None && comparison =!= "outcome_not_invariant";
  unconditional = classResolved && variantResolved && magnitudeDetermined;
  If[illegalPrecedence && qbc === "UNDETERMINED_ANALYTICALLY",
    Return["R1_REQUIRED(bc_selection)"]];
  If[unconditional,
    Which[
      comparison === "POSITIVE_R2",
        Return[If[magnitude === "magnitude_no_free_factor", "SIGN_EARNED", "R1_REQUIRED(magnitude)"]],
      comparison === "NEGATIVE_R2", Return["MECHANISM_FALSIFIED(wrong_sign)"],
      MemberQ[{"POSITIVE_WRONG_RANGE", "NEGATIVE_WRONG_RANGE"}, comparison],
        Return["MECHANISM_FALSIFIED(wrong_range)"],
      comparison === "NULL", Return["R1_REQUIRED(subleading)"]
    ]
  ];
  If[qbc === "UNDETERMINED_ANALYTICALLY" && !TrueQ[classesAgree],
    Return["R1_REQUIRED(bc_selection)"]];
  If[qbc === "MIXED" && !TrueQ[mixedInvariant],
    Return["R1_REQUIRED(mixed_bc_parameters)"]];
  If[MemberQ[{"both", "variant_unresolved"}, variant] && ro =!= ao,
    Return["R1_REQUIRED(variant_selection)"]];
  If[tier === "tier_A_conditional" || comparison === "outcome_not_invariant",
    Return["R1_REQUIRED(sign_and_magnitude)"]];
  "R1_REQUIRED(unclassified)"
];

declarativeOracle[case_Association] := Module[
  {qbc, ro, ao, variant, magnitude, tier, internal, classesAgree, mixedInvariant,
   common, variantOK, classOK},
  {qbc, ro, ao, variant, magnitude, tier, internal, classesAgree, mixedInvariant} =
    Lookup[case, {"qbc", "replace_outcome", "add_outcome", "variant", "magnitude",
      "tier", "internal", "all_classes_agree", "mixed_range_invariant"}];
  If[TrueQ[internal], Return["NO_GO(sector)"]];
  common = Which[variant === "replace", ro, variant === "add", ao, ro === ao, ro, True, None];
  variantOK = MemberQ[{"replace", "add"}, variant] || ro === ao;
  classOK = Which[
    qbc === "UNDETERMINED_ANALYTICALLY", TrueQ[classesAgree],
    qbc === "MIXED", tier === "tier_A_gaps_closed" && TrueQ[mixedInvariant],
    True, tier === "tier_A_gaps_closed"
  ];
  If[classOK && variantOK && common =!= "outcome_not_invariant",
    Return@Switch[common,
      "POSITIVE_R2", If[magnitude === "magnitude_no_free_factor", "SIGN_EARNED", "R1_REQUIRED(magnitude)"],
      "NEGATIVE_R2", "MECHANISM_FALSIFIED(wrong_sign)",
      "POSITIVE_WRONG_RANGE" | "NEGATIVE_WRONG_RANGE", "MECHANISM_FALSIFIED(wrong_range)",
      "NULL", "R1_REQUIRED(subleading)"
    ]
  ];
  FirstCase[{
      {qbc === "UNDETERMINED_ANALYTICALLY" && !TrueQ[classesAgree], "R1_REQUIRED(bc_selection)"},
      {qbc === "MIXED" && !TrueQ[mixedInvariant], "R1_REQUIRED(mixed_bc_parameters)"},
      {MemberQ[{"both", "variant_unresolved"}, variant] && ro =!= ao, "R1_REQUIRED(variant_selection)"},
      {tier === "tier_A_conditional" || common === "outcome_not_invariant", "R1_REQUIRED(sign_and_magnitude)"}
    }, {True, landing_} :> landing, "R1_REQUIRED(unclassified)"]
];

truthTable[] := Module[
  {tuples, keys, cases, landings, oracleLandings, rows, counts, digest},
  tuples = Tuples[{bcDomain, outcomeDomain, outcomeDomain, variantDomain,
    magnitudeDomain, tierDomain, {False, True}, {False, True}, {False, True}}];
  keys = {"qbc", "replace_outcome", "add_outcome", "variant", "magnitude",
    "tier", "internal", "all_classes_agree", "mixed_range_invariant"};
  cases = AssociationThread[keys, #] & /@ tuples;
  landings = adjudicate /@ cases;
  oracleLandings = declarativeOracle /@ cases;
  rows = MapThread[
    StringRiffle[(If[BooleanQ[#], ToString[#, InputForm], ToString[#]] &) /@ #1, "|"] <> "|" <> #2 &,
    {tuples, landings}];
  counts = KeySort@Counts[landings];
  digest = IntegerString[Hash[StringRiffle[rows, "\n"], "SHA256"], 16, 64];
  <|"Exact" -> (landings === oracleLandings &&
      !AnyTrue[landings, StringContainsQ[#, "unclassified"] &]),
    "Total" -> Length[rows], "Counts" -> counts, "Digest" -> digest, "Rows" -> rows|>
];

productionFacts = LandingFacts[
  classifyBC[actualEvidence], combineFacts, variantRealization, magnitudeFact, tierFact,
  amendmentSectors[replaceAmendment, addAmendment]];
productionCase = productionTuple[productionFacts];
productionLanding = adjudicate[productionCase];

landingOwnershipGuard[facts_LandingFacts, emitted_, carried_: Automatic] := Module[
  {expected, tokens, freshCase, freshLanding},
  If[!typedLandingFactsQ[facts], Return[False]];
  expected = neutralTokens[facts];
  tokens = If[carried === Automatic, expected, carried];
  If[Length[tokens] =!= Length[expected] || !And @@ MapThread[SameQ, {tokens, expected}], Return[False]];
  If[AnyTrue[tokens, StringQ[#] && (StringContainsQ[#, "R1_REQUIRED"] || StringContainsQ[#, "SIGN_EARNED"]) &],
    Return[False]];
  freshCase = productionTuple[facts];
  freshLanding = adjudicate[freshCase];
  emitted === freshLanding
];

(* ---------------------------------------------------------------------- *)
(* Canonical source-to-stage manifest.                                     *)
(* ---------------------------------------------------------------------- *)

sourceToothIDs = {
  "FIELD_IDENTITY_UNITS", "FIELD_PARENT_MAP", "FIELD_LIVE_QH", "ACTION_TRANSCRIPTION",
  "AMEND_REPLACE", "AMEND_ADD", "SHOLD_SCOPE", "ZERO_LEDGER", "MATRIX_DETERMINANT",
  "BC_ACTUAL_GAP", "BC_FREE_CONTROL", "BC_VALUE_CONTROL", "BC_MONOPOLE_CONTROL",
  "BC_MIXED_CONTROL", "FORCE_V_FUNCTIONAL", "FORCE_M_HWORK", "FORCE_J_FUNCTIONAL",
  "FORCE_MIXED_FUNCTIONAL", "MIXED_FULL_RANGE", "FALLOFF", "UNITS_RESTORED",
  "COMBINE_REPLACE", "COMBINE_ADD", "NO_DOUBLE_COUNT", "RANGE_SIGN_FLIP",
  "RANGE_TOUCH_ZERO", "RANGE_SUBDOMINANT", "MAG_FREE_FACTOR", "DENSITY_HOOK",
  "MONOPOLE_HOOK", "MODULUS_HOOK", "VERDICT_TOTALITY", "VERDICT_PRECEDENCE",
  "LANDING_OWNERSHIP", "TARGET_BLINDNESS", "DUAL_ENGINE_TERMS"
};
extraSourceClaims = {
  "EXTRA_VARIANT_UNRESOLVED", "EXTRA_QMAG_CLOSE_RANGE", "EXTRA_SGG_SELF_RESPONSE",
  "EXTRA_G0_RESPONSE_WITNESS", "EXTRA_RESPONSE_MATRIX_ESCAPE_FACTORS",
  "EXTRA_EXTERIOR_CORE_GAP", "EXTRA_SECTION4_DIGEST_COUNTS", "EXTRA_APP_E_ACCEPTANCE"
};

sourceManifest = {
  {"FIELD_IDENTITY_UNITS", "SCOPED_OUT", "STAGE031_MECHANISM"},
  {"FIELD_PARENT_MAP", "SCOPED_OUT", "STAGE031_MECHANISM"},
  {"FIELD_LIVE_QH", "SCOPED_OUT", "STAGE031_MECHANISM"},
  {"ACTION_TRANSCRIPTION", "SCOPED_OUT", "STAGE030_031_CITED"},
  {"AMEND_REPLACE", "PRESERVED", "STAGE032_Q_AMEND"},
  {"AMEND_ADD", "PRESERVED", "STAGE032_Q_AMEND"},
  {"SHOLD_SCOPE", "PRESERVED", "STAGE032_Q_AMEND"},
  {"ZERO_LEDGER", "PRESERVED", "STAGE032_Q_AMEND"},
  {"MATRIX_DETERMINANT", "SCOPED_OUT", "STAGE031_RESPONSE"},
  {"BC_ACTUAL_GAP", "REPLACED_BY_STRONGER", "STAGE032_TYPED_Q_BC"},
  {"BC_FREE_CONTROL", "PRESERVED", "STAGE032_TYPED_Q_BC"},
  {"BC_VALUE_CONTROL", "PRESERVED", "STAGE032_TYPED_Q_BC"},
  {"BC_MONOPOLE_CONTROL", "PRESERVED", "STAGE032_TYPED_Q_BC"},
  {"BC_MIXED_CONTROL", "PRESERVED", "STAGE032_TYPED_Q_BC"},
  {"FORCE_V_FUNCTIONAL", "REPLACED_BY_STRONGER", "STAGE032_FORMULA_SIGN"},
  {"FORCE_M_HWORK", "REPLACED_BY_STRONGER", "STAGE032_FORMULA_SIGN"},
  {"FORCE_J_FUNCTIONAL", "REPLACED_BY_STRONGER", "STAGE032_FORMULA_SIGN"},
  {"FORCE_MIXED_FUNCTIONAL", "REPLACED_BY_STRONGER", "STAGE032_FORMULA_SIGN"},
  {"MIXED_FULL_RANGE", "REPLACED_BY_STRONGER", "STAGE032_THREE_REGIMES"},
  {"FALLOFF", "SCOPED_OUT", "STAGE031_NEUTRAL_SHELL"},
  {"UNITS_RESTORED", "REPLACED_BY_STRONGER", "STAGE032_UNITS_FIREWALL"},
  {"COMBINE_REPLACE", "PRESERVED", "STAGE032_Q_COMBINE"},
  {"COMBINE_ADD", "PRESERVED", "STAGE032_Q_COMBINE"},
  {"NO_DOUBLE_COUNT", "PRESERVED", "STAGE032_Q_COMBINE"},
  {"RANGE_SIGN_FLIP", "PRESERVED", "STAGE032_RANGE"},
  {"RANGE_TOUCH_ZERO", "PRESERVED", "STAGE032_RANGE"},
  {"RANGE_SUBDOMINANT", "PRESERVED", "STAGE032_RANGE"},
  {"MAG_FREE_FACTOR", "PRESERVED", "STAGE032_Q_MAG"},
  {"DENSITY_HOOK", "PRESERVED", "STAGE032_Q_MAG"},
  {"MONOPOLE_HOOK", "PRESERVED", "STAGE032_Q_MAG"},
  {"MODULUS_HOOK", "PRESERVED", "STAGE032_Q_MAG"},
  {"VERDICT_TOTALITY", "REPLACED_BY_STRONGER", "STAGE032_GRID_COUNTS_DIGEST"},
  {"VERDICT_PRECEDENCE", "PRESERVED", "STAGE032_LADDER"},
  {"LANDING_OWNERSHIP", "REPLACED_BY_STRONGER", "STAGE032_UPSTREAM_RELANDING"},
  {"TARGET_BLINDNESS", "PRESERVED", "STAGE032_NEUTRAL_UPSTREAM"},
  {"DUAL_ENGINE_TERMS", "REPLACED_BY_STRONGER", "STAGE032_INDEPENDENT_ENGINES"},
  {"EXTRA_VARIANT_UNRESOLVED", "PRESERVED", "STAGE032_Q_COMBINE"},
  {"EXTRA_QMAG_CLOSE_RANGE", "PRESERVED", "STAGE032_Q_MAG"},
  {"EXTRA_SGG_SELF_RESPONSE", "SCOPED_OUT", "STAGE031_CONSUMED"},
  {"EXTRA_G0_RESPONSE_WITNESS", "SCOPED_OUT", "STAGE031_RESPONSE"},
  {"EXTRA_RESPONSE_MATRIX_ESCAPE_FACTORS", "SCOPED_OUT", "STAGE031_RESPONSE"},
  {"EXTRA_EXTERIOR_CORE_GAP", "SCOPED_OUT", "STAGE031_CONSUMED_NAMED_FACT"},
  {"EXTRA_SECTION4_DIGEST_COUNTS", "REPLACED_BY_STRONGER", "STAGE032_COMMITTED_LITERAL"},
  {"EXTRA_APP_E_ACCEPTANCE", "REPLACED_BY_STRONGER", "STAGE032_STANDALONE_ABLATIONS"}
};
lexicographicCodeLess[left_List, right_List] := Module[{limit, index},
  limit = Min[Length[left], Length[right]];
  index = SelectFirst[Range[limit], left[[#]] =!= right[[#]] &, Missing["NotFound"]];
  If[MissingQ[index], Length[left] < Length[right], left[[index]] < right[[index]]]
];
canonicalManifestText[manifest_] := Module[{rows},
  rows = StringRiffle[#, "|"] & /@ manifest;
  StringRiffle[
    Sort[rows, lexicographicCodeLess[ToCharacterCode[#1], ToCharacterCode[#2]] &], "\n"]
];
manifestSHA256[manifest_] :=
  IntegerString[Hash[canonicalManifestText[manifest], "SHA256"], 16, 64];

(* ---------------------------------------------------------------------- *)
(* Assertions.                                                            *)
(* ---------------------------------------------------------------------- *)

ok = Catch[
  If[activeMutation =!= "" && !MemberQ[mutationOrder, activeMutation],
    Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
    Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
    raise["UNKNOWN_MUTATION"]
  ];

  Print["ledger_stage032_electric_sign_bc_ensembles_landing Mathematica audit"];
  Print["CONSUMES_STAGE031={m,m_gg=B_eff*z_g^2/D,detm=z_g^2/D,S_gg,A=m_gg*C,1/R^2,NONZERO_HA_REQUIRES_CORE_HOLDER}"];
  Print["CONSUMES_STAGE030_TRANSITIVELY={D,D*=7/4}"];
  If[activeMutation =!= "",
    Print["ACTIVE_MUTATION=", activeMutation];
    Print["MUTATED_PRIMITIVE=", ablationDescriptions[activeMutation]]
  ];

  heading["Consumed upstream bundle and ownership firewall"];
  consumed = mut["PASS_STAGE031_CONSUMED_INPUTS", consumedStage031,
    DeleteCases[consumedStage031, "S_gg"]];
  expectBool["PASS_STAGE031_CONSUMED_INPUTS",
    consumed === consumedStage031 && Sort[consumed] === Sort[{
      "m", "m_gg=B_eff*z_g^2/D", "det(m)=z_g^2/D", "S_gg", "A=m_gg*C",
      "NONZERO_HA_REQUIRES_CORE_HOLDER", "s1*s2/(4*pi*R^2)"}], consumed];
  dStar = mut["PASS_STAGE030_TRANSITIVE_D_WITNESS", 7/4, 5/4];
  expectZero["PASS_STAGE030_TRANSITIVE_D_WITNESS", dStar - 7/4];
  owned = mut["PASS_SCOPE_FIREWALL", ownedStage032,
    Append[ownedStage032, "response_matrix_derivation"]];
  expectBool["PASS_SCOPE_FIREWALL", Intersection[owned, forbiddenRederivations] === {}, owned];
  coreFact = mut["PASS_NONZERO_HA_REQUIRES_CORE_HOLDER",
    "NONZERO_HA_REQUIRES_CORE_HOLDER", ""];
  expectBool["PASS_NONZERO_HA_REQUIRES_CORE_HOLDER",
    coreFact === "NONZERO_HA_REQUIRES_CORE_HOLDER"];
  Print["      stage031 mechanism cited; stage030 D*=7/4 cited transitively; no mechanism re-derived"];

  heading["Four DECIDED conditional ensembles via native Solve/Reduce"];
  avTest = mut["PASS_A_V_FORMULA", derivedAV, derivedAV + mgg phi^2];
  expectZero["PASS_A_V_FORMULA", avTest - av, <|"A_V" -> avTest|>];
  ajTest = mut["PASS_A_J_FORMULA", derivedAJ, -derivedAJ];
  expectZero["PASS_A_J_FORMULA", ajTest - aj, <|"A_J" -> ajTest|>];
  amTest = mut["PASS_A_M_FORMULA", derivedAM, mgg q^2];
  expectZero["PASS_A_M_FORMULA", amTest - am, <|"A_M" -> amTest|>];
  mixedTest = mut["PASS_A_MIXED_FORMULA", derivedAMixed,
    derivedAMixed + 2 lambda q g mgg];
  expectZero["PASS_A_MIXED_FORMULA", mixedTest - amixed, <|"A_MIXED" -> mixedTest|>];

  weakOK = mut["PASS_WEAK_SIGNS_GENERAL",
    weakVRegion === False && weakJRegion === False, False];
  expectBool["PASS_WEAK_SIGNS_GENERAL", weakOK,
    <|"Reduce[A_V<0]" -> weakVRegion, "Reduce[A_J>0]" -> weakJRegion|>];

  strictPhi = mut["PASS_STRICT_SIGNS_NONDEGENERATE_WITNESS", 3, 0];
  strictRules = {mgg -> 8/7, sgg -> 2, phi -> strictPhi, j -> 2, g -> 1};
  strictActive = strictSignContext[1, 8/7, 2, strictPhi, 2, 1];
  strictV = FullSimplify[av /. strictRules];
  strictJ = FullSimplify[aj /. strictRules];
  expectBool["PASS_STRICT_SIGNS_NONDEGENERATE_WITNESS",
    strictActive && strictV > 0 && strictJ < 0,
    <|"m_gg" -> 8/7, "S_gg" -> 2, "A_V" -> strictV, "A_J" -> strictJ|>];

  degenerateZg = mut["PASS_DEGENERATE_Z_G_ZERO", 0, 1];
  degenerateMgg = beff degenerateZg^2/delta;
  degenerateCoefficients = FullSimplify[Values[coefficients] /. mgg -> degenerateMgg];
  degenerateStrict = strictSignContext[degenerateZg, degenerateMgg, 1, 1, 1, 1];
  expectBool["PASS_DEGENERATE_Z_G_ZERO",
    degenerateCoefficients === {0, 0, 0, 0} && !degenerateStrict,
    <|"z_g" -> degenerateZg, "A_X" -> degenerateCoefficients,
      "strict_asserts_active" -> degenerateStrict|>];

  negativeQ = mut["PASS_A_M_INDEFINITE", 1, 3];
  mPositive = am /. {mgg -> 1, q -> 3, g -> 1};
  mNegative = am /. {mgg -> 1, q -> negativeQ, g -> 2};
  expectBool["PASS_A_M_INDEFINITE", mPositive > 0 && mNegative < 0,
    <|"positive" -> mPositive, "negative" -> mNegative|>];

  endpointOne = mut["PASS_MIXED_ENDPOINTS", 1, 0];
  endpoints = Factor /@ {amixed /. lambda -> 0, amixed /. lambda -> endpointOne};
  expectBool["PASS_MIXED_ENDPOINTS",
    TrueQ[FullSimplify[endpoints[[1]] - am] === 0] &&
      TrueQ[FullSimplify[endpoints[[2]] + mgg (q + g)^2] === 0], endpoints];

  rootTest = mut["PASS_MIXED_INTERIOR_ZERO", mixedRoot, mixedRoot + 1/2];
  expectBool["PASS_MIXED_INTERIOR_ZERO",
    TrueQ[FullSimplify[rootTest - (q - g)/(2 q)] === 0] &&
      TrueQ[FullSimplify[amixed /. lambda -> rootTest] === 0], rootTest];

  rootGT = d/(2 (d + g0));
  positiveLambda = mut["PASS_MIXED_REGIME_Q_GT_G", rootGT/2, (1 + rootGT)/2];
  gtValues = FullSimplify[gtExpression /. lambda -> #,
      d > 0 && g0 > 0 && m0 > 0] & /@ {positiveLambda, rootGT, (1 + rootGT)/2};
  gtSigns = FullSimplify[Sign /@ gtValues, d > 0 && g0 > 0 && m0 > 0];
  gtReduceOK = TrueQ[FullSimplify[
      Equivalent[FullSimplify[gtPositiveRegion, d > 0 && g0 > 0 && m0 > 0],
        0 <= lambda < rootGT] &&
      Equivalent[FullSimplify[gtZeroRegion, d > 0 && g0 > 0 && m0 > 0],
        lambda == rootGT] &&
      Equivalent[FullSimplify[gtNegativeRegion, d > 0 && g0 > 0 && m0 > 0],
        rootGT < lambda <= 1],
      d > 0 && g0 > 0 && m0 > 0]];
  expectBool["PASS_MIXED_REGIME_Q_GT_G",
    gtReduceOK && gtSigns === {1, 0, -1},
    <|"lambda*" -> rootGT, "Reduce" -> {gtPositiveRegion, gtZeroRegion, gtNegativeRegion},
      "signs" -> gtSigns|>];

  eqZeroLambda = mut["PASS_MIXED_REGIME_Q_EQ_G", 0, 1/2];
  eqZero = FullSimplify[eqExpression /. lambda -> eqZeroLambda,
    g0 > 0 && m0 > 0];
  eqInterior = FullSimplify[eqExpression /. lambda -> 1/2,
    g0 > 0 && m0 > 0];
  eqReduceOK = TrueQ[FullSimplify[eqZeroRegion /. {g0 -> 1, m0 -> 1}] === (lambda == 0)] &&
    TrueQ[FullSimplify[eqNegativeRegion /. {g0 -> 1, m0 -> 1}] === (0 < lambda <= 1)] &&
    eqZeroRegion =!= False && eqNegativeRegion =!= False;
  expectBool["PASS_MIXED_REGIME_Q_EQ_G",
    eqReduceOK && eqZero === 0 && TrueQ[FullSimplify[eqInterior < 0, g0 > 0 && m0 > 0]],
    <|"Reduce" -> {eqZeroRegion, eqNegativeRegion}, "A(0)" -> eqZero,
      "A(1/2)" -> eqInterior|>];

  rootLTCanonical = -d/(2 g0);
  rootLT = mut["PASS_MIXED_REGIME_Q_LT_G", rootLTCanonical, -rootLTCanonical];
  ltEndpoints = FullSimplify[ltExpression /. lambda -> #,
      d > 0 && g0 > 0 && m0 > 0] & /@ {0, 1};
  expectBool["PASS_MIXED_REGIME_Q_LT_G",
    TrueQ[FullSimplify[rootLT < 0, d > 0 && g0 > 0]] &&
      ltNonnegativeRegion === False &&
      And @@ (TrueQ[FullSimplify[# < 0, d > 0 && g0 > 0 && m0 > 0]] & /@ ltEndpoints),
    <|"lambda*" -> rootLT, "Reduce[A>=0]" -> ltNonnegativeRegion,
      "endpoints" -> ltEndpoints|>];

  spanMeaning = mut["PASS_MIXED_FULL_DOMAIN_SPANS", "FULL_PARAMETER_AND_MAGNITUDE_DOMAIN",
    "EVERY_FIXED_Q_OVER_G"];
  fullSamples = amixed /. {
    {mgg -> 1, q -> 3, g -> 1, lambda -> 0},
    {mgg -> 1, q -> 1, g -> 1, lambda -> 0},
    {mgg -> 1, q -> 1, g -> 2, lambda -> 0}
  };
  expectBool["PASS_MIXED_FULL_DOMAIN_SPANS",
    spanMeaning === "FULL_PARAMETER_AND_MAGNITUDE_DOMAIN" && (Sign /@ fullSamples) === {1, 0, -1},
    <|"meaning" -> spanMeaning, "cross-regime signs" -> (Sign /@ fullSamples)|>];

  typing = mut["PASS_DECIDED_CONDITIONAL_TYPING",
    {"coefficients=DECIDED_given_class", "bc_selection=R1", "magnitude=R1",
      "mixed_parameters=R1", "variant=R1"},
    {"coefficients=DECIDED_given_class", "bc_selection=DECIDED", "magnitude=R1",
      "mixed_parameters=R1", "variant=R1"}];
  expectBool["PASS_DECIDED_CONDITIONAL_TYPING",
    typing === {"coefficients=DECIDED_given_class", "bc_selection=R1", "magnitude=R1",
      "mixed_parameters=R1", "variant=R1"}, typing];
  Print["      native Solve coefficients; native Reduce three-regime split; strict signs witness-only"];

  heading["Q-AMEND, typed Q-BC, Q-COMBINE, and Q-MAG"];
  rep = mut["PASS_AMEND_REPLACE", replaceAmendment,
    Join[replaceAmendment, <|"New" -> {"illegal"}|>]];
  expectBool["PASS_AMEND_REPLACE",
    !MemberQ[amendmentSectors[rep, addAmendment], "replace_ledger"], rep];
  add = mut["PASS_AMEND_ADD", addAmendment,
    Join[addAmendment, <|"New" -> {"row1", "row2"}|>]];
  expectBool["PASS_AMEND_ADD",
    !MemberQ[amendmentSectors[replaceAmendment, add], "add_ledger"], add];
  zeros = mut["PASS_ZERO_LEDGER_13", zeroLedger, Most[zeroLedger]];
  expectBool["PASS_ZERO_LEDGER_13", Length[zeros] === 13 && zeros === zeroLedger, Length[zeros]];
  shold = mut["PASS_S_HOLD_SCOPE", "r_B-1/2 only", "r_B-1/2 and h"];
  expectBool["PASS_S_HOLD_SCOPE", shold === "r_B-1/2 only", shold];
  sectors = mut["PASS_INTERNAL_INCONSISTENCY_NONE",
    amendmentSectors[replaceAmendment, addAmendment], {"add_ledger"}];
  expectBool["PASS_INTERNAL_INCONSISTENCY_NONE", sectors === {},
    <|"internal_inconsistency" -> If[sectors === {}, "none", sectors]|>];

  actual = mut["PASS_BC_ACTUAL_CLASSIFIER", actualEvidence,
    Join[actualEvidence, <|"Missing" -> {}|>]];
  expectBool["PASS_BC_ACTUAL_CLASSIFIER",
    classifyBC[actual] === BCStatus["UNDETERMINED_ANALYTICALLY",
      {"missing parent-throat/boundary closure"}], classifyBC[actual]];
  free = mut["PASS_BC_FREE_CONTROL", freeControl,
    Join[freeControl, <|"Value" -> True, "Variation" -> False|>]];
  expectBool["PASS_BC_FREE_CONTROL", classifyBC[free] === BCStatus["FIXED_SOURCE", {}], classifyBC[free]];
  value = mut["PASS_BC_VALUE_CONTROL", valueControl,
    Join[valueControl, <|"Value" -> False, "Variation" -> True, "Missing" -> {"lost_value"}|>]];
  expectBool["PASS_BC_VALUE_CONTROL", classifyBC[value] === BCStatus["DIRICHLET_VALUE", {}], classifyBC[value]];
  mono = mut["PASS_BC_MONOPOLE_CONTROL", monopoleControl,
    Join[monopoleControl, <|"Conormal" -> False, "Missing" -> {"lost_flux"}|>]];
  expectBool["PASS_BC_MONOPOLE_CONTROL", classifyBC[mono] === BCStatus["FIXED_MONOPOLE", {}], classifyBC[mono]];
  mixControl = mut["PASS_BC_MIXED_CONTROL", mixedControl,
    Join[mixedControl, <|"Mixed" -> False, "Missing" -> {"lost_mix"}|>]];
  expectBool["PASS_BC_MIXED_CONTROL", classifyBC[mixControl] === BCStatus["MIXED", {}], classifyBC[mixControl]];
  classes = mut["PASS_ADMISSIBLE_CLASSES", admissibleClasses, Most[admissibleClasses]];
  expectBool["PASS_ADMISSIBLE_CLASSES", classes === admissibleClasses && Length[classes] === 4, classes];

  replaceTest = Association[replaceTotals];
  If[activeMutation === "PASS_REPLACE_TOTALS",
    replaceTest["FIXED_MONOPOLE"] = replaceTest["FIXED_MONOPOLE"] + mgg g^2];
  replaceExpected = {av, mgg q^2, -mgg j^2, mgg (1 - 2 lambda) q^2};
  expectBool["PASS_REPLACE_TOTALS",
    And @@ MapThread[TrueQ[FullSimplify[#1 - #2] === 0] &, {Values[replaceTest], replaceExpected}], replaceTest];
  addTest = Association[addTotals];
  If[activeMutation === "PASS_ADD_TOTALS", addTest["FIXED_MONOPOLE"] = mgg q^2];
  expectBool["PASS_ADD_TOTALS",
    And @@ MapThread[TrueQ[FullSimplify[#1 - #2] === 0] &,
      {Values[addTest], {av, am, aj, amixed}}], addTest];
  variantTest = mut["PASS_VARIANT_UNRESOLVED", variantRealization, VariantFact["replace"]];
  expectBool["PASS_VARIANT_UNRESOLVED", variantTest === VariantFact["variant_unresolved"], variantTest];
  rebuiltM = Factor[mgg (q + g)^2 - 2 g mgg (q + g)];
  rebuiltM = mut["PASS_NO_DOUBLE_COUNT", rebuiltM, rebuiltM + mgg g^2];
  expectZero["PASS_NO_DOUBLE_COUNT", rebuiltM - am, rebuiltM];
  outcomeTest = mut["PASS_OUTCOME_NOT_INVARIANT", combineFacts,
    AssociationThread[admissibleClasses,
      ConstantArray[<|"replace" -> OutcomeFact["POSITIVE_R2"],
        "add" -> OutcomeFact["POSITIVE_R2"]|>, 4]]];
  flatOutcomes = Flatten@Table[outcomeTest[cls][which], {cls, admissibleClasses},
    {which, {"replace", "add"}}];
  overallTest = If[Length[DeleteDuplicates[flatOutcomes]] == 1, First[flatOutcomes],
    OutcomeFact["outcome_not_invariant"]];
  expectBool["PASS_OUTCOME_NOT_INVARIANT", overallTest === OutcomeFact["outcome_not_invariant"], overallTest];

  magnitudeFactors = mut["PASS_MAGNITUDE_FREE_FACTOR", {"c_a", "c_xi"}, {"c_a"}];
  expectBool["PASS_MAGNITUDE_FREE_FACTOR",
    Sort[magnitudeFactors] === {"c_a", "c_xi"} && magnitudeFact === MagnitudeFact["magnitude_free_factor"],
    magnitudeFactors];
  densityHook = mut["PASS_QMAG_DENSITY_HOOK", hooks["density"], "YES(local_prediction)"];
  expectBool["PASS_QMAG_DENSITY_HOOK", densityHook === "NO(no_local_prediction)", densityHook];
  monopoleHook = mut["PASS_QMAG_MONOPOLE_HOOK", hooks["radial_monopole"], "DECIDED(nonzero)"];
  expectBool["PASS_QMAG_MONOPOLE_HOOK", StringStartsQ[monopoleHook, "UNDETERMINED"], monopoleHook];
  modulusHook = mut["PASS_QMAG_MODULUS_HOOK", hooks["modulus"], "YES(universal)"];
  expectBool["PASS_QMAG_MODULUS_HOOK", StringStartsQ[modulusHook, "NO("], modulusHook];
  closeHook = mut["PASS_QMAG_CLOSE_RANGE_HOOK", hooks["close_range"], "DECIDED"];
  expectBool["PASS_QMAG_CLOSE_RANGE_HOOK", StringStartsQ[closeHook, "UNDETERMINED"], closeHook];
  coblocker = mut["PASS_MAGNITUDE_COBLOCKER", magnitudeCoblocker, "R1_REQUIRED(bc_selection)"];
  expectBool["PASS_MAGNITUDE_COBLOCKER", coblocker === "R1_REQUIRED(magnitude)", coblocker];

  heading["Range controls and units firewall"];
  signFlip = mut["PASS_RANGE_SIGN_FLIP", sampleOutcome[{-1, 0, 1}], "CONSTANT_OUTCOME"];
  expectBool["PASS_RANGE_SIGN_FLIP", signFlip === "outcome_not_invariant", signFlip];
  zeroTouch = mut["PASS_RANGE_ZERO_TOUCH", sampleOutcome[{1, 0, 1}], "CONSTANT_OUTCOME"];
  expectBool["PASS_RANGE_ZERO_TOUCH", zeroTouch === "outcome_not_invariant", zeroTouch];
  subdominant = mut["PASS_RANGE_SUBDOMINANT", sampleOutcome[{3/4, 1/2}], "outcome_not_invariant"];
  expectBool["PASS_RANGE_SUBDOMINANT", subdominant === "CONSTANT_OUTCOME", subdominant];

  dimL = {1, 0, 0}; dimE = {2, -2, 1}; dimA = dimE + dimL;
  aDimension = mut["PASS_UNITS_A", dimA, dimE];
  expectBool["PASS_UNITS_A", aDimension === {3, -2, 1}, <|"[A]" -> aDimension, "target" -> "E L"|>];
  uDimension = mut["PASS_UNITS_U", dimA - dimL, dimA];
  expectBool["PASS_UNITS_U", uDimension === dimE, <|"[U]" -> uDimension, "target" -> "E"|>];
  fDimension = mut["PASS_UNITS_F", dimA - 2 dimL, dimE];
  expectBool["PASS_UNITS_F", fDimension === {1, -2, 1}, <|"[F]" -> fDimension, "target" -> "E/L"|>];

  heading["Sealed 23040-cell section-4 ladder and landing ownership"];
  If[activeMutation === "", Print["PROGRESS: native exhaustive 23040-cell landing table and SHA-256 digest"]];
  truth = truthTable[];
  typedTokens = neutralTokens[productionFacts];
  If[activeMutation === "PASS_TYPED_NEUTRAL_FACTS", typedTokens = Append[typedTokens, verdict]];
  expectBool["PASS_TYPED_NEUTRAL_FACTS",
    typedLandingFactsQ[productionFacts] && Length[typedTokens] === Length[neutralTokens[productionFacts]] &&
      AllTrue[typedTokens, !(StringQ[#] && StringContainsQ[#, "R1_REQUIRED"]) &], typedTokens];

  truthExact = mut["PASS_VERDICT_TOTALITY", truth["Exact"], False];
  expectBool["PASS_VERDICT_TOTALITY",
    truthExact && AllTrue[Keys[truth["Counts"]], !StringContainsQ[#, "unclassified"] &],
    <|"exact" -> truthExact, "landings" -> Keys[truth["Counts"]]|>];
  precedenceWitness = <|
    "qbc" -> "UNDETERMINED_ANALYTICALLY", "replace_outcome" -> "POSITIVE_R2",
    "add_outcome" -> "POSITIVE_R2", "variant" -> "both",
    "magnitude" -> "magnitude_no_free_factor", "tier" -> "tier_A_conditional",
    "internal" -> False, "all_classes_agree" -> True, "mixed_range_invariant" -> True|>;
  precedenceLanding = adjudicate[precedenceWitness,
    activeMutation === "PASS_VERDICT_PRECEDENCE"];
  expectBool["PASS_VERDICT_PRECEDENCE",
    precedenceLanding === "SIGN_EARNED" && productionLanding === verdict,
    <|"invariance witness" -> precedenceLanding, "production" -> productionLanding|>];

  totalTest = mut["PASS_TRUTH_TABLE_TOTAL", truth["Total"], truth["Total"] + 1];
  expectZero["PASS_TRUTH_TABLE_TOTAL", totalTest - expectedTotal, totalTest];
  countsTest = Association[truth["Counts"]];
  If[activeMutation === "PASS_TRUTH_TABLE_COUNTS",
    countsTest[verdict] = countsTest[verdict] + 1];
  expectBool["PASS_TRUTH_TABLE_COUNTS",
    KeySort[countsTest] === expectedCounts && Total[Values[countsTest]] === expectedTotal, countsTest];
  digestRows = truth["Rows"];
  If[activeMutation === "PASS_TRUTH_TABLE_DIGEST",
    digestRows = ReplacePart[digestRows, 1 -> First[digestRows] <> "|MUTATED"]];
  digestTest = IntegerString[Hash[StringRiffle[digestRows, "\n"], "SHA256"], 16, 64];
  expectBool["PASS_TRUTH_TABLE_DIGEST", digestTest === committedDigest,
    <|"computed" -> digestTest, "committed" -> committedDigest|>];

  alternativeFacts = productionFacts /. LandingFacts[_, outcomes_, variant_, magnitude_, tier_, sectors_] :>
    LandingFacts[BCStatus["FIXED_SOURCE", {}], outcomes, variant, magnitude, tier, sectors];
  alternativeCase = productionTuple[alternativeFacts];
  alternativeLanding = adjudicate[alternativeCase];
  ownershipFacts = mut["PASS_LANDING_OWNERSHIP", productionFacts, alternativeFacts];
  ownershipCase = productionTuple[ownershipFacts];
  ownershipLanding = adjudicate[ownershipCase];
  injectedTokens = Append[neutralTokens[productionFacts], productionLanding];
  ownershipOK = landingOwnershipGuard[ownershipFacts, productionLanding] &&
    !landingOwnershipGuard[productionFacts, productionLanding, injectedTokens] &&
    ownershipCase["overall_outcome"] === productionTuple[ownershipFacts]["overall_outcome"] &&
    ownershipCase["all_classes_agree"] === productionTuple[ownershipFacts]["all_classes_agree"] &&
    alternativeLanding === "R1_REQUIRED(sign_and_magnitude)" && alternativeLanding =!= verdict &&
    ownershipLanding === productionLanding;
  expectBool["PASS_LANDING_OWNERSHIP", ownershipOK,
    <|"upstream_qbc" -> ownershipCase["qbc"],
      "recomputed_overall" -> ownershipCase["overall_outcome"],
      "recomputed_all_classes_agree" -> ownershipCase["all_classes_agree"],
      "fresh_landing" -> ownershipLanding,
      "named_different_landing" -> alternativeLanding,
      "injected_landing_rejected" -> !landingOwnershipGuard[productionFacts, productionLanding, injectedTokens]|>];

  productionFactsTest = mut["PASS_PRODUCTION_LANDING", productionFacts, alternativeFacts];
  liveCase = productionTuple[productionFactsTest];
  liveLanding = adjudicate[liveCase];
  expectBool["PASS_PRODUCTION_LANDING",
    liveCase["overall_outcome"] === "outcome_not_invariant" &&
      !liveCase["all_classes_agree"] && liveLanding === verdict && !liveCase["internal"],
    <|"tuple" -> liveCase, "landing" -> liveLanding|>];

  upstreamVocabulary = Join[bcDomain, outcomeDomain, variantDomain, magnitudeDomain, tierDomain,
    {"missing parent-throat/boundary closure"}];
  targetToken = mut["PASS_TARGET_BLINDNESS", "outcome_not_invariant", "SIGN_EARNED"];
  expectBool["PASS_TARGET_BLINDNESS",
    MemberQ[upstreamVocabulary, targetToken] &&
      FreeQ[ToLowerCase[ToString[Values[coefficients], InputForm]], "targetsign"] &&
      AllTrue[upstreamVocabulary, !StringContainsQ[#, "R1_REQUIRED"] &], targetToken];
  Print["      cells=", truth["Total"], "; digest=", truth["Digest"],
    "; first production match=", productionLanding];

  heading["Canonical source-to-stage predicate manifest"];
  manifestTest = mut["PASS_SOURCE_PREDICATE_MANIFEST", sourceManifest, Most[sourceManifest]];
  manifestTest = mut["MANIFEST_MISDISPOSITION", manifestTest,
    ReplacePart[sourceManifest, {5, 2} -> "SCOPED_OUT"]];
  identifiers = manifestTest[[All, 1]];
  sourcePart = Take[identifiers, UpTo[Length[sourceToothIDs]]];
  extraPart = Drop[identifiers, UpTo[Length[sourceToothIDs]]];
  dispositions = DeleteDuplicates[manifestTest[[All, 2]]];
  partitionCounts = KeySort@Counts[manifestTest[[All, 2]]];
  manifestDigest = manifestSHA256[manifestTest];
  expectBool["PASS_SOURCE_PREDICATE_MANIFEST",
    sourcePart === sourceToothIDs && extraPart === extraSourceClaims &&
      Length[identifiers] === Length[DeleteDuplicates[identifiers]] === 44 &&
      Sort[dispositions] === Sort[{"PRESERVED", "REPLACED_BY_STRONGER", "SCOPED_OUT"}] &&
      AllTrue[manifestTest[[All, 3]], StringStartsQ[#, "STAGE"] &] &&
      partitionCounts === expectedManifestCounts && manifestDigest === committedManifestDigest,
    <|"entries" -> Length[manifestTest], "dispositions" -> dispositions,
      "partition" -> partitionCounts, "digest" -> manifestDigest|>];
  Print["      entries=", Length[manifestTest], "; partition=", partitionCounts,
    "; digest=", manifestDigest];

  Print[""];
  Print["INTERNAL_INCONSISTENCY=none"];
  Print["CO_BLOCKER=R1_REQUIRED(magnitude)"];
  Print["RESOLVER=SIM_DEFERRED(parent-throat boundary functional/barrier/map s->h_A)"];
  Print["VERDICT_TOKEN: ", verdict];
  If[activeMutation =!= "",
    Print["FIRST_FAILURE=MUTATION_DID_NOT_FIRE"];
    raise["MUTATION_DID_NOT_FIRE"]
  ];
  True,
  "ledgerStage032Failure",
  Function[{message, tag}, False]
];

Print[""];
Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
If[TrueQ[ok],
  Print["OVERALL PASS: Mathematica verified stage032 far-field ensembles and sealed R1 landing"];
  Exit[0],
  Print["OVERALL FAIL: Mathematica stage032 audit did not close"];
  Exit[1]
]
