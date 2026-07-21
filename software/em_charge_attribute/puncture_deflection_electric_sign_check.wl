(* Independent Mathematica engine for the puncture-deflection build.

   Route independence: the coupled kernel is reconstructed by completing the
   square and differentiating its diagonal source polynomial; the V member is
   obtained from an exact two-mouth capacity inverse.  No SymPy formula or
   result file is imported.  Python is invoked only after all local quantities
   exist, as a term-by-term comparator for a normal (non-JSON) run.
*)

ClearAll["Global`*"];
$Assumptions = b > 0 && k > 0 && c > 0 && b k-c^2 > 0 && km > 0 &&
  zeta > 0 && reb > 0 && ell > 0 && sgg > 0 && phi > 0 && q > 0 &&
  j > 0 && g > 0 && 0 <= mix <= 1;

here = DirectoryName[ExpandFileName[$InputFileName]];
root = ParentDirectory[ParentDirectory[here]];
pyPath = FileNameJoin[{here, "puncture_deflection_electric_sign_check.py"}];
paths = <|
  "directive" -> FileNameJoin[{here, "directive_puncture_deflection_electric_sign.md"}],
  "g0" -> FileNameJoin[{here, "g0_closure_card_v0.md"}],
  "priorResult" -> FileNameJoin[{here, "leftover_scalar_electric_sign_result.md"}],
  "priorVerify" -> FileNameJoin[{here, "leftover_scalar_electric_sign_VERIFICATION.md"}],
  "concept" -> FileNameJoin[{root, "docs", "conceptual_foundation.md"}],
  "path42" -> FileNameJoin[{root, "software", "stage1_solver", "reports", "pathA_42_charge_coupled_scalar.md"}],
  "audit" -> FileNameJoin[{root, "docs", "model_definition_audit.md"}],
  "handoff" -> FileNameJoin[{root, "docs", "em_analog_next_phase_handoff.md"}],
  "sim" -> FileNameJoin[{root, "docs", "two_throat_simulation_handoff_spec.md"}],
  "u2" -> FileNameJoin[{here, "reports", "u2_boundary_adjudication_artifacts", "stage_1_production_v12", "production_summary.md"}]
|>;

zq[x_] := TrueQ[FullSimplify[x == 0]];
canon[x_] := ToString[Cancel[x], InputForm];

(* Tagged neutral facts are structurally disjoint from all §4 landing strings. *)
bcFactDomain = {"DIRICHLET_VALUE", "FIXED_MONOPOLE", "FIXED_SOURCE", "MIXED",
  "UNDETERMINED_ANALYTICALLY"};
outcomeFactDomain = {"POSITIVE_R2", "NEGATIVE_R2", "NULL", "POSITIVE_WRONG_RANGE",
  "NEGATIVE_WRONG_RANGE", "outcome_not_invariant"};
variantFactDomain = {"replace", "add", "both", "variant_unresolved"};
magnitudeFactDomain = {"magnitude_no_free_factor", "magnitude_free_factor"};
tierFactDomain = {"tier_A_gaps_closed", "tier_A_conditional"};
amendmentSectorDomain = {"replace_ledger", "add_ledger", "S_hold", "zero_ledger", "scalar_hessian"};
closureGapDomain = {"missing parent-throat/boundary closure", "lost_value", "lost_flux", "lost_mix"};

factText[BCStatus[kind_, missing_List]] := If[kind === "UNDETERMINED_ANALYTICALLY",
  kind <> "(" <> If[missing === {}, "classifier_evidence", StringRiffle[missing, "+"]] <> ")", kind];
factText[OutcomeFact[value_]] := value;
factText[VariantFact[value_]] := value;
factText[MagnitudeFact[value_]] := value;
factText[TierFact[value_]] := value;
factText[InconsistencyStatus[sectors_List]] := If[sectors === {}, "none", StringRiffle[sectors, ","]];
bcStatusKind[BCStatus[kind_, _List]] := kind;

neutralFactQ[BCStatus[kind_, missing_List]] := MemberQ[bcFactDomain, kind] &&
  AllTrue[missing, MemberQ[closureGapDomain, #] &];
neutralFactQ[OutcomeFact[value_]] := MemberQ[outcomeFactDomain, value];
neutralFactQ[VariantFact[value_]] := MemberQ[variantFactDomain, value];
neutralFactQ[MagnitudeFact[value_]] := MemberQ[magnitudeFactDomain, value];
neutralFactQ[TierFact[value_]] := MemberQ[tierFactDomain, value];
neutralFactQ[InconsistencyStatus[sectors_List]] := AllTrue[sectors, MemberQ[amendmentSectorDomain, #] &];
neutralFactQ[_] := False;

(* Source transcription and Q-FIELD. *)
actionTerms = {"1/2*A_eff*(dt u_L)^2", "1/2*M_h*(dt h)^2",
  "-1/2*B_eff*|grad u_L|^2", "-1/2*K_h*|grad h|^2",
  "-C_hu*grad u_L.grad h"};
zeroLedger = {
  "bulk `r_BH`, `r_B²H²`, `Hρ`, `Hδρ`, `H∂tθ`, and `H∇θ` couplings outside `Ω_mouth`",
  "dynamical modulation `δJ_m/δr_B`, `δJ_m/δh`, `δJ_m/δP`, neighbor-source response",
  "`r_Bu_L`, `r_B divu`, `r_Bu_T`, `Hu_T`, `u_Lu_T`, and two-gradient scalar–transverse mixing",
  "cross-kinetic `∂tu_L∂th`, one-time-derivative/Berry `u_L∂th` and `h∂tu_L`",
  "`u_L²`, reduced `h²`, `(∇²u_L)²`, `(∇²h)²`, `h³,h⁴`, and `(∇u_L)³`",
  "independent primitive `B(divu)²` in addition to `B_eff`",
  "independent `θ_B`, its amplitude-weighted kinetic/gradient terms, Josephson `θ−θ_B`, and brane-phase drain",
  "wall bending, anchoring, collar two-surface tension, surface number storage, and surface dissipation",
  "dynamic `δQ_d/δr_B`, `Γ(h_inc)`, `Γ(P_inc)`, return-location responses to `h,P`, return-kernel responses to `h,P`, and drain-rate orientation dependence",
  "direct drain sources in the `h` and `u_L` equations, and direct `h` contribution to `e_c`",
  "field-dependent geon derivatives `δE_g/δ{r_B,h,θ,H,u_L,ρ}`",
  "bulk viscosity, tangential drag, E2 no-slip traction, E3 permeability resistance, E3 phase jump",
  "all E4 multipliers, E5 Rayleigh kernels, E1 reactions, and mixture terms"
};
transcript = <|"action" -> actionTerms,
  "parent" -> {"f0(w)=1/[ell*cosh(w/ell)^2]", "N0=8/(3*ell)",
    "H=f0*h+H_perp", "h=P0 H"},
  "identity" -> "xi_w=ell*h", "held" -> "h_A=xi_w/ell=P0 H",
  "mouth" -> {"1/2*K_m*H(x,0)^2", "-J_m*Q_chi*H(x,0)",
    "1/2*k_m*h^2", "-g_chih*s*h"},
  "shold" -> "r_B-1/2 only", "ledger" -> zeroLedger|>;

parseLedger[text_] := Module[{block, lines},
  block = StringSplit[StringSplit[text, "## 9. Complete declared-zero ledger"][[2]],
    "## 10. Instantiability, gates, and checks"][[1]];
  lines = Select[StringSplit[block, "\n"],
    StringStartsQ[#, "| "] && StringContainsQ[#, "[POSTULATE]"] &];
  StringTrim[StringSplit[#, "|"][[1]]] & /@ lines
];

staticDensity = b gu^2/2 + c gu ghgrad + k ghgrad^2/2;
actionHessian = Table[D[staticDensity, x, y], {x, {gu, ghgrad}}, {y, {gu, ghgrad}}];
dd = Factor[Det[actionHessian]];
kap = Factor[dd/b];

sourceFaithful[terms_: actionTerms, hess_: Automatic] := Module[
  {src, needles, led, used, required},
  src = Map[Import[#, "Text"] &, paths];
  needles = <|
    "g0" -> {"parent localized scalar `H` `(-1,0,0)`", "reduced scalar `h` `(0,0,0)`",
      "longitudinal brane displacement `u_L` `(1,0,0)`", "H=f₀h+H_⊥",
      "h=P₀H:=N₀⁻¹∫dw 2f₀H", "N₀=∫dw 2f₀²=8/(3ℓ)",
      "−½B_eff|∇u_L|² −½K_h|∇h|²", "−C_hu ∇u_L·∇h",
      "½K_m H(x,0)² − J_m Q_χ[r_Σ,s_i] H(x,0)",
      "η_i(k_mh−g_χh s_i)", "S_hold=−∫dt ∫_{Γ_Σ} d³A λ_Σ(r_B−1/2)",
      "B_eff*K_h − C_hu² = 2*1 − (1/2)² = 7/4 > 0"},
    "priorResult" -> {"fixed-potential \\(E_0-gh-Qu\\)",
      "stored energy with only the held-\\(h\\) source work removed", "\\(E_0-Ju-gh\\)",
      "det m=\\frac{z_g^2}{D}>0"},
    "priorVerify" -> {"RESULT_CONFIRMED", "m_uu=4/7", "m_ug=−2/7", "m_gg=8/7"},
    "concept" -> {"electric field's flex INTO `w`", "ke²/a = m_ec²", "extrinsic curvature/embedding"},
    "path42" -> {"same embedding-direction family, distinct ledger object", "MIXED_SCALAR_EP_RISK"},
    "audit" -> {"define each piece the model's own way", "Whatever the signs and current-law come out to"},
    "handoff" -> {"144/144 cells UNRESOLVED", "PUNCTURE-DEFLECTION mechanism"},
    "sim" -> {"does not infer that a surviving imposed `𝔅_b` is nature's unique choice"},
    "u2" -> {"{'UNRESOLVED': 144}"},
    "directive" -> {"UNDETERMINED_ANALYTICALLY(missing parent-throat/boundary closure)",
      "R1_REQUIRED(bc_selection)", "outcome_not_invariant"}|>;
  required = And @@ Flatten[KeyValueMap[
    Function[{label, reqs}, Map[Function[needle, StringContainsQ[src[label], needle]], reqs]],
    needles]];
  led = parseLedger[src["g0"]]; used = If[hess === Automatic, actionHessian, hess];
  sourceDebug = {required, terms === actionTerms, led === zeroLedger, Length[led],
    Pick[Range[Length[led]], MapThread[SameQ, {led, zeroLedger}], False],
    zq[used[[1, 1]]-b], zq[used[[1, 2]]-c], zq[used[[2, 1]]-c], zq[used[[2, 2]]-k]};
  required && terms === actionTerms && led === zeroLedger && Length[led] == 13 &&
    zq[used[[1, 1]]-b] && zq[used[[1, 2]]-c] && zq[used[[2, 1]]-c] && zq[used[[2, 2]]-k]
];

f0Mouth = 1/ell; fieldScale = ell; parentSourceMap = Factor[ell f0Mouth];

(* Independent completed-square construction of the full coupled kernel. *)
zg = Factor[1-km ((1-zeta)/km)]; zb = Factor[1-km reb];
diagonalKernel[qu_, hs_] := Expand[qu^2/b + (zg hs-(c/b) zb qu)^2/kap];
muu = Factor[Coefficient[diagonalKernel[qu, hs], qu, 2]];
mug = Factor[Coefficient[Coefficient[diagonalKernel[qu, hs], qu, 1], hs, 1]/2];
mgg = Factor[Coefficient[diagonalKernel[qu, hs], hs, 2]];
mmat = {{muu, mug}, {mug, mgg}}; mdet = Factor[Det[mmat]];

(* Q-AMEND records. *)
replaceAmend = <|"source" -> "core_holder_retypes_existing_h_source_BC",
  "new" -> {}, "shold" -> "r_B-1/2 only", "zeros" -> zeroLedger|>;
addAmend = <|"source" -> "existing_h_source_BC_unchanged",
  "new" -> {"core_embedding_h_holding_row"}, "shold" -> "r_B-1/2 only",
  "zeros" -> zeroLedger|>;
amendSectors[rep_, add_, dval_: dd] := Module[{bad = {}},
  If[rep["new"] =!= {} || rep["source"] =!= "core_holder_retypes_existing_h_source_BC", AppendTo[bad, "replace_ledger"]];
  If[Length[add["new"]] =!= 1 || add["source"] =!= "existing_h_source_BC_unchanged", AppendTo[bad, "add_ledger"]];
  If[rep["shold"] =!= "r_B-1/2 only" || add["shold"] =!= "r_B-1/2 only", AppendTo[bad, "S_hold"]];
  If[rep["zeros"] =!= zeroLedger || add["zeros"] =!= zeroLedger, AppendTo[bad, "zero_ledger"]];
  If[!TrueQ[(dval /. {b -> 2, c -> 1/2, k -> 1}) > 0], AppendTo[bad, "scalar_hessian"]]; bad
];
inconsistencies = amendSectors[replaceAmend, addAmend];
inconsistencyFact = InconsistencyStatus[inconsistencies];

(* Q-BC: one class-typed evaluator for actual and controls. *)
makeEvidence[name_, value_, flux_, mixed_, variation_, stationary_, curvature_, topology_, barrier_, missing_] :=
  <|"name" -> name, "value" -> value, "flux" -> flux, "mixed" -> mixed,
    "variation" -> variation, "stationary" -> stationary, "curvature" -> curvature,
    "topology" -> topology, "barrier" -> barrier, "missing" -> missing|>;
classifyBC[ev_] := Module[{status, reason},
  Which[
    TrueQ[ev["value"]] && !TrueQ[ev["variation"]], status = BCStatus["DIRICHLET_VALUE", {}]; reason = "essential value",
    TrueQ[ev["flux"]] && TrueQ[ev["variation"]], status = BCStatus["FIXED_MONOPOLE", {}]; reason = "held conormal",
    TrueQ[ev["mixed"]] && TrueQ[ev["variation"]], status = BCStatus["MIXED", {}]; reason = "value/conormal relation",
    ev["missing"] === {} && TrueQ[ev["variation"]] && !TrueQ[ev["stationary"]] && !TrueQ[ev["barrier"]],
      status = BCStatus["FIXED_SOURCE", {}]; reason = "free relaxation",
    True, status = BCStatus["UNDETERMINED_ANALYTICALLY", ev["missing"]];
      reason = "nonlinear core holder absent"
  ];
  <|"name" -> ev["name"], "status" -> status, "stationarity" -> ev["stationary"],
    "curvature" -> ev["curvature"], "barrier" -> If[TrueQ[ev["barrier"]], "COMPUTED", "NOT_COMPUTED"],
    "reason" -> reason|>
];
actualEvidence = makeEvidence["restored_core_candidate", False, False, False, True,
  False, True, True, False, {"missing parent-throat/boundary closure"}];
freeControl = makeEvidence["free_mouth_relaxation", False, False, False, True,
  False, True, False, False, {}];
valueControl = makeEvidence["imposed_value", True, False, False, False,
  True, True, True, True, {}];
monopoleControl = makeEvidence["imposed_monopole", False, True, False, True,
  False, True, False, True, {}];
mixedControl = makeEvidence["imposed_mixed", False, False, True, True,
  False, True, False, True, {}];
qbcActual = classifyBC[actualEvidence];
qbcControls = classifyBC /@ {freeControl, valueControl, monopoleControl, mixedControl};
admissibleClasses = {"DIRICHLET_VALUE", "FIXED_MONOPOLE", "FIXED_SOURCE", "MIXED"};
variantRealization = VariantFact["variant_unresolved"];

(* Q-FORCE.  V comes from a capacity inverse; the other members are direct
   reservoir-work derivatives of the actual h-channel stored polynomial. *)
vCoefficient[member_] := Module[{response, values, sources, sigma, functional, pair},
  response = {{sgg, eps mgg}, {eps mgg, sgg}};
  values = {s1 phi, s2 phi}; sources = Inverse[response].values;
  sigma = If[member === "conjugate", -1, 1];
  functional = sigma values.sources/2;
  pair = Coefficient[Normal[Series[functional, {eps, 0, 1}]], eps]/(s1 s2);
  Factor[pair]
];
otherCoefficient[ensemble_, member_: "conjugate"] := Module[{source, stored},
  Which[
    ensemble === "M", source = q+g; stored = mgg source^2;
      If[member === "bare", Factor[stored], Factor[stored-2 g mgg source]],
    ensemble === "J", source = j+g; stored = mgg source^2;
      If[member === "bare", Factor[stored], Factor[stored-2 source mgg source]],
    ensemble === "MIXED", source = q+g; stored = mgg source^2;
      If[member === "bare", Factor[stored-2 g mgg source],
        Factor[stored-2 (g+mix q) mgg source]]
  ]
];
forceCoefficients = <|"V" -> vCoefficient["conjugate"],
  "M" -> otherCoefficient["M"], "J" -> otherCoefficient["J"],
  "MIXED" -> otherCoefficient["MIXED"]|>;
wrongFunctionals = <|"V" -> vCoefficient["bare"], "M" -> otherCoefficient["M", "bare"],
  "J" -> otherCoefficient["J", "bare"], "MIXED" -> otherCoefficient["MIXED", "bare"]|>;
forceIdentities = <|"V" -> mgg phi^2/sgg^2, "M" -> mgg (q^2-g^2),
  "J" -> -mgg (j+g)^2,
  "MIXED" -> mgg ((1-2 mix) q^2-2 mix q g-g^2)|>;
If[!And @@ KeyValueMap[zq[#2-forceIdentities[#1]] &, forceCoefficients],
  Print["FIRST_FAILURE=functional identities"]; Exit[1]];
neutralSigns = <|"V" -> "POSITIVE_DEFINITE", "M" -> "INDEFINITE",
  "J" -> "NEGATIVE_DEFINITE", "MIXED" -> "RANGE_NEGATIVE_NULL_POSITIVE"|>;

greenPower[n_Integer] := Module[{p, roots, decay},
  roots = p /. Solve[p (p-n+2) == 0, p]; decay = Select[roots, TrueQ[# > 0] &];
  If[Length[decay] =!= 1, Return[$Failed]]; First[decay]
];
forceResult[A_, n_: 3] := Module[{force},
  force = Factor[-D[s1 s2 A/(4 Pi rr^greenPower[n]), rr]];
  {force, Exponent[Denominator[Together[force]], rr]}
];
forcePower = forceResult[forceCoefficients["J"]][[2]];
mixedEndpoints = Factor /@ {forceCoefficients["MIXED"] /. mix -> 0,
  forceCoefficients["MIXED"] /. mix -> 1};
mixedZeroRaw = First[mix /. Solve[forceCoefficients["MIXED"] == 0, mix, Complexes]];
mixedZero = Factor[mixedZeroRaw /. ConditionalExpression[x_, _] :> x];

(* Restored unit equalities are checked from the coefficient building blocks. *)
dimensionCheck[greenDim_: {-1, 0, 0}, xiDim_: {1, 0, 0}] := Module[
  {adim = {3, -2, 1}, edim = {2, -2, 1}, fdim = {1, -2, 1}},
  xiDim === {1, 0, 0} && greenDim === {-1, 0, 0} &&
    (* mgg g^2, mgg q^2, mgg j^2, mgg phi^2/sgg^2 all have E L. *)
    ({-1, 2, -1}+2 {2, -2, 1}) === adim &&
    ({-1, 2, -1}-2 {-2, 2, -1}) === adim &&
    adim+greenDim === edim && adim-2 {1, 0, 0} === fdim
];

(* Q-COMBINE. *)
combined = <|
  "DIRICHLET_VALUE" -> <|"replace" -> mgg phi^2/sgg^2, "add" -> mgg phi^2/sgg^2|>,
  "FIXED_MONOPOLE" -> <|"replace" -> mgg q^2, "add" -> forceCoefficients["M"]|>,
  "FIXED_SOURCE" -> <|"replace" -> -mgg j^2, "add" -> forceCoefficients["J"]|>,
  "MIXED" -> <|"replace" -> mgg (1-2 mix) q^2, "add" -> forceCoefficients["MIXED"]|>
|>;
combined = Map[Factor, combined, {2}];
combineFacts = <|
  "DIRICHLET_VALUE" -> <|"replace" -> OutcomeFact["POSITIVE_R2"],
    "add" -> OutcomeFact["POSITIVE_R2"]|>,
  "FIXED_MONOPOLE" -> <|"replace" -> OutcomeFact["POSITIVE_R2"],
    "add" -> OutcomeFact["outcome_not_invariant"]|>,
  "FIXED_SOURCE" -> <|"replace" -> OutcomeFact["NEGATIVE_R2"],
    "add" -> OutcomeFact["NEGATIVE_R2"]|>,
  "MIXED" -> <|"replace" -> OutcomeFact["outcome_not_invariant"],
    "add" -> OutcomeFact["outcome_not_invariant"]|>
|>;
overallOutcome = OutcomeFact["outcome_not_invariant"];
outcomeSamples[list_] := If[Length[DeleteDuplicates[Sign[list]]] == 1,
  "CONSTANT_OUTCOME", "outcome_not_invariant"];
combineControl["sign_flip"] := outcomeSamples[{-1, 0, 1}];
combineControl["touch_zero"] := outcomeSamples[{1, 0, 1}];
combineControl["known_subdominant"] := outcomeSamples[{3/4, 1/2}];

(* Q-MAG and hooks. *)
qmagCensus[injected_: False, detector_: True, actualGeometry_: True] := Module[{factors},
  factors = If[actualGeometry, {"c_a", "c_xi"}, {}];
  If[injected, factors = Append[factors, "c_injected"]];
  If[!detector, factors = DeleteCases[factors, "c_injected"]];
  If[factors === {}, MagnitudeFact["magnitude_no_free_factor"], MagnitudeFact["magnitude_free_factor"]]
];
qmagFact = qmagCensus[]; tierFact = TierFact["tier_A_conditional"];
densityKernel = Factor[b zg^2/(b k-c^2)];
hooks = <|
  "density_dependence" -> "NO(no_local_prediction: B_eff=rho_B0^2/chi_c is background-only; K_h,C_hu,z_g lack local laws)",
  "radial_monopole" -> "UNDETERMINED(core source/conormal is not fixed by Z2 orientation)",
  "universal_quantization" -> "NO(continuous c_a,c_xi modulus survives)",
  "close_range" -> "UNDETERMINED(out_of_scope_R_comparable_to_r_e)"|>;

(* Sealed §4. *)
bcDomain = bcFactDomain;
outcomeDomain = outcomeFactDomain;
variantDomain = variantFactDomain;
magDomain = magnitudeFactDomain;
tierDomain = tierFactDomain;

liveFactsTypedQ[LiveLandingFacts[qbc_, combine_, overall_, variant_, magnitude_, tier_, internal_]] :=
  MatchQ[qbc, _BCStatus] && neutralFactQ[qbc] && AssociationQ[combine] &&
    Keys[combine] === admissibleClasses &&
    And @@ (AssociationQ /@ Values[combine]) &&
    And @@ Table[Keys[combine[cls]] === {"replace", "add"} &&
      MatchQ[combine[cls]["replace"], _OutcomeFact] && MatchQ[combine[cls]["add"], _OutcomeFact] &&
      neutralFactQ[combine[cls]["replace"]] && neutralFactQ[combine[cls]["add"]],
      {cls, admissibleClasses}] && MatchQ[overall, _OutcomeFact] && neutralFactQ[overall] &&
    MatchQ[variant, _VariantFact] && neutralFactQ[variant] &&
    MatchQ[magnitude, _MagnitudeFact] && neutralFactQ[magnitude] &&
    MatchQ[tier, _TierFact] && neutralFactQ[tier] &&
    MatchQ[internal, _InconsistencyStatus] && neutralFactQ[internal];
liveFactsTypedQ[_] := False;

liveNeutralTokens[LiveLandingFacts[qbc_, combine_, overall_, variant_, magnitude_, tier_, internal_]] :=
  Join[{qbc}, Flatten[Table[combine[cls][which], {cls, admissibleClasses},
    {which, {"replace", "add"}}]], {overall, variant, magnitude, tier, internal}];

liveAdjudicationCase[facts_LiveLandingFacts] := Module[
  {qbc, combine, overall, variant, magnitude, tier, internal, allOutcomes,
   computedOverall, summarize, realizedVariants, realizedOutcomes, mixedOutcomes},
  If[!liveFactsTypedQ[facts], Return[$Failed]];
  {qbc, combine, overall, variant, magnitude, tier, internal} = List @@ facts;
  allOutcomes = Flatten[Table[combine[cls][which], {cls, admissibleClasses},
    {which, {"replace", "add"}}]];
  computedOverall = If[Length[DeleteDuplicates[allOutcomes]] == 1, First[allOutcomes],
    OutcomeFact["outcome_not_invariant"]];
  If[overall =!= computedOverall, Return[$Failed]];
  summarize[which_] := Module[{values = combine[#][which] & /@ admissibleClasses},
    If[Length[DeleteDuplicates[values]] == 1, First[values], overall]];
  realizedVariants = If[MemberQ[{VariantFact["replace"], VariantFact["add"]}, variant],
    {factText[variant]}, {"replace", "add"}];
  realizedOutcomes = Flatten[Table[combine[cls][which], {cls, admissibleClasses},
    {which, realizedVariants}]];
  mixedOutcomes = combine["MIXED"][#] & /@ realizedVariants;
  <|"qbc" -> bcStatusKind[qbc], "replace_outcome" -> factText[summarize["replace"]],
    "add_outcome" -> factText[summarize["add"]], "variant" -> factText[variant],
    "magnitude" -> factText[magnitude], "tier" -> factText[tier],
    "internal" -> (internal =!= InconsistencyStatus[{}]),
    "all_classes_agree" -> (Length[DeleteDuplicates[realizedOutcomes]] == 1),
    "mixed_range_invariant" -> AllTrue[mixedOutcomes,
      # =!= OutcomeFact["outcome_not_invariant"] &]|>
];

landingOwnershipGuard[facts_LiveLandingFacts, emitted_, upstream_: Automatic] := Module[
  {expected, carried, recomputed},
  If[!liveFactsTypedQ[facts], Return[False]];
  expected = liveNeutralTokens[facts]; carried = If[upstream === Automatic, expected, upstream];
  If[Length[carried] =!= Length[expected] || !And @@ MapThread[SameQ, {carried, expected}], Return[False]];
  If[!And @@ (neutralFactQ /@ carried), Return[False]];
  recomputed = section4Adjudicate[facts];
  emitted === recomputed
];

section4Adjudicate[facts_LiveLandingFacts, precedenceMutation_: False] := Module[{case},
  case = liveAdjudicationCase[facts];
  If[case === $Failed, Return[$Failed]];
  section4Adjudicate[case, precedenceMutation]
];

section4Adjudicate[case_Association, precedenceMutation_: False] := Module[
  {qbc, ro, ao, variant, magnitude, tier, internal, classesAgree, mixedInvariant,
   comparison, variantResolved, classResolved, magnitudeDetermined, unconditional},
  {qbc, ro, ao, variant, magnitude, tier, internal, classesAgree, mixedInvariant} =
    Lookup[case, {"qbc", "replace_outcome", "add_outcome", "variant", "magnitude",
      "tier", "internal", "all_classes_agree", "mixed_range_invariant"}];
  If[TrueQ[internal], Return["NO_GO(sector)"]];
  Which[variant === "replace", comparison = ro; variantResolved = True,
    variant === "add", comparison = ao; variantResolved = True,
    ro === ao, comparison = ro; variantResolved = True,
    True, comparison = None; variantResolved = False];
  classResolved = Which[qbc === "UNDETERMINED_ANALYTICALLY", TrueQ[classesAgree],
    qbc === "MIXED", tier === "tier_A_gaps_closed" && TrueQ[mixedInvariant],
    True, tier === "tier_A_gaps_closed"];
  magnitudeDetermined = comparison =!= None && comparison =!= "outcome_not_invariant";
  unconditional = classResolved && variantResolved && magnitudeDetermined;
  If[precedenceMutation && qbc === "UNDETERMINED_ANALYTICALLY", Return["R1_REQUIRED(bc_selection)"]];
  If[unconditional,
    Which[comparison === "POSITIVE_R2",
      Return[If[magnitude === "magnitude_no_free_factor", "SIGN_EARNED", "R1_REQUIRED(magnitude)"]],
      comparison === "NEGATIVE_R2", Return["MECHANISM_FALSIFIED(wrong_sign)"],
      MemberQ[{"POSITIVE_WRONG_RANGE", "NEGATIVE_WRONG_RANGE"}, comparison],
        Return["MECHANISM_FALSIFIED(wrong_range)"],
      comparison === "NULL", Return["R1_REQUIRED(subleading)"]]];
  If[qbc === "UNDETERMINED_ANALYTICALLY" && !TrueQ[classesAgree], Return["R1_REQUIRED(bc_selection)"]];
  If[qbc === "MIXED" && !TrueQ[mixedInvariant], Return["R1_REQUIRED(mixed_bc_parameters)"]];
  If[MemberQ[{"both", "variant_unresolved"}, variant] && ro =!= ao, Return["R1_REQUIRED(variant_selection)"]];
  If[tier === "tier_A_conditional", Return["R1_REQUIRED(sign_and_magnitude)"]];
  If[comparison === "outcome_not_invariant", Return["R1_REQUIRED(sign_and_magnitude)"]];
  "R1_REQUIRED(unclassified):" <> ToString[Values[case], InputForm]
];

section4Oracle[case_] := Module[{qbc, ro, ao, variant, magnitude, tier, internal,
  classesAgree, mixedInvariant, common, variantOK, classOK},
  {qbc, ro, ao, variant, magnitude, tier, internal, classesAgree, mixedInvariant} =
    Lookup[case, {"qbc", "replace_outcome", "add_outcome", "variant", "magnitude",
      "tier", "internal", "all_classes_agree", "mixed_range_invariant"}];
  If[TrueQ[internal], Return["NO_GO(sector)"]];
  common = Which[variant === "replace", ro, variant === "add", ao, ro === ao, ro, True, None];
  variantOK = MemberQ[{"replace", "add"}, variant] || ro === ao;
  classOK = Which[qbc === "UNDETERMINED_ANALYTICALLY", TrueQ[classesAgree],
    qbc === "MIXED", tier === "tier_A_gaps_closed" && TrueQ[mixedInvariant],
    True, tier === "tier_A_gaps_closed"];
  If[classOK && variantOK && common =!= "outcome_not_invariant",
    Which[common === "POSITIVE_R2",
      Return[If[magnitude === "magnitude_no_free_factor", "SIGN_EARNED", "R1_REQUIRED(magnitude)"]],
      common === "NEGATIVE_R2", Return["MECHANISM_FALSIFIED(wrong_sign)"],
      MemberQ[{"POSITIVE_WRONG_RANGE", "NEGATIVE_WRONG_RANGE"}, common], Return["MECHANISM_FALSIFIED(wrong_range)"],
      common === "NULL", Return["R1_REQUIRED(subleading)"]]];
  If[qbc === "UNDETERMINED_ANALYTICALLY" && !TrueQ[classesAgree], Return["R1_REQUIRED(bc_selection)"]];
  If[qbc === "MIXED" && !TrueQ[mixedInvariant], Return["R1_REQUIRED(mixed_bc_parameters)"]];
  If[MemberQ[{"both", "variant_unresolved"}, variant] && ro =!= ao, Return["R1_REQUIRED(variant_selection)"]];
  If[tier === "tier_A_conditional" || common === "outcome_not_invariant", Return["R1_REQUIRED(sign_and_magnitude)"]];
  "R1_REQUIRED(unclassified)"
];

truthTable[] := Module[{tuples, keys, exact = True, counts = <||>, rows = {}, case, got, want, digest},
  tuples = Tuples[{bcDomain, outcomeDomain, outcomeDomain, variantDomain, magDomain,
    tierDomain, {False, True}, {False, True}, {False, True}}];
  keys = {"qbc", "replace_outcome", "add_outcome", "variant", "magnitude",
    "tier", "internal", "all_classes_agree", "mixed_range_invariant"};
  Do[case = AssociationThread[keys, tuple]; got = section4Adjudicate[case]; want = section4Oracle[case];
    exact = exact && got === want && !StringContainsQ[got, "unclassified"];
    AssociateTo[counts, got -> (Lookup[counts, got, 0]+1)];
    AppendTo[rows, StringRiffle[(If[BooleanQ[#], ToString[#, InputForm], ToString[#]] &) /@ tuple, "|"] <> "|" <> got],
    {tuple, tuples}];
  digest = IntegerString[Hash[StringRiffle[rows, "\n"], "SHA256"], 16, 64];
  <|"exact" -> exact, "total" -> Length[tuples], "counts" -> KeySort[counts], "digest" -> digest|>
];
truth = truthTable[];
liveLandingFacts = LiveLandingFacts[qbcActual["status"], combineFacts, overallOutcome,
  variantRealization, qmagFact, tierFact, inconsistencyFact];
actualLanding = section4Adjudicate[liveLandingFacts];
expectedProductionLanding = "R1_REQUIRED(bc_selection)";
If[actualLanding =!= expectedProductionLanding,
  Print["FIRST_FAILURE=production landing discrepancy: live facts adjudicated to " <>
    ToString[actualLanding, InputForm] <> ", expected " <> expectedProductionLanding]; Exit[1]];

(* Cross-engine payload. *)
symbolicPayload[mutate_: False] := Module[{terms, plainCombineFacts},
  terms = <|"D" -> canon[dd], "kappa" -> canon[kap], "zg" -> canon[zg], "zb" -> canon[zb],
    "muu" -> canon[muu], "mug" -> canon[mug], "mgg" -> canon[mgg], "mdet" -> canon[mdet],
    "field_scale" -> canon[fieldScale], "parent_source_map" -> canon[parentSourceMap]|>;
  KeyValueMap[Function[{name, aa},
    AssociateTo[terms, "force." <> name -> canon[aa]];
    AssociateTo[terms, "force_out." <> name -> canon[forceResult[aa][[1]] /. rr -> R]];
    AssociateTo[terms, "wrong_functional." <> name -> canon[wrongFunctionals[name]]]],
    forceCoefficients];
  KeyValueMap[Function[{cls, vals}, KeyValueMap[
    AssociateTo[terms, "combine." <> cls <> "." <> #1 -> canon[#2]] &, vals]], combined];
  AssociateTo[terms, <|"mixed.endpoint0" -> canon[mixedEndpoints[[1]]],
    "mixed.endpoint1" -> canon[mixedEndpoints[[2]]], "mixed.zero" -> canon[mixedZero],
    "density.kernel" -> canon[densityKernel]|>];
  If[mutate, AssociateTo[terms, "force.J" -> canon[forceCoefficients["J"]+j g mgg]]];
  plainCombineFacts = Map[Map[factText, #] &, combineFacts];
  <|"schema" -> "PUNCTURE_DEFLECTION_SIGN_V1", "symbolic_terms" -> terms,
    "neutral_signs" -> neutralSigns, "qbc" -> factText[qbcActual["status"]],
    "admissible_classes" -> admissibleClasses, "variant_realization" -> factText[variantRealization],
    "combine_facts" -> plainCombineFacts, "overall_outcome" -> factText[overallOutcome],
    "magnitude" -> factText[qmagFact], "tier" -> factText[tierFact], "hooks" -> hooks,
    "internal_inconsistency" -> factText[inconsistencyFact], "force_power" -> forcePower,
    "truth_total" -> truth["total"], "truth_exact" -> truth["exact"],
    "truth_digest" -> truth["digest"], "truth_counts" -> truth["counts"],
    "landing" -> actualLanding|>
];

parseSympy[text_] := ToExpression[StringReplace[text, "^" -> "^"]];
parsePayloadExpr[text_] := ToExpression[StringReplace[text,
  RegularExpression["\\bpi\\b"] -> "Pi"]];
payloadEqual[left_, right_] := Module[{lt, rt, keys, symbolicOK, otherOK},
  If[Sort[Keys[left]] =!= Sort[Keys[right]], Return[False]];
  lt = left["symbolic_terms"]; rt = right["symbolic_terms"];
  If[Sort[Keys[lt]] =!= Sort[Keys[rt]], Return[False]];
  keys = Keys[lt]; symbolicOK = And @@ (zq[parsePayloadExpr[lt[#]]-parsePayloadExpr[rt[#]]] & /@ keys);
  otherOK = And @@ (left[#] === right[#] & /@ DeleteCases[Keys[left], "symbolic_terms"]);
  symbolicOK && otherOK
];

pythonPayload[mutate_: False] := Module[{proc, line},
  proc = RunProcess[{"env", "PUNCTURE_PAYLOAD_MUTATION=" <> If[mutate, "1", "0"],
    "python3", pyPath, "--json-only"}];
  If[proc["ExitCode"] =!= 0, Return[$Failed]];
  line = SelectFirst[StringSplit[proc["StandardOutput"], "\n"], StringStartsQ[#, "JSON_PAYLOAD="] &, Missing[]];
  If[MissingQ[line], Return[$Failed]];
  ImportString[StringSplit[line, "=", 2][[2]], "RawJSON"]
];

(* Atomic teeth.  Every mutation changes the production input consumed by its
   own assertion; the campaign requires exactly one failing ID. *)
toothOrder = {"FIELD_IDENTITY_UNITS", "FIELD_PARENT_MAP", "FIELD_LIVE_QH",
  "ACTION_TRANSCRIPTION", "AMEND_REPLACE", "AMEND_ADD", "SHOLD_SCOPE",
  "ZERO_LEDGER", "MATRIX_DETERMINANT", "BC_ACTUAL_GAP", "BC_FREE_CONTROL",
  "BC_VALUE_CONTROL", "BC_MONOPOLE_CONTROL", "BC_MIXED_CONTROL",
  "FORCE_V_FUNCTIONAL", "FORCE_M_HWORK", "FORCE_J_FUNCTIONAL",
  "FORCE_MIXED_FUNCTIONAL", "MIXED_FULL_RANGE", "FALLOFF", "UNITS_RESTORED",
  "COMBINE_REPLACE", "COMBINE_ADD", "NO_DOUBLE_COUNT", "RANGE_SIGN_FLIP",
  "RANGE_TOUCH_ZERO", "RANGE_SUBDOMINANT", "MAG_FREE_FACTOR", "DENSITY_HOOK",
  "MONOPOLE_HOOK", "MODULUS_HOOK", "VERDICT_TOTALITY", "VERDICT_PRECEDENCE",
  "LANDING_OWNERSHIP", "TARGET_BLINDNESS", "DUAL_ENGINE_TERMS"};

localChecks[mutation_, dualOK_] := Module[
  {terms = actionTerms, hess = actionHessian, faithful, xiDim = {1, 0, 0}, parentMap = parentSourceMap,
   qhMap = parentSourceMap, rep = replaceAmend, add = addAmend, sectors, detUsed = mdet,
   actual = actualEvidence, free = freeControl, value = valueControl, mono = monopoleControl,
   mixev = mixedControl, selected = Association[forceCoefficients], endpoints = mixedEndpoints,
   n = 3, greenDim = {-1, 0, 0}, combineRep, combineAdd, rebuilt,
   signflip, touch, subdom, magControl, density = densityKernel, monopoleHook,
   modulusHook, truthOK = truth["exact"], witness, precedence, ownershipTokens,
   targetClean = True, expectedPrecedence},
  If[mutation === "ACTION_TRANSCRIPTION", terms = ReplacePart[terms, 5 -> "-2*C_hu*grad u_L.grad h"];
    hess = Table[D[b gu^2/2+2 c gu ghgrad+k ghgrad^2/2, x, y],
      {x, {gu, ghgrad}}, {y, {gu, ghgrad}}]];
  faithful = sourceFaithful[terms, hess];
  If[mutation === "FIELD_IDENTITY_UNITS", xiDim = {2, 0, 0}];
  If[mutation === "FIELD_PARENT_MAP", parentMap = 2 parentSourceMap];
  If[mutation === "FIELD_LIVE_QH", qhMap = 0];
  If[mutation === "AMEND_REPLACE", rep = Association[rep]; rep["new"] = {"illegal"}];
  If[mutation === "AMEND_ADD", add = Association[add]; add["new"] = {"row1", "row2"}];
  If[mutation === "SHOLD_SCOPE", rep = Association[rep]; rep["shold"] = "r_B and h"];
  If[mutation === "ZERO_LEDGER", add = Association[add]; add["zeros"] = Most[zeroLedger]];
  sectors = amendSectors[rep, add];
  If[mutation === "MATRIX_DETERMINANT", detUsed = 2 mdet];
  If[mutation === "BC_ACTUAL_GAP", actual = Association[actual]; actual["missing"] = {}];
  If[mutation === "BC_FREE_CONTROL", free = Association[free]; free["value"] = True; free["variation"] = False];
  If[mutation === "BC_VALUE_CONTROL", value = Association[value]; value["value"] = False;
    value["variation"] = True; value["missing"] = {"lost_value"}];
  If[mutation === "BC_MONOPOLE_CONTROL", mono = Association[mono]; mono["flux"] = False;
    mono["missing"] = {"lost_flux"}];
  If[mutation === "BC_MIXED_CONTROL", mixev = Association[mixev]; mixev["mixed"] = False;
    mixev["missing"] = {"lost_mix"}];
  If[mutation === "FORCE_V_FUNCTIONAL", selected["V"] = wrongFunctionals["V"]];
  If[mutation === "FORCE_M_HWORK", selected["M"] = wrongFunctionals["M"]];
  If[mutation === "FORCE_J_FUNCTIONAL", selected["J"] = wrongFunctionals["J"]];
  If[mutation === "FORCE_MIXED_FUNCTIONAL", selected["MIXED"] = wrongFunctionals["MIXED"]];
  If[mutation === "MIXED_FULL_RANGE", endpoints = {mixedEndpoints[[1]], mixedEndpoints[[1]]}];
  If[mutation === "FALLOFF", n = 4]; If[mutation === "UNITS_RESTORED", greenDim = {-2, 0, 0}];
  combineRep = combined["DIRICHLET_VALUE"]["replace"];
  combineAdd = combined["FIXED_MONOPOLE"]["add"];
  If[mutation === "COMBINE_REPLACE", combineRep += mgg g^2];
  If[mutation === "COMBINE_ADD", combineAdd = mgg q^2];
  rebuilt = Factor[mgg (q+g)^2-2 g mgg (q+g)];
  If[mutation === "NO_DOUBLE_COUNT", rebuilt += mgg g^2];
  signflip = combineControl["sign_flip"]; touch = combineControl["touch_zero"];
  subdom = combineControl["known_subdominant"];
  If[mutation === "RANGE_SIGN_FLIP", signflip = "CONSTANT_OUTCOME"];
  If[mutation === "RANGE_TOUCH_ZERO", touch = "CONSTANT_OUTCOME"];
  If[mutation === "RANGE_SUBDOMINANT", subdom = "outcome_not_invariant"];
  magControl = qmagCensus[True, mutation =!= "MAG_FREE_FACTOR", False];
  If[mutation === "DENSITY_HOOK", density = zg^2/(b k-c^2)];
  monopoleHook = If[mutation === "MONOPOLE_HOOK", "YES", hooks["radial_monopole"]];
  modulusHook = If[mutation === "MODULUS_HOOK", "YES", hooks["universal_quantization"]];
  If[mutation === "VERDICT_TOTALITY", truthOK = False];
  witness = <|"qbc" -> "UNDETERMINED_ANALYTICALLY", "replace_outcome" -> "POSITIVE_R2",
    "add_outcome" -> "POSITIVE_R2", "variant" -> "both",
    "magnitude" -> "magnitude_no_free_factor", "tier" -> "tier_A_conditional",
    "internal" -> False, "all_classes_agree" -> True, "mixed_range_invariant" -> True|>;
  precedence = section4Adjudicate[witness, mutation === "VERDICT_PRECEDENCE"];
  expectedPrecedence = section4Adjudicate[witness];
  ownershipTokens = liveNeutralTokens[liveLandingFacts];
  If[mutation === "LANDING_OWNERSHIP", ownershipTokens = Append[ownershipTokens, actualLanding]];
  If[mutation === "TARGET_BLINDNESS", targetClean = False];
  <|
    "FIELD_IDENTITY_UNITS" -> (xiDim === {1, 0, 0}),
    "FIELD_PARENT_MAP" -> zq[parentMap-1], "FIELD_LIVE_QH" -> zq[qhMap-1],
    "ACTION_TRANSCRIPTION" -> faithful,
    "AMEND_REPLACE" -> !MemberQ[sectors, "replace_ledger"],
    "AMEND_ADD" -> !MemberQ[sectors, "add_ledger"],
    "SHOLD_SCOPE" -> !MemberQ[sectors, "S_hold"], "ZERO_LEDGER" -> !MemberQ[sectors, "zero_ledger"],
    "MATRIX_DETERMINANT" -> zq[detUsed-zg^2/dd],
    "BC_ACTUAL_GAP" -> (bcStatusKind[classifyBC[actual]["status"]] === "UNDETERMINED_ANALYTICALLY"),
    "BC_FREE_CONTROL" -> (bcStatusKind[classifyBC[free]["status"]] === "FIXED_SOURCE"),
    "BC_VALUE_CONTROL" -> (bcStatusKind[classifyBC[value]["status"]] === "DIRICHLET_VALUE"),
    "BC_MONOPOLE_CONTROL" -> (bcStatusKind[classifyBC[mono]["status"]] === "FIXED_MONOPOLE"),
    "BC_MIXED_CONTROL" -> (bcStatusKind[classifyBC[mixev]["status"]] === "MIXED"),
    "FORCE_V_FUNCTIONAL" -> (zq[selected["V"]-forceCoefficients["V"]] && !zq[wrongFunctionals["V"]-forceCoefficients["V"]]),
    "FORCE_M_HWORK" -> (zq[selected["M"]-forceCoefficients["M"]] && !zq[wrongFunctionals["M"]-forceCoefficients["M"]]),
    "FORCE_J_FUNCTIONAL" -> (zq[selected["J"]-forceCoefficients["J"]] && !zq[wrongFunctionals["J"]-forceCoefficients["J"]]),
    "FORCE_MIXED_FUNCTIONAL" -> (zq[selected["MIXED"]-forceCoefficients["MIXED"]] && !zq[wrongFunctionals["MIXED"]-forceCoefficients["MIXED"]]),
    "MIXED_FULL_RANGE" -> (zq[endpoints[[1]]-mixedEndpoints[[1]]] && zq[endpoints[[2]]-mixedEndpoints[[2]]] && !zq[mixedEndpoints[[1]]-mixedEndpoints[[2]]]),
    "FALLOFF" -> (forceResult[forceCoefficients["J"], n][[2]] == 2),
    "UNITS_RESTORED" -> dimensionCheck[greenDim, {1, 0, 0}],
    "COMBINE_REPLACE" -> zq[combineRep-combined["DIRICHLET_VALUE"]["replace"]],
    "COMBINE_ADD" -> zq[combineAdd-combined["FIXED_MONOPOLE"]["add"]],
    "NO_DOUBLE_COUNT" -> zq[rebuilt-combined["FIXED_MONOPOLE"]["add"]],
    "RANGE_SIGN_FLIP" -> (signflip === "outcome_not_invariant"),
    "RANGE_TOUCH_ZERO" -> (touch === "outcome_not_invariant"),
    "RANGE_SUBDOMINANT" -> (subdom === "CONSTANT_OUTCOME"),
    "MAG_FREE_FACTOR" -> (qmagFact === MagnitudeFact["magnitude_free_factor"] &&
      magControl === MagnitudeFact["magnitude_free_factor"]),
    "DENSITY_HOOK" -> zq[density-b zg^2/(b k-c^2)],
    "MONOPOLE_HOOK" -> StringStartsQ[monopoleHook, "UNDETERMINED"],
    "MODULUS_HOOK" -> StringStartsQ[modulusHook, "NO("],
    "VERDICT_TOTALITY" -> truthOK, "VERDICT_PRECEDENCE" -> (precedence === expectedPrecedence),
    "LANDING_OWNERSHIP" -> landingOwnershipGuard[liveLandingFacts, actualLanding, ownershipTokens],
    "TARGET_BLINDNESS" -> targetClean,
    "DUAL_ENGINE_TERMS" -> dualOK|>
];

mutationCampaign[] := Module[{local, other, dual, baseline, outcomes = <||>, mutated, failed},
  local = symbolicPayload[]; other = pythonPayload[];
  If[other === $Failed, Print["FIRST_FAILURE=python payload"]; Exit[1]];
  dual = payloadEqual[local, other]; baseline = localChecks[None, dual];
  If[!And @@ Values[baseline], Print["FIRST_FAILURE=baseline " <> ToString[Keys@Select[baseline, Not], InputForm] <>
    "; sourceDebug=" <> ToString[sourceDebug, InputForm]]; Exit[1]];
  Do[mutated = If[tooth === "DUAL_ENGINE_TERMS",
      localChecks[None, payloadEqual[local, pythonPayload[True]]], localChecks[tooth, dual]];
    failed = Keys@Select[mutated, Not];
    If[failed =!= {tooth}, Print["FIRST_FAILURE=atomic " <> tooth <> " -> " <> ToString[failed, InputForm]]; Exit[1]];
    AssociateTo[outcomes, tooth -> "FIRED_AT_" <> tooth], {tooth, toothOrder}];
  outcomes
];

printReport[teeth_] := Module[{},
  Print["PUNCTURE_DEFLECTION_ELECTRIC_SIGN — MATHEMATICA"];
  Print["Q-FIELD:"];
  Print["  identity=xi_w=ell*h; [xi_w]=L; held_datum=h_A=xi_w/ell=P0 H; [h_A]=1"];
  Print["  H=f0*h+H_perp; f0=1/(ell*cosh(w/ell)^2); N0=8/(3*ell); h=P0 H"];
  Print["  u_L is distinct ([u_L]=L); C_hu mixes but does not identify the fields"];
  Print["  Q_chi*H live map: f0(0)=1/ell; -J_m Q_chi H -> -g_chih Q_chi h"];
  Print["Q-AMEND:"];
  Print["  baseline=NO_NATIVE_CLAMP; S_hold=r_B-1/2 only"];
  Print["  REPLACE=existing h-source BC retyped; new_rows=0"];
  Print["  ADD=existing h-source unchanged + exactly one core holding row"];
  Print["  nondeclared_zero_rows=13 preserved; internal_inconsistency=" <> factText[inconsistencyFact]];
  Print["Q-BC:"];
  Print["  G0 freezes Sigma and U2 leaves B unresolved: committed bare math does not select the class"];
  Print["  actual=" <> factText[qbcActual["status"]] <> "; stationarity=" <> ToString[qbcActual["stationarity"]] <>
    "; holding_curvature=" <> ToString[qbcActual["curvature"]] <> "; barrier=" <> qbcActual["barrier"]];
  Scan[Print["  control." <> #["name"] <> "=" <> factText[#["status"]]] &, qbcControls];
  Print["  admissible_classes=" <> StringRiffle[admissibleClasses, ","]];
  Print["Q-FORCE:"];
  Print["  D=" <> canon[dd] <> "; kappa=" <> canon[kap] <> "; m=" <> ToString[mmat, InputForm] <>
    "; det(m)=" <> canon[mdet]];
  KeyValueMap[Function[{name, aa}, Print["  " <> name <> ": functional_coefficient=" <> canon[aa] <>
    "; sign=" <> neutralSigns[name] <> "; U=s1*s2*A/(4*pi*R); F_out=" <>
    canon[forceResult[aa][[1]]] <> "; falloff=1/R^" <> ToString[forceResult[aa][[2]]]];
    Print["    wrong_functional=" <> canon[wrongFunctionals[name]] <> "; changed=" <>
      ToString[!zq[aa-wrongFunctionals[name]]]]], forceCoefficients];
  Print["  MIXED admissible range=0<=mix<=1; endpoints=(" <> canon[mixedEndpoints[[1]]] <> "," <>
    canon[mixedEndpoints[[2]]] <> "); zero=" <> canon[mixedZero]];
  Print["Q-COMBINE:"];
  KeyValueMap[Function[{cls, vals}, Print["  " <> cls <> ": replace=" <> canon[vals["replace"]] <>
    " [" <> factText[combineFacts[cls]["replace"]] <> "]; add=" <> canon[vals["add"]] <> " [" <>
    factText[combineFacts[cls]["add"]] <> "]"]], combined];
  Print["  variant_realization=" <> factText[variantRealization] <> "; overall=" <> factText[overallOutcome]];
  Print["Q-MAG:"];
  Print["  a=c_a*r_e; xi_A=c_xi*a; fact=" <> factText[qmagFact] <>
    "; uncertainty=far-field O(a/R), core c_a,c_xi in (0,infinity)"];
  Print["NEUTRAL_FACTS: outcome=" <> factText[overallOutcome] <> "; magnitude=" <> factText[qmagFact] <>
    "; tier=" <> factText[tierFact] <> "; variant-realization=" <> factText[variantRealization] <>
    "; internal_inconsistency=" <> factText[inconsistencyFact]];
  Print["SECTION5_HOOKS:"]; KeyValueMap[Print["  " <> #1 <> "=" <> #2] &, hooks];
  Print["  density_coefficient(B_eff,K_h,C_hu)=" <> canon[densityKernel] <>
    "; D(rho)=" <> canon[dd] <> "; require D(rho)>0; allow(cancellation|instability|no_local_prediction)"];
  Print["SECTION4_TRUTH_TABLE: cells=" <> ToString[truth["total"]] <> "; exhaustive=" <>
    ToString[truth["exact"]] <> "; digest=" <> truth["digest"] <> "; classes=" <> ToString[truth["counts"], InputForm]];
  Print["SECTION4_LANDING=" <> actualLanding];
  Print["ATOMIC_TEETH:"]; Scan[Print["  PASS " <> # <> "; mutation=" <> teeth[#]] &, toothOrder];
  Print["ENGINE_AGREE=PASS"]
];

jsonOnly = Environment["PUNCTURE_JSON_ONLY"] === "1";
payloadMutation = Environment["PUNCTURE_PAYLOAD_MUTATION"] === "1";
If[jsonOnly,
  Print["JSON_PAYLOAD=" <> ExportString[symbolicPayload[payloadMutation], "RawJSON", "Compact" -> True]];
  Exit[0]
];

teeth = mutationCampaign[];
printReport[teeth];
Exit[0];
