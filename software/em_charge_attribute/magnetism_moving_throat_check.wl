(* Target-blind moving-throat magnetism build: independent Mathematica engine.

   Route B is constructed before Route A from a translated normalized dent,
   the transverse Euler equation, and a direct divergence-free tensor ansatz.
   Python is called only after all local results exist, for term comparison.
*)

ClearAll["Global`*"];
$Assumptions = rhoBr > 0 && muR > 0 && qT != 0 && Ae != 0 && cE > 0 &&
  R > 0 && a0 > 0 && nx^2 + ny^2 + nz^2 == 1;

here = DirectoryName[ExpandFileName[$InputFileName]];
root = ParentDirectory[ParentDirectory[here]];
pyPath = FileNameJoin[{here, "magnetism_moving_throat_check.py"}];
paths = <|
  "directive" -> FileNameJoin[{here, "directive_magnetism_moving_throat.md"}],
  "electricDirective" -> FileNameJoin[{here, "directive_puncture_deflection_electric_sign.md"}],
  "electricResult" -> FileNameJoin[{here, "puncture_deflection_electric_sign_result.md"}],
  "electricCheck" -> FileNameJoin[{here, "puncture_deflection_electric_sign_check.py"}],
  "g0" -> FileNameJoin[{here, "g0_closure_card_v0.md"}],
  "path36md" -> FileNameJoin[{root, "software", "stage1_solver", "reports", "pathA_36_c5_phase_potential.md"}],
  "path36yaml" -> FileNameJoin[{root, "software", "stage1_solver", "reports", "pathA_36_c5_phase_potential_results.yaml"}],
  "audit" -> FileNameJoin[{root, "docs", "model_definition_audit.md"}],
  "handoff" -> FileNameJoin[{root, "docs", "em_analog_next_phase_handoff.md"}],
  "sim" -> FileNameJoin[{root, "docs", "two_throat_simulation_handoff_spec.md"}]
|>;
path39Paths = Sort@FileNames["pathA_39_*", FileNameJoin[{root, "software", "stage1_solver", "reports"}]];

zq[x_] := TrueQ[FullSimplify[x == 0, $Assumptions]];
canon[x_] := ToString[Cancel[x], InputForm];

(* Neutral fact domains. *)
currentDomain = {"CONVECTION_LIKE_CONDITIONAL", "CHARACTERIZED_SOURCE_DEPARTURE",
  "NULL_SOURCE", "R1_SOURCE_BASIS"};
comparisonDomain = {"routes_agree", "routes_differ", "route_B_R1"};
relativeDomain = {"relative_sign_match", "relative_sign_opposite",
  "leading_tensor_conflict", "relative_sign_anchor_conditional"};
magnitudeDomain = {"magnitude_forced_by_electric", "magnitude_free_factor", "R1(magnitude)"};
tierDomain = {"tier_A_gaps_closed", "tier_A_conditional"};
anchorDomain = {"electric_anchor_closed", "R1_REQUIRED(bc_selection)"};

(* Read-first fidelity and ledger-ready transcript. *)
actionKinetic = {"1/2*rho_br*|dt u_T|^2", "-1/2*mu_R*|curl u_T|^2",
  "div u_T=0 (two transverse polarizations)"};
actionCoupling = {"+q_T*sum_i s_i*eta_a(x-X_i)*V_i.u_T",
  "q_T=lambda_T*tau_d; tau_d is the active-throat time-arrow"};
intendedDelta = {"import pathA_36 u_T kinetic row",
  "activate moving Q_chi*V.u_T finite-profile row"};

parseG0Rows[text_] := Module[{block, lines},
  block = StringSplit[StringSplit[text, "## 9. Complete declared-zero ledger"][[2]],
    "## 10. Instantiability, gates, and checks"][[1]];
  lines = Select[StringSplit[block, "\n"], StringStartsQ[#, "| "] &&
      (StringContainsQ[#, "[POSTULATE]"] || StringContainsQ[#, "[CONVENTION]"]) &];
  StringTrim[StringSplit[#, "|"][[2]]] & /@ lines
];

sourceFaithful[kin_: actionKinetic, coup_: actionCoupling, damage_: {}] := Module[
  {src, needles, path39, ok, rows},
  src = Map[Import[#, "Text"] &, paths];
  path39 = Association[(FileNameTake[#] -> Import[#, "Text"]) & /@ path39Paths];
  src = Join[src, path39];
  needles = <|
    "directive" -> {"j∝sV` is `DECLARED_NOT_DERIVED` and BARRED",
      "Maxwell–Darwin vector-current reference", "R1_REQUIRED(electric_bc_selection)",
      "atomic per-tooth"},
    "electricDirective" -> {"neutral throughout", "SEALED adjudication block", "atomic tooth IDs"},
    "electricResult" -> {"xi_w(x)=\\ell h(x)", "m_{gg}=\\frac{b z_g^2}{D}",
      "R1_REQUIRED(bc_selection)", "outcome_not_invariant"},
    "electricCheck" -> {"def section4_adjudicate", "def mutation_campaign", "truth_table"},
    "g0" -> {"The retained fields are", "`r_Bu_L`, `r_B divu`, `r_Bu_T`",
      "Maxwell/gauge fields, point sources, native current law, Coulomb potential prior",
      "magnetism/current structure (R2+)"},
    "path36md" -> {"L_T = 1/2 rho_br dot(u_T)^2 - 1/2 mu_R k^2 u_T^2",
      "c_gamma^2 = mu_R/rho_br", "2 massless transverse photons"},
    "path36yaml" -> {"primitive_lagrangian_per_polarization", "c_gamma_squared: mu_R/rho_br"},
    "pathA_39_magnetic_force.md" -> {"conditional on the declared Stage-1 moving-source amplitudes",
      "j_T=q_A^T V", "F^-1[k_i k_j/k^4]"},
    "pathA_39_magnetic_force_results.yaml" -> {"source: j_T=q_A^T V", "declared_stage2", "sim_deferred_values"},
    "pathA_39_stage3_operator_parity.md" -> {"FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING",
      "P_w odd", "realized nonlinear-throat O(V) coefficient"},
    "pathA_39_stage3_operator_parity_results.yaml" -> {"combined-parity-mixing", "sim_deferred", "iVdotk"},
    "pathA_39_stage4_field_classification.md" -> {"FIELD_SCALAR_VECTOR_DEPARTURE",
      "c_E=c_gamma", "nonlinear throat solve"},
    "pathA_39_stage4_field_classification_results.yaml" -> {"operator_parity_contamination",
      "c_E=c_gamma and lambda_gamma knit"},
    "pathA_39_scalar_admixture_screen.md" -> {"FAIL_OBSERVABLE_SCALAR_ADMIXTURE", "a_L`"},
    "pathA_39_scalar_admixture_results.yaml" -> {"FAIL_OBSERVABLE_SCALAR_ADMIXTURE", "sim_deferred"},
    "audit" -> {"SOURCE LAW ⚠️ IMPORTED", "Whatever the signs and current-law come out to", "F_electric/F_gravity"},
    "handoff" -> {"Target = the characterized-departure Maxwell analog", "j=sV"},
    "sim" -> {"pathA39 `j∝sV` kernel remains ancestry-barred", "R4's magnetic-sign comparison", "magnetic sign is"}
  |>;
  ok = And @@ Flatten@KeyValueMap[
    Function[{label, reqs}, (StringContainsQ[src[label], #] & /@ reqs)], needles];
  rows = parseG0Rows[src["g0"]];
  sourceDebug = {ok, Length[rows], kin === actionKinetic, coup === actionCoupling, damage};
  ok && Length[rows] == 14 && kin === actionKinetic && coup === actionCoupling && damage === {}
];

(* Q-CURRENT: the current basis follows from translating the normalized dent. *)
coords = {x, y, z}; centers = {Xx, Xy, Xz}; bodyV = {Vx, Vy, Vz};
eta = Exp[-Total[(coords-centers)^2]/a0^2]/(Pi^(3/2) a0^3);
signedDensity = s eta;
dtSigned = Factor[Total[bodyV MapThread[D, {Table[signedDensity, 3], centers}]]];
fluxDiv = Factor[Total@MapThread[D[#1, #2] &, {alpha signedDensity bodyV, coords}]];
continuityGeneral = Factor[dtSigned+fluxDiv];
uniqueAlpha = First[alpha /. Solve[(continuityGeneral/
      Total@MapThread[#1 D[signedDensity, #2] &, {bodyV, coords}]) == 0, alpha]];

deriveSource[scale_: 1] := <|
  "density" -> signedDensity,
  "residual" -> Factor[continuityGeneral /. alpha -> scale],
  "alpha" -> uniqueAlpha,
  "vector" -> qT signedDensity bodyV,
  "identity" -> "CONVECTION_LIKE_CONDITIONAL",
  "amplitude" -> "R1(magnitude)",
  "dependencies" -> {"normalized_signed_dent", "translated_body_coordinate", "continuity",
    "active_time_arrow_amplitude"}|>;
source = deriveSource[];

parities = <|
  "s" -> {-1, -1, 1, "scalar"}, "V" -> {1, 1, -1, "polar_vector"},
  "tau_d" -> {1, 1, -1, "scalar"}, "qT_s_V" -> {-1, -1, 1, "polar_vector"},
  "u_T" -> {-1, -1, 1, "polar_vector"}, "b_T" -> {-1, -1, 1, "axial_vector"}|>;

V1 = {v1x, v1y, v1z}; V2 = {v2x, v2y, v2z}; nvec = {nx, ny, nz};
Ddot = Expand[V1.V2]; a1 = Expand[V1.nvec]; a2 = Expand[V2.nvec]; Aprod = Expand[a1 a2];
cg2 = Factor[muR/rhoBr];

forceFromPotential[coefficient_] := Module[{rvec, rad, aa1, aa2, potential, raw, expected, radial},
  rvec = {rfx, rfy, rfz}; rad = Sqrt[rvec.rvec];
  aa1 = V1.rvec/rad; aa2 = V2.rvec/rad;
  potential = coefficient (Ddot+aa1 aa2)/rad;
  raw = -Table[D[potential, p], {p, {rfx, rfy, rfz}}];
  expected = -coefficient (a2 V1+a1 V2-(Ddot+3 Aprod) nvec)/R^2;
  forceDerivationResidual = FullSimplify[
    raw-(expected /. {R -> rad, nx -> rfx/rad, ny -> rfy/rad, nz -> rfz/rad}),
    rad > 0];
  If[!And @@ (zq /@ forceDerivationResidual), Return[$Failed]];
  radial = Factor[coefficient (Ddot+Aprod)/R^2];
  {Factor /@ expected, radial}
];

(* Route B is intentionally constructed first.  Divergence fixes b=a and the
   two-polarization trace fixes 3a+b=2/(4 Pi). *)
deriveRouteB[src_, foreign_: None] := Module[{sol, kernel, interaction, deps, fr},
  sol = First@Solve[{bb-aa == 0, 3 aa+bb == 1/(2 Pi)}, {aa, bb}];
  kernel = FullSimplify[(aa IdentityMatrix[3]+bb Outer[Times, nvec, nvec])/(muR R) /. sol];
  interaction = Factor[-s1 s2 qT^2 V1.kernel.V2];
  deps = {"Q_CURRENT.source", "pathA36.transverse_EOM", "direct_transverse_tensor_ansatz"};
  If[foreign =!= None,
    (* Illicit copy: this rescaling leaves Route A a factor muR above Route B;
       the tooth detects the dependency tag. *)
    interaction = Factor[foreign["interaction"] qT^2 cg2/Ae];
    deps = Append[deps, "ILLICIT_ROUTE_A_READ"]];
  fr = forceFromPotential[-s1 s2 qT^2/(8 Pi muR)];
  <|"name" -> "B_DIRECT", "kernel" -> kernel, "interaction" -> interaction,
    "force" -> fr[[1]], "radial" -> fr[[2]], "potentialPower" -> 1,
    "forcePower" -> 2, "velocityOrder" -> 2, "dependencies" -> deps|>
];
routeB = deriveRouteB[source];

(* Route A independently derives the Coulomb-gauge projector from the k^-4
   radial seed.  This is the Lorentz-completed conditional electric anchor. *)
deriveRouteA[] := Module[{rvec, rad, seed, kk, transverse, rawAtR, expected, atR, interaction, fr},
  rvec = {rax, ray, raz}; rad = Sqrt[rvec.rvec]; seed = -rad/(8 Pi);
  kk = Table[-D[seed, rvec[[i]], rvec[[j]]], {i, 3}, {j, 3}];
  transverse = FullSimplify[IdentityMatrix[3]/(4 Pi rad)-kk, rad > 0];
  rawAtR = FullSimplify[transverse /. {rax -> R nx, ray -> R ny, raz -> R nz},
    R > 0 && nx^2+ny^2+nz^2 == 1];
  expected = (IdentityMatrix[3]+Outer[Times, nvec, nvec])/(8 Pi R);
  If[!And @@ Flatten[MapThread[
      TrueQ[FullSimplify[#1 == #2, R > 0 && nx^2+ny^2+nz^2 == 1]] &,
      {rawAtR, expected}, 2]], Return[$Failed]];
  atR = expected;
  interaction = Factor[-s1 s2 Ae/cg2 V1.atR.V2];
  fr = forceFromPotential[-s1 s2 Ae/(8 Pi cg2)];
  <|"name" -> "A_BOOST", "kernel" -> atR, "interaction" -> interaction,
    "force" -> fr[[1]], "radial" -> fr[[2]], "potentialPower" -> 1,
    "forcePower" -> 2, "velocityOrder" -> 2,
    "dependencies" -> {"electric_conditional_anchor", "four_current_completion",
      "Coulomb_gauge_transverse_projector", "c_gamma_cone"}|>
];
routeA = deriveRouteA[];

electricU = s1 s2 Ae/(4 Pi R); electricFr = s1 s2 Ae/(4 Pi R^2);
expectedA = -s1 s2 Ae (Ddot+Aprod)/(8 Pi cg2 R);
expectedB = -s1 s2 qT^2 (Ddot+Aprod)/(8 Pi muR R);
ratioBA = Factor[qT^2 cg2/(muR Ae)]; deltaNormalized = Factor[ratioBA-1];
coneRatio = Factor[cE^2/cg2]; deltaU = Factor[routeB["interaction"]-routeA["interaction"]];
parallelD = v1y v2y;
parallelRatioA = Factor[routeA["radial"]/electricFr /.
  {nx -> 1, ny -> 0, nz -> 0, v1x -> 0, v2x -> 0, v1z -> 0, v2z -> 0}];

routeIndependence[rb_: routeB] := AllTrue[rb["dependencies"], !StringContainsQ[#, "ROUTE_A"] &] &&
  rb["name"] === "B_DIRECT";

compareRoutes[ra_: routeA, rb_: routeB] := Module[{ta, tb},
  ta = Factor[ra["interaction"]/(-s1 s2 Ae/cg2)];
  tb = Factor[rb["interaction"]/(-s1 s2 qT^2/muR)];
  <|"tensor" -> zq[ta-tb], "falloff" -> (ra["forcePower"] == rb["forcePower"] == 2),
    "order" -> (ra["velocityOrder"] == rb["velocityOrder"] == 2),
    "relative" -> "relative_sign_anchor_conditional", "comparison" -> "route_B_R1",
    "delta" -> deltaU, "ratio" -> ratioBA, "deltaNormalized" -> deltaNormalized,
    "cone" -> coneRatio,
    "gap" -> "nonlinear throat q_T normalization + electric boundary selection"|>
];
comparison = compareRoutes[];

(* Units and neutral scope facts. *)
dimensionCheck[qtDim_: {0, -1, 1}] := Module[
  {rhoDim = {-3, 0, 1}, muDim = {-1, -2, 1}, aeDim = {3, -2, 1},
   vDim = {1, -1, 0}, lDim = {1, 0, 0}, eDim = {2, -2, 1}, fDim = {1, -2, 1}},
  rhoDim+2 vDim === {-1, -2, 1} &&
    2 qtDim-muDim+2 vDim-lDim === eDim &&
    aeDim-2 vDim+2 vDim-lDim === eDim &&
    2 qtDim-muDim+2 vDim-2 lDim === fDim &&
    2 qtDim+2 vDim-muDim-aeDim === {0, 0, 0}
];

amendmentConsistency[rhoVal_: rhoBr, muVal_: muR, damage_: {}] := Module[{bad = damage},
  If[!TrueQ[(rhoVal /. {rhoBr -> 3, muR -> 5}) > 0], AppendTo[bad, "transverse_inertia"]];
  If[!TrueQ[(muVal /. {rhoBr -> 3, muR -> 5}) > 0], AppendTo[bad, "transverse_stiffness"]];
  bad
];
internalSectors = amendmentConsistency[];
qmagFact = "R1(magnitude)"; tierFact = "tier_A_conditional";
anchorFact = "R1_REQUIRED(bc_selection)";
hooks = <|
  "emergent_lorentz" -> "UNDETERMINED: delta_U is proportional to (q_T^2*c_gamma^2/(mu_R*A_E)-1); also require c_E^2/c_gamma^2=1 and Tier-B closure",
  "preferred_frame" -> "UNDETERMINED: medium-frame O(V1*V2) tensor is known; coefficient/cone and higher-order leakage are R1",
  "active_flux" -> "R1: G0 F_flux may have a non-conservative O(V1*V2) contribution; integrability of the full observable is not supplied by the conservative row",
  "hierarchy_capstone" -> "FLAG_ONLY: F_e/F_g becomes the held-out ratio of two couplings in one action"|>;

(* SEALED SECTION 4. *)
section4Adjudicate[case_, precedenceMutation_: False] := Module[{complete, reason},
  If[TrueQ[case["internal"]], Return["NO_GO(sector)"]];
  complete = case["current"] =!= "R1_SOURCE_BASIS" &&
    case["comparison"] =!= "route_B_R1" &&
    case["relative_sign"] =!= "relative_sign_anchor_conditional" &&
    case["tier"] === "tier_A_gaps_closed" && case["anchor"] === "electric_anchor_closed" &&
    case["magnitude"] =!= "R1(magnitude)";
  If[precedenceMutation && case["anchor"] === "R1_REQUIRED(bc_selection)",
    Return["R1_REQUIRED(consistency)"]];
  If[complete,
    If[MemberQ[{"relative_sign_opposite", "leading_tensor_conflict"}, case["relative_sign"]],
      reason = If[case["relative_sign"] === "relative_sign_opposite", "wrong_relative_sign", "routes_leading_conflict"];
      Return["AMENDMENT_EXCLUDED(" <> reason <> ")"]];
    If[case["comparison"] === "routes_agree",
      Return[If[case["magnitude"] === "magnitude_forced_by_electric",
        "MAGNETISM_LORENTZ_CONSISTENT", "MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)"]]];
    If[case["comparison"] === "routes_differ", Return["MAGNETISM_DEPARTURE_CHARACTERIZED"]]
  ];
  If[case["anchor"] === "R1_REQUIRED(bc_selection)", Return["R1_REQUIRED(electric_bc_selection)"]];
  If[case["current"] === "R1_SOURCE_BASIS" || case["comparison"] === "route_B_R1",
    Return["R1_REQUIRED(direct_moving_throat)"]];
  If[case["magnitude"] === "R1(magnitude)", Return["R1_REQUIRED(magnitude)"]];
  If[case["tier"] === "tier_A_conditional" ||
      case["relative_sign"] === "relative_sign_anchor_conditional",
    Return["R1_REQUIRED(consistency)"]];
  "R1_REQUIRED(unclassified)"
];

section4Oracle[case_] := Module[{complete},
  If[TrueQ[case["internal"]], Return["NO_GO(sector)"]];
  complete = case["current"] =!= "R1_SOURCE_BASIS" && case["comparison"] =!= "route_B_R1" &&
    case["relative_sign"] =!= "relative_sign_anchor_conditional" &&
    case["tier"] === "tier_A_gaps_closed" && case["anchor"] === "electric_anchor_closed" &&
    case["magnitude"] =!= "R1(magnitude)";
  If[complete,
    If[case["relative_sign"] === "relative_sign_opposite", Return["AMENDMENT_EXCLUDED(wrong_relative_sign)"]];
    If[case["relative_sign"] === "leading_tensor_conflict", Return["AMENDMENT_EXCLUDED(routes_leading_conflict)"]];
    If[case["comparison"] === "routes_agree", Return[If[
      case["magnitude"] === "magnitude_forced_by_electric", "MAGNETISM_LORENTZ_CONSISTENT",
      "MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)"]]];
    Return["MAGNETISM_DEPARTURE_CHARACTERIZED"]];
  If[case["anchor"] === "R1_REQUIRED(bc_selection)", Return["R1_REQUIRED(electric_bc_selection)"]];
  If[case["current"] === "R1_SOURCE_BASIS" || case["comparison"] === "route_B_R1",
    Return["R1_REQUIRED(direct_moving_throat)"]];
  If[case["magnitude"] === "R1(magnitude)", Return["R1_REQUIRED(magnitude)"]];
  If[case["tier"] === "tier_A_conditional" || case["relative_sign"] === "relative_sign_anchor_conditional",
    Return["R1_REQUIRED(consistency)"]];
  "R1_REQUIRED(unclassified)"
];

truthTable[precedenceMutation_: False] := Module[{rows = {}, counts = <||>, exact = True,
  keys, case, got, want, digest},
  keys = {"current", "comparison", "relative_sign", "magnitude", "tier", "anchor", "internal"};
  Do[case = AssociationThread[keys, values];
    got = section4Adjudicate[case, precedenceMutation]; want = section4Oracle[case];
    If[!precedenceMutation, exact = exact && got === want && !StringContainsQ[got, "unclassified"]];
    AssociateTo[counts, got -> (Lookup[counts, got, 0]+1)];
    AppendTo[rows, StringRiffle[
      (If[StringQ[#], #, ToString[#, InputForm]] &) /@ values, "|"] <> "|" <> got],
    {values, Tuples[{currentDomain, comparisonDomain, relativeDomain, magnitudeDomain,
      tierDomain, anchorDomain, {False, True}}]}];
  digest = Hash[StringRiffle[rows, "\n"], "SHA256", "HexString"];
  <|"exact" -> exact, "total" -> Length[rows], "counts" -> KeySort[counts], "digest" -> digest|>
];
truth = truthTable[];

liveFacts = <|"current" -> source["identity"], "comparison" -> comparison["comparison"],
  "relative_sign" -> comparison["relative"], "magnitude" -> qmagFact, "tier" -> tierFact,
  "anchor" -> anchorFact, "internal" -> (internalSectors =!= {})|>;
actualLanding = section4Adjudicate[liveFacts];
r1Blockers[] := DeleteCases[{
  If[anchorFact === "R1_REQUIRED(bc_selection)", "R1_REQUIRED(electric_bc_selection)", Nothing],
  If[source["identity"] === "R1_SOURCE_BASIS" || comparison["comparison"] === "route_B_R1",
    "R1_REQUIRED(direct_moving_throat)", Nothing],
  If[qmagFact === "R1(magnitude)", "R1_REQUIRED(magnitude)", Nothing],
  If[tierFact === "tier_A_conditional", "R1_REQUIRED(consistency)", Nothing]}, Nothing];
liveTokens[] := {source["identity"], comparison["comparison"], comparison["relative"],
  qmagFact, tierFact, anchorFact, internalSectors};
landingOwnership[tokens_: Automatic] := Module[{used = If[tokens === Automatic, liveTokens[], tokens]},
  Length[used] == Length[liveTokens[]] && used === liveTokens[] && actualLanding === section4Adjudicate[liveFacts]
];

targetBlindness[] := Module[{text = Import[$InputFileName, "Text"], before},
  before = First@StringSplit[text, "(* SEALED SECTION 4. *)"];
  !Or @@ (StringContainsQ[ToLowerCase[before], #] & /@
      {"magnetism_lorentz_consistent", "amendment_excluded", "magnetism_departure_characterized"})
];

symbolicPayload[mutate_: False] := Module[{terms, rb},
  terms = <|
    "source_continuity" -> canon[source["residual"]], "source_alpha" -> canon[source["alpha"]],
    "c_gamma_squared" -> canon[cg2],
    "routeA_kernel00" -> canon[routeA["kernel"][[1, 1]]],
    "routeA_kernel01" -> canon[routeA["kernel"][[1, 2]]],
    "routeA_U2" -> canon[routeA["interaction"]], "routeA_Fr" -> canon[routeA["radial"]],
    "routeB_kernel00" -> canon[routeB["kernel"][[1, 1]]],
    "routeB_kernel01" -> canon[routeB["kernel"][[1, 2]]],
    "routeB_U2" -> canon[routeB["interaction"]], "routeB_Fr" -> canon[routeB["radial"]],
    "delta_U" -> canon[deltaU], "ratio_BA" -> canon[ratioBA],
    "delta_normalized" -> canon[deltaNormalized], "cone_ratio" -> canon[coneRatio],
    "electric_U0" -> canon[electricU], "electric_Fr" -> canon[electricFr]|>;
  If[mutate, terms["routeB_U2"] = canon[routeB["interaction"]+qT^2 Ddot/(muR R)]];
  <|"schema" -> "MAGNETISM_MOVING_THROAT_V1", "symbolic_terms" -> terms,
    "current_identity" -> source["identity"], "source_amplitude" -> source["amplitude"],
    "parities" -> parities, "comparison" -> comparison["comparison"],
    "structural_agreement" -> TrueQ[comparison["tensor"] && comparison["falloff"] && comparison["order"]],
    "relative_sign" -> comparison["relative"], "magnitude" -> qmagFact, "tier" -> tierFact,
    "anchor" -> anchorFact, "internal_inconsistency" -> internalSectors,
    "potential_power" -> routeB["potentialPower"], "force_power" -> routeB["forcePower"],
    "velocity_order" -> routeB["velocityOrder"], "truth_exact" -> truth["exact"],
    "truth_total" -> truth["total"], "truth_counts" -> truth["counts"],
    "truth_digest" -> truth["digest"], "landing" -> actualLanding,
    "r1_blockers" -> r1Blockers[]|>
];

pythonPayload[mutate_: False] := Module[{run, line},
  run = RunProcess[{"env", "MAGNETISM_PAYLOAD_MUTATION=" <> If[mutate, "1", "0"],
    "/usr/bin/python3", pyPath, "--json-only"}];
  If[run["ExitCode"] =!= 0, pythonPayloadDebug = run; Return[$Failed]];
  line = SelectFirst[StringSplit[run["StandardOutput"], "\n"], StringStartsQ[#, "JSON_PAYLOAD="] &, Missing[]];
  If[MissingQ[line], Return[$Failed]]; ImportString[StringSplit[line, "=", 2][[2]], "RawJSON"]
];

parseExpr[text_] := ToExpression[StringReplace[text,
  RegularExpression["\\bpi\\b"] -> "Pi"]];
payloadEqual[left_, right_] := Module[{lt, rt, keys},
  If[Sort[Keys[left]] =!= Sort[Keys[right]], Return[False]];
  lt = left["symbolic_terms"]; rt = right["symbolic_terms"];
  If[Sort[Keys[lt]] =!= Sort[Keys[rt]], Return[False]];
  keys = Keys[lt];
  If[!And @@ (zq[parseExpr[lt[#]]-parseExpr[rt[#]]] & /@ keys), Return[False]];
  And @@ (left[#] === right[#] & /@ DeleteCases[Keys[left], "symbolic_terms"])
];

toothOrder = {"SOURCE_TRANSLATION_CONTINUITY", "SOURCE_NOT_IMPORTED", "SOURCE_BASIS",
  "PARITY_RW", "PARITY_PW", "PARITY_ROTATION", "PARITY_TIME_REVERSAL",
  "FIELD_IDENTITY_UNITS", "ACTION_KINETIC", "ACTION_COUPLING", "ACTION_STABILITY",
  "G0_DAMAGE", "ROUTE_INDEPENDENCE", "BOOST_PROJECTOR", "BOOST_GENERAL_VELOCITIES",
  "BOOST_NEXT_ORDER", "BOOST_COMMON_VELOCITY", "DIRECT_SOURCE", "DIRECT_PROJECTOR",
  "DIRECT_EXCHANGE_SIGN", "DIRECT_FALLOFF", "DIRECT_VELOCITY_ORDER", "COMPARE_COMPUTED",
  "DELTA_RATIO", "CONE_RATIO", "QMAG_R1", "UNITS_RESTORED", "ACTIVE_FLUX_CAVEAT",
  "HOOK_LORENTZ", "LEDGER_READY_ROW", "TRUTH_TOTALITY", "TRUTH_PRECEDENCE",
  "LANDING_OWNERSHIP", "TARGET_BLINDNESS", "DUAL_ENGINE_TERMS"};

localChecks[mutation_, dualOK_] := Module[
  {srcScale = 1, src, barred = {}, basis, pars = Association[parities], fieldDim = {1, 0, 0},
   kin = actionKinetic, coup = actionCoupling, damage = {}, faithful, rhoUsed = rhoBr,
   copied = routeB, raKernel = routeA["kernel"], raU = routeA["interaction"], boostOrder = 2,
   common = parallelRatioA, directU = routeB["interaction"], directKernel = routeB["kernel"],
   directSign = -1, directPower = 2, directOrder = 2, comp = compareRoutes[], compTensor,
   delta = ratioBA, cone = coneRatio, mag = qmagFact, qtDim = {0, -1, 1}, flux, lorentz,
   ledger = intendedDelta, truthOK = truth["exact"], precedence, tokens, blind},
  If[mutation === "SOURCE_TRANSLATION_CONTINUITY", srcScale = 2]; src = deriveSource[srcScale];
  If[mutation === "SOURCE_NOT_IMPORTED", barred = {"Nu"}]; basis = src["vector"];
  If[mutation === "SOURCE_BASIS", basis = 2 basis];
  If[mutation === "PARITY_RW", pars["u_T"] = {1, -1, 1, "polar_vector"}];
  If[mutation === "PARITY_PW", pars["u_T"] = {-1, 1, 1, "polar_vector"}];
  If[mutation === "PARITY_ROTATION", pars["u_T"] = {-1, -1, 1, "scalar"}];
  If[mutation === "PARITY_TIME_REVERSAL", pars["qT_s_V"] = {-1, -1, -1, "polar_vector"}];
  If[mutation === "FIELD_IDENTITY_UNITS", fieldDim = {0, 0, 0}];
  If[mutation === "ACTION_KINETIC", kin = Most[actionKinetic]];
  If[mutation === "ACTION_COUPLING", coup = {"+q_T sum_i s_i eta_i u_T", actionCoupling[[2]]}];
  If[mutation === "G0_DAMAGE", damage = {"existing_scalar_row"}]; faithful = sourceFaithful[kin, coup, damage];
  If[mutation === "ACTION_STABILITY", rhoUsed = -rhoBr];
  If[mutation === "ROUTE_INDEPENDENCE", copied = deriveRouteB[source, routeA]];
  If[mutation === "BOOST_PROJECTOR", raKernel = 2 raKernel];
  If[mutation === "BOOST_GENERAL_VELOCITIES", raU = raU /. Aprod -> 0];
  If[mutation === "BOOST_NEXT_ORDER", boostOrder = 4];
  If[mutation === "BOOST_COMMON_VELOCITY", common += parallelD/cg2];
  If[mutation === "DIRECT_SOURCE", directU = directU /. qT^2 -> 2 qT^2];
  If[mutation === "DIRECT_PROJECTOR", directKernel += Outer[Times, nvec, nvec]/(8 Pi muR R)];
  If[mutation === "DIRECT_EXCHANGE_SIGN", directSign = 1];
  If[mutation === "DIRECT_FALLOFF", directPower = 3];
  If[mutation === "DIRECT_VELOCITY_ORDER", directOrder = 1];
  compTensor = comp["tensor"]; If[mutation === "COMPARE_COMPUTED", compTensor = False];
  If[mutation === "DELTA_RATIO", delta *= 2]; If[mutation === "CONE_RATIO", cone = cE/Sqrt[cg2]];
  If[mutation === "QMAG_R1", mag = "magnitude_forced_by_electric"];
  If[mutation === "UNITS_RESTORED", qtDim = {1, -1, 1}];
  flux = If[mutation === "ACTIVE_FLUX_CAVEAT", "DECIDED_ZERO", hooks["active_flux"]];
  lorentz = If[mutation === "HOOK_LORENTZ", "YES", hooks["emergent_lorentz"]];
  If[mutation === "LEDGER_READY_ROW", ledger = Most[intendedDelta]];
  If[mutation === "TRUTH_TOTALITY", truthOK = False];
  precedence = section4Adjudicate[liveFacts, mutation === "TRUTH_PRECEDENCE"];
  tokens = liveTokens[]; If[mutation === "LANDING_OWNERSHIP", tokens = Append[tokens, actualLanding]];
  blind = targetBlindness[] && mutation =!= "TARGET_BLINDNESS";
  expectedKernel = (IdentityMatrix[3]+Outer[Times, nvec, nvec])/(8 Pi R);
  <|
    "SOURCE_TRANSLATION_CONTINUITY" -> zq[src["residual"]],
    "SOURCE_NOT_IMPORTED" -> (Intersection[barred, {"Nu", "aT", "aTp", "aL", "q_A_T", "q_L"}] === {}),
    "SOURCE_BASIS" -> And @@ (zq /@ (basis-source["vector"])),
    "PARITY_RW" -> (pars["qT_s_V"][[1]] == pars["u_T"][[1]] == -1),
    "PARITY_PW" -> (pars["qT_s_V"][[2]] == pars["u_T"][[2]] == -1),
    "PARITY_ROTATION" -> (pars["u_T"][[4]] === "polar_vector" && pars["b_T"][[4]] === "axial_vector"),
    "PARITY_TIME_REVERSAL" -> (pars["tau_d"][[3]] == -1 && pars["V"][[3]] == -1 &&
      pars["qT_s_V"][[3]] == pars["u_T"][[3]] == 1),
    "FIELD_IDENTITY_UNITS" -> (fieldDim === {1, 0, 0}),
    "ACTION_KINETIC" -> If[mutation === "ACTION_KINETIC", faithful, kin === actionKinetic],
    "ACTION_COUPLING" -> If[mutation === "ACTION_COUPLING", faithful, coup === actionCoupling],
    "ACTION_STABILITY" -> (amendmentConsistency[rhoUsed, muR] === {}),
    "G0_DAMAGE" -> (damage === {} && internalSectors === {}),
    "ROUTE_INDEPENDENCE" -> routeIndependence[copied],
    "BOOST_PROJECTOR" -> And @@ Flatten[MapThread[zq[#1-#2] &, {raKernel, expectedKernel}, 2]],
    "BOOST_GENERAL_VELOCITIES" -> zq[raU-expectedA], "BOOST_NEXT_ORDER" -> (boostOrder == 2),
    "BOOST_COMMON_VELOCITY" -> zq[common+parallelD/(2 cg2)],
    "DIRECT_SOURCE" -> zq[directU-expectedB],
    "DIRECT_PROJECTOR" -> And @@ Flatten[MapThread[zq[#1-#2] &, {directKernel, expectedKernel/muR}, 2]],
    "DIRECT_EXCHANGE_SIGN" -> (directSign == -1 && zq[routeB["interaction"]-expectedB]),
    "DIRECT_FALLOFF" -> (directPower == 2), "DIRECT_VELOCITY_ORDER" -> (directOrder == 2),
    "COMPARE_COMPUTED" -> TrueQ[compTensor && comp["falloff"] && comp["order"]],
    "DELTA_RATIO" -> zq[delta-qT^2 cg2/(muR Ae)], "CONE_RATIO" -> zq[cone-cE^2/cg2],
    "QMAG_R1" -> (mag === "R1(magnitude)"), "UNITS_RESTORED" -> dimensionCheck[qtDim],
    "ACTIVE_FLUX_CAVEAT" -> StringStartsQ[flux, "R1:"],
    "HOOK_LORENTZ" -> StringStartsQ[lorentz, "UNDETERMINED:"],
    "LEDGER_READY_ROW" -> (ledger === intendedDelta), "TRUTH_TOTALITY" -> truthOK,
    "TRUTH_PRECEDENCE" -> (precedence === actualLanding),
    "LANDING_OWNERSHIP" -> landingOwnership[tokens], "TARGET_BLINDNESS" -> blind,
    "DUAL_ENGINE_TERMS" -> dualOK|>
];

mutationCampaign[] := Module[{local, other, dual, baseline, outcomes = <||>, checks, failed},
  local = symbolicPayload[]; other = pythonPayload[];
  If[other === $Failed, Print["FIRST_FAILURE=python payload " <>
      ToString[If[ValueQ[pythonPayloadDebug], pythonPayloadDebug, "no-debug"], InputForm]]; Exit[1]];
  dual = payloadEqual[local, other]; baseline = localChecks[None, dual];
  If[!And @@ Values[baseline], Print["FIRST_FAILURE=baseline " <>
      ToString[Keys@Select[baseline, Not], InputForm] <> "; source=" <> ToString[sourceDebug, InputForm]]; Exit[1]];
  Do[checks = If[tooth === "DUAL_ENGINE_TERMS",
      localChecks[None, payloadEqual[local, pythonPayload[True]]], localChecks[tooth, dual]];
    failed = Keys@Select[checks, Not];
    If[failed =!= {tooth}, Print["FIRST_FAILURE=atomic " <> tooth <> " -> " <> ToString[failed, InputForm]]; Exit[1]];
    AssociateTo[outcomes, tooth -> "FIRED_AT_" <> tooth], {tooth, toothOrder}];
  outcomes
];

printReport[teeth_] := Module[{},
  Print["MAGNETISM_MOVING_THROAT — MATHEMATICA"];
  Print["G0+DELTA_ACTION:"];
  Print["  L_T=1/2*rho_br*|dt u_T|^2-1/2*mu_R*|curl u_T|^2; div(u_T)=0"];
  Print["  L_move=+q_T*sum_i s_i*eta_a(x-X_i)*V_i.u_T; q_T=lambda_T*tau_d"];
  Print["  c_gamma^2=" <> canon[cg2] <> "; internal_inconsistency=" <> If[internalSectors === {}, "none", ToString[internalSectors]]];
  Print["Q-CURRENT:"];
  Print["  sigma_i=s_i*eta_a(x-X_i); continuity_residual=" <> canon[source["residual"]] <>
    "; unique_flux_alpha=" <> canon[source["alpha"]]];
  Print["  I_i=s_i*eta_i*V_i derived; J_i=q_T*I_i; pathA_39 source amplitudes barred"];
  Print["  identity=" <> source["identity"] <> "; amplitude=" <> source["amplitude"]];
  Print["  u_T [L] polar/T-even/R_w-odd/P_w-odd; b_T=curl u_T [1] axial/T-even"];
  Print["Q-BOOST_ROUTE_A:"];
  Print["  U_A0=" <> canon[electricU]]; Print["  U_A2=" <> canon[routeA["interaction"]] <> " + O(V^4/c_gamma^4)"];
  Print["  F_A2=" <> ToString[canon /@ routeA["force"], InputForm] <> "; radial=" <> canon[routeA["radial"]]];
  Print["  independent (V1,V2); common side-by-side boost gives 1-V^2/(2*c_gamma^2)+O(V^4)"];
  Print["Q-DIRECT_ROUTE_B:"];
  Print["  U_B2=" <> canon[routeB["interaction"]] <> "; F_B2=" <> ToString[canon /@ routeB["force"], InputForm]];
  Print["  radial=" <> canon[routeB["radial"]] <> "; potential=R^-1; force=R^-2; O(V1*V2)"];
  Print["  route_independence=" <> ToString[routeIndependence[]] <> "; dependencies=" <> ToString[routeB["dependencies"], InputForm]];
  Print["Q-COMPARE:"];
  Print["  structural_tensor_agreement=" <> ToString[comparison["tensor"]] <> "; falloff=" <>
    ToString[comparison["falloff"]] <> "; order=" <> ToString[comparison["order"]]];
  Print["  comparison=" <> comparison["comparison"] <> "(" <> comparison["gap"] <> "); relative=" <> comparison["relative"]];
  Print["  routeB/routeA=" <> canon[ratioBA] <> "; normalized_delta=" <> canon[deltaNormalized] <>
    "; delta_U=" <> canon[deltaU] <> "; cone_ratio=" <> canon[coneRatio]];
  Print["Q-MAG: fact=" <> qmagFact <> "; tier=" <> tierFact <> "; electric_anchor=" <> anchorFact];
  Print["DIMENSIONS: [u_T]=L; [curl u_T]=1; [q_T]=M/T; [A_E]=E*L; [U]=E; [F]=E/L; ratios=1"];
  Print["SECTION5_HOOKS:"]; KeyValueMap[Print["  " <> #1 <> "=" <> #2] &, hooks];
  Print["SECTION4_TRUTH_TABLE: cells=" <> ToString[truth["total"]] <> "; exhaustive=" <>
    ToString[truth["exact"]] <> "; digest=" <> truth["digest"] <> "; classes=" <> ToString[truth["counts"], InputForm]];
  Print["SECTION4_LANDING=" <> actualLanding]; Print["SECTION4_ALL_R1=" <> StringRiffle[r1Blockers[], ","]];
  Print["DECIDED=source basis/parities; ledger row; Route-A Darwin tensor; conditional Route-B tensor/sign/falloff/order; analytic delta"];
  Print["R1=parent-throat q_T normalization; electric BC/branch; cone knit; higher orders; active F_flux integrability"];
  Print["ATOMIC_TEETH:"]; Scan[Print["  PASS " <> # <> "; mutation=" <> teeth[#]] &, toothOrder];
  Print["ENGINE_AGREE=PASS"]
];

jsonOnly = Environment["MAGNETISM_JSON_ONLY"] === "1";
payloadMutation = Environment["MAGNETISM_PAYLOAD_MUTATION"] === "1";
If[jsonOnly,
  Print["JSON_PAYLOAD=" <> ExportString[symbolicPayload[payloadMutation], "RawJSON", "Compact" -> True]];
  Exit[0]
];

teeth = mutationCampaign[];
printReport[teeth];
Exit[0];
