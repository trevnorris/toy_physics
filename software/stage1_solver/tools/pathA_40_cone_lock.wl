(* pathA_40 cone-lock gate, Mathematica engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_40_cone_lock.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
reportsDir = FileNameJoin[{stage1Root, "reports"}];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
sympyJson = FileNameJoin[{scratchDir, "pathA_40_cone_lock_sympy.json"}];
jsonOut = FileNameJoin[{scratchDir, "pathA_40_cone_lock_mathematica.json"}];

Scan[If[! FileExistsQ[#], fail["missing input: " <> #]] &, {
  FileNameJoin[{reportsDir, "pathA_35_G0_freeze.md"}],
  FileNameJoin[{reportsDir, "pathA_36_c5_phase_potential_results.yaml"}],
  FileNameJoin[{reportsDir, "pathA_38_throat_body_electric_localization.md"}],
  FileNameJoin[{reportsDir, "pathA_38_results.yaml"}],
  FileNameJoin[{reportsDir, "pathA_39_scalar_admixture_screen.md"}],
  FileNameJoin[{reportsDir, "pathA_39_stage4_field_classification_results.yaml"}],
  sympyJson
}];

p35Text = ReadString[FileNameJoin[{reportsDir, "pathA_35_G0_freeze.md"}]];
p36Text = ReadString[FileNameJoin[{reportsDir, "pathA_36_c5_phase_potential_results.yaml"}]];
p38Text = ReadString[FileNameJoin[{reportsDir, "pathA_38_throat_body_electric_localization.md"}]];
p38YamlText = ReadString[FileNameJoin[{reportsDir, "pathA_38_results.yaml"}]];
p39ScalarText = ReadString[FileNameJoin[{reportsDir, "pathA_39_scalar_admixture_screen.md"}]];
p39Stage4Text = ReadString[FileNameJoin[{reportsDir, "pathA_39_stage4_field_classification_results.yaml"}]];

assertContains[text_, token_, label_] := If[! StringContainsQ[text, token], fail["source token missing: " <> label]];
assertContains[p35Text, "postulated-ingredient", "pathA_35 postulated tags"];
assertContains[p35Text, "`rho_br` | `M L^-3` | postulated-ingredient", "rho_br postulated"];
assertContains[p35Text, "`mu_R` | `M L^-1 T^-2` | postulated-ingredient", "mu_R postulated"];
assertContains[p35Text, "The in-plane displacement `u^a", "in-plane shear DOF"];
assertContains[p38Text, "full nonlinear throat source compactness", "deferred throat program"];
assertContains[p38YamlText, "sech(w/ell)^2/ell", "h-Goldstone throat profile"];
assertContains[p36Text, "c_gamma_squared: mu_R/rho_br", "c_gamma import"];
assertContains[p36Text, "classification: SECOND_CLASS_PAIR", "real branch class"];
assertContains[p39ScalarText, "C_hu scalar mixing", "C_hu scan parameter"];
assertContains[p39ScalarText, "report does not give a closed profile", "missing closed f_u"];
assertContains[p39Stage4Text, "primary: FIELD_SCALAR_VECTOR_DEPARTURE", "field content"];
assertContains[p39Stage4Text, "status: ENGINE_AGREE", "stage-4 engine agreement"];
assertContains[p39Stage4Text, "real: 4", "real DOF count"];
assertContains[p39Stage4Text, "maxwell: 2", "Maxwell counterfactual DOF count"];

rho = Unique["rho"]; rhoBr = Unique["rhoBr"]; rhoB0 = Unique["rhoB0"];
chiC = Unique["chiC"]; muR = Unique["muR"]; kk = Unique["bigK"];
mm = Unique["massm"]; Mh = Unique["Mh"]; cE = Unique["cE"];
Chu = Unique["Chu"]; Kh = Unique["Kh"]; Beff = Unique["Beff"]; sigma = Unique["sigma"];
tauForced = Unique["tauForced"]; tauA = Unique["tauA"]; tauB = Unique["tauB"];
etaOver = Unique["etaOver"]; qhSq = Unique["qhSq"];
kSym = Unique["k"]; omegaSym = Unique["omega"];

baseVars = {rho, rhoBr, rhoB0, chiC, muR, kk, mm, Mh, cE, Chu, Kh, Beff, sigma};
basePositive = {rho, rhoBr, chiC, muR, kk, mm, Mh, cE, Kh, Beff, sigma};
baseEqs = {
  Beff*chiC - rhoB0^2,
  Kh - Mh*cE^2,
  Beff*Kh - Chu^2 - sigma
};
lockA = mm*muR - 5*kk*rho^4*rhoBr;
lockB = cE^2*rhoBr - muR;

has[kinds_, kind_] := MemberQ[kinds, kind];

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

caseKinds[case_] := Which[
  case === "production", baselineKinds[],
  case === "well_posed", {
    "candidate_bridge", "nonlinear_action_with_gnls_and_u", "in_plane_shear_profile_fu",
    "common_normalization_rho_mu", "reduction_mu_over_rho_equals_cs", "no_residual_geometric_factor",
    "route_b_distinct_dof", "over_import_guard", "freedom_free_parameter:C_hu", "freedom_free_parameter:rho_br"
  },
  case === "absent", {
    "postulated_mu_rho", "route_b_distinct_dof", "over_import_guard",
    "freedom_free_parameter:C_hu", "freedom_free_parameter:rho_br"
  },
  case === "partial_inventory", {
    "candidate_bridge", "nonlinear_action_with_gnls_and_u", "in_plane_shear_profile_fu", "postulated_mu_rho",
    "route_b_distinct_dof", "over_import_guard", "freedom_free_parameter:C_hu", "freedom_free_parameter:rho_br"
  },
  case === "forced_lock", Append[baselineKinds[], "forced_lock_fake_relation"],
  case === "a_only_partial", Append[baselineKinds[], "force_A_fake_relation"],
  case === "b_only_partial", Append[baselineKinds[], "force_B_fake_relation"],
  case === "over_constrained", Append[baselineKinds[], "over_stability_relation"],
  case === "freedom_tie", Append[baselineKinds[], "freedom_tie_relation"],
  True, fail["unknown case " <> ToString[case]]
];

routeA[kinds_] := Module[{r1, r2, r3, r4, r5, missing, grade, candidate},
  r1 = If[has[kinds, "nonlinear_action_with_gnls_and_u"], "present", "absent"];
  r2 = If[has[kinds, "in_plane_shear_profile_fu"], "present", "absent"];
  r3 = If[has[kinds, "common_normalization_rho_mu"], "present", "absent"];
  r4 = If[has[kinds, "reduction_mu_over_rho_equals_cs"], "present", "absent"];
  r5 = If[r4 === "present", If[has[kinds, "no_residual_geometric_factor"], "present", "absent"], "not_applicable"];
  missing = Cases[
    MapThread[If[#2 =!= "present", #1, Nothing] &, {{"R1", "R2", "R3", "R4", "R5"}, {r1, r2, r3, r4, r5}}],
    _String
  ];
  candidate = has[kinds, "candidate_bridge"];
  grade = Which[
    Length[missing] === 0, "ROUTE_A_WELL_POSED",
    candidate, "ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT",
    True, "ROUTE_A_ABSENT"
  ];
  <|"grade" -> grade, "missing" -> missing|>
];

routeB[kinds_] := If[has[kinds, "route_b_distinct_dof"] && ! has[kinds, "thin_plate_over_import_relation"],
  "ROUTE_B_CLOSED_CHECKED_NEGATIVE",
  "ROUTE_B_GUARD_FAIL"
];

freedomStatus[kinds_] := If[has[kinds, "freedom_tie_relation"], "FREEDOM_TIED", "FREEDOM_UNCONSTRAINED"];

sourceRelations[kinds_] := Module[{eqs = {}, vars = {}, pos = {}, rels = {}},
  If[has[kinds, "forced_lock_fake_relation"],
    eqs = Join[eqs, {muR - rhoBr*tauForced, cE^2 - tauForced, mm*tauForced - 5*kk*rho^4}];
    vars = Append[vars, tauForced]; pos = Append[pos, tauForced]; rels = Append[rels, "forced_lock_fake_relation"];
  ];
  If[has[kinds, "force_A_fake_relation"],
    eqs = Join[eqs, {muR - rhoBr*tauA, mm*tauA - 5*kk*rho^4}];
    vars = Append[vars, tauA]; pos = Append[pos, tauA]; rels = Append[rels, "force_A_fake_relation"];
  ];
  If[has[kinds, "force_B_fake_relation"],
    eqs = Join[eqs, {muR - rhoBr*tauB, cE^2 - tauB}];
    vars = Append[vars, tauB]; pos = Append[pos, tauB]; rels = Append[rels, "force_B_fake_relation"];
  ];
  If[has[kinds, "over_stability_relation"],
    eqs = Append[eqs, Chu^2 - Beff*Kh - etaOver];
    vars = Append[vars, etaOver]; pos = Append[pos, etaOver]; rels = Append[rels, "over_stability_relation"];
  ];
  If[has[kinds, "freedom_tie_relation"],
    eqs = Join[eqs, {Chu^2 - qhSq, qhSq*rhoBr - 2*Beff*Mh*muR}];
    vars = Append[vars, qhSq]; pos = Append[pos, qhSq]; rels = Append[rels, "freedom_tie_relation"];
  ];
  <|"eqs" -> eqs, "vars" -> vars, "positive" -> pos, "relations" -> rels|>
];

varsFor[kinds_] := Join[baseVars, sourceRelations[kinds]["vars"]];
positiveFor[kinds_] := Join[basePositive, sourceRelations[kinds]["positive"]];
eqsFor[kinds_, includeLocks_] := If[TrueQ[includeLocks],
  Join[baseEqs, sourceRelations[kinds]["eqs"], {lockA, lockB}],
  Join[baseEqs, sourceRelations[kinds]["eqs"]]
];
domainFor[kinds_] := (And @@ Thread[positiveFor[kinds] > 0]) && rhoB0 != 0 && Element[Chu, Reals];
formulaFor[kinds_, includeLocks_] := (And @@ Thread[eqsFor[kinds, includeLocks] == 0]) && domainFor[kinds];

satStatus[kinds_, includeLocks_] := Module[{vars, res},
  vars = varsFor[kinds];
  res = TimeConstrained[Resolve[Exists[Evaluate[vars], formulaFor[kinds, includeLocks]], Reals], 60, "UNKNOWN"];
  Which[res === True, "SAT", res === False, "UNSAT", True, "SAT_UNCERTIFIED"]
];

entailStatus[kinds_, lock_] := Module[{vars, baseFormula, lockExpr, universal, inst},
  vars = varsFor[kinds];
  baseFormula = formulaFor[kinds, False];
  lockExpr = If[lock === "A", lockA, lockB];
  universal = TimeConstrained[Resolve[ForAll[Evaluate[vars], Implies[baseFormula, lockExpr == 0]], Reals], 60, "UNKNOWN"];
  If[universal === True, Return["ENTAILED"]];
  inst = TimeConstrained[FindInstance[baseFormula && lockExpr != 0, Evaluate[vars], Reals, 1], 60, $Failed];
  If[ListQ[inst] && Length[inst] > 0, "WITNESSED", "INCONCLUSIVE"]
];

regionDim[kinds_, includeLocks_] := Module[{vars, reg, d},
  vars = varsFor[kinds];
  reg = ImplicitRegion[formulaFor[kinds, includeLocks], Evaluate[vars]];
  d = TimeConstrained[RegionDimension[reg], 60, "DIMENSION_UNCERTIFIED"];
  d
];

dimensionRecord[kinds_, sector_, locked_] := Module[{before, after},
  If[sector =!= "SAT" || locked =!= "SAT", Return[<|"status" -> "NOT_RUN"|>]];
  before = regionDim[kinds, False];
  after = regionDim[kinds, True];
  If[IntegerQ[before] && IntegerQ[after],
    <|"status" -> "CERTIFIED", "dim_before" -> before, "dim_after" -> after, "delta_r" -> before - after|>,
    <|"status" -> "DIMENSION_UNCERTIFIED", "dim_before" -> before, "dim_after" -> after, "delta_r" -> Null|>
  ]
];

decideCase[raGrade_, sector_, locked_, provA_, provB_, dimStatus_, delta_] := Which[
  sector === "UNSAT", "NO_GO(sector-ledger)",
  raGrade === "ROUTE_A_WELL_POSED", "HALT_ROUTE_A_WELL_POSED",
  locked === "UNSAT", "NO_GO(cone-lock)",
  locked =!= "SAT" || dimStatus =!= "CERTIFIED" || MemberQ[{provA, provB}, "INCONCLUSIVE"], "CONE_LOCK_PROVENANCE_INCONCLUSIVE",
  provA === "ENTAILED" && provB === "ENTAILED" && delta === 0, "CONE_LOCK_DERIVED",
  provA === "ENTAILED" && provB === "WITNESSED" && delta === 1, "CONE_LOCK_PARTIAL(derived=A, calibrated=B)",
  provB === "ENTAILED" && provA === "WITNESSED" && delta === 1, "CONE_LOCK_PARTIAL(derived=B, calibrated=A)",
  provA === "WITNESSED" && provB === "WITNESSED" && delta === 2, "CONE_LOCK_CALIBRATED",
  provA === "WITNESSED" && provB === "WITNESSED" && delta === 1, "CONE_LOCK_SHARED_CALIBRATION(delta_r=1, derived=none)",
  True, "CONE_LOCK_PROVENANCE_INCONCLUSIVE"
];

compactCase[case_] := Module[{kinds, ra, rb, freedom, sector, locked, provA = Null, provB = Null, dim, verdict},
  kinds = caseKinds[case];
  ra = routeA[kinds]; rb = routeB[kinds]; freedom = freedomStatus[kinds];
  If[ra["grade"] === "ROUTE_A_WELL_POSED",
    sector = "NOT_RUN"; locked = "NOT_RUN"; dim = <|"status" -> "NOT_RUN"|>; verdict = "HALT_ROUTE_A_WELL_POSED",
    sector = satStatus[kinds, False];
    locked = satStatus[kinds, True];
    If[sector === "SAT" && locked === "SAT",
      provA = entailStatus[kinds, "A"];
      provB = entailStatus[kinds, "B"];
    ];
    dim = dimensionRecord[kinds, sector, locked];
    verdict = decideCase[ra["grade"], sector, locked, provA, provB, dim["status"], Lookup[dim, "delta_r", Null]];
  ];
  <|
    "verdict" -> verdict,
    "route_a_grade" -> ra["grade"],
    "route_a_missing" -> ra["missing"],
    "route_b_status" -> rb,
    "freedom_status" -> freedom,
    "sector_sat" -> sector,
    "lock_sat" -> locked,
    "prov_A" -> provA,
    "prov_B" -> provB,
    "dimension_status" -> dim["status"],
    "delta_r" -> Lookup[dim, "delta_r", Null],
    "dim_before" -> Lookup[dim, "dim_before", Null],
    "dim_after" -> Lookup[dim, "dim_after", Null],
    "riders" -> {}
  |>
];

detM = (rhoBr*omegaSym^2 - Beff*kSym^2)*(Mh*omegaSym^2 - Kh*kSym^2) - Chu^2*kSym^4;
detOnCone = FullSimplify[detM /. {omegaSym^2 -> (muR/rhoBr)*kSym^2, Kh -> Mh*muR/rhoBr}];
If[! TrueQ[FullSimplify[detOnCone == -Chu^2*kSym^4]], fail["det-on-cone residual failed"]];

trace = rhoBr*Kh + Mh*Beff;
discr = (rhoBr*Kh - Mh*Beff)^2 + 4*rhoBr*Mh*Chu^2;
vMinus = FullSimplify[(trace - Sqrt[discr])/(2*rhoBr*Mh)];
vPlus = FullSimplify[(trace + Sqrt[discr])/(2*rhoBr*Mh)];
deltaMinusUnderB = FullSimplify[(vMinus - muR/rhoBr) /. muR -> Kh*rhoBr/Mh];
deltaPlusUnderB = FullSimplify[(vPlus - muR/rhoBr) /. muR -> Kh*rhoBr/Mh];
If[! TrueQ[FullSimplify[deltaMinusUnderB ==
    (Beff*Mh - Kh*rhoBr - Sqrt[Beff^2*Mh^2 - 2*Beff*Kh*Mh*rhoBr + 4*Chu^2*Mh*rhoBr + Kh^2*rhoBr^2])/(2*Mh*rhoBr)]],
  fail["Delta minus derivation failed"]
];
If[! TrueQ[FullSimplify[deltaPlusUnderB ==
    (Beff*Mh - Kh*rhoBr + Sqrt[Beff^2*Mh^2 - 2*Beff*Kh*Mh*rhoBr + 4*Chu^2*Mh*rhoBr + Kh^2*rhoBr^2])/(2*Mh*rhoBr)]],
  fail["Delta plus derivation failed"]
];

caseNames = {"production", "well_posed", "absent", "partial_inventory", "forced_lock", "a_only_partial",
  "b_only_partial", "over_constrained", "freedom_tie"};
casePayload = AssociationThread[caseNames, compactCase /@ caseNames];

payload = <|
  "production_and_controls" -> casePayload,
  "field_content" -> "FIELD_SCALAR_VECTOR_DEPARTURE",
  "det_on_light_cone_under_B" -> "-C_hu**2*k**4",
  "Delta_pole_minus_under_B" -> "(B_eff*M_h - K_h*rho_br - sqrt(B_eff**2*M_h**2 - 2*B_eff*K_h*M_h*rho_br + 4*C_hu**2*M_h*rho_br + K_h**2*rho_br**2))/(2*M_h*rho_br)",
  "Delta_pole_plus_under_B" -> "(B_eff*M_h - K_h*rho_br + sqrt(B_eff**2*M_h**2 - 2*B_eff*K_h*M_h*rho_br + 4*C_hu**2*M_h*rho_br + K_h**2*rho_br**2))/(2*M_h*rho_br)"
|>;

digest = IntegerString[Hash[ExportString[payload, "RawJSON"], "SHA256"], 16, 64];

sympyResults = Import[sympyJson, "RawJSON"];
If[! AssociationQ[sympyResults], fail["could not read SymPy JSON"]];
canon[x_Association] := SortBy[({#[[1]], canon[#[[2]]]} & /@ Normal[x]), First];
canon[x_List] := canon /@ x;
canon[x_] := x;
If[canon[sympyResults["comparison_payload"]] =!= canon[payload],
  Export[FileNameJoin[{scratchDir, "pathA_40_cone_lock_mathematica_debug.json"}],
    <|"sympy" -> sympyResults["comparison_payload"], "mathematica" -> payload|>, "RawJSON"];
  fail["Mathematica independently computed payload does not match SymPy payload"]
];

out = <|
  "schema" -> "pathA_40_cone_lock_mathematica/v1",
  "engine" -> "mathematica",
  "status" -> "OK",
  "method" -> <|
    "sat" -> "Resolve[Exists[..., semialgebraic formula], Reals]",
    "entailment" -> "Resolve[ForAll[..., E implies L_i], Reals] or FindInstance witness",
    "dimension" -> "CAD-backed RegionDimension[ImplicitRegion[...]]",
    "field_overlay" -> "symbolic determinant/root derivation"
  |>,
  "comparison_payload" -> payload,
  "comparison_digest" -> digest
|>;

Export[jsonOut, out, "RawJSON"];
Print["OK pathA_40_cone_lock_mathematica -> " <> jsonOut];
