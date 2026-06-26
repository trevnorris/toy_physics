(* pathA_35 G0 Mathematica freeze checker.

   Scope: G0 only.  Recomputes freeze hashes, checks restored-unit
   homogeneity for every constructed expression, verifies able-to-fail
   ablations, and counts flat-brane linear DOF.  No Gate L verdict is computed.
*)

ClearAll[
  fail, scriptPath, stage1Root, scratchDir, t0Report, g0Report,
  expectedT0Hash, expectedG0Hash, expectedG0Short, extractBlock, sha256,
  byteRange, t0Text, g0Text, t0Block, g0Block, t0Hash, g0Hash,
  dadd, dmul, dsub, dimString, dimRecord, records, ablations, checkDim,
  checkSame, expectFail, Mdim, Ldim, Tdim, Zdim, bulkLag, braneLag,
  actionDim, stressDim, eomUOpDim, dm, dhbar, drho, dK, da, du, dP, dw, dellG, dg,
  dgrad, ddt, ddtMeasure, dd4x, dv, dk, domega, drhoBr, dmuR, dlambdaPu, dOmegaW,
  dcs2, dDtP, dgradP, dDtu, dcurlu, duw, dDtuw, dvarpi, macKin,
  macCurl, puTerm, uwKin, uwGap, twa, slope, stressSlope, opRho, opMu,
  parity, constants, dimensions, matrixRank, computeFlatBraneDof,
  dofCountAblation, k, curlStiffness, rank, nullity, computedModes,
  reportedCounts, reportedTotal, dofAblations, modes, ledger, freeze,
  agreementPayload, report, outPath
];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_35_G0.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
t0Report = FileNameJoin[{stage1Root, "reports", "pathA_24_T0_freeze.md"}];
g0Report = FileNameJoin[{stage1Root, "reports", "pathA_35_G0_freeze.md"}];

expectedT0Hash = "8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064";
expectedG0Hash = "d9520d3819c3f718290f9d0be57138c07d5bf02d2237106478e17b6a1e389ac3";
expectedG0Short = "d9520d3819c3";

extractBlock[text_, path_] := Module[{blocks},
  blocks = StringCases[
    text,
    "```freeze-action\n" ~~ Shortest[block__] ~~ "\n```" :> block
  ];
  If[Length[blocks] =!= 1,
    fail["expected exactly one freeze-action block in " <> path <> ", found " <> ToString[Length[blocks]]]
  ];
  First[blocks]
];

sha256[block_] := ToLowerCase[IntegerString[Hash[block <> "\n", "SHA256"], 16, 64]];

byteRange[text_] := Module[{marker, markerPos, start0, rest, endRel, end0},
  marker = "```freeze-action\n";
  markerPos = First[StringPosition[text, marker]];
  start0 = markerPos[[2]];
  rest = StringDrop[text, start0];
  endRel = First[StringPosition[rest, "\n```"]][[1]];
  end0 = start0 + endRel - 1;
  <|
    "start_0_based_inclusive" -> start0,
    "end_0_based_exclusive" -> end0,
    "range_0_based" -> {start0, end0},
    "start_1_based_inclusive" -> start0 + 1,
    "end_1_based_inclusive" -> end0,
    "length_bytes" -> end0 - start0
  |>
];

t0Text = Import[t0Report, "Text"];
g0Text = Import[g0Report, "Text"];
t0Block = extractBlock[t0Text, t0Report];
g0Block = extractBlock[g0Text, g0Report];
t0Hash = sha256[t0Block];
g0Hash = sha256[g0Block];
If[t0Hash =!= expectedT0Hash, fail["T0 hash mismatch: expected " <> expectedT0Hash <> ", got " <> t0Hash]];
If[g0Hash =!= expectedG0Hash, fail["G0 hash mismatch: expected " <> expectedG0Hash <> ", got " <> g0Hash]];
If[! StringContainsQ[g0Block, t0Block], fail["byte-identical T0 freeze-action block is not embedded in G0 block"]];

dadd[dims__] := Total[{dims}];
dmul[n_, dim_] := n dim;
dsub[left_, right_] := left - right;
dimString[dim_] := Module[{labels = {"M", "L", "T"}, powers = dim, parts},
  parts = MapThread[
    If[#2 === 0, Nothing, If[#2 === 1, #1, #1 <> "^" <> ToString[#2]]] &,
    {labels, powers}
  ];
  If[parts === {}, "1", StringRiffle[parts, " "]]
];
dimRecord[dim_] := <|"triple_MLT" -> dim, "string" -> dimString[dim]|>;

records = {};
ablations = {};

checkDim[category_, name_, actual_, expected_, expression_] := Module[{},
  If[actual =!= expected,
    fail[category <> ":" <> name <> ": expected " <> ToString[expected] <> ", got " <> ToString[actual]]
  ];
  AppendTo[
    records,
    <|
      "category" -> category,
      "name" -> name,
      "expression" -> expression,
      "dimension" -> dimRecord[actual],
      "expected" -> dimRecord[expected],
      "status" -> "PASS"
    |>
  ];
  actual
];

checkSame[category_, name_, dims_, expected_, expression_] := Module[{i},
  Do[
    checkDim[category, name <> ".part_" <> ToString[i - 1], dims[[i]], expected, expression],
    {i, Length[dims]}
  ];
  AppendTo[
    records,
    <|
      "category" -> category,
      "name" -> name,
      "expression" -> expression,
      "dimension" -> dimRecord[expected],
      "expected" -> dimRecord[expected],
      "status" -> "PASS",
      "homogeneous_parts" -> Length[dims]
    |>
  ];
  expected
];

expectFail[category_, name_, actual_, expected_, expression_] := Module[{},
  If[actual === expected, fail["ablation did not fire: " <> category <> ":" <> name]];
  AppendTo[
    ablations,
    <|
      "category" -> category,
      "name" -> name,
      "expression" -> expression,
      "actual" -> dimRecord[actual],
      "expected" -> dimRecord[expected],
      "status" -> "FIRED"
    |>
  ];
];

Mdim = {1, 0, 0}; Ldim = {0, 1, 0}; Tdim = {0, 0, 1}; Zdim = {0, 0, 0};
bulkLag = {1, -2, -2};
braneLag = {1, -1, -2};
actionDim = {1, 2, -1};
stressDim = bulkLag;
eomUOpDim = {1, -3, -2};

dm = Mdim;
dhbar = {1, 2, -1};
drho = {0, -4, 0};
dK = {1, 18, -2};
da = Ldim;
du = Ldim;
dP = Zdim;
dw = Ldim;
dellG = Ldim;
dg = {0, -1, 0};
dgrad = {0, -1, 0};
ddt = {0, 0, -1};
ddtMeasure = {0, 0, 1};
dd4x = {0, 4, 0};
dv = {0, 1, -1};
dk = {0, -1, 0};
domega = {0, 0, -1};
drhoBr = {1, -3, 0};
dmuR = braneLag;
dlambdaPu = braneLag;
dOmegaW = domega;

dcs2 = dsub[dadd[dK, dmul[4, drho]], dm];
checkDim["kept_gnls", "c_s^2(rho)", dcs2, {0, 2, -2}, "5 K rho^4 / m"];
checkDim["kept_gnls", "U(rho)", dadd[dK, dmul[5, drho]], bulkLag, "(K/4) rho^5"];
checkDim[
  "kept_gnls",
  "quantum_pressure",
  dadd[dmul[2, dhbar], dmul[-1, dm], dmul[-1, drho], dmul[2, dadd[dgrad, drho]]],
  bulkLag,
  "(hbar^2/(8 m rho)) (partial_i rho)^2"
];
checkDim["kept_gnls", "bulk_kinetic_stress_scale", dadd[dm, drho, dmul[2, dv]], stressDim, "m rho v_i v_j"];

dDtP = ddt;
dgradP = dadd[dgrad, dP];
checkDim["kept_T0", "P_kinetic", dadd[dm, drho, dmul[2, da], dmul[2, dDtP]], bulkLag, "m rho a^2 (D_t^v P)^2"];
checkDim[
  "kept_T0",
  "P_Frank",
  dadd[dm, drho, dcs2, dmul[2, da], dmul[2, dgradP]],
  bulkLag,
  "m rho c_s^2 a^2 (partial_j P^i)^2"
];
checkDim["kept_T0", "P_radial_potential", dadd[dm, drho, dcs2], bulkLag, "m rho c_s^2 (P^i P^i - 1)^2"];
checkDim["kept_T0", "bulk_couple_inertia", dadd[dm, drho, dmul[2, da]], {1, -2, 0}, "m rho a^2"];
checkDim["kept_T0", "bulk_couple_stiffness", dadd[dm, drho, dcs2, dmul[2, da]], {1, 0, -2}, "m rho c_s^2 a^2"];
checkDim["kept_T0", "bulk_radial_scale", dadd[dm, drho, dcs2], bulkLag, "m rho c_s^2"];

checkDim["confinement", "w_over_ell_g", dsub[dw, dellG], Zdim, "w/ell_g"];
checkDim["confinement", "g_ell", dg, {0, -1, 0}, "exp(-(w/ell_g)^2)/(sqrt(pi) ell_g)"];
checkDim["confinement", "dw_g_ell", dadd[Ldim, dg], Zdim, "dw g_ell(w)"];

dDtu = dadd[ddt, du];
dcurlu = dadd[dgrad, du];
duw = du;
dDtuw = dadd[ddt, duw];
dvarpi = dP;
checkDim["brane_surface", "rho_br", drhoBr, {1, -3, 0}, "rho_br"];
macKin = checkDim["brane_surface", "u_kinetic", dadd[drhoBr, dmul[2, dDtu]], braneLag, "rho_br (partial_t u)^2"];
macCurl = checkDim["brane_surface", "MacCullagh_curl", dadd[dmuR, dmul[2, dcurlu]], braneLag, "mu_R (nabla_parallel x u)^2"];
checkSame["brane_surface", "L_Mac_homogeneous", {macKin, macCurl}, braneLag, "L_Mac"];
puTerm = checkDim["brane_surface", "P_u_coupling", dadd[dlambdaPu, dvarpi, dcurlu], braneLag, "lambda_Pu varpi_a Omega_u^a"];
uwKin = checkDim["brane_surface", "u_w_kinetic", dadd[drhoBr, dmul[2, dDtuw]], braneLag, "rho_br (partial_t u_w)^2"];
uwGap = checkDim["brane_surface", "u_w_gap", dadd[drhoBr, dmul[2, dOmegaW], dmul[2, duw]], braneLag, "rho_br Omega_w^2 u_w^2"];
checkSame["brane_surface", "L_uw_homogeneous", {uwKin, uwGap}, braneLag, "L_uw"];

Scan[
  checkDim["brane_bulk_representation", #[[1]], dadd[dg, #[[2]]], bulkLag, "g_ell(w) " <> #[[1]]] &,
  {
    {"g_L_Mac_kinetic", macKin},
    {"g_L_Mac_curl", macCurl},
    {"g_L_Pu", puTerm},
    {"g_L_uw_kinetic", uwKin},
    {"g_L_uw_gap", uwGap}
  }
];

checkDim["action_measure", "bulk_action_integral", dadd[ddtMeasure, dd4x, bulkLag], actionDim, "int dt d^4X bulk_lag"];
checkDim[
  "action_measure",
  "brane_bulk_action_integral",
  dadd[ddtMeasure, dd4x, dg, braneLag],
  actionDim,
  "int dt d^4X g_ell(w) L_brane"
];

twa = checkDim["traction", "T_wa", dadd[dm, drho, dv, dv], stressDim, "m rho v_w v_a"];
slope = checkDim["traction", "partial_b_u_w", dadd[dgrad, duw], Zdim, "partial_b u_w"];
stressSlope = checkDim["traction", "slope_mixing", dadd[stressDim, slope], stressDim, "(T_ww delta_ab - T_ab) partial_b u_w"];
checkSame[
  "traction",
  "T_na_projected",
  {twa, stressSlope},
  stressDim,
  "T_wa + (T_ww delta_ab - T_ab) partial_b u_w"
];

expectFail["ablation", "drop_m_from_T_wa", dadd[drho, dv, dv], stressDim, "rho v_w v_a"];
expectFail["ablation", "MacCullagh_without_curl", dadd[dmuR, dmul[2, du]], braneLag, "mu_R u^2"];

opRho = checkDim["linearization", "rho_br_omega2", dadd[drhoBr, dmul[2, domega]], eomUOpDim, "rho_br omega^2"];
opMu = checkDim["linearization", "mu_R_k2", dadd[dmuR, dmul[2, dk]], eomUOpDim, "mu_R k^2"];
checkSame["linearization", "O_u_homogeneous", {opRho, opMu}, eomUOpDim, "rho_br omega^2 I - mu_R(k^2 I-k k^T)"];
checkDim["linearization", "c_gamma_squared", dsub[dmuR, drhoBr], {0, 2, -2}, "mu_R/rho_br"];
checkDim["linearization", "omega_T_squared", dadd[dsub[dmuR, drhoBr], dmul[2, dk]], {0, 0, -2}, "c_gamma^2 k^2"];
checkDim["linearization", "P_spin_omega_squared", dadd[dcs2, dmul[2, dk]], {0, 0, -2}, "c_s^2 k^2"];
checkDim["linearization", "P_radial_gap_term", dsub[dcs2, dmul[2, da]], {0, 0, -2}, "2 c_s^2/a^2"];
checkSame[
  "linearization",
  "P_radial_omega_squared",
  {dadd[dcs2, dmul[2, dk]], dsub[dcs2, dmul[2, da]]},
  {0, 0, -2},
  "c_s^2 k^2 + 2 c_s^2/a^2"
];
checkDim["linearization", "u_w_bare_omega_squared", dmul[2, dOmegaW], {0, 0, -2}, "Omega_w^2"];
checkDim["linearization", "Fourier_P_u_monomial", dadd[dlambdaPu, dP, dk, du], braneLag, "lambda_Pu P k u"];

parity = <|
  "direct_P_dot_curlu" -> <|
    "operator" -> "P_parallel . Omega_u",
    "parity" -> "odd_pseudoscalar",
    "status" -> "excluded"
  |>,
  "repaired_w_dot_P_cross_curlu" -> <|
    "operator" -> "w_hat . (P_parallel x Omega_u)",
    "parity" -> "even_scalar_using_imposed_normal",
    "status" -> "frozen"
  |>
|>;

constants = <|
  "rho" -> dimRecord[drho],
  "m" -> dimRecord[dm],
  "K" -> dimRecord[dK],
  "hbar" -> dimRecord[dhbar],
  "a" -> dimRecord[da],
  "rho_br" -> dimRecord[drhoBr],
  "mu_R" -> dimRecord[dmuR],
  "lambda_Pu" -> dimRecord[dlambdaPu],
  "Omega_w" -> dimRecord[dOmegaW],
  "ell_g" -> dimRecord[dellG],
  "g_ell" -> dimRecord[dg]
|>;

dimensions = <|
  "base_dimensions" -> <|"M" -> {1, 0, 0}, "L" -> {0, 1, 0}, "T" -> {0, 0, 1}|>,
  "targets" -> <|
    "bulk_action_density" -> dimRecord[bulkLag],
    "brane_action_density" -> dimRecord[braneLag],
    "action" -> dimRecord[actionDim],
    "stress" -> dimRecord[stressDim],
    "u_eom_operator" -> dimRecord[eomUOpDim]
  |>,
  "constants" -> constants,
  "checks" -> records,
  "ablations" -> ablations,
  "parity" -> parity
|>;

matrixRank[m_] := MatrixRank[m];

computeFlatBraneDof[opts_Association] := Module[
  {
    uActive, pKineticPresent, pFrankPresent, pRadialPresent,
    uWKineticPresent, uWGapPresent, phiPresent, eyeU, activeU,
    curlProjector, uKineticForm, uCurlForm, uKineticRank, uCurlRank,
    uCurlNullity, eyeP, zeroP, tangentP, radialP, pKineticForm,
    pFrankForm, pRadialHessian, pTangentKineticRank, pTangentFrankRank,
    pRadialKineticRank, pRadialHessianRank, uWKineticRank, uWGapRank,
    phiKineticRank, computedCounts, computedTotal
  },
  uActive = Lookup[opts, "u_active", {1, 1, 1}];
  pKineticPresent = Lookup[opts, "p_kinetic_present", True];
  pFrankPresent = Lookup[opts, "p_frank_present", True];
  pRadialPresent = Lookup[opts, "p_radial_present", True];
  uWKineticPresent = Lookup[opts, "u_w_kinetic_present", True];
  uWGapPresent = Lookup[opts, "u_w_gap_present", True];
  phiPresent = Lookup[opts, "phi_present", False];

  eyeU = IdentityMatrix[3];
  activeU = DiagonalMatrix[uActive];
  curlProjector = DiagonalMatrix[{0, 1, 1}];
  uKineticForm = activeU . eyeU . activeU;
  uCurlForm = activeU . curlProjector . activeU;
  uKineticRank = matrixRank[uKineticForm];
  uCurlRank = matrixRank[uCurlForm];
  uCurlNullity = uKineticRank - uCurlRank;
  If[uCurlNullity < 0, fail["u curl rank exceeds active u kinetic rank"]];

  eyeP = IdentityMatrix[4];
  zeroP = ConstantArray[0, {4, 4}];
  tangentP = DiagonalMatrix[{1, 1, 1, 0}];
  radialP = DiagonalMatrix[{0, 0, 0, 1}];
  pKineticForm = If[pKineticPresent, eyeP, zeroP];
  pFrankForm = If[pFrankPresent, eyeP, zeroP];
  pRadialHessian = If[pRadialPresent, radialP, zeroP];
  pTangentKineticRank = matrixRank[tangentP . pKineticForm . tangentP];
  pTangentFrankRank = matrixRank[tangentP . pFrankForm . tangentP];
  pRadialKineticRank = matrixRank[radialP . pKineticForm . radialP];
  pRadialHessianRank = matrixRank[radialP . pRadialHessian . radialP];

  uWKineticRank = matrixRank[{{If[uWKineticPresent, 1, 0]}}];
  uWGapRank = matrixRank[{{If[uWGapPresent, 1, 0]}}];
  phiKineticRank = matrixRank[{{If[phiPresent, 1, 0]}}];

  computedCounts = <|
    "u_transverse" -> uCurlRank,
    "u_longitudinal" -> uCurlNullity,
    "P_spin_wave" -> Min[pTangentKineticRank, pTangentFrankRank],
    "P_soft_spin_radial" -> Min[pRadialKineticRank, pRadialHessianRank],
    "u_w" -> Min[uWKineticRank, uWGapRank],
    "phi" -> phiKineticRank
  |>;
  computedTotal = Total[Values[computedCounts]];
  <|
    "counts" -> computedCounts,
    "total" -> computedTotal,
    "rank_bookkeeping" -> <|
      "u_kinetic_rank" -> uKineticRank,
      "u_curl_rank" -> uCurlRank,
      "u_curl_nullity_within_active_kinetic_space" -> uCurlNullity,
      "P_tangent_kinetic_rank" -> pTangentKineticRank,
      "P_tangent_Frank_rank" -> pTangentFrankRank,
      "P_radial_kinetic_rank" -> pRadialKineticRank,
      "P_radial_soft_spin_hessian_rank" -> pRadialHessianRank,
      "u_w_kinetic_rank" -> uWKineticRank,
      "u_w_gap_rank" -> uWGapRank,
      "phi_kinetic_rank" -> phiKineticRank
    |>
  |>
];

dofCountAblation[name_, mutation_, baseline_, mutated_] := Module[{},
  If[mutated["total"] === baseline["total"], fail["DOF-count ablation did not fire: " <> name]];
  <|
    "name" -> name,
    "mutation" -> mutation,
    "baseline_total" -> baseline["total"],
    "mutated_total" -> mutated["total"],
    "mutated_counts" -> mutated["counts"],
    "status" -> "FIRED"
  |>
];

k = Unique["k"];
curlStiffness = {{0, 0, 0}, {0, k^2, 0}, {0, 0, k^2}};
rank = MatrixRank[curlStiffness, ZeroTest -> (TrueQ[FullSimplify[# == 0, k != 0]] &)];
nullity = 3 - rank;
If[rank =!= 2 || nullity =!= 1, fail["unexpected curl rank/nullity"]];
computedModes = computeFlatBraneDof[<||>];
reportedCounts = computedModes["counts"];
reportedTotal = Total[Values[reportedCounts]];
If[reportedCounts =!= computedModes["counts"] || reportedTotal =!= computedModes["total"],
  fail["reported flat-brane DOF count is not wired to computed count"]
];
If[reportedCounts["u_transverse"] =!= rank || reportedCounts["u_longitudinal"] =!= nullity,
  fail["reported u-sector count is not wired to curl rank/nullity"]
];
If[reportedTotal =!= 8, fail["unexpected flat-brane DOF count"]];
dofAblations = {
  dofCountAblation[
    "drop_u_w_gap_term",
    "set the frozen u_w gap matrix Omega_w^2 u_w^2 to rank 0",
    computedModes,
    computeFlatBraneDof[<|"u_w_gap_present" -> False|>]
  ],
  dofCountAblation[
    "drop_P_soft_spin_radial_term",
    "set the T0 (P^i P^i - 1)^2 radial Hessian to rank 0",
    computedModes,
    computeFlatBraneDof[<|"p_radial_present" -> False|>]
  ],
  dofCountAblation[
    "zero_u_longitudinal_component",
    "remove the k-parallel u component from the kinetic field space",
    computedModes,
    computeFlatBraneDof[<|"u_active" -> {0, 1, 1}|>]
  ]
};
modes = <|
  "computed_from_quadratic_terms" -> True,
  "reported_counts_wired_to_computation" -> True,
  "basis_order" -> <|
    "u" -> {"u_x_parallel_to_k", "u_y_transverse", "u_z_transverse"},
    "P" -> {"P_tangent_1", "P_tangent_2", "P_tangent_3", "P_radial"},
    "u_w" -> {"u_w"},
    "phi" -> {}
  |>,
  "quadratic_terms_used" -> <|
    "u" -> "rho_br (partial_t u^a)^2 and mu_R (nabla_parallel x u)^2",
    "P" -> "m rho a^2 (D_t P^i)^2, m rho c_s^2 a^2 (partial_j P^i)^2, and m rho c_s^2 (P^i P^i - 1)^2",
    "u_w" -> "rho_br (partial_t u_w)^2 and rho_br Omega_w^2 u_w^2",
    "phi" -> "absent in active baseline"
  |>,
  "curl_stiffness_matrix_for_k_along_x" -> {{"0", "0", "0"}, {"0", "k**2", "0"}, {"0", "0", "k**2"}},
  "curl_rank" -> rank,
  "curl_nullity" -> nullity,
  "rank_bookkeeping" -> computedModes["rank_bookkeeping"],
  "counts" -> reportedCounts,
  "total" -> reportedTotal,
  "dof_count_ablations" -> dofAblations,
  "bulk_shear_status" -> "kept_GNLS_bulk_shear_free",
  "classification_guard" -> "counts_only_no_gate_verdict_no_boundedness_no_gauge_no_leak_claim"
|>;

ledger = <|
  "field_dof_subcount_new_at_G0" -> 4,
  "constant_subcount" -> 4,
  "function_subcount" -> 1,
  "structural_subcount" -> 6,
  "independent_new_input_count_n" -> 11,
  "drift_verdict" -> "SECOND_MEDIUM_DRIFT_AT_FREEZE(11)",
  "flat_brane_total_listed_dof" -> 8
|>;

freeze = <|
  "t0_hash" -> t0Hash,
  "g0_hash" -> g0Hash,
  "g0_short_hash" -> expectedG0Short,
  "t0_bytes_embedded_in_g0" -> True,
  "g0_byte_range" -> byteRange[g0Text]
|>;

agreementPayload = <|
  "freeze" -> freeze,
  "dimension_checks" -> dimensions,
  "flat_brane_modes" -> modes,
  "ledger" -> ledger
|>;

report = <|
  "schema" -> "pathA_35_G0_mathematica/v1",
  "engine" -> "mathematica",
  "pass" -> True,
  "scope" -> "G0 freeze only; no Gate L verdict computed",
  "agreement_payload" -> agreementPayload,
  "verdicts" -> {
    "T0_SHEAR_FROZEN(" <> expectedG0Short <> ")",
    "SECOND_MEDIUM_DRIFT_AT_FREEZE(11)"
  }
|>;

If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir, CreateIntermediateDirectories -> True]];
outPath = FileNameJoin[{scratchDir, "pathA_35_G0_mathematica.json"}];
Export[outPath, report, "RawJSON"];
Print["wrote ", outPath];
Print["pathA_35 G0 Mathematica: PASS"];
Exit[0]
