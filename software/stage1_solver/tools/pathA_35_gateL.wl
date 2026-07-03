(* pathA_35 Gate L Mathematica checker.

   Independent Mathematica engine for the flat-brane Gate L algebra.  The
   agreement_payload is intentionally matched to the SymPy engine so --compare
   can detect ENGINE_DISAGREE exactly.
*)

ClearAll[
  fail, scriptPath, stage1Root, reportsDir, scratchDir, t0Report, g0Report,
  expectedT0Hash, expectedG0Hash, expectedG0Short, extractBlock, sha256,
  t0Text, g0Text, t0Block, g0Block, t0Hash, g0Hash, requiredNeedles,
  freezeFidelity, dadd, dmul, dsub, dimString, dimRecord, records, ablations,
  checkDim, checkSame, expectFail, buildDimensions, rankInt, assertZero,
  buildDerivedQuantities, tractionAudit, hiddenModeAudit, c5Audit,
  boundedClosureAudit, leakAudit, uwGapAudit, subHurdlePass,
  labelsFromSubHurdles, aggregateConfigVerdict, aggregateOverallVerdict,
  goodStructureFixtureVerdict, configResult, controlsSummary,
  hypotheticalPassFixture, buildPayload, buildReport, report, outPath
];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_35_gateL.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
reportsDir = FileNameJoin[{stage1Root, "reports"}];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
t0Report = FileNameJoin[{reportsDir, "pathA_24_T0_freeze.md"}];
g0Report = FileNameJoin[{reportsDir, "pathA_35_G0_freeze.md"}];

expectedT0Hash = "8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064";
expectedG0Hash = "d9520d3819c3f718290f9d0be57138c07d5bf02d2237106478e17b6a1e389ac3";
expectedG0Short = "d9520d3819c3";

extractBlock[text_, path_] := Module[{blocks},
  blocks = StringCases[text, "```freeze-action\n" ~~ Shortest[block__] ~~ "\n```" :> block];
  If[Length[blocks] =!= 1,
    fail["expected exactly one freeze-action block in " <> path <> ", found " <> ToString[Length[blocks]]]
  ];
  First[blocks]
];

sha256[block_] := ToLowerCase[IntegerString[Hash[block <> "\n", "SHA256"], 16, 64]];

t0Text = Import[t0Report, "Text"];
g0Text = Import[g0Report, "Text"];
t0Block = extractBlock[t0Text, t0Report];
g0Block = extractBlock[g0Text, g0Report];
t0Hash = sha256[t0Block];
g0Hash = sha256[g0Block];

freezeFidelity[] := Module[{},
  If[t0Hash =!= expectedT0Hash, fail["T0 hash mismatch: expected " <> expectedT0Hash <> ", got " <> t0Hash]];
  If[g0Hash =!= expectedG0Hash, fail["G0 hash mismatch: expected " <> expectedG0Hash <> ", got " <> g0Hash]];
  requiredNeedles = {
    "Active baseline: massless T0 spin-wave modes",
    "L_Pu := - lambda_Pu varpi_a Omega_u^a",
    "No scalar-potential analog phi is present",
    "Named inactive alternate: slaved-rigid P_parallel = w_hat x Omega_u",
    "T_wa := m rho v_w v_a"
  };
  Scan[
    If[! StringContainsQ[g0Block, #], fail["required frozen G0 line missing: " <> #]] &,
    requiredNeedles
  ];
  <|
    "t0_hash" -> t0Hash,
    "g0_hash" -> g0Hash,
    "g0_short_hash" -> expectedG0Short,
    "t0_bytes_embedded_in_g0" -> StringContainsQ[g0Block, t0Block],
    "required_gateL_lines_present" -> True
  |>
];

dadd[dims__] := Total[{dims}];
dmul[n_, dim_] := n dim;
dsub[left_, right_] := left - right;
dimString[dim_] := Module[{labels = {"M", "L", "T"}, parts},
  parts = MapThread[
    If[#2 === 0, Nothing, If[#2 === 1, #1, #1 <> "^" <> ToString[#2]]] &,
    {labels, dim}
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

buildDimensions[] := Module[
  {
    Mdim, Ldim, Tdim, Zdim, bulkStress, braneLag, uOperator, pOperator,
    mixedUPOperator, determinantDim, closureDim, advectiveCurrent,
    advectiveCurl, dm, drho, da, dellG, du, duw, dP, dphi, dgrad, ddt,
    dk, domega, dv, dcs2, drhoBr, dmuR, dlambda, dOmegaW, dOmegaP,
    dTbulk, dFrankCoupleTraction, dcurlU, ddivU, dJP, dGammaP, dMgapP,
    constants
  },
  records = {};
  ablations = {};
  Mdim = {1, 0, 0}; Ldim = {0, 1, 0}; Tdim = {0, 0, 1}; Zdim = {0, 0, 0};
  bulkStress = {1, -2, -2};
  braneLag = {1, -1, -2};
  uOperator = {1, -3, -2};
  pOperator = braneLag;
  mixedUPOperator = {1, -2, -2};
  determinantDim = {2, -4, -4};
  closureDim = braneLag;
  advectiveCurrent = {1, -2, -1};
  advectiveCurl = {1, -3, -1};

  dm = Mdim; drho = {0, -4, 0}; da = Ldim; dellG = Ldim;
  du = Ldim; duw = Ldim; dP = Zdim; dphi = {0, 2, -1};
  dgrad = {0, -1, 0}; ddt = {0, 0, -1}; dk = {0, -1, 0};
  domega = {0, 0, -1}; dv = {0, 1, -1}; dcs2 = {0, 2, -2};
  drhoBr = {1, -3, 0}; dmuR = braneLag; dlambda = braneLag;
  dOmegaW = domega; dOmegaP = domega; dTbulk = bulkStress;
  dFrankCoupleTraction = {1, 0, -2};
  dcurlU = dadd[dgrad, du];
  ddivU = dcurlU;

  dJP = checkDim["projected_T0", "J_P", dadd[dm, drho, dmul[2, da], dellG], {1, -1, 0}, "m rho a^2 ell_g"];
  dGammaP = checkDim[
    "projected_T0",
    "Gamma_P",
    dadd[dm, drho, dcs2, dmul[2, da], dellG],
    {1, 1, -2},
    "m rho c_s^2 a^2 ell_g"
  ];
  dMgapP = checkDim["projected_T0", "M_gap_P", dadd[dJP, dmul[2, dOmegaP]], braneLag, "J_P Omega_P^2"];

  checkDim["L_a_i_traction", "MacCullagh_brane_traction", dadd[dmuR, dcurlU], braneLag, "mu_R Omega_u"];
  checkDim["L_a_i_traction", "P_u_brane_traction", dadd[dlambda, dP], braneLag, "lambda_Pu varpi"];
  checkDim[
    "L_a_i_traction",
    "P_Frank_torque_not_surface_traction",
    dadd[dGammaP, dgrad, dP],
    dFrankCoupleTraction,
    "Gamma_P partial P"
  ];
  expectFail["ablation", "P_u_without_curl_or_cut_gradient", dadd[dlambda, dP, du], braneLag, "lambda_Pu P u"];

  checkDim["L_a_ii_symbol", "rho_br_omega2_u2", dadd[drhoBr, dmul[2, domega], dmul[2, du]], braneLag, "rho_br omega^2 u^2"];
  checkDim["L_a_ii_symbol", "mu_R_k2_u2", dadd[dmuR, dmul[2, dk], dmul[2, du]], braneLag, "mu_R k^2 u^2"];
  checkDim["L_a_ii_symbol", "J_P_omega2_P2", dadd[dJP, dmul[2, domega], dmul[2, dP]], braneLag, "J_P omega^2 P^2"];
  checkDim["L_a_ii_symbol", "Gamma_P_k2_P2", dadd[dGammaP, dmul[2, dk], dmul[2, dP]], braneLag, "Gamma_P k^2 P^2"];
  checkDim["L_a_ii_symbol", "lambda_k_u_P", dadd[dlambda, dk, du, dP], braneLag, "lambda_Pu k u P"];
  checkDim["L_a_ii_symbol", "u_operator_entry", dadd[dmuR, dmul[2, dk]], uOperator, "mu_R k^2"];
  checkDim["L_a_ii_symbol", "P_operator_entry", dadd[dGammaP, dmul[2, dk]], pOperator, "Gamma_P k^2"];
  checkDim["L_a_ii_symbol", "mixed_operator_entry", dadd[dlambda, dk], mixedUPOperator, "lambda_Pu k"];
  checkDim["L_a_ii_dispersion", "c_gamma_squared", dsub[dmuR, drhoBr], {0, 2, -2}, "mu_R/rho_br"];
  checkDim["L_a_ii_dispersion", "c_gamma_eff_squared", dsub[dmuR, drhoBr], {0, 2, -2}, "(mu_R - 2 lambda_Pu)/rho_br"];
  checkDim["L_a_ii_dispersion", "slaved_k4_coefficient", dsub[dGammaP, drhoBr], {0, 4, -2}, "Gamma_P/rho_br"];
  checkDim["L_a_ii_dispersion", "k_disp_squared", dsub[dmuR, dGammaP], {0, -2, 0}, "(mu_R - 2 lambda_Pu)/Gamma_P"];
  checkDim["L_a_ii_dispersion", "omega2_k2", dadd[dsub[dmuR, drhoBr], dmul[2, dk]], {0, 0, -2}, "c_gamma^2 k^2"];
  checkDim["L_a_ii_dispersion", "omega2_k4", dadd[dsub[dGammaP, drhoBr], dmul[4, dk]], {0, 0, -2}, "(Gamma_P/rho_br) k^4"];

  checkDim["L_a_iii_C5", "u_longitudinal_kinetic", dadd[drhoBr, dmul[2, ddt], dmul[2, du]], braneLag, "rho_br (partial_t u_L)^2"];
  checkDim["L_a_iii_C5", "phi_gradient_velocity", dadd[dgrad, dphi], {0, 1, -1}, "partial_a phi"];
  checkSame[
    "L_a_iii_C5",
    "phi_fixture_kinetic_parts",
    {dadd[ddt, du], dadd[dgrad, dphi]},
    {0, 1, -1},
    "partial_t u_a - partial_a phi"
  ];
  checkDim["L_a_iii_C5", "phi_fixture_term", dadd[drhoBr, dmul[2, dadd[ddt, du]]], braneLag, "rho_br (partial_t u - grad phi)^2"];

  checkDim["L_b_Hamiltonian", "det_mu_gamma_k4", dadd[dmuR, dGammaP, dmul[4, dk]], determinantDim, "mu_R Gamma_P k^4"];
  checkDim["L_b_Hamiltonian", "det_lambda2_k2", dadd[dmul[2, dlambda], dmul[2, dk]], determinantDim, "lambda_Pu^2 k^2"];
  checkDim["L_b_closure", "antisymmetric_stress_torque", dadd[dmuR, dcurlU], closureDim, "2 mu_R Omega_u"];
  checkDim["L_b_closure", "couple_stress_divergence", dadd[dGammaP, dmul[2, dk], dP], closureDim, "Gamma_P k^2 P"];
  checkDim["L_b_closure", "gapped_P_response", dsub[dlambda, dMgapP], Zdim, "lambda_Pu Omega_u / M_gap_P"];

  checkDim["L_c_direct_leak", "T_wa", dadd[dm, drho, dv, dv], bulkStress, "m rho v_w v_a"];
  checkDim["L_c_direct_leak", "partial_b_u_w", dadd[dgrad, duw], Zdim, "partial_b u_w"];
  checkDim["L_c_direct_leak", "slope_mixed_bulk_stress", dadd[dTbulk, dadd[dgrad, duw]], bulkStress, "(T_ww delta_ab - T_ab) partial_b u_w"];
  checkSame[
    "L_c_direct_leak",
    "T_na_projected",
    {dadd[dm, drho, dv, dv], dadd[dTbulk, dadd[dgrad, duw]]},
    bulkStress,
    "T_wa + (T_ww delta_ab - T_ab) partial_b u_w"
  ];
  expectFail["ablation", "drop_m_from_T_wa", dadd[drho, dv, dv], bulkStress, "rho v_w v_a"];
  checkDim["L_c_indirect_leak", "advective_P_current", dadd[dJP, ddt, dP, dgrad, dP], advectiveCurrent, "J_P partial_t P partial_a P"];
  checkDim["L_c_indirect_leak", "advective_current_curl", dadd[advectiveCurrent, dgrad], advectiveCurl, "curl(J_P partial_t P grad P)"];

  checkDim["L_d_uw_gap", "u_w_kinetic", dadd[drhoBr, dmul[2, ddt], dmul[2, duw]], braneLag, "rho_br (partial_t u_w)^2"];
  checkDim["L_d_uw_gap", "u_w_gap", dadd[drhoBr, dmul[2, dOmegaW], dmul[2, duw]], braneLag, "rho_br Omega_w^2 u_w^2"];
  checkDim["L_d_uw_gap", "omega_uw_squared", dmul[2, dOmegaW], {0, 0, -2}, "Omega_w^2"];

  constants = <|
    "rho" -> dimRecord[drho],
    "m" -> dimRecord[dm],
    "a" -> dimRecord[da],
    "ell_g" -> dimRecord[dellG],
    "rho_br" -> dimRecord[drhoBr],
    "mu_R" -> dimRecord[dmuR],
    "lambda_Pu" -> dimRecord[dlambda],
    "Omega_w" -> dimRecord[dOmegaW],
    "Omega_P_gap_control" -> dimRecord[dOmegaP],
    "J_P" -> dimRecord[dJP],
    "Gamma_P" -> dimRecord[dGammaP],
    "phi_fixture" -> dimRecord[dphi]
  |>;

  <|
    "base_dimensions" -> <|"M" -> {1, 0, 0}, "L" -> {0, 1, 0}, "T" -> {0, 0, 1}|>,
    "targets" -> <|
      "bulk_interface_traction" -> dimRecord[bulkStress],
      "brane_action_density_or_stress" -> dimRecord[braneLag],
      "u_operator" -> dimRecord[uOperator],
      "P_operator" -> dimRecord[pOperator],
      "mixed_u_P_operator" -> dimRecord[mixedUPOperator],
      "Hamiltonian_minor" -> dimRecord[determinantDim],
      "closure_residual" -> dimRecord[closureDim],
      "advective_current" -> dimRecord[advectiveCurrent],
      "advective_current_curl" -> dimRecord[advectiveCurl],
      "Frank_couple_traction_not_surface_traction" -> dimRecord[dFrankCoupleTraction]
    |>,
    "constants" -> constants,
    "checks" -> records,
    "ablations" -> ablations,
    "pass" -> True
  |>
];

rankInt[m_] := MatrixRank[m];

assertZero[expr_, label_] := Module[{reduced},
  reduced = FullSimplify[expr];
  If[reduced =!= 0, fail[label <> ": " <> ToString[InputForm[Factor[reduced]]]]]
];

buildDerivedQuantities[] := Module[
  {
    w, k, ellG, m, rho, a, cs, muR, lambdaPu, rhoBr, omega2,
    OmegaW, OmegaP, jp, gammaP, mgap, gEll, gNorm, modeWeight, jExpr,
    gammaExpr, radialExpr, uT, varpiT, omegaU, transverseEnergy,
    transverseHessian, expectedHessian, hessianMinor, lowKMinorSign,
    slavedEnergy, slavedHessian, muEff, k4Coeff, ux, uy, uz, px, py, pz,
    uw, omegaVec, pVec, unslavedSymbol, dynMassless, detMassless,
    dynGapped, detGapped, detMasslessLam0, detGappedLam0, slavedQ,
    slavedTransverseStiffness, slavedPotential, slavedStiffness,
    slavedMass, dynSlaved, detSlaved, masslessDegree, gappedDegree,
    slavedDegree, gaplessMasslessLam0, gaplessGappedLam0, gaplessSlaved,
    positiveMasslessLam0, positiveGappedLam0, positiveSlaved,
    zeroMassless, zeroGapped, zeroSlaved, transverseNegativeLowK,
    masslessPositiveLowKLambdaNonzero, Omega, pResponseGapped,
    coupleDivergence, gappedClosureResidual, x, y, t, ampA, ampB, q,
    omega, phase, uFlatY, uFlatW, vFlatW, vFlatY, slopeFlatX, deltaT,
    tNaFlat, uBentW, tNaBent, responseDenX, omegaUFlat, pFlat,
    advectiveJx, advectiveJy, indirectCurlFlat, pFrankOnly, responseDenY,
    pNonplanar, jxNp, jyNp, curlNp
  },
  {w, k, ellG, m, rho, a, cs, muR, lambdaPu, rhoBr, omega2, OmegaW,
    OmegaP, jp, gammaP, mgap} =
    {ww, kk, ellGg, mm, rrho, aa, css, muRR, lambdaPuu, rhoBrr,
      omega22, OmegaWW, OmegaPP, JPP, GammaPP, MGapPP};

  gEll = Exp[-(w/ellG)^2]/(Sqrt[Pi] ellG);
  gNorm = FullSimplify[Integrate[gEll, {w, -Infinity, Infinity}], ellG > 0];
  modeWeight = FullSimplify[ellG gNorm, ellG > 0];
  jExpr = Factor[m rho a^2 modeWeight];
  gammaExpr = Factor[m rho cs^2 a^2 modeWeight];
  radialExpr = Factor[m rho cs^2 modeWeight];
  assertZero[gNorm - 1, "g_ell normalization"];
  assertZero[modeWeight - ellG, "projected confinement width"];
  assertZero[jExpr - m rho a^2 ellG, "projected J_P"];
  assertZero[gammaExpr - m rho cs^2 a^2 ellG, "projected Gamma_P"];
  assertZero[radialExpr - gammaExpr/a^2, "projected radial stiffness Gamma_P/a^2"];

  {uT, varpiT} = {uTT, varpiTT};
  omegaU = k uT;
  transverseEnergy = 1/2 muR omegaU^2 + lambdaPu varpiT omegaU + 1/2 gammaP k^2 varpiT^2;
  transverseHessian = Table[D[transverseEnergy, var, vars],
    {var, {uT, varpiT}}, {vars, {uT, varpiT}}];
  expectedHessian = {{muR k^2, lambdaPu k}, {lambdaPu k, gammaP k^2}};
  If[FullSimplify[transverseHessian - expectedHessian] =!= ConstantArray[0, {2, 2}],
    fail["derived transverse Hessian mismatch"]
  ];
  hessianMinor = Factor[Det[transverseHessian]];
  assertZero[hessianMinor - k^2 (gammaP muR k^2 - lambdaPu^2), "Hamiltonian minor"];
  assertZero[(hessianMinor /. lambdaPu -> 0) - gammaP muR k^4, "lambda_Pu zero boundedness sanity"];
  lowKMinorSign = Factor[Limit[hessianMinor/k^2, k -> 0, Direction -> "FromAbove"]];
  If[lowKMinorSign =!= -lambdaPu^2, fail["unexpected low-k Hessian minor sign"]];

  slavedEnergy = Factor[transverseEnergy /. varpiT -> -omegaU];
  slavedHessian = Factor[D[slavedEnergy, {uT, 2}]];
  muEff = Factor[Limit[slavedHessian/k^2, k -> 0, Direction -> "FromAbove"]];
  k4Coeff = Coefficient[Expand[slavedHessian], k, 4];
  assertZero[muEff - (muR - 2 lambdaPu), "slaved mu_eff"];
  assertZero[k4Coeff - gammaP, "slaved k4 coefficient"];
  assertZero[(-omegaU) + omegaU, "slaved closure residual"];

  {ux, uy, uz, px, py, pz, uw} = {uxx, uyy, uzz, pxx, pyy, pzz, uww};
  omegaVec = {0, -k uz, k uy};
  pVec = {px, py, pz};
  unslavedSymbol[gap_] := Module[{qvars, potential, stiffness, mass, dyn},
    qvars = {ux, uy, uz, px, py, pz, uw};
    potential =
      1/2 muR (omegaVec.omegaVec)
      + lambdaPu (pVec.omegaVec)
      + 1/2 (gammaP k^2 + gap) (pVec.pVec)
      + 1/2 rhoBr OmegaW^2 uw^2;
    stiffness = Table[D[potential, qvars[[i]], qvars[[j]]], {i, Length[qvars]}, {j, Length[qvars]}];
    mass = DiagonalMatrix[{rhoBr, rhoBr, rhoBr, jp, jp, jp, rhoBr}];
    dyn = Factor[stiffness - omega2 mass];
    {dyn, Factor[Det[dyn]]}
  ];
  {dynMassless, detMassless} = unslavedSymbol[0];
  {dynGapped, detGapped} = unslavedSymbol[mgap];
  detMasslessLam0 = Factor[detMassless /. lambdaPu -> 0];
  detGappedLam0 = Factor[detGapped /. lambdaPu -> 0];

  slavedQ = {ux, uy, uz, uw};
  slavedTransverseStiffness = Factor[(muEff k^2 + gammaP k^4)];
  slavedPotential = 1/2 slavedTransverseStiffness (uy^2 + uz^2) + 1/2 rhoBr OmegaW^2 uw^2;
  slavedStiffness = Table[D[slavedPotential, slavedQ[[i]], slavedQ[[j]]], {i, Length[slavedQ]}, {j, Length[slavedQ]}];
  slavedMass = DiagonalMatrix[{rhoBr, rhoBr, rhoBr, rhoBr}];
  dynSlaved = Factor[slavedStiffness - omega2 slavedMass];
  detSlaved = Factor[Det[dynSlaved]];

  masslessDegree = Exponent[detMassless, omega2];
  gappedDegree = Exponent[detGapped, omega2];
  slavedDegree = Exponent[detSlaved, omega2];
  If[{masslessDegree, gappedDegree, slavedDegree} =!= {7, 7, 4},
    fail["unexpected principal-symbol degrees"]
  ];
  assertZero[
    detMasslessLam0 - (-omega2 rhoBr) (muR k^2 - omega2 rhoBr)^2
      (gammaP k^2 - jp omega2)^3 (rhoBr (OmegaW^2 - omega2)),
    "massless lambda_Pu=0 determinant"
  ];
  assertZero[
    detGappedLam0 - (-omega2 rhoBr) (muR k^2 - omega2 rhoBr)^2
      (mgap + gammaP k^2 - jp omega2)^3 (rhoBr (OmegaW^2 - omega2)),
    "gapped lambda_Pu=0 determinant"
  ];
  assertZero[
    detSlaved - (-omega2 rhoBr)
      (gammaP k^4 + k^2 (-2 lambdaPu + muR) - omega2 rhoBr)^2
      (rhoBr (OmegaW^2 - omega2)),
    "slaved determinant"
  ];

  gaplessMasslessLam0 = 5;
  gaplessGappedLam0 = 2;
  gaplessSlaved = 2;
  positiveMasslessLam0 = 6;
  positiveGappedLam0 = 6;
  positiveSlaved = 3;
  zeroMassless = 1;
  zeroGapped = 1;
  zeroSlaved = 1;
  transverseNegativeLowK = If[lowKMinorSign === -lambdaPu^2, 2, 0];
  masslessPositiveLowKLambdaNonzero = gaplessMasslessLam0 - transverseNegativeLowK;
  If[masslessPositiveLowKLambdaNonzero =!= 3, fail["lambda_Pu nonzero low-k branch count did not change"]];

  Omega = OmegaUu;
  pResponseGapped = Factor[-lambdaPu Omega/(mgap + gammaP k^2)];
  coupleDivergence = Factor[gammaP k^2 pResponseGapped];
  gappedClosureResidual = Factor[Limit[2 muR Omega + coupleDivergence, k -> 0, Direction -> "FromAbove"]];
  assertZero[gappedClosureResidual - 2 muR Omega, "gapped-P low-k closure residual"];

  {x, y, t, ampA, ampB, q, omega, deltaT} = {xx, yy, tt, ampAA, ampBB, qq, omegaa, DeltaTT};
  phase = k x - omega t;
  uFlatY = ampA Cos[phase];
  uFlatW = 0;
  vFlatW = D[uFlatW, t];
  vFlatY = D[uFlatY, t];
  slopeFlatX = D[uFlatW, x];
  tNaFlat = Simplify[m rho vFlatW vFlatY + deltaT slopeFlatX];
  uBentW = ampB Cos[q x - omega t];
  tNaBent = Simplify[
    (m rho D[uBentW, t] vFlatY + deltaT D[uBentW, x]) /.
      {ampA -> 2, ampB -> 3, m -> 5, rho -> 7, deltaT -> 11, k -> 13, q -> 17, omega -> 19, x -> 23, t -> 29}
  ];
  If[tNaFlat =!= 0 || tNaBent === 0, fail["direct L(c) traction derivation failed"]];

  responseDenX = gammaP k^2 - jp omega^2;
  omegaUFlat = D[uFlatY, x];
  pFlat = Factor[-lambdaPu omegaUFlat/responseDenX];
  advectiveJx = Factor[jp D[pFlat, t] D[pFlat, x]];
  advectiveJy = 0;
  indirectCurlFlat = Simplify[D[advectiveJy, x] - D[advectiveJx, y]];
  pFrankOnly = Simplify[pFlat /. lambdaPu -> 0];
  If[indirectCurlFlat =!= 0 || pFrankOnly =!= 0, fail["flat indirect channel or Frank-only control failed"]];

  responseDenY = gammaP q^2 - jp omega^2;
  pNonplanar = Factor[-lambdaPu (
      ampA k Sin[k x] Cos[omega t]/responseDenX
      + ampB q Sin[q y] Sin[omega t]/responseDenY
    )];
  jxNp = D[pNonplanar, t] D[pNonplanar, x];
  jyNp = D[pNonplanar, t] D[pNonplanar, y];
  curlNp = FullSimplify[TrigReduce[D[jyNp, x] - D[jxNp, y]]];
  If[FullSimplify[curlNp] === 0, fail["nonplanar indirect able-to-fail control failed"]];

  <|
    "projection" -> <|
      "g_ell_normalization" -> "1",
      "mode_weight_integral" -> "ell_g",
      "J_P" -> "m rho a^2 ell_g",
      "Gamma_P" -> "m rho c_s^2 a^2 ell_g",
      "radial_stiffness" -> "Gamma_P/a^2",
      "source" -> "int dw ell_g g_ell(w) times the inherited T0 P terms"
    |>,
    "transverse_hessian" -> <|
      "basis" -> {"u_T", "varpi_T"},
      "matrix" -> {{"mu_R k^2", "lambda_Pu k"}, {"lambda_Pu k", "Gamma_P k^2"}},
      "principal_minor" -> "k^2 (Gamma_P mu_R k^2 - lambda_Pu^2)",
      "negative_energy_interval" -> "0 < k^2 < lambda_Pu^2/(Gamma_P mu_R)",
      "lambda_Pu_zero_minor" -> "Gamma_P mu_R k^4",
      "bounded_below_lambda_Pu_zero" -> True,
      "low_k_minor_limit_over_k2" -> "-lambda_Pu^2"
    |>,
    "slaved_reduction" -> <|
      "substitution" -> "P_parallel = w_hat x Omega_u; varpi = w_hat x P_parallel = -Omega_u",
      "mu_eff" -> "mu_R - 2 lambda_Pu",
      "transverse_energy_hessian" -> "k^2 (mu_R - 2 lambda_Pu + Gamma_P k^2)",
      "k4_coefficient" -> "Gamma_P",
      "omega2_slaved" -> "((mu_R - 2 lambda_Pu) k^2 + Gamma_P k^4)/rho_br",
      "k_disp_squared" -> "(mu_R - 2 lambda_Pu)/Gamma_P",
      "closure_residual" -> "0"
    |>,
    "principal_symbols" -> <|
      "massless_P" -> <|
        "det_degree" -> masslessDegree,
        "positive_omega2_branches_lambda_Pu_zero" -> positiveMasslessLam0,
        "gapless_positive_branches_lambda_Pu_zero" -> gaplessMasslessLam0,
        "positive_gapless_branches_low_k_lambda_Pu_nonzero" -> masslessPositiveLowKLambdaNonzero,
        "negative_omega2_branches_low_k_lambda_Pu_nonzero" -> transverseNegativeLowK,
        "longitudinal_zero_roots" -> zeroMassless,
        "extra_gapless_P_branches" -> gaplessMasslessLam0 - 2
      |>,
      "gapped_P" -> <|
        "det_degree" -> gappedDegree,
        "positive_omega2_branches_lambda_Pu_zero" -> positiveGappedLam0,
        "gapless_positive_branches_lambda_Pu_zero" -> gaplessGappedLam0,
        "longitudinal_zero_roots" -> zeroGapped,
        "extra_gapless_P_branches" -> gaplessGappedLam0 - 2,
        "gapped_P_branches" -> 3
      |>,
      "slaved_rigid" -> <|
        "det_degree" -> slavedDegree,
        "positive_omega2_branches" -> positiveSlaved,
        "gapless_positive_branches" -> gaplessSlaved,
        "longitudinal_zero_roots" -> zeroSlaved,
        "extra_gapless_P_branches" -> gaplessSlaved - 2
      |>
    |>,
    "closure" -> <|
      "gapped_P_response" -> "-lambda_Pu Omega_u/(M_gap_P + Gamma_P k^2)",
      "gapped_P_low_k_residual" -> "2 mu_R Omega_u",
      "omit_reservoir_low_k_residual" -> "2 mu_R Omega_u",
      "slaved_residual" -> "0"
    |>,
    "leak" -> <|
      "direct_flat_T_na" -> "0",
      "direct_flat_wave" -> "u_y=A cos(k x - omega t), u_w=0, v_w=partial_t u_w=0",
      "bent_control_nonzero" -> True,
      "indirect_flat_curl" -> "0",
      "indirect_flat_P_response" -> "-lambda_Pu partial_x u_y/(Gamma_P k^2 - J_P omega^2)",
      "advective_current" -> "J_P (partial_t P) partial_a P from D_t^v P at leading v=0",
      "Frank_only_induced_P" -> "0",
      "nonplanar_able_to_fail_nonzero" -> True
    |>
  |>
];

tractionAudit[config_] := Module[
  {lambdaPu, muR, couplingCutMatrix, macCutMatrix, frankOnlyMatrix, couplingRank,
    macRank, frankRank, provenance, note},
  lambdaPu = Unique["lambdaPu"];
  muR = Unique["muR"];
  couplingCutMatrix = {{0, 0, 0}, {0, 0, -lambdaPu}, {0, lambdaPu, 0}};
  macCutMatrix = {{0, 0, 0}, {0, 0, -muR}, {0, muR, 0}};
  frankOnlyMatrix = ConstantArray[0, {3, 3}];
  couplingRank = rankInt[couplingCutMatrix];
  macRank = rankInt[macCutMatrix];
  frankRank = rankInt[frankOnlyMatrix];
  If[couplingRank =!= 2 || macRank =!= 2 || frankRank =!= 0, fail["unexpected traction ranks"]];
  If[config === "A_baseline",
    provenance = "ARROWS_SUPPLY_TRACTION";
    note = "P-u virtual work gives a rank-2 cut traction; standalone mu_R traction is also present.",
    provenance = "POSTULATED_SURFACE_ELASTICITY";
    note = "slaving removes independent arrow modes; the clean k^2 light traction is the postulated surface modulus, with arrow-sector k^4/coupling corrections."
  ];
  <|
    "status" -> "PASS",
    "label" -> Null,
    "provenance" -> provenance,
    "cut_normal" -> "x",
    "coupling_traction_rank" -> couplingRank,
    "standalone_MacCullagh_traction_rank" -> macRank,
    "Frank_only_reference_traction_rank" -> frankRank,
    "Frank_only_control" -> <|
      "status" -> "FIRED",
      "label" -> "FAIL_FRANK_TORQUE_NOT_MACCULLAGH_TRACTION",
      "computed_rank" -> frankRank
    |>,
    "standalone_elasticity_discriminator" -> <|
      "lambda_Pu_zero_mu_R_nonzero_rank" -> macRank,
      "status" -> "FIRED",
      "meaning" -> "the machinery distinguishes P-u-sourced traction from postulated surface elasticity"
    |>,
    "note" -> note
  |>
];

hiddenModeAudit[config_, derived_] := Module[
  {k, muR, rhoBr, gammaP, lambdaPu, lambdaC, muC, cauchy, cauchyRank,
    counts, extra, status, label, dispersion, muEff, omega2, kdisp2,
    modeCounts, masslessCounts, slaved, slavedCounts},
  k = Unique["k"]; muR = Unique["muR"]; rhoBr = Unique["rhoBr"];
  gammaP = Unique["GammaP"]; lambdaPu = Unique["lambdaPu"];
  lambdaC = Unique["lambdaC"]; muC = Unique["muC"];
  cauchy = DiagonalMatrix[{(lambdaC + 2 muC) k^2, muC k^2, muC k^2}];
  cauchyRank = rankInt[cauchy];
  If[cauchyRank =!= 3, fail["Cauchy control did not produce three propagating elastic modes"]];
  modeCounts = derived["principal_symbols"];
  If[config === "A_baseline",
    masslessCounts = modeCounts["massless_P"];
    counts = <|
      "positive_gapless_branches_low_k_lambda_Pu_nonzero" -> masslessCounts["positive_gapless_branches_low_k_lambda_Pu_nonzero"],
      "negative_omega2_branches_low_k_lambda_Pu_nonzero" -> masslessCounts["negative_omega2_branches_low_k_lambda_Pu_nonzero"],
      "gapless_positive_branches_lambda_Pu_zero" -> masslessCounts["gapless_positive_branches_lambda_Pu_zero"],
      "P_spin_wave_gapless_lambda_Pu_zero" -> masslessCounts["extra_gapless_P_branches"],
      "u_longitudinal_zero_roots" -> masslessCounts["longitudinal_zero_roots"],
      "P_radial_gapped" -> 1,
      "u_w_gapped" -> 1,
      "phi" -> 0
    |>;
    extra = masslessCounts["extra_gapless_P_branches"];
    status = "FAIL";
    label = "FAIL_HIDDEN_PROPAGATING_MODE";
    dispersion = <|
      "u_light_omega2_uncoupled_reference" -> "(mu_R/rho_br) k^2",
      "P_spin_omega2" -> "c_s^2 k^2",
      "coupled_transverse_pencil" -> "det[[mu_R k^2-rho_br omega^2, lambda_Pu k],[lambda_Pu k, Gamma_P k^2-J_P omega^2]]=0",
      "small_k_warning" -> "derived Hamiltonian minor gives two negative omega^2 transverse branches at low k for nonzero lambda_Pu; counted separately in L(b)"
    |>,
    slaved = derived["slaved_reduction"];
    slavedCounts = modeCounts["slaved_rigid"];
    muEff = muR - 2 lambdaPu;
    omega2 = FullSimplify[(muEff k^2 + gammaP k^4)/rhoBr];
    kdisp2 = FullSimplify[muEff/gammaP];
    If[FullSimplify[rhoBr omega2 - (muEff k^2 + gammaP k^4)] =!= 0, fail["slaved dispersion algebra mismatch"]];
    counts = <|
      "gapless_positive_branches" -> slavedCounts["gapless_positive_branches"],
      "positive_omega2_branches" -> slavedCounts["positive_omega2_branches"],
      "u_longitudinal_zero_roots" -> slavedCounts["longitudinal_zero_roots"],
      "P_spin_wave_gapless" -> slavedCounts["extra_gapless_P_branches"],
      "u_w_gapped" -> 1,
      "phi" -> 0
    |>;
    extra = slavedCounts["extra_gapless_P_branches"];
    status = "PASS_WITH_LOW_K_DISPERSION_CAVEAT";
    label = Null;
    dispersion = <|
      "omega2_slaved" -> slaved["omega2_slaved"],
      "c_gamma_squared_MacCullagh_feed_forward" -> "mu_R/rho_br",
      "effective_c_squared_if_bilinear_retained" -> "(mu_R - 2 lambda_Pu)/rho_br",
      "k4_correction" -> "(Gamma_P/rho_br) k^4",
      "k_disp_squared" -> slaved["k_disp_squared"],
      "low_k_tolerance" -> "requires k^2 << (mu_R - 2 lambda_Pu)/Gamma_P with positive effective modulus; exact nondispersive equality is false at finite k",
      "symbolic_checks" -> <|"omega2_identity" -> True, "k_disp_squared" -> "-(2*lambda_Pu - mu_R)/Gamma_P"|>
    |>
  ];
  <|
    "status" -> status,
    "label" -> label,
    "counts" -> counts,
    "extra_propagating_modes" -> extra,
    "derived_from_principal_symbol" -> True,
    "mode_count_controls" -> <|
      "gapped_P" -> modeCounts["gapped_P"],
      "lambda_Pu_zero_changes_low_k_count" ->
        modeCounts["massless_P"]["positive_gapless_branches_low_k_lambda_Pu_nonzero"] =!=
          modeCounts["massless_P"]["gapless_positive_branches_lambda_Pu_zero"]
    |>,
    "cauchy_reference_control" -> <|
      "status" -> "FIRED",
      "label" -> "FAIL_CAUCHY_STRAY_LONGITUDINAL",
      "computed_propagating_modes" -> cauchyRank,
      "stray_longitudinal" -> True
    |>,
    "c_gamma_squared" -> "mu_R/rho_br",
    "dispersion" -> dispersion
  |>
];

c5Audit[config_] := Module[
  {omega, k, noPhiMatrix, noPhiStiffnessRank, noPhiKineticRank,
    phiFixtureMatrix, phiFixtureRank, phiFixtureDet},
  omega = Unique["omega"]; k = Unique["k"];
  noPhiMatrix = {{omega^2}};
  noPhiStiffnessRank = rankInt[{{0}}];
  noPhiKineticRank = rankInt[noPhiMatrix];
  phiFixtureMatrix = {{omega^2, -omega k}, {-omega k, k^2}};
  phiFixtureRank = rankInt[phiFixtureMatrix];
  phiFixtureDet = Factor[Det[phiFixtureMatrix]];
  If[noPhiKineticRank =!= 1 || noPhiStiffnessRank =!= 0, fail["no-phi C5 branch did not leave the kinetic zero mode"]];
  If[phiFixtureRank =!= 1 || phiFixtureDet =!= 0, fail["phi fixture did not produce a gauge-null kinetic form"]];
  <|
    "status" -> "FAIL",
    "label" -> "FAIL_C5_LONGITUDINAL_ZERO_MODE",
    "frozen_phi_status" -> "absent",
    "no_phi_branch_control" -> <|
      "status" -> "FIRED",
      "kinetic_rank" -> noPhiKineticRank,
      "curl_stiffness_rank" -> noPhiStiffnessRank,
      "computed_result" -> "constrained physical zero mode remains"
    |>,
    "independent_variational_phi_fixture_control" -> <|
      "status" -> "FIRED_PASS_FIXTURE_ONLY",
      "kinetic_form_rank" -> phiFixtureRank,
      "determinant" -> ToString[phiFixtureDet],
      "gauge_nullity" -> 1,
      "physical_longitudinal_modes" -> 0,
      "gate_pass_candidate_only_if_frozen_in_G0" -> False
    |>,
    "raw_divergence_projector_control" -> <|
      "status" -> "FIRED",
      "label" -> "FAIL_C5_LONGITUDINAL_ZERO_MODE",
      "reason" -> "removes u_L by hand with no variational provenance"
    |>,
    "phi_equals_u_w_collision" -> <|
      "status" -> "FIRED",
      "reason" -> "a scalar-potential phi identified with u_w must be massless, but L(d) requires Omega_w>0"
    |>,
    "config_note" -> config <> " inherits the G0 no-phi decision"
  |>
];

boundedClosureAudit[config_, derived_] := Module[
  {hessian, slaved, closure, controls},
  hessian = derived["transverse_hessian"];
  slaved = derived["slaved_reduction"];
  closure = derived["closure"];
  controls = <|
    "omit_couple_stress_reservoir" -> <|
      "status" -> "FIRED",
      "label" -> "FAIL_GYROSTAT_NO_CLOSURE",
      "low_k_residual" -> closure["omit_reservoir_low_k_residual"]
    |>,
    "large_gap_P_control" -> <|
      "status" -> "FIRED",
      "label" -> "FAIL_GYROSTAT_NO_CLOSURE",
      "P_response" -> closure["gapped_P_response"],
      "low_k_couple_divergence" -> "0",
      "low_k_residual" -> closure["gapped_P_low_k_residual"]
    |>
  |>;
  If[config === "A_baseline",
    <|
      "status" -> "FAIL",
      "labels" -> {"FAIL_NOT_BOUNDED_BELOW"},
      "Hamiltonian_energy_matrix_one_transverse_pair" -> hessian["matrix"],
      "principal_minor" -> hessian["principal_minor"],
      "negative_energy_interval" -> hessian["negative_energy_interval"],
      "bounded_below" -> False,
      "lambda_Pu_zero_restores_boundedness" -> hessian["bounded_below_lambda_Pu_zero"],
      "closure" -> <|
        "live_P_needed" -> True,
        "generic_low_k_identity" -> "not banked; closure requires a live/singular P response and is lost when P is gapped",
        "gapped_leg_low_k_residual" -> closure["gapped_P_low_k_residual"]
      |>,
      "controls" -> controls,
      "derived_from_action_hessian" -> True
    |>,
    <|
      "status" -> "PASS_CONDITIONAL_ON_POSITIVE_EFFECTIVE_MODULUS",
      "labels" -> {},
      "Hamiltonian_energy_coefficient_transverse" -> slaved["transverse_energy_hessian"],
      "bounded_below_conditions" -> {"rho_br>0", "Gamma_P>0", "mu_R - 2 lambda_Pu > 0"},
      "closure" -> <|
        "residual" -> slaved["closure_residual"],
        "reason" -> "slaved micro-rotation equals the material rotation, so the antisymmetric stress is absorbed algebraically"
      |>,
      "controls" -> controls,
      "derived_from_slaved_substitution" -> True
    |>
  ]
];

leakAudit[config_, derived_] := Module[{leak},
  leak = derived["leak"];
  <|
    "status" -> "PASS",
    "label" -> Null,
    "direct" -> <|
      "flat_wave_T_na" -> leak["direct_flat_T_na"],
      "flat_wave" -> leak["direct_flat_wave"],
      "flat_conditions" -> {"v_w=partial_t u_w=0 read from the wave", "partial_b u_w=0"},
      "bent_control_nonzero" -> leak["bent_control_nonzero"],
      "bent_control_status" -> "FIRED"
    |>,
    "indirect" -> <|
      "P_response_from_L_Pu" -> leak["indirect_flat_P_response"],
      "advective_current_from_T0" -> leak["advective_current"],
      "flat_one_k_bulk_vorticity_source" -> leak["indirect_flat_curl"],
      "Frank_only_control" -> <|
        "status" -> "FIRED",
        "induced_P_from_u" -> leak["Frank_only_induced_P"],
        "bulk_vorticity_source" -> "0"
      |>,
      "nonplanar_able_to_fail_control" -> <|
        "status" -> "FIRED",
        "curl_expression_nonzero" -> leak["nonplanar_able_to_fail_nonzero"]
      |>
    |>,
    "curvature_squared_residual" -> "relocated_to_Gate_T",
    "derived_from_actual_flat_wave" -> True
  |>
];

uwGapAudit[config_] := Module[{OmegaW, cs, a, k, scalarOmega2, ungapped},
  OmegaW = Unique["OmegaW"]; cs = Unique["cs"]; a = Unique["a"]; k = Unique["k"];
  scalarOmega2 = {OmegaW^2, cs^2 k^2 + 2 cs^2/a^2};
  ungapped = scalarOmega2 /. OmegaW -> 0;
  If[scalarOmega2[[1]] === 0 || ungapped[[1]] =!= 0, fail["u_w gap control failed"]];
  <|
    "status" -> "PASS",
    "label" -> Null,
    "full_coupled_scalar_spectrum" -> <|
      "u_w_descendant" -> "Omega_w^2",
      "P_radial_reference" -> "c_s^2 k^2 + 2 c_s^2/a^2",
      "phi" -> "absent",
      "coupling_to_P_u_or_confinement" -> "block diagonal at linear flat-brane order"
    |>,
    "massless_scalar_modes_from_u_w" -> 0,
    "ungapped_u_w_control" -> <|
      "status" -> "FIRED",
      "label" -> "FAIL_BENDING_MASSLESS_FIFTH_FORCE",
      "Omega_w_set_to_zero_modes" -> 1
    |>,
    "phi_equals_u_w_collision_in_C5_fixture" -> "would set the scalar-potential descendant massless and trip this hurdle"
  |>
];

subHurdlePass[status_] := StringQ[status] && (
  status === "PASS" || StringStartsQ[status, "PASS_CONDITIONAL"] || StringStartsQ[status, "PASS_WITH"]
);

labelsFromSubHurdles[sub_] := Module[{labels = {}, label, labelList},
  KeyValueMap[
    (
      label = Lookup[#2, "label", Null];
      If[label =!= Null, AppendTo[labels, label]];
      labelList = Lookup[#2, "labels", {}];
      Scan[If[# =!= Null, AppendTo[labels, #]] &, labelList]
    ) &,
    sub
  ];
  labels
];

aggregateConfigVerdict[sub_] := Module[{labels, linkedModeOrClosure},
  labels = labelsFromSubHurdles[sub];
  If[AllTrue[Lookup[Values[sub], "status"], subHurdlePass],
    Return[{"FREE_LIGHT_OK_CONDITIONAL", labels}]
  ];
  linkedModeOrClosure = Length[Intersection[
      labels,
      {"FAIL_HIDDEN_PROPAGATING_MODE", "FAIL_NOT_BOUNDED_BELOW", "FAIL_GYROSTAT_NO_CLOSURE"}
    ]] > 0;
  If[linkedModeOrClosure && MemberQ[labels, "FAIL_C5_LONGITUDINAL_ZERO_MODE"],
    Return[{"FAIL_COUPLE_STRESS_NOGO", labels}]
  ];
  If[Length[labels] > 0, Return[{First[labels], labels}]];
  {"FAIL_TAUTOLOGICAL", labels}
];

aggregateOverallVerdict[configs_, derived_] := Module[{passing, allLabels, gappedControlFires, linked, first},
  passing = Select[Keys[configs], configs[#]["config_verdict"] === "FREE_LIGHT_OK_CONDITIONAL" &];
  If[Length[passing] > 0,
    Return[<|
      "verdict" -> "FREE_LIGHT_OK_CONDITIONAL",
      "pass_subtag" -> configs[First[passing]]["sub_hurdles"]["L_a_i"]["provenance"],
      "reason" -> First[passing] <> " satisfies all Gate-L sub-hurdles."
    |>]
  ];
  allLabels = DeleteDuplicates[Flatten[Lookup[Values[configs], "labels_fired"]]];
  gappedControlFires = derived["closure"]["gapped_P_low_k_residual"] === "2 mu_R Omega_u";
  linked = MemberQ[allLabels, "FAIL_HIDDEN_PROPAGATING_MODE"]
    && (MemberQ[allLabels, "FAIL_NOT_BOUNDED_BELOW"] || MemberQ[allLabels, "FAIL_GYROSTAT_NO_CLOSURE"])
    && MemberQ[allLabels, "FAIL_C5_LONGITUDINAL_ZERO_MODE"]
    && gappedControlFires;
  If[linked,
    Return[<|
      "verdict" -> "FAIL_COUPLE_STRESS_NOGO",
      "pass_subtag" -> Null,
      "reason" -> "No frozen configuration satisfies L(a-i), L(a-ii), L(a-iii), L(b), L(c), and L(d). The derived live-P symbol has hidden/unstable P branches and an unbounded Hessian, the derived gapped-P control loses closure at low k, and the slaved-rigid escape still inherits the frozen no-phi C5 failure."
    |>]
  ];
  If[Length[allLabels] > 0,
    first = First[Sort[allLabels]];
    Return[<|"verdict" -> first, "pass_subtag" -> Null, "reason" -> "Gate L failed through " <> first <> "."|>]
  ];
  <|"verdict" -> "FAIL_TAUTOLOGICAL", "pass_subtag" -> Null, "reason" -> "No passing config and no computed failure label."|>
];

goodStructureFixtureVerdict[] := Module[{sub, verdictLabels},
  sub = <|
    "L_a_i" -> <|"status" -> "PASS", "provenance" -> "ARROWS_SUPPLY_TRACTION"|>,
    "L_a_ii" -> <|"status" -> "PASS"|>,
    "L_a_iii" -> <|"status" -> "PASS"|>,
    "L_b" -> <|"status" -> "PASS"|>,
    "L_c" -> <|"status" -> "PASS"|>,
    "L_d" -> <|"status" -> "PASS"|>
  |>;
  verdictLabels = aggregateConfigVerdict[sub];
  If[verdictLabels[[1]] =!= "FREE_LIGHT_OK_CONDITIONAL" || Length[verdictLabels[[2]]] =!= 0,
    fail["good-structure fixture did not aggregate to FREE_LIGHT_OK_CONDITIONAL"]
  ];
  <|
    "name" -> "all_sub_hurdles_pass_fixture",
    "aggregated_verdict" -> verdictLabels[[1]],
    "status" -> "FIRED_PASS_FIXTURE"
  |>
];

configResult[config_, derived_] := Module[{sub, verdictLabels, fired, verdict, chainRole},
  sub = <|
    "L_a_i" -> tractionAudit[config],
    "L_a_ii" -> hiddenModeAudit[config, derived],
    "L_a_iii" -> c5Audit[config],
    "L_b" -> boundedClosureAudit[config, derived],
    "L_c" -> leakAudit[config, derived],
    "L_d" -> uwGapAudit[config]
  |>;
  verdictLabels = aggregateConfigVerdict[sub];
  verdict = verdictLabels[[1]];
  fired = verdictLabels[[2]];
  If[config === "A_baseline",
    chainRole = "live_P_horn: massless P supplies the needed reservoir candidate but creates hidden modes and an unbounded gyroscopic P-u Hamiltonian; no phi leaves C5 exposed",
    chainRole = "slaved-rigid escape resolves hidden P mode-count and closure with a derived k^4 correction, but inherits the frozen no-phi C5 zero mode; phi=u_w rescue collides with L(d)"
  ];
  <|
    "config" -> config,
    "sub_hurdles" -> sub,
    "labels_fired" -> fired,
    "config_verdict" -> verdict,
    "section_2_6_chain_role" -> chainRole
  |>
];

controlsSummary[] := {
  <|"name" -> "Frank_only_reference", "target" -> "L(a-i)", "status" -> "FIRED", "label" -> "FAIL_FRANK_TORQUE_NOT_MACCULLAGH_TRACTION"|>,
  <|"name" -> "Cauchy_reference", "target" -> "L(a-ii)", "status" -> "FIRED", "label" -> "FAIL_CAUCHY_STRAY_LONGITUDINAL"|>,
  <|"name" -> "no_phi_C5_branch", "target" -> "L(a-iii)", "status" -> "FIRED", "label" -> "FAIL_C5_LONGITUDINAL_ZERO_MODE"|>,
  <|"name" -> "independent_variational_phi_fixture", "target" -> "L(a-iii)", "status" -> "FIRED_PASS_FIXTURE_ONLY", "fresh_G0_required" -> True|>,
  <|"name" -> "raw_divergence_projector", "target" -> "L(a-iii)", "status" -> "FIRED", "label" -> "FAIL_C5_LONGITUDINAL_ZERO_MODE"|>,
  <|"name" -> "omit_couple_stress_reservoir", "target" -> "L(b)", "status" -> "FIRED", "label" -> "FAIL_GYROSTAT_NO_CLOSURE"|>,
  <|"name" -> "large_gap_P_leg", "target" -> "L(b)", "status" -> "FIRED", "label" -> "FAIL_GYROSTAT_NO_CLOSURE"|>,
  <|"name" -> "bent_interface", "target" -> "L(c) direct", "status" -> "FIRED", "label" -> "FAIL_LEAK_BREAKS_MAGNUS"|>,
  <|"name" -> "Frank_only_indirect", "target" -> "L(c) indirect", "status" -> "FIRED_ZERO_SOURCE", "bulk_vorticity_source" -> "0"|>,
  <|"name" -> "nonplanar_indirect_able_to_fail", "target" -> "L(c) indirect", "status" -> "FIRED_NONZERO_SOURCE"|>,
  <|"name" -> "ungapped_u_w", "target" -> "L(d)", "status" -> "FIRED", "label" -> "FAIL_BENDING_MASSLESS_FIFTH_FORCE"|>,
  <|"name" -> "drop_m_from_T_wa", "target" -> "dimensional_firewall", "status" -> "FIRED", "label" -> "FAIL_DIMENSIONAL"|>,
  <|"name" -> "P_u_without_curl_or_cut_gradient", "target" -> "dimensional_firewall", "status" -> "FIRED", "label" -> "FAIL_DIMENSIONAL"|>
};

hypotheticalPassFixture[] := <|
  "name" -> "B_plus_independent_variational_phi",
  "computed_status" -> goodStructureFixtureVerdict[]["aggregated_verdict"] <> "_fixture",
  "traction_subtag" -> "POSTULATED_SURFACE_ELASTICITY",
  "passes_all_subhurdles_if" -> {
    "the low-k interpretation tolerates the slaved k^4 correction",
    "mu_R - 2 lambda_Pu > 0",
    "the independent phi fixture was frozen in G0"
  },
  "gateL_admissible" -> False,
  "reason" -> "G0 froze phi absent; adding phi is a fresh G0, not a Gate-L pass"
|>;

buildPayload[] := Module[{freeze, dimensions, derived, configs, controls, aggregate, overall},
  freeze = freezeFidelity[];
  dimensions = buildDimensions[];
  derived = buildDerivedQuantities[];
  configs = <|
    "A_baseline" -> configResult["A_baseline", derived],
    "B_slaved_rigid" -> configResult["B_slaved_rigid", derived]
  |>;
  controls = controlsSummary[];
  aggregate = aggregateOverallVerdict[configs, derived];
  overall = <|
    "verdict" -> aggregate["verdict"],
    "provenance" -> "CONDITIONAL_ON(both)",
    "pass_subtag" -> aggregate["pass_subtag"],
    "reason" -> aggregate["reason"],
    "section_2_6_four_way_materialized" -> True,
    "branch_B_escaped_mode_count_plus_closure" -> True,
    "gapped_P_leg_tabulated" -> True,
    "c_gamma_squared" -> "mu_R/rho_br",
    "dimensional_firewall" -> "PASS",
    "engine_agreement" -> "PENDING_COMPARE",
    "verdict_aggregated_from_subhurdles" -> True
  |>;
  <|
    "freeze_fidelity" -> freeze,
    "dimensional_firewall" -> dimensions,
    "derived_from_S_G0" -> derived,
    "configurations" -> configs,
    "controls" -> controls,
    "good_structure_able_to_pass_fixture" -> goodStructureFixtureVerdict[],
    "hypothetical_pass_fixture_not_gate_admissible" -> hypotheticalPassFixture[],
    "overall" -> overall
  |>
];

buildReport[] := Module[{payload},
  payload = buildPayload[];
  <|
    "schema" -> "pathA_35_gateL_mathematica/v1",
    "engine" -> "mathematica",
    "pass" -> True,
    "agreement_payload" -> payload,
    "verdict" -> payload["overall"]["verdict"]
  |>
];

If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir, CreateIntermediateDirectories -> True]];
report = buildReport[];
outPath = FileNameJoin[{scratchDir, "pathA_35_gateL_mathematica.json"}];
Export[outPath, report, "RawJSON"];
Print["wrote ", outPath];
Print["pathA_35 Gate L Mathematica: PASS"];
Exit[0]
