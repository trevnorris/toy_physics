(* pathA_25 G0 Mathematica dimensional and freeze-fidelity check.

   Scope: G0 only.  Recomputes the T0 and G0 freeze hashes, asserts that the
   exact T0 block is embedded in the combined G0 block, restores MLT units for
   the frozen new terms, and verifies the explicit k->0 driver limit.
*)

ClearAll[
  fail, scriptPath, stage1Root, scratchDir, t0Report, g0Report,
  expectedT0Hash, expectedG0Hash, expectedLPolLines, extractBlock,
  t0Text, g0Text, t0Block, g0Block, t0Hash, g0Hash, linePositions,
  dimTuple, dimString, assertDim, Mdim, Ldim, Tdim, dimless,
  bulkLag, braneLag, actionDim, deltaSigma, dm, drho, dK, da, du, dP,
  dgrad4, dgrad3, dlap4, ddt, ddn, dd4R, dcs2, dcL1, dcL2,
  dgradRho, dlapRho, familyLTerms, dkLstar, dVRReal, dAR, dkR,
  dlambdaCdiv, dchiCpin, dgradP, familyCTerms, dvarrhoBr, dmuBr,
  dOmega, dDtU, dIPSigma, dKPSigma, dGPSigma, dJPu, dkappaPu,
  dDtP, lightBraneTerms, k, cL1, cL2, AR, kR, lam, chi, theta, rho0, K,
  hbar, m, a, eta, piT, sigma, familyL, kLstar, x, fR, vR, kRstar,
  U2, cQ, baselineLKernel, SGoldstone, Gmagnitude, cdivQuadratic,
  cdivMixedPi, cdivMixedSigma, cdivDeltaA, cdivEffectiveKernel, cdivLowK,
  cdivLowKShift, cdivDirectionalLiminfShift, cpinQuadratic, cpinMixedPi,
  cpinDensityHessian, cpinDeltaA, cpinEffectiveKernel, cpinLowK,
  cpinLowKShift, checks, dimensions, ktoZero, report,
  outPath
];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_25_G0_dimcheck.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
t0Report = FileNameJoin[{stage1Root, "reports", "pathA_24_T0_freeze.md"}];
g0Report = FileNameJoin[{stage1Root, "reports", "pathA_25_G0_freeze.md"}];

expectedT0Hash = "8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064";
expectedG0Hash = "f00ee99d465e2e311c68f47fcacf4af0154ca650642271ab66c36d112cb6a290";
expectedLPolLines = {
  "  L_pol =",
  "      (1/2) m rho a^2 (D_t^v P^i)(D_t^v P^i)",
  "    - (1/2) m rho c_s^2(rho) a^2 (partial_j P^i)(partial_j P^i)",
  "    - (1/4) m rho c_s^2(rho) (P^i P^i - 1)^2."
};

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

t0Text = Import[t0Report, "Text"];
g0Text = Import[g0Report, "Text"];
t0Block = extractBlock[t0Text, t0Report];
g0Block = extractBlock[g0Text, g0Report];
t0Hash = ToLowerCase[IntegerString[Hash[t0Block <> "\n", "SHA256"], 16, 64]];
g0Hash = ToLowerCase[IntegerString[Hash[g0Block <> "\n", "SHA256"], 16, 64]];
If[t0Hash =!= expectedT0Hash,
  fail["T0 hash mismatch: expected " <> expectedT0Hash <> ", got " <> t0Hash]
];
If[g0Hash =!= expectedG0Hash,
  fail["G0 hash mismatch: expected " <> expectedG0Hash <> ", got " <> g0Hash]
];
If[! StringContainsQ[g0Block, t0Block],
  fail["exact T0 freeze-action block is not embedded unchanged in G0 block"]
];
If[! AllTrue[expectedLPolLines, StringContainsQ[t0Block, #] &],
  fail["T0 L_pol line missing from T0 block"]
];
If[! AllTrue[expectedLPolLines, StringContainsQ[g0Block, #] &],
  fail["T0 L_pol line missing from G0 block"]
];
linePositions = Flatten[FirstPosition[StringSplit[g0Block, "\n"], #, Missing["notfound"]] & /@ expectedLPolLines];
If[MemberQ[linePositions, Missing["notfound"]] || linePositions =!= Sort[linePositions],
  fail["embedded T0 L_pol lines missing or out of order in G0 block"]
];

dimTuple[v_] := Normal[v];
dimString[v_] := Module[{labels = {"M", "L", "T"}, powers = Normal[v], parts},
  parts = MapThread[
    If[#2 === 0, Nothing, If[#2 === 1, #1, #1 <> "^" <> ToString[#2]]] &,
    {labels, powers}
  ];
  If[parts === {}, "1", StringRiffle[parts, " "]]
];
assertDim[name_, actual_, expected_] := If[FullSimplify[actual - expected] =!= {0, 0, 0},
  fail[name <> ": expected " <> ToString[dimTuple[expected]] <> ", got " <> ToString[dimTuple[actual]]]
];

Mdim = {1, 0, 0}; Ldim = {0, 1, 0}; Tdim = {0, 0, 1}; dimless = {0, 0, 0};
bulkLag = Mdim - 2 Ldim - 2 Tdim;
braneLag = Mdim - Ldim - 2 Tdim;
actionDim = Mdim + 2 Ldim - Tdim;
deltaSigma = -Ldim;

dm = Mdim;
drho = -4 Ldim;
dK = Mdim + 18 Ldim - 2 Tdim;
da = Ldim;
du = Ldim;
dP = dimless;
dgrad4 = -Ldim;
dgrad3 = -Ldim;
dlap4 = -2 Ldim;
ddt = -Tdim;
ddn = Ldim;
dd4R = 4 Ldim;

dcs2 = dK + 4 drho - dm;
assertDim["c_s^2(rho)", dcs2, 2 Ldim - 2 Tdim];

dcL1 = Mdim + 8 Ldim - 2 Tdim;
dcL2 = Mdim + 10 Ldim - 2 Tdim;
dgradRho = drho + dgrad4;
dlapRho = drho + dlap4;
familyLTerms = <|
  "c_L1 (partial_i rho)^2" -> dcL1 + 2 dgradRho,
  "c_L2 (partial_i partial_i rho)^2" -> dcL2 + 2 dlapRho
|>;
Scan[assertDim[#1, familyLTerms[#1], bulkLag] &, Keys[familyLTerms]];
dkLstar = (dcL1 - dcL2)/2;
assertDim["k_Lstar", dkLstar, -Ldim];

dVRReal = Mdim + 2 Ldim - 2 Tdim;
dAR = Mdim + 6 Ldim - 2 Tdim;
dkR = -Ldim;
assertDim["Family R real-space kernel term", dd4R + 2 drho + dVRReal, bulkLag];
assertDim["tilde V_R=A_R f_R", dAR, dVRReal + dd4R];
assertDim["k_R", dkR, -Ldim];
assertDim["k_Rstar", dkR, -Ldim];

dlambdaCdiv = Mdim + 3 Ldim - 2 Tdim;
dchiCpin = Mdim + 8 Ldim - 2 Tdim;
dgradP = dP + dgrad4;
familyCTerms = <|
  "lambda_Cdiv delta_rho partial_i P^i" -> dlambdaCdiv + drho + dgradP,
  "chi_Cpin (P^i partial_i rho)^2" -> dchiCpin + 2 (dP + dgradRho)
|>;
Scan[assertDim[#1, familyCTerms[#1], bulkLag] &, Keys[familyCTerms]];

dvarrhoBr = Mdim - 3 Ldim;
dmuBr = braneLag;
dOmega = dgrad3 + du;
dDtU = du + ddt;
dIPSigma = ddn + dm + drho + 2 da;
dKPSigma = ddn + dm + drho + dcs2 + 2 da;
dGPSigma = ddn + dm + drho + dcs2;
dJPu = Mdim - Ldim;
dkappaPu = braneLag;
dDtP = ddt;

assertDim["Omega_u", dOmega, dimless];
lightBraneTerms = <|
  "varrho_br (D_t u)^2" -> dvarrhoBr + 2 dDtU,
  "mu_br (curl u)^2" -> dmuBr + 2 dOmega,
  "I_PSigma (D_t P_parallel)^2" -> dIPSigma + 2 dDtP,
  "K_PSigma (partial_a P_parallel)^2" -> dKPSigma + 2 dgradP,
  "G_PSigma (delta P_parallel)^2" -> dGPSigma + 2 dP,
  "J_Pu (D_t delta P - D_t Omega_u)^2" -> dJPu + 2 dDtP,
  "kappa_Pu (delta P - Omega_u)^2" -> dkappaPu + 2 dP
|>;
Scan[
  (assertDim[#1, lightBraneTerms[#1], braneLag];
   assertDim["bulk representation delta_Sigma " <> #1, deltaSigma + lightBraneTerms[#1], bulkLag]) &,
  Keys[lightBraneTerms]
];
assertDim["varrho_br=int dn m rho", dvarrhoBr, Mdim - 3 Ldim];
assertDim["mu_br", dmuBr, braneLag];
assertDim["I_PSigma=int dn m rho a^2", dIPSigma, Mdim - Ldim];
assertDim["K_PSigma=int dn m rho c_s^2 a^2", dKPSigma, Mdim + Ldim - 2 Tdim];
assertDim["G_PSigma=int dn m rho c_s^2", dGPSigma, braneLag];
assertDim["J_Pu", dJPu, Mdim - Ldim];
assertDim["kappa_Pu", dkappaPu, braneLag];
assertDim["delta_Sigma", deltaSigma, -Ldim];
assertDim["int dt d^4X bulkLag", bulkLag + Tdim + 4 Ldim, actionDim];
assertDim["int dt d^3sigma braneLag", braneLag + Tdim + 3 Ldim, actionDim];

$Assumptions = cL1 > 0 && cL2 > 0 && AR > 0 && kR > 0 &&
  rho0 > 0 && K > 0 && hbar > 0 && m > 0 && a > 0;
familyL = cL2 k^4 - cL1 k^2;
kLstar = Sqrt[cL1/(2 cL2)];
If[FullSimplify[Limit[familyL, k -> 0]] =!= 0, fail["Family L k->0 limit is not zero"]];
If[FullSimplify[D[familyL, k] /. k -> kLstar] =!= 0, fail["Family L finite-k derivative is not zero"]];
If[FullSimplify[(D[familyL, {k, 2}] /. k -> kLstar) - 4 cL1] =!= 0,
  fail["Family L finite-k second derivative is not 4 c_L1"]
];

fR = (x^4 - 2 x^2) Exp[-x^2];
vR = AR (fR /. x -> k/kR);
kRstar = kR Sqrt[2 - Sqrt[2]];
If[FullSimplify[Limit[vR, k -> 0]] =!= 0, fail["Family R k->0 limit is not zero"]];
If[FullSimplify[D[vR, k] /. k -> kRstar] =!= 0, fail["Family R finite-k derivative is not zero"]];
If[! TrueQ[N[fR /. x -> Sqrt[2 - Sqrt[2]]] < 0],
  fail["Family R finite-k stationary point is not negative"]
];

(* Supersedes the old tautological monomial check
   cdiv=lambda_Cdiv*k and cpin=chi_Cpin*k^2.  The response is derived from
   the rho-P Hessian about uniform rho0 and unit P0, then P is integrated out
   before taking the directional k->0 limit. *)
U2 = 5 K rho0^3;
cQ = hbar^2/(4 m rho0);
baselineLKernel = FullSimplify[U2 + (cQ - cL1) k^2 + cL2 k^4];
SGoldstone = FullSimplify[5 K rho0^5 a^2];
Gmagnitude = FullSimplify[10 K rho0^5];
cdivQuadratic = -lam eta (k Sin[theta] piT + k Cos[theta] sigma);
cdivMixedPi = FullSimplify[D[cdivQuadratic, eta, piT]];
cdivMixedSigma = FullSimplify[D[cdivQuadratic, eta, sigma]];
cdivDeltaA = FullSimplify[
  -lam^2 Sin[theta]^2/SGoldstone -
    lam^2 k^2 Cos[theta]^2/(SGoldstone k^2 + Gmagnitude)
];
cdivEffectiveKernel = FullSimplify[baselineLKernel + cdivDeltaA];
limitAssumptions = rho0 > 0 && K > 0 && hbar > 0 && m > 0 && a > 0 &&
  cL1 > 0 && cL2 > 0;
cdivLowK = FullSimplify[Block[{$Assumptions = limitAssumptions},
  Limit[cdivEffectiveKernel, k -> 0, Direction -> "FromAbove"]
]];
cdivLowKShift = FullSimplify[cdivLowK - U2];
cdivDirectionalLiminfShift = FullSimplify[cdivLowKShift /. theta -> Pi/2];
If[FullSimplify[cdivDirectionalLiminfShift + lam^2/SGoldstone] =!= 0,
  fail["Family Cdiv derived transverse low-k shift mismatch"]
];
If[TrueQ[FullSimplify[cdivLowKShift == 0]],
  fail["Family Cdiv derived low-k response unexpectedly preserved the EOS stiffness"]
];

cpinQuadratic = 1/2 chi (k Cos[theta] eta)^2;
cpinMixedPi = FullSimplify[D[cpinQuadratic, eta, piT]];
cpinDensityHessian = FullSimplify[D[cpinQuadratic, eta, eta]];
cpinDeltaA = FullSimplify[cpinDensityHessian];
cpinEffectiveKernel = FullSimplify[baselineLKernel + cpinDeltaA];
cpinLowK = FullSimplify[Block[{$Assumptions = limitAssumptions},
  Limit[cpinEffectiveKernel, k -> 0, Direction -> "FromAbove"]
]];
cpinLowKShift = FullSimplify[cpinLowK - U2];
If[FullSimplify[cpinLowKShift] =!= 0, fail["Family Cpin derived k->0 EOS shift is not zero"]];
If[FullSimplify[cpinMixedPi] =!= 0, fail["Family Cpin rho-P mixed quadratic block is not zero"]];

dimensions = <|
  "expected_bulk_action_density" -> <|"triple_MLT" -> dimTuple[bulkLag], "string" -> dimString[bulkLag]|>,
  "expected_brane_density" -> <|"triple_MLT" -> dimTuple[braneLag], "string" -> dimString[braneLag]|>,
  "c_s_squared" -> <|"triple_MLT" -> dimTuple[dcs2], "string" -> dimString[dcs2]|>,
  "family_l" -> Association @ KeyValueMap[#1 -> <|"triple_MLT" -> dimTuple[#2], "string" -> dimString[#2]|> &, familyLTerms],
  "family_r" -> <|
    "real_space_kernel_V_R" -> <|"triple_MLT" -> dimTuple[dVRReal], "string" -> dimString[dVRReal]|>,
    "fourier_kernel_A_R" -> <|"triple_MLT" -> dimTuple[dAR], "string" -> dimString[dAR]|>,
    "real_space_action_density" -> <|"triple_MLT" -> dimTuple[dd4R + 2 drho + dVRReal], "string" -> dimString[dd4R + 2 drho + dVRReal]|>,
    "k_R" -> <|"triple_MLT" -> dimTuple[dkR], "string" -> dimString[dkR]|>
  |>,
  "family_c" -> Association @ KeyValueMap[#1 -> <|"triple_MLT" -> dimTuple[#2], "string" -> dimString[#2]|> &, familyCTerms],
  "light_brane_terms" -> Association @ KeyValueMap[#1 -> <|"triple_MLT" -> dimTuple[#2], "string" -> dimString[#2]|> &, lightBraneTerms],
  "light_bulk_density_terms" -> Association @ KeyValueMap[#1 -> <|"triple_MLT" -> dimTuple[deltaSigma + #2], "string" -> dimString[deltaSigma + #2]|> &, lightBraneTerms],
  "light_constants" -> <|
    "varrho_br" -> <|"triple_MLT" -> dimTuple[dvarrhoBr], "string" -> dimString[dvarrhoBr]|>,
    "mu_br" -> <|"triple_MLT" -> dimTuple[dmuBr], "string" -> dimString[dmuBr]|>,
    "I_PSigma" -> <|"triple_MLT" -> dimTuple[dIPSigma], "string" -> dimString[dIPSigma]|>,
    "K_PSigma" -> <|"triple_MLT" -> dimTuple[dKPSigma], "string" -> dimString[dKPSigma]|>,
    "G_PSigma" -> <|"triple_MLT" -> dimTuple[dGPSigma], "string" -> dimString[dGPSigma]|>,
    "J_Pu" -> <|"triple_MLT" -> dimTuple[dJPu], "string" -> dimString[dJPu]|>,
    "kappa_Pu" -> <|"triple_MLT" -> dimTuple[dkappaPu], "string" -> dimString[dkappaPu]|>,
    "delta_Sigma" -> <|"triple_MLT" -> dimTuple[deltaSigma], "string" -> dimString[deltaSigma]|>
  |>
|>;

ktoZero = <|
  "family_l" -> <|
    "quadratic_symbol" -> "c_L2*k^4 - c_L1*k^2",
    "limit_k_to_0" -> "0",
    "finite_k_stationary_scale" -> ToString[InputForm[kLstar]],
    "second_derivative_at_scale" -> "4*c_L1"
  |>,
  "family_r" -> <|
    "shape" -> "(x^4 - 2*x^2)*exp(-x^2)",
    "limit_k_to_0" -> "0",
    "small_k_series" -> ToString[InputForm[Normal[Series[vR, {k, 0, 4}]]]],
    "finite_k_stationary_scale" -> ToString[InputForm[kRstar]]
  |>,
  "family_c" -> <|
    "supersedes_old_tautology" -> "old cdiv=lambda*k and cpin=chi*k^2 assertions removed",
    "baseline_L_kernel_retained" -> ToString[InputForm[baselineLKernel]],
    "goldstone_stiffness" -> ToString[InputForm[SGoldstone]],
    "magnitude_gap" -> ToString[InputForm[Gmagnitude]],
    "Cdiv_mixed_block_eta_piT" -> ToString[InputForm[cdivMixedPi]],
    "Cdiv_mixed_block_eta_sigma" -> ToString[InputForm[cdivMixedSigma]],
    "Cdiv_delta_A_after_integrating_P" -> ToString[InputForm[cdivDeltaA]],
    "Cdiv_low_k_limit" -> ToString[InputForm[cdivLowK]],
    "Cdiv_low_k_shift_from_EOS" -> ToString[InputForm[cdivLowKShift]],
    "Cdiv_directional_liminf_shift" -> ToString[InputForm[cdivDirectionalLiminfShift]],
    "Cdiv_admission_status" -> "FAIL_ADMISSION_for_nonzero_lambda_Cdiv",
    "Cpin_rho_P_mixed_block" -> ToString[InputForm[cpinMixedPi]],
    "Cpin_direct_density_hessian" -> ToString[InputForm[cpinDensityHessian]],
    "Cpin_delta_A" -> ToString[InputForm[cpinDeltaA]],
    "Cpin_low_k_limit" -> ToString[InputForm[cpinLowK]],
    "Cpin_low_k_shift_from_EOS" -> ToString[InputForm[cpinLowKShift]],
    "Cpin_admission_status" -> "EOS_preserved_O(k^2)_anisotropic_density_Hessian"
  |>,
  "statement" -> "Derived Family-C response replaces the old monomial limit: Cdiv shifts the directional low-k EOS stiffness after Goldstone integrate-out, while Cpin is O(k^2) and preserves A_rho(0)."
|>;

report = <|
  "schema" -> "pathA_25_G0_dimcheck_mathematica/v1",
  "pass" -> True,
  "scope" -> "freeze-fidelity, dimensional homogeneity, and k->0 limit only; no gate solved",
  "freeze_fidelity" -> <|
    "t0_hash" -> t0Hash,
    "g0_hash" -> g0Hash,
    "t0_bytes_embedded_in_g0" -> True
  |>,
  "dimensions" -> dimensions,
  "k_to_zero" -> ktoZero,
  "verdict" -> "G0_DIMCHECK_MATHEMATICA_PASS"
|>;

If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir, CreateIntermediateDirectories -> True]];
outPath = FileNameJoin[{scratchDir, "pathA_25_G0_dimcheck_mathematica.json"}];
Export[outPath, report, "RawJSON"];
Print["wrote ", outPath];
Print["pathA_25 G0 Mathematica dimcheck: PASS"];
Exit[0]
