(* Gate B4 Mathematica engine: baseline B0_Lifshitz smectic test. *)

ClearAll[
  fail, assertDim, assertZero, assertSame, extractBlock, scriptPath, stage1Root,
  stripConditional,
  reportsDir, scratchDir, t0Report, g0Report, sympyJson, jsonOut, expectedT0Hash,
  expectedG0Hash, t0Text, g0Text, t0Block, g0Block, t0Hash, g0Hash, requiredLines,
  Mdim, Ldim, Tdim, bulk, drho, dK, dhbar, dm, dk, dcL1, dcL2, da, dU2, dCQ,
  dA, domega2, dcs2, dB, dKbend, k, kk, rho, rhoField, rho0, K, hbar, m, cL1, cL2, a,
  U2, cs2, cQ, AkKk, Ak, omega2, kstar2, Amin, thresholdRoots, threshold,
  kstar2AtThreshold, thetaStiffness, pInertia, pTransStiffness, pLongStiffness,
  pTransOmega2, pLongOmega2, hessian, eta, theta, sigma, pi, pgrad, pdot, coeff,
  magDensity, frankDensity, inertiaDensity, thetaDensity, totalEnergyDensity,
  rhoPZeros, thetaPZeros, dot4, isZeroVec, cosModeItems, avgPowerCoeff,
  avgEtaPGrad2Coeff, minimizedPolynomial, shellDirections, lamellar, triad,
  lamAvg2, lamAvg3, lamAvg4, triAvg2, triAvg3, triAvg4, triEtaGrad2, U3,
  cubicCoeff, cubicCoeffKk, cubicZeroK2, r, u, Aamp, positiveControl,
  positiveControlDerived, dirs, avg2, avg4, energySH, fmin, amin, candidates,
  expected, tauVal, gammaVal, uVal, Flam, Ftri, lamNCMin, lamNCAmp, lamNCCands,
  triNCMin, triNCAmp, triNCCands, eps, beta, gammaExact, gammaAbs, FlamLandau,
  FtriLandau, FlamOnset, FtriOnset, lamLandauMin, lamLandauAmp, lamLandauCands,
  triLandauMin, triLandauAmp, triLandauCands, lamOnsetMin, lamOnsetAmp,
  lamOnsetCands, triOnsetMin, triOnsetAmp, triOnsetCands, triadBeatsLamella, openNeighborhood, q, U, USecond,
  gradientSymbol, qMin, gradientLowerBound, muBr, kappaPu, curlU, mismatch,
  deltaPParallel, OmegaU, lightIntegrand, lightStationary, lightMinimum,
  supportIntegrand, supportMinimum, sympyResults, sympyExprs, assertEngine,
  engineKeys, verdict, reason, conditions, results
];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

assertDim[name_, actual_, expected_] := If[FullSimplify[actual - expected] =!= {0, 0, 0},
  fail[name <> ": expected " <> ToString[expected, InputForm] <> ", got " <> ToString[actual, InputForm]]
];

assertZero[name_, expr_] := If[! TrueQ[FullSimplify[expr == 0]], fail[name <> " did not vanish: " <> ToString[expr, InputForm]]];

assertSame[name_, actual_, expected_] := If[! TrueQ[FullSimplify[actual == expected]],
  fail[name <> " mismatch: got " <> ToString[actual, InputForm] <> ", expected " <> ToString[expected, InputForm]]
];

stripConditional[expr_] := expr /. ConditionalExpression[val_, _] :> val;

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_25_gateB4_smectic.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
reportsDir = FileNameJoin[{stage1Root, "reports"}];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir]];
t0Report = FileNameJoin[{reportsDir, "pathA_24_T0_freeze.md"}];
g0Report = FileNameJoin[{reportsDir, "pathA_25_G0_freeze.md"}];
sympyJson = FileNameJoin[{scratchDir, "pathA_25_gateB4_smectic_sympy.json"}];
jsonOut = FileNameJoin[{scratchDir, "pathA_25_gateB4_smectic_mathematica.json"}];

expectedT0Hash = "8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064";
expectedG0Hash = "f00ee99d465e2e311c68f47fcacf4af0154ca650642271ab66c36d112cb6a290";

extractBlock[text_, path_] := Module[{blocks},
  blocks = StringCases[text, "```freeze-action\n" ~~ Shortest[block__] ~~ "\n```" :> block];
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
If[t0Hash =!= expectedT0Hash, fail["T0 hash mismatch: expected " <> expectedT0Hash <> ", got " <> t0Hash]];
If[g0Hash =!= expectedG0Hash, fail["G0 hash mismatch: expected " <> expectedG0Hash <> ", got " <> g0Hash]];
If[! StringContainsQ[g0Block, t0Block], fail["exact T0 freeze-action block is not embedded unchanged in G0 block"]];
requiredLines = {
  "Family L baseline B0_Lifshitz:",
  "lambda_Cdiv=0; chi_Cpin=0",
  "L_Mac =",
  "L_Pu =",
  "No new G0 term uses V_conf, a prescribed layer profile, a fixed layer normal"
};
If[! AllTrue[requiredLines, StringContainsQ[g0Block, #] &], fail["required frozen G0 line missing"]];

Mdim = {1, 0, 0}; Ldim = {0, 1, 0}; Tdim = {0, 0, 1};
bulk = Mdim - 2 Ldim - 2 Tdim;
drho = -4 Ldim;
dK = Mdim + 18 Ldim - 2 Tdim;
dhbar = Mdim + 2 Ldim - Tdim;
dm = Mdim;
dk = -Ldim;
dcL1 = Mdim + 8 Ldim - 2 Tdim;
dcL2 = Mdim + 10 Ldim - 2 Tdim;
da = Ldim;
dU2 = dK + 3 drho;
dCQ = 2 dhbar - dm - drho;
dA = bulk - 2 drho;
domega2 = -2 Tdim;
dcs2 = dK + 4 drho - dm;
dB = bulk;
dKbend = bulk + 2 Ldim;

assertDim["U''", dU2, dA];
assertDim["hbar^2/(m rho0)", dCQ, dcL1];
assertDim["c_L1", dcL1, dCQ];
assertDim["c_L2 k^2", dcL2 + 2 dk, dcL1];
assertDim["omega^2", drho - dm + 2 dk + dA, domega2];
assertDim["c_s^2", dcs2, 2 Ldim - 2 Tdim];
assertDim["P Goldstone omega^2", dcs2 + 2 dk, domega2];
assertDim["P magnitude gap omega^2", dcs2 - 2 da, domega2];

$Assumptions = k > 0 && rho0 > 0 && rho >= 0 && K > 0 && hbar > 0 &&
  m > 0 && cL1 > 0 && cL2 > 0 && a > 0 && r > 0 && u > 0;
U2 = 5 K rho0^3;
cs2 = 5 K rho0^4/m;
cQ = hbar^2/(4 m rho0);
AkKk = FullSimplify[U2 + (cQ - cL1) kk + cL2 kk^2];
kstar2 = FullSimplify[stripConditional[kk /. First[Block[{$Assumptions = True}, Solve[D[AkKk, kk] == 0, kk]]]]];
Ak = FullSimplify[AkKk /. kk -> k^2];
omega2 = FullSimplify[rho0 k^2 Ak/m];
Amin = FullSimplify[AkKk /. kk -> kstar2];
thresholdRoots = stripConditional /@ (cL1 /. Block[{$Assumptions = True}, Solve[Amin == 0, cL1]]);
threshold = Last[SortBy[thresholdRoots, N[# /. {K -> 1, rho0 -> 1, hbar -> 1, m -> 1, cL2 -> 1}, 30] &]];
kstar2AtThreshold = FullSimplify[kstar2 /. cL1 -> threshold];

thetaStiffness = FullSimplify[rho0 hbar^2 k^2/m];
pInertia = m rho0 a^2;
pTransStiffness = FullSimplify[m rho0 cs2 a^2 k^2];
pLongStiffness = FullSimplify[m rho0 cs2 a^2 k^2 + 2 m rho0 cs2];
pTransOmega2 = FullSimplify[pTransStiffness/pInertia];
pLongOmega2 = FullSimplify[pLongStiffness/pInertia];
hessian = DiagonalMatrix[{Ak, thetaStiffness, pTransStiffness, pTransStiffness, pTransStiffness, pLongStiffness}];
If[! AllTrue[Flatten[hessian - DiagonalMatrix[Diagonal[hessian]]], # === 0 &],
  fail["unexpected off-diagonal Hessian block"]
];

rhoField = rho0 + eta;
coeff = m rhoField (5 K rhoField^4/m);
magDensity = Expand[1/4 coeff (((1 + sigma)^2 + pi^2 - 1)^2)];
frankDensity = Expand[1/2 coeff a^2 pgrad^2];
inertiaDensity = Expand[1/2 m rho a^2 pdot^2];
thetaDensity = Expand[1/2 rho0 hbar^2 k^2/m theta^2];
totalEnergyDensity = magDensity + frankDensity + inertiaDensity + thetaDensity;
rhoPZeros = {
  D[magDensity, eta, sigma] /. {eta -> 0, sigma -> 0, pi -> 0},
  D[magDensity, eta, pi] /. {eta -> 0, sigma -> 0, pi -> 0},
  D[frankDensity, eta, pgrad] /. {eta -> 0, pgrad -> 0},
  D[inertiaDensity, eta, pdot] /. {eta -> 0, pdot -> 0}
};
thetaPZeros = {
  D[totalEnergyDensity, theta, sigma] /. {theta -> 0, sigma -> 0, pi -> 0, pgrad -> 0, pdot -> 0},
  D[totalEnergyDensity, theta, pi] /. {theta -> 0, sigma -> 0, pi -> 0, pgrad -> 0, pdot -> 0},
  D[totalEnergyDensity, theta, pgrad] /. {theta -> 0, sigma -> 0, pi -> 0, pgrad -> 0, pdot -> 0},
  D[totalEnergyDensity, theta, pdot] /. {theta -> 0, sigma -> 0, pi -> 0, pgrad -> 0, pdot -> 0}
};
assertZero["rho-P quadratic block", Total[rhoPZeros]];
assertZero["theta-P quadratic block", Total[thetaPZeros]];

dot4[uvec_, vvec_] := FullSimplify[uvec.vvec];
isZeroVec[v_] := And @@ (TrueQ[FullSimplify[# == 0]] & /@ v);
cosModeItems[base_] := Join[({#, 1/2} & /@ base), ({-#, 1/2} & /@ base)];
avgPowerCoeff[base_, power_] := Module[{items = cosModeItems[base], total = 0},
  Do[
    If[isZeroVec[Total[combo[[All, 1]]]], total += Times @@ combo[[All, 2]]],
    {combo, Tuples[items, power]}
  ];
  FullSimplify[total]
];
avgEtaPGrad2Coeff[base_, etaPower_] := Module[{items = cosModeItems[base], total = 0, etaVecs, etaCoeff},
  Do[
    etaVecs = etaCombo[[All, 1]];
    etaCoeff = Times @@ etaCombo[[All, 2]];
    Do[
      If[isZeroVec[Total[Join[etaVecs, {p[[1]], qv[[1]]}]]],
        total += etaCoeff p[[2]] qv[[2]] (-dot4[p[[1]], qv[[1]]])
      ],
      {p, items}, {qv, items}
    ],
    {etaCombo, Tuples[items, etaPower]}
  ];
  FullSimplify[total]
];

minimizedPolynomial[energy_, amp_, nonnegative_, rankRules_: {}] := Module[
  {roots, reals, values, order},
  roots = DeleteDuplicates[Join[amp /. Solve[D[energy, amp] == 0, amp], {0}]];
  reals = Select[roots,
    Abs[Im[N[# /. rankRules, 30]]] < 10^-20 &&
      (! nonnegative || Re[N[# /. rankRules, 30]] >= -10^-20) &
  ];
  If[Length[reals] == 0, fail["no real stationary candidates for " <> ToString[energy, InputForm]]];
  values = FullSimplify[energy /. amp -> #] & /@ reals;
  order = First[Ordering[N[values /. rankRules, 30], 1]];
  {FullSimplify[values[[order]]], FullSimplify[reals[[order]]], reals}
];

shellDirections[mm_] := Table[UnitVector[4, ii], {ii, 1, mm}];
lamellar = {{1, 0, 0, 0}};
triad = {{1, 0, 0, 0}, {-1/2, Sqrt[3]/2, 0, 0}, {-1/2, -Sqrt[3]/2, 0, 0}};
lamAvg2 = avgPowerCoeff[lamellar, 2];
lamAvg3 = avgPowerCoeff[lamellar, 3];
lamAvg4 = avgPowerCoeff[lamellar, 4];
triAvg2 = avgPowerCoeff[triad, 2];
triAvg3 = avgPowerCoeff[triad, 3];
triAvg4 = avgPowerCoeff[triad, 4];
triEtaGrad2 = avgEtaPGrad2Coeff[triad, 1];
assertSame["single-k lamellar cubic average", lamAvg3, 0];
assertSame["triad <phi^2>", triAvg2, 3/2];
assertSame["triad <phi^3>", triAvg3, 3/2];
assertSame["triad <phi |grad phi|^2>", triEtaGrad2, 3/4];

U3 = 15 K rho0^2;
cubicCoeff = FullSimplify[U3/6 triAvg3 - hbar^2/(8 m rho0^2) k^2 triEtaGrad2];
cubicCoeffKk = cubicCoeff /. k^2 -> kk;
cubicZeroK2 = FullSimplify[stripConditional[kk /. First[Block[{$Assumptions = True}, Solve[cubicCoeffKk == 0, kk]]]]];
assertSame["triad cubic coefficient", cubicCoeff, 15 K rho0^2/4 - 3 hbar^2 k^2/(32 m rho0^2)];
assertSame["triad cubic zero condition", cubicZeroK2, 40 K m rho0^4/hbar^2];

positiveControl = <||>;
positiveControlDerived = <||>;
Do[
  dirs = shellDirections[mm];
  avg2 = avgPowerCoeff[dirs, 2];
  avg4 = avgPowerCoeff[dirs, 4];
  energySH = FullSimplify[-1/2 r avg2 Aamp^2 + 1/4 u avg4 Aamp^4];
  {fmin, amin, candidates} = minimizedPolynomial[energySH, Aamp, True, {r -> 1, u -> 1}];
  expected = FullSimplify[-r^2 mm/(6 u (2 mm - 1))];
  assertSame["Swift-Hohenberg M=" <> ToString[mm] <> " minimized energy", fmin, expected];
  positiveControl[ToString[mm] <> "_directions"] = ToString[fmin, InputForm];
  positiveControlDerived[ToString[mm] <> "_directions"] = <|
    "avg_phi2" -> ToString[avg2, InputForm],
    "avg_phi4" -> ToString[avg4, InputForm],
    "energy" -> ToString[energySH, InputForm],
    "minimizing_amplitude" -> ToString[amin, InputForm],
    "stationary_candidates" -> (ToString[#, InputForm] & /@ candidates)
  |>,
  {mm, 1, 4}
];

tauVal = 1/100; gammaVal = 1; uVal = 1;
Flam = tauVal Aamp^2/4 + 3 uVal Aamp^4/32;
Ftri = 3 tauVal Aamp^2/4 - gammaVal Aamp^3/4 + 45 uVal Aamp^4/32;
{lamNCMin, lamNCAmp, lamNCCands} = minimizedPolynomial[Flam, Aamp, False];
{triNCMin, triNCAmp, triNCCands} = minimizedPolynomial[Ftri, Aamp, True];
If[! TrueQ[N[triNCMin, 30] < 0], fail["NC-B4d benchmark did not prefer the triad"]];

eps = 1/100;
beta = 1;
gammaExact = FullSimplify[cubicCoeff /. {K -> 1, rho0 -> 1, hbar -> 1, m -> 1, k^2 -> 23/8}];
If[TrueQ[gammaExact == 0], fail["generic baseline sample accidentally lies on Gamma_T=0"]];
gammaAbs = FullSimplify[Abs[gammaExact]];
FlamLandau = FullSimplify[-eps lamAvg2 Aamp^2 + beta lamAvg4 Aamp^4];
FtriLandau = FullSimplify[-eps triAvg2 Aamp^2 - gammaAbs Aamp^3 + beta triAvg4 Aamp^4];
FlamOnset = FullSimplify[beta lamAvg4 Aamp^4];
FtriOnset = FullSimplify[-gammaAbs Aamp^3 + beta triAvg4 Aamp^4];
{lamLandauMin, lamLandauAmp, lamLandauCands} = minimizedPolynomial[FlamLandau, Aamp, False];
{triLandauMin, triLandauAmp, triLandauCands} = minimizedPolynomial[FtriLandau, Aamp, True];
{lamOnsetMin, lamOnsetAmp, lamOnsetCands} = minimizedPolynomial[FlamOnset, Aamp, False];
{triOnsetMin, triOnsetAmp, triOnsetCands} = minimizedPolynomial[FtriOnset, Aamp, True];
triadBeatsLamella = TrueQ[N[triLandauMin - lamLandauMin, 30] < 0];
openNeighborhood = TrueQ[N[triOnsetMin - lamOnsetMin, 30] < 0];
If[! triadBeatsLamella || ! openNeighborhood, fail["computed baseline Landau comparison did not prefer the triad"]];

U = K rho^5/4;
USecond = FullSimplify[D[U, {rho, 2}]];
gradientSymbol = FullSimplify[1/2 cL2 q^2 - 1/2 cL1 q];
qMin = FullSimplify[q /. First[Solve[D[gradientSymbol, q] == 0, q]]];
gradientLowerBound = FullSimplify[gradientSymbol /. q -> qMin];
assertSame["bounded-below Fourier-symbol minimum", gradientLowerBound, -cL1^2/(8 cL2)];

lightIntegrand = 1/2 muBr curlU^2 + 1/2 kappaPu mismatch^2;
lightStationary = First[Solve[{D[lightIntegrand, curlU] == 0, D[lightIntegrand, mismatch] == 0}, {curlU, mismatch}]];
lightMinimum = FullSimplify[lightIntegrand /. lightStationary];
supportIntegrand = 1/2 muBr curlU^2 + 1/2 kappaPu (deltaPParallel - OmegaU)^2;
supportMinimum = FullSimplify[supportIntegrand /. {curlU -> 0, deltaPParallel -> OmegaU}];
assertSame["light-sector integrand minimum", lightMinimum, 0];
assertSame["light-sector support integrand at minimum", supportMinimum, 0];

If[! FileExistsQ[sympyJson], fail["missing SymPy JSON for engine agreement: " <> sympyJson]];
sympyResults = Import[sympyJson, "RawJSON"];
sympyExprs = sympyResults["engine_agreement"]["mathematica_exprs"];
assertEngine[name_, actual_] := Module[{expectedText, expectedExpr},
  expectedText = sympyExprs[name];
  If[! StringQ[expectedText], fail["missing SymPy engine expression for " <> name]];
  expectedExpr = ToExpression[expectedText, InputForm];
  If[! TrueQ[FullSimplify[actual == expectedExpr]],
    fail["engine disagreement " <> name <> ": Mathematica got " <> ToString[actual, InputForm] <>
      ", SymPy exported " <> expectedText]
  ];
];

engineKeys = {
  "density_kernel_A_of_k" -> Ak,
  "omega2_density_phase" -> omega2,
  "finite_kstar_squared" -> kstar2,
  "A_min" -> Amin,
  "softening_threshold_c_L1" -> threshold,
  "kstar_squared_at_threshold" -> kstar2AtThreshold,
  "P_transverse_omega2" -> pTransOmega2,
  "P_longitudinal_omega2" -> pLongOmega2,
  "rho_P_zero_sum" -> Total[rhoPZeros],
  "theta_P_zero_sum" -> Total[thetaPZeros],
  "single_k_lamellar_cubic_average" -> lamAvg3,
  "rank2_triad_average_phi2" -> triAvg2,
  "rank2_triad_average_phi3" -> triAvg3,
  "rank2_triad_average_phi_grad_phi_sq" -> triEtaGrad2,
  "baseline_rank2_triad_cubic_coefficient" -> cubicCoeff,
  "baseline_cubic_zero_kstar_squared" -> cubicZeroK2,
  "baseline_lamella_min" -> lamLandauMin,
  "baseline_triad_min" -> triLandauMin,
  "baseline_onset_triad_min" -> triOnsetMin,
  "negative_control_NC_B4d_lamellar_min" -> lamNCMin,
  "negative_control_NC_B4d_triad_min" -> triNCMin,
  "positive_control_SH_M1_min" -> ToExpression[positiveControl["1_directions"], InputForm],
  "positive_control_SH_M2_min" -> ToExpression[positiveControl["2_directions"], InputForm],
  "positive_control_SH_M3_min" -> ToExpression[positiveControl["3_directions"], InputForm],
  "positive_control_SH_M4_min" -> ToExpression[positiveControl["4_directions"], InputForm],
  "bounded_symbol_minimum" -> gradientLowerBound,
  "bounded_symbol_minimizing_k_squared" -> qMin,
  "phase_convexity_second_derivative" -> USecond,
  "k0_compressibility" -> U2,
  "light_integrand_minimum" -> lightMinimum,
  "light_support_integrand_at_minimum" -> supportMinimum
};
Scan[assertEngine[#[[1]], #[[2]]] &, engineKeys];

conditions = <|
  "finite_k_window_open" -> TrueQ[(threshold - cQ /. {K -> 1, rho0 -> 1, hbar -> 1, m -> 1, cL2 -> 1}) > 0],
  "phase_separation_excluded_by_convexity" -> TrueQ[FullSimplify[USecond >= 0, K > 0 && rho >= 0]],
  "bounded_below_symbol_minimum" -> TrueQ[FullSimplify[gradientLowerBound == -cL1^2/(8 cL2)]],
  "static_light_sector_nonnegative" -> TrueQ[lightMinimum == 0 && supportMinimum == 0],
  "Gamma_T_nonzero_on_generic_sample" -> TrueQ[gammaExact != 0],
  "triad_minimum_below_lamella_minimum" -> triadBeatsLamella,
  "open_neighborhood_certificate" -> openNeighborhood
|>;

If[! conditions["finite_k_window_open"],
  verdict = "FAIL_NO_MODULATION"; reason = "computed finite-k softening window is empty",
  If[! conditions["phase_separation_excluded_by_convexity"],
    verdict = "FAIL_PHASE_SEPARATION"; reason = "computed local density convexity failed at fixed mean",
    If[! conditions["bounded_below_symbol_minimum"] || ! conditions["static_light_sector_nonnegative"],
      verdict = "FAIL_STABILITY"; reason = "computed boundedness or static light-sector positivity failed",
      If[conditions["Gamma_T_nonzero_on_generic_sample"] &&
          conditions["triad_minimum_below_lamella_minimum"] &&
          conditions["open_neighborhood_certificate"],
        verdict = "FAIL_NOT_CODIM1";
        reason = "computed rank-2 triad Landau minimum is below the single-k lamella in an open neighborhood above onset",
        verdict = "SMECTIC_CONDITIONAL";
        reason = "computed checks did not find a lower multi-k state off the codimension-one cubic-zero surface"
      ]
    ]
  ]
];
If[verdict =!= sympyResults["verdict"], fail["verdict disagreement with SymPy JSON"]];

results = <|
  "verdict" -> verdict,
  "reason" -> reason,
  "freeze_fidelity" -> <|
    "t0_hash" -> t0Hash,
    "g0_hash" -> g0Hash,
    "t0_bytes_embedded_in_g0" -> True,
    "required_frozen_lines_present" -> True
  |>,
  "computed_verdict" -> <|
    "conditions" -> conditions,
    "computed_baseline_comparison" -> <|
      "lamellar_min_exact" -> ToString[lamLandauMin, InputForm],
      "lamellar_min_numeric" -> N[lamLandauMin, 20],
      "triad_min_exact" -> ToString[triLandauMin, InputForm],
      "triad_min_numeric" -> N[triLandauMin, 20]
    |>
  |>,
  "coupled_spectrum" -> <|
    "density_kernel_A_of_k" -> ToString[Ak, InputForm],
    "omega2_density_phase" -> ToString[omega2, InputForm],
    "finite_kstar_squared" -> ToString[kstar2, InputForm],
    "A_min" -> ToString[Amin, InputForm],
    "softening_threshold_c_L1" -> ToString[threshold, InputForm],
    "kstar_squared_at_threshold" -> ToString[kstar2AtThreshold, InputForm],
    "P_transverse_omega2" -> ToString[pTransOmega2, InputForm],
    "P_longitudinal_omega2" -> ToString[pLongOmega2, InputForm],
    "rho_P_quadratic_block" -> "computed_zero",
    "theta_P_quadratic_block" -> "computed_zero"
  |>,
  "morphology" -> <|
    "single_k_lamellar_cubic_average" -> ToString[lamAvg3, InputForm],
    "rank2_triad_average_phi2" -> ToString[triAvg2, InputForm],
    "rank2_triad_average_phi3" -> ToString[triAvg3, InputForm],
    "rank2_triad_average_phi_grad_phi_sq" -> ToString[triEtaGrad2, InputForm],
    "baseline_rank2_triad_cubic_coefficient" -> ToString[cubicCoeff, InputForm],
    "baseline_cubic_zero_kstar_squared" -> ToString[cubicZeroK2, InputForm],
    "positive_control_NC_B4b" -> positiveControl,
    "negative_control_NC_B4d" -> <|
      "lamellar_min" -> ToString[lamNCMin, InputForm],
      "triad_min_exact" -> ToString[triNCMin, InputForm],
      "triad_min_numeric" -> N[triNCMin, 20]
    |>
  |>,
  "finite_k_vs_k0" -> <|
    "bounded_symbol_minimum" -> ToString[gradientLowerBound, InputForm],
    "bounded_symbol_minimizing_k_squared" -> ToString[qMin, InputForm],
    "phase_convexity_second_derivative" -> ToString[USecond, InputForm]
  |>,
  "light_sector" -> <|
    "integrand_minimum" -> ToString[lightMinimum, InputForm],
    "support_integrand_at_minimum" -> ToString[supportMinimum, InputForm],
    "uniform_state_support" -> "analytic_argument_not_computed"
  |>,
  "engine_agreement" -> <|
    "sympy_json" -> sympyJson,
    "shared_expression_count" -> Length[engineKeys],
    "status" -> "PASS"
  |>
|>;

Export[jsonOut, results, "JSON"];
Print["PASS pathA_25_gateB4_smectic_mathematica"];
Print[ExportString[<|
  "verdict" -> verdict,
  "engine_agreement" -> "PASS",
  "baseline_lamella_min" -> N[lamLandauMin, 20],
  "baseline_triad_min" -> N[triLandauMin, 20],
  "json" -> jsonOut
|>, "JSON"]];
Exit[0];
