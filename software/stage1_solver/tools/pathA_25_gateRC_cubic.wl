(* Gate R/C Mathematica engine: Family R/C cubic and anisotropic-shell test. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

assertSame[name_, actual_, expected_] := If[! TrueQ[FullSimplify[actual == expected]],
  fail[name <> " mismatch: got " <> ToString[actual, InputForm] <> ", expected " <> ToString[expected, InputForm]]
];

stripConditional[expr_] := expr /. ConditionalExpression[val_, _] :> val;

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_25_gateRC_cubic.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
reportsDir = FileNameJoin[{stage1Root, "reports"}];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
t0Report = FileNameJoin[{reportsDir, "pathA_24_T0_freeze.md"}];
g0Report = FileNameJoin[{reportsDir, "pathA_25_G0_freeze.md"}];
sympyJson = FileNameJoin[{scratchDir, "pathA_25_gateRC_cubic_sympy.json"}];
jsonOut = FileNameJoin[{scratchDir, "pathA_25_gateRC_cubic_mathematica.json"}];

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
  "Family R sensitivity S_R_kernel:",
  "Family C sensitivities:",
  "lambda_Cdiv",
  "chi_Cpin",
  "L_Pu =",
  "No new G0 term uses V_conf, a prescribed layer profile, a fixed layer normal"
};
If[! AllTrue[requiredLines, StringContainsQ[g0Block, #] &], fail["required frozen G0 line missing"]];

Mdim = {1, 0, 0}; Ldim = {0, 1, 0}; Tdim = {0, 0, 1};
bulk = Mdim - 2 Ldim - 2 Tdim;
drho = -4 Ldim; dK = Mdim + 18 Ldim - 2 Tdim; dhbar = Mdim + 2 Ldim - Tdim;
dm = Mdim; da = Ldim; dk = -Ldim;
dA = bulk - 2 drho; dlambda = Mdim + 3 Ldim - 2 Tdim; dchi = Mdim + 8 Ldim - 2 Tdim;
dSG = Mdim - 2 Tdim; dgap = bulk; dgamma = Mdim + 10 Ldim - 2 Tdim;
dmuBr = Mdim - Ldim - 2 Tdim; dkappaPu = Mdim - Ldim - 2 Tdim;
dKPSigma = Mdim + Ldim - 2 Tdim; dqParallel = -Ldim; dgRhoSqSigma = -9 Ldim;
assertDim[name_, actual_, expected_] := If[FullSimplify[actual - expected] =!= {0, 0, 0},
  fail[name <> ": expected " <> ToString[expected, InputForm] <> ", got " <> ToString[actual, InputForm]]
];
assertDim["A_rho", dA, Mdim + 6 Ldim - 2 Tdim];
assertDim["Delta A_Cdiv Goldstone", 2 dlambda - dSG, dA];
assertDim["Delta A_Cdiv magnitude", 2 dlambda + 2 dk - dgap, dA];
assertDim["Delta A_Cpin", dchi + 2 dk, dA];
assertDim["Gamma cubic coefficient", dgamma, Mdim + 10 Ldim - 2 Tdim];
assertDim["Goldstone stiffness S_G", dSG, Mdim - 2 Tdim];
assertDim["magnitude gap", dgap, bulk];
assertDim["hbar^2 k^2/(m rho0)", 2 dhbar + 2 dk - dm - drho, dA];
assertDim["K rho0^3", dK + 3 drho, dA];
assertDim["K rho0^5 a^2", dK + 5 drho + 2 da, dSG];
assertDim["MacCullagh mu_br", dmuBr, Mdim - Ldim - 2 Tdim];
assertDim["P-u anchoring kappa_Pu", dkappaPu, dmuBr];
assertDim["surface Frank K_PSigma", dKPSigma, Mdim + Ldim - 2 Tdim];
assertDim["K_PSigma q_parallel^2", dKPSigma + 2 dqParallel, dkappaPu];
assertDim["chi_Cpin <|grad rho|^2>_Sigma", dchi + dgRhoSqSigma, dkappaPu];

$Assumptions = k > 0 && x > 0 && rho0 > 0 && K > 0 && hbar > 0 && m > 0 &&
  a > 0 && cL1 > 0 && cL2 > 0 && AR > 0 && kR > 0 && lam > 0 &&
  eps > 0 && beta > 0 && gamma > 0 && zeta > 0 && muBr > 0 && kappaPu > 0 &&
  chiAbs > 0 && gRhoSqSigma > 0 && KPSigma > 0 && qParallel > 0 &&
  avg2O > 0 && c3O > 0 && u4O > 0;

U2 = 5 K rho0^3;
cQ = hbar^2/(4 m rho0);
AL = FullSimplify[U2 + (cQ - cL1) k^2 + cL2 k^4];
SG = FullSimplify[5 K rho0^5 a^2];
Gmag = FullSimplify[10 K rho0^5];
kperp2 = k^2 Sin[theta]^2;
kpar2 = k^2 Cos[theta]^2;
cdivTransverse = FullSimplify[-lam^2 kperp2/(SG k^2)];
cdivMagnitude = FullSimplify[-lam^2 kpar2/(SG k^2 + Gmag)];
cdivDelta = FullSimplify[cdivTransverse + cdivMagnitude];
ACdiv = FullSimplify[AL + cdivDelta];
omega2Cdiv = FullSimplify[rho0 k^2 ACdiv/m];
limitAssumptions = rho0 > 0 && K > 0 && hbar > 0 && m > 0 && a > 0 &&
  cL1 > 0 && cL2 > 0;
cdivLowk = FullSimplify[Block[{$Assumptions = limitAssumptions},
  Limit[ACdiv, k -> 0, Direction -> "FromAbove"]
]];
cdivShift = FullSimplify[cdivLowk - U2];
cdivLiminfShift = FullSimplify[cdivShift /. theta -> Pi/2];

EcdivQuad = -lam eta (k Sin[theta] piK + k Cos[theta] sigma);
cdivMixedPi = FullSimplify[D[EcdivQuad, eta, piK]];
cdivMixedSigma = FullSimplify[D[EcdivQuad, eta, sigma]];
EcpinQuad = 1/2 chi (k Cos[theta] eta)^2;
cpinDensityHessian = FullSimplify[D[EcpinQuad, eta, eta]];
cpinMixedPi = FullSimplify[D[EcpinQuad, eta, piK]];
cpinDelta = FullSimplify[chi k^2 Cos[theta]^2];
ACpin = FullSimplify[AL + cpinDelta];
omega2Cpin = FullSimplify[rho0 k^2 ACpin/m];
AcpinKk = FullSimplify[U2 + (cQ - cL1 + chi Cos[theta]^2) kk + cL2 kk^2];
cpinKstar2Angle = FullSimplify[(cL1 - cQ - chi Cos[theta]^2)/(2 cL2)];
cpinAminAngle = FullSimplify[AcpinKk /. kk -> cpinKstar2Angle];
cpinThresholdC1Angle = FullSimplify[cQ + chi Cos[theta]^2 + 2 Sqrt[U2 cL2]];
cpinKstar2Parallel = FullSimplify[cpinKstar2Angle /. theta -> 0];
cpinKstar2Perp = FullSimplify[cpinKstar2Angle /. theta -> Pi/2];
cpinLowk = FullSimplify[Block[{$Assumptions = limitAssumptions},
  Limit[ACpin, k -> 0, Direction -> "FromAbove"]
]];
cpinShift = FullSimplify[cpinLowk - U2];

fR = FullSimplify[(x^4 - 2 x^2) Exp[-x^2]];
fRprime = FullSimplify[D[fR, x]];
ARofX = FullSimplify[U2 + cQ kR^2 x^2 + AR fR];
rStationarity = FullSimplify[D[ARofX, x]];
rThresholdAR = FullSimplify[-(U2 + cQ kR^2 x^2)/fR];
rThresholdStationarity = FullSimplify[Together[D[rThresholdAR, x] fR^2/(2 x)]];
ARofK = FullSimplify[U2 + cQ k^2 + AR (fR /. x -> k/kR)];
omega2R = FullSimplify[rho0 k^2 ARofK/m];

dot4[u_, v_] := FullSimplify[u.v];
isZeroVec[v_] := And @@ (TrueQ[FullSimplify[# == 0]] & /@ v);
cosModeItems[base_] := Join[({#, 1/2} & /@ base), ({-#, 1/2} & /@ base)];
avgPowerCoeff[base_, power_] := Module[{items = cosModeItems[base], total = 0},
  Do[
    If[isZeroVec[Total[combo[[All, 1]]]], total += Times @@ combo[[All, 2]]],
    {combo, Tuples[items, power]}
  ];
  FullSimplify[total]
];
avgEtaGrad2Coeff[base_] := Module[{items = cosModeItems[base], total = 0},
  Do[
    If[isZeroVec[Total[{etaItem[[1]], p[[1]], qv[[1]]}]],
      total += etaItem[[2]] p[[2]] qv[[2]] (-dot4[p[[1]], qv[[1]]])
    ],
    {etaItem, items}, {p, items}, {qv, items}
  ];
  FullSimplify[total]
];
minimizedPolynomial[energy_, amp_, nonnegative_, rankRules_: {}] := Module[
  {roots, reals, values, order},
  roots = DeleteDuplicates[Join[amp /. Solve[D[energy, amp] == 0, amp], {0}]];
  reals = Select[roots,
    Abs[Im[N[# /. rankRules, 40]]] < 10^-18 &&
      (! nonnegative || Re[N[# /. rankRules, 40]] >= -10^-18) &
  ];
  If[Length[reals] == 0, fail["no real stationary candidates for " <> ToString[energy, InputForm]]];
  values = FullSimplify[energy /. amp -> #] & /@ reals;
  order = First[Ordering[N[values /. rankRules, 40], 1]];
  {FullSimplify[values[[order]]], FullSimplify[reals[[order]]], reals}
];

lamella = {{1, 0, 0, 0}};
triad = {{1, 0, 0, 0}, {-1/2, Sqrt[3]/2, 0, 0}, {-1/2, -Sqrt[3]/2, 0, 0}};
lamAvg2 = avgPowerCoeff[lamella, 2];
lamAvg4 = avgPowerCoeff[lamella, 4];
triAvg2 = avgPowerCoeff[triad, 2];
triAvg3 = avgPowerCoeff[triad, 3];
triAvg4 = avgPowerCoeff[triad, 4];
triEtaGrad2 = avgEtaGrad2Coeff[triad];
assertSame["lamella avg2", lamAvg2, 1/2];
assertSame["lamella avg4", lamAvg4, 3/8];
assertSame["triad avg2", triAvg2, 3/2];
assertSame["triad avg3", triAvg3, 3/2];
assertSame["triad avg4", triAvg4, 45/8];
assertSame["triad eta grad2", triEtaGrad2, 3/4];

FRquad = 1/2 AR fR triAvg2 A^2;
RCubicContribution = FullSimplify[D[FRquad, {A, 3}]];
Ubulk = K rhoSym^5/4;
U3AtRho0 = FullSimplify[D[Ubulk, {rhoSym, 3}] /. rhoSym -> rho0];
GammaREos = FullSimplify[U3AtRho0 triAvg3/6];
GammaRQuantum = FullSimplify[-hbar^2 k^2 triEtaGrad2/(8 m rho0^2)];
GammaR = FullSimplify[GammaREos + GammaRQuantum];
GammaRkk = GammaR /. k^2 -> kk;
GammaRZero = FullSimplify[stripConditional[kk /. First[Block[{$Assumptions = True}, Solve[GammaRkk == 0, kk]]]]];
GammaT = GammaR;
GammaTZero = GammaRZero;
finiteKRegion = TrueQ[FullSimplify[(fR /. x -> 1) == -Exp[-1]]] &&
  TrueQ[FullSimplify[(rThresholdAR /. x -> 1) == E (U2 + cQ kR^2)]];
isotropicShell = TrueQ[FullSimplify[D[ARofK, theta] == 0]];
rAddsCubic = ! TrueQ[FullSimplify[RCubicContribution == 0]];
genericGammaRNonzero = ! TrueQ[FullSimplify[GammaR /. {K -> 1, rho0 -> 1, hbar -> 1, m -> 1, k -> 1}] == 0];

FLam = FullSimplify[-eps lamAvg2 A^2 + beta lamAvg4 A^4];
{FLamMin, FLamAmp, FLamCandidates} = minimizedPolynomial[FLam, A, False, {eps -> 1, beta -> 1}];
qTriNeg = FullSimplify[-eps triAvg2 + zeta (3/2)/2];
uTri = FullSimplify[beta triAvg4];
FTriNeg = FullSimplify[qTriNeg A^2 - gamma A^3 + uTri A^4];
FRank3Neg = FullSimplify[(-3 eps/2 + zeta) A^2 + (45 beta/8) A^4];
FRank4Neg = FullSimplify[(-2 eps + 3 zeta/2) A^2 + (21 beta/2) A^4];
triNonnegativeThreshold = FullSimplify[gamma^2/(4 uTri)];
lamellaOpenRegionZeta = FullSimplify[2 eps + 8 gamma^2/(135 beta)];
FTriPos = FullSimplify[-eps triAvg2 A^2 - gamma A^3 + uTri A^4];
sampleRules = {eps -> 1/100, beta -> 1, gamma -> 1, zeta -> 3};
{sampleLamMin, sampleLamAmp, sampleLamCands} = minimizedPolynomial[FLam, A, False, sampleRules];
{sampleTriNegMin, sampleTriNegAmp, sampleTriNegCands} = minimizedPolynomial[FTriNeg, A, True, sampleRules];
{sampleTriPosMin, sampleTriPosAmp, sampleTriPosCands} = minimizedPolynomial[FTriPos, A, True, sampleRules];
z2Rules = {eps -> 1, beta -> 1, gamma -> 0, zeta -> 0};
baselineTriad = FullSimplify[-eps triAvg2 A^2 - gamma A^3 + uTri A^4];
FTriZ2 = FullSimplify[-eps triAvg2 A^2 + uTri A^4];
FRank3Z2 = FullSimplify[-3 eps A^2/2 + (45 beta/8) A^4];
FRank4Z2 = FullSimplify[-2 eps A^2 + (21 beta/2) A^4];
{triZ2Min, triZ2Amp, triZ2Cands} = minimizedPolynomial[FTriZ2, A, True, z2Rules];
{rank3Z2Min, rank3Z2Amp, rank3Z2Cands} = minimizedPolynomial[FRank3Z2, A, True, z2Rules];
{rank4Z2Min, rank4Z2Amp, rank4Z2Cands} = minimizedPolynomial[FRank4Z2, A, True, z2Rules];
cnullDegrades = TrueQ[FullSimplify[(FTriNeg /. zeta -> 0) - baselineTriad == 0]];
ncPosPass = TrueQ[FullSimplify[(sampleLamMin /. z2Rules) < (triZ2Min /. z2Rules)]] &&
  TrueQ[FullSimplify[(sampleLamMin /. z2Rules) < (rank3Z2Min /. z2Rules)]] &&
  TrueQ[FullSimplify[(sampleLamMin /. z2Rules) < (rank4Z2Min /. z2Rules)]];
ncRtriadPass = TrueQ[triAvg3 == 3/2 && triEtaGrad2 == 3/4];
If[! cnullDegrades, fail["NC-RC-Cnull failed: FTriNeg(zeta=0) did not restore baseline triad"]];
If[! ncPosPass, fail["NC-RC-pos failed: Z2 shell minimizer did not prefer single-k lamella"]];
If[! ncRtriadPass, fail["NC-RC-Rtriad failed: triad cubic/gradient averages changed"]];
omittedGapLower = 3/2;
omittedQLower = FullSimplify[-avg2O eps + zeta omittedGapLower/2];
omittedCoverThresholdZeta = FullSimplify[2 (eps avg2O + gamma^2 c3O^2/(4 u4O))/omittedGapLower];
omittedGapGrowth = TrueQ[Block[{$Assumptions = eps > 0 && avg2O > 0},
  Limit[omittedQLower, zeta -> Infinity]
] == Infinity];

orientationEnergy = FullSimplify[
  1/2 chiAbs gRhoSqSigma pParallel^2 +
  1/2 KPSigma qParallel^2 pParallel^2 +
  1/2 kappaPu (pParallel - OmegaU)^2 +
  2 muBr OmegaU^2
];
orientationStationarity = FullSimplify[D[orientationEnergy, pParallel]];
pParallelDriven = FullSimplify[pParallel /. First[Solve[orientationStationarity == 0, pParallel]]];
pParallelStatic = FullSimplify[pParallelDriven /. OmegaU -> 0];
pParallelStaticSq = FullSimplify[pParallelStatic^2];
pParallelDrivenSample = FullSimplify[
  pParallelDriven /. {kappaPu -> 1, chiAbs -> 1, gRhoSqSigma -> 1, KPSigma -> 1, qParallel -> 0, OmegaU -> 1}
];
quotientEnergy = FullSimplify[1/2 muBr curlU^2 + 1/2 kappaPu mismatch^2];
quotientHessian = FullSimplify[D[quotientEnergy, {{curlU, mismatch}, 2}]];
c6Det = FullSimplify[Det[quotientHessian]];
c6Leading = FullSimplify[quotientHessian[[1, 1]]];
substratePD = TrueQ[c6Det == muBr kappaPu && c6Leading == muBr];
If[! TrueQ[pParallelStatic == 0], fail["static Omega_u=0 minimization unexpectedly left P_parallel nonzero"]];
If[TrueQ[pParallelDrivenSample == 0], fail["driven C.6 fork control did not produce nonzero P_parallel"]];
If[! substratePD, fail["C.6 quotient substrate block is not positive-definite under mu_br,kappa_Pu>0"]];

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
  "A_L_kernel" -> AL,
  "SGoldstone" -> SG,
  "Gmagnitude" -> Gmag,
  "Cdiv_mixed_eta_piT" -> cdivMixedPi,
  "Cdiv_mixed_eta_sigma" -> cdivMixedSigma,
  "Cdiv_delta_A" -> cdivDelta,
  "Cdiv_omega2" -> omega2Cdiv,
  "Cdiv_lowk_limit" -> cdivLowk,
  "Cdiv_lowk_shift" -> cdivShift,
  "Cdiv_liminf_shift" -> cdivLiminfShift,
  "Cpin_mixed_eta_piT" -> cpinMixedPi,
  "Cpin_direct_density_hessian" -> cpinDensityHessian,
  "Cpin_delta_A" -> cpinDelta,
  "Cpin_omega2" -> omega2Cpin,
  "Cpin_kstar2_angle" -> cpinKstar2Angle,
  "Cpin_Amin_angle" -> cpinAminAngle,
  "Cpin_threshold_cL1_angle" -> cpinThresholdC1Angle,
  "Cpin_kstar2_parallel" -> cpinKstar2Parallel,
  "Cpin_kstar2_perp" -> cpinKstar2Perp,
  "Cpin_lowk_limit" -> cpinLowk,
  "Cpin_lowk_shift" -> cpinShift,
  "R_fR" -> fR,
  "R_fRprime" -> fRprime,
  "R_A_of_x" -> ARofX,
  "R_stationarity" -> rStationarity,
  "R_threshold_AR" -> rThresholdAR,
  "R_threshold_stationarity" -> rThresholdStationarity,
  "R_A_of_k" -> ARofK,
  "R_omega2" -> omega2R,
  "R_cubic_contribution" -> RCubicContribution,
  "Gamma_R_EOS_cubic" -> GammaREos,
  "Gamma_R_quantum_pressure_cubic" -> GammaRQuantum,
  "Gamma_R" -> GammaR,
  "Gamma_R_zero_k2" -> GammaRZero,
  "lamella_avg2" -> lamAvg2,
  "lamella_avg4" -> lamAvg4,
  "triad_avg2" -> triAvg2,
  "triad_avg3" -> triAvg3,
  "triad_avg4" -> triAvg4,
  "triad_eta_grad2" -> triEtaGrad2,
  "Cpin_F_lam" -> FLam,
  "Cpin_F_lam_min" -> FLamMin,
  "Cpin_q_tri_neg" -> qTriNeg,
  "Cpin_F_tri_neg" -> FTriNeg,
  "Cpin_F_rank3_neg" -> FRank3Neg,
  "Cpin_F_rank4_neg" -> FRank4Neg,
  "Cpin_tri_nonnegative_threshold" -> triNonnegativeThreshold,
  "Cpin_lamella_open_region_zeta" -> lamellaOpenRegionZeta,
  "Cpin_F_tri_pos" -> FTriPos,
  "Cescape_sample_lamella_min" -> sampleLamMin,
  "Cescape_sample_tri_neg_min" -> sampleTriNegMin,
  "Cpin_positive_sample_tri_min" -> sampleTriPosMin,
  "NC_pos_tri_z2_min" -> triZ2Min,
  "NC_pos_rank3_z2_min" -> rank3Z2Min,
  "NC_pos_rank4_z2_min" -> rank4Z2Min,
  "NC_Cnull_degradation" -> FullSimplify[(FTriNeg /. zeta -> 0) - baselineTriad],
  "Omitted_gap_lower_bound" -> omittedGapLower,
  "Omitted_q_lower_bound" -> omittedQLower,
  "Omitted_cover_threshold_zeta" -> omittedCoverThresholdZeta,
  "Gamma_T" -> GammaT,
  "Gamma_T_zero_k2" -> GammaTZero,
  "C6_pinning_strength" -> chiAbs gRhoSqSigma,
  "C6_frank_stiffness" -> KPSigma qParallel^2,
  "C6_total_orientation_energy" -> orientationEnergy,
  "C6_stationarity_dE_dP_parallel" -> orientationStationarity,
  "C6_P_parallel_driven" -> pParallelDriven,
  "C6_P_parallel_static" -> pParallelStatic,
  "C6_P_parallel_static_sq" -> pParallelStaticSq,
  "C6_P_parallel_driven_sample" -> pParallelDrivenSample,
  "C6_quotient_energy" -> quotientEnergy,
  "C6_hessian_det" -> c6Det,
  "C6_hessian_leading_minor" -> c6Leading
};
Scan[assertEngine[#[[1]], #[[2]]] &, engineKeys];

cdivAdmissionFail = TrueQ[FullSimplify[cdivLiminfShift == -lam^2/(5 K a^2 rho0^5)]];
rNotCodim1 = TrueQ[
  finiteKRegion && isotropicShell && (! rAddsCubic) && genericGammaRNonzero &&
  GammaRZero == 40 K m rho0^4/hbar^2
];
cpinDensityEscape = TrueQ[FullSimplify[qTriNeg /. zeta -> lamellaOpenRegionZeta] == triNonnegativeThreshold];
cpinLightFail = TrueQ[pParallelStaticSq == 0];
cpinSubstrateViable = TrueQ[pParallelStaticSq != 0 && substratePD];
cpinVerdict = Which[
  cpinDensityEscape && cpinSubstrateViable, "RC_S_L_plus_Cpin_CODIM1_CONDITIONAL",
  cpinDensityEscape && cpinLightFail, "RC_S_L_plus_Cpin_FAIL_LIGHT_STARVED",
  True, "INCONCLUSIVE"
];

branchVerdicts = <|
  "S_R_kernel" -> If[rNotCodim1, "RC_S_R_kernel_FAIL_NOT_CODIM1", "INCONCLUSIVE"],
  "S_L_plus_Cdiv" -> If[cdivAdmissionFail, "RC_S_L_plus_Cdiv_FAIL_ADMISSION", "INCONCLUSIVE"],
  "S_L_plus_Cpin" -> cpinVerdict
|>;
verdictValues = Values[branchVerdicts];
topLineVerdict = Which[
  MemberQ[verdictValues, "INCONCLUSIVE"], "INCONCLUSIVE",
  AnyTrue[verdictValues, StringEndsQ[#, "_CODIM1_CONDITIONAL"] &], "RC_DENSITY_SMECTIC_CONDITIONAL",
  AnyTrue[verdictValues, StringEndsQ[#, "_FAIL_LIGHT_STARVED"] &], "RC_DENSITY_SMECTIC_LIGHT_NOGO",
  True, "RC_DENSITY_SMECTIC_NO_GO"
];
If[topLineVerdict =!= sympyResults["top_line_verdict"], fail["top-line verdict disagreement with SymPy JSON"]];

results = <|
  "schema" -> "pathA_25_gateRC_cubic_mathematica/v1",
  "top_line_verdict" -> topLineVerdict,
  "branch_verdicts" -> branchVerdicts,
  "freeze_fidelity" -> <|
    "t0_hash" -> t0Hash,
    "g0_hash" -> g0Hash,
    "t0_bytes_embedded_in_g0" -> True,
    "required_frozen_lines_present" -> True
  |>,
  "checks" -> <|
    "Cdiv_admission_fail" -> cdivAdmissionFail,
    "R_not_codim1" -> rNotCodim1,
    "Cpin_density_escape" -> cpinDensityEscape,
    "Cpin_light_fail" -> cpinLightFail,
    "Cpin_substrate_viable" -> cpinSubstrateViable,
    "NC_RC_Cnull" -> cnullDegrades,
    "NC_RC_pos" -> ncPosPass,
    "NC_RC_Rtriad" -> ncRtriadPass,
    "NC_RC_C6_static_fork" -> cpinLightFail,
    "NC_RC_C6_driven_conditional_fork" -> TrueQ[pParallelDrivenSample != 0 && substratePD],
    "omitted_gap_growth" -> omittedGapGrowth
  |>,
  "engine_agreement" -> <|
    "sympy_json" -> sympyJson,
    "shared_expression_count" -> Length[engineKeys],
    "status" -> "PASS"
  |>
|>;

If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir, CreateIntermediateDirectories -> True]];
Export[jsonOut, results, "JSON"];
Print["PASS pathA_25_gateRC_cubic_mathematica"];
Print[ExportString[<|
  "top_line_verdict" -> topLineVerdict,
  "engine_agreement" -> "PASS",
  "shared_expression_count" -> Length[engineKeys],
  "json" -> jsonOut
|>, "JSON"]];
Exit[0];
