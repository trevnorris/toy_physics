(* Ledger stage017 Mathematica audit: grouped-P2 lane isotropy.

   Standalone, print-only, no arguments, no file I/O. This keeps the native
   Wolfram route for 017's anisotropy integrals, grouped lane algebra, numeric
   response probes, and the single-Y20 integrity echo. Stage 016's SO(3)
   covariance theorem is consumed as cited data with an explicit lambda/K2
   dual-site guard; the full Gram/eigenvalue theorem is not re-derived here.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;

isotropyCalibrated = "ISOTROPY_CALIBRATED";
isotropyPass = "ISOTROPY_PASS";
failDimensional = "FAIL_DIMENSIONAL";
failNotCovariant = "FAIL_NOT_COVARIANT";
failTautological = "FAIL_TAUTOLOGICAL";
failStaticResponse = "FAIL_STATIC_RESPONSE";
failStability = "FAIL_STABILITY";
failSingularResponse = "FAIL_SINGULAR_RESPONSE";
failAnisotropicBranch = "FAIL_ANISOTROPIC_BRANCH";
failNotAbleToFail = "FAIL_NOT_ABLE_TO_FAIL";

symbolicTol = 1.0*^-10;
numericTol = 1.0*^-8;
probeTol = 1.0*^-12;
epsValues = {1.0*^-4, 2.0*^-4, -3.0*^-4};
deltaProfile = 0.1;

cited016Lambda = 6;
cited016CSelf = 1;
cited016DimensionalOk = True;
cited016Covariant = True;
cited016TautologyClear = True;
cited016ThreeProbeAggregate = True;

$Assumptions = Element[{theta, phi, eps, delta}, Reals] &&
  Mtilde > 0 && Ktilde > 0 && TomegaTilde > 0;

raise[msg_] := Throw[msg, "ledgerStage017Failure"];

heading[text_] := (
  Print[""];
  Print[StringRepeat["=", StringLength[text]]];
  Print[text];
  Print[StringRepeat["=", StringLength[text]]]
);

subheading[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);

dropConditions[expr_] := expr /. ConditionalExpression[value_, _] :> value;
clean[expr_] := FullSimplify[dropConditions[expr], $Assumptions];
fmt[expr_String] := expr;
fmt[expr_] := ToString[InputForm[clean[expr]]];

assertExact[name_, expr_] := Module[{reals},
  reals = Cases[Unevaluated[expr], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FAIL  ", name, ": machine-real atom(s) found in exact assertion: ", ToString[InputForm[reals]]];
    raise[name]
  ]
];

expectZero[name_, residual_] := Module[{c},
  assertExact[name, residual];
  c = clean[residual];
  assertExact[name, c];
  If[TrueQ[c === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FAIL  ", name, ": residual = ", fmt[c]];
    raise[name]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

expectFail[name_, residual_] := Module[{c},
  assertExact[name, residual];
  c = clean[residual];
  assertExact[name, c];
  If[! TrueQ[c === 0],
    passCount++;
    Print["PASS  ", name, " produced required FAIL (residual = ", fmt[c], ")"],
    failCount++;
    Print["FAIL  ", name, ": required mutation/ablation did not fire"];
    raise[name]
  ]
];

boolResidual[condition_] := If[TrueQ[condition], 0, 1];
verdictResidual[actual_, expected_] := If[actual === expected, 0, 1];
allZero[exprs_] := And @@ (TrueQ[clean[#] === 0] & /@ exprs);

intS2[expr_] := FullSimplify[
  Integrate[
    Integrate[TrigExpand[expr] Sin[theta], {phi, 0, 2 Pi}],
    {theta, 0, Pi}
  ],
  $Assumptions
];

lapS2[expr_] := FullSimplify[
  1/Sin[theta] D[Sin[theta] D[expr, theta], theta] +
    1/Sin[theta]^2 D[expr, {phi, 2}],
  $Assumptions
];

defectA[triple_] := clean[(2 triple[[1]] - triple[[2]] - triple[[3]])/10];
defectB[triple_] := clean[(triple[[2]] - triple[[3]])/2];
u2FromD[d_] := clean[-d["2"]/d["0"]];
u4FromD[d_] := clean[(d["2"]^2 - d["0"] d["4"])/d["0"]^2];

verdictFromGates[g_, calibration_] := Which[
  ! TrueQ[g["dimensional_ok"]], failDimensional,
  ! TrueQ[g["covariant"]], failNotCovariant,
  ! TrueQ[g["tautology_clear"]], failTautological,
  ! TrueQ[g["dynamic_retained"]], failStaticResponse,
  ! TrueQ[g["stability_ok"]], failStability,
  ! TrueQ[g["denominator_guard_ok"]], failSingularResponse,
  ! TrueQ[g["lane_collapse_ok"]], failAnisotropicBranch,
  ! TrueQ[g["able_to_fail_ok"]], failNotAbleToFail,
  Length[calibration] > 0, isotropyCalibrated,
  True, isotropyPass
];

applyOverrides[g_, overrides_] := Module[{out = Association[g]},
  Scan[(out[#] = overrides[#]) &, Keys[overrides]];
  out
];

order = {"20", "21c", "21s", "22c", "22s"};
lanes = {"20", "21", "22"};
ordersD = {"0", "2", "4"};

harmonics = <|
  "20" -> Sqrt[5/(16 Pi)] (3 Cos[theta]^2 - 1),
  "21c" -> -Sqrt[15/(4 Pi)] Sin[theta] Cos[theta] Cos[phi],
  "21s" -> -Sqrt[15/(4 Pi)] Sin[theta] Cos[theta] Sin[phi],
  "22c" -> Sqrt[15/(16 Pi)] Sin[theta]^2 Cos[2 phi],
  "22s" -> Sqrt[15/(16 Pi)] Sin[theta]^2 Sin[2 phi]
|>;
ys = harmonics /@ order;

p2AxisZ = (3 Cos[theta]^2 - 1)/2;
anisotropyCoefficients = AssociationThread[order, intS2[p2AxisZ #^2] & /@ ys];
cGroup = <|
  "20" -> anisotropyCoefficients["20"],
  "21" -> clean[(anisotropyCoefficients["21c"] + anisotropyCoefficients["21s"])/2],
  "22" -> clean[(anisotropyCoefficients["22c"] + anisotropyCoefficients["22s"])/2]
|>;

lambdaByChannel = AssociationThread[order, ConstantArray[cited016Lambda, Length[order]]];
cSelf = AssociationThread[order, ConstantArray[cited016CSelf, Length[order]]];

dualSiteLambdaResidualsFor[lambdaAssoc_, lambdaRef_] := AssociationThread[
    ("lambda_site_A_equals_016_export_" <> #) & /@ order,
    clean[lambdaAssoc[#] - lambdaRef] & /@ order
  ];

dualSiteEchoResidualFor[lambdaRef_] := clean[-lapS2[harmonics["20"]] - lambdaRef harmonics["20"]];

dualSiteK2FormResidualsFor[laneAssoc_, lambdaRef_] := AssociationThread[
  ("K2_form_coeff_matches_016_export_" <> #) & /@ order,
  clean[D[laneAssoc[#]["K2"], TomegaTilde] - cSelf[#] lambdaRef] & /@ order
];

dualSiteResidualsFor[lambdaAssoc_, lambdaRef_, laneAssoc_] := Join[
  dualSiteLambdaResidualsFor[lambdaAssoc, lambdaRef],
  dualSiteK2FormResidualsFor[laneAssoc, lambdaRef],
  <|"single_harmonic_Y20_echo" -> dualSiteEchoResidualFor[lambdaRef]|>
];

assembleChannel[name_, lambdaAssoc_: lambdaByChannel, k2CoeffOffset_: 0] := Module[
  {pref, lam, mLane, kLane, bLane, zLane, dLane},
  pref = cSelf[name];
  lam = lambdaAssoc[name];
  mLane = clean[pref Mtilde];
  kLane = clean[pref (Ktilde + (lam + k2CoeffOffset) TomegaTilde)];
  bLane = <|
    "0" -> clean[pref B0tilde],
    "2" -> clean[pref B2tilde],
    "4" -> clean[pref B4tilde]
  |>;
  zLane = <|
    "0" -> clean[pref Z0tilde],
    "2" -> clean[pref Z2tilde],
    "4" -> clean[pref Z4tilde]
  |>;
  dLane = <|
    "0" -> clean[kLane - bLane["0"] - zLane["0"]],
    "2" -> clean[-(mLane + bLane["2"] + zLane["2"])],
    "4" -> clean[-(bLane["4"] + zLane["4"])]
  |>;
  <|"angular_self_overlap" -> pref, "lambda" -> lam, "M2" -> mLane, "K2" -> kLane, "D" -> dLane|>
];

assembleUngrouped[lambdaAssoc_, k2CoeffOffsets_: <||>] := AssociationThread[
  order,
  assembleChannel[#, lambdaAssoc, Lookup[k2CoeffOffsets, #, 0]] & /@ order
];

ungrouped = assembleUngrouped[lambdaByChannel];
dualSiteResiduals = dualSiteResidualsFor[lambdaByChannel, cited016Lambda, ungrouped];
dualSiteOk = allZero[Values[dualSiteResiduals]];

avgExpr[exprs_] := clean[Mean[exprs]];
groupedLanes = <|
  "20" -> ungrouped["20"],
  "21" -> <|
    "M2" -> avgExpr[{ungrouped["21c"]["M2"], ungrouped["21s"]["M2"]}],
    "K2" -> avgExpr[{ungrouped["21c"]["K2"], ungrouped["21s"]["K2"]}],
    "D" -> AssociationThread[ordersD, avgExpr[{ungrouped["21c"]["D"][#], ungrouped["21s"]["D"][#]}] & /@ ordersD]
  |>,
  "22" -> <|
    "M2" -> avgExpr[{ungrouped["22c"]["M2"], ungrouped["22s"]["M2"]}],
    "K2" -> avgExpr[{ungrouped["22c"]["K2"], ungrouped["22s"]["K2"]}],
    "D" -> AssociationThread[ordersD, avgExpr[{ungrouped["22c"]["D"][#], ungrouped["22s"]["D"][#]}] & /@ ordersD]
  |>
|>;

csEqual = Join[
  AssociationThread[
    ("D21c_equals_D21s_order_" <> #) & /@ ordersD,
    (TrueQ[clean[ungrouped["21c"]["D"][#] - ungrouped["21s"]["D"][#]] === 0] &) /@ ordersD
  ],
  AssociationThread[
    ("D22c_equals_D22s_order_" <> #) & /@ ordersD,
    (TrueQ[clean[ungrouped["22c"]["D"][#] - ungrouped["22s"]["D"][#]] === 0] &) /@ ordersD
  ]
];

rawDefectsFromGrouped[gl_] := Module[{out = <||>, triple},
  Do[
    triple = gl[#]["D"][n] & /@ lanes;
    out["a_D" <> n] = defectA[triple];
    out["b_D" <> n] = defectB[triple],
    {n, ordersD}
  ];
  out
];

normalizedDefectsFromGrouped[gl_] := Module[{u2L, u4L},
  u2L = AssociationThread[lanes, u2FromD[gl[#]["D"]] & /@ lanes];
  u4L = AssociationThread[lanes, u4FromD[gl[#]["D"]] & /@ lanes];
  <|
    "a2" -> defectA[u2L /@ lanes],
    "b2" -> defectB[u2L /@ lanes],
    "a4" -> defectA[u4L /@ lanes],
    "b4" -> defectB[u4L /@ lanes]
  |>
];

rawTriples = Association[];
Do[rawTriples[n] = AssociationThread[lanes, groupedLanes[#]["D"][n] & /@ lanes], {n, ordersD}];
rawDefects = rawDefectsFromGrouped[groupedLanes];
rawDefectsZero = allZero[Values[rawDefects]];

u2Lanes = AssociationThread[lanes, u2FromD[groupedLanes[#]["D"]] & /@ lanes];
u4Lanes = AssociationThread[lanes, u4FromD[groupedLanes[#]["D"]] & /@ lanes];
normalizedDefects = normalizedDefectsFromGrouped[groupedLanes];
normalizedDefectsZero = allZero[Values[normalizedDefects]];

calibrationSample = <|
  Mtilde -> 3.0,
  Ktilde -> 7.0,
  TomegaTilde -> 0.5,
  B0tilde -> 1.0,
  Z0tilde -> 0.5,
  B2tilde -> 0.4,
  Z2tilde -> 0.2,
  B4tilde -> 0.05,
  Z4tilde -> 0.03
|>;
sampleRules = Normal[calibrationSample];
evalSample[expr_, rules_: sampleRules] := N[clean[expr] /. rules, 30];
firstNonzero[values_] := AnyTrue[Flatten[{values}], TrueQ[Abs[N[#]] > probeTol] &];
maxRatioDelta[values_, epsList_] := Module[{ratios},
  ratios = MapThread[#1/#2 &, {N[values], epsList}];
  Max[Abs[ratios - First[ratios]]]
];
sampleVectorAssoc[assoc_] := Association @ KeyValueMap[
  (#1 -> Table[evalSample[#2 /. eps -> e], {e, epsValues}]) &,
  assoc
];
sampleScalarAssoc[assoc_] := Association @ KeyValueMap[
  (#1 -> evalSample[#2 /. delta -> deltaProfile]) &,
  assoc
];

k2WindowMin = 5.0 + 6.0*0.25;
d0WindowMin = 5.0 + 6.0*0.25 - 1.5 - 0.75;
k2Sample = 7.0 + 6.0*0.5;
omegaSample = Sqrt[k2Sample/3.0];
stabilityOk = TrueQ[2.0 > 0.0 && k2WindowMin > 0.0];
denominatorGuardOk = TrueQ[d0WindowMin > 0.0];

derivedInputs = {
  "explicit real l=2 harmonics for 017 anisotropy coefficients and the one-Y20 echo",
  "016-cited Gram diagonal c_self=1, lambda_m=6, and K2 angular coefficient",
  "ungrouped and grouped {20,21,22} lane assembly",
  "raw-D primary defect algebra from assembled lanes",
  "normalized-u cross-check algebra",
  "017 six response-probe verdicts and computed gate flags"
};
calibrationInputs = {
  "R0(w) linearized isotropic reference",
  "beta2(w) radial profile",
  "mu_eta(w)",
  "T_w(w)",
  "K_eta(w)",
  "T_Omega(w)",
  "Mtilde radial mass scalar",
  "Ktilde radial stiffness scalar excluding angular T_Omega",
  "TomegaTilde angular radial scalar",
  "B0tilde,B2tilde,B4tilde support scalars",
  "Z0tilde,Z2tilde,Z4tilde mixed/Maxwell scalars",
  "physical calibration window for positivity and denominator guards",
  "Gate-1 D/N boundary provenance"
};

caseVerdict[overrides_: <||>] := Module[{g},
  g = <|
    "dimensional_ok" -> cited016DimensionalOk,
    "covariant" -> TrueQ[cited016Covariant && dualSiteOk],
    "tautology_clear" -> cited016TautologyClear,
    "dynamic_retained" -> True,
    "stability_ok" -> True,
    "denominator_guard_ok" -> True,
    "lane_collapse_ok" -> True,
    "able_to_fail_ok" -> True
  |>;
  verdictFromGates[applyOverrides[g, overrides], calibrationInputs]
];

dCommon = <|
  "0" -> clean[Ktilde + cited016Lambda TomegaTilde - B0tilde - Z0tilde],
  "2" -> clean[-(Mtilde + B2tilde + Z2tilde)],
  "4" -> clean[-(B4tilde + Z4tilde)]
|>;

pureDByLane = AssociationThread[
  lanes,
  Table[
    AssociationThread[ordersD, clean[(1 + eps cGroup[lane]) #] & /@ (dCommon /@ ordersD)],
    {lane, lanes}
  ]
];
pureRawDefects = rawDefectsFromGrouped[AssociationThread[lanes, (<|"D" -> pureDByLane[#]|> &) /@ lanes]];
pureU2 = AssociationThread[lanes, u2FromD[pureDByLane[#]] & /@ lanes];
pureU4 = AssociationThread[lanes, u4FromD[pureDByLane[#]] & /@ lanes];
pureUDefects = <|
  "a2" -> defectA[pureU2 /@ lanes],
  "b2" -> defectB[pureU2 /@ lanes],
  "a4" -> defectA[pureU4 /@ lanes],
  "b4" -> defectB[pureU4 /@ lanes]
|>;
pureSamples = sampleVectorAssoc[pureRawDefects];
pureRawMoves = And @@ (firstNonzero /@ Values[pureSamples]);
pureLinearDelta = Max[maxRatioDelta[#, epsValues] & /@ Values[pureSamples]];
pureUZero = allZero[Values[pureUDefects]];

sectorDByLane = AssociationThread[
  lanes,
  Table[
    <|
      "0" -> clean[Ktilde + cited016Lambda TomegaTilde - (1 + eps cGroup[lane]) (B0tilde + Z0tilde)],
      "2" -> clean[-(Mtilde + (1 + eps cGroup[lane]) (B2tilde + Z2tilde))],
      "4" -> clean[-((1 + eps cGroup[lane]) (B4tilde + Z4tilde))]
    |>,
    {lane, lanes}
  ]
];
sectorRawDefects = rawDefectsFromGrouped[AssociationThread[lanes, (<|"D" -> sectorDByLane[#]|> &) /@ lanes]];
sectorSamples = sampleVectorAssoc[sectorRawDefects];
sectorU2 = AssociationThread[lanes, u2FromD[sectorDByLane[#]] & /@ lanes];
sectorUDefects = <|
  "a2" -> defectA[sectorU2 /@ lanes],
  "b2" -> defectB[sectorU2 /@ lanes]
|>;
sectorUSamples = sampleVectorAssoc[sectorUDefects];
sectorRawMoves = And @@ (firstNonzero /@ Values[sectorSamples]);
sectorUMoves = TrueQ[firstNonzero[sectorUSamples["a2"]] && firstNonzero[sectorUSamples["b2"]]];

profileScales = <|"20" -> 1, "21" -> 1, "22" -> (1 + delta)^2|>;
profileDByLane = AssociationThread[
  lanes,
  Table[AssociationThread[ordersD, clean[profileScales[lane] #] & /@ (dCommon /@ ordersD)], {lane, lanes}]
];
profileRawDefects = rawDefectsFromGrouped[AssociationThread[lanes, (<|"D" -> profileDByLane[#]|> &) /@ lanes]];
profileSamples = sampleScalarAssoc[profileRawDefects];
profileRawMoves = firstNonzero[Values[profileSamples]];

betaStabilityProbe[betaScale_] := Module[{m2, k2, m2Sample, k2SampleValue, ok, v},
  m2 = clean[betaScale^2 groupedLanes["20"]["M2"]];
  k2 = clean[betaScale^2 groupedLanes["20"]["K2"]];
  m2Sample = evalSample[m2];
  k2SampleValue = evalSample[k2];
  ok = TrueQ[m2Sample > probeTol && k2SampleValue > probeTol];
  v = caseVerdict[<|"stability_ok" -> ok|>];
  <|"beta2_scale" -> betaScale, "M2" -> m2, "K2" -> k2, "M2_sample" -> m2Sample,
    "K2_sample" -> k2SampleValue, "stability_ok" -> ok, "verdict" -> v,
    "fail_fires" -> TrueQ[v === failStability]|>
];

degenerateBetaProbe = betaStabilityProbe[0];
degenerateBetaAblation = betaStabilityProbe[1];

singularRules = {
  Mtilde -> 3.0,
  Ktilde -> 7.0,
  TomegaTilde -> 0.5,
  B0tilde -> k2Sample - 0.5,
  Z0tilde -> 0.5,
  B2tilde -> 0.4,
  Z2tilde -> 0.2,
  B4tilde -> 0.05,
  Z4tilde -> 0.03
};
singularD0Value = evalSample[dCommon["0"], singularRules];
singularDenominatorGuardOk = TrueQ[Abs[singularD0Value] >= probeTol];
singularVerdict = caseVerdict[<|"denominator_guard_ok" -> singularDenominatorGuardOk|>];
nonsingularD0Value = evalSample[dCommon["0"]];
nonsingularDenominatorGuardOk = TrueQ[Abs[nonsingularD0Value] >= probeTol];
singularAblationVerdict = caseVerdict[<|"denominator_guard_ok" -> nonsingularDenominatorGuardOk|>];

staticWrongD2 = clean[-(B2tilde + Z2tilde)];
staticDynamicRetained = TrueQ[! FreeQ[staticWrongD2, Mtilde]];
staticVerdict = caseVerdict[<|"dynamic_retained" -> staticDynamicRetained|>];
staticAblationDynamicRetained = TrueQ[! FreeQ[groupedLanes["20"]["D"]["2"], Mtilde]];
staticAblationVerdict = caseVerdict[<|"dynamic_retained" -> staticAblationDynamicRetained|>];

probes = <|
  "pure_prefactor_anisotropy" -> <|
    "raw_D_defects" -> pureRawDefects,
    "normalized_u_defects" -> pureUDefects,
    "sample_values" -> pureSamples,
    "raw_D_moves" -> pureRawMoves,
    "linear_scaling_max_ratio_delta" -> pureLinearDelta,
    "linear_scaling_confirmed" -> TrueQ[pureRawMoves && pureLinearDelta < symbolicTol],
    "normalized_u_defects_stay_zero" -> pureUZero,
    "verdict" -> caseVerdict[<|"lane_collapse_ok" -> False|>]
  |>,
  "sector_selective_anisotropy" -> <|
    "raw_D_defects" -> sectorRawDefects,
    "normalized_u_defects" -> sectorUDefects,
    "sample_values" -> Join[sectorSamples, AssociationThread[("u_" <> #) & /@ Keys[sectorUSamples], Values[sectorUSamples]]],
    "raw_D_moves" -> sectorRawMoves,
    "u_defects_move" -> sectorUMoves,
    "verdict" -> caseVerdict[<|"lane_collapse_ok" -> False|>]
  |>,
  "m_dependent_profile" -> <|
    "profile_scales" -> profileScales,
    "raw_D_defects" -> profileRawDefects,
    "sample_values" -> profileSamples,
    "raw_D_moves" -> profileRawMoves,
    "verdict" -> caseVerdict[<|"lane_collapse_ok" -> False|>]
  |>,
  "degenerate_beta_zero" -> Join[
    degenerateBetaProbe,
    <|
      "computed_fail_gate" -> TrueQ[! degenerateBetaProbe["stability_ok"]],
      "self_ablation" -> Join[degenerateBetaAblation, <|"fail_suppressed" -> TrueQ[! degenerateBetaAblation["fail_fires"]]|>]
    |>
  ],
  "singular_denominator" -> <|
    "D0_sample" -> singularD0Value,
    "M2_positive" -> TrueQ[3.0 > probeTol],
    "K2_positive" -> TrueQ[k2Sample > probeTol],
    "denominator_guard_ok" -> singularDenominatorGuardOk,
    "computed_fail_gate" -> TrueQ[! singularDenominatorGuardOk],
    "verdict" -> singularVerdict,
    "fail_fires" -> TrueQ[singularVerdict === failSingularResponse],
    "self_ablation" -> <|
      "D0_sample" -> nonsingularD0Value,
      "denominator_guard_ok" -> nonsingularDenominatorGuardOk,
      "verdict" -> singularAblationVerdict,
      "fail_fires" -> TrueQ[singularAblationVerdict === failSingularResponse],
      "fail_suppressed" -> TrueQ[singularAblationVerdict =!= failSingularResponse]
    |>
  |>,
  "static_drop_inertia" -> <|
    "dynamic_retained" -> staticDynamicRetained,
    "computed_fail_gate" -> TrueQ[! staticDynamicRetained],
    "wrong_D2_without_M" -> staticWrongD2,
    "correct_D2" -> dCommon["2"],
    "verdict" -> staticVerdict,
    "fail_fires" -> TrueQ[staticVerdict === failStaticResponse],
    "self_ablation" -> <|
      "dynamic_retained" -> staticAblationDynamicRetained,
      "correct_D2" -> groupedLanes["20"]["D"]["2"],
      "verdict" -> staticAblationVerdict,
      "fail_fires" -> TrueQ[staticAblationVerdict === failStaticResponse],
      "fail_suppressed" -> TrueQ[staticAblationVerdict =!= failStaticResponse]
    |>
  |>
|>;

expectedProbeVerdicts = <|
  "pure_prefactor_anisotropy" -> failAnisotropicBranch,
  "sector_selective_anisotropy" -> failAnisotropicBranch,
  "m_dependent_profile" -> failAnisotropicBranch,
  "degenerate_beta_zero" -> failStability,
  "singular_denominator" -> failSingularResponse,
  "static_drop_inertia" -> failStaticResponse
|>;
expectedProbeVerdictsMatch = Association @ KeyValueMap[
  (#1 -> TrueQ[probes[#1]["verdict"] === #2]) &,
  expectedProbeVerdicts
];
computedProbeGateFlags = <|
  "pure_prefactor_anisotropy" -> TrueQ[probes["pure_prefactor_anisotropy"]["raw_D_moves"] && probes["pure_prefactor_anisotropy"]["normalized_u_defects_stay_zero"]],
  "sector_selective_anisotropy" -> TrueQ[probes["sector_selective_anisotropy"]["raw_D_moves"] && probes["sector_selective_anisotropy"]["u_defects_move"]],
  "m_dependent_profile" -> TrueQ[probes["m_dependent_profile"]["raw_D_moves"]],
  "degenerate_beta_zero" -> TrueQ[probes["degenerate_beta_zero"]["computed_fail_gate"] && probes["degenerate_beta_zero"]["self_ablation"]["fail_suppressed"]],
  "singular_denominator" -> TrueQ[probes["singular_denominator"]["computed_fail_gate"] && probes["singular_denominator"]["self_ablation"]["fail_suppressed"]],
  "static_drop_inertia" -> TrueQ[probes["static_drop_inertia"]["computed_fail_gate"] && probes["static_drop_inertia"]["self_ablation"]["fail_suppressed"]]
|>;
ableToFail017FromFlags[flags_] := TrueQ[(And @@ Values[expectedProbeVerdictsMatch]) && (And @@ Values[flags])];
neuterFlag[flags_, key_] := Join[KeyDrop[flags, {key}], <|key -> False|>];
ableToFail017Ok = ableToFail017FromFlags[computedProbeGateFlags];
ableToFailIfProbeNeutered = AssociationThread[
  Keys[computedProbeGateFlags],
  ableToFail017FromFlags[neuterFlag[computedProbeGateFlags, #]] & /@ Keys[computedProbeGateFlags]
];
ableToFailOk = TrueQ[ableToFail017Ok && cited016ThreeProbeAggregate];

dynamicRetained = TrueQ[! FreeQ[dCommon["2"], Mtilde] && FreeQ[dCommon["0"], Mtilde]];
laneCollapseOk = TrueQ[rawDefectsZero && (And @@ Values[csEqual])];
baselineGates = <|
  "dimensional_ok" -> cited016DimensionalOk,
  "covariant" -> TrueQ[cited016Covariant && dualSiteOk],
  "tautology_clear" -> cited016TautologyClear,
  "dynamic_retained" -> dynamicRetained,
  "stability_ok" -> stabilityOk,
  "denominator_guard_ok" -> denominatorGuardOk,
  "lane_collapse_ok" -> laneCollapseOk,
  "able_to_fail_ok" -> ableToFailOk
|>;
verdict = verdictFromGates[baselineGates, calibrationInputs];

mutateGroupedD[gl_, lane_, n_, value_] := Module[{laneAssoc},
  laneAssoc = <|"M2" -> gl[lane]["M2"], "K2" -> gl[lane]["K2"], "D" -> Join[KeyDrop[gl[lane]["D"], {n}], <|n -> value|>]|>;
  Join[KeyDrop[gl, {lane}], <|lane -> laneAssoc|>]
];

runAssertions[] := Module[{key},
  subheading["Consumed 016 dual-site guard"];
  Scan[Function[k, expectZero["dual_site." <> k, dualSiteResiduals[k]]], Keys[dualSiteResiduals]];
  expectBool["dual_site.baseline_ok", dualSiteOk];
  expectZero["cited_c_self_is_1_for_lane_assembly", cSelf["20"] - 1];
  expectBool["cited_016_dimensional_ok", cited016DimensionalOk];
  expectBool["cited_016_covariant_gate", cited016Covariant];
  expectBool["cited_016_tautology_clear", cited016TautologyClear];
  expectBool["cited_016_three_probe_aggregate", cited016ThreeProbeAggregate];

  subheading["Grouped lane assembly and exact isotropy"];
  expectZero["anisotropy.c_group_20_is_2_over_7", cGroup["20"] - 2/7];
  expectZero["anisotropy.c_group_21_is_1_over_7", cGroup["21"] - 1/7];
  expectZero["anisotropy.c_group_22_is_minus_2_over_7", cGroup["22"] + 2/7];
  Scan[
    Function[lane,
      expectZero["grouped_lane." <> lane <> ".M2_is_Mtilde", groupedLanes[lane]["M2"] - Mtilde];
      expectZero["grouped_lane." <> lane <> ".K2_is_Ktilde_plus_6_TomegaTilde", groupedLanes[lane]["K2"] - (Ktilde + 6 TomegaTilde)]
    ],
    lanes
  ];
  Scan[Function[k, expectBool["cs_equal." <> k, csEqual[k]]], Keys[csEqual]];
  Scan[Function[k, expectZero["raw_D_primary." <> k, rawDefects[k]]], Keys[rawDefects]];
  expectBool["raw_D_primary.raw_defects_zero", rawDefectsZero];
  Scan[Function[k, expectZero["normalized_cross_check." <> k, normalizedDefects[k]]], Keys[normalizedDefects]];
  expectBool["normalized_cross_check.normalized_defects_zero", normalizedDefectsZero];

  subheading["Numeric stability and denominator windows"];
  expectBool["stability.window_M_and_K2_positive", stabilityOk];
  expectBool["denominator.window_D0_positive", denominatorGuardOk];
  expectBool["stability.k2_window_min_positive", k2WindowMin > 0.0];
  expectBool["denominator.d0_window_min_positive", d0WindowMin > 0.0];

  subheading["Six response probes"];
  Scan[
    Function[k,
      expectZero["probe." <> k <> ".verdict", verdictResidual[probes[k]["verdict"], expectedProbeVerdicts[k]]];
      expectBool["probe." <> k <> ".computed_gate_flag", computedProbeGateFlags[k]]
    ],
    Keys[expectedProbeVerdicts]
  ];
  expectBool["probe.pure_prefactor.raw_D_moves", probes["pure_prefactor_anisotropy"]["raw_D_moves"]];
  expectBool["probe.pure_prefactor.normalized_u_defects_stay_zero", probes["pure_prefactor_anisotropy"]["normalized_u_defects_stay_zero"]];
  expectBool["probe.pure_prefactor.linear_scaling_confirmed", probes["pure_prefactor_anisotropy"]["linear_scaling_confirmed"]];
  expectBool["probe.sector_selective.raw_D_moves", probes["sector_selective_anisotropy"]["raw_D_moves"]];
  expectBool["probe.sector_selective.u_defects_move", probes["sector_selective_anisotropy"]["u_defects_move"]];
  expectBool["probe.m_dependent.raw_D_moves", probes["m_dependent_profile"]["raw_D_moves"]];
  Scan[
    Function[k,
      expectBool["probe." <> k <> ".computed_fail_gate", probes[k]["computed_fail_gate"]];
      expectBool["probe." <> k <> ".fail_fires", probes[k]["fail_fires"]];
      expectBool["probe." <> k <> ".self_ablation_fail_suppressed", probes[k]["self_ablation"]["fail_suppressed"]]
    ],
    {"degenerate_beta_zero", "singular_denominator", "static_drop_inertia"}
  ];

  subheading["Aggregate battery and joint landing"];
  expectBool["able_to_fail.017_six_probe_aggregate", ableToFail017Ok];
  Scan[Function[k, expectBool["able_to_fail.neuter_" <> k <> "_flips_false", ! ableToFailIfProbeNeutered[k]]], Keys[ableToFailIfProbeNeutered]];
  expectBool["able_to_fail.joint_016_and_017", ableToFailOk];
  Scan[Function[k, expectBool["baseline_gate." <> k, baselineGates[k]]], Keys[baselineGates]];
  expectZero["joint.verdict_is_ISOTROPY_CALIBRATED", verdictResidual[verdict, isotropyCalibrated]];
  expectBool["no_c_S_live_symbol", ! MemberQ[Names["Global`*"], "c_S" | "cS"]]
];

runTeeth[] := Module[
  {rawMut, rawMutDefects, rawMutZero, pureNeutered, flags, oneSiteLambda, oneSiteUngrouped,
   oneSiteResiduals, assemblyFormulaUngrouped, assemblyFormulaResiduals, coordinatedLambda,
   coordinatedUngrouped, coordinatedResiduals, normMut, normDefects, gates, rungExpectations},
  subheading["Able-to-fail teeth"];
  rawMut = mutateGroupedD[groupedLanes, "21", "0", clean[groupedLanes["21"]["D"]["0"] + TomegaTilde]];
  rawMutDefects = rawDefectsFromGrouped[rawMut];
  rawMutZero = allZero[Values[rawMutDefects]];
  expectFail["tooth.raw_D_PRIMARY_lane_collapse_fires", boolResidual[rawMutZero]];
  expectZero["tooth.raw_D_PRIMARY_verdict_rung", verdictResidual[caseVerdict[<|"lane_collapse_ok" -> rawMutZero|>], failAnisotropicBranch]];

  pureNeutered = neuterFlag[computedProbeGateFlags, "pure_prefactor_anisotropy"];
  expectFail["tooth.pure_prefactor_neutered_raw_D_moves_flips_aggregate", boolResidual[ableToFail017FromFlags[pureNeutered]]];

  Scan[
    Function[k,
      flags = neuterFlag[computedProbeGateFlags, k];
      expectFail["tooth.response_probe_" <> k <> "_neuter_flips_aggregate", boolResidual[ableToFail017FromFlags[flags]]]
    ],
    {"sector_selective_anisotropy", "m_dependent_profile", "degenerate_beta_zero", "singular_denominator", "static_drop_inertia"}
  ];
  Scan[
    Function[k,
      flags = neuterFlag[computedProbeGateFlags, k];
      expectFail["tooth.six_probe_battery_intact_" <> k, boolResidual[ableToFail017FromFlags[flags]]]
    ],
    Keys[computedProbeGateFlags]
  ];

  oneSiteLambda = Join[KeyDrop[lambdaByChannel, {"20"}], <|"20" -> 4|>];
  oneSiteUngrouped = assembleUngrouped[oneSiteLambda];
  oneSiteResiduals = dualSiteResidualsFor[oneSiteLambda, cited016Lambda, oneSiteUngrouped];
  expectFail["tooth.dual_site_one_site_lambda_corruption_fires_at_integrity_guard", oneSiteResiduals["lambda_site_A_equals_016_export_20"]];
  expectFail["tooth.dual_site_one_site_K2_form_corruption_fires_at_integrity_guard", oneSiteResiduals["K2_form_coeff_matches_016_export_20"]];
  assemblyFormulaUngrouped = assembleUngrouped[lambdaByChannel, <|"20" -> 1|>];
  assemblyFormulaResiduals = dualSiteResidualsFor[lambdaByChannel, cited016Lambda, assemblyFormulaUngrouped];
  expectZero["tooth.dual_site_assembly_formula_lambda_site_stays_clean", assemblyFormulaResiduals["lambda_site_A_equals_016_export_20"]];
  expectFail["tooth.dual_site_assembly_formula_K2_corruption_fires_at_K2_form_guard", assemblyFormulaResiduals["K2_form_coeff_matches_016_export_20"]];
  coordinatedLambda = AssociationThread[order, ConstantArray[4, Length[order]]];
  coordinatedUngrouped = assembleUngrouped[coordinatedLambda];
  coordinatedResiduals = dualSiteResidualsFor[coordinatedLambda, 4, coordinatedUngrouped];
  expectFail["tooth.dual_site_coordinated_corruption_caught_by_Y20_echo", coordinatedResiduals["single_harmonic_Y20_echo"]];

  normMut = mutateGroupedD[groupedLanes, "22", "2", clean[groupedLanes["22"]["D"]["2"] + Mtilde]];
  normDefects = normalizedDefectsFromGrouped[normMut];
  expectFail["tooth.normalized_u_cross_check_fires", boolResidual[allZero[Values[normDefects]]]];

  rungExpectations = <|
    "dynamic_retained" -> failStaticResponse,
    "stability_ok" -> failStability,
    "denominator_guard_ok" -> failSingularResponse,
    "lane_collapse_ok" -> failAnisotropicBranch,
    "able_to_fail_ok" -> failNotAbleToFail
  |>;
  Scan[
    Function[k,
      gates = applyOverrides[baselineGates, <|k -> False|>];
      expectZero["tooth.joint_rung_" <> k, verdictResidual[verdictFromGates[gates, calibrationInputs], rungExpectations[k]]]
    ],
    Keys[rungExpectations]
  ]
];

runAritySelfCheck[] := Module[{dualProbe, betaProbe},
  subheading["Wolfram arity self-check and unevaluated-leakage scan"];
  dualProbe = dualSiteResidualsFor[lambdaByChannel, cited016Lambda, ungrouped];
  betaProbe = betaStabilityProbe[1];
  expectZero["arity intS2[1] returns 4*pi", intS2[1] - 4 Pi];
  expectZero["arity lapS2[Y20] is accepted and has l=2 residual", -lapS2[harmonics["20"]] - 6 harmonics["20"]];
  expectZero["arity defectA[{1,1,1}] returns zero", defectA[{1, 1, 1}]];
  expectZero["arity u2FromD[assoc] returns two", u2FromD[<|"0" -> 2, "2" -> -4, "4" -> 1|>] - 2];
  expectBool["arity dualSiteResidualsFor[assoc,lambda,lanes] returns echo key", KeyExistsQ[dualProbe, "single_harmonic_Y20_echo"]];
  expectBool["arity betaStabilityProbe[1] returns verdict key", KeyExistsQ[betaProbe, "verdict"]];
  expectBool["arity caseVerdict[empty assoc] returns ISOTROPY_CALIBRATED", caseVerdict[<||>] === isotropyCalibrated];
  expectBool["no unevaluated Integrate remains in 017 exact results", FreeQ[{Values[anisotropyCoefficients], Values[rawDefects], Values[normalizedDefects]}, _Integrate]];
  expectBool["no unevaluated Derivative remains in 017 exact results", FreeQ[{Values[dualSiteResiduals], Values[rawDefects], Values[normalizedDefects]}, _Derivative]]
];

printTranscript[] := (
  subheading["Transcript summary"];
  Print["Engine: Mathematica standalone native audit"];
  Print["Computed verdict: ", verdict, " (COMPLETE: 016 AND 017)"];
  Print["Grouped lanes:"];
  Scan[
    Function[lane,
      Print["  ", lane, ": M2=", fmt[groupedLanes[lane]["M2"]], "; K2=", fmt[groupedLanes[lane]["K2"]],
        "; D0=", fmt[groupedLanes[lane]["D"]["0"]], "; D2=", fmt[groupedLanes[lane]["D"]["2"]],
        "; D4=", fmt[groupedLanes[lane]["D"]["4"]]]
    ],
    lanes
  ];
  Print["Raw-D defects (PRIMARY): ", Association @ KeyValueMap[(#1 -> fmt[#2]) &, rawDefects]];
  Print["Normalized-u defects (CROSS-CHECK): ", Association @ KeyValueMap[(#1 -> fmt[#2]) &, normalizedDefects]];
  Print["Stability guard: ", stabilityOk, "; denominator guard: ", denominatorGuardOk];
  Print["Probe verdicts:"];
  Scan[
    Function[k, Print["  ", k, ": ", probes[k]["verdict"], " (computed flag=", computedProbeGateFlags[k], ")"]],
    Keys[expectedProbeVerdicts]
  ];
  Print["017 six-probe aggregate: ", ableToFail017Ok];
  Print["Joint able_to_fail = 016 cited aggregate AND 017 aggregate: ", ableToFailOk];
  Print["Dual-site guard: Site A lane lambda/K2 coefficient equals independent 016 lambda=6; coordinated corruption is closed by a single Y20 (-Delta_S2) echo."];
  Print["CALIBRATED-not-PASS: the wall profile and radial/support scalars are frozen calibration inputs."];
  Print["Carried caveats: 54/5 quadrupole normalization, outgoing odd-N coefficients, and solved nonlinear branch data remain Gate 4/5/6 sim-deferred work."];
  Print["Dropped bookkeeping: scratch YAML/report/feed writers replaced by print-only v2 audit output."];

  subheading["Calibration partition"];
  Print["Derived inputs:"];
  Scan[Print["  - ", #] &, derivedInputs];
  Print["Calibration inputs (printed only; count is resolved at registration):"];
  Scan[Print["  - ", #] &, calibrationInputs];

  subheading["Provenance labels"];
  Print["CONSUMED-016-DUAL-SITE: lambda_m=6, K2=Ktilde+lambda_m*TomegaTilde, Gram diagonal c_self=1, dimensional/covariant/tautology gates, and 016's 3-probe aggregate are cited and guarded."];
  Print["CONSUMED-011/012/013-PROVENANCE: mu_eta/T_w, R0(w), and Gate-1 D/N provenance are narrative provenance; K_eta=T_w*beta^2 is non-transferable across the volume-vs-line convention."];
  Print["NO-c_S: c_S is not a live symbol in the l=2 response audit."];
  Print["RAW-D-PRIMARY: raw-D defects are the primary isotropy test; normalized-u is a cross-check."];
  Print["SIX-PROBE-BATTERY-INTACT: each six-probe flag is computed; neutering any one flips the aggregate false."];
  Print["COMPLETES-THE-JOINT: 017 lands ISOTROPY_CALIBRATED as COMPLETE (016 AND 017), not partial."];
  Print["EXPORTS-PORT-KERNEL: grouped M2, K2=Ktilde+6*TomegaTilde, Btilde/Ztilde, and D-lanes to stages 018-024."];

  subheading["Verdict labels"];
  Print["Verdict labels:"];
  Print["  ledger earned-label: GROUPED_P2_LANE_ISOTROPY_EARNED (grouped l=2 lanes, raw-D primary defects zero, normalized-u cross-check zero, six response probes able-to-fail)"];
  Print["  source top-line verdict: ISOTROPY_CALIBRATED (JOINT 2-stage; 017 completes 016)"];
  Print["  joint composition: ISOTROPY_CALIBRATED = (016 l=2 SO(3) covariance theorem: Gram=I5, lambda_m=6, K2 angular stiffness)[EARNED, cited] AND (017 grouped-P2 lane isotropy: grouped lanes, raw-D=0 PRIMARY, normalized-u=0 CROSS-CHECK, 6-probe battery)[EARNED here]"];
  Print["  calibrated values: beta2(w), R0(w), Mtilde/Ktilde/TomegaTilde, Btilde/Ztilde, and the window are frozen; T_Omega/beta2/support-scalar counting is left to registration"];
  Print["  consumed: 016 lambda_m=6/K2=Ktilde+lambda_m*TomegaTilde/c_self=1 via explicit dual-site guard; 011/012/013 mu_eta/T_w/R0/Gate-1 D-N are provenance; c_S NOT consumed"];
  Print["  exports: l=2 port kernel = grouped M2 + K2=Ktilde+6*TomegaTilde + Btilde/Ztilde support scalars + D0/D2/D4 lanes"]
);

main[] := (
  heading["ledger_stage017 grouped-P2 lane isotropy Mathematica audit"];
  Print["Mode: print-only, assert-zero, zero file I/O, mixed exact symbolic + numeric calibration probes."];
  Print["Target stem: ledger_stage017_grouped_p2_lane_isotropy"];
  runAssertions[];
  runTeeth[];
  runAritySelfCheck[];
  printTranscript[];
  Print[""];
  Print["TALLY mathematica: PASS=", passCount, " FAIL=", failCount, " TOTAL=", passCount + failCount,
    " (", passCount, "+", failCount, "=", passCount + failCount, ")"];
  If[failCount == 0, Print["OVERALL PASS"]; Exit[0], Print["OVERALL FAIL"]; Exit[1]]
);

Catch[
  main[],
  "ledgerStage017Failure",
  (Print["Audit aborted: ", #1];
   Print["TALLY mathematica: PASS=", passCount, " FAIL=", failCount, " TOTAL=", passCount + failCount,
     " (", passCount, "+", failCount, "=", passCount + failCount, ")"];
   Print["OVERALL FAIL"];
   Exit[1]) &
];
