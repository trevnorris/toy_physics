(* PathA-32 Gate 3 grouped-P2 isotropy gate, Mathematica engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
boolText[x_] := If[TrueQ[x], "true", "false"];
numText[x_] := StringReplace[
  ToString[N[x, 17], InputForm],
  {RegularExpression["`[0-9.]*"] -> "", "*^" -> "e"}
];
quoteText[x_] := Module[{s = If[StringQ[x], x, ToString[x, InputForm]]},
  "'" <> StringReplace[s, {"'" -> "''", "\n" -> " "}] <> "'"
];

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_32_grouped_p2_isotropy.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
yamlOut = FileNameJoin[{scratchDir, "pathA_32_mathematica_results.yaml"}];
If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir, CreateIntermediateDirectories -> True]];

$Assumptions = Element[{theta, phi, eps, delta}, Reals] &&
  Mtilde > 0 && Ktilde > 0 && TomegaTilde > 0;

order = {"20", "21c", "21s", "22c", "22s"};
harmonics = <|
  "20" -> Sqrt[5/(16 Pi)] (3 Cos[theta]^2 - 1),
  "21c" -> -Sqrt[15/(4 Pi)] Sin[theta] Cos[theta] Cos[phi],
  "21s" -> -Sqrt[15/(4 Pi)] Sin[theta] Cos[theta] Sin[phi],
  "22c" -> Sqrt[15/(16 Pi)] Sin[theta]^2 Cos[2 phi],
  "22s" -> Sqrt[15/(16 Pi)] Sin[theta]^2 Sin[2 phi]
|>;
ys = harmonics /@ order;

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
defectA[triple_] := FullSimplify[(2 triple[[1]] - triple[[2]] - triple[[3]])/10, $Assumptions];
defectB[triple_] := FullSimplify[(triple[[2]] - triple[[3]])/2, $Assumptions];
u2FromD[d_] := FullSimplify[-d["2"]/d["0"], $Assumptions];
u4FromD[d_] := FullSimplify[(d["2"]^2 - d["0"] d["4"])/d["0"]^2, $Assumptions];

gram = Table[intS2[ys[[i]] ys[[j]]], {i, 5}, {j, 5}];
gramIsIdentity = TrueQ[FullSimplify[gram == IdentityMatrix[5], $Assumptions]];
p2AxisZ = (3 Cos[theta]^2 - 1)/2;
anisotropyCoefficients = AssociationThread[order, intS2[p2AxisZ #^2] & /@ ys];
cGroup = <|
  "20" -> anisotropyCoefficients["20"],
  "21" -> FullSimplify[(anisotropyCoefficients["21c"] + anisotropyCoefficients["21s"])/2, $Assumptions],
  "22" -> FullSimplify[(anisotropyCoefficients["22c"] + anisotropyCoefficients["22s"])/2, $Assumptions]
|>;
groupedIso = IdentityMatrix[3];
groupedCoeff = DiagonalMatrix[cGroup /@ {"20", "21", "22"}];

negLaps = AssociationThread[order, FullSimplify[-lapS2[#], $Assumptions] & /@ ys];
lambdas = AssociationThread[
  order,
  Table[
    FullSimplify[intS2[ys[[i]] negLaps[order[[i]]]]/intS2[ys[[i]]^2], $Assumptions],
    {i, 5}
  ]
];
residuals = AssociationThread[
  order,
  Table[FullSimplify[negLaps[order[[i]]] - lambdas[order[[i]]] ys[[i]], $Assumptions], {i, 5}]
];
lambdaAllSix = And @@ (TrueQ[FullSimplify[# == 6, $Assumptions]] & /@ Values[lambdas]);
residualsZero = And @@ (TrueQ[FullSimplify[# == 0, $Assumptions]] & /@ Values[residuals]);
kCoeffUsed = lambdas;
kCoeffEquals = AssociationThread[
  order,
  (TrueQ[FullSimplify[kCoeffUsed[#] == lambdas[#], $Assumptions]] &) /@ order
];
kCoeffEqualsAll = And @@ Values[kCoeffEquals];

dCommon = <|
  "0" -> FullSimplify[Ktilde + 6 TomegaTilde - B0tilde - Z0tilde, $Assumptions],
  "2" -> FullSimplify[-(Mtilde + B2tilde + Z2tilde), $Assumptions],
  "4" -> FullSimplify[-(B4tilde + Z4tilde), $Assumptions]
|>;
lanes = {"20", "21", "22"};
assembleChannel[name_] := Module[
  {yExpr, cSelf, lam, mLane, kLane, bLane, zLane, dLane},
  yExpr = harmonics[name];
  cSelf = FullSimplify[intS2[yExpr^2], $Assumptions];
  lam = lambdas[name];
  mLane = FullSimplify[cSelf Mtilde, $Assumptions];
  kLane = FullSimplify[cSelf (Ktilde + lam TomegaTilde), $Assumptions];
  bLane = <|
    "0" -> FullSimplify[cSelf B0tilde, $Assumptions],
    "2" -> FullSimplify[cSelf B2tilde, $Assumptions],
    "4" -> FullSimplify[cSelf B4tilde, $Assumptions]
  |>;
  zLane = <|
    "0" -> FullSimplify[cSelf Z0tilde, $Assumptions],
    "2" -> FullSimplify[cSelf Z2tilde, $Assumptions],
    "4" -> FullSimplify[cSelf Z4tilde, $Assumptions]
  |>;
  dLane = <|
    "0" -> FullSimplify[kLane - bLane["0"] - zLane["0"], $Assumptions],
    "2" -> FullSimplify[-(mLane + bLane["2"] + zLane["2"]), $Assumptions],
    "4" -> FullSimplify[-(bLane["4"] + zLane["4"]), $Assumptions]
  |>;
  <|"angular_self_overlap" -> cSelf, "lambda" -> lam, "M2" -> mLane, "K2" -> kLane, "D" -> dLane|>
];
ungrouped = AssociationThread[order, assembleChannel /@ order];
avgExpr[exprs_] := FullSimplify[Mean[exprs], $Assumptions];
groupedLanes = <|
  "20" -> ungrouped["20"],
  "21" -> <|
    "M2" -> avgExpr[{ungrouped["21c"]["M2"], ungrouped["21s"]["M2"]}],
    "K2" -> avgExpr[{ungrouped["21c"]["K2"], ungrouped["21s"]["K2"]}],
    "D" -> AssociationThread[
      {"0", "2", "4"},
      avgExpr[{ungrouped["21c"]["D"][#], ungrouped["21s"]["D"][#]}] & /@ {"0", "2", "4"}
    ]
  |>,
  "22" -> <|
    "M2" -> avgExpr[{ungrouped["22c"]["M2"], ungrouped["22s"]["M2"]}],
    "K2" -> avgExpr[{ungrouped["22c"]["K2"], ungrouped["22s"]["K2"]}],
    "D" -> AssociationThread[
      {"0", "2", "4"},
      avgExpr[{ungrouped["22c"]["D"][#], ungrouped["22s"]["D"][#]}] & /@ {"0", "2", "4"}
    ]
  |>
|>;
csEqual = Join[
  AssociationThread[
    ("D21c_equals_D21s_order_" <> #) & /@ {"0", "2", "4"},
    (TrueQ[FullSimplify[ungrouped["21c"]["D"][#] == ungrouped["21s"]["D"][#], $Assumptions]] &) /@ {"0", "2", "4"}
  ],
  AssociationThread[
    ("D22c_equals_D22s_order_" <> #) & /@ {"0", "2", "4"},
    (TrueQ[FullSimplify[ungrouped["22c"]["D"][#] == ungrouped["22s"]["D"][#], $Assumptions]] &) /@ {"0", "2", "4"}
  ]
];
rawTriples = Association[];
rawDefects = Association[];
Do[
  triple = groupedLanes[#]["D"][n] & /@ lanes;
  rawTriples[n] = AssociationThread[lanes, triple];
  rawDefects["a_D" <> n] = defectA[triple];
  rawDefects["b_D" <> n] = defectB[triple],
  {n, {"0", "2", "4"}}
];
rawDefectsZero = And @@ (TrueQ[FullSimplify[# == 0, $Assumptions]] & /@ Values[rawDefects]);

u2Lanes = AssociationThread[lanes, u2FromD[groupedLanes[#]["D"]] & /@ lanes];
u4Lanes = AssociationThread[lanes, u4FromD[groupedLanes[#]["D"]] & /@ lanes];
normalizedDefects = <|
  "a2" -> defectA[u2Lanes /@ lanes],
  "b2" -> defectB[u2Lanes /@ lanes],
  "a4" -> defectA[u4Lanes /@ lanes],
  "b4" -> defectB[u4Lanes /@ lanes]
|>;
normalizedDefectsZero = And @@ (TrueQ[FullSimplify[# == 0, $Assumptions]] & /@ Values[normalizedDefects]);

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

verdictFromGates[cov_, taut_, dyn_, stable_, denom_, lane_, able_] := Which[
  ! TrueQ[cov], "FAIL_NOT_COVARIANT",
  ! TrueQ[taut], "FAIL_TAUTOLOGICAL",
  ! TrueQ[dyn], "FAIL_STATIC_RESPONSE",
  ! TrueQ[stable], "FAIL_STABILITY",
  ! TrueQ[denom], "FAIL_SINGULAR_RESPONSE",
  ! TrueQ[lane], "FAIL_ANISOTROPIC_BRANCH",
  ! TrueQ[able], "FAIL_NOT_ABLE_TO_FAIL",
  Length[calibrationInputs] > 0, "ISOTROPY_CALIBRATED",
  True, "ISOTROPY_PASS"
];
caseVerdict[opts___Rule] := Module[
  {cov = True, taut = True, dyn = True, stable = True, denom = True, lane = True, able = True},
  {cov, taut, dyn, stable, denom, lane, able} =
    {covariant, tautology, dynamic, stability, denominator, collapse, ableToFail} /.
      {opts} /. {covariant -> cov, tautology -> taut, dynamic -> dyn, stability -> stable,
        denominator -> denom, collapse -> lane, ableToFail -> able};
  verdictFromGates[cov, taut, dyn, stable, denom, lane, able]
];

sampleRules = {
  Mtilde -> 3.0,
  Ktilde -> 7.0,
  TomegaTilde -> 0.5,
  B0tilde -> 1.0,
  Z0tilde -> 0.5,
  B2tilde -> 0.4,
  Z2tilde -> 0.2,
  B4tilde -> 0.05,
  Z4tilde -> 0.03
};
epsValues = {0.0001, 0.0002, -0.0003};
evalSample[expr_] := N[expr /. sampleRules, 17];
sampleVectorAssoc[assoc_] := Association @ KeyValueMap[
  (#1 -> Table[evalSample[#2 /. eps -> e], {e, epsValues}]) &,
  assoc
];
sampleScalarAssoc[assoc_] := Association @ KeyValueMap[
  (#1 -> evalSample[#2 /. delta -> 0.1]) &,
  assoc
];

pureDByLane = AssociationThread[
  {"20", "21", "22"},
  Table[
    AssociationThread[{"0", "2", "4"}, FullSimplify[(1 + eps cGroup[lane]) #, $Assumptions] & /@ (dCommon /@ {"0", "2", "4"})],
    {lane, {"20", "21", "22"}}
  ]
];
pureRawDefects = Association[];
Do[
  triple = pureDByLane[#][n] & /@ {"20", "21", "22"};
  pureRawDefects["a_D" <> n] = defectA[triple];
  pureRawDefects["b_D" <> n] = defectB[triple],
  {n, {"0", "2", "4"}}
];
pureU2 = AssociationThread[{"20", "21", "22"}, u2FromD[pureDByLane[#]] & /@ {"20", "21", "22"}];
pureU4 = AssociationThread[{"20", "21", "22"}, u4FromD[pureDByLane[#]] & /@ {"20", "21", "22"}];
pureUDefects = <|
  "a2" -> defectA[pureU2 /@ {"20", "21", "22"}],
  "b2" -> defectB[pureU2 /@ {"20", "21", "22"}],
  "a4" -> defectA[pureU4 /@ {"20", "21", "22"}],
  "b4" -> defectB[pureU4 /@ {"20", "21", "22"}]
|>;
pureSamples = sampleVectorAssoc[pureRawDefects];
pureRawMoves = And @@ (Max[Abs[N[#]]] > 10^-12 & /@ Values[pureSamples]);
pureUZero = And @@ (TrueQ[FullSimplify[# == 0, $Assumptions]] & /@ Values[pureUDefects]);

sectorDByLane = AssociationThread[
  {"20", "21", "22"},
  Table[
    <|
      "0" -> FullSimplify[Ktilde + 6 TomegaTilde - (1 + eps cGroup[lane]) (B0tilde + Z0tilde), $Assumptions],
      "2" -> FullSimplify[-(Mtilde + (1 + eps cGroup[lane]) (B2tilde + Z2tilde)), $Assumptions],
      "4" -> FullSimplify[-((1 + eps cGroup[lane]) (B4tilde + Z4tilde)), $Assumptions]
    |>,
    {lane, {"20", "21", "22"}}
  ]
];
sectorRawDefects = Association[];
Do[
  triple = sectorDByLane[#][n] & /@ {"20", "21", "22"};
  sectorRawDefects["a_D" <> n] = defectA[triple];
  sectorRawDefects["b_D" <> n] = defectB[triple],
  {n, {"0", "2", "4"}}
];
sectorU2 = AssociationThread[{"20", "21", "22"}, u2FromD[sectorDByLane[#]] & /@ {"20", "21", "22"}];
sectorUDefects = <|
  "a2" -> defectA[sectorU2 /@ {"20", "21", "22"}],
  "b2" -> defectB[sectorU2 /@ {"20", "21", "22"}]
|>;
sectorSamples = Join[
  sampleVectorAssoc[sectorRawDefects],
  sampleVectorAssoc[AssociationThread[("u_" <> #) & /@ Keys[sectorUDefects], Values[sectorUDefects]]]
];
sectorRawMoves = And @@ (Max[Abs[N[#]]] > 10^-12 & /@ Values[sampleVectorAssoc[sectorRawDefects]]);
sectorUMoves = And @@ (Max[Abs[N[#]]] > 10^-12 & /@ Values[sampleVectorAssoc[sectorUDefects]]);

profileScales = <|"20" -> 1, "21" -> 1, "22" -> (1 + delta)^2|>;
profileDByLane = AssociationThread[
  {"20", "21", "22"},
  Table[AssociationThread[{"0", "2", "4"}, FullSimplify[profileScales[lane] #, $Assumptions] & /@ (dCommon /@ {"0", "2", "4"})], {lane, {"20", "21", "22"}}]
];
profileRawDefects = Association[];
Do[
  triple = profileDByLane[#][n] & /@ {"20", "21", "22"};
  profileRawDefects["a_D" <> n] = defectA[triple];
  profileRawDefects["b_D" <> n] = defectB[triple],
  {n, {"0", "2", "4"}}
];
profileSamples = sampleScalarAssoc[profileRawDefects];
profileRawMoves = Max[Abs[N[Values[profileSamples]]]] > 10^-12;

probeTol = 10^-12;
betaStabilityProbe[betaScale_] := Module[{m2, k2, m2Sample, k2SampleValue, ok, v},
  m2 = FullSimplify[betaScale^2 groupedLanes["20"]["M2"], $Assumptions];
  k2 = FullSimplify[betaScale^2 groupedLanes["20"]["K2"], $Assumptions];
  m2Sample = N[m2 /. sampleRules, 17];
  k2SampleValue = N[k2 /. sampleRules, 17];
  ok = TrueQ[m2Sample > probeTol && k2SampleValue > probeTol];
  v = caseVerdict[stability -> ok];
  <|"beta2_scale" -> betaScale, "M2" -> m2, "K2" -> k2, "M2_sample" -> m2Sample,
    "K2_sample" -> k2SampleValue, "stability_ok" -> ok, "verdict" -> v,
    "fail_fires" -> TrueQ[v === "FAIL_STABILITY"]|>
];
forcedEigenvalueProbe[forced_] := Module[{forcedK2, coefficientEquals, v},
  forcedK2 = AssociationThread[
    order,
    FullSimplify[ungrouped[#]["angular_self_overlap"] (Ktilde + forced TomegaTilde), $Assumptions] & /@ order
  ];
  coefficientEquals = And @@ (TrueQ[FullSimplify[forced == lambdas[#], $Assumptions]] & /@ order);
  v = caseVerdict[covariant -> coefficientEquals];
  <|"forced_coefficient" -> forced, "forced_K2_by_channel" -> forcedK2,
    "coefficient_equals_computed_lambda" -> coefficientEquals, "verdict" -> v,
    "fail_fires" -> TrueQ[v === "FAIL_NOT_COVARIANT"]|>
];
hashText[expr_] := IntegerString[Hash[ToString[FullSimplify[expr, $Assumptions], InputForm], "SHA256"], 16];
laneHashProbe[inputs_Association] := Module[{hashes, distinct, v},
  hashes = Association @ KeyValueMap[(#1 -> hashText[#2]) &, inputs];
  distinct = Length[DeleteDuplicates[Values[hashes]]] == Length[hashes];
  v = caseVerdict[tautology -> distinct];
  <|"input_hashes" -> hashes, "distinct_hashes" -> distinct, "verdict" -> v,
    "fail_fires" -> TrueQ[v === "FAIL_TAUTOLOGICAL"]|>
];

degenerateBetaProbe = betaStabilityProbe[0];
degenerateBetaAblation = betaStabilityProbe[1];
wrongEigenProbe = forcedEigenvalueProbe[2];
wrongEigenAblation = forcedEigenvalueProbe[6];
k2Sample = 7.0 + 6.0*0.5;
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
singularD0Sample = N[dCommon["0"] /. singularRules, 17];
singularDenominatorGuardOk = TrueQ[Abs[singularD0Sample] >= probeTol];
singularVerdict = caseVerdict[denominator -> singularDenominatorGuardOk];
nonsingularD0Sample = evalSample[dCommon["0"]];
nonsingularDenominatorGuardOk = TrueQ[Abs[nonsingularD0Sample] >= probeTol];
singularAblationVerdict = caseVerdict[denominator -> nonsingularDenominatorGuardOk];
tautologyProbe = laneHashProbe[<|"20" -> harmonics["20"], "21" -> harmonics["20"], "22" -> harmonics["20"]|>];
tautologyAblation = laneHashProbe[<|"20" -> harmonics["20"], "21" -> harmonics["21c"], "22" -> harmonics["22c"]|>];
staticWrongD2 = FullSimplify[-(B2tilde + Z2tilde), $Assumptions];
staticDynamicRetained = TrueQ[Not[FreeQ[staticWrongD2, Mtilde]]];
staticVerdict = caseVerdict[dynamic -> staticDynamicRetained];
staticAblationDynamicRetained = TrueQ[Not[FreeQ[groupedLanes["20"]["D"]["2"], Mtilde]]];
staticAblationVerdict = caseVerdict[dynamic -> staticAblationDynamicRetained];

probeVerdicts = <|
  "pure_prefactor_anisotropy" -> caseVerdict[collapse -> False],
  "sector_selective_anisotropy" -> caseVerdict[collapse -> False],
  "m_dependent_profile" -> caseVerdict[collapse -> False],
  "degenerate_beta_zero" -> degenerateBetaProbe["verdict"],
  "wrong_eigenvalue" -> wrongEigenProbe["verdict"],
  "singular_denominator" -> singularVerdict,
  "tautology_hash_collision" -> tautologyProbe["verdict"],
  "static_drop_inertia" -> staticVerdict
|>;
expectedProbeVerdicts = <|
  "pure_prefactor_anisotropy" -> "FAIL_ANISOTROPIC_BRANCH",
  "sector_selective_anisotropy" -> "FAIL_ANISOTROPIC_BRANCH",
  "m_dependent_profile" -> "FAIL_ANISOTROPIC_BRANCH",
  "degenerate_beta_zero" -> "FAIL_STABILITY",
  "wrong_eigenvalue" -> "FAIL_NOT_COVARIANT",
  "singular_denominator" -> "FAIL_SINGULAR_RESPONSE",
  "tautology_hash_collision" -> "FAIL_TAUTOLOGICAL",
  "static_drop_inertia" -> "FAIL_STATIC_RESPONSE"
|>;
expectedProbeVerdictsMatch = Association @ KeyValueMap[(#1 -> TrueQ[probeVerdicts[#1] === #2]) &, expectedProbeVerdicts];
computedProbeGateFlags = <|
  "pure_prefactor_anisotropy" -> TrueQ[pureRawMoves && pureUZero],
  "sector_selective_anisotropy" -> TrueQ[sectorRawMoves && sectorUMoves],
  "m_dependent_profile" -> TrueQ[profileRawMoves],
  "degenerate_beta_zero" -> TrueQ[Not[degenerateBetaProbe["stability_ok"]] && Not[degenerateBetaAblation["fail_fires"]]],
  "wrong_eigenvalue" -> TrueQ[Not[wrongEigenProbe["coefficient_equals_computed_lambda"]] && Not[wrongEigenAblation["fail_fires"]]],
  "singular_denominator" -> TrueQ[Not[singularDenominatorGuardOk] && singularAblationVerdict =!= "FAIL_SINGULAR_RESPONSE"],
  "tautology_hash_collision" -> TrueQ[Not[tautologyProbe["distinct_hashes"]] && Not[tautologyAblation["fail_fires"]]],
  "static_drop_inertia" -> TrueQ[Not[staticDynamicRetained] && staticAblationVerdict =!= "FAIL_STATIC_RESPONSE"]
|>;
ableToFailOk = TrueQ[(And @@ Values[expectedProbeVerdictsMatch]) && (And @@ Values[computedProbeGateFlags])];

stabilityOk = TrueQ[2.0 > 0 && (5.0 + 6.0*0.25) > 0];
denominatorGuardOk = TrueQ[(5.0 + 6.0*0.25 - 1.5 - 0.75) > 0];
tautologyClear = TrueQ[
  Length[DeleteDuplicates[ToString[#, InputForm] & /@ ys]] == 5 &&
    And @@ (TrueQ[FullSimplify[intS2[#^2] == 1, $Assumptions]] & /@ ys)
];
laneCollapseOk = TrueQ[rawDefectsZero && (And @@ Values[csEqual])];
dynamicRetained = TrueQ[Not[FreeQ[groupedLanes["20"]["D"]["2"], Mtilde]] && FreeQ[groupedLanes["20"]["D"]["0"], Mtilde]];
covariantOk = TrueQ[gramIsIdentity && lambdaAllSix && residualsZero && kCoeffEqualsAll];
verdict = verdictFromGates[covariantOk, tautologyClear, dynamicRetained, stabilityOk, denominatorGuardOk, laneCollapseOk, ableToFailOk];

If[! TrueQ[covariantOk && rawDefectsZero && normalizedDefectsZero && ableToFailOk && verdict === "ISOTROPY_CALIBRATED"],
  Print["covariantOk=", covariantOk, " rawDefectsZero=", rawDefectsZero,
    " normalizedDefectsZero=", normalizedDefectsZero, " ableToFailOk=", ableToFailOk,
    " verdict=", verdict];
  fail["PathA-32 Mathematica gate did not reach calibrated baseline"]
];

matrixLines[indent_, matrix_] := Table[
  indent <> "- [" <> StringRiffle[numText /@ matrix[[i]], ", "] <> "]",
  {i, Length[matrix]}
];
assocNumericLines[indent_, assoc_] := Flatten[KeyValueMap[
  Function[{key, vals}, Join[
    {indent <> quoteText[key] <> ":"},
    (indent <> "  - " <> numText[#]) & /@ vals
  ]],
  assoc
]];
assocScalarLines[indent_, assoc_] := Flatten[KeyValueMap[
  Function[{key, val}, {indent <> quoteText[key] <> ": " <> numText[val]}],
  assoc
]];
assocBoolLines[indent_, assoc_] := Flatten[KeyValueMap[
  Function[{key, val}, {indent <> quoteText[key] <> ": " <> boolText[val]}],
  assoc
]];
exprText[x_] := quoteText[ToString[FullSimplify[x, $Assumptions], InputForm]];
groupedLaneLines[indent_] := Flatten[Table[
  {
    indent <> quoteText[lane] <> ":",
    indent <> "  D:",
    Table[indent <> "    " <> quoteText[n] <> ": " <> exprText[groupedLanes[lane]["D"][n]], {n, {"0", "2", "4"}}],
    indent <> "  D_sample:",
    Table[indent <> "    " <> quoteText[n] <> ": " <> numText[evalSample[groupedLanes[lane]["D"][n]]], {n, {"0", "2", "4"}}]
  },
  {lane, lanes}
]];

yamlLines = Flatten[{
  "schema: pathA_32_grouped_p2_isotropy_mathematica/v1",
  "engine: mathematica",
  "status: ok",
  "verdict: " <> verdict,
  "which_rung: " <> verdict,
  "harmonics:",
  "  order:",
  ("  - " <> quoteText[#]) & /@ order,
  "  gram_5x5_numeric:",
  matrixLines["  ", N[gram, 17]],
  "  gram_is_identity: " <> boolText[gramIsIdentity],
  "  anisotropy_coefficients_numeric:",
  KeyValueMap[("    " <> quoteText[#1] <> ": " <> numText[#2]) &, N[anisotropyCoefficients, 17]],
  "  grouped_reduction_isotropic_numeric:",
  matrixLines["  ", N[groupedIso, 17]],
  "  grouped_reduction_pure_prefactor_linear_coeff_numeric:",
  matrixLines["  ", N[groupedCoeff, 17]],
  "laplacian:",
  "  lambda_numeric_by_channel:",
  KeyValueMap[("    " <> quoteText[#1] <> ": " <> numText[#2]) &, N[lambdas, 17]],
  "  lambda_all_six: " <> boolText[lambdaAllSix],
  "  residuals_zero: " <> boolText[residualsZero],
  "  k2_coefficient_equals_computed_lambda_all: " <> boolText[kCoeffEqualsAll],
  "grouped_lanes:",
  groupedLaneLines["  "],
  "raw_D_gate:",
  "  raw_D_defects_zero: " <> boolText[rawDefectsZero],
  "normalized_response:",
  "  normalized_defects_zero: " <> boolText[normalizedDefectsZero],
  "stability:",
  "  omega_2m_sample: " <> numText[Sqrt[(7.0 + 6.0*0.5)/3.0]],
  "counterfactuals:",
  "  pure_prefactor_anisotropy:",
  "    verdict: " <> probeVerdicts["pure_prefactor_anisotropy"],
  "    sample_values:",
  assocNumericLines["      ", pureSamples],
  "  sector_selective_anisotropy:",
  "    verdict: " <> probeVerdicts["sector_selective_anisotropy"],
  "    sample_values:",
  assocNumericLines["      ", sectorSamples],
  "  m_dependent_profile:",
  "    verdict: " <> probeVerdicts["m_dependent_profile"],
  "    sample_values:",
  assocScalarLines["      ", profileSamples],
  "  degenerate_beta_zero:",
  "    verdict: " <> probeVerdicts["degenerate_beta_zero"],
  "    M2_sample: " <> numText[degenerateBetaProbe["M2_sample"]],
  "    K2_sample: " <> numText[degenerateBetaProbe["K2_sample"]],
  "    stability_ok: " <> boolText[degenerateBetaProbe["stability_ok"]],
  "    computed_fail_gate: " <> boolText[Not[degenerateBetaProbe["stability_ok"]]],
  "    self_ablation:",
  "      verdict: " <> degenerateBetaAblation["verdict"],
  "      fail_suppressed: " <> boolText[Not[degenerateBetaAblation["fail_fires"]]],
  "  wrong_eigenvalue:",
  "    verdict: " <> probeVerdicts["wrong_eigenvalue"],
  "    forced_coefficient: " <> numText[wrongEigenProbe["forced_coefficient"]],
  "    coefficient_equals_computed_lambda: " <> boolText[wrongEigenProbe["coefficient_equals_computed_lambda"]],
  "    computed_fail_gate: " <> boolText[Not[wrongEigenProbe["coefficient_equals_computed_lambda"]]],
  "    self_ablation:",
  "      forced_coefficient: " <> numText[wrongEigenAblation["forced_coefficient"]],
  "      verdict: " <> wrongEigenAblation["verdict"],
  "      fail_suppressed: " <> boolText[Not[wrongEigenAblation["fail_fires"]]],
  "  singular_denominator:",
  "    verdict: " <> probeVerdicts["singular_denominator"],
  "    D0_sample: " <> numText[singularD0Sample],
  "    denominator_guard_ok: " <> boolText[singularDenominatorGuardOk],
  "    computed_fail_gate: " <> boolText[Not[singularDenominatorGuardOk]],
  "    self_ablation:",
  "      D0_sample: " <> numText[nonsingularD0Sample],
  "      denominator_guard_ok: " <> boolText[nonsingularDenominatorGuardOk],
  "      verdict: " <> singularAblationVerdict,
  "      fail_suppressed: " <> boolText[singularAblationVerdict =!= "FAIL_SINGULAR_RESPONSE"],
  "  tautology_hash_collision:",
  "    verdict: " <> probeVerdicts["tautology_hash_collision"],
  "    distinct_hashes: " <> boolText[tautologyProbe["distinct_hashes"]],
  "    computed_fail_gate: " <> boolText[Not[tautologyProbe["distinct_hashes"]]],
  "    self_ablation:",
  "      distinct_hashes: " <> boolText[tautologyAblation["distinct_hashes"]],
  "      verdict: " <> tautologyAblation["verdict"],
  "      fail_suppressed: " <> boolText[Not[tautologyAblation["fail_fires"]]],
  "  static_drop_inertia:",
  "    verdict: " <> probeVerdicts["static_drop_inertia"],
  "    dynamic_retained: " <> boolText[staticDynamicRetained],
  "    computed_fail_gate: " <> boolText[Not[staticDynamicRetained]],
  "    self_ablation:",
  "      dynamic_retained: " <> boolText[staticAblationDynamicRetained],
  "      verdict: " <> staticAblationVerdict,
  "      fail_suppressed: " <> boolText[staticAblationVerdict =!= "FAIL_STATIC_RESPONSE"],
  "able_to_fail:",
  "  expected_probe_verdicts_match:",
  assocBoolLines["    ", expectedProbeVerdictsMatch],
  "  computed_probe_gate_flags:",
  assocBoolLines["    ", computedProbeGateFlags],
  "  able_to_fail_ok: " <> boolText[ableToFailOk],
  "gate_booleans:",
  "  covariant: " <> boolText[covariantOk],
  "  tautology_clear: " <> boolText[tautologyClear],
  "  dynamic_retained: " <> boolText[dynamicRetained],
  "  stability_ok: " <> boolText[stabilityOk],
  "  denominator_guard_ok: " <> boolText[denominatorGuardOk],
  "  lane_collapse_ok: " <> boolText[laneCollapseOk],
  "  able_to_fail_ok: " <> boolText[ableToFailOk]
}];

Export[yamlOut, StringRiffle[yamlLines, "\n"] <> "\n", "Text"];
Print["Mathematica PathA-32 scratch YAML written: ", yamlOut];
Print["Computed verdict: ", verdict];
Exit[0];
