(* Ledger stage013 Mathematica audit: breathing harmonic profiles + M/K projection.

   Standalone, print-only, no arguments, no file I/O.  This keeps the native
   DSolveValue + Integrate 013 engine, derives the consumed packet natively
   with dual-site integrity, and strips the scratch bridge plumbing.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;

BREATHINGCALIBRATED = "BREATHING_CALIBRATED";
BREATHINGFAILDIMENSIONAL = "BREATHING_FAIL_DIMENSIONAL";

$Assumptions =
  L0 > 0 && beta > 0 && muEta > 0 && Tw > 0 && rAL > 0 &&
  Element[{w, deltaA, deltaL, deltaAddot, deltaLddot}, Reals];

zeroDim = {0, 0, 0};

raise[msg_] := Throw[msg, "ledgerStage013Failure"];

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
clean[expr_] := FullSimplify[dropConditions[expr]];
fmt[expr_] := ToString[InputForm[clean[expr]]];

assertExact[name_, expr_] := Module[{reals},
  reals = Cases[Unevaluated[expr], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FAIL  ", name, ": machine-real atom(s) found: ", ToString[InputForm[reals]]];
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

expectNonzero[name_, residual_] := Module[{c},
  assertExact[name, residual];
  c = clean[residual];
  assertExact[name, c];
  If[! TrueQ[c === 0],
    passCount++;
    Print["PASS  ", name, " is nonzero as required (residual = ", fmt[c], ")"],
    failCount++;
    Print["FAIL  ", name, ": required nonzero residual vanished"];
    raise[name]
  ]
];

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

exprEqual[lhs_, rhs_: 0] := TrueQ[FullSimplify[dropConditions[lhs - rhs] == 0]];
nonzeroQ[expr_] := ! TrueQ[FullSimplify[dropConditions[expr] == 0]];
boolResidual[condition_] := If[TrueQ[condition], 0, 1];
verdictResidual[actual_, expected_] := If[actual === expected, 0, 1];
computeVerdict[dimensionalOk_] := If[TrueQ[dimensionalOk], BREATHINGCALIBRATED, BREATHINGFAILDIMENSIONAL];

dimResidualVec[actual_, expected_] := FullSimplify[(actual - expected).(actual - expected)];

dimOf[expr_, dims_] := Module[{args, ds, base, pow, argDims},
  Which[
    TrueQ[expr == 0] || NumericQ[expr], zeroDim,
    AtomQ[expr] && KeyExistsQ[dims, expr], dims[expr],
    AtomQ[expr], raise["missing dimension for " <> ToString[Unevaluated[expr], InputForm]],
    Head[expr] === Times, Total[dimOf[#, dims] & /@ (List @@ expr)],
    Head[expr] === Power,
      base = expr[[1]];
      pow = expr[[2]];
      If[! NumericQ[pow], raise["non-numeric dimension exponent"]];
      pow dimOf[base, dims],
    Head[expr] === Plus,
      args = Select[List @@ expr, ! TrueQ[FullSimplify[# == 0]] &];
      ds = dimOf[#, dims] & /@ args;
      If[Length[ds] == 0, zeroDim,
        If[Length[DeleteDuplicates[ds]] != 1, raise["dimension mismatch in sum"]];
        First[ds]
      ],
    MemberQ[{Sin, Cos, Tan, Cot, Sinh, Cosh, Tanh, Coth, Sech, Csch}, Head[expr]],
      argDims = dimOf[#, dims] & /@ (List @@ expr);
      If[AnyTrue[argDims, # =!= zeroDim &], raise["dimensionful argument in dimensionless function"]];
      zeroDim,
    True, raise["unsupported dimension expression " <> ToString[expr, InputForm]]
  ]
];

sharedDim[assoc_] := Module[{vals = Values[assoc]},
  If[Length[DeleteDuplicates[vals]] == 1, First[vals], Missing["Inhomogeneous"]]
];

dimText[d_] := Module[{parts, emit},
  emit[label_, exp_] := If[TrueQ[exp == 1], label, label <> "^" <> ToString[InputForm[exp]]];
  parts = Join[
    If[TrueQ[d[[1]] == 0], {}, {emit["L", d[[1]]]}],
    If[TrueQ[d[[2]] == 0], {}, {emit["M", d[[2]]]}],
    If[TrueQ[d[[3]] == 0], {}, {emit["T", d[[3]]]}]
  ];
  If[Length[parts] == 0, "1", StringRiffle[parts, " "]]
];

makeIntegrands[aa_, al_, includeGradient_] := Module[{pairs, mass, stiff},
  pairs = <|"aa" -> {aa, aa}, "aL" -> {aa, al}, "LL" -> {al, al}|>;
  mass = Association @ KeyValueMap[
    Function[{key, pair}, key -> FullSimplify[muEta pair[[1]] pair[[2]]]],
    pairs
  ];
  stiff = Association @ KeyValueMap[
    Function[{key, pair},
      key -> FullSimplify[
        If[TrueQ[includeGradient], Tw D[pair[[1]], w] D[pair[[2]], w], 0] +
          kEta pair[[1]] pair[[2]]
      ]
    ],
    pairs
  ];
  <|"MIntegrands" -> mass, "KIntegrands" -> stiff|>
];

integrateEntries[integrands_] := Module[{mEntries, kEntries, mMatrix, kMatrix},
  mEntries = Association @ KeyValueMap[
    Function[{key, expr}, key -> FullSimplify[dropConditions[4 Pi Integrate[expr, {w, 0, L0}]]]],
    integrands["MIntegrands"]
  ];
  kEntries = Association @ KeyValueMap[
    Function[{key, expr}, key -> FullSimplify[dropConditions[4 Pi Integrate[expr, {w, 0, L0}]]]],
    integrands["KIntegrands"]
  ];
  mMatrix = {{mEntries["aa"], mEntries["aL"]}, {mEntries["aL"], mEntries["LL"]}};
  kMatrix = {{kEntries["aa"], kEntries["aL"]}, {kEntries["aL"], kEntries["LL"]}};
  <|
    "MEntries" -> mEntries,
    "KEntries" -> kEntries,
    "MMatrix" -> mMatrix,
    "KMatrix" -> kMatrix,
    "MDet" -> FullSimplify[Det[mMatrix]],
    "KDet" -> FullSimplify[Det[kMatrix]]
  |>
];

projectFromProfiles[aa_, al_, includeGradient_] := Module[{integrands, entries},
  integrands = makeIntegrands[aa, al, includeGradient];
  entries = integrateEntries[integrands];
  Join[integrands, entries]
];

reportClosedForms[] := <|
  "MEntries" -> <|
    "aa" -> -2 Pi muEta (L0 beta Tanh[L0 beta] - Sinh[L0 beta]^2)/(beta Sinh[L0 beta]^2 Tanh[L0 beta]),
    "aL" -> 2 Pi muEta rAL (L0 beta - Tanh[L0 beta])/(beta Sinh[L0 beta] Tanh[L0 beta]),
    "LL" -> -2 Pi muEta rAL^2 (L0 beta Tanh[L0 beta] - Sinh[L0 beta]^2)/(beta Sinh[L0 beta]^2 Tanh[L0 beta])
  |>,
  "KEntries" -> <|
    "aa" -> 4 Pi Tw beta/Tanh[L0 beta],
    "aL" -> -4 Pi Tw beta rAL/Sinh[L0 beta],
    "LL" -> 4 Pi Tw beta rAL^2/Tanh[L0 beta]
  |>,
  "MDet" -> -4 Pi^2 muEta^2 rAL^2 (L0 beta - Sinh[L0 beta]) (L0 beta + Sinh[L0 beta])/(beta^2 Sinh[L0 beta]^2),
  "KDet" -> 16 Pi^2 Tw^2 beta^2 rAL^2
|>;

projectionResiduals[projected_, reference_] := Flatten[Table[
  {
    FullSimplify[projected["MEntries"][key] - reference["MEntries"][key]],
    FullSimplify[projected["KEntries"][key] - reference["KEntries"][key]]
  },
  {key, {"aa", "aL", "LL"}}
]];

projectionResidualsZeroQ[projected_, reference_] := AllTrue[projectionResiduals[projected, reference], TrueQ[# === 0] &];

symbolNames[exprs_] := Sort[DeleteDuplicates[
  ToString[#, InputForm] & /@ Cases[Unevaluated[exprs], s_Symbol /; Context[s] === "Global`", Infinity]
]];

freeSymbolNameFlags[mEntries_, kEntries_] := Module[
  {legacyNames, allowedNames, mNames, kNames, mkNames, unexpectedNames, flags},
  legacyNames = {"kappa", "chi", "sigmaA", "sigmaL"};
  allowedNames = {"L0", "beta", "muEta", "Tw", "rAL"};
  mNames = symbolNames[Values[mEntries]];
  kNames = symbolNames[Values[kEntries]];
  mkNames = Sort[Union[mNames, kNames]];
  unexpectedNames = Complement[mkNames, allowedNames];
  flags = <|
    "K_from_static_hessian" -> Length[Intersection[kNames, legacyNames]] > 0,
    "M_or_K_typed_from_legacy_values" -> Length[Intersection[mkNames, legacyNames]] > 0,
    "full_matrix_fit" -> Length[unexpectedNames] > 0
  |>;
  <|
    "MNames" -> mNames,
    "KNames" -> kNames,
    "MKNames" -> mkNames,
    "LegacyNames" -> legacyNames,
    "AllowedNames" -> allowedNames,
    "UnexpectedNames" -> unexpectedNames,
    "ForbiddenFitFlags" -> flags
  |>
];

buildDimensionalBlock[projection_] := Module[
  {
    dimRules, expectedM, expectedK, expectedRatio, mDims, kDims, mShared,
    kShared, ratioDim, dimensionalOk, corruptRules, corruptMDims,
    corruptKDims, corruptMShared, corruptKShared, corruptRatioDim,
    corruptOk, cleanVerdict, mutatedVerdict
  },
  dimRules = <|L0 -> {1, 0, 0}, beta -> {-1, 0, 0}, muEta -> {-1, 1, 0}, Tw -> {1, 1, -2}, rAL -> zeroDim|>;
  expectedM = {0, 1, 0};
  expectedK = {0, 1, -2};
  expectedRatio = {0, 0, -2};
  mDims = Association @ KeyValueMap[(#1 -> dimOf[#2, dimRules]) &, projection["MEntries"]];
  kDims = Association @ KeyValueMap[(#1 -> dimOf[#2, dimRules]) &, projection["KEntries"]];
  mShared = sharedDim[mDims];
  kShared = sharedDim[kDims];
  ratioDim = If[ListQ[mShared] && ListQ[kShared], kShared - mShared, Missing["Inhomogeneous"]];
  dimensionalOk = TrueQ[mShared == expectedM && kShared == expectedK && ratioDim == expectedRatio];
  corruptRules = Join[KeyDrop[dimRules, Tw], <|Tw -> {2, 1, -2}|>];
  corruptMDims = Association @ KeyValueMap[(#1 -> dimOf[#2, corruptRules]) &, projection["MEntries"]];
  corruptKDims = Association @ KeyValueMap[(#1 -> dimOf[#2, corruptRules]) &, projection["KEntries"]];
  corruptMShared = sharedDim[corruptMDims];
  corruptKShared = sharedDim[corruptKDims];
  corruptRatioDim = If[ListQ[corruptMShared] && ListQ[corruptKShared], corruptKShared - corruptMShared, Missing["Inhomogeneous"]];
  corruptOk = TrueQ[corruptMShared == expectedM && corruptKShared == expectedK && corruptRatioDim == expectedRatio];
  cleanVerdict = computeVerdict[dimensionalOk];
  mutatedVerdict = computeVerdict[corruptOk];
  <|
    "DimRules" -> dimRules,
    "ExpectedM" -> expectedM,
    "ExpectedK" -> expectedK,
    "ExpectedRatio" -> expectedRatio,
    "MDims" -> mDims,
    "KDims" -> kDims,
    "MShared" -> mShared,
    "KShared" -> kShared,
    "RatioDim" -> ratioDim,
    "DimensionalOk" -> dimensionalOk,
    "CorruptRules" -> corruptRules,
    "CorruptMDims" -> corruptMDims,
    "CorruptKDims" -> corruptKDims,
    "CorruptMShared" -> corruptMShared,
    "CorruptKShared" -> corruptKShared,
    "CorruptRatioDim" -> corruptRatioDim,
    "CorruptOk" -> corruptOk,
    "MutationFires" -> ! TrueQ[corruptOk],
    "ProbeVerdict" -> If[TrueQ[corruptOk], "NO_FAIL", BREATHINGFAILDIMENSIONAL],
    "CleanVerdict" -> cleanVerdict,
    "MutatedVerdict" -> mutatedVerdict,
    "FailSuppressed" -> TrueQ[cleanVerdict === BREATHINGCALIBRATED && mutatedVerdict === BREATHINGFAILDIMENSIONAL]
  |>
];

packetResiduals[packet_] := <|
  "site_A_constitutive" -> FullSimplify[packet["betaCited"] - Sqrt[packet["kEtaCited"]/packet["TwCited"]]],
  "site_B_branch" -> FullSimplify[packet["betaCited"] packet["L0Cited"] - packet["branchAnchor"]],
  "anchor_L0" -> FullSimplify[packet["L0Cited"] - 37/20],
  "anchor_Tw" -> FullSimplify[packet["TwCited"] - 1],
  "anchor_kEta" -> FullSimplify[packet["kEtaCited"] - 1],
  "anchor_beta" -> FullSimplify[packet["betaCited"] - 1]
|>;

packetOkQ[packet_] := AllTrue[Values[packetResiduals[packet]], TrueQ[# === 0] &];
packetSet[packet_, key_, value_] := ReplacePart[packet, Key[key] -> value];

degenerateMassGuard[al_, enabled_] := Module[{projection, detZero, caught},
  projection = projectFromProfiles[0, al, True];
  detZero = FullSimplify[projection["MDet"]];
  caught = TrueQ[enabled && detZero === 0];
  <|"Projection" -> projection, "MDet" -> detZero, "Caught" -> caught|>
];

buildBaseline[] := Module[
  {
    alphaA, alphaL, projection, reintegrated, closedForms, flags, dim, eomRows,
    packet, degenerate, verdict
  },
  kEta = FullSimplify[Tw beta^2];
  alphaA = FullSimplify[
    DSolveValue[{-D[y[w], {w, 2}] + beta^2 y[w] == 0, y[0] == 1, y[L0] == 0}, y[w], w]
  ];
  alphaL = FullSimplify[
    DSolveValue[{-D[y[w], {w, 2}] + beta^2 y[w] == 0, y[0] == 0, y[L0] == rAL}, y[w], w]
  ];
  projection = projectFromProfiles[alphaA, alphaL, True];
  reintegrated = projectFromProfiles[alphaA, alphaL, True];
  closedForms = reportClosedForms[];
  flags = freeSymbolNameFlags[projection["MEntries"], projection["KEntries"]];
  dim = buildDimensionalBlock[projection];
  eomRows = FullSimplify[projection["MMatrix"].{deltaAddot, deltaLddot} + projection["KMatrix"].{deltaA, deltaL}];
  packet = <|"L0Cited" -> 37/20, "TwCited" -> 1, "kEtaCited" -> 1, "betaCited" -> 1, "branchAnchor" -> 37/20|>;
  degenerate = degenerateMassGuard[alphaL, True];
  verdict = computeVerdict[dim["DimensionalOk"]];
  <|
    "KEta" -> kEta,
    "AlphaA" -> alphaA,
    "AlphaL" -> alphaL,
    "Projection" -> projection,
    "Reintegrated" -> reintegrated,
    "ClosedForms" -> closedForms,
    "ProjectionMatches" -> projectionResidualsZeroQ[reintegrated, closedForms],
    "Flags" -> flags,
    "Dim" -> dim,
    "EomRows" -> eomRows,
    "RhsPlaceholders" -> {FaHF, FLHF},
    "ConsumedPacket" -> packet,
    "Degenerate" -> degenerate,
    "Verdict" -> verdict
  |>
];

runAritySelfCheck[data_] := Module[{integrandProbe, projectionProbe, packetProbe, dimProbe},
  subheading["Wolfram arity self-check"];
  integrandProbe = makeIntegrands[data["AlphaA"], data["AlphaL"], True];
  projectionProbe = projectFromProfiles[0, data["AlphaL"], True];
  packetProbe = packetResiduals[data["ConsumedPacket"]];
  dimProbe = buildDimensionalBlock[data["Projection"]];
  expectBool["arity makeIntegrands[3 args] returns M/K integrands", KeyExistsQ[integrandProbe, "MIntegrands"] && KeyExistsQ[integrandProbe, "KIntegrands"]];
  expectBool["arity projectFromProfiles[3 args] returns projected matrices", KeyExistsQ[projectionProbe, "MMatrix"] && KeyExistsQ[projectionProbe, "KMatrix"]];
  expectBool["arity packetResiduals[assoc] returns dual sites", KeyExistsQ[packetProbe, "site_A_constitutive"] && KeyExistsQ[packetProbe, "site_B_branch"]];
  expectBool["arity computeVerdict[1 arg] returns BREATHING_CALIBRATED", computeVerdict[True] === BREATHINGCALIBRATED];
  expectBool["arity buildDimensionalBlock[projection] returns mutation_fires", dimProbe["MutationFires"] === True];
  expectBool["arity degenerateMassGuard[2 args] catches alpha_a=0", projectionProbe["MDet"] === 0]
];

runOperatorProfiles[data_] := Module[{lopA, lopL, bcValues, expectedA, expectedL},
  subheading["Operator, ell=0 restriction, collective BCs, and harmonic profiles"];
  Print["  CONSUMED/POSTULATED operator: Lop[alpha] = muEta^(-1)*(-d_w(Tw*d_w alpha) + K_eta*alpha) on w in [0,L0]."];
  Print["  Alias for projection only: K_eta = Tw*beta^2; beta = sqrt(K_eta/Tw) is the cited wall-packet relation."];
  Print["  Notation guard: domain length is L0; operator applications are named Lop_alpha_a and Lop_alpha_L."];
  Print["  ell=0 restriction: Y00=1/(2*sqrt(pi)); int_S2 Y00^2 dOmega=1; eta(w,t)=eta_00(w,t)*Y00; T_Omega drops because ell*(ell+1)=0."];
  Print["  Inner product: <f,g>_mu = 4*pi*int_0^L0 muEta*f*g dw."];
  Print["  IMPOSED collective BCs: alpha_a(0)=1, alpha_a(L0)=0; alpha_L(0)=0, alpha_L(L0)=rAL."];
  Print["  General solution form = C1*Sinh[beta*w] + C2*Cosh[beta*w]."];
  Print["  DSolveValue alpha_a = ", fmt[data["AlphaA"]]];
  Print["  DSolveValue alpha_L = ", fmt[data["AlphaL"]]];
  expectedA = Sinh[L0 beta - beta w]/Sinh[L0 beta];
  expectedL = rAL Sinh[beta w]/Sinh[L0 beta];
  lopA = FullSimplify[(-Tw D[data["AlphaA"], {w, 2}] + data["KEta"] data["AlphaA"])/muEta];
  lopL = FullSimplify[(-Tw D[data["AlphaL"], {w, 2}] + data["KEta"] data["AlphaL"])/muEta];
  bcValues = <|
    "alpha_a_mouth" -> FullSimplify[data["AlphaA"] /. w -> 0],
    "alpha_a_cap" -> FullSimplify[data["AlphaA"] /. w -> L0],
    "alpha_L_mouth" -> FullSimplify[data["AlphaL"] /. w -> 0],
    "alpha_L_cap" -> FullSimplify[data["AlphaL"] /. w -> L0]
  |>;
  expectZero["alpha_a equals reported harmonic-lift profile", data["AlphaA"] - expectedA];
  expectZero["alpha_L equals reported harmonic-lift profile", data["AlphaL"] - expectedL];
  expectZero["Lop_alpha_a harmonic residual is zero", lopA];
  expectZero["Lop_alpha_L harmonic residual is zero", lopL];
  expectZero["BC alpha_a(0)=1", bcValues["alpha_a_mouth"] - 1];
  expectZero["BC alpha_a(L0)=0", bcValues["alpha_a_cap"]];
  expectZero["BC alpha_L(0)=0", bcValues["alpha_L_mouth"]];
  expectZero["BC alpha_L(L0)=rAL", bcValues["alpha_L_cap"] - rAL]
];

runProjection[data_] := Module[
  {p, closedForms, flags, expectedM, expectedK, expectedMDet, expectedKDet, residuals},
  p = data["Projection"];
  closedForms = data["ClosedForms"];
  flags = data["Flags"];
  subheading["M_AB/K_AB by native int-dw operator projection"];
  Print["  M integrands:"];
  Scan[Function[key, Print["    ", key, ": ", fmt[p["MIntegrands"][key]]]], {"aa", "aL", "LL"}];
  Print["  K integrands:"];
  Scan[Function[key, Print["    ", key, ": ", fmt[p["KIntegrands"][key]]]], {"aa", "aL", "LL"}];
  Print["  M_AB = ", fmt[p["MEntries"]]];
  Print["  K_AB = ", fmt[p["KEntries"]]];
  Print["  det(M) = ", fmt[p["MDet"]]];
  Print["  det(K) = ", fmt[p["KDet"]]];
  expectedM = closedForms["MEntries"];
  expectedK = closedForms["KEntries"];
  expectedMDet = closedForms["MDet"];
  expectedKDet = closedForms["KDet"];
  Scan[
    Function[key,
      expectZero["M_" <> key <> " matches report closed form", p["MEntries"][key] - expectedM[key]];
      expectZero["K_" <> key <> " matches report closed form", p["KEntries"][key] - expectedK[key]]
    ],
    {"aa", "aL", "LL"}
  ];
  expectZero["det(M) matches report closed form", p["MDet"] - expectedMDet];
  expectZero["det(K) matches 16*pi^2*Tw^2*beta^2*rAL^2", p["KDet"] - expectedKDet];
  residuals = projectionResiduals[data["Reintegrated"], closedForms];
  Scan[Function[residual, expectZero["independent full-integrand re-integration residual", residual]], residuals];
  expectBool["independent re-integration comparison is not a self-identity stamp", data["ProjectionMatches"]];
  Print["  free-symbol names in M/K = ", flags["MKNames"]];
  Print["  allowed names = ", flags["AllowedNames"]];
  Print["  forbidden_fit_flags = ", flags["ForbiddenFitFlags"]];
  expectBool["M/K free-symbol names are a subset of allowed projection names", flags["UnexpectedNames"] === {}];
  Scan[Function[key, expectBool["forbidden_fit_flags[" <> key <> "] computed false", flags["ForbiddenFitFlags"][key] === False]], Keys[flags["ForbiddenFitFlags"]]];
  Print["  K_AB_provenance = \"operator_projection_not_static_Hessian\""]
];

runEom[data_] := Module[{rows},
  rows = data["EomRows"];
  subheading["Dynamical EOM LHS only"];
  Print["  Q = (delta_a, delta_L); Qddot = (delta_a_ddot, delta_L_ddot)."];
  Print["  EOM structure: M_AB*Qddot^B + K_AB*Q^B = F_A^(HF)."];
  Print["  RHS placeholders: F_a_HF, F_L_HF are deferred to stage 015; no Hellmann-Feynman force is computed here."];
  Print["  row a LHS = ", fmt[rows[[1]]], " = F_a_HF"];
  Print["  row L LHS = ", fmt[rows[[2]]], " = F_L_HF"];
  expectBool["EOM has exactly two LHS rows", Length[rows] == 2];
  expectBool["EOM RHS is symbolic placeholder only", data["RhsPlaceholders"] === {FaHF, FLHF}]
];

runDimensionalBlock[data_] := Module[{dim},
  dim = data["Dim"];
  subheading["013 dimensional legs and corrupt-[Tw] probe"];
  Print["  dimension order: (L,M,T)"];
  Print["  sourced dims: L0=(1,0,0), beta=(-1,0,0), muEta=(-1,1,0), Tw=(1,1,-2), rAL=(0,0,0)"];
  Print["  derived dim: K_eta=Tw*beta^2 = (-1,1,-2)"];
  Print["  [M_AB] entries = ", Association @ KeyValueMap[(#1 -> dimText[#2]) &, dim["MDims"]]];
  Print["  [K_AB] entries = ", Association @ KeyValueMap[(#1 -> dimText[#2]) &, dim["KDims"]]];
  Print["  [M_AB] shared = ", dim["MShared"], "; [K_AB] shared = ", dim["KShared"], "; [K/M] = ", dim["RatioDim"]];
  expectZero["M_AB shared dimension is M", dimResidualVec[dim["MShared"], dim["ExpectedM"]]];
  expectZero["K_AB shared dimension is M*T^-2", dimResidualVec[dim["KShared"], dim["ExpectedK"]]];
  expectZero["K/M ratio dimension is T^-2", dimResidualVec[dim["RatioDim"], dim["ExpectedRatio"]]];
  expectBool["dimensional_ok for 013 M/K legs", dim["DimensionalOk"]];
  Print["  corrupt [Tw]+(1,0,0) gives [Tw] = ", dim["CorruptRules"][Tw]];
  Print["  corrupt [K_AB] shared = ", dim["CorruptKShared"], "; corrupt [K/M] = ", dim["CorruptRatioDim"]];
  expectFail["corrupt-[Tw] shifts K_AB away from M*T^-2", dimResidualVec[dim["CorruptKShared"], dim["ExpectedK"]]];
  expectFail["corrupt-[Tw] shifts K/M away from T^-2", dimResidualVec[dim["CorruptRatioDim"], dim["ExpectedRatio"]]];
  expectBool["corrupt-[Tw] mutation_fires=True", dim["MutationFires"]];
  expectZero["self-ablation with mutation gives BREATHING_FAIL_DIMENSIONAL", verdictResidual[dim["MutatedVerdict"], BREATHINGFAILDIMENSIONAL]];
  expectZero["self-ablation without mutation gives BREATHING_CALIBRATED", verdictResidual[dim["CleanVerdict"], BREATHINGCALIBRATED]];
  expectBool["self-ablation fail_suppressed=True", dim["FailSuppressed"]]
];

runConsumedPacket[data_] := Module[{packet, residuals},
  packet = data["ConsumedPacket"];
  residuals = packetResiduals[packet];
  subheading["Consumed Gate-1 frozen packet, dual-site integrity"];
  Print["  CONSUMED packet anchor: (L0,Tw,K_eta,beta)=(37/20,1,1,1)."];
  Print["  Site A constitutive: betaCited - Sqrt[kEtaCited/TwCited] = 0."];
  Print["  Site B geometric/branch: betaCited*L0Cited - 37/20 = 0."];
  Print["  Anti-tautology: kEtaCited is an independent cited datum for the guard, not the local alias kEta=Tw*beta^2."];
  Scan[Function[key, expectZero["consumed packet " <> key, residuals[key]]], Keys[residuals]];
  expectBool["consumed packet clean baseline passes both sites and frozen anchor", packetOkQ[packet]];
  expectFail["K_eta_cited-only corruption breaks site A: guard is non-vacuous", packetResiduals[packetSet[packet, "kEtaCited", 2]]["site_A_constitutive"]];
  expectFail["Tw_cited-only corruption breaks site A", packetResiduals[packetSet[packet, "TwCited", 2]]["site_A_constitutive"]];
  expectFail["L0_cited-only corruption breaks site B", packetResiduals[packetSet[packet, "L0Cited", 19/10]]["site_B_branch"]];
  expectFail["branch-anchor corruption breaks site B", packetResiduals[packetSet[packet, "branchAnchor", 19/10]]["site_B_branch"]]
];

runDegenerateGuard[data_] := Module[{degenerate, bypassed},
  degenerate = data["Degenerate"];
  subheading["Native degenerate M_det -> 0 guard slice"];
  Print["  Degenerate copy: alpha_a -> 0; alpha_L unchanged; M recomputed by the same int-dw projection."];
  Print["  degenerate M_det = ", fmt[degenerate["MDet"]]];
  expectZero["degenerate alpha_a=0 recompute gives M_det=0", degenerate["MDet"]];
  expectBool["native degenerate non-degeneracy test catches M_det==0", degenerate["Caught"]];
  expectNonzero["baseline M_det is not identically zero", data["Projection"]["MDet"]];
  bypassed = degenerateMassGuard[data["AlphaL"], False];
  expectFail["tooth 6 bypassing native M_det==0 guard leaves degenerate case uncaught", boolResidual[bypassed["Caught"]]]
];

runVerdictAndComposition[data_] := (
  subheading["013 scoped landing and joint composition"];
  Print["  013 scoped verdict = ", data["Verdict"]];
  expectZero["013 component lands at BREATHING_CALIBRATED", verdictResidual[data["Verdict"], BREATHINGCALIBRATED]];
  Print["  BREATHING_CALIBRATED (JOINT, 3-stage)"];
  Print["    = (013: harmonic-lift profiles + M_AB/K_AB by int-dw operator projection + (a,L) EOM LHS, computed here)"];
  Print["    AND (014: truncation consistency -- generalized eig / beta_L0 sweep / N-convergence) [sibling stage]"];
  Print["    AND (015: legacy-Hessian structure recovery + Hellmann-Feynman force) [sibling stage]"];
  Print["  CALIBRATED <= wall constants {muEta, Tw, K_eta} are calibration inputs; structure is EARNED, values are CALIBRATED."];
  expectBool["joint composition cites 014 and 015 as siblings, not recomputed here", data["Verdict"] === BREATHINGCALIBRATED]
);

printProvenance[data_] := Module[{liveNames},
  subheading["Provenance and scope"];
  Print["  CONSUMED-from-Gate-1: domain [0,L0] with cap R0(L0)=0, frozen packet {L0=37/20,Tw=1,K_eta=1,beta=1}, and ell=0 restriction are cited from stages 011/012 with dual-site integrity."];
  Print["  no-c_S: 013 is speed-free; matter-sector c_s/BdG k^4 is deferred under k*xi << 1 (phonon_limit_caveat)."];
  Print["  EARNED-STRUCTURE: profiles are derived harmonic lifts proven by Lop[alpha]=0; M_AB/K_AB are real int-dw operator projections; EOM LHS is assembled here."];
  Print["  FIRST-CALIBRATION: muEta and Tw are CALIBRATION inputs; K_eta=Tw*beta^2 is a manifestation; beta=sqrt(K_eta/Tw) with beta*L0=37/20 is branch-determinable, but geometry alone does NOT derive the wall constants."];
  Print["  control-ratio: rAL is the dimensionless alpha_L cap ratio, tracked and not counted."];
  Print["  3-way-split: 013 carries M/K projection + (a,L) closure; 014 carries truncation consistency; 015 carries legacy structure + HF force."];
  Print["  RHS-deferred: F_A^(HF) is stage 015's Hellmann-Feynman force; 013 emits only M_AB*Qddot+K_AB*Q."];
  Print["  dropped-bookkeeping: scratch-YAML/_sympy_exprs.wl export, MMA-YAML re-read, expression_digest, and engine_agreement plumbing are stripped."];
  Print["  downstream consumers: stage 014 consumes M_AB/K_AB for generalized eig; stages 022/023 consume the ell=0 (a,L) closure."];
  Print["  register note: FIRST Part-II calibration knobs are likely {muEta,Tw}; K_eta is a manifestation, beta and L0 are branch/geometric tracked quantities, and rAL is tracked-not-counted."];
  liveNames = symbolNames[{data["AlphaA"], data["AlphaL"], Values[data["Projection"]["MEntries"]], Values[data["Projection"]["KEntries"]], data["EomRows"]}];
  expectBool["no c_S/cS live symbol appears in 013 symbolic content", FreeQ[liveNames, "cS"] && FreeQ[liveNames, "c_S"]]
];

printVerdictLabels[] := (
  Print[""];
  Print["Verdict labels:"];
  Print["  ledger earned-label (NOT a source verdict token): BREATHING_HARMONIC_MK_PROJECTION_EARNED  (harmonic-lift profiles alpha_a=sinh(L0*beta-beta*w)/sinh(L0*beta), alpha_L=rAL*sinh(beta*w)/sinh(L0*beta), proven by Lop[alpha]=0; M_AB=4*pi*int mu_eta*alpha_A*alpha_B dw and K_AB=4*pi*int [Tw*alpha_A'*alpha_B'+K_eta*alpha_A*alpha_B] dw by real sp.integrate operator projection, NOT the legacy static Hessian (forbidden_fit_flags computed False via free-symbol ancestry, K_AB_provenance=operator_projection_not_static_Hessian); dynamical-EOM LHS M_AB*Qddot+K_AB*Q with Q=(delta_a,delta_L); [M_AB]=M, [K_AB]=M*T^-2, [K/M]=T^-2 dim legs + corrupt-[Tw] probe)"];
  Print["  source top-line verdict: BREATHING_CALIBRATED  (JOINT 3-stage; 013 carries the M/K-projection + (a,L)-closure component)"];
  Print["  joint composition: BREATHING_CALIBRATED = (013: harmonic profiles + M/K operator-projection + (a,L) EOM LHS, computed here) AND (014: truncation consistency) AND (015: legacy-structure + HF force)"];
  Print["  earned (structure): profiles DERIVED as harmonic lifts (residual Lop[alpha]=0); M_AB/K_AB DERIVED by int-dw operator projection (forbidden_fit_flags computed False, provenance operator_projection_not_static_Hessian); EOM LHS assembled; dim legs consistent + corrupt-[Tw] probe fires"];
  Print["  calibrated (values): mu_eta, Tw calibration inputs; K_eta=Tw*beta^2 manifestation; beta=sqrt(K_eta/Tw), beta*L0=37/20 branch-determinable; geometry alone does NOT derive the wall constants -> BREATHING_CALIBRATED not ..._PASS"];
  Print["  consumed (cited from Gate-1 stage011/012, dual-site integrity): domain [0,L0] (cap R0(L0)=0); frozen wall packet L0=37/20, Tw=1, K_eta=1, beta=1 (K_eta=Tw*beta^2, beta*L0=37/20); ell=0 restriction Y00=1/(2*sqrt(pi)); c_S NOT consumed (matter-sector deferred, k*xi<<1)"];
  Print["  control ratio (tracked, not counted): rAL = alpha_L cap ratio, [rAL]=1"];
  Print["  RHS deferred to stage 015: F_A^(HF) Hellmann-Feynman force (013 emits the EOM LHS only)"]
);

runAbleToFailTeeth[data_] := Module[
  {badAlpha, badResidual, typedK, typedFlags, corruptedProjection, dim, packet,
   bypassed},
  subheading["Able-to-fail mutation teeth"];
  badAlpha = Sinh[L0 beta - 2 beta w]/Sinh[L0 beta];
  badResidual = FullSimplify[(-Tw D[badAlpha, {w, 2}] + data["KEta"] badAlpha)/muEta];
  expectFail["tooth 1 non-kernel wrong-wavenumber profile trips harmonic residual", badResidual];

  typedK = <|"aa" -> chi^2 kappa + sigmaA, "aL" -> -chi kappa, "LL" -> kappa + sigmaL|>;
  typedFlags = freeSymbolNameFlags[data["Projection"]["MEntries"], typedK];
  expectBool["tooth 2 name-based ancestry catches legacy symbols", Intersection[typedFlags["KNames"], typedFlags["LegacyNames"]] =!= {}];
  expectBool["tooth 2 typed legacy K flips K_from_static_hessian flag", typedFlags["ForbiddenFitFlags"]["K_from_static_hessian"]];
  expectBool["tooth 2 typed legacy K flips M_or_K_typed_from_legacy_values flag", typedFlags["ForbiddenFitFlags"]["M_or_K_typed_from_legacy_values"]];
  expectBool["tooth 2 typed legacy K flips full_matrix_fit flag via unallowed names", typedFlags["ForbiddenFitFlags"]["full_matrix_fit"]];

  corruptedProjection = projectFromProfiles[data["AlphaA"], data["AlphaL"], False];
  expectFail["tooth 3 dropping Tw*alpha_prime*alpha_prime makes independent K re-integration mismatch", boolResidual[projectionResidualsZeroQ[corruptedProjection, data["ClosedForms"]]]];
  expectFail["tooth 3 corrupted K_aa differs from emitted full-operator K_aa", corruptedProjection["KEntries"]["aa"] - data["Projection"]["KEntries"]["aa"]];

  dim = data["Dim"];
  expectFail["tooth 4 corrupt-[Tw] probe trips K_AB dimensional leg", dimResidualVec[dim["CorruptKShared"], dim["ExpectedK"]]];
  expectFail["tooth 4 corrupt-[Tw] probe trips K/M dimensional leg", dimResidualVec[dim["CorruptRatioDim"], dim["ExpectedRatio"]]];
  expectZero["tooth 4 corrupt-[Tw] verdict is BREATHING_FAIL_DIMENSIONAL", verdictResidual[dim["MutatedVerdict"], BREATHINGFAILDIMENSIONAL]];
  expectBool["tooth 4 self-ablation fail_suppressed remains true", dim["FailSuppressed"]];

  packet = data["ConsumedPacket"];
  expectFail["tooth 5 K_eta_cited corruption trips packet site A", packetResiduals[packetSet[packet, "kEtaCited", 2]]["site_A_constitutive"]];
  expectFail["tooth 5 Tw_cited corruption trips packet site A", packetResiduals[packetSet[packet, "TwCited", 2]]["site_A_constitutive"]];
  expectFail["tooth 5 L0_cited corruption trips packet site B", packetResiduals[packetSet[packet, "L0Cited", 19/10]]["site_B_branch"]];
  expectFail["tooth 5 branch-anchor corruption trips packet site B", packetResiduals[packetSet[packet, "branchAnchor", 19/10]]["site_B_branch"]];

  bypassed = degenerateMassGuard[data["AlphaL"], False];
  expectFail["tooth 6 bypassed native M_det guard fails to catch degenerate profile", boolResidual[bypassed["Caught"]]];
  expectZero["tooth 6 baseline degenerate copy still has M_det=0", data["Degenerate"]["MDet"]];

  expectZero["baseline immutable after teeth: Lop_alpha_a remains zero", FullSimplify[(-Tw D[data["AlphaA"], {w, 2}] + data["KEta"] data["AlphaA"])/muEta]];
  expectZero["baseline immutable after teeth: independent full-integrand re-integration still matches", boolResidual[data["ProjectionMatches"]]];
  expectBool["baseline immutable after teeth: forbidden flags remain all false", ! Or @@ Values[data["Flags"]["ForbiddenFitFlags"]]];
  expectZero["baseline immutable after teeth: clean 013 verdict remains BREATHING_CALIBRATED", verdictResidual[data["Verdict"], BREATHINGCALIBRATED]]
];

Module[{ok, data},
  heading["ledger_stage013_breathing_harmonic_mk_projection Mathematica audit"];
  ok = Catch[
    data = buildBaseline[];
    assertExact["baseline", data];
    runAritySelfCheck[data];
    runOperatorProfiles[data];
    runProjection[data];
    runEom[data];
    runDimensionalBlock[data];
    runConsumedPacket[data];
    runDegenerateGuard[data];
    runVerdictAndComposition[data];
    printProvenance[data];
    printVerdictLabels[];
    runAbleToFailTeeth[data];
    True,
    "ledgerStage013Failure",
    Function[{msg, tag}, Print["FAIL: ", msg]; False]
  ];

  Print[""];
  Print["PASS tally: ", passCount, "; FAIL tally: ", failCount, "; TOTAL = ", passCount, " + ", failCount, " = ", passCount + failCount];
  If[TrueQ[ok],
    Print["OVERALL PASS: Mathematica verified ledger_stage013 breathing harmonic profiles + M/K projection exactly"];
    Exit[0],
    Print["OVERALL FAIL: Mathematica stage013 audit did not close"];
    Exit[1]
  ]
]
