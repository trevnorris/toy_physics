(* Ledger stage008 Mathematica audit: monopole/dipole return constraint spec.

   Standalone, print-only, no arguments, no imports, no exports.  This keeps
   the native pathA_28 Wolfram route: closed-form spherical Hankels, a
   Coefficient/ComplexExpand radiation scan, FullSimplify residuals, Integrate
   leakage anchors, and Which-based verdict machinery.
*)

ClearAll[
  raise, heading, subheading, cleanZero, assertExact, expectZero,
  expectBool, expectNonzero, expectFail, nonzeroQ, nonzeroAssertResidual,
  sameStringResidual, fmt, hankelComponents, radiatingPower, hankelData,
  dtnData, dtnAt, expectedPower, expectedCoeff, expectedKernel,
  dtnLadderResidual, kernelResidual, copiedWithCorruptEll2Kernel,
  sourceData, rawAmplitudeResidual, conditionWorks, residualWithCorruptReturn,
  runAritySelfCheck, runDtnLadder, runSourcesAndBookkeeping, stageAnchors,
  computeVerdict, runVerdict, runGuards, runFrozenAndConsumedStage005,
  dimResidualVec, expressionDim, kernelDimensionResidual, runDimensionalBlock, printScopeAndProvenance,
  printVerdictLabels, runAbleToFailTeeth, passCount, failCount
];

passCount = 0;
failCount = 0;

raise[msg_] := Throw[msg, "ledgerStage008Failure"];

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

cleanZero[expr_] := FullSimplify[expr] /. ConditionalExpression[0, _] -> 0;

assertExact[name_, expr_] := Module[{reals},
  reals = Cases[Unevaluated[expr], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FAIL  ", name, ": machine-real atom(s) found: ", ToString[InputForm[reals]]];
    raise[name]
  ]
];

expectZero[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FAIL  ", name, ": residual = ", ToString[InputForm[clean]]];
    raise[name]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

expectNonzero[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[! TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name, " is nonzero as required (residual = ", ToString[InputForm[clean]], ")"],
    failCount++;
    Print["FAIL  ", name, ": required nonzero residual vanished"];
    raise[name]
  ]
];

expectFail[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[! TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name, " produced required FAIL (residual = ", ToString[InputForm[clean]], ")"],
    failCount++;
    Print["FAIL  ", name, ": required mutation/ablation did not fire"];
    raise[name]
  ]
];

nonzeroQ[expr_] := ! TrueQ[FullSimplify[expr == 0]];
nonzeroAssertResidual[expr_] := If[nonzeroQ[expr], 0, 1];
sameStringResidual[actual_, expected_] := If[actual === expected, 0, 1];
fmt[expr_] := ToString[InputForm[Factor[Cancel[FullSimplify[expr]]]]];

$Assumptions =
  a > 0 && k > 0 && omega > 0 && cS > 0 && lambda > 0 &&
  muW > 0 && rho0 > 0 &&
  Element[{M0, D1, Q2, R0, R1, eta0, eta1, ellw, j0, E0}, Reals];

hankelComponents[ell_, zz_, corrupt_] := Switch[
  ell,
  0, {Sin[zz]/zz, -Cos[zz]/zz},
  1, {
    If[TrueQ[corrupt], Sin[zz]/zz^2, Sin[zz]/zz^2 - Cos[zz]/zz],
    -Cos[zz]/zz^2 - Sin[zz]/zz
  },
  2, {
    ((3/zz^3) - 1/zz) Sin[zz] - 3 Cos[zz]/zz^2,
    -((3/zz^3) - 1/zz) Cos[zz] - 3 Sin[zz]/zz^2
  },
  _, raise["unsupported ell"]
];

radiatingPower[y_, max_] := Module[{n, coeff, realCoeff, imagCoeff},
  For[n = 1, n <= max, n++,
    coeff = FullSimplify[Coefficient[Expand[y], k, n]];
    realCoeff = FullSimplify[ComplexExpand[Re[coeff]]];
    imagCoeff = FullSimplify[ComplexExpand[Im[coeff]]];
    If[! TrueQ[imagCoeff == 0] && TrueQ[realCoeff == 0],
      Return[{n, imagCoeff}]
    ];
  ];
  raise["no radiation-phase coefficient found in " <> ToString[InputForm[y]]]
];

hankelData[ell_, corrupt_] := Module[
  {zz, j, y, h, lambdaEll, lambdaSeries, admittanceSeries, staticValue,
   normalized, rp, pRaw, imagCoeff, kernel},
  zz = Unique["z"];
  {j, y} = hankelComponents[ell, zz, TrueQ[corrupt && ell == 1]];
  h = FullSimplify[j + I y];
  lambdaEll = FullSimplify[(k D[h, zz]/h) /. zz -> k a];
  lambdaSeries = Expand[Normal[Series[lambdaEll, {k, 0, 6}]]];
  admittanceSeries = Expand[Normal[Series[1/lambdaSeries, {k, 0, 6}]]];
  staticValue = FullSimplify[admittanceSeries /. k -> 0];
  normalized = Expand[Normal[Series[admittanceSeries/staticValue, {k, 0, 6}]]];
  rp = radiatingPower[normalized, 8];
  pRaw = rp[[1]];
  imagCoeff = FullSimplify[rp[[2]]];
  kernel = FullSimplify[I imagCoeff (omega/cS)^pRaw];
  <|
    "h" -> h,
    "LambdaSeries" -> lambdaSeries,
    "AdmittanceSeries" -> admittanceSeries,
    "StaticValue" -> staticValue,
    "YNormSeries" -> normalized,
    "pRaw" -> pRaw,
    "ImagCoeffK" -> imagCoeff,
    "RadiationKernel" -> kernel
  |>
];

dtnData[corrupt_] := AssociationThread[{0, 1, 2}, hankelData[#, corrupt] & /@ {0, 1, 2}];
dtnAt[dtn_, ell_] := Lookup[dtn, ell];

expectedPower[ell_] := Switch[ell, 0, 1, 1, 3, 2, 5, _, raise["bad expectedPower ell"]];
expectedCoeff[ell_] := Switch[ell, 0, a, 1, a^3/2, 2, a^5/27, _, raise["bad expectedCoeff ell"]];
expectedKernel[ell_] := Switch[
  ell,
  0, I a omega/cS,
  1, I a^3 omega^3/(2 cS^3),
  2, I a^5 omega^5/(27 cS^5),
  _, raise["bad expectedKernel ell"]
];

dtnLadderResidual[dtn_] := FullSimplify[
  Total[
    Table[
      (dtnAt[dtn, ell]["pRaw"] - expectedPower[ell])^2 +
        FullSimplify[dtnAt[dtn, ell]["ImagCoeffK"] - expectedCoeff[ell]]^2,
      {ell, {0, 1, 2}}
    ]
  ]
];

kernelResidual[dtn_] := FullSimplify[
  Total[
    Table[
      FullSimplify[dtnAt[dtn, ell]["RadiationKernel"] - expectedKernel[ell]]^2,
      {ell, {0, 1, 2}}
    ]
  ]
];

copiedWithCorruptEll2Kernel[dtn_] := Join[
  dtn,
  <|2 -> Join[dtnAt[dtn, 2], <|"RadiationKernel" -> I a^5 omega^5/(9 cS^5)|>]|>
];

sourceData[dtn_] := Module[{sources, raw, residual, without, withReturn, derivativeVertex},
  sources = <|0 -> M0, 1 -> D1, 2 -> Q2|>;
  raw = AssociationThread[
    {0, 1, 2},
    FullSimplify[dtnAt[dtn, #]["RadiationKernel"] sources[#]] & /@ {0, 1, 2}
  ];
  residual = <|
    0 -> FullSimplify[dtnAt[dtn, 0]["RadiationKernel"] (M0 + R0)],
    1 -> FullSimplify[dtnAt[dtn, 1]["RadiationKernel"] (D1 + R1)]
  |>;
  without = <|
    0 -> FullSimplify[residual[0] /. R0 -> 0],
    1 -> FullSimplify[residual[1] /. R1 -> 0]
  |>;
  withReturn = <|
    0 -> FullSimplify[residual[0] /. R0 -> -M0],
    1 -> FullSimplify[residual[1] /. R1 -> -D1]
  |>;
  derivativeVertex = <|
    0 -> FullSimplify[eta0^2 omega^2 raw[0]],
    1 -> FullSimplify[eta1^2 omega^2 raw[1]]
  |>;
  <|
    "Sources" -> sources,
    "Raw" -> raw,
    "Residual" -> residual,
    "Without" -> without,
    "With" -> withReturn,
    "DerivativeVertex" -> derivativeVertex
  |>
];

rawAmplitudeResidual[src_] := FullSimplify[
  Total[
    Table[
      FullSimplify[src["Raw"][ell] - expectedKernel[ell] src["Sources"][ell]]^2,
      {ell, {0, 1, 2}}
    ]
  ]
];

conditionWorks[src_] := TrueQ[
  nonzeroQ[src["Without"][0]] && nonzeroQ[src["Without"][1]] &&
    src["With"][0] === 0 && src["With"][1] === 0
];

residualWithCorruptReturn[src_] := FullSimplify[
  (src["Residual"][0] /. R0 -> -2 M0)^2 +
    (src["Residual"][1] /. R1 -> -D1)^2
];

computeVerdict[raw_, works_, possible_] := Which[
  TrueQ[raw] && ! TrueQ[possible], "MONOPOLE_RADIATION_UNAVOIDABLE",
  TrueQ[raw] && TrueQ[works] && TrueQ[possible], "MONOPOLE_DIPOLE_RETURN_CONDITIONAL",
  True, "INCONCLUSIVE"
];

runAritySelfCheck[] := Module[{probeH, probeR},
  subheading["Wolfram arity self-check"];
  probeH = hankelComponents[0, Unique["zprobe"], False];
  probeR = radiatingPower[1 + I k, 2];
  expectBool["arity hankelComponents[ell,z,corrupt] returns a pair", MatchQ[probeH, {_, _}]];
  expectBool["arity radiatingPower[y,max] returns {power,coefficient}", probeR === {1, 1}];
  expectBool[
    "arity computeVerdict[raw,works,possible] returns a label",
    computeVerdict[True, True, True] === "MONOPOLE_DIPOLE_RETURN_CONDITIONAL"
  ]
];

runDtnLadder[] := Module[{dtn},
  subheading["Earned outgoing Hankel DtN ladder"];
  dtn = dtnData[False];
  Do[
    Print["  ell=", ell, ": h_ell = ", fmt[dtnAt[dtn, ell]["h"]]];
    Print["    Lambda_ell series = ", fmt[dtnAt[dtn, ell]["LambdaSeries"]]];
    Print["    static-normalized admittance = ", fmt[dtnAt[dtn, ell]["YNormSeries"]]];
    Print[
      "    coefficient-scan first radiation phase: p_raw=",
      dtnAt[dtn, ell]["pRaw"], ", coeff=", fmt[dtnAt[dtn, ell]["ImagCoeffK"]]
    ];
    expectZero[
      "ell=" <> ToString[ell] <> " scanned p_raw and radiation-phase coefficient match ladder",
      (dtnAt[dtn, ell]["pRaw"] - expectedPower[ell])^2 +
        FullSimplify[dtnAt[dtn, ell]["ImagCoeffK"] - expectedCoeff[ell]]^2
    ];
    expectZero[
      "ell=" <> ToString[ell] <> " kernel is i*(omega/c_S)^p times scanned coefficient",
      FullSimplify[dtnAt[dtn, ell]["RadiationKernel"] - expectedKernel[ell]]
    ],
    {ell, {0, 1, 2}}
  ];
  expectZero["DtN ladder residual p=1/3/5 and coefficients a,a^3/2,a^5/27", dtnLadderResidual[dtn]];
  expectZero["radiation kernels match i*a*w/cS, i*a^3*w^3/(2*cS^3), i*a^5*w^5/(27*cS^5)", kernelResidual[dtn]];
  dtn
];

runSourcesAndBookkeeping[dtn_] := Module[{src},
  subheading["Raw amplitudes and return-target bookkeeping"];
  src = sourceData[dtn];
  Print["  Moment definitions printed as target definitions:"];
  Print["    M0(omega)=int_brane S_leak(omega,x) d^3x"];
  Print["    D1_i(omega)=int_brane x_i S_leak(omega,x) d^3x + int_brane j_i(omega,x) d^3x, including the carried odd wake"];
  Print["    Q2 is a FREE ANCHOR symbol for the downstream quadrupole consumer, not derived ell=2 physics here"];
  Do[
    Print[
      "  ell=", ell, ": A_raw = kernel*",
      Switch[ell, 0, "M0", 1, "D1", 2, "Q2"],
      " = ", fmt[src["Raw"][ell]]
    ];
    expectNonzero["ell=" <> ToString[ell] <> " raw amplitude is present", src["Raw"][ell]],
    {ell, {0, 1, 2}}
  ];
  expectZero["raw amplitude closed forms match report lines 17-19", rawAmplitudeResidual[src]];
  Print["  Residual label: x−x bookkeeping identity, NOT an earned cancellation result"];
  Do[
    Print["    ell=", ell, " without return: ", fmt[src["Without"][ell]]];
    Print["    ell=", ell, " with target return: ", fmt[src["With"][ell]]];
    expectNonzero["ell=" <> ToString[ell] <> " residual without return condition is nonzero", src["Without"][ell]];
    expectZero["ell=" <> ToString[ell] <> " residual with bookkeeping target is exactly zero", src["With"][ell]],
    {ell, {0, 1}}
  ];
  Print["  Required targets: R0(omega)=-M0(omega); R1_i(omega)=-D1_i(omega)"];
  Print["  Raw-vs-vertex note: derivative outlet vertex g_W0(omega)=eta*omega is branch_assumption; two vertices add eta^2*omega^2 and are NOT verdict-bearing."];
  Do[
    Print["    derivative-vertex branch ell=", ell, ": ", fmt[src["DerivativeVertex"][ell]]];
    expectNonzero[
      "branch_assumption derivative-vertex ell=" <> ToString[ell] <> " remains nonzero",
      src["DerivativeVertex"][ell]
    ],
    {ell, {0, 1}}
  ];
  expectZero[
    "steady control lim_{omega->0} raw monopole amplitude is zero",
    Block[
      {$Assumptions = a > 0 && cS > 0 && Element[M0, Reals]},
      FullSimplify[Limit[src["Raw"][0], omega -> 0]]
    ]
  ];
  expectBool["dominance p(ell0)<p(ell2) from scanned powers", dtnAt[dtn, 0]["pRaw"] < dtnAt[dtn, 2]["pRaw"]];
  src
];

stageAnchors[] := Module[
  {w, W243, jw243, sleak243, target243, W244, phi, Ew, jw244, sleak244, target244},
  subheading["Stage-243/244 S_leak anchors by direct integration"];
  W243 = Exp[-w^2]/Sqrt[Pi];
  jw243 = ellw j0 w Exp[-w^2];
  sleak243 = FullSimplify[Integrate[D[W243, w] jw243, {w, -Infinity, Infinity}]];
  target243 = -Sqrt[2] ellw j0/4;
  W244 = Exp[-w^2/lambda^2]/(lambda Sqrt[Pi]);
  phi = 2 w Exp[-w^2/lambda^2]/(Sqrt[Pi] lambda^3);
  Ew = -E0 phi;
  jw244 = muW rho0 Ew;
  sleak244 = FullSimplify[
    Integrate[D[W244, w] jw244, {w, -Infinity, Infinity}],
    lambda > 0 && muW > 0 && rho0 > 0 && Element[E0, Reals]
  ];
  target244 = Sqrt[2] muW rho0 E0/(2 Sqrt[Pi] lambda^3);
  Print["  stage-243 direct integral int W' j^w dw = ", fmt[sleak243]];
  Print["  stage-244 direct integral S_leak = ", fmt[sleak244]];
  Print["  ledger_stage006 (I-3) OWNS the recovery-reduction derivation of the projected law; here the two closed forms serve only as live-source consistency anchors."];
  expectZero["stage-243 anchor equals -sqrt(2)*ell_w*j0/4", sleak243 - target243];
  expectZero["stage-244 anchor equals sqrt(2)*mu_w*rho0*E0/(2*sqrt(pi)*lambda^3)", sleak244 - target244];
  expectNonzero["stage-243 leakage lane is symbolically nonzero before recovery shortcut", sleak243];
  expectNonzero["stage-244 leakage lane is symbolically nonzero before recovery shortcut", sleak244];
  <|
    "stage243" -> sleak243,
    "target243" -> target243,
    "stage244" -> sleak244,
    "target244" -> target244,
    "ellw" -> ellw,
    "E0" -> E0,
    "muW" -> muW,
    "rho0" -> rho0,
    "lambda" -> lambda
  |>
];

runVerdict[src_] := Module[{rawPresent, works, cancellationPossible, verdict, synthetic, inconclusive},
  subheading["Computed verdict with honest literal flag"];
  rawPresent = AllTrue[{0, 1, 2}, nonzeroQ[src["Raw"][#]] &];
  works = conditionWorks[src];
  cancellationPossible = True;
  Print["  raw_present computed from amplitudes = ", rawPresent];
  Print["  condition_works computed from residual pair = ", works];
  Print["  SCOPE: parameter, not computed — track-3 decides: cancellation_possible=True"];
  verdict = computeVerdict[rawPresent, works, cancellationPossible];
  synthetic = computeVerdict[True, False, False];
  inconclusive = computeVerdict[False, False, True];
  Print["  baseline verdict = ", verdict];
  Print["  synthetic no-return control = ", synthetic];
  expectBool["baseline computed verdict is MONOPOLE_DIPOLE_RETURN_CONDITIONAL", verdict === "MONOPOLE_DIPOLE_RETURN_CONDITIONAL"];
  expectBool["synthetic (True,False,False) reaches MONOPOLE_RADIATION_UNAVOIDABLE", synthetic === "MONOPOLE_RADIATION_UNAVOIDABLE"];
  expectBool["INCONCLUSIVE verdict rung remains reachable", inconclusive === "INCONCLUSIVE"];
  {verdict, rawPresent, works}
];

runGuards[dtn_, src_, anchors_] := Module[
  {rawPresent, steadyLimit, recovery243, recovery244, derivativeNotBasis, quadrupoleSurvives, noTrack3BulkKill},
  subheading["Nine source-live guards"];
  Print["  These controls do NOT test whether suppression occurs. What they confirm is narrow: the source moments M0/D1 are kept live (no S_leak=0, no strict-recovery basis, no projection-locking that would zero them out by construction). Beyond keeping the source live, the controls pass by construction — they are not able-to-fail probes of the physical question. Treat them as guards against the obvious tautologies, not as evidence of suppression-vs-unavoidable."];
  rawPresent = AllTrue[{0, 1, 2}, nonzeroQ[src["Raw"][#]] &];
  steadyLimit = Block[
    {$Assumptions = a > 0 && cS > 0 && Element[M0, Reals]},
    FullSimplify[Limit[src["Raw"][0], omega -> 0]]
  ];
  recovery243 = FullSimplify[anchors["stage243"] /. anchors["ellw"] -> 0];
  recovery244 = FullSimplify[anchors["stage244"] /. anchors["E0"] -> 0];
  derivativeNotBasis = nonzeroQ[src["DerivativeVertex"][0]];
  quadrupoleSurvives = TrueQ[dtnAt[dtn, 2]["pRaw"] == 5] && nonzeroQ[src["Raw"][2]];
  noTrack3BulkKill = nonzeroQ[src["Without"][0]] && nonzeroQ[src["Without"][1]];
  expectBool["guard raw_monopole_present: raw monopole is live in same pipeline", nonzeroQ[src["Raw"][0]] && rawPresent];
  expectZero["guard steady_no_radiation: lim_{omega->0} raw0=0", steadyLimit];
  expectBool["guard quadrupole_survives: scanned p(ell2)=5 and raw2 nonzero", quadrupoleSurvives];
  expectBool["guard return_necessity: without nonzero and with exactly zero", conditionWorks[src]];
  Print["    declaration anti_tautology_no_S_leak_zero: M0->0 killing raw0 is pass-by-construction after the raw-amplitude closed-form checks, not a counted able-to-fail tooth."];
  Print["    declaration anti_tautology_no_S_leak_zero: no S_leak=0 shortcut is used as a verdict basis."];
  Print["    declaration anti_tautology_no_strict_recovery_basis: ell_w->0 and E0->0 recovery slices are not used as verdict bases."];
  Print["    observed_but_not_used — the strict-recovery limit exists but is NOT taken as a basis: sleak243|_{ell_w->0} = 0"];
  expectZero["observed_but_not_used strict-recovery slice sleak243|_{ell_w->0}=0", recovery243];
  Print["    observed_but_not_used — the strict-recovery limit exists but is NOT taken as a basis: sleak244|_{E0->0} = 0"];
  expectZero["observed_but_not_used strict-recovery slice sleak244|_{E0->0}=0", recovery244];
  Print["    declaration anti_tautology_no_projection_locking: M0 and D1 stay unconstrained source moments."];
  expectBool["guard anti_tautology_derivative_vertex_not_basis: branch nonzero but not verdict-bearing", derivativeNotBasis];
  expectBool["guard anti_tautology_no_track3_bulk_kill: without-return residuals stay nonzero", noTrack3BulkKill];
  Print["    declaration anti_tautology_no_track3_bulk_kill: no track-3 bulk return kill is imported here."]
];

runFrozenAndConsumedStage005[] := Module[
  {gConst, cConst, cSConst, kEos, rho, mGNLS, aStar, lStar, consumedResidual},
  subheading["Frozen slice and consumed ledger_stage005 law"];
  gConst = 1;
  cConst = 1;
  cSConst = 1;
  kEos = 1/500;
  rho = Sqrt[10];
  mGNLS = 1;
  aStar = 4731/2500;
  lStar = 18121/10000;
  Print["  Frozen slice (CALIBRATED, exact): G=c=c_s=1; K_eos=1/500; rho=sqrt(10); (a*,L*)=(4731/2500,18121/10000)"];
  Print["  Consumed from ledger_stage005 (I-2)"];
  Print["    cited symbolic law: c_s^2 = 5*K*rho^4/m_GNLS"];
  Print["    no EOS re-derivation is performed in stage008; this is citation integrity only."];
  consumedResidual = FullSimplify[5 kEos rho^4/mGNLS - 1];
  Print["    frozen slice (CALIBRATED, cited): G=", fmt[gConst]];
  Print["    frozen slice (CALIBRATED, cited): c=", fmt[cConst]];
  Print["    frozen slice (CALIBRATED, cited): c_s=", fmt[cSConst]];
  Print["    frozen slice (CALIBRATED, cited): K_eos=", fmt[kEos]];
  expectZero["frozen rho=sqrt(10) symbolic exact", rho^2 - 10];
  Print["    frozen slice (CALIBRATED, cited): a*=", fmt[aStar]];
  Print["    frozen slice (CALIBRATED, cited): L*=", fmt[lStar]];
  expectZero["consumed stage005 law exact-value citation-integrity 5*(1/500)*(sqrt(10))^4/1 - 1", consumedResidual];
  <|"K_eos" -> kEos, "rho" -> rho, "m_GNLS" -> mGNLS|>
];

dimResidualVec[actual_, expected_] := FullSimplify[(actual - expected).(actual - expected)];

expressionDim[expr_, symbolDims_] := Module[{e, dims},
  e = FullSimplify[expr];
  Which[
    NumberQ[e] || e === I || e === Pi, {0, 0},
    KeyExistsQ[symbolDims, e], symbolDims[e],
    AtomQ[e], raise["missing sourced dimension for symbol " <> ToString[InputForm[e]]],
    Head[e] === Times, Total[expressionDim[#, symbolDims] & /@ (List @@ e)],
    Head[e] === Power,
      If[
        MatchQ[e[[2]], _Integer | _Rational],
        e[[2]] expressionDim[e[[1]], symbolDims],
        raise["non-numeric exponent in dimension walk: " <> ToString[InputForm[e]]]
      ],
    Head[e] === Plus,
      dims = expressionDim[#, symbolDims] & /@ (List @@ e);
      If[
        AllTrue[Rest[dims], # === First[dims] &],
        First[dims],
        raise["dimensionally inhomogeneous expression: " <> ToString[InputForm[e]]]
      ],
    True, raise["unsupported expression in dimension walk: " <> ToString[InputForm[e]]]
  ]
];

kernelDimensionResidual[kernel_, symbolDims_] := dimResidualVec[expressionDim[kernel, symbolDims], {0, 0}];

runDimensionalBlock[dtn_] := Module[
  {zero, dimA, dimCS, dimOmega, kernelSymbolDims, corruptKernelSymbolDims, corruptDimResidual},
  subheading["Modest dimensional block"];
  zero = {0, 0};
  dimA = {1, 0};
  dimCS = {1, -1};
  dimOmega = {0, -1};
  kernelSymbolDims = <|a -> dimA, cS -> dimCS, omega -> dimOmega|>;
  corruptKernelSymbolDims = Join[kernelSymbolDims, <|a -> dimA + {1, 0}|>];
  Print["  structural note (not counted): z=k*a is dimensionless by k=omega/c_S and [a]=L."];
  Do[
    expectZero[
      "ell=" <> ToString[ell] <> " actual radiation_kernel dimension from sourced a,c_S,omega is dimensionless",
      kernelDimensionResidual[dtnAt[dtn, ell]["RadiationKernel"], kernelSymbolDims]
    ],
    {ell, {0, 1, 2}}
  ];
  corruptDimResidual = FullSimplify[
    Total[
      Table[
        kernelDimensionResidual[dtnAt[dtn, ell]["RadiationKernel"], corruptKernelSymbolDims],
        {ell, {0, 1, 2}}
      ]
    ]
  ];
  expectFail["dimensional source corruption [a]->L^2 makes actual radiation kernels dimensionful", corruptDimResidual];
  expectBool[
    "scanned ladder powers are strictly increasing",
    dtnAt[dtn, 0]["pRaw"] < dtnAt[dtn, 1]["pRaw"] < dtnAt[dtn, 2]["pRaw"]
  ];
  Print["  No dimensions are invented here for M0/D1/Q2; they remain frozen-slice normalized source moments/anchor."]
];

printScopeAndProvenance[] := (
  subheading["Scope caveat and provenance"];
  Print["  This is a VERIFIED CONSTRAINT-SPECIFICATION, not a falsifiable test of no-monopole-radiation. What is verified-solid: the raw ell=0/1/2 outgoing amplitudes and their radiation orders, and the exact moments/orders at which any brane<->bulk return must cancel. What this report does NOT establish: that monopole/dipole radiation is suppressed or unavoidable."];
  Print["  The cancellation condition (R0=-M0, R1=-D1) is the algebraic bookkeeping of what the return must cancel. Its symbolic derivation is an identity (x - x), not a deep load-bearing result."];
  Print["  The UNAVOIDABLE rung is NOT decidable in-scope: cancellation_possible is a parameter (a literal flag), not a computed quantity. Deciding suppression-vs-unavoidable requires the track-3 brane<->bulk return, which is out of scope here."];
  Print["  The real falsification lives in track 3: whether an admissible return can actually deliver R0=-M0, R1=-D1. This audit only specifies the target it must hit."];
  Print[""];
  Print["  Exports:"];
  Print["    R0=-M0 and R1=-D1 targets + raw amplitudes/kernels -> ledger stages 009/010 (pathA_29) and 022/023 (pathA_34)."];
  Print["    Q2 = FREE ANCHOR (not derived ell=2 physics) -> ledger stage 026."];
  Print["  Provenance:"];
  Print["    DtN/outgoing expansions cite reuse of research/4d_2_5pn spherical-Hankel machinery; stage029 is the later formal DOI-cite home."];
  Print["    M0/D1 are this gate's own constructions consistent with old-ledger Part-VIII projected continuity, NOT verbatim Part-VIII objects."];
  Print["    ledger_stage006 (I-3) OWNS the recovery-reduction derivation of the projected law; stage008 uses stage-243/244 only as live-source consistency anchors."];
  Print["    Raw-vs-vertex: derivative outlet vertex is branch_assumption, recorded but not used for the verdict."];
  Print["    consumed from ledger_stage005 (I-2): c_s^2 = 5*K*rho^4/m_GNLS (cited, integrity-checked)"]
);

printVerdictLabels[verdict_] := (
  Print[""];
  Print["Verdict labels:"];
  Print["  ledger earned-label (NOT a source verdict token): DTN_LADDER_DERIVED_RETURN_TARGETS_SPECIFIED  (p_raw=1/3/5 scanned from the Hankel DtN series; raw amplitudes nonzero; cancellation targets bookkept)"];
  Print["  source top-line verdict: ", verdict];
  Print["  honest scope (verbatim-class): VERIFIED CONSTRAINT-SPECIFICATION, not a falsifiable suppression test; R0=-M0 / R1=-D1 is x-x bookkeeping; cancellation_possible is a literal parameter flag (track-3 decides)"];
  Print["  earned: DtN ladder p=1/3/5 + kernels i*ak, i*a^3k^3/2, i*a^5k^5/27; steady limit; dominance ordering; stage-243/244 S_leak anchors (live-source consistency)"];
  Print["  exports: R0=-M0, R1=-D1 targets -> stages 009/010/022/023; Q2 = FREE ANCHOR (not derived l=2) -> stage 026"];
  Print["  guards: 9 carried (computed where they compute; declarations printed as declarations; pass-by-construction framing kept)"];
  Print["  frozen slice (CALIBRATED, cited): G=c=c_s=1; K_eos=1/500; (a*,L*)=(4731/2500,18121/10000)"];
  Print["  consumed from ledger_stage005 (I-2): c_s^2 = 5*K*rho^4/m_GNLS (cited, integrity-checked)"]
);

runAbleToFailTeeth[dtn_, src_, anchors_, consumed_, rawPresent_] := Module[
  {corruptDtn, corruptKernel, corruptVerdict, wrong244, steadyCorruptDtn, steadyCorruptSrc, strippedKernel0, corruptConsumed},
  subheading["Able-to-fail mutation teeth"];
  corruptDtn = dtnData[True];
  expectFail["tooth 1 hankel-form corruption changes scanned ladder", dtnLadderResidual[corruptDtn]];
  corruptKernel = copiedWithCorruptEll2Kernel[dtn];
  expectFail["tooth 2 kernel coefficient corruption a^5/27 -> a^5/9 trips kernel assert", kernelResidual[corruptKernel]];
  expectFail["tooth 3 M0->0 kills raw monopole and trips nonzero assert", nonzeroAssertResidual[src["Raw"][0] /. M0 -> 0]];
  expectFail["tooth 4 R0->-2*M0 breaks x-x bookkeeping identity", residualWithCorruptReturn[src]];
  corruptVerdict = computeVerdict[rawPresent, False, True];
  expectFail[
    "tooth 5 verdict machinery corruption falls to INCONCLUSIVE and trips baseline equality",
    sameStringResidual[corruptVerdict, "MONOPOLE_DIPOLE_RETURN_CONDITIONAL"]
  ];
  wrong244 = Sqrt[2] anchors["muW"] anchors["rho0"] anchors["E0"]/(Sqrt[Pi] anchors["lambda"]^3);
  expectFail["tooth 6 stage-244 prefactor corruption trips anchor assert", anchors["stage244"] - wrong244];
  strippedKernel0 = FullSimplify[
    dtnAt[dtn, 0]["RadiationKernel"] /. Power[omega, dtnAt[dtn, 0]["pRaw"]] -> 1
  ];
  steadyCorruptDtn = Join[
    dtn,
    <|0 -> Join[dtnAt[dtn, 0], <|"RadiationKernel" -> strippedKernel0|>]|>
  ];
  steadyCorruptSrc = sourceData[steadyCorruptDtn];
  expectFail[
    "tooth 7 steady-limit corruption strips omega from scanned kernel0 and trips limit assert",
    Block[
      {$Assumptions = a > 0 && cS > 0 && Element[M0, Reals]},
      FullSimplify[Limit[steadyCorruptSrc["Raw"][0], omega -> 0]]
    ]
  ];
  corruptConsumed = FullSimplify[4 consumed["K_eos"] consumed["rho"]^3/consumed["m_GNLS"] - 1];
  expectFail["tooth 8 consumed-law corruption 5*K*rho^4 -> 4*K*rho^3 trips citation integrity", corruptConsumed];
  expectZero["baseline immutable after copy-mutation teeth: DtN ladder unchanged", dtnLadderResidual[dtn]];
  expectZero["baseline immutable after copy-mutation teeth: kernels unchanged", kernelResidual[dtn]];
  expectZero[
    "baseline immutable after copy-mutation teeth: bookkeeping with-target residuals unchanged",
    src["With"][0]^2 + src["With"][1]^2
  ]
];

Module[{ok, dtn, src, anchors, verdictPack, verdict, rawPresent, consumed},
  heading["ledger_stage008_monopole_dipole_return_spec Mathematica audit"];
  ok = Catch[
    runAritySelfCheck[];
    dtn = runDtnLadder[];
    src = runSourcesAndBookkeeping[dtn];
    anchors = stageAnchors[];
    verdictPack = runVerdict[src];
    verdict = verdictPack[[1]];
    rawPresent = verdictPack[[2]];
    runGuards[dtn, src, anchors];
    consumed = runFrozenAndConsumedStage005[];
    runDimensionalBlock[dtn];
    printScopeAndProvenance[];
    printVerdictLabels[verdict];
    runAbleToFailTeeth[dtn, src, anchors, consumed, rawPresent];
    True,
    "ledgerStage008Failure",
    Function[{msg, tag}, Print["FAIL: ", msg]; False]
  ];

  Print[""];
  Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
  If[TrueQ[ok],
    Print["OVERALL PASS: Mathematica verified ledger_stage008 monopole/dipole return constraint spec exactly"];
    Exit[0],
    Print["OVERALL FAIL: Mathematica stage008 audit did not close"];
    Exit[1]
  ]
]
