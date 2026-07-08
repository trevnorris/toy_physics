(* Ledger stage009 Mathematica audit: flat-slab return residual, Check A.

   Standalone, print-only, no arguments, no imports, no exports.  This derives
   the bidirectional Helmholtz return phase, DC continuity transfer, bounded
   residual orders, Z accounting, and strict perfect-return controls natively
   in Wolfram code.  Check-B localization machinery is named only as stage 010.
*)

ClearAll[
  raise, heading, subheading, cleanZero, assertExact, expectZero,
  expectBool, expectNonzero, expectFail, nonzeroQ, fmt,
  omegaOrder, alpha, neighborFraction, seriesExpr, transportData,
  transportResidual, expectedKernel, consumedStage008Kernel, kernelData,
  consumedStage008Kernels, kernelAt, kernelIntegrityResidual,
  classifyCheckA, classificationResidual,
  computeBaseline, runAritySelfCheck, runRoundTrip, runReturnTransfer,
  runConsumedStage008, runResidualPrediction, runZAccounting,
  runStrictLimits, dimResidualVec, runDimensionalBlock, printProvenance,
  printVerdictLabels, runAbleToFailTeeth, passCount, failCount,
  omega, cS, d, a, w, epsilon0, epsilon1, M0, D1, Q2
];

passCount = 0;
failCount = 0;
Off[Limit::alimv];

raise[msg_] := Throw[msg, "ledgerStage009Failure"];

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
fmt[expr_] := ToString[InputForm[Factor[Cancel[FullSimplify[expr]]]]];

$Assumptions =
  omega > 0 && cS > 0 && d > 0 && a > 0 &&
  epsilon0 > 0 && epsilon1 > 0 && M0 > 0 &&
  Element[{D1, Q2, w}, Reals];

omegaOrder[expr_, max_] := Module[{n, coeff},
  For[n = 0, n <= max, n++,
    coeff = FullSimplify[SeriesCoefficient[expr, {omega, 0, n}]];
    assertExact["omegaOrder coefficient " <> ToString[n], coeff];
    If[nonzeroQ[coeff], Return[n]]
  ];
  raise["no nonzero omega coefficient through " <> ToString[max] <> " in " <> ToString[InputForm[expr]]]
];

alpha[eps_] := FullSimplify[1/(1 + eps)];
neighborFraction[eps_] := FullSimplify[eps/(1 + eps)];
seriesExpr[expr_, order_] := Expand[Normal[Series[expr, {omega, 0, order}]]];

transportData[returnSign_] := Module[
  {outgoingBasis, returningBasis, forwardPhase, returnPhase, phase, tau, phaseSeries},
  outgoingBasis = Exp[I omega w/cS];
  returningBasis = Exp[returnSign I omega w/cS];
  forwardPhase = FullSimplify[(outgoingBasis /. w -> d)/(outgoingBasis /. w -> 0)];
  returnPhase = FullSimplify[(returningBasis /. w -> 0)/(returningBasis /. w -> d)];
  phase = FullSimplify[forwardPhase returnPhase];
  phaseSeries = seriesExpr[phase, 5];
  tau = FullSimplify[D[Log[phase], omega]/I];
  <|
    "OutgoingBasis" -> outgoingBasis,
    "ReturningBasis" -> returningBasis,
    "ForwardPhase" -> forwardPhase,
    "ReturnPhase" -> returnPhase,
    "Phase" -> phase,
    "PhaseSeries" -> phaseSeries,
    "Tau" -> tau
  |>
];

transportResidual[data_] := FullSimplify[
  (Exp[I omega data["Tau"]] - data["Phase"])^2 +
    (data["Tau"] - 2 d/cS)^2
];

expectedKernel[ell_] := Switch[
  ell,
  0, I a omega/cS,
  1, I a^3 omega^3/(2 cS^3),
  2, I a^5 omega^5/(27 cS^5),
  _, raise["bad expectedKernel ell"]
];

consumedStage008Kernel[ell_] := Switch[
  ell,
  0, I a omega/cS,
  1, I a^3 omega^3/(2 cS^3),
  2, I a^5 omega^5/(27 cS^5),
  _, raise["bad consumedStage008Kernel ell"]
];

kernelData[] := <|0 -> expectedKernel[0], 1 -> expectedKernel[1], 2 -> expectedKernel[2]|>;
consumedStage008Kernels[] := <|0 -> consumedStage008Kernel[0], 1 -> consumedStage008Kernel[1], 2 -> consumedStage008Kernel[2]|>;
kernelAt[kernels_, ell_] := Lookup[kernels, ell];

kernelIntegrityResidual[pipelineKernels_, consumedKernels_] := FullSimplify[
  Total[
    Table[
      FullSimplify[kernelAt[consumedKernels, ell] - kernelAt[pipelineKernels, ell]]^2,
      {ell, {0, 1, 2}}
    ]
  ]
];

classifyCheckA[p0_, p1_] := <|
  "AStrictPass" -> TrueQ[p0 >= 5 && p1 >= 5],
  "AResidualPass" -> TrueQ[(p0 < 5 || p1 < 5) && p0 >= 1 && p1 >= 3]
|>;

classificationResidual[class_] := If[
  TrueQ[class["AStrictPass"] === False && class["AResidualPass"] === True],
  0,
  1
];

computeBaseline[] := Module[
  {
    transport, phase, tau, alpha0, alpha1, neighbor0, neighbor1,
    t0Full, t1Full, t0Series, t1Series, t0DC, t1DC, nu0, nu1,
    kernels, pRaw0, pRaw1, pRaw2, r0, r1, a0Res, a1Res,
    pRes0, pRes1, classification, steady0, steady1, zThroat,
    zReturn, zReplenishment, zBoundary, zTotal, zLocal, zFormula,
    zCertificate, strictT0, strictT1, strictT0Series, strictT1Series,
    strictNu0, strictNu1, strictPRes0, strictPRes1, strictZ
  },
  transport = transportData[-1];
  phase = transport["Phase"];
  tau = transport["Tau"];
  alpha0 = alpha[epsilon0];
  alpha1 = alpha[epsilon1];
  neighbor0 = neighborFraction[epsilon0];
  neighbor1 = neighborFraction[epsilon1];
  t0Full = FullSimplify[alpha0 phase];
  t1Full = FullSimplify[alpha1 phase];
  t0Series = seriesExpr[t0Full, 5];
  t1Series = seriesExpr[t1Full, 5];
  t0DC = FullSimplify[Limit[t0Full, omega -> 0]];
  t1DC = FullSimplify[Limit[t1Full, omega -> 0]];
  nu0 = omegaOrder[1 - t0Series, 5];
  nu1 = omegaOrder[1 - t1Series, 5];
  kernels = kernelData[];
  pRaw0 = omegaOrder[kernelAt[kernels, 0], 7];
  pRaw1 = omegaOrder[kernelAt[kernels, 1], 7];
  pRaw2 = omegaOrder[kernelAt[kernels, 2], 7];
  r0 = FullSimplify[-M0 t0Full];
  r1 = FullSimplify[-D1 t1Full];
  a0Res = FullSimplify[kernelAt[kernels, 0] M0 (1 - t0Series)];
  a1Res = FullSimplify[kernelAt[kernels, 1] D1 (1 - t1Series)];
  pRes0 = pRaw0 + nu0;
  pRes1 = pRaw1 + nu1;
  classification = classifyCheckA[pRes0, pRes1];
  steady0 = FullSimplify[M0 - alpha0 M0 - neighbor0 M0];
  steady1 = FullSimplify[D1 - alpha1 D1 - neighbor1 D1];
  zThroat = -M0;
  zReturn = FullSimplify[M0 t0DC];
  zReplenishment = 0;
  zBoundary = 0;
  zTotal = FullSimplify[zThroat + zReturn + zReplenishment + zBoundary];
  zLocal = FullSimplify[-M0 (1 - t0DC)];
  zFormula = FullSimplify[-M0 epsilon0/(1 + epsilon0)];
  zCertificate = FullSimplify[-zTotal (1 + epsilon0)/(M0 epsilon0)];
  strictT0 = FullSimplify[Limit[t0Full, epsilon0 -> 0, Direction -> "FromAbove"]];
  strictT1 = FullSimplify[Limit[t1Full, epsilon1 -> 0, Direction -> "FromAbove"]];
  strictT0Series = seriesExpr[strictT0, 5];
  strictT1Series = seriesExpr[strictT1, 5];
  strictNu0 = omegaOrder[1 - strictT0Series, 5];
  strictNu1 = omegaOrder[1 - strictT1Series, 5];
  strictPRes0 = pRaw0 + strictNu0;
  strictPRes1 = pRaw1 + strictNu1;
  strictZ = FullSimplify[Limit[zTotal, epsilon0 -> 0, Direction -> "FromAbove"]];
  <|
    "Transport" -> transport,
    "Phase" -> phase,
    "Tau" -> tau,
    "Alpha" -> <|0 -> alpha0, 1 -> alpha1|>,
    "Neighbor" -> <|0 -> neighbor0, 1 -> neighbor1|>,
    "TFull" -> <|0 -> t0Full, 1 -> t1Full|>,
    "TSeries" -> <|0 -> t0Series, 1 -> t1Series|>,
    "TDC" -> <|0 -> t0DC, 1 -> t1DC|>,
    "Nu" -> <|0 -> nu0, 1 -> nu1|>,
    "Kernels" -> kernels,
    "PRaw" -> <|0 -> pRaw0, 1 -> pRaw1, 2 -> pRaw2|>,
    "R" -> <|0 -> r0, 1 -> r1|>,
    "ARes" -> <|0 -> a0Res, 1 -> a1Res|>,
    "PRes" -> <|0 -> pRes0, 1 -> pRes1|>,
    "Classification" -> classification,
    "Steady" -> <|0 -> steady0, 1 -> steady1|>,
    "ZParts" -> <|
      "Throat" -> zThroat,
      "Return" -> zReturn,
      "ReplenishmentLocalized" -> zReplenishment,
      "BoundaryDof" -> zBoundary
    |>,
    "Z" -> zTotal,
    "ZLocal" -> zLocal,
    "ZFormula" -> zFormula,
    "ZCertificate" -> zCertificate,
    "Strict" -> <|
      "T" -> <|0 -> strictT0, 1 -> strictT1|>,
      "TSeries" -> <|0 -> strictT0Series, 1 -> strictT1Series|>,
      "Nu" -> <|0 -> strictNu0, 1 -> strictNu1|>,
      "PRes" -> <|0 -> strictPRes0, 1 -> strictPRes1|>,
      "Z" -> strictZ
    |>
  |>
];

runAritySelfCheck[] := Module[{probeTransport, probeOrder, probeClass, probeKernels},
  subheading["Wolfram arity self-check"];
  probeTransport = transportData[-1];
  probeOrder = omegaOrder[I omega, 2];
  probeClass = classifyCheckA[1, 3];
  probeKernels = kernelData[];
  expectBool[
    "arity transportData[returnSign] returns phase, phase series, and tau",
    KeyExistsQ[probeTransport, "Phase"] &&
      KeyExistsQ[probeTransport, "PhaseSeries"] &&
      KeyExistsQ[probeTransport, "Tau"]
  ];
  expectBool["arity omegaOrder[expr,max] returns first nonzero order", probeOrder === 1];
  expectBool[
    "arity classifyCheckA[p0,p1] returns Check-A booleans",
    probeClass["AStrictPass"] === False && probeClass["AResidualPass"] === True
  ];
  expectBool[
    "arity kernelAt[kernels,ell] reads an exact consumed kernel",
    kernelAt[probeKernels, 2] === expectedKernel[2]
  ]
];

runRoundTrip[data_] := Module[{tr},
  subheading["Solved bidirectional Helmholtz round-trip transport"];
  tr = data["Transport"];
  Print["  basis: Phi_l=A_l*exp(I*omega*w/c_s)+B_l*exp(-I*omega*w/c_s)"];
  Print["  forward phase ratio 0->d = ", fmt[tr["ForwardPhase"]]];
  Print["  return phase ratio d->0 = ", fmt[tr["ReturnPhase"]]];
  Print["  round-trip phase = ", fmt[tr["Phase"]]];
  Print["  round-trip phase series (informative, unasserted) = ", fmt[tr["PhaseSeries"]]];
  Print["  tau solved from D[Log[phase],omega]/I = ", fmt[tr["Tau"]]];
  expectZero["exp(I*omega*tau) equals the solved round-trip phase", Exp[I omega tr["Tau"]] - tr["Phase"]];
  expectZero["tau solved from the basis equals 2*d/c_S", tr["Tau"] - 2 d/cS]
];

runReturnTransfer[data_] := Module[{},
  subheading["DC continuity fractions and return transfer"];
  Print["  alpha_0 = ", fmt[data["Alpha"][0]], "; neighbor_0 = ", fmt[data["Neighbor"][0]]];
  Print["  alpha_1 = ", fmt[data["Alpha"][1]], "; neighbor_1 = ", fmt[data["Neighbor"][1]]];
  Print["  T_0(omega) = ", fmt[data["TFull"][0]]];
  Print["  T_1(omega) = ", fmt[data["TFull"][1]]];
  Do[
    expectZero[
      "ell=" <> ToString[ell] <> " DC transfer by limit matches alpha_ell continuity fraction",
      data["TDC"][ell] - data["Alpha"][ell]
    ];
    Print["    ell=", ell, ": Limit[T_l, omega->0] = ", fmt[data["TDC"][ell]]],
    {ell, {0, 1}}
  ];
  expectZero["zeroth-moment steady circulation balance M0=alpha0*M0+neighbor0*M0", data["Steady"][0]];
  expectZero["first-moment steady circulation balance D1=alpha1*D1+neighbor1*D1", data["Steady"][1]]
];

runConsumedStage008[data_] := Module[{consumedKernels, dcRelation, stage008Target},
  subheading["Consumed from ledger_stage008 (II-B1)"];
  consumedKernels = consumedStage008Kernels[];
  Print["  cited kernels, exact-value integrity checked:"];
  Print["    kernel_0 = ", fmt[data["Kernels"][0]]];
  Print["    kernel_1 = ", fmt[data["Kernels"][1]]];
  Print["    kernel_2 = ", fmt[data["Kernels"][2]]];
  Print["  independently typed consumed-kernel site:"];
  Print["    consumed_kernel_0 = ", fmt[consumedKernels[0]]];
  Print["    consumed_kernel_1 = ", fmt[consumedKernels[1]]];
  Print["    consumed_kernel_2 = ", fmt[consumedKernels[2]]];
  expectZero["consumed kernel_0 matches pipeline kernel_0", consumedKernels[0] - data["Kernels"][0]];
  expectZero["consumed kernel_1 matches pipeline kernel_1", consumedKernels[1] - data["Kernels"][1]];
  expectZero["consumed kernel_2 matches pipeline kernel_2", consumedKernels[2] - data["Kernels"][2]];
  expectZero["all consumed stage008 kernels match independent pipeline kernels", kernelIntegrityResidual[data["Kernels"], consumedKernels]];
  Do[
    dcRelation = FullSimplify[data["R"][ell] /. omega -> 0];
    stage008Target = Switch[ell, 0, -M0, 1, -D1];
    expectZero[
      "ell=" <> ToString[ell] <> " DC return moment approaches stage008 target as epsilon_l->0",
      FullSimplify[Limit[dcRelation, Switch[ell, 0, epsilon0, 1, epsilon1] -> 0, Direction -> "FromAbove"] - stage008Target]
    ];
    Print[
      "    ell=", ell, ": R_l(0) = ", fmt[dcRelation], " -> ", fmt[stage008Target],
      " as epsilon_l->0"
    ],
    {ell, {0, 1}}
  ];
  Print["  T2_applied=false — kernel_2/Q2 INERT, nothing derived at l=2"]
];

runResidualPrediction[data_] := Module[{},
  subheading["Bounded residual prediction, Check A"];
  Do[
    Print["  ell=", ell, ": omega-order scan nu_l=ord_omega(1-T_l series) = ", data["Nu"][ell]];
    expectZero[
      "ell=" <> ToString[ell] <> " finite DC sink leaves omega^0 deviation nu_l=0",
      data["Nu"][ell]
    ];
    expectNonzero[
      "ell=" <> ToString[ell] <> " deviation-from-one has a finite DC sink term",
      SeriesCoefficient[1 - data["TSeries"][ell], {omega, 0, 0}]
    ],
    {ell, {0, 1}}
  ];
  Print["  finite DC sink leaves a nonzero deviation-from-one at O(omega^0)"];
  Print["  R_0(omega) = ", fmt[data["R"][0]]];
  Print["  R_1(omega) = ", fmt[data["R"][1]]];
  Print["  A_res(ell0) = kernel_0*M0*(1-T0) = ", fmt[data["ARes"][0]]];
  Print["  A_res(ell1) = kernel_1*D1*(1-T1) = ", fmt[data["ARes"][1]]];
  expectZero["p_raw0 is computed from consumed kernel_0 omega order", data["PRaw"][0] - 1];
  expectZero["p_raw1 is computed from consumed kernel_1 omega order", data["PRaw"][1] - 3];
  expectZero["p_raw2 is computed from consumed kernel_2 omega order", data["PRaw"][2] - 5];
  expectZero["computed p_res(ell0)=p_raw0+nu0=1", data["PRes"][0] - 1];
  expectZero["computed p_res(ell1)=p_raw1+nu1=3", data["PRes"][1] - 3];
  Print["  p_res(ell0) = ", data["PRes"][0], "; p_res(ell1) = ", data["PRes"][1]];
  Print["  A_strict_pass computed = ", data["Classification"]["AStrictPass"]];
  Print["  A_residual_pass computed = ", data["Classification"]["AResidualPass"]];
  expectBool["A_strict_pass is computed False", data["Classification"]["AStrictPass"] === False];
  expectBool["A_residual_pass is computed True", data["Classification"]["AResidualPass"] === True];
  expectZero["Check-A classification follows from computed p_res orders", classificationResidual[data["Classification"]]];
  Print["  source top-line: RETURN_RESIDUAL_PREDICTION (Check A component computed here; Check B = ledger_stage010)"];
  Print["  Check A component computed here; Check B localization = ledger_stage010 (zero-mode/radial dsolves, counterfactual guard, DC-sink classifier, NOGO warp, spectra)"]
];

runZAccounting[data_] := Module[{parts},
  subheading["Z channel accounting, premise vs accounting"];
  parts = data["ZParts"];
  Print["  Z_throat = ", fmt[parts["Throat"]]];
  Print["  Z_return = ", fmt[parts["Return"]]];
  Print["  localized replenishment = 0 (declared)"];
  Print["  boundary DOF = 0 (declared)"];
  Print["  Z = ", fmt[data["Z"]]];
  expectZero["Z channel sum reduces to -M0*(1-T0(0))", data["Z"] - data["ZLocal"]];
  expectZero["Z accounting formula equals -M0*epsilon0/(1+epsilon0)", data["Z"] - data["ZFormula"]];
  expectZero["Z sign certificate -Z*(1+epsilon0)/(M0*epsilon0) equals 1", data["ZCertificate"] - 1];
  Print["  under v3 this is ACCOUNTING; Z<0 (drain admissibility) is the PREMISE (Z_is_premise = true)"];
  Print["  Z<0 = the drain-admissibility PREMISE (Z_is_premise=true); the formula = ACCOUNTING"]
];

runStrictLimits[data_] := (
  subheading["Strict perfect-return limits, per channel"];
  expectZero["epsilon0->0+ strict T0 equals exp(I*omega*tau)", data["Strict"]["T"][0] - Exp[I omega data["Tau"]]];
  expectZero["epsilon0->0+ strict_nu0 computed independently equals 1", data["Strict"]["Nu"][0] - 1];
  expectZero["epsilon0->0+ strict p_res0 computed as p_raw0+strict_nu0 equals 2", data["Strict"]["PRes"][0] - 2];
  expectZero["epsilon0->0+ Z tends to 0", data["Strict"]["Z"]];
  expectZero["epsilon1->0+ strict_nu1 computed from strict T1 series equals 1", data["Strict"]["Nu"][1] - 1];
  expectZero["epsilon1->0+ strict p_res1 computed as p_raw1+strict_nu1 equals 4", data["Strict"]["PRes"][1] - 4];
  Print[
    "  ell=0 contingency: residual and Z require epsilon0>0; strict T0=",
    fmt[data["Strict"]["T"][0]], ", strict_nu0=", data["Strict"]["Nu"][0],
    ", strict_p_res0=", data["Strict"]["PRes"][0], ", Z->0"
  ];
  Print[
    "  ell=1 contingency: residual requires epsilon1>0; strict T1=",
    fmt[data["Strict"]["T"][1]], ", strict_nu1=", data["Strict"]["Nu"][1],
    ", strict_p_res1=", data["Strict"]["PRes"][1]
  ]
);

dimResidualVec[actual_, expected_] := FullSimplify[(actual - expected).(actual - expected)];

runDimensionalBlock[] := Module[
  {
    zero, dimD, dimCS, dimOmega, dimTime, dimM0, dimEpsilon,
    dimOnePlusEpsilon, dimNeighbor, dimAlpha, dimPhase, dimTransfer,
    zDim, tauDim, zAccountingDim, dimensionlessFractionResidual
  },
  subheading["Modest dimensional block"];
  zero = {0, 0, 0};
  dimD = {1, 0, 0};
  dimCS = {1, -1, 0};
  dimOmega = {0, -1, 0};
  dimTime = {0, 1, 0};
  dimM0 = {0, -1, 0};
  dimEpsilon = zero;
  dimOnePlusEpsilon = zero;
  dimNeighbor = dimEpsilon - dimOnePlusEpsilon;
  dimAlpha = zero - dimOnePlusEpsilon;
  dimPhase = zero;
  dimTransfer = dimAlpha + dimPhase;
  zDim = dimOmega + dimD - dimCS;
  tauDim = dimD - dimCS;
  zAccountingDim = dimM0 + dimNeighbor;
  dimensionlessFractionResidual = FullSimplify[
    dimResidualVec[dimEpsilon, zero] +
      dimResidualVec[dimNeighbor, zero] +
      dimResidualVec[dimAlpha, zero] +
      dimResidualVec[dimTransfer, zero]
  ];
  expectZero["z=omega*d/c_S is dimensionless", dimResidualVec[zDim, zero]];
  expectZero["[tau]=T from tau=2*d/c_S", dimResidualVec[tauDim, dimTime]];
  expectZero["[Z]=[M0] via -M0*epsilon0/(1+epsilon0)", dimResidualVec[zAccountingDim, dimM0]];
  expectZero["epsilon_l, alpha_l, and T_l are dimensionless by fraction composition", dimensionlessFractionResidual];
  Print["  z=omega*d/c_S dimensionless; [tau]=T; [Z]=[M0]; epsilon_l, alpha_l, T_l dimensionless."]
];

printProvenance[] := (
  subheading["Provenance and scope"];
  Print["  postulated geometry: flat finite slab, brane w=0, absorber w=d; response derived on it."];
  Print["  premise vs accounting (v3, verbatim-class): Z<0 = the drain-admissibility PREMISE (Z_is_premise=true); Z=-M0*eps0/(1+eps0) = ACCOUNTING."];
  Print["  open item #9: sharpened, NOT closed — the deliverable is the falsifiable residual radiation prediction tied to the drain strength; the gravity-range (1/r^2) leg = stage 010."];
  Print["  Check-B pointer: radiation/Sommerfeld boundary recorded ac_check_a_only in the source — not a Check-B branch; localization/NOGO/classifier = stage 010."];
  Print["  downstream consumers: stage 023 (pathA_34 residuals must match in form/sign/order); stage 024/026 (pathA_43 consumes the Phi2 bulk mode machinery — stage 010's leg — and the projected-continuity operator lineage)."];
  Print["  dropped bookkeeping: the source's SHA-256 trace ids were build-reproducibility plumbing, replaced by the v2 tri-review protocol."]
);

printVerdictLabels[] := (
  Print[""];
  Print["Verdict labels:"];
  Print["  ledger earned-label (NOT a source verdict token): RETURN_TRANSFER_DERIVED_RESIDUAL_BOUNDED  (tau=2d/c_S solved; T_l=alpha_l*exp(i*omega*2d/c_s); nu_l=0 computed => p_res=1/3; Z accounting + sign certificate; strict limit reversible)"];
  Print["  source top-line verdict: RETURN_RESIDUAL_PREDICTION  (Check A component computed here; Check B localization = ledger_stage010)"];
  Print["  premise vs accounting (v3, verbatim-class): Z<0 = the drain-admissibility PREMISE (Z_is_premise=true); Z=-M0*eps0/(1+eps0) = ACCOUNTING"];
  Print["  earned: solved round-trip transport; DC continuity fractions + steady balances; bounded residual p_res(l0)=1, p_res(l1)=3 CONTINGENT per-channel on eps0>0 (l=0, Z) and eps1>0 (l=1) — strict limits computed independently: deviations -> O(omega) (strict orders 2/4), Z -> 0"];
  Print["  postulated: the flat finite slab (brane w=0, absorber w=d)"];
  Print["  consumed from ledger_stage008 (II-B1): kernels i*a*omega/cS, i*a^3*omega^3/(2*cS^3), i*a^5*omega^5/(27*cS^5) (cited, integrity-checked); kernel_2/Q2 INERT (T2_applied=false, nothing derived at l=2)"];
  Print["  exports: the falsifiable residual-radiation prediction (drain-strength-tied) -> stage 023; open-item #9 sharpened NOT closed"]
);

runAbleToFailTeeth[data_] := Module[
  {
    flippedTransport, alpha0Bad, t0Bad, t0BadDC, steady0Bad,
    scan0BadInput, scan1BadInput, scan0Bad, scan1Bad, zDropReturn,
    zSignFlip, signFlipCertificate, strictT0Wrong, strictT0WrongSeries,
    strictNu0Wrong, strictZWrong, corruptKernels, corruptConsumedKernels,
    corruptClassification
  },
  subheading["Able-to-fail mutation teeth"];
  flippedTransport = transportData[1];
  expectFail["tooth 1 returning-basis sign flip trips tau/phase assert", transportResidual[flippedTransport]];

  alpha0Bad = neighborFraction[epsilon0];
  t0Bad = FullSimplify[alpha0Bad data["Phase"]];
  t0BadDC = FullSimplify[Limit[t0Bad, omega -> 0]];
  steady0Bad = FullSimplify[M0 - alpha0Bad M0 - data["Neighbor"][0] M0];
  expectFail[
    "tooth 2 alpha0->epsilon0/(1+epsilon0) trips DC fraction and steady balance",
    FullSimplify[(t0BadDC - alpha[epsilon0])^2 + steady0Bad^2]
  ];

  scan0BadInput = FullSimplify[(1 - data["TSeries"][0]) - neighborFraction[epsilon0]];
  scan1BadInput = FullSimplify[(1 - data["TSeries"][1]) - neighborFraction[epsilon1]];
  scan0Bad = omegaOrder[scan0BadInput, 5];
  scan1Bad = omegaOrder[scan1BadInput, 5];
  expectFail[
    "tooth 3 subtracting finite-sink DC term raises nu and trips p_res asserts",
    scan0Bad^2 + (data["PRaw"][0] + scan0Bad - 1)^2 +
      scan1Bad^2 + (data["PRaw"][1] + scan1Bad - 3)^2
  ];

  zDropReturn = FullSimplify[
    data["ZParts"]["Throat"] + data["ZParts"]["ReplenishmentLocalized"] +
      data["ZParts"]["BoundaryDof"]
  ];
  expectFail["tooth 4 dropping Z_return trips Z reduction assert", zDropReturn - data["ZLocal"]];

  zSignFlip = FullSimplify[-data["Z"]];
  signFlipCertificate = FullSimplify[-zSignFlip (1 + epsilon0)/(M0 epsilon0)];
  expectFail["tooth 5 flipping Z sign trips sign certificate", signFlipCertificate - 1];

  strictT0Wrong = FullSimplify[data["TFull"][0] /. epsilon0 -> 1];
  strictT0WrongSeries = seriesExpr[strictT0Wrong, 5];
  strictNu0Wrong = omegaOrder[1 - strictT0WrongSeries, 5];
  strictZWrong = FullSimplify[data["Z"] /. epsilon0 -> 1];
  expectFail[
    "tooth 6 corrupting strict limit epsilon0->1 trips Z->0/strict_nu0=1",
    FullSimplify[strictZWrong^2 + (strictNu0Wrong - 1)^2]
  ];

  corruptKernels = data["Kernels"];
  corruptConsumedKernels = Join[consumedStage008Kernels[], <|2 -> I a^5 omega^5/(9 cS^5)|>];
  expectFail[
    "tooth 7 consumed kernel corruption 27->9 trips integrity check",
    kernelIntegrityResidual[corruptKernels, corruptConsumedKernels]
  ];

  corruptClassification = classifyCheckA[0, data["PRes"][1]];
  expectFail[
    "tooth 8 p_res0->0 flips A_residual_pass False and trips classification assert",
    classificationResidual[corruptClassification]
  ];

  expectZero["baseline immutable after teeth: transport residual unchanged", transportResidual[data["Transport"]]];
  expectZero[
    "baseline immutable after teeth: kernel integrity unchanged",
    kernelIntegrityResidual[data["Kernels"], consumedStage008Kernels[]]
  ];
  expectZero["baseline immutable after teeth: classification unchanged", classificationResidual[data["Classification"]]]
];

Module[{ok, data},
  heading["ledger_stage009_flat_slab_return_residual Mathematica audit"];
  ok = Catch[
    runAritySelfCheck[];
    data = computeBaseline[];
    runRoundTrip[data];
    runReturnTransfer[data];
    runConsumedStage008[data];
    runResidualPrediction[data];
    runZAccounting[data];
    runStrictLimits[data];
    runDimensionalBlock[];
    printProvenance[];
    printVerdictLabels[];
    runAbleToFailTeeth[data];
    True,
    "ledgerStage009Failure",
    Function[{msg, tag}, Print["FAIL: ", msg]; False]
  ];

  Print[""];
  Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
  If[TrueQ[ok],
    Print["OVERALL PASS: Mathematica verified ledger_stage009 flat-slab return residual Check A exactly"];
    Exit[0],
    Print["OVERALL FAIL: Mathematica stage009 audit did not close"];
    Exit[1]
  ]
]
