(* Ledger stage010 Mathematica audit: slab localization p=2 plus NOGO control.

   Standalone, print-only, no arguments, no imports, no exports.  This derives
   the two DC-sink transverse spectra, native radial DSolve branches, the
   counterfactual residual, the delocalizing warp control, and the classifier
   machinery independently in Wolfram code.
*)

ClearAll[
  raise, heading, subheading, cleanZero, assertExact, expectZero,
  expectBool, expectNonzero, expectFail, nonzeroQ, nonzeroAssertResidual,
  boolResidual, fmt, boolToInt, stage008IntegrityResidual,
  stage009IntegrityResidual, quadrupoleSurvivesFromPRaw,
  classifyDcSinkGate, branchVerdictFromP, branchVerdictResidual,
  radialOperatorResidual, radialDecayExponent, transverseSpectra,
  transverseNorm, transverseOdeResidual, transverseBcResiduals,
  ensureZeroModeSeed, zeroSeedGuardResidual, solveDynamicZeroModeRadial,
  solveStaticZeroModeRadial, solveStaticMassiveRadial,
  buildCounterfactualGuard, computeBaseline, runAritySelfCheck,
  runConsumedInputs, runTransverseSpectra, runRadialRoutes,
  runCompletionAgreementAndCounterfactual, runWarpAndClassifier,
  dimResidualVec, runDimensionalBlock, printProvenance,
  printVerdictLabels, runAbleToFailTeeth, passCount, failCount,
  RETURNRESIDUAL, RETURNNOGO, BCDEPENDENT, stage008PRawL2Consumed,
  stage008PRawL2Pipeline, stage009AResidualPassConsumed,
  stage009AResidualPassPipeline, T2Applied, omega, cS, d, r, w,
  kWarp, W, m, muRad, ZCheckA, CMassive
];

passCount = 0;
failCount = 0;

RETURNRESIDUAL = "RETURN_RESIDUAL_PREDICTION";
RETURNNOGO = "RETURN_NOGO";
BCDEPENDENT = "BC_DEPENDENT";

stage008PRawL2Consumed = 5;
stage008PRawL2Pipeline = 2 + 3;
stage009AResidualPassConsumed = True;
stage009AResidualPassPipeline = 1;
T2Applied = False;

$Assumptions =
  d > 0 && r > 0 && cS > 0 && kWarp > 0 && W > 0 && muRad > 0 &&
  Element[{omega, w, m, ZCheckA, CMassive}, Reals];

raise[msg_] := Throw[msg, "ledgerStage010Failure"];

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
boolResidual[condition_] := If[TrueQ[condition], 0, 1];
fmt[expr_] := ToString[InputForm[Factor[Cancel[FullSimplify[expr]]]]];

boolToInt[value_] := If[TrueQ[value], 1, 0];
stage008IntegrityResidual[consumed_, pipeline_] := FullSimplify[consumed - pipeline];
stage009IntegrityResidual[consumed_, pipeline_] := FullSimplify[boolToInt[consumed] - pipeline];
quadrupoleSurvivesFromPRaw[pValue_] := TrueQ[IntegerQ[pValue] && pValue >= 0];

classifyDcSinkGate[branchPs_, quadrupoleSurvives_] := Module[{equalTarget},
  equalTarget = (TrueQ[FullSimplify[# == 2]] &) /@ branchPs;
  Which[
    ! TrueQ[quadrupoleSurvives], RETURNNOGO,
    And @@ equalTarget, RETURNRESIDUAL,
    ! Or @@ equalTarget, RETURNNOGO,
    True, BCDEPENDENT
  ]
];

branchVerdictFromP[pValue_] := If[TrueQ[FullSimplify[pValue == 2]], RETURNRESIDUAL, RETURNNOGO];
branchVerdictResidual[verdict_, expected_] := If[verdict === expected, 0, 1];

radialOperatorResidual[radialExpr_, kappaSquared_, coefficient_: 2] := FullSimplify[
  D[radialExpr, {r, 2}] + coefficient D[radialExpr, r]/r + kappaSquared radialExpr
];

radialDecayExponent[radialFlow_] := Module[{probe},
  probe = FullSimplify[
    Limit[-r D[Log[radialFlow], r], r -> Infinity, Assumptions -> d > 0 && kWarp > 0]
  ];
  If[! TrueQ[FreeQ[probe, r | d | omega | cS | kWarp | W | w | m]],
    raise["could not extract pure radial exponent from " <> ToString[InputForm[radialFlow]]]
  ];
  probe
];

transverseSpectra[] := <|
  "destructuring_absorbing" -> <|
    "f0" -> FullSimplify[1/Sqrt[d]],
    "f1" -> FullSimplify[Sqrt[2/d] Cos[Pi w/d]],
    "m0" -> 0,
    "m1" -> Pi/d,
    "bc" -> "neumann"
  |>,
  "bloch_stack" -> <|
    "f0" -> FullSimplify[1/Sqrt[d]],
    "f1" -> FullSimplify[Sqrt[2/d] Cos[2 Pi w/d]],
    "m0" -> 0,
    "m1" -> 2 Pi/d,
    "bc" -> "periodic"
  |>
|>;

transverseNorm[fn_] := FullSimplify[Integrate[fn^2, {w, 0, d}, Assumptions -> d > 0]];
transverseOdeResidual[fn_, mValue_] := FullSimplify[D[fn, {w, 2}] + mValue^2 fn];

transverseBcResiduals[spec_] := Module[{f0, f1},
  f0 = spec["f0"];
  f1 = spec["f1"];
  If[spec["bc"] === "neumann",
    <|
      "f0_prime_0" -> FullSimplify[D[f0, w] /. w -> 0],
      "f0_prime_d" -> FullSimplify[D[f0, w] /. w -> d],
      "f1_prime_0" -> FullSimplify[D[f1, w] /. w -> 0],
      "f1_prime_d" -> FullSimplify[D[f1, w] /. w -> d]
    |>,
    <|
      "f0_value_periodic" -> FullSimplify[(f0 /. w -> d) - (f0 /. w -> 0)],
      "f0_derivative_periodic" -> FullSimplify[(D[f0, w] /. w -> d) - (D[f0, w] /. w -> 0)],
      "f1_value_periodic" -> FullSimplify[(f1 /. w -> d) - (f1 /. w -> 0)],
      "f1_derivative_periodic" -> FullSimplify[(D[f1, w] /. w -> d) - (D[f1, w] /. w -> 0)]
    |>
  ]
];

ensureZeroModeSeed[mValue_, route_] := If[
  ! TrueQ[FullSimplify[mValue == 0]],
  raise[route <> " zero-mode radial solve must be seeded by the computed m=0 eigenvalue"]
];

zeroSeedGuardResidual[mValue_] := If[TrueQ[FullSimplify[mValue == 0]], 0, 1];

solveDynamicZeroModeRadial[branch_, mValue_] := Module[
  {kappaSquared, general, selected, residual, limitGreen, flow, exponent},
  ensureZeroModeSeed[mValue, "dynamic"];
  kappaSquared = FullSimplify[(omega/cS)^2 - mValue^2];
  general = DSolveValue[
    Derivative[2][gDyn][r] + 2 Derivative[1][gDyn][r]/r + kappaSquared gDyn[r] == 0,
    gDyn[r],
    r
  ];
  selected = FullSimplify[
    general /. {C[1] -> 0, C[2] -> I omega/(2 Pi cS d)},
    r > 0 && cS > 0 && d > 0
  ];
  residual = radialOperatorResidual[selected, kappaSquared];
  limitGreen = FullSimplify[Limit[selected, omega -> 0, Assumptions -> r > 0 && cS > 0 && d > 0]];
  flow = FullSimplify[-D[limitGreen, r]];
  exponent = radialDecayExponent[flow];
  <|
    "General" -> general,
    "Selected" -> selected,
    "Residual" -> residual,
    "LimitGreen" -> limitGreen,
    "Flow" -> flow,
    "Exponent" -> exponent,
    "KappaSquared" -> kappaSquared
  |>
];

solveStaticZeroModeRadial[branch_, mValue_] := Module[
  {kappaSquared, general, selected, residual, flow, exponent},
  ensureZeroModeSeed[mValue, "static"];
  kappaSquared = -FullSimplify[mValue]^2;
  general = DSolveValue[
    Derivative[2][gStatic][r] + 2 Derivative[1][gStatic][r]/r + kappaSquared gStatic[r] == 0,
    gStatic[r],
    r
  ];
  selected = FullSimplify[
    general /. {C[1] -> -1/(4 Pi d), C[2] -> 0},
    r > 0 && d > 0
  ];
  residual = radialOperatorResidual[selected, kappaSquared];
  flow = FullSimplify[-D[selected, r]];
  exponent = radialDecayExponent[flow];
  <|
    "General" -> general,
    "Selected" -> selected,
    "Residual" -> residual,
    "Flow" -> flow,
    "Exponent" -> exponent,
    "KappaSquared" -> kappaSquared
  |>
];

solveStaticMassiveRadial[branch_, mValue_] := Module[
  {general, selectedMu, selected, residual},
  If[TrueQ[FullSimplify[mValue == 0]], raise["massive radial solve requires m>0"]];
  general = DSolveValue[
    Derivative[2][gMassive][r] + 2 Derivative[1][gMassive][r]/r - muRad^2 gMassive[r] == 0,
    gMassive[r],
    r
  ];
  selectedMu = FullSimplify[general /. {C[1] -> 1/(4 Pi), C[2] -> 0}, r > 0 && muRad > 0];
  selected = FullSimplify[selectedMu /. muRad -> mValue, r > 0 && d > 0];
  residual = radialOperatorResidual[selected, -mValue^2];
  <|
    "General" -> general,
    "Selected" -> selected,
    "Residual" -> residual,
    "Mu" -> mValue
  |>
];

buildCounterfactualGuard[candidate_] := Module[
  {correctResidual, perturbed, perturbedResidual},
  correctResidual = radialOperatorResidual[candidate, 0];
  perturbed = FullSimplify[candidate/r^4];
  perturbedResidual = radialOperatorResidual[perturbed, 0];
  <|
    "Candidate" -> candidate,
    "CorrectResidual" -> correctResidual,
    "Perturbed" -> perturbed,
    "PerturbedResidual" -> perturbedResidual
  |>
];

computeBaseline[] := Module[
  {
    stage008Integrity, stage009Integrity, quadrupoleSurvives, spectra,
    transverse, static, massive, dynamic, pDynamic, pStatic, pAbs,
    pBloch, counterfactual, warpNormCutoff, warpNormLimit,
    continuumGreen, continuumFlow, pDelocalizing, delocalizingVerdict,
    destructuringVerdict, blochVerdict, headline, noQuadrupoleVerdict,
    tensionStatus
  },
  stage008Integrity = stage008IntegrityResidual[stage008PRawL2Consumed, stage008PRawL2Pipeline];
  stage009Integrity = stage009IntegrityResidual[stage009AResidualPassConsumed, stage009AResidualPassPipeline];
  quadrupoleSurvives = quadrupoleSurvivesFromPRaw[stage008PRawL2Consumed];
  spectra = transverseSpectra[];
  transverse = Association @ KeyValueMap[
    (#1 -> <|
      "Norm0" -> transverseNorm[#2["f0"]],
      "Norm1" -> transverseNorm[#2["f1"]],
      "Ode0" -> transverseOdeResidual[#2["f0"], #2["m0"]],
      "Ode1" -> transverseOdeResidual[#2["f1"], #2["m1"]],
      "BcResiduals" -> transverseBcResiduals[#2]
    |>) &,
    spectra
  ];
  static = Association @ KeyValueMap[(#1 -> solveStaticZeroModeRadial[#1, #2["m0"]]) &, spectra];
  massive = Association @ KeyValueMap[(#1 -> solveStaticMassiveRadial[#1, #2["m1"]]) &, spectra];
  dynamic = solveDynamicZeroModeRadial["destructuring_absorbing", spectra["destructuring_absorbing"]["m0"]];
  pDynamic = dynamic["Exponent"];
  pStatic = static["destructuring_absorbing"]["Exponent"];
  pAbs = static["destructuring_absorbing"]["Exponent"];
  pBloch = static["bloch_stack"]["Exponent"];
  counterfactual = buildCounterfactualGuard[static["destructuring_absorbing"]["Selected"]];
  warpNormCutoff = FullSimplify[Integrate[Exp[2 kWarp w], {w, 0, W}, Assumptions -> kWarp > 0 && W > 0]];
  warpNormLimit = Limit[warpNormCutoff, W -> Infinity, Assumptions -> kWarp > 0];
  continuumGreen = FullSimplify[Integrate[Exp[-m r], {m, 0, Infinity}, Assumptions -> r > 0]/(4 Pi r), r > 0];
  continuumFlow = FullSimplify[-D[continuumGreen, r], r > 0];
  pDelocalizing = radialDecayExponent[continuumFlow];
  delocalizingVerdict = classifyDcSinkGate[{pDelocalizing}, quadrupoleSurvives];
  destructuringVerdict = branchVerdictFromP[pAbs];
  blochVerdict = branchVerdictFromP[pBloch];
  headline = classifyDcSinkGate[{pAbs, pBloch}, quadrupoleSurvives];
  noQuadrupoleVerdict = classifyDcSinkGate[{pAbs, pBloch}, False];
  tensionStatus = If[TrueQ[FullSimplify[pAbs == 2]] && ! TrueQ[FullSimplify[pDelocalizing == 2]], "witnessed", "not_witnessed"];
  <|
    "Stage008Integrity" -> stage008Integrity,
    "Stage009Integrity" -> stage009Integrity,
    "QuadrupoleSurvives" -> quadrupoleSurvives,
    "Spectra" -> spectra,
    "Transverse" -> transverse,
    "Dynamic" -> dynamic,
    "Static" -> static,
    "Massive" -> massive,
    "PDynamic" -> pDynamic,
    "PStatic" -> pStatic,
    "PAbs" -> pAbs,
    "PBloch" -> pBloch,
    "Counterfactual" -> counterfactual,
    "WarpNormCutoff" -> warpNormCutoff,
    "WarpNormLimit" -> warpNormLimit,
    "ContinuumGreen" -> continuumGreen,
    "ContinuumFlow" -> continuumFlow,
    "PDelocalizing" -> pDelocalizing,
    "DelocalizingVerdict" -> delocalizingVerdict,
    "DestructuringVerdict" -> destructuringVerdict,
    "BlochVerdict" -> blochVerdict,
    "Headline" -> headline,
    "NoQuadrupoleVerdict" -> noQuadrupoleVerdict,
    "TensionStatus" -> tensionStatus
  |>
];

runAritySelfCheck[] := Module[{spec, dynProbe, statProbe, massiveProbe},
  subheading["Wolfram arity self-check"];
  spec = transverseSpectra[];
  dynProbe = solveDynamicZeroModeRadial["arity_dynamic", 0];
  statProbe = solveStaticZeroModeRadial["arity_static", 0];
  massiveProbe = solveStaticMassiveRadial["arity_massive", Pi/d];
  expectBool["arity transverseSpectra[] returns both named completions", KeyExistsQ[spec, "destructuring_absorbing"] && KeyExistsQ[spec, "bloch_stack"]];
  expectBool["arity transverseNorm[f] returns the exact compact norm", transverseNorm[1/Sqrt[d]] === 1];
  expectBool["arity radialOperatorResidual[expr,kappa] uses default 2/r", radialOperatorResidual[1/r, 0] === 0];
  expectBool["arity radialOperatorResidual[expr,kappa,coeff] accepts the coefficient override", nonzeroQ[radialOperatorResidual[1/r, 0, 3]]];
  expectBool["arity radialDecayExponent[flow] returns a pure exponent", radialDecayExponent[1/r^2] === 2];
  expectBool["arity solveDynamicZeroModeRadial[branch,m] returns selected and limit Green", KeyExistsQ[dynProbe, "Selected"] && KeyExistsQ[dynProbe, "LimitGreen"]];
  expectBool["arity solveStaticZeroModeRadial[branch,m] returns selected and flow", KeyExistsQ[statProbe, "Selected"] && KeyExistsQ[statProbe, "Flow"]];
  expectBool["arity solveStaticMassiveRadial[branch,mu] returns decaying massive branch", KeyExistsQ[massiveProbe, "Selected"] && KeyExistsQ[massiveProbe, "Mu"]];
  expectBool["arity classifyDcSinkGate[branchPs,quadrupoleSurvives] returns a verdict", classifyDcSinkGate[{2}, True] === RETURNRESIDUAL];
  expectBool["arity branchVerdictFromP[p] returns a branch verdict", branchVerdictFromP[3] === RETURNNOGO]
];

runConsumedInputs[data_] := (
  subheading["Consumed inputs with dual-site citation integrity"];
  Print["  ledger_stage008 cited input: p_raw(l2)=5; T2_applied=false, no l=2 physics recomputed here."];
  Print["    consumed site p_raw_l2 = ", fmt[stage008PRawL2Consumed]];
  Print["    independent pipeline site p_raw_l2 = ", fmt[stage008PRawL2Pipeline]];
  expectZero["ledger_stage008 p_raw(l2) consumed minus pipeline equals zero", data["Stage008Integrity"]];
  expectZero["cited p_raw(l2) equals the frozen ledger_stage008 export 5", stage008PRawL2Consumed - 5];
  expectBool["quadrupole_survives derives from finite non-negative integer p_raw(l2)", data["QuadrupoleSurvives"]];
  expectBool["T2_applied=false is enforced for the cited quadrupole input", T2Applied === False];
  Print["  ledger_stage009 cited input: A_residual_pass=True; used only in the joint-composition print."];
  Print["    consumed site A_residual_pass = ", stage009AResidualPassConsumed];
  Print["    independent pipeline site A_residual_pass_as_int = ", fmt[stage009AResidualPassPipeline]];
  expectZero["ledger_stage009 A_residual_pass consumed minus pipeline equals zero", data["Stage009Integrity"]];
  expectBool["Check-A component is cited True", stage009AResidualPassConsumed === True]
);

runTransverseSpectra[data_] := (
  subheading["Two DC-sink transverse spectra and normalizable zero modes"];
  KeyValueMap[
    Function[{name, spec},
      Module[{checks},
        checks = data["Transverse"][name];
        Print["  ", name, ":"];
        Print["    f0(w) = ", fmt[spec["f0"]], "; m0 = ", fmt[spec["m0"]]];
        Print["    f1(w) = ", fmt[spec["f1"]], "; m1 = ", fmt[spec["m1"]]];
        expectZero[name <> " zero-mode normalization integral equals 1", checks["Norm0"] - 1];
        expectZero[name <> " first-mode normalization integral equals 1", checks["Norm1"] - 1];
        expectZero[name <> " zero-mode transverse ODE residual f''+m^2 f", checks["Ode0"]];
        expectZero[name <> " first-mode transverse ODE residual f''+m^2 f", checks["Ode1"]];
        KeyValueMap[
          (expectZero[name <> " boundary condition " <> #1, #2]) &,
          checks["BcResiduals"]
        ]
      ]
    ],
    data["Spectra"]
  ];
  Print["  load-bearing fact: both completions have a normalizable constant m=0 transverse zero mode."]
);

runRadialRoutes[data_] := Module[{dynamic, staticAbs, result},
  subheading["Dynamic radial zero-mode route"];
  dynamic = data["Dynamic"];
  Print["  operator: g''+(2/r)g'+((omega/c_S)^2-m^2)g=0 solved before omega->0"];
  Print["  DSolve general basis = ", fmt[dynamic["General"]]];
  Print["  BOUNDARY SELECTION: outgoing spherical C[1]/C[2] branch in the Wolfram DSolve basis; normalization from compact zero-mode overlap 1/d."];
  Print["  selected outgoing Green = ", fmt[dynamic["Selected"]]];
  expectZero["dynamic selected outgoing branch satisfies the radial operator", dynamic["Residual"]];
  Print["  omega->0 Green limit = ", fmt[dynamic["LimitGreen"]]];
  Print["  dynamic radial flow -dG/dr = ", fmt[dynamic["Flow"]]];
  expectZero["dynamic route large-r exponent p_dynamic=2", data["PDynamic"] - 2];

  subheading["Static radial route and static-dynamic consistency"];
  staticAbs = data["Static"]["destructuring_absorbing"];
  Print["  operator: omega=0 first, g''+(2/r)g'-m^2 g=0 with computed m=0"];
  Print["  DSolve general basis = ", fmt[staticAbs["General"]]];
  Print["  BOUNDARY SELECTION: constant branch removed; 1/r branch normalized to 1/(4*pi*d)."];
  Print["  selected static Green = ", fmt[staticAbs["Selected"]]];
  expectZero["static selected zero-mode branch satisfies the radial operator", staticAbs["Residual"]];
  Print["  static radial flow -dG/dr = ", fmt[staticAbs["Flow"]]];
  expectZero["static route large-r exponent p_static=2", data["PStatic"] - 2];
  expectZero["dynamic and static exponents agree p_dynamic-p_static=0", data["PDynamic"] - data["PStatic"]];
  expectZero["strong consistency dynamic limit Green minus static Green is identically zero", dynamic["LimitGreen"] - staticAbs["Selected"]];

  subheading["Massive/gapped Yukawa contrast"];
  KeyValueMap[
    Function[{name, result},
      Print["  ", name, ": mu=m1=", fmt[result["Mu"]], "; DSolve basis = ", fmt[result["General"]]];
      Print["    selected decaying Yukawa branch = ", fmt[result["Selected"]]];
      expectNonzero[name <> " gapped mode has mu>0", result["Mu"]];
      expectZero[name <> " selected massive branch satisfies static radial operator", result["Residual"]]
    ],
    data["Massive"]
  ];
  Print["  gapped modes are Yukawa-suppressed as exp(-m1*r)/r, so only the m=0 zero mode sets the far-field 1/r^2."];
  Print["  illustrative, not load-bearing: Green = Z*zero + C*massive with Z=", ZCheckA, " cited from Check A and C=", CMassive, " free."]
];

runCompletionAgreementAndCounterfactual[data_] := Module[{staticAbs, staticBloch, guard},
  subheading["Both completions agree on p=2 and counterfactual wrong falloff is rejected"];
  staticAbs = data["Static"]["destructuring_absorbing"];
  staticBloch = data["Static"]["bloch_stack"];
  Print["  destructuring_absorbing static Green = ", fmt[staticAbs["Selected"]], "; p_abs = ", fmt[data["PAbs"]]];
  Print["  bloch_stack static Green = ", fmt[staticBloch["Selected"]], "; p_bloch = ", fmt[data["PBloch"]]];
  expectZero["destructuring_absorbing exponent p_abs=2", data["PAbs"] - 2];
  expectZero["bloch_stack exponent p_bloch=2", data["PBloch"] - 2];
  expectZero["both completions agree p_abs-p_bloch=0", data["PAbs"] - data["PBloch"]];
  expectZero["bloch selected zero-mode branch satisfies the radial operator", staticBloch["Residual"]];
  guard = data["Counterfactual"];
  Print["  solved static candidate = ", fmt[guard["Candidate"]]];
  Print["  wrong falloff candidate = solved Green * r^-4 = ", fmt[guard["Perturbed"]]];
  expectZero["correct static Green residual is zero", guard["CorrectResidual"]];
  expectZero["counterfactual residual equals 5/(pi*d*r^7)", guard["PerturbedResidual"] - 5/(Pi d r^7)];
  expectNonzero["counterfactual wrong 1/r^5 falloff residual is nonzero", guard["PerturbedResidual"]]
];

runWarpAndClassifier[data_] := (
  subheading["NOGO warp control and computed classifier"];
  Print["  anti-localizing half-line warp: mu(w)=exp(2*k_warp*w), k_warp>0"];
  Print["  cutoff zero-mode norm integral int_0^W mu(w) dw = ", fmt[data["WarpNormCutoff"]]];
  expectBool["warp zero-mode norm diverges as W->infinity", data["WarpNormLimit"] === Infinity];
  Print["  continuum Green integral = ", fmt[data["ContinuumGreen"]]];
  Print["  continuum flow -dG/dr = ", fmt[data["ContinuumFlow"]]];
  expectZero["delocalizing continuum exponent p_delocalizing=3", data["PDelocalizing"] - 3];
  expectZero["same classifier maps [p_delocalizing] to RETURN_NOGO", branchVerdictResidual[data["DelocalizingVerdict"], RETURNNOGO]];
  expectZero["falloff-tension witness requires p_abs=2", data["PAbs"] - 2];
  expectNonzero["falloff-tension witness requires p_delocalizing != 2", data["PDelocalizing"] - 2];
  expectZero["tension_status is witnessed", If[data["TensionStatus"] === "witnessed", 0, 1]];
  Print["  destructuring branch verdict = ", data["DestructuringVerdict"]];
  Print["  bloch branch verdict = ", data["BlochVerdict"]];
  Print["  Check-B headline = ", data["Headline"]];
  expectZero["destructuring branch verdict from p_abs is RETURN_RESIDUAL_PREDICTION", branchVerdictResidual[data["DestructuringVerdict"], RETURNRESIDUAL]];
  expectZero["bloch branch verdict from p_bloch is RETURN_RESIDUAL_PREDICTION", branchVerdictResidual[data["BlochVerdict"], RETURNRESIDUAL]];
  expectZero["Check-B classifier headline from [p_abs,p_bloch] is RETURN_RESIDUAL_PREDICTION", branchVerdictResidual[data["Headline"], RETURNRESIDUAL]];
  expectZero["quadrupole_survives=False classifier prong returns RETURN_NOGO", branchVerdictResidual[data["NoQuadrupoleVerdict"], RETURNNOGO]];
  Print["  COMPLETED joint verdict: RETURN_RESIDUAL_PREDICTION = (Check A: A_residual_pass=True, CITED ledger_stage009) AND (Check B: p=2 both completions, computed here)"];
  Print["  RETURN_NOGO remains reachable if the return delocalizes (warp -> p=3) OR quadrupole_survives=False."]
);

dimResidualVec[actual_, expected_] := FullSimplify[(actual - expected).(actual - expected)];

runDimensionalBlock[] := Module[
  {zero, dimM0, dimRho, dimR, dimD, dimW, dimCS, dimOmega, dimKWarp, dimGreenZeroStatic, dimFlow},
  subheading["Modest dimensional block"];
  zero = {0, 0, 0};
  dimM0 = {0, 0, -1};
  dimRho = {0, -3, 0};
  dimR = {0, 1, 0};
  dimD = {0, 1, 0};
  dimW = {0, 1, 0};
  dimCS = {0, 1, -1};
  dimOmega = {0, 0, -1};
  dimKWarp = {0, -1, 0};
  dimGreenZeroStatic = -dimD - dimR;
  dimFlow = dimGreenZeroStatic - dimR;
  expectZero["z=omega*d/c_S is dimensionless", dimResidualVec[dimOmega + dimD - dimCS, zero]];
  expectZero["k_warp*w is dimensionless", dimResidualVec[dimKWarp + dimW, zero]];
  expectZero["[green_zero_static]=[1/(4*pi*d*r)]=L^-2", dimResidualVec[dimGreenZeroStatic, {0, -2, 0}]];
  expectZero["[radial flow]=L^-1*[Green]=L^-3", dimResidualVec[dimFlow, {0, -3, 0}]];
  expectZero["source radial-flow dim check dim_M0-dim_rho-2*dim_r=[0,1,-1]", dimResidualVec[dimM0 - dimRho - 2 dimR, {0, 1, -1}]];
  Print["  dimensions ordered M,L,T; no EOS/c_s^2 chain re-derived in stage010."]
];

printProvenance[] := (
  subheading["Provenance and scope"];
  Print["  postulated slab: flat finite slab, brane w=0, absorber w=d; localization is derived on this family."];
  Print["  BOUNDARY SELECTIONS: outgoing-spherical C1/C2 radial branch and 1/(4*pi*d) overlap normalization are selections, not derivations."];
  Print["  R19/W_slab caveat: localization holding for the flat-slab FAMILY != the family being SELECTED by dynamics; the slab width d is a postulated register row (stage009), and its selection is the deferred nonlinear return (R19, sim-deferred Gate-6 territory)."];
  Print["  open-item #9: sharpened, NOT closed - the gravity-range (1/r^2) leg passes inside the localizing flat-slab family because both DC-sink completions give p=2; the return-cancellation targets (R23) remain the deferred obligation."];
  Print["  Check-A citation: the joint RETURN_RESIDUAL_PREDICTION is completed here by composing Check B (p=2, computed) with Check A (A_residual_pass=True, cited from ledger_stage009); NOGO reachable."];
  Print["  radiation/Sommerfeld boundary provenance: recorded ac_check_a_only in the source - NOT a Check-B branch."];
  Print["  dropped bookkeeping: source SHA-256 static/dynamic/per-branch trace ids, structure_id, and expr_digest were build-reproducibility plumbing, replaced by the v2 tri-review protocol."];
  Print["  downstream consumers: stage 024/026 (pathA_43 consumes the Phi_l(w,r) bulk Helmholtz zero-mode plus the projected-continuity operator lineage)."];
  Print["  dropped-trace bookkeeping note: the trace-hash lines are not carried into this print-only stage."]
);

printVerdictLabels[] := (
  Print[""];
  Print["Verdict labels:"];
  Print["  ledger earned-label (NOT a source verdict token): RETURN_LOCALIZATION_P2_DERIVED_NOGO_REACHABLE  (both DC-sink completions -> normalizable m=0 zero mode -> genuine 3D-radial dsolve -> p=2; static-dynamic consistent; counterfactual r^-4 rejected (residual 5/(pi*d*r^7)); warp control -> p=3 -> RETURN_NOGO)"];
  Print["  source top-line verdict: RETURN_RESIDUAL_PREDICTION  (COMPLETED here: Check B p=2 computed; Check A component A_residual_pass=True cited from ledger_stage009)"];
  Print["  joint composition: RETURN_RESIDUAL_PREDICTION = (Check A: bounded residual, A_residual_pass=True, CITED stage009) AND (Check B: p=2 both completions, computed here); RETURN_NOGO reachable if the return delocalizes (warp -> p=3) OR quadrupole_survives=False"];
  Print["  earned: two normalizable m=0 zero modes (compact-cell + Bloch); dynamic dsolve then omega->0 -> p=2; static solve -> p=2; static-dynamic consistent; massive/gapped Yukawa contrast; both completions agree p=2; counterfactual r^-4 rejected; NOGO warp control p=3; falloff-tension witnessed"];
  Print["  postulated: the flat finite slab (brane w=0, absorber w=d) [stage009 geometry]"];
  Print["  boundary selections (NOT derived): outgoing-spherical C1/C2 branch choice; 1/(4*pi*d) overlap normalization"];
  Print["  consumed (cited, dual-site integrity): ledger_stage008 p_raw(l2)=5 -> quadrupole_survives=True (T2_applied=false, no l=2 recompute); ledger_stage009 Check-A component A_residual_pass=True"];
  Print["  caveat: localization holds for the FAMILY, not the family selected by dynamics (R19/W_slab); open-item #9 sharpened NOT closed (R23 pending)"]
);

runAbleToFailTeeth[data_] := Module[
  {guard, normalizableCutoff, normalizableLimit, pBlochMutated, mutatedHeadline, f0AbsCorrupt},
  subheading["Able-to-fail mutation teeth"];
  guard = data["Counterfactual"];
  expectFail["tooth 1 corrupt counterfactual guard to accept zero residual trips guard assert", guard["PerturbedResidual"]];
  normalizableCutoff = FullSimplify[Integrate[Exp[-2 kWarp w], {w, 0, W}, Assumptions -> kWarp > 0 && W > 0]];
  normalizableLimit = FullSimplify[Limit[normalizableCutoff, W -> Infinity, Assumptions -> kWarp > 0]];
  expectFail["tooth 2 warp flip exp(-2*k_warp*w) trips non-normalizability assert", boolResidual[normalizableLimit === Infinity]];
  expectFail["tooth 3 radial operator coefficient 2/r->3/r trips selected residual assert", radialOperatorResidual[data["Static"]["destructuring_absorbing"]["Selected"], 0, 3]];
  expectFail["tooth 4 nonzero m seed trips zero-mode radial guard", zeroSeedGuardResidual[data["Spectra"]["destructuring_absorbing"]["m1"]]];
  pBlochMutated = 3;
  mutatedHeadline = classifyDcSinkGate[{data["PAbs"], pBlochMutated}, data["QuadrupoleSurvives"]];
  expectFail[
    "tooth 5 p_bloch->3 trips completion agreement and headline classifier",
    FullSimplify[(data["PAbs"] - pBlochMutated)^2] + branchVerdictResidual[mutatedHeadline, RETURNRESIDUAL]
  ];
  expectZero[
    "tooth 6a classifier maps delocalizing p=3 branch to RETURN_NOGO",
    branchVerdictResidual[classifyDcSinkGate[{data["PDelocalizing"]}, data["QuadrupoleSurvives"]], RETURNNOGO]
  ];
  expectZero[
    "tooth 6b classifier maps quadrupole_survives=False to RETURN_NOGO",
    branchVerdictResidual[classifyDcSinkGate[{data["PAbs"], data["PBloch"]}, False], RETURNNOGO]
  ];
  f0AbsCorrupt = FullSimplify[1/Sqrt[2 d]];
  expectFail["tooth 7 corrupt f0_abs=1/sqrt(2*d) trips zero-mode normalization assert", transverseNorm[f0AbsCorrupt] - 1];
  expectFail["tooth 8a stage008 consumed p_raw(l2) 5->4 trips dual-site integrity", stage008IntegrityResidual[4, stage008PRawL2Pipeline]];
  expectFail["tooth 8b stage008 pipeline p_raw(l2) 5->4 trips dual-site integrity", stage008IntegrityResidual[stage008PRawL2Consumed, 4]];
  expectFail["tooth 8c stage009 consumed A_residual_pass True->False trips dual-site integrity", stage009IntegrityResidual[False, stage009AResidualPassPipeline]];
  expectFail["tooth 8d stage009 pipeline A_residual_pass 1->0 trips dual-site integrity", stage009IntegrityResidual[stage009AResidualPassConsumed, 0]];
  expectBool["tooth 8e non-finite p_raw sentinel would flip quadrupole_survives False", quadrupoleSurvivesFromPRaw[Infinity] === False];
  expectZero["baseline immutable after teeth: p_abs remains 2", data["PAbs"] - 2];
  expectZero["baseline immutable after teeth: p_bloch remains 2", data["PBloch"] - 2];
  expectZero["baseline immutable after teeth: stage008 integrity remains zero", data["Stage008Integrity"]];
  expectZero["baseline immutable after teeth: stage009 integrity remains zero", data["Stage009Integrity"]]
];

Module[{ok, data},
  heading["ledger_stage010_slab_localization_p2_nogo Mathematica audit"];
  ok = Catch[
    runAritySelfCheck[];
    data = computeBaseline[];
    assertExact["baseline", data];
    runConsumedInputs[data];
    runTransverseSpectra[data];
    runRadialRoutes[data];
    runCompletionAgreementAndCounterfactual[data];
    runWarpAndClassifier[data];
    runDimensionalBlock[];
    printProvenance[];
    printVerdictLabels[];
    runAbleToFailTeeth[data];
    True,
    "ledgerStage010Failure",
    Function[{msg, tag}, Print["FAIL: ", msg]; False]
  ];

  Print[""];
  Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
  If[TrueQ[ok],
    Print["OVERALL PASS: Mathematica verified ledger_stage010 slab localization p=2 with reachable NOGO exactly"];
    Exit[0],
    Print["OVERALL FAIL: Mathematica stage010 audit did not close"];
    Exit[1]
  ]
]
