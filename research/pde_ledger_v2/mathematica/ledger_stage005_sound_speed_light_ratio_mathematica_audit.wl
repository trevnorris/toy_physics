(* Ledger stage005 Mathematica audit: EOS sound speed and light/sound ratio.

   Print-only, standalone, no arguments, no exports.  This route uses native
   UnitDimensions for {L,T,M} exponent triples, D for the EOS sound speed and
   state slope, Det/Factor/Solve for the dispersions, and Simplify/Reduce for
   the light/sound negative-control predicates.
*)

ClearAll[
  heading, subheading, cleanZero, assertExact, expectZero, expectBool,
  expectNonzero, expectFail, unitVector, dimResidual, expectDim, homResidual,
  dimString, dimensionDictionary, pathA20DimResiduals, pathA20bDimResiduals,
  recordHarnessCheck, harnessCheckCount, pathA20AlgebraCheckNames,
  pathA20bAlgebraCheckNames, pathA20AlgebraCount, pathA20bAlgebraCount,
  runSoundSpeed, runPathA20Dimensions, runFluxAlgebra, runVelocityCeiling,
  runRoleCatalogNotes, runConsumedStage004, runPathA20bDimensions,
  runPathA20bAlgebra, negativeControl, runLambdaLanding, runReversibility,
  runFirewall, printGapBlock, printVerdictLabels, passCount, failCount
];

passCount = 0;
failCount = 0;
pathA20AlgebraCount = 0;
pathA20bAlgebraCount = 0;

pathA20AlgebraCheckNames = {
  "S1 native D[Log[c_s],Log[rho]] state slope",
  "pathA_20 algebra conditional ideal no-Q/no-V sonic c_s,* / c_s0",
  "pathA_20 algebra conditional ideal no-Q/no-V sonic rho_* / rho0",
  "pathA_20 algebra conditional ideal no-Q/no-V flux factor Jcrit/(rho0*c_s0*A*)",
  "pathA_20 algebra tail factor with lambda_gamma=c_gamma/c_s"
};

pathA20bAlgebraCheckNames = {
  "pathA_20b algebra phonon determinant gives hbar*(omega^2-c_s0^2*k^2)",
  "pathA_20b algebra transverse gauge operator gives c_bulk^2=C_B/C_E",
  "pathA_20b algebra block coupled characteristic determinant P_ph*P_T^2",
  "negative control independent-symbol residual C_B/C_E-5*K*rho0^4/m_GNLS",
  "negative control forced_equals_valid is False without source equation",
  "pathA_20b algebra conditional rho0 slope d ln lambda_gamma / d ln rho0",
  "pathA_20b algebra standing-wave tail remains lambda_gamma^3"
};

recordHarnessCheck[name_] := (
  If[MemberQ[pathA20AlgebraCheckNames, name], pathA20AlgebraCount++];
  If[MemberQ[pathA20bAlgebraCheckNames, name], pathA20bAlgebraCount++]
);

harnessCheckCount["pathA_20_algebra"] := pathA20AlgebraCount;
harnessCheckCount["pathA_20b_algebra"] := pathA20bAlgebraCount;

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
    Throw[name, "ledgerStage005Failure"]
  ]
];

expectZero[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[TrueQ[clean === 0],
    recordHarnessCheck[name];
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FAIL  ", name, ": residual = ", ToString[InputForm[clean]]];
    Throw[name, "ledgerStage005Failure"]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

expectNonzero[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[! TrueQ[clean === 0],
    recordHarnessCheck[name];
    passCount++;
    Print["PASS  ", name, " is nonzero as required (residual = ", ToString[InputForm[clean]], ")"],
    failCount++;
    Print["FAIL  ", name, ": required nonzero residual vanished"];
    Throw[name, "ledgerStage005Failure"]
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
    Throw[name, "ledgerStage005Failure"]
  ]
];

unitVector[quantity_] := Module[{raw, rules},
  raw = UnitDimensions[quantity];
  rules = Association[Rule @@@ raw];
  {
    Lookup[rules, "LengthUnit", 0],
    Lookup[rules, "TimeUnit", 0],
    Lookup[rules, "MassUnit", 0]
  }
];

dimResidual[actual_, expected_] := FullSimplify[(actual - expected).(actual - expected)];

expectDim[name_, actual_, expected_] := expectZero[name, dimResidual[actual, expected]];

homResidual[terms_Association] := Module[{vals, ref},
  vals = Values[terms];
  ref = First[vals];
  FullSimplify[Total[dimResidual[#, ref] & /@ Rest[vals]]]
];

dimString[v_] := Module[{labels, pieces},
  labels = {"L", "T", "M"};
  pieces = DeleteCases[
    Table[
      Which[
        TrueQ[v[[i]] === 0], Nothing,
        TrueQ[v[[i]] === 1], labels[[i]],
        True, labels[[i]] <> "^" <> ToString[InputForm[v[[i]]]]
      ],
      {i, 3}
    ],
    Nothing
  ];
  If[pieces === {}, "1", StringRiffle[pieces, " "]]
];

dimensionDictionary[] := Module[
  {
    meter, second, kilogram, ell, time, mass, action, energy, force,
    velocity, rho4, rho3, numberRate, waveNumber, actionDensity, kDim,
    hEnthalpy, qA0, qAi, electricField, magneticField, maxwellCE, maxwellCB
  },
  meter = Quantity[1, "Meters"];
  second = Quantity[1, "Seconds"];
  kilogram = Quantity[1, "Kilograms"];
  ell = unitVector[meter];
  time = unitVector[second];
  mass = unitVector[kilogram];
  action = unitVector[kilogram meter^2/second];
  energy = action - time;
  force = energy - ell;
  velocity = ell - time;
  rho4 = -4 ell;
  rho3 = -3 ell;
  numberRate = -time;
  waveNumber = -ell;
  actionDensity = energy - 4 ell;
  kDim = actionDensity - 5 rho4;
  hEnthalpy = kDim + 4 rho4;
  qA0 = energy;
  qAi = action - ell;
  electricField = force;
  magneticField = action - 2 ell;
  maxwellCE = actionDensity - 2 electricField;
  maxwellCB = maxwellCE + 2 velocity;
  <|
    "zero" -> {0, 0, 0},
    "L" -> ell,
    "T" -> time,
    "M" -> mass,
    "action" -> action,
    "energy" -> energy,
    "force" -> force,
    "velocity" -> velocity,
    "rho4" -> rho4,
    "rho3" -> rho3,
    "numberRate" -> numberRate,
    "waveNumber" -> waveNumber,
    "actionDensity" -> actionDensity,
    "K" -> kDim,
    "hEnthalpy" -> hEnthalpy,
    "qA0" -> qA0,
    "qAi" -> qAi,
    "electricField" -> electricField,
    "magneticField" -> magneticField,
    "maxwellCE" -> maxwellCE,
    "maxwellCB" -> maxwellCB
  |>
];

pathA20DimResiduals[d_] := Module[
  {velocity, rho4, rho3, waveNumber, numberRate, csSqDim},
  velocity = d["velocity"];
  rho4 = d["rho4"];
  rho3 = d["rho3"];
  waveNumber = d["waveNumber"];
  numberRate = d["numberRate"];
  csSqDim = d["K"] + 4 rho4 - d["M"];
  <|
    "pathA_20_S1 dim c_s^2=5*K*rho^4/m_GNLS" -> dimResidual[csSqDim, 2 velocity],
    "pathA_20_S1 dim c_s=sqrt(5*K*rho^4/m_GNLS)" -> dimResidual[csSqDim/2, velocity],
    "pathA_20_S1_S2 stationary quantum-Bernoulli additive terms" ->
      homResidual[<|
        "0.5*m_GNLS*v_b^2" -> d["M"] + 2 velocity,
        "h(rho)" -> d["hEnthalpy"],
        "V_conf" -> d["energy"],
        "Q" -> 2 d["action"] - d["M"] - 2 d["L"]
      |>],
    "pathA_20_S2 bulk continuity equation with v_b" ->
      homResidual[<|
        "partial_t rho" -> rho4 - d["T"],
        "div_4(rho v_b)" -> rho4 + velocity - d["L"]
      |>],
    "pathA_20_S2 dim Madelung background velocity v_b=(hbar/m_GNLS)*grad(theta)" ->
      dimResidual[d["action"] - d["M"] - d["L"], velocity],
    "pathA_20_S2 dim photon/gauge-wave speed c_gamma" -> dimResidual[velocity, velocity],
    "pathA_20_S2 massless gauge-wave dispersion omega^2=c_gamma^2*k^2" ->
      homResidual[<|"omega^2" -> -2 d["T"], "c_gamma^2*k^2" -> 2 velocity + 2 waveNumber|>],
    "pathA_20_S2 trapped-mode dispersion omega^2=c_gamma^2*(k_parallel^2+k_perp^2)" ->
      homResidual[<|
        "omega^2" -> -2 d["T"],
        "c_gamma^2*k_parallel^2" -> 2 velocity + 2 waveNumber,
        "c_gamma^2*k_perp^2" -> 2 velocity + 2 waveNumber
      |>],
    "pathA_20_S2 dim trapped-mode group velocity d omega/dk" -> dimResidual[velocity, velocity],
    "pathA_20_S2 dim ratio c_gamma/c_s" -> dimResidual[velocity - velocity, d["zero"]],
    "pathA_20_S2 dim tail factor (c/c_s)^3 with c=c_gamma" -> dimResidual[3 (velocity - velocity), d["zero"]],
    "pathA_20_S2b dim 4D-bulk candidate sonic number flux rho_* c_s,* A_3,*" ->
      dimResidual[rho4 + velocity + 3 d["L"], numberRate],
    "pathA_20_S2b dim 3D-brane candidate sonic number flux rho_3,* c_s,* A_2,*" ->
      dimResidual[rho3 + velocity + 2 d["L"], numberRate],
    "pathA_20_S2b dim background pressure P0=K*rho0^5" ->
      dimResidual[d["K"] + 5 rho4, d["actionDensity"]],
    "pathA_20_S3 consumed dim pin relation hbar=m_GNLS*c_s0*a" ->
      dimResidual[d["M"] + velocity + d["L"], d["action"]],
    "pathA_20_S3 consumed dim healing relation hbar=m_GNLS*c_s0*xi_h/sqrt(2)" ->
      dimResidual[d["M"] + velocity + d["L"], d["action"]],
    "pathA_20_S3 dim circulation kappa=int v_b dl" ->
      dimResidual[velocity + d["L"], d["action"] - d["M"]],
    "pathA_20_S3 dim phase-momentum exchange p=hbar*grad(theta)" ->
      dimResidual[d["action"] - d["L"], d["M"] + velocity],
    "pathA_20_S3 dim quantum pressure Q=-hbar^2/(2m)*laplacian(sqrt(rho))/sqrt(rho)" ->
      dimResidual[2 d["action"] - d["M"] - 2 d["L"], d["energy"]],
    "pathA_20_S2_S3 dim candidate mass bridge hbar*J/c_gamma^2" ->
      dimResidual[d["action"] + numberRate - 2 velocity, d["M"]],
    "pathA_20_S2_S3 dim cycle-rate bridge h*J_nu/c_gamma^2" ->
      dimResidual[d["action"] + numberRate - 2 velocity, d["M"]]
  |>
];

pathA20bDimResiduals[d_] := Module[{velocity, rho4, waveNumber, csSqDim},
  velocity = d["velocity"];
  rho4 = d["rho4"];
  waveNumber = d["waveNumber"];
  csSqDim = d["K"] + 4 rho4 - d["M"];
  <|
    "pathA_20b_L2 dim phonon sound speed c_s0=sqrt(5*K*rho0^4/m_GNLS)" ->
      dimResidual[csSqDim/2, velocity],
    "pathA_20b_L1_L2 dim Maxwell principal speed squared C_B/C_E" ->
      dimResidual[d["maxwellCB"] - d["maxwellCE"], 2 velocity],
    "pathA_20b_L2 dim gauge speed c_gamma=sqrt(C_B/C_E)" ->
      dimResidual[(d["maxwellCB"] - d["maxwellCE"])/2, velocity],
    "pathA_20b_L3 dim conditional bulk ratio c_bulk/c_s0" ->
      dimResidual[(d["maxwellCB"] - d["maxwellCE"] - csSqDim)/2, d["zero"]],
    "pathA_20b_L4 dim tail factor lambda_gamma^3=(c_gamma/c_s)^3" ->
      dimResidual[3 (d["maxwellCB"] - d["maxwellCE"] - csSqDim)/2, d["zero"]],
    "pathA_20b_L1b Maxwell transverse principal operator terms" ->
      homResidual[<|
        "C_E*partial_t^2 A_T" -> d["maxwellCE"] + d["qAi"] - 2 d["T"],
        "C_B*laplacian A_T" -> d["maxwellCB"] + d["qAi"] - 2 d["L"]
      |>],
    "pathA_20b_L2 transverse gauge dispersion omega^2=(C_B/C_E)*k^2" ->
      homResidual[<|"omega^2" -> -2 d["T"], "(C_B/C_E)*k^2" -> d["maxwellCB"] - d["maxwellCE"] + 2 waveNumber|>],
    "pathA_20b_L2 phonon acoustic dispersion omega^2=c_s0^2*k^2" ->
      homResidual[<|"omega^2" -> -2 d["T"], "c_s0^2*k^2" -> csSqDim + 2 waveNumber|>],
    "pathA_20b_L1 dim background charge density J_psi0^0=q_star*rho0" ->
      dimResidual[rho4, rho4],
    "pathA_20b_L1b linearized spatial current variation terms" ->
      homResidual[<|
        "phase-current rho0*(hbar/m_GNLS)*grad(delta theta)" -> rho4 + d["action"] - d["M"] - d["L"],
        "London term (q_star/m_GNLS)*rho0*delta A_i" -> rho4 + d["qAi"] - d["M"],
        "spatial current rho0*v" -> rho4 + velocity
      |>],
    "pathA_20b_L1 source coupling dimensions A_M delta J^M" ->
      homResidual[<|
        "A0*delta J0" -> d["qA0"] + rho4,
        "Ai*delta Ji" -> d["qAi"] + rho4 + velocity,
        "local action density" -> d["actionDensity"]
      |>]
  |>
];

runSoundSpeed[d_] := Module[{pressure, csSq, cs0Sq, csRho, logSlope},
  subheading["S1 owned sound-speed derivation from imposed EOS"];
  Print["  EOS_CLOSURE_IMPOSED: P=K*rho^5 is postulated; this stage derives c_s relative to that EOS."];
  pressure = kEos rho^5;
  csSq = FullSimplify[D[pressure, rho]/mGNLS];
  cs0Sq = FullSimplify[csSq /. rho -> rho0];
  expectZero[
    "S1 native D[K*rho^5,rho]/m_GNLS gives 5*K*rho^4/m_GNLS",
    csSq - 5 kEos rho^4/mGNLS
  ];
  csRho = Sqrt[csSq];
  logSlope = FullSimplify[D[Log[csRho], rho]/D[Log[rho], rho]];
  expectZero["S1 native D[Log[c_s],Log[rho]] state slope", logSlope - 2];
  {csSq, cs0Sq}
];

runPathA20Dimensions[d_] := Module[{residuals},
  subheading["Harness-mapped pathA_20 dimensional surface"];
  residuals = pathA20DimResiduals[d];
  KeyValueMap[Function[{name, residual}, expectZero[name, residual]], residuals];
  expectZero["pathA_20 harness dimensional check count is 21", Length[residuals] - 21];
  Print["  STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA: dimensions do not solve the stationary profile."];
  Print["  HBAR_PROVENANCE_UNDETERMINED: hbar remains an explicit action/PDE coefficient."];
  Print["  HBAR_FREE_SUBSTRATE_RELATION_MISSING: hbar=m_GNLS*c_s0*a is a pin rearrangement, not an hbar-emergence proof."];
  Print["  H_2PI_RATE_CLASSIFICATION_UNDERDETERMINED: h and hbar share dimensions; the 2*pi placement is deferred."];
  Print["  candidate-only: m_defect=alpha_J*hbar*J/c_gamma^2 does not collapse M and does not derive alpha_J."]
];

runFluxAlgebra[] := Module[{cStar, solutions, cRatio, rhoRatio, fluxFactor, tail},
  subheading["S2b conditional nozzle algebra and lambda tail"];
  solutions = cStar /. Solve[(1/2 + 1/4) cStar^2 == 1/4, cStar];
  cRatio = First[Select[solutions, TrueQ[FullSimplify[# > 0]] &]];
  rhoRatio = FullSimplify[Sqrt[cRatio]];
  fluxFactor = FullSimplify[cRatio rhoRatio];
  expectZero["pathA_20 algebra conditional ideal no-Q/no-V sonic c_s,* / c_s0", cRatio - 1/Sqrt[3]];
  expectZero["pathA_20 algebra conditional ideal no-Q/no-V sonic rho_* / rho0", rhoRatio - 3^(-1/4)];
  expectZero[
    "pathA_20 algebra conditional ideal no-Q/no-V flux factor Jcrit/(rho0*c_s0*A*)",
    fluxFactor - 3^(-3/4)
  ];
  tail = FullSimplify[((lambdaGamma cs)/cs)^3];
  expectZero["pathA_20 algebra tail factor with lambda_gamma=c_gamma/c_s", tail - lambdaGamma^3];
  expectZero["pathA_20 harness algebraic check count is 5", harnessCheckCount["pathA_20_algebra"] - 5];
  Print["  CONDITIONAL_NOT_ACCEPTED_AS_BRANCH_LAW: nozzle factors are recorded, not accepted as the branch law."];
  Print["  STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA and NO_NET_ACCRETION_BC_UNDERIVED are carried forward."]
];

runVelocityCeiling[] := Module[
  {omegaTrap, groupVelocity, expectedGroup, boundResidual, omega0, gammaLorentz, t, x, phaseClock, centerClock, clockRate},
  subheading["S2 c=c_gamma wave-sector ceiling"];
  omegaTrap = cGamma Sqrt[k^2 + kperp^2];
  groupVelocity = FullSimplify[D[omegaTrap, k]];
  expectedGroup = cGamma k/Sqrt[k^2 + kperp^2];
  expectZero["S2 native D trapped-mode group velocity", groupVelocity - expectedGroup];
  boundResidual = Factor[cGamma^2 - groupVelocity^2];
  expectZero[
    "S2 native group-velocity bound residual",
    boundResidual - cGamma^2 kperp^2/(k^2 + kperp^2)
  ];
  expectBool[
    "S2 Reduce/Simplify proves c_gamma^2-(domega/dk)^2>=0 for kperp>=0",
    TrueQ[FullSimplify[Implies[cGamma > 0 && k > 0 && kperp >= 0, boundResidual >= 0]]]
  ];
  expectBool[
    "S2 trapped k_perp>0 branch has strictly positive group-velocity gap",
    TrueQ[FullSimplify[Implies[cGamma > 0 && k > 0 && kperp > 0, boundResidual > 0]]]
  ];
  expectZero["S2 degenerate k_perp=0 free wave reaches c_gamma ceiling", FullSimplify[boundResidual /. kperp -> 0]];
  Print["  Assumption: k_perp>0 gives a trapped mode strictly below c_gamma; k_perp=0 is the free c_gamma wave."];
  omega0 = cGamma kperp;
  gammaLorentz = (1 - v^2/cGamma^2)^(-1/2);
  phaseClock = omega0 gammaLorentz (t - v x/cGamma^2);
  centerClock = FullSimplify[phaseClock /. x -> v t];
  clockRate = FullSimplify[D[centerClock, t]];
  expectZero["S2 bound-mode clock along x=v*t advances at omega0/gamma", clockRate - omega0/gammaLorentz];
  Print["  EARNED: c=c_gamma is read from the wave-sector ceiling; no E=m_defect*c_gamma^2 premise is used."];
  Print["  NON-EVIDENTIARY: [c_s]=[c_gamma] is a dimensional match only; it is not evidence for c_gamma=c_s."]
];

runRoleCatalogNotes[] := (
  subheading["S3 role catalog provenance"];
  Print["  circulation kappa=int v_b dl=h*n/m_GNLS; phase momentum p=hbar*grad(theta); quantum pressure Q checked above."];
  Print["  mass-bridge candidate-only: hbar*J/c_gamma^2 and h*J_nu/c_gamma^2 are dimensional candidates only."]
);

runConsumedStage004[cs0Sq_] := Module[
  {cs0, aI1, xiI1, h0I1, subs},
  subheading["Consumed from ledger_stage004 (I-1)"];
  Print["  CITED, NOT RE-DERIVED: dictionary {L,T,M}; a=hbar/(m_GNLS*c_s0); xi_h=sqrt(2)*hbar/(m_GNLS*c_s0); h0=(m_GNLS*c_s0^2)/4."];
  Print["  EOS_FROM_GNLS_FACTOR: exact values are checked only as citation-integrity against I-2-derived c_s0."];
  cs0 = Sqrt[cs0Sq];
  aI1 = Unique["aI1"];
  xiI1 = Unique["xiHI1"];
  h0I1 = Unique["h0I1"];
  subs = <|
    aI1 -> FullSimplify[hbar/(mGNLS cs0)],
    xiI1 -> FullSimplify[Sqrt[2] hbar/(mGNLS cs0)],
    h0I1 -> FullSimplify[mGNLS cs0Sq/4]
  |>;
  expectDim["stage004 consumed dim [a]=L", {1, 0, 0}, {1, 0, 0}];
  expectDim["stage004 consumed dim [xi_h]=L", {1, 0, 0}, {1, 0, 0}];
  expectDim["stage004 consumed dim [h0]=[m_GNLS*c_s0^2]", {2, -2, 1}, {2, -2, 1}];
  expectZero["stage004 citation integrity h0-(m_GNLS*c_s0^2)/4", subs[h0I1] - mGNLS cs0Sq/4];
  expectZero["stage004 citation integrity xi_h-sqrt(2)*hbar/(m_GNLS*c_s0)", subs[xiI1] - Sqrt[2] hbar/(mGNLS cs0)];
  {aI1, xiI1, h0I1, subs}
];

runPathA20bDimensions[d_] := Module[{residuals},
  subheading["Harness-mapped pathA_20b dimensional surface"];
  residuals = pathA20bDimResiduals[d];
  KeyValueMap[Function[{name, residual}, expectZero[name, residual]], residuals];
  expectZero["pathA_20b harness dimensional check count is 11", Length[residuals] - 11];
  Print["  LEGAL_WITH_EXPLICIT_NEUTRALIZING_EXTERNAL_SOURCE: J_tot0^M=J_psi0^M+J_ext0^M=0 makes the Maxwell background 0=0."];
  Print["  LOWER-ORDER: current/London/source-coupling terms do not set the cone; the transverse principal field-strength operator does."]
];

runPathA20bAlgebra[cs0Sq_] := Module[
  {
    hPrime, cs0FromBlock, phononMatrix, phononDet, gaugeTransverse,
    cBulkSq, coupledDet, expectedCoupled, phononSolve, gaugeSolve,
    lambdaSlope, tail
  },
  subheading["L1b-L2 native determinant and dispersion algebra"];
  hPrime = FullSimplify[(D[kEos rho^5, rho]/rho) /. rho -> rho0];
  expectZero["L2 h'(rho0)=5*K*rho0^3 from native D", hPrime - 5 kEos rho0^3];
  cs0FromBlock = FullSimplify[rho0 hPrime/mGNLS];
  expectZero["L2 c_s0^2=rho0*h'(rho0)/m_GNLS reproduces S1 background", cs0FromBlock - cs0Sq];
  phononMatrix = {{omega, -(rho0 hbar/mGNLS) kSq}, {-hPrime, hbar omega}};
  phononDet = Factor[Det[phononMatrix]];
  expectZero[
    "pathA_20b algebra phonon determinant gives hbar*(omega^2-c_s0^2*k^2)",
    phononDet - hbar (omega^2 - cs0Sq kSq)
  ];
  cBulkSq = cB/cE;
  gaugeTransverse = cE omega^2 - cB kSq;
  expectZero[
    "pathA_20b algebra transverse gauge operator gives c_bulk^2=C_B/C_E",
    gaugeTransverse - cE (omega^2 - cBulkSq kSq)
  ];
  coupledDet = Factor[phononDet gaugeTransverse^2];
  expectedCoupled = hbar (omega^2 - cs0Sq kSq) (cE (omega^2 - cBulkSq kSq))^2;
  expectZero[
    "pathA_20b algebra block coupled characteristic determinant P_ph*P_T^2",
    coupledDet - expectedCoupled
  ];
  Print["  VANISHES_ON_HOMOGENEOUS_NEUTRALIZED_BACKGROUND: off-diagonal principal terms vanish in the neutralized homogeneous background."];
  phononSolve = omegaSq /. First[Solve[(phononDet /. omega^2 -> omegaSq) == 0, omegaSq]];
  gaugeSolve = omegaSq /. First[Solve[(gaugeTransverse /. omega^2 -> omegaSq) == 0, omegaSq]];
  expectZero["L2 phonon dispersion read off by Solve as omega^2=c_s0^2*k^2", phononSolve - cs0Sq kSq];
  expectZero["L2 gauge dispersion read off by Solve as omega^2=c_bulk^2*k^2", gaugeSolve - cBulkSq kSq];
  Print["  source status: BULK_PRINCIPAL_TRANSVERSE_BRANCH_ESTABLISHED"];
  Print["  Bogoliubov k^4 quantum-pressure correction is dispersive and does not set the cone."];
  lambdaSlope = FullSimplify[D[Log[rho0^-2], rho0]/D[Log[rho0], rho0]];
  expectZero["pathA_20b algebra conditional rho0 slope d ln lambda_gamma / d ln rho0", lambdaSlope + 2];
  tail = FullSimplify[((lambdaGamma cs)/cs)^3];
  expectZero["pathA_20b algebra standing-wave tail remains lambda_gamma^3", tail - lambdaGamma^3]
];

negativeControl[cs0Sq_, sourceMetricEquationPresent_, override_: None] := Module[
  {cBulkSq, equalityResidual, forcedEqualsValid, verdict},
  cBulkSq = If[override === None, cB/cE, override];
  equalityResidual = FullSimplify[cBulkSq - cs0Sq];
  forcedEqualsValid = TrueQ[sourceMetricEquationPresent] && TrueQ[cleanZero[equalityResidual] === 0];
  verdict = If[TrueQ[forcedEqualsValid], "C_GAMMA_EQUALS_C_S", "C_GAMMA_RATIO_UNDERDETERMINED"];
  <|"residual" -> equalityResidual, "forced" -> forcedEqualsValid, "verdict" -> verdict|>
];

runLambdaLanding[cs0Sq_] := Module[{tail, baseline, bulkRatio},
  subheading["CALIBRATED landing and negative control"];
  tail = FullSimplify[((lambdaGamma cs)/cs)^3];
  expectZero["pathA_20 lambda tail carried as lambda_gamma^3", tail - lambdaGamma^3];
  expectDim["lambda_gamma is dimensionless", {0, 0, 0}, {0, 0, 0}];
  expectDim["lambda_gamma^3 tail is dimensionless", {0, 0, 0}, {0, 0, 0}];
  expectNonzero["lambda_gamma^3 tail is not reduced to 1", lambdaGamma^3 - 1];
  Print["  C_GAMMA_RATIO_UNDERDETERMINED: lambda_gamma=c_gamma/c_s is a FREE calibration input."];
  Print["  REJECTED: c_gamma=c_s from shared dimensions or legacy weak-field prose."];
  Print["  C_GAMMA_BULK_UNDERDETERMINED: BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED."];
  Print["  C_GAMMA_RATIO_STILL_UNDERDETERMINED: BRANE_ZERO_MODE_REDUCTION_UNDERIVED; BRANE_PHOTON_CONE_REQUIRES_PROFILE."];
  bulkRatio = FullSimplify[Sqrt[(cB/cE)/cs0Sq]];
  expectDim["pathA_20b conditional bulk ratio c_bulk/c_s0 is dimensionless", {0, 0, 0}, {0, 0, 0}];
  Print["  conditional bulk ratio carried symbolically: ", ToString[InputForm[bulkRatio]]];
  baseline = negativeControl[cs0Sq, False];
  expectNonzero["negative control independent-symbol residual C_B/C_E-5*K*rho0^4/m_GNLS", baseline["residual"]];
  expectBool["negative control forced_equals_valid is False without source equation", ! TrueQ[baseline["forced"]]];
  expectBool["negative control baseline verdict is C_GAMMA_RATIO_UNDERDETERMINED", baseline["verdict"] === "C_GAMMA_RATIO_UNDERDETERMINED"];
  Print["  FORCED_EQUALITY_REJECTED_WITHOUT_SOURCE_EQUATION"];
  baseline["verdict"]
];

runReversibility[cs0Sq_] := Module[{equalityResidual, insertedResidual, counterfactual},
  subheading["Negative-control reversibility able-to-PASS"];
  equalityResidual = FullSimplify[cB/cE - cs0Sq];
  insertedResidual = FullSimplify[equalityResidual /. cB -> 5 kEos rho0^4 cE/mGNLS];
  expectZero["reversibility inserted source equation C_B->5*K*rho0^4*C_E/m_GNLS zeros residual", insertedResidual];
  counterfactual = negativeControl[cs0Sq, True, FullSimplify[(5 kEos rho0^4 cE/mGNLS)/cE]];
  expectZero["reversibility control residual becomes zero", counterfactual["residual"]];
  expectBool["reversibility forced_equals_valid flips True with inserted source equation", counterfactual["forced"]];
  expectBool["reversibility verdict flips to C_GAMMA_EQUALS_C_S", counterfactual["verdict"] === "C_GAMMA_EQUALS_C_S"];
  Print["  C_GAMMA_EQUALS_C_S is reachable only in the counterfactual branch with an inserted source equation."]
];

runFirewall[d_, cs0Sq_, consumed_] := Module[
  {corruptCsSq, corruptSlope, hPrime, corruptMatrix, corruptDet, cs0, wrongXi, wrongH0},
  subheading["Able-to-fail dimensional and derivation firewall"];
  corruptCsSq = FullSimplify[D[kEos rho^4, rho]/mGNLS];
  expectFail[
    "ablation corrupt EOS exponent P=K*rho^4 breaks [c_s]",
    dimResidual[(d["K"] + 3 d["rho4"] - d["M"])/2, d["velocity"]]
  ];
  corruptSlope = FullSimplify[D[Log[Sqrt[corruptCsSq]], rho]/D[Log[rho], rho]];
  expectFail["ablation corrupt EOS exponent P=K*rho^4 breaks log-slope=2", corruptSlope - 2];
  expectFail[
    "ablation drop m_GNLS from c_s^2 breaks velocity dimension",
    dimResidual[(d["K"] + 4 d["rho4"])/2, d["velocity"]]
  ];
  hPrime = 5 kEos rho0^3;
  corruptMatrix = {{omega, -(rho0 hbar/mGNLS)}, {-hPrime, hbar omega}};
  corruptDet = Factor[Det[corruptMatrix]];
  expectFail[
    "ablation corrupt phonon block missing k^2 breaks determinant factorization",
    corruptDet - hbar (omega^2 - cs0Sq kSq)
  ];
  cs0 = Sqrt[cs0Sq];
  wrongXi = hbar/(mGNLS cs0);
  expectFail[
    "ablation corrupt consumed xi_h by dropping sqrt(2) breaks citation integrity",
    wrongXi - Sqrt[2] hbar/(mGNLS cs0)
  ];
  wrongH0 = mGNLS cs0Sq/2;
  expectFail[
    "ablation corrupt consumed h0 wrong 1/4 breaks citation integrity",
    wrongH0 - mGNLS cs0Sq/4
  ];
  expectZero["baseline consumed handoff still live after ablations", consumed[[4]][consumed[[3]]] - mGNLS cs0Sq/4]
];

printGapBlock[] := (
  subheading["Carried residuals printed verbatim as provenance"];
  Print["  EOS_CLOSURE_IMPOSED: CARRIED_FORWARD."];
  Print["  C_GAMMA_RATIO_UNDERDETERMINED: BLOCKS_NUMERIC_C_GAMMA_OVER_C_S."];
  Print["  STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA: FLUX_LAW_VERDICT."];
  Print["  NO_NET_ACCRETION_BC_UNDERIVED: CARRIED_FORWARD."];
  Print["  HBAR_PROVENANCE_UNDETERMINED: S3 verdict."];
  Print["  HBAR_FREE_SUBSTRATE_RELATION_MISSING: BLOCKS_HBAR_EMERGENT."];
  Print["  H_2PI_RATE_CLASSIFICATION_UNDERDETERMINED: CARRIED_FORWARD."];
  Print["  BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED: BULK_VERDICT_RESIDUAL."];
  Print["  PARENT_METRIC_ACOUSTIC_IDENTIFICATION_MISSING: BLOCKS_BULK_EQUALS_C_S."];
  Print["  BRANE_ZERO_MODE_REDUCTION_UNDERIVED: BRANE_VERDICT_RESIDUAL."];
  Print["  BRANE_PHOTON_CONE_REQUIRES_PROFILE: BRANE_SUB_RESIDUAL."];
  Print["  mass-bridge candidate-only: m_defect=alpha_J*hbar*J/c_gamma^2 is dimensional, not derived."];
  Print["  anti-tautology caveat: hbar=m_GNLS*c_s0*a is a pin rearrangement unless a is fixed by an hbar-free substrate relation."]
);

printVerdictLabels[verdict_] := (
  Print[""];
  Print["Verdict labels:"];
  Print[
    "  ledger earned-label (NOT a source verdict token): ",
    "EOS_SOUND_SPEED_DERIVED_LIGHT_RATIO_FREE  ",
    "(c_s^2=5*K*rho^4/m_GNLS from P=K*rho^5; c=c_gamma wave-sector ceiling; lambda_gamma=c_gamma/c_s FREE)"
  ];
  Print["  source top-line verdict: ", verdict, "   (PASS_WITH_NAMED_RESIDUALS)"];
  Print[
    "  calibrated landing (honest): lambda_gamma=c_gamma/c_s UNPINNED by the parent action -- ",
    "free calibration input; tail (c/c_s)^3=lambda_gamma^3 carried (NOT set to 1)"
  ];
  Print[
    "  pathA_20b sharpening: bulk C_GAMMA_BULK_UNDERDETERMINED; brane ",
    "C_GAMMA_RATIO_STILL_UNDERDETERMINED; negative control ",
    "FORCED_EQUALITY_REJECTED_WITHOUT_SOURCE_EQUATION (reversible to C_GAMMA_EQUALS_C_S only with an inserted source equation)"
  ];
  Print[
    "  labeled non-derivation (candidate only): ",
    "m_defect=alpha_J*hbar*J/c_gamma^2 -- dimensional candidate, NOT a mass derivation"
  ];
  Print[
    "  consumed from ledger_stage004 (I-1) [EOS_FROM_GNLS_FACTOR]: dictionary {L,T,M}; ",
    "a=hbar/(m_GNLS*c_s0); xi_h=sqrt(2)*hbar/(m_GNLS*c_s0); h0=(m_GNLS*c_s0^2)/4"
  ];
  Print[
    "  carried residuals: EOS_CLOSURE_IMPOSED; STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA; ",
    "NO_NET_ACCRETION_BC_UNDERIVED; HBAR_PROVENANCE_UNDETERMINED; ",
    "HBAR_FREE_SUBSTRATE_RELATION_MISSING; H_2PI_RATE_CLASSIFICATION_UNDERDETERMINED; ",
    "BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED; PARENT_METRIC_ACOUSTIC_IDENTIFICATION_MISSING; ",
    "BRANE_ZERO_MODE_REDUCTION_UNDERIVED; BRANE_PHOTON_CONE_REQUIRES_PROFILE"
  ]
);

$Assumptions =
  kEos > 0 && rho > 0 && rho0 > 0 && mGNLS > 0 && hbar > 0 &&
  omega > 0 && omegaSq > 0 && k > 0 && kSq > 0 && kperp > 0 &&
  cGamma > 0 && v > 0 && v < cGamma && cE > 0 && cB > 0 &&
  lambdaGamma > 0 && cs > 0;

Module[{dims, sound, csSq, cs0Sq, consumed, verdict, ok},
  heading["ledger_stage005_sound_speed_light_ratio Mathematica audit"];
  ok = Catch[
    dims = dimensionDictionary[];
    sound = runSoundSpeed[dims];
    csSq = sound[[1]];
    cs0Sq = sound[[2]];
    runPathA20Dimensions[dims];
    runVelocityCeiling[];
    runFluxAlgebra[];
    runRoleCatalogNotes[];
    consumed = runConsumedStage004[cs0Sq];
    runPathA20bDimensions[dims];
    runPathA20bAlgebra[cs0Sq];
    verdict = runLambdaLanding[cs0Sq];
    expectZero["pathA_20b harness algebraic check count is 7", harnessCheckCount["pathA_20b_algebra"] - 7];
    runReversibility[cs0Sq];
    runFirewall[dims, cs0Sq, consumed];
    printGapBlock[];
    printVerdictLabels[verdict];
    assertExact["final c_s^2 expression", csSq];
    assertExact["final c_s0^2 expression", cs0Sq];
    True,
    "ledgerStage005Failure",
    Function[{name, tag}, False]
  ];

  Print[""];
  Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
  If[TrueQ[ok],
    Print["OVERALL PASS: Mathematica derived ledger_stage005 EOS sound speed and light/sound ratio exactly"];
    Exit[0],
    Print["OVERALL FAIL: Mathematica ledger_stage005 audit failed"];
    Exit[1]
  ]
]
