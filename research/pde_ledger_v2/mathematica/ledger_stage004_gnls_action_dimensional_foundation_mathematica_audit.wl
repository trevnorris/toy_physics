(* Ledger stage004 Mathematica audit: GNLS action dimensional foundation.

   Print-only, standalone, no arguments, no exports.  This route uses native
   Quantity/UnitDimensions to obtain primitive exponent triples, Solve for the
   derived action dictionary and healing scale, and NullSpace/RowReduce for
   the pin relation.
*)

ClearAll[
  heading, subheading, cleanZero, assertExact, expectZero, expectBool,
  expectNonzero, expectFail, unitVector, dimResidual, expectDim, dimString,
  homResidual, factorToReach, dropM, primitiveDimensions, deriveDictionary,
  runTwoTierDictionary, runSymbolicCore, normalizeNullVector,
  runPinAnalysis, foundationResiduals, ltRepresentationResiduals,
  flaggedResiduals, allHarnessPassQ, ltRejectionGate, verdictFromPredicate,
  runHarnessChecks, runFlaggedResiduals, runVerdict, runFirewall,
  printGapBlock, printVerdictLabels, passCount, failCount
];

passCount = 0;
failCount = 0;

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
    Throw[name, "ledgerStage004Failure"]
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
    Throw[name, "ledgerStage004Failure"]
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
    Throw[name, "ledgerStage004Failure"]
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
    Throw[name, "ledgerStage004Failure"]
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

homResidual[terms_Association] := Module[{vals, ref},
  vals = Values[terms];
  ref = First[vals];
  FullSimplify[Total[dimResidual[#, ref] & /@ Rest[vals]]]
];

factorToReach[expected_, actual_] := expected - actual;

dropM[v_] := {v[[1]], v[[2]], 0};

primitiveDimensions[] := Module[{meter, second, kilogram},
  meter = Quantity[1, "Meters"];
  second = Quantity[1, "Seconds"];
  kilogram = Quantity[1, "Kilograms"];
  <|
    "hbar" -> unitVector[kilogram meter^2/second],
    "mGNLS" -> unitVector[kilogram],
    "rho0" -> unitVector[meter^-4]
  |>
];

deriveDictionary[] := Module[
  {
    prim, zero, ell, time, mass, hbar, mGNLS, rho0, action, energy,
    velocity, actionDensity, pL, pT, pM, psi, rho, kEos, cs0, h0,
    xiH, force, electricField, magneticField, maxwellCoeff
  },

  prim = primitiveDimensions[];
  zero = {0, 0, 0};
  ell = unitVector[Quantity[1, "Meters"]];
  time = unitVector[Quantity[1, "Seconds"]];
  mass = unitVector[Quantity[1, "Kilograms"]];
  hbar = prim["hbar"];
  mGNLS = prim["mGNLS"];
  rho0 = prim["rho0"];

  action = hbar;
  energy = action - time;
  velocity = ell - time;
  actionDensity = energy - 4 ell;

  psi = {pL, pT, pM} /. First[Solve[
      Thread[hbar - time + 2 {pL, pT, pM} == actionDensity],
      {pL, pT, pM}
    ]];
  rho = 2 psi;
  kEos = actionDensity - 5 rho;
  cs0 = (kEos + 4 rho0 - mGNLS)/2;
  h0 = kEos + 4 rho0;
  xiH = (2 hbar - mGNLS - h0)/2;
  force = energy - ell;
  electricField = force;
  magneticField = action - 2 ell;
  maxwellCoeff = actionDensity - 2 electricField;

  <|
    "zero" -> zero,
    "L" -> ell,
    "T" -> time,
    "M" -> mass,
    "hbar" -> hbar,
    "mGNLS" -> mGNLS,
    "rho0" -> rho0,
    "action" -> action,
    "energy" -> energy,
    "force" -> force,
    "velocity" -> velocity,
    "psi" -> psi,
    "rho" -> rho,
    "rho3" -> -3 ell,
    "K" -> kEos,
    "cS0" -> cs0,
    "h0" -> h0,
    "xiH" -> xiH,
    "lagrangianDensity" -> actionDensity,
    "qA0" -> energy,
    "qAi" -> action - ell,
    "electricField" -> electricField,
    "magneticField" -> magneticField,
    "maxwellCoeff" -> maxwellCoeff,
    "muWall" -> force - 2 velocity,
    "Tw" -> force,
    "USigmaRR" -> force - 2 ell,
    "G3" -> {3, -2, -1},
    "G4" -> {4, -2, -1}
  |>
];

runTwoTierDictionary[d_] := (
  subheading["Two-tier dictionary derivation"];
  Print["Primitive inputs (posted through UnitDimensions):"];
  Print["  [hbar] = ", dimString[d["hbar"]], "  PRIMITIVE INPUT"];
  Print["  [m_GNLS] = ", dimString[d["mGNLS"]], "  PRIMITIVE INPUT"];
  Print["  [rho0] = ", dimString[d["rho0"]], "  PRIMITIVE INPUT"];
  Print["Derived-by-composition dimensions:"];
  Print["  [psi] = ", dimString[d["psi"]], "  from kinetic density Solve"];
  Print["  [rho] = ", dimString[d["rho"]], "  from psi^2"];
  Print["  [K] = ", dimString[d["K"]], "  from EOS density"];
  Print["  [c_s0] = ", dimString[d["cS0"]], "  from 5*K*rho0^4/m_GNLS"];
  Print["  [h0] = ", dimString[d["h0"]], "  from (5K/4)*rho0^4"];
  Print["  [xi_h] = ", dimString[d["xiH"]], "  from core balance"];

  expectDim["primitive hbar used as action in kinetic density", d["hbar"], d["action"]];
  expectDim["primitive m_GNLS used in Madelung velocity hbar/(m*L)", d["hbar"] - d["mGNLS"] - d["L"], d["velocity"]];
  expectDim["primitive rho0 checked against derived rho=psi^2", d["rho"], d["rho0"]];
  expectDim["[psi] derived from action-density equation matches sqrt(rho0)", d["psi"], d["rho0"]/2];
  expectDim["[K] derived from EOS matches m*c_s0^2/rho0^4 target", d["K"], d["mGNLS"] + 2 d["velocity"] - 4 d["rho0"]];
  expectDim["[c_s0] derived from sound-speed composition", d["cS0"], d["velocity"]];
  expectDim["[h0] EOS scale matches m_GNLS*c_s0^2", d["h0"], d["mGNLS"] + 2 d["cS0"]];
  expectDim["[xi_h] core-balance dimension", d["xiH"], d["L"]]
);

runSymbolicCore[] := Module[
  {hbar, m, kEos, rho0, cs0, h0, xi, csSq, kRule, h0Core, xiRules, xiCore},
  subheading["Sound-speed and healing algebra"];
  $Assumptions = hbar > 0 && m > 0 && kEos > 0 && rho0 > 0 && cs0 > 0 && h0 > 0 && xi > 0;
  csSq = cs02 /. First[Solve[cs02 == 5 kEos rho0^4/m, cs02]];
  expectZero["algebraic sound-speed law c_s0^2=5*K*rho0^4/m_GNLS", csSq - 5 kEos rho0^4/m];
  kRule = First[Solve[cs0^2 == 5 kEos rho0^4/m, kEos]];
  h0Core = FullSimplify[(5/4) kEos rho0^4 /. kRule];
  expectZero["GNLS core balance h0=(5K/4)*rho0^4=(m_GNLS*c_s0^2)/4", h0Core - m cs0^2/4];
  xiRules = Solve[xi^2 == hbar^2/(2 m h0), xi];
  xiCore = FullSimplify[(xi /. Last[xiRules]) /. h0 -> m cs0^2/4];
  expectZero["healing length xi_h=sqrt(2)*hbar/(m_GNLS*c_s0)", xiCore - Sqrt[2] hbar/(m cs0)]
];

normalizeNullVector[v_] := Module[{denoms, scale, ints, gcd},
  denoms = Denominator /@ v;
  scale = Apply[LCM, denoms];
  ints = scale v;
  gcd = Apply[GCD, Abs[ints]];
  ints = ints/gcd;
  If[First[ints] < 0, -ints, ints]
];

runPinAnalysis[d_] := Module[
  {pinDims, mat, rowReduced, rank, nullity, relation, relationDim, a, cs0, hbar, m, aRule, xi},
  subheading["Pin null-relation"];
  pinDims = {d["L"], d["cS0"], d["hbar"], d["mGNLS"]};
  mat = Transpose[pinDims];
  rowReduced = RowReduce[mat];
  rank = MatrixRank[mat];
  nullity = Length[pinDims] - rank;
  relation = normalizeNullVector[First[NullSpace[mat]]];
  relationDim = relation.pinDims;
  expectZero["pin matrix rank is 3", rank - 3];
  expectZero["four pins on three bases leave nullity 1", nullity - 1];
  expectZero["pin RowReduce has three pivot rows", MatrixRank[rowReduced] - 3];
  expectZero["pin null vector is a*c_s0*hbar^-1*m_GNLS", Total[(relation - {1, 1, -1, 1})^2]];
  expectDim["derived pin relation a*c_s0*m_GNLS/hbar is dimensionless", relationDim, d["zero"]];

  $Assumptions = a > 0 && cs0 > 0 && hbar > 0 && m > 0;
  aRule = a /. First[Solve[a cs0 m/hbar == 1, a]];
  expectZero["derived pin relation a=hbar/(m_GNLS*c_s0)", aRule - hbar/(m cs0)];
  xi = Sqrt[2] hbar/(m cs0);
  expectZero["raw four pins give a/xi_h=1/sqrt(2)", aRule/xi - 1/Sqrt[2]];
  Print["  A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT: a is a branch collective moment, not a base invariant."]
];

(* Harness mapping: these 14 entries preserve _patha19_foundation_checks in
   order; ltRepresentationResiduals preserves the 3 LT representation gates. *)
foundationResiduals[d_] := Module[
  {
    rho4, rho3, psi, velocity, numberRate, energy, action, lpsiDensity,
    electricField, magneticField, maxwellCoeff
  },
  rho4 = d["rho0"];
  rho3 = d["rho3"];
  psi = d["psi"];
  velocity = d["velocity"];
  numberRate = -d["T"];
  energy = d["energy"];
  action = d["action"];
  lpsiDensity = d["lagrangianDensity"];
  electricField = d["electricField"];
  magneticField = d["magneticField"];
  maxwellCoeff = d["maxwellCoeff"];
  <|
    "pathA_19_F2: 4D-bulk closed-3-surface number flux J_bulk" ->
      dimResidual[rho4 + velocity + 3 d["L"], numberRate],
    "pathA_19_F2: 3D-brane reduced 2-sphere number flux J_brane" ->
      dimResidual[rho3 + velocity + 2 d["L"], numberRate],
    "pathA_19_F2: 4D-bulk volumetric flux Q_vol=rho^-1 J" ->
      dimResidual[-rho4 + numberRate, 4 d["L"] - d["T"]],
    "pathA_19_F2: 3D-brane volumetric flux Q_vol=rho_3^-1 J" ->
      dimResidual[-rho3 + numberRate, 3 d["L"] - d["T"]],
    "pathA_19_F2: constituent mass flux m_GNLS*J" ->
      dimResidual[d["mGNLS"] + numberRate, d["M"] - d["T"]],
    "pathA_19_F1: conditional defect rest-frequency conversion hbar*J/c_gamma^2" ->
      dimResidual[action + numberRate - 2 velocity, d["M"]],
    "pathA_19_F2: bulk continuity equation" ->
      homResidual[<|"partial_t rho" -> rho4 - d["T"], "div_4(rho v)" -> rho4 + velocity - d["L"]|>],
    "pathA_19_F3: sound-speed law 5*K*rho^4/m" ->
      dimResidual[d["K"] + 4 rho4 - d["mGNLS"], 2 velocity],
    "pathA_19_F3: EOS enthalpy scale h0=(m_GNLS*c_s0^2)/4" ->
      dimResidual[d["mGNLS"] + 2 velocity, energy],
    "pathA_19_F3: GNLS healing length sqrt(hbar^2/(2*m_GNLS*h0))" ->
      dimResidual[(2 action - d["mGNLS"] - energy)/2, d["L"]],
    "pathA_19_F3: parent GNLS Lagrangian density terms" ->
      homResidual[<|
        "i*hbar*psi*partial_t psi" -> action - d["T"] + 2 psi,
        "hbar^2/(2m)*|D_i psi|^2" -> 2 action - d["mGNLS"] + 2 (psi - d["L"]),
        "V_conf*rho" -> energy + rho4,
        "U=K*rho^5/4" -> d["K"] + 5 rho4
      |>],
    "pathA_19_F3: spatial gauge minimal-coupling dimension q*A_i/hbar" ->
      dimResidual[d["qAi"] - action, -d["L"]],
    "pathA_19_F3: localized Maxwell sector with explicit c factors" ->
      homResidual[<|
        "(Z/mu0)*E_i^2" -> maxwellCoeff + 2 electricField,
        "(Z/mu0)*c^2*B_ij^2" -> maxwellCoeff + 2 velocity + 2 magneticField,
        "A0*J0_ext" -> d["qA0"] + rho4,
        "Ai*Ji_ext" -> d["qAi"] + rho4 + velocity
      |>],
    "pathA_19_F3: wall action density before dt*dw integration" ->
      homResidual[<|
        "mu_eta*(partial_t eta)^2" -> d["muWall"] + 2 velocity,
        "T_w*(partial_w eta)^2" -> d["Tw"],
        "K_eta*eta^2" -> d["USigmaRR"] + 2 d["L"]
      |>]
  |>
];

ltRepresentationResiduals[] := Module[
  {
    ell, time, zero, rho, psi, velocity, actionLT, energyLT, massLT,
    forceLT, kLT, densityLT, eFieldLT, bFieldLT, maxwellCoeffLT
  },
  ell = {1, 0, 0};
  time = {0, 1, 0};
  zero = {0, 0, 0};
  rho = -4 ell;
  psi = -2 ell;
  velocity = ell - time;
  actionLT = {2, -1, 0};
  energyLT = {2, -2, 0};
  massLT = zero;
  forceLT = {1, -2, 0};
  kLT = energyLT - 4 rho;
  densityLT = energyLT - 4 ell;
  eFieldLT = forceLT;
  bFieldLT = actionLT - 2 ell;
  maxwellCoeffLT = densityLT - 2 eFieldLT;
  <|
    "pathA_19_LT_representation: local GNLS terms after projecting m_GNLS to dimensionless" ->
      homResidual[<|
        "i*hbar*psi*partial_t psi" -> actionLT - time + 2 psi,
        "hbar^2/(2m)*|D_i psi|^2" -> 2 actionLT - massLT + 2 (psi - ell),
        "V_conf*rho" -> energyLT + rho,
        "U=K*rho^5/4" -> kLT + 5 rho
      |>],
    "pathA_19_LT_representation: local Maxwell terms after M projection" ->
      homResidual[<|
        "(Z/mu0)*E_i^2" -> maxwellCoeffLT + 2 eFieldLT,
        "(Z/mu0)*c^2*B_ij^2" -> maxwellCoeffLT + 2 velocity + 2 bFieldLT
      |>],
    "pathA_19_LT_representation: local wall terms after M projection" ->
      homResidual[<|
        "mu_eta*(partial_t eta)^2" -> {-1, 0, 0} + 2 velocity,
        "T_w*(partial_w eta)^2" -> forceLT,
        "K_eta*eta^2" -> {-1, -2, 0} + 2 ell
      |>]
  |>
];

flaggedResiduals[d_] := Module[{velocity, formal, observed, ltGate},
  velocity = d["velocity"];
  formal = d["G4"] + 5 velocity - 5 d["L"] - 5 velocity;
  observed = d["G3"] + 5 velocity - 5 d["L"] - 5 velocity;
  ltGate = {3, -2, 0} + 5 velocity - 5 d["L"] - 5 velocity;
  {
    <|
      "token" -> "formal_4D_R_norm_target_not_dimensionless_without_conversion",
      "actual" -> formal,
      "factor" -> factorToReach[d["zero"], formal]
    |>,
    <|
      "token" -> "observed_3D_GR_target_not_dimensionless_without_conversion",
      "actual" -> observed,
      "factor" -> factorToReach[d["zero"], observed]
    |>,
    <|
      "token" -> "LT_R_norm_gate_fails_without_new_conversion_factor",
      "actual" -> ltGate,
      "factor" -> factorToReach[d["zero"], ltGate]
    |>
  }
];

allHarnessPassQ[d_] := AllTrue[
  Join[Values[foundationResiduals[d]], Values[ltRepresentationResiduals[]]],
  TrueQ[cleanZero[#] === 0] &
];

ltRejectionGate[d_, repairResiduals_: False] := Module[{residualGateFails, localPasses},
  residualGateFails = If[
    TrueQ[repairResiduals],
    False,
    AnyTrue[flaggedResiduals[d], ! TrueQ[cleanZero[dimResidual[dropM[#["actual"]], d["zero"]]] === 0] &]
  ];
  localPasses = AllTrue[Values[ltRepresentationResiduals[]], TrueQ[cleanZero[#] === 0] &];
  TrueQ[allHarnessPassQ[d] && localPasses && residualGateFails]
];

verdictFromPredicate[d_, mDefectDerivedHere_, repairLTResiduals_: False] := If[
  TrueQ[ltRejectionGate[d, repairLTResiduals] && ! TrueQ[mDefectDerivedHere]],
  "RETAIN_L_T_M",
  "NOT_RETAIN_L_T_M"
];

runHarnessChecks[d_] := (
  subheading["Harness-mapped 17 dimensional checks"];
  KeyValueMap[Function[{name, residual}, expectZero[name, residual]], foundationResiduals[d]];
  KeyValueMap[Function[{name, residual}, expectZero[name, residual]], ltRepresentationResiduals[]];
  expectZero[
    "harness check count is 14 foundation + 3 LT representation",
    Length[foundationResiduals[d]] + Length[ltRepresentationResiduals[]] - 17
  ];
  Print["  hbar*J/c_gamma^2 = M check above is dimensional-only, not a mass derivation."]
);

runFlaggedResiduals[d_] := Module[{flags},
  subheading["Carried flagged residual gates"];
  flags = flaggedResiduals[d];
  Scan[
    Function[flag,
      Print["  ", flag["token"], ": actual ", dimString[flag["actual"]], "; factor needed ", dimString[flag["factor"]]];
      expectNonzero[
        flag["token"] <> " remains non-dimensionless after M drop",
        dimResidual[dropM[flag["actual"]], d["zero"]]
      ];
      expectDim[
        flag["token"] <> " exact conversion factor",
        flag["factor"] + flag["actual"],
        d["zero"]
      ]
    ],
    flags
  ]
];

runVerdict[d_] := Module[{mDefectDerivedHere, verdict},
  subheading["Mass fork and computed verdict"];
  mDefectDerivedHere = False;
  Print["  m_defect_derived_here=", mDefectDerivedHere];
  Print["  INFLOW_MASS_SOURCE_MISSING: m_defect is NOT emergent at this foundation gate."];
  Print["  hbar*J/c_gamma^2 = M is a dimensional conversion only, not a mass theorem."];
  verdict = verdictFromPredicate[d, mDefectDerivedHere];
  expectZero["computed verdict is RETAIN_L_T_M", If[verdict === "RETAIN_L_T_M", 0, 1]];
  expectBool[
    "verdict flips if LT residual gates are counterfactually repaired",
    verdictFromPredicate[d, False, True] =!= "RETAIN_L_T_M"
  ];
  expectBool[
    "verdict flips if m_defect_derived_here=True",
    verdictFromPredicate[d, True] =!= "RETAIN_L_T_M"
  ];
  verdict
];

runFirewall[d_] := Module[
  {psi, action, energy, rho4, velocity, wrongK, pinDims, baselineNullity, dropMat, corruptRelation},
  subheading["Able-to-fail dimensional firewall"];
  psi = d["psi"];
  action = d["action"];
  energy = d["energy"];
  rho4 = d["rho0"];
  velocity = d["velocity"];
  expectFail[
    "drop hbar from kinetic term breaks GNLS density homogeneity",
    homResidual[<|
      "kinetic_without_hbar" -> -d["T"] + 2 psi,
      "gradient" -> 2 action - d["mGNLS"] + 2 (psi - d["L"]),
      "V_conf*rho" -> energy + rho4
    |>]
  ];
  wrongK = {17, -2, 1};
  expectFail[
    "wrong K exponent M L^17 T^-2 breaks EOS density",
    homResidual[<|
      "kinetic" -> action - d["T"] + 2 psi,
      "U_wrong=K*rho^5/4" -> wrongK + 5 rho4
    |>]
  ];
  expectFail[
    "drop rho0^4 from c_s0^2=5*K*rho0^4/m_GNLS breaks velocity dimension",
    dimResidual[(d["K"] - d["mGNLS"])/2, velocity]
  ];
  pinDims = {d["L"], d["cS0"], d["hbar"], d["mGNLS"]};
  baselineNullity = Length[pinDims] - MatrixRank[Transpose[pinDims]];
  dropMat = Transpose[{d["L"], d["hbar"], d["mGNLS"]}];
  expectBool[
    "pin corruption by dropping c_s0 pin changes nullity",
    3 - MatrixRank[dropMat] =!= baselineNullity
  ];
  corruptRelation = d["L"] + d["zero"] - d["hbar"] + d["mGNLS"];
  expectFail[
    "pin corruption c_s0 dimensionless breaks a=hbar/(m_GNLS*c_s0)",
    dimResidual[corruptRelation, d["zero"]]
  ]
];

printGapBlock[] := (
  subheading["Carried gaps printed as provenance"];
  Print["  A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT: a is branch geometry, not a base invariant."];
  Print["  EOS_FROM_GNLS_FACTOR: h0=(m_GNLS*c_s0^2)/4 and xi_h=sqrt(2)*hbar/(m_GNLS*c_s0) carried to pathA_20."];
  Print["  NO_NET_ACCRETION_BC_UNDERIVED: no-net-accretion is a boundary condition, not derived here."];
  Print["  M_TO_G_UNIFICATION: defect-mass/back-reaction relation is deferred to pathA_21."];
  Print["  SCALE_MAP_INPUTS: pathA_22 consumes J, a(branch), rho0, K, m_GNLS, hbar, 3D-reduction factors."]
);

printVerdictLabels[verdict_] := (
  Print[""];
  Print["Verdict labels:"];
  Print[
    "  ledger earned-label (NOT a source verdict token): ",
    "DIMENSIONAL_FOUNDATION_LTM_RETAINED  ",
    "(base {L,T,M}; a=hbar/(m*c_s0); xi_h=sqrt(2)*hbar/(m*c_s0))"
  ];
  Print["  source top-line verdict: ", verdict, "   (PASS_WITH_NAMED_RESIDUALS)"];
  Print[
    "  labeled non-derivation (carried gap): ",
    "m_defect NOT emergent -- INFLOW_MASS_SOURCE_MISSING (m_defect_derived_here=False)"
  ];
  Print[
    "  carried residuals: LT_R_norm_gate_fails_without_new_conversion_factor ",
    "(REJECTS_TRUE_LT_BASE); ",
    "formal_4D_R_norm_target_not_dimensionless_without_conversion; ",
    "observed_3D_GR_target_not_dimensionless_without_conversion"
  ];
  Print[
    "  carried forward (provenance, §3c): ",
    "A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT; EOS_FROM_GNLS_FACTOR; ",
    "NO_NET_ACCRETION_BC_UNDERIVED; M_TO_G_UNIFICATION; SCALE_MAP_INPUTS"
  ]
);

Module[{dims, verdict, ok},
  heading["ledger_stage004_gnls_action_dimensional_foundation Mathematica audit"];
  ok = Catch[
    dims = deriveDictionary[];
    runTwoTierDictionary[dims];
    runSymbolicCore[];
    runPinAnalysis[dims];
    runHarnessChecks[dims];
    runFlaggedResiduals[dims];
    verdict = runVerdict[dims];
    runFirewall[dims];
    printGapBlock[];
    printVerdictLabels[verdict];
    True,
    "ledgerStage004Failure",
    Function[{name, tag}, False]
  ];

  Print[""];
  Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
  If[TrueQ[ok],
    Print["OVERALL PASS: Mathematica derived ledger_stage004 GNLS dimensional foundation exactly"];
    Exit[0],
    Print["OVERALL FAIL: Mathematica stage004 audit did not close"];
    Exit[1]
  ]
]
