(* Ledger stage007 Mathematica audit: pathA_35 G0 shear-surface freeze.

   Standalone, print-only, no arguments, no exports.  This route uses a native
   byte scanner over BinaryReadList data, exponent associations in {L,T,M},
   native Association set partitions, projector/NullSpace rank constructions,
   and Count/GroupBy tallies for the historical and operative enumerations.
*)

ClearAll[
  raise, heading, subheading, cleanZero, assertExact, expectZero, expectBool,
  expectNonzero, expectFail, dim, dimVector, dimFromVector, dimAdd, dimMul,
  dimSub, dimDiv, dimPow, dimResidual, expectDim, dimString, homResidual,
  readReportBytes, extractFenceBytes, hashBytes, hashResidual, byteContainsQ,
  missingFenceResidual, runFreezeFidelity, historicalActionSummands,
  actionPartitionResidual, runPostD16ActionPartition, runDimensionalFirewall,
  checkSameDims, computeFlatBraneDof, runFlatBraneDof, computePostD16Dof,
  runPostD16Dof, structuralPostulates, driftTable, categoryCounts,
  driftExpectedN, setResidual, historicalDriftResidual,
  postD16DriftResidual, driftDimResidual,
  runComputedDriftLedger, distinctnessFailureResidual, runMuRFirewall,
  printStructuralPostulates, printProvenanceBlocks, printVerdictLabels,
  passCount, failCount
];

passCount = 0;
failCount = 0;

raise[msg_] := Throw[msg, "ledgerStage007Failure"];

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

dim[l_, t_, m_] := <|"L" -> l, "T" -> t, "M" -> m|>;
dimVector[d_Association] := Lookup[d, {"L", "T", "M"}];
dimFromVector[v_] := AssociationThread[{"L", "T", "M"}, v];
dimAdd[dims__] := dimFromVector[Total[dimVector /@ {dims}]];
dimMul[n_, d_Association] := dimFromVector[n dimVector[d]];
dimSub[left_Association, right_Association] := dimAdd[left, dimMul[-1, right]];
dimDiv[left_Association, right_Association] := dimSub[left, right];
dimPow[d_Association, n_] := dimMul[n, d];
dimResidual[actual_Association, expected_Association] := FullSimplify[(dimVector[actual] - dimVector[expected]).(dimVector[actual] - dimVector[expected])];
expectDim[name_, actual_Association, expected_Association] := expectZero[name, dimResidual[actual, expected]];

dimString[d_Association] := Module[{labels, powers, pieces},
  labels = {"L", "T", "M"};
  powers = dimVector[d];
  pieces = DeleteCases[
    Table[
      Which[
        TrueQ[powers[[i]] === 0], Nothing,
        TrueQ[powers[[i]] === 1], labels[[i]],
        True, labels[[i]] <> "^" <> ToString[InputForm[powers[[i]]]]
      ],
      {i, 3}
    ],
    Nothing
  ];
  If[pieces === {}, "1", StringRiffle[pieces, " "]]
];

homResidual[terms_Association] := Module[{vals, ref},
  vals = Values[terms];
  If[Length[vals] == 0, raise["homogeneity check requires at least one term"]];
  ref = First[vals];
  FullSimplify[Total[dimResidual[#, ref] & /@ Rest[vals]]]
];

scriptPath = ExpandFileName[
  If[StringQ[$InputFileName] && $InputFileName =!= "",
    $InputFileName,
    FileNameJoin[{"research", "pde_ledger_v2", "mathematica", "ledger_stage007_shear_surface_g0_freeze_mathematica_audit.wl"}]
  ]
];
repoRoot = Nest[ParentDirectory, DirectoryName[scriptPath], 3];
reportRoot = FileNameJoin[{repoRoot, "software", "stage1_solver", "reports"}];
t0Report = FileNameJoin[{reportRoot, "pathA_24_T0_freeze.md"}];
g0Report = FileNameJoin[{reportRoot, "pathA_35_G0_freeze.md"}];

expectedT0Hash = "8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064";
expectedG0Hash = "d9520d3819c3f718290f9d0be57138c07d5bf02d2237106478e17b6a1e389ac3";
expectedG0Short = "d9520d3819c3";
expectedDriftToken = "SECOND_MEDIUM_DRIFT_AT_FREEZE(11)";
expectedFreezeToken = "T0_SHEAR_FROZEN(" <> expectedG0Short <> ")";
expectedPostD16ActionToken = "POST_D16_ACTION{S_GNLS,gL_Mac,gL_uw}_OF(" <> expectedG0Short <> ")";
expectedPostD16DriftToken = "POST_D16_DRIFT(7)";
historicalActionNames = {"S_GNLS", "L_pol", "gL_Mac", "gL_Pu", "gL_uw"};
retiredActionNames = {"L_pol", "gL_Pu"};
postD16ActionNames = {"S_GNLS", "gL_Mac", "gL_uw"};
driftDimTargets = <|
  "rho_br" -> "rho_br",
  "mu_R" -> "mu_R",
  "lambda_Pu" -> "lambda_Pu",
  "Omega_w" -> "Omega_w",
  "g_ell(w)" -> "g_ell"
|>;

readReportBytes[path_] := Module[{},
  If[! FileExistsQ[path], raise["missing report: " <> path]];
  BinaryReadList[path, "Byte"]
];

extractFenceBytes[path_, label_] := Module[
  {bytes, opening, closing, opens, contentStart, tail, closes, closingStart, block},
  bytes = readReportBytes[path];
  opening = ToCharacterCode["```" <> label <> "\n", "UTF-8"];
  closing = ToCharacterCode["```\n", "UTF-8"];
  opens = SequencePosition[bytes, opening];
  If[Length[opens] == 0, raise["missing `" <> label <> "` fence in " <> path]];
  If[Length[opens] > 1, raise["ambiguous `" <> label <> "` fence in " <> path <> ": found " <> ToString[Length[opens]] <> " blocks"]];
  contentStart = opens[[1, 2]] + 1;
  tail = Drop[bytes, opens[[1, 2]]];
  closes = SequencePosition[tail, closing];
  If[Length[closes] == 0, raise["unterminated `" <> label <> "` fence in " <> path]];
  closingStart = opens[[1, 2]] + closes[[1, 1]];
  block = If[closingStart > contentStart, Take[bytes, {contentStart, closingStart - 1}], {}];
  <|
    "Path" -> path,
    "Label" -> label,
    "Bytes" -> block,
    "Start0" -> contentStart - 1,
    "End0" -> closingStart - 1,
    "Start1" -> contentStart,
    "End1" -> closingStart - 1,
    "Length" -> Length[block]
  |>
];

hashBytes[bytes_] := Hash[ByteArray[bytes], "SHA256", "HexString"];
hashResidual[bytes_, expected_] := If[hashBytes[bytes] === expected, 0, 1];
byteContainsQ[haystack_, needle_] := SequencePosition[haystack, needle] =!= {};

missingFenceResidual[] := Module[{result},
  result = Catch[
    extractFenceBytes[g0Report, "nonexistent-freeze-action"];
    "NO_FAILURE",
    "ledgerStage007Failure",
    Function[{msg, tag}, "EXPECTED_FAILURE"]
  ];
  If[result === "EXPECTED_FAILURE", 1, 0]
];

runFreezeFidelity[] := Module[{t0, g0, corrupted},
  subheading["HISTORICAL freeze (pre-Decision-16, freeze-as-run): fidelity byte audit"];
  t0 = extractFenceBytes[t0Report, "freeze-action"];
  g0 = extractFenceBytes[g0Report, "freeze-action"];
  Print["  T0 report: ", t0["Path"]];
  Print["  G0 report: ", g0["Path"]];
  Print[
    "  G0 freeze-action byte range (informative only): 0-based [",
    g0["Start0"], ",", g0["End0"], "), 1-based ", g0["Start1"], "-",
    g0["End1"], ", length ", g0["Length"]
  ];
  Print[
    "  T0 freeze-action byte range (informative only): 0-based [",
    t0["Start0"], ",", t0["End0"], "), 1-based ", t0["Start1"], "-",
    t0["End1"], ", length ", t0["Length"]
  ];
  expectZero["T0 freeze-action SHA-256 matches frozen report", hashResidual[t0["Bytes"], expectedT0Hash]];
  expectBool["G0 short hash is prefix of frozen full SHA-256 constant", StringStartsQ[expectedG0Hash, expectedG0Short]];
  expectZero["G0 freeze-action SHA-256 matches frozen report", hashResidual[g0["Bytes"], expectedG0Hash]];
  expectBool["byte-identical T0 freeze-action block is embedded inside G0 block", byteContainsQ[g0["Bytes"], t0["Bytes"]]];
  corrupted = ReplacePart[g0["Bytes"], 1 -> BitXor[First[g0["Bytes"]], 1]];
  expectFail["hash tooth: one-byte in-memory G0 corruption trips SHA-256 mismatch", hashResidual[corrupted, expectedG0Hash]];
  expectFail["hash tooth: nonexistent fence tag trips extractor missing-fence path", missingFenceResidual[]];
  {t0, g0}
];

historicalActionSummands[] := <|
  "S_GNLS" -> HoldComplete[sGnlsExisting],
  "L_pol" -> HoldComplete[
    (1/2) mSym rhoSym aSym^2 dtPSq -
    (1/2) mSym rhoSym csSq aSym^2 gradPSq -
    (1/4) mSym rhoSym csSq radialPSq
  ],
  "gL_Mac" -> HoldComplete[gEll ((1/2) rhoBrSym dtUSq - (1/2) muRSym omegaUSq)],
  "gL_Pu" -> HoldComplete[-gEll lambdaPuSym varpiDotOmegaU],
  "gL_uw" -> HoldComplete[gEll ((1/2) rhoBrSym dtUwSq - (1/2) rhoBrSym omegaWSym^2 uwSq)]
|>;

actionPartitionResidual[historical_Association, operative_Association, retired_Association] := Module[
  {historicalKeys, operativeKeys, retiredKeys, residual},
  historicalKeys = Keys[historical];
  operativeKeys = Keys[operative];
  retiredKeys = Keys[retired];
  residual = 0;
  residual += If[historicalKeys === historicalActionNames, 0, 1];
  residual += If[Sort[retiredKeys] === Sort[retiredActionNames], 0, 1];
  residual += If[Intersection[operativeKeys, retiredKeys] === {}, 0, 1];
  residual += If[Sort[Union[operativeKeys, retiredKeys]] === Sort[historicalKeys], 0, 1];
  residual += If[Sort[operativeKeys] === Sort[Complement[historicalKeys, retiredActionNames]], 0, 1];
  residual += If[Intersection[operativeKeys, retiredActionNames] === {}, 0, 1];
  Do[residual += If[SameQ[operative[key], historical[key]], 0, 1], {key, Intersection[operativeKeys, historicalKeys]}];
  Do[residual += If[SameQ[retired[key], historical[key]], 0, 1], {key, Intersection[retiredKeys, historicalKeys]}];
  residual
];

runPostD16ActionPartition[] := Module[
  {historical, retired, operative, token, altered, movedRetired, droppedSurvivor},
  subheading["OPERATIVE post-Decision-16 action: native symbolic summand set partition"];
  historical = historicalActionSummands[];
  retired = KeyTake[historical, retiredActionNames];
  operative = KeyDrop[historical, retiredActionNames];
  expectZero["historical symbolic action summands match the five-summand frozen S_G0 grammar", If[Keys[historical] === historicalActionNames, 0, 1]];
  expectZero["post-D16 operative and retired action sets form an exact disjoint partition", actionPartitionResidual[historical, operative, retired]];
  expectZero["post-D16 operative action names are exactly the ordered survivor set", If[Keys[operative] === postD16ActionNames, 0, 1]];
  token = "POST_D16_ACTION{" <> StringRiffle[Keys[operative], ","] <> "}_OF(" <> expectedG0Short <> ")";
  expectZero["post-D16 action token is assembled from computed survivor names and historical hash", If[token === expectedPostD16ActionToken, 0, 1]];
  Print["  historical summands = {", StringRiffle[Keys[historical], ","], "}"];
  Print["  retired complement = {", StringRiffle[Keys[retired], ","], "}"];
  Print["  operative token = ", token];
  Print["  Route A is a symbolic action-summand partition over the immutable hash anchor; no byte-substring surgery is performed."];

  altered = Join[operative, <|"gL_Mac" -> HoldComplete[gEll lMacHistorical + deltaSurvivorMutation]|>];
  expectFail["action-partition tooth: alter survivor gL_Mac trips definition fidelity", actionPartitionResidual[historical, altered, retired]];
  movedRetired = Join[operative, KeyTake[retired, {"L_pol"}]];
  expectFail["action-partition tooth: move retired L_pol into operative set trips disjointness and absence", actionPartitionResidual[historical, movedRetired, retired]];
  droppedSurvivor = KeyDrop[operative, {"gL_uw"}];
  expectFail["action-partition tooth: drop survivor gL_uw trips union-equals-historical", actionPartitionResidual[historical, droppedSurvivor, retired]];
  expectZero["baseline action partition remains valid after copy-mutation teeth", actionPartitionResidual[historical, operative, retired]];
  operative
];

checkSameDims[name_, terms_Association, expected_] := (
  KeyValueMap[expectDim[name <> " part " <> #1, #2, expected] &, terms];
  expectZero[name <> " homogeneous", homResidual[terms]]
);

runDimensionalFirewall[] := Module[
  {
    Z, L, T, M, bulkLag, braneLag, actionDim, stress, eomUOp,
    dm, dhbar, drho, dK, da, du, dP, dw, dellG, dg, dgrad, ddt,
    ddtMeasure, dd4x, dv, dk, domega, drhoBr, dmuR, dlambdaPu,
    dOmegaW, dcs2, dDtP, dgradP, w, ellG, gExpr, normalization,
    dDtu, dcurlu, duw, dDtuw, dvarpi, macKin, macCurl, puTerm,
    uwKin, uwGap, twa, slope, stressSlope, opRho, opMu, cGammaSq
  },
  subheading["Dimensional firewalls: OPERATIVE live surface + RETIRED-HISTORICAL P surface"];
  Z = dim[0, 0, 0]; L = dim[1, 0, 0]; T = dim[0, 1, 0]; M = dim[0, 0, 1];
  bulkLag = dimDiv[M, dimAdd[dimPow[L, 2], dimPow[T, 2]]];
  braneLag = dimDiv[M, dimAdd[L, dimPow[T, 2]]];
  actionDim = dimDiv[dimAdd[M, dimPow[L, 2]], T];
  stress = bulkLag;
  eomUOp = dimDiv[M, dimAdd[dimPow[L, 3], dimPow[T, 2]]];

  dm = M;
  dhbar = actionDim;
  drho = dimPow[L, -4];
  dK = dimAdd[M, dimPow[L, 18], dimPow[T, -2]];
  da = L; du = L; dP = Z; dw = L; dellG = L; dg = dimPow[L, -1];
  dgrad = dimPow[L, -1]; ddt = dimPow[T, -1]; ddtMeasure = T; dd4x = dimPow[L, 4];
  dv = dimDiv[L, T]; dk = dimPow[L, -1]; domega = dimPow[T, -1];
  drhoBr = dimAdd[M, dimPow[L, -3]];
  dmuR = braneLag; dlambdaPu = braneLag; dOmegaW = domega;

  Print["  OPERATIVE LIVE firewall: kept GNLS parent, L_Mac, L_uw, g_ell, projected traction, O_u, c_gamma^2, and omega_uw,bare^2."];
  dcs2 = dimDiv[dimAdd[dK, dimMul[4, drho]], dm];
  expectDim["OPERATIVE kept GNLS c_s^2(rho)=5 K rho^4/m", dcs2, dimDiv[dimPow[L, 2], dimPow[T, 2]]];
  expectDim["OPERATIVE kept GNLS U(rho)=(K/4)rho^5", dimAdd[dK, dimMul[5, drho]], bulkLag];
  expectDim[
    "OPERATIVE kept GNLS quantum pressure (hbar^2/(8 m rho))(partial_i rho)^2",
    dimAdd[dimMul[2, dhbar], dimMul[-1, dm], dimMul[-1, drho], dimMul[2, dimAdd[dgrad, drho]]],
    bulkLag
  ];
  expectDim["OPERATIVE kept GNLS bulk kinetic stress scale m rho v_i v_j", dimAdd[dm, drho, dimMul[2, dv]], stress];

  Print["  RETIRED-HISTORICAL firewall: L_pol and inherited P-sector coefficients remain dimensionally audited as freeze-as-run records."];
  dDtP = ddt;
  dgradP = dimAdd[dgrad, dP];
  expectDim["RETIRED-HISTORICAL T0 P^i dimensionless", dP, Z];
  expectDim["RETIRED-HISTORICAL T0 L_pol kinetic term", dimAdd[dm, drho, dimMul[2, da], dimMul[2, dDtP]], bulkLag];
  expectDim["RETIRED-HISTORICAL T0 L_pol Frank term", dimAdd[dm, drho, dcs2, dimMul[2, da], dimMul[2, dgradP]], bulkLag];
  expectDim["RETIRED-HISTORICAL T0 L_pol radial term", dimAdd[dm, drho, dcs2], bulkLag];
  expectZero[
    "RETIRED-HISTORICAL T0 L_pol all terms homogeneous",
    homResidual[<|
      "P_kinetic" -> dimAdd[dm, drho, dimMul[2, da], dimMul[2, dDtP]],
      "P_Frank" -> dimAdd[dm, drho, dcs2, dimMul[2, da], dimMul[2, dgradP]],
      "P_radial" -> dimAdd[dm, drho, dcs2]
    |>]
  ];
  expectDim["RETIRED-HISTORICAL T0 couple-stress inertia m rho a^2", dimAdd[dm, drho, dimMul[2, da]], dimAdd[M, dimPow[L, -2]]];
  expectDim["RETIRED-HISTORICAL T0 couple-stress stiffness m rho c_s^2 a^2", dimAdd[dm, drho, dcs2, dimMul[2, da]], dimDiv[M, dimPow[T, 2]]];
  expectDim["RETIRED-HISTORICAL T0 bulk radial scale m rho c_s^2", dimAdd[dm, drho, dcs2], bulkLag];

  expectDim["profile ratio w/ell_g is dimensionless", dimDiv[dw, dellG], Z];
  expectDim["profile g_ell(w)=exp(-(w/ell_g)^2)/(sqrt(pi) ell_g)", dg, dimPow[L, -1]];
  expectDim["profile measure dw g_ell(w) is dimensionless", dimAdd[L, dg], Z];
  Clear[w, ellG];
  gExpr = Exp[-(w/ellG)^2]/(Sqrt[Pi] ellG);
  normalization = FullSimplify[Integrate[gExpr, {w, -Infinity, Infinity}, Assumptions -> ellG > 0], ellG > 0];
  expectZero["derived Gaussian normalization integral int g_ell(w) dw = 1", normalization - 1];

  dDtu = dimAdd[ddt, du];
  dcurlu = dimAdd[dgrad, du];
  duw = du;
  dDtuw = dimAdd[ddt, duw];
  dvarpi = dP;
  expectDim["target [u^a]=L", du, L];
  expectDim["target [u_w]=L", duw, L];
  expectDim["target [rho_br]=M L^-3", drhoBr, dimAdd[M, dimPow[L, -3]]];
  expectDim["target [mu_R]=M L^-1 T^-2", dmuR, braneLag];
  expectDim["RETIRED-HISTORICAL target [lambda_Pu]=M L^-1 T^-2", dlambdaPu, braneLag];
  expectDim["target [Omega_w]=T^-1", dOmegaW, dimPow[T, -1]];
  expectDim["target [g_ell]=L^-1", dg, dimPow[L, -1]];
  expectDim["target [ell_g]=L", dellG, L];

  macKin = dimAdd[drhoBr, dimMul[2, dDtu]];
  macCurl = dimAdd[dmuR, dimMul[2, dcurlu]];
  expectDim["L_Mac kinetic rho_br (partial_t u^a)^2", macKin, braneLag];
  expectDim["L_Mac MacCullagh curl mu_R Omega_u^a Omega_u^a", macCurl, braneLag];
  expectZero["L_Mac homogeneous", homResidual[<|"kinetic" -> macKin, "curl" -> macCurl|>]];
  puTerm = dimAdd[dlambdaPu, dvarpi, dcurlu];
  expectDim["RETIRED-HISTORICAL L_Pu parity-repaired lambda_Pu varpi_a Omega_u^a", puTerm, braneLag];
  uwKin = dimAdd[drhoBr, dimMul[2, dDtuw]];
  uwGap = dimAdd[drhoBr, dimMul[2, dOmegaW], dimMul[2, duw]];
  expectDim["L_uw kinetic rho_br (partial_t u_w)^2", uwKin, braneLag];
  expectDim["L_uw bare gap rho_br Omega_w^2 u_w^2", uwGap, braneLag];
  expectZero["L_uw homogeneous", homResidual[<|"kinetic" -> uwKin, "gap" -> uwGap|>]];

  KeyValueMap[
    expectDim["OPERATIVE brane bulk representation " <> #1, dimAdd[dg, #2], bulkLag] &,
    <|
      "g_ell L_Mac kinetic" -> macKin,
      "g_ell L_Mac curl" -> macCurl,
      "g_ell L_uw kinetic" -> uwKin,
      "g_ell L_uw gap" -> uwGap
    |>
  ];
  expectDim["RETIRED-HISTORICAL brane bulk representation g_ell L_Pu", dimAdd[dg, puTerm], bulkLag];

  expectDim["action measure int dt d^4X bulk_lag", dimAdd[ddtMeasure, dd4x, bulkLag], actionDim];
  expectDim["action measure int dt d^4X g_ell(w) L_brane", dimAdd[ddtMeasure, dd4x, dg, braneLag], actionDim];

  twa = dimAdd[dm, drho, dv, dv];
  slope = dimAdd[dgrad, duw];
  stressSlope = dimAdd[stress, slope];
  expectDim["projected traction T_wa=m rho v_w v_a", twa, stress];
  expectDim["projected traction partial_b u_w is dimensionless", slope, Z];
  expectDim["projected traction slope mixing", stressSlope, stress];
  expectZero["full projected traction T_na homogeneous", homResidual[<|"T_wa" -> twa, "slope" -> stressSlope|>]];

  expectFail["dim ablation drop_m_from_T_wa", dimResidual[dimAdd[drho, dv, dv], stress]];
  expectFail["dim ablation MacCullagh_without_curl", dimResidual[dimAdd[dmuR, dimMul[2, du]], braneLag]];

  opRho = dimAdd[drhoBr, dimMul[2, domega]];
  opMu = dimAdd[dmuR, dimMul[2, dk]];
  expectDim["OPERATIVE linearization O_u rho_br omega^2 term", opRho, eomUOp];
  expectDim["OPERATIVE linearization O_u mu_R k^2 term", opMu, eomUOp];
  expectZero["OPERATIVE linearization O_u homogeneous", homResidual[<|"rho_br omega^2" -> opRho, "mu_R k^2" -> opMu|>]];
  cGammaSq = dimDiv[dmuR, drhoBr];
  expectDim["OPERATIVE target [c_gamma^2=mu_R/rho_br]=L^2 T^-2", cGammaSq, dimDiv[dimPow[L, 2], dimPow[T, 2]]];
  expectDim["OPERATIVE linearization omega_T^2=c_gamma^2 k^2", dimAdd[cGammaSq, dimMul[2, dk]], dimPow[T, -2]];
  expectDim["OPERATIVE linearization omega_uw,bare^2=Omega_w^2", dimMul[2, dOmegaW], dimPow[T, -2]];
  Print["  RETIRED-HISTORICAL linearization block: omega_P^2 and omega_radial^2 are audited only as removed P-sector modes, never as live survivors."];
  expectDim["RETIRED-HISTORICAL linearization omega_P^2=c_s^2 k^2", dimAdd[dcs2, dimMul[2, dk]], dimPow[T, -2]];
  expectDim["RETIRED-HISTORICAL linearization P radial gap 2 c_s^2/a^2", dimSub[dcs2, dimMul[2, da]], dimPow[T, -2]];
  expectZero[
    "RETIRED-HISTORICAL linearization omega_radial^2 homogeneous",
    homResidual[<|"c_s^2 k^2" -> dimAdd[dcs2, dimMul[2, dk]], "2 c_s^2/a^2" -> dimSub[dcs2, dimMul[2, da]]|>]
  ];
  expectDim["RETIRED-HISTORICAL linearization Fourier P-u monomial lambda_Pu P k u", dimAdd[dlambdaPu, dP, dk, du], braneLag];

  Print["  retired-historical parity record: direct P_parallel.Omega_u was parity-ODD and excluded."];
  Print["  retired-historical parity record: w_hat.(P_parallel x Omega_u) was parity-EVEN and frozen before Decision 16 retired the P-u block."];

  <|
    "P" -> dP,
    "u" -> du,
    "u_w" -> duw,
    "rho_br" -> drhoBr,
    "mu_R" -> dmuR,
    "lambda_Pu" -> dlambdaPu,
    "Omega_w" -> dOmegaW,
    "g_ell" -> dg,
    "ell_g" -> dellG,
    "c_gamma_sq" -> cGammaSq,
    "mu_R_4D" -> dimAdd[M, dimPow[L, -2], dimPow[T, -2]]
  |>
];

computeFlatBraneDof[opts_Association] := Module[
  {
    removeLongitudinal, pKineticPresent, pFrankPresent, pRadialPresent,
    uWKineticPresent, uWGapPresent, phiPresent, kVec, k2, pLong, pTrans,
    uKinetic, curlStiffness, uKineticRank, uCurlRank, uCurlNullity,
    eyeP, zeroP, tangentP, radialP, pKineticForm, pFrankForm,
    pRadialHessian, pTangentKineticRank, pTangentFrankRank,
    pRadialKineticRank, pRadialHessianRank, uWKineticRank, uWGapRank,
    phiKineticRank, counts
  },
  removeLongitudinal = Lookup[opts, "remove_u_longitudinal", False];
  pKineticPresent = Lookup[opts, "p_kinetic_present", True];
  pFrankPresent = Lookup[opts, "p_frank_present", True];
  pRadialPresent = Lookup[opts, "p_radial_present", True];
  uWKineticPresent = Lookup[opts, "u_w_kinetic_present", True];
  uWGapPresent = Lookup[opts, "u_w_gap_present", True];
  phiPresent = Lookup[opts, "phi_present", False];

  kVec = {1, 2, 3};
  k2 = kVec.kVec;
  pLong = Outer[Times, kVec, kVec]/k2;
  pTrans = IdentityMatrix[3] - pLong;
  uKinetic = If[TrueQ[removeLongitudinal], pTrans, IdentityMatrix[3]];
  curlStiffness = k2 pTrans;
  uKineticRank = MatrixRank[uKinetic];
  uCurlRank = MatrixRank[uKinetic . curlStiffness . uKinetic];
  uCurlNullity = uKineticRank - uCurlRank;
  If[uCurlNullity < 0, raise["u curl rank exceeds active u kinetic rank"]];

  eyeP = IdentityMatrix[4];
  zeroP = ConstantArray[0, {4, 4}];
  tangentP = DiagonalMatrix[{1, 1, 1, 0}];
  radialP = DiagonalMatrix[{0, 0, 0, 1}];
  pKineticForm = If[TrueQ[pKineticPresent], eyeP, zeroP];
  pFrankForm = If[TrueQ[pFrankPresent], eyeP, zeroP];
  pRadialHessian = If[TrueQ[pRadialPresent], radialP, zeroP];
  pTangentKineticRank = MatrixRank[tangentP . pKineticForm . tangentP];
  pTangentFrankRank = MatrixRank[tangentP . pFrankForm . tangentP];
  pRadialKineticRank = MatrixRank[radialP . pKineticForm . radialP];
  pRadialHessianRank = MatrixRank[radialP . pRadialHessian . radialP];
  uWKineticRank = MatrixRank[{{If[TrueQ[uWKineticPresent], 1, 0]}}];
  uWGapRank = MatrixRank[{{If[TrueQ[uWGapPresent], 1, 0]}}];
  phiKineticRank = MatrixRank[{{If[TrueQ[phiPresent], 1, 0]}}];

  counts = <|
    "u_transverse" -> uCurlRank,
    "u_longitudinal" -> uCurlNullity,
    "P_spin_wave" -> Min[pTangentKineticRank, pTangentFrankRank],
    "P_soft_spin_radial" -> Min[pRadialKineticRank, pRadialHessianRank],
    "u_w" -> Min[uWKineticRank, uWGapRank],
    "phi" -> phiKineticRank
  |>;
  <|
    "Counts" -> counts,
    "Total" -> Total[Values[counts]],
    "Bookkeeping" -> <|
      "u_kinetic_rank" -> uKineticRank,
      "u_curl_rank" -> uCurlRank,
      "u_curl_nullity_within_active_kinetic_space" -> uCurlNullity,
      "P_tangent_kinetic_rank" -> pTangentKineticRank,
      "P_tangent_Frank_rank" -> pTangentFrankRank,
      "P_radial_kinetic_rank" -> pRadialKineticRank,
      "P_radial_soft_spin_hessian_rank" -> pRadialHessianRank,
      "u_w_kinetic_rank" -> uWKineticRank,
      "u_w_gap_rank" -> uWGapRank,
      "phi_kinetic_rank" -> phiKineticRank
    |>
  |>
];

runFlatBraneDof[] := Module[{computed, reportedCounts, expected, pMutated},
  subheading["HISTORICAL freeze-as-run flat-brane DOF rank computation"];
  computed = computeFlatBraneDof[<||>];
  reportedCounts = computed["Counts"];
  expectZero["projector curl stiffness rank transverse to generic k_a is 2", computed["Bookkeeping"]["u_curl_rank"] - 2];
  expectZero["projector curl stiffness nullity is one k-parallel u component", computed["Bookkeeping"]["u_curl_nullity_within_active_kinetic_space"] - 1];
  expectZero["flat-brane total DOF rank-computes to frozen value 8", computed["Total"] - 8];
  expected = <|"u_transverse" -> 2, "u_longitudinal" -> 1, "P_spin_wave" -> 3, "P_soft_spin_radial" -> 1, "u_w" -> 1, "phi" -> 0|>;
  KeyValueMap[expectZero["G0.5 breakdown " <> #1 <> " reported->computed", reportedCounts[#1] - #2] &, expected];
  Print["  structural-postulate-6 fact: C5 phi is absent in the active baseline, so phi contributes 0 DOF."];
  Print["  rank bookkeeping: ", ToString[InputForm[computed["Bookkeeping"]]]];
  Print["  computed G0.5 counts: ", ToString[InputForm[reportedCounts]], "; total=", computed["Total"]];

  pMutated = computeFlatBraneDof[<|"p_radial_present" -> False|>];
  expectZero["HISTORICAL DOF ablation drop_P_soft_spin_radial_term computes total 7", pMutated["Total"] - 7];
  expectFail["HISTORICAL DOF ablation drop_P_soft_spin_radial_term changes computed total away from 8", pMutated["Total"] - computed["Total"]];
  expectFail["historical-integrity tooth: freeze-as-run DOF cannot be falsified downward to operative 4", computed["Total"] - 4];
  computed
];

computePostD16Dof[opts_Association] := Module[
  {
    removeLongitudinal, uWGapPresent, uWKineticPresent, phiPresent, pReinject,
    kVec, k2, pLong, pTrans, uKinetic, curlForm, uKineticRank, uCurlRank,
    uNullSpaceDimension, uCurlNullity, uWKineticRank, uWGapRank, phiRank,
    pRank, counts
  },
  removeLongitudinal = Lookup[opts, "remove_u_longitudinal", False];
  uWKineticPresent = Lookup[opts, "u_w_kinetic_present", True];
  uWGapPresent = Lookup[opts, "u_w_gap_present", True];
  phiPresent = Lookup[opts, "phi_present", False];
  pReinject = Lookup[opts, "p_reinject_form", None];
  kVec = {1, 2, 3};
  k2 = kVec.kVec;
  pLong = Outer[Times, kVec, kVec]/k2;
  pTrans = IdentityMatrix[3] - pLong;
  uKinetic = If[TrueQ[removeLongitudinal], pTrans, IdentityMatrix[3]];
  curlForm = FullSimplify[uKinetic . (k2 pTrans) . uKinetic];
  uKineticRank = MatrixRank[uKinetic];
  uCurlRank = MatrixRank[curlForm];
  uNullSpaceDimension = Length[NullSpace[curlForm]];
  uCurlNullity = uKineticRank - uCurlRank;
  If[uCurlNullity < 0, raise["operative u curl rank exceeds active kinetic rank"]];
  uWKineticRank = MatrixRank[{{If[TrueQ[uWKineticPresent], 1, 0]}}];
  uWGapRank = MatrixRank[{{If[TrueQ[uWGapPresent], 1, 0]}}];
  phiRank = MatrixRank[{{If[TrueQ[phiPresent], 1, 0]}}];
  pRank = If[pReinject === None, 0, MatrixRank[pReinject]];
  counts = <|
    "u_transverse" -> uCurlRank,
    "u_longitudinal" -> uCurlNullity,
    "u_w" -> Min[uWKineticRank, uWGapRank],
    "phi" -> phiRank,
    "reinjected_P" -> pRank
  |>;
  <|
    "Counts" -> counts,
    "Total" -> Total[Values[counts]],
    "Bookkeeping" -> <|
      "u_kinetic_rank" -> uKineticRank,
      "u_curl_rank" -> uCurlRank,
      "ambient_curl_nullspace_dimension" -> uNullSpaceDimension,
      "u_curl_nullity_within_active_kinetic_space" -> uCurlNullity,
      "u_w_kinetic_rank" -> uWKineticRank,
      "u_w_gap_rank" -> uWGapRank,
      "phi_rank" -> phiRank,
      "reinjected_P_rank" -> pRank
    |>
  |>
];

runPostD16Dof[] := Module[
  {computed, expected, operativeTarget, oneMode, fullBlock, survivorMutations},
  subheading["OPERATIVE post-Decision-16 flat-brane DOF rank computation (P field removed)"];
  computed = computePostD16Dof[<||>];
  expected = <|"u_transverse" -> 2, "u_longitudinal" -> 1, "u_w" -> 1, "phi" -> 0, "reinjected_P" -> 0|>;
  KeyValueMap[expectZero["post-D16 breakdown " <> #1 <> " reported->rank-computed", computed["Counts"][#1] - #2] &, expected];
  operativeTarget = Total[Values[expected]];
  expectZero["post-D16 operative DOF rank-computes to 4", computed["Total"] - operativeTarget];
  expectZero["post-D16 computed target is the required operative DOF=4", operativeTarget - 4];
  Print["  operative rank bookkeeping: ", ToString[InputForm[computed["Bookkeeping"]]]];
  Print["  operative counts: ", ToString[InputForm[computed["Counts"]]], "; operative DOF=", computed["Total"]];
  Print["  removed P block = historical tangent 3 + radial 1 = 4 DOF; historical 8 -> operative 4."];

  oneMode = computePostD16Dof[<|"p_reinject_form" -> IdentityMatrix[1]|>];
  expectZero["operative DOF tooth fixture: re-inject one retired P mode computes 5", oneMode["Total"] - 5];
  expectFail["operative DOF tooth: one retired P mode changes operative total away from 4", oneMode["Total"] - computed["Total"]];
  fullBlock = computePostD16Dof[<|"p_reinject_form" -> IdentityMatrix[4]|>];
  expectZero["operative DOF tooth fixture: re-inject full retired P block computes 8", fullBlock["Total"] - 8];
  expectFail["operative DOF tooth: full retired P block changes operative total away from 4", fullBlock["Total"] - computed["Total"]];

  survivorMutations = {
    {"drop_u_w_gap_term", computePostD16Dof[<|"u_w_gap_present" -> False|>]},
    {"zero_u_longitudinal_component", computePostD16Dof[<|"remove_u_longitudinal" -> True|>]}
  };
  Do[
    expectZero["OPERATIVE DOF ablation " <> mutation[[1]] <> " computes total 3", mutation[[2]]["Total"] - 3];
    expectFail["OPERATIVE DOF ablation " <> mutation[[1]] <> " changes computed total away from 4", mutation[[2]]["Total"] - computed["Total"]],
    {mutation, survivorMutations}
  ];
  expectZero["baseline operative DOF remains 4 after copy-mutation teeth", computed["Total"] - 4];
  computed
];

structuralPostulates[] := {
  <|"Key" -> "postulate_1", "Text" -> "imposed `ŵ` axis + `w=0` surface (conceded-wall)"|>,
  <|"Key" -> "postulate_2", "Text" -> "`uᵃ` same-medium surface collective, tangentially free-slip (`u̇ᵃ ≠ vᵃ`)"|>,
  <|"Key" -> "postulate_3", "Text" -> "T0 `Pⁱ` reused as the Cosserat micro-rotation reservoir (0 new DOF)"|>,
  <|"Key" -> "postulate_4", "Text" -> "baseline `Pⁱ` spin-wave status = `massless` (alternates `gapped`/`slaved-rigid` named-inactive)"|>,
  <|"Key" -> "postulate_5", "Text" -> "the `ŵ`-dependent parity-EVEN P–u operator re-admits the ε-contracted/chiral class excluded by T0 and REQUIRES the conceded axis `ŵ` (a structural-postulate cost, not a free choice; the direct `P_∥·Ω_u` is parity-ODD, excluded)"|>,
  <|"Key" -> "postulate_6", "Text" -> "no C5 `φ` analog / no longitudinal constraint"|>
};

constantKeys = {"rho_br", "mu_R", "lambda_Pu", "Omega_w"};
functionKeys = {"g_ell"};
historicalPostulateKeys = Lookup[structuralPostulates[], "Key"];
retiredDriftKeys = {"lambda_Pu", "postulate_3", "postulate_4", "postulate_5"};
antiAbsorptionNames = {"rho_B0", "chi_c", "C_hu"};
validCategories = {"constant", "function", "structural_postulate"};

driftTable[] := Module[{L, T, M, posts},
  L = dim[1, 0, 0]; T = dim[0, 1, 0]; M = dim[0, 0, 1];
  posts = structuralPostulates[];
  Join[
    {
      <|"Key" -> "rho_br", "Name" -> "rho_br", "Category" -> "constant", "Dim" -> dimAdd[M, dimPow[L, -3]], "Note" -> "surface inertia"|>,
      <|"Key" -> "mu_R", "Name" -> "mu_R", "Category" -> "constant", "Dim" -> dimAdd[M, dimPow[L, -1], dimPow[T, -2]], "Note" -> "MacCullagh modulus"|>,
      <|"Key" -> "lambda_Pu", "Name" -> "lambda_Pu", "Category" -> "constant", "Dim" -> dimAdd[M, dimPow[L, -1], dimPow[T, -2]], "Note" -> "parity-repaired P-u coupling"|>,
      <|"Key" -> "Omega_w", "Name" -> "Omega_w", "Category" -> "constant", "Dim" -> dimPow[T, -1], "Note" -> "bare u_w gap scale"|>,
      <|"Key" -> "g_ell", "Name" -> "g_ell(w)", "Category" -> "function", "Dim" -> dimPow[L, -1], "Note" -> "fixed Gaussian shape, ONE width knob ell_g; no free-form profile"|>
    },
    (<|"Key" -> #["Key"], "Name" -> #["Text"], "Category" -> "structural_postulate", "Note" -> "verbatim historical structural postulate"|> &) /@ posts
  ]
];

categoryCounts[table_] := Counts[Lookup[table, "Category"]];
driftExpectedN[] := ToExpression[StringSplit[expectedDriftToken, {"(", ")"}][[2]]];
setResidual[actual_, expected_] := If[Sort[actual] === Sort[expected], 0, 1];

historicalDriftResidual[table_, verdictNDelta_: 0] := Module[
  {counts, keys, names, categories, residual, n, verdict, constantActual, functionActual, postulateActual},
  counts = categoryCounts[table];
  keys = Lookup[table, "Key"];
  names = Lookup[table, "Name"];
  categories = Lookup[table, "Category"];
  constantActual = Lookup[Select[table, #["Category"] === "constant" &], "Key"];
  functionActual = Lookup[Select[table, #["Category"] === "function" &], "Key"];
  postulateActual = Lookup[Select[table, #["Category"] === "structural_postulate" &], "Key"];
  residual = 0;
  residual += If[DuplicateFreeQ[keys], 0, 1];
  residual += If[Complement[categories, validCategories] === {}, 0, 1];
  residual += (Lookup[counts, "constant", 0] - Length[constantKeys])^2;
  residual += (Lookup[counts, "function", 0] - Length[functionKeys])^2;
  residual += (Lookup[counts, "structural_postulate", 0] - Length[structuralPostulates[]])^2;
  residual += setResidual[constantActual, constantKeys];
  residual += setResidual[functionActual, functionKeys];
  residual += setResidual[postulateActual, historicalPostulateKeys];
  residual += If[Intersection[names, antiAbsorptionNames] === {}, 0, 1];
  n = Total[Lookup[counts, validCategories, 0]];
  residual += (n - driftExpectedN[])^2;
  verdict = "SECOND_MEDIUM_DRIFT_AT_FREEZE(" <> ToString[n + verdictNDelta] <> ")";
  residual += If[verdict === expectedDriftToken, 0, 1];
  FullSimplify[residual]
];

postD16DriftResidual[historical_, operative_, retired_, verdictNDelta_: 0] := Module[
  {
    historicalMap, operativeMap, retiredMap, historicalKeys, operativeKeys,
    retiredKeys, expectedOperativeKeys, expectedEntries, expectedGroups,
    operativeGroups, categories, names, residual, n, verdict
  },
  historicalMap = AssociationThread[Lookup[historical, "Key"], historical];
  operativeMap = AssociationThread[Lookup[operative, "Key"], operative];
  retiredMap = AssociationThread[Lookup[retired, "Key"], retired];
  historicalKeys = Keys[historicalMap];
  operativeKeys = Keys[operativeMap];
  retiredKeys = Keys[retiredMap];
  expectedOperativeKeys = Complement[historicalKeys, retiredDriftKeys];
  expectedEntries = Lookup[historicalMap, expectedOperativeKeys];
  expectedGroups = GroupBy[expectedEntries, #["Category"] &, Length];
  operativeGroups = GroupBy[operative, #["Category"] &, Length];
  categories = Lookup[operative, "Category"];
  names = Lookup[operative, "Name"];
  residual = 0;
  residual += If[DuplicateFreeQ[Lookup[historical, "Key"]], 0, 1];
  residual += If[DuplicateFreeQ[Lookup[operative, "Key"]], 0, 1];
  residual += If[DuplicateFreeQ[Lookup[retired, "Key"]], 0, 1];
  residual += If[Complement[categories, validCategories] === {}, 0, 1];
  residual += setResidual[retiredKeys, retiredDriftKeys];
  residual += If[Intersection[operativeKeys, retiredKeys] === {}, 0, 1];
  residual += setResidual[Union[operativeKeys, retiredKeys], historicalKeys];
  residual += setResidual[operativeKeys, expectedOperativeKeys];
  Do[residual += If[SameQ[operativeMap[key], historicalMap[key]], 0, 1], {key, Intersection[operativeKeys, historicalKeys]}];
  Do[residual += If[SameQ[retiredMap[key], historicalMap[key]], 0, 1], {key, Intersection[retiredKeys, historicalKeys]}];
  Do[
    residual += (Lookup[operativeGroups, category, 0] - Lookup[expectedGroups, category, 0])^2,
    {category, validCategories}
  ];
  residual += If[Intersection[names, antiAbsorptionNames] === {}, 0, 1];
  n = Total[Lookup[operativeGroups, validCategories, 0]];
  residual += (n - Length[expectedOperativeKeys])^2;
  verdict = "POST_D16_DRIFT(" <> ToString[n + verdictNDelta] <> ")";
  residual += If[verdict === expectedPostD16DriftToken, 0, 1];
  FullSimplify[residual]
];

driftDimResidual[table_, dims_] := Module[{residual, targetKey},
  residual = 0;
  Do[
    If[MemberQ[{"constant", "function"}, entry["Category"]],
      targetKey = Lookup[driftDimTargets, entry["Name"], None];
      residual += If[
        targetKey === None || ! KeyExistsQ[entry, "Dim"] || ! KeyExistsQ[dims, targetKey],
        1,
        dimResidual[entry["Dim"], dims[targetKey]]
      ]
    ],
    {entry, table}
  ];
  FullSimplify[residual]
];

runComputedDriftLedger[dims_] := Module[
  {
    table, counts, constantSubcount, functionSubcount, structuralSubcount,
    n, verdict, fields, fieldSubcount, keptT0, t0NewCount, dropped,
    miscategorized, injected, dimCorrupted, retired, operative, operativeGroups,
    expectedOperative, expectedGroups, operativeN, operativeVerdict, note,
    lambdaEntry, lambdaLeftLive, postLeftLive, sameCardinalitySwap,
    operativeInjected, operativeDropped, operativeMiscategorized,
    operativeDimCorrupted
  },
  subheading["HISTORICAL freeze-as-run drift ledger enumeration"];
  table = driftTable[];
  counts = categoryCounts[table];
  constantSubcount = Lookup[counts, "constant", 0];
  functionSubcount = Lookup[counts, "function", 0];
  structuralSubcount = Lookup[counts, "structural_postulate", 0];
  n = Total[Lookup[counts, validCategories, 0]];
  verdict = "SECOND_MEDIUM_DRIFT_AT_FREEZE(" <> ToString[n] <> ")";
  Print["  Enumerated drift members:"];
  Do[
    Print[
      "    - ", entry["Name"], ": ", entry["Category"],
      If[KeyExistsQ[entry, "Dim"], ", [" <> dimString[entry["Dim"]] <> "]", ""],
      "; ", entry["Note"]
    ],
    {entry, table}
  ];
  expectZero["historical constant subcount computed from enumeration", constantSubcount - Length[constantKeys]];
  expectZero["historical function subcount computed from enumeration", functionSubcount - Length[functionKeys]];
  expectZero["historical structural-postulate subcount computed from enumeration", structuralSubcount - Length[structuralPostulates[]]];
  expectZero["independent new input count n computed from subcounts", n - driftExpectedN[]];
  expectZero["verdict string built from computed n equals frozen token", If[verdict === expectedDriftToken, 0, 1]];
  expectZero["historical drift table anti-absorption and exact-enumeration guard", historicalDriftResidual[table]];
  expectZero["historical drift table Dim fields match dimensional-firewall targets", driftDimResidual[table, dims]];
  Print["  computed n = ", constantSubcount, "+", functionSubcount, "+", structuralSubcount, " = ", n];
  Print["  computed verdict = ", verdict];

  fields = {
    <|"Name" -> "u^a", "Components" -> {"u_x", "u_y", "u_z"}|>,
    <|"Name" -> "u_w", "Components" -> {"u_w"}|>
  };
  fieldSubcount = Total[Length /@ Lookup[fields, "Components"]];
  expectZero["new-field subcount computed separately from field list", fieldSubcount - 4];
  expectBool["new-field names are kept out of the 11-input drift table", Intersection[Lookup[fields, "Name"], Lookup[table, "Name"]] === {}];
  Print["  new-field subcount = ", fieldSubcount, " from u^a (3) + u_w (1), separate from the 11-input drift count."];

  keptT0 = {
    <|"Name" -> "m rho a^2", "New" -> False|>,
    <|"Name" -> "m rho c_s^2 a^2", "New" -> False|>,
    <|"Name" -> "m rho c_s^2", "New" -> False|>
  };
  t0NewCount = Count[Lookup[keptT0, "New"], True];
  expectZero["T0 couple-stress coefficients contribute 0 new entries", t0NewCount];
  Print["  T0 couple-stress coefficients are kept-not-new: m rho a^2; m rho c_s^2 a^2; m rho c_s^2."];

  dropped = Most[table];
  expectFail["HISTORICAL enumeration tooth: drop one entry gives n=10 and trips drift validation", historicalDriftResidual[dropped]];
  miscategorized = (If[#["Name"] === "Omega_w", Join[KeyDrop[#, "Category"], <|"Category" -> "structural_postulate"|>], #] &) /@ table;
  expectFail["HISTORICAL enumeration tooth: miscategorize Omega_w trips subcount assertions", historicalDriftResidual[miscategorized]];
  injected = Append[table, <|"Key" -> "rho_B0", "Name" -> "rho_B0", "Category" -> "constant", "Dim" -> table[[1, "Dim"]], "Note" -> "forbidden Part-VI injection"|>];
  expectFail["HISTORICAL enumeration tooth: inject rho_B0 trips anti-absorption guard", historicalDriftResidual[injected]];
  expectFail["HISTORICAL enumeration tooth: corrupt computed n before verdict assembly trips token equality", historicalDriftResidual[table, 1]];
  dimCorrupted = (If[#["Name"] === "rho_br", Join[KeyDrop[#, "Dim"], <|"Dim" -> dimAdd[#["Dim"], dim[1, 0, 0]]|>], #] &) /@ table;
  expectFail["HISTORICAL enumeration tooth: corrupt rho_br table Dim trips firewall consistency", driftDimResidual[dimCorrupted, dims]];
  expectFail["historical-integrity tooth: freeze-as-run drift cannot be falsified downward to 7", n - 7];
  expectZero["baseline historical drift table remains valid after copy-mutation teeth", historicalDriftResidual[table]];

  subheading["OPERATIVE post-Decision-16 drift derived as historical minus exact retired set"];
  retired = Select[table, MemberQ[retiredDriftKeys, #["Key"]] &];
  operative = Select[table, ! MemberQ[retiredDriftKeys, #["Key"]] &];
  operativeGroups = GroupBy[operative, #["Category"] &, Length];
  expectedOperative = Select[table, ! MemberQ[retiredDriftKeys, #["Key"]] &];
  expectedGroups = GroupBy[expectedOperative, #["Category"] &, Length];
  operativeN = Total[Lookup[operativeGroups, validCategories, 0]];
  operativeVerdict = "POST_D16_DRIFT(" <> ToString[operativeN] <> ")";
  Print["  Enumerated operative survivors:"];
  Do[
    note = If[
      entry["Key"] === "postulate_1",
      "live annotation softened: intrinsic retained wall normal and w=0 geometry; no longer an extra axis conceded for the retired P-u operator",
      entry["Note"]
    ];
    Print[
      "    - ", entry["Key"], ": ", entry["Name"], ": ", entry["Category"],
      If[KeyExistsQ[entry, "Dim"], ", [" <> dimString[entry["Dim"]] <> "]", ""],
      "; ", note
    ],
    {entry, operative}
  ];
  Print["  Postulate-(1) adjudication: RETAIN the wall-normal geometry, but soften 'conceded-wall' to intrinsic wall normal because L_Pu is gone; this changes annotation, not membership or count."];
  Print["  exact retired keys = {", StringRiffle[Lookup[retired, "Key"], ", "], "}"];
  Do[
    expectZero[
      "post-D16 " <> category <> " subcount computed from survivor enumeration",
      Lookup[operativeGroups, category, 0] - Lookup[expectedGroups, category, 0]
    ],
    {category, validCategories}
  ];
  expectZero["post-D16 operative n computes as 3 constants + 1 function + 3 structural postulates", operativeN - 7];
  expectZero["post-D16 drift is exact historical-minus-retired set partition", postD16DriftResidual[table, operative, retired]];
  expectZero["post-D16 drift table Dim fields match live firewall survivors", driftDimResidual[operative, dims]];
  expectZero["post-D16 verdict string built from computed n", If[operativeVerdict === expectedPostD16DriftToken, 0, 1]];
  Print["  computed n = ", Lookup[operativeGroups, "constant", 0], "+", Lookup[operativeGroups, "function", 0], "+", Lookup[operativeGroups, "structural_postulate", 0], " = ", operativeN];
  Print["  computed operative verdict = ", operativeVerdict];

  lambdaEntry = First[Select[retired, #["Key"] === "lambda_Pu" &]];
  lambdaLeftLive = Append[operative, lambdaEntry];
  expectZero["post-D16 drift tooth fixture: leave lambda_Pu live computes n=8", Length[lambdaLeftLive] - 8];
  expectFail["post-D16 drift tooth: leave lambda_Pu live trips DRIFT(7)", postD16DriftResidual[table, lambdaLeftLive, retired]];
  Do[
    postLeftLive = Append[operative, First[Select[retired, #["Key"] === postKey &]]];
    expectZero["post-D16 drift tooth fixture: leave " <> postKey <> " live computes n=8", Length[postLeftLive] - 8];
    expectFail["post-D16 drift tooth: leave " <> postKey <> " live trips DRIFT(7)", postD16DriftResidual[table, postLeftLive, retired]],
    {postKey, {"postulate_3", "postulate_4", "postulate_5"}}
  ];
  sameCardinalitySwap = Append[Select[operative, #["Key"] =!= "Omega_w" &], lambdaEntry];
  expectZero["post-D16 drift swap fixture retains cardinality n=7", Length[sameCardinalitySwap] - 7];
  expectFail["post-D16 drift tooth: same-cardinality Omega_w/lambda_Pu swap trips exact set partition", postD16DriftResidual[table, sameCardinalitySwap, retired]];
  operativeInjected = Append[operative, <|"Key" -> "rho_B0", "Name" -> "rho_B0", "Category" -> "constant", "Dim" -> table[[1, "Dim"]], "Note" -> "forbidden Part-VI injection"|>];
  expectFail["post-D16 drift tooth: inject rho_B0 trips operative anti-absorption guard", postD16DriftResidual[table, operativeInjected, retired]];
  expectFail["post-D16 drift tooth: corrupt n before token assembly trips POST_D16_DRIFT equality", postD16DriftResidual[table, operative, retired, 1]];

  operativeDropped = Most[operative];
  expectZero["OPERATIVE enumeration tooth fixture: drop one survivor computes n=6", Length[operativeDropped] - 6];
  expectFail["OPERATIVE enumeration tooth: drop one survivor trips n=7 and set partition", postD16DriftResidual[table, operativeDropped, retired]];
  operativeMiscategorized = (If[#["Key"] === "Omega_w", Join[KeyDrop[#, "Category"], <|"Category" -> "structural_postulate"|>], #] &) /@ operative;
  expectFail["OPERATIVE enumeration tooth: miscategorize Omega_w trips survivor subcount", postD16DriftResidual[table, operativeMiscategorized, retired]];
  operativeDimCorrupted = (If[#["Key"] === "rho_br", Join[KeyDrop[#, "Dim"], <|"Dim" -> dimAdd[#["Dim"], dim[1, 0, 0]]|>], #] &) /@ operative;
  expectFail["OPERATIVE enumeration tooth: corrupt rho_br table Dim trips live firewall consistency", driftDimResidual[operativeDimCorrupted, dims]];
  expectZero["baseline post-D16 drift remains valid after copy-mutation teeth", postD16DriftResidual[table, operative, retired]];
  {table, operative}
];

distinctnessFailureResidual[muR_, muR4D_] := If[dimResidual[muR, muR4D] === 0, 1, 0];

runMuRFirewall[dims_] := Module[{muR, muR4D, L, forcedEqualDims},
  subheading["mu_R notational firewall and R17 pending edge"];
  L = dim[1, 0, 0];
  muR = dims["mu_R"];
  muR4D = dims["mu_R_4D"];
  Print["  [mu_R] 3D brane modulus = ", dimString[muR]];
  Print["  [mu_R^(4)] 4D shear-stiffness density = ", dimString[muR4D]];
  expectNonzero["[mu_R] != [mu_R^(4)] as exponent triples", dimResidual[muR, muR4D]];
  expectDim["R17 dim consistency [mu_R^(4)]*L = [mu_R]", dimAdd[muR4D, L], muR];
  Print["  R17 status: PENDING (mu_R = int chi_B mu_R^(4) dw; deferred nonlinear throat/projection)."];
  forcedEqualDims = Association[dims];
  forcedEqualDims["mu_R_4D"] = forcedEqualDims["mu_R"];
  expectFail[
    "mu_R firewall tooth: forced equality trips distinctness inequality",
    distinctnessFailureResidual[forcedEqualDims["mu_R"], forcedEqualDims["mu_R_4D"]]
  ]
];

printStructuralPostulates[] := (
  subheading["HISTORICAL six structural postulates printed verbatim"];
  Do[
    Print["  ", i, ". ", structuralPostulates[][[i, "Key"]], ": ", structuralPostulates[][[i, "Text"]]],
    {i, Length[structuralPostulates[]]}
  ];
  Print["  OPERATIVE survivors: postulate_1 (intrinsic wall-normal annotation), postulate_2, postulate_6; postulates 3/4/5 retired with P."]
);

printProvenanceBlocks[] := (
  subheading["Erratum, supersession, methodology, Gate-L scope, and carried debt"];
  Print["  ERRATUM (2026-07-04): the SECOND_MEDIUM_DRIFT_AT_FREEZE(11) count is NOT inflated by a rho_br overcount; this freeze's count STANDS."];
  Print["  pathA_25 varrho_br[rho] belongs to the CLOSED density-smectic candidate, FAIL_NOT_CODIM1, OUT_OF_ACTIVE_NG5."];
  Print["  This rho_br/mu_R is genuine postulated shear-surface inertia/modulus with registered-pending pathA_40 Route-A reduction."];
  Print["  Corroboration token: NO_OVERCOUNT_ROUTE_A_PENDING."];
  Print["  Honest cross-sector drift {rho_B0, chi_c, C_hu} is a Part-VI item, not absorbed into the historical 11 or operative 7."];
  Print[""];
  Print["  Supersession fact 1: stage006 chi_B order-field wall superseded fixed-shape g_ell(w) as the MATERIAL-STATE closure."];
  Print["  Supersession fact 2: G0 REMAINS the light-sector CONSTITUTIVE freeze; stage003 consumes L_Mac as-is."];
  Print[""];
  Print["  Methodology: postulating an ingredient is allowed; postulating an outcome is not."];
  Print["  Methodology: late ingredient = AD_HOC_RESCUE -> fresh G0; every knob counted; >=2 new inputs => drift reported plainly."];
  Print["  Methodology: g(w) admitted on locality/minimality grounds ONLY, target-blind."];
  Print["  Methodology anti-impose: G0 freezes TERMS, not gate answers; no bounded-below, traction, or longitudinal-is-gauge claims."];
  Print["  Methodology: a clean all-pass is suspicious, so able-to-fail teeth are live."];
  Print[""];
  Print["  Scope: G0 freeze only; no Gate L verdict computed."];
  Print["  classification_guard: counts only; no gate verdict, no boundedness, no gauge, no leak claim."];
  Print["  Gate-L exposure names are prose provenance only: FAIL_HIDDEN_PROPAGATING_MODE; FAIL_GYROSTAT_NO_CLOSURE; FAIL_NOT_BOUNDED_BELOW; linked FAIL_COUPLE_STRESS_NOGO chain remains able to fire."];
  Print["  Gate-L artifacts are not imported; exposure strings are not used as computed predicates."];
  Print[""];
  Print["  Route-A carried debt: {rho_br, mu_R} reduction = Route-A PENDING (R10), ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT, free-unreduced brane constants on the deferred nonlinear throat."];
  Print["  Historical labeled postulated constants: {lambda_Pu, Omega_w, ell_g}; Decision 16 retires lambda_Pu, leaving live {Omega_w, ell_g} alongside survivor {rho_br, mu_R}."];
  Print[""];
  Print["  Downstream consumers: ledger_stage003 consumes mu_R, rho_br, c_gamma^2=mu_R/rho_br, and frozen L_Mac."];
  Print["  Downstream consumers: ledger_stage006 cites c_gamma^2=mu_R/rho_br plus the chi_B/G0 supersession relationship."]
);

printVerdictLabels[] := (
  Print[""];
  Print["Verdict labels:"];
  Print["  ledger earned-label (NOT a source verdict token): G0_FREEZE_FIDELITY_PLUS_POST_D16_LAYER_VERIFIED  (historical hash/DOF/drift preserved; operative action partition, DOF, drift, and tiered firewall computed)"];
  Print["  HISTORICAL freeze-as-run immutable: ", expectedFreezeToken, " + ", expectedDriftToken, " + historical DOF=8"];
  Print["  OPERATIVE post-Decision-16 live: ", expectedPostD16ActionToken, " + ", expectedPostD16DriftToken, " + operative DOF=4"];
  Print["  operative action core token: POST_D16_ACTION{S_GNLS,gL_Mac,gL_uw}; exact partition of historical summands with retired complement {L_pol,gL_Pu}"];
  Print["  historical landing: 11 = 4 constants {rho_br, mu_R, lambda_Pu, Omega_w} + 1 function g_l(w; l_g) + 6 structural postulates; operative 7 = 3 constants {rho_br,mu_R,Omega_w} + 1 function + survivor postulates {1,2,6}"];
  Print["  erratum (2026-07-04): historical 11 STANDS; operative 7 is Decision-16 retirement, not an overcount correction; {rho_B0, chi_c, C_hu} remains Part-VI and is excluded from both tables [NO_OVERCOUNT_ROUTE_A_PENDING]"];
  Print["  DECISION16_PROVENANCE retired={L_pol,L_Pu,lambda_Pu,postulates_3/4/5}; reason=P_RETIRED_ALL_PAYOFFS_FAILED_PLUS_LIFSHITZ_INSTABILITY"];
  Print["  supersession: stage006 chi_B wall = the MATERIAL-STATE closure (supersedes fixed-shape g_l(w) as material wall); G0 REMAINS the light-sector CONSTITUTIVE freeze (stage003 consumes L_Mac as-is)"];
  Print["  notational firewall: mu_R (3D brane, M L^-1 T^-2) != mu_R_4D (4D density, M L^-2 T^-2); related only by PENDING R17 projection"];
  Print["  Gate-L: EXCLUDED — no gate verdict computed or imported; exposure names printed as provenance only"];
  Print["  carried: Route-A reduction PENDING (pathA_40, R10) for {rho_br, mu_R}; live postulated Omega_w and l_g remain; lambda_Pu is retired-historical"]
);

Module[{ok, dims},
  heading["ledger_stage007_shear_surface_g0_freeze Mathematica audit"];
  ok = Catch[
    runFreezeFidelity[];
    runPostD16ActionPartition[];
    dims = runDimensionalFirewall[];
    runFlatBraneDof[];
    runPostD16Dof[];
    runComputedDriftLedger[dims];
    runMuRFirewall[dims];
    printStructuralPostulates[];
    printProvenanceBlocks[];
    printVerdictLabels[];
    True,
    "ledgerStage007Failure",
    Function[{msg, tag}, Print["FAIL: ", msg]; False]
  ];

  Print[""];
  Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
  If[TrueQ[ok],
    Print["OVERALL PASS: Mathematica verified historical stage007 freeze fidelity plus the computed post-Decision-16 operative layer"];
    Exit[0],
    Print["OVERALL FAIL: Mathematica stage007 audit did not close"];
    Exit[1]
  ]
]
