(* Ledger stage015 Mathematica audit: breathing structure + HF force.

   Standalone, print-only, no arguments, no file I/O. This is the pathA_31
   II-G2c slice only. Stage 013 profiles and M/K closed forms are cited as
   native literals with dual-site integrity; stage 014's truncation certificate
   is cited. The native 015 route uses D, Integrate, Det, MatrixRank, and its
   own raw-tree comparison. No Get/Import/Export bridge is present.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;

BREATHINGCALIBRATED = "BREATHING_CALIBRATED";
BREATHINGFAILSTRUCTURE = "BREATHING_FAIL_STRUCTURE";
BREATHINGFAILHFFORCE = "BREATHING_FAIL_HF_FORCE";

$Assumptions =
  L0 > 0 && beta > 0 && muEta > 0 && Tw > 0 && rAL > 0 &&
  rhoStar > 0 && Vp0 > 0 && ellC > 0 &&
  kappa > 0 && chi > 0 && sigmaA > 0 && sigmaL > 0 && B > 0 &&
  Element[{w, deltaA, deltaL, qa, qL, r}, Reals];

raise[msg_] := Throw[msg, "ledgerStage015Failure"];

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
rawString[expr_] := ToString[Unevaluated[expr], InputForm];

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

boolResidual[condition_] := If[TrueQ[condition], 0, 1];
verdictResidual[actual_, expected_] := If[actual === expected, 0, 1];
zeroPattern[m_] := Map[TrueQ[clean[#] === 0] &, m, {2}];

citedProfiles[] := <|
  "AlphaA" -> Sinh[L0 beta - beta w]/Sinh[L0 beta],
  "AlphaL" -> rAL Sinh[beta w]/Sinh[L0 beta]
|>;

citedMK[] := Module[{mEntries, kEntries, m, k},
  mEntries = <|
    "aa" -> -2 Pi muEta (L0 beta Tanh[L0 beta] - Sinh[L0 beta]^2)/(beta Sinh[L0 beta]^2 Tanh[L0 beta]),
    "aL" -> 2 Pi muEta rAL (L0 beta - Tanh[L0 beta])/(beta Sinh[L0 beta] Tanh[L0 beta]),
    "LL" -> -2 Pi muEta rAL^2 (L0 beta Tanh[L0 beta] - Sinh[L0 beta]^2)/(beta Sinh[L0 beta]^2 Tanh[L0 beta])
  |>;
  kEntries = <|
    "aa" -> 4 Pi Tw beta/Tanh[L0 beta],
    "aL" -> -4 Pi Tw beta rAL/Sinh[L0 beta],
    "LL" -> 4 Pi Tw beta rAL^2/Tanh[L0 beta]
  |>;
  m = {{mEntries["aa"], mEntries["aL"]}, {mEntries["aL"], mEntries["LL"]}};
  k = {{kEntries["aa"], kEntries["aL"]}, {kEntries["aL"], kEntries["LL"]}};
  <|"MEntries" -> mEntries, "KEntries" -> kEntries, "M" -> m, "K" -> k|>
];

detM013Form[] := 4 Pi^2 muEta^2 rAL^2 (Sinh[L0 beta]^2 - (L0 beta)^2)/(beta^2 Sinh[L0 beta]^2);
detK013Form[] := 16 Pi^2 Tw^2 beta^2 rAL^2;

positivityCertificateIdentities[] := Module[{f1, g, f2, h},
  f1 = Sinh[B] Cosh[B] - B;
  g = Sinh[B] - B;
  f2 = Sinh[B]^2 - B^2;
  h = B - Tanh[B];
  <|
    "f1_at_0" -> (f1 /. B -> 0),
    "f1_prime_square" -> D[f1, B] - 2 Sinh[B]^2,
    "g_at_0" -> (g /. B -> 0),
    "g_prime_square" -> D[g, B] - 2 Sinh[B/2]^2,
    "f2_factorization" -> f2 - g (Sinh[B] + B),
    "h_at_0" -> (h /. B -> 0),
    "h_prime_square" -> D[h, B] - Tanh[B]^2
  |>
];

entryCertificateResiduals[m_, k_] := Module[{localB, f1, f2, h},
  localB = L0 beta;
  f1 = Sinh[localB] Cosh[localB] - localB;
  f2 = Sinh[localB]^2 - localB^2;
  h = localB - Tanh[localB];
  <|
    "M_aa_core_is_f1" -> clean[m[[1, 1]] beta Sinh[localB]^2/(2 Pi muEta) - f1],
    "det_M_core_is_f2" -> clean[Det[m] beta^2 Sinh[localB]^2/(4 Pi^2 muEta^2 rAL^2) - f2],
    "M_aL_core_is_h" -> clean[m[[1, 2]] beta Sinh[localB] Tanh[localB]/(2 Pi muEta rAL) - h],
    "K_aL_negative_core" -> clean[k[[1, 2]] Sinh[localB]/(4 Pi Tw beta rAL) + 1]
  |>
];

mPosdefByCertificate[m_] := Module[{ids, residuals, checks},
  ids = positivityCertificateIdentities[];
  residuals = entryCertificateResiduals[m, IdentityMatrix[2]];
  checks = {ids["f1_at_0"], ids["f1_prime_square"], ids["g_at_0"], ids["g_prime_square"], ids["f2_factorization"], residuals["M_aa_core_is_f1"], residuals["det_M_core_is_f2"]};
  AllTrue[checks, TrueQ[clean[#] === 0] &]
];

mALPositiveByCertificate[m_] := Module[{ids, residuals, checks},
  ids = positivityCertificateIdentities[];
  residuals = entryCertificateResiduals[m, IdentityMatrix[2]];
  checks = {ids["h_at_0"], ids["h_prime_square"], residuals["M_aL_core_is_h"]};
  AllTrue[checks, TrueQ[clean[#] === 0] &]
];

kOffdiagNegativeByCertificate[k_] := Module[{residuals},
  residuals = entryCertificateResiduals[IdentityMatrix[2], k];
  TrueQ[clean[residuals["K_aL_negative_core"]] === 0]
];

mkCitationIntegrity[m_, k_] := Module[{detMResidual, detKResidual, malPositive, kalNegative},
  detMResidual = clean[Det[m] - detM013Form[]];
  detKResidual = clean[Det[k] - detK013Form[]];
  malPositive = mALPositiveByCertificate[m];
  kalNegative = kOffdiagNegativeByCertificate[k];
  <|
    "DetMResidual" -> detMResidual,
    "DetKResidual" -> detKResidual,
    "MALPositive" -> malPositive,
    "KALNegative" -> kalNegative,
    "Ok" -> TrueQ[detMResidual === 0 && detKResidual === 0 && malPositive && kalNegative]
  |>
];

structureGate[m_, k_, h_] := Module[
  {
    legacyDetPositiveForm, legacySymmetric, legacyOffdiagNegative,
    legacyDetPositive, legacyRank, legacyZeros, mPosdef, kSymmetric,
    kOffdiagNegative, kRank, kRankMatchesLegacy, kZeroPattern,
    kZeroPatternMatchesLegacy, kStructureOk
  },
  legacyDetPositiveForm = kappa sigmaA + kappa chi^2 sigmaL + sigmaA sigmaL;
  legacySymmetric = TrueQ[h === Transpose[h]];
  legacyOffdiagNegative = TrueQ[clean[h[[1, 2]] + chi kappa] === 0];
  legacyDetPositive = TrueQ[clean[Det[h] - legacyDetPositiveForm] === 0];
  legacyRank = MatrixRank[h];
  legacyZeros = zeroPattern[h];
  mPosdef = mPosdefByCertificate[m];
  kSymmetric = TrueQ[k === Transpose[k]];
  kOffdiagNegative = kOffdiagNegativeByCertificate[k];
  kRank = MatrixRank[k];
  kRankMatchesLegacy = TrueQ[kRank === legacyRank && legacyRank === 2];
  kZeroPattern = zeroPattern[k];
  kZeroPatternMatchesLegacy = TrueQ[kZeroPattern === legacyZeros];
  kStructureOk = TrueQ[
    kSymmetric && kOffdiagNegative && kRankMatchesLegacy &&
    kZeroPatternMatchesLegacy && legacySymmetric && legacyOffdiagNegative &&
    legacyDetPositive
  ];
  <|
    "MPosdef" -> mPosdef,
    "KSymmetric" -> kSymmetric,
    "KOffdiagNegative" -> kOffdiagNegative,
    "KRank" -> kRank,
    "KRankMatchesLegacy" -> kRankMatchesLegacy,
    "KZeroPattern" -> kZeroPattern,
    "KZeroPatternMatchesLegacy" -> kZeroPatternMatchesLegacy,
    "KStructureOk" -> kStructureOk,
    "LegacySymmetric" -> legacySymmetric,
    "LegacyOffdiagNegative" -> legacyOffdiagNegative,
    "LegacyDetPositive" -> legacyDetPositive,
    "LegacyRank" -> legacyRank,
    "LegacyZeroPattern" -> legacyZeros,
    "StructureFromComputedMatrices" -> True,
    "FullMatrixFit" -> False
  |>
];

structureProbes[m_, k_, h_] := Module[{mProbe, kProbe},
  mProbe = m;
  mProbe[[1, 1]] = -mProbe[[1, 1]];
  kProbe = k;
  kProbe[[1, 2]] = -kProbe[[1, 2]];
  kProbe[[2, 1]] = -kProbe[[2, 1]];
  <|
    "NonPosdefMProbe" -> structureGate[mProbe, k, h],
    "SignFlippedKProbe" -> structureGate[m, kProbe, h]
  |>
];

computeVerdict[gate_, hfForceReduces_, unsimplifiedRoutesIdentical_] := Which[
  ! TrueQ[gate["MPosdef"] && gate["KSymmetric"] && gate["KStructureOk"]], BREATHINGFAILSTRUCTURE,
  ! TrueQ[hfForceReduces] || TrueQ[unsimplifiedRoutesIdentical], BREATHINGFAILHFFORCE,
  True, BREATHINGCALIBRATED
];

packetResiduals[packet_] := <|
  "site_A_betaL0_anchor" -> clean[packet["betaL0Cited"] - 37/20],
  "site_A_packet_product" -> clean[packet["betaCited"] packet["L0Cited"] - packet["betaL0Cited"]],
  "anchor_L0" -> clean[packet["L0Cited"] - 37/20],
  "anchor_beta" -> clean[packet["betaCited"] - 1]
|>;

packetOkQ[packet_] := AllTrue[Values[packetResiduals[packet]], TrueQ[# === 0] &];
packetSet[packet_, key_, value_] := ReplacePart[packet, Key[key] -> value];

profileSiteBResiduals[aa_, al_] := <|
  "residual_alpha_a" -> clean[-D[aa, {w, 2}] + beta^2 aa],
  "residual_alpha_L" -> clean[-D[al, {w, 2}] + beta^2 al],
  "bc_alpha_a_0" -> clean[(aa /. w -> 0) - 1],
  "bc_alpha_a_L0" -> clean[aa /. w -> L0],
  "bc_alpha_L_0" -> clean[al /. w -> 0],
  "bc_alpha_L_L0" -> clean[(al /. w -> L0) - rAL]
|>;

profileSiteBOkQ[aa_, al_] := AllTrue[Values[profileSiteBResiduals[aa, al]], TrueQ[# === 0] &];

hessianEntryResiduals[h_] := <|
  "aa" -> clean[h[[1, 1]] - (chi^2 kappa + sigmaA)],
  "aL" -> clean[h[[1, 2]] + chi kappa],
  "LL" -> clean[h[[2, 2]] - (kappa + sigmaL)],
  "det" -> clean[Det[h] - (kappa sigmaA + kappa chi^2 sigmaL + sigmaA sigmaL)]
|>;

buildBaseline[] := Module[
  {
    profiles, alphaA, alphaL, mk, m, k, egeom, legacyH, gate, probes,
    src, forceDistAUns, forceDistLUns, forceDistA, forceDistL, rParam,
    vConf, forceLegacyAUns, forceLegacyLUns, forceLegacyA, forceLegacyL,
    hfAOk, hfLOk, unsimplifiedRoutesIdentical, wrongZeroFA,
    wrongConstFA, degenerateMismatch, constantMismatch, verdict
  },
  profiles = citedProfiles[];
  alphaA = profiles["AlphaA"];
  alphaL = profiles["AlphaL"];
  mk = citedMK[];
  m = mk["M"];
  k = mk["K"];
  egeom = 1/2 kappa (deltaL - chi deltaA)^2 + 1/2 sigmaA deltaA^2 + 1/2 sigmaL deltaL^2;
  legacyH = D[egeom, {{deltaA, deltaL}, 2}];
  gate = structureGate[m, k, legacyH];
  probes = structureProbes[m, k, legacyH];
  src = clean[rhoStar Vp0/ellC];
  forceDistAUns = 4 Pi Inactive[Integrate][src alphaA, {w, 0, L0}];
  forceDistLUns = 4 Pi Inactive[Integrate][src alphaL, {w, 0, L0}];
  forceDistA = clean[Activate[forceDistAUns]];
  forceDistL = clean[Activate[forceDistLUns]];
  rParam = r0[w] + qa alphaA + qL alphaL;
  vConf = (Vp0/ellC) (r - rParam);
  forceLegacyAUns = -4 Pi Inactive[Integrate][rhoStar (D[vConf, qa] /. {qa -> 0, qL -> 0}), {w, 0, L0}];
  forceLegacyLUns = -4 Pi Inactive[Integrate][rhoStar (D[vConf, qL] /. {qa -> 0, qL -> 0}), {w, 0, L0}];
  forceLegacyA = clean[Activate[forceLegacyAUns]];
  forceLegacyL = clean[Activate[forceLegacyLUns]];
  hfAOk = TrueQ[clean[forceDistA - forceLegacyA] === 0];
  hfLOk = TrueQ[clean[forceDistL - forceLegacyL] === 0];
  unsimplifiedRoutesIdentical = rawString[forceDistAUns] === rawString[forceLegacyAUns];
  wrongZeroFA = clean[4 Pi Integrate[src 0, {w, 0, L0}]];
  wrongConstFA = clean[4 Pi Integrate[src 1, {w, 0, L0}]];
  degenerateMismatch = ! TrueQ[clean[wrongZeroFA - forceLegacyA] === 0];
  constantMismatch = ! TrueQ[clean[wrongConstFA - forceLegacyA] === 0];
  verdict = computeVerdict[gate, TrueQ[hfAOk && hfLOk], unsimplifiedRoutesIdentical];
  <|
    "Profiles" -> profiles,
    "MK" -> mk,
    "EGeom" -> egeom,
    "HLegacy" -> legacyH,
    "Gate" -> gate,
    "Probes" -> probes,
    "SourceDensity" -> src,
    "HF" -> <|
      "ForceDistAUns" -> forceDistAUns,
      "ForceDistLUns" -> forceDistLUns,
      "ForceDistA" -> forceDistA,
      "ForceDistL" -> forceDistL,
      "ForceLegacyAUns" -> forceLegacyAUns,
      "ForceLegacyLUns" -> forceLegacyLUns,
      "ForceLegacyA" -> forceLegacyA,
      "ForceLegacyL" -> forceLegacyL,
      "HFAOk" -> hfAOk,
      "HFLOk" -> hfLOk,
      "HFForceReduces" -> TrueQ[hfAOk && hfLOk],
      "UnsimplifiedRoutesIdentical" -> unsimplifiedRoutesIdentical
    |>,
    "Guards" -> <|
      "WrongZeroFA" -> wrongZeroFA,
      "WrongConstFA" -> wrongConstFA,
      "DegenerateHFMismatch" -> degenerateMismatch,
      "ConstantHFMismatch" -> constantMismatch,
      "ConstantOverlapPasses014" -> True
    |>,
    "Packet" -> <|"L0Cited" -> 37/20, "betaCited" -> 1, "betaL0Cited" -> 37/20|>,
    "MKIntegrity" -> mkCitationIntegrity[m, k],
    "Verdict" -> verdict
  |>
];

runAritySelfCheck[data_] := Module[{profiles, mk, certs, packet, gate, integ},
  subheading["Wolfram arity self-check"];
  profiles = citedProfiles[];
  mk = citedMK[];
  certs = positivityCertificateIdentities[];
  packet = packetResiduals[data["Packet"]];
  gate = structureGate[mk["M"], mk["K"], data["HLegacy"]];
  integ = mkCitationIntegrity[mk["M"], mk["K"]];
  expectBool["arity citedProfiles[0 args] returns AlphaA/AlphaL", KeyExistsQ[profiles, "AlphaA"] && KeyExistsQ[profiles, "AlphaL"]];
  expectBool["arity citedMK[0 args] returns M/K matrices", KeyExistsQ[mk, "M"] && KeyExistsQ[mk, "K"]];
  expectBool["arity positivityCertificateIdentities[0 args] returns f1/f2/h identities", KeyExistsQ[certs, "f1_prime_square"] && KeyExistsQ[certs, "h_prime_square"]];
  expectBool["arity packetResiduals[assoc] returns site A keys", KeyExistsQ[packet, "site_A_betaL0_anchor"] && KeyExistsQ[packet, "site_A_packet_product"]];
  expectBool["arity structureGate[3 args] returns structure booleans", KeyExistsQ[gate, "MPosdef"] && KeyExistsQ[gate, "KStructureOk"]];
  expectBool["arity mkCitationIntegrity[2 args] returns det/sign checks", KeyExistsQ[integ, "DetMResidual"] && KeyExistsQ[integ, "MALPositive"]];
  expectBool["arity computeVerdict[3 args] returns BREATHING_CALIBRATED", computeVerdict[data["Gate"], True, False] === BREATHINGCALIBRATED]
];

runConsumedInputs[data_] := Module[
  {packet, alphaA, alphaL, m, k, integrity, residuals, kernelBad, nonkernelBad,
   kernelResiduals, mBad, kBad, signFlipIntegrity, kSignFlipIntegrity},
  subheading["CONSUMED-from-013/014 citation integrity"];
  packet = data["Packet"];
  alphaA = data["Profiles"]["AlphaA"];
  alphaL = data["Profiles"]["AlphaL"];
  m = data["MK"]["M"];
  k = data["MK"]["K"];
  integrity = data["MKIntegrity"];
  Print["  CONSUMED-from-013: profiles alpha_a,alpha_L; operator-projected M_AB/K_AB; frozen packet {L0=37/20,beta=1,beta*L0=37/20}."];
  Print["  CONSUMED-from-014: truncation certificate beta_L0 in [0.1,3.0] / K_eta/Tw <= ~2.6; cited only, no eigensolve here."];
  Print["  Site A reads betaL0Cited as an independently corruptible cited datum, not local beta*L0."];
  residuals = packetResiduals[packet];
  Scan[Function[key, expectZero["site A consumed packet " <> key, residuals[key]]], Keys[residuals]];
  expectBool["site A clean packet passes", packetOkQ[packet]];
  expectFail["site A betaL0Cited-only corruption is non-vacuous", boolResidual[packetOkQ[packetSet[packet, "betaL0Cited", 19/10]]]];
  Print["  Site B verifies cited profile residual AND endpoint BCs; it does not solve the BVP."];
  residuals = profileSiteBResiduals[alphaA, alphaL];
  Scan[Function[key, expectZero["site B profile " <> key, residuals[key]]], Keys[residuals]];
  expectBool["site B residual-and-BC guard passes cited profiles", profileSiteBOkQ[alphaA, alphaL]];
  kernelBad = Sinh[L0 beta - beta w]/Cosh[L0 beta];
  nonkernelBad = alphaA + 1;
  kernelResiduals = profileSiteBResiduals[kernelBad, alphaL];
  expectBool["kernel-preserving corruption has zero residual but wrong BC", kernelResiduals["residual_alpha_a"] === 0 && kernelResiduals["bc_alpha_a_0"] =!= 0];
  expectFail["kernel-preserving profile corruption breaks site B BC leg", boolResidual[profileSiteBOkQ[kernelBad, alphaL]]];
  expectFail["non-kernel profile corruption breaks site B residual leg", boolResidual[profileSiteBOkQ[nonkernelBad, alphaL]]];
  Print["  Consumed M/K det-identities plus off-diagonal sign checks cover every cited entry, including det-blind sign flips."];
  expectZero["consumed det(M) identity vs independent 013 closed form", integrity["DetMResidual"]];
  expectZero["consumed det(K) identity vs independent 013 closed form", integrity["DetKResidual"]];
  expectBool["consumed M_aL > 0 via B-tanh(B) certificate", integrity["MALPositive"]];
  expectBool["consumed K_aL < 0 via -1/sinh(B) certificate", integrity["KALNegative"]];
  expectBool["consumed M/K citation integrity clean baseline passes", integrity["Ok"]];
  mBad = m; mBad[[1, 1]] = 2 mBad[[1, 1]];
  expectFail["M_aa diagonal corruption breaks det(M) identity", mkCitationIntegrity[mBad, k]["DetMResidual"]];
  mBad = m; mBad[[2, 2]] = 2 mBad[[2, 2]];
  expectFail["M_LL diagonal corruption breaks det(M) identity", mkCitationIntegrity[mBad, k]["DetMResidual"]];
  mBad = m; mBad[[1, 2]] = 2 mBad[[1, 2]]; mBad[[2, 1]] = 2 mBad[[2, 1]];
  expectFail["M_aL magnitude corruption breaks det(M) identity", mkCitationIntegrity[mBad, k]["DetMResidual"]];
  mBad = m; mBad[[1, 2]] = -mBad[[1, 2]]; mBad[[2, 1]] = -mBad[[2, 1]];
  signFlipIntegrity = mkCitationIntegrity[mBad, k];
  expectZero["M_aL sign flip leaves det(M) identity blind as expected", signFlipIntegrity["DetMResidual"]];
  expectFail["M_aL sign flip breaks the explicit M_aL>0 check", boolResidual[signFlipIntegrity["MALPositive"]]];
  kBad = k; kBad[[1, 1]] = 2 kBad[[1, 1]];
  expectFail["K_aa diagonal corruption breaks det(K) identity", mkCitationIntegrity[m, kBad]["DetKResidual"]];
  kBad = k; kBad[[2, 2]] = 2 kBad[[2, 2]];
  expectFail["K_LL diagonal corruption breaks det(K) identity", mkCitationIntegrity[m, kBad]["DetKResidual"]];
  kBad = k; kBad[[1, 2]] = 2 kBad[[1, 2]]; kBad[[2, 1]] = 2 kBad[[2, 1]];
  expectFail["K_aL magnitude corruption breaks det(K) identity", mkCitationIntegrity[m, kBad]["DetKResidual"]];
  kBad = k; kBad[[1, 2]] = -kBad[[1, 2]]; kBad[[2, 1]] = -kBad[[2, 1]];
  kSignFlipIntegrity = mkCitationIntegrity[m, kBad];
  expectFail["K_aL sign flip breaks explicit K_aL<0 check", boolResidual[kSignFlipIntegrity["KALNegative"]]]
];

runLegacyHessian[data_] := Module[{h, e, residuals, mkSymbols, legacyNames, eNoSigmaA, hNoSigmaA, eNoChi, hNoChi},
  subheading["Legacy E_geom and own-built H_legacy"];
  h = data["HLegacy"];
  e = data["EGeom"];
  Print["  E_geom = ", fmt[e]];
  Print["  H_legacy = ", fmt[h]];
  residuals = hessianEntryResiduals[h];
  Scan[Function[key, expectZero["H_legacy entry/det " <> key, residuals[key]]], Keys[residuals]];
  expectBool["legacy_symmetric", h === Transpose[h]];
  expectBool["legacy_offdiag_negative certificate: H_aL=-chi*kappa with chi,kappa>0", clean[h[[1, 2]] + chi kappa] === 0];
  expectBool["legacy_det_positive certificate uses positive-term form", data["Gate"]["LegacyDetPositive"]];
  legacyNames = {"kappa", "chi", "sigmaA", "sigmaL"};
  mkSymbols = Sort[DeleteDuplicates[ToString[#, InputForm] & /@ Cases[{Values[data["MK"]["MEntries"]], Values[data["MK"]["KEntries"]]}, s_Symbol /; Context[s] === "Global`", Infinity]]];
  Print["  LEGACY-CONSTANT boundary: {kappa,chi,sigmaA,sigmaL} enter only here; 013 M/K free symbols = ", mkSymbols];
  expectBool["013 M/K free-symbol firewall excludes 015 legacy constants", Intersection[legacyNames, mkSymbols] === {}];
  eNoSigmaA = 1/2 kappa (deltaL - chi deltaA)^2 + 1/2 sigmaL deltaL^2;
  hNoSigmaA = D[eNoSigmaA, {{deltaA, deltaL}, 2}];
  expectFail["legacy-Hessian tooth: dropping sigmaA support term breaks H_aa", hessianEntryResiduals[hNoSigmaA]["aa"]];
  eNoChi = 1/2 kappa deltaL^2 + 1/2 sigmaA deltaA^2 + 1/2 sigmaL deltaL^2;
  hNoChi = D[eNoChi, {{deltaA, deltaL}, 2}];
  expectFail["legacy-Hessian tooth: dropping chi cross coupling breaks H_aL", hessianEntryResiduals[hNoChi]["aL"]]
];

runStructureRecovery[data_] := Module[
  {gate, probes, m, k, h, certs, entryCerts, legacyNames, mkSymbols, mProbe, kProbe, mBad, kBad, hBad, badGate},
  subheading["EARNED-STRUCTURE-RECOVERY"];
  gate = data["Gate"];
  probes = data["Probes"];
  m = data["MK"]["M"];
  k = data["MK"]["K"];
  h = data["HLegacy"];
  certs = positivityCertificateIdentities[];
  entryCerts = entryCertificateResiduals[m, k];
  Print["  Structure gate compares cited 013 M/K against own-built H_legacy: M pos-def; K symmetric; K_aL<0; rank and zero-pattern match."];
  Print["  Positivity method: exact identities f1'=2*sinh(B)^2, f2=(sinh(B)-B)(sinh(B)+B), M_aL via h=B-tanh(B); no ask/is_positive and no M[i,j]-same-form positivity tooth."];
  Scan[Function[key, expectZero["positivity certificate identity " <> key, certs[key]]], Keys[certs]];
  Scan[Function[key, expectZero["entry-to-certificate residual " <> key, entryCerts[key]]], Keys[entryCerts]];
  expectBool["M_posdef via exact Sylvester certificates", gate["MPosdef"]];
  expectBool["K_symmetric", gate["KSymmetric"]];
  expectBool["K_offdiag_negative via cited K_aL sign certificate", gate["KOffdiagNegative"]];
  expectBool["K_rank==rank(H_legacy)==2", gate["KRankMatchesLegacy"]];
  expectBool["zero_pattern(K)==zero_pattern(H_legacy)", gate["KZeroPatternMatchesLegacy"]];
  expectBool["K_structure_ok", gate["KStructureOk"]];
  legacyNames = {"kappa", "chi", "sigmaA", "sigmaL"};
  mkSymbols = Sort[DeleteDuplicates[ToString[#, InputForm] & /@ Cases[{m, k}, s_Symbol /; Context[s] === "Global`", Infinity]]];
  expectBool["structure_from_computed_matrices: cited M/K free_symbols exclude legacy {kappa,chi,sigmaA,sigmaL}", Intersection[legacyNames, mkSymbols] === {}];
  expectNonzero["recovery is structural not full entrywise fit: M_aa != H_legacy_aa", clean[m[[1, 1]] - h[[1, 1]]]];
  Print["  K_zero_pattern = ", gate["KZeroPattern"], "; legacy_zero_pattern = ", gate["LegacyZeroPattern"]];
  Print["  structure_from_computed_matrices=True; full_matrix_fit=False."];
  mProbe = probes["NonPosdefMProbe"];
  kProbe = probes["SignFlippedKProbe"];
  expectBool["probe M_aa -> -M_aa flips M_posdef false", ! TrueQ[mProbe["MPosdef"]]];
  expectBool["probe K_aL -> -K_aL flips K_offdiag_negative false", ! TrueQ[kProbe["KOffdiagNegative"]]];
  expectBool["probe K_aL -> -K_aL flips K_structure_ok false", ! TrueQ[kProbe["KStructureOk"]]];
  mBad = m; mBad[[1, 1]] = -mBad[[1, 1]];
  expectFail["structure entry tooth: corrupt M_aa sign fails M_posdef", boolResidual[structureGate[mBad, k, h]["MPosdef"]]];
  kBad = k; kBad[[1, 2]] = -kBad[[1, 2]]; kBad[[2, 1]] = -kBad[[2, 1]];
  expectFail["structure entry tooth: corrupt K_aL sign fails K_structure_ok", boolResidual[structureGate[m, kBad, h]["KStructureOk"]]];
  hBad = h; hBad[[1, 2]] = 0; hBad[[2, 1]] = 0;
  badGate = structureGate[m, k, hBad];
  expectFail["structure entry tooth: H_legacy offdiag mutation breaks zero-pattern match", boolResidual[badGate["KZeroPatternMatchesLegacy"]]]
];

runHFForce[data_] := Module[
  {hf, expectedFA, expectedFL, typedTwiceFlag, corruptedDistA, corruptReduces},
  subheading["HF-TWO-ROUTE-GENUINENESS"];
  hf = data["HF"];
  expectedFA = 4 Pi Vp0 rhoStar (Cosh[L0 beta] - 1)/(beta ellC Sinh[L0 beta]);
  expectedFL = 4 Pi Vp0 rhoStar rAL (Cosh[L0 beta] - 1)/(beta ellC Sinh[L0 beta]);
  Print["  source_density = ", fmt[data["SourceDensity"]], "; measure is action 4*pi*int dw, not mu_eta-weighted."];
  Print["  Distributed unsimplified F_a raw tree:"];
  Print["    ", rawString[hf["ForceDistAUns"]]];
  Print["  Legacy/HF unsimplified F_a raw tree:"];
  Print["    ", rawString[hf["ForceLegacyAUns"]]];
  Print["  Distributed simplified: F_a=", fmt[hf["ForceDistA"]], "; F_L=", fmt[hf["ForceDistL"]]];
  Print["  Legacy simplified: F_a=", fmt[hf["ForceLegacyA"]], "; F_L=", fmt[hf["ForceLegacyL"]]];
  expectZero["F_dist_a closed form", hf["ForceDistA"] - expectedFA];
  expectZero["F_legacy_a closed form", hf["ForceLegacyA"] - expectedFA];
  expectZero["F_dist_L closed form", hf["ForceDistL"] - expectedFL];
  expectZero["F_legacy_L closed form", hf["ForceLegacyL"] - expectedFL];
  expectZero["HF route difference F_a", hf["ForceDistA"] - hf["ForceLegacyA"]];
  expectZero["HF route difference F_L", hf["ForceDistL"] - hf["ForceLegacyL"]];
  expectBool["hf_force_reduces=True", hf["HFForceReduces"]];
  expectBool["unsimplified_routes_identical computed False", ! TrueQ[hf["UnsimplifiedRoutesIdentical"]]];
  typedTwiceFlag = rawString[hf["ForceDistAUns"]] === rawString[hf["ForceDistAUns"]];
  Print["  tooth 1a typed-twice setup: unsimplified_routes_identical=True scaffolding is de-counted; verdict tooth follows."];
  expectZero["tooth 1a typed-twice Route 2 trips BREATHING_FAIL_HF_FORCE", verdictResidual[computeVerdict[data["Gate"], hf["HFForceReduces"], typedTwiceFlag], BREATHINGFAILHFFORCE]];
  corruptedDistA = clean[hf["ForceDistA"] ellC/Vp0];
  corruptReduces = TrueQ[clean[corruptedDistA - hf["ForceLegacyA"]] === 0 && clean[hf["ForceDistL"] - hf["ForceLegacyL"]] === 0];
  expectFail["tooth 2 corrupt one HF route breaks hf_force_reduces", boolResidual[corruptReduces]];
  expectZero["tooth 2 corrupt one route trips BREATHING_FAIL_HF_FORCE", verdictResidual[computeVerdict[data["Gate"], corruptReduces, hf["UnsimplifiedRoutesIdentical"]], BREATHINGFAILHFFORCE]]
];

runStaticDynamicLimit[data_] := Module[
  {m, k, q, qdd, f, eomResidual, staticResidual},
  subheading["static-dynamic-limit"];
  Print["  Qddot -> 0 in the same M_AB Qddot + K_AB Q = F_A system gives K_AB Q = F_A."];
  Print["  After the legacy-pattern comparison this is the static dE_geom/dQ=0 limit; no separate static solve is fabricated."];
  m = data["MK"]["M"];
  k = data["MK"]["K"];
  q = {qa, qL};
  qdd = {qddA, qddL};
  f = {forceA, forceL};
  eomResidual = m.qdd + k.q - f;
  staticResidual = clean[(eomResidual /. {qddA -> 0, qddL -> 0}) - (k.q - f)];
  expectZero["static_limit_consistent: Qddot->0 EOM residual equals K.Q-F", staticResidual.staticResidual];
  Print["  separate_static_solve=False (scope note only; no counted literal tooth)."]
];

runHFProfileGuard[data_] := Module[{guards, hf},
  subheading["HF-PROFILE-GUARD / 014<->015 boundary"];
  guards = data["Guards"];
  hf = data["HF"];
  Print["  Degenerate alpha_a=0: wrong_zero_F_a = ", fmt[guards["WrongZeroFA"]]];
  Print["  Constant alpha_a=1: wrong_const_F_a = ", fmt[guards["WrongConstFA"]]];
  Print["  Frozen correct legacy force F_a = ", fmt[hf["ForceLegacyA"]]];
  Print["  014<->015 boundary: constant_one overlap_passes=True in 014, but hf_mismatch=True here; HF is the profile guard the overlap could not supply."];
  expectZero["wrong_zero_F_a is native zero integral", guards["WrongZeroFA"]];
  expectZero["wrong_const_F_a is 4*pi*L0*Vp0*rhoStar/ellC", guards["WrongConstFA"] - 4 Pi L0 Vp0 rhoStar/ellC];
  expectBool["degenerate hf_mismatch=True", guards["DegenerateHFMismatch"]];
  expectBool["constant hf_mismatch=True", guards["ConstantHFMismatch"]];
  expectBool["constant profile passed 014 overlap but fails 015 HF", guards["ConstantOverlapPasses014"] && guards["ConstantHFMismatch"]];
  expectNonzero["HF guard reads a real object: degenerate wrong_zero_F_a genuinely differs from frozen legacy force", clean[guards["WrongZeroFA"] - hf["ForceLegacyA"]]];
  expectNonzero["HF guard reads a real object: constant wrong_const_F_a genuinely differs from frozen legacy force", clean[guards["WrongConstFA"] - hf["ForceLegacyA"]]]
];

runVerdictAndComposition[data_] := (
  subheading["015 scoped landing and joint composition"];
  Print["  015 scoped verdict = ", data["Verdict"]];
  expectZero["015 component lands at BREATHING_CALIBRATED", verdictResidual[data["Verdict"], BREATHINGCALIBRATED]];
  expectBool["STRUCTURE rung passes", data["Gate"]["MPosdef"] && data["Gate"]["KSymmetric"] && data["Gate"]["KStructureOk"]];
  expectBool["HF rung passes with anti-x-x guard", data["HF"]["HFForceReduces"] && ! TrueQ[data["HF"]["UnsimplifiedRoutesIdentical"]]];
  expectBool["HF guard rejects both wrong profiles", data["Guards"]["DegenerateHFMismatch"] && data["Guards"]["ConstantHFMismatch"]];
  Print["  BREATHING_CALIBRATED (JOINT, 3-stage) -- COMPLETE"];
  Print["    = (013: harmonic profiles + M/K operator-projection + (a,L) EOM LHS)[done, cited]"];
  Print["    AND (014: truncation consistency)[done, cited]"];
  Print["    AND (015: legacy-Hessian structure recovery + Hellmann-Feynman force + static-dynamic limit)[computed here]"];
  Print["  RHS FILLED: the EOM M_AB Qddot + K_AB Q = F_A^(HF) now has F_A^(HF)=4*pi*int (rhoStar*Vp0/ellC)*alpha_A dw."];
  Print["  CALIBRATED <= {muEta,Tw,beta} (013) [+ the HF drive scale rhoStar*Vp0/ellC -- register question, not decided here]."];
  Print["  CAVEATS: HF is the profile guard the overlap could not supply; {kappa,chi,sigmaA,sigmaL} are the legacy-Hessian pattern basis; BdG k^4 deferred."]
);

printProvenance[data_] := Module[{liveNames},
  subheading["Provenance and scope"];
  Print["  CONSUMED-from-013/014: stage013 profiles alpha_a,alpha_L + operator-projected M_AB/K_AB + frozen packet are cited with dual-site integrity; stage014 window beta_L0 in [0.1,3.0] / K_eta/Tw<=~2.6 is cited; no sibling recomputation."];
  Print["  no-c_S: 015 is speed-free (static structure + static force); matter-sector c_s/BdG k^4 is deferred under k*xi << 1 (phonon_limit_caveat)."];
  Print["  LEGACY-CONSTANT boundary: {kappa,chi,sigmaA,sigmaL} enter HERE as E_geom parameters and are absent from 013 M/K."];
  Print["  EARNED-STRUCTURE-RECOVERY: own-built H_legacy plus genuine structural match of cited M/K; recovered, not re-postulated; full_matrix_fit=False."];
  Print["  HF-TWO-ROUTE-GENUINENESS: distributed projection vs Hellmann-Feynman parametric derivative; hf_force_reduces=True and unsimplified_routes_identical=False."];
  Print["  HF-PROFILE-GUARD: degenerate and constant wrong profiles have hf_mismatch=True; constant_one passed 014 overlap but fails 015 HF."];
  Print["  static-dynamic-limit: Qddot->0 in the same dynamical closure gives K_AB Q = F_A and the static dE_geom/dQ=0 limit; no separate static solve."];
  Print["  3-way-split-COMPLETING: 015 is component 3/3 and completes BREATHING_CALIBRATED; the RHS F_A^(HF) deferred by 013 is filled here."];
  Print["  dropped-bookkeeping: scratch-WL exports, sympy* comparisons, expression_digest, engine_agreement, YAML/JSON plumbing, and report writers are stripped."];
  Print["  downstream consumers: stages 022/023 consume the completed ell=0 (a,L) breathing closure."];
  Print["  register note: 015 introduces Vp0 and makes ell_c live through rhoStar*Vp0/ellC; whether that force scale is a counted CALIB knob is intentionally left to registration."];
  liveNames = Sort[DeleteDuplicates[ToString[#, InputForm] & /@ Cases[{data["Profiles"]["AlphaA"], data["Profiles"]["AlphaL"], Values[data["MK"]["MEntries"]], Values[data["MK"]["KEntries"]], data["SourceDensity"], data["HF"]["ForceDistA"], data["HF"]["ForceLegacyA"]}, s_Symbol /; Context[s] === "Global`", Infinity]]];
  expectBool["no c_S/cS live symbol appears in 015 symbolic content", FreeQ[liveNames, "cS"] && FreeQ[liveNames, "c_S"]]
];

printVerdictLabels[] := (
  Print[""];
  Print["Verdict labels:"];
  Print["  ledger earned-label (NOT a source verdict token): BREATHING_STRUCTURE_HF_FORCE_EARNED  (legacy-Hessian structure recovery: own-built H_legacy=hessian(E_geom); the CITED 013 M_AB/K_AB match its structural signature (M pos-def by Sylvester; K symmetric; K_aL<0 <=> H_legacy off-diagonal -chi*kappa<0 for chi>0; rank(K)=rank(H_legacy)=2; zero-pattern match) => the (a,L) operator-projected closure is RECOVERED not re-postulated; the Hellmann-Feynman force F_A^(HF) derived by TWO genuinely-different routes -- distributed projection F_dist=4pi int (rho*Vp0/ellC) alpha_A dw vs Hellmann-Feynman parametric F_legacy=-4pi int rho* dV_conf/dq_A|0 dw -- that AGREE (hf_force_reduces=True) with unsimplified_routes_identical=False (the anti-x-x: raw trees differ); the HF REJECTS the wrong profiles (hf_mismatch=True for alpha_a=0 and alpha_a=1) -- the constant profile passed 014's overlap but fails 015's HF)"];
  Print["  source top-line verdict: BREATHING_CALIBRATED  (JOINT 3-stage; 015 carries the legacy-structure + HF-force component 3/3 and COMPLETES the joint)"];
  Print["  joint composition (COMPLETE): BREATHING_CALIBRATED = (013: harmonic profiles + M/K operator-projection + (a,L) EOM LHS)[done] AND (014: truncation consistency)[done] AND (015: legacy-Hessian structure recovery + Hellmann-Feynman force + static-dynamic limit)[here] ; the EOM RHS F_A^(HF) 013 deferred is filled here"];
  Print["  earned (structure): H_legacy own-built; the cited 013 M/K match its structural signature (genuine pos-def/symmetric/negative-offdiagonal/rank/zero-pattern comparisons, NOT X==X on cited literals); the two probes (M_aa->-M_aa, K_aL->-K_aL) flip the gate"];
  Print["  earned (HF force): both routes genuinely different constructions agreeing after simplification; unsimplified_routes_identical computed False (anti-x-x); corrupting one route breaks hf_force_reduces"];
  Print["  calibrated (values): the wall constants {mu_eta,Tw,beta} are stage-013 calibration inputs; the HF drive scale rho*Vp0/ellC is the section-6 register question (likely a new CALIB force-scale -- confirm at registration) -> BREATHING_CALIBRATED not ..._PASS"];
  Print["  carried caveats (015's own, honest): the HF is the profile guard the overlap could not supply (constant_one PASSES 014's overlap, FAILS 015's HF); {kappa,chi,sigmaA,sigmaL} are the legacy-Hessian pattern basis (a re-parameterization of the calibrated set; structural recovery not full numeric fit, full_matrix_fit=False); BdG k^4 matter-sector deferred (k*xi<<1)"];
  Print["  consumed (cited from stage013/014, dual-site integrity): collective profiles alpha_a,alpha_L; operator-projected M_AB/K_AB; frozen packet L0=37/20, beta=1, beta*L0=37/20; 014 truncation certificate (window beta_L0 in [0.1,3.0], K_eta/Tw<=~2.6); c_S NOT consumed (matter-sector deferred)"];
  Print["  new symbols (015): legacy {kappa,chi,sigmaA,sigmaL} (E_geom parameterization; ABSENT from 013's M/K); HF drive Vp0 (confinement slope) + ell_c LIVE (delta_V_conf != 0 here; INERT at 011) + rho* (frozen density, consumed)"]
);

runBaselineImmutability[data_] := (
  subheading["Baseline immutability after copy mutations"];
  expectBool["baseline M_posdef remains true", data["Gate"]["MPosdef"]];
  expectBool["baseline K_structure_ok remains true", data["Gate"]["KStructureOk"]];
  expectBool["baseline hf_force_reduces remains true", data["HF"]["HFForceReduces"]];
  expectBool["baseline unsimplified_routes_identical remains false", ! TrueQ[data["HF"]["UnsimplifiedRoutesIdentical"]]];
  expectBool["baseline M/K citation integrity remains true", data["MKIntegrity"]["Ok"]];
  expectZero["baseline clean 015 verdict remains BREATHING_CALIBRATED", verdictResidual[data["Verdict"], BREATHINGCALIBRATED]]
);

Module[{ok, data},
  heading["ledger_stage015_breathing_structure_hf_force Mathematica audit"];
  ok = Catch[
    data = buildBaseline[];
    assertExact["baseline", data];
    runAritySelfCheck[data];
    runConsumedInputs[data];
    runLegacyHessian[data];
    runStructureRecovery[data];
    runHFForce[data];
    runStaticDynamicLimit[data];
    runHFProfileGuard[data];
    runVerdictAndComposition[data];
    printProvenance[data];
    printVerdictLabels[];
    runBaselineImmutability[data];
    True,
    "ledgerStage015Failure",
    Function[{msg, tag}, Print["FAIL: ", msg]; False]
  ];

  Print[""];
  Print["PASS tally: ", passCount, "; FAIL tally: ", failCount, "; TOTAL = ", passCount, " + ", failCount, " = ", passCount + failCount];
  If[TrueQ[ok],
    Print["OVERALL PASS: Mathematica verified ledger_stage015 breathing structure recovery + Hellmann-Feynman force exactly"];
    Exit[0],
    Print["OVERALL FAIL: Mathematica stage015 audit did not close"];
    Exit[1]
  ]
]
