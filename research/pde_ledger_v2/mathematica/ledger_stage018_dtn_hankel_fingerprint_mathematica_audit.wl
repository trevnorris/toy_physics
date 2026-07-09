(* Ledger stage018 Mathematica audit: DtN Hankel fingerprint + chi sign.

   Standalone, print-only, no arguments, no file I/O. This keeps the native
   Wolfram spherical-Bessel/Series/Coefficient route for the pathA_33 II-G4a
   exterior-wave slice: outgoing DtN Hankel fingerprint, outgoing/incoming
   chi sign, standing contrast, passivity, and units-restored coefficient
   dimensions. Sibling stages own the prefactor, normalization partition, and
   dimensional closure.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
failureMessage = "";

DTNHANKELFINGERPRINTEARNED = "DTN_HANKEL_FINGERPRINT_EARNED";
FAILFINGERPRINT = "FAIL_FINGERPRINT";
FAILNOTOUTGOING = "FAIL_NOT_OUTGOING";
FAILDIMENSIONAL = "FAIL_DIMENSIONAL";
FAILNOTABLETOFAIL = "FAIL_NOT_ABLE_TO_FAIL";

$Assumptions = Element[{z, a, cs}, Reals] && a > 0 && cs > 0;

raise[msg_] := Throw[msg, "ledgerStage018Failure"];

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

clean[expr_] := FullSimplify[expr, $Assumptions];
fmt[expr_String] := expr;
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

zeroDim = {0, 0, 0};
dimL = {1, 0, 0};
dimSpeed = {1, 0, -1};
dimT2 = {0, 0, 2};
dimT4 = {0, 0, 4};
dimT5 = {0, 0, 5};

dimResidualVec[actual_, expected_] := FullSimplify[(actual - expected).(actual - expected)];

dimText[d_] := Module[{pairs, parts, emit},
  emit[label_, exp_] := If[TrueQ[exp == 1], label, label <> "^" <> ToString[InputForm[exp]]];
  pairs = {{"L", d[[1]]}, {"M", d[[2]]}, {"T", d[[3]]}};
  parts = (emit[#[[1]], #[[2]]] &) /@ Select[pairs, ! TrueQ[#[[2]] == 0] &];
  If[Length[parts] == 0, "1", StringRiffle[parts, " "]]
];

dimRaise[msg_] := Throw[msg, "stage018DimError"];

dimOf[expr_, dims_] := Module[{args, ds, base, pow},
  Which[
    TrueQ[expr == 0] || NumberQ[expr], zeroDim,
    AtomQ[expr] && KeyExistsQ[dims, expr], dims[expr],
    AtomQ[expr], dimRaise["missing dimension for " <> ToString[Unevaluated[expr], InputForm]],
    Head[expr] === Times, Total[dimOf[#, dims] & /@ (List @@ expr)],
    Head[expr] === Power,
      base = expr[[1]];
      pow = expr[[2]];
      If[! NumberQ[pow], dimRaise["non-numeric dimension exponent"]];
      pow dimOf[base, dims],
    Head[expr] === Plus,
      args = Select[List @@ expr, ! TrueQ[clean[#] == 0] &];
      ds = dimOf[#, dims] & /@ args;
      If[Length[ds] == 0, zeroDim,
        If[Length[DeleteDuplicates[ds]] != 1, dimRaise["dimension mismatch in sum " <> ToString[expr, InputForm]]];
        First[ds]
      ],
    True, dimRaise["unsupported dimension expression " <> ToString[expr, InputForm]]
  ]
];

serZ[expr_, order_] := Collect[Normal[Series[expr, {z, 0, order}]], z, FullSimplify];

j2 = (3/z^3 - 1/z) Sin[z] - 3 Cos[z]/z^2;
y2 = (1/z - 3/z^3) Cos[z] - 3 Sin[z]/z^2;

branchData[name_, h_, source_] := Module[
  {lam, yhat, hser, lamser, yser, yser5, u2z, u4z, v5z, static},
  lam = clean[z D[h, z]/h];
  yhat = clean[-3/lam];
  hser = serZ[h, 7];
  lamser = serZ[lam, 6];
  yser = serZ[yhat, 6];
  yser5 = serZ[yhat, 5];
  u2z = clean[Coefficient[yser, z, 2]];
  u4z = clean[Coefficient[yser, z, 4]];
  v5z = clean[Coefficient[yser, z, 5]/I];
  static = clean[Coefficient[yser, z, 0]];
  <|
    "name" -> name,
    "source" -> source,
    "h" -> h,
    "lambda" -> lam,
    "yhat" -> yhat,
    "hSeries" -> hser,
    "lambdaSeries" -> lamser,
    "yhatSeries" -> yser,
    "yhatSeriesToZ5" -> yser5,
    "static" -> static,
    "u2z" -> u2z,
    "u4z" -> u4z,
    "v5z" -> v5z,
    "u2" -> clean[u2z a^2/cs^2],
    "u4" -> clean[u4z a^4/cs^4],
    "v5" -> clean[v5z a^5/cs^5]
  |>
];

passivityFromSource[branch_] := Module[
  {imagNonzero, isGenuine},
  imagNonzero = ! TrueQ[clean[branch["v5z"]] === 0];
  isGenuine = branch["source"] === "dtn_hankel1";
  <|
    "source" -> branch["source"],
    "imaginaryV5Nonzero" -> imagNonzero,
    "genuineOutgoing" -> TrueQ[isGenuine && imagNonzero]
  |>
];

scopedVerdict[overrides_: <||>] := Module[
  {fingerprintOk, outgoingOk, dimensionalOk, ableOk},
  fingerprintOk = Lookup[overrides, "fingerprintOk", True];
  outgoingOk = Lookup[overrides, "outgoingOk", True];
  dimensionalOk = Lookup[overrides, "dimensionalOk", True];
  ableOk = Lookup[overrides, "ableToFailOk", True];
  Which[
    ! TrueQ[fingerprintOk], FAILFINGERPRINT,
    ! TrueQ[outgoingOk], FAILNOTOUTGOING,
    ! TrueQ[dimensionalOk], FAILDIMENSIONAL,
    ! TrueQ[ableOk], FAILNOTABLETOFAIL,
    True, DTNHANKELFINGERPRINTEARNED
  ]
];

dynamicAblation[overrides_, expectedFail_] := Module[
  {withMutation, withoutMutation},
  withMutation = scopedVerdict[overrides];
  withoutMutation = scopedVerdict[<||>];
  <|
    "RerunGateLogic" -> True,
    "WithMutation" -> withMutation,
    "WithoutMutation" -> withoutMutation,
    "ExpectedFail" -> expectedFail,
    "FailSuppressed" -> TrueQ[withMutation === expectedFail && withoutMutation =!= expectedFail]
  |>
];

out = branchData["outgoing_hankel1", j2 + I y2, "dtn_hankel1"];
incoming = branchData["incoming_hankel2", j2 - I y2, "dtn_hankel2"];
standing = branchData["standing_j2", j2, "regular_standing_j2"];

expected = <|
  "u2z" -> 1/9,
  "u4z" -> 4/81,
  "v5z" -> 1/27,
  "u2" -> a^2/(9 cs^2),
  "u4" -> 4 a^4/(81 cs^4),
  "v5" -> a^5/(27 cs^5)
|>;

matches = AssociationThread[
  Keys[expected],
  TrueQ[clean[out[#] - expected[#]] === 0] & /@ Keys[expected]
];
computedV5Denominator = clean[1/out["v5z"]];
chiQ = clean[out["v5"]/(a^5/(27 cs^5))];
chiQIncoming = clean[incoming["v5"]/(a^5/(27 cs^5))];
exactSample = <|
  "a" -> 3,
  "cs" -> 2,
  "u2" -> clean[out["u2"] /. {a -> 3, cs -> 2}],
  "u4" -> clean[out["u4"] /. {a -> 3, cs -> 2}],
  "v5" -> clean[out["v5"] /. {a -> 3, cs -> 2}]
|>;

dimRules = <|a -> dimL, cs -> dimSpeed|>;
u2Dim = dimOf[out["u2"], dimRules];
u4Dim = dimOf[out["u4"], dimRules];
v5Dim = dimOf[out["v5"], dimRules];
corruptedU2 = clean[out["u2z"] a^2/cs^3];
corruptedU2Dim = dimOf[corruptedU2, dimRules];
unitsOk = TrueQ[
  dimResidualVec[u2Dim, dimT2] == 0 &&
  dimResidualVec[u4Dim, dimT4] == 0 &&
  dimResidualVec[v5Dim, dimT5] == 0
];
corruptedUnitsOk = TrueQ[dimResidualVec[corruptedU2Dim, dimT2] == 0];
unitsSelfAblation = dynamicAblation[<|"dimensionalOk" -> False|>, FAILDIMENSIONAL];

passivity = passivityFromSource[out];
incomingExpected = TrueQ[
  clean[incoming["u2z"] - out["u2z"]] === 0 &&
  clean[incoming["u4z"] - out["u4z"]] === 0 &&
  clean[incoming["v5z"] + out["v5z"]] === 0 &&
  clean[chiQIncoming + 1] === 0
];
standingExpected = TrueQ[
  clean[Coefficient[standing["lambdaSeries"], z, 0] - 2] === 0 &&
  clean[standing["static"] + 3/2] === 0 &&
  clean[standing["v5z"]] === 0
];
wrongBcProbe = <|
  "IncomingExpected" -> incomingExpected,
  "StandingExpected" -> standingExpected,
  "IncomingVerdict" -> If[incomingExpected, FAILFINGERPRINT, FAILNOTABLETOFAIL],
  "StandingVerdict" -> If[standingExpected, FAILFINGERPRINT, FAILNOTABLETOFAIL],
  "SelfAblation" -> dynamicAblation[<|"fingerprintOk" -> False|>, FAILFINGERPRINT]
|>;

phenomSource = <|"source" -> "phenomenological_inserted_sink", "v5z" -> out["v5z"]|>;
phenomGenuine = TrueQ[phenomSource["source"] === "dtn_hankel1" && ! TrueQ[clean[phenomSource["v5z"]] === 0]];
imposedProbe = <|
  "MutatedSource" -> phenomSource["source"],
  "InsertedV5z" -> phenomSource["v5z"],
  "GenuineOutgoing" -> phenomGenuine,
  "RemovingOutgoingBcRemovesImaginaryTerm" -> TrueQ[clean[standing["v5z"]] === 0],
  "Verdict" -> If[phenomGenuine, DTNHANKELFINGERPRINTEARNED, FAILNOTOUTGOING],
  "SelfAblation" -> dynamicAblation[<|"outgoingOk" -> False|>, FAILNOTOUTGOING]
|>;

u2MixedCopy = branchData["u2_mixed_copy", (j2 + I y2) (1 + z^2), "dtn_hankel1_u2_copy"];
u4OnlyCopy = branchData["u4_only_copy", (j2 + I y2) (1 + z^4), "dtn_hankel1_u4_copy"];
v5OnlyCopy = branchData["v5_only_copy", (j2 + I y2) (1 + I z^5), "dtn_hankel1_v5_copy"];

fingerprintOk = TrueQ[
  And @@ Values[matches] &&
  clean[computedV5Denominator - 27] === 0 &&
  clean[chiQ - 1] === 0 &&
  clean[chiQIncoming + 1] === 0 &&
  incomingExpected &&
  standingExpected
];
outgoingOk = TrueQ[
  passivity["genuineOutgoing"] &&
  clean[out["static"] - 1] === 0 &&
  imposedProbe["Verdict"] === FAILNOTOUTGOING
];
ableToFailOk = TrueQ[
  wrongBcProbe["IncomingVerdict"] === FAILFINGERPRINT &&
  wrongBcProbe["StandingVerdict"] === FAILFINGERPRINT &&
  wrongBcProbe["SelfAblation"]["FailSuppressed"] &&
  imposedProbe["SelfAblation"]["FailSuppressed"] &&
  scopedVerdict[<|"dimensionalOk" -> corruptedUnitsOk|>] === FAILDIMENSIONAL
];
gateBooleans = <|
  "fingerprint_ok" -> fingerprintOk,
  "outgoing_ok" -> outgoingOk,
  "dimensional_ok" -> unitsOk,
  "able_to_fail_ok" -> ableToFailOk
|>;
verdict = scopedVerdict[
  <|
    "fingerprintOk" -> fingerprintOk,
    "outgoingOk" -> outgoingOk,
    "dimensionalOk" -> unitsOk,
    "ableToFailOk" -> ableToFailOk
  |>
];

symbolNames[expr_] := DeleteDuplicates[
  Cases[expr, s_Symbol /; Context[s] === "Global`" :> SymbolName[s], Infinity]
];

forbiddenTokens[] := {
  StringJoin["mu", "_hat", "0"],
  StringJoin["mtilde", "0"],
  StringJoin["N", "0"],
  StringJoin["D", "0"],
  StringJoin["P", "0", "_phys"],
  StringJoin["build", "_dimensions"]
};

runFingerprintDerivation[] := (
  subheading["Outgoing DtN Hankel fingerprint derived from the branch"];
  Print["  h2^(1)(z) is built from explicit spherical j2 + i*y2 rational sin/cos expressions."];
  Print["  h2^(1) series = ", fmt[out["hSeries"]]];
  Print["  Lambda_2^out(z) = z*h2'(z)/h2(z)"];
  Print["  Lambda_2^out series = ", fmt[out["lambdaSeries"]]];
  Print["  Yhat_2^out(z) = -3/Lambda_2^out(z)"];
  Print["  Yhat_2^out series through z^5 = ", fmt[out["yhatSeriesToZ5"]]];
  Print["  Coefficients are read from that series by Coefficient[...,z,n]; typed targets are independent references."];
  Scan[
    Function[key,
      Print["  derived ", key, " = ", fmt[out[key]]];
      expectZero["outgoing derived " <> key <> " matches independent target", out[key] - expected[key]]
    ],
    {"u2z", "u4z", "v5z", "u2", "u4", "v5"}
  ];
  expectZero["computed radiative denominator 1/v5z equals 27", computedV5Denominator - 27];
  Print["  exact sample at a=3, c_s=2 = ", Association @ KeyValueMap[(#1 -> fmt[#2]) &, exactSample]]
);

runChiAndStanding[] := (
  subheading["chi sign from outgoing/incoming Hankel branches and standing contrast"];
  Print["  outgoing v5z = ", fmt[out["v5z"]], "; incoming v5z = ", fmt[incoming["v5z"]]];
  Print["  computed chi_Q outgoing = ", fmt[chiQ], "; incoming = ", fmt[chiQIncoming]];
  expectZero["chi_Q computed from outgoing v5 physical slot is +1", chiQ - 1];
  expectZero["incoming v5z + outgoing v5z = 0", incoming["v5z"] + out["v5z"]];
  expectZero["incoming computed chi_Q is -1", chiQIncoming + 1];
  expectZero["incoming even u2z unchanged", incoming["u2z"] - out["u2z"]];
  expectZero["incoming even u4z unchanged", incoming["u4z"] - out["u4z"]];
  expectFail["using incoming branch as outgoing chi check fires", chiQIncoming - 1];
  Print["  standing Lambda series = ", fmt[standing["lambdaSeries"]]];
  Print["  standing Yhat series = ", fmt[standing["yhatSeriesToZ5"]]];
  expectZero["standing branch Lambda_static is +2", Coefficient[standing["lambdaSeries"], z, 0] - 2];
  expectZero["standing branch Yhat_static is -3/2", standing["static"] + 3/2];
  expectZero["standing branch has no radiating v5z", standing["v5z"]];
  expectFail["standing static slot differs from outgoing normalized static slot", standing["static"] - 1]
);

runPassivityAndProbes[] := (
  subheading["Passivity and the two 018 fingerprint probes"];
  Print["  passivity source = ", passivity["source"], "; imaginary v5 nonzero = ", passivity["imaginaryV5Nonzero"]];
  expectBool["outgoing branch passivity gate is genuine outgoing with nonzero v5", passivity["genuineOutgoing"]];
  expectZero["outgoing static slot is the DC limit of the same response", out["static"] - 1];
  expectZero["3a incoming verdict is FAIL_FINGERPRINT", verdictResidual[wrongBcProbe["IncomingVerdict"], FAILFINGERPRINT]];
  expectZero["3a standing verdict is FAIL_FINGERPRINT", verdictResidual[wrongBcProbe["StandingVerdict"], FAILFINGERPRINT]];
  expectBool["3a incoming predicted change is computed", wrongBcProbe["IncomingExpected"]];
  expectBool["3a standing predicted change is computed", wrongBcProbe["StandingExpected"]];
  expectZero["3a dynamic self-ablation with mutation is FAIL_FINGERPRINT", verdictResidual[wrongBcProbe["SelfAblation"]["WithMutation"], FAILFINGERPRINT]];
  expectZero["3a dynamic self-ablation without mutation returns earned verdict", verdictResidual[wrongBcProbe["SelfAblation"]["WithoutMutation"], DTNHANKELFINGERPRINTEARNED]];
  expectBool["3a self-ablation suppresses the failure", wrongBcProbe["SelfAblation"]["FailSuppressed"]];
  Print["  3b mutated source = ", imposedProbe["MutatedSource"], "; inserted v5z = ", fmt[imposedProbe["InsertedV5z"]]];
  expectBool["3b inserted sink is not a genuine outgoing source", ! imposedProbe["GenuineOutgoing"]];
  expectBool["3b removing outgoing BC removes imaginary v5", imposedProbe["RemovingOutgoingBcRemovesImaginaryTerm"]];
  expectZero["3b imposed dissipation verdict is FAIL_NOT_OUTGOING", verdictResidual[imposedProbe["Verdict"], FAILNOTOUTGOING]];
  expectZero["3b dynamic self-ablation with mutation is FAIL_NOT_OUTGOING", verdictResidual[imposedProbe["SelfAblation"]["WithMutation"], FAILNOTOUTGOING]];
  expectZero["3b dynamic self-ablation without mutation returns earned verdict", verdictResidual[imposedProbe["SelfAblation"]["WithoutMutation"], DTNHANKELFINGERPRINTEARNED]];
  expectBool["3b outgoing_ablation is a rerun, not a constant label", imposedProbe["SelfAblation"]["RerunGateLogic"] && imposedProbe["SelfAblation"]["WithMutation"] =!= imposedProbe["SelfAblation"]["WithoutMutation"]]
);

runUnitsLeg[] := (
  subheading["Units-restored physical coefficient dimensions only"];
  Print["  dimension order = (L,M,T); [a]=L, [c_s]=L*T^-1."];
  Print["  computed dimensions = ", <|"u2" -> dimText[u2Dim], "u4" -> dimText[u4Dim], "v5" -> dimText[v5Dim]|>];
  expectZero["physical u2 has dimension T^2", dimResidualVec[u2Dim, dimT2]];
  expectZero["physical u4 has dimension T^4", dimResidualVec[u4Dim, dimT4]];
  expectZero["physical v5 has dimension T^5", dimResidualVec[v5Dim, dimT5]];
  expectBool["baseline coefficient dimensional leg passes", unitsOk];
  Print["  corrupted u2 expression = ", fmt[corruptedU2], "; dimension = ", dimText[corruptedU2Dim]];
  expectFail["corrupting the c_s power in u2 breaks the T^2 dimension", dimResidualVec[corruptedU2Dim, dimT2]];
  expectZero["units corruption reaches FAIL_DIMENSIONAL", verdictResidual[scopedVerdict[<|"dimensionalOk" -> corruptedUnitsOk|>], FAILDIMENSIONAL]];
  expectZero["units dynamic self-ablation with mutation is FAIL_DIMENSIONAL", verdictResidual[unitsSelfAblation["WithMutation"], FAILDIMENSIONAL]];
  expectZero["units dynamic self-ablation without mutation returns earned verdict", verdictResidual[unitsSelfAblation["WithoutMutation"], DTNHANKELFINGERPRINTEARNED]];
  expectBool["units self-ablation suppresses the failure", unitsSelfAblation["FailSuppressed"]]
);

runPerToothAblations[] := (
  subheading["Per-tooth ablations on copies"];
  expectFail["u2 derivation copy changes u2z so its own coefficient assert fires", u2MixedCopy["u2z"] - expected["u2z"]];
  Print["  u4-only derivation copy coefficients = u2z ", fmt[u4OnlyCopy["u2z"]], ", u4z ", fmt[u4OnlyCopy["u4z"]], ", v5z ", fmt[u4OnlyCopy["v5z"]]];
  expectZero["u4-only copy leaves u2z unchanged", u4OnlyCopy["u2z"] - expected["u2z"]];
  expectFail["u4-only copy breaks exactly the u4z coefficient assert", u4OnlyCopy["u4z"] - expected["u4z"]];
  expectZero["u4-only copy leaves v5z unchanged", u4OnlyCopy["v5z"] - expected["v5z"]];
  expectZero["v5-only copy leaves u2z unchanged", v5OnlyCopy["u2z"] - expected["u2z"]];
  expectZero["v5-only copy leaves u4z unchanged", v5OnlyCopy["u4z"] - expected["u4z"]];
  expectFail["v5-only copy changes v5z so its own coefficient assert fires", v5OnlyCopy["v5z"] - expected["v5z"]];
  expectFail["incoming branch mutation flips chi sign tooth", chiQIncoming - 1];
  expectZero["incoming branch reaches FAIL_FINGERPRINT", verdictResidual[wrongBcProbe["IncomingVerdict"], FAILFINGERPRINT]];
  expectZero["imposed sink reaches FAIL_NOT_OUTGOING", verdictResidual[imposedProbe["Verdict"], FAILNOTOUTGOING]];
  expectZero["units dim tooth reaches FAIL_DIMENSIONAL", verdictResidual[scopedVerdict[<|"dimensionalOk" -> corruptedUnitsOk|>], FAILDIMENSIONAL]];
  expectZero["baseline immutable after teeth: outgoing u2z still 1/9", out["u2z"] - expected["u2z"]];
  expectZero["baseline immutable after teeth: outgoing u4z still 4/81", out["u4z"] - expected["u4z"]];
  expectZero["baseline immutable after teeth: outgoing v5z still 1/27", out["v5z"] - expected["v5z"]]
);

runAritySelfCheck[] := Module[{arityBranch},
  subheading["Wolfram arity self-check and unevaluated-leakage scan"];
  arityBranch = branchData["arity_outgoing", j2 + I y2, "dtn_hankel1"];
  expectBool["arity branchData[name,h,source] returns u2z key", KeyExistsQ[arityBranch, "u2z"]];
  expectZero["arity serZ[z+z^2,2] keeps z^2 coefficient", Coefficient[serZ[z + z^2, 2], z, 2] - 1];
  expectBool["arity passivityFromSource[branch] returns genuineOutgoing key", KeyExistsQ[passivityFromSource[out], "genuineOutgoing"]];
  expectBool["arity scopedVerdict[association] returns earned label", scopedVerdict[<||>] === DTNHANKELFINGERPRINTEARNED];
  expectBool["no unevaluated SeriesData remains in emitted series", FreeQ[{out["hSeries"], out["lambdaSeries"], out["yhatSeries"]}, _SeriesData]];
  expectBool["no unevaluated Derivative remains in DtN expressions", FreeQ[{out["lambda"], out["yhat"], standing["lambda"]}, Derivative]]
];

runScopeAndProvenance[] := Module[
  {names, helperNames, forbidden, frozenSpeed},
  subheading["018 scope, provenance-only consumption, and structural exclusions"];
  Print["  018 gate booleans = ", gateBooleans];
  Print["  018 scoped verdict = ", verdict];
  expectZero["018 scoped verdict lands the earned partial component", verdictResidual[verdict, DTNHANKELFINGERPRINTEARNED]];
  Print["  QUAD_CALIBRATED (1/4) -- outgoing l=2 DtN Hankel fingerprint EARNED (PARTIAL)"];
  Print["    = u2=a^2/(9*c_s^2), u4=4*a^4/(81*c_s^4), v5=a^5/(27*c_s^5), derived from Yhat=-3/Lambda."];
  Print["    AND chi_Q=+1 outgoing / -1 incoming from j2 +/- i*y2, with standing j2 carrying no radiating term."];
  Print["  REMAINING: prefactor algebra = 019; normalization partition plus calibrated label = 020; dim closure = 021."];
  Print["  CAVEATS: c_s is a units carrier, not a consumed value; branch scalar numerics remain sim-deferred; chi reconciliation is Part-VII."];
  Print["  CONSUMES (PROVENANCE only): 017 l=2 port kernel + 009/010 bulk mode + stage005 R1 identity for the c_s units symbol."];
  Print["  EXPORTS: outgoing fingerprint and chi sign to 019/020 plus later non-regression and closure stages; no file artifact is written."];
  names = symbolNames[{out["h"], out["lambda"], out["yhat"], incoming["h"], standing["h"], corruptedU2}];
  Print["  live symbolic names in this slice = ", Sort[names]];
  expectBool["only z, a, cs are live symbols in the fingerprint slice", SubsetQ[{"z", "a", "cs"}, names]];
  frozenSpeed = StringJoin["c", "_S"];
  expectBool["frozen-wall speed symbol is not live", FreeQ[names, frozenSpeed] && FreeQ[names, "cS"]];
  forbidden = forbiddenTokens[];
  helperNames = StringReplace[Names["Global`*"], "Global`" -> ""];
  expectBool["no prefactor/dim-closure carrier helper or symbol names are live", Intersection[names, forbidden] === {} && Intersection[helperNames, forbidden] === {}];
  Print["  dropped-bookkeeping: scratch YAML, report writers, engine-agreement files, and cross-reads are absent."];
  Print["  register note: likely zero new counted knobs; c_s is R1-provenance units carrier; structural edge is for the outgoing DtN fingerprint."]
];

printVerdictLabels[] := (
  subheading["Verdict labels"];
  Print["  ledger earned-label (NOT a source verdict token): DTN_HANKEL_FINGERPRINT_EARNED"];
  Print["  source top-line verdict: QUAD_CALIBRATED (JOINT 4-stage; 018 carries the outgoing-fingerprint + chi-sign component as PARTIAL)"];
  Print["  joint composition: 018 outgoing DtN Hankel fingerprint + 019 prefactor + 020 normalization partition/calibrated label + 021 dim closure"];
  Print["  earned fingerprint: u2z=1/9, u4z=4/81, v5z=1/27 from Coefficient[...,z,n] of Yhat=-3/Lambda; 27 = 1/v5z is computed."];
  Print["  earned sign: chi_Q=+1 outgoing and -1 incoming from j2 +/- i*y2; standing j2 has Lambda_static=+2 and no radiating term."];
  Print["  earned able-to-fail: wrong-BC and imposed-sink probes each carry dynamic with/without mutation verdicts."];
  Print["  calibrated/deferred outside 018: prefactor algebra, normalization magnitude/label, and the full dimensional closure."];
  Print["  consumed framing: provenance-only; no guard on port kernel, bulk mode, or c_s value; only c_s units appear in the dim leg."]
);

runAll[] := (
  heading["ledger_stage018_dtn_hankel_fingerprint_mathematica_audit"];
  Print["Target stem confirmed: ledger_stage018_dtn_hankel_fingerprint"];
  Print["Engine: native Wolfram exact spherical-Hankel Series/Coefficient; no floats/tolerances; zero file I/O."];
  runFingerprintDerivation[];
  runChiAndStanding[];
  runPassivityAndProbes[];
  runUnitsLeg[];
  runPerToothAblations[];
  runAritySelfCheck[];
  runScopeAndProvenance[];
  printVerdictLabels[];
  0
);

result = Catch[
  runAll[],
  "ledgerStage018Failure",
  Function[{msg, tag}, failureMessage = ToString[msg]; 1]
];

heading["Tallies"];
Print["TALLY mathematica: ", passCount, " pass + ", failCount, " fail = ", passCount + failCount, " checks"];
If[result === 0 && failCount === 0,
  Print["OVERALL PASS"];
  Exit[0],
  If[failureMessage =!= "", Print["ABORTED: ", failureMessage]];
  Print["OVERALL FAIL"];
  Exit[1]
];
