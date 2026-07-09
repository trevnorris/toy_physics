(* Ledger stage019 Mathematica audit: squared-denominator prefactor algebra.

   Standalone, print-only, no arguments, no file I/O. This keeps the native
   Wolfram Series/Coefficient/FullSimplify route for the pathA_33 II-G4b
   abstract-algebra slice, including the plain-N/D discriminator and the
   dynamically rerun 019-local 3g failure. Adjacent products are provenance.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
failureMessage = "";

PREFACTORALGEBRAEARNED = "PREFACTOR_ALGEBRA_EARNED";
FAILPREFACTORALGEBRA = "FAIL_PREFACTOR_ALGEBRA";
NOFAIL = "NO_FAIL";

$Assumptions = Element[{omega, D0, D2, D4, N0, N2, N4}, Reals] &&
  D0 != 0 && D2 != 0 && D4 != 0 && N0 != 0 && N2 != 0 && N4 != 0;

fail[msg_] := Throw[msg, "ledgerStage019Failure"];

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
    fail[name]
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
    fail[name]
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
    fail[name]
  ]
];

verdictResidual[actual_, expected_] := If[actual === expected, 0, 1];

serW[expr_, order_] := Collect[
  Normal[Series[expr, {omega, 0, order}]],
  omega,
  FullSimplify
];

extractCoefficients[obj_] := Module[{derivedSeries},
  derivedSeries = serW[obj, 6];
  <|
    "Series" -> derivedSeries,
    "P0" -> clean[Coefficient[derivedSeries, omega, 0]],
    "P2" -> clean[Coefficient[derivedSeries, omega, 2]],
    "P4" -> clean[Coefficient[derivedSeries, omega, 4]]
  |>
];

scopedVerdict[prefactorOk_, ndSelfCheckOk_] := If[
  TrueQ[prefactorOk && ndSelfCheckOk],
  PREFACTORALGEBRAEARNED,
  FAILPREFACTORALGEBRA
];

dynamicAblation[baselinePrefactorOk_, mutatedPrefactorOk_, ndSelfCheckOk_] := Module[
  {withMutation, withoutMutation},
  withMutation = scopedVerdict[mutatedPrefactorOk, ndSelfCheckOk];
  withoutMutation = scopedVerdict[baselinePrefactorOk, ndSelfCheckOk];
  <|
    "RerunGateLogic" -> True,
    "WithMutation" -> withMutation,
    "WithoutMutation" -> withoutMutation,
    "ExpectedFail" -> FAILPREFACTORALGEBRA,
    "FailSuppressed" -> TrueQ[
      withMutation === FAILPREFACTORALGEBRA &&
      withoutMutation =!= FAILPREFACTORALGEBRA
    ]
  |>
];

Nomega = N0 + N2 omega^2 + N4 omega^4;
Dcons = D0 + D2 omega^2 + D4 omega^4;
prefObj = D0 Nomega/Dcons^2;
plainObj = Nomega/Dcons;

prefData = extractCoefficients[prefObj];
plainData = extractCoefficients[plainObj];
prefSeries = prefData["Series"];
plainSeries = plainData["Series"];
p0 = prefData["P0"];
p2 = prefData["P2"];
p4 = prefData["P4"];
plainP0 = plainData["P0"];
plainP2 = plainData["P2"];
plainP4 = plainData["P4"];

expectedP0 = N0/D0;
expectedP2 = (D0 N2 - 2 D2 N0)/D0^2;
expectedP4 = (D0^2 N4 - 2 D0 (D2 N2 + D4 N0) + 3 D2^2 N0)/D0^3;
expected = <|"P0" -> expectedP0, "P2" -> expectedP2, "P4" -> expectedP4|>;
coefficients = <|"P0" -> p0, "P2" -> p2, "P4" -> p4|>;
plainCoefficients = <|"P0" -> plainP0, "P2" -> plainP2, "P4" -> plainP4|>;

resP0 = clean[p0 - expectedP0];
resP2 = clean[p2 - expectedP2];
resP4 = clean[p4 - expectedP4];
residuals = <|"P0" -> resP0, "P2" -> resP2, "P4" -> resP4|>;
p0Match = TrueQ[resP0 === 0];
p2Match = TrueQ[resP2 === 0];
p4Match = TrueQ[resP4 === 0];
prefactorOk = TrueQ[p0Match && p2Match && p4Match];

differencePlainMinusCorrectP2 = clean[plainP2 - expectedP2];
plainEqualsCorrectP2 = TrueQ[clean[plainP2 == expectedP2]];
ndSelfCheckOk = ! plainEqualsCorrectP2;
correctEqualsCorrectP2 = TrueQ[clean[p2 == expectedP2]];

mutatedMatches = AssociationMap[
  TrueQ[clean[plainCoefficients[#] - expected[#]] === 0] &,
  {"P0", "P2", "P4"}
];
mutatedPrefactorOk = TrueQ[And @@ Values[mutatedMatches]];
probeVerdict = If[! plainEqualsCorrectP2, FAILPREFACTORALGEBRA, NOFAIL];
correctObjectVerdict = If[correctEqualsCorrectP2, NOFAIL, FAILPREFACTORALGEBRA];
prefactorAblation = dynamicAblation[prefactorOk, mutatedPrefactorOk, ndSelfCheckOk];

gateBooleans = <|
  "prefactor_match" -> prefactorOk,
  "N_over_D_self_check" -> ndSelfCheckOk
|>;
verdict = scopedVerdict[prefactorOk, ndSelfCheckOk];

sampleRules = {
  D0 -> 19, D2 -> 23, D4 -> 29,
  N0 -> 11, N2 -> 13, N4 -> 17
};
sampleValues = <|
  "P0" -> clean[p0 /. sampleRules],
  "P2" -> clean[p2 /. sampleRules],
  "P4" -> clean[p4 /. sampleRules],
  "plainP2" -> clean[plainP2 /. sampleRules],
  "plainMinusCorrectP2" -> clean[differencePlainMinusCorrectP2 /. sampleRules]
|>;

symbolNames[expr_] := DeleteDuplicates[
  Cases[expr, s_Symbol /; Context[s] === "Global`" :> SymbolName[s], Infinity]
];

runPrefactorDerivation[] := (
  subheading["Squared-denominator object and series-extracted coefficients"];
  Print["  N(omega) = ", fmt[Nomega]];
  Print["  Dcons(omega) = ", fmt[Dcons]];
  Print["  correctObj = ", fmt[prefObj]];
  Print["  correctSeries = ", fmt[prefSeries]];
  Print["  P0/P2/P4 below are read by Coefficient[correctSeries,omega,n]; targets are independent typed references."];
  Scan[
    Function[name,
      Print["  ", name, "(series-extracted) = ", fmt[coefficients[name]]];
      expectZero["series-extracted " <> name <> " matches independent expected " <> name, residuals[name]]
    ],
    {"P0", "P2", "P4"}
  ];
  Print["  exact sample controls (D0,D2,D4,N0,N2,N4 only) = ", sampleRules];
  Print["  exact sample values (transcript only, non-tooth) = ", sampleValues]
);

runNdSelfCheck[] := (
  subheading["Plain N/D factor-of-2 discriminator"];
  Print["  plainObj = ", fmt[plainObj]];
  Print["  plainSeries = ", fmt[plainSeries]];
  Print["  plainP2(series-extracted) = ", fmt[plainP2]];
  Print["  differencePlainMinusCorrectP2(computed) = ", fmt[differencePlainMinusCorrectP2]];
  Print["  plainEqualsCorrectP2 = ", plainEqualsCorrectP2];
  expectZero[
    "plain N/D P2 is series-extracted with the single D2*N0 term",
    plainP2 - (D0 N2 - D2 N0)/D0^2
  ];
  expectZero[
    "computed plain-minus-correct P2 is D2*N0/D0^2",
    differencePlainMinusCorrectP2 - D2 N0/D0^2
  ];
  expectBool[
    "plainEqualsCorrectP2 is False (squared-denominator factor-of-2 detected)",
    ! plainEqualsCorrectP2
  ]
);

runProbe[] := (
  subheading["Probe 3g wrong-prefactor object and dynamic 019-local self-ablation"];
  Print["  3g mutated object = ", fmt[plainObj]];
  Print["  3g plainEqualsCorrectP2 = ", plainEqualsCorrectP2];
  Print["  3g wrong-object verdict = ", probeVerdict];
  Print["  3g correct-object verdict = ", correctObjectVerdict];
  Print["  3g dynamic self-ablation = ", prefactorAblation];
  expectZero[
    "3g plain N/D reaches FAIL_PREFACTOR_ALGEBRA",
    verdictResidual[probeVerdict, FAILPREFACTORALGEBRA]
  ];
  expectZero[
    "3g correct squared-denominator object is NO_FAIL",
    verdictResidual[correctObjectVerdict, NOFAIL]
  ];
  expectZero[
    "3g dynamic rerun with plain-N/D mutation fails locally",
    verdictResidual[prefactorAblation["WithMutation"], FAILPREFACTORALGEBRA]
  ];
  expectZero[
    "3g dynamic rerun without mutation returns the earned local verdict",
    verdictResidual[prefactorAblation["WithoutMutation"], PREFACTORALGEBRAEARNED]
  ];
  expectBool["3g self-ablation suppresses the failure", prefactorAblation["FailSuppressed"]]
);

runPerToothAblations[] := Module[{mutated, correctControl, correctControlEquals},
  subheading["Per-tooth derivation ablations on independent series copies"];
  Scan[
    Function[pair,
      mutated = extractCoefficients[prefObj + omega^pair[[2]]];
      expectFail[
        pair[[1]] <> " derivation copy changes its own series-extracted coefficient",
        mutated[pair[[1]]] - expected[pair[[1]]]
      ]
    ],
    {{"P0", 0}, {"P2", 2}, {"P4", 4}}
  ];
  correctControl = extractCoefficients[prefObj];
  correctControlEquals = TrueQ[clean[correctControl["P2"] == expectedP2]];
  expectBool[
    "N/D discriminator comparison flips True when its plain control is replaced by the correct object",
    correctControlEquals
  ];
  expectZero[
    "baseline remains immutable after derivation-copy ablations",
    Total[Values[residuals]]
  ]
];

runAritySelfCheck[] := Module[{aritySeries, arityCoefficients, arityAblation},
  subheading["Wolfram arity self-check and unevaluated-leakage scan"];
  aritySeries = serW[1 + omega^2, 4];
  arityCoefficients = extractCoefficients[prefObj];
  arityAblation = dynamicAblation[True, False, True];
  expectZero[
    "arity serW[expr,order] returns the requested omega^2 coefficient",
    Coefficient[aritySeries, omega, 2] - 1
  ];
  expectBool[
    "arity extractCoefficients[obj] returns P0/P2/P4 keys",
    And @@ (KeyExistsQ[arityCoefficients, #] & /@ {"P0", "P2", "P4"})
  ];
  expectBool[
    "arity scopedVerdict[prefactorOk,ndSelfCheckOk] returns earned label",
    scopedVerdict[True, True] === PREFACTORALGEBRAEARNED
  ];
  expectBool[
    "arity dynamicAblation[baseline,mutated,selfCheck] returns changed verdicts",
    arityAblation["WithMutation"] =!= arityAblation["WithoutMutation"]
  ];
  expectBool[
    "no unevaluated SeriesData remains in emitted series",
    FreeQ[{prefSeries, plainSeries}, _SeriesData]
  ];
  expectBool[
    "no unevaluated Series/Coefficient/FullSimplify call remains",
    FreeQ[{prefSeries, plainSeries, p0, p2, p4, plainP2}, Series | Coefficient | FullSimplify]
  ]
];

runScopeAndProvenance[] := Module[{names, allowedNames},
  subheading["019 scope, provenance-only consumption, and PARTIAL landing"];
  Print["  019 gate booleans = ", gateBooleans];
  Print["  019 scoped verdict = ", verdict];
  expectZero[
    "019 scoped verdict lands the earned partial component",
    verdictResidual[verdict, PREFACTORALGEBRAEARNED]
  ];
  names = symbolNames[{
    Nomega, Dcons, prefObj, plainObj, prefSeries, plainSeries,
    Values[coefficients], Values[plainCoefficients]
  }];
  allowedNames = {"omega", "D0", "D2", "D4", "N0", "N2", "N4"};
  Print["  live symbolic names in the earned slice = ", Sort[names]];
  expectBool[
    "earned algebra has exactly omega and D0..N4 as live symbols",
    Sort[names] === Sort[allowedNames]
  ];
  Print["  QUAD_CALIBRATED (2/4) -- squared-denominator prefactor algebra EARNED (PARTIAL)"];
  Print["    = P0=N0/D0, P2=(D0*N2-2*D2*N0)/D0^2, P4=(D0^2*N4-2*D0*(D2*N2+D4*N0)+3*D2^2*N0)/D0^3."];
  Print["    SERIES-EXTRACTED from P(omega)=D0*N/Dcons^2; the factor-of-2 is the squared-denominator signature."];
  Print["    AND plain N/D reaches FAIL_PREFACTOR_ALGEBRA while the correct object is NO_FAIL."];
  Print["  REMAINING: fingerprint = 018 (DONE); magnitude partition plus calibrated label = 020; dimension closure = 021."];
  Print["  CAVEATS: D-lane and N-moment numerical branch scalars remain symbolic/Gate-6 deferred; the magnitude is calibrated in 020."];
  Print["  CONSUMES (PROVENANCE only): 017 D-lanes + deferred port N-moments + 018 fingerprint context; no guard/dual-site."];
  Print["  EXPORTS: the series-extracted prefactor algebra and N/D self-check to 020/021 and 027; no file artifact is written."];
  Print["  PORT-MOMENTS PROVENANCE EXPORT -- LABELED NON-CHECK: concrete N-moments are deferred Gate-6 branch data, asserted against nothing."];
  Print["  reduction certificate: FROZEN-INPUT symbolic D_n lanes and outgoing port-moment inputs/N_n; COMPUTED P0/P2/P4 plus N/D self-check; DEFERRED numerical branch scalars."];
  Print["  dropped-bookkeeping: scratch-YAML agreement handoff and report/feed writers are absent; agreement is transcript-level."];
  Print["  register note: likely zero new counted knobs; registration decides the structural edge."]
];

printVerdictLabels[] := (
  subheading["Verdict labels"];
  Print["  ledger earned-label (NOT a source verdict token): PREFACTOR_ALGEBRA_EARNED  (the squared-denominator prefactor P(omega)=D0*N(omega)/Dcons(omega)^2, N=N0+N2*omega^2+N4*omega^4, Dcons=D0+D2*omega^2+D4*omega^4, series-expands to the coefficients P0=N0/D0, P2=(D0*N2-2*D2*N0)/D0^2, P4=(D0^2*N4-2*D0*(D2*N2+D4*N0)+3*D2^2*N0)/D0^3 -- SERIES-EXTRACTED, NOT typed; the -2*D2*N0 factor-of-2 is the squared-denominator signature (plain N/D gives only -D2*N0, provably wrong -> FAIL_PREFACTOR_ALGEBRA))"];
  Print["  source top-line verdict: QUAD_CALIBRATED  (JOINT 4-stage; 019 carries the squared-denominator prefactor-algebra component 2/4 and lands the token as a PARTIAL)"];
  Print["  joint composition (019 = the SECOND leg; 018 DONE, 020/021 remaining): QUAD_CALIBRATED = (018: outgoing DtN Hankel fingerprint + chi_Q sign)[EARNED, PARTIAL, DONE] AND (019: squared-denominator prefactor algebra P(omega)=D0*N/Dcons^2)[EARNED here, PARTIAL] AND (020: 54/5=2*27/5 provenance partition + the CALIBRATED label, G=GENUINE_BLOCKED) AND (021: mu_hat0-free [P0_phys]=1 dim closure)"];
  Print["  earned (the prefactor algebra): P0=N0/D0, P2=(D0*N2-2*D2*N0)/D0^2, P4=(D0^2*N4-2*D0*(D2*N2+D4*N0)+3*D2^2*N0)/D0^3 SERIES-EXTRACTED from P(omega)=D0*N/Dcons^2 (NOT typed); the N/D self-check (plain N/D gives -D2*N0 vs the correct -2*D2*N0 -> FAIL_PREFACTOR_ALGEBRA; the correct object does NOT fire)"];
  Print["  earned (able-to-fail): 3g_wrong_prefactor_object (plain N/D -> missing factor of 2 -> FAIL_PREFACTOR_ALGEBRA; correct object NO_FAIL), with a DYNAMIC self-ablation (re-run, NOT a constant); the prefactor coefficients are SERIES-EXTRACTED off the actual series (not a hardcode-and-compare, the section-4 firewall)"];
  Print["  calibrated / deferred (NOT 019): the fingerprint (018, DONE); the 54/5=2*27/5 magnitude + G (020, external_bridge_input, G=GENUINE_BLOCKED); the mu_hat0-free dim closure (021); the numerical (D_n,N_n) port scalars (Gate-6 sim-deferred, report :49)"];
  Print["  consumed (PROVENANCE only -- NO guard/dual-site; abstract port-agnostic algebra): 017's l=2 port kernel (D-lanes D0/D2/D4 + {Btilde,Ztilde}, carried as abstract symbols) + build_port_moments' concrete N-moments (deferred Gate-6 branch data, asserted-against-nothing) + 018's fingerprint context; NO c_s/a/G/mu_hat0"];
  Print["  exports (the prefactor algebra): P0=N0/D0 (-> 021's [P0_phys]=1 dim closure) + P2/P4 + the N/D self-check (-> 020's 54/5=2*27/5 partition context) => 020/021 + 027 (pathA_43 closure)"];
  Print["  new symbols first-appearing (019): none new-counted (the D0/D2/D4 are 017's exported D-lanes + the N0/N2/N4 are build_port_moments' deferred Gate-6 port N-moments, NOT new knobs); no units symbol (units-free abstract algebra); no counted knobs expected (an EARNED prefactor-algebra slice, like 018)"]
);

runAll[] := (
  heading["ledger_stage019_prefactor_algebra_mathematica_audit"];
  Print["Target stem confirmed: ledger_stage019_prefactor_algebra"];
  Print["Engine: native Wolfram exact Series/Coefficient/FullSimplify derivation; no floats/tolerances; zero file I/O."];
  runPrefactorDerivation[];
  runNdSelfCheck[];
  runProbe[];
  runPerToothAblations[];
  runAritySelfCheck[];
  runScopeAndProvenance[];
  printVerdictLabels[];
  0
);

result = Catch[
  runAll[],
  "ledgerStage019Failure",
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
