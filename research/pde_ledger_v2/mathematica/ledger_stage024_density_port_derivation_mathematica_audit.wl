(* Ledger stage024 Wolfram audit: density-native l=2 two-port derivation.

   Standalone, print-only, no arguments, and no file I/O.  The load-bearing
   route stays native: eliminate q2 from the wall row, substitute it into the
   bulk row, Solve for Phi2, and only then build N0_den from that response.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
failureMessage = "";

localVerdict = "DENSITY_TWO_PORT_DERIVED";
jointPartial = "DENSITY_PORT_HOSTED (1/4, DERIVATION — the density-native two-port N0_den; 025 = vector-freedom taint, 026 = continuity lineage, 027 = port checks + closure)";

$Assumptions =
  Element[{a, cS, rhoEff, I25, xiQ, etaQ, etaPhi, varpiQ2, varpiPhi2, lambdaC}, Reals] &&
  a > 0 && cS > 0 && rhoEff > 0 && varpiQ2 > 0 && varpiPhi2 > 0 &&
  I25 != 0 && xiQ != 0 && etaQ != 0 && etaPhi != 0;

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

clean[expr_] := Factor[Cancel[Together[FullSimplify[expr, $Assumptions]]]];
fmt[expr_String] := expr;
externalSymbolNames = <|
  "cS" -> "c_s", "rhoEff" -> "rho_eff", "xiQ" -> "Xi_Q",
  "etaQ" -> "eta_q", "etaPhi" -> "eta_phi",
  "varpiQ2" -> "varpi_q2", "varpiPhi2" -> "varpi_Phi2",
  "lambdaC" -> "lambda_c", "phi2" -> "Phi2",
  "OmegaU" -> "Omega_U", "OmegaW" -> "Omega_W",
  "Rmix" -> "R_mix", "gU" -> "g_U", "gW" -> "g_W"
|>;
fmt[expr_] := StringReplace[
  ToString[InputForm[clean[expr]]],
  Normal[externalSymbolNames]
];

fail[] := Throw[failureMessage, "ledgerStage024Failure"];

expectZero[name_, residual_] := Module[{c},
  c = clean[residual];
  If[TrueQ[c === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    failureMessage = name <> ": residual = " <> fmt[c];
    Print["FAIL  ", failureMessage];
    fail[]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

makeN0[responseArgument_] := clean[responseArgument^2];

densityHostNames = {
  "q2", "Phi2", "c_s", "a", "rho_eff", "I25", "Xi_Q", "eta_q",
  "eta_phi", "varpi_q2", "varpi_Phi2", "lambda_c"
};
forbiddenVectorNames = {
  "A_w", "F_muw", "Jw", "U", "W", "Omega_U", "Omega_W",
  "R_mix", "g_U", "g_W"
};

globalSymbolNames[expr_] := Sort[DeleteDuplicates[
  Lookup[externalSymbolNames, SymbolName[#], SymbolName[#]] & /@
  Cases[Unevaluated[expr], symbol_Symbol /; Context[symbol] === "Global`", Infinity]
]];

guardNonsingular[delta_, context_String] := expectBool[
  "B Delta!=0 nonsingular guard before Solve (" <> context <> ")",
  ! TrueQ[clean[delta] === 0]
];

(* Keeping lambdaOp separate provides the singular-operator ablation point. *)
lambdaOp = lambdaC;

densityExpression[gBaseInput_, context_String] := Module[
  {gBase, gqDen, gphiDen, delta, qRule, phiEq, phiResponse, pden, n0},
  gBase = clean[gBaseInput];
  gqDen = clean[gBase etaQ];
  gphiDen = clean[gBase etaPhi];
  delta = clean[varpiQ2 varpiPhi2 - lambdaOp^2];

  (* Tooth B is deliberately before qRule and Solve. *)
  guardNonsingular[delta, context];

  (* Genuine Green-DtN route: eliminate the wall coordinate first. *)
  qRule = q2 -> clean[(lambdaOp phi2 + gqDen)/varpiQ2];
  phiEq = clean[varpiPhi2 phi2 - lambdaOp q2 == gphiDen /. qRule];
  phiResponse = clean[phi2 /. First[Solve[phiEq, phi2]]];

  (* Independent adjugate/determinant oracle; it does not construct N0_den. *)
  pden = clean[varpiQ2 gphiDen + lambdaOp gqDen];
  n0 = makeN0[phiResponse];
  <|
    "StaticOperator" -> {{varpiQ2, -lambdaOp}, {-lambdaOp, varpiPhi2}},
    "Source" -> {gqDen, gphiDen},
    "QRule" -> qRule,
    "PhiEquation" -> phiEq,
    "Phi2Response" -> phiResponse,
    "DeltaDen" -> delta,
    "PDen" -> pden,
    "N0Den" -> n0,
    "GBase" -> gBase
  |>
];

definitionArity[function_Symbol] := Module[{definitions, lhs},
  definitions = DownValues[function];
  If[Length[definitions] =!= 1, Return[-1]];
  lhs = Extract[definitions, {1, 1}, HoldComplete];
  Extract[
    lhs,
    {1, 1},
    Function[call, Length[Unevaluated[call]], HoldAllComplete]
  ]
];

heldCallArity[held_HoldComplete] := Extract[
  held,
  {1},
  Function[call, Length[Unevaluated[call]], HoldAllComplete]
];

runAudit[] := Module[
  {gBase, baselineCall, makeN0Call, guardCall, baseline, deltaProbe,
   zeroCall, zeroControl, liveNames, transcriptObjects, expectedKeys},
  gBase = clean[Sqrt[rhoEff] cS^2 I25 xiQ/a^(7/2)];

  subheading["E. Wolfram arity scan"];
  baselineCall = HoldComplete[densityExpression[gBase, "baseline"]];
  makeN0Call = HoldComplete[makeN0[phi2]];
  guardCall = HoldComplete[guardNonsingular[varpiQ2 varpiPhi2 - lambdaC^2, "arity probe"]];
  expectBool[
    "E definition/call arity scan matches densityExpression, makeN0, and guardNonsingular",
    definitionArity[densityExpression] === 2 && heldCallArity[baselineCall] === 2 &&
    definitionArity[makeN0] === 1 && heldCallArity[makeN0Call] === 1 &&
    definitionArity[guardNonsingular] === 2 && heldCallArity[guardCall] === 2
  ];
  baseline = ReleaseHold[baselineCall];
  expectedKeys = {
    "StaticOperator", "Source", "QRule", "PhiEquation", "Phi2Response",
    "DeltaDen", "PDen", "N0Den", "GBase"
  };
  expectBool[
    "E densityExpression Module call evaluates to the complete result association",
    AssociationQ[baseline] && Sort[Keys[baseline]] === Sort[expectedKeys]
  ];

  subheading["A. elimination factorization and response-to-N0 dataflow"];
  expectZero[
    "A eliminated response equals independent P_den/Delta oracle",
    baseline["Phi2Response"] - baseline["PDen"]/baseline["DeltaDen"]
  ];
  expectBool[
    "A make_N0 runtime dataflow probe reads its response argument",
    ! TrueQ[
      clean[
        makeN0[baseline["Phi2Response"] + deltaProbe] -
        makeN0[baseline["Phi2Response"]]
      ] === 0
    ]
  ];

  subheading["C. coupling-vanishes boundary"];
  expectBool[
    "C baseline N0_den is nonzero for symbolic nonzero coupling",
    ! TrueQ[clean[baseline["N0Den"]] === 0]
  ];
  zeroCall = HoldComplete[densityExpression[0, "g_base=0 control"]];
  zeroControl = ReleaseHold[zeroCall];
  expectZero[
    "C zero-control N0_den|g_base=0 is computed zero",
    zeroControl["N0Den"]
  ];

  subheading["D. density-only host set"];
  liveNames = globalSymbolNames[baseline["N0Den"]];
  expectBool[
    "D N0_den host-set membership is density-only and vector-free",
    SubsetQ[densityHostNames, liveNames] &&
    Intersection[liveNames, forbiddenVectorNames] === {}
  ];

  subheading["E. Wolfram unevaluated-leakage transcript scan"];
  transcriptObjects = {
    baseline, zeroControl, baseline["StaticOperator"], baseline["Phi2Response"],
    baseline["N0Den"], liveNames
  };
  expectBool[
    "E no unevaluated authored helper or Solve leakage remains in transcript objects",
    FreeQ[
      transcriptObjects,
      densityExpression | guardNonsingular | makeN0 | Solve | FullSimplify |
      Simplify | Together | Cancel | Factor
    ]
  ];

  Append[baseline, "LiveNames" -> liveNames]
];

emit[data_Association] := (
  subheading["Derived density two-port transcript"];
  Print["port_picture: ii two-port(q2,Phi2)"];
  Print["method: Wolfram Green-DtN eliminate-q2 then Solve Phi2"];
  Print["static_operator: ", fmt[data["StaticOperator"]]];
  Print["g_base: ", fmt[data["GBase"]]];
  Print["Phi2_response: ", fmt[data["Phi2Response"]]];
  Print["Delta_den: ", fmt[data["DeltaDen"]]];
  Print["P_den: ", fmt[data["PDen"]]];
  Print["N0_den (canonical factored): ", fmt[data["N0Den"]]];
  Print["Delta!=0 guard status: PASS (checked before each Solve)"];
  Print["density-only host-set: {", StringRiffle[densityHostNames, ", "], "}"];
  Print["N0_den live symbols: {", StringRiffle[data["LiveNames"], ", "], "}"];

  subheading["physical_relations provenance"];
  Print["varpi_q2 <- pathA_32 wall: K2/M2 = (c_s/a)^2*kappa_q from the wall angular l=2 operator (stages016/017)"];
  Print["varpi_Phi2 <- pathA_29 bulk: (c_s/a)^2*(6+(m*a)^2) from the bulk Helmholtz l=2 mode (stages009/010/016)"];
  Print["lambda_c <- projected continuity/interface: (c_s/a)^2*lambda_hat_Q"];
  Print["I25: typed l=2 continuity-moment input; lineage validation is forward-referenced to stage026"];
  Print["rho_eff: effective reduced-3D mass density; its literal reduction is SIM_DEFERRED/GAP and it is not stage005 rho0"];

  subheading["Verdict labels"];
  Print["LOCAL_AUDIT_VERDICT: ", localVerdict];
  Print["JOINT_LANDING_LABEL (PARTIAL): ", jointPartial]
);

runAll[] := Module[{data},
  heading["ledger_stage024_density_port_derivation_mathematica_audit"];
  Print["Target stem confirmed: ledger_stage024_density_port_derivation"];
  Print["Engine: exact Wolfram Green-DtN elimination; zero file I/O."];
  data = runAudit[];
  emit[data];
  0
];

result = Catch[
  runAll[],
  "ledgerStage024Failure",
  Function[{msg, tag}, 1]
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
