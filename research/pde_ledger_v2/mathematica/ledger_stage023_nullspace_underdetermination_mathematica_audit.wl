(* Ledger stage023 Mathematica audit: native nullspace underdetermination.

   Standalone, print-only, exact, and zero-file-I/O.  This is a genuinely
   re-authored Wolfram route: it constructs NullSpace bases and measures the
   return-map image on those bases.  It does not use the SymPy engine's
   augmented-rank subtraction and does not rebuild stage022's fingerprints or
   Gate-4 non-regression.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
failureMessage = "";

CROSSLRESIDUALPREDICTION = "CROSS_L_RESIDUAL_PREDICTION";
FAILUNDERDETERMINED = "FAIL_UNDERDETERMINED_NOT_PREDICTIVE";
FAILDECOUPLED = "FAIL_DECOUPLED";
FAILTAUTOLOGICAL = "FAIL_TAUTOLOGICAL";
FAILDIMENSIONAL = "FAIL_DIMENSIONAL";
FAILNOCONSISTENTRETURN = "FAIL_NO_CONSISTENT_RETURN";
FAILOVERCANCEL = "FAIL_OVERCANCEL";
FAILEPSILONMISMATCH = "FAIL_EPSILON_MISMATCH";
NOFAIL = "NO_FAIL";

$Assumptions = Element[
    {
      omega, a, cs, M0, D1, R0, R1, D0, K0c, Keta, TOmega,
      Z0ret, Z1ret, OmegaU, OmegaW, Rmix, gU, gW,
      etaNull, gain0, gain1, qfree
    },
    Reals
  ] && a > 0 && cs > 0 && D0 != 0 && M0 != 0 && D1 != 0 &&
  K0c > 0 && Keta > 0 && TOmega > 0 && Z0ret > 0 && Z1ret > 0 &&
  OmegaU != 0 && OmegaW != 0 && Rmix != 0 && gU != 0 && gW != 0 &&
  etaNull > 0 && gain0 > 0 && gain1 > 0;

fail[msg_] := Throw[msg, "ledgerStage023Failure"];

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

clean[expr_] := FullSimplify[Cancel[Together[expr]], $Assumptions];
fmt[expr_String] := expr;
fmt[expr_] := ToString[InputForm[clean[expr]]];
boolText[value_] := If[TrueQ[value], "True", "False"];

assertExact[name_, expr_] := Module[{machineReals},
  machineReals = Cases[Unevaluated[expr], _Real, Infinity];
  If[machineReals =!= {},
    failCount++;
    Print["FAIL  ", name, ": machine-real atom(s) found: ", InputForm[machineReals]];
    fail[name]
  ]
];

expectZero[name_, residual_] := Module[{c0},
  assertExact[name, residual];
  c0 = clean[residual];
  assertExact[name, c0];
  If[TrueQ[c0 === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FAIL  ", name, ": residual = ", fmt[c0]];
    fail[name]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

expectFail[name_, residual_] := Module[{c0},
  assertExact[name, residual];
  c0 = clean[residual];
  assertExact[name, c0];
  If[! TrueQ[c0 === 0],
    passCount++;
    Print["PASS  ", name, " produced required FAIL (residual = ", fmt[c0], ")"],
    failCount++;
    Print["FAIL  ", name, ": required mutation/ablation did not fire"];
    fail[name]
  ]
];

selfConsistency[name_String, condition_] := If[
  TrueQ[condition],
  Print["SELF-CONSISTENCY OK (not earned; not counted)  ", name],
  Print["SELF-CONSISTENCY FAIL (not earned; not counted)  ", name];
  fail["self-consistency check failed: " <> name]
];

verdictResidual[actual_, expected_] := If[actual === expected, 0, 1];

(* Stage017 grouped-P2 port kernel retained; no stage019 prefactor is built. *)
portKernel[factor_] := Module[{ou, pport, dport, nport, p0},
  ou = factor OmegaU;
  pport = ou^2 gW + Rmix gU;
  dport = ou^2 OmegaW^2 - Rmix^2;
  nport = clean[pport^2/dport^2];
  p0 = clean[nport/D0];
  <|
    "PPort" -> pport,
    "DeltaPort" -> dport,
    "N0Port" -> nport,
    "P0Raw" -> p0,
    "P0Physical" -> clean[(cs/a)^2 p0]
  |>
];

baselinePort = portKernel[1];
P0port = baselinePort["P0Raw"];
K1dc = Keta + 2 TOmega;
rankDofs = {OmegaU, OmegaW, Rmix, gU, gW, D0, K0c, Keta, TOmega, Z0ret, Z1ret};
baseConstraints = {P0port, K0c, K1dc};

constraintJacobian[constraints_List, dofs_List] := Table[
  clean[D[constraints[[row]], dofs[[column]]]],
  {row, Length[constraints]},
  {column, Length[dofs]}
];

constructNullBasis[constraints_List, dofs_List] := NullSpace[constraintJacobian[constraints, dofs]];

returnGradients[t0_, t1_, dofs_List] := constraintJacobian[{t0, t1}, dofs];

returnImageRank[basis_List, gradients_List] := MatrixRank[basis . Transpose[gradients]];

T0dc = clean[K0c/(K0c + Z0ret)];
T1dc = clean[K1dc/(K1dc + Z1ret)];
eps0 = clean[Z0ret/K0c];
eps1 = clean[Z1ret/K1dc];

Jbase = constraintJacobian[baseConstraints, rankDofs];
basis = constructNullBasis[baseConstraints, rankDofs];
Greturn = returnGradients[T0dc, T1dc, rankDofs];
rank0 = MatrixRank[Jbase];
nativeNullity = Length[basis];
returnMovingNullity = returnImageRank[basis, Greturn];

z0Unit = UnitVector[Length[rankDofs], First[First[Position[rankDofs, Z0ret]]]];
z1Unit = UnitVector[Length[rankDofs], First[First[Position[rankDofs, Z1ret]]]];

selectorResiduals = {Z0ret - K0c, Z1ret - K1dc};
selectorConstraints = Join[baseConstraints, selectorResiduals];
Jselector = constraintJacobian[selectorConstraints, rankDofs];
selectorBasis = constructNullBasis[selectorConstraints, rankDofs];
selectorReturnImage = clean[Greturn . Transpose[selectorBasis]];
selectorRank = MatrixRank[Jselector];
selectorNativeNullity = Length[selectorBasis];
selectorReturnMovingNullity = MatrixRank[selectorReturnImage];

selectorMetadata[mode_String] := Switch[mode,
  "selector",
    <|
      "Present" -> True, "DerivedFromNamedPDE" -> False,
      "ControlOnly" -> True, "Tautological" -> False,
      "Equations" -> selectorResiduals
    |>,
  "asserted_selector",
    <|
      "Present" -> True, "DerivedFromNamedPDE" -> False,
      "ControlOnly" -> False, "Tautological" -> True,
      "Equations" -> selectorResiduals
    |>,
  _,
    <|
      "Present" -> False, "DerivedFromNamedPDE" -> False,
      "ControlOnly" -> False, "Tautological" -> False,
      "Equations" -> {}
    |>
];

pathA29Premise = <|"ZIsPremise" -> True, "BoundaryDof" -> "none"|>;

oneSelectorConstraints = Append[baseConstraints, Z0ret - K0c];
oneSelectorBasis = constructNullBasis[oneSelectorConstraints, rankDofs];
oneSelectorReturnMovingNullity = returnImageRank[oneSelectorBasis, Greturn];

injectDofs = Append[rankDofs, etaNull];
injectBasis = constructNullBasis[baseConstraints, injectDofs];
injectT0 = clean[K0c/(K0c + Z0ret + etaNull K0c)];
injectT1 = clean[K1dc/(K1dc + Z1ret + etaNull K1dc)];
injectGradients = returnGradients[injectT0, injectT1, injectDofs];
etaUnit = UnitVector[Length[injectDofs], Length[injectDofs]];

(* stage022's {1,1/2} are cited typed inputs, not re-derived fingerprints. *)
v0 = 1;
v1 = 1/2;
A0lead = clean[I v0 (a omega/cs) M0 (1 - T0dc)];
A1lead = clean[I v1 (a omega/cs)^3 D1 (1 - T1dc)];
expectedA0 = clean[I a omega M0 eps0/(cs (1 + eps0))];
expectedA1 = clean[I a^3 omega^3 D1 eps1/(2 cs^3 (1 + eps1))];
resA0 = clean[A0lead - expectedA0];
resA1 = clean[A1lead - expectedA1];
corruptV1A1 = clean[I (1/3) (a omega/cs)^3 D1 (1 - T1dc)];
corruptV1Residual = clean[corruptV1A1 - expectedA1];
corruptOrderA1 = clean[I v1 (a omega/cs)^2 D1 (1 - T1dc)];

zeroDim = {0, 0, 0};
dimScale[vector_, factor_] := factor vector;

dimensionAxisSlots = {{"L", 1}, {"M", 2}, {"T", 3}};
dimensionAxesLabel[] := StringRiffle[dimensionAxisSlots[[All, 1]], ","];
dimensionComponents[dimension_] := (dimension[[#[[2]]]] &) /@ dimensionAxisSlots;

printDimRecord[name_String, binding_] := Print[
  "DIM|axes=", dimensionAxesLabel[],
  "|name=", name,
  "|exponents=", ToString[InputForm[dimensionComponents[binding]]]
];

dimOf[expr_, dims_Association] := Module[{arguments, dimensions, base, power},
  Which[
    TrueQ[expr == 0] || NumericQ[expr], zeroDim,
    AtomQ[expr] && KeyExistsQ[dims, expr], dims[expr],
    AtomQ[expr], fail["missing dimension for " <> ToString[Unevaluated[expr], InputForm]],
    Head[expr] === Times, Total[dimOf[#, dims] & /@ (List @@ expr)],
    Head[expr] === Power,
      base = expr[[1]];
      power = expr[[2]];
      If[! NumericQ[power], fail["non-numeric dimension exponent"]];
      dimScale[dimOf[base, dims], power],
    Head[expr] === Plus,
      arguments = Select[List @@ expr, ! TrueQ[# == 0] &];
      dimensions = dimOf[#, dims] & /@ arguments;
      If[dimensions === {},
        zeroDim,
        If[Length[DeleteDuplicates[dimensions]] =!= 1,
          fail["dimension mismatch in sum " <> ToString[expr, InputForm]]
        ];
        First[dimensions]
      ],
    True, fail["unsupported dimension expression " <> ToString[expr, InputForm]]
  ]
];

baseDims = <|
  a -> {1, 0, 0}, cs -> {1, 0, -1}, omega -> {0, 0, -1},
  M0 -> {0, 1, -1}, D1 -> {1, 1, -1}, R0 -> {0, 1, -1}, R1 -> {1, 1, -1},
  D0 -> {-1, 1, -2}, K0c -> {0, 1, -2}, Keta -> {0, 1, -2},
  TOmega -> {0, 1, -2}, Z0ret -> {0, 1, -2}, Z1ret -> {0, 1, -2},
  OmegaU -> {0, 0, -1}, OmegaW -> {0, 0, -1}, Rmix -> {0, 0, -2},
  gU -> {-1/2, 1/2, -2}, gW -> {-1/2, 1/2, -2},
  etaNull -> zeroDim, gain0 -> zeroDim, gain1 -> zeroDim, qfree -> zeroDim
|>;

expectedDims = <|
  "A0" -> {0, 1, -1}, "A1" -> {1, 1, -1},
  "T0" -> zeroDim, "T1" -> zeroDim,
  "epsilon0" -> zeroDim, "epsilon1" -> zeroDim,
  "P0Physical" -> zeroDim
|>;

dimensionAudit[dims_Association, a0_, a1_, t0_, t1_, e0_, e1_, p0physical_] := Module[
  {expressions, computed, ok},
  expressions = <|
    "A0" -> a0, "A1" -> a1, "T0" -> t0, "T1" -> t1,
    "epsilon0" -> e0, "epsilon1" -> e1, "P0Physical" -> p0physical
  |>;
  computed = AssociationMap[
    Function[name,
      If[TrueQ[expressions[name] == 0], expectedDims[name], dimOf[expressions[name], dims]]
    ],
    Keys[expressions]
  ];
  ok = And @@ KeyValueMap[TrueQ[computed[#1] === #2] &, expectedDims];
  <|
    "Computed" -> computed,
    "Expected" -> expectedDims,
    "DimensionalOk" -> TrueQ[ok],
    "Verdict" -> If[TrueQ[ok], NOFAIL, FAILDIMENSIONAL],
    "MentionsQFree" -> ! FreeQ[Values[expressions], qfree]
  |>
];

baselineDimAudit = dimensionAudit[
  baseDims, A0lead, A1lead, T0dc, T1dc, eps0, eps1, baselinePort["P0Physical"]
];
corruptSourcedDims = Join[KeyDrop[baseDims, {M0}], <|M0 -> {1, 1, -1}|>];
corruptSourcedDimAudit = dimensionAudit[
  corruptSourcedDims, A0lead, A1lead, T0dc, T1dc, eps0, eps1, baselinePort["P0Physical"]
];
corruptFreeDims = Join[KeyDrop[baseDims, {qfree}], <|qfree -> {7, 0, 0}|>];
corruptFreeDimAudit = dimensionAudit[
  corruptFreeDims, A0lead, A1lead, T0dc, T1dc, eps0, eps1, baselinePort["P0Physical"]
];

emitDimensionRecords[] := (
  Print["DIMENSIONS|axes=", dimensionAxesLabel[]];
  printDimRecord["sourced_dims.a", baseDims[a]];
  printDimRecord["sourced_dims.c_s", baseDims[cs]];
  printDimRecord["sourced_dims.omega", baseDims[omega]];
  printDimRecord["sourced_dims.M0", baseDims[M0]];
  printDimRecord["sourced_dims.D1", baseDims[D1]];
  printDimRecord["sourced_dims.R0", baseDims[R0]];
  printDimRecord["sourced_dims.R1", baseDims[R1]];
  printDimRecord["sourced_dims.D0", baseDims[D0]];
  printDimRecord["sourced_dims.K0c", baseDims[K0c]];
  printDimRecord["sourced_dims.K_eta", baseDims[Keta]];
  printDimRecord["sourced_dims.T_Omega", baseDims[TOmega]];
  printDimRecord["sourced_dims.Z0_ret", baseDims[Z0ret]];
  printDimRecord["sourced_dims.Z1_ret", baseDims[Z1ret]];
  printDimRecord["sourced_dims.Omega_U", baseDims[OmegaU]];
  printDimRecord["sourced_dims.Omega_W", baseDims[OmegaW]];
  printDimRecord["sourced_dims.R_mix", baseDims[Rmix]];
  printDimRecord["sourced_dims.g_U", baseDims[gU]];
  printDimRecord["sourced_dims.g_W", baseDims[gW]];
  printDimRecord["sourced_dims.eta_null", baseDims[etaNull]];
  printDimRecord["sourced_dims.gain0", baseDims[gain0]];
  printDimRecord["sourced_dims.gain1", baseDims[gain1]];
  printDimRecord["sourced_dims.q_free", baseDims[qfree]];
  printDimRecord["computed_dims.A0", baselineDimAudit["Computed"]["A0"]];
  printDimRecord["computed_dims.A1", baselineDimAudit["Computed"]["A1"]];
  printDimRecord["computed_dims.T0", baselineDimAudit["Computed"]["T0"]];
  printDimRecord["computed_dims.T1", baselineDimAudit["Computed"]["T1"]];
  printDimRecord["computed_dims.epsilon0", baselineDimAudit["Computed"]["epsilon0"]];
  printDimRecord["computed_dims.epsilon1", baselineDimAudit["Computed"]["epsilon1"]];
  printDimRecord[
    "computed_dims.P0_physical",
    baselineDimAudit["Computed"]["P0Physical"]
  ]
);

makeItem[tags_List, computedClass_String] := <|
  "Tags" -> tags,
  "ComputedClass" -> computedClass
|>;

buildProvenance[mode_String, underdetermined_] := Module[{items, epsilonClass},
  epsilonClass = If[TrueQ[underdetermined], "deferred_branch_data", "derived_in_gate"];
  items = <|
    "stage022_v0_v1_coefficients" -> makeItem[
      {"stage022_cross_l_fingerprints", "checkable_typed_input"}, "cited_earned_input"
    ],
    "forward_T0_T1_transfer_map" -> makeItem[
      {"continuity_partition", "forward_from_Z_over_K"}, "derived_in_gate"
    ],
    "forward_A0_A1_residual_map" -> makeItem[
      {"stage022_coefficients", "pathA29_independent_form"}, "derived_in_gate"
    ],
    "ell2_P0_map" -> makeItem[{"stage017_grouped_port_kernel"}, "derived_in_gate"],
    "epsilon_eff_magnitude" -> makeItem[{"native_constructive_nullspace"}, epsilonClass],
    "pathA29_premise" -> makeItem[{"Z_is_premise", "boundary_dof_none"}, "cited_earned_input"],
    "pathA28_side_conditions" -> makeItem[
      {"stage008_R0_minus_M0", "stage008_R1_minus_D1"}, "external_bridge_input"
    ],
    "stage022_gate4_result" -> makeItem[
      {"quad_regression_false", "stage022_done"}, "cited_earned_input"
    ],
    "time_convention" -> makeItem[{"exp_minus_i_omega_t"}, "convention"]
  |>;
  Switch[mode,
    "assert_not_derive",
      items["forward_T0_T1_transfer_map"] = Join[
        items["forward_T0_T1_transfer_map"], <|"EmittedClass" -> "asserted_literal"|>
      ],
    "emit_epsilon",
      items["epsilon_eff_magnitude"] = Join[
        items["epsilon_eff_magnitude"],
        <|"EmittedClass" -> "derived_in_gate", "ConcreteValue" -> 1/2|>
      ],
    "asserted_selector",
      items["branch_selector"] = <|
        "Tags" -> {"selector_without_named_pde"},
        "ComputedClass" -> "counterfactual_control",
        "EmittedClass" -> "asserted_literal"
      |>,
    _, Null
  ];
  items = Map[
    Function[item,
      With[{emitted = Lookup[item, "EmittedClass", item["ComputedClass"]]},
        Join[
          item,
          <|
            "EmittedClass" -> emitted,
            "ClassMatchesComputed" -> TrueQ[emitted === item["ComputedClass"]]
          |>
        ]
      ]
    ],
    items
  ];
  <|
    "Items" -> items,
    "Ok" -> TrueQ[And @@ Lookup[Values[items], "ClassMatchesComputed"]],
    "NoConcreteEpsilonMagnitude" -> ! KeyExistsQ[items["epsilon_eff_magnitude"], "ConcreteValue"]
  |>
];

failurePriority = {
  {"Decoupled", FAILDECOUPLED},
  {"Tautological", FAILTAUTOLOGICAL},
  {"Dimensional", FAILDIMENSIONAL},
  {"NoConsistentReturn", FAILNOCONSISTENTRETURN},
  {"Overcancel", FAILOVERCANCEL},
  {"EpsilonMismatch", FAILEPSILONMISMATCH},
  {"Underdetermined", FAILUNDERDETERMINED}
};

verdictFromFlags[flags_Association] := Module[{firstFailure},
  firstFailure = SelectFirst[
    failurePriority,
    TrueQ[Lookup[flags, #[[1]], False]] &,
    Missing["NoFailure"]
  ];
  If[MissingQ[firstFailure], CROSSLRESIDUALPREDICTION, firstFailure[[2]]]
];

positiveBoundedPair[t0_, t1_, e0_, e1_] := TrueQ[
  FullSimplify[e0 > 0 && e1 > 0 && 0 < t0 < 1 && 0 < t1 < 1, $Assumptions]
];

caseFor[mode_String] := caseFor[mode] = Module[
  {
    z0, z1, t0, t1, e0, e1, coeff1, power1, a0, a1, expected0, expected1,
    form0, form1, order0, order1, caseBasis, moving, dims, dimAudit,
    provenance, decoupled, positive, overcancel, noConsistent, flags, verdict
  },
  z0 = Z0ret;
  z1 = Z1ret;
  coeff1 = 1/2;
  power1 = 3;
  caseBasis = basis;
  dims = baseDims;
  Switch[mode,
    "selector", z0 = K0c; z1 = K1dc; caseBasis = selectorBasis,
    "asserted_selector", z0 = K0c; z1 = K1dc; caseBasis = selectorBasis,
    "wrong_sign", z0 = -Z0ret; z1 = -Z1ret,
    "perfect", z0 = 0; z1 = 0,
    "no_consistent", z0 = -2 K0c; z1 = -2 K1dc,
    "corrupt_v1", coeff1 = 1/3,
    "corrupt_order", power1 = 2,
    "corrupt_dimension", dims = corruptSourcedDims,
    _, Null
  ];
  t0 = clean[K0c/(K0c + z0)];
  t1 = clean[K1dc/(K1dc + z1)];
  e0 = clean[z0/K0c];
  e1 = clean[z1/K1dc];
  If[mode === "decouple", t0 = clean[gain0 t0]; t1 = clean[gain1 t1]];
  a0 = clean[I (a omega/cs) M0 (1 - t0)];
  a1 = clean[I coeff1 (a omega/cs)^power1 D1 (1 - t1)];
  expected0 = clean[I a omega M0 e0/(cs (1 + e0))];
  expected1 = clean[I a^3 omega^3 D1 e1/(2 cs^3 (1 + e1))];
  form0 = TrueQ[clean[a0 - expected0] === 0];
  form1 = TrueQ[clean[a1 - expected1] === 0];
  order0 = Exponent[Cancel[Together[a0]], omega];
  order1 = Exponent[Cancel[Together[a1]], omega];
  moving = returnImageRank[caseBasis, returnGradients[t0, t1, rankDofs]];
  dimAudit = dimensionAudit[
    dims, a0, a1, t0, t1, e0, e1, baselinePort["P0Physical"]
  ];
  provenance = buildProvenance[mode, moving > 0];
  decoupled = TrueQ[
    ! FreeQ[t0, gain0] && FreeQ[t0, gain1] &&
    ! FreeQ[t1, gain1] && FreeQ[t1, gain0] &&
    FreeQ[P0port, gain0 | gain1] &&
    clean[D[t0, gain0]] =!= 0 && clean[D[t1, gain1]] =!= 0
  ];
  positive = positiveBoundedPair[t0, t1, e0, e1];
  overcancel = TrueQ[clean[e0] === 0 && clean[e1] === 0];
  noConsistent = TrueQ[clean[e0 + 2] === 0 && clean[e1 + 2] === 0 && ! positive];
  flags = <|
    "Decoupled" -> decoupled,
    "Tautological" -> ! provenance["Ok"],
    "Dimensional" -> ! dimAudit["DimensionalOk"],
    "NoConsistentReturn" -> noConsistent,
    "Overcancel" -> overcancel,
    "EpsilonMismatch" -> (! positive || ! form0 || ! form1 || order0 =!= 1 || order1 =!= 3),
    "Underdetermined" -> TrueQ[moving > 0]
  |>;
  verdict = verdictFromFlags[flags];
  <|
    "Mode" -> mode, "Z0" -> z0, "Z1" -> z1, "T0" -> t0, "T1" -> t1,
    "Epsilon0" -> e0, "Epsilon1" -> e1, "A0" -> a0, "A1" -> a1,
    "ExpectedA0" -> expected0, "ExpectedA1" -> expected1,
    "ResidualA0" -> clean[a0 - expected0], "ResidualA1" -> clean[a1 - expected1],
    "A0Order" -> order0, "A1Order" -> order1,
    "ReturnMovingNullity" -> moving, "Dimension" -> dimAudit,
    "Provenance" -> provenance, "SelectorProvenance" -> selectorMetadata[mode],
    "PathA29Premise" -> pathA29Premise, "Flags" -> flags, "Verdict" -> verdict
  |>
];

dynamicAblation[mode_String, expected_String] := Module[
  {withCase, withoutCase, neutralMode, neutralCase, withVerdict, withoutVerdict},
  withCase = caseFor[mode];
  withoutCase = caseFor["baseline"];
  neutralMode = "neutralized_independent_rerun";
  neutralCase = caseFor[neutralMode];
  withVerdict = withCase["Verdict"];
  withoutVerdict = withoutCase["Verdict"];
  <|
    "VerdictTrace" -> {withVerdict, withoutVerdict},
    "WithMutation" -> withVerdict,
    "WithoutMutation" -> withoutVerdict,
    "ComputedVerdictsDiffer" -> TrueQ[withVerdict =!= withoutVerdict],
    "ExpectedFailFires" -> TrueQ[withVerdict === expected && withoutVerdict =!= expected],
    "NeutralizedVerdict" -> neutralCase["Verdict"],
    "NeutralizedMutationFires" -> TrueQ[neutralCase["Verdict"] === expected],
    "NeutralizedIndependentRerun" -> TrueQ[
      neutralMode =!= "baseline" && neutralCase["Mode"] === neutralMode && neutralCase =!= withoutCase
    ]
  |>
];

checkFailureAblation[label_String, mode_String, expected_String] := Module[{audit},
  audit = dynamicAblation[mode, expected];
  Print["  ", label, ": verdict_trace=", audit["VerdictTrace"]];
  expectZero[
    label <> " dynamic rerun with mutation reaches " <> expected,
    verdictResidual[audit["WithMutation"], expected]
  ];
  expectZero[
    label <> " dynamic rerun without mutation returns earned native-nullspace FAIL",
    verdictResidual[audit["WithoutMutation"], FAILUNDERDETERMINED]
  ];
  expectBool[label <> " two computed verdicts differ", audit["ComputedVerdictsDiffer"]];
  expectBool[label <> " own expected token fires", audit["ExpectedFailFires"]];
  expectBool[
    label <> " independently rerun neutralized mutation does not fire",
    audit["NeutralizedIndependentRerun"] && ! audit["NeutralizedMutationFires"] &&
      audit["NeutralizedVerdict"] === audit["WithoutMutation"]
  ];
  TrueQ[
    audit["ExpectedFailFires"] && audit["NeutralizedIndependentRerun"] &&
      ! audit["NeutralizedMutationFires"]
  ]
];

portDependency[corrupt_] := Module[{candidate, baselineRow, candidateRow},
  candidate = portKernel[If[TrueQ[corrupt], 2, 1]];
  baselineRow = First[constraintJacobian[{baselinePort["P0Raw"]}, rankDofs]];
  candidateRow = First[constraintJacobian[{candidate["P0Raw"]}, rankDofs]];
  <|
    "P0Changes" -> ! TrueQ[clean[candidate["P0Raw"] - baselinePort["P0Raw"]] === 0],
    "Ell2RowChanges" -> ! TrueQ[clean[candidateRow - baselineRow] === ConstantArray[0, Length[rankDofs]]]
  |>
];

runNativeRankAndSelector[] := Module[{witnesses},
  subheading["Constructive native NullSpace and return-map image"];
  Print["  GENERATOR_DOFS=", rankDofs];
  Print["  symbolic constraint Jacobian dimensions=", Dimensions[Jbase]];
  Print["  rank0=", rank0];
  Print["  native_nullspace_dimension=", nativeNullity];
  Print["  return_moving_nullity=", returnMovingNullity];
  Print["  native_null_moves_return=", boolText[returnMovingNullity > 0]];
  expectBool["constraint Jacobian has exactly 3 genuine rows and 11 columns", Dimensions[Jbase] === {3, 11}];
  expectZero["native MatrixRank[Jbase] gives rank0=3", rank0 - 3];
  expectZero["constructive Length[NullSpace[Jbase]] gives native nullity 8", nativeNullity - 8];
  expectZero["constructive return image MatrixRank[basis.Transpose[Greturn]] is 2", returnMovingNullity - 2];
  expectBool["native constructive nullspace moves return", returnMovingNullity > 0];
  witnesses = {
    <|"Dof" -> "Z0_ret", "Vector" -> z0Unit|>,
    <|"Dof" -> "Z1_ret", "Vector" -> z1Unit|>
  };
  Scan[
    Function[witness,
      With[
        {
          constraintResiduals = clean[Jbase . witness["Vector"]],
          deltaReturn = clean[Greturn . witness["Vector"]]
        },
        Print[
          "  witness ", witness["Dof"], ": preserves_every_constraint=",
          boolText[constraintResiduals === ConstantArray[0, Length[Jbase]]],
          ", delta_return=", deltaReturn
        ];
        expectBool[
          witness["Dof"] <> " appended unit vector preserves every computed Jacobian row",
          constraintResiduals === ConstantArray[0, Length[Jbase]]
        ];
        expectBool[witness["Dof"] <> " unit witness moves T0/T1", deltaReturn =!= {0, 0}]
      ]
    ],
    witnesses
  ];

  subheading["Isolated constructive-nullspace teeth and selector control"];
  Print[
    "  isolated nullspaces: inject=", Length[injectBasis],
    "; one-row return-moving=", oneSelectorReturnMovingNullity,
    "; selector rank=", selectorRank, ", native=", selectorNativeNullity,
    ", return-moving=", selectorReturnMovingNullity
  ];
  Print["  selector_rank0=", selectorRank];
  Print["  selector_native_nullspace_dimension=", selectorNativeNullity];
  Print["  selector_return_moving_nullity=", selectorReturnMovingNullity];
  expectZero["inject_null adds etaNull and changes native nullity 8->9", Length[injectBasis] - 9];
  expectBool["injected etaNull unit direction preserves every base constraint", clean[constraintJacobian[baseConstraints, injectDofs] . etaUnit] === ConstantArray[0, 3]];
  expectBool["injected etaNull direction moves the injected return map", clean[injectGradients . etaUnit] =!= {0, 0}];
  expectZero["one selector row gives native nullity 7", Length[oneSelectorBasis] - 7];
  expectZero["one selector row changes return-moving nullity 2->1", oneSelectorReturnMovingNullity - 1];
  expectBool["one selector row leaves native_moves True", oneSelectorReturnMovingNullity > 0];
  expectZero["two selector rows raise rank 3->5", selectorRank - 5];
  expectZero["Length[NullSpace[Jselector]] changes native nullity 8->6", selectorNativeNullity - 6];
  expectBool["original Greturn annihilates the selector nullspace", selectorReturnImage === ConstantArray[0, {2, selectorNativeNullity}]];
  expectZero["two selector rows collapse return-moving nullity 2->0", selectorReturnMovingNullity];
  expectZero["selector makes T0=1/2", caseFor["selector"]["T0"] - 1/2];
  expectZero["selector makes T1=1/2", caseFor["selector"]["T1"] - 1/2];
  expectZero["baseline verdict is earned underdetermination FAIL", verdictResidual[caseFor["baseline"]["Verdict"], FAILUNDERDETERMINED]];
  expectZero["selector control verdict is CROSS_L_RESIDUAL_PREDICTION", verdictResidual[caseFor["selector"]["Verdict"], CROSSLRESIDUALPREDICTION]];
  expectBool["selector changes the two computed verdicts", caseFor["selector"]["Verdict"] =!= caseFor["baseline"]["Verdict"]];
  expectBool["neutered selector has no flip", caseFor["neutral"]["Verdict"] === caseFor["baseline"]["Verdict"]];
  expectBool[
    "selector provenance is control-only, not derived and not tautological",
    caseFor["selector"]["SelectorProvenance"]["ControlOnly"] &&
      ! caseFor["selector"]["SelectorProvenance"]["DerivedFromNamedPDE"] &&
      ! caseFor["selector"]["SelectorProvenance"]["Tautological"]
  ];
  Print[
    "  PROVENANCE (documentation; not counted): pathA_29 Z_is_premise=",
    boolText[caseFor["baseline"]["PathA29Premise"]["ZIsPremise"]],
    ", boundary_dof=", caseFor["baseline"]["PathA29Premise"]["BoundaryDof"]
  ]
];

runResidualTeeth[] := Module[{v1Flag, orderFlag},
  subheading["Forward A0/A1 residuals versus independent pathA_29 form"];
  Print["  cited stage022 coefficients: v0=", v0, ", v1=", v1];
  Print["  A0=", fmt[A0lead]];
  Print["  A1=", fmt[A1lead]];
  Print["  A0_residual_to_pathA29_form=", resA0];
  Print["  A1_residual_to_pathA29_form=", resA1];
  Print["  A0_order=omega^", Exponent[A0lead, omega], "; A1_order=omega^", Exponent[A1lead, omega]];
  expectZero["forward A0 matches independently built pathA_29 form", resA0];
  expectZero["forward A1 matches independently built pathA_29 form", resA1];
  expectZero["A0 realized omega order is 1", Exponent[A0lead, omega] - 1];
  expectZero["A1 realized omega order is 3", Exponent[A1lead, omega] - 3];
  selfConsistency[
    "T/epsilon identities are forward from Z/K",
    clean[T0dc - 1/(1 + eps0)] === 0 && clean[T1dc - 1/(1 + eps1)] === 0
  ];
  expectFail["stage022 v1 consumption tooth 1/2->1/3 fires A1 form", corruptV1Residual];
  expectFail["A1 order-realization corruption omega^3->omega^2 fires", Exponent[corruptOrderA1, omega] - 3];
  v1Flag = checkFailureAblation["v1 consumption", "corrupt_v1", FAILEPSILONMISMATCH];
  orderFlag = checkFailureAblation["A1 order", "corrupt_order", FAILEPSILONMISMATCH];
  <|"v1_consumption" -> v1Flag, "A1_order" -> orderFlag|>
];

runDimensionalGate[] := Module[{dimFlag},
  subheading["Exact (L,M,T) dimensional gate and free-carrier control"];
  KeyValueMap[
    Function[{name, expected},
      Print["  [", name, "]=", baselineDimAudit["Computed"][name], " expected=", expected];
      expectBool["[" <> name <> "] matches its sourced expected dimension", baselineDimAudit["Computed"][name] === expected]
    ],
    expectedDims
  ];
  Print["  dimensional_ok=", boolText[baselineDimAudit["DimensionalOk"]]];
  Print["  corrupt_sourced_verdict=", corruptSourcedDimAudit["Verdict"]];
  Print["  corrupt_free_carrier_verdict=", corruptFreeDimAudit["Verdict"]];
  expectBool["baseline dimensional_ok is True", baselineDimAudit["DimensionalOk"]];
  expectZero["corrupt sourced [M0]+=L reaches FAIL_DIMENSIONAL", verdictResidual[corruptSourcedDimAudit["Verdict"], FAILDIMENSIONAL]];
  expectZero["corrupt free qfree dimension reaches NO_FAIL locally", verdictResidual[corruptFreeDimAudit["Verdict"], NOFAIL]];
  expectBool["qfree appears in no checked expression", ! corruptFreeDimAudit["MentionsQFree"]];
  dimFlag = checkFailureAblation["3f corrupt sourced dimension", "corrupt_dimension", FAILDIMENSIONAL];
  <|"3f" -> dimFlag|>
];

runFirewall[] := Module[{base, asserted, epsilonEmission, assertedSelector, flags},
  subheading["Strengthened provenance firewall and anti-back-solve teeth"];
  base = caseFor["baseline"]["Provenance"];
  expectBool["baseline class_matches_computed holds for every item", base["Ok"]];
  Print[
    "  PROVENANCE (documentation; not counted): stage022 {1,1/2} class=",
    base["Items"]["stage022_v0_v1_coefficients"]["ComputedClass"]
  ];
  Print[
    "  PROVENANCE (documentation; not counted): ell2_P0_map tags=",
    base["Items"]["ell2_P0_map"]["Tags"]
  ];
  expectBool["epsilon magnitude computed class is deferred_branch_data", base["Items"]["epsilon_eff_magnitude"]["ComputedClass"] === "deferred_branch_data"];
  expectBool["baseline structurally emits no concrete epsilon magnitude/value", base["NoConcreteEpsilonMagnitude"]];
  asserted = caseFor["assert_not_derive"]["Provenance"];
  expectBool[
    "assert_not_derive is rewired to genuinely-023-derived transfer map",
    ! asserted["Items"]["forward_T0_T1_transfer_map"]["ClassMatchesComputed"] &&
      asserted["Items"]["stage022_v0_v1_coefficients"]["ClassMatchesComputed"]
  ];
  epsilonEmission = caseFor["emit_epsilon"]["Provenance"];
  expectBool[
    "epsilon-emission keeps computed deferred class but emits a concrete value",
    epsilonEmission["Items"]["epsilon_eff_magnitude"]["ComputedClass"] === "deferred_branch_data" &&
      KeyExistsQ[epsilonEmission["Items"]["epsilon_eff_magnitude"], "ConcreteValue"]
  ];
  assertedSelector = caseFor["asserted_selector"]["Provenance"];
  expectBool["asserted selector is separately class-mismatched", ! assertedSelector["Ok"]];
  flags = <|
    "assert_not_derive" -> checkFailureAblation["3g assert_not_derive", "assert_not_derive", FAILTAUTOLOGICAL],
    "asserted_selector" -> checkFailureAblation["3g asserted-selector", "asserted_selector", FAILTAUTOLOGICAL],
    "emit_epsilon_magnitude_as_derived" -> checkFailureAblation["3g emit_epsilon_magnitude_as_derived", "emit_epsilon", FAILTAUTOLOGICAL]
  |>;
  flags
];

runTransferAndPortProbes[] := Module[{flags, r1Mutated, r1Neutral},
  subheading["Dynamic transfer/return-family, decoupling, and R1 probes"];
  flags = <|
    "3a" -> checkFailureAblation["3a decouple_knobs", "decouple", FAILDECOUPLED],
    "3c" -> checkFailureAblation["3c wrong_sign_return", "wrong_sign", FAILEPSILONMISMATCH],
    "3d" -> checkFailureAblation["3d perfect_return", "perfect", FAILOVERCANCEL],
    "3h" -> checkFailureAblation["3h no_consistent_return", "no_consistent", FAILNOCONSISTENTRETURN]
  |>;
  r1Mutated = portDependency[True];
  r1Neutral = portDependency[False];
  Print["  R1 corrupt-port signals=", r1Mutated, "; neutralized=", r1Neutral];
  expectBool["R1 corrupt OmegaU->2*OmegaU changes P0_raw", r1Mutated["P0Changes"]];
  expectBool["R1 corrupt OmegaU->2*OmegaU changes ell=2 rank row", r1Mutated["Ell2RowChanges"]];
  expectBool["R1 neutralized port mutation changes neither object", ! Or @@ Values[r1Neutral]];
  Join[flags, <|"R1" -> TrueQ[And @@ Values[r1Mutated] && ! Or @@ Values[r1Neutral]]|>]
];

runAritySelfCheck[] := Module[{coreResults},
  subheading["Wolfram arity self-check and unevaluated-leakage scan"];
  expectBool["arity constraintJacobian[constraints,dofs] returns a 3x11 Jacobian", Dimensions[constraintJacobian[baseConstraints, rankDofs]] === {3, 11}];
  expectBool["arity constructNullBasis[constraints,dofs] returns eight basis vectors", Length[constructNullBasis[baseConstraints, rankDofs]] === 8];
  expectBool["arity returnGradients[t0,t1,dofs] returns a 2x11 map", Dimensions[returnGradients[T0dc, T1dc, rankDofs]] === {2, 11}];
  expectZero["arity returnImageRank[basis,gradients] returns 2", returnImageRank[basis, Greturn] - 2];
  expectBool["arity dimensionAudit[dims,A0,A1,T0,T1,e0,e1,P0] returns dimensional status", KeyExistsQ[baselineDimAudit, "DimensionalOk"]];
  expectBool["arity buildProvenance[mode,underdetermined] returns Items/Ok", And @@ (KeyExistsQ[buildProvenance["baseline", True], #] & /@ {"Items", "Ok"})];
  expectZero["arity verdictFromFlags[association] returns baseline FAIL", verdictResidual[verdictFromFlags[caseFor["baseline"]["Flags"]], FAILUNDERDETERMINED]];
  expectBool["arity dynamicAblation[mode,label] returns two distinct verdicts", dynamicAblation["wrong_sign", FAILEPSILONMISMATCH]["ComputedVerdictsDiffer"]];
  coreResults = {
    Jbase, basis, Greturn, selectorBasis, selectorReturnImage, injectBasis,
    A0lead, A1lead, expectedA0, expectedA1, resA0, resA1,
    baselineDimAudit, corruptSourcedDimAudit, corruptFreeDimAudit,
    caseFor["baseline"], caseFor["selector"], caseFor["wrong_sign"],
    portDependency[True]
  };
  expectBool[
    "no unevaluated authored helper, NullSpace, or MatrixRank call remains in core results",
    FreeQ[
      coreResults,
      constraintJacobian | constructNullBasis | returnGradients | returnImageRank |
      dimOf | dimensionAudit | buildProvenance | verdictFromFlags | caseFor |
      dynamicAblation | portDependency | NullSpace | MatrixRank
    ]
  ]
];

runVerdictSeamAndScope[probeFlags_Association] := Module[{baseline, selector},
  subheading["022/023 verdict seam and PROVENANCE consumption"];
  baseline = caseFor["baseline"];
  selector = caseFor["selector"];
  Print["  baseline flags=", baseline["Flags"]];
  expectBool["stage023 failure priority has no quad_regression rung", FreeQ[failurePriority, "QuadRegression" | "FAIL_QUAD_REGRESSION"]];
  expectBool["stage023 failure priority has no able_to_fail_bad rung", FreeQ[failurePriority, "AbleToFailBad"]];
  expectBool["actual probe results non-recursively establish able-to-fail", And @@ Values[probeFlags]];
  expectZero["baseline physics verdict is earned underdetermination FAIL", verdictResidual[baseline["Verdict"], FAILUNDERDETERMINED]];
  expectZero["selector control physics verdict is residual prediction", verdictResidual[selector["Verdict"], CROSSLRESIDUALPREDICTION]];
  expectBool["vestigial Delta and Sport symbols are dropped", FreeQ[Names["Global`*"], "Global`Delta" | "Global`Sport" | "Global`S_port"]];
  Print["  CHECKABLE: stage022 cited {v0=1,v1=1/2} feeds the forward A0/A1 map; v1 corruption fires."];
  Print["  CHECKABLE: stage009/010 pathA_29 epsilon/(1+epsilon) form is built independently."];
  Print["  PROVENANCE: pathA_29 Z_is_premise=True, boundary_dof=none is the keystone return premise."];
  Print["  PROVENANCE: stage008 M0/D1 and R0=-M0/R1=-D1; stage017 grouped-P2 P0_raw port kernel."];
  Print["  PROVENANCE CONTEXT: pathA_34-convention K0c/K1; stage013/017 scalar reduction remains pending."];
  Print["  PROVENANCE: c_s from R1 and a from CONV; stage022 Gate-4 quad_regression=False is cited, not rebuilt."];
  Print["  ZERO file I/O: standalone print-only audit; no scratch bridge or publication artifact is written."]
];

printVerdictLabels[] := (
  subheading["Verdict labels:"];
  Print["  earned physics: the Gate-5 named constraints leave two return-moving directions."];
  Print["  able-to-fail control: two counterfactual return equations annihilate the return map on the selector nullspace."];
  Print["  honest export: Gate 6 must supply two independent return equations or an equivalent return law."];
  Print["baseline_verdict=", caseFor["baseline"]["Verdict"]];
  Print["selector_control_verdict=", caseFor["selector"]["Verdict"]];
  Print["AUDIT_STATUS=PASS"];
  Print["PHYSICS_VERDICT=FAIL_UNDERDETERMINED_NOT_PREDICTIVE (2/2, COMPLETING)"]
);

runAll[] := Module[{residualFlags, dimensionalFlags, firewallFlags, transferFlags, allFlags},
  heading["ledger_stage023_nullspace_underdetermination_mathematica_audit"];
  Print["Target stem confirmed: ledger_stage023_nullspace_underdetermination"];
  Print["Engine: re-authored constructive NullSpace basis + return-image rank; exact, float-free, standalone, zero file I/O."];
  emitDimensionRecords[];
  runNativeRankAndSelector[];
  residualFlags = runResidualTeeth[];
  dimensionalFlags = runDimensionalGate[];
  firewallFlags = runFirewall[];
  transferFlags = runTransferAndPortProbes[];
  allFlags = Join[residualFlags, dimensionalFlags, firewallFlags, transferFlags];
  runAritySelfCheck[];
  runVerdictSeamAndScope[allFlags];
  printVerdictLabels[];
  0
];

result = Catch[
  runAll[],
  "ledgerStage023Failure",
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
