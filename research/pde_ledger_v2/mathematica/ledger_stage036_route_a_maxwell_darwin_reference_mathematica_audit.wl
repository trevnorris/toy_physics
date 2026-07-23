(* Ledger stage036 Mathematica audit: Route-A Maxwell--Darwin reference.

   Standalone, print-only, assert-zero, exact, and file-I/O-free.

   Genuinely independent route: this engine evaluates the analytically
   continued Riesz transform for k^-4 and direct Schwinger/Gaussian momentum
   moments for k_i k_j/k^4.  It assembles the transverse projector from those
   moment integrals, without differentiating a supplied radial seed.  The
   velocity anchor is organized by parallel/perpendicular decomposition, and
   the force follows from the differential geometry of R and n rather than a
   Cartesian Hessian/gradient port.

   Tooth-local runtime ablation uses LEDGER_STAGE036_MUTATION.
*)

ClearAll["Global`*"];
$HistoryLength = 0;

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE036_MUTATION";
activeMutation = Environment[mutationEnvironment];
If[!StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];

verdictToken = "MAXWELL_DARWIN_REFERENCE";
uncertifiedToken = "MAXWELL_DARWIN_REFERENCE_UNCERTIFIED";
tierToken = "tier_A_conditional";
electricR1 = "R1_REQUIRED(bc_selection)";
manifestDigestExpected =
  "f8b1569834c1c5cfd404fee4fba7bc49d3d7263f332e551270ad499fd3585cac";

toothOrder = {
  "BOOST_PROJECTOR",
  "BOOST_GENERAL_VELOCITIES",
  "BOOST_NEXT_ORDER",
  "TARGET_BLINDNESS",
  "DUAL_ENGINE_TERMS",
  "UNITS_RESTORED",
  "VERDICT_REDERIVATION",
  "SOURCE_TO_STAGE_MANIFEST"
};

ablationDescriptions = <|
  "BOOST_PROJECTOR" ->
    "double only the Schwinger-moment kernel checked by BOOST_PROJECTOR",
  "BOOST_GENERAL_VELOCITIES" ->
    "drop A_V only from the live parallel/perpendicular anchor",
  "BOOST_NEXT_ORDER" ->
    "raise only the claimed computed velocity order from 2 to 4",
  "TARGET_BLINDNESS" ->
    "inject a Route-B q_T^2/mu_R term only into the dependency object",
  "DUAL_ENGINE_TERMS" ->
    "drop radial_force_A2 only from the canonical computed-term inventory",
  "UNITS_RESTORED" ->
    "make R dimensionless only inside the formal unit rescaling",
  "VERDICT_REDERIVATION" ->
    "double a computed Schwinger projector input consumed only by verdict derivation",
  "SOURCE_TO_STAGE_MANIFEST" ->
    "drop one scoped-out row only from the canonical 35-tooth manifest"
|>;

raise[name_] := Throw[name, "ledgerStage036Failure"];

assertExact[name_, expression_] := Module[{reals},
  reals = Cases[Unevaluated[expression], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FIRST_FAILURE=", name];
    Print["FAIL  ", name, ": machine-real atom(s): ", InputForm[reals]];
    raise[name]
  ]
];

cleanZero[expression_] :=
  FullSimplify[expression] /. ConditionalExpression[0, _] -> 0;

expectZero[name_, residual_, evidence_: None] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FIRST_FAILURE=", name];
    If[activeMutation === name, Print["FIRED_AT_OWN_ASSERT=", name]];
    Print["FAIL  ", name, ": residual = ", InputForm[clean]];
    If[evidence =!= None, Print["      evidence = ", InputForm[evidence]]];
    raise[name]
  ]
];

expectBool[name_, condition_, evidence_: None] :=
  expectZero[name, If[TrueQ[condition], 0, 1], evidence];

section[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);

zeroQ[expression_] := TrueQ[Factor[Together[expression]] === 0];

unitSphereReduce[expression_] := Factor[
  FullSimplify[
    Factor[expression] /. nz^2 -> 1 - nx^2 - ny^2,
    R > 0 &&
      Element[
        {nx, ny, nz, v1x, v1y, v1z, v2x, v2y, v2z,
          s1, s2, aE, cGamma}, Reals] &&
      cGamma > 0
  ]
];

(* ---------------------------------------------------------------------- *)
(* Direct Riesz/Schwinger momentum-space construction.                     *)
(* ---------------------------------------------------------------------- *)

nvec = {nx, ny, nz};
velocity1 = {v1x, v1y, v1z};
velocity2 = {v2x, v2y, v2z};
velocitySymbols = Join[velocity1, velocity2];
dV = Expand[velocity1.velocity2];
v1n = Expand[velocity1.nvec];
v2n = Expand[velocity2.nvec];
aV = Expand[v1n v2n];

rieszPotential[alpha_] :=
  Gamma[3/2 - alpha] R^(2 alpha - 3)/
    (4^alpha Pi^(3/2) Gamma[alpha]);

coulombScalar = FullSimplify[rieszPotential[1], R > 0];
seedK4 = FullSimplify[rieszPotential[2], R > 0];

schwingerGaussian =
  Exp[-R^2/(4 tau)]/(4 Pi tau)^(3/2);
deltaMomentIntegrand = schwingerGaussian/2;
longitudinalMomentIntegrand = schwingerGaussian/(4 tau);

Print["PROGRESS MATHEMATICA STAGE036 SCHWINGER INTEGRALS START"];
deltaMoment = FullSimplify[
  Integrate[
    deltaMomentIntegrand,
    {tau, 0, Infinity},
    Assumptions -> R > 0,
    GenerateConditions -> False],
  R > 0];
longitudinalMoment = FullSimplify[
  Integrate[
    longitudinalMomentIntegrand,
    {tau, 0, Infinity},
    Assumptions -> R > 0,
    GenerateConditions -> False],
  R > 0];
Print["PROGRESS MATHEMATICA STAGE036 SCHWINGER INTEGRALS COMPLETE"];

kkMomentAtR =
  deltaMoment IdentityMatrix[3] -
    longitudinalMoment R^2 Outer[Times, nvec, nvec];
darwinKernel = Map[
  unitSphereReduce,
  coulombScalar IdentityMatrix[3] - kkMomentAtR,
  {2}];
expectedKernel =
  (IdentityMatrix[3] + Outer[Times, nvec, nvec])/(8 Pi R);
kernelResiduals =
  unitSphereReduce /@ Flatten[darwinKernel - expectedKernel];

(* ---------------------------------------------------------------------- *)
(* Independent parallel/perpendicular anchor and geometric force route.    *)
(* ---------------------------------------------------------------------- *)

velocity1Perpendicular = Expand[velocity1 - v1n nvec];
velocity2Perpendicular = Expand[velocity2 - v2n nvec];
velocityBilinear = unitSphereReduce[
  Expand[
    velocity1Perpendicular.velocity2Perpendicular + 2 aV]];
bilinearResidual = unitSphereReduce[velocityBilinear - (dV + aV)];

electricU0 = s1 s2 aE/(4 Pi R);
routeAU2 =
  Factor[-s1 s2 aE velocityBilinear/(8 Pi cGamma^2 R)];
expectedA2 =
  -s1 s2 aE (dV + aV)/(8 Pi cGamma^2 R);
fullAnchor = Factor[electricU0 + routeAU2];
expectedFullAnchor =
  Factor[electricU0 (1 - (dV + aV)/(2 cGamma^2))];

(* grad_n(V.n)=(V-(V.n)n)/R and grad_R(1/R)=-n/R^2. *)
gradientA = Expand[
  (v2n (velocity1 - v1n nvec) +
    v1n (velocity2 - v2n nvec))/R];
gradientShape = Map[
  unitSphereReduce,
  gradientA/R - (dV + aV) nvec/R^2];
forceA2 = Map[
  unitSphereReduce,
  s1 s2 aE/(8 Pi cGamma^2) gradientShape];
expectedForceA2 =
  s1 s2 aE/(8 Pi cGamma^2 R^2) *
    (v2n velocity1 + v1n velocity2 - (dV + 3 aV) nvec);
forceResiduals =
  unitSphereReduce /@ (forceA2 - expectedForceA2);
radialForceA2 = unitSphereReduce[nvec.forceA2];
expectedRadialForceA2 =
  -s1 s2 aE (dV + aV)/(8 Pi cGamma^2 R^2);
radialResidual =
  unitSphereReduce[radialForceA2 - expectedRadialForceA2];

(* ---------------------------------------------------------------------- *)
(* Runtime velocity-order object.                                          *)
(* ---------------------------------------------------------------------- *)

velocityScaledAnchor = Expand[
  fullAnchor /.
    Thread[velocitySymbols -> epsilonV velocitySymbols]];
velocityCoefficientRules =
  CoefficientRules[velocityScaledAnchor, {epsilonV}];
computedVelocityOrders = Sort[
  DeleteDuplicates[
    Cases[
      velocityCoefficientRules,
      Rule[{power_}, coefficient_] /; !zeroQ[coefficient] :> power]]];
computedVelocityOrder = Max[computedVelocityOrders];
nextUncomputedOrder = 4;

(* ---------------------------------------------------------------------- *)
(* Native formal-rescaling dimensional firewall on the real expressions.   *)
(* ---------------------------------------------------------------------- *)

unitRules[badRadius_] := Join[
  {
    R -> If[TrueQ[badRadius], R, lengthScale R],
    aE -> massScale lengthScale^3 timeScale^-2 aE,
    cGamma -> lengthScale timeScale^-1 cGamma,
    s1 -> s1,
    s2 -> s2,
    epsilonV -> epsilonV
  },
  Thread[
    velocitySymbols ->
      (lengthScale timeScale^-1 velocitySymbols)]
];

scalingResidual[
  expression_, requiredScale_, badRadius_] :=
  Factor[Together[
    (expression /. unitRules[badRadius]) -
      requiredScale expression]];

dimensionResiduals[badRadius_] := Join[
  scalingResidual[
    #, lengthScale^-1, badRadius] & /@ Flatten[darwinKernel],
  {
    scalingResidual[
      aE,
      massScale lengthScale^3 timeScale^-2,
      badRadius],
    scalingResidual[
      cGamma^2,
      lengthScale^2 timeScale^-2,
      badRadius]
  },
  scalingResidual[
    #, lengthScale timeScale^-1, badRadius] & /@ velocitySymbols,
  {
    scalingResidual[
      dV, lengthScale^2 timeScale^-2, badRadius],
    scalingResidual[
      aV, lengthScale^2 timeScale^-2, badRadius],
    scalingResidual[
      electricU0,
      massScale lengthScale^2 timeScale^-2,
      badRadius],
    scalingResidual[
      routeAU2,
      massScale lengthScale^2 timeScale^-2,
      badRadius],
    scalingResidual[
      fullAnchor,
      massScale lengthScale^2 timeScale^-2,
      badRadius]
  },
  scalingResidual[
    #,
    massScale lengthScale timeScale^-2,
    badRadius] & /@ forceA2,
  {
    scalingResidual[
      radialForceA2,
      massScale lengthScale timeScale^-2,
      badRadius],
    scalingResidual[s1, 1, badRadius],
    scalingResidual[s2, 1, badRadius],
    scalingResidual[epsilonV^2, 1, badRadius]
  }
];

(* ---------------------------------------------------------------------- *)
(* Runtime reference-only guard, term inventory, and verdict precedence.   *)
(* ---------------------------------------------------------------------- *)

forbiddenTargetSymbols =
  {qT, muR, rBA, rCone, deltaU, sealedLanding};

expectedTermInventory = {
  "darwin_kernel_9_components",
  "electric_U0",
  "route_A_U2",
  "full_anchor_UA",
  "force_A2_3_components",
  "radial_force_A2",
  "velocity_orders_0_2_next_4",
  "dimension_firewall",
  "verdict_precedence"
};

computedTermInventory[includeRadialForce_: True] := Module[{terms = {}},
  If[Dimensions[darwinKernel] === {3, 3},
    AppendTo[terms, "darwin_kernel_9_components"]];
  If[!FreeQ[electricU0, aE] && !FreeQ[electricU0, R],
    AppendTo[terms, "electric_U0"]];
  If[
    !FreeQ[routeAU2, aE] &&
    !FreeQ[routeAU2, cGamma] &&
    !FreeQ[routeAU2, R],
    AppendTo[terms, "route_A_U2"]];
  If[
    !FreeQ[fullAnchor, aE] &&
    zeroQ[unitSphereReduce[fullAnchor - expectedFullAnchor]],
    AppendTo[terms, "full_anchor_UA"]];
  If[Length[forceA2] === 3,
    AppendTo[terms, "force_A2_3_components"]];
  If[TrueQ[includeRadialForce] && !FreeQ[radialForceA2, R],
    AppendTo[terms, "radial_force_A2"]];
  If[computedVelocityOrders === {0, 2},
    AppendTo[terms, "velocity_orders_0_2_next_4"]];
  AppendTo[terms, "dimension_firewall"];
  AppendTo[terms, "verdict_precedence"];
  terms
];

deriveVerdict[
  projectorResiduals_List,
  anchorResidual_,
  orderIsTwo_] :=
  If[
    And @@ (zeroQ /@ projectorResiduals) &&
      zeroQ[anchorResidual] &&
      TrueQ[orderIsTwo],
    verdictToken,
    uncertifiedToken
  ];

(* Exact source-build order: all 35 teeth, no wildcard families. *)
sourceToothIDs = {
  "SOURCE_TRANSLATION_CONTINUITY",
  "SOURCE_NOT_IMPORTED",
  "SOURCE_BASIS",
  "PARITY_RW",
  "PARITY_PW",
  "PARITY_ROTATION",
  "PARITY_TIME_REVERSAL",
  "FIELD_IDENTITY_UNITS",
  "ACTION_KINETIC",
  "ACTION_COUPLING",
  "ACTION_STABILITY",
  "G0_DAMAGE",
  "ROUTE_INDEPENDENCE",
  "BOOST_PROJECTOR",
  "BOOST_GENERAL_VELOCITIES",
  "BOOST_NEXT_ORDER",
  "BOOST_COMMON_VELOCITY",
  "DIRECT_SOURCE",
  "DIRECT_PROJECTOR",
  "DIRECT_EXCHANGE_SIGN",
  "DIRECT_FALLOFF",
  "DIRECT_VELOCITY_ORDER",
  "COMPARE_COMPUTED",
  "DELTA_RATIO",
  "CONE_RATIO",
  "QMAG_R1",
  "UNITS_RESTORED",
  "ACTIVE_FLUX_CAVEAT",
  "HOOK_LORENTZ",
  "LEDGER_READY_ROW",
  "TRUTH_TOTALITY",
  "TRUTH_PRECEDENCE",
  "LANDING_OWNERSHIP",
  "TARGET_BLINDNESS",
  "DUAL_ENGINE_TERMS"
};

sourceManifest = {
  {"SOURCE_TRANSLATION_CONTINUITY", "SCOPED_OUT", "STAGE035_V2_DONE"},
  {"SOURCE_NOT_IMPORTED", "SCOPED_OUT", "STAGE035_V2_DONE"},
  {"SOURCE_BASIS", "SCOPED_OUT", "STAGE035_V2_DONE"},
  {"PARITY_RW", "SCOPED_OUT", "STAGE035_V2_DONE"},
  {"PARITY_PW", "SCOPED_OUT", "STAGE035_V2_DONE"},
  {"PARITY_ROTATION", "SCOPED_OUT", "STAGE035_V2_DONE"},
  {"PARITY_TIME_REVERSAL", "SCOPED_OUT", "STAGE035_V2_DONE"},
  {"FIELD_IDENTITY_UNITS", "SCOPED_OUT", "STAGE034_V1_DONE"},
  {"ACTION_KINETIC", "SCOPED_OUT", "STAGE034_V1_DONE"},
  {"ACTION_COUPLING", "SCOPED_OUT", "STAGE034_V1_DONE"},
  {"ACTION_STABILITY", "SCOPED_OUT", "STAGE034_V1_DONE"},
  {"G0_DAMAGE", "SCOPED_OUT", "STAGE034_V1_DONE"},
  {"ROUTE_INDEPENDENCE", "SCOPED_OUT", "STAGE037_V4"},
  {"BOOST_PROJECTOR", "REPLACED_BY_STRONGER",
    "STAGE036_INDEPENDENT_PROJECTOR_RECONSTRUCTION"},
  {"BOOST_GENERAL_VELOCITIES", "REPLACED_BY_STRONGER",
    "STAGE036_ANCHOR_FORCE_RADIAL_RECONSTRUCTION"},
  {"BOOST_NEXT_ORDER", "PRESERVED",
    "STAGE036_RUNTIME_VELOCITY_ORDER"},
  {"BOOST_COMMON_VELOCITY", "SCOPED_OUT", "STAGE037_V4"},
  {"DIRECT_SOURCE", "SCOPED_OUT", "STAGE037_V4"},
  {"DIRECT_PROJECTOR", "SCOPED_OUT", "STAGE037_V4"},
  {"DIRECT_EXCHANGE_SIGN", "SCOPED_OUT", "STAGE037_V4"},
  {"DIRECT_FALLOFF", "SCOPED_OUT", "STAGE037_V4"},
  {"DIRECT_VELOCITY_ORDER", "SCOPED_OUT", "STAGE037_V4"},
  {"COMPARE_COMPUTED", "SCOPED_OUT", "STAGE037_V4"},
  {"DELTA_RATIO", "SCOPED_OUT", "STAGE037_V4"},
  {"CONE_RATIO", "SCOPED_OUT", "STAGE037_V4"},
  {"QMAG_R1", "SCOPED_OUT", "STAGE037_V4"},
  {"UNITS_RESTORED", "REPLACED_BY_STRONGER",
    "STAGE036_EXPRESSION_DIMENSION_FIREWALL"},
  {"ACTIVE_FLUX_CAVEAT", "SCOPED_OUT", "STAGE038_V5"},
  {"HOOK_LORENTZ", "SCOPED_OUT", "STAGE038_V5"},
  {"LEDGER_READY_ROW", "SCOPED_OUT", "STAGE034_V1_DONE"},
  {"TRUTH_TOTALITY", "SCOPED_OUT", "STAGE038_V5"},
  {"TRUTH_PRECEDENCE", "SCOPED_OUT", "STAGE038_V5"},
  {"LANDING_OWNERSHIP", "SCOPED_OUT", "STAGE038_V5"},
  {"TARGET_BLINDNESS", "PRESERVED",
    "STAGE036_REFERENCE_DEPENDENCY_OBJECT"},
  {"DUAL_ENGINE_TERMS", "REPLACED_BY_STRONGER",
    "STAGE036_CANONICAL_TERM_INVENTORY"}
};

expectedManifestCounts = KeySort@<|
  "PRESERVED" -> 2,
  "REPLACED_BY_STRONGER" -> 4,
  "SCOPED_OUT" -> 29
|>;

expectedInScope = {
  "BOOST_PROJECTOR",
  "BOOST_GENERAL_VELOCITIES",
  "BOOST_NEXT_ORDER",
  "UNITS_RESTORED",
  "TARGET_BLINDNESS",
  "DUAL_ENGINE_TERMS"
};

expectedScopedOut = {
  "SOURCE_TRANSLATION_CONTINUITY",
  "SOURCE_NOT_IMPORTED",
  "SOURCE_BASIS",
  "PARITY_RW",
  "PARITY_PW",
  "PARITY_ROTATION",
  "PARITY_TIME_REVERSAL",
  "FIELD_IDENTITY_UNITS",
  "ACTION_KINETIC",
  "ACTION_COUPLING",
  "ACTION_STABILITY",
  "G0_DAMAGE",
  "ROUTE_INDEPENDENCE",
  "BOOST_COMMON_VELOCITY",
  "DIRECT_SOURCE",
  "DIRECT_PROJECTOR",
  "DIRECT_EXCHANGE_SIGN",
  "DIRECT_FALLOFF",
  "DIRECT_VELOCITY_ORDER",
  "COMPARE_COMPUTED",
  "DELTA_RATIO",
  "CONE_RATIO",
  "QMAG_R1",
  "ACTIVE_FLUX_CAVEAT",
  "HOOK_LORENTZ",
  "LEDGER_READY_ROW",
  "TRUTH_TOTALITY",
  "TRUTH_PRECEDENCE",
  "LANDING_OWNERSHIP"
};

lexicographicCodeLess[left_, right_] := Module[{index, limit},
  limit = Min[Length[left], Length[right]];
  index = SelectFirst[
    Range[limit], left[[#]] =!= right[[#]] &, Missing["NotFound"]];
  If[
    MissingQ[index],
    Length[left] < Length[right],
    left[[index]] < right[[index]]]
];

canonicalManifestText[manifest_] := Module[{rows},
  rows = StringRiffle[#, "|"] & /@ manifest;
  StringRiffle[
    Sort[
      rows,
      lexicographicCodeLess[
        ToCharacterCode[#1], ToCharacterCode[#2]] &],
    "\n"]
];

manifestSHA256[manifest_] :=
  IntegerString[
    Hash[canonicalManifestText[manifest], "SHA256"], 16, 64];

ok = Catch[
  If[
    activeMutation =!= "" &&
      !MemberQ[toothOrder, activeMutation],
    Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
    Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
    raise["UNKNOWN_MUTATION"]
  ];

  Print[
    "ledger_stage036_route_a_maxwell_darwin_reference Mathematica audit"];
  Print[
    "ROUTE=Riesz analytic continuation + direct Schwinger/Gaussian k-space moments + perpendicular/parallel anchor"];
  Print[
    "INDEPENDENCE=no seeded Hessian; no Cartesian potential differentiation"];
  Print["FILE_IO=none; CROSS_ENGINE_COMPARE=none"];
  Print["PROGRESS MATHEMATICA STAGE036 ASSERTIONS START"];
  If[activeMutation =!= "",
    Print["ACTIVE_MUTATION=", activeMutation];
    Print[
      "MUTATED_PRIMITIVE=",
      ablationDescriptions[activeMutation]]
  ];

  section[
    "Inverse-FT transverse projector from Riesz/Schwinger momentum moments"];
  liveKernel = If[
    activeMutation === "BOOST_PROJECTOR",
    2 darwinKernel,
    darwinKernel];
  liveProjectorResiduals =
    unitSphereReduce /@
      Flatten[liveKernel - expectedKernel];
  projectorOK =
    zeroQ[coulombScalar - 1/(4 Pi R)] &&
    zeroQ[seedK4 + R/(8 Pi)] &&
    zeroQ[deltaMoment - 1/(8 Pi R)] &&
    zeroQ[longitudinalMoment - 1/(8 Pi R^3)] &&
    And @@ (zeroQ /@ liveProjectorResiduals);
  expectBool[
    "BOOST_PROJECTOR",
    projectorOK,
    <|
      "RieszK4Seed" -> seedK4,
      "DeltaMoment" -> deltaMoment,
      "LongitudinalMoment" -> longitudinalMoment,
      "ComponentResiduals" -> liveProjectorResiduals
    |>];
  Print[
    "      Riesz F^-1[k^-4]=-R/(8*pi); Schwinger moments=1/(8*pi*R),1/(8*pi*R^3)"];
  Print[
    "      I_ij=(delta_ij+n_i*n_j)/(8*pi*R), nine component residuals zero"];

  section[
    "Independent-velocity anchor by perpendicular/parallel decomposition"];
  liveU2 = If[
    activeMutation === "BOOST_GENERAL_VELOCITIES",
    -s1 s2 aE dV/(8 Pi cGamma^2 R),
    routeAU2];
  liveFullAnchor = Factor[electricU0 + liveU2];
  anchorResidual =
    unitSphereReduce[liveU2 - expectedA2];
  fullAnchorResidual =
    unitSphereReduce[liveFullAnchor - expectedFullAnchor];
  generalVelocityOK =
    zeroQ[bilinearResidual] &&
    zeroQ[anchorResidual] &&
    zeroQ[fullAnchorResidual] &&
    And @@ (zeroQ /@ forceResiduals) &&
    zeroQ[radialResidual];
  expectBool[
    "BOOST_GENERAL_VELOCITIES",
    generalVelocityOK,
    <|
      "BilinearResidual" -> bilinearResidual,
      "AnchorResidual" -> anchorResidual,
      "FullAnchorResidual" -> fullAnchorResidual,
      "ForceResiduals" -> forceResiduals,
      "RadialResidual" -> radialResidual
    |>];
  Print["      D_V=V1.V2; A_V=(V1.n)(V2.n)"];
  Print[
    "      U_A=(s1*s2*A_E/(4*pi*R))*(1-(D_V+A_V)/(2*c_gamma^2))"];
  Print[
    "      F_A2=(s1*s2*A_E/(8*pi*c_gamma^2*R^2))*[(V2.n)V1+(V1.n)V2-(D_V+3A_V)n]"];
  Print[
    "      F_A2,r=-s1*s2*A_E*(D_V+A_V)/(8*pi*c_gamma^2*R^2)"];

  section[
    "Explicit computed orders and named uncomputed remainder"];
  claimedOrder = If[
    activeMutation === "BOOST_NEXT_ORDER",
    4,
    computedVelocityOrder];
  orderOK =
    computedVelocityOrders === {0, 2} &&
    claimedOrder === 2 &&
    zeroQ[Coefficient[velocityScaledAnchor, epsilonV, 4]];
  expectBool[
    "BOOST_NEXT_ORDER",
    orderOK,
    <|
      "ComputedOrders" -> computedVelocityOrders,
      "ClaimedOrder" -> claimedOrder,
      "NextUncomputed" -> nextUncomputedOrder,
      "Epsilon4Coefficient" ->
        Coefficient[velocityScaledAnchor, epsilonV, 4]
    |>];
  Print[
    "      explicit=O(1)+O(v^2/c_gamma^2); next_uncomputed=O(v^4/c_gamma^4)"];

  section[
    "Build-global reference-only dependency and term coverage"];
  dependencyExpression =
    Total[Flatten[darwinKernel]] +
    electricU0 + routeAU2 + fullAnchor +
    Total[forceA2] + radialForceA2;
  If[
    activeMutation === "TARGET_BLINDNESS",
    dependencyExpression +=
      qT^2 (dV + aV)/(muR R)];
  forbiddenLive = Select[
    forbiddenTargetSymbols,
    !FreeQ[dependencyExpression, #] &];
  expectBool[
    "TARGET_BLINDNESS",
    forbiddenLive === {},
    <|"ForbiddenLive" -> forbiddenLive|>];
  Print[
    "      reference depends on cited A_E,c_gamma only; no Route-B/comparison/ratio/landing input"];

  inventory = computedTermInventory[
    activeMutation =!= "DUAL_ENGINE_TERMS"];
  expectBool[
    "DUAL_ENGINE_TERMS",
    inventory === expectedTermInventory,
    <|"Computed" -> inventory, "Required" -> expectedTermInventory|>];
  Print[
    "      TERM_INVENTORY=",
    StringRiffle[inventory, ","]];

  section[
    "Whole-stage restored-unit firewall on live expressions"];
  restoredResiduals = dimensionResiduals[
    activeMutation === "UNITS_RESTORED"];
  expectBool[
    "UNITS_RESTORED",
    And @@ (zeroQ /@ restoredResiduals),
    <|"Residuals" -> restoredResiduals|>];
  Print[
    "      [I_ij]=L^-1; [A_E]=M L^3 T^-2; [c_gamma^2]=L^2 T^-2"];
  Print[
    "      [D_V]=[A_V]=L^2 T^-2; [U_A]=M L^2 T^-2; [F_A2]=M L T^-2"];

  section[
    "Verdict re-derived from projector, anchor, and order objects"];
  verdictKernel = If[
    activeMutation === "VERDICT_REDERIVATION",
    2 darwinKernel,
    darwinKernel];
  verdictProjectorResiduals =
    unitSphereReduce /@
      Flatten[verdictKernel - expectedKernel];
  liveVerdict = deriveVerdict[
    verdictProjectorResiduals,
    unitSphereReduce[routeAU2 - expectedA2],
    computedVelocityOrder === 2];
  projectorNegative = deriveVerdict[
    unitSphereReduce /@
      Flatten[2 darwinKernel - expectedKernel],
    0,
    True];
  anchorNegative = deriveVerdict[
    kernelResiduals,
    unitSphereReduce[
      -s1 s2 aE dV/(8 Pi cGamma^2 R) -
        expectedA2],
    True];
  orderNegative = deriveVerdict[
    kernelResiduals,
    0,
    False];
  Print["      REDERIVED_TOKEN=", liveVerdict];
  expectBool[
    "VERDICT_REDERIVATION",
    liveVerdict === verdictToken &&
    projectorNegative === uncertifiedToken &&
    anchorNegative === uncertifiedToken &&
    orderNegative === uncertifiedToken,
    <|
      "Live" -> liveVerdict,
      "ProjectorNegative" -> projectorNegative,
      "AnchorNegative" -> anchorNegative,
      "OrderNegative" -> orderNegative
    |>];
  Print["      FAILURE_TOKEN=", uncertifiedToken];
  Print[
    "      OPEN_ELECTRIC_R1=", electricR1,
    " is a scope tag, not a certification failure"];

  section["Canonical source-to-stage predicate manifest"];
  manifestTest = If[
    activeMutation === "SOURCE_TO_STAGE_MANIFEST",
    Most[sourceManifest],
    sourceManifest];
  identifiers = manifestTest[[All, 1]];
  partitionCounts =
    KeySort@Counts[manifestTest[[All, 2]]];
  scopedOut = Cases[
    manifestTest,
    {identifier_, "SCOPED_OUT", _} :> identifier];
  inScope = Cases[
    manifestTest,
    {identifier_, disposition_, _} /;
      disposition =!= "SCOPED_OUT" :> identifier];
  manifestDigest = manifestSHA256[manifestTest];
  manifestOK =
    identifiers === sourceToothIDs &&
    Length[identifiers] ===
      Length[DeleteDuplicates[identifiers]] === 35 &&
    partitionCounts === expectedManifestCounts &&
    scopedOut === expectedScopedOut &&
    Length[scopedOut] === 29 &&
    inScope === expectedInScope &&
    manifestDigest === manifestDigestExpected &&
    AllTrue[
      manifestTest[[All, 3]],
      StringStartsQ[
        #,
        Alternatives[
          "STAGE034_", "STAGE035_", "STAGE036_",
          "STAGE037_", "STAGE038_"]] &];
  expectBool[
    "SOURCE_TO_STAGE_MANIFEST",
    manifestOK,
    <|
      "Entries" -> Length[manifestTest],
      "Partition" -> partitionCounts,
      "InScope" -> inScope,
      "ScopedOut" -> scopedOut,
      "Digest" -> manifestDigest
    |>];
  Print[
    "      entries=", Length[manifestTest],
    "; partition=", partitionCounts,
    "; scoped_out=", Length[scopedOut],
    "; digest=", manifestDigest];
  Print[
    "      IN_SCOPE=", StringRiffle[inScope, ","]];
  Print[
    "      SCOPED_OUT=",
    StringRiffle[scopedOut, ","]];

  Print[""];
  Print[
    "REFERENCE_ONLY=Route-A electric boost; Route-B/comparison/ratios deferred to stage037"];
  Print[
    "PROVENANCE=A_E carries R1_REQUIRED(bc_selection); c_gamma cited from stage003"];
  Print[
    "SCOPE=EARNED reference kernel at tier_A_conditional; electric sector not re-derived"];
  Print[
    "REMAINDER=O(v^4/c_gamma^4) named but not computed"];
  Print["VERDICT_TOKEN: ", liveVerdict];
  Print["PROGRESS MATHEMATICA STAGE036 COMPLETE"];

  If[activeMutation =!= "",
    Print["FIRST_FAILURE=MUTATION_DID_NOT_FIRE"];
    raise["MUTATION_DID_NOT_FIRE"]
  ];
  True,
  "ledgerStage036Failure",
  Function[{message, tag}, False]
];

Print[""];
Print["TOOTH_COUNT=", Length[toothOrder]];
Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
If[TrueQ[ok],
  Print[
    "OVERALL PASS: Mathematica independently reached ",
    verdictToken];
  Exit[0],
  Print[
    "OVERALL FAIL: Mathematica stage036 audit did not close"];
  Exit[1]
]
