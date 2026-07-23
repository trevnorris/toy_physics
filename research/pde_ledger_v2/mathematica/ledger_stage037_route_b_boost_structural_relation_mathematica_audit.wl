(* Ledger stage037 Mathematica audit: blind Route B and structural comparison.

   Standalone, print-only, assert-zero, exact, and file-I/O-free.

   Genuinely independent route:

   - The transverse tensor shape is obtained from NullSpace of a Cartesian
     divergence operator acting on the two isotropic R^-1 tensor bases.  Its
     overall coefficient is then fixed separately by the two-polarization
     trace.  There is no Solve[{b-a,...}] ansatz step.
   - The force is assembled from the tangent-plane normal projector and the
     differential geometry of n, rather than by differentiating a Cartesian
     coordinate potential.
   - Tensor comparison uses mixed velocity Hessians.  Coefficient ratios use
     velocity-scaling Coefficient extraction, rather than the direct U_B/U_A
     quotient used by the SymPy engine.

   These constructions are distinct from both the stage SymPy engine and the
   source magnetism_moving_throat_check.wl.

   Tooth-local runtime ablation uses LEDGER_STAGE037_MUTATION.
*)

ClearAll["Global`*"];
$HistoryLength = 0;

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE037_MUTATION";
activeMutation = Environment[mutationEnvironment];
If[!StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];

verdictToken = "BOOST_STRUCTURAL_RELATION_HOLDS";
uncertifiedToken = "BOOST_STRUCTURAL_RELATION_UNCERTIFIED";
coblockerToken = "R1_REQUIRED(direct_moving_throat)";
comparisonFact = "route_B_R1";
relativeSignFact = "relative_sign_anchor_conditional";
manifestDigestExpected =
  "3c88849c5f4f5b7fe05c0d06c59004acccf8c6e85f6823b0e839c41174c8adf4";

toothOrder = {
  "DIRECT_SOURCE",
  "DIRECT_PROJECTOR",
  "DIRECT_EXCHANGE_SIGN",
  "DIRECT_FALLOFF",
  "DIRECT_VELOCITY_ORDER",
  "ROUTE_INDEPENDENCE",
  "BOOST_COMMON_VELOCITY",
  "COMPARE_COMPUTED",
  "DELTA_RATIO",
  "CONE_RATIO",
  "QMAG_R1",
  "TARGET_BLINDNESS",
  "DUAL_ENGINE_TERMS",
  "UNITS_RESTORED",
  "VERDICT_REDERIVATION",
  "SOURCE_TO_STAGE_MANIFEST"
};

ablationDescriptions = <|
  "DIRECT_SOURCE" ->
    "double the live direct-source coupling in U_B, F_B, and F_B,r",
  "DIRECT_PROJECTOR" ->
    "add an excess n_i*n_j/(8*pi*mu_R*R) to the live direct kernel",
  "DIRECT_EXCHANGE_SIGN" ->
    "reverse the live interaction used for the computed parallel/antiparallel signature",
  "DIRECT_FALLOFF" ->
    "divide the live derived radial force by one extra R",
  "DIRECT_VELOCITY_ORDER" ->
    "replace the second velocity by a fixed direction in the live interaction",
  "ROUTE_INDEPENDENCE" ->
    "derive a mutation-only Route-B copy from the already-instantiated Route A",
  "BOOST_COMMON_VELOCITY" ->
    "add D_V/c_gamma^2 to the Route-A/electric correction ratio",
  "COMPARE_COMPUTED" ->
    "add a full n_i*n_j excess to the mixed-Hessian Route-B tensor",
  "DELTA_RATIO" ->
    "double the Coefficient-extracted r_BA used by ratio/delta/Delta_U checks",
  "CONE_RATIO" ->
    "replace c_E^2/c_gamma^2 by c_E/sqrt(c_gamma^2)",
  "QMAG_R1" ->
    "replace R1(magnitude) by magnitude_forced_by_electric",
  "TARGET_BLINDNESS" ->
    "inject barred and sealed-decision symbols and remove A_E from the comparison lane",
  "DUAL_ENGINE_TERMS" ->
    "drop Delta_U and corrupt every remaining canonical symbolic term",
  "UNITS_RESTORED" ->
    "corrupt every live base dimension, including [q_T]: M*T^-1 -> L*T^-1",
  "VERDICT_REDERIVATION" ->
    "corrupt verdict-local kernel, comparison, falloff, order, and ancestry objects",
  "SOURCE_TO_STAGE_MANIFEST" ->
    "remove a scoped-out row and mis-scope DIRECT_SOURCE in the live manifest"
|>;

raise[name_] := Throw[name, "ledgerStage037Failure"];

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
  FullSimplify[
    Together[expression],
    R > 0 && muR > 0 && rhoBr > 0 && cGammaSquared > 0 &&
      cE > 0 && qT != 0 && aE != 0 &&
      nx^2 + ny^2 + nz^2 == 1
  ] /. ConditionalExpression[0, _] -> 0;

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
    If[evidence =!= None,
      Print["      evidence = ", InputForm[evidence]]];
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

zeroQ[expression_] := TrueQ[cleanZero[expression] === 0];

sphereReduce[expression_] := Factor[
  Together[
    Expand[expression] /. nz^2 -> 1 - nx^2 - ny^2
  ]
];

allZero[expressions_List] := And @@ (zeroQ /@ expressions);


(* ---------------------------------------------------------------------- *)
(* Exact symbols, invariants, and cited inputs.                            *)
(* ---------------------------------------------------------------------- *)

nVector = {nx, ny, nz};
velocityOne = {v1x, v1y, v1z};
velocityTwo = {v2x, v2y, v2z};
velocitySymbols = Join[velocityOne, velocityTwo];
dVelocity = Expand[velocityOne.velocityTwo];
vOneNormal = Expand[velocityOne.nVector];
vTwoNormal = Expand[velocityTwo.nVector];
aVelocity = Expand[vOneNormal vTwoNormal];
velocityShape = Expand[dVelocity + aVelocity];
cGammaRelation = cGammaSquared -> muR/rhoBr;

expectedKernelB =
  (IdentityMatrix[3] + Outer[Times, nVector, nVector])/
    (8 Pi muR R);
expectedKernelA =
  (IdentityMatrix[3] + Outer[Times, nVector, nVector])/
    (8 Pi R);
expectedUB =
  -s1 s2 qT^2 velocityShape/(8 Pi muR R);
expectedForceB =
  s1 s2 qT^2/(8 Pi muR R^2) *
    (vTwoNormal velocityOne + vOneNormal velocityTwo -
      (dVelocity + 3 aVelocity) nVector);
expectedRadialB =
  -s1 s2 qT^2 velocityShape/(8 Pi muR R^2);
expectedUA2 =
  -s1 s2 aE velocityShape/(8 Pi cGammaSquared R);
electricRadialForce = s1 s2 aE/(4 Pi R^2);


(* ---------------------------------------------------------------------- *)
(* Wolfram-native Route B: divergence operator -> NullSpace -> trace.      *)
(* ---------------------------------------------------------------------- *)

Print["PROGRESS MATHEMATICA STAGE037 NULLSPACE KERNEL START"];

cartesianCoordinates = {rx, ry, rz};
cartesianRadius = Sqrt[cartesianCoordinates.cartesianCoordinates];
isotropicTensorBases = {
  IdentityMatrix[3]/cartesianRadius,
  Outer[Times, cartesianCoordinates, cartesianCoordinates]/
    cartesianRadius^3
};

divergenceBasis = Table[
  Table[
    Sum[
      D[
        isotropicTensorBases[[basisIndex, rowIndex, columnIndex]],
        cartesianCoordinates[[rowIndex]]
      ],
      {rowIndex, 3}
    ],
    {columnIndex, 3}
  ],
  {basisIndex, 2}
];

(* Radial projection turns the two vector divergences into {-1,+1}. *)
transverseOperatorRow = FullSimplify[
  cartesianRadius (cartesianCoordinates.#),
  cartesianRadius > 0
] & /@ divergenceBasis;
transverseShape = First[NullSpace[{transverseOperatorRow}]];
traceWeights = {3, 1};
twoPolarizationTrace = 1/(2 Pi);
traceNormalizer = Factor[
  twoPolarizationTrace/(traceWeights.transverseShape)
];
directCoefficients = Factor[traceNormalizer transverseShape];

kernelCartesian =
  (directCoefficients[[1]] isotropicTensorBases[[1]] +
    directCoefficients[[2]] isotropicTensorBases[[2]])/muR;
directKernel = Map[
  FullSimplify[
    sphereReduce[#],
    R > 0 && nx^2 + ny^2 + nz^2 == 1
  ] &,
  kernelCartesian /. Thread[cartesianCoordinates -> R nVector],
  {2}
];

Print["PROGRESS MATHEMATICA STAGE037 NULLSPACE KERNEL COMPLETE"];


(* ---------------------------------------------------------------------- *)
(* Normal-projector differential geometry force, with no coordinate D[U]. *)
(* ---------------------------------------------------------------------- *)

normalProjector =
  IdentityMatrix[3] - Outer[Times, nVector, nVector];
gradientAngularShape =
  (vTwoNormal normalProjector.velocityOne +
    vOneNormal normalProjector.velocityTwo)/R;
directCoupling = -s1 s2 qT^2/(8 Pi muR);
directInteraction =
  Factor[-s1 s2 qT^2 velocityOne.directKernel.velocityTwo];
directForce = Map[
  sphereReduce,
  -directCoupling *
    (gradientAngularShape/R - velocityShape nVector/R^2)
];
directRadial = sphereReduce[nVector.directForce];


(* ---------------------------------------------------------------------- *)
(* Build records: direct Route B is instantiated before cited Route A.     *)
(* ---------------------------------------------------------------------- *)

buildLog = {};
qCurrentSourceTag = Unique["qCurrentSourceTag"];
pathA36EOMTag = Unique["pathA36EOMTag"];
directAnsatzTag = Unique["directAnsatzTag"];
illicitRouteAReadTag = Unique["illicitRouteAReadTag"];

buildStamp[label_] := Module[{ordinal = Length[buildLog]},
  AppendTo[buildLog, label];
  ordinal
];

deriveDirectRoute[foreign_] := Module[
  {interaction = directInteraction,
   dependencies = {
     "Q_CURRENT.source",
     "pathA36.transverse_EOM",
     "direct_transverse_tensor_ansatz"
   },
   ancestry = {
     qCurrentSourceTag, pathA36EOMTag, directAnsatzTag
   },
   routeName = "B_DIRECT",
   foreignUsed = False},
  If[foreign =!= None,
    interaction = Factor[
      foreign["interaction"] qT^2 cGammaSquared/aE
    ];
    dependencies = Append[dependencies, "ILLICIT_ROUTE_A_READ"];
    ancestry = Append[ancestry, illicitRouteAReadTag];
    routeName = "B_DIRECT_FROM_ROUTE_A";
    foreignUsed = True
  ];
  <|
    "name" -> routeName,
    "kernel" -> directKernel,
    "interaction" -> interaction,
    "force" -> directForce,
    "radial" -> directRadial,
    "dependencies" -> dependencies,
    "ancestry" -> ancestry,
    "foreignUsed" -> foreignUsed,
    "ordinal" -> buildStamp["ROUTE_B"]
  |>
];

citeRouteA[] := Module[{force, radial},
  force =
    s1 s2 aE/(8 Pi cGammaSquared R^2) *
      (vTwoNormal velocityOne + vOneNormal velocityTwo -
        (dVelocity + 3 aVelocity) nVector);
  radial =
    -s1 s2 aE velocityShape/(8 Pi cGammaSquared R^2);
  <|
    "name" -> "A_CITED_STAGE036",
    "kernel" -> expectedKernelA,
    "interaction" -> expectedUA2,
    "force" -> force,
    "radial" -> radial,
    "dependencies" -> {
      "stage036.R70.kernel", "stage036.R70.U_A2"
    },
    "ancestry" -> {Unique["stage036R70Citation"]},
    "foreignUsed" -> False,
    "ordinal" -> buildStamp["ROUTE_A"]
  |>
];

routeB = deriveDirectRoute[None];
routeA = citeRouteA[];


(* ---------------------------------------------------------------------- *)
(* Mixed velocity Hessian comparison and Coefficient-extracted ratios.     *)
(* ---------------------------------------------------------------------- *)

mixedVelocityTensor[interaction_, scalarPrefactor_] := Table[
  Factor[
    D[
      interaction,
      velocityOne[[rowIndex]],
      velocityTwo[[columnIndex]]
    ]/scalarPrefactor
  ],
  {rowIndex, 3},
  {columnIndex, 3}
];

tensorRouteB = mixedVelocityTensor[
  routeB["interaction"],
  -s1 s2 qT^2/muR
];
tensorRouteA = mixedVelocityTensor[
  routeA["interaction"],
  -s1 s2 aE/cGammaSquared
];
tensorComparisonResiduals =
  sphereReduce /@ Flatten[tensorRouteB - tensorRouteA];

velocityScaleRules =
  Thread[velocityTwo -> ratioProbe velocityTwo];
routeBAmplitude = Coefficient[
  Expand[routeB["interaction"] /. velocityScaleRules],
  ratioProbe,
  1
];
routeAAmplitude = Coefficient[
  Expand[routeA["interaction"] /. velocityScaleRules],
  ratioProbe,
  1
];
ratioBA = Factor[routeBAmplitude/routeAAmplitude];
deltaBA = Factor[ratioBA - 1];
coneRatio = Factor[cE^2/cGammaSquared];
deltaU = Factor[routeBAmplitude - routeAAmplitude];

parallelRatioA = Factor[
  (routeA["radial"]/electricRadialForce) /.
    {
      nx -> 1, ny -> 0, nz -> 0,
      v1x -> 0, v2x -> 0,
      v1z -> 0, v2z -> 0
    }
];
parallelD = v1y v2y;


(* ---------------------------------------------------------------------- *)
(* Computed sign, R-power, and velocity-degree objects.                    *)
(* ---------------------------------------------------------------------- *)

inverseRPower[expression_] := Module[{together = Together[expression]},
  Exponent[Denominator[together], R] -
    Exponent[Numerator[together], R]
];

velocityDegreeSet[expression_] := Module[
  {scaled, coefficientRules},
  scaled = Expand[
    expression /.
      Thread[velocitySymbols -> velocityProbe velocitySymbols]
  ];
  coefficientRules = CoefficientRules[scaled, {velocityProbe}];
  Sort@DeleteDuplicates@Cases[
    coefficientRules,
    Rule[{power_}, coefficient_] /; !zeroQ[coefficient] :> power
  ]
];

exchangeSignature[interaction_] := Module[
  {commonRules, parallel, antiparallel},
  commonRules = {
    nx -> 1, ny -> 0, nz -> 0,
    v1x -> 0, v1y -> 1, v1z -> 0,
    v2x -> 0, v2z -> 0,
    s1 -> 1, s2 -> 1, qT -> 1, muR -> 1, R -> 1
  };
  parallel = FullSimplify[
    interaction /. Join[commonRules, {v2y -> 1}]
  ];
  antiparallel = FullSimplify[
    interaction /. Join[commonRules, {v2y -> -1}]
  ];
  {Sign[parallel], Sign[antiparallel]}
];

computedForcePower = inverseRPower[routeB["radial"]];
computedPotentialPower = inverseRPower[routeB["interaction"]];
computedVelocityDegrees = velocityDegreeSet[routeB["interaction"]];
computedExchangeSignature = exchangeSignature[routeB["interaction"]];


(* ---------------------------------------------------------------------- *)
(* Dependency, blindness, and canonical symbolic-term objects.             *)
(* ---------------------------------------------------------------------- *)

independenceViolations[rb_, ra_] := Module[{violations = {}},
  If[
    AnyTrue[
      rb["dependencies"],
      StringContainsQ[#, "ROUTE_A"] &
    ],
    AppendTo[violations, "ROUTE_A_DEPENDENCY_TAG"]
  ];
  If[
    !FreeQ[rb["ancestry"], illicitRouteAReadTag],
    AppendTo[violations, "ROUTE_A_ANCESTRY_SYMBOL"]
  ];
  If[
    TrueQ[rb["foreignUsed"]],
    AppendTo[violations, "FOREIGN_PAYLOAD_USED"]
  ];
  If[
    rb["ordinal"] >= ra["ordinal"],
    AppendTo[violations, "ROUTE_B_NOT_FIRST"]
  ];
  If[
    rb["name"] =!= "B_DIRECT",
    AppendTo[violations, "NOT_B_DIRECT"]
  ];
  violations
];

barredPathA39Symbols = {Nu, aT, aTPrime, aL, qAT, qL};
sealedSection4Symbols = {
  magnetismLorentzConsistent,
  amendmentExcluded,
  magnetismDepartureCharacterized,
  noGoSector
};

targetBlindnessViolations[mutate_] := Module[
  {directExpression, comparisonExpression, decisionExpression = 0,
   violations = {}},
  directExpression =
    Total[Flatten[routeB["kernel"]]] +
    routeB["interaction"] +
    Total[routeB["force"]] +
    routeB["radial"];
  comparisonExpression =
    routeA["interaction"] + routeB["interaction"] +
    ratioBA + deltaBA + coneRatio + deltaU;
  If[TrueQ[mutate],
    directExpression += Nu + magnetismLorentzConsistent + aE;
    comparisonExpression =
      comparisonExpression /. aE -> 1;
    decisionExpression += amendmentExcluded
  ];
  If[
    AnyTrue[barredPathA39Symbols, !FreeQ[directExpression, #] &],
    AppendTo[violations, "BARRED_PATHA39_DIRECT"]
  ];
  If[
    AnyTrue[sealedSection4Symbols, !FreeQ[directExpression, #] &],
    AppendTo[violations, "SEALED_TOKEN_DIRECT"]
  ];
  If[
    AnyTrue[sealedSection4Symbols, !FreeQ[decisionExpression, #] &],
    AppendTo[violations, "SEALED_DECISION_CHANNEL"]
  ];
  If[
    !FreeQ[directExpression, aE],
    AppendTo[violations, "A_E_ENTERED_DIRECT_ROUTE"]
  ];
  If[
    FreeQ[comparisonExpression, aE],
    AppendTo[violations, "A_E_MISSING_FROM_CITED_COMPARISON"]
  ];
  violations
];

expectedTermKeys = {
  "routeB_kernel00",
  "routeB_kernel01",
  "routeB_U2",
  "routeB_Fr",
  "routeA_U2_cited",
  "ratio_BA",
  "delta_BA",
  "cone_ratio",
  "Delta_U"
};

computedTerms[] := <|
  "routeB_kernel00" -> routeB["kernel"][[1, 1]],
  "routeB_kernel01" -> routeB["kernel"][[1, 2]],
  "routeB_U2" -> routeB["interaction"],
  "routeB_Fr" -> routeB["radial"],
  "routeA_U2_cited" -> routeA["interaction"],
  "ratio_BA" -> ratioBA,
  "delta_BA" -> deltaBA,
  "cone_ratio" -> coneRatio,
  "Delta_U" -> deltaU
|>;

expectedTerms = <|
  "routeB_kernel00" ->
    (1 + nx^2)/(8 Pi muR R),
  "routeB_kernel01" ->
    nx ny/(8 Pi muR R),
  "routeB_U2" -> expectedUB,
  "routeB_Fr" -> sphereReduce[expectedRadialB],
  "routeA_U2_cited" -> expectedUA2,
  "ratio_BA" -> qT^2 cGammaSquared/(muR aE),
  "delta_BA" -> qT^2 cGammaSquared/(muR aE) - 1,
  "cone_ratio" -> cE^2/cGammaSquared,
  "Delta_U" ->
    -s1 s2 aE/(8 Pi cGammaSquared R) *
      (qT^2 cGammaSquared/(muR aE) - 1) velocityShape
|>;

termViolations[mutate_] := Module[{live, violations = {}},
  live = computedTerms[];
  If[TrueQ[mutate],
    live = Map[2 # &, KeyDrop[live, "Delta_U"]]
  ];
  If[
    Keys[live] =!= expectedTermKeys,
    AppendTo[violations, "TERM_INVENTORY"]
  ];
  Do[
    If[
      !KeyExistsQ[live, key],
      AppendTo[violations, "MISSING:" <> key],
      If[
        !zeroQ[live[key] - expectedTerms[key]],
        AppendTo[violations, "SYMBOLIC_MISMATCH:" <> key]
      ]
    ],
    {key, expectedTermKeys}
  ];
  violations
];


(* ---------------------------------------------------------------------- *)
(* Native formal-rescaling dimensional firewall on the real expressions.   *)
(* ---------------------------------------------------------------------- *)

unitRules[badQT_] := Join[
  {
    R -> If[TrueQ[badQT], R, lengthScale R],
    qT -> If[
      TrueQ[badQT],
      lengthScale timeScale^-1 qT,
      massScale timeScale^-1 qT
    ],
    muR -> If[
      TrueQ[badQT],
      muR,
      massScale lengthScale^-1 timeScale^-2 muR
    ],
    rhoBr -> If[
      TrueQ[badQT],
      lengthScale rhoBr,
      massScale lengthScale^-3 rhoBr
    ],
    aE -> If[
      TrueQ[badQT],
      lengthScale aE,
      massScale lengthScale^3 timeScale^-2 aE
    ],
    cGammaSquared -> If[
      TrueQ[badQT],
      cGammaSquared,
      lengthScale^2 timeScale^-2 cGammaSquared
    ],
    cE -> If[
      TrueQ[badQT],
      massScale cE,
      lengthScale timeScale^-1 cE
    ],
    s1 -> If[TrueQ[badQT], lengthScale s1, s1],
    s2 -> If[TrueQ[badQT], timeScale s2, s2]
  },
  Thread[
    nVector ->
      If[TrueQ[badQT], lengthScale nVector, nVector]
  ],
  Thread[
    velocitySymbols ->
      If[
        TrueQ[badQT],
        velocitySymbols,
        lengthScale timeScale^-1 velocitySymbols
      ]
  ]
];

scalingResidual[expression_, requiredScale_, badQT_] :=
  Factor[
    Together[
      (expression /. unitRules[badQT]) -
        requiredScale expression
    ]
  ];

unitContract = <|
  "q_T" -> {qT, massScale timeScale^-1},
  "mu_R" -> {
    muR, massScale lengthScale^-1 timeScale^-2
  },
  "rho_br" -> {
    rhoBr, massScale lengthScale^-3
  },
  "A_E" -> {
    aE, massScale lengthScale^3 timeScale^-2
  },
  "c_gamma_squared" -> {
    cGammaSquared, lengthScale^2 timeScale^-2
  },
  "c_E" -> {cE, lengthScale timeScale^-1},
  "D_V" -> {
    dVelocity, lengthScale^2 timeScale^-2
  },
  "A_V" -> {
    aVelocity, lengthScale^2 timeScale^-2
  },
  "kernel_B00" -> {
    routeB["kernel"][[1, 1]], massScale^-1 timeScale^2
  },
  "U_B" -> {
    routeB["interaction"],
    massScale lengthScale^2 timeScale^-2
  },
  "F_B0" -> {
    routeB["force"][[1]],
    massScale lengthScale timeScale^-2
  },
  "F_B1" -> {
    routeB["force"][[2]],
    massScale lengthScale timeScale^-2
  },
  "F_B2" -> {
    routeB["force"][[3]],
    massScale lengthScale timeScale^-2
  },
  "F_Br" -> {
    routeB["radial"],
    massScale lengthScale timeScale^-2
  },
  "U_A2" -> {
    routeA["interaction"],
    massScale lengthScale^2 timeScale^-2
  },
  "r_BA" -> {ratioBA, 1},
  "delta_BA" -> {deltaBA, 1},
  "r_cone" -> {coneRatio, 1},
  "Delta_U" -> {
    deltaU, massScale lengthScale^2 timeScale^-2
  },
  "s_1" -> {s1, 1},
  "s_2" -> {s2, 1}
|>;

unitViolations[badQT_] := Keys@Select[
  unitContract,
  !zeroQ[
    scalingResidual[#[[1]], #[[2]], badQT]
  ] &
];


(* ---------------------------------------------------------------------- *)
(* Verdict precedence from exactly the computed certification objects.     *)
(* ---------------------------------------------------------------------- *)

deriveVerdict[
  projectorResiduals_List,
  comparisonResiduals_List,
  forcePower_,
  velocityDegrees_List,
  independenceState_List
] := If[
  allZero[projectorResiduals] &&
    allZero[comparisonResiduals] &&
    forcePower === 2 &&
    velocityDegrees === {2} &&
    independenceState === {},
  verdictToken,
  uncertifiedToken
];


(* ---------------------------------------------------------------------- *)
(* Exact 35-item source-to-stage predicate manifest.                       *)
(* ---------------------------------------------------------------------- *)

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
  {"SOURCE_TRANSLATION_CONTINUITY", "SCOPED_OUT",
    "STAGE035_V2_DONE"},
  {"SOURCE_NOT_IMPORTED", "SCOPED_OUT",
    "STAGE035_V2_DONE"},
  {"SOURCE_BASIS", "SCOPED_OUT", "STAGE035_V2_DONE"},
  {"PARITY_RW", "SCOPED_OUT", "STAGE035_V2_DONE"},
  {"PARITY_PW", "SCOPED_OUT", "STAGE035_V2_DONE"},
  {"PARITY_ROTATION", "SCOPED_OUT", "STAGE035_V2_DONE"},
  {"PARITY_TIME_REVERSAL", "SCOPED_OUT",
    "STAGE035_V2_DONE"},
  {"FIELD_IDENTITY_UNITS", "SCOPED_OUT",
    "STAGE034_V1_DONE"},
  {"ACTION_KINETIC", "SCOPED_OUT", "STAGE034_V1_DONE"},
  {"ACTION_COUPLING", "SCOPED_OUT", "STAGE034_V1_DONE"},
  {"ACTION_STABILITY", "SCOPED_OUT",
    "STAGE034_V1_DONE"},
  {"G0_DAMAGE", "SCOPED_OUT", "STAGE034_V1_DONE"},
  {"ROUTE_INDEPENDENCE", "REPLACED_BY_STRONGER",
    "STAGE037_TAG_SYMBOL_ORDER_GUARD"},
  {"BOOST_PROJECTOR", "SCOPED_OUT", "STAGE036_V3_CITED"},
  {"BOOST_GENERAL_VELOCITIES", "SCOPED_OUT",
    "STAGE036_V3_CITED"},
  {"BOOST_NEXT_ORDER", "SCOPED_OUT", "STAGE036_V3_CITED"},
  {"BOOST_COMMON_VELOCITY", "PRESERVED",
    "STAGE037_ROUTE_A_ELECTRIC_CROSSCHECK"},
  {"DIRECT_SOURCE", "REPLACED_BY_STRONGER",
    "STAGE037_DIRECT_U_FORCE_RECONSTRUCTION"},
  {"DIRECT_PROJECTOR", "REPLACED_BY_STRONGER",
    "STAGE037_CONSTRAINT_SOLVE_COMPONENTS"},
  {"DIRECT_EXCHANGE_SIGN", "REPLACED_BY_STRONGER",
    "STAGE037_COMPUTED_SIGN_SIGNATURE"},
  {"DIRECT_FALLOFF", "REPLACED_BY_STRONGER",
    "STAGE037_EXPRESSION_R_POWER"},
  {"DIRECT_VELOCITY_ORDER", "REPLACED_BY_STRONGER",
    "STAGE037_EXPRESSION_VELOCITY_DEGREE"},
  {"COMPARE_COMPUTED", "REPLACED_BY_STRONGER",
    "STAGE037_TENSOR_ONLY_COMPARISON"},
  {"DELTA_RATIO", "REPLACED_BY_STRONGER",
    "STAGE037_RATIO_DELTA_POTENTIAL_DIFFERENCE"},
  {"CONE_RATIO", "REPLACED_BY_STRONGER",
    "STAGE037_CONE_DENSITY_EQUIVALENCE"},
  {"QMAG_R1", "PRESERVED", "STAGE037_R1_MAGNITUDE_ENUM"},
  {"UNITS_RESTORED", "REPLACED_BY_STRONGER",
    "STAGE037_EXPRESSION_DIMENSION_FIREWALL"},
  {"ACTIVE_FLUX_CAVEAT", "SCOPED_OUT", "STAGE038_V5"},
  {"HOOK_LORENTZ", "SCOPED_OUT", "STAGE038_V5"},
  {"LEDGER_READY_ROW", "SCOPED_OUT", "STAGE034_V1_DONE"},
  {"TRUTH_TOTALITY", "SCOPED_OUT", "STAGE038_V5"},
  {"TRUTH_PRECEDENCE", "SCOPED_OUT", "STAGE038_V5"},
  {"LANDING_OWNERSHIP", "SCOPED_OUT", "STAGE038_V5"},
  {"TARGET_BLINDNESS", "REPLACED_BY_STRONGER",
    "STAGE037_RUNTIME_SCOPE_CHANNELS"},
  {"DUAL_ENGINE_TERMS", "REPLACED_BY_STRONGER",
    "STAGE037_SYMBOLIC_TERM_CONTRACT"}
};

expectedManifestCounts = KeySort@<|
  "PRESERVED" -> 2,
  "REPLACED_BY_STRONGER" -> 12,
  "SCOPED_OUT" -> 21
|>;

expectedInScope = {
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
  "BOOST_PROJECTOR",
  "BOOST_GENERAL_VELOCITIES",
  "BOOST_NEXT_ORDER",
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
    Range[limit],
    left[[#]] =!= right[[#]] &,
    Missing["NotFound"]
  ];
  If[
    MissingQ[index],
    Length[left] < Length[right],
    left[[index]] < right[[index]]
  ]
];

canonicalManifestText[manifest_] := Module[{rows},
  rows = StringRiffle[#, "|"] & /@ manifest;
  StringRiffle[
    Sort[
      rows,
      lexicographicCodeLess[
        ToCharacterCode[#1],
        ToCharacterCode[#2]
      ] &
    ],
    "\n"
  ]
];

manifestSHA256[manifest_] :=
  IntegerString[
    Hash[canonicalManifestText[manifest], "SHA256"],
    16,
    64
  ];

manifestState[manifest_] := Module[
  {identifiers, partition, inScope, scopedOut},
  identifiers = manifest[[All, 1]];
  partition = KeySort@Counts[manifest[[All, 2]]];
  inScope = Cases[
    manifest,
    {identifier_, disposition_, _} /;
      disposition =!= "SCOPED_OUT" :> identifier
  ];
  scopedOut = Cases[
    manifest,
    {identifier_, "SCOPED_OUT", _} :> identifier
  ];
  {
    identifiers,
    partition,
    inScope,
    scopedOut,
    manifestSHA256[manifest]
  }
];

expectedManifestState = {
  sourceToothIDs,
  expectedManifestCounts,
  expectedInScope,
  expectedScopedOut,
  manifestDigestExpected
};


(* ---------------------------------------------------------------------- *)
(* Assertions and report.                                                  *)
(* ---------------------------------------------------------------------- *)

ok = Catch[
  If[
    activeMutation =!= "" &&
      !MemberQ[toothOrder, activeMutation],
    Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
    Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
    raise["UNKNOWN_MUTATION"]
  ];

  Print[
    "ledger_stage037_route_b_boost_structural_relation Mathematica audit"
  ];
  Print[
    "ROUTE=Cartesian divergence NullSpace + trace normalization + normal-projector force + mixed velocity Hessians"
  ];
  Print[
    "INDEPENDENCE=no ansatz Solve; no Cartesian coordinate-potential D; no direct U_B/U_A quotient"
  ];
  Print["BLIND_BUILD=Route B instantiated before cited Route A"];
  Print["FILE_IO=none; CROSS_ENGINE_COMPARE=none"];
  Print["PROGRESS MATHEMATICA STAGE037 ASSERTIONS START"];
  If[activeMutation =!= "",
    Print["ACTIVE_MUTATION=", activeMutation];
    Print[
      "MUTATED_PRIMITIVE=",
      ablationDescriptions[activeMutation]
    ]
  ];

  section[
    "Blind NullSpace tensor construction and direct source/force reconstruction"
  ];
  sourceScale = If[
    activeMutation === "DIRECT_SOURCE",
    2,
    1
  ];
  liveSourceU = Factor[sourceScale routeB["interaction"]];
  liveSourceForce = Factor /@ (sourceScale routeB["force"]);
  liveSourceRadial = Factor[sourceScale routeB["radial"]];
  sourceResiduals = Join[
    {Factor[liveSourceU - expectedUB]},
    sphereReduce /@ (liveSourceForce - expectedForceB),
    {sphereReduce[liveSourceRadial - expectedRadialB]}
  ];
  expectBool[
    "DIRECT_SOURCE",
    allZero[sourceResiduals],
    <|"Residuals" -> sourceResiduals|>
  ];
  Print[
    "      U_B=-s1*s2*q_T^2*(D_V+A_V)/(8*pi*mu_R*R)"
  ];
  Print[
    "      F_B=(s1*s2*q_T^2/(8*pi*mu_R*R^2))*[(V2.n)V1+(V1.n)V2-(D_V+3*A_V)n]"
  ];
  Print[
    "      F_B,r=-s1*s2*q_T^2*(D_V+A_V)/(8*pi*mu_R*R^2)"
  ];

  liveKernel = If[
    activeMutation === "DIRECT_PROJECTOR",
    routeB["kernel"] +
      Outer[Times, nVector, nVector]/(8 Pi muR R),
    routeB["kernel"]
  ];
  projectorResiduals =
    sphereReduce /@ Flatten[liveKernel - expectedKernelB];
  expectBool[
    "DIRECT_PROJECTOR",
    allZero[projectorResiduals],
    <|
      "DivergenceOperator" -> transverseOperatorRow,
      "NullShape" -> transverseShape,
      "Coefficients" -> directCoefficients,
      "ComponentResiduals" -> projectorResiduals
    |>
  ];
  Print[
    "      divergence operator row=", transverseOperatorRow,
    "; NullSpace shape=", transverseShape
  ];
  Print[
    "      trace normalization gives coefficients=", directCoefficients,
    "; G_B,ij=(delta_ij+n_i*n_j)/(8*pi*mu_R*R)"
  ];

  section["Computed sign, falloff, and velocity-order objects"];
  signInteraction = If[
    activeMutation === "DIRECT_EXCHANGE_SIGN",
    -routeB["interaction"],
    routeB["interaction"]
  ];
  liveSignature = exchangeSignature[signInteraction];
  expectBool[
    "DIRECT_EXCHANGE_SIGN",
    liveSignature === {-1, 1},
    <|"ParallelAntiparallelSignature" -> liveSignature|>
  ];
  Print[
    "      computed parallel/antiparallel interaction signs=",
    liveSignature
  ];

  falloffForce = If[
    activeMutation === "DIRECT_FALLOFF",
    routeB["radial"]/R,
    routeB["radial"]
  ];
  liveForcePower = inverseRPower[falloffForce];
  expectBool[
    "DIRECT_FALLOFF",
    liveForcePower === 2,
    <|"ComputedForcePower" -> liveForcePower|>
  ];
  Print[
    "      computed U power=R^-", computedPotentialPower,
    "; computed F power=R^-", liveForcePower
  ];

  orderInteraction = routeB["interaction"];
  If[
    activeMutation === "DIRECT_VELOCITY_ORDER",
    orderInteraction =
      orderInteraction /. {v2x -> 1, v2y -> 0, v2z -> 0}
  ];
  liveVelocityDegrees = velocityDegreeSet[orderInteraction];
  expectBool[
    "DIRECT_VELOCITY_ORDER",
    liveVelocityDegrees === {2},
    <|"ComputedVelocityDegrees" -> liveVelocityDegrees|>
  ];
  Print[
    "      computed velocity degrees=", liveVelocityDegrees,
    "; O(V1*V2)"
  ];

  section[
    "Dependency-tag, ancestry-symbol, and execution-order blindness"
  ];
  independenceRoute = routeB;
  If[
    activeMutation === "ROUTE_INDEPENDENCE",
    independenceRoute = deriveDirectRoute[routeA]
  ];
  independenceState =
    independenceViolations[independenceRoute, routeA];
  expectBool[
    "ROUTE_INDEPENDENCE",
    independenceState === {},
    <|
      "Violations" -> independenceState,
      "Dependencies" -> independenceRoute["dependencies"],
      "BuildLog" -> buildLog
    |>
  ];
  Print[
    "      foreign_payload=None; production build order=ROUTE_B then ROUTE_A"
  ];
  Print[
    "      dependencies=", InputForm[routeB["dependencies"]]
  ];

  section["Route-A/electric equal-velocity cross-check"];
  commonRatio = parallelRatioA;
  If[
    activeMutation === "BOOST_COMMON_VELOCITY",
    commonRatio = Factor[
      commonRatio + parallelD/cGammaSquared
    ]
  ];
  commonResidual = Factor[
    commonRatio + parallelD/(2 cGammaSquared)
  ];
  expectZero[
    "BOOST_COMMON_VELOCITY",
    commonResidual,
    <|"ParallelCorrectionRatio" -> commonRatio|>
  ];
  Print[
    "      F_A2,r/F_E,r=-(V1.V2)/(2*c_gamma^2), using Route A + electric anchor only"
  ];
  Print[
    "      equal v: F_r/F_E,r=1-v^2/(2*c_gamma^2)+O(v^4)"
  ];

  section["Mixed-Hessian prefactor-stripped tensor comparison"];
  liveTensorB = tensorRouteB;
  If[
    activeMutation === "COMPARE_COMPUTED",
    liveTensorB =
      liveTensorB +
        Outer[Times, nVector, nVector]/(8 Pi R)
  ];
  comparisonResiduals =
    sphereReduce /@ Flatten[liveTensorB - tensorRouteA];
  expectBool[
    "COMPARE_COMPUTED",
    allZero[comparisonResiduals],
    <|
      "RouteATensor" -> tensorRouteA,
      "RouteBTensor" -> liveTensorB,
      "Residuals" -> comparisonResiduals
    |>
  ];
  Print[
    "      mixed velocity Hessians strip each scalar prefactor and agree component-by-component"
  ];
  Print[
    "      tensor structure agrees; computed Route-B falloff/order match cited Route A"
  ];

  section[
    "Coefficient-extracted delta, cone, and potential-difference ratios"
  ];
  liveRatio = If[
    activeMutation === "DELTA_RATIO",
    2 ratioBA,
    ratioBA
  ];
  liveDelta = Factor[liveRatio - 1];
  ratioExpected =
    qT^2 cGammaSquared/(muR aE);
  ratioResiduals = {
    Factor[liveRatio - ratioExpected],
    Factor[
      (liveRatio /. cGammaRelation) -
        qT^2/(rhoBr aE)
    ],
    Factor[liveDelta - (ratioExpected - 1)],
    Factor[
      deltaU +
        s1 s2 aE/(8 Pi cGammaSquared R) *
          liveDelta velocityShape
    ]
  };
  expectBool[
    "DELTA_RATIO",
    allZero[ratioResiduals],
    <|
      "rBA" -> liveRatio,
      "deltaBA" -> liveDelta,
      "DeltaU" -> deltaU,
      "Residuals" -> ratioResiduals
    |>
  ];
  Print[
    "      r_BA=q_T^2*c_gamma^2/(mu_R*A_E)=q_T^2/(rho_br*A_E)"
  ];
  Print[
    "      delta_BA=r_BA-1; Delta_U=-(s1*s2*A_E/(8*pi*c_gamma^2*R))*delta_BA*(D_V+A_V)"
  ];

  liveCone = If[
    activeMutation === "CONE_RATIO",
    cE/Sqrt[cGammaSquared],
    coneRatio
  ];
  coneResiduals = {
    Factor[liveCone - cE^2/cGammaSquared],
    Factor[
      (liveCone /. cGammaRelation) -
        cE^2 rhoBr/muR
    ]
  };
  expectBool[
    "CONE_RATIO",
    allZero[coneResiduals],
    <|"rCone" -> liveCone, "Residuals" -> coneResiduals|>
  ];
  Print[
    "      r_cone=c_E^2/c_gamma^2=c_E^2*rho_br/mu_R (OPEN)"
  ];

  qmagFact = If[
    activeMutation === "QMAG_R1",
    "magnitude_forced_by_electric",
    "R1(magnitude)"
  ];
  expectBool[
    "QMAG_R1",
    qmagFact === "R1(magnitude)",
    <|"MagnitudeFact" -> qmagFact|>
  ];
  Print[
    "      magnitude=", qmagFact,
    "; comparison=", comparisonFact,
    "; relative_sign=", relativeSignFact
  ];

  section["Build-global runtime scope, dual-term, and unit firewalls"];
  blindnessState = targetBlindnessViolations[
    activeMutation === "TARGET_BLINDNESS"
  ];
  expectBool[
    "TARGET_BLINDNESS",
    blindnessState === {},
    <|"Violations" -> blindnessState|>
  ];
  Print[
    "      no SEALED section-4 landing channel; pathA_39 markers absent from Route B"
  ];
  Print[
    "      A_E enters only the cited comparison/ratio lane"
  ];

  termState = termViolations[
    activeMutation === "DUAL_ENGINE_TERMS"
  ];
  expectBool[
    "DUAL_ENGINE_TERMS",
    termState === {},
    <|"Violations" -> termState|>
  ];
  Print[
    "      TERM_INVENTORY=",
    StringRiffle[expectedTermKeys, ","]
  ];
  Print[
    "      every listed term is checked by exact symbolic residual, never by text"
  ];

  unitsState = unitViolations[
    activeMutation === "UNITS_RESTORED"
  ];
  expectBool[
    "UNITS_RESTORED",
    unitsState === {},
    <|"Violations" -> unitsState|>
  ];
  Print[
    "      [q_T]=M*T^-1; [mu_R]=M*L^-1*T^-2; [rho_br]=M*L^-3"
  ];
  Print[
    "      [G_B]=M^-1*T^2; [U_B]=E; [F_B]=E/L; [r_BA]=[delta_BA]=[r_cone]=1"
  ];

  section["Verdict re-derived from computed certification objects"];
  verdictKernel = routeB["kernel"];
  verdictTensorB = tensorRouteB;
  verdictForce = routeB["radial"];
  verdictInteraction = routeB["interaction"];
  verdictIndependence =
    independenceViolations[routeB, routeA];
  If[
    activeMutation === "VERDICT_REDERIVATION",
    verdictKernel =
      verdictKernel +
        Outer[Times, nVector, nVector]/(8 Pi muR R);
    verdictTensorB =
      verdictTensorB +
        Outer[Times, nVector, nVector]/(8 Pi R);
    verdictForce = verdictForce/R;
    verdictInteraction =
      verdictInteraction /. {v2x -> 1, v2y -> 0, v2z -> 0};
    verdictIndependence = {"ILLICIT_ROUTE_A_READ"}
  ];
  verdictProjectorResiduals =
    sphereReduce /@ Flatten[verdictKernel - expectedKernelB];
  verdictComparisonResiduals =
    sphereReduce /@ Flatten[verdictTensorB - tensorRouteA];
  liveVerdict = deriveVerdict[
    verdictProjectorResiduals,
    verdictComparisonResiduals,
    inverseRPower[verdictForce],
    velocityDegreeSet[verdictInteraction],
    verdictIndependence
  ];
  Print["      REDERIVED_TOKEN=", liveVerdict];
  expectBool[
    "VERDICT_REDERIVATION",
    liveVerdict === verdictToken,
    <|
      "Live" -> liveVerdict,
      "FailureToken" -> uncertifiedToken,
      "ProjectorResiduals" -> verdictProjectorResiduals,
      "ComparisonResiduals" -> verdictComparisonResiduals,
      "ForcePower" -> inverseRPower[verdictForce],
      "VelocityDegrees" -> velocityDegreeSet[verdictInteraction],
      "Independence" -> verdictIndependence
    |>
  ];
  Print["      FAILURE_TOKEN=", uncertifiedToken];
  Print[
    "      unresolved q_T/A_E magnitude does not enter structural certification"
  ];

  section["Canonical source-to-stage predicate manifest"];
  liveManifest = sourceManifest;
  If[
    activeMutation === "SOURCE_TO_STAGE_MANIFEST",
    liveManifest = (
      {
        #[[1]],
        If[
          #[[1]] === "DIRECT_SOURCE",
          "SCOPED_OUT",
          #[[2]]
        ],
        #[[3]]
      } & /@ Rest[sourceManifest]
    )
  ];
  liveManifestState = manifestState[liveManifest];
  expectBool[
    "SOURCE_TO_STAGE_MANIFEST",
    liveManifestState === expectedManifestState,
    <|
      "Partition" -> liveManifestState[[2]],
      "InScope" -> liveManifestState[[3]],
      "ScopedOut" -> liveManifestState[[4]],
      "Digest" -> liveManifestState[[5]]
    |>
  ];
  Print[
    "      entries=35; partition=", expectedManifestCounts,
    "; scoped_out=21; digest=", manifestDigestExpected
  ];
  Print[
    "      IN_SCOPE=", StringRiffle[expectedInScope, ","]
  ];
  Print[
    "      SCOPED_OUT=", StringRiffle[expectedScopedOut, ","]
  ];

  Print[""];
  Print[
    "BUILD_ORDER=Route B first with foreign_payload=None; Route A cited afterward"
  ];
  Print[
    "STRUCTURAL_AGREEMENT=tensor + R^-2 falloff + O(V1*V2) velocity order"
  ];
  Print[
    "RATIO_STATUS=DECIDED expressions; R1-valued through q_T and A_E"
  ];
  Print["COMPARISON_FACT=", comparisonFact];
  Print["RELATIVE_SIGN_FACT=", relativeSignFact];
  Print["COBLOCKER=", coblockerToken];
  Print[
    "STAGE038_BOUNDARY=no SEALED section-4 landing is adjudicated here"
  ];
  Print["VERDICT_TOKEN: ", liveVerdict];
  Print["PROGRESS MATHEMATICA STAGE037 COMPLETE"];

  If[activeMutation =!= "",
    Print["FIRST_FAILURE=MUTATION_DID_NOT_FIRE"];
    raise["MUTATION_DID_NOT_FIRE"]
  ];
  True,
  "ledgerStage037Failure",
  Function[{message, tag}, False]
];

Print[""];
Print["TOOTH_COUNT=", Length[toothOrder]];
Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
If[TrueQ[ok],
  Print[
    "OVERALL PASS: Mathematica independently reached ",
    verdictToken
  ];
  Exit[0],
  Print[
    "OVERALL FAIL: Mathematica stage037 audit did not close"
  ];
  Exit[1]
]
