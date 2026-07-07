(* Ledger stage001 Mathematica audit: solid angles and second moments.
   Print-only, no exports.  This route builds its own parametrized normals,
   derives the induced surface element from the Gram determinant, and then
   integrates the normalized moments with native Integrate. *)

ClearAll[
  heading, subheading, cleanZero, expectZero, jacobianByColumns,
  inducedSurfaceDensity, integrateOnChart, momentOnChart, baselineResiduals,
  mutationRows, mutationMustFail
];

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

expectZero[name_, residual_] := Module[{clean = cleanZero[residual]},
  If[TrueQ[clean === 0],
    Print["PASS  ", name],
    Print["FAIL  ", name, ": residual = ", ToString[InputForm[clean]]];
    Throw[name, "ledgerStage001Failure"]
  ]
];

jacobianByColumns[map_, coordinates_] := Table[
  D[map[[row]], coordinates[[col]]],
  {row, Length[map]},
  {col, Length[coordinates]}
];

inducedSurfaceDensity[map_, coordinates_, assumptions_] := Module[{jac},
  jac = jacobianByColumns[map, coordinates];
  FullSimplify[Sqrt[Det[Transpose[jac].jac]], assumptions]
];

integrateOnChart[expr_, bounds_] := FullSimplify[Apply[Integrate, Prepend[bounds, expr]]];

momentOnChart[normal_, density_, bounds_, i_, j_, omega_] := FullSimplify[
  integrateOnChart[normal[[i]] normal[[j]] density, bounds]/omega
];

baselineResiduals[targets_] := Module[
  {
    u, v, r, s, t, patchTwo, patchThree, boundsTwo, boundsThree,
    densityTwo, densityThree, areaTwo, areaThree, built
  },

  patchTwo = {Sin[u] Cos[v], Sin[u] Sin[v], Cos[u]};
  boundsTwo = {{u, 0, Pi}, {v, 0, 2 Pi}};
  densityTwo = inducedSurfaceDensity[
    patchTwo,
    {u, v},
    0 <= u <= Pi && 0 <= v <= 2 Pi
  ];
  areaTwo = integrateOnChart[densityTwo, boundsTwo];

  patchThree = {
    Sin[r] Sin[s] Cos[t],
    Sin[r] Sin[s] Sin[t],
    Sin[r] Cos[s],
    Cos[r]
  };
  boundsThree = {{r, 0, Pi}, {s, 0, Pi}, {t, 0, 2 Pi}};
  densityThree = inducedSurfaceDensity[
    patchThree,
    {r, s, t},
    0 <= r <= Pi && 0 <= s <= Pi && 0 <= t <= 2 Pi
  ];
  areaThree = integrateOnChart[densityThree, boundsThree];

  built = <|
    "Omega_2 from native S2 Integrate" -> areaTwo,
    "<n1^2>_S2 native normalized moment" -> momentOnChart[patchTwo, densityTwo, boundsTwo, 1, 1, areaTwo],
    "<n1 n2>_S2 native normalized moment" -> momentOnChart[patchTwo, densityTwo, boundsTwo, 1, 2, areaTwo],
    "Omega_3 from native S3 Integrate" -> areaThree,
    "<n1^2>_S3 native normalized moment" -> momentOnChart[patchThree, densityThree, boundsThree, 1, 1, areaThree],
    "<n1 n2>_S3 native normalized moment" -> momentOnChart[patchThree, densityThree, boundsThree, 1, 2, areaThree]
  |>;

  Association @ KeyValueMap[#1 -> cleanZero[built[#1] - targets[#1]] &, built]
];

mutationRows[baseTargets_] := {
  <|
    "label" -> "perturb Omega_2 target 4 Pi -> 2 Pi",
    "check" -> "Omega_2 from native S2 Integrate",
    "targets" -> Join[baseTargets, <|"Omega_2 from native S2 Integrate" -> 2 Pi|>]
  |>,
  <|
    "label" -> "perturb Omega_3 target 2 Pi^2 -> 4 Pi^2",
    "check" -> "Omega_3 from native S3 Integrate",
    "targets" -> Join[baseTargets, <|"Omega_3 from native S3 Integrate" -> 4 Pi^2|>]
  |>,
  <|
    "label" -> "perturb S2 <n1^2> target 1/3 -> 1/2",
    "check" -> "<n1^2>_S2 native normalized moment",
    "targets" -> Join[baseTargets, <|"<n1^2>_S2 native normalized moment" -> 1/2|>]
  |>,
  <|
    "label" -> "perturb S3 <n1^2> target 1/4 -> 1/2",
    "check" -> "<n1^2>_S3 native normalized moment",
    "targets" -> Join[baseTargets, <|"<n1^2>_S3 native normalized moment" -> 1/2|>]
  |>,
  <|
    "label" -> "cross-moment control: perturb S2 <n1 n2> target 0 -> 1/5",
    "check" -> "<n1 n2>_S2 native normalized moment",
    "targets" -> Join[baseTargets, <|"<n1 n2>_S2 native normalized moment" -> 1/5|>]
  |>
};

mutationMustFail[row_] := Module[{mutated, outcome},
  mutated = baselineResiduals[row["targets"]];
  outcome = Catch[
    expectZero[row["check"], mutated[row["check"]]];
    "SURVIVED",
    "ledgerStage001Failure",
    Function[{name, tag}, "FAILED"]
  ];
  If[TrueQ[outcome === "FAILED"],
    Print["PASS  mutation probe: ", row["label"], " produced the required FAIL"];
    True,
    Print["FAIL  mutation probe: ", row["label"], " survived"];
    False
  ]
];

Module[{targets, baseline, baselineOutcome, probeOutcome},
  heading["ledger_stage001_solid_angle_second_moment_primitives Mathematica audit"];

  targets = <|
    "Omega_2 from native S2 Integrate" -> 4 Pi,
    "<n1^2>_S2 native normalized moment" -> 1/3,
    "<n1 n2>_S2 native normalized moment" -> 0,
    "Omega_3 from native S3 Integrate" -> 2 Pi^2,
    "<n1^2>_S3 native normalized moment" -> 1/4,
    "<n1 n2>_S3 native normalized moment" -> 0
  |>;

  subheading["Baseline exact residuals"];
  baseline = baselineResiduals[targets];
  baselineOutcome = Catch[
    Scan[Function[key, expectZero[key, baseline[key]]], Keys[baseline]];
    True,
    "ledgerStage001Failure",
    Function[{name, tag}, False]
  ];

  subheading["Able-to-fail mutation probe"];
  probeOutcome = AllTrue[mutationRows[targets], mutationMustFail];

  If[TrueQ[baselineOutcome && probeOutcome],
    Print[""];
    Print["OVERALL PASS: Mathematica derived all stage001 geometry primitives exactly"];
    Exit[0],
    Print[""];
    Print["OVERALL FAIL: Mathematica stage001 audit did not close"];
    Exit[1]
  ]
]
