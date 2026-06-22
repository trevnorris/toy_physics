(* PathA 22b Gate 0 cross-check.
   Scope: dimensional/algebraic verification only; no branch solve and no
   target-value comparison. *)

ClearAll[
  dimString, checkDim, homogeneous, L, T, M, dim0, velocity, g3,
  gamma5, btRhs, mhat0, factorizedRhs, naturalMap, naturalCorrection,
  IZ, Kstress, Ksource, factoredRatio, w, v, negStress, negSource,
  negCondition, propStress, propSource, propCondition, propRouteBCancels,
  propWeightedRatio, z0, z1, checks, report,
  outDir, jsonPath, allPass
];

dim0 = {0, 0, 0};
L = {1, 0, 0};
T = {0, 1, 0};
M = {0, 0, 1};
velocity = L - T;
g3 = 3 L - 2 T - M;
gamma5 = 5 T;
btRhs = g3 - 5 velocity;
mhat0 = -L - T - M/2;
factorizedRhs = g3 + 5 velocity - 5 L - 5 velocity;
naturalMap = dim0;
naturalCorrection = 2 L - 2 L;

dimString[d_] := Module[
  {labels = {"L", "T", "M"}, pairs, nonzero},
  pairs = Transpose[{labels, d}];
  nonzero = Select[pairs, #[[2]] =!= 0 &];
  If[
    nonzero === {},
    "1",
    StringRiffle[
      Map[
        If[#[[2]] === 1, #[[1]], #[[1]] <> "^" <> ToString[InputForm[#[[2]]]]] &,
        nonzero
      ],
      " "
    ]
  ]
];

checkDim[name_, actual_, expected_, expectedNegative_: False] := <|
  "name" -> name,
  "pass" -> If[TrueQ[expectedNegative], TrueQ[actual =!= expected], TrueQ[actual === expected]],
  "raw_equal" -> TrueQ[actual === expected],
  "expected_negative" -> expectedNegative,
  "expected" -> dimString[expected],
  "actual" -> dimString[actual]
|>;

homogeneous[name_, terms_Association] := Module[
  {values = Values[terms], expected, rawPass},
  expected = First[values];
  rawPass = AllTrue[values, TrueQ[# === expected] &];
  <|
    "name" -> name,
    "pass" -> rawPass,
    "expected" -> dimString[expected],
    "terms" -> Association @ KeyValueMap[#1 -> dimString[#2] &, terms]
  |>
];

factoredRatio = FullSimplify[(IZ*Kstress)/(IZ*Ksource)];
negStress = 1 + w;
negSource = 1 + w^2;
negCondition = Expand[negStress*(1 + v^2) - (1 + v)*(1 + w^2)];
propStress = 2*(1 + w^2);
propSource = 1 + w^2;
propCondition = Expand[propStress*(1 + v^2) - 2*(1 + v^2)*propSource];
propRouteBCancels = TrueQ[propCondition === 0];
propWeightedRatio = FullSimplify[
  (z0*propStress /. w -> 0) + (z1*propStress /. w -> 1)
];
propWeightedRatio = FullSimplify[
  propWeightedRatio / ((z0*propSource /. w -> 0) + (z1*propSource /. w -> 1))
];

checks = {
  checkDim["Gamma5", gamma5, 5 T],
  checkDim["G/c^5", btRhs, -2 L + 3 T - M],
  checkDim["mhat0 from odd law", mhat0, (btRhs - gamma5)/2],
  checkDim["factorized normalization right side", factorizedRhs, 2 mhat0],
  homogeneous[
    "dimensionful mhat0 law",
    <|"mhat0^2*Gamma5" -> 2 mhat0 + gamma5, "G/c^5" -> btRhs|>
  ],
  checkDim[
    "direct dimensionless mhat reading fails odd law",
    2 naturalMap + gamma5,
    btRhs,
    True
  ],
  checkDim["natural correction a^2/r^2", naturalCorrection, dim0],
  <|
    "name" -> "common scalar I_Z cancellation",
    "pass" -> TrueQ[factoredRatio === Kstress/Ksource],
    "actual" -> ToString[InputForm[factoredRatio]],
    "expected" -> ToString[InputForm[Kstress/Ksource]]
  |>,
  <|
    "name" -> "non-factorizing weighted negative control",
    "pass" -> TrueQ[negCondition =!= 0],
    "condition_residual" -> ToString[InputForm[Factor[negCondition]]],
    "expected" -> "nonzero residual"
  |>,
  <|
    "name" -> "weighted proportional control",
    "pass" -> TrueQ[propRouteBCancels && propWeightedRatio === 2],
    "outcome" -> If[TrueQ[propRouteBCancels && propWeightedRatio === 2], "CANCELS", "DOES_NOT_CANCEL"],
    "condition_residual" -> ToString[InputForm[propCondition]],
    "ratio_for_arbitrary_two_point_Z" -> ToString[InputForm[propWeightedRatio]],
    "expected" -> "CANCELS by route (b): Z-independent ratio"
  |>
};

allPass = AllTrue[checks, TrueQ[#["pass"]] &];
report = <|
  "schema" -> "stage1_pathA_22b_gate0_mathematica_crosscheck/v1",
  "scope" -> "dimensional and algebraic cross-check only",
  "pass" -> allPass,
  "checks" -> checks,
  "outcomes" -> <|
    "0a" -> "MHAT_DIMENSIONFUL_CONFIRMED",
    "0b_source_assessment" -> "DOES_NOT_CANCEL (NOT_ESTABLISHED — sources do not establish either cancellation route; a later Gate-4 action-level derivation could still find one)"
  |>
|>;

outDir = FileNameJoin[{"software", "stage1_solver", "_scratch"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
jsonPath = FileNameJoin[{outDir, "pathA_22b_gate0_mathematica_crosscheck.json"}];
Export[jsonPath, report, "RawJSON"];
Print["wrote ", jsonPath];
Print["pathA_22b Gate 0 Mathematica cross-check: ", Count[checks[[All, "pass"]], True], "/", Length[checks], " checks"];
If[TrueQ[report["pass"]], Exit[0], Exit[1]]
