(* PathA 22b Gate 4 cross-check.
   Scope: stress-lane algebra, source-map blocker propagation, alpha_J
   non-cancellation, and comparator controls. *)

ClearAll[
  dimString, checkDim, checkExpr, checkBool, condition,
  dim0, L, T, M, velocity, actionDim, rho4, dwDim, n3Dim,
  g3Dim, mhat0Dim, stressDim, stressIntDim, gCondDim, gGDim,
  w, v, eps, chiN, rhoInf4, kStress, mutatedStress,
  negativeResidual, mutatedResidual, theta1, theta2, alpha1, alpha2,
  if12, j1, j2, n3, mG, cGamma, hbar, a, cS, q1, q2, cF,
  mass1, mass2, gCond, expectedGCond, residualAlpha, checks,
  algebra, report, outDir, jsonPath, allPass
];

dim0 = {0, 0, 0};
L = {1, 0, 0};
T = {0, 1, 0};
M = {0, 0, 1};
velocity = L - T;
actionDim = 2 L - T + M;
rho4 = -4 L;
dwDim = L;
n3Dim = -3 L;
g3Dim = 3 L - 2 T - M;
mhat0Dim = -L - T - M/2;
stressDim = rho4;
stressIntDim = stressDim + dwDim;
gCondDim = 4 velocity + M - n3Dim - 2 actionDim;
gGDim = gCondDim + M - L - 2 velocity;

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

checkDim[name_, actual_, expected_, note_: ""] := <|
  "name" -> name,
  "pass" -> TrueQ[actual === expected],
  "expected" -> dimString[expected],
  "actual" -> dimString[actual],
  "note" -> note
|>;

checkExpr[name_, actual_, expected_, note_: ""] := Module[
  {residual = FullSimplify[actual - expected]},
  <|
    "name" -> name,
    "pass" -> TrueQ[residual === 0],
    "expected" -> ToString[InputForm[expected]],
    "actual" -> ToString[InputForm[FullSimplify[actual]]],
    "residual" -> ToString[InputForm[residual]],
    "note" -> note
  |>
];

checkBool[name_, actual_, expected_, note_: ""] := <|
  "name" -> name,
  "pass" -> TrueQ[actual === expected],
  "expected" -> ToString[InputForm[expected]],
  "actual" -> ToString[InputForm[actual]],
  "note" -> note
|>;

condition[ks_, kq_] := FullSimplify[ks*(kq /. w -> v) - (ks /. w -> v)*kq];

kStress = chiN[w]*rhoInf4[w];
mutatedStress = kStress + eps*w;
negativeResidual = condition[1 + w, 1 + w^2];
mutatedResidual = condition[mutatedStress, kStress];

q1 = theta1*j1/n3;
q2 = theta2*j2/n3;
cF = FullSimplify[mG*n3*q1*q2*if12/(4*Pi)];
mass1 = alpha1*hbar*j1/cGamma^2;
mass2 = alpha2*hbar*j2/cGamma^2;
gCond = FullSimplify[cF/(mass1*mass2)];
expectedGCond = FullSimplify[cGamma^4*mG*theta1*theta2*if12/(4*Pi*n3*alpha1*alpha2*hbar^2)];
residualAlpha = FullSimplify[1/(1/(alpha1*alpha2))];

checks = {
  checkDim["K_stress(w)", stressDim, rho4],
  checkDim["int K_stress dw", stressIntDim, n3Dim],
  checkDim["mhat0 dimensional carrier", mhat0Dim, -L - T - M/2],
  checkDim["conditional G dimension", gCondDim, g3Dim],
  checkDim["g_G dimension", gGDim, dim0]
};

algebra = {
  checkExpr[
    "C_F substitution",
    cF,
    mG*theta1*theta2*j1*j2*if12/(4*Pi*n3),
    "Retains both defect flux factors."
  ],
  checkExpr[
    "G_cond pair expression",
    gCond,
    expectedGCond,
    "No square-root reduction."
  ],
  checkBool[
    "negative control does not cancel",
    ! TrueQ[FullSimplify[negativeResidual] === 0],
    True
  ],
  checkBool[
    "mutated real-stress control does not cancel",
    ! TrueQ[FullSimplify[mutatedResidual] === 0],
    True
  ],
  checkBool[
    "alpha residual remains when source lane is alpha-free",
    ! TrueQ[residualAlpha === 1],
    True
  ]
};

allPass = AllTrue[Join[checks, algebra], TrueQ[#["pass"]] &];

report = <|
  "schema" -> "stage1_pathA_22b_gate4_mathematica_crosscheck/v1",
  "pass" -> allPass,
  "K_stress" -> ToString[InputForm[kStress]],
  "K_source_status" -> "BLOCKED_NEEDS_SOURCE_MAP_PROVENANCE",
  "G_cond" -> ToString[InputForm[gCond]],
  "negative_control_residual" -> ToString[InputForm[negativeResidual]],
  "mutated_kernel_residual" -> ToString[InputForm[mutatedResidual]],
  "alpha_residual_if_source_alpha_free" -> ToString[InputForm[residualAlpha]],
  "overall_verdict" -> "BLOCKED_NEEDS_SOURCE_MAP_PROVENANCE",
  "checks" -> checks,
  "algebra" -> algebra
|>;

outDir = FileNameJoin[{"software", "stage1_solver", "_scratch"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
jsonPath = FileNameJoin[{outDir, "pathA_22b_gate4_mathematica_crosscheck.json"}];
Export[jsonPath, report, "RawJSON"];
Print["wrote ", jsonPath];
Print["pathA_22b Gate 4 Mathematica cross-check: ", Count[Join[checks, algebra][[All, "pass"]], True], "/", Length[Join[checks, algebra]], " checks"];
If[TrueQ[report["pass"]], Exit[0], Exit[1]]
