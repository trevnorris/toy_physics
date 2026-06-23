(* Machine check for pathA_24 T0 frozen polar-OP bulk modes.

   Mathematica is used only as the independent engine for the derivation of
   the quadratic OP spectrum.  The frozen action block itself is not edited.
*)

ClearAll[
  scriptPath, stage1Root, scratchDir, freezeReport, expectedLPolLines,
  freezeText, freezeBlock, linePresentQ, eps, m, rho0, K, a, omega,
  kx, ky, kz, kw, sigma, pi1, pi2, pi3, sigmaT, pi1T, pi2T, pi3T,
  q, qT, grad, cs20, iP, kP, bDepth, p, dtvP, gradP, lPol, l2,
  hessian, kineticMatrix, stiffnessMatrices, massMatrix, expectedKinetic,
  expectedStiffness, expectedMass, k2, inertia, stiffness, massDiag,
  waveOperators, omega2, transverseIndices, longitudinalIndices, cOpSq,
  gapSq, exprString, checks, checkNames, failures, report, outPath
];

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "", $InputFileName, "tools/pathA_24_T0_op_modes.wl"];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
freezeReport = FileNameJoin[{stage1Root, "reports", "pathA_24_T0_freeze.md"}];

expectedLPolLines = {
  "  L_pol =",
  "      (1/2) m rho a^2 (D_t^v P^i)(D_t^v P^i)",
  "    - (1/2) m rho c_s^2(rho) a^2 (partial_j P^i)(partial_j P^i)",
  "    - (1/4) m rho c_s^2(rho) (P^i P^i - 1)^2."
};

freezeText = Import[freezeReport, "Text"];
freezeBlock = StringCases[
  freezeText,
  "```freeze-action\n" ~~ block__ ~~ "\n```" :> block
];
If[Length[freezeBlock] =!= 1,
  Print["FAIL: freeze-action block not found uniquely in ", freezeReport];
  Exit[1]
];
linePresentQ[line_] := StringContainsQ[First[freezeBlock], line];
If[! And @@ (linePresentQ /@ expectedLPolLines),
  Print["FAIL: frozen L_pol line missing or changed"];
  Exit[1]
];

$Assumptions = m > 0 && rho0 > 0 && K > 0 && a > 0;

q = {sigma, pi1, pi2, pi3};
qT = {sigmaT, pi1T, pi2T, pi3T};
grad = Array[g, {4, 4}];

cs20 = 5 K rho0^4/m;
iP = m rho0 a^2;
kP = m rho0 cs20 a^2;
bDepth = m rho0 cs20;

(* Uniform background: rho=rho0, v=0, and an O(4)-rotated representative
   P0=(1,0,0,0).  This transcribes the frozen L_pol terms after those
   background substitutions. *)
p = {1 + eps sigma, eps pi1, eps pi2, eps pi3};
dtvP = eps qT;
gradP = eps grad;
lPol = (
  1/2 m rho0 a^2 Total[dtvP^2]
  - 1/2 m rho0 cs20 a^2 Total[Flatten[gradP]^2]
  - 1/4 m rho0 cs20 (p.p - 1)^2
);
l2 = Expand[Coefficient[Normal[Series[lPol, {eps, 0, 2}]], eps, 2]];

hessian[expr_, vars_] := Outer[D[expr, #1, #2] &, vars, vars];
kineticMatrix = FullSimplify[hessian[l2, qT]];
stiffnessMatrices = Table[
  FullSimplify[-hessian[l2, grad[[All, direction]]]],
  {direction, 1, 4}
];
massMatrix = FullSimplify[-hessian[l2, q]];

expectedKinetic = iP IdentityMatrix[4];
expectedStiffness = kP IdentityMatrix[4];
expectedMass = DiagonalMatrix[{2 bDepth, 0, 0, 0}];

k2 = kx^2 + ky^2 + kz^2 + kw^2;
inertia = Diagonal[kineticMatrix];
stiffness = Diagonal[First[stiffnessMatrices]];
massDiag = Diagonal[massMatrix];
waveOperators = FullSimplify[
  Table[-inertia[[i]] omega^2 + stiffness[[i]] k2 + massDiag[[i]], {i, 1, 4}]
];
omega2 = FullSimplify[
  Table[(stiffness[[i]] k2 + massDiag[[i]])/inertia[[i]], {i, 1, 4}]
];

transverseIndices = Flatten[Position[FullSimplify[massDiag], 0]];
longitudinalIndices = Complement[Range[4], transverseIndices];
cOpSq = FullSimplify[stiffness[[2]]/inertia[[2]]];
gapSq = FullSimplify[massDiag[[1]]/inertia[[1]]];

checks = {
  kineticMatrix == expectedKinetic,
  And @@ (FullSimplify[# == expectedStiffness] & /@ stiffnessMatrices),
  massMatrix == expectedMass,
  transverseIndices == {2, 3, 4},
  longitudinalIndices == {1},
  FullSimplify[cOpSq == cs20],
  FullSimplify[Sqrt[cOpSq] == Sqrt[cs20]],
  FullSimplify[gapSq == 2 cs20/a^2],
  FullSimplify[omega2[[2]] == cs20 k2],
  FullSimplify[omega2[[3]] == cs20 k2],
  FullSimplify[omega2[[4]] == cs20 k2]
};
checkNames = {
  "kinetic matrix",
  "stiffness matrices",
  "mass matrix",
  "transverse count/indices",
  "longitudinal count/index",
  "c_OP^2 equals c_s0^2",
  "c_OP equals c_s0",
  "longitudinal gap",
  "pi1 dispersion",
  "pi2 dispersion",
  "pi3 dispersion"
};
failures = Pick[checkNames, Not /@ (TrueQ /@ checks), True];
If[failures =!= {},
  Print["FAIL: ", failures];
  Print["checks: ", InputForm[checks]];
  Print["truths: ", InputForm[TrueQ /@ checks]];
  Print["lPol: ", InputForm[lPol]];
  Print["l2: ", InputForm[l2]];
  Print["massDiag: ", InputForm[massDiag]];
  Print["transverseIndices: ", InputForm[transverseIndices]];
  Print["longitudinalIndices: ", InputForm[longitudinalIndices]];
  Print["cOpSq: ", InputForm[cOpSq]];
  Print["cs20: ", InputForm[cs20]];
  Print["gapSq: ", InputForm[gapSq]];
  Print["omega2: ", InputForm[omega2]];
  Exit[1]
];

exprString[x_] := ToString[InputForm[FullSimplify[x]]];

report = <|
  "schema" -> "stage1_pathA_24_T0_op_modes_mathematica/v1",
  "engine" -> "mathematica",
  "pass" -> True,
  "freeze_fidelity_guard" -> <|
    "path" -> freezeReport,
    "checked_lines" -> expectedLPolLines
  |>,
  "modes" -> <|
    "quadratic_lagrangian" -> exprString[l2],
    "I_P" -> exprString[iP],
    "K_P" -> exprString[kP],
    "c_s0_squared" -> exprString[cs20],
    "c_OP_squared" -> exprString[cOpSq],
    "c_OP" -> ("Sqrt[" <> exprString[cOpSq] <> "]"),
    "transverse_count" -> Length[transverseIndices],
    "longitudinal_count" -> Length[longitudinalIndices],
    "longitudinal_gap_squared" -> exprString[gapSq],
    "longitudinal_gap" -> ("Sqrt[" <> exprString[gapSq] <> "]"),
    "wave_operators" -> (exprString /@ waveOperators),
    "dispersion_omega_squared" -> (exprString /@ omega2),
    "mass_matrix_diag" -> (exprString /@ massDiag)
  |>
|>;

If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir, CreateIntermediateDirectories -> True]];
outPath = FileNameJoin[{scratchDir, "pathA_24_T0_op_modes_mathematica.json"}];
Export[outPath, report, "RawJSON"];

Print["T0_DIMCHECK_MATHEMATICA: PASS"];
Print["c_OP_mathematica: ", report["modes", "c_OP"]];
Print["longitudinal_gap_mathematica: ", report["modes", "longitudinal_gap"]];
Print["transverse_count_mathematica: ", report["modes", "transverse_count"]];
Print["wrote: ", outPath];
Exit[0]
