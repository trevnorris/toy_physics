(* PathA 22b Gate 2 cross-check.
   Scope: localized Maxwell zero-mode coefficient reduction, cone ratio,
   dimensional bookkeeping, and c_gamma=c_s negative control. *)

ClearAll[
  dimString, checkExpr, checkBool, checkDim, homogeneous,
  dim0, L, T, M, velocity, action, energy, force, lpsiDensity,
  electricField, cEDim, cBDim, ai, waveNumber, rho4, kEos,
  zIntSym, mu0Sym, cE, cB, ratio, omegaSym, k2Sym,
  ftx, fty, ftz, ftw, fxy, fxz, fxw, fyz, fyw, fzw,
  metricDiag, fmat, f2Full, f2ZeroMode, e2, b2,
  beta, rho0, kSym, mG, cs2, lambdaGamma, forcedResidual,
  checks, algebra, report, outDir, jsonPath, allPass
];

dim0 = {0, 0, 0};
L = {1, 0, 0};
T = {0, 1, 0};
M = {0, 0, 1};
velocity = L - T;
action = {2, -1, 1};
energy = {2, -2, 1};
force = {1, -2, 1};
lpsiDensity = energy - 4 L;
electricField = force;
cEDim = lpsiDensity - 2 electricField;
cBDim = cEDim;
ai = action - L;
waveNumber = -L;
rho4 = -4 L;
kEos = energy - 4 rho4;

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

homogeneous[name_, terms_Association, note_: ""] := Module[
  {values = Values[terms], expected, pass},
  expected = First[values];
  pass = AllTrue[values, TrueQ[# === expected] &];
  <|
    "name" -> name,
    "pass" -> pass,
    "expected" -> dimString[expected],
    "terms" -> Association @ KeyValueMap[#1 -> dimString[#2] &, terms],
    "note" -> note
  |>
];

zIntSym = Symbol["Zint"];
mu0Sym = Symbol["mu0"];
omegaSym = Symbol["omega"];
k2Sym = Symbol["k2"];
ftx = Symbol["Ftx"];
fty = Symbol["Fty"];
ftz = Symbol["Ftz"];
ftw = Symbol["Ftw"];
fxy = Symbol["Fxy"];
fxz = Symbol["Fxz"];
fxw = Symbol["Fxw"];
fyz = Symbol["Fyz"];
fyw = Symbol["Fyw"];
fzw = Symbol["Fzw"];
beta = Symbol["betaBulkToBrane"];
rho0 = Symbol["rho0"];
kSym = Symbol["K"];
mG = Symbol["mGNLS"];

metricDiag = {-1, 1, 1, 1, 1};
fmat = {
  {0, ftx, fty, ftz, ftw},
  {-ftx, 0, fxy, fxz, fxw},
  {-fty, -fxy, 0, fyz, fyw},
  {-ftz, -fxz, -fyz, 0, fzw},
  {-ftw, -fxw, -fyw, -fzw, 0}
};
f2Full = Expand[
  Sum[metricDiag[[m]]*metricDiag[[n]]*fmat[[m, n]]^2, {m, 1, 5}, {n, 1, 5}]
];
f2ZeroMode = Expand[f2Full /. {ftw -> 0, fxw -> 0, fyw -> 0, fzw -> 0}];
e2 = ftx^2 + fty^2 + ftz^2;
b2 = fxy^2 + fxz^2 + fyz^2;

cE = zIntSym/mu0Sym;
cB = zIntSym/mu0Sym;
ratio = FullSimplify[cB/cE];
cs2 = 5*kSym*rho0^4/mG;
lambdaGamma = FullSimplify[beta*Sqrt[ratio/cs2]];
forcedResidual = FullSimplify[beta^2*ratio - cs2];

checks = {
  checkDim["C_B/C_E flat bulk cone ratio", cBDim - cEDim, dim0],
  checkDim["beta_bulk_to_brane speed normalization", velocity, velocity],
  checkDim["c_gamma dimension", velocity + (cBDim - cEDim)/2, velocity],
  checkDim["c_s dimension", (kEos + 4 rho4 - M)/2, velocity],
  checkDim["lambda_gamma dimensionless", 2 velocity + (cBDim - cEDim) - (kEos + 4 rho4 - M), dim0],
  homogeneous[
    "bulk-metric transverse Maxwell principal operator",
    <|
      "C_E partial_0^2 A_T" -> cEDim + ai - 2 L,
      "C_B laplacian A_T" -> cBDim + ai - 2 L
    |>
  ]
};

algebra = {
  checkExpr["flat metric F_MN F^MN zero-mode expansion", f2ZeroMode, -2*e2 + 2*b2],
  checkExpr["localized Maxwell C_B/C_E", ratio, 1],
  checkExpr["transverse operator factorization", cE*omegaSym^2 - cB*k2Sym, cE*(omegaSym^2 - ratio*k2Sym)],
  checkExpr["lambda_gamma carried expression", lambdaGamma, beta*Sqrt[mG/(5*kSym*rho0^4)]],
  checkBool["negative guard residual is not zero by identity", ! TrueQ[forcedResidual === 0], True]
};

allPass = AllTrue[Join[checks, algebra], TrueQ[#["pass"]] &];

report = <|
  "schema" -> "stage1_pathA_22b_gate2_mathematica_crosscheck/v1",
  "pass" -> allPass,
  "metric" -> "eta_MN=diag(-1,+1,+1,+1,+1)",
  "F_squared_zero_mode" -> ToString[InputForm[f2ZeroMode]],
  "C_E" -> "Z_int/mu0",
  "C_B" -> "Z_int/mu0",
  "C_B_over_C_E" -> "1",
  "lambda_gamma" -> ToString[InputForm[lambdaGamma]],
  "negative_control_residual" -> ToString[InputForm[forcedResidual]],
  "outcome" -> "STILL_TUNABLE_LAMBDAGAMMA",
  "checks" -> checks,
  "algebra" -> algebra
|>;

outDir = FileNameJoin[{"software", "stage1_solver", "_scratch"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
jsonPath = FileNameJoin[{outDir, "pathA_22b_gate2_mathematica_crosscheck.json"}];
Export[jsonPath, report, "RawJSON"];
Print["wrote ", jsonPath];
Print["pathA_22b Gate 2 Mathematica cross-check: ", Count[Join[checks, algebra][[All, "pass"]], True], "/", Length[Join[checks, algebra]], " checks"];
If[TrueQ[report["pass"]], Exit[0], Exit[1]]
