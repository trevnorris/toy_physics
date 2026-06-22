(* Narrow algebraic/dimensional cross-check for pathA_22a.
   Scope: verifies the restored-unit homogeneity and symbolic xi reduction under
   the same assumptions as the Python harness. It does not prove that chi_Q,
   g_G, g_mhat, lambda_gamma, or branch forms are physically fixed. *)

ClearAll[
  dimString, factorToReach, checkDim, homogeneous,
  dim0, L, T, M, velocity, g3, k0, mass0, omegaU, omegaW, rMix,
  gPort, delta, sMix, qMix, pMix, z0, z2, z4, n0, n2, n4,
  p0Faithful, freqNorm, mhat, target, p0, sPort,
  lambdaGamma, gMhat, gG, chiQ, p0sym, cs, a, m, gMonomial,
  mhatMonomial, denominator, xiTimes, checks, algebra, report,
  outDir, jsonPath, allDimPass, allAlgPass
];

dim0 = {0, 0, 0};
L = {1, 0, 0};
T = {0, 1, 0};
M = {0, 0, 1};
velocity = L - T;
g3 = 3 L - 2 T - M;
k0 = M - L - 2 T;
mass0 = M - L;
omegaU = -T;
omegaW = -T;
rMix = -2 T;
gPort = (k0 - 2 T)/2;
delta = 2 omegaU + 2 omegaW;
sMix = 2 omegaU;
qMix = 2 gPort + 2 omegaW;
pMix = 2 omegaU + gPort;
z0 = qMix - delta;
z2 = qMix + sMix - 2 delta;
z4 = qMix + 2 sMix - 3 delta;
n0 = 2 pMix - 2 delta;
n2 = pMix + pMix + sMix - 3 delta;
n4 = 2 delta + 2 gPort - 4 delta;
p0Faithful = n0 - k0;
freqNorm = 2 velocity - 2 L;
mhat = -L - T - M/2;
target = g3 + 5 velocity - 5 L - 5 velocity;
p0 = p0Faithful + freqNorm;
sPort = dim0;

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

factorToReach[expected_, actual_] := expected - actual;

checkDim[name_, actual_, expected_, note_: ""] := <|
  "name" -> name,
  "pass" -> TrueQ[actual === expected],
  "expected" -> dimString[expected],
  "actual" -> dimString[actual],
  "factor_needed_to_reach_expected" -> dimString[factorToReach[expected, actual]],
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

checks = {
  homogeneous[
    "D0=K-B0-Z0",
    <|"K" -> k0, "B0" -> k0, "Z0" -> k0|>
  ],
  homogeneous[
    "Delta=Omega_U^2*Omega_W^2-R^2",
    <|"Omega_U^2*Omega_W^2" -> delta, "R^2" -> 2 rMix|>
  ],
  homogeneous[
    "P=Omega_U^2*g_W+R*g_U",
    <|"Omega_U^2*g_W" -> pMix, "R*g_U" -> rMix + gPort|>
  ],
  checkDim["Z0=Q/Delta", z0, k0],
  checkDim["Z2 formula", z2, mass0],
  checkDim["Z4 formula", z4, mass0 + 2 T],
  checkDim["N0=P^2/Delta^2", n0, mass0],
  checkDim["N2 formula", n2, mass0 + 2 T],
  checkDim["N4 formula", n4, mass0 + 4 T],
  checkDim["faithful P0=N0/D0 before normalization", p0Faithful, 2 T],
  checkDim["frequency normalization (c_s/a)^2", freqNorm, -2 T],
  checkDim["normalized P0=(c_s/a)^2*N0/D0", p0, dim0],
  checkDim["mhat0", mhat, -L - T - M/2],
  checkDim["mhat0^2", 2 mhat, target],
  checkDim["S_port/chi_Q", sPort, dim0],
  homogeneous[
    "R_norm",
    <|"mhat0^2*S_port*P0" -> 2 mhat + sPort + p0, "GR_target" -> target|>
  ],
  checkDim[
    "planted missing a^5 denominator control",
    target + 5 L,
    target,
    "This must fail; the pass flag is expected to be False."
  ],
  checkDim[
    "wrong N0 reduced stiffness assertion control",
    n0,
    k0,
    "This must fail; N0 is derived from P^2/Delta^2 as reduced mass."
  ]
};

gMonomial = a*cs^2/m;
mhatMonomial = cs/(a^2*Sqrt[m]);
denominator = gMonomial*gG*cs^5/(a^5*(lambdaGamma*cs)^5);
xiTimes = FullSimplify[(mhatMonomial^2*gMhat^2*chiQ*p0sym)/denominator];
algebra = {
  <|
    "name" -> "xi*S_port*P0 keeps dimensionless residual multipliers",
    "pass" -> TrueQ[xiTimes === p0sym*chiQ*gMhat^2*lambdaGamma^5/gG],
    "actual" -> ToString[InputForm[xiTimes]],
    "expected" -> "p0sym*chiQ*gMhat^2*lambdaGamma^5/gG"
  |>
};

allDimPass =
  AllTrue[Drop[checks, -2], TrueQ[#["pass"]] &] &&
  TrueQ[checks[[-2]]["pass"] === False] &&
  TrueQ[checks[[-1]]["pass"] === False];
allAlgPass = AllTrue[algebra, TrueQ[#["pass"]] &];
report = <|
  "schema" -> "stage1_pathA_22a_dimensional_skeleton_mathematica_crosscheck/v1",
  "scope" -> "dimensional/algebraic only under Python harness assumptions",
  "pass" -> TrueQ[allDimPass && allAlgPass],
  "checks" -> checks,
  "algebraic_checks" -> algebra,
  "summary" -> <|
    "dimensional_total" -> Length[checks],
    "dimensional_passed_including_expected_negative" -> Count[Drop[checks, -2][[All, "pass"]], True] + If[TrueQ[checks[[-2]]["pass"] === False], 1, 0] + If[TrueQ[checks[[-1]]["pass"] === False], 1, 0],
    "algebraic_total" -> Length[algebra],
    "algebraic_passed" -> Count[algebra[[All, "pass"]], True]
  |>
|>;

outDir = FileNameJoin[{"software", "stage1_solver", "_scratch"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
jsonPath = FileNameJoin[{outDir, "pathA_22a_dimensional_skeleton_mathematica_crosscheck.json"}];
Export[jsonPath, report, "RawJSON"];
Print["wrote ", jsonPath];
Print[
  "pathA_22a Mathematica cross-check: ",
  report["summary", "dimensional_passed_including_expected_negative"], "/",
  report["summary", "dimensional_total"], " dimensional including expected negative; ",
  report["summary", "algebraic_passed"], "/", report["summary", "algebraic_total"], " algebraic"
];
If[TrueQ[report["pass"]], Exit[0], Exit[1]]
