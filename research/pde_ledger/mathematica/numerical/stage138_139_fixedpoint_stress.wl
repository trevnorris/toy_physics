(* Multi-seed numerical stress harness for stages 138-139. *)

ClearAll["Global`*"];

rootDir = DirectoryName[$InputFileName];
configPath = ExpandFileName @ If[
  ValueQ[stage138139FixedpointConfigOverride],
  stage138139FixedpointConfigOverride,
  FileNameJoin[{
    rootDir, "..", "..", "..", "scripts", "numerical",
    "stage138_139_fixedpoint_samples.json"
  }]
];
config = Import[ExpandFileName[configPath], "RawJSON"];

If[config["schema"] =!= "moving_throat_numerical_stage138_139_v1",
  Print["Unexpected config schema."];
  Exit[1];
];

fmt[x_] := ToString @ NumberForm[N[x, 14], {Infinity, 12}, ExponentFunction -> (Null &)];
nearQ[lhs_, rhs_, tol_] := Abs[lhs - rhs] <= tol (1 + Abs[rhs]);

require[label_, condition_, detail_] := Module[{status},
  status = If[TrueQ[condition], "PASS", "FAIL"];
  Print["[", status, "] ", label, ": ", detail];
  If[!TrueQ[condition], Throw[$Failed]]
];

const = config["constants"];
rF1 = const["rF1"];
gStar = const["g_star"];
piStar = const["Pi_star"];
sigma0Star = const["Sigma0_star"];
tHatStar = const["T_hat_star"];
sigma0CanExpected = const["Sigma0_can_expected"];
sCanExpected = const["S_can_expected"];
piCanExpected = const["Pi_can_expected"];
tHatCanExpected = const["T_hat_can_expected"];

n = config["grid_points"];
maxIterations = config["max_iterations"];
fixedPointTol = config["fixed_point_tol"];
profileTol = config["profile_tolerance"];
momentTol = config["moment_tolerance"];
rootTol = config["root_tolerance"];
bisectionIterations = config["bisection_iterations"];

xGrid = Subdivide[0.0, 1.0, n - 1];
dx = xGrid[[2]] - xGrid[[1]];
weights = ConstantArray[dx, n];
weights[[1]] = dx/2;
weights[[-1]] = dx/2;
kappa = N[Pi/2];
cGrid = Cos[(Pi*xGrid)/2];
kqGrid = Cosh[kappa*(1 - xGrid)]/Cosh[kappa];

normalize[sig_List] := sig/Total[sig*weights];

seedProfile[seed_] := Module[{kind, lam},
  kind = seed["kind"];
  Switch[kind,
    "canonical_exponential",
      normalize[piStar*Exp[-piStar*xGrid]],
    "uniform",
      normalize[ConstantArray[1.0, n]],
    "derivative_match",
      normalize[(Pi/2)*Cos[(Pi*xGrid)/2]],
    "convex_blend",
      lam = seed["lambda"];
      normalize[(1 - lam)*ConstantArray[1.0, n] + lam*(Pi/2)*Cos[(Pi*xGrid)/2]],
    _,
      Print["Unknown seed kind: ", kind];
      Exit[1]
  ]
];

canonicalSeed = seedProfile[<|"kind" -> "canonical_exponential"|>];

tsOperator[sig_List] := Module[{cumSig, cumY},
  cumSig = Accumulate[sig*weights];
  cumY = Accumulate[xGrid*sig*weights];
  cumY + xGrid*(Last[cumSig] - cumSig)
];

tqOperator[sig_List] := Module[{aTerm, bTerm},
  aTerm = Accumulate[Sinh[kappa*xGrid]*sig*weights];
  bTerm = Reverse[Accumulate[Reverse[Cosh[kappa*(1 - xGrid)]*sig*weights]]];
  (Cosh[kappa*(1 - xGrid)]*aTerm + Sinh[kappa*xGrid]*bTerm)/(kappa*Cosh[kappa])
];

gMoment[sig_List] := Total[cGrid*sig*weights];
sMoment[sig_List] := Total[kqGrid*sig*weights];
rMoment[sig_List] := ((gMoment[sig] - rF1)^2)/(1 + rF1^2);
phi[sig_List, sigma0_?NumericQ] := sigma0*(tsOperator[sig] - rMoment[sig]*tqOperator[sig]);
nextSigma[sig_List, sigma0_?NumericQ] := Module[{ph, phShift},
  ph = phi[sig, sigma0];
  phShift = ph - Min[ph];
  normalize[Exp[-phShift]]
];

solveFixedPoint[sigma0_?NumericQ, initial_List, maxIt_: maxIterations, tol_: fixedPointTol] := Module[
  {sig, sigNew, err = Infinity, iter},
  sig = initial;
  For[iter = 1, iter <= maxIt, iter++,
    sigNew = nextSigma[sig, sigma0];
    err = Max[Abs[sigNew - sig]];
    sig = sigNew;
    If[err < tol,
      Return[<|"sigma" -> sig, "iterations" -> iter, "error" -> err|>]
    ];
  ];
  Print["[FAIL] fixed-point convergence: sigma0=", fmt[sigma0], ", error=", fmt[err]];
  Throw[$Failed]
];

summaryAt[sigma0_?NumericQ, initial_List] := Module[{fp, sig, gFp, sFp, rFp, piFp, tHat},
  fp = solveFixedPoint[sigma0, initial];
  sig = fp["sigma"];
  gFp = gMoment[sig];
  sFp = sMoment[sig];
  rFp = rMoment[sig];
  piFp = sigma0*(1 - rFp*sFp);
  tHat = Sqrt[9*sigma0/20];
  <|"sigma" -> sig, "iterations" -> fp["iterations"], "error" -> fp["error"],
    "g" -> gFp, "S" -> sFp, "R" -> rFp, "Pi" -> piFp, "T_hat" -> tHat|>
];

gCache = <||>;
gFpDefault[sigma0_?NumericQ] := Module[{key},
  key = ToString[NumberForm[sigma0, {Infinity, 14}, ExponentFunction -> (Null &)]];
  If[!KeyExistsQ[gCache, key],
    gCache[key] = summaryAt[sigma0, canonicalSeed]["g"];
  ];
  gCache[key]
];

checkSeedConsistency[case_] := Module[{sigma0, summaries, refName, refSummary, dg, pred, direct},
  sigma0 = case["sigma0"];
  Print["\n=== ", case["name"], " (", case["kind"], ") ==="];
  Print["sigma0 = ", fmt[sigma0]];
  If[KeyExistsQ[case, "assumptions"],
    Print["assumptions:"];
    Scan[Print["  - ", #] &, case["assumptions"]];
  ];

  summaries = Table[
    With[{seed = seedAssoc},
      Module[{initial, summary},
        initial = seedProfile[seed];
        summary = summaryAt[sigma0, initial];
        Print[
          "seed ", seed["name"], ": iterations=", summary["iterations"],
          ", max_residual=", fmt[summary["error"]],
          ", g=", fmt[summary["g"]],
          ", S=", fmt[summary["S"]],
          ", R=", fmt[summary["R"]],
          ", Pi=", fmt[summary["Pi"]]
        ];
        {seed["name"], summary}
      ]
    ],
    {seedAssoc, config["seeds"]}
  ];

  {refName, refSummary} = First[summaries];
  Scan[
    Function[{pair},
      Module[{name, summary, profileDiff},
        {name, summary} = pair;
        If[name === refName, Return[]];
        profileDiff = Max[Abs[summary["sigma"] - refSummary["sigma"]]];
        require["profile agreement vs " <> refName <> " for " <> name, profileDiff <= profileTol,
          "L_inf=" <> fmt[profileDiff]];
        require["g agreement vs " <> refName <> " for " <> name, nearQ[summary["g"], refSummary["g"], momentTol],
          "g=" <> fmt[summary["g"]] <> ", ref=" <> fmt[refSummary["g"]]];
        require["S agreement vs " <> refName <> " for " <> name, nearQ[summary["S"], refSummary["S"], momentTol],
          "S=" <> fmt[summary["S"]] <> ", ref=" <> fmt[refSummary["S"]]];
        require["R agreement vs " <> refName <> " for " <> name, nearQ[summary["R"], refSummary["R"], momentTol],
          "R=" <> fmt[summary["R"]] <> ", ref=" <> fmt[refSummary["R"]]];
        require["Pi agreement vs " <> refName <> " for " <> name, nearQ[summary["Pi"], refSummary["Pi"], momentTol],
          "Pi=" <> fmt[summary["Pi"]] <> ", ref=" <> fmt[refSummary["Pi"]]];
      ]
    ],
    Rest[summaries]
  ];

  If[nearQ[sigma0, sigma0Star, rootTol],
    dg = refSummary["g"] - gStar;
    pred = -dg/Sqrt[1 + rF1^2] + dg^2/(1 + rF1^2);
    direct = refSummary["R"] - 1/4;
    require["Stage-138 transport-law consistency", nearQ[pred, direct, momentTol],
      "pred=" <> fmt[pred] <> ", direct=" <> fmt[direct]];
  ];
];

bracketRoot[lo_?NumericQ, hi_?NumericQ] := Module[{flo, fhi},
  flo = gFpDefault[lo] - gStar;
  fhi = gFpDefault[hi] - gStar;
  require["root bracket sign change", flo == 0 || fhi == 0 || flo*fhi < 0,
    "[" <> fmt[lo] <> ", " <> fmt[hi] <> "] -> flo=" <> fmt[flo] <> ", fhi=" <> fmt[fhi]];
  {flo, fhi}
];

bisectRoot[lo_?NumericQ, hi_?NumericQ] := Module[{left, right, flo, fhi, mid, fm},
  {flo, fhi} = bracketRoot[lo, hi];
  left = lo; right = hi;
  Do[
    mid = (left + right)/2;
    fm = gFpDefault[mid] - gStar;
    If[flo == 0, Return[left]];
    If[fhi == 0, Return[right]];
    If[flo*fm <= 0,
      right = mid;
      fhi = fm;,
      left = mid;
      flo = fm;
    ],
    {bisectionIterations}
  ];
  (left + right)/2
];

checkStage139[] := Module[{scanPoints, scanValues, diffs, roots, refRoot, final},
  Print["\n=== stage139_scan_and_root ==="];
  scanPoints = config["stage139_scan_points"];
  scanValues = gFpDefault /@ scanPoints;
  Do[
    Print["scan sigma0=", fmt[scanPoints[[i]]], " -> g_fp=", fmt[scanValues[[i]]]],
    {i, Length[scanPoints]}
  ];
  diffs = Differences[scanValues];
  require["monotone increase on analyzed scan grid", Min[diffs] > 0,
    "min diff=" <> fmt[Min[diffs]]];

  roots = Table[
    Module[{lo = bracket[[1]], hi = bracket[[2]], root},
      root = bisectRoot[lo, hi];
      Print["root from bracket [", fmt[lo], ", ", fmt[hi], "] = ", fmt[root]];
      root
    ],
    {bracket, config["stage139_root_brackets"]}
  ];

  refRoot = First[roots];
  Do[
    require["root consistency vs bracket 0 for bracket " <> ToString[i - 1],
      nearQ[roots[[i]], refRoot, rootTol],
      "root=" <> fmt[roots[[i]]] <> ", ref=" <> fmt[refRoot]],
    {i, 2, Length[roots]}
  ];

  final = summaryAt[refRoot, canonicalSeed];
  require["g(root) = g_star", nearQ[final["g"], gStar, momentTol],
    "g=" <> fmt[final["g"]] <> ", g_star=" <> fmt[gStar]];
  require["R(root) = 1/4", nearQ[final["R"], 0.25, momentTol],
    "R=" <> fmt[final["R"]]];
  require["Sigma0_can matches expected audit value", nearQ[refRoot, sigma0CanExpected, rootTol],
    "root=" <> fmt[refRoot] <> ", expected=" <> fmt[sigma0CanExpected]];
  require["S_can matches expected audit value", nearQ[final["S"], sCanExpected, momentTol],
    "S=" <> fmt[final["S"]] <> ", expected=" <> fmt[sCanExpected]];
  require["Pi_can matches expected audit value", nearQ[final["Pi"], piCanExpected, momentTol],
    "Pi=" <> fmt[final["Pi"]] <> ", expected=" <> fmt[piCanExpected]];
  require["T_hat_can matches expected audit value", nearQ[final["T_hat"], tHatCanExpected, momentTol],
    "T_hat=" <> fmt[final["T_hat"]] <> ", expected=" <> fmt[tHatCanExpected]];

  Print["Sigma0_can = ", fmt[refRoot]];
  Print["g_can      = ", fmt[final["g"]]];
  Print["S_can      = ", fmt[final["S"]]];
  Print["R_can      = ", fmt[final["R"]]];
  Print["Pi_can     = ", fmt[final["Pi"]]];
  Print["T_hat_can  = ", fmt[final["T_hat"]]];
  Print[""];
  Print["Relative shifts from original canonical point:"];
  Print["Sigma0 ratio - 1 = ", fmt[refRoot/sigma0Star - 1]];
  Print["Pi ratio - 1     = ", fmt[final["Pi"]/piStar - 1]];
  Print["T_hat ratio - 1  = ", fmt[final["T_hat"]/tHatStar - 1]];

  Print["\n=== stage139_negative_controls ==="];
  Scan[
    Function[{control},
      Module[{sigma0, summary, side, gGap, rGap},
        sigma0 = control["sigma0"];
        summary = summaryAt[sigma0, canonicalSeed];
        side = control["side"];
        gGap = Abs[summary["g"] - gStar];
        rGap = Abs[summary["R"] - 1/4];
        Print[
          control["name"], ": sigma0=", fmt[sigma0],
          ", g=", fmt[summary["g"]],
          ", R=", fmt[summary["R"]],
          ", Pi=", fmt[summary["Pi"]]
        ];
        require[control["name"] <> " stays away from the compensated moment",
          gGap >= control["min_g_gap"],
          "|g-g_star|=" <> fmt[gGap] <> ", min gap=" <> fmt[control["min_g_gap"]]];
        require[control["name"] <> " stays away from R=1/4",
          rGap >= control["min_R_gap"],
          "|R-1/4|=" <> fmt[rGap] <> ", min gap=" <> fmt[control["min_R_gap"]]];
        If[side === "below",
          require[control["name"] <> " remains on the pre-compensation side",
            summary["g"] < gStar && summary["R"] > 1/4,
            "g=" <> fmt[summary["g"]] <> ", g_star=" <> fmt[gStar] <> ", R=" <> fmt[summary["R"]]],
          require[control["name"] <> " remains on the post-compensation side",
            summary["g"] > gStar && summary["R"] < 1/4,
            "g=" <> fmt[summary["g"]] <> ", g_star=" <> fmt[gStar] <> ", R=" <> fmt[summary["R"]]]
        ]
      ]
    ],
    config["stage139_negative_controls"]
  ];
];

If[
  Catch[
    Print["Loaded config from ", ExpandFileName[configPath]];
    Scan[checkSeedConsistency, config["seed_check_sigmas"]];
    checkStage139[];
    Print["\nAll stage-138/139 fixed-point stress checks passed."];
    "ok"
  ] === "ok",
  Null,
  Print["\nStage-138/139 stress harness failed."];
  Exit[1]
];
