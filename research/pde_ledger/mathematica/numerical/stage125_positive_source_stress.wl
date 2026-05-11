(* Positive-source numerical stress harness for Stage 125. *)

ClearAll["Global`*"];
$HistoryLength = 0;

rootDir = DirectoryName[$InputFileName];
configPath = FileNameJoin[{
   rootDir, "..", "..", "scripts", "numerical",
   "stage125_positive_source_samples.json"
}];
config = Import[ExpandFileName[configPath], "RawJSON"];

If[config["schema"] =!= "moving_throat_numerical_stage125_v1",
  Print["Unexpected config schema."];
  Exit[1];
];

fmt[x_] := ToString @ NumberForm[N[x, 14], {Infinity, 12}, ExponentFunction -> (Null &)];

require[label_, condition_, detail_] := Module[{status},
  status = If[TrueQ[condition], "PASS", "FAIL"];
  Print["[", status, "] ", label, ": ", detail];
  If[!TrueQ[condition], Throw[$Failed]]
];

stage125Block[] := Module[
  {n, iters, tol, gMinus, gPlus, xGrid, dx, weights, kernel, trapz, normalizedProfile, gValue, values = <||>, bracket, lo, hi, flo, fhi, mid, fm, betaStar, gStar},
  Print["=== Stage 125: positive local mouth-source stress ==="];
  n = config["quadrature_points"];
  iters = config["bisection_iterations"];
  tol = config["tolerances"];
  gMinus = N[config["g_minus"]];
  gPlus = N[config["g_plus"]];
  xGrid = Subdivide[0.0, 1.0, n - 1];
  dx = xGrid[[2]] - xGrid[[1]];
  weights = ConstantArray[dx, n];
  weights[[1]] = dx/2;
  weights[[-1]] = dx/2;
  kernel = Cos[Pi xGrid/2];

  trapz[v_List] := Total[v weights];

  normalizedProfile[profile_] := Module[{kind, raw, p, beta},
    kind = profile["kind"];
    raw = Switch[kind,
      "uniform", ConstantArray[1.0, n],
      "left_linear", 2.0 (1.0 - xGrid),
      "left_beta", p = N[profile["p"]]; (p + 1.0) (1.0 - xGrid)^p,
      "right_beta", p = N[profile["p"]]; (p + 1.0) xGrid^p,
      "exp_left", beta = N[profile["beta"]]; beta Exp[-beta xGrid]/(1.0 - Exp[-beta]),
      _, Print["Unknown profile kind: ", kind]; Throw[$Failed]
    ];
    raw / trapz[raw]
  ];

  gValue[profile_List] := trapz[profile kernel];

  Scan[
    Function[{profile},
      Module[{sigma, norm, gSigma},
        sigma = normalizedProfile[profile];
        norm = trapz[sigma];
        gSigma = gValue[sigma];
        values[profile["name"]] = gSigma;
        Print[""];
        Print[profile["name"], " (", profile["kind"], "): g[sigma]=", fmt[gSigma]];
        require[profile["name"] <> " stays nonnegative",
          Min[sigma] >= -10^-12,
          "min sigma=" <> fmt[Min[sigma]]];
        require[profile["name"] <> " normalization",
          Abs[norm - 1] <= N[tol["normalization_tol"]],
          "integral=" <> fmt[norm]];
        require[profile["name"] <> " moment stays in [0,1]",
          -N[tol["moment_tol"]] <= gSigma <= 1 + N[tol["moment_tol"]],
          "g[sigma]=" <> fmt[gSigma]];
      ]
    ],
    config["profiles"]
  ];

  require["g_- lies in (0,1)", 0 < gMinus < 1, "g_-=" <> fmt[gMinus]];
  require["g_+ lies above 1", gPlus > 1, "g_+=" <> fmt[gPlus]];
  require["uniform < g_- < left_linear",
    values["uniform"] < gMinus < values["left_linear"],
    "uniform=" <> fmt[values["uniform"]] <> ", g_-=" <> fmt[gMinus] <> ", left_linear=" <> fmt[values["left_linear"]]];
  require["right-localized profiles drive g toward 0",
    values["right_beta_20"] < values["right_beta_8"] < values["uniform"],
    "right_beta_20=" <> fmt[values["right_beta_20"]] <> ", right_beta_8=" <> fmt[values["right_beta_8"]] <> ", uniform=" <> fmt[values["uniform"]]];
  require["left-localized profiles drive g toward 1",
    values["left_beta_20"] > values["left_beta_8"] > values["left_linear"] > gMinus,
    "left_beta_20=" <> fmt[values["left_beta_20"]] <> ", left_beta_8=" <> fmt[values["left_beta_8"]] <> ", left_linear=" <> fmt[values["left_linear"]] <> ", g_-=" <> fmt[gMinus]];

  bracket = config["gminus_realization"]["beta_bracket"];
  lo = N[bracket[[1]]];
  hi = N[bracket[[2]]];

  expProfile[beta_] := Module[{raw},
    raw = beta Exp[-beta xGrid]/(1.0 - Exp[-beta]);
    raw / trapz[raw]
  ];
  f[beta_] := gValue[expProfile[beta]] - gMinus;

  flo = f[lo];
  fhi = f[hi];
  require["g_- realization bracket",
    flo <= 0 && fhi >= 0,
    "f(lo)=" <> fmt[flo] <> ", f(hi)=" <> fmt[fhi]];
  Do[
    mid = (lo + hi)/2;
    fm = f[mid];
    If[fm <= 0, lo = mid, hi = mid],
    {iters}
  ];
  betaStar = (lo + hi)/2;
  gStar = gValue[expProfile[betaStar]];
  require["exp_left family realizes g_-",
    Abs[gStar - gMinus] <= N[tol["root_tol"]],
    "beta*=" <> fmt[betaStar] <> ", g(beta*)=" <> fmt[gStar] <> ", g_-=" <> fmt[gMinus]];

  Print[""];
  Print["Solved exp_left beta* = ", fmt[betaStar], " with g[beta*] = ", fmt[gStar]];
  Print[""];
  Print["All stage-125 positive-source stress checks passed."];
];

If[
  Catch[
    Print["Loaded config from ", ExpandFileName[configPath]];
    stage125Block[];
    "ok"
  ] === "ok",
  Null,
  Print[""];
  Print["Stage 125 positive-source stress harness failed."];
  Exit[1]
];
