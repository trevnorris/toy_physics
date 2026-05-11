(* Foundational numerical spot checks for Stages 041 and 153. *)

ClearAll["Global`*"];
$HistoryLength = 0;

rootDir = DirectoryName[$InputFileName];
configPath = FileNameJoin[{
   rootDir, "..", "..", "scripts", "numerical",
   "stage058_170_foundational_samples.json"
}];
config = Import[ExpandFileName[configPath], "RawJSON"];

If[config["schema"] =!= "moving_throat_numerical_stage058_170_v1",
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

trapz[y_List, x_List] := Module[{dx},
  dx = Differences[x];
  Total[dx*(Most[y] + Rest[y])/2]
];

stage058Block[] := Module[
  {block, tol, quadraturePoints, bisectionIterations, xGrid, w, kernel, sigma, deltaFormula, deltaQuadrature, delta0, deltaInf, bisect},
  Print[""];
  Print["=== Stage 058: coupled support/source fixed point ==="];
  block = config["stage058"];
  tol = block["tolerances"];
  quadraturePoints = block["quadrature_points"];
  bisectionIterations = block["bisection_iterations"];
  xGrid = Subdivide[0.0, 1.0, quadraturePoints - 1];

  w[alpha_, eta_] := alpha Sinh[alpha] + eta Cosh[alpha];
  kernel[alpha_, eta_, xvals_List] :=
    (Cosh[alpha xvals] + (eta/alpha) Sinh[alpha xvals] - Cosh[alpha (1 - xvals)]) / w[alpha, eta];
  sigma[pe_, xvals_List] := pe Exp[pe xvals] / (Exp[pe] - 1);
  deltaFormula[pe_, alpha_, eta_] := Module[{ic, is},
    ic = (Exp[pe] (pe Cosh[alpha] - alpha Sinh[alpha]) - pe)/(pe^2 - alpha^2);
    is = (Exp[pe] (pe Sinh[alpha] - alpha Cosh[alpha]) + alpha)/(pe^2 - alpha^2);
    pe/(Exp[pe] - 1) (((1 - Cosh[alpha]) ic + (eta/alpha + Sinh[alpha]) is)/w[alpha, eta])
  ];
  deltaQuadrature[pe_, alpha_, eta_] := trapz[kernel[alpha, eta, xGrid] sigma[pe, xGrid], xGrid];
  delta0[alpha_, eta_] := eta (Cosh[alpha] - 1)/(alpha^2 w[alpha, eta]);
  deltaInf[alpha_, eta_] := (Cosh[alpha] + (eta/alpha) Sinh[alpha] - 1)/w[alpha, eta];
  bisect[alpha_, eta_, xi_] := Module[{lo, hi, flo, fhi, mid, fm, root, residual},
    lo = xi delta0[alpha, eta];
    hi = xi deltaInf[alpha, eta];
    flo = lo - xi deltaFormula[lo, alpha, eta];
    fhi = hi - xi deltaFormula[hi, alpha, eta];
    require["stage058 bracket sign", flo <= tol["bracket_tol"] && fhi >= -tol["bracket_tol"],
      "flo=" <> fmt[flo] <> ", fhi=" <> fmt[fhi]];
    Do[
      mid = (lo + hi)/2;
      fm = mid - xi deltaFormula[mid, alpha, eta];
      If[fm <= 0, lo = mid, hi = mid],
      {bisectionIterations}
    ];
    root = (lo + hi)/2;
    residual = root - xi deltaFormula[root, alpha, eta];
    <|"root" -> root, "PeLo" -> xi delta0[alpha, eta], "PeHi" -> xi deltaInf[alpha, eta], "residual" -> residual|>
  ];

  Scan[
    Function[{case},
      Module[{alpha, eta, xi, d0Formula, d0Quad, dInf, rootData, root, deltaRootFormula, deltaRootQuad},
        alpha = N[case["alpha"]];
        eta = N[case["eta"]];
        xi = N[case["Xi"]];
        Print[""];
        Print[case["name"], " (", case["kind"], "): alpha=", fmt[alpha], ", eta=", fmt[eta], ", Xi=", fmt[xi]];
        d0Formula = delta0[alpha, eta];
        d0Quad = trapz[kernel[alpha, eta, xGrid], xGrid];
        dInf = deltaInf[alpha, eta];
        rootData = bisect[alpha, eta, xi];
        root = rootData["root"];
        deltaRootFormula = deltaFormula[root, alpha, eta];
        deltaRootQuad = deltaQuadrature[root, alpha, eta];
        require[case["name"] <> " Delta0 quadrature check",
          nearQ[d0Quad, d0Formula, tol["quadrature_tol"]],
          "quadrature=" <> fmt[d0Quad] <> ", formula=" <> fmt[d0Formula]];
        require[case["name"] <> " Delta(root) quadrature check",
          nearQ[deltaRootQuad, deltaRootFormula, tol["quadrature_tol"]],
          "quadrature=" <> fmt[deltaRootQuad] <> ", formula=" <> fmt[deltaRootFormula]];
        require[case["name"] <> " fixed point residual",
          Abs[rootData["residual"]] <= tol["residual_tol"],
          "residual=" <> fmt[rootData["residual"]]];
        require[case["name"] <> " Pe bracket",
          rootData["PeLo"] - tol["bracket_tol"] <= root <= rootData["PeHi"] + tol["bracket_tol"],
          "Pe_lo=" <> fmt[rootData["PeLo"]] <> ", root=" <> fmt[root] <> ", Pe_hi=" <> fmt[rootData["PeHi"]]];
        Print[
          "  Delta0=", fmt[d0Formula], ", Delta_inf=", fmt[dInf],
          ", Pe_*=", fmt[root], ", Xi Delta(Pe_*)=", fmt[xi deltaRootFormula]
        ];
      ]
    ],
    block["cases"]
  ];
];

stage170Block[] := Module[
  {block, tol, eps, baseline, d0, n0, sigma, p0, d2, d4, finiteDerivatives, expectedMap, pairResults = <||>, first, second},
  Print[""];
  Print["=== Stage 170: grouped outlet map ==="];
  block = config["stage170"];
  tol = block["tolerances"];
  eps = N[block["epsilon_step"]];
  baseline = block["baseline"];
  d0 = N[baseline["D0"]];
  n0 = N[baseline["N0"]];
  sigma = N[baseline["sigma"]];
  p0 = n0/d0;
  d2 = -d0/9;
  d4 = -d0/27;

  finiteDerivatives[dD0_, dD2_, dD4_, dN0_] := Module[{u2Full, u4Full, p0Full, du2, du4, dp0},
    u2Full = -(d2 + eps dD2)/(d0 + eps dD0);
    u4Full = ((d2 + eps dD2)^2 - (d0 + eps dD0) (d4 + eps dD4))/(d0 + eps dD0)^2;
    p0Full = (n0 + eps dN0)/(d0 + eps dD0);
    du2 = (u2Full - 1/9)/eps;
    du4 = (u4Full - 4/81)/eps;
    dp0 = (p0Full - p0)/eps;
    <|
      "du2" -> du2,
      "du4" -> du4,
      "dP0" -> dp0,
      "dkappa" -> -3 (1 - sigma) du2/sigma,
      "dgamma" -> -(1 - sigma) dp0/(9 sigma p0)
    |>
  ];

  expectedMap[dD0_, dD2_, dN0_] := Module[{kA, gA},
    kA = dD2 + dD0/9;
    gA = dN0 - p0 dD0;
    <|
      "K_A" -> kA,
      "G_A" -> gA,
      "dkappa" -> 3 (1 - sigma) kA/(sigma d0),
      "dgamma" -> -(1 - sigma) gA/(9 sigma n0)
    |>
  ];

  Scan[
    Function[{case},
      Module[{dD0, dD2, dN0, dD4, exact, expected},
        dD0 = N[case["dD0"]];
        dD2 = N[case["dD2"]];
        dN0 = N[case["dN0"]];
        dD4 = 2 dD2/3 + dD0/27;
        exact = finiteDerivatives[dD0, dD2, dD4, dN0];
        expected = expectedMap[dD0, dD2, dN0];
        Print[""];
        Print[case["name"], ": dD0=", fmt[dD0], ", dD2=", fmt[dD2], ", dD4=", fmt[dD4], ", dN0=", fmt[dN0]];
        require[case["name"] <> " direct even consistency",
          nearQ[exact["du4"], (8/9) exact["du2"], tol["pair_tol"]],
          "du4=" <> fmt[exact["du4"]] <> ", 8/9 du2=" <> fmt[(8/9) exact["du2"]]];
        require[case["name"] <> " delta kappa map",
          nearQ[exact["dkappa"], expected["dkappa"], tol["linear_tol"]],
          "direct=" <> fmt[exact["dkappa"]] <> ", expected=" <> fmt[expected["dkappa"]]];
        require[case["name"] <> " delta gamma map",
          nearQ[exact["dgamma"], expected["dgamma"], tol["linear_tol"]],
          "direct=" <> fmt[exact["dgamma"]] <> ", expected=" <> fmt[expected["dgamma"]]];
        pairResults[case["name"]] = Join[exact, expected];
        Print[
          "  K_A=", fmt[expected["K_A"]], ", G_A=", fmt[expected["G_A"]],
          ", delta_kappa=", fmt[exact["dkappa"]], ", delta_gamma=", fmt[exact["dgamma"]]
        ];
      ]
    ],
    block["consistent_pairs"]
  ];

  first = pairResults[block["consistent_pairs"][[1]]["name"]];
  second = pairResults[block["consistent_pairs"][[2]]["name"]];
  require["consistent pair K_A collapse",
    nearQ[first["K_A"], second["K_A"], tol["pair_tol"]],
    "K_A(1)=" <> fmt[first["K_A"]] <> ", K_A(2)=" <> fmt[second["K_A"]]];
  require["consistent pair G_A collapse",
    nearQ[first["G_A"], second["G_A"], tol["pair_tol"]],
    "G_A(1)=" <> fmt[first["G_A"]] <> ", G_A(2)=" <> fmt[second["G_A"]]];
  require["consistent pair delta kappa equality",
    nearQ[first["dkappa"], second["dkappa"], tol["pair_tol"]],
    "delta_kappa(1)=" <> fmt[first["dkappa"]] <> ", delta_kappa(2)=" <> fmt[second["dkappa"]]];
  require["consistent pair delta gamma equality",
    nearQ[first["dgamma"], second["dgamma"], tol["pair_tol"]],
    "delta_gamma(1)=" <> fmt[first["dgamma"]] <> ", delta_gamma(2)=" <> fmt[second["dgamma"]]];

  Scan[
    Function[{case},
      Module[{dD0, dD2, dD4, dN0, exact, mismatch},
        dD0 = N[case["dD0"]];
        dD2 = N[case["dD2"]];
        dD4 = N[case["dD4"]];
        dN0 = N[case["dN0"]];
        exact = finiteDerivatives[dD0, dD2, dD4, dN0];
        mismatch = exact["du4"] - (8/9) exact["du2"];
        Print[""];
        Print[case["name"], ": mismatch=", fmt[mismatch]];
        require[case["name"] <> " violates hidden-even consistency",
          Abs[mismatch] >= tol["inconsistency_floor"],
          "du4 - 8/9 du2=" <> fmt[mismatch]];
      ]
    ],
    block["inconsistent_cases"]
  ];
];

If[
  Catch[
    Print["Loaded config from ", ExpandFileName[configPath]];
    stage058Block[];
    stage170Block[];
    Print[""];
    Print["All stage-058/170 foundational stress checks passed."];
    "ok"
  ] === "ok",
  Null,
  Print[""];
  Print["Stage 058/170 foundational stress harness failed."];
  Exit[1]
];
