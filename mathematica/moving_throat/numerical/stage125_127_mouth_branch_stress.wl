(* Numerical stress harness for the explicit Family-1 mouth branch (Stages 125-127). *)

ClearAll["Global`*"];
$HistoryLength = 0;

rootDir = DirectoryName[$InputFileName];
configPath = FileNameJoin[{
   rootDir, "..", "..", "..", "scripts", "moving_throat", "numerical",
   "stage125_127_mouth_branch_samples.json"
}];
config = Import[ExpandFileName[configPath], "RawJSON"];

If[config["schema"] =!= "moving_throat_numerical_stage125_127_v1",
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
tolerances = config["tolerances"];

rF1 = N[const["rF1"]];
gMinusExpected = N[const["g_minus_expected"]];
gPlusExpected = N[const["g_plus_expected"]];
PiStarExpected = N[const["Pi_star_expected"]];
PiMatchExpected = N[const["Pi_match_expected"]];
THatStarExpected = N[const["T_hat_star_expected"]];
THatMatchExpected = N[const["T_hat_match_expected"]];
RInfExpected = N[const["R_inf_expected"]];
denomInfExpected = N[const["denom_inf_expected"]];
tHatSqrtPiLimitExpected = N[const["that_sqrtPi_limit_expected"]];

quadratureTol = N[tolerances["quadrature_tol"]];
normalizationTol = N[tolerances["normalization_tol"]];
residualTol = N[tolerances["residual_tol"]];
rootTol = N[tolerances["root_tol"]];
ratioTol = N[tolerances["ratio_tol"]];
bisectionIterations = config["bisection_iterations"];

n = config["quadrature_points"];
xGrid = Subdivide[0.0, 1.0, n - 1];
dx = xGrid[[2]] - xGrid[[1]];
weights = ConstantArray[dx, n];
weights[[1]] = dx/2;
weights[[-1]] = dx/2;
cGrid = Cos[Pi xGrid/2];
KqGrid = Cosh[(Pi/2) (1 - xGrid)]/Cosh[Pi/2];

trapz[v_List] := Total[v weights];
sigmaProfile[piValue_] := Module[{raw},
  raw = piValue Exp[-piValue xGrid]/(1 - Exp[-piValue]);
  raw / trapz[raw]
];
gFormula[piValue_] := 2 piValue (2 piValue Exp[piValue] + Pi)/((4 piValue^2 + Pi^2) (Exp[piValue] - 1));
SFormula[piValue_] := Module[{num, denom},
  num = piValue ((Pi/2) Tanh[Pi/2] + piValue (Exp[-piValue]/Cosh[Pi/2] - 1));
  denom = (1 - Exp[-piValue]) (((Pi/2)^2) - piValue^2);
  num/denom
];
RFromG[gValue_] := (gValue - rF1)^2/(1 + rF1^2);
sigma0Formula[piValue_] := Module[{gv, sv},
  gv = gFormula[piValue];
  sv = SFormula[piValue];
  piValue/(1 - RFromG[gv] sv)
];
tHatFormula[piValue_] := Sqrt[(9 sigma0Formula[piValue])/20];

bisectionRoot[target_, loIn_, hiIn_] := Module[{lo, hi, flo, fhi, mid, fm},
  lo = N[loIn];
  hi = N[hiIn];
  flo = gFormula[lo] - target;
  fhi = gFormula[hi] - target;
  require["root bracket sign change on [" <> fmt[lo] <> ", " <> fmt[hi] <> "]",
    flo == 0 || fhi == 0 || flo fhi < 0,
    "f(lo)=" <> fmt[flo] <> ", f(hi)=" <> fmt[fhi]];
  Do[
    mid = (lo + hi)/2;
    fm = gFormula[mid] - target;
    If[flo == 0, Return[lo]];
    If[fhi == 0, Return[hi]];
    If[flo fm <= 0,
      hi = mid; fhi = fm,
      lo = mid; flo = fm
    ],
    {bisectionIterations}
  ];
  (lo + hi)/2
];

stage125Samples[] := Module[{samplePoints},
  Print["=== Stages 125-126: explicit positive-mouth branch samples ==="];
  samplePoints = config["sample_points"];
  Scan[
    Function[{case},
      Module[{piValue, sigma, norm, minSigma, gQuad, SQuad, gExact, SExact, rQuad, sigma0, denomExact, tHat, residual},
        piValue = N[case["Pi"]];
        sigma = sigmaProfile[piValue];
        norm = trapz[sigma];
        minSigma = Min[sigma];
        gQuad = trapz[sigma cGrid];
        SQuad = trapz[sigma KqGrid];
        gExact = gFormula[piValue];
        SExact = SFormula[piValue];
        rQuad = RFromG[gQuad];
        sigma0 = sigma0Formula[piValue];
        denomExact = 1 - RFromG[gExact] SExact;
        tHat = tHatFormula[piValue];
        residual = piValue - sigma0 (1 - rQuad SQuad);

        Print[""];
        Print[case["name"], " (", case["kind"], "): Pi=", fmt[piValue]];
        If[KeyExistsQ[case, "assumptions"],
          Print["assumptions:"];
          Scan[Function[{item}, Print["  - ", item]], case["assumptions"]];
        ];
        Print[
          "  g_quad=", fmt[gQuad], ", g_exact=", fmt[gExact],
          ", S_quad=", fmt[SQuad], ", S_exact=", fmt[SExact]
        ];
        Print[
          "  R=", fmt[rQuad], ", denom=", fmt[denomExact],
          ", Sigma0=", fmt[sigma0], ", T_hat=", fmt[tHat]
        ];

        require[case["name"] <> " positive source normalization",
          Abs[norm - 1] <= normalizationTol,
          "integral=" <> fmt[norm]];
        require[case["name"] <> " source stays nonnegative",
          minSigma >= -10^-15,
          "min sigma=" <> fmt[minSigma]];
        require[case["name"] <> " g quadrature matches formula",
          nearQ[gQuad, gExact, quadratureTol],
          "quadrature=" <> fmt[gQuad] <> ", formula=" <> fmt[gExact]];
        require[case["name"] <> " S quadrature matches formula",
          nearQ[SQuad, SExact, quadratureTol],
          "quadrature=" <> fmt[SQuad] <> ", formula=" <> fmt[SExact]];
        require[case["name"] <> " self-consistent residual",
          Abs[residual] <= residualTol,
          "residual=" <> fmt[residual]];
        require[case["name"] <> " stays in 2/pi < g < 1",
          2/Pi < gExact < 1,
          "2/pi=" <> fmt[2/Pi] <> ", g=" <> fmt[gExact]];
        require[case["name"] <> " denominator stays positive",
          denomExact > 0.8,
          "1-R_q S_q=" <> fmt[denomExact]];
      ]
    ],
    samplePoints
  ];

  require["canonical T_hat matches audited value",
    nearQ[tHatFormula[PiStarExpected], THatStarExpected, rootTol],
    "T_hat(Pi_*)=" <> fmt[tHatFormula[PiStarExpected]] <> ", expected=" <> fmt[THatStarExpected]];
  require["derivative-match T_hat matches audited value",
    nearQ[tHatFormula[PiMatchExpected], THatMatchExpected, rootTol],
    "T_hat(Pi_match)=" <> fmt[tHatFormula[PiMatchExpected]] <> ", expected=" <> fmt[THatMatchExpected]];
];

stage127BranchSelection[] := Module[{scanValues, diffs, roots = <||>, target, blockRoots, refRoot, exclusionValues},
  Print[""];
  Print["=== Stage 127: branch selection and finite-bias exclusion ==="];
  scanValues = gFormula /@ N[config["monotone_scan_points"]];
  Scan[
    Function[{pair},
      Print["scan Pi=", fmt[pair[[1]]], " -> g(Pi)=", fmt[pair[[2]]]]
    ],
    Transpose[{N[config["monotone_scan_points"]], scanValues}]
  ];

  diffs = Differences[scanValues];
  require["g(Pi) increases on the explicit branch scan",
    AllTrue[diffs, # > 0 &],
    "min diff=" <> fmt[Min[diffs]]];

  Scan[
    Function[{block},
      target = If[block["target_key"] === "pi_over_4", N[Pi/4], N[const[block["target_key"]]]];
      blockRoots = Table[
        Module[{lo, hi, root},
          lo = N[pair[[1]]];
          hi = N[pair[[2]]];
          root = bisectionRoot[target, lo, hi];
          Print[block["name"], " root from [", fmt[lo], ", ", fmt[hi], "] = ", fmt[root]];
          root
        ],
        {pair, block["brackets"]}
      ];
      refRoot = First[blockRoots];
      Do[
        require[block["name"] <> " bracket " <> ToString[idx] <> " agrees with bracket 0",
          nearQ[blockRoots[[idx + 1]], refRoot, rootTol],
          "root=" <> fmt[blockRoots[[idx + 1]]] <> ", ref=" <> fmt[refRoot]],
        {idx, 1, Length[blockRoots] - 1}
      ];
      require[block["name"] <> " root matches expected value",
        nearQ[refRoot, N[block["expected_pi"]], rootTol],
        "root=" <> fmt[refRoot] <> ", expected=" <> fmt[N[block["expected_pi"]]]];
      roots[block["name"]] = refRoot;
    ],
    config["root_checks"]
  ];

  require["derivative-match point lies above canonical point",
    roots["derivative_match_branch"] > roots["lower_compensated_branch"],
    "Pi_match=" <> fmt[roots["derivative_match_branch"]] <> ", Pi_*=" <> fmt[roots["lower_compensated_branch"]]];

  exclusionValues = gFormula /@ N[config["finite_exclusion_points"]];
  require["finite explicit branch never reaches g=1",
    AllTrue[exclusionValues, # < 1 &],
    "max g(Pi) on finite window=" <> fmt[Max[exclusionValues]]];
  require["finite explicit branch never reaches g_+",
    AllTrue[exclusionValues, # < gPlusExpected &],
    "max g(Pi)=" <> fmt[Max[exclusionValues]] <> ", g_+=" <> fmt[gPlusExpected]];
  require["upper compensated branch remains outside the positive-source range",
    gPlusExpected > 1,
    "g_+=" <> fmt[gPlusExpected]];
];

stage126Asymptotics[] := Module[{denomErrors = {}, ratioErrors = {}, gTail = {}, samples, gValue, sValue, rValue, denom, ratio},
  Print[""];
  Print["=== Stage 126: singular-limit asymptotics ==="];
  samples = N[config["asymptotic_samples"]];
  Scan[
    Function[{piValue},
      gValue = gFormula[piValue];
      sValue = SFormula[piValue];
      rValue = RFromG[gValue];
      denom = 1 - rValue sValue;
      ratio = tHatFormula[piValue]/Sqrt[piValue];
      AppendTo[gTail, gValue];
      AppendTo[denomErrors, Abs[denom - denomInfExpected]];
      AppendTo[ratioErrors, Abs[ratio - tHatSqrtPiLimitExpected]];
      Print[
        "Pi=", fmt[piValue], ": 1-g=", fmt[1 - gValue],
        ", denom=", fmt[denom], ", T_hat/sqrt(Pi)=", fmt[ratio]
      ];
      require["Pi=" <> fmt[piValue] <> " keeps g(Pi)<1",
        gValue < 1,
        "g(Pi)=" <> fmt[gValue]];
      require["Pi=" <> fmt[piValue] <> " asymptotic T_hat ratio stays near limit",
        Abs[ratio - tHatSqrtPiLimitExpected] <= ratioTol,
        "ratio=" <> fmt[ratio] <> ", limit=" <> fmt[tHatSqrtPiLimitExpected]];
    ],
    samples
  ];

  require["1-g(Pi) decreases toward the singular limit",
    AllTrue[Differences[1 - # & /@ gTail], # < 0 &],
    "tail gaps=" <> StringRiffle[fmt /@ (1 - # & /@ gTail), ", "]];
  require["denominator approaches positive large-Pi limit",
    AllTrue[Differences[denomErrors], # < 0 &],
    "errors=" <> StringRiffle[fmt /@ denomErrors, ", "]];
  require["T_hat/sqrt(Pi) approaches the predicted asymptotic coefficient",
    AllTrue[Differences[ratioErrors], # < 0 &],
    "errors=" <> StringRiffle[fmt /@ ratioErrors, ", "]];
  require["large-Pi denominator limit matches Stage 126",
    nearQ[1 - RInfExpected, denomInfExpected, rootTol],
    "1-R_inf=" <> fmt[1 - RInfExpected] <> ", expected=" <> fmt[denomInfExpected]];
];

If[
  Catch[
    Print["Loaded config from ", ExpandFileName[configPath]];
    stage125Samples[];
    stage127BranchSelection[];
    stage126Asymptotics[];
    Print[""];
    Print["All stage-125/127 mouth-branch stress checks passed."];
    "ok"
  ] === "ok",
  Null,
  Print[""];
  Print["Stage-125/127 mouth-branch stress harness failed."];
  Exit[1]
];
