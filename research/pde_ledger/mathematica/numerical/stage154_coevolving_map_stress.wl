(* Numerical bridge harness for the Stage 154 co-evolving core-mouth map. *)

ClearAll["Global`*"];

rootDir = DirectoryName[$InputFileName];
configPath = FileNameJoin[{
   rootDir, "..", "..", "scripts", "numerical",
   "stage154_coevolving_map_samples.json"
}];
config = Import[ExpandFileName[configPath], "RawJSON"];

If[config["schema"] =!= "moving_throat_numerical_stage154_v1",
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
tol = config["tolerances"];

rF1 = const["rF1"];
gStar = const["g_star"];
sStar = const["S_star"];
piStar = const["Pi_star"];
sigma0Star = const["Sigma0_star"];
tHatStar = const["T_hat_star"];
sigma0Can = const["Sigma0_can"];
sCan = const["S_can"];
piCan = const["Pi_can"];
tHatCan = const["T_hat_can"];
gFpStar = const["g_fp_star"];
sFpStar = const["S_fp_star"];
rFpStar = const["R_fp_star"];
piFpStar = const["Pi_fp_star"];

normalizationTol = tol["normalization_tol"];
positivityTol = tol["positivity_tol"];
momentTol = tol["moment_tol"];
fixedPointResidualTol = tol["fixed_point_residual_tol"];
kernelSlopeTol = tol["kernel_slope_tol"];
phiSlopeTol = tol["phi_slope_tol"];

n = config["grid_points"];
maxIterations = config["max_iterations"];
fixedPointTol = config["fixed_point_tol"];
fdSteps = config["fd_steps"];

xGrid = Subdivide[0.0, 1.0, n - 1];
dx = xGrid[[2]] - xGrid[[1]];
weights = ConstantArray[dx, n];
weights[[1]] = dx/2;
weights[[-1]] = dx/2;
kappa = N[Pi/2];
cGrid = Cos[(Pi*xGrid)/2];
kqGrid = Cosh[kappa*(1 - xGrid)]/Cosh[kappa];
sqrtTerm = Sqrt[1 + rF1^2];
gStarExpected = rF1 - sqrtTerm/2;

normalize[sig_List] := sig/Total[sig*weights];

sigmaVar[lam_?NumericQ] := normalize[(1 - lam) ConstantArray[1.0, n] + lam*(Pi/2)*Cos[(Pi*xGrid)/2]];

seedProfile[kind_String, lam_: Missing["NotAvailable"]] := Switch[kind,
  "canonical_exponential",
    normalize[piStar*Exp[-piStar*xGrid]],
  "uniform",
    normalize[ConstantArray[1.0, n]],
  "derivative_match",
    normalize[(Pi/2)*Cos[(Pi*xGrid)/2]],
  "convex_blend",
    sigmaVar[lam],
  _,
    Print["Unknown seed kind: ", kind];
    Throw[$Failed]
];

canonicalProfile = seedProfile["canonical_exponential"];

gMoment[sig_List] := Total[cGrid*sig*weights];
sMoment[sig_List] := Total[kqGrid*sig*weights];
rMoment[sig_List] := ((gMoment[sig] - rF1)^2)/(1 + rF1^2);

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

phiGrid[sig_List, sigma0_?NumericQ] := sigma0*(tsOperator[sig] - rMoment[sig]*tqOperator[sig]);
nextSigma[sig_List, sigma0_?NumericQ] := Module[{ph, phShift},
  ph = phiGrid[sig, sigma0];
  phShift = ph - Min[ph];
  normalize[Exp[-phShift]]
];

solveFixedPoint[sigma0_?NumericQ, initial_List] := Module[{sig, sigNew, err = Infinity, iter},
  sig = initial;
  For[iter = 1, iter <= maxIterations, iter++,
    sigNew = nextSigma[sig, sigma0];
    err = Max[Abs[sigNew - sig]];
    sig = sigNew;
    If[err < fixedPointTol,
      Return[<|"sigma" -> sig, "iterations" -> iter, "error" -> err|>]
    ];
  ];
  Print["[FAIL] fixed-point convergence: sigma0=", fmt[sigma0], ", error=", fmt[err]];
  Throw[$Failed]
];

tsAt[sig_List, x_?NumericQ] := Total[MapThread[Min[x, #1] #2 #3 &, {xGrid, sig, weights}]];

tqAt[sig_List, x_?NumericQ] := Module[{kernel},
  kernel = Map[
    (Sinh[kappa*Min[x, #]]*Cosh[kappa*(1 - Max[x, #])])/(kappa*Cosh[kappa]) &,
    xGrid
  ];
  Total[kernel*sig*weights]
];

phiAt[sig_List, sigma0_?NumericQ, x_?NumericQ] := sigma0*(tsAt[sig, x] - rMoment[sig]*tqAt[sig, x]);

forwardSlope[fn_, h_?NumericQ] := (4 fn[h] - fn[2 h])/(2 h);

checkSlopes[sig_List, sigma0_?NumericQ, prefix_String] := Module[
  {sValue, rValue, piExpected, phiSlopes = {}, tsSlope, tqSlope, phiSlope, h, avgPhi},
  sValue = sMoment[sig];
  rValue = rMoment[sig];
  piExpected = sigma0*(1 - rValue*sValue);
  Do[
    tsSlope = forwardSlope[Function[z, tsAt[sig, z]], h];
    tqSlope = forwardSlope[Function[z, tqAt[sig, z]], h];
    phiSlope = forwardSlope[Function[z, phiAt[sig, sigma0, z]], h];
    AppendTo[phiSlopes, phiSlope];
    require[prefix <> " T_s'(0) at h=" <> fmt[h], nearQ[tsSlope, 1.0, kernelSlopeTol],
      "slope=" <> fmt[tsSlope]];
    require[prefix <> " T_q'(0)=S at h=" <> fmt[h], nearQ[tqSlope, sValue, kernelSlopeTol],
      "slope=" <> fmt[tqSlope] <> ", S=" <> fmt[sValue]];
    require[prefix <> " Phi'(0)=Sigma0(1-RS) at h=" <> fmt[h], nearQ[phiSlope, piExpected, phiSlopeTol],
      "slope=" <> fmt[phiSlope] <> ", expected=" <> fmt[piExpected]];
    ,
    {h, fdSteps}
  ];
  avgPhi = Mean[phiSlopes];
  require[prefix <> " forward-slope self-consistency", nearQ[avgPhi, piExpected, phiSlopeTol],
    "avg slope=" <> fmt[avgPhi] <> ", expected=" <> fmt[piExpected]];
  avgPhi
];

representativeProfile[case_] := Module[{kind, seedKind, initial},
  kind = case["kind"];
  Switch[kind,
    "explicit_canonical",
      canonicalProfile,
    "derivative_match",
      seedProfile["derivative_match"],
    "convex_blend",
      seedProfile["convex_blend", case["lambda"]],
    "fixed_point",
      seedKind = Lookup[case, "seed_kind", "canonical_exponential"];
      initial = seedProfile[seedKind];
      solveFixedPoint[case["sigma0"], initial]["sigma"],
    _,
      Print["Unknown representative profile kind: ", kind];
      Throw[$Failed]
  ]
];

checkRepresentativeProfiles[] := Module[
  {sigma0, sigma, gValue, sValue, rValue, norm, minSigma, piDirect, residual},
  Print["=== Stage 154: representative profile checks ==="];
  require["loaded g_* matches the exact compensation formula", nearQ[gStar, gStarExpected, momentTol],
    "g_*=" <> fmt[gStar] <> ", expected=" <> fmt[gStarExpected]];
  Scan[
    Function[case,
      sigma0 = case["sigma0"];
      sigma = representativeProfile[case];
      gValue = gMoment[sigma];
      sValue = sMoment[sigma];
      rValue = rMoment[sigma];
      norm = Total[sigma*weights];
      minSigma = Min[sigma];
      Print[""];
      Print["=== ", case["name"], " (", case["kind"], ") ==="];
      Print["sigma0 = ", fmt[sigma0]];
      If[KeyExistsQ[case, "assumptions"],
        Print["assumptions:"];
        Scan[Print["  - ", #] &, case["assumptions"]];
      ];
      piDirect = checkSlopes[sigma, sigma0, case["name"]];
      Print[
        "g=", fmt[gValue], ", S=", fmt[sValue], ", R=", fmt[rValue],
        ", Pi_direct=", fmt[piDirect], ", T_hat=", fmt[Sqrt[9 sigma0/20]]
      ];
      require[case["name"] <> " normalization", Abs[norm - 1.0] <= normalizationTol,
        "integral=" <> fmt[norm]];
      require[case["name"] <> " positivity", minSigma >= -positivityTol,
        "min sigma=" <> fmt[minSigma]];
      require[case["name"] <> " cosine moment stays in the positive-source range",
        -positivityTol <= gValue <= 1.0 + positivityTol,
        "g=" <> fmt[gValue]];

      If[KeyExistsQ[case, "expected_g"],
        require[case["name"] <> " expected g", nearQ[gValue, case["expected_g"], momentTol],
          "g=" <> fmt[gValue] <> ", expected=" <> fmt[case["expected_g"]]]
      ];
      If[KeyExistsQ[case, "expected_S"],
        require[case["name"] <> " expected S", nearQ[sValue, case["expected_S"], momentTol],
          "S=" <> fmt[sValue] <> ", expected=" <> fmt[case["expected_S"]]]
      ];
      If[KeyExistsQ[case, "expected_R"],
        require[case["name"] <> " expected R", nearQ[rValue, case["expected_R"], momentTol],
          "R=" <> fmt[rValue] <> ", expected=" <> fmt[case["expected_R"]]]
      ];
      If[KeyExistsQ[case, "expected_Pi"],
        require[case["name"] <> " expected Pi", nearQ[piDirect, case["expected_Pi"], phiSlopeTol],
          "Pi=" <> fmt[piDirect] <> ", expected=" <> fmt[case["expected_Pi"]]]
      ];
      If[TrueQ @ Lookup[case, "expect_compensated", False],
        require[case["name"] <> " hits g_*", nearQ[gValue, gStar, momentTol],
          "g=" <> fmt[gValue] <> ", g_*=" <> fmt[gStar]];
        require[case["name"] <> " hits R=1/4", nearQ[rValue, 1/4, momentTol],
          "R=" <> fmt[rValue]];
      ];
      If[case["kind"] === "fixed_point",
        residual = Max[Abs[nextSigma[sigma, sigma0] - sigma]];
        require[case["name"] <> " fixed-point residual", residual <= fixedPointResidualTol,
          "L_inf residual=" <> fmt[residual]];
      ];
    ],
    config["representative_profiles"]
  ];
];

checkTransport[] := Module[
  {sigmaStar, piStarDirect, sigmaVarCase, lam, sigma0, gVar, sVar, dgRate, dSRate, dRRate, dPiRate,
   eps, sigmaEps, norm, minSigma, gEps, sEps, rEps, piEps, deltaPiExact, deltaPiLinear, deltaPiRel},
  Print[""];
  Print["=== Stage 154: first-order transport from the canonical branch ==="];
  sigmaStar = canonicalProfile;
  piStarDirect = checkSlopes[sigmaStar, sigma0Star, "canonical_transport_base"];
  require["canonical transport base Pi matches Pi_*", nearQ[piStarDirect, piStar, phiSlopeTol],
    "Pi_direct=" <> fmt[piStarDirect] <> ", Pi_*=" <> fmt[piStar]];

  Scan[
    Function[case,
      lam = case["lambda"];
      sigma0 = case["sigma0"];
      sigmaVarCase = sigmaVar[lam];
      gVar = gMoment[sigmaVarCase];
      sVar = sMoment[sigmaVarCase];
      dgRate = gVar - gStar;
      dSRate = sVar - sStar;
      dRRate = -dgRate/sqrtTerm;
      dPiRate = -sigma0*(1/4 dSRate + sStar*dRRate);

      Print[""];
      Print["=== ", case["name"], " (", case["kind"], ") ==="];
      Print["lambda = ", fmt[lam]];
      Print[
        "direction rates: dg/dε=", fmt[dgRate], ", dS/dε=", fmt[dSRate],
        ", dR/dε=", fmt[dRRate], ", dPi/dε=", fmt[dPiRate]
      ];

      Scan[
        Function[sample,
          eps = sample["epsilon"];
          sigmaEps = normalize[(1 - eps) sigmaStar + eps sigmaVarCase];
          norm = Total[sigmaEps*weights];
          minSigma = Min[sigmaEps];
          gEps = gMoment[sigmaEps];
          sEps = sMoment[sigmaEps];
          rEps = rMoment[sigmaEps];
          piEps = checkSlopes[sigmaEps, sigma0, case["name"] <> " ε=" <> fmt[eps]];
          deltaPiExact = piEps - piStarDirect;
          deltaPiLinear = eps*dPiRate;
          deltaPiRel = Abs[deltaPiExact - deltaPiLinear]/Abs[deltaPiLinear];

          Print["epsilon = ", fmt[eps]];
          Print[
            "  g=", fmt[gEps], ", S=", fmt[sEps], ", R=", fmt[rEps],
            ", deltaPi_exact=", fmt[deltaPiExact], ", deltaPi_linear=", fmt[deltaPiLinear]
          ];
          require[case["name"] <> " ε=" <> fmt[eps] <> " normalization", Abs[norm - 1.0] <= normalizationTol,
            "integral=" <> fmt[norm]];
          require[case["name"] <> " ε=" <> fmt[eps] <> " positivity", minSigma >= -positivityTol,
            "min sigma=" <> fmt[minSigma]];
          require[case["name"] <> " ε=" <> fmt[eps] <> " cosine moment stays in the positive-source range",
            -positivityTol <= gEps <= 1.0 + positivityTol,
            "g=" <> fmt[gEps]];
          require[case["name"] <> " ε=" <> fmt[eps] <> " first-order deltaPi envelope",
            deltaPiRel <= sample["deltaPi_rel_envelope"],
            "relative error=" <> fmt[deltaPiRel] <> ", envelope=" <> fmt[sample["deltaPi_rel_envelope"]]];
        ],
        case["samples"]
      ];
    ],
    config["transport_cases"]
  ];
];

Catch[
  checkRepresentativeProfiles[];
  checkTransport[];
  Print[""];
  Print["Stage 154 co-evolving map stress harness passed."];
  Exit[0];
];

Exit[1];
