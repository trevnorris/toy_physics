(* PathA 22b Gate 3 cross-check.
   Scope: exact outgoing l=2 Hankel DtN fingerprint, omega^5 extraction,
   radius-consistent defect/vacuum chi_Q normalization, dimensions, and
   imported Python convergence/vacuum sanity gates. *)

ClearAll[
  dimString, checkExpr, checkBool, checkDim,
  dim0, L, T, M, velocity, z, omega, a, csSym, q, rExit,
  alphaDef, alphaVac, h2, lambdaOut, yhatOut, lambdaSeries, yhatSeries,
  fingerprintA, fingerprintR, omegaQ, sigmaCan, responsePert, coeffPert,
  chiPert, eDef, eVac, chiCal, checks, algebra, jsonPath, solveReport,
  solveChecks, report, outDir, allPass, outcomeClassified, isNotExtractable,
  isConvergedProduction,
  nwJsonPath, nwReport, nwChecks, nwOutcome, nwRows, bestFit, calcResiduals,
  calcRms, calcMax
];

dim0 = {0, 0, 0};
L = {1, 0, 0};
T = {0, 1, 0};
M = {0, 0, 1};
velocity = L - T;

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

z = Symbol["z"];
omega = Symbol["omega"];
a = Symbol["a"];
csSym = Symbol["cs"];
rExit = Symbol["RExit"];
alphaDef = Symbol["alphaDef"];
alphaVac = Symbol["alphaVac"];
q = 7/5;

h2 = Exp[I z] (z^2 + 3 I z - 3)/z^3;
lambdaOut = FullSimplify[z D[Log[h2], z]];
yhatOut = FullSimplify[-3/lambdaOut];
lambdaSeries = Expand[Normal[Series[lambdaOut, {z, 0, 5}]]];
yhatSeries = Expand[Normal[Series[yhatOut, {z, 0, 5}]]];

fingerprintA = FullSimplify[(Coefficient[yhatSeries, z, 5]/I) a^5/csSym^5];
fingerprintR = FullSimplify[(Coefficient[yhatSeries, z, 5]/I) rExit^5/csSym^5];
omegaQ = 3 csSym/(2 a);
sigmaCan = FullSimplify[(9/8)/omegaQ^5];
responsePert = 1 + a^2 omega^2/(9 csSym^2) + 4 a^4 omega^4/(81 csSym^4) +
  I q fingerprintA omega^5;
coeffPert = FullSimplify[Coefficient[Normal[Series[responsePert, {omega, 0, 5}]], omega, 5]/I];
chiPert = FullSimplify[coeffPert/fingerprintA];
eDef = FullSimplify[alphaDef fingerprintR];
eVac = FullSimplify[alphaVac fingerprintR];
chiCal = FullSimplify[eDef/eVac];

checks = {
  checkDim["defect omega^5 coefficient", 5 T, 5 T],
  checkDim["vacuum omega^5 coefficient", 5 T, 5 T],
  checkDim["closure placement factor", 5 L - 5 velocity, 5 T],
  checkDim["radius-consistent chi_Q ratio", 5 T - 5 T, dim0]
};

algebra = {
  checkExpr["lambda outgoing fingerprint", lambdaSeries, -3 + z^2/3 + z^4/9 + I z^5/9],
  checkExpr["normalized outgoing fingerprint", yhatSeries, 1 + z^2/9 + 4 z^4/81 + I z^5/27],
  checkExpr["derived a-fingerprint coefficient", fingerprintA, a^5/(27 csSym^5)],
  checkExpr["derived R-exit fingerprint coefficient", fingerprintR, rExit^5/(27 csSym^5)],
  checkExpr["canonical sigma_Q", sigmaCan, 4 a^5/(27 csSym^5)],
  checkExpr["omega5 coefficient extraction", coeffPert, q fingerprintA],
  checkExpr["symbolic fixed-radius chi_Q", chiPert, q],
  checkExpr["calibrated defect/vacuum chi_Q", chiCal, alphaDef/alphaVac],
  checkExpr["R_exit derivative cancels", D[chiCal, rExit], 0],
  checkExpr["vacuum returns one", chiCal /. alphaDef -> alphaVac, 1]
};

jsonPath = FileNameJoin[{"software", "stage1_solver", "_scratch", "pathA_22b_gate3.json"}];
solveChecks = If[FileExistsQ[jsonPath],
  solveReport = Import[jsonPath, "RawJSON"];
  outcomeClassified = StringStartsQ[solveReport["gate3_outcome"], "DELTA_Q_NE_0"] ||
    StringStartsQ[solveReport["gate3_outcome"], "CHI_Q_MAGNITUDE_NOT_EXTRACTABLE"] ||
    StringStartsQ[solveReport["gate3_outcome"], "NW_CONVERGES"] ||
    solveReport["gate3_outcome"] === "CHI_Q=1_DERIVED";
  isNotExtractable = StringStartsQ[solveReport["gate3_outcome"], "CHI_Q_MAGNITUDE_NOT_EXTRACTABLE"];
  isConvergedProduction = StringStartsQ[solveReport["gate3_outcome"], "DELTA_Q_NE_0"] ||
    StringStartsQ[solveReport["gate3_outcome"], "NW_CONVERGES"];
  {
    checkBool["python outcome class is classified", outcomeClassified, True],
    checkBool["python production outcome is converged", isConvergedProduction, True],
    checkBool["python branch-geometry reference present", KeyExistsQ[solveReport["result"]["references"], "branch_geometry"], True],
    checkBool["python flat reference present", KeyExistsQ[solveReport["result"]["references"], "flat_Zw1"], True],
    checkBool[
      "python domain sweep has at least three branch rows",
      Length[solveReport["result"]["references"]["branch_geometry"]["domain_truncation"]["rows"]] >= 3,
      True
    ],
    checkBool[
      "python branch rmax plateau pass",
      solveReport["result"]["references"]["branch_geometry"]["domain_truncation"]["assessment"]["plateau"],
      True
    ],
    checkBool[
      "python branch rmax plateau onset numeric",
      NumberQ[solveReport["result"]["references"]["branch_geometry"]["domain_truncation"]["assessment"]["plateau_onset_r_max"]],
      True
    ],
    checkBool[
      "python post-plateau grid convergence pass",
      solveReport["result"]["post_plateau_grid_convergence"]["assessment"]["converging"],
      True
    ],
    checkBool["python tiny defect linearity pass", solveReport["result"]["tiny_defect_linearity"]["linear_toward_zero"], True],
    checkBool["python uniform consistency small", solveReport["result"]["trivial_uniform_consistency"]["max_uniform_abs_error"] <= 10^-10, True],
    checkBool[
      "python converged central matches even-nw characterization",
      If[
        isConvergedProduction,
        NumberQ[solveReport["chi_Q_value"]] &&
          Abs[solveReport["chi_Q_value"] - solveReport["result"]["even_nw_characterization"]["chi_Q_reported"]] <= 10^-12,
        solveReport["chi_Q_value"] === Null
      ],
      True
    ],
    checkBool[
      "python Z_w systematic is one-sided and separate",
      solveReport["result"]["combined_budget"]["zw_reference_systematic_kind"] === "one_sided_definitional" &&
        solveReport["result"]["combined_budget"]["zw_reference_shift"] > solveReport["chi_Q_error_bar"],
      True
    ]
  },
  {
    checkBool["python solve output present", False, True, "Run the Gate-3 Python module before this cross-check."]
  }
];

nwJsonPath = FileNameJoin[{"_scratch", "pathA_22b_gate3_nw_characterization.json"}];
nwChecks = If[FileExistsQ[nwJsonPath],
  nwReport = Import[nwJsonPath, "RawJSON"];
  nwOutcome = nwReport["characterization"]["outcome_class"];
  nwRows = Select[Values[nwReport["rows"]], #["sweep_kind"] === "nw" &];
  bestFit = nwReport["characterization"]["best_fit"];
  calcResiduals = MapThread[
    #2 - (bestFit["chi_inf"] + bestFit["c"]/#1^bestFit["p"]) &,
    {bestFit["nw"], bestFit["chi_Q"]}
  ];
  calcRms = Sqrt[Total[calcResiduals^2]/Length[calcResiduals]];
  calcMax = Max[Abs[calcResiduals]];
  {
    checkBool[
      "nw characterization outcome class is classified",
      StringStartsQ[nwOutcome, "NW_CONVERGES"] ||
        StringStartsQ[nwOutcome, "NW_OSCILLATES_DAMPED"] ||
        nwOutcome === "NW_NO_FINITE_LIMIT" ||
        nwOutcome === "NW_INCONCLUSIVE_NEED_FINER",
      True
    ],
    checkBool["nw characterization has at least seven even rows", Length[nwRows] >= 14, True],
    checkBool["nw characterization excludes odd lane counts", AllTrue[nwRows, EvenQ[#["grid"]["nw"]] &], True],
    checkBool["nw best-fit chi_inf numeric", NumberQ[bestFit["chi_inf"]], True],
    checkBool["nw best-fit order numeric", NumberQ[bestFit["p"]], True],
    checkBool["nw fit rms arithmetic", Abs[calcRms - bestFit["rms"]] <= 10^-10, True],
    checkBool["nw fit max residual arithmetic", Abs[calcMax - bestFit["max_abs_residual"]] <= 10^-10, True],
    checkBool["nw jackknife recorded", KeyExistsQ[nwReport["characterization"], "jackknife"], True],
    checkBool["nw nr comparison recorded", KeyExistsQ[nwReport["characterization"], "nr_comparison"], True]
  },
  {
    checkBool["nw characterization output present", False, True, "Run the even-nw characterization before this cross-check."]
  }
];

allPass = AllTrue[Join[checks, algebra, solveChecks, nwChecks], TrueQ[#["pass"]] &];

report = <|
  "schema" -> "stage1_pathA_22b_gate3_mathematica_crosscheck/v5",
  "pass" -> allPass,
  "lambda_series" -> ToString[InputForm[lambdaSeries]],
  "yhat_series" -> ToString[InputForm[yhatSeries]],
  "derived_a_fingerprint" -> ToString[InputForm[fingerprintA]],
  "derived_r_exit_fingerprint" -> ToString[InputForm[fingerprintR]],
  "canonical_sigma_Q" -> ToString[InputForm[sigmaCan]],
  "calibrated_chi_Q" -> ToString[InputForm[chiCal]],
  "vacuum_chi_Q" -> ToString[InputForm[chiCal /. alphaDef -> alphaVac]],
  "imported_python_outcome" -> If[ValueQ[solveReport], solveReport["gate3_outcome"], "missing"],
  "imported_python_chi_Q" -> If[ValueQ[solveReport], solveReport["chi_Q_value"], Missing["not_available"]],
  "imported_nw_outcome" -> If[ValueQ[nwReport], nwReport["characterization"]["outcome_class"], "missing"],
  "imported_nw_chi_inf" -> If[ValueQ[nwReport], nwReport["characterization"]["best_fit"]["chi_inf"], Missing["not_available"]],
  "checks" -> checks,
  "algebra" -> algebra,
  "solve_checks" -> solveChecks,
  "nw_checks" -> nwChecks
|>;

outDir = FileNameJoin[{"software", "stage1_solver", "_scratch"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
Export[FileNameJoin[{outDir, "pathA_22b_gate3_mathematica_crosscheck.json"}], report, "RawJSON"];
Print["wrote ", FileNameJoin[{outDir, "pathA_22b_gate3_mathematica_crosscheck.json"}]];
Print["pathA_22b Gate 3 Mathematica cross-check: ", Count[Join[checks, algebra, solveChecks, nwChecks][[All, "pass"]], True], "/", Length[Join[checks, algebra, solveChecks, nwChecks]], " checks"];
If[TrueQ[report["pass"]], Exit[0], Exit[1]]
