(* ===================================================================== *)
(* dynamicSlabScalingTest.wl                                              *)
(* Version v4.2 — robust sweep/exports + regime split fits + bend plots    *)
(* ===================================================================== *)

Print["\n[dynamicSlabScalingTest] version v4.2"];

ClearAll["Global`*"];

(* ----------------------------- *)
(* 0. Configuration              *)
(* ----------------------------- *)

baseConfig = <|
  "G" -> 1.0,

  "rho0" -> 1.0,
  "DeltaPcrit" -> 10^-6,

  (* Mass grid *)
  "Mlist" -> Table[10.^x, {x, 9, 12, 0.5}],

  (* Q(M)=kappa*M^alphaQ *)
  "kappa" -> 2.0,
  "alphaQ" -> 1.0,

  (* EOS *)
  "EOS" -> "Incompressible",          (* "Incompressible" or "Polytrope" *)
  "gamma" -> 5/3,
  "K" -> 1.0,

  (* Flow model *)
  "FlowModel" -> "Crossover",         (* "3D" | "2D" | "Crossover" *)
  "Hgeom" -> 10^6,
  "CrossoverSharpness" -> 4,

  (* Regime thresholds by ratio heff/hgeom *)
  "rLow" -> 0.3,
  "rHigh" -> 3.0,

  (* Numerics *)
  "WorkingPrecision" -> 50,
  "MaxIterations" -> 200,

  (* Output *)
  "OutDir" -> "slabScalingOutputs",
  "MakePlots" -> True,
  "ExportResults" -> True,

  (* Logging control *)
  "PrintBaselineSummary" -> True,
  "PrintCaseSummaries" -> False,
  "PrintProgress" -> True,
  "ProgressEvery" -> 25,

  (* Sweeps *)
  "RunSweeps" -> True,
  "alphaQList" -> {0.5, 1.0, 1.5},
  "kappaList" -> {1.0, 2.0},
  "EOSList" -> {"Incompressible", "Polytrope"},
  "DeltaPcritList" -> {10^-7, 10^-6, 10^-5},
  "FlowModelList" -> {"3D", "Crossover"},
  "HgeomList" -> {10^6, 10^7}
|>;

(* ----------------------------- *)
(* Helpers: safe config lookups  *)
(* ----------------------------- *)

lookupKey[assoc_Association, k_String, kSym_Symbol, default_] := Module[{v},
  v = Lookup[assoc, k, Missing["no"]];
  If[v =!= Missing["no"], Return[v]];
  v = Lookup[assoc, kSym, Missing["no"]];
  If[v =!= Missing["no"], Return[v]];
  default
];

outDirFrom[cfg_Association] := Module[{d},
  d = lookupKey[cfg, "OutDir", OutDir, "slabScalingOutputs"];
  If[!StringQ[d] || StringLength[d] == 0, d = "slabScalingOutputs"];
  d
];

ensureDir[dir_] := Module[{d = dir},
  If[!StringQ[d] || StringLength[d] == 0, d = "slabScalingOutputs"];
  If[!DirectoryQ[d], CreateDirectory[d, CreateIntermediateDirectories -> True]];
  d
];

shortNum[x_] := If[NumericQ[x], NumberForm[N[x], {6, 3}], x];

(* ----------------------------- *)
(* 1) Derived flow fields         *)
(* ----------------------------- *)

u3D[R_, Q_] := Q/(4 Pi R^2);
u2D[R_, Q_, H_] := Q/(4 Pi H R);
wCrossover[R_, H_, n_] := 1/(1 + (R/H)^n);

uModel[R_, Q_, cfg_Association] := Module[
  {fm = cfg["FlowModel"], H = cfg["Hgeom"], n = cfg["CrossoverSharpness"], w},
  Switch[fm,
    "3D", u3D[R, Q],
    "2D", u2D[R, Q, H],
    "Crossover",
      (w = wCrossover[R, H, n];
       w*u3D[R, Q] + (1 - w)*u2D[R, Q, H]),
    _, u3D[R, Q]
  ]
];

qFromM[M_, cfg_Association] := cfg["kappa"]*M^cfg["alphaQ"];

(* ----------------------------- *)
(* 2) Bernoulli / EOS closure     *)
(* ----------------------------- *)

deltaPIncompressible[u_, rho0_] := (1/2) rho0 u^2;

polytropeRho[u_, rho0_, gamma_, K_] := Module[{term},
  term = rho0^(gamma - 1) - ((gamma - 1)/(2 gamma K)) u^2;
  If[term <= 0, 0, term^(1/(gamma - 1))]
];

deltaPPolytrope[u_, rho0_, gamma_, K_] := Module[{rho, p0, p},
  rho = polytropeRho[u, rho0, gamma, K];
  p0 = K rho0^gamma;
  p = K rho^gamma;
  p0 - p
];

deltaP[u_, cfg_Association] := Module[{eos = cfg["EOS"]},
  Which[
    eos === "Incompressible", deltaPIncompressible[u, cfg["rho0"]],
    eos === "Polytrope", deltaPPolytrope[u, cfg["rho0"], cfg["gamma"], cfg["K"]],
    True, $Failed
  ]
];

uStarIncompressible[cfg_Association] := Sqrt[2 cfg["DeltaPcrit"]/cfg["rho0"]];
rbGuess3D[Q_, cfg_Association] := Sqrt[Q/(4 Pi uStarIncompressible[cfg])];
rbGuess2D[Q_, cfg_Association] := Q/(4 Pi cfg["Hgeom"] uStarIncompressible[cfg]);

bracketRoot[f_, R0_] := Module[
  {a = R0/10., b = R0*10., fa, fb, k = 0, grow = 3., maxK = 20},
  fa = Quiet@Check[f[a], Indeterminate];
  fb = Quiet@Check[f[b], Indeterminate];
  While[
    (k < maxK) && (Not[NumericQ[fa]] || Not[NumericQ[fb]] || fa*fb > 0),
    a = a/grow; b = b*grow;
    fa = Quiet@Check[f[a], Indeterminate];
    fb = Quiet@Check[f[b], Indeterminate];
    k++
  ];
  If[NumericQ[fa] && NumericQ[fb] && fa*fb <= 0, {a, b}, $Failed]
];

rbNumeric[Q_, cfg_Association] := Module[
  {dpcrit = cfg["DeltaPcrit"], wp = cfg["WorkingPrecision"], it = cfg["MaxIterations"], f, R0, br, sol},
  f[R_?NumericQ] := deltaP[uModel[R, Q, cfg], cfg] - dpcrit;

  R0 = Max[10^-30, rbGuess3D[Q, cfg]];
  br = bracketRoot[f, R0];
  If[br === $Failed,
    R0 = Max[10^-30, rbGuess2D[Q, cfg]];
    br = bracketRoot[f, R0];
  ];
  If[br === $Failed, Return[$Failed]];

  sol = Quiet@Check[
    FindRoot[f[R] == 0, {R, br[[1]], br[[2]]}, WorkingPrecision -> wp, MaxIterations -> it],
    $Failed
  ];
  If[sol === $Failed, $Failed, R /. sol]
];

heffFromM[M_, cfg_Association] := Module[{Q, rb},
  Q = qFromM[M, cfg];
  rb = rbNumeric[Q, cfg];
  rb
];

vinf[M_, heff_, cfg_Association] := Sqrt[cfg["G"] M/heff];

(* ----------------------------- *)
(* 3) Fits + regime diagnostics   *)
(* ----------------------------- *)

fitLine[data_, x_Symbol] := Module[{fit, expr, slope, intercept, r2},
  If[Length[data] < 2,
    Return[<|"Slope" -> Missing["TooFewPoints"], "Intercept" -> Missing["TooFewPoints"], "RSq" -> Missing["TooFewPoints"]|>]
  ];
  fit = LinearModelFit[data, {1, x}, x];
  expr = fit["BestFit"];
  slope = Coefficient[expr, x, 1];
  intercept = Coefficient[expr, x, 0];
  r2 = fit["RSquared"];
  <|"Fit" -> fit, "Slope" -> slope, "Intercept" -> intercept, "RSq" -> r2|>
];

transitionMass[goodRows_] := Module[
  {M = goodRows[[All, 1]], ratio = goodRows[[All, 5]], lr, lm, i, x1, x2, y1, y2, xt},
  lr = Log10[ratio]; lm = Log10[M];
  i = FirstPosition[
    Partition[lr, 2, 1],
    {a_, b_} /; (NumericQ[a] && NumericQ[b] && a <= 0 && b >= 0),
    Missing["NoCross"]
  ];
  If[i === Missing["NoCross"], Return[Missing["NoCross"]]];
  i = i[[1]];
  x1 = lm[[i]]; x2 = lm[[i + 1]];
  y1 = lr[[i]]; y2 = lr[[i + 1]];
  xt = x1 + (0 - y1)*(x2 - x1)/(y2 - y1);
  10^xt
];

localSlopeLogMLogV[goodRows_] := Module[
  {M = goodRows[[All, 1]], v = goodRows[[All, 4]], x, y, mids, slopes},
  x = Log10[v]; y = Log10[M];
  slopes = Table[
    If[x[[i + 1]] == x[[i]], Indeterminate, (y[[i + 1]] - y[[i]])/(x[[i + 1]] - x[[i]])],
    {i, 1, Length[x] - 1}
  ];
  mids = Sqrt[M[[;; -2]]*M[[2 ;;]]];
  Select[Transpose[{mids, slopes}], NumericQ[Last[#]] &]
];

(* ----------------------------- *)
(* 4) Run one case                *)
(* ----------------------------- *)

runCase[cfg_Association] := Module[
  {Mlist = cfg["Mlist"], rows, failures = 0, Q, heff, v, ratio,
   goodRows, lowRows, highRows,
   fitHAll, fitTFAll, fitTFLow, fitTFHigh,
   ratioStats, mTrans, slopeLocal, pltRatio, pltTF, pltSlope},

  rows = Table[
    Q = qFromM[M, cfg];
    heff = heffFromM[M, cfg];
    If[heff === $Failed || !NumericQ[heff] || heff <= 0,
      failures++;
      {M, Q, Missing["noRoot"], Missing["noRoot"], Missing["noRoot"]},
      v = vinf[M, heff, cfg];
      ratio = heff/cfg["Hgeom"];
      {M, Q, heff, v, ratio}
    ],
    {M, Mlist}
  ];

  goodRows = Select[rows, NumericQ[#[[3]]] && NumericQ[#[[4]]] && NumericQ[#[[5]]] &];
  lowRows  = Select[goodRows, #[[5]] < cfg["rLow"] &];
  highRows = Select[goodRows, #[[5]] > cfg["rHigh"] &];

  fitHAll  = fitLine[Log10 /@ goodRows[[All, {1, 3}]], x];
  fitTFAll = fitLine[Log10 /@ goodRows[[All, {4, 1}]], x];

  fitTFLow  = fitLine[Log10 /@ lowRows[[All, {4, 1}]], x];
  fitTFHigh = fitLine[Log10 /@ highRows[[All, {4, 1}]], x];

  ratioStats = <|
    "min" -> If[Length[goodRows] > 0, Min[goodRows[[All, 5]]], Missing["noData"]],
    "max" -> If[Length[goodRows] > 0, Max[goodRows[[All, 5]]], Missing["noData"]],
    "med" -> If[Length[goodRows] > 0, Median[goodRows[[All, 5]]], Missing["noData"]]
  |>;

  mTrans = If[Length[goodRows] > 1, transitionMass[goodRows], Missing["NoData"]];
  slopeLocal = If[Length[goodRows] >= 2, localSlopeLogMLogV[goodRows], {}];

  pltRatio = ListLogLogPlot[
    goodRows[[All, {1, 5}]],
    Frame -> True,
    FrameLabel -> {"M", "heff/hgeom"},
    GridLines -> {None, {{1}}},
    PlotRange -> All,
    ImageSize -> Large,
    PlotLabel -> "Ratio (heff/hgeom) vs mass"
  ];

  pltTF = ListLogLogPlot[
    goodRows[[All, {1, 4}]],
    Frame -> True,
    FrameLabel -> {"M", "vinf"},
    GridLines -> Automatic,
    PlotRange -> All,
    ImageSize -> Large,
    PlotLabel -> "vinf vs mass (shows flattening/bends)"
  ];

  pltSlope = ListLogLogPlot[
    slopeLocal,
    Frame -> True,
    FrameLabel -> {"M midpoints", "local slope dlogM/dlogv"},
    GridLines -> {None, {{4}}},
    PlotRange -> All,
    ImageSize -> Large,
    PlotLabel -> "Local TF slope vs mass"
  ];

  <|
    "Config" -> cfg,
    "Rows" -> rows,
    "GoodRows" -> goodRows,
    "Failures" -> failures,
    "Counts" -> <|"nTotal" -> Length[rows], "nGood" -> Length[goodRows], "nLow" -> Length[lowRows], "nHigh" -> Length[highRows]|>,
    "RatioStats" -> ratioStats,
    "TransitionMass" -> mTrans,
    "Fits" -> <|"bH" -> fitHAll, "bTF" -> fitTFAll, "bTFlow" -> fitTFLow, "bTFhigh" -> fitTFHigh|>,
    "Plots" -> <|"ratio" -> pltRatio, "tf" -> pltTF, "slope" -> pltSlope|>
  |>
];

(* ----------------------------- *)
(* 5) Sweep + export              *)
(* ----------------------------- *)

numKey[x_] := StringReplace[ToString[FortranForm@N[x]], {"." -> "p", "-" -> "m", "+" -> ""}];

safeKey[x_] := Module[{t},
  t = If[StringQ[x], x, ToString[x, InputForm]];
  StringReplace[t, {"\"" -> "", " " -> "", "/" -> "-", "." -> "p", "-" -> "m", "+" -> ""}]
];

caseId[cfg_Association] := StringRiffle[
  {
    "aQ" <> safeKey[cfg["alphaQ"]],
    "k" <> safeKey[cfg["kappa"]],
    "EOS" <> safeKey[cfg["EOS"]],
    "DP" <> numKey[cfg["DeltaPcrit"]],
    "FM" <> safeKey[cfg["FlowModel"]],
    "Hg" <> numKey[cfg["Hgeom"]]
  },
  "-"
];

summaryAssoc[out_] := Module[{cfg = out["Config"], fits = out["Fits"], c = out["Counts"], rs = out["RatioStats"]},
  <|
    "alphaQ" -> cfg["alphaQ"],
    "kappa" -> cfg["kappa"],
    "EOS" -> cfg["EOS"],
    "DeltaPcrit" -> cfg["DeltaPcrit"],
    "FlowModel" -> cfg["FlowModel"],
    "Hgeom" -> cfg["Hgeom"],
    "nGood" -> c["nGood"],
    "nLow" -> c["nLow"],
    "nHigh" -> c["nHigh"],
    "ratioMin" -> rs["min"],
    "ratioMed" -> rs["med"],
    "ratioMax" -> rs["max"],
    "Mtrans" -> out["TransitionMass"],
    "bH" -> fits["bH"]["Slope"],
    "R2H" -> fits["bH"]["RSq"],
    "bTF" -> fits["bTF"]["Slope"],
    "R2TF" -> fits["bTF"]["RSq"],
    "bTFlow" -> fits["bTFlow"]["Slope"],
    "R2TFlow" -> fits["bTFlow"]["RSq"],
    "bTFhigh" -> fits["bTFhigh"]["Slope"],
    "R2TFhigh" -> fits["bTFhigh"]["RSq"]
  |>
];

exportCase[out_] := Module[{cfg = out["Config"], dir, id, grid, rowsFile, pdfFile},
  dir = ensureDir[outDirFrom[cfg]];
  id = caseId[cfg];

  grid = GraphicsGrid[
    {
      {out["Plots"]["ratio"], out["Plots"]["tf"]},
      {out["Plots"]["slope"], Graphics[{}]}
    },
    Spacings -> {0.5, 0.5}
  ];

  pdfFile = FileNameJoin[{dir, id <> ".pdf"}];
  rowsFile = FileNameJoin[{dir, id <> "-rows.csv"}];

  Export[pdfFile, grid];
  Export[rowsFile, Prepend[out["GoodRows"], {"M", "Q", "heff", "vinf", "heffOverHgeom"}]];
];

runSweeps[base_Association] := Module[
  {dir = ensureDir[outDirFrom[base]], combos, outs, summaries, cfg, i, n, pe = base["ProgressEvery"]},

  combos = Flatten@Table[
    <|
      base,
      "alphaQ" -> aQ,
      "kappa" -> kap,
      "EOS" -> eos,
      "DeltaPcrit" -> dp,
      "FlowModel" -> fm,
      "Hgeom" -> hg
    |>,
    {aQ, base["alphaQList"]},
    {kap, base["kappaList"]},
    {eos, base["EOSList"]},
    {dp, base["DeltaPcritList"]},
    {fm, base["FlowModelList"]},
    {hg, base["HgeomList"]}
  , 5];

  n = Length[combos];

  outs = Table[
    cfg = combos[[i]];

    If[TrueQ[base["PrintProgress"]] && (i == 1 || i == n || Mod[i, pe] == 0),
      Print["Running case ", i, " of ", n, ": ", caseId[cfg]];
    ];

    With[{o = runCase[cfg]},
      If[TrueQ[base["PrintCaseSummaries"]],
        Print["  Summary: bTF=", shortNum[o["Fits"]["bTF"]["Slope"]],
              "  R2TF=", shortNum[o["Fits"]["bTF"]["RSq"]],
              "  ratioMed=", shortNum[o["RatioStats"]["med"]],
              "  Mtrans=", shortNum[o["TransitionMass"]]
        ];
      ];
      If[TrueQ[base["ExportResults"]], exportCase[o]];
      o
    ],
    {i, 1, n}
  ];

  summaries = summaryAssoc /@ outs;

  If[TrueQ[base["ExportResults"]],
    Export[FileNameJoin[{dir, "sweepSummary.csv"}], Dataset[summaries] // Normal];
  ];

  Dataset[summaries]
];

(* ----------------------------- *)
(* 6) Main execution              *)
(* ----------------------------- *)

baselineOut = runCase[baseConfig];

If[TrueQ[baseConfig["PrintBaselineSummary"]],
  Print["\n================================================="];
  Print["BASELINE SUMMARY"];
  Print["================================================="];
  Print["EOS=", baseConfig["EOS"], "  FM=", baseConfig["FlowModel"], "  alphaQ=", baseConfig["alphaQ"], "  kappa=", baseConfig["kappa"],
        "  DP=", baseConfig["DeltaPcrit"], "  Hgeom=", baseConfig["Hgeom"]];
  Print["Failures: ", baselineOut["Failures"], " of ", baselineOut["Counts"]["nTotal"]];
  Print["RatioStats: ", baselineOut["RatioStats"], "  Mtrans: ", baselineOut["TransitionMass"]];
  Print["bH=", baselineOut["Fits"]["bH"]["Slope"], "  R2H=", baselineOut["Fits"]["bH"]["RSq"]];
  Print["bTF=", baselineOut["Fits"]["bTF"]["Slope"], "  R2TF=", baselineOut["Fits"]["bTF"]["RSq"]];
  Print["bTFlow=", baselineOut["Fits"]["bTFlow"]["Slope"], "  bTFhigh=", baselineOut["Fits"]["bTFhigh"]["Slope"]];
];

Print["\n--- Preview table (M, Q(M), heff(M), vinf(M), heff/hgeom) ---"];
Print[Grid[
  Prepend[
    (If[NumericQ[#], NumberForm[#, {10, 3}], #] & /@ #) & /@ baselineOut["Rows"],
    {"M", "Q", "heff", "vinf", "heffOverHgeom"}
  ],
  Frame -> All
]];

If[TrueQ[baseConfig["MakePlots"]],
  Show[baselineOut["Plots"]["ratio"]];
  Show[baselineOut["Plots"]["tf"]];
  Show[baselineOut["Plots"]["slope"]];
];

If[TrueQ[baseConfig["RunSweeps"]],
  Print["\n================================================="];
  Print["RUNNING FULL SWEEP (all levers)"];
  Print["================================================="];
  sweepSummary = runSweeps[baseConfig];
  Print["\nSweep summary created (and exported if enabled)."];
  sweepSummary
];
