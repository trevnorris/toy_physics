(* ===================================================================== *)
(* DYNAMIC SLAB SCALING TEST (DERIVATION-DRIVEN, PATCHED v3)              *)
(* ===================================================================== *)
(* What this script does
   1) Derives u(R) for 3D sink from continuity + spherical symmetry:
        u3D(R)=Q/(4 Pi R^2)
   2) Derives u(R) for slab/2D radial flow (uniform across thickness 2H):
        u2D(R)=Q/(4 Pi H R)
   3) Defines boundary by Bernoulli back-reaction:
        DeltaP(u(Rb)) = DeltaPcrit
   4) Sets H_eff(M) = Rb(M) and predicts plateau v_inf via:
        v_inf^2 = G M / H_eff(M)
   5) Automatically reports:
        - global fits
        - 3D-dominated subset fits (H_eff/Hgeom < rLow)
        - 2D-dominated subset fits (H_eff/Hgeom > rHigh)
        - transition mass where H_eff/Hgeom ~ 1 (if crossed)
        - local TF slope vs mass (reveals bends cleanly)
   6) Sweeps across all levers + exports plots and CSVs per case.
*)

ClearAll["Global`*"];

(* ----------------------------- *)
(* 0. Configuration              *)
(* ----------------------------- *)

baseConfig = <|
   "G" -> 1.0,

   "rho0" -> 1.0,
   "p0" -> 1.0,

   "DeltaPcrit" -> 10^-6,

   (* Mass grid (edit freely) *)
   "Mlist" -> Table[10.^x, {x, 9, 12, 0.5}],

   (* Q(M)=kappa*M^alphaQ *)
   "kappa" -> 2.0,
   "alphaQ" -> 1.0,

   (* EOS: "Incompressible" or "Polytrope" *)
   "EOS" -> "Incompressible",
   "gamma" -> 5/3,
   "K" -> 1.0,

   (* Flow model: "3D" | "2D" | "Crossover" *)
   "FlowModel" -> "Crossover",
   "Hgeom" -> 10^6,
   "CrossoverSharpness" -> 4,

   (* Regime thresholds by ratio H_eff/Hgeom *)
   "rLow" -> 0.3,
   "rHigh" -> 3.0,

   (* Numerics *)
   "WorkingPrecision" -> 50,
   "MaxIterations" -> 200,

   (* Output *)
   "Verbose" -> True,
   "MakePlots" -> True,
   "ExportResults" -> True,
   "OutDir" -> "slab_scaling_outputs",

   (* Sweeps (edit lists to broaden/target) *)
   "RunSweeps" -> True,
   "alphaQList" -> {0.5, 1.0, 1.5},
   "kappaList" -> {1.0, 2.0},
   "EOSList" -> {"Incompressible", "Polytrope"},
   "DeltaPcritList" -> {10^-7, 10^-6, 10^-5},
   "FlowModelList" -> {"3D", "Crossover"},
   "HgeomList" -> {10^6, 10^7}
|>;

(* ----------------------------------------------------------- *)
(* 1. DERIVATIONS: u(R) for 3D and slab (2D radial)             *)
(* ----------------------------------------------------------- *)

If[TrueQ[baseConfig["Verbose"]],
  Print["\n================================================="];
  Print["DERIVATION CHECK: u(R) from div u=0 + symmetry (3D and slab)"];
  Print["================================================="];

  Module[{u, sol},
    sol = DSolveValue[{D[R^2 u[R], R] == 0}, u, R];
    Print["3D: d/dR (R^2 u)=0 => u(R) = ", sol];
    Print["3D flux: Q = 4 Pi R^2 u(R) => u3D(R)=Q/(4 Pi R^2)."];
  ];

  Print["\nSlab/2D radial (uniform across thickness 2H):"];
  Print["Flux through cylindrical surface: Q = (2H)(2PiR) u(R) = 4 Pi H R u(R)"];
  Print["=> u2D(R)=Q/(4 Pi H R)."];
];

u3D[R_, Q_] := Q/(4 Pi R^2);
u2D[R_, Q_, H_] := Q/(4 Pi H R);

wCrossover[R_, H_, n_] := 1/(1 + (R/H)^n);

uModel[R_, Q_, cfg_Association] := Module[
  {fm = cfg["FlowModel"], H = cfg["Hgeom"], n = cfg["CrossoverSharpness"], w},
  Switch[fm,
    "3D", u3D[R, Q],
    "2D", u2D[R, Q, H],
    "Crossover", (w = wCrossover[R, H, n];
                  w*u3D[R, Q] + (1 - w)*u2D[R, Q, H]),
    _, u3D[R, Q]
  ]
];

QfromM[M_, cfg_Association] := cfg["kappa"]*M^cfg["alphaQ"];

(* ----------------------------------------------------------- *)
(* 2. Bernoulli + EOS: DeltaP(u)                                *)
(* ----------------------------------------------------------- *)

DeltaPIncompressible[u_, rho0_] := (1/2) rho0 u^2;

PolytropeRho[u_, rho0_, gamma_, K_] := Module[{term},
  term = rho0^(gamma - 1) - ((gamma - 1)/(2 gamma K)) u^2;
  If[term <= 0, 0, term^(1/(gamma - 1))]
];

DeltaPPolytrope[u_, rho0_, gamma_, K_] := Module[{rho, p0, p},
  rho = PolytropeRho[u, rho0, gamma, K];
  p0 = K rho0^gamma;
  p  = K rho^gamma;
  p0 - p
];

DeltaP[u_, cfg_Association] := Module[{eos = cfg["EOS"]},
  Which[
    eos === "Incompressible", DeltaPIncompressible[u, cfg["rho0"]],
    eos === "Polytrope",      DeltaPPolytrope[u, cfg["rho0"], cfg["gamma"], cfg["K"]],
    True, $Failed
  ]
];

(* ----------------------------------------------------------- *)
(* 3. Solve boundary: DeltaP(u(R))=DeltaPcrit => Rb(M)=H_eff(M) *)
(* ----------------------------------------------------------- *)

uStarIncompressible[cfg_] := Sqrt[2 cfg["DeltaPcrit"]/cfg["rho0"]];
RbGuess3D[Q_, cfg_] := Sqrt[Q/(4 Pi uStarIncompressible[cfg])];
RbGuess2D[Q_, cfg_] := Q/(4 Pi cfg["Hgeom"] uStarIncompressible[cfg]);

BracketRoot[f_, R0_] := Module[
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

RbNumeric[Q_, cfg_Association] := Module[
  {DeltaPcrit = cfg["DeltaPcrit"], wp = cfg["WorkingPrecision"],
   it = cfg["MaxIterations"], f, R0, br, sol},

  f[R_?NumericQ] := DeltaP[uModel[R, Q, cfg], cfg] - DeltaPcrit;

  R0 = Max[10^-30, RbGuess3D[Q, cfg]];
  br = BracketRoot[f, R0];

  If[br === $Failed,
    R0 = Max[10^-30, RbGuess2D[Q, cfg]];
    br = BracketRoot[f, R0];
  ];

  If[br === $Failed, Return[$Failed]];

  sol = Quiet@Check[
     FindRoot[f[R] == 0, {R, br[[1]], br[[2]]},
              WorkingPrecision -> wp, MaxIterations -> it],
     $Failed
  ];

  If[sol === $Failed, $Failed, R /. sol]
];

HeffFromM[M_, cfg_Association] := Module[{Q, Rb},
  Q = QfromM[M, cfg];
  Rb = RbNumeric[Q, cfg];
  Rb
];

vInfinity[M_, Heff_, cfg_Association] := Sqrt[cfg["G"] M/Heff];

(* ----------------------------------------------------------- *)
(* 4. Regimes, fits, transition mass, local slope               *)
(* ----------------------------------------------------------- *)

RegimeFromRatio[ratio_, cfg_Association] := Which[
  Not[NumericQ[ratio]], "no-root",
  ratio < cfg["rLow"], "3D-dominated",
  ratio > cfg["rHigh"], "2D-dominated",
  True, "transition"
];

FitLine[data_, x_Symbol] := Module[{fit, expr, slope, intercept, r2},
  If[Length[data] < 2,
    Return[<|"Slope" -> Missing["TooFewPoints"], "Intercept" -> Missing["TooFewPoints"], "RSquared" -> Missing["TooFewPoints"]|>]
  ];
  fit = LinearModelFit[data, {1, x}, x];
  expr = fit["BestFit"];
  slope = Coefficient[expr, x, 1];
  intercept = Coefficient[expr, x, 0];
  r2 = fit["RSquared"];
  <|"Fit" -> fit, "Slope" -> slope, "Intercept" -> intercept, "RSquared" -> r2|>
];

TransitionMass[goodRows_, cfg_Association] := Module[
  {M = goodRows[[All, 1]], ratio = goodRows[[All, 5]], lr, lm, i, x1, x2, y1, y2, xt},
  lr = Log10[ratio];
  lm = Log10[M];
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

LocalSlopeLogMLogv[goodRows_] := Module[
  {M = goodRows[[All, 1]], v = goodRows[[All, 4]], x, y, mids, slopes},
  x = Log10[v]; y = Log10[M];
  slopes = Table[
    If[x[[i + 1]] == x[[i]], Indeterminate, (y[[i + 1]] - y[[i]])/(x[[i + 1]] - x[[i]])],
    {i, 1, Length[x] - 1}
  ];
  mids = Sqrt[M[[;; -2]]*M[[2 ;;]]];
  Select[Transpose[{mids, slopes}], NumericQ[Last[#]] &]
];

(* ----------------------------------------------------------- *)
(* 5. Run one case                                               *)
(* ----------------------------------------------------------- *)

RunCase[cfg_Association] := Module[
  {Mlist = cfg["Mlist"], rows, failures = 0, Q, Heff, v, ratio, regime,
   goodRows, lowRows, highRows, transRows,
   fitHAll, fitTFAll, fitHLow, fitTFLow, fitHHigh, fitTFHigh,
   Mtrans, ratioStats, slopeLocal, pltRatio, pltTF, pltSlope, out},

  rows = Table[
    Q = QfromM[M, cfg];
    Heff = HeffFromM[M, cfg];
    If[Heff === $Failed || !NumericQ[Heff] || Heff <= 0,
      failures++;
      {M, Q, Missing["no-root"], Missing["no-root"], Missing["no-root"], "no-root"},
      v = vInfinity[M, Heff, cfg];
      ratio = Heff/cfg["Hgeom"];
      regime = RegimeFromRatio[ratio, cfg];
      {M, Q, Heff, v, ratio, regime}
    ],
    {M, Mlist}
  ];

  goodRows = Select[rows, NumericQ[#[[3]]] && NumericQ[#[[4]]] && NumericQ[#[[5]]] &];
  lowRows  = Select[goodRows, #[[5]] < cfg["rLow"] &];
  highRows = Select[goodRows, #[[5]] > cfg["rHigh"] &];
  transRows= Select[goodRows, cfg["rLow"] <= #[[5]] <= cfg["rHigh"] &];

  fitHAll  = FitLine[Log10 /@ goodRows[[All, {1, 3}]], x];
  fitTFAll = FitLine[Log10 /@ goodRows[[All, {4, 1}]], x];

  fitHLow  = FitLine[Log10 /@ lowRows[[All, {1, 3}]], x];
  fitTFLow = FitLine[Log10 /@ lowRows[[All, {4, 1}]], x];

  fitHHigh  = FitLine[Log10 /@ highRows[[All, {1, 3}]], x];
  fitTFHigh = FitLine[Log10 /@ highRows[[All, {4, 1}]], x];

  Mtrans = TransitionMass[goodRows, cfg];

  ratioStats = <|
    "min" -> If[Length[goodRows] > 0, Min[goodRows[[All, 5]]], Missing["no-data"]],
    "max" -> If[Length[goodRows] > 0, Max[goodRows[[All, 5]]], Missing["no-data"]],
    "median" -> If[Length[goodRows] > 0, Median[goodRows[[All, 5]]], Missing["no-data"]]
  |>;

  slopeLocal = If[Length[goodRows] >= 2, LocalSlopeLogMLogv[goodRows], {}];

  pltRatio = ListLogLogPlot[
    goodRows[[All, {1, 5}]],
    Frame -> True,
    FrameLabel -> {"M", "H_eff/H_geom"},
    PlotLabel -> Row[{
      "Ratio vs mass (FM=", cfg["FlowModel"], ", EOS=", cfg["EOS"], ", alphaQ=", cfg["alphaQ"],
      ", kappa=", cfg["kappa"], ", DPcrit=", cfg["DeltaPcrit"], ", Hgeom=", cfg["Hgeom"], ")"
    }],
    GridLines -> {None, {{1}}},
    PlotRange -> All,
    ImageSize -> Large
  ];

  pltTF = ListLogLogPlot[
    goodRows[[All, {1, 4}]],
    Frame -> True,
    FrameLabel -> {"M", "v_inf"},
    PlotLabel -> "Plateau speed vs mass (bends/flattening show crossover)",
    GridLines -> Automatic,
    PlotRange -> All,
    ImageSize -> Large
  ];

  pltSlope = ListLogLogPlot[
    slopeLocal,
    Frame -> True,
    FrameLabel -> {"M (midpoints)", "local slope d log M / d log v"},
    PlotLabel -> "Local TF slope vs mass (BTFR ~4 in 3D-dominated regime)",
    GridLines -> {None, {{4}}},
    PlotRange -> All,
    ImageSize -> Large
  ];

  out = <|
    "Config" -> cfg,
    "Rows" -> rows,
    "GoodRows" -> goodRows,
    "Failures" -> failures,
    "Counts" -> <|
      "nTotal" -> Length[rows],
      "nGood" -> Length[goodRows],
      "nLow" -> Length[lowRows],
      "nTrans" -> Length[transRows],
      "nHigh" -> Length[highRows]
    |>,
    "RatioStats" -> ratioStats,
    "TransitionMass" -> Mtrans,
    "Fits" -> <|
      "H_all" -> fitHAll, "TF_all" -> fitTFAll,
      "H_low" -> fitHLow, "TF_low" -> fitTFLow,
      "H_high" -> fitHHigh, "TF_high" -> fitTFHigh
    |>,
    "Plots" -> <|"Ratio" -> pltRatio, "TF" -> pltTF, "LocalSlope" -> pltSlope|>
  |>;

  If[TrueQ[cfg["Verbose"]],
    Print["\n================================================="];
    Print["CASE SUMMARY"];
    Print["================================================="];
    Print["EOS: ", cfg["EOS"], " | FM: ", cfg["FlowModel"], " | alphaQ=", cfg["alphaQ"], " | kappa=", cfg["kappa"],
          " | DPcrit=", cfg["DeltaPcrit"], " | Hgeom=", cfg["Hgeom"]];
    Print["Failures: ", failures, " of ", Length[rows]];
    Print["Counts: ", out["Counts"]];
    Print["Ratio stats: ", ratioStats];
    Print["Transition mass (ratio~1): ", out["TransitionMass"]];
    Print["Global bH=", out["Fits"]["H_all"]["Slope"], "  R2=", out["Fits"]["H_all"]["RSquared"]];
    Print["Global bTF=", out["Fits"]["TF_all"]["Slope"], "  R2=", out["Fits"]["TF_all"]["RSquared"]];
    Print["Low  (ratio<", cfg["rLow"], ") bTF_low=", out["Fits"]["TF_low"]["Slope"], "  R2=", out["Fits"]["TF_low"]["RSquared"]];
    Print["High (ratio>", cfg["rHigh"], ") bTF_high=", out["Fits"]["TF_high"]["Slope"], "  R2=", out["Fits"]["TF_high"]["RSquared"]];
  ];

  out
];

(* ----------------------------------------------------------- *)
(* 6. Sweeps + exports                                           *)
(* ----------------------------------------------------------- *)

EnsureDir[path_] := If[!DirectoryQ[path], CreateDirectory[path, CreateIntermediateDirectories -> True]];

NumKey[x_] := StringReplace[ToString[FortranForm@N[x]], {"." -> "p", "-" -> "m", "+" -> ""}];
SafeKey[s_] := StringReplace[ToString[s, InputForm], {" " -> "", "/" -> "_", "." -> "p", "-" -> "m"}];

CaseID[cfg_] := StringJoin[
  "aQ", SafeKey[cfg["alphaQ"]],
  "_k",  SafeKey[cfg["kappa"]],
  "_EOS", SafeKey[cfg["EOS"]],
  "_DP", NumKey[cfg["DeltaPcrit"]],
  "_FM", SafeKey[cfg["FlowModel"]],
  "_Hg", NumKey[cfg["Hgeom"]]
];

SummaryRow[out_] := Module[{cfg = out["Config"], fits = out["Fits"], c = out["Counts"], rs = out["RatioStats"]},
  <|
    "alphaQ" -> cfg["alphaQ"],
    "kappa" -> cfg["kappa"],
    "EOS" -> cfg["EOS"],
    "DeltaPcrit" -> cfg["DeltaPcrit"],
    "FlowModel" -> cfg["FlowModel"],
    "Hgeom" -> cfg["Hgeom"],
    "nGood" -> c["nGood"],
    "nLow" -> c["nLow"],
    "nTrans" -> c["nTrans"],
    "nHigh" -> c["nHigh"],
    "ratioMin" -> rs["min"],
    "ratioMed" -> rs["median"],
    "ratioMax" -> rs["max"],
    "Mtrans" -> out["TransitionMass"],
    "bH_all" -> fits["H_all"]["Slope"],
    "R2_H_all" -> fits["H_all"]["RSquared"],
    "bTF_all" -> fits["TF_all"]["Slope"],
    "R2_TF_all" -> fits["TF_all"]["RSquared"],
    "bTF_low" -> fits["TF_low"]["Slope"],
    "R2_TF_low" -> fits["TF_low"]["RSquared"],
    "bTF_high" -> fits["TF_high"]["Slope"],
    "R2_TF_high" -> fits["TF_high"]["RSquared"]
  |>
];

ExportCase[out_] := Module[{cfg = out["Config"], dir = cfg["OutDir"], id, grid},
  EnsureDir[dir];
  id = CaseID[cfg];
  grid = GraphicsGrid[{
      {out["Plots"]["Ratio"], out["Plots"]["TF"]},
      {out["Plots"]["LocalSlope"], Graphics[{}]}
    }, Spacings -> {0.5, 0.5}];
  Export[FileNameJoin[{dir, id <> ".pdf"}], grid];
  Export[FileNameJoin[{dir, id <> "_rows.csv"}],
    Prepend[out["GoodRows"], {"M", "Q", "H_eff", "v_inf", "Heff_over_Hgeom", "Regime"}]
  ];
];

RunSweeps[base_Association] := Module[
  {dir = base["OutDir"], combos, outs, summaries, cfg},

  If[TrueQ[base["ExportResults"]], EnsureDir[dir]];

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

  outs = Table[
    cfg = combos[[i]];
    If[TrueQ[base["Verbose"]],
      Print["\n--- Running case ", i, " of ", Length[combos], ": ", CaseID[cfg], " ---"];
    ];
    With[{o = RunCase[cfg]},
      If[TrueQ[base["ExportResults"]], ExportCase[o]];
      o
    ],
    {i, 1, Length[combos]}
  ];

  summaries = SummaryRow /@ outs;

  If[TrueQ[base["ExportResults"]],
    Export[FileNameJoin[{dir, "sweep_summary.csv"}], Dataset[summaries] // Normal];
  ];

  Dataset[summaries]
];

(* ----------------------------------------------------------- *)
(* 7. Main execution                                             *)
(* ----------------------------------------------------------- *)

baselineOut = RunCase[baseConfig];

Print["\n--- Preview table (M, Q(M), H_eff(M), v_inf(M), H_eff/Hgeom, Regime) ---"];
Print[Grid[
  Prepend[
    (If[NumericQ[#], NumberForm[#, {10, 3}], #] & /@ #) & /@ baselineOut["Rows"],
    {"M", "Q", "H_eff", "v_inf", "Heff/Hgeom", "Regime"}
  ],
  Frame -> All
]];

If[TrueQ[baseConfig["MakePlots"]],
  Show[baselineOut["Plots"]["Ratio"]];
  Show[baselineOut["Plots"]["TF"]];
  Show[baselineOut["Plots"]["LocalSlope"]];
];

If[TrueQ[baseConfig["RunSweeps"]],
  Print["\n================================================="];
  Print["RUNNING FULL SWEEP (all levers)"];
  Print["================================================="];
  sweepSummary = RunSweeps[baseConfig];
  Print["\nSweep summary Dataset created (also exported to CSV if enabled)."];
  sweepSummary
];

