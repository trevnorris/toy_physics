(* Stress harness for the full-profile mouth correction chain (Stages 133-135). *)

ClearAll["Global`*"];
$HistoryLength = 0;

rootDir = DirectoryName[$InputFileName];
configPath = FileNameJoin[{
   rootDir, "..", "..", "..", "scripts", "moving_throat", "numerical",
   "stage133_135_profile_correction_samples.json"
}];
config = Import[ExpandFileName[configPath], "RawJSON"];

If[config["schema"] =!= "moving_throat_numerical_stage133_135_v1",
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

valueTol = N[tolerances["value_tol"]];
derivativeTol = N[tolerances["derivative_tol"]];
curvatureTol = N[tolerances["curvature_tol"]];

PiStarExpected = N[const["Pi_star_expected"]];
SigmaMStarExpected = N[const["Sigma_m_star_expected"]];
gStarExpected = N[const["g_star_expected"]];
SStarExpected = N[const["S_star_expected"]];
gPrimeStarExpected = N[const["gprime_star_expected"]];
ATExpected = N[const["A_T_expected"]];
BTExpected = N[const["B_T_expected"]];
THatStarExpected = N[const["T_hat_star_expected"]];
covCRExpected = N[const["cov_cR_expected"]];
covKRExpected = N[const["cov_KR_expected"]];
deltaGExpected = N[const["delta_g_expected"]];
deltaSExpected = N[const["delta_S_expected"]];
deltaPiExpected = N[const["deltaPi_expected"]];
deltaTExpected = N[const["deltaT_expected"]];
g1Expected = N[const["g1_expected"]];
S1Expected = N[const["S1_expected"]];
Pi1Expected = N[const["Pi1_expected"]];
T1Expected = N[const["T1_expected"]];

fdSteps = N /@ config["finite_difference_steps"];
curvatureSteps = N /@ config["curvature_steps"];
lambdaSamples = config["lambda_samples"];

kappa = Pi/2;
rF1 = Sqrt[4107 - 100 Pi^2]/(10 Pi);
gMinus = rF1 - Sqrt[1 + rF1^2]/2;

gFormula[p_] := 2 p (2 p Exp[p] + Pi)/((4 p^2 + Pi^2) (Exp[p] - 1));
sFormula[p_] := p (kappa Tanh[kappa] + p (Exp[-p] Sech[kappa] - 1))/((1 - Exp[-p]) (kappa^2 - p^2));

PiStar = p /. FindRoot[gFormula[p] == gMinus, {p, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
SStar = N[sFormula[PiStar], 50];
gStar = N[gFormula[PiStar], 50];
gPrimeStar = N[D[gFormula[p], p] /. p -> PiStar, 50];
sPrimeStar = N[D[sFormula[p], p] /. p -> PiStar, 50];
SigmaMStar = N[PiStar/(4 - SStar), 50];
THatStar = N[Sqrt[9 (PiStar/(1 - SStar/4))/20], 50];
AT = N[
  -(9/(40 THatStar)) (1/(gPrimeStar (1 - SStar/4)) + PiStar sPrimeStar/(4 gPrimeStar (1 - SStar/4)^2)),
  50
];
BT = N[(9/(40 THatStar)) PiStar/(4 (1 - SStar/4)^2), 50];

tS[x_] := (1 - Exp[-PiStar x])/(PiStar (1 - Exp[-PiStar])) - x Exp[-PiStar]/(1 - Exp[-PiStar]);
tQ[x_] := Module[{cQ, aQ},
  cQ = PiStar/((1 - Exp[-PiStar]) (kappa^2 - PiStar^2));
  aQ = cQ (kappa Sinh[kappa] + PiStar Exp[-PiStar])/(kappa Cosh[kappa]);
  aQ Sinh[kappa x] - cQ Cosh[kappa x] + cQ Exp[-PiStar x]
];
sigmaStar[x_] := PiStar Exp[-PiStar x]/(1 - Exp[-PiStar]);
rStar[x_] := SigmaMStar (4 tS[x] - tQ[x]) - PiStar x;
cKernel[x_] := Cos[Pi x/2];
sKernel[x_] := Cosh[(Pi/2) (1 - x)]/Cosh[Pi/2];

nint[expr_] := Quiet[
  NIntegrate[
    Evaluate[expr],
    {t, 0, 1},
    WorkingPrecision -> 60,
    AccuracyGoal -> 20,
    PrecisionGoal -> 20,
    Method -> {"GlobalAdaptive", "MaxErrorIncreases" -> 100}
  ],
  NIntegrate::precw
];

checkStage133[] := Module[{targetR2, forward, second},
  Print["=== Stage 133: full-profile residual geometry ==="];
  require["Pi_* matches expected canonical bias",
    nearQ[PiStar, PiStarExpected, valueTol],
    "Pi_*=" <> fmt[PiStar] <> ", expected=" <> fmt[PiStarExpected]];
  require["Sigma_m* matches expected canonical gain",
    nearQ[SigmaMStar, SigmaMStarExpected, valueTol],
    "Sigma_m*=" <> fmt[SigmaMStar] <> ", expected=" <> fmt[SigmaMStarExpected]];
  require["g_* matches expected canonical overlap",
    nearQ[gStar, gStarExpected, valueTol],
    "g_*=" <> fmt[gStar] <> ", expected=" <> fmt[gStarExpected]];
  require["S_* matches expected canonical response",
    nearQ[SStar, SStarExpected, valueTol],
    "S_*=" <> fmt[SStar] <> ", expected=" <> fmt[SStarExpected]];
  require["R_*(0)=0",
    Abs[rStar[0]] <= valueTol,
    "R(0)=" <> fmt[rStar[0]]];

  Scan[
    Function[{h},
      forward = (rStar[h] - rStar[0])/h;
      Print["forward derivative at h=", fmt[h], " -> ", fmt[forward]];
      require["R_*'(0) small at h=" <> fmt[h],
        Abs[forward] <= derivativeTol,
        "forward diff=" <> fmt[forward]];
      require["R_*(h) is negative at h=" <> fmt[h],
        rStar[h] < 0,
        "R(h)=" <> fmt[rStar[h]]];
    ],
    fdSteps
  ];

  targetR2 = -3 SigmaMStar PiStar/(1 - Exp[-PiStar]);
  Scan[
    Function[{h},
      second = (rStar[h] - 2 rStar[0] + rStar[-h])/h^2;
      Print["central second derivative at h=", fmt[h], " -> ", fmt[second]];
      require["R_*''(0) matches target at h=" <> fmt[h],
        Abs[second - targetR2] <= curvatureTol,
        "finite diff=" <> fmt[second] <> ", target=" <> fmt[targetR2]];
      require["R_*''(0) stays negative at h=" <> fmt[h],
        second < 0,
        "finite diff=" <> fmt[second]];
    ],
    curvatureSteps
  ];
];

checkStage135[] := Module[
  {avgR, avgC, avgS, covCR, covKR, deltaG, deltaS, deltaPi, deltaT, piErrors = {}, tErrors = {}, sigmaLam, zLam, gLam, sLam, deltaGActual, deltaSActual, deltaPiActual, deltaTActual, deltaGLinear, deltaSLinear, deltaPiLinear, deltaTLinear, piRel, tRel, lam},
  Print[""];
  Print["=== Stage 135: covariance correction and nonlinear stress ==="];

  avgR = nint[sigmaStar[t] rStar[t]];
  avgC = nint[sigmaStar[t] cKernel[t]];
  avgS = nint[sigmaStar[t] sKernel[t]];

  covCR = nint[sigmaStar[t] (cKernel[t] - avgC) (rStar[t] - avgR)];
  covKR = nint[sigmaStar[t] (sKernel[t] - avgS) (rStar[t] - avgR)];

  deltaG = -covCR;
  deltaS = -covKR;
  deltaPi = covCR/gPrimeStar;
  deltaT = AT deltaG + BT deltaS;

  require["Cov(c,R_*) matches Stage 135",
    nearQ[covCR, covCRExpected, valueTol],
    "cov_cR=" <> fmt[covCR] <> ", expected=" <> fmt[covCRExpected]];
  require["Cov(K_q,R_*) matches Stage 135",
    nearQ[covKR, covKRExpected, valueTol],
    "cov_KR=" <> fmt[covKR] <> ", expected=" <> fmt[covKRExpected]];
  require["delta g_act matches Stage 135",
    nearQ[deltaG, deltaGExpected, valueTol],
    "delta_g=" <> fmt[deltaG] <> ", expected=" <> fmt[deltaGExpected]];
  require["delta S_act matches Stage 135",
    nearQ[deltaS, deltaSExpected, valueTol],
    "delta_S=" <> fmt[deltaS] <> ", expected=" <> fmt[deltaSExpected]];
  require["delta Pi_act matches Stage 135",
    nearQ[deltaPi, deltaPiExpected, valueTol],
    "deltaPi=" <> fmt[deltaPi] <> ", expected=" <> fmt[deltaPiExpected]];
  require["delta T_act matches Stage 135",
    nearQ[deltaT, deltaTExpected, valueTol],
    "deltaT=" <> fmt[deltaT] <> ", expected=" <> fmt[deltaTExpected]];

  Scan[
    Function[{case},
      lam = N[case["lambda"]];
      sigmaLam[x_] := Exp[-PiStar x - lam rStar[x]];
      zLam = nint[sigmaLam[t]];
      gLam = nint[(sigmaLam[t]/zLam) cKernel[t]];
      sLam = nint[(sigmaLam[t]/zLam) sKernel[t]];
      deltaGActual = gLam - gStar;
      deltaSActual = sLam - SStar;
      deltaPiActual = -deltaGActual/gPrimeStar;
      deltaTActual = AT deltaGActual + BT deltaSActual;

      deltaGLinear = lam deltaG;
      deltaSLinear = lam deltaS;
      deltaPiLinear = lam deltaPi;
      deltaTLinear = lam deltaT;

      piRel = Abs[deltaPiActual - deltaPiLinear]/Abs[deltaPiLinear];
      tRel = Abs[deltaTActual - deltaTLinear]/Abs[deltaTLinear];
      AppendTo[piErrors, piRel];
      AppendTo[tErrors, tRel];

      Print[""];
      Print[case["name"], " (", case["kind"], "): lambda=", fmt[lam]];
      Print[
        "  delta_g_actual=", fmt[deltaGActual], ", delta_g_linear=", fmt[deltaGLinear],
        ", delta_S_actual=", fmt[deltaSActual], ", delta_S_linear=", fmt[deltaSLinear]
      ];
      Print[
        "  deltaPi_actual=", fmt[deltaPiActual], ", deltaPi_linear=", fmt[deltaPiLinear],
        ", deltaT_actual=", fmt[deltaTActual], ", deltaT_linear=", fmt[deltaTLinear]
      ];

      require[case["name"] <> " keeps the correction direction",
        deltaGActual < 0 && deltaSActual < 0 && deltaPiActual > 0 && deltaTActual > 0,
        "delta_g=" <> fmt[deltaGActual] <> ", delta_S=" <> fmt[deltaSActual] <> ", deltaPi=" <> fmt[deltaPiActual] <> ", deltaT=" <> fmt[deltaTActual]];
      require[case["name"] <> " deltaPi stays within the linear envelope",
        piRel <= N[case["relative_envelope"]],
        "relative error=" <> fmt[piRel] <> ", envelope=" <> fmt[N[case["relative_envelope"]]]];
      require[case["name"] <> " deltaT stays within the linear envelope",
        tRel <= N[case["relative_envelope"]],
        "relative error=" <> fmt[tRel] <> ", envelope=" <> fmt[N[case["relative_envelope"]]]];

      If[lam == 1,
        require["lambda=1 g_1 matches Stage 135",
          nearQ[gLam, g1Expected, valueTol],
          "g_1=" <> fmt[gLam] <> ", expected=" <> fmt[g1Expected]];
        require["lambda=1 S_1 matches Stage 135",
          nearQ[sLam, S1Expected, valueTol],
          "S_1=" <> fmt[sLam] <> ", expected=" <> fmt[S1Expected]];
        require["lambda=1 Pi_1 matches Stage 135",
          nearQ[PiStar + deltaPiActual, Pi1Expected, valueTol],
          "Pi_1=" <> fmt[PiStar + deltaPiActual] <> ", expected=" <> fmt[Pi1Expected]];
        require["lambda=1 T_1 matches Stage 135",
          nearQ[THatStar + deltaTActual, T1Expected, valueTol],
          "T_1=" <> fmt[THatStar + deltaTActual] <> ", expected=" <> fmt[T1Expected]];
      ];
    ],
    lambdaSamples
  ];

  require["deltaPi linearization error grows with lambda",
    AllTrue[Differences[piErrors], # > 0 &],
    "errors=" <> StringRiffle[fmt /@ piErrors, ", "]];
  require["deltaT linearization error grows with lambda",
    AllTrue[Differences[tErrors], # > 0 &],
    "errors=" <> StringRiffle[fmt /@ tErrors, ", "]];
];

If[
  Catch[
    Print["Loaded config from ", ExpandFileName[configPath]];
    checkStage133[];
    checkStage135[];
    Print[""];
    Print["All stage-133/135 profile-correction stress checks passed."];
    "ok"
  ] === "ok",
  Null,
  Print[""];
  Print["Stage-133/135 profile-correction stress harness failed."];
  Exit[1]
];
