(* Finite-deformation stress harness for the Stage 146-148 rigidity bridge. *)

ClearAll["Global`*"];
$HistoryLength = 0;

rootDir = DirectoryName[$InputFileName];
configPath = FileNameJoin[{
   rootDir, "..", "..", "scripts", "numerical",
   "stage146_148_rigidity_samples.json"
}];
config = Import[ExpandFileName[configPath], "RawJSON"];

If[config["schema"] =!= "moving_throat_numerical_stage146_148_v1",
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
integralTol = N[tolerances["integral_tol"]];
compTol = N[tolerances["compensation_residual_tol"]];

PiStarExpected = N[const["Pi_star_expected"]];
gStarExpected = N[const["g_star_expected"]];
SStarExpected = N[const["S_star_expected"]];
gPrimeStarExpected = N[const["gprime_star_expected"]];
SPrimeStarExpected = N[const["Sprime_star_expected"]];
THatStarExpected = N[const["T_hat_star_expected"]];
ATExpected = N[const["A_T_expected"]];
BTExpected = N[const["B_T_expected"]];
lambdaPiZeroExpected = N[const["lambda_pi_zero_expected"]];
lambdaTZeroExpected = N[const["lambda_T_zero_expected"]];

kappa = Pi/2;
rF1 = Sqrt[4107 - 100 Pi^2]/(10 Pi);
gMinus = rF1 - Sqrt[1 + rF1^2]/2;
gFormula[p_] := 2 p (2 p Exp[p] + Pi)/((4 p^2 + Pi^2) (Exp[p] - 1));
SFormula[p_] := p (kappa Tanh[kappa] + p (Exp[-p] Sech[kappa] - 1))/((1 - Exp[-p]) (kappa^2 - p^2));

PiStar = p /. FindRoot[gFormula[p] == gMinus, {p, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
gStar = N[gFormula[PiStar], 50];
SStar = N[SFormula[PiStar], 50];
gPrimeStar = N[D[gFormula[p], p] /. p -> PiStar, 50];
SPrimeStar = N[D[SFormula[p], p] /. p -> PiStar, 50];
THatStar = N[Sqrt[9 (PiStar/(1 - SStar/4))/20], 50];
AT = N[
  -(9/(40 THatStar)) (1/(gPrimeStar (1 - SStar/4)) + PiStar SPrimeStar/(4 gPrimeStar (1 - SStar/4)^2)),
  50
];
BT = N[(9/(40 THatStar)) PiStar/(4 (1 - SStar/4)^2), 50];

gUniform = N[2/Pi, 50];
SUniform = N[2 Tanh[Pi/2]/Pi, 50];
gDerivative = N[Pi/4, 50];
SDerivative = N[(1 + Sinh[Pi/2])/(2 Cosh[Pi/2]), 50];

sigmaExp[x_, p_] := p Exp[-p x]/(1 - Exp[-p]);
varsigmaLambda[x_, lam_] := (1 - lam) + lam (Pi/2) Cos[Pi x/2];
gLambdaFormula[lam_] := (1 - lam) gUniform + lam gDerivative;
SLambdaFormula[lam_] := (1 - lam) SUniform + lam SDerivative;
sigmaTotal[x_, p_, eps_, lam_] := (1 - eps) sigmaExp[x, p] + eps varsigmaLambda[x, lam];

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

checkBaseConstants[] := Module[{},
  Print["=== Stage 146-148: base constants ==="];
  require["Pi_* matches audited value",
    nearQ[PiStar, PiStarExpected, valueTol],
    "Pi_*=" <> fmt[PiStar] <> ", expected=" <> fmt[PiStarExpected]];
  require["g_* matches audited value",
    nearQ[gStar, gStarExpected, valueTol],
    "g_*=" <> fmt[gStar] <> ", expected=" <> fmt[gStarExpected]];
  require["S_* matches audited value",
    nearQ[SStar, SStarExpected, valueTol],
    "S_*=" <> fmt[SStar] <> ", expected=" <> fmt[SStarExpected]];
  require["g_*' matches audited value",
    nearQ[gPrimeStar, gPrimeStarExpected, valueTol],
    "g_*'=" <> fmt[gPrimeStar] <> ", expected=" <> fmt[gPrimeStarExpected]];
  require["S_*' matches audited value",
    nearQ[SPrimeStar, SPrimeStarExpected, valueTol],
    "S_*'=" <> fmt[SPrimeStar] <> ", expected=" <> fmt[SPrimeStarExpected]];
  require["T_* matches audited value",
    nearQ[THatStar, THatStarExpected, valueTol],
    "T_*=" <> fmt[THatStar] <> ", expected=" <> fmt[THatStarExpected]];
  require["A_T matches audited value",
    nearQ[AT, ATExpected, valueTol],
    "A_T=" <> fmt[AT] <> ", expected=" <> fmt[ATExpected]];
  require["B_T matches audited value",
    nearQ[BT, BTExpected, valueTol],
    "B_T=" <> fmt[BT] <> ", expected=" <> fmt[BTExpected]];
];

runCase[case_] := Module[
  {name, lam, normVar, gVarInt, SVarInt, gVar, SVar, dPiPerEps, dSPerEps, dTPerEps, piErrors = {}, sErrors = {}, tErrors = {}, quadCoeffs = {}, eps, predictedDPi, predictedDS, predictedDT, pExact, sigmaNorm, gInt, SInt, compensationResidual, tExact, deltaPiExact, deltaSExact, deltaTExact, piRel, sRel, tRel, quadCoeff},
  name = case["name"];
  lam = N[Rationalize[case["lambda"], 0], 50];
  Print[""];
  Print["=== ", name, " (", case["kind"], ") ==="];
  Print["lambda = ", fmt[lam]];
  If[KeyExistsQ[case, "assumptions"],
    Print["assumptions:"];
    Scan[Function[{item}, Print["  - ", item]], case["assumptions"]];
  ];

  If[name === "bias_neutral_direction",
    require["bias-neutral lambda matches Stage 148",
      nearQ[lam, lambdaPiZeroExpected, valueTol],
      "lambda=" <> fmt[lam] <> ", expected=" <> fmt[lambdaPiZeroExpected]];
  ];
  If[name === "traction_neutral_direction",
    require["traction-neutral lambda matches Stage 148",
      nearQ[lam, lambdaTZeroExpected, valueTol],
      "lambda=" <> fmt[lam] <> ", expected=" <> fmt[lambdaTZeroExpected]];
  ];

  normVar = nint[varsigmaLambda[t, lam]];
  gVarInt = nint[varsigmaLambda[t, lam] Cos[Pi t/2]];
  SVarInt = nint[varsigmaLambda[t, lam] Cosh[(Pi/2) (1 - t)]/Cosh[Pi/2]];
  gVar = gLambdaFormula[lam];
  SVar = SLambdaFormula[lam];

  require["varsigma_lambda normalization",
    nearQ[normVar, 1, integralTol],
    "integral=" <> fmt[normVar]];
  require["varsigma_lambda is nonnegative at x=1",
    varsigmaLambda[1, lam] >= -10^-30,
    "varsigma(1)=" <> fmt[varsigmaLambda[1, lam]]];
  require["g_lambda integral matches formula",
    nearQ[gVarInt, gVar, integralTol],
    "integral=" <> fmt[gVarInt] <> ", formula=" <> fmt[gVar]];
  require["S_lambda integral matches formula",
    nearQ[SVarInt, SVar, integralTol],
    "integral=" <> fmt[SVarInt] <> ", formula=" <> fmt[SVar]];

  dPiPerEps = -(gVar - gStar)/gPrimeStar;
  dSPerEps = (SVar - SStar) - (SPrimeStar/gPrimeStar) (gVar - gStar);
  dTPerEps = AT (gVar - gStar) + BT (SVar - SStar);
  Print["linear dPi/eps=", fmt[dPiPerEps], ", linear dS/eps=", fmt[dSPerEps], ", linear dT/eps=", fmt[dTPerEps]];

  Scan[
    Function[{sample},
      eps = N[Rationalize[sample["epsilon"], 0], 50];
      predictedDPi = eps dPiPerEps;
      predictedDS = eps dSPerEps;
      predictedDT = eps dTPerEps;

      pExact = p /. Quiet[
        FindRoot[
          (1 - eps) gFormula[p] + eps gVar == gStar,
          {p, PiStar + predictedDPi},
          WorkingPrecision -> 50,
          AccuracyGoal -> 20,
          PrecisionGoal -> 20,
          MaxIterations -> 100
        ],
        FindRoot::precw
      ];
      sigmaNorm = nint[sigmaTotal[t, pExact, eps, lam]];
      gInt = nint[sigmaTotal[t, pExact, eps, lam] Cos[Pi t/2]];
      SInt = nint[sigmaTotal[t, pExact, eps, lam] Cosh[(Pi/2) (1 - t)]/Cosh[Pi/2]];
      compensationResidual = N[(1 - eps) gFormula[pExact] + eps gVar - gStar, 50];
      tExact = Sqrt[9 pExact/(20 (1 - SInt/4))];

      deltaPiExact = pExact - PiStar;
      deltaSExact = SInt - SStar;
      deltaTExact = tExact - THatStar;

      Print[""];
      Print["epsilon = ", fmt[eps]];
      Print["  Pi_exact=", fmt[pExact], ", deltaPi_exact=", fmt[deltaPiExact], ", deltaPi_linear=", fmt[predictedDPi]];
      Print["  deltaS_exact=", fmt[deltaSExact], ", deltaS_linear=", fmt[predictedDS], ", deltaT_exact=", fmt[deltaTExact], ", deltaT_linear=", fmt[predictedDT]];

      require["sigma_total normalization",
        nearQ[sigmaNorm, 1, integralTol],
        "integral=" <> fmt[sigmaNorm]];
      require["sigma_total stays nonnegative at x=1",
        sigmaTotal[1, pExact, eps, lam] >= -10^-30,
        "sigma_total(1)=" <> fmt[sigmaTotal[1, pExact, eps, lam]]];
      require["g_total integral matches g_*",
        nearQ[gInt, gStar, integralTol],
        "g_int=" <> fmt[gInt] <> ", g_*=" <> fmt[gStar]];
      require["S_total integral matches affine formula",
        nearQ[SInt, (1 - eps) SFormula[pExact] + eps SVar, integralTol],
        "S_int=" <> fmt[SInt]];
      require["exact compensation residual",
        TrueQ[Abs[compensationResidual] <= compTol],
        "residual=" <> fmt[compensationResidual]];

      If[KeyExistsQ[sample, "deltaPi_abs_tol"],
        require["epsilon=" <> fmt[eps] <> " bias shift stays neutral",
          Abs[deltaPiExact] <= N[sample["deltaPi_abs_tol"]],
          "deltaPi_exact=" <> fmt[deltaPiExact]],
        piRel = Abs[deltaPiExact - predictedDPi]/Abs[predictedDPi];
        AppendTo[piErrors, piRel];
        require["epsilon=" <> fmt[eps] <> " deltaPi stays within first-order envelope",
          piRel <= N[sample["deltaPi_rel_envelope"]],
          "relative error=" <> fmt[piRel] <> ", envelope=" <> fmt[N[sample["deltaPi_rel_envelope"]]]]
      ];

      sRel = Abs[deltaSExact - predictedDS]/Abs[predictedDS];
      AppendTo[sErrors, sRel];
      require["epsilon=" <> fmt[eps] <> " deltaS stays within first-order envelope",
        sRel <= N[sample["deltaS_rel_envelope"]],
        "relative error=" <> fmt[sRel] <> ", envelope=" <> fmt[N[sample["deltaS_rel_envelope"]]]];

      If[KeyExistsQ[sample, "deltaT_quadratic_coeff_max"],
        quadCoeff = Abs[deltaTExact]/eps^2;
        AppendTo[quadCoeffs, quadCoeff];
        require["epsilon=" <> fmt[eps] <> " deltaT remains quadratic-small",
          quadCoeff <= N[sample["deltaT_quadratic_coeff_max"]],
          "|deltaT|/eps^2=" <> fmt[quadCoeff] <> ", max=" <> fmt[N[sample["deltaT_quadratic_coeff_max"]]]],
        tRel = Abs[deltaTExact - predictedDT]/Abs[predictedDT];
        AppendTo[tErrors, tRel];
        require["epsilon=" <> fmt[eps] <> " deltaT stays within first-order envelope",
          tRel <= N[sample["deltaT_rel_envelope"]],
          "relative error=" <> fmt[tRel] <> ", envelope=" <> fmt[N[sample["deltaT_rel_envelope"]]]]
      ];
    ],
    case["samples"]
  ];

  If[Length[piErrors] > 1,
    require["deltaPi relative error grows with epsilon",
      AllTrue[Differences[piErrors], # > 0 &],
      "errors=" <> StringRiffle[fmt /@ piErrors, ", "]];
  ];
  If[Length[sErrors] > 1,
    require["deltaS relative error grows with epsilon",
      AllTrue[Differences[sErrors], # > 0 &],
      "errors=" <> StringRiffle[fmt /@ sErrors, ", "]];
  ];
  If[Length[tErrors] > 1,
    require["deltaT relative error grows with epsilon",
      AllTrue[Differences[tErrors], # > 0 &],
      "errors=" <> StringRiffle[fmt /@ tErrors, ", "]];
  ];
  If[Length[quadCoeffs] > 1,
    require["quadratic traction coefficient grows mildly with epsilon",
      AllTrue[Differences[quadCoeffs], # >= 0 &],
      "coeffs=" <> StringRiffle[fmt /@ quadCoeffs, ", "]];
  ];
];

If[
  Catch[
    Print["Loaded config from ", ExpandFileName[configPath]];
    checkBaseConstants[];
    Scan[runCase, config["cases"]];
    Print[""];
    Print["All stage-146/148 rigidity stress checks passed."];
    "ok"
  ] === "ok",
  Null,
  Print[""];
  Print["Stage 146/148 rigidity stress harness failed."];
  Exit[1]
];
