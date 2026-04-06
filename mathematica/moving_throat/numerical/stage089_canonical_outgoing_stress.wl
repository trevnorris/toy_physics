(* Branch-sensitivity stress harness for Stage 089. *)

ClearAll["Global`*"];
$HistoryLength = 0;

rootDir = DirectoryName[$InputFileName];
configPath = FileNameJoin[{
   rootDir, "..", "..", "..", "scripts", "moving_throat", "numerical",
   "stage089_canonical_outgoing_samples.json"
}];
config = Import[ExpandFileName[configPath], "RawJSON"];

If[config["schema"] =!= "moving_throat_numerical_stage089_v1",
  Print["Unexpected config schema."];
  Exit[1];
];

fmt[x_] := ToString @ NumberForm[N[x, 14], {Infinity, 12}, ExponentFunction -> (Null &)];

require[label_, condition_, detail_] := Module[{status},
  status = If[TrueQ[condition], "PASS", "FAIL"];
  Print["[", status, "] ", label, ": ", detail];
  If[!TrueQ[condition], Throw[$Failed]]
];

sourceMapModel[aOverR_] := 1 + aOverR^2;
nqValue[mhat0_, chiQ_] := 1/(mhat0^2 chiQ);

stage089Block[] := Module[{canonicalErrors = {}, aOverR, mhat0, chiQ, nq, gammaEffRatio, coeffRelError, nameSmall, aSmall, errSmall, nameLarge, aLarge, errLarge, separation},
  Print["=== Stage 089: canonical outgoing closure stress ==="];
  Print["source-map model: ", config["source_map_model"]];

  Scan[
    Function[{case},
      aOverR = N[case["a_over_r"]];
      mhat0 = sourceMapModel[aOverR];
      chiQ = 1.0;
      nq = nqValue[mhat0, chiQ];
      gammaEffRatio = mhat0^2 nq;
      coeffRelError = Abs[nq - 1];
      AppendTo[canonicalErrors, {case["name"], aOverR, coeffRelError}];
      Print[""];
      Print[case["name"], " (", case["kind"], "): a/r=", fmt[aOverR], ", mhat0=", fmt[mhat0]];
      Print["  N_Q=", fmt[nq], ", coefficient ratio=", fmt[nq], ", gamma_eff ratio=", fmt[gammaEffRatio]];
      require[case["name"] <> " canonical gamma_eff hits target",
        Abs[gammaEffRatio - 1] <= 10^-12,
        "gamma_eff/target=" <> fmt[gammaEffRatio]];
      require[case["name"] <> " canonical branch keeps N_Q <= 1",
        nq <= 1,
        "N_Q=" <> fmt[nq]];
    ],
    config["canonical_samples"]
  ];

  Scan[
    Function[{pair},
      {nameSmall, aSmall, errSmall} = pair[[1]];
      {nameLarge, aLarge, errLarge} = pair[[2]];
      require[nameSmall <> " is closer to target than " <> nameLarge,
        errSmall < errLarge,
        "|N_Q-1|(" <> fmt[aSmall] <> ")=" <> fmt[errSmall] <> ", |N_Q-1|(" <> fmt[aLarge] <> ")=" <> fmt[errLarge]];
    ],
    Partition[canonicalErrors, 2, 1]
  ];

  Scan[
    Function[{case},
      aOverR = N[case["a_over_r"]];
      mhat0 = sourceMapModel[aOverR];
      chiQ = N[case["chi_Q"]];
      nq = nqValue[mhat0, chiQ];
      gammaEffRatio = mhat0^2 nq;
      separation = Abs[gammaEffRatio - 1];
      Print[""];
      Print[case["name"], " (", case["kind"], "): a/r=", fmt[aOverR], ", mhat0=", fmt[mhat0], ", chi_Q=", fmt[chiQ]];
      Print["  N_Q=", fmt[nq], ", coefficient ratio=", fmt[nq], ", gamma_eff ratio=", fmt[gammaEffRatio]];
      require[case["name"] <> " branch deformation moves gamma_eff off target",
        separation >= N[case["separation_floor"]],
        "|gamma_eff/target - 1|=" <> fmt[separation] <> ", floor=" <> fmt[N[case["separation_floor"]]]];
      require[case["name"] <> " gamma_eff follows 1/chi_Q",
        Abs[gammaEffRatio - 1/chiQ] <= 10^-12,
        "gamma_eff/target=" <> fmt[gammaEffRatio] <> ", 1/chi_Q=" <> fmt[1/chiQ]];
    ],
    config["branch_perturbations"]
  ];

  Print[""];
  Print["All stage-089 canonical outgoing stress checks passed."];
];

If[
  Catch[
    Print["Loaded config from ", ExpandFileName[configPath]];
    stage089Block[];
    "ok"
  ] === "ok",
  Null,
  Print[""];
  Print["Stage-089 canonical outgoing stress harness failed."];
  Exit[1]
];
