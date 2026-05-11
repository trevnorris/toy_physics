(* Branch-sensitivity stress harness for Stage 106. *)

ClearAll["Global`*"];
$HistoryLength = 0;

rootDir = DirectoryName[$InputFileName];
configPath = FileNameJoin[{
   rootDir, "..", "..", "scripts", "numerical",
   "stage106_canonical_outgoing_samples.json"
}];
config = Import[ExpandFileName[configPath], "RawJSON"];

If[config["schema"] =!= "moving_throat_numerical_stage106_v1",
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
canonicalPointParticleApprox[aOverR_] := 1 - 2 aOverR^2;
outgoingProbeZ[aOverR_] := Min[2.0*10^-2, Max[5.0*10^-4, 10.0*aOverR]];
outgoingYhatExact[z_] := Module[{lambdaOut},
  lambdaOut = z D[SphericalHankelH1[2, zz], zz]/SphericalHankelH1[2, zz] /. zz -> z;
  -3/lambdaOut
];
canonicalChiFromOutgoingDtN[z_] := Module[{yhat, evenPart},
  yhat = outgoingYhatExact[z];
  evenPart = 1 + z^2/9 + 4 z^4/81;
  27 Im[yhat - evenPart]/z^5
];

stage106Block[] := Module[
  {
    canonicalErrors = {}, perturbationPairs = <||>, aOverR, mhat0, chiQ, zProbe,
    chiExact, nq, ppApprox, remainder, scaledRemainder, coeffRelError,
    nameSmall, aSmall, errSmall, nameLarge, aLarge, errLarge,
    nqCanonical, branchKey, chiPlus, chiExactPlus, nqPlus, chiMinus, chiExactMinus, nqMinus,
    slope, curvature
  },
  Print["=== Stage 106: canonical outgoing closure stress ==="];
  Print["source-map model: ", config["source_map_model"]];

  Scan[
    Function[{case},
      aOverR = N[case["a_over_r"]];
      mhat0 = sourceMapModel[aOverR];
      chiQ = 1.0;
      zProbe = outgoingProbeZ[aOverR];
      chiExact = canonicalChiFromOutgoingDtN[zProbe];
      nq = nqValue[mhat0, chiQ];
      ppApprox = canonicalPointParticleApprox[aOverR];
      remainder = nq - ppApprox;
      scaledRemainder = remainder/(aOverR^4);
      coeffRelError = Abs[nq - 1];
      AppendTo[canonicalErrors, {case["name"], aOverR, coeffRelError}];

      Print[""];
      Print[
        case["name"], " (", case["kind"], "): a/r=", fmt[aOverR],
        ", mhat0=", fmt[mhat0], ", z_probe=", fmt[zProbe]
      ];
      Print[
        "  N_Q=", fmt[nq], ", point-particle approx=", fmt[ppApprox],
        ", chi_Q^out(exact)=", fmt[chiExact]
      ];

      require[case["name"] <> " canonical branch keeps N_Q <= 1",
        nq <= 1.0,
        "N_Q=" <> fmt[nq]];
      require[case["name"] <> " exact outgoing DtN fixes chi_Q = 1",
        Abs[chiExact - 1.0] <= 5.0*10^-5,
        "chi_Q^out(exact)=" <> fmt[chiExact]];
      require[case["name"] <> " outgoing DtN and source-map canonical N_Q agree",
        Abs[nqValue[mhat0, chiExact] - nq] <= 5.0*10^-5,
        "N_Q(chi_exact)=" <> fmt[nqValue[mhat0, chiExact]] <> ", N_Q(chi=1)=" <> fmt[nq]];
      require[case["name"] <> " exact N_Q stays above the O((a/r)^2) point-particle truncation",
        remainder > 0.0,
        "N_Q - [1 - 2(a/r)^2]=" <> fmt[remainder]];
      If[aOverR >= 10^-3,
        require[case["name"] <> " point-particle remainder is O((a/r)^4)",
          2.9 <= scaledRemainder <= 3.05,
          "(N_Q - [1 - 2(a/r)^2])/(a/r)^4=" <> fmt[scaledRemainder]],
        require[case["name"] <> " point-particle remainder stays below floating-point noise",
          remainder <= 10^-14,
          "N_Q - [1 - 2(a/r)^2]=" <> fmt[remainder]]
      ];
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
      zProbe = outgoingProbeZ[aOverR];
      chiExact = canonicalChiFromOutgoingDtN[zProbe];
      nq = nqValue[mhat0, chiQ];
      nqCanonical = nqValue[mhat0, chiExact];

      Print[""];
      Print[
        case["name"], " (", case["kind"], "): a/r=", fmt[aOverR],
        ", mhat0=", fmt[mhat0], ", chi_Q=", fmt[chiQ],
        ", chi_Q^out(exact)=", fmt[chiExact]
      ];
      Print[
        "  N_Q=", fmt[nq], ", canonical N_Q=", fmt[nqCanonical],
        ", chi gap=", fmt[chiQ - chiExact]
      ];

      If[chiQ > chiExact,
        branchKey = "positive";
        require[case["name"] <> " positive chi_Q shift lowers N_Q",
          nq < nqCanonical,
          "N_Q=" <> fmt[nq] <> ", canonical=" <> fmt[nqCanonical]];
        require[case["name"] <> " positive chi_Q shift stays above the outgoing DtN branch",
          chiQ > chiExact,
          "chi_Q=" <> fmt[chiQ] <> ", chi_Q^out(exact)=" <> fmt[chiExact]],
        branchKey = "negative";
        require[case["name"] <> " negative chi_Q shift raises N_Q",
          nq > nqCanonical,
          "N_Q=" <> fmt[nq] <> ", canonical=" <> fmt[nqCanonical]];
        require[case["name"] <> " negative chi_Q shift stays below the outgoing DtN branch",
          chiQ < chiExact,
          "chi_Q=" <> fmt[chiQ] <> ", chi_Q^out(exact)=" <> fmt[chiExact]]
      ];

      If[!KeyExistsQ[perturbationPairs, aOverR], perturbationPairs[aOverR] = <||>];
      perturbationPairs[aOverR][branchKey] = {chiQ, chiExact, nq};
    ],
    config["branch_perturbations"]
  ];

  Scan[
    Function[{aKey},
      If[Sort[Keys[perturbationPairs[aKey]]] =!= {"negative", "positive"},
        Throw[$Failed]
      ];
      {chiPlus, chiExactPlus, nqPlus} = perturbationPairs[aKey]["positive"];
      {chiMinus, chiExactMinus, nqMinus} = perturbationPairs[aKey]["negative"];
      chiExact = (chiExactPlus + chiExactMinus)/2;
      nqCanonical = nqValue[sourceMapModel[aKey], chiExact];
      slope = (nqPlus - nqMinus)/(chiPlus - chiMinus);
      curvature = nqPlus + nqMinus - 2.0*nqCanonical;

      require["a/r=" <> fmt[aKey] <> " symmetric branch slope around the outgoing DtN branch is negative",
        slope < 0.0,
        "slope=" <> fmt[slope]];
      require["a/r=" <> fmt[aKey] <> " symmetric branch response is centered on chi_Q^out",
        Abs[chiExactPlus - chiExactMinus] <= 1.0*10^-4,
        "chi_Q^out(+)= " <> fmt[chiExactPlus] <> ", chi_Q^out(-)= " <> fmt[chiExactMinus]];
      require["a/r=" <> fmt[aKey] <> " reciprocal conservative response is convex",
        curvature > 0.0,
        "N_Q(chi_+) + N_Q(chi_-) - 2 N_Q(chi_out)=" <> fmt[curvature]];
      require["a/r=" <> fmt[aKey] <> " symmetric chi_Q shifts stay on opposite sides of chi_Q^out",
        chiMinus < chiExact < chiPlus,
        "chi_-=" <> fmt[chiMinus] <> ", chi_out=" <> fmt[chiExact] <> ", chi_+=" <> fmt[chiPlus]];
      require["a/r=" <> fmt[aKey] <> " symmetric chi_Q shifts separate the conservative coefficient",
        nqMinus - nqPlus > 0.0,
        "N_Q(chi_-)-N_Q(chi_+)=" <> fmt[nqMinus - nqPlus]];
    ],
    Sort[Keys[perturbationPairs]]
  ];

  Print[""];
  Print["All stage-106 canonical outgoing stress checks passed."];
];

If[
  Catch[
    Print["Loaded config from ", ExpandFileName[configPath]];
    stage106Block[];
    "ok"
  ] === "ok",
  Null,
  Print[""];
  Print["Stage 106 canonical outgoing stress harness failed."];
  Exit[1]
];
