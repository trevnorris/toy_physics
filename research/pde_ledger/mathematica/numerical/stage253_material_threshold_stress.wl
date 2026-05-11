(* Numerical stress harness for the Stage 253 material-threshold companion.

   This harness starts from primitive host inputs and verifies the Stage 253
   screening stack forward. It deliberately avoids back-solving host values
   from target Pi margins. *)

ClearAll["Global`*"];
$HistoryLength = 0;

rootDir = DirectoryName[$InputFileName];
configPath = FileNameJoin[{
   rootDir, "..", "..", "scripts", "numerical",
   "stage253_material_threshold_samples.json"
}];
config = Import[ExpandFileName[configPath], "RawJSON"];

If[config["schema"] =!= "moving_throat_numerical_stage253_v2",
  Print["Unexpected config schema."];
  Exit[1];
];

flagOrder = {"Pi_ep", "Pi_chi", "Pi_k", "Pi_t"};

fmt[x_] := ToString @ NumberForm[N[x, 16], {Infinity, 12}, ExponentFunction -> (Null &)];
nearQ[lhs_, rhs_, tol_] := Abs[N[lhs - rhs, 20]] <= tol (1 + Abs[N[rhs, 20]]);

require[label_, condition_, detail_] := Module[{status},
  status = If[TrueQ[condition], "PASS", "FAIL"];
  Print["[", status, "] ", label, ": ", detail];
  If[!TrueQ[condition], Throw[$Failed]]
];

resolveUpsilon[case_, gammaLatSafeEq_] := Module[{mode},
  mode = case["mode"];
  Which[
    mode === "micro", 1.,
    mode === "legacy", gammaLatSafeEq/N[case["gamma_lattice_legacy"]],
    mode === "custom", N[case["Upsilon_lat"]],
    True, Throw[$Failed]
  ]
];

orderedFlags[assoc_] := Association @ Map[# -> TrueQ[assoc[#]] &, flagOrder];
passesThreshold[value_, tol_] := N[value] >= 1.0 - tol;

stage253Block[] := Module[{ratioTol, consistencyTol},
  Print[""];
  Print["=== Stage 253: material-threshold stress ==="];
  ratioTol = N[config["tolerances"]["ratio_tol"]];
  consistencyTol = N[config["tolerances"]["consistency_tol"]];

  Scan[
    Function[{case},
      Module[
        {
          sc, s0, fLat, muEta, zetaEp, tStar, lambdaPhys, lambdaRef, eStar,
          kTurn, kCorr, host, lambdaEpOmegaD, rTurn, rTurnPhys, vprimeTurnAbs,
          kEff, temperature, expected, gammaLatSafeEq, upsilonLat,
          thresholdLambda, tCrossPhys, thresholdProduct, tMax,
          rTurnPhysExpected, kTurnFromVprime, kEffReqDirect, kEffReqReduced,
          piEpDirect, piEpRatio, piChiPhys, piChiReduced, piKDirect, piKRatio,
          piTDirect, piTRatio, actual, gammaLegacy
        },
        sc = N[case["s_c"]];
        s0 = N[case["s_0"]];
        fLat = N[case["f_lat"]];
        muEta = N[case["mu_eta"]];
        zetaEp = N[case["zeta_ep"]];
        tStar = N[case["t_star"]];
        lambdaPhys = N[case["lambda_phys"]];
        lambdaRef = N[case["lambda_ref"]];
        eStar = N[case["E_star"]];
        kTurn = N[case["K_turn"]];
        kCorr = N[case["K_corr"]];

        host = case["host"];
        lambdaEpOmegaD = N[host["lambda_ep_omega_D"]];
        rTurn = N[host["r_turn"]];
        rTurnPhys = N[host["r_turn_phys"]];
        vprimeTurnAbs = N[host["Vprime_turn_abs"]];
        kEff = N[host["k_eff"]];
        temperature = N[host["temperature"]];

        expected = orderedFlags[case["expected_flags"]];

        gammaLatSafeEq = fLat muEta (s0^2 - sc^2)/sc;
        upsilonLat = resolveUpsilon[case, gammaLatSafeEq];
        thresholdLambda = gammaLatSafeEq/(upsilonLat zetaEp tStar);
        tCrossPhys = tStar/sc;
        thresholdProduct = thresholdLambda tCrossPhys;
        tMax = kCorr/tCrossPhys;

        rTurnPhysExpected = lambdaPhys rTurn/lambdaRef;
        kTurnFromVprime = lambdaRef^2 vprimeTurnAbs/rTurn;
        kEffReqDirect = eStar lambdaRef vprimeTurnAbs/(lambdaPhys rTurnPhys);
        kEffReqReduced = kTurn eStar/lambdaPhys^2;

        piEpDirect = upsilonLat zetaEp lambdaEpOmegaD tStar/gammaLatSafeEq;
        piEpRatio = lambdaEpOmegaD/thresholdLambda;
        piChiPhys = 2 lambdaPhys/rTurnPhys;
        piChiReduced = 2 lambdaRef/rTurn;
        piKDirect = kEff lambdaPhys^2/(kTurn eStar);
        piKRatio = kEff/kEffReqReduced;
        piTDirect = kCorr/(temperature tCrossPhys);
        piTRatio = tMax/temperature;

        actual = Association[
          "Pi_ep" -> passesThreshold[piEpDirect, ratioTol],
          "Pi_chi" -> passesThreshold[piChiPhys, ratioTol],
          "Pi_k" -> passesThreshold[piKDirect, ratioTol],
          "Pi_t" -> passesThreshold[piTDirect, ratioTol]
        ];

        Print[""];
        Print[
          case["name"], " (", case["mode"], "): gamma_safe=", fmt[gammaLatSafeEq],
          ", Upsilon=", fmt[upsilonLat], ", lambda_min=", fmt[thresholdLambda]
        ];
        Print[
          "  host inputs: lambda_ep*omega_D=", fmt[lambdaEpOmegaD],
          ", r_turn=", fmt[rTurn], ", r_turn_phys=", fmt[rTurnPhys],
          ", |V'|=", fmt[vprimeTurnAbs], ", k_eff=", fmt[kEff],
          ", T=", fmt[temperature]
        ];
        Print[
          "  computed margins: Pi_ep=", fmt[piEpDirect],
          ", Pi_chi=", fmt[piChiPhys], ", Pi_k=", fmt[piKDirect],
          ", Pi_t=", fmt[piTDirect]
        ];

        require[
          case["name"] <> " primitive admissibility",
          s0 > sc > 0 && fLat > 0 && muEta > 0 && zetaEp > 0 && upsilonLat > 0,
          "s0=" <> fmt[s0] <> ", s_c=" <> fmt[sc] <> ", f_lat=" <> fmt[fLat] <>
            ", mu_eta=" <> fmt[muEta] <> ", zeta_ep=" <> fmt[zetaEp] <>
            ", Upsilon=" <> fmt[upsilonLat]
        ];
        require[
          case["name"] <> " host-input positivity",
          lambdaEpOmegaD > 0 && rTurn > 0 && rTurnPhys > 0 && vprimeTurnAbs > 0 &&
            kEff > 0 && temperature > 0,
          "lambda_ep*omega_D=" <> fmt[lambdaEpOmegaD] <> ", r_turn=" <> fmt[rTurn] <>
            ", r_turn_phys=" <> fmt[rTurnPhys] <> ", |V'|=" <> fmt[vprimeTurnAbs] <>
            ", k_eff=" <> fmt[kEff] <> ", T=" <> fmt[temperature]
        ];
        require[
          case["name"] <> " turnover threshold positivity",
          gammaLatSafeEq > 0 && thresholdLambda > 0 && thresholdProduct > 0,
          "gamma_safe=" <> fmt[gammaLatSafeEq] <> ", lambda_min=" <> fmt[thresholdLambda] <>
            ", product_min=" <> fmt[thresholdProduct]
        ];
        require[
          case["name"] <> " turning-radius dictionary",
          nearQ[rTurnPhys, rTurnPhysExpected, consistencyTol],
          "r_turn_phys=" <> fmt[rTurnPhys] <> ", expected=" <> fmt[rTurnPhysExpected]
        ];
        require[
          case["name"] <> " K_turn / Vprime consistency",
          nearQ[kTurnFromVprime, kTurn, consistencyTol],
          "K_turn_from_Vprime=" <> fmt[kTurnFromVprime] <> ", K_turn=" <> fmt[kTurn]
        ];
        require[
          case["name"] <> " stiffness compiler consistency",
          nearQ[kEffReqDirect, kEffReqReduced, consistencyTol],
          "k_req_direct=" <> fmt[kEffReqDirect] <> ", k_req_reduced=" <> fmt[kEffReqReduced]
        ];
        require[
          case["name"] <> " product-threshold factorization",
          nearQ[thresholdProduct, gammaLatSafeEq/(upsilonLat zetaEp sc), consistencyTol],
          "product=" <> fmt[thresholdProduct] <> ", reduced=" <>
            fmt[gammaLatSafeEq/(upsilonLat zetaEp sc)]
        ];
        require[
          case["name"] <> " Pi_ep direct vs threshold ratio",
          nearQ[piEpDirect, piEpRatio, ratioTol],
          "direct=" <> fmt[piEpDirect] <> ", ratio=" <> fmt[piEpRatio]
        ];
        require[
          case["name"] <> " Pi_chi physical vs reduced ratio",
          nearQ[piChiPhys, piChiReduced, ratioTol],
          "physical=" <> fmt[piChiPhys] <> ", reduced=" <> fmt[piChiReduced]
        ];
        require[
          case["name"] <> " Pi_k direct vs reduced ratio",
          nearQ[piKDirect, piKRatio, ratioTol],
          "direct=" <> fmt[piKDirect] <> ", ratio=" <> fmt[piKRatio]
        ];
        require[
          case["name"] <> " Pi_t direct vs Tmax ratio",
          nearQ[piTDirect, piTRatio, ratioTol],
          "direct=" <> fmt[piTDirect] <> ", ratio=" <> fmt[piTRatio]
        ];
        require[
          case["name"] <> " screening classification",
          actual === expected,
          "actual=" <> ToString[actual, InputForm] <> ", expected=" <> ToString[expected, InputForm]
        ];
        require[
          case["name"] <> " overall screen verdict",
          (And @@ Values[actual]) === (And @@ Values[expected]),
          "actual=" <> ToString[actual, InputForm] <> ", expected=" <> ToString[expected, InputForm]
        ];

        If[case["mode"] === "legacy",
          gammaLegacy = N[case["gamma_lattice_legacy"]];
          require[
            case["name"] <> " legacy Upsilon recovery",
            nearQ[upsilonLat, gammaLatSafeEq/gammaLegacy, consistencyTol],
            "Upsilon=" <> fmt[upsilonLat] <> ", gamma_safe/gamma_legacy=" <> fmt[gammaLatSafeEq/gammaLegacy]
          ];
          require[
            case["name"] <> " legacy lambda threshold",
            nearQ[thresholdLambda, gammaLegacy/(zetaEp tStar), consistencyTol],
            "threshold=" <> fmt[thresholdLambda] <> ", legacy=" <> fmt[gammaLegacy/(zetaEp tStar)]
          ];
          require[
            case["name"] <> " legacy product threshold",
            nearQ[thresholdProduct, gammaLegacy/(zetaEp sc), consistencyTol],
            "product=" <> fmt[thresholdProduct] <> ", legacy=" <> fmt[gammaLegacy/(zetaEp sc)]
          ];
        ];
      ]
    ],
    config["cases"]
  ];
];

If[
  Catch[
    Print["Loaded config from ", ExpandFileName[configPath]];
    stage253Block[];
    Print[""];
    Print["All Stage 253 numerical stress checks passed."];
    "ok"
  ] === "ok",
  Null,
  Print[""];
  Print["Stage 253 numerical stress harness failed."];
  Exit[1]
];
