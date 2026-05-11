(* Numerical stress harness for the invariant/orbit closure chain (Stages 168-170). *)

ClearAll["Global`*"];
$HistoryLength = 0;

rootDir = DirectoryName[$InputFileName];
configPath = FileNameJoin[{
   rootDir, "..", "..", "scripts", "numerical",
   "stage185_187_orbit_samples.json"
}];
config = Import[ExpandFileName[configPath], "RawJSON"];

If[config["schema"] =!= "moving_throat_numerical_stage185_187_v1",
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

varOrder = {"lambda_W", "c_etaU", "gamma", "K_U", "K_eta_eff", "K_W_eff", "mu_W", "T_U"};
freeOrder = {"Delta_lambda", "Delta_c", "Delta_gamma", "Delta_U", "Delta_W"};
fullOrder = {"Delta_lambda", "Delta_c", "Delta_gamma", "Delta_U", "Delta_eta", "Delta_W", "Delta_mu", "Delta_T"};

vectorFromMapping[m_, order_List] := N[Lookup[m, order]];

tolerances = config["tolerances"];
stage185Step = N[config["stage185_step"]];
stage186Step = N[config["stage186_step"]];
monomialLinearTol = N[tolerances["monomial_linear_tol"]];
observableLinearTol = N[tolerances["observable_linear_tol"]];
tangentZeroTol = N[tolerances["tangent_zero_tol"]];
kernelTol = N[tolerances["kernel_tol"]];
invariantTol = N[tolerances["invariant_tol"]];
transverseFloor = N[tolerances["transverse_floor"]];

stateVector[case_] := vectorFromMapping[case["reference_state"], varOrder];

derived[state_List] := Module[
  {lam, ceta, gamma, kU, kEta, kW, mu, tU, chi, delta, epsW, eps, zratio, epsEta,
   eStar, fStar, cTr, aTr, bStar, cStar, rTr, tSq},
  {lam, ceta, gamma, kU, kEta, kW, mu, tU} = state;
  chi = gamma ceta/kU;
  delta = tU/kU;
  epsW = gamma^2 lam^2/(kU kW);
  eps = epsW (1 - (2/11) delta/(1 + delta));
  zratio = lam^2 mu/(kEta kW^2);
  epsEta = ceta^2/(kU kEta);
  eStar = 2 epsW/(1 - eps) (11 + 9 delta)/(11 (1 + delta));
  fStar = 2 chi/(1 + delta) + 4 epsW delta/(11 (1 - eps) (1 + delta)^2);
  cTr = chi delta/((1 + chi) (1 + delta) (1 + chi + delta));
  aTr = 2 chi/((1 + chi) (1 + delta));
  bStar = 2 (1 + chi + delta)/delta;
  cStar = ((1 + chi) (1 + delta) (1 + chi + delta))/(chi delta);
  rTr = (1 + chi/(1 + delta))/(1 + chi);
  tSq = zratio (1 + chi)^2/(1 - eps)^2;
  <|
    "chi" -> chi,
    "delta" -> delta,
    "eps_w" -> epsW,
    "eps" -> eps,
    "zratio" -> zratio,
    "eps_eta" -> epsEta,
    "E" -> eStar,
    "F" -> fStar,
    "C_tr" -> cTr,
    "A_tr" -> aTr,
    "B_star" -> bStar,
    "C_star" -> cStar,
    "R_tr" -> rTr,
    "T_sq" -> tSq,
    "T_star" -> rTr^(-cStar),
    "N_star" -> tSq rTr^bStar
  |>
];

monomialInvariants[state_List, ref_Association] := Module[{d, ctr, cnt},
  d = derived[state];
  ctr = d["chi"]^(1 + ref["delta"]) d["delta"]^(1 + ref["chi"]);
  cnt = d["zratio"] d["eps_w"]^ref["E"] d["delta"]^(-ref["F"]);
  <|"Ctr" -> ctr, "Cnt" -> cnt, "eps_eta" -> d["eps_eta"]|>
];

observableInvariants[state_List, ref_Association] := Module[{d},
  d = derived[state];
  <|
    "T_star" -> d["R_tr"]^(-ref["C_star"]),
    "N_star" -> d["T_sq"] d["R_tr"]^ref["B_star"],
    "eps_eta" -> d["eps_eta"]
  |>
];

matrixAt[ref_Association] := {
  {0., 1 + ref["delta"], 1 + ref["delta"], -(2 + ref["chi"] + ref["delta"]), 0., 0., 0., 1 + ref["chi"]},
  {2 (1 + ref["E"]), 0., 2 ref["E"], ref["F"] - ref["E"], -1., -(2 + ref["E"]), 1., -ref["F"]},
  {0., 2., 0., -1., -1., 0., 0., 0.}
};

selectedMinor[ref_Association] := {
  {0., 0., 1 + ref["chi"]},
  {-1., 1., -ref["F"]},
  {-1., 0., 0.}
};

fullTangentFromFree[ref_Association, free_] := Module[
  {lam, ceta, gamma, kU, kW, deta, dt, dmu},
  {lam, ceta, gamma, kU, kW} = vectorFromMapping[free, freeOrder];
  deta = 2 ceta - kU;
  dt = kU - ((1 + ref["delta"])/(1 + ref["chi"])) (gamma + ceta - kU);
  dmu =
    deta
    + 2 kW
    - 2 lam
    - ref["E"] (2 gamma + 2 lam - kU - kW)
    - ref["F"] ((1 + ref["delta"])/(1 + ref["chi"])) (gamma + ceta - kU);
  N[{lam, ceta, gamma, kU, deta, kW, dmu, dt}]
];

linearMonomialDrifts[ref_Association, drift_List] := Module[
  {lam, ceta, gamma, kU, kEta, kW, mu, tU, sigmaChi, sigmaDelta, sigmaEps, sigmaZ, sigmaEta},
  {lam, ceta, gamma, kU, kEta, kW, mu, tU} = drift;
  sigmaChi = gamma + ceta - kU;
  sigmaDelta = tU - kU;
  sigmaEps = 2 gamma + 2 lam - kU - kW;
  sigmaZ = 2 lam + mu - kEta - 2 kW;
  sigmaEta = 2 ceta - kU - kEta;
  <|
    "Sigma_tr" -> (1 + ref["chi"]) sigmaDelta + (1 + ref["delta"]) sigmaChi,
    "Sigma_nt" -> sigmaZ + ref["E"] sigmaEps - ref["F"] sigmaDelta,
    "Sigma_eta" -> sigmaEta
  |>
];

linearObservableDefect[ref_Association, drifts_Association] := <|
  "Theta_1" -> -ref["C_tr"] drifts["Sigma_tr"],
  "Xi_1" -> ref["A_tr"] drifts["Sigma_tr"] + drifts["Sigma_nt"],
  "R1_plus_Xi1" -> -(ref["eps_eta"]/(1 - ref["eps_eta"])) drifts["Sigma_eta"]
|>;

steppedState[base_List, drift_List, step_?NumericQ] := base Exp[step drift];

stage185Compare[case_, refState_List, ref_Association] := Module[
  {refObs, refMon, tangentDrifts, transverseDrifts, stepState, obs, mon, obsDrifts, monDrifts, lin, defects, maxDrift},
  Print[""];
  Print["-- Stage 185: observable vs microscopic drifts --"];
  refObs = observableInvariants[refState, ref];
  refMon = monomialInvariants[refState, ref];

  tangentDrifts = Table[{item["name"], fullTangentFromFree[ref, item]}, {item, case["tangent_free_vectors"]}];
  transverseDrifts = Table[{item["name"], vectorFromMapping[item, fullOrder]}, {item, case["transverse_full_vectors"]}];

  Scan[
    Function[{pair},
      Module[{name, drift},
        {name, drift} = pair;
        stepState = steppedState[refState, drift, stage185Step];
        obs = observableInvariants[stepState, ref];
        mon = monomialInvariants[stepState, ref];
        obsDrifts = <|
          "Sigma_tr" -> Log[obs["T_star"]/refObs["T_star"]]/stage185Step,
          "Sigma_nt" -> Log[obs["N_star"]/refObs["N_star"]]/stage185Step,
          "Sigma_eta" -> Log[obs["eps_eta"]/refObs["eps_eta"]]/stage185Step
        |>;
        monDrifts = <|
          "Sigma_tr" -> Log[mon["Ctr"]/refMon["Ctr"]]/stage185Step,
          "Sigma_nt" -> Log[mon["Cnt"]/refMon["Cnt"]]/stage185Step,
          "Sigma_eta" -> Log[mon["eps_eta"]/refMon["eps_eta"]]/stage185Step
        |>;
        lin = linearMonomialDrifts[ref, drift];
        Print[
          "tangent ", name, ": obs=(", fmt[obsDrifts["Sigma_tr"]], ", ", fmt[obsDrifts["Sigma_nt"]], ", ", fmt[obsDrifts["Sigma_eta"]],
          ") mon=(", fmt[monDrifts["Sigma_tr"]], ", ", fmt[monDrifts["Sigma_nt"]], ", ", fmt[monDrifts["Sigma_eta"]], ")"
        ];
        Scan[
          Function[{key},
            require[name <> " monomial drift matches linear " <> key,
              nearQ[monDrifts[key], lin[key], monomialLinearTol],
              "direct=" <> fmt[monDrifts[key]] <> ", linear=" <> fmt[lin[key]]]
          ],
          {"Sigma_tr", "Sigma_nt", "Sigma_eta"}
        ];
        Scan[
          Function[{key},
            require[name <> " observable drift matches monomial " <> key,
              nearQ[obsDrifts[key], monDrifts[key], observableLinearTol],
              "observable=" <> fmt[obsDrifts[key]] <> ", monomial=" <> fmt[monDrifts[key]]]
          ],
          {"Sigma_tr", "Sigma_nt", "Sigma_eta"}
        ];
        defects = linearObservableDefect[ref, monDrifts];
        Print[
          "  defects: Theta_1=", fmt[defects["Theta_1"]],
          ", Xi_1=", fmt[defects["Xi_1"]],
          ", R_1+Xi_1=", fmt[defects["R1_plus_Xi1"]]
        ];
      ]
    ],
    tangentDrifts
  ];

  Scan[
    Function[{pair},
      Module[{name, drift},
        {name, drift} = pair;
        stepState = steppedState[refState, drift, stage185Step];
        obs = observableInvariants[stepState, ref];
        mon = monomialInvariants[stepState, ref];
        obsDrifts = <|
          "Sigma_tr" -> Log[obs["T_star"]/refObs["T_star"]]/stage185Step,
          "Sigma_nt" -> Log[obs["N_star"]/refObs["N_star"]]/stage185Step,
          "Sigma_eta" -> Log[obs["eps_eta"]/refObs["eps_eta"]]/stage185Step
        |>;
        monDrifts = <|
          "Sigma_tr" -> Log[mon["Ctr"]/refMon["Ctr"]]/stage185Step,
          "Sigma_nt" -> Log[mon["Cnt"]/refMon["Cnt"]]/stage185Step,
          "Sigma_eta" -> Log[mon["eps_eta"]/refMon["eps_eta"]]/stage185Step
        |>;
        lin = linearMonomialDrifts[ref, drift];
        Print[
          "transverse ", name, ": obs=(", fmt[obsDrifts["Sigma_tr"]], ", ", fmt[obsDrifts["Sigma_nt"]], ", ", fmt[obsDrifts["Sigma_eta"]],
          ") mon=(", fmt[monDrifts["Sigma_tr"]], ", ", fmt[monDrifts["Sigma_nt"]], ", ", fmt[monDrifts["Sigma_eta"]], ")"
        ];
        Scan[
          Function[{key},
            require[name <> " monomial drift matches linear " <> key,
              nearQ[monDrifts[key], lin[key], monomialLinearTol],
              "direct=" <> fmt[monDrifts[key]] <> ", linear=" <> fmt[lin[key]]]
          ],
          {"Sigma_tr", "Sigma_nt", "Sigma_eta"}
        ];
        Scan[
          Function[{key},
            require[name <> " observable drift matches monomial " <> key,
              nearQ[obsDrifts[key], monDrifts[key], observableLinearTol],
              "observable=" <> fmt[obsDrifts[key]] <> ", monomial=" <> fmt[monDrifts[key]]]
          ],
          {"Sigma_tr", "Sigma_nt", "Sigma_eta"}
        ];
        defects = linearObservableDefect[ref, monDrifts];
        Print[
          "  defects: Theta_1=", fmt[defects["Theta_1"]],
          ", Xi_1=", fmt[defects["Xi_1"]],
          ", R_1+Xi_1=", fmt[defects["R1_plus_Xi1"]]
        ];
        maxDrift = Max[Abs /@ Lookup[monDrifts, {"Sigma_tr", "Sigma_nt", "Sigma_eta"}]];
        require[name <> " is genuinely transverse", maxDrift > transverseFloor,
          "max monomial drift=" <> fmt[maxDrift]];
      ]
    ],
    transverseDrifts
  ];
];

stage186Check[case_, refState_List, ref_Association] := Module[
  {m, minor, detMinor, svals, minorSvals, condMatrix, condMinor, refMon, drift, residual, maxResidual, stepState, mon, monDrifts, maxDrift},
  Print[""];
  Print["-- Stage 186: tangent-space kernel checks --"];
  m = matrixAt[ref];
  minor = selectedMinor[ref];
  detMinor = N[Det[minor]];
  svals = SingularValueList[N[m]];
  minorSvals = SingularValueList[N[minor]];
  condMatrix = First[svals]/Last[svals];
  condMinor = First[minorSvals]/Last[minorSvals];
  Print["det selected minor = ", fmt[detMinor]];
  Print["cond(M_*) = ", fmt[condMatrix]];
  Print["cond(selected minor) = ", fmt[condMinor]];
  require["selected minor determinant = 1 + chi_*",
    nearQ[detMinor, 1 + ref["chi"], kernelTol],
    "det=" <> fmt[detMinor] <> ", expected=" <> fmt[1 + ref["chi"]]];

  refMon = monomialInvariants[refState, ref];
  Scan[
    Function[{item},
      drift = fullTangentFromFree[ref, item];
      residual = m . drift;
      maxResidual = Max[Abs[residual]];
      require[item["name"] <> " lies in ker(M_*)", maxResidual <= kernelTol,
        "max residual=" <> fmt[maxResidual]];
      stepState = steppedState[refState, drift, stage186Step];
      mon = monomialInvariants[stepState, ref];
      monDrifts = {
        Log[mon["Ctr"]/refMon["Ctr"]]/stage186Step,
        Log[mon["Cnt"]/refMon["Cnt"]]/stage186Step,
        Log[mon["eps_eta"]/refMon["eps_eta"]]/stage186Step
      };
      maxDrift = Max[Abs[monDrifts]];
      require[item["name"] <> " preserves the invariant triple to first order",
        maxDrift <= tangentZeroTol,
        "max drift=" <> fmt[maxDrift]];
    ],
    case["tangent_free_vectors"]
  ];
];

stage187Check[case_, refState_List, ref_Association] := Module[
  {refMon, delta, pairedState, recovered, maxRecover, lam, ceta, gamma, kU, kW,
   detaVal, muVal, dtVal, detaExpected, dtExpected, dmuExpected, pairMon, ratioLogs},
  Print[""];
  Print["-- Stage 187: finite orbit/quotient checks --"];
  refMon = monomialInvariants[refState, ref];
  Scan[
    Function[{item},
      delta = fullTangentFromFree[ref, item];
      pairedState = steppedState[refState, delta, 1.0];
      recovered = Log[pairedState/refState];
      maxRecover = Max[Abs[recovered - delta]];
      require[item["name"] <> " recovers its finite log-ratio vector",
        maxRecover <= invariantTol,
        "max recovery error=" <> fmt[maxRecover]];

      {lam, ceta, gamma, kU, detaVal, kW, muVal, dtVal} = delta;
      detaExpected = 2 ceta - kU;
      dtExpected = kU - ((1 + ref["delta"])/(1 + ref["chi"])) (gamma + ceta - kU);
      dmuExpected =
        detaExpected
        + 2 kW
        - 2 lam
        - ref["E"] (2 gamma + 2 lam - kU - kW)
        - ref["F"] ((1 + ref["delta"])/(1 + ref["chi"])) (gamma + ceta - kU);

      require[item["name"] <> " Delta_eta finite law",
        nearQ[detaVal, detaExpected, invariantTol],
        "direct=" <> fmt[detaVal] <> ", expected=" <> fmt[detaExpected]];
      require[item["name"] <> " Delta_T finite law",
        nearQ[dtVal, dtExpected, invariantTol],
        "direct=" <> fmt[dtVal] <> ", expected=" <> fmt[dtExpected]];
      require[item["name"] <> " Delta_mu finite law",
        nearQ[muVal, dmuExpected, invariantTol],
        "direct=" <> fmt[muVal] <> ", expected=" <> fmt[dmuExpected]];

      pairMon = monomialInvariants[pairedState, ref];
      ratioLogs = <|
        "Ctr" -> Log[pairMon["Ctr"]/refMon["Ctr"]],
        "Cnt" -> Log[pairMon["Cnt"]/refMon["Cnt"]],
        "eps_eta" -> Log[pairMon["eps_eta"]/refMon["eps_eta"]]
      |>;
      Scan[
        Function[{key},
          require[item["name"] <> " preserves " <> key,
            Abs[ratioLogs[key]] <= invariantTol,
            "log ratio=" <> fmt[ratioLogs[key]]]
        ],
        {"Ctr", "Cnt", "eps_eta"}
      ];
    ],
    case["finite_orbit_pairs"]
  ];
];

runCase[case_] := Module[{refState, ref},
  Print[""];
  Print["=== ", case["name"], " ==="];
  refState = stateVector[case];
  ref = derived[refState];
  Print[
    "reference coefficients: chi0*=", fmt[ref["chi"]],
    ", deltaU*=", fmt[ref["delta"]],
    ", epsilon_W*=", fmt[ref["eps_w"]],
    ", epsilon*=", fmt[ref["eps"]],
    ", E*=", fmt[ref["E"]],
    ", F*=", fmt[ref["F"]],
    ", epsilon_eta*=", fmt[ref["eps_eta"]]
  ];
  stage185Compare[case, refState, ref];
  stage186Check[case, refState, ref];
  stage187Check[case, refState, ref];
];

If[
  Catch[
    Print["Loaded config from ", ExpandFileName[configPath]];
    Scan[runCase, config["cases"]];
    Print[""];
    Print["All stage-185/187 orbit-conditioning stress checks passed."];
    "ok"
  ] === "ok",
  Null,
  Print[""];
  Print["Stage 185/187 orbit-conditioning stress harness failed."];
  Exit[1]
];
