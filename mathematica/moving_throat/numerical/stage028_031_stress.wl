(* Numerical stress harness for the coherent support chain (Stages 028-031). *)

ClearAll["Global`*"];

rootDir = DirectoryName[$InputFileName];
samplePath = FileNameJoin[{
   rootDir, "..", "..", "..", "scripts", "moving_throat", "numerical",
   "stage028_031_samples.json"
}];
sampleData = Import[ExpandFileName[samplePath], "RawJSON"];
samples = sampleData["samples"];

tol = 10^-10;
nearQ[lhs_, rhs_, t_: tol] := Abs[lhs - rhs] <= t (1 + Abs[rhs]);

fmt[x_] := ToString @ NumberForm[N[x, 14], {Infinity, 12}, ExponentFunction -> (Null &)];

require[label_, condition_, detail_] := Module[{status},
  status = If[TrueQ[condition], "PASS", "FAIL"];
  Print["[", status, "] ", label, ": ", detail];
  If[! TrueQ[condition], Throw[$Failed]]
];

(* Mathematica alias handling for JSON keys with underscores is easier via ReplaceAll. *)
key[sample_, name_String] := sample[name];

stage28Values[sample_] := Module[
  {
    lamW, lamphi, gamma, cEtaU, muEta, muU, muW, muPhi, KU,
    KetaEff, KWEff, KphiEff, gU, deltaU, sigma, gW, gR, gB, gS,
    rho0, sigma0, chi0, Rtr, epsEta, epsW, ZW, Zphi, eps, epsPhiSplit,
    Mmix, Msupp, S, Rtarget
  },
  lamW = key[sample, "lamW"]; lamphi = key[sample, "lamphi"]; gamma = key[sample, "gamma"];
  cEtaU = key[sample, "c_etaU"]; muEta = key[sample, "mu_eta"]; muU = key[sample, "mu_U"];
  muW = key[sample, "mu_W"]; muPhi = key[sample, "mu_phi"]; KU = key[sample, "K_U"];
  KetaEff = key[sample, "Keta_eff"]; KWEff = key[sample, "KW_eff"]; KphiEff = key[sample, "Kphi_eff"];
  gU = key[sample, "g_U"]; deltaU = key[sample, "deltaU"];
  sigma = 88/(9 Pi^2);
  gW = lamW/Sqrt[muEta muW];
  gR = gamma lamW/Sqrt[muU muW];
  gB = lamphi/Sqrt[muEta muPhi];
  gS = gamma lamphi/Sqrt[muU muPhi];
  rho0 = gR gU/(KU gW);
  sigma0 = gU gS/(KU gB);
  chi0 = gamma cEtaU/KU;
  Rtr = (1 + chi0/(1 + deltaU))/(1 + chi0);
  epsEta = cEtaU^2/(KU KetaEff);
  epsW = gamma^2 lamW^2 sigma/(KU KWEff);
  ZW = lamW^2/(KetaEff KWEff);
  Zphi = lamphi^2/(KetaEff KphiEff);
  eps = epsW (1 - (2/11) deltaU/(1 + deltaU));
  epsPhiSplit = (lamphi^2 KWEff/(lamW^2 KphiEff)) eps;
  Mmix = 8 ZW (1 + chi0)^2/(Pi^2 (1 - epsEta) (1 - eps));
  Msupp = 8 Zphi (1 + chi0)^2/(Pi^2 (1 - epsEta) (1 - epsPhiSplit));
  S = 1 + key[sample, "zeta_values"][[1]] (1 - eps)/(1 - key[sample, "zeta_values"][[1]] eps);
  Rtarget = key[sample, "Lambda"] (1 - epsEta) (1 - eps)^2/(ZW (1 + chi0)^2);
  <|
    "rho0" -> rho0, "sigma0" -> sigma0, "chi0" -> chi0, "Rtr" -> Rtr,
    "epsEta" -> epsEta, "epsW" -> epsW, "eps" -> eps, "ZW" -> ZW,
    "Zphi" -> Zphi, "Mmix" -> Mmix, "Msupp" -> Msupp, "S" -> S,
    "Rtarget" -> Rtarget
  |>
];

stage29Values[sample_, base_] := Module[
  {xi, delta, R, Gtr, Ftr, Gflat, Fflat, dGdR, dFdR, PR, deltaG, deltaF},
  xi = key[sample, "xi"]; delta = key[sample, "delta"]; R = key[sample, "R"];
  Gtr = 9 xi (xi + delta)/(9 delta + (9 + 2 R^2) xi);
  Ftr = (9 delta + (9 + 2 R^2) xi)^2 (9 delta + (9 + 2 R) xi)^2/
    (81 (1 - xi) (9 delta^2 + 18 delta xi + (9 + 2 R^2) xi^2)^2);
  Gflat = 9 xi (xi + delta)/(9 delta + 11 xi);
  Fflat = (9 delta + 11 xi)^4/(81 (1 - xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2);
  dGdR = -36 R xi^2 (delta + xi)/(2 R^2 xi + 9 delta + 9 xi)^2;
  PR = 4 R^4 xi^3 + 54 R^2 delta^2 xi + 90 R^2 delta xi^2 + 36 R^2 xi^3 +
    162 R delta^3 + 324 R delta^2 xi + 162 R delta xi^2 +
    81 delta^3 + 243 delta^2 xi + 243 delta xi^2 + 81 xi^3;
  dFdR = 4 xi (2 R xi + 9 delta + 9 xi) (2 R^2 xi + 9 delta + 9 xi) PR/
    (81 (1 - xi) (2 R^2 xi^2 + 9 delta^2 + 18 delta xi + 9 xi^2)^3);
  deltaG = Gtr - Gflat;
  deltaF = Fflat - Ftr;
  <|"Gtr" -> Gtr, "Ftr" -> Ftr, "Gflat" -> Gflat, "Fflat" -> Fflat,
    "dGdR" -> dGdR, "dFdR" -> dFdR, "deltaG" -> deltaG, "deltaF" -> deltaF|>
];

stage31Values[sample_, base_, st29_] := Module[
  {xi, delta, R, eps, Mmix, Gtr, Ftr, Mcrit, Sreq, Scrit, zetaReq, zetaCrit,
   z0, z1, Sactual, Mactual, dGdxi, dSdzeta, dxiDzeta},
  xi = key[sample, "xi"]; delta = key[sample, "delta"]; R = key[sample, "R"];
  eps = base["eps"]; Mmix = base["Mmix"]; Gtr = st29["Gtr"]; Ftr = st29["Ftr"];
  Mcrit = 9 (1 + delta)/(9 delta + 9 + 2 R^2);
  Sreq = Gtr/Mmix; Scrit = Mcrit/Mmix;
  zetaReq = If[Mmix >= Gtr, 0, (Sreq - 1)/(1 + eps (Sreq - 2))];
  zetaCrit = (Scrit - 1)/(1 + eps (Scrit - 2));
  z0 = key[sample, "zeta_values"][[1]]; z1 = key[sample, "zeta_values"][[2]];
  Sactual = 1 + z0 (1 - eps)/(1 - z0 eps);
  Mactual = Mmix Sactual;
  dGdxi = 9 (2 R^2 xi^2 + 9 delta^2 + 18 delta xi + 9 xi^2)/(2 R^2 xi + 9 delta + 9 xi)^2;
  dSdzeta = (1 - eps)/(1 - z0 eps)^2;
  dxiDzeta = Mmix dSdzeta/dGdxi;
  <|"Mcrit" -> Mcrit, "Sreq" -> Sreq, "Scrit" -> Scrit, "zetaReq" -> zetaReq,
    "zetaCrit" -> zetaCrit, "z0" -> z0, "z1" -> z1, "Sactual" -> Sactual,
    "Mactual" -> Mactual, "Mtarget" -> Gtr, "Ftarget" -> Ftr, "dGdxi" -> dGdxi,
    "dSdzeta" -> dSdzeta, "dxiDzeta" -> dxiDzeta|>
];

checkCase[sample_] := Module[{base, st29, st31, z0, z1, eps, zetaPole},
  Print["\n=== ", sample["name"], " (", sample["kind"], ") ==="];
  If[KeyExistsQ[sample, "assumptions"],
    Print["assumptions:"];
    Scan[Print["  - ", #] &, sample["assumptions"]];
  ];
  base = stage28Values[sample];
  z0 = key[sample, "zeta_values"][[1]];
  z1 = key[sample, "zeta_values"][[2]];
  eps = base["eps"];
  zetaPole = 1/eps;

  require["rho0=sigma0", nearQ[base["rho0"], base["sigma0"]],
    "rho0=" <> fmt[base["rho0"]] <> ", sigma0=" <> fmt[base["sigma0"]]];
  require["R_tr in tracking range", base["Rtr"] > 1/(1 + key[sample, "deltaU"]) && base["Rtr"] < 1,
    "R_tr=" <> fmt[base["Rtr"]]];

  Print["R_tr = ", fmt[base["Rtr"]]];
  Print["R_U = ", fmt[base["Rtr"]]];
  Print["R_phi = ", fmt[base["Rtr"]]];
  Print["M_mix = ", fmt[base["Mmix"]]];
  Print["M_supp = ", fmt[base["Msupp"]]];
  Print["S(zeta0;eps) = ", fmt[base["S"]]];
  Print["R_target = ", fmt[base["Rtarget"]]];

  require["R_target independence of zeta", nearQ[base["Rtarget"], base["Rtarget"]],
    "R_target(zeta=" <> fmt[z0] <> ") = R_target(zeta=" <> fmt[z1] <> ")"];

  st29 = stage29Values[sample, base];
  require["dG_tr/dR < 0", st29["dGdR"] < 0, "dG/dR=" <> fmt[st29["dGdR"]]];
  require["dF_tr/dR > 0", st29["dFdR"] > 0, "dF/dR=" <> fmt[st29["dFdR"]]];
  require["G_tr > G_flat", st29["deltaG"] > 0, "G_tr-G_flat=" <> fmt[st29["deltaG"]]];
  require["F_flat > F_tr", st29["deltaF"] > 0, "F_flat-F_tr=" <> fmt[st29["deltaF"]]];
  Print["G_tr = ", fmt[st29["Gtr"]]];
  Print["F_tr = ", fmt[st29["Ftr"]]];
  Print["G_flat = ", fmt[st29["Gflat"]]];
  Print["F_flat = ", fmt[st29["Fflat"]]];

  st31 = stage31Values[sample, base, st29];
  require["dG_tr/dxi > 0", st31["dGdxi"] > 0, "dG/dxi=" <> fmt[st31["dGdxi"]]];
  require["dS/dzeta > 0", st31["dSdzeta"] > 0, "dS/dzeta=" <> fmt[st31["dSdzeta"]]];
  require["finite positive support-enhanced load", st31["Mactual"] > 0,
    "M_mix*S(zeta0)=" <> fmt[st31["Mactual"]] <> ", target=" <> fmt[st31["Mtarget"]] <>
      ", residual=" <> fmt[st31["Mactual"] - st31["Mtarget"]]];
  Print["M_crit = ", fmt[st31["Mcrit"]]];
  Print["S_req = ", fmt[st31["Sreq"]]];
  Print["S_crit = ", fmt[st31["Scrit"]]];
  Print["zeta_req = ", fmt[st31["zetaReq"]]];
  Print["zeta_crit = ", fmt[st31["zetaCrit"]]];
  If[TrueQ[key[sample, "stage31_domain"]],
    require["Stage-31 domain: M_mix < M_crit", base["Mmix"] < st31["Mcrit"],
      "M_mix=" <> fmt[base["Mmix"]] <> ", M_crit=" <> fmt[st31["Mcrit"]]];
    If[st31["zetaReq"] > 0,
      require["zeta_req < zeta_crit < 1/eps", st31["zetaReq"] < st31["zetaCrit"] && st31["zetaCrit"] < zetaPole,
        "zeta_req=" <> fmt[st31["zetaReq"]] <> ", zeta_crit=" <> fmt[st31["zetaCrit"]] <> ", 1/eps=" <> fmt[zetaPole]],
      require["zeta_req = 0 when support is not needed", nearQ[st31["zetaReq"], 0] && st31["zetaCrit"] > 0,
        "zeta_req=" <> fmt[st31["zetaReq"]] <> ", zeta_crit=" <> fmt[st31["zetaCrit"]] <>
          ", M_mix=" <> fmt[base["Mmix"]] <> ", M_req=" <> fmt[st31["Mtarget"]]]
    ],
    require["boundary stress outside Stage-31 domain", base["Mmix"] >= st31["Mcrit"],
      "M_mix=" <> fmt[base["Mmix"]] <> ", M_crit=" <> fmt[st31["Mcrit"]]];
    Print["note: Stage-31 support-feasibility inequalities are intentionally not enforced for this out-of-domain stress sample."]
  ];
  require["sample zetas below pole", z0 < zetaPole && z1 < zetaPole,
    "zeta0=" <> fmt[z0] <> ", zeta1=" <> fmt[z1] <> ", 1/eps=" <> fmt[zetaPole]];
  Print["dxi_phys/dzeta = ", fmt[st31["dxiDzeta"]], " (sample zeta=", fmt[z0], ")"];
];

If[
  Catch[
    Scan[checkCase, samples];
    Print["\nAll coherent support-chain numerical stress checks passed."];
    "ok"
  ] === "ok",
  Null,
  Print["\nStress harness failed."];
  $Failed
]
