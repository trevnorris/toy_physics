(* Foundational resonance/outgoing stress checks for Stages 003 and 021. *)

ClearAll["Global`*"];
$HistoryLength = 0;

rootDir = DirectoryName[$InputFileName];
configPath = FileNameJoin[{
   rootDir, "..", "..", "scripts", "numerical",
   "stage003_021_foundational_samples.json"
}];
config = Import[ExpandFileName[configPath], "RawJSON"];

If[config["schema"] =!= "moving_throat_numerical_stage003_021_v1",
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

stage003Block[] := Module[{tol, cases},
  Print[""];
  Print["=== Stage 003: one-mode resonance stress ==="];
  tol = config["stage003"]["tolerances"];
  cases = config["stage003"]["cases"];
  Scan[
    Function[{case},
      Module[{m, om2, delta, g, varpi2, disc, minus, plus, approxMinus, approxPlus, wallShift, matterShift, lambdaRatio, wallRel, matterRel, tolRel, floor},
        m = N[case["M"]];
        om2 = N[case["Omega_eta2"]];
        delta = N[case["delta"]];
        g = N[case["g"]];
        varpi2 = om2 + delta;
        disc = Sqrt[(om2 - varpi2)^2 + 4 g^2/m];
        minus = (om2 + varpi2 - disc)/2;
        plus = (om2 + varpi2 + disc)/2;
        approxMinus = om2 - g^2/(m delta);
        approxPlus = varpi2 + g^2/(m delta);
        wallShift = minus - om2;
        matterShift = plus - varpi2;
        lambdaRatio = g^2/(m delta^2);

        Print[""];
        Print[case["name"], " (", case["kind"], "): M=", fmt[m], ", Omega_eta^2=", fmt[om2], ", delta=", fmt[delta], ", g=", fmt[g]];
        Print["  resonance parameter g^2/(M delta^2) = ", fmt[lambdaRatio]];

        require[case["name"] <> " exact root ordering",
          minus < om2 < varpi2 < plus,
          "omega_-^2=" <> fmt[minus] <> ", Omega_eta^2=" <> fmt[om2] <> ", varpi^2=" <> fmt[varpi2] <> ", omega_+^2=" <> fmt[plus]];
        require[case["name"] <> " exact sum of roots",
          nearQ[minus + plus, om2 + varpi2, tol["sum_product_tol"]],
          "sum=" <> fmt[minus + plus] <> ", target=" <> fmt[om2 + varpi2]];
        require[case["name"] <> " exact product of roots",
          nearQ[minus plus, om2 varpi2 - g^2/m, tol["sum_product_tol"]],
          "product=" <> fmt[minus plus] <> ", target=" <> fmt[om2 varpi2 - g^2/m]];
        require[case["name"] <> " wall/matter shift signs",
          wallShift < 0 && matterShift > 0,
          "wall_shift=" <> fmt[wallShift] <> ", matter_shift=" <> fmt[matterShift]];

        wallRel = Abs[wallShift - (approxMinus - om2)]/Abs[approxMinus - om2];
        matterRel = Abs[matterShift - (approxPlus - varpi2)]/Abs[approxPlus - varpi2];
        Print[
          "  exact shifts: wall=", fmt[wallShift], ", matter=", fmt[matterShift],
          "; perturbative: wall=", fmt[approxMinus - om2], ", matter=", fmt[approxPlus - varpi2]
        ];
        If[case["regime"] === "perturbative_valid",
          tolRel = N[case["perturbative_rel_tol"]];
          require[case["name"] <> " wall-like perturbative validity",
            wallRel <= tolRel,
            "relative error=" <> fmt[wallRel] <> ", tol=" <> fmt[tolRel]];
          require[case["name"] <> " matter-like perturbative validity",
            matterRel <= tolRel,
            "relative error=" <> fmt[matterRel] <> ", tol=" <> fmt[tolRel]];
          ,
          floor = N[case["breakdown_floor"]];
          require[case["name"] <> " perturbative breakdown is visible",
            wallRel >= floor && matterRel >= floor,
            "wall_rel=" <> fmt[wallRel] <> ", matter_rel=" <> fmt[matterRel] <> ", floor=" <> fmt[floor]]
        ];
      ]
    ],
    cases
  ];
];

stage021Block[] := Module[{cases},
  Print[""];
  Print["=== Stage 021: outgoing mixed-sector stress ==="];

  h1[z_] := (Sin[z]/z^2 - Cos[z]/z) + I (-Cos[z]/z^2 - Sin[z]/z);
  h2[z_] := (((3/z^3) - 1/z) Sin[z] - 3 Cos[z]/z^2) + I (-((3/z^3) - 1/z) Cos[z] - 3 Sin[z]/z^2);
  y2Hat[radius_, cS_, omega_] := Module[{z, hh2, hh2p, lam, y2},
    z = omega radius/cS;
    hh2 = h2[z];
    hh2p = h1[z] - 3 hh2/z;
    lam = (omega/cS) hh2p/hh2;
    y2 = 1/lam;
    y2/(-radius/3)
  ];

  cases = config["stage021"]["cases"];
  Scan[
    Function[{case},
      Module[{oA, oW, r, gA, gW, radius, cS, denom0, n0, gamma5, aker, wker, delta, nExact, sigmaCons, piOut, sigmaFull, dCorr, y2Exact, y2Even, nRel, dcorrRel, y2ImagRel, y2RealAbs, omega},
        oA = N[case["Omega_A"]];
        oW = N[case["Omega_W"]];
        r = N[case["R"]];
        gA = N[case["g_A"]];
        gW = N[case["g_W"]];
        radius = N[case["radius"]];
        cS = N[case["c_s"]];
        denom0 = oA^2 oW^2 - r^2;
        n0 = (oA^2 gW + r gA)^2/denom0^2;
        gamma5 = radius^5/(27 cS^5);

        Print[""];
        Print[case["name"], " (", case["kind"], "): Omega_A=", fmt[oA], ", Omega_W=", fmt[oW], ", R=", fmt[r]];
        Print["  denominator Omega_A^2 Omega_W^2 - R^2 = ", fmt[denom0]];
        require[case["name"] <> " positive mixed denominator", denom0 > 0, "denominator=" <> fmt[denom0]];

        Scan[
          Function[{omegaRaw},
            omega = N[omegaRaw];
            aker = oA^2 - omega^2;
            wker = oW^2 - omega^2;
            delta = aker wker - r^2;
            nExact = (aker gW + r gA)^2/delta^2;
            sigmaCons = (gA^2 wker + 2 gA gW r + gW^2 aker)/delta;
            piOut = I gamma5 omega^5;
            sigmaFull = (gA^2 (wker - piOut) + 2 gA gW r + gW^2 aker)/(aker (wker - piOut) - r^2);
            dCorr = -(sigmaFull - sigmaCons);
            y2Exact = y2Hat[radius, cS, omega];
            y2Even = 1 + radius^2 omega^2/(9 cS^2) + 4 radius^4 omega^4/(81 cS^4);
            nRel = Abs[nExact - n0]/(1 + Abs[n0]);
            dcorrRel = Abs[(Im[dCorr]/omega^5) - (-gamma5 n0)]/(1 + Abs[gamma5 n0]);
            y2ImagRel = Abs[(Im[y2Exact]/omega^5) - gamma5]/(1 + Abs[gamma5]);
            y2RealAbs = Abs[Re[y2Exact] - y2Even];

            Print[
              "  omega=", fmt[omega], ": N(omega)=", fmt[nExact],
              ", Im(delta D)/omega^5=", fmt[Im[dCorr]/omega^5],
              ", Im(Y2hat)/omega^5=", fmt[Im[y2Exact]/omega^5]
            ];
            require[case["name"] <> " N0 convergence at omega=" <> fmt[omega],
              nRel <= N[case["N0_rel_tol"]],
              "relative error=" <> fmt[nRel] <> ", tol=" <> fmt[N[case["N0_rel_tol"]]]];
            require[case["name"] <> " wall odd coefficient at omega=" <> fmt[omega],
              dcorrRel <= N[case["Dcorr_rel_tol"]],
              "relative error=" <> fmt[dcorrRel] <> ", tol=" <> fmt[N[case["Dcorr_rel_tol"]]]];
            require[case["name"] <> " Y2hat odd coefficient at omega=" <> fmt[omega],
              y2ImagRel <= N[case["Y2_imag_rel_tol"]],
              "relative error=" <> fmt[y2ImagRel] <> ", tol=" <> fmt[N[case["Y2_imag_rel_tol"]]]];
            require[case["name"] <> " Y2hat even branch at omega=" <> fmt[omega],
              y2RealAbs <= N[case["Y2_real_abs_tol"]],
              "abs error=" <> fmt[y2RealAbs] <> ", tol=" <> fmt[N[case["Y2_real_abs_tol"]]]];
          ],
          case["omegas"]
        ];
      ]
    ],
    cases
  ];
];

If[
  Catch[
    Print["Loaded config from ", ExpandFileName[configPath]];
    stage003Block[];
    stage021Block[];
    Print[""];
    Print["All stage-003/021 foundational stress checks passed."];
    "ok"
  ] === "ok",
  Null,
  Print[""];
  Print["Stage-003/021 foundational stress harness failed."];
  Exit[1]
];
