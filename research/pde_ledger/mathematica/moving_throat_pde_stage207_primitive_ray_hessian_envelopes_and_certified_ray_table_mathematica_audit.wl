ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

cleanScalar[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

cleanTensor[expr_] := If[
  ListQ[expr],
  Map[cleanScalar, expr, {ArrayDepth[expr]}],
  cleanScalar[expr]
];

zeroTensorQ[expr_] := And @@ (TrueQ[# === 0] & /@ Flatten[{expr}]);
pretty[expr_] := If[ListQ[expr], MatrixForm[expr], fmt[expr]];

expectZero[name_String, expr_] := Module[{res},
  res = cleanTensor[expr];
  If[ListQ[res],
    Print[name, " ="];
    Print[pretty[res]],
    Print[name, " = ", fmt[res]]
  ];
  If[zeroTensorQ[res], pass[name], fail[name, res]];
];

banner["STAGE 207 -- PRIMITIVE-RAY HESSIAN ENVELOPES AND CERTIFIED RAY TABLE"];

Clear[
  Hll, Hcc, Hgg, HUU, HWW,
  Hlc, Hlg, HlU, HlW, Hcg, HcU, HcW, HgU, HgW, HUW,
  a, b, Gam, H0, k, kLambda, kC, kGamma, kU, kW,
  kappaLo, kappaHi, kappaTurnLo, kappaTurnHi, tau,
  chi0Star, deltaUStar, eStar, fStar,
  sLambda, sC, sGamma, sU, sW,
  epsLambda, epsC, epsGamma, epsU, epsW
];

$Assumptions = (
  Element[
    {
      Hll, Hcc, Hgg, HUU, HWW,
      Hlc, Hlg, HlU, HlW, Hcg, HcU, HcW, HgU, HgW, HUW,
      a, b, Gam, kappaLo, kappaHi, tau,
      sLambda, sC, sGamma, sU, sW,
      epsLambda, epsC, epsGamma, epsU, epsW
    },
    Reals
  ] &&
  H0 > 0 && k > 0 &&
  kLambda > 0 && kC > 0 && kGamma > 0 && kU > 0 && kW > 0 &&
  k^2 - 2 kappaLo H0 > 0 && k^2 - 2 kappaHi H0 > 0 &&
  kappaTurnLo < 0 && kappaTurnHi < 0 &&
  chi0Star > 0 && deltaUStar > 0 && eStar > 0 && fStar > 0 &&
  Gam != 0
);

labels = {"lambda", "c", "gamma", "U", "W"};
H = {
  {Hll, Hlc, Hlg, HlU, HlW},
  {Hlc, Hcc, Hcg, HcU, HcW},
  {Hlg, Hcg, Hgg, HgU, HgW},
  {HlU, HcU, HgU, HUU, HUW},
  {HlW, HcW, HgW, HUW, HWW}
};

subbanner["M1. Diagonal Hessian reduction on primitive axes"];
Do[
  Module[{e = UnitVector[5, i], diag = H[[i, i]], label = labels[[i]]},
    expectZero["M1 +e_" <> label <> " diagonal residual", e . H . e - diag];
    expectZero["M1 -e_" <> label <> " diagonal residual", (-e) . H . (-e) - diag];
  ],
  {i, 1, 5}
];

subbanner["M2. Mixed-ray quadratic form and off-diagonal coefficient"];
Do[
  Module[{eI, eJ, ray, form, expected, label},
    eI = UnitVector[5, i];
    eJ = UnitVector[5, j];
    ray = a eI + b eJ;
    form = Expand[ray . H . ray];
    expected = a^2 H[[i, i]] + 2 a b H[[i, j]] + b^2 H[[j, j]];
    label = labels[[i]] <> "," <> labels[[j]];
    expectZero["M2 quadratic form residual (" <> label <> ")", form - expected];
    expectZero[
      "M2 cross coefficient residual (" <> label <> ")",
      Coefficient[form, a b] - 2 H[[i, j]]
    ];
  ],
  {i, 1, 4}, {j, i + 1, 5}
];

subbanner["M3. Canonical orientation law"];
epsilonCanonical = -Sign[Gam];
orientedSlope = epsilonCanonical Gam;
expectZero["M3 K + Abs[Gamma]", orientedSlope + Abs[Gam]];

subbanner["M4. Certified monotone bracket root map"];
comparisonPolynomial[c_, x_] := H0 - k x + 1/2 c x^2;
closedMonotoneRoot[c_] := 2 H0/(k + Sqrt[k^2 - 2 c H0]);

monotoneFirstRoot[label_String, cSym_] := Module[
  {ass, roots, positiveBranch, negativeBranch, linearBranch, positiveSet},
  ass = $Assumptions && Element[cSym, Reals] && k^2 - 2 cSym H0 > 0;
  roots = stripConditional[tau /. Solve[comparisonPolynomial[cSym, tau] == 0, tau]];
  positiveSet = FullSimplify[
    Reduce[comparisonPolynomial[cSym, tau] == 0 && tau > 0, tau, Reals],
    Assumptions -> ass
  ];
  Print["M4 positive root set (" <> label <> ") = ", fmt[positiveSet]];
  positiveBranch = SelectFirst[
    roots,
    Function[root,
      Module[{other = First[DeleteCases[roots, root]]},
        TrueQ[
          FullSimplify[
            root > 0 && other > 0 && root < other,
            Assumptions -> ass && cSym > 0
          ]
        ]
      ]
    ]
  ];
  negativeBranch = SelectFirst[
    roots,
    TrueQ[FullSimplify[# > 0, Assumptions -> ass && cSym < 0]] &
  ];
  linearBranch = tau /. First[Solve[comparisonPolynomial[0, tau] == 0, tau]];
  If[MissingQ[positiveBranch], fail["M4 " <> label <> " c>0 root selection", roots]];
  If[MissingQ[negativeBranch], fail["M4 " <> label <> " c<0 root selection", roots]];
  Piecewise[{{positiveBranch, cSym > 0}, {negativeBranch, cSym < 0}}, linearBranch]
];

expectZero[
  "M4 slope convention coefficient for tau is -k",
  Coefficient[comparisonPolynomial[kappaLo, tau], tau] + k
];

tauLo = monotoneFirstRoot["lower envelope", kappaLo];
tauHi = monotoneFirstRoot["upper envelope", kappaHi];

expectZero["M4 tau_lo selected root - closed form", tauLo - closedMonotoneRoot[kappaLo]];
expectZero["M4 tau_lo quadratic residual", comparisonPolynomial[kappaLo, tauLo]];
expectZero["M4 tau_hi selected root - closed form", tauHi - closedMonotoneRoot[kappaHi]];
expectZero["M4 tau_hi quadratic residual", comparisonPolynomial[kappaHi, tauHi]];

subbanner["M5. Certified turning bracket"];
turningPolynomial[curvature_, x_] := H0 + 1/2 curvature x^2;

turningPositiveRoot[label_String, curvature_] := Module[
  {roots, root},
  roots = stripConditional[
    tau /. Solve[H0 - 1/2 (-curvature) tau^2 == 0, tau]
  ];
  root = SelectFirst[
    roots,
    TrueQ[FullSimplify[# > 0, Assumptions -> $Assumptions && curvature < 0]] &
  ];
  If[MissingQ[root], fail["M5 " <> label <> " positive turning root selection", roots]];
  root
];

tauTurnLo = turningPositiveRoot["lower", kappaTurnLo];
tauTurnHi = turningPositiveRoot["upper", kappaTurnHi];

expectZero["M5 tau_lo^(tp) - Sqrt[-2H0/kappa_lo]", tauTurnLo - Sqrt[-2 H0/kappaTurnLo]];
expectZero["M5 tau_lo^(tp) turning residual", turningPolynomial[kappaTurnLo, tauTurnLo]];
expectZero["M5 tau_hi^(tp) - Sqrt[-2H0/kappa_hi]", tauTurnHi - Sqrt[-2 H0/kappaTurnHi]];
expectZero["M5 tau_hi^(tp) turning residual", turningPolynomial[kappaTurnHi, tauTurnHi]];

subbanner["M6. Sign-adapted primitive drift table"];
aStar = (1 + deltaUStar)/(1 + chi0Star);

primitiveDrift[dir_List] := Module[
  {sd, st, sk, sm},
  sd = -aStar (dir[[3]] + dir[[2]] - dir[[4]]);
  st = dir[[4]] + sd;
  sk = 2 dir[[2]] - dir[[4]];
  sm = 2 dir[[2]] - dir[[4]] + 2 dir[[5]] - 2 dir[[1]] -
    eStar (2 dir[[3]] + 2 dir[[1]] - dir[[4]] - dir[[5]]) +
    fStar sd;
  {sd, st, sk, sm}
];

epsList = {epsLambda, epsC, epsGamma, epsU, epsW};
primitiveRows = {
  {0, 0, 0, -2 - 2 eStar},
  {-aStar, -aStar, 2, 2 - fStar aStar},
  {-aStar, -aStar, 0, -2 eStar - fStar aStar},
  {aStar, 1 + aStar, -1, -1 + eStar + fStar aStar},
  {0, 0, 0, 2 + eStar}
};

Do[
  Module[{dir = epsList[[i]] UnitVector[5, i], expected = epsList[[i]] primitiveRows[[i]]},
    expectZero[
      "M6 primitive exponent row (" <> labels[[i]] <> ")",
      primitiveDrift[dir] - expected
    ];
  ],
  {i, 1, 5}
];

banner["STAGE 207 MATHEMATICA AUDIT PASSED"];
Exit[0];
