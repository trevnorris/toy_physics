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

expectZero[name_String, expr_] := Module[{res},
  res = cleanScalar[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, expr_] := Module[{res},
  res = FullSimplify[expr, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

selectByLimit[name_String, roots_List, var_Symbol, target_, assumptions_] := Module[
  {choice},
  choice = SelectFirst[
    roots,
    Function[root,
      TrueQ[
        FullSimplify[
          stripConditional[
            FullSimplify[
              Limit[root, var -> 0, Assumptions -> assumptions] - target,
              Assumptions -> assumptions
            ]
          ] == 0,
          Assumptions -> assumptions
        ]
      ]
    ],
    Missing["NoBranch"]
  ];
  If[MissingQ[choice], fail[name <> " branch selection", roots]];
  FullSimplify[choice, Assumptions -> assumptions]
];

banner["STAGE 205 -- DIRECTIONAL HESSIAN AND QUADRATIC ROOT REFINEMENT"];

Clear[
  h11, h12, h13, h14, h15, h22, h23, h24, h25, h33, h34, h35,
  h44, h45, h55, g1, g2, g3, g4, g5, d1, d2, d3, d4, d5,
  x1, x2, x3, x4, x5, u, tau, phiBase,
  p0ap, p1ap, p2ap, p0an, p1anAbs, p2an,
  p0lp, l0lp, l1lp, p0ln, l0lnAbs, l1ln,
  phi0turn, phi2turn, tauRoot,
  phi0tan, phi1tan, phi2tan,
  eps, l0s, l1s
];

$Assumptions = (
  Element[
    {
      h11, h12, h13, h14, h15, h22, h23, h24, h25, h33, h34, h35,
      h44, h45, h55, g1, g2, g3, g4, g5, d1, d2, d3, d4, d5,
      x1, x2, x3, x4, x5, u, tau, phiBase,
      p0ap, p1ap, p2ap, p0an, p1anAbs, p2an,
      p0lp, l0lp, l1lp, p0ln, l0lnAbs, l1ln,
      phi0turn, phi2turn, tauRoot,
      phi0tan, phi1tan, phi2tan,
      eps, l0s, l1s
    },
    Reals
  ] &&
  phiBase > 0 &&
  p0ap > 0 && p1ap > 0 &&
  p0an > 0 && p1anAbs > 0 &&
  p0lp > 0 && l0lp > 0 &&
  p0ln > 0 && l0lnAbs > 0 &&
  l0s > 0
);
baseAssumptions = $Assumptions;

ordinaryModel[p0_, p1_, p2_, z_] := (p0 - 1) + p1*z + 1/2*p2*z^2;
logModel[p0_, ell0_, ell1_, z_] := Log[p0] + ell0*z + 1/2*ell1*z^2;

ordinaryBranch[p0_, p1_, p2_, assumptions_] := Module[{roots, target},
  roots = Block[
    {$Assumptions = True},
    tau /. Solve[ordinaryModel[p0, p1, p2, tau] == 0, tau]
  ];
  target = (1 - p0)/p1;
  selectByLimit["ordinary quadratic", roots, p2, target, assumptions]
];

logBranch[p0_, ell0_, ell1_, assumptions_] := Module[{roots, target},
  roots = Block[
    {$Assumptions = True},
    tau /. Solve[logModel[p0, ell0, ell1, tau] == 0, tau]
  ];
  target = -Log[p0]/ell0;
  selectByLimit["logarithmic quadratic", roots, ell1, target, assumptions]
];

nearZeroBranch[name_String, roots_List, assumptions_] :=
  selectByLimit[name, roots, eps, 0, assumptions];

subbanner["M1-M2. Log-Hessian identities from second derivatives"];
vars = {x1, x2, x3, x4, x5};
gradVec = {g1, g2, g3, g4, g5};
rayDir = {d1, d2, d3, d4, d5};
hessMat = {
  {h11, h12, h13, h14, h15},
  {h12, h22, h23, h24, h25},
  {h13, h23, h33, h34, h35},
  {h14, h24, h34, h44, h45},
  {h15, h25, h35, h45, h55}
};

chiLocal = (
  phiBase + gradVec . vars + 1/2*vars . hessMat . vars
);
chiAlong = chiLocal /. Thread[vars -> u*rayDir];
ordinarySlope = D[chiAlong, u] /. u -> 0;
ordinaryCurve = D[chiAlong, {u, 2}] /. u -> 0;
logHessianByD = Table[
  D[Log[chiLocal], vars[[i]], vars[[j]]] /. Thread[vars -> ConstantArray[0, 5]],
  {i, 1, 5}, {j, 1, 5}
];
logSlope = ordinarySlope/phiBase;
logCurve = rayDir . logHessianByD . rayDir;

expectZero[
  "M1 log-curvature identity",
  logCurve - (ordinaryCurve/phiBase - ordinarySlope^2/phiBase^2)
];
expectZero[
  "M2 ordinary/log bridge",
  ordinaryCurve - phiBase*(logCurve + logSlope^2)
];

subbanner["M3-M4. Affine quadratic predictor branches from Solve"];
tauAffPos = ordinaryBranch[
  p0ap,
  p1ap,
  p2ap,
  Element[{p0ap, p1ap}, Reals] && p0ap > 0 && p1ap > 0
];
expectZero[
  "M3 affine residual positive slope",
  ordinaryModel[p0ap, p1ap, p2ap, tauAffPos]
];
expectZero[
  "M4 affine zero-curvature limit positive slope",
  Limit[
    tauAffPos,
    p2ap -> 0,
    Assumptions -> Element[{p0ap, p1ap}, Reals] && p0ap > 0 && p1ap > 0
  ] - (1 - p0ap)/p1ap
];

tauAffNeg = ordinaryBranch[
  p0an,
  -p1anAbs,
  p2an,
  Element[{p0an, p1anAbs}, Reals] && p0an > 0 && p1anAbs > 0
];
expectZero[
  "M3 affine residual negative slope",
  ordinaryModel[p0an, -p1anAbs, p2an, tauAffNeg]
];
expectZero[
  "M4 affine zero-curvature limit negative slope",
  Limit[
    tauAffNeg,
    p2an -> 0,
    Assumptions -> Element[{p0an, p1anAbs}, Reals] && p0an > 0 && p1anAbs > 0
  ] - (1 - p0an)/(-p1anAbs)
];

subbanner["M5-M6. Logarithmic quadratic predictor branches from Solve"];
tauLogPos = logBranch[
  p0lp,
  l0lp,
  l1lp,
  Element[{p0lp, l0lp}, Reals] && p0lp > 0 && l0lp > 0
];
expectZero[
  "M5 log residual positive slope",
  logModel[p0lp, l0lp, l1lp, tauLogPos]
];
expectZero[
  "M6 log zero-curvature limit positive slope",
  Limit[
    tauLogPos,
    l1lp -> 0,
    Assumptions -> Element[{p0lp, l0lp}, Reals] && p0lp > 0 && l0lp > 0
  ] + Log[p0lp]/l0lp
];

tauLogNeg = logBranch[
  p0ln,
  -l0lnAbs,
  l1ln,
  Element[{p0ln, l0lnAbs}, Reals] && p0ln > 0 && l0lnAbs > 0
];
expectZero[
  "M5 log residual negative slope",
  logModel[p0ln, -l0lnAbs, l1ln, tauLogNeg]
];
expectZero[
  "M6 log zero-curvature limit negative slope",
  Limit[
    tauLogNeg,
    l1ln -> 0,
    Assumptions -> Element[{p0ln, l0lnAbs}, Reals] && p0ln > 0 && l0lnAbs > 0
  ] + Log[p0ln]/(-l0lnAbs)
];

subbanner["M7. Turning-point reality criterion"];
turnRadicand = 2*(1 - phi0turn)/phi2turn;
turnProduct = (1 - phi0turn)*phi2turn;
radicandRegion = Reduce[
  turnRadicand >= 0 && phi0turn != 1 && phi2turn != 0,
  {phi0turn, phi2turn},
  Reals
];
criterionRegion = Reduce[
  turnProduct > 0 && phi0turn != 1 && phi2turn != 0,
  {phi0turn, phi2turn},
  Reals
];
negativeRegion = Reduce[
  turnRadicand >= 0 && turnProduct < 0 && phi0turn != 1 && phi2turn != 0,
  {phi0turn, phi2turn},
  Reals
];
rootExistRegion = Reduce[
  Exists[
    tauRoot,
    Element[tauRoot, Reals] &&
      ordinaryModel[phi0turn, 0, phi2turn, tauRoot] == 0
  ] && phi0turn != 1 && phi2turn != 0,
  {phi0turn, phi2turn},
  Reals
];

expectTrue["M7 radicand criterion via Reduce", Equivalent[radicandRegion, criterionRegion]];
expectTrue["M7 no real root on negative product", negativeRegion === False];
expectTrue["M7 real-root existence region", Equivalent[rootExistRegion, criterionRegion]];

Block[
  {$Assumptions = baseAssumptions && turnProduct > 0 && phi0turn != 1 && phi2turn != 0},
  expectZero[
    "M7 tau-plus residual on criterion region",
    ordinaryModel[phi0turn, 0, phi2turn, Sqrt[turnRadicand]]
  ];
  expectZero[
    "M7 tau-minus residual on criterion region",
    ordinaryModel[phi0turn, 0, phi2turn, -Sqrt[turnRadicand]]
  ];
];

subbanner["M8. Tangency model"];
tangentResidual = (
  ordinaryModel[phi0tan, phi1tan, phi2tan, tau] /. {phi0tan -> 1, phi1tan -> 0}
) - 1/2*phi2tan*tau^2;
expectZero["M8 tangency model from quadratic closure", tangentResidual];

subbanner["M9-M11. Local curvature expansions"];
p0Series = 1 + eps;
p1Series = p0Series*l0s;
p2Series = p0Series*(l1s + l0s^2);
quadNearRoots = Block[
  {$Assumptions = True},
  tau /. Solve[ordinaryModel[p0Series, p1Series, p2Series, tau] == 0, tau]
];
logNearRoots = Block[
  {$Assumptions = True},
  tau /. Solve[logModel[p0Series, l0s, l1s, tau] == 0, tau]
];
tauQuadNear = nearZeroBranch[
  "near ordinary predictor",
  quadNearRoots,
  Element[{l0s, l1s}, Reals] && l0s > 0
];
tauLog2Near = nearZeroBranch[
  "near log predictor",
  logNearRoots,
  Element[{l0s, l1s}, Reals] && l0s > 0
];
tauAffNear = (1 - p0Series)/p1Series;
tauLogNear = -Log[p0Series]/l0s;

ordinaryCorrection = Normal[
  Series[
    tauQuadNear - tauAffNear + (p2Series/(2*p1Series))*tauAffNear^2,
    {eps, 0, 2}
  ]
];
expectZero["M9 ordinary curvature correction", ordinaryCorrection];

logCorrection = Normal[
  Series[
    tauLog2Near - tauLogNear + (l1s/(2*l0s))*tauLogNear^2,
    {eps, 0, 2}
  ]
];
expectZero["M10 logarithmic curvature correction", logCorrection];

gapSeries = Normal[Series[tauLog2Near - tauQuadNear, {eps, 0, 3}]];
expectZero["M11 gap coefficient eps^0", Coefficient[gapSeries, eps, 0]];
expectZero["M11 gap coefficient eps^1", Coefficient[gapSeries, eps, 1]];
expectZero["M11 gap coefficient eps^2", Coefficient[gapSeries, eps, 2]];
expectZero[
  "M11 gap coefficient eps^3",
  Coefficient[gapSeries, eps, 3] - (l0s^2 + 3*l1s)/(6*l0s^3)
];

banner["STAGE 205 MATHEMATICA AUDIT PASSED"];
Exit[0];
