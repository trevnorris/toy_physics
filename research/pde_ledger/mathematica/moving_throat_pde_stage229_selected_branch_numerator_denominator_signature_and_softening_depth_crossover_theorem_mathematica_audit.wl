ClearAll["Global`*"];
$HistoryLength = 0;

Print["Stage 229 Mathematica audit: selected-branch numerator/denominator signature"];

fmt[expr_] := ToString[InputForm[expr]];

stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL ", name];
  If[!MissingQ[detail], Print["  detail = ", fmt[detail]]];
  Exit[1];
);

pass[name_String] := Print["PASS ", name];

cleanScalar[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Cancel[Together[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

expectZero[expr_] := expectZero["unnamed zero check", expr];

expectZero[name_String, expr_] := Module[{res},
  res = cleanScalar[expr];
  Print[name, " residual = ", fmt[res]];
  If[Simplify[res, Assumptions -> $Assumptions] =!= 0, fail[name, res], pass[name]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = stripConditional[FullSimplify[cond, Assumptions -> $Assumptions]];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

expectNumeric[name_String, actual_, expected_, tol_] := Module[{err},
  err = N[Abs[actual - expected], 30];
  Print[name, " error = ", fmt[err]];
  If[TrueQ[err < tol], pass[name], fail[name, err]];
];

stableRoot[poly_, deltaValue_] := Module[{sols, roots},
  sols = xi /. NSolve[(poly /. delta -> deltaValue) == 0, xi, Reals];
  roots = Select[N[sols, 30], TrueQ[0 < # < 1] &];
  Print["M10 roots for delta = ", fmt[deltaValue], " -> ", fmt[roots]];
  If[Length[roots] == 1, First[roots], fail["M10 unique root in (0,1)", roots]]
];

domainAssumptions = (
  Element[{xi, delta, A, beta0, x, DeltaKax}, Reals]
    && xi > 0 && xi < 1 && delta > 0 && A > 0 && beta0 > 0
);
$Assumptions = domainAssumptions;

Assuming[
  xi > 0 && xi < 1 && delta > 0 && A > 0 && beta0 > 0,
  (
    kappa0sq = 8/Pi^2;
    kappa1sq = 16/(9 Pi^2);

    sMinus = (
      (kappa0sq*(x + DeltaKax) + kappa1sq*x)^2
        /(kappa0sq*(x + DeltaKax)^2 + kappa1sq*x^2)
    );
    nMinus = (
      beta0*sMinus^2/(kappa0sq*(A - x))
    );

    selectedN = Cancel[Together[nMinus /. {x -> A*xi, DeltaKax -> A*delta}]];
    scaleN = 8*beta0/(Pi^2*A);
    F = Cancel[Together[selectedN/scaleN]];
    FClaim = (
      (9*delta + 11*xi)^4
        /(81*(1 - xi)*(9*delta^2 + 18*delta*xi + 11*xi^2)^2)
    );

    Print["M1 sMinus = ", fmt[Factor[sMinus]]];
    Print["M1 derived F = ", fmt[Factor[F]]];
    expectZero["M1 selected N_- reduction", selectedN - scaleN*FClaim];
    expectZero["M1 dimensionless F from constants", F - FClaim];

    FNum = (
      (9*delta + 11*xi)^4
        /(81*(9*delta^2 + 18*delta*xi + 11*xi^2)^2)
    );
    FDen = 1/(1 - xi);
    expectZero["M2 F factorization", F - FNum*FDen];

    LNum = FullSimplify[D[Log[FNum], xi], Assumptions -> $Assumptions];
    LDen = FullSimplify[D[Log[FDen], xi], Assumptions -> $Assumptions];
    LNumClaim = (
      72*delta^2/((9*delta + 11*xi)*(9*delta^2 + 18*delta*xi + 11*xi^2))
    );
    LDenClaim = 1/(1 - xi);
    Print["M3 L_num = ", fmt[Factor[LNum]]];
    Print["M3 L_den = ", fmt[Factor[LDen]]];
    expectZero["M3 numerator log derivative", LNum - LNumClaim];
    expectZero["M3 denominator log derivative", LDen - LDenClaim];

    RND = Cancel[Together[LNum/LDen]];
    RNDClaim = (
      72*delta^2*(1 - xi)
        /((9*delta + 11*xi)*(9*delta^2 + 18*delta*xi + 11*xi^2))
    );
    Print["M4 R_ND = ", fmt[Factor[RND]]];
    expectZero["M4 classifier reduction", RND - RNDClaim];

    onset = FullSimplify[RND /. xi -> 0, Assumptions -> delta > 0];
    expectZero["M5 onset value", onset - 8/(9*delta)];

    RNDSoft = Block[
      {$Assumptions = delta > 0},
      Limit[RND, xi -> 1, Direction -> "FromBelow"]
    ];
    LNumSoft = FullSimplify[
      Block[
        {$Assumptions = delta > 0},
        Limit[LNum, xi -> 1, Direction -> "FromBelow"]
      ],
      Assumptions -> delta > 0
    ];
    LDenSoft = Block[
      {$Assumptions = True},
      Limit[LDen, xi -> 1, Direction -> "FromBelow"]
    ];
    LNumSoftClaim = 72*delta^2/((9*delta + 11)*(9*delta^2 + 18*delta + 11));
    expectZero["M6 R_ND softening limit", RNDSoft];
    expectZero["M6 L_num softening limit", LNumSoft - LNumSoftClaim];
    expectTrue[
      "M6 L_den pole reciprocal",
      Simplify[1/LDenSoft == 0, Assumptions -> delta > 0]
    ];

    rawP = Expand[Numerator[Together[RND - 1]]];
    crossoverDen = Factor[Denominator[Together[RND - 1]]];
    P = Expand[If[TrueQ[Coefficient[rawP, xi, 3] < 0], -rawP, rawP]];
    pCoeffs = CoefficientList[P, xi];
    targetCoeffs = {
      81*delta^3 - 72*delta^2,
      333*delta^2,
      297*delta,
      121
    };
    Print["M7 raw numerator = ", fmt[Factor[rawP]]];
    Print["M7 crossover denominator = ", fmt[crossoverDen]];
    Print["M7 sign-normalized P = ", fmt[Factor[P]]];
    Print["M7 coefficient list = ", fmt[pCoeffs]];
    Do[
      expectZero[
        "M7 coefficient " <> ToString[j - 1],
        pCoeffs[[j]] - targetCoeffs[[j]]
      ],
      {j, 1, Length[targetCoeffs]}
    ];

    dP = Expand[D[P, xi]];
    dPCoeffs = CoefficientList[dP, xi];
    targetDPCoeffs = {333*delta^2, 594*delta, 363};
    Print["M8 dP/dxi = ", fmt[Factor[dP]]];
    Do[
      expectZero[
        "M8 derivative coefficient " <> ToString[j - 1],
        dPCoeffs[[j]] - targetDPCoeffs[[j]]
      ],
      {j, 1, Length[targetDPCoeffs]}
    ];
    positiveD = Resolve[
      ForAll[
        {xiPos, deltaPos},
        Implies[
          xiPos >= 0 && deltaPos > 0,
          (dP /. {xi -> xiPos, delta -> deltaPos}) > 0
        ]
      ],
      Reals
    ];
    expectTrue["M8 derivative positivity", positiveD];

    P0 = FullSimplify[P /. xi -> 0, Assumptions -> delta > 0];
    threshold = delta /. First[Solve[9*delta - 8 == 0, delta, Reals]];
    Print["M9 P(0,delta) = ", fmt[Factor[P0]]];
    expectZero["M9 onset polynomial factor", P0 - 9*delta^2*(9*delta - 8)];
    expectZero["M9 threshold root", threshold - 8/9];

    sampleData = {
      {1/4, 0.107223051105697},
      {1/2, 0.081847937860074},
      {3/4, 0.032505121082825}
    };
    Do[
      Module[{delta0, expected, root, leftProbe, rightProbe, leftVal, rightVal},
        {delta0, expected} = sample;
        root = stableRoot[P, delta0];
        expectNumeric[
          "M10 crossover root delta = " <> fmt[delta0],
          root,
          expected,
          5*10^-13
        ];
        leftProbe = Max[10^-6, root/2];
        rightProbe = Min[1 - 10^-6, (root + 1)/2];
        leftVal = N[RND /. {xi -> leftProbe, delta -> delta0}, 30];
        rightVal = N[RND /. {xi -> rightProbe, delta -> delta0}, 30];
        Print[
          "M10 side values delta = ", fmt[delta0],
          ": left ", fmt[leftVal], ", right ", fmt[rightVal]
        ];
        If[!TrueQ[leftVal > 1], fail["M10 left-of-root numerator-like", leftVal]];
        If[!TrueQ[rightVal < 1], fail["M10 right-of-root denominator-like", rightVal]];
        pass["M10 side classification delta = " <> fmt[delta0]];
      ],
      {sample, sampleData}
    ];

    denominatorProbes = {1/100, 1/5, 3/5, 9/10};
    denominatorValues = N[RND /. delta -> 1 /. xi -> #, 30] & /@ denominatorProbes;
    Print["M10 delta=1 probe values = ", fmt[denominatorValues]];
    expectTrue["M10 delta=1 always denominator slice", AllTrue[denominatorValues, # < 1 &]];

    Print["Stage 229 Mathematica audit passed."];
  )
];

Exit[0];
