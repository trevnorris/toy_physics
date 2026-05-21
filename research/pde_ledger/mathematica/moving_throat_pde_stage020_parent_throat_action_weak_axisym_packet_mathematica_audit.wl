(* Unit 020 Mathematica audit: parent throat action weak-axisymmetric packet. *)

ClearAll["Global`*"];
$HistoryLength = 0;

ClearAll[GauntIntegral];
SetAttributes[GauntIntegral, Listable];

Module[
  {
    show, gauntWeight, overlapBase, lambda0, lambda1, lambda2,
    crossOne, crossTwo, angularResiduals,
    KSigma, B0, Z0, N0, N01, dKSigma, dMSigma,
    B01, B21, B41, Z01, Z21, Z41,
    D0, D01, D21, D41, K1, HEven, Xi1,
    gateJacobian, branchSolve, branchRules, expectedRules,
    determinantResidual, solveResiduals, deficitResiduals, xiResidual
  },

  show[expr_] := ToString[InputForm[expr]];
  gauntWeight[la_, lb_, lc_] := Sqrt[(2 la + 1) (2 lb + 1) (2 lc + 1)/(4 Pi)];
  GauntIntegral[la_Integer, lb_Integer, lc_Integer, ma_Integer, mb_Integer, mc_Integer] :=
    Quiet[
      FullSimplify[
        gauntWeight[la, lb, lc]*
          ThreeJSymbol[{la, 0}, {lb, 0}, {lc, 0}]*
          ThreeJSymbol[{la, ma}, {lb, mb}, {lc, mc}]
      ],
      ClebschGordan::phy
    ];

  Print["STAGE 020 PARENT THROAT ACTION WEAK-AXISYMMETRIC PACKET MATHEMATICA AUDIT"];

  overlapBase = GauntIntegral[2, 2, 2, 0, 0, 0];
  lambda0 = FullSimplify[(-1)^0 GauntIntegral[2, 2, 2, 0, 0, 0]/overlapBase];
  lambda1 = FullSimplify[(-1)^1 GauntIntegral[2, 2, 2, 0, 1, -1]/overlapBase];
  lambda2 = FullSimplify[(-1)^2 GauntIntegral[2, 2, 2, 0, 2, -2]/overlapBase];
  crossOne = GauntIntegral[2, 2, 2, 0, 1, 1];
  crossTwo = GauntIntegral[2, 2, 2, 0, 2, 2];
  angularResiduals = {
    FullSimplify[lambda0 - 1],
    FullSimplify[lambda1 - 1/2],
    FullSimplify[lambda2 + 1],
    FullSimplify[crossOne],
    FullSimplify[crossTwo]
  };

  Print["M1 lambda_0 residual = ", show[angularResiduals[[1]]]];
  If[!TrueQ[FullSimplify[lambda0 - 1] === 0],
    Print["FAIL: M1 lambda_0"]; Exit[1]
  ];
  Print["M1 lambda_1 residual = ", show[angularResiduals[[2]]]];
  If[!TrueQ[FullSimplify[lambda1 - 1/2] === 0],
    Print["FAIL: M1 lambda_1"]; Exit[1]
  ];
  Print["M1 lambda_2 residual = ", show[angularResiduals[[3]]]];
  If[!TrueQ[FullSimplify[lambda2 - (-1)] === 0],
    Print["FAIL: M1 lambda_2"]; Exit[1]
  ];
  Print["M1 same-sign m=1 residual = ", show[angularResiduals[[4]]]];
  If[!TrueQ[FullSimplify[crossOne - 0] === 0],
    Print["FAIL: M1 same-sign m=1"]; Exit[1]
  ];
  Print["M1 same-sign m=2 residual = ", show[angularResiduals[[5]]]];
  If[!TrueQ[FullSimplify[crossTwo - 0] === 0],
    Print["FAIL: M1 same-sign m=2"]; Exit[1]
  ];
  Print["M1 OK"];

  Clear[KSigma, B0, Z0, N0, N01, dKSigma, dMSigma, B01, B21, B41, Z01, Z21, Z41];
  $Assumptions =
    Element[
      {KSigma, B0, Z0, N0, N01, dKSigma, dMSigma, B01, B21, B41, Z01, Z21, Z41},
      Reals
    ] && N0 != 0 && KSigma - B0 - Z0 != 0;

  D0 = KSigma - B0 - Z0;
  D01 = dKSigma - B01 - Z01;
  D21 = -(dMSigma + B21 + Z21);
  D41 = -(B41 + Z41);
  K1 = D21 + D01/9;
  HEven = D41 - (2/3) D21 - D01/27;
  Xi1 = N01/N0 - D01/D0;

  gateJacobian = {
    {D[K1, dKSigma], D[K1, dMSigma]},
    {D[HEven, dKSigma], D[HEven, dMSigma]}
  };
  determinantResidual = FullSimplify[Det[gateJacobian] - 1/27];
  Print["M2 even-gate Jacobian determinant residual = ", show[determinantResidual]];
  If[!TrueQ[FullSimplify[Det[gateJacobian] - 1/27] === 0],
    Print["FAIL: M2 even-gate Jacobian determinant"]; Exit[1]
  ];
  Print["M2 OK"];

  branchSolve = Solve[{K1 == 0, HEven == 0}, {dKSigma, dMSigma}];
  Print["M3 solution count residual = ", show[Length[branchSolve] - 1]];
  If[!TrueQ[FullSimplify[Length[branchSolve] - 1] === 0],
    Print["FAIL: M3 unique solution"]; Exit[1]
  ];
  branchRules = First[branchSolve];
  expectedRules = {
    dKSigma -> B01 + Z01 + 27 (B41 + Z41),
    dMSigma -> -(B21 + Z21) + 3 (B41 + Z41)
  };
  solveResiduals = FullSimplify[
    ({dKSigma, dMSigma} /. branchRules) - ({dKSigma, dMSigma} /. expectedRules)
  ];
  Print["M3 dKSigma residual = ", show[solveResiduals[[1]]]];
  If[!TrueQ[FullSimplify[(dKSigma /. branchRules) - (B01 + Z01 + 27 (B41 + Z41))] === 0],
    Print["FAIL: M3 dKSigma"]; Exit[1]
  ];
  Print["M3 dMSigma residual = ", show[solveResiduals[[2]]]];
  If[!TrueQ[FullSimplify[(dMSigma /. branchRules) - (-(B21 + Z21) + 3 (B41 + Z41))] === 0],
    Print["FAIL: M3 dMSigma"]; Exit[1]
  ];
  Print["M3 OK"];

  deficitResiduals = FullSimplify[
    {
      (D01 /. branchRules) - 27 (B41 + Z41),
      (D21 /. branchRules) - (-3 (B41 + Z41)),
      (D41 /. branchRules) - (-(B41 + Z41))
    }
  ];
  Print["M4 D01 residual = ", show[deficitResiduals[[1]]]];
  If[!TrueQ[FullSimplify[(D01 /. branchRules) - 27 (B41 + Z41)] === 0],
    Print["FAIL: M4 D01"]; Exit[1]
  ];
  Print["M4 D21 residual = ", show[deficitResiduals[[2]]]];
  If[!TrueQ[FullSimplify[(D21 /. branchRules) - (-3 (B41 + Z41))] === 0],
    Print["FAIL: M4 D21"]; Exit[1]
  ];
  Print["M4 D41 residual = ", show[deficitResiduals[[3]]]];
  If[!TrueQ[FullSimplify[(D41 /. branchRules) - (-(B41 + Z41))] === 0],
    Print["FAIL: M4 D41"]; Exit[1]
  ];
  Print["M4 OK"];

  xiResidual =
    FullSimplify[
      (Xi1 /. branchRules) -
        (N01/N0 - 27 (B41 + Z41)/(KSigma - B0 - Z0))
    ];
  Print["M5 Xi1 residual = ", show[xiResidual]];
  If[!TrueQ[FullSimplify[(Xi1 /. branchRules) -
      (N01/N0 - 27 (B41 + Z41)/(KSigma - B0 - Z0))] === 0],
    Print["FAIL: M5 Xi1"]; Exit[1]
  ];
  Print["M5 OK"];

  Print["STATUS: PASS"];
  Exit[0]
];
