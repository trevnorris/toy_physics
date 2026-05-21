(* moving_throat_pde_stage017_parent_throat_action_weak_axisym_mathematica_audit.wl *)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];
reduce[expr_] := FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];

checkZero[label_String, expr_] := Module[{res},
  res = reduce[expr];
  Print[label, " residual = ", fmt[res]];
  If[TrueQ[res == 0],
    Print["PASS: ", label],
    Print["FAIL: ", label]; Exit[1]
  ];
];

Print["STAGE 017 PARENT THROAT ACTION WEAK AXISYM MATHEMATICA AUDIT"];

Clear[theta, phi];
$Assumptions = Element[{theta, phi}, Reals];

tripleOnSphere[l1_Integer, l2_Integer, l3_Integer, q1_Integer, q2_Integer, q3_Integer] :=
  tripleOnSphere[l1, l2, l3, q1, q2, q3] = Module[{expanded},
    expanded =
      FullSimplify[
        ComplexExpand[
          FunctionExpand[
            SphericalHarmonicY[l1, q1, theta, phi]*
              SphericalHarmonicY[l2, q2, theta, phi]*
              SphericalHarmonicY[l3, q3, theta, phi]
          ]
        ],
        Assumptions -> 0 <= theta <= Pi && 0 <= phi <= 2*Pi
      ];
    FullSimplify[
      Integrate[
        Sin[theta]*expanded,
        {phi, 0, 2*Pi},
        {theta, 0, Pi},
        GenerateConditions -> False
      ]
    ]
  ];

baseTriple = tripleOnSphere[2, 2, 2, 0, 0, 0];
laneFactor[m_Integer] :=
  laneFactor[m] =
    reduce[(-1)^m*tripleOnSphere[2, 2, 2, 0, m, -m]/baseTriple];

checkZero["M1 lane ratio m=0", laneFactor[0] - 1];
checkZero["M2 lane ratio m=1", laneFactor[1] - 1/2];
checkZero["M3 lane ratio m=2", laneFactor[2] + 1];
checkZero["M4 same-sign cross term m=1", tripleOnSphere[2, 2, 2, 0, 1, 1]];
checkZero["M4 same-sign cross term m=2", tripleOnSphere[2, 2, 2, 0, 2, 2]];

Clear[dMsym, dKsym, M1, K1w, eps, D0, N0];
$Assumptions =
  Element[{dMsym, dKsym, M1, K1w, eps, D0, N0}, Reals] && D0 != 0;

kWall = -dMsym + dKsym/9;
hWall = 2*dMsym/3 - dKsym/27;
wallJacobian = {{D[kWall, dKsym], D[kWall, dMsym]}, {D[hWall, dKsym], D[hWall, dMsym]}};
checkZero["M9 wall Jacobian determinant", Det[wallJacobian] - 1/27];
wallSolve = Solve[{kWall == 0, hWall == 0}, {dKsym, dMsym}];
wallSolutionVector = If[Length[wallSolve] == 1, {dKsym, dMsym} /. First[wallSolve], {1, 1}];
checkZero["M10 wall-only solution count", Length[wallSolve] - 1];
checkZero["M10 wall-only dK solution", wallSolutionVector[[1]]];
checkZero["M10 wall-only dM solution", wallSolutionVector[[2]]];

massLane[m_Integer] := eps*laneFactor[m]*M1;
stiffLane[m_Integer] := eps*laneFactor[m]*K1w;
d01Lane[m_Integer] := stiffLane[m];
d21Lane[m_Integer] := -massLane[m];
d41Lane[m_Integer] := 0;
kGate[m_Integer] := d21Lane[m] + d01Lane[m]/9;
hGate[m_Integer] := d41Lane[m] - (2/3)*d21Lane[m] - d01Lane[m]/27;

Do[
  checkZero[
    "M7 K1 gate lane m=" <> ToString[m],
    kGate[m] - eps*laneFactor[m]*(-M1 + K1w/9)
  ],
  {m, {0, 1, 2}}
];
Do[
  checkZero[
    "M8 H_even gate lane m=" <> ToString[m],
    hGate[m] - eps*laneFactor[m]*(2*M1/3 - K1w/27)
  ],
  {m, {0, 1, 2}}
];

weightedMean[vals_List] := reduce[(vals[[1]] + 2*vals[[2]] + 2*vals[[3]])/5];
defectA[vals_List] := reduce[(2*vals[[1]] - vals[[2]] - vals[[3]])/10];
defectB[vals_List] := reduce[(vals[[2]] - vals[[3]])/2];

massValues = {massLane[0], massLane[1], massLane[2]};
stiffValues = {stiffLane[0], stiffLane[1], stiffLane[2]};
xiValues = {-d01Lane[0]/D0, -d01Lane[1]/D0, -d01Lane[2]/D0};
prefactorValues = {
  -N0*d01Lane[0]/D0^2,
  -N0*d01Lane[1]/D0^2,
  -N0*d01Lane[2]/D0^2
};

checkZero["M5 wall inertia grouped trace", weightedMean[massValues]];
checkZero["M5 wall inertia grouped line", defectB[massValues] - 3*defectA[massValues]];
checkZero["M6 wall stiffness grouped trace", weightedMean[stiffValues]];
checkZero["M6 wall stiffness grouped line", defectB[stiffValues] - 3*defectA[stiffValues]];
checkZero["M11 Xi load grouped trace", weightedMean[xiValues]];
checkZero["M11 Xi load grouped line", defectB[xiValues] - 3*defectA[xiValues]];
checkZero["M12 prefactor grouped trace", weightedMean[prefactorValues]];
checkZero["M12 prefactor grouped line", defectB[prefactorValues] - 3*defectA[prefactorValues]];

Print["STATUS: PASS"];
Exit[0];
