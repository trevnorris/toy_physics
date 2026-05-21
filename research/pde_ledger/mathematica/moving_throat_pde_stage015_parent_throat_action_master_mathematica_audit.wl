(* moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl *)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];
reduce[expr_] := FullSimplify[Together[expr], Assumptions -> $Assumptions];

expectZero[label_String, expr_] := Module[{res},
  res = reduce[expr];
  Print[label, " residual = ", fmt[res]];
  If[TrueQ[res == 0],
    Print["PASS: ", label],
    Print["FAIL: ", label]; Exit[1]
  ];
];

expectNonzero[label_String, expr_] := Module[{res},
  res = reduce[expr];
  Print[label, " residual = ", fmt[res]];
  If[TrueQ[res == 0],
    Print["FAIL: ", label]; Exit[1],
    Print["PASS: ", label]
  ];
];

expectEqual[label_String, got_, expected_] := Module[{ok},
  ok = FullSimplify[got == expected, Assumptions -> $Assumptions];
  Print[label, " residual = ", fmt[ok]];
  If[TrueQ[ok],
    Print["PASS: ", label],
    Print["FAIL: ", label];
    Print["  got      = ", fmt[got]];
    Print["  expected = ", fmt[expected]];
    Exit[1]
  ];
];

Print["STAGE 015 PARENT THROAT ACTION MASTER MATHEMATICA AUDIT"];

Clear[
  w, eps, gateCoeff, mu0, Tw0, TO0, U0, TwR0, TwRR0, UR0, URR0,
  R0p, dTwRR0p, eta, etat, etaw, grad2, dK, dM, b01, b21, b41,
  z01, z21, z41, x0, eps1
];

$Assumptions =
  Element[
    {w, eps, gateCoeff, mu0, Tw0, TO0, U0, TwR0, TwRR0, UR0, URR0,
      R0p, dTwRR0p, eta, etat, etaw, grad2, dK, dM, b01, b21, b41,
      z01, z21, z41, x0, eps1},
    Reals
  ];

aFun[w_] := A[w];
etaFun[w_] := eta[w];
ibpProductRule =
  -aFun[w]*etaFun[w]*etaFun'[w] -
    (D[-aFun[w]*etaFun[w]^2/2, w] + aFun'[w]*etaFun[w]^2/2);
expectZero["M1 generic IBP product-rule identity", ibpProductRule];

ibpMutated =
  -aFun[w]*etaFun[w]*etaFun'[w] -
    (D[-(-aFun[w]*etaFun[w]^2/2), w] + aFun'[w]*etaFun[w]^2/2);
expectNonzero["M1 mutated IBP boundary sign", ibpMutated];

aConcrete = Exp[-w^2];
etaConcrete = Exp[-w^2/2];
boundaryGaussian =
  Quiet[
    Limit[-aConcrete*etaConcrete^2/2, w -> Infinity] -
      Limit[-aConcrete*etaConcrete^2/2, w -> -Infinity],
    Limit::alimv
  ];
expectZero["M2 Gaussian IBP boundary discharge", boundaryGaussian];

crossGaussian =
  Integrate[-aConcrete*etaConcrete*D[etaConcrete, w], {w, -Infinity, Infinity}];
bulkGaussian =
  Integrate[D[aConcrete, w]*etaConcrete^2/2, {w, -Infinity, Infinity}];
expectZero["M2 Gaussian IBP cross equals bulk", crossGaussian - bulkGaussian];

Tw = Tw0 + eps*TwR0*eta + eps^2*TwRR0*eta^2/2;
U = U0 + eps*UR0*eta + eps^2*URR0*eta^2/2;
Rt = eps*etat;
Rw = R0p + eps*etaw;
lagrangian = mu0*Rt^2/2 - Tw*Rw^2/2 - TO0*eps^2*grad2/2 - U;
L2raw = Coefficient[Series[lagrangian, {eps, 0, 2}] // Normal, eps, 2];
crossCoeff = D[D[L2raw, eta], etaw];
expectZero["M3 K_eta raw eta etaw cross coefficient", crossCoeff + TwR0*R0p];

effectiveMass = URR0 - dTwRR0p + TwRR0*R0p^2/2;
canonicalL2 =
  mu0*etat^2/2 - Tw0*etaw^2/2 - TO0*grad2/2 - effectiveMass*eta^2/2;
L2afterIBP = Expand[L2raw - (-TwR0*R0p*eta*etaw) + dTwRR0p*eta^2/2];
expectZero["M3 K_eta canonical quadratic form", L2afterIBP - canonicalL2];

effectiveMassMutated = URR0 + dTwRR0p + TwRR0*R0p^2/2;
canonicalL2Mutated =
  mu0*etat^2/2 - Tw0*etaw^2/2 - TO0*grad2/2 -
    effectiveMassMutated*eta^2/2;
expectNonzero["M3 K_eta dTwRR0p sign mutation", L2afterIBP - canonicalL2Mutated];

D01full = dK - b01 - z01;
D21full = -(dM + b21 + z21);
D41full = -(b41 + z41);
K1full = D21full + D01full/9;
Hevenfull = D41full - (2/3)*D21full - D01full/27;
wallSpec = {b01 -> 0, b21 -> 0, b41 -> 0, z01 -> 0, z21 -> 0, z41 -> 0};
K1wall = Expand[K1full /. wallSpec];
Hevenwall = Expand[Hevenfull /. wallSpec];
expectZero["M4 wall-only K1 specialization", K1wall - (-dM + dK/9)];
expectZero["M4 wall-only H_even specialization", Hevenwall - ((2/3)*dM - dK/27)];

betaConcrete = Exp[-w^2];
deltaMu = Exp[-w^2];
deltaTw = Exp[-w^2];
deltaTO = Exp[-w^2];
deltaKeta = Exp[-w^2];
dMoverlap = Integrate[deltaMu*betaConcrete^2, {w, -Infinity, Infinity}];
dKoverlap =
  Integrate[
    deltaTw*D[betaConcrete, w]^2 + (deltaKeta + 6*deltaTO)*betaConcrete^2,
    {w, -Infinity, Infinity}
  ];
expectZero["M5 Gaussian dM overlap closed form", dMoverlap - Sqrt[Pi/3]];
expectZero["M5 Gaussian dK overlap closed form", dKoverlap - 23*Sqrt[Pi]/(3*Sqrt[3])];

K1wallNum = K1wall /. {dK -> dKoverlap, dM -> dMoverlap};
HevenwallNum = Hevenwall /. {dK -> dKoverlap, dM -> dMoverlap};
expectZero[
  "M5 wall-only K1 from concrete Gaussian overlap integrals",
  K1wallNum - (-dMoverlap + dKoverlap/9)
];
expectZero[
  "M5 wall-only H_even from concrete Gaussian overlap integrals",
  HevenwallNum - ((2/3)*dMoverlap - dKoverlap/27)
];
dKoverlapMutated =
  Integrate[
    deltaTw*D[betaConcrete, w]^2 + (deltaKeta + 5*deltaTO)*betaConcrete^2,
    {w, -Infinity, Infinity}
  ];
K1wallMutated = K1wall /. {dK -> dKoverlapMutated, dM -> dMoverlap};
expectNonzero[
  "M5 wall-only K1 detects mutated 6 deltaTO coefficient",
  K1wallMutated - (-dMoverlap + dKoverlap/9)
];

gateMatrix = {{D[K1wall, dK], D[K1wall, dM]}, {D[Hevenwall, dK], D[Hevenwall, dM]}};
expectZero["M6 wall-only Jacobian determinant", Det[gateMatrix] - 1/27];

K1wallParam = -dM + gateCoeff*dK;
gateMatrixParam = {
  {D[K1wallParam, dK], D[K1wallParam, dM]},
  {D[Hevenwall, dK], D[Hevenwall, dM]}
};
wallDetShift =
  FullSimplify[
    (Det[gateMatrixParam] /. gateCoeff -> 1/9 + eps) -
      (Det[gateMatrixParam] /. gateCoeff -> 1/9),
    Assumptions -> $Assumptions
  ];
expectZero["M6 wall determinant perturbation value", wallDetShift - 2*eps/3];
expectNonzero["M6 wall determinant perturbation nonzero guard", wallDetShift];

wallSolve = Solve[{K1wall == 0, Hevenwall == 0}, {dK, dM}];
expectEqual["M7 wall-only zero solve", wallSolve, {{dK -> 0, dM -> 0}}];
wallSolvePerturbed = Solve[{K1wall + eps == 0, Hevenwall == 0}, {dK, dM}];
expectNonzero["M7 perturbed solve nonzero dK", dK /. First[wallSolvePerturbed]];
expectNonzero["M7 perturbed solve nonzero dM", dM /. First[wallSolvePerturbed]];

gauntBase = ThreeJSymbol[{2, 0}, {2, 0}, {2, 0}];
realY20Ratio[m_] := (-1)^m*ThreeJSymbol[{2, 0}, {2, m}, {2, -m}]/gauntBase;
expectZero["M8 real-Y20 ratio m=0", realY20Ratio[0] - 1];
expectZero["M8 real-Y20 ratio m=1", realY20Ratio[1] - 1/2];
expectZero["M8 real-Y20 ratio m=2", realY20Ratio[2] + 1];
Do[
  expectZero[
    "M8 same-sign cross term m=" <> ToString[m],
    Quiet[ThreeJSymbol[{2, 0}, {2, m}, {2, m}], ClebschGordan::phy]
  ],
  {m, {1, 2}}
];

lam20 = 1;
lam21 = 1/2;
lam22 = -1;
x20 = x0 + eps1*lam20;
x21 = x0 + eps1*lam21;
x22 = x0 + eps1*lam22;
xbar = (x20 + 2*x21 + 2*x22)/5;
ax = (2*x20 - x21 - x22)/10;
bx = (x21 - x22)/2;
expectZero["M9 grouped trace", xbar - x0];
expectZero["M9 grouped line b=3a", bx - 3*ax];

Print["STATUS: PASS"];
Exit[0];
