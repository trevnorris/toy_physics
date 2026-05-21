(* Projected Maxwell -> grouped P2 bridge

   Independent Mathematica audit of the projection-first Maxwell corrections
   against the grouped real-P2 / full-bundle algebra used in the moving-throat
   PDE program.
*)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

check[label_, expr_] := If[
  FullSimplify[expr] =!= 0,
  Print["FAIL: ", label, " residual = ", fmt[FullSimplify[expr]]]; Exit[1],
  Print[label, " residual = ", fmt[FullSimplify[expr]]]
];

Print["STAGE 011 PROJECTED MAXWELL P2 BRIDGE MATHEMATICA AUDIT"];

Clear[
  D0, D2, D4, N0, N2, N4, z0, z2, z4, n0, n2, n4, eps,
  B4, Z4, M, B2, Z2, K, B0, Z0slot, Ptarget, S, T, D0target,
  ea, x0, x1, lam, za, na, z2a
];

$Assumptions =
  Element[
    {D0, D2, D4, N0, N2, N4, z0, z2, z4, n0, n2, n4, eps,
      B4, Z4, M, B2, Z2, K, B0, Z0slot, Ptarget, S, T, D0target,
      ea, x0, x1, lam, za, na, z2a},
    Reals
  ] && D0 != 0 && N0 != 0 && T != 0 && Ptarget != 0 &&
    D0target != 0 && lam != 0;

momentD[0, e_] := D0 - e z0;
momentD[2, e_] := D2 - e z2;
momentD[4, e_] := D4 - e z4;
momentN[0, e_] := N0 + e n0;
momentN[2, e_] := N2 + e n2;
momentN[4, e_] := N4 + e n4;

u2Bundle[e_] := -momentD[2, e]/momentD[0, e];
u4Bundle[e_] := (momentD[2, e]^2 - momentD[0, e] momentD[4, e])/momentD[0, e]^2;
p0Bundle[e_] := momentN[0, e]/momentD[0, e];
p2Bundle[e_] :=
  (momentD[0, e] momentN[2, e] - 2 momentD[2, e] momentN[0, e])/
    momentD[0, e]^2;
p4Bundle[e_] :=
  (momentD[0, e]^2 momentN[4, e] -
      2 momentD[0, e] (momentD[2, e] momentN[2, e] +
        momentD[4, e] momentN[0, e]) +
      3 momentD[2, e]^2 momentN[0, e])/momentD[0, e]^3;

firstShift[expr_, var_] :=
  Coefficient[Normal[Series[expr - (expr /. var -> 0), {var, 0, 1}]], var, 1];

deltaU2 = firstShift[u2Bundle[eps], eps];
check["M1 delta u2", deltaU2 - (D0 z2 - D2 z0)/D0^2];

deltaP0 = firstShift[p0Bundle[eps], eps];
check["M2 delta P0", deltaP0 - (D0 n0 + N0 z0)/D0^2];

poleDefect[e_] :=
  (D0 - e z0) (B4 + Z4 + e z4) -
    3 (M + B2 + Z2 + e z2)^2;
deltaPole = firstShift[poleDefect[eps], eps];
check[
  "M3 one-pole first variation",
  deltaPole - (D0 z4 - z0 (B4 + Z4) - 6 z2 (M + B2 + Z2))
];

kPole =
  K /. First[
    Solve[(K - B0 - (Z0slot + eps z0)) (T + eps z4) ==
      3 (S + eps z2)^2, K]
  ];
kNorm =
  K /. First[
    Solve[(N0 + eps n0)/(K - B0 - (Z0slot + eps z0)) == Ptarget, K]
  ];
compatFixed = FullSimplify[kNorm - kPole];
compatFixedClosed =
  (N0 + eps n0)/Ptarget - 3 (S + eps z2)^2/(T + eps z4);
check[
  "M4 fixed-target K-eliminated surface",
  compatFixed - compatFixedClosed
];

compatFixedShift = firstShift[compatFixed, eps];
check[
  "M5 fixed-target compatibility first variation",
  compatFixedShift - (n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2)
];
check["M5 fixed-target z0 independence", D[compatFixedShift, z0]];

pTargetTransport = (N0 + eps n0)/D0target;
kNormTransport =
  K /. First[
    Solve[
      (N0 + eps n0)/(K - B0 - (Z0slot + eps z0)) == pTargetTransport,
      K
    ]
  ];
check[
  "M6 transported-target normalization K surface",
  kNormTransport - (B0 + Z0slot + eps z0 + D0target)
];

compatTransport = FullSimplify[kNormTransport - kPole];
compatTransportShift = firstShift[compatTransport, eps];
check[
  "M7 transported-target compatibility first variation",
  compatTransportShift - (-6 S z2/T + 3 S^2 z4/T^2)
];
check["M7 transported-target z0 independence", D[compatTransportShift, z0]];

g[m1_, m2_, m3_] :=
  Quiet[
    Sqrt[5/(4 Pi)] *
      ThreeJSymbol[{2, 0}, {2, 0}, {2, 0}] *
      ThreeJSymbol[{2, m1}, {2, m2}, {2, m3}],
    ClebschGordan::phy
  ];

laneLambda[0] := 1;
laneLambda[m_Integer] := FullSimplify[(-1)^m g[0, m, -m]/g[0, 0, 0]];

check["M8 real-Y20 lambda(0)", laneLambda[0] - 1];
check["M8 real-Y20 lambda(1)", laneLambda[1] - 1/2];
check["M8 real-Y20 lambda(2)", laneLambda[2] + 1];
check["M8 same-sign selection m=1", g[0, 1, 1]];
check["M8 same-sign selection m=2", g[0, 2, 2]];

x20 = x0 + ea laneLambda[0] x1;
x21 = x0 + ea laneLambda[1] x1;
x22 = x0 + ea laneLambda[2] x1;
xbar = FullSimplify[(x20 + 2 x21 + 2 x22)/5];
ax = FullSimplify[(2 x20 - x21 - x22)/10];
bx = FullSimplify[(x21 - x22)/2];
check["M9 weak-axisymmetric grouped trace", xbar - x0];
check["M9 weak-axisymmetric b=3a", bx - 3 ax];

pLane = (N0 + ea lam na)/(D0 - ea lam za);
xiSlope = FullSimplify[((D[pLane, ea] /. ea -> 0)/lam)/(N0/D0)];
check["M10 static Xi1 slope", xiSlope - (na/N0 + za/D0)];

u2Lane = -(D2 - ea lam z2a)/(D0 - ea lam za);
u2LaneSlope = FullSimplify[(D[u2Lane, ea] /. ea -> 0)/lam];
check[
  "M11 u2 projected-Maxwell lane slope",
  u2LaneSlope - (D0 z2a - D2 za)/D0^2
];

Print["STATUS: PASS"];
Exit[0];
