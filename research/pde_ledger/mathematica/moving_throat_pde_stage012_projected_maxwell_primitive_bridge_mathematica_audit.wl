(* moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.wl *)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[expr], Assumptions -> $Assumptions];
  Print[name, " residual = ", fmt[res]];
  If[TrueQ[res === 0],
    Print["PASS: ", name],
    Print["FAIL: ", name]; Exit[1]
  ];
];

expectNonZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[expr], Assumptions -> $Assumptions];
  Print[name, " residual = ", fmt[res]];
  If[TrueQ[res === 0],
    Print["FAIL: ", name]; Exit[1],
    Print["PASS: ", name]
  ];
];

Print["STAGE 012 PROJECTED MAXWELL PRIMITIVE BRIDGE MATHEMATICA AUDIT"];

Clear[
  Q, S2, H, Delta, P, Gw, q1, s1, h1, d1, p1, g1, ell,
  D0, kVar, bBare, zSlot, nSeed, dTarget, pGoal, shapeS, shapeT
];

$Assumptions =
  Element[
    {Q, S2, H, Delta, P, Gw, q1, s1, h1, d1, p1, g1, ell,
      D0, kVar, bBare, zSlot, nSeed, dTarget, pGoal, shapeS, shapeT},
    Reals
  ] && Delta != 0 && P != 0 && D0 != 0 && pGoal != 0 &&
    dTarget != 0 && shapeT != 0;

Z0form = Q/Delta;
Z2form = (Q S2 - H Delta)/Delta^2;
Z4form = (Q (S2^2 - Delta) - S2 H Delta)/Delta^3;
N0form = P^2/Delta^2;
N2form = 2 P (P S2 - Delta Gw)/Delta^3;
N4form =
  (Delta^2 Gw^2 - 2 Delta P^2 - 4 Delta P S2 Gw + 3 P^2 S2^2)/
    Delta^4;

expectZero["M1 Z0 primitive one-port", Z0form - Q/Delta];
expectZero["M1 Z2 primitive one-port", Z2form - (Q S2 - H Delta)/Delta^2];
expectZero[
  "M1 Z4 primitive one-port",
  Z4form - (Q (S2^2 - Delta) - S2 H Delta)/Delta^3
];
expectZero["M1 N0 primitive one-port", N0form - P^2/Delta^2];
expectZero[
  "M1 N2 primitive one-port",
  N2form - 2 P (P S2 - Delta Gw)/Delta^3
];
expectZero[
  "M1 N4 primitive one-port",
  N4form -
    (Delta^2 Gw^2 - 2 Delta P^2 - 4 Delta P S2 Gw +
        3 P^2 S2^2)/Delta^4
];

primitiveShift = {
  Q -> Q + ell q1,
  S2 -> S2 + ell s1,
  H -> H + ell h1,
  Delta -> Delta + ell d1,
  P -> P + ell p1,
  Gw -> Gw + ell g1
};

z0Series =
  FullSimplify[
    Coefficient[Normal[Series[Z0form /. primitiveShift, {ell, 0, 1}]], ell, 1],
    Assumptions -> $Assumptions
  ];
z2Series =
  FullSimplify[
    Coefficient[Normal[Series[Z2form /. primitiveShift, {ell, 0, 1}]], ell, 1],
    Assumptions -> $Assumptions
  ];
z4Series =
  FullSimplify[
    Coefficient[Normal[Series[Z4form /. primitiveShift, {ell, 0, 1}]], ell, 1],
    Assumptions -> $Assumptions
  ];
n0Series =
  FullSimplify[
    Coefficient[Normal[Series[N0form /. primitiveShift, {ell, 0, 1}]], ell, 1],
    Assumptions -> $Assumptions
  ];
n2Series =
  FullSimplify[
    Coefficient[Normal[Series[N2form /. primitiveShift, {ell, 0, 1}]], ell, 1],
    Assumptions -> $Assumptions
  ];
n4Series =
  FullSimplify[
    Coefficient[Normal[Series[N4form /. primitiveShift, {ell, 0, 1}]], ell, 1],
    Assumptions -> $Assumptions
  ];

z0Partial =
  FullSimplify[
    D[Z0form, Q] q1 + D[Z0form, S2] s1 + D[Z0form, H] h1 +
      D[Z0form, Delta] d1 + D[Z0form, P] p1 + D[Z0form, Gw] g1,
    Assumptions -> $Assumptions
  ];
z2Partial =
  FullSimplify[
    D[Z2form, Q] q1 + D[Z2form, S2] s1 + D[Z2form, H] h1 +
      D[Z2form, Delta] d1 + D[Z2form, P] p1 + D[Z2form, Gw] g1,
    Assumptions -> $Assumptions
  ];
z4Partial =
  FullSimplify[
    D[Z4form, Q] q1 + D[Z4form, S2] s1 + D[Z4form, H] h1 +
      D[Z4form, Delta] d1 + D[Z4form, P] p1 + D[Z4form, Gw] g1,
    Assumptions -> $Assumptions
  ];
n0Partial =
  FullSimplify[
    D[N0form, Q] q1 + D[N0form, S2] s1 + D[N0form, H] h1 +
      D[N0form, Delta] d1 + D[N0form, P] p1 + D[N0form, Gw] g1,
    Assumptions -> $Assumptions
  ];
n2Partial =
  FullSimplify[
    D[N2form, Q] q1 + D[N2form, S2] s1 + D[N2form, H] h1 +
      D[N2form, Delta] d1 + D[N2form, P] p1 + D[N2form, Gw] g1,
    Assumptions -> $Assumptions
  ];
n4Partial =
  FullSimplify[
    D[N4form, Q] q1 + D[N4form, S2] s1 + D[N4form, H] h1 +
      D[N4form, Delta] d1 + D[N4form, P] p1 + D[N4form, Gw] g1,
    Assumptions -> $Assumptions
  ];

z0Expected = (Delta q1 - Q d1)/Delta^2;
z2Expected =
  (-Delta^2 h1 + Delta (H d1 + Q s1 + S2 q1) - 2 Q S2 d1)/
    Delta^3;
z4Expected =
  (-Delta^2 H s1 - Delta^2 S2 h1 - Delta^2 q1 +
      2 Delta H S2 d1 + 2 Delta Q S2 s1 + 2 Delta Q d1 +
      Delta S2^2 q1 - 3 Q S2^2 d1)/Delta^4;
n0Expected = 2 P (Delta p1 - P d1)/Delta^3;
n2Expected =
  -(2 Delta^2 (Gw p1 + P g1) -
      2 Delta P (2 Gw d1 + P s1 + 2 S2 p1) +
      6 P^2 S2 d1)/Delta^4;
n4Expected =
  2 (Delta^3 Gw g1 - Delta^2 Gw^2 d1 -
      2 Delta^2 Gw P s1 - 2 Delta^2 Gw S2 p1 -
      2 Delta^2 P S2 g1 - 2 Delta^2 P p1 +
      6 Delta Gw P S2 d1 + 3 Delta P^2 S2 s1 +
      3 Delta P^2 d1 + 3 Delta P S2^2 p1 -
      6 P^2 S2^2 d1)/Delta^5;

expectZero["M2 z0 series closed form", z0Series - z0Expected];
expectZero["M2 z2 series closed form", z2Series - z2Expected];
expectZero["M2 z4 series closed form", z4Series - z4Expected];
expectZero["M2 n0 series closed form", n0Series - n0Expected];
expectZero["M2 n2 series closed form", n2Series - n2Expected];
expectZero["M2 n4 series closed form", n4Series - n4Expected];
expectZero["M2 z0 partial route", z0Partial - z0Expected];
expectZero["M2 z2 partial route", z2Partial - z2Expected];
expectZero["M2 z4 partial route", z4Partial - z4Expected];
expectZero["M2 n0 partial route", n0Partial - n0Expected];
expectZero["M2 n2 partial route", n2Partial - n2Expected];
expectZero["M2 n4 partial route", n4Partial - n4Expected];

z0 = z0Partial;
z2 = z2Partial;
z4 = z4Partial;
n0 = n0Partial;
n2 = n2Partial;
n4 = n4Partial;

xiStatic = FullSimplify[n0/N0form + z0/D0, Assumptions -> $Assumptions];
expectZero[
  "M3 static Xi1 closed form",
  xiStatic -
    (2 p1/P - 2 d1/Delta + (Delta q1 - Q d1)/(D0 Delta^2))
];

oneEquation =
  (kVar - bBare - (zSlot + ell z0)) (shapeT + ell z4) ==
    3 (shapeS + ell z2)^2;
oneSolution = kVar /. First[Solve[oneEquation, kVar]];
normEquation =
  (nSeed + ell n0)/(kVar - bBare - (zSlot + ell z0)) == pGoal;
normSolution = kVar /. First[Solve[normEquation, kVar]];

expectZero[
  "M4 one-pole K round-trip",
  ((oneSolution - bBare - (zSlot + ell z0)) (shapeT + ell z4) -
    3 (shapeS + ell z2)^2)
];
expectZero[
  "M4 normalization K round-trip",
  (nSeed + ell n0)/(normSolution - bBare - (zSlot + ell z0)) - pGoal
];

compatDirect =
  FullSimplify[
    (nSeed + ell n0)/pGoal -
      3 (shapeS + ell z2)^2/(shapeT + ell z4),
    Assumptions -> $Assumptions
  ];
compatShift =
  FullSimplify[
    Coefficient[Normal[Series[compatDirect, {ell, 0, 1}]], ell, 1],
    Assumptions -> $Assumptions
  ];
expectZero[
  "M5 fixed-target linear compatibility shift",
  compatShift - (n0/pGoal - 6 shapeS z2/shapeT + 3 shapeS^2 z4/shapeT^2)
];

transportTarget = FullSimplify[(nSeed + ell n0)/dTarget, Assumptions -> $Assumptions];
transportEquation =
  (nSeed + ell n0)/(kVar - bBare - (zSlot + ell z0)) == transportTarget;
transportSolution = kVar /. First[Solve[transportEquation, kVar]];
transportClosed = bBare + zSlot + ell z0 + dTarget;
expectZero[
  "M6 transported normalization K surface",
  transportSolution - transportClosed
];
expectZero[
  "M6 transported normalization round-trip",
  (nSeed + ell n0)/(transportSolution - bBare - (zSlot + ell z0)) -
    transportTarget
];

transportCompatDirect =
  FullSimplify[
    dTarget - 3 (shapeS + ell z2)^2/(shapeT + ell z4),
    Assumptions -> $Assumptions
  ];
transportCompatFromK =
  FullSimplify[transportSolution - oneSolution, Assumptions -> $Assumptions];
transportShift =
  FullSimplify[
    Coefficient[Normal[Series[transportCompatDirect, {ell, 0, 1}]], ell, 1],
    Assumptions -> $Assumptions
  ];
expectZero[
  "M7 transported compatibility surface",
  transportCompatFromK - transportCompatDirect
];
expectZero[
  "M7 transported linear compatibility shift",
  transportShift - (-6 shapeS z2/shapeT + 3 shapeS^2 z4/shapeT^2)
];

expectZero[
  "M8 fixed-target no z0 channel in q1",
  D[normSolution - oneSolution, q1] - D[compatDirect, q1]
];
expectZero[
  "M8 fixed-target no z0 channel in d1",
  D[normSolution - oneSolution, d1] - D[compatDirect, d1]
];
expectZero[
  "M8 transported no z0 channel in q1",
  D[transportSolution - oneSolution, q1] - D[transportCompatDirect, q1]
];
expectZero[
  "M8 transported no z0 channel in d1",
  D[transportSolution - oneSolution, d1] - D[transportCompatDirect, d1]
];
expectNonZero[
  "M8 normalization K keeps q1 z0 channel",
  D[normSolution, q1]
];
expectNonZero[
  "M8 normalization K keeps d1 z0 channel",
  D[normSolution, d1]
];

expectNonZero[
  "M9 fixed-target mutation flips z4 sign",
  compatShift - (n0/pGoal - 6 shapeS z2/shapeT - 3 shapeS^2 z4/shapeT^2)
];
expectNonZero[
  "M9 transported mutation flips z4 sign",
  transportShift - (-6 shapeS z2/shapeT - 3 shapeS^2 z4/shapeT^2)
];

Print["STATUS: PASS"];
Exit[0];
