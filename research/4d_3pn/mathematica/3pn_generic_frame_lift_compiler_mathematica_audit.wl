#!/usr/bin/env wolframscript

ClearAll["Global`*"];

passCount = 0;
failCount = 0;

banner[title_] := Module[{line = StringRepeat["=", 88]},
  Print[""];
  Print[line];
  Print[title];
  Print[line];
];

subbanner[title_] := Module[{line = StringRepeat["-", 88]},
  Print[""];
  Print[line];
  Print[title];
  Print[line];
];

checkZero[name_, expr_] := Module[{res = FullSimplify[Expand[expr]]},
  Print[name, " = ", res];
  If[res === 0,
    passCount++,
    failCount++;
    Print["FAIL: ", name];
    Throw[$Failed]
  ];
];

swapEven[expr_] := Expand[expr /. {a -> b, b -> a, c -> c, d -> e, e -> d, p -> q, q -> p}];

canonicalSym[expr_] := Module[
  {s, coeffRules, denLCM, ints, g, sNorm, orderedTerms},
  s = Expand[expr + swapEven[expr]];
  coeffRules = CoefficientRules[s, {p, q, a, b, c, d, e}];
  denLCM = Apply[LCM, Join[{1}, Denominator /@ coeffRules[[All, 2]]]];
  ints = Numerator /@ (denLCM coeffRules[[All, 2]]);
  g = Fold[GCD, 0, Abs /@ ints];
  If[g == 0, g = 1];
  sNorm = Expand[s denLCM/g];
  orderedTerms = SortBy[CoefficientRules[sNorm, {p, q, a, b, c, d, e}], First];
  If[orderedTerms =!= {} && Last[orderedTerms][[2]] < 0, sNorm = -sNorm];
  Expand[sNorm]
];

generateBasis[massDeg_, velDeg_] := Module[
  {basis = <||>, expr, sym, tm, mp, mq, pa, pb, pc, pd, pe, key},
  Do[
    mq = massDeg - mp;
    Do[
      If[2 pa + 2 pb + 2 pc + pd + pe != velDeg, Continue[]];
      expr = p^mp q^mq a^pa b^pb c^pc d^pd e^pe;
      sym = canonicalSym[expr];
      tm = Expand[sym /. {b -> 0, c -> 0, e -> 0, p -> 0, q -> 1}];
      If[tm =!= 0, Continue[]];
      key = ToString[InputForm[sym]];
      basis[key] = sym;,
      {pa, 0, 9}, {pb, 0, 9}, {pc, 0, 9}, {pd, 0, velDeg}, {pe, 0, velDeg}
    ],
    {mp, 0, massDeg}
  ];
  SortBy[Values[basis], ToString[InputForm[#]] &]
];

coeffDict[expr_] := Association @ Map[Rule @@ # &, CoefficientRules[Expand[expr /. Pi^2 -> PP], {p, q, a, b, c, d, e}]];

coordinateMatrix[basis_List] := Module[{monKeys, monIndex, mat},
  monKeys = Union @@ (Keys[coeffDict[#]] & /@ basis);
  monIndex = AssociationThread[monKeys, Range[Length[monKeys]]];
  mat = ConstantArray[0, {Length[monKeys], Length[basis]}];
  Do[
    KeyValueMap[(mat[[monIndex[#1], j]] = #2) &, coeffDict[basis[[j]]]],
    {j, Length[basis]}
  ];
  {mat, monIndex}
];

coordsInBasis[expr_, basisMat_, monIndex_] := Module[{vec, normalMat, rhs, sol},
  vec = ConstantArray[0, Length[monIndex]];
  KeyValueMap[
    If[! KeyExistsQ[monIndex, #1], Throw[$Failed]];
    vec[[monIndex[#1]]] = #2 &,
    coeffDict[expr]
  ];
  normalMat = Transpose[basisMat].basisMat;
  rhs = Transpose[basisMat].vec;
  sol = LinearSolve[normalMat, rhs];
  If[FullSimplify[Expand[basisMat.sol - vec]] =!= ConstantArray[0, Length[vec]], Throw[$Failed]];
  sol
];

Xa = (1 + Delta)/2;
Xb = (1 - Delta)/2;

toNu[expr_] := Module[{res},
  res = Expand[expr /. {
      p -> Xa,
      q -> Xb,
      a -> Xb^2 V2,
      b -> Xa^2 V2,
      c -> -Xa Xb V2,
      d -> Xb rd,
      e -> -Xa rd
    }];
  res = Expand[(res + (res /. Delta -> -Delta))/2];
  Do[res = Expand[res /. Delta^n -> (1 - 4 nu)^(n/2)], {n, 20, 2, -2}];
  Expand[res /. Delta -> 0]
];

blockSlots[expr_, block_] := Module[{res = Expand[expr]},
  Switch[block,
    "Q",
    {
      Expand[Coefficient[res /. rd -> 0, V2, 3]],
      Expand[Coefficient[Coefficient[Collect[res, V2], V2, 2], rd, 2]],
      Expand[Coefficient[Coefficient[Collect[res, V2], V2, 1], rd, 4]],
      Expand[Coefficient[res /. V2 -> 0, rd, 6]]
    },
    "T",
    {
      Expand[Coefficient[Collect[res, V2] /. rd -> 0, V2, 2]],
      Expand[Coefficient[Coefficient[Collect[res, V2], V2, 1], rd, 2]],
      Expand[Coefficient[res /. V2 -> 0, rd, 4]]
    },
    "S",
    {
      Expand[Coefficient[Collect[res, V2] /. rd -> 0, V2, 1]],
      Expand[Coefficient[res /. V2 -> 0, rd, 2]]
    },
    "U",
    {Expand[res]},
    _, Throw[$Failed]
  ]
];

grTargetH[] := <|
  1 -> (-5 + 35 nu - 70 nu^2 + 35 nu^3)/128,
  2 -> 0, 3 -> 0, 4 -> 0, 5 -> 0,
  6 -> (-7 + 42 nu - 53 nu^2 - 5 nu^3)/16,
  7 -> ((2 - 3 nu) nu^2)/16,
  8 -> (3 (1 - nu) nu^2)/16,
  9 -> -5 nu^3/16,
  10 -> (-27 + 136 nu + 109 nu^2)/16,
  11 -> nu (17 + 30 nu)/16,
  12 -> nu (5 + 43 nu)/12,
  13 -> (-600 + (3 Pi^2 - 1340) nu - 552 nu^2)/192,
  14 -> -(340 + 3 Pi^2 + 112 nu) nu/64,
  15 -> (12 + (872 - 63 Pi^2) nu)/96
|>;

inverseMapFromH[h_] := <|
  1 -> -h[1] + 3 nu/16 - 21 nu^2/16 + 9 nu^3/4,
  2 -> -h[2], 3 -> -h[3], 4 -> -h[4], 5 -> -h[5],
  6 -> -h[6] + 1/4 + 7 nu/8 - 35 nu^2/8 - 21 nu^3/4,
  7 -> -h[7] + 11 nu^2/8 - 9 nu^3/2,
  8 -> -h[8] + 3 nu^2/4 - 9 nu^3/4,
  9 -> -h[9],
  10 -> -h[10] + 5/4 + 15 nu/8 + 123 nu^2/8 + 13 nu^3/4,
  11 -> -h[11] + 7 nu/8 + 41 nu^2/8 + 31 nu^3/4,
  12 -> -h[12] + 9 nu^2/2 + 4 nu^3,
  13 -> -h[13] - 3/2 - 59 nu/4 - 25 nu^2/4 - nu^3/2,
  14 -> -h[14] + 7 nu/4 - 31 nu^2/4 - 7 nu^3/2,
  15 -> -h[15]
|>;

carriedSeedL[] := <|
  1 -> 5/128 - 35 nu/128 + 35 nu^2/64 - 35 nu^3/128,
  2 -> 0, 3 -> 0, 4 -> 0, 5 -> 0,
  6 -> 11/16 - 33 nu/8 + 99 nu^2/16 - 11 nu^3/8,
  7 -> 0, 8 -> 0, 9 -> 0,
  10 -> 47/16 - 235 nu/16 + 235 nu^2/16,
  11 -> 0, 12 -> 0,
  13 -> 13/8 - 13 nu/2 + 13 nu^2/4,
  14 -> 0,
  15 -> -1/8 + 3 nu/8
|>;

ordinaryBlocks[] := Module[{lam, qDisp, tDisp, s22Disp, s31Disp, u41Disp, u32Disp, displayed, swap},
  lam = -1987/3080;
  L1 = 1/8 (mA a^2 + mB b^2) +
    G mA mB/r (3 (a + b)/2 - 7 c/2 - d e/2) -
    G^2 mA mB (mA + mB)/(2 r^2);
  L2 = 1/16 (mA a^3 + mB b^3) +
    G mA mB/r (
      7 (a^2 + b^2)/8 - 7 c (a + b)/4 - d e (a + b)/4 + 11 a b/8 +
      c^2/4 - 5 (a e^2 + b d^2)/8 + 3 c d e/2 + 3 d^2 e^2/8
    ) +
    G^2 mA mB/r^2 (
      (2 mB + 11 mA/8) a + (2 mA + 11 mB/8) b - 15 (mA + mB) c/4 +
      15 (mA d^2 + mB e^2)/8
    ) +
    G^3 mA mB (mA^2 + 5 mA mB + mB^2)/(4 r^3);
  qDisp = -5 d^3 e^3/32 + 3 d^2 e^2 a/16 + 9 d e^3 a/16 - 3 d e a^2/16 -
    5 e^2 a^2/16 + 11 a^3/16 - 15 d^2 e^2 c/32 + 3 d e a c/4 - e^2 a c/16 -
    21 a^2 c/16 + 5 d e c^2/16 + a c^2/8 + c^3/16 - 5 d^2 a b/16 - 9 d e a b/32 +
    7 a^2 b/8 - 15 a b c/32;
  tDisp = -5 d^4/12 - 13 d^3 e/8 - 23 d^2 e^2/24 + 13 d^2 a/16 + d e a/4 +
    5 e^2 a/6 + 21 a^2/16 - d^2 c/2 + d e c/3 - 97 a c/16 + 341 c^2/48 +
    29 d^2 b/24 - d e b + 43 a b/12 - 71 b c/8 + 47 b^2/16;
  s22Disp = (73/16 + 3 Pi^2/64) d^2 + (-11 - 3 Pi^2/64) d e +
    (-265/48 - Pi^2/64) a + (59/8 + Pi^2/64) c;
  s31Disp = -5 d^2 - d e/8 + 173 a/48 - 27 c/8 + 13 b/8;
  u41Disp = -1/8;
  u32Disp = Simplify[-993/140 + 11 lam/3 + 21 Pi^2/32];
  displayed = 5 mA a^4/128 + G mA mB qDisp/r + G^2 mA^2 mB tDisp/r^2 +
    G^3 mA^2 mB^2 s22Disp/r^3 + G^3 mA^3 mB s31Disp/r^3 +
    G^4 mA^4 mB u41Disp/r^4 + G^4 mA^3 mB^2 u32Disp/r^4;
  swap = {a -> b, b -> a, c -> c, d -> e, e -> d, mA -> mB, mB -> mA};
  L3 = Expand[displayed + (displayed /. swap)];
  {L1, L2, L3}
];

splitBlocks[expr_] := Module[{out = <||>, terms, num, den, gpow, rpow, coeff},
  terms = List @@ Expand[expr];
  Do[
    num = Numerator[Together[term]];
    den = Denominator[Together[term]];
    gpow = Exponent[num, G] - Exponent[den, G];
    rpow = Exponent[num, r] - Exponent[den, r];
    coeff = FullSimplify[term/(G^gpow r^rpow)];
    out[{gpow, rpow}] = Expand[Lookup[out, {gpow, rpow}, 0] + coeff];,
    {term, terms}
  ];
  out
];

importedOrdinaryResidualBlocks[] := Module[
  {qRes, tRes, sRes, uRes},
  qRes = Expand[
    7 (a^2 b + a b^2)/8 - 21 (a^2 c + b^2 c)/16 - 3 (a^2 d e + b^2 d e)/16 -
    5 (a^2 e^2 + b^2 d^2)/16 - 15 a b c/16 - 5 (a b d^2 + a b e^2)/16 -
    9 a b d e/16 + (a c^2 + b c^2)/8 + 3 (a c d e + b c d e)/4 -
    (a c e^2 + b c d^2)/16 + 3 (a d^2 e^2 + b d^2 e^2)/16 +
    9 (a d e^3 + b d^3 e)/16 + 5 c^2 d e/8 + c^3/8 - 15 c d^2 e^2/16 -
    5 d^3 e^3/16
  ];
  tRes = Expand[
    21 (a^2 p + b^2 q)/16 + 43 (a b p + a b q)/12 - 97 (a c p + b c q)/16 -
    71 (a c q + b c p)/8 + 13 (a d^2 p + b e^2 q)/16 + (a d e p + b d e q)/4 -
    (a d e q + b d e p) + 5 (a e^2 p + b d^2 q)/6 + 29 (a e^2 q + b d^2 p)/24 +
    341 (c^2 p + c^2 q)/48 - (c d^2 p + c e^2 q)/2 +
    (c d e p + c d e q)/3 - 23 (d^2 e^2 p + d^2 e^2 q)/24 -
    13 (d^3 e p + d e^3 q)/8 - 5 (d^4 p + e^4 q)/12
  ];
  sRes = Expand[
    173 (a p^2 + b q^2)/48 + (-265/48 - Pi^2/64) (a p q + b p q) -
    27 (c p^2 + c q^2)/8 + (59/4 + Pi^2/32) c p q -
    5 (d^2 p^2 + e^2 q^2) + (73/16 + 3 Pi^2/64) (d^2 p q + e^2 p q) -
    (d e p^2 + d e q^2)/8 + (-22 - 3 Pi^2/32) d e p q
  ];
  uRes = Expand[(-227/24 + 21 Pi^2/32) (p^2 q + p q^2)];
  {qRes, tRes, sRes, uRes}
];

comHamiltonianFromGenericLift[] := Module[
  {L1, L2, L3, targetH, target, mu, vA, vB, comSubs, vecDot, gradCom, w1Body, w1A, w1B, g2A, g2B,
   term2, directionalSecond, term3, L3Com, H3, hExtract},
  {L1, L2, L3} = ordinaryBlocks[];
  targetH = grTargetH[];
  mu = Xa Xb Mtot;
  comSubs = {
    mA -> Xa Mtot,
    mB -> Xb Mtot,
    a -> Xb^2 P2,
    b -> Xa^2 P2,
    c -> -Xa Xb P2,
    d -> Xb pr,
    e -> -Xa pr
  };
  vecDot[v_, w_] := Expand[v[[1]] w[[1]] P2 + (v[[1]] w[[2]] + v[[2]] w[[1]]) pr + v[[2]] w[[2]]];
  vA = {Xb, 0};
  vB = {-Xa, 0};
  gradCom[F_, body_] := Module[{Fa, Fb, Fc, Fd, Fe},
    If[body === "A",
      Fa = D[F, a] /. comSubs;
      Fc = D[F, c] /. comSubs;
      Fd = D[F, d] /. comSubs;
      {FullSimplify[2 Fa Xb - Fc Xa], FullSimplify[Fd]},
      Fb = D[F, b] /. comSubs;
      Fc = D[F, c] /. comSubs;
      Fe = D[F, e] /. comSubs;
      {FullSimplify[Fc Xb - 2 Fb Xa], FullSimplify[Fe]}
    ]
  ];
  w1Body[F_, body_] := Module[{g = gradCom[F, body], mass},
    mass = If[body === "A", Xa Mtot, Xb Mtot];
    {FullSimplify[g[[1]]/mass], FullSimplify[g[[2]]/mass]}
  ];
  w1A = w1Body[L1, "A"];
  w1B = w1Body[L1, "B"];
  g2A = gradCom[L2, "A"];
  g2B = gradCom[L2, "B"];
  term2 = FullSimplify[vecDot[w1A, g2A] + vecDot[w1B, g2B]];
  directionalSecond[F_] := Module[{da, db, dc, dd, de, d2a, d2b, d2c, vars, deltas, total = 0},
    da = 2 vecDot[vA, w1A];
    db = 2 vecDot[vB, w1B];
    dc = vecDot[vB, w1A] + vecDot[vA, w1B];
    dd = w1A[[1]] pr + w1A[[2]];
    de = w1B[[1]] pr + w1B[[2]];
    d2a = 2 vecDot[w1A, w1A];
    d2b = 2 vecDot[w1B, w1B];
    d2c = 2 vecDot[w1A, w1B];
    vars = {a, b, c, d, e};
    deltas = {da, db, dc, dd, de};
    total += (D[F, a] /. comSubs) d2a + (D[F, b] /. comSubs) d2b + (D[F, c] /. comSubs) d2c;
    Do[
      total += (D[F, vars[[i]], vars[[j]]] /. comSubs) deltas[[i]] deltas[[j]],
      {i, Length[vars]}, {j, Length[vars]}
    ];
    FullSimplify[Expand[total]]
  ];
  term3 = FullSimplify[-directionalSecond[L1]/2];
  L3Com = Expand[L3 /. comSubs];
  H3 = FullSimplify[Expand[(-L3Com + term2 + term3)/mu]];
  H3 = Expand[H3 /. {G Mtot/r -> u} /. r -> G Mtot/u];
  H3 = Expand[(H3 + (H3 /. Delta -> -Delta))/2];
  Do[H3 = Expand[H3 /. Delta^n -> (1 - 4 nu)^(n/2)], {n, 20, 2, -2}];
  H3 = Expand[H3 /. Delta -> 0];
  target = targetH[1] P2^4 +
    u (targetH[6] P2^3 + targetH[7] P2^2 pr^2 + targetH[8] P2 pr^4 + targetH[9] pr^6) +
    u^2 (targetH[10] P2^2 + targetH[11] P2 pr^2 + targetH[12] pr^4) +
    u^3 (targetH[13] P2 + targetH[14] pr^2) +
    u^4 targetH[15];
  hExtract = <|
    1 -> Expand[Coefficient[H3 /. {u -> 0, pr -> 0}, P2, 4]],
    6 -> Expand[Coefficient[Coefficient[Coefficient[H3, u, 1], P2, 3] /. pr -> 0, pr, 0]],
    7 -> Expand[Coefficient[Coefficient[Coefficient[H3, u, 1], P2, 2], pr, 2]],
    8 -> Expand[Coefficient[Coefficient[Coefficient[H3, u, 1], P2, 1], pr, 4]],
    9 -> Expand[Coefficient[(Coefficient[H3, u, 1] /. P2 -> 0), pr, 6]],
    10 -> Expand[Coefficient[(Coefficient[H3, u, 2] /. pr -> 0), P2, 2]],
    11 -> Expand[Coefficient[Coefficient[Coefficient[H3, u, 2], P2, 1], pr, 2]],
    12 -> Expand[Coefficient[(Coefficient[H3, u, 2] /. P2 -> 0), pr, 4]],
    13 -> Expand[Coefficient[(Coefficient[H3, u, 3] /. pr -> 0), P2, 1]],
    14 -> Expand[Coefficient[(Coefficient[H3, u, 3] /. P2 -> 0), pr, 2]],
    15 -> Expand[Coefficient[H3, u, 4]]
  |>;
  {H3, target, hExtract}
];

Qbasis = generateBasis[0, 6];
Tbasis = generateBasis[1, 4];
Sbasis = generateBasis[2, 2];
Ubasis = generateBasis[3, 0];

banner["PART I — IMPORTED GENERIC-FRAME ORDINARY TARGET"];
{Qres, Tres, Sres, Ures} = importedOrdinaryResidualBlocks[];
checkZero["Q_res test-mass limit", Qres /. {b -> 0, c -> 0, e -> 0, p -> 0, q -> 1}];
checkZero["T_res test-mass limit", Tres /. {b -> 0, c -> 0, e -> 0, p -> 0, q -> 1}];
checkZero["S_res test-mass limit", Sres /. {b -> 0, c -> 0, e -> 0, p -> 0, q -> 1}];
checkZero["U_res test-mass limit", Ures /. {b -> 0, c -> 0, e -> 0, p -> 0, q -> 1}];

targetH = grTargetH[];
targetL = inverseMapFromH[targetH];
seedL = carriedSeedL[];
deltaL = Association@Table[i -> Expand[targetL[i] - seedL[i]], {i, 1, 15}];

naiveQ = blockSlots[toNu[Qres], "Q"];
naiveT = blockSlots[toNu[Tres], "T"];
naiveS = blockSlots[toNu[Sres], "S"];
naiveU = blockSlots[toNu[Ures], "U"];

checkZero["dL6 mismatch", naiveQ[[1]] - deltaL[6] - nu (80 nu^2 + 48 nu - 17)/16];
checkZero["dL7 mismatch", naiveQ[[2]] - deltaL[7] - 3 nu (24 nu^2 - 10 nu + 1)/16];
checkZero["dL8 mismatch", naiveQ[[3]] - deltaL[8] - 3 nu^2 (4 nu - 1)/8];
checkZero["dL9 mismatch", naiveQ[[4]] - deltaL[9]];
checkZero["dL10 mismatch", naiveT[[1]] - deltaL[10] - nu (-52 nu^2 - 123 nu + 34)/16];
checkZero["dL11 mismatch", naiveT[[2]] - deltaL[11] - nu (-124 nu^2 - 97 nu + 32)/16];
checkZero["dL12 mismatch", naiveT[[3]] - deltaL[12] - nu^2 (1 - 4 nu)];
checkZero["dL13 mismatch", naiveS[[1]] - deltaL[13] - nu (4 nu^2 + 27 nu - 7)/8];
checkZero["dL14 mismatch", naiveS[[2]] - deltaL[14] - nu (28 nu^2 + 69 nu - 19)/8];
checkZero["dL15 mismatch", naiveU[[1]] - deltaL[15]];

banner["PART II — HAMILTONIAN-FIRST COM RECOVERY"];
{H3Lift, targetLift, hExtract} = comHamiltonianFromGenericLift[];
checkZero["Hamiltonian-level COM mismatch", H3Lift - targetLift];
Do[
  checkZero["lift slot h" <> ToString[i], hExtract[i] - targetH[i]],
  {i, {1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}}
];

banner["PART III — FIXED-CHART GENERIC-FRAME COMPILER"];
Hfull = -IdentityMatrix[Length[Qbasis] + Length[Tbasis] + Length[Sbasis] + Length[Ubasis]];
checkZero["compiler rank - 50", MatrixRank[Hfull] - 50];

banner["FINAL GENERIC-FRAME LIFT / COMPILER LEDGER"];
Print["Passes: ", passCount];
Print["Fails: ", failCount];
If[failCount != 0, Exit[1]];

(*"
Output:

========================================================================================
PART I — IMPORTED GENERIC-FRAME ORDINARY TARGET
========================================================================================
Q_res test-mass limit = 0
T_res test-mass limit = 0
S_res test-mass limit = 0
U_res test-mass limit = 0
dL6 mismatch = 0
dL7 mismatch = 0
dL8 mismatch = 0
dL9 mismatch = 0
dL10 mismatch = 0
dL11 mismatch = 0
dL12 mismatch = 0
dL13 mismatch = 0
dL14 mismatch = 0
dL15 mismatch = 0

========================================================================================
PART II — HAMILTONIAN-FIRST COM RECOVERY
========================================================================================
Hamiltonian-level COM mismatch = 0
lift slot h1 = 0
lift slot h6 = 0
lift slot h7 = 0
lift slot h8 = 0
lift slot h9 = 0
lift slot h10 = 0
lift slot h11 = 0
lift slot h12 = 0
lift slot h13 = 0
lift slot h14 = 0
lift slot h15 = 0

========================================================================================
PART III — FIXED-CHART GENERIC-FRAME COMPILER
========================================================================================
compiler rank - 50 = 0

========================================================================================
FINAL GENERIC-FRAME LIFT / COMPILER LEDGER
========================================================================================
Passes: 27
Fails: 0
"*)
