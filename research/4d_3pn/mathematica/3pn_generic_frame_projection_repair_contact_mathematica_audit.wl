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

checkEqual[name_, lhs_, rhs_] := checkZero[name, lhs - rhs];

swapEven[expr_] := Expand[expr /. {a -> b, b -> a, c -> c, d -> e, e -> d, p -> q, q -> p}];
swapOdd[expr_] := Expand[expr /. {a -> b, b -> a, c -> c, d -> -e, e -> -d, p -> q, q -> p}];

canonicalSym[expr_, odd_] := Module[
  {s, coeffRules, denLCM, ints, g, sNorm, orderedTerms},
  s = Expand[expr + If[odd, swapOdd[expr], swapEven[expr]]];
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

generateBasis[massDeg_, velDeg_, odd_: False] := Module[
  {basis = <||>, expr, sym, tm, mp, mq, pa, pb, pc, pd, pe, key, paMax = 10, pbMax = 10, pcMax = 10},
  Do[
    mq = massDeg - mp;
    Do[
      If[2 pa + 2 pb + 2 pc + pd + pe != velDeg, Continue[]];
      expr = p^mp q^mq a^pa b^pb c^pc d^pd e^pe;
      sym = canonicalSym[expr, odd];
      tm = Expand[sym /. {b -> 0, c -> 0, e -> 0, p -> 0, q -> 1}];
      If[tm =!= 0, Continue[]];
      key = ToString[InputForm[sym]];
      basis[key] = sym;,
      {pa, 0, paMax}, {pb, 0, pbMax}, {pc, 0, pcMax}, {pd, 0, velDeg}, {pe, 0, velDeg}
    ],
    {mp, 0, massDeg}
  ];
  SortBy[Values[basis], ToString[InputForm[#]] &]
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
  res = Expand[res /. Delta -> 0];
  res
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

coeffRowFromSlots[block_, slots_] := Module[{rows = {}, poly},
  Switch[block,
    "Q",
    Do[
      poly = CoefficientRules[Expand[slots[[i]]], nu];
      rows = Join[rows, {
        Coefficient[slots[[i]], nu, 1],
        Coefficient[slots[[i]], nu, 2],
        Coefficient[slots[[i]], nu, 3]
      }],
      {i, Length[slots]}
    ],
    "T" | "S",
    Do[
      rows = Join[rows, {
        Coefficient[slots[[i]], nu, 1],
        Coefficient[slots[[i]], nu, 2]
      }],
      {i, Length[slots]}
    ],
    "U",
    rows = {Coefficient[slots[[1]], nu, 1]}
  ];
  rows
];

imageMatrixPolynomial[block_, exprs_List] := Transpose[coeffRowFromSlots[block, blockSlots[toNu[#], block]] & /@ exprs];

termsCoeffDict[expr_] := Association @ Map[Rule @@ # &, CoefficientRules[Expand[expr], {p, q, a, b, c, d, e}]];

coordinateMatrix[basis_List] := Module[{monKeys, monIndex, mat},
  monKeys = Union @@ (Keys[termsCoeffDict[#]] & /@ basis);
  monIndex = AssociationThread[monKeys, Range[Length[monKeys]]];
  mat = ConstantArray[0, {Length[monKeys], Length[basis]}];
  Do[
    KeyValueMap[(mat[[monIndex[#1], j]] = #2) &, termsCoeffDict[basis[[j]]]],
    {j, Length[basis]}
  ];
  {mat, monIndex}
];

coordsInBasis[expr_, basisMat_, monIndex_] := Module[{vec, normalMat, rhs, sol},
  vec = ConstantArray[0, Length[monIndex]];
  KeyValueMap[
    If[! KeyExistsQ[monIndex, #1], Throw[$Failed]];
    vec[[monIndex[#1]]] = #2 &,
    termsCoeffDict[expr]
  ];
  normalMat = Transpose[basisMat].basisMat;
  rhs = Transpose[basisMat].vec;
  sol = LinearSolve[normalMat, rhs];
  If[FullSimplify[Expand[basisMat.sol - vec]] =!= ConstantArray[0, Length[vec]], Throw[$Failed]];
  sol
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

rdot = d - e;
adot = -2 G q d/r^2;
bdot = 2 G p e/r^2;
cdot = G (p d - q e)/r^2;
ddot = -G q/r^2 + (a - c - d^2 + d e)/r;
edot = G p/r^2 + (c - b - d e + e^2)/r;

dtExpr[expr_] := Expand[
  D[expr, a] adot +
  D[expr, b] bdot +
  D[expr, c] cdot +
  D[expr, d] ddot +
  D[expr, e] edot +
  D[expr, r] rdot
];

DL6 = Expand[38 nu/16 - 116 nu^2/16 - 57 nu^3/16];
DL7 = Expand[20 nu^2/16 - 69 nu^3/16];
DL8 = Expand[9 nu^2/16 - 33 nu^3/16];
DL9 = Expand[5 nu^3/16];
DL10 = Expand[129 nu/16 - 98 nu^2/16 + 52 nu^3/16];
DL11 = Expand[-3 nu/16 + 52 nu^2/16 + 124 nu^3/16];
DL12 = Expand[-5 nu/12 + 11 nu^2/12 + 4 nu^3];
DL13 = Expand[-244 nu/192 - 3 Pi^2 nu/192 - 1272 nu^2/192 - 96 nu^3/192];
DL14 = Expand[452 nu/64 + 3 Pi^2 nu/64 - 6 nu^2 - 7 nu^3/2];
DL15 = Expand[(-908 + 63 Pi^2) nu/96];

targetQ = {DL6, DL7, DL8, DL9};
targetT = {DL10, DL11, DL12};
targetS = {DL13, DL14};
targetU = {DL15};

tObs = {13 nu^3/4, 31 nu^3/4, 4 nu^3};
sObs = {-nu^3/2, -7 nu^3/2};

deltaTNu = Expand[(p q/(p + q)^2) (13 a b/4 + 31 c d e/4 + 4 d^2 e^2)];
deltaSNu = Expand[(p q/(p + q)^2)^2 (c + 7 d e)/2];

QCan = Expand[(72 a^2 b - 76 a^2 c + 40 a^2 e^2 + 72 a b^2 + 122 a b c + 58 a b d e + 18 a d^2 e^2 + 15 a d e^3 - 76 b^2 c + 40 b^2 d^2 + 15 b d^3 e + 18 b d^2 e^2 - 10 d^3 e^3)/32];
TCan = Expand[(387 a^2 p^2 + 387 a^2 p q + 867 a b p^2 + 1890 a b p q + 867 a b q^2 - 9 a d^2 p^2 - 9 a d^2 p q - 129 a d e p^2 - 129 a d e p q + 387 b^2 p q + 387 b^2 q^2 - 129 b d e p q - 129 b d e q^2 - 9 b e^2 p q - 9 b e^2 q^2 + 372 c d e p q + 20 d^3 e p q + 20 d^3 e q^2 - 16 d^2 e^2 p^2 + 160 d^2 e^2 p q - 16 d^2 e^2 q^2 + 20 d e^3 p^2 + 20 d e^3 p q)/(48 (p + q))];
SCan = Expand[-(3 Pi^2 a p^4 + 880 a p^4 + 9 Pi^2 a p^3 q + 2004 a p^3 q + 9 Pi^2 a p^2 q^2 + 1368 a p^2 q^2 + 3 Pi^2 a p q^3 + 244 a p q^3 + 3 Pi^2 b p^3 q + 244 b p^3 q + 9 Pi^2 b p^2 q^2 + 1368 b p^2 q^2 + 9 Pi^2 b p q^3 + 2004 b p q^3 + 3 Pi^2 b q^4 + 880 b q^4 - 96 c p^2 q^2 - 780 d^2 p^4 - 9 Pi^2 d^2 p^4 - 2916 d^2 p^3 q - 27 Pi^2 d^2 p^3 q - 3492 d^2 p^2 q^2 - 27 Pi^2 d^2 p^2 q^2 - 1356 d^2 p q^3 - 9 Pi^2 d^2 p q^3 - 672 d e p^2 q^2 - 1356 e^2 p^3 q - 9 Pi^2 e^2 p^3 q - 3492 e^2 p^2 q^2 - 27 Pi^2 e^2 p^2 q^2 - 2916 e^2 p q^3 - 27 Pi^2 e^2 p q^3 - 780 e^2 q^4 - 9 Pi^2 e^2 q^4)/(192 (p + q)^2)];
UCan = Expand[p q (-908 + 63 Pi^2) (p + q)/96];

Qbasis = generateBasis[0, 6];
Tbasis = generateBasis[1, 4];
Sbasis = generateBasis[2, 2];
Ubasis = generateBasis[3, 0];

FQRaw = generateBasis[0, 7, True];
FTRaw = generateBasis[0, 5, True];
FSRaw = generateBasis[1, 3, True];
FURaw = generateBasis[2, 1, True];

FQBasis = r # & /@ FQRaw;
FTBasis = G # & /@ FTRaw;
FSBasis = (G^2/r) # & /@ FSRaw;
FUBasis = (G^3/r^2) # & /@ FURaw;
allF = Join[FQBasis, FTBasis, FSBasis, FUBasis];

MQ = imageMatrixPolynomial["Q", Qbasis];
MT = imageMatrixPolynomial["T", Tbasis];
MS = imageMatrixPolynomial["S", Sbasis];
MU = imageMatrixPolynomial["U", Ubasis];

{QB, Qmap} = coordinateMatrix[Qbasis];
{TB, Tmap} = coordinateMatrix[Tbasis];
{SB, Smap} = coordinateMatrix[Sbasis];
{UB, Umap} = coordinateMatrix[Ubasis];

C3 = d p + e q;
C4 = a b - c^2;
C5 = a e - c d;
C6 = b d - c e;
comNullGroebner = GroebnerBasis[{C3, C4, C5, C6}, {a, b, c, d, e, p, q}];

gammaList = {
  G (-a C5 + b C6),
  G (b C5 - a C6),
  G (d + e) (e C5 - d C6),
  -G (d - e) C4,
  -G d e (C5 - C6),
  G (-d^2 C5 + e^2 C6),
  G^2/r (q C5 - p C6),
  G^2/r (a - b) C3,
  G^2/r (-p C5 + q C6),
  G^2/r (d^2 - e^2) C3,
  -G^3/r^2 (p - q) C3
};

banner["PART I — GENERIC-FRAME COM PROJECTION RANKS"];
checkZero["Q basis count - 24", Length[Qbasis] - 24];
checkZero["T basis count - 17", Length[Tbasis] - 17];
checkZero["S basis count - 8", Length[Sbasis] - 8];
checkZero["U basis count - 1", Length[Ubasis] - 1];
checkZero["Q image rank - 12", MatrixRank[MQ] - 12];
checkZero["T image rank - 6", MatrixRank[MT] - 6];
checkZero["S image rank - 4", MatrixRank[MS] - 4];
checkZero["U image rank - 1", MatrixRank[MU] - 1];

banner["PART II — EXACT MINIMAL SEED REPAIR"];
Print["DeltaT_nu = ", Factor[deltaTNu]];
Print["DeltaS_nu = ", Factor[deltaSNu]];
checkZero["strict test-mass limit of DeltaT_nu", deltaTNu /. {b -> 0, c -> 0, e -> 0, p -> 0, q -> 1}];
checkZero["strict test-mass limit of DeltaS_nu", deltaSNu /. {b -> 0, c -> 0, e -> 0, p -> 0, q -> 1}];
Do[
  checkEqual["DeltaT_nu COM slot dL" <> ToString[i + 9], blockSlots[toNu[deltaTNu], "T"][[i]], tObs[[i]]],
  {i, 1, 3}
];
Do[
  checkEqual["DeltaS_nu COM slot dL" <> ToString[i + 12], blockSlots[toNu[deltaSNu], "S"][[i]], sObs[[i]]],
  {i, 1, 2}
];

banner["PART III — CANONICAL GENERIC-FRAME REPRESENTATIVE"];
Do[
  checkEqual["Q_can COM slot dL" <> ToString[i + 5], blockSlots[toNu[QCan], "Q"][[i]], targetQ[[i]]],
  {i, 1, 4}
];
Do[
  checkEqual["T_can COM slot dL" <> ToString[i + 9], blockSlots[toNu[TCan], "T"][[i]], targetT[[i]]],
  {i, 1, 3}
];
Do[
  checkEqual["S_can COM slot dL" <> ToString[i + 12], blockSlots[toNu[SCan], "S"][[i]], targetS[[i]]],
  {i, 1, 2}
];
checkEqual["U_can COM slot dL15", blockSlots[toNu[UCan], "U"][[1]], targetU[[1]]];

banner["PART IV — RAW ODD GENERATOR COUNTS"];
checkZero["FQ count - 31", Length[FQBasis] - 31];
checkZero["FT count - 12", Length[FTBasis] - 12];
checkZero["FS count - 8", Length[FSBasis] - 8];
checkZero["FU count - 2", Length[FUBasis] - 2];
checkZero["total odd generators - 53", Length[allF] - 53];

banner["PART V — SPARSE 11-GENERATOR CONTACT FAMILY"];
Do[
  checkZero[
    "Gamma_" <> ToString[i] <> " ideal remainder",
    Last @ PolynomialReduce[
      Numerator[Together[Expand[gammaList[[i]] /. {G -> 1, r -> 1}]]],
      comNullGroebner,
      {a, b, c, d, e, p, q}
    ]
  ],
  {i, 1, Length[gammaList]}
];

contactImages = Table[
  Module[{blk, qv, tv, sv, uv},
    blk = splitBlocks[dtExpr[gammaList[[i]]]];
    qv = ConstantArray[0, Length[Qbasis]];
    tv = ConstantArray[0, Length[Tbasis]];
    sv = ConstantArray[0, Length[Sbasis]];
    uv = ConstantArray[0, Length[Ubasis]];
    If[KeyExistsQ[blk, {1, -1}] && blk[{1, -1}] =!= 0, qv = coordsInBasis[blk[{1, -1}], QB, Qmap]];
    If[KeyExistsQ[blk, {2, -2}] && blk[{2, -2}] =!= 0, tv = coordsInBasis[blk[{2, -2}], TB, Tmap]];
    If[KeyExistsQ[blk, {3, -3}] && blk[{3, -3}] =!= 0, sv = coordsInBasis[blk[{3, -3}], SB, Smap]];
    If[KeyExistsQ[blk, {4, -4}] && blk[{4, -4}] =!= 0, uv = coordsInBasis[blk[{4, -4}], UB, Umap]];
    Join[qv, tv, sv, uv]
  ],
  {i, 1, Length[gammaList]}
];

contactImageMatrix = Transpose[contactImages];
nullFull = ArrayFlatten[{
  {MQ, ConstantArray[0, {Length[MQ], Length[Tbasis] + Length[Sbasis] + Length[Ubasis]}]},
  {ConstantArray[0, {Length[MT], Length[Qbasis]}], MT, ConstantArray[0, {Length[MT], Length[Sbasis] + Length[Ubasis]}]},
  {ConstantArray[0, {Length[MS], Length[Qbasis] + Length[Tbasis]}], MS, ConstantArray[0, {Length[MS], Length[Ubasis]}]},
  {ConstantArray[0, {Length[MU], Length[Qbasis] + Length[Tbasis] + Length[Sbasis]}], MU}
}];

checkZero["contact image rank - 11", MatrixRank[contactImageMatrix] - 11];
checkZero["full COM-null rank - 27", MatrixRank[NullSpace[nullFull]] - 27];
checkZero["contact image is COM-null", nullFull.contactImageMatrix];
checkZero["Q contact block rank - 6", MatrixRank[contactImageMatrix[[1 ;; Length[Qbasis], All]]] - 6];
checkZero["T contact block rank - 9", MatrixRank[contactImageMatrix[[Length[Qbasis] + 1 ;; Length[Qbasis] + Length[Tbasis], All]]] - 9];
checkZero["S contact block rank - 4", MatrixRank[contactImageMatrix[[Length[Qbasis] + Length[Tbasis] + 1 ;; Length[Qbasis] + Length[Tbasis] + Length[Sbasis], All]]] - 4];
checkZero["U contact block rank - 0", MatrixRank[contactImageMatrix[[-Length[Ubasis] ;;, All]]] - 0];

banner["FINAL GENERIC-FRAME PROJECTION / CONTACT LEDGER"];
Print["Passes: ", passCount];
Print["Fails: ", failCount];
If[failCount != 0, Exit[1]];

(*"
Output:

========================================================================================
PART I — GENERIC-FRAME COM PROJECTION RANKS
========================================================================================
Q basis count - 24 = 0
T basis count - 17 = 0
S basis count - 8 = 0
U basis count - 1 = 0
Q image rank - 12 = 0
T image rank - 6 = 0
S image rank - 4 = 0
U image rank - 1 = 0

========================================================================================
PART II — EXACT MINIMAL SEED REPAIR
========================================================================================
DeltaT_nu = ((13*a*b + 31*c*d*e + 16*d^2*e^2)*p*q)/(4*(p + q)^2)
DeltaS_nu = ((c + 7*d*e)*p^2*q^2)/(2*(p + q)^4)
strict test-mass limit of DeltaT_nu = 0
strict test-mass limit of DeltaS_nu = 0
DeltaT_nu COM slot dL10 = 0
DeltaT_nu COM slot dL11 = 0
DeltaT_nu COM slot dL12 = 0
DeltaS_nu COM slot dL13 = 0
DeltaS_nu COM slot dL14 = 0

========================================================================================
PART III — CANONICAL GENERIC-FRAME REPRESENTATIVE
========================================================================================
Q_can COM slot dL6 = 0
Q_can COM slot dL7 = 0
Q_can COM slot dL8 = 0
Q_can COM slot dL9 = 0
T_can COM slot dL10 = 0
T_can COM slot dL11 = 0
T_can COM slot dL12 = 0
S_can COM slot dL13 = 0
S_can COM slot dL14 = 0
U_can COM slot dL15 = 0

========================================================================================
PART IV — RAW ODD GENERATOR COUNTS
========================================================================================
FQ count - 31 = 0
FT count - 12 = 0
FS count - 8 = 0
FU count - 2 = 0
total odd generators - 53 = 0

========================================================================================
PART V — SPARSE 11-GENERATOR CONTACT FAMILY
========================================================================================
Gamma_1 ideal remainder = 0
Gamma_2 ideal remainder = 0
Gamma_3 ideal remainder = 0
Gamma_4 ideal remainder = 0
Gamma_5 ideal remainder = 0
Gamma_6 ideal remainder = 0
Gamma_7 ideal remainder = 0
Gamma_8 ideal remainder = 0
Gamma_9 ideal remainder = 0
Gamma_10 ideal remainder = 0
Gamma_11 ideal remainder = 0
"*)
