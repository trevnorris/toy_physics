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
fmt[expr_] := StringReplace[ToString[InputForm[expr]], "Global`" -> ""];
stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  detail -> ", fmt[detail]]];
  Exit[1];
);

cleanExpr[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Cancel[Together[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

zeroQ[expr_] := And @@ (TrueQ[# === 0] & /@ Flatten[{expr}]);

expectZero[name_String, expr_] := Module[{res},
  res = cleanExpr[expr];
  Print[name, " residual = ", fmt[res]];
  If[zeroQ[res], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

expectNonZero[name_String, expr_] := Module[{res},
  res = cleanExpr[expr];
  Print[name, " residual = ", fmt[res]];
  If[zeroQ[res], fail[name, res], pass[name]];
];

num50[s_String] := ToExpression[s <> "`50"];

expectClose[name_String, actual_, expected_String, tol_: 10^-12] := Module[
  {a, e, diff},
  a = N[actual, 16];
  e = N[num50[expected], 16];
  diff = N[Abs[a - e], 16];
  Print[name, " actual = ", fmt[N[a, 16]], ", residual = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

expectVectorClose[name_String, actual_List, expected_List, tol_: 10^-12] := Module[
  {diffs},
  diffs = MapThread[N[Abs[N[#1, 16] - N[num50[#2], 16]], 16] &, {actual, expected}];
  Print[name, " residuals = ", fmt[diffs]];
  If[And @@ (TrueQ[# <= tol] & /@ diffs), pass[name], fail[name, diffs]];
];

exactDecimal[s_String] := Module[
  {txt, sign, parts, whole, frac, scale, value},
  txt = s;
  sign = 1;
  If[StringStartsQ[txt, "-"],
    sign = -1;
    txt = StringDrop[txt, 1]
  ];
  parts = StringSplit[txt, "."];
  whole = ToExpression[parts[[1]]];
  If[Length[parts] == 1,
    value = whole,
    frac = parts[[2]];
    scale = 10^StringLength[frac];
    value = (whole scale + ToExpression[frac])/scale
  ];
  sign Rationalize[value, 0]
];

epsSlope[expr_] := Coefficient[Normal[Series[expr, {eps, 0, 1}]], eps];

Clear[
  eps, lambdaA, D0, D2, D4, N0, D01, D21, D41, N01,
  u2Base, u4Base, p0Base, d0A, d2A, d4A, n0A, u2A, u4A, p0A,
  kappa, Kwall, mass, lamB, lamU, lamW, lamR, omU, omW, varpi,
  xK, xM, xLB, xV, xLU, xLW, xLR, xOU, xOW
];

$Assumptions = (
  Element[
    {
      eps, lambdaA, D0, D2, D4, N0, D01, D21, D41, N01,
      kappa, Kwall, mass, lamB, lamU, lamW, lamR, omU, omW, varpi,
      xK, xM, xLB, xV, xLU, xLW, xLR, xOU, xOW
    },
    Reals
  ]
  && D0 != 0 && N0 != 0 && lambdaA != 0
  && kappa > 0 && Kwall > 0 && mass > 0
  && lamB > 0 && lamU > 0 && lamW > 0 && lamR > 0
  && omU > 0 && omW > 0 && varpi > 0
  && omU^2 omW^2 - (kappa lamR)^2 != 0
);

banner["STAGE 225 -- MICROSCOPIC XI1 COMPILER MATHEMATICA AUDIT"];

subbanner["M1. Arbitrary-base first-order slopes from Series coefficients"];

d0A = D0 + eps lambdaA D01;
d2A = D2 + eps lambdaA D21;
d4A = D4 + eps lambdaA D41;
n0A = N0 + eps lambdaA N01;

u2Base = -D2/D0;
u4Base = (D2^2 - D0 D4)/D0^2;
p0Base = N0/D0;

u2A = -d2A/d0A;
u4A = (d2A^2 - d0A d4A)/d0A^2;
p0A = n0A/d0A;

u2Slope = cleanExpr[epsSlope[u2A]/lambdaA];
u4Slope = cleanExpr[epsSlope[u4A]/lambdaA];
pSlope = cleanExpr[epsSlope[p0A]/lambdaA];
xi1Arb = cleanExpr[pSlope/p0Base];

expectZero[
  "M1 u2^(1)",
  u2Slope - (-D0 D21 + D2 D01)/D0^2
];
expectZero[
  "M1 u4^(1)",
  u4Slope
    - (
      -D0 (D0 D41 + D01 D4 - 2 D2 D21)
      + 2 D01 (D0 D4 - D2^2)
    )/D0^3
];
expectZero["M1 Xi1", xi1Arb - (N01/N0 - D01/D0)];

subbanner["M2. Conservative compensation surface and one-pole reduction"];

d21Comp = D21 /. First[Solve[u2Slope == 0, D21]];
d41Comp = D41 /. First[Solve[cleanExpr[u4Slope /. D21 -> d21Comp] == 0, D41]];

expectZero["M2 D21 compensation", d21Comp + u2Base D01];
expectZero["M2 D41 compensation", d41Comp - (D4/D0) D01];
expectZero["M2 D41 as (u2^2-u4)D01", d41Comp - (u2Base^2 - u4Base) D01];

d4OnePole = -3 D0 u2Base^2;
onePoleD41 = cleanExpr[d41Comp /. D4 -> d4OnePole];
expectZero["M2 one-pole D41 reduction", onePoleD41 - (-3 u2Base^2) D01];
expectNonZero[
  "M2 wrong -2 coefficient is not an identity",
  onePoleD41 - (-2 u2Base^2) D01
];

subbanner["M3. Primitive logarithmic-slope compiler"];

cBase = kappa lamB;
guBase = lamU;
gwBase = kappa lamW;
rBase = kappa lamR;

deltaBase = omU^2 omW^2 - rBase^2;
s2Base = omU^2 + omW^2;
hBase = guBase^2 + gwBase^2;
qBase = guBase^2 omW^2 + 2 guBase gwBase rBase + gwBase^2 omU^2;
pBase = omU^2 gwBase + rBase guBase;

b0Base = cBase^2/varpi^2;
b2Base = cBase^2/varpi^4;
b4Base = cBase^2/varpi^6;
z0Base = qBase/deltaBase;
z2Base = (qBase s2Base - hBase deltaBase)/deltaBase^2;
z4Base = (qBase (s2Base^2 - deltaBase) - s2Base hBase deltaBase)/deltaBase^3;
n0Base = pBase^2/deltaBase^2;

kDressed = Kwall Exp[eps xK];
mDressed = mass Exp[eps xM];
lamBDressed = lamB Exp[eps xLB];
varpiDressed = varpi Exp[eps xV];
lamUDressed = lamU Exp[eps xLU];
lamWDressed = lamW Exp[eps xLW];
lamRDressed = lamR Exp[eps xLR];
omUDressed = omU Exp[eps xOU];
omWDressed = omW Exp[eps xOW];

cDressed = kappa lamBDressed;
guDressed = lamUDressed;
gwDressed = kappa lamWDressed;
rDressed = kappa lamRDressed;

deltaDressed = omUDressed^2 omWDressed^2 - rDressed^2;
s2Dressed = omUDressed^2 + omWDressed^2;
hDressed = guDressed^2 + gwDressed^2;
qDressed = (
  guDressed^2 omWDressed^2
  + 2 guDressed gwDressed rDressed
  + gwDressed^2 omUDressed^2
);
pDressed = omUDressed^2 gwDressed + rDressed guDressed;

b0Dressed = cDressed^2/varpiDressed^2;
b2Dressed = cDressed^2/varpiDressed^4;
b4Dressed = cDressed^2/varpiDressed^6;
z0Dressed = qDressed/deltaDressed;
z2Dressed = (qDressed s2Dressed - hDressed deltaDressed)/deltaDressed^2;
z4Dressed = (
  qDressed (s2Dressed^2 - deltaDressed)
  - s2Dressed hDressed deltaDressed
)/deltaDressed^3;
n0Dressed = pDressed^2/deltaDressed^2;

b0Slope = cleanExpr[epsSlope[b0Dressed]];
b2Slope = cleanExpr[epsSlope[b2Dressed]];
b4Slope = cleanExpr[epsSlope[b4Dressed]];
deltaSlope = cleanExpr[epsSlope[deltaDressed]];
s2Slope = cleanExpr[epsSlope[s2Dressed]];
hSlope = cleanExpr[epsSlope[hDressed]];
qSlope = cleanExpr[epsSlope[qDressed]];
pRawSlope = cleanExpr[epsSlope[pDressed]];
z0Slope = cleanExpr[epsSlope[z0Dressed]];
z2Slope = cleanExpr[epsSlope[z2Dressed]];
z4Slope = cleanExpr[epsSlope[z4Dressed]];
n0Slope = cleanExpr[epsSlope[n0Dressed]];

deltaClosed = 2 omU^2 omW^2 (xOU + xOW) - 2 rBase^2 xLR;
s2Closed = 2 omU^2 xOU + 2 omW^2 xOW;
hClosed = 2 guBase^2 xLU + 2 gwBase^2 xLW;
qClosed = (
  2 guBase^2 omW^2 (xLU + xOW)
  + 2 guBase gwBase rBase (xLU + xLW + xLR)
  + 2 gwBase^2 omU^2 (xLW + xOU)
);
pRawClosed = omU^2 gwBase (2 xOU + xLW) + rBase guBase (xLR + xLU);

z0Closed = (qClosed deltaBase - qBase deltaClosed)/deltaBase^2;
z2Closed = (
  deltaBase (-deltaBase hClosed - hBase deltaClosed + qBase s2Closed + s2Base qClosed)
  + 2 deltaClosed (deltaBase hBase - qBase s2Base)
)/deltaBase^3;
z4Closed = -(
  deltaBase^2 hBase s2Closed
  + deltaBase^2 s2Base hClosed
  + deltaBase^2 qClosed
  - 2 deltaBase hBase s2Base deltaClosed
  - 2 deltaBase qBase s2Base s2Closed
  - 2 deltaBase qBase deltaClosed
  - deltaBase s2Base^2 qClosed
  + 3 qBase s2Base^2 deltaClosed
)/deltaBase^4;
n0Closed = 2 pBase pRawClosed/deltaBase^2 - 2 pBase^2 deltaClosed/deltaBase^3;

expectZero["M3 B0 drift", b0Slope - b0Base (2 xLB - 2 xV)];
expectZero["M3 B2 drift", b2Slope - b2Base (2 xLB - 4 xV)];
expectZero["M3 B4 drift", b4Slope - b4Base (2 xLB - 6 xV)];
expectZero["M3 Delta drift", deltaSlope - deltaClosed];
expectZero["M3 S2 drift", s2Slope - s2Closed];
expectZero["M3 H drift", hSlope - hClosed];
expectZero["M3 Q drift", qSlope - qClosed];
expectZero["M3 P raw drift", pRawSlope - pRawClosed];
expectZero["M3 Z0 drift", z0Slope - z0Closed];
expectZero["M3 Z2 drift", z2Slope - z2Closed];
expectZero["M3 Z4 drift", z4Slope - z4Closed];
expectZero["M3 N0 drift", n0Slope - n0Closed];

d01Compiler = cleanExpr[epsSlope[kDressed - b0Dressed - z0Dressed]];
d21Compiler = cleanExpr[epsSlope[-(mDressed + b2Dressed + z2Dressed)]];
d41Compiler = cleanExpr[epsSlope[-(b4Dressed + z4Dressed)]];
n01Compiler = n0Slope;

expectZero["M3 D01 compiler", d01Compiler - (Kwall xK - b0Slope - z0Slope)];
expectZero["M3 D21 compiler", d21Compiler + (mass xM + b2Slope + z2Slope)];
expectZero["M3 D41 compiler", d41Compiler + (b4Slope + z4Slope)];
expectZero["M3 N01 compiler", n01Compiler - n0Closed];

subbanner["M4. Concrete Stage 223 compatibility sample"];

sampleRules = {
  kappa -> 2 Sqrt[2]/Pi,
  lamB -> 1/2,
  lamU -> 3/10,
  lamW -> 2/5,
  lamR -> 1/4,
  omU -> 1,
  omW -> 7/5,
  varpi -> 2,
  mass -> 1
};

b0Sample = cleanExpr[b0Base /. sampleRules];
b2Sample = cleanExpr[b2Base /. sampleRules];
b4Sample = cleanExpr[b4Base /. sampleRules];
z0Sample = cleanExpr[z0Base /. sampleRules];
z2Sample = cleanExpr[z2Base /. sampleRules];
z4Sample = cleanExpr[z4Base /. sampleRules];
n0Sample = cleanExpr[n0Base /. sampleRules];

d0Compat = cleanExpr[(3 (mass + b2Base + z2Base)^2/(b4Base + z4Base)) /. sampleRules];
kCompat = cleanExpr[b0Sample + z0Sample + d0Compat];
sampleFullRules = Join[sampleRules, {Kwall -> kCompat}];

d0Sample = cleanExpr[(Kwall - b0Base - z0Base) /. sampleFullRules];
d2Sample = cleanExpr[-(mass + b2Base + z2Base) /. sampleFullRules];
d4Sample = cleanExpr[-(b4Base + z4Base) /. sampleFullRules];
u2Sample = cleanExpr[-d2Sample/d0Sample];
u4Sample = cleanExpr[(d2Sample^2 - d0Sample d4Sample)/d0Sample^2];
p0TargetCompat = cleanExpr[n0Sample/d0Compat];

expectClose["M4 D0", d0Sample, "24.2373099886223"];
expectClose["M4 D2", d2Sample, "-1.18562046858190"];
expectClose["M4 D4", d4Sample, "-0.173991572849491"];
expectClose["M4 u2", u2Sample, "0.0489171640391802"];
expectClose["M4 u4", u4Sample, "0.00957155575054425"];
expectClose["M4 D4/D0", d4Sample/d0Sample, "-0.00717866681290820"];
expectClose["M4 P0 target", p0TargetCompat, "0.00206979231806289"];
expectZero["M4 exact one-pole identity", u4Sample - 4 u2Sample^2];

subbanner["M5. Concrete Xi1 coefficients at the sample"];

d01Sample = cleanExpr[d01Compiler /. sampleFullRules];
d21Sample = cleanExpr[d21Compiler /. sampleFullRules];
d41Sample = cleanExpr[d41Compiler /. sampleFullRules];
n01Sample = cleanExpr[n01Compiler /. sampleFullRules];
xi1Sample = cleanExpr[Expand[n01Sample/n0Sample - d01Sample/d0Sample]];

expectClose["M5 coeff xK", Coefficient[xi1Sample, xK], "-1.00975540977030"];
expectZero["M5 coeff xM", Coefficient[xi1Sample, xM]];
expectClose["M5 coeff xLB", Coefficient[xi1Sample, xLB], "0.00418038073077834"];
expectClose["M5 coeff xV", Coefficient[xi1Sample, xV], "-0.00418038073077834"];
expectClose["M5 coeff xLU", Coefficient[xi1Sample, xLU], "0.324464020216766"];
expectClose["M5 coeff xLW", Coefficient[xi1Sample, xLW], "1.69086641859305"];
expectClose["M5 coeff xLR", Coefficient[xi1Sample, xLR], "0.423379354382463"];
expectClose["M5 coeff xOU", Coefficient[xi1Sample, xOU], "-0.747843374599229"];
expectClose["M5 coeff xOW", Coefficient[xi1Sample, xOW], "-4.11424577297551"];

subbanner["M6. Wall-only and pure-BdG no-go checks"];

wallRules = {xLB -> 0, xV -> 0, xLU -> 0, xLW -> 0, xLR -> 0, xOU -> 0, xOW -> 0};
wallEquations = cleanExpr[
  {
    d21Sample + u2Sample d01Sample,
    d41Sample - (d4Sample/d0Sample) d01Sample
  } /. wallRules
];
wallMatrix = Table[Coefficient[wallEquations[[i]], v], {i, 1, 2}, {v, {xK, xM}}];
expectTrue["M6 wall matrix rank two", MatrixRank[wallMatrix] == 2];
expectTrue["M6 wall-only nullspace trivial", NullSpace[wallMatrix] === {}];

Clear[b0g, b2g, b4g, d0g, d4g, u2g];
bdgMatrix = {
  {-(b2g + u2g b0g), 2 b2g + u2g b0g},
  {-(b4g - (d4g/d0g) b0g), 3 b4g - (d4g/d0g) b0g}
};
bdgDet = Expand[Det[bdgMatrix]];
bdgExpected = -b0g b2g (d4g/d0g) - 2 b0g b4g u2g - b2g b4g;
expectZero["M6 pure-BdG determinant formula", bdgDet - bdgExpected];

bdgSampleDet = cleanExpr[
  bdgExpected /. {
    b0g -> b0Sample,
    b2g -> b2Sample,
    b4g -> b4Sample,
    d0g -> d0Sample,
    d4g -> d4Sample,
    u2g -> u2Sample
  }
];
expectClose["M6 pure-BdG sample determinant", bdgSampleDet, "-0.0000511886996120011", 10^-15];
expectNonZero["M6 pure-BdG determinant nonzero at sample", bdgSampleDet];

subbanner["M7. Mixed/U survivor from direct NullSpace"];

mixedRules = {xK -> 0, xM -> 0, xLB -> 0, xV -> 0};
mixedVars = {xLU, xLW, xLR, xOU, xOW};
mixedEquations = cleanExpr[
  {
    d21Sample + u2Sample d01Sample,
    d41Sample - (d4Sample/d0Sample) d01Sample
  } /. mixedRules
];
mixedMatrix = Table[Coefficient[mixedEquations[[i]], v], {i, 1, 2}, {v, mixedVars}];

expectTrue["M7 mixed matrix rank two", MatrixRank[mixedMatrix] == 2];
mixedNullRaw = NullSpace[mixedMatrix];
expectTrue["M7 mixed nullity three", Length[mixedNullRaw] == 3];

freeBlock = mixedNullRaw[[All, 3 ;; 5]];
mixedNullNotes = Table[
  Module[{coeffs},
    coeffs = LinearSolve[Transpose[freeBlock], UnitVector[3, j]];
    cleanExpr[coeffs . mixedNullRaw]
  ],
  {j, 1, 3}
];

Do[
  expectZero[
    "M7 normalized basis v" <> ToString[j] <> " lies in nullspace",
    mixedMatrix . mixedNullNotes[[j]]
  ],
  {j, 1, 3}
];

expectVectorClose[
  "M7 null vector v1",
  mixedNullNotes[[1]],
  {"-0.610255553634424", "0.671187016268095", "1.0", "0.0", "0.0"}
];
expectVectorClose[
  "M7 null vector v2",
  mixedNullNotes[[2]],
  {"7.05469842496522", "-9.44615143817664", "0.0", "1.0", "0.0"}
];
expectVectorClose[
  "M7 null vector v3",
  mixedNullNotes[[3]],
  {"1.61486053113911", "-0.839860892848583", "0.0", "0.0", "1.0"}
];

xi1Mixed = cleanExpr[xi1Sample /. mixedRules];
xi1NullValues = Table[
  cleanExpr[xi1Mixed /. Thread[mixedVars -> mixedNullNotes[[j]]]],
  {j, 1, 3}
];

expectClose["M7 Xi1(v1)", xi1NullValues[[1]], "1.36026097049402"];
expectClose["M7 Xi1(v2)", xi1NullValues[[2]], "-14.4310278139755"];
expectClose["M7 Xi1(v3)", xi1NullValues[[3]], "-5.01037421295998"];

subbanner["M8. Transported Stage 224 amplitude windows"];

sigma1Window = exactDecimal["1.36026097049402"];
budgetData = {
  {"both_10", exactDecimal["0.367930328492646"], "0.270485102839510"},
  {"one_10", exactDecimal["0.737619063660757"], "0.542262903708006"},
  {"both_30", exactDecimal["2.94889585703134"], "2.16788978070904"},
  {"one_30", exactDecimal["4.63505472371892"], "3.40747461278373"}
};

Do[
  Module[{label, budget, expectedWindow, window},
    {label, budget, expectedWindow} = entry;
    window = budget/sigma1Window;
    expectClose["M8 " <> label <> " window", window, expectedWindow];
  ],
  {entry, budgetData}
];

Print[""];
Print["All Stage 225 Mathematica checks passed."];
