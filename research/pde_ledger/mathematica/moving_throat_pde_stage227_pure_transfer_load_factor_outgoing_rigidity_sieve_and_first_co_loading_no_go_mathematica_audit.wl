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

fmt[expr_] := StringReplace[ToString[InputForm[expr]], "Global`" -> ""];
num[s_String] := ToExpression[s <> "`60"];
stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

expectZero[tag_String, expr_] := Module[{res},
  res = FullSimplify[Cancel[Together[expr]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[tag, " residual = ", fmt[res]];
  If[TrueQ[res === 0],
    Print["OK  ", tag],
    Print["FAIL ", tag, " -> ", fmt[res]];
    Exit[1]
  ];
];

expectTrue[tag_String, cond_] := (
  Print[tag, " condition = ", fmt[cond]];
  If[TrueQ[cond],
    Print["OK  ", tag],
    Print["FAIL ", tag];
    Exit[1]
  ];
);

expectClose[tag_String, actual_, expected_, tol_] := Module[{aa, ee, diff},
  aa = N[actual, 50];
  ee = N[expected, 50];
  diff = Abs[aa - ee];
  Print[tag, " actual = ", fmt[aa], " expected = ", fmt[ee], " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol],
    Print["OK  ", tag],
    Print["FAIL ", tag];
    Exit[1]
  ];
];

expectVectorClose[tag_String, actual_List, expected_List, tol_] := Module[
  {aa, ee, diff},
  aa = N[actual, 50];
  ee = SetPrecision[expected, 50];
  diff = Max[Abs[aa - ee]];
  Print[tag, " actual = ", fmt[N[aa, 16]], " maxdiff = ", fmt[N[diff, 16]]];
  If[TrueQ[diff <= tol],
    Print["OK  ", tag],
    Print["FAIL ", tag, " expected = ", fmt[ee]];
    Exit[1]
  ];
];

expectDirectionClose[tag_String, actual_List, expected_List, tol_] := Module[
  {aa, ee, diff},
  aa = N[actual, 50];
  ee = SetPrecision[expected, 50];
  diff = Min[Max[Abs[aa - ee]], Max[Abs[aa + ee]]];
  Print[tag, " actual = ", fmt[N[aa, 16]], " signed maxdiff = ", fmt[N[diff, 16]]];
  If[TrueQ[diff <= tol],
    Print["OK  ", tag],
    Print["FAIL ", tag, " expected direction = ", fmt[ee]];
    Exit[1]
  ];
];

coeffRow[expr_, vars_List] := (Coefficient[Expand[expr], #] & /@ vars);
unit[v_List] := FullSimplify[v/Sqrt[v.v], Assumptions -> $Assumptions];

Clear[
  kap, kStiff, mass, lamB, lamU, lamW, lamR, omegaU, omegaW, varpi,
  xlu, xlw, xlr, xou, xow,
  gUSym, gWSym, rSym, omUSym, omWSym
];

$Assumptions = (
  Element[
    {
      kap, kStiff, mass, lamB, lamU, lamW, lamR, omegaU, omegaW, varpi,
      xlu, xlw, xlr, xou, xow,
      gUSym, gWSym, rSym, omUSym, omWSym
    },
    Reals
  ]
  && kap > 0 && kStiff > 0 && mass > 0
  && lamB > 0 && lamU > 0 && lamW > 0 && lamR > 0
  && omegaU > 0 && omegaW > 0 && varpi > 0
  && gUSym > 0 && gWSym > 0 && rSym > 0 && omUSym > 0 && omWSym > 0
  && omegaU^2 omegaW^2 - (kap lamR)^2 != 0
  && omUSym^2 omWSym^2 - rSym^2 != 0
);

banner["STAGE 227 -- PURE-TRANSFER LOAD FACTOR AND CO-LOADING NO-GO MATHEMATICA AUDIT"];

mixedVars = {xlu, xlw, xlr, xou, xow};
logPairs = {{xlu, lamU}, {xlw, lamW}, {xlr, lamR}, {xou, omegaU}, {xow, omegaW}};
deltaMixed[expr_] := FullSimplify[
  Total[(#[[1]] #[[2]] D[expr, #[[2]]]) & /@ logPairs],
  Assumptions -> $Assumptions
];

sampleRules = {
  kap -> 2 Sqrt[2]/Pi,
  lamB -> 1/2,
  lamU -> 3/10,
  lamW -> 2/5,
  lamR -> 1/4,
  omegaU -> 1,
  omegaW -> 7/5,
  varpi -> 2,
  mass -> 1
};

cWall = kap lamB;
gU = lamU;
gW = kap lamW;
rMix = kap lamR;
delta = omegaU^2 omegaW^2 - rMix^2;
s2 = omegaU^2 + omegaW^2;
hMix = gU^2 + gW^2;
qMix = gU^2 omegaW^2 + 2 gU gW rMix + gW^2 omegaU^2;
pPort = omegaU^2 gW + rMix gU;
n0Same = pPort^2/delta^2;

b0 = cWall^2/varpi^2;
b2 = cWall^2/varpi^4;
b4 = cWall^2/varpi^6;
z0 = qMix/delta;
z2 = (qMix s2 - hMix delta)/delta^2;
z4 = (qMix (s2^2 - delta) - s2 hMix delta)/delta^3;
d0Bundle = kStiff - b0 - z0;
d2Bundle = -(mass + b2 + z2);
d4Bundle = -(b4 + z4);

d0Compat = FullSimplify[
  3 (1 + b2 + z2)^2/(b4 + z4) /. sampleRules,
  Assumptions -> $Assumptions
];
kCompat = FullSimplify[(b0 + z0 /. sampleRules) + d0Compat, Assumptions -> $Assumptions];
Print["Stage-225 K compatibility value = ", fmt[N[kCompat, 20]]];

p01Mixed = Expand[FullSimplify[deltaMixed[pPort] /. sampleRules, Assumptions -> $Assumptions]];
delta01Mixed = Expand[FullSimplify[deltaMixed[delta] /. sampleRules, Assumptions -> $Assumptions]];
n01Mixed = Expand[FullSimplify[deltaMixed[n0Same] /. sampleRules, Assumptions -> $Assumptions]];
pSample = FullSimplify[pPort /. sampleRules, Assumptions -> $Assumptions];
deltaSample = FullSimplify[delta /. sampleRules, Assumptions -> $Assumptions];
n0Sample = FullSimplify[n0Same /. sampleRules, Assumptions -> $Assumptions];
xi1Mixed = Expand[FullSimplify[n01Mixed/n0Sample, Assumptions -> $Assumptions]];

subbanner["M1. Pure-transfer theorem"];

loadLawResidual = xi1Mixed - 2 (p01Mixed/pSample - delta01Mixed/deltaSample);
Print["Xi1 from direct delta N0/N0 = ", fmt[xi1Mixed]];
Print["2(delta P/P - delta Delta/Delta) = ", fmt[
  Expand[FullSimplify[2 (p01Mixed/pSample - delta01Mixed/deltaSample), Assumptions -> $Assumptions]]
]];
expectZero["M1 Xi1 direct N0 variation equals load-factor variation", loadLawResidual];

subbanner["M2. One-port factorization"];

lambdaSym = (omUSym^2 gWSym + rSym gUSym)/(omUSym^2 omWSym^2 - rSym^2);
iSym = rSym gUSym/(omUSym^2 gWSym);
hSym = rSym^2/(omUSym^2 omWSym^2);
lambdaFactored = (gWSym/omWSym^2) (1 + iSym)/(1 - hSym);
Print["Lambda symbolic residual = ", fmt[
  FullSimplify[Cancel[Together[lambdaSym - lambdaFactored]], Assumptions -> $Assumptions]
]];
expectZero["M2 Lambda equals one-port factorization", lambdaSym - lambdaFactored];

subbanner["M3. Microscopic slope law and sample coefficients"];

iPrimitive = rMix gU/(omegaU^2 gW);
hPrimitive = rMix^2/(omegaU^2 omegaW^2);
mForm = Expand[FullSimplify[deltaMixed[Log[gW/omegaW^2]], Assumptions -> $Assumptions]];
iForm = Expand[FullSimplify[deltaMixed[Log[iPrimitive]], Assumptions -> $Assumptions]];
hForm = Expand[FullSimplify[deltaMixed[Log[hPrimitive]], Assumptions -> $Assumptions]];
iSample = FullSimplify[iPrimitive /. sampleRules, Assumptions -> $Assumptions];
hSample = FullSimplify[hPrimitive /. sampleRules, Assumptions -> $Assumptions];
xiSlopeLaw = Expand[FullSimplify[
  2 (mForm + iSample/(1 + iSample) iForm + hSample/(1 - hSample) hForm),
  Assumptions -> $Assumptions
]];
xiSpecial = Expand[FullSimplify[
  2 mForm + (6/19) iForm + 50 hForm/(98 Pi^2 - 25),
  Assumptions -> $Assumptions
]];

Print["m = ", fmt[mForm]];
Print["i = ", fmt[iForm]];
Print["h = ", fmt[hForm]];
Print["I_sample = ", fmt[iSample]];
Print["H_sample = ", fmt[hSample]];
Print["Xi1 sample slope law = ", fmt[xiSpecial]];
expectZero["M3 m log-form", mForm - (xlw - 2 xow)];
expectZero["M3 i log-form", iForm - (xlr + xlu - xlw - 2 xou)];
expectZero["M3 h log-form", hForm - (2 xlr - 2 xou - 2 xow)];
expectZero["M3 I sample equals 3/16", iSample - 3/16];
expectZero["M3 H sample equals 25/(98 Pi^2)", hSample - 25/(98 Pi^2)];
expectZero["M3 Xi1 equals microscopic slope law", xi1Mixed - xiSlopeLaw];
expectZero["M3 Xi1 equals specialized sample law", xi1Mixed - xiSpecial];

subbanner["M4. Combined i=h rigidity determinant"];

d01Mixed = Expand[FullSimplify[deltaMixed[d0Bundle] /. sampleRules, Assumptions -> $Assumptions]];
d21Mixed = Expand[FullSimplify[deltaMixed[d2Bundle] /. sampleRules, Assumptions -> $Assumptions]];
d41Mixed = Expand[FullSimplify[deltaMixed[d4Bundle] /. sampleRules, Assumptions -> $Assumptions]];
transferMatrix = {
  coeffRow[d01Mixed, mixedVars],
  coeffRow[d21Mixed, mixedVars],
  coeffRow[d41Mixed, mixedVars]
};

expectTrue["M4 transfer matrix has rank 3", MatrixRank[transferMatrix] === 3];
transferKernel = NullSpace[transferMatrix];
expectTrue["M4 transfer nullity is 2", Length[transferKernel] === 2];

kernelColumns = Transpose[transferKernel];
tailSolveMatrix = kernelColumns[[{4, 5}, All]];
expectTrue["M4 nullspace tail coordinates are independent", Det[tailSolveMatrix] =!= 0];
t1 = FullSimplify[kernelColumns.LinearSolve[tailSolveMatrix, {1, 0}], Assumptions -> $Assumptions];
t2 = FullSimplify[kernelColumns.LinearSolve[tailSolveMatrix, {0, 1}], Assumptions -> $Assumptions];
transferBasis = {t1, t2};
transferBasisColumns = Transpose[transferBasis];

expectedT1 = {-4.359222794718, 3.107402039105, 18.703510605854, 1.0, 0.0};
expectedT2 = {1.909256655687, -1.163651238154, -0.482414494705, 0.0, 1.0};
Print["t1 = ", fmt[N[t1, 16]]];
Print["t2 = ", fmt[N[t2, 16]]];
expectVectorClose["M4 Stage-226 transfer basis t1", t1, expectedT1, 10^-11];
expectVectorClose["M4 Stage-226 transfer basis t2", t2, expectedT2, 10^-11];

iRow = coeffRow[iForm, mixedVars];
hRow = coeffRow[hForm, mixedVars];
mRow = coeffRow[mForm, mixedVars];
reducedIH = {iRow, hRow}.transferBasisColumns;
detIH = FullSimplify[Det[reducedIH], Assumptions -> $Assumptions];
Print["det[(i,h)|pure transfer] = ", fmt[detIH]];
expectTrue["M4 combined i=h determinant is structurally nonzero", FullSimplify[detIH, Assumptions -> $Assumptions] =!= 0];
expectTrue["M4 combined i=h reduced matrix has rank 2", MatrixRank[reducedIH] === 2];

subbanner["M5. One-dimensional rigid survivors"];

reducedI = {iRow}.transferBasisColumns;
reducedH = {hRow}.transferBasisColumns;
reducedM = {mRow}.transferBasisColumns;
expectTrue["M5 i=0 reduced nullity is 1", Length[NullSpace[reducedI]] === 1];
expectTrue["M5 h=0 reduced nullity is 1", Length[NullSpace[reducedH]] === 1];
expectTrue["M5 m=0 reduced nullity is 1", Length[NullSpace[reducedM]] === 1];

vi = unit[FullSimplify[transferBasisColumns.NullSpace[reducedI][[1]], Assumptions -> $Assumptions]];
vh = unit[FullSimplify[transferBasisColumns.NullSpace[reducedH][[1]], Assumptions -> $Assumptions]];
vm = unit[FullSimplify[transferBasisColumns.NullSpace[reducedM][[1]], Assumptions -> $Assumptions]];

expectedVi = {0.45280825, -0.29424612, -0.82815170, -0.04054866, 0.14458380};
expectedVh = {0.66561963, -0.38941932, 0.46712837, 0.03609301, 0.43103536};
expectedVm = {0.13386239, -0.10586713, -0.98242900, -0.05389175, -0.05293356};

expectDirectionClose["M5 i=0 unit survivor direction", vi, expectedVi, 10^-8];
expectDirectionClose["M5 h=0 unit survivor direction", vh, expectedVh, 10^-8];
expectDirectionClose["M5 m=0 unit survivor direction", vm, expectedVm, 10^-8];

subbanner["M6. Same-charge gain scales"];

xiCoeff = coeffRow[xi1Mixed, mixedVars];
sigmaI = Abs[N[xiCoeff.vi, 50]];
sigmaH = Abs[N[xiCoeff.vh, 50]];
sigmaM = Abs[N[xiCoeff.vm, 50]];
expectClose["M6 |Xi1(v_i)|", sigmaI, num["1.26576248"], 10^-8];
expectClose["M6 |Xi1(v_h)|", sigmaH, num["2.04509123"], 10^-8];
expectClose["M6 |Xi1(v_m)|", sigmaM, num["0.29342952"], 10^-8];

subbanner["M7. Corridor norm and transported 10%-loss ceilings"];

gramTransfer = FullSimplify[Transpose[transferBasisColumns].transferBasisColumns, Assumptions -> $Assumptions];
projectorTransfer = FullSimplify[
  transferBasisColumns.Inverse[gramTransfer].Transpose[transferBasisColumns],
  Assumptions -> $Assumptions
];
sigmaTransfer = Sqrt[N[xiCoeff.projectorTransfer.xiCoeff, 70]];
budgetBoth = num["0.367930328492646"];
budgetNonempty = num["0.737619063660757"];

ceilTransfer = {budgetBoth/sigmaTransfer, budgetNonempty/sigmaTransfer};
ceilI = {budgetBoth/sigmaI, budgetNonempty/sigmaI};
ceilH = {budgetBoth/sigmaH, budgetNonempty/sigmaH};
ceilM = {budgetBoth/sigmaM, budgetNonempty/sigmaM};

Print["sigma_transfer = ", fmt[N[sigmaTransfer, 20]]];
Print["transfer ceilings = ", fmt[N[ceilTransfer, 16]]];
Print["i ceilings = ", fmt[N[ceilI, 16]]];
Print["h ceilings = ", fmt[N[ceilH, 16]]];
Print["m ceilings = ", fmt[N[ceilM, 16]]];

expectClose["M7 sigma_transfer", sigmaTransfer, num["2.31561904386057"], 10^-11];
expectClose["M7 transfer both ceiling", ceilTransfer[[1]], num["0.15889070"], 10^-8];
expectClose["M7 transfer nonempty ceiling", ceilTransfer[[2]], num["0.31854077"], 10^-8];
expectClose["M7 i both ceiling", ceilI[[1]], num["0.29067881"], 10^-8];
expectClose["M7 i nonempty ceiling", ceilI[[2]], num["0.58274682"], 10^-8];
expectClose["M7 h both ceiling", ceilH[[1]], num["0.17990900"], 10^-8];
expectClose["M7 h nonempty ceiling", ceilH[[2]], num["0.36067783"], 10^-8];
expectClose["M7 m both ceiling", ceilM[[1]], num["1.25389678"], 10^-8];
expectClose["M7 m nonempty ceiling", ceilM[[2]], num["2.51378617"], 10^-8];

Print[""];
Print["All Stage 227 Mathematica checks passed."];
Exit[0];
