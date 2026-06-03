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
num[s_String] := ToExpression[s <> "`50"];
stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

expectZero[expr_, tag_] := Module[{res, nres, ok},
  res = stripConditional[FullSimplify[Cancel[Together[expr]], Assumptions -> $Assumptions]];
  res = stripConditional[res];
  nres = Quiet[N[res, 30]];
  ok = TrueQ[PossibleZeroQ[Chop[nres]]] ||
    If[NumberQ[nres], TrueQ[Abs[nres] < 10^-12], False];
  If[ok,
    Print["OK  ", tag],
    Print["FAIL ", tag, " -> ", nres];
    Exit[1]
  ];
];

expectClose[a_, b_, tol_: 10^-11, tag_] := Module[{aa, bb, diff},
  aa = N[a, 30];
  bb = N[b, 30];
  diff = Abs[aa - bb];
  If[TrueQ[diff <= tol],
    Print["OK  ", tag],
    Print["FAIL ", tag, " -> ", aa, " vs ", bb];
    Exit[1]
  ];
];

coeffRow[expr_, vars_List] := Normal[CoefficientArrays[Expand[expr], vars][[2]]];

projectedNorm[coeff_List, basis_List] := Module[{qRows, projector},
  qRows = Orthogonalize[N[basis, 80]];
  projector = Transpose[qRows].qRows;
  Sqrt[N[coeff, 80].projector.N[coeff, 80]]
];

Clear[
  eps, lam, d0, d2, d4, n0, d01, d21, d41, n01,
  u2var, u4var, kap, kStiff, mass, lamB, lamU, lamW, lamR,
  omegaU, omegaW, varpi, xk, xm, xlb, xv, xlu, xlw, xlr, xou, xow
];

$Assumptions = (
  Element[
    {
      eps, lam, d0, d2, d4, n0, d01, d21, d41, n01,
      u2var, u4var, kap, kStiff, mass, lamB, lamU, lamW, lamR,
      omegaU, omegaW, varpi, xk, xm, xlb, xv, xlu, xlw, xlr, xou, xow
    },
    Reals
  ]
  && lam != 0 && d0 != 0 && n0 != 0
  && kap > 0 && kStiff > 0 && mass > 0
  && lamB > 0 && lamU > 0 && lamW > 0 && lamR > 0
  && omegaU > 0 && omegaW > 0 && varpi > 0
);

banner["STAGE 226 -- STRICT 5PN EVEN-GATE PACKAGE MATHEMATICA AUDIT"];

subbanner["M1. Load-defect bridge from first-order series"];

n0A = n0 + eps lam n01;
d0A = d0 + eps lam d01;
p0 = n0/d0;
p0A = n0A/d0A;
p0Series = Normal[Series[p0A, {eps, 0, 1}]];
xi1Series = FullSimplify[Coefficient[p0Series, eps]/(lam p0), Assumptions -> $Assumptions];
xiLoad = n01/n0 - d01/d0;
Print["Xi1(series) = ", fmt[xi1Series]];
expectZero[xi1Series - xiLoad, "M1 Xi1 equals load-defect bridge"];

subbanner["M2. Even-gate closed forms on the compensation surface"];

u2 = -d2/d0;
k1Gate = d21 + d01/9;
hEvenGate = d41 - (2/3) d21 - d01/27;
k1Comp = FullSimplify[k1Gate /. d21 -> -u2 d01, Assumptions -> $Assumptions];
hEvenComp = FullSimplify[
  hEvenGate /. {d21 -> -u2 d01, d41 -> (d4/d0) d01},
  Assumptions -> $Assumptions
];
hOnePole = FullSimplify[
  ((u2var^2 - u4var) + (2/3) u2var - 1/27) d01 /. u4var -> 4 u2var^2,
  Assumptions -> $Assumptions
];

Print["K1(comp) = ", fmt[k1Comp]];
Print["H_even(comp) = ", fmt[hEvenComp]];
Print["H_even(one-pole) = ", fmt[hOnePole]];
expectZero[k1Comp - (1/9 - u2) d01, "M2 K1 compensation formula"];
expectZero[
  hEvenComp - (d4/d0 + (2/3) u2 - 1/27) d01,
  "M2 H_even compensation formula"
];
expectZero[
  hOnePole - (-3 u2var^2 + (2/3) u2var - 1/27) d01,
  "M2 H_even one-pole formula"
];

subbanner["M3. Compatibility branch from physical primitives"];

cWall = kap lamB;
gU = lamU;
gW = kap lamW;
rMix = kap lamR;
delta = omegaU^2 omegaW^2 - rMix^2;
s2 = omegaU^2 + omegaW^2;
h2 = gU^2 + gW^2;
qMix = gU^2 omegaW^2 + 2 gU gW rMix + gW^2 omegaU^2;
pTrans = omegaU^2 gW + rMix gU;

b0 = cWall^2/varpi^2;
b2 = cWall^2/varpi^4;
b4 = cWall^2/varpi^6;
z0 = qMix/delta;
z2 = (qMix s2 - h2 delta)/delta^2;
z4 = (qMix (s2^2 - delta) - s2 h2 delta)/delta^3;
n0Bundle = pTrans^2/delta^2;

d0Bundle = kStiff - b0 - z0;
d2Bundle = -(mass + b2 + z2);
d4Bundle = -(b4 + z4);

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

d0Compat = FullSimplify[
  3 (1 + b2 + z2)^2/(b4 + z4) /. sampleRules,
  Assumptions -> $Assumptions
];
kCompat = FullSimplify[(b0 + z0 /. sampleRules) + d0Compat, Assumptions -> $Assumptions];
sampleFull = Join[sampleRules, {kStiff -> kCompat}];

d0Value = FullSimplify[d0Bundle /. sampleFull, Assumptions -> $Assumptions];
d2Value = FullSimplify[d2Bundle /. sampleFull, Assumptions -> $Assumptions];
d4Value = FullSimplify[d4Bundle /. sampleFull, Assumptions -> $Assumptions];
n0Value = FullSimplify[n0Bundle /. sampleFull, Assumptions -> $Assumptions];
u2Value = FullSimplify[-d2Value/d0Value, Assumptions -> $Assumptions];
u4Value = FullSimplify[(d2Value^2 - d0Value d4Value)/d0Value^2, Assumptions -> $Assumptions];
p0Value = FullSimplify[n0Value/d0Value, Assumptions -> $Assumptions];

Print["D0 = ", fmt[N[d0Value, 20]]];
Print["D2 = ", fmt[N[d2Value, 20]]];
Print["D4 = ", fmt[N[d4Value, 20]]];
Print["u2 = ", fmt[N[u2Value, 20]]];
Print["u4 = ", fmt[N[u4Value, 20]]];
Print["P0 = ", fmt[N[p0Value, 20]]];
expectClose[d0Value, num["24.2373099886223"], 10^-12, "M3 D0 compatibility value"];
expectClose[d2Value, num["-1.18562046858190"], 10^-12, "M3 D2 compatibility value"];
expectClose[d4Value, num["-0.173991572849491"], 10^-12, "M3 D4 compatibility value"];
expectClose[u2Value, num["0.0489171640391802"], 10^-12, "M3 u2 compatibility value"];
expectClose[u4Value, num["0.00957155575054425"], 10^-12, "M3 u4 compatibility value"];
expectClose[d4Value/d0Value, num["-0.00717866681290820"], 10^-12, "M3 D4/D0 compatibility value"];
expectClose[p0Value, num["0.00206979231806289"], 10^-12, "M3 P0 compatibility value"];
expectZero[u4Value - 4 u2Value^2, "M3 exact one-pole identity"];

subbanner["M4. Concrete even-gate coefficients on the branch"];

k1Coeff = FullSimplify[1/9 - u2Value, Assumptions -> $Assumptions];
hEvenCoeff = FullSimplify[d4Value/d0Value + (2/3) u2Value - 1/27, Assumptions -> $Assumptions];
Print["K1 coefficient = ", fmt[N[k1Coeff, 20]]];
Print["H_even coefficient = ", fmt[N[hEvenCoeff, 20]]];
expectClose[k1Coeff, num["0.0621939470719309"], 10^-12, "M4 K1 coefficient"];
expectClose[hEvenCoeff, num["-0.0116042611571584"], 10^-12, "M4 H_even coefficient"];

subbanner["M5/M6. Mixed-sector primitive compiler and corridor matrices"];

dressRules = {
  kStiff -> kStiff Exp[eps xk],
  mass -> mass Exp[eps xm],
  lamB -> lamB Exp[eps xlb],
  varpi -> varpi Exp[eps xv],
  lamU -> lamU Exp[eps xlu],
  lamW -> lamW Exp[eps xlw],
  lamR -> lamR Exp[eps xlr],
  omegaU -> omegaU Exp[eps xou],
  omegaW -> omegaW Exp[eps xow]
};

d01Expr = D[d0Bundle /. dressRules, eps] /. eps -> 0;
d21Expr = D[d2Bundle /. dressRules, eps] /. eps -> 0;
d41Expr = D[d4Bundle /. dressRules, eps] /. eps -> 0;
n01Expr = D[n0Bundle /. dressRules, eps] /. eps -> 0;

nonMixedZero = {xk -> 0, xm -> 0, xlb -> 0, xv -> 0};
mixedVars = {xlu, xlw, xlr, xou, xow};

d01Mixed = Expand[FullSimplify[d01Expr /. sampleFull /. nonMixedZero, Assumptions -> $Assumptions]];
d21Mixed = Expand[FullSimplify[d21Expr /. sampleFull /. nonMixedZero, Assumptions -> $Assumptions]];
d41Mixed = Expand[FullSimplify[d41Expr /. sampleFull /. nonMixedZero, Assumptions -> $Assumptions]];
n01Mixed = Expand[FullSimplify[n01Expr /. sampleFull /. nonMixedZero, Assumptions -> $Assumptions]];
xi1Mixed = Expand[FullSimplify[n01Mixed/n0Value - d01Mixed/d0Value, Assumptions -> $Assumptions]];

k1Mixed = Expand[FullSimplify[d21Mixed + d01Mixed/9, Assumptions -> $Assumptions]];
hEvenMixed = Expand[FullSimplify[d41Mixed - (2/3) d21Mixed - d01Mixed/27, Assumptions -> $Assumptions]];

evenMatrixExact = {coeffRow[k1Mixed, mixedVars], coeffRow[hEvenMixed, mixedVars]};
transferMatrixExact = {
  coeffRow[d01Mixed, mixedVars],
  coeffRow[d21Mixed, mixedVars],
  coeffRow[d41Mixed, mixedVars]
};
xiCoeffExact = coeffRow[xi1Mixed, mixedVars];

Print["Strict even-gate matrix = ", fmt[N[evenMatrixExact, 14]]];
expectZero[MatrixRank[evenMatrixExact] - 2, "M5 strict even-gate rank 2"];
evenNull = NullSpace[evenMatrixExact];
expectZero[Length[evenNull] - 3, "M5 strict even-gate nullity 3"];

sigmaEven = projectedNorm[xiCoeffExact, evenNull];
Print["sigma_even = ", fmt[N[sigmaEven, 20]]];
expectClose[sigmaEven, num["2.67386816837173"], 10^-11, "M5 sigma_even"];

subbanner["M6. Pure-transfer subcorridor"];

Print["Pure-transfer matrix = ", fmt[N[transferMatrixExact, 14]]];
expectZero[MatrixRank[transferMatrixExact] - 3, "M6 transfer matrix rank 3"];
transferNull = NullSpace[transferMatrixExact];
expectZero[Length[transferNull] - 2, "M6 transfer matrix nullity 2"];

Do[
  Module[{residuals},
    residuals = Map[
      FullSimplify[Cancel[Together[#]], Assumptions -> $Assumptions] &,
      transferMatrixExact.transferNull[[i]]
    ];
    Print["transfer null vector ", i, " constraint residuals = ", fmt[N[residuals, 20]]];
    Do[
      expectClose[residuals[[j]], 0, 10^-12, "M6c transfer null residual " <> ToString[i] <> "." <> ToString[j]],
      {j, Length[residuals]}
    ];
  ],
  {i, Length[transferNull]}
];

sigmaTransfer = projectedNorm[xiCoeffExact, transferNull];
Print["sigma_transfer = ", fmt[N[sigmaTransfer, 20]]];
expectClose[sigmaTransfer, num["2.31561904386057"], 10^-11, "M6 sigma_transfer"];

subbanner["M7. Transported ceiling budgets"];

baseBudgets = num /@ {
  "0.367930328492646",
  "0.737619063660757",
  "2.94889585703134",
  "4.63505472371892"
};
expectedEvenBudgets = num /@ {
  "0.137602269567650",
  "0.275862165676603",
  "1.10285760977778",
  "1.73346419189450"
};
expectedTransferBudgets = num /@ {
  "0.158890698998242",
  "0.318540765855427",
  "1.27348056877049",
  "2.00164821411704"
};

evenBudgets = baseBudgets/sigmaEven;
transferBudgets = baseBudgets/sigmaTransfer;
Print["even budgets = ", fmt[N[evenBudgets, 18]]];
Print["transfer budgets = ", fmt[N[transferBudgets, 18]]];

Do[
  expectClose[
    evenBudgets[[i]],
    expectedEvenBudgets[[i]],
    10^-11,
    "M7 even budget " <> ToString[i]
  ],
  {i, Length[baseBudgets]}
];

Do[
  expectClose[
    transferBudgets[[i]],
    expectedTransferBudgets[[i]],
    10^-11,
    "M7 transfer budget " <> ToString[i]
  ],
  {i, Length[baseBudgets]}
];

Print[""];
Print["All Stage 226 Mathematica checks passed."];
Exit[0];
