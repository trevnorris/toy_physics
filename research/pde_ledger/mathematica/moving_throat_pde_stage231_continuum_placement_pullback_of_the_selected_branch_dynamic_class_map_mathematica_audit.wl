ClearAll["Global`*"];
$HistoryLength = 0;
$MaxExtraPrecision = 1000;

fmt[expr_] := ToString[InputForm[expr]];
pass[name_String] := Print["PASS: ", name];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  detail -> ", fmt[detail]]];
  Exit[1];
);

stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

cleanScalar[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanScalar[expr];
  Print[name, " residual = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = stripConditional[FullSimplify[cond, Assumptions -> $Assumptions]];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

expectNumericTrue[name_String, cond_] := (
  Print[name, " = ", fmt[cond]];
  If[TrueQ[cond], pass[name], fail[name, cond]];
);

expectClose[name_String, actual_, expected_, tol_] := Module[{diff},
  diff = N[Abs[actual - expected], 50];
  Print[name, " actual = ", fmt[N[actual, 18]], " expected = ", fmt[N[expected, 18]], " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

expectInfinityPole[name_String, actual_] := Module[{inv},
  inv = FullSimplify[1/actual, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[actual], "; reciprocal = ", fmt[inv]];
  If[TrueQ[inv === 0], pass[name], fail[name, actual]];
];

thresholdXi[deltaValue_, capValue_] := Module[{onset, poly, sols, roots},
  onset = N[rnd /. {x -> 0, d -> deltaValue}, 80];
  If[TrueQ[onset <= N[capValue, 80] + 10^-30],
    0,
    poly = Together[capValue (9 deltaValue + 11 x) (9 deltaValue^2 + 18 deltaValue x + 11 x^2) - 72 deltaValue^2 (1 - x)];
    sols = x /. NSolve[poly == 0, x, Reals, WorkingPrecision -> 80];
    roots = Select[N[sols, 50], TrueQ[0 < # < 1] &];
    If[Length[roots] == 1, First[roots], fail["M5 unique threshold root", roots]]
  ]
];

thresholdR[deltaValue_, capValue_] := Module[{xc},
  xc = thresholdXi[deltaValue, capValue];
  N[f /. {d -> deltaValue, x -> xc}, 50]
];

Print["=== Stage 231 Mathematica audit: continuum pullback of selected-branch dynamic-class map ==="];

$Assumptions = Element[{x, d, c}, Reals] && 0 <= x < 1 && d > 0 && c > 0;

f = (9 d + 11 x)^4/(81 (1 - x) (9 d^2 + 18 d x + 11 x^2)^2);
g = 9 x (x + d)/(9 d + 11 x);
rnd = 72 d^2 (1 - x)/((9 d + 11 x) (9 d^2 + 18 d x + 11 x^2));

dF = Factor[D[f, x]];
dG = Factor[D[g, x]];
dRnd = Factor[D[rnd, x]];
dFNumerator = Factor[Numerator[Together[dF]]];
dFPolyFactor = Factor[dFNumerator/(9 d + 11 x)^3];
dFPublishedPoly = 81 d^3 + 189 d^2 x + 72 d^2 + 297 d x^2 + 121 x^3;
dGPublished = 9 (9 d^2 + 18 d x + 11 x^2)/(9 d + 11 x)^2;

Print["M1 Mathematica Factor[D[F,x]] = ", fmt[dF]];
Print["M1 Mathematica D[F,x] numerator polynomial factor = ", fmt[dFPolyFactor]];
Print["M1 Mathematica Factor[D[G,x]] = ", fmt[dG]];
Print["M1 Mathematica Factor[D[RND,x]] = ", fmt[dRnd]];
expectZero["M1 D[G,x] closed form", dG - dGPublished];
expectZero["M1 D[F,x] numerator polynomial factor", dFPolyFactor - dFPublishedPoly];
expectTrue["M1 D[F,x] > 0 on stable strip", Resolve[ForAll[{x, d}, Implies[0 <= x < 1 && d > 0, dF > 0]], Reals]];
expectTrue["M1 D[G,x] > 0 on stable strip", Resolve[ForAll[{x, d}, Implies[0 <= x < 1 && d > 0, dG > 0]], Reals]];
expectTrue["M1 D[RND,x] < 0 on stable strip", Resolve[ForAll[{x, d}, Implies[0 <= x < 1 && d > 0, dRnd < 0]], Reals]];

fOnset = FullSimplify[f /. x -> 0, Assumptions -> d > 0];
gOnset = FullSimplify[g /. x -> 0, Assumptions -> d > 0];
rndOnset = FullSimplify[rnd /. x -> 0, Assumptions -> d > 0];
fSoftLimit = Block[{$Assumptions = Element[d, Reals] && d > 0}, Limit[f, x -> 1, Direction -> "FromBelow"]];

expectZero["M2 F(0,d) - 1", fOnset - 1];
expectZero["M2 G(0,d)", gOnset];
expectZero["M2 RND(0,d) - 8/(9 d)", rndOnset - 8/(9 d)];
expectInfinityPole["M2 Limit F as x -> 1 from below", fSoftLimit];

pullbackSlope = Together[dRnd/dF];
probeGrid = {{1/100, 1/5}, {1/20, 1/4}, {1/10, 1/2}, {1/5, 3/4}, {2/5, 5/4}, {4/5, 2}};

Do[
  Module[{xv, dv, slopeVal},
    {xv, dv} = pt;
    slopeVal = N[pullbackSlope /. {x -> xv, d -> dv}, 50];
    Print["M3 dRphys/dRtarget at {x,d} = ", fmt[pt], " -> ", fmt[slopeVal]];
    expectNumericTrue["M3 sampled pullback slope is negative", slopeVal < 0];
  ],
  {pt, probeGrid}
];

pc = Expand[c (9 d + 11 x) (9 d^2 + 18 d x + 11 x^2) - 72 d^2 (1 - x)];
dPc = Factor[D[pc, x]];
dPcPublished = 3 (87 c d^2 + 198 c d x + 121 c x^2 + 24 d^2);
pcOnset = Factor[pc /. x -> 0];
deltaC = 8/(9 c);

Print["M4 P_c = ", fmt[pc]];
Print["M4 Factor[D[P_c,x]] = ", fmt[dPc]];
expectZero["M4 D[P_c,x] closed form", dPc - dPcPublished];
expectTrue["M4 D[P_c,x] > 0 on stable strip", Resolve[ForAll[{x, d, c}, Implies[0 <= x < 1 && d > 0 && c > 0, dPc > 0]], Reals]];
expectZero["M4 P_c(0,d) closed form", pcOnset - 9 d^2 (9 c d - 8)];
expectZero["M4 onset law P_c(0,delta_c)", pcOnset /. d -> deltaC];

rStar = (411024574532864/10^15)/(334368725711457/10^15);
rFlip080 = thresholdR[4/5, rStar];
rDen100 = thresholdR[1, 1];

thresholdRows = {
  {1/4, 1330868539/10^9, 1393832566/10^9},
  {1/2, 1139956630/10^9, 1221087062/10^9},
  {3/4, 1, 1071471867/10^9}
};

Do[
  Module[{dv, expectedFlip, expectedDen, xFlip, rFlip, xDen, rDen},
    {dv, expectedFlip, expectedDen} = row;
    xFlip = thresholdXi[dv, rStar];
    rFlip = N[f /. {x -> xFlip, d -> dv}, 50];
    xDen = thresholdXi[dv, 1];
    rDen = N[f /. {x -> xDen, d -> dv}, 50];
    Print[
      "M5 delta = ", fmt[N[dv, 18]],
      " x_flip = ", fmt[N[xFlip, 18]],
      " R_flip = ", fmt[N[rFlip, 18]],
      " x_den = ", fmt[N[xDen, 18]],
      " R_den = ", fmt[N[rDen, 18]]
    ];
    expectClose["M5 R_flip sample", rFlip, expectedFlip, 10^-9];
    expectClose["M5 R_den sample", rDen, expectedDen, 10^-9];
    expectNumericTrue["M5 R_flip <= R_den", N[rFlip, 50] <= N[rDen, 50] + 10^-30];
  ],
  {row, thresholdRows}
];

expectClose["M5 R_flip(0.80) collapse", rFlip080, 1, 10^-14];
expectClose["M5 R_den(1.00) collapse", rDen100, 1, 10^-14];

$Assumptions = Element[{epsEta, epsW, rho, zW, delta0, lam}, Reals] && epsEta != 1 && epsW != 1 && zW != 0 && 1 + rho != 0;
deltaPlacement = delta0/(1 - epsEta);
mMix = 8 zW (1 + rho)^2/(Pi^2 (1 - epsEta) (1 - epsW));
rTargetPlacement = lam (1 - epsEta) (1 - epsW)^2/(zW (1 + rho)^2);
productLawResidual = rTargetPlacement mMix - 8 lam (1 - epsW)/Pi^2;

Print["M6 delta placement = ", fmt[deltaPlacement]];
Print["M6 M_mix = ", fmt[mMix]];
Print["M6 R_target = ", fmt[rTargetPlacement]];
expectZero["M6 continuum placement product law", productLawResidual];

$Assumptions = True;
bDynBothInf = 967282389363822/10^15;
bDynNonemptyInf = 990581810705233/10^15;
bStatBoth = 367930328492646/10^15;
bStatNonempty = 737619063660757/10^15;

expectTrue["M7 dynamic both ceiling exceeds static both budget", bDynBothInf > bStatBoth];
expectTrue["M7 dynamic nonempty ceiling exceeds static nonempty budget", bDynNonemptyInf > bStatNonempty];

Do[
  Module[{xv, dv, sample, onset},
    {xv, dv} = pt;
    sample = N[rnd /. {x -> xv, d -> dv}, 50];
    onset = N[rnd /. {x -> 0, d -> dv}, 50];
    Print["M7 RND subset sample {x,d} = ", fmt[pt], " sample = ", fmt[sample], " onset = ", fmt[onset]];
    expectNumericTrue["M7 sampled RND lies in [0,RND(0,d)]", 0 <= sample <= onset];
  ],
  {pt, probeGrid}
];

Print["All Stage 231 Mathematica audits passed."];
Exit[0];
