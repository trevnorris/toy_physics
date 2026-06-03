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
num50[s_String] := ToExpression[s <> "`50"];
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

expectClose[name_String, actual_, expected_String, tol_: 10^-12] := Module[
  {a, e, diff},
  a = N[actual, 60];
  e = num50[expected];
  diff = N[Abs[a - e], 30];
  Print[name, " actual = ", fmt[N[a, 18]], ", residual = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

expectVectorClose[name_String, actual_List, expected_List, tol_: 10^-10] := Module[
  {diffs},
  diffs = MapThread[N[Abs[N[#1, 60] - num50[#2]], 30] &, {actual, expected}];
  Print[name, " residuals = ", fmt[diffs]];
  If[And @@ (TrueQ[# <= tol] & /@ diffs), pass[name], fail[name, diffs]];
];

positiveFrequencyRoots[poly_] := Module[{roots, yRoots},
  roots = y /. NSolve[poly == 0, y, WorkingPrecision -> 50];
  yRoots = Sort[
    N[
      Re /@ Select[
        roots,
        TrueQ[Abs[Im[N[#, 60]]] < 10^-35 && Re[N[#, 60]] > 0] &
      ],
      50
    ]
  ];
  N[Sqrt /@ yRoots, 40]
];

strictIncreasingQ[vals_List] := And @@ (TrueQ[#2 > #1] & @@@ Partition[vals, 2, 1]);
strictDecreasingQ[vals_List] := And @@ (TrueQ[#2 < #1] & @@@ Partition[vals, 2, 1]);

Clear[
  s, len, y, omega, kap, lamB, lamU, lamW, lamR,
  omegaU, omegaW, varpi, kWall, mass, aScale, cs,
  d0p, omegaStar, nStar, eta, x, deltaVreq, rqRatio, rDet
];

$Assumptions = (
  Element[
    {
      s, len, y, omega, kap, lamB, lamU, lamW, lamR,
      omegaU, omegaW, varpi, kWall, mass, aScale, cs,
      d0p, omegaStar, nStar, eta, x, deltaVreq, rqRatio, rDet
    },
    Reals
  ]
  && len > 0 && kap > 0
  && lamB > 0 && lamU > 0 && lamW > 0 && lamR > 0
  && omegaU > 0 && omegaW > 0 && varpi > 0
  && kWall > 0 && mass > 0
  && aScale > 0 && cs > 0
  && d0p > 0 && omegaStar > 0 && nStar > 0
  && eta > 0 && x > 0 && deltaVreq > 0 && rqRatio > 0
);

banner["STAGE 222 -- FINITE-THROAT PRIMITIVE BRANCH MATHEMATICA AUDIT"];

subbanner["M1. Primitive finite-throat overlap"];

u0 = 1/Sqrt[len];
f0 = Sqrt[2/len] Sin[Pi s/(2 len)];
kapIntegral = FullSimplify[Integrate[u0 f0, {s, 0, len}], Assumptions -> len > 0];
expectZero["M1 overlap constant", kapIntegral - 2 Sqrt[2]/Pi];
Print["kappa = ", fmt[kapIntegral]];

subbanner["M2. Quartic pole polynomial from native denominator"];

cCoef = kap lamB;
gU = lamU;
gW = kap lamW;
rMix = kap lamR;

kbExpr = kWall - mass omega^2 - cCoef^2/(varpi^2 - omega^2);
uBranch = omegaU^2 - omega^2;
wBranch = omegaW^2 - omega^2;
deltaExpr = uBranch wBranch - rMix^2;
qExpr = gU^2 wBranch + 2 gU gW rMix + gW^2 uBranch;
dBranch = Together[kbExpr - qExpr/deltaExpr];

deltaY = (omegaU^2 - y) (omegaW^2 - y) - rMix^2;
qY = gU^2 (omegaW^2 - y) + 2 gU gW rMix + gW^2 (omegaU^2 - y);
fProductY = (
  ((kWall - mass y) (varpi^2 - y) - cCoef^2) deltaY
  - (varpi^2 - y) qY
);

nativeNumerator = Numerator[Together[dBranch]];
quarticFromD = Collect[
  Expand[PowerExpand[nativeNumerator /. omega -> Sqrt[y]]],
  y,
  FullSimplify
];
dBranchY = Together[
  Collect[Expand[PowerExpand[dBranch /. omega -> Sqrt[y]]], y, FullSimplify]
];
quarticDegree = Exponent[quarticFromD, y];
quarticCoeffs = CoefficientList[quarticFromD, y];

Print["native quartic degree = ", quarticDegree];
Print["native quartic coefficients = ", fmt[quarticCoeffs]];

expectZero[
  "M2 numerator-over-denominator reconstructs Together[D]",
  quarticFromD/((varpi^2 - y) deltaY) - dBranchY
];
expectZero["M2 native numerator equals product F[y]", quarticFromD - fProductY];
expectZero[
  "M2 D[omega] equals F[omega^2]/denominator",
  dBranch - (fProductY /. y -> omega^2)/((varpi^2 - omega^2) deltaExpr)
];
expectTrue["M2 F[y] has degree four", quarticDegree == 4];
expectTrue["M2 product-form F[y] has degree four", Exponent[fProductY, y] == 4];

subbanner["M3. Residue/linewidth and low-loss survival algebra"];

gamma5 = aScale^5/(27 cs^5);
aqqStar = 1/d0p;
gammaStar = gamma5 omegaStar^5 nStar/d0p;
rqDerived = FullSimplify[aqqStar/gammaStar, Assumptions -> $Assumptions];
rqExpected = 27 cs^5/(aScale^5 omegaStar^5 nStar);
expectZero["M3 residue/linewidth cancellation", rqDerived - rqExpected];

lineShapeAtBoundary = FullSimplify[
  (1/2) rqRatio (rDet/(1 + rDet^2))/x^6 /. rDet -> eta,
  Assumptions -> $Assumptions
];
stagePeak = (1/2) rqRatio eta/(1 + eta^2)/x^6;
thresholdFromPeak = stripConditional[rqRatio /. First[Solve[stagePeak == deltaVreq, rqRatio]]];
thresholdExpected = 2 deltaVreq (1 + eta^2) x^6/eta;

expectZero["M3 low-loss peak form", lineShapeAtBoundary - stagePeak];
expectZero["M3 survival threshold form", thresholdFromPeak - thresholdExpected];
Print["R_Q,* = ", fmt[rqDerived]];
Print["|Re V_Q|_max = ", fmt[stagePeak]];
Print["required R_Q,* = ", fmt[thresholdExpected]];

subbanner["M4. Static sample-slice data"];

sampleRulesFor[lw_] := {
  kap -> 2 Sqrt[2]/Pi,
  lamB -> 1/2,
  lamU -> 3/10,
  lamW -> lw,
  lamR -> 1/4,
  omegaU -> 1,
  omegaW -> 7/5,
  varpi -> 2,
  kWall -> 3,
  mass -> 1,
  aScale -> 1,
  cs -> 1
};
sampleRules = sampleRulesFor[2/5];

pPrimitive = uBranch gW + rMix gU;
nFun = Together[pPrimitive^2/deltaExpr^2];
rqExpr = Together[27 cs^5/(aScale^5 omega^5 nFun)];

delta0 = deltaExpr /. omega -> 0;
d0 = dBranch /. omega -> 0;
n0 = nFun /. omega -> 0;
p0 = n0/d0;

expectClose["M4 Delta0", delta0 /. sampleRules, "1.9093394081788311", 10^-12];
expectClose["M4 D0", d0 /. sampleRules, "2.7635551093312736", 10^-12];
expectClose["M4 N0", n0 /. sampleRules, "0.05016619802495911", 10^-12];
expectClose["M4 P0", p0 /. sampleRules, "0.018152776420332848", 10^-12];

subbanner["M5. Sample-slice pole census and pure-Q residues"];

fSampleY = fProductY /. sampleRules;
wallPolyY = ((kWall - mass y) (varpi^2 - y) - cCoef^2) /. sampleRules;
internalPolyY = deltaY /. sampleRules;

wallRoots = positiveFrequencyRoots[wallPolyY];
internalRoots = positiveFrequencyRoots[internalPolyY];
coupledRoots = positiveFrequencyRoots[fSampleY];
rqValues = N[(rqExpr /. sampleRules /. omega -> #) & /@ coupledRoots, 40];

Print["uncoupled wall/BdG roots = ", fmt[N[wallRoots, 18]]];
Print["uncoupled internal U/W roots = ", fmt[N[internalRoots, 18]]];
Print["coupled omega roots = ", fmt[N[coupledRoots, 18]]];
Print["R_Q,* values = ", fmt[N[rqValues, 18]]];

expectVectorClose[
  "M5 uncoupled wall/BdG roots",
  wallRoots,
  {"1.6814318259147836", "2.0427400751933362"},
  10^-12
];
expectVectorClose[
  "M5 uncoupled internal U/W roots",
  internalRoots,
  {"0.9746017237463136", "1.417798109771174"},
  10^-12
];
expectVectorClose[
  "M5 coupled roots",
  coupledRoots,
  {"0.9382727417467537", "1.3914108765380409", "1.7204537104800286", "2.045399487836587"},
  10^-10
];
expectVectorClose[
  "M5 pure-Q residue/linewidth values",
  rqValues,
  {"18.7069287828307", "0.380740659074003", "16.0250330226177", "32.0025481088465"},
  10^-10
];

subbanner["M6. Low-loss thresholds and lambda_W scan"];

vKnown = num50["1.181909222592"];
epsilon = 1/10;
deltaVreqNum = vKnown - epsilon;
thresholdEta01 = N[2 deltaVreqNum (1 + (1/10)^2)/(1/10), 50];
thresholdEta03 = N[2 deltaVreqNum (1 + (3/10)^2)/(3/10), 50];

expectClose["M6 threshold eta=0.1", thresholdEta01, "21.854566296358396", 10^-12];
expectClose["M6 threshold eta=0.3", thresholdEta03, "7.8618736841685335", 10^-12];

scanRow[lw_] := Module[{rules, roots, upperRoot, rqUpper},
  rules = sampleRulesFor[lw];
  roots = positiveFrequencyRoots[fProductY /. rules];
  upperRoot = Max[roots];
  rqUpper = N[rqExpr /. rules /. omega -> upperRoot, 50];
  N[{lw, p0 /. rules, d0 /. rules, upperRoot, rqUpper}, 40]
];

lambdaWScan = {1/5, 2/5, 3/5, 4/5, 1};
scanRows = scanRow /@ lambdaWScan;

Print["lambda_W scan rows {lambda_W, P0, D0, upper omega_*, upper R_Q,*}:"];
Scan[Print[fmt[N[#, 18]]] &, scanRows];
Print["lambda_W=0.2 upper-wall R_Q,* = ", fmt[N[scanRows[[1, 5]], 18]]];

expectTrue["M6 P0 increases with lambda_W", strictIncreasingQ[scanRows[[All, 2]]]];
expectTrue["M6 upper-wall R_Q,* decreases with lambda_W", strictDecreasingQ[scanRows[[All, 5]]]];

expectedScanRows = {
  {"0.2", "0.005947405317693074", "2.8272344215844973", "2.04402272302752", "145.4838586578641"},
  {"0.4", "0.018152776420332847", "2.7635551093312736", "2.045399487836587", "32.00254810884649"},
  {"0.6", "0.03800016314040507", "2.665913497209664", "2.04793277506821", "13.688535635680836"},
  {"0.8", "0.06717078268091271", "2.5343095852196678", "2.0519066888921147", "7.580971267465815"},
  {"1.0", "0.10847330811047598", "2.368743373361286", "2.0577833903550973", "4.827389255647022"}
};

Do[
  expectVectorClose[
    "M6 lambda_W scan row " <> ToString[i],
    scanRows[[i]],
    expectedScanRows[[i]],
    10^-10
  ],
  {i, Length[scanRows]}
];

Print[""];
Print["All Stage 222 Mathematica checks passed."];
