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
  Print[name, " actual = ", fmt[N[a, 30]], ", residual = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

expectVectorClose[name_String, actual_List, expected_List, tol_: 10^-10] := Module[
  {diffs},
  diffs = MapThread[N[Abs[N[#1, 60] - num50[#2]], 30] &, {actual, expected}];
  Print[name, " residuals = ", fmt[diffs]];
  If[And @@ (TrueQ[# <= tol] & /@ diffs), pass[name], fail[name, diffs]];
];

positiveFrequencyRoots[poly_] := Module[{roots, yRoots, polyWork},
  polyWork = SetPrecision[N[poly, 70], 70];
  roots = y /. NSolve[polyWork == 0, y, WorkingPrecision -> 60];
  yRoots = Sort[
    N[
      Re /@ Select[
        roots,
        TrueQ[Abs[Im[N[#, 70]]] < 10^-35 && Re[N[#, 70]] > 0] &
      ],
      60
    ]
  ];
  N[Sqrt /@ yRoots, 50]
];

wallLikeRoots[roots_List] := Sort[TakeLargest[roots, 2]];

bisectPcompat[func_, left_, right_, target_, index_Integer, steps_Integer: 44] := Module[
  {a, b, fa, fb, mid, fm, root},
  a = N[left, 60];
  b = N[right, 60];
  fa = N[func[a][[index]] - target, 60];
  fb = N[func[b][[index]] - target, 60];
  If[TrueQ[fa == 0], Return[{a, func[a][[1]]}]];
  If[TrueQ[fb == 0], Return[{b, func[b][[1]]}]];
  If[!TrueQ[fa fb < 0], fail["bisection bracket straddles target", {left, right, fa, fb}]];
  Do[
    mid = N[(a + b)/2, 60];
    fm = N[func[mid][[index]] - target, 60];
    If[TrueQ[fa fm > 0],
      a = mid; fa = fm,
      b = mid
    ],
    {steps}
  ];
  root = N[(a + b)/2, 60];
  {root, func[root][[1]]}
];

Clear[
  s, len, y, omega, kap, lamB, lamU, lamW, lamR,
  omegaU, omegaW, varpi, kWall, mass, p0Target,
  aScale, cs, b0, b2, b4, z0, z2, z4, nZero
];

$Assumptions = (
  Element[
    {
      s, len, y, omega, kap, lamB, lamU, lamW, lamR,
      omegaU, omegaW, varpi, kWall, mass, p0Target,
      aScale, cs, b0, b2, b4, z0, z2, z4, nZero
    },
    Reals
  ]
  && len > 0 && kap > 0
  && lamB > 0 && lamU > 0 && lamW > 0 && lamR > 0
  && omegaU > 0 && omegaW > 0 && varpi > 0
  && mass > 0 && p0Target > 0
  && aScale > 0 && cs > 0 && nZero > 0
  && kWall - b0 - z0 != 0
  && b4 + z4 != 0
  && mass + b2 + z2 != 0
);

banner["STAGE 223 -- ISOTROPIC TARGET-SURFACE COMPATIBILITY"];

subbanner["M1. Primitive overlap constant"];

uConst = 1/Sqrt[len];
fProfile = Sqrt[2/len] Sin[Pi s/(2 len)];
kapIntegral = FullSimplify[Integrate[uConst fProfile, {s, 0, len}], Assumptions -> len > 0];
expectZero["M1 kappa integral", kapIntegral - 2 Sqrt[2]/Pi];
Print["kappa = ", fmt[kapIntegral]];

subbanner["M2. One-pole numerator identity"];

d0Iso = kWall - b0 - z0;
d2Iso = -(mass + b2 + z2);
d4Iso = -(b4 + z4);
u2Iso = -d2Iso/d0Iso;
u4Iso = (d2Iso^2 - d0Iso d4Iso)/d0Iso^2;
onePoleResidual = (
  u4Iso - 4 u2Iso^2
  - (d0Iso (b4 + z4) - 3 (mass + b2 + z2)^2)/d0Iso^2
);
expectZero["M2 u4 - 4 u2^2 numerator identity", onePoleResidual];

subbanner["M3. Compatibility surface from solved stiffness equality"];

kPoleIso = 3 (mass + b2 + z2)^2/(b4 + z4) + b0 + z0;
kNormIso = nZero/p0Target + b0 + z0;
solvedTarget = stripConditional[p0Target /. First[Solve[kNormIso == kPoleIso, p0Target]]];
boxedTarget = nZero (b4 + z4)/(3 (mass + b2 + z2)^2);

expectZero["M3 solved target minus boxed target", solvedTarget - boxedTarget];
expectZero[
  "M3 equivalent N0/P0 compatibility equation",
  nZero/boxedTarget - 3 (mass + b2 + z2)^2/(b4 + z4)
];

subbanner["M4. Primitive specialization"];

cCoef = kap lamB;
gU = lamU;
gW = kap lamW;
rMix = kap lamR;

delta0 = omegaU^2 omegaW^2 - rMix^2;
s2 = omegaU^2 + omegaW^2;
h2 = gU^2 + gW^2;
q0 = gU^2 omegaW^2 + 2 gU gW rMix + gW^2 omegaU^2;
pPrimitive = omegaU^2 gW + rMix gU;

b0Primitive = cCoef^2/varpi^2;
b2Primitive = cCoef^2/varpi^4;
b4Primitive = cCoef^2/varpi^6;
z0Primitive = q0/delta0;
z2Primitive = (q0 s2 - h2 delta0)/delta0^2;
z4Primitive = (q0 (s2^2 - delta0) - s2 h2 delta0)/delta0^3;
n0Primitive = pPrimitive^2/delta0^2;

pCompatPrimitive = FullSimplify[
  boxedTarget /. {
    nZero -> n0Primitive,
    b4 -> b4Primitive,
    b2 -> b2Primitive,
    z4 -> z4Primitive,
    z2 -> z2Primitive
  },
  Assumptions -> $Assumptions
];
primitiveClaim = FullSimplify[
  (pPrimitive^2/delta0^2) (cCoef^2/varpi^6 + z4Primitive)
    /(3 (mass + cCoef^2/varpi^4 + z2Primitive)^2),
  Assumptions -> $Assumptions
];
expectZero["M4 primitive boxed surface", pCompatPrimitive - primitiveClaim];

subbanner["M5. Quartic determinant/product form"];

deltaY = (omegaU^2 - y) (omegaW^2 - y) - rMix^2;
sourceY = gU^2 (omegaW^2 - y) + 2 gU gW rMix + gW^2 (omegaU^2 - y);
fy = (
  ((kWall - mass y) (varpi^2 - y) - cCoef^2) deltaY
  - (varpi^2 - y) sourceY
);
fyExpanded = Expand[fy];
fyCoefficients = CoefficientList[fyExpanded, y];

Print["F[y] degree = ", Exponent[fyExpanded, y]];
Print["F[y] coefficient count = ", Length[fyCoefficients]];
expectTrue["M5 determinant/product polynomial is quartic", Exponent[fyExpanded, y] == 4];
expectZero["M5 leading quartic coefficient is mass", Last[fyCoefficients] - mass];

subbanner["M6. Concrete compatibility-branch slice"];

sampleRulesFor[lw_] := {
  kap -> 2 Sqrt[2]/Pi,
  lamB -> 1/2,
  lamU -> 3/10,
  lamW -> lw,
  lamR -> 1/4,
  omegaU -> 1,
  omegaW -> 7/5,
  varpi -> 2,
  mass -> 1,
  aScale -> 1,
  cs -> 1
};
sampleRules = sampleRulesFor[2/5];

kCompatPrimitive = (
  3 (mass + b2Primitive + z2Primitive)^2/(b4Primitive + z4Primitive)
  + b0Primitive + z0Primitive
);
d0CompatPrimitive = kCompatPrimitive - b0Primitive - z0Primitive;

Print["P0_target_compat = ", fmt[N[pCompatPrimitive /. sampleRules, 30]]];
Print["K_compat = ", fmt[N[kCompatPrimitive /. sampleRules, 30]]];
Print["D0_compat = ", fmt[N[d0CompatPrimitive /. sampleRules, 30]]];

expectClose["M6 P0_target_compat", pCompatPrimitive /. sampleRules, "0.0020697923180628827", 10^-15];
expectClose["M6 K_compat", kCompatPrimitive /. sampleRules, "24.473754879290977", 10^-12];
expectClose["M6 D0_compat", d0CompatPrimitive /. sampleRules, "24.23730998862225", 10^-12];

subbanner["M7. Compatibility-branch pole census"];

fyCompatSample = fyExpanded /. sampleRules /. kWall -> (kCompatPrimitive /. sampleRules);
coupledRoots = positiveFrequencyRoots[fyCompatSample];
wallRoots = wallLikeRoots[coupledRoots];
internalRoots = Complement[coupledRoots, wallRoots];

Print["positive omega roots = ", fmt[N[coupledRoots, 18]]];
Print["internal-like roots = ", fmt[N[internalRoots, 18]]];
Print["wall-like roots = ", fmt[N[wallRoots, 18]]];

expectVectorClose[
  "M7 positive omega roots",
  coupledRoots,
  {"0.971575315129468", "1.41651290122561", "1.99753567893361", "4.94905432364313"},
  2 10^-12
];
expectTrue["M7 wall-like roots are the largest two roots", Min[wallRoots] > Max[internalRoots]];

subbanner["M8. Residue/linewidth values"];

nOmega = Together[
  ((omegaU^2 - omega^2) gW + rMix gU)^2
    /(((omegaU^2 - omega^2) (omegaW^2 - omega^2) - rMix^2)^2)
];
rqOmega = Together[27 cs^5/(aScale^5 omega^5 nOmega)];
rqValues = N[(rqOmega /. sampleRules /. omega -> #) & /@ coupledRoots, 50];

Print["R_Q,* values = ", fmt[N[rqValues, 18]]];
expectVectorClose[
  "M8 R_Q,* values",
  rqValues,
  {"0.159888393135835", "0.000806281535937178", "30.1999075602499", "36.1711864832695"},
  2 10^-10
];

subbanner["M9. Compatibility-family scan and survival windows"];

compatibilityData[lw_?NumericQ] := Module[
  {rules, pVal, kVal, poly, roots, walls, rqWalls},
  rules = sampleRulesFor[SetPrecision[lw, 60]];
  pVal = N[pCompatPrimitive /. rules, 60];
  kVal = N[kCompatPrimitive /. rules, 60];
  poly = N[fyExpanded /. rules /. kWall -> kVal, 70];
  roots = positiveFrequencyRoots[poly];
  walls = wallLikeRoots[roots];
  rqWalls = N[(rqOmega /. rules /. omega -> #) & /@ walls, 60];
  {pVal, kVal, rqWalls[[1]], rqWalls[[2]]}
];

scanLambdas = {2/5, 3/5, 4/5, 1};
scanRows = compatibilityData /@ scanLambdas;
Print["scan rows (lambda_W, P0, K, lower wall R_Q, upper wall R_Q):"];
Do[Print[fmt[N[Prepend[scanRows[[i]], scanLambdas[[i]]], 18]]], {i, Length[scanLambdas]}];

expectVectorClose[
  "M9 scan P0 values for lambda_W = 0.4..1.0",
  scanRows[[All, 1]],
  {"0.002069792318063", "0.004865681200486", "0.009169913681573", "0.014981190324091"},
  5 10^-13
];
expectVectorClose[
  "M9 scan K values for lambda_W = 0.4..1.0",
  scanRows[[All, 2]],
  {"24.4737548792910", "21.1544287401845", "19.0298300900561", "17.7824591822917"},
  2 10^-11
];
expectVectorClose[
  "M9 scan lower wall R_Q values for lambda_W = 0.4..1.0",
  scanRows[[All, 3]],
  {"30.1999075602499", "12.8348600273988", "7.06074242207991", "4.45922850098679"},
  3 10^-10
];
expectVectorClose[
  "M9 scan upper wall R_Q values for lambda_W = 0.4..1.0",
  scanRows[[All, 4]],
  {"36.1711864832695", "16.7575510327116", "9.69035785242054", "6.30111094469551"},
  3 10^-10
];

vKnown = num50["1.181909222592"];
deltaVreq = vKnown - 1/10;
threshold01 = N[2 deltaVreq (1 + (1/10)^2)/(1/10), 50];
threshold03 = N[2 deltaVreq (1 + (3/10)^2)/(3/10), 50];

expectClose["M9 threshold eta = 0.1", threshold01, "21.854566296358396", 10^-14];
expectClose["M9 threshold eta = 0.3", threshold03, "7.8618736841685335", 10^-15];

lowerStrict = bisectPcompat[compatibilityData, 2/5, 3/5, threshold01, 3];
upperStrict = bisectPcompat[compatibilityData, 2/5, 3/5, threshold01, 4];
lowerLoose = bisectPcompat[compatibilityData, 3/5, 4/5, threshold03, 3];
upperLoose = bisectPcompat[compatibilityData, 4/5, 1, threshold03, 4];

Print["10% lower wall ceiling = ", fmt[N[lowerStrict[[2]], 30]]];
Print["10% upper wall ceiling = ", fmt[N[upperStrict[[2]], 30]]];
Print["30% lower wall ceiling = ", fmt[N[lowerLoose[[2]], 30]]];
Print["30% upper wall ceiling = ", fmt[N[upperLoose[[2]], 30]]];

expectClose["M9 lower wall 10% P0 ceiling", lowerStrict[[2]], "0.00283133168555932", 5 10^-13];
expectClose["M9 upper wall 10% P0 ceiling", upperStrict[[2]], "0.00359651058968466", 5 10^-13];
expectClose["M9 lower wall 30% P0 ceiling", lowerLoose[[2]], "0.00817339430971383", 5 10^-13];
expectClose["M9 upper wall 30% P0 ceiling", upperLoose[[2]], "0.0116633929790174", 5 10^-13];

Print[""];
Print["All Stage 223 Mathematica checks M1-M9 completed successfully."];
