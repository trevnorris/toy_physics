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
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

normalizeExpr[expr_] := Module[{res},
  res = If[
    MatrixQ[expr],
    Map[FullSimplify[#, Assumptions -> $Assumptions] &, expr, {2}],
    If[
      VectorQ[expr],
      Map[FullSimplify[#, Assumptions -> $Assumptions] &, expr],
      FullSimplify[expr, Assumptions -> $Assumptions]
    ]
  ];
  res
];

allZeroQ[expr_] := If[
  ListQ[expr],
  And @@ Flatten[Map[TrueQ[# == 0] &, expr, {ArrayDepth[expr]}]],
  TrueQ[expr == 0]
];

prettyArray[arr_] := If[VectorQ[arr], MatrixForm[{arr}], MatrixForm[arr]];

expectZero[name_String, expr_] := Module[{res},
  res = normalizeExpr[expr];
  If[ListQ[res],
    Print[name, " ="];
    Print[prettyArray[res]];
    If[allZeroQ[res], pass[name], fail[name, res]],
    Print[name, " = ", fmt[res]];
    If[allZeroQ[res], pass[name], fail[name, res]]
  ];
];

eulerLagrange2D[L_, field_, tVar_, wVar_] := Module[
  {qt, qw},
  qt = D[field, tVar];
  qw = D[field, wVar];
  FullSimplify[
    D[D[L, qt], tVar] + D[D[L, qw], wVar] - D[L, field],
    Assumptions -> $Assumptions
  ]
];

banner["STAGE 001 — GEOMETRY LIFT AND LINEARIZED PDE SKELETON"];

subbanner["I. Real-harmonic bookkeeping"];

Clear[theta, phi];
$Assumptions = Element[{theta, phi}, Reals] && 0 < theta < Pi;

dOmega = Sin[theta];
y00 = 1/(2 Sqrt[Pi]);
y20 = Sqrt[5] (3 Cos[theta]^2 - 1)/(4 Sqrt[Pi]);
y21c = Sqrt[15] Sin[theta] Cos[theta] Cos[phi]/(2 Sqrt[Pi]);
y21s = Sqrt[15] Sin[theta] Cos[theta] Sin[phi]/(2 Sqrt[Pi]);
y22c = Sqrt[15] Sin[theta]^2 Cos[2 phi]/(4 Sqrt[Pi]);
y22s = Sqrt[15] Sin[theta]^2 Sin[2 phi]/(4 Sqrt[Pi]);

basis = <|
  "Y00" -> y00,
  "Y20" -> y20,
  "Y21c" -> y21c,
  "Y21s" -> y21s,
  "Y22c" -> y22c,
  "Y22s" -> y22s
|>;

sphereIntegral[expr_] := FullSimplify[
  Integrate[expr dOmega, {phi, 0, 2 Pi}, {theta, 0, Pi}],
  Assumptions -> $Assumptions
];

Do[
  avg = sphereIntegral[basis[name]];
  norm = sphereIntegral[basis[name]^2];
  Print[name, ": average = ", fmt[avg], ", norm = ", fmt[norm]],
  {name, Keys[basis]}
];

expectZero["average(Y20)", sphereIntegral[y20]];
expectZero["average(Y21c)", sphereIntegral[y21c]];
expectZero["average(Y21s)", sphereIntegral[y21s]];
expectZero["average(Y22c)", sphereIntegral[y22c]];
expectZero["average(Y22s)", sphereIntegral[y22s]];
expectZero["norm(Y00)-1", sphereIntegral[y00^2] - 1];
expectZero["norm(Y20)-1", sphereIntegral[y20^2] - 1];
expectZero["cross(Y00,Y20)", sphereIntegral[y00 y20]];

subbanner["I.2 Mouth-average extraction"];

Clear[q00, q20, q21c, q22c];
$Assumptions = Element[{theta, phi, q00, q20, q21c, q22c}, Reals] && 0 < theta < Pi;

eta = q00 y00 + q20 y20 + q21c y21c + q22c y22c;
mouthAvg = FullSimplify[sphereIntegral[eta]/(4 Pi), Assumptions -> $Assumptions];
expectZero["mouth average - q00/(2 sqrt(pi))", mouthAvg - q00/(2 Sqrt[Pi])];

subbanner["I.3 Spherical Laplacian eigenvalues"];

lapS2[expr_] := FullSimplify[
  (1/Sin[theta]) D[Sin[theta] D[expr, theta], theta] + D[expr, {phi, 2}]/Sin[theta]^2,
  Assumptions -> $Assumptions
];

expectZero["-Delta_S2 Y00", -lapS2[y00]];
expectZero["-Delta_S2 Y20 - 6 Y20", -lapS2[y20] - 6 y20];
expectZero["-Delta_S2 Y21c - 6 Y21c", -lapS2[y21c] - 6 y21c];
expectZero["-Delta_S2 Y21s - 6 Y21s", -lapS2[y21s] - 6 y21s];
expectZero["-Delta_S2 Y22c - 6 Y22c", -lapS2[y22c] - 6 y22c];
expectZero["-Delta_S2 Y22s - 6 Y22s", -lapS2[y22s] - 6 y22s];

subbanner["II. Level-set confinement linearization"];

Clear[sigma0, etaMode, ellc, eps];
$Assumptions = Element[{sigma0, etaMode, ellc, eps}, Reals] && ellc > 0;

expr = Vwall[(sigma0 - eps etaMode)/ellc];
firstVar = FullSimplify[D[expr, eps] /. eps -> 0, Assumptions -> $Assumptions];
targetVar = -etaMode D[Vwall[sigma0/ellc], sigma0];
Print["First variation = ", fmt[firstVar]];
expectZero["linearized confinement variation", firstVar - targetVar];

subbanner["III. Modal wall action"];

Clear[t, w, ell];
$Assumptions =
  Element[{t, w}, Reals] && Element[ell, Integers] && ell >= 0;

muEta = MuEta[w];
tW = TW[w];
tOmega = TOmega[w];
kEta = KEta[w];
g = G[w];
qField = q[t, w];
kEll = kEta + ell (ell + 1) tOmega;

ldens = (1/2) muEta D[qField, t]^2 - (1/2) tW D[qField, w]^2 - (1/2) kEll qField^2;
elDens = FullSimplify[-eulerLagrange2D[ldens, qField, t, w], Assumptions -> $Assumptions];
targetDens = -muEta D[qField, {t, 2}] + D[tW D[qField, w], w] - kEll qField;
expectZero["densitized Euler-Lagrange equation", elDens - targetDens];

lweighted = g ldens;
elWeighted = FullSimplify[-eulerLagrange2D[lweighted, qField, t, w], Assumptions -> $Assumptions];
targetWeighted = -g muEta D[qField, {t, 2}] + D[g tW D[qField, w], w] - g kEll qField;
expectZero["weighted Euler-Lagrange equation", elWeighted - targetWeighted];

expectZero["K_l at ell = 0", (kEll /. ell -> 0) - kEta];
expectZero["K_l at ell = 2", (kEll /. ell -> 2) - (kEta + 6 tOmega)];

subbanner["III.4 Sourced modal wall equation"];
Clear[Slm, fext];
$Assumptions = Element[{t, w}, Reals] && Element[ell, Integers] && ell >= 0;
sourceTotal = Slm[t, w] + fext[t, w];
ldensForced = ldens - qField*sourceTotal;
elForced = FullSimplify[-eulerLagrange2D[ldensForced, qField, t, w], Assumptions -> $Assumptions];
targetForced = targetDens - sourceTotal;
expectZero["sourced densitized Euler-Lagrange equation", elForced - targetForced];

subbanner["IV. Representative localized-Maxwell linearization"];
Clear[x, gaugeXi, mu0, zloc, axField, awField, jxField, jwField];
$Assumptions =
  Element[{x, w, gaugeXi, mu0}, Reals] && gaugeXi > 0 && mu0 > 0;
zloc = Z[w];
axField = Ax[x, w];
awField = Aw[x, w];
jxField = Jx[x, w];
jwField = Jw[x, w];
fwx = D[axField, w] - D[awField, x];
divA = D[axField, x] + D[awField, w];
lmax = (1/2) zloc fwx^2 - divA^2/(2 gaugeXi) + mu0 (jxField axField + jwField awField);
elAx = FullSimplify[
  D[D[lmax, D[axField, x]], x] + D[D[lmax, D[axField, w]], w] - D[lmax, axField],
  Assumptions -> $Assumptions
];
elAw = FullSimplify[
  D[D[lmax, D[awField, x]], x] + D[D[lmax, D[awField, w]], w] - D[lmax, awField],
  Assumptions -> $Assumptions
];
targetAx = D[zloc fwx, w] - D[divA, x]/gaugeXi - mu0 jxField;
targetAw = -D[zloc fwx, x] - D[divA, w]/gaugeXi - mu0 jwField;
expectZero["localized-Maxwell x-component", elAx - targetAx];
expectZero["localized-Maxwell w-component", elAw - targetAw];

Print[""];
Print["All Stage-001 symbolic checks passed."];
