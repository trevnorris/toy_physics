ClearAll["Global`*"];
$HistoryLength = 0;
Needs["VariationalMethods`"];

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
  If[! MissingQ[detail], Print["  residual -> ", fmt[detail]]];
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

eulerLagrange1D[L_, field_, tVar_] := FullSimplify[
  D[D[L, D[field, tVar]], tVar] - D[L, field],
  Assumptions -> $Assumptions
];

lapS2[expr_, theta_, phi_] := FullSimplify[
  (1/Sin[theta]) D[Sin[theta] D[expr, theta], theta] + D[expr, {phi, 2}]/Sin[theta]^2,
  Assumptions -> $Assumptions
];

gradS2Inner[y1_, y2_, theta_, phi_] := FullSimplify[
  D[y1, theta] D[y2, theta] + D[y1, phi] D[y2, phi]/Sin[theta]^2,
  Assumptions -> $Assumptions
];

banner["STAGE 002 — BREATHING REDUCTION BACK TO THE OLD (a, L) CLOSURE"];

subbanner["I. Monopole normalization bridge"];

Clear[theta, phi, deltaA, q00];
$Assumptions = Element[{theta, phi, deltaA, q00}, Reals] && 0 < theta < Pi;

dOmega = Sin[theta];
y00 = FullSimplify[SphericalHarmonicY[0, 0, theta, phi]];

avgY00 = FullSimplify[
  Integrate[y00 dOmega, {phi, 0, 2 Pi}, {theta, 0, Pi}]/(4 Pi),
  Assumptions -> $Assumptions
];
normY00 = FullSimplify[
  Integrate[y00^2 dOmega, {phi, 0, 2 Pi}, {theta, 0, Pi}],
  Assumptions -> $Assumptions
];

expectZero["Y00 mouth average - 1/(2 sqrt(pi))", avgY00 - 1/(2 Sqrt[Pi])];
expectZero["norm(Y00) - 1", normY00 - 1];

mouthAvg = FullSimplify[q00 avgY00, Assumptions -> $Assumptions];
expectZero["mouth average from q00 Y00 - q00/(2 sqrt(pi))", mouthAvg - q00/(2 Sqrt[Pi])];
Print["Identifying the physical mouth-average shift with delta_a then gives q00 = 2 sqrt(pi) delta_a."];

angularPrefactor = FullSimplify[
  Integrate[(2 Sqrt[Pi] y00)^2 dOmega, {phi, 0, 2 Pi}, {theta, 0, Pi}],
  Assumptions -> $Assumptions
];
expectZero["angular prefactor from (2 sqrt(pi) Y00)^2 - 4 pi", angularPrefactor - 4 Pi];

subbanner["II. Axisymmetric two-mode reduction from the Stage 001 action"];

Clear[w, wL, wR, da, dL, dadt, dLdt];
$Assumptions = Element[{w, wL, wR, da, dL, dadt, dLdt}, Reals] && wL < wR;

muEta = muEtaF[w];
Tw = TwF[w];
K0 = K0F[w];
alphaA = alphaAF[w];
alphaL = alphaLF[w];

axisym = alphaA da + alphaL dL;
axisymT = alphaA dadt + alphaL dLdt;
axisymW = D[alphaA, w] da + D[alphaL, w] dL;

eta = 2 Sqrt[Pi] axisym y00;
etaT = 2 Sqrt[Pi] axisymT y00;
etaW = 2 Sqrt[Pi] axisymW y00;

lw = FullSimplify[
  (1/2) Integrate[(muEta etaT^2 - Tw etaW^2 - K0 eta^2) dOmega, {phi, 0, 2 Pi}, {theta, 0, Pi}],
  Assumptions -> $Assumptions
];

Q = {da, dL};
Qdot = {dadt, dLdt};

MaaExt = 2 Coefficient[lw, dadt, 2];
MLLExt = 2 Coefficient[lw, dLdt, 2];
MaLExt = Coefficient[Coefficient[lw, dadt], dLdt];
KaaExt = -2 Coefficient[lw, da, 2];
KLLExt = -2 Coefficient[lw, dL, 2];
KaLExt = -Coefficient[Coefficient[lw, da], dL];
MintegrandExtracted = {{MaaExt, MaLExt}, {MaLExt, MLLExt}};
KintegrandExtracted = {{KaaExt, KaLExt}, {KaLExt, KLLExt}};
MintegrandBoxed = 4 Pi {
  {muEta alphaA^2, muEta alphaA alphaL},
  {muEta alphaA alphaL, muEta alphaL^2}
};
KintegrandBoxed = 4 Pi {
  {Tw D[alphaA, w]^2 + K0 alphaA^2, Tw D[alphaA, w] D[alphaL, w] + K0 alphaA alphaL},
  {Tw D[alphaA, w] D[alphaL, w] + K0 alphaA alphaL, Tw D[alphaL, w]^2 + K0 alphaL^2}
};
expectZero["extracted M matches boxed M (4 Pi overlap form)", MintegrandExtracted - MintegrandBoxed];
expectZero["extracted K matches boxed K (4 Pi overlap form)", KintegrandExtracted - KintegrandBoxed];

Mintegrand = MintegrandExtracted;
Kintegrand = KintegrandExtracted;
lwTarget = FullSimplify[
  (1/2) (Qdot . MintegrandBoxed . Qdot - Q . KintegrandBoxed . Q),
  Assumptions -> $Assumptions
];

Mmat = {
  {Integrate[Mintegrand[[1, 1]], {w, wL, wR}], Integrate[Mintegrand[[1, 2]], {w, wL, wR}]},
  {Integrate[Mintegrand[[2, 1]], {w, wL, wR}], Integrate[Mintegrand[[2, 2]], {w, wL, wR}]}
};
Kmat = {
  {Integrate[Kintegrand[[1, 1]], {w, wL, wR}], Integrate[Kintegrand[[1, 2]], {w, wL, wR}]},
  {Integrate[Kintegrand[[2, 1]], {w, wL, wR}], Integrate[Kintegrand[[2, 2]], {w, wL, wR}]}
};
Lred = Integrate[lw, {w, wL, wR}];
LredTarget = Integrate[lwTarget, {w, wL, wR}];
expectZero["formal reduced Lagrangian - boxed matrix form", Lred - LredTarget];

Print["Recovered overlap matrices:"];
Print[MatrixForm[Mmat]];
Print[MatrixForm[Kmat]];

Clear[t];
$Assumptions = Element[{t, wL, wR}, Reals] && wL < wR;

qa = qA[t];
qLfun = qLfunF[t];
Maa = Mmat[[1, 1]]; MaL = Mmat[[1, 2]]; MLL = Mmat[[2, 2]];
Kaa = Kmat[[1, 1]]; KaL = Kmat[[1, 2]]; KLL = Kmat[[2, 2]];

lredTime = (1/2) (
  Maa D[qa, t]^2 + 2 MaL D[qa, t] D[qLfun, t] + MLL D[qLfun, t]^2 -
  Kaa qa^2 - 2 KaL qa qLfun - KLL qLfun^2
);

{elAEq, elLEq} = EulerEquations[lredTime, {qa, qLfun}, t];
elA = elAEq[[2]] - elAEq[[1]];
elL = elLEq[[2]] - elLEq[[1]];

expectZero[
  "Euler-Lagrange equation for q_a",
  elA - (Maa D[qa, {t, 2}] + MaL D[qLfun, {t, 2}] + Kaa qa + KaL qLfun)
];
expectZero[
  "Euler-Lagrange equation for q_L",
  elL - (MaL D[qa, {t, 2}] + MLL D[qLfun, {t, 2}] + KaL qa + KLL qLfun)
];

subbanner["III. Grouped real P2 orthonormality and degeneracy"];

Clear[theta, phi, w, wL, wR];
$Assumptions = Element[{theta, phi, w, wL, wR}, Reals] && 0 < theta < Pi && wL < wR;

dOmega = Sin[theta];
y20 = Sqrt[5]/(4 Sqrt[Pi]) (3 Cos[theta]^2 - 1);
y21c = Sqrt[15]/(2 Sqrt[Pi]) Sin[theta] Cos[theta] Cos[phi];
y21s = Sqrt[15]/(2 Sqrt[Pi]) Sin[theta] Cos[theta] Sin[phi];
y22c = Sqrt[15]/(4 Sqrt[Pi]) Sin[theta]^2 Cos[2 phi];
y22s = Sqrt[15]/(4 Sqrt[Pi]) Sin[theta]^2 Sin[2 phi];
basis = {y20, y21c, y21s, y22c, y22s};
basisNames = {"Y20", "Y21c", "Y21s", "Y22c", "Y22s"};

expectZero["phase shift: Y21s(theta,phi) - Y21c(theta,phi-Pi/2)", y21s - (y21c /. phi -> phi - Pi/2)];
expectZero["phase shift: Y22s(theta,phi) - Y22c(theta,phi-Pi/4)", y22s - (y22c /. phi -> phi - Pi/4)];

Do[
  normI = FullSimplify[
    Integrate[basis[[i]]^2 dOmega, {phi, 0, 2 Pi}, {theta, 0, Pi}],
    Assumptions -> $Assumptions
  ];
  angI = FullSimplify[
    Integrate[gradS2Inner[basis[[i]], basis[[i]], theta, phi] dOmega, {phi, 0, 2 Pi}, {theta, 0, Pi}],
    Assumptions -> $Assumptions
  ];
  expectZero[StringJoin["norm(", basisNames[[i]], ") - 1"], normI - 1];
  expectZero[StringJoin["angular energy(", basisNames[[i]], ") - 6"], angI - 6];
  ,
  {i, 1, Length[basis]}
];

normMatrix5 = Table[
  FullSimplify[
    Integrate[basis[[i]] basis[[j]] dOmega, {phi, 0, 2 Pi}, {theta, 0, Pi}],
    Assumptions -> $Assumptions
  ],
  {i, 1, Length[basis]}, {j, 1, Length[basis]}
];
gradMatrix5 = Table[
  FullSimplify[
    Integrate[gradS2Inner[basis[[i]], basis[[j]], theta, phi] dOmega, {phi, 0, 2 Pi}, {theta, 0, Pi}],
    Assumptions -> $Assumptions
  ],
  {i, 1, Length[basis]}, {j, 1, Length[basis]}
];
expectZero["real P2 norm matrix - IdentityMatrix[5]", normMatrix5 - IdentityMatrix[5]];
expectZero["real P2 angular stiffness matrix - 6 IdentityMatrix[5]", gradMatrix5 - 6 IdentityMatrix[5]];

Do[
  expectZero[
    StringJoin["-Delta_S2 ", basisNames[[i]], " - 6 ", basisNames[[i]]],
    -lapS2[basis[[i]], theta, phi] - 6 basis[[i]]
  ],
  {i, 1, Length[basis]}
];

muEta = muEtaF[w];
Tw = TwF[w];
TOmega = TOmegaF[w];
KEta = KEtaF[w];
beta2 = beta2F[w];

qvec = {q20, q21c, q21s, q22c, q22s};
qdotvec = {q20d, q21cd, q21sd, q22cd, q22sd};

Do[
  yi = basis[[i]];
  qi = qvec[[i]];
  qdi = qdotvec[[i]];
  lwP2i = FullSimplify[
    (1/2) Integrate[
      (
        muEta beta2^2 qdi^2 yi^2 -
        Tw D[beta2, w]^2 qi^2 yi^2 -
        KEta beta2^2 qi^2 yi^2 -
        TOmega beta2^2 qi^2 gradS2Inner[yi, yi, theta, phi]
      ) dOmega,
      {phi, 0, 2 Pi},
      {theta, 0, Pi}
    ],
    Assumptions -> $Assumptions
  ];
  lwP2Target = FullSimplify[
    (1/2) (
      muEta beta2^2 qdi^2 -
      (Tw D[beta2, w]^2 + (KEta + 6 TOmega) beta2^2) qi^2
    ),
    Assumptions -> $Assumptions
  ];
  expectZero[
    StringJoin["single-component reduced density for ", basisNames[[i]]],
    lwP2i - lwP2Target
  ],
  {i, 1, Length[basis]}
];

Clear[t];
$Assumptions = Element[{t, wL, wR}, Reals] && wL < wR;

q2 = q2fun[t];
M2 = Integrate[muEtaF[w] beta2F[w]^2, {w, wL, wR}];
K2 = Integrate[TwF[w] D[beta2F[w], w]^2 + (KEtaF[w] + 6 TOmegaF[w]) beta2F[w]^2, {w, wL, wR}];
l2 = (1/2) (M2 D[q2, t]^2 - K2 q2^2);
el2 = eulerLagrange1D[l2, q2, t];
expectZero["single-component P2 Euler-Lagrange equation", el2 - (M2 D[q2, {t, 2}] + K2 q2)];

banner["STAGE 002 SUMMARY"];
Print["Verified with Mathematica:"];
Print["  • the Y00 normalization bridge behind q00 = 2 sqrt(pi) delta a;"];
Print["  • insertion of the two-mode ansatz into the Stage 001 wall action, producing"];
Print["    the boxed overlap-integral formulas for M_AB and K_AB;"];
Print["  • the conservative breathing-reduction matrix system M_AB Qdd^B + K_AB Q^B"];
Print["    using those overlap integrals as coefficients;"];
Print["  • real P2 orthonormality, common angular stiffness 6, and the resulting"];
Print["    isotropic grouped-real P2 degeneracy before coupling."];

Exit[0];
