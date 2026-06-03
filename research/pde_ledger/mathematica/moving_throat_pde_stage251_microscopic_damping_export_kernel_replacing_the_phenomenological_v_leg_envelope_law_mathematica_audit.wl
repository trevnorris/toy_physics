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

stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[stripConditional[expr]]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[stripConditional[cond], Assumptions -> $Assumptions];
  res = stripConditional[res];
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 40] - N[target, 40]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 251 - MICROSCOPIC DAMPING/EXPORT KERNEL"];

Clear[
  omega, omu0, omw0, r0, eta0, gamma1, piV0, piVm, a, cs, beta0,
  sMinus, lamMinus, kProj, s, s0, sc, g3, g5, muEta, kappaV, t,
  v2, v3, ss, gg3, gg5, mm, xRoot, yRoot, gh3, gh5, gBook, eps1,
  rootVar, sInf
];

$Assumptions = (
  Element[
    {
      omega, omu0, omw0, r0, eta0, gamma1, piV0, piVm, a, cs, beta0,
      sMinus, lamMinus, kProj, s, s0, sc, g3, g5, muEta, kappaV,
      t, v2, v3, ss, gg3, gg5, mm, xRoot, yRoot, gh3, gh5,
      gBook, eps1, rootVar
    },
    Reals
  ] &&
  omu0 > 0 && omw0 > 0 && r0 > 0 && eta0 > 0 && gamma1 > 0 &&
  a > 0 && cs > 0 && beta0 > 0 && sMinus > 0 && lamMinus > 0 &&
  kProj > 0 && s > 0 && s0 > 0 && sc > 0 && g3 > 0 && g5 > 0 &&
  muEta > 0 && kappaV > 0 && gBook > 0
);

subbanner["M1. Cubic scalar outlet from a Series coefficient"];

a0[x_] := omu0^2 - x^2;
w0[x_] := omw0^2 - x^2;
gw0[x_] := eta0 x;
delta0 = omu0^2 omw0^2 - r0^2;
n0[x_] := (a0[x] gw0[x])^2/(a0[x] w0[x] - r0^2)^2;
n0Series = Normal[Series[n0[omega], {omega, 0, 4}]];
n0Coeff2 = FullSimplify[Coefficient[n0Series, omega, 2], Assumptions -> $Assumptions];
gamma30 = FullSimplify[gamma1 n0Coeff2, Assumptions -> $Assumptions];
gamma3Projected = FullSimplify[piV0^2 gamma30, Assumptions -> $Assumptions];

Print["N0 series = ", fmt[n0Series]];
Print["omega^2 coefficient = ", fmt[n0Coeff2]];
expectZero["M1 omega^2 coefficient", n0Coeff2 - eta0^2 omu0^4/delta0^2];
expectZero["M1 Gamma_{3,0}", gamma30 - gamma1 eta0^2 omu0^4/delta0^2];
expectZero["M1 projected Gamma3", gamma3Projected - piV0^2 gamma1 eta0^2 omu0^4/delta0^2];

subbanner["M2. Quintic projection structure"];

p0Minus = beta0 sMinus/lamMinus;
gamma5Minus = a^5 p0Minus/(27 cs^5);
gamma5Projected = FullSimplify[piVm^2 gamma5Minus, Assumptions -> $Assumptions];
kExpProjected = Expand[gamma3Projected s^3 + gamma5Projected s^5];

Print["P0_minus = ", fmt[p0Minus]];
Print["Gamma_{5,-} = ", fmt[gamma5Minus]];
Print["Gamma5 = ", fmt[gamma5Projected]];
expectZero[
  "M2 projection homogeneity",
  (gamma5Projected /. piVm -> kProj piVm) - kProj^2 gamma5Projected
];
expectZero["M2 projected kernel s^5 coefficient", Coefficient[kExpProjected, s, 5] - gamma5Projected];
expectZero["M2 projected kernel s^3 coefficient", Coefficient[kExpProjected, s, 3] - gamma3Projected];

subbanner["M3. Odd export kernel"];

kExp = Expand[g3 s^3 + g5 s^5];
sigmaExp = -I g3 omega^3 - I g5 omega^5;

Print["K_exp(s) = ", fmt[kExp]];
Print["Sigma_exp(omega) = ", fmt[sigmaExp]];
expectZero["M3 K s^3 coefficient", Coefficient[kExp, s, 3] - g3];
expectZero["M3 K s^5 coefficient", Coefficient[kExp, s, 5] - g5];
expectZero["M3 K s^2 coefficient absent", Coefficient[kExp, s, 2]];
expectZero["M3 K s^4 coefficient absent", Coefficient[kExp, s, 4]];
expectZero["M3 Sigma omega^3 coefficient", Coefficient[sigmaExp, omega, 3] + I g3];
expectZero["M3 Sigma omega^5 coefficient", Coefficient[sigmaExp, omega, 5] + I g5];
expectZero["M3 Sigma omega^2 coefficient absent", Coefficient[sigmaExp, omega, 2]];
expectZero["M3 Sigma omega^4 coefficient absent", Coefficient[sigmaExp, omega, 4]];

subbanner["M4. Schott identity"];

vfun = v[t];
fOdd = g3 D[vfun, {t, 3}] - g5 D[vfun, {t, 5}];
sOdd = (
  g3 D[vfun, t] D[vfun, {t, 2}]
  - g5 (D[vfun, t] D[vfun, {t, 4}] - D[vfun, {t, 2}] D[vfun, {t, 3}])
);
schottResidual = Expand[
  D[vfun, t] fOdd - D[sOdd, t] + g3 D[vfun, {t, 2}]^2 + g5 D[vfun, {t, 3}]^2
];

Print["F_odd = ", fmt[fOdd]];
Print["S_odd = ", fmt[sOdd]];
expectZero["M4 Schott residual", schottResidual];
expectTrue[
  "M4 P_exp nonnegative for real derivative data",
  Resolve[
    ForAll[{v2, v3, gg3, gg5},
      Implies[gg3 > 0 && gg5 > 0 && Element[{v2, v3}, Reals],
        gg3 v2^2 + gg5 v3^2 >= 0]
    ],
    Reals
  ]
];

subbanner["M5. Characteristic polynomial and positive-root certificate"];

fChar[x_] := g5 x^5 + g3 x^3 + muEta x^2 - kappaV;
fPrime = D[fChar[s], s];
diffQuot = (
  g5 (xRoot^4 + xRoot^3 yRoot + xRoot^2 yRoot^2 + xRoot yRoot^3 + yRoot^4)
  + g3 (xRoot^2 + xRoot yRoot + yRoot^2)
  + muEta (xRoot + yRoot)
);
m5DerivativePositive = Resolve[
  ForAll[{ss, gg3, gg5, mm},
    Implies[gg3 > 0 && gg5 > 0 && mm > 0 && ss > 0,
      5 gg5 ss^4 + 3 gg3 ss^2 + 2 mm ss > 0]
  ],
  Reals
];
m5QuotientPositive = Resolve[
  ForAll[{xRoot, yRoot, gg3, gg5, mm},
    Implies[gg3 > 0 && gg5 > 0 && mm > 0 && xRoot > 0 && yRoot > 0,
      gg5 (xRoot^4 + xRoot^3 yRoot + xRoot^2 yRoot^2 + xRoot yRoot^3 + yRoot^4)
      + gg3 (xRoot^2 + xRoot yRoot + yRoot^2)
      + mm (xRoot + yRoot) > 0]
  ],
  Reals
];

Print["F(s) = ", fmt[fChar[s]]];
Print["F'(s) = ", fmt[fPrime]];
expectZero["M5 derivative form", fPrime - (5 g5 s^4 + 3 g3 s^2 + 2 muEta s)];
expectTrue["M5 F(0) is negative", fChar[0] < 0];
expectTrue["M5 F(s) tends to +Infinity", Limit[fChar[sInf], sInf -> Infinity] == Infinity];
expectTrue["M5 derivative positive for all s>0", m5DerivativePositive];
expectZero["M5 difference quotient identity", fChar[xRoot] - fChar[yRoot] - (xRoot - yRoot) diffQuot];
expectTrue["M5 positive-root uniqueness quotient", m5QuotientPositive];

subbanner["M6. Small-kernel slowdown from a bookkeeping Series"];

weakPoly = Expand[
  fChar[s0 + gBook eps1] /. {kappaV -> muEta s0^2, g3 -> gBook g3, g5 -> gBook g5}
];
weakSeries = Normal[Series[weakPoly, {gBook, 0, 1}]];
weakBalance = Coefficient[weakSeries, gBook, 1];
eps1SolRaw = eps1 /. First[Solve[weakBalance == 0, eps1]];
eps1Sol = FullSimplify[stripConditional[eps1SolRaw], Assumptions -> $Assumptions];
rootShift = -(g3 s0^2 + g5 s0^4)/(2 muEta);

Print["weak balance = ", fmt[weakBalance]];
Print["delta s derived = ", fmt[eps1Sol]];
expectZero["M6 root-shift coefficient", eps1Sol - rootShift];

subbanner["M7. Event-safe surface and Session-IV benchmark"];

safeExpr = gh5 sc^5 + gh3 sc^3 + sc^2 - s0^2;
gh3SafeRaw = gh3 /. First[Solve[safeExpr == 0, gh3]];
gh3Safe = FullSimplify[stripConditional[gh3SafeRaw], Assumptions -> $Assumptions];

expectZero["M7 safe half-plane", gh3Safe + sc^2 gh5 - (s0^2 - sc^2)/sc^3];

tCrossNum = 182169718/100000000;
gammaCritNum = 694311167/100000000;
scNum = N[1/tCrossNum, 50];
s0Num = N[gammaCritNum, 50];
weightNum = N[scNum^2, 50];
rhsNum = N[(s0Num^2 - scNum^2)/scNum^3, 50];
g3HatSafeNum = rhsNum;
g5HatSafeNum = N[(s0Num^2 - scNum^2)/scNum^5, 50];

Print["s_c = ", fmt[scNum]];
Print["s_c^2 = ", fmt[weightNum]];
Print["safe rhs = ", fmt[rhsNum]];
Print["Gamma5hat safe = ", fmt[g5HatSafeNum]];
expectApprox["M7 weight sc^2", weightNum, ToExpression["0.3013336471`20"], 10^-10];
expectApprox["M7 safe rhs", rhsNum, ToExpression["289.61004918`20"], 10^-8];
expectApprox["M7 Gamma5hat safe", g5HatSafeNum, ToExpression["961.09429528`20"], 10^-8];

roots3 = rootVar /. NSolve[
  rootVar^2 + g3HatSafeNum rootVar^3 - s0Num^2 == 0,
  rootVar,
  Reals,
  WorkingPrecision -> 50
];
roots5 = rootVar /. NSolve[
  rootVar^2 + g5HatSafeNum rootVar^5 - s0Num^2 == 0,
  rootVar,
  Reals,
  WorkingPrecision -> 50
];
positiveRoots3 = Select[roots3, TrueQ[N[#] > 0] &];
positiveRoots5 = Select[roots5, TrueQ[N[#] > 0] &];

Print["positive cubic roots = ", fmt[positiveRoots3]];
Print["positive quintic roots = ", fmt[positiveRoots5]];
expectTrue["M7 cubic polynomial has one positive real root", Length[positiveRoots3] == 1];
expectTrue["M7 quintic polynomial has one positive real root", Length[positiveRoots5] == 1];
expectApprox["M7 cubic positive root equals sc", First[positiveRoots3], scNum, 10^-10];
expectApprox["M7 quintic positive root equals sc", First[positiveRoots5], scNum, 10^-10];

Print[""];
Print["STAGE 251 MATHEMATICA AUDIT PASSED"];
Exit[0];
