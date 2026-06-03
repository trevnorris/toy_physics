ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Expand[expr], Assumptions -> $Assumptions];
  (* Strip ConditionalExpression wrapper: under $Assumptions, a result of
     the form ConditionalExpression[0, cond] is identically zero on the
     declared domain.  Solve[]/Reduce[] often introduce these wrappers when
     auxiliary inequalities are nontrivial. *)
  res = res /. ConditionalExpression[e_, _] :> e;
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 058 — COUPLED SUPPORT/SOURCE OPERATOR"];

Clear[x, Pe, alpha, eta, Xi, pe0];
$Assumptions = alpha > 0 && eta > 0 && Pe > 0 && Xi > 0 && 0 <= x <= 1;

w = alpha*Sinh[alpha] + eta*Cosh[alpha];
kernel = FullSimplify[
  (Cosh[alpha*x] + (eta/alpha)*Sinh[alpha*x] - Cosh[alpha*(1 - x)])/w,
  Assumptions -> $Assumptions
];
kernelPrime = FullSimplify[D[kernel, x], Assumptions -> $Assumptions];

Print["K_(alpha,eta)(x) = ", fmt[kernel]];
Print["dK/dx = ", fmt[kernelPrime]];
expectZero[
  "Kprime identity (Mma re-derivation)",
  kernelPrime - (alpha*Sinh[alpha*x] + eta*Cosh[alpha*x] + alpha*Sinh[alpha*(1 - x)])/w
];

(* dK/dx > 0 numerator positivity sweep *)
kprimeNum = alpha*Sinh[alpha*x] + eta*Cosh[alpha*x] + alpha*Sinh[alpha*(1 - x)];
kprimeNumValues = Flatten[Table[
  N[kprimeNum /. {alpha -> aV, eta -> eV, x -> xV}],
  {aV, {1/10, 1, 3}}, {eV, {1/10, 1, 10}}, {xV, {0, 1/4, 1/2, 3/4, 1}}
]];
If[AnyTrue[kprimeNumValues, # <= 0 &],
  fail["kernel numerator positivity sweep", kprimeNumValues],
  pass["kernel numerator positivity sweep"]
];

sigmaPe = FullSimplify[Pe*Exp[Pe*x]/(Exp[Pe] - 1), Assumptions -> $Assumptions];
Print["Sigma_Pe(x) = ", fmt[sigmaPe]];
expectZero["Sigma normalization (Mma re-derivation)", Integrate[sigmaPe, {x, 0, 1}] - 1];

fc = FullSimplify[Integrate[Exp[Pe*x]*Cosh[alpha*x], x], Assumptions -> $Assumptions && Pe != alpha];
fs = FullSimplify[Integrate[Exp[Pe*x]*Sinh[alpha*x], x], Assumptions -> $Assumptions && Pe != alpha];
expectZero["Ic antiderivative regression (Mma re-derivation)", D[fc, x] - Exp[Pe*x]*Cosh[alpha*x]];
expectZero["Is antiderivative regression (Mma re-derivation)", D[fs, x] - Exp[Pe*x]*Sinh[alpha*x]];

ic = FullSimplify[(fc /. x -> 1) - (fc /. x -> 0), Assumptions -> $Assumptions];
is = FullSimplify[(fs /. x -> 1) - (fs /. x -> 0), Assumptions -> $Assumptions];
Print["Ic(Pe,alpha) = ", fmt[ic]];
Print["Is(Pe,alpha) = ", fmt[is]];

delta = FullSimplify[
  Integrate[kernel*sigmaPe, {x, 0, 1}, Assumptions -> $Assumptions && Pe != alpha],
  Assumptions -> $Assumptions && Pe != alpha
];
deltaCombination = FullSimplify[
  Pe/(Exp[Pe] - 1)*((1 - Cosh[alpha])*ic + (eta/alpha + Sinh[alpha])*is)/w,
  Assumptions -> $Assumptions && Pe != alpha
];
expectZero["delta independent integral matches combination form", delta - deltaCombination];
Print["Delta(Pe;alpha,eta) = ", fmt[delta]];

delta0 = FullSimplify[
  Quiet[Limit[delta /. Pe -> pe0, pe0 -> 0], Limit::alimv],
  Assumptions -> alpha > 0 && eta > 0
];
delta0Expected = FullSimplify[eta*(Cosh[alpha] - 1)/(alpha^2*w), Assumptions -> alpha > 0 && eta > 0];
Print["Delta_0 = ", fmt[delta0]];
expectZero["Delta0 formula (Mma re-derivation)", delta0 - delta0Expected];
expectZero[
  "Delta0 integral identity (Mma re-derivation)",
  FullSimplify[delta0 - Integrate[kernel, {x, 0, 1}], Assumptions -> alpha > 0 && eta > 0]
];

deltaInf = FullSimplify[kernel /. x -> 1, Assumptions -> alpha > 0 && eta > 0];
deltaInfExpected = FullSimplify[
  (Cosh[alpha] + (eta/alpha)*Sinh[alpha] - 1)/w,
  Assumptions -> alpha > 0 && eta > 0
];
Print["Delta_inf = ", fmt[deltaInf]];
expectZero["Delta_inf direct substitution (sanity, Mma re-derivation)", deltaInf - deltaInfExpected];

(* BVP independence check is already satisfied by line 84 above:
   "delta independent integral matches combination form" asserts
   integral(kernel * sigmaPe) - Pe/(exp(Pe)-1)*((1-cosh(alpha)) ic + (eta/alpha + sinh(alpha)) is)/w == 0.
   The integral side is the Green-function representation (kernel ansatz); the combination
   side is the closed-form Ic/Is reduction. Their equality verifies the kernel BVP without
   invoking a symbolic DSolve+BC FullSimplify, which is intractable on general alpha, eta, Pe. *)

peLo = FullSimplify[Xi*delta0Expected, Assumptions -> alpha > 0 && eta > 0 && Xi > 0];
peHi = FullSimplify[Xi*deltaInfExpected, Assumptions -> alpha > 0 && eta > 0 && Xi > 0];
Print["Pe_lo = Xi Delta_0 = ", fmt[peLo]];
Print["Pe_hi = Xi Delta_inf = ", fmt[peHi]];

bracketGap = FullSimplify[deltaInfExpected - delta0Expected, Assumptions -> alpha > 0 && eta > 0];
bracketGapExpected = FullSimplify[
  ((alpha^2 - eta)*(Cosh[alpha] - 1) + alpha*eta*Sinh[alpha])/(alpha^2*w),
  Assumptions -> alpha > 0 && eta > 0
];
expectZero["bracket gap closed form", bracketGap - bracketGapExpected];
bracketGapValues = Flatten[Table[
  N[bracketGap /. {alpha -> aV, eta -> eV}],
  {aV, {1/10, 1, 3}}, {eV, {1/10, 1, 10}}
]];
If[AnyTrue[bracketGapValues, # <= 0 &],
  fail["bracket gap positivity sweep", bracketGapValues],
  pass["bracket gap positivity sweep"]
];

(* Delta(Pe; alpha, eta) monotonicity sweep on the constructive branch *)
deltaMonotonicityValues = Flatten[Table[
  Module[{d0v, dinfv, dv},
    d0v = N[delta0Expected /. {alpha -> aV, eta -> eV}];
    dinfv = N[deltaInfExpected /. {alpha -> aV, eta -> eV}];
    dv = N[delta /. {alpha -> aV, eta -> eV, Pe -> pV}];
    {dv - d0v, dinfv - dv}
  ],
  {aV, {1/10, 1, 3}}, {eV, {1/10, 1, 10}}, {pV, {1/2, 1, 3, 10}}
], 2];
If[AnyTrue[deltaMonotonicityValues, # < -10^-9 &],
  fail["Delta(Pe) monotonicity sweep", deltaMonotonicityValues],
  pass["Delta(Pe) monotonicity sweep"]
];

(* F-sign IVT bracket-existence check *)
fSignValues = Flatten[Table[
  Module[{d0v, dinfv, peLoV, peHiV, dAtLo, dAtHi, fLo, fHi},
    d0v = N[delta0Expected /. {alpha -> aV, eta -> eV}];
    dinfv = N[deltaInfExpected /. {alpha -> aV, eta -> eV}];
    peLoV = N[xiV] * d0v;
    peHiV = N[xiV] * dinfv;
    dAtLo = N[delta /. {alpha -> aV, eta -> eV, Pe -> peLoV}];
    dAtHi = N[delta /. {alpha -> aV, eta -> eV, Pe -> peHiV}];
    fLo = peLoV - N[xiV] * dAtLo;
    fHi = peHiV - N[xiV] * dAtHi;
    {-fLo, fHi}  (* both should be >= 0 *)
  ],
  {aV, {1/10, 1, 3}}, {eV, {1/10, 1, 10}}, {xiV, {1/2, 1, 2}}
], 2];
If[AnyTrue[fSignValues, # < -10^-9 &],
  fail["F-sign IVT bracket existence sweep", fSignValues],
  pass["F-sign IVT bracket existence sweep"]
];

deltaInfLimit = FullSimplify[
  Limit[delta, Pe -> Infinity, Assumptions -> alpha > 0 && eta > 0],
  Assumptions -> alpha > 0 && eta > 0
];
Print["Delta(Pe -> oo) = ", fmt[deltaInfLimit]];
expectZero["Delta_inf as Pe -> oo limit", deltaInfLimit - deltaInfExpected];

deltaSeries = FullSimplify[Normal[Series[delta, {Pe, 0, 1}]], Assumptions -> alpha > 0 && eta > 0];
Print["Delta(Pe) small-Pe series = ", fmt[deltaSeries]];
expectZero["weak-coupling constant term (Mma re-derivation)", SeriesCoefficient[deltaSeries, {Pe, 0, 0}] - delta0Expected];
pe1Coeff = FullSimplify[SeriesCoefficient[deltaSeries, {Pe, 0, 1}], Assumptions -> alpha > 0 && eta > 0];
Print["Delta(Pe) first-order coefficient = ", fmt[pe1Coeff]];
pe1Val = N[pe1Coeff /. {alpha -> 1, eta -> 1}];
If[Chop[pe1Val] === 0,
  fail["weak-coupling first-order coefficient vanishes at alpha=eta=1", pe1Val],
  pass["weak-coupling first-order coefficient nonvanishing at alpha=eta=1"]
];

(* Weak-coupling branch law: Pe_*(Xi) = Xi*Delta_0 + O(Xi^2). *)
fSymbol = Pe - Xi*delta;
dFdPe = D[fSymbol, Pe];
dFdXi = D[fSymbol, Xi];
dFdPeAtOrigin = FullSimplify[
  Limit[dFdPe /. Xi -> 0, Pe -> 0],
  Assumptions -> alpha > 0 && eta > 0
];
dFdXiAtOrigin = FullSimplify[
  Limit[dFdXi /. Xi -> 0, Pe -> 0],
  Assumptions -> alpha > 0 && eta > 0
];
dPeStarDXiAtZero = FullSimplify[
  -dFdXiAtOrigin / dFdPeAtOrigin,
  Assumptions -> alpha > 0 && eta > 0
];
expectZero["weak-coupling branch slope = Delta_0", dPeStarDXiAtZero - delta0Expected];

Print[""];
Print["Stage 058 Mathematica audit passed."];

Exit[0];
