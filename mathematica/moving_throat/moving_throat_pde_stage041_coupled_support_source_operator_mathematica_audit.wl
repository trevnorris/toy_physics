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
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 041 — COUPLED SUPPORT/SOURCE OPERATOR"];

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
  "Kprime identity",
  kernelPrime - (alpha*Sinh[alpha*x] + eta*Cosh[alpha*x] + alpha*Sinh[alpha*(1 - x)])/w
];

sigmaPe = FullSimplify[Pe*Exp[Pe*x]/(Exp[Pe] - 1), Assumptions -> $Assumptions];
Print["Sigma_Pe(x) = ", fmt[sigmaPe]];
expectZero["Sigma normalization", Integrate[sigmaPe, {x, 0, 1}] - 1];

fc = Exp[Pe*x]*(Pe*Cosh[alpha*x] - alpha*Sinh[alpha*x])/(Pe^2 - alpha^2);
fs = Exp[Pe*x]*(Pe*Sinh[alpha*x] - alpha*Cosh[alpha*x])/(Pe^2 - alpha^2);
expectZero["Ic antiderivative check", D[fc, x] - Exp[Pe*x]*Cosh[alpha*x]];
expectZero["Is antiderivative check", D[fs, x] - Exp[Pe*x]*Sinh[alpha*x]];

ic = FullSimplify[(fc /. x -> 1) - (fc /. x -> 0), Assumptions -> $Assumptions];
is = FullSimplify[(fs /. x -> 1) - (fs /. x -> 0), Assumptions -> $Assumptions];
Print["Ic(Pe,alpha) = ", fmt[ic]];
Print["Is(Pe,alpha) = ", fmt[is]];

delta = FullSimplify[
  Pe/(Exp[Pe] - 1)*((1 - Cosh[alpha])*ic + (eta/alpha + Sinh[alpha])*is)/w,
  Assumptions -> $Assumptions
];
Print["Delta(Pe;alpha,eta) = ", fmt[delta]];

delta0 = FullSimplify[
  Quiet[Limit[delta /. Pe -> pe0, pe0 -> 0], Limit::alimv],
  Assumptions -> alpha > 0 && eta > 0
];
delta0Expected = FullSimplify[eta*(Cosh[alpha] - 1)/(alpha^2*w), Assumptions -> alpha > 0 && eta > 0];
Print["Delta_0 = ", fmt[delta0]];
expectZero["Delta0 formula", delta0 - delta0Expected];
expectZero[
  "Delta0 integral identity",
  FullSimplify[delta0 - Integrate[kernel, {x, 0, 1}], Assumptions -> alpha > 0 && eta > 0]
];

deltaInf = FullSimplify[kernel /. x -> 1, Assumptions -> alpha > 0 && eta > 0];
deltaInfExpected = FullSimplify[
  (Cosh[alpha] + (eta/alpha)*Sinh[alpha] - 1)/w,
  Assumptions -> alpha > 0 && eta > 0
];
Print["Delta_inf = ", fmt[deltaInf]];
expectZero["Delta_inf formula", deltaInf - deltaInfExpected];

peLo = FullSimplify[Xi*delta0Expected, Assumptions -> alpha > 0 && eta > 0 && Xi > 0];
peHi = FullSimplify[Xi*deltaInfExpected, Assumptions -> alpha > 0 && eta > 0 && Xi > 0];
Print["Pe_lo = Xi Delta_0 = ", fmt[peLo]];
Print["Pe_hi = Xi Delta_inf = ", fmt[peHi]];

deltaSeries = FullSimplify[Normal[Series[delta, {Pe, 0, 1}]], Assumptions -> alpha > 0 && eta > 0];
Print["Delta(Pe) small-Pe series = ", fmt[deltaSeries]];
expectZero["weak-coupling constant term", SeriesCoefficient[deltaSeries, {Pe, 0, 0}] - delta0Expected];

Print[""];
Print["Stage 041 Mathematica audit passed."];

Exit[0];
