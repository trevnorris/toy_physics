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
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 116 — SCALAR D/N RESPONSE KERNEL"];

Clear[x, piM, kappa, gSrc];
$Assumptions =
  Element[{x, piM, kappa, gSrc}, Reals] &&
  0 <= x <= 1 && piM > 0 && kappa > 0 && gSrc > 0 && kappa != piM;

sigma = piM*Exp[-piM*x]/(1 - Exp[-piM]);
cCoeff = FullSimplify[gSrc*piM/((1 - Exp[-piM])*(kappa^2 - piM^2)), Assumptions -> $Assumptions];
aCoeff = FullSimplify[cCoeff*(kappa*Sinh[kappa] + piM*Exp[-piM])/(kappa*Cosh[kappa]), Assumptions -> $Assumptions];
u = FullSimplify[aCoeff*Sinh[kappa*x] - cCoeff*Cosh[kappa*x] + cCoeff*Exp[-piM*x], Assumptions -> $Assumptions];

residual = FullSimplify[-D[u, {x, 2}] + kappa^2*u - gSrc*sigma, Assumptions -> $Assumptions];
bc0 = FullSimplify[u /. x -> 0, Assumptions -> $Assumptions];
bc1 = FullSimplify[D[u, x] /. x -> 1, Assumptions -> $Assumptions];

sKernel = FullSimplify[(D[u, x] /. x -> 0)/gSrc, Assumptions -> $Assumptions];
sTarget = FullSimplify[
  piM*(kappa*Tanh[kappa] + piM*(Exp[-piM]/Cosh[kappa] - 1))/((1 - Exp[-piM])*(kappa^2 - piM^2)),
  Assumptions -> $Assumptions
];

Print["u(x) = ", fmt[u]];
expectZero["ODE residual", residual];
expectZero["u(0)", bc0];
expectZero["u'(1)", bc1];
expectZero["mouth derivative kernel", sKernel - sTarget];

Clear[kk];
$Assumptions = Element[{kk, piM}, Reals] && kk > 0 && piM > 0;
s0 = FullSimplify[Normal[Series[sTarget /. kappa -> kk, {kk, 0, 0}]], Assumptions -> $Assumptions];
Print["S(Pi,0) = ", fmt[s0]];
expectZero["static-shell limit", s0 - 1];

banner["STAGE 116 — GENERAL TWO-CHANNEL FIXED-POINT LAW"];

Print["Pi = Mplus*S(Pi,kappa_plus) + Mminus*S(Pi,kappa_minus)"];
Print["with"];
Print["  S(Pi,kappa) ="];
Print["  ", fmt[sTarget]];

Print[""];
Print["Stage 116 Mathematica audit passed."];

Exit[0];
