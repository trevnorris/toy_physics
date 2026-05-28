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

banner["STAGE 133 — SCALAR D/N RESPONSE KERNEL"];

Clear[x, piM, kappa, gSrc, uFun];
$Assumptions =
  Element[{x, piM, kappa, gSrc}, Reals] &&
  0 <= x <= 1 && piM > 0 && kappa > 0 && gSrc > 0 && kappa != piM;

sigma = piM*Exp[-piM*x]/(1 - Exp[-piM]);

(* Independent derivation: let DSolveValue solve the D/N problem from scratch. *)
uSol = DSolveValue[
  {-uFun''[x] + kappa^2*uFun[x] == gSrc*sigma, uFun[0] == 0, uFun'[1] == 0},
  uFun[x],
  x
];
u = FullSimplify[uSol, Assumptions -> $Assumptions];

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

banner["STAGE 133 — GENERAL TWO-CHANNEL FIXED-POINT LAW"];

Print["Pi = Mplus*S(Pi,kappa_plus) + Mminus*S(Pi,kappa_minus)"];
Print["with"];
Print["  S(Pi,kappa) ="];
Print["  ", fmt[sTarget]];

Print[""];
Print["Stage 133 Mathematica audit passed."];

Exit[0];
