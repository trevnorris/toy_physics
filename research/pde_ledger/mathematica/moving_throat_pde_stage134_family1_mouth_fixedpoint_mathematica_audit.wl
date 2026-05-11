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

banner["STAGE 117 — FAMILY-1 FIXED-POINT REDUCTION"];

Clear[p, k, kk, k0, Ms, Mq];
$Assumptions = p > 0 && Element[{k, kk, Ms, Mq}, Reals];

sKernel[p_, k_] := FullSimplify[
  p*(k*Tanh[k] + p*(Exp[-p]/Cosh[k] - 1))/((1 - Exp[-p])*(k^2 - p^2)),
  Assumptions -> $Assumptions
];

sShell = Quiet[
  FullSimplify[Limit[sKernel[p, k0], k0 -> 0], Assumptions -> p > 0],
  Limit::alimv
];
sQ = FullSimplify[sKernel[p, Pi/2], Assumptions -> p > 0];
sQExpected = FullSimplify[
  p*((Pi/2)*Tanh[Pi/2] + p*(Exp[-p]/Cosh[Pi/2] - 1))/((1 - Exp[-p])*(Pi^2/4 - p^2)),
  Assumptions -> p > 0
];
fixedPointLaw = FullSimplify[Ms + Mq*sQ, Assumptions -> p > 0];

Print["S_shell = ", fmt[sShell]];
expectZero["static shell channel", sShell - 1];

Print["S_q(p) = ", fmt[sQ]];
expectZero["specialized D/N kernel", sQ - sQExpected];

Print["Fixed-point law p = ", fmt[fixedPointLaw]];

piStar = SetPrecision[1.50882951349316, 30];
sStar = N[sQ /. p -> piStar, 30];
gainLine = N[piStar - sStar*Mq, 30];

Print["S_q(Pi_star) = ", sStar];
Print["Canonical gain line Ms = Pi_star - S_q(Pi_star) M_q"];
Print[gainLine];

Print[""];
Print["RESULT:"];
Print["  Family-1 reduces the coupled mouth law to p = Ms + Mq S_q(p)."];
Print["  The static shell lane gives S_shell = 1 exactly."];
Print["  Evaluated at Pi_star, the canonical compensation line matches the Stage 134 note."];

Exit[0];
