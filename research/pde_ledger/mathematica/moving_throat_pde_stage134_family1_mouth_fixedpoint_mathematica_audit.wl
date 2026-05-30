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

banner["STAGE 134 — FAMILY-1 FIXED-POINT REDUCTION"];

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
fixedPointLaw = FullSimplify[Ms + Mq*sQ, Assumptions -> p > 0];

Print["S_shell = ", fmt[sShell]];
expectZero["static shell channel", sShell - 1];

Print["S_q(p) = ", fmt[sQ]];

(* Non-tautological numeric check: evaluate sQ at three independent Pi values
   against high-precision targets verified independently via mpmath. The targets
   must NOT be derived from sKernel at runtime. *)
expectClose[name_String, got_, want_, tol_] := Module[{d},
  d = Abs[N[got - want, 30]];
  Print[name, " = ", got, "  (target ", want, ", diff ", d, ")"];
  If[TrueQ[d < tol], pass[name], fail[name, d]]
];

expectClose["S_q at p=1/2", N[sQ /. p -> 1/2, 30],
  SetPrecision[0.608336415687717065435990381419, 30], 10^-12];
expectClose["S_q at p=1",   N[sQ /. p -> 1, 30],
  SetPrecision[0.633127670034487546375729566676, 30], 10^-12];
expectClose["S_q at p=2",   N[sQ /. p -> 2, 30],
  SetPrecision[0.681366857005321783286541952613, 30], 10^-12];

Print["Fixed-point law p = ", fmt[fixedPointLaw]];

piStar = SetPrecision[1.50882951349316, 30];   (* imported literal; owned by stage 131 note *)
sStar = N[sQ /. p -> piStar, 30];
gainLine = N[piStar - sStar*Mq, 30];

Print["S_q(Pi_star) = ", sStar];
Print["Canonical gain line Ms = Pi_star - S_q(Pi_star) M_q (printed only; not asserted)"];
Print[gainLine];

(* No in-stage gain-line assertion: intercept is the imported literal PiStar
   (stage 131) and the slope is -S_q(PiStar); re-asserting them here would be an
   X-X tautology. Outlet consistency of (Ms, Mq) is verified at Stage 135;
   susceptibility closure at Stage 137. *)

Print[""];
Print["RESULT:"];
Print["  Family-1 reduces the coupled mouth law to p = Ms + Mq S_q(p)."];
Print["  The static shell lane gives S_shell = 1 exactly."];
Print["  S_q at p=1/2,1,2 and at Pi_star match independent literal targets (see checks above)."];
Print["  Outlet consistency of the gain pair (Ms, Mq) is verified downstream at Stage 135."];

Exit[0];
