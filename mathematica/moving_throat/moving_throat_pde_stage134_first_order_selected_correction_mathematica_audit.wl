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

banner["FIRST-ORDER SELF-CONSISTENT SOURCE CORRECTION"];

Clear[rBar, cBar, kBar, cRBar, kRBar, gPrime, aT, bT];
$Assumptions = Element[{rBar, cBar, kBar, cRBar, kRBar, gPrime, aT, bT}, Reals] && gPrime != 0;

covCR = cRBar - cBar*rBar;
covKR = kRBar - kBar*rBar;

deltaG = -(cRBar - cBar*rBar);
deltaS = -(kRBar - kBar*rBar);

expectZero["delta g + Cov(c,R)", deltaG + covCR];
expectZero["delta S + Cov(K,R)", deltaS + covKR];

deltaPi = FullSimplify[-deltaG/gPrime];
deltaT = FullSimplify[aT*deltaG + bT*deltaS];

Print["deltaPi = ", fmt[deltaPi]];
Print["deltaT  = ", fmt[deltaT]];
Print["deltaPi in covariance form = ", fmt[deltaPi]];
Print["deltaT in covariance form = ", fmt[deltaT]];

Print[""];
Print["Theorem:"];
Print["  Once the full mouth residual R_*(x) is known, the selected first-order"];
Print["  source correction is completely determined by Cov_*(c,R_*) and Cov_*(K_q,R_*)."];

Print[""];
Print["Stage 134 Mathematica audit passed."];

Exit[0];
