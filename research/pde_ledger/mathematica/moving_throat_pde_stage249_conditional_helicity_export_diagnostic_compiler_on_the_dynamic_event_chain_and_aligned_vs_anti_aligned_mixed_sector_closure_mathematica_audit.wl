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

fmt[expr_] := ToString[InputForm[expr]];
pass[name_String] := Print["PASS: ", name];

fail[name_String, detail_] := (
  Print["FAIL: ", name];
  Print["  residual -> ", fmt[detail]];
  Exit[1];
);

stripConditional[expr_] := expr /. ConditionalExpression[value_, _] :> value;

cleanScalar[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanScalar[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, statement_] := Module[{res},
  res = stripConditional[FullSimplify[statement, Assumptions -> $Assumptions]];
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = N[Abs[N[value, 50] - N[target, 50]], 50];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 249 - CONDITIONAL HELICITY EXPORT DIAGNOSTIC"];

Clear[
  Sfull, Sres, Scov, Gamma0, Gamma1, sig, ah, Rpksym, etah, I0, I1,
  Rintsym, abar
];

$Assumptions = (
  Element[
    {Sfull, Sres, Gamma0, Gamma1, sig, ah, Rpksym, etah, I0, I1,
      Rintsym, abar},
    Reals
  ] &&
  Gamma0 != 0 && etah != 0 && I0 != I1 && ah != 1 &&
  Rpksym != -1 && Rintsym != -1
);

subbanner["M1. Subscale transfer-law source"];

Scov = Sfull - Sres;
fullLedgerRHS = -2 Sfull;
resolvedLedgerRHS = -2 Sres;
subtractedRHS = fullLedgerRHS - resolvedLedgerRHS;

Print["Scov = ", fmt[Scov]];
Print["subtracted RHS = ", fmt[subtractedRHS]];
Print["matched -2 Scov = ", fmt[-2 Scov]];
expectZero["M1 projected-minus-resolved source", FullSimplify[subtractedRHS - (-2 Scov)]];

subbanner["M2-M3. Closure reduction and peak inverse"];

Hdot[s_] := Gamma0 + s Gamma1;
ah = Gamma1/Gamma0;
expectZero["M2 closure reduction", FullSimplify[Hdot[sig] - Gamma0 (1 + sig ah)]];
Clear[ah];

RpkExpr = FullSimplify[
  Hdot[1]/Hdot[-1] /. Gamma1 -> ah Gamma0,
  Assumptions -> $Assumptions
];
peakSolution = Solve[Rpksym == RpkExpr, ah];
ahSolved = stripConditional[ah /. First[peakSolution]];

Print["Rpk expression = ", fmt[RpkExpr]];
Print["ah solved from Rpk = ", fmt[ahSolved]];
expectZero["M3 peak Mobius inverse", FullSimplify[ahSolved - (Rpksym - 1)/(Rpksym + 1)]];

subbanner["M4. Integrated ratio and scale cancellation"];

Hplus = etah (I0 + I1);
Hminus = etah (I0 - I1);
Rint = FullSimplify[Hplus/Hminus, Assumptions -> $Assumptions];

Print["Hplus = ", fmt[Hplus]];
Print["Hminus = ", fmt[Hminus]];
Print["Rint = ", fmt[Rint]];
expectTrue["M4 eta_h cancels", FreeQ[Rint, etah]];
expectZero["M4 integrated ratio", FullSimplify[Rint - (I0 + I1)/(I0 - I1)]];

integratedSolution = Solve[Rintsym == (1 + abar)/(1 - abar), abar];
abarSolved = stripConditional[abar /. First[integratedSolution]];

Print["abar solved from Rint = ", fmt[abarSolved]];
expectZero[
  "M4 integrated Mobius inverse",
  FullSimplify[abarSolved - (Rintsym - 1)/(Rintsym + 1)]
];

subbanner["M5. Session-II benchmark packet"];

peakAligned = 281.79830789;
peakAnti = 56.96878122;
hAligned = 20.58070146;
hAnti = 5.00843357;
RintReport = 4.10920923;
XiTurn = 0.34437471;
lambdaTh = 0.42826825;
vCross = 2.59221845;

RpkBench = peakAligned/peakAnti;
alphaPeakBench = (RpkBench - 1)/(RpkBench + 1);
RfinalBench = hAligned/hAnti;
abarBench = (RintReport - 1)/(RintReport + 1);

Print["Rpk benchmark = ", fmt[N[RpkBench, 16]]];
Print["alpha_pk benchmark = ", fmt[N[alphaPeakBench, 16]]];
Print["Rfinal benchmark = ", fmt[N[RfinalBench, 16]]];
Print["abar benchmark = ", fmt[N[abarBench, 16]]];

expectApprox["M5 R_pk benchmark", RpkBench, 4.94653917, 1*^-7];
expectApprox["M5 alpha_pk benchmark", alphaPeakBench, 0.66366992, 1*^-7];
expectApprox["M5 final-load ratio", RfinalBench, RintReport, 1*^-7];
expectApprox["M5 abar_h benchmark", abarBench, 0.60854999, 1*^-7];
expectTrue[
  "M5 asymmetry ordering",
  0 < abarBench < alphaPeakBench < 1
];
expectTrue["M5 positive event-chain packet", XiTurn > 0 && lambdaTh > 0 && vCross > 0];

banner["Stage 249 Mathematica audit result"];
Print["All Mathematica checks passed."];
