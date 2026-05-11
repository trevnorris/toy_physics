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

banner["STAGE 085 — HIGHER ODD IRRELEVANCE"];

Clear[omega, omegaQ, chiQ, tauQ];
$Assumptions =
  Element[{omega, omegaQ, chiQ, tauQ}, Reals] &&
  omegaQ > 0 && chiQ > 0 && tauQ > 0;

sigmaCan = FullSimplify[(9/8)/omegaQ^5, Assumptions -> $Assumptions];
yRet = 3/4 + (1/4)/(1 - omega^2/omegaQ^2 - I*chiQ*sigmaCan*omega^5 - I*tauQ*omega^7);
ySer5 = Expand[Normal[Series[yRet, {omega, 0, 5}]]];
ySer7 = Expand[Normal[Series[yRet, {omega, 0, 7}]]];

Print["series through O(omega^5) = ", fmt[ySer5]];
Print["series through O(omega^7) = ", fmt[ySer7]];

tauCoeff5 = FullSimplify[D[Coefficient[ySer5, omega, 5]/I, tauQ], Assumptions -> $Assumptions];
tauCoeff7 = FullSimplify[D[Coefficient[ySer7, omega, 7]/I, tauQ], Assumptions -> $Assumptions];
Print["tauQ coefficient in omega^5 term = ", fmt[tauCoeff5]];
Print["tauQ coefficient in omega^7 term = ", fmt[tauCoeff7]];

expectZero["tauQ irrelevance at omega^5", tauCoeff5];
expectZero["tauQ coefficient at omega^7 - 1/4", tauCoeff7 - 1/4];
expectZero["check canonical odd coefficient", (Coefficient[ySer5, omega, 5]/I) - chiQ*(9/32)/omegaQ^5];

Print[""];
Print["Stage 102 Mathematica audit passed."];

Exit[0];
