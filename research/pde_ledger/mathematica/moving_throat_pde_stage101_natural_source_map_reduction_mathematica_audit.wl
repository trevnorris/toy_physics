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

banner["STAGE 084 — NATURAL SOURCE-MAP REDUCTION"];

Clear[mHat0, chiQ, nQ, deltaQ];
$Assumptions =
  Element[{mHat0, chiQ, nQ, deltaQ}, Reals] &&
  mHat0 > 0 && chiQ > 0;

nQSol = First[Solve[mHat0^2*chiQ*nQ == 1, nQ]];
nQExact = FullSimplify[nQ /. nQSol, Assumptions -> $Assumptions];

Print["NQ from exact factorized odd normalization = ", fmt[nQExact]];
Print["point-particle natural branch mhat0->1 gives = ", fmt[FullSimplify[nQExact /. mHat0 -> 1, Assumptions -> $Assumptions]]];
Print["canonical compact outgoing branch chiQ=1 gives = ", fmt[FullSimplify[nQExact /. {mHat0 -> 1, chiQ -> 1}, Assumptions -> $Assumptions]]];

expectZero["point-particle natural branch reduction", (nQExact /. mHat0 -> 1) - 1/chiQ];
expectZero["canonical compact outgoing branch gives NQ=1", (nQExact /. {mHat0 -> 1, chiQ -> 1}) - 1];

exprDelta = FullSimplify[1/(1 + deltaQ) - 1, Assumptions -> $Assumptions];
seriesDelta = Expand[Normal[Series[exprDelta, {deltaQ, 0, 2}]]];
Print["NQ - 1 in terms of DeltaQ = ", fmt[exprDelta]];
Print["small-DeltaQ series = ", fmt[seriesDelta]];
expectZero["exact replacement chiQ=1+DeltaQ", FullSimplify[(nQExact /. {mHat0 -> 1, chiQ -> 1 + deltaQ}) - 1/(1 + deltaQ), Assumptions -> $Assumptions]];

Print[""];
Print["Stage 101 Mathematica audit passed."];

Exit[0];
