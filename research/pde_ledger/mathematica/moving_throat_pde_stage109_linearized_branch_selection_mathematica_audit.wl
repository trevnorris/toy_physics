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

banner["STAGE 092 — LINEARIZED BRANCH-SELECTION LAW"];

Clear[eps, s, b, a0, a5];
$Assumptions = Element[{eps, s, b, a0, a5}, Reals];

sNorm = 1 + eps*s;
beta = 1 + eps*b;
sigma0 = eps*a0;
sigma5 = eps*a5;

chiQ = FullSimplify[3*(sNorm*beta^5 + 9*sigma5)/(3*sNorm - sigma0), Assumptions -> $Assumptions];
chiSeries = Expand[Normal[Series[chiQ, {eps, 0, 1}]]];
expected = 1 + eps*(5*b + a0/3 + 9*a5);

Print["chi_Q series = ", fmt[chiSeries]];
expectZero["linearized chi law", chiSeries - expected];

coeff = FullSimplify[Expand[(chiSeries - 1)/eps], Assumptions -> $Assumptions];
Print["first-order defect coefficient = ", fmt[coeff]];
expectZero["overall scale cancels", D[coeff, s]];

a5Pres = FullSimplify[a5 /. First[Solve[coeff == 0, a5, Reals]], Assumptions -> $Assumptions];
Print["a5 preservation condition = ", fmt[a5Pres]];
expectZero["a5 preservation condition + 5 b/9 + a0/27", a5Pres + 5*b/9 + a0/27];
expectZero["preservation substitution", coeff /. a5 -> a5Pres];

Print[""];
Print["Stage 109 Mathematica audit passed."];

Exit[0];
