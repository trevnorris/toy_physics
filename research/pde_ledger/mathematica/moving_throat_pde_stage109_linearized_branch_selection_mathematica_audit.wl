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

banner["STAGE 109 — LINEARIZED BRANCH-SELECTION LAW"];

(* Paper card's secondary Checks (Robin/mixed-pole limits, even-coefficient
   preservation under the compensated branch) are exercised at downstream
   stages 110 (Robin), 111 (mixed-pole), and 112 (hybrid compensation). *)

Clear[eps, s, b, a0, a5];
$Assumptions = Element[{eps, s, b, a0, a5}, Reals];

sNorm = 1 + eps*s;
beta = 1 + eps*b;
sigma0 = eps*a0;
sigma5 = eps*a5;

(* Independent derivation: expand numerator and denominator separately, *)
(* then form the ratio via series of 1/denominator. *)
num = Expand[3*(sNorm*beta^5 + 9*sigma5)];
den = Expand[3*sNorm - sigma0];
numLin = Normal[Series[num, {eps, 0, 1}]];
denLin = Normal[Series[den, {eps, 0, 1}]];
invDenLin = Normal[Series[1/denLin, {eps, 0, 1}]];
chiSeries = Expand[Normal[Series[numLin*invDenLin, {eps, 0, 1}]]];
expected = 1 + eps*(5*b + a0/3 + 9*a5);

Print["chi_Q series = ", fmt[chiSeries]];
expectZero["linearized chi law", chiSeries - expected];

(* Scale cancellation: read off the first-order eps coefficient and *)
(* check it has no s-dependence. *)
defectCoeff = FullSimplify[Coefficient[chiSeries, eps, 1], Assumptions -> $Assumptions];
Print["first-order defect coefficient = ", fmt[defectCoeff]];
expectZero["overall scale cancels", D[defectCoeff, s]];

(* Preservation condition: solve chiSeries == 1 directly for a5. *)
a5Pres = FullSimplify[a5 /. First[Solve[chiSeries - 1 == 0, a5, Reals]], Assumptions -> $Assumptions];
Print["a5 preservation condition = ", fmt[a5Pres]];
expectZero["a5 preservation condition + 5 b/9 + a0/27", a5Pres + 5*b/9 + a0/27];
expectZero["preservation substitution", (chiSeries - 1) /. a5 -> a5Pres];

Print[""];
Print["Stage 109 Mathematica audit passed."];

Exit[0];
