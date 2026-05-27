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

banner["STAGE 112 — EXACT ROBIN-MIXED COMPENSATION LAW"];

Clear[z, rho, sigma, kappa, gamma];
$Assumptions =
  Element[{z, rho, sigma, kappa, gamma}, Reals] &&
  sigma != 1 && rho - sigma != 3;

lambdaOut = -3 + z^2/3 + z^4/9 + I*z^5/9;
lambdaHyb = Expand[Normal[Series[lambdaOut + rho - sigma/(1 - kappa*z^2 - I*gamma*z^5), {z, 0, 5}]]];

l0 = FullSimplify[Coefficient[lambdaHyb, z, 0], Assumptions -> $Assumptions];
l2 = FullSimplify[Coefficient[lambdaHyb, z, 2], Assumptions -> $Assumptions];
l4 = FullSimplify[Coefficient[lambdaHyb, z, 4], Assumptions -> $Assumptions];
l5 = FullSimplify[Coefficient[lambdaHyb, z, 5]/I, Assumptions -> $Assumptions];

sols = SortBy[Solve[{-l2/l0 == 1/9, l2^2/l0^2 - l4/l0 == 4/81}, {rho, kappa}, Reals], FullSimplify[kappa /. #, Assumptions -> $Assumptions] &];
If[Length[sols] =!= 2, fail["expected two canonical-even solutions", sols]];

solA = sols[[1]];
solB = sols[[2]];

Print["Lambda_hyb(z) = ", fmt[lambdaHyb]];
Print["canonical-even solutions = ", fmt[sols]];

expectZero["branch A rho - sigma", (rho /. solA) - sigma];
expectZero["branch A kappa", kappa /. solA];
expectZero["branch B rho - 4 sigma", (rho /. solB) - 4*sigma];
expectZero["branch B kappa - 1/3", (kappa /. solB) - 1/3];

(* Independent Stage-92 (= stage 109) linearized cross-check on solB.        *)
(* The notes' Branch-selection data box gives (b, a_0, a_5) = (0, 3 sigma_W, *)
(* -sigma_W gamma_W) on the nontrivial compensated branch.  The preservation *)
(* condition a_0/3 + 9 a_5 = sigma_W (1 - 9 gamma_W) = 0 determines gamma_W   *)
(* = 1/9 by an algebraically independent route from the chi_Q-based solve.    *)
a0Def = FullSimplify[(Coefficient[lambdaHyb /. solB, z, 0]) - (-3), Assumptions -> $Assumptions];
a5Def = FullSimplify[Coefficient[lambdaHyb /. solB, z, 5]/I - 1/9, Assumptions -> $Assumptions];
expectZero["independent: a_0 - 3 sigma on solB", a0Def - 3*sigma];
expectZero["independent: a_5 + sigma gamma on solB", a5Def + sigma*gamma];
gammaFromLinear = FullSimplify[gamma /. First[Solve[a0Def/3 + 9*a5Def == 0, gamma, Reals]], Assumptions -> $Assumptions];
Print["gamma_W from linearized preservation = ", fmt[gammaFromLinear]];
expectZero["independent: gamma_W from a_0/3 + 9 a_5 = 0", gammaFromLinear - 1/9];

chiA = FullSimplify[((-l5/l0)/(1/27)) /. solA, Assumptions -> $Assumptions];
chiB = FullSimplify[((-l5/l0)/(1/27)) /. solB, Assumptions -> $Assumptions];

Print["chi_Q branch A = ", fmt[chiA]];
Print["chi_Q branch B = ", fmt[chiB]];

expectZero["chi_Q branch A - (1 - 9 sigma gamma)", chiA - (1 - 9*sigma*gamma)];
expectZero["chi_Q branch B - (1 - 9 sigma gamma)/(1 - sigma)", chiB - (1 - 9*sigma*gamma)/(1 - sigma)];
expectZero["chi_Q branch B at gamma=1/9", (chiB /. gamma -> 1/9) - 1];

scaledIdentity = Expand[((lambdaHyb /. solB) - (1 - sigma)*lambdaOut) /. gamma -> 1/9];
Print["scaled identity on branch B = ", fmt[scaledIdentity]];
expectZero["scaled identity on branch B", scaledIdentity];

Print[""];
Print["Stage 112 Mathematica audit passed."];

Exit[0];
