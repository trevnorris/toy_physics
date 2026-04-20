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

banner["PART I — EXACT SELECTED OVERLAP DERIVATIVE"];

Clear[a, dK, alpha, x0, x1, beta0];
$Assumptions = Element[{a, dK, alpha, x0, x1, beta0}, Reals] &&
  a > 0 && dK > 0 && alpha >= 0 && x0 > 0 && x1 > 0 && beta0 > 0;

sigma = x0 + x1;
deltaKappa = x0 - x1;
kappaProd = x0*x1;
r = Sqrt[(dK + alpha*deltaKappa)^2 + 4*alpha^2*kappaProd];

lamMinus = FullSimplify[(2*a + dK - alpha*sigma - r)/2, Assumptions -> $Assumptions];
lamPlus = FullSimplify[(2*a + dK - alpha*sigma + r)/2, Assumptions -> $Assumptions];
sMinus = FullSimplify[-D[lamMinus, alpha], Assumptions -> $Assumptions];
dsExact = FullSimplify[D[sMinus, alpha], Assumptions -> $Assumptions];
dsExpected = FullSimplify[2*dK^2*kappaProd/r^3, Assumptions -> $Assumptions];
expectZero["ds_-/dalpha exact formula", dsExact - dsExpected];
Print["ds_-/dalpha = ", fmt[dsExpected]];

banner["PART II — EXACT PREFATOR MONOTONICITY IDENTITY"];
p0Sel = FullSimplify[beta0*sMinus/lamMinus, Assumptions -> $Assumptions];
dPExact = FullSimplify[D[p0Sel, alpha], Assumptions -> $Assumptions];
dPExpected = FullSimplify[beta0*(dsExpected*lamMinus + sMinus^2)/lamMinus^2, Assumptions -> $Assumptions];
expectZero["exact monotonicity identity", dPExact - dPExpected];
Print["dP0_-/dalpha = ", fmt[dPExpected]];

banner["PART III — INITIAL VALUES"];
expectZero["lambda_-(0)", (lamMinus /. alpha -> 0) - a];
expectZero["s_-(0)", (sMinus /. alpha -> 0) - x0];
expectZero["P0_-(0)", (p0Sel /. alpha -> 0) - beta0*x0/a];

banner["PART IV — EXACT SOFTENING THRESHOLD"];
t0 = FullSimplify[(a + dK)*x0 + a*x1, Assumptions -> $Assumptions];
alphaCrit = FullSimplify[a*(a + dK)/t0, Assumptions -> $Assumptions];
expectZero["det factorization", Expand[lamMinus*lamPlus - (a*(a + dK) - alpha*t0)]];
expectZero["alpha_crit condition", (a*(a + dK) - alpha*t0) /. alpha -> alphaCrit];

lambdaPlusCrit = FullSimplify[lamPlus /. alpha -> alphaCrit, Assumptions -> $Assumptions];
radCrit = a^4*x0^2 + 2*a^4*x0*x1 + a^4*x1^2 + 4*a^3*dK*x0^2 + 4*a^3*dK*x0*x1 +
  6*a^2*dK^2*x0^2 + 2*a^2*dK^2*x0*x1 + 4*a*dK^3*x0^2 + dK^4*x0^2;
expectZero["threshold radical square identity", radCrit - (a^2*x1 + (a + dK)^2*x0)^2];
Print["alpha_crit = ", fmt[alphaCrit]];
Print["lambda_+^(crit) = ", fmt[lambdaPlusCrit]];

banner["PART V — STABLE-SIDE DIVERGENCE STRUCTURE"];
expectZero["lambda_- * lambda_+ - T0*(alpha_crit-alpha)", Expand[lamMinus*lamPlus - t0*(alphaCrit - alpha)]];
p0Factored = FullSimplify[beta0*sMinus*lamPlus/(t0*(alphaCrit - alpha)), Assumptions -> $Assumptions];
expectZero["P0_- factorization", p0Sel - p0Factored];
Print["P0_- = ", fmt[p0Factored]];

banner["STAGE 14 AUDIT COMPLETE"];
Print["Verified:"];
Print["  exact overlap derivative ds_-/dalpha = 2 DK^2 kappa0^2 kappa1^2 / R^3"];
Print["  exact monotonicity identity for dP0_-/dalpha"];
Print["  exact initial-value formulas at alpha = 0"];
Print["  exact softening threshold alpha_crit and determinant factorization"];
Print["  exact stable-side factorization P0_- = beta0 s_- lambda_+ / [T0 (alpha_crit-alpha)]"];

Exit[0];
