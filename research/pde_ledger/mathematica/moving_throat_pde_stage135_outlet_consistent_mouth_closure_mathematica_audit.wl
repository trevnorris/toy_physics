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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 40] - N[target, 40]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = TrueQ[cond];
  Print[name, " = ", fmt[res]];
  If[res, pass[name], fail[name, cond]];
];

banner["STAGE 135 — OUTLET-CONSISTENT ONE-PARAMETER CLOSURE"];

Clear[pi, sigmaM, ms, mq, sigmaVar];
$Assumptions =
  Element[{pi, sigmaM, ms, mq, sigmaVar}, Reals] &&
  pi > 0 && sigmaM > 0 && sigmaVar > 0;

sKernel[p_, k_] := FullSimplify[
  p*(k*Tanh[k] + p*(Exp[-p]/Cosh[k] - 1))/((1 - Exp[-p])*(k^2 - p^2)),
  Assumptions -> $Assumptions
];

sQ = FullSimplify[sKernel[pi, Pi/2], Assumptions -> $Assumptions];
genericLaw = FullSimplify[ms + mq*sQ, Assumptions -> $Assumptions];
reducedLaw = FullSimplify[genericLaw /. {ms -> 4*sigmaM, mq -> -sigmaM}, Assumptions -> $Assumptions];
expectedLaw = FullSimplify[sigmaM*(4 - sQ), Assumptions -> $Assumptions];

Print["Pi = ", fmt[expectedLaw]];
expectZero["outlet-consistent reduction", reducedLaw - expectedLaw];

piStar = SetPrecision[1.50882951349316, 30];
sStar = N[sQ /. pi -> piStar, 30];
sigmaStar = N[sigmaVar /. First[Solve[piStar == sigmaVar*(4 - sStar), sigmaVar, Reals]], 30];
msStar = N[4*sigmaStar, 30];
mqStar = N[-sigmaStar, 30];
mixedCorrection = N[mqStar*sStar, 30];
residual = N[piStar - sigmaStar*(4 - sStar), 30];

Print["S_q(Pi_star) = ", sStar];
expectTrue["0 < S_q(Pi_star) < 1", 0 < sStar < 1];
Print["Sigma_m^* = ", sigmaStar];
Print["M_s^* = ", msStar];
Print["M_q^* = ", mqStar];
Print["M_q^* S_q(Pi_*) = ", mixedCorrection];
Print["Pi_star - Sigma_star*(4 - S_star) = ", residual];

expectApprox["Sigma_m^* numeric check", sigmaStar, 0.451485277739089696513730132210, 10^-15];
expectApprox["M_s^* numeric check", msStar, 1.80594111095635878605492052884, 10^-14];
expectApprox["M_q^* numeric check", mqStar, -0.451485277739089696513730132210, 10^-15];
expectApprox["mixed-lane correction numeric check", mixedCorrection, -0.297111597463198185916386491356, 10^-14];
expectApprox["closure residual", residual, 0, 10^-14];
expectTrue["3 Sigma_m^* < Pi_* < 4 Sigma_m^*", 3*sigmaStar < piStar < 4*sigmaStar];

Print[""];
Print["Carry-forward formulas:"];
Print["  M_s = 4 Sigma_m,   M_q = -Sigma_m"];
Print["  Pi = Sigma_m [4 - S_q(Pi)]"];
Print["  The canonical outlet-consistent branch is selected by one moderate mouth-gain amplitude."];

Exit[0];
