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

expectApprox[name_String, value_, target_, tol_] := Module[{delta},
  delta = N[value - target, 40];
  Print[name, " residual = ", fmt[delta]];
  If[TrueQ[Abs[delta] < tol], pass[name], fail[name, delta]];
];

banner["STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS"];

(* r_F1 closed form from Stage 121. *)
rFClosed = Sqrt[(12/Pi^2)*(37/20)^2 - 1];
rF = N[rFClosed, 50];
(* Self-matched susceptibility closure (Sigma_0 = mS) is established at Stage 140, not here;
   Stage 139 only evaluates the gain pair on the Family-1 branch. See P4-42. *)

Clear[p, x];
$Assumptions = Element[{p, x}, Reals] && p > 0 && 0 <= x <= 1;

(* g_- is the LOWER compensated branch g_c = rF - (1/2) Sqrt[1 + rF^2] (notes stage139 lines
   91-100). R_q^comp = 1/4 is DEFINITIONAL on this branch (true for any rF) AND branch-blind.
   Compute g_- directly as the closed form; the value check below discriminates the lower branch. *)
gMinusClosed = FullSimplify[rFClosed - Sqrt[1 + rFClosed^2]/2, Assumptions -> True];
gMinus = N[gMinusClosed, 50];
gMinusRoot = N[gMinusClosed, 80];

(* Independent Pi_* route: clear the positive denominator of g(Pi) == g_- and solve the
   numerator with a bracketed high-precision root. The residual uses (Exp[p] - 1). *)
gThresholdResidual[z_] := (
  2*z*(2*z*Exp[z] + Pi) - gMinusRoot*(4*z^2 + Pi^2)*(Exp[z] - 1)
);
piStar = p /. FindRoot[gThresholdResidual[p] == 0, {p, 1.4, 1.6},
  WorkingPrecision -> 80, AccuracyGoal -> 35, PrecisionGoal -> 35, MaxIterations -> 100];

(* Independent S_q at Pi_star route: direct source-moment quadrature of Sigma*K_q. *)
kappaQ = Pi/2;
sigmaStarX = piStar*Exp[-piStar*x]/(1 - Exp[-piStar]);
kQStarX = Cosh[kappaQ*(1 - x)]/Cosh[kappaQ];
sQStar = NIntegrate[
  Evaluate[SetPrecision[sigmaStarX*kQStarX, 70]],
  {x, 0, 1},
  WorkingPrecision -> 70,
  AccuracyGoal -> 30,
  PrecisionGoal -> 30,
  MaxRecursion -> 30
];

sigmaPX = p*Exp[-p*x]/(1 - Exp[-p]);
kQPX = Cosh[kappaQ*(1 - x)]/Cosh[kappaQ];
sQDirectMoment = Integrate[sigmaPX*kQPX, {x, 0, 1},
  Assumptions -> p > 0, GenerateConditions -> False];
sQDirectAtStar = N[sQDirectMoment /. p -> piStar, 50];

rQNat = N[(1 - rFClosed)^2/(1 + rFClosed^2), 50];
mSNat = N[piStar/(1 - rQNat*sQStar), 30];
mQNat = N[-rQNat*mSNat, 30];

rQComp = N[FullSimplify[(gMinusClosed - rFClosed)^2/(1 + rFClosed^2), Assumptions -> True], 50];
mSComp = N[piStar/(1 - rQComp*sQStar), 30];
mQComp = N[-rQComp*mSComp, 30];

Print["r_F1 = ", fmt[rF]];
Print["S_q(Pi_*) = ", fmt[sQStar]];
Print["R_q^nat = ", fmt[rQNat]];
Print["M_s^nat,* = ", fmt[mSNat]];
Print["M_q^nat,* = ", fmt[mQNat]];
Print["g_-^F1 = ", fmt[gMinus]];
Print["R_q^comp = ", fmt[rQComp]];
Print["M_s^comp,* = ", fmt[mSComp]];
Print["M_q^comp,* = ", fmt[mQComp]];
Print["shell gain fractional shift = ", fmt[N[mSComp/mSNat - 1, 20]]];
Print["mixed gain magnitude ratio = ", fmt[N[Abs[mQComp]/Abs[mQNat], 20]]];

(* Numerical-deliverable assertions against the boxed values in
   notes/stages/moving_throat_pde_stage139_family1_actual_mouth_gains.md. *)
tolLit = 10^-12;
tolAlg = 10^-25;

expectApprox["r_F1", rF, 1.77799353547498, tolLit];
expectApprox["Pi_* cleared-denominator residual", gThresholdResidual[piStar], 0, 10^-25];
expectApprox["Pi_* value", piStar, 1.50882951349316, tolLit];
expectApprox["S_q(Pi_*) value", sQStar, 0.658075937605429, tolLit];
expectApprox["S_q quadrature vs direct Integrate", sQStar, sQDirectAtStar, 10^-25];
expectApprox["R_q^nat", rQNat, 0.145454452260421, tolLit];
(* (R2 anchor) g_-^F1 VALUE vs the canonical cross-stage literal (owned at 127/142/144/164/169);
   DISCRIMINATES the lower branch (g_- ~ 0.758) from the upper (g_+ ~ 2.79), which R_q = 1/4
   cannot. Falsifiable: a sign/branch/rF typo gives ~2.79 and FAILS. *)
expectApprox["g_-^F1 value", gMinus, 0.758035078944662826919680890414, tolLit];
expectApprox["M_s^nat,*", mSNat, 1.66854252965624, tolLit];
expectApprox["M_q^nat,*", mQNat, -0.242696939724365, tolLit];
expectApprox["M_s^comp,*", mSComp, 1.80594111095636, tolLit];
expectApprox["M_q^comp,*", mQComp, -0.451485277739090, tolLit];
(* Structural sanity only (true by construction, NOT a literal check): *)
Print["outlet form residual nat structural = ", fmt[N[piStar - (mSNat + mQNat*sQStar), 5]]];
Print["outlet form residual comp structural = ", fmt[N[piStar - (mSComp + mQComp*sQStar), 5]]];
expectApprox["R_q^comp - 1/4 (definitional-consistency)", rQComp, 1/4, tolAlg];

Print[""];
Print["Stage 139 Mathematica audit passed."];

Exit[0];
