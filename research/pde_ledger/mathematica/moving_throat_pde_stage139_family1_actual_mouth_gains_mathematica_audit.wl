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
rF = N[Sqrt[(12/Pi^2)*(37/20)^2 - 1], 30];
(* piStar and sQStar (= S_q evaluated at piStar) imported from Stage 134. *)
piStar = SetPrecision[1.50882951349316, 30];
sQStar = SetPrecision[0.658075937605429, 30];

rQNat = N[(1 - rF)^2/(1 + rF^2), 30];
mSNat = N[piStar/(1 - rQNat*sQStar), 30];
mQNat = N[-rQNat*mSNat, 30];

(* Derive g_c on the compensated branch by solving the defining condition
   (g_c - r_F1)^2 / (1 + r_F1^2) == 1/4, branch g_c < r_F1 (lower
   compensated branch, see notes/stages/moving_throat_pde_stage139*.md
   section "Exact compensated branch"). *)
gMinusSolutions = gc /. Solve[(gc - rF)^2 == (1 + rF^2)/4 && gc < rF, gc, Reals];
If[Length[gMinusSolutions] =!= 1,
  Print["FAIL: g_minus branch selection ambiguous, solutions = ", gMinusSolutions];
  Exit[1];
];
gMinus = N[First[gMinusSolutions], 30];
(* Cross-check the derived value against the closed form quoted in the notes. *)
expectApprox["g_minus closed form", gMinus, rF - Sqrt[1 + rF^2]/2, 10^-25];
rQComp = N[(gMinus - rF)^2/(1 + rF^2), 30];
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
expectApprox["R_q^nat", rQNat, 0.145454452260421, tolLit];
expectApprox["M_s^nat,*", mSNat, 1.66854252965624, tolLit];
expectApprox["M_q^nat,*", mQNat, -0.242696939724365, tolLit];
expectApprox["M_s^comp,*", mSComp, 1.80594111095636, tolLit];
expectApprox["M_q^comp,*", mQComp, -0.451485277739090, tolLit];
expectApprox["outlet consistency nat", piStar, mSNat + mQNat*sQStar, tolAlg];
expectApprox["outlet consistency comp", piStar, mSComp + mQComp*sQStar, tolAlg];
expectApprox["R_q^comp - 1/4", rQComp, 1/4, tolAlg];

Print[""];
Print["Stage 139 Mathematica audit passed."];

Exit[0];
