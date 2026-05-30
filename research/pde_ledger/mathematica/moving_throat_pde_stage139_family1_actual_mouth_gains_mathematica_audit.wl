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
(* Self-matched susceptibility closure (Sigma_0 = mS) is established at Stage 140, not here;
   Stage 139 only evaluates the gain pair on the Family-1 branch. See P4-42. *)

rQNat = N[(1 - rF)^2/(1 + rF^2), 30];
mSNat = N[piStar/(1 - rQNat*sQStar), 30];
mQNat = N[-rQNat*mSNat, 30];

(* g_- is the LOWER compensated branch g_c = rF - (1/2) Sqrt[1 + rF^2] (notes stage139 lines
   91-100). R_q^comp = 1/4 is DEFINITIONAL on this branch (true for any rF) AND branch-blind.
   Compute g_- DIRECTLY as the closed form (a sanctioned mirror of the SymPy route; rF is
   independently anchored at line 71 and the branch is discriminated by the g_-^F1 value check
   below) — NOT by solving (gc - rF)^2 == (1 + rF^2)/4, which would re-bake 1/4. *)
gMinus = N[rF - Sqrt[1 + rF^2]/2, 30];
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
(* (R2 anchor) g_-^F1 VALUE vs the canonical cross-stage literal (owned at 127/142/144/164/169);
   DISCRIMINATES the lower branch (g_- ~ 0.758) from the upper (g_+ ~ 2.79), which R_q = 1/4
   cannot. Falsifiable: a sign/branch/rF typo gives ~2.79 and FAILS. *)
expectApprox["g_-^F1 value", gMinus, 0.758035078944662826919680890414, tolLit];
expectApprox["M_s^nat,*", mSNat, 1.66854252965624, tolLit];
expectApprox["M_q^nat,*", mQNat, -0.242696939724365, tolLit];
expectApprox["M_s^comp,*", mSComp, 1.80594111095636, tolLit];
expectApprox["M_q^comp,*", mQComp, -0.451485277739090, tolLit];
(* R1 (was tautological outlet residual): independently reconstruct S_q at piStar
   from the Stage 134 closed-form kernel S(p, kappa) at kappa = Pi/2, and confirm it
   matches the imported Stage 134 literal sQStar via a route that does NOT reuse
   mS = piStar/(1 - rQ sQStar). ASCII-safe names only in comments/strings. *)
kappaQ = Pi/2;
sQRecon = N[piStar*(kappaQ*Tanh[kappaQ]
            + piStar*(Exp[-piStar]/Cosh[kappaQ] - 1))
            / ((1 - Exp[-piStar])*(kappaQ^2 - piStar^2)), 30];
expectApprox["S_q recon from Stage 134 kernel", sQRecon, sQStar, tolLit];
(* Structural sanity only (true by construction, NOT a literal check): *)
Print["outlet form residual nat structural = ", fmt[N[piStar - (mSNat + mQNat*sQStar), 5]]];
Print["outlet form residual comp structural = ", fmt[N[piStar - (mSComp + mQComp*sQStar), 5]]];
expectApprox["R_q^comp - 1/4 (definitional-consistency)", rQComp, 1/4, tolAlg];

Print[""];
Print["Stage 139 Mathematica audit passed."];

Exit[0];
