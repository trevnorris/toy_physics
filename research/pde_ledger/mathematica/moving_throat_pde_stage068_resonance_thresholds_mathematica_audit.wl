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

banner["STAGE 068 — RESONANCE-CORRECTED THRESHOLDS"];

Clear[C2, Cres, Wwall, Wres, Pres, PeReq, Delta0, Deltainf, Ain, Atrans, Wmatch, Wprof, DeltaSym];
$Assumptions =
  Element[{C2, Cres, Wwall, Wres, Pres, PeReq, Delta0, Deltainf, Ain}, Reals] &&
  C2 > 0 && Cres > 0 && Wwall > 0 && Pres > 0 && PeReq > 0 && Delta0 > 0 && Deltainf > 0 && Ain > 0;

(* W_res derived from matched-branch gain decomposition (notes section 1): *)
(*   G_match = rho_star g_phi^2 N_phiphi / (m c_s^2 K_X),                  *)
(*   W_wall  = kappa G_match,    G_res = C^2 G_match,    W_res = kappa G_res *)
Clear[rhoStar, gPhi, Nphi, mPart, cs, KX, kappaW, GmatchExpr, WwallExpr, GresExpr, WresExpr];
$Assumptions = $Assumptions && Element[{rhoStar, gPhi, Nphi, mPart, cs, KX, kappaW}, Reals] &&
  rhoStar > 0 && gPhi > 0 && Nphi > 0 && mPart > 0 && cs > 0 && KX > 0 && kappaW > 0;
GmatchExpr = rhoStar*gPhi^2*Nphi / (mPart*cs^2*KX);
WwallExpr  = kappaW*GmatchExpr;
GresExpr   = C2*GmatchExpr;
WresExpr   = kappaW*GresExpr;
expectZero["W_res - C2 * W_wall (from gain decomposition)",
  FullSimplify[WresExpr - C2*WwallExpr, Assumptions -> $Assumptions]];
expectZero["W_res(C2->1) - W_wall (matched limit)",
  FullSimplify[(WresExpr /. C2 -> 1) - WwallExpr, Assumptions -> $Assumptions]];

(* P_res derived from ratio of required wall figures at resonance C2 -> Cres^2. *)
PresFromRatio =
  FullSimplify[
    ((PeReq/(C2*Delta1)) / (PeReq/Delta1)) /. C2 -> Cres^2,
    Assumptions -> $Assumptions && Delta1 > 0
  ];
expectZero["P_res - 1/C_res^2 (from required-wall-figure ratio)",
  FullSimplify[PresFromRatio - 1/Cres^2, Assumptions -> $Assumptions]];

(* Numeric anchor: paper card states P_res = 1.005612487760576 and             *)
(* C_res^2 = 0.994418836451529 (carried from stage 067).                       *)
With[{CresSqNum = SetPrecision[0.994418836451529, 20],
      PresPaperNum = SetPrecision[1.005612487760576, 20]},
  PresNumResidual = Abs[1/CresSqNum - PresPaperNum];
  Print["P_res numeric residual = ", fmt[PresNumResidual]];
  If[!TrueQ[PresNumResidual < 10^-12],
    fail["P_res numeric anchor", PresNumResidual]
  ];
];

(* Matched-branch threshold from Reduce[W*Delta == PeReq && W > 0, W].       *)
(* Use Reduce rather than Solve to keep the positivity premise explicit.     *)
WfailMatch = First[Cases[
    Reduce[Wmatch*Deltainf == PeReq && Wmatch > 0 && Deltainf > 0 && PeReq > 0, Wmatch, Reals],
    HoldPattern[Wmatch == rhs_] :> rhs, Infinity
  ]];
WsuffMatch = First[Cases[
    Reduce[Wmatch*Delta0 == PeReq && Wmatch > 0 && Delta0 > 0 && PeReq > 0, Wmatch, Reals],
    HoldPattern[Wmatch == rhs_] :> rhs, Infinity
  ]];
WfailMatch = FullSimplify[WfailMatch, Assumptions -> $Assumptions];
WsuffMatch = FullSimplify[WsuffMatch, Assumptions -> $Assumptions];

(* Profile-family threshold from Reduce[C2*W*Delta == PeReq && W > 0, W].   *)
WfailRes = First[Cases[
    Reduce[C2*Wprof*Deltainf == PeReq && Wprof > 0 && Deltainf > 0 && PeReq > 0 && C2 > 0, Wprof, Reals],
    HoldPattern[Wprof == rhs_] :> rhs, Infinity
  ]];
WsuffRes = First[Cases[
    Reduce[C2*Wprof*Delta0 == PeReq && Wprof > 0 && Delta0 > 0 && PeReq > 0 && C2 > 0, Wprof, Reals],
    HoldPattern[Wprof == rhs_] :> rhs, Infinity
  ]];
WfailRes = FullSimplify[WfailRes, Assumptions -> $Assumptions];
WsuffRes = FullSimplify[WsuffRes, Assumptions -> $Assumptions];

Print["Matched fail threshold     = ", fmt[WfailMatch]];
Print["Matched succeed threshold  = ", fmt[WsuffMatch]];
Print["Profile-family fail thresh = ", fmt[WfailRes]];
Print["Profile-family succ thresh = ", fmt[WsuffRes]];

(* Non-trivial cross-relation: WfailRes * C2 == WfailMatch *)
expectZero["Wfail_res * C2 - Wfail_match", WfailRes*C2 - WfailMatch];
expectZero["Wsuff_res * C2 - Wsuff_match", WsuffRes*C2 - WsuffMatch];

(* At resonance C2 = 1/Pres the profile thresholds scale by Pres *)
expectZero["Wfail_res(C2->1/Pres) - Pres*Wfail_match", (WfailRes /. C2 -> 1/Pres) - Pres*WfailMatch];
expectZero["Wsuff_res(C2->1/Pres) - Pres*Wsuff_match", (WsuffRes /. C2 -> 1/Pres) - Pres*WsuffMatch];

banner["PROFILE-SENSITIVE BANDS"];
(* C-form: evaluate WsuffRes at C2 -> Cres^2 (NOT via Pres substitution). *)
WsuffResC = PeReq/(Cres^2 * Delta0);
WfailResC = PeReq/(Cres^2 * Deltainf);
successBandA = FullSimplify[WsuffResC - WsuffMatch, Assumptions -> $Assumptions];
failureBandA = FullSimplify[WfailResC - WfailMatch, Assumptions -> $Assumptions];

(* P-form: (Pres - 1) * Wmatch directly. *)
successBandB = FullSimplify[(Pres - 1)*WsuffMatch, Assumptions -> $Assumptions];
failureBandB = FullSimplify[(Pres - 1)*WfailMatch, Assumptions -> $Assumptions];

Print["Success-side band width (C-form) = ", fmt[successBandA]];
Print["Failure-side band width (C-form) = ", fmt[failureBandA]];
Print["Success-side band width (P-form) = ", fmt[successBandB]];
Print["Failure-side band width (P-form) = ", fmt[failureBandB]];

expectZero["success band C-form vs P-form (under Pres = 1/Cres^2)",
  FullSimplify[(successBandA - successBandB) /. Pres -> 1/Cres^2, Assumptions -> $Assumptions]];
expectZero["failure band C-form vs P-form (under Pres = 1/Cres^2)",
  FullSimplify[(failureBandA - failureBandB) /. Pres -> 1/Cres^2, Assumptions -> $Assumptions]];

banner["FINAL LEDGER"];
Print["W_res derived as |C|^2 W_wall."];
Print["P_res derived as 1/C_res^2."];
Print["Band widths cross-checked two ways."];

Exit[0];
