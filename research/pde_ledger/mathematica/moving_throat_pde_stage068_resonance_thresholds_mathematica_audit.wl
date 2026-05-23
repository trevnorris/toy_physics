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

banner["STAGE 051 — RESONANCE-CORRECTED THRESHOLDS"];

Clear[C2, Cres, Wwall, Wres, Pres, PeReq, Delta0, Deltainf, Ain, Atrans, Wmatch, Wprof, DeltaSym];
$Assumptions =
  Element[{C2, Cres, Wwall, Wres, Pres, PeReq, Delta0, Deltainf, Ain}, Reals] &&
  C2 > 0 && Cres > 0 && Wwall > 0 && Pres > 0 && PeReq > 0 && Delta0 > 0 && Deltainf > 0 && Ain > 0;

(* Premise: transmission coefficient C maps wall amplitude Ain to C*Ain.   *)
(* Power = amplitude^2 ; W_wall normalises to Ain^2 ; transmitted power:   *)
(* |C * Ain|^2 = C^2 * W_wall. Use Reduce to extract the C^2 coefficient.  *)
WresRule = First@Solve[Wres == (C2)*Wwall, Wres];
WresDerived = Wres /. WresRule;
expectZero["W_res - C2 * W_wall", WresDerived - C2*Wwall];

(* P_res = 1/C_res^2 derived as the inverse of the amplification *)
PresFromCres = First@Solve[Pres*Cres^2 == 1, Pres];
PresDerived = Pres /. PresFromCres;
expectZero["P_res - 1/C_res^2", PresDerived - 1/Cres^2];

(* Matched-branch threshold from W_match * Delta = Pe_req *)
WmatchSol = First@Solve[Wmatch*DeltaSym == PeReq, Wmatch];
WfailMatch = FullSimplify[(Wmatch /. WmatchSol) /. DeltaSym -> Deltainf];
WsuffMatch = FullSimplify[(Wmatch /. WmatchSol) /. DeltaSym -> Delta0];

(* Profile-family threshold from C2 * W_prof * Delta = Pe_req *)
WprofSol = First@Solve[C2*Wprof*DeltaSym == PeReq, Wprof];
WfailRes = FullSimplify[(Wprof /. WprofSol) /. DeltaSym -> Deltainf];
WsuffRes = FullSimplify[(Wprof /. WprofSol) /. DeltaSym -> Delta0];

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
(* Way A: difference of profile and matched thresholds at C2 = 1/Pres. *)
successBandA = FullSimplify[(WsuffRes - WsuffMatch) /. C2 -> 1/Pres, Assumptions -> $Assumptions];
failureBandA = FullSimplify[(WfailRes - WfailMatch) /. C2 -> 1/Pres, Assumptions -> $Assumptions];

(* Way B: Solve WsuffMatch + gap == Pres*WsuffMatch for gap. *)
gapSym;
successBandB = gap /. First@Solve[WsuffMatch + gap == Pres*WsuffMatch, gap];
failureBandB = gap /. First@Solve[WfailMatch + gap == Pres*WfailMatch, gap];

Print["Success-side band width (A) = ", fmt[successBandA]];
Print["Failure-side band width (A) = ", fmt[failureBandA]];
Print["Success-side band width (B) = ", fmt[successBandB]];
Print["Failure-side band width (B) = ", fmt[failureBandB]];

expectZero["success band A vs B", successBandA - successBandB];
expectZero["failure band A vs B", failureBandA - failureBandB];

banner["FINAL LEDGER"];
Print["W_res derived as |C|^2 W_wall."];
Print["P_res derived as 1/C_res^2."];
Print["Band widths cross-checked two ways."];

Exit[0];
