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

Clear[C2, Wwall, Wres, Pres, PeReq, Delta0, Deltainf];
$Assumptions =
  Element[{C2, Wwall, Wres, Pres, PeReq, Delta0, Deltainf}, Reals] &&
  C2 > 0 && Wwall > 0 && Wres > 0 && Pres > 0 && PeReq > 0 && Delta0 > 0 && Deltainf > 0;

Wres = FullSimplify[C2*Wwall, Assumptions -> $Assumptions];
Print["W_res = ", fmt[Wres]];

WfailMatch = FullSimplify[PeReq/Deltainf, Assumptions -> $Assumptions];
WsuffMatch = FullSimplify[PeReq/Delta0, Assumptions -> $Assumptions];
WfailRes = FullSimplify[PeReq/(C2*Deltainf), Assumptions -> $Assumptions];
WsuffRes = FullSimplify[PeReq/(C2*Delta0), Assumptions -> $Assumptions];

Print["Matched fail threshold     = ", fmt[WfailMatch]];
Print["Matched succeed threshold  = ", fmt[WsuffMatch]];
Print["Profile-family fail thresh = ", fmt[WfailRes]];
Print["Profile-family succ thresh = ", fmt[WsuffRes]];

expectZero["Wfail_res - P_res*Wfail_match", (WfailRes /. C2 -> 1/Pres) - Pres*WfailMatch];
expectZero["Wsuff_res - P_res*Wsuff_match", (WsuffRes /. C2 -> 1/Pres) - Pres*WsuffMatch];

banner["PROFILE-SENSITIVE BANDS"];
successBand = FullSimplify[WsuffRes - WsuffMatch, Assumptions -> $Assumptions];
failureBand = FullSimplify[WfailRes - WfailMatch, Assumptions -> $Assumptions];
Print["Success-side band width = ", fmt[successBand]];
Print["Failure-side band width = ", fmt[failureBand]];

successBandRes = FullSimplify[successBand /. C2 -> 1/Pres, Assumptions -> $Assumptions];
failureBandRes = FullSimplify[failureBand /. C2 -> 1/Pres, Assumptions -> $Assumptions];
expectZero["success width / matched threshold", successBandRes/WsuffMatch - (Pres - 1)];
expectZero["failure width / matched threshold", failureBandRes/WfailMatch - (Pres - 1)];

banner["FINAL LEDGER"];
Print["Exact formulas:"];
Print["  W_res(r) = C^2(r) W_wall"];
Print["  W_fail^(res) = Pe_req / [ C^2(r) Delta_inf ]"];
Print["  W_suff^(res) = Pe_req / [ C^2(r) Delta_0 ]"];
Print["At the resonance point C^2 = C_res^2, these become:"];
Print["  W_fail^(res,*) = P_res Pe_req / Delta_inf"];
Print["  W_suff^(res,*) = P_res Pe_req / Delta_0"];
Print["with P_res = 1/C_res^2."];
Print[""];
Print["Interpretation:"];
Print["  The explicit independent-profile benchmark modifies the matched-layer"];
Print["  threshold window only by the multiplicative factor P_res."];

Exit[0];
