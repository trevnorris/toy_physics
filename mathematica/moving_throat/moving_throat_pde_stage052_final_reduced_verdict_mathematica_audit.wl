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

banner["STAGE 052 — FINAL REDUCED SUPPORT/SOURCE VERDICT"];

Clear[PeReq, Delta0, Deltainf, Pres, Wwall];
$Assumptions =
  Element[{PeReq, Delta0, Deltainf, Pres, Wwall}, Reals] &&
  PeReq > 0 && Delta0 > 0 && Deltainf > 0 && Pres > 0 && Wwall > 0;

WfailMatch = FullSimplify[PeReq/Deltainf, Assumptions -> $Assumptions];
WsuffMatch = FullSimplify[PeReq/Delta0, Assumptions -> $Assumptions];
WfailRes = FullSimplify[Pres*WfailMatch, Assumptions -> $Assumptions];
WsuffRes = FullSimplify[Pres*WsuffMatch, Assumptions -> $Assumptions];

Print["Universal matched fail threshold    = ", fmt[WfailMatch]];
Print["Universal matched success threshold = ", fmt[WsuffMatch]];
Print["Resonance-family fail threshold     = ", fmt[WfailRes]];
Print["Resonance-family success threshold  = ", fmt[WsuffRes]];

deltaFail = FullSimplify[WfailRes - WfailMatch, Assumptions -> $Assumptions];
deltaSuff = FullSimplify[WsuffRes - WsuffMatch, Assumptions -> $Assumptions];
Print["Failure-side profile-sensitive width = ", fmt[deltaFail]];
Print["Success-side profile-sensitive width = ", fmt[deltaSuff]];

expectZero["failure relative width", deltaFail/WfailMatch - (Pres - 1)];
expectZero["success relative width", deltaSuff/WsuffMatch - (Pres - 1)];

Print["Success if W_wall - P_res*Pe_req/Delta_0 >= 0:"];
Print["  ", fmt[Wwall - WsuffRes]];
Print["Failure if Pe_req/Delta_inf - W_wall >= 0:"];
Print["  ", fmt[WfailMatch - Wwall]];

banner["FINAL LEDGER"];
Print["Exact reduced verdict:"];
Print["  Universal fail   : W_wall <= Pe_req / Delta_inf"];
Print["  Universal succeed: W_wall >= Pe_req / Delta_0"];
Print["  Resonance-family thresholds are shifted only by P_res."];
Print["Therefore profile mismatch matters only in the narrow O(P_res - 1) bands"];
Print["around the universal matched-branch thresholds."];

Exit[0];
