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
  If[! MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res == 0], pass[name], fail[name, res]];
];

expectPositive[name_String, expr_] := Module[{shown, res},
  shown = FullSimplify[expr, Assumptions -> $Assumptions];
  res = FullSimplify[expr > 0, Assumptions -> $Assumptions];
  Print[name, " > 0  -> ", fmt[shown]];
  If[TrueQ[res], pass[name], fail[name, shown]];
];

banner["STAGE 052 — FINAL REDUCED SUPPORT/SOURCE VERDICT"];

Clear[PeReq, Delta0, DeltaGap, PresGap, Deltainf, Pres, Cres2, uFail, uSucc, vFail, vSucc];
$Assumptions =
  Element[{PeReq, Delta0, DeltaGap, PresGap, vFail, vSucc}, Reals] &&
  PeReq > 0 && Delta0 > 0 && DeltaGap > 0 && PresGap > 0 &&
  vFail > 0 && vSucc > 0;

Deltainf = FullSimplify[Delta0 + DeltaGap, Assumptions -> $Assumptions];
Pres = FullSimplify[1 + PresGap, Assumptions -> $Assumptions];
Cres2 = FullSimplify[1/Pres, Assumptions -> $Assumptions];
uFail = FullSimplify[vFail/(1 + vFail), Assumptions -> $Assumptions];
uSucc = FullSimplify[vSucc/(1 + vSucc), Assumptions -> $Assumptions];

WfailMatch = FullSimplify[PeReq/Deltainf, Assumptions -> $Assumptions];
WsuffMatch = FullSimplify[PeReq/Delta0, Assumptions -> $Assumptions];
WfailRes = FullSimplify[Pres WfailMatch, Assumptions -> $Assumptions];
WsuffRes = FullSimplify[Pres WsuffMatch, Assumptions -> $Assumptions];

Print["Matched-branch fail threshold       = ", fmt[WfailMatch]];
Print["Matched-branch success threshold    = ", fmt[WsuffMatch]];
Print["Resonance-family fail threshold     = ", fmt[WfailRes]];
Print["Resonance-family success threshold  = ", fmt[WsuffRes]];

expectZero["P_res - 1/C_res^2", Pres - 1/Cres2];
expectZero[
  "matched window width",
  WsuffMatch - WfailMatch - PeReq (Deltainf - Delta0)/(Delta0 Deltainf)
];
expectPositive["Delta_inf - Delta_0", Deltainf - Delta0];
expectPositive["matched success threshold - matched fail threshold", WsuffMatch - WfailMatch];
expectZero[
  "resonance fail threshold - Pe_req/(C_res^2 Delta_inf)",
  WfailRes - PeReq/(Cres2 Deltainf)
];
expectZero[
  "resonance success threshold - Pe_req/(C_res^2 Delta_0)",
  WsuffRes - PeReq/(Cres2 Delta0)
];
expectPositive["1 - C_res^2", 1 - Cres2];
expectPositive["P_res - 1", Pres - 1];

deltaFail = FullSimplify[WfailRes - WfailMatch, Assumptions -> $Assumptions];
deltaSucc = FullSimplify[WsuffRes - WsuffMatch, Assumptions -> $Assumptions];

Print["Failure-side profile-sensitive width = ", fmt[deltaFail]];
Print["Success-side profile-sensitive width = ", fmt[deltaSucc]];

expectZero[
  "failure-side width - Pe_req(1-C_res^2)/(C_res^2 Delta_inf)",
  deltaFail - PeReq (1 - Cres2)/(Cres2 Deltainf)
];
expectZero[
  "success-side width - Pe_req(1-C_res^2)/(C_res^2 Delta_0)",
  deltaSucc - PeReq (1 - Cres2)/(Cres2 Delta0)
];
expectZero[
  "P_res - 1 - (1-C_res^2)/C_res^2",
  (Pres - 1) - (1 - Cres2)/Cres2
];
expectPositive["failure-side width", deltaFail];
expectPositive["success-side width", deltaSucc];

WfailureBand = FullSimplify[WfailMatch + uFail deltaFail, Assumptions -> $Assumptions];
WsuccessBand = FullSimplify[WsuffMatch + uSucc deltaSucc, Assumptions -> $Assumptions];

expectPositive["failure-band point - matched fail edge", WfailureBand - WfailMatch];
expectPositive["resonance fail edge - failure-band point", WfailRes - WfailureBand];
expectPositive["success-band point - matched success edge", WsuccessBand - WsuffMatch];
expectPositive["resonance success edge - success-band point", WsuffRes - WsuccessBand];

banner["FINAL LEDGER"];
Print["Stage-49 matched branch:"];
Print["  Universal fail   : W_wall <= Pe_req / Delta_inf"];
Print["  Universal succeed: W_wall >= Pe_req / Delta_0"];
Print["Stage-51 resonance family:"];
Print["  Fail threshold   : W_wall <= Pe_req / (C_res^2 Delta_inf)"];
Print["  Success threshold: W_wall >= Pe_req / (C_res^2 Delta_0)"];
Print["Therefore the only profile-sensitive regions are the exact side-bands"];
Print["  Pe_req/Delta_inf < W_wall < Pe_req/(C_res^2 Delta_inf)"];
Print["and"];
Print["  Pe_req/Delta_0 < W_wall < Pe_req/(C_res^2 Delta_0),"];
Print["whose relative widths are both P_res - 1 = (1 - C_res^2)/C_res^2."];

Exit[0];
