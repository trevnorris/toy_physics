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

banner["STAGE 121 — NORMALIZED MOUTH-GAIN FAMILY AND COMPENSATION RATIO"];

Clear[lM, thetaSigma, kS, kQ, lam, gS, gQ, rHat, gCore, sigma0];
$Assumptions =
  Element[{lM, thetaSigma, kS, kQ, lam, gS, gQ, rHat, gCore, sigma0}, Reals] &&
  lM > 0 && thetaSigma > 0 && kS > 0 && kQ > 0 && gS > 0 && sigma0 > 0;

mS = FullSimplify[lM*gS^2/(kS*thetaSigma), Assumptions -> $Assumptions];
mQ = FullSimplify[-lM*(kS*gQ - lam*gS)^2/(kS*(kS*kQ + lam^2)*thetaSigma), Assumptions -> $Assumptions];

subs = {lam -> rHat*Sqrt[kS*kQ], gQ -> gCore*gS*Sqrt[kQ/kS]};
mSNorm = FullSimplify[mS /. subs, Assumptions -> $Assumptions];
mQNorm = FullSimplify[mQ /. subs, Assumptions -> $Assumptions];
rQRaw = FullSimplify[-mQNorm/mSNorm, Assumptions -> $Assumptions];
rQ = FullSimplify[rQRaw /. mSNorm -> sigma0, Assumptions -> $Assumptions];

Print["M_s normalized = ", fmt[sigma0]];
Print["M_q normalized = ", fmt[FullSimplify[mQNorm /. mSNorm -> sigma0, Assumptions -> $Assumptions]]];
Print["R_q = ", fmt[rQ]];
expectZero["R_q exact formula", rQ - (gCore - rHat)^2/(1 + rHat^2)];

solPlus = FullSimplify[rQ /. gCore -> rHat + Sqrt[1 + rHat^2]/2, Assumptions -> $Assumptions];
solMinus = FullSimplify[rQ /. gCore -> rHat - Sqrt[1 + rHat^2]/2, Assumptions -> $Assumptions];
Print["R_q on + branch = ", fmt[solPlus]];
Print["R_q on - branch = ", fmt[solMinus]];
expectZero["R_q on + branch - 1/4", solPlus - 1/4];
expectZero["R_q on - branch - 1/4", solMinus - 1/4];

Print[""];
Print["Normalized mouth-gain family verified."];

Exit[0];
