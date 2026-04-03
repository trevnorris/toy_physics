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

banner["STAGE 122 — ACTUAL FAMILY-1 MOUTH GAINS"];

rF = N[Sqrt[(12/Pi^2)*(37/20)^2 - 1], 30];
piStar = SetPrecision[1.50882951349316, 30];
sQStar = SetPrecision[0.658075937605429, 30];

rQNat = N[(1 - rF)^2/(1 + rF^2), 30];
mSNat = N[piStar/(1 - rQNat*sQStar), 30];
mQNat = N[-rQNat*mSNat, 30];

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

expectApprox["R_q^comp - 1/4", rQComp, 1/4, 10^-25];

Print[""];
Print["Stage 122 Mathematica audit passed."];

Exit[0];
