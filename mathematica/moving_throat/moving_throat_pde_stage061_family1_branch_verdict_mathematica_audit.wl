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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 40] - N[target, 40]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, cond]];
];

banner["STAGE 061 — FAMILY-1 BRANCH VERDICT"];

Clear[lambdaMu, peReq];
$Assumptions = Element[{lambdaMu, peReq}, Reals] && lambdaMu > 0 && peReq > 0;

thetaChiCoeff = SetPrecision[4.06863235008162, 20];
thetaJCoeff = SetPrecision[0.927552032539308, 20];
thetaFailCoeff = SetPrecision[3.62605617972939*^-4, 20];
thetaSuffCoeff = SetPrecision[4.21495341569977*^-2, 20];

thetaChi = thetaChiCoeff*lambdaMu^2;
thetaJ = thetaJCoeff*lambdaMu^2;
thetaFail = thetaFailCoeff*peReq;
thetaSuff = thetaSuffCoeff*peReq;

peSuffChi = N[thetaChiCoeff/thetaSuffCoeff, 30];
peFailChi = N[thetaChiCoeff/thetaFailCoeff, 30];
peSuffJ = N[thetaJCoeff/thetaSuffCoeff, 30];
peFailJ = N[thetaJCoeff/thetaFailCoeff, 30];

Print["Pe_suff^(chi) / lambda_mu^2 = ", fmt[peSuffChi]];
Print["Pe_fail^(chi) / lambda_mu^2 = ", fmt[peFailChi]];
Print["Pe_suff^(J) / lambda_mu^2 = ", fmt[peSuffJ]];
Print["Pe_fail^(J) / lambda_mu^2 = ", fmt[peFailJ]];

expectApprox["Pe_suff^(chi) numeric check", peSuffChi, SetPrecision[96.528524726438575954, 25], 10^-12];
expectApprox["Pe_fail^(chi) numeric check", peFailChi, SetPrecision[11220.544162625905301, 25], 10^-9];
expectApprox["Pe_suff^(J) numeric check", peSuffJ, SetPrecision[22.006222633075413597, 25], 10^-12];
expectApprox["Pe_fail^(J) numeric check", peFailJ, SetPrecision[2558.0189234920526360, 25], 10^-10];
expectTrue["Pe_suff^(chi) < Pe_fail^(chi)", peSuffChi < peFailChi];
expectTrue["Pe_suff^(J) < Pe_fail^(J)", peSuffJ < peFailJ];

Print[""];
Print["Stage 061 Mathematica audit passed."];

Exit[0];
