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

expectTrue[name_String, cond_] := Module[{res},
  res = TrueQ[cond];
  Print[name, " = ", fmt[res]];
  If[res, pass[name], fail[name, cond]];
];

banner["STAGE 072 — EXPLICIT FAMILY-1 VERDICT FOR THE MINIMAL ISOTROPIC BRANCH"];

rhoMin = 4/3;
zetaMin = 1/3;
rhoSuff = SetPrecision[3.46622291347846, 20];
rhoFail = SetPrecision[3.46752913273870, 20];
rhoMax = SetPrecision[3.46752922945601, 20];
zetaMax = SetPrecision[2.46752922945601, 20];
aF1 = SetPrecision[1.00005192880220, 20];

deltaSuff = N[rhoSuff - rhoMin, 25];
deltaFail = N[rhoFail - rhoMin, 25];
deltaMax = N[rhoMax - rhoMin, 25];
deltaZeta = N[zetaMax - zetaMin, 25];
deltaAF1 = N[aF1 - zetaMin, 25];

Print["rho_min = ", fmt[N[rhoMin, 25]]];
Print["rho_suff = ", fmt[rhoSuff]];
Print["rho_fail = ", fmt[rhoFail]];
Print["rho_max = ", fmt[rhoMax]];
Print["zeta_min = ", fmt[N[zetaMin, 25]]];
Print["zeta_max = ", fmt[zetaMax]];
Print["A_F1 = ", fmt[aF1]];

Print["Delta_suff = ", fmt[deltaSuff]];
Print["Delta_fail = ", fmt[deltaFail]];
Print["Delta_max = ", fmt[deltaMax]];
Print["Delta_zeta = ", fmt[deltaZeta]];
Print["Delta_AF1 = ", fmt[deltaAF1]];

expectTrue["Family-1 loading-ratio ordering", rhoMin < rhoSuff < rhoFail < rhoMax];
expectTrue["minimal isotropic branch stays in the symmetric-lowest-twin regime", zetaMin < 1];
expectTrue["minimal isotropic branch succeeds at zero transport bias", zetaMin < aF1];
expectTrue["minimal isotropic branch stays below the Family-1 ceiling", zetaMin < zetaMax];

Print[""];
Print["Stage 072 Mathematica audit passed."];

Exit[0];
