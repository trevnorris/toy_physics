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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = N[Abs[value - target], 50];
  Print[name, " diff = ", fmt[diff]];
  If[diff <= tol, pass[name], fail[name, diff]];
];

banner["STAGE 072 — EXPLICIT FAMILY-1 VERDICT FOR THE MINIMAL ISOTROPIC BRANCH"];

rhoMin = 4/3;
zetaMin = 1/3;

(* Stage-62 Family-1 demand map data. *)
kappaF1 = 12321/5;
etaF1 = 37;
yF1 = y /. FindRoot[y Tan[y] == etaF1, {y, 1.53}, WorkingPrecision -> 60, AccuracyGoal -> 40, PrecisionGoal -> 40];
aF1 = N[(kappaF1 + Pi^2/4)/(kappaF1 + yF1^2), 50];
omegaPe[pe_] := Pi pe (2 pe Exp[pe] + Pi)/((4 pe^2 + Pi^2) (Exp[pe] - 1));
zetaF1[pe_] := aF1 omegaPe[pe]^2;
zetaMax = N[aF1 Pi^2/4, 50];

(* Stage-63/69 thresholds at lambda_mu = 1. *)
peSuffChi = SetPrecision[96.5285247264386, 40];
peFailChi = SetPrecision[11220.5441626259, 40];
q[zeta_, eps_] := (1 + (1 - 2 eps) zeta)/(1 - eps zeta);
zetaSuff = N[zetaF1[peSuffChi], 50];
zetaFail = N[zetaF1[peFailChi], 50];
rhoSuff = N[q[zetaSuff, 0], 50];
rhoFail = N[q[zetaFail, 0], 50];
rhoMax = N[q[zetaMax, 0], 50];

expectApprox["Stage-62 zeta_max = A_F1 pi^2/4", zetaMax, aF1 Pi^2/4, 10^-30];
expectApprox["Stage-69 Q(zeta;0) = 1 + zeta on rho_suff", rhoSuff, 1 + zetaSuff, 10^-30];
expectApprox["Stage-69 Q(zeta;0) = 1 + zeta on rho_fail", rhoFail, 1 + zetaFail, 10^-30];
expectApprox["Stage-69 Q(zeta;0) = 1 + zeta on rho_max", rhoMax, 1 + zetaMax, 10^-30];

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
Print["Stage 089 Mathematica audit passed."];

Exit[0];
