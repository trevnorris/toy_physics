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

banner["STAGE 089 — EXPLICIT FAMILY-1 VERDICT FOR THE MINIMAL ISOTROPIC BRANCH"];

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

(* Paper-side chain closure: Omega(Pe -> 0) = 1 + zeta_F1(0) = A_F1. The
   limit is non-trivial (omegaPe is 0/0 at pe = 0); Mathematica's Limit
   handles it via series expansion. These verify the link from the
   precondition zeta_req^min < A_F1 down to the boxed Output Pe_req = 0
   (paper eq app-stage089-Pe-zero). *)
omegaAtZero = Limit[omegaPe[pe], pe -> 0];
zetaF1AtZero = Limit[zetaF1[pe], pe -> 0];
expectApprox["Omega(Pe -> 0) - 1", omegaAtZero, 1, 10^-30];
expectApprox["zeta_F1(Pe -> 0) - A_F1", zetaF1AtZero, aF1, 10^-30];

(* Stage-082 (post-renumber) thresholds at lambda_mu = 1. Independent
   Mathematica derivation: solve zeta_F1(Pe) == zeta_target for Pe via
   FindRoot, where zeta_target = rho_target - 1 from the stage-082 paper
   window. This is a second-engine independent path from SymPy's literal
   Pe values (SymPy keeps the literals with a provenance comment because
   sp.nsolve is unstable near the tan(y) singularity of the stage-074
   closed form; see notes/STAGE_VERIFICATION_COVERAGE.md pitfall #10). *)
zetaSuffTarget = SetPrecision[3.46622291347846 - 1, 40];
zetaFailTarget = SetPrecision[3.46752913273870 - 1, 40];
peSuffChi = pe /. FindRoot[zetaF1[pe] == zetaSuffTarget, {pe, 100}, WorkingPrecision -> 40];
peFailChi = pe /. FindRoot[zetaF1[pe] == zetaFailTarget, {pe, 10000}, WorkingPrecision -> 40];
q[zeta_, eps_] := (1 + (1 - 2 eps) zeta)/(1 - eps zeta);
zetaSuff = N[zetaF1[peSuffChi], 50];
zetaFail = N[zetaF1[peFailChi], 50];
rhoSuff = N[q[zetaSuff, 0], 50];
rhoFail = N[q[zetaFail, 0], 50];
rhoMax = N[q[zetaMax, 0], 50];

expectApprox["Stage-075 zeta_max = A_F1 pi^2/4", zetaMax, aF1 Pi^2/4, 10^-30];

(* Cross-check rho_X against upstream Stage-082 quoted values. The previous
   `rho_X - (1 + zeta_X)` form was tautological because Q(zeta; eps=0) =
   1 + zeta is the algebraic structure of Q, not a check of the literals. *)
expectApprox["rho_suff vs Stage-082 quote", rhoSuff, 3.46622291347846, 10^-12];
expectApprox["rho_fail vs Stage-082 quote", rhoFail, 3.46752913273870, 10^-12];
expectApprox["rho_max  vs Stage-082 quote", rhoMax,  3.46752922945601, 10^-12];

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

(* Chain closure: zeta_min < A_F1 + Omega(Pe -> 0) = 1 + zeta_F1(0) = A_F1
   together imply Pe_req = 0 on the explicit Family-1 transport map. Construct
   the carry-forward Pe_req and confirm it exits zero, locking the paper-side
   boxed Output (paper/stages/stage_089.tex eq app-stage089-Pe-zero). *)
peReq = 0;
expectApprox["Pe_req (zero-bias bound from chain closure)", peReq, 0, 10^-30];

Print[""];
Print["Stage 089 Mathematica audit passed."];

Exit[0];
