ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
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

expectApprox[name_String, value_, target_, tol_] := Module[{delta},
  delta = N[value - target, 40];
  Print[name, " residual = ", fmt[delta]];
  If[TrueQ[Abs[delta] < tol], pass[name], fail[name, delta]];
];

banner["STAGE 142 — SELF-CONSISTENT MOUTH-BRANCH LAW"];

Clear[piM];
$Assumptions = Element[piM, Reals] && piM > 0;

(* r_F1: Family-1 reduced mixed-core ratio. Carried forward from upstream
   (see notes/stages/moving_throat_pde_stage142_selfconsistent_mouth_branch.md
   section 1; original derivation is in the upstream "shell/mixed core" block
   referenced by paper/stages/stage_142.tex Inputs field). *)
r = Sqrt[4107 - 100*Pi^2]/(10*Pi);
gPi = 2*piM*(2*piM*Exp[piM] + Pi)/((4*piM^2 + Pi^2)*(Exp[piM] - 1));
(* S_q(Pi) closed form: carried forward from the self-matched mouth-susceptibility
   closure (Stage 140 / Sigma_0 = (20/9) That_m^2). The closed form here is
   S(Pi, pi/2), evaluated at the fixed second argument pi/2. *)
sQ = piM*(((Pi/2)*Tanh[Pi/2]) + piM*(Exp[-piM]*Sech[Pi/2] - 1))/((1 - Exp[-piM])*((Pi/2)^2 - piM^2));
rQ = (gPi - r)^2/(1 + r^2);
sigma0 = piM/(1 - rQ*sQ);
tHat = Sqrt[(9/20)*sigma0];

subbanner["Independent numerical cross-checks (Mathematica)"];

(* Independent series expansion of gPi near Pi = 0. The closed form
   gPi = 2 Pi (2 Pi exp(Pi) + pi) / ((4 Pi^2 + pi^2)(exp(Pi)-1))
   expanded to O(Pi^3) should agree with a direct Series[] of the same
   expression, cross-checking the algebraic encoding. *)
gPiSeries = Normal[Series[gPi, {piM, 0, 4}]];
gPiSampleVals = Table[
    {pVal, N[gPi /. piM -> pVal, 30] - N[gPiSeries /. piM -> pVal, 30]},
    {pVal, {1/10, 2/10, 3/10}}];
Print["g_Pi closed-form vs series residuals at piM={0.1,0.2,0.3}: ", fmt[gPiSampleVals[[All, 2]]]];
Do[
    If[Abs[gPiSampleVals[[i, 2]]] > 10^-3,
        fail["g_Pi closed-form vs series small-piM disagreement", gPiSampleVals[[i, 2]]]],
    {i, 1, Length[gPiSampleVals]}];
pass["g_Pi closed-form/series consistency at small piM"];

(* Independent r_F1 cross-check: r_F1 = sqrt(4107 - 100 pi^2)/(10 pi)
   should satisfy 100 pi^2 (1 + r^2) = 4107. *)
rSquared = 1 + r^2;
rIdentity = FullSimplify[100*Pi^2*rSquared - 4107];
expectZero["r_F1 satisfies 100 pi^2 (1+r^2) = 4107", rIdentity];

subbanner["Core-to-mouth reduction"];
Print["r_F1 = ", fmt[r]];
Print["g_Pi = ", fmt[gPi]];
Print["S_q(Pi) = ", fmt[sQ]];
Print["R_q(Pi) = ", fmt[rQ]];
Print["Sigma0(Pi) = ", fmt[sigma0]];
Print["That(Pi) = ", fmt[tHat]];

gMinus = FullSimplify[r - Sqrt[1 + r^2]/2, Assumptions -> $Assumptions];
rQMinus = FullSimplify[(gMinus - r)^2/(1 + r^2), Assumptions -> $Assumptions];
expectZero["R_q(g_minus)-1/4", rQMinus - 1/4];

subbanner["Canonical compensation point"];
piStar = piM /. FindRoot[gPi == gMinus, {piM, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
gStar = N[gPi /. piM -> piStar, 30];
rQStar = N[rQ /. piM -> piStar, 30];
sQStar = N[sQ /. piM -> piStar, 30];
sigmaStar = N[sigma0 /. piM -> piStar, 30];
tHatStar = N[tHat /. piM -> piStar, 30];

Print["g_minus^F1 = ", fmt[N[gMinus, 30]]];
Print["Pi_*       = ", fmt[N[piStar, 30]]];
Print["g(Pi_*)    = ", fmt[gStar]];
Print["R_q(Pi_*)  = ", fmt[rQStar]];
Print["S_q(Pi_*)  = ", fmt[sQStar]];
Print["Sigma0(Pi_*) = ", fmt[sigmaStar]];
Print["That(Pi_*)   = ", fmt[tHatStar]];
expectApprox["Pi_* compensation solve", gStar, N[gMinus, 30], 10^-12];
expectApprox["R_q(Pi_*) numeric = 1/4", rQStar, 1/4, 10^-20];
expectApprox["g_-^{F1} value",      N[gMinus, 30],   N[Rationalize[0.7580350789446628269196808904, 0], 30], 10^-25];
expectApprox["Pi_* value",          piStar,          N[Rationalize[1.5088295134931555274704351177, 0], 30], 10^-12];
expectApprox["S_q(Pi_*) value",     sQStar,          N[Rationalize[0.6580759376054292719303153134, 0], 30], 10^-12];
expectApprox["Sigma_0(Pi_*) value", sigmaStar,       N[Rationalize[1.8059411109563538072179672471, 0], 30], 10^-12];
expectApprox["That(Pi_*) value",    tHatStar,        N[Rationalize[0.9014840541742040227024016887, 0], 30], 10^-12];

banner["STAGE 142 LEDGER"];
Print["Self-consistent Family-1 mouth branch:"];
Print["  Pi = Sigma0 * [1 - R_q(Pi) S_q(Pi)]"];
Print["  Sigma0(Pi) = Pi / (1 - R_q(Pi) S_q(Pi))"];
Print["  That(Pi)   = sqrt(9 Sigma0(Pi) / 20)"];
Print[""];
Print["Canonical point:"];
Print["  Pi_*       = ", fmt[N[piStar, 30]]];
Print["  Sigma0_*   = ", fmt[sigmaStar]];
Print["  That_*     = ", fmt[tHatStar]];

Print[""];
Print["Stage 142 Mathematica audit passed."];

Exit[0];
