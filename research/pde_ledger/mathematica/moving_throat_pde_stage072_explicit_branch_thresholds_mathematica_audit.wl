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

banner["STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES"];

Clear[chiS, lambdaEll, upsilonW, peReq, rhoW, cSw, v0, hbar];
$Assumptions =
  Element[{chiS, lambdaEll, upsilonW, peReq, rhoW, cSw, v0, hbar}, Reals] &&
  chiS > 0 && lambdaEll > 0 && upsilonW > 0 && peReq > 0 && rhoW > 0 && cSw > 0 && v0 > 0 && hbar > 0;

kappa = FullSimplify[4*chiS^2 + (4/5)*lambdaEll^2, Assumptions -> $Assumptions];
eta = lambdaEll;
alpha = FullSimplify[Sqrt[kappa], Assumptions -> $Assumptions];

delta0 = FullSimplify[
  eta*(Cosh[alpha] - 1)/(alpha^2*(alpha*Sinh[alpha] + eta*Cosh[alpha])),
  Assumptions -> $Assumptions
];
deltaInf = FullSimplify[
  (Cosh[alpha] + (eta/alpha)*Sinh[alpha] - 1)/(alpha*Sinh[alpha] + eta*Cosh[alpha]),
  Assumptions -> $Assumptions
];

upsilonFail = FullSimplify[peReq/(lambdaEll^2*deltaInf), Assumptions -> $Assumptions];
upsilonSuff = FullSimplify[peReq/(lambdaEll^2*delta0), Assumptions -> $Assumptions];

Print["kappa = ", fmt[kappa]];
Print["eta = ", fmt[eta]];
Print["Delta_0 = ", fmt[delta0]];
Print["Delta_inf = ", fmt[deltaInf]];
Print["Upsilon_fail = ", fmt[upsilonFail]];
Print["Upsilon_suff = ", fmt[upsilonSuff]];

v0FailSq = FullSimplify[hbar^2*cSw^2*upsilonFail/(4*rhoW^2), Assumptions -> $Assumptions];
v0SuffSq = FullSimplify[hbar^2*cSw^2*upsilonSuff/(4*rhoW^2), Assumptions -> $Assumptions];
Print["V0_fail^2 = ", fmt[v0FailSq]];
Print["V0_suff^2 = ", fmt[v0SuffSq]];

banner["SHELL-GRADIENT DOMINATED ASYMPTOTICS"];

cShell = 2/Sqrt[5];
upsilonFailShell = FullSimplify[2*peReq/(Sqrt[5]*lambdaEll), Assumptions -> $Assumptions];
upsilonSuffShell = FullSimplify[(4/5)*(1 + 2/Sqrt[5])*peReq, Assumptions -> $Assumptions];
delta0Shell = FullSimplify[lambdaEll/((cShell*lambdaEll)^2*(cShell*lambdaEll + lambdaEll)), Assumptions -> $Assumptions];
deltaInfShell = FullSimplify[(1 + lambdaEll/(cShell*lambdaEll))/(cShell*lambdaEll + lambdaEll), Assumptions -> $Assumptions];

(* Extract the leading-order asymptotic of the full Delta_0, Delta_inf as
   chi_s -> 0 then Lambda_ell -> infinity, and confirm the ratio to the
   hand-built shell forms tends to 1. *)
delta0ShellRatio = Limit[FullSimplify[(delta0 /. chiS -> 0)/delta0Shell,
    Assumptions -> lambdaEll > 0], lambdaEll -> Infinity];
deltaInfShellRatio = Limit[FullSimplify[(deltaInf /. chiS -> 0)/deltaInfShell,
    Assumptions -> lambdaEll > 0], lambdaEll -> Infinity];
Print["Delta0  shell leading-order ratio = ", fmt[delta0ShellRatio]];
Print["DeltaInf shell leading-order ratio = ", fmt[deltaInfShellRatio]];
expectZero["Delta0  shell leading-order matches full delta0",
  delta0ShellRatio - 1];
expectZero["DeltaInf shell leading-order matches full deltaInf",
  deltaInfShellRatio - 1];

Print["Upsilon_fail_shell = ", fmt[upsilonFailShell]];
Print["Upsilon_suff_shell = ", fmt[upsilonSuffShell]];
expectZero["shell fail asymptotic", peReq/(lambdaEll^2*deltaInfShell) - upsilonFailShell];
expectZero["shell suff asymptotic", peReq/(lambdaEll^2*delta0Shell) - upsilonSuffShell];

banner["COMPRESSION DOMINATED ASYMPTOTICS"];

delta0Comp = FullSimplify[lambdaEll/((2*chiS)^2*(2*chiS + lambdaEll)), Assumptions -> $Assumptions];
deltaInfComp = FullSimplify[(1 + lambdaEll/(2*chiS))/(2*chiS + lambdaEll), Assumptions -> $Assumptions];
upsilonFailComp = FullSimplify[peReq/(lambdaEll^2*deltaInfComp), Assumptions -> $Assumptions];
upsilonSuffComp = FullSimplify[peReq/(lambdaEll^2*delta0Comp), Assumptions -> $Assumptions];

(* Confirm the hand-built compression forms are the leading-order asymptotic of
   the full Delta_0, Delta_inf as chi_s -> infinity at fixed Lambda_ell. *)
delta0CompRatio = Limit[FullSimplify[delta0/delta0Comp,
    Assumptions -> chiS > 0 && lambdaEll > 0], chiS -> Infinity];
deltaInfCompRatio = Limit[FullSimplify[deltaInf/deltaInfComp,
    Assumptions -> chiS > 0 && lambdaEll > 0], chiS -> Infinity];
Print["Delta0  comp leading-order ratio = ", fmt[delta0CompRatio]];
Print["DeltaInf comp leading-order ratio = ", fmt[deltaInfCompRatio]];
expectZero["Delta0  comp leading-order matches full delta0",
  delta0CompRatio - 1];
expectZero["DeltaInf comp leading-order matches full deltaInf",
  deltaInfCompRatio - 1];

Print["Upsilon_fail_comp = ", fmt[upsilonFailComp]];
Print["Upsilon_suff_comp = ", fmt[upsilonSuffComp]];
expectZero["compression fail asymptotic", upsilonFailComp - 2*peReq*chiS/lambdaEll^2];
expectZero[
  "compression suff asymptotic",
  upsilonSuffComp - 4*peReq*chiS^2*(lambdaEll + 2*chiS)/lambdaEll^3
];

Print[""];
Print["Stage 072 Mathematica audit passed."];

Exit[0];
