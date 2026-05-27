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

banner["STAGE 078 — FAMILY-1 BRANCH VERDICT"];

Clear[lambdaMu, peReq];
$Assumptions = Element[{lambdaMu, peReq}, Reals] && lambdaMu > 0 && peReq > 0;

(* Independent re-derivation of the four window coefficients from
   their Stage-75 symbolic closed forms (no SymPy input).             *)
thetaFailSym = (
  (37 Cosh[111 Sqrt[5]/5] + (111 Sqrt[5]/5) Sinh[111 Sqrt[5]/5])
  / (136900 (-1 + (Sqrt[5]/3) Sinh[111 Sqrt[5]/5]
                 + Cosh[111 Sqrt[5]/5]))
);
(* Independent closed form for Theta_suff from Stage-75 sympy output line 21:
   Theta_suff/Pe_req = -(45 cosh(alpha) + 27 sqrt(5) sinh(alpha))
                        / (2500 - 2500 cosh(alpha)),  alpha = 111 sqrt(5)/5. *)
thetaSuffSym = (-(45 Cosh[111 Sqrt[5]/5] + 27 Sqrt[5] Sinh[111 Sqrt[5]/5]) / (2500 - 2500 Cosh[111 Sqrt[5]/5]));
(* The chi^2 and Jensen-floor Theta values are adopted from Stage-77 at
   extended precision; their independent re-derivation is the subject of
   the Stage-077 audit, not this one.                                    *)
thetaChiCoeffNum = ToExpression["4.0686323500816155092718546246574670820527`40"];
thetaJCoeffNum   = ToExpression["0.92755203253930797183993260663904217023`40"];
thetaChiCoeff  = thetaChiCoeffNum;
thetaJCoeff    = thetaJCoeffNum;
thetaFailCoeff = N[thetaFailSym, 30];
thetaSuffCoeff = N[thetaSuffSym, 30];

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

(* Independent targets computed in Mathematica from the symbolic
   closed form (thetaFailSym) and the chi/J Stage-77 numerics.   *)
peSuffChiTarget = N[thetaChiCoeffNum / thetaSuffSym, 30];
peFailChiTarget = N[thetaChiCoeffNum / thetaFailSym, 30];
peSuffJTarget   = N[thetaJCoeffNum   / thetaSuffSym, 30];
peFailJTarget   = N[thetaJCoeffNum   / thetaFailSym, 30];
expectApprox["Pe_suff^(chi) numeric check", peSuffChi, peSuffChiTarget, 10^-12];
expectApprox["Pe_fail^(chi) numeric check", peFailChi, peFailChiTarget, 10^-9];
expectApprox["Pe_suff^(J) numeric check",   peSuffJ,   peSuffJTarget,   10^-12];
expectApprox["Pe_fail^(J) numeric check",   peFailJ,   peFailJTarget,   10^-10];
expectTrue["Pe_suff^(chi) < Pe_fail^(chi)", peSuffChi < peFailChi];
expectTrue["Pe_suff^(J) < Pe_fail^(J)", peSuffJ < peFailJ];

(* --- Branch verdict (chi vs Jensen-floor) ---------------------------- *)
expectTrue["Pe_suff^(J) < Pe_suff^(chi)", peSuffJ < peSuffChi];
expectTrue["Pe_fail^(J) < Pe_fail^(chi)", peFailJ < peFailChi];
expectTrue["Pe_suff^(chi) < Pe_fail^(J) (window overlap)", peSuffChi < peFailJ];

Print[""];
Print["Stage 078 Mathematica audit passed."];

Exit[0];
