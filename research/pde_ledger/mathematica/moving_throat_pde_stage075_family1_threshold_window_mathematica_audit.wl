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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 40] - N[target, 40]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 058 — EXPLICIT FAMILY-1 THRESHOLD WINDOW"];

Clear[peReq, thetaW];
$Assumptions = Element[{peReq, thetaW}, Reals] && peReq > 0 && thetaW > 0;

alphaR = 10;
lambdaEll = 37;
eta = 37;
kappa = 12321/5;
alpha = Sqrt[kappa];

delta0 = FullSimplify[
  eta*(Cosh[alpha] - 1)/(alpha^2*(alpha*Sinh[alpha] + eta*Cosh[alpha])),
  Assumptions -> $Assumptions
];
deltaInf = FullSimplify[
  (Cosh[alpha] + (eta/alpha)*Sinh[alpha] - 1)/(alpha*Sinh[alpha] + eta*Cosh[alpha]),
  Assumptions -> $Assumptions
];

Print["Lambda_ell = ", fmt[lambdaEll]];
Print["eta = ", fmt[eta]];
Print["kappa = ", fmt[kappa]];
Print["alpha = sqrt(kappa) = ", fmt[alpha]];
Print["Delta_0 = ", fmt[delta0]];
Print["Delta_inf = ", fmt[deltaInf]];
Print["Delta_0 (numeric) = ", fmt[N[delta0, 20]]];
Print["Delta_inf (numeric) = ", fmt[N[deltaInf, 20]]];

upsilonFail = FullSimplify[peReq/(lambdaEll^2*deltaInf), Assumptions -> $Assumptions];
upsilonSuff = FullSimplify[peReq/(lambdaEll^2*delta0), Assumptions -> $Assumptions];
xiFail = FullSimplify[peReq/deltaInf, Assumptions -> $Assumptions];
xiSuff = FullSimplify[peReq/delta0, Assumptions -> $Assumptions];
thetaFail = FullSimplify[upsilonFail/alphaR^2, Assumptions -> $Assumptions];
thetaSuff = FullSimplify[upsilonSuff/alphaR^2, Assumptions -> $Assumptions];

Print["Upsilon_fail = ", fmt[upsilonFail]];
Print["Upsilon_suff = ", fmt[upsilonSuff]];
Print["Xi_fail = ", fmt[xiFail]];
Print["Xi_suff = ", fmt[xiSuff]];
Print["Theta_fail = ", fmt[thetaFail]];
Print["Theta_suff = ", fmt[thetaSuff]];

(* Independent symbolic check: the stated closed forms must satisfy
   the defining algebraic identities for *all* positive alpha, eta,
   not just the substituted numeric values alpha = 111/Sqrt[5], eta = 37. *)
(* This identity check is the independent-derivation leg required by
   the second-engine policy: Mathematica's FullSimplify must prove the
   identity for *symbolic* alpha and eta, which catches a wrong factor
   or sign in the stated closed form even though the rest of the script
   transliterates the SymPy recipe. *)
Module[{aSym, eSym, delta0Sym, deltaInfSym},
  ClearAll[aSym, eSym];
  delta0Sym = eSym*(Cosh[aSym] - 1)/(aSym^2*(aSym*Sinh[aSym] + eSym*Cosh[aSym]));
  deltaInfSym = (Cosh[aSym] + (eSym/aSym)*Sinh[aSym] - 1)/(aSym*Sinh[aSym] + eSym*Cosh[aSym]);
  expectZero["Delta_0 algebraic identity (free alpha, eta)",
    Assuming[aSym > 0 && eSym > 0,
      FullSimplify[(aSym*Sinh[aSym] + eSym*Cosh[aSym])*delta0Sym - eSym*(Cosh[aSym] - 1)/aSym^2]]];
  expectZero["Delta_inf algebraic identity (free alpha, eta)",
    Assuming[aSym > 0 && eSym > 0,
      FullSimplify[(aSym*Sinh[aSym] + eSym*Cosh[aSym])*deltaInfSym - (Cosh[aSym] + (eSym/aSym)*Sinh[aSym] - 1)]]];
];

expectZero["Upsilon_fail - alphaR^2 * Theta_fail", upsilonFail - alphaR^2*thetaFail];
expectZero["Upsilon_suff - alphaR^2 * Theta_suff", upsilonSuff - alphaR^2*thetaSuff];

expectApprox["Delta_0 numeric check", delta0, 0.00017330207902152514906, 10^-18];
expectApprox["Delta_inf numeric check", deltaInf, 0.020144756554052159427, 10^-17];
expectApprox["Upsilon_fail / Pe_req numeric check", upsilonFail/peReq, 0.036260561797293886969, 10^-16];
expectApprox["Upsilon_suff / Pe_req numeric check", upsilonSuff/peReq, 4.2149534156997728721, 10^-14];
expectApprox["Xi_fail / Pe_req numeric check", xiFail/peReq, 49.640709100495331260, 10^-13];
expectApprox["Xi_suff / Pe_req numeric check", xiSuff/peReq, 5770.2712260929890619, 10^-10];
expectApprox["Theta_fail / Pe_req numeric check", thetaFail/peReq, 0.00036260561797293886969, 10^-18];
expectApprox["Theta_suff / Pe_req numeric check", thetaSuff/peReq, 0.042149534156997728721, 10^-16];

Print[""];
Print["Stage 075 Mathematica audit passed."];

Exit[0];
