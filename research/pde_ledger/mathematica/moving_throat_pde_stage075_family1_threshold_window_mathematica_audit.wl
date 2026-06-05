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

banner["STAGE 075 — EXPLICIT FAMILY-1 THRESHOLD WINDOW"];

Clear[peReq, thetaW];
$Assumptions = Element[{peReq, thetaW}, Reals] && peReq > 0 && thetaW > 0;

alphaR = 10;
(* F4 (v2 paper-alignment Q2 direction (a) lock): paper Inputs line states
   Upsilon_w = alpha_r^2 Theta_w with alpha_r^2 = 100. Lock the value so any
   future drift between the paper Inputs line and the script surfaces here. *)
expectZero["alpha_r^2 - 100 (paper Inputs line lock)", alphaR^2 - 100];
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

(* Note: the algebraic identity below is structurally tautological (it follows
   from the definition of Delta_0 / Delta_inf by canceling a common factor).
   The genuine independent check is the asymptotic-limit block further below,
   which exercises a non-trivial property of the closed forms via Mathematica's
   Limit operator (computed independently from SymPy's sp.limit). *)
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

(* F1 (v2): asymptotic-limit checks for Delta_0 and Delta_inf. These are
   non-trivial consequences of the closed forms (a wrong factor would change
   the limit), and Mathematica's Limit operator computes them by an algorithm
   independent of SymPy's sp.limit. *)
Module[{aSym, eSym, delta0Sym, deltaInfSym, largeAlphaLimit, smallAlphaLimit},
  ClearAll[aSym, eSym];
  delta0Sym = eSym*(Cosh[aSym] - 1)/(aSym^2*(aSym*Sinh[aSym] + eSym*Cosh[aSym]));
  deltaInfSym = (Cosh[aSym] + (eSym/aSym)*Sinh[aSym] - 1)/(aSym*Sinh[aSym] + eSym*Cosh[aSym]);
  largeAlphaLimit = Limit[aSym*deltaInfSym, aSym -> Infinity, Assumptions -> eSym > 0];
  Print["alpha * Delta_inf large-alpha limit = ", fmt[largeAlphaLimit]];
  If[TrueQ[largeAlphaLimit === 1],
    pass["alpha * Delta_inf -> 1 (large alpha)"],
    fail["alpha * Delta_inf -> 1 (large alpha)", largeAlphaLimit]];
  smallAlphaLimit = Limit[delta0Sym, aSym -> 0, Assumptions -> eSym > 0];
  Print["Delta_0 small-alpha limit = ", fmt[smallAlphaLimit]];
  If[TrueQ[smallAlphaLimit === 1/2],
    pass["Delta_0 -> 1/2 (small alpha)"],
    fail["Delta_0 -> 1/2 (small alpha)", smallAlphaLimit]];
];

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
