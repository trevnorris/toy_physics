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

banner["FULL-PROFILE FAMILY-1 MOUTH RESIDUAL"];

Clear[x, p, sigmaM];
$Assumptions = Element[{x, p, sigmaM}, Reals] && p > 0 && sigmaM > 0 && 0 <= x <= 1 && p != Pi/2;

k = Pi/2;
ts = (1 - Exp[-p*x])/(p*(1 - Exp[-p])) - x*Exp[-p]/(1 - Exp[-p]);
cq = p/((1 - Exp[-p])*(k^2 - p^2));
aq = cq*(k*Sinh[k] + p*Exp[-p])/(k*Cosh[k]);
tq = aq*Sinh[k*x] - cq*Cosh[k*x] + cq*Exp[-p*x];

(* Hand-derived closed form: T_q'(0) = aq*k - cq*p
   (differentiate tq = aq*Sinh[k*x] - cq*Cosh[k*x] + cq*Exp[-p*x] at x=0). *)
(* Build the slope compactly from FREE coefficient symbols so the PRINTED form is provably
   the real slope; then substitute the concrete aq, cq definitions for the load-bearing checks. *)
sQsymbolic = aqS*k - cqS*p;
sQ = sQsymbolic /. {aqS -> aq, cqS -> cq};

Print["S_q(Pi) = ", fmt[sQsymbolic]];
Print["S_q(Pi) [expanded] = ", fmt[sQ]];

expectZero["T_s(0)", ts /. x -> 0];
expectZero["T_q(0)", tq /. x -> 0];
expectZero["T_s'(0)-1", (D[ts, x] /. x -> 0) - 1];
expectZero["T_q'(0)-S_q", (D[tq, x] /. x -> 0) - sQ];

r = sigmaM*(4*ts - tq - (4 - sQ)*x);
expectZero["R(0)", r /. x -> 0];
expectZero["R'(0)", D[r, x] /. x -> 0];

r2 = FullSimplify[D[r, {x, 2}] /. x -> 0, Assumptions -> $Assumptions];
targetR2 = FullSimplify[-3*sigmaM*p/(1 - Exp[-p]), Assumptions -> $Assumptions];
Print["R''(0) = ", fmt[r2]];
expectZero["R''(0) - target", r2 - targetR2];

Print[""];
Print["Theorem:"];
Print["  R(0)=0, R'(0)=0, and R''(0) = -3 Sigma Pi/(1-exp(-Pi)) < 0."];
Print["  So the full compensated mouth potential is tangent-matched but sublinear at the mouth."];

Print[""];
Print["Stage 150 Mathematica audit passed."];

Exit[0];
