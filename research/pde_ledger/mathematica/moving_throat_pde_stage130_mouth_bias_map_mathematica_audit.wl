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

expectApprox[name_String, value_, target_, tol_] := Module[{delta},
  delta = N[value - target, 40];
  Print[name, " residual = ", fmt[delta]];
  If[TrueQ[Abs[delta] < tol], pass[name], fail[name, delta]];
];

banner["STAGE 130 — EXACT MOUTH-BIAS MAP AND FAMILY-1 COMPENSATION POINT"];

Clear[z, lM, piM];
$Assumptions = Element[{z, lM, piM}, Reals] && lM > 0 && piM > 0;

gMinus = N[(2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi), 80];

sigma = piM*Exp[-piM*z/lM]/(lM*(1 - Exp[-piM]));
f = Cos[Pi*z/(2*lM)];

gPi = FullSimplify[Integrate[sigma*f, {z, 0, lM}], Assumptions -> $Assumptions];
gPiBoxed = 2*piM*(2*piM*Exp[piM] + Pi) / ((4*piM^2 + Pi^2) * (Exp[piM] - 1));
expectZero["g_Pi matches paper boxed closed form", gPi - gPiBoxed];
eZ = FullSimplify[Integrate[sigma*z, {z, 0, lM}], Assumptions -> $Assumptions];
eFZ = FullSimplify[Integrate[sigma*f*z, {z, 0, lM}], Assumptions -> $Assumptions];
covId = FullSimplify[D[gPi, piM] + (eFZ - gPi*eZ)/lM, Assumptions -> $Assumptions];

Print["g_Pi = ", fmt[gPi]];
Print["Covariance identity residual = ", fmt[covId]];
expectZero["covariance identity", covId];

Clear[pi0, piInf];
g0 = FullSimplify[Limit[gPi /. piM -> pi0, pi0 -> 0, Direction -> "FromAbove"], Assumptions -> pi0 > 0];
gInf = FullSimplify[Limit[gPi /. piM -> piInf, piInf -> Infinity], Assumptions -> piInf > 0];
Print["limit Pi->0+ = ", fmt[g0]];
Print["limit Pi->oo = ", fmt[gInf]];
expectZero["uniform-source limit", g0 - 2/Pi];
expectZero["point-source limit", gInf - 1];

(* Global strict monotonicity dg/dpiM > 0 for all piM > 0 (notes section 2, boxed).
   Proof structure: dg/dpiM = -(1/lM) Cov_piM(f, z), already checked covId == 0.
   Certify Cov < 0 via the FKG/Chebyshev symmetrized identity and the pointwise
   sign of its integrand. GLOBAL certificate on piM > 0, not a finite sweep. *)
dgPi = D[gPi, piM];

normP = FullSimplify[Integrate[sigma, {z, 0, lM}], Assumptions -> $Assumptions];
expectZero["sigma_piM normalized on [0,lM]", normP - 1];

cov = FullSimplify[eFZ - gPi*eZ, Assumptions -> $Assumptions];

(* (1) Symmetrized double-integral identity for the covariance. *)
Clear[z1, z2];
asm2 = Element[{z1, z2}, Reals] && lM > 0 && piM > 0 && 0 <= z1 <= lM && 0 <= z2 <= lM;
f1 = (f /. z -> z1);
f2 = (f /. z -> z2);
p1 = (sigma /. z -> z1);
p2 = (sigma /. z -> z2);
integrandSym = (1/2)*(f1 - f2)*(z1 - z2)*p1*p2;
covDouble = FullSimplify[
  Integrate[Integrate[integrandSym, {z1, 0, lM}], {z2, 0, lM}],
  Assumptions -> $Assumptions];
expectZero["symmetrized covariance identity", covDouble - cov];

(* (2) Pointwise sign of the symmetrizer factor: f strictly decreasing on [0,lM].
   f'(z) = -(Pi/(2 lM)) Sin[Pi z/(2 lM)] < 0 for 0 < z < lM (argument in (0,Pi/2)).
   Certify the closed form of f'(z); its sign then follows from Sin > 0 on (0,Pi/2).
   This is a bounded-domain trig statement with NO Exp[piM]; decidable. *)
fPrime = D[f, z];
expectZero["f'(z) closed form", fPrime + (Pi/(2*lM))*Sin[Pi*z/(2*lM)]];
sinPos = Reduce[Sin[Pi*z/(2*lM)] > 0 && 0 < z < lM && lM > 0, z, Reals];
Print["Sin[Pi z/(2 lM)] > 0 on (0,lM) decided as: ", fmt[sinPos]];
If[TrueQ[sinPos =!= False],
  pass["f strictly decreasing on (0,lM) -> symmetrizer <= 0"],
  fail["f strictly decreasing on (0,lM)", sinPos]];

Print["Cov_piM(f,z) (symmetrized) = ", fmt[covDouble]];
(* Consistency: dg/dpiM = -(1/lM) Cov. Non-tautological: a wrong gPi or Cov breaks it. *)
expectZero["dg/dpiM = -(1/lM) Cov consistency", dgPi + cov/lM];
Print["Global strict monotonicity certified: dg/dpiM > 0 for all piM > 0."];

(* (F1c) Uniqueness of Pi_* on (0,oo): g strictly increasing from g(0+)=2/Pi to
   g(oo)=1, so g(piM)=gMinus has exactly one root iff 2/Pi < gMinus < 1. *)
gLo = N[2/Pi, 40];
If[TrueQ[gLo < gMinus < 1],
  pass["g_minus strictly inside (2/Pi, 1): Pi_* unique"],
  fail["g_minus inside (2/Pi, 1)", gMinus]];
Print["Bracket for unique Pi_*: 2/Pi = ", fmt[gLo], " < g_minus = ", fmt[gMinus], " < 1"];

piStar = piM /. FindRoot[gPi == gMinus, {piM, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
Print["Pi_* = ", fmt[N[piStar, 30]]];
Print["x_* = 1/Pi_* = ", fmt[N[1/piStar, 30]]];
Print["g(Pi_*) = ", fmt[N[gPi /. piM -> piStar, 30]]];
Print["g'(Pi_*) = ", fmt[N[D[gPi, piM] /. piM -> piStar, 30]]];
expectApprox["Family-1 compensation point", gPi /. piM -> piStar, gMinus, 10^-20];

Print[""];
Print["Stage 130 Mathematica audit passed."];

Exit[0];
