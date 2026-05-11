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

banner["STAGE 078 — SECOND-ORDER GEOMETRY CONTAMINATION"];

Clear[chi, m0, g0, g2, g4, w, oQ, kPole];
$Assumptions =
  Element[{chi, m0, g0, g2, g4, w, oQ, kPole}, Reals] &&
  g0 != 0 && oQ != 0 && kPole != 0;

dG = g0 + g2*w^2 + g4*w^4;
corr = Expand[Normal[Series[-chi^2*m0^2/dG, {w, 0, 5}]]];
k0Corr = FullSimplify[corr /. w -> 0, Assumptions -> $Assumptions];
k2Corr = FullSimplify[Coefficient[corr, w, 2], Assumptions -> $Assumptions];
k4Corr = FullSimplify[Coefficient[corr, w, 4], Assumptions -> $Assumptions];

Print["Schur-complement correction = ", fmt[corr]];
Print["K0corr = ", fmt[k0Corr]];
Print["K2corr = ", fmt[k2Corr]];
Print["K4corr = ", fmt[k4Corr]];

expectZero["K0corr / chi^2 static factor", k0Corr/chi^2 + m0^2/g0];
expectZero["K2corr / chi^2 dynamic factor", k2Corr/chi^2 - g2*m0^2/g0^2];
expectZero["K4corr / chi^2 dynamic factor", k4Corr/chi^2 - m0^2*(g0*g4 - g2^2)/g0^3];

eps2 = FullSimplify[oQ^2*k2Corr/kPole, Assumptions -> $Assumptions];
eps4 = FullSimplify[oQ^4*k4Corr/kPole, Assumptions -> $Assumptions];
Print["eps2 = ", fmt[eps2]];
Print["eps4 = ", fmt[eps4]];

cPole = FullSimplify[(1 + eps4)/(4*(1 + eps2)^2), Assumptions -> $Assumptions];
cPoleSeries = Expand[Normal[Series[cPole, {chi, 0, 2}]]];
delta = FullSimplify[cPoleSeries - 1/4, Assumptions -> $Assumptions];

Print["c_pole = ", fmt[cPoleSeries]];
Print["delta c_pole = ", fmt[delta]];
expectZero["d c_pole / dchi |0", (D[cPole, chi] /. chi -> 0)];

Print["d^2 c_pole / dchi^2 |0 = ", fmt[FullSimplify[D[cPole, {chi, 2}] /. chi -> 0, Assumptions -> $Assumptions]]];

Print[""];
Print["Stage 095 Mathematica audit passed."];

Exit[0];
