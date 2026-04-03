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

banner["STAGE 082 — REDUCED FINISH LINE"];

Clear[G, c, cS, a, omega, omegaQ, nQ];
$Assumptions =
  Element[{G, c, cS, a, omega, omegaQ, nQ}, Reals] &&
  And @@ Thread[{G, c, cS, a, omegaQ, nQ} > 0];

yhatCons = FullSimplify[3/4 + 1/(4*(1 - omega^2/omegaQ^2)), Assumptions -> $Assumptions];
Print["Yhat_Q^cons(omega) = ", fmt[yhatCons]];

k0Target = FullSimplify[64*G*omegaQ^5/(45*c^5), Assumptions -> $Assumptions];
k0TargetGeom = FullSimplify[k0Target /. omegaQ -> 3*cS/(2*a), Assumptions -> $Assumptions];
expectZero["K0_target geometric form", k0TargetGeom - 54*G*cS^5/(5*a^5*c^5)];

kbar = FullSimplify[nQ*k0Target*yhatCons, Assumptions -> $Assumptions];
series = FullSimplify[Normal[Series[kbar, {omega, 0, 4}]], Assumptions -> $Assumptions];

k0 = FullSimplify[series /. omega -> 0, Assumptions -> $Assumptions];
k2 = FullSimplify[Coefficient[series, omega, 2], Assumptions -> $Assumptions];
k4 = FullSimplify[Coefficient[series, omega, 4], Assumptions -> $Assumptions];
gamma5 = FullSimplify[9*Sqrt[k2^5/k0^3], Assumptions -> $Assumptions];

Print["Kbar_Q^cons series = ", fmt[series]];
Print["K0 = ", fmt[k0]];
Print["K2 = ", fmt[k2]];
Print["K4 = ", fmt[k4]];
Print["Gamma5 = ", fmt[gamma5]];

k2Target = FullSimplify[k0Target/(4*omegaQ^2), Assumptions -> $Assumptions];
k4Target = FullSimplify[k0Target/(4*omegaQ^4), Assumptions -> $Assumptions];
gamma5Target = FullSimplify[9*Sqrt[k2Target^5/k0Target^3], Assumptions -> $Assumptions];

r0 = FullSimplify[k0/k0Target - 1, Assumptions -> $Assumptions];
r2 = FullSimplify[k2/k2Target - 1, Assumptions -> $Assumptions];
r4 = FullSimplify[k4/k4Target - 1, Assumptions -> $Assumptions];
r5 = FullSimplify[gamma5/gamma5Target - 1, Assumptions -> $Assumptions];

expectZero["R0 - (N_Q - 1)", r0 - (nQ - 1)];
expectZero["R2 - (N_Q - 1)", r2 - (nQ - 1)];
expectZero["R4 - (N_Q - 1)", r4 - (nQ - 1)];
expectZero["R5 - (N_Q - 1)", r5 - (nQ - 1)];

Print["R0 = ", fmt[Factor[r0]]];
Print["R2 = ", fmt[Factor[r2]]];
Print["R4 = ", fmt[Factor[r4]]];
Print["R5 = ", fmt[Factor[r5]]];
Print[""];
Print["STAGE 082 AUDIT PASSED"];

Exit[0];
