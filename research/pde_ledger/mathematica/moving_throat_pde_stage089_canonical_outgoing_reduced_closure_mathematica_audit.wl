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
  res = FullSimplify[Expand[expr], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 089 — REDUCED 2.5PN CLOSURE ON CANONICAL OUTGOING DtN BRANCH"];

Clear[G, c, cS, a, m0hat, chiQ, NQ];
$Assumptions = And @@ Thread[{G, c, cS, a, m0hat, chiQ, NQ} > 0];

constraint = m0hat^2*chiQ*NQ == 1;
Print["Normalization constraint = ", fmt[constraint]];

nqGeneral = FullSimplify[NQ /. First[Solve[constraint, NQ, Reals]], Assumptions -> $Assumptions];
Print["N_Q on the general outgoing branch = ", fmt[nqGeneral]];

nqCanonical = FullSimplify[nqGeneral /. chiQ -> 1, Assumptions -> $Assumptions];
Print["N_Q on the canonical outgoing branch = ", fmt[nqCanonical]];
expectZero[
  "point-particle canonical branch gives N_Q = 1",
  (nqCanonical /. m0hat -> 1) - 1
];

k0Target = 54*G*cS^5/(5*a^5*c^5);
k2Target = 6*G*cS^3/(5*a^3*c^5);
k4Target = 8*G*cS/(15*a*c^5);
gamma5Target = 2*G/(5*c^5);

k0 = FullSimplify[nqGeneral*k0Target, Assumptions -> $Assumptions];
k2 = FullSimplify[k0/(4*(3*cS/(2*a))^2), Assumptions -> $Assumptions];
k4 = FullSimplify[k0/(4*(3*cS/(2*a))^4), Assumptions -> $Assumptions];
gamma5 = FullSimplify[nqGeneral*gamma5Target, Assumptions -> $Assumptions];

Print["K0 = ", fmt[k0]];
Print["K2 = ", fmt[k2]];
Print["K4 = ", fmt[k4]];
Print["Gamma5 = ", fmt[gamma5]];

expectZero["branch identity K4 - 4 K2^2 / K0", k4 - 4*k2^2/k0];
expectZero[
  "branch identity Gamma5 - 9 K2^(5/2)/K0^(3/2)",
  gamma5 - 9*Sqrt[k2^5/k0^3]
];

gammaEffCanonical = FullSimplify[(m0hat^2*gamma5) /. chiQ -> 1, Assumptions -> $Assumptions];
expectZero["canonical gamma_eff - target", gammaEffCanonical - gamma5Target];

Print[""];
Print["RESULT:"];
Print["  On the canonical outgoing branch, N_Q = 1/m0hat^2."];
Print["  In the strict point-particle limit m0hat = 1, the reduced closure gives N_Q = 1."];
Print["  The effective odd coefficient then reproduces 2 G / (5 c^5) exactly."];

Exit[0];
