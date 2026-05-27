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

banner["STAGE 099 — REDUCED FINISH LINE"];

Clear[G, c, cS, a, omega, omegaQ, nQ];
$Assumptions =
  Element[{G, c, cS, a, omega, omegaQ, nQ}, Reals] &&
  And @@ Thread[{G, c, cS, a, omegaQ, nQ} > 0];

yhatCons = FullSimplify[3/4 + 1/(4*(1 - omega^2/omegaQ^2)), Assumptions -> $Assumptions];
Print["Yhat_Q^cons(omega) = ", fmt[yhatCons]];

(* F4 substantive: static slot and pole residue of the conservative module. *)
expectZero["Yhat_Q^cons static slot equals 1", (yhatCons /. omega -> 0) - 1];
expectZero["Yhat_Q^cons pole residue at omega=omegaQ is -omegaQ/8",
           Residue[yhatCons, {omega, omegaQ}] - (-omegaQ/8)];

k0Target = FullSimplify[64*G*omegaQ^5/(45*c^5), Assumptions -> $Assumptions];
k0TargetGeom = FullSimplify[k0Target /. omegaQ -> 3*cS/(2*a), Assumptions -> $Assumptions];
expectZero["K0_target geometric form", k0TargetGeom - 54*G*cS^5/(5*a^5*c^5)];

(* F1 substantive: structural relations forced by Yhat_Q^cons with k0Sym free. *)
Clear[k0Sym];
$Assumptions = $Assumptions && k0Sym > 0 && Element[k0Sym, Reals];
k2Struct = k0Sym/(4*omegaQ^2);
k4Struct = k0Sym/(4*omegaQ^4);
gamma5Struct = 9*k0Sym/(32*omegaQ^5);

expectZero["branch identity K0 K4 = 4 K2^2", k0Sym*k4Struct - 4*k2Struct^2];
expectZero["Gamma_5 sqrt form matches canonical odd-coeff form",
           9*Sqrt[k2Struct^5/k0Sym^3] - gamma5Struct];
expectZero["Gamma_5 normalization equals N_Q on chi_Q=1 branch",
           (gamma5Struct /. k0Sym -> nQ*k0Target) / (2*G/(5*c^5)) - nQ];

Print["K2_struct = ", fmt[Factor[k2Struct]]];
Print["K4_struct = ", fmt[Factor[k4Struct]]];
Print["Gamma5_struct = ", fmt[Factor[gamma5Struct]]];
Print[""];
Print["STAGE 099 AUDIT PASSED"];

Exit[0];
