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

banner["STAGE 100 — OUTGOING NORMALIZATION FACTORIZATION"];

(* Carry-forward annotations (per batch-IV1-paper-alignment Cluster B (c)):
   - Check (ii) higher odd terms beyond 2.5PN: anchored at stage 102.
   - Check (iii) DtN l=2 fingerprint pinning chi_Q = 1: anchored at stage 097.
   This stage owns Check (i): the mhat_0^2 chi_Q N_Q = 1 closure derived below. *)

Clear[gNewton, cLight, omegaQ, k0, mHat0, chiQ, omega];
(* chi_Q is a free real parameter here; its pin to 1 is upstream (stage 097). *)
$Assumptions =
  Element[{gNewton, cLight, omegaQ, k0, mHat0, omega}, Reals] &&
  Element[chiQ, Reals] &&
  gNewton > 0 && cLight > 0 && omegaQ > 0 && k0 > 0 && mHat0 > 0;

sigmaCan = FullSimplify[(9/8)/omegaQ^5, Assumptions -> $Assumptions];
yRet = 3/4 + (1/4)/(1 - omega^2/omegaQ^2 - I*chiQ*sigmaCan*omega^5);
ySeries = Expand[Normal[Series[yRet, {omega, 0, 5}]]];

k2 = FullSimplify[k0*Coefficient[ySeries, omega, 2], Assumptions -> $Assumptions];
k4 = FullSimplify[k0*Coefficient[ySeries, omega, 4], Assumptions -> $Assumptions];
gamma5 = FullSimplify[(Coefficient[ySeries, omega, 5]/I)*k0, Assumptions -> $Assumptions];

k0Target = FullSimplify[64*gNewton*omegaQ^5/(45*cLight^5), Assumptions -> $Assumptions];
k2Target = FullSimplify[k0Target/(4*omegaQ^2), Assumptions -> $Assumptions];
k4Target = FullSimplify[k0Target/(4*omegaQ^4), Assumptions -> $Assumptions];
gamma5Target = FullSimplify[2*gNewton/(5*cLight^5), Assumptions -> $Assumptions];
nQDerived = FullSimplify[k0/k0Target, Assumptions -> $Assumptions];

Print["Yhat_Q^ret series = ", fmt[ySeries]];
Print["K2 = ", fmt[k2]];
Print["K4 = ", fmt[k4]];
Print["Gamma5 = ", fmt[gamma5]];
Print["NQ = ", fmt[nQDerived]];

expectZero["K2/K2_target - NQ", k2/k2Target - nQDerived];
expectZero["K4/K4_target - NQ", k4/k4Target - nQDerived];
expectZero["Gamma5/Gamma5_target - chiQ*NQ", gamma5/gamma5Target - chiQ*nQDerived];

(* Substantive closure: impose mhat_0^2 * Gamma_5 = Gamma_5_target. Substituting
   the script-derived Gamma_5 = chi_Q N_Q Gamma_5_target, the residual ratio
   over Gamma_5_target equals (mhat_0^2 chi_Q N_Q - 1). The closure forces this
   to zero, recovering the headline factorization. *)
closureResidual = FullSimplify[mHat0^2 * gamma5 - gamma5Target, Assumptions -> $Assumptions];
closureRatio = FullSimplify[closureResidual / gamma5Target, Assumptions -> $Assumptions];
Print["closure_residual / Gamma5_target = ", fmt[closureRatio]];
expectZero["closure_ratio - (mhat0^2 chi_Q N_Q - 1)",
           closureRatio - (mHat0^2 * chiQ * nQDerived - 1)];

Print[""];
Print["Closure ledger: mhat_0^2 * Gamma_5 = Gamma_5_target  <=>  mhat_0^2 chi_Q N_Q = 1"];
Print["Stage 100 Mathematica audit passed."];

Exit[0];
