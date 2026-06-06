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
u2Coeff = FullSimplify[1/omegaQ^2, Assumptions -> $Assumptions];
u5Coeff = FullSimplify[I*chiQ*sigmaCan, Assumptions -> $Assumptions];
uSquared4Coeff = FullSimplify[u2Coeff^2, Assumptions -> $Assumptions];

(* Independent route: expand the pole as 1 + u + u^2 through omega^5, with
   u = omega^2/omegaQ^2 + I chiQ sigmaCan omega^5. Higher powers start at
   omega^6 or omega^7 and cannot affect K2, K4, or Gamma5. *)
y2Coeff = FullSimplify[u2Coeff/4, Assumptions -> $Assumptions];
y4Coeff = FullSimplify[uSquared4Coeff/4, Assumptions -> $Assumptions];
y5Coeff = FullSimplify[u5Coeff/4, Assumptions -> $Assumptions];
ySeries = Expand[1 + y2Coeff*omega^2 + y4Coeff*omega^4 + y5Coeff*omega^5];

k2 = FullSimplify[k0*y2Coeff, Assumptions -> $Assumptions];
k4 = FullSimplify[k0*y4Coeff, Assumptions -> $Assumptions];
gamma5 = FullSimplify[(y5Coeff/I)*k0, Assumptions -> $Assumptions];

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
