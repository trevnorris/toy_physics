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

banner["STAGE 097 — SINGLE NORMALIZATION DEFECT"];

Clear[gConst, cLight, cSound, aRad, k0, omegaQ, nQ, w];
$Assumptions =
  Element[{gConst, cLight, cSound, aRad, k0, omegaQ, nQ, w}, Reals] &&
  gConst > 0 && cLight > 0 && cSound > 0 && aRad > 0 && k0 > 0 && omegaQ > 0 && nQ > 0;

(* Independent derivation route (different from SymPy's direct definitions):
   build K_n by extracting series coefficients of the conservative quadrupole
   module Yhat_Q^cons(w) = 3/4 + (1/4)/(1 - w^2/omegaQ^2) carried by K_0. *)
yhatCons[wvar_] := 3/4 + (1/4)/(1 - wvar^2/omegaQ^2);
kbarCons[wvar_] := k0 * yhatCons[wvar];

k2 = FullSimplify[SeriesCoefficient[kbarCons[w], {w, 0, 2}], Assumptions -> $Assumptions];
k4 = FullSimplify[SeriesCoefficient[kbarCons[w], {w, 0, 4}], Assumptions -> $Assumptions];
expectZero["k2 series == k0/(4 omegaQ^2)", k2 - k0/(4*omegaQ^2)];
expectZero["k4 series == k0/(4 omegaQ^4)", k4 - k0/(4*omegaQ^4)];

(* gamma5 from the series-derived k2 (the load-bearing odd coefficient). *)
gamma5 = FullSimplify[9*k2^(5/2)/k0^(3/2), Assumptions -> $Assumptions];
gamma5Expected = FullSimplify[9*k0/(32*omegaQ^5), Assumptions -> $Assumptions];

Print["K2 = ", fmt[k2]];
Print["K4 = ", fmt[k4]];
Print["Gamma5 = ", fmt[gamma5Expected]];
expectZero["Gamma5 - 9 K0/(32 Omega^5)", gamma5 - gamma5Expected];

k0Target = FullSimplify[64*gConst*omegaQ^5/(45*cLight^5), Assumptions -> $Assumptions];
omegaGeom = FullSimplify[3*cSound/(2*aRad), Assumptions -> $Assumptions];
k0TargetGeom = FullSimplify[k0Target /. omegaQ -> omegaGeom, Assumptions -> $Assumptions];

Print["K0_target = ", fmt[k0Target]];
Print["K0_target (Omega=3 c_s/(2a)) = ", fmt[k0TargetGeom]];
expectZero["geometric target reduction", k0TargetGeom - 54*gConst*cSound^5/(5*aRad^5*cLight^5)];

(* k2Target, k4Target via the same series route, gamma5Target via the
   already-established 9 k0/(32 omegaQ^5) form (not 9 k2^(5/2)/k0^(3/2)). *)
k2Target = FullSimplify[SeriesCoefficient[k0Target*yhatCons[w], {w, 0, 2}], Assumptions -> $Assumptions];
k4Target = FullSimplify[SeriesCoefficient[k0Target*yhatCons[w], {w, 0, 4}], Assumptions -> $Assumptions];
gamma5Target = FullSimplify[9*k0Target/(32*omegaQ^5), Assumptions -> $Assumptions];
expectZero["Gamma5_target - 2G/(5c^5)", gamma5Target - 2*gConst/(5*cLight^5)];

(* R_i reductions via actual-branch coefficients built before simplification. *)
k0Actual = nQ * k0Target;
k2Actual = FullSimplify[SeriesCoefficient[k0Actual*yhatCons[w], {w, 0, 2}], Assumptions -> $Assumptions];
k4Actual = FullSimplify[SeriesCoefficient[k0Actual*yhatCons[w], {w, 0, 4}], Assumptions -> $Assumptions];
gamma5Actual = FullSimplify[9*k0Actual/(32*omegaQ^5), Assumptions -> $Assumptions];

r0 = FullSimplify[k0Actual/k0Target - 1, Assumptions -> $Assumptions];
r2 = FullSimplify[k2Actual/k2Target - 1, Assumptions -> $Assumptions];
r4 = FullSimplify[k4Actual/k4Target - 1, Assumptions -> $Assumptions];
r5 = FullSimplify[gamma5Actual/gamma5Target - 1, Assumptions -> $Assumptions];

Print["R0 = ", fmt[Factor[r0]]];
Print["R2 = ", fmt[Factor[r2]]];
Print["R4 = ", fmt[Factor[r4]]];
Print["R5 = ", fmt[Factor[r5]]];

expectZero["R0 - (N_Q - 1)", r0 - (nQ - 1)];
expectZero["R2 - (N_Q - 1)", r2 - (nQ - 1)];
expectZero["R4 - (N_Q - 1)", r4 - (nQ - 1)];
expectZero["R5 - (N_Q - 1)", r5 - (nQ - 1)];

Print[""];
Print["Stage 097 Mathematica audit passed."];

Exit[0];
