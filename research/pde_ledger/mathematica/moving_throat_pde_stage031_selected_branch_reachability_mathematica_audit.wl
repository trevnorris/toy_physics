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

banner["PART I — INDEPENDENT EIGEN-DECOMPOSITION AND HF IDENTITY"];

Clear[a, dK, alpha, kappa0, kappa1, beta0];
$Assumptions = Element[{a, dK, alpha, kappa0, kappa1, beta0}, Reals] &&
  a > 0 && dK > 0 && alpha >= 0 && kappa0 > 0 && kappa1 > 0 && beta0 > 0;

(* Loaded 2x2 wall matrix, built directly from the physical setup:
   diagonal stiffnesses (a, a+dK) minus rank-one loading alpha * v v^T
   with loading vector v = {kappa0, kappa1}. This is the rank-one loaded
   wall problem of Stage 028 expressed without the (sigma, deltaKappa)
   choreography of the SymPy mirror. *)
wallMatrix = {{a, 0}, {0, a + dK}} -
  alpha*Outer[Times, {kappa0, kappa1}, {kappa0, kappa1}];

(* Independent eigenvalue extraction via Eigenvalues, then identify
   the lower branch as the one that reduces to a at alpha=0. *)
eigvals = Sort[Eigenvalues[wallMatrix]];
lamMinusIndep = eigvals[[1]];
lamPlusIndep  = eigvals[[2]];
expectZero["lam_-(0) initial value", (lamMinusIndep /. alpha -> 0) - a];
expectZero["lam_+(0) initial value", (lamPlusIndep  /. alpha -> 0) - (a + dK)];

(* Independent eigenvector extraction, then loading-vector overlap.
   The selected overlap s_- = (v . e_-)^2 / (e_- . e_-) is the
   physical definition used in the paper notes. *)
eigvecs = Eigenvectors[wallMatrix];
(* Pair eigenvectors with eigenvalues so we always pick the lower one *)
pairs = Sort[Transpose[{Eigenvalues[wallMatrix], eigvecs}], #1[[1]] < #2[[1]] &];
eMinus = pairs[[1, 2]];
loadingVec = {kappa0, kappa1};
sMinusOverlap = FullSimplify[
  (loadingVec . eMinus)^2 / (eMinus . eMinus),
  Assumptions -> $Assumptions
];

(* HF identity: d lam_- / d alpha = - s_-. Independent verification,
   not a definition. *)
hfResidual = FullSimplify[
  D[lamMinusIndep, alpha] + sMinusOverlap,
  Assumptions -> $Assumptions
];
expectZero["Hellmann-Feynman d lam_-/d alpha + s_- = 0", hfResidual];

banner["PART II — EXACT OVERLAP DERIVATIVE FROM OVERLAP DEFINITION"];

dsOverlap = FullSimplify[D[sMinusOverlap, alpha], Assumptions -> $Assumptions];
rSquared = FullSimplify[(dK + alpha*(kappa0^2 - kappa1^2))^2 + 4*alpha^2*kappa0^2*kappa1^2,
                        Assumptions -> $Assumptions];
dsExpected = FullSimplify[2*dK^2*kappa0^2*kappa1^2 / rSquared^(3/2),
                          Assumptions -> $Assumptions];
expectZero["ds_-/dalpha closed form", FullSimplify[dsOverlap - dsExpected, Assumptions -> $Assumptions]];
Print["ds_-/dalpha = ", fmt[dsExpected]];

banner["PART III — EXACT PREFACTOR DERIVATIVE"];

p0SelIndep = FullSimplify[beta0*sMinusOverlap/lamMinusIndep, Assumptions -> $Assumptions];
dPDirect = FullSimplify[D[p0SelIndep, alpha], Assumptions -> $Assumptions];
dPClosed = FullSimplify[
  beta0*(dsExpected*lamMinusIndep + sMinusOverlap^2)/lamMinusIndep^2,
  Assumptions -> $Assumptions
];
expectZero["dP0_-/dalpha closed form", FullSimplify[dPDirect - dPClosed, Assumptions -> $Assumptions]];
Print["dP0_-/dalpha = ", fmt[dPClosed]];

banner["PART IV — INITIAL VALUES AT alpha=0"];

expectZero["s_-(0) = kappa0^2", (sMinusOverlap /. alpha -> 0) - kappa0^2];
expectZero["P0_-(0) = beta0 kappa0^2 / a", (p0SelIndep /. alpha -> 0) - beta0*kappa0^2/a];

banner["PART V — SOFTENING THRESHOLD FROM Det[M]=0"];

(* alpha_crit is the root of Det[M] = 0, derived from the matrix directly. *)
detM = FullSimplify[Det[wallMatrix], Assumptions -> $Assumptions];
critSolutions = Solve[detM == 0, alpha];
(* The stable-side threshold is the smaller positive root. Both solutions are
   non-negative; pick the one that equals the closed form. *)
alphaCritClosed = FullSimplify[
  a*(a + dK) / ((a + dK)*kappa0^2 + a*kappa1^2),
  Assumptions -> $Assumptions
];
alphaCritFromDet = FullSimplify[alpha /. First[critSolutions], Assumptions -> $Assumptions];
expectZero["alpha_crit from Det[M]=0", FullSimplify[alphaCritFromDet - alphaCritClosed, Assumptions -> $Assumptions]];
Print["alpha_crit = ", fmt[alphaCritClosed]];

expectZero["lam_-(alpha_crit) = 0",
  FullSimplify[lamMinusIndep /. alpha -> alphaCritClosed, Assumptions -> $Assumptions]
];

banner["PART VI — STABLE-SIDE POLE STRUCTURE"];

(* lam_- lam_+ = Det[M] (basic linear algebra), so factorization through
   alpha_crit is automatic from Det[M]. *)
detExpanded = FullSimplify[detM, Assumptions -> $Assumptions];
t0Closed = FullSimplify[(a + dK)*kappa0^2 + a*kappa1^2, Assumptions -> $Assumptions];
expectZero["Det[M] - t0(alpha_crit - alpha)",
  FullSimplify[detExpanded - t0Closed*(alphaCritClosed - alpha), Assumptions -> $Assumptions]
];
expectZero["lam_- lam_+ - Det[M]",
  FullSimplify[lamMinusIndep*lamPlusIndep - detExpanded, Assumptions -> $Assumptions]
];

p0Factored = FullSimplify[
  beta0*sMinusOverlap*lamPlusIndep / (t0Closed*(alphaCritClosed - alpha)),
  Assumptions -> $Assumptions
];
expectZero["P0_- pole factorization",
  FullSimplify[p0SelIndep - p0Factored, Assumptions -> $Assumptions]
];
Print["P0_- factored = ", fmt[p0Factored]];

banner["STAGE 031 INDEPENDENT MATHEMATICA AUDIT COMPLETE"];
Print["Verified (independently from Eigenvalues/Eigenvectors of the loaded 2x2 wall matrix):"];
Print["  HF identity d lam_-/d alpha + s_- = 0 (derived, not assumed)"];
Print["  exact overlap derivative ds_-/dalpha = 2 dK^2 kappa0^2 kappa1^2 / R^3"];
Print["  exact prefactor derivative formula"];
Print["  initial values at alpha=0"];
Print["  alpha_crit from Det[M] = 0"];
Print["  pole factorization of P0_- at alpha_crit"];

Exit[0];
