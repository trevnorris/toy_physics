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

banner["STAGE 116 — FINITE D/N MIXED-TUBE REALIZATION"];

Clear[a, lW, kS, kQ, lam, z, kSym, omegaSym, csSym, x];
$Assumptions =
  Element[{a, lW, kS, kQ, lam, z, kSym, omegaSym, csSym, x}, Reals] &&
  a > 0 && lW > 0 && kS > 0 && kQ > 0 && lam > 0 &&
  kSym > 0 && omegaSym > 0 && csSym > 0;

(* D/N half-wave eigenvalue derivation:
   Solve q'' + k^2 q = 0 with q(0)=0, then impose q'(lW)=0.
   The normalized nonzero solution branch gives the characteristic equation
   Cos[k lW] == 0; solve it for the smallest positive eigenvalue. *)
gensol = DSolve[{q''[xv] + kSym^2*q[xv] == 0, q[0] == 0, q'[0] == 1}, q, xv];
qGenExpr = FullSimplify[q[xv] /. First[gensol], Assumptions -> $Assumptions];
odeRes = FullSimplify[D[qGenExpr, {xv, 2}] + kSym^2*qGenExpr, Assumptions -> $Assumptions];
expectZero["D/N solved mode satisfies q'' + k^2 q = 0", odeRes];
bcLeft = FullSimplify[qGenExpr /. xv -> 0, Assumptions -> $Assumptions];
expectZero["D/N solved mode satisfies q(0) = 0", bcLeft];
charEq = FullSimplify[D[qGenExpr, xv] /. xv -> lW, Assumptions -> $Assumptions];
(* Solve for the product u = kW*lW on (0, Pi), avoiding symbolic division by lW. *)
uRoot = u /. First[Solve[Cos[u] == 0 && 0 < u < Pi, u, Reals]];
kWValue = FullSimplify[uRoot/lW, Assumptions -> $Assumptions];
expectZero["D/N eigenvalue solves Cos[kW lW]==0", Cos[kWValue*lW]];
expectZero["D/N eigenvalue kW = Pi/(2 lW)", kWValue - Pi/(2*lW)];

OmegaW = kWValue * csSym;
kappa0Derived = FullSimplify[(omegaSym/OmegaW)^2 / (a*omegaSym/csSym)^2,
                             Assumptions -> $Assumptions];
expectZero["kappa0 from D/N eigenvalue matches geometric expression",
           kappa0Derived - 4*lW^2/(Pi^2*a^2)];

rC = FullSimplify[lam^2/(kS*kQ), Assumptions -> $Assumptions];
kappa0FromTube = FullSimplify[kappa0Derived, Assumptions -> $Assumptions];
lWRequired = FullSimplify[lW /. First[Solve[kappa0FromTube == (1 + rC)/3, lW, Reals]], Assumptions -> $Assumptions];

Print["kappa0 from D/N half-wave tube = ", fmt[kappa0FromTube]];
Print["Required tube length L_W = ", fmt[lWRequired]];
expectZero["tube-length law", lWRequired - (Pi*a*Sqrt[(1 + rC)/3])/2];

(* --- Renormalization to canonical coefficients (REPORTED, not asserted) --- *)
(* Load-bearing physics verified above; gamma0 is an upstream-carried input    *)
(* (Stage 98), so kappa_c/gamma_c here are definitional consequences, printed  *)
(* not asserted (an expectZero here would be tautological).                    *)
kappa0Bare = FullSimplify[4*lWRequired^2/(Pi^2*a^2), Assumptions -> $Assumptions];  (* derived tube coeff at required length *)
gamma0Bare = FullSimplify[(1 + rC)/9, Assumptions -> $Assumptions];                  (* upstream-carried input (Stage 98) *)
commonScale = 1 + rC;
kappaC = FullSimplify[kappa0Bare/commonScale, Assumptions -> $Assumptions];
gammaC = FullSimplify[gamma0Bare/commonScale, Assumptions -> $Assumptions];
Print["Renormalization (definitional consequence, not an independent check):"];
Print["  kappa0_bare (derived tube coeff at lWRequired) = ", kappa0Bare];
Print["  gamma0_bare (upstream-carried input, Stage 98) = ", gamma0Bare];
Print["  kappa_c = kappa0Bare/(1+rC) = ", kappaC];
Print["  gamma_c = gamma0Bare/(1+rC) = ", gammaC];

Print[""];
Print["Summary:"];
Print["  D/N half-wave mixed tube length: ", fmt[lWRequired]];
Print["  Bare mixed coefficients: kappa0 = ", fmt[kappa0Bare], ", gamma0 = ", fmt[gamma0Bare]];
Print["  Final coefficients: kappa_c = ", fmt[kappaC], ", gamma_c = ", fmt[gammaC]];

Print[""];
Print["Stage 116 Mathematica audit passed."];

Exit[0];
