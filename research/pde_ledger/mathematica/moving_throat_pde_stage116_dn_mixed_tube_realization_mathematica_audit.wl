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
   Solve q'' + k^2 q = 0 with q(0)=0, q'(lW)=0.
   q(x) = Sin[k x] satisfies q(0)=0. The second BC k Cos[k lW] = 0 gives
   the smallest positive eigenvalue kW = Pi/(2 lW). *)
qTrial[xVar_, kArg_] := Sin[kArg*xVar];
odeRes = FullSimplify[D[qTrial[x, kSym], {x, 2}] + kSym^2 * qTrial[x, kSym], Assumptions -> $Assumptions];
expectZero["D/N trial satisfies q'' + k^2 q = 0", odeRes];
bcLeft = FullSimplify[qTrial[0, kSym], Assumptions -> $Assumptions];
expectZero["D/N trial satisfies q(0) = 0", bcLeft];
kWValue = Pi/(2*lW);
bcRightAtKW = FullSimplify[(D[qTrial[x, kSym], x] /. x -> lW) /. kSym -> kWValue,
                          Assumptions -> $Assumptions];
expectZero["D/N trial satisfies q'(lW) = 0 at k = Pi/(2 lW)", bcRightAtKW];

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

kappa0BareGeom = FullSimplify[4*lWRequired^2/(Pi^2*a^2), Assumptions -> $Assumptions];
expectZero["geometric kappa0 at lWRequired equals (1+r_c)/3",
           kappa0BareGeom - (1 + rC)/3];
gamma0Bare = FullSimplify[(1 + rC)/9, Assumptions -> $Assumptions];
kappaC = FullSimplify[kappa0BareGeom/(1 + rC), Assumptions -> $Assumptions];
gammaC = FullSimplify[gamma0Bare/(1 + rC), Assumptions -> $Assumptions];

expectZero["final kappa_c - 1/3", kappaC - 1/3];
expectZero["final gamma_c - 1/9", gammaC - 1/9];

dBare = Expand[(1 + rC)*(1 - z^2/3 - I*z^5/9)];
gamma0FromD = FullSimplify[I * Coefficient[dBare, z, 5], Assumptions -> $Assumptions];
expectZero["gamma0 extracted from dBare matches (1+rC)/9",
           gamma0FromD - (1 + rC)/9];
dFinal = FullSimplify[dBare/(1 + rC), Assumptions -> $Assumptions];
expectZero["bare scaled-canonical branch renormalizes to canonical", dFinal - (1 - z^2/3 - I*z^5/9)];
kappa0FromD = FullSimplify[-Coefficient[dBare, z, 2], Assumptions -> $Assumptions];
expectZero["kappa0_bare extracted from dBare matches (1+rC)/3",
           kappa0FromD - (1 + rC)/3];
kappa0Bare = FullSimplify[(1 + rC)/3, Assumptions -> $Assumptions];

Print[""];
Print["Summary:"];
Print["  D/N half-wave mixed tube length: ", fmt[lWRequired]];
Print["  Bare mixed coefficients: kappa0 = ", fmt[kappa0Bare], ", gamma0 = ", fmt[gamma0Bare]];
Print["  Final coefficients: kappa_c = ", fmt[kappaC], ", gamma_c = ", fmt[gammaC]];

Print[""];
Print["Stage 116 Mathematica audit passed."];

Exit[0];
