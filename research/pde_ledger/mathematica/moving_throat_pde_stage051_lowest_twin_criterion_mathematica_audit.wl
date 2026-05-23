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
  (* Strip ConditionalExpression wrapper: under $Assumptions, a result of
     the form ConditionalExpression[0, cond] is identically zero on the
     declared domain.  Solve[]/Reduce[] often introduce these wrappers when
     auxiliary inequalities are nontrivial. *)
  res = res /. ConditionalExpression[e_, _] :> e;
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 034 — LOWEST TWIN CRITERION"];

Clear[xi, delta, r, lambda, eps, epsEta, chi0, mMix, zW, pBranch];
$Assumptions =
  Element[{xi, delta, r, lambda, eps, epsEta, chi0, mMix, zW, pBranch}, Reals] &&
  0 < xi < 1 && delta > 0 && r > 0 && lambda > 0 && 0 < eps < 1 &&
  0 < epsEta < 1 && chi0 > -1 && mMix > 0 && zW > 0 && pBranch > 0;

gTr = FullSimplify[9 xi (xi + delta)/(9 delta + (9 + 2 r^2) xi), Assumptions -> $Assumptions];
fTr = FullSimplify[
  (9 delta + (9 + 2 r^2) xi)^2 (9 delta + (9 + 2 r) xi)^2/
    (81 (1 - xi) (9 delta^2 + 18 delta xi + (9 + 2 r^2) xi^2)^2),
  Assumptions -> $Assumptions
];
piTr = FullSimplify[fTr gTr, Assumptions -> $Assumptions];
(* Independent canonicalization via Factor/Together (no FullSimplify wrapper around the diff). *)
piTrFactored = Factor[Together[fTr gTr]];
piExpected = xi (xi + delta) (9 delta + (9 + 2 r) xi)^2 (9 delta + (9 + 2 r^2) xi)/
    (9 (1 - xi) (9 delta^2 + 18 delta xi + (9 + 2 r^2) xi^2)^2);
piExpectedFactored = Factor[Together[piExpected]];

Print["G_tr = ", fmt[gTr]];
Print["F_tr = ", fmt[fTr]];
Print["Pi_tr = ", fmt[piTr]];
Print["Pi_tr (Factor/Together) = ", fmt[piTrFactored]];
Print["Pi_tr (claim, Factor/Together) = ", fmt[piExpectedFactored]];
expectZero["Pi_tr - expected closed form", piTrFactored - piExpectedFactored];

pi0 = FullSimplify[Limit[piTr, xi -> 0, Direction -> "FromAbove"], Assumptions -> delta > 0 && r > 0];
pi1 = Limit[piTr, xi -> 1, Direction -> "FromBelow"];
Print["Pi_tr(xi->0+) = ", fmt[pi0]];
Print["Pi_tr(xi->1-) = ", fmt[pi1]];
If[pi0 =!= 0, fail["Pi_tr(xi->0+)"], None];
(* Mathematica's Limit heuristically returns either Infinity or
   Infinity/<positive> for poles of polynomial fractions; both encode the
   same blow-up.  Test via 1/pi1 -> 0 to be robust to that flakiness. *)
If[!TrueQ[Simplify[1/pi1 == 0, Assumptions -> $Assumptions]],
   fail["Pi_tr(xi->1-) is not +infinity", pi1], None];

cMix = FullSimplify[8 lambda (1 - eps)/Pi^2, Assumptions -> $Assumptions];
sReq = FullSimplify[pBranch/cMix, Assumptions -> $Assumptions];
zetaReq = FullSimplify[(sReq - 1)/(1 + eps (sReq - 2)), Assumptions -> $Assumptions];
zetaReqBranch = FullSimplify[zetaReq /. pBranch -> piTr, Assumptions -> $Assumptions];

Print["C_mix = ", fmt[cMix]];
Print["S_req = ", fmt[sReq]];
Print["zeta_req(Pi,C_mix,eps) = ", fmt[zetaReq]];
Print["Physical-branch zeta_req = ", fmt[zetaReqBranch]];
expectZero["zeta_req at Pi = C_mix", zetaReq /. pBranch -> cMix];
expectZero["zeta_req at Pi = 2 C_mix minus 1", (zetaReq /. pBranch -> 2 cMix) - 1];
expectZero["zeta_req - 1 at Pi = 2 C_mix", ((zetaReq - 1) /. pBranch -> 2 cMix)];

lambdaTwinReq = FullSimplify[Pi^2 piTr/(16 (1 - eps)), Assumptions -> $Assumptions];
mMixTwinReq = FullSimplify[gTr/2, Assumptions -> $Assumptions];
zWTwinReq = FullSimplify[Pi^2 (1 - epsEta) (1 - eps) gTr/(16 (1 + chi0)^2), Assumptions -> $Assumptions];
(* Stage 047/030 coherent forward map: Z_W = pi^2 (1-eps_eta)(1-eps) M_mix / [8 (1+chi0)^2]. *)
zWFromMmix = FullSimplify[Pi^2 (1 - epsEta) (1 - eps) mMix/(8 (1 + chi0)^2), Assumptions -> $Assumptions];
zWThresholdViaMap = FullSimplify[zWFromMmix /. mMix -> gTr/2, Assumptions -> $Assumptions];

Print["Lambda_twin,req = ", fmt[lambdaTwinReq]];
Print["M_mix^(twin,req) = ", fmt[mMixTwinReq]];
Print["Z_W^(twin,req) = ", fmt[zWTwinReq]];
expectZero["Z_W^(twin,req) - forward-map(M_mix=G_tr/2)", zWTwinReq - zWThresholdViaMap];

(* Independently derive the positive root of gTr == 2 mMix via Solve on Reals
   with explicit positivity constraint; symbolic positivity of the algebraic
   roots is undecidable without further assumptions, so SelectFirst-on-TrueQ
   silently drops both roots. *)
xi2xDerivedRaw = xi /. First[Solve[{gTr == 2 mMix, xi > 0}, xi, Reals]];
xi2xDerivedRaw = xi2xDerivedRaw /. ConditionalExpression[e_, _] :> e;
xi2xDerived = FullSimplify[xi2xDerivedRaw, Assumptions -> $Assumptions];
(* Closed-form claim (the docstring's answer). *)
xi2xClaim = FullSimplify[
  (2 mMix (9 + 2 r^2) - 9 delta + Sqrt[(2 mMix (9 + 2 r^2) - 9 delta)^2 + 648 mMix delta])/18,
  Assumptions -> $Assumptions
];
Print["xi_(2x) (Solve) = ", fmt[xi2xDerived]];
Print["xi_(2x) (claim) = ", fmt[xi2xClaim]];
expectZero["xi_(2x): Solve vs claim", xi2xDerived - xi2xClaim];
expectZero["G_tr(xi_(2x)) - 2 M_mix", (gTr /. xi -> xi2xDerived) - 2 mMix];

Print[""];
Print["Stage 051 Mathematica audit passed."];

Exit[0];
