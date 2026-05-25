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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 50] - N[target, 50]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS"];

Clear[zeta, piTr, cMix, epsBlk];
$Assumptions =
  Element[{zeta, piTr, cMix, epsBlk}, Reals] &&
  zeta >= 0 && piTr > 0 && cMix > 0 && epsBlk >= 0 && epsBlk < 1;

zetaExpr = FullSimplify[(piTr - cMix)/(cMix - epsBlk*(2*cMix - piTr)), Assumptions -> $Assumptions];
piOfZeta = FullSimplify[piTr /. First[Solve[zeta == zetaExpr, piTr]], Assumptions -> $Assumptions];
(* Strip ConditionalExpression wrapper Solve may have introduced; the underlying
   expression is the solution on the declared domain, and we need the bare form
   for downstream substitutions zeta -> {0, 1, ...} that would otherwise hit
   the ConditionalExpression's boundary and return Undefined. *)
piOfZeta = piOfZeta /. ConditionalExpression[e_, _] :> e;
qq = FullSimplify[piOfZeta / cMix, Assumptions -> $Assumptions];
qq = qq /. ConditionalExpression[e_, _] :> e;
expectZero["Q matches closed form",
  qq - (1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)];

Print["Pi_of_zeta = ", fmt[piOfZeta]];
Print["Q(zeta;eps_blk) = ", fmt[qq]];
expectZero["Q(0)-1", (qq /. zeta -> 0) - 1];
expectZero["Q(1)-2", (qq /. zeta -> 1) - 2];

zetaSuffChi1 = ToExpression["2.46622291347846`20"];
zetaFailChi1 = ToExpression["2.46752913273870`20"];
zetaSuffJ1 = ToExpression["2.44257571477179`20"];
zetaFailJ1 = ToExpression["2.46752736855058`20"];
zetaMaxF1 = ToExpression["2.46752922945601`20"];

piSuffChiOverC = FullSimplify[qq /. zeta -> zetaSuffChi1];
piFailChiOverC = FullSimplify[qq /. zeta -> zetaFailChi1];
piSuffJOverC = FullSimplify[qq /. zeta -> zetaSuffJ1];
piFailJOverC = FullSimplify[qq /. zeta -> zetaFailJ1];
piMaxOverC = FullSimplify[qq /. zeta -> zetaMaxF1];

Print["Pi_suff^(chi)/C_mix = ", fmt[piSuffChiOverC]];
Print["Pi_fail^(chi)/C_mix = ", fmt[piFailChiOverC]];
Print["Pi_suff^(J)/C_mix = ", fmt[piSuffJOverC]];
Print["Pi_fail^(J)/C_mix = ", fmt[piFailJOverC]];
Print["Pi_max^(F1)/C_mix = ", fmt[piMaxOverC]];

expectApprox["Pi_suff^(chi)/C_mix at eps=0 matches 1+zeta", N[(piSuffChiOverC - (1 + zetaSuffChi1)) /. epsBlk -> 0, 40], 0, 10^-14];
expectApprox["Pi_fail^(chi)/C_mix at eps=0 matches 1+zeta", N[(piFailChiOverC - (1 + zetaFailChi1)) /. epsBlk -> 0, 40], 0, 10^-14];
expectApprox["Pi_suff^(J)/C_mix at eps=0 matches 1+zeta", N[(piSuffJOverC - (1 + zetaSuffJ1)) /. epsBlk -> 0, 40], 0, 10^-14];
expectApprox["Pi_fail^(J)/C_mix at eps=0 matches 1+zeta", N[(piFailJOverC - (1 + zetaFailJ1)) /. epsBlk -> 0, 40], 0, 10^-14];
expectApprox["Pi_max^(F1)/C_mix at eps=0 matches 1+zeta", N[(piMaxOverC - (1 + zetaMaxF1)) /. epsBlk -> 0, 40], 0, 10^-14];

epsCeiling = N[1/zetaMaxF1, 40];
Print["Blocking ceiling eps_blk < ", fmt[epsCeiling]];
expectApprox["blocking ceiling reciprocal", N[epsCeiling*zetaMaxF1 - 1, 40], 0, 10^-14];

Print[""];
Print["Stage 081 Mathematica audit passed."];

Exit[0];
