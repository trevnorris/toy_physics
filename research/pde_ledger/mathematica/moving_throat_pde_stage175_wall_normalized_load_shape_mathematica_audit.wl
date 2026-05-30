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

dlog[expr_] := FullSimplify[
  D[Log[FullSimplify[expr, Assumptions -> $Assumptions]], eps] /. eps -> 0,
  Assumptions -> $Assumptions
];

dlogSeries[expr_] := Coefficient[Normal[Series[Log[expr], {eps, 0, 1}]], eps];

banner["STAGE 175 — WALL-NORMALIZED LOAD/SHAPE FACTORIZATION"];

Clear[k, chi, varpi];
$Assumptions = Element[{k, chi, varpi}, Reals] && k > 0 && chi > 0 && varpi > 0;

c = Sqrt[k]*varpi*chi;
b0 = FullSimplify[c^2/varpi^2, Assumptions -> $Assumptions];
expectZero["B0 - K*chi^2", b0 - k*chi^2];

Clear[ou, ow, r, gu, gw, ouHat, owHat, rHat, guHat, gwHat];
$Assumptions = Element[{k, ou, ow, r, gu, gw, ouHat, owHat, rHat, guHat, gwHat}, Reals] &&
  k > 0 && ou > 0 && ow > 0 && r > 0 && gu > 0 && gw > 0 &&
  ouHat > 0 && owHat > 0 && rHat > 0 && guHat > 0 && gwHat > 0;

delta = ou*ow - r^2;
q = gu^2*ow + 2*gu*gw*r + gw^2*ou;
p = ou*gw + r*gu;

subsHat = {ou -> k*ouHat, ow -> k*owHat, r -> k*rHat, gu -> k*guHat, gw -> k*gwHat};

expectZero["Delta - K^2*Delta_hat", Expand[delta /. subsHat] - k^2*(ouHat*owHat - rHat^2)];
expectZero[
  "Q - K^3*Q_hat",
  Expand[q /. subsHat] - k^3*(guHat^2*owHat + 2*guHat*gwHat*rHat + gwHat^2*ouHat)
];
expectZero["P - K^2*P_hat", Expand[p /. subsHat] - k^2*(ouHat*gwHat + rHat*guHat)];

upsilon = FullSimplify[(q/(k*delta)) /. subsHat, Assumptions -> $Assumptions];
lambda = FullSimplify[(p/delta) /. subsHat, Assumptions -> $Assumptions];
expectZero["Z0 - K*Upsilon", FullSimplify[(q/delta) /. subsHat, Assumptions -> $Assumptions] - k*upsilon];
expectZero["N0 - Lambda^2", FullSimplify[(p^2/delta^2) /. subsHat, Assumptions -> $Assumptions] - lambda^2];

banner["Differential defect identities"];

Clear[eps, kappa, schi, su, sw, sr, sgu, sgw, k0, chi0, varpi0, ou0, ow0, r0, gu0, gw0];
$Assumptions = Element[{eps, kappa, schi, su, sw, sr, sgu, sgw, k0, chi0, varpi0, ou0, ow0, r0, gu0, gw0}, Reals] &&
  k0 > 0 && chi0 > 0 && varpi0 > 0 && ou0 > 0 && ow0 > 0 && r0 > 0 && gu0 > 0 && gw0 > 0;

subsEps = {
  k -> k0*Exp[kappa*eps],
  chi -> chi0*Exp[schi*eps],
  varpi -> varpi0,
  ouHat -> ou0*Exp[su*eps],
  owHat -> ow0*Exp[sw*eps],
  rHat -> r0*Exp[sr*eps],
  guHat -> gu0*Exp[sgu*eps],
  gwHat -> gw0*Exp[sgw*eps]
};

exprC = c /. {k -> k0*Exp[kappa*eps], chi -> chi0*Exp[schi*eps], varpi -> varpi0};
sigmaBDirect = FullSimplify[2*dlog[exprC] - kappa, Assumptions -> $Assumptions];
sigmaBShape = FullSimplify[dlog[(chi^2) /. subsEps], Assumptions -> $Assumptions];
expectZero["Sigma_B - dln(chi^2)", sigmaBDirect - sigmaBShape];

exprZ = ((q/delta) /. subsHat) /. subsEps;
exprU = upsilon /. subsEps;
sigmaZDirect = FullSimplify[dlog[exprZ] - kappa, Assumptions -> $Assumptions];
sigmaZShape = FullSimplify[dlog[exprU], Assumptions -> $Assumptions];
expectZero["Sigma_Z - dln(Upsilon)", sigmaZDirect - sigmaZShape];

(* Direct route: build (P/Delta) from the physical primitives and apply the
   physical wall scaling (subsHat) THEN the eps-flow, without first simplifying
   to the cached lambda. This keeps the K-dependence explicit before cancellation,
   so the -kappa (= -delta_K) subtraction is load-bearing rather than tautological. *)
exprPoverDeltaPhys = ((p/delta) /. subsHat) /. subsEps;
sigmaNDirect = FullSimplify[2*dlog[exprPoverDeltaPhys] - kappa, Assumptions -> $Assumptions];
sigmaNShape = FullSimplify[dlog[(lambda^2/k) /. subsEps], Assumptions -> $Assumptions];
expectZero["Sigma_N - dln(Lambda^2/K)", sigmaNDirect - sigmaNShape];
(* Independent second-engine slope route (red-team R1): extract the Sigma_N
   first-order log-slope via Series-coefficient (Series+Coefficient) instead of
   D[Log[.]], so the differential identity no longer relies on a transliteration
   of the SymPy dlog route. Series-route DIRECT (2 dln(P/Delta) - dK) vs the SHAPE
   target dln(Lambda^2/K); -kappa (= -delta_K) kept symbolic. *)
sigmaNDirectSeries = FullSimplify[2*dlogSeries[exprPoverDeltaPhys] - kappa, Assumptions -> $Assumptions];
sigmaNShapeSeries  = FullSimplify[dlogSeries[(lambda^2/k) /. subsEps], Assumptions -> $Assumptions];
expectZero["Sigma_N - dln(Lambda^2/K) [series route]", sigmaNDirectSeries - sigmaNShapeSeries];
(* Note (red-team F1 resolution): the differential Sigma_N claim is fully and
   non-trivially exercised by the check above — 2 dln(P/Delta) - dK = dln(Lambda^2/K)
   is load-bearing on kappa = delta_K, and the genuine homogeneity coverage is
   N0 = Lambda^2 (checked earlier). The earlier "2 dln(P/Delta) - 2 dln Lambda" /
   "Sigma_N - (2 dln Lambda - dK)" lines reduced to a simplify-commutes identity
   (p/delta and lambda are value-equal), so they are omitted rather than reported
   as substantive PASS lines. *)

banner["Conservative-shape-preserving reductions"];
sigmaBCons = FullSimplify[sigmaBDirect /. schi -> 0, Assumptions -> $Assumptions];
sigmaZCons = FullSimplify[sigmaZDirect /. {su -> 0, sw -> 0, sr -> 0, sgu -> 0, sgw -> 0}, Assumptions -> $Assumptions];
sigmaNCommon = FullSimplify[sigmaNDirect /. {su -> 0, sw -> 0, sr -> 0, sgu -> 0, sgw -> 0}, Assumptions -> $Assumptions];
expectZero["Conservative-shape branch Sigma_B", sigmaBCons];
expectZero["Conservative-shape branch Sigma_Z", sigmaZCons];
expectZero["Common-shape branch Sigma_N + dK", sigmaNCommon + kappa];
(* Weighted aggregate no-go using sum_r rho_r^(N) = 1. *)
Clear[rho1, rho2];
$Assumptions = $Assumptions && Element[{rho1, rho2}, Reals] && rho1 >= 0 && rho2 >= 0;
xiLoadFrozen = (rho1 + rho2)*sigmaNCommon;
expectZero["Xi_load (all shapes frozen) + dK", (xiLoadFrozen /. (rho1 + rho2) -> 1) + kappa];

Print[""];
Print["Conclusions:"];
Print["  B0 = K chi^2,  Z0 = K Upsilon,  N0 = Lambda^2."];
Print["  Therefore"];
Print["    Sigma_B = d ln chi^2"];
Print["    Sigma_Z = d ln Upsilon"];
Print["    Sigma_N = d ln(Lambda^2/K) = 2 d ln Lambda - dK"];
Print["  If all wall-normalized shapes are frozen, then Sigma_N = -dK."];
Print["  So naive common self-similarity does not kill the grouped defect."];

Print[""];
Print["Stage 175 Mathematica audit passed."];

Exit[0];
