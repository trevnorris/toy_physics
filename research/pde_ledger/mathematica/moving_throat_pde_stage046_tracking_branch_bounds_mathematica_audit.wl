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

expectPositiveCoefficients[name_String, expr_, vars_List] := Module[{coeffs},
  coeffs = Last /@ CoefficientRules[Expand[expr], vars];
  Print[name, " coefficients = ", fmt[coeffs]];
  If[AllTrue[coeffs, Positive], pass[name], fail[name, coeffs]];
];

banner["STAGE 029 — TRACKING-BRANCH BOUNDS"];

Clear[xi, delta, r];
$Assumptions =
  Element[{xi, delta, r}, Reals] &&
  0 < xi < 1 && delta > 0 && r > 0;

gTr = FullSimplify[9*xi*(delta + xi)/(9*delta + (9 + 2*r^2)*xi), Assumptions -> $Assumptions];
fTr = FullSimplify[
  (9*delta + (9 + 2*r^2)*xi)^2*(9*delta + (9 + 2*r)*xi)^2/
    (81*(1 - xi)*(9*delta^2 + 18*delta*xi + (9 + 2*r^2)*xi^2)^2),
  Assumptions -> $Assumptions
];
gFlat = FullSimplify[gTr /. r -> 1, Assumptions -> $Assumptions];
fFlat = FullSimplify[fTr /. r -> 1, Assumptions -> $Assumptions];

Print["G_tr = ", fmt[gTr]];
Print["F_tr = ", fmt[fTr]];
Print["G_flat = ", fmt[gFlat]];
Print["F_flat = ", fmt[fFlat]];
expectZero["strong-split endpoint for G", FullSimplify[gTr /. r -> 0, Assumptions -> 0 < xi < 1 && delta > 0] - xi];
expectZero["strong-split endpoint for F", FullSimplify[fTr /. r -> 0, Assumptions -> 0 < xi < 1 && delta > 0] - 1/(1 - xi)];

(* Derivatives obtained directly from gTr/fTr; no hand-typed comparison polynomial *)
dGdR = Together[D[gTr, r]];
dFdR = Together[D[fTr, r]];
Print["dG_tr/dR = ", fmt[Factor[dGdR]]];
Print["dF_tr/dR = ", fmt[Factor[dFdR]]];

(* Sign of dG/dR on the bound domain: must be <= 0 on 0 <= r <= 1.
   Use Reduce on the closure (here we test strict negativity on the open interval). *)
dGSignClaim =
  Reduce[ForAll[{r, xi, delta}, 0 < r < 1 && 0 < xi < 1 && delta > 0, dGdR < 0], Reals];
Print["Reduce[dG/dR < 0 on (0,1)^3] = ", fmt[dGSignClaim]];
If[TrueQ[dGSignClaim], pass["dG/dR < 0 on bound domain"],
  fail["dG/dR < 0 on bound domain", dGSignClaim]];

dFSignClaim =
  Reduce[ForAll[{r, xi, delta}, 0 < r < 1 && 0 < xi < 1 && delta > 0, dFdR > 0], Reals];
Print["Reduce[dF/dR > 0 on (0,1)^3] = ", fmt[dFSignClaim]];
If[TrueQ[dFSignClaim], pass["dF/dR > 0 on bound domain"],
  fail["dF/dR > 0 on bound domain", dFSignClaim]];

(* Branch-difference forms derived directly; no hand-typed p1/p2/gDiffExpected/fDiffExpected *)
deltaG = Together[gTr - gFlat];
deltaF = Together[fFlat - fTr];
Print["G_tr - G_flat = ", fmt[Factor[deltaG]]];
Print["F_flat - F_tr = ", fmt[Factor[deltaF]]];

(* Verify the (1 - r^2) factor in deltaG via polynomial division *)
{deltaGQuotient, deltaGRemainder} =
  PolynomialQuotientRemainder[Together[Numerator[Factor[deltaG]]], 1 - r^2, r];
Print["deltaG numerator / (1 - r^2) quotient = ", fmt[Factor[deltaGQuotient]]];
Print["deltaG numerator / (1 - r^2) remainder = ", fmt[Factor[deltaGRemainder]]];
expectZero["(1 - r^2) divides numerator of G_tr - G_flat", deltaGRemainder];

(* Verify the (1 - r) factor in deltaF via polynomial division *)
{deltaFQuotient, deltaFRemainder} =
  PolynomialQuotientRemainder[Together[Numerator[Factor[deltaF]]], 1 - r, r];
Print["deltaF numerator / (1 - r) quotient = ", fmt[Factor[deltaFQuotient]]];
Print["deltaF numerator / (1 - r) remainder = ", fmt[Factor[deltaFRemainder]]];
expectZero["(1 - r) divides numerator of F_flat - F_tr", deltaFRemainder];

(* Sign of branch difference on the bound domain. Both should be > 0 on (0,1). *)
deltaGSignClaim =
  Reduce[ForAll[{r, xi, delta}, 0 < r < 1 && 0 < xi < 1 && delta > 0, deltaG > 0], Reals];
Print["Reduce[G_tr - G_flat > 0 on (0,1)^3] = ", fmt[deltaGSignClaim]];
If[TrueQ[deltaGSignClaim], pass["G_tr > G_flat on bound domain"],
  fail["G_tr > G_flat on bound domain", deltaGSignClaim]];

deltaFSignClaim =
  Reduce[ForAll[{r, xi, delta}, 0 < r < 1 && 0 < xi < 1 && delta > 0, deltaF > 0], Reals];
Print["Reduce[F_flat - F_tr > 0 on (0,1)^3] = ", fmt[deltaFSignClaim]];
If[TrueQ[deltaFSignClaim], pass["F_flat > F_tr on bound domain"],
  fail["F_flat > F_tr on bound domain", deltaFSignClaim]];

banner["3b. Boundary-value sign checks"];
expectZero["G_tr - G_flat vanishes at r=1", FullSimplify[(gTr - gFlat) /. r -> 1, Assumptions -> $Assumptions]];
expectZero["F_flat - F_tr vanishes at r=1", FullSimplify[(fFlat - fTr) /. r -> 1, Assumptions -> $Assumptions]];
expectZero["G_tr at r=0 equals xi", FullSimplify[(gTr /. r -> 0) - xi, Assumptions -> 0 < xi < 1 && delta > 0]];
expectZero["F_tr at r=0 equals 1/(1-xi)", FullSimplify[(fTr /. r -> 0) - 1/(1 - xi), Assumptions -> 0 < xi < 1 && delta > 0]];

samplePts = {
  {r -> 1/4, xi -> 1/3, delta -> 1/2},
  {r -> 1/2, xi -> 1/2, delta -> 1/4},
  {r -> 3/4, xi -> 1/5, delta -> 2/3}
};
Do[
  Module[{pt, gs, fs},
    pt = samplePts[[k]];
    gs = FullSimplify[(gTr - gFlat) /. pt];
    fs = FullSimplify[(fFlat - fTr) /. pt];
    Print["sample ", k, ": G_tr - G_flat = ", fmt[gs], ", F_flat - F_tr = ", fmt[fs]];
    If[!TrueQ[gs > 0], fail["sample " <> ToString[k] <> " violates G_tr > G_flat", gs]];
    If[!TrueQ[fs > 0], fail["sample " <> ToString[k] <> " violates F_tr < F_flat", fs]];
  ],
  {k, 1, Length[samplePts]}
];

Print[""];
Print["Stage 046 Mathematica audit passed."];

Exit[0];
