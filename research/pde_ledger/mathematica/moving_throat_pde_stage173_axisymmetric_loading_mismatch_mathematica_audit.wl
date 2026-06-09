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

banner["STAGE 173 — WEAK-AXISYMMETRIC PHYSICAL-SLOPE TRANSPORT"];

Clear[eps, lam, d0, d01, d2, d21, d4, d41, n0, n01, u2, u4, p0];
$Assumptions = Element[{eps, lam, d0, d01, d2, d21, d4, d41, n0, n01, u2, u4, p0}, Reals] &&
  d0 != 0 && n0 != 0;

u21 = FullSimplify[
  (D[-(d2 + eps*lam*d21)/(d0 + eps*lam*d01), eps] /. eps -> 0)/lam,
  Assumptions -> $Assumptions];
u41 = FullSimplify[
  (D[
      ((d2 + eps*lam*d21)^2 - (d0 + eps*lam*d01)*(d4 + eps*lam*d41))/
        (d0 + eps*lam*d01)^2,
      eps] /. eps -> 0)/lam,
  Assumptions -> $Assumptions];
p1 = FullSimplify[
  (D[(n0 + eps*lam*n01)/(d0 + eps*lam*d01), eps] /. eps -> 0)/lam,
  Assumptions -> $Assumptions];

Print["u2^(1) general = ", fmt[u21]];
Print["u4^(1) general = ", fmt[u41]];
Print["P1 general     = ", fmt[p1]];

expectZero["u2 slope identity", u21 - (-(d21 + u2*d01)/d0 /. u2 -> -d2/d0)];

banner["Canonical branch formulas"];
canonicalBranchRules = {d2 -> -d0/9, d4 -> -d0/27};
branchU2Slope = FullSimplify[u21 /. {d2 -> -d0/9}, Assumptions -> $Assumptions];
branchU4Slope = FullSimplify[u41 /. canonicalBranchRules, Assumptions -> $Assumptions];
staticPressureRatio = FullSimplify[p1/(n0/d0), Assumptions -> $Assumptions];

Print["u2^(1) canonical = ", fmt[branchU2Slope]];
Print["u4^(1) canonical = ", fmt[branchU4Slope]];
Print["P1/P0            = ", fmt[staticPressureRatio]];

expectZero["u4 canonical formula", branchU4Slope + (5*d01 + 18*d21 + 81*d41)/(81*d0)];
expectZero["P1/P0 formula", staticPressureRatio - (n01/n0 - d01/d0)];

banner["Hidden-even operator identity"];
hiddenEvenResidual = Expand[branchU4Slope - 8*branchU2Slope/9 - (d01/(27*d0) + 2*d21/(3*d0) - d41/d0)];
expectZero["hidden-even residual", hiddenEvenResidual];

banner["Even-preserving collapse"];
u2BalanceNumerator = Numerator[Together[branchU2Slope]];
u21ZeroD21 = FullSimplify[
  -(u2BalanceNumerator /. d21 -> 0)/Coefficient[u2BalanceNumerator, d21],
  Assumptions -> $Assumptions];
Print["D21 from u2^(1)=0 = ", fmt[u21ZeroD21]];

d41BalanceNumerator = Numerator[Together[(branchU4Slope - 8*branchU2Slope/9) /. d21 -> u21ZeroD21]];
d41Even = FullSimplify[
  -(d41BalanceNumerator /. d41 -> 0)/Coefficient[d41BalanceNumerator, d41],
  Assumptions -> $Assumptions];
Print["D41 on even-preserving branch = ", fmt[d41Even]];

expectZero["D21 + D01/9", u21ZeroD21 + d01/9];
expectZero["D41 + D01/27", d41Even + d01/27];

banner["Single static loading mismatch"];
xiLoad = FullSimplify[n01/n0 - d01/d0, Assumptions -> $Assumptions];
Print["Xi_load = ", fmt[xiLoad]];

lam20 = 1;
lam21 = 1/2;
lam22 = -1;
Print["Delta_Q^(20)/eps = ", fmt[FullSimplify[lam20*xiLoad, Assumptions -> $Assumptions]]];
Print["Delta_Q^(21)/eps = ", fmt[FullSimplify[lam21*xiLoad, Assumptions -> $Assumptions]]];
Print["Delta_Q^(22)/eps = ", fmt[FullSimplify[lam22*xiLoad, Assumptions -> $Assumptions]]];

Print[""];
Print["Mathematica summary: verified slope identities, even-preserving collapse, and Xi_load lanes above."];

Print[""];
Print["Stage 173 Mathematica audit passed."];

Exit[0];
