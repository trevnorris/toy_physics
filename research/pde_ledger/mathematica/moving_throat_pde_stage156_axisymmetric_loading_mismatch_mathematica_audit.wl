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

banner["STAGE 156 — WEAK-AXISYMMETRIC PHYSICAL-SLOPE TRANSPORT"];

Clear[eps, lam, d0, d01, d2, d21, d4, d41, n0, n01, u2, u4, p0];
$Assumptions = Element[{eps, lam, d0, d01, d2, d21, d4, d41, n0, n01, u2, u4, p0}, Reals] &&
  d0 != 0 && n0 != 0;

d0A = d0 + eps*lam*d01;
d2A = d2 + eps*lam*d21;
d4A = d4 + eps*lam*d41;
n0A = n0 + eps*lam*n01;

u2A = Expand[Normal[Series[-d2A/d0A, {eps, 0, 1}]]];
u4A = Expand[Normal[Series[(d2A^2 - d0A*d4A)/d0A^2, {eps, 0, 1}]]];
p0A = Expand[Normal[Series[n0A/d0A, {eps, 0, 1}]]];

u21 = FullSimplify[(D[u2A, eps] /. eps -> 0)/lam, Assumptions -> $Assumptions];
u41 = FullSimplify[(D[u4A, eps] /. eps -> 0)/lam, Assumptions -> $Assumptions];
p1 = FullSimplify[(D[p0A, eps] /. eps -> 0)/lam, Assumptions -> $Assumptions];

Print["u2^(1) general = ", fmt[u21]];
Print["u4^(1) general = ", fmt[u41]];
Print["P1 general     = ", fmt[p1]];

expectZero["u2 slope identity", u21 - (-(d21 + u2*d01)/d0 /. u2 -> -d2/d0)];

banner["Canonical branch formulas"];
u21Can = FullSimplify[u21 /. d2 -> -(1/9)*d0, Assumptions -> $Assumptions];
u41Can = FullSimplify[u41 /. {d2 -> -(1/9)*d0, d4 -> -(1/27)*d0}, Assumptions -> $Assumptions];
p1Ratio = FullSimplify[p1/(n0/d0), Assumptions -> $Assumptions];

Print["u2^(1) canonical = ", fmt[u21Can]];
Print["u4^(1) canonical = ", fmt[u41Can]];
Print["P1/P0            = ", fmt[p1Ratio]];

expectZero["u4 canonical formula", u41Can + (5*d01 + 18*d21 + 81*d41)/(81*d0)];
expectZero["P1/P0 formula", p1Ratio - (n01/n0 - d01/d0)];

banner["Hidden-even operator identity"];
hiddenEvenResidual = Expand[u41Can - 8*u21Can/9 - (d01/(27*d0) + 2*d21/(3*d0) - d41/d0)];
expectZero["hidden-even residual", hiddenEvenResidual];

d41Hidden = FullSimplify[d41 /. First[Solve[u41Can == 8*u21Can/9, d41]], Assumptions -> $Assumptions];
Print["D41 from hidden-even relation = ", fmt[d41Hidden]];

banner["Even-preserving collapse"];
u21ZeroD21 = FullSimplify[d21 /. First[Solve[u21Can == 0, d21]], Assumptions -> $Assumptions];
Print["D21 from u2^(1)=0 = ", fmt[u21ZeroD21]];

d41Even = FullSimplify[d41Hidden /. d21 -> u21ZeroD21, Assumptions -> $Assumptions];
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
Print["Carry-forward formulas:"];
Print["  u2^(1) = -(D21 + u2 D01)/D0"];
Print["  u4^(1) = -(5 D01 + 18 D21 + 81 D41)/(81 D0) on the canonical branch"];
Print["  P1/P0  = N01/N0 - D01/D0"];
Print["  hidden-even  <=>  D41 = 2 D21/3 + D01/27"];
Print["  if u2^(1)=0, then D21 = -D01/9 and D41 = -D01/27"];
Print["  remaining defect Xi_load = N01/N0 - D01/D0"];

Print[""];
Print["Stage 156 Mathematica audit passed."];

Exit[0];
