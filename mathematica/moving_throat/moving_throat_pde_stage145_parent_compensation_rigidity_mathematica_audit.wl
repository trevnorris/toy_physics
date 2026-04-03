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

banner["STAGE 145 — PARENT COMPENSATION-SURFACE RIGIDITY"];

Clear[r, dr];
$Assumptions = Element[{r, dr}, Reals];

gamma0 = (1 + r^2)/9;
lRatio = Pi/2*Sqrt[(1 + r^2)/3];
gLower = r - Sqrt[1 + r^2]/2;

Print["gamma0(r) = ", fmt[gamma0]];
Print["L_W/a (r) = ", fmt[lRatio]];
Print["g_lower(r) = ", fmt[gLower]];

dlogGamma = FullSimplify[D[Log[gamma0], r]*dr];
dlogL = FullSimplify[D[Log[lRatio], r]*dr];
Print["d ln gamma0 = ", fmt[dlogGamma]];
Print["d ln(L_W/a) = ", fmt[dlogL]];
expectZero["similarity identity", dlogGamma - 2*dlogL];

slope = FullSimplify[D[gLower, r]];
Print["dg_lower/dr = ", fmt[slope]];
expectZero["lower-branch differential law", slope - (1 - r/(2*Sqrt[1 + r^2]))];

slopePos = FullSimplify[(4 + 3*r^2)/(2*Sqrt[1 + r^2]*(2*Sqrt[1 + r^2] + r))];
expectZero["positive slope decomposition", slope - slopePos];
Print["The lower-branch slope is manifestly positive for all real r."];

rF1 = SetPrecision[1.77799353547498, 30];
slopeNum = N[slope /. r -> rF1, 30];
invSlopeNum = N[1/slopeNum, 30];
Print["r_F1 = ", fmt[rF1]];
Print["(dg/dr)|_F1 = ", fmt[slopeNum]];
Print["dr/dg |_F1 = ", fmt[invSlopeNum]];

Print[""];
Print["Carry-forward formulas:"];
Print["  On the exact parent family: d ln gamma0 = 2 d ln(L_W/a), so Xi_slip = 0."];
Print["  On the lower branch: delta g = (dg/dr) delta r with dg/dr > 0."];
Print["  Therefore delta g = 0 implies delta r = 0, hence the first-order"];
Print["  D/N similarity defect and the reduced 2.5PN normalization defect vanish."];

Print[""];
Print["Stage 145 Mathematica audit passed."];

Exit[0];
