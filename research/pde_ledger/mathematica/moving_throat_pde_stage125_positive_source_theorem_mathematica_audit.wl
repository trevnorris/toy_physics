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
  res = FullSimplify[Expand[expr], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 125 — POSITIVE LOCAL MOUTH-SOURCE THEOREM"];

Clear[x, z, L];
$Assumptions = L > 0 && Element[{x, z}, Reals];

k = Pi/(2*L);
kernel = Cos[Pi*x/2];

Print["D/N half-wave derivative factor:"];
Print["Cos[k z] with k = Pi/(2 L) = ", fmt[Cos[k*z]]];

kernelMin = FullSimplify[MinValue[{kernel, 0 <= x <= 1}, x], Assumptions -> $Assumptions];
kernelMax = FullSimplify[MaxValue[{kernel, 0 <= x <= 1}, x], Assumptions -> $Assumptions];
expectZero["kernel minimum on [0,1]", kernelMin];
expectZero["kernel maximum on [0,1] - 1", kernelMax - 1];

Print[""];
Print["Because the normalized kernel lies in [0,1] on x in [0,1],"];
Print["every positive normalized source law has its cosine moment in [0,1]."];

rrad = Sqrt[4107 - 100*Pi^2];
r = FullSimplify[Sqrt[12*(37/20)^2/Pi^2 - 1], Assumptions -> $Assumptions];

Clear[gSym];
branchSolutions = Solve[1 + r^2 - 4*(gSym - r)^2 == 0, gSym];
branchValues = FullSimplify[gSym /. branchSolutions, Assumptions -> $Assumptions];
gminus = First[Sort[branchValues]];
gplus  = Last[Sort[branchValues]];

expectZero["r - sqrt(4107 - 100 Pi^2)/(10 Pi)", r - rrad/(10*Pi)];
expectZero["g_- matches closed form (2*rrad - 37*Sqrt[3])/(20*Pi)",
  gminus - (2*rrad - 37*Sqrt[3])/(20*Pi)];
expectZero["g_+ matches closed form (2*rrad + 37*Sqrt[3])/(20*Pi)",
  gplus - (2*rrad + 37*Sqrt[3])/(20*Pi)];

Print[""];
Print["Explicit Family-1 compensated branches:"];
Print["g_- = ", fmt[gminus]];
Print["g_+ = ", fmt[gplus]];
Print["g_- (numeric) = ", N[gminus, 20]];
Print["g_+ (numeric) = ", N[gplus, 20]];

expectTrue["g_- > 0", gminus > 0];
expectTrue["g_- < 1", gminus < 1];
expectTrue["g_+ > 1", gplus > 1];

(* Integral-bound test on a one-parameter family of positive normalized sources.
   sigmaA(z) = (a + 1) * (z/L)^a / L is nonneg on [0,L] for a >= 0 and integrates to 1.
   Endpoint values: a = 0 -> uniform (g = 2/Pi); a -> oo -> peaked at z = L (g -> 0). *)
Clear[aSym, zSym];
$Assumptions = $Assumptions && aSym >= 0;
sigmaProfile = (aSym + 1)*(zSym/L)^aSym / L;
normA = FullSimplify[Integrate[sigmaProfile, {zSym, 0, L}], Assumptions -> $Assumptions];
expectZero["parametric family normalization", normA - 1];

gA = FullSimplify[Integrate[sigmaProfile*Cos[Pi*zSym/(2*L)], {zSym, 0, L}],
                  Assumptions -> $Assumptions];
Print["g_a (parametric moment) = ", fmt[gA]];

expectZero["moment g[uniform] - 2/Pi", (gA /. aSym -> 0) - 2/Pi];
expectZero["moment g[peaked@L] limit", Limit[gA, aSym -> Infinity]];
gApeaked = Chop[Re[Block[{$MaxExtraPrecision = 10000}, N[(gA /. aSym -> 100), 30]]]];
Print["g_a at aSym -> 100 (peaked-at-L proxy) = ", fmt[gApeaked]];
expectTrue["g[peaked@L proxy a=100] >= 0", gApeaked >= 0];
expectTrue["g[peaked@L proxy a=100] <= 1", gApeaked <= 1];

expectTrue["g[uniform] >= 0", (gA /. aSym -> 0) >= 0];
expectTrue["g[uniform] <= 1", (gA /. aSym -> 0) <= 1];

Print[""];
Print["Conclusion: under any positive localized mouth source law,"];
Print["the upper compensated branch is impossible and the lower branch is unique."];

Exit[0];
