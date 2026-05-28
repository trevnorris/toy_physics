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

expectApprox[name_String, value_, target_, tol_] := Module[{delta},
  delta = N[value - target, 40];
  Print[name, " residual = ", fmt[delta]];
  If[TrueQ[Abs[delta] < tol], pass[name], fail[name, delta]];
];

banner["STAGE 127 — GEOMETRIC MOUTH-PENETRATION FAMILIES"];

Clear[x, z];
$Assumptions = Element[x, Reals] && x > 0 && Element[z, Reals];

rDisc = Sqrt[4107 - 100*Pi^2];
gMinus = N[(2*rDisc - 37*Sqrt[3])/(20*Pi), 80];

(* Independent derivation of g_slab from the slab source integral.
   With L = 1 (the bias factor is L-independent after u = z/L). *)
gSlabFromIntegral = FullSimplify[Integrate[(1/x)*Cos[Pi*z/2], {z, 0, x}], Assumptions -> x > 0];
gSlab = FullSimplify[2*Sin[Pi*x/2]/(Pi*x), Assumptions -> $Assumptions];
slabClosedFormResidual = FullSimplify[gSlabFromIntegral - gSlab, Assumptions -> x > 0];
If[TrueQ[slabClosedFormResidual === 0],
   pass["slab closed-form matches source integral"],
   fail["slab closed-form matches source integral", slabClosedFormResidual]];

(* Independent derivation of g_exp from the exponential source integral. *)
gExpFromIntegral = FullSimplify[Integrate[(Exp[-z/x]/(x*(1 - Exp[-1/x])))*Cos[Pi*z/2], {z, 0, 1}], Assumptions -> x > 0];
gExp = FullSimplify[2*(2 + Pi*x*Exp[-1/x])/((4 + Pi^2*x^2)*(1 - Exp[-1/x])), Assumptions -> $Assumptions];
expClosedFormResidual = FullSimplify[gExpFromIntegral - gExp, Assumptions -> x > 0];
If[TrueQ[expClosedFormResidual === 0],
   pass["exp closed-form matches source integral"],
   fail["exp closed-form matches source integral", expClosedFormResidual]];

Print["g_slab(x) = ", fmt[gSlab]];
Print["g_exp(x) = ", fmt[gExp]];

xStarSlab = x /. FindRoot[gSlab == gMinus, {x, 0.8}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 200];
xStarExp = x /. FindRoot[gExp == gMinus, {x, 0.66}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 200];

Print["g_-^F1 = ", fmt[N[gMinus, 30]]];
Print["x_*^slab = ", fmt[N[xStarSlab, 30]]];
Print["x_*^exp  = ", fmt[N[xStarExp, 30]]];

expectApprox["slab compensation root", gSlab /. x -> xStarSlab, gMinus, 10^-20];
expectApprox["exponential compensation root", gExp /. x -> xStarExp, gMinus, 10^-20];

Print[""];
Print["Stage 127 Mathematica audit passed."];

Exit[0];
