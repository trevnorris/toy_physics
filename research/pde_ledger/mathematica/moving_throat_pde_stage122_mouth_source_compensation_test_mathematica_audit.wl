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

stripCE[expr_] := expr /. ConditionalExpression[v_, _] :> v;

reduceExact[expr_] := Module[{res},
  res = stripCE[expr];
  res = FullSimplify[RootReduce[Cancel[Together[res]]], Assumptions -> $Assumptions];
  stripCE[res]
];

expectZero[name_String, expr_] := Module[{res},
  res = reduceExact[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectNonzero[name_String, expr_] := Module[{res, nonzeroQ, numericAwayQ},
  res = reduceExact[expr];
  Print[name, " = ", fmt[res]];
  Print[name, " numeric = ", fmt[N[res, 20]]];
  nonzeroQ = TrueQ[FullSimplify[res != 0, Assumptions -> $Assumptions]];
  numericAwayQ = TrueQ[Abs[N[res, 40]] > 10^-6];
  If[nonzeroQ && numericAwayQ, pass[name], fail[name, res]];
];

banner["STAGE 122 - MOUTH SOURCE COMPENSATION TEST"];

Clear[g, r, gNat, cStage];

$Assumptions =
  Element[{g, r, gNat, cStage}, Reals] &&
  cStage > 0 && gNat > 0 &&
  4107 - 100*Pi^2 > 0;

rFromGeometry[radius_] := Sqrt[12*radius^2/Pi^2 - 1];
compensationResidual[rval_, gain_] := 1 + rval^2 - 4*(gain - rval)^2;

mouthRadius = 37/20;
rF1 = reduceExact[rFromGeometry[mouthRadius]];

Print["r_F1 = ", fmt[rF1]];
expectZero["r_F1 geometric relation", rF1^2 - (12*mouthRadius^2/Pi^2 - 1)];

compensationEquation = compensationResidual[rF1, g] == 0;
rootRules = Solve[compensationEquation, g, Reals];
If[Length[rootRules] =!= 2, fail["expected two compensated branches", rootRules]];

roots = reduceExact /@ (g /. rootRules);
gMinus = SelectFirst[
  roots,
  TrueQ[FullSimplify[# - rF1 < 0, Assumptions -> $Assumptions]] &
];
gPlus = SelectFirst[
  roots,
  TrueQ[FullSimplify[# - rF1 > 0, Assumptions -> $Assumptions]] &
];

If[MissingQ[gMinus] || MissingQ[gPlus], fail["branch selection by sign around r_F1", roots]];

Print["derived g_minus = ", fmt[gMinus]];
Print["derived g_plus  = ", fmt[gPlus]];
Print["numeric g_minus = ", fmt[N[gMinus, 20]]];
Print["numeric g_plus  = ", fmt[N[gPlus, 20]]];

gMinusExact = (2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi);
gPlusExact = (2*Sqrt[4107 - 100*Pi^2] + 37*Sqrt[3])/(20*Pi);

expectZero["gminus exact form", gMinus - gMinusExact];
expectZero["gplus exact form", gPlus - gPlusExact];

expectZero["compensation quadratic at gminus", compensationResidual[rF1, gMinus]];
expectZero["compensation quadratic at gplus", compensationResidual[rF1, gPlus]];

gNatAnsatz = 1;
naturalDefect = reduceExact[compensationResidual[rF1, gNat] /. gNat -> gNatAnsatz];
defectExact = (-12321 + 80*Pi*Sqrt[4107 - 100*Pi^2])/(100*Pi^2);

Print["natural branch defect = ", fmt[naturalDefect]];
Print["natural branch defect numeric = ", fmt[N[naturalDefect, 20]]];

expectZero["defect closed form", naturalDefect - defectExact];
expectNonzero["natural off compensation", naturalDefect];

tmNatRaw = cStage/gNat;
tmMinusRaw = cStage/gMinus;
tmPlusRaw = cStage/gPlus;

tractionRatioMinusRaw = reduceExact[tmMinusRaw/tmNatRaw];
tractionRatioPlusRaw = reduceExact[tmPlusRaw/tmNatRaw];
tractionResidualMinusBeforeAnsatz = reduceExact[tractionRatioMinusRaw - 1/gMinus];
tractionResidualPlusBeforeAnsatz = reduceExact[tractionRatioPlusRaw - 1/gPlus];

Print["T_m(-)/T_m(nat), before natural ansatz = ", fmt[tractionRatioMinusRaw]];
Print["T_m(+)/T_m(nat), before natural ansatz = ", fmt[tractionRatioPlusRaw]];
Print["traction residual (-), before natural ansatz = ", fmt[tractionResidualMinusBeforeAnsatz]];
Print["traction residual (+), before natural ansatz = ", fmt[tractionResidualPlusBeforeAnsatz]];

tractionRatioMinus = reduceExact[tractionRatioMinusRaw /. gNat -> gNatAnsatz];
tractionRatioPlus = reduceExact[tractionRatioPlusRaw /. gNat -> gNatAnsatz];

Print["numeric T ratio (-) = ", fmt[N[tractionRatioMinus, 20]]];
Print["numeric T ratio (+) = ", fmt[N[tractionRatioPlus, 20]]];

expectZero[
  "traction ratio (-) identity",
  tractionResidualMinusBeforeAnsatz /. gNat -> gNatAnsatz
];
expectZero[
  "traction ratio (+) identity",
  tractionResidualPlusBeforeAnsatz /. gNat -> gNatAnsatz
];

Print[""];
Print["Stage 122 Mathematica audit passed."];

Exit[0];
