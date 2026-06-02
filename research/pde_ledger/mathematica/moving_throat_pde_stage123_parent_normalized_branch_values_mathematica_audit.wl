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
stripCE[expr_] := expr /. ConditionalExpression[e_, _] :> e;

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

reduceExact[expr_] := Module[{res},
  res = stripCE[expr];
  res = FullSimplify[PowerExpand[Cancel[Together[res]]], Assumptions -> $Assumptions];
  res = FullSimplify[RootReduce[res], Assumptions -> $Assumptions];
  stripCE[res]
];

expectZero[name_String, expr_] := Module[{res},
  res = reduceExact[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

branchesFromReduce[red_, var_] := DeleteDuplicates[
  Join[
    Cases[LogicalExpand[red], HoldPattern[var == rhs_] :> rhs, Infinity],
    Cases[LogicalExpand[red], HoldPattern[lhs_ == var] :> lhs, Infinity]
  ]
];

uniqueBranch[name_String, red_, var_] := Module[{branches},
  branches = reduceExact /@ branchesFromReduce[red, var];
  If[Length[branches] =!= 1, fail[name <> " branch selection", red]];
  First[branches]
];

printNumeric[name_String, expr_] := Print[name, " numeric = ", fmt[N[RootReduce[expr], 20]]];

banner["STAGE 123 - PARENT-NORMALIZED BRANCH VALUES"];

Clear[a, L, ell, mu0, Zq, m, rho, q, hbar, cs, Tm, vw, rr, gg];

$Assumptions =
  Element[{a, L, ell, mu0, Zq, m, rho, q, hbar, cs, Tm, rr, gg}, Reals] &&
  a > 0 && L > 0 && ell > 0 && mu0 > 0 && Zq > 0 &&
  m > 0 && rho > 0 && q > 0 && hbar > 0 && cs > 0 &&
  Tm > 0 && rr > 0 && gg > 0 &&
  Element[vw, Reals];

Print["Stage 118 anchors: K_s, K_q, J_s, negative lambda, and the Xi_T healing lock."];
Print["Stage 121 and 122 anchors: r_F1 and g branches."];

Ks = 3*Pi*a^2*hbar^2/(5*m*rho*ell);
Kq = (Zq/mu0)*Pi^2*cs^2/(4*L^2);
Js = 4*Pi*a^2*ell/3;
lambda = -(8*Sqrt[2]/3)*q*vw*a^2*ell*Sqrt[L];

rParent = reduceExact[lambda/Sqrt[Ks*Kq]];
xiVDef = reduceExact[
  q*Sqrt[mu0*m*rho]*a*L^(3/2)*ell^(3/2)*vw/(hbar*Sqrt[Zq]*cs)
];

vwReduce = Reduce[rr == rParent && Element[vw, Reals], vw, Reals];
vwFromR = uniqueBranch["v_w0 from r", vwReduce, vw];
xiVFromParent = reduceExact[xiVDef /. vw -> vwFromR];

Print["v_w0(r) = ", fmt[vwFromR]];
Print["Xi_v(r) derived = ", fmt[xiVFromParent]];

xiVExpected = -(3*Sqrt[30]*Pi^(3/2)/160)*rr;
expectZero["Xi_v law", xiVFromParent - xiVExpected];

gParentUnlocked = Sqrt[2*Zq*Ks]/(Tm*Js*cs*Sqrt[mu0*L]);
gParentLocked = reduceExact[gParentUnlocked /. cs -> hbar/(2*m*ell)];
xiTDef = reduceExact[Tm*a*Sqrt[mu0*rho*L*ell]/Sqrt[Zq*m]];

tmReduce = Reduce[gg == gParentLocked && Tm > 0 && gg > 0, Tm, Reals];
tmFromG = uniqueBranch["T_m from g", tmReduce, Tm];
xiTFromParent = reduceExact[xiTDef /. Tm -> tmFromG];

Print["T_m(g) = ", fmt[tmFromG]];
Print["Xi_T(g) derived = ", fmt[xiTFromParent]];

xiTExpected = (3*Sqrt[30]/(10*Sqrt[Pi]))/gg;
expectZero["Xi_T law", xiTFromParent - xiTExpected];

R = 37/20;
rF1 = reduceExact[Sqrt[12*R^2/Pi^2 - 1]];
gMinus = reduceExact[(2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi)];
gPlus = reduceExact[(2*Sqrt[4107 - 100*Pi^2] + 37*Sqrt[3])/(20*Pi)];

xiVF1 = reduceExact[xiVFromParent /. rr -> rF1];
xiVF1Expected = reduceExact[xiVExpected /. rr -> rF1];
printNumeric["Xi_v(F1)", xiVF1];
expectZero["Xi_v(F1)", xiVF1 - xiVF1Expected];

xiTNat = reduceExact[xiTFromParent /. gg -> 1];
xiTNatExpected = reduceExact[xiTExpected /. gg -> 1];
printNumeric["Xi_T(nat)", xiTNat];
expectZero["Xi_T(nat)", xiTNat - xiTNatExpected];

xiTMinus = reduceExact[xiTFromParent /. gg -> gMinus];
xiTMinusExpected = reduceExact[xiTExpected /. gg -> gMinus];
printNumeric["Xi_T(-)", xiTMinus];
expectZero["Xi_T(-)", xiTMinus - xiTMinusExpected];

xiTPlus = reduceExact[xiTFromParent /. gg -> gPlus];
xiTPlusExpected = reduceExact[xiTExpected /. gg -> gPlus];
printNumeric["Xi_T(+)", xiTPlus];
expectZero["Xi_T(+)", xiTPlus - xiTPlusExpected];

Print[""];
Print["Stage 123 Mathematica audit passed."];

Exit[0];
