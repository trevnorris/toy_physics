ClearAll["Global`*"];
$HistoryLength = 0;
$MaxExtraPrecision = 1000;

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  detail -> ", fmt[detail]]];
  Exit[1];
);

stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

cleanScalar[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanScalar[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

expectClose[name_String, actual_, expected_, tol_] := Module[{diff},
  diff = N[Abs[actual - expected], 40];
  Print[name, " actual = ", fmt[N[actual, 18]], " expected = ", fmt[N[expected, 18]], " diff = ", fmt[diff]];
  If[TrueQ[diff < tol], pass[name], fail[name, diff]];
];

expectInfinity[name_String, actual_] := (
  Print[name, " = ", fmt[actual]];
  If[TrueQ[actual === Infinity], pass[name], fail[name, actual]];
);

Print["=== Stage 230 Mathematica audit: selected-branch classifier -> dynamic window ==="];

tolTight = 2 10^-15;
tolSample = 5 10^-13;
$Assumptions = Element[{xi, delta, R}, Reals] && delta > 0 && R >= 0;

(* M1: exact classifier and monotonicity. *)
rND = 72 delta^2 (1 - xi)/((9 delta + 11 xi) (9 delta^2 + 18 delta xi + 11 xi^2));
rNDAtOnset = rND /. xi -> 0;
rNDSoft = Block[
  {$Assumptions = Element[delta, Reals] && delta > 0},
  Limit[rND, xi -> 1, Direction -> "FromBelow"]
];
rNDMonotoneReduce = Reduce[0 <= xi < 1 && delta > 0 && D[rND, xi] < 0, {xi, delta}, Reals];

Print["M1 Reduce domain = ", fmt[rNDMonotoneReduce]];
expectZero["M1 R_ND(0,delta) - 8/(9 delta)", rNDAtOnset - 8/(9 delta)];
expectZero["M1 limit xi -> 1-", rNDSoft];
expectTrue[
  "M1 D[R_ND,xi] < 0 on stable strip by Reduce",
  Resolve[ForAll[{xi, delta}, (0 <= xi < 1 && delta > 0) \[Implies] rNDMonotoneReduce], Reals]
];
expectTrue[
  "M1 Reduce recovers stable strip",
  Resolve[ForAll[{xi, delta}, rNDMonotoneReduce \[Equivalent] (0 <= xi < 1 && delta > 0)], Reals]
];

(* M2: affine share compiler. *)
sPlusNum = -508643465308977/10^15;
sPlusDen = -301516097158113/10^15;
sMinusNum = -334368725711457/10^15;
sMinusDen = 411024574532864/10^15;

wNum = R/(1 + R);
wDen = 1/(1 + R);
sPlus = FullSimplify[wNum sPlusNum + wDen sPlusDen];
sMinus = FullSimplify[wNum sMinusNum + wDen sMinusDen];
sPlusReduce = Reduce[D[sPlus, R] < 0 && R >= 0, R, Reals];
sMinusReduce = Reduce[D[sMinus, R] < 0 && R >= 0, R, Reals];

Print["M2 S_plus(R) = ", fmt[sPlus]];
Print["M2 S_minus(R) = ", fmt[sMinus]];
Print["M2 Reduce S_plus' = ", fmt[sPlusReduce]];
Print["M2 Reduce S_minus' = ", fmt[sMinusReduce]];
expectZero["M2 S_plus affine share form", sPlus - (R sPlusNum + sPlusDen)/(1 + R)];
expectZero["M2 S_minus affine share form", sMinus - (R sMinusNum + sMinusDen)/(1 + R)];
expectTrue[
  "M2 D[S_plus,R] < 0 for R >= 0 by Reduce",
  Resolve[ForAll[R, sPlusReduce \[Equivalent] R >= 0], Reals]
];
expectTrue[
  "M2 D[S_minus,R] < 0 for R >= 0 by Reduce",
  Resolve[ForAll[R, sMinusReduce \[Equivalent] R >= 0], Reals]
];

(* M3: sign-flip threshold. *)
rStar = FullSimplify[sMinusDen/(-sMinusNum)];
rRoots = R /. Solve[sMinus == 0, R, Reals];
rRootsNonnegative = Select[rRoots, TrueQ[FullSimplify[# >= 0]] &];
sMinusPositiveReduce = Reduce[0 <= R < rStar && sMinus > 0, R, Reals];
sMinusNegativeReduce = Reduce[R > rStar && sMinus < 0, R, Reals];

Print["M3 R roots = ", fmt[rRoots]];
Print["M3 Reduce S_minus positive side = ", fmt[sMinusPositiveReduce]];
Print["M3 Reduce S_minus negative side = ", fmt[sMinusNegativeReduce]];
expectZero["M3 unique nonnegative root count", Length[rRootsNonnegative] - 1];
expectZero["M3 solved root equals R_star", First[rRootsNonnegative] - rStar];
expectClose["M3 N[R_star,30]", N[rStar, 30], 1229255438463336/10^15, 2 10^-15];
expectTrue["M3 S_plus(R_star) < 0", (sPlus /. R -> rStar) < 0];
expectTrue[
  "M3 S_minus positive for 0 <= R < R_star",
  Resolve[ForAll[R, (0 <= R < rStar) \[Implies] sMinus > 0], Reals]
];
expectTrue[
  "M3 S_minus negative for R > R_star",
  Resolve[ForAll[R, R > rStar \[Implies] sMinus < 0], Reals]
];

(* M4: onset threshold by Solve. *)
deltaDynStar = FullSimplify[8/(9 rStar)];
deltaRoots = delta /. Solve[rNDAtOnset == rStar, delta, Reals];
deltaRootsPositive = Select[deltaRoots, TrueQ[FullSimplify[# > 0]] &];

Print["M4 delta roots = ", fmt[deltaRoots]];
expectZero["M4 unique positive onset root count", Length[deltaRootsPositive] - 1];
expectZero["M4 solved onset root equals delta_dyn_star", First[deltaRootsPositive] - deltaDynStar];
expectClose["M4 N[delta_dyn_star,30]", N[deltaDynStar, 30], 723111617875019/10^15, 2 10^-15];

(* M5: dynamic ceilings in |eps Xi_1|. *)
rqMinus = 30199907560250075/10^15;
rqPlus = 36171186483269487/10^15;
rqReq = 21854566296358396/10^15;
ellMinus = Log[rqMinus/rqReq];
ellPlus = Log[rqPlus/rqReq];

bPlus = FullSimplify[ellPlus/(-sPlus)];
bMinus = FullSimplify[ellMinus/(-sMinus)];

ceilingBothValue[rval_] := Module[{sp, sm, bp, bm},
  sp = N[sPlus /. R -> rval, 50];
  sm = N[sMinus /. R -> rval, 50];
  bp = N[ellPlus/(-sp), 50];
  Which[
    TrueQ[sm < 0],
      bm = N[ellMinus/(-sm), 50];
      Min[bp, bm],
    True,
      bp
  ]
];

ceilingNonemptyValue[rval_] := Module[{sp, sm, bp, bm},
  sp = N[sPlus /. R -> rval, 50];
  sm = N[sMinus /. R -> rval, 50];
  Which[
    TrueQ[sm < 0],
      bp = N[ellPlus/(-sp), 50];
      bm = N[ellMinus/(-sm), 50];
      Max[bp, bm],
    True,
      Infinity
  ]
];

bPlusInf = Block[{$Assumptions = True}, Limit[bPlus, R -> Infinity]];
bMinusInf = Block[{$Assumptions = True}, Limit[bMinus, R -> Infinity]];
expectTrue["M5 endpoint ordering B_minus(inf) < B_plus(inf)", N[bMinusInf, 50] < N[bPlusInf, 50]];
bBothInf = bMinusInf;
bNonemptyInf = bPlusInf;

expectClose["M5 ell_minus", N[ellMinus, 30], 323428979934714/10^15, tolTight];
expectClose["M5 ell_plus", N[ellPlus, 30], 503852964869151/10^15, tolTight];
expectClose["M5 B_both(0)", ceilingBothValue[0], 1671064893775584/10^15, tolTight];
expectInfinity["M5 B_nonempty(0)", ceilingNonemptyValue[0]];
expectClose["M5 Limit B_both at infinity", N[bBothInf, 30], 967282389363822/10^15, tolTight];
expectClose["M5 Limit B_nonempty at infinity", N[bNonemptyInf, 30], 990581810705233/10^15, tolTight];

sampleRows = {
  {0, -301516097158113/10^15, 411024574532864/10^15, 1671064893775584/10^15, Infinity},
  {1, -405079781233545/10^15, 383279244107035/10^16, 1243836370541187/10^15, Infinity},
  {rStar, -415730215182002/10^15, 0, 1211971000588856/10^15, Infinity},
  {10, -489813704567990/10^15, -266605698416519/10^15, 1028662448947899/10^15, 1213136035184892/10^15}
};

Do[
  Module[{rval, expectedSp, expectedSm, expectedBoth, expectedNonempty, spVal, smVal, bothVal, nonemptyVal},
    {rval, expectedSp, expectedSm, expectedBoth, expectedNonempty} = row;
    spVal = N[sPlus /. R -> rval, 50];
    smVal = N[sMinus /. R -> rval, 50];
    bothVal = ceilingBothValue[rval];
    nonemptyVal = ceilingNonemptyValue[rval];
    Print[
      "M5 sample R = ", fmt[N[rval, 18]],
      " S_plus = ", fmt[N[spVal, 18]],
      " S_minus = ", fmt[N[smVal, 18]],
      " B_both = ", fmt[N[bothVal, 18]],
      " B_nonempty = ", If[expectedNonempty === Infinity, "+inf", fmt[N[nonemptyVal, 18]]]
    ];
    expectClose["M5 sample S_plus", spVal, expectedSp, tolSample];
    expectClose["M5 sample S_minus", smVal, expectedSm, tolSample];
    expectClose["M5 sample B_both", bothVal, expectedBoth, tolSample];
    If[expectedNonempty === Infinity,
      expectInfinity["M5 sample B_nonempty", nonemptyVal],
      expectClose["M5 sample B_nonempty", nonemptyVal, expectedNonempty, tolSample]
    ];
  ],
  {row, sampleRows}
];

(* M6: static-first theorem. *)
bPlusMonotoneReduce = Reduce[D[bPlus, R] < 0 && R >= 0, R, Reals];
bMinusMonotoneReduce = Reduce[D[bMinus, R] < 0 && R > rStar, R, Reals];
bStatBoth = 367930328492646/10^15;
bStatNonempty = 737619063660757/10^15;

Print["M6 Reduce B_plus' = ", fmt[bPlusMonotoneReduce]];
Print["M6 Reduce B_minus' = ", fmt[bMinusMonotoneReduce]];
expectTrue[
  "M6 B_plus decreases on R >= 0 by Reduce",
  Resolve[ForAll[R, bPlusMonotoneReduce \[Equivalent] R >= 0], Reals]
];
expectTrue[
  "M6 B_minus decreases on finite nonempty region by Reduce",
  Resolve[ForAll[R, bMinusMonotoneReduce \[Equivalent] R > rStar], Reals]
];
expectTrue["M6 inf B_both exceeds static both budget", N[bBothInf, 50] > N[bStatBoth, 50]];
expectTrue["M6 inf B_nonempty exceeds static nonempty budget", N[bNonemptyInf, 50] > N[bStatNonempty, 50]];

Print["R_star = ", fmt[N[rStar, 30]]];
Print["delta_dyn_star = ", fmt[N[deltaDynStar, 30]]];
Print["inf B_dyn^(both) = ", fmt[N[bBothInf, 30]]];
Print["inf B_dyn^(nonempty) = ", fmt[N[bNonemptyInf, 30]]];
Print["All Stage 230 Mathematica audits passed."];
Exit[0];
