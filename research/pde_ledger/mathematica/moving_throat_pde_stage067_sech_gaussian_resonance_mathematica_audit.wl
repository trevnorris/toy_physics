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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 50] - N[target, 50]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = TrueQ[cond];
  Print[name, " = ", fmt[res]];
  If[res, pass[name], fail[name, cond]];
];

banner["STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK"];

Clear[y, x, wf, wg, r];
$Assumptions =
  Element[{y, x, wf, wg, r}, Reals] &&
  wf > 0 && wg > 0 && r > 0;

nssDirect = FullSimplify[Integrate[Sech[y/wf]^2, {y, -Infinity, Infinity}], Assumptions -> $Assumptions];
nppDirect = FullSimplify[Integrate[Exp[-2*y^2/wg^2], {y, -Infinity, Infinity}], Assumptions -> $Assumptions];
nssExpected = 2*wf;
nppExpected = wg*Sqrt[Pi/2];

Print["N_(sigma sigma) = ", fmt[nssDirect]];
Print["N_(phi phi) = ", fmt[nppDirect]];
expectZero["N_(sigma sigma) - 2 w_f", nssDirect - nssExpected];
expectZero["N_(phi phi) - w_g sqrt(pi/2)", nppDirect - nppExpected];

overlapHead = Symbol["OverlapI"];
c2Sym = FullSimplify[overlapHead[r]^2/(r*Sqrt[2*Pi]), Assumptions -> $Assumptions];
Print["C^2(r) = ", fmt[c2Sym]];

banner["EXACT DUALITY IMPLICATION"];

dualityRhs = FullSimplify[(r/Sqrt[Pi])*overlapHead[Pi/r], Assumptions -> $Assumptions];
c2Dual = FullSimplify[dualityRhs^2/(r*Sqrt[2*Pi]), Assumptions -> $Assumptions];
c2Target = FullSimplify[overlapHead[Pi/r]^2/((Pi/r)*Sqrt[2*Pi]), Assumptions -> $Assumptions];
(* Algebraic implication only: substitutes OverlapI -> (r/Sqrt[Pi]) OverlapI[Pi/r] into
   C^2(r) and checks it equals C^2(Pi/r). Holds for ANY function OverlapI; the duality
   identity for the sech-Gaussian overlap is exercised numerically below. *)
expectZero["C^2(r) - C^2(pi/r) under duality", c2Dual - c2Target];

banner["SELF-DUAL STATIONARY POINT"];

rStar = Sqrt[Pi];
c2Fn = Symbol["C2fn"];
c2PrimeLeft = Symbol["C2PrimeLeft"];
c2PrimeDual = Symbol["C2PrimeDual"];
c2SymmetryTangent = FullSimplify[
  D[c2Fn[r] - c2Fn[Pi/r], r] /. {
    Derivative[1][c2Fn][r] -> c2PrimeLeft,
    Derivative[1][c2Fn][Pi/r] -> c2PrimeDual
  },
  Assumptions -> $Assumptions
];
c2SymmetryAtRStar = FullSimplify[
  c2SymmetryTangent /. {r -> rStar, c2PrimeDual -> c2PrimeLeft},
  Assumptions -> $Assumptions
];
Print["differentiated C^2 symmetry at r_* = ", fmt[c2SymmetryAtRStar]];
c2StationarySolve = Solve[c2SymmetryAtRStar == 0, c2PrimeLeft, Reals];
Print["Solve[c2 symmetry tangent == 0] = ", fmt[c2StationarySolve]];
(* Tautological: Solve[2*C2PrimeLeft == 0] returns C2PrimeLeft -> 0; substituting that
   back into C2PrimeLeft yields 0. This is the calculus fact that a function symmetric
   under r <-> Pi/r has zero derivative at r = Sqrt[Pi], not a sech-Gaussian-specific
   result. The numerical monotonicity scan below provides the substantive evidence. *)
expectZero[
  "self-dual C^2 stationary slope from symmetry solve",
  FullSimplify[c2PrimeLeft /. First[c2StationarySolve], Assumptions -> $Assumptions]
];

banner["NUMERICAL OVERLAP BENCHMARK"];

Clear[overlapNum, c2Num];
overlapNum[rr_?NumericQ] := overlapNum[rr] = Module[{rrp, wp = 80, ag = 32, pg = 32},
  rrp = SetPrecision[rr, wp];
  Quiet[
    2*NIntegrate[
      Evaluate@N[Sech[x]*Exp[-x^2/rrp^2], wp],
      {x, 0, Infinity},
      WorkingPrecision -> wp,
      AccuracyGoal -> ag,
      PrecisionGoal -> pg,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
    ],
    NIntegrate::precw
  ]
];

c2Num[rr_?NumericQ] := c2Num[rr] = overlapNum[rr]^2/(rr*Sqrt[2*Pi]);

rStarNum = N[rStar, 60];
c2Star = N[c2Num[rStarNum], 60];
pres = N[1/c2Star, 60];

(* c2Target / presTarget are the sympy mpmath quad results from
   scripts/output/moving_throat_pde_stage067_sech_gaussian_sympy_audit.txt.
   This block confirms cross-engine numerical agreement on the same definite
   integral, not agreement with any closed-form benchmark. *)
c2Target = ToExpression["0.994418836451529348706428351608877628170873348983716948813464`60"];
presTarget = ToExpression["1.00561248776057621695172301479763550405448504648609605997534`60"];

Print["r_* = ", fmt[rStarNum]];
Print["C_res^2 = ", fmt[c2Star]];
Print["P_res = ", fmt[pres]];
Print["1 - C_res^2 = ", fmt[N[1 - c2Star, 50]]];

expectApprox["C_res^2 numeric check", c2Star, c2Target, 10^-35];
expectApprox["P_res numeric check", pres, presTarget, 10^-34];

banner["NUMERICAL DUALITY CHECKS"];

sampleGrid = {3/4, 1, 6/5, 3/2, 2};
Do[
  lhs = N[overlapNum[N[rr, 50]], 50];
  rhs = N[(rr/Sqrt[Pi])*overlapNum[N[Pi/rr, 50]], 50];
  diff = Abs[lhs - rhs];
  Print["r = ", fmt[N[rr, 20]], ": |I(r) - r/sqrt(pi) I(pi/r)| = ", fmt[N[diff, 20]]];
  expectTrue["duality sample " <> fmt[N[rr, 20]], diff <= 10^-35],
  {rr, sampleGrid}
];

banner["NUMERICAL MONOTONICITY SCAN"];

leftGrid = {11/20, 7/10, 9/10, 11/10, 27/20, 8/5, rStar};
leftVals = N[c2Num[N[#, 50]], 40] & /@ leftGrid;
Do[
  Print["C^2(", fmt[N[leftGrid[[k]], 20]], ") = ", fmt[leftVals[[k]]]],
  {k, Length[leftGrid]}
];
expectTrue[
  "constructive-branch increase up to r_*",
  And @@ Table[leftVals[[k + 1]] > leftVals[[k]], {k, Length[leftVals] - 1}]
];

rightGrid = {rStar, 2, 5/2, 3, 4};
rightVals = N[c2Num[N[#, 50]], 40] & /@ rightGrid;
Do[
  Print["C^2(", fmt[N[rightGrid[[k]], 20]], ") = ", fmt[rightVals[[k]]]],
  {k, Length[rightGrid]}
];
expectTrue[
  "constructive-branch decrease after r_*",
  And @@ Table[rightVals[[k + 1]] < rightVals[[k]], {k, Length[rightVals] - 1}]
];

Print[""];
Print["Stage 067 Mathematica audit passed."];

Exit[0];
