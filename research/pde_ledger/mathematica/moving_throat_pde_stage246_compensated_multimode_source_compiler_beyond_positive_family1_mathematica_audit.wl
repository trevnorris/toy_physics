ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
);

fmt[expr_] := ToString[InputForm[expr]];
pass[name_String] := Print["PASS: ", name];

fail[name_String, detail_] := (
  Print["FAIL: ", name];
  Print["  residual -> ", fmt[detail]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  res = res /. ConditionalExpression[value_, _] :> value;
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, statement_] := Module[{res},
  res = FullSimplify[statement, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

expectApprox[name_String, value_, target_, tol_] := Module[{delta},
  delta = N[Abs[value - target], 50];
  Print[name, " delta = ", fmt[delta]];
  If[TrueQ[delta <= tol], pass[name], fail[name, delta]];
];

banner["STAGE 246 - COMPENSATED MULTIMODE SOURCE COMPILER BEYOND POSITIVE FAMILY-1"];

Clear[a, b, x, y, gt, St, gc, r, rsig, a0, b0];
$Assumptions = Element[{a, b, x, y, gt, St, gc, r, rsig, a0, b0}, Reals];

source = 1 + a Cos[Pi x] + b Cos[2 Pi x];

subbanner["M1-M3. Direct moment integrals"];

meanMoment = Integrate[source, {x, 0, 1}, GenerateConditions -> False];
mouthMoment = Integrate[
  source Cos[Pi x/2],
  {x, 0, 1},
  GenerateConditions -> False
];
shellKernel = Cosh[Pi (1 - x)/2]/Cosh[Pi/2];
shellMoment = Integrate[
  source shellKernel,
  {x, 0, 1},
  GenerateConditions -> False
];

mouthExpected = (2/Pi) (1 + a/3 - b/15);
shellExpected = (2 Tanh[Pi/2]/Pi) (1 + a/5 + b/17);

expectZero["M1 mean preservation", meanMoment - 1];
expectTrue[
  "M2 mouth-bias functional",
  FullSimplify[mouthMoment - mouthExpected == 0, Assumptions -> $Assumptions]
];
expectTrue[
  "M3 shell-loading functional",
  FullSimplify[shellMoment - shellExpected == 0, Assumptions -> $Assumptions]
];

subbanner["M4. Exact minimum via MinValue"];

m4Cases = {
  {{a -> 1/2, b -> -1/5}, 3/10, "M4 minimum at a=1/2,b=-1/5"},
  {{a -> 5/2, b -> 1/4}, -5/4, "M4 minimum at a=5/2,b=1/4"},
  {{a -> 1/2, b -> 1/2}, 7/16, "M4 minimum at a=1/2,b=1/2"}
};

Do[
  Module[{rules, target, name, minValue},
    {rules, target, name} = entry;
    minValue = FullSimplify[MinValue[{source /. rules, 0 <= x <= 1}, x]];
    expectApprox[name, minValue, target, 10^-20];
  ],
  {entry, m4Cases}
];

subbanner["M5. Two-moment map"];

Msrc = {{1/3, -1/15}, {1/5, 1/17}};
inverseImage = FullSimplify[
  Inverse[Msrc] . {gt, St},
  Assumptions -> $Assumptions
];
inverseExpected = {
  (85/42) St + (25/14) gt,
  (425/42) St - (85/14) gt
};

expectZero["M5 det(Msrc)", Det[Msrc] - 14/425];
expectTrue[
  "M5 inverse coefficients",
  FullSimplify[inverseImage == inverseExpected, Assumptions -> $Assumptions]
];

subbanner["M6-M7. Quarter-ratio and compensation line"];

rF1 = Sqrt[(12/Pi^2) (37/20)^2 - 1];
gMinus = rF1 - (1/2) Sqrt[1 + rF1^2];
gPlus = rF1 + (1/2) Sqrt[1 + rF1^2];

expectZero["M6 quarter-ratio lower branch", (gMinus - rF1)^2/(1 + rF1^2) - 1/4];
expectZero["M6 quarter-ratio upper branch", (gPlus - rF1)^2/(1 + rF1^2) - 1/4];

bSolved = b /. First[Solve[(2/Pi) (1 + a/3 - b/15) == gc, b]];
expectZero["M7 compensation line", bSolved - (5 a + 15 - (15 Pi/2) gc)];

subbanner["M8. Transported sign-change threshold"];

sRad = rsig^2/(r^2 + rsig^2);
aRad = a0 sRad;
bRad = b0 sRad;
orientationAssumptions = a0 > 0 && b0 < 0 && rsig > 0 && r > 0;
transportMin = FullSimplify[
  1 + bRad - Abs[aRad],
  Assumptions -> orientationAssumptions
];
thresholdReduced = FullSimplify[
  Reduce[
    orientationAssumptions && a0 - b0 > 1 && transportMin < 0,
    r,
    Reals
  ],
  Assumptions -> a0 > 0 && b0 < 0 && rsig > 0 && a0 - b0 > 1
];

expectZero["M8 transported boundary minimum", transportMin - (1 - (a0 - b0) sRad)];
expectTrue[
  "M8 threshold solves to r bound",
  FullSimplify[
    thresholdReduced \[Equivalent] (a0 > 0 && b0 < 0 && rsig > 0 &&
      a0 - b0 > 1 && 0 < r < rsig Sqrt[a0 - b0 - 1]),
    Assumptions -> a0 > 0 && b0 < 0 && rsig > 0 && a0 - b0 > 1
  ]
];

subbanner["M9. Session-I numeric readback"];

a0Num = 11/5;
b0Num = -3/5;
rsigNum = 4/5;
reval = 25054257/25000000;
sNum[rv_] := rsigNum^2/(rv^2 + rsigNum^2);
gNum[rv_] := (2/Pi) (1 + a0Num sNum[rv]/3 - b0Num sNum[rv]/15);
SNum[rv_] := (2 Tanh[Pi/2]/Pi) (1 + a0Num sNum[rv]/5 + b0Num sNum[rv]/17);
RNum[rv_] := (gNum[rv] - rF1)^2/(1 + rF1^2);
sigmaMinNum[rv_] := 1 - (a0Num - b0Num) sNum[rv];
rThrNum = rsigNum Sqrt[a0Num - b0Num - 1];

expectApprox["M9 g(r_eval)", gNum[reval], 0.82823667, 5*10^-9];
expectApprox["M9 S(r_eval)", SNum[reval], 0.67584771, 5*10^-9];
expectApprox["M9 R(r_eval)", RNum[reval], 0.21677037, 5*10^-9];
expectApprox["M9 sigmaMin(r_eval)", sigmaMinNum[reval], -0.08979545, 5*10^-9];
expectTrue["M9 g(0) > 1", gNum[0] > 1];
expectTrue["M9 rThr > r_eval", rThrNum > reval];

banner["Stage 246 Mathematica audit result"];
Print["All Mathematica checks passed."];
