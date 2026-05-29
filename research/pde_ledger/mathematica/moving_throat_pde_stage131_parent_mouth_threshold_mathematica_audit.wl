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

banner["STAGE 131 — PARENT MICRO-THRESHOLD FOR CANONICAL MOUTH COMPENSATION"];

Clear[piM, thetaSigma, lM, tM, qStar, a0Prime];
$Assumptions =
  Element[{piM, thetaSigma, lM, tM, qStar, a0Prime}, Reals] &&
  piM > 0 && thetaSigma > 0 && lM > 0;

(* F1 lower-branch value g_-^{F1}: closed form below; see notes/stages/
   moving_throat_pde_stage131_parent_mouth_threshold.md Sec. 3 for context. *)
gMinusExact = (2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi);
gMinusLiteral = 0.758035078944663`50;
expectApprox["g_-^F1 closed form vs literal",
  N[gMinusExact, 50], gMinusLiteral, 10^-14];
gMinus = N[gMinusExact, 80];
gPi = FullSimplify[2*piM*(2*piM*Exp[piM] + Pi)/((4*piM^2 + Pi^2)*(Exp[piM] - 1)), Assumptions -> $Assumptions];
(* INDEPENDENT Pi_* route (not a transliteration of SymPy's nsolve on gPi == gMinus): *)
(* clear denominators so the root equation is polynomial-in-(piM, Exp[piM]) rather than *)
(* the rational gPi == gMinus, and isolate the unique positive root with a bracketing *)
(* seed pair. g_Pi is monotone on (0, Infinity), so the bracket robustly fixes the root. *)
gThresholdResidual[p_] := 40*Pi*p*(2*p*Exp[p] + Pi) - 20*Pi*gMinus*(4*p^2 + Pi^2)*(Exp[p] - 1);
piStar = piM /. FindRoot[gThresholdResidual[piM] == 0, {piM, 1.4, 1.6},
  WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];

v1 = piM*thetaSigma/lM;
v1Star = N[piStar, 30]*thetaSigma/lM;
gPrimeStar = N[D[gPi, piM] /. piM -> piStar, 30];
thresholdResidual = FullSimplify[(tM - qStar*a0Prime) - piM*thetaSigma/lM, Assumptions -> $Assumptions];

Print["Pi_* = ", fmt[N[piStar, 30]]];
Print["V1_* = ", fmt[v1Star]];
Print["g'(Pi_*) = ", fmt[gPrimeStar]];
Print["Parent bias mismatch formula = ", fmt[thresholdResidual]];
(* Anchored checks vs notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md *)

(* Anchor 1: piStar matches notes Sec. 1 quoted value. *)
expectApprox["piStar notes Sec. 1 value",
  N[piStar, 50], 1.50882951349316`50, 10^-14];

(* Anchor 2: slope D[gPi,piM] at piStar matches notes Sec. 3 quoted value. *)
expectApprox["slope at piStar notes Sec. 3 value",
  N[D[gPi, piM] /. piM -> piStar, 50], 0.0714453558083195`50, 10^-14];

(* Anchor 4: lower-vs-singular branch discrimination, paper Checks item 3. *)
(* gPi rises from 2/Pi to a supremum of 1, so gNat = 1 is the unreachable singular *)
(* equal-normalized branch and gPlus > 1 is never attained. Wrap each Rule in parens. *)
gNat = 1;
gPlusExact = (2*Sqrt[4107 - 100*Pi^2] + 37*Sqrt[3])/(20*Pi);
gPiAtStar = N[(gPi /. piM -> piStar), 40];

(* 4a: lower-branch MEMBERSHIP. *)
lowerResidual = Abs[N[gPiAtStar - gMinus, 40]];
If[TrueQ[lowerResidual < 10^-30],
  pass["Pi_* on lower branch (membership)"],
  fail["Pi_* on lower branch (membership)", lowerResidual]
];

(* 4b: SINGULAR equal-normalized branch EXCLUDED, separation matches notes Delta g_-. *)
singSep = N[gNat - gPiAtStar, 30];
deltaGMinusNotes = 0.241964921055337`30;
If[TrueQ[Abs[N[singSep - deltaGMinusNotes, 30]] < 10^-12 && singSep > 10^-3],
  pass["singular equal-normalized branch excluded (notes Delta g_-)"],
  fail["singular equal-normalized branch excluded (notes Delta g_-)", singSep]
];

(* 4c: UPPER branch EXCLUDED. *)
upperSep = Abs[N[gPiAtStar - N[gPlusExact, 40], 30]];
If[TrueQ[upperSep > 1],
  pass["upper branch excluded"],
  fail["upper branch excluded", upperSep]
];

Print[""];
Print["Stage 131 Mathematica audit passed."];

Exit[0];
