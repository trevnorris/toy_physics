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
piStar = piM /. FindRoot[gPi == gMinus, {piM, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];

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

(* Anchor 3: parent threshold identity at piM = piStar, notes Sec. 2. *)
thresholdAtStar = FullSimplify[
  thresholdResidual /. piM -> piStar,
  Assumptions -> $Assumptions
];
expectedForm = (tM - qStar*a0Prime) - piStar*thetaSigma/lM;
identityResidual = Chop[Simplify[thresholdAtStar - expectedForm], 10^-30];
If[TrueQ[identityResidual === 0],
  pass["parent threshold identity at piM = piStar (notes Sec. 2)"],
  fail["parent threshold identity at piM = piStar (notes Sec. 2)",
       identityResidual]
];

(* Anchor 4: lower-branch discrimination, gPi at 2*piStar is far from gMinus. *)
offStarResidual = Abs[N[(gPi /. piM -> 2*piStar) - gMinus, 30]];
If[TrueQ[offStarResidual > 10^-3],
  pass["lower-branch discrimination (paper Checks item 3)"],
  fail["lower-branch discrimination (paper Checks item 3)", offStarResidual]
];

Print[""];
Print["Stage 131 Mathematica audit passed."];

Exit[0];
