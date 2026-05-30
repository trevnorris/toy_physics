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

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

banner["STAGE 144 — UNIQUE REGULAR CANONICAL MOUTH BRANCH"];

Clear[piM];
$Assumptions = Element[piM, Reals] && piM > 0;

r = Sqrt[4107 - 100*Pi^2]/(10*Pi);
gMinus = FullSimplify[r - Sqrt[1 + r^2]/2, Assumptions -> $Assumptions];
gPlus = FullSimplify[r + Sqrt[1 + r^2]/2, Assumptions -> $Assumptions];
gPi = 2*piM*(2*piM*Exp[piM] + Pi)/((4*piM^2 + Pi^2)*(Exp[piM] - 1));
sQ = piM*(((Pi/2)*Tanh[Pi/2]) + piM*(Exp[-piM]*Sech[Pi/2] - 1))/((1 - Exp[-piM])*((Pi/2)^2 - piM^2));
rQ = (gPi - r)^2/(1 + r^2);
sigma0 = piM/(1 - rQ*sQ);
tHat = Sqrt[(9/20)*sigma0];

subbanner["Compensated branches"];
Print["g_-^F1 = ", fmt[N[gMinus, 30]]];
Print["g_+^F1 = ", fmt[N[gPlus, 30]]];
Print["2/pi   = ", fmt[N[2/Pi, 30]]];

(* INDEPENDENT Pi_* route (NOT a transliteration of SymPy's nsolve on the *)
(* rational gPi == gMinus): clear the strictly-positive denominator *)
(* (4 piM^2 + Pi^2)(Exp[piM]-1) so the root equation is polynomial-in- *)
(* (piM, Exp[piM]), and isolate the unique positive root with a bracketing *)
(* seed pair. gPi is monotone on (0, Infinity), so the bracket robustly fixes *)
(* the root. Precedent: batch-4 stage 131-F3, batch-5 stage 142. *)
gThresholdResidual[p_] := 2*p*(2*p*Exp[p] + Pi) - gMinus*(4*p^2 + Pi^2)*(Exp[p] - 1);
piStar = piM /. FindRoot[gThresholdResidual[piM] == 0, {piM, 1.4, 1.6},
  WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
tHatStar = N[tHat /. piM -> piStar, 30];

piMatch = piM /. FindRoot[gPi == Pi/4, {piM, 1.9}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
tHatMatch = N[tHat /. piM -> piMatch, 30];

Print["Pi_*          = ", fmt[N[piStar, 30]]];
Print["That(Pi_*)    = ", fmt[tHatStar]];
Print["Pi_match      = ", fmt[N[piMatch, 30]]];
Print["That(Pi_match)= ", fmt[tHatMatch]];

If[!(N[piStar, 30] > 0 && N[piMatch, 30] > N[piStar, 30]), fail["unexpected ordering of canonical and matched-derivative points", {N[piStar, 20], N[piMatch, 20]}]];

sigma0Star = N[sigma0 /. piM -> piStar, 30];
sigma0Match = N[sigma0 /. piM -> piMatch, 30];
Print["Sigma0(Pi_*)  = ", fmt[sigma0Star]];
Print["Sigma0(Pi_match)= ", fmt[sigma0Match]];

If[!(N[gPlus, 30] > 1), fail["upper branch must satisfy g_+^F1 > 1", N[gPlus, 30]], pass["upper branch g_+^F1 > 1"]];
If[!(N[2/Pi, 30] < N[gMinus, 30] < 1), fail["lower branch must satisfy 2/pi < g_-^F1 < 1", N[gMinus, 30]], pass["lower branch bracket 2/pi < g_-^F1 < 1"]];
tol = 10^(-12);
(* (R1, INDEPENDENT root anchor) piStar above came from the cleared-denominator *)
(* FindRoot, NOT from gPi == gMinus; assert the residual is genuinely ~0 there, *)
(* then anchor to STAGE 131's owning Pi_* (its independent nsolve, prec 80): *)
(* scripts/output/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.txt:2 *)
piStarResidual = Chop[N[gThresholdResidual[piStar], 30], 10^-30];
If[!(piStarResidual === 0 || Abs[piStarResidual] < 10^-25),
  fail["Pi_* cleared-denominator residual not ~0", piStarResidual],
  pass["Pi_* solves cleared-denominator residual (independent root)"]];
piStarOwner = N[Rationalize[1.50882951349315558300555075595, 0], 30];
If[!(Abs[N[piStar - piStarOwner, 30]] < tol), fail["Pi_* drift vs stage 131 owner", N[piStar, 30]], pass["Pi_* matches stage-131 owning value (independent route)"]];
If[!(Abs[N[tHatStar - 0.901484054174205`30, 30]] < tol), fail["That(Pi_*) drift", tHatStar], pass["That(Pi_*) matches notes target"]];
If[!(Abs[N[sigma0Star - 1.80594111095636`30, 30]] < tol), fail["Sigma0(Pi_*) drift", sigma0Star], pass["Sigma0(Pi_*) matches notes target"]];
If[!(Abs[N[piMatch - 1.90848600654854`30, 30]] < tol), fail["Pi_match drift", N[piMatch, 30]], pass["Pi_match matches notes target"]];
If[!(Abs[N[tHatMatch - 1.01132972803599`30, 30]] < tol), fail["That(Pi_match) drift", tHatMatch], pass["That(Pi_match) matches notes target"]];

banner["STAGE 144 LEDGER"];
Print["Positive-source theorem gives 0 <= g <= 1, so the upper compensated branch is impossible."];
Print["The exponential mouth family is monotone and spans (2/pi, 1)."];
Print["Since g_-^F1 lies strictly inside that interval, it is reached at one unique finite Pi_*."];
Print["That point has finite moderate traction, whereas g=1 is only the singular Pi->infty limit."];
Print[""];
Print["Pi_*       = ", fmt[N[piStar, 30]]];
Print["That_*     = ", fmt[tHatStar]];
Print["Pi_match   = ", fmt[N[piMatch, 30]]];
Print["That_match = ", fmt[tHatMatch]];

Print[""];
Print["Stage 144 Mathematica audit passed."];

Exit[0];
