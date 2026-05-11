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

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 023 — GENERALIZED SELECTED-BRANCH NORMALIZATION"];

Clear[a0, delta, z0, q, eta, xi, rU, eps];
$Assumptions =
  Element[{a0, delta, z0, q, eta, xi, rU, eps}, Reals] &&
  a0 > 0 && delta > 0 && z0 > 0 && xi > 0 && rU > 0 && eps > 0 && q != 0;

lambda0 = 2/9;

subbanner["1. Exact 2x2 selected-branch solve with generic loading ratio"];

lamMinus = a0 (1 - xi);
alphaReq = FullSimplify[a0 xi (delta + xi)/(z0^2 (delta + xi) + (q z0)^2 xi), Assumptions -> $Assumptions];
r = FullSimplify[(a0 xi - alphaReq z0^2)/(alphaReq z0 (q z0)), Assumptions -> $Assumptions];

Print["alpha_req = ", fmt[alphaReq]];
Print["e1/e0 = ", fmt[r]];
expectZero["e1/e0 closed form", r - xi q/(delta + xi)];

subbanner["2. Exact overlap formulas and generalized F,G functions"];

zOverlapSq = FullSimplify[(1 + q r)^2/(1 + r^2), Assumptions -> $Assumptions];
sOverlapSq = FullSimplify[(1 + eta xi/(delta + xi))^2/(1 + r^2), Assumptions -> $Assumptions];
fGeneral = FullSimplify[(a0/lamMinus) zOverlapSq sOverlapSq, Assumptions -> $Assumptions];
gGeneral = FullSimplify[(z0^2/a0) alphaReq, Assumptions -> $Assumptions];

fExpected = FullSimplify[
  (delta + (1 + q^2) xi)^2 (delta + (1 + eta) xi)^2/
    ((1 - xi) ((delta + xi)^2 + q^2 xi^2)^2),
  Assumptions -> $Assumptions
];
gExpected = FullSimplify[xi (delta + xi)/(delta + (1 + q^2) xi), Assumptions -> $Assumptions];

Print["(z.e_-)^2 / z0^2 = ", fmt[zOverlapSq]];
Print["(s.e_-)^2 / s0^2 = ", fmt[sOverlapSq]];
Print["F_(q,eta) = ", fmt[fGeneral]];
Print["G_q = ", fmt[gGeneral]];
expectZero["F_general - expected", fGeneral - fExpected];
expectZero["G_general - expected", gGeneral - gExpected];

subbanner["3. Split-U specialization"];

qU = FullSimplify[-Sqrt[lambda0] rU, Assumptions -> $Assumptions];
etaU = FullSimplify[lambda0 rU, Assumptions -> $Assumptions];
fU = FullSimplify[fExpected /. {q -> qU, eta -> etaU}, Assumptions -> $Assumptions];
gU = FullSimplify[gExpected /. q -> qU, Assumptions -> $Assumptions];

fStage18 = FullSimplify[
  (9 delta + 11 xi)^4/(81 (1 - xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2),
  Assumptions -> $Assumptions
];
gStage19 = FullSimplify[9 xi (delta + xi)/(9 delta + 11 xi), Assumptions -> $Assumptions];

Print["F_U(xi,delta;R_U) = ", fmt[fU]];
Print["G_U(xi,delta;R_U) = ", fmt[gU]];
expectZero["F_U(R_U=1) - Stage18 F", (fU /. rU -> 1) - fStage18];
expectZero["G_U(R_U=1) - Stage19 G", (gU /. rU -> 1) - gStage19];

subbanner["4. Exact small-deformation expansion about the flat-U limit"];

fRatio = FullSimplify[Normal[Series[(fU /. rU -> 1 + eps)/fStage18, {eps, 0, 1}]], Assumptions -> $Assumptions];
gRatio = FullSimplify[Normal[Series[(gU /. rU -> 1 + eps)/gStage19, {eps, 0, 1}]], Assumptions -> $Assumptions];
hF = FullSimplify[(D[fU /. rU -> 1 + eps, eps] /. eps -> 0)/fStage18, Assumptions -> $Assumptions];
hG = FullSimplify[(D[gU /. rU -> 1 + eps, eps] /. eps -> 0)/gStage19, Assumptions -> $Assumptions];

Print["F_U / F_flat = ", fmt[fRatio], " + O(eps^2)"];
Print["G_U / G_flat = ", fmt[gRatio], " + O(eps^2)"];
Print["H_F = ", fmt[hF]];
Print["H_G = ", fmt[hG]];
expectZero["F_ratio - (1 + eps H_F)", Expand[fRatio - (1 + eps hF)]];
expectZero["G_ratio - (1 + eps H_G)", Expand[gRatio - (1 + eps hG)]];

Print[""];
Print["Stage 040 Mathematica audit passed."];

Exit[0];
