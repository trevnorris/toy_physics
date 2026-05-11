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

banner["I. DIMENSIONLESS REDUCTION OF THE COMPENSATION LAW"];

Clear[kS, kQ, lam, gS, gQ, rHat, gHat];
$Assumptions =
  Element[{kS, kQ, lam, gS, gQ, rHat, gHat}, Reals] &&
  kS > 0 && kQ > 0 && gS > 0;

law = Expand[gS^2*(kS*kQ + lam^2) - 4*(kS*gQ - lam*gS)^2];
lawRed = FullSimplify[
  (law /. {lam -> rHat*Sqrt[kS*kQ], gQ -> gHat*gS*Sqrt[kQ]/Sqrt[kS]})/(gS^2*kS*kQ),
  Assumptions -> $Assumptions
];

Print["Reduced law = ", fmt[lawRed]];
expectZero["dimensionless law", lawRed - (1 + rHat^2 - 4*(gHat - rHat)^2)];

banner["II. EXACT SOLUTION FOR THE MOUTH-COUPLING RATIO"];

gHatSol = Solve[1 + rHat^2 == 4*(gHat - rHat)^2, gHat, Reals];
If[Length[gHatSol] =!= 2, fail["expected two gHat branches", gHatSol]];
Print["gHat solutions = ", fmt[gHat /. gHatSol]];
expectZero["first branch check", 1 + rHat^2 - 4*((gHat /. gHatSol[[1]]) - rHat)^2];
expectZero["second branch check", 1 + rHat^2 - 4*((gHat /. gHatSol[[2]]) - rHat)^2];

banner["III. D/N-TUBE LENGTH IN TERMS OF THE SAME NORMALIZED HYBRIDIZATION"];

Clear[a, lW, rC];
$Assumptions = Element[{a, lW, rC}, Reals] && a > 0 && lW > 0 && rC > 0;

kappa0 = FullSimplify[4*lW^2/(Pi^2*a^2), Assumptions -> $Assumptions];
lSel = FullSimplify[lW /. First[Solve[kappa0 == (1 + rC)/3, lW, Reals]], Assumptions -> $Assumptions];

Print["L_W = ", fmt[lSel]];
expectZero["tube-length law", lSel - (Pi*a*Sqrt[(1 + rC)/3])/2];

banner["IV. EXPLICIT TRACTION LAW"];

Clear[zQ, mu0, cSound, tM, jS];
$Assumptions =
  Element[{kS, rHat, zQ, mu0, cSound, tM, jS, lW}, Reals] &&
  kS > 0 && zQ > 0 && mu0 > 0 && cSound > 0 && tM > 0 && jS > 0 && lW > 0;

kQExpr = zQ*Pi^2*cSound^2/(4*mu0*lW^2);
gQExpr = zQ*Pi/(Sqrt[2]*mu0*lW^(3/2));
gHatExpr = FullSimplify[gQExpr*Sqrt[kS]/(tM*jS*Sqrt[kQExpr]), Assumptions -> $Assumptions];

Print["gHat explicit = ", fmt[gHatExpr]];
expectZero[
  "gHat explicit simplification",
  gHatExpr - Sqrt[2*zQ*kS]/(tM*jS*cSound*Sqrt[mu0*lW])
];

tMPlus = FullSimplify[tM /. First[Solve[gHatExpr == rHat + Sqrt[1 + rHat^2]/2, tM, Reals]], Assumptions -> $Assumptions];
tMMinus = FullSimplify[tM /. First[Solve[gHatExpr == rHat - Sqrt[1 + rHat^2]/2, tM, Reals]], Assumptions -> $Assumptions];

Print["T_m (+ branch) = ", fmt[tMPlus]];
Print["T_m (- branch) = ", fmt[tMMinus]];

Print[""];
Print["Stage 119 Mathematica audit passed."];

Exit[0];
