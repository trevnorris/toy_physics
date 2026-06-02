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

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[stripCE[expr]]], Assumptions -> $Assumptions];
  res = FullSimplify[stripCE[res], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 121 - GEOMETRIC SELECTION OF r"];

Clear[ell, throatScale, tubeLength, radius, soundSpeed];

baseAssumptions = (
  Element[{ell, throatScale, tubeLength, radius, soundSpeed}, Reals] &&
  ell > 0 && throatScale > 0 && tubeLength > 0 && soundSpeed > 0
);
branchAssumptions = (
  baseAssumptions &&
  ell/throatScale > Pi/(2*Sqrt[3]) &&
  tubeLength/throatScale > Pi/(2*Sqrt[3]) &&
  -3*Pi^2*throatScale^2 + 36*ell^2 > 0 &&
  -3*Pi^2*throatScale^2 + 36*tubeLength^2 > 0
);
$Assumptions = branchAssumptions;

stage99TubeLength = (Pi*throatScale/2)*Sqrt[(1 + radius^2)/3];
radiusRoots = stripCE /@ (radius /. Solve[tubeLength == stage99TubeLength, radius, Reals]);
positiveRadiusRoots = Select[
  radiusRoots,
  TrueQ[FullSimplify[# > 0, Assumptions -> branchAssumptions]] &
];

If[Length[positiveRadiusRoots] =!= 1,
  fail["positive branch selection", radiusRoots]
];

derivedRadius = FullSimplify[
  First[positiveRadiusRoots] /. tubeLength -> ell,
  Assumptions -> branchAssumptions
];
Print["r_geom(L/a) = ", fmt[derivedRadius]];

expectZero[
  "r_geom closed-form (explicit)",
  derivedRadius^2 - (12*ell^2/(Pi^2*throatScale^2) - 1)
];

expectZero[
  "tube-length relation",
  (stage99TubeLength /. radius -> derivedRadius) - ell
];

ratioF1 = 37/20;
rF1Derived = FullSimplify[
  derivedRadius /. ell -> ratioF1*throatScale,
  Assumptions -> baseAssumptions
];
rF1Target = Sqrt[4107 - 100*Pi^2]/(10*Pi);
Print["r_F1 = ", fmt[rF1Derived]];
Print["r_F1 numeric = ", N[rF1Derived, 20]];

expectZero[
  "r_F1 symbolic value at L/a = 37/20",
  FullSimplify[rF1Derived - rF1Target, Assumptions -> baseAssumptions]
];

rcF1Derived = FullSimplify[rF1Derived^2, Assumptions -> baseAssumptions];
rcF1Target = 4107/(100*Pi^2) - 1;
Print["r_c(F1) = ", fmt[rcF1Derived]];
Print["r_c(F1) numeric = ", N[rcF1Derived, 20]];

expectZero[
  "r_c(F1) symbolic value",
  rcF1Derived - rcF1Target
];

halfWavePole = Pi*soundSpeed/(2*tubeLength);
omegaAtEqualTubeLength = FullSimplify[
  halfWavePole /. tubeLength -> ell,
  Assumptions -> baseAssumptions
];
Print["Omega_W(L_W = L) = ", fmt[omegaAtEqualTubeLength]];

expectZero[
  "Omega_W identification at L_W = L",
  omegaAtEqualTubeLength - Pi*soundSpeed/(2*ell)
];

thresholdAspect = Pi/(2*Sqrt[3]);
Print["existence threshold L/a >= ", fmt[thresholdAspect], " ~= ", N[thresholdAspect, 20]];

$Assumptions = baseAssumptions;
expectZero[
  "r_geom vanishes at existence threshold",
  FullSimplify[derivedRadius /. ell -> thresholdAspect*throatScale,
               Assumptions -> baseAssumptions]
];

Print[""];
Print["STAGE 121 MATHEMATICA AUDIT PASSED"];

Exit[0];
