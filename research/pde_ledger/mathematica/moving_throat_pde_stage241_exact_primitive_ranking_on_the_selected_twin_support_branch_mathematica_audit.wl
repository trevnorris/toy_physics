ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];
stripCE[expr_] := expr /. ConditionalExpression[e_, _] :> e;

pass[name_String] := Print["[ok] ", name];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

reduceResidual[expr_] := Module[{res},
  res = FullSimplify[Together[Expand[stripCE[expr]]], Assumptions -> $Assumptions];
  res = FullSimplify[stripCE[res], Assumptions -> $Assumptions];
  res
];

checkZero[name_String, expr_] := Module[{res},
  res = reduceResidual[expr];
  Print[name, " residual = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

checkPositive[name_String, expr_] := Module[{res, test},
  res = reduceResidual[expr];
  test = FullSimplify[res > 0, Assumptions -> $Assumptions];
  Print[name, " value = ", fmt[res]];
  If[TrueQ[test], pass[name], fail[name, res]];
];

takeNonzeroRoot[roots_List, name_String] := Module[{cleanRoots, nonzeroRoots},
  cleanRoots = FullSimplify[stripCE /@ roots, Assumptions -> $Assumptions];
  nonzeroRoots = Select[
    cleanRoots,
    ! TrueQ[FullSimplify[# == 0, Assumptions -> $Assumptions]] &
  ];
  If[Length[nonzeroRoots] != 1, fail[name, cleanRoots]];
  First[nonzeroRoots]
];

Print["Stage 241 - exact primitive ranking on the selected twin-support branch"];

Clear[beta, varrho, epsStar, sigma];

$Assumptions = (
  Element[{beta, varrho, epsStar, sigma}, Reals] &&
  0 < beta < 1 &&
  0 < varrho < 2/3 &&
  0 < epsStar < 1 &&
  sigma > 0
);

supportEquation = varrho == 2*(1 - epsStar)/3;
epsFromSupport = FullSimplify[
  stripCE[epsStar /. First[Solve[supportEquation, epsStar]]],
  Assumptions -> $Assumptions
];
epsSelected = epsFromSupport;
sigmaSelected = FullSimplify[
  (2*epsStar/(1 - epsStar)) /. epsStar -> epsSelected,
  Assumptions -> $Assumptions
];

checkZero[
  "M1 support law solved for epsStar",
  epsFromSupport - (1 - 3*varrho/2)
];
checkZero[
  "M1 selected sigma reduction",
  sigmaSelected - (4/(3*varrho) - 2)
];

sigmaSel = 4/(3*varrho) - 2;
checkZero[
  "M2 lower twin-window residual",
  sigmaSel - (1/varrho - 2) - 1/(3*varrho)
];
checkZero[
  "M2 upper twin-window residual",
  (2/varrho - 2) - sigmaSel - 2/(3*varrho)
];
checkPositive["M2 lower twin-window margin", sigmaSel - (1/varrho - 2)];
checkPositive["M2 upper twin-window margin", (2/varrho - 2) - sigmaSel];

den = 5*(1 - epsStar)^2 + 6*epsStar^2*(1 + beta^2);
wLambda = epsStar^2*(1 + beta^2)/den;
wZ = (1 - 2*epsStar + (2 + beta^2)*epsStar^2)/den;
wChi = 2*(1 - 2*epsStar + (2 + beta^2)*epsStar^2)/den;
wW = epsStar*(1 - epsStar)/den;
wUmag = beta*epsStar*(1 - epsStar)/den;

checkZero["M3 wChi = 2 wZ", wChi - 2*wZ];
checkZero[
  "M3 wZ - wLambda exact form",
  wZ - wLambda - (1 - epsStar)^2/den
];
checkZero[
  "M3 wZ - wW exact form",
  wZ - wW -
    (beta^2*epsStar^2 + 3*(epsStar - 1/2)^2 + 1/4)/den
];
checkZero[
  "M3 wW - wUmag exact form",
  wW - wUmag - epsStar*(1 - epsStar)*(1 - beta)/den
];

epsWRoots = epsStar /. Solve[Together[wW - wLambda] == 0, epsStar];
epsWCross = FullSimplify[
  takeNonzeroRoot[epsWRoots, "M4 wW-wLambda nonzero root"],
  Assumptions -> $Assumptions
];
checkZero[
  "M4 Solve[wW - wLambda == 0] nonzero root",
  epsWCross - 1/(2 + beta^2)
];

epsURoots = epsStar /. Solve[Together[wUmag - wLambda] == 0, epsStar];
epsUCross = FullSimplify[
  takeNonzeroRoot[epsURoots, "M4 wUmag-wLambda nonzero root"],
  Assumptions -> $Assumptions
];
checkZero[
  "M4 Solve[wUmag - wLambda == 0] nonzero root",
  epsUCross - beta/(1 + beta + beta^2)
];

varrhoWLambda = FullSimplify[
  stripCE[varrho /. First[Solve[epsSelected == epsWCross, varrho]]],
  Assumptions -> $Assumptions
];
varrhoULambda = FullSimplify[
  stripCE[varrho /. First[Solve[epsSelected == epsUCross, varrho]]],
  Assumptions -> $Assumptions
];
checkZero[
  "M4 support law gives varrhoWLambda",
  varrhoWLambda - 2*(1 + beta^2)/(3*(2 + beta^2))
];
checkZero[
  "M4 support law gives varrhoULambda",
  varrhoULambda - 2*(1 + beta^2)/(3*(1 + beta + beta^2))
];

selectedRules = {epsStar -> epsSelected};
wLambdaSel = FullSimplify[wLambda /. selectedRules, Assumptions -> $Assumptions];
wWSel = FullSimplify[wW /. selectedRules, Assumptions -> $Assumptions];
wUmagSel = FullSimplify[wUmag /. selectedRules, Assumptions -> $Assumptions];
wZSel = FullSimplify[wZ /. selectedRules, Assumptions -> $Assumptions];
wChiSel = FullSimplify[wChi /. selectedRules, Assumptions -> $Assumptions];
dPoly = (
  18*beta^2*varrho^2 - 24*beta^2*varrho + 8*beta^2
    + 33*varrho^2 - 24*varrho + 8
);

checkZero[
  "M5 factorized sign law for wLambda - wW",
  (wLambdaSel - wWSel) -
    (2 - 3*varrho)*(2 + beta^2)*(varrhoWLambda - varrho)/dPoly
];
checkZero[
  "M5 factorized sign law for wLambda - wUmag",
  (wLambdaSel - wUmagSel) -
    (2 - 3*varrho)*(1 + beta + beta^2)*(varrhoULambda - varrho)/dPoly
];

checkZero[
  "M6 varrhoULambda - varrhoWLambda exact form",
  varrhoULambda - varrhoWLambda -
    2*(1 + beta^2)*(1 - beta)/(3*(1 + beta + beta^2)*(2 + beta^2))
];
checkPositive["M6 varrhoULambda - varrhoWLambda positive", varrhoULambda - varrhoWLambda];
checkZero[
  "M6 two-thirds minus varrhoULambda exact form",
  2/3 - varrhoULambda - 2*beta/(3*(1 + beta + beta^2))
];
checkPositive["M6 two-thirds minus varrhoULambda positive", 2/3 - varrhoULambda];

betaMax = 2/11;
checkZero["M7 varrhoWLambda(beta=0)", (varrhoWLambda /. beta -> 0) - 1/3];
checkZero["M7 varrhoWLambda(beta=2/11)", (varrhoWLambda /. beta -> betaMax) - 125/369];
checkZero["M7 varrhoULambda(beta=0)", (varrhoULambda /. beta -> 0) - 2/3];
checkZero["M7 varrhoULambda(beta=2/11)", (varrhoULambda /. beta -> betaMax) - 250/441];
checkZero[
  "M7 d varrhoWLambda / d beta",
  D[varrhoWLambda, beta] - 4*beta/(3*(beta^2 + 2)^2)
];
checkZero[
  "M7 d varrhoULambda / d beta",
  D[varrhoULambda, beta] + 2*(1 - beta^2)/(3*(1 + beta + beta^2)^2)
];

betaSample = 1/10;
rhoWAtSample = FullSimplify[varrhoWLambda /. beta -> betaSample];
rhoUAtSample = FullSimplify[varrhoULambda /. beta -> betaSample];
rhoRegionI = FullSimplify[rhoWAtSample/2];
rhoRegionII = FullSimplify[(rhoWAtSample + rhoUAtSample)/2];
rhoRegionIII = FullSimplify[(rhoUAtSample + 2/3)/2];

regionIRules = {beta -> betaSample, varrho -> rhoRegionI};
checkPositive["M8 Region I wChi > wZ", wChiSel - wZSel /. regionIRules];
checkPositive["M8 Region I wZ > wLambda", wZSel - wLambdaSel /. regionIRules];
checkPositive["M8 Region I wLambda > wW", wLambdaSel - wWSel /. regionIRules];
checkPositive["M8 Region I wW > wUmag", wWSel - wUmagSel /. regionIRules];

regionIIRules = {beta -> betaSample, varrho -> rhoRegionII};
checkPositive["M8 Region II wChi > wZ", wChiSel - wZSel /. regionIIRules];
checkPositive["M8 Region II wZ > wW", wZSel - wWSel /. regionIIRules];
checkPositive["M8 Region II wW > wLambda", wWSel - wLambdaSel /. regionIIRules];
checkPositive["M8 Region II wLambda > wUmag", wLambdaSel - wUmagSel /. regionIIRules];

regionIIIRules = {beta -> betaSample, varrho -> rhoRegionIII};
checkPositive["M8 Region III wChi > wZ", wChiSel - wZSel /. regionIIIRules];
checkPositive["M8 Region III wZ > wW", wZSel - wWSel /. regionIIIRules];
checkPositive["M8 Region III wW > wUmag", wWSel - wUmagSel /. regionIIIRules];
checkPositive["M8 Region III wUmag > wLambda", wUmagSel - wLambdaSel /. regionIIIRules];

Print[""];
Print["Stage 241 Mathematica audit passed."];

Exit[0];
