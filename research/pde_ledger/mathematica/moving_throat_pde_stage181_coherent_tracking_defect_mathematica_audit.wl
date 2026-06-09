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

reduceExpr[expr_] := FullSimplify[Together[Cancel[expr]], Assumptions -> $Assumptions];

expectZero[name_String, expr_] := Module[{res},
  res = reduceExpr[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectNonZero[name_String, expr_] := Module[{res},
  res = reduceExpr[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], fail[name, res], pass[name]];
];

directional[expr_, pairs_List] := reduceExpr[
  Total[(D[expr, #[[1]]]*#[[2]]) & /@ pairs]
];

logDirectional[expr_, pairs_List] := reduceExpr[directional[expr, pairs]/expr];

banner["STAGE 181 - COHERENT TRACKING-BRANCH DEFECT LAW"];

Clear[
  grav, light, sound, radius, etaSlack, wallBlock, uSplit, joint, supportShare,
  wallOverlap, freqSq
];

$Assumptions = Element[
    {grav, light, sound, radius, etaSlack, wallBlock, uSplit, joint, supportShare,
     wallOverlap, freqSq},
    Reals
  ] && grav > 0 && light > 0 && sound > 0 && radius > 0 &&
  wallOverlap > 0 && freqSq > 0;

shapeConstant = reduceExpr[27*Pi^2*grav*sound^5/(20*radius^5*light^5)];
lambdaBranch = reduceExpr[shapeConstant*freqSq];
splitMultiplier = reduceExpr[1 - 2*uSplit/(11*(1 + uSplit))];
splitBlock = reduceExpr[wallBlock*splitMultiplier];

targetFromSelectedBranch = reduceExpr[
  lambdaBranch*(1 - etaSlack)*(1 - splitBlock)^2/(wallOverlap*(1 + joint)^2)
];

shapeFromDirectPort = reduceExpr[
  wallOverlap*(1 + joint)^2/(freqSq*(1 - splitBlock)^2)
];
shapeFromSelectedDemand = reduceExpr[
  shapeConstant*(1 - etaSlack)/targetFromSelectedBranch
];

expectZero[
  "direct-selected transfer-shape identity",
  shapeFromDirectPort - shapeFromSelectedDemand
];

mixedLegMass = reduceExpr[
  8*wallOverlap*(1 + joint)^2/(Pi^2*(1 - etaSlack)*(1 - splitBlock))
];
coherentSupportMass = reduceExpr[
  8*supportShare*wallOverlap*(1 + joint)^2/
    (Pi^2*(1 - etaSlack)*(1 - supportShare*splitBlock))
];
supportMultiplier = reduceExpr[
  1 + supportShare*(1 - splitBlock)/(1 - supportShare*splitBlock)
];

loadedMassRaw = mixedLegMass + coherentSupportMass;
loadedProductRaw = 8*lambdaBranch*(1 - splitBlock)*supportMultiplier/Pi^2;
loadedTargetRaw = loadedProductRaw/loadedMassRaw;
loadedShapeRaw = shapeConstant*(1 - etaSlack)/loadedTargetRaw;

expectZero[
  "support-loaded R_target reconstruction",
  loadedTargetRaw - targetFromSelectedBranch
];
expectZero[
  "support-loaded T^2 reconstruction",
  loadedShapeRaw - shapeFromDirectPort
];
expectZero[
  "d/dzeta ln T^2 (support-loaded route)",
  D[Log[loadedShapeRaw], supportShare]
];
expectZero[
  "d/dzeta ln R_target (support-loaded route)",
  D[Log[loadedTargetRaw], supportShare]
];

Clear[spoilerAmplitude];
spoiledMassRaw = loadedMassRaw + spoilerAmplitude*supportShare*mixedLegMass;
spoiledTargetRaw = loadedProductRaw/spoiledMassRaw;
spoiledSupportSlope = reduceExpr[
  D[Log[spoiledTargetRaw], supportShare] /. spoilerAmplitude -> 1
];
expectNonZero["spoiled d/dzeta ln R_target", spoiledSupportSlope];

banner["Weak-axisymmetric drift transport"];

Clear[
  overlapDrift, frequencyDrift, jointDrift, bareBlockDrift, uSplitDrift,
  etaSlackDrift
];

$Assumptions = Element[
    {grav, light, sound, radius, etaSlack, wallBlock, uSplit, joint, supportShare,
     wallOverlap, freqSq, overlapDrift, frequencyDrift, jointDrift, bareBlockDrift,
     uSplitDrift, etaSlackDrift},
    Reals
  ] && grav > 0 && light > 0 && sound > 0 && radius > 0 &&
  wallOverlap > 0 && freqSq > 0;

placementDrifts = {
  {wallOverlap, wallOverlap*overlapDrift},
  {freqSq, freqSq*frequencyDrift},
  {joint, jointDrift},
  {wallBlock, bareBlockDrift},
  {uSplit, uSplitDrift},
  {etaSlack, etaSlackDrift}
};

splitVariableDrifts = {
  {wallBlock, bareBlockDrift},
  {uSplit, uSplitDrift}
};

splitDriftFromFullExpression = directional[splitBlock, splitVariableDrifts];
splitDriftFromProductLedger = reduceExpr[
  bareBlockDrift*splitMultiplier +
    wallBlock*directional[splitMultiplier, {{uSplit, uSplitDrift}}]
];

expectZero[
  "split-blocking drift eps_1",
  splitDriftFromFullExpression - splitDriftFromProductLedger
];

epsilonOne = reduceExpr[splitDriftFromFullExpression];
Print["epsilon_1 = ", fmt[epsilonOne]];

xiOneFromShape = logDirectional[shapeFromDirectPort, placementDrifts];
xiOneLaw = reduceExpr[
  overlapDrift - frequencyDrift + 2*jointDrift/(1 + joint) +
    2*epsilonOne/(1 - splitBlock)
];

expectZero[
  "Xi_1 derived from transfer differential matches defect law",
  xiOneFromShape - xiOneLaw
];

rOneFromTarget = logDirectional[targetFromSelectedBranch, placementDrifts];
rOneLaw = reduceExpr[
  frequencyDrift - etaSlackDrift/(1 - etaSlack) - overlapDrift -
    2*jointDrift/(1 + joint) - 2*epsilonOne/(1 - splitBlock)
];

expectZero[
  "R_1 derived from selected demand matches defect law",
  rOneFromTarget - rOneLaw
];
expectZero[
  "selected-branch identity",
  xiOneFromShape + etaSlackDrift/(1 - etaSlack) + rOneFromTarget
];

xiOneFromLoadedRoute = logDirectional[loadedShapeRaw, placementDrifts];
expectZero[
  "support-loaded Xi_1 matches direct differential",
  xiOneFromLoadedRoute - xiOneFromShape
];
expectZero[
  "d/dzeta Xi_1 (support-loaded route)",
  D[xiOneFromLoadedRoute, supportShare]
];

Print["Xi_1 = ", fmt[xiOneFromShape]];
Print["R_1  = ", fmt[rOneFromTarget]];

banner["Tracking-factor drift"];

trackingFactorDirect = reduceExpr[(1 + joint/(1 + uSplit))/(1 + joint)];
trackingFactorFactored = reduceExpr[(1 + joint + uSplit)/((1 + joint)*(1 + uSplit))];

expectZero[
  "tracking-factor quotient form",
  trackingFactorDirect - trackingFactorFactored
];

trackingDrifts = {
  {joint, jointDrift},
  {uSplit, uSplitDrift}
};

thetaFromQuotient = logDirectional[trackingFactorDirect, trackingDrifts];
thetaFromFactorLedger = reduceExpr[
  logDirectional[1 + joint + uSplit, trackingDrifts] -
    logDirectional[1 + joint, trackingDrifts] -
    logDirectional[1 + uSplit, trackingDrifts]
];

expectZero[
  "tracking-factor drift",
  thetaFromQuotient - thetaFromFactorLedger
];

thetaOne = reduceExpr[thetaFromQuotient];
Print["Theta_1 = ", fmt[thetaOne]];

banner["Support-blindness consequence"];

xiSupportRigid = reduceExpr[xiOneFromShape /. {jointDrift -> 0, uSplitDrift -> 0}];
thetaSupportRigid = reduceExpr[thetaOne /. {jointDrift -> 0, uSplitDrift -> 0}];

Print["Xi_1 with chi1=deltaU1=0 = ", fmt[xiSupportRigid]];
Print["Theta_1 with chi1=deltaU1=0 = ", fmt[thetaSupportRigid]];

If[TrueQ[xiSupportRigid === 0], fail["support-rigid specialization unexpectedly killed Xi_1"]];

Print[""];
Print["Conclusion:"];
Print["  The coherent support ratio zeta drops out identically from T^2 and R_target."];
Print["  The grouped weak-axisymmetric defect is carried only by"];
Print["  Z_W, Omega_W^2, chi_0, eps_W, and delta_U."];
Print["  Exact tracking-factor rigidity is not sufficient to kill Xi_1."];

Print[""];
Print["Stage 181 Mathematica audit passed."];

Exit[0];
