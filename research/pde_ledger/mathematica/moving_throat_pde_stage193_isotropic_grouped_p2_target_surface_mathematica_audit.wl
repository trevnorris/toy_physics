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

stripCE[expr_] := expr /. ConditionalExpression[e_, _] :> e;

cleanScalar[expr_] := FullSimplify[
  Together[Expand[stripCE[expr]]],
  Assumptions -> $Assumptions
];

cleanTensor[expr_] := If[
  ListQ[expr],
  Map[cleanScalar, expr, {Length[Dimensions[expr]]}],
  cleanScalar[expr]
];

zeroTensorQ[expr_] := If[
  ListQ[expr],
  And @@ (TrueQ[# === 0] & /@ Flatten[expr]),
  TrueQ[expr === 0]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanTensor[expr];
  Print[name, " = ", fmt[res]];
  If[zeroTensorQ[res], pass[name], fail[name, res]];
];

expectTrue[name_String, claim_] := Module[{res},
  res = FullSimplify[stripCE[claim], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 193 -- MATHEMATICA ISOTROPIC GROUPED-P2 TARGET SURFACE AUDIT"];

Clear[
  x20, x21, x22, xMean, xAxial, xSplit,
  nuTwo, nuFour, w, d0, d2, d4,
  d0geom, d2iso, chi, c20, c21, c22
];

$Assumptions = (
  Element[
    {
      x20, x21, x22, xMean, xAxial, xSplit,
      nuTwo, nuFour, w, d0, d2, d4,
      d0geom, d2iso, chi, c20, c21, c22
    },
    Reals
  ] &&
  d0 != 0 && d2 != 0 && d0geom != 0 && d2iso != 0
);

subbanner["I. Grouped trace/anomaly projector"];

laneVector = {x20, x21, x22};
projector = {
  {1/5, 2/5, 2/5},
  {1/5, -1/10, -1/10},
  {0, 1/2, -1/2}
};

groupCoordinates = cleanTensor[projector . laneVector];
inverseProjector = cleanTensor[Inverse[projector]];
inverseTarget = {
  {1, 4, 0},
  {1, -1, 1},
  {1, -1, -1}
};

Print["projector lanes -> {xbar, ax, bx} = ", fmt[groupCoordinates]];
Print["inverse projector = ", fmt[inverseProjector]];
expectZero["inverse projector - SymPy grouped inverse target", inverseProjector - inverseTarget];
expectZero["projector inverse reconstructs lanes", inverseProjector . (projector . laneVector) - laneVector];

commonNu2 = projector . {nuTwo, nuTwo, nuTwo};
commonNu4 = projector . {nuFour, nuFour, nuFour};
expectZero["a2,b2 vanish on common-lane p2 branch", commonNu2[[2 ;; 3]]];
expectZero["a4,b4 vanish on common-lane p4 branch", commonNu4[[2 ;; 3]]];

subbanner["II. One-pole conservative identity from response derivatives"];

denominator[w_] := d0 + d2*w^2 + d4*w^4;
responseShape[w_] := d0/denominator[w];

nu2FromD = cleanScalar[(D[responseShape[w], {w, 2}]/2) /. w -> 0];
nu4FromD = cleanScalar[(D[responseShape[w], {w, 4}]/24) /. w -> 0];
deltaPole = cleanScalar[nu4FromD - 4*nu2FromD^2];

Print["nu2 from d^2 Y/dw^2 = ", fmt[nu2FromD]];
Print["nu4 from d^4 Y/dw^4 = ", fmt[nu4FromD]];
Print["Delta_pole = ", fmt[deltaPole]];
expectZero["nu2 + d2/d0", nu2FromD + d2/d0];
expectZero["nu4 - (d2^2 - d0 d4)/d0^2", nu4FromD - (d2^2 - d0*d4)/d0^2];
expectZero["Delta_pole + (3 d2^2 + d0 d4)/d0^2", deltaPole + (3*d2^2 + d0*d4)/d0^2];

onePoleD4 = cleanScalar[d4 /. First[Solve[deltaPole == 0, d4]]];
expectZero["Solve[Delta_pole==0] gives d0 d4 + 3 d2^2 = 0", onePoleD4 + 3*d2^2/d0];
expectZero["Delta_pole on one-pole surface", deltaPole /. d4 -> onePoleD4];

subbanner["III. One-parameter conservative carrier through O(w^4)"];

omegaQ2 = cleanScalar[-d0/(4*d2)];
carrier[w_] := 3/4 + (1/4)/(1 - w^2/omegaQ2);
carrierSeries = cleanScalar[Normal[Series[carrier[w], {w, 0, 4}]]];
carrierCoeff2 = cleanScalar[Coefficient[carrierSeries, w, 2]];
carrierCoeff4 = cleanScalar[Coefficient[carrierSeries, w, 4]];
targetSeries = cleanScalar[1 + nu2FromD*w^2 + (nu4FromD /. d4 -> onePoleD4)*w^4];

Print["Omega_Q^2 = ", fmt[omegaQ2]];
Print["carrier series through w^4 = ", fmt[carrierSeries]];
expectZero["carrier w^2 coefficient - nu2", carrierCoeff2 - nu2FromD];
expectZero["carrier w^4 coefficient - one-pole nu4", carrierCoeff4 - (nu4FromD /. d4 -> onePoleD4)];
expectZero["carrier series - SymPy one-pole target series", carrierSeries - targetSeries];

subbanner["IV. Scalar/geometry firewall by native block Schur complement"];

scalarBlock = {{d0geom}};
mixSeed = {{c20, c21, c22}};
l2Block = d2iso*IdentityMatrix[3];
upperMix = chi*mixSeed;
lowerMix = Transpose[upperMix];

blockOperator = ArrayFlatten[{
  {scalarBlock, upperMix},
  {lowerMix, l2Block}
}];
schurEff = cleanTensor[l2Block - lowerMix . Inverse[scalarBlock] . upperMix];
schurTarget = cleanTensor[d2iso*IdentityMatrix[3] - chi^2*Transpose[mixSeed] . mixSeed/d0geom];

Print["D_block(chi) = ", fmt[blockOperator]];
Print["D_eff,l=2(chi) from Schur complement = ", fmt[schurEff]];
expectZero["upper off-diagonal derivative is C", D[blockOperator[[1, 2 ;; 4]], chi] - mixSeed[[1]]];
expectZero["upper off-diagonal second derivative vanishes", D[blockOperator[[1, 2 ;; 4]], {chi, 2}]];
expectZero["Schur complement - (D2 I - chi^2 C^T C / D0scalar)", schurEff - schurTarget];
expectZero["det block - det scalar det Schur", Det[blockOperator] - Det[scalarBlock]*Det[schurEff]];
expectZero["d/dchi D_eff at chi=0 (linear-order firewall)", D[schurEff, chi] /. chi -> 0];
expectZero["D_eff at chi=0 - D2 I", (schurEff /. chi -> 0) - d2iso*IdentityMatrix[3]];

banner["STAGE 193 MATHEMATICA LEDGER"];
Print["1. The grouped projector has the exact inverse used by the SymPy witness, and"];
Print["   common grouped lanes force a2=b2=a4=b4=0."];
Print["2. Native response derivatives give Delta_pole = -(3 d2^2 + d0 d4)/d0^2,"];
Print["   hence the one-pole surface d0 d4 + 3 d2^2 = 0."];
Print["3. The one-parameter carrier 3/4 + 1/4 (1 - w^2/Omega_Q^2)^(-1) matches"];
Print["   the one-pole conservative series through O(w^4)."];
Print["4. A block operator with linear chi mixing gives a Schur correction"];
Print["   D_eff,l=2 = D2 I - chi^2 C^T C / D0scalar, with no O(chi) term."];

Exit[0];
