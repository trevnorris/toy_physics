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
fmt[expr_] := StringReplace[ToString[InputForm[expr]], "Global`" -> ""];
stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

cleanExpr[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Cancel[Together[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

zeroQ[expr_] := And @@ (TrueQ[# === 0] & /@ Flatten[{expr}]);

expectZero[name_String, expr_] := Module[{res},
  res = cleanExpr[expr];
  Print[name, " = ", fmt[res]];
  If[zeroQ[res], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

paperSymbol[name_String] := Quiet[Symbol["Global`" <> name], General::shdw];

{
  kWall, mWall, cBdg, varpi, OmegaU, OmegaW, rMix, gU, gW,
  omega, piOut, gammaOut, Jq, JU, JW, sq, sU, sW, sourceS,
  x, kappa, betaQ, betaU, betaW, y
} = paperSymbol /@ {
  "K", "M", "C", "varpi", "OmegaU", "OmegaW", "R", "GU", "GW",
  "omega", "Pi", "Gamma", "Jq", "JU", "JW", "sq", "sU", "sW", "S",
  "x", "kappa", "betaQ", "betaU", "betaW", "y"
};

$Assumptions = (
  Element[
    {
      kWall, mWall, cBdg, varpi, OmegaU, OmegaW, rMix, gU, gW,
      omega, gammaOut, Jq, JU, JW, sq, sU, sW, sourceS,
      betaQ, betaU, betaW
    },
    Reals
  ]
  && gammaOut >= 0
  && x > 0
  && kappa > 0
);

sampleRules = {
  kWall -> 11,
  mWall -> 2,
  cBdg -> 1,
  varpi -> 5,
  OmegaU -> 3,
  OmegaW -> 4,
  rMix -> 2,
  gU -> 1,
  gW -> 2,
  omega -> 1/2,
  Jq -> 1,
  JU -> 2,
  JW -> 1,
  gammaOut -> 1/10
};
sampleConservativeRules = Join[sampleRules, {piOut -> 0}];

banner["STAGE 220 -- DYNAMIC MIXED-PORT KERNEL MATHEMATICA AUDIT"];

subbanner["Symbol setup"];

KB = kWall - mWall omega^2 - cBdg^2/(varpi^2 - omega^2);
Aw = OmegaU^2 - omega^2;
Ww = OmegaW^2 - omega^2 - piOut;

Kdyn = {
  {KB, -gU, -gW},
  {-gU, Aw, -rMix},
  {-gW, -rMix, Ww}
};

DeltaPi = Aw Ww - rMix^2;
QPi = gU^2 Ww + 2 gU gW rMix + gW^2 Aw;
DPi = KB - QPi/DeltaPi;

Print["KB = ", fmt[KB]];
Print["Aw = ", fmt[Aw]];
Print["Ww = ", fmt[Ww]];
Print["DeltaPi = ", fmt[DeltaPi]];
Print["QPi = ", fmt[QPi]];
Print["DPi = ", fmt[DPi]];

subbanner["M1. Determinant identity"];

detNative = Det[Kdyn];
expectZero["M1 Det[Kdyn] - DeltaPi DPi", detNative - DeltaPi DPi];

subbanner["M2. Static reduction"];

Kstar = kWall - cBdg^2/varpi^2;
Delta0 = OmegaU^2 OmegaW^2 - rMix^2;
Q0 = gU^2 OmegaW^2 + 2 gU gW rMix + gW^2 OmegaU^2;
D0 = Kstar - Q0/Delta0;
staticRules = {omega -> 0, piOut -> 0};
staticResiduals = {
  (KB /. staticRules) - Kstar,
  (Aw /. staticRules) - OmegaU^2,
  (Ww /. staticRules) - OmegaW^2,
  (DeltaPi /. staticRules) - Delta0,
  (QPi /. staticRules) - Q0,
  (DPi /. staticRules) - D0
};
expectZero["M2 static Stage 219 bundle residuals", staticResiduals];

subbanner["M3. Native inverse entries"];

Kinv = Map[Cancel[Together[#]] &, Inverse[Kdyn], {2}];

chiQQ = 1/DPi;
chiQU = (gU Ww + rMix gW)/(DeltaPi DPi);
chiQW = (Aw gW + rMix gU)/(DeltaPi DPi);
chiUU = (KB Ww - gW^2)/(DeltaPi DPi);
chiUW = (KB rMix + gU gW)/(DeltaPi DPi);
chiWW = (KB Aw - gU^2)/(DeltaPi DPi);

Scan[
  Function[row, expectZero[row[[1]], row[[2]] - row[[3]]]],
  {
    {"M3 chi_qq native inverse entry", chiQQ, Kinv[[1, 1]]},
    {"M3 chi_qU native inverse entry", chiQU, Kinv[[1, 2]]},
    {"M3 chi_qW native inverse entry", chiQW, Kinv[[1, 3]]},
    {"M3 chi_UU native inverse entry", chiUU, Kinv[[2, 2]]},
    {"M3 chi_UW native inverse entry", chiUW, Kinv[[2, 3]]},
    {"M3 chi_WW native inverse entry", chiWW, Kinv[[3, 3]]}
  }
];

subbanner["M4. Susceptibility quadratic law"];

Jvec = {Jq, JU, JW};
Vmix = -1/2 (Jvec . Kinv . Jvec);
VmixExpected = -1/2 (
  chiQQ Jq^2
  + 2 chiQU Jq JU
  + 2 chiQW Jq JW
  + chiUU JU^2
  + 2 chiUW JU JW
  + chiWW JW^2
);
expectZero["M4 Vmix native matrix form - susceptibility form", Vmix - VmixExpected];

subbanner["M5. Collinear-source factorization"];

Jcol = sourceS {sq, sU, sW};
Ns = (
  DeltaPi sq^2
  + 2 (gU Ww + rMix gW) sq sU
  + 2 (Aw gW + rMix gU) sq sW
  + (KB Ww - gW^2) sU^2
  + 2 (KB rMix + gU gW) sU sW
  + (KB Aw - gU^2) sW^2
);
chiS = Ns/(DeltaPi DPi);
expectZero[
  "M5 collinear source response",
  -1/2 (Jcol . Kinv . Jcol) + 1/2 chiS sourceS^2
];

subbanner["M6. Primitive product-family theorem"];

Jprim = {
  betaQ x^-3,
  betaU Exp[-2 kappa x]/x,
  betaW Exp[-2 kappa x]/x
};
Vprim = -1/2 (Jprim . Kinv . Jprim);

C6 = chiQQ betaQ^2;
C4 = chiQU betaQ betaU + chiQW betaQ betaW;
C2 = chiUU betaU^2 + 2 chiUW betaU betaW + chiWW betaW^2;
VprimExpected = -1/2 (
  C6 x^-6
  + 2 C4 Exp[-2 kappa x]/x^4
  + C2 Exp[-4 kappa x]/x^2
);

expectZero["M6 primitive product-family formula", Vprim - VprimExpected];

VprimY = Expand[Vprim /. {Exp[-4 kappa x] -> y^2, Exp[-2 kappa x] -> y}];
VprimYCollected = Collect[VprimY, {x^-6, y/x^4, y^2/x^2}, FullSimplify];
extractedFamilyResiduals = {
  Coefficient[VprimYCollected, x^-6] + C6/2,
  Coefficient[VprimYCollected, y/x^4] + C4,
  Coefficient[VprimYCollected, y^2/x^2] + C2/2
};
expectZero["M6 extracted family coefficients", extractedFamilyResiduals];

supportPolynomial = Expand[Cancel[Together[-2 VprimY DeltaPi DPi x^6]]];
supportRules = CoefficientRules[supportPolynomial, {x, y}];
shiftedSupportSet = Sort[supportRules[[All, 1]]];
laurentSupportSet = Sort[({#[[1]] - 6, #[[2]]} & /@ shiftedSupportSet)];
expectedSupportSet = Sort[{{-6, 0}, {-4, 1}, {-2, 2}}];
Print["M6 Laurent support = ", fmt[laurentSupportSet]];
expectTrue["M6 no additional x-y monomial families", laurentSupportSet === expectedSupportSet];

subbanner["M7. Outgoing-port derivative identity"];

TJ = chiQW Jq + chiUW JU + chiWW JW;
dVdPi = D[VmixExpected, piOut];
expectTrue[
  "M7 self-test VmixExpected depends on Pi",
  (dVdPi /. sampleConservativeRules) != 0
];
expectZero["M7 dVmix/dPi + 1/2 TJ^2", dVdPi + 1/2 TJ^2];

subbanner["M8. Linear outgoing correction"];

TJ0 = TJ /. piOut -> 0;
linearOutgoing = (dVdPi /. piOut -> 0) I gammaOut;
linearExpected = -1/2 I gammaOut TJ0^2;
expectZero["M8 linear outgoing correction", linearOutgoing - linearExpected];

subbanner["M9. Phase-lag no-go and absorbed power"];

dV1 = -1/2 I gammaOut TJ0^2;
phaseReal = ComplexExpand[Re[dV1]];
phaseImag = ComplexExpand[Im[dV1]];
Pabs = -omega phaseImag;

expectTrue["M9 Re[dV1] is identically zero", phaseReal === 0];
expectZero["M9 Pabs perfect-square form", Pabs - omega gammaOut/2 TJ0^2];
expectTrue[
  "M9 off-pole sample is dissipative",
  (DeltaPi /. sampleConservativeRules) != 0
    && (DPi /. sampleConservativeRules) != 0
    && (phaseReal /. sampleConservativeRules) == 0
    && (phaseImag /. sampleConservativeRules) != 0
    && (Pabs /. sampleConservativeRules) > 0
];

Print[""];
Print["All Stage 220 Mathematica checks passed."];
