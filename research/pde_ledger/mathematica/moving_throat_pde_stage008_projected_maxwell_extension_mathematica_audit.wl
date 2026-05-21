ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

Print["STAGE 008 PROJECTED MAXWELL EXTENSION MATHEMATICA AUDIT"];

Clear[w, lambda, sigma, xi, mu0, R, iwzSym, iwhSym, zIntSym];
$Assumptions =
  lambda > 0 && sigma > 0 && xi > 0 && mu0 > 0 && R > 0 &&
    iwzSym > 0 && iwhSym > 0 && zIntSym > 0 && Element[w, Reals];

xiEffProj[iwz_, iwh_] := xi*iwz/iwh;
invXiEffProj[iwz_, iwh_] := iwh/(xi*iwz);

m1Residual =
  FullSimplify[xiEffProj[iwzSym, iwhSym]*invXiEffProj[iwzSym, iwhSym] - 1];
Print["M1 reciprocal residual = ", fmt[m1Residual]];
If[FullSimplify[xiEffProj[iwzSym, iwhSym]*invXiEffProj[iwzSym, iwhSym] - 1] =!= 0,
  Print["FAIL: M1 reciprocal identity"]; Exit[1]
];
Print["PASS: M1 reciprocal identity"];

gaussianWeight = Exp[-w^2/sigma^2]/(Sqrt[Pi]*sigma);
localizedGaussian = Exp[-w^2/lambda^2];
lorentzWeight = (sigma/Pi)/(w^2 + sigma^2);

gaussNorm =
  Integrate[gaussianWeight, {w, -Infinity, Infinity},
    Assumptions -> lambda > 0 && sigma > 0];
gaussOverlap =
  Integrate[gaussianWeight*localizedGaussian, {w, -Infinity, Infinity},
    Assumptions -> lambda > 0 && sigma > 0];
gaussGaugeWeight =
  Integrate[gaussianWeight*localizedGaussian, {w, -Infinity, Infinity},
    Assumptions -> lambda > 0 && sigma > 0];
m2aNormResidual = FullSimplify[gaussNorm - 1];
m2aGaugeResidual = FullSimplify[xi*gaussOverlap/gaussGaugeWeight - xi];
Print["M2 Pair A normalization residual = ", fmt[m2aNormResidual]];
Print["M2 Pair A H=Z residual = ", fmt[m2aGaugeResidual]];
If[FullSimplify[gaussNorm - 1] =!= 0,
  Print["FAIL: M2 Pair A normalization"]; Exit[1]
];
If[FullSimplify[xi*gaussOverlap/gaussGaugeWeight - xi] =!= 0,
  Print["FAIL: M2 Pair A H=Z gauge invariance"]; Exit[1]
];
Print["PASS: M2 Pair A Gaussian-Gaussian"];

lorentzNorm =
  Integrate[lorentzWeight, {w, -Infinity, Infinity},
    Assumptions -> lambda > 0 && sigma > 0];
lorentzOverlap =
  Integrate[lorentzWeight*localizedGaussian, {w, -Infinity, Infinity},
    Assumptions -> lambda > 0 && sigma > 0];
lorentzGaugeWeight =
  Integrate[lorentzWeight*localizedGaussian, {w, -Infinity, Infinity},
    Assumptions -> lambda > 0 && sigma > 0];
m2bNormResidual = FullSimplify[lorentzNorm - 1];
m2bGaugeResidual = FullSimplify[xi*lorentzOverlap/lorentzGaugeWeight - xi];
Print["M2 Pair B normalization residual = ", fmt[m2bNormResidual]];
Print["M2 Pair B H=Z residual = ", fmt[m2bGaugeResidual]];
If[FullSimplify[lorentzNorm - 1] =!= 0,
  Print["FAIL: M2 Pair B normalization"]; Exit[1]
];
If[FullSimplify[xi*lorentzOverlap/lorentzGaugeWeight - xi] =!= 0,
  Print["FAIL: M2 Pair B H=Z gauge invariance"]; Exit[1]
];
Print["PASS: M2 Pair B Lorentzian-Gaussian"];

zArea =
  Integrate[localizedGaussian, {w, -Infinity, Infinity},
    Assumptions -> lambda > 0];
gaussMatchedSourceOverlap =
  Integrate[gaussianWeight*(localizedGaussian/zArea), {w, -Infinity, Infinity},
    Assumptions -> lambda > 0 && sigma > 0];
m3Residual =
  FullSimplify[mu0*gaussMatchedSourceOverlap/gaussOverlap - mu0/zArea];
Print["M3 matched source residual = ", fmt[m3Residual]];
If[FullSimplify[mu0*gaussMatchedSourceOverlap/gaussOverlap - mu0/zArea] =!= 0,
  Print["FAIL: M3 mu0_eff_proj matched source"]; Exit[1]
];
Print["PASS: M3 mu0_eff_proj matched source"];

matchedWeight = localizedGaussian/(Sqrt[Pi]*lambda);
zSelfArea =
  Integrate[localizedGaussian^2, {w, -Infinity, Infinity},
    Assumptions -> lambda > 0];
matchedOverlap =
  Integrate[matchedWeight*localizedGaussian, {w, -Infinity, Infinity},
    Assumptions -> lambda > 0];
matchedNorm =
  Integrate[matchedWeight, {w, -Infinity, Infinity},
    Assumptions -> lambda > 0];
matchedDeltaSource =
  Integrate[matchedWeight*DiracDelta[w], {w, -Infinity, Infinity},
    Assumptions -> lambda > 0];
matchedDistributedSource =
  Integrate[matchedWeight*(localizedGaussian/zArea), {w, -Infinity, Infinity},
    Assumptions -> lambda > 0];
m4Residuals = {
  FullSimplify[zArea - Sqrt[Pi]*lambda],
  FullSimplify[zSelfArea - Sqrt[Pi/2]*lambda],
  FullSimplify[matchedOverlap - Sqrt[2]/2],
  FullSimplify[xi*matchedOverlap/matchedOverlap - xi],
  FullSimplify[xi*matchedOverlap/matchedNorm - xi/Sqrt[2]],
  FullSimplify[(mu0*matchedDeltaSource/matchedOverlap)/(mu0/zArea) - Sqrt[2]],
  FullSimplify[mu0*matchedDistributedSource/matchedOverlap - mu0/zArea]
};
Print["M4 residuals = ", fmt[m4Residuals]];
If[FullSimplify[zArea - Sqrt[Pi]*lambda] =!= 0,
  Print["FAIL: M4 Gaussian Z_int"]; Exit[1]
];
If[FullSimplify[zSelfArea - Sqrt[Pi/2]*lambda] =!= 0,
  Print["FAIL: M4 Gaussian Z^2 integral"]; Exit[1]
];
If[FullSimplify[matchedOverlap - Sqrt[2]/2] =!= 0,
  Print["FAIL: M4 matched I_WZ"]; Exit[1]
];
If[FullSimplify[xi*matchedOverlap/matchedOverlap - xi] =!= 0,
  Print["FAIL: M4 H=Z matched gauge"]; Exit[1]
];
If[FullSimplify[xi*matchedOverlap/matchedNorm - xi/Sqrt[2]] =!= 0,
  Print["FAIL: M4 H=1 matched gauge"]; Exit[1]
];
If[FullSimplify[(mu0*matchedDeltaSource/matchedOverlap)/(mu0/zArea) - Sqrt[2]] =!= 0,
  Print["FAIL: M4 delta source ratio"]; Exit[1]
];
If[FullSimplify[mu0*matchedDistributedSource/matchedOverlap - mu0/zArea] =!= 0,
  Print["FAIL: M4 matched distributed source"]; Exit[1]
];
Print["PASS: M4 Gaussian matched-kernel concrete values"];

xi4UnweightedReg[r_] := xi*Sqrt[Pi]*lambda/(2*r);
m5Residual =
  FullSimplify[Limit[xi4UnweightedReg[R], R -> Infinity,
      Assumptions -> lambda > 0 && xi > 0] - 0];
Print["M5 regulator limit residual = ", fmt[m5Residual]];
If[FullSimplify[Limit[xi4UnweightedReg[R], R -> Infinity,
      Assumptions -> lambda > 0 && xi > 0] - 0] =!= 0,
  Print["FAIL: M5 reduction-first H=1 regulator limit"]; Exit[1]
];
Print["PASS: M5 reduction-first H=1 regulator limit"];

xi4General[zAreaArg_, hAreaArg_] := xi*zAreaArg/hAreaArg;
m6Residual = FullSimplify[xi4General[zIntSym, zIntSym] - xi];
Print["M6 H=Z reduction residual = ", fmt[m6Residual]];
If[FullSimplify[xi4General[zIntSym, zIntSym] - xi] =!= 0,
  Print["FAIL: M6 reduction-first H=Z identity"]; Exit[1]
];
Print["PASS: M6 reduction-first H=Z identity"];

lorentzMatchedSourceOverlap =
  Integrate[lorentzWeight*(localizedGaussian/zArea), {w, -Infinity, Infinity},
    Assumptions -> lambda > 0 && sigma > 0];
m7GaugeResiduals =
  N[(xi*lorentzOverlap/lorentzGaugeWeight - xi) /.
    {{lambda -> 1, sigma -> 1/2, xi -> 7/5},
      {lambda -> 1, sigma -> 2, xi -> 7/5}}, 30];
m7SourceResiduals =
  N[(mu0*lorentzMatchedSourceOverlap/lorentzOverlap - mu0/zArea) /.
    {{lambda -> 1, sigma -> 1/2, mu0 -> 11/10},
      {lambda -> 1, sigma -> 2, mu0 -> 11/10}}, 30];
Print["M7 Lorentzian numeric H=Z residuals = ", fmt[m7GaugeResiduals]];
Print["M7 Lorentzian numeric matched-source residuals = ", fmt[m7SourceResiduals]];
If[Chop[m7GaugeResiduals] =!= {0, 0},
  Print["FAIL: M7 Lorentzian numeric H=Z gauge invariance"]; Exit[1]
];
If[Chop[m7SourceResiduals] =!= {0, 0},
  Print["FAIL: M7 Lorentzian numeric matched source"]; Exit[1]
];
Print["PASS: M7 Lorentzian numerical independent-profile checks"];

Print["STATUS: PASS"];
Exit[0];
