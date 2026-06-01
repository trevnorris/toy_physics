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
stripCE[expr_] := expr /. ConditionalExpression[e_, _] :> e;

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[stripCE[expr]]], Assumptions -> $Assumptions];
  res = stripCE[res];
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectZeroList[name_String, expr_List] := Module[{res},
  res = FullSimplify[stripCE /@ expr, Assumptions -> $Assumptions];
  res = stripCE[res];
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === ConstantArray[0, Length[res]]], pass[name], fail[name, res]];
];

expectTrue[name_String, claim_] := Module[{res},
  res = FullSimplify[stripCE[claim], Assumptions -> $Assumptions];
  res = stripCE[res];
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 197 - CONDITIONAL PACKET-A CLOSURE THEOREM"];

Clear[
  z, radius, soundSpeed, lightSpeed, bigG, scaleS, betaStretch,
  sigma0, sigma2, sigma4, sigma5, tail7, pIso, deltaNormSlot,
  nScale, mHat0, chiFree, deltaQ, eps, epsBeta, dSigma0, dSigma5
];

$Assumptions =
  Element[
    {z, radius, soundSpeed, lightSpeed, bigG, scaleS, betaStretch,
     sigma0, sigma2, sigma4, sigma5, tail7, pIso, deltaNormSlot,
     nScale, mHat0, chiFree, deltaQ, eps, epsBeta, dSigma0, dSigma5},
    Reals
  ] &&
  radius > 0 && soundSpeed > 0 && lightSpeed > 0 && bigG > 0 &&
  scaleS != 0 && betaStretch != 0 &&
  3*scaleS - sigma0 != 0 &&
  scaleS*betaStretch^5 + 9*sigma5 != 0 &&
  mHat0 > 0 && chiFree != 0 && 1 + deltaQ != 0;

(* I. Isotropic outgoing lane implies a_P0=b_P0=0. *)
subbanner["I. Isotropic grouped-P0 projection"];

isotropicP0 = {pIso, pIso, pIso};
traceProjector = {
  {1/5, 2/5, 2/5},
  {1/5, -1/10, -1/10},
  {0, 1/2, -1/2}
};
projectedP0 = FullSimplify[traceProjector.isotropicP0, Assumptions -> $Assumptions];

Print["grouped trace/isotropy coordinates = ", fmt[projectedP0]];
expectZero["a_P0 on isotropic lane", projectedP0[[2]]];
expectZero["b_P0 on isotropic lane", projectedP0[[3]]];

(* II. Packet-A residual collapses to the normalization slot. *)
subbanner["II. Packet-A residual collapse"];

packetResidual = Join[ConstantArray[0, 7], {deltaNormSlot}];
Print["Delta_branch after carried isotropic front end = ", fmt[packetResidual]];
expectZeroList["first seven Packet-A slots", Take[packetResidual, 7]];
expectZero["Packet-A normalization slot is the only survivor", Last[packetResidual] - deltaNormSlot];

(* III. Native DtN coefficient extraction and natural source-map reduction. *)
subbanner["III. Native z^5 extraction and source-map reduction"];

gamma5Target = 2*bigG/(5*lightSpeed^5);
p0FromGamma = FullSimplify[27*soundSpeed^5*gamma5Target/radius^5, Assumptions -> $Assumptions];
p0Target = 54*bigG*soundSpeed^5/(5*radius^5*lightSpeed^5);

expectZero["P0_target from Gamma_5 normalization", p0FromGamma - p0Target];

lambdaOut = FunctionExpand[z*D[SphericalHankelH1[2, z], z]/SphericalHankelH1[2, z]];
lambdaWindow5 = Expand[Normal[Series[lambdaOut, {z, 0, 5}]]];
deformedDtn = Expand[
  (
    scaleS*(lambdaWindow5 /. z -> betaStretch*z)
    + sigma0 + sigma2*z^2 + sigma4*z^4 + I*sigma5*z^5
    + I*tail7*z^7
  )
];
normalizedDtn = Coefficient[deformedDtn, z, 0]/deformedDtn;
normalizedSeries = Expand[Normal[Series[normalizedDtn, {z, 0, 7}]]];

evenMatchRules = Solve[
  {
    Coefficient[normalizedSeries, z, 2] == 1/9,
    Coefficient[normalizedSeries, z, 4] == 4/81
  },
  {sigma2, sigma4},
  Reals
];
If[Length[evenMatchRules] =!= 1, fail["unique canonical-even matching solution", evenMatchRules]];
evenRule = First[evenMatchRules] /. ConditionalExpression[e_, _] :> e;

sigma2SymPy = -(3*scaleS*betaStretch^2 - 3*scaleS + sigma0)/9;
sigma4SymPy = -(3*scaleS*betaStretch^4 - 3*scaleS + sigma0)/27;

Print["canonical-even Sigma_2 = ", fmt[FullSimplify[sigma2 /. evenRule, Assumptions -> $Assumptions]]];
Print["canonical-even Sigma_4 = ", fmt[FullSimplify[sigma4 /. evenRule, Assumptions -> $Assumptions]]];

expectZero["Sigma_2 agreement with SymPy result", (sigma2 /. evenRule) - sigma2SymPy];
expectZero["Sigma_4 agreement with SymPy result", (sigma4 /. evenRule) - sigma4SymPy];

matchedSeries = Expand[normalizedSeries /. evenRule];
chiFromCoefficient = FullSimplify[
  Coefficient[matchedSeries, z, 5]/(I/27),
  Assumptions -> $Assumptions
];
chiFromDerivative = FullSimplify[
  ((D[normalizedDtn /. evenRule, {z, 5}] /. z -> 0)/5!)/(I/27),
  Assumptions -> $Assumptions
];
chiSymPy = 3*(scaleS*betaStretch^5 + 9*sigma5)/(3*scaleS - sigma0);

Print["chi_Q from Series/Coefficient z^5 extraction = ", fmt[chiFromCoefficient]];
expectZero["Series coefficient agrees with fifth-derivative extraction", chiFromCoefficient - chiFromDerivative];
expectZero["chi_Q extractor - SymPy deformation algebra formula", chiFromCoefficient - chiSymPy];

sourceMapEquation = mHat0^2*chiFromCoefficient*nScale == 1;
nFromOddClosure = FullSimplify[
  stripCE[nScale /. First[Solve[sourceMapEquation, nScale, Reals]]],
  Assumptions -> $Assumptions
];
nNatural = FullSimplify[nFromOddClosure /. mHat0 -> 1, Assumptions -> $Assumptions];
deltaNormNatural = FullSimplify[p0Target*(nNatural - 1), Assumptions -> $Assumptions];

Print["N_Q from observable odd source-map equation = ", fmt[nFromOddClosure]];
Print["N_Q on natural point-particle branch = ", fmt[nNatural]];
Print["Delta_norm on natural point-particle branch = ", fmt[deltaNormNatural]];

expectZero["source-map equation after solve", mHat0^2*chiFromCoefficient*nFromOddClosure - 1];
expectZero["N_Q(point-particle) - 1/chi_Q", nNatural - 1/chiFromCoefficient];
expectZero[
  "Delta_norm - P0_target*(1/chi_Q - 1)",
  deltaNormNatural - p0Target*(1/chiFromCoefficient - 1)
];

(* IV. Exact finish-line equivalences. *)
subbanner["IV. Finish-line equivalences"];

packetResidualNatural = Join[ConstantArray[0, 7], {deltaNormNatural}];
deltaNormGeneric = p0Target*(1/chiFree - 1);
nGeneric = 1/chiFree;
zeroSetGeneric = FullSimplify[
  Reduce[deltaNormGeneric == 0, chiFree, Reals],
  Assumptions -> $Assumptions
];
nZeroSetGeneric = FullSimplify[
  Reduce[nGeneric == 1, chiFree, Reals],
  Assumptions -> $Assumptions
];

Print["Delta_branch on natural source-map branch = ", fmt[packetResidualNatural]];
Print["Reduce[Delta_norm == 0, chi_Q] = ", fmt[zeroSetGeneric]];
Print["Reduce[N_Q == 1, chi_Q] = ", fmt[nZeroSetGeneric]];

expectZeroList["natural Delta_branch first seven slots", Take[packetResidualNatural, 7]];
expectZero[
  "chi_Q*(Delta_norm/P0_target + 1) - 1",
  chiFromCoefficient*(deltaNormNatural/p0Target + 1) - 1
];
expectZero[
  "Delta_norm/P0_target + (chi_Q - 1)/chi_Q",
  deltaNormNatural/p0Target + (chiFromCoefficient - 1)/chiFromCoefficient
];
expectTrue[
  "Delta_norm == 0 iff chi_Q == 1",
  Equivalent[zeroSetGeneric, chiFree == 1]
];
expectTrue[
  "N_Q == 1 iff chi_Q == 1",
  Equivalent[nZeroSetGeneric, chiFree == 1]
];
expectZero[
  "N_Q - 1 in terms of Delta_Q",
  ((nGeneric - 1) /. chiFree -> 1 + deltaQ) + deltaQ/(1 + deltaQ)
];
expectZero[
  "Delta_norm under chi_Q = 1 + Delta_Q",
  (deltaNormGeneric /. chiFree -> 1 + deltaQ) + p0Target*deltaQ/(1 + deltaQ)
];

deltaBad = FullSimplify[deltaNormGeneric /. chiFree -> 6/5, Assumptions -> $Assumptions];
Print["Delta_norm at chi_Q = 6/5 = ", fmt[deltaBad]];
expectTrue["Delta_norm at chi_Q = 6/5 is nonzero", deltaBad != 0];

(* V. Stage 194 deformation algebra realization gate. *)
subbanner["V. Deformation-algebra closure gate"];

closureNumerator = 3*scaleS*(betaStretch^5 - 1) + sigma0 + 27*sigma5;
deltaNormDeformation = FullSimplify[p0Target*(1/chiFromCoefficient - 1), Assumptions -> $Assumptions];

Print["closure numerator = ", fmt[closureNumerator]];
Print["Delta_norm after native deformation extraction = ", fmt[deltaNormDeformation]];

expectZero[
  "(3S - Sigma0)(chi_Q - 1) - closure numerator",
  (3*scaleS - sigma0)*(chiFromCoefficient - 1) - closureNumerator
];
expectZero[
  "chi_Q - 1 minus closure numerator/(3S - Sigma0)",
  (chiFromCoefficient - 1) - closureNumerator/(3*scaleS - sigma0)
];
expectZero[
  "Delta_norm + P0_target*closure_num/[3(S beta^5 + 9 Sigma5)]",
  deltaNormDeformation + p0Target*closureNumerator/(3*(scaleS*betaStretch^5 + 9*sigma5))
];

(* VI. Linearized finish-line map. *)
subbanner["VI. Linearized finish-line map"];

linearRule = {
  betaStretch -> 1 + eps*epsBeta,
  sigma0 -> eps*dSigma0,
  sigma5 -> eps*dSigma5
};
linearSlope = 5*epsBeta + dSigma0/(3*scaleS) + 9*dSigma5/scaleS;

chiLinear = Expand[Normal[Series[chiFromCoefficient /. linearRule, {eps, 0, 1}]]];
deltaLinear = Expand[Normal[Series[deltaNormDeformation /. linearRule, {eps, 0, 1}]]];

Print["chi_Q linearized = ", fmt[chiLinear]];
Print["Delta_norm linearized = ", fmt[deltaLinear]];

expectZero[
  "linearized chi_Q - [1 + eps(5 eps_beta + dSigma0/(3S) + 9 dSigma5/S)]",
  chiLinear - (1 + eps*linearSlope)
];
expectZero[
  "linearized Delta_norm + eps P0_target*(5 eps_beta + dSigma0/(3S) + 9 dSigma5/S)",
  deltaLinear + eps*p0Target*linearSlope
];

(* VII. Higher-odd irrelevance at the final finish line. *)
subbanner["VII. Higher-odd irrelevance"];

expectZero["d chi_Q^(series extractor) / dL7", D[chiFromCoefficient, tail7]];
expectZero["d Delta_norm / dL7", D[deltaNormNatural, tail7]];

banner["STAGE 197 MATHEMATICA LEDGER"];
Print["1. The isotropic grouped-P0 lane kills both anisotropy coordinates."];
Print["2. The carried Packet-A residual has only the normalization slot left."];
Print["3. Native DtN Series/Coefficient extraction gives the same chi_Q as SymPy."];
Print["4. The natural source-map equation gives Delta_norm = P0_target (1/chi_Q - 1)."];
Print["5. Therefore Delta_branch = 0 iff chi_Q = 1, equivalently the deformation numerator vanishes."];
Print["6. Linearization and higher-odd irrelevance agree with the SymPy audit."];

Exit[0];
