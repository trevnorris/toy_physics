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

banner["STAGE 195 - SOURCE-MAP REDUCTION OF THE CANONICAL OUTGOING BRANCH"];

Clear[radius, soundSpeed, lightSpeed, bigG, p0Bar, chiQ, mhat0, nQ, deltaQ];
Clear[x, scaleS, betaStretch, sigma0, sigma2, sigma4, sigma5];
Clear[epsBeta, dSigma0, dSigma5, tau];

$Assumptions =
  Element[
    {radius, soundSpeed, lightSpeed, bigG, p0Bar, chiQ, mhat0, nQ, deltaQ,
     x, scaleS, betaStretch, sigma0, sigma2, sigma4, sigma5,
     epsBeta, dSigma0, dSigma5, tau},
    Reals
  ] &&
  radius > 0 && soundSpeed > 0 && lightSpeed > 0 && bigG > 0 &&
  p0Bar > 0 && chiQ > 0 && mhat0 > 0 &&
  scaleS != 0 && betaStretch != 0 &&
  3*scaleS - sigma0 != 0 && scaleS*betaStretch^5 + 9*sigma5 != 0 &&
  1 + deltaQ != 0;

(* I. Carry-forward isotropic retarded invariant tuple. *)
subbanner["I. Carry-forward isotropic retarded invariant tuple"];

omegaQ = 3*soundSpeed/(2*radius);
p0Target = 54*bigG*soundSpeed^5/(5*radius^5*lightSpeed^5);
gamma5Target = 2*bigG/(5*lightSpeed^5);

gammaByGeometry = chiQ*radius^5*p0Bar/(27*soundSpeed^5);
gammaByPole = FullSimplify[chiQ*9*p0Bar/(32*omegaQ^5), Assumptions -> $Assumptions];
gammaRatio = FullSimplify[
  (gammaByGeometry/gamma5Target) /. p0Bar -> nQ*p0Target,
  Assumptions -> $Assumptions
];

Print["Omega_Q = ", fmt[omegaQ]];
Print["P0_target = ", fmt[p0Target]];
Print["Gamma5_target = ", fmt[gamma5Target]];
Print["Gamma5 = ", fmt[gammaByGeometry]];

expectZero["Gamma5 - chi_Q*9*P0/(32 Omega_Q^5)", gammaByGeometry - gammaByPole];
expectZero["Gamma5/Gamma5_target - chi_Q*N_Q", gammaRatio - chiQ*nQ];

(* II. Observable odd closure and Packet-A collapse. *)
subbanner["II. Observable odd closure and Packet-A collapse"];

oddObservableResidual = FullSimplify[
  (mhat0^2*gammaByGeometry - gamma5Target) /. p0Bar -> nQ*p0Target,
  Assumptions -> $Assumptions
];
oddClosureEquation = mhat0^2*chiQ*nQ == 1;
nQFromOdd = FullSimplify[
  stripCE[nQ /. First[Solve[oddClosureEquation, nQ, Reals]]],
  Assumptions -> $Assumptions
];
deltaAfterOdd = FullSimplify[
  (mhat0^2*nQ*p0Target - p0Target) /. nQ -> nQFromOdd,
  Assumptions -> $Assumptions
];

Print["odd observable residual = ", fmt[oddObservableResidual]];
Print["N_Q from odd closure = ", fmt[nQFromOdd]];
Print["Delta_norm after imposing odd closure = ", fmt[deltaAfterOdd]];

expectZero[
  "observable odd condition factorizes as Gamma5_target*(mhat0^2 chi_Q N_Q - 1)",
  oddObservableResidual - gamma5Target*(mhat0^2*chiQ*nQ - 1)
];
expectZero[
  "Delta_norm from odd closure - P0_target*(1/chi_Q - 1)",
  deltaAfterOdd - p0Target*(1/chiQ - 1)
];
expectZero[
  "Delta_norm in terms of Delta_Q",
  (deltaAfterOdd /. chiQ -> 1 + deltaQ) + p0Target*deltaQ/(1 + deltaQ)
];

(* III. Natural source-map reduction. *)
subbanner["III. Natural source-map reduction"];

nQNatural = FullSimplify[nQFromOdd /. mhat0 -> 1, Assumptions -> $Assumptions];
deltaNatural = FullSimplify[
  (p0Target*nQ - p0Target) /. nQ -> nQNatural,
  Assumptions -> $Assumptions
];

Print["N_Q on the natural point-particle source-map branch = ", fmt[nQNatural]];
Print["Delta_norm on the natural point-particle source-map branch = ", fmt[deltaNatural]];

expectZero["N_Q(point-particle) - 1/chi_Q", nQNatural - 1/chiQ];
expectZero[
  "Delta_norm(point-particle) - P0_target*(1/chi_Q - 1)",
  deltaNatural - p0Target*(1/chiQ - 1)
];
expectZero[
  "N_Q(point-particle)-1 in terms of Delta_Q",
  ((nQNatural - 1) /. chiQ -> 1 + deltaQ) + deltaQ/(1 + deltaQ)
];

(* IV. DtN deformation algebra rederived from the outgoing l=2 operator. *)
subbanner["IV. Stage 194 deformation algebra after source-map reduction"];

lambdaOut = FunctionExpand[x*D[SphericalHankelH1[2, x], x]/SphericalHankelH1[2, x]];
lambdaWindow = Expand[Normal[Series[lambdaOut, {x, 0, 5}]]];
deformedDtn = Expand[
  scaleS*(lambdaWindow /. x -> betaStretch*x)
  + sigma0 + sigma2*x^2 + sigma4*x^4 + I*sigma5*x^5
];
normalizedDtn = Expand[
  Normal[Series[Coefficient[deformedDtn, x, 0]/deformedDtn, {x, 0, 5}]]
];

evenRules = Solve[
  {
    Coefficient[normalizedDtn, x, 2] == 1/9,
    Coefficient[normalizedDtn, x, 4] == 4/81
  },
  {sigma2, sigma4},
  Reals
];
If[Length[evenRules] =!= 1, fail["unique canonical-even DtN deformation", evenRules]];

chiFromDtn = FullSimplify[
  Coefficient[normalizedDtn /. First[evenRules], x, 5]/(I/27),
  Assumptions -> $Assumptions
];
nQFromDtn = FullSimplify[1/chiFromDtn, Assumptions -> $Assumptions];
nQDeformationTarget = (3*scaleS - sigma0)/(3*(scaleS*betaStretch^5 + 9*sigma5));
deltaDeformation = FullSimplify[p0Target*(nQFromDtn - 1), Assumptions -> $Assumptions];
deltaDeformationTarget = FullSimplify[
  -p0Target*(3*scaleS*(betaStretch^5 - 1) + sigma0 + 27*sigma5)/
    (3*(scaleS*betaStretch^5 + 9*sigma5)),
  Assumptions -> $Assumptions
];

Print["chi_Q from native DtN deformation = ", fmt[chiFromDtn]];
Print["N_Q after natural source-map reduction = ", fmt[nQFromDtn]];
Print["Delta_norm(point-particle) after inserting DtN deformation algebra = ",
  fmt[deltaDeformation]];

expectZero["N_Q deformation law", nQFromDtn - nQDeformationTarget];
expectZero["Delta_norm deformation law", deltaDeformation - deltaDeformationTarget];

nQLinear = FullSimplify[
  Coefficient[
    Normal[
      Series[
        (nQFromDtn - 1) /. {
          betaStretch -> 1 + tau*epsBeta,
          sigma0 -> tau*dSigma0,
          sigma5 -> tau*dSigma5
        },
        {tau, 0, 1}
      ]
    ],
    tau,
    1
  ],
  Assumptions -> $Assumptions
];
deltaLinear = FullSimplify[p0Target*nQLinear, Assumptions -> $Assumptions];

Print["linearized N_Q - 1 = ", fmt[Expand[nQLinear]]];
Print["linearized Delta_norm(point-particle) = ", fmt[Expand[deltaLinear]]];

expectZero[
  "linearized N_Q - 1",
  nQLinear + 5*epsBeta + dSigma0/(3*scaleS) + 9*dSigma5/scaleS
];
expectZero[
  "linearized Delta_norm(point-particle)",
  deltaLinear + p0Target*(5*epsBeta + dSigma0/(3*scaleS) + 9*dSigma5/scaleS)
];

(* V. Canonical compact outgoing branch. *)
subbanner["V. Canonical compact outgoing branch"];

canonicalRule = {betaStretch -> 1, sigma0 -> 0, sigma5 -> 0};

expectZero["chi_Q(canonical) - 1", (chiFromDtn /. canonicalRule) - 1];
expectZero["N_Q(canonical) - 1", (nQFromDtn /. canonicalRule) - 1];
expectZero["Delta_norm(canonical)", deltaDeformation /. canonicalRule];
expectZero[
  "P0(canonical source-map-reduced) - P0_target",
  p0Target*(nQFromDtn /. canonicalRule) - p0Target
];

banner["STAGE 195 MATHEMATICA LEDGER"];
Print["1. The retarded invariant tuple gives Gamma5/Gamma5_target = chi_Q N_Q."];
Print["2. The observable odd condition factorizes as mhat0^2 chi_Q N_Q = 1."];
Print["3. Imposing that condition collapses Delta_norm to P0_target (1/chi_Q - 1)."];
Print["4. The native DtN deformation route gives the same source-map-reduced laws."];
Print["5. The canonical compact outgoing branch closes with N_Q = 1 and Delta_norm = 0."];

Exit[0];
