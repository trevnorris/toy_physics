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

zeroQ[res_] := If[ListQ[res],
  And @@ (TrueQ[# === 0] & /@ Flatten[res]),
  TrueQ[res === 0]
];

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[zeroQ[res], pass[name], fail[name, res]];
];

banner["STAGE 189 — TRANSFER-SHAPE / OUTGOING-PREFACTOR COMPILER"];

Clear[
  bStar, etaStar, xR, xN, xEta, s, theta1, xi1, sigmaEta,
  r0Obs, n0Obs, lambdaFront, zW, rho, omegaW2, epsW, epsEta,
  gravII, cLightII, aII, csII, chi0, deltaU,
  epsAmp, lambdaAxis, zetaZ, omegaWDrift, chi1, epsW1, deltaU1,
  w, d0, d2, d4, n0, n2, n4, kbl,
  eps, lamA, d1, n1, aa, cs, grav, cLight, mhat0
];

$Assumptions =
  Element[
    {
      bStar, etaStar, xR, xN, xEta, s, theta1, xi1, sigmaEta,
      r0Obs, n0Obs, lambdaFront, zW, rho, omegaW2, epsW, epsEta,
      gravII, cLightII, aII, csII, chi0, deltaU,
      epsAmp, lambdaAxis, zetaZ, omegaWDrift, chi1, epsW1, deltaU1,
      w, d0, d2, d4, n0, n2, n4, kbl,
      eps, lamA, d1, n1, aa, cs, grav, cLight, mhat0
    },
    Reals
  ] &&
  bStar > 0 && 0 < etaStar < 1 && r0Obs > 0 && n0Obs > 0 &&
  lambdaFront > 0 && zW > 0 && rho > -1 && omegaW2 > 0 &&
  epsW < 1 && epsEta < 1 && gravII > 0 && cLightII > 0 &&
  aII > 0 && csII > 0 && chi0 > 0 && deltaU > 0 &&
  lambdaAxis != 0 &&
  d0 != 0 && n0 != 0 && kbl != 0 && aa > 0 && cs > 0 &&
  grav > 0 && cLight > 0 && mhat0 > 0 && lamA != 0;

banner["I. Observable packet -> transfer-shape packet"];

obsVars = {xR, xN, xEta};
transferLogFunctions = {
  xN - bStar*xR,
  Log[(1 - etaStar*Exp[xEta])/(1 - etaStar)],
  Log[(1 - etaStar*Exp[xEta])/(1 - etaStar)] - (xN - bStar*xR)
};
zeroObs = Thread[obsVars -> {0, 0, 0}];
cObsToTrfNative = Table[
  FullSimplify[D[transferLogFunctions[[i]], obsVars[[j]]] /. zeroObs, Assumptions -> $Assumptions],
  {i, 1, 3}, {j, 1, 3}
];
cObsToTrfExpected = {
  {-bStar, 1, 0},
  {0, 0, -etaStar/(1 - etaStar)},
  {bStar, -1, -etaStar/(1 - etaStar)}
};

Print["C_obs->trf from finite-log Jacobian = ", fmt[cObsToTrfNative]];
expectZero["C_obs->trf * Delta_obs^(1) - Delta_trf^(1)", cObsToTrfNative - cObsToTrfExpected];
expectZero["rank(C_obs->trf) - 2", MatrixRank[cObsToTrfNative] - 2];

rTrPath = r0Obs*Exp[s*theta1];
nStarPath = n0Obs*Exp[s*(xi1 + bStar*theta1)];
etaPath = etaStar*Exp[s*sigmaEta];
tSqPath = nStarPath/rTrPath^bStar;
xiFromObservable = FullSimplify[D[Log[tSqPath], s] /. s -> 0, Assumptions -> $Assumptions];
oneMinusSlope = FullSimplify[D[Log[1 - etaPath], s] /. s -> 0, Assumptions -> $Assumptions];
rTargetPath = lambdaFront*(1 - etaPath)/tSqPath;
rOneSlope = FullSimplify[D[Log[rTargetPath], s] /. s -> 0, Assumptions -> $Assumptions];

expectZero["dln T^2 - Xi_1", xiFromObservable - xi1];
expectZero["dln(1-epseta) - (R_1 + Xi_1)", oneMinusSlope - (rOneSlope + xiFromObservable)];
expectZero[
  "compatibility: dln R_target + dln T^2 - dln(1-epseta)",
  rOneSlope + xiFromObservable - oneMinusSlope
];

banner["II. One-port continuum transfer-shape identity"];

lambda0Value = 27*Pi^2*gravII*csII^5/(20*aII^5*cLightII^5);
t2OnePort = zW*(1 + rho)^2/(omegaW2*(1 - epsW)^2);
rTargetDefinition = FullSimplify[lambdaFront*(1 - epsEta)/t2OnePort, Assumptions -> $Assumptions];
Print["T_A^2 (direct continuum form) = ", fmt[t2OnePort]];
Print["R_target selected-branch definition := Lambda_0 (1-epseta) / T_A^2 = ", fmt[rTargetDefinition]];
Print["Lambda_0 value = ", fmt[lambda0Value]];

epsSplit = epsW*(1 - 2*deltaU/(11*(1 + deltaU)));
t2Coherent = t2OnePort /. {rho -> chi0, epsW -> epsSplit};
expectZero[
  "coherent local D/N specialization",
  t2Coherent - zW*(1 + chi0)^2/(omegaW2*(1 - epsSplit)^2)
];
eps1Bridge = FullSimplify[D[epsSplit, epsW]*epsW1 + D[epsSplit, deltaU]*deltaU1, Assumptions -> $Assumptions];
xi1Closed = FullSimplify[
  zetaZ - omegaWDrift + 2*chi1/(1 + chi0) + 2*eps1Bridge/(1 - epsSplit),
  Assumptions -> $Assumptions
];
epsSplitPert = (epsW + s*epsAmp*lambdaAxis*epsW1)*
  (1 - (2/11)*(deltaU + s*epsAmp*lambdaAxis*deltaU1)/(1 + deltaU + s*epsAmp*lambdaAxis*deltaU1));
t2CoherentPert = (zW*(1 + s*epsAmp*lambdaAxis*zetaZ))*
  (1 + chi0 + s*epsAmp*lambdaAxis*chi1)^2/
  ((omegaW2*(1 + s*epsAmp*lambdaAxis*omegaWDrift))*(1 - epsSplitPert)^2);
directSlope = FullSimplify[(D[Log[t2CoherentPert/t2Coherent], s] /. s -> 0), Assumptions -> $Assumptions];
expectZero[
  "direct-slope bridge dln T_A^2 - epsilon lambda_A Xi_1",
  directSlope - epsAmp*lambdaAxis*xi1Closed
];

banner["III. Exact isotropic grouped response / prefactor compiler"];

den[w_] := d0 + d2*w^2 + d4*w^4;
num[w_] := n0 + n2*w^2 + n4*w^4;
response[w_] := d0/den[w];
prefactor[w_] := d0*num[w]/den[w]^2;

u2FromDerivative = FullSimplify[(D[response[w], {w, 2}]/2) /. w -> 0, Assumptions -> $Assumptions];
u4FromDerivative = FullSimplify[(D[response[w], {w, 4}]/24) /. w -> 0, Assumptions -> $Assumptions];
u2Expected = -d2/d0;
u4Expected = (d2^2 - d0*d4)/d0^2;

expectZero["Y(omega) - (1 + u2 w^2 + u4 w^4)", {
  u2FromDerivative - u2Expected,
  u4FromDerivative - u4Expected
}];

p0Coeff = FullSimplify[prefactor[w] /. w -> 0, Assumptions -> $Assumptions];
p2Coeff = FullSimplify[(D[prefactor[w], {w, 2}]/2) /. w -> 0, Assumptions -> $Assumptions];
p4Coeff = FullSimplify[(D[prefactor[w], {w, 4}]/24) /. w -> 0, Assumptions -> $Assumptions];
p0Expected = n0/d0;
p2Expected = (d0*n2 - 2*d2*n0)/d0^2;
p4Expected = (d0^2*n4 - 2*d0*(d2*n2 + d4*n0) + 3*d2^2*n0)/d0^3;

expectZero["Pref(omega) - (P0 + P2 w^2 + P4 w^4)", {
  p0Coeff - p0Expected,
  p2Coeff - p2Expected,
  p4Coeff - p4Expected
}];

tEff2 = n0/kbl;
expectZero["P0 - (K_bl/D0) T_eff^2", p0Coeff - (kbl/d0)*tEff2];

banner["IV. Weak-axisymmetric prefactor slope"];

pLane[eps_] := (n0 + eps*lamA*n1)/(d0 + eps*lamA*d1);
p1FromDerivative = FullSimplify[(D[pLane[eps], eps] /. eps -> 0)/lamA, Assumptions -> $Assumptions];
p1Expected = (n1*d0 - n0*d1)/d0^2;
logSlope = FullSimplify[D[Log[pLane[eps]/p0Expected], eps] /. eps -> 0, Assumptions -> $Assumptions];

expectZero["P_A - (P0 + eps lambda P1)", p1FromDerivative - p1Expected];
expectZero["log(P_A/P0) - eps lambda (P1/P0)", logSlope - lamA*(p1Expected/p0Expected)];
expectZero["P1/P0 - (N1/N0 - D1/D0)", p1Expected/p0Expected - (n1/n0 - d1/d0)];

banner["V. Compact outgoing l=2 fingerprint compiler"];

aOut = aa^2/(9*cs^2);
bOut = 4*aa^4/(81*cs^4);
g5Out = aa^5/(27*cs^5);
prefPoly = p0Expected + p2Expected*w^2 + p4Expected*w^4;
yhatOut = 1 + aOut*w^2 + bOut*w^4 + I*g5Out*w^5;
outTrunc = Normal[Series[Expand[prefPoly*yhatOut], {w, 0, 5}]];
k0Coeff = Coefficient[outTrunc, w, 0];
k2Coeff = Coefficient[outTrunc, w, 2];
k4Coeff = Coefficient[outTrunc, w, 4];
gamma5Coeff = Coefficient[outTrunc, w, 5]/I;

k0Expected = p0Expected;
k2Expected = p2Expected + aOut*p0Expected;
k4Expected = p4Expected + aOut*p2Expected + bOut*p0Expected;
gamma5Expected = g5Out*p0Expected;

expectZero[
  "outgoing branch expansion - (K0 + K2 w^2 + K4 w^4 + i Gamma5 w^5)",
  outTrunc - (k0Expected + k2Expected*w^2 + k4Expected*w^4 + I*gamma5Expected*w^5)
];

banner["VI. Normalization equivalence and constant-prefactor branch"];

p0Target = 54*grav*cs^5/(5*aa^5*cLight^5*mhat0^2);
expectZero[
  "mhat0^2 Gamma5 - 2G/(5c^5) at P0_target",
  (mhat0^2*gamma5Expected - 2*grav/(5*cLight^5)) /. p0Expected -> p0Target
];
expectZero["Gamma5 - a^5 P0/(27 c_s^5)", gamma5Expected - aa^5*p0Expected/(27*cs^5)];

n2Const = FullSimplify[n2 /. First[Solve[p2Expected == 0, n2]], Assumptions -> $Assumptions];
n4Const = FullSimplify[n4 /. First[Solve[(p4Expected /. n2 -> n2Const) == 0, n4]], Assumptions -> $Assumptions];

expectZero["P2 on constant-prefactor branch", p2Expected /. n2 -> n2Const];
expectZero["P4 on constant-prefactor branch", p4Expected /. {n2 -> n2Const, n4 -> n4Const}];
expectZero["K2 on constant-prefactor branch - A P0", (k2Expected /. n2 -> n2Const) - aOut*p0Expected];
expectZero[
  "K4 on constant-prefactor branch - B P0",
  (k4Expected /. {n2 -> n2Const, n4 -> n4Const}) - bOut*p0Expected
];

banner["STAGE 189 MATHEMATICA LEDGER"];
Print["The transfer packet is recovered from finite log perturbations and the"];
Print["prefactor/outgoing compilers are recovered by derivative, coefficient, and"];
Print["Solve routes independent of the SymPy script choreography."];

Exit[0];
