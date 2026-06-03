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

stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

cleanScalar[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

cleanExpr[expr_] := If[
  ListQ[expr],
  Map[cleanScalar, expr, {ArrayDepth[expr]}],
  cleanScalar[expr]
];

allZeroQ[expr_] := If[
  ListQ[expr],
  And @@ Flatten[Map[TrueQ[# === 0] &, expr, {ArrayDepth[expr]}]],
  TrueQ[expr === 0]
];

prettyArray[arr_] := If[VectorQ[arr], MatrixForm[{arr}], MatrixForm[arr]];

expectZero[name_String, expr_] := Module[{res},
  res = cleanExpr[expr];
  If[ListQ[res],
    Print[name, " ="];
    Print[prettyArray[res]];
    If[allZeroQ[res], pass[name], fail[name, res]],
    Print[name, " = ", fmt[res]];
    If[allZeroQ[res], pass[name], fail[name, res]]
  ];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 50] - N[target, 50]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 245 -- NON-RIGID MOUTH/DRESSING PACKET AND U/V DRAIN COMPILER"];

Clear[
  aU, aV, chiUV, fU, uVar, vVar, T2ref, epsref, Lambda0, u1, v1,
  eps, Lam, varrho
];

$Assumptions = (
  Element[{aU, aV, T2ref, epsref, Lambda0}, Reals] &&
  Element[{chiUV, fU, uVar, vVar, u1, v1, eps, Lam, varrho}, Reals] &&
  aU > 0 && aV > 0 && T2ref > 0 && Lambda0 > 0 && 0 < epsref < 1
);

deltaUV = aU aV - chiUV^2;

subbanner["M1. Stationary packet"];

freeEnergy = (
  1/2 aU uVar^2 + 1/2 aV vVar^2 - chiUV uVar vVar - fU uVar
);
stationaryEquations = {
  D[freeEnergy, uVar] == 0,
  D[freeEnergy, vVar] == 0
};
solutionRule = First[Solve[stationaryEquations, {uVar, vVar}]];
uSol = cleanScalar[uVar /. solutionRule];
vSol = cleanScalar[vVar /. solutionRule];
hessian = {{aU, -chiUV}, {-chiUV, aV}};

Print["F_nr(U,V) = ", fmt[freeEnergy]];
Print["stationary solution = ", fmt[solutionRule]];
Print["Delta_UV = ", fmt[deltaUV]];
Print["U_sol = ", fmt[uSol]];
Print["V_sol = ", fmt[vSol]];
Print["Hessian = ", MatrixForm[hessian]];

expectZero["M1 U_sol closed form", uSol - aV fU/deltaUV];
expectZero["M1 V_sol closed form", vSol - chiUV fU/deltaUV];
expectZero["M1 V_sol/U_sol ratio", vSol/uSol - chiUV/aV];
expectZero["M1 Hessian determinant", Det[hessian] - deltaUV];
expectZero["M1 U_sol at f_U=0", uSol /. fU -> 0];
expectZero["M1 V_sol at f_U=0", vSol /. fU -> 0];
expectZero["M1 V_sol at chi_UV=0", vSol /. chiUV -> 0];
expectZero["M1 U_sol at chi_UV=0", (uSol /. chiUV -> 0) - fU/aU];

subbanner["M2. Drain form"];

drain = cleanScalar[chiUV uSol vSol];
drainClosed = chiUV^2 aV fU^2/deltaUV^2;
Print["D_UV = ", fmt[drain]];
expectZero["M2 drain closed form", drain - drainClosed];

subbanner["M3. Drain nonnegativity on opposite-sign branch"];

drainOppositeSign = N[
  drainClosed /. {aU -> 5/2, aV -> 3, chiUV -> -19/25, fU -> 33/100},
  30
];
Print["D_UV at opposite-sign point = ", fmt[drainOppositeSign]];
expectTrue["M3 opposite-sign drain > 0", drainOppositeSign > 0];

subbanner["M4. Finite physical compiler from selected-branch identity"];

T2 = T2ref Exp[uSol];
epseta = epsref Exp[vSol];
Rtarget = Lambda0 (1 - epseta)/T2;
Rtargetref = Lambda0 (1 - epsref)/T2ref;
RtargetRatio = cleanScalar[Rtarget/Rtargetref];
RtargetPaper = ((1 - epsref Exp[vSol])/(1 - epsref)) Exp[-uSol];

Print["T^2/T^2_ref = ", fmt[cleanScalar[T2/T2ref]]];
Print["epsilon_eta/eps_ref = ", fmt[cleanScalar[epseta/epsref]]];
Print["R_target/R_ref from identity = ", fmt[RtargetRatio]];

expectZero["M4 R_target ratio from identity", RtargetRatio - RtargetPaper];
expectZero["M4 T^2 compiler", T2/T2ref - Exp[uSol]];
expectZero["M4 epsilon_eta compiler", epseta/epsref - Exp[vSol]];

subbanner["M5. Dependent microscopic correction"];

dependentVector = {0, -vSol, uSol - vSol};
dependentExpected = {
  0,
  -chiUV fU/deltaUV,
  (aV - chiUV) fU/deltaUV
};
Print["dependent correction = ", fmt[dependentVector]];
expectZero["M5 dependent microscopic correction", dependentVector - dependentExpected];

subbanner["M6. First-order packet"];

xiCoeff = Coefficient[Normal[Series[Log[Exp[eps u1]], {eps, 0, 1}]], eps];
oneMinusCoeff = Coefficient[
  Normal[Series[Log[(1 - epsref Exp[eps v1])/(1 - epsref)], {eps, 0, 1}]],
  eps
];
ratioCoeff = Coefficient[
  Normal[
    Series[
      Log[((1 - epsref Exp[eps v1])/(1 - epsref)) Exp[-eps u1]],
      {eps, 0, 1}
    ]
  ],
  eps
];

Xi1 = xiCoeff;
R1plusXi1 = oneMinusCoeff;
R1 = ratioCoeff;

Print["Xi1 = ", fmt[cleanScalar[Xi1]]];
Print["R1 + Xi1 = ", fmt[cleanScalar[R1plusXi1]]];
Print["R1 = ", fmt[cleanScalar[R1]]];

expectZero["M6 Xi1 coefficient", xiCoeff - u1];
expectZero["M6 R1+Xi1 coefficient", oneMinusCoeff - (-epsref v1/(1 - epsref))];
expectZero["M6 R1 coefficient", ratioCoeff - (-u1 - epsref v1/(1 - epsref))];
expectZero["M6 packet consistency", R1 + Xi1 - R1plusXi1];

subbanner["M7. Orbit-side support-blindness with positive control"];

deltaKeta = -vSol;
deltaMu = uSol - vSol;
supportObjects = {
  {"U_sol", uSol},
  {"V_sol", vSol},
  {"epseta", epseta},
  {"Rtarget/Rtargetref", RtargetRatio},
  {"drain", drain},
  {"Delta_Keta", deltaKeta},
  {"Delta_mu", deltaMu}
};

Do[
  expectZero["M7 d/dLam " <> entry[[1]], D[entry[[2]], Lam]];
  expectZero["M7 d/dvarrho " <> entry[[1]], D[entry[[2]], varrho]],
  {entry, supportObjects}
];

controlDerivative = D[aV (fU + Lam)/(aU aV - chiUV^2), Lam];
Print["control d/dLam U_bad = ", fmt[controlDerivative]];
expectTrue["M7 positive control detects support contamination", controlDerivative =!= 0];
expectZero["M7 positive control exact derivative", controlDerivative - aV/(aU aV - chiUV^2)];

subbanner["M8. Session-I readback"];

uObs = 14313458/100000000;
vObs = -3619791/100000000;
epsrefObs = 3/10;
epsetaObs = epsrefObs Exp[vObs];
RratioObs = ((1 - epsrefObs Exp[vObs])/(1 - epsrefObs)) Exp[-uObs];
R1Obs = -uObs - epsrefObs vObs/(1 - epsrefObs);

Print["epsilon_eta(obs) = ", fmt[N[epsetaObs, 20]]];
Print["R_target/R_ref(obs) = ", fmt[N[RratioObs, 20]]];
Print["R1_obs = ", fmt[N[R1Obs, 20]]];

expectApprox["M8 epsilon_eta Session-I", epsetaObs, 14466741/50000000, 5*10^-9];
expectApprox["M8 R_target/R_ref Session-I", RratioObs, 87984149/100000000, 5*10^-9];
expectApprox["M8 R1 Session-I", R1Obs, -12762119/100000000, 5*10^-9];

Print[""];
Print["Stage 245 Mathematica audit passed."];

Exit[0];
