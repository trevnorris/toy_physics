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

clean[expr_] := Module[{res},
  res = expr /. ConditionalExpression[e_, _] :> e;
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = res /. ConditionalExpression[e_, _] :> e;
  FullSimplify[res, Assumptions -> $Assumptions]
];

expectZero[name_String, expr_] := Module[{res},
  res = clean[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectZeroPacket[name_String, expr_] := Module[{res, zeros},
  res = clean[expr];
  zeros = ConstantArray[0, Dimensions[res]];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === zeros], pass[name], fail[name, res]];
];

firstLogDrift[ratio_] := FullSimplify[
  Coefficient[Normal[Series[PowerExpand[Log[ratio]], {small, 0, 1}]], small, 1],
  Assumptions -> $Assumptions
];

linearCompiler[outputs_, inputs_] := FullSimplify[
  Outer[D, outputs, inputs],
  Assumptions -> $Assumptions
];

banner["STAGE 188 -- BRANCH-OBSERVABLE COMPILER (MATHEMATICA AUDIT)"];

Clear[
  chi, del, eta, small,
  rObs, nObs, etaObs,
  sigTr, sigNt, sigEta,
  theta, xi, rDef
];

$Assumptions = (
  Element[
    {chi, del, eta, small, rObs, nObs, etaObs, sigTr, sigNt, sigEta, theta, xi, rDef},
    Reals
  ] &&
  chi > 0 && del > 0 && 0 < eta < 1
);

obsVars = {rObs, nObs, etaObs};
quotVars = {sigTr, sigNt, sigEta};
defVars = {theta, xi, rDef};

(* Coherent reference coefficients from the appendix definitions. *)
cTr = FullSimplify[
  chi*del/((1 + chi)*(1 + del)*(1 + chi + del)),
  Assumptions -> $Assumptions
];
cStar = FullSimplify[
  ((1 + chi)*(1 + del)*(1 + chi + del))/(chi*del),
  Assumptions -> $Assumptions
];
bStar = FullSimplify[
  2*(1 + chi + del)/del,
  Assumptions -> $Assumptions
];
aTr = FullSimplify[
  2*chi/((1 + chi)*(1 + del)),
  Assumptions -> $Assumptions
];

banner["I. Coherent-branch coefficient identities"];
Print["C_tr,* = ", fmt[cTr]];
Print["C_* = ", fmt[cStar]];
Print["B_* = ", fmt[bStar]];
Print["A_tr,* = ", fmt[aTr]];
expectZero["C_* - 1/C_tr,*", cStar - 1/cTr];
expectZero["A_tr,* - B_* C_tr,*", aTr - bStar*cTr];

banner["II. Observable packet to tangent quotient packet"];

(* Derive the quotient drift from finite branch-composite ratios:
   tracking composite is a fixed negative power of R_tr, the nontracking
   composite is itself the second observable, and epsilon_eta is itself
   the dressing observable. *)
quotFromObsRules = First[
  Solve[
    {
      sigTr == firstLogDrift[Exp[-cStar*small*rObs]],
      sigNt == firstLogDrift[Exp[small*nObs]],
      sigEta == firstLogDrift[Exp[small*etaObs]]
    },
    quotVars,
    Reals
  ]
];
quotFromObs = FullSimplify[quotVars /. quotFromObsRules, Assumptions -> $Assumptions];
obsToQuot = linearCompiler[quotFromObs, obsVars];
sympyObsToQuot = {
  {-cStar, 0, 0},
  {0, 1, 0},
  {0, 0, 1}
};
quotToObs = FullSimplify[Inverse[obsToQuot], Assumptions -> $Assumptions];

Print["Delta_quot^(1) from finite branch drifts = ", fmt[quotFromObs]];
Print["C_obs->quot = ", fmt[obsToQuot]];
Print["det(C_obs->quot) = ", fmt[FullSimplify[Det[obsToQuot], Assumptions -> $Assumptions]]];
expectZeroPacket["observable-to-quotient compiler agrees with SymPy", obsToQuot - sympyObsToQuot];
expectZero["det(C_obs->quot) + C_*", Det[obsToQuot] + cStar];
expectZeroPacket["C_quot->obs * C_obs->quot - I", quotToObs.obsToQuot - IdentityMatrix[3]];
expectZeroPacket["C_obs->quot * C_quot->obs - I", obsToQuot.quotToObs - IdentityMatrix[3]];

banner["III. Tangent quotient packet to defect packet"];

quotToDefRules = First[
  Solve[
    {
      theta == -cTr*sigTr,
      xi == aTr*sigTr + sigNt,
      rDef + xi == -(eta/(1 - eta))*sigEta
    },
    defVars,
    Reals
  ]
];
defFromQuot = FullSimplify[defVars /. quotToDefRules, Assumptions -> $Assumptions];
quotToDef = linearCompiler[defFromQuot, quotVars];
sympyQuotToDef = {
  {-cTr, 0, 0},
  {aTr, 1, 0},
  {-aTr, -1, -eta/(1 - eta)}
};

Print["C_quot->def = ", fmt[quotToDef]];
expectZeroPacket["quotient-to-defect compiler agrees with SymPy", quotToDef - sympyQuotToDef];

defViaQuotFromObs = FullSimplify[
  defFromQuot /. Thread[quotVars -> quotFromObs],
  Assumptions -> $Assumptions
];
expectZero["Theta_from_quot - dln(R_tr)", defViaQuotFromObs[[1]] - rObs];
expectZero[
  "Xi_from_quot - (dln(N_*) - B_* dln(R_tr))",
  defViaQuotFromObs[[2]] - (nObs - bStar*rObs)
];
expectZero[
  "R_from_quot - (-(epseta_*/(1-epseta_*)) dln(epseta) - Xi)",
  defViaQuotFromObs[[3]] - (-(eta/(1 - eta))*etaObs - defViaQuotFromObs[[2]])
];

banner["IV. Direct observable packet to defect packet"];

thetaDriftFromR = firstLogDrift[Exp[small*rObs]];
nCompositeDrift = firstLogDrift[Exp[small*xi]*Exp[small*theta]^bStar];
etaComplementDrift = FullSimplify[
  SeriesCoefficient[Log[(1 - eta*Exp[small*etaObs])/(1 - eta)], {small, 0, 1}],
  Assumptions -> $Assumptions
];

directDefRules = First[
  Solve[
    {
      theta == thetaDriftFromR,
      nObs == nCompositeDrift,
      rDef + xi == etaComplementDrift
    },
    defVars,
    Reals
  ]
];
defFromObsDirect = FullSimplify[defVars /. directDefRules, Assumptions -> $Assumptions];
obsToDefDirect = linearCompiler[defFromObsDirect, obsVars];
obsToDefFactorized = FullSimplify[quotToDef.obsToQuot, Assumptions -> $Assumptions];
sympyObsToDef = {
  {1, 0, 0},
  {-bStar, 1, 0},
  {bStar, -1, -eta/(1 - eta)}
};

Print["Direct defect packet from branch drift equations = ", fmt[defFromObsDirect]];
Print["C_obs->def (direct) = ", fmt[obsToDefDirect]];
Print["det(C_obs->def) = ", fmt[FullSimplify[Det[obsToDefDirect], Assumptions -> $Assumptions]]];
expectZeroPacket["factorized compiler - direct compiler", obsToDefFactorized - obsToDefDirect];
expectZeroPacket["direct compiler agrees with SymPy expected compiler", obsToDefDirect - sympyObsToDef];
expectZero["det(C_obs->def) + epseta_*/(1-epseta_*)", Det[obsToDefDirect] + eta/(1 - eta)];
expectZero["Theta - dln(R_tr)", defFromObsDirect[[1]] - rObs];
expectZero["Xi - (dln(N_*) - B_* dln(R_tr))", defFromObsDirect[[2]] - (nObs - bStar*rObs)];
expectZero[
  "R - (B_* dln(R_tr) - dln(N_*) - epseta_*/(1-epseta_*) dln(epseta))",
  defFromObsDirect[[3]] - (bStar*rObs - nObs - eta/(1 - eta)*etaObs)
];

banner["V. Exact inverse observable compiler"];

inverseObsRules = First[
  Solve[
    {
      theta == thetaDriftFromR,
      nObs == nCompositeDrift,
      rDef + xi == etaComplementDrift
    },
    obsVars,
    Reals
  ]
];
obsFromDefSolved = FullSimplify[obsVars /. inverseObsRules, Assumptions -> $Assumptions];
defToObsSolved = linearCompiler[obsFromDefSolved, defVars];
defToObsInverse = FullSimplify[Inverse[obsToDefDirect], Assumptions -> $Assumptions];
sympyDefToObs = {
  {1, 0, 0},
  {bStar, 1, 0},
  {0, -(1 - eta)/eta, -(1 - eta)/eta}
};

Print["C_def->obs from solving branch equations = ", fmt[defToObsSolved]];
expectZeroPacket["inverse compiler - expected inverse", defToObsSolved - sympyDefToObs];
expectZeroPacket["solved inverse - matrix inverse", defToObsSolved - defToObsInverse];
expectZeroPacket["C_def->obs * C_obs->def - I", defToObsSolved.obsToDefDirect - IdentityMatrix[3]];
expectZeroPacket["C_obs->def * C_def->obs - I", obsToDefDirect.defToObsSolved - IdentityMatrix[3]];
expectZeroPacket[
  "observable roundtrip - Delta_obs^(1)",
  defToObsSolved.(obsToDefDirect.obsVars) - obsVars
];

banner["VI. Complementary selected-branch observable"];

Print["delta ln(1 - epseta) = ", fmt[etaComplementDrift]];
expectZero[
  "(R + Xi) - delta ln(1-epseta)",
  (defFromObsDirect[[3]] + defFromObsDirect[[2]]) - etaComplementDrift
];

banner["VII. Zero-set equivalence (shared zero set via invertibility)"];

expectZeroPacket[
  "C_obs->def then inverse recovers generic obs (bijection, def side)",
  defToObsSolved.(obsToDefDirect.obsVars) - obsVars
];
expectZeroPacket[
  "C_obs->quot then inverse recovers generic obs (bijection, quot side)",
  quotToObs.(obsToQuot.obsVars) - obsVars
];

detObsToDef = FullSimplify[Det[obsToDefDirect], Assumptions -> $Assumptions];
detObsToQuot = FullSimplify[Det[obsToQuot], Assumptions -> $Assumptions];
expectZero[
  "1/det(C_obs->def) well-defined (nonzero det)",
  detObsToDef*(1/detObsToDef) - 1
];
expectZero[
  "1/det(C_obs->quot) well-defined (nonzero det)",
  detObsToQuot*(1/detObsToQuot) - 1
];

Print[""];
Print["Stage 188 Mathematica audit passed."];

Exit[0];
