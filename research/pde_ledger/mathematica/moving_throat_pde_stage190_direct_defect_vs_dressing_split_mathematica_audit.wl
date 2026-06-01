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

stripCE[expr_] := expr /. ConditionalExpression[e_, _] :> e;

reduce[expr_] := Module[{res},
  res = stripCE[expr];
  res = FullSimplify[Together[res], Assumptions -> $Assumptions];
  res = stripCE[res];
  FullSimplify[Cancel[res], Assumptions -> $Assumptions]
];

expectZero[name_String, expr_] := Module[{res},
  res = reduce[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectNonzero[name_String, expr_] := Module[{res},
  res = reduce[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], fail[name, res], pass[name]];
];

banner["STAGE 190 -- DIRECT DEFECT VS DRESSING SPLIT"];

Clear[
  zeta, zW, chi, omegaW2, eps, epsEta, lambda0, deltaU, branchR, badPacket,
  lambdaW, muW, gamma, cEtaU, kU, kEta, kW, tU, ell, sig,
  lambdaOne, muOne, gammaOne, cOne, kappaU, kappaEta, kappaW, tauOne,
  sigmaZ, sigmaChi, sigmaEps, sigmaDelta, sigmaEta, epsW,
  sigmaTrVar, sigmaNtVar, sigmaEtaVar, thetaVar, xiVar, dressVar,
  epsAx, x0, x1, y0, y1
];

$Assumptions =
  Element[
    {
      zeta, zW, chi, omegaW2, eps, epsEta, lambda0, deltaU, branchR,
      badPacket, lambdaW, muW, gamma, cEtaU, kU, kEta, kW, tU, ell, sig,
      lambdaOne, muOne, gammaOne, cOne, kappaU, kappaEta, kappaW, tauOne,
      sigmaZ, sigmaChi, sigmaEps, sigmaDelta, sigmaEta, epsW,
      sigmaTrVar, sigmaNtVar, sigmaEtaVar, thetaVar, xiVar, dressVar,
      epsAx, x0, x1, y0, y1
    },
    Reals
  ] &&
  0 < zeta < 1 && zW > 0 && chi > 0 && omegaW2 > 0 &&
  0 < eps < 1 && 0 < epsEta < 1 && lambda0 > 0 && deltaU > 0 &&
  branchR > 0 && lambdaW > 0 && muW > 0 && gamma > 0 &&
  cEtaU > 0 && kU > 0 && kEta > 0 && kW > 0 && tU > 0 &&
  ell > 0 && sig > 0 && epsW > 0;

(* I. Support-blindness is checked by rebuilding the loaded demand from the
   mass packet and then cancelling it back to the direct transfer shape. *)
banner["I. Support-blind transfer and selected-branch identities"];

mixMass =
  8*zW*(1 + chi)^2/(Pi^2*(1 - epsEta)*(1 - eps));
supportMass =
  8*zeta*zW*(1 + chi)^2/(Pi^2*(1 - epsEta)*(1 - zeta*eps));
supportBoost = 1 + zeta*(1 - eps)/(1 - zeta*eps);
loadedMass = mixMass + supportMass;

loadedTarget =
  8*lambda0*omegaW2*(1 - eps)*supportBoost/(Pi^2*loadedMass);
loadedTransfer = lambda0*(1 - epsEta)/loadedTarget;

directTransfer = zW*(1 + chi)^2/(omegaW2*(1 - eps)^2);
directTarget =
  lambda0*omegaW2*(1 - epsEta)*(1 - eps)^2/(zW*(1 + chi)^2);
branchExponent = 2*(1 + chi + deltaU)/deltaU;
branchShape = (1 + chi + deltaU)/((1 + chi)*(1 + deltaU));
loadedComposite = loadedTransfer*branchR^branchExponent;
directComposite = directTransfer*branchShape^branchExponent;
selectedProduct = loadedTarget*loadedTransfer/lambda0;

expectZero["M_tr - M_mix S(zeta;eps)", loadedMass - mixMass*supportBoost];
expectZero["support-loaded T^2 reconstruction", loadedTransfer - directTransfer];
expectZero["support-loaded R_target reconstruction", loadedTarget - directTarget];
expectZero[
  "support-loaded N_* reconstruction",
  (loadedComposite /. branchR -> branchShape) - directComposite
];
expectZero["E - (1 - epsilon_eta)", selectedProduct - (1 - epsEta)];
expectZero["dln(T^2_loaded)/dzeta", D[Log[loadedTransfer], zeta]];
expectZero["dln(R_target_loaded)/dzeta", D[Log[loadedTarget], zeta]];
expectZero["dln(N_* loaded)/dzeta", D[Log[loadedComposite], zeta]];

spoiledTarget =
  8*lambda0*omegaW2*(1 - eps)*supportBoost/
    (Pi^2*(loadedMass + zeta*mixMass));
spoiledResidual = D[Log[spoiledTarget], zeta];
expectNonzero[
  "spoiled dln(R_target)/dzeta is not identically zero",
  spoiledResidual
];
expectZero[
  "spoiled exact witness at eps=1/3,zeta=1/2",
  (spoiledResidual /. {eps -> 1/3, zeta -> 1/2}) + 46/133
];

(* II. The primitive microscopic ledger is extracted with logarithmic Euler
   operators on the defining monomials, rather than by a one-parameter drift. *)
banner["II. Microscopic slippage ledger"];

primitiveVariables = {lambdaW, muW, gamma, cEtaU, kU, kEta, kW, tU};
primitiveDrifts = {
  lambdaOne, muOne, gammaOne, cOne, kappaU, kappaEta, kappaW, tauOne
};
logEuler[monomial_] := reduce[
  Total[
    MapThread[
      #2*#1*D[Log[monomial], #1] &,
      {primitiveVariables, primitiveDrifts}
    ]
  ]
];

sigmaChiFromKernel = logEuler[gamma*cEtaU/kU];
sigmaDeltaFromKernel = logEuler[Pi^2*tU/(ell^2*kU)];
sigmaEtaFromKernel = logEuler[cEtaU^2/(kU*kEta)];
sigmaEpsFromKernel = logEuler[gamma^2*lambdaW^2*sig/(kU*kW)];
sigmaZFromKernel = logEuler[lambdaW^2*muW/(kEta*kW^2)];

expectZero["Sigma_chi - (gamma1 + c1 - kappaU)", sigmaChiFromKernel - (gammaOne + cOne - kappaU)];
expectZero["Sigma_delta - (tau1 - kappaU)", sigmaDeltaFromKernel - (tauOne - kappaU)];
expectZero["Sigma_eta - (2 c1 - kappaU - kappaEta)", sigmaEtaFromKernel - (2*cOne - kappaU - kappaEta)];
expectZero[
  "Sigma_epsilon - (2 gamma1 + 2 lambda1 - kappaU - kappaW)",
  sigmaEpsFromKernel - (2*gammaOne + 2*lambdaOne - kappaU - kappaW)
];
expectZero[
  "Sigma_Z - (2 lambda1 + mu1 - kappaEta - 2 kappaW)",
  sigmaZFromKernel - (2*lambdaOne + muOne - kappaEta - 2*kappaW)
];

(* III. Remove the tracking direction by projection along Sigma_chi.  The
   residual is the nontracking coordinate and contains no Sigma_chi term. *)
banner["III. Direct-defect projection"];

xiDirect =
  sigmaZ + 2*chi/(1 + chi)*sigmaChi +
    2*epsW/(1 - eps)*
      (
        (11 + 9*deltaU)/(11*(1 + deltaU))*sigmaEps -
        2*deltaU/(11*(1 + deltaU)^2)*sigmaDelta
      );
trackingCoordinate = (1 + chi)*sigmaDelta + (1 + deltaU)*sigmaChi;
atrProjected =
  reduce[Coefficient[xiDirect, sigmaChi]/Coefficient[trackingCoordinate, sigmaChi]];
sigmaNtProjected = reduce[xiDirect - atrProjected*trackingCoordinate];

atrTarget = 2*chi/((1 + chi)*(1 + deltaU));
ctrTarget =
  chi*deltaU/((1 + chi)*(1 + deltaU)*(1 + chi + deltaU));
sigmaNtTarget =
  sigmaZ +
    2*epsW/(1 - eps)*(11 + 9*deltaU)/(11*(1 + deltaU))*sigmaEps -
    (
      2*chi/(1 + deltaU) +
      4*epsW*deltaU/(11*(1 - eps)*(1 + deltaU)^2)
    )*sigmaDelta;

expectZero["projected A_tr", atrProjected - atrTarget];
expectZero["nontracking residual has no Sigma_chi", Coefficient[sigmaNtProjected, sigmaChi]];
expectZero["projected Sigma_nt", sigmaNtProjected - sigmaNtTarget];
expectZero["Xi_direct - (A_tr Sigma_tr + Sigma_nt)", xiDirect - (atrTarget*trackingCoordinate + sigmaNtTarget)];

(* IV. The direct/dressing split is inverted by native matrix algebra. *)
banner["IV. Triangular direct-defect / dressing compiler"];

compiler = {
  {-ctrTarget, 0, 0},
  {atrTarget, 1, 0},
  {0, 0, -epsEta/(1 - epsEta)}
};
normalVector = {sigmaTrVar, sigmaNtVar, sigmaEtaVar};
observableVector = compiler.normalVector;
thetaExpr = observableVector[[1]];
xiExpr = observableVector[[2]];
dressExpr = observableVector[[3]];

Print["compiler = ", fmt[compiler]];
Print["det(compiler) = ", fmt[reduce[Det[compiler]]]];
expectZero["determinant formula", Det[compiler] - ctrTarget*epsEta/(1 - epsEta)];
expectZero["dXi/dSigma_eta", D[xiExpr, sigmaEtaVar]];
expectZero["d(R+Xi)/dSigma_tr", D[dressExpr, sigmaTrVar]];
expectZero["d(R+Xi)/dSigma_nt", D[dressExpr, sigmaNtVar]];

inverseByLinearSolve = reduce[LinearSolve[compiler, {thetaVar, xiVar, dressVar}]];
trackingRecovered =
  -((1 + chi)*(1 + deltaU)*(1 + chi + deltaU))/(chi*deltaU)*thetaVar;
nontrackingRecovered = xiVar + atrTarget/ctrTarget*thetaVar;
dressingRecovered = -(1 - epsEta)/epsEta*dressVar;

expectZero["Sigma_tr inverse by LinearSolve", inverseByLinearSolve[[1]] - trackingRecovered];
expectZero[
  "A_tr/C_tr - 2(1+chi+delta)/delta",
  atrTarget/ctrTarget - 2*(1 + chi + deltaU)/deltaU
];
expectZero["Sigma_nt inverse by LinearSolve", inverseByLinearSolve[[2]] - nontrackingRecovered];
expectZero["Sigma_eta inverse by LinearSolve", inverseByLinearSolve[[3]] - dressingRecovered];
expectZero[
  "inverse reconstructs Sigma_tr",
  (trackingRecovered /. thetaVar -> thetaExpr) - sigmaTrVar
];
expectZero[
  "inverse reconstructs Sigma_nt",
  (nontrackingRecovered /. {thetaVar -> thetaExpr, xiVar -> xiExpr}) - sigmaNtVar
];
expectZero[
  "inverse reconstructs Sigma_eta",
  (dressingRecovered /. dressVar -> dressExpr) - sigmaEtaVar
];
expectZero["tracking-rigid Xi - Sigma_nt", (xiExpr /. sigmaTrVar -> 0) - sigmaNtVar];
expectZero["grouped-defect cancellation theorem", xiExpr /. sigmaNtVar -> -atrTarget*sigmaTrVar];
expectZero["selected-branch rigidity theorem", dressExpr /. sigmaEtaVar -> 0];

(* V. Project the grouped real P2 lane into scalar and anisotropic pieces, then
   use the invariant quadratic form to show the scalar feed-down is O(eps^2). *)
banner["V. Grouped real P2 scalar no-go filter"];

p2Signature = {1, 1/2, -1};
projectors = {
  {1, 2, 2}/5,
  {2, -1, -1}/10,
  {0, 1, -1}/2
};
xLane = x0*{1, 1, 1} + epsAx*x1*p2Signature;
yLane = y0*{1, 1, 1} + epsAx*y1*p2Signature;
xProjected = reduce[projectors.xLane];
yProjected = reduce[projectors.yLane];

xBar = xProjected[[1]];
ax = xProjected[[2]];
bx = xProjected[[3]];
yBar = yProjected[[1]];
ay = yProjected[[2]];
by = yProjected[[3]];
invariantKernel = DiagonalMatrix[{4, 4/5}];
iXY = reduce[{ax, bx}.invariantKernel.{ay, by}];
iXX = reduce[{ax, bx}.invariantKernel.{ax, bx}];

expectZero["xbar - x0", xBar - x0];
expectZero["ybar - y0", yBar - y0];
expectZero["b_x - 3 a_x", bx - 3*ax];
expectZero["b_y - 3 a_y", by - 3*ay];
expectZero["I[x,y] - 7/10 eps^2 x1 y1", iXY - 7/10*epsAx^2*x1*y1];
expectZero["I[x,x] - 7/10 eps^2 x1^2", iXX - 7/10*epsAx^2*x1^2];
expectZero[
  "linear term in xbar",
  Coefficient[Normal[Series[xBar, {epsAx, 0, 1}]], epsAx, 1]
];
expectZero[
  "linear term in I[x,x] at eps=0",
  Coefficient[Normal[Series[iXX, {epsAx, 0, 1}]], epsAx, 1]
];

banner["STAGE 190 LEDGER"];
Print["1. The support-loaded branch algebra cancels exactly back to T^2, R_target, and N_*."];
Print["2. Logarithmic Euler operators reproduce the five microscopic slippage laws."];
Print["3. Projection along the tracking coordinate leaves the corrected nontracking defect."];
Print["4. The direct-defect/dressing compiler is triangular and exactly invertible."];
Print["5. Pure grouped real P2 anisotropy has no linear scalar feed-down."];

Exit[0];
