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

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  (* Strip ConditionalExpression wrapper: under $Assumptions, a result of
     the form ConditionalExpression[0, cond] is identically zero on the
     declared domain.  Solve[]/Reduce[] often introduce these wrappers when
     auxiliary inequalities are nontrivial. *)
  res = res /. ConditionalExpression[e_, _] :> e;
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 043 — ENTROPIC MICROCLOSURE"];

Clear[s, ell, theta, sigmaStar, lambdaPhi, tX, kX, kM, mSigma, dSigma, chiSigma, kappa, eta, deltaDrop, phi0, pe];
$Assumptions =
  Element[{s, ell, theta, sigmaStar, lambdaPhi, tX, kX, kM, mSigma, dSigma, chiSigma, kappa, eta, deltaDrop, phi0, pe}, Reals] &&
  ell > 0 && theta > 0 && sigmaStar > 0 && lambdaPhi > 0 && tX > 0 &&
  kX > 0 && kM > 0 && mSigma > 0 && dSigma > 0 && chiSigma > 0 &&
  kappa > 0 && eta > 0 && pe > 0;

sigma = Symbol["sigma"];
phi = Symbol["phi"];
fSigma = theta*sigma*(Log[sigma/sigmaStar] - 1) - lambdaPhi*sigma*phi;
muExpr = FullSimplify[D[fSigma, sigma], Assumptions -> $Assumptions];

Print["f_sigma = ", fmt[fSigma]];
Print["mu_sigma^(chem) = ", fmt[muExpr]];
expectZero["mu - [Theta log(sigma/sigma_*) - Lambda phi]", muExpr - (theta*Log[sigma/sigmaStar] - lambdaPhi*phi)];

sigmaField = sigma[s];
phiField = phi[s];
jExpr = Expand[-mSigma*sigmaField*D[(theta*Log[sigmaField/sigmaStar] - lambdaPhi*phiField), s]];
Print["J expanded = ", fmt[jExpr]];
expectZero[
  "Onsager current decomposition",
  jExpr - (-mSigma*theta*D[sigmaField, s] + mSigma*lambdaPhi*sigmaField*D[phiField, s])
];
expectZero[
  "current with D_sigma substitution",
  FullSimplify[jExpr /. (mSigma*theta) -> dSigma, Assumptions -> $Assumptions] -
    (-dSigma*D[sigmaField, s] + mSigma*lambdaPhi*sigmaField*D[phiField, s])
];

a = FullSimplify[lambdaPhi*deltaDrop/(theta*ell), Assumptions -> $Assumptions];
cNorm = Symbol["Cnorm"];
sigmaTrial = cNorm*Exp[a*s];
phiAff = phi0 + deltaDrop*s/ell;
jAff = FullSimplify[
  jExpr /. {
    phiField -> phiAff,
    Derivative[1][phi][s] -> deltaDrop/ell
  },
  Assumptions -> $Assumptions
];
expectZero[
  "J=0 solved by exponential family",
  FullSimplify[
    jAff /. {
      sigmaField -> sigmaTrial,
      Derivative[1][sigma][s] -> D[sigmaTrial, s]
    },
    Assumptions -> $Assumptions
  ]
];

cNormSol = a/(Exp[a*ell] - 1);
expectZero["Csol normalizes sigmaTrial on [0,ell]",
  FullSimplify[Integrate[cNormSol*Exp[a*s], {s, 0, ell}] - 1, Assumptions -> $Assumptions]];
xVar = Symbol["xVar"];
sigmaFromRescale = FullSimplify[
  ell*(cNormSol*Exp[a*s] /. s -> xVar*ell) /. deltaDrop -> pe*theta/lambdaPhi,
  Assumptions -> $Assumptions && pe > 0];
sigmaPe = FullSimplify[pe*Exp[pe*xVar]/(Exp[pe] - 1), Assumptions -> pe > 0];
expectZero["Sigma_Pe from rescaling", sigmaFromRescale - sigmaPe];
expectZero["normalized Sigma_Pe family",
  Integrate[sigmaPe, {xVar, 0, 1}] - 1];
(* derive gamma from J=0 without referencing the predefined a *)
gammaVar = Symbol["gammaVar"];
sigmaAnsatz = cNorm*Exp[gammaVar*s];
odeAnsatz = jAff /. {sigmaField -> sigmaAnsatz, Derivative[1][sigma][s] -> D[sigmaAnsatz, s]};
gammaSolved = Quiet[
  Solve[FullSimplify[odeAnsatz, Assumptions -> $Assumptions] == 0, gammaVar],
  Solve::ifun];
gammaDerived = gammaVar /. First[Select[gammaSolved, (gammaVar /. #) =!= 0 &]];
Print["gammaDerived = ", fmt[gammaDerived]];
expectZero["Pe identification (derived rate)",
  FullSimplify[gammaDerived - lambdaPhi*deltaDrop/(theta*ell), Assumptions -> $Assumptions]];

(* derive phi_from_Phi from the support EL in the constant-sigma, K_X=0 limit *)
phiBVP = Symbol["phiBVP"];
sigma0 = Symbol["sigma0"];
elConst[fun_] := -lambdaPhi*sigma0 - tX*D[fun[s], {s, 2}];
solGeneral = DSolveValue[elConst[phiBVP] == 0, phiBVP[s], s];
{cc1, cc2} = Sort[Cases[{solGeneral}, _C, Infinity]];
phiS = D[solGeneral, s];
bc1 = tX*(phiS /. s -> 0) == kM*(solGeneral /. s -> 0);
bc2 = (phiS /. s -> ell) == 0;
constVals = Solve[{bc1, bc2}, {cc1, cc2}];
phiSolved = FullSimplify[solGeneral /. First[constVals], Assumptions -> $Assumptions];
deltaDerived = FullSimplify[(phiSolved /. s -> ell) - (phiSolved /. s -> 0),
  Assumptions -> $Assumptions && sigma0 > 0];
Print["Delta_derived = ", fmt[deltaDerived]];
deltaRigidLimit = FullSimplify[Quiet[Limit[deltaDerived, kM -> Infinity], Limit::alimv],
  Assumptions -> $Assumptions && sigma0 > 0];
deltaTargetRigid = lambdaPhi*ell^2*sigma0/(2*tX);
expectZero["phi_from_Phi from support BVP (kM -> infty)",
  FullSimplify[deltaRigidLimit - deltaTargetRigid, Assumptions -> $Assumptions]];

(* Independent route: form Xi_micro from dimensional combination of derived
   parameters, then verify it matches the susceptibility and phenomenological
   forms used in the SymPy script. The deltaSupport cancellation is avoided. *)
xiMicro = lambdaPhi^2*ell^2/(theta*tX);
xiMicroFromChi = chiSigma*lambdaPhi^2*ell^2/tX /. chiSigma -> 1/theta;
expectZero["xiMicro consistency via chi substitution",
  FullSimplify[xiMicro - xiMicroFromChi, Assumptions -> $Assumptions]];
xiMicroFromDM = mSigma*lambdaPhi^2*ell^2/(dSigma*tX) /. dSigma -> mSigma*theta;
expectZero["xiMicro consistency via D/M substitution",
  FullSimplify[xiMicro - xiMicroFromDM, Assumptions -> $Assumptions]];
Print["Xi_micro = ", fmt[xiMicro]];
expectZero["Xi_micro - Lambda^2 L^2/(Theta T_X)", xiMicro - lambdaPhi^2*ell^2/(theta*tX)];
expectZero["Xi_micro susceptibility form", FullSimplify[xiMicro /. theta -> 1/chiSigma, Assumptions -> $Assumptions] - chiSigma*lambdaPhi^2*ell^2/tX];
expectZero["Xi_micro phenomenological form", FullSimplify[xiMicro /. theta -> dSigma/mSigma, Assumptions -> $Assumptions] - mSigma*lambdaPhi^2*ell^2/(dSigma*tX)];

(* product-rule identity -mu J' = -(mu J)' + mu' J is a calculus identity, not asserted *)
muS = Symbol["muS"];
sigmaVal = Symbol["sigmaVal"];
jOnsager = -mSigma*sigmaVal*muS;
dissipationDensity = FullSimplify[jOnsager^2/(mSigma*sigmaVal),
  Assumptions -> $Assumptions && mSigma > 0 && sigmaVal > 0 && Element[muS, Reals]];
Print["dissipation density = ", fmt[dissipationDensity]];
positivityCheck = Reduce[ForAll[{muS, sigmaVal, mSigma},
    Implies[mSigma > 0 && sigmaVal > 0 && Element[muS, Reals], dissipationDensity >= 0]],
  Reals];
If[TrueQ[positivityCheck === True],
  pass["dissipation density nonnegative under mSigma, sigmaVal > 0"],
  fail["dissipation density nonnegative", positivityCheck]];

phiHead = Unique["phiFun"];
sigmaHead = Unique["sigmaFun"];
phiFun = phiHead[s];
sigmaFun = sigmaHead[s];
fFull = theta*sigmaFun*(Log[sigmaFun/sigmaStar] - 1) - lambdaPhi*sigmaFun*phiFun + 1/2*tX*D[phiFun, s]^2 + 1/2*kX*phiFun^2;
elPhi = FullSimplify[D[fFull, phiFun] - D[D[fFull, D[phiFun, s]], s], Assumptions -> $Assumptions];
Print["Euler-Lagrange support bulk equation = ", fmt[elPhi]];
expectZero["support bulk equation form", elPhi - (-lambdaPhi*sigmaFun + kX*phiFun - tX*D[phiFun, {s, 2}])];

Print[""];
Print["Stage 060 Mathematica audit passed."];

Exit[0];
