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
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 186 — EXACT MICROSCOPIC SIMILARITY ORBIT"];

Clear[
  chi, delta, eParam, fParam, lam1, c1, gam1, kU, kEta, kW, mu1,
  tau1, scaleS, lambdaW0, cEtaU0, gamma0, kU0, kEta0, kW0, muW0,
  tU0, ell0, sigma0
];

driftVars = {lam1, c1, gam1, kU, kEta, kW, mu1, tau1};
microVars = {lambdaW0, cEtaU0, gamma0, kU0, kEta0, kW0, muW0, tU0};
positiveMicroAssumptions = And @@ Thread[Join[microVars, {ell0, sigma0}] > 0];
$Assumptions =
  Element[{chi, delta, eParam, fParam, lam1, c1, gam1, kU, kEta, kW, mu1, tau1, scaleS}, Reals] &&
  chi > 0 && delta > 0 && positiveMicroAssumptions;

scaleRules = Thread[microVars -> microVars*Exp[scaleS*driftVars]];

chi0Monomial = gamma0*cEtaU0/kU0;
deltaUMonomial = Pi^2*tU0/(ell0^2*kU0);
epsilonWMonomial = gamma0^2*lambdaW0^2*sigma0/(kU0*kW0);
zOverOmegaMonomial = lambdaW0^2*muW0/(kEta0*kW0^2);
epsEtaMonomial = cEtaU0^2/(kU0*kEta0);

ctrMonomial = chi0Monomial^(1 + delta)*deltaUMonomial^(1 + chi);
cntMonomial = zOverOmegaMonomial*epsilonWMonomial^eParam*deltaUMonomial^(-fParam);

logDriftFromMonomial[monomial_] := Module[{scaled, drift},
  scaled = monomial /. scaleRules;
  drift = D[PowerExpand[Log[scaled/monomial]], scaleS] /. scaleS -> 0;
  FullSimplify[drift, Assumptions -> $Assumptions]
];

monomialLogDrifts = logDriftFromMonomial /@ {
  ctrMonomial,
  cntMonomial,
  epsEtaMonomial
};

dlogCtr = monomialLogDrifts[[1]];
dlogCnt = monomialLogDrifts[[2]];
dlogEta = monomialLogDrifts[[3]];

Print["delta log C_tr = ", fmt[Expand[dlogCtr]]];
Print["delta log C_nt = ", fmt[Expand[dlogCnt]]];
Print["delta log eps_eta = ", fmt[Expand[dlogEta]]];

mMat = {
  {0, 1 + delta, 1 + delta, -(2 + chi + delta), 0, 0, 0, 1 + chi},
  {2 + 2*eParam, 0, 2*eParam, fParam - eParam, -1, -(2 + eParam), 1, -fParam},
  {0, 2, 0, -1, -1, 0, 0, 0}
};
mMatDerived = Table[
  Coefficient[Expand[monomialLogDrifts[[row]]], driftVars[[col]]],
  {row, Length[monomialLogDrifts]}, {col, Length[driftVars]}
];
Print[""];
Print["M_* = ", fmt[MatrixForm[mMat]]];
expectZero["M_* row 1 from C_tr monomial matches reference", Total[Abs[mMatDerived[[1]] - mMat[[1]]]]];
expectZero["M_* row 2 from C_nt monomial matches reference", Total[Abs[mMatDerived[[2]] - mMat[[2]]]]];
expectZero["M_* row 3 from eps_eta monomial matches reference", Total[Abs[mMatDerived[[3]] - mMat[[3]]]]];

minor = mMat[[All, {8, 5, 7}]];
detMinor = FullSimplify[Det[minor], Assumptions -> $Assumptions];
Print[""];
Print["minor det(tau1,kEta,mu1) = ", fmt[detMinor]];
If[!TrueQ[detMinor === 1 + chi], fail["Unexpected rank-3 minor determinant", detMinor]];

banner["Linear compatibility solve"];
sol = First[Solve[{dlogCtr == 0, dlogEta == 0, dlogCnt == 0}, {tau1, kEta, mu1}]];
tauFormula = FullSimplify[tau1 /. sol, Assumptions -> $Assumptions];
kEtaFormula = FullSimplify[kEta /. sol, Assumptions -> $Assumptions];
muFormula = FullSimplify[mu1 /. sol, Assumptions -> $Assumptions];

Print["tau1 = ", fmt[tauFormula]];
Print["kEta = ", fmt[kEtaFormula]];
Print["mu1  = ", fmt[muFormula]];

expectZero[
  "tau1 - [kU - ((1+delta)/(1+chi)) (gam1+c1-kU)]",
  tauFormula - (kU - (1 + delta)*(gam1 + c1 - kU)/(1 + chi))
];
expectZero["kEta - (2 c1 - kU)", kEtaFormula - (2*c1 - kU)];
expectZero[
  "mu1 - Stage 185 form",
  muFormula - (
    (2*c1 - kU) + 2*kW - 2*lam1 -
    eParam*(2*gam1 + 2*lam1 - kU - kW) -
    fParam*(1 + delta)*(gam1 + c1 - kU)/(1 + chi)
  )
];

banner["Finite five-parameter similarity orbit"];
Clear[lamSym, cSym, gamSym, uSym, wSym, kEtaScale, muScale, tauScale];
$Assumptions = Element[{chi, delta, eParam, fParam, lamSym, cSym, gamSym, uSym, wSym, kEtaScale, muScale, tauScale}, Reals] &&
  chi > 0 && delta > 0;

orbitLogVector = {lamSym, cSym, gamSym, uSym, kEtaScale, wSym, muScale, tauScale};
scaleSol = First[Solve[Thread[mMat . orbitLogVector == {0, 0, 0}], {kEtaScale, muScale, tauScale}]];
etaExp = FullSimplify[kEtaScale /. scaleSol, Assumptions -> $Assumptions];
tauExp = FullSimplify[tauScale /. scaleSol, Assumptions -> $Assumptions];
muExp = FullSimplify[muScale /. scaleSol, Assumptions -> $Assumptions];

expectZero["solved K_eta scaling matches paper", etaExp - (2*cSym - uSym)];
expectZero[
  "solved T_U scaling matches paper",
  tauExp - (uSym - (1 + delta)*(gamSym + cSym - uSym)/(1 + chi))
];
expectZero[
  "solved mu_W scaling matches paper",
  muExp - (
    (2*cSym - uSym) + 2*wSym - 2*lamSym -
    eParam*(2*gamSym + 2*lamSym - uSym - wSym) -
    fParam*(1 + delta)*(gamSym + cSym - uSym)/(1 + chi)
  )
];

Print["K_eta exponent = ", fmt[etaExp]];
Print["T_U exponent   = ", fmt[tauExp]];
Print["mu_W exponent  = ", fmt[FullSimplify[muExp, Assumptions -> $Assumptions]]];

ctrOrbit = (1 + delta)*(gamSym + cSym - uSym) + (1 + chi)*(tauExp - uSym);
cntOrbit = (2 + 2*eParam)*lamSym + 2*eParam*gamSym + muExp + (fParam - eParam)*uSym - etaExp - (2 + eParam)*wSym - fParam*tauExp;
etaOrbit = 2*cSym - uSym - etaExp;

expectZero["finite orbit preserves C_tr", ctrOrbit];
expectZero["finite orbit preserves C_nt", cntOrbit];
expectZero["finite orbit preserves eps_eta", etaOrbit];
(* Ground check: derive the K_eta^eff scaling from the physical eps_eta monomial. *)
Clear[etaScaling];
$Assumptions = Element[{chi, delta, eParam, fParam, lamSym, cSym, gamSym, uSym, wSym, etaScaling}, Reals] &&
  chi > 0 && delta > 0 && positiveMicroAssumptions;
epsEtaLogDrift = FullSimplify[
  logDriftFromMonomial[epsEtaMonomial] /. {
    lam1 -> 0, c1 -> cSym, gam1 -> 0, kU -> uSym,
    kEta -> etaScaling, kW -> 0, mu1 -> 0, tau1 -> 0
  },
  Assumptions -> $Assumptions
];
solvedEta = etaScaling /. First[Solve[epsEtaLogDrift == 0, etaScaling]];
Print["derived eps_eta log drift = ", fmt[epsEtaLogDrift]];
expectZero["K_eta preserving scaling matches paper 2C-U", solvedEta - (2*cSym - uSym)];
expectZero["chosen Eta_exp solves eps_eta preservation",
  epsEtaLogDrift /. etaScaling -> etaExp];

banner["Linearization reproduces compatibility ledger"];
basisSubs = {lamSym -> lam1, cSym -> c1, gamSym -> gam1, uSym -> kU, wSym -> kW};
expectZero["linearized tau formula", (tauExp /. basisSubs) - tauFormula];
expectZero["linearized kEta formula", (etaExp /. basisSubs) - kEtaFormula];
expectZero["linearized mu formula", (muExp /. basisSubs) - muFormula];

Print[""];
Print["Conclusion:"];
Print["  The three Stage 185 compatibility equations are exactly the tangent-space"];
Print["  equations of a finite five-parameter multiplicative similarity orbit."];
Print["  The coherent weak-axisymmetric zero-defect branch is therefore codimension 3,"];
Print["  not an isolated fine-tuning point."];

Print[""];
Print["Stage 186 Mathematica audit passed."];

Exit[0];
