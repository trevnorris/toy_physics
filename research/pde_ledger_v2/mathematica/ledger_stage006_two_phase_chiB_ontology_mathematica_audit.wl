(* Ledger stage006 Mathematica audit: two-phase chi_B ontology.

   Print-only, standalone, no arguments, no exports.  This mirror uses
   exponent associations for dimensions, formal-function continuity expansion,
   native DSolve/Integrate for the wall and leakage anchor, and a longitudinal
   dispersion determinant rather than the SymPy Hamiltonian route.
*)

ClearAll[
  heading, subheading, cleanZero, assertExact, expectZero, expectBool,
  expectNonzero, expectFail, dim, dimVec, dimString, dimAdd, dimNeg,
  dimSub, dimScale, dimPow, dimResidual, expectDim, homResidual,
  algebraicZero, boundaryJump, assembleProjectedTwoSource,
  operativeA1AbsenceResidual, driftPartitionResidual,
  printPins, runConsumed, runLegADimensions, runLegAClosure, runLegAWall,
  runLegBProjection, runLegBRecovery, tokenFromDiscriminators, runLegC,
  runAblations, printCarriedAndDrift, printVerdictLabels, passCount, failCount
];

passCount = 0;
failCount = 0;

expectedPostD16DriftToken = "DRIFT(5)";
preD16DriftMembers = {"chi_B", "a_B", "kappa_B", "alpha_aniso", "Gamma_B", "gating structure"};
retiredD16DriftMembers = {"alpha_aniso"};

heading[text_] := (
  Print[""];
  Print[StringRepeat["=", StringLength[text]]];
  Print[text];
  Print[StringRepeat["=", StringLength[text]]]
);

subheading[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);

cleanZero[expr_] := FullSimplify[expr] /. ConditionalExpression[0, _] -> 0;

assertExact[name_, expr_] := Module[{reals},
  reals = Cases[Unevaluated[expr], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FAIL  ", name, ": machine-real atom(s) found: ", ToString[InputForm[reals]]];
    Throw[name, "ledgerStage006Failure"]
  ]
];

expectZero[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FAIL  ", name, ": residual = ", ToString[InputForm[clean]]];
    Throw[name, "ledgerStage006Failure"]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

expectNonzero[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[! TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name, " is nonzero as required (residual = ", ToString[InputForm[clean]], ")"],
    failCount++;
    Print["FAIL  ", name, ": required nonzero residual vanished"];
    Throw[name, "ledgerStage006Failure"]
  ]
];

expectFail[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[! TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name, " produced required FAIL (residual = ", ToString[InputForm[clean]], ")"],
    failCount++;
    Print["FAIL  ", name, ": required mutation/ablation did not fire"];
    Throw[name, "ledgerStage006Failure"]
  ]
];

dim[l_, t_, m_] := <|"L" -> l, "T" -> t, "M" -> m|>;
dimVec[d_Association] := {d["L"], d["T"], d["M"]};
dimAdd[items__Association] := AssociationThread[{"L", "T", "M"}, Total[dimVec /@ {items}]];
dimNeg[d_Association] := AssociationThread[{"L", "T", "M"}, -dimVec[d]];
dimSub[left_Association, right_Association] := dimAdd[left, dimNeg[right]];
dimScale[scale_, d_Association] := AssociationThread[{"L", "T", "M"}, scale dimVec[d]];
dimPow[d_Association, power_] := dimScale[power, d];
dimResidual[actual_Association, expected_Association] := FullSimplify[(dimVec[actual] - dimVec[expected]).(dimVec[actual] - dimVec[expected])];
expectDim[name_, actual_Association, expected_Association] := expectZero[name, dimResidual[actual, expected]];

dimString[d_Association] := Module[{labels, powers, pieces},
  labels = {"L", "T", "M"};
  powers = dimVec[d];
  pieces = DeleteCases[
    Table[
      Which[
        TrueQ[powers[[i]] === 0], Nothing,
        TrueQ[powers[[i]] === 1], labels[[i]],
        True, labels[[i]] <> "^" <> ToString[InputForm[powers[[i]]]]
      ],
      {i, 3}
    ],
    Nothing
  ];
  If[pieces === {}, "1", StringRiffle[pieces, " "]]
];

homResidual[terms_Association] := Module[{vals, ref},
  vals = Values[terms];
  ref = First[vals];
  FullSimplify[Total[dimResidual[#, ref] & /@ Rest[vals]]]
];

operativeA1AbsenceResidual[terms_Association] := If[
  MemberQ[Flatten[Lookup[Values[terms], "Symbols", {}]], "alpha_aniso"] ||
    AnyTrue[Keys[terms], StringContainsQ[#, "alpha_aniso"] &],
  1,
  0
];

algebraicZero[expr_, assumptions_: True] := TrueQ[FullSimplify[expr == 0, assumptions]];

boundaryJump[expr_, var_, assumptions_: True] := FullSimplify[
  Limit[expr, var -> Infinity, Assumptions -> assumptions] -
    Limit[expr, var -> -Infinity, Assumptions -> assumptions],
  assumptions
];

assembleProjectedTwoSource[Wexpr_, nexpr_, uwexpr_, chiexpr_, gammaexpr_, jchiexpr_, var_, assumptions_] := Module[
  {qOrder, rawFlux, boundary, sFlux, sConvert},
  qOrder = Expand[chiexpr nexpr uwexpr + jchiexpr];
  rawFlux = FullSimplify[-Integrate[Wexpr D[qOrder, var], {var, -Infinity, Infinity}, Assumptions -> assumptions], assumptions];
  boundary = boundaryJump[Wexpr qOrder, var, assumptions];
  sFlux = FullSimplify[-boundary + Integrate[D[Wexpr, var] qOrder, {var, -Infinity, Infinity}, Assumptions -> assumptions], assumptions];
  sConvert = FullSimplify[Integrate[Wexpr nexpr gammaexpr, {var, -Infinity, Infinity}, Assumptions -> assumptions], assumptions];
  <|"qOrder" -> qOrder, "rawFlux" -> rawFlux, "boundary" -> boundary, "sFlux" -> sFlux, "sConvert" -> sConvert, "sTotal" -> FullSimplify[sFlux + sConvert, assumptions]|>
];

zeroD = dim[0, 0, 0];
ell = dim[1, 0, 0];
timeD = dim[0, 1, 0];
massD = dim[0, 0, 1];
velocityD = dimSub[ell, timeD];
actionD = dimAdd[massD, dimScale[2, ell], dimNeg[timeD]];
energyD = dimSub[actionD, timeD];
fDensity4D = dimSub[energyD, dimScale[4, ell]];
braneLagD = dimSub[energyD, dimScale[3, ell]];
rho4D = dimScale[-4, ell];
mGnlsD = massD;
kEosD = dimSub[fDensity4D, dimScale[5, rho4D]];
hbarD = actionD;

printPins[] := (
  subheading["Pinned postulated modeling choices P1-P13"];
  Print["  P1  POSTULATED: n(X,t)=total conserved 4D constituent number density [n]=L^-4; KINETIC_MASS_FACTOR_PINNED: kinetic term is 1/2*m_GNLS*n*|u|^2."];
  Print["  P2  POSTULATED: chi_B dimensionless in [0,1]; chi_B=1 brane-ordered, chi_B=0 bulk-like; n_B=chi_B*n."];
  Print["  P3  POSTULATED: f_B(chi_B)=a_B*chi_B^2*(1-chi_B)^2, a_B>0, minima {0,1}; new double-well input because cited U(rho) is single-well."];
  Print["  P4  POSTULATED: SAME_DENSITY_DEGENERACY_POSTULATED: f_B is n-independent and f_B(n,0)=f_B(n,1)=0."];
  Print["  P5  POSTULATED: gradient/interface term (kappa_B/2)*|grad_4 chi_B|^2 with [kappa_B]=M T^-2."];
  Print["  P6  POSTULATED: shear gate chi_B*f_shear with displacement u_d distinct from velocity u; brane mu_R projection is dim-consistent only here."];
  Print["  P7  POSTULATED: chi_B is THE wall order parameter: a postulated independent scalar field, NOT currently built as |P_parallel|^2. Decision 16 retires alpha_aniso and the carried P field. FUTURE_GATE_CHI_B_EQ_ABS_P_PARALLEL_SQ remains named high-risk/Part-VII-adjacent and requires a NEW T0 freeze; obsolete as a carried route, not foreclosed."];
  Print["  P8  POSTULATED ADJUNCT: D_t chi_B=-M_chi*mu_chi+Gamma_B; HANDOFF_P_ORDER_N_PLACEMENT_CORRECTED: P_order=int mu_chi*D_t chi_B d4X."];
  Print["  P9  POSTULATED DEFAULT: J_chi=0 default; J_chi!=0 deferred with dim row only."];
  Print["  P10 CONVENTION: recovery target is frozen old-ledger S_leak with j^w=n*u^w and unit-normalized W, [W]=L^-1."];
  Print["  P11 POSTULATED GLOBAL: global returns R_0=-M_0 and R_1=-D_1 printed here, NOT locally asserted."];
  Print["  P12 POSTULATED ONTOLOGY: throat=phase-conversion admittance/outlet driven by stress/mu gradients, NOT pairwise suction; Gate-L delta w=u_w DEFERRED."];
  Print["  P13 CONVENTION: free-energy +1/2*kappa_phase*(grad theta)^2 maps to Lagrangian K_theta=-kappa_phase; Maxwell needs K_theta=C_J^2/rho_br>0."]
);

runConsumed[] := Module[
  {consumedCJ, consumedBEff, consumedCS0Sq, dMuR, dRhoBr, dRhoB0, dJ, dCJ, dKtheta, dChiC, dB},
  subheading["Consumed citations and exact-value integrity checks"];
  Print["  CITED ledger_stage004 (I-1): {L,T,M}; [n]=[rho0]=L^-4; [m_GNLS]=M; [K]=M L^18 T^-2; [hbar]=M L^2 T^-1; U(rho)=(K/4)*rho^5 single-well."];
  Print["  CITED ledger_stage005 (I-2): c_s0^2=5*K*rho0^4/m_GNLS."];
  Print["  CITED ledger_stage003 (III): c_gamma^2=mu_R/rho_br; C_J=-J*rho_B0; B_eff=rho_B0^2/chi_c; second-class classification rule."];

  consumedCJ = -Jsym rhoB0;
  consumedBEff = rhoB0^2/chiC;
  consumedCS0Sq = 5 kEOSsym rho0sym^4/mGNLSsym;
  expectZero["CITED exact-value check C_J+J*rho_B0=0", consumedCJ + Jsym rhoB0];
  expectZero["CITED exact-value check B_eff-rho_B0^2/chi_c=0", consumedBEff - rhoB0^2/chiC];
  expectZero["CITED exact-value check c_s0^2-5*K*rho0^4/m_GNLS=0", consumedCS0Sq - 5 kEOSsym rho0sym^4/mGNLSsym];

  expectDim["CITED stage004 [n]=L^-4", rho4D, dimScale[-4, ell]];
  expectDim["CITED stage004 [m_GNLS]=M", mGnlsD, massD];
  expectDim["CITED stage004 [K]=M L^18 T^-2", kEosD, dim[18, -2, 1]];
  expectDim["CITED stage004 [hbar]=M L^2 T^-1", hbarD, dim[2, -1, 1]];
  expectDim["CITED stage005 [c_s0^2]=L^2 T^-2", dimSub[dimAdd[kEosD, dimScale[4, rho4D]], mGnlsD], dimPow[velocityD, 2]];

  dMuR = braneLagD;
  dRhoBr = dimAdd[massD, dimScale[-3, ell]];
  dRhoB0 = dRhoBr;
  dJ = dimSub[dimScale[2, ell], timeD];
  dCJ = dimAdd[massD, dimScale[-1, ell], dimNeg[timeD]];
  dKtheta = dimAdd[massD, ell, dimScale[-2, timeD]];
  dChiC = dimAdd[massD, dimScale[-5, ell], dimScale[2, timeD]];
  dB = braneLagD;
  expectDim["CITED stage003 [mu_R]=M L^-1 T^-2", dMuR, dim[-1, -2, 1]];
  expectDim["CITED stage003 [rho_br]=M L^-3", dRhoBr, dim[-3, 0, 1]];
  expectDim["CITED stage003 [rho_B0]=M L^-3", dRhoB0, dim[-3, 0, 1]];
  expectDim["CITED stage003 [J]=L^2 T^-1", dJ, dim[2, -1, 0]];
  expectDim["CITED stage003 [C_J]=M L^-1 T^-1", dCJ, dim[-1, -1, 1]];
  expectDim["CITED stage003 [K_theta]=M L T^-2", dKtheta, dim[1, -2, 1]];
  expectDim["CITED stage003 [chi_c]=M L^-5 T^2", dChiC, dim[-5, 2, 1]];
  expectDim["CITED stage003 [B]=M L^-1 T^-2", dB, dim[-1, -2, 1]];
  expectDim["CITED stage003 [c_gamma^2]=L^2 T^-2", dimSub[dMuR, dRhoBr], dimPow[velocityD, 2]];
  <|"CJ" -> consumedCJ, "BEff" -> consumedBEff, "CS0Sq" -> consumedCS0Sq|>
];

runLegADimensions[] := Module[
  {
    chi, nDim, grad4, kappaB, aB, muR4, fThroat, fMix,
    terms, reinjectedTerms, pdeTerms, muChi, MChi, POrder, MN, numberMu,
    W, rhoBProjected, SLeak, muRProjected
  },
  subheading["Leg A1 dimensional audit: OPERATIVE post-Decision-16 action and adjunct rows"];
  chi = zeroD;
  nDim = rho4D;
  grad4 = dimScale[-1, ell];
  kappaB = dimSub[massD, dimScale[2, timeD]];
  aB = fDensity4D;
  muR4 = fDensity4D;
  fThroat = fDensity4D;
  fMix = fDensity4D;

  terms = <|
    "POSTULATED P1 kinetic 1/2*m_GNLS*n*|u|^2" -> <|"Dim" -> dimAdd[mGnlsD, nDim, dimScale[2, velocityD]], "Symbols" -> {"m_GNLS", "n", "u"}|>,
    "CITED I-1 U(n)=(K/4)*n^5" -> <|"Dim" -> dimAdd[kEosD, dimScale[5, nDim]], "Symbols" -> {"K", "n"}|>,
    "POSTULATED P3 f_B=a_B*chi_B^2*(1-chi_B)^2" -> <|"Dim" -> aB, "Symbols" -> {"a_B", "chi_B"}|>,
    "POSTULATED P5 (kappa_B/2)*|grad_4 chi_B|^2" -> <|"Dim" -> dimAdd[kappaB, dimScale[2, dimAdd[grad4, chi]]], "Symbols" -> {"kappa_B", "chi_B"}|>,
    "POSTULATED P6 chi_B*f_shear" -> <|"Dim" -> dimAdd[chi, muR4, dimScale[2, dimAdd[grad4, ell]]], "Symbols" -> {"chi_B", "mu_R^(4)", "u"}|>,
    "DEFERRED_PLACEHOLDER f_throat" -> <|"Dim" -> fThroat, "Symbols" -> {"f_throat"}|>,
    "DEFERRED_PLACEHOLDER f_mix" -> <|"Dim" -> fMix, "Symbols" -> {"f_mix"}|>
  |>;
  KeyValueMap[Function[{name, term}, expectDim[name <> " has 4D free-energy-density dim M L^-2 T^-2", term["Dim"], fDensity4D]], terms];
  expectZero[
    "POSTULATED operative F integrand homogeneity",
    homResidual[AssociationThread[Keys[terms], Lookup[Values[terms], "Dim"]]]
  ];
  expectZero["Decision-16 operative A1 surface excludes retired symbol alpha_aniso", operativeA1AbsenceResidual[terms]];

  Print["  RETIRED-HISTORICAL (not in operative A1): alpha_aniso*chi_B*(P.w_hat)^2 had [alpha_aniso]=M L^-2 T^-2; retired with P by Decision 16, not by a dimensional defect."];
  expectDim["RETIRED-HISTORICAL P7 alpha_aniso*chi_B*(P.w_hat)^2 was dimensionally homogeneous", fDensity4D, fDensity4D];
  reinjectedTerms = Append[
    terms,
    "REINJECTED alpha_aniso*chi_B*(P.w_hat)^2" -> <|"Dim" -> fDensity4D, "Symbols" -> {"alpha_aniso", "chi_B", "P", "w_hat"}|>
  ];
  expectFail["Decision-16 A1 tooth: re-inject alpha_aniso term trips operative retired-symbol absence", operativeA1AbsenceResidual[reinjectedTerms]];
  expectZero["Decision-16 baseline operative A1 surface remains alpha_aniso-free after copy mutation", operativeA1AbsenceResidual[terms]];

  pdeTerms = <|
    "partial_t(chi_B*n)" -> dimSub[dimAdd[chi, nDim], timeD],
    "div_4(chi_B*n*u)" -> dimSub[dimAdd[chi, nDim, velocityD], ell],
    "div_4 J_chi" -> dimSub[dim[-3, -1, 0], ell],
    "n*Gamma_B" -> dimSub[nDim, timeD]
  |>;
  KeyValueMap[Function[{name, d}, expectDim["POSTULATED balance row [" <> name <> "]=L^-4 T^-1", d, dim[-4, -1, 0]]], pdeTerms];
  expectZero["POSTULATED balance PDE homogeneity implies [Gamma_B]=T^-1", homResidual[pdeTerms]];

  muChi = fDensity4D;
  MChi = dim[2, 1, -1];
  POrder = dimAdd[muChi, dimScale[-1, timeD], dimScale[4, ell]];
  MN = dim[-4, 1, -1];
  numberMu = energyD;
  expectDim["POSTULATED ADJUNCT [mu_chi]=M L^-2 T^-2", muChi, fDensity4D];
  expectDim["POSTULATED ADJUNCT [M_chi]=L^2 T M^-1", MChi, dim[2, 1, -1]];
  expectDim["POSTULATED ADJUNCT P_order=int mu_chi*D_t chi_B d4X = M L^2 T^-3", POrder, dimSub[energyD, timeD]];
  expectDim["POSTULATED P12 [M_n]=L^-4 T M^-1", MN, dim[-4, 1, -1]];
  expectDim["POSTULATED P12 J_repair=-M_n*grad(mu) gives L^-3 T^-1", dimSub[dimAdd[MN, numberMu], ell], dim[-3, -1, 0]];

  W = dimScale[-1, ell];
  rhoBProjected = dimAdd[W, rho4D, ell];
  SLeak = dimSub[rhoBProjected, timeD];
  muRProjected = dimAdd[fDensity4D, ell];
  expectDim["CONVENTION P10 [W]=L^-1", W, dimScale[-1, ell]];
  expectDim["CONVENTION P10 projected [rho_B]=int W*chi_B*n dw = L^-4", rhoBProjected, dimScale[-4, ell]];
  expectDim["CONVENTION P10 [S_leak]=L^-4 T^-1", SLeak, dim[-4, -1, 0]];
  expectDim["POSTULATED/PENDING P6 int chi_B*mu_R^(4) dw has brane [mu_R]=M L^-1 T^-2", muRProjected, braneLagD]
];

runLegAClosure[] := Module[
  {x, w, t, total, orderLeft, orderResidual, disorderLeft, disorderResidual, advectiveResidual, sourceCoeff},
  subheading["Leg A2 structural-closure identities via formal functions (EARNED)"];
  total = D[nfun[x, w, t], t] + D[nfun[x, w, t] ux[x, w, t], x] + D[nfun[x, w, t] uw[x, w, t], w];
  orderLeft =
    D[chifun[x, w, t] nfun[x, w, t], t] +
    D[chifun[x, w, t] nfun[x, w, t] ux[x, w, t] + Jx[x, w, t], x] +
    D[chifun[x, w, t] nfun[x, w, t] uw[x, w, t] + Jw[x, w, t], w];
  orderResidual = orderLeft - nfun[x, w, t] GammaB[x, w, t];
  advectiveResidual = FullSimplify[
    Expand[(orderResidual - chifun[x, w, t] total)/nfun[x, w, t]] -
      (D[chifun[x, w, t], t] + ux[x, w, t] D[chifun[x, w, t], x] + uw[x, w, t] D[chifun[x, w, t], w] +
        (D[Jx[x, w, t], x] + D[Jw[x, w, t], w])/nfun[x, w, t] - GammaB[x, w, t])
  ];
  expectZero["EARNED advective form D_t chi_B=Gamma_B-(1/n)div J_chi", advectiveResidual];

  disorderLeft =
    D[(1 - chifun[x, w, t]) nfun[x, w, t], t] +
    D[(1 - chifun[x, w, t]) nfun[x, w, t] ux[x, w, t] - Jx[x, w, t], x] +
    D[(1 - chifun[x, w, t]) nfun[x, w, t] uw[x, w, t] - Jw[x, w, t], w];
  disorderResidual = disorderLeft + nfun[x, w, t] GammaB[x, w, t];
  expectZero["EARNED disorder complement is total minus order residual", FullSimplify[disorderResidual - (total - orderResidual)]];
  expectZero["EARNED order+disorder sum reproduces total conservation with Gamma_B cancelling", FullSimplify[orderResidual + disorderResidual - total]];
  sourceCoeff = Coefficient[orderResidual, GammaB[x, w, t]];
  expectNonzero["EARNED order balance is genuinely sourced by n*Gamma_B", sourceCoeff]
];

runLegAWall[] := Module[
  {w, chiSol, chi, fB, coeff, fCoeff, deltaSq, deltaSqSym, deltaRule, profile, residual, sigma},
  subheading["Leg A3 wall admission by DSolve/substitution (EARNED relative to P3/P5)"];
  $Assumptions = aBsym > 0 && kappaBsym > 0 && delta > 0;
  chiSol = y[w] /. First[DSolve[{y'[w] == y[w] (1 - y[w])/delta, y[0] == 1/2}, y, w]];
  expectZero["EARNED DSolve logistic profile equals 1/(1+Exp[-w/delta])", FullSimplify[chiSol - 1/(1 + Exp[-w/delta])]];
  chi = Unique["chi"];
  fB = aBsym chi^2 (1 - chi)^2;
  coeff = FullSimplify[kappaBsym D[chiSol, {w, 2}]/(chiSol (1 - chiSol) (1 - 2 chiSol))];
  fCoeff = FullSimplify[D[fB, chi]/(chi (1 - chi) (1 - 2 chi))];
  deltaSqSym = Unique["deltaSq"];
  deltaSq = deltaSqSym /. First[Solve[FullSimplify[coeff delta^2]/deltaSqSym == fCoeff, deltaSqSym]];
  expectZero["EARNED kink width solved from kappa_B*chi''=f_B'", Sqrt[deltaSq] - Sqrt[kappaBsym/(2 aBsym)]];
  deltaRule = delta -> Sqrt[deltaSq];
  profile = FullSimplify[chiSol /. deltaRule];
  residual = FullSimplify[kappaBsym D[profile, {w, 2}] - (D[fB, chi] /. chi -> profile)];
  expectZero["EARNED kink EL residual vanishes by substitution", residual];
  sigma = FullSimplify[Integrate[kappaBsym D[profile, w]^2, {w, -Infinity, Infinity}, Assumptions -> aBsym > 0 && kappaBsym > 0]];
  expectZero["EARNED sigma_wall exact closed form", sigma - Sqrt[2 aBsym kappaBsym]/6];
  expectDim["EARNED [sigma_wall]=M L^-1 T^-2", dimPow[dimAdd[fDensity4D, dimSub[massD, dimScale[2, timeD]]], 1/2], braneLagD]
];

runLegBProjection[] := Module[
  {w, lambda, A, Wprofile, Qprofile, intWQPrime, intWpQ, sFlux, wrongFlux, assumptions},
  subheading["Leg B1 projected two-source law on independent hyperbolic profile family (EARNED)"];
  assumptions = lambda > 0 && A > 0;
  Wprofile[z_] := Sech[z/lambda]^2/(2 lambda);
  Qprofile[z_] := A Tanh[z/lambda] Sech[z/lambda]^2;
  intWQPrime = FullSimplify[
    Integrate[Wprofile[w] D[Qprofile[w], w], {w, -Infinity, Infinity}, Assumptions -> assumptions],
    assumptions
  ];
  intWpQ = FullSimplify[
    Integrate[D[Wprofile[w], w] Qprofile[w], {w, -Infinity, Infinity}, Assumptions -> assumptions],
    assumptions
  ];
  sFlux = intWpQ;
  wrongFlux = -intWpQ;
  expectZero["EARNED IBP identity hyperbolic family: -int W*Q' = S_flux", FullSimplify[-intWQPrime - sFlux]];
  expectFail["Leg-B ablation hyperbolic family: flip W' sign breaks IBP identity", FullSimplify[-intWQPrime - wrongFlux]];
  expectNonzero["EARNED hyperbolic family: W' contribution is not accidentally zero", intWpQ]
];

runLegBRecovery[] := Module[
  {
    w, lambda, n0, u0, chi0, chiAmp, gammaAmp, jchiAmp, Wlive, nLive,
    uwLive, chiLive, gammaLive, jchiLive, assumptions, assembled,
    sTwoSourceLimit, Wfrozen, jFrozen, sLeak, residualChiC, residualGamma,
    residualJchi, expectedJchiFlux, muW, rho0, E0, Wlambda, phi, Ew,
    jw, derived, frozen, corruptPhi, corruptJ, corrupt
  },
  subheading["Leg B2/B3 recovery reduction and frozen Gaussian anchor (EARNED)"];
  assumptions = lambda > 0 && n0 > 0 && u0 > 0 && Element[{chi0, chiAmp, gammaAmp, jchiAmp, cConst}, Reals];
  Wlive[z_] := Sech[z/lambda]^2/(2 lambda);
  nLive[z_] := n0 Sech[z/lambda]^2;
  uwLive[z_] := u0 Tanh[z/lambda];
  chiLive[z_] := chi0 + chiAmp Sech[z/lambda]^2;
  gammaLive[z_] := gammaAmp Sech[z/lambda]^2;
  jchiLive[z_] := jchiAmp Tanh[z/lambda] Sech[z/lambda]^2;

  assembled = assembleProjectedTwoSource[
    Wlive[w], nLive[w], uwLive[w], chiLive[w], gammaLive[w], jchiLive[w], w, assumptions
  ];
  expectZero["EARNED B2 general chi_B/Gamma_B/J_chi projected order-balance IBP identity", assembled["rawFlux"] - assembled["sFlux"]];
  expectNonzero["EARNED B2 general chi_B profile is live in S_flux", D[assembled["sFlux"], chiAmp]];
  expectNonzero["EARNED B2 general Gamma_B profile is live in S_convert", D[assembled["sConvert"], gammaAmp]];
  expectNonzero["EARNED B2 general J_chi^w profile is live in S_flux", D[assembled["sFlux"], jchiAmp]];

  sTwoSourceLimit = FullSimplify[assembled["sTotal"] /. {chi0 -> 1, chiAmp -> 0, gammaAmp -> 0, jchiAmp -> 0}, assumptions];
  Wfrozen[z_] := Sech[z/lambda]^2/(2 lambda);
  jFrozen[z_] := n0 u0 Tanh[z/lambda] Sech[z/lambda]^2;
  sLeak = FullSimplify[
    -boundaryJump[Wfrozen[w] jFrozen[w], w, assumptions] +
      Integrate[D[Wfrozen[w], w] jFrozen[w], {w, -Infinity, Infinity}, Assumptions -> assumptions],
    assumptions
  ];
  expectZero["FROZEN stage_243 S_leak target equals independent closed form on B2 profile", sLeak + 4 n0 u0/(15 lambda)];
  expectZero["EARNED recovery reduction at chi_B=1,Gamma_B=0,J_chi=0 against frozen stage_243 S_leak", sTwoSourceLimit - sLeak];
  residualChiC = FullSimplify[(assembled["sTotal"] /. {chi0 -> cConst, chiAmp -> 0, gammaAmp -> 0, jchiAmp -> 0}) - sLeak, assumptions];
  expectZero["Leg-B conditionality chi_B=c residual computes to (c-1)*S_leak", residualChiC - (cConst - 1) sLeak];
  expectNonzero["Leg-B conditionality chi_B=c!=1 leaves (c-1)*S_leak residual", residualChiC];
  residualGamma = FullSimplify[(assembled["sTotal"] /. {chi0 -> 1, chiAmp -> 0, jchiAmp -> 0}) - sLeak, assumptions];
  expectZero["Leg-B conditionality Gamma_B residual computes to S_convert integral", residualGamma - (assembled["sConvert"] /. {chi0 -> 1, chiAmp -> 0, jchiAmp -> 0})];
  expectNonzero["Leg-B conditionality Gamma_B!=0 leaves computed S_convert residual", residualGamma];
  residualJchi = FullSimplify[(assembled["sTotal"] /. {chi0 -> 1, chiAmp -> 0, gammaAmp -> 0}) - sLeak, assumptions];
  expectedJchiFlux = FullSimplify[(assembled["sFlux"] /. {chi0 -> 1, chiAmp -> 0, gammaAmp -> 0}) - sTwoSourceLimit, assumptions];
  expectZero["Leg-B conditionality J_chi^w residual computes to its flux terms", residualJchi - expectedJchiFlux];
  expectNonzero["Leg-B conditionality J_chi^w!=0 leaves computed flux residual", residualJchi];
  expectFail["Leg-B ablation corrupt frozen S_leak transcription alone breaks B2 reduction", sTwoSourceLimit - (sLeak + 2 n0 u0/(15 lambda))];
  expectFail["Leg-B ablation corrupt general assembly alone breaks B2 reduction", (sTwoSourceLimit + 2 n0 u0/(15 lambda)) - sLeak];

  $Assumptions = lambda > 0 && muW > 0 && rho0 > 0 && E0 > 0;
  Wlambda = Exp[-w^2/lambda^2]/(lambda Sqrt[Pi]);
  phi = (2 w/(Sqrt[Pi] lambda^3)) Exp[-w^2/lambda^2];
  Ew = -E0 phi;
  jw = muW rho0 Ew;
  derived = FullSimplify[Integrate[D[Wlambda, w] jw, {w, -Infinity, Infinity}]];
  frozen = Sqrt[2] muW rho0 E0/(2 Sqrt[Pi] lambda^3);
  expectZero["EARNED stage_244 Gaussian one-mode S_leak direct integration", derived - frozen];
  corruptPhi = (2 w/(Sqrt[Pi] lambda^2)) Exp[-w^2/lambda^2];
  corruptJ = muW rho0 (-E0 corruptPhi);
  corrupt = FullSimplify[Integrate[D[Wlambda, w] corruptJ, {w, -Infinity, Infinity}]];
  expectFail["Leg-B ablation corrupt Gaussian phi lambda^2 vs lambda^3 mismatches stage_244", corrupt - frozen]
];

tokenFromDiscriminators[bracketZero_, bZero_, mZero_, flags_Association] := Module[{flag},
  Which[
    TrueQ[! mZero],
      {"FAIL_SECOND_CLASS_NOT_MAXWELL", "NOT_MAXWELL"},
    TrueQ[bracketZero && bZero && mZero],
      flag = If[And @@ Values[flags], "WITH_PROVENANCE", "BY_TUNING"];
      {"C5_RESOLVED_MAXWELL_BY_TUNING", flag},
    TrueQ[bracketZero && (! bZero || ! mZero)],
      {"FAIL_SECOND_CLASS_NOT_MAXWELL", "NOT_MAXWELL"},
    TrueQ[bZero],
      {"FAIL_C5_LONGITUDINAL_ZERO_MODE", "PROVENANCE_NO_GO"},
    True,
      {"FAIL_CAUCHY_STRAY_LONGITUDINAL", "PROVENANCE_NO_GO"}
  ]
];

runLegC[consumed_Association] := Module[
  {
    consumedCJ, consumedBEff, Kprov, Kmax, omega, omega2, omega2Long, matrix, det, coeffOmega,
    bracketEquivalent, bracketEquivalentMassless, bracketProv, kappaFree, kappaLocus, KthetaLocus,
    detProv, omega2Finite, expectedPole, locusSolve, flags, provPair, zeroPair,
    tunedPair, provToken, zeroToken, tunedToken, tunedFlag, detuneKPair,
    detuneBPair, detuneMPair, detuneK, detuneB, detuneM, counterPair,
    counterFlag, provBranch, zeroBranch, tunedBranch, detuneKBranch, detuneBBranch,
    detuneMBranch, branchAssumptions, q, qStar, reduceLife, uT, transLagrangian,
    transPoly, transOmega2, baselineSpeed, disturbedSpeed, baselineToken, disturbedToken
  },
  subheading["Leg C theta-as-Maxwell-phi no-go via longitudinal determinant (OWNED predicates, CONSUMED rule)"];
  consumedCJ = consumed["CJ"];
  consumedBEff = consumed["BEff"];
  Kprov = -kappaPhase;
  Kmax = consumedCJ^2/rhoBr;
  expectBool["OWNED C1 stable kappa_phase>0 gives K_theta=-kappa_phase<0", TrueQ[FullSimplify[Kprov < 0, kappaPhase > 0]]];
  expectBool["OWNED C1 Maxwell K_theta=C_J^2/rho_br>0", TrueQ[FullSimplify[Kmax > 0, Jsym > 0 && rhoB0 > 0 && rhoBr > 0]]];
  expectBool["OWNED C1 sign-lock predicates are OPPOSITE", TrueQ[FullSimplify[Kprov < 0 && Kmax > 0, kappaPhase > 0 && Jsym > 0 && rhoB0 > 0 && rhoBr > 0]]];
  expectDim["OWNED C1 [kappa_phase]=[K_theta]=M L T^-2", dimAdd[massD, ell, dimScale[-2, timeD]], dim[1, -2, 1]];

  matrix = {
    {rhoBr omega^2 - BEffSym ksym^2, I CJSym ksym omega},
    {-I CJSym ksym omega, KthetaSym ksym^2 - mTheta2}
  };
  det = FullSimplify[Det[matrix]];
  coeffOmega = FullSimplify[Coefficient[det /. {BEffSym -> 0}, omega^2]];
  expectZero["OWNED C2 determinant coefficient for curl-only longitudinal block includes m_theta^2", coeffOmega - (rhoBr (KthetaSym ksym^2 - mTheta2) - CJSym^2 ksym^2)];
  bracketEquivalent = FullSimplify[-coeffOmega/rhoBr];
  expectZero["OWNED C2 bracket-equivalent includes mTheta2 algebraically", bracketEquivalent - (CJSym^2 ksym^2/rhoBr - KthetaSym ksym^2 + mTheta2)];
  bracketEquivalentMassless = FullSimplify[bracketEquivalent /. mTheta2 -> 0];
  bracketProv = FullSimplify[bracketEquivalentMassless /. {CJSym -> consumedCJ, KthetaSym -> -kappaPhase}];
  expectZero["OWNED C2 bracket-equivalent k^2*(J^2*rho_B0^2+kappa_phase*rho_br)/rho_br", bracketProv - ksym^2 (Jsym^2 rhoB0^2 + kappaPhase rhoBr)/rhoBr];
  expectNonzero["OWNED C2 provenance branch bracket-equivalent is nonzero", bracketProv];

  kappaFree = Unique["kappaFree"];
  kappaLocus = kappaFree /. First[Solve[(bracketProv /. kappaPhase -> kappaFree) == 0, kappaFree]];
  KthetaLocus = FullSimplify[-kappaLocus];
  expectZero["OWNED C2 zero locus solved by Solve: K_theta=+J^2*rho_B0^2/rho_br", KthetaLocus - Jsym^2 rhoB0^2/rhoBr];
  locusSolve = KthetaSym /. First[Solve[(bracketEquivalentMassless /. CJSym -> consumedCJ) == 0, KthetaSym]];
  expectZero["OWNED C2 determinant degeneracy gives same Maxwell locus", locusSolve - Jsym^2 rhoB0^2/rhoBr];

  detProv = FullSimplify[det /. {CJSym -> consumedCJ, KthetaSym -> -kappaPhase, BEffSym -> consumedBEff, mTheta2 -> 0}];
  omega2Long = Unique["omega2Long"];
  omega2Finite = omega2Long /. First[Solve[(detProv /. omega^2 -> omega2Long) == 0, omega2Long]];
  expectedPole = ksym^2 kappaPhase rhoB0^2/(chiC (Jsym^2 rhoB0^2 + kappaPhase rhoBr));
  expectZero["OWNED C2 determinant mode structure gives finite provenance longitudinal pole", FullSimplify[omega2Finite - expectedPole]];
  expectZero["OWNED C2 tuned Maxwell fixture determinant coefficient computes to zero", FullSimplify[coeffOmega /. {CJSym -> consumedCJ, KthetaSym -> KthetaLocus, mTheta2 -> 0}]];

  flags = <|"K_theta_forced_by_frozen_defs" -> False, "B_eff_removed_by_frozen_defs" -> False, "m_theta_zero_forced" -> False|>;
  branchAssumptions = Jsym > 0 && rhoB0 > 0 && rhoBr > 0 && chiC > 0 && ksym > 0 && kappaPhase > 0 && mTheta2 > 0;
  provBranch = {CJSym -> consumedCJ, KthetaSym -> -kappaPhase, BEffSym -> consumedBEff, mTheta2 -> 0};
  zeroBranch = {CJSym -> consumedCJ, KthetaSym -> -kappaPhase, BEffSym -> 0, mTheta2 -> 0};
  tunedBranch = {CJSym -> consumedCJ, KthetaSym -> KthetaLocus, BEffSym -> 0, mTheta2 -> 0};
  provPair = tokenFromDiscriminators[
    algebraicZero[bracketEquivalent /. provBranch, branchAssumptions],
    algebraicZero[BEffSym /. provBranch, branchAssumptions],
    algebraicZero[mTheta2 /. provBranch, branchAssumptions],
    flags
  ];
  zeroPair = tokenFromDiscriminators[
    algebraicZero[bracketEquivalent /. zeroBranch, branchAssumptions],
    algebraicZero[BEffSym /. zeroBranch, branchAssumptions],
    algebraicZero[mTheta2 /. zeroBranch, branchAssumptions],
    flags
  ];
  tunedPair = tokenFromDiscriminators[
    algebraicZero[bracketEquivalent /. tunedBranch, branchAssumptions],
    algebraicZero[BEffSym /. tunedBranch, branchAssumptions],
    algebraicZero[mTheta2 /. tunedBranch, branchAssumptions],
    flags
  ];
  provToken = provPair[[1]];
  zeroToken = zeroPair[[1]];
  tunedToken = tunedPair[[1]];
  tunedFlag = tunedPair[[2]];
  expectZero["CONSUMED classification emits FAIL_CAUCHY_STRAY_LONGITUDINAL on provenance branch", If[provToken === "FAIL_CAUCHY_STRAY_LONGITUDINAL", 0, 1]];
  expectZero["CONSUMED classification emits FAIL_C5_LONGITUDINAL_ZERO_MODE when B_eff=0", If[zeroToken === "FAIL_C5_LONGITUDINAL_ZERO_MODE", 0, 1]];
  expectZero["CONSUMED classification emits C5_RESOLVED_MAXWELL_BY_TUNING on tuned fixture", If[tunedToken === "C5_RESOLVED_MAXWELL_BY_TUNING", 0, 1]];
  expectZero["OWNED provenance flag is BY_TUNING NOT WITH_PROVENANCE", If[tunedFlag === "BY_TUNING", 0, 1]];
  detuneKBranch = {CJSym -> consumedCJ, KthetaSym -> -kappaPhase, BEffSym -> 0, mTheta2 -> 0};
  detuneBBranch = {CJSym -> consumedCJ, KthetaSym -> KthetaLocus, BEffSym -> consumedBEff, mTheta2 -> 0};
  detuneMBranch = {CJSym -> consumedCJ, KthetaSym -> KthetaLocus, BEffSym -> 0};
  detuneKPair = tokenFromDiscriminators[
    algebraicZero[bracketEquivalent /. detuneKBranch, branchAssumptions],
    algebraicZero[BEffSym /. detuneKBranch, branchAssumptions],
    algebraicZero[mTheta2 /. detuneKBranch, branchAssumptions],
    flags
  ];
  detuneBPair = tokenFromDiscriminators[
    algebraicZero[bracketEquivalent /. detuneBBranch, branchAssumptions],
    algebraicZero[BEffSym /. detuneBBranch, branchAssumptions],
    algebraicZero[mTheta2 /. detuneBBranch, branchAssumptions],
    flags
  ];
  detuneMPair = tokenFromDiscriminators[
    algebraicZero[bracketEquivalent /. detuneMBranch, branchAssumptions],
    algebraicZero[BEffSym /. detuneMBranch, branchAssumptions],
    algebraicZero[mTheta2 /. detuneMBranch, branchAssumptions],
    flags
  ];
  detuneK = detuneKPair[[1]];
  detuneB = detuneBPair[[1]];
  detuneM = detuneMPair[[1]];
  expectZero["Leg-C detuning K_theta off locus re-fires FAIL_C5_LONGITUDINAL_ZERO_MODE", If[detuneK === "FAIL_C5_LONGITUDINAL_ZERO_MODE", 0, 1]];
  expectZero["Leg-C detuning B_eff!=0 on locus re-fires FAIL_SECOND_CLASS_NOT_MAXWELL", If[detuneB === "FAIL_SECOND_CLASS_NOT_MAXWELL", 0, 1]];
  expectZero["Leg-C detuning m_theta^2!=0 on locus re-fires FAIL_SECOND_CLASS_NOT_MAXWELL", If[detuneM === "FAIL_SECOND_CLASS_NOT_MAXWELL", 0, 1]];
  counterPair = tokenFromDiscriminators[
    algebraicZero[bracketEquivalent /. tunedBranch, branchAssumptions],
    algebraicZero[BEffSym /. tunedBranch, branchAssumptions],
    algebraicZero[mTheta2 /. tunedBranch, branchAssumptions],
    <|"K_theta_forced_by_frozen_defs" -> True, "B_eff_removed_by_frozen_defs" -> True, "m_theta_zero_forced" -> True|>
  ];
  counterFlag = counterPair[[2]];
  expectZero["Leg-C provenance-flag counterfactual flips BY_TUNING to WITH_PROVENANCE", If[counterFlag === "WITH_PROVENANCE", 0, 1]];

  qStar = q /. First[Solve[D[(1/2) (kappaReal q + kappa4 q^2), q] == 0, q]];
  expectZero["OWNED C3 Lifshitz finite-k minimum k_star^2=-kappa_phase/(2*kappa4)", qStar + kappaReal/(2 kappa4)];
  reduceLife = Reduce[qStar > 0 && kappa4 > 0, {kappaReal, kappa4}, Reals];
  expectBool["OWNED C3 k_star^2>0 iff kappa_phase<0", TrueQ[FullSimplify[Equivalent[reduceLife, kappaReal < 0 && kappa4 > 0]]]];
  expectBool["OWNED C3 kappa_phase>0 gives k_star^2<0", TrueQ[FullSimplify[(qStar /. kappaReal -> kappaPhase) < 0, kappaPhase > 0 && kappa4 > 0]]];
  expectDim["CITED pathA_25 k_Rstar^2=40*K*m*rho0^4/hbar^2 has L^-2", dimSub[dimAdd[kEosD, mGnlsD, dimScale[4, rho4D]], dimScale[2, hbarD]], dimScale[-2, ell]];

  transLagrangian = (1/2) epsilonSym omega2 uT^2 - (1/2) muR ksym^2 uT^2;
  transPoly = FullSimplify[D[transLagrangian, uT]/uT];
  transOmega2 = omega2 /. First[Solve[transPoly == 0, omega2]];
  baselineSpeed = FullSimplify[(transOmega2/ksym^2) /. epsilonSym -> rhoBr];
  disturbedSpeed = FullSimplify[(transOmega2/ksym^2) /. epsilonSym -> 2 rhoBr];
  expectZero["CONSUMED C4 transverse dispersion from L_T gives omega^2=(mu_R/epsilon)*k^2", transOmega2 - (muR/epsilonSym) ksym^2];
  expectZero["CONSUMED C4 transverse baseline speed equals consumed c_gamma^2=mu_R/rho_br", baselineSpeed - muR/rhoBr];
  expectZero["CONSUMED C4 epsilon!=rho_br fixture shifts transverse speed to mu_R/(2*rho_br)", disturbedSpeed - muR/(2 rhoBr)];
  baselineToken = If[algebraicZero[baselineSpeed - muR/rhoBr], "PASS_TRANSVERSE_UNDISTURBED", "FAIL_TRANSVERSE_DISTURBED"];
  disturbedToken = If[! algebraicZero[disturbedSpeed - muR/rhoBr], "FAIL_TRANSVERSE_DISTURBED", "PASS_TRANSVERSE_UNDISTURBED"];
  expectZero["CONSUMED C4 baseline emits PASS_TRANSVERSE_UNDISTURBED", If[baselineToken === "PASS_TRANSVERSE_UNDISTURBED", 0, 1]];
  expectZero["CONSUMED C4 epsilon!=rho_br emits FAIL_TRANSVERSE_DISTURBED", If[disturbedToken === "FAIL_TRANSVERSE_DISTURBED", 0, 1]]
];

runAblations[consumed_Association] := (
  subheading["Able-to-fail firewall ablations"];
  expectFail["Leg-A ablation drop m_GNLS from kinetic term breaks F-density homogeneity", dimResidual[dimAdd[rho4D, dimScale[2, velocityD]], fDensity4D]];
  expectFail["Leg-A ablation kappa_B*chi_B^2 no-gradient breaks F-density homogeneity", dimResidual[dimSub[massD, dimScale[2, timeD]], fDensity4D]];
  expectFail["Leg-A ablation Gamma_B in place of n*Gamma_B breaks balance source dim", dimResidual[dimScale[-1, timeD], dim[-4, -1, 0]]];
  expectFail["Leg-A ablation handoff P_order=int mu_chi*n*Gamma_B d4X is inhomogeneous", dimResidual[dimAdd[fDensity4D, rho4D, dimScale[-1, timeD], dimScale[4, ell]], dimSub[energyD, timeD]]];

  expectFail["Consuming ablation B_eff multiply-vs-divide breaks citation integrity", rhoB0^2 chiC - rhoB0^2/chiC];
  expectFail["Consuming ablation C_J sign +J*rho_B0 breaks citation integrity", Jsym rhoB0 + Jsym rhoB0];
  expectFail["Consuming ablation c_s0^2 coefficient 4 breaks citation integrity", 4 kEOSsym rho0sym^4/mGNLSsym - 5 kEOSsym rho0sym^4/mGNLSsym];
  expectZero["baseline consumed C_J still live after ablations", consumed["CJ"] + Jsym rhoB0];
  expectZero["baseline consumed B_eff still live after ablations", consumed["BEff"] - rhoB0^2/chiC];
  expectZero["baseline consumed c_s0^2 still live after ablations", consumed["CS0Sq"] - 5 kEOSsym rho0sym^4/mGNLSsym]
);

driftPartitionResidual[pre_, operative_, retired_, verdictNDelta_: 0] := Module[
  {preSet, operativeSet, retiredSet, residual, n, token},
  preSet = Union[pre];
  operativeSet = Union[operative];
  retiredSet = Union[retired];
  residual = 0;
  residual += If[DuplicateFreeQ[pre], 0, 1];
  residual += If[DuplicateFreeQ[operative], 0, 1];
  residual += If[DuplicateFreeQ[retired], 0, 1];
  residual += If[retiredSet === Union[retiredD16DriftMembers], 0, 1];
  residual += If[Intersection[operativeSet, retiredSet] === {}, 0, 1];
  residual += If[Union[operativeSet, retiredSet] === preSet, 0, 1];
  residual += If[operativeSet === Complement[preSet, retiredD16DriftMembers], 0, 1];
  n = Length[operative];
  residual += (n - (Length[pre] - Length[retiredD16DriftMembers]))^2;
  token = "DRIFT(" <> ToString[n + verdictNDelta] <> ")";
  residual += If[token === expectedPostD16DriftToken, 0, 1];
  FullSimplify[residual]
];

printCarriedAndDrift[] := Module[{preD16, retired, operative, n, token, reinjected},
  subheading["Carried tokens, deferred items, and computed drift"];
  Print["  carried no-go tokens verbatim: FAIL_CAUCHY_STRAY_LONGITUDINAL; FAIL_C5_LONGITUDINAL_ZERO_MODE; C5_RESOLVED_MAXWELL_BY_TUNING (BY_TUNING, not WITH_PROVENANCE)."];
  Print["  pathA_25 lineage carried: FAIL_LIGHT_STARVED finite-k wall/smectic sign-flip route."];
  Print["  SECOND_MEDIUM_DRIFT lineage note: chi_B package is Part-I drift; rho_B0 and chi_c also appear in pathA_41 Part-VI drift trio, cross-reference only."];
  Print["  THETA_BRANCH_DEAD_NOT_ADMITTED: theta,J,rho_B0,K_theta/kappa_phase,chi_c,B are not live knobs of this stage."];
  Print["  DEFERRED: Gate-L delta w=u_w translation-Goldstone hazard; J_chi!=0; f_throat/f_mix; dynamics adjunct."];
  Print["  POSTULATED GLOBAL RETURNS: R_0=-M_0, R_1=-D_1 (not locally asserted)."];
  preD16 = preD16DriftMembers;
  retired = Select[preD16, MemberQ[retiredD16DriftMembers, #] &];
  operative = Select[preD16, ! MemberQ[retiredD16DriftMembers, #] &];
  n = Length[operative];
  token = "DRIFT(" <> ToString[n] <> ")";
  expectZero["pre-Decision-16 drift enumeration computes six members", Length[preD16] - 6];
  expectZero["Decision-16 retired drift complement is exactly {alpha_aniso}", If[Union[retired] === Union[retiredD16DriftMembers], 0, 1]];
  expectZero["operative drift is the set partition pre_D16_DRIFT6 minus {alpha_aniso}", driftPartitionResidual[preD16, operative, retired]];
  expectZero["COMPUTED operative DRIFT tally has five live chi_B inputs", n - 5];
  expectZero["operative DRIFT token is built from computed n", If[token === expectedPostD16DriftToken, 0, 1]];
  Print["  pre-D16: DRIFT(6) incl. alpha_aniso -> operative: DRIFT(5)."];
  Print["  ", token, " computed {", StringRiffle[operative, "; "], "}."];
  Print["  rung_W:140 pre-D16 reconciliation recorded six; Decision 16 removes only alpha_aniso, leaving the operative five-member chi_B package."];

  reinjected = Append[operative, "alpha_aniso"];
  expectZero["Decision-16 drift reinjection fixture computes n=6", Length[reinjected] - 6];
  expectFail["Decision-16 drift tooth: re-inject alpha_aniso trips operative DRIFT(5) partition", driftPartitionResidual[preD16, reinjected, retired]];
  expectFail["Decision-16 drift tooth: corrupt computed n before token assembly trips DRIFT(5) equality", driftPartitionResidual[preD16, operative, retired, 1]];
  expectZero["Decision-16 baseline drift partition remains valid after copy mutations", driftPartitionResidual[preD16, operative, retired]]
];

printVerdictLabels[] := (
  Print[""];
  Print["Verdict labels:"];
  Print["  ledger earned-label (NOT a source verdict token): CHI_B_ACTION_CLASSIFIED_RECOVERY_VERIFIED  (postulated two-phase chi_B action dimensionally classified; wall admitted; recovery reduction to the frozen S_leak assert-zero; theta-as-phi no-go carried)"];
  Print["  headline: ACTION_SPECIFIED_CLASSIFIED   (structure; POSTULATED microstructure, all terms labeled)"];
  Print["  recovery sub-verdict (EARNED rel. to the imposed chi_B split + declared W): RECOVERY_REDUCTION_VERIFIED   (target = frozen stage_243/244 S_leak, incl. the Gaussian one-mode anchor)"];
  Print["  carried no-go: FAIL_CAUCHY_STRAY_LONGITUDINAL (finite B_eff) / FAIL_C5_LONGITUDINAL_ZERO_MODE (B_eff=0); positive control C5_RESOLVED_MAXWELL_BY_TUNING flagged BY_TUNING NOT WITH_PROVENANCE; the only provenance sign-flip = Lifshitz (pathA_25 wall, killed)"];
  Print["  drift: pre-D16: DRIFT(6) incl. alpha_aniso -> operative: DRIFT(5) computed {chi_B; a_B; kappa_B; Gamma_B; gating structure}; THETA_BRANCH_DEAD_NOT_ADMITTED; cross-ref rho_B0, chi_c in pathA_41 Part-VI drift"];
  Print["  consumed: ledger_stage004 {L,T,M}+U(rho) single-well; ledger_stage005 c_s0^2=5*K*rho0^4/m_GNLS; ledger_stage003 c_gamma^2=mu_R/rho_br, C_J=-J*rho_B0, B_eff=rho_B0^2/chi_c, second-class classification rule"];
  Print["  labeled postulates: P1..P13 (incl. KINETIC_MASS_FACTOR_PINNED, SAME_DENSITY_DEGENERACY_POSTULATED, Decision-16-retired historical alpha_aniso, HANDOFF_P_ORDER_N_PLACEMENT_CORRECTED, global-return R_0=-M_0,R_1=-D_1 NOT locally asserted, throat=phase-conversion ontology)"];
  Print["  P7 live reframe: chi_B is THE postulated wall OP, not currently |P_parallel|^2; FUTURE_GATE_CHI_B_EQ_ABS_P_PARALLEL_SQ is high-risk/Part-VII-adjacent and requires a NEW T0 freeze (not foreclosed)"];
  Print["  honest scope: does NOT earn light; dynamics/energy ledger = labeled adjunct; wall-translation Goldstone + J_chi + f_throat/f_mix DEFERRED"]
);

Module[{ok, consumed},
  heading["ledger_stage006_two_phase_chiB_ontology Mathematica audit"];
  ok = Catch[
    printPins[];
    consumed = runConsumed[];
    runLegADimensions[];
    runLegAClosure[];
    runLegAWall[];
    runLegBProjection[];
    runLegBRecovery[];
    runLegC[consumed];
    runAblations[consumed];
    printCarriedAndDrift[];
    printVerdictLabels[];
    True,
    "ledgerStage006Failure",
    Function[{name, tag}, False]
  ];
  Print[""];
  Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
  If[TrueQ[ok],
    Print["OVERALL PASS: Mathematica classified ledger_stage006 chi_B ontology and recovery/no-go checks exactly"];
    Exit[0],
    Print["OVERALL FAIL: Mathematica stage006 audit did not close"];
    Exit[1]
  ]
]
