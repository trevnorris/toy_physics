(* Ledger stage011 Mathematica audit: frozen-reduction certificate.

   Standalone, print-only, no arguments, no imports, no exports.  This is the
   Part-II pathA_30 II-G1a slice only: native frozen operator assembly,
   c_S^2 extraction, solved pinch-off domain, reduction certificate, c_S^2
   dimensional leg, and able-to-fail teeth.  Stage 012 owns D/N, DtN, Robin,
   pole ladder, and static-series work.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;

REDUCTIONCERTIFIED = "REDUCTION_CERTIFIED";
FAILDIMENSIONAL = "FAIL_DIMENSIONAL";
FAILOPERATORINTRUSION = "FAIL_OPERATOR_INTRUSION";
FAILWRONGSPEED = "FAIL_WRONG_SPEED";
FAILWRONGDOMAIN = "FAIL_WRONG_DOMAIN";
DNUNITTESTFAILDIMENSIONAL = "DN_UNITTEST_FAIL_DIMENSIONAL";

$Assumptions =
  L0 > 0 && Rmouth > 0 && K > 0 && rhoStar > 0 && m > 0 &&
  hbar > 0 && omega > 0 && Element[{s, rho, deltaWall, epsilonRho}, Reals] &&
  deltaWall != 0 && epsilonRho != 0;

cSSquaredBulk = FullSimplify[5 K rhoStar^4/m];
bulkN = FullSimplify[omega^2/cSSquaredBulk];
linearTaper = FullSimplify[Rmouth (1 - s/L0)];

raise[msg_] := Throw[msg, "ledgerStage011Failure"];

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
    raise[name]
  ]
];

fmt[expr_] := ToString[InputForm[Factor[Cancel[FullSimplify[expr]]]]];

expectZero[name_, residual_] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FAIL  ", name, ": residual = ", ToString[InputForm[clean]]];
    raise[name]
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
    raise[name]
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
    raise[name]
  ]
];

exprEqual[lhs_, rhs_: 0] := TrueQ[FullSimplify[lhs - rhs == 0]];
nonzeroQ[expr_] := ! TrueQ[FullSimplify[expr == 0]];
boolResidual[condition_] := If[TrueQ[condition], 0, 1];
verdictResidual[actual_, expected_] := If[actual === expected, 0, 1];
applyAssembledLs[operator_, trial_] := FullSimplify[
  operator /. {
    HoldPattern[Derivative[order_Integer][psiHat][s]] :> D[trial, {s, order}],
    psiHat[s] -> trial
  }
];

r1SiteFromExponent[exponent_] := FullSimplify[exponent K rho^(exponent - 1)/m];
r1EosSiteFromExponent[exponent_] := FullSimplify[D[K rho^exponent, rho]/m];

computeVerdict[dimensionalOk_, unsuppressedOperatorIntrusion_, operatorIsHelmholtz_, speedIsCs_, domainIsL0_] :=
  Which[
    ! TrueQ[dimensionalOk], FAILDIMENSIONAL,
    TrueQ[unsuppressedOperatorIntrusion], FAILOPERATORINTRUSION,
    ! TrueQ[operatorIsHelmholtz], FAILOPERATORINTRUSION,
    ! TrueQ[speedIsCs], FAILWRONGSPEED,
    ! TrueQ[domainIsL0], FAILWRONGDOMAIN,
    True, REDUCTIONCERTIFIED
  ];

xiEllCFirewallOk[xiExpr_] := TrueQ[xiSymbol =!= ellC] && FreeQ[xiExpr, ellC];

solveCapEndpoint[taperExpr_] := Module[{roots},
  roots = s /. Solve[taperExpr == 0, s];
  If[Length[roots] < 1, raise["taper has no symbolic cap root"]];
  FullSimplify[First[roots]]
];

buildReductionCase[
  sqrtGamma0Expr_,
  rho0Expr_,
  taperExpr_,
  bdgFlag_,
  bdgDeferredBySmallness_,
  deltaVConfExpr_: 0,
  speedFactor_: 1
] := Module[
  {
    measureCoeff, rho0grad, csSquaredLocal, csgrad, nCoeff, qBackground,
    qgrad, bFlag, bdgOperatorCoeff, ls, ideal, residual,
    operatorIsHelmholtz, scalarIntrusion, bdgIntrusion,
    unsuppressedOperatorIntrusion, psiCoeff, extractedSpeedSquared,
    speedIsCs, capEndpoint, domainIsL0, verdict
  },
  measureCoeff = FullSimplify[D[Log[sqrtGamma0Expr], s]];
  rho0grad = FullSimplify[D[rho0Expr, s]];
  csSquaredLocal = FullSimplify[5 K rho0Expr^4/m];
  csgrad = FullSimplify[D[Sqrt[csSquaredLocal], s]];
  nCoeff = FullSimplify[omega^2/(csSquaredLocal speedFactor^2)];
  qBackground = FullSimplify[-hbar^2 D[Sqrt[rho0Expr], {s, 2}]/(2 m Sqrt[rho0Expr])];
  qgrad = FullSimplify[D[qBackground, s]];
  bFlag = bdgFlag;
  bdgOperatorCoeff = FullSimplify[hbar^2/(4 m^2 cSSquaredBulk)];
  ls = FullSimplify[
    D[psiHat[s], {s, 2}] +
    measureCoeff D[psiHat[s], s] +
    nCoeff psiHat[s] -
    bFlag bdgOperatorCoeff D[psiHat[s], {s, 4}]
  ];
  ideal = FullSimplify[D[psiHat[s], {s, 2}] + bulkN psiHat[s]];
  residual = FullSimplify[ls - ideal];
  operatorIsHelmholtz = exprEqual[residual, 0];
  scalarIntrusion = nonzeroQ[deltaVConfExpr] || nonzeroQ[qgrad];
  bdgIntrusion = ! TrueQ[bFlag == 0] && ! TrueQ[bdgDeferredBySmallness];
  unsuppressedOperatorIntrusion =
    nonzeroQ[measureCoeff] || nonzeroQ[nCoeff - bulkN] || bdgIntrusion || scalarIntrusion;
  psiCoeff = FullSimplify[Coefficient[Expand[ls], psiHat[s]]];
  extractedSpeedSquared = FullSimplify[omega^2/psiCoeff];
  speedIsCs =
    exprEqual[psiCoeff, bulkN] &&
    exprEqual[extractedSpeedSquared, cSSquaredBulk] &&
    exprEqual[csgrad, 0];
  capEndpoint = solveCapEndpoint[taperExpr];
  domainIsL0 = exprEqual[capEndpoint, L0];
  verdict = computeVerdict[True, unsuppressedOperatorIntrusion, operatorIsHelmholtz, speedIsCs, domainIsL0];
  <|
    "SqrtGamma0" -> sqrtGamma0Expr,
    "Rho0" -> rho0Expr,
    "M" -> measureCoeff,
    "Rho0Grad" -> rho0grad,
    "CsSquaredLocal" -> csSquaredLocal,
    "CsGrad" -> csgrad,
    "N" -> nCoeff,
    "DeltaVConf" -> FullSimplify[deltaVConfExpr],
    "Q" -> qBackground,
    "QGrad" -> qgrad,
    "B" -> bFlag,
    "BdgOperatorCoeff" -> bdgOperatorCoeff,
    "L_s" -> ls,
    "Ideal" -> ideal,
    "OperatorResidual" -> residual,
    "OperatorIsHelmholtz" -> operatorIsHelmholtz,
    "UnsuppressedOperatorIntrusion" -> unsuppressedOperatorIntrusion,
    "PsiCoeff" -> psiCoeff,
    "ExtractedSpeedSquared" -> extractedSpeedSquared,
    "SpeedIsCs" -> speedIsCs,
    "Taper" -> taperExpr,
    "CapEndpoint" -> capEndpoint,
    "DomainIsL0" -> domainIsL0,
    "Verdict" -> verdict
  |>
];

dimResidualVec[actual_, expected_] := FullSimplify[(actual - expected).(actual - expected)];

dimOf[expr_, dims_] := Module[{clean, args, argDims, base, power},
  clean = FullSimplify[expr];
  Which[
    TrueQ[clean == 0] || NumericQ[clean], {0, 0, 0},
    AtomQ[clean] && KeyExistsQ[dims, clean], dims[clean],
    AtomQ[clean], raise["missing dimension for " <> ToString[Unevaluated[clean], InputForm]],
    Head[clean] === Times, Total[dimOf[#, dims] & /@ (List @@ clean)],
    Head[clean] === Power,
      base = clean[[1]];
      power = clean[[2]];
      If[! NumericQ[power], raise["non-numeric dimension power"]];
      power dimOf[base, dims],
    Head[clean] === Plus,
      args = Select[List @@ clean, ! TrueQ[FullSimplify[# == 0]] &];
      argDims = dimOf[#, dims] & /@ args;
      If[Length[argDims] == 0, {0, 0, 0},
        If[Length[DeleteDuplicates[argDims]] != 1, raise["dimension mismatch in sum"]];
        First[argDims]
      ],
    True, raise["unsupported dimension expression " <> ToString[clean, InputForm]]
  ]
];

buildDimensionalBlock[] := Module[
  {
    lengthDim, energyDim, fourVolumeDim, pressureDim, rhoDim, kDim,
    dimRules, expectedCsSquaredDim, csSquaredDim, dimensionalOk,
    corruptRules, corruptCsSquaredDim, corruptOk, probeVerdict,
    mutationFires, cleanVerdict, mutatedVerdict, failSuppressed
  },
  lengthDim = {1, 0, 0};
  energyDim = {2, 1, -2};
  fourVolumeDim = 4 lengthDim;
  pressureDim = energyDim - fourVolumeDim;
  rhoDim = -4 lengthDim;
  kDim = pressureDim - 5 rhoDim;
  dimRules = <|
    L0 -> lengthDim,
    omega -> {0, 0, -1},
    K -> kDim,
    rhoStar -> rhoDim,
    m -> {0, 1, 0}
  |>;
  expectedCsSquaredDim = {2, 0, -2};
  csSquaredDim = dimOf[cSSquaredBulk, dimRules];
  dimensionalOk = TrueQ[csSquaredDim == expectedCsSquaredDim];
  corruptRules = Join[KeyDrop[dimRules, K], <|K -> dimRules[K] + {1, 0, 0}|>];
  corruptCsSquaredDim = dimOf[cSSquaredBulk, corruptRules];
  corruptOk = TrueQ[corruptCsSquaredDim == expectedCsSquaredDim];
  probeVerdict = If[corruptOk, "NO_FAIL", DNUNITTESTFAILDIMENSIONAL];
  mutationFires = TrueQ[probeVerdict === DNUNITTESTFAILDIMENSIONAL];
  cleanVerdict = computeVerdict[dimensionalOk, False, True, True, True];
  mutatedVerdict = computeVerdict[corruptOk, False, True, True, True];
  failSuppressed = TrueQ[
    cleanVerdict === REDUCTIONCERTIFIED &&
    mutatedVerdict === FAILDIMENSIONAL &&
    mutationFires
  ];
  <|
    "LengthDim" -> lengthDim,
    "EnergyDim" -> energyDim,
    "FourVolumeDim" -> fourVolumeDim,
    "PressureDim" -> pressureDim,
    "RhoDim" -> rhoDim,
    "KDim" -> kDim,
    "CsSquaredDim" -> csSquaredDim,
    "ExpectedCsSquaredDim" -> expectedCsSquaredDim,
    "DimensionalOk" -> dimensionalOk,
    "CorruptKDim" -> corruptRules[K],
    "CorruptCsSquaredDim" -> corruptCsSquaredDim,
    "CorruptOk" -> corruptOk,
    "ProbeVerdict" -> probeVerdict,
    "MutationFires" -> mutationFires,
    "CleanVerdict" -> cleanVerdict,
    "MutatedVerdict" -> mutatedVerdict,
    "FailSuppressed" -> failSuppressed
  |>
];

buildBaseline[] := Module[
  {
    reduction, k, xi, bdgRatioDirect, bdgRatioWindow, siteA, siteB,
    consumed, dim, verdict, firewallOk
  },
  reduction = buildReductionCase[Aperp0, rhoStar, linearTaper, 0, True, D[0, s]];
  k = FullSimplify[omega/Sqrt[cSSquaredBulk]];
  xi = FullSimplify[hbar/(m Sqrt[cSSquaredBulk])];
  bdgRatioDirect = FullSimplify[hbar^2 k^2/(4 m^2 cSSquaredBulk)];
  bdgRatioWindow = FullSimplify[(k xi/2)^2];
  siteA = r1SiteFromExponent[5];
  siteB = r1EosSiteFromExponent[5];
  consumed = FullSimplify[siteA /. rho -> rhoStar];
  dim = buildDimensionalBlock[];
  verdict = computeVerdict[
    dim["DimensionalOk"],
    reduction["UnsuppressedOperatorIntrusion"],
    reduction["OperatorIsHelmholtz"],
    reduction["SpeedIsCs"],
    reduction["DomainIsL0"]
  ];
  reduction["Verdict"] = verdict;
  firewallOk = xiEllCFirewallOk[xi];
  <|
    "Reduction" -> reduction,
    "k" -> k,
    "Xi" -> xi,
    "BdgRatioDirect" -> bdgRatioDirect,
    "BdgRatioWindow" -> bdgRatioWindow,
    "SiteA" -> siteA,
    "SiteB" -> siteB,
    "Consumed" -> consumed,
    "Dim" -> dim,
    "Verdict" -> verdict,
    "FirewallOk" -> firewallOk
  |>
];

runAritySelfCheck[] := Module[{case, dim},
  subheading["Wolfram arity self-check"];
  case = buildReductionCase[Aperp0, rhoStar, linearTaper, 0, True];
  dim = buildDimensionalBlock[];
  expectBool["arity buildReductionCase[5 args] returns L_s and cap endpoint", KeyExistsQ[case, "L_s"] && KeyExistsQ[case, "CapEndpoint"]];
  expectBool["arity buildReductionCase optional delta/speed args accepts 7 args", KeyExistsQ[buildReductionCase[Aperp0, rhoStar, linearTaper, 0, True, 0, 1], "PsiCoeff"]];
  expectBool["arity exprEqual[lhs,rhs] handles two args", exprEqual[1, 1]];
  expectBool["arity solveCapEndpoint[taper] returns L0", exprEqual[solveCapEndpoint[linearTaper], L0]];
  expectBool["arity computeVerdict[5 booleans] returns REDUCTION_CERTIFIED", computeVerdict[True, False, True, True, True] === REDUCTIONCERTIFIED];
  expectBool["arity r1SiteFromExponent[e] returns exact literal site", exprEqual[r1SiteFromExponent[5], 5 K rho^4/m]];
  expectBool["arity r1EosSiteFromExponent[e] returns exact EOS route", exprEqual[r1EosSiteFromExponent[5], 5 K rho^4/m]];
  expectBool["arity buildDimensionalBlock[] returns mutation_fires", dim["MutationFires"] === True]
];

runReductionCertificate[data_] := Module[{red},
  red = data["Reduction"];
  subheading["Frozen background and reduction certificate"];
  Print["  POSTULATED geometry: straight finite throat, frozen eta=0, brane/mouth s=0, cap s=L0."];
  Print["  POSTULATED fields: rho0(s)=rho_star, A_M0=0, sqrt(gamma0)=A_perp0, matter perturbation exp(-I*omega*t)."];
  Print["  CITED speed: c_S^2 = 5*K*rho_star^4/m from R1 at rho_star; EOS exponent-5 P=K*rho^5 IMPOSED."];
  Print["  M=d_s log sqrt(gamma0) = ", fmt[red["M"]]];
  expectZero["projection measure coefficient M is computed zero", red["M"]];
  Print["  rho0grad=d_s rho_star = ", fmt[red["Rho0Grad"]]];
  expectZero["rho0grad is computed zero", red["Rho0Grad"]];
  Print["  c_s(s)^2 = ", fmt[red["CsSquaredLocal"]]];
  Print["  N=(omega/c_s(s))^2 = ", fmt[red["N"]]];
  expectZero["frozen N equals (omega/c_S)^2 with R1 c_S^2", red["N"] - bulkN];
  Print["  csgrad=d_s sqrt(c_s^2) = ", fmt[red["CsGrad"]]];
  expectZero["csgrad is computed zero", red["CsGrad"]];
  Print["  delta_V_conf witness = ", fmt[red["DeltaVConf"]], " (ell_c inert in frozen eta=0 test)"];
  expectZero["delta_V_conf witness is computed zero", red["DeltaVConf"]];
  Print["  Q(rho0) = ", fmt[red["Q"]], "; Qgrad=d_s Q = ", fmt[red["QGrad"]]];
  expectZero["Qgrad is computed by differentiation as zero", red["QGrad"]];
  Print["  BdG ratio direct = ", fmt[data["BdgRatioDirect"]]];
  Print["  BdG ratio window form (k*xi/2)^2 = ", fmt[data["BdgRatioWindow"]]];
  expectZero["BdG smallness witness hbar^2*k^2/(4*m^2*c_S^2) equals (k*xi/2)^2", data["BdgRatioDirect"] - data["BdgRatioWindow"]];
  expectZero["baseline BdG inclusion flag B=0 under k*xi<<1 deferral", red["B"]]
];

runOperatorAndSpeed[data_] := Module[{red, kOp, trigSinResidual, trigCosResidual},
  red = data["Reduction"];
  subheading["Produced operator, speed extraction, and domain solve"];
  Print["  assembled L_s = ", fmt[red["L_s"]]];
  Print["  ideal Helmholtz target = ", fmt[red["Ideal"]]];
  Print["  residual L_s - ideal = ", fmt[red["OperatorResidual"]]];
  expectZero["operator_is_helmholtz is produced by L_s assembly", boolResidual[red["OperatorIsHelmholtz"]]];
  expectZero["assembled L_s minus ideal Helmholtz operator is zero", red["OperatorResidual"]];
  expectBool["unsuppressed_operator_intrusion is computed false", red["UnsuppressedOperatorIntrusion"] === False];
  kOp = FullSimplify[omega/Sqrt[cSSquaredBulk]];
  trigSinResidual = applyAssembledLs[red["L_s"], Sin[kOp s]];
  trigCosResidual = applyAssembledLs[red["L_s"], Cos[kOp s]];
  expectZero["native trig-basis confirmation sin branch satisfies produced operator", trigSinResidual];
  expectZero["native trig-basis confirmation cos branch satisfies produced operator", trigCosResidual];
  Print["  ODE artifact: psi''(s) + (omega/c_S)^2 psi(s) = 0"];
  Print["  extracted psi_hat coefficient = ", fmt[red["PsiCoeff"]]];
  Print["  extracted speed squared omega^2/coeff = ", fmt[red["ExtractedSpeedSquared"]]];
  Print["  bulk c_S^2 trace = ", fmt[cSSquaredBulk]];
  expectZero["extracted psi_hat coefficient equals (omega/c_S)^2", red["PsiCoeff"] - bulkN];
  expectZero["extracted speed squared equals bulk R1 c_S^2", red["ExtractedSpeedSquared"] - cSSquaredBulk];
  expectBool["speed_is_cs is extracted from L_s and csgrad", red["SpeedIsCs"]];
  Print["  POSTULATED reference taper: R0(s)=R_mouth*(1-s/L0), monotone pinch with R0(0)>0 and R0(L0)=0."];
  Print["  solved cap endpoint solve(R0(s)=0,s) = ", fmt[red["CapEndpoint"]]];
  Print["  domain = [0, ", fmt[red["CapEndpoint"]], "]"];
  expectZero["domain_is_L0 is solved from the taper root", boolResidual[red["DomainIsL0"]]];
  expectZero["cap endpoint minus L0 is zero", red["CapEndpoint"] - L0]
];

runConsumedR1[data_] := (
  subheading["Consumed R1 input with dual-site integrity"];
  Print["  CITED ledger_stage005 R1 at rho_star; pathA_30 bare m is the stage004 m_GNLS ACTION primitive."];
  Print["  site A literal = ", fmt[data["SiteA"]]];
  Print["  site B EOS route d(K*rho^5)/d rho / m = ", fmt[data["SiteB"]]];
  expectZero["R1 site A minus site B equals zero", data["SiteA"] - data["SiteB"]];
  expectZero["R1 evaluated at rho_star equals c_S^2 bulk", data["Consumed"] - cSSquaredBulk];
  expectZero["explicit frozen-export anchor consumed - 5*K*rho_star^4/m equals zero", data["Consumed"] - 5 K rhoStar^4/m]
);

runDimensionalBlock[data_] := Module[{dim},
  dim = data["Dim"];
  subheading["Stage011 c_S^2 dimensional leg and corrupt-[K] probe"];
  Print["  dimension order: (L,M,T)"];
  Print["DIMENSIONS|axes=L,M,T"];
  Print["DIM|axes=L,M,T|name=LengthDim|exponents=", ToString[InputForm[dim["LengthDim"]]]];
  Print["DIM|axes=L,M,T|name=ExpectedCsSquaredDim|exponents=", ToString[InputForm[dim["ExpectedCsSquaredDim"]]]];
  Print["  [energy] = ", dim["EnergyDim"], "; [four-volume] = ", dim["FourVolumeDim"], "; [P] = ", dim["PressureDim"]];
  Print["  [rho] = ", dim["RhoDim"], "; [K]=[P]-5[rho] = ", dim["KDim"]];
  Print["  [c_S^2=5*K*rho_star^4/m] = ", dim["CsSquaredDim"]];
  expectZero["c_S^2 dimensional leg equals (2,0,-2)", dimResidualVec[dim["CsSquaredDim"], dim["ExpectedCsSquaredDim"]]];
  expectBool["dimensional_ok for the 011 c_S^2 leg", dim["DimensionalOk"]];
  Print["  corrupt [K]+(1,0,0) gives [K] = ", dim["CorruptKDim"]];
  Print["  corrupt [c_S^2] = ", dim["CorruptCsSquaredDim"], " -> ", dim["ProbeVerdict"]];
  expectZero["corrupt-[K] mutated c_S^2 dimension is exactly (3,0,-2)", dimResidualVec[dim["CorruptCsSquaredDim"], {3, 0, -2}]];
  expectBool["corrupt-[K] mutation_fires=True", dim["MutationFires"]];
  expectZero["self-ablation with mutation gives FAIL_DIMENSIONAL", verdictResidual[dim["MutatedVerdict"], FAILDIMENSIONAL]];
  expectZero["self-ablation without mutation gives REDUCTION_CERTIFIED", verdictResidual[dim["CleanVerdict"], REDUCTIONCERTIFIED]];
  expectBool["self-ablation fail_suppressed=True", dim["FailSuppressed"]]
];

runFirewallAndVerdict[data_] := Module[{red},
  red = data["Reduction"];
  subheading["xi != ell_c firewall and 011 verdict"];
  Print["  xi = hbar/(m*c_s) = ", fmt[data["Xi"]]];
  Print["  ell_c is the confinement length in V_wall(Sigma/ell_c), inert here because delta_V_conf=0."];
  expectBool["xi and ell_c are distinct symbols and never substituted", data["FirewallOk"]];
  Print["  011 scoped verdict = ", data["Verdict"]];
  expectZero["011 verdict is REDUCTION_CERTIFIED", verdictResidual[data["Verdict"], REDUCTIONCERTIFIED]];
  Print["  DN_UNITTEST_BC_DEPENDENT (JOINT) = (011: REDUCTION_CERTIFIED, computed here) AND (012: D/N ladder + bc_derivation_emitted=False -> BC_DEPENDENT landing, stage 012)"];
  expectBool["operator/speed/domain booleans are all computed true in baseline", red["OperatorIsHelmholtz"] && red["SpeedIsCs"] && red["DomainIsL0"]]
];

printProvenance[] := (
  subheading["Provenance and scope"];
  Print["  POSTULATED geometry: straight finite throat + frozen eta=0; L0 = ACTION-geometry throat depth, NOT a medium constant; R0(s) taper POSTULATED."];
  Print["  CITED-speed: c_S^2=5*K*rho_star^4/m is Part I edge R1 (stage005) evaluated at rho_star; EOS exponent-5 IMPOSED."];
  Print["  de-rig: operator PRODUCED by assembly, speed EXTRACTED from L_s, domain SOLVED from the pinch-off; replaces X==X/L0==L0/literal-True checks; the pinch-off domain is a labeled SELECTION."];
  Print["  validity-window: L_s is const-coeff Helmholtz conditional on {rho0'/rho0=0, sqrt(gamma0) const, delta_V_conf=0, grad_Q=0, k*xi<<1}; BdG k^4 is DEFERRED, not dropped unconditionally."];
  Print["  firewall: xi=hbar/(m*c_s) healing length is distinct from ell_c confinement length; ell_c is INERT here because delta_V_conf=0."];
  Print["  split: 011 carries the reduction-certificate component; D/N ladder + Robin + BC_DEPENDENT landing are stage 012; bc_derivation_emitted=False is 012's rung."];
  Print["  dropped-bookkeeping: scratch-YAML/_sympy_exprs.wl export, MMA-YAML re-read, expression_digest, and engine_agreement plumbing are stripped."];
  Print["  downstream consumers: stage 013 (harmonic beta lift) + stage 017 (calibration input) consume frozen L_s, domain, c_S, and validity window."];
  Print["  register note: zero new counted knobs; L0 is POSTULATED ACTION-geometry; ell_c INERT; xi DERIVED; validity record + firewall are structural edge candidates."]
);

printVerdictLabels[] := (
  Print[""];
  Print["Verdict labels:"];
  Print["  ledger earned-label (NOT a source verdict token): FROZEN_REDUCTION_HELMHOLTZ_CERTIFIED  (frozen wall eta=0 -> L_s assembled from the reduction (projection measure + every intruding coeff computed with its vanishing/deferral condition) -> const-coeff Helmholtz psi''+(omega/c_S)^2 psi=0; c_S^2=5*K*rho_star^4/m bulk; domain [0,L0] solved from the pinch-off R0(L0)=0; validity window {rho0'/rho0=0, sqrt(gamma0) const, delta_V_conf=0, grad_Q=0, k*xi<<1})"];
  Print["  source top-line verdict: DN_UNITTEST_BC_DEPENDENT  (JOINT; 011 carries the reduction-certificate component REDUCTION_CERTIFIED; the D/N ladder + bc_derivation_emitted=False -> BC_DEPENDENT landing = stage 012)"];
  Print["  joint composition: DN_UNITTEST_BC_DEPENDENT = (011: REDUCTION_CERTIFIED, computed here) AND (012: D/N pole ladder + Robin + BC_DEPENDENT landing, stage 012)"];
  Print["  earned (de-rig): operator_is_helmholtz PRODUCED by assembly (not X==X); speed_is_cs EXTRACTED from L_s (not literal True); domain_is_L0 SOLVED from R0(L0)=0 pinch-off (not L0==L0); unsuppressed_operator_intrusion COMPUTED; c_S^2 dim leg (2,0,-2) via [K]=[P]-5[rho] + corrupt-[K] probe fires"];
  Print["  postulated: straight finite throat, brane s=0, cap s=L0 (R0(L0)=0), rho0=rho_star, A_M0=0, frozen eta=0; L0 = ACTION-geometry throat depth (NOT a medium constant); R0(s) taper postulated"];
  Print["  cited (R1, stage005, dual-site integrity): c_S^2 = 5*K*rho^4/m evaluated at rho_star; EOS exponent-5 P=K*rho^5 IMPOSED"];
  Print["  deferred (validity window): BdG k^4 term hbar^2*k^4/(4*m^2), ratio hbar^2*k^2/(4*m^2*c_S^2)=(k*xi/2)^2, deferred only under k*xi<<1"];
  Print["  firewall: xi=hbar/(m*c_s) (healing length) != ell_c (confinement length) -- distinct symbols; ell_c INERT here (delta_V_conf=0)"]
);

runAbleToFailTeeth[data_] := Module[
  {
    baseline, measureMut, bdgMut, nonuniformRho0, rhoMut,
    deltaVConfMut,
    bdgFourthDerivative, baselineBdgFourthDerivativeCoeff,
    retainedBdgFourthDerivativeCoeff, speedMutVerdict,
    taperMut, domainMut, xiConflated,
    conflatedFirewallOk, dim
  },
  subheading["Able-to-fail mutation teeth"];
  baseline = data["Reduction"];

  measureMut = buildReductionCase[AperpMut[s], rhoStar, linearTaper, 0, True];
  expectFail["tooth 1 nonconstant sqrt(gamma0) makes M nonzero", measureMut["M"]];
  expectFail["tooth 1 operator_is_helmholtz boolean flips false", boolResidual[measureMut["OperatorIsHelmholtz"]]];
  expectZero["tooth 1 verdict is FAIL_OPERATOR_INTRUSION", verdictResidual[measureMut["Verdict"], FAILOPERATORINTRUSION]];

  bdgMut = buildReductionCase[Aperp0, rhoStar, linearTaper, 1, False];
  bdgFourthDerivative = D[psiHat[s], {s, 4}];
  baselineBdgFourthDerivativeCoeff = FullSimplify[Coefficient[Expand[baseline["L_s"]], bdgFourthDerivative]];
  retainedBdgFourthDerivativeCoeff = FullSimplify[Coefficient[Expand[bdgMut["L_s"]], bdgFourthDerivative]];
  expectZero[
    "tooth 2 deferred BdG flag leaves fourth-derivative term absent in baseline",
    baselineBdgFourthDerivativeCoeff
  ];
  expectNonzero[
    "tooth 2 retained BdG flag injects fourth-derivative term into L_s",
    retainedBdgFourthDerivativeCoeff
  ];
  expectFail["tooth 2 operator_is_helmholtz boolean flips false", boolResidual[bdgMut["OperatorIsHelmholtz"]]];
  expectZero["tooth 2 verdict is FAIL_OPERATOR_INTRUSION", verdictResidual[bdgMut["Verdict"], FAILOPERATORINTRUSION]];

  nonuniformRho0 = FullSimplify[rhoStar (1 + epsilonRho s/L0)];
  rhoMut = buildReductionCase[Aperp0, nonuniformRho0, linearTaper, 0, True];
  expectFail["tooth 3 nonuniform rho0 makes N-(omega/c_S)^2 nonzero", rhoMut["N"] - bulkN];
  expectFail["tooth 3 operator_is_helmholtz boolean flips false", boolResidual[rhoMut["OperatorIsHelmholtz"]]];
  expectZero["tooth 3 verdict is FAIL_OPERATOR_INTRUSION", verdictResidual[rhoMut["Verdict"], FAILOPERATORINTRUSION]];

  deltaVConfMut = buildReductionCase[Aperp0, rhoStar, linearTaper, 0, True, deltaWall];
  expectZero[
    "tooth 3b verdict is FAIL_OPERATOR_INTRUSION",
    verdictResidual[deltaVConfMut["Verdict"], FAILOPERATORINTRUSION]
  ];
  expectZero[
    "tooth 3b nonzero delta_V_conf witness makes unsuppressed_operator_intrusion true",
    boolResidual[deltaVConfMut["UnsuppressedOperatorIntrusion"]]
  ];
  expectZero[
    "tooth 3b operator_is_helmholtz stays true because the witness is not in L_s",
    boolResidual[deltaVConfMut["OperatorIsHelmholtz"]]
  ];

  Print["  note: FAIL_WRONG_SPEED is a defensive verdict branch not reachable by real operator corruption; intrusion dominates."];
  speedMutVerdict = computeVerdict[
    True,
    baseline["UnsuppressedOperatorIntrusion"],
    baseline["OperatorIsHelmholtz"],
    False,
    baseline["DomainIsL0"]
  ];
  expectZero[
    "compute_verdict logic branch returns FAIL_WRONG_SPEED when only speed_is_cs is false",
    verdictResidual[speedMutVerdict, FAILWRONGSPEED]
  ];

  taperMut = FullSimplify[Rmouth (1 - s/(2 L0))];
  domainMut = buildReductionCase[Aperp0, rhoStar, taperMut, 0, True];
  expectFail["tooth 5 corrupt taper root differs from L0", domainMut["CapEndpoint"] - L0];
  expectFail["tooth 5 domain_is_L0 boolean flips false", boolResidual[domainMut["DomainIsL0"]]];
  expectZero["tooth 5 verdict is FAIL_WRONG_DOMAIN", verdictResidual[domainMut["Verdict"], FAILWRONGDOMAIN]];

  xiConflated = FullSimplify[data["Xi"] /. hbar -> ellC hbar];
  conflatedFirewallOk = xiEllCFirewallOk[xiConflated];
  expectFail["tooth 6 ell_c -> xi conflation trips distinct-symbol firewall", boolResidual[conflatedFirewallOk]];

  Scan[
    Function[exponent,
      expectFail[
        "tooth 7 site A exponent 5->" <> ToString[exponent] <> " trips R1 dual-site integrity",
        r1SiteFromExponent[exponent] - data["SiteB"]
      ];
      expectFail[
        "tooth 7 site B exponent 5->" <> ToString[exponent] <> " trips R1 dual-site integrity",
        data["SiteA"] - r1EosSiteFromExponent[exponent]
      ]
    ],
    {4, 6}
  ];
  expectFail[
    "tooth 7 coordinated both-site exponent drift still trips frozen-export anchor",
    (r1SiteFromExponent[6] /. rho -> rhoStar) - 5 K rhoStar^4/m
  ];

  dim = data["Dim"];
  expectFail[
    "tooth 8 corrupt-[K] probe trips dimensional gate",
    dimResidualVec[dim["CorruptCsSquaredDim"], dim["ExpectedCsSquaredDim"]]
  ];
  expectZero["tooth 8 corrupt-[K] verdict is FAIL_DIMENSIONAL", verdictResidual[dim["MutatedVerdict"], FAILDIMENSIONAL]];
  expectBool["tooth 8 self-ablation fail_suppressed remains true", dim["FailSuppressed"]];

  expectZero["baseline immutable after teeth: operator residual remains zero", baseline["OperatorResidual"]];
  expectZero["baseline immutable after teeth: speed extraction remains bulk", baseline["ExtractedSpeedSquared"] - cSSquaredBulk];
  expectZero["baseline immutable after teeth: cap endpoint remains L0", baseline["CapEndpoint"] - L0];
  expectZero["baseline immutable after teeth: R1 site integrity remains zero", data["SiteA"] - data["SiteB"]];
  expectZero["baseline immutable after teeth: clean 011 verdict remains REDUCTION_CERTIFIED", verdictResidual[data["Verdict"], REDUCTIONCERTIFIED]]
];

Module[{ok, data},
  heading["ledger_stage011_frozen_reduction_certificate Mathematica audit"];
  ok = Catch[
    runAritySelfCheck[];
    data = buildBaseline[];
    assertExact["baseline", data];
    runReductionCertificate[data];
    runOperatorAndSpeed[data];
    runConsumedR1[data];
    runDimensionalBlock[data];
    runFirewallAndVerdict[data];
    printProvenance[];
    printVerdictLabels[];
    runAbleToFailTeeth[data];
    True,
    "ledgerStage011Failure",
    Function[{msg, tag}, Print["FAIL: ", msg]; False]
  ];

  Print[""];
  Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
  If[TrueQ[ok],
    Print["OVERALL PASS: Mathematica verified ledger_stage011 frozen reduction certificate exactly"];
    Exit[0],
    Print["OVERALL FAIL: Mathematica stage011 audit did not close"];
    Exit[1]
  ]
]
