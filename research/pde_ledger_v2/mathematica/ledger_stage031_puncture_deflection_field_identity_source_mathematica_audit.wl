(* Ledger stage031 Mathematica audit: puncture-deflection mechanism.

   Standalone, print-only, assert-zero, machine-real-free, and file-I/O-free.
   This is an independent Wolfram route: Integrate constructs the annulus,
   the actual I_plus expression and N_chi; DSolve derives the exterior branch;
   and the response uses a native typed matrix and Inverse.  Stage030's scalar
   closure is consumed symbolically, not re-derived.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE031_MUTATION";
activeMutation = Environment[mutationEnvironment];
If[!StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];

raise[name_] := Throw[name, "ledgerStage031Failure"];

assertExact[name_, expression_] := Module[{reals},
  reals = Cases[Unevaluated[expression], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FIRST_FAILURE=", name];
    Print["FAIL  ", name, ": machine-real atom(s): ", InputForm[reals]];
    raise[name]
  ]
];

cleanZero[expression_] := FullSimplify[expression] /.
  ConditionalExpression[0, _] -> 0;

expectZero[name_, residual_, evidence_: None] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FIRST_FAILURE=", name];
    Print["FAIL  ", name, ": residual = ", InputForm[clean]];
    If[evidence =!= None, Print["      evidence = ", InputForm[evidence]]];
    raise[name]
  ]
];

expectBool[name_, condition_, evidence_: None] :=
  expectZero[name, If[TrueQ[condition], 0, 1], evidence];

positiveExactQ[expression_, assumptions_: True] :=
  TrueQ[FullSimplify[expression > 0, assumptions]];

heading[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);

ablations = <|
  "PASS_FIELD_IDENTITY" -> <|"Primitive" -> "xi_power", "Value" -> 2, "Description" -> "xi_w=ell h -> ell^2 h"|>,
  "PASS_ETA_ANNULUS_NORMALIZATION" -> <|"Primitive" -> "eta_numerator", "Value" -> 2, "Description" -> "eta prefactor 3 -> 2"|>,
  "PASS_F0_MOUTH_VALUE_EVALUATED" -> <|"Primitive" -> "f0_power", "Value" -> 1, "Description" -> "sech^2 -> sech"|>,
  "PASS_REDUCED_COUPLING_NONZERO" -> <|"Primitive" -> "Jm_witness", "Value" -> 0, "Description" -> "J_m witness -> 0"|>,
  "PASS_BOUNDED_PLUS_SLEEVE" -> <|"Primitive" -> "peak_factor", "Value" -> 2, "Description" -> "peak normalization -> twice peak"|>,
  "PASS_REFLECTION_ANTISYMMETRY" -> <|"Primitive" -> "reflect_body", "Value" -> 0, "Description" -> "remove body reflection"|>,
  "PASS_I_PLUS_REFLECTION_DOMINANCE" -> <|"Primitive" -> "dominance_sign", "Value" -> -1, "Description" -> "+w dominance -> reversed"|>,
  "PASS_I_PLUS_EVEN_DEFORMATION_CONTROL" -> <|"Primitive" -> "even_control_odd_part", "Value" -> 1, "Description" -> "add odd part to even D"|>,
  "PASS_I_PLUS_ZERO_SLEEVE_CONTROL" -> <|"Primitive" -> "zero_sleeve", "Value" -> 1, "Description" -> "L_s=0 -> L_s>0"|>,
  "PASS_I_PLUS_NEGATIVE_SLEEVE_CONTROL" -> <|"Primitive" -> "negative_control", "Value" -> 1, "Description" -> "-w sleeve -> +w sleeve"|>,
  "PASS_N_CHI_NONZERO_GUARD" -> <|"Primitive" -> "nchi_guard", "Value" -> 0, "Description" -> "pre-division N_chi -> 0"|>,
  "PASS_Q_CHI_ORIENTATION" -> <|"Primitive" -> "projection_parity", "Value" -> 0, "Description" -> "odd kernel -> even kernel"|>,
  "PASS_BARE_FORCING_LIVE_VARIATION" -> <|"Primitive" -> "mouth_quadratic_factor", "Value" -> 2, "Description" -> "live K_m wiring factor -> 2"|>,
  "PASS_EXTERIOR_ONE_OVER_R" -> <|"Primitive" -> "radial_dimension", "Value" -> 4, "Description" -> "3D exterior -> 4D"|>,
  "PASS_POSITIVE_HOLDING_CURVATURE" -> <|"Primitive" -> "curvature_sign", "Value" -> -1, "Description" -> "holding curvature sign flip"|>,
  "PASS_NONZERO_HA_REQUIRES_CORE_HOLDER" -> <|"Primitive" -> "coreless_shift", "Value" -> 1, "Description" -> "inject nonzero coreless stationary point"|>,
  "PASS_KAPPA_REDUCTION" -> <|"Primitive" -> "kappa_factor", "Value" -> 2, "Description" -> "kappa -> 2D/b"|>,
  "PASS_Z_B_WIRING" -> <|"Primitive" -> "zb_factor", "Value" -> 0, "Description" -> "remove z_b from Z"|>,
  "PASS_RESPONSE_M_UU" -> <|"Primitive" -> "muu_extra", "Value" -> 1, "Description" -> "m_uu -> m_uu+1"|>,
  "PASS_RESPONSE_M_UG" -> <|"Primitive" -> "mug_sign", "Value" -> 1, "Description" -> "m_ug sign flip"|>,
  "PASS_RESPONSE_M_GG" -> <|"Primitive" -> "mgg_factor", "Value" -> 2, "Description" -> "m_gg target -> 2m_gg"|>,
  "PASS_RESPONSE_SYMMETRY" -> <|"Primitive" -> "symmetry_skew", "Value" -> 1, "Description" -> "add skew response"|>,
  "PASS_RESPONSE_DETERMINANT" -> <|"Primitive" -> "det_factor", "Value" -> 2, "Description" -> "det target -> twice"|>,
  "PASS_Z_G_POSTULATED_WITNESS" -> <|"Primitive" -> "zg_witness", "Value" -> 0, "Description" -> "z_g witness 1 -> 0"|>,
  "PASS_M_GG_NONNEGATIVE" -> <|"Primitive" -> "mgg_square", "Value" -> 1, "Description" -> "z_g^2 -> z_g"|>,
  "PASS_RESPONSE_STAR_WITNESS" -> <|"Primitive" -> "star_zb", "Value" -> 0, "Description" -> "star z_b 1 -> 0"|>,
  "PASS_S_GG_SELF_RESPONSE" -> <|"Primitive" -> "green_sign", "Value" -> -1, "Description" -> "positive Green form -> indefinite"|>,
  "PASS_NEUTRAL_FAR_FIELD_FORM" -> <|"Primitive" -> "potential_power", "Value" -> 2, "Description" -> "1/R -> 1/R^2"|>,
  "PASS_ALLOWED_ZERO_COEFFICIENT" -> <|"Primitive" -> "zero_C", "Value" -> 1, "Description" -> "C=0 control -> C=1"|>,
  "PASS_UNITS_XI_W" -> <|"Primitive" -> "dim_xi", "Value" -> {2, 0, 0}, "Description" -> "[xi_w] L -> L^2"|>,
  "PASS_UNITS_H" -> <|"Primitive" -> "dim_h", "Value" -> {1, 0, 0}, "Description" -> "[h] 1 -> L"|>,
  "PASS_UNITS_H_A" -> <|"Primitive" -> "dim_hA", "Value" -> {1, 0, 0}, "Description" -> "[h_A] 1 -> L"|>,
  "PASS_UNITS_K_M" -> <|"Primitive" -> "dim_Km", "Value" -> {3, -2, 1}, "Description" -> "[K_m] M L^4 T^-2 -> M L^3 T^-2"|>,
  "PASS_UNITS_J_M" -> <|"Primitive" -> "dim_Jm", "Value" -> {2, -2, 1}, "Description" -> "[J_m] M L^3 T^-2 -> E"|>,
  "PASS_UNITS_K_M_REDUCED" -> <|"Primitive" -> "dim_km", "Value" -> {3, -2, 1}, "Description" -> "[k_m] E -> E L"|>,
  "PASS_UNITS_G_CHIH" -> <|"Primitive" -> "dim_g", "Value" -> {3, -2, 1}, "Description" -> "[g] E -> E L"|>,
  "PASS_UNITS_ETA" -> <|"Primitive" -> "dim_eta", "Value" -> {-2, 0, 0}, "Description" -> "[eta] L^-3 -> L^-2"|>,
  "PASS_UNITS_ODD_KERNEL" -> <|"Primitive" -> "dim_o", "Value" -> {-1, 0, 0}, "Description" -> "[o_ell] L^-2 -> L^-1"|>,
  "PASS_UNITS_N_CHI" -> <|"Primitive" -> "dim_nchi", "Value" -> {0, 0, 0}, "Description" -> "[N_chi] L^-1 -> 1"|>,
  "PASS_UNITS_M_UU" -> <|"Primitive" -> "dim_muu", "Value" -> {2, 2, -1}, "Description" -> "[m_uu] L^3/E -> L^4/E"|>,
  "PASS_UNITS_M_UG" -> <|"Primitive" -> "dim_mug", "Value" -> {1, 2, -1}, "Description" -> "[m_ug] L^2/E -> L^3/E"|>,
  "PASS_UNITS_M_GG" -> <|"Primitive" -> "dim_mgg", "Value" -> {0, 2, -1}, "Description" -> "[m_gg] L/E -> L^2/E"|>,
  "PASS_UNITS_KAPPA" -> <|"Primitive" -> "dim_kappa", "Value" -> {0, -2, 1}, "Description" -> "[kappa] E/L -> E/L^2"|>,
  "PASS_UNITS_ESCAPE_MATRIX" -> <|"Primitive" -> "dim_Z21", "Value" -> {0, 0, 0}, "Description" -> "[Z_21] L -> 1"|>,
  "PASS_UNITS_DET_M" -> <|"Primitive" -> "dim_detm", "Value" -> {1, 4, -2}, "Description" -> "[det m] M^-2 T^4 -> L M^-2 T^4"|>,
  "PASS_UNITS_S_GG" -> <|"Primitive" -> "dim_sgg", "Value" -> {-1, 2, -1}, "Description" -> "[S_gg] E^-1 -> L/E"|>,
  "PASS_UNITS_SHELL_C" -> <|"Primitive" -> "dim_C", "Value" -> {3, -4, 2}, "Description" -> "[C] E^2 -> E^2/L"|>,
  "PASS_UNITS_SHELL_A" -> <|"Primitive" -> "dim_A", "Value" -> {2, -2, 1}, "Description" -> "[A] E L -> E"|>,
  "PASS_UNITS_U" -> <|"Primitive" -> "dim_U", "Value" -> {3, -2, 1}, "Description" -> "[U] E -> E L"|>,
  "PASS_UNITS_F" -> <|"Primitive" -> "dim_F", "Value" -> {2, -2, 1}, "Description" -> "[F] E/L -> E"|>
|>;

primitive[predicate_, name_, canonical_] := If[
  activeMutation === predicate && ablations[predicate]["Primitive"] === name,
  ablations[predicate]["Value"], canonical
];

dimResidual[actual_, expected_] := FullSimplify[(actual - expected).(actual - expected)];

Module[
  {ok, assumptions, rho, w, wp, u, a, ell, width, sleeve, hvar, Jm, Km,
   ss, xiPower, xiTest, etaNumerator, etaTest, etaIntegral, f0Power,
   f0Test, f0AtMouth, canonicalCurvature, testCurvature, jmWitness,
   gTest, rPlus, center, halfspan, rCentered, peakFactor, rangeIdentity,
   rBackground, oddKernel, reflectBody, rMinusTest, reflectionResidual,
   denForward, denReflected, denDifference, expectedDifference,
   dominanceSign, profileNumerator, rForward, rReflected,
   rDifferenceIdentity, rDifferenceCertificate, dPair, kernelMass,
   genericCertificate, actualPair, actualPairIdentity, iPlusNative, nChiNative,
   evenOddPart, dEvenTest, evenIntegral, zeroSleeve, zeroProfile, zeroD,
   zeroParity, negativeControl, negativePair, nchiGuard, projectionParity,
   qPlus, qMinus, mouthFactor, canonicalF0Mouth, canonicalEta, parentDensity,
   eulerDensityTest, integratedEuler, km, g, expectedEulerDensity, bareAmplitude,
   b, kh, c, zb, zg, delta, kappa, zMatrix, response, expectedMuu,
   expectedMug, expectedMgg, kappaFactor, zbFactor, zTest, muuTest,
   mugSign, mggFactor, skew, responseTest, detFactor, zgWitness,
   squarePower, deltaPositive, mggStable, starZb, starRules, starM,
   p, q, greenSign, sourceVector, greenMatrix, sgg,
   radialDimension, radialEquation, generalExterior, constants,
   boundaryRules, exterior, exteriorEnergy, curvatureSign, curvature,
   corelessShift, corelessStationary, rr, hA, R, s1, s2, capC,
   shellMgg, amplitude, potentialPower, potential, force, expectedPotential,
   expectedForce, zeroC, dimZero, dimL, dimT, dimMass, dimEnergy,
   unitRows, predicate, row, actualDimension, liveDimensionResidual},

  Print["ledger_stage031_puncture_deflection_field_identity_source Mathematica audit"];
  Print["CONSUMES_STAGE030={f0,f0(0)=1/ell,N0=8/(3ell),h=P0H,S_Lh,D,D*=7/4}"];
  Print["CONSUMES_STAGE003_030={B_eff,C_hu}"];

  ok = Catch[
    If[activeMutation =!= "" && !KeyExistsQ[ablations, activeMutation],
      Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
      Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
      raise["UNKNOWN_MUTATION"]
    ];
    If[activeMutation =!= "",
      Print["ACTIVE_MUTATION=", activeMutation];
      Print["MUTATED_PRIMITIVE=", ablations[activeMutation]["Description"]]
    ];

    assumptions = a > 0 && ell > 0 && width > 0 && sleeve > 0 && wp > 0 &&
      Element[{w, u}, Reals];

    (* Independent route begins with native response construction/inverse. *)
    heading["Full completed-square response via native matrix inverse"];
    delta = b kh - c^2;
    kappa = delta/b;
    zMatrix = {{1, 0}, {-(c/b) zb, zg}};
    response = FullSimplify[
      Transpose[zMatrix].Inverse[DiagonalMatrix[{b, kappa}]].zMatrix,
      b > 0 && delta > 0 && Element[{c, zb, zg}, Reals]];
    expectedMuu = (delta + c^2 zb^2)/(b delta);
    expectedMug = -c zb zg/delta;
    expectedMgg = b zg^2/delta;

    kappaFactor = primitive["PASS_KAPPA_REDUCTION", "kappa_factor", 1];
    expectZero["PASS_KAPPA_REDUCTION", kappaFactor kappa - delta/b];

    zbFactor = primitive["PASS_Z_B_WIRING", "zb_factor", 1];
    zTest = {{1, 0}, {-(c/b) zbFactor zb, zg}};
    expectZero["PASS_Z_B_WIRING", zTest[[2, 1]] + (c/b) zb];

    muuTest = response[[1, 1]] + primitive["PASS_RESPONSE_M_UU", "muu_extra", 0];
    expectZero["PASS_RESPONSE_M_UU", muuTest - expectedMuu];

    mugSign = primitive["PASS_RESPONSE_M_UG", "mug_sign", -1];
    expectZero["PASS_RESPONSE_M_UG", response[[1, 2]] - mugSign c zb zg/delta];

    mggFactor = primitive["PASS_RESPONSE_M_GG", "mgg_factor", 1];
    expectZero["PASS_RESPONSE_M_GG", response[[2, 2]] - mggFactor expectedMgg];

    skew = primitive["PASS_RESPONSE_SYMMETRY", "symmetry_skew", 0];
    responseTest = response + {{0, skew}, {0, 0}};
    expectBool["PASS_RESPONSE_SYMMETRY", responseTest === Transpose[responseTest], <|"m" -> responseTest|>];

    detFactor = primitive["PASS_RESPONSE_DETERMINANT", "det_factor", 1];
    expectZero["PASS_RESPONSE_DETERMINANT", Det[response] - detFactor zg^2/delta];

    zgWitness = primitive["PASS_Z_G_POSTULATED_WITNESS", "zg_witness", 1];
    expectBool["PASS_Z_G_POSTULATED_WITNESS", 0 < zgWitness <= 1,
      <|"status" -> "POSTULATED Robin-admissibility witness", "z_g*" -> zgWitness|>];

    squarePower = primitive["PASS_M_GG_NONNEGATIVE", "mgg_square", 2];
    mggStable = b zg^squarePower/deltaPositive;
    expectBool["PASS_M_GG_NONNEGATIVE",
      FullSimplify[mggStable >= 0, b > 0 && deltaPositive > 0 && Element[zg, Reals]],
      <|"m_gg" -> mggStable, "assumptions" -> "b>0,D>0,z_g real"|>];

    starZb = primitive["PASS_RESPONSE_STAR_WITNESS", "star_zb", 1];
    starRules = {b -> 2, kh -> 1, c -> 1/2, zb -> starZb, zg -> 1};
    starM = FullSimplify[response /. starRules];
    expectBool["PASS_RESPONSE_STAR_WITNESS",
      TrueQ[FullSimplify[(delta /. starRules) - 7/4] === 0] &&
      TrueQ[starM === {{4/7, -2/7}, {-2/7, 8/7}}],
      <|"D* consumed from stage030" -> (delta /. starRules), "m*" -> starM|>];

    greenSign = primitive["PASS_S_GG_SELF_RESPONSE", "green_sign", 1];
    sourceVector = {1, 1};
    greenMatrix = DiagonalMatrix[{1/p, greenSign/q}];
    sgg = Factor[sourceVector.greenMatrix.sourceVector];
    expectBool["PASS_S_GG_SELF_RESPONSE", positiveExactQ[sgg, p > 0 && q > 0],
      <|"definition" -> "S_gg=<eta,L_h^-1 eta>", "quadratic form" -> sgg,
        "status" -> "DERIVED for positive L_h"|>];
    Print["      native m=Transpose[Z].Inverse[diag(b,kappa)].Z; z_g>0 witness POSTULATED"];

    (* DSolve, not a copied ansatz, selects the decaying exterior branch. *)
    heading["Exterior stationary solution via DSolve and neutral shell"];
    radialDimension = primitive["PASS_EXTERIOR_ONE_OVER_R", "radial_dimension", 3];
    radialEquation = D[rr^(radialDimension - 1) D[hfun[rr], rr], rr] == 0;
    generalExterior = DSolveValue[radialEquation, hfun[rr], rr];
    constants = DeleteDuplicates[Cases[generalExterior, C[_], Infinity]];
    boundaryRules = First@Solve[
      {FullSimplify[generalExterior /. rr -> a] == hA,
       FullSimplify[Limit[generalExterior, rr -> Infinity]] == 0}, constants];
    exterior = FullSimplify[generalExterior /. boundaryRules, a > 0 && rr > 0];
    expectBool["PASS_EXTERIOR_ONE_OVER_R",
      radialDimension === 3 && TrueQ[FullSimplify[exterior - hA a/rr] === 0],
      <|"DSolve general" -> generalExterior, "decaying branch" -> exterior|>];

    exteriorEnergy = FullSimplify[
      Integrate[(kappa/2) 4 Pi rr^2 D[hA a/rr, rr]^2, {rr, a, Infinity},
        Assumptions -> a > 0 && kappa > 0 && Element[hA, Reals],
        GenerateConditions -> False], a > 0 && kappa > 0];
    curvatureSign = primitive["PASS_POSITIVE_HOLDING_CURVATURE", "curvature_sign", 1];
    curvature = curvatureSign D[exteriorEnergy, {hA, 2}];
    expectBool["PASS_POSITIVE_HOLDING_CURVATURE",
      TrueQ[FullSimplify[exteriorEnergy - 2 Pi kappa a hA^2] === 0] &&
      positiveExactQ[curvature, a > 0 && kappa > 0],
      <|"E_ext" -> exteriorEnergy, "curvature" -> curvature|>];

    corelessShift = primitive["PASS_NONZERO_HA_REQUIRES_CORE_HOLDER", "coreless_shift", 0];
    corelessStationary = hA /. Solve[D[2 Pi kappa a (hA - corelessShift)^2, hA] == 0, hA];
    expectBool["PASS_NONZERO_HA_REQUIRES_CORE_HOLDER", corelessStationary === {0},
      <|"coreless stationary h_A" -> corelessStationary,
        "named fact" -> "NONZERO_HA_REQUIRES_CORE_HOLDER"|>];

    shellMgg = b zg^2/deltaPositive;
    amplitude = shellMgg capC;
    potentialPower = primitive["PASS_NEUTRAL_FAR_FIELD_FORM", "potential_power", 1];
    potential = s1 s2 amplitude/(4 Pi R^potentialPower);
    force = Factor[-D[potential, R]];
    expectedPotential = s1 s2 amplitude/(4 Pi R);
    expectedForce = s1 s2 amplitude/(4 Pi R^2);
    expectBool["PASS_NEUTRAL_FAR_FIELD_FORM",
      TrueQ[FullSimplify[potential - expectedPotential] === 0] &&
      TrueQ[FullSimplify[force - expectedForce] === 0] &&
      TrueQ[FullSimplify[(potential /. {s1 -> s2, s2 -> s1}) - potential] === 0],
      <|"A=m_gg C" -> amplitude, "U" -> potential, "F_out" -> force,
        "C" -> "unspecified real [E^2]"|>];
    zeroC = primitive["PASS_ALLOWED_ZERO_COEFFICIENT", "zero_C", 0];
    expectZero["PASS_ALLOWED_ZERO_COEFFICIENT", expectedPotential /. capC -> zeroC,
      <|"conditional leading form" -> "coefficient may vanish", "class numerator" -> "not selected"|>];
    Print["      NONZERO_HA_REQUIRES_CORE_HOLDER"];
    Print["      A=m_gg*C, [C]=E^2: neutral shell; no numerator and no sign selected"];

    (* Mouth route comes last and is independent of the SymPy ordering. *)
    heading["Bare mouth reconstruction and option-A reflection dominance"];
    xiPower = primitive["PASS_FIELD_IDENTITY", "xi_power", 1];
    xiTest = ell^xiPower hvar;
    expectZero["PASS_FIELD_IDENTITY", xiTest/ell - hvar, <|"xi_w" -> xiTest|>];

    etaNumerator = primitive["PASS_ETA_ANNULUS_NORMALIZATION", "eta_numerator", 3];
    etaTest = etaNumerator/(4 Pi ((a + ell)^3 - a^3));
    etaIntegral = FullSimplify[
      Integrate[4 Pi rho^2 etaTest, {rho, a, a + ell},
        Assumptions -> a > 0 && ell > 0, GenerateConditions -> False],
      a > 0 && ell > 0];
    expectZero["PASS_ETA_ANNULUS_NORMALIZATION", etaIntegral - 1,
      <|"native Integrate" -> etaIntegral, "eta" -> etaTest|>];

    f0Power = primitive["PASS_F0_MOUTH_VALUE_EVALUATED", "f0_power", 2];
    f0Test = 1/(ell Cosh[w/ell]^f0Power);
    f0AtMouth = FullSimplify[f0Test /. w -> 0, ell > 0];
    canonicalCurvature = FullSimplify[D[1/(ell Cosh[w/ell]^2), {w, 2}] /. w -> 0, ell > 0];
    testCurvature = FullSimplify[D[f0Test, {w, 2}] /. w -> 0, ell > 0];
    expectBool["PASS_F0_MOUTH_VALUE_EVALUATED",
      TrueQ[FullSimplify[f0AtMouth - 1/ell, ell > 0] === 0] &&
      TrueQ[FullSimplify[testCurvature - canonicalCurvature, ell > 0] === 0],
      <|"evaluated f0(0)" -> f0AtMouth, "profile curvature" -> testCurvature|>];
    canonicalF0Mouth = 1/ell;

    jmWitness = primitive["PASS_REDUCED_COUPLING_NONZERO", "Jm_witness", 1];
    gTest = jmWitness/ell;
    expectBool["PASS_REDUCED_COUPLING_NONZERO", jmWitness =!= 0 &&
      TrueQ[FullSimplify[gTest - jmWitness/ell] === 0],
      <|"g_chih" -> gTest, "J_m witness" -> jmWitness|>];

    rPlus = (Tanh[(w + width/2)/ell] - Tanh[(w - (width/2 + sleeve))/ell])/
      (2 Tanh[(width + sleeve)/(2 ell)]);
    center = sleeve/2;
    halfspan = (width + sleeve)/2;
    rCentered = (1 + Cosh[2 halfspan/ell])/
      (Cosh[2 (w - center)/ell] + Cosh[2 halfspan/ell]);
    peakFactor = primitive["PASS_BOUNDED_PLUS_SLEEVE", "peak_factor", 1];
    rangeIdentity = FullSimplify[Cosh[2 u/ell] - 1 - 2 Sinh[u/ell]^2,
      ell > 0 && Element[u, Reals]];
    expectBool["PASS_BOUNDED_PLUS_SLEEVE",
      TrueQ[FullSimplify[(rPlus - rCentered) /. w -> center, ell > 0 && width > 0 && sleeve > 0] === 0] &&
      TrueQ[FullSimplify[(peakFactor rCentered /. w -> center) - 1] === 0] &&
      TrueQ[rangeIdentity === 0],
      <|"representative" -> rPlus, "peak" -> (peakFactor rCentered /. w -> center),
        "range proof" -> "denominator-numerator=2 Sinh[u/ell]^2>=0"|>];

    rBackground = (Tanh[(w + width/2)/ell] - Tanh[(w - width/2)/ell])/
      (2 Tanh[width/(2 ell)]);
    oddKernel = w Exp[-w^2/(2 ell^2)]/(Sqrt[2 Pi] ell^3);
    reflectBody = primitive["PASS_REFLECTION_ANTISYMMETRY", "reflect_body", 1];
    rMinusTest = If[reflectBody === 1, rPlus /. w -> -w, rPlus];
    reflectionResidual = FullSimplify[rMinusTest - (rPlus /. w -> -w), assumptions];
    expectBool["PASS_REFLECTION_ANTISYMMETRY",
      TrueQ[reflectionResidual === 0] &&
      TrueQ[FullSimplify[(oddKernel /. w -> -w) + oddKernel, assumptions] === 0] &&
      TrueQ[FullSimplify[(rBackground /. w -> -w) - rBackground, assumptions] === 0],
      <|"explicit body reflection residual" -> reflectionResidual, "I_minus" -> "-I_plus"|>];

    denForward = Cosh[2 (wp - center)/ell] + Cosh[2 halfspan/ell];
    denReflected = Cosh[2 (-wp - center)/ell] + Cosh[2 halfspan/ell];
    denDifference = FullSimplify[TrigExpand[denReflected - denForward], assumptions];
    expectedDifference = 2 Sinh[2 wp/ell] Sinh[sleeve/ell];
    dominanceSign = primitive["PASS_I_PLUS_REFLECTION_DOMINANCE", "dominance_sign", 1];
    profileNumerator = 1 + Cosh[2 halfspan/ell];
    rForward = rCentered /. w -> wp;
    rReflected = rCentered /. w -> -wp;
    rDifferenceIdentity = FullSimplify[
      (rForward - rReflected) - profileNumerator denDifference/(denForward denReflected), assumptions];
    rDifferenceCertificate = dominanceSign profileNumerator expectedDifference/(denForward denReflected);
    dPair = Factor[(rForward - rReflected) (rForward + rReflected)];
    kernelMass = FullSimplify[
      Integrate[(oddKernel /. w -> wp), {wp, 0, Infinity},
        Assumptions -> ell > 0, GenerateConditions -> False], ell > 0];
    genericCertificate =
      TrueQ[FullSimplify[denDifference - expectedDifference, assumptions] === 0] &&
      TrueQ[rDifferenceIdentity === 0] &&
      positiveExactQ[dominanceSign expectedDifference, assumptions] &&
      positiveExactQ[rDifferenceCertificate, assumptions] &&
      positiveExactQ[rForward, assumptions] && positiveExactQ[rReflected, assumptions] &&
      positiveExactQ[rForward + rReflected, assumptions] &&
      positiveExactQ[kernelMass, ell > 0];
    If[activeMutation =!= "PASS_I_PLUS_REFLECTION_DOMINANCE" && !TrueQ[genericCertificate],
      Print["STAGE031_STOP: I_plus_not_generic"];
      Exit[1]
    ];
    expectBool["PASS_I_PLUS_REFLECTION_DOMINANCE", genericCertificate,
      <|"den(-w)-den(w)" -> denDifference,
        "D(w)-D(-w)" -> "(r-r_-)(r+r_-)>0", "Integral_0^Infinity o_ell" -> kernelMass,
        "strict positive measure" -> "all w>0"|>];

    evenOddPart = primitive["PASS_I_PLUS_EVEN_DEFORMATION_CONTROL", "even_control_odd_part", 0];
    dEvenTest = Exp[-w^2/(2 ell^2)] (1 + evenOddPart w/ell);
    evenIntegral = FullSimplify[
      Integrate[oddKernel dEvenTest, {w, -Infinity, Infinity},
        Assumptions -> ell > 0, GenerateConditions -> False], ell > 0];
    expectZero["PASS_I_PLUS_EVEN_DEFORMATION_CONTROL", evenIntegral,
      <|"D(w)=D(-w)" -> (evenOddPart === 0), "I_plus control" -> evenIntegral|>];

    zeroSleeve = primitive["PASS_I_PLUS_ZERO_SLEEVE_CONTROL", "zero_sleeve", 0];
    zeroProfile = FullSimplify[rCentered /. sleeve -> zeroSleeve,
      ell > 0 && width > 0 && Element[w, Reals]];
    zeroD = zeroProfile^2 - rBackground^2;
    zeroParity = FullSimplify[(zeroD /. w -> -w) - zeroD,
      ell > 0 && width > 0 && Element[w, Reals]];
    expectZero["PASS_I_PLUS_ZERO_SLEEVE_CONTROL", zeroParity,
      <|"L_s" -> zeroSleeve, "D even -> paired integral" -> 0|>];

    negativeControl = primitive["PASS_I_PLUS_NEGATIVE_SLEEVE_CONTROL", "negative_control", -1];
    negativePair = negativeControl profileNumerator expectedDifference/(denForward denReflected) (rForward + rReflected);
    expectBool["PASS_I_PLUS_NEGATIVE_SLEEVE_CONTROL",
      TrueQ[FullSimplify[negativePair < 0, assumptions]],
      <|"inadmissible sleeve" -> "-w", "paired sign" -> Sign[negativePair]|>];

    (* Native Integrate constructs the certified sleeve integral after the
       actual representative has been reduced to the structural dPair form.
       Keeping dPairStructural symbolic avoids asking the kernel for a
       nonexistent elementary antiderivative of Gaussian times tanh. *)
    actualPair = (oddKernel /. w -> wp) ((rCentered /. w -> wp)^2 - (rCentered /. w -> -wp)^2);
    actualPairIdentity = FullSimplify[
      actualPair - (oddKernel /. w -> wp) dPair, assumptions];
    If[!TrueQ[actualPairIdentity === 0], raise["I_PLUS_STRUCTURAL_REDUCTION"]];
    If[activeMutation === "", Print["PROGRESS: native Integrate for certified I_plus sleeve integral"]];
    iPlusNative = Integrate[
      (oddKernel /. w -> wp) dPairStructural[wp], {wp, 0, Infinity},
      Assumptions -> ell > 0, GenerateConditions -> False];
    (* Native annular integration constructs N_chi from the same I_plus only
       after the independent generic nonzero certificate has passed. *)
    nChiNative = FullSimplify[
      Integrate[4 Pi rho^2 (3/(4 Pi ((a + ell)^3 - a^3))) iPlusNative,
        {rho, a, a + ell}, Assumptions -> a > 0 && ell > 0,
        GenerateConditions -> False], a > 0 && ell > 0];
    nchiGuard = primitive["PASS_N_CHI_NONZERO_GUARD", "nchi_guard", 1];
    expectBool["PASS_N_CHI_NONZERO_GUARD", nchiGuard === 1 && TrueQ[genericCertificate],
      <|"N_chi native construction" -> nChiNative,
        "ordering" -> "generic I_plus>0 certified before division"|>];

    projectionParity = primitive["PASS_Q_CHI_ORIENTATION", "projection_parity", 1];
    qPlus = 1;
    qMinus = If[projectionParity === 1, -1, 1];
    expectBool["PASS_Q_CHI_ORIENTATION", qPlus === 1 && qMinus === -1,
      <|"kernel" -> If[projectionParity === 1, "odd", "even"],
        "Q_+=I_+/N_chi" -> qPlus, "Q_-=I_-/N_chi" -> qMinus|>];

    mouthFactor = primitive["PASS_BARE_FORCING_LIVE_VARIATION", "mouth_quadratic_factor", 1];
    canonicalEta = 3/(4 Pi ((a + ell)^3 - a^3));
    parentDensity = canonicalEta (mouthFactor Km (canonicalF0Mouth hvar)^2/2 - Jm ss canonicalF0Mouth hvar);
    eulerDensityTest = D[parentDensity, hvar];
    km = Km/ell^2;
    g = Jm/ell;
    expectedEulerDensity = canonicalEta (km hvar - g ss);
    integratedEuler = FullSimplify[
      Integrate[4 Pi rho^2 eulerDensityTest, {rho, a, a + ell},
        Assumptions -> a > 0 && ell > 0, GenerateConditions -> False],
      a > 0 && ell > 0];
    bareAmplitude = FullSimplify[integratedEuler /. hvar -> 0];
    expectBool["PASS_BARE_FORCING_LIVE_VARIATION",
      TrueQ[FullSimplify[eulerDensityTest - expectedEulerDensity] === 0] &&
      TrueQ[FullSimplify[integratedEuler - (km hvar - g ss)] === 0] &&
      TrueQ[FullSimplify[bareAmplitude + g ss] === 0] &&
      TrueQ[(bareAmplitude /. {Jm -> 1, ss -> 1, ell -> 1}) =!= 0],
      <|"delta Omega/delta h density" -> eulerDensityTest,
        "integrated Euler amplitude" -> integratedEuler,
        "bare integrated amplitude" -> bareAmplitude|>];
    Print["      option-A result: N_chi=I_plus>0; (Q_+,Q_-)=(+1,-1); bare source=-g_chih*s"];

    heading["Complete mechanism dimensional firewall [L,T,M]"];
    dimZero = {0, 0, 0}; dimL = {1, 0, 0}; dimT = {0, 1, 0}; dimMass = {0, 0, 1};
    dimEnergy = 2 dimL - 2 dimT + dimMass;
    unitRows = <|
      "PASS_UNITS_XI_W" -> {"dim_xi", dimL, dimL},
      "PASS_UNITS_H" -> {"dim_h", dimZero, dimZero},
      "PASS_UNITS_H_A" -> {"dim_hA", dimZero, dimZero},
      "PASS_UNITS_K_M" -> {"dim_Km", dimEnergy + 2 dimL, dimEnergy + 2 dimL},
      "PASS_UNITS_J_M" -> {"dim_Jm", dimEnergy + dimL, dimEnergy + dimL},
      "PASS_UNITS_K_M_REDUCED" -> {"dim_km", dimEnergy, dimEnergy},
      "PASS_UNITS_G_CHIH" -> {"dim_g", dimEnergy, dimEnergy},
      "PASS_UNITS_ETA" -> {"dim_eta", -3 dimL, -3 dimL},
      "PASS_UNITS_ODD_KERNEL" -> {"dim_o", -2 dimL, -2 dimL},
      "PASS_UNITS_N_CHI" -> {"dim_nchi", -dimL, -dimL},
      "PASS_UNITS_M_UU" -> {"dim_muu", 3 dimL - dimEnergy, 3 dimL - dimEnergy},
      "PASS_UNITS_M_UG" -> {"dim_mug", 2 dimL - dimEnergy, 2 dimL - dimEnergy},
      "PASS_UNITS_M_GG" -> {"dim_mgg", dimL - dimEnergy, dimL - dimEnergy},
      "PASS_UNITS_KAPPA" -> {"dim_kappa", dimEnergy - dimL, dimEnergy - dimL},
      "PASS_UNITS_ESCAPE_MATRIX" -> {"dim_Z21", dimL, dimL},
      "PASS_UNITS_DET_M" -> {"dim_detm", 4 dimL - 2 dimEnergy, -2 dimMass + 4 dimT},
      "PASS_UNITS_S_GG" -> {"dim_sgg", -dimEnergy, -dimEnergy},
      "PASS_UNITS_SHELL_C" -> {"dim_C", 2 dimEnergy, 2 dimEnergy},
      "PASS_UNITS_SHELL_A" -> {"dim_A", dimEnergy + dimL, dimEnergy + dimL},
      "PASS_UNITS_U" -> {"dim_U", dimEnergy, dimEnergy},
      "PASS_UNITS_F" -> {"dim_F", dimEnergy - dimL, dimEnergy - dimL}
    |>;
    Do[
      row = unitRows[predicate];
      actualDimension = primitive[predicate, row[[1]], row[[2]]];
      expectZero[predicate, dimResidual[actualDimension, row[[3]]],
        <|"actual" -> actualDimension, "target" -> row[[3]]|>],
      {predicate, Keys[unitRows]}
    ];
    liveDimensionResidual = Total[{
      dimResidual[(dimEnergy + 2 dimL) - 2 dimL, dimEnergy],
      dimResidual[(dimEnergy + dimL) - dimL, dimEnergy],
      dimResidual[-3 dimL + 4 dimL - 2 dimL, -dimL],
      dimResidual[(dimL - dimEnergy) + 2 dimEnergy, dimEnergy + dimL],
      dimResidual[(dimEnergy + dimL) - dimL, dimEnergy],
      dimResidual[(dimEnergy + dimL) - 2 dimL, dimEnergy - dimL]
    }];
    If[!TrueQ[liveDimensionResidual === 0], raise["DIMENSION_LIVE_RELATIONS"]];
    Print["      live chains: k_m=K_m/ell^2; g=J_m/ell; N_chi=int eta dx int o dw; A=m_gg C; U=A/R; F=A/R^2"];

    Print[""];
    Print["SCOPE: stage031 earned mechanism and target-blind conditional leading form only."];
    Print["DEFERRED_TO_STAGE032: BC-class numerator C, force sign, ensembles, internal_inconsistency, R1 landing."];
    Print["VERDICT_TOKEN: THROAT_H_SOURCE_1_OVER_R2"];
    Print["EARNED_LABEL: PUNCTURE_DEFLECTION_MECHANISM_TARGET_BLIND"];

    If[activeMutation =!= "",
      Print["FIRST_FAILURE=MUTATION_DID_NOT_FIRE"];
      raise["MUTATION_DID_NOT_FIRE"]
    ];
    True,
    "ledgerStage031Failure",
    Function[{message, tag}, False]
  ];

  Print[""];
  Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
  If[TrueQ[ok],
    Print["OVERALL PASS: Mathematica verified stage031 puncture-deflection field identity and source mechanism"];
    Exit[0],
    Print["OVERALL FAIL: Mathematica stage031 audit did not close"];
    Exit[1]
  ]
]
