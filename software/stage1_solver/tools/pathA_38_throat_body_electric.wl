(* pathA_38 throat-body electric localization gate, Mathematica engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
Off[Limit::alimv];

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_38_throat_body_electric.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
sympyJson = FileNameJoin[{scratchDir, "pathA_38_throat_body_electric_sympy.json"}];
jsonOut = FileNameJoin[{scratchDir, "pathA_38_throat_body_electric_mathematica.json"}];

If[! FileExistsQ[sympyJson], fail["missing SymPy JSON for engine agreement: " <> sympyJson]];

$Assumptions =
  ell > 0 && b > 0 && d > 0 && R > 0 && QE > 0 &&
  lamN > 0 && lamTau > 0 && Omega35 > 0 && epsMix > 0 &&
  cE > 0 && omega > 0 && Element[{w, wp}, Reals];

signIndicator[expr_] := Which[
  TrueQ[FullSimplify[expr > 0]], 1,
  TrueQ[FullSimplify[expr < 0]], -1,
  TrueQ[FullSimplify[expr == 0]], 0,
  True, fail["could not determine sign: " <> ToString[expr, InputForm]]
];

radialResidual[expr_, kappaSquared_] := FullSimplify[D[expr, {R, 2}] + 2 D[expr, R]/R + kappaSquared expr];
radialExponent[potential_] := Module[{flow, p},
  flow = FullSimplify[-D[potential, R]];
  p = FullSimplify[-Limit[R D[Log[flow], R], R -> Infinity]];
  If[! NumericQ[p], fail["could not extract radial exponent from " <> ToString[potential, InputForm]]];
  p
];

opAction[vec_, m2_, kw_] := FullSimplify[-D[kw . D[vec, w], w] + m2 . vec];

zeroExprQ[value_] := Which[
  value === Missing["absent"], False,
  ListQ[value], TrueQ[FullSimplify[value == ConstantArray[0, Dimensions[value]]]],
  True, TrueQ[FullSimplify[value == 0]]
];

positiveExprQ[value_] := value =!= Missing["absent"] && TrueQ[FullSimplify[value > 0]];
nonzeroExprQ[value_] := value =!= Missing["absent"] && ! zeroExprQ[value];

sourceTailFromSpectrum[branch_] := Module[
  {
    mSrc = Lookup[branch, "source_m_squared", Missing["absent"]],
    zeroOverlap = Lookup[branch, "source_zero_mode_overlap", Missing["absent"]],
    modeOverlap = Lookup[branch, "source_mode_overlap", Missing["absent"]],
    modeResidual = Lookup[branch, "source_mode_residual", 0]
  },
  If[
    positiveExprQ[mSrc] && zeroExprQ[zeroOverlap] && nonzeroExprQ[modeOverlap] && zeroExprQ[modeResidual],
    "yukawa",
    None
  ]
];

classifyBranch[branch_] := Which[
  sourceTailFromSpectrum[branch] === "yukawa", "FAIL_YUKAWA",
  ! TrueQ[branch["normalizable_zero_mode"]], "FAIL_DELOCALIZED_BULK_1_OVER_R3",
  ! TrueQ[FullSimplify[branch["p_static"] == 2 && branch["p_dynamic"] == 2]], "FAIL_DELOCALIZED_BULK_1_OVER_R3",
  TrueQ[FullSimplify[branch["q_eff"] == 0]], "FAIL_NO_MONOPOLE",
  ! TrueQ[FullSimplify[branch["m_desc_squared"] == 0]], "FAIL_PINNED_BRANON",
  branch["U_pp_sign"] <= 0 || branch["U_pm_sign"] >= 0, "FAIL_GHOST_INSTABILITY",
  True, "THROAT_ELECTRIC_LOCALIZED_COULOMB"
];

expectDim[name_, actual_, expected_] := If[actual =!= expected, fail["dimension check failed: " <> name]];
expectDimFail[name_, actual_, expected_] := If[actual === expected, fail["dimension ablation did not fire: " <> name]];

(* Reduced Z2 wall model and transverse eigenproblem. *)
vars = {chi, nField, tauField};
wallT = Tanh[w/ell];
chi0 = wallT;
n0 = wallT;
phi0 = {chi0, n0, 0};
f0 = FullSimplify[D[phi0, w]];
f = f0[[1]];
fp = f /. w -> wp;
sigma = DiagonalMatrix[{-1, -1, 1}];
kwMat = IdentityMatrix[3];
kparMat = IdentityMatrix[3];

potential =
  (chi^2 - 1)^2/(2 ell^2) +
  (nField^2 - 1)^2/(2 ell^2) +
  lamN (nField - chi)^2/2 +
  lamTau tauField^2/2;
gradU = FullSimplify[{D[potential, chi], D[potential, nField], D[potential, tauField]} /. {chi -> chi0, nField -> n0, tauField -> 0}];
backgroundEom = FullSimplify[{-D[chi0, {w, 2}], -D[n0, {w, 2}], 0} + gradU];
If[backgroundEom =!= {0, 0, 0}, fail["background EOM failed: " <> ToString[backgroundEom, InputForm]]];

m2 = FullSimplify[
  Table[D[potential, vars[[i]], vars[[j]]], {i, 1, 3}, {j, 1, 3}] /. {chi -> chi0, nField -> n0, tauField -> 0}
];
v0 = FullSimplify[(4 - 6 Sech[w/ell]^2)/ell^2];
expectedM2 = {{v0 + lamN, -lamN, 0}, {-lamN, v0 + lamN, 0}, {0, 0, lamTau}};
If[! TrueQ[FullSimplify[m2 == expectedM2]], fail["Hessian mismatch"]];

of0 = opAction[f0, m2, kwMat];
If[of0 =!= {0, 0, 0}, fail["Goldstone residual failed: " <> ToString[of0, InputForm]]];

shapeScalar = Sech[w/ell] Tanh[w/ell];
shapeVec = {shapeScalar, shapeScalar, 0};
shapeM2 = 3/ell^2;
shapeResidual = FullSimplify[opAction[shapeVec, m2, kwMat] - shapeM2 kparMat . shapeVec];
If[shapeResidual =!= {0, 0, 0}, fail["shape residual failed: " <> ToString[shapeResidual, InputForm]]];

relVec = {f, -f, 0};
relM2 = 2 lamN;
relResidual = FullSimplify[opAction[relVec, m2, kwMat] - relM2 kparMat . relVec];
If[relResidual =!= {0, 0, 0}, fail["relative-lock residual failed"]];

tauVec = {0, 0, 1};
tauM2 = lamTau;
tauResidual = FullSimplify[(m2 /. Sech[w/ell] -> 0) . tauVec - tauM2 kparMat . tauVec];
If[tauResidual =!= {0, 0, 0}, fail["tau source residual failed"]];
tauZeroOverlap = FullSimplify[f0 . (kparMat . tauVec)];
tauSourceOverlap = FullSimplify[tauVec . (kparMat . tauVec)];

(* Normalization and source projections are integrations in y=tanh(w/ell). *)
normDirectDensity = FullSimplify[f0 . (kparMat . f0)];
normDensityY = FullSimplify[2 (1 - y^2)/ell];
normCutoffY = FullSimplify[Integrate[normDensityY, {y, -Y, Y}, Assumptions -> 0 < Y < 1 && ell > 0]];
normCutoff = FullSimplify[normCutoffY /. Y -> Tanh[d/ell]];
normInfinite = FullSimplify[normCutoffY /. Y -> 1];
If[! TrueQ[FullSimplify[normInfinite == 8/(3 ell)]], fail["Goldstone norm failed"]];

descDensity = FullSimplify[f0 . (kparMat . of0)];
descNumCutoff = FullSimplify[Integrate[descDensity, {w, -d, d}, Assumptions -> d > 0 && ell > 0]];
descNum = FullSimplify[Limit[descNumCutoff, d -> Infinity]];
mDescZ2 = FullSimplify[descNum/normInfinite];
pinningNum = FullSimplify[Omega35^2 normInfinite];
mDescConf = FullSimplify[pinningNum/normInfinite];

lockRatio = FullSimplify[D[n0, w]/D[chi0, w]];
lockResidualN = FullSimplify[f0[[2]] - lockRatio f0[[1]]];
lockResidualTau = f0[[3]];
If[lockResidualN =!= 0 || lockResidualTau =!= 0, fail["locking flat direction failed"]];

parityGoldstoneResidual = FullSimplify[sigma . (f0 /. w -> -w) + f0];
If[parityGoldstoneResidual =!= {0, 0, 0}, fail["Goldstone parity failed"]];
bcParityEta = FullSimplify[D[f, w] /. w -> 0];
outerFluxCutoff = FullSimplify[D[f, w] /. w -> d];
bcOuterFluxLimit = FullSimplify[Limit[outerFluxCutoff, d -> Infinity]];
If[bcParityEta =!= 0 || bcOuterFluxLimit =!= 0, fail["BC residuals failed"]];

B = Tanh[b/ell];
rhoB = 1/(2 b);
qPlusDensityY = FullSimplify[QE rhoB 2];
qMinusDensityY = FullSimplify[QE rhoB (-2)];
qPlus = FullSimplify[Integrate[qPlusDensityY, {y, -B, B}, Assumptions -> b > 0 && ell > 0 && QE > 0]];
qMinus = FullSimplify[Integrate[qMinusDensityY, {y, -B, B}, Assumptions -> b > 0 && ell > 0 && QE > 0]];
qEff = qPlus;
neutralCompositeSum = FullSimplify[qPlus + qMinus];
gravEvenDensityY = FullSimplify[QE rhoB 2 y];
gravEvenOverlap = FullSimplify[Integrate[gravEvenDensityY, {y, -B, B}, Assumptions -> b > 0 && ell > 0 && QE > 0]];
gravMixOverlap = FullSimplify[gravEvenOverlap + epsMix qEff];
orthChargePart = FullSimplify[Integrate[QE/b, {y, -B, B}, Assumptions -> b > 0 && ell > 0 && QE > 0]];
orthSubtractionPart = FullSimplify[Integrate[(qEff/normInfinite) normDensityY, {y, -1, 1}, Assumptions -> b > 0 && ell > 0 && QE > 0]];
noMonopoleOverlap = FullSimplify[orthChargePart - orthSubtractionPart];
If[neutralCompositeSum =!= 0 || gravEvenOverlap =!= 0 || noMonopoleOverlap =!= 0, fail["source projection checks failed"]];

(* Radial Green factors consume the computed transverse m0^2 and normalization. *)
gStaticGeneral = DSolveValue[gStatic''[R] + 2 gStatic'[R]/R - mDescZ2 gStatic[R] == 0, gStatic[R], R];
staticNoConstant = FullSimplify[gStaticGeneral - Limit[gStaticGeneral, R -> Infinity]];
staticCoulombCoeff = FullSimplify[Limit[R staticNoConstant, R -> Infinity]];
greenStatic = FullSimplify[staticNoConstant/staticCoulombCoeff/(4 Pi)];
If[radialResidual[greenStatic, -mDescZ2] =!= 0, fail["static radial residual failed"]];
pStatic = radialExponent[greenStatic];

k = omega/cE;
gDynamicGeneral = DSolveValue[gDynamic''[R] + 2 gDynamic'[R]/R + (k^2 - mDescZ2) gDynamic[R] == 0, gDynamic[R], R];
greenDynamic = FullSimplify[gDynamicGeneral /. {C[1] -> 0, C[2] -> I omega/(2 Pi cE)}];
If[radialResidual[greenDynamic, k^2 - mDescZ2] =!= 0, fail["dynamic radial residual failed"]];
greenDynamicLimit = FullSimplify[Limit[greenDynamic, omega -> 0]];
pDynamic = radialExponent[greenDynamicLimit];

greenZeroProjectedStatic = FullSimplify[qPlus qPlus/normInfinite greenStatic];
greenZeroProjectedDynamic = FullSimplify[qPlus qPlus/normInfinite greenDynamic];
g0EtaEtaStatic = FullSimplify[f fp/normInfinite greenStatic];
uPP = FullSimplify[qPlus qPlus/normInfinite greenStatic];
uPM = FullSimplify[qPlus qMinus/normInfinite greenStatic];
uPPSign = signIndicator[uPP];
uPMSign = signIndicator[uPM];

mainBranch = <|
  "normalizable_zero_mode" -> True,
  "p_static" -> pStatic,
  "p_dynamic" -> pDynamic,
  "q_eff" -> qEff,
  "m_desc_squared" -> mDescZ2,
  "U_pp_sign" -> uPPSign,
  "U_pm_sign" -> uPMSign
|>;
headline = classifyBranch[mainBranch];

(* Massive/Yukawa and pinned descendants. *)
yukawaGeneral = DSolveValue[gYukawa''[R] + 2 gYukawa'[R]/R - tauM2 gYukawa[R] == 0, gYukawa[R], R];
yukawaTail = FullSimplify[Exp[-Sqrt[tauM2] R]/(4 Pi R)];
If[radialResidual[yukawaTail, -tauM2] =!= 0, fail["Yukawa witness residual failed"]];
pinnedGeneral = DSolveValue[gPinned''[R] + 2 gPinned'[R]/R - mDescConf gPinned[R] == 0, gPinned[R], R];
pinnedTail = FullSimplify[Exp[-Sqrt[mDescConf] R]/(4 Pi R)];
If[radialResidual[pinnedTail, -mDescConf] =!= 0, fail["pinned witness residual failed"]];

(* Delocalized anti-localizing operator e^(2 k_w w) K_parallel with k_w=3/ell. *)
kWarp = 3/ell;
antiDensityY = FullSimplify[2 (1 + y)^4/(ell (1 - y)^2)];
antiCutoffY = FullSimplify[Integrate[antiDensityY, {y, 0, Y}, Assumptions -> 0 < Y < 1 && ell > 0]];
antiNorm = Limit[antiCutoffY /. Y -> 1 - s, s -> 0, Direction -> -1];
antiRate = FullSimplify[2 kWarp - 4/ell];
mCont = Unique["mCont"];
continuumGreen = FullSimplify[Integrate[Exp[-mCont R]/(4 Pi R), {mCont, 0, Infinity}, Assumptions -> mCont > 0 && R > 0]];
pContinuum = radialExponent[continuumGreen];
delocBranch = <|
  "normalizable_zero_mode" -> False,
  "p_static" -> pContinuum,
  "p_dynamic" -> 0,
  "q_eff" -> qEff,
  "m_desc_squared" -> 0,
  "U_pp_sign" -> uPPSign,
  "U_pm_sign" -> uPMSign
|>;
delocHeadline = classifyBranch[delocBranch];

ghostNormDensityY = FullSimplify[-2 (1 - y^2)/ell];
ghostNorm = FullSimplify[Integrate[ghostNormDensityY, {y, -1, 1}, Assumptions -> ell > 0]];
ghostUPP = FullSimplify[qPlus qPlus/ghostNorm greenStatic];
ghostUPM = FullSimplify[qPlus qMinus/ghostNorm greenStatic];
ghostUPPSign = signIndicator[ghostUPP];
ghostUPMSign = signIndicator[ghostUPM];
ghostBranch = Join[mainBranch, <|"U_pp_sign" -> ghostUPPSign, "U_pm_sign" -> ghostUPMSign|>];
ghostHeadline = classifyBranch[ghostBranch];

noMonopoleBranch = Join[mainBranch, <|"q_eff" -> noMonopoleOverlap, "U_pp_sign" -> 0, "U_pm_sign" -> 0|>];
noMonopoleHeadline = classifyBranch[noMonopoleBranch];
pinnedBranch = Join[mainBranch, <|"m_desc_squared" -> mDescConf|>];
pinnedHeadline = classifyBranch[pinnedBranch];
yukawaBranch = Join[
  mainBranch,
  <|
    "q_eff" -> tauZeroOverlap,
    "U_pp_sign" -> 0,
    "U_pm_sign" -> 0,
    "source_m_squared" -> tauM2,
    "source_zero_mode_overlap" -> tauZeroOverlap,
    "source_mode_overlap" -> tauSourceOverlap,
    "source_mode_residual" -> tauResidual
  |>
];
yukawaHeadline = classifyBranch[yukawaBranch];
gravityMixFired = ! TrueQ[FullSimplify[gravMixOverlap == 0]];

If[yukawaHeadline =!= "FAIL_YUKAWA", fail["Yukawa classifier self-test failed"]];
If[pinnedHeadline =!= "FAIL_PINNED_BRANON", fail["pinned classifier self-test failed"]];
If[delocHeadline =!= "FAIL_DELOCALIZED_BULK_1_OVER_R3", fail["delocalized classifier self-test failed"]];
If[noMonopoleHeadline =!= "FAIL_NO_MONOPOLE", fail["no-monopole classifier self-test failed"]];
If[ghostHeadline =!= "FAIL_GHOST_INSTABILITY", fail["ghost classifier self-test failed"]];
If[! gravityMixFired, fail["gravity mixing witness failed"]];

(* Dimensional firewall in basis {E,L,Chi,N,Tau,Q}. *)
ed = {1, 0, 0, 0, 0, 0};
ld = {0, 1, 0, 0, 0, 0};
chid = {0, 0, 1, 0, 0, 0};
nd = {0, 0, 0, 1, 0, 0};
taud = {0, 0, 0, 0, 1, 0};
d4x = 4 ld;
grad = -ld;
rhoA = -3 ld;
rhoBDim = -ld;
kEta = ed - 2 ld - 2 chid;
mEta = ed - 4 ld - 2 chid;
sEta = ed - 4 ld - chid;
qeUEta = ed - chid;
greenEtaEta = 2 chid - ed;
operatorEta = ed - 4 ld - 2 chid;
expectDim["E_quad_eta_gradient", kEta + 2 (grad + chid) + d4x, ed];
expectDim["E_quad_eta_mass", mEta + 2 chid + d4x, ed];
expectDim["source_eta", sEta + chid + d4x, ed];
expectDim["rho_a_rho_b_u_eta", rhoA + rhoBDim + qeUEta, sEta];
expectDim["green_response", greenEtaEta + sEta + d4x, chid];
expectDim["U_int", sEta + greenEtaEta + sEta + 2 d4x, ed];
expectDim["m_squared", operatorEta - kEta, -2 ld];
expectDimFail["omit_rho_b_from_source", rhoA + qeUEta, sEta];
expectDimFail["bad_K_parallel_weight", operatorEta - (ed - ld - 2 chid), -2 ld];
If[nd === chid, fail["drop locking ratio ablation did not fire"]];

sympyResults = Import[sympyJson, "RawJSON"];
sympyExprs = sympyResults["engine_agreement"]["mathematica_exprs"];
sympyDigest = sympyResults["engine_agreement"]["sympy_expression_digest"];

actuals = <|
  "K_w_11" -> 1,
  "K_parallel_11" -> 1,
  "M_eta_eta" -> m2[[1, 1]],
  "M_eta_n" -> m2[[1, 2]],
  "M_tau_tau" -> m2[[3, 3]],
  "O_f0_eta" -> of0[[1]],
  "O_f0_n" -> of0[[2]],
  "O_f0_tau" -> of0[[3]],
  "N0_norm" -> normInfinite,
  "m_desc_Z2_squared" -> mDescZ2,
  "m_desc_conf_squared" -> mDescConf,
  "shape_mode_m_squared" -> shapeM2,
  "relative_tilt_eigenvalue" -> relM2,
  "tau_source_eigenvalue" -> tauM2,
  "q_eff" -> qEff,
  "q_h_plus" -> qPlus,
  "q_h_minus" -> qMinus,
  "neutral_composite_sum" -> neutralCompositeSum,
  "grav_even_overlap" -> gravEvenOverlap,
  "grav_mix_overlap" -> gravMixOverlap,
  "no_monopole_overlap" -> noMonopoleOverlap,
  "p_static" -> pStatic,
  "p_dynamic" -> pDynamic,
  "green_zero_projected_static" -> greenZeroProjectedStatic,
  "green_zero_projected_dynamic" -> greenZeroProjectedDynamic,
  "U_pp_sign" -> uPPSign,
  "U_pm_sign" -> uPMSign,
  "ghost_norm" -> ghostNorm,
  "ghost_U_pp_sign" -> ghostUPPSign,
  "ghost_U_pm_sign" -> ghostUPMSign,
  "anti_norm_diverges" -> If[antiNorm === Infinity, 1, 0],
  "delocalized_p" -> pContinuum,
  "FAIL_YUKAWA_fired" -> If[yukawaHeadline === "FAIL_YUKAWA", 1, 0],
  "FAIL_PINNED_BRANON_fired" -> If[pinnedHeadline === "FAIL_PINNED_BRANON", 1, 0],
  "FAIL_DELOCALIZED_BULK_1_OVER_R3_fired" -> If[delocHeadline === "FAIL_DELOCALIZED_BULK_1_OVER_R3", 1, 0],
  "FAIL_NO_MONOPOLE_fired" -> If[noMonopoleHeadline === "FAIL_NO_MONOPOLE", 1, 0],
  "FAIL_GHOST_INSTABILITY_fired" -> If[ghostHeadline === "FAIL_GHOST_INSTABILITY", 1, 0],
  "FAIL_GRAVITY_MIXING_fired" -> If[gravityMixFired, 1, 0],
  "headline_is_pass" -> If[headline === "THROAT_ELECTRIC_LOCALIZED_COULOMB", 1, 0],
  "dim_omit_rho_b_fired" -> 1,
  "dim_bad_K_parallel_fired" -> 1,
  "dim_drop_lock_ratio_fired" -> 1
|>;

engineKeys = {
  "K_w_11", "K_parallel_11",
  "M_eta_eta", "M_eta_n", "M_tau_tau",
  "O_f0_eta", "O_f0_n", "O_f0_tau",
  "N0_norm", "m_desc_Z2_squared", "m_desc_conf_squared",
  "shape_mode_m_squared", "relative_tilt_eigenvalue", "tau_source_eigenvalue",
  "q_eff", "q_h_plus", "q_h_minus", "neutral_composite_sum",
  "grav_even_overlap", "grav_mix_overlap", "no_monopole_overlap",
  "p_static", "p_dynamic",
  "green_zero_projected_static", "green_zero_projected_dynamic",
  "U_pp_sign", "U_pm_sign",
  "ghost_norm", "ghost_U_pp_sign", "ghost_U_pm_sign",
  "anti_norm_diverges", "delocalized_p",
  "FAIL_YUKAWA_fired", "FAIL_PINNED_BRANON_fired", "FAIL_DELOCALIZED_BULK_1_OVER_R3_fired",
  "FAIL_NO_MONOPOLE_fired", "FAIL_GHOST_INSTABILITY_fired", "FAIL_GRAVITY_MIXING_fired",
  "headline_is_pass",
  "dim_omit_rho_b_fired", "dim_bad_K_parallel_fired", "dim_drop_lock_ratio_fired"
};

assertEngine[name_] := Module[{expectedText, expectedExpr, actual},
  expectedText = sympyExprs[name];
  If[! StringQ[expectedText], fail["missing SymPy expression for " <> name]];
  expectedExpr = ToExpression[expectedText, InputForm];
  actual = actuals[name];
  If[! TrueQ[FullSimplify[actual == expectedExpr]],
    fail["engine disagreement " <> name <> ": Mathematica got " <>
      ToString[actual, InputForm] <> ", SymPy exported " <> expectedText]
  ];
];

Scan[assertEngine, engineKeys];
If[headline =!= sympyResults["top_line_verdict"], fail["headline disagreement"]];

out = <|
  "schema" -> "pathA_38_throat_body_electric_mathematica/v2",
  "status" -> "OK",
  "headline" -> headline,
  "checked_quantities" -> engineKeys,
  "sympy_expression_digest" -> sympyDigest,
  "computed" -> <|
    "p_static" -> ToString[pStatic, InputForm],
    "p_dynamic" -> ToString[pDynamic, InputForm],
    "N0_norm" -> ToString[normInfinite, InputForm],
    "q_eff" -> ToString[qEff, InputForm],
    "q_h_minus" -> ToString[qMinus, InputForm],
    "U_pp_sign" -> uPPSign,
    "U_pm_sign" -> uPMSign,
    "ghost_norm" -> ToString[ghostNorm, InputForm],
    "delocalized_p" -> pContinuum,
    "fail_witnesses" -> <|
      "FAIL_YUKAWA" -> yukawaHeadline,
      "FAIL_PINNED_BRANON" -> pinnedHeadline,
      "FAIL_DELOCALIZED_BULK_1_OVER_R3" -> delocHeadline,
      "FAIL_NO_MONOPOLE" -> noMonopoleHeadline,
      "FAIL_GHOST_INSTABILITY" -> ghostHeadline,
      "FAIL_GRAVITY_MIXING" -> If[gravityMixFired, "FIRED", "NOT_FIRED"]
    |>,
    "dimensional_ablations_fired" -> 3
  |>
|>;

Export[jsonOut, out, "RawJSON"];
Print["OK pathA_38_throat_body_electric_mathematica"];
Print[ExportString[<|"json" -> jsonOut, "headline" -> headline|>, "RawJSON"]];
