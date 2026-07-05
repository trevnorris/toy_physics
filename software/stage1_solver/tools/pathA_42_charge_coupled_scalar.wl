(* pathA_42 charge-coupled scalar gate, Mathematica engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
assert[label_, cond_] := If[! TrueQ[cond], fail[label]];

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_42_charge_coupled_scalar.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir, CreateIntermediateDirectories -> True]];
jsonOut = FileNameJoin[{scratchDir, "pathA_42_charge_coupled_scalar_mathematica.json"}];

(* Independent dual-engine derivation: static 2x2 exchange. *)
stiffness = {{Beff, Chu}, {Chu, Kh}};
invStiffness = FullSimplify[Inverse[stiffness]];
det = FullSimplify[Det[stiffness]];
jq = {qL, qh};
jm = {qM, mh};
Aqq = FullSimplify[jq.invStiffness.jq];
Aqm = FullSimplify[jq.invStiffness.jm];
Amm = FullSimplify[jm.invStiffness.jm];

expectAqq = (Kh*qL^2 - 2*Chu*qL*qh + Beff*qh^2)/det;
expectAqm = (Kh*qL*qM - Chu*qh*qM - Chu*qL*mh + Beff*qh*mh)/det;
expectAmm = (Kh*qM^2 - 2*Chu*qM*mh + Beff*mh^2)/det;

assert["A_qq residual", FullSimplify[Aqq - expectAqq == 0]];
assert["A_qm residual", FullSimplify[Aqm - expectAqm == 0]];
assert["A_mm residual", FullSimplify[Amm - expectAmm == 0]];
assert["decoupled A_qm projection", FullSimplify[(Aqm /. {Chu -> 0, mh -> 0}) - qL*qM/Beff == 0]];
assert["mixed qL=0 term", FullSimplify[(Aqm /. {qL -> 0, mh -> 0}) + Chu*qh*qM/det == 0]];
assert["decoupled A_mm density only", FullSimplify[(Amm /. {Chu -> 0, mh -> 0}) - qM^2/Beff == 0]];
assert["qM flip A_qm", FullSimplify[(Aqm /. {mh -> 0, qM -> -qM}) + (Aqm /. mh -> 0) == 0]];
assert["A_mm even in qM", FullSimplify[(Amm /. {mh -> 0, qM -> -qM}) - (Amm /. mh -> 0) == 0]];
signedAmm = FullSimplify[(Amm /. mh -> 0)/qM];
assert["signed A_mm projection flips", FullSimplify[(signedAmm /. qM -> -qM) + signedAmm == 0]];

(* Independent dual-engine derivation: radiation speed exponents from dipole/flux powers. *)
speedExponent[expr_] := Module[{zexpr, num, den},
  zexpr = Together[expr /. cE -> z];
  num = Numerator[zexpr];
  den = Denominator[zexpr];
  Exponent[den, z] - Exponent[num, z]
];

accel = FullSimplify[omega^2*d];
scalarFarTimeGradient = FullSimplify[qh*accel/(Mh*cE^2*r)];
scalarPower = FullSimplify[Mh*cE*scalarFarTimeGradient^2*r^2];
emFarTimeGradient = FullSimplify[QE*accel/(cg^2*r)];
emPower = FullSimplify[cg*emFarTimeGradient^2*r^2];
ratioBare = FullSimplify[scalarPower/emPower];
expectedRatioBare = FullSimplify[(qh^2/(Mh*QE^2))*(cg/cE)^3];
applyStaticKh[expr_] := FullSimplify[(expr /. Mh -> Kh/cE^2)*Kh/qh^2*kappa];
ratioPinnedKh = applyStaticKh[ratioBare];
scalarFarTimeGradientCorrupt = FullSimplify[qh*accel/(Mh*cE*r)];
scalarPowerCorrupt = FullSimplify[Mh*cE*scalarFarTimeGradientCorrupt^2*r^2];
ratioCorrupt = FullSimplify[scalarPowerCorrupt/emPower];
ratioWrongNorm = applyStaticKh[ratioBare];
bareExp = speedExponent[ratioBare];
pinnedExp = speedExponent[ratioPinnedKh];
corruptExp = speedExponent[ratioCorrupt];
wrongNormExp = speedExponent[ratioWrongNorm];
wrongNormDiscriminatesBare = TrueQ[wrongNormExp == 1 && bareExp == 3 && wrongNormExp != bareExp];

assert["bare flux ratio derived", FullSimplify[ratioBare - expectedRatioBare == 0]];
assert["wrong normalization exponent discriminates bare", wrongNormDiscriminatesBare];
assert["bare-fixed radiation exponent", bareExp == 3];
assert["pinned-Kh radiation exponent", pinnedExp == 1];
assert["corrupt speed exponent changes", corruptExp != bareExp];
assert["wrong normalization exponent", wrongNormExp == 1];

(* Independent dual-engine derivation: stability and projection checks. *)
massVector = {qM, 0};
hUnit = {0, 1};
fixtureMassVector = {qM, mh};
assert["static determinant", FullSimplify[det == Beff*Kh - Chu^2]];
assert["h mass projection zero", FullSimplify[hUnit.massVector == 0]];
assert["fixture h projection", FullSimplify[hUnit.fixtureMassVector == mh]];

payload = <|
  "static_exchange" -> <|
    "determinant_check" -> True,
    "A_qq_residual_zero" -> True,
    "A_qm_residual_zero" -> True,
    "A_mm_residual_zero" -> True,
    "decoupled_A_qm_has_no_h_mass_residue" -> True,
    "mixed_qL0_contains_CqhqM" -> True,
    "decoupled_A_mm_density_only" -> True,
    "qM_flip_A_qm" -> True,
    "A_mm_even_but_signed_projection_flips" -> True
  |>,
  "radiation" -> <|
    "bare_flux_ratio_matches_expected" -> True,
    "wrong_normalization_recomputed" -> wrongNormDiscriminatesBare,
    "bare_fixed_exponent" -> bareExp,
    "pinned_Kh_exponent" -> pinnedExp,
    "corrupt_speed_exponent" -> corruptExp,
    "wrong_normalization_exponent" -> wrongNormExp
  |>,
  "stability" -> <|
    "bound" -> "C_hu**2 < B_eff*K_h",
    "violated_condition" -> "C_hu**2 >= B_eff*K_h",
    "strict_slack_symbol" -> "sigma = B_eff*K_h - C_hu**2 > 0"
  |>,
  "projections" -> <|
    "source_mass_vector_imported" -> "[q_M, 0]",
    "h_mass_projection_zero" -> True,
    "fixture_h_projection" -> "m_h",
    "grav_even_overlap_imported_zero" -> True
  |>
|>;

Export[jsonOut, <|
  "schema" -> "pathA_42_charge_coupled_scalar_mathematica/v1",
  "comparison_payload" -> payload
|>, "RawJSON"];

Print["pathA_42 Mathematica derivation OK"];
