# Measurements backing `S11c_b_adjudication_build_directive.md`
Regenerated from the commands below (rule 2). NOTE: the reconcile run was regenerated against the
time-order-fixed comparator; the family tally here is post-fix (the directive is re-folded to v4 against
these corrected artifacts). The Bridge A identity and protected-set citations are source-stable.

## Claim (Bridge A): B_rho_3 = bRho * W_0, from the theta^2 energy coefficient in each engine
$ grep -n "B_rho_3" research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
133:B_rho_3 = inherited_symbol("B_rho_3", "KNOB", "uniform integrated compression modulus", DIM_ENERGY)
1135:        (btheta**2, B_rho_3 * W_bg / (2 * W0)),
$ sed -n "471,479p" research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl   # WL: uniformCoefficients[[3]] = bRho anchoredWidth, uniformFactors[[3]] = 1/2
  uniformCoefficients = {
    anchoredModulus, modulusDivU, bRho anchoredWidth,
    kW anchoredWidth^2, cCross anchoredWidth,
    gradientThetaCoefficient,
    kappaW anchoredWidth^4,
    gradientThetaEwCoefficient,
    thetaDivUCoefficient,
    ewDivUCoefficient anchoredWidth/WZero};
  uniformFactors = {1/2, 1/2, 1/2, 1/2, 1, 1/2, 1/2, 1, 1, 1};
# alignment: WL (1/2)*bRho*W_bg  ==  PY B_rho_3*W_bg/(2*W0)  <=>  B_rho_3 = bRho*W_0

## Claim: bRho is the load-bearing residue Bridge A touches (families with bRho in the residual)
$ grep "^A_minus_B " /home/trevnorris/.s11_build/S11c_b_reconcile_run.out | grep bRho | grep -oE "A_minus_B [A-Z_]+ \([^,)]*" | sed "s/gamma[^ ]*//" | sort | uniq -c
      1 A_minus_B COUPLING_KERNEL (BRANCH=LAB_HELD
      2 A_minus_B COUPLING_KERNEL (BRANCH=MATERIAL_ADVECTED
      2 A_minus_B MU_THETA_OPERATOR (BRANCH=LAB_HELD
      2 A_minus_B MU_THETA_OPERATOR (BRANCH=MATERIAL_ADVECTED
      4 A_minus_B SLAB_OPERATOR (OBJECT=MASS_EVOLUTION_ROW
      4 A_minus_B SLAB_OPERATOR (OBJECT=MU_THETA
      4 A_minus_B SLAB_OPERATOR (OBJECT=THICKNESS_ROW
      4 A_minus_B SLAB_OPERATOR (OBJECT=U_MOMENTUM_ROWS

## Claim: CORE-family verdict tally (post-fix reconcile)
$ grep -cE "^MATCH|^FLAG|^NAMESPACE_INCOMPLETE" /home/trevnorris/.s11_build/S11c_b_reconcile_run.out  (per keyword below)
MATCH   2
FLAG    14
NAMESPC 12

## Claim (protected one-engine): WL gamma-DivGrad heads = PY omitted candidates 08/11 (kept unfolded)
$ grep -nE "DivGrad(Theta|Ew)" research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl | head
392:    gammaWidthDivGradTheta, gammaWidthShearGradTheta,
393:    gammaWidthDivGradEw, gammaWidthShearGradEw,
396:    gammaModulusUDiv, gammaModulusDivGradTheta,
397:    gammaModulusShearGradTheta, gammaModulusDivGradEw,
# selected reps {1,4,5,6,7,9,10,13} / omitted {2,3,8,11,12,14,15} recomputed independently by both
# round-2 & round-3 decision-list legs (see _measurements/S11c_b_adjudication_directive_review.md);
# PY-only 07/10 and WL-only DivGrad 08/11 remain unmapped in H.WL_TO_PY_RENAME.
