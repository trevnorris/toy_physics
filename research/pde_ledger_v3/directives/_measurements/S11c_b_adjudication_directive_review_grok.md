# Measurements — S11c-b adjudication directive review (Grok leg)

## Bridge A coefficient alignment
```
WL_theta2_energy_coeff = W_bg*bRho/2
PY_theta2_energy_coeff = B_rho_3*W_bg/(2*W0)
solve_B_rho_3 = [W0*bRho]
residual_under_B_rho_3_is_bRho_W0 = 0
WL_hom = W0*bRho/2 PY_hom = B_rho_3/2
hom_residual_under_ident = 0
```

## Commands
- Cited lines read: WL:472,479,1625-1629; PY:1135; spec SHARED_PHYSICS:99-102 (`B_ρ⁽³⁾ ≡ B_ρ W₀`).
- WL transcript tag scan (Derivative of thetaWave / widthBackground / pressure):
```
ENERGY_BASIS_VARIABLE {'d_theta': True, 'd_ew': True, 'd_widthBg': False, 'w1JetOne': True, 'deriv_pressure': False}
SLAB_OPERATOR {'d_theta': True, 'd_ew': True, 'd_widthBg': False, 'w1JetOne': True, 'deriv_pressure': False}
MU_THETA_OPERATOR {'d_theta': True, 'd_ew': True, 'd_widthBg': False, 'w1JetOne': True, 'deriv_pressure': False}
COUPLING_KERNEL {'d_theta': False, 'd_ew': False, 'd_widthBg': False, 'w1JetOne': True, 'deriv_pressure': False}
```

## Post-C.parse_wl_value forms (gate evidence object)
```
Derivative[1,0,0,0][thetaWave][xOne,xTwo,xThree,time] -> theta_d1
thetaWave[xOne,xTwo,xThree,time] -> theta
Derivative[1,0,0][widthBackground][xOne,xTwo,xThree] -> W_bg_d1
widthBackground[xOne,xTwo,xThree] -> W_bg
Derivative[1,0,0,0][longitudinalTestPotential][xOne,xTwo,xThree,time] -> psi_L_s11cb_d1
longitudinalTestPotential[xOne,xTwo,xThree,time] -> psi_L_s11cb
pressureUpper[xOne,xTwo,xThree,time] -> pressureUpper(xOne, xTwo, xThree, time)
```

## v1 FLAG family tally (from ~/.s11_build/S11c_b_reconcile_run.out)
```
FLAG COUPLING_KERNEL
FLAG COUPLING_KERNEL_TERM_ORIGINS
FLAG SLAB_OPERATOR
FLAG SLAB_OPERATOR_TERM_ORIGINS
FLAG MU_THETA_OPERATOR
FLAG ADMISSIBILITY_OPERATOR_OPERAND
FLAG ADMISSIBILITY_RESIDUAL
MATCH ADMISSIBILITY_SUPPORT_OPERAND
FLAG ENERGY_BASIS_VARIABLE
MATCH ENERGY_BASIS_COUNT
FLAG ENERGY_BASIS_NEW_INVARIANTS
FLAG ENERGY_BASIS_OMISSIONS
FLAG DIMENSIONS
FLAG HOMOGENEITY_BASE_OPERAND
FLAG HOMOGENEITY_CONTROL_OPERAND
FLAG HOMOGENEITY_RESIDUAL
NAMESPACE_INCOMPLETE CONTROL_FORM_ABLATED_OPERAND
NAMESPACE_INCOMPLETE CONTROL_FORM_BASE_OPERAND
NAMESPACE_INCOMPLETE CONTROL_FORM_RESIDUAL
NAMESPACE_INCOMPLETE CONTROL_INDEPENDENCE_BASE_OPERAND
NAMESPACE_INCOMPLETE CONTROL_INDEPENDENCE_CORRUPTED_OPERAND
NAMESPACE_INCOMPLETE CONTROL_INDEPENDENCE_RESIDUAL
NAMESPACE_INCOMPLETE REP_INVARIANCE_EULERIAN_OPERAND
NAMESPACE_INCOMPLETE REP_INVARIANCE_MATERIAL_OPERAND
NAMESPACE_INCOMPLETE REP_INVARIANCE_RESIDUAL
NAMESPACE_INCOMPLETE UNIFORM_LIMIT_RESIDUAL
NAMESPACE_INCOMPLETE UNIFORM_LIMIT_S11B_OPERAND
NAMESPACE_INCOMPLETE UNIFORM_LIMIT_S11CB_OPERAND
Function heads remaining after v1: ['HeldDiv']
```
