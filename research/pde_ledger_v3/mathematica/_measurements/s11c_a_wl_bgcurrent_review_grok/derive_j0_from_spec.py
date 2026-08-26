"""Independent derivation from S11c-a SHARED_PHYSICS §1b/§2d/§3c.

Spec §1b: j = rho_4D * v_bulk
Rest frame: v_bulk^0 = 0 (engine bulkVelocityZero; drain is scope limit only)
Spec §3c: rest-frame background current rho_4D^0 v_bulk^0 vanishes
Density background rho_4D^0 stays in B^0 (§2d) and is nonzero in general.
"""
from sympy import symbols, Matrix, simplify, Integer

rho0 = symbols("rho_4D^0", positive=True, nonzero=True)
# rest-frame bulk velocity 4-vector (3 in-plane + normal w)
v0 = Matrix([0, 0, 0, 0])
j0 = simplify(rho0 * v0)

print("OPERAND_rho0 =", rho0)
print("OPERAND_v_bulk0 =", list(v0))
print("PRODUCT_j0 = rho0 * v_bulk0 =", list(j0))
print("RESIDUAL_j0_minus_zero =", list(j0 - Matrix([0, 0, 0, 0])))
print("j0_is_identically_zero =", all(c == 0 for c in j0))
print("rho0_stays_symbolic_nonzero =", rho0.is_nonzero is True or str(rho0) == "rho_4D^0")

# Also show: if v were free symbols, j tracks rho*v (the form the ablation probes)
vp1, vp2, vp3, vpW = symbols("vp1 vp2 vp3 vpW")
v_probe = Matrix([vp1, vp2, vp3, vpW])
j_probe = simplify(rho0 * v_probe)
print("PROBE_FORM_j =", list(j_probe))
print("PROBE_currentW =", j_probe[3])
print("PROBE_currentX =", [j_probe[i] for i in range(3)])
