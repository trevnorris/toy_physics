"""Retarded tangential translation dipole of the homogeneous curl-elastic brane."""

import sympy as sp


def main():
    a, r, k, s, rho, wave_speed = sp.symbols("a r k s rho c_T", positive=True)
    theta = sp.Symbol("theta", real=True)
    d, amplitude = sp.symbols("d C", real=True)
    modulus = rho*wave_speed**2

    # Polar axis along the drive. The curl of an azimuthal vector potential
    # is divergence-free; decaying modified spherical waves select no influx.
    f = amplitude*sp.exp(-k*(r-a))*(1+k*r)/r**2
    ur = 2*f*sp.cos(theta)/r
    ut = -sp.diff(r*f, r)*sp.sin(theta)/r
    div = sp.simplify(sp.diff(r**2*ur, r)/r**2
                      + sp.diff(sp.sin(theta)*ut, theta)/(r*sp.sin(theta)))
    curl = sp.simplify((sp.diff(r*ut, r) - sp.diff(ur, theta))/r)
    if div != 0 or sp.simplify(curl+k**2*f*sp.sin(theta)) != 0:
        raise ArithmeticError("Exterior transverse wave residual")
    tangent = (ut/(-sp.sin(theta))).subs(r, a)
    normalization, = sp.solve(tangent-d, amplitude)
    # Traction on the exterior is mu*(curl u) cross (-e_r).
    traction_theta = -modulus*curl.subs(r, a)
    angular_weight = 2*sp.pi*sp.sin(theta)
    surface_force = sp.integrate(-traction_theta*sp.sin(theta)*a**2*angular_weight,
                                 (theta, 0, sp.pi))
    # Independent momentum inventory. Its radial integrand is an exact
    # derivative and r^2*f vanishes at infinity when k>0.
    radial = sp.integrate((ur*sp.cos(theta)-ut*sp.sin(theta))*r**2*angular_weight,
                          (theta, 0, sp.pi))
    moment_primitive = 8*sp.pi*r**2*f/3
    if sp.simplify(radial-sp.diff(moment_primitive, r)) != 0:
        raise ArithmeticError("Momentum integrand disagrees with boundary work")
    force_from_storage = rho*s**2*(sp.limit(moment_primitive, r, sp.oo)
                                    - moment_primitive.subs(r, a))
    balance = sp.simplify((surface_force-force_from_storage).subs(k, s/wave_speed))
    if balance != 0:
        raise ArithmeticError("Traction disagrees with the exterior momentum rate")
    ratio = sp.factor((surface_force/(d*s**2)).subs(amplitude, normalization)
                      .subs(k, s/wave_speed))
    slow = sp.limit(ratio, s, 0)
    tau, t = sp.symbols("tau t", positive=True)
    response = sp.factor((ratio/slow).subs(wave_speed, a/tau))
    ramp = sp.simplify(sp.inverse_laplace_transform(response/s, s, t))
    print("divergence and bulk wave residual =", div, sp.simplify(curl+k**2*f*sp.sin(theta)))
    print("projected traction minus exterior momentum rate =", balance)
    print("force per acceleration of prescribed field trace =", ratio)
    print("slow coefficient and normalized causal transfer =", slow, response)
    print("normalized force for a field acceleration step =", ramp)
    print("onset value, initial slope, late limit =", sp.limit(ramp,t,0,dir="+"),
          sp.limit(sp.diff(ramp,t),t,0,dir="+"), sp.limit(ramp,t,sp.oo))
    # At the rim the transverse field's normal displacement is solved, not
    # fixed by a rigid-body constraint on the solid angle.
    normal_ratio = sp.factor((ur/(d*sp.cos(theta))).subs(r,a).subs(amplitude,normalization))
    print("normal displacement divided by d cos(theta) =", normal_ratio, flush=True)
    if (sp.simplify(slow-8*sp.pi*rho*a**3/3) != 0
            or sp.limit(ramp,t,0,dir="+") != 0 or sp.limit(ramp,t,sp.oo) != 1
            or sp.simplify(sp.limit(sp.diff(ramp,t),t,0,dir="+")-1/tau) != 0):
        raise ArithmeticError("Limits of the ideal shear response failed")


if __name__ == "__main__":
    main()
