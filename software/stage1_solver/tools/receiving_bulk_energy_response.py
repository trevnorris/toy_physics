"""Fixed 4D cut in a potential receiver with a separation-dependent dipole."""
import sympy as sp


def check(label, computed, expected):
    reduced = sp.factor(sp.simplify(computed))
    residual = sp.simplify(reduced - expected)
    print(f"{label}: {reduced} | expected {expected} | residual {residual}", flush=True)
    if residual != 0:
        raise ArithmeticError(label)


def main():
    print("SymPy version =", sp.__version__, flush=True)
    r, d, rho, q1, q2 = sp.symbols("R d rho q1 q2", positive=True)
    z, s, slope = sp.symbols("z s k", real=True)
    omega = 2 * sp.pi**2
    # z is the cosine toward the companion. Outside the mouth the zonal
    # angular Laplacian on the unit 3-sphere is (1-z^2) f_zz - 3 z f_z.
    h = z / r**3
    lap = (sp.diff(r**3 * sp.diff(h, r), r) / r**3
           + ((1-z**2)*sp.diff(h, z, 2)-3*z*sp.diff(h, z)) / r**2)
    check("exterior dipole Laplacian", lap, 0)
    weight = 4 * sp.pi * sp.sqrt(1-z**2)
    check("change in total flux per unit dipole amplitude",
          r**3 * sp.integrate(weight * sp.diff(h, r), (z, -1, 1)), 0)
    check("direction cosine squared, 4D surface mean",
          sp.integrate(weight*z**2, (z, -1, 1))/omega, sp.Rational(1, 4))

    # By divergence, integral(phi2 * normal_x dA) = integral(d_x phi2 dV).
    # The harmonic volume mean gives the center derivative throughout R<d.
    # Differentiate its *angular* moment at fixed directions for the flux.
    phi2 = -q2 / (2*rho*omega*(d**2 - 2*d*r*z + r**2))
    axial_gradient = sp.diff(phi2, r).subs({r: 0, z: 1})
    ball_volume = sp.integrate(omega*r**3, (r, 0, r))
    surface_velocity_moment = sp.diff(ball_volume*axial_gradient/r**3, r)
    check("integral of h times companion normal velocity",
          surface_velocity_moment, omega*axial_gradient/4)

    # Compare Bernoulli energy current in two boundary responses at b(d*)=0.
    # At fixed x relative to mouth 1, phi_t=-s phi_d. The same velocities
    # occur at the comparison event for any b'(d*)=k.
    distance = sp.Symbol("distance", positive=True)
    b = sp.Function("b")(distance)
    dt_added = (-s*sp.diff(b*h, distance)).subs(sp.diff(b, distance), slope)
    normal, kinetic, base_dt, head = sp.symbols("vN E phi_t head", real=True)
    p = head - kinetic - rho*base_dt
    change = sp.expand(((kinetic+p-rho*dt_added) - (kinetic+p))*normal)
    check("change in local power density", change, rho*s*slope*h*normal)
    radial_self = q1/(rho*omega*r**3)
    self_moment = r**3*sp.integrate(weight*h*radial_self, (z, -1, 1))
    check("self-flow projection in power change", self_moment, 0)
    delta_slope = rho*slope*(self_moment+surface_velocity_moment)
    check("change in four-state odd slope", delta_slope, -slope*q2/(4*d**3))
    mono_slope = -q1*axial_gradient
    check("pure monopole four-state odd slope", mono_slope, q1*q2/(rho*omega*d**3))
    check("ratio for k=q1/(rho Omega4)",
          ((mono_slope+delta_slope)/mono_slope).subs(slope, q1/(rho*omega)), sp.Rational(3, 4))
    # For the added potential only, extract the amplitude from normal velocity.
    response = sp.Symbol("B", real=True)
    measured = sp.integrate(weight*z*sp.diff(response*h, r), (z,-1,1))/omega
    check("amplitude recovered from normal-velocity moment", -4*r**4*measured/3, response)


if __name__ == "__main__":
    main()
