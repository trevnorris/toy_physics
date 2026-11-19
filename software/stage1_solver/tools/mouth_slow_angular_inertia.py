"""Isolated signed spherical cut. Slow translational inertia of specific sectors.

Compare the coefficient of the outgoing moving-cut force after subtracting
its static prestress, with the kinetic energy of the translated exterior
angular field. Inputs are the material flux, angular wave speed and modulus.
The open-throat carrier, brane deformation and interior are absent.
"""
import mpmath as mp
import sympy as sp


def integrals():
    r, R, C, c = sp.symbols("r R C c", positive=True)
    theta, q = sp.symbols("theta q", real=True)
    field = q/(4*sp.pi*r**2)
    def outside(density):
        radial = sp.integrate(density*r**2, (r, R, sp.oo))
        return sp.simplify(2*sp.pi*sp.integrate(radial*sp.sin(theta), (theta, 0, sp.pi)))
    energy = outside(C*field**2/2)
    # The curl of the azimuthal translation derivative sin(theta)*D(r)
    # has exactly the first spatial translation of the monopolar field.
    angle_per_X = field*sp.sin(theta)
    radial_curl = sp.diff(sp.sin(theta)*angle_per_X, theta)/(r*sp.sin(theta))
    polar_curl = -sp.diff(r*angle_per_X,r)/r
    curl_residuals = [sp.simplify(radial_curl + sp.diff(field,r)*sp.cos(theta)),
                      sp.simplify(polar_curl-field*sp.sin(theta)/r)]
    mass = outside(C/c**2*angle_per_X**2)
    z = sp.symbols("z", real=True)
    outgoing = -sp.exp(sp.I*z)*(1+sp.I/z)/z
    traction = sp.factor(-(1+z*sp.diff(outgoing,z)/outgoing))
    reactive_limit = sp.limit((1-traction)/z**2,z,0)
    from_force = (8*sp.pi*R**2/3)*(C/R)*(q/(4*sp.pi*R**2))**2*(R/c)**2*reactive_limit
    print(f"curl_translation_residuals={curl_residuals} exterior_energy={energy} "
          f"exterior_inertia={mass} traction_limit={reactive_limit} "
          f"ratio_Mc2_over_energy={sp.simplify(mass*c**2/energy)} "
          f"traction_minus_volume_mass={sp.simplify(from_force-mass)}")
    assert curl_residuals==[0,0] and sp.simplify(from_force-mass)==0


def direct(radius, ratio, sign):
    # Solve the full tagged material, spin slip and outgoing wave at finite rate.
    # This time retain the conservative static prestress in the force.
    C, speed = mp.mpf(4), 1/mp.sqrt(2)
    mass_density, area_modulus = mp.mpf(10), mp.mpf(".06")
    pin, exchange = mp.mpf(".1"), mp.mpf(".1")
    half_width = mp.sqrt(pin**2+(exchange/2)**2)
    gap = half_width+exchange/2
    omega = ratio*gap/10
    kR = omega*radius/speed
    b = sign/(2*radius**2)
    h, B = 4*radius**2*pin**2/half_width, 8
    inertia = h*(10/gap)**2
    wave = lambda z: mp.hankel1(mp.mpf("1.5"),z)/mp.sqrt(z)
    Z = -C/radius*(1+kR*mp.diff(wave,kR)/wave(kR))
    mass = mp.matrix([[Z+B,-B,0],[-B,B+h-inertia*omega**2,inertia*b*omega**2],
                      [0,inertia*b*omega**2,2*area_modulus/radius**2-(mass_density+inertia*b**2)*omega**2]])
    a,s,u = mp.lu_solve(mass,[b*(C/radius+B),-B*b,mass_density*omega**2])
    area = 4*mp.pi*radius**2
    torque = B*(a-b-s)
    force_field = -2*area/3*b*(C/radius*(a-b)+torque)
    force_layer = -2*area/3*omega**2*mass_density*(mp.mpf("1.5")+u)
    return -mp.re(force_field)/omega**2, -mp.re(force_layer)/omega**2


def main():
    integrals()
    mp.mp.dps = 55
    for R in (mp.mpf(5),mp.mpf(10)):
        field_limit = 8*(2*mp.pi)**2/(6*mp.pi*R)
        layer_limit = 4*mp.pi*R**2*10
        for ratio in map(mp.mpf, (".01", ".001", ".0001")):
            masses = direct(R, ratio, 1)
            mirror = direct(R, ratio, -1)
            deviations = [abs(m/benchmark-1) for m, benchmark in zip(masses,(field_limit,layer_limit))]
            print(f"R={float(R):g} freq_over_gap={float(ratio):g} "
                  f"M_field={float(masses[0]):.11g} limit_field={float(field_limit):.11g} "
                  f"M_layer={float(masses[1]):.11g} limit_layer={float(layer_limit):.11g} "
                  f"max_mirror_diff={float(max(abs(a-b) for a,b in zip(masses,mirror))):.3g}")
            assert masses==mirror
            if ratio<mp.mpf(".0002"):
                assert max(deviations)<mp.mpf("1e-7")


if __name__=="__main__":
    main()
