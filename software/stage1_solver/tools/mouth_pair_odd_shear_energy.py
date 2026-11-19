"""Odd thickness displacement and opposite-face curl wave, fixed positive moduli.

Time-independent quadratic energy and unit-volume, lossless plane waves in a
homogeneous slab. The spring (v - lambda curl A)^2 is counted once. The term
"ordinary mechanical" refers only to rho vdot^2 and G |grad v|^2.
"""

import sympy as sp


rho, G, spring, b, K = sp.symbols("rho G kappa b K", positive=True)
lam = sp.symbols("lambda", real=True)


def real_space_current():
    # Keeping both transverse Cartesian components in a real wave makes the
    # current sign and the work done by the spring independent of helicity.
    x, t = sp.symbols("x t", real=True)
    vy, vz, ay, az = [sp.Function(name)(x, t) for name in ("vy", "vz", "Ay", "Az")]
    fields = (vy, vz, ay, az)
    by, bz = -az.diff(x), ay.diff(x)
    v_grad = G * (vy.diff(x)**2 + vz.diff(x)**2) / 2
    layer = spring * ((vy - lam*by)**2 + (vz - lam*bz)**2) / 2
    potential = v_grad + layer + K * (by**2 + bz**2) / 2
    inertia = (rho, rho, b, b)
    kinetic = sum(m * f.diff(t)**2 / 2 for m, f in zip(inertia, fields))
    acceleration = {f.diff(t, 2): (sp.diff(potential, f.diff(x)).diff(x)
                                    - sp.diff(potential, f)) / m
                    for m, f in zip(inertia, fields)}
    current = -sum(f.diff(t) * sp.diff(potential, f.diff(x)) for f in fields)
    ordinary = -sum(G * f.diff(x) * f.diff(t) for f in (vy, vz))
    hy, hz = [K * B - spring*lam*(v - lam*B) for v, B in ((vy,by), (vz,bz))]
    residual = ((kinetic + potential).diff(t) + current.diff(x)).subs(acceleration).expand()
    flux_residual = (current - ordinary - (hy*az.diff(t) - hz*ay.diff(t))).expand()
    print("local continuity residual and J - Jv - (H x Adot)_x:", residual, flux_residual)
    assert residual == flux_residual == 0


def wave():
    k = sp.symbols("k", positive=True)
    z = sp.symbols("z", real=True)          # omega^2, includes zero at k=0
    v, a = sp.symbols("v A", real=True)      # normalized helicity polarization
    mass = sp.diag(rho, b)
    helicities = {}
    for h in (1, -1):
        flux = h*k*a
        potential2 = G*k**2*v**2 + K*flux**2 + spring*(v - lam*flux)**2
        stiffness = sp.hessian(potential2/2, (v, a))
        dynamic = stiffness - z*mass
        determinant = dynamic.det().expand()
        ratio = sp.factor(-dynamic[0,1] / (dynamic[0,0]*h*k))
        # Eigenmode with A=1. r = v/B, so v = h k r.
        vec = sp.Matrix([h*k*ratio, 1])
        # Energy transport from the real-space current in a circular plane wave
        # of phase exp[i(kx-omega t)]. Jv and J are cycle averages divided by
        # omega*k/2; energy terms below are four times their cycle averages.
        jv = G*vec[0]**2
        jt = jv + K - spring*lam*(ratio - lam)
        flux_check = sp.factor(jt - (vec.T*stiffness.diff(k)*vec)[0]/(2*k))
        print(f"helicity {h:+d} v/B, flux derivative residual:", ratio, flux_check)
        assert flux_check == 0
        helicities[h] = determinant, ratio, vec
    determinant, ratio, vec = helicities[1]
    assert determinant == helicities[-1][0] and ratio == helicities[-1][1]

    # The determinant has two simple roots at k=0. Expand the gapless one,
    # solving the k^2 and k^4 coefficients before taking energy/current ratios.
    c2, c4 = sp.symbols("c2 c4")
    trial = determinant.subs(z, c2*k**2 + c4*k**4).expand()
    speed2 = sp.solve(trial.coeff(k,2), c2)[0]
    correction = sp.factor(sp.solve(trial.coeff(k,4).subs(c2,speed2), c4)[0])
    low = speed2*k**2 + correction*k**4
    print("squared frequencies at k=0, gapless omega^2/k^2 and k^4 coefficient:",
          sp.solve(determinant.subs(k,0),z), speed2, correction)
    print("mechanical displacement v/B at k=0 and its k^2 correction:",
          sp.limit(ratio.subs(z,low), k,0),
          sp.factor(sp.limit((ratio.subs(z,low)-lam)/k**2, k,0)))

    em = (rho*z + G*k**2)*vec[0]**2
    ef = b*z + K*k**2
    es = spring*(vec[0] - lam*k)**2
    et = em + ef + es
    jv = G*vec[0]**2
    jt = jv + K - spring*lam*(ratio - lam)
    # For a normal mode, four times the time-averaged energy equals
    # 2 z (rho v^2 + b A^2). Off shell the mismatch is det(D)/D_vv.
    energy_residual = sp.factor(et - 2*z*(rho*vec[0]**2 + b)
                                - determinant/(spring + G*k**2 - rho*z))
    print("kinetic-potential equipartition residual after shell subtraction:", energy_residual)
    e_coeff, j_coeff = [sp.factor(sp.limit((part/whole/k**2).subs(z,low), k,0))
                        for part, whole in ((em,et), (jv,jt))]
    s_coeff = sp.factor(sp.limit((es/et/k**4).subs(z,low), k,0))
    print("lim ordinary energy fraction/k^2, ordinary current fraction/k^2,"
          " spring energy fraction/k^4:", e_coeff, j_coeff, s_coeff)
    equal = {G: rho*K/b, z: K*k**2/b}
    matched = [sp.factor(f.subs(equal)) for f in (determinant, ratio, em/et, jv/jt, es)]
    print("matched speeds: dispersion residual, v/B, energy fraction, current fraction,"
          " spring energy:", matched)
    assert energy_residual == matched[0] == matched[-1] == 0
    assert sp.factor(matched[2] - matched[3]) == 0

    # Solve the exact lower generalized eigenvalue at finite wavenumber, with
    # real positive A on the + helicity at x=t=0 defining the phase reference.
    sample = {rho:3, G:5, spring:4, lam:sp.Rational(1,2), b:2, K:3}
    roots = sp.solve(determinant.subs(sample), z)
    print("sample rho,G,kappa,lambda,b,K=3,5,4,1/2,2,3; + helicity, A=1 real")
    print("k omega v/B ordinary_E/total_E ordinary_J/total_J scaled_E scaled_J"
          " Jtotal/Etotal group_velocity")
    for kval in (sp.Rational(1,2), sp.Rational(1,8), sp.Rational(1,32), sp.Rational(1,128)):
        root = min(roots, key=lambda r: r.subs(k, kval).evalf())
        at = {**sample, k:kval, z:root.subs(k,kval)}
        omega = sp.sqrt(at[z])
        er, jr = [(p/q).subs(at) for p,q in ((em,et), (jv,jt))]
        velocity = (2*k*sp.sqrt(z)*jt/et).subs(at)
        group_velocity = sp.diff(sp.sqrt(root),k).subs(k,kval)
        values = (kval, omega, ratio.subs(at), er, jr, er/kval**2, jr/kval**2,
                  velocity, group_velocity)
        print(*(f"{float(value):.12g}" for value in values))
        assert abs(float(velocity-group_velocity)) < 1e-12
    print("sample scaled fraction limits:", e_coeff.subs(sample), j_coeff.subs(sample))


if __name__ == "__main__":
    real_space_current()
    wave()
