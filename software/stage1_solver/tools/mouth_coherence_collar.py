"""Sector comparison for a spherical material collar with fixed exterior flux.

Radial gauge for the prescribed connection is A_r=0. The angular connection
has uniform degree m=n-C_h and angular Laplacian floor |m|/2 on the unit
sphere. Neumann material edges are free. The exterior energy k_e*n**2/R is
counted once. Only the trial fields are restricted to be constant in r.
"""

import sympy as sp


def angular_moments(f, z, degree):
    # In the north chart A_phi=m(1-z), z=(1+cos(theta))/2, dOmega/(4pi)=dz
    # after azimuthal integration. South-chart transition removes the phase
    # singularity at the zero of these sections.
    kinetic = z*(1-z)*sp.diff(f, z)**2 + degree**2*(1-z)/(4*z)*f**2
    return tuple(sp.factor(sp.integrate(sp.cancel(p), (z, 0, 1)))
                 for p in (kinetic, f**2, f**4))


def main():
    r, R, lam, K, u, g, k_e, eta, mu = sp.symbols(
        "r R Lambda K u g k_e eta mu0", positive=True)
    n, h, s = sp.symbols("n h s", integer=True)
    q, gain = sp.symbols("q gain", nonnegative=True)
    z = sp.Symbol("z", positive=True)
    checks = []
    radial_length, radial_volume, radial_inv = (
        sp.factor(sp.integrate(p, (r, R, lam*R))) for p in (1, r**2, r**-2)
    )
    radial_factor = sp.factor(radial_length**2/(radial_volume*radial_inv))
    energy_unit = 2*sp.pi*K**2*g**2*radial_inv/u
    y = sp.factor(k_e/(R*energy_unit))
    # q=<Σ|Psi_h|²> over angles at a given r and may vary freely with r.
    # Angular spectrum gives sum_h (|m_h|/2 - g_h)<|Psi_h|²> >= -gain*q.
    # Jensen gives <S²> >= q².  Discard the nonnegative radial kinetic term.
    # gain=max(0,max_h(g_h-|m_h|/2)).  All mixed fields obey the same bound.
    shell = 4*sp.pi*(-K*gain*q + u*r**2*q**2/2)
    shell_q = sp.solve(sp.diff(shell, q), q)[0]
    shell_min = sp.factor(shell.subs(q, shell_q))
    floor = sp.factor(k_e*n**2/R + sp.integrate(shell_min, (r, R, lam*R)))
    print("radial moments (length, volume, inverse square):",
          radial_length, radial_volume, radial_inv)
    print("all-profile shell floor polynomial, minimizing mass, minimum:",
          shell, shell_q, shell_min)
    print("sector floor/B, y, residual to y*n^2-(gain/g)^2:",
          sp.factor(floor/energy_unit), y, sp.factor(floor/energy_unit-y*n**2+(gain/g)**2))
    checks.append(sp.factor(floor/energy_unit-y*n**2+(gain/g)**2) == 0)

    # F is the squared amplitude of a radially constant angular section.
    T, M, P, F = sp.symbols("T M P F", positive=True)
    trial_mat = 4*sp.pi*(K*(T - g*M)*radial_length*F +
                         u*P*radial_volume*F**2/2)
    stationary = sp.factor(sp.solve(sp.diff(trial_mat, F), F)[0])
    trial_lowering = sp.factor(-trial_mat.subs(F, stationary)/energy_unit)
    print("constant-r trial optimal F for g*M>T:", stationary)
    print("constant-r trial lowering/B, radial factor, factorization residual:",
          trial_lowering, radial_factor,
          sp.factor(trial_lowering - radial_factor*(M-T/g)**2/P))
    checks.append(sp.factor(trial_lowering - radial_factor*(M-T/g)**2/P) == 0)

    aspect = sp.Rational(5, 4)
    strength = sp.Rational(3, 2)
    ratio = sp.Rational(4, 5)
    y0 = sp.Rational(2, 9)
    material = {g: strength, lam: aspect}
    calibration = sp.factor(sp.solve(y.subs(material)-y0, k_e)[0])
    reference = {g: strength, lam: aspect, k_e: calibration}
    y_hom = sp.factor(y.subs(reference))
    weight = radial_factor.subs(lam, aspect)
    print("fixed g, eta in split control, Lambda, k_e, B, y, radial factor:",
          strength, ratio, aspect, calibration, sp.factor(energy_unit.subs(material)),
          y_hom, weight)
    checks.append(y_hom == y0)

    # The fixed shape has a single regular zero at the south pole in its
    # regular gauge. Complex conjugation gives its opposite degree.
    shape = sp.sqrt(z)*(1 + sp.Rational(2, 5)*(1-z))
    unit_moments = angular_moments(shape, z, 1)
    angular_values = dict(zip((T, M, P), unit_moments))
    unit_drop = trial_lowering.subs(material).subs(angular_values)
    unit_upper = y_hom - unit_drop
    zero_floor = -((strength - sp.Integer(2)/2)/strength)**2
    two_floor = 4*y_hom - 1  # n²>=4, no local field can lower F/B below -1.
    print("unit trial T M P, positive g*M-T, R^2 F:",
          *unit_moments, strength*unit_moments[1]-unit_moments[0],
          sp.factor(R**2*stationary.subs(material).subs(angular_values)))
    print("all R (unit upper, zero lower, |n|>=2 lower, normal lower):",
          unit_upper, zero_floor, two_floor, 0)
    margins = (zero_floor-unit_upper, two_floor-unit_upper, -unit_upper)
    print("all R dimensionless margins (zero, |n|>=2, normal):", *margins)
    checks.extend(margin > 0 for margin in margins)
    checks.append((strength*unit_moments[1] - unit_moments[0]) > 0)

    degree = n + 2*h*s
    same_s_conjugate = degree.subs({n: -n, h: -h}, simultaneous=True)
    mirror_swap = degree.subs({s: -s, h: -h}, simultaneous=True)
    print("equal g, degree (original, conjugation+swap, mirror+swap):",
          degree, same_s_conjugate, mirror_swap)
    checks.extend((same_s_conjugate == -degree, mirror_swap == degree))
    for orientation, favored in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
        favored_c = -2*orientation*favored
        selected = sp.sign(favored_c)
        favored_gain = strength - abs(-selected - favored_c)/2
        weaker_gain = ratio*strength - abs(-selected + favored_c)/2
        gain_max = sp.Max(0, favored_gain, weaker_gain)
        opposite_floor = y_hom - (gain_max/strength)**2
        gap = opposite_floor - unit_upper
        print(f"split s={orientation} favored_h={favored} (trial n, opposite n,"
              " opposite gains, upper, opposite lower, gap):",
              selected, -selected, favored_gain, weaker_gain, unit_upper,
              opposite_floor, gap)
        checks.append(gap > 0)

    # C_h=0 admits a constant flat angular section at n=0. Its radial trial
    # need not be exact. The rough shell floor for every charged sector is
    # already above that trial, with either choice of gains.
    flat_trial = -trial_lowering.subs(material).subs({T: 0, M: 1, P: 1})
    charge_floor_flat = y_hom - 1
    print("C_h=0 at all R (n=0 upper, |n|>=1 lower, gap):",
          flat_trial, charge_floor_flat, charge_floor_flat-flat_trial)
    checks.append(charge_floor_flat > flat_trial)

    # Fixed physical L=R0/4, same k_e and all material parameters.
    # For neutral flux at C_h=+-2, f=z is a smooth degree-two section.
    # Its angular gradient is small enough to give a negative trial energy.
    radius_ratio = sp.Symbol("a", positive=True)  # R/R0
    fixed_aspect = 1 + 1/(4*radius_ratio)
    y_fixed = sp.factor(y.subs({lam: fixed_aspect, g: strength, k_e: calibration}))
    fixed_weight = sp.factor(radial_factor.subs(lam, fixed_aspect))
    zero_moments = angular_moments(z, z, 2)
    zero_upper = -trial_lowering.subs({lam: fixed_aspect, g: strength}).subs(
        dict(zip((T, M, P), zero_moments)))
    single_gain = (strength-sp.Rational(1, 2))/strength
    single_floor = sp.factor(y_fixed - single_gain**2)
    double_floor = sp.factor(4*y_fixed - 1)
    threshold = sp.solve(single_floor, radius_ratio)[0]
    # Express the inequalities on the entire half-line a=threshold+delta.
    delta = sp.Symbol("delta", nonnegative=True)
    shifted = {radius_ratio: threshold+delta}
    print("fixed L=R0/4, a=R/R0 (y, radial factor):", y_fixed, fixed_weight)
    print("fixed L n=0 T M P, positive g*M-T, n=0 upper/B:",
          *zero_moments, strength*zero_moments[1]-zero_moments[0],
          sp.factor(zero_upper))
    print("fixed L sufficient threshold, a=threshold+delta (unit lower, |n|>=2 lower):",
          threshold, sp.expand(single_floor.subs(shifted)), sp.expand(double_floor.subs(shifted)))
    checks.extend((strength*zero_moments[1] > zero_moments[0], zero_upper.is_negative,
                   sp.expand(single_floor.subs(shifted)).is_nonnegative,
                   sp.expand(double_floor.subs(shifted)).is_positive))

    # Positive bulk detuning +mu0*S suppresses even a matched flat component.
    # For arbitrary g_h<=g and r>=R, write r=R+dr and g_h=g-dg.
    dr, dg = sp.symbols("dr dg", nonnegative=True)
    mass_excess = sp.factor(mu - K*(g-dg)/(R+dr)**2 - (mu-K*g/R**2))
    mass_threshold = sp.sqrt(K*g/mu)
    print("bulk detuning (sufficient R threshold, worst coefficient at threshold,"
          " coefficient surplus for r=R+dr, g_h=g-dg):",
          mass_threshold, sp.simplify((mu-K*g/R**2).subs(R, mass_threshold)), mass_excess)
    checks.extend((mass_excess.is_nonnegative,
                   sp.simplify((mu-K*g/R**2).subs(R, mass_threshold)) == 0))
    print("g_h=0 both branches (shell floor, total floor, n=0 normal energy):",
          shell_min.subs(gain, 0), floor.subs(gain, 0).factor(), 0)

    if not all(checks):
        raise ArithmeticError("An emitted shell bound, calibration, control or margin failed.")


if __name__ == "__main__":
    main()
