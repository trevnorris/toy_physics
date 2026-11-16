"""Leading radial acceleration of a free draining mouth in a single velocity monopole."""

import sympy as sp


def main():
    r, density = sp.symbols("r rho", positive=True)
    mass_a, mass_b, mass_s = sp.symbols("M_A M_B M_s", positive=True)
    rate_a, rate_b, rate_s = sp.symbols("D_A D_B D_s", positive=True)
    far_a, far_b, far_s = sp.symbols("Q_A Q_B Q_s", positive=True)
    transfer_a, transfer_b = sp.symbols("e_A e_B", real=True, nonzero=True)

    def acceleration(transfer, target_rate, target_mass, far_source_rate):
        speed = far_source_rate / (4 * sp.pi * r**2)
        retained_force = transfer * target_rate * speed
        return retained_force / target_mass

    a_from_s = acceleration(transfer_a, rate_a, mass_a, far_s)
    b_from_s = acceleration(transfer_b, rate_b, mass_b, far_s)
    # The third body's target properties cancel, so any nonzero target response suffices.
    reference_response = sp.symbols("C_over_M", positive=True)
    at_reference_from_a = reference_response * far_a / (4 * sp.pi * r**2)
    at_reference_from_b = reference_response * far_b / (4 * sp.pi * r**2)
    target_ratio = sp.factor(a_from_s / b_from_s)
    source_ratio = sp.factor((at_reference_from_a / mass_a) / (at_reference_from_b / mass_b))
    g_pair = sp.factor(r**2 * a_from_s / mass_s)
    g_matched = g_pair.subs(far_s, rate_s / density)
    flux_per_mass = sp.symbols("lambda", positive=True)
    g_universal = sp.factor(g_matched.subs({rate_a: flux_per_mass * mass_a,
                                           rate_s: flux_per_mass * mass_s}))
    matched_source_ratio = source_ratio.subs({far_a: rate_a / density, far_b: rate_b / density})
    transfer_ratio = sp.factor(target_ratio / matched_source_ratio)

    for label, value in (
        ("acceleration A toward s", a_from_s),
        ("target acceleration ratio A/B at one source", target_ratio),
        ("source acceleration per unit source mass ratio A/B at one target", source_ratio),
        ("pair value of r^2 a / M_s", g_pair),
        ("target ratio divided by source ratio with matched flux", transfer_ratio),
        ("G with matched flux and common intake per mass", g_universal),
    ):
        print(f"{label} = {value}", flush=True)

    assert sp.simplify(target_ratio - (transfer_a * rate_a / mass_a) / (transfer_b * rate_b / mass_b)) == 0
    assert sp.simplify(source_ratio - (far_a / mass_a) / (far_b / mass_b)) == 0
    assert sp.simplify(transfer_ratio - transfer_a / transfer_b) == 0
    assert sp.simplify(g_universal - transfer_a * flux_per_mass**2 / (4 * sp.pi * density)) == 0


if __name__ == "__main__":
    main()
