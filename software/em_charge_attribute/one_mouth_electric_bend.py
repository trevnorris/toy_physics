"""Static response to a distant, spherically uniform change in the electric scalar."""

import sympy as sp


def check(label, actual, expected):
    result = sp.factor(actual)
    print(f"{label} = {result}", flush=True)
    if sp.simplify(result - expected) != 0:
        raise ArithmeticError(label)


def main():
    a, gap, r, kappa = sp.symbols("a gap r kappa", positive=True)
    outer = a + gap
    mouth, bath, const, tail = sp.symbols("h b c d", real=True)
    profile = const + tail/r
    ends = sp.solve([profile.subs(r, a) - mouth, profile.subs(r, outer) - bath],
                    [const, tail], dict=True)[0]
    field = profile.subs(ends)
    energy = sp.integrate(2*sp.pi*kappa*r**2*sp.diff(field, r)**2, (r, a, outer))
    force = sp.factor(sp.diff(energy, mouth))
    stiffness = sp.diff(force, mouth)
    check("shell stiffness", stiffness, 4*sp.pi*kappa/(1/a - 1/outer))
    check("distant bath limit", sp.limit(stiffness, gap, sp.oo), 4*sp.pi*kappa*a)
    check("virtual work minus inner conormal force", force + 4*sp.pi*kappa*a**2
          * sp.diff(field, r).subs(r, a), 0)

    K, k = sp.symbols("K k", positive=True)
    offset = sp.symbols("h_ref")
    # Baseline-subtracted equations. All source strengths and the reference
    # for the core are kept fixed during the signed change in the bath.
    load = force.subs(kappa, K / sp.diff(stiffness, kappa)).factor()
    fixed_value = sp.Integer(0)
    fixed_load = sp.solve(load, mouth)[0]
    elastic = sp.solve(load + k*mouth, mouth)[0]
    for label, h, h_slope, load_slope in [
        ("fixed value", fixed_value, 0, -K),
        ("fixed net load", fixed_load, 1, 0),
        ("finite restoring slope", elastic, K/(K+k), -K*k/(K+k)),
    ]:
        reaction = load.subs(mouth, h)
        check(label + " bend slope", sp.diff(h, bath), h_slope)
        check(label + " load slope", sp.diff(reaction, bath), load_slope)
        for reading in (h, reaction):
            signed = (reading - reading.subs(bath, -bath))/2
            check(label + " half-pair residual",
                  signed.subs(bath, bath/2) - reading/2, 0)

    check("stiffness inferred from the nonzero bend",
          -load.subs(mouth, elastic)/elastic, k)
    moving_reference = sp.solve(load + k*(mouth - offset), mouth)[0]
    check("common reference shift bend", moving_reference.subs(offset, bath), bath)
    check("common reference shift load", load.subs(mouth, moving_reference).subs(offset, bath), 0)


if __name__ == "__main__":
    main()
