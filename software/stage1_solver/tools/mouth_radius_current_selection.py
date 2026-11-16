"""Check which weak inlet currents are compatible with the fixed-bath carrier exterior."""

import sympy as s


def main():
    eps, eta, phase, lap_eta, lap_phase = s.symbols("eps eta phi lap_eta lap_phi", real=True)
    g, speed, mass = s.symbols("g V m", positive=True)
    field = 1 + eps * (eta + s.I * phase)
    density = s.expand_complex(field * field.conjugate()).simplify()
    laplacian = eps * (lap_eta + s.I * lap_phase)
    equation = -laplacian / 2 + (density**4 - 1) * field / 4
    # Keep the real restoring boundary law. Supply a trial normal phase slope
    # of eps*V*g only to ask whether the exterior accepts that slope.
    normal = (-g * (1 - 1 / s.sqrt(density)) + s.I * eps * speed * g) * field
    current = s.im(field.conjugate() * normal)

    def linear(expression):
        return s.simplify(s.diff(expression, eps).subs(eps, 0))

    k, w, depth = s.symbols("k w H", positive=True)
    amplitude, alpha, beta = s.symbols("q alpha beta", real=True)
    mode = alpha * s.cosh(k*w) + beta * s.sinh(k*w)
    slope = s.diff(mode, w)
    constants = s.solve([slope.subs(w, 0), slope.subs(w, depth) - amplitude], [alpha, beta], dict=True)
    profile = mode.subs(constants[0])
    recovered = s.simplify(s.diff(profile, w).subs(w, depth))

    r, radius, reference = s.symbols("r a a_0", positive=True)
    weight = (1 - (r/radius)**2)**3
    throughput = s.factor(s.integrate(4*s.pi*r**2 * mass * linear(current).subs(g, weight),
                                      (r, 0, radius)))
    candidates = [throughput, throughput.subs(speed, speed * reference / radius)]
    reference_rates = [s.factor(rate.subs(radius, reference)) for rate in candidates]
    doubling = [s.factor(rate.subs(radius, 2*reference) / initial)
                for rate, initial in zip(candidates, reference_rates, strict=True)]

    print("linear carrier equation (real, imaginary) =", linear(equation).as_real_imag(), flush=True)
    print("linear normal derivative (real, imaginary) =", linear(normal).as_real_imag(), flush=True)
    print("lateral mode phi(w) =", profile, flush=True)
    print("recovered trial phase slope at the mouth =", recovered, flush=True)
    print("linear inward mass rate at unit ambient number density =", throughput, flush=True)
    print("mass rates at calibration radius =", reference_rates, flush=True)
    print("mass-rate doubling factors (constant V, V*a_0/a) =", doubling, flush=True)

    assert linear(equation).as_real_imag() == (2*eta - lap_eta/2, -lap_phase/2)
    assert linear(normal).as_real_imag() == (-g*eta, g*speed)
    assert s.simplify(s.diff(profile, w, 2) - k**2*profile) == 0
    assert s.diff(profile, w).subs(w, 0) == 0 and recovered == amplitude
    assert reference_rates[0] == reference_rates[1] and doubling == [8, 4]


if __name__ == "__main__":
    main()
