"""Cycle and surface projections for a conditional trapped-shear mouth model."""

import sympy as sp


def emit(label, value, control):
    value = sp.factor(value)
    print(f"{label} = {value}", flush=True)
    if sp.simplify(value - control) != 0:
        raise ArithmeticError(label)


def main():
    phase, mu, flow, speed = sp.symbols("phase mu U V", real=True)
    a, kappa, rho, omega = sp.symbols("a kappa rho omega", positive=True)
    amplitude, drive, detuning, loss = sp.symbols("S h Delta L", real=True)
    cosine, sine = sp.symbols("bc bs", real=True)

    def average(expr):
        return sp.integrate(sp.expand(expr), (phase, 0, 2 * sp.pi)) / (2 * sp.pi)

    # First harmonic of y'' + 2 gamma y' + omega1^2 y; L = 2 gamma omega.
    y = cosine * sp.cos(phase) + sine * sp.sin(phase)
    residual = (omega**2 * sp.diff(y, phase, 2) + loss * sp.diff(y, phase)
                + (detuning + omega**2) * y - drive * amplitude * sp.cos(phase))
    harmonic = [average(2 * residual * basis) for basis in (sp.cos(phase), sp.sin(phase))]
    coeff = sp.solve(harmonic, (cosine, sine), dict=True)[0]
    emit("in-phase response per speed", coeff[cosine], drive * amplitude * detuning / (detuning**2 + loss**2))
    emit("quadrature response per speed", coeff[sine], drive * amplitude * loss / (detuning**2 + loss**2))

    pressure = average(kappa * (amplitude * sp.cos(phase) + flow * mu * y)**2)
    force = 2 * sp.pi * a**2 * sp.integrate(mu * pressure, (mu, -1, 1))
    emit("radial supporting pressure at zero flow", pressure.subs(flow, 0), kappa * amplitude**2 / 2)
    emit("axial load at zero flow", force.subs(flow, 0), 0)
    slope = sp.diff(force, flow)
    emit("axial load per speed", slope, 4 * sp.pi * a**2 * kappa * amplitude * cosine / 3)
    emit("load with quadrature alone", force.subs(cosine, 0), 0)
    emit("detuned axial slope", slope.subs(coeff),
         4 * sp.pi * a**2 * kappa * drive * amplitude**2 * detuning / (3 * (detuning**2 + loss**2)))

    # One polarized displacement component; components add for a vector field.
    f, fx = sp.symbols("f fx", real=True)
    velocity = -omega * f * sp.sin(phase) - speed * fx * sp.cos(phase)
    kinetic = average(rho * velocity**2 / 2)
    emit("translation kinetic Hessian per volume", sp.diff(kinetic, speed, 2), rho * fx**2 / 2)


if __name__ == "__main__":
    main()
