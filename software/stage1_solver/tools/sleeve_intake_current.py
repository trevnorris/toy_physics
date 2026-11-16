"""Material intake through a flat patch of a smooth 0PN carrier.

Hold the surface density, tangential phase and its time evolution fixed.
Allow the normal phase derivative in a local stationary solution to vary.
The far end and wall of a physical core are outside this calculation.
"""
import sympy as sp


def main():
    w, x, t = sp.symbols("w x t", real=True)
    m, hbar, K = sp.symbols("m hbar K", positive=True)
    k, mu, V = sp.symbols("k mu V", real=True)
    amplitude = sp.Symbol("a", positive=True)
    a1, a2, p1, p2 = sp.symbols("a1 a2 p1 p2", real=True)
    # The normal Taylor coefficients at an arbitrary smooth point suffice
    # because the carrier equation contains at most two normal derivatives.
    a = amplitude + a1*w + a2*w**2/2
    phase = p1*w + p2*w**2/2
    carrier_phase = sp.exp(sp.I*(k*x + phase - mu*t/hbar))
    psi = a*carrier_phase
    laplacian = sp.diff(psi, w, 2) + sp.diff(psi, x, 2)
    # h(n) = d[(K/4)n^5]/dn, with n=a^2.
    residual = (sp.I*hbar*sp.diff(psi, t) + hbar**2*laplacian/(2*m)
                - (V + 5*K*a**8/4)*psi)
    at_surface = {w: 0, x: 0, t: 0}
    real_eq, imag_eq = sp.expand(residual.subs(at_surface)).as_real_imag()
    complex_flux = sp.conjugate(psi)*sp.diff(psi, w)
    normal_current = sp.simplify(hbar*sp.im(complex_flux.subs(at_surface))/m)
    current_derivative = sp.simplify(hbar*sp.im(sp.diff(complex_flux, w).subs(at_surface))/m)
    q = sp.Symbol("j", positive=True)
    phase_slope, = sp.solve(normal_current-q, p1)
    curvature, = sp.solve(real_eq, a2)
    curvature = sp.expand(curvature.subs(p1, phase_slope))
    print("normal current =", normal_current)
    print("imaginary equation =", imag_eq)
    print("hbar/(2*a) times current derivative =", sp.expand(hbar*current_derivative/(2*amplitude)))
    print("difference =", sp.simplify(imag_eq - hbar*current_derivative/(2*amplitude)))
    print("normal phase gradient at fixed current =", phase_slope)
    print("normal amplitude curvature =", curvature)

    # Dimensionless examples, with a=1 and zero amplitude slope at the intake.
    at_intake = {hbar: 1, m: 1, K: sp.Rational(4, 5), mu: sp.Rational(3, 2), V: 0,
                 amplitude: 1, a1: 0}
    density_curvature = sp.diff(a**2, w, 2).subs(w, 0).subs(a2, curvature)
    cases = [("isolated", 0, sp.Rational(1, 2)),
             ("incident, continuation A", sp.Rational(1, 4), sp.Rational(1, 2)),
             ("incident, continuation B", sp.Rational(1, 4), sp.Rational(9, 16))]
    for name, tangential_k, j in cases:
        values = {**at_intake, k: tangential_k, q: j}
        slope = phase_slope.subs(values)
        measured_current = normal_current.subs(p1, slope).subs(at_intake)
        print(name, "[tangential k, inward current, density curvature] =",
              [tangential_k, measured_current, density_curvature.subs(values)])


if __name__ == "__main__":
    main()
