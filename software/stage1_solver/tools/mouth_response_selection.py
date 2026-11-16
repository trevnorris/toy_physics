"""Linear mouth traces of the existing flat, steady 0PN carrier and passive layer."""
import sympy as sp

from sleeve_outlet_boundary_data import beta, c2, linear_equations, v, w, x


def main():
    B, D, U, J, H = sp.symbols("B D U J H", real=True)
    z2 = sp.exp(-2*w)*sp.sin(x)
    z3 = sp.exp(-3*w)*sp.sin(x)
    dn = B*z2 + D*(8*z3 - 3*z2)
    dphi = U*x + 2*B*z2/3 + D*(3*z3 - 2*z2)
    params = {v: 1, c2: sp.Rational(35, 24), beta: sp.Rational(1, 24)}
    bulk = linear_equations(dn, dphi).subs(params).applyfunc(sp.simplify)
    dux, duw = sp.diff(dphi, x), sp.diff(dphi, w)
    print("bulk [continuity, Bernoulli] =", list(bulk))

    # Y=log[(1-chi)/chi], Y0=w-b, material derivative of Y equals one.
    dy = sp.integrate(-duw, w)
    dy -= sp.limit(dy, w, sp.oo)
    dy_w = sp.simplify(sp.diff(dy, w))
    print("order [delta u_w, d_w delta Y, sum] =", [duw, dy_w, sp.simplify(duw + dy_w)])
    print("deep [dn, du_x, du_w, dY] =", [sp.limit(f, w, sp.oo) for f in (dn, dux, duw, dy)])

    current = (dn + duw).subs(w, 0)
    sine_amplitudes = sp.Matrix([current, dn.subs(w, 0)]).subs(x, sp.pi/2)
    matching = sine_amplitudes.jacobian([B, D])
    print("mouth [J,H] =", list(sine_amplitudes))
    print("matching matrix, determinant, current-only nullspace =",
          matching.tolist(), matching.det(), [list(col) for col in matching[:1, :].nullspace()])
    trace_data = sp.Matrix([J, H])
    amplitudes = matching.inv()*trace_data
    chosen = {B: amplitudes[0], D: amplitudes[1]}
    print("matched [B,D] =", list(amplitudes))

    # Cauchy traction on +w is the negative of the stress momentum flux.
    traction = params[beta]*sp.diff(dn, x, w)
    lateral, stress = [sp.expand(f.subs(w, 0).subs(chosen)).collect(sp.cos(x))
                       for f in (dux, traction)]
    print("mouth [delta u_x, delta carrier traction_x] =", [lateral, stress])
    vel_cos = (lateral-U).coeff(sp.cos(x))
    flux_cos = -stress.coeff(sp.cos(x))
    relation = sp.simplify(flux_cos.subs(H, sp.solve(sp.Symbol("V")-vel_cos, H)[0]))
    print("stress momentum flux cosine coefficient in terms of V and J =", relation)

    # Conditional square-gradient term evaluated on this passive profile.
    # chi'_0 = d chi_0/dw = d chi/dY at Y0=w-b. No mechanical feedback here.
    slope, kappa = sp.symbols("slope kappa", real=True)
    order_flux = kappa*slope**2*sp.diff(dy, x).subs(w, 0)
    print("mouth conditional order stress momentum flux =",
          sp.factor(order_flux.subs(chosen)))
    check_match = (sine_amplitudes.subs(chosen) - trace_data).applyfunc(sp.simplify)
    print("matched [current,density] residuals =", list(check_match))
    if any(bulk) or sp.simplify(duw + dy_w) != 0 or any(check_match):
        raise ValueError("Bulk, order or entrance matching failed.")


if __name__ == "__main__":
    main()
