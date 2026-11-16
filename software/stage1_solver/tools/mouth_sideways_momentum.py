"""First-order x momentum through the w=0 entrance of the existing bulk test.

The two spectator directions in the strip have unit area. A ball uses the
full three-dimensional aperture measure. m=n0=u_w0=1 in this test.
"""
import sympy as sp

from sleeve_outlet_boundary_data import beta, c2, linear_equations, v, w, x


def main():
    eps, U, B, D = sp.symbols("eps U B D", real=True)
    a = sp.Symbol("a", positive=True)
    mode2 = sp.exp(-2*w)*sp.sin(x)
    mode3 = sp.exp(-3*w)*sp.sin(x)
    # B is the previous conversion example. D is the earlier bulk mode
    # that preserves local normal current at the reference section.
    dn = B*mode2 + D*(8*mode3 - 3*mode2)
    dphi = U*x + sp.Rational(2, 3)*B*mode2 + D*(3*mode3 - 2*mode2)
    params = {v: 1, c2: sp.Rational(35, 24), beta: sp.Rational(1, 24)}
    residuals = linear_equations(dn, dphi).subs(params).applyfunc(sp.simplify)
    print("bulk [continuity, phase] =", list(residuals))
    n = 1 + eps*dn
    phi = w + eps*dphi
    ux, uw = sp.diff(phi, x), sp.diff(phi, w)
    j = n*uw
    p = (params[c2]/5)*n**5

    def first(field):
        return sp.simplify(sp.diff(field, eps).subs(eps, 0))

    def gradient_stress(row, col):
        return -params[beta]*n*sp.diff(sp.log(n), row, col)

    flow = first(j*ux)
    stress = first(gradient_stress(x, w))
    total = sp.simplify(flow + stress)
    diagonal = first(n*ux**2 + p + gradient_stress(x, x))
    print("x-momentum [d_x Pi_xx, d_w Pi_xw] =",
          [sp.simplify(sp.diff(diagonal, x)), sp.simplify(sp.diff(total, w))])
    balance = sp.simplify(sp.diff(diagonal, x) + sp.diff(total, w))
    print("x-momentum divergence =", balance)
    readings = [sp.simplify(field.subs(w, 0)) for field in (first(j), flow, stress, total)]
    print("mouth [number, carried, density gradient, total x momentum] =", readings)
    print("deep [number, total x momentum] =", [sp.limit(field, w, sp.oo) for field in (first(j),total)])
    weights = {"strip": 1, "ball": sp.pi*(a**2-x**2)}
    for shape, weight in weights.items():
        integrals = [sp.simplify(sp.integrate(weight*field, (x, -a, a))) for field in readings]
        print(shape, "[number, carried, density gradient, total] =", integrals)
        if shape == "strip":
            print("prior full strip at a=pi =", [value.subs(a, sp.pi) for value in integrals])

    # Radial baseline and pure dipole redistribution on a spherical shell
    # within the entrance ball. mu is the cosine to the companion.
    mu, j0, j1, ur = sp.symbols("mu j0 j1 ur", real=True)
    weight_change = first((j0+eps*j1*mu)*ur*mu)
    print("angular means [intake change, redistribution times existing velocity] =",
          [sp.integrate(term,(mu,-1,1))/2 for term in (j1*mu,weight_change)])

    # Conditional reversible stress of the postulated square-gradient order
    # energy under a virtual shear (x,w) -> (x+shear*w,w). The order
    # field is carried as a scalar. Minus the energy derivative is the flux.
    shear, kappa, chi_x, chi_w = sp.symbols("shear kappa chi_x chi_w", real=True)
    new_gradient = sp.Matrix([[1,shear],[0,1]]).inv().T * sp.Matrix([chi_x,chi_w])
    energy = kappa*new_gradient.dot(new_gradient)/2
    print("reversible order-gradient mixed flux =",
          -sp.diff(energy,shear).subs(shear,0))
    if any(residuals) or balance != 0:
        raise ValueError("Nonzero bulk equation or momentum-balance residual.")


if __name__ == "__main__":
    main()
