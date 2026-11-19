"""Project the bare order-gradient stress on a radial cut of the feeding slab."""

import sympy as sp


def check(label, value, expected):
    value = sp.factor(value)
    print(f"{label} = {value}", flush=True)
    if sp.simplify(value - expected) != 0:
        raise ArithmeticError(label)


def main():
    kappa, delta = sp.symbols("kappa delta", positive=True)
    c = sp.symbols("c", real=True)
    m_r, h_r, g_r = sp.symbols("m_r h_r g_r", real=True)

    # Arbitrary mirror-symmetric slab profile C(w-m,h). At reflected points
    # C_w has opposite signs, C_h the same. Cauchy traction is -k C_w C_r.
    cw, ch = sp.symbols("Cw Ch", real=True)
    def traction(dw):
        return -kappa * dw * (-m_r*dw + h_r*ch)

    check("paired depth elements of common traction",
          traction(cw) + traction(-cw), 2*kappa*cw**2*m_r)

    # Isolated logistic interface: dc/dw=c(1-c)/delta on the rising face.
    # Use c as the depth coordinate; reversing a face preserves its squared norm.
    gradient = c*(1-c)/delta
    sigma = kappa*sp.integrate(gradient, (c, 0, 1))
    check("single-face translation stiffness", sigma, kappa/(6*delta))
    check("two separated faces, common translation stiffness", 2*sigma, kappa/(3*delta))

    # Let each order level shift by m + g (2c-1)^2. The half-order crossing
    # is insensitive to g. Work to first order about the undeformed profile.
    level_weight = (2*c-1)**2
    check("shape weight at the face marker", level_weight.subs(c, sp.Rational(1, 2)), 0)
    shape_stiffness = kappa*sp.integrate(gradient*level_weight, (c, 0, 1))
    check("shape to translation traction ratio", shape_stiffness/sigma, sp.Rational(1, 5))
    check("two-face traction per area for translation and shape",
          2*(sigma*m_r + shape_stiffness*g_r), 2*sigma*(m_r+g_r/5))


if __name__ == "__main__":
    main()
