"""Reduce a linear translation probe and check its fixed-time release readout.

The Laplace-domain entries are unspecified physical response operators. The
two-row check has one other field coordinate and assumes an invertible field
block and a regular coupled solve. No entry is a measured mouth coefficient.
"""

import sympy as sp


def main():
    s = sp.Symbol("s", positive=True)
    f, hold = sp.symbols("f H")
    xx, xy, yx, yy = sp.symbols("Dxx Dxy Dyx Dyy")
    # Subtract the continued held solution from the released solution.
    # The independent initial fields and maintained flow cancel; the
    # calibrated external load replaces the holder in the translation row.
    operator = sp.Matrix([[xx, xy], [yx, yy]])
    readout = operator.inv() * sp.Matrix([(f - hold) / s, 0])
    response = sp.factor(sp.diff(readout[0], f))
    stiffness = xx - xy * yx / yy
    print("center displacement per force step in Laplace space =", response)
    print("reapplied holder displacement =", readout[0].subs(f, hold))
    residual = sp.factor(s * stiffness * readout[0] - (f - hold))
    print("reduced center equation residual =", residual)

    # At a fixed elapsed time after release, differentiate the condition
    # X(t; H(U), U) = X0 directly, without using the operator inverse.
    u = sp.Symbol("U")
    t = sp.Symbol("t", positive=True)
    force_slope = sp.Function("S")(t)
    speed_slope = sp.Function("flow_slope")(t)
    c, h2, ff, fu, uu = sp.symbols("C h2 ff fu uu")
    travel = force_slope*f + speed_slope*u + ff*f**2 + fu*f*u + uu*u**2
    held_curve = travel.subs(f, -c*u + h2*u**2)
    flow_slope, = sp.solve(sp.diff(held_curve, u).subs(u, 0), speed_slope)
    print("fixed-time displacement derivative in maintained flow =", flow_slope)
    if residual != 0 or flow_slope != c * force_slope:
        raise ArithmeticError("The coupled solve and held-state differentiation must be consistent")


if __name__ == "__main__":
    main()
