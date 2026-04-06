#!/usr/bin/env python3

import sympy as sp


def grad(f, coords):
    return sp.Matrix([sp.diff(f, c) for c in coords])


def div(v, coords):
    return sum(sp.diff(v[i], coords[i]) for i in range(3))


def curl(v, coords):
    x, y, z = coords
    return sp.Matrix(
        [
            sp.diff(v[2], y) - sp.diff(v[1], z),
            sp.diff(v[0], z) - sp.diff(v[2], x),
            sp.diff(v[1], x) - sp.diff(v[0], y),
        ]
    )


def laplacian_vec(v, coords):
    x, y, z = coords
    return sp.Matrix(
        [sp.diff(comp, x, 2) + sp.diff(comp, y, 2) + sp.diff(comp, z, 2) for comp in v]
    )


def main():
    print("=== PART 1: Cavity Aspect Ratio ===")
    x01, a, L, chi, beta_vac = sp.symbols("x01 a L chi beta_vac", positive=True)
    H = L * x01**2 + a**2 * sp.pi**2 / L + beta_vac * a**2 * L
    radial_condition = sp.simplify(sp.diff(H, a) / (2 * a))
    beta_from_radial = sp.simplify(-sp.pi**2 / L**2)
    chi_equation = sp.simplify(
        sp.diff(H, L).subs({beta_vac: beta_from_radial.subs(L, chi * a), L: chi * a})
    )
    aspect_ratio = sp.sqrt(2) * sp.pi / x01
    print("beta_vac from dH/da = 0:", beta_from_radial)
    print("Aspect ratio from dH/dL = 0:", aspect_ratio)
    print("Axial residual at that aspect ratio:", sp.simplify(chi_equation.subs(chi, aspect_ratio)))

    print("\n=== PART 2: Hierarchy Scaling ===")
    Gamma, rho0, kappa_m, kappa_q = sp.symbols(
        "Gamma rho0 kappa_m kappa_q", positive=True
    )
    m_g = kappa_m * rho0 * sp.pi * a**2 * L
    q = kappa_q * rho0 * sp.pi * a**2 * Gamma
    q2_over_m2 = sp.simplify((q**2 / m_g**2).subs(L, aspect_ratio * a))
    print("q^2 / m_G^2:", q2_over_m2)
    print("Exponent of a:", sp.degree(sp.together(q2_over_m2).as_numer_denom()[1], a) * -1)

    print("\n=== PART 3: Coulomb / Gauss ===")
    x, y, z, eps0 = sp.symbols("x y z eps0", positive=True)
    r = sp.sqrt(x**2 + y**2 + z**2)
    phi_c = q / (4 * sp.pi * eps0 * r)
    e_field = -grad(phi_c, (x, y, z))
    div_e = sp.simplify(div(e_field, (x, y, z)))
    print("phi_C:", phi_c)
    print("div(E) away from origin:", div_e)

    print("\n=== PART 4: Ampere Identity ===")
    t, c = sp.symbols("t c", positive=True)
    ax = sp.Function("Ax")(t, x, y, z)
    ay = sp.Function("Ay")(t, x, y, z)
    az = sp.Function("Az")(t, x, y, z)
    phi = sp.Function("phi")(t, x, y, z)
    Avec = sp.Matrix([ax, ay, az])
    coords = (x, y, z)
    B = curl(Avec, coords)
    E = -grad(phi, coords) - sp.diff(Avec, t)
    lhs = sp.simplify(curl(B, coords) - sp.diff(E, t) / c**2)
    rhs = sp.simplify(
        -(-sp.diff(Avec, t, 2) / c**2 + laplacian_vec(Avec, coords))
        + grad(div(Avec, coords) + sp.diff(phi, t) / c**2, coords)
    )
    component_checks = [sp.simplify(lhs[i] - rhs[i]) == 0 for i in range(3)]
    print("Componentwise identity holds:", component_checks)

    print("\n=== PART 5: Magnus / Lorentz Matching ===")
    ux, uy, uz, v0, b0 = sp.symbols("ux uy uz v0 b0", real=True)
    uvec = sp.Matrix([ux, uy, uz])
    khat = sp.Matrix([0, 0, 1])
    v_inf = sp.Matrix([v0, 0, 0])
    f_m = rho0 * Gamma * khat.cross(uvec - v_inf)
    f_m_u = rho0 * Gamma * khat.cross(uvec)
    f_l_mag = (q / L) * uvec.cross(b0 * khat)
    b0_solution = sp.solve(sp.Eq(f_m_u[0], f_l_mag[0]), b0)[0]
    print("u-dependent Magnus term:", f_m_u)
    print("Lorentz magnetic term:", f_l_mag)
    print("B0 required for matching:", b0_solution)
    print("Background Magnus term:", sp.simplify(f_m - f_m_u))


if __name__ == "__main__":
    main()
