import sympy as sp


def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


mu = sp.symbols("mu", real=True)
P2 = sp.legendre(2, mu)

# ----------------------------------------------------------------------------
# 1) Target support data carried forward from the passed 2PN support-sector fit
# ----------------------------------------------------------------------------
section("1) Carried-forward support targets")

targets = {
    (0, 0): sp.Rational(4, 45),
    (1, 1): sp.Rational(2, 7),
    (1, 0): sp.Rational(1, 4),
    (2, 0): sp.Rational(4, 9),
    (2, 1): sp.Rational(2, 3),
    (2, 2): sp.Rational(8, 3),
}

for key in targets:
    print(f"K_{key} = {targets[key]}")

print("\nWe will keep K_(0,0) = 4/45 as the separate monopole wall channel, and solve")
print("the exact local traction-balance completion on the l=1,2 support sector.")

# ----------------------------------------------------------------------------
# 2) Minimal local traction-balance ansatz on the l=1,2 support sector
# ----------------------------------------------------------------------------
section("2) Minimal local traction-balance ansatz")

b0, b2, t0, t2, t4 = sp.symbols("b0 b2 t0 t2 t4", real=True)

z_base = b0 + b2 * mu**2
# Tangential wall-stress / curvature-compliance profile.
t_prof = t0 + t2 * mu**2 + t4 * mu**4

print("Ansatz:")
print("  z_base(mu) = b0 + b2 mu^2")
print("  t(mu)      = t0 + t2 mu^2 + t4 mu^4")
print("")
print("Channel formula on the l=1,2 support sector:")
print("  K_(l,m) = <z_base>_(l,m) + (l(l+1)-2) <t>_(l,m)")
print("This is the diagonal matrix-element form of the previously solved")
print("generalized Robin / support operator, rewritten as a local")
print("traction-balance completion.")


def density_lm(ell: int, m: int) -> sp.Expr:
    P = sp.assoc_legendre(ell, m, mu)
    norm = sp.Rational(2 * ell + 1, 2) * sp.factorial(ell - m) / sp.factorial(ell + m)
    return sp.simplify(norm * P**2)


eqs = []
for ell, m in [(1, 1), (1, 0), (2, 0), (2, 1), (2, 2)]:
    w = density_lm(ell, m)
    eq = sp.Eq(
        sp.integrate(sp.expand(z_base * w), (mu, -1, 1))
        + (ell * (ell + 1) - 2) * sp.integrate(sp.expand(t_prof * w), (mu, -1, 1)),
        targets[(ell, m)],
    )
    eqs.append(eq)

sol_list = sp.solve(eqs, [b0, b2, t0, t2, t4], dict=True)
assert len(sol_list) == 1
sol = sol_list[0]

z_base_sol = sp.expand(z_base.subs(sol))
t_prof_sol = sp.expand(t_prof.subs(sol))

print("Exact solution:")
for sym in [b0, b2, t0, t2, t4]:
    print(f"  {sym} = {sp.simplify(sol[sym])}")

print("\nSolved profiles:")
print(f"  z_base(mu) = {z_base_sol}")
print(f"  t(mu)      = {t_prof_sol}")

print("\nLegendre-basis form:")
a0, a2 = sp.symbols("a0 a2")
sol_base_leg = sp.solve(
    sp.Poly(sp.expand(a0 + a2 * P2 - z_base_sol), mu).coeffs(), [a0, a2], dict=True
)[0]
print(
    "  z_base(mu) = "
    f"{sp.simplify(sol_base_leg[a0])} + ({sp.simplify(sol_base_leg[a2])}) P2(mu)"
)

c0, c2, c4 = sp.symbols("c0 c2 c4")
sol_t_leg = sp.solve(
    sp.Poly(sp.expand(c0 + c2 * P2 + c4 * P2**2 - t_prof_sol), mu).coeffs(),
    [c0, c2, c4],
    dict=True,
)[0]
print(
    "  t(mu)      = "
    f"{sp.simplify(sol_t_leg[c0])} + ({sp.simplify(sol_t_leg[c2])}) P2(mu) + ({sp.simplify(sol_t_leg[c4])}) P2(mu)^2"
)

# ----------------------------------------------------------------------------
# 3) Equivalent local wall-energy form
# ----------------------------------------------------------------------------
section("3) Equivalent local wall-energy form")

print("For the operator z_base + 1/2{L-2, t}, with L = -Delta_S, the local")
print("quadratic form on the l=1,2 support sector is")
print("  E = 1/2 ∫ dΩ [ p(mu) Psi^2 + t(mu) |∇_S Psi|^2 ]")
print("with")
print("  p(mu) = z_base(mu) - 2 t(mu) - 1/2 Delta_S t(mu)")

Delta_t = sp.expand(sp.diff((1 - mu**2) * sp.diff(t_prof_sol, mu), mu))
p_prof_sol = sp.expand(z_base_sol - 2 * t_prof_sol - sp.Rational(1, 2) * Delta_t)
print(f"  Delta_S t(mu) = {Delta_t}")
print(f"  p(mu)         = {p_prof_sol}")

sol_p_leg = sp.solve(
    sp.Poly(sp.expand(c0 + c2 * P2 + c4 * P2**2 - p_prof_sol), mu).coeffs(),
    [c0, c2, c4],
    dict=True,
)[0]
print(
    "  p(mu)         = "
    f"{sp.simplify(sol_p_leg[c0])} + ({sp.simplify(sol_p_leg[c2])}) P2(mu) + ({sp.simplify(sol_p_leg[c4])}) P2(mu)^2"
)

# ----------------------------------------------------------------------------
# 4) Verification against the passed support data
# ----------------------------------------------------------------------------
section("4) Verification")


def grad_theta_density(ell: int, m: int) -> sp.Expr:
    P = sp.assoc_legendre(ell, m, mu)
    dP = sp.diff(P, mu)
    norm = sp.Rational(2 * ell + 1, 2) * sp.factorial(ell - m) / sp.factorial(ell + m)
    return sp.simplify(norm * (1 - mu**2) * dP**2)

for ell, m in [(1, 1), (1, 0), (2, 0), (2, 1), (2, 2)]:
    w = density_lm(ell, m)
    gt = grad_theta_density(ell, m)
    gp = sp.expand(m**2 / (1 - mu**2) * w)
    val_local = sp.simplify(sp.integrate(p_prof_sol * w + t_prof_sol * (gt + gp), (mu, -1, 1)))
    val_operator = sp.simplify(
        sp.integrate(z_base_sol * w, (mu, -1, 1))
        + (ell * (ell + 1) - 2) * sp.integrate(t_prof_sol * w, (mu, -1, 1))
    )
    print(f"({ell},{m}) local residual    = {sp.simplify(val_local - targets[(ell, m)])}")
    print(f"({ell},{m}) operator residual = {sp.simplify(val_operator - targets[(ell, m)])}")

print("\nMonopole channel:")
print(f"  K_(0,0) is kept separate and fixed to {targets[(0,0)]}.")
print("  The local traction-balance completion derived here is the exact")
print("  l=1,2 support-sector completion that supplements that monopole wall mode.")

# ----------------------------------------------------------------------------
# 5) Sign structure / support interpretation
# ----------------------------------------------------------------------------
section("5) Sign structure")

roots_t = [sp.N(r, 16) for r in sp.nroots(sp.Poly(t_prof_sol, mu)) if abs(sp.im(r)) < 1e-12]
roots_p = [sp.N(r, 16) for r in sp.nroots(sp.Poly(p_prof_sol, mu)) if abs(sp.im(r)) < 1e-12]
print(f"Real roots of t(mu): {roots_t}")
print(f"Real roots of p(mu): {roots_p}")
print("")
print(f"t(0)  = {sp.N(t_prof_sol.subs(mu, 0), 16)}")
print(f"t(1)  = {sp.N(t_prof_sol.subs(mu, 1), 16)}")
print(f"p(0)  = {sp.N(p_prof_sol.subs(mu, 0), 16)}")
print(f"p(1)  = {sp.N(p_prof_sol.subs(mu, 1), 16)}")

# ----------------------------------------------------------------------------
# 6) Compact conclusion
# ----------------------------------------------------------------------------
section("6) Conclusion")

print("Main result:")
print("- The passed l=1,2 support-sector data admit an exact minimal local")
print("  traction-balance completion with one base pressure profile z_base(mu)")
print("  and one tangential wall-stress profile t(mu), both in the Family-1")
print("  low-order basis {1, mu^2, mu^4}.")
print("- The resulting local wall-energy form reproduces all dipole and")
print("  quadrupole support channels exactly, with the monopole K_(0,0)=4/45")
print("  kept as the separate isotropic wall channel already identified in")
print("  the earlier PDE steps.")
print("- This is the first explicit traction-balance / wall-stress completion")
print("  of the strict isotropic layer that is exact on the solved 2PN")
print("  support sector.")
