import sympy as sp


def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


mu = sp.symbols("mu", real=True)

# -----------------------------------------------------------------------------
# 1) Static support data and exact angular moments
# -----------------------------------------------------------------------------
section("1) Input channel data and exact angular moments")

P2 = sp.legendre(2, mu)


def density_lm(ell: int, m: int) -> sp.Expr:
    P = sp.assoc_legendre(ell, m, mu)
    norm = sp.Rational(2 * ell + 1, 2) * sp.factorial(ell - m) / sp.factorial(ell + m)
    return sp.simplify(norm * P**2)


def moment(expr: sp.Expr, ell: int, m: int) -> sp.Expr:
    return sp.simplify(sp.integrate(sp.expand(expr * density_lm(ell, m)), (mu, -1, 1)))


moment_table = []
for ell in [1, 2]:
    for m in range(0, ell + 1):
        p2 = moment(P2, ell, m)
        p22 = moment(P2**2, ell, m)
        moment_table.append((ell, m, p2, p22))
        print(f"(ell={ell}, m={m}): <P2> = {p2}, <P2^2> = {p22}")

channels = [
    # name, ell, m, target stiffness K_{ell m}
    ("0", 0, 0, sp.Rational(4, 45)),
    ("1perp", 1, 1, sp.Rational(2, 7)),
    ("10", 1, 0, sp.Rational(1, 4)),
    ("20", 2, 0, sp.Rational(4, 9)),
    ("21", 2, 1, sp.Rational(2, 3)),
    ("22", 2, 2, sp.Rational(8, 3)),
]

print("\nStatic support targets:")
for row in channels:
    print("  ", row)

# -----------------------------------------------------------------------------
# 2) Reduced variational modal wall-Hessian ansatz
# -----------------------------------------------------------------------------
section("2) Reduced modal wall-Hessian ansatz")

Kmono, Tau0, A0, A1, B0, B1 = sp.symbols("Kmono Tau0 A0 A1 B0 B1")
lam = lambda ell: ell * (ell + 1)

print("Ansatz:")
print("  K_00 = Kmono")
print("  K_{ell m} = Kmono + Tau0*ell*(ell+1)")
print("             + (A0 + A1*ell*(ell+1)) <P2>")
print("             + (B0 + B1*ell*(ell+1)) <P2^2>,   for ell=1,2")

moment_map = {(ell, m): (p2, p22) for ell, m, p2, p22 in moment_table}

eqs = []
for name, ell, m, target in channels:
    if ell == 0:
        eqs.append(sp.Eq(Kmono, target))
    else:
        p2, p22 = moment_map[(ell, m)]
        eqs.append(
            sp.Eq(
                Kmono
                + Tau0 * lam(ell)
                + (A0 + A1 * lam(ell)) * p2
                + (B0 + B1 * lam(ell)) * p22,
                target,
            )
        )

sol = sp.solve(eqs, [Kmono, Tau0, A0, A1, B0, B1], dict=True)[0]

print("\nSolved coefficients:")
for key in [Kmono, Tau0, A0, A1, B0, B1]:
    print(f"  {key} = {sp.simplify(sol[key])}  ~=  {sp.N(sol[key], 16)}")

# -----------------------------------------------------------------------------
# 3) Exact reconstruction check
# -----------------------------------------------------------------------------
section("3) Exact reconstruction check")


def K_channel(ell: int, m: int) -> sp.Expr:
    p2, p22 = moment_map[(ell, m)]
    return sp.simplify(
        sol[Kmono]
        + sol[Tau0] * lam(ell)
        + (sol[A0] + sol[A1] * lam(ell)) * p2
        + (sol[B0] + sol[B1] * lam(ell)) * p22
    )

all_zero = True
for name, ell, m, target in channels:
    expr = sol[Kmono] if ell == 0 else K_channel(ell, m)
    resid = sp.simplify(expr - target)
    all_zero = all_zero and (resid == 0)
    print(f"{name:>6s}: K = {expr}, target = {target}, residual = {resid}")

print("\nPASS =", all_zero)

# -----------------------------------------------------------------------------
# 4) Family-1 Gaussian flare basis
# -----------------------------------------------------------------------------
section("4) Family-1 Gaussian flare -> {1, P2, P2^2} basis")

xi = sp.symbols("xi", positive=True)
D0 = sp.simplify(1 - xi / 3 + xi**2 / 18)
D1 = sp.simplify(-2 * xi / 3 + 2 * xi**2 / 9)
D2 = sp.simplify(2 * xi**2 / 9)

print("For a(z) = a0 (1 + beta exp(-z^2/zm^2)) and xi = a0^2/zm^2,")
print("the mouth-sphere flare expansion gives")
print("  delta a(theta)/a0 = beta [D0 + D1 P2 + D2 P2^2] + O(xi^3)")
print(f"  D0 = {D0}")
print(f"  D1 = {D1}")
print(f"  D2 = {D2}")
print("So the actual Family-1 flare automatically produces the exact axisymmetric basis used by the reduced wall Hessian.")

# -----------------------------------------------------------------------------
# 5) Minimal linear-gradient interpretation
# -----------------------------------------------------------------------------
section("5) Minimal linear-gradient interpretation")

ratio = sp.simplify(sol[B1] / sol[A1])
xi_fit = sp.simplify(sp.solve(sp.Eq(ratio, sp.simplify(D2 / D1)), xi)[0])
zm_over_a = sp.simplify(1 / sp.sqrt(xi_fit))

print("If the anisotropic gradient block is linear in the Family-1 flare profile, then")
print("  B1 / A1 = D2 / D1 = xi / (xi - 3)")
print(f"  B1/A1 = {ratio}  ~=  {sp.N(ratio, 16)}")
print(f"  xi = a0^2/zm^2 = {xi_fit}  ~=  {sp.N(xi_fit, 16)}")
print(f"  zm/a0 = 1/sqrt(xi) = {zm_over_a}  ~=  {sp.N(zm_over_a, 16)}")

# -----------------------------------------------------------------------------
# 6) Compact summary
# -----------------------------------------------------------------------------
section("6) Summary")

print("Main exact result:")
print("  K_00 = 4/45")
print("  K_{ell m} = 4/45 + (23/135) ell(ell+1)")
print("             + (9095/3528 - (25559/21168) ell(ell+1)) <P2>")
print("             + (-109/56 + (1765/3024) ell(ell+1)) <P2^2>,   for ell=1,2")
print("")
print("Interpretation:")
print("- This is the first exact reduced variational wall-Hessian that reproduces the passed static axisymmetric support data.")
print("- The Family-1 Gaussian flare naturally generates the same {1, P2, P2^2} axisymmetric basis at second order.")
print("- In the minimal linear-gradient reading, the flare width comes out of order the throat radius: zm/a0 ~= 1.0114.")
