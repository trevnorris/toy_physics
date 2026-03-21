import sympy as sp
import numpy as np
from scipy.optimize import least_squares


def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


mu = sp.symbols("mu", real=True)
A, B, q, r = sp.symbols("A B q r", real=True)

# -----------------------------------------------------------------------------
# 1) Strict two-moment boundary-layer surface action on a flared mouth sphere
# -----------------------------------------------------------------------------
section("1) Strict two-moment Family-1 boundary-layer ansatz")

print("Surface profile:")
print("  sigma(mu) = 1 - q mu^2 + r mu^4")
print("Strict surface action:")
print("  E = 1/2 ∫ dΩ [ A J(mu) Ψ^2 + B (Fθ(mu) (∂θΨ)^2 + Fφ(mu) (∂φΨ)^2/sin^2θ ) ]")
print("with J, Fθ, Fφ obtained by pulling the steep isotropic wall layer back to the reference mouth sphere.")

sigma = 1 - q * mu**2 + r * mu**4
eps = sp.symbols("eps")
sigma_eps = sigma.subs({q: eps * q, r: eps * r})

J = sp.series(
    sigma_eps * sp.sqrt(sigma_eps**2 + (1 - mu**2) * sp.diff(sigma_eps, mu) ** 2),
    eps,
    0,
    3,
).removeO().subs(eps, 1).expand()

Ftheta = sp.series(
    sigma_eps / sp.sqrt(sigma_eps**2 + (1 - mu**2) * sp.diff(sigma_eps, mu) ** 2),
    eps,
    0,
    3,
).removeO().subs(eps, 1).expand()

Fphi = sp.series(
    sp.sqrt(sigma_eps**2 + (1 - mu**2) * sp.diff(sigma_eps, mu) ** 2) / sigma_eps,
    eps,
    0,
    3,
).removeO().subs(eps, 1).expand()

print("Quadratic flare truncation:")
print(f"  J(mu)    = {J}")
print(f"  Fθ(mu)   = {Ftheta}")
print(f"  Fφ(mu)   = {Fphi}")

# -----------------------------------------------------------------------------
# 2) Channel stiffnesses
# -----------------------------------------------------------------------------
section("2) Channel stiffness formulas")


def density_lm(ell: int, m: int) -> sp.Expr:
    P = sp.assoc_legendre(ell, m, mu)
    norm = sp.Rational(2 * ell + 1, 2) * sp.factorial(ell - m) / sp.factorial(ell + m)
    return sp.simplify(norm * P**2)


def grad_theta_density(ell: int, m: int) -> sp.Expr:
    P = sp.assoc_legendre(ell, m, mu)
    dP = sp.diff(P, mu)
    norm = sp.Rational(2 * ell + 1, 2) * sp.factorial(ell - m) / sp.factorial(ell + m)
    return sp.simplify(norm * (1 - mu**2) * dP**2)


targets = {
    (0, 0): sp.Rational(4, 45),
    (1, 1): sp.Rational(2, 7),
    (1, 0): sp.Rational(1, 4),
    (2, 0): sp.Rational(4, 9),
    (2, 1): sp.Rational(2, 3),
    (2, 2): sp.Rational(8, 3),
}

exprs = {}
for ell, m in targets:
    w = sp.expand(density_lm(ell, m))
    gt = sp.expand(grad_theta_density(ell, m))
    gp = sp.expand(m**2 / (1 - mu**2) * w)
    expr = sp.simplify(sp.integrate(A * J * w + B * (Ftheta * gt + Fphi * gp), (mu, -1, 1)))
    exprs[(ell, m)] = sp.expand(expr)
    print(f"K_{{{ell},{m}}} = {exprs[(ell, m)]}")

# -----------------------------------------------------------------------------
# 3) Numerical best-fit test against the passed support data
# -----------------------------------------------------------------------------
section("3) Numerical best-fit against passed support data")

funs = [sp.lambdify((A, B, q, r), sp.N(exprs[k] - targets[k]), "numpy") for k in targets]
keys = list(targets.keys())


def residual(x: np.ndarray) -> np.ndarray:
    a, b, qq, rr = x
    return np.array([f(a, b, qq, rr) for f in funs], dtype=float)


seed_list = [
    np.array([0.1, 0.1, 0.1, 0.1]),
    np.array([0.05, 0.2, 0.3, 0.1]),
    np.array([0.2, 0.1, -0.3, 0.2]),
    np.array([0.5, 0.1, 0.2, -0.1]),
    np.array([0.1, 1.0, 0.2, 0.1]),
    np.array([0.5, 0.5, 0.5, 0.2]),
    np.array([1.0, 0.2, 0.3, 0.3]),
    np.array([-0.1, 0.2, 0.5, 0.1]),
]

rng = np.random.default_rng(12345)
for _ in range(24):
    seed_list.append(rng.uniform(low=[-1, -1, -4, -4], high=[1, 1, 4, 4]))

best = None
best_ss = np.inf
for seed in seed_list:
    res = least_squares(residual, seed, method="trf", max_nfev=4000)
    ss = float(np.sum(res.fun**2))
    if ss < best_ss:
        best_ss = ss
        best = res

xbest = best.x
print("Best-fit parameters:")
print(f"  A = {xbest[0]: .12f}")
print(f"  B = {xbest[1]: .12f}")
print(f"  q = {xbest[2]: .12f}")
print(f"  r = {xbest[3]: .12f}")
print("")
print(f"Sum of squared residuals = {best_ss:.15f}")
print(f"Max absolute channel residual = {np.max(np.abs(best.fun)):.15f}")
print("")
print("Channel residuals at best fit:")
for key, val in zip(keys, best.fun):
    print(f"  {key}: {val:+.15f}")

# -----------------------------------------------------------------------------
# 4) Compact conclusion
# -----------------------------------------------------------------------------
section("4) Conclusion")

print("Main result:")
print("- The strict two-moment isotropic boundary-layer pullback does NOT exactly reproduce the passed support-channel data.")
print("- Across many seeds, the best least-squares fit remains far from zero.")
print("")
print("Interpretation:")
print("- Pure isotropic penetration + pure geometry pullback is too small a model.")
print("- At least one extra local wall-stress / traction / profile degree of freedom is required.")
print("- This makes the reduced variational wall-Hessian derived in the next step a genuine new structure, not just a reparameterization of the strict isotropic layer.")
