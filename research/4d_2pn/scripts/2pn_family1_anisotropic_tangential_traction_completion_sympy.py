import sympy as sp
import numpy as np
from scipy.optimize import least_squares


def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


mu = sp.symbols("mu", real=True)
A, Btan0, Btan2, q, r = sp.symbols("A Btan0 Btan2 q r", real=True)

# -----------------------------------------------------------------------------
# 1) Minimal strict Family-1 completion
# -----------------------------------------------------------------------------
section("1) Minimal anisotropic tangential-traction completion")

print("Start from the strict Family-1 steep-wall pullback, but keep only the")
print("smallest new local wall degree of freedom that the previous no-go forced:")
print("")
print("  isotropic normal penetration moment   A")
print("  axisymmetric tangential wall moment   B_tan(mu) = Btan0 + Btan2 mu^2")
print("  flared mouth profile                  sigma(mu) = 1 - q mu^2 + r mu^4")
print("")
print("This is the smallest local PDE-side completion of the strict isotropic")
print("boundary-layer model that can possibly generate the exact l=1,2 support")
print("sector found in the earlier operator steps.")

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

Btan = Btan0 + Btan2 * mu**2

print("Quadratic flare truncation:")
print(f"  J(mu)    = {J}")
print(f"  Ftheta   = {Ftheta}")
print(f"  Fphi     = {Fphi}")
print(f"  B_tan(mu)= {Btan}")

# -----------------------------------------------------------------------------
# 2) Channel formulas
# -----------------------------------------------------------------------------
section("2) Exact channel formulas")


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
    (1, 1): sp.Rational(2, 7),
    (1, 0): sp.Rational(1, 4),
    (2, 0): sp.Rational(4, 9),
    (2, 1): sp.Rational(2, 3),
    (2, 2): sp.Rational(8, 3),
}

all_channels = [(0, 0)] + list(targets.keys())
exprs = {}

for ell, m in all_channels:
    w = sp.expand(density_lm(ell, m))
    gt = sp.expand(grad_theta_density(ell, m))
    gp = sp.expand(sp.simplify(m**2 / (1 - mu**2) * w))
    expr = sp.simplify(sp.integrate(A * J * w + Btan * (Ftheta * gt + Fphi * gp), (mu, -1, 1)))
    exprs[(ell, m)] = sp.expand(expr)
    print(f"K_{{{ell},{m}}} = {exprs[(ell, m)]}")

# -----------------------------------------------------------------------------
# 3) Multi-seed exact solve on the l=1,2 support sector
# -----------------------------------------------------------------------------
section("3) Multi-seed exact solve on the l=1,2 support sector")

keys = list(targets.keys())
f_syms = [exprs[k] - targets[k] for k in keys]
jac_syms = sp.Matrix(f_syms).jacobian([A, Btan0, Btan2, q, r])

f_num = sp.lambdify((A, Btan0, Btan2, q, r), f_syms, "numpy")
jac_num = sp.lambdify((A, Btan0, Btan2, q, r), jac_syms, "numpy")


def residual(x: np.ndarray) -> np.ndarray:
    return np.array(f_num(*x), dtype=float)


def jacobian(x: np.ndarray) -> np.ndarray:
    return np.array(jac_num(*x), dtype=float)


seed_list = [
    np.array([-0.28, 0.65, -1.09, 2.37, 2.76]),
    np.array([-0.15, 0.68, -1.18, 0.72, -1.62]),
    np.array([-0.15, 0.64, -1.05, -2.27, -2.30]),
    np.array([-0.086, 0.709, -1.266, -0.734, 1.250]),
]

rng = np.random.default_rng(20260321)
for _ in range(28):
    seed_list.append(rng.uniform(low=[-1.0, -0.5, -2.0, -3.5, -3.5], high=[0.5, 1.5, 0.5, 3.5, 3.5]))

raw_solutions = []
for seed in seed_list:
    res = least_squares(
        residual,
        seed,
        jac=jacobian,
        method="trf",
        max_nfev=12000,
        xtol=1e-14,
        ftol=1e-14,
        gtol=1e-14,
    )
    ss = float(np.sum(res.fun**2))
    if ss < 1e-22:
        raw_solutions.append(res.x)

# cluster near-duplicate roots
clusters = []
for x in raw_solutions:
    if not any(np.linalg.norm(x - y) < 1e-7 for y in clusters):
        clusters.append(x)

print(f"Found {len(clusters)} exact real branches (within numerical tolerance).")

def sigma_profile(x: np.ndarray, mu_vals: np.ndarray) -> np.ndarray:
    return 1 - x[3] * mu_vals**2 + x[4] * mu_vals**4

def bprof_profile(x: np.ndarray, mu_vals: np.ndarray) -> np.ndarray:
    return x[1] + x[2] * mu_vals**2

mus = np.linspace(-1.0, 1.0, 2001)
for i, x in enumerate(clusters):
    sig = sigma_profile(x, mus)
    bpr = bprof_profile(x, mus)
    print("")
    print(f"Branch {i}:")
    print(f"  A      = {x[0]: .15f}")
    print(f"  Btan0  = {x[1]: .15f}")
    print(f"  Btan2  = {x[2]: .15f}")
    print(f"  q      = {x[3]: .15f}")
    print(f"  r      = {x[4]: .15f}")
    print(f"  sigma min/max = {sig.min(): .12f}, {sig.max(): .12f}")
    print(f"  B_tan min/max = {bpr.min(): .12f}, {bpr.max(): .12f}")

# choose the Family-1-like positive-flare branch: q>0, r>0, sigma(mu)>0
physical = None
for x in clusters:
    sig = sigma_profile(x, mus)
    if x[3] > 0 and x[4] > 0 and sig.min() > 0:
        physical = x
        break

if physical is None:
    raise RuntimeError("No positive-flare Family-1-like branch found.")

print("")
print("Selected Family-1-like physical branch:")
for name, val in zip(["A", "Btan0", "Btan2", "q", "r"], physical):
    print(f"  {name} = {val: .15f}")

print("")
print("Residuals on selected branch:")
for key, val in zip(keys, residual(physical)):
    print(f"  {key}: {val:+.18e}")

beta0 = physical[1] + physical[2] / 3.0
beta2 = 2.0 * physical[2] / 3.0
print("")
print("Legendre form of the tangential wall profile on the selected branch:")
print("  B_tan(mu) = beta0 + beta2 P2(mu)")
print(f"  beta0 = {beta0: .15f}")
print(f"  beta2 = {beta2: .15f}")

# -----------------------------------------------------------------------------
# 4) Universal monopole prediction
# -----------------------------------------------------------------------------
section("4) Universal monopole prediction from the support match")

identity = sp.simplify(
    exprs[(0, 0)]
    - (
        exprs[(1, 1)]
        + sp.Rational(1, 2) * exprs[(1, 0)]
        - sp.Rational(1, 10) * exprs[(2, 0)]
        - sp.Rational(1, 5) * exprs[(2, 1)]
        - sp.Rational(1, 5) * exprs[(2, 2)]
    )
)
print("Exact structural identity:")
print("  K_(0,0) - [ K_(1,1) + 1/2 K_(1,0) - 1/10 K_(2,0) - 1/5 K_(2,1) - 1/5 K_(2,2) ]")
print(f"  = {identity}")

K00_from_support = sp.simplify(
    targets[(1, 1)]
    + sp.Rational(1, 2) * targets[(1, 0)]
    - sp.Rational(1, 10) * targets[(2, 0)]
    - sp.Rational(1, 5) * targets[(2, 1)]
    - sp.Rational(1, 5) * targets[(2, 2)]
)
monopole_target = sp.Rational(4, 45)
monopole_add = sp.simplify(monopole_target - K00_from_support)

print("")
print("So once the l=1,2 support targets are matched exactly, the model predicts")
print(f"  K_(0,0)^BL = {K00_from_support}")
print("independently of which exact support branch is chosen.")
print("")
print("To recover the carried-forward monopole target K_(0,0)=4/45,")
print(f"the separate monopole wall add-on must be {monopole_add}.")

k00_numeric = float(sp.N(K00_from_support))
print("")
print(f"Numerically: K_(0,0)^BL = {k00_numeric:.15f}")
print(f"            monopole add = {float(sp.N(monopole_add)):.15f}")

# -----------------------------------------------------------------------------
# 5) Conclusion
# -----------------------------------------------------------------------------
section("5) Conclusion")

print("Main result:")
print("- The smallest genuine PDE-side completion of the strict isotropic layer that")
print("  works is to keep the normal penetration moment isotropic but promote the")
print("  tangential wall moment to the axisymmetric profile B_tan(mu)=Btan0+Btan2 mu^2.")
print("- That single extra tangential wall-stress degree of freedom is enough to")
print("  reproduce the entire l=1,2 support sector exactly.")
print("- The monopole remains separate, but now with an exact universal prediction")
print("  for the raw boundary-layer value K_(0,0)^BL = -757/2520 and therefore an")
print("  exact required monopole add-on 109/280.")
print("")
print("This is the first strict boundary-layer / soft-wall PDE completion that")
print("actually reproduces the solved 2PN support sector, rather than fitting the")
print("support operator directly.")
