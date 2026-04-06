import numpy as np
from scipy.integrate import quad

G = 1.0
M = 1.0
c0 = 1.0
mu = G * M

CASES = (
    (100.0, 10000.0, 10000.0),
    (200.0, 10000.0, 10000.0),
    (100.0, 20000.0, 20000.0),
)


def N_of_r(r):
    return 1.0 + 2.0 * mu / (c0**2 * r)


def numerical_delay(b, z_e, z_r):
    integrand = lambda z: N_of_r(np.hypot(b, z)) / c0
    t_num, _ = quad(integrand, -z_e, z_r, epsabs=1e-12, epsrel=1e-12, limit=200)
    t_flat = (z_e + z_r) / c0
    return t_num - t_flat


def exact_delay(b, z_e, z_r):
    return 2.0 * mu / c0**3 * (np.arcsinh(z_e / b) + np.arcsinh(z_r / b))


def asymptotic_delay(b, z_e, z_r):
    return 2.0 * mu / c0**3 * np.log(4.0 * z_e * z_r / b**2)


print("Shapiro-delay check for the preferred n=5 pure-refraction branch")
print("Using N(r) ≈ 1 + 2 GM/(r c^2) and the direct line integral from the paper.")
print(
    f"{'b':>8} | {'Z_E':>8} | {'Z_R':>8} | "
    f"{'num':>12} | {'exact':>12} | {'asymptotic':>12} | {'num/asym':>10}"
)
print("-" * 92)

for b, z_e, z_r in CASES:
    delay_num = numerical_delay(b, z_e, z_r)
    delay_exact = exact_delay(b, z_e, z_r)
    delay_asym = asymptotic_delay(b, z_e, z_r)
    print(
        f"{b:8.1f} | {z_e:8.1f} | {z_r:8.1f} | "
        f"{delay_num:12.6f} | {delay_exact:12.6f} | {delay_asym:12.6f} | "
        f"{delay_num / delay_asym:10.6f}"
    )

print("-" * 92)
print("Conclusion: the direct integral matches the exact one-way delay and stays close")
print("to the weak-field asymptotic formula 2GM c^-3 ln(4 r_E r_R / b^2).")
