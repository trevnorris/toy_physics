
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr: sp.Expr) -> None:
    s = sp.simplify(sp.expand(expr))
    print(f"{name} = {s}")
    if s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 116 — FINITE D/N MIXED-TUBE REALIZATION")

a, L_W = sp.symbols("a L_W", positive=True, real=True)
K_s, K_q, lam = sp.symbols("K_s K_q lam", positive=True, real=True)
z = sp.symbols("z", real=True)

# D/N half-wave eigenvalue derivation: solve q'' + k^2 q = 0 with q(0)=0, q'(L_W)=0.
# q(x) = sin(k*x) satisfies q(0)=0; the second BC requires cos(k*L_W) = 0, so the
# smallest positive eigenvalue is k_W = pi/(2*L_W).
x_var = sp.symbols("x", real=True)
k_sym = sp.symbols("k_sym", positive=True, real=True)
q_trial = sp.sin(k_sym * x_var)
ode_residual = sp.simplify(sp.diff(q_trial, x_var, 2) + k_sym**2 * q_trial)
expect_zero("D/N trial satisfies q'' + k^2 q = 0", ode_residual)
bc_left = sp.simplify(q_trial.subs(x_var, 0))
expect_zero("D/N trial satisfies q(0) = 0", bc_left)
bc_right = sp.simplify(sp.diff(q_trial, x_var).subs(x_var, L_W))
k_W_value = sp.pi / (2 * L_W)
expect_zero(
    "D/N trial satisfies q'(L_W) = 0 at k = pi/(2 L_W)",
    sp.simplify(bc_right.subs(k_sym, k_W_value)),
)

omega_sym, c_s_sym = sp.symbols("omega c_s", positive=True, real=True)
Omega_W = k_W_value * c_s_sym
kappa0_derived = sp.simplify((omega_sym / Omega_W)**2 / (a * omega_sym / c_s_sym)**2)
expect_zero(
    "kappa0 from D/N eigenvalue matches geometric expression",
    kappa0_derived - 4 * L_W**2 / (sp.pi**2 * a**2),
)

r_c = lam**2 / (K_s * K_q)

kappa0_from_tube = sp.simplify(kappa0_derived)
L_W_required = sp.solve(
    sp.Eq(kappa0_from_tube, (1 + r_c) / 3),
    L_W,
)[0]

print("kappa0 from D/N half-wave tube =", kappa0_from_tube)
print("Required tube length L_W =", sp.simplify(L_W_required))
expect_zero(
    "tube-length law: L_W = pi a sqrt((1+r_c)/3) / 2",
    sp.simplify(L_W_required - sp.pi * a * sp.sqrt((1 + r_c) / 3) / 2),
)

# --- Renormalization to canonical coefficients (REPORTED, not asserted) ---
# The load-bearing, falsifiable physics is verified above (eigenvalue, kappa0
# collapse, tube-length law). The renormalization below carries no independent
# falsifiable content in this stage: kappa0 is already established above and
# gamma0 is an upstream-carried input (Stage 115 compensation requirement) with
# no in-stage derivation. Dividing out the common (1+r_c) hybridization factor
# is therefore a definitional consequence, so these values are PRINTED, not
# asserted (an expect_zero here would be tautological).
kappa0_bare = sp.simplify(kappa0_from_tube.subs(L_W, L_W_required))  # derived tube coeff at required length
gamma0_bare = sp.simplify((1 + r_c) / 9)                            # upstream-carried input (Stage 115), not derived in-stage
common_scale = 1 + r_c
kappa_c = sp.simplify(kappa0_bare / common_scale)
gamma_c = sp.simplify(gamma0_bare / common_scale)
print("Renormalization (definitional consequence, not an independent check):")
print("  kappa0_bare (derived tube coeff at L_W_required) =", kappa0_bare)
print("  gamma0_bare (upstream-carried input, Stage 115)   =", gamma0_bare)
print("  kappa_c = kappa0_bare/(1+r_c) =", kappa_c)
print("  gamma_c = gamma0_bare/(1+r_c) =", gamma_c)

print("\nSummary:")
print("  D/N half-wave mixed tube length:", sp.simplify(L_W_required))
print("  Bare mixed coefficients: kappa0 =", kappa0_bare, ", gamma0 =", gamma0_bare)
print("  Final coefficients: kappa_c =", kappa_c, ", gamma_c =", gamma_c)
