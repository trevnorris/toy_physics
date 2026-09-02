"""Independent (engine-free) term-level derivation of the maximum
background-jet order carried by the bulk momentum EL row

    row_a = dL/du_a - d_i (dL / d(d_i u_a))

of a quadratic first-order brane energy density whose only background
spurions are first jets of W(x).  Prints operands and residuals; states
no conclusion.

Terms:
  A  L = gamma * (grad W · grad theta) * div u     [DIVU x GRADTHETA x spurion]
  B  L = (kappa/2) * W * |grad u|^2                [uniform coeff of W, no spurion]
  C  L = gamma * (d_x W) * d_x(W0 e/W) * d_x u_1   [spurion x GRADEW x GRADU]
  D  L = gamma * (grad W · u) * theta              [undifferentiated-u spurion]

For A, also form the CONSTRAINT-MIXED object  EL_u - grad(EL_theta),
which is what a virtual-work elimination of theta through linearized
continuity  theta ~ -div u  injects into a U-row.
"""
from __future__ import annotations

import sympy as sp

x, y, z = sp.symbols("x y z")
COORDS = (x, y, z)
W = sp.Function("W")
theta = sp.Function("theta")
eW = sp.Function("eW")
u = [sp.Function("u1"), sp.Function("u2"), sp.Function("u3")]
gamma, kappa, W0 = sp.symbols("gamma kappa W0", nonzero=True)
POS = COORDS


def Wpos():
    return W(*POS)


def upos(a):
    return u[a](*POS)


def thpos():
    return theta(*POS)


def epos():
    return eW(*POS)


def grad(fexpr):
    return [sp.diff(fexpr, c) for c in COORDS]


def div_u():
    return sum(sp.diff(upos(i), COORDS[i]) for i in range(3))


def w_orders(expr):
    """Spatial-derivative orders of W appearing in expr (0 = undifferentiated W)."""
    found = set()
    expr = sp.expand(expr)
    for node in sp.preorder_traversal(expr):
        if isinstance(node, sp.Derivative) and node.expr == Wpos():
            found.add(sum(n for _v, n in node.variable_count))
        elif node == Wpos():
            found.add(0)
    return sorted(found)


def flux_u(L, a):
    return [sp.diff(L, sp.diff(upos(a), COORDS[i])) for i in range(3)]


def local_u(L, a):
    return sp.diff(L, upos(a))


def el_u(L, a):
    f = flux_u(L, a)
    return sp.expand(local_u(L, a) - sum(sp.diff(f[i], COORDS[i]) for i in range(3)))


def el_theta(L):
    flux = [sp.diff(L, sp.diff(thpos(), COORDS[i])) for i in range(3)]
    return sp.expand(sp.diff(L, thpos()) - sum(sp.diff(flux[i], COORDS[i]) for i in range(3)))


def report(name, L):
    print(f"\n=== {name} ===")
    p_flux = set()
    p_local = set()
    p_row = set()
    flux_orders_per_comp = []
    for a in range(3):
        f = flux_u(L, a)
        loc = local_u(L, a)
        row = el_u(L, a)
        fo = sorted({o for comp in f for o in w_orders(comp)})
        lo = w_orders(loc)
        ro = w_orders(row)
        flux_orders_per_comp.append(fo)
        p_flux.update(fo)
        p_local.update(lo)
        p_row.update(ro)
    dens_o = w_orders(L)
    print(f"{name}_DENSITY_W_ORDERS: {dens_o}")
    print(f"{name}_FLUX_dL_d(du)_W_ORDERS: {sorted(p_flux)}")
    print(f"{name}_LOCAL_dL_du_W_ORDERS: {sorted(p_local)}")
    print(f"{name}_UNCONSTRAINED_EL_U_W_ORDERS: {sorted(p_row)}")
    print(f"{name}_P_MAX_FLUX_ORDER: {max(p_flux) if p_flux else 0}")
    print(f"{name}_ROW_MAX_ORDER: {max(p_row) if p_row else 0}")
    return {
        "density": dens_o,
        "flux": sorted(p_flux),
        "local": sorted(p_local),
        "row": sorted(p_row),
        "L": L,
    }


# --- KEY FACT: a single jet of any order is sigma_W^1 ---
print("=== KEY_FACT_SIGMA_W_GRADE ===")
eta, sigma, L_W, W0s = sp.symbols("eta_bg sigma_W L_W W0", nonzero=True)
w1 = sp.Function("w1")
# W = W0 (1 + eta * w1) is the background; the engine's profile map is
# d^n W -> sigma_W * (profile jet) / L_W**(n-1)  for n>=1.
# Verify: W = W0 + (W-W0), sigma_W parametrizes the amplitude of EVERY spatial
# derivative of the profile, not a new power per extra derivative.
W_profile = W0s * (1 + eta * w1(*POS))
for n, multi in [(1, (x,)), (2, (x, x)), (3, (x, x, x)), (3, (x, y, z))]:
    dW = sp.diff(W_profile, *multi)
    # dW = W0 * eta * d^n w1 ; one power of the shape amplitude eta.
    amp = sp.Poly(sp.expand(dW), eta).degree()
    print(f"KEY_FACT_dW_order_{n}_{multi}_ETA_DEGREE: {amp}")
    print(f"KEY_FACT_dW_order_{n}_{multi}_EXPR: {dW}")
# product of two first jets:
prod = sp.diff(W_profile, x) * sp.diff(W_profile, y)
print(f"KEY_FACT_PRODUCT_TWO_FIRST_JETS_ETA_DEGREE: {sp.Poly(sp.expand(prod), eta).degree()}")


# Term A
gW = grad(Wpos())
gth = grad(thpos())
L_A = gamma * sum(gW[k] * gth[k] for k in range(3)) * div_u()
resA = report("TERM_A_DIVU_GRADTHETA_SPURION", L_A)

# Term B
L_B = (kappa / 2) * Wpos() * sum(
    sp.diff(upos(a), COORDS[i]) ** 2 for a in range(3) for i in range(3)
)
resB = report("TERM_B_UNIFORM_W_TIMES_GRADU_SQ", L_B)

# Term C: spurion_x * d_x(W0 e/W) * d_x u1
local_e = W0 * epos() / Wpos()
L_C = gamma * sp.diff(Wpos(), x) * sp.diff(local_e, x) * sp.diff(upos(0), x)
resC = report("TERM_C_SPURION_GRADEW_GRADU", L_C)
# After expand, TERM_C has (dW)^2 * e/W^2 pieces (two-jet, sigma_W^2) and
# (dW)(de) pieces (single jet).  Split them.
L_C_exp = sp.expand(L_C)
two_jet = sp.Integer(0)
one_jet = sp.Integer(0)
for term in sp.Add.make_args(L_C_exp):
    orders = w_orders(term)
    # a product of two first jets has two Derivative(W,...) factors of order 1,
    # which w_orders reports as {1} still.  Count Derivative factors.
    n_w_deriv = sum(
        1
        for node in sp.preorder_traversal(term)
        if isinstance(node, sp.Derivative) and node.expr == Wpos()
    )
    n_w_undiff = sum(1 for node in sp.preorder_traversal(term) if node == Wpos())
    # Use degree in a dummy replacement: W -> 1 + eps * W, count eps in the
    # derivatives-only sense.  Simpler: substitute a scaled W.
print("TERM_C_DENSITY_EXPANDED_W_ORDERS: ", w_orders(L_C_exp))

# Term D
L_D = gamma * sum(gW[k] * upos(k) for k in range(3)) * thpos()
resD = report("TERM_D_U_THETA_SPURION", L_D)


# Constraint mixing on A: EL_u - grad(EL_theta)
print("\n=== CONSTRAINT_MIX_A_EL_U_MINUS_GRAD_MU_THETA ===")
mu = el_theta(L_A)
print(f"MU_THETA_A_W_ORDERS: {w_orders(mu)}")
mixed_orders = set()
unconst_orders = set()
for a in range(3):
    un = el_u(L_A, a)
    mixed = sp.expand(un - sp.diff(mu, COORDS[a]))
    unconst_orders.update(w_orders(un))
    mixed_orders.update(w_orders(mixed))
    print(f"MIXED_ROW_{a}_W_ORDERS: {w_orders(mixed)}")
    print(f"UNCONSTRAINED_ROW_{a}_W_ORDERS: {w_orders(un)}")
print(f"CONSTRAINT_MIX_MAX_W_ORDER: {max(mixed_orders) if mixed_orders else 0}")
print(f"UNCONSTRAINED_MAX_W_ORDER: {max(unconst_orders) if unconst_orders else 0}")

# Residual: mixed - unconstrained = -grad mu_theta.  Report its orders.
grad_mu_orders = set()
for a in range(3):
    grad_mu_orders.update(w_orders(sp.diff(mu, COORDS[a])))
print(f"NEG_GRAD_MU_THETA_W_ORDERS: {sorted(grad_mu_orders)}")

# Continuity-style extra lower term mu * (grad W)/W  (LAB_HELD, RHO4 ~ W)
print("\n=== CONSTRAINT_MIX_A_PLUS_MU_GRADLOGW ===")
extra_orders = set()
for a in range(3):
    extra = sp.expand(el_u(L_A, a) - sp.diff(mu, COORDS[a]) - mu * sp.diff(Wpos(), COORDS[a]) / Wpos())
    extra_orders.update(w_orders(extra))
    print(f"FULL_CONSTRAINED_ROW_{a}_W_ORDERS: {w_orders(extra)}")
print(f"FULL_CONSTRAINED_MAX_W_ORDER: {max(extra_orders) if extra_orders else 0}")

print("\nDONE")
