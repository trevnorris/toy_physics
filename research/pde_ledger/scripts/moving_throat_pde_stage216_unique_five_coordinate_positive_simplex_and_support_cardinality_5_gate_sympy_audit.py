import sympy as sp
from pathlib import Path


def main() -> str:
    out = []
    out.append("Stage 216 SymPy audit — unique five-coordinate positive simplex and support-cardinality-5 gate")
    out.append("")

    # ------------------------------------------------------------------
    # Primitive data
    # ------------------------------------------------------------------
    k_lam, k_c, k_gam, k_U, k_W = sp.symbols("k_lambda k_c k_gamma k_U k_W", positive=True)
    ks = [k_lam, k_c, k_gam, k_U, k_W]
    k_norm_sq = sum(k**2 for k in ks)
    k_norm = sp.sqrt(k_norm_sq)

    out.append("1. Primitive five-coordinate simplex ledger")
    out.append("   primitive axes = (lambda, c, gamma, U, W)")
    out.append("   number of codimension-one faces = 5")
    out.append("   boundary faces = Q_hatlambda, Q_hatc, Q_hatgamma, Q_hatU, Q_hatW")
    out.append("")

    # ------------------------------------------------------------------
    # Gradient-optimal interior point
    # ------------------------------------------------------------------
    a_grad = [sp.simplify(k / k_norm) for k in ks]
    grad_norm_sq = sp.simplify(sum(a**2 for a in a_grad))
    k_grad = sp.simplify(sum(a * k for a, k in zip(a_grad, ks)))

    out.append("2. Gradient-optimal interior ray")
    out.append(f"   ||a_grad||^2 = {grad_norm_sq}")
    out.append(f"   k_grad = {k_grad}")

    face_labels = ["Q_hatlambda", "Q_hatc", "Q_hatgamma", "Q_hatU", "Q_hatW"]
    face_max_sq = [
        k_c**2 + k_gam**2 + k_U**2 + k_W**2,
        k_lam**2 + k_gam**2 + k_U**2 + k_W**2,
        k_lam**2 + k_c**2 + k_U**2 + k_W**2,
        k_lam**2 + k_c**2 + k_gam**2 + k_W**2,
        k_lam**2 + k_c**2 + k_gam**2 + k_U**2,
    ]
    diffs = [sp.simplify(k_norm_sq - s) for s in face_max_sq]
    for label, diff in zip(face_labels, diffs):
        out.append(f"   (k_5^grad)^2 - (k_{label}^grad)^2 = {diff}")
    out.append("")

    # Ratio patch a_lambda > 0
    r_grad = sp.simplify(k_c / k_lam)
    s_grad = sp.simplify(k_gam / k_lam)
    t_grad = sp.simplify(k_U / k_lam)
    u_grad = sp.simplify(k_W / k_lam)
    out.append("   interior ratio patch (a_lambda > 0):")
    out.append(f"   r_grad = {r_grad}")
    out.append(f"   s_grad = {s_grad}")
    out.append(f"   t_grad = {t_grad}")
    out.append(f"   u_grad = {u_grad}")
    out.append("")

    # ------------------------------------------------------------------
    # Equal-mix barycenter and off-diagonal leverage
    # ------------------------------------------------------------------
    a1, a2, a3, a4, a5 = sp.symbols("a1 a2 a3 a4 a5", real=True)
    asum = a1 + a2 + a3 + a4 + a5
    norm_sq = a1**2 + a2**2 + a3**2 + a4**2 + a5**2
    w_sigma = 2 * (
        a1*a2 + a1*a3 + a1*a4 + a1*a5 + a2*a3 + a2*a4 + a2*a5 + a3*a4 + a3*a5 + a4*a5
    )
    identity_check = sp.expand(w_sigma - (asum**2 - norm_sq))
    slack = sp.expand(5 * norm_sq - asum**2)
    pair_sum = sum((x - y)**2 for i, x in enumerate([a1, a2, a3, a4, a5]) for y in [a1, a2, a3, a4, a5][i+1:])
    slack_check = sp.expand(slack - pair_sum)

    out.append("3. Total ten-way off-diagonal leverage")
    out.append(f"   w_sigma - ((sum a_i)^2 - sum a_i^2) = {identity_check}")
    out.append(f"   5*sum a_i^2 - (sum a_i)^2 - sum_(p<q)(a_p-a_q)^2 = {slack_check}")

    a_eq = [sp.Integer(1) / sp.sqrt(5)] * 5
    w_eq = sp.simplify(w_sigma.subs({a1: a_eq[0], a2: a_eq[1], a3: a_eq[2], a4: a_eq[3], a5: a_eq[4]}))
    w_quad_eq = sp.Integer(3)  # imported from Stage 213
    w_triple_eq = sp.Integer(2)  # imported from Stage 210/194 chain
    w_pair_eq = sp.Integer(1)
    out.append(f"   equal-mix barycenter a_eq = (1,1,1,1,1)/sqrt(5)")
    out.append(f"   w_sigma(a_eq) = {w_eq}")
    out.append(f"   quadruple equal-mix leverage = {w_quad_eq}")
    out.append(f"   triple equal-mix leverage = {w_triple_eq}")
    out.append(f"   pairwise equal-mix leverage = {w_pair_eq}")
    out.append("")

    # ------------------------------------------------------------------
    # Fixed-point certified bracket
    # ------------------------------------------------------------------
    H0, kappa = sp.symbols("H0 kappa", positive=True)
    k = sp.symbols("k", positive=True)
    tau = sp.simplify(2 * H0 / (k + sp.sqrt(k**2 - 2 * H0 * kappa)))
    out.append("4. Fixed-point certified bracket")
    out.append(f"   tau(H0,k,kappa) = {tau}")
    out.append("   with k = k_5(a) and kappa = a^T H_(5,star) a")
    out.append("")

    # ------------------------------------------------------------------
    # Canonical screen packet and support-cardinality ceiling
    # ------------------------------------------------------------------
    out.append("5. Canonical five-way screen packet")
    out.append("   imported faces  = 5")
    out.append("   interior screens = {gradient-optimal, equal-mix}")
    out.append("   total canonical rows = 7")
    out.append("")

    out.append("6. Support-cardinality ceiling")
    out.append("   free-quintuple dimension = 5")
    out.append("   therefore no local mixed ray can have support-cardinality > 5")
    out.append("")

    return "\n".join(out)


if __name__ == "__main__":
    text = main()
    print(text)
