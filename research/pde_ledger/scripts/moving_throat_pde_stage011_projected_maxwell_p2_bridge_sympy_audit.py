#!/usr/bin/env python3
"""
Projected Maxwell -> P2 compatibility bridge.

This Stage 011 audit checks only the paper-card output:

1. the one-pole and fixed-target normalization K-surfaces,
2. the K-eliminated compatibility surface C = N0/Ptarget - 3 S^2/T,
3. its first variation,
4. the z0-cancellation in that first variation.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> str:
    return "\n" + "=" * 96 + f"\n{title}\n" + "=" * 96 + "\n"


def line(label: str, expr: sp.Expr | str) -> str:
    return f"{label}\n  {expr}\n"


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def main() -> None:
    eps = sp.symbols("eps")
    K, B0, Z0slot, N0, Ptarget, S, T = sp.symbols("K B0 Z0slot N0 Ptarget S T", nonzero=True)
    z0, z2, z4, n0 = sp.symbols("z0 z2 z4 n0")

    def lin(expr: sp.Expr) -> sp.Expr:
        return sp.expand(sp.series(sp.expand(expr), eps, 0, 2).removeO())

    K_one_pole_p = sp.solve(
        sp.Eq((K - B0 - (Z0slot + eps * z0)) * (T + eps * z4), 3 * (S + eps * z2) ** 2),
        K,
    )[0]
    K_norm_p = sp.solve(
        sp.Eq((N0 + eps * n0) / (K - B0 - (Z0slot + eps * z0)), Ptarget),
        K,
    )[0]

    K_one_pole = K_one_pole_p.subs(eps, 0)
    K_norm = K_norm_p.subs(eps, 0)
    dK_one_pole = sp.simplify((lin(K_one_pole_p) - K_one_pole) / eps)
    dK_norm = sp.simplify((lin(K_norm_p) - K_norm) / eps)

    compat_direct_p = sp.simplify((N0 + eps * n0) / Ptarget - 3 * (S + eps * z2) ** 2 / (T + eps * z4))
    compat_direct = compat_direct_p.subs(eps, 0)
    d_compat_direct = sp.simplify((lin(compat_direct_p) - compat_direct) / eps)
    d_compat = sp.simplify((lin(K_norm_p - K_one_pole_p) - (K_norm - K_one_pole)) / eps)

    assert_zero("compatibility surface after eliminating K", (K_norm_p - K_one_pole_p) - compat_direct_p)
    assert_zero("one-pole K shift", dK_one_pole - (z0 + 6 * S * z2 / T - 3 * S**2 * z4 / T**2))
    assert_zero("normalization K shift", dK_norm - (z0 + n0 / Ptarget))
    assert_zero("compatibility shift from competing K surfaces", d_compat - (dK_norm - dK_one_pole))
    assert_zero("compatibility first variation from eliminated surface", d_compat - d_compat_direct)
    assert_zero("compatibility first variation", d_compat_direct - (n0 / Ptarget - 6 * S * z2 / T + 3 * S**2 * z4 / T**2))
    assert_zero("fixed-target compatibility z0 cancellation", sp.diff(d_compat_direct, z0))

    lines = []
    lines.append(banner("1) K-eliminated isotropic compatibility surface"))
    lines.append("One-pole and fixed-target normalization surfaces:\n")
    lines.append(line("K_pole(eps) =", K_one_pole_p))
    lines.append(line("K_norm(eps) =", K_norm_p))
    lines.append(line("Compat =", compat_direct))
    lines.append(line("delta Compat =", d_compat_direct))
    lines.append(
        "Structural result:\n"
        "  z0 drops out of delta Compat after exact elimination of K and D0.\n"
        "  The static conservative shift changes P0 directly, but not this fixed-target compatibility surface.\n"
    )

    print("".join(lines))
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
