#!/usr/bin/env python3
"""
Projected Maxwell primitive bridge

This script starts one step closer to the Stage-4 / Stage-5 one-port moving-throat
bundle. Instead of inserting abstract bundle corrections z_n and n_n directly,
it perturbs the primitive Maxwell/mixed one-port data:

    Q, S2, H, Delta, P, Gw

by small mouth-local projected-Maxwell slippages:

    q1, s1, h1, d1, p1, g1

and derives the induced bundle corrections

    z0, z2, z4, n0, n2, n4

exactly at first order.

These are the quantities that then feed the full-bundle formulas for
D0, D2, D4, P0, P2, P4, Xi1, and the isotropic compatibility surface.
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


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    value = sp.factor(sp.together(sp.simplify(expr)))
    if value == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


def main() -> None:
    Q, S2, H, Delta, P, Gw = sp.symbols("Q S2 H Delta P Gw", nonzero=True)
    q1, s1, h1, d1, p1, g1 = sp.symbols("q1 s1 h1 d1 p1 g1")
    ell = sp.symbols("ell")

    # Stage-4 / Stage-5 primitive one-port formulas
    Z0 = Q / Delta
    Z2 = (Q * S2 - H * Delta) / Delta**2
    Z4 = (Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3

    N0 = P**2 / Delta**2
    N2 = 2 * P * (P * S2 - Delta * Gw) / Delta**3
    N4 = (Delta**2 * Gw**2 - 2 * Delta * P**2 - 4 * Delta * P * S2 * Gw + 3 * P**2 * S2**2) / Delta**4

    subs = {
        Q: Q + ell * q1,
        S2: S2 + ell * s1,
        H: H + ell * h1,
        Delta: Delta + ell * d1,
        P: P + ell * p1,
        Gw: Gw + ell * g1,
    }

    def dlin(expr: sp.Expr) -> sp.Expr:
        return sp.simplify((sp.series(sp.expand(expr.subs(subs)), ell, 0, 2).removeO() - expr) / ell)

    z0 = dlin(Z0)
    z2 = dlin(Z2)
    z4 = dlin(Z4)
    n0 = dlin(N0)
    n2 = dlin(N2)
    n4 = dlin(N4)

    assert_zero("z0 closed form", z0 - (Delta * q1 - Q * d1) / Delta**2)
    assert_zero(
        "z2 closed form",
        z2 - (-Delta**2 * h1 + Delta * (H * d1 + Q * s1 + S2 * q1) - 2 * Q * S2 * d1) / Delta**3,
    )
    assert_zero(
        "z4 closed form",
        z4
        - (
            -Delta**2 * H * s1
            - Delta**2 * S2 * h1
            - Delta**2 * q1
            + 2 * Delta * H * S2 * d1
            + 2 * Delta * Q * S2 * s1
            + 2 * Delta * Q * d1
            + Delta * S2**2 * q1
            - 3 * Q * S2**2 * d1
        )
        / Delta**4,
    )
    assert_zero("n0 closed form", n0 - 2 * P * (Delta * p1 - P * d1) / Delta**3)
    assert_zero(
        "n2 closed form",
        n2
        - (
            -(
                2 * Delta**2 * (Gw * p1 + P * g1)
                - 2 * Delta * P * (2 * Gw * d1 + P * s1 + 2 * S2 * p1)
                + 6 * P**2 * S2 * d1
            )
            / Delta**4
        ),
    )
    assert_zero(
        "n4 closed form",
        n4
        - 2
        * (
            Delta**3 * Gw * g1
            - Delta**2 * Gw**2 * d1
            - 2 * Delta**2 * Gw * P * s1
            - 2 * Delta**2 * Gw * S2 * p1
            - 2 * Delta**2 * P * S2 * g1
            - 2 * Delta**2 * P * p1
            + 6 * Delta * Gw * P * S2 * d1
            + 3 * Delta * P**2 * S2 * s1
            + 3 * Delta * P**2 * d1
            + 3 * Delta * P * S2**2 * p1
            - 6 * P**2 * S2**2 * d1
        )
        / Delta**5,
    )

    slips = {Q: q1, S2: s1, H: h1, Delta: d1, P: p1, Gw: g1}

    def frechet(expr: sp.Expr) -> sp.Expr:
        return sp.simplify(sum(sp.diff(expr, var) * slip for var, slip in slips.items()))

    assert_zero("z0 Frechet derivative", z0 - frechet(Z0))
    assert_zero("z2 Frechet derivative", z2 - frechet(Z2))
    assert_zero("z4 Frechet derivative", z4 - frechet(Z4))
    assert_zero("n0 Frechet derivative", n0 - frechet(N0))
    assert_zero("n2 Frechet derivative", n2 - frechet(N2))
    assert_zero("n4 Frechet derivative", n4 - frechet(N4))

    D0sym, N0sym = sp.symbols("D0 N0", nonzero=True)
    Xi1_static = sp.simplify(n0 / N0 + z0 / D0sym)
    assert_zero(
        "primitive Xi static formula",
        Xi1_static - (2*p1/P - 2*d1/Delta + (Delta*q1 - Q*d1)/(D0sym*Delta**2)),
    )

    K, B0, Z0slot, N0base, D0target, Ptarget, S, T = sp.symbols("K B0 Z0slot N0base D0target Ptarget S T", nonzero=True)
    K_one_p = sp.solve(
        sp.Eq((K - B0 - (Z0slot + ell * z0)) * (T + ell * z4), 3 * (S + ell * z2) ** 2),
        K,
    )[0]
    K_norm_p = sp.solve(
        sp.Eq((N0base + ell * n0) / (K - B0 - (Z0slot + ell * z0)), Ptarget),
        K,
    )[0]
    compat_direct_p = sp.simplify((N0base + ell * n0) / Ptarget - 3 * (S + ell * z2) ** 2 / (T + ell * z4))
    Ptarget_transport_p = sp.simplify((N0base + ell * n0) / D0target)
    K_norm_transport_p = sp.solve(
        sp.Eq((N0base + ell * n0) / (K - B0 - (Z0slot + ell * z0)), Ptarget_transport_p),
        K,
    )[0]
    compat_transport_p = sp.simplify(K_norm_transport_p - K_one_p)
    compat_direct = compat_direct_p.subs(ell, 0)
    compat_transport = compat_transport_p.subs(ell, 0)
    compat = K_norm_p.subs(ell, 0) - K_one_p.subs(ell, 0)
    dCompat = sp.simplify((sp.series(sp.expand(K_norm_p - K_one_p), ell, 0, 2).removeO() - compat) / ell)
    dCompat_direct = sp.simplify((sp.series(sp.expand(compat_direct_p), ell, 0, 2).removeO() - compat_direct) / ell)
    dCompat_transport = sp.simplify((sp.series(sp.expand(compat_transport_p), ell, 0, 2).removeO() - compat_transport) / ell)
    assert_zero(
        "K_one solve round-trip",
        (K_one_p - B0 - (Z0slot + ell * z0)) * (T + ell * z4) - 3 * (S + ell * z2) ** 2,
    )
    assert_zero(
        "K_norm solve round-trip",
        (N0base + ell * n0) / (K_norm_p - B0 - (Z0slot + ell * z0)) - Ptarget,
    )
    assert_zero(
        "primitive compatibility shift from competing K surfaces",
        dCompat_direct - (n0 / Ptarget - 6 * S * z2 / T + 3 * S**2 * z4 / T**2),
    )
    assert_zero("primitive transported-target normalization K surface", K_norm_transport_p - (B0 + Z0slot + ell * z0 + D0target))
    assert_zero(
        "primitive transported-target compatibility surface",
        compat_transport_p - (D0target - 3 * (S + ell * z2) ** 2 / (T + ell * z4)),
    )
    assert_zero(
        "primitive transported-target compatibility shift",
        dCompat_transport - (-6 * S * z2 / T + 3 * S**2 * z4 / T**2),
    )
    # The derived z0 = (Delta*q1 - Q*d1)/Delta**2 enters K_norm_p and K_one_p
    # through identical +ell*z0 terms, so it must cancel from their difference.
    # The substantive check is that the q1 and d1 partials of the compatibility
    # difference receive NO contribution from the z0 channel: they must equal
    # the q1 and d1 partials of compat_direct_p, which depends only on n0, z2, z4.
    for slip in (q1, d1):
        assert_zero(
            f"primitive fixed-target compatibility has no z0 channel in {slip}",
            sp.diff(K_norm_p - K_one_p, slip) - sp.diff(compat_direct_p, slip),
        )
        assert_zero(
            f"primitive transported-target compatibility has no z0 channel in {slip}",
            sp.diff(K_norm_transport_p - K_one_p, slip) - sp.diff(compat_transport_p, slip),
        )
    # Conversely, the normalization K surface alone DOES retain a q1/d1 channel
    # through the derived z0.
    assert_nonzero(
        "primitive normalization K surface retains q1 channel from z0",
        sp.diff(K_norm_p, q1),
    )
    assert_nonzero(
        "primitive normalization K surface retains d1 channel from z0",
        sp.diff(K_norm_p, d1),
    )
    assert_nonzero(
        "mutated primitive compatibility transport should fail",
        dCompat_direct - (n0 / Ptarget - 6 * S * z2 / T - 3 * S**2 * z4 / T**2),
    )
    assert_nonzero(
        "mutated primitive transported-target compatibility should fail",
        dCompat_transport - (-6 * S * z2 / T - 3 * S**2 * z4 / T**2),
    )

    lines = []
    lines.append(banner("1) Primitive one-port Maxwell/mixed formulas"))
    lines.append(line("Z0 =", Z0))
    lines.append(line("Z2 =", Z2))
    lines.append(line("Z4 =", Z4))
    lines.append(line("N0 =", N0))
    lines.append(line("N2 =", N2))
    lines.append(line("N4 =", N4))
    lines.append(
        "We now perturb the primitive one-port data by small mouth-local projected-Maxwell slippages:\n"
        "  Q -> Q + ell*q1,\n"
        "  S2 -> S2 + ell*s1,\n"
        "  H  -> H  + ell*h1,\n"
        "  Delta -> Delta + ell*d1,\n"
        "  P -> P + ell*p1,\n"
        "  Gw -> Gw + ell*g1.\n"
    )

    lines.append(banner("2) Induced bundle corrections z0, z2, z4, n0, n2, n4"))
    lines.append(line("z0 =", z0))
    lines.append(line("z2 =", z2))
    lines.append(line("z4 =", z4))
    lines.append(line("n0 =", n0))
    lines.append(line("n2 =", n2))
    lines.append(line("n4 =", n4))

    lines.append(banner("3) Reading the primitive slippages"))
    lines.append(
        "Three useful structural facts are immediate:\n"
        "  1. z0 depends only on (q1, d1).\n"
        "  2. z2 depends only on (q1, s1, h1, d1).\n"
        "  3. n0 depends only on (p1, d1).\n"
        "So the first isotropic normalization shift is controlled by the primitive transfer-leg slippage p1,\n"
        "the primitive denominator slippage d1, and the primitive conservative source slippage q1.\n"
    )

    lines.append(banner("4) Direct primitive contribution to Xi1 = P1/P0"))
    lines.append(line("Xi1^(primitive,static) =", Xi1_static))
    lines.append(
        "So at the one-port primitive level the weak-axisymmetric static prefactor slope receives contributions from:\n"
        "  - transfer-leg slippage p1,\n"
        "  - denominator slippage d1,\n"
        "  - and direct conservative source slippage q1 through z0/D0.\n"
    )

    lines.append(banner("5) Direct primitive contribution to the isotropic compatibility shift"))
    lines.append(line("delta Compat^(primitive) =", dCompat))
    lines.append(
        "Because delta Compat depends on z2, z4, and n0 only, the primitive contributors that matter here are\n"
        "  q1, s1, h1, d1, p1,\n"
        "while the Gw slope g1 first appears only through n2 and n4, i.e. in the constant-prefactor transport.\n"
        "The script also checks that z0 cancels from both the fixed-target and transported-target compatibility shifts,\n"
        "even though the normalization K surface still carries z0 before compatibility elimination.\n"
    )

    lines.append(banner("6) Practical reading"))
    lines.append(
        "At this one-port primitive level the projected-Maxwell near-throat correction has a clean hierarchy:\n"
        "  - (q1, d1, p1) are the first static normalization movers,\n"
        "  - (s1, h1) first matter through z2 and z4, hence through one-pole / compatibility transport,\n"
        "  - g1 first matters for N2 and N4, hence for whether the branch stays on the constant-prefactor outgoing surface.\n"
        "\n"
        "This is a useful way to decide which concrete throat-local projected ansatz to try next.\n"
    )

    print("".join(lines))
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
