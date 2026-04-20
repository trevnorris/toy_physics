#!/usr/bin/env python3
"""
moving_throat_pde_stage001_geometry_lift_sympy_audit.py

Stage 001 SymPy audit for the geometry lift and the first linearized wall PDE
bookkeeping.

Scope
-----
This script verifies the explicitly reusable algebra in Stage 001:

  • normalized real-harmonic bookkeeping for the monopole and grouped real P2,
  • the mouth-average extraction rule q_00 = 2 sqrt(pi) delta a,
  • the confinement chain-rule sign in the promoted level-set coupling,
  • the modal wall equation in the densitized convention used by the ledger,
  • the corresponding weighted-surface form before densitization,
  • the sourced modal-wall RHS forcing,
  • and a representative localized-Maxwell component variation with
    gauge-fixing and source current.

The stage is foundational and intentionally does not attempt to solve the full
coupled GNLS + Maxwell + moving-wall PDE.
"""

from __future__ import annotations

import sympy as sp
from sympy.calculus.euler import euler_equations


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        simplified = expr.applyfunc(lambda z: sp.simplify(sp.trigsimp(sp.expand(z))))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
        return

    simplified = sp.simplify(sp.trigsimp(sp.expand(expr)))
    print(f"{name} = {simplified}")
    if simplified != 0:
        raise AssertionError(f"{name} is not zero")


def harmonic_bookkeeping_audit() -> None:
    banner("SECTION I — REAL-HARMONIC BOOKKEEPING")

    theta, phi = sp.symbols("theta phi", real=True)
    dOmega = sp.sin(theta)

    y00 = sp.Rational(1, 2) / sp.sqrt(sp.pi)
    y20 = sp.sqrt(5) / (4 * sp.sqrt(sp.pi)) * (3 * sp.cos(theta) ** 2 - 1)
    y21c = sp.sqrt(15) / (2 * sp.sqrt(sp.pi)) * sp.sin(theta) * sp.cos(theta) * sp.cos(phi)
    y21s = sp.sqrt(15) / (2 * sp.sqrt(sp.pi)) * sp.sin(theta) * sp.cos(theta) * sp.sin(phi)
    y22c = sp.sqrt(15) / (4 * sp.sqrt(sp.pi)) * sp.sin(theta) ** 2 * sp.cos(2 * phi)
    y22s = sp.sqrt(15) / (4 * sp.sqrt(sp.pi)) * sp.sin(theta) ** 2 * sp.sin(2 * phi)

    basis = {
        "Y00": y00,
        "Y20": y20,
        "Y21c": y21c,
        "Y21s": y21s,
        "Y22c": y22c,
        "Y22s": y22s,
    }

    subbanner("I.1 — Norms and zero averages")
    for name, y in basis.items():
        avg = sp.simplify(sp.integrate(sp.integrate(y * dOmega, (phi, 0, 2 * sp.pi)), (theta, 0, sp.pi)))
        norm = sp.simplify(sp.integrate(sp.integrate(y * y * dOmega, (phi, 0, 2 * sp.pi)), (theta, 0, sp.pi)))
        print(f"{name}: average = {avg}, norm = {norm}")

    for name in ["Y20", "Y21c", "Y21s", "Y22c", "Y22s"]:
        y = basis[name]
        expect_zero(
            f"average({name})",
            sp.integrate(sp.integrate(y * dOmega, (phi, 0, 2 * sp.pi)), (theta, 0, sp.pi)),
        )

    expect_zero(
        "norm(Y00)-1",
        sp.integrate(sp.integrate(y00 * y00 * dOmega, (phi, 0, 2 * sp.pi)), (theta, 0, sp.pi)) - 1,
    )
    expect_zero(
        "norm(Y20)-1",
        sp.integrate(sp.integrate(y20 * y20 * dOmega, (phi, 0, 2 * sp.pi)), (theta, 0, sp.pi)) - 1,
    )
    expect_zero(
        "cross(Y00,Y20)",
        sp.integrate(sp.integrate(y00 * y20 * dOmega, (phi, 0, 2 * sp.pi)), (theta, 0, sp.pi)),
    )

    subbanner("I.2 — Mouth-average extraction")
    q00, q20, q21c, q22c = sp.symbols("q00 q20 q21c q22c", real=True)
    eta = q00 * y00 + q20 * y20 + q21c * y21c + q22c * y22c
    mouth_avg = sp.simplify(
        sp.integrate(sp.integrate(eta * dOmega, (phi, 0, 2 * sp.pi)), (theta, 0, sp.pi)) / (4 * sp.pi)
    )
    expect_zero("mouth average - q00/(2 sqrt(pi))", mouth_avg - q00 / (2 * sp.sqrt(sp.pi)))
    print("Therefore the physical mouth-average shift delta a satisfies q00 = 2 sqrt(pi) delta a.")

    subbanner("I.3 — Spherical Laplacian eigenvalues")

    def lap_s2(expr: sp.Expr) -> sp.Expr:
        return sp.simplify(
            (1 / sp.sin(theta)) * sp.diff(sp.sin(theta) * sp.diff(expr, theta), theta)
            + sp.diff(expr, phi, 2) / sp.sin(theta) ** 2
        )

    expect_zero("-Delta_S2 Y00", -lap_s2(y00))
    for name in ["Y20", "Y21c", "Y21s", "Y22c", "Y22s"]:
        expect_zero(f"-Delta_S2 {name} - 6 {name}", -lap_s2(basis[name]) - 6 * basis[name])


def confinement_chain_rule_audit() -> None:
    banner("SECTION II — LEVEL-SET CONFINEMENT LINEARIZATION")

    sigma0, eta, ell_c, eps = sp.symbols("Sigma0 eta ell_c eps", nonzero=True)
    vwall = sp.Function("Vwall")

    expr = vwall((sigma0 - eps * eta) / ell_c)
    first_var = sp.diff(expr, eps).subs(eps, 0)
    target = -eta * sp.diff(vwall(sigma0 / ell_c), sigma0)
    print("First variation d/dε Vwall((Sigma0 - ε eta)/ell_c)|_{ε=0} =")
    sp.pprint(first_var)
    expect_zero("linearized confinement variation", first_var - target)
    print("Equivalent note form: delta V_conf = -(V'_wall(Sigma0/ell_c) / ell_c) eta")


def modal_wall_action_audit() -> None:
    banner("SECTION III — MODAL WALL ACTION")

    t, w = sp.symbols("t w", real=True)
    ell = sp.symbols("ell", integer=True, nonnegative=True)
    mu_eta = sp.Function("mu_eta")(w)
    t_w = sp.Function("T_w")(w)
    t_omega = sp.Function("T_Omega")(w)
    k_eta = sp.Function("K_eta")(w)
    q = sp.Function("q")
    g = sp.Function("g")(w)

    k_l = k_eta + ell * (ell + 1) * t_omega

    subbanner("III.1 — Densitized convention used by the ledger")
    ldens = (
        sp.Rational(1, 2) * mu_eta * sp.diff(q(t, w), t) ** 2
        - sp.Rational(1, 2) * t_w * sp.diff(q(t, w), w) ** 2
        - sp.Rational(1, 2) * k_l * q(t, w) ** 2
    )
    el_dens = euler_equations(ldens, q(t, w), [t, w])[0]
    target_dens = -mu_eta * sp.diff(q(t, w), t, 2) + sp.diff(t_w * sp.diff(q(t, w), w), w) - k_l * q(t, w)
    expect_zero("densitized Euler-Lagrange equation", el_dens.lhs - target_dens)

    subbanner("III.2 — Weighted surface form before densitization")
    lweighted = g * ldens
    el_weighted = euler_equations(lweighted, q(t, w), [t, w])[0]
    target_weighted = (
        -g * mu_eta * sp.diff(q(t, w), t, 2)
        + sp.diff(g * t_w * sp.diff(q(t, w), w), w)
        - g * k_l * q(t, w)
    )
    expect_zero("weighted Euler-Lagrange equation", el_weighted.lhs - target_weighted)
    print("Dividing by g(w) gives the extra first-derivative term that disappears in the densitized convention.")

    subbanner("III.3 — Grouped real P2 restoring shift")
    expect_zero("K_l at ell=0", sp.simplify(k_l.subs(ell, 0) - k_eta))
    expect_zero("K_l at ell=2", sp.simplify(k_l.subs(ell, 2) - (k_eta + 6 * t_omega)))
    print("The constant 6 is not an unexplained literal; it is the specialization ell(ell+1) at ell = 2.")

    subbanner("III.4 — Explicit source forcing on the modal wall equation")
    S_lm = sp.Function("S_lm")
    f_ext = sp.Function("f_ext")
    source_total = S_lm(t, w) + f_ext(t, w)
    ldens_forced = ldens - q(t, w) * source_total
    el_forced = euler_equations(ldens_forced, q(t, w), [t, w])[0]
    target_forced = target_dens - source_total
    expect_zero("sourced densitized Euler-Lagrange equation", el_forced.lhs - target_forced)


def linearized_maxwell_audit() -> None:
    banner("SECTION IV — REPRESENTATIVE LOCALIZED-MAXWELL LINEARIZATION")

    x, w = sp.symbols("x w", real=True)
    gauge_xi, mu0 = sp.symbols("xi mu0", positive=True, real=True)
    Zloc = sp.Function("Z")
    Ax = sp.Function("A_x")
    Aw = sp.Function("A_w")
    Jx = sp.Function("J_x")
    Jw = sp.Function("J_w")

    Fwx = sp.diff(Ax(x, w), w) - sp.diff(Aw(x, w), x)
    divA = sp.diff(Ax(x, w), x) + sp.diff(Aw(x, w), w)

    lmax = (
        sp.Rational(1, 2) * Zloc(w) * Fwx**2
        - sp.Rational(1, 2) * divA**2 / gauge_xi
        + mu0 * (Jx(x, w) * Ax(x, w) + Jw(x, w) * Aw(x, w))
    )
    el_Ax = (
        sp.diff(sp.diff(lmax, sp.diff(Ax(x, w), x)), x)
        + sp.diff(sp.diff(lmax, sp.diff(Ax(x, w), w)), w)
        - sp.diff(lmax, Ax(x, w))
    )
    el_Aw = (
        sp.diff(sp.diff(lmax, sp.diff(Aw(x, w), x)), x)
        + sp.diff(sp.diff(lmax, sp.diff(Aw(x, w), w)), w)
        - sp.diff(lmax, Aw(x, w))
    )

    target_Ax = sp.diff(Zloc(w) * Fwx, w) - sp.diff(divA, x) / gauge_xi - mu0 * Jx(x, w)
    target_Aw = -sp.diff(Zloc(w) * Fwx, x) - sp.diff(divA, w) / gauge_xi - mu0 * Jw(x, w)

    expect_zero("localized-Maxwell x-component", el_Ax - target_Ax)
    expect_zero("localized-Maxwell w-component", el_Aw - target_Aw)


def main() -> None:
    harmonic_bookkeeping_audit()
    confinement_chain_rule_audit()
    modal_wall_action_audit()
    linearized_maxwell_audit()

    banner("STAGE 001 SUMMARY")
    print("Verified with SymPy:")
    print("  • normalized real-harmonic bookkeeping for the monopole and grouped real P2 basis;")
    print("  • the mouth-average extraction rule q00 = 2 sqrt(pi) delta a;")
    print("  • the confinement-chain-rule sign delta V_conf = -(V'_wall/ell_c) eta;")
    print("  • the Stage 001 modal wall equation in the densitized convention used by the ledger;")
    print("  • the sourced modal-wall RHS forcing S_lm^(psi,A) + f_lm^ext;")
    print("  • the weighted-surface variant that explains the previously noted surface-measure caveat;")
    print("  • and representative localized-Maxwell component equations with gauge-fixing and source current.")


if __name__ == "__main__":
    main()
