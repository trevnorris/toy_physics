#!/usr/bin/env python3
"""
4pn_local_target_import_audit.py

Fixed-chart local 4PN target import audit.

What this script does
---------------------
1. Derives the strict one-body 4PN Hamiltonian gate by applying the exact quartic
   Legendre compiler to the previously frozen one-body Schwarzschild Lagrangian.
2. Imports the reduced center-of-mass *local* 4PN ADM Hamiltonian target from the
   fixed-chart source used for the local solve.
3. Organizes the imported target into the natural local 4PN COM Hamiltonian slots:
      - free kinetic: 1 slot,
      - G/r block   : 5 slots,
      - G^2/r^2     : 4 slots,
      - G^3/r^3     : 3 slots,
      - G^4/r^4     : 2 slots,
      - G^5/r^5     : 1 slot.
4. Verifies that the imported target reduces exactly to the one-body Hamiltonian gate
   in the strict test-mass limit nu -> 0.
5. Extracts the highest power of nu appearing in each local interaction block.
6. Answers the first exact target-import theorem gate:

      Do the upper local blocks G^4/r^4 or G^5/r^5 contain nu^3 or nu^4 tails?

The answer is no for the imported fixed-chart local ADM target.

Primary fixed-chart source
--------------------------
P. Jaranowski and G. Schaefer,
"Derivation of local-in-time fourth post-Newtonian ADM Hamiltonian for spinless compact binaries",
arXiv:1508.01016, Eq. (8.41e).

Scope note
----------
This script intentionally imports only the *local-in-time* part of the 4PN target.
The nonlocal/tail sector remains separate and is not folded into this audit.
"""

from __future__ import annotations

import sympy as sp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


def max_nu_degree(expr: sp.Expr, nu: sp.Symbol) -> int:
    poly = sp.Poly(sp.expand(expr), nu)
    return int(poly.degree())


# ---------------------------------------------------------------------------
# Part I. Strict one-body 4PN Hamiltonian gate from the quartic compiler
# ---------------------------------------------------------------------------

def onebody_hamiltonian_gate() -> dict[str, sp.Expr]:
    banner("PART I — STRICT ONE-BODY 4PN HAMILTONIAN GATE")

    U, v, p = sp.symbols("U v p", real=True)

    # Exact one-body Lagrangian coefficients carried from the Schwarzschild audit.
    L1 = sp.Rational(1, 8) * v**4 + sp.Rational(3, 2) * U * v**2 - sp.Rational(1, 2) * U**2
    L2 = sp.Rational(1, 16) * v**6 + sp.Rational(7, 8) * U * v**4 + 2 * U**2 * v**2 + sp.Rational(1, 4) * U**3
    L3 = (
        sp.Rational(5, 128) * v**8
        + sp.Rational(11, 16) * U * v**6
        + sp.Rational(47, 16) * U**2 * v**4
        + sp.Rational(13, 8) * U**3 * v**2
        - sp.Rational(1, 8) * U**4
    )
    L4 = (
        sp.Rational(7, 256) * v**10
        + sp.Rational(75, 128) * U * v**8
        + sp.Rational(59, 16) * U**2 * v**6
        + sp.Rational(203, 32) * U**3 * v**4
        + sp.Rational(31, 32) * U**4 * v**2
        + sp.Rational(1, 16) * U**5
    )

    # Exact quartic compiler in the unit-mass isotropic one-body case.
    v0 = p
    A0 = sp.diff(L1, v).subs(v, v0)
    B0 = sp.diff(L2, v).subs(v, v0)
    D0 = sp.diff(L3, v).subs(v, v0)
    C0 = sp.diff(L1, v, 2).subs(v, v0)
    E0 = sp.diff(L2, v, 2).subs(v, v0)
    T0 = sp.diff(L1, v, 3).subs(v, v0)

    H4 = sp.expand(
        -L4.subs(v, v0)
        + A0 * D0
        + sp.Rational(1, 2) * B0**2
        - B0 * C0 * A0
        - sp.Rational(1, 2) * A0**2 * E0
        + sp.Rational(1, 2) * A0**2 * C0**2
        + sp.Rational(1, 6) * A0**3 * T0
    )

    # Natural one-body Hamiltonian slots.
    gate = {
        "K": sp.simplify(H4.coeff(p, 10)),
        "Q1": sp.simplify(H4.coeff(U, 1).coeff(p, 8)),
        "T1": sp.simplify(H4.coeff(U, 2).coeff(p, 6)),
        "S1": sp.simplify(H4.coeff(U, 3).coeff(p, 4)),
        "U1": sp.simplify(H4.coeff(U, 4).coeff(p, 2)),
        "W1": sp.simplify(H4.coeff(U, 5).coeff(p, 0)),
    }

    print("H4^(1-body) =")
    sp.pprint(H4)
    print("\nSlot coefficients:")
    for key in ["K", "Q1", "T1", "S1", "U1", "W1"]:
        print(f"  {key} = {gate[key]}")

    return gate


# ---------------------------------------------------------------------------
# Part II. Imported reduced COM local 4PN ADM Hamiltonian target
# ---------------------------------------------------------------------------

def local_adm_target() -> dict[str, sp.Expr]:
    banner("PART II — IMPORTED REDUCED COM LOCAL 4PN ADM TARGET")

    nu = sp.symbols("nu", real=True)
    pi = sp.pi

    target = {
        # free kinetic slot
        "K": sp.Rational(7, 256) - sp.Rational(63, 256) * nu + sp.Rational(189, 256) * nu**2 - sp.Rational(105, 128) * nu**3 + sp.Rational(63, 256) * nu**4,

        # G/r block: (p^2)^4, (n.p)^2(p^2)^3, (n.p)^4(p^2)^2, (n.p)^6 p^2, (n.p)^8
        "Q1": sp.Rational(45, 128) - sp.Rational(45, 16) * nu + sp.Rational(423, 64) * nu**2 - sp.Rational(1013, 256) * nu**3 - sp.Rational(35, 128) * nu**4,
        "Q2": -sp.Rational(3, 32) * nu**2 + sp.Rational(23, 64) * nu**3 - sp.Rational(5, 32) * nu**4,
        "Q3": -sp.Rational(9, 64) * nu**2 + sp.Rational(69, 128) * nu**3 - sp.Rational(9, 64) * nu**4,
        "Q4": -sp.Rational(5, 64) * nu**3 - sp.Rational(5, 32) * nu**4,
        "Q5": sp.Rational(35, 256) * nu**3 - sp.Rational(35, 128) * nu**4,

        # G^2/r^2 block: (p^2)^3, (n.p)^2(p^2)^2, (n.p)^4 p^2, (n.p)^6
        "T1": sp.Rational(13, 8) - sp.Rational(791, 64) * nu + sp.Rational(4857, 256) * nu**2 + sp.Rational(2335, 256) * nu**3,
        "T2": sp.Rational(49, 16) * nu - sp.Rational(545, 64) * nu**2 + sp.Rational(1135, 256) * nu**3,
        "T3": -sp.Rational(889, 192) * nu + sp.Rational(9475, 768) * nu**2 - sp.Rational(1649, 768) * nu**3,
        "T4": sp.Rational(369, 160) * nu - sp.Rational(1151, 128) * nu**2 + sp.Rational(10353, 1280) * nu**3,

        # G^3/r^3 block: (p^2)^2, (n.p)^2 p^2, (n.p)^4
        "S1": sp.Rational(105, 32) + (sp.Rational(2749, 8192) * pi**2 - sp.Rational(589189, 19200)) * nu + (sp.Rational(18491, 16384) * pi**2 - sp.Rational(1189789, 28800)) * nu**2 - sp.Rational(553, 128) * nu**3,
        "S2": (sp.Rational(63347, 1600) - sp.Rational(1059, 1024) * pi**2) * nu + (-sp.Rational(127, 3) - sp.Rational(4035, 2048) * pi**2) * nu**2 - sp.Rational(225, 64) * nu**3,
        "S3": (sp.Rational(375, 8192) * pi**2 - sp.Rational(23533, 1280)) * nu + (sp.Rational(57563, 1920) - sp.Rational(38655, 16384) * pi**2) * nu**2 - sp.Rational(381, 128) * nu**3,

        # G^4/r^4 block: p^2, (n.p)^2
        "U1": sp.Rational(105, 32) + (sp.Rational(185761, 19200) - sp.Rational(21837, 8192) * pi**2) * nu + (sp.Rational(672811, 19200) - sp.Rational(158177, 49152) * pi**2) * nu**2,
        "U2": (sp.Rational(3401779, 57600) - sp.Rational(28691, 24576) * pi**2) * nu + (sp.Rational(110099, 49152) * pi**2 - sp.Rational(21827, 3840)) * nu**2,

        # G^5/r^5 static block
        "W1": -sp.Rational(1, 16) + (sp.Rational(6237, 1024) * pi**2 - sp.Rational(169199, 2400)) * nu + (sp.Rational(7403, 3072) * pi**2 - sp.Rational(1256, 45)) * nu**2,
    }

    print("Imported slot coefficients from the fixed local ADM chart:")
    for key in ["K", "Q1", "Q2", "Q3", "Q4", "Q5", "T1", "T2", "T3", "T4", "S1", "S2", "S3", "U1", "U2", "W1"]:
        print(f"  {key} = {sp.factor(target[key])}")

    return target


# ---------------------------------------------------------------------------
# Part III. One-body limit checks
# ---------------------------------------------------------------------------

def onebody_checks(gate: dict[str, sp.Expr], target: dict[str, sp.Expr]) -> None:
    banner("PART III — STRICT ONE-BODY LIMIT OF THE IMPORTED TARGET")

    nu = sp.symbols("nu", real=True)

    expect_zero("free slot K at nu->0", target["K"].subs(nu, 0) - gate["K"])
    expect_zero("G/r leading slot Q1 at nu->0", target["Q1"].subs(nu, 0) - gate["Q1"])
    expect_zero("G^2/r^2 leading slot T1 at nu->0", target["T1"].subs(nu, 0) - gate["T1"])
    expect_zero("G^3/r^3 leading slot S1 at nu->0", target["S1"].subs(nu, 0) - gate["S1"])
    expect_zero("G^4/r^4 leading slot U1 at nu->0", target["U1"].subs(nu, 0) - gate["U1"])
    expect_zero("G^5/r^5 leading slot W1 at nu->0", target["W1"].subs(nu, 0) - gate["W1"])

    for key in ["Q2", "Q3", "Q4", "Q5", "T2", "T3", "T4", "S2", "S3", "U2"]:
        expect_zero(f"subleading slot {key} at nu->0", target[key].subs(nu, 0))


# ---------------------------------------------------------------------------
# Part IV. nu-degree theorem gate by block
# ---------------------------------------------------------------------------

def block_degree_checks(target: dict[str, sp.Expr]) -> None:
    banner("PART IV — HIGHEST nu-DEGREE BY LOCAL 4PN BLOCK")

    nu = sp.symbols("nu", real=True)

    blocks = {
        "free": ["K"],
        "G/r": ["Q1", "Q2", "Q3", "Q4", "Q5"],
        "G^2/r^2": ["T1", "T2", "T3", "T4"],
        "G^3/r^3": ["S1", "S2", "S3"],
        "G^4/r^4": ["U1", "U2"],
        "G^5/r^5": ["W1"],
    }

    degree_by_block: dict[str, int] = {}
    for block, keys in blocks.items():
        degs = [max_nu_degree(target[key], nu) for key in keys]
        degree_by_block[block] = max(degs)
        print(f"{block:8s} max degree in nu = {degree_by_block[block]}  from {dict(zip(keys, degs))}")

    print("\nExact theorem gate:")
    print("  Do the upper local blocks G^4/r^4 or G^5/r^5 contain nu^3 or nu^4 tails?")
    print("  -> No.")

    if degree_by_block["G^4/r^4"] > 2:
        raise AssertionError("Imported G^4/r^4 block unexpectedly contains nu^3 or higher tails.")
    if degree_by_block["G^5/r^5"] > 2:
        raise AssertionError("Imported G^5/r^5 block unexpectedly contains nu^3 or higher tails.")


# ---------------------------------------------------------------------------
# Part V. Readout in natural slot language
# ---------------------------------------------------------------------------

def final_readout(target: dict[str, sp.Expr]) -> None:
    banner("PART V — FINAL LOCAL TARGET-IMPORT LEDGER")

    print("Imported fixed-chart local 4PN COM Hamiltonian slot counts:")
    print("  free kinetic : 1 slot")
    print("  G/r          : 5 slots")
    print("  G^2/r^2      : 4 slots")
    print("  G^3/r^3      : 3 slots")
    print("  G^4/r^4      : 2 slots")
    print("  G^5/r^5      : 1 slot")
    print("  total local  : 16 slots")
    print()
    print("Most important structural outcome:")
    print("  - the fixed local 4PN target is now frozen in the COM Hamiltonian chart;")
    print("  - it matches the strict one-body gate exactly;")
    print("  - the upper local interaction blocks G^4/r^4 and G^5/r^5 stop at nu^2;")
    print("  - so the imported target does not force any nu^3 / nu^4 upper-block tails.")
    print()
    print("That means the next exact local theorem gate is no longer source selection.")
    print("It is the Hamiltonian-to-ordinary translation of these fixed target slots against")
    print("the local 4PN scaffold and the quartic compiler.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    gate = onebody_hamiltonian_gate()
    target = local_adm_target()
    onebody_checks(gate, target)
    block_degree_checks(target)
    final_readout(target)
