#!/usr/bin/env python3
"""
4pn_generic_frame_ordinary_translation_audit.py

Translate the canonical local 4PN generic-frame Hamiltonian representative back to the
ordinary chart, while keeping the tail bridge separate.

What this script does
---------------------
1. Reuses the exact quartic-compiler facts already frozen earlier and states the key
   residual theorem:

       if the lower-order ordinary ledger (through 3PN) is frozen, then for any 4PN
       ordinary residual ΔL4 written in the same formal generic-frame scaffold,

           ΔH4 = -ΔL4(v0)

   once the 4PN seed is chosen in the Hamiltonian-aligned convention.

2. Rebuilds the sparse canonical generic-frame Hamiltonian representative found in the
   previous audit:

       ΔH4,loc^can = (Gpq/r) Q_can + (G^2 pq/r^2) T_can + (G^3 pq/r^3) S_can
                    + (G^4/r^4) U_can + (G^5/r^5) W_can.

3. Translates that representative to the ordinary chart by the exact quartic residual
   sign flip:

       ΔL4,loc^can = - ΔH4,loc^can.

4. Verifies directly that the COM reduction of this canonical ordinary residual matches
   the exact fixed-chart *ordinary* residual once the latter is measured relative to the
   ordinary seed aligned with the chosen Hamiltonian self/static seed.

5. Computes the aligned ordinary COM seed and compares it to the earlier natural
   one-body ordinary seed, showing that the previous ordinary-chart "nu^4 obstruction"
   is entirely a seed-alignment issue rather than a failure of the Hamiltonian generic-
   frame lift.

Main structural result
----------------------
The canonical local 4PN generic-frame Hamiltonian representative does translate cleanly
back to a canonical generic-frame ordinary *residual* representative.  The exact local
ordinary target splits as

    L4,target = L4,seed^(aligned) + ΔL4,loc^can,

where ΔL4,loc^can is the translated residual above.  The earlier ordinary-chart span
problem is therefore not a failure of the comparable-mass lift; it is the separate chart
bookkeeping problem of lifting the Hamiltonian-aligned ordinary seed (or equivalently the
seed-alignment correction beyond the natural one-body ordinary seed).
"""

from __future__ import annotations

import contextlib
import io
import os
import runpy
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


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        simplified = expr.applyfunc(lambda z: sp.expand(sp.simplify(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.expand(sp.simplify(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


def max_nu_degree(expr: sp.Expr, nu: sp.Symbol) -> int:
    if sp.expand(expr) == 0:
        return -1
    return int(sp.Poly(sp.expand(expr), nu).degree())


def block_degree_vector(block: list[sp.Expr], nu: sp.Symbol) -> list[int]:
    return [max_nu_degree(x, nu) for x in block]


def nonzero_slot_dict(slot_map: dict[int, sp.Expr]) -> dict[int, sp.Expr]:
    return {i: sp.factor(sp.simplify(v)) for i, v in slot_map.items() if sp.simplify(v) != 0}


# ---------------------------------------------------------------------------
# Load prior exact audits without executing their main blocks
# ---------------------------------------------------------------------------

BASE = os.path.dirname(os.path.abspath(__file__))
mod_can = runpy.run_path(os.path.join(BASE, "4pn_hamiltonian_chart_canonical_slice_audit.py"))
mod_ord = runpy.run_path(os.path.join(BASE, "4pn_local_hamiltonian_to_ordinary_audit.py"))

nu = mod_can["nu"]
sp_pi = sp.pi


# ---------------------------------------------------------------------------
# Part I. Exact quartic residual compiler theorem
# ---------------------------------------------------------------------------

def residual_compiler_theorem() -> None:
    banner("PART I — EXACT QUARTIC RESIDUAL COMPILER THEOREM")

    d1, d2, s1, s2, v0, F0 = sp.symbols("d1 d2 s1 s2 v0 F0", real=True)

    # Once the lower-order ledger is frozen, the quartic compiler has the form
    #   H4 = F[L1,L2,L3](v0) - L4(v0),
    # where F depends only on the lower orders.
    L4_seed = s1 * v0 + s2 * v0**2
    dL4 = d1 * v0 + d2 * v0**2

    H4_seed = F0 - L4_seed
    H4_full = F0 - (L4_seed + dL4)

    expect_zero("H4(full) - H4(seed) + dL4(v0)", sp.expand(H4_full - H4_seed + dL4))

    print("\nTherefore, once the lower-order 1PN/2PN/3PN ordinary ledger is frozen and the 4PN seed")
    print("is chosen in the Hamiltonian-aligned convention, the quartic residual compiler is")
    print("exactly the same sign flip as at 3PN:")
    print("\n    ΔH4 = -ΔL4(v0)\n")
    print("on the common formal generic-frame scaffold, with the ordinary velocity invariants")
    print("reinterpreted on the Hamiltonian side as Newtonian-order momentum invariants.")


# ---------------------------------------------------------------------------
# Part II. Rebuild the sparse canonical Hamiltonian representative
# ---------------------------------------------------------------------------

def canonical_hamiltonian_blocks() -> tuple[dict[str, sp.Expr], dict[str, list[sp.Expr]]]:
    banner("PART II — REBUILD THE CANONICAL GENERIC-FRAME HAMILTONIAN RESIDUAL")

    with contextlib.redirect_stdout(io.StringIO()):
        Qbasis = mod_can["generate_basis"](0, 8)
        Tbasis = mod_can["generate_basis"](1, 6)
        Sbasis = mod_can["generate_basis"](2, 4)
        Ubasis = mod_can["generate_basis"](3, 2)
        Wbasis = mod_can["generate_basis"](4, 0)

        QM = mod_can["image_matrix_polynomial"]("Q", Qbasis, 4)
        TM = mod_can["image_matrix_polynomial"]("T", Tbasis, 3)
        SM = mod_can["image_matrix_polynomial"]("S", Sbasis, 3)
        UM = mod_can["image_matrix_polynomial"]("U", Ubasis, 2)
        WM = mod_can["image_matrix_polynomial"]("W", Wbasis, 2)

        target = mod_can["local_adm_target"]()
        seedH = mod_can["hamiltonian_seed"]()
        residual = mod_can["residual_slots"](target, seedH)

        Qcoords = mod_can["particular_solution"](QM, mod_can["target_vector"](residual["Q"], 4))
        Tcoords = mod_can["particular_solution"](TM, mod_can["target_vector"](residual["T"], 3))
        Scoords = mod_can["particular_solution"](SM, mod_can["target_vector"](residual["S"], 3))
        Ucoords = mod_can["particular_solution"](UM, mod_can["target_vector"](residual["U"], 2))
        Wcoords = mod_can["particular_solution"](WM, mod_can["target_vector"](residual["W"], 2))

    Qcan = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Qcoords, Qbasis)))
    Tcan = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Tcoords, Tbasis)))
    Scan = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Scoords, Sbasis)))
    Ucan = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Ucoords, Ubasis)))
    Wcan = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Wcoords, Wbasis)))

    print("Nonzero scaffold directions in the carried canonical Hamiltonian representative:")
    print(f"  Q : {len(mod_can['nonzero_terms'](Qcoords, Qbasis))}")
    print(f"  T : {len(mod_can['nonzero_terms'](Tcoords, Tbasis))}")
    print(f"  S : {len(mod_can['nonzero_terms'](Scoords, Sbasis))}")
    print(f"  U : {len(mod_can['nonzero_terms'](Ucoords, Ubasis))}")
    print(f"  W : {len(mod_can['nonzero_terms'](Wcoords, Wbasis))}")

    Hblocks = {"Q": Qcan, "T": Tcan, "S": Scan, "U": Ucan, "W": Wcan}
    residual_slots = {
        "Q": residual["Q"],
        "T": residual["T"],
        "S": residual["S"],
        "U": residual["U"],
        "W": residual["W"],
    }
    return Hblocks, residual_slots


# ---------------------------------------------------------------------------
# Part III. Translate to the ordinary residual representative
# ---------------------------------------------------------------------------

def ordinary_residual_from_hamiltonian(Hblocks: dict[str, sp.Expr], residual_slots_H: dict[str, list[sp.Expr]]) -> tuple[dict[str, sp.Expr], dict[str, list[sp.Expr]], dict[int, sp.Expr], dict[int, sp.Expr], dict[int, sp.Expr]]:
    banner("PART III — CANONICAL GENERIC-FRAME ORDINARY RESIDUAL")

    Lblocks = {k: sp.expand(-v) for k, v in Hblocks.items()}

    print("The translated canonical ordinary residual is obtained by the exact quartic sign flip:")
    print("\n    Q_can^(L) = - Q_can^(H)")
    print("    T_can^(L) = - T_can^(H)")
    print("    S_can^(L) = - S_can^(H)")
    print("    U_can^(L) = - U_can^(H)")
    print("    W_can^(L) = - W_can^(H)\n")

    # Extract the exact fixed-chart ordinary target and the Hamiltonian-aligned ordinary seed.
    with contextlib.redirect_stdout(io.StringIO()):
        hmap, coeffs, _ = mod_ord["quartic_com_map"]()
        target_h = mod_ord["local_hamiltonian_target"]()
        target_l = mod_ord["ordinary_target_from_hamiltonian"](hmap, coeffs, target_h)
        natural_seed = mod_ord["natural_seed_local_ordinary"]()
        seedH_named = mod_can["hamiltonian_seed"]()

    feedback = {i: sp.simplify(hmap[i].subs({coeffs[j]: 0 for j in range(21)})) for i in range(1, 22)}
    seedH_idx = {
        1: seedH_named["K"],
        2: sp.Integer(0), 3: sp.Integer(0), 4: sp.Integer(0), 5: sp.Integer(0), 6: sp.Integer(0),
        7: seedH_named["Q1"], 8: seedH_named["Q2"], 9: seedH_named["Q3"], 10: seedH_named["Q4"], 11: seedH_named["Q5"],
        12: seedH_named["T1"], 13: seedH_named["T2"], 14: seedH_named["T3"], 15: seedH_named["T4"],
        16: seedH_named["S1"], 17: seedH_named["S2"], 18: seedH_named["S3"],
        19: seedH_named["U1"], 20: seedH_named["U2"],
        21: seedH_named["W1"],
    }
    aligned_seed = {i: sp.simplify(feedback[i] - seedH_idx[i]) for i in range(1, 22)}

    ordinary_residual_slots_aligned = {
        "Q": [sp.simplify(target_l[i] - aligned_seed[i]) for i in range(7, 12)],
        "T": [sp.simplify(target_l[i] - aligned_seed[i]) for i in range(12, 16)],
        "S": [sp.simplify(target_l[i] - aligned_seed[i]) for i in range(16, 19)],
        "U": [sp.simplify(target_l[i] - aligned_seed[i]) for i in range(19, 21)],
        "W": [sp.simplify(target_l[21] - aligned_seed[21])],
    }

    subbanner("III.1 — Direct COM verification against the aligned ordinary residual")
    for label in ["Q", "T", "S", "U", "W"]:
        lhs = mod_can["block_slots"](mod_can["to_nu"](Lblocks[label]), label)
        rhs = ordinary_residual_slots_aligned[label]
        for i, (x, y) in enumerate(zip(lhs, rhs), start=1):
            expect_zero(f"{label} slot {i}", x - y)

    subbanner("III.2 — Residual nu-degree ceilings are unchanged by the chart translation")
    for label in ["Q", "T", "S", "U", "W"]:
        print(f"{label} Hamiltonian residual degrees = {block_degree_vector(residual_slots_H[label], nu)}")
        print(f"{label} ordinary residual degrees    = {block_degree_vector(ordinary_residual_slots_aligned[label], nu)}")

    return Lblocks, ordinary_residual_slots_aligned, aligned_seed, natural_seed, target_l


# ---------------------------------------------------------------------------
# Part IV. Seed alignment and the earlier ordinary-chart obstruction
# ---------------------------------------------------------------------------

def seed_alignment_analysis(aligned_seed: dict[int, sp.Expr], natural_seed: dict[int, sp.Expr], target_l: dict[int, sp.Expr], ordinary_residual_slots_aligned: dict[str, list[sp.Expr]]) -> None:
    banner("PART IV — THE ORDINARY-SIDE SEED-ALIGNMENT CORRECTION")

    misalign = {i: sp.simplify(aligned_seed[i] - natural_seed[i]) for i in range(1, 22)}
    nz_align = nonzero_slot_dict({i: aligned_seed[i] for i in range(1, 22)})
    nz_mis = nonzero_slot_dict({i: misalign[i] for i in range(1, 22)})

    subbanner("IV.1 — Hamiltonian-aligned ordinary COM seed coefficients")
    for i, expr in nz_align.items():
        print(f"L{i}^(seed,aligned) = {expr}")

    subbanner("IV.2 — Seed-alignment correction relative to the natural one-body ordinary seed")
    for i, expr in nz_mis.items():
        print(f"delta_seed[{i}] = {expr}")

    mis_blocks = {
        "Q": [misalign[i] for i in range(7, 12)],
        "T": [misalign[i] for i in range(12, 16)],
        "S": [misalign[i] for i in range(16, 19)],
        "U": [misalign[i] for i in range(19, 21)],
        "W": [misalign[21]],
    }

    subbanner("IV.3 — Degree ceilings of the seed-alignment correction")
    for label in ["Q", "T", "S", "U", "W"]:
        print(f"{label} correction degrees = {block_degree_vector(mis_blocks[label], nu)}")

    subbanner("IV.4 — Exact decomposition of the full ordinary target")
    # Rebuild target blocks.
    target_blocks = {
        "Q": [target_l[i] for i in range(7, 12)],
        "T": [target_l[i] for i in range(12, 16)],
        "S": [target_l[i] for i in range(16, 19)],
        "U": [target_l[i] for i in range(19, 21)],
        "W": [target_l[21]],
    }
    natural_blocks = {
        "Q": [natural_seed[i] for i in range(7, 12)],
        "T": [natural_seed[i] for i in range(12, 16)],
        "S": [natural_seed[i] for i in range(16, 19)],
        "U": [natural_seed[i] for i in range(19, 21)],
        "W": [natural_seed[21]],
    }

    for label in ["Q", "T", "S", "U", "W"]:
        for i, tgt in enumerate(target_blocks[label], start=1):
            rhs = natural_blocks[label][i - 1] + mis_blocks[label][i - 1] + ordinary_residual_slots_aligned[label][i - 1]
            expect_zero(f"{label} full-target split slot {i}", tgt - rhs)

    print("\nSo the exact fixed-chart local ordinary target splits as")
    print("\n    L4,target = L4,seed^(natural) + delta_seed + ΔL4,loc^can,\n")
    print("with delta_seed = L4,seed^(aligned) - L4,seed^(natural).")
    print("The earlier ordinary-chart nu^4 obstruction is therefore not a failure of the")
    print("Hamiltonian generic-frame lift.  It is the separate chart/seed object delta_seed,")
    print("whose decisive new nu^4 content sits in the G^2/r^2, G^3/r^3, and G^4/r^4 blocks.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    residual_compiler_theorem()
    Hblocks, residual_slots_H = canonical_hamiltonian_blocks()
    Lblocks, ordinary_residual_slots_aligned, aligned_seed, natural_seed, target_l = ordinary_residual_from_hamiltonian(Hblocks, residual_slots_H)
    seed_alignment_analysis(aligned_seed, natural_seed, target_l, ordinary_residual_slots_aligned)

    banner("PART V — FINAL LEDGER")
    print("1. The quartic generic-frame residual compiler is exactly a sign flip once the lower-order")
    print("   ordinary ledger is frozen and the 4PN seed is chosen in the Hamiltonian-aligned convention:")
    print("       ΔH4 = -ΔL4.")
    print("2. Therefore the sparse canonical Hamiltonian representative from the previous audit")
    print("   translates immediately to a sparse canonical generic-frame ordinary residual representative.")
    print("3. Its COM reduction matches the exact fixed-chart ordinary residual measured relative to")
    print("   the Hamiltonian-aligned ordinary seed.")
    print("4. The earlier ordinary-chart nu^4 obstruction is not a failure of the comparable-mass")
    print("   Hamiltonian lift.  It is the separate seed-alignment correction beyond the natural")
    print("   one-body ordinary seed.")
    print("5. In that precise sense, the next local chart problem is no longer residual existence.")
    print("   It is the explicit generic-frame realization of the aligned ordinary seed / seed-")
    print("   alignment correction, while the nonlocal tail bridge remains separate.")


if __name__ == "__main__":
    main()
