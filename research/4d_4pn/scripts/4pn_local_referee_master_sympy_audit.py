#!/usr/bin/env python3
"""
4pn_local_referee_master_sympy_audit.py

End-to-end referee audit for the *local* 4PN derivation chain.

What this script does
---------------------
1. Rebuilds the exact reduced COM local 4PN ordinary target from the imported
   local ADM Hamiltonian target using the exact quartic Legendre compiler.
2. Verifies the strict one-body gate and the natural local self/static seed.
3. Rebuilds the canonical generic-frame Hamiltonian residual, translates it to the
   ordinary chart, and verifies the exact sign-flip theorem

       Delta L4^(can) = - Delta H4^(can).

4. Rebuilds the aligned ordinary seed and the seed-alignment correction

       delta_seed = L4,seed^(aligned) - L4,seed^(natural).

5. Solves the generic-frame lift of delta_seed in the minimal structured seed spaces

       Q,
       T ⊕ (pq)T,
       S ⊕ (pq)S,
       U ⊕ (p^2 q^2)U,
       W = 0.

6. Assembles one exact generic-frame ordinary local 4PN candidate

       L4,loc^(gen)
         = L4,seed^(natural,gen)
         + delta_seed^(gen)
         + Delta L4,loc^(can,gen),

   and verifies that its COM reduction reproduces the full fixed-chart local
   ordinary 4PN target slot by slot.

Scope note
----------
This is a referee audit for the *local* 4PN sector only.  The hereditary/tail bridge
remains separate by construction.
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


def union_exprs(*lists: list[sp.Expr]) -> list[sp.Expr]:
    seen: set[sp.Expr] = set()
    out: list[sp.Expr] = []
    for lst in lists:
        for expr in lst:
            expr = sp.expand(expr)
            if expr not in seen:
                seen.add(expr)
                out.append(expr)
    return sorted(out, key=str)


def canonical_expr(coords: sp.Matrix, basis: list[sp.Expr]) -> sp.Expr:
    return sp.expand(sum(sp.expand(c) * e for c, e in zip(coords, basis)))


def nonzero_count(vec: sp.Matrix) -> int:
    return sum(1 for x in vec if sp.simplify(x) != 0)


# ---------------------------------------------------------------------------
# Load prior exact audits quietly
# ---------------------------------------------------------------------------

BASE = os.path.dirname(os.path.abspath(__file__))
with contextlib.redirect_stdout(io.StringIO()):
    mod_ord = runpy.run_path(os.path.join(BASE, "4pn_local_hamiltonian_to_ordinary_audit.py"))
    mod_can = runpy.run_path(os.path.join(BASE, "4pn_hamiltonian_chart_generic_frame_lift_audit.py"))
    mod_tr = runpy.run_path(os.path.join(BASE, "4pn_generic_frame_ordinary_translation_audit.py"))

nu = mod_can["nu"]
a = mod_can["a"]
b = mod_can["b"]
c = mod_can["c"]
d = mod_can["d"]
e = mod_can["e"]
p = mod_can["p"]
q = mod_can["q"]


# ---------------------------------------------------------------------------
# Slot helpers
# ---------------------------------------------------------------------------

BLOCK_RANGES = {
    "Q": range(7, 12),
    "T": range(12, 16),
    "S": range(16, 19),
    "U": range(19, 21),
    "W": range(21, 22),
}


def block_slots_from_target(target_l: dict[int, sp.Expr], block: str) -> list[sp.Expr]:
    return [target_l[i] for i in BLOCK_RANGES[block]]


# ---------------------------------------------------------------------------
# Rebuild exact local target data
# ---------------------------------------------------------------------------

def rebuild_exact_targets() -> tuple[
    dict[int, sp.Expr],
    dict[int, sp.Expr],
    dict[int, sp.Expr],
    dict[str, sp.Expr],
    dict[str, list[sp.Expr]],
    dict[int, sp.Expr],
]:
    hmap, coeffs, _ = mod_ord["quartic_com_map"]()
    target_h = mod_ord["local_hamiltonian_target"]()
    target_l = mod_ord["ordinary_target_from_hamiltonian"](hmap, coeffs, target_h)
    natural_seed = mod_ord["natural_seed_local_ordinary"]()

    Hblocks, residual_slots_H = mod_tr["canonical_hamiltonian_blocks"]()
    Lblocks_res, ordinary_residual_slots_aligned, aligned_seed, natural_seed_again, target_l_again = mod_tr["ordinary_residual_from_hamiltonian"](Hblocks, residual_slots_H)

    for i in range(1, 22):
        if i in natural_seed_again:
            if sp.expand(sp.simplify(natural_seed_again[i] - natural_seed[i])) != 0:
                raise AssertionError("natural seed mismatch across audits")
        if i in target_l_again:
            if sp.expand(sp.simplify(target_l_again[i] - target_l[i])) != 0:
                raise AssertionError("target ordinary slot mismatch across audits")

    return target_h, target_l, natural_seed, Lblocks_res, ordinary_residual_slots_aligned, aligned_seed


# ---------------------------------------------------------------------------
# Part I. Exact local target, one-body gate, and aligned seed free block
# ---------------------------------------------------------------------------

def target_and_seed_audit(target_h: dict[int, sp.Expr], target_l: dict[int, sp.Expr], natural_seed: dict[int, sp.Expr], aligned_seed: dict[int, sp.Expr]) -> None:
    banner("PART I — EXACT LOCAL 4PN TARGET, ONE-BODY GATE, AND ALIGNED SEED")

    subbanner("I.1 — strict one-body local ordinary gate")
    expect_zero("L1 one-body gate", target_l[1].subs(nu, 0) - sp.Rational(7, 256))
    expect_zero("L7 one-body gate", target_l[7].subs(nu, 0) - sp.Rational(75, 128))
    expect_zero("L12 one-body gate", target_l[12].subs(nu, 0) - sp.Rational(59, 16))
    expect_zero("L16 one-body gate", target_l[16].subs(nu, 0) - sp.Rational(203, 32))
    expect_zero("L19 one-body gate", target_l[19].subs(nu, 0) - sp.Rational(31, 32))
    expect_zero("L21 one-body gate", target_l[21].subs(nu, 0) - sp.Rational(1, 16))
    for i in [2, 3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15, 17, 18, 20]:
        expect_zero(f"subleading one-body slot {i}", sp.simplify(target_l[i].subs(nu, 0)))

    subbanner("I.2 — the natural seed is not yet the full target")
    delta_nat = {i: sp.simplify(target_l[i] - natural_seed[i]) for i in range(1, 22)}
    print("Nonzero residual slots beyond the natural seed:")
    for i in range(1, 22):
        if delta_nat[i] != 0:
            print(f"  Delta_nat[{i}] = {sp.factor(delta_nat[i])}")

    subbanner("I.3 — the aligned seed already fixes the full free block")
    for i in range(1, 7):
        expect_zero(f"aligned seed free slot {i}", aligned_seed[i] - target_l[i])

    print("\nLocal interaction target slots 7–21 remain to be built generically;")
    print("the free slot is already fixed exactly by the Hamiltonian-aligned seed.")


# ---------------------------------------------------------------------------
# Part II. Canonical residual sign-flip theorem and COM slots
# ---------------------------------------------------------------------------

def residual_audit(Lblocks_res: dict[str, sp.Expr], residual_slots_aligned: dict[str, list[sp.Expr]]) -> None:
    banner("PART II — CANONICAL LOCAL RESIDUAL AND THE ORDINARY SIGN-FLIP THEOREM")

    Hblocks, residual_slots_H = mod_tr["canonical_hamiltonian_blocks"]()

    subbanner("II.1 — blockwise sign-flip theorem")
    for block in ["Q", "T", "S", "U", "W"]:
        expect_zero(f"{block} block sign flip", Lblocks_res[block] + Hblocks[block])

    subbanner("II.2 — direct COM slot verification of the canonical ordinary residual")
    for block in ["Q", "T", "S", "U", "W"]:
        lhs = mod_can["block_slots"](mod_can["to_nu"](Lblocks_res[block]), block)
        rhs = residual_slots_aligned[block]
        for i, (x, y) in enumerate(zip(lhs, rhs), start=1):
            expect_zero(f"{block} residual slot {i}", x - y)

    print("\nSo the Hamiltonian-derived generic-frame comparable-mass residual translates back")
    print("to the ordinary chart by exact sign flip, with no remaining local ambiguity.")


# ---------------------------------------------------------------------------
# Part III. Rebuild the exact seed-alignment correction in COM slots
# ---------------------------------------------------------------------------

def seed_alignment_data(target_l: dict[int, sp.Expr], natural_seed: dict[int, sp.Expr], aligned_seed: dict[int, sp.Expr]) -> tuple[dict[int, sp.Expr], dict[str, list[sp.Expr]]]:
    delta_seed = {i: sp.simplify(aligned_seed[i] - natural_seed[i]) for i in range(1, 22)}
    mis_blocks = {block: [delta_seed[i] for i in BLOCK_RANGES[block]] for block in ["Q", "T", "S", "U", "W"]}
    return delta_seed, mis_blocks


def print_seed_alignment_ledger(delta_seed: dict[int, sp.Expr], mis_blocks: dict[str, list[sp.Expr]]) -> None:
    banner("PART III — EXACT SEED-ALIGNMENT CORRECTION")

    print("Free correction slot:")
    print(f"  delta_seed[1] = {sp.factor(delta_seed[1])}")

    print("\nLocal interaction correction slots:")
    for block in ["Q", "T", "S", "U", "W"]:
        print(f"  {block}: {[sp.factor(x) for x in mis_blocks[block]]}")

    print("\nNu-degree ceilings:")
    for block in ["Q", "T", "S", "U", "W"]:
        degs = [mod_tr["max_nu_degree"](x, nu) for x in mis_blocks[block]]
        print(f"  {block}: {degs}")


# ---------------------------------------------------------------------------
# Part IV. Minimal structured generic-frame lift of delta_seed
# ---------------------------------------------------------------------------

def structured_seed_lift(mis_blocks: dict[str, list[sp.Expr]]) -> tuple[dict[str, list[sp.Expr]], dict[str, sp.Matrix], dict[str, sp.Expr]]:
    banner("PART IV — STRUCTURED GENERIC-FRAME LIFT OF THE SEED-ALIGNMENT CORRECTION")

    old_bases = {
        "Q": mod_can["generate_basis"](0, 8),
        "T": mod_can["generate_basis"](1, 6),
        "S": mod_can["generate_basis"](2, 4),
        "U": mod_can["generate_basis"](3, 2),
        "W": mod_can["generate_basis"](4, 0),
    }

    pq = p * q
    seed_bases = {
        "Q": old_bases["Q"],
        "T": union_exprs(old_bases["T"], [sp.expand(pq * expr) for expr in old_bases["T"]]),
        "S": union_exprs(old_bases["S"], [sp.expand(pq * expr) for expr in old_bases["S"]]),
        "U": union_exprs(old_bases["U"], [sp.expand((pq ** 2) * expr) for expr in old_bases["U"]]),
        "W": [],
    }

    print("Lift spaces:")
    print("  Q_delta in Q")
    print("  T_delta in T ⊕ (pq)T")
    print("  S_delta in S ⊕ (pq)S")
    print("  U_delta in U ⊕ (p^2 q^2)U")
    print("  W_delta = 0")

    coords: dict[str, sp.Matrix] = {}
    exprs: dict[str, sp.Expr] = {"W": sp.Integer(0)}

    for block, maxdeg in [("Q", 4), ("T", 4), ("S", 4), ("U", 4)]:
        basis = seed_bases[block]
        M = mod_can["image_matrix_polynomial"](block, basis, maxdeg)
        vec = mod_can["target_vector"](mis_blocks[block], maxdeg)
        sol = mod_can["particular_solution"](M, vec)
        coords[block] = sol
        exprs[block] = canonical_expr(sol, basis)
        expect_zero(f"{block} structured solve", M * sol - vec)
        print(f"{block}: basis size = {len(basis)}, nonzero canonical directions = {nonzero_count(sol)}")

    subbanner("IV.1 — direct COM verification of delta_seed^(gen)")
    for block in ["Q", "T", "S", "U"]:
        lhs = mod_can["block_slots"](mod_can["to_nu"](exprs[block]), block)
        rhs = mis_blocks[block]
        for i, (x, y) in enumerate(zip(lhs, rhs), start=1):
            expect_zero(f"{block} delta slot {i}", x - y)
    expect_zero("W delta is identically zero", mis_blocks["W"][0])

    return seed_bases, coords, exprs


# ---------------------------------------------------------------------------
# Part V. Assemble one exact generic-frame ordinary local candidate
# ---------------------------------------------------------------------------

def assemble_full_candidate(target_l: dict[int, sp.Expr], natural_seed: dict[int, sp.Expr], aligned_seed: dict[int, sp.Expr], delta_exprs: dict[str, sp.Expr], residual_blocks: dict[str, sp.Expr]) -> dict[str, sp.Expr]:
    banner("PART V — ASSEMBLE ONE EXACT GENERIC-FRAME ORDINARY LOCAL 4PN CANDIDATE")

    Qnat = sp.Rational(75, 128) * (a**4 + b**4)
    Tnat = sp.Rational(59, 16) * (q * a**3 + p * b**3)
    Snat = sp.Rational(203, 32) * (q**2 * a**2 + p**2 * b**2)
    Unat = sp.Rational(31, 32) * (q**3 * a + p**3 * b)
    Wnat = sp.Rational(1, 16) * (q**4 + p**4)

    aligned_blocks = {
        "Q": sp.expand(Qnat + delta_exprs["Q"]),
        "T": sp.expand(Tnat + delta_exprs["T"]),
        "S": sp.expand(Snat + delta_exprs["S"]),
        "U": sp.expand(Unat + delta_exprs["U"]),
        "W": sp.expand(Wnat),
    }

    subbanner("V.1 — direct COM check of the aligned seed blocks")
    for block, slots in [
        ("Q", [aligned_seed[i] for i in BLOCK_RANGES["Q"]]),
        ("T", [aligned_seed[i] for i in BLOCK_RANGES["T"]]),
        ("S", [aligned_seed[i] for i in BLOCK_RANGES["S"]]),
        ("U", [aligned_seed[i] for i in BLOCK_RANGES["U"]]),
        ("W", [aligned_seed[21]]),
    ]:
        lhs = mod_can["block_slots"](mod_can["to_nu"](aligned_blocks[block]), block)
        for i, (x, y) in enumerate(zip(lhs, slots), start=1):
            expect_zero(f"aligned {block} slot {i}", x - y)

    full_blocks = {
        "Q": sp.expand(aligned_blocks["Q"] + residual_blocks["Q"]),
        "T": sp.expand(aligned_blocks["T"] + residual_blocks["T"]),
        "S": sp.expand(aligned_blocks["S"] + residual_blocks["S"]),
        "U": sp.expand(aligned_blocks["U"] + residual_blocks["U"]),
        "W": sp.expand(aligned_blocks["W"] + residual_blocks["W"]),
    }

    subbanner("V.2 — full direct COM verification against the exact local ordinary target")
    for block in ["Q", "T", "S", "U", "W"]:
        lhs = mod_can["block_slots"](mod_can["to_nu"](full_blocks[block]), block)
        rhs = block_slots_from_target(target_l, block)
        for i, (x, y) in enumerate(zip(lhs, rhs), start=1):
            expect_zero(f"full {block} slot {i}", x - y)

    print("\nFull local generic-frame ordinary candidate:")
    print("  L4,loc^(gen) = L4,seed^(natural,gen) + delta_seed^(gen) + Delta L4,loc^(can,gen)")
    print("with the free block already fixed by the aligned seed and the interaction blocks")
    print("reproducing the exact local ordinary target after COM reduction.")

    return full_blocks


# ---------------------------------------------------------------------------
# Part VI. Final theorem ledger
# ---------------------------------------------------------------------------

def final_ledger(seed_bases: dict[str, list[sp.Expr]], coords: dict[str, sp.Matrix], full_blocks: dict[str, sp.Expr]) -> None:
    banner("PART VI — FINAL LOCAL 4PN THEOREM LEDGER")

    print("Structured lift sizes / canonical sparsities:")
    for block in ["Q", "T", "S", "U"]:
        print(f"  {block}: basis size {len(seed_bases[block])}, nonzero directions {nonzero_count(coords[block])}")
    print("  W: no lift required")

    print("\nMain local result:")
    print("  1. The full fixed-chart local 4PN ordinary target is reconstructed exactly from the")
    print("     quartic Hamiltonian target plus the exact Hamiltonian->ordinary compiler.")
    print("  2. The canonical comparable-mass local residual lifts cleanly in the Hamiltonian chart")
    print("     and translates back by exact sign flip.")
    print("  3. The only ordinary-chart obstruction is the aligned-seed correction, and that")
    print("     correction is exactly generic-frame liftable in the minimal structured spaces")
    print("        Q,  T ⊕ (pq)T,  S ⊕ (pq)S,  U ⊕ (p^2 q^2)U,  W=0.")
    print("  4. Therefore one exact generic-frame ordinary representative of the *entire local* 4PN")
    print("     sector exists and has now been assembled explicitly.")
    print("  5. What remains open is not local existence.  It is the sharper uniqueness/")
    print("     interpretation question for the aligned-seed sector, and the separate tail bridge.")

    print("\nFinal assembled interaction blocks:")
    for block in ["Q", "T", "S", "U", "W"]:
        print(f"  {block}(gen) =")
        sp.pprint(sp.factor(full_blocks[block]))
        print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    target_h, target_l, natural_seed, residual_blocks, residual_slots_aligned, aligned_seed = rebuild_exact_targets()
    target_and_seed_audit(target_h, target_l, natural_seed, aligned_seed)
    residual_audit(residual_blocks, residual_slots_aligned)
    delta_seed, mis_blocks = seed_alignment_data(target_l, natural_seed, aligned_seed)
    print_seed_alignment_ledger(delta_seed, mis_blocks)
    seed_bases, coords, delta_exprs = structured_seed_lift(mis_blocks)
    full_blocks = assemble_full_candidate(target_l, natural_seed, aligned_seed, delta_exprs, residual_blocks)
    final_ledger(seed_bases, coords, full_blocks)


if __name__ == "__main__":
    main()
