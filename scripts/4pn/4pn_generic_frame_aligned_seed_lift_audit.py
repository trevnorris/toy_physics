#!/usr/bin/env python3
"""
4pn_generic_frame_aligned_seed_lift_audit.py

Lift the aligned ordinary local 4PN seed / seed-alignment correction to the generic frame.

What this script does
---------------------
1. Rebuilds the exact COM seed-alignment correction

       delta_seed = L4,seed^(aligned) - L4,seed^(natural)

   from the earlier Hamiltonian->ordinary translation audit.

2. Shows that the old constant-coefficient *comparable-mass interaction* scaffold does
   NOT lift this correction in the ordinary chart:
      - Q block is solvable already,
      - T/S/U blocks are not,
      - W block vanishes.

3. Identifies a minimal structured seed-sector enlargement that *does* lift the local
   correction exactly:

      Q_delta in Q,
      T_delta in T ⊕ (pq) T,
      S_delta in S ⊕ (pq) S,
      U_delta in U ⊕ (p^2 q^2) U,
      W_delta = 0.

   Here Q,T,S,U are the original local generic-frame block families and p,q are the
   mass-fraction variables used in the existing scaffold audits.

4. Constructs one exact canonical generic-frame representative of the local seed-alignment
   correction in those structured enlarged seed spaces.

5. Combines that lift with the natural generic-frame local seed and the already solved
   canonical ordinary residual to produce one exact generic-frame local ordinary candidate
   whose COM reduction reproduces the fixed-chart local 4PN ordinary target exactly.

Main structural result
----------------------
The earlier ordinary-chart obstruction is not a failure of the local generic-frame lift.
It is an *aligned-seed* issue, and the aligned ordinary local seed is generically liftable
once the seed sector is enlarged in the specific structured way above.  In particular,
no enlargement is needed in Q, a single pq-dressed copy suffices in T and S, a p^2 q^2-
 dressed copy suffices in U, and W is unchanged.
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


def nonzero_terms(coords: sp.Matrix, basis: list[sp.Expr]) -> list[tuple[int, sp.Expr, sp.Expr]]:
    out: list[tuple[int, sp.Expr, sp.Expr]] = []
    for i, (coef, expr) in enumerate(zip(coords, basis)):
        coef = sp.expand(sp.simplify(coef))
        if coef != 0:
            out.append((i, sp.factor(coef), expr))
    return out


# ---------------------------------------------------------------------------
# Load prior exact audits quietly
# ---------------------------------------------------------------------------

BASE = "/mnt/data"
with contextlib.redirect_stdout(io.StringIO()):
    mod_can = runpy.run_path(os.path.join(BASE, "4pn_hamiltonian_chart_generic_frame_lift_audit.py"))
    mod_ord = runpy.run_path(os.path.join(BASE, "4pn_local_hamiltonian_to_ordinary_audit.py"))
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
# Rebuild the exact COM seed-alignment data
# ---------------------------------------------------------------------------

def rebuild_seed_alignment_data() -> tuple[dict[int, sp.Expr], dict[int, sp.Expr], dict[int, sp.Expr], dict[str, list[sp.Expr]], dict[str, sp.Expr]]:
    with contextlib.redirect_stdout(io.StringIO()):
        hmap, coeffs, _ = mod_ord["quartic_com_map"]()
        target_h = mod_ord["local_hamiltonian_target"]()
        target_l = mod_ord["ordinary_target_from_hamiltonian"](hmap, coeffs, target_h)
        natural_seed = mod_ord["natural_seed_local_ordinary"]()
        Hblocks, residual_slots_H = mod_tr["canonical_hamiltonian_blocks"]()
        Lblocks_res, ordinary_residual_slots_aligned, aligned_seed, _, _ = mod_tr["ordinary_residual_from_hamiltonian"](Hblocks, residual_slots_H)

    delta_seed = {i: sp.simplify(aligned_seed[i] - natural_seed[i]) for i in range(1, 22)}
    mis_blocks = {
        "Q": [delta_seed[i] for i in range(7, 12)],
        "T": [delta_seed[i] for i in range(12, 16)],
        "S": [delta_seed[i] for i in range(16, 19)],
        "U": [delta_seed[i] for i in range(19, 21)],
        "W": [delta_seed[21]],
    }
    return target_l, natural_seed, aligned_seed, mis_blocks, Lblocks_res


# ---------------------------------------------------------------------------
# Part I. Exact COM correction ledger
# ---------------------------------------------------------------------------

def print_com_correction_ledger(delta_seed: dict[int, sp.Expr], mis_blocks: dict[str, list[sp.Expr]]) -> None:
    banner("PART I — EXACT COM SEED-ALIGNMENT CORRECTION LEDGER")

    print("Pure-kinetic COM correction (separate from the local interaction lift):")
    print(f"  delta_seed[1] = {sp.factor(delta_seed[1])}")

    subbanner("I.1 — Local interaction correction slots")
    for block, idxs in [("Q", range(7, 12)), ("T", range(12, 16)), ("S", range(16, 19)), ("U", range(19, 21))]:
        print(f"\n{block}-block:")
        for i in idxs:
            print(f"  delta_seed[{i}] = {sp.factor(delta_seed[i])}")
    print(f"\nW-block: delta_seed[21] = {delta_seed[21]}")

    print("\nDegree ceilings in nu:")
    for block in ["Q", "T", "S", "U", "W"]:
        degs = [mod_tr["max_nu_degree"](x, nu) for x in mis_blocks[block]]
        print(f"  {block}: {degs}")

    print("\nInterpretation:")
    print("  - The local seed-alignment correction is nontrivial in Q/T/S/U and vanishes in W.")
    print("  - The T/S/U corrections carry higher nu-degree than the original ordinary")
    print("    comparable-mass interaction scaffold allows, so a direct lift in that old")
    print("    scaffold is not expected to work.")


# ---------------------------------------------------------------------------
# Part II. Old scaffold failure
# ---------------------------------------------------------------------------

def old_scaffold_failure(mis_blocks: dict[str, list[sp.Expr]]) -> tuple[dict[str, list[sp.Expr]], dict[str, sp.Matrix]]:
    banner("PART II — THE OLD ORDINARY COMPARABLE-MASS SCAFFOLD FAILS ON delta_seed")

    Qbasis = mod_can["generate_basis"](0, 8)
    Tbasis = mod_can["generate_basis"](1, 6)
    Sbasis = mod_can["generate_basis"](2, 4)
    Ubasis = mod_can["generate_basis"](3, 2)
    Wbasis = mod_can["generate_basis"](4, 0)

    old_bases = {"Q": Qbasis, "T": Tbasis, "S": Sbasis, "U": Ubasis, "W": Wbasis}
    old_coords: dict[str, sp.Matrix] = {}

    for block, maxdeg in [("Q", 4), ("T", 4), ("S", 4), ("U", 4), ("W", 2)]:
        basis = old_bases[block]
        M = mod_can["image_matrix_polynomial"](block, basis, maxdeg)
        vec = mod_can["target_vector"](mis_blocks[block], maxdeg)
        rank = M.rank()
        rows, cols = M.shape
        print(f"{block}: matrix {rows} x {cols}, rank = {rank}, nullity = {cols-rank}")
        try:
            sol = mod_can["particular_solution"](M, vec)
            old_coords[block] = sol
            expect_zero(f"{block} old-image solve", M * sol - vec)
            print(f"  -> solvable in the old scaffold")
        except Exception as exc:
            print(f"  -> NOT solvable in the old scaffold ({exc.__class__.__name__})")

    print("\nConclusion:")
    print("  - Q delta already lies in the old Q family.")
    print("  - W delta vanishes identically.")
    print("  - T/S/U delta do NOT lie in the old constant-coefficient comparable-mass")
    print("    ordinary interaction scaffold.")

    return old_bases, old_coords


# ---------------------------------------------------------------------------
# Part III. Minimal structured seed-sector enlargement
# ---------------------------------------------------------------------------

def structured_seed_spaces(mis_blocks: dict[str, list[sp.Expr]], old_bases: dict[str, list[sp.Expr]]) -> tuple[dict[str, list[sp.Expr]], dict[str, sp.Matrix], dict[str, sp.Expr]]:
    banner("PART III — MINIMAL STRUCTURED GENERIC-FRAME SEED-SPACES")

    pq = p * q

    seed_bases = {
        "Q": old_bases["Q"],
        "T": union_exprs(old_bases["T"], [sp.expand(pq * expr) for expr in old_bases["T"]]),
        "S": union_exprs(old_bases["S"], [sp.expand(pq * expr) for expr in old_bases["S"]]),
        "U": union_exprs(old_bases["U"], [sp.expand((pq ** 2) * expr) for expr in old_bases["U"]]),
        "W": [],
    }

    print("Structured lift spaces:")
    print("  Q_delta ∈ Q")
    print("  T_delta ∈ T ⊕ (pq) T")
    print("  S_delta ∈ S ⊕ (pq) S")
    print("  U_delta ∈ U ⊕ (p^2 q^2) U")
    print("  W_delta = 0")

    coords: dict[str, sp.Matrix] = {}
    exprs: dict[str, sp.Expr] = {"W": sp.Integer(0)}

    for block, maxdeg in [("Q", 4), ("T", 4), ("S", 4), ("U", 4)]:
        basis = seed_bases[block]
        M = mod_can["image_matrix_polynomial"](block, basis, maxdeg)
        vec = mod_can["target_vector"](mis_blocks[block], maxdeg)
        rank = M.rank()
        rows, cols = M.shape
        print(f"\n{block}: matrix {rows} x {cols}, rank = {rank}, nullity = {cols-rank}")
        sol = mod_can["particular_solution"](M, vec)
        expect_zero(f"{block} structured-image solve", M * sol - vec)
        coords[block] = sol
        exprs[block] = canonical_expr(sol, basis)
        print(f"  nonzero canonical directions = {sum(1 for x in sol if sp.simplify(x) != 0)}")

    subbanner("III.1 — Direct COM verification")
    for block, expr in [("Q", exprs["Q"]), ("T", exprs["T"]), ("S", exprs["S"]), ("U", exprs["U"]), ("W", exprs["W"] )]:
        if block == "W":
            expect_zero("W delta direct COM verification", mis_blocks["W"][0])
            continue
        lhs = mod_can["block_slots"](mod_can["to_nu"](expr), block)
        rhs = mis_blocks[block]
        for i, (x, y) in enumerate(zip(lhs, rhs), start=1):
            expect_zero(f"{block} slot {i}", x - y)

    return seed_bases, coords, exprs


# ---------------------------------------------------------------------------
# Part IV. Canonical representative of the local seed-alignment correction
# ---------------------------------------------------------------------------

def print_canonical_representative(seed_bases: dict[str, list[sp.Expr]], coords: dict[str, sp.Matrix], exprs: dict[str, sp.Expr]) -> None:
    banner("PART IV — CANONICAL GENERIC-FRAME REPRESENTATIVE OF delta_seed")

    for block in ["Q", "T", "S", "U"]:
        subbanner(f"IV.{['Q','T','S','U'].index(block)+1} — {block}-block nonzero directions")
        nz = nonzero_terms(coords[block], seed_bases[block])
        print(f"nonzero directions = {len(nz)}")
        for i, coef, expr in nz:
            print(f"  [{i:02d}] {coef}   *   {expr}")
        print("\nCompact block formula:")
        print(sp.factor(exprs[block]))

    subbanner("IV.5 — W block")
    print("W_delta = 0")


# ---------------------------------------------------------------------------
# Part V. Assemble the aligned local ordinary seed and one full local candidate
# ---------------------------------------------------------------------------

def assemble_and_verify(target_l: dict[int, sp.Expr], natural_seed: dict[int, sp.Expr], aligned_seed: dict[int, sp.Expr], exprs: dict[str, sp.Expr], residual_blocks: dict[str, sp.Expr]) -> None:
    banner("PART V — ASSEMBLE THE ALIGNED LOCAL ORDINARY SEED AND THE FULL LOCAL CANDIDATE")

    # Natural generic-frame local seed (interaction blocks only).
    Qnat = sp.Rational(75, 128) * (a**4 + b**4)
    Tnat = sp.Rational(59, 16) * (q * a**3 + p * b**3)
    Snat = sp.Rational(203, 32) * (q**2 * a**2 + p**2 * b**2)
    Unat = sp.Rational(31, 32) * (q**3 * a + p**3 * b)
    Wnat = sp.Rational(1, 16) * (q**4 + p**4)

    seed_aligned_blocks = {
        "Q": sp.expand(Qnat + exprs["Q"]),
        "T": sp.expand(Tnat + exprs["T"]),
        "S": sp.expand(Snat + exprs["S"]),
        "U": sp.expand(Unat + exprs["U"]),
        "W": sp.expand(Wnat),
    }

    subbanner("V.1 — Natural seed direct COM check")
    for block, expr, slots in [
        ("Q", Qnat, [natural_seed[i] for i in range(7, 12)]),
        ("T", Tnat, [natural_seed[i] for i in range(12, 16)]),
        ("S", Snat, [natural_seed[i] for i in range(16, 19)]),
        ("U", Unat, [natural_seed[i] for i in range(19, 21)]),
        ("W", Wnat, [natural_seed[21]]),
    ]:
        lhs = mod_can["block_slots"](mod_can["to_nu"](expr), block)
        for i, (x, y) in enumerate(zip(lhs, slots), start=1):
            expect_zero(f"natural {block} slot {i}", x - y)

    subbanner("V.2 — Aligned seed direct COM check")
    for block, expr, slots in [
        ("Q", seed_aligned_blocks["Q"], [aligned_seed[i] for i in range(7, 12)]),
        ("T", seed_aligned_blocks["T"], [aligned_seed[i] for i in range(12, 16)]),
        ("S", seed_aligned_blocks["S"], [aligned_seed[i] for i in range(16, 19)]),
        ("U", seed_aligned_blocks["U"], [aligned_seed[i] for i in range(19, 21)]),
        ("W", seed_aligned_blocks["W"], [aligned_seed[21]]),
    ]:
        lhs = mod_can["block_slots"](mod_can["to_nu"](expr), block)
        for i, (x, y) in enumerate(zip(lhs, slots), start=1):
            expect_zero(f"aligned {block} slot {i}", x - y)

    subbanner("V.3 — One full local ordinary generic-frame candidate")
    full_blocks = {
        "Q": sp.expand(seed_aligned_blocks["Q"] + residual_blocks["Q"]),
        "T": sp.expand(seed_aligned_blocks["T"] + residual_blocks["T"]),
        "S": sp.expand(seed_aligned_blocks["S"] + residual_blocks["S"]),
        "U": sp.expand(seed_aligned_blocks["U"] + residual_blocks["U"]),
        "W": sp.expand(seed_aligned_blocks["W"] + residual_blocks["W"]),
    }

    for block, expr, slots in [
        ("Q", full_blocks["Q"], [target_l[i] for i in range(7, 12)]),
        ("T", full_blocks["T"], [target_l[i] for i in range(12, 16)]),
        ("S", full_blocks["S"], [target_l[i] for i in range(16, 19)]),
        ("U", full_blocks["U"], [target_l[i] for i in range(19, 21)]),
        ("W", full_blocks["W"], [target_l[21]]),
    ]:
        lhs = mod_can["block_slots"](mod_can["to_nu"](expr), block)
        for i, (x, y) in enumerate(zip(lhs, slots), start=1):
            expect_zero(f"full local {block} slot {i}", x - y)

    print("\nSo one exact local generic-frame ordinary candidate is")
    print("\n    L4,loc^(gen) = L4,seed^(natural,gen) + delta_seed^(gen) + ΔL4,loc^(can),\n")
    print("with delta_seed^(gen) given by the structured lifts above and ΔL4,loc^(can)")
    print("the already solved canonical comparable-mass residual from the previous audit.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    target_l, natural_seed, aligned_seed, mis_blocks, residual_blocks = rebuild_seed_alignment_data()
    delta_seed = {i: sp.simplify(aligned_seed[i] - natural_seed[i]) for i in range(1, 22)}

    print_com_correction_ledger(delta_seed, mis_blocks)
    old_bases, _ = old_scaffold_failure(mis_blocks)
    seed_bases, coords, exprs = structured_seed_spaces(mis_blocks, old_bases)
    print_canonical_representative(seed_bases, coords, exprs)
    assemble_and_verify(target_l, natural_seed, aligned_seed, exprs, residual_blocks)

    banner("PART VI — FINAL LEDGER")
    print("1. The local ordinary seed-alignment correction is a real generic-frame object, not a")
    print("   failure of the Hamiltonian-derived comparable-mass residual lift.")
    print("2. The old constant-coefficient comparable-mass interaction scaffold lifts Q_delta and")
    print("   trivially handles W_delta = 0, but it does not lift T_delta, S_delta, or U_delta.")
    print("3. The correction becomes exactly liftable once the seed sector is enlarged in the")
    print("   structured way")
    print("        Q,   T ⊕ (pq)T,   S ⊕ (pq)S,   U ⊕ (p^2 q^2)U,   W=0.")
    print("4. With that seed-sector lift in hand, the natural generic-frame local seed plus")
    print("   delta_seed^(gen) plus the already solved canonical ordinary residual reproduces the")
    print("   exact fixed-chart local 4PN ordinary target block by block after COM reduction.")
    print("5. So the local-first 4PN program now has one exact generic-frame ordinary candidate for")
    print("   the whole local interaction sector.  The remaining open issue is no longer local")
    print("   existence; it is the sharper uniqueness/interpretation question for this aligned-seed")
    print("   sector, with the nonlocal tail bridge still separate.")


if __name__ == "__main__":
    main()
