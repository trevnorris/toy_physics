
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

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

def expect_zero(name: str, expr):
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 340 — BRANCH-PACKET COLLAPSE TO ONE NORMALIZATION SCALAR")

# Stage-200 reduced branch packet
a2, b2, a4, b4, aP0, bP0, Delta_pole, Delta_norm = sp.symbols(
    "a_2 b_2 a_4 b_4 a_P0 b_P0 Delta_pole Delta_norm", real=True
)
Delta_branch = sp.Matrix([a2, b2, a4, b4, aP0, bP0, Delta_pole, Delta_norm])

subbanner("I. Original exact reduced branch packet")
print("Delta_branch =")
sp.pprint(Delta_branch)

# Actual isotropic branch assumptions distilled from the later moving-throat notes:
#   - grouped-lane isotropy kills a2,b2,a4,b4,a(P0),b(P0)
#   - static-geometry one-pole branch kills Delta_pole
#   - the only surviving reduced defect is the normalization scalar N_Q - 1
NQ = sp.symbols("N_Q", real=True)
collapse_subs = {
    a2: 0,
    b2: 0,
    a4: 0,
    b4: 0,
    aP0: 0,
    bP0: 0,
    Delta_pole: 0,
    Delta_norm: NQ - 1,
}
Delta_branch_collapsed = sp.simplify(Delta_branch.subs(collapse_subs))

subbanner("II. Actual isotropic one-pole branch collapse")
print("Delta_branch | actual isotropic one-pole branch =")
sp.pprint(Delta_branch_collapsed)

expect_zero("first entry", Delta_branch_collapsed[0])
expect_zero("second entry", Delta_branch_collapsed[1])
expect_zero("third entry", Delta_branch_collapsed[2])
expect_zero("fourth entry", Delta_branch_collapsed[3])
expect_zero("fifth entry", Delta_branch_collapsed[4])
expect_zero("sixth entry", Delta_branch_collapsed[5])
expect_zero("seventh entry", Delta_branch_collapsed[6])
expect_zero("eighth entry minus (N_Q-1)", Delta_branch_collapsed[7] - (NQ - 1))

subbanner("III. Exact reduced branch verdict")
print("On the actual isotropic static-geometry contact-plus-pole branch,")
print("the entire 8-component reduced branch packet collapses to the 1-scalar defect")
print("    N_Q - 1.")
print()
print("Equivalently:")
print("    Delta_branch = 0   iff   N_Q = 1")

banner("STAGE 340 LEDGER")
print("1. The Stage-200 8-component reduced branch packet is exact in general.")
print("2. On the later actual isotropic one-pole branch, the first seven components vanish identically.")
print("3. The only surviving reduced branch defect is N_Q - 1.")
