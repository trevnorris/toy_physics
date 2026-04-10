
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.factor(sp.together(sp.expand(expr))))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

"""
5pn_stage104_family1_core_geometry_map.py

Stage 104 audit: explicit Family-1 geometric fixing of the normalized
hybridization parameter.
"""

banner("STAGE 104 — FAMILY-1 CORE GEOMETRY MAP")

a = sp.symbols("a", positive=True, real=True)
r = sp.symbols("mathfrak_r", positive=True, real=True)

Lambda_star = sp.Rational(37, 20)  # carried Family-1 reference geometry L/a
L = sp.simplify(Lambda_star * a)

# Compensation-selected D/N mixed-tube length from Stage 102
L_W = sp.simplify(sp.pi * a * sp.sqrt((1 + r**2) / 3) / 2)

# Family-1 closure: identify the auxiliary D/N tube length with the actual throat span.
r_F1_solutions = sp.solve(sp.Eq(L_W, L), r)
r_F1 = sp.simplify([sol for sol in r_F1_solutions if sol.is_real][0])

print("L/a =", Lambda_star)
print("L_W/a =")
sp.pprint(sp.simplify(L_W / a))
print("mathfrak_r_F1 =")
sp.pprint(r_F1)
print("mathfrak_r_F1 (numeric) =", sp.N(r_F1, 20))

expect_zero(
    "Family-1 closure identity",
    L_W.subs(r, r_F1) - L
)

banner("STAGE 104 FINAL LEDGER")
print("The explicit Family-1 geometry fixes the normalized hybridization to")
print("  mathfrak_r_F1 = sqrt(12/pi^2 * (37/20)^2 - 1)")
print("               ≈", sp.N(r_F1, 15))
print("This is the unique positive hybridization selected by")
print("  L_W = L = (37/20) a")
print("on the first D/N mixed-tube branch.")
