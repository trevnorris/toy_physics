#!/usr/bin/env python3
"""Stage V2-22C lightweight SymPy formula audit.

This audit checks only the orchestration identities for the V2-22C end-to-end
pipeline.  The heavier component audits remain in V2-21, V2-22A, and V2-22B.
"""
from __future__ import annotations

import sympy as sp


def main() -> int:
    checks = []

    x20, x21, x22 = sp.symbols("x20 x21 x22")
    xbar = (x20 + 2*x21 + 2*x22)/5
    ax = (2*x20 - x21 - x22)/10
    bx = (x21 - x22)/2
    checks.append(("group_inverse_20", sp.simplify(xbar + 4*ax - x20) == 0))
    checks.append(("group_inverse_21", sp.simplify(xbar - ax + bx - x21) == 0))
    checks.append(("group_inverse_22", sp.simplify(xbar - ax - bx - x22) == 0))

    eps, x0, x1 = sp.symbols("eps x0 x1")
    y20 = x0 + eps*x1
    y21 = x0 + eps*sp.Rational(1, 2)*x1
    y22 = x0 - eps*x1
    trace = sp.simplify((y20 + 2*y21 + 2*y22)/5)
    ay = sp.simplify((2*y20 - y21 - y22)/10)
    by = sp.simplify((y21 - y22)/2)
    checks.append(("axisym_trace_unchanged", sp.simplify(trace - x0) == 0))
    checks.append(("axisym_b_equals_3a", sp.simplify(by - 3*ay) == 0))

    D0, D2, D4, N0 = sp.symbols("D0 D2 D4 N0", nonzero=True)
    N2_const = 2*D2*N0/D0
    N4_const = N0*(D2**2 + 2*D0*D4)/D0**2
    P2 = (D0*N2_const - 2*D2*N0)/D0**2
    P4 = (D0**2*N4_const - 2*D0*(D2*N2_const + D4*N0) + 3*D2**2*N0)/D0**3
    checks.append(("constant_prefactor_P2_zero", sp.simplify(P2) == 0))
    checks.append(("constant_prefactor_P4_zero", sp.simplify(P4) == 0))

    G, cs, a, c, mhat, Sport = sp.symbols("G cs a c mhat Sport", positive=True)
    Ptarget = 54*G*cs**5/(5*a**5*c**5)
    gamma_eff = mhat**2*Sport*Ptarget*a**5/(27*cs**5)
    gamma_gr = Sport*mhat**2*2*G/(5*c**5)
    checks.append(("target_gamma_equivalence", sp.simplify(gamma_eff - gamma_gr) == 0))

    passed = sum(1 for _, ok in checks if ok)
    print("STAGE V2-22C FORMULA SYMPY AUDIT")
    print(f"checks_passed: {passed}/{len(checks)}")
    for name, ok in checks:
        print(f"{name}: {'PASS' if ok else 'FAIL'}")
    return 0 if passed == len(checks) else 2


if __name__ == "__main__":
    raise SystemExit(main())
