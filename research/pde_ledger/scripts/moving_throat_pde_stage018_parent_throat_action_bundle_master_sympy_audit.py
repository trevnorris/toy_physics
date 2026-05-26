#!/usr/bin/env python3
"""Master-note audit for step_16_parent_throat_action_bundle_master_notes.md."""
from __future__ import annotations

import sympy as sp


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


def main() -> None:
    w = sp.symbols("w", real=True)
    beta = sp.exp(-w**2 / 2)
    MSigma_example = sp.integrate(beta**2, (w, -sp.oo, sp.oo))
    KSigma_example = sp.integrate(sp.diff(beta, w)**2 + beta**2, (w, -sp.oo, sp.oo))
    assert_zero("example wall inertia integral", MSigma_example - sp.sqrt(sp.pi))
    assert_zero("example wall stiffness integral", KSigma_example - 3*sp.sqrt(sp.pi)/2)

    print("STEP 16 PARENT THROAT ACTION BUNDLE MASTER AUDIT")
    print("Checked concrete Gaussian parent-wall inertia and stiffness branch integrals.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
