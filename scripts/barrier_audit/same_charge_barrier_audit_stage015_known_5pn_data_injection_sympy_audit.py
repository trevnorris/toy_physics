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


def fmt(name: str, expr) -> None:
    print(f"{name} = {sp.N(expr, 20)}")


def expect_true(name: str, condition: bool) -> None:
    print(f"{name} = {condition}")
    if not condition:
        raise AssertionError(f"{name} failed")


banner("SAME-CHARGE BARRIER AUDIT — STAGE 015")
print("Known 5PN support/source data injection into the same-charge audit chain")

subbanner("I. Exact numerical inputs from the current 5PN stack")
Lambda_ell = sp.Float("36.94973154240256")
kappa = sp.Float("2457.5087899001137")
zeta_max = sp.Float("2.4675297457259358")

Pe_chi = sp.Float("11155.7265863205869")
zeta_phys_chi = sp.Float("2.4675296478814376")
rho_alpha_chi = sp.Float("3.4675296478814376")

Pe_J = sp.Float("2504.9703142859238")
zeta_phys_J = sp.Float("2.4675278051675084")
rho_alpha_J = sp.Float("3.4675278051675084")

zeta_req = sp.Rational(1, 3)
rho_alpha_req = sp.Rational(4, 3)

fmt("Lambda_ell", Lambda_ell)
fmt("kappa", kappa)
fmt("zeta_max", zeta_max)
fmt("Pe_chi", Pe_chi)
fmt("zeta_phys_chi", zeta_phys_chi)
fmt("rho_alpha_chi", rho_alpha_chi)
fmt("Pe_J", Pe_J)
fmt("zeta_phys_J", zeta_phys_J)
fmt("rho_alpha_J", rho_alpha_J)
fmt("zeta_req", zeta_req)
fmt("rho_alpha_req", rho_alpha_req)

subbanner("II. Exact support/source margins")
zeta_margin_chi = sp.simplify(zeta_phys_chi - zeta_req)
rho_margin_chi = sp.simplify(rho_alpha_chi - rho_alpha_req)

zeta_margin_J = sp.simplify(zeta_phys_J - zeta_req)
rho_margin_J = sp.simplify(rho_alpha_J - rho_alpha_req)

fmt("zeta_margin_chi", zeta_margin_chi)
fmt("rho_margin_chi", rho_margin_chi)
fmt("zeta_margin_J", zeta_margin_J)
fmt("rho_margin_J", rho_margin_J)

subbanner("III. Ratio form")
zeta_ratio_chi = sp.simplify(zeta_phys_chi / zeta_req)
rho_ratio_chi = sp.simplify(rho_alpha_chi / rho_alpha_req)

zeta_ratio_J = sp.simplify(zeta_phys_J / zeta_req)
rho_ratio_J = sp.simplify(rho_alpha_J / rho_alpha_req)

fmt("zeta_ratio_chi", zeta_ratio_chi)
fmt("rho_ratio_chi", rho_ratio_chi)
fmt("zeta_ratio_J", zeta_ratio_J)
fmt("rho_ratio_J", rho_ratio_J)

subbanner("IV. Distance to the Family-1 ceiling")
ceiling_headroom_chi = sp.simplify(zeta_max - zeta_phys_chi)
ceiling_headroom_J = sp.simplify(zeta_max - zeta_phys_J)
fmt("ceiling_headroom_chi", ceiling_headroom_chi)
fmt("ceiling_headroom_J", ceiling_headroom_J)

subbanner("V. Basic theorem-gate booleans")
expect_true("chi branch clears zeta requirement", bool(zeta_phys_chi > zeta_req))
expect_true("chi branch clears rho_alpha requirement", bool(rho_alpha_chi > rho_alpha_req))
expect_true("J branch clears zeta requirement", bool(zeta_phys_J > zeta_req))
expect_true("J branch clears rho_alpha requirement", bool(rho_alpha_J > rho_alpha_req))

expect_true("chi branch lies below Family-1 ceiling", bool(zeta_phys_chi < zeta_max))
expect_true("J branch lies below Family-1 ceiling", bool(zeta_phys_J < zeta_max))

subbanner("VI. Transported same-charge verdict compiler")
# Carried from the earlier same-charge audit chain:
B_stat_both = sp.Float("0.367930328492646")
B_stat_nonempty = sp.Float("0.737619063660757")
B_dyn_both = sp.Float("0.967282389363822")
B_dyn_nonempty = sp.Float("0.990581810705233")

fmt("B_stat_both", B_stat_both)
fmt("B_stat_nonempty", B_stat_nonempty)
fmt("B_dyn_both", B_dyn_both)
fmt("B_dyn_nonempty", B_dyn_nonempty)

print("\nLogical consequences:")
print("1. The dynamic window is looser than the static window (carried theorem from stages 013–014).")
print("2. The support/source side is numerically safe on both explicit Family-1 branches.")
print("3. Therefore the first unresolved kill condition is still the actual PDE-selected")
print("   orbit-lock / coherent placement point, not the support/source branch itself.")

orbit_lock_point_numerically_located = False
print(f"orbit_lock_point_numerically_located = {orbit_lock_point_numerically_located}")

print("\nFinal verdict:")
print("same_charge_route_status = ALIVE_BUT_UNRESOLVED")
print("resolved_bottlenecks = {dynamic_window, support_source_side}")
print("remaining_unresolved_object = actual_PDE_selected_orbit_lock_point")
