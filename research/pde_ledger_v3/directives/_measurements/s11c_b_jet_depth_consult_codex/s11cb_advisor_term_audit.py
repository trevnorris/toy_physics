#!/usr/bin/env python3
"""Term-level audit of raw versus mass-constrained S11c-b bulk U rows."""

from __future__ import annotations

import importlib.util
import sys

import sympy as sp


ENGINE = "/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py"
spec = importlib.util.spec_from_file_location("s11cb_engine_for_advisor", ENGINE)
eng = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = eng
spec.loader.exec_module(eng)


def jet_atoms_by_order(maximum: int = 4) -> dict[int, set[sp.Symbol]]:
    result: dict[int, set[sp.Symbol]] = {order: set() for order in range(maximum + 1)}
    result[0].update((eng.W_bg, eng.mu_R_bg))
    for source in ("W_BG", "MU_R_BG"):
        result[1].update(eng.BACKGROUND_FIRST_JETS[source])
        for order, table in eng.BACKGROUND_PROFILE_JETS[source].items():
            if order <= maximum:
                result[order].update(table.values())
    return result


ORDER_ATOMS = jet_atoms_by_order()


def report(label: str, expressions) -> None:
    expressions = tuple(sp.expand(expression) for expression in expressions)
    symbols = set().union(*(expression.free_symbols for expression in expressions))
    orders = sorted(order for order, atoms in ORDER_ATOMS.items() if symbols & atoms)
    maximum = max(orders, default=0)
    order3 = sorted(str(atom) for atom in symbols & ORDER_ATOMS[3])
    print(f"{label}_JET_ORDERS_PRESENT: {orders}")
    print(f"{label}_MAX_JET_ORDER: {maximum}")
    print(f"{label}_ORDER3_ATOMS ({len(order3)}): {order3}")


build = eng.construct_energy("LAB_HELD", background_depth=3)
label = "W_BG_FIRST_JET_CONTRACTION_06"
term = dict(build.terms)[label]

# This selected energy term is gamma W_i (d_j u_i) (d_j theta).
u_flux = tuple(
    sp.diff(term, eng.grad_u[a][j]) for a in eng.DIRECTIONS for j in eng.DIRECTIONS
)
raw_u = tuple(
    eng.euler_derivative(term, eng.u[a], eng.grad_u[a], background_depth=3)
    for a in eng.DIRECTIONS
)
mu_theta = eng.euler_derivative(
    term, eng.theta, eng.grad_theta, background_depth=3
)

# For EULERIAN/LAB_HELD/RHO4_CONSTANT, the normalized virtual mass constraint is
# delta theta = -div(delta u) - (d_a W/W) delta u_a - (W0/W) delta e.
# Hence the constrained U row is E_u,a + d_a(mu_theta) - (d_a W/W) mu_theta.
constrained_u = tuple(
    raw_u[a]
    + eng.total_derivative(mu_theta, a, background_depth=3)
    - eng.grad_W[a] * mu_theta / eng.W_bg
    for a in eng.DIRECTIONS
)

print(f"SELECTED_TERM_LABEL: {label}")
print("SELECTED_TERM_INDEXED_FORM: gamma * W_i * (d_j u_i) * (d_j theta)")
report("ENERGY_DENSITY", (term,))
report("DU_COEFFICIENT", u_flux)
report("RAW_U_EL_UNREDUCED", raw_u)
report("THETA_EL_UNREDUCED", (mu_theta,))
report("CONSTRAINED_U_UNREDUCED", constrained_u)

raw_retained = tuple(eng.retained_grade(expression) for expression in raw_u)
constrained_retained = tuple(
    eng.retained_grade(expression) for expression in constrained_u
)
report("RAW_U_EL_RETAINED", raw_retained)
report("CONSTRAINED_U_RETAINED", constrained_retained)

# Print the exact coefficient of one decisive third jet after retained reduction.
third_atom = eng.BACKGROUND_PROFILE_JETS["W_BG"][3][(0, 0, 0)]
coefficient = sp.expand(constrained_retained[0]).coeff(third_atom)
print(f"DECISIVE_ATOM: {third_atom}")
print(f"DECISIVE_ATOM_COEFFICIENT: {coefficient}")
print(f"DECISIVE_COEFFICIENT_NONZERO: {coefficient != 0}")

# Global check: no selected energy term or dL/d(d_i u_a) contains a Hessian.
all_terms = tuple(expression for _, expression in build.terms)
all_u_flux = tuple(
    sp.diff(build.density, eng.grad_u[a][j])
    for a in eng.DIRECTIONS
    for j in eng.DIRECTIONS
)
report("ALL_SELECTED_ENERGY_TERMS", all_terms)
report("ALL_DU_COEFFICIENTS", all_u_flux)
print("DONE")
