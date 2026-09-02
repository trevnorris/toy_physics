"""PY-engine measurement: jet orders in the ENERGY DENSITY, the u-FLUX
dL/d(d_i u_a), the local dL/du_a, and the unconstrained bulk EL row.

Prints operands; states no conclusion.  Case = LAB_HELD / RHO4_CONSTANT /
EULERIAN.  Builds construct_energy at background_depth=3 so that if the
density already carried a Hessian, the EL divergence would be allowed to
raise it to order 3.
"""
from __future__ import annotations

import sys
import importlib.util

import sympy as sp

ENGINE = "/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py"
spec = importlib.util.spec_from_file_location("s11cb_engine", ENGINE)
eng = importlib.util.module_from_spec(spec)
sys.modules["s11cb_engine"] = eng
spec.loader.exec_module(eng)

BRANCH, REP, ROUTE = "LAB_HELD", "RHO4_CONSTANT", "EULERIAN"
DEPTH = 3

order_atoms = {}
for order in range(0, 5):
    atoms = set()
    if order == 0:
        atoms.add(eng.W_bg)
        atoms.add(eng.mu_R_bg)
    for src in ("W_BG", "MU_R_BG"):
        table = eng.BACKGROUND_PROFILE_JETS.get(src, {})
        if order in table:
            atoms |= set(table[order].values())
        if order == 1:
            atoms |= set(eng.BACKGROUND_FIRST_JETS.get(src, ()))
    order_atoms[order] = atoms


def orders_in(expr):
    fs = sp.expand(expr).free_symbols
    return sorted(o for o, atoms in order_atoms.items() if o >= 1 and fs & atoms)


def atoms_of_order(expr, order):
    fs = sp.expand(expr).free_symbols
    return sorted(str(a) for a in (fs & order_atoms.get(order, set())))


def sigma_powers(expr):
    powers = set()
    for term in sp.Add.make_args(sp.expand(expr)):
        p = term.as_powers_dict().get(eng.sigma_W, 0)
        if p != 0:
            powers.add(int(p) if p == int(p) else p)
    return sorted(powers)


print("Building construct_energy(LAB_HELD, background_depth=3) ...", flush=True)
build = eng.construct_energy(BRANCH, background_depth=DEPTH)
density = build.density
print(f"ENERGY_TERM_COUNT: {len(build.terms)}", flush=True)
print(f"ENERGY_COUNT_FIELD: {build.count}", flush=True)

print("\n=== ENERGY_DENSITY ===")
print(f"DENSITY_BG_JET_ORDERS: {orders_in(density)}")
print(f"DENSITY_ORDER>=2_ATOMS ({len(atoms_of_order(density, 2))}): {atoms_of_order(density, 2)[:20]}")
print(f"DENSITY_ORDER>=3_ATOMS ({len(atoms_of_order(density, 3))}): {atoms_of_order(density, 3)}")
print(f"DENSITY_SIGMA_W_POWERS_PRESENT: {sigma_powers(density)}")

# Per-term density orders (identify any term that already has Hessian).
print("\n=== PER_TERM_DENSITY_JET_ORDERS ===")
hessian_terms = []
for label, term in build.terms:
    o = orders_in(term)
    o3 = atoms_of_order(term, 3)
    o2 = atoms_of_order(term, 2)
    if o2 or o3 or (o and max(o) >= 2):
        hessian_terms.append((label, o, len(o2), len(o3)))
        print(f"TERM {label} orders={o} n_hess={len(o2)} n_o3={len(o3)}")
if not hessian_terms:
    print("NO_TERM_HAS_ORDER_GE_2")

# Flux dL/d(grad u) and local dL/du
print("\n=== U_FLUX_AND_LOCAL ===")
flux_orders = set()
local_orders = set()
flux_o3 = []
local_o3 = []
flux_o2_count = 0
for a in eng.DIRECTIONS:
    loc = sp.diff(density, eng.u[a])
    lo = orders_in(loc)
    local_orders.update(lo)
    local_o3.extend(atoms_of_order(loc, 3))
    print(f"LOCAL_dL_du[{a}]_ORDERS: {lo}")
    for i in eng.DIRECTIONS:
        fl = sp.diff(density, eng.grad_u[a][i])
        fo = orders_in(fl)
        flux_orders.update(fo)
        flux_o3.extend(atoms_of_order(fl, 3))
        flux_o2_count += len(atoms_of_order(fl, 2))
        print(f"FLUX_dL_d(du)[{a},{i}]_ORDERS: {fo}")
print(f"FLUX_MAX_ORDER_P: {max(flux_orders) if flux_orders else 0}")
print(f"FLUX_ALL_ORDERS: {sorted(flux_orders)}")
print(f"FLUX_N_HESSIAN_ATOM_HITS: {flux_o2_count}")
print(f"FLUX_ORDER3_ATOMS: {sorted(set(flux_o3))}")
print(f"LOCAL_ALL_ORDERS: {sorted(local_orders)}")
print(f"LOCAL_ORDER3_ATOMS: {sorted(set(local_o3))}")

# Unconstrained EL via the engine's live_strong_rows at depth 3
print("\n=== UNCONSTRAINED_LIVE_STRONG_ROWS depth=3 ===")
rows = eng.live_strong_rows(density, BRANCH, REP, ROUTE, DEPTH)
u_rows = eng.named_tuple_row(eng.casify(rows), "U")
theta_row = eng.named_tuple_row(eng.casify(rows), "THETA")
ew_row = eng.named_tuple_row(eng.casify(rows), "E_W")

def report_row(name, expr):
    if isinstance(expr, (tuple, list, sp.Tuple)):
        ords = set()
        o3 = []
        o2 = []
        for item in expr:
            ords.update(orders_in(item))
            o3.extend(atoms_of_order(item, 3))
            o2.extend(atoms_of_order(item, 2))
        print(f"{name}_BG_JET_ORDERS: {sorted(ords)}")
        print(f"{name}_N_HESSIAN_ATOMS: {len(set(o2))}")
        print(f"{name}_ORDER3_ATOMS ({len(set(o3))}): {sorted(set(o3))}")
        return sorted(ords)
    ords = orders_in(expr)
    print(f"{name}_BG_JET_ORDERS: {ords}")
    print(f"{name}_ORDER3_ATOMS: {atoms_of_order(expr, 3)}")
    return ords

u_ords = report_row("U_EL_UNCONSTRAINED", u_rows)
report_row("THETA_EL", theta_row)
report_row("EW_EL", ew_row)

# Retained-grade U rows
print("\n=== RETAINED_GRADE_U_EL ===")
rg = eng.retained_grade(u_rows)
report_row("RETAINED_U_EL", rg)

# KEY FACT on engine map
print("\n=== ENGINE_BACKGROUND_JET_EXPRESSION_SIGMA_POWER ===")
for order, index in [(1, (0,)), (2, (0, 0)), (3, (0, 0, 0)), (3, (0, 1, 2))]:
    if order == 1:
        expr = eng.background_jet_expression("W_BG", 1, index)
    else:
        expr = eng.background_jet_expression("W_BG", order, index)
    print(f"JET_ORDER_{order}_{index}: {expr}")
    print(f"JET_ORDER_{order}_{index}_SIGMA_POWER: {sp.expand(expr).as_powers_dict().get(eng.sigma_W, 0)}")

print("\nDONE", flush=True)
