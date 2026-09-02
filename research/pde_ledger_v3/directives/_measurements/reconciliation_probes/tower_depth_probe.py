"""Decisive PY-side reconciliation probe: does forming a strong U-momentum row
generate nonzero ORDER-3 background jets that SURVIVE the retained grade?

PRINTS operands + residuals; states no conclusion (rule 2). The step record interprets.

Compares retained-grade strong rows at background_depth in {1,2,3,4} for one
representative case (LAB_HELD, RHO4_CONSTANT, EULERIAN). Key objects:
  - residual_32 = rows(depth3) - rows(depth2)   <- the WL(3) vs PY-default(2) question
  - residual_43 = rows(depth4) - rows(depth3)   <- does the tower TERMINATE at 3?
  - residual_21 = rows(depth2) - rows(depth1)   <- sanity: the Hessian must matter (live instrument)
For any surviving residual it prints the max sigma_W power (retained grade keeps <=1)
and which background-jet ORDER atoms appear.
"""
import sys
import importlib.util
import sympy as sp

ENGINE = "/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py"
spec = importlib.util.spec_from_file_location("s11cb_engine", ENGINE)
eng = importlib.util.module_from_spec(spec)
sys.modules["s11cb_engine"] = eng  # dataclasses resolve cls.__module__ via sys.modules
spec.loader.exec_module(eng)

BRANCH, REP, ROUTE = "LAB_HELD", "RHO4_CONSTANT", "EULERIAN"
sigma_W = eng.sigma_W

# background-jet atoms by ORDER (order 1 = grad, 2 = Hessian, 3, 4 ...)
order_atoms = {}
for order in range(1, 5):
    atoms = set()
    for src in ("W_BG", "MU_R_BG"):
        table = eng.BACKGROUND_PROFILE_JETS.get(src, {})
        if order in table:
            atoms |= set(table[order].values())
        if order == 1:
            atoms |= set(eng.BACKGROUND_FIRST_JETS.get(src, ()))
    order_atoms[order] = atoms


def rows_at(depth):
    build = eng.construct_energy(BRANCH, background_depth=depth)
    return eng.retained_grade(
        eng.live_strong_rows(build.density, BRANCH, REP, ROUTE, depth)
    )


def leaves(value):
    return [sp.expand(s) for s in eng.object_scalars(value)]


def report_pair(name, hi, lo):
    assert len(hi) == len(lo), f"{name}: operand shape mismatch {len(hi)} vs {len(lo)}"
    diffs = [sp.expand(a - b) for a, b in zip(hi, lo)]
    nz = [d for d in diffs if d != 0]
    print(f"\n=== {name} ===")
    print(f"{name}_NONZERO_LEAF_COUNT: {len(nz)} of {len(diffs)}")
    if not nz:
        print(f"{name}_RESIDUAL: 0  (identical after retained grade)")
        return
    total = sp.Add(*nz)
    maxsig = 0
    present = set()
    for term in sp.Add.make_args(sp.expand(total)):
        p = term.as_powers_dict().get(sigma_W, 0)
        maxsig = max(maxsig, int(p) if p == int(p) else p)
        fs = term.free_symbols
        for order, atoms in order_atoms.items():
            if fs & atoms:
                present.add(order)
    print(f"{name}_MAX_SIGMA_W_POWER: {maxsig}   (retained grade keeps <=1)")
    print(f"{name}_BG_JET_ORDERS_PRESENT: {sorted(present)}")
    sample = nz[0]
    print(f"{name}_FIRST_NONZERO_LEAF (truncated 400 chars):\n  {str(sample)[:400]}")


print("Building retained-grade strong rows at depths 1..4 ...", flush=True)
r1 = leaves(rows_at(1)); print(f"depth1 leaves={len(r1)}", flush=True)
r2 = leaves(rows_at(2)); print(f"depth2 leaves={len(r2)}", flush=True)
r3 = leaves(rows_at(3)); print(f"depth3 leaves={len(r3)}", flush=True)
r4 = leaves(rows_at(4)); print(f"depth4 leaves={len(r4)}", flush=True)

report_pair("RESIDUAL_21_HESSIAN_SANITY", r2, r1)
report_pair("RESIDUAL_32_DECISIVE", r3, r2)
report_pair("RESIDUAL_43_TERMINATION", r4, r3)

# Distinguish "order-3 never generated in strong rows" from "generated then cancels".
# Scan the UN-reduced depth-3 strong rows for order>=3 background-jet atoms.
print("\n=== UNREDUCED DEPTH-3 STRONG ROWS: background-jet orders present ===")
build3 = eng.construct_energy(BRANCH, background_depth=3)
unreduced3 = eng.live_strong_rows(build3.density, BRANCH, REP, ROUTE, 3)
u3_syms = set()
for s in eng.object_scalars(unreduced3):
    u3_syms |= sp.expand(s).free_symbols
orders_unreduced = sorted(o for o, atoms in order_atoms.items() if u3_syms & atoms)
print(f"UNREDUCED_DEPTH3_BG_JET_ORDERS_PRESENT: {orders_unreduced}")
o3 = order_atoms.get(3, set()) | order_atoms.get(4, set())
hit = sorted(str(a) for a in (u3_syms & o3))
print(f"UNREDUCED_DEPTH3_ORDER>=3_ATOMS ({len(hit)}): {hit[:12]}")
print("\nDONE", flush=True)
