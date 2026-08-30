#!/usr/bin/env python3
"""Adjudication reporter: parse the committed row-residual instrument's stdout and,
per case, print family/key + whether ROW_RESIDUAL / IN_SCOPE_WEAK_REMAINDER are zero,
their simplified form, and the coupling FULL_PREBRIDGE_ROUTE. Reporting only — no physics."""
import sys
from pathlib import Path
SD = Path("/var/projects/toy_physics/research/pde_ledger_v3/scripts")
sys.path.insert(0, str(SD))
import sympy as sp
import S11c_b_cross_engine_comparator as C  # has TextAtom, Association, serialise

# eval namespace: sympy names + the comparator container classes
NS = {name: getattr(sp, name) for name in dir(sp)}
NS.update({"TextAtom": C.TextAtom, "Association": C.Association})

def ev(rhs: str):
    return eval(rhs, {"__builtins__": {}}, NS)

def get(assoc, key):
    if isinstance(assoc, C.Association):
        return assoc.entries.get(key)
    return None

def key_str(case):
    k = get(case, "KEY")
    fam = get(case, "FAMILY")
    parts = []
    if isinstance(k, C.Association):
        parts = [f"{kk}={getattr(vv,'value',vv)}" for kk, vv in k.entries.items()]
    fam = getattr(fam, "value", fam)
    return f"{fam} [{', '.join(parts)}]"

def zero_and_form(val):
    """val is an Association with VALUE; return (is_zero, short_form, multigrade)."""
    v = get(val, "VALUE") if isinstance(val, C.Association) else val
    mg = get(val, "ETA_SIGMAW_MULTIGRADE_SUPPORT") if isinstance(val, C.Association) else None
    if not isinstance(v, sp.Basic):
        return (v == 0, repr(v), mg)
    s = sp.simplify(v)
    return (s == 0, sp.sstr(s)[:240], mg)

def main(path):
    lines = Path(path).read_text().splitlines()
    cur = None
    rows = []
    pending = {}
    for ln in lines:
        if ln.startswith("CASE = "):
            if cur is not None:
                rows.append((cur, pending))
            cur = ev(ln[len("CASE = "):]); pending = {}
        else:
            for tag in ("ROW_RESIDUAL", "IN_SCOPE_WEAK_REMAINDER", "FULL_PREBRIDGE_ROUTE",
                        "NO_CLEAN_QUOTIENT", "ROW_OPERAND_WL", "ROW_OPERAND_PY_TRUNC"):
                if ln.startswith(tag + " = "):
                    try:
                        pending[tag] = ev(ln[len(tag) + 3:])
                    except Exception as e:
                        pending[tag] = f"<parse-error {e}>"
    if cur is not None:
        rows.append((cur, pending))

    for case, p in rows:
        ks = key_str(case)
        eq = getattr(get(case, "EQUIVALENCE"), "value", "?")
        out = [f"\n=== {ks}  [{eq}]"]
        if "ROW_RESIDUAL" in p:
            z, form, mg = zero_and_form(p["ROW_RESIDUAL"])
            mgs = f"  grades={C.serialise(mg)}" if mg is not None else ""
            out.append(f"  ROW_RESIDUAL {'== 0' if z else 'NONZERO'}{mgs}: {form}")
        if "FULL_PREBRIDGE_ROUTE" in p:
            r = p["FULL_PREBRIDGE_ROUTE"]
            out.append(f"  FULL_PREBRIDGE_ROUTE = {getattr(get(r,'VALUE') if isinstance(r,C.Association) else r,'value', r)}")
        if "IN_SCOPE_WEAK_REMAINDER" in p:
            z, form, mg = zero_and_form(p["IN_SCOPE_WEAK_REMAINDER"])
            out.append(f"  IN_SCOPE_WEAK_REMAINDER {'== 0' if z else 'NONZERO'}: {form}")
        if "NO_CLEAN_QUOTIENT" in p:
            r = p["NO_CLEAN_QUOTIENT"]
            out.append(f"  NO_CLEAN_QUOTIENT = {getattr(get(r,'VALUE') if isinstance(r,C.Association) else r,'value', r)}")
        print("\n".join(out))

if __name__ == "__main__":
    main(sys.argv[1])
