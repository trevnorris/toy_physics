#!/usr/bin/env python3
"""S11c-a T7 step-1 DENSITY-AXIS decidable test (rule 2: prints computed objects; asserts nothing).

For each law where ONE engine emits BOTH density representatives, pair the two rep cases by their
(branch, face, dof) coordinates and compare the stored expression. DIFFER ⇒ the law's shape-derivative
depends on the density representative (the engine that omits the axis is MISSING cases). IDENTICAL ⇒ the
axis is redundant for that law. ⛔ Does not decide the normative requirement (that is the §4/§5 read).

  WL emits both reps for KINEMATIC_BALANCE, RELATIVE_FLUX  -> test there.
  PY emits both reps for the PROJECTION cluster            -> test there.
"""
import re

PY = "/home/trevnorris/.s11_build/S11c_a_sympy_engine.out"
WL = "/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out"

def load(path, prefix, tag):
    for line in open(path):
        if line.startswith(f"{prefix}_{tag}:"):
            return line.split(": ", 1)[1].rstrip("\n")
    return None

def split_top(s, opens, closes, sep=','):
    args, depth, i, start, n, instr = [], 0, 0, 0, len(s), None
    while i < n:
        c = s[i]
        if instr:
            if c == instr: instr = None
            i += 1; continue
        if c in ("'", '"'): instr = c; i += 1; continue
        two = s[i:i+2]
        if two in opens: depth += 1; i += 2; continue
        if two in closes: depth -= 1; i += 2; continue
        if c in opens: depth += 1; i += 1; continue
        if c in closes: depth -= 1; i += 1; continue
        if c == sep and depth == 0: args.append(s[start:i]); start = i+1; i += 1; continue
        i += 1
    args.append(s[start:])
    return [a.strip() for a in args if a.strip()]

WLB, WLC = {'<|', '[', '{', '('}, {'|>', ']', '}', ')'}

def wl_entries(payload):
    p = payload.strip()
    inner = p[2:-2] if p.startswith("<|") and p.endswith("|>") else p[2:]
    out = []
    for e in split_top(inner, WLB, WLC):
        m = re.match(r'^"([^"]*)"\s*->\s*(.*)$', e.strip(), re.S)
        if m: out.append((m.group(1), m.group(2)))
    return out

def wl_field(value, field):
    """Extract a named field's value from a WL <|...|> value blob, else the whole blob."""
    for k, v in wl_entries(value):
        if k == field: return v
    return value

def wl_density_test(tag):
    entries = wl_entries(load(WL, "WL_S11CA", tag))
    groups = {}
    for key, val in entries:
        parts = key.split("|")
        dens = [p for p in parts if p in ("RHO4_CONSTANT", "RHOBR_CONSTANT")]
        if not dens: continue
        coord = tuple(p for p in parts if p not in ("RHO4_CONSTANT", "RHOBR_CONSTANT"))
        groups.setdefault(coord, {})[dens[0]] = wl_field(val, "EXPRESSION")
    print(f"--- WL {tag} ---")
    ndiff = nsame = 0
    for coord, d in sorted(groups.items()):
        if "RHO4_CONSTANT" in d and "RHOBR_CONSTANT" in d:
            same = d["RHO4_CONSTANT"] == d["RHOBR_CONSTANT"]
            nsame += same; ndiff += (not same)
            if not same:
                a, b = d["RHO4_CONSTANT"], d["RHOBR_CONSTANT"]
                # first differing char window
                j = next((k for k in range(min(len(a), len(b))) if a[k] != b[k]), min(len(a), len(b)))
                print(f"  {'|'.join(coord):40s} DIFFER")
    print(f"  => groups: {len(groups)}, DIFFER={ndiff}, IDENTICAL={nsame}")

def py_entries_top(payload):
    p = payload.strip()
    inner = p[len("Tuple("):-1]
    return split_top(inner, {'(', '[', '{'}, {')', ']', '}'})

def py_key_tokens(keystr):
    return [a or b for a, b in re.findall(r"Str\('([^']*)'\)|Integer\((-?\d+)\)", keystr)]

def py_density_test(tag):
    pairs = py_entries_top(load(PY, "PY_S11CA", tag))
    groups = {}
    for pr in pairs:
        kv = split_top(pr[len("Tuple("):-1], {'(', '[', '{'}, {')', ']', '}'})
        key, val = kv[0], kv[1]
        toks = py_key_tokens(key)
        dens = [t for t in toks if t in ("RHO4_CONSTANT", "RHOBR_CONSTANT")]
        if not dens: continue
        coord = tuple(t for t in toks if t not in ("RHO4_CONSTANT", "RHOBR_CONSTANT"))
        groups.setdefault(coord, {})[dens[0]] = val
    print(f"--- PY {tag} ---")
    ndiff = nsame = 0
    for coord, d in sorted(groups.items()):
        if "RHO4_CONSTANT" in d and "RHOBR_CONSTANT" in d:
            same = d["RHO4_CONSTANT"] == d["RHOBR_CONSTANT"]
            nsame += same; ndiff += (not same)
            if not same:
                print(f"  {'|'.join(coord):24s} DIFFER")
    print(f"  => groups: {len(groups)}, DIFFER={ndiff}, IDENTICAL={nsame}")

for t in ("KINEMATIC_BALANCE", "RELATIVE_FLUX"):
    wl_density_test(t)
for t in ("PROJECTION_SHAPE_DERIV", "PROJECTION_STATIC_OPERAND", "PROJECTION_DYNAMIC_OPERAND", "PROJECTION_RESIDUAL"):
    py_density_test(t)
