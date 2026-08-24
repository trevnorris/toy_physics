#!/usr/bin/env python3
"""S11c-a T7 step-0 AXIS-SET CENSUS v2 (rule 2: prints computed objects; asserts nothing).

Per joined tag, per engine: fully-flattened LEAF-case count + the set of case AXES used.
PY nests quantity-grouped tags one level (bare Str('QUANT') key); leaf cases are positional
Tuple keys. WL flattens every axis into one pipe-string, so top-level keys are already leaves.
The decisive object is the AXIS-SET DIFF, not the raw count.
"""
import re

PY = "/home/trevnorris/.s11_build/S11c_a_sympy_engine.out"
WL = "/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out"
OUTKEYS = "/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/s11ca_axis_keys.txt"

def load(path, prefix):
    out = {}
    with open(path) as f:
        for line in f:
            m = re.match(r'^(' + prefix + r'_S11CA_[A-Z0-9_]+):\s?(.*)$', line.rstrip('\n'))
            if m:
                out[m.group(1)[len(prefix)+1:]] = m.group(2)
    return out

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
    return [a.strip() for a in args if a.strip() != ""]

PYB, PYC = {'(', '[', '{'}, {')', ']', '}'}

def py_key_tokens(keystr):
    toks = re.findall(r"Str\('([^']*)'\)|Integer\((-?\d+)\)", keystr)
    return [a or b for a, b in toks]

def py_flatten(payload):
    """Return list of leaf keys (each a tuple of tokens incl. any accumulated quantity)."""
    p = payload.strip()
    if not p.startswith("Tuple("):
        return None
    def recurse(s, prefix):
        inner = s[len("Tuple("):-1] if s.startswith("Tuple(") and s.endswith(")") else None
        if inner is None:
            return []
        pairs = split_top(inner, PYB, PYC)
        leaves = []
        for pr in pairs:
            kv = split_top(pr[len("Tuple("):-1], PYB, PYC) if pr.startswith("Tuple(") and pr.endswith(")") else None
            if not kv or len(kv) < 2:
                return None  # non-standard -> signal
            key, val = kv[0].strip(), kv[1].strip()
            if key.startswith("Str(") and not key.startswith("Str('") == False and key.count(",") == 0 and "Integer(" not in key and "Tuple(" not in key:
                # bare Str -> quantity grouping level: descend one level
                q = py_key_tokens(key)
                sub = recurse(val, prefix + q)
                if sub is None: return None
                leaves.extend(sub)
            elif key.startswith("Tuple("):
                # positional leaf case key
                leaves.append(tuple(prefix + py_key_tokens(key)))
            else:
                return None
        return leaves
    return recurse(p, [])

def wl_keys(payload):
    p = payload.strip()
    if not p.startswith("<|"): return None
    inner = p[2:-2] if p.endswith("|>") else p[2:]
    entries = split_top(inner, {'<|', '[', '{', '('}, {'|>', ']', '}', ')'})
    keys = []
    for e in entries:
        m = re.match(r'^"([^"]*)"\s*->', e.strip())
        keys.append(tuple(m.group(1).split("|")) if m else ("RAW",))
    return keys

QUANT = {  # physics quantity stems that can appear as a grouping / leading axis
 'RELATIVE_FLUX','TRACTION','CLOSURE_SHAPE_DERIV','VIRTUAL_CONSTRAINT','KINEMATIC_BALANCE',
 'FACE_NORMAL','CONORMAL_DERIV','FACE_MEASURE_SHAPE_DERIV','FACE_VELOCITY','EVOLUTION_MASS_BALANCE',
 'PROJECTION_SHAPE_DERIV','PROJECTION_STATIC_OPERAND','PROJECTION_DYNAMIC_OPERAND','PROJECTION_RESIDUAL',
 'VIRTUAL_WORK_SHAPE_DERIV','FACE_SHIFT','FACE_MAP'}

def axis_of(tok):
    if tok in ('LAB_HELD','MATERIAL_ADVECTED'): return 'BRANCH'
    if tok in ('RHO4_CONSTANT','RHOBR_CONSTANT'): return 'DENSITY'
    if tok in ('FACE_PLUS','FACE_MINUS','1','-1','BOTH_FACES'): return 'FACE'
    if tok in ('DELTA_W','ZETA_C') or tok.startswith('DOF_'): return 'DOF'
    if tok.startswith('VIRTUAL_DOF_'): return 'VIRTUAL_DOF'
    if tok.startswith('DIRECTION'): return 'DIRECTION'
    if tok.startswith('FIELD'): return 'FIELD'
    if tok.startswith('ORIGIN'): return 'ORIGIN'
    if tok in QUANT: return 'QUANTITY'
    return 'OTHER:'+tok

def axes_of(keys):
    ax = set()
    for k in keys:
        for t in k:
            ax.add(axis_of(t))
    return ax

py = load(PY, "PY"); wl = load(WL, "WL")
joined = sorted(set(py) & set(wl))
fk = open(OUTKEYS, "w")
AX_ORDER = ['QUANTITY','BRANCH','DENSITY','FACE','DOF','VIRTUAL_DOF','DIRECTION','FIELD','ORIGIN']
def fmt(ax):
    known = [a for a in AX_ORDER if a in ax]
    other = sorted(a for a in ax if a.startswith('OTHER'))
    return "+".join([a[:4] for a in known] + other) or "-"

print(f"{'TAG':44s} {'PYn':>4s} {'WLn':>4s}  {'PY_axes':<28s} {'WL_axes':<28s} DIFF")
print("-"*130)
for tag in joined:
    pk = py_flatten(py[tag]); wk = wl_keys(wl[tag])
    if pk is None:
        pn, pax = "BESPOKE", set()
    else:
        pn, pax = len(pk), axes_of(pk)
    wn, wax = (len(wk), axes_of(wk)) if wk else ("?", set())
    diff = ""
    if pk is not None and wk:
        only_py = pax - wax; only_wl = wax - pax
        parts = []
        if only_wl: parts.append("WL+{" + ",".join(sorted(only_wl)) + "}")
        if only_py: parts.append("PY+{" + ",".join(sorted(only_py)) + "}")
        if pn != wn and not parts: parts.append("count %s/%s" % (pn, wn))
        diff = " ".join(parts)
    print(f"{tag:44s} {str(pn):>4s} {str(wn):>4s}  {fmt(pax):<28s} {fmt(wax):<28s} {diff}")
    fk.write(f"\n=== {tag} (PYn={pn} WLn={wn}) PY_axes={sorted(pax)} WL_axes={sorted(wax)} ===\n")
    for k in (pk or [])[:40]: fk.write("  PY " + "|".join(map(str,k)) + "\n")
    for k in (wk or [])[:40]: fk.write("  WL " + "|".join(map(str,k)) + "\n")
fk.close()
print("\n(full leaf keys -> scratchpad/s11ca_axis_keys.txt)")
