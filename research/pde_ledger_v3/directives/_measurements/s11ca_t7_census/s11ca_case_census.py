#!/usr/bin/env python3
"""S11c-a T7 step-0 CASE-STRUCTURE CENSUS (rule 2: prints computed objects, states no conclusion).

For every joined tag it emits, per engine, the top-level CASE-KEY set and count + an axis
signature. It does NOT decide which engine is correct (that is step 1) and asserts nothing.

PY payload  : sympy-srepr   Tuple( Tuple(<key>, <payload>), ... ), <key>=Tuple(Str('..'),Integer(n),..)
WL payload  : Association    <| "a|b|c" -> <| "FIELD" -> .. |> , .. |>
"""
import sys, re

PY = "/home/trevnorris/.s11_build/S11c_a_sympy_engine.out"
WL = "/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out"

def load(path, prefix):
    out = {}
    with open(path) as f:
        for line in f:
            m = re.match(r'^(' + prefix + r'_S11CA_[A-Z0-9_]+):\s?(.*)$', line.rstrip('\n'))
            if m:
                out[m.group(1)[len(prefix)+1:]] = m.group(2)   # strip "PY_"/"WL_"
    return out

# ---- generic bracket-depth top-level splitter (skips string interiors) ----
def split_top(s, opens, closes, sep=','):
    """Split s on sep at bracket-depth 0. opens/closes are dicts of multichar tokens."""
    args, depth, i, start, n = [], 0, 0, 0, len(s)
    instr = None
    while i < n:
        c = s[i]
        if instr:
            if c == instr:
                instr = None
            i += 1; continue
        if c in ("'", '"'):
            instr = c; i += 1; continue
        # multichar tokens first (WL <| |>)
        two = s[i:i+2]
        if two in opens:
            depth += 1; i += 2; continue
        if two in closes:
            depth -= 1; i += 2; continue
        if c in opens:
            depth += 1; i += 1; continue
        if c in closes:
            depth -= 1; i += 1; continue
        if c == sep and depth == 0:
            args.append(s[start:i]); start = i+1; i += 1; continue
        i += 1
    args.append(s[start:])
    return [a.strip() for a in args if a.strip() != ""]

# ---------- PY ----------
def py_cases(payload):
    """Return list of normalized case-key strings for a PY tag payload, or None if non-standard."""
    p = payload.strip()
    if not p.startswith("Tuple("):
        return None, "non-Tuple-top"
    inner = p[len("Tuple("):]
    if inner.endswith(")"):
        inner = inner[:-1]
    pairs = split_top(inner, {'(', '[', '{'}, {')', ']', '}'})
    keys = []
    for pr in pairs:
        pr = pr.strip()
        if not pr.startswith("Tuple("):
            return None, "pair-not-Tuple"
        pin = pr[len("Tuple("):]
        if pin.endswith(")"):
            pin = pin[:-1]
        kv = split_top(pin, {'(', '[', '{'}, {')', ']', '}'})
        if len(kv) < 2:
            return None, "pair-arity<2"
        keystr = kv[0].strip()
        # key = Tuple(Str('A'), Integer(1), Str('B'))  OR sometimes a bare Str/other
        toks = re.findall(r"Str\('([^']*)'\)|Integer\((-?\d+)\)|Symbol\('([^']*)'", keystr)
        norm = []
        for a, b, c in toks:
            norm.append(a or b or c)
        keys.append("|".join(norm) if norm else keystr[:40])
    return keys, "ok"

# ---------- WL ----------
def wl_cases(payload):
    p = payload.strip()
    if not p.startswith("<|"):
        return None, "non-Assoc-top"
    inner = p[2:]
    if inner.endswith("|>"):
        inner = inner[:-2]
    # top-level entries split on comma at depth 0
    entries = split_top(inner, {'<|', '[', '{', '('}, {'|>', ']', '}', ')'})
    keys = []
    for e in entries:
        e = e.strip()
        m = re.match(r'^"([^"]*)"\s*->', e)
        if not m:
            # key may be an integer or symbol (e.g. <|1 -> ..|>)
            m2 = re.match(r'^([^\->]+?)\s*->', e)
            keys.append(("RAW:" + m2.group(1).strip()) if m2 else ("?" + e[:30]))
            continue
        keys.append(m.group(1))
    return keys, "ok"

def axis_sig_py(keys):
    if not keys: return "-"
    arities = sorted(set(k.count("|")+1 for k in keys))
    return "arity=" + ",".join(map(str, arities))

def axis_sig_wl(keys):
    if not keys: return "-"
    arities = sorted(set(k.count("|")+1 for k in keys if not k.startswith("RAW:")))
    prefixes = set()
    for k in keys:
        for f in k.split("|"):
            mm = re.match(r'^([A-Z]+[A-Z_]*?)(?:_[A-Z0-9]+.*)?$', f)
            # capture leading axis word up to first meaningful underscore-group
            prefixes.add(f.split("_")[0] if "_" not in f else "_".join(f.split("_")[:2]))
    return "arity=" + ",".join(map(str, arities))

py = load(PY, "PY")
wl = load(WL, "WL")
joined = sorted(set(py) & set(wl))

sidecar = open("/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/s11ca_census_keys.txt", "w")
print(f"{'TAG':44s} {'PY_N':>5s} {'WL_N':>5s}  {'PY_sig':<10s} {'WL_sig':<10s}  DIVERGE?")
print("-"*100)
for tag in joined:
    pk, ps = py_cases(py[tag])
    wk, wsx = wl_cases(wl[tag])
    pn = len(pk) if pk is not None else -1
    wn = len(wk) if wk is not None else -1
    psig = axis_sig_py(pk) if pk else ps
    wsig = axis_sig_wl(wk) if wk else wsx
    div = "***" if (pn != wn) else ""
    print(f"{tag:44s} {pn:5d} {wn:5d}  {psig:<10s} {wsig:<10s}  {div}")
    sidecar.write(f"\n=== {tag}  (PY_N={pn}, WL_N={wn}) ===\n")
    sidecar.write("  PY keys:\n")
    for k in (pk or ["<%s>" % ps]): sidecar.write("    " + k + "\n")
    sidecar.write("  WL keys:\n")
    for k in (wk or ["<%s>" % wsx]): sidecar.write("    " + k + "\n")
sidecar.close()
print("\n(full key sets -> scratchpad/s11ca_census_keys.txt)")
