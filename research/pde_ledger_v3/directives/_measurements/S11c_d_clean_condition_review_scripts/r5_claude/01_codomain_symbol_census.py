#!/usr/bin/env python3
"""Text-level census of the free Symbol names in the committed S11c-b exports
`slab_operator` and `mu_theta_operator` values (no sympy evaluation), and their
registry class in the same exports ledger.  Mechanical: regex over the stored
srepr text.  Reports names absent from the ledger and, per substrate key, which
coordinate families occur."""
import re, sys, collections
P = "/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_exports.py"
lines = open(P, encoding="utf-8").read().split("\n")
# ledger: key -> class (parse 'key': { ... 'class': 'X' })
ledger_class = {}
key = None
for ln in lines:
    m = re.match(r"^    '([^']+)':\s+\{", ln)
    if m:
        key = m.group(1); continue
    m = re.match(r"^        'class': '([A-Z_]+)',", ln)
    if m and key is not None:
        ledger_class[key] = m.group(1)
# symbol display names may differ from ledger keys; build value-name -> class map
value_class = {}
key = None
for ln in lines:
    m = re.match(r"^    '([^']+)':\s+\{", ln)
    if m:
        key = m.group(1); continue
    m = re.match(r"^        'value': _restore\(\"Symbol\('([^']+)'", ln)
    if m and key is not None:
        value_class[m.group(1)] = ledger_class.get(key, "?")
def value_line(name):
    for i, ln in enumerate(lines):
        if ln.startswith(f"    '{name}':"):
            for j in range(i + 1, i + 4):
                if lines[j].lstrip().startswith("'value':"):
                    return lines[j]
    raise KeyError(name)
for obj in ("slab_operator", "mu_theta_operator"):
    text = value_line(obj)
    names = collections.Counter(re.findall(r"Symbol\('([^']+)'", text))
    print("OBJECT", obj, "VALUE_CHARS", len(text), "DISTINCT_SYMBOLS", len(names))
    by_class = collections.defaultdict(list)
    for n in sorted(names):
        by_class[value_class.get(n, "NOT_IN_LEDGER")].append(n)
    for c in sorted(by_class):
        print(" CLASS", c, len(by_class[c]))
        print("   ", " ".join(by_class[c]))
