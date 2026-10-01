#!/usr/bin/env python3
"""Check every repo line-citation in A and B against the grounding files and sources."""
import re
from pathlib import Path

A = Path("/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_clean_condition.md").read_text()
B = Path("/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_zinvariant_operator_blocks_directive.md").read_text()
GA = Path("/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition.md").read_text()
GB = Path("/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11c_d_zinvariant_operator_blocks_directive.md").read_text()

# Citations of the form path:start-end or `:start-end` or spec `:n-m`
cite_re = re.compile(
    r"((?:[\w./-]+\.(?:md|py|wl|txt))?)\s*`:(\d+)(?:[–-](\d+))?`"
)

def cites(text):
    out = []
    for m in cite_re.finditer(text):
        path, a, b = m.group(1), int(m.group(2)), m.group(3)
        b = int(b) if b else a
        out.append((path, a, b, m.group(0)))
    return out

print("A_CITE_COUNT", len(cites(A)))
print("B_CITE_COUNT", len(cites(B)))

# Grounding files should contain each numeric range as "sed -n start,endp" or "sed -n startp"
def grounded(gtext, start, end):
    if start == end:
        pats = [f"sed -n {start}p", f"sed -n '{start}p'", f":{start}", f"{start}p"]
    else:
        pats = [
            f"sed -n {start},{end}p",
            f"sed -n '{start},{end}p'",
            f"{start},{end}p",
            f":{start}–{end}",
            f":{start}-{end}",
            f"`:{start}–{end}`",
            f"{start}–{end}",
        ]
    return any(p in gtext for p in pats)

print("--- A citations vs A grounding ---")
for path, a, b, raw in cites(A):
    ok = grounded(GA, a, b)
    # also accept if the range is a subset of a dumped range in the grounding
    print(("OK" if ok else "MISSING_IN_GROUNDING"), raw, "path=", path or "(bare)")

print("--- B citations vs B grounding ---")
for path, a, b, raw in cites(B):
    ok = grounded(GB, a, b)
    print(("OK" if ok else "MISSING_IN_GROUNDING"), raw, "path=", path or "(bare)")

# Specific known-risk checks against the files
p2b = Path("/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_p2b_gamma_bridge_directive.md")
print("P2B_IN_B_GROUNDING", "S11c_b_p2b_gamma_bridge_directive.md" in GB)
print("P2B_EXISTS", p2b.exists())
if p2b.exists():
    lines = p2b.read_text().splitlines()
    print("P2B_LINE_10", lines[9][:120] if len(lines) >= 10 else "SHORT")
    print("P2B_LINE_13", lines[12][:160] if len(lines) >= 13 else "SHORT")
