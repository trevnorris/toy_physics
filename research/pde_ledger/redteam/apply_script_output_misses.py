#!/usr/bin/env python3
"""
Phase-D miss completion: stale cross-ref LABELS the per-line classifier flagged on
SOME lines but missed on others (same token, additional occurrence), or where a
capital `Stage18` LABEL coexists with the lowercase/camelCase identifier
`F_stage18`/`fStage18` (VAR-LEAVE) on one line and only the identifier was logged.

Per the user rule: LEAVE code identifiers, FIX adjacent labels. Each edit is digit-only
(strip-digits guard) and matched with a negative LOOK-BEHIND so `Stage18` never matches
inside `fStage18`/`gStage19`, and a trailing non-digit guard.

Owners (content-verified): Stage18 F = canon 035 (F dimensionless shape), Stage19 G =
canon 036 (G support-feasibility), Stage 22 split-U continuum = canon 039,
Stage-27 branch equation = canon 044. SOURCE only; outputs refresh on re-run.
"""
import os, re, sys
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DRY = "--apply" not in sys.argv

EDITS = [  # (relpath, line, old, new)
    ("scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py", 101, "Stage 22", "Stage 039"),
    ("scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py", 121, "Stage18", "Stage035"),
    ("scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py", 122, "Stage19", "Stage036"),
    ("mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl", 128, "Stage18", "Stage035"),
    ("mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl", 129, "Stage19", "Stage036"),
    ("scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py", 128, "Stage-27", "Stage-044"),
]
sd = lambda s: re.sub(r"\d", "", s)
applied, skips = 0, []
from collections import defaultdict
byfile = defaultdict(list)
for rp, ln, o, n in EDITS:
    byfile[rp].append((ln, o, n))
for rp, edits in sorted(byfile.items()):
    path = os.path.join(ROOT, rp)
    lines = open(path, encoding="utf-8").read().split("\n")
    changed = False
    for ln, old, new in edits:
        if sd(old) != sd(new):
            skips.append(f"{rp}:{ln} NON-LABEL-ONLY {old!r}->{new!r}"); continue
        cur = lines[ln-1]
        pat = r"(?<![A-Za-z0-9_])" + re.escape(old) + r"(?!\d)"
        cur2, k = re.subn(pat, new, cur, count=1)
        if k != 1:
            if old not in cur and new in cur: continue
            skips.append(f"{rp}:{ln} `{old}` matched {k}x: {cur.strip()[:70]!r}"); continue
        lines[ln-1] = cur2; applied += 1; changed = True
    if changed and not DRY:
        open(path, "w", encoding="utf-8").write("\n".join(lines))
print(f"{'DRY-RUN' if DRY else 'APPLIED'}: {applied}/{len(EDITS)} miss fixes")
for s in skips: print("  SKIP", s)
