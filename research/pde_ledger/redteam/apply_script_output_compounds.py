#!/usr/bin/env python3
"""
Phase-D compound-ref completion for the SCRIPT/OUTPUT band.

The Phase-A scan regex captured only the FIRST number of a compound stage ref
(`Stage-X/Y`, `Stages X-Y`), so the per-token applier left the trailing number(s)
stale (e.g. `Stage-3/4` -> half-fixed `Stage-003/4`). This script completes those
compounds with EXPLICIT, content-verified old->new strings. Each trailing number maps
to its canonical owner exactly as its already-fixed partner did (old K -> 021 if K=4,
else K+17), confirmed against the cited deliverable. Digit-only (strip-digits guard);
each old substring is matched once on its exact line with a trailing non-digit guard.

Trailing-owner content notes:
  /4->021 (old-4 maxwell/mixed), /5->022 (grouped P2), /9->026 (concrete-axial branch),
  /19->036 (G support-feasibility, partner of F@035; stage037 note says "034/036"),
  /25->042 (rank2 normalization law), /46->063 (parent thresholds, adj. to 062),
  /63->080 (zeta thresholds), /64->081 (Pi thresholds / support ceiling),
  stage036 /019->036 (F@035 / G@036 share xi,delta,M_mix — old-19=self 036).
"""
import os, re, sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DRY = "--apply" not in sys.argv

# (relpath, line, old, new)
EDITS = [
    ("mathematica/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.wl", 109, "STAGE-025/9", "STAGE-025/026"),
    ("mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl", 191, "Stage-078/63/64", "Stage-078/080/081"),
    ("scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py", 12, "Stage-003/4", "Stage-003/021"),
    ("scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py", 338, "Stage-003/4", "Stage-003/021"),
    ("scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py", 345, "Stage-021/5", "Stage-021/022"),
    ("scripts/moving_throat_pde_stage027_nonconstant_axial_family_sympy_audit.py", 18, "Stage-025/9", "Stage-025/026"),
    ("scripts/moving_throat_pde_stage027_nonconstant_axial_family_sympy_audit.py", 179, "Stage-025/9", "Stage-025/026"),
    ("scripts/moving_throat_pde_stage027_nonconstant_axial_family_sympy_audit.py", 183, "STAGE-025/9", "STAGE-025/026"),
    ("scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py", 12, "Stage 035/019", "Stage 035/036"),
    ("scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py", 6, "Stage-034/19", "Stage-034/036"),
    ("scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py", 224, "Stage-034/19", "Stage-034/036"),
    ("scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py", 14, "Stage-035/19", "Stage-035/036"),
    ("scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py", 160, "Stage-035/19", "Stage-035/036"),
    ("scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py", 121, "Stage 041/25", "Stage 041/042"),
    ("scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py", 25, "Stage-062/46", "Stage-062/063"),
    ("scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py", 180, "Stage-062/46", "Stage-062/063"),
    ("scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py", 171, "Stage-078/63/64", "Stage-078/080/081"),
    ("scripts/moving_throat_pde_stage086_family1_loading_ratio_window_sympy_audit.py", 37, "Stages 080-64", "Stages 080-081"),
    ("scripts/moving_throat_pde_stage090_updated_reduced_status_sympy_audit.py", 18, "Stage 080/64", "Stage 080/081"),
    ("scripts/moving_throat_pde_stage090_updated_reduced_status_sympy_audit.py", 69, "Stages 080/64", "Stages 080/081"),
]

def strip_digits(s):
    return re.sub(r"\d", "", s)

applied = 0
skips = []
from collections import defaultdict
byfile = defaultdict(list)
for rp, line, old, new in EDITS:
    byfile[rp].append((line, old, new))

for rp, edits in sorted(byfile.items()):
    path = os.path.join(ROOT, rp)
    lines = open(path, encoding="utf-8").read().split("\n")
    changed = False
    for line, old, new in edits:
        if strip_digits(old) != strip_digits(new):
            skips.append(f"{rp}:{line} NON-LABEL-ONLY {old!r}->{new!r}"); continue
        cur = lines[line-1]
        cur2, n = re.subn(re.escape(old) + r"(?!\d)", new, cur, count=1)
        if n != 1:
            if old not in cur and new in cur:
                continue  # idempotent
            skips.append(f"{rp}:{line} `{old}` matched {n}x (need 1): {cur.strip()[:70]!r}"); continue
        lines[line-1] = cur2
        applied += 1
        changed = True
    if changed and not DRY:
        open(path, "w", encoding="utf-8").write("\n".join(lines))

print(f"{'DRY-RUN' if DRY else 'APPLIED'}: {applied}/{len(EDITS)} compound trailing-number fixes")
for s in skips:
    print("  SKIP", s)
