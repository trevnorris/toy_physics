#!/usr/bin/env python3
"""
Band 091-180 SCRIPT/OUTPUT numbering fixes — explicit, orchestrator-vetted edit list.

Unlike the band-1 map-driven applier, band 091-180 has only a handful of content-stale
cross-refs, each individually verified by the orchestrator against the cited stage's
filename-stem + paper card + corroborating siblings. Source-only (.py/.wl); .txt outputs
are refreshed by re-running afterward. Every edit is digit-only (strip-all-digits guard),
applied on its exact line with a negative look-ahead so `Stage 98` never matches inside
`Stage 980`. Idempotent. Run with --apply to write.

Vetted edits (see redteam/script_output_maps/BAND_091-180_RESULT.md for the evidence):
  #2  stage112.wl:55   Stage-92  -> Stage-109   (old-epoch +17; stage109=linearized_branch_selection
                                                 owns the (b,a0,a5) cross-check w/ a0/3+9a5; stage092=
                                                 dynamic_geometry_obstruction does not)
  #3  stage116 x6      Stage 98  -> Stage 115    (old-epoch +17; stage115=core_balance_compensation,
                                                 line 47 gamma0_can=(1+r_c)/9, owns the "compensation
                                                 requirement"; stage098=family1_support does not)
  #4  stage134 x2      Stage 137 -> Stage 140    (susceptibility closure owned by stage140=
                                                 selfmatched_mouth_susceptibility, corroborated by
                                                 stage139:10 "established at Stage 140, not here" +
                                                 stage142:37; stage137=core_to_mouth_gain_map does not)
"""
import os, re, sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DRY = "--apply" not in sys.argv

# (relpath, line, old_token, new_token)
EDITS = [
    # #2
    ("mathematica/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl", 55, "Stage-92", "Stage-109"),
    # #3  (.py)
    ("scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py", 68, "Stage 98", "Stage 115"),
    ("scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py", 73, "Stage 98", "Stage 115"),
    ("scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py", 79, "Stage 98", "Stage 115"),
    # #3  (.wl)
    ("mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl", 67, "Stage 98", "Stage 115"),
    ("mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl", 70, "Stage 98", "Stage 115"),
    ("mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl", 76, "Stage 98", "Stage 115"),
    # #4
    ("scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py", 90, "Stage 137", "Stage 140"),
    ("mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl", 101, "Stage 137", "Stage 140"),
]

def strip_digits(s):
    return re.sub(r"\d", "", s)

# group by file
byfile = {}
for relpath, line, old, new in EDITS:
    byfile.setdefault(relpath, []).append((line, old, new))

applied = skipped = files_changed = 0
skips = []
for relpath in sorted(byfile):
    path = os.path.join(ROOT, relpath)
    if not os.path.exists(path):
        skips.append(f"{relpath}: FILE NOT FOUND"); continue
    lines = open(path, encoding="utf-8").read().split("\n")
    changed = False
    for line, old, new in sorted(byfile[relpath]):
        if line < 1 or line > len(lines):
            skips.append(f"{relpath}:{line} `{old}` out of range"); skipped += 1; continue
        if strip_digits(old) != strip_digits(new):
            skips.append(f"{relpath}:{line} NON-DIGIT-ONLY `{old}`->`{new}`"); skipped += 1; continue
        cur = lines[line-1]
        pat = re.escape(old) + r"(?!\d)"
        cur2, n = re.subn(pat, new, cur, count=1)
        if n == 0:
            if old not in cur and new in cur:
                skipped += 1  # idempotent (already applied)
            else:
                skips.append(f"{relpath}:{line} `{old}` NOT FOUND (got: {cur.strip()[:70]!r})"); skipped += 1
            continue
        lines[line-1] = cur2
        applied += 1
        changed = True
    if changed:
        files_changed += 1
        if not DRY:
            open(path, "w", encoding="utf-8").write("\n".join(lines))

print("DRY-RUN (pass --apply to write)" if DRY else "APPLIED")
print(f"  edits applied: {applied}   skipped: {skipped}   files changed: {files_changed}   total: {len(EDITS)}")
for s in skips:
    print("    -", s)
