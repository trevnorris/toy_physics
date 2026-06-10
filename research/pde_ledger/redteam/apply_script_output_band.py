#!/usr/bin/env python3
"""
Apply SCRIPT/OUTPUT-band numbering label fixes (deterministic, label-only).

Consumes the Phase-A classifier maps redteam/script_output_maps/map_*.md. Applies ONLY
rows whose CLASS column is in APPLY_CLASSES, minus an explicit EXCLUSIONS set (escalation
items the user/orchestrator decided to LEAVE). SOURCE files only (.py/.wl); .txt outputs
are refreshed separately by re-running. Every edit is a digit-only label change verified
by a strip-all-digits guard; tokens are replaced on their exact line with a negative
look-ahead so `STAGE-4` never matches inside `STAGE-44`. Idempotent. Anything ambiguous
is SKIPPED and reported -- nothing is guessed.

Width rules (empirically confirmed against already-canonical files):
  - .wl  : 3-digit zero-padded everywhere (STAGE 001, FINAL STAGE-0NN, Stage 0NN ...).
  - .py SELF prose/banner/ledger/audit (A2-SELF, non-filename): UNPADDED (STAGE 30, FINAL STAGE-21).
  - .py docstring FILENAME self-ref (lowercase `stageN`, no sep): 3-digit (stage022_...).
  - .py CROSS-ref (A3-CROSS / A4-AMBIG): 3-digit zero-padded (Stage-021, Stage 049).
"""
import os, re, glob, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
MAPS = sorted(glob.glob(os.path.join(HERE, "script_output_maps", "map_*.md")))

APPLY_CLASSES = {"A2-SELF", "A3-CROSS", "A4-AMBIG"}

# (relpath, line, token) the user/orchestrator decided to LEAVE (do NOT apply)
EXCLUSIONS = {
    ("scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py", 3, "stage4"),      # docstring old FILENAME, stem mismatch -> bogus name
    ("scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py", 121, "Stage-30"),           # M_supp defined upstream (043); non-printed comment; prior-pass ambiguous
    ("scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py", 12, "stage 65"),  # intentional "(post-renumber) former stage" history
    ("scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py", 13, "stage 69"),  # intentional history
}

DRY = "--apply" not in sys.argv

def strip_digits(s):
    return re.sub(r"\d", "", s)

def fmt_new(token, owner, engine, klass):
    m = re.match(r"(?i)(stages?)([ \-_~]?)(\d{1,3})$", token)
    if not m:
        return None
    word, sep, _ = m.groups()
    if engine == "wl":
        num = f"{owner:03d}"                              # .wl -> 3-digit everywhere
    elif token == token.lower() and sep == "":          # .py lowercase filename stem `stageN`
        num = f"{owner:03d}"
    elif klass in ("A3-CROSS", "A4-AMBIG"):             # .py cross-ref -> 3-digit
        num = f"{owner:03d}"
    else:                                                # .py A2-SELF prose/banner/ledger -> unpadded
        num = str(owner)
    return f"{word}{sep}{num}"

# parse maps -> rows (relpath, line, token, klass, owner)
rows = []
row_re = re.compile(r"^\|\s*([^|]+?)\s*\|\s*(\d+)\s*\|\s*`([^`]+)`\s*\|\s*([A-Z0-9\-]+)\s*\|\s*([0-9]+|-)\s*\|")
for mf in MAPS:
    for ln in open(mf, encoding="utf-8"):
        m = row_re.match(ln)
        if not m:
            continue
        relpath, line, token, klass, owner = m.group(1), int(m.group(2)), m.group(3), m.group(4), m.group(5)
        if klass not in APPLY_CLASSES:
            continue
        if (relpath, line, token) in EXCLUSIONS:
            continue
        if owner == "-":
            print(f"  ?? no owner for apply-row {relpath}:{line} `{token}` ({klass}) — SKIP")
            continue
        owner = int(owner)
        rows.append((relpath, line, token, klass, owner))

# group by file, then line
byfile = defaultdict(lambda: defaultdict(list))
for relpath, line, token, klass, owner in rows:
    byfile[relpath][line].append((token, klass, owner))

applied = skipped = files_changed = 0
skips = []
for relpath in sorted(byfile):
    path = os.path.join(ROOT, relpath)
    engine = "wl" if relpath.endswith(".wl") else "py"
    if not os.path.exists(path):
        skips.append(f"{relpath}: FILE NOT FOUND"); continue
    lines = open(path, encoding="utf-8").read().split("\n")
    changed = False
    for lineno, edits in sorted(byfile[relpath].items()):
        if lineno < 1 or lineno > len(lines):
            for t, k, o in edits:
                skips.append(f"{relpath}:{lineno} `{t}` line out of range")
            continue
        orig = lines[lineno-1]
        cur = orig
        for token, klass, owner in edits:
            new = fmt_new(token, owner, engine, klass)
            if new is None:
                skips.append(f"{relpath}:{lineno} `{token}` unparsable token"); skipped += 1; continue
            if strip_digits(token) != strip_digits(new):
                skips.append(f"{relpath}:{lineno} NON-LABEL-ONLY `{token}`->`{new}`"); skipped += 1; continue
            if new == token:
                skipped += 1; continue
            pat = re.escape(token) + r"(?!\d)"
            cur2, n = re.subn(pat, new, cur, count=1)
            if n == 0:
                # idempotent? already applied
                if token not in cur and new in cur:
                    pass
                else:
                    skips.append(f"{relpath}:{lineno} `{token}` NOT FOUND on line (got: {cur.strip()[:70]!r})"); skipped += 1
                continue
            cur = cur2
            applied += 1
        if cur != orig:
            lines[lineno-1] = cur
            changed = True
    if changed:
        files_changed += 1
        if not DRY:
            open(path, "w", encoding="utf-8").write("\n".join(lines))

print(f"\n{'DRY-RUN (no writes; pass --apply to write)' if DRY else 'APPLIED'}")
print(f"  edits applied: {applied}   skipped: {skipped}   files changed: {files_changed}")
print(f"  total apply-rows parsed: {len(rows)}")
if skips:
    print("  SKIPS / WARNINGS:")
    for s in skips:
        print("    -", s)
