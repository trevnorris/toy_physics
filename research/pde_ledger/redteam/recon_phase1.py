#!/usr/bin/env python3
import os, re, glob, sys
ROOT = "/var/projects/toy_physics/research/pde_ledger"
APPLY = "--apply" in sys.argv

# ---- stem -> canonical (from scripts/ filenames; stems shared across .py/.wl/.md) ----
stem2canon = {}
for p in glob.glob(os.path.join(ROOT, "scripts", "moving_throat_pde_stage*_sympy_audit.py")):
    m = re.match(r"moving_throat_pde_stage(\d{3})_(.+)_sympy_audit\.py$", os.path.basename(p))
    if m:
        stem2canon[m.group(2)] = m.group(1)
# Full-stem A citation matching. Capture stageNNN_ then the ENTIRE trailing
# identifier run, then resolve to a stem by testing the full run and the run
# minus a trailing audit-artifact suffix. This avoids the prefix false positive
# where a shorter real stem (e.g. grouped_p2 / mouth_boundary_layer) sits inside
# a longer real-but-scriptless stem (grouped_p2_status_update /
# mouth_boundary_layer_status) — those scriptless tex-only stages must NOT match.
A_re = re.compile(r"stage(\d{3})_([a-z0-9_]+)")
A_suffix_re = re.compile(r"_(sympy|mathematica)_audit(_output)?$")

def resolve_A_stem(tail):
    """Return the matched stem if `tail` (the run after stageNNN_) names a real
    stem either in full or after stripping an audit-artifact suffix; else None."""
    if tail in stem2canon:
        return tail
    m = A_suffix_re.search(tail)
    if m and tail[:m.start()] in stem2canon:
        return tail[:m.start()]
    return None

def render(canon3, stale_str):
    if stale_str.startswith('0') or len(stale_str) == 3:
        return f"{int(canon3):03d}"
    return str(int(canon3))

def own_canon(path):
    b = os.path.basename(path)
    m = re.search(r"stage(\d{3})_", b)
    if m: return m.group(1)
    m = re.search(r"stage_(\d{3})\.tex$", b)
    if m: return m.group(1)
    return None

edits = []  # (path, lineno1, cls, old, new)

def add(path, i, cls, old, new):
    if old != new:
        edits.append((path, i + 1, cls, old, new))

# ---- C.2: paper cards 091-124, line 1: Stage~MMM (both occurrences) -> own ----
for n in range(91, 125):
    path = os.path.join(ROOT, "paper", "stages", f"stage_{n:03d}.tex")
    if not os.path.exists(path): continue
    lines = open(path).read().split("\n")
    l0 = lines[0]
    new0 = re.sub(r"Stage~(\d+)", lambda m: "Stage~" + render(f"{n:03d}", m.group(1)), l0)
    add(path, 0, "C2", l0, new0)

# ---- B: notes/stages H1 self-title (first '# ' line within first 6 lines) ----
for path in glob.glob(os.path.join(ROOT, "notes", "stages", "*.md")):
    canon = own_canon(path)
    if not canon: continue
    lines = open(path).read().split("\n")
    for i, l in enumerate(lines[:6]):
        if l.startswith("# ") and re.search(r"Stage\s+\d+", l):
            new = re.sub(r"(Stage\s+)(\d+)", lambda m: m.group(1) + render(canon, m.group(2)), l, count=1)
            add(path, i, "B", l, new)
            break

# ---- C.1: .py + .wl self-banners + docstrings ----
def fix_c1(path):
    canon = own_canon(path)
    if not canon: return
    lines = open(path).read().split("\n")
    for i, l in enumerate(lines):
        nl = l
        # banner("STAGE n ... / banner('STAGE n / banner["STAGE n
        nl = re.sub(r'(banner\([\'"]STAGE )(\d+)', lambda m: m.group(1) + render(canon, m.group(2)), nl)
        nl = re.sub(r'(banner\[[\'"]STAGE )(\d+)', lambda m: m.group(1) + render(canon, m.group(2)), nl)
        # docstring: 'audit for Stage n of the moving-throat PDE'
        nl = re.sub(r'(audit for Stage )(\d+)( of the moving-throat PDE)',
                    lambda m: m.group(1) + render(canon, m.group(2)) + m.group(3), nl)
        add(path, i, "C1", l, nl)

for path in glob.glob(os.path.join(ROOT, "scripts", "*.py")):
    fix_c1(path)
for path in glob.glob(os.path.join(ROOT, "mathematica", "*.wl")):
    fix_c1(path)

# ---- A: stageMMM_<stem> citations (stem-keyed), in notes/**/*.md ----
for path in glob.glob(os.path.join(ROOT, "notes", "**", "*.md"), recursive=True):
    lines = open(path).read().split("\n")
    for i, l in enumerate(lines):
        def repl(m):
            num, tail = m.group(1), m.group(2)
            stem = resolve_A_stem(tail)
            if not stem:
                return m.group(0)
            canon = stem2canon.get(stem)
            if canon and canon != num:
                # rewrite ONLY the stageNNN_ number; preserve the full trailing run
                return "stage" + canon + "_" + tail
            return m.group(0)
        nl = A_re.sub(repl, l)
        add(path, i, "A", l, nl)

# ---- report ----
from collections import Counter
cnt = Counter(e[2] for e in edits)
print(f"PROPOSED EDITS: total={len(edits)}  by class={dict(cnt)}  (report target: A=25 B=19 C1=155 C2=34 => 233)")
print("files touched:", len(set(e[0] for e in edits)))
print()
for cls in ("C2", "B", "A", "C1"):
    rows = [e for e in edits if e[2] == cls]
    print(f"===== CLASS {cls} ({len(rows)}) =====")
    for path, ln, _, old, new in rows:
        rel = path.replace(ROOT + "/", "")
        # show only the changed token context, trimmed
        print(f"  {rel}:{ln}")
        print(f"      - {old.strip()[:140]}")
        print(f"      + {new.strip()[:140]}")
    print()

if APPLY:
    byfile = {}
    for path, ln, cls, old, new in edits:
        byfile.setdefault(path, {})[ln - 1] = new
    for path, changes in byfile.items():
        lines = open(path).read().split("\n")
        for idx, new in changes.items():
            lines[idx] = new
        open(path, "w").write("\n".join(lines))
    print(f"APPLIED {len(edits)} edits to {len(byfile)} files.")
else:
    print("DRY RUN (no files changed). Re-run with --apply to write.")
