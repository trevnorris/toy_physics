#!/usr/bin/env python3
"""
Apply SECONDARY-pass Cluster-3 label-only numbering maps (compound/forward tail).

Sibling of apply_secondary_maps.py with the SAME deterministic guards but globbing
redteam/secondary_maps/c3_*.md and NO inline EXTRA list. For every (note, OLD, NEW):
  - confirm OLD occurs EXACTLY once in the note's CURRENT content,
  - confirm the edit is label-only (OLD == NEW after stripping digits),
  - replace that single occurrence.
Idempotent: a re-run where OLD is gone and NEW present is skipped silently.
Any non-unique / not-found / non-label-only edit is SKIPPED and reported; nothing guessed.
"""
import os, re, glob
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
MAPS = os.path.join(HERE, "secondary_maps")
NOTES = os.path.join(HERE, "..", "notes", "stages")

def strip_digits(s):
    return re.sub(r"\d", "", s)

def resolve(note):
    if note.endswith(".md"):
        p = os.path.join(NOTES, note)
        return p if os.path.exists(p) else None
    hits = glob.glob(os.path.join(NOTES, f"moving_throat_pde_stage{note}_*.md"))
    return hits[0] if len(hits) == 1 else None

edits = defaultdict(list)
parse_errors = []

for mf in sorted(glob.glob(os.path.join(MAPS, "c3_*.md"))):
    cur = None
    pending_old = None
    for raw in open(mf, encoding="utf-8"):
        line = raw.rstrip("\n")
        if line.startswith("### "):
            cur = line[4:].strip()
            pending_old = None
        elif line.lstrip().startswith("- OLD:"):
            v = line.split(":", 1)[1].strip()
            pending_old = v[1:-1] if v.startswith("`") and v.endswith("`") else v
        elif line.lstrip().startswith("NEW:"):
            v = line.split(":", 1)[1].strip()
            new = v[1:-1] if v.startswith("`") and v.endswith("`") else v
            if cur is None or pending_old is None:
                parse_errors.append(f"{os.path.basename(mf)}: NEW without OLD/### ({new!r})")
            else:
                p = resolve(cur)
                if not p:
                    parse_errors.append(f"{os.path.basename(mf)}: cannot resolve note {cur!r}")
                else:
                    edits[p].append((pending_old, new))
            pending_old = None

applied = 0
skipped = []
total = sum(len(v) for v in edits.values())
for p, lst in sorted(edits.items()):
    content = open(p, encoding="utf-8").read()
    changed = False
    for old, new in lst:
        if strip_digits(old) != strip_digits(new):
            skipped.append(f"{os.path.basename(p)}: NON-LABEL-ONLY {old!r} -> {new!r}")
            continue
        c = content.count(old)
        if c != 1:
            if c == 0 and content.count(new) >= 1:
                continue  # already applied (idempotent re-run)
            skipped.append(f"{os.path.basename(p)}: count={c} (need 1) for {old!r}")
            continue
        content = content.replace(old, new, 1)
        applied += 1
        changed = True
    if changed:
        open(p, "w", encoding="utf-8").write(content)

print(f"edits parsed: {total}  notes touched: {len(edits)}")
print(f"APPLIED: {applied}")
print(f"SKIPPED: {len(skipped)}")
for s in skipped:
    print("  SKIP", s)
for e in parse_errors:
    print("  PARSE", e)
