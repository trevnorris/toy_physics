#!/usr/bin/env python3
"""
Apply SECONDARY-pass Cluster-2 label-only numbering maps (band 001-090 + 164-178).

Consumes redteam/secondary_maps/c2_*.md (each block: `### <note>` then
`- OLD: \`...\`` / `  NEW: \`...\``) plus an inline EXTRA list of verifier-found
coverage-gap fixes. For every (note, OLD, NEW):
  - confirm OLD occurs EXACTLY once in the note's CURRENT content (deterministic),
  - confirm the edit is label-only (OLD and NEW identical after stripping digits),
  - replace that single occurrence.
Edits are grouped per note and written once. Any non-unique / not-found / non-label-only
edit is SKIPPED and reported; nothing is guessed.
"""
import os, re, glob
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
MAPS = os.path.join(HERE, "secondary_maps")
NOTES = os.path.join(HERE, "..", "notes", "stages")

# Verifier-found coverage-gap fixes (stage_num, OLD, NEW) — 3-digit padded, label-only.
EXTRA = [
    ("024", "the Stage-6 radial/axial overlap objects", "the Stage-023 radial/axial overlap objects"),
    ("024", "the full Stage-6 coefficients", "the full Stage-023 coefficients"),
    ("024", "just the Stage-6 combined gauge/mixed coupling strength", "just the Stage-023 combined gauge/mixed coupling strength"),
    ("026", "at the end of Stage 8.", "at the end of Stage 025."),
    ("028", "into the Stage-10 normalization equation", "into the Stage-027 normalization equation"),
    ("036", "Stages 18 and 19 together now give", "Stages 035 and 036 together now give"),
    ("064", "## 6. What Stage 47 changes", "## 6. What Stage 064 changes"),
    ("070", "## 5. What Stage 53 changes", "## 5. What Stage 070 changes"),
    ("071", "## 5. What Stage 54 changes", "## 5. What Stage 071 changes"),
    ("072", "## 4. What Stage 55 changes", "## 4. What Stage 072 changes"),
    ("073", "## 5. What Stage 56 changes", "## 5. What Stage 073 changes"),
    ("074", "## 6. What Stage 57 changes", "## 6. What Stage 074 changes"),
    ("075", "## 5. What Stage 58 changes", "## 5. What Stage 075 changes"),
    ("076", "## 5. What Stage 59 changes", "## 5. What Stage 076 changes"),
    ("077", "## 6. What Stage 60 changes", "## 6. What Stage 077 changes"),
    ("078", "## 4. What remains open after Stage 61", "## 4. What remains open after Stage 078"),
    ("079", "## 5. What Stage 62 changes", "## 5. What Stage 079 changes"),
]

def strip_digits(s):
    return re.sub(r"\d", "", s)

def resolve(note):
    """note is either a full basename or a 3-digit stage number."""
    if note.endswith(".md"):
        p = os.path.join(NOTES, note)
        return p if os.path.exists(p) else None
    hits = glob.glob(os.path.join(NOTES, f"moving_throat_pde_stage{note}_*.md"))
    return hits[0] if len(hits) == 1 else None

# (path, OLD, NEW) edits, grouped per path
edits = defaultdict(list)
parse_errors = []

for mf in sorted(glob.glob(os.path.join(MAPS, "c2_*.md"))):
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

for num, old, new in EXTRA:
    p = resolve(num)
    if not p:
        parse_errors.append(f"EXTRA: cannot resolve stage {num}")
    else:
        edits[p].append((old, new))

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
                continue  # already applied (idempotent re-run): old gone, new present
            skipped.append(f"{os.path.basename(p)}: count={c} (need 1) for {old!r}")
            continue
        content = content.replace(old, new, 1)
        applied += 1
        changed = True
    if changed:
        open(p, "w", encoding="utf-8").write(content)

print(f"edits parsed: {total} (+{len(EXTRA)} extra incl.)  notes touched: {sum(1 for p in edits if True)}")
print(f"APPLIED: {applied}")
print(f"SKIPPED: {len(skipped)}")
for s in skipped:
    print("  SKIP", s)
for e in parse_errors:
    print("  PARSE", e)
