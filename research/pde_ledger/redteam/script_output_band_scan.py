#!/usr/bin/env python3
"""
SCRIPT / OUTPUT-band numbering-reconciliation INVENTORY generator (read-only).

Scope (this band ONLY): scripts/*.py, mathematica/*.wl, scripts/output/*.txt,
mathematica/output/*.txt. The notes/.tex band is DONE (NUMBERING_BROAD_SWEEP_PLAN.md)
and is NOT touched here.

Ground truth = the canonical stage number in each file's FILENAME (stage(\\d{3})_).
Enumerates every in-CONTENT stage-number token whose number != that file's canonical
number -- the candidate pool -- and buckets each by a PRIORITY heuristic to drive the
Phase-A classifier agents. NEVER decides stale-vs-canonical mechanically; content work
is the agents' job (a cited number != canonical may be a legitimate cross-ref). The
buckets are triage, not verdicts. NEVER offset-sweep.

Token forms captured (broad, per the plan): Stage/STAGE/stage/Stages + any single
separator (space, -, _, ~, none) + number + optional sub-stage `.k`; including embedded
identifiers (chiFromStage180, J1_stage48, F_stage18) and old-filename stems (stage8_).

Priority buckets (each candidate assigned to exactly one, top-down):
  VARIABLE   : ident-form (token embedded in a code identifier / variable name), Δ≠0
               -> CODE identifiers; renaming brushes the IRON RULE "zero variable bytes"
               -> USER DECISION (rename-consistently vs leave). Includes chiFromStage180.
  SELF       : line matches a self-label pattern (banner/pass-line/ledger/AUDIT/docstring), Δ≠0
               -> A1 (output stale) or A2 (source self-label). The dominant actionable class.
  FORWARD    : cited > canonical (a script citing a LATER stage), not self/ident
               -> A3 high-signal (stale forward-label) OR genuine forward pointer. Content.
  OLD-BACK   : cited < canonical, cited token is 1-2 digit (unpadded), not self/ident
               -> A3/A4 old-epoch candidate (e.g. `Stage-4`=old-4). Content + padding signal.
  CANON-BACK : cited < canonical, cited token is 3-digit (zero-padded canonical form)
               -> mostly GENUINE upstream citations; default LEAVE, agents spot-check.
"""
import os, re, glob
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(HERE, "NUMBERING_SCRIPT_OUTPUT_BAND_WORKLIST.md")
AGENTDIR = os.path.join(HERE, "script_output_maps")

DIRS = [
    ("py",     os.path.join(ROOT, "scripts"),                "*.py"),
    ("wl",     os.path.join(ROOT, "mathematica"),            "*.wl"),
    ("py-out", os.path.join(ROOT, "scripts", "output"),      "*.txt"),
    ("wl-out", os.path.join(ROOT, "mathematica", "output"),  "*.txt"),
]

re_canon = re.compile(r"stage(\d{3})_")
re_tok = re.compile(r"(?i)stages?[ \-_~]?(\d{1,3})(?:\.(\d{1,3}))?")
re_self_line = re.compile(
    r"(?i)(AUDIT COMPLETE|LEDGER|checks? passed|all .*pass|^\s*\(?\*?\s*STAGE[ \-]?\d"
    r"|_audit\.(py|wl|txt)|^\s*STATUS:|banner\()")

def band(n):
    for lo, hi in [(1, 90), (91, 180), (181, 253)]:
        if lo <= n <= hi:
            return f"{lo:03d}-{hi:03d}"
    return "other"

def bucket(r):
    if r["ident"] and r["delta"] != 0:
        return "VARIABLE"
    if r["looks_self"] and r["delta"] != 0:
        return "SELF"
    if r["delta"] > 0:
        return "FORWARD"
    if r["delta"] < 0 and r["digw"] <= 2:
        return "OLD-BACK"
    if r["delta"] < 0 and r["digw"] == 3:
        return "CANON-BACK"
    return "OTHER"

rows = []
masters = []
selfmatch = 0
for kind, d, pat in DIRS:
    for path in sorted(glob.glob(os.path.join(d, pat))):
        fn = os.path.basename(path)
        m = re_canon.search(fn)
        canon = int(m.group(1)) if m else None
        for i, line in enumerate(open(path, encoding="utf-8", errors="replace").readlines(), 1):
            for mm in re_tok.finditer(line):
                cited = int(mm.group(1))
                if not (1 <= cited <= 253):
                    continue
                start = mm.start()
                prev = line[start-1] if start else ""
                ident = bool(re.match(r"[A-Za-z0-9_]", prev))
                if canon is not None and cited == canon and not ident:
                    selfmatch += 1
                    continue
                tok = line[start:mm.end()].strip()
                after = line[mm.end():mm.end()+2]
                rec = dict(kind=kind, fn=fn, canon=canon, line=i, cited=cited, tok=tok,
                           digw=len(mm.group(1)), ident=ident,
                           fileform=bool(re.match(r"_[a-z]", after)),
                           is_out=kind.endswith("out"),
                           looks_self=bool(re_self_line.search(line)),
                           delta=(cited - canon) if canon is not None else 0,
                           ctx=line.rstrip("\n"))
                if canon is None:
                    masters.append(rec)
                else:
                    rows.append(rec)

# de-dup
seen = set(); uniq = []
for r in rows:
    k = (r["fn"], r["line"], r["cited"], r["tok"])
    if k in seen:
        continue
    seen.add(k); uniq.append(r)
rows = uniq
for r in rows:
    r["bucket"] = bucket(r)

def stem(fn):
    return (fn.replace("moving_throat_pde_", "")
              .replace("_sympy_audit", "").replace("_mathematica_audit", "")
              .replace(".py", "").replace(".wl", "").replace(".txt", ""))

BUCKETS = ["SELF", "FORWARD", "OLD-BACK", "VARIABLE", "CANON-BACK", "OTHER"]
byband = defaultdict(list)
for r in rows:
    byband[band(r["canon"])].append(r)

def fmt(r):
    snip = r["ctx"].strip()
    snip = (snip[:120] + "…") if len(snip) > 120 else snip
    snip = snip.replace("|", "\\|")
    fl = []
    if r["is_out"]: fl.append("OUT")
    if r["fileform"]: fl.append("stem")
    return (f"| {r['canon']:03d} | {r['kind']} | {r['line']} | `{r['tok']}` | "
            f"{r['cited']} | {r['delta']:+d} | {','.join(fl) or '-'} | {snip} |")

# ---- main worklist ----
with open(OUT, "w", encoding="utf-8") as fh:
    fh.write("# Numbering — SCRIPT/OUTPUT band candidate INVENTORY (read-only, machine-generated)\n\n")
    fh.write("Generated by `redteam/script_output_band_scan.py`. Ground truth = filename "
             "`stageNNN_`; NEVER offset-sweep. `Δ` = cited−canonical DISTANCE (NOT a staleness "
             "offset). Buckets are TRIAGE for the Phase-A classifier agents, not verdicts.\n\n")
    fh.write(f"- total candidate tokens (number ≠ canonical): **{len(rows)}**\n")
    fh.write(f"- distinct files with ≥1 candidate: **{len(set(r['fn'] for r in rows))}**\n")
    fh.write(f"- canonical self-matches skipped: {selfmatch}; master-file tokens: {len(masters)}\n\n")
    fh.write("## counts by bucket × band\n\n")
    fh.write("| bucket | 001-090 | 091-180 | 181-253 | total |\n| --- | ---: | ---: | ---: | ---: |\n")
    for b in BUCKETS:
        cells = [sum(1 for r in byband[bd] if r["bucket"] == b) for bd in ("001-090","091-180","181-253")]
        tot = sum(cells)
        if tot:
            fh.write(f"| {b} | {cells[0]} | {cells[1]} | {cells[2]} | {tot} |\n")
    fh.write(f"| **TOTAL** | {len(byband['001-090'])} | {len(byband['091-180'])} | {len(byband['181-253'])} | {len(rows)} |\n\n")
    fh.write("## bucket meaning\n\n"
             "- **SELF** = self-label (banner/pass-line/ledger/AUDIT/docstring), Δ≠0 → A1 (output stale, source canonical → refresh) or A2 (source self-label stale → fix+rerun).\n"
             "- **FORWARD** = cites a LATER stage → A3 (stale forward-label) or genuine forward pointer (content).\n"
             "- **OLD-BACK** = cites earlier stage by a 1-2-digit (unpadded) number → A3/A4 old-epoch candidate (content; padding is a signal, not a verdict).\n"
             "- **CANON-BACK** = cites earlier stage by a 3-digit padded number → mostly GENUINE upstream citation; default LEAVE, spot-check.\n"
             "- **VARIABLE** = stage number embedded in a CODE identifier (`A_stage4`, `chiFromStage180`) → IRON-RULE-sensitive (zero variable bytes) → USER DECISION.\n\n")
    fh.write("---\n\n")
    for bd in ("001-090", "091-180", "181-253"):
        fh.write(f"# BAND {bd}\n\n")
        bybuck = defaultdict(list)
        for r in byband[bd]:
            bybuck[r["bucket"]].append(r)
        for b in BUCKETS:
            if not bybuck[b]:
                continue
            fh.write(f"## {bd} · {b}  ({len(bybuck[b])})\n\n")
            fh.write("| canon | kind | line | token | cited | Δ | flags | context |\n")
            fh.write("| ---: | --- | ---: | --- | ---: | ---: | --- | --- |\n")
            for r in sorted(bybuck[b], key=lambda r: (r["canon"], r["kind"], r["line"])):
                fh.write(fmt(r) + "\n")
            fh.write("\n")

# ---- per-band agent input files (actionable buckets only; CANON-BACK collapsed) ----
os.makedirs(AGENTDIR, exist_ok=True)
for bd in ("001-090", "091-180", "181-253"):
    p = os.path.join(AGENTDIR, f"candidates_{bd}.md")
    bybuck = defaultdict(list)
    for r in byband[bd]:
        bybuck[r["bucket"]].append(r)
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(f"# Phase-A classifier input — BAND {bd} (read-only)\n\n")
        fh.write("Adjudicate each ACTIONABLE row by CONTENT (which stage owns the cited deliverable?); "
                 "NEVER offset-sweep. Ground truth = filename `stageNNN_`. Report A1/A2/A3/A4 + the "
                 "per-ref owning-stage map (or LEAVE/cannot-confirm).\n\n")
        for b in ["SELF", "FORWARD", "OLD-BACK", "VARIABLE"]:
            fh.write(f"## {b}  ({len(bybuck[b])})\n\n")
            fh.write("| canon | kind | line | token | cited | Δ | flags | context |\n")
            fh.write("| ---: | --- | ---: | --- | ---: | ---: | --- | --- |\n")
            for r in sorted(bybuck[b], key=lambda r: (r["canon"], r["kind"], r["line"])):
                fh.write(fmt(r) + "\n")
            fh.write("\n")
        fh.write(f"## CANON-BACK ({len(bybuck['CANON-BACK'])}) — 3-digit padded backward refs; default LEAVE, spot-check only\n\n")
        cb = defaultdict(int)
        for r in bybuck["CANON-BACK"]:
            cb[(r["canon"], stem(r["fn"]))] += 1
        fh.write("| canon | stem | count |\n| ---: | --- | ---: |\n")
        for (c, st), n in sorted(cb.items()):
            fh.write(f"| {c:03d} | {st} | {n} |\n")

print(f"candidates={len(rows)}  files={len(set(r['fn'] for r in rows))}  masters={len(masters)}")
gb = defaultdict(lambda: defaultdict(int))
for r in rows:
    gb[band(r["canon"])][r["bucket"]] += 1
for bd in ("001-090", "091-180", "181-253"):
    parts = "  ".join(f"{b}={gb[bd][b]}" for b in BUCKETS if gb[bd][b])
    print(f"  {bd}: total={len(byband[bd])}  {parts}")
print("ACTIONABLE (SELF+FORWARD+OLD-BACK+VARIABLE) totals:")
for b in ["SELF","FORWARD","OLD-BACK","VARIABLE","CANON-BACK"]:
    print(f"   {b}: {sum(1 for r in rows if r['bucket']==b)}")
print(f"written -> {os.path.relpath(OUT, ROOT)} + per-band files in {os.path.relpath(AGENTDIR, ROOT)}/")
