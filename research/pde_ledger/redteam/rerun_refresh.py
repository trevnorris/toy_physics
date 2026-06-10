#!/usr/bin/env python3
"""
Re-run scripts and refresh their committed .txt outputs (Phase B/C output refresh).
Usage: rerun_refresh.py <py|wl> stageNNN [stageNNN ...]
Runs each stage's source for the given engine (timeout 600), and on exit 0 overwrites
the committed output .txt with fresh stdout. Reports exit codes; does NOT write on failure.
"""
import sys, os, glob, subprocess

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
engine = sys.argv[1]
stages = sys.argv[2:]
SRCDIR = "scripts" if engine == "py" else "mathematica"
EXT = ".py" if engine == "py" else ".wl"
OUTDIR = os.path.join(SRCDIR, "output")

fail = []
for s in stages:
    src = glob.glob(os.path.join(ROOT, SRCDIR, f"moving_throat_pde_stage{s}_*{EXT}"))
    if len(src) != 1:
        print(f"  stage{s} {engine}: SRC AMBIGUOUS/MISSING ({len(src)})"); fail.append(s); continue
    src = src[0]
    base = os.path.basename(src)[:-len(EXT)]
    out = os.path.join(ROOT, OUTDIR, base + ".txt")
    cmd = ["python3", src] if engine == "py" else ["math", "-script", src]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd=os.path.join(ROOT, SRCDIR))
    except subprocess.TimeoutExpired:
        print(f"  stage{s} {engine}: TIMEOUT (>600s)"); fail.append(s); continue
    if r.returncode != 0:
        print(f"  stage{s} {engine}: EXIT {r.returncode}  stderr: {r.stderr.strip()[:160]}"); fail.append(s); continue
    if not os.path.exists(out):
        print(f"  stage{s} {engine}: OUTPUT PATH MISSING {os.path.relpath(out, ROOT)}"); fail.append(s); continue
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(r.stdout)
    print(f"  stage{s} {engine}: OK -> {os.path.relpath(out, ROOT)}")

print(f"\n{engine}: {len(stages)-len(fail)}/{len(stages)} refreshed; failures: {fail}")
sys.exit(1 if fail else 0)
