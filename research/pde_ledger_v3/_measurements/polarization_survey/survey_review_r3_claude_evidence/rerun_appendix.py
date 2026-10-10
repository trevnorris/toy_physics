#!/usr/bin/env python3
# Re-run every Appendix A command in the survey and diff against its recorded literal output.
import re, subprocess, sys
path = "/var/projects/toy_physics/_scratch/polarization/POLARIZATION_SURVEY.md"
lines = open(path, encoding="utf-8").read().split("\n")
start = next(i for i,l in enumerate(lines) if l.startswith("## Appendix A"))
blocks = []
i = start
cur = None
while i < len(lines):
    l = lines[i]
    m = re.match(r'<a id="([a-z0-9]+)"></a>', l)
    if m:
        cur = {"id": m.group(1), "cmd": None, "out": None, "exit": None}
        blocks.append(cur)
    if l.startswith("``````bash") and cur is not None:
        j = i + 1; buf = []
        while not lines[j].startswith("``````"):
            buf.append(lines[j]); j += 1
        cur["cmd"] = "\n".join(buf); i = j
    elif l.startswith("``````text") and cur is not None:
        j = i + 1; buf = []
        while not lines[j].startswith("``````"):
            buf.append(lines[j]); j += 1
        cur["out"] = "\n".join(buf); i = j
    elif l.startswith("Exit code:") and cur is not None:
        cur["exit"] = l.strip()
    i += 1
for b in blocks:
    r = subprocess.run(["bash", "-c", b["cmd"]], cwd="/var/projects/toy_physics", capture_output=True, text=True)
    got = r.stdout.rstrip("\n")
    want = (b["out"] or "").rstrip("\n")
    same = (got == want)
    print(f'{b["id"]}: recorded_exit="{b["exit"]}" rerun_exit={r.returncode} recorded_lines={len(want.splitlines())} rerun_lines={len(got.splitlines())} identical={same}')
    if not same:
        gl, wl = got.splitlines(), want.splitlines()
        for k in range(max(len(gl), len(wl))):
            a = wl[k] if k < len(wl) else "<none>"
            c = gl[k] if k < len(gl) else "<none>"
            if a != c:
                print("   first diff at output line", k+1)
                print("   recorded:", a[:160]); print("   rerun:   ", c[:160]); break
