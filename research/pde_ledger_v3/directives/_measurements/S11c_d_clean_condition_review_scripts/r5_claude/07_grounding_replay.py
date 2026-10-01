#!/usr/bin/env python3
"""Replay every `$ <command>` block in both grounding files against the current working tree and report
blocks whose literal recorded output differs from the command's present output.  Mechanical only."""
import re, subprocess, sys
ROOT = "/var/projects/toy_physics"
FILES = [
    "research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition.md",
    "research/pde_ledger_v3/directives/_measurements/S11c_d_zinvariant_operator_blocks_directive.md",
]
SAFE = re.compile(r"^(sed -n|grep|wc|git log|git show|ls|awk|head|tail|cat|rg|nl)\b")
for f in FILES:
    text = open(f"{ROOT}/{f}", encoding="utf-8").read()
    blocks = re.findall(r"^(`{3,})\n\$ (.+?)\n(.*?)^\1$", text, flags=re.M | re.S)
    n_run = n_same = n_diff = n_skip = 0
    print(f"=== {f}: {len(blocks)} command blocks")
    for fence, cmd, out in blocks:
        if not SAFE.match(cmd) or "python" in cmd:
            n_skip += 1
            print(f"SKIP_NOT_A_LOOKUP {cmd[:110]}")
            continue
        r = subprocess.run(cmd, shell=True, cwd=ROOT, capture_output=True, text=True, timeout=120)
        now = r.stdout
        n_run += 1
        if now.rstrip("\n") == out.rstrip("\n"):
            n_same += 1
        else:
            n_diff += 1
            print(f"DIFFERS {cmd}")
            a, b = out.rstrip("\n").splitlines(), now.rstrip("\n").splitlines()
            for i in range(max(len(a), len(b))):
                x = a[i] if i < len(a) else "<none>"
                y = b[i] if i < len(b) else "<none>"
                if x != y:
                    print(f"   recorded[{i}]: {x[:160]}")
                    print(f"   current [{i}]: {y[:160]}")
                    break
    print(f"RUN {n_run} SAME {n_same} DIFFERS {n_diff} SKIPPED {n_skip}")
