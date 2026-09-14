"""Run one bounded native Wolfram audit with durable, source-pinned evidence.

No producer edits, automatic restart, parallel kernels, or publication occur here.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
SOURCES = [
    "mathematica/S11c_b_brane_operator_mathematica_audit.wl",
    "mathematica/S11c_c1_bulk_closure_mathematica_audit.wl",
    "mathematica/S11c_c2_N6_mathematica_audit.wl",
    "directives/S11b_SHARED_PHYSICS.md",
    *[f"directives/S11c_{part}_SHARED_PHYSICS.md" for part in ("a", "b", "c1", "c2")],
    "_measurements/S11c_wolfram_repair_audit_plan.md",
    "_measurements/S11c_wolfram_repair_audit_run.py",
]


def pin(path: Path) -> dict:
    return {"bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def main() -> None:
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("script", type=Path)
    parser.add_argument("--run-directory", type=Path, required=True)
    args = parser.parse_args()
    script = (ROOT / args.script).resolve()
    run = args.run_directory.resolve()
    run.relative_to(REPO / "_scratch" / "s11c")
    script.relative_to(ROOT)
    run.mkdir(parents=True, exist_ok=False)
    # Inventory all host processes, including other sessions, before launching.
    processes = subprocess.check_output(["ps", "-eo", "pid,comm,args"], text=True)
    (run / "processes-before.txt").write_text(processes)
    occupied = []
    for line in processes.splitlines()[1:]:
        parts = line.strip().split(None, 2)
        if len(parts) < 3 or int(parts[0]) == os.getpid():
            continue
        _, name, command = parts
        if name in {"MathKernel", "WolframKernel", "math", "Mathematica"}:
            occupied.append(line)
        if name.startswith("python") and any(token in command for token in
                ("S11c_d_", "S11c_thickness_coordinate_run_stage", "S11c_mechanical_repair_run_stage")):
            occupied.append(line)
    if occupied:
        (run / "occupied-processes.json").write_text(json.dumps(occupied, indent=2) + "\n")
        raise SystemExit("CAS process present; audit not launched. See occupied-processes.json")
    paths = list(dict.fromkeys([*SOURCES, str(script.relative_to(ROOT))]))
    source_pins = {}
    for relative in paths:
        source = ROOT / relative
        target = run / "sources" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        source_pins[relative] = pin(source)
    prior_outputs = {}
    for relative in SOURCES[:3]:
        output = ROOT / "mathematica/out" / (Path(relative).stem + ".out")
        prior_outputs[str(output.relative_to(ROOT))] = pin(output)
    manifest = {
        "status": "prepared", "head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "startedUtc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "sources": source_pins, "priorOutputs": prior_outputs,
        "command": ["math", "-script", str(run / "sources" / script.relative_to(ROOT))],
        "cwd": str(ROOT), "runDirectory": str(run),
    }
    manifest_path = run / "manifest.json"

    def save():
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    save()
    env = dict(os.environ, S11C_WOLFRAM_AUDIT_ROOT=str(run))
    started = time.monotonic()
    with (run / "full.out").open("wb") as stdout, (run / "stderr.txt").open("wb") as stderr:
        child = subprocess.Popen(manifest["command"], cwd=ROOT, env=env,
                                 stdout=stdout, stderr=stderr)
        manifest.update(status="running", childPid=child.pid)
        save()
        print(json.dumps({"pid": child.pid, "runDirectory": str(run)}), flush=True)
        returncode = child.wait()
    manifest.update(status="completed", exitCode=returncode,
                    wallSeconds=time.monotonic() - started,
                    peakRssKiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                    finishedUtc=dt.datetime.now(dt.timezone.utc).isoformat(),
                    output=pin(run / "full.out"), stderr=pin(run / "stderr.txt"),
                    sourceStability={p: pin(ROOT / p) == original for p, original in source_pins.items()},
                    priorOutputStability={p: pin(ROOT / p) == original for p, original in prior_outputs.items()})
    save()
    print(json.dumps(manifest), flush=True)
    raise SystemExit(returncode)


if __name__ == "__main__":
    main()
