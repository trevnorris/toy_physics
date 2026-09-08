#!/usr/bin/env python3
"""Serialize the WL harness workers and record their literal transcript.

Manifest (orchestrator-owned):
K_A: faceSources flux assignment; drop lambdaAResponse affinity.
K_T: faceSources virtualWork product assignment; replace the whole product by 0.
K_W: pressureField[-1]; replace the applied lower pressure by the applied upper.
The executable source patches and extractor self-tests live in the WL file.

1. The harness may **PRINT** computed objects (carriers, diffs). It may ⛔ NOT state conclusions — no `PASS`, no
   verdict, no "bites"/"fails". Interpretation is the orchestrator's, from the printed triples.
2. **Print operand and residual, then guard.** Emit `baseline`, `corrupted`, and their `diff`; a residual asserted
   zero/nonzero carries no information.
3. Interpretation belongs to the review / step record, ⛔ not the script.
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    root = Path(__file__).resolve().parents[3]
    ledger = root / "research/pde_ledger_v3"
    engine = ledger / "mathematica/S11c_b_brane_operator_mathematica_audit.wl"
    harness = ledger / "mathematica/S11c_b_carrier_ablation_harness.wl"
    directive = ledger / "directives/S11c_b_carrier_ablation_harness_directive.md"
    record = ledger / "_measurements/S11c_b_carrier_ablation_harness_wl.md"
    runner = shutil.which("wolframscript")
    timeout = shutil.which("timeout")
    if runner is None or timeout is None:
        raise SystemExit("Required executables: wolframscript and timeout")
    sources = {path: digest(path) for path in (engine, harness, Path(__file__).resolve(), directive)}
    transcript: list[bytes] = []

    def emit(data: bytes) -> None:
        transcript.append(data)
        sys.stdout.buffer.write(data)
        sys.stdout.buffer.flush()

    modes = ("BASELINE", "UNABLATED_COPY", "K_A", "K_T", "K_W", "DEAD_PATH",
             "RESCALE_K_A", "RESCALE_K_T", "RESCALE_K_W", "REPORT")
    status = 0
    with tempfile.TemporaryDirectory(prefix="s11cb_carrier_wl_") as temporary:
        for mode in modes:
            environment = os.environ.copy()
            configured = {
                "S11CB_CARRIER_MODE": mode,
                "S11CB_CARRIER_WORKDIR": temporary,
                "S11CB_CARRIER_ENGINE": str(engine),
            }
            environment.update(configured)
            command = [timeout, "600", runner, "-file", str(harness)]
            invocation = " ".join(f"{key}={shlex.quote(value)}" for key, value in configured.items())
            emit(("$ " + invocation + " " + shlex.join(command) + "\n").encode())
            process = subprocess.Popen(command, cwd=root, env=environment,
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            assert process.stdout is not None
            for chunk in iter(process.stdout.readline, b""):
                emit(chunk)
            status = process.wait()
            emit(f"PROCESS_EXIT_CODE: {status}\n".encode())
            if status:
                break
    for path, original_digest in sources.items():
        current_digest = digest(path)
        if current_digest != original_digest:
            emit((f"SOURCE_DRIFT: {path.relative_to(root)} "
                  f"{original_digest} {current_digest}\n").encode())
            status = 2
    output = b"".join(transcript)
    invocation = shlex.join([sys.executable, *sys.argv])
    lines = ["# S11c-b WL carrier ablation harness run", "",
             f"Working directory: `{root}`", "", "Exact invocation:", "",
             "```sh", invocation, "```", "", "Source SHA-256 digests:", ""]
    lines.extend(f"- `{path.relative_to(root)}`: `{value}`" for path, value in sources.items())
    lines.extend(["", f"Output SHA-256: `{hashlib.sha256(output).hexdigest()}`", "",
                  f"Process exit code: `{status}`", "", "Literal transcript:", "", "```text"])
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_bytes(("\n".join(lines) + "\n").encode() + output + b"```\n")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
