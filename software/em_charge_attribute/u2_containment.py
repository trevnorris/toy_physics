#!/usr/bin/env python3
"""Parse syscall evidence and enforce U2 whole-process network/write containment."""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path
from typing import Any

import yaml


class ContainmentFailure(RuntimeError):
    def __init__(self, assert_id: str, detail: str):
        super().__init__(detail); self.assert_id = assert_id; self.detail = detail


def require(condition: bool, assert_id: str, detail: str) -> None:
    if not condition: raise ContainmentFailure(assert_id, detail)


def sha256_path(path: Path) -> str: return hashlib.sha256(path.read_bytes()).hexdigest()


def is_under(path: Path, root: Path) -> bool:
    try: path.relative_to(root); return True
    except ValueError: return False


def quoted_values(line: str) -> list[str]:
    return [bytes(value, "utf-8").decode("unicode_escape") for value in re.compile(r'"((?:[^"\\]|\\.)*)"').findall(line)]


def syscall_name(line: str) -> str | None:
    match = re.search(r"(?:^|\s)([a-zA-Z0-9_]+)\(", line); return match.group(1) if match else None


def parse_trace(trace: Path, cwd: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    networks = []; sockets = []; writes = []
    write_flags = ("O_WRONLY", "O_RDWR", "O_CREAT", "O_TRUNC", "O_APPEND", "O_TMPFILE")
    mutators = {"creat", "mkdir", "mkdirat", "mknod", "mknodat", "rename", "renameat", "renameat2", "rmdir", "symlink", "symlinkat", "truncate", "unlink", "unlinkat"}
    for number, line in enumerate(trace.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        name = syscall_name(line)
        if not name: continue
        result = line.rsplit("=", 1)[-1].strip() if "=" in line else "unknown"
        if name == "socket" and ("AF_INET" in line or "AF_INET6" in line):
            sockets.append({"trace": str(trace), "line": number, "result": result})
        if name in {"connect", "bind", "sendto", "sendmsg", "sendmmsg"} and ("AF_INET" in line or "AF_INET6" in line):
            networks.append({"trace": str(trace), "line": number, "syscall": name, "result": result})
        values = quoted_values(line); targets = []
        if name in {"open", "openat", "openat2"} and values and any(flag in line for flag in write_flags): targets = [values[0]]
        elif name in mutators and values: targets = values
        for raw in targets:
            if raw.startswith(("/proc/", "/dev/")): continue
            if name == "openat" and raw in {"uid_map", "gid_map", "setgroups"}:
                target = Path("/__u2_namespace_control__/proc/self") / raw
            elif name == "mkdir" and raw in {"newroot", "oldroot"}:
                target = Path("/__u2_namespace_control__/bwrap") / raw
            else:
                target = Path(raw); target = target if target.is_absolute() else cwd / target; target = target.resolve(strict=False)
            writes.append({"trace": str(trace), "line": number, "syscall": name, "path": str(target), "result": result})
    return networks, writes, sockets


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("--trace", action="append", required=True)
    parser.add_argument("--cwd", required=True); parser.add_argument("--allow-write-root", action="append", default=[])
    parser.add_argument("--mapped-write-root", action="append", default=[]); parser.add_argument("--output")
    parser.add_argument("--probe-network", action="store_true"); parser.add_argument("--probe-write", action="store_true")
    args = parser.parse_args(); cwd = Path(args.cwd).resolve()
    allowed = [Path(value).resolve(strict=False) for value in args.allow_write_root]; mapped = []
    for value in args.mapped_write_root:
        virtual, separator, physical = value.partition("=")
        if not separator: parser.error("mapped root requires VIRTUAL=PHYSICAL")
        mapped.append({"virtual_root": str(Path(virtual).resolve(strict=False)), "physical_scratch_root": str(Path(physical).resolve(strict=False))})
        allowed.append(Path(virtual).resolve(strict=False))
    traces = [Path(value).resolve() for value in args.trace]
    try:
        networks = []; writes = []; sockets = []
        for trace in traces:
            new_networks, new_writes, new_sockets = parse_trace(trace, cwd)
            networks.extend(new_networks); writes.extend(new_writes); sockets.extend(new_sockets)
        require(not networks and not args.probe_network, "ASSERT_CONTAINMENT_NETWORK", f"observed {len(networks)} network attempts")
        forbidden = [row for row in writes if not any(is_under(Path(row["path"]), root) for root in allowed)]
        require(not forbidden and not args.probe_write, "ASSERT_CONTAINMENT_WRITES", f"observed {len(forbidden)} forbidden writes; first={forbidden[:1]}")
        evidence = {
            "schema_version": "U2_CONTAINMENT_EVIDENCE_V1", "status": "PASS",
            "measurement": "strace whole descendant tree including failed attempts",
            "trace_records": [{"path": str(trace), "sha256": sha256_path(trace)} for trace in traces],
            "allowed_write_roots": [str(root) for root in allowed], "mapped_write_roots": mapped,
            "network_attempt_count": len(networks), "inet_socket_creation_count": len(sockets),
            "inet_socket_creations": sockets, "write_attempt_count": len(writes), "write_attempts": writes,
            "forbidden_write_attempt_count": len(forbidden),
            "assertions": [{"assert_id": "ASSERT_CONTAINMENT_NETWORK", "status": "PASS"}, {"assert_id": "ASSERT_CONTAINMENT_WRITES", "status": "PASS"}],
        }
        require(bool(args.output), "ASSERT_CONTAINMENT_OUTPUT", "missing output")
        output = Path(args.output); output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(yaml.safe_dump(evidence, sort_keys=False, allow_unicode=True, width=140), encoding="utf-8")
        print(f"U2_CONTAINMENT_PASS traces={len(traces)} writes={len(writes)} network=0"); return 0
    except ContainmentFailure as failure:
        print(f"ASSERTION_FAILED {failure.assert_id}: {failure.detail}", file=sys.stderr); return 1


if __name__ == "__main__": raise SystemExit(main())
