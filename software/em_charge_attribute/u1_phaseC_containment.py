#!/usr/bin/env python3
"""Parse syscall evidence and enforce Phase-C network/write containment."""

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
        super().__init__(detail)
        self.assert_id = assert_id
        self.detail = detail


def require(condition: bool, assert_id: str, detail: str) -> None:
    if not condition:
        raise ContainmentFailure(assert_id, detail)


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def is_under(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def quoted_values(line: str) -> list[str]:
    quoted = re.compile(r'"((?:[^"\\]|\\.)*)"')
    return [bytes(value, "utf-8").decode("unicode_escape") for value in quoted.findall(line)]


def syscall_name(line: str) -> str | None:
    match = re.search(r"(?:^|\s)([a-zA-Z0-9_]+)\(", line)
    return match.group(1) if match else None


def path_from_raw(raw: str, cwd: Path) -> Path:
    path = Path(raw)
    if not path.is_absolute():
        path = cwd / path
    return path.resolve(strict=False)


def parse_trace(
    trace: Path, cwd: Path
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    networks: list[dict[str, Any]] = []
    inet_socket_creations: list[dict[str, Any]] = []
    writes: list[dict[str, Any]] = []
    write_open_flags = ("O_WRONLY", "O_RDWR", "O_CREAT", "O_TRUNC", "O_APPEND", "O_TMPFILE")
    path_mutators = {
        "creat",
        "mkdir",
        "mkdirat",
        "mknod",
        "mknodat",
        "rename",
        "renameat",
        "renameat2",
        "rmdir",
        "symlink",
        "symlinkat",
        "truncate",
        "unlink",
        "unlinkat",
    }
    for line_number, line in enumerate(
        trace.read_text(encoding="utf-8", errors="replace").splitlines(), start=1
    ):
        name = syscall_name(line)
        if not name:
            continue
        if name == "socket" and ("AF_INET" in line or "AF_INET6" in line):
            inet_socket_creations.append(
                {
                    "trace": str(trace),
                    "line": line_number,
                    "address_family": "AF_INET6" if "AF_INET6" in line else "AF_INET",
                    "result": line.rsplit("=", 1)[-1].strip() if "=" in line else "unknown",
                    "classification": "socket_creation_without_network_traffic",
                }
            )
        if name in {"connect", "bind", "sendto", "sendmsg", "sendmmsg"} and (
            "AF_INET" in line or "AF_INET6" in line
        ):
            networks.append(
                {
                    "trace": str(trace),
                    "line": line_number,
                    "syscall": name,
                    "address_family": "AF_INET6" if "AF_INET6" in line else "AF_INET",
                    "result": line.rsplit("=", 1)[-1].strip() if "=" in line else "unknown",
                }
            )
        values = quoted_values(line)
        targets: list[str] = []
        if name in {"open", "openat", "openat2"}:
            if values and any(flag in line for flag in write_open_flags):
                targets = [values[0]]
        elif name in path_mutators and values:
            # For rename/symlink families every quoted pathname is a write target.
            targets = values
        for raw in targets:
            if raw.startswith(("/proc/", "/dev/")):
                continue
            if name == "openat" and raw in {"uid_map", "gid_map", "setgroups"}:
                resolved_target = Path("/__phaseC_namespace_control__/proc/self") / raw
            elif name == "mkdir" and raw in {"newroot", "oldroot"}:
                resolved_target = Path("/__phaseC_namespace_control__/bwrap") / raw
            else:
                resolved_target = path_from_raw(raw, cwd)
            writes.append(
                {
                    "trace": str(trace),
                    "line": line_number,
                    "syscall": name,
                    "path": str(resolved_target),
                    "result": line.rsplit("=", 1)[-1].strip() if "=" in line else "unknown",
                }
            )
    return networks, writes, inet_socket_creations


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trace", action="append", required=True)
    parser.add_argument("--cwd", required=True)
    parser.add_argument("--allow-write-root", action="append", default=[])
    parser.add_argument("--mapped-write-root", action="append", default=[])
    parser.add_argument("--output")
    parser.add_argument("--probe-network", action="store_true")
    parser.add_argument("--probe-write", action="store_true")
    args = parser.parse_args()
    cwd = Path(args.cwd).resolve()
    allowed = [Path(value).resolve(strict=False) for value in args.allow_write_root]
    mapped_roots = []
    for value in args.mapped_write_root:
        virtual, separator, physical = value.partition("=")
        if not separator:
            parser.error("--mapped-write-root requires VIRTUAL=PHYSICAL")
        record = {
            "virtual_root": str(Path(virtual).resolve(strict=False)),
            "physical_scratch_root": str(Path(physical).resolve(strict=False)),
        }
        mapped_roots.append(record)
        allowed.append(Path(record["virtual_root"]))
    traces = [Path(value).resolve() for value in args.trace]
    try:
        networks: list[dict[str, Any]] = []
        inet_socket_creations: list[dict[str, Any]] = []
        writes: list[dict[str, Any]] = []
        for trace in traces:
            new_networks, new_writes, new_socket_creations = parse_trace(trace, cwd)
            networks.extend(new_networks)
            writes.extend(new_writes)
            inet_socket_creations.extend(new_socket_creations)
        require(
            not networks,
            "ASSERT_CONTAINMENT_NETWORK",
            f"observed {len(networks)} AF_INET/AF_INET6 syscall attempts",
        )
        forbidden = [
            row
            for row in writes
            if not any(is_under(Path(row["path"]), root) for root in allowed)
        ]
        require(
            not forbidden,
            "ASSERT_CONTAINMENT_WRITES",
            f"observed {len(forbidden)} writes outside declared roots; first={forbidden[:1]}",
        )
        require(
            not args.probe_network,
            "ASSERT_CONTAINMENT_NETWORK",
            "network probe trace escaped the network tooth",
        )
        require(
            not args.probe_write,
            "ASSERT_CONTAINMENT_WRITES",
            "write probe trace escaped the write tooth",
        )
        evidence = {
            "schema_version": "U1_PHASE_C_CONTAINMENT_EVIDENCE_V1",
            "status": "PASS",
            "measurement": "strace syscall observation, including failed attempts",
            "trace_records": [
                {"path": str(trace), "sha256": sha256_path(trace)} for trace in traces
            ],
            "allowed_write_roots": [str(root) for root in allowed],
            "scratch_bind_mapped_write_roots": mapped_roots,
            "network_attempt_count": len(networks),
            "inet_socket_creation_count": len(inet_socket_creations),
            "inet_socket_creations": inet_socket_creations,
            "write_attempt_count": len(writes),
            "write_attempts": writes,
            "forbidden_write_attempt_count": len(forbidden),
            "assertions": [
                {"assert_id": "ASSERT_CONTAINMENT_NETWORK", "status": "PASS"},
                {"assert_id": "ASSERT_CONTAINMENT_WRITES", "status": "PASS"},
            ],
        }
        require(bool(args.output), "ASSERT_CONTAINMENT_OUTPUT", "missing output")
        output = Path(args.output).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            yaml.safe_dump(evidence, sort_keys=False, allow_unicode=True, width=120),
            encoding="utf-8",
        )
        print(
            f"PHASEC_CONTAINMENT_PASS traces={len(traces)} writes={len(writes)} network=0"
        )
        return 0
    except ContainmentFailure as failure:
        print(f"ASSERTION_FAILED {failure.assert_id}: {failure.detail}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
