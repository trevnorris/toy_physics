#!/usr/bin/env python3
"""Execute a manifest-authenticated source through its held file descriptor."""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
from pathlib import Path

import yaml

from u1_body_b2_common import ROOT, UniqueKeyLoader, dump_yaml, load_yaml, require


def held_blob(path: Path) -> tuple[int, bytes, os.stat_result]:
    fd = os.open(path, os.O_RDONLY)
    chunks = []
    while True:
        block = os.read(fd, 1 << 20)
        if not block:
            break
        chunks.append(block)
    os.lseek(fd, 0, os.SEEK_SET)
    return fd, b"".join(chunks), os.fstat(fd)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--stage0", type=Path, required=True)
    p.add_argument("--stage0-contract-digest", required=True)
    p.add_argument("--consumer", required=True)
    p.add_argument("--target", type=Path, required=True)
    p.add_argument("--evidence", type=Path, required=True)
    p.add_argument("--producer-record", type=Path, action="append", default=[])
    p.add_argument("command", nargs=argparse.REMAINDER)
    a = p.parse_args()
    stage_fd = target_fd = -1
    bound_fds: list[int] = []
    try:
        stage_fd, stage_blob, stage_stat = held_blob(a.stage0.resolve())
        stage_sha = hashlib.sha256(stage_blob).hexdigest()
        require(stage_sha == a.stage0_contract_digest, "B2_A1_STAGE0_CONTAINER", a.consumer)
        stage0 = yaml.load(stage_blob.decode("utf-8"), Loader=UniqueKeyLoader)
        target = a.target.resolve()
        target_fd, target_blob, target_stat = held_blob(target)
        target_sha = hashlib.sha256(target_blob).hexdigest()
        rel = target.relative_to(ROOT).as_posix()
        expected = stage0["observation_contract"]["Obs_B2_manifest"].get(rel)
        require(expected == target_sha, "B2_A1_OBS_B2_FIRST_USE", f"{a.consumer}:{rel}")
        command = list(a.command)
        if command and command[0] == "--":
            command = command[1:]
        require(bool(command), "B2_A1_OBS_B2_FIRST_USE", "missing authenticated command")
        held_path = f"/proc/self/fd/{target_fd}"
        replaced = [held_path if Path(token).resolve() == target else token for token in command]
        require(replaced != command, "B2_A1_OBS_B2_FIRST_USE", f"target absent from command:{target}")
        bound = []
        for record_path in a.producer_record:
            record = load_yaml(record_path.resolve())
            require(record.get("status") == "PASS", "B2_A1_PRODUCER_DIGEST", f"{a.consumer}:{record_path}")
            for row in record["outputs"]:
                produced = Path(row["path"]).resolve()
                fd, blob, stat = held_blob(produced)
                observed = hashlib.sha256(blob).hexdigest()
                require(observed == row["sha256"], "B2_A1_PRODUCER_DIGEST", f"{a.consumer}:{produced}")
                bound_fds.append(fd)
                fd_path = f"/proc/self/fd/{fd}"
                before = list(replaced)
                replaced = [fd_path if Path(token).resolve() == produced else token for token in replaced]
                bound.append({"path": str(produced), "producer": row["producer"], "expected_sha256": row["sha256"], "consumed_fd_sha256": observed, "device": stat.st_dev, "inode": stat.st_ino, "passed_as_immutable_blob": fd_path, "appeared_in_command": before != replaced})
        if a.producer_record:
            require(any(row["appeared_in_command"] for row in bound), "B2_A1_PRODUCER_DIGEST", f"{a.consumer}:no producer-bound input in command")
        proc = subprocess.run(replaced, cwd=target.parent, check=False, pass_fds=tuple([target_fd, *bound_fds]))
        # The child consumed immutable descriptors.  Also fail if a protected
        # path was swapped during the run, rather than silently tolerating it.
        require(hashlib.sha256(target.read_bytes()).hexdigest() == target_sha, "B2_A1_PROTECTED_REWRITE", f"{a.consumer}:{target}")
        for row in bound:
            require(hashlib.sha256(Path(row["path"]).read_bytes()).hexdigest() == row["expected_sha256"], "B2_A1_PROTECTED_REWRITE", f"{a.consumer}:{row['path']}")
        dump_yaml(a.evidence, {"schema_version": "U1_PHASE_B2_FIRST_USE_EXEC_V2", "status": "PASS" if proc.returncode == 0 else "CHILD_EXIT", "consumer": a.consumer, "stage0": {"consumed_fd_sha256": stage_sha, "device": stage_stat.st_dev, "inode": stage_stat.st_ino}, "target": {"path": rel, "expected_sha256": expected, "consumed_fd_sha256": target_sha, "device": target_stat.st_dev, "inode": target_stat.st_ino, "executed_via_held_descriptor": held_path}, "producer_bound_inputs": bound, "protected_path_rewrite_check": "PASS", "child_exit_code": proc.returncode})
        return proc.returncode
    except Exception as exc:
        print(str(exc))
        return 1
    finally:
        if stage_fd >= 0:
            os.close(stage_fd)
        if target_fd >= 0:
            os.close(target_fd)
        for fd in bound_fds:
            os.close(fd)


if __name__ == "__main__":
    raise SystemExit(main())
