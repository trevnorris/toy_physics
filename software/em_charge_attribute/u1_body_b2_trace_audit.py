#!/usr/bin/env python3
"""Summarize traced process-tree file/network events into YAML closure evidence."""

from __future__ import annotations

import argparse
import ipaddress
import re
from pathlib import Path
from typing import Any

from u1_body_b2_common import HERE, ROOT, digest, dump_yaml, load_yaml_authenticated, manifest_absolute, parse_manifest_sets, rel_repo, require, sha256_file

CERTIFICATE_SHA256 = "650656fd2ef8a87884161825d25eced2620d8602099efbb172a960073480373c"


OPEN_RE = re.compile(r'^(?P<time>\d+\.\d+)\s+(?P<op>openat|openat2)\((?P<dirfd>[^,]+),\s+"(?P<raw>(?:\\.|[^"])*)",\s+(?P<flags>.*?)\)\s+=\s+(?P<result>-?\d+)(?:<(?P<resolved>[^>]+)>)?')
EXEC_RE = re.compile(r'^(?P<time>\d+\.\d+)\s+execve\("(?P<path>(?:\\.|[^"])*)".*?=\s+(?P<result>-?\d+)')
PATH_OP_RE = re.compile(r'^(?P<time>\d+\.\d+)\s+(?P<op>unlink|unlinkat|rename|renameat|renameat2|truncate)\((?P<body>.*)\)\s+=\s+(?P<result>-?\d+)')
CHDIR_RE = re.compile(r'^(?P<time>\d+\.\d+)\s+chdir\("(?P<path>(?:\\.|[^"])*)"\)\s+=\s+(?P<result>-?\d+)')
FCHDIR_RE = re.compile(r'^(?P<time>\d+\.\d+)\s+fchdir\((?P<fd>[^)]+)\)\s+=\s+(?P<result>-?\d+)')
CLOSE_RE = re.compile(r'^(?P<time>\d+\.\d+)\s+close\((?P<fd>\d+)(?:<[^>]+>)?\)\s+=\s+(?P<result>-?\d+)')
DUP_RE = re.compile(r'^(?P<time>\d+\.\d+)\s+(?P<op>dup|dup2|dup3)\((?P<body>.*)\)\s+=\s+(?P<result>-?\d+)')
CLONE_RE = re.compile(r'^(?P<time>\d+\.\d+)\s+(?P<op>clone|clone3|fork|vfork)\(.*\)\s+=\s+(?P<result>\d+)')
IP_RE = re.compile(r'(?:inet_addr\("(?P<v4>[0-9.]+)"\)|inet_pton\(AF_INET6,\s+"(?P<v6>[0-9a-fA-F:]+)")')
NETWORK_RE = re.compile(r'^(?P<time>\d+\.\d+)\s+(?P<op>connect|sendto|sendmmsg|sendmsg)\((?P<body>.*)\)\s+=\s+(?P<result>-?\d+)')
QUOTED_RE = re.compile(r'"((?:\\.|[^"])*)"')
FIRST_USE_RE = re.compile(r'^(?P<pid>\d+)\t(?P<kind>EXEC|LOAD|OPEN)\t(?P<sha>[0-9a-f]{64})\t(?P<device>\d+)\t(?P<inode>\d+)\t(?P<path>.*)$')


def unescape(raw: str) -> str:
    return bytes(raw, "utf-8").decode("unicode_escape")


def descriptor_path(raw: str, cwd: Path, dirfds: dict[int, Path]) -> Path | None:
    decorated = re.search(r'<([^>]+)>', raw)
    if decorated and decorated.group(1).startswith("/"):
        return Path(decorated.group(1).split("->", 1)[0]).resolve()
    token = raw.split("<", 1)[0].strip()
    if token == "AT_FDCWD":
        return cwd
    if token.lstrip("-").isdigit():
        return dirfds.get(int(token))
    return None


def resolve_operand(raw: str, cwd: Path, dirfds: dict[int, Path], dirfd: str | None = None) -> Path | None:
    value = unescape(raw)
    if value.startswith("/"):
        return Path(value).resolve()
    base = descriptor_path(dirfd, cwd, dirfds) if dirfd is not None else cwd
    return (base / value).resolve() if base is not None else None


def path_from_open(match: re.Match[str], cwd: Path, dirfds: dict[int, Path]) -> Path | None:
    resolved = match.group("resolved")
    if resolved and resolved.startswith("/"):
        return Path(resolved.split("->", 1)[0]).resolve()
    return resolve_operand(match.group("raw"), cwd, dirfds, match.group("dirfd"))


def repository_local(path: Path) -> bool:
    try:
        path.relative_to(ROOT)
        return True
    except ValueError:
        return False


def under(path: Path, prefixes: list[Path]) -> bool:
    for prefix in prefixes:
        try:
            path.relative_to(prefix)
            return True
        except ValueError:
            pass
    return False


def trace_pid(path: Path) -> int:
    match = re.search(r'\.(\d+)$', path.name)
    return int(match.group(1)) if match else 0


def parse_path_operation(op: str, body: str, cwd: Path, dirfds: dict[int, Path]) -> list[tuple[Path, bool]]:
    quoted = QUOTED_RE.findall(body)
    resolved: list[tuple[Path, bool]] = []
    if op in {"unlink", "truncate"} and quoted:
        path = resolve_operand(quoted[0], cwd, dirfds)
        if path is not None:
            resolved.append((path, False))
    elif op == "unlinkat" and quoted:
        dirfd = body.split(",", 1)[0]
        path = resolve_operand(quoted[0], cwd, dirfds, dirfd)
        if path is not None:
            resolved.append((path, False))
    elif op == "rename" and len(quoted) >= 2:
        old, new = resolve_operand(quoted[0], cwd, dirfds), resolve_operand(quoted[1], cwd, dirfds)
        if old is not None:
            resolved.append((old, False))
        if new is not None:
            resolved.append((new, True))
    elif op in {"renameat", "renameat2"} and len(quoted) >= 2:
        match = re.match(r'(?P<oldfd>[^,]+),\s*"(?:\\.|[^"])*",\s*(?P<newfd>[^,]+),\s*"', body)
        if match:
            old = resolve_operand(quoted[0], cwd, dirfds, match.group("oldfd"))
            new = resolve_operand(quoted[1], cwd, dirfds, match.group("newfd"))
            if old is not None:
                resolved.append((old, False))
            if new is not None:
                resolved.append((new, True))
    return resolved


def parse_events(trace_paths: list[Path], generated_prefixes: list[Path], initially_existing: set[Path]) -> dict[str, Any]:
    raw_rows: list[dict[str, Any]] = []
    for trace in trace_paths:
        pid = trace_pid(trace)
        for line_number, line in enumerate(trace.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            line = line.strip()
            time_match = re.match(r'^(\d+\.\d+)', line)
            if time_match:
                raw_rows.append({"time": float(time_match.group(1)), "pid": pid, "trace": trace.name, "line": line_number, "raw": line})
    raw_rows.sort(key=lambda row: (row["time"], row["trace"], row["line"]))
    states: dict[int, dict[str, Any]] = {}
    events: list[dict[str, Any]] = []
    network: list[dict[str, Any]] = []
    for raw_row in raw_rows:
            line = raw_row["raw"]
            pid = raw_row["pid"]
            state = states.setdefault(pid, {"cwd": ROOT, "dirfds": {}})
            trace, line_number = Path(raw_row["trace"]), raw_row["line"]
            match = CLONE_RE.match(line)
            if match:
                child = int(match.group("result"))
                states[child] = {"cwd": state["cwd"], "dirfds": dict(state["dirfds"])}
                continue
            match = CHDIR_RE.match(line)
            if match and int(match.group("result")) == 0:
                path = resolve_operand(match.group("path"), state["cwd"], state["dirfds"])
                if path is not None:
                    state["cwd"] = path
                continue
            match = FCHDIR_RE.match(line)
            if match and int(match.group("result")) == 0:
                path = descriptor_path(match.group("fd"), state["cwd"], state["dirfds"])
                if path is not None:
                    state["cwd"] = path
                continue
            match = CLOSE_RE.match(line)
            if match and int(match.group("result")) == 0:
                state["dirfds"].pop(int(match.group("fd")), None)
                continue
            match = DUP_RE.match(line)
            if match and int(match.group("result")) >= 0:
                fds = [int(value) for value in re.findall(r'\b\d+\b', match.group("body"))]
                if fds and fds[0] in state["dirfds"]:
                    state["dirfds"][int(match.group("result"))] = state["dirfds"][fds[0]]
                continue
            match = OPEN_RE.match(line)
            if match and int(match.group("result")) >= 0:
                path = path_from_open(match, state["cwd"], state["dirfds"])
                if path is not None and ("O_DIRECTORY" in match.group("flags") or path.is_dir()):
                    state["dirfds"][int(match.group("result"))] = path
                if path is not None and repository_local(path) and "/.git/" not in str(path):
                    flags = match.group("flags")
                    if "O_DIRECTORY" in flags:
                        continue
                    read = "O_WRONLY" not in flags
                    write = any(flag in flags for flag in ["O_WRONLY", "O_RDWR", "O_CREAT", "O_TRUNC", "O_APPEND"])
                    events.append({"time": float(match.group("time")), "pid": pid, "kind": "open", "path": path, "read": read, "write": write, "create": "O_CREAT" in flags, "truncate": "O_TRUNC" in flags, "resolved_against": str(descriptor_path(match.group("dirfd"), state["cwd"], state["dirfds"]) or state["cwd"]), "trace": trace.name, "line": line_number})
                continue
            match = EXEC_RE.match(line)
            if match and int(match.group("result")) == 0:
                path = resolve_operand(match.group("path"), state["cwd"], state["dirfds"])
                if path is not None and repository_local(path) and "/.git/" not in str(path):
                    events.append({"time": float(match.group("time")), "pid": pid, "kind": "exec", "path": path, "read": True, "write": False, "create": False, "truncate": False, "trace": trace.name, "line": line_number})
                continue
            match = PATH_OP_RE.match(line)
            if match and int(match.group("result")) == 0:
                for path, creates_destination in parse_path_operation(match.group("op"), match.group("body"), state["cwd"], state["dirfds"]):
                    if repository_local(path) and "/.git/" not in str(path):
                        events.append({"time": float(match.group("time")), "pid": pid, "kind": match.group("op"), "path": path, "read": False, "write": True, "create": creates_destination, "truncate": match.group("op") == "truncate", "resolved_against": str(state["cwd"]), "trace": trace.name, "line": line_number})
            match = NETWORK_RE.match(line)
            if match:
                addresses = []
                for ipmatch in IP_RE.finditer(match.group("body")):
                    value = ipmatch.group("v4") or ipmatch.group("v6")
                    addresses.append(value)
                network.append({"time": float(match.group("time")), "pid": pid, "operation": match.group("op"), "result": int(match.group("result")), "addresses": addresses, "trace": trace.name, "line": line_number})
    event_paths = {row["path"] for row in events}
    directory_paths = {path for path in event_paths if path.is_dir()}
    # A traced path is directory-like when it is an ancestor of another traced
    # path, including after the directory itself has been removed.  Building
    # the ancestor set is equivalent to the former all-pairs containment test
    # but scales linearly with path depth for large descendant trace sets.
    traced_parents: set[Path] = set()
    for path in event_paths:
        for candidate in path.parents:
            if candidate == ROOT.parent:
                break
            traced_parents.add(candidate)
    directory_paths.update(event_paths & traced_parents)
    events = [row for row in events if row["path"] not in directory_paths]
    events.sort(key=lambda row: (row["time"], row["trace"], row["line"]))
    first: dict[Path, dict[str, Any]] = {}
    reads, writes, generated, read_before_write = set(), set(), set(), set()
    for event in events:
        path = event["path"]
        first.setdefault(path, event)
        was_read = path in reads
        if event["read"]:
            reads.add(path)
        if event["write"]:
            writes.add(path)
            if was_read:
                read_before_write.add(path)
            if under(path, generated_prefixes) and path not in initially_existing and event["create"] and first[path] is event and not was_read:
                generated.add(path)
    wpure = {path for path in writes if under(path, generated_prefixes) and path not in read_before_write}
    preexisting_reads = reads - generated
    external_network = []
    for event in network:
        for raw in event["addresses"]:
            address = ipaddress.ip_address(raw)
            if not (address.is_loopback or address.is_unspecified):
                external_network.append({**event, "external_address": raw})
    return {
        "events": events,
        "reads": reads,
        "writes": writes,
        "generated": generated,
        "read_before_write": read_before_write,
        "W_pure": wpure,
        "preexisting_reads": preexisting_reads,
        "initially_existing": initially_existing,
        "network_events": network,
        "external_network_events": external_network,
    }


def load_first_use(trace_dir: Path) -> tuple[dict[Path, dict[str, Any]], dict[str, Any]]:
    first: dict[Path, dict[str, Any]] = {}
    environment: dict[Path, dict[str, Any]] = {}
    logs = sorted(path for path in trace_dir.rglob("*.firstuse") if path.is_file())
    require(bool(logs), "B2_FIRST_USE_TRACE", "no first-use audit logs")
    # Opening one completed audit log can append that log's binding to the
    # aggregate first-use log.  Warm every log first, then parse a second-pass
    # snapshot so authentication is independent of lexicographic log order.
    for log in logs:
        log.read_bytes()
    snapshots = {
        log: log.read_text(encoding="utf-8", errors="replace").splitlines()
        for log in logs
    }
    sequence = 0
    for log in logs:
        for line_number, line in enumerate(snapshots[log], 1):
            match = FIRST_USE_RE.match(line)
            require(match is not None, "B2_FIRST_USE_TRACE", f"malformed:{log.name}:{line_number}")
            if not match.group("path").startswith("/"):
                continue
            sequence += 1
            path = Path(match.group("path")).resolve()
            row = {"path": str(path), "sha256": match.group("sha"), "device": int(match.group("device")), "inode": int(match.group("inode")), "event_kind": match.group("kind"), "pid": int(match.group("pid")), "first_use_sequence": sequence, "audit_log": log.name, "audit_line": line_number, "binding_semantics": "regular-file bytes SHA-256 before consumer use: open/openat/openat2 descriptors before return, Wolfram regular-file arguments before kernel dispatch, and executable/loader objects before application main or dlopen return"}
            first.setdefault(path, row)
            categories = []
            raw = str(path)
            kind = match.group("kind")
            if kind == "EXEC":
                categories.append("executed_binary")
            if kind == "LOAD":
                categories.append("loaded_shared_library")
                if path.name.startswith(("ld-linux", "ld-")):
                    categories.append("loader")
            python_tree = raw.startswith("/usr/lib/python") or raw.startswith("/usr/local/lib/python") or "/.local/lib/python" in raw
            if python_tree and (kind == "OPEN" or path.suffix in {".so", ".py", ".pyc"}):
                categories.append("imported_python_module")
            if raw.startswith(("/opt/Wolfram/", "/usr/local/Wolfram/")):
                categories.append("wolfram_installation_file")
            if categories:
                current = environment.get(path)
                if current is None:
                    environment[path] = {**row, "categories": sorted(set(categories)), "event_kinds": [kind]}
                else:
                    require(current["sha256"] == row["sha256"], "B2_ENVIRONMENT_FIRST_USE", f"environment bytes changed:{path}")
                    current["categories"] = sorted(set(current["categories"]) | set(categories))
                    current["event_kinds"] = sorted(set(current["event_kinds"]) | {kind})
    environment_rows = sorted(environment.values(), key=lambda row: row["path"])
    closure = {
        "schema_version": "U1_PHASE_B2_ENVIRONMENT_CLOSURE_V1",
        "first_use_mechanism": "runner hashes successful open/openat/openat2 descriptors before return, pre-reads Wolfram regular-file arguments before kernel dispatch, and snapshots executable/loader objects before application main or dlopen return",
        "first_use_record_count": sum(len(lines) for lines in snapshots.values()),
        "entry_count": len(environment_rows),
        "entries": environment_rows,
        "category_counts": {category: sum(category in row["categories"] for row in environment_rows) for category in ["executed_binary", "loaded_shared_library", "loader", "imported_python_module", "wolfram_installation_file"]},
    }
    closure["environment_closure_digest"] = digest({"entries": [{"path": row["path"], "sha256": row["sha256"], "categories": row["categories"]} for row in environment_rows]})
    require(all(closure["category_counts"][category] > 0 for category in closure["category_counts"]), "B2_ENVIRONMENT_FIRST_USE", f"incomplete categories:{closure['category_counts']}")
    return first, closure


def serial_paths(paths: set[Path], first_use: dict[Path, dict[str, Any]], with_hash: bool = True) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(paths):
        row: dict[str, Any] = {"path": rel_repo(path)}
        if with_hash:
            authentication = first_use.get(path.resolve())
            require(authentication is not None, "B2_FIRST_USE_TRACE", f"observed dependency lacks first-use binding:{path}")
            row["sha256"] = authentication["sha256"]
            row["first_use"] = authentication
        rows.append(row)
    return rows


def run_mutation_fixture(case: str, fixture_root: Path) -> int:
    fixture_root.mkdir(parents=True, exist_ok=True)
    trace = fixture_root / "fixture.strace.41001"
    protected = fixture_root / "arbitrary_protected.bin"
    replacement = fixture_root / "replacement.bin"
    protected.write_text("protected\n", encoding="utf-8")
    replacement.write_text("replacement\n", encoding="utf-8")
    timestamp = "1000.000001"
    if case == "relative_dirfd_protected_rewrite":
        trace.write_text(
            f'{timestamp} chdir("{fixture_root}") = 0\n'
            f'1000.000002 rename("replacement.bin", "arbitrary_protected.bin") = 0\n'
            f'1000.000003 renameat(7<{fixture_root}>, "replacement.bin", 7<{fixture_root}>, "arbitrary_protected.bin") = 0\n',
            encoding="utf-8",
        )
        state = parse_events([trace], [fixture_root], {protected.resolve(), replacement.resolve()})
        require(protected.resolve() not in state["writes"], "B2_A1_PROTECTED_REWRITE", "mutation evasion: arbitrary protected relative rename resolved through cwd/dirfd")
    elif case == "preexisting_generated_laundering":
        trace.write_text(f'{timestamp} openat(AT_FDCWD<{fixture_root}>, "arbitrary_protected.bin", O_WRONLY|O_CREAT, 0666) = 8<{protected}>\n', encoding="utf-8")
        state = parse_events([trace], [fixture_root], {protected.resolve(), replacement.resolve()})
        require(protected.resolve() in state["generated"], "B2_TRACE_GENERATED_CREATE", "mutation evasion: O_CREAT on initially-existing path cannot be generated")
    else:
        raise ValueError(case)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutation-fixture", choices=["relative_dirfd_protected_rewrite", "preexisting_generated_laundering"])
    parser.add_argument("--fixture-root", type=Path)
    parser.add_argument("--trace-dir", type=Path)
    parser.add_argument("--scope", choices=["full_b1", "targeted_b1", "stage0"])
    parser.add_argument("--certificate", type=Path)
    parser.add_argument("--initial-existence", type=Path)
    parser.add_argument("--generated-prefix", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.mutation_fixture:
            require(args.fixture_root is not None, "B2_TRACE_RECIPE", "fixture root")
            return run_mutation_fixture(args.mutation_fixture, args.fixture_root.resolve())
        require(all(value is not None for value in [args.trace_dir, args.scope, args.certificate, args.initial_existence, args.output]), "B2_TRACE_RECIPE", "normal audit arguments")
        traces = sorted(path for path in args.trace_dir.rglob("*") if path.is_file() and ".strace" in path.name)
        require(bool(traces), "B2_TRACE_RECIPE", "no trace logs found")
        require(args.initial_existence.is_file(), "B2_TRACE_INITIAL_EXISTENCE", "scope-start existence snapshot missing")
        initially_existing = {Path(line).resolve() for line in args.initial_existence.read_text(encoding="utf-8").splitlines() if line.startswith("/")}
        state = parse_events(traces, [path.resolve() for path in args.generated_prefix], initially_existing)
        first_use, environment_closure = load_first_use(args.trace_dir)
        certificate, certificate_auth = load_yaml_authenticated(args.certificate, CERTIFICATE_SHA256, f"trace_audit:{args.scope}:certificate")
        sets = parse_manifest_sets(certificate)
        allowed = {manifest_absolute(raw).resolve() for rows in sets.values() for raw in rows}
        mutable = {HERE / "reports/u1_body_dynamics_results.yaml", HERE / "reports/u1_body_dynamics.md"}
        if args.scope in {"full_b1", "targeted_b1"}:
            unexpected = state["preexisting_reads"] - allowed - {path.resolve() for path in mutable}
            require(not unexpected, "B2_A1_OBS_B1_CONTAINMENT", f"unmanifested B1 reads: {[rel_repo(p) for p in sorted(unexpected)]}")
        require(not state["external_network_events"], "B2_NO_NETWORK", "external network syscall observed inside network namespace")
        artifact = {
            "schema_version": "U1_PHASE_B2_TRACE_CLOSURE_V1",
            "scope": args.scope,
            "status": "PASS",
            "first_use_authentication": [certificate_auth],
            "trace_recipe": {
                "mechanism": "strace -ff descendant recipe plus runner first-use content digests",
                "parser_syscall_coverage": ["openat", "openat2", "execve", "chdir", "fchdir", "close", "dup", "dup2", "dup3", "connect", "sendto", "sendmsg", "sendmmsg", "unlink", "unlinkat", "rename", "renameat", "renameat2", "truncate"],
                "process_tree_event_coverage": ["clone", "clone3", "fork", "vfork"],
                "network_enforcement_only_syscalls": ["socket", "bind", "listen"],
                "first_use_binding": "open/openat/openat2 descriptors are hashed before return; Wolfram regular-file arguments are hashed before kernel dispatch; executable and loaded-object snapshots occur before application main or dlopen returns",
                "process_tree": "all descendants",
                "repository_root": str(ROOT),
                "dot_git_exclusion": "git implementation metadata excluded after git-show blob content binding",
                "relative_path_rule": "per-process cwd inherited over clone/fork/vfork and updated by chdir/fchdir; every *at operand resolves against its traced dirfd (including -yy descriptor decorations)",
                "initial_existence_snapshot": str(args.initial_existence.resolve()),
                "initial_existence_count": len(initially_existing),
                "generated_output_rule": "path absent from the scope-start snapshot whose actual first file event is a successful CREATE before any read; a mere first write/O_TRUNC is never generated",
                "directory_open_rule": "O_DIRECTORY events excluded; repository-local regular-file universe only",
            },
            "trace_file_count": len(traces),
            "observed_reads": serial_paths(state["reads"], first_use),
            "observed_preexisting_reads": serial_paths(state["preexisting_reads"], first_use),
            "observed_writes": serial_paths(state["writes"], first_use, with_hash=False),
            "generated_output": serial_paths(state["generated"], first_use, with_hash=False),
            "read_before_write": serial_paths(state["read_before_write"], first_use, with_hash=False),
            "W_pure": serial_paths(state["W_pure"], first_use, with_hash=False),
            "environment_closure": environment_closure,
            "network": {
                "network_namespace": "bubblewrap --unshare-net",
                "event_count_including_local_IPC_and_netlink": len(state["network_events"]),
                "external_event_count": len(state["external_network_events"]),
                "external_events": state["external_network_events"],
            },
            "containment": {
                "equation": "Obs_B1-W_pure subset P union S union Dep_B1 union mutable_aggregate" if args.scope != "stage0" else "record_only_until_Obs_B2_manifest_is_frozen",
                "passed": True,
                "P_count": len(sets["P"]), "S_count": len(sets["S"]), "Dep_B1_count": len(sets["Dep_B1"]),
            },
            "set_sha256": digest(sorted(rel_repo(path) for path in state["preexisting_reads"])),
        }
        dump_yaml(args.output, artifact)
        print(f"B2_TRACE_AUDIT: PASS scope={args.scope} reads={len(state['reads'])} external_network=0")
        return 0
    except Exception as exc:
        print(str(exc))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
