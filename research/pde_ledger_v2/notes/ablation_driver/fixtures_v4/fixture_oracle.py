#!/usr/bin/env python3
"""Independent readers and relational discriminators for the v4 fixtures."""

from __future__ import annotations

import copy
import hashlib
import json
import pathlib
import re
from collections import Counter
from typing import Any, Iterable


RESULT_COLUMNS = [
    "schema",
    "list_sha256",
    "row_number",
    "row_id",
    "target",
    "axis",
    "key",
    "record",
    "line",
    "old_text",
    "new_text",
    "include_extra",
    "match_count",
    "applied",
    "target_entry_sha256",
    "target_mutant_sha256",
    "producer_exit",
    "producer_tally",
    "first_failed_assertion",
    "artifacts_emitted",
    "moved_values",
    "checker_exits",
    "checker_verdicts",
    "captures",
    "outcome",
]
RESTORE_COLUMNS = [
    "event",
    "operation",
    "kind",
    "name",
    "entry_state",
    "entry_sha256",
    "pre_restore_state",
    "pre_restore_sha256",
    "restored_state",
    "restored_sha256",
    "restored",
]
OUTCOMES = {
    "NOT_APPLIED_ZERO",
    "NOT_APPLIED_MULTIPLE",
    "PRODUCER_NONZERO",
    "CHECKER_NONZERO",
    "ALL_ZERO",
}
INAPPLICABLE_EXECUTION_FIELDS = [
    "target_mutant_sha256",
    "producer_exit",
    "producer_tally",
    "first_failed_assertion",
    "artifacts_emitted",
    "moved_values",
    "checker_exits",
    "checker_verdicts",
    "captures",
]
TALLY_RE = re.compile(
    r"^TALLY (?P<label>.+?): (?P<pass>[0-9]+) pass \+ "
    r"(?P<fail>[0-9]+) fail = (?P<total>[0-9]+) checks$"
)
RESULT_FIELD_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
SHA_RE = re.compile(r"[0-9a-f]{64}")
EVIDENCE_MD_SHA256 = "bee409ce92a59d48fa1fd85eabf9c80b6fc228b59afadf6342931948c332ba95"
A7_LIST_SHA256 = "081ad6ab482adff9ac0be5c3d9355685ec07e52df051878becfd97a1264e762b"
A7_RESULTS_SHA256 = "28807ff0910d98c59526539423f566962dd1faeef79a1ea70110f0158873cd13"
A7_PROJECTION_SHA256 = "e498caf729bb6570610423b50859a12bc04c7ac810004ebd3078dc41535a408a"
LEGACY_RESULT_COLUMNS = [
    "axis",
    "key",
    "record",
    "stage_exit",
    "pass_count",
    "fail_count",
    "first_fail",
    "sidecar_written",
    "record_moved",
    "emitted_value",
    "cmp_exit",
    "cmp_status",
    "mismatch_names",
]


class FixtureFailure(AssertionError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FixtureFailure(message)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


EntrySnapshot = dict[tuple[str, str], bytes | None]


def capture_entry_snapshot(
    project_root: pathlib.Path,
    config: dict[str, Any],
) -> EntrySnapshot:
    """Capture restoration-set bytes before the driver is invoked."""
    snapshot: EntrySnapshot = {}
    for kind, specifications in (
        ("target", config["targets"]),
        ("artifact", config["artifacts"]),
    ):
        for specification in specifications:
            path = project_root / specification["path"]
            if path.is_file():
                snapshot[(kind, specification["name"])] = path.read_bytes()
            else:
                require(kind == "artifact" and not path.exists(), f"invalid entry path {path}")
                snapshot[(kind, specification["name"])] = None
    return snapshot


def canonical_json(value: Any) -> str:
    return json.dumps(
        value, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":")
    )


def parse_standalone_jcs(data: bytes, label: str) -> Any:
    require(data.endswith(b"\n") and not data.endswith(b"\n\n"), f"{label}: one LF required")
    require(b"\r" not in data and data[:3] != b"\xef\xbb\xbf", f"{label}: text encoding")
    try:
        text = data[:-1].decode("utf-8")
        value = json.loads(text)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise FixtureFailure(f"{label}: invalid JSON: {exc}") from exc
    require(canonical_json(value) == text, f"{label}: JSON is not canonical")
    return value


def parse_jcs_cell(cell: str, label: str) -> Any:
    try:
        value = json.loads(cell)
    except json.JSONDecodeError as exc:
        raise FixtureFailure(f"{label}: invalid JSON cell: {exc}") from exc
    require(canonical_json(value) == cell, f"{label}: JSON cell is not canonical")
    return value


def read_tsv_bytes(data: bytes, label: str) -> tuple[list[str], list[list[str]]]:
    require(data.endswith(b"\n"), f"{label}: final LF required")
    require(b"\r" not in data and data[:3] != b"\xef\xbb\xbf", f"{label}: encoding")
    try:
        text = data.decode("utf-8")
    except UnicodeError as exc:
        raise FixtureFailure(f"{label}: invalid UTF-8: {exc}") from exc
    rows = [line.split("\t") for line in text[:-1].split("\n")]
    require(bool(rows), f"{label}: empty TSV")
    width = len(rows[0])
    require(width > 0 and all(len(row) == width for row in rows), f"{label}: row width")
    require(all(cell != "" for cell in rows[0]), f"{label}: empty header")
    require(len(set(rows[0])) == width, f"{label}: duplicate header")
    return rows[0], rows[1:]


def read_include(data: bytes) -> tuple[list[str], list[dict[str, str]]]:
    header, body = read_tsv_bytes(data, "include_list.tsv")
    required = {"axis", "key", "record", "old_text", "new_text"}
    require(required <= set(header), "include_list.tsv: required header")
    rows = [dict(zip(header, fields, strict=True)) for fields in body]
    return header, rows


def occurrence_count(text: str, old: str, line: int | None) -> int:
    if line is None:
        return text.count(old)
    logical = text.splitlines(keepends=True)
    if line > len(logical):
        return 0
    return logical[line - 1].count(old)


def exact_mutant_bytes(entry: bytes, old: str, new: str, line: int | None) -> bytes:
    try:
        text = entry.decode("utf-8")
    except UnicodeError as exc:
        raise FixtureFailure("target entry is not UTF-8") from exc
    require(occurrence_count(text, old, line) == 1, "exact mutant requires one match")
    if line is None:
        mutant = text.replace(old, new, 1)
    else:
        logical = text.splitlines(keepends=True)
        logical[line - 1] = logical[line - 1].replace(old, new, 1)
        mutant = "".join(logical)
    return mutant.encode("utf-8")


def parse_result_line(line: str) -> str:
    if not line.startswith("RESULT|"):
        raise FixtureFailure("not a RESULT line")
    fields: dict[str, str] = {}
    for raw in line.split("|")[1:]:
        require("=" in raw, "malformed RESULT field")
        name, value = raw.split("=", 1)
        require(bool(RESULT_FIELD_RE.fullmatch(name)), "malformed RESULT name")
        require(name not in fields and value != "" and "|" not in value, "malformed RESULT value")
        fields[name] = value
    require("status" in fields, "RESULT has no status")
    return fields["status"]


def checker_observation(data: bytes, report_format: str, name: str, exit_code: int) -> dict[str, Any]:
    if report_format == "exit-only":
        return {
            "mismatch_names": [],
            "name": name,
            "status": "EXIT_ZERO" if exit_code == 0 else "EXIT_NONZERO",
        }
    require(report_format == "ledger-result-v1", "unsupported checker report format")
    try:
        lines = data.decode("utf-8").splitlines()
    except UnicodeError as exc:
        raise FixtureFailure(f"checker {name}: stdout is not UTF-8") from exc
    result_lines = [line for line in lines if line.startswith("RESULT|")]
    require(len(result_lines) <= 1, f"checker {name}: multiple RESULT lines")
    mismatches = [
        line[len("MISMATCH ") :].split(":", 1)[0]
        for line in lines
        if line.startswith("MISMATCH ") and ":" in line and line[len("MISMATCH ") :].split(":", 1)[0]
    ]
    status = "NO_RESULT_LINE"
    if result_lines:
        status = parse_result_line(result_lines[0])
    return {"mismatch_names": mismatches, "name": name, "status": status}


def producer_observation(data: bytes, report_format: str) -> tuple[Any, Any]:
    if report_format == "exit-only":
        return "-", "-"
    require(report_format == "ledger-tally-v1", "unsupported producer report format")
    try:
        lines = data.decode("utf-8").splitlines()
    except UnicodeError as exc:
        raise FixtureFailure("producer stdout is not UTF-8") from exc
    tally_positions = [
        (position, match)
        for position, line in enumerate(lines)
        if (match := TALLY_RE.fullmatch(line)) is not None
    ]
    require(len(tally_positions) <= 1, "multiple producer tallies")
    malformed = [line for line in lines if line.startswith("TALLY ") and not TALLY_RE.fullmatch(line)]
    require(not malformed, "malformed producer tally")
    if not tally_positions:
        return "-", "-"
    tally_position, match = tally_positions[0]
    passed = int(match.group("pass"))
    failed = int(match.group("fail"))
    total = int(match.group("total"))
    require(total == passed + failed, "producer tally arithmetic")
    first = next(
        (line for line in lines[:tally_position] if re.match(r"^FAIL +", line)),
        None,
    )
    require((failed == 0 and first is None) or (failed > 0 and first is not None), "first failure grammar")
    return {"fail": failed, "pass": passed}, first if first is not None else "-"


def dimension_records(data: bytes) -> dict[str, str]:
    try:
        lines = data.decode("utf-8").splitlines()
    except UnicodeError as exc:
        raise FixtureFailure("dimension artifact is not UTF-8") from exc
    records: dict[str, str] = {}
    for line in lines:
        if not line.startswith("DIM|"):
            continue
        parts = line.split("|")
        require(len(parts) == 4, "malformed DIM record")
        fields: dict[str, str] = {}
        for part in parts[1:]:
            require("=" in part, "malformed DIM field")
            name, value = part.split("=", 1)
            require(name not in fields and value != "", "malformed DIM value")
            fields[name] = value
        require(set(fields) == {"axes", "name", "exponents"}, "DIM field set")
        require(fields["name"] not in records, "duplicate DIM name")
        records[fields["name"]] = fields["exponents"]
    return records


def outcome_for(row: dict[str, Any]) -> str:
    match_count = int(row["match_count"])
    applied = row["applied"]
    require(applied == (match_count == 1), "applied does not match match_count")
    if match_count == 0:
        return "NOT_APPLIED_ZERO"
    if match_count >= 2:
        return "NOT_APPLIED_MULTIPLE"
    producer_exit = int(row["producer_exit"])
    if producer_exit != 0:
        return "PRODUCER_NONZERO"
    checker_exits = row["checker_exits"]
    require(bool(checker_exits), "applied row has no checker")
    if any(int(item["exit"]) != 0 for item in checker_exits):
        return "CHECKER_NONZERO"
    return "ALL_ZERO"


def validate_outcome(row: dict[str, Any]) -> None:
    require(row["outcome"] in OUTCOMES, "unknown outcome")
    if not row["applied"]:
        for field in INAPPLICABLE_EXECUTION_FIELDS:
            require(row[field] == "-", f"non-applied row has {field}")
    required = outcome_for(row)
    require(row["outcome"] == required, f"outcome {row['outcome']} is unsupported")


def _target_map(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {item["name"]: item for item in config["targets"]}


def _artifact_map(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {item["name"]: item for item in config["artifacts"]}


def parse_and_audit_results(
    project_root: pathlib.Path,
    evidence: pathlib.Path,
    config: dict[str, Any],
    include_data: bytes,
    entry_snapshot: EntrySnapshot,
    *,
    result_data: bytes | None = None,
    require_complete: bool = True,
) -> list[dict[str, Any]]:
    header, include_rows = read_include(include_data)
    if result_data is None:
        result_data = (evidence / "results.tsv").read_bytes()
    result_header, result_body = read_tsv_bytes(result_data, "results.tsv")
    require(result_header == RESULT_COLUMNS, "results.tsv: exact header")
    if require_complete:
        require(len(result_body) == len(include_rows), "results.tsv: one row per include row")
    else:
        require(len(result_body) <= len(include_rows), "results.tsv: prefix exceeds include-list")
    require(
        all(all(cell != "" for cell in cells) for cells in result_body),
        "results.tsv: empty data cell",
    )
    list_sha = sha256_bytes(include_data)
    targets = _target_map(config)
    artifacts = _artifact_map(config)
    target_entry: dict[str, bytes] = {}
    for name in targets:
        data = entry_snapshot[("target", name)]
        require(data is not None, f"target {name}: fixture entry snapshot is absent")
        target_entry[name] = data
    artifact_entry = {
        name: entry_snapshot[("artifact", name)]
        for name in artifacts
    }
    decoded: list[dict[str, Any]] = []
    audited_include_rows = include_rows[: len(result_body)]
    for position, (include, cells) in enumerate(
        zip(audited_include_rows, result_body, strict=True),
        start=1,
    ):
        raw = dict(zip(RESULT_COLUMNS, cells, strict=True))
        row: dict[str, Any] = copy.deepcopy(raw)
        resolved = include.get("target")
        if resolved is None:
            require(len(targets) == 1, "list without target in multi-target config")
            resolved = next(iter(targets))
        expected_id = f"{resolved}:{include['axis']}:{include['key']}"
        require(raw["schema"] == "ablation-result-v1", "result schema")
        require(raw["list_sha256"] == list_sha, "result list digest")
        require(raw["row_number"] == str(position), "result row number/order")
        require(raw["row_id"] == expected_id, "result row identity")
        require(raw["target"] == resolved, "result target")
        require(raw["axis"] == include["axis"] and raw["key"] == include["key"], "result key")
        for name in ("record", "old_text", "new_text"):
            row[name] = parse_jcs_cell(raw[name], f"{expected_id}:{name}")
            require(row[name] == include[name], f"{expected_id}: list text preservation")
        line_value = include.get("line", "")
        line = int(line_value) if line_value else None
        require(raw["line"] == (str(line) if line is not None else "-"), f"{expected_id}: line")
        extras = {
            name: include[name]
            for name in header
            if name not in {"axis", "key", "record", "old_text", "new_text", "target", "line"}
        }
        row["include_extra"] = parse_jcs_cell(raw["include_extra"], f"{expected_id}:include_extra")
        require(row["include_extra"] == extras, f"{expected_id}: include extras")
        entry_text = target_entry[resolved].decode("utf-8")
        count = occurrence_count(entry_text, include["old_text"], line)
        require(raw["match_count"] == str(count), f"{expected_id}: match count")
        row["match_count"] = count
        require(raw["applied"] in {"true", "false"}, f"{expected_id}: applied spelling")
        row["applied"] = raw["applied"] == "true"
        require(row["applied"] == (count == 1), f"{expected_id}: applied invariant")
        entry_digest = sha256_bytes(target_entry[resolved])
        require(raw["target_entry_sha256"] == entry_digest, f"{expected_id}: entry digest")
        if not row["applied"]:
            for field in INAPPLICABLE_EXECUTION_FIELDS:
                require(raw[field] == "-", f"{expected_id}: {field} must be inapplicable")
            decoded.append(row)
            validate_outcome(row)
            continue
        exact_mutant = exact_mutant_bytes(
            target_entry[resolved],
            include["old_text"],
            include["new_text"],
            line,
        )
        require(
            raw["target_mutant_sha256"] == sha256_bytes(exact_mutant),
            f"{expected_id}: exact C-3 mutant digest",
        )
        row["producer_exit"] = int(raw["producer_exit"])
        require(0 <= row["producer_exit"] <= 255, f"{expected_id}: producer exit")
        for field in (
            "artifacts_emitted",
            "moved_values",
            "checker_exits",
            "checker_verdicts",
            "captures",
        ):
            row[field] = parse_jcs_cell(raw[field], f"{expected_id}:{field}")
        emitted = row["artifacts_emitted"]
        require(
            emitted == [item["name"] for item in config["artifacts"] if item["name"] in emitted],
            f"{expected_id}: artifact order/identity",
        )
        checker_exits = row["checker_exits"]
        expected_checker_names = [item["name"] for item in config["checkers"]]
        require(
            [item.get("name") for item in checker_exits] == expected_checker_names,
            f"{expected_id}: checker exit sequence",
        )
        for item in checker_exits:
            require(set(item) == {"exit", "name"} and isinstance(item["exit"], int), "checker exit grammar")
            require(0 <= item["exit"] <= 255, "checker exit range")
        expected_capture_specs: list[tuple[str, str]] = [
            ("producer.stdout", f"captures/{position:06d}/producer.stdout"),
            ("producer.stderr", f"captures/{position:06d}/producer.stderr"),
        ]
        for checker_position, checker in enumerate(config["checkers"], start=1):
            expected_capture_specs.extend(
                [
                    (
                        f"checker.{checker['name']}.stdout",
                        f"captures/{position:06d}/checker-{checker_position:02d}-{checker['name']}.stdout",
                    ),
                    (
                        f"checker.{checker['name']}.stderr",
                        f"captures/{position:06d}/checker-{checker_position:02d}-{checker['name']}.stderr",
                    ),
                ]
            )
        for artifact in config["artifacts"]:
            if artifact["name"] in emitted:
                expected_capture_specs.append(
                    (
                        f"artifact.{artifact['name']}",
                        f"captures/{position:06d}/artifact-{artifact['name']}",
                    )
                )
        require(len(row["captures"]) == len(expected_capture_specs), f"{expected_id}: capture count")
        capture_data: dict[str, bytes] = {}
        for reference, (role, path) in zip(row["captures"], expected_capture_specs, strict=True):
            require(set(reference) == {"bytes", "path", "role", "sha256"}, "capture reference fields")
            require(reference["role"] == role and reference["path"] == path, "capture role/path")
            capture_path = evidence / path
            require(capture_path.is_file(), f"missing claimed capture {path}")
            data = capture_path.read_bytes()
            require(reference["bytes"] == len(data), f"capture size {path}")
            require(reference["sha256"] == sha256_bytes(data), f"capture digest {path}")
            capture_data[role] = data
        tally, first = producer_observation(
            capture_data["producer.stdout"], config["producer"]["report_format"]
        )
        row["producer_tally"] = (
            "-" if raw["producer_tally"] == "-" else parse_jcs_cell(raw["producer_tally"], "producer_tally")
        )
        row["first_failed_assertion"] = (
            "-"
            if raw["first_failed_assertion"] == "-"
            else parse_jcs_cell(raw["first_failed_assertion"], "first_failed_assertion")
        )
        require(row["producer_tally"] == tally, f"{expected_id}: producer tally derivation")
        require(row["first_failed_assertion"] == first, f"{expected_id}: failure derivation")
        expected_verdicts = []
        for checker, exit_item in zip(config["checkers"], checker_exits, strict=True):
            expected_verdicts.append(
                checker_observation(
                    capture_data[f"checker.{checker['name']}.stdout"],
                    checker["report_format"],
                    checker["name"],
                    exit_item["exit"],
                )
            )
        require(row["checker_verdicts"] == expected_verdicts, f"{expected_id}: checker derivation")
        expected_moved: list[dict[str, Any]] = []
        for artifact in config["artifacts"]:
            name = artifact["name"]
            if name not in emitted or artifact["observation_format"] == "none":
                continue
            require(artifact["observation_format"] == "ledger-dimension-v1", "artifact format")
            before_data = artifact_entry[name]
            require(before_data is not None, f"{expected_id}: dimension baseline absent")
            before = dimension_records(before_data)
            after = dimension_records(capture_data[f"artifact.{name}"])
            record = include["record"]
            require(record in before and record in after, f"{expected_id}: selected DIM record")
            expected_moved.append(
                {
                    "after": after[record],
                    "artifact": name,
                    "before": before[record],
                    "moved": before[record] != after[record],
                    "name": record,
                }
            )
        require(row["moved_values"] == expected_moved, f"{expected_id}: moved values")
        decoded.append(row)
        validate_outcome(row)
    return decoded


def audit_restore(
    project_root: pathlib.Path,
    evidence: pathlib.Path,
    config: dict[str, Any],
    entry_snapshot: EntrySnapshot,
    required_operation: str | None = None,
) -> list[dict[str, str]]:
    header, body = read_tsv_bytes((evidence / "restore.tsv").read_bytes(), "restore.tsv")
    require(header == RESTORE_COLUMNS, "restore.tsv: exact header")
    require(all(all(cell != "" for cell in cells) for cells in body), "restore.tsv: empty cell")
    rows = [dict(zip(header, cells, strict=True)) for cells in body]
    members = [
        ("target", item["name"], item["path"]) for item in config["targets"]
    ] + [("artifact", item["name"], item["path"]) for item in config["artifacts"]]
    require(bool(rows) and len(rows) % len(members) == 0, "restore event cardinality")
    event_count = len(rows) // len(members)
    for event_index in range(event_count):
        group = rows[event_index * len(members) : (event_index + 1) * len(members)]
        require(all(row["event"] == str(event_index + 1) for row in group), "restore event numbering")
        require(len({row["operation"] for row in group}) == 1, "restore operation grouping")
        for row, (kind, name, path_text) in zip(group, members, strict=True):
            require(row["kind"] == kind and row["name"] == name, "restore member order")
            require(row["operation"] in {"run", "signal", "error", "repair"}, "restore operation")
            require(row["restored"] in {"true", "false"}, "restore boolean")
            for state_name, digest_name in (
                ("entry_state", "entry_sha256"),
                ("pre_restore_state", "pre_restore_sha256"),
                ("restored_state", "restored_sha256"),
            ):
                require(row[state_name] in {"file", "absent"}, "restore state")
                if row[state_name] == "absent":
                    require(row[digest_name] == "-", "absent restore digest")
                else:
                    require(bool(SHA_RE.fullmatch(row[digest_name])), "file restore digest")
            entry = entry_snapshot[(kind, name)]
            entry_state = "file" if entry is not None else "absent"
            entry_digest = sha256_bytes(entry) if entry is not None else "-"
            require(
                row["entry_state"] == entry_state
                and row["entry_sha256"] == entry_digest,
                "entry state claim differs from fixture snapshot",
            )
            path = project_root / path_text
            restored_state = "file" if path.is_file() else "absent"
            restored_digest = sha256_bytes(path.read_bytes()) if restored_state == "file" else "-"
            require(
                row["restored_state"] == restored_state
                and row["restored_sha256"] == restored_digest,
                "restored state",
            )
            require(
                restored_state == entry_state and restored_digest == entry_digest,
                "restored bytes differ from fixture snapshot",
            )
            require(row["restored"] == "true", "restoration proof is false")
    if required_operation is not None:
        require(rows[-1]["operation"] == required_operation, f"final restore is not {required_operation}")
    return rows


def grade_a8_inventory(claimed_paths: Iterable[str], committed_paths: Iterable[str]) -> None:
    missing = set(claimed_paths) - set(committed_paths)
    require(not missing, f"committed claim has no file: {sorted(missing)}")


def grade_a8_clean_restore_claims(
    evidence: pathlib.Path,
    rows: list[dict[str, Any]],
    restore_rows: list[dict[str, str]],
    config: dict[str, Any],
) -> None:
    """Bind a clean fixture's final pre-restore claims to committed result/capture bytes."""
    require(bool(rows), "A8 clean run has no result row")
    member_count = len(config["targets"]) + len(config["artifacts"])
    final_event = restore_rows[-member_count:]
    require(all(row["operation"] == "run" for row in final_event), "A8 is not a clean run event")
    last = rows[-1]
    for claim in final_event:
        if claim["kind"] == "target":
            if claim["name"] == last["target"] and last["applied"]:
                expected_state = "file"
                expected_digest = last["target_mutant_sha256"]
            else:
                expected_state = claim["entry_state"]
                expected_digest = claim["entry_sha256"]
        else:
            role = f"artifact.{claim['name']}"
            reference = next(
                (
                    item
                    for item in last["captures"]
                    if last["applied"] and item["role"] == role
                ),
                None,
            )
            if reference is None:
                expected_state = "absent"
                expected_digest = "-"
            else:
                capture = evidence / reference["path"]
                expected_state = "file"
                expected_digest = sha256_bytes(capture.read_bytes())
        require(
            claim["pre_restore_state"] == expected_state
            and claim["pre_restore_sha256"] == expected_digest,
            f"A8 unsupported pre-restore claim for {claim['kind']} {claim['name']}",
        )


def audit_completed_evidence(
    project_root: pathlib.Path,
    evidence: pathlib.Path,
    contract_path: pathlib.Path,
    entry_snapshot: EntrySnapshot,
    *,
    invocation_config: bytes | None,
    invocation_include: bytes | None,
    required_restore_operation: str = "run",
) -> list[dict[str, Any]]:
    required = {
        "EVIDENCE.md",
        "config.json",
        "include_list.tsv",
        "results.tsv",
        "restore.tsv",
        "captures",
        "COMPLETE",
    }
    committed = {path.name for path in evidence.iterdir()}
    grade_a8_inventory(required, committed)
    require((evidence / "COMPLETE").read_bytes() == b"ablation-evidence-v1\n", "COMPLETE marker")
    config_data = (evidence / "config.json").read_bytes()
    include_data = (evidence / "include_list.tsv").read_bytes()
    if invocation_config is not None:
        require(config_data == invocation_config, "config.json is not the invocation input")
    if invocation_include is not None:
        require(include_data == invocation_include, "include_list.tsv is not the invocation input")
    config = parse_standalone_jcs(config_data, "config.json")
    rows = parse_and_audit_results(
        project_root,
        evidence,
        config,
        include_data,
        entry_snapshot,
    )
    audit_restore(
        project_root,
        evidence,
        config,
        entry_snapshot,
        required_operation=required_restore_operation,
    )
    command_entries = list(
        dict.fromkeys(command["argv"][0] for command in [config["producer"], *config["checkers"]])
    )
    evidence_lines = [
        "# Ablation evidence scope",
        "",
        "This directory records one ablation run's observations; it is not a self-contained replay package.",
        "The claims are the fields in `results.tsv`, the bytes and digests referenced under `captures/`, and the restoration events in `restore.tsv`.",
        "The accepted run configuration and include-list are included as `config.json` and `include_list.tsv`.",
        "Audit also requires the repository tree at the commit containing this directory: each target and any pre-existing artifact must match the final restored state.",
        "Reproducing a row additionally requires that repository tree, compatible external toolchain behavior, POSIX filesystem/process/signal semantics, and the child environment recorded in `config.json`; stable-input entry digests and external executable, library, OS, kernel, or hardware bytes are not captured here.",
        "Accordingly this evidence does not claim that the commit alone can re-run a row or that a new run will repeat an outcome.",
        "",
        "Command entry points required for reproduction:",
        *(f"- {canonical_json(entry)}" for entry in command_entries),
    ]
    evidence_data = ("\n".join(evidence_lines) + "\n").encode("utf-8")
    require(
        (evidence / "EVIDENCE.md").read_bytes() == evidence_data,
        "EVIDENCE.md fixed scope/prerequisites",
    )
    if command_entries == ["python3"]:
        require(
            sha256_bytes(evidence_data) == EVIDENCE_MD_SHA256,
            "EVIDENCE.md additional digest binding",
        )
    claimed = [
        reference["path"]
        for row in rows
        if row["applied"]
        for reference in row["captures"]
    ]
    committed_capture_paths = [
        path.relative_to(evidence).as_posix()
        for path in (evidence / "captures").rglob("*")
        if path.is_file()
    ]
    grade_a8_inventory(claimed, committed_capture_paths)
    return rows


def grade_a1(rows: list[dict[str, Any]]) -> None:
    require(len(rows) == 2, "A1 row cardinality")
    first, second = rows
    require(first["artifacts_emitted"] == ["probe"], "A1 emitting row has no fresh artifact")
    require(second["artifacts_emitted"] == [], "A1 non-emitting row observed residual artifact")
    require(
        int(second["producer_exit"]) == 0,
        "A1 non-emitting row observed stale reset state",
    )
    require(
        all(int(item["exit"]) == 0 for item in second["checker_exits"]),
        "A1 residual state changed the checker result",
    )


def grade_a2(rows: list[dict[str, Any]]) -> None:
    require(len(rows) == 2, "A2 row cardinality")
    counts = {int(row["match_count"]) for row in rows}
    require(0 in counts and any(count >= 2 for count in counts), "A2 does not exercise both floors")
    for row in rows:
        validate_outcome(row)
        require(not row["applied"], "A2 unusable row was executed")
    require(rows[0]["outcome"] != rows[1]["outcome"], "A2 floors collapsed")


def grade_a3(facts: dict[str, Any]) -> None:
    require(bool(facts["banked_prefix"]), "A3 interruption has no banked prefix")
    require(
        set(facts["banked_prefix"]).isdisjoint(facts["resume_observed"]),
        "A3 resume re-ran a banked row",
    )
    require(facts["resumed_results"] == facts["uninterrupted_results"], "A3 result bytes differ")
    require(facts["complete"], "A3 did not complete all unbanked rows")


def grade_a4(facts: dict[str, dict[str, Any]]) -> None:
    for path_name, observation in facts.items():
        require(observation["target_restored"], f"A4 {path_name}: target not restored")
        require(observation["artifact_restored"], f"A4 {path_name}: artifact not restored")


def grade_a5(include_data: bytes, evidence: pathlib.Path) -> None:
    _, rows = read_include(include_data)
    require(len(rows) >= 2, "A5 needs multiple mutants")
    for position, row in enumerate(rows, start=1):
        token = row["new_text"].encode("utf-8")
        capture_root = evidence / f"captures/{position:06d}"
        expected = [
            capture_root / "producer.stdout",
            capture_root / "producer.stderr",
            capture_root / "checker-01-probe_checker.stdout",
            capture_root / "checker-01-probe_checker.stderr",
            capture_root / "artifact-probe",
        ]
        require(
            {path.name for path in capture_root.iterdir()} == {path.name for path in expected},
            f"A5 capture set for row {position}",
        )
        for path in expected:
            require(path.is_file(), f"A5 missing capture for row {position}")
            require(token in path.read_bytes(), f"A5 capture does not correspond to row {position}")


def grade_a6(
    include_data: bytes,
    validation: dict[str, Any],
    named_selectors: list[tuple[str, str]],
    required_shape: dict[str, int] | None = None,
) -> None:
    header, rows = read_include(include_data)
    require(
        set(validation)
        == {"axis_counts", "columns", "list_sha256", "row_count", "rows", "schema"},
        "A6 validation member set",
    )
    require(validation.get("schema") == "ablation-list-validation-v1", "A6 validation schema")
    require(validation.get("columns") == header, "A6 column order")
    require(validation.get("list_sha256") == sha256_bytes(include_data), "A6 list digest")
    require(validation.get("row_count") == len(rows), "A6 row count")
    counts = dict(Counter(row["axis"] for row in rows))
    if required_shape is not None:
        require(len(rows) == sum(required_shape.values()), "A6 contract row count")
        require(counts == required_shape, "A6 contract axis counts")
    require(validation.get("axis_counts") == counts, "A6 axis counts")
    observed = validation.get("rows")
    require(isinstance(observed, list) and len(observed) == len(rows), "A6 accept-and-ignore parse")
    for source, item in zip(rows, observed, strict=True):
        resolved_target = item.get("target")
        require(
            item == {
                "cells": [source[name] for name in header],
                "row_id": f"{resolved_target}:{source['axis']}:{source['key']}",
                "target": resolved_target,
            },
            "A6 lossless row reproduction",
        )
    source_by_key = {(row["axis"], row["key"]): row for row in rows}
    output_by_key = {
        (item["cells"][header.index("axis")], item["cells"][header.index("key")]): item
        for item in observed
    }
    for selector in named_selectors:
        require(selector in source_by_key and selector in output_by_key, f"A6 named row {selector}")
        source = source_by_key[selector]
        item = output_by_key[selector]
        require(
            item["cells"][header.index("old_text")] == source["old_text"]
            and item["cells"][header.index("new_text")] == source["new_text"],
            f"A6 named mutation text {selector}",
        )


def new_projection(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    projection: list[dict[str, Any]] = []
    for row in rows:
        moved = [
            {
                "after": item["after"],
                "artifact": item["artifact"],
                "moved": item["moved"],
                "name": item["name"],
            }
            for item in row["moved_values"]
        ]
        verdict = row["checker_verdicts"][0]
        projection.append(
            {
                "axis": row["axis"],
                "key": row["key"],
                "record": row["record"],
                "old_text": row["old_text"],
                "new_text": row["new_text"],
                "producer_exit": row["producer_exit"],
                "producer_tally": row["producer_tally"],
                "first_failed_assertion": row["first_failed_assertion"],
                "artifacts_emitted": row["artifacts_emitted"],
                "moved_values": moved,
                "checker_exit": row["checker_exits"][0]["exit"],
                "checker_status": verdict["status"],
                "mismatch_names": verdict["mismatch_names"],
                "outcome": row["outcome"],
                "match_count": row["match_count"],
                "applied": row["applied"],
            }
        )
    return projection


def grade_a7(
    legacy_list_data: bytes,
    legacy_result_data: bytes,
    candidate_projection: list[dict[str, Any]],
) -> None:
    require(sha256_bytes(legacy_list_data) == A7_LIST_SHA256, "A7 legacy list authority")
    require(sha256_bytes(legacy_result_data) == A7_RESULTS_SHA256, "A7 legacy result authority")
    _, legacy_include_rows = read_include(legacy_list_data)
    legacy_header, legacy_result_cells = read_tsv_bytes(
        legacy_result_data,
        "legacy results.tsv",
    )
    require(legacy_header == LEGACY_RESULT_COLUMNS, "A7 legacy result header")
    legacy_results = [
        dict(zip(legacy_header, cells, strict=True))
        for cells in legacy_result_cells
    ]
    include_by_key = {
        (row["axis"], row["key"]): row
        for row in legacy_include_rows
    }
    legacy_by_key = {
        (row["axis"], row["key"]): row
        for row in legacy_results
    }
    candidate_by_key = {(row["axis"], row["key"]): row for row in candidate_projection}
    require(
        len(include_by_key) == len(legacy_include_rows),
        "A7 legacy include keys are not unique",
    )
    require(len(legacy_by_key) == len(legacy_results), "A7 legacy result keys are not unique")
    require(len(candidate_by_key) == len(candidate_projection), "A7 new keys are not unique")
    require(
        set(candidate_by_key) == set(legacy_by_key) == set(include_by_key),
        "A7 joined row set",
    )

    def legacy_integer(value: str, field: str) -> int:
        require(bool(re.fullmatch(r"0|[1-9][0-9]*", value)), f"A7 legacy {field} integer")
        return int(value)

    for key, legacy in legacy_by_key.items():
        include = include_by_key[key]
        candidate = candidate_by_key[key]
        require(candidate["axis"] == legacy["axis"], f"A7 {key}: axis")
        require(candidate["key"] == legacy["key"], f"A7 {key}: key")
        require(
            candidate["record"] == legacy["record"] == include["record"],
            f"A7 {key}: record",
        )
        require(
            candidate["old_text"] == include["old_text"]
            and candidate["new_text"] == include["new_text"],
            f"A7 {key}: mutation text",
        )
        stage_exit = legacy_integer(legacy["stage_exit"], "stage_exit")
        pass_count = legacy_integer(legacy["pass_count"], "pass_count")
        fail_count = legacy_integer(legacy["fail_count"], "fail_count")
        cmp_exit = legacy_integer(legacy["cmp_exit"], "cmp_exit")
        require(candidate["producer_exit"] == stage_exit, f"A7 {key}: producer exit")
        require(
            candidate["producer_tally"] == {"pass": pass_count, "fail": fail_count},
            f"A7 {key}: producer tally",
        )
        expected_failure: str = legacy["first_fail"] if legacy["first_fail"] else "-"
        require(
            candidate["first_failed_assertion"] == expected_failure,
            f"A7 {key}: first failure",
        )

        sidecar_written = legacy["sidecar_written"]
        require(sidecar_written in {"yes", "no"}, f"A7 {key}: sidecar spelling")
        require(
            ("stage023_sidecar" in candidate["artifacts_emitted"])
            == (sidecar_written == "yes"),
            f"A7 {key}: sidecar membership",
        )
        moved = [
            item
            for item in candidate["moved_values"]
            if item["artifact"] == "stage023_sidecar" and item["name"] == legacy["record"]
        ]
        moved_spelling = legacy["record_moved"]
        emitted_value = legacy["emitted_value"]
        if moved_spelling == "no_sidecar":
            require(not emitted_value and not moved, f"A7 {key}: no-sidecar mapping")
        else:
            require(moved_spelling in {"yes", "no"}, f"A7 {key}: moved spelling")
            require(len(moved) == 1, f"A7 {key}: moved-value cardinality")
            require(
                moved[0]["moved"] is (moved_spelling == "yes"),
                f"A7 {key}: moved boolean",
            )
            require(bool(emitted_value), f"A7 {key}: emitted value is blank")
            require(moved[0]["after"] == emitted_value, f"A7 {key}: emitted value")

        require(candidate["checker_exit"] == cmp_exit, f"A7 {key}: checker exit")
        require(bool(legacy["cmp_status"]), f"A7 {key}: checker status is blank")
        require(candidate["checker_status"] == legacy["cmp_status"], f"A7 {key}: checker status")
        require(bool(legacy["mismatch_names"]), f"A7 {key}: mismatch names are blank")
        expected_mismatches = (
            [] if legacy["mismatch_names"] == "none" else [legacy["mismatch_names"]]
        )
        require(
            candidate["mismatch_names"] == expected_mismatches,
            f"A7 {key}: mismatch names",
        )
        expected_outcome = (
            "PRODUCER_NONZERO"
            if stage_exit != 0
            else ("CHECKER_NONZERO" if cmp_exit != 0 else "ALL_ZERO")
        )
        require(candidate["outcome"] == expected_outcome, f"A7 {key}: outcome")
        require(
            candidate["match_count"] == 1 and candidate["applied"] is True,
            f"A7 {key}: application",
        )

    require(
        sha256_bytes(canonical_json(candidate_projection).encode("utf-8"))
        == A7_PROJECTION_SHA256,
        "A7 projection additional digest binding",
    )


def grade_a9(rows: list[dict[str, Any]]) -> None:
    require(bool(rows), "A9 has no driver rows")
    for row in rows:
        validate_outcome(row)
    require(
        {outcome_for(row) for row in rows}
        == {"PRODUCER_NONZERO", "CHECKER_NONZERO", "ALL_ZERO"},
        "A9 applied truth-table branches were not all observed",
    )
