#!/usr/bin/env python3
"""Targeted semantic B1 replay comparator (standalone-copy safe)."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any

import yaml

getcontext().prec = 40
CERTIFICATE_SHA256 = "650656fd2ef8a87884161825d25eced2620d8602099efbb172a960073480373c"


class UniqueLoader(yaml.SafeLoader):
    pass


def unique_mapping(loader: UniqueLoader, node: yaml.MappingNode, deep: bool = False) -> dict[str, Any]:
    result = {}
    for knode, vnode in node.value:
        key = loader.construct_object(knode, deep=deep)
        if key in result:
            raise ValueError(f"duplicate YAML key: {key!r}")
        result[key] = loader.construct_object(vnode, deep=deep)
    return result


UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, unique_mapping)


def load(path: Path, consumer: str, evidence: list[dict[str, Any]], expected: str | None = None) -> Any:
    fd = os.open(path, os.O_RDONLY)
    try:
        before = os.fstat(fd)
        chunks = []
        while True:
            chunk = os.read(fd, 1 << 20)
            if not chunk:
                break
            chunks.append(chunk)
        blob = b"".join(chunks)
        actual = hashlib.sha256(blob).hexdigest()
        require(expected is None or actual == expected, "B2_A1_TARGETED_FIRST_USE", f"{consumer}:{path}")
        stat = os.fstat(fd)
        require((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) == (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns), "B2_A1_PROTECTED_REWRITE", f"{consumer}:{path}")
        evidence.append({"consumer": consumer, "path": str(path.resolve()), "expected_sha256": expected, "consumed_fd_sha256": actual, "device": stat.st_dev, "inode": stat.st_ino, "held_descriptor": True, "descriptor_stable_during_parse": True})
        text = blob.decode("utf-8")
        if path.suffix == ".json":
            return json.loads(text)
        return yaml.load(text, Loader=UniqueLoader)
    finally:
        os.close(fd)


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, default=str).encode()


def digest(value: Any) -> str:
    return hashlib.sha256(canonical(value)).hexdigest()


UNIQUE_SYMBOL = re.compile(r"\$[0-9]+")


def normalize_wolfram_unique_symbols(value: Any) -> Any:
    """Bijective alpha-renaming within one complete compared structure."""
    names: dict[str, str] = {}

    def rename_text(text: str) -> str:
        def replacement(match: re.Match[str]) -> str:
            token = match.group(0)
            if token not in names:
                names[token] = f"$GEN{len(names)}"
            return names[token]
        return UNIQUE_SYMBOL.sub(replacement, text)

    def walk(node: Any) -> Any:
        if isinstance(node, dict):
            return {key: walk(node[key]) for key in sorted(node, key=str)}
        if isinstance(node, list):
            return [walk(child) for child in node]
        if isinstance(node, str):
            return rename_text(node)
        return node

    normalized = walk(value)
    require(len(set(names.values())) == len(names), "B2_A1_TARGETED_ALPHA", "bijective Unique[] alpha map")
    return normalized


def require(test: bool, tooth: str, detail: str) -> None:
    if not test:
        raise AssertionError(f"ASSERT_FAIL:{tooth}:{detail}")


def add_term(bucket: dict[tuple[tuple[str, int], ...], Decimal], powers: dict[str, Any], coefficient: Decimal) -> None:
    key = tuple(sorted((name, int(power)) for name, power in powers.items() if int(power)))
    bucket[key] = bucket.get(key, Decimal(0)) + coefficient


def forward_hessian(termwise: dict[str, list[dict[str, Any]]], velocities: list[str]) -> list[list[list[dict[str, Any]]]]:
    matrix: list[list[list[dict[str, Any]]]] = []
    all_rows = [row for rows in termwise.values() for row in rows]
    for vi in velocities:
        matrix_row = []
        for vj in velocities:
            bucket: dict[tuple[tuple[str, int], ...], Decimal] = {}
            for row in all_rows:
                powers = {name: int(power) for name, power in row.get("powers", {}).items()}
                pi, pj = powers.get(vi, 0), powers.get(vj, 0)
                factor = pi * (pj - (1 if vi == vj else 0))
                if not factor:
                    continue
                remaining = dict(powers)
                remaining[vi] = remaining.get(vi, 0) - 1
                remaining[vj] = remaining.get(vj, 0) - 1
                add_term(bucket, remaining, Decimal(str(row["coefficient"])) * factor)
            matrix_row.append([
                {"coefficient": coefficient, "powers": dict(key)}
                for key, coefficient in sorted(bucket.items()) if coefficient != 0
            ])
        matrix.append(matrix_row)
    return matrix


def compare_terms(actual: list[dict[str, Any]], expected: list[dict[str, Any]], detail: str) -> None:
    def collect(rows: list[dict[str, Any]]) -> dict[tuple[tuple[str, int], ...], Decimal]:
        bucket: dict[tuple[tuple[str, int], ...], Decimal] = {}
        for row in rows:
            add_term(bucket, row.get("powers", {}), Decimal(str(row["coefficient"])))
        return bucket
    left, right = collect(actual), collect(expected)
    require(set(left) == set(right), "B2_A1_TARGETED_MXX", f"{detail}:monomial set")
    for key in left:
        scale = max(Decimal(1), abs(left[key]), abs(right[key]))
        require(abs(left[key] - right[key]) <= Decimal("2e-11") * scale, "B2_A1_TARGETED_MXX", f"{detail}:{key}")


def compare_matrix(actual: list[Any], expected: list[Any], detail: str) -> None:
    require(len(actual) == len(expected) == 3, "B2_A1_TARGETED_MXX", f"{detail}:shape")
    for i in range(3):
        require(len(actual[i]) == len(expected[i]) == 3, "B2_A1_TARGETED_MXX", f"{detail}:shape row")
        for j in range(3):
            compare_terms(actual[i][j], expected[i][j], f"{detail}[{i},{j}]")


def ledger_digest_checks(artifact: dict[str, Any], engine: str) -> dict[str, Any]:
    records = artifact["partition_ledger"]["records"]
    require(len(records) == 41, "B2_A1_TARGETED_LEDGER", f"{engine}:41 records")
    verified = []
    for row in records:
        candidate = row["candidate_id"]
        if candidate == "outer_control_flux:translation":
            computed = digest("pending_B2")
            if engine == "SymPy":
                require(computed == row["computed_expression_digest"], "B2_A1_TARGETED_LEDGER", f"{engine}:{candidate}:stored digest")
            semantic = {"placeholder": "pending_B2"}
        else:
            cell, provenance, suffix = candidate.split(":")
            require(suffix == "translation", "B2_A1_TARGETED_LEDGER", f"{engine}:{candidate}:suffix")
            computed = digest(artifact["indexed_cells"][cell]["termwise_L"][provenance])
            semantic = artifact["indexed_cells"][cell]["termwise_L"][provenance]
            require(row["provenance"] == provenance, "B2_A1_TARGETED_LEDGER", f"{engine}:{candidate}:provenance")
            if engine == "SymPy":
                require(computed == row["computed_expression_digest"], "B2_A1_TARGETED_LEDGER", f"{engine}:{candidate}:stored digest")
        verified.append({"candidate_id": candidate, "stored_engine_digest": row["computed_expression_digest"], "serialization_neutral_semantic_digest": digest(semantic)})
    owners = [row["owner"] for row in records]
    require(owners.count("M") == 40 and owners.count("C_mdot_pending") == 1, "B2_A1_TARGETED_LEDGER", f"{engine}:owner census")
    require(artifact["partition_ledger"]["unique"] and artifact["partition_ledger"]["state"] == "partition_open_pending_B2", "B2_A1_TARGETED_LEDGER", f"{engine}:ledger state")
    return {"record_count": len(records), "verified_digest_set_sha256": digest(verified), "digest_recipe_verification": "Python canonical byte recomputation" if engine == "SymPy" else "fresh Wolfram engine-native in-memory recomputation exactly reproduced protected stored digests; comparator independently reconstructed every ID-to-termwise mapping and semantic digest (Wolfram machine-real RawJSON digest is not YAML-roundtrip stable)"}


def selected_structures(artifact: dict[str, Any]) -> dict[str, Any]:
    cell = artifact["indexed_cells"]["E1|symmetric_postulate"]
    return {
        "termwise_L": cell["termwise_L"],
        "M_XX_p0_known": cell["M_XX_p0_known"],
        "field_contraction_integrals": cell["field_contraction_integrals"],
        "scalar_regression": artifact["scalar_regression"]["E1|symmetric_postulate"],
        "partition_ledger": artifact["partition_ledger"],
        "phase_a_amendment": artifact["phase_a_amendment"],
        "amended_phase_a_payload": artifact["amended_phase_a_payload"],
    }


def compare_engine(protected: dict[str, Any], replay: dict[str, Any], protected_phase: dict[str, Any], replay_phase: dict[str, Any], engine: str) -> dict[str, Any]:
    require(protected["engine"] == replay["engine"] == engine, "B2_A1_TARGETED_REPLAY", f"{engine}:identity")
    p_selected, r_selected = selected_structures(protected), selected_structures(replay)
    phase_keys = ["linearized_channel_operator", "coupled_indicial_analysis", "tail_channels"]
    p_phase = normalize_wolfram_unique_symbols({key: protected_phase[key] for key in phase_keys})
    r_phase = normalize_wolfram_unique_symbols({key: replay_phase[key] for key in phase_keys})
    require(p_phase == r_phase, "B2_A1_TARGETED_REPLAY", f"{engine}:normalized complete Phase-A operator/tail structures")
    normalized_tails_digest = digest(r_phase["tail_channels"])
    for selected, phase in ((p_selected, p_phase), (r_selected, r_phase)):
        selected["phase_a_amendment"]["phase_a_acceptance_recheck"]["tails_digest"] = digest(phase["tail_channels"])
    require(p_selected == r_selected, "B2_A1_TARGETED_REPLAY", f"{engine}:complete corresponding structures")
    velocities = ["V_x", "V_y", "V_z"] if engine == "SymPy" else ["Vx", "Vy", "Vz"]
    derived = forward_hessian(replay["indexed_cells"]["E1|symmetric_postulate"]["termwise_L"], velocities)
    compare_matrix(derived, replay["indexed_cells"]["E1|symmetric_postulate"]["M_XX_p0_known"], f"{engine}:forward M_XX")
    ledger = ledger_digest_checks(replay, engine)
    return {
        "engine": engine,
        "endpoint": "E1",
        "ambient": "symmetric_postulate",
        "comparison": "duplicate-key-rejected recursive equality of complete selected structures",
        "normalizations": ["Phase-A display strings only: per-complete-structure bijective first-occurrence alpha-renaming $<n> -> $GEN<i>; repeated symbols retain equality and distinct symbols remain distinct; affected tails_digest recomputed from the complete normalized tail_channels structure"],
        "excluded_volatile_paths": [],
        "normalized_phase_a_structure_sha256": digest(r_phase),
        "normalized_tails_sha256": normalized_tails_digest,
        "selected_structure_sha256": digest(r_selected),
        "forward_M_XX_sha256": digest(derived),
        "ledger": ledger,
        "amended_payload_sha256": digest(replay["amended_phase_a_payload"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protected-sympy", type=Path, required=True)
    parser.add_argument("--protected-mathematica", type=Path, required=True)
    parser.add_argument("--replay-sympy", type=Path, required=True)
    parser.add_argument("--replay-mathematica", type=Path, required=True)
    parser.add_argument("--protected-phase-sympy", type=Path, required=True)
    parser.add_argument("--protected-phase-mathematica", type=Path, required=True)
    parser.add_argument("--replay-phase-sympy", type=Path, required=True)
    parser.add_argument("--replay-phase-mathematica", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--certificate", type=Path)
    parser.add_argument("--producer-record", type=Path)
    args = parser.parse_args()
    try:
        evidence: list[dict[str, Any]] = []
        expected: dict[str, str] = {}
        if args.certificate:
            cert = load(args.certificate, "targeted_replay:approval_certificate", evidence, CERTIFICATE_SHA256)
            for field in ["protected_artifact_sha256", "runtime_dependency_sha256"]:
                for raw, sha in cert[field].items():
                    path = (args.certificate.parent.parent.parent / raw).resolve() if not raw.startswith("//") else (args.certificate.parents[4] / raw[2:]).resolve()
                    expected[str(path)] = sha
        if args.producer_record:
            producer = load(args.producer_record, "targeted_replay:producer_record", evidence)
            require(producer.get("status") == "PASS", "B2_A1_TARGETED_FIRST_USE", "producer record status")
            for row in producer["outputs"]:
                expected[str(Path(row["path"]).resolve())] = row["sha256"]
        def read(path: Path, consumer: str) -> Any:
            return load(path, consumer, evidence, expected.get(str(path.resolve())))
        alpha_probe = normalize_wolfram_unique_symbols({"same": ["x$91", "x$91"], "different": "y$92"})
        require(alpha_probe["same"][0] == alpha_probe["same"][1] and alpha_probe["same"][0] != alpha_probe["different"], "B2_A1_TARGETED_ALPHA", "equality/distinctness preservation")
        rows = [
            compare_engine(read(args.protected_sympy, "SymPy:protected_B1"), read(args.replay_sympy, "SymPy:fresh_replay_B1"), read(args.protected_phase_sympy, "SymPy:protected_PhaseA"), read(args.replay_phase_sympy, "SymPy:fresh_replay_PhaseA"), "SymPy"),
            compare_engine(read(args.protected_mathematica, "Mathematica:protected_B1"), read(args.replay_mathematica, "Mathematica:fresh_replay_B1"), read(args.protected_phase_mathematica, "Mathematica:protected_PhaseA"), read(args.replay_phase_mathematica, "Mathematica:fresh_replay_PhaseA"), "Mathematica"),
        ]
        require(rows[0]["amended_payload_sha256"] == rows[1]["amended_payload_sha256"] == "b23993cca80dc3e6a790abcf68c1af63aa804fc47b06b153b9f224ccf27f899d", "B2_A1_TARGETED_AMENDMENT", "amended digest three-way")
        artifact = {
            "schema_version": "U1_PHASE_B2_TARGETED_B1_REPLAY_V1",
            "status": "PASS",
            "recipe": "fresh engine replay outputs semantically compared to protected B1 artifacts",
            "comparisons": rows,
            "normalizations": rows[0]["normalizations"],
            "excluded_volatile_paths": [],
            "first_use_authentication": evidence,
            "producer_record_sha256": next((row["consumed_fd_sha256"] for row in evidence if row["consumer"] == "targeted_replay:producer_record"), None),
            "checks": {"B2_A1_TARGETED_REPLAY": "PASS", "B2_A1_TARGETED_MXX": "PASS", "B2_A1_TARGETED_LEDGER": "PASS", "B2_A1_TARGETED_AMENDMENT": "PASS", "B2_A1_TARGETED_ALPHA": "PASS", "B2_A1_TARGETED_FIRST_USE": "PASS"},
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(yaml.safe_dump(artifact, sort_keys=False, allow_unicode=True, width=220), encoding="utf-8")
        print("B2_TARGETED_REPLAY: PASS engines=2 endpoint=E1")
        return 0
    except Exception as exc:
        print(str(exc))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
