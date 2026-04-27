#!/usr/bin/env python3
"""Validate the Fluid Universe derivation atlas graph artifacts."""

from __future__ import annotations

import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import yaml


ROOT = Path(__file__).resolve().parents[1]
GRAPH_DIR = Path(__file__).resolve().parent
GRAPH_JSON = GRAPH_DIR / "fluid_universe_derivation_atlas_graph.json"
GRAPH_YAML = GRAPH_DIR / "fluid_universe_derivation_atlas_graph.yaml"
MANIFEST_JSON = GRAPH_DIR / "fluid_universe_derivation_atlas_paper_insertion_manifest.json"


def fail_with(message: str) -> None:
    print(f"FAIL: {message}", file=sys.stderr)
    raise SystemExit(1)


def load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        fail_with(f"{path} is not valid JSON: {exc}")


def load_yaml(path: Path) -> Any:
    try:
        return yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        fail_with(f"{path} is not valid YAML: {exc}")


def source_path(path: str) -> str:
    return re.sub(r":V2.*$", "", str(path))


def repo_path_exists(path: str) -> bool:
    return (ROOT / source_path(path)).exists()


def main() -> None:
    graph = load_yaml(GRAPH_YAML)
    generated_json = load_json(GRAPH_JSON)
    manifest = load_json(MANIFEST_JSON)

    errors: list[str] = []
    warnings: list[str] = []

    if graph != generated_json:
        errors.append("generated graph JSON differs from YAML; run python3 graph/sync_graph_formats.py")
    if graph["node_count"] != len(graph["nodes"]):
        errors.append("node_count does not match nodes length")
    if graph["edge_count"] != len(graph["edges"]):
        errors.append("edge_count does not match edges length")

    node_ids = [node["id"] for node in graph["nodes"]]
    duplicate_nodes = [node_id for node_id, count in Counter(node_ids).items() if count > 1]
    if duplicate_nodes:
        errors.append(f"duplicate node ids: {', '.join(duplicate_nodes)}")

    node_id_set = set(node_ids)
    for index, edge in enumerate(graph["edges"]):
        source = edge.get("source")
        target = edge.get("target")
        relation = edge.get("relation")
        if not source:
            errors.append(f"edge {index} missing source")
        if not target:
            errors.append(f"edge {index} missing target")
        if not relation:
            errors.append(f"edge {index} missing relation")
        if source and source not in node_id_set:
            errors.append(f"edge {index} source not found: {source}")
        if target and target not in node_id_set:
            errors.append(f"edge {index} target not found: {target}")

    source_resolution = graph.get("source_resolution", [])
    if not source_resolution:
        errors.append("source_resolution is empty")
    for entry in source_resolution:
        canonical = entry.get("canonical_ref")
        kind = entry.get("source_kind")
        if not entry.get("legacy_ref"):
            errors.append(f"source_resolution entry missing legacy_ref: {entry}")
        if not canonical:
            errors.append(f"source_resolution entry missing canonical_ref: {entry}")
        if not kind:
            errors.append(f"source_resolution entry missing source_kind: {entry}")
        if canonical and kind != "synthetic_provenance" and not repo_path_exists(canonical):
            warnings.append(f"resolved source path missing: {canonical}")

    for node in graph["nodes"]:
        for path in node.get("sources", []):
            if path == "derived_from_current_task":
                continue
            if not repo_path_exists(path):
                warnings.append(f"node {node['id']} source path missing: {path}")

        for path in node.get("files", []):
            if not repo_path_exists(path):
                warnings.append(f"node {node['id']} file path missing: {path}")

        if node.get("future_paper_needed") and node.get("source_kind") != "future_paper_note":
            errors.append(
                f"node {node['id']} has future_paper_needed without "
                "source_kind=future_paper_note"
            )

        tex_anchor = node.get("tex_anchor")
        if tex_anchor:
            tex_file = tex_anchor.get("file")
            tex_line = tex_anchor.get("line")
            if not tex_file:
                errors.append(f"node {node['id']} tex_anchor missing file")
            elif not repo_path_exists(tex_file):
                warnings.append(f"node {node['id']} tex_anchor file missing: {tex_file}")
            if not isinstance(tex_line, int) or tex_line < 1:
                errors.append(f"node {node['id']} tex_anchor has invalid line: {tex_line}")
            if not tex_anchor.get("match_basis"):
                errors.append(f"node {node['id']} tex_anchor missing match_basis")
            if not tex_anchor.get("confidence"):
                errors.append(f"node {node['id']} tex_anchor missing confidence")

        if (
            node.get("layer") == "source_section_anchor"
            and node.get("source_kind") == "paper"
            and not tex_anchor
        ):
            warnings.append(f"paper section anchor missing tex_anchor: {node['id']}")

    manifest_entries = manifest.get("entries", [])
    if len(manifest_entries) != 16:
        errors.append(f"manifest should have 16 entries, found {len(manifest_entries)}")

    for entry in manifest_entries:
        block_id = entry.get("block_id", "(missing block_id)")
        target = entry.get("resolved_target_path")
        if not target:
            errors.append(f"{block_id} missing resolved_target_path")
        elif not repo_path_exists(target):
            warnings.append(f"{block_id} target path missing: {target}")

        if entry.get("future_paper_needed") and entry.get("resolved_target_kind") != "future_paper_note":
            errors.append(
                f"{block_id} has future_paper_needed without "
                "resolved_target_kind=future_paper_note"
            )

        for candidate in entry.get("candidate_full_draft_paths_or_globs", []):
            if not repo_path_exists(candidate):
                warnings.append(f"{block_id} candidate path missing: {candidate}")

    for message in warnings:
        print(f"WARN: {message}", file=sys.stderr)

    if errors:
        for message in errors:
            print(f"ERROR: {message}", file=sys.stderr)
        fail_with(f"{len(errors)} validation error(s)")

    print("OK: graph validation passed")
    print(f"nodes: {len(graph['nodes'])}")
    print(f"edges: {len(graph['edges'])}")
    print(f"manifest entries: {len(manifest_entries)}")
    print(f"warnings: {len(warnings)}")


if __name__ == "__main__":
    main()
