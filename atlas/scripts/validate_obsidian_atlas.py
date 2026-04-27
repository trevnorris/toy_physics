#!/usr/bin/env python3
"""Validate the generated Obsidian-facing atlas."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

import yaml

import generate_obsidian_atlas as gen


ROOT = Path(__file__).resolve().parents[2]
ATLAS_DIR = ROOT / "atlas"
GRAPH_YAML = ROOT / "graph" / "fluid_universe_derivation_atlas_graph.yaml"
REPORT = ATLAS_DIR / "exports" / "validation_report.md"
WIKILINK_RE = re.compile(r"!?\[\[([^\]]+)\]\]")


def load_yaml(path: Path) -> Any:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def split_frontmatter(path: Path) -> tuple[dict[str, Any], str]:
    text = path.read_text(encoding="utf-8")
    if not text.startswith("---\n"):
        raise ValueError("missing YAML frontmatter")
    end = text.find("\n---", 4)
    if end < 0:
        raise ValueError("unterminated YAML frontmatter")
    frontmatter = yaml.safe_load(text[4:end]) or {}
    body = text[end + 4 :]
    return frontmatter, body


def note_path_for(atlas_dir: Path, node: dict[str, Any]) -> Path:
    return atlas_dir / "nodes" / gen.node_folder(node) / f"{node['id']}.md"


def link_target(raw: str) -> str:
    target = raw.split("|", 1)[0].split("#", 1)[0].strip()
    return Path(target).stem if target.endswith(".md") else target


def validate_notes(
    graph: dict[str, Any],
    atlas_dir: Path,
    errors: list[str],
    warnings: list[str],
) -> int:
    node_ids = {node["id"] for node in graph["nodes"]}
    count = 0
    for node in graph["nodes"]:
        path = note_path_for(atlas_dir, node)
        if not path.exists():
            errors.append(f"missing node note: {path.relative_to(atlas_dir)}")
            continue
        count += 1
        if path.stem != node["id"]:
            errors.append(f"{path}: filename stem does not match node id {node['id']}")
        try:
            frontmatter, body = split_frontmatter(path)
        except ValueError as exc:
            errors.append(f"{path}: {exc}")
            continue

        for key in ("id", "title", "type", "layer", "status", "summary_short", "source_graph_version"):
            if key not in frontmatter:
                errors.append(f"{path.relative_to(atlas_dir)} missing frontmatter key: {key}")
        if frontmatter.get("id") != node["id"]:
            errors.append(f"{path.relative_to(atlas_dir)} frontmatter id mismatch")
        if frontmatter.get("status") != node.get("status"):
            errors.append(f"{node['id']} status differs from source graph")
        if gen.GENERATED_HEADER not in body:
            errors.append(f"{path.relative_to(atlas_dir)} missing generated header")

        if node.get("layer") == "equation_anchor" or node.get("type") == "equation":
            if "## Equation" not in body:
                errors.append(f"{node['id']} equation node missing ## Equation")
        if node.get("layer") == "open_gate" or node["id"].startswith("OPEN_"):
            if "## What Remains Open" not in body or "## What Would Close It" not in body:
                errors.append(f"{node['id']} open gate missing open/closure sections")
        if node["id"].startswith("FIREWALL_") or node.get("type") in gen.STATUS_FIREWALL_TYPES:
            if "## Invalid Inference" not in body or "## Corrected Inference" not in body:
                errors.append(f"{node['id']} status firewall missing invalid/corrected inference sections")

        for raw in WIKILINK_RE.findall(body):
            target = link_target(raw)
            if target.endswith(".base"):
                continue
            if target and target not in node_ids:
                warnings.append(f"{node['id']} has unresolved wikilink: {raw}")
    return count


def validate_exports(graph: dict[str, Any], atlas_dir: Path, errors: list[str]) -> None:
    export_dir = atlas_dir / "exports"
    yaml_path = export_dir / "atlas_graph.yaml"
    json_path = export_dir / "atlas_graph.json"
    nodes_tsv = export_dir / "atlas_nodes.tsv"
    edges_tsv = export_dir / "atlas_edges.tsv"
    manifest_path = export_dir / gen.GENERATED_MANIFEST
    for path in (yaml_path, json_path, nodes_tsv, edges_tsv, manifest_path):
        if not path.exists():
            errors.append(f"missing export: {path.relative_to(atlas_dir)}")
    if yaml_path.exists():
        exported = load_yaml(yaml_path)
        if len(exported.get("nodes", [])) != len(graph["nodes"]):
            errors.append("atlas_graph.yaml node count mismatch")
        if len(exported.get("edges", [])) != len(graph["edges"]):
            errors.append("atlas_graph.yaml edge count mismatch")
    if json_path.exists():
        exported = json.loads(json_path.read_text(encoding="utf-8"))
        if len(exported.get("nodes", [])) != len(graph["nodes"]):
            errors.append("atlas_graph.json node count mismatch")
        if len(exported.get("edges", [])) != len(graph["edges"]):
            errors.append("atlas_graph.json edge count mismatch")
    if manifest_path.exists():
        for raw_line in manifest_path.read_text(encoding="utf-8").splitlines():
            rel_path = raw_line.strip()
            if not rel_path or rel_path.startswith("#"):
                continue
            if not (atlas_dir / rel_path).exists():
                errors.append(f"generated manifest points to missing file: {rel_path}")


def validate_bases(atlas_dir: Path, errors: list[str]) -> int:
    count = 0
    for path in sorted((atlas_dir / "bases").glob("*.base")):
        count += 1
        try:
            data = load_yaml(path)
        except yaml.YAMLError as exc:
            errors.append(f"{path.relative_to(atlas_dir)} invalid YAML: {exc}")
            continue
        if "views" not in data:
            errors.append(f"{path.relative_to(atlas_dir)} missing views")
    return count


def validate_canvas(atlas_dir: Path, errors: list[str]) -> int:
    count = 0
    for path in sorted((atlas_dir / "canvas").glob("*.canvas")):
        count += 1
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            errors.append(f"{path.relative_to(atlas_dir)} invalid JSON: {exc}")
            continue
        positions: dict[tuple[Any, Any], str] = {}
        for node in data.get("nodes", []):
            if node.get("type") != "file":
                continue
            coords = (node.get("x"), node.get("y"))
            if coords in positions:
                errors.append(
                    f"{path.relative_to(atlas_dir)} overlapping file cards at {coords}: "
                    f"{positions[coords]} and {node.get('id')}"
                )
            else:
                positions[coords] = str(node.get("id"))
            file_path = atlas_dir / node.get("file", "")
            if not file_path.exists():
                errors.append(f"{path.relative_to(atlas_dir)} file card points to missing file: {node.get('file')}")
    return count


def write_report(
    graph: dict[str, Any],
    atlas_dir: Path,
    note_count: int,
    base_count: int,
    canvas_count: int,
    errors: list[str],
    warnings: list[str],
) -> None:
    report = atlas_dir / "exports" / "validation_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Obsidian Atlas Validation Report",
        "",
        f"- source nodes: `{len(graph['nodes'])}`",
        f"- source edges: `{len(graph['edges'])}`",
        f"- generated node notes: `{note_count}`",
        f"- bases: `{base_count}`",
        f"- canvas files: `{canvas_count}`",
        f"- errors: `{len(errors)}`",
        f"- warnings: `{len(warnings)}`",
        "",
        "## Errors",
        "",
    ]
    lines.extend(f"- {error}" for error in errors) if errors else lines.append("- none")
    lines.extend(["", "## Warnings", ""])
    lines.extend(f"- {warning}" for warning in warnings) if warnings else lines.append("- none")
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=GRAPH_YAML, help="Input graph YAML")
    parser.add_argument("--atlas", type=Path, default=ATLAS_DIR, help="Atlas vault root")
    args = parser.parse_args()

    graph = load_yaml(args.graph)
    errors: list[str] = []
    warnings: list[str] = []

    note_count = validate_notes(graph, args.atlas, errors, warnings)
    validate_exports(graph, args.atlas, errors)
    base_count = validate_bases(args.atlas, errors)
    canvas_count = validate_canvas(args.atlas, errors)
    write_report(graph, args.atlas, note_count, base_count, canvas_count, errors, warnings)

    if errors:
        print(f"FAIL: {len(errors)} error(s); see {REPORT}")
        raise SystemExit(1)
    print(f"OK: Obsidian atlas validation passed ({note_count} notes, {base_count} bases, {canvas_count} canvas files)")
    if warnings:
        print(f"warnings: {len(warnings)}")


if __name__ == "__main__":
    main()
