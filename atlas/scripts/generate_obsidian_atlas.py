#!/usr/bin/env python3
"""Generate the Obsidian-facing Fluid Universe atlas from graph YAML."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml


ROOT = Path(__file__).resolve().parents[2]
ATLAS_DIR = ROOT / "atlas"
GRAPH_YAML = ROOT / "graph" / "fluid_universe_derivation_atlas_graph.yaml"

NODE_DIR = ATLAS_DIR / "nodes"
VIEW_DIR = ATLAS_DIR / "views"
BASE_DIR = ATLAS_DIR / "bases"
CANVAS_DIR = ATLAS_DIR / "canvas"
EXPORT_DIR = ATLAS_DIR / "exports"
TEMPLATE_DIR = ATLAS_DIR / "templates"

GENERATED_HEADER = "<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->"
GENERATED_MANIFEST = "generated_manifest.txt"
ATLAS_ID_RE = re.compile(r"^[A-Z][A-Z0-9]*(?:_[A-Z0-9]+)+$")
BASE_FILENAMES = (
    "claims.base",
    "open_gates.base",
    "physical_objects.base",
    "status_firewalls.base",
    "source_crosswalk.base",
    "equations.base",
)
EXPORT_FILENAMES = (
    "atlas_graph.yaml",
    "atlas_graph.json",
    "atlas_nodes.tsv",
    "atlas_edges.tsv",
)

NODE_FOLDER_BY_LAYER = {
    "physical_ontology": "physical",
    "physical_register": "physical",
    "math_object": "math",
    "equation_anchor": "equations",
    "claim_theorem": "claims",
    "open_gate": "open_gates",
    "file_anchor": "sources",
    "source_section_anchor": "sources",
    "derivation": "derivations",
    "query_validation": "audits",
    "status_audit": "audits",
    "extension_workflow_step": "extensions",
    "atlas_meta": "meta",
    "atlas_completion_phase": "meta",
    "paper_backlink": "meta",
}

STATUS_FIREWALL_TYPES = {
    "status_firewall_rule",
    "ontology_firewall",
}

CANVAS_SPECS = {
    "00_backbone.canvas": {
        "title": "Backbone",
        "seeds": [
            "CLAIM_PARENT_ACTION_CURRENT_EXACT",
            "CLAIM_PROJECTION_NOT_REDUCTION",
            "CLAIM_1PN_EIH_WITHIN_HIERARCHY",
            "CLAIM_2PN_ADM_WITHIN_HIERARCHY",
            "CLAIM_25PN_QUAD_NARROWING",
            "CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY",
            "CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL",
            "OPEN_QUAD_NORMALIZATION",
        ],
        "terms": ["parent action", "projection", "PN"],
    },
    "01_physical_ontology.canvas": {
        "title": "Physical Ontology",
        "terms": ["physical throat brane bulk outgoing support mouth charge"],
        "layers": {"physical_ontology", "physical_register"},
    },
    "02_moving_throat_pde.canvas": {
        "title": "Moving Throat PDE",
        "terms": ["moving throat PDE wall branch BdG Maxwell mixed grouped"],
    },
    "03_quadrupole_normalization.canvas": {
        "title": "Quadrupole Normalization",
        "terms": ["quadrupole normalization P0 outgoing 2.5PN 4PN tail"],
    },
    "04_pn_chain.canvas": {
        "title": "PN Chain",
        "terms": ["1PN 2PN 2.5PN 3PN 4PN 5PN"],
    },
    "05_open_gates.canvas": {
        "title": "Open Gates",
        "layers": {"open_gate"},
        "terms": ["open gate"],
    },
    "06_charge_ontology_firewall.canvas": {
        "title": "Charge Ontology Firewall",
        "terms": ["charge q_eff q_star eta_Q circulation kappa_rho Maxwell"],
    },
    "07_status_firewalls.canvas": {
        "title": "Status Firewalls",
        "terms": ["firewall invalid inference corrected inference"],
        "types": STATUS_FIREWALL_TYPES | {"invalid_inference"},
    },
    "08_lepton_atom_anomaly.canvas": {
        "title": "Lepton Atom Anomaly",
        "terms": ["lepton atom anomaly g-2 hydrogen"],
    },
}


def load_graph(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def ensure_dirs() -> None:
    for path in [
        NODE_DIR / "physical",
        NODE_DIR / "math",
        NODE_DIR / "equations",
        NODE_DIR / "claims",
        NODE_DIR / "open_gates",
        NODE_DIR / "status_firewalls",
        NODE_DIR / "sources",
        NODE_DIR / "derivations",
        NODE_DIR / "audits",
        NODE_DIR / "extensions",
        NODE_DIR / "meta",
        VIEW_DIR,
        BASE_DIR,
        CANVAS_DIR,
        EXPORT_DIR,
        TEMPLATE_DIR,
    ]:
        path.mkdir(parents=True, exist_ok=True)


def manifest_path() -> Path:
    return EXPORT_DIR / GENERATED_MANIFEST


def clean_previous_generated_files() -> None:
    path = manifest_path()
    if not path.exists():
        return
    atlas_root = ATLAS_DIR.resolve()
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        rel_path = raw_line.strip()
        if not rel_path or rel_path.startswith("#"):
            continue
        target = (ATLAS_DIR / rel_path).resolve()
        try:
            target.relative_to(atlas_root)
        except ValueError:
            continue
        if target.is_file():
            target.unlink()


def generated_markdown_files(directory: Path) -> list[Path]:
    paths: list[Path] = []
    if not directory.exists():
        return paths
    for path in sorted(directory.rglob("*.md")):
        if GENERATED_HEADER in path.read_text(encoding="utf-8"):
            paths.append(path)
    return paths


def write_generated_manifest() -> None:
    paths: set[Path] = set()
    for directory in (NODE_DIR, VIEW_DIR, TEMPLATE_DIR):
        paths.update(generated_markdown_files(directory))
    for filename in BASE_FILENAMES:
        paths.add(BASE_DIR / filename)
    for filename in CANVAS_SPECS:
        paths.add(CANVAS_DIR / filename)
    for filename in EXPORT_FILENAMES:
        paths.add(EXPORT_DIR / filename)

    entries = sorted(path.relative_to(ATLAS_DIR).as_posix() for path in paths if path.exists())
    content = "# Generated by atlas/scripts/generate_obsidian_atlas.py\n" + "\n".join(entries) + "\n"
    manifest_path().write_text(content, encoding="utf-8")


def compact(value: Any, limit: int = 140) -> str:
    text = " ".join(str(value or "").split())
    return text if len(text) <= limit else f"{text[:limit - 3]}..."


def title_for(node: dict[str, Any]) -> str:
    return str(node.get("label") or node["id"]).strip()


def sanitize_tag(value: Any) -> str:
    text = str(value or "unknown").lower()
    text = re.sub(r"[^a-z0-9_/-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_/")
    return text or "unknown"


def wiki(node_id: str) -> str:
    return f"[[{node_id}]]"


def looks_like_atlas_id(value: str) -> bool:
    return "/" not in value and "." not in value and bool(ATLAS_ID_RE.fullmatch(value))


def node_folder(node: dict[str, Any]) -> str:
    if node["id"].startswith("FIREWALL_") or node.get("type") in STATUS_FIREWALL_TYPES:
        return "status_firewalls"
    return NODE_FOLDER_BY_LAYER.get(node.get("layer"), "meta")


def node_path(node: dict[str, Any]) -> Path:
    return NODE_DIR / node_folder(node) / f"{node['id']}.md"


def vault_path(path: Path) -> str:
    return str(path.relative_to(ATLAS_DIR))


def graph_text(node: dict[str, Any]) -> str:
    return " ".join(flatten(node)).lower()


def flatten(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, dict):
        out = []
        for key, child in value.items():
            out.append(str(key))
            out.extend(flatten(child))
        return out
    if isinstance(value, list):
        out = []
        for child in value:
            out.extend(flatten(child))
        return out
    return [str(value)]


def source_paths(node: dict[str, Any]) -> list[str]:
    paths: list[str] = []
    for key in ("file", "source_note_file"):
        value = node.get(key)
        if isinstance(value, str) and not looks_like_atlas_id(value):
            paths.append(value)
    for key in (
        "files",
        "sources",
        "source_files",
        "source_note_files",
        "canonical_target_files",
        "legacy_sources",
        "legacy_files",
    ):
        for value in node.get(key, []) or []:
            if isinstance(value, str) and not looks_like_atlas_id(value):
                paths.append(value)
    anchor = node.get("tex_anchor") or {}
    if anchor.get("file"):
        paths.append(anchor["file"])
    seen = set()
    result = []
    for path in paths:
        clean = re.sub(r":V2.*$", "", path)
        if clean not in seen:
            seen.add(clean)
            result.append(clean)
    return result


def legacy_source_paths(node: dict[str, Any]) -> list[str]:
    paths = []
    for key in ("legacy_sources", "legacy_files"):
        for value in node.get(key, []) or []:
            if isinstance(value, str) and not looks_like_atlas_id(value):
                paths.append(value)
    seen = set()
    result = []
    for path in paths:
        if path not in seen:
            seen.add(path)
            result.append(path)
    return result


def edge_ref(edge: dict[str, Any], direction: str) -> dict[str, Any]:
    other_key = "target" if direction == "out" else "source"
    ref = {
        other_key: edge[other_key],
        "relation": edge.get("relation"),
    }
    if edge.get("status"):
        ref["status"] = edge["status"]
    if edge.get("note"):
        ref["note"] = compact(edge["note"], 120)
    return ref


def related_ids(
    node: dict[str, Any],
    incoming: list[dict[str, Any]],
    outgoing: list[dict[str, Any]],
    node_by_id: dict[str, dict[str, Any]],
) -> dict[str, list[str]]:
    candidates = set()
    for key in (
        "physical_anchors",
        "linked_objects",
        "equation_anchors",
        "anchors_claims",
        "anchors",
        "open_gates",
        "protects",
        "protected_by",
        "guards_against",
        "section_anchors",
    ):
        candidates.update(value for value in node.get(key, []) or [] if isinstance(value, str))
    if isinstance(node.get("parent_node"), str):
        candidates.add(node["parent_node"])
    if isinstance(node.get("file_id"), str):
        candidates.add(node["file_id"])
    for edge in incoming:
        candidates.add(edge["source"])
    for edge in outgoing:
        candidates.add(edge["target"])
    candidates.discard(node["id"])

    buckets = {
        "physical_ids": [],
        "math_ids": [],
        "equation_ids": [],
        "claim_ids": [],
        "open_gate_ids": [],
        "status_firewall_ids": [],
        "source_ids": [],
    }
    for node_id in sorted(candidates):
        other = node_by_id.get(node_id)
        if not other:
            continue
        layer = other.get("layer", "")
        node_type = other.get("type", "")
        if node_id.startswith("PHYS_") or layer in {"physical_ontology", "physical_register"}:
            buckets["physical_ids"].append(node_id)
        if node_id.startswith("MATH_") or layer == "math_object":
            buckets["math_ids"].append(node_id)
        if node_id.startswith("EQ_") or layer == "equation_anchor":
            buckets["equation_ids"].append(node_id)
        if node_id.startswith("CLAIM_") or layer == "claim_theorem":
            buckets["claim_ids"].append(node_id)
        if node_id.startswith("OPEN_") or layer == "open_gate":
            buckets["open_gate_ids"].append(node_id)
        if node_id.startswith("FIREWALL_") or node_type in STATUS_FIREWALL_TYPES:
            buckets["status_firewall_ids"].append(node_id)
        if layer in {"file_anchor", "source_section_anchor"}:
            buckets["source_ids"].append(node_id)
    return buckets


def topic_tags(node: dict[str, Any]) -> list[str]:
    text = graph_text(node)
    topics = []
    checks = {
        "moving_throat": ["moving throat", "pde", "wall", "bdg"],
        "pn_chain": ["1pn", "2pn", "2.5pn", "3pn", "4pn", "5pn"],
        "quadrupole": ["quadrupole", "p0", "stf"],
        "charge": ["charge", "q_eff", "qstar", "q_star", "eta_q"],
        "maxwell": ["maxwell", "zero mode", "mixed"],
        "projection": ["projection", "reduction"],
        "atom": ["atom", "hydrogen"],
        "lepton": ["lepton"],
        "g2": ["g-2", "anomaly"],
    }
    for topic, needles in checks.items():
        if any(needle in text for needle in needles):
            topics.append(f"topic/{topic}")
    return topics


def frontmatter_for(
    graph: dict[str, Any],
    node: dict[str, Any],
    incoming: list[dict[str, Any]],
    outgoing: list[dict[str, Any]],
    node_by_id: dict[str, dict[str, Any]],
    generated_utc: str,
) -> dict[str, Any]:
    title = title_for(node)
    related = related_ids(node, incoming, outgoing, node_by_id)
    tags = [
        "atlas/node",
        f"atlas/{node_folder(node)}",
        f"layer/{sanitize_tag(node.get('layer'))}",
        f"type/{sanitize_tag(node.get('type'))}",
        f"status/{sanitize_tag(node.get('status'))}",
    ]
    tags.extend(topic_tags(node))
    data = {
        "id": node["id"],
        "title": title,
        "type": node.get("type"),
        "layer": node.get("layer"),
        "status": node.get("status"),
        "atlas_version": "obsidian-v0.1",
        "source_graph_version": graph.get("schema_version") or graph.get("release_version"),
        "source_graph_file": "graph/fluid_universe_derivation_atlas_graph.yaml",
        "generated_by": "codex",
        "generated": True,
        "last_generated_utc": generated_utc,
        "summary_short": compact(node.get("meaning") or title, 180),
        "source_kind": node.get("source_kind"),
        "future_paper_needed": bool(node.get("future_paper_needed", False)),
        "source_files": source_paths(node),
        "legacy_sources": legacy_source_paths(node),
        "source_links": [wiki(node_id) for node_id in related["source_ids"]],
        "tex_anchor": node.get("tex_anchor"),
        **related,
        "outgoing_edges": [edge_ref(edge, "out") for edge in outgoing],
        "incoming_edges": [edge_ref(edge, "in") for edge in incoming],
        "tags": sorted(set(tags)),
    }
    return {key: value for key, value in data.items() if value not in (None, [], {})}


def yaml_frontmatter(data: dict[str, Any]) -> str:
    return "---\n" + yaml.safe_dump(data, sort_keys=False, allow_unicode=True, width=1000) + "---"


def math_blocks(node: dict[str, Any]) -> list[str]:
    blocks = []
    if node.get("expression"):
        blocks.append(str(node["expression"]))
    for value in node.get("math", []) or []:
        if isinstance(value, str):
            blocks.append(value)
    return blocks


def display_math(value: str) -> str:
    text = value.strip()
    if text.startswith("$$"):
        return text
    text = text.strip("$")
    return f"$$\n{text}\n$$"


def link_list(ids: list[str]) -> str:
    if not ids:
        return "- none"
    return "\n".join(f"- {wiki(node_id)}" for node_id in ids)


def edge_table(edges: list[dict[str, Any]], direction: str) -> str:
    if not edges:
        return "- none"
    other = "target" if direction == "out" else "source"
    lines = ["| Relation | Node | Note |", "|---|---|---|"]
    for edge in edges:
        lines.append(
            f"| `{edge.get('relation', '')}` | {wiki(edge[other])} | {compact(edge.get('note'), 110)} |"
        )
    return "\n".join(lines)


def source_anchor_section(node: dict[str, Any], related: dict[str, list[str]]) -> str:
    paths = source_paths(node)
    lines = []
    if related["source_ids"]:
        lines.append("### Source anchor notes")
        lines.extend(f"- {wiki(node_id)}" for node_id in related["source_ids"])
    else:
        lines.append("### Source anchor notes")
        lines.append("- No source anchor note recorded.")
    lines.append("")
    if paths:
        lines.append("### Source files")
        lines.extend(f"- `{path}`" for path in paths)
    else:
        lines.append("### Source files")
        lines.append("- No source file path recorded.")

    anchor = node.get("tex_anchor")
    if anchor:
        lines.append("\n### TeX anchor")
        label = anchor.get("nearest_label") or {}
        lines.append(f"- File: `{anchor.get('file')}`")
        lines.append(f"- Line: `{anchor.get('line')}`")
        if anchor.get("heading"):
            lines.append(f"- Heading: {anchor.get('heading')}")
        if label.get("name"):
            lines.append(f"- Nearest label: `{label.get('name')}` at line `{label.get('line')}`")
        lines.append(f"- Match basis: `{anchor.get('match_basis')}`")
        lines.append(f"- Confidence: `{anchor.get('confidence')}`")
        if anchor.get("note"):
            lines.append(f"- Note: {anchor.get('note')}")
    return "\n".join(lines)


def note_body(
    node: dict[str, Any],
    incoming: list[dict[str, Any]],
    outgoing: list[dict[str, Any]],
    related: dict[str, list[str]],
) -> str:
    title = title_for(node)
    blocks = math_blocks(node)
    lines = [
        GENERATED_HEADER,
        "",
        f"# {title}",
        "",
        "> [!warning] Generated note",
        "> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.",
        "",
        f"> **Atlas ID:** `{node['id']}`  ",
        f"> **Status:** `{node.get('status', '')}`  ",
        f"> **Layer:** `{node.get('layer', '')}`  ",
        f"> **Type:** `{node.get('type', '')}`",
        "",
        "## Summary",
        "",
        node.get("meaning") or title,
        "",
    ]

    if node.get("future_paper_needed"):
        lines.extend(
            [
                "> [!note] Future paper needed",
                "> This node is intentionally anchored to a notes/derivation source until a maintained paper exists.",
                "",
            ]
        )

    if node.get("layer") == "claim_theorem" or node.get("type") == "claim":
        lines.extend(["## Claim", "", node.get("meaning") or "Claim text is represented by the graph summary.", ""])
        if "open" in str(node.get("status", "")).lower() or "conditional" in str(node.get("status", "")).lower():
            lines.extend(
                [
                    "## What It Does Not Claim",
                    "",
                    "This generated note preserves the graph status. It should not be read as closing an open gate or promoting a conditional theorem.",
                    "",
                ]
            )

    if node.get("layer") == "open_gate" or node["id"].startswith("OPEN_"):
        lines.extend(
            [
                "## What Remains Open",
                "",
                node.get("meaning") or "The graph marks this as an open gate.",
                "",
                "## What Would Close It",
                "",
                "A source-backed derivation, branch computation, theorem, or paper update must change the graph source of truth before this note can change status.",
                "",
            ]
        )

    if node["id"].startswith("FIREWALL_") or node.get("type") in STATUS_FIREWALL_TYPES:
        lines.extend(
            [
                "## Invalid Inference",
                "",
                node.get("forbidden_target") or node.get("guards_against") or "See the status-firewall rule in the graph metadata.",
                "",
                "## Corrected Inference",
                "",
                node.get("rule") or node.get("meaning") or "Use the graph status and source anchors as the allowed reading.",
                "",
            ]
        )

    lines.extend(
        [
            "## Physical Meaning",
            "",
            node.get("meaning") or "No additional physical interpretation is recorded beyond the graph metadata.",
            "",
            "## Mathematical Role",
            "",
            f"- Layer: `{node.get('layer')}`",
            f"- Type: `{node.get('type')}`",
            f"- Status: `{node.get('status')}`",
        ]
    )
    if node.get("parent_node"):
        related_parent_ids = set(related.get("math_ids", []) + related.get("equation_ids", []) + related.get("claim_ids", []))
        parent_text = wiki(node["parent_node"]) if node["parent_node"] in related_parent_ids else f"`{node['parent_node']}`"
        lines.append(f"- Parent node: {parent_text}")
    if node.get("target"):
        lines.append(f"- Target: `{node['target']}`")
    if node.get("outputs"):
        lines.append("- Outputs: " + ", ".join(f"`{value}`" for value in node["outputs"]))
    lines.append("")

    if blocks:
        lines.extend(["## Equation", ""])
        for block in blocks:
            lines.extend([display_math(block), ""])
        lines.extend(
            [
                "## Variable Dictionary",
                "",
                "The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.",
                "",
            ]
        )

    lines.extend(
        [
            "## Atlas Links",
            "",
            "### Related physical nodes",
            link_list(related["physical_ids"]),
            "",
            "### Related math nodes",
            link_list(related["math_ids"]),
            "",
            "### Related equations",
            link_list(related["equation_ids"]),
            "",
            "### Related claims",
            link_list(related["claim_ids"]),
            "",
            "### Open gates",
            link_list(related["open_gate_ids"]),
            "",
            "### Status firewalls",
            link_list(related["status_firewall_ids"]),
            "",
            "### Source anchors",
            link_list(related["source_ids"]),
            "",
            "## Outgoing Edges",
            "",
            edge_table(outgoing, "out"),
            "",
            "## Incoming Edges",
            "",
            edge_table(incoming, "in"),
            "",
            "## Source Anchors",
            "",
            source_anchor_section(node, related),
            "",
            "## AI Maintenance Notes",
            "",
            "- Treat this file as generated read-only presentation material.",
            "- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.",
            "- Do not use this note to upgrade statuses or weaken firewalls.",
            "",
        ]
    )
    return "\n".join(lines)


def build_indexes(graph: dict[str, Any]) -> tuple[dict[str, dict[str, Any]], dict[str, list[dict[str, Any]]], dict[str, list[dict[str, Any]]]]:
    node_by_id = {node["id"]: node for node in graph["nodes"]}
    incoming: dict[str, list[dict[str, Any]]] = defaultdict(list)
    outgoing: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for edge in graph["edges"]:
        outgoing[edge["source"]].append(edge)
        incoming[edge["target"]].append(edge)
    for edges in incoming.values():
        edges.sort(key=lambda edge: (edge.get("relation", ""), edge.get("source", "")))
    for edges in outgoing.values():
        edges.sort(key=lambda edge: (edge.get("relation", ""), edge.get("target", "")))
    return node_by_id, incoming, outgoing


def write_node_notes(graph: dict[str, Any], generated_utc: str) -> dict[str, str]:
    node_by_id, incoming, outgoing = build_indexes(graph)
    paths: dict[str, str] = {}
    for node in graph["nodes"]:
        inc = incoming.get(node["id"], [])
        out = outgoing.get(node["id"], [])
        related = related_ids(node, inc, out, node_by_id)
        frontmatter = frontmatter_for(graph, node, inc, out, node_by_id, generated_utc)
        content = yaml_frontmatter(frontmatter) + "\n\n" + note_body(node, inc, out, related)
        path = node_path(node)
        path.write_text(content, encoding="utf-8")
        paths[node["id"]] = vault_path(path)
    return paths


def write_exports(graph: dict[str, Any], note_paths: dict[str, str]) -> None:
    graph_export = dict(graph)
    for node in graph_export["nodes"]:
        node["obsidian_note"] = note_paths[node["id"]]
    (EXPORT_DIR / "atlas_graph.yaml").write_text(
        yaml.safe_dump(graph_export, sort_keys=False, allow_unicode=True, width=1000),
        encoding="utf-8",
    )
    (EXPORT_DIR / "atlas_graph.json").write_text(
        json.dumps(graph_export, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    with (EXPORT_DIR / "atlas_nodes.tsv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(["id", "label", "layer", "type", "status", "note", "source_kind", "tex_file", "tex_line"])
        for node in graph["nodes"]:
            anchor = node.get("tex_anchor") or {}
            writer.writerow(
                [
                    node["id"],
                    node.get("label", ""),
                    node.get("layer", ""),
                    node.get("type", ""),
                    node.get("status", ""),
                    note_paths[node["id"]],
                    node.get("source_kind", ""),
                    anchor.get("file", ""),
                    anchor.get("line", ""),
                ]
            )

    with (EXPORT_DIR / "atlas_edges.tsv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(["source", "relation", "target", "status", "note"])
        for edge in graph["edges"]:
            writer.writerow([edge["source"], edge.get("relation", ""), edge["target"], edge.get("status", ""), edge.get("note", "")])


def nodes_matching_canvas_spec(
    graph: dict[str, Any],
    node_by_id: dict[str, dict[str, Any]],
    incoming: dict[str, list[dict[str, Any]]],
    outgoing: dict[str, list[dict[str, Any]]],
    spec: dict[str, Any],
    max_nodes: int = 34,
) -> list[str]:
    selected = set(node_id for node_id in spec.get("seeds", []) if node_id in node_by_id)
    layers = set(spec.get("layers", []) or [])
    types = set(spec.get("types", []) or [])
    terms = [term.lower() for term in spec.get("terms", [])]
    for node in graph["nodes"]:
        text = graph_text(node)
        if layers and node.get("layer") in layers:
            selected.add(node["id"])
        if types and node.get("type") in types:
            selected.add(node["id"])
        if terms and any(all(part in text for part in term.split()) for term in terms):
            selected.add(node["id"])

    queue = deque((node_id, 0) for node_id in list(selected))
    while queue and len(selected) < max_nodes:
        node_id, depth = queue.popleft()
        if depth >= 1:
            continue
        adjacent = outgoing.get(node_id, []) + incoming.get(node_id, [])
        for edge in adjacent:
            other = edge["target"] if edge["source"] == node_id else edge["source"]
            if other not in selected:
                selected.add(other)
                queue.append((other, depth + 1))
            if len(selected) >= max_nodes:
                break

    return sorted(selected, key=lambda node_id: (node_by_id[node_id].get("layer", ""), node_id))[:max_nodes]


def write_canvas_files(graph: dict[str, Any], note_paths: dict[str, str]) -> None:
    node_by_id, incoming, outgoing = build_indexes(graph)
    layer_order = {layer: index for index, layer in enumerate(sorted({node.get("layer", "") for node in graph["nodes"]}))}
    for filename, spec in CANVAS_SPECS.items():
        node_ids = nodes_matching_canvas_spec(graph, node_by_id, incoming, outgoing, spec)
        canvas_layers = sorted(
            {node_by_id[node_id].get("layer", "") for node_id in node_ids},
            key=lambda layer: layer_order.get(layer, 999),
        )
        column_by_layer = {layer: index for index, layer in enumerate(canvas_layers)}
        column_counts: defaultdict[int, int] = defaultdict(int)
        nodes = []
        for node_id in node_ids:
            node = node_by_id[node_id]
            column = column_by_layer[node.get("layer", "")]
            row = column_counts[column]
            column_counts[column] += 1
            nodes.append(
                {
                    "id": node_id,
                    "type": "file",
                    "file": note_paths[node_id],
                    "x": column * 380,
                    "y": row * 240,
                    "width": 330,
                    "height": 190,
                }
            )
        selected = set(node_ids)
        edges = []
        edge_index = 0
        for edge in graph["edges"]:
            if edge["source"] not in selected or edge["target"] not in selected:
                continue
            edges.append(
                {
                    "id": f"e{edge_index}",
                    "fromNode": edge["source"],
                    "fromSide": "right",
                    "toNode": edge["target"],
                    "toSide": "left",
                    "label": edge.get("relation", ""),
                }
            )
            edge_index += 1
        canvas = {"nodes": nodes, "edges": edges}
        (CANVAS_DIR / filename).write_text(json.dumps(canvas, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_base(path: Path, folder: str, view_name: str, order: list[str]) -> None:
    data = {
        "filters": {"and": [f'file.inFolder("{folder}")']},
        "views": [
            {
                "type": "table",
                "name": view_name,
                "order": order,
            }
        ],
    }
    path.write_text(yaml.safe_dump(data, sort_keys=False, allow_unicode=True), encoding="utf-8")


def write_bases() -> None:
    common = ["file.name", "status", "summary_short", "source_links", "source_files"]
    write_base(BASE_DIR / "claims.base", "nodes/claims", "Claims", ["file.name", "layer", "status", "claim_strength", "summary_short", "source_links", "source_files"])
    write_base(BASE_DIR / "open_gates.base", "nodes/open_gates", "Open gates", common + ["claim_ids", "equation_ids"])
    write_base(BASE_DIR / "physical_objects.base", "nodes/physical", "Physical objects", ["file.name", "type", "status", "summary_short", "math_ids", "equation_ids"])
    write_base(BASE_DIR / "status_firewalls.base", "nodes/status_firewalls", "Status firewalls", common + ["claim_ids", "open_gate_ids"])
    write_base(BASE_DIR / "source_crosswalk.base", "nodes/sources", "Source crosswalk", ["file.name", "source_kind", "status", "summary_short", "source_links", "source_files", "tex_anchor"])
    write_base(BASE_DIR / "equations.base", "nodes/equations", "Equations", ["file.name", "status", "summary_short", "claim_ids", "source_links", "source_files", "tex_anchor"])


def static_node_list(nodes: list[dict[str, Any]], limit: int = 40) -> str:
    if not nodes:
        return "- none"
    lines = []
    for node in nodes[:limit]:
        lines.append(f"- {wiki(node['id'])} — `{node.get('status', '')}` — {compact(node.get('meaning') or node.get('label'), 110)}")
    if len(nodes) > limit:
        lines.append(f"- ... {len(nodes) - limit} more")
    return "\n".join(lines)


def write_views(graph: dict[str, Any]) -> None:
    nodes = graph["nodes"]
    by_layer = defaultdict(list)
    for node in nodes:
        by_layer[node.get("layer", "")].append(node)
    for layer_nodes in by_layer.values():
        layer_nodes.sort(key=lambda node: node["id"])

    open_nodes = [node for node in nodes if node.get("layer") == "open_gate" or "open" in str(node.get("status", "")).lower() or node.get("future_paper_needed")]
    claim_nodes = [node for node in nodes if node.get("layer") == "claim_theorem"]
    source_nodes = [node for node in nodes if node.get("layer") in {"file_anchor", "source_section_anchor"}]
    firewall_nodes = [node for node in nodes if node["id"].startswith("FIREWALL_") or node.get("type") in STATUS_FIREWALL_TYPES]

    views = {
        "00_home.md": [
            "# Fluid Universe Atlas",
            "",
            GENERATED_HEADER,
            "",
            "This vault is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`.",
            "",
            f"- Nodes: `{len(graph['nodes'])}`",
            f"- Edges: `{len(graph['edges'])}`",
            f"- Source graph version: `{graph.get('schema_version')}`",
            "",
            "## Dashboards",
            "",
            "- [[01_open_gates]]",
            "- [[02_claim_register]]",
            "- [[03_physics_to_math]]",
            "- [[04_moving_throat_pde]]",
            "- [[05_pn_chain]]",
            "- [[06_charge_ontology_firewall]]",
            "- [[07_status_firewalls]]",
            "- [[08_coverage_audit]]",
            "- [[09_source_crosswalk]]",
            "- [[10_query_validation]]",
            "",
            "## Bases",
            "",
            "- ![[claims.base]]",
            "- ![[open_gates.base]]",
            "- ![[source_crosswalk.base]]",
        ],
        "01_open_gates.md": ["# Open Gates", "", GENERATED_HEADER, "", "![[open_gates.base]]", "", static_node_list(open_nodes, 120)],
        "02_claim_register.md": ["# Claim Register", "", GENERATED_HEADER, "", "![[claims.base]]", "", static_node_list(claim_nodes, 120)],
        "03_physics_to_math.md": ["# Physics To Math", "", GENERATED_HEADER, "", static_node_list(by_layer["physical_ontology"] + by_layer["math_object"] + by_layer["equation_anchor"], 140)],
        "04_moving_throat_pde.md": ["# Moving Throat PDE", "", GENERATED_HEADER, "", static_node_list([node for node in nodes if "moving_throat" in " ".join(topic_tags(node)) or "pde" in graph_text(node)], 120)],
        "05_pn_chain.md": ["# PN Chain", "", GENERATED_HEADER, "", static_node_list([node for node in nodes if any(term in graph_text(node) for term in ["1pn", "2pn", "2.5pn", "3pn", "4pn", "5pn"])], 140)],
        "06_charge_ontology_firewall.md": ["# Charge Ontology Firewall", "", GENERATED_HEADER, "", static_node_list([node for node in nodes if any(term in graph_text(node) for term in ["charge", "q_eff", "qstar", "eta_q", "circulation", "kappa_rho"])], 120)],
        "07_status_firewalls.md": ["# Status Firewalls", "", GENERATED_HEADER, "", "![[status_firewalls.base]]", "", static_node_list(firewall_nodes, 120)],
        "08_coverage_audit.md": ["# Coverage Audit", "", GENERATED_HEADER, "", "## Nodes by layer", ""],
        "09_source_crosswalk.md": ["# Source Crosswalk", "", GENERATED_HEADER, "", "![[source_crosswalk.base]]", "", static_node_list(source_nodes, 140)],
        "10_query_validation.md": ["# Query Validation", "", GENERATED_HEADER, "", static_node_list(by_layer["query_validation"], 120)],
    }
    layer_counts = Counter(node.get("layer", "unknown") for node in nodes)
    views["08_coverage_audit.md"].extend(f"- `{layer}`: {count}" for layer, count in sorted(layer_counts.items()))
    anchored = sum(1 for node in nodes if node.get("tex_anchor"))
    views["08_coverage_audit.md"].extend(["", f"TeX-anchored nodes: `{anchored}`"])

    for filename, parts in views.items():
        (VIEW_DIR / filename).write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_templates(generated_utc: str) -> None:
    templates = {
        "node_template.md": "Generic generated node template",
        "claim_node_template.md": "Generated claim node template",
        "equation_node_template.md": "Generated equation node template",
        "open_gate_template.md": "Generated open gate template",
        "status_firewall_template.md": "Generated status firewall template",
        "source_anchor_template.md": "Generated source anchor template",
    }
    for filename, title in templates.items():
        content = f"""---
id: NODE_ID
title: {title}
generated: true
generated_by: codex
last_generated_utc: {generated_utc}
---

# {title}

{GENERATED_HEADER}

Templates are generated placeholders. The node generator writes concrete notes from the graph YAML.
"""
        (TEMPLATE_DIR / filename).write_text(content, encoding="utf-8")


def main() -> None:
    global ATLAS_DIR, NODE_DIR, VIEW_DIR, BASE_DIR, CANVAS_DIR, EXPORT_DIR, TEMPLATE_DIR

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=GRAPH_YAML, help="Input graph YAML")
    parser.add_argument("--atlas", type=Path, default=ATLAS_DIR, help="Atlas vault root")
    args = parser.parse_args()

    ATLAS_DIR = args.atlas
    NODE_DIR = ATLAS_DIR / "nodes"
    VIEW_DIR = ATLAS_DIR / "views"
    BASE_DIR = ATLAS_DIR / "bases"
    CANVAS_DIR = ATLAS_DIR / "canvas"
    EXPORT_DIR = ATLAS_DIR / "exports"
    TEMPLATE_DIR = ATLAS_DIR / "templates"

    clean_previous_generated_files()
    ensure_dirs()
    graph = load_graph(args.graph)
    generated_utc = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")

    note_paths = write_node_notes(graph, generated_utc)
    write_exports(graph, note_paths)
    write_bases()
    write_views(graph)
    write_canvas_files(graph, note_paths)
    write_templates(generated_utc)
    write_generated_manifest()

    print(f"generated notes: {len(note_paths)}")
    print(f"generated edges: {len(graph['edges'])}")
    print(f"atlas root: {ATLAS_DIR}")


if __name__ == "__main__":
    main()
