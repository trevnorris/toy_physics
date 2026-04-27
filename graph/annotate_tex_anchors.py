#!/usr/bin/env python3
"""Annotate graph nodes with TeX-local anchors.

This script reads the YAML graph, parses TeX headings and labels, then writes
`tex_anchor` metadata back to YAML. It does not edit TeX papers.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any

import yaml


ROOT = Path(__file__).resolve().parents[1]
GRAPH_DIR = Path(__file__).resolve().parent
DEFAULT_YAML = GRAPH_DIR / "fluid_universe_derivation_atlas_graph.yaml"

HEADING_COMMANDS = {
    "part": 0,
    "chapter": 0,
    "section": 1,
    "subsection": 2,
    "subsubsection": 3,
    "paragraph": 4,
}
LABEL_RE = re.compile(r"\\label\{([^}]+)\}")
TOKEN_RE = re.compile(r"[a-z0-9]+")
MANUAL_EQUIVALENT_HEADINGS = {
    "SEC_1PN_BRIDGE_DOWNSTREAM": {
        "file": "research/4d_1pn_bridge/paper/4d_1pn_bridge.tex",
        "heading_contains": "Conclusions",
        "note": "No exact TeX heading matches the legacy summary heading `Takeaway for downstream reasoning`; the paper-level downstream takeaway is anchored to Conclusions.",
    }
}


@dataclass
class Label:
    name: str
    line: int
    heading_index: int | None = None


@dataclass
class Heading:
    command: str
    level: int
    title: str
    clean_title: str
    line: int
    labels: list[Label] = field(default_factory=list)


@dataclass
class TeXFile:
    path: str
    headings: list[Heading]
    labels: list[Label]
    line_count: int


def fail(message: str) -> None:
    print(f"FAIL: {message}", file=sys.stderr)
    raise SystemExit(1)


def load_yaml(path: Path) -> dict[str, Any]:
    try:
        return yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        fail(f"{path} is not valid YAML: {exc}")


def write_yaml(path: Path, graph: dict[str, Any]) -> None:
    path.write_text(
        yaml.safe_dump(graph, sort_keys=False, allow_unicode=True, width=1000),
        encoding="utf-8",
    )


def command_arg(line: str, command: str) -> str | None:
    start = line.find(f"\\{command}")
    if start < 0:
        return None
    index = start + len(command) + 1
    if index < len(line) and line[index] == "*":
        index += 1
    while index < len(line) and line[index].isspace():
        index += 1
    if index < len(line) and line[index] == "[":
        depth = 1
        index += 1
        while index < len(line) and depth:
            if line[index] == "[":
                depth += 1
            elif line[index] == "]":
                depth -= 1
            index += 1
    while index < len(line) and line[index].isspace():
        index += 1
    if index >= len(line) or line[index] != "{":
        return None

    depth = 1
    index += 1
    out: list[str] = []
    while index < len(line) and depth:
        char = line[index]
        if char == "{":
            depth += 1
            out.append(char)
        elif char == "}":
            depth -= 1
            if depth:
                out.append(char)
        else:
            out.append(char)
        index += 1
    return "".join(out) if depth == 0 else None


def strip_tex(value: Any) -> str:
    text = str(value or "")
    previous = None
    while previous != text:
        previous = text
        text = re.sub(r"\\texorpdfstring\{([^{}]*)\}\{([^{}]*)\}", r"\2", text)
    replacements = {
        r"\oplus": " plus ",
        r"\mathrm": " ",
        r"\texttt": " ",
        r"\text": " ",
        r"\emph": " ",
        r"\rm": " ",
        r"\left": " ",
        r"\right": " ",
        r"\big": " ",
        r"\quad": " ",
        r"\,": " ",
        r"\_": "_",
        "--": " ",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    text = re.sub(r"\\[a-zA-Z]+", " ", text)
    text = re.sub(r"[{}$^_~\\]", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def normalized(value: Any) -> str:
    return " ".join(TOKEN_RE.findall(strip_tex(value).lower()))


def tokens(value: Any) -> set[str]:
    return set(TOKEN_RE.findall(strip_tex(value).lower()))


def parse_tex(path: str) -> TeXFile:
    full_path = ROOT / path
    lines = full_path.read_text(encoding="utf-8").splitlines()
    headings: list[Heading] = []
    labels: list[Label] = []
    current_heading_index: int | None = None

    for index, line in enumerate(lines, 1):
        for command, level in HEADING_COMMANDS.items():
            title = command_arg(line, command)
            if title is None:
                continue
            current_heading_index = len(headings)
            headings.append(
                Heading(
                    command=command,
                    level=level,
                    title=title,
                    clean_title=strip_tex(title),
                    line=index,
                )
            )
            break

        for label_name in LABEL_RE.findall(line):
            label = Label(label_name, index, current_heading_index)
            labels.append(label)
            if current_heading_index is not None:
                headings[current_heading_index].labels.append(label)

    return TeXFile(path=path, headings=headings, labels=labels, line_count=len(lines))


def tex_paths_for_node(node: dict[str, Any]) -> list[str]:
    paths: list[str] = []
    for key in ("file",):
        value = node.get(key)
        if isinstance(value, str):
            paths.append(value)
    for key in ("files", "sources", "canonical_target_files", "source_note_files"):
        for value in node.get(key, []) or []:
            if isinstance(value, str):
                paths.append(value)
    seen: set[str] = set()
    tex_paths = []
    for path in paths:
        clean = re.sub(r":V2.*$", "", path)
        if clean.endswith(".tex") and (ROOT / clean).exists() and clean not in seen:
            seen.add(clean)
            tex_paths.append(clean)
    return tex_paths


def score_text(query: str, candidate: str) -> float:
    query_norm = normalized(query)
    candidate_norm = normalized(candidate)
    if not query_norm or not candidate_norm:
        return 0.0
    ratio = SequenceMatcher(None, query_norm, candidate_norm).ratio()
    query_tokens = set(query_norm.split())
    candidate_tokens = set(candidate_norm.split())
    overlap = len(query_tokens & candidate_tokens) / max(1, len(query_tokens))
    containment_bonus = 0.15 if query_tokens and query_tokens <= candidate_tokens else 0.0
    return min(1.0, 0.58 * ratio + 0.42 * overlap + containment_bonus)


def heading_anchor(tex_file: TeXFile, heading: Heading, *, basis: str, score: float) -> dict[str, Any]:
    nearest = None
    nearby = []
    if heading.labels:
        nearest_label = min(heading.labels, key=lambda label: abs(label.line - heading.line))
        nearest = {"name": nearest_label.name, "line": nearest_label.line}
        nearby = [{"name": label.name, "line": label.line} for label in heading.labels[:10]]
    confidence = "high" if score >= 0.72 else "medium" if score >= 0.50 else "low"
    return {
        "file": tex_file.path,
        "line": heading.line,
        "heading_level": heading.command,
        "heading": heading.clean_title,
        "nearest_label": nearest,
        "nearby_labels": nearby,
        "match_basis": basis,
        "match_score": round(score, 3),
        "confidence": confidence,
    }


def label_anchor(tex_file: TeXFile, label: Label, *, basis: str, score: float) -> dict[str, Any]:
    heading = tex_file.headings[label.heading_index] if label.heading_index is not None else None
    anchor = {
        "file": tex_file.path,
        "line": label.line,
        "heading_level": heading.command if heading else None,
        "heading": heading.clean_title if heading else None,
        "nearest_label": {"name": label.name, "line": label.line},
        "nearby_labels": [{"name": label.name, "line": label.line}],
        "match_basis": basis,
        "match_score": round(score, 3),
        "confidence": "high" if score >= 0.72 else "medium" if score >= 0.50 else "low",
    }
    return {key: value for key, value in anchor.items() if value is not None}


def best_heading_match(tex_file: TeXFile, queries: list[str]) -> tuple[float, Heading | None, str]:
    best_score = 0.0
    best_heading = None
    best_query = ""
    for query in queries:
        if not normalized(query):
            continue
        for heading in tex_file.headings:
            score = score_text(query, heading.clean_title)
            if score > best_score:
                best_score = score
                best_heading = heading
                best_query = query
    return best_score, best_heading, best_query


def best_label_match(tex_file: TeXFile, queries: list[str]) -> tuple[float, Label | None, str]:
    best_score = 0.0
    best_label = None
    best_query = ""
    for query in queries:
        query_terms = tokens(query)
        if not query_terms:
            continue
        for label in tex_file.labels:
            label_terms = tokens(label.name.replace(":", " ").replace("-", " "))
            if not label_terms:
                continue
            overlap = len(query_terms & label_terms) / max(1, len(query_terms))
            ratio = SequenceMatcher(None, normalized(query), normalized(label.name)).ratio()
            score = 0.45 * ratio + 0.55 * overlap
            if score > best_score:
                best_score = score
                best_label = label
                best_query = query
    return best_score, best_label, best_query


def graph_edges_by_target(graph: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    incoming: dict[str, list[dict[str, Any]]] = {}
    for edge in graph["edges"]:
        incoming.setdefault(edge["target"], []).append(edge)
    return incoming


def node_queries(node: dict[str, Any]) -> list[str]:
    queries = [
        node.get("label"),
        node.get("query"),
        node.get("meaning"),
        node.get("insert_after"),
        node.get("expression"),
    ]
    return [str(query) for query in queries if query]


def annotate_source_section(
    node: dict[str, Any],
    tex_files: dict[str, TeXFile],
) -> bool:
    override = MANUAL_EQUIVALENT_HEADINGS.get(node["id"])
    if override:
        tex_file = tex_files.get(override["file"])
        if tex_file:
            needle = normalized(override["heading_contains"])
            for heading in tex_file.headings:
                if needle and needle in normalized(heading.clean_title):
                    anchor = heading_anchor(
                        tex_file,
                        heading,
                        basis="manual_equivalent_heading",
                        score=0.5,
                    )
                    anchor["confidence"] = "low"
                    anchor["note"] = override["note"]
                    node["tex_anchor"] = anchor
                    return True

    paths = tex_paths_for_node(node)
    if not paths:
        return False
    best: tuple[float, dict[str, Any] | None] = (0.0, None)
    for path in paths:
        tex_file = tex_files[path]
        score, heading, _query = best_heading_match(tex_file, node_queries(node))
        if heading and score > best[0]:
            best = (score, heading_anchor(tex_file, heading, basis="semantic_heading_match", score=score))
    if best[1] and best[0] >= 0.38:
        node["tex_anchor"] = best[1]
        return True
    return False


def inherited_anchor(
    node: dict[str, Any],
    incoming: dict[str, list[dict[str, Any]]],
    node_by_id: dict[str, dict[str, Any]],
) -> dict[str, Any] | None:
    preferred_relations = {
        "ANCHORS_CLAIM_SECTION",
        "SUPPORTS_CLAIM",
        "ANCHORED_IN",
        "CONTAINS_EQUATION",
    }
    candidates = []
    for edge in incoming.get(node["id"], []):
        if edge.get("relation") not in preferred_relations:
            continue
        source = node_by_id.get(edge["source"])
        if not source or "tex_anchor" not in source:
            continue
        anchor = dict(source["tex_anchor"])
        anchor["match_basis"] = f"graph_edge:{edge['relation']}"
        anchor["source_anchor_node"] = source["id"]
        anchor["confidence"] = "medium" if anchor.get("confidence") == "high" else anchor.get("confidence", "medium")
        candidates.append((source.get("layer") == "source_section_anchor", anchor))
    if not candidates:
        return None
    candidates.sort(key=lambda item: (not item[0], item[1].get("line") or 0))
    return candidates[0][1]


def annotate_direct_match(
    node: dict[str, Any],
    tex_files: dict[str, TeXFile],
) -> bool:
    if node.get("layer") in {"file_anchor", "paper_backlink"}:
        return False
    paths = tex_paths_for_node(node)
    if not paths:
        return False
    queries = node_queries(node)
    if not queries:
        return False

    best_anchor: dict[str, Any] | None = None
    best_score = 0.0
    for path in paths:
        tex_file = tex_files[path]
        heading_score, heading, _query = best_heading_match(tex_file, queries)
        label_score, label, _label_query = best_label_match(tex_file, queries)
        if label and label_score >= heading_score and label_score > best_score:
            best_score = label_score
            best_anchor = label_anchor(tex_file, label, basis="semantic_label_match", score=label_score)
        if heading and heading_score > best_score:
            best_score = heading_score
            best_anchor = heading_anchor(tex_file, heading, basis="semantic_heading_match", score=heading_score)

    if best_anchor and best_score >= 0.50:
        node["tex_anchor"] = best_anchor
        return True
    return False


def annotate_paper_backlink(
    node: dict[str, Any],
    node_by_id: dict[str, dict[str, Any]],
) -> bool:
    if node.get("layer") != "paper_backlink":
        return False
    for anchor_id in reversed(node.get("section_anchors", []) or []):
        source = node_by_id.get(anchor_id)
        if not source or "tex_anchor" not in source:
            continue
        anchor = dict(source["tex_anchor"])
        anchor["match_basis"] = "backlink_section_anchor"
        anchor["source_anchor_node"] = source["id"]
        anchor["confidence"] = "medium" if anchor.get("confidence") == "high" else anchor.get("confidence", "medium")
        node["tex_anchor"] = anchor
        return True
    return False


def annotate_file_anchor(node: dict[str, Any], tex_files: dict[str, TeXFile]) -> bool:
    if node.get("layer") != "file_anchor":
        return False
    paths = tex_paths_for_node(node)
    if not paths:
        return False
    tex_file = tex_files[paths[0]]
    first_heading = tex_file.headings[0] if tex_file.headings else None
    if first_heading:
        node["tex_anchor"] = heading_anchor(tex_file, first_heading, basis="file_first_heading", score=1.0)
    else:
        node["tex_anchor"] = {
            "file": tex_file.path,
            "line": 1,
            "heading": None,
            "nearest_label": None,
            "nearby_labels": [{"name": label.name, "line": label.line} for label in tex_file.labels[:10]],
            "match_basis": "file_start",
            "match_score": 1.0,
            "confidence": "high",
        }
    return True


def annotate_nodes(graph: dict[str, Any], tex_files: dict[str, TeXFile]) -> dict[str, int]:
    for node in graph["nodes"]:
        node.pop("tex_anchor", None)

    stats = {
        "source_section": 0,
        "inherited": 0,
        "backlink": 0,
        "direct": 0,
        "file_anchor": 0,
        "missing_source_section": 0,
    }

    for node in graph["nodes"]:
        if node.get("layer") == "source_section_anchor" and node.get("source_kind") == "paper":
            if annotate_source_section(node, tex_files):
                stats["source_section"] += 1
            else:
                stats["missing_source_section"] += 1

    incoming = graph_edges_by_target(graph)
    node_by_id = {node["id"]: node for node in graph["nodes"]}

    for node in graph["nodes"]:
        if "tex_anchor" in node:
            continue
        anchor = inherited_anchor(node, incoming, node_by_id)
        if anchor:
            node["tex_anchor"] = anchor
            stats["inherited"] += 1

    for node in graph["nodes"]:
        if "tex_anchor" in node:
            continue
        if annotate_paper_backlink(node, node_by_id):
            stats["backlink"] += 1

    for node in graph["nodes"]:
        if "tex_anchor" in node:
            continue
        if annotate_file_anchor(node, tex_files):
            stats["file_anchor"] += 1

    for node in graph["nodes"]:
        if "tex_anchor" in node:
            continue
        if annotate_direct_match(node, tex_files):
            stats["direct"] += 1

    return stats


def used_tex_paths(graph: dict[str, Any]) -> list[str]:
    paths: set[str] = set()
    for node in graph["nodes"]:
        paths.update(tex_paths_for_node(node))
    return sorted(paths)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-yaml", type=Path, default=DEFAULT_YAML, help="Primary graph YAML path")
    parser.add_argument("--check", action="store_true", help="Fail if anchors would change")
    parser.add_argument("--summary-json", type=Path, help="Optional path to write summary JSON")
    args = parser.parse_args()

    graph = load_yaml(args.graph_yaml)
    paths = used_tex_paths(graph)
    tex_files = {path: parse_tex(path) for path in paths}

    before = json.dumps(graph, sort_keys=True, ensure_ascii=False)
    stats = annotate_nodes(graph, tex_files)
    after = json.dumps(graph, sort_keys=True, ensure_ascii=False)

    anchored = sum(1 for node in graph["nodes"] if node.get("tex_anchor"))
    summary = {
        "tex_files_scanned": len(tex_files),
        "nodes_with_tex_anchor": anchored,
        **stats,
    }

    if args.summary_json:
        args.summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if args.check:
        if before != after:
            fail("TeX anchors are stale; run python3 graph/annotate_tex_anchors.py")
        print("OK: TeX anchors are current")
        return

    write_yaml(args.graph_yaml, graph)
    print(f"wrote {args.graph_yaml}")
    for key, value in summary.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
