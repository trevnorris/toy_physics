#!/usr/bin/env python3
"""Query focused subsets of the Fluid Universe derivation atlas.

The tool is intentionally read-only unless --write-subgraph is provided. Subset
JSON files use the same basic shape as the main graph and can be passed to
render_graph.py for a focused visualization.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict, deque
from pathlib import Path
from typing import Any, Iterable

import yaml


GRAPH_DIR = Path(__file__).resolve().parent
DEFAULT_GRAPH = GRAPH_DIR / "fluid_universe_derivation_atlas_graph.yaml"


OPEN_RE = re.compile(r"(^|[_\W])open([_\W]|$)", re.IGNORECASE)


def flatten_text(value: Any) -> Iterable[str]:
    if value is None:
        return
    if isinstance(value, dict):
        for key, child in value.items():
            yield str(key)
            yield from flatten_text(child)
        return
    if isinstance(value, list):
        for child in value:
            yield from flatten_text(child)
        return
    yield str(value)


def compact(value: Any, limit: int = 130) -> str:
    text = " ".join(str(value or "").split())
    return text if len(text) <= limit else f"{text[:limit - 3]}..."


def normalize_terms(query: str) -> list[str]:
    return [term.lower() for term in re.split(r"\s+", query.strip()) if term]


def arg_text(value: str | list[str] | None) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        return " ".join(value)
    return value


def edge_key(edge: dict[str, Any]) -> tuple[str, str, str, str]:
    return (
        edge.get("source", ""),
        edge.get("target", ""),
        edge.get("relation", ""),
        edge.get("note", ""),
    )


class GraphIndex:
    def __init__(self, graph: dict[str, Any]) -> None:
        self.graph = graph
        self.nodes = graph["nodes"]
        self.edges = graph["edges"]
        self.node_by_id = {node["id"]: node for node in self.nodes}
        self.out_edges: dict[str, list[dict[str, Any]]] = defaultdict(list)
        self.in_edges: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for edge in self.edges:
            self.out_edges[edge["source"]].append(edge)
            self.in_edges[edge["target"]].append(edge)
        self.node_text = {node["id"]: " ".join(flatten_text(node)).lower() for node in self.nodes}

    @classmethod
    def load(cls, path: Path) -> "GraphIndex":
        text = path.read_text(encoding="utf-8")
        if path.suffix.lower() in {".yaml", ".yml"}:
            return cls(yaml.safe_load(text))
        return cls(json.loads(text))

    def require_node(self, node_id: str) -> dict[str, Any]:
        node = self.node_by_id.get(node_id)
        if not node:
            fail(f"unknown node id: {node_id}")
        return node

    def score_node(self, node: dict[str, Any], terms: list[str]) -> int:
        text = self.node_text[node["id"]]
        if not all(term in text for term in terms):
            return 0
        score = 1
        node_id = node["id"].lower()
        label = str(node.get("label", "")).lower()
        path_text = " ".join(
            flatten_text(
                {
                    "file": node.get("file"),
                    "files": node.get("files"),
                    "sources": node.get("sources"),
                    "legacy_sources": node.get("legacy_sources"),
                    "legacy_file": node.get("legacy_file"),
                }
            )
        ).lower()
        for term in terms:
            if term == node_id:
                score += 100
            if term in node_id:
                score += 25
            if term in label:
                score += 15
            if term in str(node.get("layer", "")).lower():
                score += 8
            if term in str(node.get("type", "")).lower():
                score += 8
            if term in str(node.get("status", "")).lower():
                score += 8
            if term in path_text:
                score += 7
        return score

    def search(
        self,
        query: str,
        *,
        layer: str | None = None,
        node_type: str | None = None,
        status: str | None = None,
        source_kind: str | None = None,
        limit: int = 20,
    ) -> list[tuple[int, dict[str, Any]]]:
        terms = normalize_terms(query)
        if not terms:
            return []
        matches: list[tuple[int, dict[str, Any]]] = []
        for node in self.nodes:
            if layer and node.get("layer") != layer:
                continue
            if node_type and node.get("type") != node_type:
                continue
            if status and status.lower() not in str(node.get("status", "")).lower():
                continue
            if source_kind and node.get("source_kind") != source_kind:
                continue
            score = self.score_node(node, terms)
            if score:
                matches.append((score, node))
        matches.sort(key=lambda item: (-item[0], item[1]["id"]))
        return matches[:limit]

    def source_resolution_matches(self, query: str) -> list[dict[str, Any]]:
        terms = normalize_terms(query)
        results = []
        for entry in self.graph.get("source_resolution", []):
            text = " ".join(flatten_text(entry)).lower()
            if all(term in text for term in terms):
                results.append(entry)
        return results

    def source_nodes(self, query: str) -> list[dict[str, Any]]:
        terms = normalize_terms(query)
        results = []
        for node in self.nodes:
            text = " ".join(
                flatten_text(
                    {
                        "file": node.get("file"),
                        "files": node.get("files"),
                        "sources": node.get("sources"),
                        "legacy_sources": node.get("legacy_sources"),
                        "legacy_file": node.get("legacy_file"),
                        "source_note_files": node.get("source_note_files"),
                        "canonical_target_files": node.get("canonical_target_files"),
                    }
                )
            ).lower()
            if text and all(term in text for term in terms):
                results.append(node)
        results.sort(key=lambda node: (node.get("layer", ""), node["id"]))
        return results

    def neighborhood(
        self,
        seeds: Iterable[str],
        *,
        depth: int,
        direction: str,
        relations: set[str] | None = None,
        max_nodes: int = 200,
    ) -> tuple[set[str], list[dict[str, Any]]]:
        seen = set(seeds)
        queue = deque((seed, 0) for seed in seen)
        selected_edges: dict[tuple[str, str, str, str], dict[str, Any]] = {}

        while queue:
            node_id, current_depth = queue.popleft()
            if current_depth >= depth:
                continue

            adjacent: list[tuple[str, dict[str, Any]]] = []
            if direction in ("down", "both"):
                adjacent.extend((edge["target"], edge) for edge in self.out_edges.get(node_id, []))
            if direction in ("up", "both"):
                adjacent.extend((edge["source"], edge) for edge in self.in_edges.get(node_id, []))

            for neighbor, edge in adjacent:
                if relations and edge.get("relation") not in relations:
                    continue
                selected_edges[edge_key(edge)] = edge
                if neighbor not in seen:
                    if len(seen) >= max_nodes:
                        continue
                    seen.add(neighbor)
                    queue.append((neighbor, current_depth + 1))

        return seen, sorted(selected_edges.values(), key=edge_key)

    def edges_within(self, node_ids: set[str]) -> list[dict[str, Any]]:
        return sorted(
            [
                edge
                for edge in self.edges
                if edge["source"] in node_ids and edge["target"] in node_ids
            ],
            key=edge_key,
        )

    def open_candidates(self) -> list[dict[str, Any]]:
        candidates = []
        for node in self.nodes:
            status = str(node.get("status", ""))
            if (
                node.get("layer") == "open_gate"
                or node.get("future_paper_needed")
                or OPEN_RE.search(status)
            ):
                candidates.append(node)
        candidates.sort(key=lambda node: (node.get("layer", ""), node["id"]))
        return candidates

    def shortest_paths(
        self,
        start: str,
        goal: str,
        *,
        max_depth: int,
        direction: str,
        limit: int,
    ) -> list[list[tuple[str, dict[str, Any] | None, str | None]]]:
        self.require_node(start)
        self.require_node(goal)
        queue = deque([[(start, None, None)]])
        found: list[list[tuple[str, dict[str, Any] | None, str | None]]] = []

        while queue and len(found) < limit:
            path = queue.popleft()
            current = path[-1][0]
            if current == goal:
                found.append(path)
                continue
            if len(path) - 1 >= max_depth:
                continue

            used = {step[0] for step in path}
            adjacent: list[tuple[str, dict[str, Any], str]] = []
            if direction in ("out", "both"):
                adjacent.extend((edge["target"], edge, "out") for edge in self.out_edges.get(current, []))
            if direction in ("in", "both"):
                adjacent.extend((edge["source"], edge, "in") for edge in self.in_edges.get(current, []))

            adjacent.sort(key=lambda item: (item[0], item[1].get("relation", "")))
            for neighbor, edge, edge_direction in adjacent:
                if neighbor in used:
                    continue
                queue.append(path + [(neighbor, edge, edge_direction)])

        return found


def fail(message: str) -> None:
    print(f"FAIL: {message}", file=sys.stderr)
    raise SystemExit(1)


def node_summary(node: dict[str, Any]) -> str:
    parts = [node["id"], f"[{node.get('layer')}/{node.get('type')}/{node.get('status')}]"]
    if node.get("source_kind"):
        parts.append(f"source_kind={node['source_kind']}")
    if node.get("future_paper_needed"):
        parts.append("future_paper_needed=true")
    if node.get("tex_anchor"):
        anchor = node["tex_anchor"]
        label = anchor.get("nearest_label") or {}
        label_text = f"#{label.get('name')}" if label.get("name") else ""
        parts.append(f"tex={anchor.get('file')}:{anchor.get('line')}{label_text}")
    label = compact(node.get("label") or node.get("meaning") or "")
    if label and label != node["id"]:
        parts.append(f"- {label}")
    return " ".join(parts)


def edge_summary(edge: dict[str, Any]) -> str:
    note = compact(edge.get("note"), 90)
    suffix = f" - {note}" if note else ""
    return f"{edge['source']} -[{edge.get('relation', '?')}]-> {edge['target']}{suffix}"


def make_subgraph(
    index: GraphIndex,
    node_ids: set[str],
    edges: list[dict[str, Any]] | None,
    description: str,
) -> dict[str, Any]:
    selected_nodes = [node for node in index.nodes if node["id"] in node_ids]
    selected_edges = edges if edges is not None else index.edges_within(node_ids)
    return {
        "atlas_name": index.graph.get("atlas_name"),
        "release_version": index.graph.get("release_version"),
        "schema_version": index.graph.get("schema_version"),
        "subset_description": description,
        "node_count": len(selected_nodes),
        "edge_count": len(selected_edges),
        "nodes": selected_nodes,
        "edges": selected_edges,
    }


def maybe_write_subgraph(path: Path | None, subgraph: dict[str, Any] | None) -> None:
    if not path:
        return
    if not subgraph:
        fail("--write-subgraph was requested but no subgraph was produced")
    path.write_text(json.dumps(subgraph, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {path}", file=sys.stderr)


def print_nodes(title: str, nodes: list[dict[str, Any]], limit: int | None = None) -> None:
    if title:
        print(f"## {title}")
    shown = nodes[:limit] if limit else nodes
    if not shown:
        print("- none")
        return
    for node in shown:
        print(f"- {node_summary(node)}")
    if limit and len(nodes) > limit:
        print(f"- ... {len(nodes) - limit} more")


def print_edges(title: str, edges: list[dict[str, Any]], limit: int = 40) -> None:
    if title:
        print(f"## {title}")
    if not edges:
        print("- none")
        return
    for edge in edges[:limit]:
        print(f"- {edge_summary(edge)}")
    if len(edges) > limit:
        print(f"- ... {len(edges) - limit} more")


def output_json(payload: Any) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))


def add_common_filters(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--layer", help="Only match nodes in this layer")
    parser.add_argument("--type", dest="node_type", help="Only match nodes with this type")
    parser.add_argument("--status", help="Only match nodes whose status contains this text")
    parser.add_argument("--source-kind", help="Only match nodes with this source_kind")


def command_search(index: GraphIndex, args: argparse.Namespace) -> None:
    query = arg_text(args.query)
    matches = index.search(
        query,
        layer=args.layer,
        node_type=args.node_type,
        status=args.status,
        source_kind=args.source_kind,
        limit=args.limit,
    )
    node_ids = {node["id"] for _score, node in matches}
    subgraph = make_subgraph(index, node_ids, index.edges_within(node_ids), f"search: {query}")
    maybe_write_subgraph(args.write_subgraph, subgraph)

    if args.format == "json":
        output_json({"query": query, "matches": [{"score": score, "node": node} for score, node in matches]})
    elif args.format == "ids":
        for _score, node in matches:
            print(node["id"])
    else:
        print(f"# Search: {query}")
        if not matches:
            print("- none")
        for score, node in matches:
            print(f"- score={score} {node_summary(node)}")


def command_show(index: GraphIndex, args: argparse.Namespace) -> None:
    node = index.require_node(args.node_id)
    incoming = sorted(index.in_edges.get(node["id"], []), key=edge_key)
    outgoing = sorted(index.out_edges.get(node["id"], []), key=edge_key)
    if args.format == "json":
        output_json({"node": node, "incoming": incoming, "outgoing": outgoing})
        return
    print(f"# {node['id']}")
    for key in sorted(node):
        value = node[key]
        if isinstance(value, (list, dict)):
            print(f"- {key}: {json.dumps(value, ensure_ascii=False)}")
        else:
            print(f"- {key}: {value}")
    print_edges("Incoming", incoming, args.edge_limit)
    print_edges("Outgoing", outgoing, args.edge_limit)


def command_neighborhood(index: GraphIndex, args: argparse.Namespace) -> None:
    index.require_node(args.node_id)
    relations = set(args.relation or []) or None
    node_ids, edges = index.neighborhood(
        [args.node_id],
        depth=args.depth,
        direction=args.direction,
        relations=relations,
        max_nodes=args.max_nodes,
    )
    subgraph = make_subgraph(index, node_ids, edges, f"neighborhood: {args.node_id}")
    maybe_write_subgraph(args.write_subgraph, subgraph)

    if args.format == "json":
        output_json(subgraph)
        return
    print(f"# Neighborhood: {args.node_id}")
    print(f"- depth: {args.depth}")
    print(f"- direction: {args.direction}")
    print(f"- nodes: {len(subgraph['nodes'])}")
    print(f"- edges: {len(edges)}")
    print_nodes("Nodes", subgraph["nodes"], args.node_limit)
    print_edges("Edges", edges, args.edge_limit)


def command_paths(index: GraphIndex, args: argparse.Namespace) -> None:
    paths = index.shortest_paths(
        args.source,
        args.target,
        max_depth=args.max_depth,
        direction=args.direction,
        limit=args.limit,
    )
    if args.format == "json":
        payload = []
        for path in paths:
            payload.append(
                [
                    {"node": node_id, "edge": edge, "direction": edge_direction}
                    for node_id, edge, edge_direction in path
                ]
            )
        output_json({"source": args.source, "target": args.target, "paths": payload})
        return

    print(f"# Paths: {args.source} to {args.target}")
    if not paths:
        print("- none")
        return
    for index_path, path in enumerate(paths, 1):
        parts = [path[0][0]]
        for node_id, edge, edge_direction in path[1:]:
            relation = edge.get("relation", "?") if edge else "?"
            if edge_direction == "out":
                parts.append(f"-[{relation}]-> {node_id}")
            else:
                parts.append(f"<-[{relation}]- {node_id}")
        print(f"{index_path}. {' '.join(parts)}")


def command_source(index: GraphIndex, args: argparse.Namespace) -> None:
    query = arg_text(args.query)
    resolutions = index.source_resolution_matches(query)
    nodes = index.source_nodes(query)
    node_ids = {node["id"] for node in nodes}
    subgraph = make_subgraph(index, node_ids, index.edges_within(node_ids), f"source: {query}")
    maybe_write_subgraph(args.write_subgraph, subgraph)

    if args.format == "json":
        output_json({"query": query, "source_resolution": resolutions, "nodes": nodes})
        return
    print(f"# Source: {query}")
    print("## Source Resolution")
    if not resolutions:
        print("- none")
    for entry in resolutions:
        print(
            "- "
            f"{entry.get('legacy_ref')} -> {entry.get('canonical_ref')} "
            f"({entry.get('source_kind')})"
        )
    print_nodes("Nodes", nodes, args.node_limit)


def command_open_gates(index: GraphIndex, args: argparse.Namespace) -> None:
    candidates = index.open_candidates()
    query = arg_text(args.query)
    if query:
        seed_matches = index.search(query, limit=args.seed_limit)
        seeds = [node["id"] for _score, node in seed_matches]
        if seeds:
            context_nodes, _edges = index.neighborhood(
                seeds,
                depth=args.depth,
                direction="both",
                max_nodes=args.max_nodes,
            )
        else:
            context_nodes = set()
        terms = normalize_terms(query)
        filtered = []
        for node in candidates:
            text_match = all(term in index.node_text[node["id"]] for term in terms)
            if text_match or node["id"] in context_nodes:
                filtered.append(node)
        candidates = filtered

    node_ids = {node["id"] for node in candidates}
    subgraph = make_subgraph(index, node_ids, index.edges_within(node_ids), f"open-gates: {query or 'all'}")
    maybe_write_subgraph(args.write_subgraph, subgraph)

    if args.format == "json":
        output_json({"query": query, "nodes": candidates})
        return
    title = "Open Gates" if not query else f"Open Gates: {query}"
    print_nodes(title, candidates, args.limit)


def classify_context_nodes(nodes: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    buckets: dict[str, list[dict[str, Any]]] = {
        "claims": [],
        "equations": [],
        "open_gates": [],
        "future_paper_notes": [],
        "sources": [],
        "paper_backlinks": [],
        "other": [],
    }
    for node in nodes:
        layer = node.get("layer")
        node_type = node.get("type")
        if node.get("future_paper_needed") or node.get("source_kind") == "future_paper_note":
            buckets["future_paper_notes"].append(node)
        elif layer == "open_gate" or OPEN_RE.search(str(node.get("status", ""))):
            buckets["open_gates"].append(node)
        elif layer == "claim_theorem" or node_type == "claim":
            buckets["claims"].append(node)
        elif layer == "equation_anchor" or node_type == "equation":
            buckets["equations"].append(node)
        elif layer in ("file_anchor", "source_section_anchor") or node_type in ("source_file", "section_anchor"):
            buckets["sources"].append(node)
        elif layer == "paper_backlink":
            buckets["paper_backlinks"].append(node)
        else:
            buckets["other"].append(node)
    for bucket in buckets.values():
        bucket.sort(key=lambda node: (node.get("layer", ""), node["id"]))
    return buckets


def command_context(index: GraphIndex, args: argparse.Namespace) -> None:
    topic = arg_text(args.topic)
    seed_matches = index.search(topic, limit=args.seed_limit)
    seeds = [node["id"] for _score, node in seed_matches]
    if seeds:
        node_ids, edges = index.neighborhood(
            seeds,
            depth=args.depth,
            direction="both",
            max_nodes=args.max_nodes,
        )
    else:
        node_ids, edges = set(), []

    subgraph = make_subgraph(index, node_ids, edges, f"context: {topic}")
    maybe_write_subgraph(args.write_subgraph, subgraph)

    if args.format == "json":
        output_json({"topic": topic, "seeds": seed_matches, "subgraph": subgraph})
        return

    buckets = classify_context_nodes(subgraph["nodes"])
    print(f"# Context: {topic}")
    print(f"- seed matches: {len(seed_matches)}")
    print(f"- packet nodes: {len(subgraph['nodes'])}")
    print(f"- packet edges: {len(edges)}")
    print_nodes("Seed Matches", [node for _score, node in seed_matches], args.seed_limit)
    print_nodes("Claims", buckets["claims"], args.section_limit)
    print_nodes("Equations", buckets["equations"], args.section_limit)
    print_nodes("Open Gates", buckets["open_gates"], args.section_limit)
    print_nodes("Future Paper Notes", buckets["future_paper_notes"], args.section_limit)
    print_nodes("Sources", buckets["sources"], args.section_limit)
    print_nodes("Paper Backlinks", buckets["paper_backlinks"], args.section_limit)
    print_edges("Key Edges", edges, args.edge_limit)


def command_stats(index: GraphIndex, args: argparse.Namespace) -> None:
    layers = defaultdict(int)
    statuses = defaultdict(int)
    relations = defaultdict(int)
    source_kinds = defaultdict(int)
    for node in index.nodes:
        layers[node.get("layer", "unknown")] += 1
        statuses[node.get("status", "unknown")] += 1
        source_kinds[node.get("source_kind", "none")] += 1
    for edge in index.edges:
        relations[edge.get("relation", "unknown")] += 1

    payload = {
        "nodes": len(index.nodes),
        "edges": len(index.edges),
        "layers": dict(sorted(layers.items())),
        "source_kinds": dict(sorted(source_kinds.items())),
        "top_statuses": dict(sorted(statuses.items(), key=lambda item: (-item[1], item[0]))[: args.limit]),
        "top_relations": dict(sorted(relations.items(), key=lambda item: (-item[1], item[0]))[: args.limit]),
    }
    if args.format == "json":
        output_json(payload)
        return
    print("# Graph Stats")
    print(f"- nodes: {payload['nodes']}")
    print(f"- edges: {payload['edges']}")
    print("## Layers")
    for key, count in payload["layers"].items():
        print(f"- {key}: {count}")
    print("## Source Kinds")
    for key, count in payload["source_kinds"].items():
        print(f"- {key}: {count}")
    print("## Top Statuses")
    for key, count in payload["top_statuses"].items():
        print(f"- {key}: {count}")
    print("## Top Relations")
    for key, count in payload["top_relations"].items():
        print(f"- {key}: {count}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=DEFAULT_GRAPH, help="Input graph JSON")
    subparsers = parser.add_subparsers(dest="command", required=True)

    search = subparsers.add_parser("search", help="Find graph nodes by text")
    search.add_argument("query", nargs="+")
    search.add_argument("--limit", type=int, default=20)
    search.add_argument("--format", choices=("md", "json", "ids"), default="md")
    search.add_argument("--write-subgraph", type=Path)
    add_common_filters(search)
    search.set_defaults(func=command_search)

    show = subparsers.add_parser("show", help="Show one node and its direct edges")
    show.add_argument("node_id")
    show.add_argument("--format", choices=("md", "json"), default="md")
    show.add_argument("--edge-limit", type=int, default=40)
    show.set_defaults(func=command_show)

    neighborhood = subparsers.add_parser("neighborhood", help="Extract a bounded neighborhood")
    neighborhood.add_argument("node_id")
    neighborhood.add_argument("--depth", type=int, default=1)
    neighborhood.add_argument("--direction", choices=("up", "down", "both"), default="both")
    neighborhood.add_argument("--relation", action="append", help="Restrict to this relation; repeatable")
    neighborhood.add_argument("--max-nodes", type=int, default=200)
    neighborhood.add_argument("--node-limit", type=int, default=80)
    neighborhood.add_argument("--edge-limit", type=int, default=80)
    neighborhood.add_argument("--format", choices=("md", "json"), default="md")
    neighborhood.add_argument("--write-subgraph", type=Path)
    neighborhood.set_defaults(func=command_neighborhood)

    paths = subparsers.add_parser("paths", help="Find short graph paths between two nodes")
    paths.add_argument("source")
    paths.add_argument("target")
    paths.add_argument("--max-depth", type=int, default=4)
    paths.add_argument("--direction", choices=("out", "in", "both"), default="out")
    paths.add_argument("--limit", type=int, default=5)
    paths.add_argument("--format", choices=("md", "json"), default="md")
    paths.set_defaults(func=command_paths)

    source = subparsers.add_parser("source", help="Find nodes tied to a source path or legacy source name")
    source.add_argument("query", nargs="+")
    source.add_argument("--node-limit", type=int, default=80)
    source.add_argument("--format", choices=("md", "json"), default="md")
    source.add_argument("--write-subgraph", type=Path)
    source.set_defaults(func=command_source)

    open_gates = subparsers.add_parser("open-gates", help="List open gates and future-paper markers")
    open_gates.add_argument("query", nargs="*")
    open_gates.add_argument("--depth", type=int, default=1)
    open_gates.add_argument("--seed-limit", type=int, default=16)
    open_gates.add_argument("--max-nodes", type=int, default=220)
    open_gates.add_argument("--limit", type=int, default=120)
    open_gates.add_argument("--format", choices=("md", "json"), default="md")
    open_gates.add_argument("--write-subgraph", type=Path)
    open_gates.set_defaults(func=command_open_gates)

    context = subparsers.add_parser("context", help="Build a compact topic packet")
    context.add_argument("topic", nargs="+")
    context.add_argument("--depth", type=int, default=1)
    context.add_argument("--seed-limit", type=int, default=12)
    context.add_argument("--max-nodes", type=int, default=180)
    context.add_argument("--section-limit", type=int, default=30)
    context.add_argument("--edge-limit", type=int, default=80)
    context.add_argument("--format", choices=("md", "json"), default="md")
    context.add_argument("--write-subgraph", type=Path)
    context.set_defaults(func=command_context)

    stats = subparsers.add_parser("stats", help="Summarize graph counts")
    stats.add_argument("--limit", type=int, default=16)
    stats.add_argument("--format", choices=("md", "json"), default="md")
    stats.set_defaults(func=command_stats)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    index = GraphIndex.load(args.graph)
    args.func(index, args)


if __name__ == "__main__":
    main()
