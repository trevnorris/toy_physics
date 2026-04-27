#!/usr/bin/env python3
"""Validate the generated human-readable reader site."""

from __future__ import annotations

import argparse
import html
import re
import sys
from pathlib import Path
from typing import Any

import yaml

import generate_reader_site as reader


ROOT = Path(__file__).resolve().parents[2]
ATLAS_DIR = ROOT / "atlas"
GRAPH_YAML = ROOT / "graph" / "fluid_universe_derivation_atlas_graph.yaml"
TOPICS_YAML = ATLAS_DIR / "topics.yaml"
README = ROOT / "README.md"
SITE_DIR = ATLAS_DIR / "site"
EXPORT_DIR = ATLAS_DIR / "exports"
HREF_RE = re.compile(r'href="([^"]+)"')
MATH_BLOCK_RE = re.compile(r"\\\[(.*?)\\\]", re.S)
BARE_GROUP_RE = re.compile(r"[_^](?![{\\])(?:[A-Za-z0-9]{2,}|\([^)]+\))")
BAD_MATH_SHORTHAND_RE = re.compile(r"\bmhat(?:[_A-Za-z0-9]|\b)")
UNICODE_SUBSUP_RE = re.compile(r"[⁰¹²³⁴⁵⁶⁷⁸⁹₀₁₂₃₄₅₆₇₈₉ᵀ]")
BARE_MATH_TEXT_RE = re.compile(r"(?<![{A-Za-z\\])(?:equivalently|rank)(?![A-Za-z])")


def load_yaml(path: Path) -> Any:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def validate_topics(graph: dict[str, Any], topics: dict[str, Any], errors: list[str], warnings: list[str]) -> None:
    node_ids = {node["id"] for node in graph["nodes"]}
    slugs = {topic["slug"] for topic in topics["topics"]}
    for topic in topics["topics"]:
        if not topic.get("summary"):
            errors.append(f"{topic['slug']} missing summary")
        sections = topic.get("sections") or {}
        for key in ("physical_picture", "mathematical_core", "status_note"):
            if not sections.get(key):
                warnings.append(f"{topic['slug']} missing {key} prose")
            for index, item in enumerate(sections.get(key, []) or []):
                if not isinstance(item, str):
                    errors.append(f"{topic['slug']} {key}[{index}] must be a string, got {type(item).__name__}")
        for node_id in topic.get("node_ids", []) or []:
            if node_id not in node_ids:
                errors.append(f"{topic['slug']} references missing node {node_id}")
        for slug in topic.get("related_topics", []) or []:
            if slug not in slugs:
                errors.append(f"{topic['slug']} references missing related topic {slug}")


def validate_site_files(site_dir: Path, topics: dict[str, Any], errors: list[str]) -> None:
    required = [
        site_dir / "index.html",
        site_dir / "future-paper-worklist.html",
        site_dir / "assets" / "reader.css",
        site_dir / "assets" / "reader.js",
        site_dir / "data" / "topics.json",
        site_dir / "data" / "future_paper_backfill.yaml",
        site_dir / reader.MANIFEST,
    ]
    required.extend(site_dir / reader.slug_to_file(topic["slug"]) for topic in topics["topics"])
    for path in required:
        if not path.exists():
            errors.append(f"missing generated site file: {path.relative_to(site_dir)}")


def validate_links(site_dir: Path, errors: list[str]) -> None:
    for path in sorted(site_dir.rglob("*.html")):
        text = path.read_text(encoding="utf-8")
        for href in HREF_RE.findall(text):
            if href.startswith(("http://", "https://", "mailto:", "#")):
                continue
            clean = href.split("#", 1)[0]
            if not clean:
                continue
            target = (path.parent / clean).resolve()
            try:
                target.relative_to(site_dir.resolve())
            except ValueError:
                errors.append(f"{path.relative_to(site_dir)} links outside site: {href}")
                continue
            if not target.exists():
                errors.append(f"{path.relative_to(site_dir)} broken local link: {href}")


def validate_math_blocks(site_dir: Path, errors: list[str]) -> None:
    for path in sorted(site_dir.rglob("*.html")):
        text = path.read_text(encoding="utf-8")
        for raw_block in MATH_BLOCK_RE.findall(text):
            block = html.unescape(raw_block)
            shorthand = BAD_MATH_SHORTHAND_RE.search(block)
            if shorthand:
                errors.append(
                    f"{path.relative_to(site_dir)} has unexpanded graph math shorthand "
                    f"`{shorthand.group(0)}` in `{block.strip()}`"
                )
            unicode_subsup = UNICODE_SUBSUP_RE.search(block)
            if unicode_subsup:
                errors.append(
                    f"{path.relative_to(site_dir)} has Unicode sub/superscript shorthand "
                    f"`{unicode_subsup.group(0)}` in `{block.strip()}`"
                )
            bare_text = BARE_MATH_TEXT_RE.search(block)
            if bare_text:
                errors.append(
                    f"{path.relative_to(site_dir)} has bare text/operator word "
                    f"`{bare_text.group(0)}` in `{block.strip()}`"
                )
            if len(reader.split_top_level_semicolons(block)) > 1:
                errors.append(
                    f"{path.relative_to(site_dir)} has a top-level semicolon in a MathJax block: "
                    f"`{block.strip()}`"
                )
            match = BARE_GROUP_RE.search(block)
            if match:
                errors.append(
                    f"{path.relative_to(site_dir)} has unbraced multi-character TeX group "
                    f"`{match.group(0)}` in `{block.strip()}`"
                )


def validate_backfill(site_dir: Path, export_dir: Path, errors: list[str]) -> int:
    site_backfill = site_dir / "data" / "future_paper_backfill.yaml"
    export_backfill = export_dir / "future_paper_backfill.yaml"
    if not site_backfill.exists() or not export_backfill.exists():
        errors.append("future paper backfill YAML is missing")
        return 0
    site_data = load_yaml(site_backfill) or []
    export_data = load_yaml(export_backfill) or []
    if site_data != export_data:
        errors.append("site and export future_paper_backfill.yaml differ")
    for index, item in enumerate(site_data):
        for key in ("topic_slug", "source", "node_ids", "tag"):
            if key not in item:
                errors.append(f"backfill entry {index} missing {key}")
        if item.get("tag") != "future_paper_needed":
            errors.append(f"backfill entry {index} has unexpected tag {item.get('tag')}")
    return len(site_data)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=GRAPH_YAML)
    parser.add_argument("--topics", type=Path, default=TOPICS_YAML)
    parser.add_argument("--readme", type=Path, default=README)
    parser.add_argument("--site", type=Path, default=SITE_DIR)
    parser.add_argument("--exports", type=Path, default=EXPORT_DIR)
    args = parser.parse_args()

    errors: list[str] = []
    warnings: list[str] = []
    graph = load_yaml(args.graph)
    topics = load_yaml(args.topics)
    doi_entries = reader.parse_readme_dois(args.readme)
    if len(doi_entries) < 10:
        errors.append(f"expected README DOI metadata; found only {len(doi_entries)} DOI entries")

    validate_topics(graph, topics, errors, warnings)
    validate_site_files(args.site, topics, errors)
    validate_links(args.site, errors)
    validate_math_blocks(args.site, errors)
    backfill_count = validate_backfill(args.site, args.exports, errors)

    if errors:
        print("Reader site validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        if warnings:
            print("Warnings:", file=sys.stderr)
            for warning in warnings:
                print(f"- {warning}", file=sys.stderr)
        sys.exit(1)

    print(f"OK: reader site validation passed ({len(topics['topics'])} topics, {backfill_count} backfill entries)")
    if warnings:
        print(f"warnings: {len(warnings)}")


if __name__ == "__main__":
    main()
