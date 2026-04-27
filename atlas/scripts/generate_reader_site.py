#!/usr/bin/env python3
"""Generate a human-readable static reader site from the atlas graph."""

from __future__ import annotations

import argparse
import html
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml


ROOT = Path(__file__).resolve().parents[2]
ATLAS_DIR = ROOT / "atlas"
GRAPH_YAML = ROOT / "graph" / "fluid_universe_derivation_atlas_graph.yaml"
TOPICS_YAML = ATLAS_DIR / "topics.yaml"
README = ROOT / "README.md"
SITE_DIR = ATLAS_DIR / "site"
EXPORT_DIR = ATLAS_DIR / "exports"
MANIFEST = "generated_manifest.txt"

GENERATED_NOTICE = "Generated from atlas/topics.yaml and graph/fluid_universe_derivation_atlas_graph.yaml."
ATLAS_ID_RE = re.compile(r"^[A-Z][A-Z0-9]*(?:_[A-Z0-9]+)+$")
MATH_UNICODE_TRANSLATION = str.maketrans(
    {
        "⁰": "^{0}",
        "¹": "^{1}",
        "²": "^{2}",
        "³": "^{3}",
        "⁴": "^{4}",
        "⁵": "^{5}",
        "⁶": "^{6}",
        "⁷": "^{7}",
        "⁸": "^{8}",
        "⁹": "^{9}",
        "₀": "_{0}",
        "₁": "_{1}",
        "₂": "_{2}",
        "₃": "_{3}",
        "₄": "_{4}",
        "₅": "_{5}",
        "₆": "_{6}",
        "₇": "_{7}",
        "₈": "_{8}",
        "₉": "_{9}",
        "ᵀ": "^{T}",
    }
)

SOURCE_KEYS = (
    "file",
    "source_note_file",
    "legacy_file",
)
SOURCE_LIST_KEYS = (
    "files",
    "sources",
    "source_files",
    "source_note_files",
    "canonical_target_files",
    "legacy_sources",
    "legacy_files",
)

FOLDER_DOI_MATCHERS = {
    "pde_ledger": ("moving-throat pde derivation companion",),
    "4d": ("4d toy model", "action", "projections"),
    "4d_1pn_bridge": ("deriving key post-newtonian coefficients",),
    "4d_em_fields": ("maxwell from the unified 4d toy model",),
    "4d_plasma": ("plasma dynamics from the unified 4d",),
    "4d_1pn_full": ("full first post-newtonian sector",),
    "4d_2pn": ("full conservative second post-newtonian",),
    "4d_2_5pn": ("point-particle 2.5pn",),
    "4d_3pn": ("two-body 3pn",),
    "4d_4pn": ("two-body 4pn",),
    "1pn_orbital_dynamics": ("newtonian and 1pn orbital dynamics",),
    "1pn_optics": ("gravitational optics",),
    "1pn_spin_and_nbody": ("spin, vorticity, and n-body",),
    "em_fields": ("electromagnetic fields and charged defects",),
    "brane_bulk_ontology": ("brane-bulk throat ontology",),
    "1pn_hybrid": ("hybrid 1pn dynamics",),
}

STATUS_EXPLAINERS = {
    "published_backed": "paper-backed",
    "mixed_published_and_open": "paper-backed with open gates",
    "notes_backed_pending_papers": "paper pending",
    "graph_backed": "graph guardrail",
    "generated_backfill": "generated worklist",
}


def load_yaml(path: Path) -> Any:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def generated_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def clean_previous_site(site_dir: Path) -> None:
    manifest = site_dir / MANIFEST
    if not manifest.exists():
        return
    root = site_dir.resolve()
    for raw in manifest.read_text(encoding="utf-8").splitlines():
        rel = raw.strip()
        if not rel or rel.startswith("#"):
            continue
        target = (site_dir / rel).resolve()
        try:
            target.relative_to(root)
        except ValueError:
            continue
        if target.is_file():
            target.unlink()


def write_file(path: Path, content: str, generated_files: set[Path]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    generated_files.add(path)


def parse_readme_dois(path: Path) -> list[dict[str, str]]:
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(r"^- \*(?P<title>.+?)\.\* \[(?P<doi>10\.5281/zenodo\.\d+)\]\((?P<url>https://doi\.org/[^)]+)\)", re.M)
    return [match.groupdict() for match in pattern.finditer(text)]


def doi_for_folder(folder: str, doi_entries: list[dict[str, str]]) -> dict[str, str] | None:
    needles = FOLDER_DOI_MATCHERS.get(folder)
    if not needles:
        return None
    for entry in doi_entries:
        title = entry["title"].lower()
        if all(needle in title for needle in needles):
            return entry
    return None


def node_title(node: dict[str, Any]) -> str:
    return str(node.get("label") or node["id"])


def compact(text: Any, limit: int = 180) -> str:
    value = " ".join(str(text or "").split())
    return value if len(value) <= limit else value[: limit - 3] + "..."


def escape(text: Any) -> str:
    return html.escape(str(text or ""), quote=True)


def looks_like_atlas_id(value: str) -> bool:
    return "/" not in value and "." not in value and bool(ATLAS_ID_RE.fullmatch(value))


def expand_source_value(value: str) -> list[str]:
    if looks_like_atlas_id(value):
        return []
    if " / " in value:
        return [part.strip() for part in value.split(" / ") if part.strip()]
    return [value]


def slug_to_file(slug: str) -> str:
    return f"topics/{slug}.html"


def topic_href(slug: str, from_topic: bool = False) -> str:
    prefix = "../" if from_topic else ""
    return f"{prefix}{slug_to_file(slug)}"


def css_href(from_topic: bool = False) -> str:
    return "../assets/reader.css" if from_topic else "assets/reader.css"


def js_href(from_topic: bool = False) -> str:
    return "../assets/reader.js" if from_topic else "assets/reader.js"


def normalize_source(raw: str, resolution: dict[str, dict[str, Any]]) -> tuple[str, dict[str, Any] | None]:
    if raw in resolution:
        entry = resolution[raw]
        return str(entry.get("canonical_ref") or raw), entry
    if ".md:" in raw:
        base, suffix = raw.split(".md:", 1)
        legacy = base + ".md"
        if legacy in resolution:
            entry = resolution[legacy]
            canonical = str(entry.get("canonical_ref") or legacy)
            return canonical + ":" + suffix, entry
    return raw, None


def iter_node_sources(node: dict[str, Any], resolution: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    found: list[dict[str, Any]] = []
    for key in SOURCE_KEYS:
        value = node.get(key)
        if isinstance(value, str):
            for expanded in expand_source_value(value):
                path, entry = normalize_source(expanded, resolution)
                found.append({"path": path, "source_key": key, "resolution": entry})
    for key in SOURCE_LIST_KEYS:
        for value in node.get(key, []) or []:
            if isinstance(value, str):
                for expanded in expand_source_value(value):
                    path, entry = normalize_source(expanded, resolution)
                    found.append({"path": path, "source_key": key, "resolution": entry})
    anchor = node.get("tex_anchor") or {}
    if isinstance(anchor.get("file"), str):
        found.append({"path": anchor["file"], "source_key": "tex_anchor", "resolution": None})
    return found


def classify_source(source: dict[str, Any], node: dict[str, Any], doi_entries: list[dict[str, str]]) -> dict[str, Any]:
    path = source["path"]
    resolution = source.get("resolution") or {}
    source_kind = resolution.get("source_kind") or node.get("source_kind") or ""
    future_needed = bool(resolution.get("future_paper_needed") or node.get("future_paper_needed"))
    info = {
        "path": path,
        "source_kind": source_kind or "unknown",
        "status": "uncategorized",
        "title": "",
        "doi": "",
        "url": "",
        "reason": "",
    }
    if future_needed or source_kind == "future_paper_note" or path.startswith("notes/"):
        info["status"] = "paper_pending"
        info["reason"] = resolution.get("note") or "No formal paper source is recorded yet."
        return info
    if path.startswith("research/") and "/notes/" in path:
        parts = path.split("/")
        folder = parts[1] if len(parts) > 1 else ""
        entry = doi_for_folder(folder, doi_entries)
        info.update(
            {
                "status": "support_note",
                "title": entry["title"] if entry else folder,
                "doi": entry["doi"] if entry else "",
                "url": entry["url"] if entry else "",
                "reason": "Supporting research note associated with a paper; cite the paper DOI in public prose.",
            }
        )
        return info
    if path.startswith("research/"):
        parts = path.split("/")
        folder = parts[1] if len(parts) > 1 else ""
        entry = doi_for_folder(folder, doi_entries)
        if entry:
            info.update(
                {
                    "status": "published",
                    "title": entry["title"],
                    "doi": entry["doi"],
                    "url": entry["url"],
                    "source_kind": source_kind or "paper",
                }
            )
        else:
            info.update(
                {
                    "status": "research_no_doi",
                    "title": folder,
                    "reason": "Research source exists but no DOI was matched from README.md.",
                }
            )
        return info
    if path == "derived_from_current_task":
        info["status"] = "synthetic"
        info["reason"] = "Synthetic graph provenance."
        return info
    info["status"] = "paper_pending" if path.endswith(".md") else "uncategorized"
    info["reason"] = "Markdown or legacy source should be replaced by a formal paper citation when available."
    return info


def unique_sources(sources: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen = set()
    result = []
    for source in sources:
        key = (source["path"], source.get("status"), source.get("doi"))
        if key in seen:
            continue
        seen.add(key)
        result.append(source)
    return sorted(result, key=lambda item: (item["status"], item["path"]))


def build_indexes(graph: dict[str, Any]) -> tuple[dict[str, dict[str, Any]], dict[str, list[dict[str, Any]]], dict[str, list[dict[str, Any]]]]:
    node_by_id = {node["id"]: node for node in graph["nodes"]}
    incoming: dict[str, list[dict[str, Any]]] = defaultdict(list)
    outgoing: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for edge in graph["edges"]:
        outgoing[edge["source"]].append(edge)
        incoming[edge["target"]].append(edge)
    return node_by_id, incoming, outgoing


def topic_nodes(topic: dict[str, Any], node_by_id: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    nodes = []
    for node_id in topic.get("node_ids", []) or []:
        node = node_by_id.get(node_id)
        if node:
            nodes.append(node)
    return nodes


def source_anchor_nodes(
    nodes: list[dict[str, Any]],
    node_by_id: dict[str, dict[str, Any]],
    incoming: dict[str, list[dict[str, Any]]],
    outgoing: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    selected = {node["id"] for node in nodes}
    source_ids = set()
    for node in nodes:
        for edge in incoming[node["id"]] + outgoing[node["id"]]:
            other_id = edge["source"] if edge["source"] != node["id"] else edge["target"]
            other = node_by_id.get(other_id)
            if not other:
                continue
            if other.get("layer") in {"file_anchor", "source_section_anchor", "paper_backlink"}:
                source_ids.add(other_id)
        if node.get("layer") in {"file_anchor", "source_section_anchor", "paper_backlink"}:
            source_ids.add(node["id"])
    return [node_by_id[node_id] for node_id in sorted(source_ids - selected) if node_id in node_by_id]


def collect_sources_for_topic(
    topic: dict[str, Any],
    node_by_id: dict[str, dict[str, Any]],
    incoming: dict[str, list[dict[str, Any]]],
    outgoing: dict[str, list[dict[str, Any]]],
    resolution: dict[str, dict[str, Any]],
    doi_entries: list[dict[str, str]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    nodes = topic_nodes(topic, node_by_id)
    source_nodes = source_anchor_nodes(nodes, node_by_id, incoming, outgoing)
    sources = []
    for node in nodes + source_nodes:
        for source in iter_node_sources(node, resolution):
            classified = classify_source(source, node, doi_entries)
            classified["node_id"] = node["id"]
            classified["node_title"] = node_title(node)
            sources.append(classified)
    return unique_sources(sources), source_nodes


def node_group(nodes: list[dict[str, Any]], predicate) -> list[dict[str, Any]]:
    return [node for node in nodes if predicate(node)]


def is_equation(node: dict[str, Any]) -> bool:
    return node.get("layer") == "equation_anchor" or node.get("type") == "equation"


def is_claim(node: dict[str, Any]) -> bool:
    return node.get("layer") == "claim_theorem" or node.get("type") == "claim"


def is_open(node: dict[str, Any]) -> bool:
    return node.get("layer") == "open_gate" or "open" in str(node.get("status", "")).lower()


def is_firewall(node: dict[str, Any]) -> bool:
    return node["id"].startswith("FIREWALL_") or "firewall" in str(node.get("type", "")).lower()


def expression_for(node: dict[str, Any]) -> str:
    for key in ("expression", "math", "equation"):
        value = node.get(key)
        if isinstance(value, str) and value.strip():
            return normalize_reader_math(value.strip())
        if isinstance(value, list) and value:
            expression = "; ".join(normalize_tex_expression(str(item)) for item in value)
            return format_display_expression(expression)
    return ""


def normalize_reader_math(expression: str) -> str:
    return format_display_expression(normalize_tex_expression(expression))


def format_display_expression(expression: str) -> str:
    clauses = split_top_level_semicolons(expression)
    if len(clauses) < 2 or not all("=" in clause for clause in clauses):
        return expression

    rows = []
    for clause in clauses:
        lhs, rhs = clause.split("=", 1)
        rows.append(f"{lhs.strip()} &= {rhs.strip()}")
    return r"\begin{aligned}" + r" \\ ".join(rows) + r"\end{aligned}"


def split_top_level_semicolons(expression: str) -> list[str]:
    parts: list[str] = []
    start = 0
    depth = 0
    for index, char in enumerate(expression):
        if char in "([{":
            depth += 1
        elif char in ")]}":
            depth = max(0, depth - 1)
        elif char == ";" and depth == 0 and (index == 0 or expression[index - 1] != "\\"):
            parts.append(expression[start:index].strip())
            start = index + 1
    parts.append(expression[start:].strip())
    return [part for part in parts if part]


def normalize_tex_expression(expression: str) -> str:
    """Normalize compact graph math shorthand for correct MathJax output."""
    expression = expand_math_shorthand(expression)
    out: list[str] = []
    index = 0
    while index < len(expression):
        char = expression[index]
        if char not in "_^":
            out.append(char)
            index += 1
            continue

        out.append(char)
        index += 1
        if index >= len(expression):
            continue

        next_char = expression[index]
        if next_char in "{\\":
            continue
        if next_char.isspace():
            continue
        if next_char == "(":
            end = expression.find(")", index + 1)
            if end != -1:
                out.append("{" + expression[index : end + 1] + "}")
                index = end + 1
                continue

        start = index
        while index < len(expression):
            token_char = expression[index]
            if index > start and token_char.isupper() and expression[index - 1].islower():
                break
            if index > start and expression[start].isdigit() and token_char.isalpha() and not token_char.isupper():
                break
            if token_char.isalnum() or token_char in "'′":
                index += 1
                continue
            break
        token = expression[start:index]
        if not token:
            out.append(expression[index])
            index += 1
        elif len(token) > 1:
            out.append(tex_group(token))
        else:
            out.append(token)
    return "".join(out)


def expand_math_shorthand(expression: str) -> str:
    text = expression.translate(MATH_UNICODE_TRANSLATION)
    text = text.replace("μ0_eff", r"\mu_{0,\mathrm{eff}}")
    text = text.replace("μ0", r"\mu_0")
    text = text.replace("μ", r"\mu")
    text = text.replace("ν", r"\nu")
    text = text.replace("η", r"\eta")
    text = text.replace("Ω", r"\Omega")
    text = text.replace("Ξ", r"\Xi")
    text = text.replace("Λ", r"\Lambda")
    text = text.replace("α", r"\alpha")
    text = text.replace("ε", r"\epsilon")
    text = text.replace("λ", r"\lambda")
    text = text.replace("∈", r"\in")
    text = text.replace("≈", r"\approx")
    text = text.replace("∂", r"\partial")
    text = text.replace("Ŷ", r"\hat{Y}")
    text = text.replace("δx", r"\delta x")
    text = text.replace("δln", r"\delta\ln")
    text = re.sub(r"(?<=\\in )Z\b", r"\\mathbb{Z}", text)
    text = re.sub(r"\bq_lm,([A-Za-z]+)", r"q_{\\ell m,\1}", text)
    text = re.sub(r"\b([qS])_lm\b", r"\1_{\\ell m}", text)
    text = text.replace("l(l+1)", r"\ell(\ell+1)")
    text = re.sub(r"sqrt\(([^()]+)\)", r"\\sqrt{\1}", text)
    text = re.sub(r"(?<=_)psi\b", r"\\psi", text)
    text = re.sub(r"(?<![A-Za-z_\\])psi(?![A-Za-z])", r"\\psi", text)
    text = re.sub(r"(?<![A-Za-z_\\])Sigma(?![A-Za-z])", r"\\Sigma", text)
    text = re.sub(r"\s+~\s+", lambda _: r" \sim ", text)
    text = re.sub(r",\s*equivalently\s+", lambda _: r",\;\mathrm{equivalently}\; ", text)
    text = re.sub(r"(?<!\\)\brank(?=\s*\()", r"\\operatorname{rank}", text)
    text = re.sub(r"\bmhat([A-Za-z0-9]+)\b", r"\\hat{m}_\1", text)
    text = re.sub(r"\bmhat(?=_|\b)", r"\\hat{m}", text)
    text = re.sub(r"\^([0-9]+)", r"^{\1}", text)
    text = re.sub(r"\bXi([0-9]+)\b", r"\\Xi_\1", text)
    text = re.sub(r"(?<!\\)\bDelta\b", r"\\Delta", text)
    text = re.sub(r"(?<!\\)\bvarpi\b", r"\\varpi", text)
    text = re.sub(r"\b([A-Z])([0-9]+)\b", r"\1_\2", text)
    return text


def tex_group(token: str) -> str:
    if any(char.isalpha() for char in token):
        return "{\\mathrm{" + token + "}}"
    return "{" + token + "}"


def render_paragraphs(items: list[str]) -> str:
    return "\n".join(f"<p>{escape(item)}</p>" for item in items)


def render_node_list(nodes: list[dict[str, Any]], empty: str = "No graph nodes selected for this section.") -> str:
    if not nodes:
        return f"<p class=\"muted\">{escape(empty)}</p>"
    lines = ["<ul class=\"node-list\">"]
    for node in nodes:
        lines.append(
            "<li>"
            f"<code>{escape(node['id'])}</code> "
            f"<strong>{escape(node_title(node))}</strong> "
            f"<span class=\"badge small\">{escape(node.get('status', ''))}</span>"
            f"<p>{escape(compact(node.get('meaning') or node_title(node), 220))}</p>"
            "</li>"
        )
    lines.append("</ul>")
    return "\n".join(lines)


def render_equations(nodes: list[dict[str, Any]]) -> str:
    equations = [node for node in nodes if is_equation(node)]
    if not equations:
        return "<p class=\"muted\">No equation anchors selected for this topic yet.</p>"
    cards = []
    for node in equations:
        expr = expression_for(node)
        math = f"<div class=\"math-block\">\\[{escape(expr)}\\]</div>" if expr else ""
        cards.append(
            "<article class=\"evidence-card\">"
            f"<h3>{escape(node_title(node))}</h3>"
            f"<p class=\"node-id\">{escape(node['id'])} · {escape(node.get('status', ''))}</p>"
            f"{math}"
            f"<p>{escape(compact(node.get('meaning') or '', 260))}</p>"
            "</article>"
        )
    return "\n".join(cards)


def render_sources(sources: list[dict[str, Any]]) -> str:
    if not sources:
        return "<p class=\"muted\">No source metadata was found for this topic.</p>"
    published = [source for source in sources if source["status"] == "published"]
    pending = [source for source in sources if source["status"] == "paper_pending"]
    support = [source for source in sources if source["status"] == "support_note"]
    other = [source for source in sources if source["status"] not in {"published", "paper_pending", "support_note"}]
    chunks = []
    if published:
        chunks.append("<h3>Published sources</h3><ul class=\"source-list\">")
        for source in published:
            chunks.append(
                "<li>"
                f"<strong>{escape(source['title'])}</strong> "
                f"<a href=\"{escape(source['url'])}\">{escape(source['doi'])}</a>"
                f"<br><code>{escape(source['path'])}</code>"
                "</li>"
            )
        chunks.append("</ul>")
    if pending:
        chunks.append("<h3>Paper pending sources</h3><ul class=\"source-list pending\">")
        for source in pending:
            chunks.append(
                "<li>"
                "<span class=\"badge pending-badge\">paper pending</span> "
                f"<code>{escape(source['path'])}</code>"
                f"<p>{escape(source.get('reason') or 'A formal paper citation should be added when available.')}</p>"
                "</li>"
            )
        chunks.append("</ul>")
    if support:
        chunks.append("<h3>Supporting notes</h3><ul class=\"source-list\">")
        for source in support:
            paper = (
                f" Associated paper: <a href=\"{escape(source['url'])}\">{escape(source['doi'])}</a>."
                if source.get("url")
                else ""
            )
            chunks.append(
                "<li>"
                f"<span class=\"badge small\">support note</span> <code>{escape(source['path'])}</code>"
                f"<p>{escape(source.get('reason') or '')}{paper}</p>"
                "</li>"
            )
        chunks.append("</ul>")
    if other:
        chunks.append("<h3>Other source metadata</h3><ul class=\"source-list\">")
        for source in other:
            chunks.append(
                "<li>"
                f"<span class=\"badge small\">{escape(source['status'])}</span> "
                f"<code>{escape(source['path'])}</code>"
                f"<p>{escape(source.get('reason') or source.get('source_kind') or '')}</p>"
                "</li>"
            )
        chunks.append("</ul>")
    return "\n".join(chunks)


def page_shell(title: str, body: str, from_topic: bool = False) -> str:
    prefix = "../" if from_topic else ""
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)}</title>
  <link rel="stylesheet" href="{css_href(from_topic)}">
  <script defer src="{js_href(from_topic)}"></script>
  <script>
    window.MathJax = {{ tex: {{ inlineMath: [['$', '$'], ['\\\\(', '\\\\)']] }} }};
  </script>
  <script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
</head>
<body>
  <header class="site-header">
    <a class="brand" href="{prefix}index.html">Fluid Universe Reader</a>
    <nav>
      <a href="{prefix}index.html#topics">Topics</a>
      <a href="{prefix}future-paper-worklist.html">Paper Backfill</a>
      <a href="{prefix}data/topics.json">Data</a>
    </nav>
  </header>
  {body}
  <footer class="site-footer">
    <p>{escape(GENERATED_NOTICE)} Reader pages are additive commentary; published research papers remain the frozen references.</p>
  </footer>
</body>
</html>
"""


def render_topic_graph(topics: list[dict[str, Any]], from_topic: bool = False) -> str:
    tracks = {track: idx for idx, track in enumerate(sorted({topic.get("track", "Other") for topic in topics}))}
    positions: dict[str, tuple[int, int]] = {}
    by_track: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for topic in topics:
        by_track[topic.get("track", "Other")].append(topic)
    for track, grouped in by_track.items():
        for row, topic in enumerate(grouped):
            positions[topic["slug"]] = (80 + tracks[track] * 220, 70 + row * 95)
    width = max(760, 160 + max((x for x, _ in positions.values()), default=0))
    height = max(360, 120 + max((y for _, y in positions.values()), default=0))
    lines = [f'<svg class="topic-graph" viewBox="0 0 {width} {height}" role="img" aria-label="Readable topic graph">']
    for topic in topics:
        x1, y1 = positions[topic["slug"]]
        for rel in topic.get("related_topics", []) or []:
            if rel not in positions:
                continue
            x2, y2 = positions[rel]
            lines.append(f'<line x1="{x1 + 70}" y1="{y1 + 22}" x2="{x2 + 70}" y2="{y2 + 22}" />')
    for topic in topics:
        x, y = positions[topic["slug"]]
        href = topic_href(topic["slug"], from_topic)
        status = topic.get("status", "")
        label = compact(topic["title"], 30)
        lines.append(f'<a href="{escape(href)}">')
        lines.append(f'<g class="topic-node {escape(status)}">')
        lines.append(f'<rect x="{x}" y="{y}" width="160" height="54" rx="6" />')
        lines.append(f'<text x="{x + 12}" y="{y + 22}">{escape(label)}</text>')
        lines.append(f'<text x="{x + 12}" y="{y + 40}" class="track-label">{escape(topic.get("track", ""))}</text>')
        lines.append("</g></a>")
    lines.append("</svg>")
    return "\n".join(lines)


def render_index(
    topics_config: dict[str, Any],
    graph: dict[str, Any],
    topic_summaries: dict[str, dict[str, Any]],
    generated_at: str,
) -> str:
    topics = topics_config["topics"]
    pending_count = sum(len(summary["pending_sources"]) for summary in topic_summaries.values())
    status_counts = Counter(topic.get("status", "unknown") for topic in topics)
    cards = []
    for topic in topics:
        summary = topic_summaries[topic["slug"]]
        cards.append(
            "<article class=\"topic-card\" data-topic-card "
            f"data-track=\"{escape(topic.get('track', ''))}\" data-status=\"{escape(topic.get('status', ''))}\">"
            f"<a href=\"{escape(topic_href(topic['slug']))}\"><h3>{escape(topic['title'])}</h3></a>"
            f"<p>{escape(topic['summary'])}</p>"
            f"<p><span class=\"badge\">{escape(STATUS_EXPLAINERS.get(topic.get('status'), topic.get('status', '')))}</span> "
            f"<span class=\"muted\">{summary['node_count']} graph nodes · {len(summary['published_sources'])} published sources · {len(summary['pending_sources'])} pending sources</span></p>"
            "</article>"
        )
    status_items = "".join(
        f"<li><span>{escape(STATUS_EXPLAINERS.get(status, status))}</span><strong>{count}</strong></li>"
        for status, count in sorted(status_counts.items())
    )
    body = f"""
<main>
  <section class="reader-hero">
    <div>
      <p class="eyebrow">Generated {escape(generated_at)}</p>
      <h1>{escape(topics_config['site']['title'])}</h1>
      <p class="lede">{escape(topics_config['site']['subtitle'])}</p>
      <p>{escape(topics_config['site']['status_policy'])}</p>
    </div>
  </section>

  <section class="stats-band">
    <ul>
      <li><span>Graph nodes</span><strong>{len(graph['nodes'])}</strong></li>
      <li><span>Graph edges</span><strong>{len(graph['edges'])}</strong></li>
      <li><span>Reader topics</span><strong>{len(topics)}</strong></li>
      <li><span>Pending citations</span><strong>{pending_count}</strong></li>
    </ul>
  </section>

  <section class="layout-two">
    <div>
      <h2>Topic Graph</h2>
      <p>The public site is organized around readable topic pages, not raw graph nodes. Edges show suggested reading flow and conceptual dependencies.</p>
      {render_topic_graph(topics)}
    </div>
    <aside>
      <h2>Status Mix</h2>
      <ul class="status-list">{status_items}</ul>
      <p class="muted">Material sourced only from notes is tagged as paper pending and exported for later backfill.</p>
    </aside>
  </section>

  <section id="topics">
    <div class="section-heading">
      <h2>Readable Topics</h2>
      <input type="search" data-topic-filter placeholder="Filter topics">
    </div>
    <div class="topic-grid">
      {''.join(cards)}
    </div>
  </section>
</main>
"""
    return page_shell(topics_config["site"]["title"], body)


def render_topic_page(
    topic: dict[str, Any],
    topics: list[dict[str, Any]],
    nodes: list[dict[str, Any]],
    sources: list[dict[str, Any]],
    source_nodes: list[dict[str, Any]],
    generated_at: str,
) -> str:
    claims = node_group(nodes, is_claim)
    open_nodes = node_group(nodes, is_open)
    firewalls = node_group(nodes, is_firewall)
    physical = [node for node in nodes if str(node.get("layer", "")).startswith("physical")]
    math_nodes = [node for node in nodes if node.get("layer") == "math_object"]
    pending = [source for source in sources if source["status"] == "paper_pending"]
    related_links = " ".join(
        f'<a class="pill-link" href="{escape(topic_href(slug, from_topic=True))}">{escape(next((t["title"] for t in topics if t["slug"] == slug), slug))}</a>'
        for slug in topic.get("related_topics", []) or []
    )
    sections = topic.get("sections", {})
    body = f"""
<main>
  <div class="breadcrumb"><a href="../index.html">Reader</a> / {escape(topic['title'])}</div>
  <article class="topic-page">
    <header class="topic-header">
      <p class="eyebrow">{escape(topic.get('track', ''))} · generated {escape(generated_at)}</p>
      <h1>{escape(topic['title'])}</h1>
      <p class="lede">{escape(topic['summary'])}</p>
      <p>
        <span class="badge {escape(topic.get('status', ''))}">{escape(STATUS_EXPLAINERS.get(topic.get('status'), topic.get('status', '')))}</span>
        <span class="muted">{len(nodes)} selected graph nodes · {len(pending)} paper-pending source entries</span>
      </p>
    </header>

    <section>
      <h2>Physical Picture</h2>
      {render_paragraphs(sections.get('physical_picture', []))}
    </section>

    <section>
      <h2>Mathematical Core</h2>
      {render_paragraphs(sections.get('mathematical_core', []))}
    </section>

    <section>
      <h2>Status</h2>
      {render_paragraphs(sections.get('status_note', []))}
    </section>

    <section>
      <h2>Key Equations</h2>
      <div class="evidence-grid">{render_equations(nodes)}</div>
    </section>

    <section>
      <h2>Physical and Math Objects</h2>
      <div class="split">
        <div><h3>Physical objects</h3>{render_node_list(physical, 'No physical object nodes selected yet.')}</div>
        <div><h3>Math objects</h3>{render_node_list(math_nodes, 'No math object nodes selected yet.')}</div>
      </div>
    </section>

    <section>
      <h2>Claims</h2>
      {render_node_list(claims, 'No claim nodes selected yet.')}
    </section>

    <section>
      <h2>Open Gates and Firewalls</h2>
      <div class="split">
        <div><h3>Open or conditional gates</h3>{render_node_list(open_nodes, 'No open gates selected for this topic.')}</div>
        <div><h3>Firewalls</h3>{render_node_list(firewalls, 'No firewalls selected for this topic.')}</div>
      </div>
    </section>

    <section>
      <h2>Sources</h2>
      {render_sources(sources)}
    </section>

    <section>
      <h2>Source Anchor Nodes</h2>
      <details class="provenance-details">
        <summary>Show graph source-anchor nodes ({len(source_nodes)})</summary>
        {render_node_list(source_nodes, 'No adjacent source-anchor nodes found.')}
      </details>
    </section>

    <section>
      <h2>Related Topics</h2>
      <p>{related_links or '<span class="muted">No related topics configured.</span>'}</p>
    </section>
  </article>
</main>
"""
    return page_shell(topic["title"], body, from_topic=True)


def render_backfill_page(topics: list[dict[str, Any]], backfill: list[dict[str, Any]], generated_at: str) -> str:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for item in backfill:
        grouped[item["topic_slug"]].append(item)
    sections = []
    for topic in topics:
        items = grouped.get(topic["slug"], [])
        if not items:
            continue
        sections.append(f"<section><h2>{escape(topic['title'])}</h2><ul class=\"source-list pending\">")
        for item in items:
            nodes = ", ".join(item.get("node_ids", []))
            sections.append(
                "<li>"
                f"<span class=\"badge pending-badge\">paper pending</span> <code>{escape(item['source'])}</code>"
                f"<p>{escape(item.get('reason') or 'Add a formal paper citation when available.')}</p>"
                f"<p class=\"muted\">Nodes: {escape(nodes)}</p>"
                "</li>"
            )
        sections.append("</ul></section>")
    body = f"""
<main>
  <article class="topic-page">
    <header class="topic-header">
      <p class="eyebrow">Generated {escape(generated_at)}</p>
      <h1>Future Paper Backfill</h1>
      <p class="lede">Notes-backed material is visible by design. This page lists the sources that should be replaced or supplemented once formal papers are published.</p>
      <p><a href="data/future_paper_backfill.yaml">Download YAML worklist</a></p>
    </header>
    {''.join(sections) if sections else '<p>No paper-pending sources found.</p>'}
  </article>
</main>
"""
    return page_shell("Future Paper Backfill", body)


def render_css() -> str:
    return """
:root {
  color-scheme: light;
  --bg: #f7f8f5;
  --paper: #ffffff;
  --ink: #17201b;
  --muted: #5f6b62;
  --line: #d9dfd6;
  --accent: #235f7a;
  --accent-2: #6b4e16;
  --pending: #9a4b13;
  --good: #2f6b45;
  --soft: #edf2ee;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--bg);
  color: var(--ink);
  font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  line-height: 1.55;
}
a { color: var(--accent); text-decoration-thickness: 1px; text-underline-offset: 3px; }
code {
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: 0.92em;
  background: #eef1ed;
  padding: 0.08rem 0.22rem;
  border-radius: 4px;
}
.site-header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 1rem;
  padding: 0.85rem clamp(1rem, 4vw, 3rem);
  border-bottom: 1px solid var(--line);
  background: rgba(255, 255, 255, 0.9);
  position: sticky;
  top: 0;
  z-index: 10;
  backdrop-filter: blur(8px);
}
.brand { font-weight: 700; color: var(--ink); text-decoration: none; }
nav { display: flex; gap: 1rem; flex-wrap: wrap; font-size: 0.95rem; }
main { width: min(1180px, calc(100% - 2rem)); margin: 0 auto; }
.reader-hero {
  min-height: 42vh;
  display: grid;
  align-items: end;
  padding: 5rem 0 3rem;
  border-bottom: 1px solid var(--line);
  background:
    linear-gradient(90deg, rgba(247,248,245,0.92), rgba(247,248,245,0.68)),
    radial-gradient(circle at 85% 30%, rgba(35,95,122,0.18), transparent 34%),
    linear-gradient(135deg, rgba(47,107,69,0.12), rgba(107,78,22,0.10));
}
.reader-hero > div { width: min(820px, 100%); }
.eyebrow { color: var(--accent-2); text-transform: uppercase; letter-spacing: 0.08em; font-size: 0.78rem; font-weight: 700; }
h1 { font-size: clamp(2.2rem, 5vw, 4.8rem); line-height: 1.02; margin: 0.35rem 0 1rem; }
h2 { font-size: 1.45rem; margin-top: 2rem; }
h3 { margin-bottom: 0.35rem; }
.lede { font-size: 1.2rem; max-width: 760px; color: #2e3932; }
.stats-band ul, .status-list {
  list-style: none;
  padding: 0;
  margin: 0;
}
.stats-band ul {
  display: grid;
  grid-template-columns: repeat(4, minmax(0, 1fr));
  gap: 1px;
  border: 1px solid var(--line);
  background: var(--line);
  margin: 1.5rem 0 2rem;
}
.stats-band li, .status-list li {
  background: var(--paper);
  padding: 1rem;
  display: flex;
  justify-content: space-between;
  gap: 1rem;
}
.stats-band strong { font-size: 1.45rem; }
.layout-two {
  display: grid;
  grid-template-columns: minmax(0, 1fr) 280px;
  gap: 2rem;
  align-items: start;
}
aside {
  background: var(--paper);
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 1rem;
}
.section-heading {
  display: flex;
  justify-content: space-between;
  align-items: center;
  gap: 1rem;
  margin-top: 2rem;
}
input[type="search"] {
  width: min(360px, 100%);
  border: 1px solid var(--line);
  border-radius: 6px;
  padding: 0.7rem 0.85rem;
  font: inherit;
}
.topic-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
  gap: 1rem;
  margin-bottom: 3rem;
}
.topic-card, .evidence-card {
  background: var(--paper);
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 1rem;
}
.topic-card h3 { font-size: 1.08rem; }
.badge {
  display: inline-flex;
  align-items: center;
  border: 1px solid var(--line);
  border-radius: 999px;
  padding: 0.16rem 0.5rem;
  background: var(--soft);
  color: var(--ink);
  font-size: 0.82rem;
  font-weight: 650;
}
.small { font-size: 0.72rem; }
.pending-badge, .notes_backed_pending_papers { border-color: #e1b78d; background: #fff2e4; color: var(--pending); }
.published_backed { border-color: #b9d5c3; background: #eef8f1; color: var(--good); }
.mixed_published_and_open { border-color: #d8c891; background: #fff9df; color: #6b4e16; }
.muted { color: var(--muted); }
.topic-graph {
  width: 100%;
  height: auto;
  background: var(--paper);
  border: 1px solid var(--line);
  border-radius: 8px;
}
.topic-graph line { stroke: #b8c2b7; stroke-width: 1.3; }
.topic-graph rect { fill: #f9fbf8; stroke: #9aa79d; stroke-width: 1.1; }
.topic-graph text { fill: var(--ink); font-size: 12px; font-weight: 650; }
.topic-graph .track-label { fill: var(--muted); font-size: 10px; font-weight: 500; }
.breadcrumb { margin: 1rem 0; color: var(--muted); }
.topic-page {
  background: var(--paper);
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: clamp(1rem, 3vw, 2.4rem);
  margin-bottom: 3rem;
}
.topic-header {
  border-bottom: 1px solid var(--line);
  margin-bottom: 1.5rem;
  padding-bottom: 1.5rem;
}
.split {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 1rem;
}
.node-list, .source-list {
  list-style: none;
  padding: 0;
  margin: 0;
}
.node-list li, .source-list li {
  border-top: 1px solid var(--line);
  padding: 0.85rem 0;
}
.node-list p, .source-list p { margin: 0.35rem 0 0; }
.evidence-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
  gap: 1rem;
}
.math-block {
  overflow-x: auto;
  padding: 0.8rem;
  background: #f4f6f3;
  border: 1px solid var(--line);
  border-radius: 6px;
}
.node-id {
  color: var(--muted);
  font-size: 0.9rem;
}
.pill-link {
  display: inline-block;
  margin: 0.2rem 0.35rem 0.2rem 0;
  padding: 0.32rem 0.55rem;
  border: 1px solid var(--line);
  border-radius: 999px;
  background: #f4f6f3;
}
.provenance-details {
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 0.85rem 1rem;
  background: #fafbf9;
}
.provenance-details summary {
  cursor: pointer;
  font-weight: 700;
  color: var(--accent);
}
.site-footer {
  border-top: 1px solid var(--line);
  padding: 1.2rem clamp(1rem, 4vw, 3rem);
  color: var(--muted);
}
@media (max-width: 760px) {
  .layout-two, .split { grid-template-columns: 1fr; }
  .stats-band ul { grid-template-columns: repeat(2, minmax(0, 1fr)); }
  .site-header { align-items: flex-start; flex-direction: column; }
}
"""


def render_js() -> str:
    return """
const filter = document.querySelector('[data-topic-filter]');
if (filter) {
  filter.addEventListener('input', () => {
    const query = filter.value.trim().toLowerCase();
    document.querySelectorAll('[data-topic-card]').forEach((card) => {
      card.hidden = query && !card.textContent.toLowerCase().includes(query);
    });
  });
}
"""


def build_topic_summaries(
    topics: list[dict[str, Any]],
    node_by_id: dict[str, dict[str, Any]],
    incoming: dict[str, list[dict[str, Any]]],
    outgoing: dict[str, list[dict[str, Any]]],
    resolution: dict[str, dict[str, Any]],
    doi_entries: list[dict[str, str]],
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    summaries = {}
    backfill: list[dict[str, Any]] = []
    for topic in topics:
        nodes = topic_nodes(topic, node_by_id)
        sources, _source_nodes = collect_sources_for_topic(topic, node_by_id, incoming, outgoing, resolution, doi_entries)
        published = [source for source in sources if source["status"] == "published"]
        pending = [source for source in sources if source["status"] == "paper_pending"]
        summaries[topic["slug"]] = {
            "node_count": len(nodes),
            "published_sources": published,
            "pending_sources": pending,
        }
        pending_by_path: dict[str, set[str]] = defaultdict(set)
        reason_by_path: dict[str, str] = {}
        for source in pending:
            pending_by_path[source["path"]].add(source["node_id"])
            reason_by_path[source["path"]] = source.get("reason", "")
        for path, node_ids in sorted(pending_by_path.items()):
            backfill.append(
                {
                    "topic_slug": topic["slug"],
                    "topic_title": topic["title"],
                    "source": path,
                    "reason": reason_by_path.get(path, ""),
                    "node_ids": sorted(node_ids),
                    "tag": "future_paper_needed",
                }
            )
    return summaries, backfill


def write_manifest(site_dir: Path, generated_files: set[Path]) -> None:
    entries = sorted(path.relative_to(site_dir).as_posix() for path in generated_files if path.exists())
    (site_dir / MANIFEST).write_text("# Generated by atlas/scripts/generate_reader_site.py\n" + "\n".join(entries) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=GRAPH_YAML)
    parser.add_argument("--topics", type=Path, default=TOPICS_YAML)
    parser.add_argument("--readme", type=Path, default=README)
    parser.add_argument("--site", type=Path, default=SITE_DIR)
    parser.add_argument("--exports", type=Path, default=EXPORT_DIR)
    args = parser.parse_args()

    graph = load_yaml(args.graph)
    topics_config = load_yaml(args.topics)
    doi_entries = parse_readme_dois(args.readme)
    resolution = {entry["legacy_ref"]: entry for entry in graph.get("source_resolution", []) if entry.get("legacy_ref")}
    node_by_id, incoming, outgoing = build_indexes(graph)
    timestamp = generated_utc()

    clean_previous_site(args.site)
    args.site.mkdir(parents=True, exist_ok=True)
    args.exports.mkdir(parents=True, exist_ok=True)
    generated_files: set[Path] = set()

    summaries, backfill = build_topic_summaries(
        topics_config["topics"],
        node_by_id,
        incoming,
        outgoing,
        resolution,
        doi_entries,
    )

    write_file(args.site / "assets" / "reader.css", render_css(), generated_files)
    write_file(args.site / "assets" / "reader.js", render_js(), generated_files)
    write_file(
        args.site / "data" / "topics.json",
        json.dumps({"topics": topics_config["topics"], "generated_utc": timestamp}, indent=2) + "\n",
        generated_files,
    )
    write_file(
        args.site / "data" / "future_paper_backfill.yaml",
        yaml.safe_dump(backfill, sort_keys=False, allow_unicode=True),
        generated_files,
    )
    (args.exports / "future_paper_backfill.yaml").write_text(
        yaml.safe_dump(backfill, sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )

    write_file(args.site / "index.html", render_index(topics_config, graph, summaries, timestamp), generated_files)
    write_file(args.site / "future-paper-worklist.html", render_backfill_page(topics_config["topics"], backfill, timestamp), generated_files)

    for topic in topics_config["topics"]:
        nodes = topic_nodes(topic, node_by_id)
        sources, source_nodes = collect_sources_for_topic(topic, node_by_id, incoming, outgoing, resolution, doi_entries)
        write_file(
            args.site / slug_to_file(topic["slug"]),
            render_topic_page(topic, topics_config["topics"], nodes, sources, source_nodes, timestamp),
            generated_files,
        )

    write_manifest(args.site, generated_files)
    print(f"generated reader topics: {len(topics_config['topics'])}")
    print(f"paper-pending backfill entries: {len(backfill)}")
    print(f"site root: {args.site}")


if __name__ == "__main__":
    main()
