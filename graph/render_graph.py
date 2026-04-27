#!/usr/bin/env python3
"""Render the Fluid Universe derivation atlas as DOT and standalone HTML.

The HTML output has no external dependencies. It embeds the graph JSON and uses
Canvas plus a small force layout so it can be opened directly in a browser.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml


GRAPH_DIR = Path(__file__).resolve().parent
DEFAULT_GRAPH = GRAPH_DIR / "fluid_universe_derivation_atlas_graph.yaml"
DEFAULT_DOT = GRAPH_DIR / "fluid_universe_derivation_atlas_graph.dot"
DEFAULT_HTML = GRAPH_DIR / "fluid_universe_derivation_atlas_graph.html"
DEFAULT_SVG = GRAPH_DIR / "fluid_universe_derivation_atlas_graph.svg"
PALETTE = [
    "#0f766e", "#2563eb", "#9333ea", "#c2410c", "#15803d", "#be123c",
    "#0891b2", "#7c3aed", "#b45309", "#4f46e5", "#047857", "#a21caf",
    "#0369a1", "#65a30d", "#dc2626", "#475569",
]


def dot_quote(value: object) -> str:
    text = str(value if value is not None else "")
    text = text.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")
    return f'"{text}"'


def xml_escape(value: object) -> str:
    text = str(value if value is not None else "")
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#39;")
    )


def truncate(value: object, limit: int = 36) -> str:
    text = str(value if value is not None else "")
    return text if len(text) <= limit else f"{text[:limit - 1]}..."


def node_label(node: dict) -> str:
    label = node.get("label") or node.get("id")
    node_id = node.get("id", "")
    if label == node_id:
        return str(label)
    return f"{label}\\n{node_id}"


def dot_node_label(node: dict) -> str:
    return str(node.get("id", ""))


def write_dot(graph: dict, output_path: Path) -> None:
    layers = sorted({node.get("layer", "unknown") for node in graph["nodes"]})
    layer_index = {layer: index for index, layer in enumerate(layers)}

    lines = [
        "digraph FluidUniverseDerivationAtlas {",
        "  graph [layout=sfdp, overlap=false, splines=true, bgcolor=\"white\", outputorder=edgesfirst];",
        "  node [shape=box, style=\"rounded,filled\", fontname=\"Helvetica\", fontsize=10];",
        "  edge [fontname=\"Helvetica\", fontsize=8, color=\"#64748b\", arrowsize=0.6];",
    ]

    for node in graph["nodes"]:
        layer = node.get("layer", "unknown")
        attrs = {
            "label": dot_node_label(node),
            "fillcolor": PALETTE[layer_index[layer] % len(PALETTE)],
            "tooltip": node_label(node),
            "group": layer,
        }
        attr_text = ", ".join(f"{key}={dot_quote(value)}" for key, value in attrs.items())
        lines.append(f"  {dot_quote(node['id'])} [{attr_text}];")

    for edge in graph["edges"]:
        tooltip = edge.get("relation", "")
        if edge.get("note"):
            tooltip = f"{tooltip}: {edge['note']}"
        lines.append(
            f"  {dot_quote(edge['source'])} -> {dot_quote(edge['target'])} "
            f"[tooltip={dot_quote(tooltip)}];"
        )

    lines.append("}")
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def node_tooltip(node: dict) -> str:
    fields = [
        ("id", node.get("id")),
        ("label", node.get("label")),
        ("layer", node.get("layer")),
        ("type", node.get("type")),
        ("status", node.get("status")),
        ("source_kind", node.get("source_kind")),
        ("file", node.get("file")),
        ("meaning", node.get("meaning")),
        ("tex_anchor", tex_anchor_summary(node)),
        ("resolution_note", node.get("resolution_note")),
    ]
    return "\n".join(f"{key}: {value}" for key, value in fields if value)


def tex_anchor_summary(node: dict) -> str:
    anchor = node.get("tex_anchor")
    if not anchor:
        return ""
    label = anchor.get("nearest_label") or {}
    label_text = f" label={label.get('name')}" if label.get("name") else ""
    heading = anchor.get("heading") or ""
    heading_text = f" heading={heading}" if heading else ""
    return f"{anchor.get('file')}:{anchor.get('line')}{label_text}{heading_text}"


def write_svg(graph: dict, output_path: Path) -> None:
    layers = sorted({node.get("layer", "unknown") for node in graph["nodes"]})
    layer_index = {layer: index for index, layer in enumerate(layers)}
    by_layer = {layer: [] for layer in layers}
    for node in graph["nodes"]:
        by_layer[node.get("layer", "unknown")].append(node)
    for group in by_layer.values():
        group.sort(key=lambda node: node.get("id", ""))

    margin_x = 70
    margin_y = 86
    column_gap = 270
    row_gap = 28
    label_width = 178
    width = max(900, margin_x * 2 + column_gap * max(0, len(layers) - 1) + label_width)
    height = max(760, margin_y * 2 + row_gap * max(len(group) for group in by_layer.values()))

    positions = {}
    for x_index, layer in enumerate(layers):
        x = margin_x + x_index * column_gap
        group = by_layer[layer]
        for y_index, node in enumerate(group):
            positions[node["id"]] = (x, margin_y + y_index * row_gap)

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "  <style>",
        "    .bg { fill: #fbfdff; }",
        "    .layer { fill: #0f172a; font: 700 13px system-ui, sans-serif; }",
        "    .node-label { fill: #334155; font: 9px system-ui, sans-serif; }",
        "    .edge { stroke: #94a3b8; stroke-opacity: 0.16; stroke-width: 0.8; }",
        "    .node { stroke: #ffffff; stroke-width: 1.4; }",
        "    .future { stroke: #be123c; stroke-width: 2; }",
        "  </style>",
        '  <rect class="bg" width="100%" height="100%"/>',
        "  <title>Fluid Universe Derivation Atlas static overview</title>",
    ]

    for layer in layers:
        x = positions[by_layer[layer][0]["id"]][0] if by_layer[layer] else margin_x
        lines.append(
            f'  <text class="layer" x="{x}" y="42">{xml_escape(truncate(layer, 28))}</text>'
        )

    for edge in graph["edges"]:
        source = positions.get(edge.get("source"))
        target = positions.get(edge.get("target"))
        if not source or not target:
            continue
        tooltip = edge.get("relation", "")
        if edge.get("note"):
            tooltip = f"{tooltip}: {edge['note']}"
        lines.extend(
            [
                f'  <line class="edge" x1="{source[0]}" y1="{source[1]}" x2="{target[0]}" y2="{target[1]}">',
                f"    <title>{xml_escape(tooltip)}</title>",
                "  </line>",
            ]
        )

    for layer in layers:
        fill = PALETTE[layer_index[layer] % len(PALETTE)]
        for node in by_layer[layer]:
            x, y = positions[node["id"]]
            classes = "node future" if node.get("source_kind") == "future_paper_note" else "node"
            lines.extend(
                [
                    f'  <g><title>{xml_escape(node_tooltip(node))}</title>',
                    f'    <circle class="{classes}" cx="{x}" cy="{y}" r="5.2" fill="{fill}"/>',
                    f'    <text class="node-label" x="{x + 9}" y="{y + 3}">{xml_escape(truncate(node["id"]))}</text>',
                    "  </g>",
                ]
            )

    lines.append("</svg>")
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Fluid Universe Derivation Atlas</title>
  <style>
    :root {
      color-scheme: light;
      --bg: #f8fafc;
      --panel: #ffffff;
      --text: #0f172a;
      --muted: #64748b;
      --line: #dbe3ee;
      --accent: #0f766e;
    }
    * { box-sizing: border-box; }
    html, body { height: 100%; margin: 0; }
    body {
      font: 14px/1.4 system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: var(--text);
      background: var(--bg);
      overflow: hidden;
    }
    #app {
      display: grid;
      grid-template-columns: 340px 1fr;
      height: 100vh;
      min-height: 0;
    }
    aside {
      border-right: 1px solid var(--line);
      background: var(--panel);
      padding: 14px;
      overflow: auto;
    }
    main { position: relative; min-width: 0; min-height: 0; }
    canvas { display: block; width: 100%; height: 100%; background: #fbfdff; }
    h1 { margin: 0 0 10px; font-size: 18px; line-height: 1.2; }
    label { display: block; margin: 11px 0 5px; font-size: 12px; color: var(--muted); }
    input, select, button {
      width: 100%;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: white;
      color: var(--text);
      padding: 8px 9px;
      font: inherit;
    }
    button {
      margin-top: 12px;
      cursor: pointer;
      background: var(--accent);
      color: white;
      border-color: var(--accent);
    }
    .stats {
      display: grid;
      grid-template-columns: repeat(2, 1fr);
      gap: 8px;
      margin: 12px 0;
    }
    .stat {
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 8px;
      background: #f8fafc;
    }
    .stat strong { display: block; font-size: 18px; }
    #details {
      margin-top: 12px;
      border-top: 1px solid var(--line);
      padding-top: 12px;
      white-space: pre-wrap;
      word-break: break-word;
    }
    #legend {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      margin-top: 10px;
    }
    .chip {
      display: inline-flex;
      align-items: center;
      gap: 5px;
      border: 1px solid var(--line);
      border-radius: 999px;
      padding: 3px 7px;
      font-size: 12px;
      color: var(--muted);
      background: #fff;
    }
    .swatch { width: 9px; height: 9px; border-radius: 50%; }
    #hint {
      position: absolute;
      left: 12px;
      bottom: 10px;
      color: var(--muted);
      background: rgba(255,255,255,0.88);
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 6px 8px;
      pointer-events: none;
    }
  </style>
</head>
<body>
<div id="app">
  <aside>
    <h1>Derivation Atlas</h1>
    <div class="stats">
      <div class="stat"><strong id="visibleNodes">0</strong> nodes</div>
      <div class="stat"><strong id="visibleEdges">0</strong> edges</div>
    </div>

    <label for="search">Search ID, label, status, or source</label>
    <input id="search" placeholder="e.g. 4PN, OPEN_QUAD, Maxwell">

    <label for="layer">Layer</label>
    <select id="layer"></select>

    <label for="kind">Source Kind</label>
    <select id="kind"></select>

    <label for="relation">Relation</label>
    <select id="relation"></select>

    <button id="reset">Reset View</button>

    <div id="legend"></div>
    <div id="details">Hover or click a node to inspect it. Drag to pan, wheel to zoom.</div>
  </aside>
  <main>
    <canvas id="graph"></canvas>
    <div id="hint">Click a node to lock its neighborhood highlight.</div>
  </main>
</div>

<script>
const GRAPH = __GRAPH_DATA__;

const canvas = document.getElementById("graph");
const ctx = canvas.getContext("2d");
const searchInput = document.getElementById("search");
const layerSelect = document.getElementById("layer");
const kindSelect = document.getElementById("kind");
const relationSelect = document.getElementById("relation");
const details = document.getElementById("details");
const visibleNodesEl = document.getElementById("visibleNodes");
const visibleEdgesEl = document.getElementById("visibleEdges");
const legend = document.getElementById("legend");

const nodes = GRAPH.nodes.map((node, index) => ({
  ...node,
  index,
  x: 0,
  y: 0,
  vx: 0,
  vy: 0,
  r: 4 + Math.min(7, Math.sqrt((node.anchors_claims || node.sources || []).length + 1))
}));
const nodeById = new Map(nodes.map(node => [node.id, node]));
const edges = GRAPH.edges
  .map(edge => ({...edge, sourceNode: nodeById.get(edge.source), targetNode: nodeById.get(edge.target)}))
  .filter(edge => edge.sourceNode && edge.targetNode);

const layers = unique(nodes.map(node => node.layer || "unknown"));
const kinds = unique(nodes.map(node => node.source_kind || "none"));
const relations = unique(edges.map(edge => edge.relation || "unknown"));
const palette = [
  "#0f766e", "#2563eb", "#9333ea", "#c2410c", "#15803d", "#be123c",
  "#0891b2", "#7c3aed", "#b45309", "#4f46e5", "#047857", "#a21caf",
  "#0369a1", "#65a30d", "#dc2626", "#475569"
];
const layerColor = new Map(layers.map((layer, index) => [layer, palette[index % palette.length]]));

let transform = {x: 0, y: 0, k: 1};
let dragging = false;
let lastMouse = {x: 0, y: 0};
let hoverNode = null;
let selectedNode = null;
let visibleNodes = [];
let visibleEdges = [];
let running = true;

function unique(values) {
  return [...new Set(values)].sort((a, b) => a.localeCompare(b));
}

function fillSelect(select, values, allLabel) {
  select.innerHTML = "";
  const all = document.createElement("option");
  all.value = "";
  all.textContent = allLabel;
  select.appendChild(all);
  values.forEach(value => {
    const option = document.createElement("option");
    option.value = value;
    option.textContent = value;
    select.appendChild(option);
  });
}

function initializePositions() {
  const byLayer = new Map();
  nodes.forEach(node => {
    const layer = node.layer || "unknown";
    if (!byLayer.has(layer)) byLayer.set(layer, []);
    byLayer.get(layer).push(node);
  });

  const ringRadius = 245;
  const layerRadius = 48;
  layers.forEach((layer, layerIndex) => {
    const group = byLayer.get(layer) || [];
    const centerAngle = (Math.PI * 2 * layerIndex) / Math.max(1, layers.length);
    const cx = Math.cos(centerAngle) * ringRadius;
    const cy = Math.sin(centerAngle) * ringRadius;
    group.forEach((node, index) => {
      const angle = (Math.PI * 2 * index) / Math.max(1, group.length);
      const radius = layerRadius + 7 * Math.sqrt(group.length);
      node.x = cx + Math.cos(angle) * radius;
      node.y = cy + Math.sin(angle) * radius;
    });
  });
}

function rebuildLegend() {
  legend.innerHTML = "";
  layers.forEach(layer => {
    const chip = document.createElement("span");
    chip.className = "chip";
    chip.innerHTML = `<span class="swatch" style="background:${layerColor.get(layer)}"></span>${escapeHtml(layer)}`;
    legend.appendChild(chip);
  });
}

function escapeHtml(value) {
  return String(value || "").replace(/[&<>"']/g, ch => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;", "'": "&#39;"
  }[ch]));
}

function nodeSearchText(node) {
  return [
    node.id, node.label, node.layer, node.type, node.status, node.source_kind,
    node.meaning, node.file, ...(node.sources || []), ...(node.files || [])
  ].filter(Boolean).join(" ").toLowerCase();
}

function applyFilters() {
  const query = searchInput.value.trim().toLowerCase();
  const layer = layerSelect.value;
  const kind = kindSelect.value;
  const relation = relationSelect.value;

  visibleNodes = nodes.filter(node => {
    if (layer && (node.layer || "unknown") !== layer) return false;
    if (kind && (node.source_kind || "none") !== kind) return false;
    if (query && !nodeSearchText(node).includes(query)) return false;
    return true;
  });
  const visibleSet = new Set(visibleNodes.map(node => node.id));
  visibleEdges = edges.filter(edge => {
    if (!visibleSet.has(edge.source) || !visibleSet.has(edge.target)) return false;
    if (relation && edge.relation !== relation) return false;
    return true;
  });

  visibleNodesEl.textContent = visibleNodes.length;
  visibleEdgesEl.textContent = visibleEdges.length;
  running = true;
}

function resize() {
  const rect = canvas.getBoundingClientRect();
  const scale = window.devicePixelRatio || 1;
  canvas.width = Math.max(1, Math.floor(rect.width * scale));
  canvas.height = Math.max(1, Math.floor(rect.height * scale));
  ctx.setTransform(scale, 0, 0, scale, 0, 0);
  draw();
}

function tick() {
  if (!running || visibleNodes.length === 0) return;
  const visibleSet = new Set(visibleNodes.map(node => node.id));
  const maxNodes = Math.min(visibleNodes.length, 450);

  for (let i = 0; i < maxNodes; i++) {
    const a = visibleNodes[i];
    for (let j = i + 1; j < maxNodes; j++) {
      const b = visibleNodes[j];
      let dx = a.x - b.x;
      let dy = a.y - b.y;
      let dist2 = dx * dx + dy * dy + 0.01;
      const force = Math.min(2.4, 1500 / dist2);
      const dist = Math.sqrt(dist2);
      dx /= dist;
      dy /= dist;
      a.vx += dx * force;
      a.vy += dy * force;
      b.vx -= dx * force;
      b.vy -= dy * force;
    }
  }

  visibleEdges.forEach(edge => {
    if (!visibleSet.has(edge.source) || !visibleSet.has(edge.target)) return;
    const a = edge.sourceNode;
    const b = edge.targetNode;
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const dist = Math.sqrt(dx * dx + dy * dy) || 1;
    const target = 95;
    const force = (dist - target) * 0.0025;
    const fx = (dx / dist) * force;
    const fy = (dy / dist) * force;
    a.vx += fx;
    a.vy += fy;
    b.vx -= fx;
    b.vy -= fy;
  });

  let totalMotion = 0;
  visibleNodes.forEach(node => {
    node.vx += -node.x * 0.0008;
    node.vy += -node.y * 0.0008;
    node.vx *= 0.82;
    node.vy *= 0.82;
    node.x += node.vx;
    node.y += node.vy;
    totalMotion += Math.abs(node.vx) + Math.abs(node.vy);
  });
  if (totalMotion < 0.05) running = false;
}

function worldToScreen(node) {
  const rect = canvas.getBoundingClientRect();
  return {
    x: rect.width / 2 + transform.x + node.x * transform.k,
    y: rect.height / 2 + transform.y + node.y * transform.k
  };
}

function screenToWorld(x, y) {
  const rect = canvas.getBoundingClientRect();
  return {
    x: (x - rect.width / 2 - transform.x) / transform.k,
    y: (y - rect.height / 2 - transform.y) / transform.k
  };
}

function neighborSet(node) {
  if (!node) return new Set();
  const set = new Set([node.id]);
  edges.forEach(edge => {
    if (edge.source === node.id) set.add(edge.target);
    if (edge.target === node.id) set.add(edge.source);
  });
  return set;
}

function draw() {
  const rect = canvas.getBoundingClientRect();
  ctx.clearRect(0, 0, rect.width, rect.height);
  const focus = selectedNode || hoverNode;
  const focusedNeighbors = neighborSet(focus);

  ctx.save();
  ctx.translate(rect.width / 2 + transform.x, rect.height / 2 + transform.y);
  ctx.scale(transform.k, transform.k);

  visibleEdges.forEach(edge => {
    const focused = !focus || (focusedNeighbors.has(edge.source) && focusedNeighbors.has(edge.target));
    ctx.beginPath();
    ctx.moveTo(edge.sourceNode.x, edge.sourceNode.y);
    ctx.lineTo(edge.targetNode.x, edge.targetNode.y);
    ctx.lineWidth = focused ? 1.1 / transform.k : 0.45 / transform.k;
    ctx.strokeStyle = focused ? "rgba(15, 23, 42, 0.32)" : "rgba(100, 116, 139, 0.08)";
    ctx.stroke();
  });

  visibleNodes.forEach(node => {
    const focused = !focus || focusedNeighbors.has(node.id);
    ctx.beginPath();
    ctx.arc(node.x, node.y, node.r / Math.sqrt(transform.k), 0, Math.PI * 2);
    ctx.fillStyle = focused ? layerColor.get(node.layer || "unknown") || "#64748b" : "rgba(148, 163, 184, 0.35)";
    ctx.fill();
    if (node === selectedNode) {
      ctx.lineWidth = 3 / transform.k;
      ctx.strokeStyle = "#0f172a";
      ctx.stroke();
    } else if (node === hoverNode) {
      ctx.lineWidth = 2 / transform.k;
      ctx.strokeStyle = "#0f766e";
      ctx.stroke();
    }
  });

  if (transform.k > 0.82) {
    ctx.font = `${11 / transform.k}px system-ui, sans-serif`;
    ctx.fillStyle = "#334155";
    visibleNodes.forEach(node => {
      if (focus && !focusedNeighbors.has(node.id)) return;
      ctx.fillText(node.id, node.x + 7 / transform.k, node.y - 7 / transform.k);
    });
  }
  ctx.restore();
}

function nearestNode(clientX, clientY) {
  const world = screenToWorld(clientX, clientY);
  let best = null;
  let bestDist = Infinity;
  visibleNodes.forEach(node => {
    const dx = node.x - world.x;
    const dy = node.y - world.y;
    const dist = Math.sqrt(dx * dx + dy * dy);
    if (dist < bestDist && dist < 12 / Math.sqrt(transform.k)) {
      best = node;
      bestDist = dist;
    }
  });
  return best;
}

function describeNode(node) {
  if (!node) {
    details.textContent = "Hover or click a node to inspect it. Drag to pan, wheel to zoom.";
    return;
  }
  const texAnchor = node.tex_anchor ? [
    `${node.tex_anchor.file}:${node.tex_anchor.line}`,
    node.tex_anchor.heading ? `heading=${node.tex_anchor.heading}` : "",
    node.tex_anchor.nearest_label ? `label=${node.tex_anchor.nearest_label.name}` : "",
    node.tex_anchor.match_basis ? `basis=${node.tex_anchor.match_basis}` : "",
    node.tex_anchor.confidence ? `confidence=${node.tex_anchor.confidence}` : ""
  ].filter(Boolean).join(" | ") : "";
  const fields = [
    ["id", node.id],
    ["label", node.label],
    ["layer", node.layer],
    ["type", node.type],
    ["status", node.status],
    ["source_kind", node.source_kind],
    ["file", node.file],
    ["meaning", node.meaning],
    ["tex_anchor", texAnchor],
    ["sources", (node.sources || []).join(", ")],
    ["legacy_file", node.legacy_file],
    ["resolution_note", node.resolution_note]
  ].filter(([, value]) => value);
  details.textContent = fields.map(([key, value]) => `${key}: ${value}`).join("\n");
}

function animate() {
  tick();
  draw();
  requestAnimationFrame(animate);
}

canvas.addEventListener("mousedown", event => {
  dragging = true;
  lastMouse = {x: event.clientX, y: event.clientY};
});
window.addEventListener("mouseup", () => { dragging = false; });
window.addEventListener("mousemove", event => {
  if (dragging) {
    transform.x += event.clientX - lastMouse.x;
    transform.y += event.clientY - lastMouse.y;
    lastMouse = {x: event.clientX, y: event.clientY};
    draw();
    return;
  }
  hoverNode = nearestNode(event.clientX, event.clientY);
  if (!selectedNode) describeNode(hoverNode);
  draw();
});
canvas.addEventListener("click", event => {
  selectedNode = nearestNode(event.clientX, event.clientY);
  describeNode(selectedNode || hoverNode);
  draw();
});
canvas.addEventListener("wheel", event => {
  event.preventDefault();
  const before = screenToWorld(event.clientX, event.clientY);
  const scale = Math.exp(-event.deltaY * 0.001);
  transform.k = Math.max(0.08, Math.min(5, transform.k * scale));
  const after = screenToWorld(event.clientX, event.clientY);
  transform.x += (after.x - before.x) * transform.k;
  transform.y += (after.y - before.y) * transform.k;
  draw();
}, {passive: false});

document.getElementById("reset").addEventListener("click", () => {
  transform = {x: 0, y: 0, k: 1};
  selectedNode = null;
  hoverNode = null;
  running = true;
  describeNode(null);
});

[searchInput, layerSelect, kindSelect, relationSelect].forEach(control => {
  control.addEventListener("input", () => {
    selectedNode = null;
    applyFilters();
    describeNode(null);
  });
});

fillSelect(layerSelect, layers, "All layers");
fillSelect(kindSelect, kinds, "All source kinds");
fillSelect(relationSelect, relations, "All relations");
initializePositions();
rebuildLegend();
applyFilters();
resize();
window.addEventListener("resize", resize);
animate();
</script>
</body>
</html>
"""


def write_html(graph: dict, output_path: Path) -> None:
    graph_data = json.dumps(graph, ensure_ascii=False, separators=(",", ":"))
    graph_data = graph_data.replace("</", "<\\/")
    output_path.write_text(HTML_TEMPLATE.replace("__GRAPH_DATA__", graph_data), encoding="utf-8")


def load_graph(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() in {".yaml", ".yml"}:
        return yaml.safe_load(text)
    return json.loads(text)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=DEFAULT_GRAPH, help="Input graph JSON path")
    parser.add_argument("--dot", type=Path, default=DEFAULT_DOT, help="Output Graphviz DOT path")
    parser.add_argument("--html", type=Path, default=DEFAULT_HTML, help="Output standalone HTML path")
    parser.add_argument("--svg", type=Path, default=DEFAULT_SVG, help="Output static SVG path")
    parser.add_argument("--skip-dot", action="store_true", help="Do not write DOT output")
    parser.add_argument("--skip-html", action="store_true", help="Do not write HTML output")
    parser.add_argument("--skip-svg", action="store_true", help="Do not write SVG output")
    args = parser.parse_args()

    graph = load_graph(args.graph)

    if not args.skip_dot:
        write_dot(graph, args.dot)
        print(f"wrote {args.dot}")

    if not args.skip_html:
        write_html(graph, args.html)
        print(f"wrote {args.html}")

    if not args.skip_svg:
        write_svg(graph, args.svg)
        print(f"wrote {args.svg}")


if __name__ == "__main__":
    main()
