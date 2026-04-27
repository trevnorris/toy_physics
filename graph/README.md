# Fluid Universe Derivation Atlas

This directory contains the graph-level provenance layer for the Fluid Universe
research program. The stable graph artifacts are versionless filenames; version
lineage is kept inside the files when useful.

## Core Files

- `fluid_universe_derivation_atlas_graph.yaml` - primary maintained graph.
- `fluid_universe_derivation_atlas_graph.json` - generated JSON mirror for tools and browser embedding.
- `fluid_universe_derivation_atlas_paper_insertion_manifest.json` - backlink insertion manifest.
- `fluid_universe_derivation_atlas_source_resolution.md` - map from legacy markdown references to canonical TeX papers or future-paper notes.
- `fluid_universe_derivation_atlas_paper_backlink_blocks.md` - paste-ready backlink blocks.
- `fluid_universe_derivation_atlas_status_firewall_register.md` - false-upgrade/status-firewall register.

## Format Workflow

Edit the YAML graph, then regenerate JSON:

```sh
python3 graph/sync_graph_formats.py
```

Check that JSON is current:

```sh
python3 graph/sync_graph_formats.py --check
```

The JSON file is intentionally kept because it is convenient for `jq`, browser
visualization, and external graph tooling. It should not be edited by hand.

## TeX Anchors

Graph nodes that point into maintained TeX papers carry `tex_anchor` metadata:

```yaml
tex_anchor:
  file: research/4d_4pn/paper/4d_4pn.tex
  line: 2453
  heading_level: section
  heading: Final conditional conservative 4PN theorem
  nearest_label:
    name: sec:final-conditional-4pn-theorem
    line: 2454
  match_basis: semantic_heading_match
  confidence: high
```

Regenerate those anchors after TeX headings or labels move:

```sh
python3 graph/annotate_tex_anchors.py
python3 graph/sync_graph_formats.py
```

The annotator parses TeX headings and `\label{...}` entries. Exact semantic
matches are marked high-confidence; inherited and equivalent-heading anchors
carry their match basis so callers can tell how direct the anchor is.

## Validate

Run:

```sh
python3 graph/validate_graph.py
```

The validator checks:

- generated JSON/YAML equivalence.
- node and edge counts.
- duplicate node IDs.
- missing edge endpoints or missing `relation`.
- missing source-resolution metadata.
- missing resolved source or manifest target paths.
- future-paper markers on note-only sources.
- malformed TeX anchors.

## Query

Use the query CLI for small context packets while working on a specific physics
thread:

```sh
python3 graph/query_graph.py search 4PN
python3 graph/query_graph.py show OPEN_QUAD_NORMALIZATION
python3 graph/query_graph.py neighborhood CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL --depth 1
python3 graph/query_graph.py paths FILE_4PN_FULL OPEN_QUAD_NORMALIZATION --max-depth 3
python3 graph/query_graph.py source research/4d_4pn/paper/4d_4pn.tex
python3 graph/query_graph.py open-gates 5PN
python3 graph/query_graph.py context moving throat --depth 1
```

The query modes are:

- `search` - find nodes by ID, label, status, source path, type, or metadata.
- `show` - inspect one node and its incoming/outgoing edges.
- `neighborhood` - extract bounded upstream/downstream graph context.
- `paths` - find short dependency paths between two nodes.
- `source` - list graph nodes tied to a TeX paper, note file, or legacy markdown name.
- `open-gates` - list open gates, open-status claims, and future-paper markers near a topic.
- `context` - produce a compact Markdown or JSON packet grouped by claims, equations, open gates, future-paper notes, sources, and backlink blocks.
- `stats` - summarize graph counts by layer, source kind, status, and relation.

Subset-producing commands accept `--write-subgraph path.json`. The resulting
JSON can be rendered with the visualization generator:

```sh
python3 graph/query_graph.py context 4PN --write-subgraph /tmp/atlas_4pn_subset.json
python3 graph/render_graph.py --graph /tmp/atlas_4pn_subset.json --html /tmp/atlas_4pn_subset.html --svg /tmp/atlas_4pn_subset.svg --dot /tmp/atlas_4pn_subset.dot
```

## Visualize

Generate the visualization artifacts:

```sh
python3 graph/render_graph.py
```

This writes:

- `fluid_universe_derivation_atlas_graph.html` - standalone interactive Canvas view.
- `fluid_universe_derivation_atlas_graph.svg` - static dependency overview.
- `fluid_universe_derivation_atlas_graph.dot` - compact Graphviz DOT export.

The HTML file is the primary full-graph view. It has no external dependencies
and can be opened directly in a browser. It supports search, layer filtering,
source-kind filtering, relation filtering, hover inspection, click-to-highlight
neighborhoods, panning, zooming, and TeX-anchor inspection.

The DOT export uses compact node IDs and relation tooltips so static Graphviz
renders stay tractable. Use the HTML view when you need full labels and
metadata inspection.
The SVG output is generated directly by the renderer, so Graphviz is not
required for the default visualization workflow.

## Source Policy

Existing maintained TeX papers are canonical graph targets. Markdown-only
derivation ledgers remain linked deliberately and carry `future_paper_needed`
metadata so they continue to track papers that still need to be written.
