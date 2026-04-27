# Fluid Universe Obsidian Atlas

This directory is the Obsidian vault for browsing the Fluid Universe derivation
atlas.

The vault is read-only for human use. Durable content here should be generated
from the maintained graph source:

```text
graph/fluid_universe_derivation_atlas_graph.yaml
```

The graph YAML, TeX papers, and notes remain the source of truth. Obsidian
Markdown notes, Canvas files, Bases, and exports are presentation/query
artifacts that Codex regenerates and validates.

The public-facing reader site is generated separately from `topics.yaml`. It
turns the raw graph into curated topic pages with DOI-backed paper citations,
paper-pending badges for notes-only material, and a future-paper backfill list.

## Open In Obsidian

Open this folder as the vault:

```text
/var/projects/toy_physics/atlas
```

Core plugins currently enabled include Graph, Backlinks, Canvas, Properties,
Bases, Search, and File Explorer. Obsidian Sync is disabled for this generated
vault; Git remains the durable sync/versioning mechanism.

## Maintenance Flow

When the source graph or TeX papers change:

```sh
python3 graph/annotate_tex_anchors.py
python3 graph/sync_graph_formats.py
python3 graph/validate_graph.py
python3 atlas/scripts/generate_obsidian_atlas.py
python3 atlas/scripts/validate_obsidian_atlas.py
python3 atlas/scripts/build_reader_site.py
```

The generator rewrites `nodes/`, `views/`, `bases/`, `canvas/`, `templates/`,
and `exports/` from the graph YAML. The validator checks generated node
frontmatter, wikilinks, Bases, Canvas files, and export counts.
`exports/generated_manifest.txt` records generator-owned files so stale
generated outputs can be cleaned on later runs without touching hand-authored
vault docs.

The reader-site generator writes `site/`, `site/data/`, and
`exports/future_paper_backfill.yaml`. Open `site/index.html` directly in a
browser to inspect the generated static website. Reader math regression cases
live in `fixtures/reader_math_cases.yaml` and cover graph shorthand that should
render as proper MathJax. The `build_reader_site.py` pipeline runs those
fixtures, regenerates the reader site, and validates the output.

## Reference

The original GPT-Pro implementation notes are preserved in `reference/`.
