# Atlas Agent Rules

This vault is AI-maintained and human-read-only.

- Do not make physics-status changes directly in generated Obsidian notes.
- Update `graph/fluid_universe_derivation_atlas_graph.yaml` first, then regenerate derived artifacts.
- Do not weaken status firewalls or promote open/conditional claims from Obsidian text.
- Keep stable atlas IDs exactly as written.
- Put equations in Markdown bodies, not long YAML frontmatter.
- Treat Canvas and Bases as views, not graph truth.
- Ignore volatile `.obsidian` workspace/cache/plugin state.
- Regenerate generated files with `python3 atlas/scripts/generate_obsidian_atlas.py`.
- Validate generated files with `python3 atlas/scripts/validate_obsidian_atlas.py`.
- Keep public reader prose in `atlas/topics.yaml`; generated files in `atlas/site/` are output.
- Build and validate the reader site with `python3 atlas/scripts/build_reader_site.py`.
- Treat `exports/future_paper_backfill.yaml` as the notes-to-paper citation cleanup queue.

Generated files must keep the generated-file notice and should not be edited by
hand.
