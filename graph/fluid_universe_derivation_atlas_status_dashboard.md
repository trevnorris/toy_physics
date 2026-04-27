# Fluid Universe Derivation Atlas v0.8 — Codex Paper-Linking Handoff

## What changed in v0.8

v0.8 corrects the source-side backlink plan based on the current repository reality. The generated atlas originally used summary and note filenames where the maintained TeX papers were not visible. This cleanup resolves the maintained paper targets where they exist and leaves markdown-only derivation ledgers marked as future-paper sources.

This pass turns the v0.6/v0.7 backlink layer into a **Codex-ready paper sweep package**:

```text
Backlink entries: 16
Canonical source-side blocks: 16
Status firewalls carried forward: 16
Primary manifest: fluid_universe_derivation_atlas_paper_insertion_manifest.json
Codex prompt: fluid_universe_derivation_atlas_codex_sweep_prompt.md
Source-resolution register: fluid_universe_derivation_atlas_source_resolution.md
Primary graph source: fluid_universe_derivation_atlas_graph.yaml
Generated graph mirror: fluid_universe_derivation_atlas_graph.json
```

## Corrected completion status

| Atlas pass | Status after v0.8 | Notes |
|---|---|---|
| Curated dependency graph | Complete | v0.1–v0.4 established the graph backbone. |
| Status gates / open problems | Complete | v0.5–v0.7 preserve gates and false-upgrade firewalls. |
| Equation / claim / source anchors | Working complete | Enough for project navigation and future extension. |
| Paper backlink register | Complete as canonical register | v0.6 blocks remain the source of truth. |
| Source resolution to maintained TeX papers | Complete for existing papers | Markdown-only derivation ledgers remain marked as future-paper sources. |
| TeX-local graph anchors | Complete for paper-backed section anchors | Generated `tex_anchor` metadata points to paper lines, headings, and nearby labels. |
| Physical insertion into full paper drafts | Delegated / pending | Backlink blocks are ready, but this cleanup did not insert them into paper prose. |
| Codex handoff for paper sweep | **Complete in v0.8** | This pass supplies manifest, prompt, validation policy, and source-resolution register. |

## Operational rule

The backlink register is now treated as a **bidirectional linking contract**:

```text
atlas node -> resolved source target -> Codex insertion block -> full paper section
```

Codex should update the full papers by inserting compact backlink blocks, but the atlas should remain the canonical source of anchor IDs and status labels.

For graph maintenance, edit the YAML graph and regenerate the JSON mirror with
`python3 graph/sync_graph_formats.py`. If TeX headings or labels move,
regenerate `tex_anchor` metadata with `python3 graph/annotate_tex_anchors.py`
before syncing JSON.

## What Codex should do

1. Load `fluid_universe_derivation_atlas_paper_insertion_manifest.json`.
2. For each entry, find the corresponding full paper draft in the local repository.
3. Locate the recommended insertion point by section title or semantic match, not by stale summary-file line number.
4. Insert the block only if the same `Backlink block: ...` is not already present.
5. Preserve every atlas ID exactly.
6. Preserve every status note exactly or with stronger caution, never weaker caution.
7. Produce an application report with patched files, skipped files, missing full drafts, and conflicts.

## Core status firewalls Codex must preserve

- `FIREWALL_CHARGE_NOT_CIRCULATION` — Do not identify circulation, throat geometry, or magnetic/vortical winding with electric charge sign.
- `FIREWALL_QEFF_THICKNESS_NOT_BREATHING` — Keep microscopic branch charge q_star and observed brane coupling q_eff separate from throat radius/length dynamics.
- `FIREWALL_KAPPARHO_NOT_CHARGE` — Do not recycle the historical gravity-side q=1 notation as electromagnetic charge.
- `FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE` — Treat A_w, J^w, F_{mu w}, E_w, and C_a as suppressed in the brane limit, not deleted.
- `FIREWALL_PROJECTION_NOT_REDUCTION` — Projection is an exact measurement definition; reduction is a controlled dynamical simplification.
- `FIREWALL_PARENT_WALL_NOT_STRICT` — Do not upgrade the moving-wall PDE to strict parent-level status unless S_eta/S_Sigma is included.
- `FIREWALL_WALL_COEFFS_BRANCH_DATA` — Treat wall constitutive coefficients as computed/constrained branch data, not as adjustable rescue parameters.
- `FIREWALL_4PN_LOCAL_NOT_FULL_TAIL` — Keep local instantaneous 4PN closure separate from the inherited passive/outgoing quadrupole normalization gate.
- `FIREWALL_25PN_CONDITIONAL` — Preserve the conditional theorem status of 2.5PN until actual branch normalization is realized.
- `FIREWALL_ANGULAR_NOT_RADIAL` — Do not convert the exact angular identity mhat_ang=1 into a solved radial/axial normalization theorem.
- `FIREWALL_READOUTS_NOT_THROAT` — Keep response packets downstream of the physical throat, mouth, support, mixed, and outgoing branch data.
- `FIREWALL_ATOM_REDUCED_SECTOR` — Label atomic hydrogen as a reduced-sector consequence, not a completed full-defect PDE theorem.
- `FIREWALL_LEPTON_CONDITIONAL` — Keep the lepton same-charge quantizer as a conditional theorem until the autonomous/recirculation closure is file-grounded.
- `FIREWALL_SIMILARITY_NOT_FULL_5PN` — Treat the similarity/orbit theorem as a normalization-defect theorem, not as closure of the separate conservative even gates.
- `FIREWALL_G2_COMMON_CONDITIONAL` — Preserve the common-layer status as a sharp conditional target until the actual moving-throat branch tangent is computed.
- `FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED` — Carry the weighted-gauge-fixing caveat/patch in any derivation that uses localized Maxwell equations beyond the safe zero-mode reading.

## v1.0 interpretation

If the project accepts the backlink register as the canonical source-side backlink layer, then v0.7 was already a v1.0 release candidate. If the project requires the actual papers to contain backlinks, then v0.8 is the handoff that enables Codex to complete that final source-editing pass.

## Deliverables in this pass

- `fluid_universe_derivation_atlas_paper_insertion_manifest.json`
- `fluid_universe_derivation_atlas_codex_sweep_prompt.md`
- `fluid_universe_derivation_atlas_linking_policy.md`
- `fluid_universe_derivation_atlas_graph_patch.md`
- `sync_graph_formats.py`
- `annotate_tex_anchors.py`
- `query_graph.py`
- `render_graph.py`
- `validate_graph.py`
- `fluid_universe_derivation_atlas_status_dashboard.md`
- `fluid_universe_derivation_atlas_validation.md`
- `fluid_universe_derivation_atlas_source_resolution.md`
