# Codex sweep prompt — Fluid Universe Atlas paper backlinks v0.8

You are working in the local Fluid Universe repository. The generated atlas and paper-backlink register now include a source-resolution layer: existing maintained TeX papers are canonical targets, while markdown-only derivation ledgers are explicit future-paper backlog. Your job is to apply the backlink layer to maintained papers where they exist and report note-only entries without inventing paper targets.

## Inputs

Use these files from the graph package:

- `fluid_universe_derivation_atlas_paper_insertion_manifest.json`
- `fluid_universe_derivation_atlas_paper_backlink_blocks.md` if present
- `fluid_universe_derivation_atlas_status_firewall_register.md` if present
- `fluid_universe_derivation_atlas_source_resolution.md` if present
- `fluid_universe_derivation_atlas_graph.yaml` for `tex_anchor` hints

## Required behavior

For each manifest entry:

1. Search the repository for the candidate full draft paths or globs listed in `candidate_full_draft_paths_or_globs`.
2. Prefer the maintained full paper over a summary file when both exist.
3. If `future_paper_needed` is true, do not invent a TeX target; keep the markdown note as the source and report it as paper backlog.
4. Find the `recommended_insertion_point` using section headings and nearby content. Use graph-node `tex_anchor` fields as hints when present. Do not rely on summary-file line numbers as exact positions inside full drafts.
5. Insert the `block_markdown` as a compact subsection near the insertion point.
6. If the exact `Backlink block: <ID>` already exists, do not duplicate it. Instead validate that anchors and status notes are current.
7. Preserve all atlas anchor IDs exactly.
8. Preserve status notes exactly unless making them more conservative.
9. Never use backlinks to upgrade claim status.
10. Write a final report named `codex_atlas_backlink_application_report.md` with:
   - patched files,
   - already-present blocks,
   - missing target files,
   - ambiguous insertion points,
   - conflicts or needed human review,
   - any changed atlas version tags.

## Status firewalls that must not be violated

- `FIREWALL_CHARGE_NOT_CIRCULATION`: Do not identify circulation, throat geometry, or magnetic/vortical winding with electric charge sign.
- `FIREWALL_QEFF_THICKNESS_NOT_BREATHING`: Keep microscopic branch charge q_star and observed brane coupling q_eff separate from throat radius/length dynamics.
- `FIREWALL_KAPPARHO_NOT_CHARGE`: Do not recycle the historical gravity-side q=1 notation as electromagnetic charge.
- `FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE`: Treat A_w, J^w, F_{mu w}, E_w, and C_a as suppressed in the brane limit, not deleted.
- `FIREWALL_PROJECTION_NOT_REDUCTION`: Projection is an exact measurement definition; reduction is a controlled dynamical simplification.
- `FIREWALL_PARENT_WALL_NOT_STRICT`: Do not upgrade the moving-wall PDE to strict parent-level status unless S_eta/S_Sigma is included.
- `FIREWALL_WALL_COEFFS_BRANCH_DATA`: Treat wall constitutive coefficients as computed/constrained branch data, not as adjustable rescue parameters.
- `FIREWALL_4PN_LOCAL_NOT_FULL_TAIL`: Keep local instantaneous 4PN closure separate from the inherited passive/outgoing quadrupole normalization gate.
- `FIREWALL_25PN_CONDITIONAL`: Preserve the conditional theorem status of 2.5PN until actual branch normalization is realized.
- `FIREWALL_ANGULAR_NOT_RADIAL`: Do not convert the exact angular identity mhat_ang=1 into a solved radial/axial normalization theorem.
- `FIREWALL_READOUTS_NOT_THROAT`: Keep response packets downstream of the physical throat, mouth, support, mixed, and outgoing branch data.
- `FIREWALL_ATOM_REDUCED_SECTOR`: Label atomic hydrogen as a reduced-sector consequence, not a completed full-defect PDE theorem.
- `FIREWALL_LEPTON_CONDITIONAL`: Keep the lepton same-charge quantizer as a conditional theorem until the autonomous/recirculation closure is file-grounded.
- `FIREWALL_SIMILARITY_NOT_FULL_5PN`: Treat the similarity/orbit theorem as a normalization-defect theorem, not as closure of the separate conservative even gates.
- `FIREWALL_G2_COMMON_CONDITIONAL`: Preserve the common-layer status as a sharp conditional target until the actual moving-throat branch tangent is computed.
- `FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED`: Carry the weighted-gauge-fixing caveat/patch in any derivation that uses localized Maxwell equations beyond the safe zero-mode reading.

## Special caution notes

- The 4D parent action/projection/continuity identities can be exact, but Poisson and zero-mode Maxwell outputs remain controlled reductions.
- The localized Maxwell sector must not erase the mixed core; `A_w`, `J^w`, `F_{\mu w}`, `E_w`, and `C_a` remain microscopic channels outside the zero-mode brane limit.
- The moving-wall PDE must not be described as strict parent-level unless `S_eta` or `S_Sigma` is explicitly included in the total parent action.
- The 2.5PN and 4PN quadrupole/tail results remain conditional on the same passive/outgoing quadrupole normalization gate.
- The atom, lepton, anomaly, and 5PN branches must keep their reduced/conditional/open statuses.

## Output format

After applying patches, produce:

```markdown
# Codex Atlas Backlink Application Report

## Summary
- total manifest entries:
- patched:
- already present:
- skipped / missing files:
- ambiguous / human review:

## Patched files
| file | backlink block | insertion heading | notes |

## Already present
| file | backlink block | validation result |

## Missing or ambiguous
| manifest entry | candidate paths checked | issue | recommended action |

## Status-firewall check
| firewall | pass/fail | notes |
```
