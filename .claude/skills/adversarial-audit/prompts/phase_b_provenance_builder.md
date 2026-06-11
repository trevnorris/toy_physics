# Phase B Provenance Builder

Build a parameter-value provenance slice for candidate `{CANDIDATE_ID}`.

Binding rules:
- No attribution without opening the per-stage `notes/` source. Do not commit an "introduced at stage N" claim from cards, scripts, outputs, or graph nodes alone.
- A self-contradictory origin claim is itself a finding, never silently resolved.
- If the atlas graph has no node for a candidate source or dependency, log `graph_gap` with the source path and continue with cited source evidence.

Taxonomy:
- File/source provenance tracks where text or scripts live.
- Result/stage provenance tracks which stage established a result.
- Parameter-value provenance tracks where a numerical value or parameter assignment was introduced and what constrained it.

This ledger is parameter-value provenance only. Do not conflate it with file lineage or stage attribution alone.

Inputs:

```yaml
{CANDIDATE_YAML}
```

Emit only YAML with `parameter_name`, `source_evidence`, `origin_claims`, `constraints`, `graph_context`, and `provenance_findings`.
