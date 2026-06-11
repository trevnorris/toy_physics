# Phase A Graph Scan

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128`
- Graph wrapper command: `timeout 600 python3 graph/query_graph.py --graph graph/fluid_universe_derivation_atlas_graph.yaml`

Task:
Use only the atlas graph and graph wrapper outputs to identify nodes carrying numerical commitments, parameter fixes, external-match claims, or status claims that may be fit-insertion points. If a stage source has no graph node, emit a `graph_gap` entry instead of inventing graph evidence.

Emit only YAML:

```yaml
modality: graph
candidates:
  - candidate_key:
    anchor_stage:
    parameter_names: []
    graph_node:
    citation:
      path:
      line:
      excerpt:
    reason:
graph_gaps:
  - stage:
    source:
    reason:
```
