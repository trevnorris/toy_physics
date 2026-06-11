# Phase A Graph Scan

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224`
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
