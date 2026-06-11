# Phase A Graph Scan

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160`
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
