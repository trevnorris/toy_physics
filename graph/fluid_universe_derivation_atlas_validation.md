# Fluid Universe Derivation Atlas v0.8 — Validation Report

## Generated artifacts

```text
paper insertion manifest: yes
codex sweep prompt: yes
linking policy: yes
graph patch: yes
graph v0.8 YAML primary: yes
generated JSON mirror: yes
source-resolution register: yes
TeX-anchor annotations: yes
zip packet: no
```

## Manifest statistics

```text
entries: 16
entries with open gates: 12
entries with primary anchors: 16
entries with source-section anchors: 16
```

## Non-claim-upgrade check

v0.8 adds only repository-operation metadata and source-resolution metadata. It does not modify any physics claim status, open gate, equation anchor, or theorem status.

## Format and anchor check

The maintained graph source is YAML. JSON is generated from YAML with
`python3 graph/sync_graph_formats.py`, and validation fails if the mirror is
stale. TeX anchors are generated with `python3 graph/annotate_tex_anchors.py`;
they point graph nodes to paper-local lines, headings, and nearby `\label{...}`
entries without editing the TeX papers.

## Dangling operational nodes check

The v0.8 operation nodes are connected by new operation edges. Register artifacts `BACKLINK_REGISTER_V06` and `STATUS_FIREWALL_REGISTER_V07` are represented as atlas-meta nodes, so the v0.8 operation edges have no intentional dangling targets.

## Final graph size

```text
nodes: 373
edges: 1325
```
