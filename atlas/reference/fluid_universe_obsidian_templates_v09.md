# Fluid Universe Obsidian Templates v0.9

## Generic node template

```markdown
---
id: NODE_ID
title: Human-readable title
type: claim
layer: derivation
status: exact_within_closure
claim_strength: exact_algebra_within_declared_hierarchy
atlas_version: v0.9
source_graph_version: v0.8
generated_by: codex
last_generated_utc: 2026-04-25T00:00:00Z
summary_short: One-sentence summary.
source_files: []
source_sections: []
physical_ids: []
math_ids: []
equation_ids: []
claim_ids: []
open_gate_ids: []
status_firewall_ids: []
outgoing_edges: []
incoming_edges: []
tags:
  - atlas/node
---

# Human-readable title

> **Status:** exact_within_closure  
> **Layer:** derivation  
> **Type:** claim

## Summary

A readable summary goes here. Keep this useful on its own.

## Physical meaning

Explain what this means physically in the Fluid Universe ontology.

## Mathematical role

Explain the formal object, reduction, equation, or derivation role.

## Key equations

$$
...
$$

## Atlas links

- Physical: [[PHYS_...]]
- Math: [[MATH_...]]
- Equations: [[EQ_...]]
- Open gates: [[OPEN_...]]
- Status firewalls: [[FIREWALL_...]]

## Source anchors

- `source_file.md`

## AI maintenance notes

Explain what must be updated if this node changes.
```

## Claim node template

```markdown
---
id: CLAIM_EXAMPLE
title: Example claim
type: claim
layer: derivation
status: controlled_reduction
claim_strength: controlled_reduction
atlas_version: v0.9
source_graph_version: v0.8
generated_by: codex
last_generated_utc: 2026-04-25T00:00:00Z
summary_short: Example claim summary.
source_files: []
source_sections: []
physical_ids: []
math_ids: []
equation_ids: []
claim_ids: []
open_gate_ids: []
status_firewall_ids: []
outgoing_edges: []
incoming_edges: []
tags:
  - atlas/claim
---

# Example claim

## Claim

State what is claimed.

## What it does not claim

State non-claims explicitly.

## Physical meaning

State the physical interpretation.

## Mathematical route

State the derivation/reduction path.

## Status discipline

Explain why the status label is correct.

## Links

- Depends on: [[...]]
- Feeds: [[...]]
- Blocked by: [[...]]
```

## Physical node template

```markdown
---
id: PHYS_EXAMPLE
title: Example physical object
type: physical_object
layer: physical_ontology
status: physical_ontology
atlas_version: v0.9
source_graph_version: v0.8
generated_by: codex
last_generated_utc: 2026-04-25T00:00:00Z
summary_short: Example physical object summary.
source_files: []
source_sections: []
physical_ids: []
math_ids: []
equation_ids: []
claim_ids: []
open_gate_ids: []
status_firewall_ids: []
outgoing_edges: []
incoming_edges: []
tags:
  - atlas/physical
---

# Example physical object

## Physical picture

Describe the object in ordinary physical language.

## Mathematical representation

List the symbols, fields, modes, or reductions representing it.

## What a brane observer sees

Separate projected brane observables from hidden bulk/mixed structure.

## Downstream uses

List major claims, equations, and open gates that use this object.
```

## Equation node template

```markdown
---
id: EQ_EXAMPLE
title: Example equation
type: equation
layer: math_object
status: exact
atlas_version: v0.9
source_graph_version: v0.8
generated_by: codex
last_generated_utc: 2026-04-25T00:00:00Z
summary_short: Example equation summary.
source_files: []
source_sections: []
physical_ids: []
math_ids: []
equation_ids: []
claim_ids: []
open_gate_ids: []
status_firewall_ids: []
outgoing_edges: []
incoming_edges: []
tags:
  - atlas/equation
---

# Example equation

## Equation

$$
...
$$

## Variable dictionary

- $x$: meaning.
- $y$: meaning.

## Derivation role

Explain where the equation enters the derivation graph.

## Status

Explain exact / reduced / closure status.
```

## Open gate template

```markdown
---
id: OPEN_EXAMPLE
title: Example open gate
type: open_gate
layer: open_gate
status: open
atlas_version: v0.9
source_graph_version: v0.8
generated_by: codex
last_generated_utc: 2026-04-25T00:00:00Z
summary_short: Example open theorem gate.
source_files: []
source_sections: []
physical_ids: []
math_ids: []
equation_ids: []
claim_ids: []
open_gate_ids: []
status_firewall_ids: []
outgoing_edges: []
incoming_edges: []
tags:
  - atlas/open_gate
---

# Example open gate

## What remains open

State the unresolved theorem/data/normalization problem.

## What would close it

Define the required branch data, proof, numerical result, or parent-action upgrade.

## Downstream consequences

List claims affected by closing or failing this gate.

## Invalid shortcuts

List status firewalls that protect this gate.
```

## Status firewall template

```markdown
---
id: FIREWALL_EXAMPLE
title: Example status firewall
type: status_firewall
layer: status_audit
status: active_firewall
atlas_version: v0.9
source_graph_version: v0.8
generated_by: codex
last_generated_utc: 2026-04-25T00:00:00Z
summary_short: Prevents an invalid inference.
source_files: []
source_sections: []
physical_ids: []
math_ids: []
equation_ids: []
claim_ids: []
open_gate_ids: []
status_firewall_ids: []
outgoing_edges: []
incoming_edges: []
tags:
  - atlas/status_firewall
---

# Example status firewall

## Invalid inference

State the false shortcut.

## Corrected inference

State the permitted reading.

## Why it matters

Explain what would be overclaimed if this firewall failed.

## Affected nodes

- [[...]]
```

## Dataview dashboard snippets

### Open gates

```dataview
TABLE status, summary_short, source_files
FROM "atlas/nodes/open_gates"
SORT status ASC, file.name ASC
```

### Claims by status

```dataview
TABLE layer, status, claim_strength, summary_short, source_files
FROM "atlas/nodes/claims"
SORT status ASC, layer ASC, file.name ASC
```

### Status firewalls

```dataview
TABLE status, summary_short, source_files
FROM "atlas/nodes/status_firewalls"
SORT file.name ASC
```

### Moving-throat PDE

```dataview
TABLE type, layer, status, summary_short, open_gate_ids
FROM "atlas/nodes"
WHERE contains(tags, "topic/moving_throat")
SORT layer ASC, status ASC, file.name ASC
```
