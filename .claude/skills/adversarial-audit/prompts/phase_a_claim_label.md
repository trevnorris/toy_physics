# Phase A Claim-Label Scan

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `{STAGES}`
- Source files: `{SOURCE_FILES}`

Task:
Find every claim label or prose claim that says or implies a value is derived, forced, exact, non-tunable, matched, canonical, not free, or fixed. The target is the claim, not whether the claim is true.

Emit only YAML:

```yaml
modality: claim_label
candidates:
  - candidate_key:
    anchor_stage:
    parameter_names: []
    citation:
      path:
      line:
      excerpt:
    reason:
```
