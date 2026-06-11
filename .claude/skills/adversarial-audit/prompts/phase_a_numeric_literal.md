# Phase A Numeric-Literal Scan

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `{STAGES}`
- Source files: `{SOURCE_FILES}`

Task:
Find numerical literals, rational coefficients, closed-form constants, or parameter assignments that could be fit-insertion-point candidates. Exclude pure identities such as harmonic normalizations, dimension labels, residual tolerances, line numbers, and pass/fail counters unless the source itself claims the value is matched, fixed, canonical, forced, or derived.

Emit only YAML:

```yaml
modality: numeric_literal
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
