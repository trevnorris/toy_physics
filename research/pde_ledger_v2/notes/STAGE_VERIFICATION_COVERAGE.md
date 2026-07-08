# Stage Verification Coverage Baseline

This document is the stage-coverage control sheet for the rebuilt PDE ledger.

```yaml
canonical_stage_count: 11
verified_stage_count: 11
sympy_audit_count: 11
mathematica_audit_count: 11
numerical_stress_count: 0
reviewed_stage_count: 11
```

## Coverage By Part

| Part | Stage Range | Total | SymPy | Mathematica | Numerical | Review |
|---|---|---:|---:|---:|---:|---:|
| I | 004-007 | 4 | 4 | 4 | 0 | 4 |
| II | 001-002, 008-011 | 6 | 6 | 6 | 0 | 6 |
| III | 003 | 1 | 1 | 1 | 0 | 1 |

## Coverage Classes

| Coverage class | Count | Stage ranges |
|---|---:|---|
| SymPy + Mathematica | 11 | 001-011 |
| SymPy only | 0 | none |
| Mathematica only | 0 | none |
| No executable audit | 0 | none |

## Constant Provenance Rule

Coverage counts are not trust grades. Any unexplained literal in a future stage
audit should be treated as a verification defect until its provenance is
classified.
