# Stage Verification Coverage Baseline

This document is the stage-coverage control sheet for the rebuilt PDE ledger.

```yaml
canonical_stage_count: 13
verified_stage_count: 13
sympy_audit_count: 13
mathematica_audit_count: 13
numerical_stress_count: 0
reviewed_stage_count: 13
```

## Coverage By Part

| Part | Stage Range | Total | SymPy | Mathematica | Numerical | Review |
|---|---|---:|---:|---:|---:|---:|
| I | 004-007 | 4 | 4 | 4 | 0 | 4 |
| II | 001-002, 008-013 | 8 | 8 | 8 | 0 | 8 |
| III | 003 | 1 | 1 | 1 | 0 | 1 |

## Coverage Classes

| Coverage class | Count | Stage ranges |
|---|---:|---|
| SymPy + Mathematica | 13 | 001-013 |
| SymPy only | 0 | none |
| Mathematica only | 0 | none |
| No executable audit | 0 | none |

## Constant Provenance Rule

Coverage counts are not trust grades. Any unexplained literal in a future stage
audit should be treated as a verification defect until its provenance is
classified.
