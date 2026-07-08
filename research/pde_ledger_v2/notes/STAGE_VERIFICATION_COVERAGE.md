# Stage Verification Coverage Baseline

This document is the stage-coverage control sheet for the rebuilt PDE ledger.

```yaml
canonical_stage_count: 6
verified_stage_count: 6
sympy_audit_count: 6
mathematica_audit_count: 6
numerical_stress_count: 0
reviewed_stage_count: 6
```

## Coverage By Part

| Part | Stage Range | Total | SymPy | Mathematica | Numerical | Review |
|---|---|---:|---:|---:|---:|---:|
| I | 004-006 | 3 | 3 | 3 | 0 | 3 |
| II | 001-002 | 2 | 2 | 2 | 0 | 2 |
| III | 003 | 1 | 1 | 1 | 0 | 1 |

## Coverage Classes

| Coverage class | Count | Stage ranges |
|---|---:|---|
| SymPy + Mathematica | 6 | 001-006 |
| SymPy only | 0 | none |
| Mathematica only | 0 | none |
| No executable audit | 0 | none |

## Constant Provenance Rule

Coverage counts are not trust grades. Any unexplained literal in a future stage
audit should be treated as a verification defect until its provenance is
classified.
