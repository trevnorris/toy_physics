# Stage Verification Coverage Baseline

This document is the stage-coverage control sheet for the rebuilt PDE ledger.

```yaml
canonical_stage_count: 32
verified_stage_count: 31
sympy_audit_count: 31
mathematica_audit_count: 31
numerical_stress_count: 0
reviewed_stage_count: 32
```

> **Stage 029 (PN corpus DOI-cite) is CITE-only** — a documentary provenance stage with **no executable audit**
> (no `scripts/`/`mathematica/` pair). It increments `canonical_stage_count` and `reviewed_stage_count` (it gets the
> adapted doc-fidelity review) but NOT `sympy_audit_count`/`mathematica_audit_count`/`verified_stage_count`
> (`verified` counts executable-verified stages only). **Note (028-lag reconciliation):** at Stage 028 only
> `canonical`/`reviewed` had been bumped (27→28) while the audit trio was left at 27, even though Stage 028 is
> dual-engine; the By-Part and Coverage-Classes tables already counted 028. The trio is corrected to 28 here (Stage
> 028's two audits) alongside the Stage 029 bumps, so the header counts, the tables, and the manifest agree.

## Coverage By Part

| Part | Stage Range | Total | SymPy | Mathematica | Numerical | Review |
|---|---|---:|---:|---:|---:|---:|
| I | 004-007 | 4 | 4 | 4 | 0 | 4 |
| II | 001-002, 008-029 | 24 | 23 | 23 | 0 | 24 |
| III | 003 | 1 | 1 | 1 | 0 | 1 |
| IV | 030-032 | 3 | 3 | 3 | 0 | 3 |

## Coverage Classes

| Coverage class | Count | Stage ranges |
|---|---:|---|
| SymPy + Mathematica | 31 | 001-028, 030-032 |
| SymPy only | 0 | none |
| Mathematica only | 0 | none |
| No executable audit | 1 | 029 |

## Constant Provenance Rule

Coverage counts are not trust grades. Any unexplained literal in a future stage
audit should be treated as a verification defect until its provenance is
classified.
