# Stage Provenance Index

This repo-local index tracks raw derivation-note sources and executable audit
artifacts for canonical stage files after they are added to the rebuilt ledger.

Current scope: Stages 001.

| Stage | Canonical TeX | Note Source | SymPy Audit | Mathematica Audit | Numerical Stress |
|---|---|---|---|---|---|
| 001 | paper/stages/stage_001.tex | notes/stages/ledger_stage001_solid_angle_second_moment_primitives.md | scripts/ledger_stage001_solid_angle_second_moment_primitives_sympy_audit.py (+ scripts/output/*.txt) | mathematica/ledger_stage001_solid_angle_second_moment_primitives_mathematica_audit.wl (+ mathematica/output/*.txt) | — |
