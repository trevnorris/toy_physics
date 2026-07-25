# Composite manifest report

Headline: **FAIL**

Coverage: `citations=17/20 claims=106/106 unresolved_producers=3 c7_edges=1/17 closure=stage030,stage031,stage032 mathematica_outputs_checked=0`

Dimensional consistency scope: Verifies manifest-internal dimensional algebra: relation homogeneity, declared-vs-recovered agreement for symbols whose dimensions are recoverable from the stage's audit script, and cross-stage agreement of shared quantities. It does NOT independently certify that a stage's dimensions are physically correct — that is owned by the stage's dual-engine unit audit.

Declared dimensional bases (legacy omission is reported as the implicit `LMT [L, M, T]` basis):

| Stage | Basis |
|---|---|
| stage030 | `LMT [L, M, T]` |
| stage031 | `LMT [L, M, T]` |
| stage032 | `LMT [L, M, T]` |

| Check | Outcome |
|---|---|
| SCHEMA | PASS |
| EVIDENCE | PASS |
| ADJUDICATION | PASS |
| IMPORT | PASS |
| C1 | PASS |
| C2 | PASS |
| C3 | PARTIAL |
| DIMENSIONAL_CONSISTENCY | FAIL |
| C5 | PARTIAL |
| C6 | PASS |
| C7 | PARTIAL |
| C8 | PASS |

## SCHEMA

No findings.

## EVIDENCE

No findings.

## ADJUDICATION

No findings.

## IMPORT

No findings.

## C1

No findings.

## C2

No findings.

## C3

- `ABSENT_PRODUCER` — stage030:stage003/b_eff_definition
- `ABSENT_PRODUCER` — stage030:stage003/brane_inertia
- `ABSENT_PRODUCER` — stage031:stage003/b_eff_definition

## DIMENSIONAL_CONSISTENCY

Verifies manifest-internal dimensional algebra: relation homogeneity, declared-vs-recovered agreement for symbols whose dimensions are recoverable from the stage's audit script, and cross-stage agreement of shared quantities. It does NOT independently certify that a stage's dimensions are physically correct — that is owned by the stage's dual-engine unit audit.

- `DIM_SOURCE_NOT_REGISTERED` — stage032:scripts/ledger_stage031_puncture_deflection_field_identity_source_sympy_audit.py

## C5

- `REGISTER_COVERAGE_PARTIAL` — 44 knob/GAP register rows await classification or their owning stages
- `LIFECYCLE_PRODUCER_UNEXTRACTED` — stage030:k.brane.effective_inertia:stage003
- `LIFECYCLE_PRODUCER_UNEXTRACTED` — stage030:k.brane.effective_stiffness:stage003
- `RANGE_CLAIM_UNRESOLVED` — stage043/count_range

## C6

No findings.

## C7

- `C7_EDGE_UNCOVERED` — 1/17 edges have executable C7 metadata

## C8

No findings.
