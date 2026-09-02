# Derived Registry Review Record

Date: 2026-08-24
Superseded for current corpus selection: 2026-09-01 by
[`pde_ledger_exclusion_validation.md`](pde_ledger_exclusion_validation.md).
The historical review below remains as recorded.
Artifact: `output_budgets`, `derived_pages`, mixed-page bases, and Atlas routes
Method: two independent read-only review legs, followed by one fold

## Convergent blockers

Both legs found that unit-wide derived inputs were technically below an 8 MB
ceiling but operationally far too broad, reaching 5.9–7.3 MB and roughly
one-to-two million tokens. They also found missing Atlas source coverage,
undeclared benchmark destinations, incomplete conflict dependencies, missing
actual PDE-audit engines, advisory rather than enforced per-entry budgets, and
no deterministic migration-completeness contract.

The operational leg additionally required full dependency freshness metadata,
derived-orphan handling, and a citation allowlist narrowed to actually readable
direct inputs.

## Fold

The fold:

- added narrow committed support records from PDE audit, the legacy/v2 lineage,
  and Stage-1;
- created dedicated v3-status, Stage-1, emergent-charge, and
  emergent-magnetism topics;
- established the benchmark's one-time topic-name crosswalk;
- made all 15 topics inputs to the conflict register;
- selected representative Python/Wolfram audit engines;
- added deterministic `direct_sources` projections with literal-verified
  excerpt anchors;
- lowered the derived semantic packet cap to 750,000 bytes;
- added and enforced maximum entry counts and per-entry word limits;
- normalized broad Atlas directory references to exact committed files;
- hashed and routed migration requirements by exact derived target.

Validation after projection resolved 33 source units and 22 derived pages. All
21 semantic derived tasks were below 147 KB; the index has no direct source
payload. Every selected excerpt anchor existed literally in the committed HEAD
blob, and every direct source belonged to a declared input unit.
