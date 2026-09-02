# Config Resolution Validation

Date: 2026-08-24
Superseded for current corpus selection: 2026-09-01 by
[`pde_ledger_exclusion_validation.md`](pde_ledger_exclusion_validation.md).
The historical validation below remains as recorded.
Config: `memory/_meta/config.yaml` schema version 2
Committed snapshot: `ddecdbc249a4d9d860dd25ecd4531726a01d26c6`

## Purpose

This read-only check independently resolved every configured unit against one
NUL-delimited enumeration of the committed `research/` and `software/` tree.
It applied the selector contract locally, checked required exact members and
selectors, enforced unique unit IDs and capsule targets, rejected duplicate
resolved members and primary ownership overlaps, and checked semantic member
and unit byte limits. It did not read the mutable index or working-tree source
bytes.

Git enumeration:

```bash
git ls-tree -rz --long --full-tree HEAD -- research software
```

Literal summary:

```text
HEAD ddecdbc249a4d9d860dd25ecd4531726a01d26c6
TREE_BLOBS 7840
UNITS 32
UNIQUE_RESOLVED_PATHS 168
PRIMARY_OWNED_PATHS 168
SEMANTIC_BYTES 9762954
ERRORS 0
```

The 168 resolved paths include semantic and identity-only members. The
9,762,954-byte figure is the sum of semantic members across all units, not a
single transaction read. The largest configured unit is
`pde-ledger-v3-s10` at 3,657,813 semantic bytes, below the configured
8,000,000-byte per-unit limit. Every individual semantic member is below the
configured 2,000,000-byte member limit.

## Result

The normalized Phase-0 configuration resolves without missing required
members, empty required selectors, duplicate targets, ownership conflicts, or
read-limit violations at the recorded snapshot. This is a configuration
fixture result, not a substitute for the production scanner's own validation
and tests.
