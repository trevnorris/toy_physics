# Ablation-driver v4 conformance fixtures

This suite grades `REQUIREMENTS.md` A1--A9 against the public executable fixed by
`CONTRACT.md`. It contains no driver. `fixture_oracle.py` is an adoptable
derivation library spanning all items, but it is not a complete driver: it has no
child execution, recovery, resume, or orchestration. The separation rests on
independent authorship and the freeze, not on secrecy.

The byte authority is
`research/pde_ledger_v2/notes/ablation_driver/fixtures_v4.accepted.sha256`.
The authority is intentionally outside this directory. The verifier requires its
committed `HEAD` bytes and a clean suite/authority worktree before checking the
complete inventory. Until the user commits the authority and suite, freeze status
is PENDING-COMMIT.

`run_conformance.py` runs that verifier first and refuses conformance when the
freeze is not verified. `FREEZE_LIMITS.md` states the protection this provides
and both the byte-authority and execution-environment boundaries it does not
cover.

From the project root, verify the committed freeze:

```text
PYTHONDONTWRITEBYTECODE=1 python3 research/pde_ledger_v2/notes/ablation_driver/fixtures_v4/verify_freeze.py
```

Run the full ladder:

```text
PYTHONDONTWRITEBYTECODE=1 python3 research/pde_ledger_v2/notes/ablation_driver/fixtures_v4/run_conformance.py
```

For diagnosis, append one or more names from `A1` through `A9`. Exit zero means
all selected items passed, exit one means an item failed, and exit two means the
contract-path driver is absent. Exit 124 is always a failed run.

A7 is one uninterrupted invocation in a disposable non-Git mirror. It reports
each fully audited, durably captured row as progress and uses a per-row inactivity
guard rather than an aggregate wall-clock deadline. Successful completion is
audited through its final `run` restoration event; `repair` is exercised only
after the uncatchable-kill boundary in A4.

Every item applies the C-4/C-5 result and capture audit to completed evidence. A6
also submits a contract-invalid list and requires the validation failure path, in
addition to checking both the committed and runtime-generated valid lists.

`NEGATIVE_DEMONSTRATIONS.md` records the pre-driver able-to-fail execution. The
demonstration used disposable probes; no probe, expected table, or driver stub is
retained.

## Known and accepted coverage gaps

This gate does not attempt exhaustive failure-code coverage: it exercises success
and a catchable-signal driver exit, while the uncatchable-kill case proceeds to
repair. Its invalid-list floor is not exhaustive; other invalid-input grammar,
timeout or child-start failure, and stable-input drift cases are not covered.

Multi-target runs and restoration, multiple checkers or artifacts and their
combined behavior, and the truth of every restore event's pre-restore claim are
not covered. The direct-argv probe distinguishes a quoting shell pass from direct
execution, and the child requires exact configured-environment membership so a
selectively leaked ambient whitelist also fails.

Performance, resistance to hardcoding against the public oracles, and external
trust or signing are also outside this gate. The caller-tree guard covers `.git`,
but `FREEZE_LIMITS.md` explains why that does not authenticate the environment or
repository-local configuration used by the earlier freeze check.
