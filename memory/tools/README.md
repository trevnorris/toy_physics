# Research-memory deterministic tool

`memory.py` is the file/transaction layer for the research memory. It does not
summarize sources or call a model. Its normal boundary is the committed Git
tree under the candidate roots declared in `_meta/config.yaml`.

Run it from anywhere inside the repository:

```bash
python3 memory/tools/memory.py init
python3 memory/tools/memory.py status
python3 memory/tools/memory.py update --prepare --units paper-1pn-hybrid
python3 memory/tools/run_isolated.py <transaction-id> --task paper-1pn-hybrid --profile codex
python3 memory/tools/review_isolated.py <transaction-id> --task paper-1pn-hybrid
python3 memory/tools/memory.py lint --staged <transaction-id> --record
# Repair an exact recorded machine-lint failure after the prior Grok PASS:
python3 memory/tools/run_isolated.py <new-transaction-id> --task paper-1pn-hybrid --profile codex \
  --revise-lint-from <lint-failed-transaction-id>
# After an incomplete reviewer infrastructure/network failure only:
python3 memory/tools/review_isolated.py <transaction-id> --task paper-1pn-hybrid --retry-incomplete
python3 memory/tools/memory.py lint --staged <transaction-id>
python3 memory/tools/memory.py update --finalize <transaction-id>
# Retry from an attested Codex candidate after policy/editorial changes at the same commit:
python3 memory/tools/run_isolated.py <new-transaction-id> --task paper-1pn-hybrid --profile codex \
  --revise-from <prior-transaction-id>
# Model-free reuse after an attested Grok PASS review:
python3 memory/tools/run_isolated.py <new-transaction-id> --task paper-1pn-hybrid \
  --reuse-reviewed-from <reviewed-transaction-id>
# Model-free reuse of an original Codex candidate that needs a fresh review:
python3 memory/tools/run_isolated.py <new-transaction-id> --task paper-1pn-hybrid \
  --reuse-candidate-from <candidate-transaction-id>
python3 memory/tools/review_isolated.py <new-transaction-id> --task paper-1pn-hybrid
python3 memory/tools/memory.py lint --served
python3 memory/tools/memory.py query where is the latest treatment of shear
```

## Transaction storage lifecycle

Transactions contain the complete evidence needed for revisions and reviewed
reuse, so cleanup is deliberately lossless and dry-run-first. Inspect the
current protection graph and proposed operations with:

```bash
python3 memory/tools/transaction_gc.py status
python3 memory/tools/transaction_gc.py archive <transaction-id>
python3 memory/tools/transaction_gc.py prune <transaction-id>
python3 memory/tools/transaction_gc.py restore <transaction-id>
```

Repeat a command with `--apply` only after reviewing its dry run.
`archive --apply` writes and fully verifies a compressed archive plus hash receipt but
never removes the source. `prune --apply` can remove that source only when the
archive still reproduces it exactly and the transaction is not the current
state generation, a current served-page producer, or part of an incoming
revision/candidate-reuse/reviewed-reuse chain. Restore also verifies the complete tree before
publication back under `memory/transactions/`.

Immediately after a successful finalize, while it is still the state's named
generation, record durable proof with:

```bash
python3 memory/tools/transaction_gc.py record-finalized <transaction-id>
python3 memory/tools/transaction_gc.py record-finalized <transaction-id> --apply
```

That receipt lets the tool later distinguish a finalized-but-superseded
transaction from an unclassified current-target work item. Transactions whose
target is no longer `HEAD`, explicitly failed transactions, and finalized
transactions fully superseded in state become prune-eligible only after a
verified archive exists. Archives and lifecycle receipts live in the already
gitignored `memory/transactions/` tree as `_archive-<transaction-id>.*`
artifacts; they are files, not directories that older status commands could
mistake for live transactions.

`update --prepare` prints a directory under `memory/transactions/`. Its
`manifest.json` contains explicit `writer_tasks` and `task_order` lists. Each
task has a sealed, self-contained packet. Initial packet seals are anchored in
the immutable transaction seal. The isolation runner supports only the named
`codex` writer profile, clears the inherited environment, mounts exact runtime
and credential files rather than a home directory, hides the live
workspace with bubblewrap, provides one persistent `/output` directory,
rejects undeclared outputs, and imports the declared page into staging.
Codex authoring has a hard 20-minute wall-clock limit. A timeout preserves any
partial stdout/stderr as `failed-stdout.bin` and `failed-stderr.bin` in the task
output directory, but never imports a candidate into staging.
Finalize requires the matching versioned-runtime isolation attestation.
Identity-only members have identities in the packet but no staged contents.
Normal Codex publication also requires an exact current Grok `PASS`; finalize
rechecks the writer packet/attestation, candidate, review packet/contract/model,
report, reviewer implementation, and review attestation. The normal order is
therefore Codex authoring → Grok review → staged lint → finalize.
Grok review has a hard 15-minute wall-clock limit. A timeout preserves the same
two diagnostic artifacts and leaves an incomplete review eligible for the
normal `--retry-incomplete` archive-and-retry path.
`--retry-incomplete` refuses any attempt containing a recognizable `PASS` or
`FAIL` attestation, report, or structured verdict. It atomically preserves an
eligible failed attempt under the transaction's `rejected-reviews/` directory
before creating a clean active review; finalization considers only `reviews/`.

`--revise-from` never reads a served/live page. It admits only the same task's
staged page from another transaction after verifying a Codex isolation
attestation, output mapping/source-member identities, and the candidate,
packet, and attestation hashes. A source-capsule revision may cross unrelated
commits only when its task/output/source-kind/member identities and unit digest
remain exact; both commits stay recorded and machine frontmatter is normalized
to the current sealed packet. Any changed source blob or unit contract rejects
the revision. Derived revisions remain target-commit-strict because their
prerequisite snapshots are commit-bound. The candidate and its
verification record become read-only, append-only inputs in the new sealed
packet. Prior and current unit digests may differ so current `editorial_focus`
or frozen-policy corrections can be applied only for same-commit revisions;
the revised output receives a normal new Codex writer attestation and must earn
a fresh current Grok `PASS`. A revision also requires the prior
hash-complete Grok `FAIL` review; its report and attestation are sealed into the
retry packet. The exact `--task` is the retry selector, and derived retries bind
their direct-source/dependency selector plus prior/current prerequisite hashes.
A failed candidate created by `--reuse-candidate-from` is also eligible: before
admission the launcher independently re-verifies its complete original
Codex/bubblewrap provenance chain, then verifies the candidate-reuse
transaction's Grok `FAIL` is hash-bound to that exact candidate and writer
attestation. Missing or tampered provenance is rejected. Output from
`--reuse-reviewed-from` is never an eligible revision input.

`--reuse-reviewed-from` is source-unit-only and never invokes a model. It
requires the same task/kind/output/member identities, a hash-bound prior Codex
bubblewrap attestation, and a hash-bound Grok `PASS` review. The repository
commit may advance for unrelated changes only when the complete source-unit
digest also remains identical; changed source blobs or unit contracts reject
cross-commit reuse. It preserves the reviewed body, normalizes machine-owned
frontmatter to the current sealed task, records both commits, and records
`runtime_profile: codex-reviewed-reuse`,
`isolation: deterministic-reviewed-reuse`, and `model_invoked: false`.
No fresh Grok review is required because the prose and source bytes are
unchanged. Finalize independently rechecks the complete prior writer/reviewer
chain and the cross-commit identity constraints.

`--reuse-candidate-from` is the source-unit-only recheck path and never invokes
a model. It accepts only an original `runtime_profile: codex`, bubblewrap-bound
candidate with the same task kind, output mapping, source-member identities,
source kind, and unit digest. The repository commit may advance when unrelated
files change; both prior and current commits are recorded, and machine-owned
frontmatter is normalized against the current sealed packet. Any source blob or
unit-contract change still rejects reuse. Derived tasks are rejected because a
safe exact-prerequisite reuse contract is not yet implemented. The body is
preserved. The new attestation records `runtime_profile: codex-candidate-reuse`,
`isolation: deterministic-candidate-reuse`, `model_invoked: false`, and the
prior manifest, packet seal, packet, candidate, writer-attestation, runner,
runtime, unit, and member-identity hashes. This mode deliberately does not
inherit any prior review: a fresh current Grok `PASS` is mandatory before
finalization, and finalize independently rechecks the full prior writer chain
plus the new review chain. It is mutually exclusive with `--profile`,
`--revise-from`, and `--reuse-reviewed-from`.

`memory.py lint --staged <transaction> --record` is the source-only bridge for
a candidate that passed Grok but failed deterministic staged lint. Recording is
allowed only for a single writer task after its writer provenance and exact
current Grok `PASS` verify. Lint must fail, and every recorded error must begin
with that task's exact output path; passing lint, missing review evidence,
multi-task transactions, and global/ambiguous errors are refused. The command
writes deterministic `record.json` and `report.md` artifacts under
`lint-failures/<task-id>/`, outside staged pages, binding the transaction,
policy/tool identities, task/output, candidate, writer/review attestations, and
exact ordered errors.

`--revise-lint-from` is mutually exclusive with the other revision/reuse modes
and requires `--profile codex`. It accepts only a source task with unchanged
task/output/source-member/unit identities (unrelated commits may differ), an
eligible direct Codex or fully verified candidate-reuse writer chain, a
hash-bound Grok `PASS`, and the trusted recorded lint `FAIL`. The candidate,
PASS review, lint JSON, and lint report are copied read-only into the new sealed
packet. The lint-specific prompt delegates only those exact errors. The revised
page is normalized to the current packet and must earn a fresh Grok `PASS`,
pass staged lint, and pass normal finalization. Derived lint revisions are not
supported.

Prepare defaults to bounded unit, writer-task, and total semantic-byte batches.
Select fewer `--units` for routine work. `memory.py update --prepare
--allow-large-batch` is the explicit override and is recorded in the transaction
manifest. Review prompt bytes are independently capped by frozen configuration.

Derived tasks run only after their declared prerequisites. The trusted runner
admits prerequisite pages append-only into the next sealed packet only when
their staged hashes and isolation attestations verify. Hydrated page and
attestation identities are written into `task.json` and the packet is resealed
as a chain from its transaction-anchored initial seal. Derived direct source
inputs are deduplicated and bounded by `read_limits.derived_task_max_bytes`.

Finalize refuses HEAD, state, policy, prompt, tool, snapshot, membership, or
page-lint drift. It journals page replacement, writes state last, and restores
the prior pages/state after a handled failure. Every command checks for and
recovers an interrupted publication before doing other work.

Run the isolated temporary-repository tests with:

```bash
python3 -m unittest discover -s memory/tools/tests -v
```

Deliberate boundaries: `memory.py` never invokes a model; only an explicit
`run_isolated.py --profile codex` command does. Codex is the only authoring
runtime. The tooling does not
parse PDFs and requires Linux bubblewrap for semantic writers. A configured derived page
removed after publication is reported as an orphan and blocks preparation;
automatic derived-page retirement is intentionally unsupported until an
explicit lifecycle/cutover task is defined. Missing source or derived inputs
remain visible as blocked tasks and prevent a full-baseline advance. Atlas
migration accounting is distinct from cutover: cutover additionally requires
full freshness, served lint/integrity, no legacy references, and recorded
benchmark plus fresh-agent retrieval passes. The tooling does not currently
produce those benchmark records, so cutover remains false until they are
supplied by a separate reviewed workflow.
