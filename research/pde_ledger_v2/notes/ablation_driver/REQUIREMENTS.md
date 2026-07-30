# Ablation driver — requirements and acceptance intent (v5)

**Status: USER DECISION, 2026-07-29, RE-SCOPED SMALL 2026-07-29/30**
(`manifests/DIMENSION_REWRITE.md` §12b(b)). Leads the build queue; existing hand-rolled harnesses switch
over to it; every subsequent stage reuses it. It is built before stage027 begins, and 027 is the first
conversion to use it — that is the build queue, not the conversion order (§8).

⛔⛔ **v5 — USER DECISION 2026-07-29/30: BUILD IT SMALL. This is the whole specification.**
**Mutate a declaration, confirm the declared assert fires, record it. Reviewed by one fresh agent.**
⛔ **No contract, no frozen fixtures, no three-session shape, no freeze tower, no byte authority.** What was
cut, and why: exact **resume**, crash **repair**, the outcome **truth table**, the **digest restore proof**,
the evidence-prose infrastructure, and the whole acceptance-tooling separation. None of them can catch a
wrong derivation — they catch ways the *tooling* could be wrong, and that is the layer this project stopped
buying (`docs/development_pipeline.md`, *THE POSTURE*).
⚠ **The prior artefacts are DELETED (2026-07-30), not superseded-in-place:** `CONTRACT.md`,
`CONTRACT_NOTES.md`, `CONTRACT_NOTES_V4.md`, `FIXTURES_REPORT.md`, and the whole `fixtures_v4/` suite
with its `fixtures_v4.accepted.sha256`. They are recoverable from git history and nothing needs them.
⇒ **What was worth keeping out of them is in the Appendix at the foot of this file** — the stage023
legacy mapping (was `CONTRACT.md` §C-9), the builder hazards, and the able-to-fail shape.
⛔ Nothing in that appendix is a gate.

⚠ **What the earlier rounds correctly established, and how the surviving set is shaped.**
v1 and v2 failed design review for trying to fix an **interface** — argv, wire format, sentinels, outcome
predicates, column mappings — which is design, and design of deliverable code belongs to the builder
(`docs/development_pipeline.md`, Phase 1 step 1: the orchestrator *"never pre-designs the computational
route"*). v3/v4 then answered that by **sequencing** the work across three sessions with frozen fixtures in
the middle, which is what grew the tower. ⇒ **The interface is the builder's to choose and state.** What
the orchestrator specifies is the property set below.
⚠ **Two counts live here and must not be collapsed — an earlier draft said "the four surviving
requirements", which contradicts the set it then enumerates:**
- **SIX requirements survive** — **R1, R2, R3, R4, R6, R7** (stated at the foot of §R). That is the live set.
- **FOUR of them were *measured* to decide real verdicts** — **R2, R3, R6** plus acceptance item **A6**.
  Those four are the reason the driver exists at all; R1, R4 and R7 survive because they are what makes the
  loop honest and its output readable, not because a measurement forced them.

**What it is.** A per-tooth ablation runner: it takes an orchestrator-supplied list of mutations, applies
them one at a time to the real file at its real path, runs a producer and one or more checkers after
each, and records one result row plus retained captures per mutation.

**What it is not.** It never decides *what* to mutate. Tools whose `init` enumerates the sites were
rejected partly for that: the target list belongs to the orchestrator (`development_pipeline.md`
per-gate checklist item 6).

---

## R — requirements (properties the driver must have)

**R1 — the list is an input.** The driver never scans for candidate sites and never adds a row.

**R2 — the mutation is visible at the real file's real path while the producer runs.** ⭐ **KEEP —
load-bearing.** This workstream's controls hash the file's own bytes and derive artifact paths from
`__file__`; a mutant living anywhere else makes those controls pass vacuously or fail uniformly.

**R3 — no state from a previous mutant is observable by the next** — no artifact it emitted, no stale
bytecode. ⭐ **KEEP — load-bearing:** its absence let a residual artifact, rather than a comparison, decide
stage023 verdicts. ⚠⚠ **The "16 of 22" figure is SELF-REPORTED and repository-unverifiable** — the defective
run was gitignored and is not in the tree (see *Evidential note*); cite the failure class, not the count.

**R4 — captures are per mutant and are never overwritten.** (This is the "record it" half.)

**R5 — ⛔ CUT 2026-07-29/30.** ~~resume re-runs no already-banked row, completes every unbanked row, and
yields the same result file as an uninterrupted run.~~ Exact resume is tooling robustness; re-run the list.

**R6 — a mutation that does not apply exactly once is a recorded outcome**, distinguishing zero from
multiple matches, never recorded as a mutant that ran. ⭐ **KEEP — load-bearing.**

**R7 — restore the file after every mutation**, on normal end and on error. ⛔ Never by `git checkout`,
`restore` or `stash` — this workstream ablates **uncommitted** work and those commands restore from HEAD and
destroy it [`never-checkout-to-restore-uncommitted`]. Use `cp`. ⚠ **CUT 2026-07-29/30:** the **digest
proof** of the restore, and the separately invocable **repair** operation for the uncatchable-kill class.
The `cp` restore stays; it is data safety, not ceremony.

**R8 — ⛔ CUT 2026-07-29/30 as a requirement; it survives as a PROSE rule** (see Process constraints).
~~the evidence claims exactly what it supports~~ — the defect it was written against was a claim that
outlived its evidence, and that gets a prose fix, not an output-format property.

**R9 — ⛔ CUT 2026-07-29/30.** ~~no row's outcome asserts more than its own fields support~~ — this
required the contract's outcome truth table to exist. Record what you observed; do not invent an outcome
vocabulary that needs a schema to police it.

⚠ **The numbering is deliberately UNCHANGED** so that every existing `R2`/`R3`/`R6` citation elsewhere in
the corpus still resolves. **The live set is R1, R2, R3, R4, R6, R7.**

---

## A — acceptance intent

Each item states **what must be demonstrated by execution**, and each must be able to fail: say which
implementation deficiency it catches, and show it catching it. ⚠ The **mechanics** are the builder's to
choose and to state in the build report — there is no contract for them to bind to.

**A1** — R3 is load-bearing: a residual artifact that would otherwise change a mutant's outcome does not.
**A2** — R6's zero-match and multi-match floors fire, distinguishably.
**A3** — ⛔ **CUT** (resume, old R5).
**A4** — restore holds after a clean run and after an error (R7). ⚠ The uncatchable-kill/repair leg is
**cut**.
**A5** — every mutant's captures survive the run and correspond to that mutant (R4).
**A6** — ⭐ **KEEP — the real committed include-list validates, with a floor that an accept-and-ignore
parser fails:** all 51 rows, the two axis counts, and the exact mutation text of named rows reproduced.
⛔ **A parser must be tested against the REAL committed file, not a synthetic one**
[`parsers-need-real-input`].
**A7** — ⭐ **THE RETROFIT — a cross-check worth running, DOWNGRADED from a blocking gate (2026-07-29/30).**
Re-run stage023's two ablation axes and compare against the committed
`notes/stage023_step_h_evidence/results.tsv` (22 `A1_DECLARATION` + 29 `A2_BINDING` rows). ⭐ **The property
worth having, and it is cheap:** the oracle is **already committed**, so nobody gets to decide the right
answer afterwards, and a driver that reproduces the old verdicts has demonstrably replaced the harness.
⛔ **What is NO LONGER required:** exact agreement on every field of the 13-column legacy schema via the
old mapping table. That is **tooling-replay fidelity** — it cannot catch a wrong derivation, and it is
beyond the four *measured* properties (R2, R3, R6, A6) this driver exists for. ⇒ Run the comparison,
**report every disagreement**, and ⛔ never adjust either side to close one; a remaining disagreement is a
finding, not a blocker. **Appendix A.1** keeps the mapping as a **reference** for which legacy column
corresponds to which new observation — ⛔ not an authority, and nothing in it is binding.
**A8** — ⛔ **CUT** (old R8's evidence audit; the prose rule replaces it).
**A9** — ⛔ **CUT** (old R9's truth table).

⚠ **Numbering unchanged, for the same citation reason as R.** **The live set is A1, A2, A4, A5, A6, A7.**
⛔ Also cut: the whole notion of a separately authored and separately graded fixture suite.

---

## Evidential note — what motivates R3, and what survives

⚠ Measured in a gitignored scratch tree. State it as motivation, never as reproducible evidence:
- **R3** — in the pre-fix run of stage023's declaration axis, most checker verdicts were freshness
  failures decided by a residual artifact rather than comparisons. The committed `results.tsv` is the
  **corrected** run; the defective one is not in git.

---

## Process constraints

⭐ **One builder, one fresh reviewer** (`docs/development_pipeline.md`, Roles table). Whoever writes the
driver does not review it; the review runs on a **fresh** agent. ⛔ That is the whole separation — there is
no three-session sequence, no independently frozen acceptance suite, and no "may neither author nor weaken"
clause.

⛔ **The orchestrator owns the mutation target list**, never the builder
(`docs/development_pipeline.md` per-gate checklist item 6, [`per-tooth-ablation`]). This is the one
ownership rule the driver's design must not quietly take over.

⛔ **A directive to the builder still states no expected answer, reason or premise** — it states what to
determine and requires the evidence for that determination [`never-supply-the-expected-reason`]. ⚠ The
acceptance oracles here are ordinary committed files and the builder may read them; what is forbidden is
*supplying the conclusion* in the directive.

⚠ **A PROSE rule, deliberately not an engineering requirement:** whatever a run leaves behind must claim
only what the **committed** files support, not what was true of the scratch directory it ran in. Four claims
in stage023's hand-written evidence summary were true of the run and false of the commit and were corrected
under review. ⛔ **Do not convert this into an output-format requirement again** — that conversion is what
pulled in transitive read sets, replay capsules and platform prerequisites and took a ~100-line tool to a
full interface contract. **A prose defect gets a prose fix** [`fix-defect-at-its-own-level`].

⚠ **Recorded cost, kept because it is honest:** an implementer who reads the oracles can special-case every
disclosed input and pass every item without establishing general behaviour. ⛔ **That is a motivated-builder
threat, and this project does not harden against it** — the operative risk is drift and honest error
(`docs/development_pipeline.md`, *THE POSTURE*). The fresh review leg is the mitigation, and it is enough
for a two-person project.

⛔ **Scripts get `timeout 600`; a 124 is a failure to reformulate, never a raised cap**
[`script-timeout-policy`]. ⚠ The bespoke suite-wide wait-allowance bound that used to live here is **cut**
(2026-07-29/30) along with the suite it bounded.

---

# Appendix — EXTRACTED HISTORY

⚠ **What this is.** The superseded contract, its two notes files, the fixtures report and the `fixtures_v4/`
suite were **deleted 2026-07-30** (recoverable from git history). This appendix keeps only what would
otherwise have been lost. ⛔ **None of it is a requirement, a gate, or an authority** — it is reference for
whoever builds the driver. ⛔ The freeze, the external byte authority, the three-session separation and the
outcome truth table are **deliberately not reproduced**; that layer is what was retired.

## A.1 — the stage023 legacy oracle and the A7 mapping (was `CONTRACT.md` §C-9)

**The two committed legacy tables**, both under `notes/stage023_step_h_evidence/`:
- `include_list.tsv` — 5 columns `axis · key · record · old_text · new_text`, **51 data rows**,
  **22 `A1_DECLARATION` + 29 `A2_BINDING`**.
- `results.tsv` — 13 columns `axis · key · record · stage_exit · pass_count · fail_count · first_fail ·
  sidecar_written · record_moved · emitted_value · cmp_exit · cmp_status · mismatch_names`, same 51 rows.

Join the two on `(axis,key)`; keys are unique and `record` agrees across both.

**Mapping — legacy column → the observation a new driver would make.** Right-hand side is stated as an
observation, not as a field name, because the schema it used to map into is gone with the contract.

| Legacy field / value | Corresponding new observation |
|---|---|
| `axis`, `key`, `record` | the same three, exact strings |
| `stage_exit` | the producer's exit status |
| `pass_count`, `fail_count` | the producer's PASS/FAIL tally |
| blank `first_fail` | no failed assertion was named |
| nonblank `first_fail` | the name of the first failed assertion, exact string |
| `sidecar_written=yes` / `=no` | the `stage023_sidecar` artifact was / was not emitted |
| `record_moved=yes` / `=no` | the record's value did / did not move |
| `record_moved=no_sidecar` (with blank `emitted_value`) | no sidecar, so no moved-value observation exists |
| nonblank `emitted_value` | the value after the move, exact string |
| `cmp_exit` | the comparator checker's exit status |
| `cmp_status` | that checker's verdict status, exact string |
| `mismatch_names=none` | the empty mismatch list |
| any other `mismatch_names` cell | a one-element mismatch list holding that exact cell |

⚠ **Legacy sentinels, so they are not misread:** `yes`/`no` are booleans; `none` is the empty-list sentinel,
not a missing value; a blank cell is meaningful only where the table above says so. `no_sidecar` together
with a nonblank `emitted_value` is an inconsistent combination and would be a comparison error, not a row.

**Measured distribution of the committed `results.tsv`** — derived from the file, re-verified 2026-07-30:
**16** rows with `stage_exit != 0`; **35** with `stage_exit = 0` and `cmp_exit != 0`; **0** with both zero.
The same file has **35** `sidecar_written=yes` / `record_moved=yes` and **16** `no` / `no_sidecar`.
⇒ Useful as a shape check on any re-run before comparing row by row.

⚠ **Why A7 is a mapped comparison and never reproduction.** The legacy table has **no field that
establishes a mutation applied** (no match count) and **no outcome column**. Any successor with a richer
schema is structurally incapable of byte equality with it, which is why the "exactly" in earlier drafts had
to be read as *on mapped fields* — and why v5 dropped even that and asks only that disagreements be
reported.

## A.2 — hazards a builder needs (from the two notes files and the fixtures report)

1. ⭐ **The real include-list has no target-file column.** R1 says the list is an input; it does **not**
   mean each row names its own target — the committed 5-column list cannot. The target comes from
   run-level configuration. If a per-row target is ever wanted, the list format changes and the stage023
   oracle stops parsing.
2. ⭐ **Fix who may write what, or attribution breaks.** If the producer may rewrite its own mutant, or a
   checker may create the artifact, a row's evidence no longer attributes the output to the producer.
   Observe the artifact **before** running any checker.
3. ⭐ **R3 has a boundary — state it, do not claim hermeticity.** "No state from a previous mutant is
   observable" is decidable only over state the run *declares*: the target, the artifact, reset paths, and
   bytecode caches. An unconstrained child can read host state no driver can enumerate.
4. ⭐ **R7 has a boundary, and cutting `repair` has a consequence.** Cleanup after a catchable signal is
   executable; cleanup after `SIGKILL` or machine loss is not, because the terminated process cannot run
   it. ⚠ **A hard-killed run therefore leaves the target MUTATED on disk**, and this workstream ablates
   **uncommitted** work. Recovery is a manual `cp` back from the copy R7 requires — ⛔ **never**
   `git checkout`/`restore`/`stash`, which restore from HEAD and destroy it
   [`never-checkout-to-restore-uncommitted`].
5. **R3 acceptance cannot prove internal ordering.** It can show reset happened before the producer ran; it
   cannot distinguish two internal orders with identical observable bytes. Do not write an acceptance item
   that claims to.
6. **A6's concrete failure mode, measured:** a whitespace-tokenizing parser reads the real 51-row list into
   rows of width 7/9/11/13 — `record`, `old_text` and `new_text` contain spaces. ⇒ The floor is **51
   losslessly parsed rows each preserving all five TSV fields**, not a successful open and not a line count.
   ⚠ A6 is **parse-only**; actually applying the 51 mutations needs a production target, producer and
   checker, and that is A7's job, not A6's.

## A.3 — the able-to-fail shape that worked (was `NEGATIVE_DEMONSTRATIONS.md`)

⭐ **Three steps per acceptance item, in this order:** accept a good observation → **reject a twin deficient
in exactly one property** → restore that one property and watch the same twin become acceptable. The third
step is what separates a real discriminator from one that rejects everything put in front of it.
⚠ Running that sweep found several harness predicates that were **entailed by their own harness** and so
could not fail; they were deleted rather than kept as reassurance
[`decisive-test-not-tautological`, `script-review-depth`].
