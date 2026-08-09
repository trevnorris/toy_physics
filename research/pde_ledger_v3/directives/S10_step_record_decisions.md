# S10 step record — decision list

**Artifact:** `research/pde_ledger_v3/steps/S10_two_transverse_photons.md`.
**No other file is edited.** Reading anything is in scope; editing is not.

---

## 0. What the record is for — every item below is subordinate to this

The record states **what was measured** and **what may not be claimed on the strength of it**. A reader
who has read only this record must not be able to overstate S10's result.

⇒ If satisfying any item below would make the record **easier** to overstate, that is a finding to
report, not an edit to make.

---

## 1. Every claim is sourced from an artifact that exists, **at a locus that supports it**

The record was written against machinery that has since been deleted. For each claim whose source no
longer exists, exactly one of these becomes true:

- the same quantity is **re-measured** from an artifact that exists at HEAD, and the claim is re-cited to it;
- the claim is **deleted**;
- the measurement is real but no longer reproducible ⇒ the record says so **in its own sentence**, names
  what a reader would have to run to re-establish it, and does not present it as current evidence.

⛔ **No claim is kept by annotating its dead citation.**

⛔⛔ **A resolving path is not a source.** Every cited **locus** — file *and* line range — contains an
operand that supports the claim made about it. ⚠ A live file cited at the wrong lines passes a path check
and is exactly as wrong as a dead one.

## 2. Every claim of cross-engine agreement names the instrument that compared it, and the family it compared

`scripts/out/S10_cross_engine_comparator.out` is the live cross-engine instrument.

- ⛔⛔ **A claim of cross-engine agreement for a family of objects the join does not contain says so.**
  ⚠ Agreement between two engines that were never joined is the defect class this rebuild exists to
  remove; a record can carry it just as an engine can.
- ⛔ **Rows that are unparsed, shape- or type-mismatched, or emitted by one engine only are NOT
  ESTABLISHED** — ⛔ not "a harmless disagreement", ⛔ not "accounted for".
- ⛔ **No agreement count appears without how much of it is bare integers and empty containers.**
  ⛔ **No disagreement count appears without how many carry a genuine numeric or algebraic residual, and
  whether those sit in `MAIN` or in a control package.** ⛔ **No shared-name count appears without how many
  are value-comparable at all.**
- ⚠ **Omitting the counts does not discharge this item** — a qualitative sentence (*"no computational
  conflict"*, *"every one is accounted for"*) carries the same claim and gets the same treatment.
- The record reports the comparator's **final guard state** and what it does and does not signal.
- The record states **which of the step's own load-bearing quantities are integers**, so a reader can see
  for themselves what an integer agreeing across two engines does and does not establish.

## 3. The common-mode limit sits beside **every** independently quotable cross-engine headline

Both engines take the action from **one spec, written by the orchestrator**. What a cross-engine agreement
corroborates is the route **from** the action **to** the spectrum. A defect in the action moves both
engines together and no cross-engine row can see it.

⛔ Not one chosen location. A reader quoting any single agreement claim must find this beside it.

## 4. Every bullet under *WHAT THIS STEP STILL DOES NOT ESTABLISH* is re-checked against the live instrument

For each, the record says whether it is **still open**. ⚠ A limit that merely *looks* closed because the
instrument that measured it was deleted **stays open**.

⛔⛔ **A limit the live instrument closes is restated as EXACTLY the property that instrument established,
⛔ never as the limit's headline.** ⚠ Where any part of a comparison runs an engine's object against that
same engine's other object, that part is **internal consistency** and is named as such.

## 5. The export chain is evidence this record does not have at all

S10 **imports** S9's ledger instead of re-declaring it, and **overwrites** some of what it imports.

- The record **inventories the actual overwrite operands** and **bounds the claim to those objects**.
  ⛔ It does not escalate to the imported ledger, the action, or the spectrum.
- The record states what a **form** change and a **coefficient** change to S9's action each do to whether
  S10 publishes ⇒ ⭐ **run both on scratch copies under `/tmp` and report the literal output.** The
  committed bytes do not carry this; ⛔ do not transcribe it from this list or from the front door.
- It is a **cross-step** comparison between two SymPy engines, independent because they are built from
  actions at different component counts. ⛔ **It is not a second cross-engine leg.**
- ⚠ The front door records **measured limits on the export chain**. The record carries the ones that bound
  what a refusal to publish, and what inherited assumptions, may be claimed to mean.

Sources: `scripts/S9_exports.py`, `scripts/S10_exports.py`, both engines.

## 6. The physics does not move — and that fence covers content, not evidence

The walk, the nullity table, the controls table, the dimension identity, the prior-art reading and
*Departure* are the step's result. ⛔ If a live artifact contradicts one of them, that is a **finding to
report**, not an edit to make.

⛔⛔ **The fence is the ALGEBRAIC CONTENT only.** Evidence, counts, provenance and citations sitting beside
frozen content are **always in scope and are re-measured** — ⚠ a false strength claim adjacent to a true
result is not protected by this item.

⚠ **A generic result stays adjacent to its exceptional-stratum result.** Neither may be quoted without the
other, and the record says whether the stratum layer reaches a cross-engine verdict at all.

## 7. The registry section describes a registry that no longer exists

Whatever S10 introduced there, the record says where it now lives — or that it lives nowhere and is
**owed**. ⚠ The provenance defect recorded in that section is filed against a deleted file.

---

## Inputs

- `research/pde_ledger_v3/REBUILD_HANDOFF.md` — the front door. It carries the provenance of measurements
  no committed artifact holds, and its list of measured limits on the export chain.
  ⛔⛔ **Numbers come from the artifacts, not from the handoff.** Where the handoff and a committed
  artifact disagree, **the artifact wins and the disagreement is reported.**
  ⚠ **Its S10 ledger inventory is known stale** — that one is already measured; expect others.
- `scripts/out/S10_brane_mode_spectrum_sympy_audit.out`, `mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out`,
  `scripts/out/S10_cross_engine_comparator.out` — committed inputs. ⛔ **Not regenerated.**
- `scripts/S9_exports.py`, `scripts/S10_exports.py`, both engine sources, `directives/S10_SHARED_PHYSICS.md`,
  `steps/S10_PREREGISTERED_PREDICTION.md`, `steps/S9_light_requires_shear.md` (the rebuilt record, for form).
- `reduction/` was deleted with the registry design. Anything under it is gone.

## Constraints

- ⛔ **No new checks, no new residuals, no new scripts in the tree.** Extraction scripts for reading the
  transcripts go under `/tmp` and are reported with their literal stdout.
- ⛔ **Do not modify the comparator or either engine in the tree.** Both are closed.
  ⭐ **Item 5's ablations are run on scratch COPIES under `/tmp`** — SymPy only, and the literal output is
  reported. ⛔ Never `git checkout` to restore an ablated file; copy first.
- ⛔ Do not start `wolframscript`. The Wolfram transcript is committed and is read, never regenerated.

## Acceptance

⭐ **A reader holding only the committed artifacts can locate every number in the record** — or the record
says of that number that it is not reproducible from them and names what would re-establish it.
⚠ The committed artifacts are the three `.out` files **and** the two exports modules; a number from the
export chain is sourced to the module, ⛔ not disclaimed as unreproducible.

⇒ **and no sentence claims cross-engine agreement for a family of objects the join does not contain.**

⚠ The mechanical half — every path resolves at HEAD — is **necessary and not sufficient**: it passes a
live file cited at the wrong lines, which is why item 1 puts the burden on the **locus**.
