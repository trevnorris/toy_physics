# Mathematica `DIM|` emitter — requirements and acceptance intent (v2)

**Status: USER DECISION, 2026-07-29** — a committed deliverable, second in the **build queue** behind the
ablation driver. ⚠ The build queue is not the conversion order: 027 remains the next *conversion*
(`manifests/DIMENSION_REWRITE.md` §8). "Correctness is king, regardless of cost" applies here with the
same force as D1–D5 (§1b).

⚠ **v2 after design review, which found that v1 avoided the driver's interface sprawl but reproduced its
other failure.** The review returned **five** blocking findings; the three structural ones were: **a no-op template conformed to every
requirement** (nothing required the shared machinery to be the *cause* of any record, and the retrofit
passed by changing nothing); **R7's anti-mirroring half did not bind** (a stage-local shadow table
carrying Python's names, order and hand-copied vectors satisfies every literal clause); and a process
constraint **broadened a narrow rule into an unsatisfiable one**, which is the same prose-to-infrastructure
over-correction that cost the driver three spec versions. The other two — a converted-stage count
contradiction and six under-specified contract questions — are fixed at R9 and questions 7–12. Reviews:
`research/pde_ledger_v2/_scratch/codex_wl_emitter_review{,2}.log`.
⚠ **A confirm-pass on v2 then found the anti-shadow rule was itself an over-correction** (it forbade
D2-authorized derivations) **and that v2's own "only 016 and 023" measurement was false.** Both are fixed
below. That is three over-corrections of the same shape in one night — the tell each time was a rule
written to forbid a defect rather than to state the property that distinguishes it from the legitimate
case next door.

⚠ **This document states properties and acceptance intent**, not an invocation signature or a deliverable
layout — those are Codex's, settled by a **contract** authored before fixtures are frozen. ⚠ Stated
honestly, it is not free of design consequence: R1's single-structure property and R2's before-exit
property both constrain how an emitter can be built. They are kept because each is a *measured* defect
class, not a preference. Sequence, per
`docs/development_pipeline.md` checklist 1b: **contract → fixtures frozen by a session that will not
build → build by a third session that may neither author nor weaken them.**

## Why it exists

23 stages remain, and every one currently hand-rolls its `.wl` emission block. Each has re-derived the
same hazards, and each has cost review effort proportional to that. The measured hazards, all from this
workstream's own record:
- ⛔ **43 of 43 stage-audit `.wl` files terminate in `Exit[]`** (the `mathematica/` directory holds 44
  `.wl` in total; the extra is `midway_knob_audit_codimension_mathematica.wl`, which also exits) — a print appended at end-of-file is dead code and
  emits nothing (§8).
- ⛔ **Deriving only the `axes=` label is insufficient** — stage016's first emitter printed the mapped
  label beside raw storage order, so permuting the map would have relabelled all 21 records while every
  standing gate stayed green (§4-b, `5b29f400`).
- ⛔ **The emission holder runs many times** — measured 6× (016), 13× (021), 17× (023), 19× (027),
  several under deliberately corrupted bases. An in-body print emits duplicate *and* corrupted records
  (§9).
- ⛔ **Hardcoded exponent literals in a `Print`** are vacuous against the `.py`'s constant, and two such
  defects already exist in the corpus (§5).
- ⛔ **Fractional values occur** (023 `g_U`/`g_W`, 027, 021) and must serialise as exact rationals; no
  validator currently rejects a decimal (§9, §12).
- Axis order is **per stage** and is not the same everywhere — `(L,M,T)`, `(L,T,M)`, `(M,L,T)`, a 4-axis
  `(M,L,T,E)`, and a non-physical `(stiffness,L,T)` all occur.

## R — requirements (properties)

**R1 — one structure feeds both.** The emitted `axes=` label and every emitted exponent vector derive
from a single label→slot structure in that stage. ⛔ A second, parallel axis list used only by the records
is a defect, not a fix.

**R2 — the emission is reached.** It executes in the stage's real run, before any terminal exit.

**R3 — exactly one clean set.** However many times the emitting holder is invoked, and under whatever
corrupted or mutated bases the stage exercises, the artifact contains exactly one record per emitted
quantity, from the stage's clean baseline state.

**R4 — values are live, never transcribed.** Each record's exponents come from the binding the record
names, not from a literal typed beside it. ⛔ Source-grepping for literals is dodgeable; the property is
that perturbing the binding moves the printed value.

**R5 — exactness is preserved.** Rationals serialise exactly; no float representation can reach the
artifact.

**R6 — the stage's own content stays the stage's.** The template may carry *machinery* — structure,
rendering, emission locus. ⛔ It may not carry, default, or infer any stage's basis, axis order, quantity
set or values (§10: shared machinery imported, per-stage declarations not shared, a generated table
output-only).

**R7 — the charter is untouched.** ⛔ Nothing here may import, read, or receive anything from the Python
side; Mathematica stays an independent engine (§1). D3 is not relaxed by convenience: if the template
makes it *easier* for a `.wl` to mirror its `.py`, that is a defect in the template.

**R8 — coverage is checkable against the stage's own enumeration.** The count of emitted records is
comparable to the count marked *emitted* in that stage's §4-a table, so a dead or partial emission is
detectable without anyone remembering to look.

**R9 — adoption is per stage and reversible.** Converting a stage to the template is a change a reviewer
can read as a diff, and the **seven** already-converted stages are not silently re-plumbed by it.
⚠ **Measured, and corrected twice — take the numbers from here, not from earlier drafts.** By label
construction, **004/013/018** render a literal label while **011/012/016/023** compute one. By R1's
property as written — *one structure feeding both label and vectors* — the compliant set is **at least
011/012/016/023**: 011 and 012 build their label from `Keys[dimensionAxes]` and every emitted vector from
lookups into that same association, with no second parallel axis list. ⚠ An earlier draft claimed "only
016 and 023"; that was false under R1's own wording. If the intended property is narrower — every vector
re-projected through explicit slots at serialization time, as 016 and 023 do — **R1 must say so**, and it
does not today. Bringing 004/013/018 across is a **separate, per-stage decision** with its own review,
never a side effect of adoption.

**R10 — the template must be the CAUSE of the records.** ⛔ v1's fatal gap: an empty file satisfied
everything, because each stage's existing hand-written block already emitted. Disabling or removing the
template must remove the records from a stage that has adopted it, and a stage with no prior emitter must
be able to emit **through the template alone**. Conformance is a statement about participation, not about
coexistence.

**R11 — a record's inputs are derived on the Mathematica side, not transcribed.** ⛔ What must fail: a
table assembled beside the emitter whose entries are **typed in** — names, ordering and vectors carried
across from the Python side or from anywhere else — because that satisfies every literal clause of R6/R7
while being exactly the mirroring D3 forbids. ✅ **What must remain admissible:** a computation added to
the `.wl` specifically to expose an otherwise unreachable value, provided it is its own route — **D2
expressly authorizes this** (§1b), and stage013's `kEtaDim = dimOf[data["KEta"], …]`, added inside its
emission holder as the corpus's first D2 case, is the working precedent. ⚠ The discriminator is
**derivation versus transcription**, not when the code was added or why.
⚠ **Stated honestly: this is a reviewable property, not a
mechanizable gate** — the manifest already concedes that no current mechanism establishes independent
authorship (§1). What the contract must settle is therefore *what evidence a reviewer reads* to judge
binding lineage, not a check that decides it.

## A — acceptance intent (each item able to fail)

**A1** — permuting the label→slot structure changes the emitted records, and a variant that derives only
the label fails this. *(This is stage016's measured defect; the fixture reproduces it.)*
**A2** — an emission placed after a terminal exit produces zero records and is caught, not silently
accepted.
**A3** — a stage whose holder is invoked many times, including under corrupted bases, yields exactly one
clean record set; a variant that emits in-body fails.
**A4** — perturbing a bound quantity moves its emitted exponents; a variant that renders a typed literal
fails.
**A5** — a float-valued exponent cannot reach the artifact; a variant using approximate rendering fails.
**A6** — a mismatch between emitted count and the stage's §4-a emitted count is detected.
**A7** — the retrofit oracle: applying the template to **stage023's `.wl`** reproduces its committed
`.out` `DIM|` block exactly — 29 records, `axes=L,M,T`, byte-identical after kernel-ID normalization.
⛔ A disagreement is a finding to report, never a difference to adjust either side into.
**A8** — a template variant that reads any Python-side artifact, or that supplies a default basis or axis
order, fails (R6/R7).
**A9** — **R10's attribution:** a fixture stage with no hand-written emitter emits its records through the
template, and disabling the template removes them. ⛔ A no-op template fails this, which v1's ladder did
not catch.
**A10** — **R11's lineage, as far as it is testable:** a fixture in which the emitter's inputs are a table
assembled solely for emission is distinguishable, in the evidence a reviewer reads, from one whose inputs
are the stage's own live bindings. ⚠ If the fixture author concludes this cannot be made decidable
without reading the implementation, that is a finding to report, not a gap to paper over.

## What the CONTRACT must settle — questions, not answers

1. **Distribution.** Is this a file each stage loads, or a block each stage carries? ⚠ If it is a loaded
   file it becomes shared machinery with the same trust problem the Python module had — say what pins it
   and where that authority lives, because an unpinned shared emitter is a single point of failure for
   every stage's reference half (§9, the module-pin history).
2. **Invocation and inputs** — how a stage supplies its basis, axis order, names and bindings, given R6.
3. **Rendering** — the exact record grammar it produces, and how it guarantees R5 across the axis sets
   in use, including the 4-axis and non-physical bases already cleared as non-blocking (§8).
4. **The exactly-once property** — what makes the clean state identifiable in a stage whose holder runs
   under mutated bases, without the template knowing anything stage-specific (R6).
5. **Adoption** — what converting one stage looks like as a diff, and what the **seven** converted
   stages' position is (R9). ⚠ Seven, not six: 004, 011, 012, 013, 016, 018, 023 all emit `DIM|` records
   today, and stage023 is both a converted stage and the A7 oracle.
6. **How A6 is checked**, and by which existing tool, if any — including **what authority holds the
   §4-a emitted set** and how a stage note's prose table is consumed without becoming a second source of
   truth (§10: the derived artifact is never a second source of truth).
7. **Attribution (R10).** What makes the template demonstrably the cause of a stage's records, such that
   disabling it removes them — and how a fixture exercises that without the implementation's cooperation.
8. **Binding lineage (R11).** What evidence a reviewer reads to distinguish a **derived** input — including
   a D2-authorized computation added to expose an unreachable value — from a **transcribed** one carried
   across from another engine. ⚠ If your conclusion is that no evidence short of reading the `.wl` settles it, say
   so — that is a finding, and it belongs in the contract rather than in a reviewer's head.
9. **Input validity.** Duplicate names, duplicate or out-of-range slots, vector arity mismatches, empty
   input, symbolic exponents, approximate numbers, zero denominators — what is rejected, and how.
10. **Failure behaviour.** Whether an invalid record may leave earlier records already written: fail-closed
   or partial. ⚠ A partially written artifact that still parses is the shape this workstream has been
   bitten by before.
11. **Artifact shape.** Destination, header, record ordering, line endings, and the exact extraction and
   normalization operation A7 compares under — the `.out` comparison already depends on a kernel-ID
   normalization, and A7 is not decidable until that operation is named.
12. **Trust closure.** If a shared file is loaded (question 1), what pins it, and what pins the A6 checker
   and the A7 normalizer — ⛔ otherwise the trust root is not removed, only moved into the checking
   machinery, which is the mistake §9 records for the Python module. If instead a block is copied per
   stage, what detects **copy drift** between stages, and how a template version is recorded.

⚠ Questions 7–12 exist because design review found six of them missing from v1 and two of them
load-bearing enough to make the whole ladder passable by a no-op.

## Process constraints

⛔ Acceptance tooling cannot grade itself — hence the three-session sequence above.
⛔ No expected **scientific answer or rationale** stated in a directive — ask for the determination plus
its evidence (`docs/development_pipeline.md`, checklist item 0). ⚠ **v1 broadened this to "any file a
build session can read", which is unsatisfiable and was itself the over-correction pattern:** A7's oracle,
the frozen fixtures' expected outcomes, and every stage's committed expected transcript are all builder-
readable by design, and must be.
Scripts run under `timeout 600`; a 124 is a failure to report, never a reason to raise the cap.
⚠ **≤2 concurrent `math -script` seats**, and per this session's measurement the cap must count agents,
not only Codex sessions.
