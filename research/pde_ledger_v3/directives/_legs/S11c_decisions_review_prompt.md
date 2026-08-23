# Independent review — the S11c decision list (RE-VALIDATION, not transcription check)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md`

## What this artifact is, and why this review is load-bearing
S11c is the **next** step of a dual-engine physics ledger: the **non-uniform, variable-coefficient
transverse coupling** at a brane–bulk interface — *is light's confinement unconditional?* The immediately
prior step **S11b** (the uniform interface coupling law) is CLOSED; it proved the uniform transverse
coupling is **identically zero**, so the non-uniform, gradient-driven channel is exactly what S11c must
compute.

This artifact is S11c's **decision list** (orchestrator-written): it settles structure, naming, provenance,
and scope for the S11c spec and build directives that follow. ⚠ **It is the ONE artifact the eventual
builder trusts** — everything downstream of it is reviewed twice, but the decision list itself is reviewed
only here (these two legs), then folded once and used. An error you miss propagates into both engines.

## The CORE task — RE-VALIDATE the five requirements, do NOT merely check transcription
The decision list's central job is to re-validate five requirements for the non-uniform problem (items
`N3`–`N7`). Those requirements were named by earlier reviews and preserved as **historical input** in
`steps/S11c_SCOPE.md` and `steps/S11b_HANDOFF.md` — ⚠ **that package is explicitly NOT ratified** (the
unified S11b decision list ratified only C's *scope*, not this requirement set;
`steps/S11b_HANDOFF.md:12-15`). ⛔ **"The list faithfully copied the scope doc" is NOT a pass** — the scope
doc is itself unratified. Your job is to independently judge the **physics**.

⭐ **Derive your own view first.** For each requirement, decide from first principles whether it is (a)
**correct physics** for a non-uniform slab, (b) **load-bearing** on the transverse–thickness coupling, and
(c) **complete** — is any necessary requirement MISSING (a sixth ingredient the five omit)? A requirement
you cannot independently justify — or one you would state differently — is a finding. Say, in your own
words, *why* each requirement you accept is correct for the non-uniform problem; a bare "looks right" is not
a review.

## What you are handed (read the sources FIRST, then the artifact)
Read these and form your own view of what S11c must require **before** opening the decision list:
- `research/pde_ledger_v3/steps/S11c_SCOPE.md` — the consolidated scope (source of the five requirements,
  the carry-ins, the scope boundary, the split-finer lesson). Unratified historical input.
- `research/pde_ledger_v3/steps/S11b_HANDOFF.md` — its "WHAT S11b-C MUST HANDLE" section (lines ~66–89) is
  the earlier-review origin of the requirements and carry-ins.
- `research/pde_ledger_v3/steps/S11b_interface_coupling_law.md` — the CLOSED S11b step record: the uniform
  result S11c builds on (esp. the transverse mode identically-zero coupling, and the `O(v₀|q_n|/ω)` limit).
- `research/pde_ledger_v3/directives/S11b_unified_decisions.md` — the S11b decision list; `G14` is the scope
  boundary S11c inherits, `G12(a)` is the "presuppose the FORM of an answer" defect, `G8(a)` the T7
  comparator contract, `G3` the `v_0` false-equal naming hazard.
- `research/pde_ledger_v3/directives/_measurements/S11c_decisions.md` — the rule-2 twin: every artifact-claim
  in the decision list with the command that produced it. Spot-check that the citations resolve to what the
  list says they do.
There is **no do-not-read list** — read anything in the repo that helps you judge the physics.

## Required method — DOCUMENT branch
Read the sources of truth first, form your own view of what the non-uniform physics demands, and **only
then** read the decision list. Quote **both sides** (source and decision-list line) for every finding.
⚠ This reading order is a method instruction, not a blindness control — do not describe yourself as "blind."

## Four traps to hunt by name (each has cost a real defect)
1. **Form-presupposition (`N5`).** Does any decision presuppose the FORM of an answer the engines are meant
   to COMPUTE — above all, does anything demand a **global dispersion relation**? That exact error killed
   **two** S11b directive revisions. (It is the `G12(a)` defect in a new place.)
2. **A vacuous replacement control (`N6`).** The uniform-limit control is known-vacuous (the uniform
   coupling is identically zero). The list requires a REPLACEMENT control that can actually **fail** — a
   FORM/gradient-order change or a one-sided corruption of the gradient channel. Is what the list requires a
   genuine independent-route / form check, or a coefficient-only / definition-vs-its-own-substitution
   tautology that passes vacuously?
3. **A leaked bound (`N7`).** S11c's coupling coefficient is bounded by a **bench-top-optics** measurement;
   that numeric bound is a withheld acceptance criterion (rule 5) — the builder computes the coefficient and
   we diff on our side. Does ANY phrasing in the decision list state, narrow, or hint at that bound (a
   number, an order of magnitude, a threshold, an "expected" value)? Conversely, is the falsification still
   stated clearly enough to be actionable?
4. **The split (`N2`).** The list does NOT pre-commit the sub-step decomposition — it asks YOU to recommend
   one. S11b took **eleven** directive revisions because its surface was too large. **Return a recommended
   split** of S11c into tightly-scoped, independently-spec-able sub-steps (or a reasoned finding that a
   single pass is tractable). Are the candidate seams the list names genuinely finer than S11b's, and does
   each produce an independently spec-able piece — or does a seam re-merge the surface?

Plus the naming trap inherited from `G3`: any place a **new** non-uniform quantity reuses an imported S11b
key such that the export chain's `F9` object-comparison would prove a **false equal** and silently merge two
different quantities.

## Physics filter
Report a finding only if it catches a way the **physics or the method** could go wrong — not "this could be
wrong on a different input," and not stylistic preference. But ⛔ a leg that returns "nothing survives the
filter" is weak evidence: state explicitly what you checked, what you derived, and where you looked hardest.

## Output
1. **Findings**, most-serious first — each with the source quote, the decision-list line, and the concrete
   way the physics/method goes wrong if uncorrected.
2. **Your independent verdict on each of the five requirements** `N3`–`N7`: correct / needs-correction /
   wrong / incomplete — with your first-principles reason.
3. **Your recommended S11c sub-step split** (`N2`).
4. Whether any requirement is **MISSING** from the five.
