# Prior-art brief — how should a two-person project verify CAS-derived physics?

**Status:** open question, commissioned 2026-08-06. ⛔ Nothing here is a decision.
**Read first:** `CLAUDE.md` (how we work), `docs/development_pipeline.md` (the posture).
**Run this from a fresh session.** ⛔ It does not depend on the S11 session's context.

---

## 1 · What we are actually trying to do

A toy-physics ledger derives results symbolically, step by step. Each step must end with a defensible
answer to: **"is this derivation right?"**

The current answer is **dual-engine**: two independently built CAS scripts — one Wolfram, one SymPy — are
written from **one shared specification**, run, and their emitted objects compared by a harness. The
disagreement is the measurement.

⭐ **The question this brief asks: is that the right architecture, and is there a known better one?**

⚠ **We are not committed to it.** Redoing earlier steps (S9, S10) under a better method is explicitly
acceptable if it makes a step faster. ⛔ Do not treat anything currently in the repository as a constraint.

## 2 · The measured problem — why this is being asked

⚠ **Roughly one full day per step goes into making the scripts and the comparison work, and almost none
into physics.** We are about a week in with near-zero physics output. The numbers below are from this
repository and from the S11 session; ⭐ verify any you intend to rely on.

**The defects we find are not physics defects.** In S11, about **25** defects were found and fixed across
the specification, both engines, the comparator and the registry. Essentially **none** were wrong physics:

- Engine 1: **9** defects over 3 review rounds and 3 fix rounds — all *declaration wiring* except one.
- Engine 2: **7** defects — one structural, the rest declarations, missing premise tags, conditional
  emission.
- The shared spec: **9** defects over 3 repair rounds — all naming, pinning and section interactions.
  **3 of those 9 were introduced by the repair itself.**
- ⭐ Meanwhile an independent reviewer reconstructed the dynamical matrices, determinant and roots from the
  supplied action and **reproduced both engines' values**. ⇒ **the physics was right in both engines on the
  first build.**

**The comparison layer does not converge.** `research/pde_ledger_v3/reduction/checks_S10.yaml` is **3,121
lines** and its own header concedes it compares a **subset**; coverage was **545 / 690** rows, **145 rows
never reached a verdict**, and **11 reported DISAGREE purely on nullspace-basis representation** — two
correct answers written differently.

**Naming is the recurring failure.** The two engines emitted **17 vs 22** premise families with **5 names
shared**. Three rounds of pinning names in prose twice reintroduced divergence while closing it. The S10
config already recorded this as a known spec gap; S11 inherited it unchanged and it recurred.

⇒ ⭐⭐ **The suspicion to test: most of the cost is in making two independently-written programs agree on
how to SAY things, not on what is true.**

## 3 · What to look for

⭐ Search each area. For each, we want: **what the established practice is, whether it applies to a
two-person project with no infrastructure team, and what it would cost us.**

**A · Differential testing and N-version programming.** Running independent implementations and comparing
outputs is an old, named practice. ⭐ Find what the literature says about **correlated failure between
"independent" versions**, and specifically about **the shared specification as the correlation channel** —
that is exactly our architecture, and our spec generated 9 defects while the physics was right. ⭐ What
mitigations are known? Does the literature recommend specification *diversity* rather than more careful
specifications? *(Starting points believed relevant, ⛔ verify rather than assume: McKeeman on differential
testing; Knight & Leveson 1986 on correlated N-version failures.)*

**B · Cross-CAS verification in practice.** ⭐ How do people actually compare results across Mathematica,
SymPy, Maple, Maxima? Are there established test suites, benchmark corpora, or published methodologies?
*(e.g. integration test suites of the Rubi kind; CAS comparison studies.)* ⭐ What do they do about
representation differences?

**C · ⭐⭐ Canonical semantic interchange.** ⭐ Formats designed so two systems can exchange mathematical
objects **without agreeing on names** — OpenMath, content MathML, and any practical toolchain between
SymPy and Wolfram. ⭐ **This is the most directly relevant area**: we have been pinning tag spellings in
prose to solve a problem that may have a standard representational answer. ⭐ Report whether a usable
round-trip actually exists today, or whether it is a paper standard with no working bridge.

**D · ⭐⭐ Deciding whether two symbolic expressions are the same object.** ⭐ This is the irreducible core
of what we need. ⭐ Look for practical zero-testing and equivalence methods — normal forms, randomized
evaluation at rational or finite-field points (Schwartz–Zippel style probabilistic identity testing),
certified simplification, SMT/polynomial-identity approaches. ⭐ **Specifically: is random numeric probing a
recognised substitute for symbolic simplification when comparing two CAS?** If so it would let us compare
*values* and stop negotiating *representations* — and would make basis-representation disagreements
disappear.

**E · Formal methods, honestly costed.** ⭐ Proof assistants (Lean/mathlib, Coq, Isabelle) and SMT solvers.
⭐ We need a realistic verdict on cost for a two-person toy-physics project, ⛔ not an endorsement. Is there
a lightweight subset — e.g. checking one algebraic identity per step — that buys real confidence cheaply?

**F · How physics actually validates symbolic derivations.** ⭐ Look at communities that derive by CAS —
GR tensor algebra (xAct/xTensor, Cadabra), particle physics (FeynCalc, FORM). ⭐ **What do they do instead
of dual-engine?** Known limits, special cases, dimensional analysis, comparison against published results,
independent re-derivation by a person? ⚠ **This may be the most valuable answer in the brief**: the
alternative to a better comparator might be a different validation strategy entirely.

**G · Provenance.** ⭐ We maintain a registry mapping declared quantities to the code that produced them,
by file path and line range — which breaks whenever a script is rewritten. ⭐ Is there an established
pattern for computational provenance that does not key on line numbers? *(W3C PROV, workflow-manager
provenance, research-object patterns.)*

**H · The whole activity.** ⭐ Has anyone published on **agent-driven symbolic derivation with automated
verification**? Thin literature expected, ⛔ but worth one pass — including negative results and
cautionary write-ups.

## 4 · What to bring back

⭐ **A single document, under ~4 pages**, containing:

1. ⭐ **For each area A–H: what exists, whether it is real and usable today, and what it would cost us.**
   ⛔ Not a literature review — a practitioner's assessment.
2. ⭐⭐ **A direct answer to three questions:**
   - Is **dual-engine** a recognised practice, a known anti-pattern, or an unusual choice? Does anything in
     the literature predict our specific failure — two implementations agreeing on the physics while
     diverging on presentation?
   - Is there a **standard way to compare two CAS results** that does not require agreeing on names? If
     yes, name it concretely and say what it would take to adopt.
   - Is there a **cheaper validation strategy** that gives comparable confidence — such that dual-engine
     could be reduced or dropped for most steps?
3. ⭐ **A recommendation on comparator architecture**: one universal harness across all steps, versus a
   per-step comparator over a small shared library. ⚠ Our evidence points to the latter; ⭐ say whether
   prior art agrees or disagrees, and why.
4. ⭐ **Anything that suggests we are solving a problem that does not need solving.**
5. ⭐ **What you could not find**, and where you looked. ⚠ **NOT-FOUND is not evidence of originality** —
   report it as absence of evidence.

## 5 · How to treat the findings — ⛔ this constrains the whole brief

⛔⛔ **PRIOR ART IS AN ORACLE, NEVER A PREMISE** (`CLAUDE.md` rule 16). ⭐ Check our approach against what
you find; ⛔ never assume a published result holds for our object, and ⛔ never adopt a method because it is
standard. Its conditions may not be ours.

⭐ **Apply the governing test to every candidate**, and label each explicitly:
- **Does it catch a way the PHYSICS could be wrong?** → high value.
- **Does it only reduce TOOLING friction?** → ⭐ still valuable here, that is the whole point of this
  brief — ⛔ but say which it is, and never present the second as the first.

⛔ **Do not build anything. Do not modify the repository.** ⭐ This brief produces a document and nothing
else. Adoption is a separate decision.

⚠ **Constraints any recommendation must respect:** two people, no infrastructure team; agents do the
building; the reviewer must not be the builder; falsification is a first-class outcome, so nothing may
require the answer in advance; and a method that takes a week to stand up must save more than a week.

## 6 · Where this came from

⭐ Commissioned after the S11 session, in which the physics was correct in both engines from the first
build and the entire day went to making the two engines and the shared specification agree on how to
report it. ⭐ The open question is whether that is intrinsic to the architecture or an accident of how we
have implemented it.

⇒ ⭐ Findings feed a **post-S11 method assessment**. ⚠ Going back and redoing S9 and S10 under a better
method is explicitly on the table.
