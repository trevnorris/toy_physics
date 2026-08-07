# Prior-art findings — how should a two-person project verify CAS-derived physics?

**Answers the brief at `docs/method_prior_art_brief.md`.** Written 2026-08-06 from a fresh session.
⛔ **Nothing here is a decision, and nothing was built or modified.** Adoption is a separate call.

⚠ **Provenance of the numbers below.** Repo figures were re-measured, not taken from the brief. Papers
marked ⓕ were fetched and read; papers marked ⓢ are from search summaries only and should be opened before
being leaned on. Several 2026 citations are recent preprints.

---

## 0 · Two measurements taken here first, because they change the reading

**① The comparator compares NAMES, then SIMPLIFIES.** `reduction/engine_output_checks.py:728` is
`sp.simplify(lhs - rhs) == 0`, reached only after a naming-negotiation layer — `_snake_to_lower_camel`
(`:799`), `_build_naming_rule` (`:815`), `_build_symbol_identities` (`:867`), `_symbol_bijection`
(`:1157`). Confirmed: `reduction/` is **8,081 lines**; `checks_S10.yaml` is **3,121** and its own header
concedes a subset.

⇒ ⭐⭐ **Our architecture is the two things prior art independently identifies as the expensive and
unreliable path: a nominal join key, and symbolic simplification as the equality test.** Everything below
follows from that.

**② The granularity is unprecedented.** S11 engine 1 emits **3,750 tags** (`REBUILD_HANDOFF.md`). ⭐ **No
community found in this survey compares two CAS at the level of intermediate tagged emissions.** Every
corpus compares a *small number of final objects* against a reference. This is the single biggest
divergence between us and everyone else, and it is a cost driver nobody else pays.

---

## 1 · Areas A–H

### A · N-version programming and differential testing — **real, named, and it predicts us**

⭐ **Dual-engine is a recognised practice**, not an invention and not an anti-pattern. It is N-version
programming (Avizienis) / differential testing (McKeeman).

- **Knight & Leveson 1986** ⓢ: 27 independently written versions, 1M random tests. Individually very
  reliable; **coincident failures far exceeded the independence assumption.** The named mechanism is ours:
  *because all versions are developed from the same specification, and specifications are finite and
  sometimes ambiguous, there exists a nonzero set of inputs on which they are commonly misinterpreted.*
- **Hatton 1997** ⓢ argued the opposite economics (N ordinary versions beat 1 excellent one); later work
  refutes/qualifies it. ⇒ The literature is **divided on cost-effectiveness** and **unanimous on the
  mechanism**. ⛔ Do not read "recognised practice" as "settled good idea."
- ⭐⭐ **The mitigation literature does NOT recommend more careful specifications.** It recommends moving
  the oracle out of the specification — see **metamorphic testing** (§F/§4).
- ⭐⭐ **Brilliant, Knight & Leveson, *The Consistent Comparison Problem in N-Version Software*
  (1986/IEEE TSE 1989)** ⓢ — the closest published thing to our failure, and it is **forty years old**:
  *an N-version system may be unable to reach consensus even when none of its versions has failed.*
  ⇒ **The literature predicts our failure at the COMPARATOR, not at the physics.** Our 11 nullspace-basis
  DISAGREE rows are this paper's phenomenon exactly.
- **2026, and directly on us — "N-Version Programming with Coding Agents"** ⓕ (arXiv 2606.20158):
  48 agent-generated implementations, 1M inputs. **Substantial common-mode failure**, and the co-occurring
  failures **trace to where the specification is hard or ambiguous**. Still worth doing: single version
  387.44 failures → three-version majority **130.99**.
- **MAST multi-agent failure taxonomy** ⓢ: **Specification Issues 41.77%** — the largest of three
  categories, ahead of inter-agent misalignment (36.94%) and task verification (21.30%).

**Verdict.** Keep dual-engine as a *concept*; it catches **PHYSICS**. ⛔ But our spec-repair loop (9 defects,
**3 of them introduced by the repair**) is the exact failure mode the 2026 agent study measures, and it is
**TOOLING**. Prior art says you cannot write your way out of it. Cost of doing nothing: recurring, per step.

### B · Cross-CAS verification in practice — **established, and pairs by INPUT not by NAME**

- **Rubi** ⓢ: 72,000+ integration problems, translated into several CAS syntaxes, each stored with an
  **optimal antiderivative**. Abbasi's *Computer Algebra Independent Integration Tests* ⓢ run Rubi,
  Mathematica, Maple, Maxima, Fricas, SymPy, Giac, Sage against it and grade A/B/C/F.
- ⭐⭐ **The structural point: there is no naming negotiation anywhere in this, because the join key is the
  PROBLEM, not the output label.** Two systems are compared because they were handed the same integrand.
- ⚠ Their reported reason for adding numeric testing: **CAS simplification verification is not effective
  enough on its own** (§D).

**Verdict.** Usable pattern, near-zero cost: **anchor every cross-engine comparison on a shared input, and
let the output labels be whatever each engine wants.** **TOOLING** — but it deletes our recurring cost.

### C · Canonical semantic interchange (OpenMath / content MathML) — ⛔ **paper standard, no bridge**

- OpenMath is real: XML, 200+ content dictionaries, formal properties. Content MathML likewise.
- ⭐ **There is no maintained OpenMath↔SymPy converter.** The SymPy community answer is *"there isn't an
  out-of-the-box way; parse the XML yourself and write the transform."* Same for Wolfram.
- SymPy↔Wolfram in practice is **one-way and lexical**: `parse_mathematica` (WL→SymPy, maintained) and the
  `mathematica_code` **printer** (SymPy→WL). Not a semantic round-trip. ⚠ We already do the WL→SymPy leg
  by hand — `engine_output_checks.py` parses WL FullForm into SymPy trees.

**Verdict.** ⛔ **Do not adopt OpenMath.** Answering the brief's question directly: it is a paper standard
with no working bridge for our purpose. **TOOLING**, and unavailable.

### D · ⭐⭐ Deciding whether two symbolic expressions are the same object — **THE payoff area**

**Yes — random numeric probing is a recognised substitute for symbolic simplification, and the people who
tried symbolic-first published that it was not enough.**

- **Schwartz–Zippel / polynomial identity testing** ⓢ: evaluating a nonzero polynomial at random points
  from a large set almost certainly gives nonzero. It was **the first randomized polynomial-time algorithm
  proven correct** for identity testing. This is textbook, not a hack.
- ⭐⭐ **Measured, and the most useful number in this document — DLMF/LaCASt** ⓕ (arXiv 2109.08899):
  of **2,405** translated formulae, **symbolic** simplification verified **660 (27.4%)**; **numeric**
  testing on the remainder verified **418 more**. Their own conclusion in print: *symbolic simplification
  verification is not extremely effective.* Pipeline: symbolic first, numeric second, 300 s timeout on the
  numeric stage only.
- ⚠ **And their failure mode, which is ours to inherit:** 892 cases (~51.1%) landed above threshold, many
  classified **false positives** caused by missing constraint semantics and **branch-cut** differences —
  not by real errors. ⇒ Numeric probing is weakest exactly where special functions and branch cuts live.
  ⭐ **S9–S11 objects are Lagrangians, matrices, determinants and polynomial dispersion relations — the
  rational/algebraic good case where Schwartz–Zippel actually applies.**
- ⭐⭐ **The most demanding symbolic-physics community already abandoned expression comparison.** Multiloop
  amplitude reduction (**FiniteFlow**, finite-field reconstruction, arXiv 1608.01902 / 1905.08019) ⓢ does
  *all* the algebra by **evaluating at random rational points over a finite field**, precisely to sidestep
  intermediate expression swell. Exact arithmetic, no floats.
- **Fingerprinting** ⓢ: evaluate at N random points at high precision, round, hash (MD5) → a compact
  fingerprint; equal fingerprints are treated as equivalent. Used in symbolic regression and, in a
  structural form, as hash consing (arXiv 2509.20534).
- **e-graphs / equality saturation (`egg`)** ⓢ: the modern *symbolic* answer. Real and usable, but you must
  supply the rewrite theory, and it decides equivalence only modulo the rules you give it. Higher cost than
  probing, lower confidence per unit effort for our objects.
- ⛔⛔ **The caution, from the same 1989 paper as §A: never probe in floating point.** The consistent
  comparison problem is *about* finite-precision comparisons producing disagreement between two correct
  versions. ⇒ **Probe at exact rationals or in a finite field** — which is what FiniteFlow does and what
  Schwartz–Zippel assumes.

**On the 11 nullspace DISAGREE rows specifically:** ⭐ that is not an equivalence problem at all. It is a
**subspace** comparison performed as a **basis** comparison. Standard, textbook, no negotiation: RREF is a
canonical form for a subspace; or check each engine's basis vectors lie in the other's kernel and that the
ranks match. This class of disagreement disappears. **TOOLING** — those were two correct answers.

**On naming:** fingerprint every emitted object in both engines and **join on the fingerprint, not the
tag**. Content-addressed matching. Names become a human label, not a join key. **TOOLING**, and it is the
single change with the largest measured cost saving available.

⚠ **Honest limits.** Probing is one-sided (it *proves different*; it can only *fail to disprove same*); it
requires agreement on **inputs** (a symbol→value map), which is a far smaller and more checkable agreement
than agreement on output names; and it says nothing about whether the object being compared was *computed*
(see §4).

### E · Formal methods, honestly costed — ⛔ **full formalisation no; result-checking yes**

- Cost is real and large: **de Bruijn factor ≈ 4** typical (5.6 called high) ⓢ; a serious formalisation
  effort cited at **251 days**. mathlib is ~232,000 theorems. **PhysLean / PhysLib** exist but are early.
  ⇒ ⛔ **Full formalisation is out of reach for two people and would not save more than a week.**
- ⭐⭐ **The lightweight subset that is real, and it is not proof — it is CERTIFICATION / RESULT CHECKING.**
  Blum & Kannan program checking; McConnell et al., *Certifying algorithms* (Computer Science Review 2011)
  ⓢ: an algorithm emits a **witness** with each output, and *the checker is much simpler than the
  algorithm — yet it is all the user has to trust.* Applied to CAS directly: *Certifying Differential
  Equation Solutions from Computer Algebra Systems in Isabelle/HOL* (arXiv 2102.02679) ⓕ — **do not verify
  the derivation; substitute the answer back.**
- ⭐ **Costed for us, in hours not weeks, needing no proof assistant:** each engine emits a witness the
  other side (or a five-line script) checks by substitution. *"Here is a root"* → check `det M(root) ≡ 0`.
  *"Here is the nullvector"* → check `M·v ≡ 0`. *"Here is the factored determinant"* → expand and subtract.

**Verdict.** ⭐⭐ **This catches PHYSICS**, it is cheap, and — importantly — **it is single-engine.** A wrong
root fails its own substitution check with no second engine present.

### F · How physics actually validates symbolic derivations — ⚠ **the brief was right that this matters most**

⭐ **The answer is: mostly not by dual-engine.** The named practices are:

- **Limiting cases and known-solution reduction**, **symmetry constraints**, **dimensional analysis**, and
  **high-precision numerical evaluation**. Contemporary and explicit: verifiable-AI-physicist work ⓢ
  predefines exactly three assertion types — *limiting-case tests* (exact diagonalisation at small size),
  *symmetry tests*, and *analytical tests* in regimes with known results.
- ⚠ Dimensional analysis is described in this literature as *powerful but not able to spot all mistakes* —
  a fast check, not a verdict.
- **Where multi-system cross-checking IS used it exists and is respectable**: one HEP workflow ⓢ verifies by
  **triple CAS (SymPy + cadabra2 + FORM)**, **≥100-digit numerical evaluation**, and **Lean 4 for the
  coefficient arithmetic**. FormLink embeds FORM inside Mathematica/FeynCalc.
- ⭐⭐ **But look at the shape.** In every case the redundancy lands on a **small number of final, physically
  meaningful quantities** — an amplitude, a decay width, a coefficient — cross-checked end-to-end, plus
  numeric evaluation. **Not on hundreds or thousands of tagged intermediate emissions.**

**Verdict.** ⭐⭐ **Limits, symmetry and ablation are single-engine, cheap, and catch PHYSICS.** Our existing
control matrix (`XFORM_FULLGRAD`, `XFORM_DIVONLY`) already *is* this — it is a **metamorphic relation set**
in the software-testing sense, which is the recognised answer to the oracle problem. ⭐ We built the right
instrument and then spent the budget elsewhere.

### G · Provenance — **the standard pattern exists and we are on the wrong side of it**

- W3C **PROV / PROV-O**, extensions **ProvONE**, **RO-Crate** (arXiv 2312.07852) ⓢ. The distinction that
  matters: **prospective** provenance (the plan — a hand-written specification of what will run) versus
  **retrospective** provenance (**recorded by the run itself**).
- ⭐ **Our registry, keyed on file path + line range, is hand-maintained prospective provenance.** It breaks
  on every rewrite **by construction** — that is not a bug in our registry, it is the known cost of the
  form.
- ⭐ **Standard practice: the producer emits its own provenance at run time, keyed by a stable entity
  identifier.** We already have the stable identifier — **the tag**.

**Verdict.** ⛔ Do not adopt PROV-O/RDF (weeks of ontology work, buys nothing here). ⭐ Adopt the *idea* in
our existing YAML: attach provenance to the emitted record; delete the line ranges. **TOOLING**, low cost,
removes recurring rework.

### H · Agent-driven symbolic derivation with automated verification — **thin in 2024, active in 2026**

Three independent 2026 groups, converging:

- ⭐⭐ **Agentic Diagrammatica** (Alabama/Fermilab, arXiv 2603.26990) ⓕ — **publishes our exact diagnosis**:
  *correctness is governed by implicit mathematical conventions that are not encoded in a form that can be
  easily checked in the computational backend.* Their remedy is **not** two engines. It is
  **tool-constrained computation**: the agent emits a *compact, human-auditable specification*, and a
  **trusted backend performs the manipulations exactly**.
- **Grounded autonomous research** (arXiv 2607.02329) ⓕ — the operative grounding mechanism is
  **numerical confrontation at calibration checkpoints**: reproduce a known published result before
  attempting a novel computation.
- **N-Version Programming with Coding Agents** (arXiv 2606.20158) ⓕ — §A above.

**Verdict.** ⭐ **The convergent 2026 answer is: narrow trusted backend + numeric confrontation, not N
versions of the whole derivation.** ⚠ It is a young literature and none of these papers is about a ledger
of the kind we keep. ⛔ Treat as oracle, not premise.

---

## 2 · The three direct answers

**Q1 · Is dual-engine recognised, an anti-pattern, or unusual? Does anything predict our failure?**
⭐ **Recognised** — N-version programming / differential testing, and as of 2026 actively studied with
coding agents specifically. ⛔ Not an anti-pattern, but its economics are *contested* (Hatton vs. his
critics), and its central assumption — independent development gives independent failure — was **measured
false in 1986 and measured false again for LLM agents in 2026**, with the shared specification named as the
correlation channel both times.
⭐⭐ **Yes, our exact failure is predicted, and by name**: *the consistent comparison problem* (Brilliant,
Knight & Leveson, 1986/1989) — an N-version system unable to reach consensus **when no version has failed**.
⇒ ⭐ **The literature locates our cost where we measured it: in the comparator and the spec, not in the
physics.** That is a strong external confirmation that we are not doing something idiosyncratically wrong;
it is also a warning that no amount of spec care will close it.

**Q2 · Is there a standard way to compare two CAS results without agreeing on names?**
⭐⭐ **Yes. Evaluate both sides at a shared set of random EXACT points (rationals or a finite field) and
compare values** — Schwartz–Zippel identity testing, operationalised as *fingerprint and join*: hash each
engine's evaluation vector, match on the hash, and let each engine name things however it likes.
**What it takes to adopt:** (i) agree a **symbol→value map** per step — this is the only agreement left, and
it is small, mechanical, and checkable, unlike agreement on output names; (ii) each engine emits a
numerically evaluable form (both already emit parseable expressions); (iii) compare subspaces as subspaces
(RREF or kernel-membership), not as bases. ⚠ Probing is one-sided and weak at branch cuts — our objects are
the good case, but this must be stated in the record, not assumed away.
**Precedent:** Rubi and the DLMF suites both do symbolic-then-numeric; HEP's hardest calculations do
finite-field evaluation *only*.

**Q3 · Is there a cheaper validation strategy of comparable confidence?**
⭐ **Yes, and there are two, both single-engine.**
**(a) Result checking / witnesses** (§E): substitute the answer back. Catches wrong physics; needs no second
engine; costs hours.
**(b) Metamorphic relations** (§F): limits, symmetry, dimensional reduction, and **form ablations** — which
we already build.
⇒ ⭐⭐ **The defensible reduction is by GRANULARITY, not by dropping the second engine:** keep dual-engine on
the **handful of headline objects per step** (the action, the system matrix, the determinant, the roots),
and stop cross-comparing the thousands of intermediates. ⚠ **Evidence from our own S11 session supports
this more strongly than anything in the literature:** an independent reviewer reconstructed the matrices,
determinant and roots from the supplied action and reproduced both engines' values — **one human-scale
re-derivation of the headline objects delivered the verdict that 3,750 tags and 25 defect-fixes did not.**

---

## 3 · Comparator architecture — universal harness vs. per-step

⭐⭐ **The brief poses the wrong axis, and prior art shows why.** Rubi and the DLMF suites run **one universal
harness** across tens of thousands of problems with no per-problem configuration — but they can only do that
because their comparison predicate is **semantic and universal** (evaluate; compare to a reference), so
there is nothing per-problem to configure.

⇒ **The variable that decides is not universal-vs-per-step. It is NOMINAL-vs-SEMANTIC comparison.**

- With **nominal** matching, per-step is correct — the naming negotiation is irreducibly per-step, which is
  why `checks_S10.yaml` reached 3,121 lines and still covered a subset.
- With **semantic** (value-based) matching, one universal harness is cheap, and the per-step file collapses
  to three things: **the symbol→value map, the probe point set, and the list of headline objects.** Tens of
  lines, not thousands.

⭐ **Recommendation: choose semantic first; the harness question then answers itself.** ⚠ So prior art
**disagrees** with the brief's stated leaning toward per-step comparators — but only conditionally, and the
condition is the one thing we would change anyway.

---

## 4 · ⭐⭐ Where we are solving a problem that does not need solving

1. **The 11 nullspace-basis DISAGREE rows.** Not a disagreement. A subspace compared as a basis. Textbook
   fix, zero negotiation.
2. **Pinning tag spellings in prose.** A **join-key problem being solved in natural language** — three
   rounds, two of which reintroduced divergence while closing it. Content-address it and the class is gone.
3. **Comparing 3,750 intermediate emissions.** No community in this survey compares at that granularity.
   The physics content of a step lives in a handful of objects; the rest is the engine narrating itself.

⛔⛔ **THE COUNTERWEIGHT, AND IT IS THE MOST IMPORTANT PARAGRAPH HERE — do not adopt §4 without it.**
The 3,750 tags exist because of a *measured* defect: engines emitting physics conclusions as **typed
sentences with no computation behind them**, missed by eight review legs. ⛔ **Value-based comparison does
not fix that. A typed literal evaluates to a number perfectly well, and will fingerprint-match a typed
literal in the other engine.** Numeric probing solves the **naming** problem and the **representation**
problem. It does **not** solve the **fabrication** problem.

⭐⭐ **What solves fabrication is mutation, and we already measured it: of engine 1's nine defects, NONE was
visible by reading and EVERY one was visible by mutation; the ablating leg found 8 of 9.** That is rule 14
with a number attached, and it is consistent with everything in §F. ⇒ **If any control is expanded, expand
that one.** ⛔ If a value-based comparator were adopted as a *replacement* for tag inventory rather than for
tag *matching*, it would re-open the exact hole the rebuild exists to close.

---

## 5 · What I could not find, and where I looked

Searched: ACM/IEEE via web index, arXiv (fetched full abstracts and one full text), the OpenMath standard
and CD registry, SymPy documentation and mailing-list threads, the Rubi and 12000.org test-suite sites, and
the xAct/Cadabra/FeynCalc/FORM documentation and papers.

- ⛔ **No published methodology for dual-CAS comparison of *intermediate tagged emissions* at scale.**
  Every corpus found compares final results against a stored reference. ⚠ **This is absence of evidence,
  not evidence of originality** — the most likely explanation is that nobody needs to, because the
  granularity is a choice we made and they did not.
- ⛔ **No maintained OpenMath↔SymPy or OpenMath↔Wolfram round-trip.**
- ⚠ **Abbasi's exact equivalence procedure could not be retrieved** — 12000.org returned HTTP 403. The
  description of its symbolic/numeric grading here is **secondary**; open it directly before relying on it.
- ⛔ **No cost figures for adopting metamorphic testing in a two-person project.** The MT literature is
  large but reports bug counts, not adoption cost.
- ⛔ **No prior art on "how many objects per derivation step should be cross-checked."** The §2 Q3
  recommendation rests on our own S11 measurement, not on a published result.
- ⚠ **Papers marked ⓢ were not read in full.** In particular the Knight & Leveson and Brilliant/Knight/
  Leveson primary sources, the MAST percentages, and the FiniteFlow methodology are second-hand here.
  ⛔ Anything load-bearing should be read before it becomes a premise.

⛔⛔ **Rule 16 applies to every line above: this is an oracle, never a premise.** The DLMF conditions are
special functions with branch cuts, not our rational algebra; the N-version conditions are safety-critical
software with a reference implementation, which we do not have; the HEP conditions have published results to
confront, which for our object do not exist. ⭐ **Check our computed result against these; never assume
their result holds for our object.**
