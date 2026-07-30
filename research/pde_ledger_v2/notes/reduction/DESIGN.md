# The reduction workstream — design, after two independent reviews

**Question:** how many independent physical inputs does the model actually take?
**Why it matters:** the calibrate-predict surplus (held-out matches − tuned knobs) has no denominator
until this is known.

Reviewed independently by Codex (xhigh) and Grok-4.5, neither seeing the other's answer. ⭐ **They
converged on every major point below.** Where only one raised something, it is attributed.

---

## 1. ⛔ The premise this workstream started from was wrong, and the correction is the finding

**Started from:** *"stage043's [40,49] is an inventory of knobs; nothing composes `c_s = f(P,ρ,m)` and
eliminates intermediates, so a real reduction should collapse it to single digits."*

**Both reviewers rejected this, and it checks out against the corpus:**

- stage043's counting rule already subtracts executed reductions —
  `count = nominal − DERIVED-and-EXECUTED − CONV − external-bridge`
  (`notes/stages/ledger_stage043_irreducible_count_range.md:69`, reinforced `:78`).
- The composition machinery **exists and runs**:
  `scripts/midway_knob_audit_codimension_sympy.py` computes Krull dimension from the leading-monomial
  ideal of a grevlex Gröbner basis (`:28-50`), certifies it with an exact **positive smooth real
  witness** plus a **Jacobian corank** check (`:53-73`, `assert len(variables) - rank == dimension`),
  and asserts an `EXPECTED_PAYLOAD` (`:161`, `:185`).
- The relations it composes are declared `DERIVED` with explicit collapse effects:
  `R1` `c_s0 → {K, ρ0, m_GNLS}`, `R2` `ξ_h, a, h0 → primitives`, `R3` `λγ` not independent,
  `R20` `δ, σ_wall → {a_B, κ_B}`, `R29` `K_η → {T_w, β}` (`notes/parameter_register.md:268-296`).

⇒ **The easy algebraic chain has already been reduced.** The residue is not "variables nobody wrote an
equation for."

### What the residue actually is (Grok's decomposition; ⚠ relayed, not independently verified)

| | approx |
|---|---|
| medium primitives `{ħ, m, K, ρ0}` | **~4** |
| route-ful reduction debt, still counted | ~15–21 |
| route-less postulates / structure | ~14–20 |
| force-magnitude + Parts III–VI continuous | ~13 |
| open convention C1/C2 | ±9 |

⭐ **Both halves of the owner's intuition are right about different things.** The medium primitives
really are a handful. The ~40 is primitives **plus** un-executed debt **plus** route-less postulates
**plus** sectors outside the Path-A constant set. Many of the debt items *do* have equation-shaped prose
— `R30` says a nonlinear throat solve *would* derive `{μ_η, T_w, β}` — and stay free because the
equation **is not executed**, not because nobody composed it.

⚠ Grok's diagnosis of where "it can't be 40" comes from: `decisions/14`'s ~8 genuine free parameters is
a **narrower scope**. Comparing it to `[40,49]` confuses **scope expansion** with **failed reduction**.

---

## 2. ⭐⭐ The universe is the hard part, not the algebra (Codex, and it is the biggest single point)

> *"Every named variable" is the wrong universe.*

It mixes action/constitutive parameters · fields · coordinates · initial and boundary data · derived
observables · intermediate expressions · numerical controls · discrete branch choices · functions and
profiles.

⛔ **If fields and arbitrary initial data are admitted, the input space is infinite-dimensional.** If
they are excluded, there must be a *precise stated rule* for why.

⚠ **Companion failure — universe incompleteness:** *"the most impeccable rank computation still gives a
falsely small number if an entire parameter, profile or boundary condition never entered the
registry."* The register is measured at **105 of 226** distinct dimension-bearing quantities
(`notes/measure_register_sufficiency.md:27-29`). ⛔ Carry that figure with its own limit: the two
scripts that file names as producing it are not retained (`:142-145`), so it is not re-derivable.
Completeness cannot be assumed either way.

---

## 3. The method both reviewers converged on

⛔ **Not** a definition digraph with substitution. **A typed constraint system with an admission gate,
and dimension by rank.**

Graph reduction is the *triangular special case* — valid only where admitted relations are uniquely
solvable for one head and acyclic. Keep it for human explanation; **rank is authoritative**.

1. **Freeze the counting contract** — model version, regimes, sectors, whether state/boundary data
   count, how functional inputs are represented.
2. **Canonical quantity IDs** — kind, scope, regime, live/retired, dimension convention, aliases.
   ⛔ Same spelling is not identity; different spelling is not non-identity.
3. **Every candidate relation as typed data** — residual expression, designated output, input QIDs and
   applied functions, provenance status, regime and assumptions, domain/measure/BCs/denominator guards,
   source and execution loci, external benchmark refs.
4. **Admission gate.** Only `DERIVED-EXECUTED` enters the earned constraint set. A relation **fails** if
   it is prose not a parseable expression · references an unknown QID · applies a function lacking a
   definition or an explicit functional-input classification · its output recurs without implicit-block
   treatment · it is pending / calibrated / uncommitted / conventional / merely structural ·
   ⭐ **the implementation independently freezes the alleged output rather than computing it** · its
   regime assumptions are incompatible with the current branch.
5. **Alias-canonicalise, then a bipartite factor graph** (variables | constraints), computed by
   connected component. A separate directed graph records dataflow and SCCs.
6. **Substitute acyclic explicit definitions** — an optimisation only. Adding a derived symbol plus its
   definition must leave total dimension unchanged.
7. **Rank each finite scalar block** — maximum bipartite matching for a structural upper bound; symbolic
   Jacobian at exact constraint-satisfying generic witnesses; a nonzero minor reaching the structural
   bound is a cheap exact generic-rank certificate; Gröbner/Krull for small polynomial blocks, singular
   cases, or where the Jacobian is inconclusive. ⛔ Do not fake-polynomialise rational, transcendental,
   integral or PDE blocks.
8. **Quotient conventions; compute calibration separately.** Report per regime: continuous model
   dimension · convention-orbit dimension · calibrated rank · residual post-calibration dimension ·
   discrete choices · functional inputs / unresolved functional dimension · unresolved candidates.

**Report a bracket, never a scalar:** `[N_optimistic, N_certified]`. ⭐ **A pending relation may lower
the optimistic endpoint; it must never lower the certified headline.**

### Traps in the algebra, each of which silently returns a wrong dimension

- **Singular / non-reduced constraints** — a Jacobian at a singular point gives the wrong dimension;
  `x² = 0` is the standard counterexample.
- **Denominator clearing** adds spurious branches unless nonzero conditions are retained or the ideal is
  saturated.
- **Component choice** — different connected components can have different dimensions.
- **Existence and uniqueness** — a PDE "defines" a profile only with its domain, boundary/initial data,
  gauge fixing **and** a uniqueness/nullspace result. ⭐ A nullspace coefficient or branch selector is
  **another input**.
- **Cycles are not auto-zero.** `dim = |vars_SCC| − rank(J)`. `A=B, B=A` → **1**. `A=2B, B=2A` → **0**.
  An inconsistent cycle → **empty variety**, not "negative inputs".
- **Inactive / observationally redundant parameters** need a *demonstrated reparameterisation symmetry*,
  ⛔ not merely the present inability to compute an observable from them.

---

## 4. ⭐ The flattering-undercount failure mode, and the gate that catches it

Both reviewers named the same one, and the same concrete template.

> A pending or formal equation is admitted as an earned constraint, while its undefined profiles,
> missing boundary data, or calibration anchors disappear from the universe.

**`R35` is the template, not an anecdote.** `M̃ = ∫ μ_η β₂² dV` reads as three integral equations
deleting three scalars — while no functional form of `β₂(w)` or `μ_η(w)` exists in the ledger, and the
scripts still freeze the outputs independently (stage017 carries them as `calibration_inputs`).

⭐⭐ **And the deeper version:** if `β₂(w)` is a free function not fixed by an in-corpus PDE/BVP, the
"reduction" collapses a scalar into a **function space**. That is ⛔ not `dim − 1` — it trades one
unknown for infinitely many unless the function is expanded in a finite basis and its coefficients
counted.

**The mechanical gate — transitive closure:**
1. Parse the RHS.
2. Enumerate every scalar, function application, derivative, operator, measure, domain and boundary datum.
3. Require each leaf to terminate in exactly one of: an admitted primitive scalar input · a fully
   defined finite-parameter function · an explicitly counted functional input · a cited, well-posed
   field-equation solution · a convention or external datum handled on its own axis.
4. **Any unassigned leaf ⇒ the relation cannot reduce the certified count.**
5. **Mutate away each functional definition in turn; the count must stay unchanged or increase, never
   decrease.**

**The standing check on any published number** (Grok): a decrease relative to stage043's low endpoint
requires a named admitted relation, a dual-engine dimension drop, able-to-fail controls, and an explicit
**before/after residue-set diff** a hostile reader can re-run. ⛔ If N drops without that diff, reject
the number.

---

## 5. Acceptance tests that can fail

duplicate/entailed equation → dimension unchanged · genuinely independent equation → drops by one ·
vacuous relation (`x−x=0`) → dimension **increases** · pending/calibration promoted into the earned
ideal → schema failure · undefined profile inserted into an integral → certified dimension **cannot**
decrease · orientation reversed → dimension unchanged, provenance output changes · alias split into two
names plus an identity → physical dimension unchanged · incompatible regimes combined → schema failure ·
`A=B,B=A` → 1 · `A=2B,B=2A` → 0 · inconsistent cycle → empty variety · denominator-zero branch from
clearing fractions → caught · gauge/unit reparameterisation → quotient dimension unchanged ·
independently-frozen output substituted for a computed definition → dataflow failure ·
⭐ **top-down inventory from the assembled action disagrees with bottom-up stage dependencies →
completeness failure.**

⚠ The existing medium and wall blocks are usable fixtures (dimensions **5** and **7**), but ⛔ passing
them cannot certify the untested blocks.

---

## 6. The 96-edge relations ledger: triage queue, ⛔ not a constraint file

`notes/parameter_register.md:264` contains fundamentally different objects: real algebraic candidates
(`R1`–`R5`, `R15`–`R16`, `R24`, `R29`) · pending prose routes (`R6`, `R10`, `R12`, `R13`, `R19`, `R30`,
`R33`, `R36`) · calibrations **not in force** (`R7`–`R8`) · closed-negative routes · structural facts
explicitly marked "not a reduction" · count sensitivities and bookkeeping (`R85`–`R97`).

⛔ **"Edge" there does not mean "independent equality constraint."** Feeding all 96 into a digraph and
eliminating every node with an edge over-reduces catastrophically. Each row is promoted **individually**
through the admission gate.

---

## 7. Where it lives — ⚠ both reviewers say NOT inside `ledger_dimensions.py`

⭐ **The owner's goal stands and is achieved; the mechanism moves.** The goal — *one place holding every
input and what it derives from, so the leftovers are visible* — is exactly right. Both reviewers say
that place should be a **dedicated relation registry**, not an enlargement of the dimensions module:

- dimensions are **necessary but not sufficient** — a dimension algebra is not a derivation graph, and
  merging them is a category error (Grok);
- only 7 of 30 dimension-bearing stages import the shared module today, and making every stage import
  the global relation truth **risks circular audits** (Codex).

**Architecture both endorse:**
1. a dedicated machine-readable **QID + relation registry** as the single semantic input;
2. stage artifacts **emitting or citing** relation records;
3. a **standalone whole-model reducer** consuming that registry;
4. **independent rank implementations** over the same normalised inputs;
5. CI checks for stale source loci, unresolved closures, and count mutations.

---

## 8. What this changes about the deliverable

⛔ **The product is not a single integer.** It is:
1. an **admitted-collapse certificate** — block dimensions, deltas, able-to-fail controls;
2. an **equation-shaped-hole list** — `R35`-class, machine-refused rather than reader-noticed;
3. a **residue partition with names** — primitives / route-ful debt / route-less postulate / calibrated
   — reported as a **range**;
4. a **debt-discharge roadmap** — which single future execution (one throat solve, `R0` support, the
   Gate-6 return law) moves which residue members into the earned set.

⭐ **(3) is the same three-way split the census was already asked to produce** (debt · structural ·
postulate), arrived at independently from the algebraic side. That convergence is the strongest
argument that the partition — not the integer — is the real deliverable.

⚠ **Will the certified number be single digits soon?** Both say almost certainly not, unless the
deferred throat/return/support solves are discharged, or "input" is redefined to exclude debt and
calibration — which ⛔ breaks surplus honesty, because debt knobs are free or tuned until discharged.
