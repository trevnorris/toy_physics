# THE DERIVATION WALKTHROUGH — plan

**Status: IN USE.** User decision 2026-07-30 to change method. Reviewed independently by Codex and Grok;
their findings are folded in (§1.0, §3's retraction, §7.4's cut-back claim, §7a, checks 10–13).
**Phase 0 is running — two steps banked**, `research/pde_ledger_v2/walkthrough/`.

⭐ **Read alongside this file: `research/pde_ledger_v2/walkthrough/DECISIONS.md`** — the user calls that
shape what gets counted. ⛔ A decision recorded there can **overturn a ruling recorded here**, and one
already has (**D-01**, the `a`-pin — see §1.1 and §8). Where the two disagree, `DECISIONS.md` is the later
document and governs; this file is not self-sufficient on any question it touches.

---

## 0. Why the method changed

The audit route worked **backward**: take 43 stages of finished artifacts and infer what was derived.
It produced real findings but at a cost that made the physics unfollowable — and the physics has exactly
one reviewer.

⚠ **This is the second occurrence of the same failure.** `docs/development_pipeline.md` already records
apparatus growing above the physics once, and the user stopping it. This time the shape was: census
schema → two review legs → pilot → blocking review → fix → re-run → blocking review → registry. Eleven
commits, and **no physics verified**. Every genuinely useful output was a finding about a *specific
artifact*, obtained by opening a file.

⭐ **The walkthrough works forward.** Start at the medium's defining properties, take one derivation step
at a time, and record at each: **what it is · what it does · what's new.** The irreducible input count
then falls out **by construction** — "what's new" accumulated down the chain — rather than being
inferred backward from artifacts. It is also followable, which the audit route was not.

⚠ Convergent evidence: Codex and Grok, reviewing the reduction design independently, both concluded the
deliverable is a **residue partition plus a debt-discharge roadmap** — *which future execution moves
which member into the earned set*. That is a forward-order work list. They reached the user's framing
from the algebra side.

---

## 1. What each step records

One step = one derivation move. Its record is short and fixed:

| field | meaning |
|---|---|
| **what it is** | the statement, in one or two sentences |
| **what it does** | what it buys — what becomes computable that was not |
| ⭐ **what's new** | every quantity or assumption entering here that did not exist upstream. ⛔ **The introduction INVENTORY — not the count.** Only §7a certifies a count. |
| **inputs** | each upstream quantity consumed, with the step that produced it |
| **the equation(s)** | canonical form, in the registry (`prefix-v1`), evaluable |
| **class per new item** | `derived` · `calibrated` · `postulated` · `debt` (a named-but-unexecuted route) |
| **regime** | the assumptions under which it holds |
| **departure** | where this differs from standard GR/EM, if it does |

⭐ **"What's new" is the introduction inventory.** A step that introduces nothing new is pure consequence. A step
that introduces a postulate is tier-1 material, recorded where it enters rather than reconstructed later.

⛔ **But "what's new" is NOT the count.** See §1.0 — it is the readable *introduction inventory*, and it
requires a closing certification step before any number is quoted.

### 1.0 ⭐⭐ What the count IS (user definition, 2026-07-30)

> **How many variables must I define for the simulation to run?** That is the minimal set.

⭐ This is concrete, checkable, and it absorbs the three cases the algebraic framing handles badly:
- a **free function** (`β₂(w)`) ⇒ you must supply a *profile* — one entry, not one scalar;
- a **boundary condition / ensemble class** ⇒ you must supply a *BC*;
- a **nullspace coefficient or branch selector** ⇒ you must *pick one*.

⇒ Report the set **partitioned by kind**: scalars · profiles/functions · BCs and domains · discrete
branch choices. ⛔ Never as one integer, because the kinds are not interchangeable.

**Its relationship to the algebraic residual dimension — they are a matched pair, not rivals:**
- the **sim-input set** is what you must *supply*;
- the **residual dimension after admitted constraints** is the algebraic *minimum*;
- ⭐ **where they disagree, that is a finding.** Either something being supplied is actually derivable,
  or a constraint is hidden. ⚠ Real template: `R9`, where two locks read as one story until codimension
  shows `Δr = 2`. Introduction accounting alone cannot discover independence; only rank can.

⇒ ⛔ **The forward walkthrough does not by itself certify the count.** It produces the inventory and
makes it followable. Certification needs the closing step in §7a.

### 1.1 ⭐⭐ The classification test (user decision, 2026-07-30)

Apply in order. It is mechanical, which is why it is usable at every step:

1. **Does it have a defining equation — an expression in terms of other model quantities?**
   → **`derived`**. ⚠ This holds *even where the quantity is later used as a unit*: `c_s0` is one of the
   four natural-unit pins, yet `c_s² = nKρⁿ⁻¹/m` defines it from `K`, `ρ0`, `m` and the EOS exponent
   **independently of any unit choice** — so `c_s0` is `derived`, and being set to 1 in a pin does not
   change that.
   ⭐ **Clause (from D-01, and it OUTLIVES the case that produced it):** *a relation arising from
   imposing unit pins is not a defining equation* — the residue of a units choice expresses the choice,
   not the model's content. ✅ `a = ħ/(m c_s0)` was that case; it is **RESOLVED BY REMOVAL**
   (`407eed94`) — the relation and the quantity are deleted, and `a` now means the throat radius. ⛔ There
   is no open class. The clause stands for **future** pin-shaped relations (§8).
2. **Was it chosen because the calibration was necessary for the model to work?** (`n = 5`, `β = 3`)
   → **`calibrated`**.
3. **Is there a named route to a defining equation that nobody has executed?** → **`debt`**.
4. **None of the above** → **`postulated`** — genuinely irreducible, tier-1 core.

⛔ **`convention` is NOT a class here.** It survives only as the narrow case of an **assignment to a
number** that carries no content — `c = 1`, `ℓ_P = 1` — which has no defining equation in terms of model
quantities and so never reaches test 1.

✅ **CLOSED 2026-07-31 by removal.** This paragraph once argued that `a` *"does have a defining
equation"*, then D-01 reopened it. Both readings are moot: the user retired the pin, `R2.a_pin` and
`Q.medium.a_pin` are **deleted**, and no class needs picking. ⛔ Do not re-open it. ⚠ The durable lesson
is the **kind-test**, not the verdict: ask *"is this one number for the whole medium, or one per
particle?"* — `ħ/(m c_s0)` **is** a healing length in a standard convention, so a value-check passes it
and only a kind-check catches it.

⇒ The four classes map straight onto the tiers: `derived` → **tier 3** · `calibrated` → **tier 2** ·
`postulated` and `debt` → **tier 1**, split as the user specified.

---

## 2. Checks per step

⛔ **Governing test applies to every check below**: does it catch a way the **physics** could be wrong?
If it only catches tooling or a motivated adversary, it is not here.

### Mandatory — every step
1. **Dimensional homogeneity.** Every equation checked from declared dimensions.
   *Catches:* a dropped Jacobian, a wrong power, a volume-vs-line confusion.
   **Built and working** (`reduction/dimensional_homogeneity_gate.py`, demonstrated able to fail).
2. **What's-new accounting.** Each new quantity classified `derived`/`calibrated`/`postulated`/`debt`.
   *Catches:* a fit entering disguised as a derivation.
3. **Input provenance.** Each input cites the step that produced it, by **locus**, not by name.
   *Catches:* the measured corpus-wide defect — **zero cross-artifact citations resolve to a locus**.
4. **Step-to-step identity.** Inputs must be the same *quantity* as the prior step's output, adjudicated,
   not name-matched. *Catches:* the same-name-different-quantity merge, and its converse.
   ⚠ Real here: eight known-same cross-stage pairs that name-matching misses.

### Conditional — only when the step has the feature
5. **Dual-engine, with teeth** — when the step contains real algebra. One canonical expression, ⛔ never
   two hand-typed copies. ⚠ **But parsing the same tree twice only tests transport tooling.** Each engine
   must **independently derive or evaluate the result from upstream primitives**; two result-emitters
   agreeing on a literal is vacuous.
6. **Able-to-fail with a PLAUSIBLE PHYSICS MUTATION** — when the step asserts something. ⚠ "Perturb
   something and an assertion fires" proves liveness, not correctness. The mutation must be a real
   physics error: a wrong-sign source, an omitted term, a wrong branch, a singular denominator, a bad
   boundary condition. ⛔ Mutating a label or a report string does not count.
7. **Existence / uniqueness** — when a PDE or BVP appears. A profile is defined only with its domain,
   boundary data, gauge fixing **and** a uniqueness result. ⭐ A nullspace coefficient or branch selector
   is **another input** and must appear under "what's new".
8. **Held-out vs calibration-consuming** — when a number meets data. Tag which. This is the surplus
   ledger, kept forward instead of reconstructed.
9. **Regime compatibility** — when an assumption is introduced: does it conflict with one already live?
10. ⭐ **Term-by-term fidelity** — when a step derives from an action or a balance law. Derive the
    Euler–Lagrange equation, its source, sign and coefficients from the upstream object **term by term**.
    ⚠ Dimensions and two parsers ⛔ cannot catch a *faithful transcription of the wrong operator*.
11. ⭐ **Dynamical health** — when a dynamical block appears: conservation/balance laws, physical DOF and
    constraint count, stability/ghost check, characteristic/hyperbolic structure. ⚠ Measured relevance:
    both current drain candidates are **dimensionally valid yet different physics** — dimensions alone
    cannot separate them.
12. **Approximation validity**, not merely regime compatibility — for any reduction or truncation, show
    which discarded terms vanish or are controlled, and what error order remains. ⚠ "No assumption
    conflicts" is strictly weaker than "the approximation is valid".
13. ⭐ **A second independent physics review leg** on physics-bearing steps. ⛔ Do **not** cut this as
    ceremony: `docs/development_pipeline.md` records a tautological "decisive test" that **only a third
    independent reader caught**. Orchestrator + user alone re-opens that hole.

### Considered and rejected
Byte-level artifact custody · freeze authorities · citation-integrity tooling beyond check 3 · a manifest
per step. ⛔ All fail the governing test — they catch tooling or an adversary, not a wrong derivation.

---

## 3. Step order

⚠ **Dependency order, not the ledger's part numbering** — the existing order is a build history.
⛔ The physical sequence is the user's call; this is a proposal to correct.

| phase | content | expected "what's new" |
|---|---|---|
| **0** | **Brane + bulk defining properties** — what the medium *is*, before any derivation | ⭐ almost everything; this is tier 1's expected core |
| **1** | **Medium mechanics** — EoS, sound speed, length scales | few; mostly consequence |
| **2** | **Excitations** — brane shear, signal speed (light) | the light-sector anchors |
| **3** | **Defects** — the ±w puncture, its topology and conservation | the charge-sector postulates |
| **4** | **Flow** — the drain, Newtonian limit, the PN ladder | calibration meets data here |
| **5** | **Motion** — the moving throat (magnetism) | boost-consistency assumptions |
| **6** | **Knit + integration** — cross-sector consistency, the surplus tally | ⭐ new cross-sector consequences EXPECTED; ⛔ no new input that revises an earlier sector — see the rule below |

⛔ **RETRACTED — "phase 6 introduces nothing new" was already false when written.** `stage044` states
outright that it is **not** a "no new knob" stage and introduces `Z_χ`
(`notes/stages/ledger_stage044_parent_action_reconciliation.md:30`, `:358`). ⚠ It is also **gameable in
the other direction**: declare every future coupling in phase 0 and the test passes by construction.

⭐ **Replacement (corrected 2026-07-31 — the first replacement was still too strong):**

> ⛔ **The knit may not introduce a new *input*** — action, constitutive, source or BC — **that revises what
> an earlier sector already derived.** Such an input falsifies **completeness of the proposed substrate**
> — ⛔ not, by itself, the one-medium hypothesis.
>
> ⭐ **It is expected to produce new *consequences*.** Cross-sector combination yielding constraints that no
> single sector gives — ⭐ including constraints on the throat interior — is the knit's **purpose**, not a
> violation. **C10 is the worked example**: the inertial sector wants a compact throat, the spin sector
> wants an extended one, and neither produces that tension alone.

⚠ New inputs at integration admit several readings — an incomplete phase-0 substrate, an unchosen
inertial-vs-dissipative wall branch, an omitted cross-sector coupling, or two historical descriptions
that are not dynamically identical. ⛔ Do not collapse them to "different media". The corpus itself says
the knit establishes **bookkeeping/labelling** closure, not a dynamical one-medium proof
(`docs/model_map.md:134`, `:139`).

⚠ **v3 deliberately diverges from this ordering: it runs gravity BEFORE charge.**
`research/pde_ledger_v3/V3_STEP_PLAN.md` puts PHASE 4 (gravity) ahead of PHASE 4b (charge/magnetism,
`Q1`–`Q7`) — gravity is the more solid sector and the user prioritized it, which is the call ⛔ *"The
physical sequence is the user's call; this is a proposal to correct"* (above) reserves. ⇒ ⭐ A **RECORDED**
divergence, ⛔ not a contradiction.

---

## 4. What carries forward from the audit work

⛔ **Extract before archiving.** These are findings, not apparatus, and must survive:

- **Zero cross-artifact citations resolve to a locus** — measured on both pilot stages, both engines.
  The single highest-leverage mechanical fix available, and check 3 above enforces it going forward.
- **A false provenance attribution** — stage016's engines assert `M̃`/`K̃`/`T̃_Ω` are
  `CONSUMED-from-011/012/013`; those stages contain none of those symbols.
- **A live contradiction between tracked documents** — `parameter_register.md:170` calls
  `K0c`/`K_eta`/`T_Omega` `FREE-UNREDUCED`, PENDING, counted debt, *"NOT identified with the raw 013/017
  densities"*; `stage023_pathA34_..._source_map.md:250-253` states that identification as performed and
  calls them *"likely DERIVED manifestations"*. ⭐ A tier-1-vs-tier-3 disagreement, unresolved.
- **A wrong locus in four tracked files** — stage016's dimension literals are at `:314-325`; `:355-366`
  is cited in `parameter_register.md:182/:183/:184`, the stage016 note `:194`,
  `rewrite_reference_table.md:205` and `measure_register_sufficiency.md:100` (en-dash spelling).
- **An off-locus Jacobian understated rank** — the helper reported 4 where the exact satisfying-witness
  rank is 5. ✅ **FIXED 2026-07-30 (`6e79e2ec`)** — rank is now evaluated at exact constraint-satisfying
  witnesses, `acceptance_check.py` returns `MATCH` exit 0 on all four medium cases, and two independent
  reviewers confirmed the four numbers. ⛔ **No longer a live latent defect.** ⭐ It stays on this list
  because the *finding* carries forward: an off-locus Jacobian biases **flattering** — a too-low
  `dim_after` makes the model look more determined than it is.
- **The easy algebraic reduction is already done** — `midway_knob_audit_codimension_sympy.py` composes
  and certifies block dimensions. ⛔ The residue is *not* uncomposed inventory.

**Kept as live infrastructure, not archived:**
`reduction/` (registry, reader, homogeneity gate, able-to-fail harness) — it is the recording surface for
§1 and already evaluable for the propagation check · `notes/reduction/DESIGN.md` · the dimension rewrite
and its canonical doc.

---

## 5. Cleanup — ✅ DONE 2026-07-30

Archived by `git mv` (history preserved, nothing deleted): `archive/census/` and
`archive/manifests/`. ⭐ `DIMENSION_REWRITE.md` and `DIM_ORDER_DECISION.md` deliberately **kept** at
`research/pde_ledger_v2/manifests/` with unchanged paths; `notes/ablation_driver/` and
`notes/wl_emitter/` kept as the specs for checks 6 and 5.

⇒ **What was archived, why, and where the census's findings went: `archive/README.md`.**

## 5a. ⭐⭐ The two scripts — a standing deliverable, grown step by step

⭐ **These remain the end product** (user, 2026-07-30), and they are ⛔ **not** built up front. They grow
as the walkthrough proceeds, so kinks surface on real steps instead of being designed around in advance.

1. **The one script every other SymPy script imports** — the key to all variables: each quantity, its
   defining equation if it has one, and what it derives from. ⭐ This is what makes the leftovers
   *readable* rather than reconstructed — the original reason the shared import was asked for.
   Today: `research/pde_ledger_v2/reduction/` (`quantities.yaml`, `relations.yaml`, `registry_read.py`).
2. **The dimensions-check script** — verifies every equation is dimensionally homogeneous from those
   declared dimensions. Today: `reduction/dimensional_homogeneity_gate.py`, working and demonstrated
   able to fail.

**How they grow:** each walkthrough step adds its quantities and its equation to (1), and must pass (2)
before the step is banked. ⇒ By the last step, (1) *is* the grand check and (2) has been exercised on
every equation in it — rather than both being written against a corpus nobody walked.

⚠ Both must end up readable by **SymPy and Mathematica**. The `prefix-v1` form already stores one
canonical expression both engines parse; ⛔ never two hand-typed copies.

---

## 6. Working discipline

- **One step at a time**, with a stop for the user at each. ⛔ No fan-out across steps.
- ⛔ **No agent parallelism** unless a single step genuinely needs it. The pace problem was caused by
  running three and four agents and delivering dense synthesis.
- **The step record is the deliverable**, not a report about the step.
- Commit per step. Commit before anything destructive.
- ⭐ **If a step cannot be derived, that is the result.** Record it as `debt` or `postulated` with what is
  missing, and move on. ⛔ Do not grind.
- ⛔⛔ **THE TRAP THIS PLAN IS MOST LIKELY TO FALL INTO, named by review:** building the registry schema
  and its CI as a *precondition* for banking a step, instead of adding one quantity and one residual when
  a step needs them. ⚠ **If execution recreates census-grade apparatus before phase-1 physics, it fails
  the same governing test that killed the audit route.** §7a is a **closing** step, ⛔ not a gate to build
  first.

---

## 7. Acceptance

The walkthrough is working if, at any point:
1. the accumulated **"what's new"** list is the irreducible input count, with each member traceable to
   the step that introduced it;
2. every step passes its mandatory checks, with failures **recorded, not fixed by adjustment**;
3. a reader can follow the chain from phase 0 to any step without consulting the audit apparatus;
4. the registry can take values for the sim-input set and propagate them forward — the
   **explicit-scalar dataflow smoke test**.

⚠ **(4) renamed, and its claim cut back.** It is ⛔ **not** "the real test of the whole program". A green
run proves, for the traversed explicit branch: every required scalar leaf was supplied or recursively
computed · denominators and declared assumptions held · derived values were **recomputed rather than
independently frozen** · the stored equations are internally evaluable. That is real and worth having.

⛔ **It does not prove:** registry completeness · correctness of the stored physics · PDE/BVP existence
or uniqueness · implicit or cyclic block solvability · conservation, stability or constraint consistency
· correct calibration/held-out separation · agreement with data · **the irreducible count**.

⚠ Report PDE solves and functional inputs **separately** — not every sim-input member has a numeric
value, and law-form postulates, discrete choices and profiles do not fit that interface.

---

## 7a. ⭐⭐ The closing certification — required before any count is quoted

The walkthrough produces the **inventory**. It does ⛔ not certify it. Before a number leaves this
project:

1. **Admission gate** — only a `derived` relation whose output is actually computed (⛔ not
   independently frozen elsewhere) may reduce the count.
2. **Transitive leaf closure** — enumerate every scalar, function, derivative, operator, measure, domain
   and boundary datum on the RHS; each must terminate in a supplied input or a well-posed computed
   solution. ⛔ Any unassigned leaf ⇒ the relation does not reduce the count. ⭐ This is the `R35` gate.
3. **Block rank** — because introduction accounting cannot discover independence (`R9`, `Δr = 2`).
4. **Top-down reconciliation** — map the assembled action, its non-variational sources, BCs, domains and
   profiles back onto the forward steps. ⭐ This is the only check that catches a **universe hole**: a
   parameter that no step ever named.
5. **Sim-input vs residual-dimension diff** (§1.0) — report it; ⛔ do not reconcile it silently.

⇒ Report as a **range with a residue partition**, never a scalar.

---

## 8. Open for the user

1. **The phase order in §3** — proposed on dependency grounds; the physical picture is the user's.
2. ✅ **CLOSED 2026-07-31 — the `R2` `a`-pin question no longer exists.** It was reopened by D-01, then
   **resolved by removal**: the user retired the pin entirely (commit `407eed94`). `R2.a_pin` and
   `Q.medium.a_pin` are **deleted** from the registry; `a` now denotes the **throat radius**, a
   defect-sector quantity with no defining equation, absent from the medium block until a defect step
   introduces it. ⛔ There is no open class to pick. ⇒ `research/pde_ledger_v3/` is the live workstream.


✅ **Resolved 2026-07-30:** the archive items (both kept, §5) · Zenodo (user-owned, no constraint).
