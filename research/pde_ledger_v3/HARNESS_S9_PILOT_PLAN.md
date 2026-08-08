# S9 pilot — phase 1 smoke test for the harness standard

**Status: PROPOSAL, 2026-08-06. ⛔ Nothing built. ⛔ Not reviewed — this plan needs its own two legs
(this revision was applied by Codex from an orchestrator decision list ⇒ fresh Claude agent + Grok)
before any of it is built.**
Evidence base: `docs/method_prior_art_findings.md` **as corrected by its two review legs** (Grok + Codex,
2026-08-06) plus the S9/S10/S11 measurements below.
⚠⚠ **Both legs found material defects in the findings document, and one was a FALSE UNIVERSAL that this
plan's first draft was built on** (§3b). ⛔ Read §3b before §3. ⭐ The corrections are folded in here and
back into the findings document; ⛔ do not work from an uncorrected copy of either.

⛔⛔ **THIS CHANGES NO PHYSICS.** If any computed value differs from the committed output, **STOP and
report it — ⛔ do not adjust either engine or the comparator.**

---

## 1 · Why S9 and not S10 — the measurement

| | **S9** | **S10** | **S11 (in flight)** |
|---|---|---|---|
| `checks_S<n>.yaml` | **194 lines** | **3,121** | — |
| declared cross-engine pairs | **12** | **690** | — |
| tags emitted (py / wl) | 635 / 1,559 | 4,227 / 2,983 | 3,539 / **3,750** |
| engine scripts (py / wl) | 640 / 782 | 1,738 / 1,807 | — |

⭐ **91% of `checks_S10.yaml` (lines 273–3114, 690 entries) is a hand-written name→name pair list.**
⚠ `REBUILD_HANDOFF.md` currently designates **S10** as the template. ⛔ On this evidence that is the wrong
template: S10 is where the apparatus outgrew the physics. **S9 is the Phase 1 smoke-test target.**

⭐ **And S9 has a compact adjudication inventory**: 12 declared pairs, plus
`measurements/declaration_load_ablation.py` has already isolated a load-bearing subset of today's **six
naming exceptions and one separate symbol identity**. ⇒ a rebuild can be **falsified** here, cheaply,
against a pinned, builder-inaccessible differential target that involves **no expected physics value**.

### The two-phase boundary

⭐ **Phase 1 is S9, and S9 alone is a smoke test.** It is the right cheap falsifier for the naming
question and nothing more: it has 12 pairs, one unparsed conditional sign, and
exercises few object kinds.

⛔ **Phase 2 is required before any method adoption or coverage reduction.** It adds representative rows
for an S10 **matrix**, an **action/derivative** object, a **relation**, a **mapping**, a **multiset with
multiplicity**, and a **parameter-dependent subspace**; and for S11b **integral**, **transcendental**,
**pole**, **branch-condition**, **sheet**, and **textual** cases. These cases are live today at
`scripts/S11bA_interface_response_sympy_audit.py:37` (unevaluated integrals), `:336`
(transcendentals), `scripts/S11bB_interface_assembly_sympy_audit.py:29` (sheet-sensitive algebra), and
`mathematica/S11bB_interface_assembly_mathematica_audit.wl:166` (monodromy and non-fixed multiplicity).

⭐ Phase 2 also includes these synthetic adversaries by name: **probe-set vanisher · pole hit · branch
collision · wrong operator with the same kernel · missing duplicate · symbol relabel · typed literal ·
inert premise.** ⛔ **No method adoption and no comparison reduction follows from S9 green alone.**

---

## 2 · What actually costs the days — four negotiation layers

Measured in `checks_S9.yaml` and `engine_output_checks.py`:

| layer | what it is | S9 | S10 | verdict |
|---|---|---|---|---|
| **L1 · tag-name join** | `wl: WL_S9_DETERMINANT` ↔ `py: PY_S9_MAIN_DET_M_FACTORED` | 12 rows | **690 rows** | ⭐ **KEEP** — declaring 12 pairs is free; declaring 690 is the disease. Fix by **granularity**, not by automation |
| **L2 · symbol naming** | `registry_snake_case_to_lower_camel` + 6 hand-written exceptions; `_map_tree`, `_apply_symbol_identities` run **before** comparison; `_symbol_bijection:1157` runs **after** equality fails | 6 naming exceptions + 1 symbol identity; load-bearing subset measured | 119 lines | ⛔ **DELETE** — this is the layer value-comparison removes outright |
| **L3 · shape negotiation** | `wl_select: scalar` / `py_select: last_pair_second` / `list_of_pairs_second` / `sequence_third` | 8 selectors | pervasive | ⛔ **REPLACE** with a declared **object kind** — see §3 L0 |
| **semantic binding / domain** | coordinate/ambient order · coefficient field · parameter premises · exceptional rank strata · branch conditions and sheet metadata · pole avoidance and degree bounds · root multiplicity · defining polynomial | already required by L6 modes | cross-step and new-kind dependent | ⭐ **MEASURE SEPARATELY** — these are not shape properties and do not vanish with selectors |

⚠ **L3 is already costing coverage.** `checks_S9.yaml`'s own comment: *"The main action's residual is
per-root in the `.wl` and a per-root list in the `.py`; see the step record for why that pair is compared
by hand rather than mapped."* ⇒ a pair **dropped from automated comparison because the shapes would not
line up.** ⛔ Nobody has named this layer before; it is not a naming problem and value-comparison alone
does not fix it.

⭐ **These are not shape properties.** L0's object `kind` does not remove them: **L0 must also
carry semantic symbol roles, domains, coordinate bases, branch data, and multiplicity semantics.** A
squared-variable identity of the form `x2 = x**2` is an **algebraic identity, not a spelling exception**. Deleting the
symbol-identity machinery does not remove it; it relocates into the symbol→value map or the object schema.
⇒ L0 may trade selector negotiation for **kind-taxonomy negotiation**. The fourth layer's cost must be
**measured separately**, ⛔ not assumed to vanish.

---

## 3 · The check stack — ⭐⭐ ranked by what it catches, cheapest first

⛔⛔ **The reordering IS the plan.** Cross-engine comparison — where nearly all the days went — sits
**below four cheaper checks** and is **blind** to the shared-spec defect class that S9 already
measured (`steps/S9_light_requires_shear.md:188`: the wrong homogeneity test came from the shared
directive, so *both engines computed it the same wrong way and agreed*).

### ⛔ Declared standing hole · shared-spec blindness

Both review legs walked L0–L7 against a defect introduced identically into both engines from the shared
specification, and **no layer fires**. The exhibited class is a **dimensionless coefficient error in the
shared construction**: it preserves declared kinds, shapes, cardinality, multiplicity, and structural
scaling behaviour; the dimension and unit-covariance checks cannot distinguish it; every emitted root
remains a genuine root of the wrong constructed object, leaving L5 self-consistent; and cross-engine
comparison sees identical objects. This is a declared standing hole, ⛔ not a possible later reporting
outcome.

⭐ The harness requires **at least one physics-bearing oracle whose construction does not come from the
shared specification**: an independent re-derivation from the action, a separately commissioned
metamorphic invariant, or an external published result consulted **only after** computation. ⛔ **An
uncaught priority shared-spec mutant blocks any coverage reduction.**

### L0 · Emission contract — every emitted object declares a KIND

⭐ **This is the keystone; everything else is cheap once it exists.** Each engine emits, per declared
object, a record carrying a declared kind and a canonical serialisation for that kind:

| kind | serialisation | S9 examples |
|---|---|---|
| `scalar` | one expression | `factored_determinant`, `transverse_speed_squared` |
| `multiset` | unordered, multiplicities preserved | `full_root_set` |
| `intvec` | integer exponent tuple | the **6** dimension rows |
| `matrix` | declared shape, row-major | dynamical matrix |
| `subspace` | ⭐ a **subspace**, not a basis | nullspaces |
| `branchy` | ⛔ a **quarantined kind**: carries `Piecewise` / `sign` / `ConditionalExpression` / `Exp[I·]` | S9's anisotropic third-root sign |
| `boolean` | a predicate value | control flags |

⇒ **shape becomes a property of the declared object, not a per-engine accident** ⇒ **L3 dies**, and the
dropped scaling-residual pair becomes comparable.

⭐ In addition to `kind`, every L0 record carries the semantic symbol roles, domains, coordinate bases,
branch data, and multiplicity semantics required by its declared comparison mode.

### L1 · Structural gate — free, no CAS

Kind matches · shape matches · cardinality matches · both engines emitted it. Fails loudly and early.
⚠ This is where `no_comparison` currently hides: S9's dimension gate reached **312 of 1,559** wl tags.
⭐ A structural gate makes non-comparison **a reported number**, ⛔ not silence.

### L2 · Integer-vector dimension equality — free, no CAS

⭐ **Dimensions are an abelian group** — this is why Osprey (Jiang & Su, ICSE 2006) checks units with type
inference plus **Gaussian elimination** rather than a type checker. ⇒ dimensional consistency is **linear
algebra over integer exponent vectors**, ⛔ not a `sp.simplify` problem.
⇒ **6 of S9's 12 rows** become exact integer-tuple equality plus a rank/consistency solve. Instant.

⚠⚠ **But count them as ONE check applied six times, ⛔ not six independent checks.** Dimensional analysis
is a **homomorphism** check. It catches unlike-dimension sums and dropped dimensioned factors — *exactly*
the S9 defect the second engine caught. It is **structurally blind to every dimension-preserving error**:
a wrong dimensionless coefficient, a sign flip, a transposed index, a compensated power.

### L3 · unit-covariance check — ⛔ DEMOTED 2026-08-07, and the demotion is measured

For each emitted expression object, substitute `t → α·t`, `x → β·x`, and each dimensioned parameter by its
declared scaling; require the object to scale as `α^a β^b γ^c` with the exponents it declared.

⛔⛔ **RETRACTED — this section used to call L3 "the fabrication-resistant one" and claim "a typed literal
fails it."** ⭐ A review leg exhibited counterexamples with a script and literal stdout, and they are
correct **by inspection** — a dimensionless factor changes nothing about how an expression scales:

```
DECLARED A -> gamma*beta**-1*alpha**-2*A; B -> gamma*beta**-3*B
=== WRONG dimensionless coefficient 17*A/B ===  L3 VERDICT: PASS
=== WRONG SIGN -A/B ===                       L3 VERDICT: PASS
=== FABRICATED dimensionless Q=7 ===           L3 VERDICT: PASS
```

⇒ ⭐ **What L3 actually catches, stated at the width the measurement supports — TWO things:**
1. a **bare number** claimed to carry a non-zero dimension;
2. an **inhomogeneous sum** (`A/B + 1`).

⇒ ⛔ **What it does NOT catch:** any **dimensionally correct** fabrication · any **dimensionless**
fabrication (multiplicities, flags, counts, pure-number residuals) · wrong coefficients · wrong signs.

⚠ **And it is partly circular:** it compares the observed scaling of a body against **declared** exponents
under **declared** parameter scalings. A co-fabricated body-plus-declaration passes.

⛔⛔ **L3 IS NOT A FABRICATION DEFENCE. Do not let it back into R4 as one.** The plan's own §6 already
names **mutation** as the fabrication instrument, and elevating L3 contradicted it. ⭐ L3 is a cheap
**monomial scaling consistency** check — worth keeping at that width, ⛔ worth nothing beyond it.

⚠ Further counterexample classes bound this check. **`0` satisfies the scaling law for every declared
exponent.** A matrix or coefficient list may legitimately carry different dimensions per entry, so one
scalar exponent is the wrong covariance law. A `Piecewise`/conditional object requires branchwise
dimensions plus dimensionally valid conditions, ⛔ not one global test.

### L4 · Metamorphic relations — S9 already has one; name the family

⭐ The scientific-computing MT literature elicits MRs from **monotonicity, conservation, scaling** and
domain expectations. **Scaling is a named MR family**, and S9's `k → λk` test
(`ω²(λk) − λ²ω²(k)`) is one. ⇒ ⭐ **the existing 9-action control matrix IS a metamorphic relation set**
and stays exactly as it is.

⭐ **Also name what the controls already are: manufactured solutions.** MMS (Roache; term coined by
Oberkampf & Blottner) is *the* code-verification standard in computational physics — an answer known **by
construction**, pushed through the whole machinery.

⛔⛔ **LEAK REPAIRED 2026-08-07 — this paragraph used to STATE S9's answer** (a speed formula and a mode
count) while the next sentence forbade exactly that. ⭐ A review leg caught the self-contradiction. ⛔ It
matters because step 1 turns this plan into `HARNESS_STANDARD.md`, **the artifact both engines read**, and
rule 5 says *a builder iterating to exit 0 converges on any target it can see.*

⇒ ⭐ **Every control states its expectation STRUCTURALLY and only structurally.** The admissible form:

> *"the gradient-elastic control must emit a **longitudinal root in addition to** the transverse family;
> the curl-only action must not"* · *"the anisotropic control must **reduce** the transverse multiplicity"*
> · *"the sign control must **change the sign** of at least one emitted `ω²`"*

⛔ **Inadmissible in this document or anything downstream of it:** a speed, a ratio, a formula, a specific
multiplicity, or any number the builder could iterate toward. ⭐ **A structural expectation is still
able-to-fail** — it names *what must move*, ⛔ never *what it must equal*.

### L5 · Witness / result checking — ⭐ single-engine, catches PHYSICS, costs hours

Blum–Kannan result checking, and the Isabelle/HOL CAS-certification pattern (arXiv 2102.02679): ⛔ **do
not verify the derivation — substitute the answer back.** *"The checker is much simpler than the
algorithm, yet it is all the user has to trust."*

| the engine EMITS (an object) | the HARNESS DECIDES (a verdict) |
|---|---|
| the root `ω²₀` **and** `det M` evaluated at it | is that residual zero? |
| the nullvector `v` **and** `M·v` **and** `v` itself | is `M·v` zero and `v` non-zero? |
| the factored determinant **and** its expansion **and** the original | is the difference zero? |
| the corank of `[M; kᵀ]` recomputed at the emitted root | does it match the multiplicity emitted elsewhere? |

⛔⛔ **THE SPLIT IN THAT TABLE IS THE POINT, AND IT IS A CORRECTION TO THIS PLAN'S FIRST DRAFT.** The draft
wrote witnesses as `det M(ω²₀) ≡ 0` — i.e. **as an assertion inside the engine.** ⛔ That is the defect
class this whole rebuild exists to remove, reintroduced through the back door: an engine that asserts
prints `0` whether or not anything computed it, and ⛔ **a residual asserted zero carries no information**
(`CLAUDE.md` rule 2). ⚠ Worse, the review-legs standard already records that **an `assert` before an emit
blinds the ablating leg** — a perturbation strong enough to flip the check kills the process, so the leg
sees only PASS-or-crash.

⇒ ⭐⭐ **Engines emit both operands and the residual. The harness decides.** No engine emits a verdict, a
`PASS`, a boolean conclusion, or a value it asserted rather than computed. ⭐ Same rule applies to **L3**:
the engine emits the **rescaled object**; the harness compares it to the declared scaling.

⚠ **And the tautology bar still applies:** a witness must not be evaluated with the predicate that
produced the value (`c ≔ √x` then reporting `c² − x`). ⇒ ⭐ **each witness names which INDEPENDENT object
it substitutes into**, and that name is checked at review, ⛔ not assumed.

⭐ **A wrong root fails this with no second engine present.**

### L6 · Cross-engine value comparison — ⭐ the expensive one, now scoped

⛔⛔ **INVERTED 2026-08-07 — the first draft had this backwards, and a review leg PROVED it with a mutant.**
The draft replaced `sp.simplify(lhs - rhs) == 0` with evaluation at a **recorded** point set. ⛔ Recording
the points makes them **targetable**, and a wrong expression can be built to vanish on exactly them:

```
mutant = true * (1 + ∏(kx/ky − i), i = 1..8)      # ⇒ agrees at all 8 recorded probes
 rescaling_residual = 0        ← passes L3 (same dimension AND same scaling)
 PROBE_VECTORS_EQUAL True      ← passes the proposed L6
 CURRENT_SYMBOLIC_EQUAL False  ← ⭐ TODAY'S COMPARATOR CATCHES IT
```

⇒ ⭐⭐ **As drafted, L6 was a REGRESSION from the comparator we already have.**

### ⇒ The corrected rule

1. ⭐ **Symbolic residual equality stays PRIMARY** for algebraic scalars — it is already a strong oracle and
   ⛔ nothing may be adopted that is weaker than it on any row.
2. ⭐ **Exact-point evaluation is a SECONDARY, DIFFERENTIAL technique** — it decides rows symbolic equality
   cannot close, and it cross-checks rows it can. ⛔ It is not a replacement.
3. ⛔⛔ **Probe points are NOT fixed and recorded in the config.** Generate them **after** the comparator is
   implemented, and draw **fresh exact points on every adjudication run**. ⭐ A seed recorded for
   reproducibility must be re-drawn for adjudication.
4. ⭐ **Degree/pole/domain bounds are derived and enforced INDEPENDENTLY** — ⛔ not chosen by the party
   whose comparator they gate, and ⛔ never widened to fit a case.

For each **declared pair**, evaluate both sides at exact points — rationals or a finite field, ⛔ **never
floating point** — under those four conditions.

⛔⛔ **Fingerprints are an EQUALITY TEST on intentionally paired objects. They are ⛔ NOT auto-discovery of
which tag matches which.** (Review-leg correction, accepted: content-addressed joining can collide, and
two physically different objects can coincide.) L1's declared pairing stays. What dies is L2 naming.

⛔⛔ **THE COMPARATOR MODE IS CHOSEN BY OBJECT KIND, AND EACH MODE HAS CONDITIONS.** Schwartz–Zippel is
stated for a **low-degree multivariate polynomial over a field** — ⛔ not for arbitrary CAS expressions.
Our own artifacts already leave that class (`sp.Integral`, `sech²`, `tanh`, `Abs`, natural-language
strings, deliberately-unsimplified square roots, sheet filtering and monodromy in S11bA/S11bB).

| kind | mode | required conditions ⛔ none optional |
|---|---|---|
| `scalar` polynomial/rational | exact PIT at the point set | denominator clearing · **pole avoidance** · explicit **degree bound** · points/primes reported after the run · fresh adjudication draw |
| `matrix` | entrywise PIT | + declared **coordinate order** · structural (shape, rank) checks |
| `multiset` (roots) | multiset **of values** per point ⇒ ⭐ ordering stops being a disagreement | + defining polynomial · residual witness · **multiplicity** · sheet metadata |
| `subspace` | RREF **of the transposed** basis matrix (⛔ naïve RREF of a column-basis matrix does not canonicalise its column space), **or** mutual cross-annihilation + equal nullity | + fixed ambient coordinate order · declared coefficient field and parameter premises · **exceptional strata where rank changes** · each basis must annihilate **its own** operator · ⭐ **keep comparing the operators too — different operators can share a kernel** |
| `intvec` | exact integer-tuple equality | decided at L2; ⛔ never routed through a CAS |
| `branchy` | ⛔⛔ **never silently probed** — per-branch with the branch condition as part of the object, else symbolic | ⭐ **row marked symbolically-decided in the report** |
| transcendental / `Integral` | symbolic identity, metamorphic relation, or **certified interval/ball evaluation over a declared domain** | ⛔ PIT does not apply |
| textual / natural-language | ⛔ **not comparable** — flag as an emission defect | ⚠ these exist today in S11bA and are a finding in their own right |

⭐ **Compare evaluation vectors directly; hash them only for storage and indexing.**

⭐ **Every row reports WHICH LAYER DECIDED IT.** The count decided by symbolic fallback is a published
number, ⛔ not silence.

### ⛔ Named finding · textual emissions are engine defects

⭐ **An emission that is not a parseable CAS object is an ENGINE DEFECT.** It is comparable by no mode,
and it lands on the S11b rebuilds as work — ⛔ not on the comparator as a gap to accommodate.

### L7 · Mutation adequacy — ⭐ the meta-check, and it is not optional

⛔⛔ **A layer with no mutation that flips it is not a check.** Measured basis: of S11 engine 1's nine
defects, **none was visible by reading and the ablating leg found 8 of 9**; and there is now published work
on **mutation-based adequacy metrics for metamorphic relations** in scientific computing (arXiv 2605.17437),
i.e. this is a recognised way to ask *"is my control set enough?"*.

⇒ Extend `harness_ablation.py` so **every layer L0–L6** has at least one mutation that flips it, with
**literal stdout recorded**. Output-representable defects may still use its engine-free in-memory
mappings; the fabrication and inert-premise classes obey §3b's source-level requirement.

### L8 · ⭐⭐ CROSS-STEP CONSISTENCY — the existing oracle OUTSIDE a single step's spec

⚠ **Missing from the first draft entirely.** The whole stack above is **within-step**: does engine A agree
with engine B at step N. ⛔ It never asks whether step N's `c_γ` is step N+1's `c_γ`.

⭐ **The instrument already exists and is the right shape.** `registry_residual` binds a relation from
`reduction/relations.yaml` to tags **this step's engine emitted**. Both S9 and S10 declare **R4**:

```
S9:  rho_br, mu_R, c_gamma  ←  ALL THREE from one tag, WL_S9_CANDIDATE_SPEED_SQUARED1
S10: rho_br ← Q6_INERTIAL_COEFFICIENTS · mu_R ← Q6_STIFFNESS_COEFFICIENTS
     c_gamma ← Q3_DISTINCT_ROOTS
```

⇒ ⭐⭐ **Why this matters more than its size suggests: it is currently the ONE implemented check whose
oracle does not come from the step's own shared specification.** Every layer L0–L7 is downstream of the
directive both engines read, which is why L0–L7 are blind end-to-end to a shared-spec defect (both legs,
independently, with a passing mutant). ⭐ Step N+1 has its **own** directive ⇒ a defect in step N's spec can
surface as a **cross-step residual** that neither of step N's engines could see. ⛔ It does not close that
hole; §3's separately constructed physics-bearing oracle is still required.

### ⛔ Two defects in the instrument as it stands, both live

1. ⛔⛔ **S9's R4 residual is NONZERO** and has been recorded rather than resolved — S9's record:
   *"REGISTRY_RESIDUAL: nonzero=1 (R4 — it compares a squared speed against a speed)."* ⇒ the one
   cross-step check we own is reporting a **category mismatch** and is being carried as a note.
2. ⛔ **S9's binding is near-self-referential** — all three quantities are read out of **one** emitted
   expression, so within S9 it largely checks that expression against itself. ⚠ That is the
   *control-inside-the-thing-it-polices* failure. ⭐ S10's three-tag binding is the strong form.

### ⇒ What L8 requires

- ⭐ Every step **declares which registry quantities it PRODUCES and which it CONSUMES.**
- ⛔ A binding must read from **independent tags**, ⛔ never several qids out of one expression.
- ⭐ A nonzero cross-step residual is a **BLOCKING finding**, ⛔ not a recorded note.
- ⚠ **Coverage is the honest limit:** `relations.yaml` holds **5** relations (R1, R2.h0, R2.xi_h, R4, R5)
  over **14** quantities in **2** sectors — against ~25 remaining steps. ⇒ ⭐ **growing the registry is
  part of every step's cost**, ⛔ not a separate project, and a step that produces a quantity no later step
  consumes gets **no cross-step check at all** — say so in its record.

---

## 3b · ⛔⛔ WHAT MAY BE CUT IS DECIDED BY MEASUREMENT — ⛔ NOT BY THIS PLAN

⚠⚠ **This section exists because the first draft got it wrong, and the correction is load-bearing.**

The draft argued for cutting from ~3,750 emissions to a handful of headline objects, on the grounds that
*"no community compares two CAS at the level of intermediate tagged emissions."* ⛔⛔ **That is FALSE.**
Both review legs caught it and the primary source confirms it:

> Knight & Leveson's versions produced **241 results**; *"The launch condition is the only true output…
> **The other results are really intermediate**"*; failure was declared on *"any discrepancy between the
> 241 results"*; and the stated reason is *"The intermediate results were compared because, **in a
> practical multi-version system, votes would be taken on intermediate results**."*

⇒ ⭐ **The canonical N-version experiment compared intermediates on purpose.** Granularity is ⛔ not a
reason to cut.

⭐⭐ **The REAL difference, and it is the one to act on: they had a GOLD PROGRAM.** Each of their 241
intermediates was compared against a **reference**, so each cost nothing. ⛔ We have no reference, so every
intermediate is a **peer negotiation between two engines**. ⇒ **the cost is in the negotiation, not the
count** — which is exactly what §2's four layers measure, and it is why L0–L3 (kinds, structure, integer
vectors, unit covariance) plus semantic binding/domain are the right target: ⭐ **they make an
intermediate cheap to compare rather than fewer intermediates to compare.**

### ⇒ The reduction gate

⛔ **No emission is dropped from comparison until a MUTATION COVERAGE STUDY says which contract catches
which defect class.** The reduction gate has all of these requirements:

1. Seed **multiple mutants per class across different object kinds**, including the A6 catalogue, the
   shared-spec class declared in §3, and S11b-like object kinds. Include **coupled/common-mode mutants**.
2. ⛔ The **fabrication** and **inert-premise** classes are tested by mutating **engine source** or an
   instrumented derivation IR, ⛔ **never output mappings**. A committed output cannot distinguish a value
   computed from the model, the identical payload typed literally, a premise used in construction, and
   the same premise merely printed afterward. The existing honest instrument is
   `reduction/derived_or_declared.py:318`: it labels a tag derived when a selected perturbation changes
   its text. Its own recorded caveat at `:332-337` says constant status is not proof of literal source.
   The converse also fails: a typed **symbolic** literal can change under perturbation and be
   misclassified as derived.
3. Hold back mutants that are **not visible to the builder**, and run them only after the candidate
   comparator is frozen.
4. Run **deliberately weaker comparator implementations as negative controls**.
5. The stated **kill threshold** is that every priority mutant is killed; the orchestrator-only
   adjudication artifact pins any class-specific thresholds before the run. Then record which contracts
   fire for the visible and held-out catalogues.
6. ⛔ **Any uncaught priority mutant blocks reduction.** In particular, an uncaught priority shared-spec
   mutant blocks any coverage reduction; it is not merely published.
7. The catalogue must cover the shared-spec class and S11b-like object kinds, or carry the explicit
   result **"not yet tested ⇒ no reduction for those steps."**
8. Only after those gates may a minimal hitting set be considered, and it still retains every provenance
   and inventory check regardless of whether it fired.
9. ⭐ Publish additional classes that no contract catches as a **report**, ⛔ never as a substitute for
   the priority-mutant gate.

⚠ **This study has not been made.** ⛔ Until it is, the honest position is: *we do not yet know which of
the 3,750 can go*, and any number in this plan describing a smaller comparison surface is a **hypothesis**,
⛔ not a target.

---

## 3c · ⛔⛔ THE BUILDER MUST NOT SEE THE ORACLE — leak inventory and separation

⚠⚠ **Both legs independently ruled the plan leaks, and the inventory is larger than the one sentence I
asked about.** Predeclared outcomes reachable by whoever builds the comparator:

| leaked | why it is a target |
|---|---|
| the former speed expression and mode count (L4) | ⛔ a **physics answer**, in the paragraph that forbids exactly that |
| the former expected multiplicity in the witness table (L5) | an expected outcome |
| S9's baseline aggregate (§1, old A1) | a **predeclared verdict**, ⛔ not merely an inventory count |
| the former named load-bearing declarations and concrete identity (§1/§2) | names the answer to A2 |
| the root-set cardinality in the L0 table | a physics cardinality |
| the former computed-value outcome (§6) | tells the builder what a green run looks like |

⭐ **Rule 5's whole content: a builder iterating to exit 0 converges on any target it can see.** ⛔ Calling
this a "method" document does not undo the information flow — step 1 turns it into the artifact both
engines and the comparator read.

⇒ ⭐⭐ **STRUCTURAL FIX — separate the two artifacts, and freeze in this order:**

1. **Builder-visible:** the emission contract, the object kinds, the layer definitions, the *generic*
   per-kind comparator logic. ⛔ No S9 values, no verdicts, no cardinalities, no named declarations.
2. ⛔ **Orchestrator-only, unreachable from the repository sandbox, and applied ONLY AFTER the generic
   comparator is frozen:** the adjudication artifact exists at
   `steps/S9_PILOT_ADJUDICATION.md` (in-tree, carrying a warning header — ⛔ NOT a custody control).
   It holds expected per-row verdicts, the load-bearing declaration list, held-out-mutant and probe-point
   stubs, removed result values, and structural control expectations. ⭐ **Held-out mutants and fresh probe
   points are the part that makes freezing meaningful.**
3. ⭐ **Explicitly prohibited and tested for: any S9 tag-name or S9-value special case in the comparator.**

⚠ **This is not theoretical.** `REBUILD_HANDOFF.md:218`: *"Three times a half-fix passed the new test, the
whole suite, the full battery and produced byte-identical comparator output."*

---

## 4 · ⭐⭐ Acceptance gates and reporting requirements — no new physics

⛔ **No expected physics value appears anywhere in this list.** Every criterion is an **invariant** — a
count, a coverage, a which-layer-decided — ⛔ never a sampled outcome.

| # | criterion | why it is the right test |
|---|---|---|
| **A1** | ⛔⛔ **REPLACED 2026-08-07 — the old A1 was DEFEATED BY EXECUTION.** A leg deleted the naming sections, replaced both equality functions with `return True`, ran the real S9 machinery, and got `FAKE_COMPARATOR_COUNTS {'AGREE': 12}` · `COVERAGE 12 12 GAPS ()` · `OPERANDS_ACTUALLY_TESTED_FOR_EQUALITY 0`. ⇒ the advertised outcome was reachable **with zero equality checks**. ⭐ **New A1: PER-ROW verdict equivalence against the CURRENT symbolic comparator, over the baseline AND an adversarial mutant set — ⛔ not an aggregate count** — with the row inventory pinned by an artifact **independent of the candidate comparator**. **Runnable criterion:** run both comparators on that pinned baseline and mutant inventory and compare each row's verdict. | ⛔ an aggregate count is not a bar; a per-row differential against a known-strong oracle is |
| **A2** | Run `declaration_load_ablation.py` against the **declaration set as it stands today**, pinned independently before the candidate is built. Deleting a naming section must not shrink the declaration set the ablator iterates; any declaration that remains load-bearing is an A2 failure. | prevents deletion of the measured input from satisfying its own ablation |
| **A3** | Every row names its deciding layer; an **independent confirmation verifies that the claimed layer actually decided the row**; the symbolic-fallback count is printed | ⛔ prevents an unverified label from hiding the actual deciding predicate |
| **A4** | The pair currently dropped for shape mismatch (main-action per-root scaling residual) is compared **by its kind's declared mode** | excludes an always-true or projection-only predicate from satisfying “is compared” |
| **A5** | ⭐ Every layer L0–L6 has ≥1 mutation that flips it, with literal stdout at a named path | ⛔ without this the stack is a demo |
| **A6** | Seeded checks — **must FIRE**: a sign flip in one engine's determinant · a coefficient change · a **form** change · a dropped dimensioned factor · wrong root multiplicity · a typed bare literal carrying a non-zero claimed dimension · one-sided fabrication against a computed counterpart · a **wrong operator that shares a kernel**, firing the operator comparison · a pole/branch error · a compensated error · a missing duplicate. **Must NOT fire**: a permuted root order · a pure alpha-renaming **only after semantic symbol roles exist** · a differently-written nullspace basis. Until semantic roles exist, a relabel legitimately fails loudly under R3. | ⭐ tests both directions; ⛔ a check that only fires is not calibrated |
| **A7** | `checks_S9.yaml` under **~120 lines** with the naming section gone; the S10 config is **mechanically re-derived** and measured against its 3,121; and the **compared row set is preserved** | excludes minification or row-discarding from satisfying the size measurement |
| **A8 · REPORTING, NOT A GATE** | ⭐ Total wall-clock for a full S9 re-verify from committed outputs is **recorded** | disclosure cannot pass an adoption decision |
| **A9** | ⭐⭐ The **mutation coverage study** of §3b passes its stated threshold, negative controls, held-out catalogue, priority-mutant bar, and untested-step rule | ⛔ the **gate on any reduction**; without it a smaller comparison surface is a guess |
| **A10 · REPORTING, NOT A GATE** | ⭐⭐ Additional defect classes that no contract catches are published | disclosure cannot pass an adoption decision, and the standing shared-spec hole remains governed by A9 |

⛔⛔ **THE STOP CONDITION — and it is the point of a pilot.** If **A1 fails per-row verdict equivalence**,
the proposal is **dead**, we learned it on the smallest step, and we go back to the current method with
one measurement more than we had. ⭐ Say this out loud before building, so a failing A1 is a **result**,
⛔ not a thing to iterate around.

---

## 5 · Sequencing — who builds, who reviews

⛔ **Rule 7 governs: whatever writes does not review.** ⛔ No commit before both legs report.

| # | step | author | legs |
|---|---|---|---|
| 0 | **this plan** | **Codex applying an orchestrator decision list** | **fresh Claude agent + Grok** ⛔ before anything is built |
| 1 | **Emission contract + check-stack spec** (`HARNESS_STANDARD.md`) — the artifact **both engines and the comparator read** ⇒ physics-bearing by rule 7 | **assigned non-orchestrator author** | **Codex + Grok** |
| 2 | New comparator, built against S9's **committed** outputs — ⛔ engines untouched | **a different Codex session from the session that edited this plan** | fresh Claude + Grok |
| 3 | Run acceptance gates **A1–A7 and A9**; record A8 and A10 separately; preserve literal stdout | orchestrator | — (measurement, not a claim) |
| 4 | ⚠ **Only if** the acceptance gates need engine changes (new kinds, rescaling hooks): **one engine at a time** | Codex | fresh Claude + Grok, **per engine** |
| 5 | **Phase 2:** add the required S10/S11b representative rows and named adversaries; mechanically re-derive the S10 config; rerun the gates before adoption or reduction | Codex | fresh Claude + Grok |
| 6 | Fold the standard into `REBUILD_HANDOFF.md`; ⛔ retire "S10 is the template" | orchestrator | Codex + Grok |

⭐ **Step 2 is deliberately engine-free.** Rebuilding the comparator against frozen outputs means the
engines cannot be tuned to make the comparator green — which is the co-authorship failure the project has
already recorded.

⭐ **Authorship change:** this revision was applied by Codex from an orchestrator decision list, so rule 7
requires a **fresh Claude agent + Grok** as the two review legs for this document, ⛔ not Codex + Grok.
⛔ **Step 1's `HARNESS_STANDARD.md` is not authored by the orchestrator:** three consecutive
orchestrator-written artifacts each carried a central defect (rule 15). ⛔ **The Codex session that edits
this plan must not be the Codex session that later builds the comparator.**

---

## 6 · ⛔ What does NOT change

- ⛔ **The physics.** If any computed value differs from the committed output, **STOP and report it — ⛔ do
  not adjust either engine or the comparator.**
- ⛔ **The 9-action control matrix** — it is the metamorphic-relation set and it stays.
- ⛔ **`derived_or_declared.py` and the `.premises` sidecar** — fabrication triage, ⛔ not source proof.
- ⛔ **The mutation legs.** ⭐⭐ **The instrument that catches fabrication is source/derivation-IR
  mutation, ⛔ not comparison or output mutation.**
- ⛔ **Two engines.** ⭐ S9 measured the second engine catching a wrong dimension that two review legs and a
  full ablation suite missed. ⛔ This plan scopes dual-engine; it does **not** drop it.

---

## 7 · Risks, and what each is mitigated by

| # | risk | mitigation | residual |
|---|---|---|---|
| **R1** | ⛔ **Branch-carrying objects.** 17 hits of `ConditionalExpression`/`ComplexExpand`/`Sign[`/`Reduce[`/`Exp[I` in S11's `.wl` alone; S9 already has one unparsed `Piecewise` root sign | the `branchy` kind: ⛔ never silently probed; per-branch or symbolic-with-a-label | ⚠ a step dominated by branchy objects gains little. **Record that per step**, ⛔ do not average it away |
| **R2** | Accidental agreement / fingerprint collision | ⭐ pairing stays **declared**; ≥8 independent exact points; point count reported | one-sided by nature — probing proves *different*, never *same* |
| **R3** | ⛔ The symbol→value map becomes the new negotiation | derive from the registry where possible; both engines must accept the **same** map or the step fails loudly; pure alpha-renaming may be transparent **only after semantic symbol roles exist** | ⚠ until those roles exist, a relabel legitimately fails loudly; this is the honest cost, ⛔ not zero |
| **R4** | ⛔⛔ **Shrinking the comparator shrinks fabrication coverage** | ⭐ L5 plus source/derivation-IR mutation at L7 are the fabrication defences and they **GROW**. ⛔ A typed literal fingerprint-matches another typed literal — value comparison does **not** fix fabrication | ⛔ **If witness, provenance, inventory, or source/IR mutation coverage is descoped, this plan is not safe to run** |
| **R5** | "It works on S9 because S9 is small" | **Phase 2** exercises the required S10/S11b kinds and adversaries before adoption or reduction | ⚠ S9 green alone authorises neither |
| **R6** | The plan itself is a shared spec ⇒ rule 7's failure class | step 0 and step 1 both get two legs; the spec is the highest-value review target | ⚠ Knight & Leveson says this cannot be fully closed ⇒ ⭐ **and their answer is to USE the correlated failure as a spec-repair signal**, which is what our three S11 repair rounds actually were |
| **R7** | ⛔⛔ **The first draft's reduction rationale was a FALSE UNIVERSAL** ("nobody compares intermediates"), caught by both legs | §3b replaces it with a **measurement gate** (A9); A10 is reporting only | ⚠ ⭐ **The measurement has not been made.** ⛔ If A9 is descoped, this plan reverts to the same unjustified cut, and the fact that it *feels* obviously right is precisely what made it survive one draft and two of my own readings |

---

## 8 · What this is worth

Measured: **220 commits over 7 days → 2.5 closed steps.** Remaining: S1–S8 re-banking · S11 finish ·
S11b-A/B rebuild · S11b-C · S12–S23 ≈ **25 units**. At the current rate that is **8–10 weeks**, nearly all
of it comparator, spec and naming.

The brief's own bar: *a method that takes a week to stand up must save more than a week.*
⭐ Phase 1 is bounded — frozen outputs, an existing answer key, a falsifiable A1 — and either advances to
the required Phase 2 or fails **cheap and visibly**. ⛔ It does not clear the adoption bar by itself.

⚠ ⛔ **What it does NOT buy, stated so nobody reads it in:** it does not make the physics faster to *find*,
it does not close S9's open debts (P2 uncomputed; the longitudinal absence **assumed, not derived**; the
curl-only conditionality), and it does not remove the shared spec as a correlation channel — ⭐ nothing
does. It buys back the **days per step** currently spent making two programs agree on how to *say* things.
