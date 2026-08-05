# Repair `S10_SHARED_PHYSICS.md` — ⭐ SEVEN targeted edits, ⛔ NOT a rewrite

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md`
⛔ **Edit in place. ⛔ Do not commit. ⛔ Do not modify any other file.**
⚠ The as-built text is pinned at tag `s10-as-built`; ⭐ the step record cites that, ⛔ not your output.

⭐⭐ **SCOPE, and it is deliberately narrow (user, 2026-08-05): keep S10 to S10.**
⛔ **A full rewrite into a schema-checked contract with a tag registry is DEFERRED** —
⇒ `S10_spec_rewrite_directive_v2.md` + `S10_contract_DEFERRED_findings.md`. ⭐ That work buys **future
cross-engine comparability**, ⛔ not S10's correctness. ⛔ **Do not do any of it here.**

⛔⛔ **EXPLICITLY OUT OF SCOPE — ⛔ do not touch, even if you see it is wrong:**
the `§8` tag-name grammar and any quantity registry · the payload serialiser · the sign-diagnostic
classifier · the `Q6d` unknown-coefficient count · `report_contract`-style machinery.
⚠ **All are known-defective and all are recorded elsewhere.** ⭐ If you spot another, ⭐ **report it, ⛔ do
not fix it.**

---

## ⛔⛔ E1 — `§Q7`: rewrite the section. ⭐ It is the reason this directive exists.

⚠ **The section currently instructs TWO DIFFERENT OBJECTS.** One sentence names a fixed curl density; a
later amendment says to use the package's **own** stiffness density. ⇒ ⛔ **the two engines each followed a
different sentence and both exited clean.** ⭐ **Delete, ⛔ do not annotate.**

**⭐ What Q7 IS — state this in the section, ⛔ and note what it does NOT say:**
> The stiffness is written in general `D`, because the 3-vector cross product exists only at `D = 3`.
> `§Q7` **COMPARES this package's stiffness density against the ordinary curl-squared at `D = 3`**, and
> emits both operands and their difference. ⭐ It is a form-and-normalisation comparison bearing on the
> coefficient that carries the wave speed. ⛔ It is **not** evidence for the mode count, which is computed
> in every `D` and does not depend on it.

⛔⛔ **DO NOT WRITE THAT THE GENERAL-`D` FORM "BECOMES" THE CURL-SQUARED. ⚠ An earlier draft of this
directive did, and BOTH legs caught it — twice over.** ⭐ It states Q7's **outcome**, in a section whose
entire content is *whether* that holds; ⇒ ⛔ **a leak.** ⚠ **And it names a SECOND compared object** — the
fixed general-`D` form — alongside the package's own density. ⛔ Those coincide for the baseline package
and **differ** for the form controls, ⇒ **it is the very divergence this edit exists to remove.**
⭐ **Name exactly ONE compared pair, and state no outcome.**

### ⭐⭐⭐ NAME THE OBJECT. ⛔ DO NOT SPECIFY HOW TO OBTAIN IT.

⛔⛔ **THIS IS THE CORRECTION THAT MATTERS, and every earlier draft of this directive got it wrong.**
⚠ Successive drafts specified a **recipe** — zero the velocities, divide by a weight, prove the weight
nonzero, guard the division. ⇒ **Five rounds of review then argued about the RECIPE**: is the weight
unique, is the quotient well-defined, is the residual a tautology. ⛔ **None of that is a question about
the physics.** ⭐ **Every one of those questions was manufactured by specifying a derivation path the spec
never needed to specify.**

⭐⭐ **The object already exists in both engines.** Each action constructor **already returns the stiffness
density it used**, as an object, alongside the Lagrangian. ⇒ ⛔ **nothing needs to be divided by anything**,
and ⭐ **the weight never appears, so it never needs pinning.**

**⭐ The comparison — ⛔ exactly one pair, ⭐ stated as OBJECTS:**
```
g_ij      := INDEPENDENT symbols standing for ∂_i u_j     ⛔ never a k×a amplitude curl
S_pkg     := the stiffness density THIS package's action USES,
             with the g_ij substituted in                 ⭐ the object the action was built from
c_i       := Σ_{j,k} ε_ijk g_jk                           COMPUTED by the CAS
emit:       S_pkg  ·  c·c  ·  (S_pkg − c·c)               all three
```

⛔⛔ **DO NOT write an extraction procedure. ⛔ No quotient, ⛔ no weight, ⛔ no velocity-zeroing step.**
⭐ **Ask for the object; ⛔ let the engine hand over what it built.**

⛔⛔ **`c·c` MUST BE REACHED BY COMPUTATION from the Levi-Civita definition. ⛔ Do NOT write out the
expanded polynomial.** ⭐ That computed side is the **only** reason this comparison can fail on physics —
the density side is whatever the engine used, so a hand-expanded reference would compare a typed object
against a typed object.

### ⭐⭐ AND LET THE ENGINES ANSWER — ⛔ do not pre-empt them in prose

⭐⭐ **If the two engines emit DIFFERENT densities, that is a FINDING, ⛔ not a spec defect to be prevented
by better wording.** ⚠ It is exactly what happened: one engine emits a fixed curl density for every
package, the other emits the package's own. ⇒ ⭐ **a cross-engine comparison of this tag shows it.**
⛔ **Do not try to make divergence impossible by writing more precise prose about how to derive the
object.** ⭐ Name the object, emit it from both engines, and **compare afterwards** — ⇒ that is what the
two engines are FOR.

⭐ **The one requirement that stays, because it is about DATA DEPENDENCY and not about a recipe:**
> Mutating **the action alone**, with the package selector held FIXED, must move the emitted `S_pkg`.

⚠ **Measured:** an implementation keyed on the **selector**, taking no action object at all, passes a
*"change the form in one package and watch it move"* test. ⇒ ⛔ **re-deriving from the selector is not
re-entering at the action**, and ⭐ it needs no extraction algebra to state.

⛔⛔ **AND SCOPE IT HONESTLY — ⚠ an earlier draft claimed this requirement "catches an engine reporting a
density it did not use". ⛔ IT DOES NOT, and a leg proved it.** ⭐ Mutation sensitivity establishes **some**
dependence on the action, ⛔ **not identity with the density the action used.** ⚠ **The counterexample,
computed:** an engine whose action uses a density while its reporter emits a **rescaled** multiple of it
passes the reference residual **and** moves under mutation — ⇒ ⛔ **both surviving checks pass while the
reported object is the wrong one**, and the rescaling is exactly the normalisation error `§Q7` exists to
catch.
⇒ ⭐ **State in the spec what actually polices this: the SECOND ENGINE.** An engine that uses one object
and reports another is **internally inconsistent**, and ⭐ the cross-engine comparison of this tag is what
surfaces it. ⛔ A within-engine algebraic identity cannot — ⚠ a leg verified the removed one was
**identically zero on every conforming build**, so it fired only on non-conforming ones, and ⛔ it never
caught a **shared** wrong formula either.
⭐ **Require the emitted object to BE the one the action was assembled from** — ⛔ not an equal-valued
rebuild — and ⭐ **say plainly that provenance rests on cross-engine comparison, ⛔ not on this residual.**

⛔⛔ **THE "Measured:" RATIONALE ABOVE IS FOR YOU, ⛔ NOT FOR THE SPEC.** ⭐ It explains why the
requirement is worded as it is; ⛔ **do not copy it into `S10_SHARED_PHYSICS.md`** — ⚠ a measured
failure shape steers an implementer, and `E6` forbids exactly that shape in the target file.

⛔ **DELETE from `§Q7`:** the sentence naming a fixed curl density as the compared object · the sentence
stating what the residual is **expected** to be and for which packages · the sentence describing the
**measured** failure shape under a corrupted action · the "curl vector of the **amplitude**" wording, which
contradicts the derivative formula on the next line.

⚠ ⭐ **Say plainly in the spec what this comparison does NOT do:** it does ⛔ **not** establish that the
emitted density is the one the action was assembled from. ⭐ **Nothing inside a single engine does.**
⇒ the action-mutation requirement establishes **dependence**, and ⭐ **cross-engine comparison of this tag
is what establishes provenance.**

## ⛔⛔ E2 — add the DISTINCTNESS premises. ⭐ Without them a control can police nothing.

⚠⚠ **CORRECTED — an earlier draft of this directive got the facts wrong and both legs caught it.**
⭐ **What is actually true, verified in both engines:**

| condition | in `§7`'s package row? | in the joint premise set? | in the engines? |
|---|---|---|---|
| anisotropy-scale distinctness | ⭐ **yes, already** | ⛔ no | ⭐ **yes, both** |
| coefficient-scale distinctness | ⛔ **no** | ⛔ no | ⛔ **neither** |

⇒ ⭐ **The anisotropy condition is being CENTRALISED** (it exists, in the package row, and both engines
picked it up); ⛔ **the coefficient-scale condition is genuinely MISSING everywhere.**
⭐ **Add both to the joint premise set, as package-domain premises.**

⛔⛔ **STATE THE CONSEQUENCE FOR THE COEFFICIENT-SCALE PACKAGE ONLY. ⚠ Two earlier drafts got this wrong in
OPPOSITE directions** — the first said both controls were "dead"; the second said **each** admitted domain
contains its collapse locus. ⛔ **Both false.**

⭐ **What is true, and it is asymmetric:**
- **Anisotropy scale:** its distinctness condition is **already in its `§7` row** and **both engines
  enforce it** ⇒ ⛔ **its collapse point is NOT admitted, and its control is sound.** ⭐ Moving the
  condition into the joint premise set is **centralisation**, ⛔ not a repair.
- **Coefficient scale:** its distinctness condition is **nowhere** ⇒ ⭐ **its unit-scale collapse point IS
  admitted.** ⇒ an implementation may sit on it, satisfy **every** stated premise, and emit the baseline —
  ⛔ a control that polices nothing. ⭐ **This is the actual defect.**

⛔ **Do not write into the spec that every control's domain admits a collapse point.** ⚠ It is false for the
anisotropy package, and a builder made to write it would either state a falsehood or **weaken that
package's premise** — ⛔ which would break a control that currently works.
⚠ ⭐ **Also: there are six packages but FIVE controls** — the baseline is explicitly not a control.
⚠ ⭐ **And say why they are legitimate**, because the existing instruction *"do not add a premise to force a
solver to decide"* is being read as licence to omit them: ⭐ **a premise that keeps a CONTROL DISTINCT is
not a premise that forces a decision.** ⛔ Do not conflate them.

## ⛔ E3 — `§Q2`: delete the stated claim about the two matrix routes

It asserts a **value** for the difference between the two routes; ⚠ **both engines refute it**, and it
contradicts the two lines below it in its own section.
⛔ **Delete the VALUE CLAIM only.**

⛔⛔ **DO NOT delete the neighbouring LIMITATION** — ⚠ an earlier draft's blanket "compute and stop" would
have taken it out with the rest, and a leg caught that. ⭐ **The true and load-bearing statement is that the
two routes come from the SAME action, so their agreement is a coding-consistency check and ⛔ NOT
independent physical evidence.** ⇒ ⭐ **that is a disclosure, ⛔ not a stated result — it stays.**

## ⛔ E4 — `§Q6`: write the FOUR unwritten equations

The coefficient dimensions are said to follow from *"requiring `L` to be an energy density on a
`D`-dimensional brane"* — ⛔ **and nothing else is written.** ⇒ both engines had to **assume** the dimension
of energy, the volume element, and the two derivative operators. ⚠ They happened to assume the same four.

⛔⛔ **AN EARLIER DRAFT SAID "write all four as equations" AND DID NOT WRITE THEM.** ⚠ Both legs flagged it:
⇒ a builder would either invent them or raid a deferred document. ⭐ **They are supplied PREMISES, ⛔ not
results, so here they are — write them into `§Q6` verbatim, in the file's own `(length, time, mass)`
convention:**

```
[energy] = ( 2, −2,  1)
[L]      = [energy] · length^(−D)  = ( 2−D, −2, 1)      an energy DENSITY on a D-dimensional brane
[∂_i]    = (−1,  0,  0)
[∂_t]    = ( 0, −1,  0)
```

⚠ `[u]` is **already** supplied and stays as it is. ⛔ Nothing else in `§Q6` changes — ⚠ **except** that
`E6` still applies to any stated **result** inside it.

## ⛔ E5 — `§7`: write all six actions IN FULL

⚠ `§7` insists **every control re-enters at the ACTION** and then gives replacement **fragments**, with one
package's row being **prose about a sign**. ⇒ a reader must reconstruct the action mentally.
⭐ **Write each package's Lagrangian in full.** ⚠ It is six lines, and it is the object everything else is
computed from.

⛔⛔ **AN EARLIER DRAFT ALSO ORDERED A DECLARED SPATIAL WEIGHT HERE. ⛔ THAT IS WITHDRAWN.**
⚠ It existed **only** to serve an extraction recipe in `§Q7` that has now been removed. ⇒ ⛔ with no
quotient there is **no weight to declare**, and the factorisation question disappears with it.
⭐ **Write each action in full. ⛔ Nothing more is required of `§7` by `§Q7`.**

## ⛔ E6 — sweep the file for STATED RESULTS and delete them

⛔⛔ **DEFINITION — ⚠ an earlier draft said "what something EQUALS goes", which literally deletes the
supplied equations and the six actions this very directive ORDERS ADDED. Both legs caught the collision.**

⭐ **A STATED RESULT is an asserted OUTCOME of a computation the spec asks for** — what a residual, root,
rank, count, sign, locus, or comparison **came out to**, what is **expected**, what "the control working"
looks like, what was **measured**, or a runtime. ⛔ **Those go.**

⭐ **EXEMPT, ⛔ and these STAY:** supplied premises and definitions · the dimension equations · the six
package actions · instructions for **how** to compute · and ⭐ **disclosures that a check is weak or
vacuous by construction** — ⚠ a disclosure is the opposite of a stated result and ⛔ deleting one makes the
file *less* honest.

⭐ **`E6` OVERRIDES every other edit's "nothing else changes"** — ⭐ deleting a stated result is in scope in
**any** section, including ones otherwise out of scope. ⛔ Deleting anything else there is not.
⛔ **Delete them; ⛔ do NOT relocate them.**

⚠ **Known loci to check, ⛔ and this list is NOT exhaustive** — ⭐ read the whole file: the two-route matrix
claim; a root/sign shape early in the file; an identically-zero-root example; an exceptional-stratum answer
shape; a CAS timing figure; and a prior step's tag-parity count.

---

## ⛔ E7 — the ONE piece of `§8` that S10 actually needs. ⭐ A minimal exception, ⛔ not the rewrite.

⚠⚠ **My out-of-scope list wrongly excluded this, and a leg caught it.** `§8`'s tag grammar has **no scope
token for a stratum**, while `Q8` requires the mode count **recomputed on each allowed stratum**. ⇒ the two
engines emit **different raw names** for those tags and they match only after ad-hoc canonicalisation.
⇒ ⛔ **S10 cannot claim automatic cross-engine parity on the stratum mode counts.**

⭐ **Do the minimum that removes the ambiguity, and nothing more:** add a **stratum scope token** to the
grammar, and say where it sits relative to the other scope tokens. ⛔ **Do NOT** touch the quantity
registry, the payload serialiser, or any other part of `§8`.
⭐ **If you judge even that too wide, the acceptable alternative is a DISCLOSURE** — state in the file that
stratum tag names are not aligned across engines and that their comparison is manual. ⛔ What is **not**
acceptable is leaving it unstated, ⚠ because an unaligned name is **never compared at all** — ⭐ the
harness records it as missing rather than as a disagreement, ⇒ ⛔ the pair silently drops out of the
comparison it was supposed to be in.

## ⭐ ACCEPTANCE — ⛔ run these and paste literal output

1. ⭐ **The two-reader test on `§Q7`:** read your rewritten section as an implementer and list **every**
   object it could tell you to compare. ⛔ If the list has more than one entry, ⛔ **the rewrite has
   failed.**
2. ⭐⭐ **THE ACTION-MUTATION TEST — ⛔ this one was missing and a leg caught it.** Acceptance previously
   checked only prose. ⭐ **Take your rewritten `§Q7` and answer, in writing:** if the ACTION of one package
   is mutated while its **package selector is held FIXED**, does your text require the emitted `S_pkg` to
   move? ⭐ Quote the sentence that forces it. ⛔ If an implementation keyed on the **selector**, taking no
   action object at all, could satisfy your text, ⛔ **the rewrite has failed.**
   ⚠ ⭐ **And check your own text for a RECIPE:** if it tells the engine *how* to obtain the density rather
   than *which object* to emit, ⛔ **you have reintroduced the defect this edit removed.**
3. ⭐ **Grep the whole file** for sentences beginning "an earlier version", and for `expected`, `measured`,
   `it turns out`, `structurally zero`, `will return`, `comes back`, `seconds`, and `gaps`.
   ⛔⛔ **A grep is a floor, ⛔ not the test** — ⚠ a leg found stated results carrying **none** of these
   cues. ⭐ **Also read every section end to end** and paste the line number of every stated result you
   deleted, cued or not.
4. ⭐ **Confirm the six package actions** are each written as a complete Lagrangian, and that the
   distinctness premises appear in the supplied premise set.
5. ⭐ **Old and new line counts**, and a list of every section you touched. ⛔ Any section outside `E1`–`E7`
   appearing in that list is a scope violation — ⭐ report it.

## Report back — ⛔ under 25 lines

1. One line per `E1`–`E7`: done / partially / not.
2. The acceptance output.
3. ⭐ Anything you saw that is wrong and that you did **not** fix because it was out of scope. ⭐ **This is
   wanted** — it is how the deferred list stays honest.
4. ⛔ Do not report what any engine's values came out to be.
