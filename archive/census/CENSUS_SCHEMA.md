# The DERIVED-vs-DECLARED census — schema

**Status:** SPEC. This file defines *how* the census is taken. It does not take it.
⛔ **This file contains no census results** — no tier counts, no per-quantity verdicts, no worked
example presented as settled. Every locus quoted below is quoted as *evidence about the substrate*
(what a taxonomy already says, where a defect already sits), never as a census finding.

---

## 1. What the census is for, and the constraint that shapes it

The census answers one question: **how much of this ledger computes something rather than asserting
it.**

The result must **size tier 1** (§5). This is not bookkeeping for its own sake. The project's headline
measure is the calibrate-predict surplus — held-out matches minus tuned knobs — and that ratio has **no
denominator until tier 1 is sized**. A surplus quoted against an unsized irreducible set is not a
number.

⭐ **The governing design constraint: the schema must be able to return "almost nothing here is
derived."**

A schema whose defaults produce a flattering number is a broken schema, and it is broken in a way that
is invisible from inside the result — a high derived-fraction reads the same whether it is true or
whether the classifier simply leans that way. So the schema is built asymmetrically on purpose:

> ⭐⭐ **Every default falls AWAY FROM THE FLATTERING ANSWER — on both reported quantities, and the
> flattering direction is opposite on the two.**
>
> - **Numerator (§4 `DERIVED`): flattering is LARGE.** Where a row's status cannot be established from
>   evidence, the row is excluded from the derived set, not admitted to it.
> - **Denominator (§5 TIER 1): flattering is SMALL.** ⛔ The same defaults must therefore *not* drain
>   tier 1. Where a row's status cannot be established, the row is **not** removed from tier 1's
>   account: it lands in the tier-1 **upper bound** (§5, the reported range), never silently outside
>   the tier system.

**Why the second half is stated separately.** Unresolved states in this schema route rows out
of tiers — `A-UNADJUDICATED` has no tier (§5.3), a failed evidence requirement demotes to it (§9), an
undemonstrated convention claim is `UNADJUDICATED` (§3.4), an out-of-scope row leaves the denominator
(§7.3). Read one-directionally, that machinery is conservative for the numerator and **anti**-conservative
for the denominator: it makes the irreducible set look smaller than the evidence supports. ⭐ The
**tier-1 range** (§5) is the mechanism that closes this, and it is required, not optional.

### 1.1 ⭐⭐ The asymmetry binds ASSERTIONS, not only defaults

Everything above is written about **defaults** — what happens when the classifier stops. But a builder
does not set defaults; a builder **asserts**, and an assertion in the flattering direction walks past a
rule written about defaults without ever contradicting it. So the constraint is stated once, generally,
and it governs every classification in this schema **including ones added to it later**:

> ⛔⛔ **WHERE A CLASSIFICATION HAS AN OPTIMISTIC BRANCH AND A CONSERVATIVE ONE, THE OPTIMISTIC BRANCH
> CARRIES THE EVIDENCE BURDEN, AND AN UNMET BURDEN FALLS TO THE CONSERVATIVE BRANCH — NEVER THE
> REVERSE.**

⭐ **Reading it for a branch this schema does not name.** The optimistic branch is whichever reading
makes the ledger look **more reducible or more derived**: `A-REDUCED` over `A-UNADJUDICATED`,
`tier1-debt` over `tier1-structural`, `PHYSICS-FED` over `C-UNRESOLVED`, and — since small is the
flattering direction for a denominator (above) — **any disposition that removes a row from tier 1's
account over one that keeps it in**. ⛔ It is never settled by which reading is more likely, more
charitable to the ledger, or easier to write down; the burden sits on the optimistic side whether or
not anyone contests it.

⇒ **How to write a new evidence requirement.** Any requirement added to §9 later is written to this
rule: **state what the optimistic branch must cite, and name the conservative branch the row falls to
when it cannot.** A requirement that names only the first half is not a requirement — it is a
preference, and §9.0's fallback has nothing to fire into.

⚠ **This is why §9's requirements are not symmetric, and that is deliberate.** A rule that made both
branches cost the same evidence would leave the classifier free to pick either when neither is
evidenced, which is the state §1 exists to remove.

---

## 2. Why this is not a fourth taxonomy

Three provenance taxonomies already exist — `decisions/14`'s value classes, `parameter_register.md`'s
provenance column, stage043's `ratified_category` constants — and ⛔ **each fuses three independent
questions into one label**: where a value came from, whether a computation executed, and whether
physics entered from outside the artifact. ⇒ ⛔ **No census row inherits a substrate class, and no
precedence rule is written between the three** — a precedence rule preserves the fusion, picking a
winner among labels that each already conflate the same three questions. The census adds **orthogonal
axes** (§3) instead, and defines the user-facing tiers as a **projection** of them (§5).

---

## 3. The three axes

Every occurrence carries a value on **each** of the three axes. The axes are independent: **no value
on one implies any value on another.** A schema user who finds themselves inferring axis B from axis A
has re-created the fusion §2 exists to break.

> ⛔⛔ **NO AXIS VALUE MAY BE INFERRED FROM ANOTHER AXIS'S VALUE.** Specifically: **`B-EXECUTED` is NOT
> evidence for `A-REDUCED`.** That code ran is axis B's claim; that the quantity reduces to something
> more fundamental is axis A's, and they are **established separately or not at all**.

⛔ **This bars the blanket rule as well as the individual assertion.** *"The reduction is performed at
the binding site itself"* is not an axis-A finding — it is axis B restated in axis-A words. Applied
across a file it assigns `A-REDUCED` to every computed row at once, which makes `A-UNADJUDICATED`
**unreachable for anything the artifact executes**, and an axis value no row can take reports an
artifact of the rule rather than a measurement. ⚠ The prohibition runs in **every** pair and both
directions, not only this one: axis A implies nothing about axis B or C, axis C implies nothing about
axis A (§3.3, and §9.0 is written to match), and the `CONVENTION-LADEN` flag (§3.4) is inferred from
none of them. §9.2(1) carries the evidence consequence.

⚠ **The shape this is written against, recorded so a later reader sees the stakes.** ⛔ **Nothing here
is a census verdict on these loci** — they are cited as evidence about the substrate, in the manner of
§3.4.1 and §3.5. A `.py` artifact obtains an ordering by re-listing a dict's keys,
`order = list(harmonics)`
(`research/pde_ledger_v2/scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py:279`, over the basis
typed out at `:268-275`), while its `.wl` twin **types the identical ordering in** as a literal,
`order = {"20", "21c", "21s", "22c", "22s"}`
(`research/pde_ledger_v2/mathematica/ledger_stage016_l2_so3_covariance_mathematica_audit.wl:164`). The
quantity is the same and its content — the 5-ness — comes from the same typed basis in both engines.
⇒ Under the blanket rule the two engines land on **different axis-A values**, and the divergence is
produced **by the rule, not by the artifacts**. ⚠ A re-listing of keys is exactly the case the rule
flatters: real execution, no reduction.

### 3.1 Axis A — REDUCIBILITY

*Can this be obtained from something more fundamental?*

| value | meaning |
|---|---|
| `A-REDUCED` | it **is** so obtained; the reduction exists and is recorded |
| `A-REDUCIBLE-UNDERIVED` | a **named** route exists but has not been executed — this is reduction **debt** |
| `A-IRREDUCIBLE-STRUCTURAL` | no route exists *within this framework* |
| `A-IRREDUCIBLE-POSTULATE` | no route exists *in principle*; a defining property of the medium |
| `A-CALIBRATED` | fixed against an external physical number |
| `A-INDEPENDENT-VARIABLE` | a **coordinate**, not a knob — an independent variable of the problem (`ω`, `r`, `t`) that the artifact sweeps, integrates over, or differentiates with respect to. Nothing is tuned by choosing it and there is nothing to reduce |
| `A-UNADJUDICATED` | not established — the default (§1, §9) |

⭐ **Convention is deliberately absent from this table.** It is not an axis-A value; it is an
independent flag carried by every occurrence (§3.4). Axis A asks **only** about reducibility.

⛔ **`A-INDEPENDENT-VARIABLE` costs evidence, because it drains the denominator.** It is the one
axis-A value that removes a row from tier 1's account *without* the row being convention-laden, so
by §1's asymmetry it may not be reachable by assertion. §9 requires the loci at which the ledger
**uses** the quantity as an independent variable — swept, integrated over, or differentiated with
respect to. ⚠ Where that cannot be cited the value is `A-UNADJUDICATED` (§9.0), which keeps the row
in the tier-1 upper bound. ⛔ A symbol is not an independent variable because it *looks* like one.

#### 3.1.1 ⭐ STRUCTURAL versus POSTULATE — the test

"No route *within this framework*" and "no route *in principle*" are the two most consequential
irreducibility claims the census makes, and without a test they are distinguished by nothing but the
classifier's mood. The test is a single question, asked in this form:

> **Would an extension of the framework — one that leaves the medium's defining properties unchanged —
> open a route?**

- **Yes** ⇒ `A-IRREDUCIBLE-STRUCTURAL`. The block is a property of the framework *as currently
  written*, not of the medium. ⛔ **Name that property** (the sector not modelled, the solve deferred,
  the structure the ledger does not posit) — §9 requires it as evidence. Such a row is a standing
  statement that a specific piece of machinery is missing.
- **No, because the quantity is part of what the medium IS** ⇒ `A-IRREDUCIBLE-POSTULATE`. Deriving it
  would mean deriving the medium from something more fundamental, which this model does not posit.
  §9 requires the locus where the postulate is stated as such.
- **Cannot be decided from evidence** ⇒ `A-UNADJUDICATED`. ⚠ Note what this costs and does not cost:
  it does not remove the row from tier 1's account — it places it in the tier-1 **upper** bound (§1,
  §5). ⛔ It is never resolved by picking whichever of the two reads better.

⛔ **The two are not interchangeable, and `A-IRREDUCIBLE-STRUCTURAL` is not the safe default.** A
framework-level block is a repairable finding; a postulate is a permanent cost of the model. Recording
one as the other misstates what the ledger owes.

#### 3.1.2 ⭐ A NAMED route that requires a framework extension

§3.1's `A-REDUCIBLE-UNDERIVED` reads a named-but-unexecuted route as **debt**; §3.1.1's test reads a
route no current machinery can execute as **structural**. One named route can answer to both
descriptions at once, and they land in different tier-1 sub-buckets. The reading is fixed here:

- a named route **executable within this framework** — every piece of machinery it needs already
  exists in the ledger and what is missing is that nobody ran it ⇒ **`A-REDUCIBLE-UNDERIVED`**;
- a named route that **requires a framework extension** — it runs through a sector the ledger does not
  model, or a structure the ledger does not posit ⇒ **`A-IRREDUCIBLE-STRUCTURAL`**, and ⛔ **§9
  requires the extension itself to be named**, not merely the route.

⭐ **Why that second requirement is the point of this rule.** It makes `tier1-structural` an
**evidenced** claim rather than an assertable one. Without it, "there is a route, but it needs an
extension" is the cheapest sentence available and would absorb every awkward row — reinstating exactly
the asymmetry §9 exists to remove, where *"cannot be derived"* costs less than *"not derived yet"*.

⚠ The test is about **machinery, not effort**. "Nobody has had time to run it" is debt. "The ledger
does not model the sector the route runs through" is structural.

##### ⭐⭐ The bar for "executable within this framework"

Without one, the test above is answered by intuition, and `debt` becomes the cheap assertable claim on
exactly the split §9 already made `structural` pay for.

> ⛔ **A route is `executable within this framework` ONLY IF every quantity it consumes has a
> FUNCTIONAL FORM present in the ledger, cited by locus.**

⛔ **None of the following is a functional form**, and each has to be named because each *looks* like
one at the citation:

- a **name in a calibration or input list** — that is a label, and a label integrates nothing;
- a **dimension rule** — an `(L,M,T)` exponent vector is a statement about units, not about a
  function's shape (§3.5 is the whole reason this needs saying);
- a **free symbol** the artifact declares and never assigns (`B-DECLARED-UNASSIGNED`, §3.2) — an
  unfilled slot;
- a **literal scalar** standing where a profile is needed — one number is not the function it samples;
- an **expression assembled from any of the above** — a product of free symbols is a product of free
  symbols, however much it is shaped like an integrand.

⇒ **Where any factor's form is absent, the route is NOT executable-here.** The row goes to
`A-IRREDUCIBLE-STRUCTURAL` **with the missing form named** — and that missing form **is** the framework
extension §9 requires (§3.1.1, §9): "the ledger nowhere gives `f(·)` a functional form" is a property of
the framework as currently written, stated concretely enough to be repaired. ⚠ Where even that cannot be
established, the row is `A-UNADJUDICATED` (§3.1.1's third branch), which keeps it in the tier-1 upper
bound.

⛔ **And a record that nobody ran the route is not evidence that anybody could.** *Un-run*, *pending*,
*deferred* and *owed* are statements about **execution**; `executable within this framework` is a
statement about **machinery**. §9.2(2) states this generally, because it is not a rule about this
section.

⚠ **Why the bar sits on `debt` rather than on `structural`.** §9's requirements made *"cannot be
derived"* cost evidence; nothing yet made *"not derived yet"* cost any, and `debt` is the **optimistic**
branch — it says the reduction is available to whoever does the work. By §1.1 the optimistic branch is
the one that carries the burden, so the burden moves here.

### 3.2 Axis B — EXECUTION

*Did a computation actually run, in this artifact, to produce this value?*

| value | meaning |
|---|---|
| `B-EXECUTED` | code ran and produced the value |
| `B-DERIVED-IN-FORM-UNEXECUTED` | a symbolic derivation exists; no code path evaluates it |
| `B-DECLARED-LITERAL` | typed in |
| `B-ASSERTED-TARGET` | a value the artifact checks *against*, not one it produces |
| `B-POSTULATED` | **asserted as a premise**: no computation produces it, and none is owed (§7.1 route 2) |
| `B-DECLARED-UNASSIGNED` | the artifact **declares the symbol and never assigns it** — a free-symbol declaration |
| `B-UNADJUDICATED` | execution state not established from evidence — the axis-B default (§9) |

⚠ `B-UNADJUDICATED` exists so that a missing **code** locus demotes **axis B only**. Without it, §9's
fallback had no axis-B target and had to demote axis A instead, dropping a well-evidenced reduction
out of the tier system over a citation defect on a different axis (§9).

⛔ **`B-POSTULATED` is NOT a near-miss, and must never be read as one.** A postulate is not a failed
derivation: nothing was supposed to compute it, so there is no unexecuted route and no missing code
locus. It therefore never enters §4's `derived-in-form-but-unexecuted`, and ⛔ the absence of a code
locus on such a row is not an evidence failure under §9.0 — there is no computation to cite. ⚠ What a
`B-POSTULATED` row *does* owe is its axis-A evidence: where the premise is **stated** as one (§3.1.1,
§9).

⛔ **`B-DECLARED-UNASSIGNED` and `B-UNADJUDICATED` are opposite findings and must not be merged.**
`B-DECLARED-UNASSIGNED` is **established**: the census opened the artifact and found a declaration with
no assignment. `B-UNADJUDICATED` is **not established**: the census could not determine the execution
state at all. ⚠ Collapsing the first into the second would report a measurement as a gap.

### 3.3 Axis C — EXTERNAL INPUT

⭐⭐ **REQUIRED FIELD. Never optional, never satisfied by a prose caveat.** A row without axis C is not
a row.

**The input closure.** For an occurrence, walk back the expression that produced its value,
transitively, until every branch terminates. The set of terminals is the occurrence's **input
closure**. Tag each leaf:

| leaf tag | meaning |
|---|---|
| `C-SELF` | a literal declared in the same artifact |
| `C-PEER` | a value imported from another ledger artifact, **at a cited source locus** |
| `C-FIELDEQ` | an **equation** the value is derived *from* — a field equation, a variational principle, or the operator a solver solves (§3.3.1) |
| `C-STATED` | a value a document merely **asserts** — a stated model postulate quoted as a number, a value read out of a note's table |
| `C-EXTERNAL` | a measured physical number |
| `C-MATH` | a mathematical constant or a pure-mathematical identity |
| `C-CONVENTION` | a **convention the artifact adopts** — a time arrow, a sign convention, an orientation or a normalisation choice — which fixes part of the value, characteristically its **sign** |
| `C-FREE` | a **free symbol carrying no value** — the artifact never assigns it, and nothing it stands for has been fixed |

⛔ **`C-STATED` is split out of `C-FIELDEQ` on purpose, and it does NOT confer `PHYSICS-FED`.** An
equation is something a value follows *from*; a stated value is something a document *says*. Folding
the second into the first would make **a document asserting a number a valid physics leaf** — which is
exactly the assertion-counted-as-derivation this census exists to detect. ⚠ This does not weaken
`A-IRREDUCIBLE-POSTULATE`: where a postulate is *stated*, that statement is evidence on **axis A**
(§9), not a physics leaf on axis C.

⛔ **`C-PEER` requires a cited source locus.** A number that matches another artifact's number, with no
citation, is **hand transcription**, and the leaf is `C-SELF`. Without this rule a chain of
transcription propagates `PHYSICS-FED` transitively and re-admits §4's near-miss 3 one hop downstream.

⭐ **Whose citation establishes a `C-PEER` leaf.** It may live in **either** of two places, and ⛔ **the
row records which**:

- **`peer-cited-in-artifact`** — the locus is carried by the artifact that consumes the value;
- **`peer-cited-in-stage-note`** — the locus is carried by **that artifact's own stage note**.

⛔ **Three conditions, none of them relaxed by this widening.** (i) It must be a **locus** — a file and
a line. ⛔ **Naming a source *stage* is not a citation**, in either place; that is the case the rule
above sends to `C-SELF`. (ii) The note must be **that artifact's** stage note, not any note in the
corpus that happens to mention the quantity — otherwise the leaf is established by whatever document a
classifier chooses to go looking in. (iii) §9.1 applies unchanged: the locus is **opened and read**,
and a note-carried attribution that does not check out is a finding, not a pass.

⭐ **The two populations are reported SEPARATELY, and may not be summed into one `C-PEER` count.** A
locus carried by the consuming artifact is evidence that *the artifact* imported the value; a locus
carried only by its note is evidence that *a document about the artifact* says where the value came
from — weaker, and weaker in a way that is invisible once the two are added together. ⚠ Where a
`PHYSICS-FED` verdict rests transitively on a chain of `C-PEER` hops, the row records the **weakest**
hop kind in the chain: one note-carried hop makes the whole chain note-carried.

⛔ **`C-CONVENTION` never confers `PHYSICS-FED`, and it does NOT set the `CONVENTION-LADEN` flag.** A
convention is a choice the artifact makes, not physics entering from outside, so it is a terminal like
any other declaration. ⚠ **And it is not a shortcut into §3.4**: the flag is set by §3.4's clauses
(a)/(b)/(c) and by nothing else, so a `C-CONVENTION` leaf makes an occurrence a *candidate* the §3.4
test must be run on — never a convention by leaf tag. Reading it the other way would hand back exactly
the author-labelling §3.4 replaced. ⭐ **What the leaf is for:** where a claimed **sign** data-depends
on one, the row records it, because a sign fixed by a convention is not a sign the model predicts —
and §3.4's clause (b) is then answerable rather than notional, since flipping the convention is a
transformation whose effect on the claim can actually be checked.

⚠ **`C-SELF` versus `C-MATH`, when a leaf reads as both.** A hand-typed literal that happens to equal a
mathematical constant satisfies both descriptions. The test: **could the artifact have declared a
different value and still be internally consistent?**

- **Yes** ⇒ `C-SELF`. The artifact *chose* the number; that it coincides with a mathematical constant
  is not provenance.
- **No — the value is forced by a mathematical identity the artifact cites or executes** ⇒ `C-MATH`.
- **Cannot be decided** ⇒ `C-SELF`, the conservative direction (§1): it keeps §3.5's shape — code
  running over the artifact's own declarations — visible rather than dissolving it into mathematics.

⚠ Neither tag confers `PHYSICS-FED`, so this precedence does not move the headline; it decides what the
`self-referential` and near-miss buckets are able to show.

⛔ **`C-FREE` never confers `PHYSICS-FED`.** A symbol nothing assigns has fed nothing in, and ⛔ it is
not `C-SELF` — nothing was declared. ⚠ **What a closure containing one is:** where no other leaf
establishes physics-feeding, the closure is **`C-UNRESOLVED`**, not `PHYSICS-FED = false`. The walk
terminated but the terminal carries no evidence either way, and that is not the positive finding
`unclassified-nonfed` records (§5.1); by §1's asymmetry the row must stay in **tier 1's account**
rather than leave the tier system — at whichever bound its **axis A** places it (§5.6, and the
`C-UNRESOLVED` rule below). Where another leaf *does* establish physics-feeding, a `C-FREE` leaf
changes nothing.

**Two derived fields, computed from the closure:**

- **`PHYSICS-FED`** — three-valued: `true` / `false` / `C-UNRESOLVED`. It is `true` iff the closure
  contains a `C-FIELDEQ` or a `C-EXTERNAL` leaf, **or** contains a `C-PEER` leaf that is itself
  `PHYSICS-FED`. Evaluate transitively. A `C-PEER` cycle yields `C-UNRESOLVED`.
- **`SELF-REFERENTIAL`** — true iff the walk reaches the occurrence's own declaration **through at least
  one intermediate step**: a walk that leaves the occurrence and returns to its origin. ⛔ A closure that
  terminates *immediately* at the occurrence's own declaration is a plain literal, and is
  `SELF-REFERENTIAL = false`.
  ⚠ **What counts as an intermediate step.** The round trip must pass through **at least one other
  occurrence** — a different `(artifact, binding-site, quantity)` triple (§7.2) — and the row **records
  which**. ⛔ A sub-expression *within the same binding site* is not an intermediate step: rearranging a
  literal algebraically on one line is still a plain literal, and reading it as a round trip would fill
  the `self-referential` bucket with formatting. ⭐ The consumer is §3.5's shape — a value handed back
  to itself **through the ledger's own machinery** — and machinery means another binding site.

⛔ **Rule, not recommendation: if the closure cannot be determined, the field is `C-UNRESOLVED` and
the row is NOT counted as derived.** Exclusion is the default. An undetermined closure is never
resolved by assuming the favourable branch. ⚠ **Excluded from the DERIVED set — never from tier 1**;
which of the two an unresolved closure touches is settled immediately below.

(`C-UNRESOLVED` is a value of `PHYSICS-FED`, the closure-level field — not a leaf tag. A closure
containing an untraceable leaf is `C-UNRESOLVED` as a whole, and `SELF-REFERENTIAL` is then undefined
rather than false.)

⭐ **Where `C-UNRESOLVED` lands.** For every tier test **and every §4 bucket test** it counts as
**¬`PHYSICS-FED`** — such a row can never reach TIER 3 (§5), and ⛔ it is never asserted into §4's
`executed-but-not-physics-fed`, which is a positive claim the evidence does not support. It does ⛔
**not** go to `unclassified-nonfed` (§5.1) either: that bucket is a
**positive** finding (the closure *was* determined and contained no physics), and an unresolved closure
is not that finding.

> ⛔⛔ **AXIS C GATES TIER 3 AND THE `DERIVED` HEADLINE. IT NEVER REMOVES A ROW FROM TIER 1.**

An unresolved closure blocks `tier3-emergent` and blocks `DERIVED` (§4). It ⛔ does **not** overturn an
**evidenced** axis-A value: a row carrying an evidenced `A-REDUCIBLE-UNDERIVED`,
`A-IRREDUCIBLE-STRUCTURAL`, `A-IRREDUCIBLE-POSTULATE` or `A-CALIBRATED` keeps that tier and its
sub-bucket, whatever axis C could not resolve. ⇒ `C-UNRESOLVED` sends a row to `unadjudicated` (§5.3)
**only where axis A is itself unestablished** — and for `A-REDUCED`, whose only tier is tier 3, that is
the entirety of what the gate does.

**Two independent reasons. Either alone is sufficient:**

1. **The axes answer different questions (§3).** Axis A asks whether the quantity can be reduced; axis
   C asks whether physics entered *this occurrence*. An unresolved closure says nothing about
   reducibility — a named-but-unexecuted route is still a named-but-unexecuted route, and a stated
   postulate is still a stated postulate. ⛔ Letting C overturn A re-creates exactly the fusion §2
   exists to break.
2. ⭐ **§1's two-directional asymmetry requires it.** Tier 1 is a **denominator**, where *small* is the
   flattering direction. Routing an **evidenced** debt row out of tier 1 and into `unadjudicated`
   **understates** tier 1 — it spends an axis-C gap on shrinking the very count the census exists to
   size. The conservative direction, the one §1 mandates, keeps the row in.

⇒ §5.7's enum and §9.0's fallback are written to match. ⚠ §5.6 is where the consequence shows up: this
rule moves evidenced rows out of the range's **upper span** and into its **lower bound**. It narrows
the span by *establishing* rows, never by assuming them.

⭐ **What `SELF-REFERENTIAL` is for — its consumer.** ⛔ An occurrence whose `SELF-REFERENTIAL` is
**true** is **never `A-REDUCED`**: the closure returned the value it was handed. Its axis-A value is
whatever a non-self branch supports, and where there is none it is `A-UNADJUDICATED`. ⚠ A bare literal
is not blocked by this rule — it is `SELF-REFERENTIAL = false` (above), and its axis A is adjudicated
normally. Self-referential occurrences are additionally reported in a named bucket
**`self-referential`** with its own count (§5). ⚠ §3.5's worked shape is what this field is built to
catch.

#### 3.3.1 Axis C for values no expression produced

"Walk back the expression" presumes an expression. Four production modes have none, and each of them
flips `PHYSICS-FED` — and hence tier 3 — depending on which reading a classifier picks. The readings
are fixed here.

1. ⭐ **Solver / numerical-integration / root-find output.** The **operator** — the equation the routine
   solves or integrates — is a `C-FIELDEQ` leaf of the closure, **in addition to** the literal inputs,
   ⛔ **iff the artifact cites the locus of the equation being solved.** Where no model equation is
   cited, the operator is a generic numerical kernel and contributes no leaf: the closure is then just
   the literal inputs, and it is typically `{C-SELF, …}`. ⚠ Test that this can fail: a solve assembled
   entirely from the artifact's own constants, with no cited equation, is **not** `PHYSICS-FED`.
   ⚠ **A generic linear-algebra routine is NOT an operator.** `rank()`, `NullSpace`, a determinant, an
   eigen-solve and their kin compute a **property of a matrix**; they solve no model equation, so ⛔ the
   routine contributes **no leaf** under this rule, cited or not. ⭐ What can be a `C-FIELDEQ` leaf is
   **what was assembled into the matrix** — the constraint rows — each judged individually and each
   under the same citation condition. ⚠ This is not pedantry: a rank is a fact *about the rows*, so the
   physics content of the answer is exactly the physics content of the rows and no more, and attaching
   the leaf to the routine would let an uncited assembly inherit `PHYSICS-FED` from a library call.
2. **Fit output.** The fitted **target** is the leaf: `C-EXTERNAL` if it is a measured number,
   `C-PEER` if it is another artifact's value at a cited locus, `C-SELF` otherwise. The fit *form* is a
   `C-FIELDEQ` leaf only under the same citation condition as (1). ⚠ A fit is also almost always
   `A-CALIBRATED` on axis A — the axes are independent and both must be recorded (§3).
3. **Transcription.** Covered by the `C-PEER` citation rule above: cited ⇒ `C-PEER`, uncited ⇒
   `C-SELF`.
4. ⭐ **A constitutive proposition (§7.1 route 2).** Nothing produced it — it is a premise. Its closure
   terminates at **where the model posits it**: that terminal is `C-STATED` where a document states it,
   or `C-SELF` where the artifact itself declares it (§3.3). ⛔ Either way the row is **not**
   `PHYSICS-FED`. ⚠ That is not a demotion: the statement of a
   premise is evidence on **axis A** (§9), which is where a postulate is meant to score, and §3.3's
   `C-STATED` split is exactly what keeps it off axis C. Its axis B is `B-POSTULATED` (§3.2).
   ⭐ **The dependence hop.** Where a reported result *depends on* such a proposition and no expression
   consumes a value of it, the walk takes a **dependence hop**: the proposition is a terminal of that
   result's closure, tagged as above, and the row records the hop kind as `premise-dependence` in its
   reachability witness (§7.1.1). ⛔ **The hop confers no `PHYSICS-FED` on anything it passes through.**
   It puts the premise **in the universe**, which is the whole of what route 2 claims.

### 3.4 The `CONVENTION-LADEN` flag — an independent field, not a fourth axis

⭐ **REQUIRED FIELD, boolean-plus-evidence, carried by every occurrence.** Values: `true`, `false`,
`UNADJUDICATED`.

**Why it is a flag and not an axis-A value.** A value can be reducible **and** convention-laden at
once, so making the two mutually exclusive loses information. `ℓ_P = √(ħG/c³)` is *reduced from
calibrations*; `ℓ_P = 1` in Planck units is *convention*. Those are different occurrences of the same
symbol, and axis A must be able to say the first is reduced without being blocked by the second.

With `A-CONVENTION` gone there is no axis-A home for a bare units pin. Such an occurrence is
adjudicated on axis A like any other — and where no route is established it is `A-UNADJUDICATED`
(§5.3), not silently exempt. It is **this flag, never axis A**, that removes an occurrence from the
tier totals (§5.2) and quotients it out of the rank denominator (§10.1).

**The operational test — this replaces author labelling.** An occurrence is:

- **CONVENTION** ⇒ flag `true`, iff **(a)** a documented transformation group exists for the units /
  gauge / coordinates / normalisation / pinning at issue; **(b)** an alternative admissible value,
  applied consistently to every affected quantity, leaves the action, the equations, the boundary data
  and **every dimensionless observable the MODEL states** invariant; and **(c)** no external datum is
  consumed.
- **CALIBRATION** ⇒ flag **`false`**, iff an external benchmark selects a coordinate on the quotient
  space, and changing it at fixed convention **changes a dimensionless physical relation**. (A
  calibration is not convention-laden; it is a spent tuning, and it stays counted — §5.2, §10.2.)
- **NOT A CANDIDATE CONVENTION** ⇒ flag **`false`**, the ordinary case: no transformation group is
  claimed for the occurrence at all. This — not `UNADJUDICATED` — is the value for any quantity nobody
  proposes is a units / gauge / coordinate / normalisation choice.
- ⛔ **`UNADJUDICATED` is reserved for a row where a convention claim IS made or implied and clauses
  (a)/(b) are not demonstrated — never convention by intuition, and never by author label.**

⭐⭐ **Clause (b)'s quantifier is over the MODEL's observables, not over the corpus's current reach.**
⛔ **Vacuous invariance is not invariance.** If the reason nothing varies is that **nothing computes
anything** from the quantity — a sector whose simulation is deferred, a branch not yet solved — clause
(b) is **untested, not satisfied**. Concretely:

> ⛔ **A quantity from which no dimensionless observable of the model is currently computable is
> `UNADJUDICATED`. It is NEVER `CONVENTION`.**

Without this the flag is **easiest to earn exactly where the model is least developed**, and the
laundering hazard of §3.4.1 — the one clause (b) exists to close — walks straight through the clause
built to close it. ⚠ Clause (b) demands a demonstrated invariance of *something*; trivial invariance of
an empty set is not a demonstration.

⛔ **An `UNADJUDICATED` convention claim buys no exclusion.** The row keeps whatever tier axis A gives
it, and is additionally reported in a named bucket **`convention-unadjudicated`** with its own count.
Exclusion shrinks tier 1, so by §1's asymmetry the undemonstrated case must fall toward *not*
excluded. ⛔ **That count reads `is_tier` only** (§6.1). ⚠ §9's evidence requirement demotes **the
flag**, never axis A — see §9.0.

#### 3.4.1 The laundering hazard the test exists to block

⚠ Recorded as a hazard about the substrate. ⛔ **Nothing here is a census verdict on those rows.**

In `software/em_charge_attribute/g0_closure_card_v0.md`, the collar rounding `δ/a` = `1/20` is tagged
`[CONVENTION]`, dimension `1` (`:44`), while `ℓ/a = ε_r` = `1/20` — two rows above it in the **same
six-row table** — is tagged `[CALIBRATED]`, dimension `1` (`:42`). Two dimensionless `1/20` ratios in
one table, classed oppositely, the distinction resting on the prose "reuse `ℓ`, so no new scale".

⚠ A **dimensionless** ratio is what the standing rule says *tests* the model — so mislabelling one as
convention would quotient away a testable prediction. That is the laundering route, and clause **(b)**
is what closes it.

### 3.5 Why axis C exists — stated as method

Consider a dimension-checking walk of the following shape. A recursive `dim_of` descends an
expression; at every symbol it terminates in a single lookup into a dimension table; and that table is
built by copying a module-level map of hand-typed exponent triples. Every leaf of every closure it
produces is then `C-SELF`.

Such a walk is `B-EXECUTED` — code genuinely ran, and it can genuinely fail. But its closure is
`{C-SELF}`, so `PHYSICS-FED = false`. **Code ran; no physics entered.** Without axis C the two are
indistinguishable, and every executed check reads as a derivation.

The mechanism is real and present in the corpus, which is why the axis is required rather than
advisory. In `research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py`,
`dim_of`'s leaf lookup is `symbol_dims[expr]` (`:436-439`), and the dimension map handed to it is built
as `dims = dict(SOURCED_DIMS)` (`:505`), `SOURCED_DIMS` being a typed table of exponent triples.

⛔ **Nothing in this paragraph is a census verdict about stage023.** Classifying that stage's
occurrences is the census's work, not this spec's. The loci are cited to show the *shape* the axis
detects.

---

## 4. The headline quantity

There is **exactly one** headline definition, and no looser one is permitted anywhere in the census
output, its summaries, or any document that quotes it:

> **DERIVED := `B-EXECUTED` ∧ `PHYSICS-FED` ∧ `A-REDUCED`.**

Three near-misses fail this conjunction. ⛔ **Each is reported as its own named bucket and is never
folded into DERIVED**, in a summary, a headline, or a rounding:

1. **executed-but-not-physics-fed** — `B-EXECUTED ∧ ¬PHYSICS-FED`. Code ran over the artifact's own
   declarations. Real work; not a derivation from the model.
2. **derived-in-form-but-unexecuted** — `B-DERIVED-IN-FORM-UNEXECUTED`. The symbolic route exists and
   nothing evaluates it.
3. **physics-fed-but-declared-literal** — `PHYSICS-FED ∧ B-DECLARED-LITERAL`. A value derived from
   physics *elsewhere* and **stored here as a literal**. The derivation is real and belongs to
   whichever artifact performed it; it is ⛔ not this artifact's derivation, and counting it here would
   double-count the corpus's derived fraction.

⛔ **A `B-POSTULATED` row is in NONE of the three** (§3.2). Each near-miss is a **derivation that fell
short** — one that ran on the wrong inputs, one that was never evaluated, one performed elsewhere. A
premise is not a derivation that fell short; nothing was owed. ⚠ Filing postulates under near-miss 2
would report the model's inputs as unfinished work, which is the opposite of what tier 1 records.

---

## 5. The tier projection — the output shape

The tiers are a projection of axis A (with one stated axis-C guard on tier 3, §5.1), taken over the
occurrences the `CONVENTION-LADEN` flag does not remove (§3.4, §5.2). Every occurrence in this list
carries the projection **twice** — once as `is_tier` and once as `should_be_tier` (§6); ⛔ the counts
below are always `is_tier` counts. The field's permitted values are fixed at §5.7.

⚠ **The tier list is not the census's whole output.** Buckets defined in §3.3, §3.4, §4, §5.1–§5.3,
§7.1.1, §7.3, §8.4 and §10.3 are equally required. ⭐ **The complete required output set is enumerated
at §5.8**; a census that reports only the three tiers is incomplete.

- **TIER 1 — IRREDUCIBLE.** ⭐ **Reported as a RANGE, never a scalar** (§5.6). It **splits three ways,
  and the split is the point of the deliverable** —
  "irreducible" collapses three very different situations that score differently:
  - `tier1-debt` ← `A-REDUCIBLE-UNDERIVED` — a reduction someone could actually go and do.
  - `tier1-structural` ← `A-IRREDUCIBLE-STRUCTURAL` — not derivable within this framework.
  - `tier1-postulate` ← `A-IRREDUCIBLE-POSTULATE` — a genuine postulate about the medium.
- **TIER 2 — CALIBRATED** ← `A-CALIBRATED`.
- **TIER 3 — EMERGENT** ← `A-REDUCED` ∧ `PHYSICS-FED`. **It splits two ways, and the split is
  required (§5.5):**
  - `tier3-calibration-propagated` — follows algebraically from tier-1 and tier-2 values.
  - `tier3-held-out` — follows, **and** is not implied by any calibrated input, so it can be compared
    against a number nobody tuned.

### 5.1 Consequence 1 — reduction without physics is not emergence

⛔ `A-REDUCED ∧ ¬PHYSICS-FED` is **not** tier 3. It goes to a named holding bucket,
**`unclassified-nonfed`**, which is reported with its own count.

A value that "follows" only from other declarations has been shown to follow from *the artifact*. It
has not been shown to follow from *the model*. Tier 3 is a claim about the model.

### 5.2 Consequence 2 — convention sits outside the tiers

**Convention-laden occurrences (§3.4) sit outside all three tiers**, in their own reported bucket,
**whatever their axis-A value**. They are excluded from every tier total and reported as their own
count.

✅ **RATIFIED (user, 2026-07-30): convention outside all three tiers.**

⭐ The decisive argument is that **the corpus already counts this way**. Folding convention into
tier 2 would not be a schema choice; it would silently change an existing counting rule:

- `notes/parameter_register.md:33` gives `` `CONV` | convention / units / pin choice — never free |
  no ``, while `:35` gives `` `CALIB` / `CALIB-ANCHOR` | calibrated to a benchmark / anchor
  (magnitude tuned) | **yes** (a tuned knob) `` — the third column being "Counts toward irreducible?"
  (header `:29`).
- `notes/stages/ledger_stage043_irreducible_count_range.md:69` states the coherent counting rule as
  `count = nominal − DERIVED-and-EXECUTED − CONV − external-bridge`, reinforced at `:78`:
  "**IMPOSED calibrations STAY COUNTED** — a tuning is NOT a reduction".

Folding convention into tier 2 would additionally inflate tier 2 by exactly the set of quantities that
carry no physical content, which is the opposite of what tier 2 is for.

⚠ **Reconciliation with the standing rule — which this does not overturn.** The standing project rule
phrases the absolute constants `G, c, ħ, ℓ_P` as "convention **and** calibration". That is true **at
the parameter level**, where one symbol carries both roles. The census is at the **occurrence** level
(§7.2), where the two separate cleanly: the measured anchor and the unit choice are different
occurrences of the same symbol. ⇒ The census **refines** the standing rule at finer grain; it ⛔ does
**not** overturn it.

### 5.3 Rows with no tier

`A-UNADJUDICATED` has no tier by construction. Such rows go to a named bucket **`unadjudicated`**,
reported with its count alongside the tiers, and are excluded from every derived-fraction numerator.
A large `unadjudicated` bucket is an honest census result, not a failure of the census.

⭐⭐ **But they are NOT excluded from tier 1's account.** By §1's second half, an unresolved row may not
drain the denominator: `unadjudicated` **is** the tier-1 upper bound's span (§5.6). ⛔ There is no
state in this schema in which a row's status "could not be established" and the row silently leaves
the irreducible count. The demonstrated case the range exists to catch is a bare literal knob with no
named route, no benchmark and no stated postulate — an object that satisfies **no** §9 evidence row,
lands here, and under a point estimate would have vanished from the very count it belongs to.

### 5.4 ✅ RESOLVED — tier 3 versus near-miss 2

§4 requires **derived-in-form-but-unexecuted** to be reported as its own bucket and never folded in.
§5's tier-3 rule (`A-REDUCED ∧ PHYSICS-FED`) does not mention axis B, so a row that is
`A-REDUCED ∧ PHYSICS-FED ∧ B-DERIVED-IN-FORM-UNEXECUTED` satisfies tier 3 while also being a §4
near-miss. Tier 3 is therefore a **superset** of DERIVED, and the unexecuted rows sit inside it.

✅ **RATIFIED (user, 2026-07-30): TIER 3 is defined WITHOUT the `B-EXECUTED` conjunct.** It stays
`A-REDUCED ∧ PHYSICS-FED`. The tier definition is *"everything that follows from tiers 1 and 2"* — a
claim about **logical consequence**, not about whether a code path ran. A symbolic derivation that
nothing evaluates still follows.

Two consequences, both rules:

1. ⭐ **`DERIVED` (§4) is a strict sub-count of TIER 3, never a synonym.** The difference between them
   is exactly the unexecuted rows.
2. ⛔ **TIER 3 is always reported split by axis B.** This was a mitigation while the question was
   open; it is now a **rule**, so the unexecuted subset can never be invisible inside a tier-3
   headline.

### 5.5 ⭐ The TIER 3 split — calibration-propagated vs held-out

**Required.** Every tier-3 report carries the two sub-counts named in §5:

- `tier3-calibration-propagated` — follows algebraically from tier-1 and tier-2 values;
- `tier3-held-out` — follows, **and** is not implied by any calibrated input, so it can be compared
  against a number nobody tuned.

Both "emerge" algebraically, but ⭐ **only the held-out subset earns predictive surplus**. A tier-3
total that mixes them cannot support the surplus claim the census exists to enable (§1).

⛔ **The undivided tier-3 total may not be quoted as evidence of prediction.**

⇒ Tier 3 therefore carries **two** required splits at once — by axis B (§5.4) and by
calibration-propagation (here). They are independent; neither substitutes for the other.

### 5.6 ⭐⭐ TIER 1 is reported as a RANGE. A bare tier-1 scalar may not be quoted.

**Required.**

> **TIER 1 = `[ confirmed-tier-1 , confirmed-tier-1 + unadjudicated ]`**
>
> - **lower bound** = occurrences (resp. QIDs) whose tier-1 assignment is *established from evidence*;
> - **upper bound** = the lower bound plus everything the census could **not** adjudicate —
>   `unadjudicated` (§5.3), including rows demoted there by an **axis-A** evidence failure (§9.0) or an
>   undecidable STRUCTURAL/POSTULATE call (§3.1.1).

⚠ **`unadjudicated` remains the upper span; what changed is which rows are in it.** A `C-UNRESOLVED`
closure **alone no longer widens the span** (§3.3): axis C gates tier 3 and the `DERIVED` headline,
never tier-1 membership, so a row carrying an **evidenced** axis-A value together with an unresolved
closure sits in the **lower bound**. Such a row reaches the upper span only where its **axis A** is
itself unestablished — including `A-REDUCED ∧ C-UNRESOLVED`, whose only tier is tier 3 and which
therefore has no tier at all (§5.7). ⭐ The span narrows by rows being **established**, never by rows
being assumed.

⛔ **And nothing enters either bound by assertion.** A row the pass did **not classify** is in neither
bound. A claim that such a row *"would land in tier 1"* is a classification claim without a
classification (§7.3.1) and ⛔ does not place it — at either bound, or in the span between them. ⚠ The
range moves only by rows the census actually adjudicated, and the count of rows it did not reach is
reported separately (§5.8, outputs 15 and 16) rather than folded into the span as a guess.

⛔ **A bare tier-1 scalar may not be quoted** — not in the census, not in a summary, not in anything
that cites the census.

**Why a range and not a point.** Every unresolved state **on axis A** routes a row *out* of a tier
(§1), and small is the flattering direction for a denominator. A point estimate would therefore report the *most
flattering* reading of its own unresolved rows as if it were the measurement. The range makes the
unresolved set **visible in the headline instead of absorbed by it**: a wide span is a statement about
how much the census could not establish, which is itself a result (§5.3).

⭐ **This is the corpus's own existing idiom, not a new convention.** The stage043 note
(`notes/stages/ledger_stage043_irreducible_count_range.md:14-18`) states its Decision A as "**the count
is REPORTED AS A RANGE, not a single headline**", carries the continuous codimension as `[40, 49]` with
its spread stated, and reports a separate discrete count alongside it rather than folding the two.

⚠ **Level consistency.** §10.1 requires every reported fraction to carry both an occurrence and a QID
denominator. The range is computed **at each level from that level's own units**, with the §10.2
aggregation rule applied first — so a QID that already has a confirmed tier-1 occurrence sits in the
QID-level *lower* bound and is ⛔ not added again to the QID-level upper bound by its unadjudicated
occurrences.

### 5.7 The `is_tier` / `should_be_tier` value enum

Both fields (§6) take values from **this list and no other**. It is total: every occurrence gets one.

| value | assigned when |
|---|---|
| `tier1-debt` | `A-REDUCIBLE-UNDERIVED` |
| `tier1-structural` | `A-IRREDUCIBLE-STRUCTURAL` |
| `tier1-postulate` | `A-IRREDUCIBLE-POSTULATE` |
| `tier2-calibrated` | `A-CALIBRATED` |
| `tier3-emergent` | `A-REDUCED` ∧ `PHYSICS-FED` |
| `no-tier:unclassified-nonfed` | `A-REDUCED` ∧ `PHYSICS-FED = false` (§5.1) |
| `no-tier:independent-variable` | `A-INDEPENDENT-VARIABLE` (§3.1) — a coordinate, not a physical input |
| `no-tier:convention` | `CONVENTION-LADEN = true` (§5.2) — whatever axis A says |
| `no-tier:unadjudicated` | `A-UNADJUDICATED` — including where §9.0's per-axis fallback demoted **axis A**, and where `A-REDUCED` fails tier 3's `PHYSICS-FED` conjunct by an unresolved closure (§3.3) |

⛔ **`PHYSICS-FED = C-UNRESOLVED` is NOT, by itself, a line in this table.** It gates `tier3-emergent`
and it gates `DERIVED` (§4); it does not place a row. A row with an evidenced axis-A value and an
unresolved closure takes that axis-A value's line — `tier1-debt`, `tier1-structural`,
`tier1-postulate` or `tier2-calibrated` — exactly as it would with a resolved one (§3.3, §9.0). The one
case where an unresolved closure decides the line is `A-REDUCED`, whose only tier is tier 3: it fails
that tier's conjunct, and having no other tier to fall to, lands at `no-tier:unadjudicated` (⛔ never at
`no-tier:unclassified-nonfed`, which is a positive finding — §5.1, §3.3).

⚠ **Precedence, since a row can match two lines.** `CONVENTION-LADEN = true` is evaluated **first** and
wins: a convention-laden occurrence is `no-tier:convention` whatever axis A says (§5.2, ratified). Only
where the flag is `false` or `UNADJUDICATED` does the axis-A projection run — and an `UNADJUDICATED`
flag ⛔ buys no exclusion, so such a row takes its axis-A value normally (§3.4).

⭐ **The tier-1 sub-bucket is part of the value, never a separate optional column** — §10.2's
aggregation rule needs it present to detect an intra-tier-1 conflict.

⚠ **Tier 3's two required splits are NOT enum values.** By axis B (§5.4) and by calibration-propagation
(§5.5) — they are independent of each other, so neither can be folded into a single tier token. They
are carried as their own required fields on every `tier3-emergent` row.

⛔ **The bare words `IRREDUCIBLE` / `CALIBRATED` / `EMERGENT` are section headings, not field values.**
Writing a bare `IRREDUCIBLE` into `should_be_tier` discards the sub-bucket the deliverable is about.

### 5.8 ⭐ The complete required output set

The tier list above is **not** the whole output. All of the following are required, each reported with
its own count, and ⛔ none may be omitted or folded into a neighbour:

| # | output | defined at |
|---|---|---|
| 1 | **TIER 1**, as a range, split three ways | §5, §5.6 |
| 2 | **TIER 2** — at QID level, the **calibrated** count, which overlaps the others by design | §5, §10.2.1 |
| 3 | **TIER 3**, split by axis B **and** by calibration-propagation | §5, §5.4, §5.5 |
| 4 | **`DERIVED`** — the strict sub-count of tier 3 | §4, §5.4 |
| 5 | near-miss **`executed-but-not-physics-fed`** | §4 |
| 6 | near-miss **`derived-in-form-but-unexecuted`** | §4 |
| 7 | near-miss **`physics-fed-but-declared-literal`** | §4 |
| 8 | **`unclassified-nonfed`** | §5.1 |
| 9 | the **convention** bucket | §5.2 |
| 10 | **`unadjudicated`** | §5.3 |
| 11 | **`convention-unadjudicated`** | §3.4 |
| 12 | **`self-referential`** | §3.3 |
| 13 | the **conflict set** | §10.3 |
| 14 | the **stage043 delta**, all three directions: census QIDs with no `REG:` ID, `REG:` IDs reconciled as out-of-scope, and `REG:` IDs reconciled as **`unvalued-in-universe`** | §8.4 |
| 15 | the **out-of-scope list**, each row with its reason — with the **provisional** (paused / pending) exclusions carrying their own sub-count | §7.3 |
| 16 | the **reached-by-no-reported-result** count | §7.1.1, §7.3 |
| 17 | the enumerated **reported-result set** itself, with loci | §7.1.1 |
| 18 | **`unvalued-in-universe`** — reachable, valued by no artifact; ⚠ QID-level only, zero occurrences by construction | §8.4 |
| 19 | the **`C-PEER` populations, reported SEPARATELY** — `peer-cited-in-artifact` and `peer-cited-in-stage-note`; ⛔ never summed | §3.3, §9 |
| 20 | **`no-tier:independent-variable`** — reported because it is the one axis-A value that leaves tier 1's account without the convention flag | §3.1, §5.7 |

⚠ **These do not partition anything.** Buckets 5–7 overlap the tiers by construction, and after §10.2
a QID can be counted both in tier 3 and in the calibrated total. ⛔ **Do not sum this table.**

---

## 6. ⭐⭐ Every occurrence carries IS **and** SHOULD-BE

§5 says how one tier assignment is made. **Every occurrence carries two of them**, and the pair — not
either column alone — is what makes the census answerable (§6.2).

| field | meaning |
|---|---|
| `is_tier` | what the evidence currently supports. **This is the measurement.** Values: §5.7. |
| `should_be_tier` | what the physical picture says it *ought* to be. Values: §5.7. |
| `delta` | flagged wherever the two differ |

⭐ **The set of deltas IS the revisit list**, and that is the point of the field: it says where to go
back and look.

**`should_be_basis` is ALWAYS REQUIRED** — on every occurrence, whether or not the two tiers differ —
and must be one of:

| basis | meaning |
|---|---|
| `named-route` | a specific reduction is identified and could be executed — **actionable debt** |
| `physical-picture-expectation` | the model's picture implies it should reduce, but no route is named yet. A hunch, and ⛔ it must stay legible as one. |
| `convention-candidate` | the picture says the occurrence is a units / gauge / coordinate / normalisation choice, but §3.4's clauses (a)/(b) are **not demonstrated**. Then `should_be_tier` = `no-tier:convention` while `is_tier` stays whatever axis A gives. |
| `none` | **the default.** No basis is claimed — and then `should_be_tier` **must** equal `is_tier`. ⛔ A should-be that departs from `is_tier` on no basis is not recorded. |

⛔⛔ **`convention-candidate` sets nothing and excludes nothing.** It is a `should_be_tier` value, so
§6.1's guard already bars it from every numerator and every denominator — and that is precisely why it
is safe to record. The `CONVENTION-LADEN` flag is set by §3.4's clauses and by nothing else; an
undemonstrated convention claim leaves the flag `UNADJUDICATED`, buys **no** exclusion, and the row
keeps its axis-A tier (§3.4, §9.0). ⭐ **The value exists so the belief is written down where it can be
audited instead of acted on** — §3.4.1's laundering route runs through exactly this belief, and a
schema with no place to put it invites the classifier to express it by setting the flag.

⚠ **Which of X11's two options this is.** The field was previously required *only* when the tiers
differ, which made `none` — whose own definition requires them equal — unreachable. It is resolved by
making the field **always required with `none` as its default**, rather than by deleting the value: an
explicit `none` on an agreeing row is the difference between "considered and no basis claimed" and
"never looked at", and only the first is a census row.

### 6.1 ⛔⛔ THE GUARD

> **Every headline count uses `is_tier`. NEVER `should_be_tier`.**

Without this the census becomes quotable as *"the real figure is much lower, once you count what we
expect to derive"* — which is exactly the flattering number §1's asymmetry rule exists to prevent.
⚠ The guard is stated without naming any expected magnitude; what the picture expects is at §14, and
⛔ it is not for the classifier.

⭐ **`should_be_tier` is a work list. It is never a discount, never a subtraction, and never appears in
a numerator or a denominator.**

### 6.2 Why the two columns are what make the census answerable

A tier-1 count on its own **cannot** distinguish two opposite conclusions. Whatever tier 1 comes in
at, either:

- most of those rows carry a `should_be_tier` of `tier3-emergent` ⇒ the ledger has a large **reduction
  debt**, and the deltas are the work list; or
- most carry a `should_be_tier` still in **tier 1** ⇒ **the medium requires more specification than
  the physical picture suggests**, which is a finding about the model itself.

⛔ **Do not assume which.** ⭐ **That gap IS the result**, and the `is`/`should-be` pair is the
instrument that reads it. A census reporting only `is_tier` leaves the question unanswered.

⚠ **No statement of the expected magnitude appears in §3–§9.** What the physical picture expects tier 1
to be is deliberately **not** stated anywhere a classifier reads while assigning rows; it is addressed
to the **reader of the finished result**, at §14. A census-taker who knows the expected answer before
classifying is not measuring.

---

## 7. Universe, occurrence, inclusion

### 7.1 Universe

The census ranges over what the ledger **asserts** and the physics **depends on**. There are **two
admission routes**, and a row admitted by either is in the universe:

> **Route 1 — valued quantities.** Any named quantity to which the ledger assigns a value — numeric or
> closed-form symbolic — on which the physics depends.
>
> **Route 2 — constitutive propositions.** A proposition the model *posits* about the medium or the
> framework, which is **not derived from anything else in the ledger**, and on which a reported result
> depends. Admitted **whether or not it carries a numeric or symbolic value**.

⚠ **"On which the physics depends" is a description, not a test.** The operational test that decides
membership is at **§7.1.1**, and it is the only one — it decides **both** routes, route 2 by the
`premise-dependence` hop of §3.3.1(4). Nothing below §7.1.1's rule may be used to admit or exclude a
quantity.

#### ⭐⭐ Why route 2 exists — and what it is NOT

A structural postulate is **label-valued**: a boundary-condition class, a time arrow, an existence
claim. Under route 1 alone no expression walk can ever reach one, because there is no value to walk
back to. ⇒ Under route 1 alone, `tier1-postulate` (§5) is a bucket **no row can enter** — and a bucket
nothing can enter reports an artifact of the admission rule, never a measurement. ⛔ That is not an
acceptable state for the one sub-bucket the model's own picture expects to carry its inputs (§3.1's
`A-IRREDUCIBLE-POSTULATE` exists for precisely these), so the rule is widened rather than the finding
suppressed.

**Route 2's three bounds, and they are not optional** (§7.2, §7.3). The proposition must be:

1. **depended on by a reported result** (§7.1.1);
2. **not derived elsewhere in the ledger** — if something derives it, it is a conclusion, not a premise;
   and
3. ⭐⭐ **not a restatement of a freedom the universe already carries via route 1.**

⛔ An assertion failing test 1 or test 2 is **out of scope**, recorded with its reason (§7.3). Without
bounds 1 and 2 route 2 would admit every sentence in the corpus, and the universe would become a reading
of prose rather than a closure. Bound 3 is different in kind and is stated in full below: a proposition
failing it is not out of scope — it is **already in the universe**, under another row.

#### ⭐⭐ Bound 3 — route 2 may not re-count a route-1 freedom

> ⛔ **A proposition admitted by route 2 must not restate a freedom already counted via route 1.**

**Before admitting, check** whether the proposition's **content** is already carried by a **valued
quantity in the universe**. Where it is, the assertion is **that occurrence's evidence, not a new
occurrence** — it is recorded on the route-1 row (⛔ never silently dropped) and route 2 admits nothing.

⚠ **The check is on content, not on wording.** A proposition and a valued quantity restate **one**
freedom whenever **fixing either fixes the other**. A premise that a degree of freedom exists, and the
value the ledger assigns to that degree of freedom, are two statements of one thing; admitting both
counts one input twice.

⭐ **One source literal asserting several things is ONE occurrence** — ⛔ **unless each part is
independently a freedom the model takes, and that independence is EVIDENCED, not assumed.** How many
keys an author typed into one literal is a formatting decision; reading each as its own premise
multiplies `tier1-postulate` by that decision. ⚠ Where independence cannot be shown, the parts are one
occurrence and the surplus parts are that occurrence's **evidence**.

⚠ **Reconciling bound 3 with §1, which is the objection to read first.** De-duplication **shrinks**
tier 1, and §7.3 says in as many words that shrinking the denominator is exactly where a census goes
wrong without anyone quoting a number. ⇒ **So the fold costs evidence too, in the same shape as
everything else in §9:** a fold must **name the route-1 quantity** whose freedom the proposition
restates, **at its locus**, and show that fixing either fixes the other. ⛔ **A fold on suspicion, on
similarity of wording, or on "that sounds like the same thing" is not a fold** — the proposition stays
admitted as its own occurrence. ⭐ What §1.1 forbids is an **unresolved** state falling toward the
flattering answer; a **demonstrated** identity between two rows is not an unresolved state, and §10.1
already states that two rows can encode one degree of freedom and that the rank, not the row count, is
the answer.

> ⛔⛔ **ROUTE 2 WIDENS THE UNIVERSE. IT RE-TAGS NO AXIS-C LEAF AND CONFERS NO `PHYSICS-FED`.**

⚠ **The §3.3 split between `C-FIELDEQ` and `C-STATED` stands exactly as written, and route 2 must not
be read as loosening it.** They answer different questions:

- **`C-FIELDEQ`** is a **premise the model posits**, in the form of an equation, a variational
  principle or a solver's operator — something a value is derived *from*. A value that follows from one
  is `PHYSICS-FED`.
- **`C-STATED`** is a **document asserting a value**. It confers **no** `PHYSICS-FED`, ⛔ however
  central the document, however the document labels the value, and ⛔ whether or not the proposition the
  document is discussing is itself admitted by route 2.

⇒ **Admitting a proposition does not convert a number into a physics leaf.** A note stating `x = 1/20`
is `C-STATED` even where the note calls `x` a postulate; what route 2 admits is the *proposition* the
model posits, as a **row of its own** to be classified (axis A by §3.1.1, axis B `B-POSTULATED`, axis C
by §3.3.1(4)) — not a promotion of anybody's stated number. ⛔ Collapsing the two would make a document
asserting a value a valid physics leaf, which is the assertion-counted-as-derivation this census exists
to detect (§3.3, §4 near-miss 3).

⭐ **Two live instances in the corpus, named so this rule is checkable.** ⚠ Cited as evidence that
route 2 has referents — ⛔ **neither is classified here; the census does that.**

- stage043's **11** `discrete-structural-postulate` IDs: the category constant at
  `research/pde_ledger_v2/scripts/ledger_stage043_irreducible_count_range_sympy_audit.py:143`, its
  expected count at `:172`, reached from the status mapping `"STRUCTURAL-POSTULATE"` at `:158`.
- stage023's **`pathA29_premise`**, carried at
  `research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:425`
  (`Z_is_premise`, `boundary_dof`) and `:600-602`, and described by the artifact itself at `:1051` as
  "the keystone return premise" — a proposition, consumed by no expression.

⛔ **The universe is NOT the dimension-record universe.**

This has to be said explicitly, because the corpus already offers three tempting candidate
enumerations, all inside a single stage. In
`research/pde_ledger_v2/notes/stages/ledger_stage023_nullspace_underdetermination.md`:

- **42** source objects under the stage's declared membership rule (rule stated `§1.6`, `:116-124`,
  whose criterion is that the object's *value is an `(L,M,T)` exponent vector*; the count is carried
  machine-readably at `:352`, `<!-- ENUM_COUNT|row_count|42 -->`, and referenced in prose at
  `:205-207`);
- **29** emitted records (the `dimension_records` map,
  `research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:539-571`);
- **34** physics-routed objects (`§1.7(1)`, `:379-390`).

**All three enumerate objects whose value is an `(L,M,T)` exponent vector.** A derived-vs-declared
census over exponent vectors would measure the dimension bookkeeping, not the physics.

⇒ Those three counts are evidence that **the universe question is hard**. They are ⛔ **not a menu of
candidate universes** for this census, and the census must not adopt any of them as its scope.

⚠ Everything above rules candidates **out** — and §7.5's survey figure is likewise ruled out, as an
estimate of effort rather than a scope. What follows is the rule that decides what is **in**. Without
it this section excludes every enumeration the corpus offers and names no replacement, and a
census-taker cannot select row 1.

#### 7.1.1 ⭐⭐ The positive membership rule

> **A named quantity — or a constitutive proposition (§7.1 route 2) — is IN THE CENSUS UNIVERSE iff it
> lies in the transitive input closure (§3.3) of at least one REPORTED PHYSICAL RESULT of the ledger.**

The closure machinery is the same one axis C already defines — including §3.3.1's rules for values no
expression produced. Nothing new is invented here; the universe is the axis-C closure taken from a
fixed starting set. ⚠ **Route 2 adds a hop, not a second test:** a proposition enters the closure by
the `premise-dependence` hop of §3.3.1(4), so it is admitted by **this** rule like everything else, and
⛔ it carries the same reachability witness every other in-universe row carries.

**What counts as a reported physical result.** A value or relation the ledger presents as a statement
about the **modelled world**, rather than about its own bookkeeping. An entry qualifies iff all three
hold:

- **(a)** it is presented as an **outcome**, in a stage note, a paper section or a script emission —
  not as an intermediate on the way to one;
- **(b)** it is a claim about the medium or its phenomenology — a dimensionless ratio compared or
  comparable to an observation, the form / sign / power of a force law, an equation of motion or one of
  its coefficients, a spectrum, a stability statement about the medium, a predicted relation between
  physical quantities;
- **(c)** ⭐ **the discriminator: it would still be a claim about the world if every bookkeeping
  artifact of this ledger were deleted.**

⛔ **Not reported physical results** — these are the ledger's statements about *itself*: dimension
exponent vectors and dimensional CORRECT/WRONG/UNDETERMINED verdicts (§11), coverage figures and row
counts, gate pass/fail flags, ID manifests and reconciliation deltas, category tallies, and ⛔ the
census's own outputs. ⚠ Test (c) is what separates these: delete the bookkeeping and they say nothing.

> ⛔⛔ **A RESULT'S LOCUS IS THE SITE OF THE CLAIM'S VALUE-BEARING CONTENT.**

A result is entered with the locus at which the claim's **content** lives — the expression, assertion
or emission carrying the claimed value, form, sign, power or relation. ⛔ **A verdict token, a gate
flag, a PASS/FAIL tally or a status string is NEVER a result locus.** By this section's own test it is
a claim about the ledger rather than about the world, and the exclusion list above already names it.

⭐ **Why this is a rule and not a matter of style.** The result set is the **starting set** of every
closure walk (step 2 below), so each locus decides what the walk pulls in. A verdict token's value
characteristically data-depends on **every guard that feeds it** — so entering one as a locus admits an
artifact's whole gate apparatus to the universe on the strength of a formatting decision about which
string the script prints last. The universe would then be fixed by how a script *reports*, not by what
it *claims*, and it would grow by whatever a later author adds to the ladder.

⚠ **What a verdict token may still be used for.** It remains legitimate evidence under test **(a)** —
that a claim is *presented as an outcome*, and that the artifact separates its own status flag from its
physics. ⛔ What it may not be is the closure walk's **entry point**.

**⭐⭐ The ordering, stated because it is the whole point.**

1. ⛔ **FIRST**, enumerate the reported-result set — **as a listed artifact, one entry per result, each
   with its locus.** It is a required census output in its own right (§5.8, output 17).
2. ⛔ **THEN** derive the universe from it, by walking each result's closure and **marking** every
   object reached, together with **which result reached it and by what hop** — the row's
   *reachability witness*, a required field on every in-universe row.
3. Only then classify.

> ⛔ **The universe is derived from the result set. The result set is NEVER derived from the universe,
> and it is never amended once classification has begun.** A result set that can be adjusted after the
> rows are seen is a knob on the denominator.

⚠ Step 2's walk and the axis-C determination are the **same walk, done once** — the closure is what
both need. What the ordering constrains is that the **starting set is fixed before the walk begins**.

**⭐ How someone holding one stage's script decides membership.** By the row's **reachability witness**,
which step 2 records per binding site — so the marking pass emits, per artifact, that artifact's slice
of the universe: the in-universe binding sites in that file and the result each answers to. A holder of
one script consults its slice.

⛔ **Membership is NOT decidable by reading one script in isolation, and this rule does not pretend
otherwise.** Whether a quantity feeds a reported result is a question about what *consumes* it, and
nothing in the file that declares it can answer that. ⚠ A locally-decidable membership criterion is
precisely what produced the three enumerations above — "is this object's value an `(L,M,T)` exponent
vector?" is answerable from one file, and it enumerates the bookkeeping. **The cost of a correct
universe is one corpus-wide pass; the schema pays it rather than accepting a local proxy.**

**⚠ The consequence, stated honestly.** An object reachable from **no** reported result is **out of
scope**. It is recorded as such — with its locus and the reason `reached-by-no-reported-result` — and
⛔ **never quietly dropped**: the bucket carries its own count (§5.8, output 16, and §7.3).

⛔ **If that set is large, that is itself a census finding and must be reported as one.** A ledger most
of whose value-bearing objects reach no reported result is a statement about the ledger, and the census
is the instrument that would show it. Suppressing the count would delete the finding.

### 7.2 Occurrence

**One occurrence = one `(artifact, binding-site, quantity-or-proposition)` triple** at which an
artifact assigns or asserts a value, **or asserts a proposition constitutive of the model** (§7.1
route 2).

- The same quantity at **two binding sites in one artifact** is **two occurrences**.
- The same quantity in **two artifacts** is **two occurrences**.
- The same proposition asserted in **two artifacts** is likewise **two occurrences** — a premise
  restated in a second place is a second assertion, and whether the two agree is exactly what §10.3 is
  for.

⛔ **The proposition half is BOUNDED, and the bound is checked per occurrence, not per proposition.**
An assertion is an occurrence only if it satisfies **all three** of §7.1 route 2's tests: the
proposition is **depended on by a reported result**, it is **not derived elsewhere in the ledger**, and
it **does not restate a freedom already counted via route 1**. ⚠ An assertion failing either of the
first two is **out of scope under §7.3**, recorded with its reason — ⛔ never silently dropped and ⛔
never admitted "because it sounded foundational". Without this bound the occurrence definition would
swallow every assertion in the corpus, and §1's denominator would be set by prose volume.

⚠ **Bound 3 narrows this triple-rule for propositions, and the narrowing is deliberate.** A single
literal asserting several propositions would otherwise be as many occurrences as it has parts, by the
same reading that makes one quantity at two binding sites two occurrences. §7.1's bound 3 makes it
**one** occurrence unless each part's independence is evidenced — ⛔ and a proposition folded under
bound 3 is not an out-of-scope exclusion: it is recorded as evidence on the route-1 occurrence that
already carries the freedom.

*Operational definitions (needed to execute; see §13 for their status):*

- **artifact** = one committed file the ledger treats as a source of record — a stage's `.py`, a
  stage's `.wl`, a stage note, `parameter_register.md`, a `paper/stages/*.tex`. The `.py` and `.wl` of
  one stage are **two artifacts**, so a dual-engine value has at least two occurrences. This matters:
  it is exactly the case where one engine may compute what the other types in.
  ⚠ **Locus style.** The engines are `scripts/ledger_stage0NN_<slug>_sympy_audit.py` and
  `mathematica/ledger_stage0NN_<slug>_mathematica_audit.wl`; ⛔ every locus in a census row carries its
  **own full path**, never a bare `:NNN`.
- **binding-site** = the file-and-line of the assignment or assertion expression. A multi-line
  expression takes the line of its head. **One markdown table row is one binding site.**
  ⚠ **Wolfram globals**, since a `.wl` artifact binds differently from a `.py` and the rule must be
  decidable without judgement. In a file under `ClearAll["Global`*"]`
  (`research/pde_ledger_v2/mathematica/ledger_stage023_nullspace_underdetermination_mathematica_audit.wl:10`),
  a global's binding site is the line of the **top-level `Set` / `SetDelayed` head** that assigns it:
  1. a symbol assigned at **several** top-level lines has **one binding site per assignment**, exactly
     as in a `.py` — ⛔ the last assignment does not absorb the earlier ones;
  2. a symbol assigned only **inside** a `Module` / `Block` / `With` is **local**: it is not a global
     binding site, and it belongs to the enclosing definition;
  3. a **pattern definition** (`f[x_] := …`, e.g. `:38`) binds a rule, not a value — its binding site
     is the definition line, and it is an occurrence only where the artifact names its result as a
     quantity;
  4. a **system global the artifact sets** is a binding site like any other — it is set by this
     artifact and it can carry model content. `$Assumptions` at `:26-36` of that file is one, and its
     head line is its binding site under the multi-line rule above.

### 7.3 Inclusion

**In scope:** values the physics depends on.

**Out of scope** — each exclusion **recorded with its reason so it is auditable**, ⛔ never silently
dropped:

- **`reached-by-no-reported-result`** — the §7.1.1 universe rule did not reach it. ⭐ This exclusion
  carries its **own reported count**, separately from the ones below (§5.8, output 16), because it is
  the one exclusion that is a finding about the ledger rather than a housekeeping decision;
- **`derived-proposition-not-constitutive`** — an asserted proposition a reported result depends on
  which **is** derived elsewhere in the ledger, so it fails §7.1 route 2's second bound. It is a
  conclusion, not a premise. ⚠ Route 2's *first* bound needs no separate reason: an assertion no
  reported result depends on is already `reached-by-no-reported-result` above;
- controls and deliberate negatives (constructions built to prove a check able-to-fail);
- test scaffolding;
- display and formatting constants;
- loop bounds and tolerances;
- **`retired-or-excluded`** — values belonging to an approach the ledger has **retired or excluded**;
- ⭐ **`paused-or-pending`** — values belonging to an approach the ledger records as **PAUSED**, or an
  obligation it records as **PENDING**. ⛔ **This is a DIFFERENT reason from `retired-or-excluded` and
  must not be recorded as one.** The corpus keeps the two apart and this list previously did not: the
  register carries the R49 field-profile rows as living on "the **PAUSED** EM/U1 Phase-C parallel
  track" and discharged by a BVP marked "(R49, **PENDING**)"
  (`notes/parameter_register.md:355`, and `:887`), and the stage043 note uses the same word,
  "the **PAUSED** EM/U1 Phase-C parallel track"
  (`notes/stages/ledger_stage043_irreducible_count_range.md:222`).
  ⭐ **A paused item is expected to RETURN.** Its exclusion is therefore **provisional**: it is marked
  as such, it carries its own sub-count inside the out-of-scope list (§5.8, output 15), and ⛔ it is
  **re-checked** on any later pass rather than treated as settled. ⚠ A retired approach is a decision
  the ledger has taken; a paused one is a decision it has deferred, and reporting the second as the
  first would retire work nobody retired.

⛔ **Out of scope removes a row from the denominator (§10.1), so this list is exactly where a tier-1
denominator can be shrunk without anyone quoting a number.** Every exclusion carries its reason and is
auditable; ⚠ an exclusion reason that cannot be checked against the row is not a reason. ⭐ The
provisional reason above is the sharpest case: it shrinks the denominator on the strength of a
deferral, so it is the one exclusion that must be able to come back.

#### 7.3.1 ⛔⛔ A claim about where an UNCLASSIFIED row would land is a classification claim

> ⛔ **Saying of a suppressed, deferred or out-of-scope row that it *"would land in tier 1"* — or in any
> other tier or bucket — asserts an axis-A value for a row the pass did not classify, and it costs the
> same §9 evidence as making the classification.**

⛔ **It may not be asserted** — not in a row file, not in a summary, not in a note about either, and not
as a reassurance that an exclusion did not cost anything. ⚠ It is available in exactly **one** form:
**classify the row**, with its evidence, under §9. Until then the honest statement is the **count** of
rows the pass did not reach and the **reason** it did not reach them — which §7.3 and §5.8 (outputs 15
and 16) already require.

⭐ **Why this is not pedantry about phrasing.** Such a claim is read as a statement about the tier-1
range while carrying **none of the range's evidence**: it moves §5.6's span in the reader's head without
appearing in either bound, in whichever direction the writer already believes. ⚠ And nothing downstream
would catch it being wrong, because it is attached to **no row** — there is no locus to open (§9.1), no
axis to demote (§9.0) and no conflict to detect (§10.3). ⇒ It is an assertion in precisely the sense
§1.1 governs, made where the schema's whole apparatus cannot reach it.

### 7.4 ⭐ Every row carries a LIVE / RETIRED flag

**Required field.** The reason is measured, not hypothetical: `parameter_register.md` carries retired
rows — `λ_Pu` (`:139`) and `α_aniso` (`:159`), both retired by Decision 16 — which still list their
historical dimensions, and the register has **no live/retired field** (see
`research/pde_ledger_v2/notes/measure_register_sufficiency.md:131`). The rows are
not unmarked — each says "**RETIRED (Decision 16)**" in its name cell and "dropped from the live knob
set" in its notes — but the marking is **prose inside cells, not a field**. A census that ingested
those rows mechanically would count retired physics into a live denominator.

### 7.5 Expected scale — an estimate, flagged as one

A prior survey reports the corpus at **226** distinct dimension-bearing quantities across the **30 of
43** scripts that carry dimension machinery, against **145** register-defined quantities spanning
**31 of 43** stages (`notes/measure_register_sufficiency.md:27-29`, coverage line `:32`).

⛔ **It bounds expected effort and is not re-derivable** (the scripts named as producing it are absent
from the tree and from git history); it is **not an input to any census number**.

---

## 8. Identity, IDs, and the join

### 8.1 The join rule

Occurrences join into one quantity by **adjudicated identity**. ⛔ **Never by candidate name.**

Both halves of the evidence for this are measured:

- `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md` §"GROUPING LIMITATIONS" (`:190-207`) lists **eight**
  known-same cross-stage pairs that its general rule does **not** group, because their initial case
  differs after separator-boundary normalization (`CsSquaredDim`/`cs_squared_dim`,
  `CorruptKDim`/`corrupt_K_dim`, `EnergyDim`/`energy_dim`, `FourVolumeDim`/`four_volume_dim`,
  `MassDim`/`mass_dim`, `OmegaDim`/`omega_dim`, `PressureDim`/`pressure_dim`, `RhoDim`/`rho_dim`).
  The same file shows the converse hazard: `T_w`→`TW` at stage016 (`:126`) never meets `Tw`→`Tw` at
  stage004 (`:51`) and stage013 (`:108`), and the grouped row at `:186` accordingly carries only the
  `Tw` pair. ⚠ Those two are *not* a missed merge — they carry different triples under different
  conventions — which is precisely why the decision is an adjudication and not a string match.
- The parameter register has **no machine key at all**: "no `quantity_id`; names collide and one
  quantity has several spellings" (`notes/measure_register_sufficiency.md:119-131`, gap 3 at `:124`), over "one flat
  namespace over a corpus that is not" (gap 5, `:127`).

### 8.2 IDs

**Extend stage043's existing scheme rather than inventing one.** That scheme is
`REG:<bucket-or-class>[:<part>]:<name>` — e.g. `REG:a:hbar`, `REG:E:IV:c_E` and `REG:C2:Ztilde0`, at
`scripts/ledger_stage043_irreducible_count_range_sympy_audit.py:216`, `:271` and `:377`.

Census IDs:

- **quantity ID:** `QID:<adjudicated-name>`
- **occurrence ID:** `QID:<adjudicated-name>@<artifact>#<binding-site>`

### 8.3 The registry

**One file:** `research/pde_ledger_v2/notes/census/QID_REGISTRY.md`. One row per QID, carrying its
adjudicated members and **the adjudication's evidence**.

- The **builder mints** IDs.
- ⭐ The **physics review leg adjudicates identity**, because a naming decision is a physics decision:
  saying two symbols are the same quantity is a claim about the model, not about strings.

⚠ **The two roles are the project's, not this spec's** — defined in `docs/development_pipeline.md`
(Roles table, `:45`, `:49`, and the review leg at `:109`). In short: the **builder** is the single party
that writes the deliverable — here, the census rows and this registry; the **physics review leg** is an
independent **fresh** agent that reviews it and whose verdict is **blocking**. ⛔ **They are never the
same party**, and ⛔ the builder does not adjudicate its own merges — that is the whole content of the
rule above.
- ⛔ **An unadjudicated merge is not a merge.** It stays two QIDs, and it is reported as two, flagged
  as an open identity question.

### 8.4 Reconciliation with stage043 — required, and mechanical

- **Every one of the script's 152 `REG:` IDs must be RECONCILED**, in exactly one of **three** ways:
  1. it **maps to ≥ 1 census occurrence**; or
  2. it is **recorded as out of scope**, carrying **its stage043 category** (the `ratified_category`
     constants, `scripts/ledger_stage043_irreducible_count_range_sympy_audit.py:137-146`) **and a
     reason** referred to §7.1.1 or §7.3; or
  3. ⭐ **Route C — `unvalued-in-universe`**: it is **reachable from a reported result** and **assigned
     a value by no artifact in the ledger**.

  ⛔ **Route C is neither of the other two, which is why the binary was not exhaustive.** It is not out
  of scope — §7.1.1's walk reached it. And it maps to no occurrence — §7.2 requires a binding site that
  assigns or asserts, and no artifact ever provides one. Under a two-route rule such an ID has nowhere
  to go and is either forced into a wrong bucket or left as an unlisted remainder.

  ⭐ **The consequence is substantive rather than clerical, and must be stated as one:** **a quantity
  the model uses and never pins is a free parameter the model takes as input.** ⇒ It **counts toward
  TIER 1**. Its axis-A value is decided by §3.1.1's STRUCTURAL / POSTULATE test, or is
  `A-REDUCIBLE-UNDERIVED` where a route is named — carrying §9's evidence in either case. ⛔ **It is not
  `unadjudicated` on account of having no value.** That nothing in the corpus pins it is a **positive**
  finding, established by the same corpus-wide pass §7.1.1 already pays for, and it is exactly that
  finding which places the quantity in tier 1.

  ⚠ **Two mechanics, so the route neither double-counts nor disappears.**
  - **A Route-C entry is QID-level with ZERO occurrences**, by construction. It enters the QID-level
    counts and the QID-level tier-1 range (§5.6, §10.2.2) and contributes **nothing** to the
    occurrence-level ones. ⛔ Its absence from the occurrence denominator is not a licence to drop it
    from the QID one — that is the silent shrink §7.3 warns about. §10.1's both-denominators rule is
    what keeps the asymmetry visible instead of hiding it.
  - **The tier-1 sub-bucket still costs evidence.** Where §3.1.1's test cannot be decided, its third
    branch fires as written: the sub-bucket is `A-UNADJUDICATED` and the entry sits in the tier-1
    **upper** bound. ⚠ That does not contradict the rule above — Route C establishes that the quantity
    is inside **tier 1's account**, which holds at either bound; what the evidence decides is *which*.
- **All three directions are reported as named deltas with counts** (§5.8, outputs 14 and 18): census
  QIDs mapping to no `REG:` ID — the census's genuine extension — the option-2 out-of-scope
  reconciliations, and the Route-C `unvalued-in-universe` set. ⛔ None is absorbed silently or left as
  an unlisted remainder.
- ⛔ **Do not re-derive the 152-ID manifest.** It exists; reconcile with it.

⚠ Stage043's own "engine-qualified" manifest guard is `all(identifier.startswith("REG:") …)` (script
`:437`, `:444`), which is **always true** for these IDs. It is not a check the census may rely on.

---

## 9. Evidence every row must cite

The census must be **auditable rather than an opinion table**. Per row, the following evidence is
**required**, by axis value:

| when the row asserts | it must cite |
|---|---|
| axis B `B-EXECUTED` | the **code locus of the computation**, and the loci of its **input leaves** |
| axis C `PHYSICS-FED` | the locus of the **field equation or external anchor** actually reached; for a solver or fit, the locus of **the equation being solved or the form being fitted** (§3.3.1) — ⛔ **not the locus of the call that solves it**; a routine is not an operator (§3.3.1, §9.2(1)) |
| axis C `C-PEER` leaf | the **source locus** the value was imported from, **and which of `peer-cited-in-artifact` / `peer-cited-in-stage-note` carries it** (§3.3) |
| axis A `A-REDUCED` | ⭐⭐ **a DEFINING RELATION expressing this quantity in terms of other NAMED quantities**, where that relation is recorded, **and the loci of those quantities** — ⛔ **pointing at the computation that produced the value does NOT satisfy this; that is axis B** (§3, §9.2(1)). Direction-checked per §9.2(3) |
| axis A `A-REDUCIBLE-UNDERIVED` | the **named route**, where that route is recorded, **direction-checked** (§9.2(3)), **and that it is executable within this framework** — every quantity the route consumes carrying a **functional form present in the ledger, cited by locus** (§3.1.2). ⛔ A record that the route is un-run/pending/deferred is not evidence of this (§9.2(2)) |
| axis A `A-IRREDUCIBLE-STRUCTURAL` | ⭐ **which framework property forecloses the route** (§3.1.1) — **and where a route IS named, the framework EXTENSION it requires** (§3.1.2), the route **direction-checked** (§9.2(3)). ⚠ A **missing functional form** named under §3.1.2 satisfies the extension requirement |
| axis A `A-IRREDUCIBLE-POSTULATE` | where the postulate is **stated** as a defining property (§3.1.1) |
| axis A `A-CALIBRATED` | the **benchmark** it is calibrated against — the external value and **its** locus; ⛔ the fit or solve that consumed it is axis B (§9.2(1)) |
| axis A `A-INDEPENDENT-VARIABLE` | ⭐ the loci at which the ledger **VARIES** it — a sweep over **more than one value**, an integration over it, or a differentiation with respect to it (§3.1). ⛔ **That an executed expression consumes the symbol is not this**; consumption is axis B, and every input of every computation is consumed (§9.2(1)) |
| axis B `B-POSTULATED` | where the model **posits** the proposition (§7.1 route 2); ⛔ no code locus is owed — there is no computation |
| axis B `B-DECLARED-UNASSIGNED` | the **declaration** locus, and that the artifact contains **no assignment** of it |
| a route-2 row | ⭐ **all three** of §7.1's bounds: **which reported result depends on it**, that it is **not derived elsewhere** in the ledger, and that it **does not restate a freedom already counted via route 1** (§7.1 bound 3). ⚠ A **fold** under bound 3 costs its own evidence — the route-1 quantity **named at its locus** |
| `CONVENTION-LADEN = true` | the **transformation group** (§3.4 clause a) and where the **invariance** (clause b) is demonstrated over the model's stated observables |
| every in-universe row | its **reachability witness** (§7.1.1) |
| `should_be_tier ≠ is_tier` | the `should_be_basis` (§6), and for `named-route`, where the route is recorded — **direction-checked** (§9.2(3)) |

⭐⭐ **Why `A-REDUCED` and `A-IRREDUCIBLE-STRUCTURAL` carry requirements.** Both were previously
assertable bare, and both decide the deliverable:

- `DERIVED` (§4) is a three-way conjunction. With `A-REDUCED` unevidenced it had two evidenced
  conjuncts and one free one — the fallback below could never fire on it, and §9.1's open-the-locus
  rule had no locus to open.
- Calling a row `tier1-debt` cost a named route; calling it `tier1-structural` cost nothing. ⛔ That
  made **"cannot be derived" the CHEAPER claim than "not derived yet"** — on precisely the tier-1 split
  the census exists to produce. Both now cost evidence.
- ⭐ **And an evidence requirement is only as strong as the object it names.** `A-REDUCED`'s requirement
  was first written as *"the reduction — where it is performed"*, which any executed binding site
  satisfies: the row could be evidenced by **pointing at the computation**, so the requirement demoted
  nothing and axis A was decided by axis B. It now names the **relation**, not the run (§9.2(1)).

### 9.0 ⭐ The fallback is PER-AXIS

**A row that fails to carry required evidence is demoted — but only on the axis whose evidence is
missing.** The three axes and the flag are independent (§3); a citation defect on one is not evidence
about another.

| what is missing | what it demotes | what it does NOT touch |
|---|---|---|
| axis-A evidence | axis A → `A-UNADJUDICATED` | axes B, C; the flag |
| axis-B evidence | axis B → `B-UNADJUDICATED` | axes A, C; the flag |
| axis-C evidence | `PHYSICS-FED` → `C-UNRESOLVED` | axes A, B; the flag; ⭐ **and the row's TIER, where axis A is evidenced** |
| convention evidence | the flag → `UNADJUDICATED` | **axis A — untouched** |
| the reachability witness | the row is out of universe (§7.1.1), recorded with its reason | — |

⛔ **A whole-row demotion is not permitted.** A row with a well-evidenced `A-REDUCED` and an uncited
code locus keeps its axis-A value and its tier; it loses only its execution claim.

⭐⭐ **The axis-C case, stated in the same form, because it is the one that decides tier 1.** A row with
a well-evidenced `A-REDUCIBLE-UNDERIVED` and an unresolved closure **keeps `tier1-debt`**; a row with a
well-evidenced `A-IRREDUCIBLE-STRUCTURAL`, `A-IRREDUCIBLE-POSTULATE` or `A-CALIBRATED` likewise keeps
its tier. What such a row loses is **tier 3 and the `DERIVED` headline** — and that is the whole of what
axis C gates (§3.3, §5.7). ⛔ **An axis-C failure never routes an evidenced row into `unadjudicated`.**
Doing so would let a citation defect on the provenance axis shrink the irreducible **denominator**,
which is the anti-conservative direction §1 forbids — and it would restate axis C's finding as if it
were a finding about reducibility, which is the fusion §2 breaks.

⚠ The single exception is not an exception to that rule but a consequence of tier 3's definition:
`A-REDUCED` projects **only** to tier 3 (§5.7), so a row that is `A-REDUCED` with an unresolved closure
has no tier left to keep and lands at `no-tier:unadjudicated`.

⚠ **§3.4 governs the convention flag, not §9.** An undemonstrated convention claim demotes **the flag**
to `UNADJUDICATED` and leaves axis A alone — the row keeps whatever tier axis A gives it and is
additionally reported in `convention-unadjudicated` (§3.4). ⛔ It does **not** wipe axis A. Exclusion
shrinks tier 1, so by §1 the undemonstrated case must fall toward *not* excluded, and a rule that
answered "we could not demonstrate this is convention" by deleting the row's reducibility would do the
opposite.

⛔ A demoted value is not a guess, not a "probably", and not a provisional classification to be firmed
up later by the same person who could not evidence it.

### 9.1 ⭐⭐ A cited locus must be where the value is COMPUTED

Not where a document attributes it. The reason is measured.

`notes/measure_register_sufficiency.md:98-102` records the register's rows `:182`–`:185` as attributing their
dimensions to "stage 017, dual-engine", while **stage017 has no dimension machinery in either
engine** — `research/pde_ledger_v2/scripts/ledger_stage017_grouped_p2_lane_isotropy_sympy_audit.py:62`
carries only `CITED_016_DIMENSIONAL_OK = True`, mirrored as `cited016DimensionalOk = True` at
`research/pde_ledger_v2/mathematica/ledger_stage017_grouped_p2_lane_isotropy_mathematica_audit.wl:34`
— and the values actually live in stage016's `make_dim_rules` at
`research/pde_ledger_v2/scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py:314-325`, the twelve
literal exponent triples. **The values agree; the pointer does not.**

⚠ **Current state of that defect, which the census must check for itself rather than inherit from the
survey:** register rows `:182`, `:183` and `:184` have since been corrected in place — each now names
the stage016 script and states outright that "stage017 has no dimension machinery in either engine".
Row `:185` still reads "(stage 017, dual-engine: …)". The survey's `:182–:185` framing is therefore
stale for three of the four rows.

⇒ **Two rules follow:**

1. **Census evidence is verified by opening the cited locus.** Never by trusting an attribution — not
   the register's, not a stage note's, and not this spec's.
2. Where a row's own substrate attributes it to a stage that does not compute it, the census records
   **both** loci — the attributed one and the computing one — and **marks the attribution defect**.

### 9.2 ⭐⭐ Three things that are NOT evidence

Each of these has carried a classification past the table above, and each fails for a **different**
reason. They are stated as **general rules over the whole table**, not as notes on the rows that
exposed them: a rule attached to one row is walked around by the next row that needs it, and §1.1
requires these to govern requirements added later as well.

**(1) ⛔⛔ AN EXECUTION IS NOT AXIS-A EVIDENCE.** Restating §3 in the form the table has to be read
against: **no axis value may be inferred from another axis's value**, and specifically **`B-EXECUTED`
is not evidence for `A-REDUCED`**. That code ran is axis B's claim; that the quantity reduces to
something more fundamental is axis A's.

⇒ **What `A-REDUCED` must cite is a DEFINING RELATION** — one that expresses this quantity **in terms
of other named quantities** — together with the loci of those quantities. ⛔ The locus of the
computation that produced the value does **not** satisfy it, and this holds whether the pointer is made
**per row** or asserted **once as a file-wide rule** ("the reduction is performed at the binding site
itself"). ⚠ A file-wide version is the worse of the two: it decides every computed row at once, without
any row's evidence ever being opened.

⭐ **The general form, because `A-REDUCED` is not the only row it reaches.** Wherever an evidence row's
locus could be satisfied by a **call** rather than the **object**, it is the object that is required:
`PHYSICS-FED` wants the **equation being solved**, not the solve (§3.3.1 already says a routine is not
an operator); `A-CALIBRATED` wants the **benchmark value**, not the fit that consumed it;
`A-INDEPENDENT-VARIABLE` wants the ledger **varying** the quantity, not an executed expression
consuming it. ⚠ That last one costs the most if it slips: it is the one axis-A value that drains the
denominator without the convention flag (§3.1), and **every** input of **every** computation is
consumed by something.

**(2) ⭐⭐ ABSENCE OF A DENIAL IS NOT EVIDENCE.** A record that a route is **un-executed**, **pending**,
**deferred**, **not yet attempted** or **owed** says nothing about whether it is **buildable**. ⛔ It
may not be cited in support of `executable within this framework` (§3.1.2), nor of **any** positive
evidence requirement in the table above.

⚠ **The general form:** a source that is **silent** on a question is not a source that **answers** it —
and under §1.1 an unanswered question falls to the conservative branch, so the silence cannot be read
in the optimistic direction either. ⭐ The tell is a citation whose content is a negation or a status
flag (`…_executed=False`, `PENDING`, `DEFERRED`, `TODO`): those are axis-B facts about what ran, and
§9.2(1) already says what they cannot decide.

**(3) ⭐ A CITED ROUTE MUST BE DIRECTION-CHECKED.** Every evidence row that accepts a **route** —
`A-REDUCIBLE-UNDERIVED`, `A-IRREDUCIBLE-STRUCTURAL` where a route is named, `A-REDUCED`'s defining
relation, and the `named-route` `should_be_basis` — requires that the relation **express the quantity IN
TERMS OF other quantities**.

> ⛔ **A relation that CONSUMES the quantity as an input does not reduce it**, and satisfies neither the
> debt nor the structural evidence requirement.

⇒ **The check, stated so it can actually be run:** **open the relation and confirm the quantity appears
on the DEFINED side, not only among the inputs.** ⚠ A route *out of* a quantity and a route *to* it are
indistinguishable from the citation alone — both look like "the route at `X:nnn` relates these symbols"
— which is why the check is on the **relation** and never on the name of the route or the note that
records it. ⭐ Where the direction fails, the citation is struck and the row is re-adjudicated on
whatever evidence remains; ⛔ it is not repaired by re-describing the same relation.

---

## 10. Denominator, aggregation, conflict

### 10.1 Denominator — a RANK, not a row count

⭐ The census's ultimate denominator is:

> **`N_physical_inputs` = dim( parameter variety, after the derived relations are imposed, modulo the
> convention group )**

**Occurrence counts and QID counts are audit receipts, not the number of independent inputs.** Two
rows can encode one degree of freedom, and a convention-laden row (§3.4) encodes none.

The both-denominators rule stands, subordinate to that:

⛔ **Never report a bare fraction.** Every *reported fraction* carries **both** denominators:

- the **occurrence count**, and
- the **adjudicated-quantity (QID) count**.

They differ, and which one is meant changes the answer. A quantity derived once and declared in six
places scores very differently under the two.

**Excluded from the denominator:** the out-of-scope rows of §7.3, listed with their reasons; and the
convention-laden occurrences of §3.4, quotiented out.

⚠ **`A-INDEPENDENT-VARIABLE` is a third case and is handled differently from both.** A coordinate is
not a parameter, so it contributes **nothing to the rank**. But it is in-universe, and it stays in the
occurrence and QID **audit receipts** under `no-tier:independent-variable` (§5.7, §5.8 output 20)
rather than being deleted from them. ⭐ The distinction is load-bearing while the rank is not yet
computable (below): the receipts are what an auditor can actually check, so a row that leaves them
leaves the census's only checkable account.

⇒ **The fraction is the audit trail. The rank is the answer.**

⚠ **Substrate note — this is not a new invention.** `notes/parameter_register.md`'s own stated purpose
is already the "**dimension of the parameter variety after the derived relations are imposed**, not
the number of symbols" (`:3-12`, the phrase at `:10-11`), and its "Relations & reductions ledger (the
edges)" (`:264`) — **96 edges, `R1`–`R97` less `R49`** — is the existing machinery for imposing those
relations. ⇒ The census **connects to** that machinery rather than replacing it.

⚠ **Flagged as a limit, not a defect:** computing a rank requires the relations to be independent, and
**independence of that edge set is not currently established.** Marked `OPEN-METHOD` (§13). The rank
is the target quantity; the row counts are what the census can deliver first.

### 10.2 Aggregation

**A QID's reduction tier = the strongest reduction achieved at any of its occurrences.**

The ordering of "strongest" is `TIER 3` (emergent) `>` `TIER 1` (irreducible), and **that is the whole
ladder**. ⚠ Within tier 1 there is **no reduction ordering** — `tier1-debt`, `tier1-structural` and
`tier1-postulate` are not more or less reduced than one another, and a disagreement among a QID's
occurrences *inside* tier 1 is a conflict (§10.3), not something the max resolves.

#### 10.2.1 ⭐⭐ Calibration is NOT a rung on that ladder

⛔ **TIER 2 is deliberately absent from the ordering above.** A calibration is not a weaker reduction;
it is a **different axis** — a tuning that was *spent*, not a derivation that was *achieved*.

> **A QID with ANY `A-CALIBRATED` occurrence is counted as calibrated — a spent tuning — REGARDLESS of
> any other occurrence's reduction. Its reduction tier is reported alongside, never instead.**

**Why.** Under a single max ordering with `TIER 3 > TIER 2`, a quantity calibrated to a benchmark in
one artifact and reduced-and-physics-fed in another maxes to tier 3 and **leaves the tuned-knob count
entirely** — the tuning disappears from the ledger's accounts because the quantity is also derivable
somewhere. That contradicts the counting rule §5.2 cites as decisive:
"**IMPOSED calibrations STAY COUNTED** — a tuning is NOT a reduction"
(`notes/stages/ledger_stage043_irreducible_count_range.md:78`). ⭐ **A benchmark that was consumed stays
consumed**, whatever else the quantity can also be shown to follow from.

⚠ A QID whose occurrences are *all* calibrated has no reduction tier and is reported in the calibrated
count only (§5.8).

#### 10.2.2 Occurrences with no tier

`no-tier:unclassified-nonfed` (§5.1), `no-tier:convention` (§5.2) and `no-tier:unadjudicated` (§5.3,
§3.4) carry no tier, so the max above does not place them. The rule:

1. ⛔ **Tier-less occurrences never enter the max.** They cannot raise or lower a QID's reduction tier.
2. If a QID has **≥ 1 tiered occurrence**, its reduction tier is the max over the tiered ones; its
   tier-less occurrences are still reported, as a per-QID count, and the QID is flagged
   **mixed-adjudication**. ⛔ A tiered occurrence does not license discarding an unadjudicated one.
3. If **all** of a QID's occurrences are tier-less, the QID has no tier and is reported in the bucket
   of its tier-less kind. ⚠ Where its tier-less occurrences disagree **in kind** — one convention, one
   unadjudicated — that is a conflict (§10.3), not a precedence question.
4. ⭐ **For the tier-1 range (§5.6):** at QID level, a QID already in the *lower* bound by rule 2 is ⛔
   **not** added again to the upper bound by its unadjudicated occurrences. A QID widens the QID-level
   span iff it has **no tiered occurrence** and **at least one** `no-tier:unadjudicated` occurrence —
   which covers both the wholly-unadjudicated QIDs of rule 3 and the mixes of rule 5. ⛔ A QID with no
   tiered occurrence does not leave the range merely because *some* of its tier-less occurrences are of
   another kind; that would let the span be narrowed by rows that established nothing, against §1. At
   occurrence level the range is computed over occurrences, where no such deduplication applies.
5. ⭐ **A QID mixing `no-tier:unadjudicated` with `no-tier:unclassified-nonfed`.** It is reported in
   **both** buckets, its tier-less occurrences are counted per-QID as under rule 2, and ⛔ it is
   **not** entered in the conflict set (§10.3).
   ⚠ **Why this mix is not a conflict, while rule 3's convention-vs-unadjudicated mix is.** A conflict
   is a disagreement between two claims **about the quantity**. `CONVENTION-LADEN = true` is such a
   claim: §3.4's transformation group and invariance are properties of the *quantity*, so a convention
   verdict at one binding site collides with any other disposition of the same QID.
   `no-tier:unclassified-nonfed` claims nothing about the quantity — it records that **this
   occurrence's** closure was determined and contained no physics — and `no-tier:unadjudicated` records
   that another occurrence could not be established at all. One says something local, the other says
   nothing; ⛔ two statements cannot disagree when one of them is silence. ⚠ Entering the mix as a
   conflict would fill a first-class output with non-findings and dilute the ones §10.3 exists to
   surface.

### 10.3 ⭐ Conflicts are a first-class output

**Where a QID's occurrences disagree, that is a CONFLICT.** Report the **conflict set with its count**
as a first-class census output, next to the tier list.

⛔ **The conflict set is computed from `is_tier` ONLY** (§6.1). It is not filtered by, ranked by, or
cross-checked against `should_be_tier` — a conflict is a disagreement between two pieces of evidence,
never a disagreement between evidence and expectation.

⛔ **Do not silently resolve conflicts by precedence.** A quantity derived in one stage and declared in
another is *exactly* the kind of thing this census exists to surface; a precedence rule would delete
the finding and leave a clean number in its place.

⭐ **Conflicts WITHIN one occurrence — the intra-occurrence case.** The rule above compares a QID's
occurrences with each other. A second shape exists and needs its own line: **one** occurrence whose
substrate carries **two incompatible claims about it** — the same value classified in two classes at
once (§13 records `decisions/14` doing this). The rule:

- it is recorded **on the row itself**, carrying **both** claims and **both** loci;
- it is reported in the **same conflict set**, marked **intra-occurrence** so the two shapes stay
  distinguishable in the count;
- ⛔ it is **never** resolved by picking one claim. §2 already forbids inheriting a substrate class at
  all, so there is nothing to pick *between*: the census adjudicates its own axis values from evidence,
  and the substrate's disagreement is recorded as a **finding about the substrate**.

⚠ An intra-occurrence conflict does **not** by itself demote any axis. Where the census's own evidence
settles the axes, the row keeps them; where it does not, §9.0's per-axis fallback applies as usual — ⛔
the substrate's confusion is not a separate demotion mechanism.

⭐ **The intra-occurrence mark carries its KIND**, because two different disagreements now land in it
and §10.3's count is only useful if they stay apart:

- **`substrate`** — the shape immediately above: the artifact or a substrate taxonomy carries two
  incompatible claims about one occurrence;
- **`census-rows`** — the shape immediately below: the census's **own rows** disagree.

#### ⭐ Intra-file consistency — where the census disagrees with itself

The shapes above are disagreements the census **finds**. This one is a disagreement the census
**makes**: two rows that cite the **same expression at the same locus** and read it incompatibly — one
treating the machinery as present, another treating the same machinery as absent.

> ⛔ **They must AGREE.** Both rows answer to one artifact, so at most one reading is right, and the
> census is the party that adjudicates — the reading is settled from the evidence **before the rows
> ship**.

⚠ **This is not in tension with the "never resolved by picking one claim" rule above, and the
difference is exactly whose claims they are.** The substrate case is two claims made by **others**,
which the census records rather than arbitrates (§2 forbids inheriting a substrate class at all, so
there is nothing to pick between). This case is two claims made by **the census**, on evidence the
census gathered, and leaving them unresolved is not neutrality — it is shipping a finding the pass
already had everything needed to settle.

⛔ **Where it genuinely cannot be resolved**, it is recorded in the same conflict set, marked
**intra-occurrence** with kind **`census-rows`**, carrying **both** rows and the **shared locus**. ⚠ It
is then a finding **about the census**, not about the substrate, and it is a defect in the pass that
produced it. ⭐ A pass that ships one has **established neither reading**, so by §1.1 both rows fall to
their conservative branch until one of them is evidenced — ⛔ the optimistic reading is not the one that
survives on the strength of having been written twice.

The §10.2 max-rule and this rule are not in tension because the max is **reported alongside** the
conflict set, never instead of it. A QID that appears in the conflict set is legible as such wherever
its tier is quoted.

---

## 11. What this schema explicitly does NOT do

Stated so that scope creep is visible later.

- **It does not rule on dimensional correctness.** ⚠ The per-quantity tallies already in the corpus —
  stage023's `24 CORRECT / 0 WRONG / 10 UNDETERMINED` (`ledger_stage023_…md:394`) and
  `27 / 0 / 7` (`:396`), and stage016's 21-of-21 CORRECT (`ledger_stage016_l2_so3_covariance.md:175`,
  with 12 of 21 declared literals in both engines at `:194`) — are **dimensional
  CORRECT/WRONG/UNDETERMINED verdicts, a different axis from provenance.**
  The census **may reuse their row universes and their correctness fields**. ⛔ It must **not** be
  expected to reproduce those distributions, and they are ⛔ **not an oracle for it**.
- **It does not establish that any value is right.** Axis A records reducibility, not truth. A
  `tier1-postulate` may be a bad postulate; a `TIER 3` value may be emergently wrong.
- **It does not re-run or re-verify any stage.**

---

## 12. Overlap with Part VII

`research/pde_ledger_v2/notes/part7_integration_atomic_split.md:67` requires of stage **046** a
calibration map over **every constant**; this census is over **occurrences** (§7.2) and is broader.
⇒ ⛔ **State it as overlap, not equivalence** — neither discharges the other.

---

## 13. Open items carried by this spec

| tag | item | where |
|---|---|---|
| `OPEN-METHOD` | the rank denominator needs `parameter_register.md`'s relation edges to be independent; **independence of that edge set is not established**, so the rank is not yet computable | §10.1 |
| `OPEN-METHOD` | the universe rule requires **one corpus-wide closure pass** to emit per-artifact slices and reachability witnesses; ⛔ membership is not decidable from a single artifact and this spec does not pretend it is | §7.1.1 |
| `PARTIAL` | the **reported-result set** exists as `notes/census/REPORTED_RESULTS.md` for stages **016, 023 and 043 only**; it is ⛔ not a partial census and ⛔ not amendable once classification begins. The set for the remaining stages is still owed, and is required *before* any of their rows is classified | §7.1.1 |
| `OPEN-SUBSTRATE` | the wrong stage016 range `:355-366` is carried by five further corpus files, including the record of the previous repair; ⛔ **not fixed here** — another workstream | §9.1 |
| `OPEN-SUBSTRATE` | stage043 note-vs-script disagreement: roll-up sub-tags vs disjoint peer categories, "≈ 152" vs exactly 152; the census cites the **script** as authoritative for counts | §2, §8.4 |
| `OPEN-SUBSTRATE` | `decisions/14` classifying `λγ` and `g_mhat` in two classes at once — now handled by §10.3's **intra-occurrence** rule: recorded on the row with both loci, reported in the conflict set, ⛔ never inherited as a class | §2, §10.3 |
| `OPEN-SUBSTRATE` | `part7_integration_atomic_split.md` internally inconsistent on stage 046 (VII-3 vs VII-4; computable vs documentary) | §12 |
| `OPEN-SUBSTRATE` | the `paused-or-pending` exclusion is **provisional by construction** — a paused approach is expected to return, so its rows re-enter the denominator when it does; ⛔ the sub-count must be re-checked on every pass, never inherited | §7.3 |
| `BUILDER-DECISION` | operational definitions of **artifact** and **binding-site**, including "one stage's `.py` and `.wl` are two artifacts" | §7.2 |
| `BUILDER-DECISION` | the **Wolfram-global** binding-site convention: one site per top-level `Set`/`SetDelayed`, `Module`/`Block` locals excluded, pattern definitions bind a rule not a value, artifact-set system globals included | §7.2 |
| `BUILDER-DECISION` | an `UNADJUDICATED` convention claim buys no exclusion — the row keeps its axis-A tier and is reported in `convention-unadjudicated` | §3.4 |
| `BUILDER-DECISION` | `C-SELF` is the tie-break where a leaf reads as both `C-SELF` and `C-MATH`; it is the conservative direction and it moves no headline | §3.3 |
| `BUILDER-DECISION` | the two `C-PEER` populations (artifact-carried vs stage-note-carried loci) are reported **separately and never summed**, and a chain records its weakest hop | §3.3, §5.8 |
| `BUILDER-DECISION` | ⭐⭐ the asymmetry binds **assertions**, not only defaults: the optimistic branch carries the evidence burden and an unmet burden falls to the conservative branch. Every later addition to §9 is written to this shape — what the optimistic branch cites, and the branch it falls to | §1.1, §9 |
| `BUILDER-DECISION` | ⭐⭐ **no axis value is inferred from another axis's value**; `A-REDUCED` requires a **defining relation**, not the locus of the computation. Swept across §9: `PHYSICS-FED`, `A-CALIBRATED` and `A-INDEPENDENT-VARIABLE` were rewritten to name the **object** rather than the **call** | §3, §9, §9.2 |
| `BUILDER-DECISION` | `executable within this framework` requires every consumed quantity to have a **functional form present in the ledger, cited by locus**; a name, a dimension rule, a free symbol and a literal scalar are not forms. Missing form ⇒ `A-IRREDUCIBLE-STRUCTURAL` with the form named, else `A-UNADJUDICATED` | §3.1.2, §9 |
| `BUILDER-DECISION` | **absence of a denial is not evidence** — un-run / pending / deferred says nothing about buildable; stated generally over §9's whole table | §9.2, §3.1.2 |
| `BUILDER-DECISION` | a cited route is **direction-checked**: the quantity must appear on the **defined** side, not only among the inputs | §9.2, §9 |
| `BUILDER-DECISION` | route 2's **third bound** — no restatement of a route-1 freedom; one source literal is **one** occurrence unless each part's independence is evidenced. ⚠ The **fold** costs its own evidence (it shrinks tier 1), so a fold names the route-1 quantity at its locus | §7.1, §7.2, §9 |
| `BUILDER-DECISION` | a *"would land in tier N"* claim about a row the pass did **not classify** is barred; only the count of unreached rows and the reason may be stated | §5.6, §7.3.1 |
| `BUILDER-DECISION` | the conflict set's **intra-occurrence** mark carries a **kind**: `substrate` (the artifact disagrees with itself, recorded not arbitrated) vs `census-rows` (the census's own rows disagree about one expression at one locus — must be resolved, and if unresolvable both rows fall conservative) | §10.3 |
| `NOT-YET-CREATED` | `QID_REGISTRY.md` — required by §8.3, not authored by this spec | §8.3 |

---

## 14. ⭐⭐ FOR THE READER OF THE RESULT — the reassessment trigger

> ⛔⛔ **ADDRESSED TO WHOEVER READS THE FINISHED CENSUS. NOT to whoever takes it.**
>
> ⛔ **If you are classifying rows, STOP HERE.** Nothing in this section is an input to §3–§9, and
> §3–§9 deliberately contain no statement of what the tier-1 answer is expected to be. Reading an
> expected magnitude before assigning rows contaminates the measurement — a classifier who knows the
> target is no longer measuring, and cannot tell from the inside that they have stopped.

**The standing expectation.** The model is expected to take **few physical inputs** — the brane/bulk
defining properties plus a small number of extras.

**The trigger, and it fires in two directions:**

- ⇒ **If the tier-1 RANGE (§5.6) comes in substantially above that**, the correct response is to
  reassess what those rows are and why there are so many — **not** to report the number and move on.
- ⇒ ⭐ **If the UNADJUDICATED SPAN is wide** — that is, the range's upper bound sits far above its
  lower bound — the correct response is the same, and it is triggered independently. A wide span says
  the census could not establish a large part of its own denominator. ⛔ A narrow, comfortable lower
  bound underneath a wide span is **not** a reassuring result; it is an unread one.

A large tier 1 is a signal that something is missing or miscounted. A wide span is a signal that the
evidence was not there to look. In both cases the `should_be` column (§6) and the conflict set (§10.3)
are where the reassessment starts.

⛔ This is a **reporting obligation**. ⛔ It is **not** a licence to adjust rows until the count looks
right.

⭐ **The honest outcome to hold open.** §6.2 states the two readings a large tier 1 could carry — a
reduction debt, or a medium that needs more specification than the picture suggests — and ⛔ **does not
say which.** Neither does this section. That gap is the result.
