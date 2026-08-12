# Script rebuild — state and what's next

Read `CLAUDE.md` first (how we work), then this (where we are).
Process detail: `docs/development_pipeline.md`. Defects: `DEFECT_REGISTER.md`.

Last updated 2026-08-09. ⛔ The registry design is ABANDONED — read `S9_REWRITE_PLAN.md` before anything else.

---

## The scope — three steps need rebuilding, and the debt is exactly these

The measured defect was in **three independently-built steps**: engines emitting physics conclusions as
typed sentences with no CAS object behind them, missed by eight review legs.

| step | state |
|---|---|
| S9 | rebuilt |
| S10 | ▶ **CHAIN + COMPARATOR + RECORD + PAPER CARD ALL DONE** — export `e644876c`+`c84263ed`, comparator `82443c95`, **record `e167b07f`** (664 → 828 lines, 4 builds / **14 legs**), **card rebuilt `2998029f`+ fix round 1** (299 → 477 lines, 2 builds / **8 legs**). ⛔ **No number in either was ever found wrong.** ⭐⭐ **NEXT: D12 → registers.** `.wl` untouched |
| S11 `stray_longitudinal` | ▶ **IN PROGRESS** — spec done `f49a1684`; engines not yet rebuilt |
| S11b-A `interface_response` | built under the old pattern → **rebuild** |
| S11b-B `interface_assembly` | built under the old pattern → **rebuild** |
| S11b-C non-uniform coupling | **never built** — the MacCullagh differentiator |

S12 onward does not exist yet, so every remaining step in every remaining sector gets built under the new
pattern the first time. This tax is paid once.

### ⛔⛔ THE PATTERN CHANGED 2026-08-08 — ⭐ read `S9_REWRITE_PLAN.md`

⛔ **The old pattern told you to write a `reduction/checks_S<n>.yaml` and register quantities.** ⛔ That
design is abandoned. ⛔ Do not build it for any step.

⭐ **The pattern now:**

1. Both engines emit **CAS objects**. ⛔ No typed conclusions; a residual asserted zero always prints `0`
   and carries no information — ⭐ emit both operands and the residual, then guard.
2. ⭐ **Both engines emit the STANDARD NAME of each object** (below). ⇒ comparison is a lookup.
3. ⭐ The **SymPy** engine writes a flat `LEDGER` the next step imports; ⭐ the **Wolfram** engine imports
   nothing and re-derives independently. ⇒ consistency by dataflow, contention py-vs-wl within a step.
4. ⭐ Every declared symbol carries a **class tag and an English description** —
   `KNOB · STRUCTURAL · COORDINATE · CONTROL · PREMISE · DERIVED`.
   ⛔⛔ **NEVER annotate what the computation will produce** — no value, dimension, count or sign.
5. ⭐ The step record **cites what it measured**, ⛔ not what the author concluded.
6. ⭐ Requirements-register entries fall out of the same read ⇒ `SUBSTRATE_REQUIREMENTS.md`.

---

## ⭐⭐ S9 — CLOSED `5d6b56d3`, RE-KEYED `83eace28`, STRUCTURALS EXPORTED `b23c5334`. ⭐ NEXT IS S10

### ⭐ The ledger is **145** entries — KNOB 2 · STRUCTURAL 4 · PREMISE 7 · DERIVED 132

⭐ **`D`, `length_dimension`, `time_dimension`, `mass_dimension` now cross directly** ⇒ a consumer **binds**
them instead of declaring look-alikes. ⚠ `Symbol("D")` ≠ `Symbol("D", integer=True, positive=True)`; their
difference **prints as something that looks like zero and is not** ⇒
[[project-ledger-symbol-identity-is-the-point]].

⛔⛔ **FOUR MEASURED LIMITS — ⭐ they bound what S9 may claim, and ⛔ none is fixable inside the export:**

1. ⛔⛔ **The assumptions the export exists to carry are INVISIBLE in the record.** The assumption set
   **hand-types** `Q.integer(D)` beside the declaration instead of reading it, and both spellings print
   identically ⇒ dropping `integer=True` at the declaration leaves stdout **byte-identical** while **10**
   ledger records change. ⭐ **S10 binding `D` is the test** — ⛔ not another in-run check.
2. ⛔ **The export boundary is a hand-typed comment.** Restoring one pre-round tag ships the
   **ablation-package table as a structural input at exit 0**; the round-trip rejects it only when the
   value is not a SymPy object.
3. ⚠ **The unit-basis convention is unpoliced inside S9** — every dimension is built *from* the markers, so
   permuting L↔T moves **no** residual. ⭐ The `.wl` declares the same basis independently and constrains
   the derived vectors, ⛔ but **neither engine emits the markers**, so the four new records have **no
   cross-engine row**.
4. ⚠ `field_dimension` is an **alias** of `length_dimension` and both are now keys. ⚠ `sp.eye(3)` is still
   typed at `:89`; `:95`'s read resolves a **typed triple**.

### ⛔⛔ THREE IN-RUN CHECKS ON THE EXPORT HAVE NOW BEEN BUILT AND DELETED. ⭐ STOP BUILDING THEM.

⚠ Measured across three rounds: **every** guard added inside the export writer compared two operands
descending from one source, and each was proven unable to fail — one of them by **reverting the change it
policed** and still exiting 0. ⇒ ⭐⭐ **The export writer cannot verify itself: the classification is its own
input.** ⭐ What establishes the export is a **diff against the committed baseline**, from outside the run
⇒ [[feedback-a-check-cannot-audit-its-own-input]].

⭐ **CONTROL now carries two meanings** — an ablation coefficient **and** engine machinery (user decision,
2026-08-08) ⇒ ⚠ the knob inventory mixes them, and **the tag decides what the next step imports.**

### ⭐⭐ THE D-KEY — decided by the user 2026-08-08, and it is a STANDING convention

⭐ **A ledger key records the component count of the CONSTRUCTION that produced the object.** Built from
the fixed-component action ⇒ carries the count (`transverse_speed_squared_d3`). Built from unit algebra or
from the dimension solve's symbolic output ⇒ carries none (`dim_solution`, `field_dimension`).
⭐ Now **130 suffixed · 11 unsuffixed**; 126 keys moved, **0** payload fields did.

⛔⛔ **The key does NOT record whether the value depends on `D`, and ⛔ you cannot make it.** ⚠ Measured:
`transverse_multiplicity` is the bare value `2` and is `D−1`; `lagrangian` has no `D` and is D=3 by its
`t,x,y,z` coordinates. ⇒ **free-symbol inspection is wrong in BOTH directions** (a leg measured 8 of 141
misclassified). Deciding D-dependence needs a run at a **second component count**, which S9 does not have.
⭐ Four equal D-slots at S10 are the measurement instead.

⚠ **Two limits the step record must carry:**
1. ⛔ **Nothing inside the artifact distinguishes the live count readout from a literal** — the engine runs
   at one count, so `f"_d{n}"` and `"_d3"` produce byte-identical output. ⭐ Only a **cross-count ablation**
   separates them, and it lives in the review, ⛔ not in the artifact.
2. ⛔ **The `.out` record name for the four D-specialisations is still a hand-typed literal** (`DIM_SOLUTION_D3`
   …), **36 lines at HEAD and 36 now** ⇒ ⚠ **pre-existing, ⛔ not introduced by the re-key.** Under mutation
   the ledger key follows the computation and the **emitted tag name does not** ⇒ this is a **D11 naming**
   question ⇒ [[feedback-name-binding-is-unpoliced]].

⭐ `EXPORT_D_PARTITION` is **printed once and guards nothing** — deliberately ⇒
[[feedback-a-check-cannot-audit-its-own-input]].

⭐ Cleanup `67dd3ce2` + `fb29bba2`; build `753ae7b1`; **seven** fix rounds through `890a359d`.
⭐ Each round: Codex built, **two independent legs** reviewed (fresh agent + Grok), commit after both.
⭐ Deliverable: `scripts/S9_exports.py` — flat **141** entries, `{display, value, dim?, class, step}`,
**KNOB 2 · PREMISE 7 · DERIVED 132**, 3 with a computed `dim`.

### ⛔⛔ READ THIS BEFORE OPENING ANOTHER FIX ROUND ON ANYTHING

⚠ **Rounds 1–5 repaired the ARTIFACT. Rounds 6–7 repaired MY OWN EDITS.** Round 6 was opened in answer to
*"any smaller changes?"*, bundled a **schema change**, and introduced two defects; round 7 repaired those
and introduced one; the close-out was **one deleted line**. ⇒ ⭐ **after a round that removes code and leaves
output byte-identical, STOP** ⇒ [[feedback-fix-rounds-converging-on-my-own-edits]].
⚠ Across all seven rounds the **derivation never moved** and every leg reproduced every value.

## ⭐⭐ S10 — THE CHAIN IS REAL. ⚠ Three open defects; ⛔ do NOT treat S10 as closed.

⭐⭐ **`scripts/S10_exports.py`: 617 entries** — KNOB 2 · STRUCTURAL 4 · COORDINATE 27 · CONTROL 9 ·
PREMISE 78 · DERIVED 497; **162 carried from S9, 452 added by S10** (S9's own ledger is **165**).
⚠ **Re-measured 2026-08-09 from the modules; the previous figures here (574/145/142) were STALE** — a
review leg caught them ⇒ ⭐ **numbers in this file are a convenience, ⛔ the artifact is the authority.**
⭐ **3 keys overwritten and AGREEING · ZERO drift in anything untouched.**

⛔⛔ **AND THE THREE OVERWRITES ARE TWO INDEPENDENT VALUES.** `coefficient_dimension_difference` is
built as `stiffness − inertia` at `S10:1790-1795` ⇒ ⛔ it **cannot** disagree once the other two agree.
⛔⛔ **All three are DIMENSION VECTORS ⇒ the guard is blind to any DIMENSIONLESS mutation BY
CONSTRUCTION** — measured: a ×7 rescale moves the speed to `7μ_R/ρ_br` and S10 **exits 0 and publishes**,
while a FORM change makes it refuse. ⚠ Reproduced independently by a review leg from scratch.

⭐⭐ **What this bought — the first cross-step corroboration the ledger has ever had.** S10 re-derives
`inertia_coefficient_dimension`, `stiffness_coefficient_dimension` and `coefficient_dimension_difference`
from its own general-`D` action and lands on S9's keys with **matching values**: two independently built
engines, at different component counts, writing the same slot.

⭐⭐ **The physics is INDEPENDENTLY CONFIRMED.** A review leg wrote its own derivation from
`S11_SHARED_PHYSICS`-style first principles **before opening the artifact** and reproduced **26 of 26**
committed tags — root count, ordering, stacked ranks, scale exponent, sign-flip control, Q7.
⇒ **every finding below is the export chain, ⛔ not the derivation.**

⭐⭐ **The key-naming problem is SETTLED BY MEASUREMENT.** Under a FORM control that moved the root count
`2 → 3` and produced a stratum where there was none, the **exported key set came back IDENTICAL** (empty
symmetric difference both ways). Runtime root/stratum/coefficient indices are payload **inside** fixed
`indexed_derivations` collections ⇒ ⛔ no key is minted by the solver's answer.

⭐ **Assumption inheritance did what the user predicted.** Binding S9's `positive`/`integer` objects moved
exactly two control readouts — `XFORM_ANISO_D{3,4}_Q8_STRATUM1_ROOT2_Q3_SIGN`:
`undecided_under_joint_assumptions → 1`. ⚠ Not a changed answer — a **newly decidable** one. A leg derived
`+1` independently.

### ⭐⭐ ROUNDS 2–4 DONE — committed `e644876c`. ⚠ Round 5 is a REMOVAL round, in flight

⭐ **Eight legs across four rounds. ⛔ No computed physics value moved in any of them** — every transcript
change is a rename, a container retype, or export bookkeeping.

⭐⭐ **What the chain now has, ⭐ each DEMONSTRATED ABLE TO FAIL:**
- ⭐ **Every symbol in an exported value is a ledger entry** — S9 `3/23 → 23/23`, S10 `3/158 → 64/64`.
- ⭐⭐ **The S9→S10 overwrite residual is a GENUINE cross-step measurement.** A FORM change
  (`curl·curl → ∇²u·∇²u`) makes S9 publish `transverse_speed_squared = mu_R*q/rho_br` — ⛔ no longer a
  speed — and S10, re-deriving from its own solve, reports `value_residual = 1` and **refuses to publish**
  (`S10:2110`). ⭐ A **coefficient** rescale moves nothing. ⚠ Reproduced in two independent sandboxes.
- ⭐ `corroborated_steps` records the S9/S10 agreement **in the ledger**; `dimension_key` makes units
  discoverable; `LEDGER` is a `mappingproxy` of `mappingproxy` with **0** mutable values; both engines are
  **byte-deterministic across four independent rebuilds**; `BUILD_INPUT_DIGESTS` makes a stale module
  detectable.
- ⭐ **The symbol-identity scan is the ONE in-run check that works** — it walks the **written module**, and
  it rejects its weakest change.

### ⛔⛔⛔ THE DEFECT OF ROUNDS 2–4 WAS MINE — ⭐ I commissioned FIVE in-run checks this file already forbade

⚠ **This document said, in bold, "three in-run checks have been built and deleted — STOP BUILDING THEM."**
⭐ I sanctioned one (the identity scan) and then wrote three decision lists commissioning **five more inside
the export writer.** ⇒ ⛔ each round's own review prompt warned about the defect while the round built it.

⛔ **`value_kind`** is written by `export_value_kind(value)` and checked against `export_value_kind(value)` —
**the same function**; falsifying it to one constant leaves every residual 0 and publishes.
⛔ **`BUILD_INPUT_DIGEST_RESIDUAL`** compares a `str→str` dict with itself across `repr`→`eval`.
⛔ **The assumption-channel residual scores 0 when the ENTIRE channel is deleted.**

⇒ ⭐⭐ **The tell to use: does this residual live inside the thing that produced its operands?** If yes it is
decoration, whatever its shape ⇒ [[feedback-a-check-cannot-audit-its-own-input]],
[[feedback-control-outside-the-thing-it-polices]].

### ⛔ MEASURED LIMITS — ⭐ record them, ⛔ do NOT build checks for them

⚠ Every repair proposed for these is **either a denylist or another self-comparison.**

1. ⛔ **`value_kind` marks ONE carrier of five.** A sentence reaches the ledger tagged `COMPUTED_OBJECT` as
   `Symbol('…')` (⚠ which **mints a ledger key equal to the sentence**), `Function('…')(t)`,
   `Tuple(Str('…'))`, `Eq(Str,Str)`, or in **`display`** — a raw `str` on every record, in no residual.
2. ⛔ **The authored-word count does not bound the population**: emitted **24**, records carrying authored
   text **36**, distinct authored strings **104**.
3. ⛔ **The assumption channel is one-way and half-covered** — `11/23` S9 symbols, `14/46` S10, `MAIN` only.
   ⚠ My decision list stated both denominators **one too low**; the builder measured them. Conclusion unchanged.
   ⚠ Also measured: **S9's `value_kind` field takes ONE value across its whole ledger** ⇒ upstream it
   distinguishes nothing.
   ⭐ **Revisit AFTER the comparator**: the `.wl` declares assumptions independently ⇒ an **outside oracle**,
   which is the only kind that has worked here.
4. ⛔ **A parse error leaves a stale module** — Python compiles before executing, so ⛔ no module-scope
   statement can prevent it. ⭐ Detectable via `BUILD_INPUT_DIGESTS`; ⛔ nothing checks it.
5. ⛔ **Every guard is an `assert`** ⇒ `PYTHONOPTIMIZE=1` publishes anyway.
6. ⛔ **`BUILD_INPUT_DIGESTS` omits the CAS** — no sympy/Python version.
7. ⛔ **`symbol_binding_residual` cannot see one quantity under two names** — it counts variants per name.
   ⚠ The ledger is clean of splits **by inspection, ⛔ not by instrument**; **S11 adding names is unpoliced**
   ⇒ [[feedback-name-binding-is-unpoliced]].

⚠ **Recorded, ⛔ not a defect:** the relational round-trip repair is **correct but inert** — `MAIN` produces
no allowed stratum at any count, so nothing on the export path exercises it.

### ⚠ WHAT S10 STILL OWES

| owed | note |
|---|---|
| ~~round 5~~ | ⭐ **DONE `c84263ed`** — removal round, author changed per rule 15. **Both legs CLEAR.** ⭐ The cross-step residual was re-established on the committed bytes by **three parties using two different FORM changes**, with a ×7 coefficient control that correctly moves nothing |
| ~~the comparator (F-2)~~ | ⭐⭐ **DONE `82443c95`** — 4 rounds, 8 legs. ⭐ See the F-2 section below for what it measures, the **deflations the step record must carry**, and the injectivity check D12 must run |
| ~~the S10 STEP RECORD~~ | ⭐⭐ **DONE** — 664 → 828 lines. **4 builds, 14 legs**; every number and locus verified by **three** independent passes, ⛔ none wrong. ⭐ The trigger fired on all four decision lists and **blocked every one** ⇒ see below |
| **the D12 naming pass** | ⭐⭐ **DRIVEN BY THE COMPARATOR'S OUTPUT, ⛔ not by hand.** User decision 2026-08-09: **adjudicate BOTH engines**, ⛔ not wl-only — the py tag suffix **IS the ledger key**, so py-as-default-authority is the unpoliced binding the light-cone re-point exploited |
| `.wl` rename | ⭐ engine EXISTS (1807 lines, 660 KB output) and is correct — ⛔ **do NOT rebuild it**. Emitted **name strings only** |
| **step record** `steps/S10_two_transverse_photons.md` | ⛔ **STALE — 664 lines, 23 citations to deleted machinery** |
| ~~paper card~~ | ⭐⭐ **DONE — 299 → 477 lines, 2 builds / 8 legs.** ⛔ The old *"2 stale refs"* estimate was wrong by an order of magnitude: it carried **six** defect classes, including an entire verbatim block quoting the **deleted harness** (⛔ no live artifact emits `CROSS_ENGINE_COVERAGE` / `CONTROL_RESPONSE` / `TAG_PARITY`) and a claim the record **contradicts** (*"nullspace bases are never compared for span"* vs **26 PASS / 0 FAIL**). ⭐ See the card section below for what the round measured |
| **D12 naming pass** | 2 verified pairs; ⛔ **injectivity across the worklist FIRST**. Touches **both** engines ⇒ physics-bearing ⇒ full gate |
| **defect register** | ⭐ `C16` added `e167b07f`+ — `Q8a`'s stratum enumeration. ⚠ Remaining: sweep the S10 record's 15 dispositions for anything not yet carrying a row |
| requirements register | ⭐ already populated — **11 entries, 7 naming S10**. ⚠ Check only that the record's new limits are reflected |

### ⭐⭐⭐ WHAT THE S10 RECORD ROUNDS MEASURED — ⚠ physics, ⛔ not bookkeeping

⭐⭐ **`D − 1` rests on TWO independent structural premises, and the sector's controls only ever probed ONE.**

1. ⭐ **The stiffness functional** — and it buys **less** than the sector's slogan claimed. `FULLGRAD` (full
   gradient) returns the **same root expression AND the same transverse count** as curl-only at `D = 3`
   and `D = 4`. ⇒ ⛔ *"curl-only gives the two transverse modes"* is **refuted by its own control**;
   ⭐ what curl-only determines is that the **longitudinal does not propagate**.
2. ⭐⭐ **ISOTROPIC INERTIA — and nothing in the headline mentioned it.** Breaking it on one axis
   (`ANISO`) splits the degenerate pair into two speeds and one branch becomes **generically OBLIQUE**:

   | | per-root `N2` (**modes**) | per-root `N3` (**exactly transverse**) | propagating totals |
   |---|---|---|---|
   | `MAIN` D3 | 1, **2** | 0, **2** | 2 / 2 |
   | `ANISO` D3 generic | 1, **1, 1** | 0, **1**, 0 | **2 / 1** |
   | `ANISO` D3 on stratum | 1, **2** | 0, **2** | 2 / 2 |

   ⛔⛔ **THE PROPAGATING COUNT NEVER MOVES.** ⛔ No mode is removed; ⭐ one is **tilted**.
   ⇒ ⭐⭐ **isotropic inertia is what makes the transverse modes EXACTLY transverse, ⛔ not what makes them
   exist.** ⚠ Both engines agree on every cell; **20 generic `N2`/`N3` comparator rows, 20 PASS.**

⇒ ⭐ **Scope the ablation maxim:** *"only a form change moves physics"* is **overbroad** — `XCOEF_SCALE`
moves the **root value**, `SIGNFLIP` turns two waves into two **growing modes**. ⭐ What a coefficient
rescale cannot move is a **COUNT** or a **MODE STRUCTURE** ⇒ [[feedback-per-tooth-ablation]].

### ⛔⛔ TWO OPEN FINDINGS FROM THE FINAL LEG — ⭐ recorded, ⛔ NOT folded into the record

1. ⛔⛔ **`Q8a`'s STRATUM ENUMERATION IS STRUCTURALLY BLIND TO THE LOCUS WHERE TRANSVERSALITY CHANGES.**
   ⚠ It searches **minors of `M_r`** and **root-coincidence loci** (`S10_SHARED_PHYSICS.md:390-396`);
   ⛔ **neither looks at the STACKED matrix `[M_r; kᵀ]`**, which is what governs `nu_T`. ⇒ a locus where the
   **transverse count alone** moves **cannot be found by construction.**
   ⭐ **One exists and was exhibited** by a leg that re-derived the `ANISO` matrix from the spec: at
   `k₁ = 0` (`k` ⊥ the anisotropy axis) the exactly-transverse total returns to `D − 1`, with the **root
   count unchanged** and **every `rank(M_r)` unchanged**. ⛔ Neither engine enumerates it.
   ⚠ **This falsifies nothing in the record** (which says *"generically"* and disclaims completeness);
   ⭐ it is a **SPEC-LEVEL** defect — the engines implement `Q8a` faithfully. ⇒ **red-team / S11 spec item.**
   ⭐ `MAIN` is **not** exposed: a 748-wavevector sweep found **0** deviations.
2. ⚠ **A correct `MAIN` limit was DOWNGRADED relative to HEAD** — *"there is nowhere in the allowed region
   for `D−1` to fail"* is still **decided by both engines** (`Q8_ALLOWED_STRATA` empty at every `D`, by two
   independent routes: PY `locus_conflicts_with_positive_wavevector_norm`, WL `decidedEmpty`).
   ⭐ The new record keeps it only as a comparator-table row. ⛔ **Under-claiming, not overstatement** — so
   it was not worth a fourth round; ⚠ fix it if the record is opened again.

### ⭐⭐ THE PAPER CARD ROUND — ⚠ what it measured, and ⛔ the queue it leaves

⭐ **`2998029f` (rebuild) + fix round 1. 2 builds, 8 legs, 4 gate cycles.** ⛔ Both legs re-parsed the
comparator and re-read every cited engine locus; ⛔ **no number was found wrong in either build.** ⭐ Every
defect was a **claim, a strength, or an evidence trail.**

⭐⭐⭐ **THE PHYSICS DEFECT, and it is the whole reason the round was worth running.** The rebuilt card's
headline listed **nine** conditions — and **`XFORM_ANISO` satisfied every one of them** (positive
coefficients, cosine ansatz, unstrained rest, no dissipation, linear, nonzero `k`, generic, **and
`S_curl` stiffness**, since `L_XFORM_ANISO` keeps `S_curl` exactly) **while measuring `D − 2`**
exactly-transverse. ⇒ ⛔⛔ **the second premise S10's record was reopened to establish was missing from the
one sentence a reader quotes.** ⚠ Only the words *"MAIN case"* excluded it, ⭐ and a reader reads those as
naming the sweep, ⛔ not as a physics premise. ⚠ **The builder's own claim ledger missed it too.**
⭐ **Closed at the headline sentence itself**, and re-verified case-by-case: ⛔ no measured package now
satisfies every listed condition and returns a different count.

⭐ **Also removed:** the analytic passage stated at **unrestricted general `D`** where the record licenses
*"for the cases measured"* (⚠ and `S10_SHARED_PHYSICS.md:44-47` says *"typing its reduced form is the
defect this whole rebuild removes"* — ⛔ the card was reintroducing it); a load-bearing citation pointing
at `V3_STEP_PLAN.md:354-380`, which **asserted what the record refutes**; two of the export guard's five
hard limits; and the absent prior-art attribution.

⛔⛔ **AND THE GATE CAUGHT ME COMMISSIONING DAMAGE, ⛔ not just missing things.** ⚠ One fix item would have
required the dimensional paragraph to carry the engines' action-term vacuity markers. ⭐ **That was not a
defect** — the card never quotes those booleans, and the record scopes the vacuity to *"the solved
action-term homogeneity booleans only"*, expressly **not** downgrading the dimension solution.
⇒ [[feedback-decision-list-length-is-the-defect-rate]]. ⚠ **Two more of my items were narrowed** — one
would have attributed the supplied split (the record attributes only the out-of-plane **object**), one
named one citation site when there were two.

### ⚠ BANKED FROM THE CARD'S FINAL LEGS — ⛔ recorded, ⛔ NOT fixed

⚠ **One leg reported these; the other returned nothing surviving its filter. ⛔ Neither changes what may be
claimed** ⇒ stopped per rule 10, ⛔ not iterated to green.

1. ⚠ **`"non-vacuous"` is a term of art the card inherits without its defining paragraph.** ⭐ Both the
   card's *"This is a homogeneity result…"* (record `:390`) and *"…other non-vacuous dimensional objects
   retain their stated evidential status"* (record `:511`) are **verbatim from the record**, so ⛔ there is
   no over-claim — ⚠ but a **card-only reader cannot tell what partition the qualifier excludes.**
   ⛔ Deleting the qualifier would be worse: it would assert that **all** dimensional objects retain status.
2. ⚠ Record item 10 (**whole-layer dimensional agreement open**) and the three `−3` Q6 unknown-count
   residuals have **no card counterpart**. ⭐ The card makes no whole-layer claim, so this is omission,
   ⛔ not over-claim.

### ⛔⛔ THE RECORD-REOPENING QUEUE — ⭐ THREE items, ⛔ do NOT open a round for any one alone

⚠ **`steps/S10_two_transverse_photons.md` is closed at `e167b07f`.** ⭐ If it is **ever** reopened, fix all
three in one pass:

1. ⚠ **A correct `MAIN` limit was DOWNGRADED** — *"nowhere in the allowed region for `D−1` to fail"* is
   still decided by both engines. ⭐ Under-claiming, ⛔ not overstatement.
2. ⛔⛔ **`:127-128` cites the cosine ansatz to `:13-28, :30-47, :82-107`. ⚠ The ansatz is at `:49-56` — in
   NONE of them.** ⭐ Found 2026-08-09 by a card leg; ⛔ the card fixed **its** copy at both sites.
3. ⛔ **`:386-389` cites the dimension DIFFERENCE `(2,−2,0)` to loci that do not contain it** (it is at
   `:4215`).

⇒ ⭐⭐ **2 and 3 falsify a sentence this handoff carried:** *"every number AND every cited locus verified by
three independent passes, none wrong."* ⚠ **True of the NUMBERS. ⛔ FALSE of the loci** — three passes
checked values and did not open every range ⇒ [[feedback-a-check-cannot-audit-its-own-input]].

### ⚠ S9's CARD HAS THE SAME DISEASE — ⛔ out of scope, ⛔ do not fold in

⚠ `paper/steps/S9_light_requires_shear.tex` quotes the **deleted harness**: *"312 comparisons"* / *"329"*,
*"3 per engine"*, *"exit 2 on the unparsed `Piecewise` tag"*, *"12 cross-engine quantities agreeing"*.
⛔ **No live artifact produces any of them**, and the card predates S9's own rebuild.
⚠ `paper/steps/S11_stray_longitudinal.tex` carries **4** references to the deleted registry; ⭐ that one
travels with S11's rebuild.

### ⛔⛔ THE "GRANULARITY DECISION" DOES NOT EXIST — ⚠ it was inherited from the S11 survey

⛔ **"py emits per-root stdout tags, wl emits lists" is FALSE for S10.** ⭐ Measured: both engines emit
**per-package, per-`D`, per-root** tags, and **562 fully-qualified names already match byte-for-byte** —
including `Q3_DETERMINANT`, `Q2_MATRIX_A`/`_B`, `ROOT#_N1_MATRIX`, `Q3_ROOT_COUNT`. ⇒ the difference is
**vocabulary**, which is exactly what D12 is for.

### ⭐ (superseded) the orchestrator prototype

⛔ **Its numbers are SUPERSEDED by the committed comparator** — see the F-2 section below. ⚠ The prototype
read `348` agreements where the committed instrument reads `388`; ⭐ the difference is categorisation
(non-canonical and unparsed rows), ⛔ not a disagreement about any value.
⭐ **What the prototype was for, and it worked:** it was withheld from the builder so a differing tally
would be a finding. ⭐ The join agreed exactly — `562` shared names on both — which is the corroboration
that mattered.

---

## ⭐⭐ F-2 IS BUILT AND COMMITTED `82443c95` — the cross-engine join EXISTS

⭐ `scripts/S10_cross_engine_comparator.py` + committed `scripts/out/S10_cross_engine_comparator.out`.
⭐ **4 rounds, 8 independent legs.** ⛔ No config, ⛔ no YAML, ⛔ no pair table — it reads exactly two `.out`
files and joins on the standard name.

⭐ **Why it was needed, and it is not hypothetical:** a leg re-pointed ONE standard name so the light cone
became an `ω²` instead of a speed squared, and **every check in the repository passed.**

### ⭐⭐ WHAT IT MEASURES — ⚠ use THESE numbers, ⛔ not the headline ones

⭐ `4233` py names · `2983` wl names · `562` shared. ⭐⭐ **Every shared name carries exactly ONE verdict** —
⛔ nothing falls out of the join uncounted, which is the property that matters most, because a row the
comparator could not parse **reads exactly like a row that agreed.**

⭐⭐ **`Q3_DETERMINANT` AGREES at `D = 2,3,4,5` AND under all four FORM controls** — py **expanded**, wl
**factored**, so `simplify` does real work. ⭐ Also `Q2_MATRIX_A`/`_B`, `ROOT_ORDERING`, `Q3_ROOT_COUNT`,
`N1_MATRIX`, `N2_NULLITY`, `N3_STACKED_MATRIX`.
⭐ **Three independent implementations agree row for row**; ⛔ no leg found a row where the artifact reports
agreement and an independent route disputes it.

⛔⛔ **THE DEFLATIONS — the step record MUST carry these:**
- ⭐ Of **388** agreements: **215 bare integers**, **14 empty containers**, only **159** symbolic/structured.
- ⭐ Of **164** disagreements: 23 naming-only · 13 representational · 128 content — and of those 128, only
  **SIX** carry a genuine numeric/algebraic residual, **five of the six in CONTROL packages, ⛔ not MAIN.**
  ⭐ The one `MAIN` row is a sign one engine evaluates under a premise **both engines declare**, and it is
  **zero under that premise**.
- ⚠ **109 shape/type mismatches** ⇒ only **~390 of 562** shared names are value-comparable at all.
- ⛔⛔ **COMMON-MODE IS UNTOUCHED** — the action is authored identically into both engines from the shared
  spec. ⭐ What is corroborated is the route **from** the action **to** the spectrum.

### ⛔ OPEN ON THE COMPARATOR — ⭐ recorded, ⛔ deliberately NOT fixed

1. ⛔⛔ **A repointed symbol still classifies NAMING_ONLY when the substituted spelling is emitted by only
   ONE engine.** ⚠ The transcripts DO refute it — as a **non-injective worklist**: two distinct py symbols
   claiming one wl name, both marked unrefuted. ⇒ ⭐⭐ **THE D12 PASS MUST CHECK INJECTIVITY ACROSS THE
   WORKLIST BEFORE CONSUMING ANY PAIR.** ⚠ The row still FAILS its guard; ⛔ nothing is silently accepted.
2. ⛔ The constant/shadow rule is applied to **one** engine's payload parser; the other fails **loudly**
   (not-computed) rather than agreeing falsely.
3. ⛔ Two sub-population selectors test the emitted **name suffix**, ⛔ not the residual. ⭐ Clean on this input.
4. ⛔ The process guard **already fails at baseline** ⇒ the exit code carries **no** ablation signal; ⭐ only
   per-row residuals do.

### ⭐⭐ WHAT THE COMPARATOR CHANGED ABOUT THE PLAN

⭐⭐ **It SCOPES D12 rather than depending on it.** ⛔ I was about to adjudicate ~778 names by hand; the join
reports which differences actually **block** a comparison — **two** verified symbol pairs
(`D ↔ braneDimension`, `s ↔ coefficientScale`) — and it **caught a third that was MIS-TYPED**:
`M_B ↔ quadraticFormRoute` is a **route-tag** divergence, and consuming it as a rename would have bound two
different descriptions ⇒ [[feedback-name-binding-is-unpoliced]].

⭐ **It is also the OUTSIDE ORACLE for the assumption channel.** The `Q3_SIGN` rows are exactly where the two
engines' assumption strengths differ — one decides what the other calls undecidable. ⇒ ⭐ that is real
evidence, ⛔ not another self-comparison, and it is the seed of the check that could close the export
chain's open assumption limit.

⭐⭐ **AND IT PUTS S9'S BIGGEST OPEN LIMIT WITHIN REACH — ⚠ measured 2026-08-09, ⛔ not yet acted on.**
S9's record closes with *"nothing in the repository performs the cross-engine standard-name lookup"*, and
that is the limit the **light-cone re-point** exploited: one re-pointed name turned the cone into an `ω²`
and **every check in the repository passed**. ⚠ **That sentence is still literally TRUE** — the only
comparator is S10's, and ⛔ nothing reads `mathematica/out/S9_light_requires_shear_mathematica_audit.out`.
⇒ ⭐ **But the instrument now exists and BOTH S9 transcripts are already committed**, so the join S9 says
cannot be done needs **no new instrument** — only a second pair of inputs. ⭐ **This is the cheapest open
item in the rebuild and the one that closes the largest hole.** ⛔ Not on the S10 queue; do not fold it in.

## ⭐⭐⭐ CURRENT STATE — 2026-08-12. ⛔ READ THIS BLOCK FIRST; everything after it is earlier.

### ⭐⭐ THE NAMING QUESTION IS SETTLED. ⭐⭐ NEXT: BUILD THE S11 PY ENGINE.

⭐⭐ **THE USER'S DECISION:** *"Only add `s11_` if it would override something it shouldn't. The point was to
override the times when it should."* ⇒ `directives/S11_export_chain_decisions_v2.md`, **`F9`** — three
outcomes, decided per key by one comparison against the imported `LEDGER`:

| | condition | outcome |
|---|---|---|
| `F9a` | key absent from the import | ⭐ bare key |
| `F9b` | present, objects **proved equal** | ⭐ bare key — ⚠ **this is override-when-it-should** |
| `F9c` | present, equality **not proved** (different **or undecided**) | ⭐ `s11_`-prefixed; ⛔ imported row untouched |

⭐⭐ **Decided under the assumptions the two live objects themselves carry** — ⛔ never a premise this step
asserts on top. ⚠ Measured: the same pair is `UNDECIDED` under the objects' own assumptions and `EQUAL`
under S11's joint premise set, and ⛔ the second reading **overwrites a predecessor whose step never
asserted that premise.**

⭐⭐ **EVERY STEP WRITES ITS LEDGER, S11 INCLUDED.** ⛔⛔ *"S11 writes no ledger"* is **REVERSED** — it removed
S11 from the accumulating record, ⚠ and that record **is** the point of the chain: the list of everything
the model defines and the true list of knobs.

### ⛔⛔ WHAT THIS COST — ⭐ the most useful thing on this page

⚠⚠ **8 review legs, 4 revisions, ~15 findings, and ⛔ ZERO LINES OF ENGINE.** ⭐ The physics — does the
compressional term lift the zero root, and what happens to the mode count — has ⛔ **not been computed
once.**

⇒ ⭐⭐⭐ **RULE 4, verbatim: many turns reasoning toward an answer a script would settle in one.**
⇒ ⭐⭐⭐ **THE LEDGER EXPORT IS NOT THE INSTRUMENT.** ⭐ The instrument is **PY vs WL**, and the WL engine
**imports nothing**. ⛔ A control on the export catches **nothing** about S11's physics — it governs where a
value is **filed**.
⇒ ⭐⭐ **`F9`'s rule survived all four rounds untouched.** ⛔ What broke, every round, was guard machinery
written against **a script that does not exist**, where ⛔ nothing can be ablated ⇒ ⭐ those findings are now
a **BUILD-REVIEW CHECKLIST** in the same file. ⛔ Do not fold them back into prose.
⇒ ⛔ **A leak prohibition was written and CUT** — a denylist means the architecture is wrong; ⭐ blindness is
enforced by **absence**, by bounding what the builder is handed.

### ⚠ Still open, and small

- ⚠⚠ **`S11_SHARED_PHYSICS.md` carries a +13-line UNCOMMITTED amendment, BLOCKED on 2 of 3 items.**
  ⭐ **Keep** the four tag names (`LAGRANGIAN`, `EULER_LAGRANGE_SYSTEM`, `M_RESIDUAL`, `M_RATIO`) — ⛔ without
  them the two engines name the step's primary objects differently and the comparator yields **zero rows**
  for them. ⛔ **Withdraw** the decomposition rule (⚠ it merges a per-pair locus family into one system and
  computes a **different locus**) and the boolean-rendering rule (⚠ defeated by a wrapper, ⛔ penalises the
  stronger zero-test, and ⛔ manufactures a false cross-engine disagreement on **every** boolean row).
- ⚠ `chain_accumulate_or_generate_decision.md` is **RECORDED, ⛔ NOT ADOPTED** — both legs blocked it.
- ⛔ **QUARANTINE:** several leg workspaces computed real S11 physics (a `D=5` determinant, root
  multiplicities) ⇒ `/tmp/s11_fold_leg/`, `/tmp/f9*_leg_*/`. ⛔ They must **not** enter the repo and ⛔ must
  not reach a builder.

---


### ⭐⭐ (superseded 2026-08-12) BUILD S11's ENGINES. ⛔ THE EXPORT CHAIN WAS DEFERRED

⛔⛔ **SUPERSEDED — read the block above.** ⚠ This said *"the S11 SymPy engine writes no ledger"*; ⭐ the user
reversed it: **every step writes its output so the next can import it**, and ⛔ deferring dropped S11 out of
the accumulating record. ⭐ `F9` settles the collision it was avoiding.

| | |
|---|---|
| spec | ⭐⭐ **DONE, `cf4a21a4`** — 1149 lines, two legs. ⛔ It was never the problem |
| export chain | ⭐⭐ **SETTLED — `F9`.** ⚠ Four earlier designs blocked ⇒ `DEFECT_REGISTER.md#c20` carries **why each died**; ⛔ do not re-propose one. ⭐ **`F9` is the fifth and it is the user's, ⛔ not a design** |
| S10 rename | ⚠ **DEFERRED, ⛔ not abandoned** — ⭐ `F9c` prefixes rather than re-points, so ⛔ **no rename is needed to build S11** |

⭐⭐ **WHY THE RENAME WAITS — ⛔ the reason, ⛔ not the difficulty:** the only validation of a rename is that
**the chain still works afterwards**, and ⛔ there is no chain until something writes into it. ⇒ ⛔
Re-pointing **429** committed rows with no consumer that would notice an error is the failure where **every
repo check passes and the physics moves silently.**

⚠ ⭐ **S10 owes FOUR items, and they are ONE pass, ⛔ not four** — the `F3` regeneration, the naming defect,
**12 `dimension_key` refs crossing a `D` boundary**, and **one spectrum written in two vocabularies**
(`k1,k2,k3` in one row, `kx,ky,kz` in the next). ⇒ ⭐ done when a consumer exists that would catch an error.

⛔⛔ **THE LESSON, and it outranks the rest:** ⚠ **not one item in any of the four designs read a PHYSICS
residual** — while `ROOT_DEGREE_RESIDUAL`, `N7_RESIDUAL`, `V7_RESIDUAL`, `Q7_RESIDUAL` and `M_A − M_B` sit
in the spec already, each with **two operands from different routes**, unread. ⇒ ⭐⭐ **Four rounds of
machinery that cannot catch wrong physics.** ⛔ Check the artifact against the physics, ⛔ not against
itself.

---

## ⭐⭐ (superseded) CURRENT STATE — 2026-08-11 (midday)

### ⭐⭐ S11's SPEC IS DONE — `cf4a21a4`, 1149 lines. ⭐⭐ NEXT IS THE **PY DECISION LIST**.

| | state |
|---|---|
| `C17` | ⭐⭐ **CLOSED.** Counts computed **on the component**; each carries `CONSTANT` + certificate · `VARIES` + the **sub-locus where it moves** · `UNDECIDED`; a point survives only under `POINT_EVIDENCE_`. ⭐ Both legs re-ran the register's witness and reproduced it |
| `C18` | ⚠ **PARTIALLY CLOSED — ⛔ say so whenever citing it.** The split is now `ADMISSIBLE` vs `UNDECIDED` **with operands** ⇒ a **coverage gap**, ⛔ not a phantom disagreement. ⛔ **No computation resolves it**, and `§9` states that |
| `T3` witness exchange | ⛔⛔ **CUT** (user, 2026-08-11), struck through in the decisions doc |
| `T7` comparator half | ⚠ **OPEN — ⛔ belongs to the COMPARATOR CONTRACT**, which does not exist yet. ⛔ Not a spec defect |

⭐⭐ **WHY `T3` WAS CUT — ⛔ this is a FINDING, ⛔ not unfinished business, and ⛔ do NOT re-open it without
a new way to establish the correspondence:** ⚠ handing a point from one engine to the other requires the
two to **share a coordinate vocabulary**, ⛔ and the blind build exists precisely so they do not. ⭐ Every
route violated a rule of the same file — a cross-engine **name map** is what this rebuild abolished, each
engine **sorts its own `COEFFICIENT_ORDERING`** so positions do not correspond, and the coefficient list
may ⛔ **never be hardcoded** (`:517-520`).

⭐⭐ **THE TWO MEASURED FACTS THE SPEC NOW RESTS ON — ⛔ do not re-derive:**
1. ⭐ A component's **counts are invariant** under which variable is eliminated (verified both CASes, both
   charts) ⇒ **counts are the comparison rows.**
2. ⛔⛔ **Counts CANNOT see a dispersion relation wrong by a factor of two** — measured: counts stay `1`
   vs `1`, and only the difference **reduced modulo each engine's own defining equations** goes nonzero.
   ⇒ ⭐ that reduction is why component-scoped symbolic payloads are compared at all.
⚠ ⛔ **Which variable an engine eliminates is deliberately NOT pinned** — a round pinned it and **deleted a
branch** (`x·y = 0` keeps `y = 0`, loses `x = 0`).

⚠⚠ **WHAT NINE ROUNDS COST, and the lesson is one line:** ⭐ **thirteen of the sixteen bred defects lived
in a mechanism invented to stop two engines describing a locus differently.** ⇒ ⭐⭐ `CLAUDE.md` **rule 3**
(*name the object, do not specify the recipe*) and **rule 6** (*do not make divergence impossible*) — ⛔ a
review arguing whether a construction is well-defined is answering a question the construction manufactured.

---

## ⭐⭐ (superseded) CURRENT STATE — 2026-08-10

### ⭐⭐ DONE, and gated

| what | commit | gate paid |
|---|---|---|
| **S11 spec REPAIRED** — `directives/S11_SHARED_PHYSICS.md`, **914 → 1005** lines | ⭐ **`ab8cb50e`** | decision list 2 legs · repaired artifact 2 legs · fold verified at source |
| `C17` · `C18` registered | `f87132ef` | — |
| `C19` · `C20` registered, `C19` amended | `62c12f36` · `20f57adf` | — |
| **naming + chain PLAN**, revised after 2 legs | ⭐ **`e60419bc`** | 2 legs · adjudicated by measurement |
| S11 PY decision list — ⛔ **BLOCKED, ⛔ not a directive** | `11bf8e05` | 2 legs; **5 findings, 2 critical** |

### ⛔⛔⛔ THE ROUTE — ⭐ `directives/S11_export_chain_decisions_v2.md` (`4d81e9de`)

⛔⛔ **`S11_naming_and_chain_plan.md` IS SUPERSEDED — ⛔ do not follow it.** ⚠ Its retrofit-first ordering,
its rule-2 charge and its producer-scoping were **all three refuted** (`36589024`, then blocked by 2 legs).
⭐ **There is NO S10 retrofit. ⭐ S10 is REGENERATED, ⛔ not re-keyed, and its tag names do not move.**

1. ⭐ **Storage keys stay FLAT; `D5` is unchanged.** A later step re-deriving an object writes **the same
   key** ⇒ ⭐⭐ **two steps deriving one object must be able to MEET. ⛔ The collision IS the measurement.**
2. ⭐ **Before writing an imported key, compare the OBJECT** (`DEFECT_REGISTER.md:675`, which already said
   so). Same object ⇒ re-derivation: both operands + residual, then guard. Different object ⇒ ⛔ **fails
   loudly.** ⚠ The *"same object"* predicate belongs to the S11 PY list.
3. ⭐ **A re-derived row carries its evidence IN THE ROW.** ⛔ `corroborated_steps` alone is a claim with no
   operands, and `S11:527-529` specifies `Q6r` to **propagate** it.
4. ⭐ **S10's export is REGENERATED under 3** — ⚠ every value byte-identical apart from the new fields, ⛔ or
   it was not a regeneration.
5. ⭐ `C19` is a **REAL deviation** (`S10:197` orders *"Emit `M_A`, `M_B`"*; engines emit `Q2_MATRIX_A/B`).
   ⛔ The record must **disclose** it; ⭐ the rename is its **own gate**, and ⛔ S11's build does **not**
   depend on it — S11's spec fixes its own quantity names.
6. Then S11 PY → WL → comparator → record → card → registers. ⛔ Every stage two legs.

### ⭐⭐ WHY THE PLAN DIED — ⛔ and it is the same lesson three times

⚠⚠ **The user's decision to put the unreviewed plan IN SCOPE in the next review prompt is what caught it.**
⭐ All three unreviewed pieces failed:

| piece | whose | outcome |
|---|---|---|
| **option E** — producer-scoping | ⚠ Codex's own proposal, written inside a review | ⛔ **refuted by Codex**, in the review that put it in scope |
| **retrofit-first ordering** | ⛔ mine, invented during a fold | ⛔ refuted by measurement **before** the legs ran |
| **the 1-vs-8 adjudication** | ⛔ mine, and I overruled a leg | ⭐ **8 reproduces** under a stated proxy — ⛔ but `E2` made the count irrelevant, and that was the wrong way to make it irrelevant |

⛔⛔ **Producer-scoping fails for the reason rule 6 names:** it *"tries to make divergence impossible"*.
⭐ Under `D5` the three overwritten rows meet on one key and are compared; scope them by producer and they
write different keys, **nothing compares them**, and the detector never fires because there is no longer a
collision to catch.

### ⭐⭐ WHAT S10 STILL OWES — ⛔ and the PHYSICS filter applied to each (user directive, 2026-08-10)

> ⭐ *"make sure all changes we pursue are physics related and not process related"*

⇒ ⭐⭐ **Keep an item only if it catches a way the PHYSICS could be wrong or a CLAIM could be unsupported.**
⛔ Nothing below blocks S11's build. ⚠ **Documented so it is not silently dropped** — ⛔ not scheduled.

| owed | verdict | why |
|---|---|---|
| **3 sentences corrected in the S10 record after the last leg reported** (the eight-of-twelve verdict count + 2 companions) | ⭐⭐ **KEEP — CLAIM** | ⛔ A number in a signed record that **no leg has seen**. ⚠ This is the only owed item that could make S10 *say something false*. |
| **`C19` — the emitted names are the engines'** (`S10:197` orders *"Emit `M_A`, `M_B`"*; engines emit `Q2_MATRIX_A/B`) | ⚠ **SPLIT — ⭐ cut STANDS, and ⛔ BOTH LEGS WERE WRONG ABOUT WHY IT DOESN'T** | ⭐ The **forward** discipline is PHYSICS — a name is a binding. ⛔ But S11's spec fixes its own quantity names ⇒ ⭐ **S11 is covered without touching S10.** ⚠⚠ **Both legs claimed the record carries an UNSUPPORTED CLAIM** — *"`M_A`/`M_B` PASS rows while comparator rows with those names are 0."* ⛔ **MEASURED FALSE.** `steps/S10_two_transverse_photons.md:658` says *"each have 12 PASS rows"* and **names the one failing row**; the comparator has **12 PASS + 1 FAIL** for `Q2_MATRIX_A`, `_B` and `_RESIDUAL` alike. ⇒ ⭐ **the number is RIGHT; the deviation is real but produced NO false claim** — it is a **traversal gap**, and the disclosure is cheap ⛔ but not physics. |
| **`F3`/`F4` — regenerate S10's export so re-derived rows carry their operands** | ⛔ **CEREMONY** | ⭐ The residual is already computed, emitted and committed (`out/…:4215`), and `Q6r` emits **its own** derived/imported/difference triple (`S11:527-529`). ⇒ nothing physics-bearing changes. |
| **`F6` — no export published from a partial `MAIN` sweep** | ⚠ **CUT, ⛔ but ONLY with an import obligation** — ⭐ both legs moved this | ⛔ Codex's counterexample: `MAIN, D=5` fails **after** the imported ledger is copied ⇒ the merged export publishes a **predecessor** `D=5` row that a record or card reads as S11 output. ⇒ ⭐⭐ **an importer MUST read `RUN_PAIRS`/`SKIPPED_PAIRS` before using a row** (`S11:895`). ⛔ Without that, the cut is wrong. |
| **S10's requirements registers** (`SUBSTRATE_REQUIREMENTS.md`) | ⚠ **OWED CONTENT, ⛔ not a defect** | ⭐ It captures what the sector obliges; ⛔ no claim is wrong without it. |
| **`DEFECT_REGISTER#F7`** — the kernel equates a boolean with a number | ⛔⛔ **NOT DISCHARGED — ⚠ MY VERDICT HERE WAS WRONG.** ⭐ **KEEP; it moves to the comparator contract (`T7`)** | ⛔⛔ **Measured in the live comparator:** both parsers turn the token `False` into a **native Python `bool`**, `residual_is_zero` is `value == 0`, and ⇒ **`False` vs `0` scores as AGREEMENT**; so does `True` vs `1`. ⭐ Lowercase `false` correctly mismatches ⇒ ⚠ the defect **depends on the spelling in the transcript**, and uppercase is what `str()` emits. ⭐⭐ **NOT REALIZED:** S10's join has **0** boolean payloads in 562 keys; S11's old engines have **400** bool/bool + **11** bool/other in 3042, all scored correctly. ⛔ **But S11's locus protocol requires 3 booleans PER LOCUS in BOTH engines** ⇒ on the critical path. ⚠⚠ **How I got it wrong: I tested `sp.false == 0` — the SymPy object — ⛔ not what the PARSER returns for the transcript token.** ⇒ a component of the path, generalised to the path. |
| **`C20` detector** | ⭐ **KEEP — moves to S11's PY list as `F2`** | ⛔ Not S10 work. |

⇒ ⭐⭐ **NET after both legs: TWO survive** — the three unlegged record sentences, and **`F7`**, which
⛔ **I wrongly discharged** and which now lands in the comparator contract as `T7`.
⚠ **`F6`'s cut is conditional**: an importer must honour `RUN_PAIRS`/`SKIPPED_PAIRS`.
⚠ **The registers row needs re-checking** — a leg reports S9/S10 pass 1 is already marked complete, ⛔ which
would make the row stale rather than owed.

⛔⛔ **AND THE LEGS SPLIT ON FOUR OF THESE ROWS — ⭐ the split, ⛔ not either verdict, is what settled them.**
⚠ On `F7` **both** were wrong in **opposite** directions (Codex *"reproduces"*, Grok *"correctly treated,
closed"*) ⇒ ⭐ neither had measured whether it was **realized**. On `C19` **both** alleged the same
unsupported claim and **both were refuted by opening the record and the comparator** ⇒ ⭐⭐ **a finding
reproduced by two independent legs is still not a measurement** (rule 13).

### ⛔ THE FIVE FINDINGS THAT BLOCKED THE S11 PY LIST — ⭐ the rewrite must answer all five

1. ⛔⛔ **naming rule corrupted the chain 3 ways** — `root_ordering_d3` collision; S10's spec never defines
   `ROOT_ORDERING` (its engine coined it); and **lowercasing contradicts `Q6r:518`**, which pins the map
   case-sensitively (`mu_R`, `B_comp`).
2. ⛔⛔ **schema fork** — the list cited `S9_REWRITE_PLAN#D5`'s sketch `{value, dim, class, step}`.
   ⭐ Measured: **ZERO** entries carry `dim`; real schema is `{class, display, step, value, value_kind}`
   + `dimension_key` (50) + `corroborated_steps` (3).
3. ⛔ `_LOCAL_` split named only `Q6r` — `PREMISE_INVENTORY`, the local name-list tag and solver-condition
   tags fall in **neither** bucket.
4. ⛔ **partial-run semantics undecided** — a skipped `MAIN` cell can leave predecessor rows in a ledger
   that looks valid.
5. ⛔ **the chain-integrity guard is TAUTOLOGICAL** — it compares against the engine's own record of what it
   touched, so an unintended write is classified as touched and passes.

### ⭐⭐ MEASURED, and ⛔ do not re-derive these

| measurement | result |
|---|---|
| `S10_SHARED_PHYSICS.md:195-199` names the route matrices | `M_A`, `M_B` |
| what **both** S10 engines emit | ⛔ `Q2_MATRIX_A` / `Q2_MATRIX_B` — ⭐ and `Q2_MATRIX` is in **neither** build directive |
| S10 spec mentions of the `N2`/`N3` families | ⛔ **0** — engines coined those too |
| S10 py · wl tag suffixes · shared | 4233 · 2983 · **562** |
| S10 `MAIN` `D3` · **S11 engine** emissions · shared keys | 304 · 211 · ⛔ **8** |
| S10 py · wl `emit` CALL SITES | **106** · **137** ⚠ ⭐ *this* is the rename work, ⛔ not 4233+2983 |
| S10 record `:658` claims PASS rows for `M_A`/`M_B`; comparator rows ending `_M_A`/`_M_B` | ⛔ **0** ⇒ a **traceability** gap |

### ⛔ OPEN, ⛔ repaired by nothing above

`C17` (stratum reruns point-dependent) · `C18` (locus protocol pins no construction) ·
S10 registers · the `C20` collision detector (⭐ still needed even under producer-scoping).

### ⚠⚠ PROCESS — what cost the most this round

- ⛔⛔ **THE LEAK TEST IS DERIVABILITY, ⛔ NOT LITERAL PRESENCE.** A decision list stated the *reason* a
  control scale was a `Q6` unknown; subtracting two dimension equations from that reason **fixes the
  answer**. ⚠ One leg asked *"does it print the value?"* → cleared it. The other asked *"can a builder
  derive it?"* → blocked it. ⇒ ⭐ **justification belongs in the REVIEW PROMPT; the decision list states the
  DECISION.** The builder never reads the prompt.
- ⛔⛔ **I CITED DOCUMENTS OVER ARTIFACTS I HAD ALREADY MEASURED** — twice, both critical (the `dim` field;
  the `1` collision). ⇒ ⭐ **every claim in the next list carries the command that produced it.**
- ⛔⛔ **MY ALTERNATIVE COSTINGS WERE BIASED TOWARD MY OWN PROPOSAL** — both errors inflated the options I
  did not choose. ⇒ ⭐ when a plan compares routes, **cost each one with a measurement**, ⛔ not an estimate.
- ⭐⭐ **THE LEGS DISAGREED TWICE AND BOTH TIMES THE DISAGREEMENT WAS THE FINDING.** ⛔ Do not treat either
  verdict as authority — ⭐ adjudicate by running the measurement yourself.
- ⚠ **A fold leaves stale prose ABOVE the fix.** It happened in the spec (`§2`) and then in my own plan,
  one document later ⇒ ⭐ **after folding, re-read the artifact's OPENING sections.**

---

## ⭐⭐⭐ S11 — THE CURRENT STATE, MEASURED 2026-08-09. ⛔ Everything below this block is PROVENANCE.

⚠ **Read this block. ⛔ The older S11 prose further down predates BOTH the S9 and S10 rebuilds and is kept
only for history.**

### ⭐ What EXISTS on disk, with sizes

| artifact | state |
|---|---|
| `directives/S11_SHARED_PHYSICS.md` | **914 lines** = ⭐⭐ **`fc920079`**, which is HEAD. ⛔⛔ **NOT `f49a1684` — that is the 877-line PREDECESSOR.** ⚠ An earlier draft of this block named it and a leg caught it: **repairing `f49a1684` resurrects the tautological point residual, the undefined stratum Jacobian, and the old per-premise emission scheme** that `fc920079` deleted. ⭐ **Pin `fc920079` or HEAD as the repair base.** ⛔ **CLOSED ≠ CORRECT — five defects below** |
| `mathematica/S11_stray_longitudinal_mathematica_audit.wl` | 44,804 B. ⭐ Runs: 3750 tags, 84 s, exit 0. ⚠ **Fix round 3 (`F1`) is UNREVIEWED** |
| `scripts/S11_stray_longitudinal_sympy_audit.py` | 68,039 B. ⛔⛔ **WILL NOT RUN** — `from registry_read import …` at **`:21`**, and `reduction/` is deleted |
| `mathematica/out/S11_…out` · `scripts/out/S11_…out` | 972 KB · 1.15 MB, both committed |
| `_asbuilt/S11_stray_longitudinal_sympy_audit.py` + `README.md` | frozen provenance for two `Q6r` rows |
| `scripts/S11_stray_longitudinal_sympy_audit.premises` | committed |
| S11b-A / S11b-B engines (`.py` + `.wl`) | exist, **old pattern**, rebuilt after S11 |

⛔ **The PY engine does NOT import `S10_exports`** — measured, zero references. ⇒ the dataflow pattern that
defines this rebuild is **absent from S11 entirely.**

### ⛔⛔⛔ FIVE SPEC DEFECTS — ⭐ REPAIR THE SPEC BEFORE EITHER ENGINE IS TOUCHED

⚠ **This list was TWO items until two review legs ran on it (2026-08-09). ⭐ Both confirmed the two;
⛔ between them they added three more and showed both original repairs were UNDER-SPECIFIED.** ⚠ One leg
built a **SymPy witness** for defect 2 rather than arguing it.

⚠ **A spec both engines read is physics-bearing: an error there makes both engines agree on the same wrong
thing** (`CLAUDE.md` rule 7). ⚠ ⛔ **A CLOSED spec is not a CORRECT spec** ⇒
[[feedback-one-engine-fix-is-a-spec-question]].

**1. ⛔⛔ S11 HAS NO INERTIA CONTROL — and this is precisely the hole S10 fell into.**
⭐ Measured: all **seven** packages (`MAIN`, `XFORM_CURLONLY`, `XFORM_DIVONLY`, `XFORM_TRACELESS`,
`XFORM_EXTRA`, `XCOEF_BSCALE`, `XCOEF_BSIGN`) vary **`W`, the stiffness functional, and nothing else**.
The kinetic term `(ρ_br/2)·Σ_j(∂_t u_j)²` is **byte-identical across every package** (spec §2, §7).
⚠ Exact-token search of all 914 lines: **zero** occurrences of `aniso`, `isotropic inertia`, `s_rho`, or a
kinetic-form control.
⇒ ⛔⛔ **S10's headline needed TWO structural premises and its stiffness controls probed ONE.** ⭐ S10 only
found the second because it happened to carry `XFORM_ANISO`. **S11 does not carry its analogue**, so any
S11 result about mode structure inherits an unprobed premise **with no instrument that could reveal it.**
⇒ ⭐ **The object to add is a control that breaks inertial isotropy on one axis while holding `W` fixed.**
⛔ Do not specify its expected effect ⇒ [[project-s10-record-closed]], [[feedback-per-tooth-ablation]].

⛔⛔ **AND THAT SENTENCE IS NOT ENOUGH TO BUILD FROM — a leg showed why.** *"Hold `W` fixed"* leaves
**seven** choices, and holding `W_XFORM_CURLONLY` merely reconstructs S10's ANISO system ⇒ ⛔ it would
**not probe the inertia dependence of S11's own primary `MAIN` claim** (spec `:807`). ⭐ **The repair
decision list must name, as INPUTS:** which `W` is held (⭐ `W_MAIN` if it is meant to police the primary
result); the scale **symbol** and its **domain**; the package's **`D` sweep**; whether the scale is
declared dimensionless or is a **`Q6` unknown**; how it enters **`COEFFICIENT_ORDERING`**; and whether the
coefficient-space and `Q11` loci solve over it.
⭐⭐ **These are ACTION and DOMAIN definitions, ⛔ not expected effects — supplying them leaks nothing.**
⚠ S10 supplied exactly this set (`s_ρ > 0`, `s_ρ ≠ 1`, and `s_ρ` as a `Q6` unknown rather than declared
dimensionless).
⚠ **Watch `:406`, which says "the inertial coefficient" SINGULAR** ⇒ ⛔ a partial repair leaves the new
scale out of `Q6`, `Q10`, the coefficient loci and `Q11`.

**2. ⛔ `C16` APPLIES TO S11 VERBATIM — ⭐ and S11 is where it was always going to matter.**
⭐ S11 **does** specify the right object: `ROOT<r>_N3_STACKED_RANK` / `_N3_TRANSVERSE_NULLITY` at spec
`:321-340`, with an explicit warning against the two-way parallel/perpendicular classification.
⛔⛔ **But §Q8a/§Q8b enumerate strata from the minors of `M_r` and from rank-drop loci of `M_r`** (`:534`,
`:549`) — ⛔ **never from the stacked matrix `[M_r; kᵀ]`, which is what governs `nu_T`.** ⇒ **a locus where
the transverse count ALONE moves cannot be found by construction.** ⚠ Identical to `C16`; ⭐ both S10
engines implemented it faithfully, so ⛔ this is **not** an engine defect and cannot be fixed in one.

⭐⭐ **A leg BUILT THE WITNESS rather than arguing it** — `M(t) = [[t²,−t],[−t,1]]`, `kᵀ = (0,1)`:
generic `(nu, nu_T) = (1, 0)`; at `t = 0`, `(1, 1)`. ⛔ **`rank(M)` never drops**, so `M`'s minors
`[t², −t, −t, 1]` cannot find it; the **stacked** minors `[0, t², −t]` can.

⚠ **SCOPE IT HONESTLY:** `C16`'s exhibited witness was under **ANISO**, and `MAIN` was **not** exposed by a
748-wavevector sweep. ⇒ ⛔ this is **not** evidence that S11's current stiffness-only packages already hide
a `nu_T`-only locus. ⭐ It becomes live **as a COUPLING with defect 1** — once an inertia control exists.

⛔⛔ **AND "ADD THE STACKED MINORS" IS NOT THE REPAIR — a leg measured why.** §Q8b says *"emit the
components your CAS returns"* (`:549`), and **that wording has ALREADY produced a physics-bearing
divergence**: at `XFORM_EXTRA, D = 2`, Wolfram emits `STRATUM_ORDERING: {{beta == 0, muR == bComp}}`
(`mathematica/out/…:561`) while SymPy emits `()` (`scripts/out/…:524`) — ⚠ and the Wolfram point satisfies
the stated domain and reruns `Q3`/`Q4`. ⇒ ⭐ **the repair must define the admissible UNION of root-matrix
and stacked-matrix degeneracy components, their deduplication, and either REQUIRE completeness or
explicitly BOUND the claim.** ⛔ *"Whatever the CAS returns"* cannot support a complete exceptional-mode
statement.
⚠ **`DEFECT_REGISTER.md#C16` allows two legitimate fixes** — add the stacked source, **or** declare `Q8`
silent on pure transversality. ⭐ **The decision list must PICK ONE**, ⛔ or a builder invents a third.

**3. ⛔⛔ §Q6r REQUIRES A DIRECTORY THAT NO LONGER EXISTS — ⭐ found by BOTH legs, and it is the most
concrete of the five.** Spec `:473` still reads: *"The registry is the YAML quantity register under
`research/pde_ledger_v3/reduction/`; read it with that directory's own reader."*
⛔ `reduction/` was **deleted (63 files)**, and `S9_REWRITE_PLAN.md:6` says *"There is no YAML in this
design. No registry."* ⚠ A leg **ran the engine**: `ModuleNotFoundError: No module named 'registry_read'`.
⇒ ⛔ a faithful builder must **restore prohibited machinery, omit required `Q6r` output, or build another
engine that cannot run.** ⭐ **Retire or replace `Q6r` for the export chain. This is a SPEC repair, ⛔ not
an engine task.**

**4. ⛔ §Q3's MULTIPLICITY REQUIREMENT IS SELF-CONTRADICTORY.** `:283` asks for *"the solution set as
returned, retaining multiplicity"* — ⚠ **a solution set does not retain multiplicity.** ⭐ A leg
demonstrated it: for `(x−1)²(x+2)`, `solve` → `[-2, 1]`, `solveset` → `{-2, 1}`, `roots` → `{-2:1, 1:2}`.
⇒ ⛔ the two engines **already invented different constructions** (PY `sp.roots` at `:691-700`; WL factors
and repeats at `:351-361`), so **`ROOT_COUNT_ALL` can differ while both builders follow the words.**
⭐ The spec must name **two separate objects**: a polynomial-root/multiplicity object, **and** the distinct
solution set.

**5. ⚠ TWO MORE FIXED-METADATA TOKENS ARE FABRICATION-FORCING**, in the same class as
`PREMISE_INVENTORY`: `M_ROUTE_RESIDUAL_SCOPE = CODING_CONSISTENCY_ONLY` (`:268`) and `FAILURE_TOKEN:
MISSING_TANGENT_COORDINATES_AND_OFF_STRATUM_EXTENSION` (`:661`). ⛔ **Neither can be read out of a CAS
expression** — the first is an epistemic classification, the second names a *missing* construction.
⇒ ⛔ a literal implementation violates the live-read rule, and manufacturing a holder object to read the
token back violates corollary 5's second sentence. ⭐ **Extend the `:899-901` exemption to cover fixed
metadata explicitly** ⇒ [[feedback-no-fabrication-forcing-rules]].

### ⭐⭐ THE BUILD ORDER — ⛔ each arrow is a gate, ⛔ not a suggestion

⛔ **Rule 7's TRIGGER fires on EVERY decision list below, including the spec-repair list.**

1. ⭐⭐ **SPEC REPAIR FIRST** — the five defects above. ⛔ Orchestrator-written ⇒ decision list gets
   **Codex + Grok** ⇒ repair ⇒ **two legs on the repaired spec.** ⛔ No engine is touched until this closes.
   ⭐⭐ **DECISION LIST DONE AND FOLDED — `dfd95dcb`**, `directives/S11_spec_repair_decisions_v2.md`.
   ⚠ Two legs returned **ten** findings, incl. an **answer leak** (the list stated the *reason* `s_ρ` is a
   `Q6` unknown, and the reason is derivable into the answer) and an item that **would have opened the hole
   it was meant to close** (`FAILURE_TOKEN` is a **field** of a five-field object whose other four are live
   operands; a tag-level exemption licenses typing all five). ⭐ **The fold came out SHORTER.**
   ⚠ ⭐ **The legs DISAGREED twice, and that was the most informative output** — one blocked item 2c with a
   computation the other had cleared. ⇒ `C17`/`C18` registered `f87132ef` rather than absorbed.
2. **PY engine rebuilt** under `S9_REWRITE_PLAN.md`'s pattern — ⭐ imports `scripts/S10_exports.py`, binds
   its objects, writes `scripts/S11_exports.py`. ⛔ Codex-written ⇒ **fresh Claude agent + Grok.**
3. **WL engine rebuilt** — ⭐ imports **nothing**, re-derives independently. ⛔ Same two legs.

⛔⛔⛔ **BOTH ENGINES ARE A REWRITE FROM THE REPAIRED SPEC, ⛔ NOT A REPAIR OF WHAT IS THERE.**
⚠ **User, 2026-08-09:** *"we'll need to rewrite the sympy and wl scripts for S11. It's probably too big for
a simple repair job."*
⭐ **The spec already forces it**: `§5`'s NO VERDICT section says *"Both as-built S11 engines end in a
`VERDICT` tag, and one renders a symbolic boolean as the typed word `TRUE`/`FALSE` with the residual
discarded. ⛔ Neither survives here."* ⇒ **the ending both engines were built toward is deleted.**
⚠⚠ **THE WL ENGINE IS THE TRAP: it RUNS — 3750 tags, 84 s, exit 0 — and that is exactly what will tempt a
patch.** ⛔ Running is not the bar; **being built from the repaired spec is.** ⭐ It also predates
`§8`'s one-tag-per-named-object rule, and the two as-built engines share **exactly one** tag suffix
(`VERDICT`) ⇒ ⛔ **there is no tag namespace worth preserving in either.**
⚠ The PY engine will not even start (`registry_read`, `:21`), so only the WL side invites this mistake.
4. Run both into committed `out/`; commit `S11_exports.py`.
5. **Cross-engine comparator** — ⭐ join on the standard name, as `S10_cross_engine_comparator.py` does.
   ⛔ No pair table, ⛔ no YAML, ⛔ no harness.
   ⛔⛔ **BUT A SAME-NAME JOIN IS NOT BY ITSELF A COMPARISON POLICY — ⭐ FREEZE AND REVIEW THE CONTRACT
   BEFORE THE COMPARATOR SEES EITHER OUTPUT.** ⚠ Otherwise it gets **fitted to the disagreements it
   finds** ⇒ [[feedback-control-outside-the-thing-it-polices]]. The contract must settle: `N6_BASIS`
   (⛔ the spec says at `:342` it is **not cross-engine comparable**) · arbitrary `R0_MATRIX` (⚠ the two
   engines may pick different `R₀`, `:620`) · how **independently chosen stratum points** and their
   point-evaluated reruns are compared (⚠ `:569` — a single point does not characterise a
   positive-dimensional component) · representation normalisation · unparsed-row treatment · the
   shared-vs-`_LOCAL_` tag inventory.
6. Step record → paper card → registers.
   ⛔⛔ **EACH OF THESE GETS ITS OWN TWO LEGS — ⚠ a leg caught that the plan gated the spec and the engines
   and then left this a bare ungated tail.** ⭐ The comparator, the generated `S11_exports.py`, the step
   record, the card and the register edits are **each physics-bearing** and each can propagate a wrong
   interpretation ⇒ rule 7.

⚠ **`W_XFORM_CURLONLY` IS LITERALLY S10's `MAIN` ACTION** — same density, ansatz, phase average, `D` sweep.
⇒ ⭐⭐ **the same physical system computed by FOUR independently built engines**, which the ledger has never
had. ⛔ **Do NOT tell either builder this** — it is a target they could converge on. ⭐ The comparison is
the orchestrator's to run **after** both engines exist.

### ⭐⭐⭐ THE TWO ENGINES HAVE DIFFERENT JOBS — ⚠ user decision, restated 2026-08-09

> ⭐ **"PY should import, the wl scripts derive from scratch. That's the tension that allows us to double
> check that each step is in sync with the others along with verifying one blind build."**

| engine | job | what it establishes |
|---|---|---|
| **SymPy** | ⭐ **imports** the previous step's `LEDGER`, binds, overwrites what it re-derives | ⭐⭐ **CROSS-STEP SYNC** — that this step is consistent with the last one |
| **Wolfram** | ⭐ imports **NOTHING**, re-derives from the spec alone | ⭐⭐ **INDEPENDENCE** — ⛔ this is the **only** blind build |

⭐⭐ **S11's PY ENGINE SEEING S10 IS FINE, AND IT IS THE POINT.** ⛔ Do NOT scope, fence, or filter the
import; ⛔ do not build any apparatus to keep S11 from seeing S10.

⚠ **A review leg proposed fencing it** (2026-08-09), because `S10_exports.py` carries
`transverse_speed_squared_d3` and friends. ⛔ **Rejected — fencing breaks the cross-step check the import
exists for** ⇒ rule 13.

⛔⛔⛔ **AND THE CONTAMINATION WORRY BEHIND IT IS RULE 12, WHICH THIS REPO HAS ALREADY PAID FOR IN DAYS.**

> **"Don't build blindness apparatus. The measured failure is absence of computation, not anchoring."**

⚠ **Measured cost, user 2026-08-09:** *"we've lost days with you doing crazy crap like making full copies
of the repo so there's no contamination."* ⇒ ⭐ **A script that does not compute is caught by §4's
structural rule and by rule 14's ablation — ⛔ both already exist and ⛔ neither needs a new control.**
⛔ Do not invent an acceptance criterion to police a contamination that is not the measured failure mode.

### ⛔⛔⛔ MODE SWITCH — ⚠ S10's card was a DOCUMENT round; ⭐ S11 IS A SCRIPT ROUND

⚠ **The last four gate cycles reviewed prose, where blindness comes from READING ORDER.** ⛔ That muscle
memory is **wrong** for what comes next.

⇒ ⭐⭐⭐ **FOR A SCRIPT, A FORM ABLATION IS MANDATORY IN EVERY REVIEW LEG** (`CLAUDE.md` rule 14;
`.claude/skills/review-legs`). ⛔ Not optional, ⛔ not "if time permits."
- ⭐ **Change the STRUCTURE of a load-bearing object** and report the **literal** diff. ⛔ A **coefficient**
  rescale tests arithmetic; only a **FORM** change tests physics.
- ⚠ **Measured on S11's OWN engine 1: all nine defects were one class — a tag DECLARING what the run used,
  assembled from a literal beside the computation. ⛔ NONE was visible by reading. The ablating leg found
  8 of 9** ⇒ [[feedback-wl-reader-vs-verifier]].
- ⭐ **Demand a script and its literal stdout.** ⛔ A prose re-derivation is the same defect relocated into
  the review ⇒ [[feedback-review-agents]].
- ⛔ **Serialize the legs if both ablate the `.wl`** — the licence has **two** seats and a leg has already
  died mid-run (exit 144) from contention. ⛔ Wrap every kernel run in `timeout 600`.
- ⛔⛔ **DO NOT move engines, directives or transcripts out of the tree to manufacture blindness.** ⚠ An
  earlier draft of this block said to, and **both a review leg and the user rejected it**: it contradicts
  `CLAUDE.md` rule 12 (*"Don't build blindness apparatus… Quarantine is cut"*), and ⛔ it **cannot** work
  anyway, since the PY engine is **required** to import `S10_exports.py` ⇒ [[project-export-chain-pivot]].
  ⭐ **Withhold the explicit target sentence from builder-visible inputs; ⛔ do not CLAIM target blindness
  for the PY engine, and ⛔ do not build machinery to enforce it.** ⚠ The four-engine comparison is still
  valuable — ⭐ as a cross-step consistency check, ⛔ not as four blinded predictions.
- ⚠ **Keep the *rationale* out of builder-visible inputs too** — a leg noted that *"do not state its
  expected effect"* itself signals that an effect is anticipated.

### ⚠ S11's OWN OPEN ITEMS, carried forward

⚠ **Line cites CORRECTED by both legs; ⭐ two rows added; ⭐ one row was wrong.**

| owed | state |
|---|---|
| `PREMISE_INVENTORY` — one tag per `(package, D)` | ⛔ absent in both; WL has **17** `PREMISE_*` suffixes, PY **23**. ⛔⛔ **EXEMPT from corollary 5's live-read at `:899-901`** (⚠ the section opens at `:888`; the earlier cite `:887-901` was loose) |
| the pinned failure object | ⛔ neither emits it. ⚠ **The pinned form is at `:661-669`, ⛔ not `:647`** (`:647` is an ambiguity warning). ⭐ Both emit wrong-shaped stand-ins: WL 11 `Failure[…]` objects, PY 10 zero matrices |
| `STRATUM1_POINT_RESIDUAL` | ⛔ PY emits **10**, WL **0**. ⚠ **The prohibition is at `:559-561`, ⛔ not `:549`.** A **specification artifact**, ⛔ not a physics disagreement |
| ~~`c_s0` unconfirmed~~ | ⭐ **CONFIRMED PRESENT, ⛔ the old row was wrong** — PY at source `:423`, passed to KW sign and zero-locus admissibility at `:1271-1275`; WL builds `bulkAssumptions` at `:709`, used at `:861-866` |
| ⭐ **`M_ROUTE_RESIDUAL_SCOPE`** | ⛔ **ZERO tags in BOTH outputs**, despite the requirement at `:268-271`. ⚠ New row |
| ⭐ **sign payload ORDER** | ⛔ both emit `(operand, token)`; the current spec pins **`SIGN_TOKEN` then `OPERAND`** at `:294-301`. ⚠ New row |
| WL fix round 3 (`F1`) | ⛔ **UNREVIEWED.** ⚠ Recorded in a commit message only; ⛔ absence of review cannot be proven from the tree |
| ⚠ WL *"3750 tags, 84 s, exit 0"* | ⭐ **3750 tags CONFIRMED** in the committed transcript, all 18 cells, none skipped. ⛔ **The 84 s / exit 0 are NOT independently reproduced** — a leg's rerun exited 255 before stdout. ⚠ Treat the runtime as provenance, ⛔ not measurement |

## S11 — where it actually is (⚠ PROVENANCE, superseded by the block above)

### ⭐ Done

| | |
|---|---|
| **Shared physics spec** | `directives/S11_SHARED_PHYSICS.md`, committed **`f49a1684`**. Three review rounds, four legs. |
| **As-built baseline** | tag **`s11-as-built`** (annotated, at `e7658d3a`) — both old engines and the two directives they were built from. |

### ⚠ Where the 2026-08-06 S11 push actually stopped — ⛔ all of it PREDATES the pattern change

| | |
|---|---|
| **Comparator derivative atoms** | ⭐ **CLOSED `681925bb`** — five fix rounds, eight legs. ⚠ The comparator itself lived in `reduction/` and is **deleted**; ⭐ the *finding* is banked below. |
| **Engine directives** | `S11_wl_rebuild_directive.md`, `S11_py_rebuild_directive.md`. Two review rounds, four legs, **seven** blocking defects. ⚠ Both predate `S9_REWRITE_PLAN.md` and ⛔ neither carries the export/annotation pattern. |
| **Registry repoint** | ⚠ **Moot** — the registry is deleted. `_asbuilt/README.md` survives as provenance. |
| **Engine 1 (WL)** | Rebuilt · **3750 tags, 351 `_LOCAL_`, 84 s, exit 0, all 18 cells, 0 skipped**. **3 review rounds / 6 legs · 3 fix rounds · 9 defects.** ⚠ Fix round 3 (`F1`) is **UNREVIEWED**. |
| **Engine 2 (PY)** | Built; output committed (1.15 MB). ⚠ It **imports `registry_read`**, which no longer exists ⇒ ⛔ it will not re-run as written. ⭐ That is expected and is not a defect to fix — it is one of the engines being rebuilt. |

### ⛔⛔ ENGINE 1's NINE DEFECTS WERE ALL ONE CLASS — ⭐ carry this to every remaining engine

⚠ **Every one was §5 corollary 5 or corollary 3: a tag that DECLARES what the run used, assembled from a
literal beside the computation instead of read out of it** — the stripped factor, the simplifier, the sort
key, the bulk premises, `[s]`'s dimension, the supplied action premises, and finally `LOCAL_TAG_NAMES`
itself, the tag whose only job is to inventory the others.

⭐⭐ **NONE was visible by reading. EVERY one was visible by mutation.** ⇒ rule 14, measured: a leg that
ablates finds these; a leg that reasons does not. ⚠ Across three rounds the ablating leg found **8 of 9**.

⇒ ⭐ **For S11b-A / S11b-B / S11b-C: audit every `PREMISE_*` and every declaration tag by mutation before
the first review round**, ⛔ not after. ⚠ And a premise stating an ABSENCE (`v₀ = 0`, no dissipation,
frozen wall width) cannot drive a construction ⇒ corollary 5's second branch is the honest outcome — ⭐ mark
it explicitly, ⛔ do not manufacture a consumer for it.

### ⭐⭐ S11 IS UNHELD — 2026-08-07 (supersedes the hold below)

⛔⛔ **THE PILOT IS DEAD. S11 was held on something that no longer exists.** The S9 comparison-method
pilot (`HARNESS_S9_PILOT_PLAN.md`) was reviewed by four legs across two rounds and the verdict was **do
not build from it** — its central move and its own acceptance criterion were measured mutually
unsatisfiable, and the criterion was passed by a comparator that never compares.

⭐ **S11's engines have both run and their outputs are committed:**
`mathematica/out/S11_stray_longitudinal_mathematica_audit.out` (972 KB) and
`scripts/out/S11_stray_longitudinal_sympy_audit.out` (1.15 MB).

⛔⛔ **S11 does NOT get a `checks_S11.yaml`.** ⚠ Earlier versions of this file said write one, S9-shaped.
⛔ That instruction is dead with the design that produced it — ⭐ S11 is rebuilt under
`S9_REWRITE_PLAN.md`'s pattern like every other step, and it comes **after** S9 and S10 have been rebuilt,
because it imports their `LEDGER`.
⭐ S11 already emits `rho_br`, `mu_R` and `B_comp` dimensions symbolically, so it starts closer to the new
pattern than S9 or S10 did.

⚠ **The fingerprint/semantic-comparison decision recorded below was CORRECTED by review.** ⛔ A fingerprint
is an equality oracle for **already-paired** objects — ⛔ it is **not** a join key and cannot discover
which object corresponds to which. Equal content is not identity: this tree emits many rows with
identical payloads. ⭐ **Declared pairing stays; what dies is the symbol-spelling negotiation.**
⇒ `docs/method_prior_art_findings.md` carries the retractions.

---

### ⛔⛔ THE INSTRUMENT TRACK IS ABANDONED — 2026-08-08, user decision

⛔⛔ **`reduction/` is DELETED — 63 files, commit below.** ⛔ Do not resume it, ⛔ do not "finish"
W0/W1/W2, ⛔ do not restore any of it. ⭐ Git has every byte if it is ever wanted.

⭐⭐ **The replacement is `S9_REWRITE_PLAN.md`.** Read that.

**Why it was abandoned.** The registry reconciles a derivation against a **hand-copied restatement of
itself**. ⛔ No quantity has two independent sources, so every defect this track fought — a declaration
drifting from a derivation, a wrong `D` coefficient the gate cannot see, a witness needing a witness — was
**created by keeping the second copy.** ⚠ Measured: five build rounds, ten review legs, four recurrences of
the same defect one level up, and **no physics.**

⭐ **The requirement is unchanged**: each step consistent with every other, each internally dimensionally
consistent, and every undetermined value tracked. ⭐ The new answer is **dataflow, not reconciliation** —
each step's SymPy engine writes a flat `LEDGER` the next step imports; each Wolfram engine stays siloed and
re-derives, which is where contention lives.

### ⭐⭐ STANDING CONVENTION — emitted names are the STANDARD NAME OF THE OBJECT

⭐⭐ **Every engine emits the name the mathematics uses** — `rho_br`, `mu_R`, `B_comp`, `c_gamma` — ⛔ not a
per-engine spelling. ⭐ A name belongs to the **object**, ⛔ not to the engine that computed it, and ⛔ neither
engine is the authority.

⚠ **Mathematica cannot use underscores in variable names — ⭐ irrelevant.** ⭐ Internal variables stay
whatever WL needs (`rhoBr`); ⭐ only the **string in the emit call** carries the standard name.

⚠ **What the alternative cost, measured:** `checks_S10.yaml` is **3,121 lines**, ~690 rows, **91% a
hand-written name→name pair list**, existing for one reason — the engines named the same objects
differently. ⇒ ⭐ a full day of 2026-08-07/08 went into cross-checking names rather than physics.
⇒ ⛔ **Never hand-author a cross-engine pair table again.** ⭐ Fix the names instead.

### ⭐ CLEANUP — ⭐ DONE 2026-08-08, ⛔ nothing to resume

⭐ **63 tracked files deleted in one commit**: all of `reduction/` (52), the seven root-level
`_review_prompt*.md`, `HARNESS_S9_PILOT_PLAN.md`, `CROSS_STEP_DIMENSION_SCOPE.md`, `INTEGRATION_TODO.md`,
`EXPORT_CHAIN_ARCHITECTURE.md`. ⭐ W0/W1/W2/W3/W4 went with them and ⛔ none of it is owed.

⭐ **Nothing outside `reduction/` imported it at runtime** — measured: every remaining reference is prose in
`directives/` or `docs/`. The `.wl` engines carry **zero** registry references and are untouched.

⭐ **The `out/` files all survive.** ⚠ The plan's sentence *"`scripts/out/` and `mathematica/out/` hold
outputs from the abandoned design"* is **wrong about these two directories** — measured, they hold exactly
six files and all six are engine transcripts (S9/S10/S11 × py/wl). ⛔ Deleting the three SymPy ones would
have destroyed the SymPy side of comparisons whose Wolfram side the same table keeps.

#### ⭐ Retrieving the pair lists after the deletion — **D11 needs them, git has them**

```
git show 67dd3ce2:research/pde_ledger_v3/reduction/checks_S9.yaml    #  12 cross_engine pairs
git show 67dd3ce2:research/pde_ledger_v3/reduction/checks_S10.yaml   # 690 cross_engine pairs
```

⭐⭐ **Read them for ONE thing: the `quantity:` field is already the standard name of the object**
(`factored_determinant`, `transverse_speed_squared`, `inertia_coefficient_dimension`, …). ⇒ **D11's naming
decision is largely already written down** — under the new pattern both engines emit that name directly and
the pair list ceases to exist. ⛔ Do not restore the files; ⛔ do not rebuild a pair table from them.

⚠ Also there: `symbol_naming.rule: registry_snake_case_to_lower_camel` plus six exceptions. ⭐ That is the
**payload** symbol correspondence (`muR` ↔ `mu_R`), which **D12 leaves in place** — WL keeps its internal
spellings and only the emitted tag becomes standard. ⭐ Keep it as a *rule* if a comparator ever needs it;
⛔ never as a hardcoded list.

---

### ⭐ What the abandoned track measured about the physics — ⭐ KEEP, ⛔ this part is not superseded

⛔ **No physics defect was ever found by it.** Three independent parties derived the same verdict on the
`D` question: the engines are correct per branch and the declared `D_brane: 3` matched them.
⛔⛔ **Two blindnesses are MEASURED, not hypothetical, and ⛔ nothing in the new design touches either:**
- **Common-mode:** shifting the length exponent of **all three** brane constituents together leaves **all
  five relations `HOMOGENEOUS`**. ⇒ the gate anchors **differences**, ⛔ never absolute brane dimensions.
- **Shared-spec:** a defect introduced identically into both engines is caught by **no layer**. Both legs
  walked the stack. ⇒ independent re-derivation stays mandatory.

---

### ⚠ SUPERSEDED — the hold, kept for provenance

⭐ **The method decision, taken by the user 2026-08-07** ⇒ `docs/method_prior_art_findings.md`,
[[project-method-prior-art-verdict]]:

- ⭐ **Comparison becomes SEMANTIC, not nominal.** Evaluate both engines' headline objects at shared
  **exact** points (rationals / finite field, ⛔ **never floats** — that is the 1989 caution) and **join on
  a fingerprint, ⛔ not on a tag name.**
- ⭐ The per-step file collapses to **three things**: the symbol→exact-value map, the probe point set, and
  the list of headline objects. ⚠ Tens of lines, against `checks_S10.yaml`'s **3,121**.
- ⛔⛔ **IT DOES NOT FIX FABRICATION.** A typed literal fingerprint-matches another typed literal. ⇒ the
  fingerprint replaces tag **MATCHING**, ⛔ **never** tag **INVENTORY**, and the mutation/ablation controls
  stay exactly as they are. ⚠ Measured: **8 of engine 1's 9 defects were visible ONLY by mutation.**
- ⛔ OpenMath is a paper standard with no working SymPy/Wolfram bridge — ⛔ do not adopt.

⚠ **Our failure has a name**: the **consistent comparison problem** (1986) — an N-version system cannot
reach consensus **when no version has failed**. ⇒ ⛔ **no amount of spec care closes it**, which is why
four spec repair rounds could not.

### ⭐ WHAT S11 STILL OWES — ⚠ REMEASURED 2026-08-07, ⛔ the previous version of this list was WRONG

⛔⛔ **The earlier list said the engines were built against the pre-repair spec. They were not.** Both
engines were **rebuilt and the spec repaired in the same commit** — `fc920079`, which supersedes
`f49a1684`. ⇒ ⭐ **check this list against the spec before using it**; ⛔ a Codex build correctly refused a
directive derived from the old version.

⭐ Measured from the committed outputs and `S11_SHARED_PHYSICS.md`:

| owed | state |
|---|---|
| `PREMISE_INVENTORY` — **one** tag per `(package, D)` | ⛔ absent in **both**. WL carries 17 distinct `PREMISE_*` suffixes, PY carries **23**. |
| the pinned failure object (spec `:647`) | ⛔ **neither** engine emits it — WL has 11 non-pinned `Failure` payloads, PY has 10 zero matrices. |
| `STRATUM1_POINT_RESIDUAL` — spec **forbids** it at `:549` | ⛔ PY emits **10**; WL emits **0**. ⚠ A specification artifact, ⛔ not a physics disagreement. |
| `c_s0` in both admissibility sets (spec `:691-704`) | ⚠ in the spec; engine state **unconfirmed**. |

⛔⛔ **NOT owed — ⛔ do not reinstate it:** *"emit fingerprints at shared probe points"* appears in **neither
the spec nor either output**. ⚠ The fingerprint decision was **corrected by review** — a fingerprint is an
equality oracle for **already-paired** objects, ⛔ never a join key.

### ⛔⛔ `PREMISE_INVENTORY` is EXEMPT from corollary 5 — ⚠ the spec says so at `:887-901`

> *"Corollary 5's live-read requirement does NOT apply to this tag. Several supplied premises are
> qualitative or assert an **absence** — there is no live CAS object to read them from, and one must not be
> manufactured."*

⛔⛔ **Do not require it to be live-read.** ⚠ **Measured 2026-08-07:** the orchestrator wrote a decision list
demanding exactly that, one paragraph after warning that a one-engine fix is a spec question first. ⇒ its
only satisfiable outcome was an **invented** value ⇒ [[feedback-no-fabrication-forcing-rules]].
⭐ Its entries are **declarations**, listed in whatever form the engine holds them.

⚠ **The two deferred engine defects DISSOLVED** — WL's missing §9 density premises and PY's 8 disconnected
premise cells both vanish once premises are one engine-local inventory tag. ⛔ Do not fix them.

### ⚠ The probe set is the next artifact, and it is the one that can go quietly wrong

⛔ A probe point where a headline object degenerates makes the comparison **vacuous while reporting green**.
⭐ Requirements: several **distinct exact rational** points; values satisfying §3's positivity premises so
the point is admissible; per-`D` tables for `D ∈ {2,3,4,5}`; and a comparator check that **flags any
headline object evaluating to zero at every probe point.**

---

### ⛔ SPEC DEFECTS — ⭐ ALL REPAIRED, spec CLOSED at 914 lines (4 rounds, 8 legs)

⚠ **`S11_SHARED_PHYSICS.md` is closed at `f49a1684` and contradicts itself.** ⛔ Both engines read it, so
this is rule 7's failure class exactly.

1. ⛔ **§Q8b line 550 requires `STRATUM<s>_POINT_RESIDUAL`; §5 corollary 3 forbids it.** The point is
   solved **from** the defining equations by `FindInstance`, so the residual is structurally zero —
   measured: changing the action moved the stratum **and** the point, and the residual stayed `{0,0}`.
   ⇒ a check that **cannot fail**.
2. ⛔ **§Q10 does not define differentiation along a stratum** — no tangent coordinates, no off-stratum
   extension. Both legs agree.

⭐ **Authorised (user, 2026-08-06): repair BOTH in one spec pass, then align the engines.** ⇒ remove the
residual from §Q8b; make §Q10 either define the construction or sanction a failure object as the expected
emission. ⚠ Spec repair is physics-bearing ⇒ **two legs before either engine is touched.**

### ⛔ AN ASYMMETRY THE ORCHESTRATOR INTRODUCED — ⚠ the lesson, not just the bug

⛔ **The WL fix directive told engine 1 to drop `POINT_RESIDUAL` on corollary-3 grounds without noticing
§Q8b names it.** ⇒ WL emits **0**, PY emits **5**. ⭐ The finding was right; ⛔ applying it to **one**
engine by directive manufactured a cross-engine disagreement that is a **specification artifact, not
physics**.
⇒ ⭐⭐ **A defect found in one engine is a SPEC question first.** ⛔ Never repair one engine's reading of a
shared clause without asking what the other engine will do with the same clause.

### ⛔ THE SURVEY CORRECTED THIS FILE'S OWN PREVIOUS CLAIMS — verify, do not inherit

The earlier version of this handoff was wrong on three of four points. Measured 2026-08-06:

- ⛔ **"The two engine directives are gitignored and not in the repository."** ⭐ **False.** Both are
  **tracked in HEAD** (`911d0af8`); gitignore never applies to already-tracked files. Nothing needed
  promoting.
- ⛔ **"The namespace is misaligned — WL emits `WL_S11_*`, PY emits `S11_*`."** ⭐ **Understated.** The two
  engines have different tag **vocabularies** and different **granularity** — WL bundles
  `{root → {nullity, orientation}}` into one payload where PY emits three. Stripping engine prefixes, the
  two share exactly **one** tag suffix: `VERDICT`, which the new spec deletes. ⇒ renaming prefixes yields a
  machine-comparable row count of **one**.
- ⛔ **"The `.wl` prints a boolean as `TRUE`/`FALSE` — a rule-2 slip."** ⭐ **It is specified by its own
  directive** (T4, T7, T8, T9, T10). Repairing the engine's line 9 would have fixed nothing.
- ✅ **"No committed outputs, no `checks_S11.yaml`."** Both correct.

⇒ ⭐ **S11 is a RECONSTRUCTION, not a repair.** Emission volume, for scale: S11's engines emit **23** and
**79** payloads; S10's emit **2983** and **4227** against a 690-row declaration.

### ⭐ The remaining order — ⚠ REWRITTEN 2026-08-08, ⛔ the old order's last five items are deleted machinery

⭐ **S11 comes after S9 and S10**, because under the new pattern its SymPy engine **imports their
`LEDGER`.** ⛔ It cannot be built first any more.

1. **S9 rebuilt** under `S9_REWRITE_PLAN.md` ⇒ `scripts/S9_exports.py`. Two legs.
2. **S10 rebuilt**, importing `S9_exports` ⇒ `scripts/S10_exports.py`. Two legs.
3. **Build both S11 engines**, each from `S11_SHARED_PHYSICS.md` plus a short engine-specific directive —
   ⭐ SymPy imports `S10_exports`, ⛔ Wolfram imports nothing. Two legs per engine.
4. Run both into committed outputs under `mathematica/out/` and `scripts/out/`, and commit
   `scripts/S11_exports.py` alongside them.
5. ⭐ **Compare py against wl by tag lookup** — ⛔ no pair table, ⛔ no checks config, ⛔ no harness.
   ⚠ Under **D11** both engines emit the standard name, so the comparison is a join on the tag.
6. Step record · paper card · requirements-register pass 2 · defect-register entries.

⛔ **Items 4–8 of the previous order are gone with `reduction/`**: no `checks_S11.yaml`, no
`harness_ablation.py` battery, no registry provenance to repoint.

### ⚠ What the next session must NOT be told, and must run itself

`W_XFORM_CURLONLY` is **literally S10's `MAIN` action** — same density, ansatz, phase average and `D`
sweep. ⇒ S11's curl-only package should reproduce S10's baseline exactly: the **same physical system
computed by four independently built engines**, which the ledger has never had.

⛔ The cross-step language was **deliberately removed from the spec** during the loosening round, because
pointing a builder at committed S10 rows is a target it can converge on. ⭐ **The comparison is the
orchestrator's to run after both engines exist.**

### ⚠ One trap that OUTLIVES the config it was written for

⛔ **Do not compare `N6_BASIS` across engines.** A nullspace basis is not canonical; S10 produced **11 rows
reporting `DISAGREE` on representation alone.** The spec already says `N6` is display-only.
⇒ ⚠ **This is now a naming decision, not a config one:** under **D11** anything both engines emit under the
same standard name reads as a claim that they should match. ⭐ A non-canonical object must therefore be
emitted under a name that says so, or not paired at all.

---

## ⭐ Banked for the S11 step record — findings, not fixes

- ⛔ **`V7_RESIDUAL = 0` does not validate `V1`.** A leg's own flawed route returned `V1_DIM` 2 against 4,
  and `V2`, `V6` and `V7` were all self-consistent on the wrong `V1`. `V7` tests `V2` against `V6`
  **within** `V1`; only the cross-engine `V1_BASIS` comparison tests `V1` itself.
- ⛔ **`Q.brane.B_comp` and `Q.brane.c_L` still carry provenance into the engine being replaced**, and `R5`
  came from that build. ⇒ the registry is **not an independent operand** for those two rows — only for
  `ρ_br` and `μ_R`. ⚠ **Repointed 2026-08-06**: the artifact is frozen byte-identical at
  `_asbuilt/S11_stray_longitudinal_sympy_audit.py` and all seven loci now address it, so the provenance
  survives the rebuild ⇒ `_asbuilt/README.md`. ⛔ **Freezing it did not make it independent** — a zero
  `Q6r` residual against these two rows still means *the new engine reproduces its predecessor*.
- ⛔ **The dimension provenance for those two rows was MIS-POINTED, and it predates the rebuild.** It named
  `632-664` — the tail of the `A2` check and the start of `A3` — which derives **no dimension**. The
  derivation is `A4`, `676-715`. Corrected in the same pass. ⚠ Two legs computed `A4` independently and
  agreed.
- ⚠⚠ **`B_comp`'s declared dimension is a `D = 3` SPECIALISATION.** `A4` gives
  `modulus_dimension = (2−D, −2, 1)`; the declared `(−1, −2, 1)` matches **only at `D = 3`** (general-`D`
  difference `(3−D, 0, 0)`). `c_L`'s `(1, −1, 0)` matches at **every** `D`. ⇒ ⛔ a `Q6r` residual against
  `B_comp` is a `D = 3` statement and the record must say so. ⭐ Found by one leg, missed by the other.
- ⛔ **The old build's `N_O = 3 for every D` is a SINGLE-ENGINE claim**: the `.wl` has zero reflection /
  `O(D)` content, inside a record section headed *"two engines… agreeing on every computed value"*.
- ⚠ **The S11 pre-registration was wrong on the invariant count and the engines were right.** Not a clean
  pass; say so.
- ⛔ **A `CanonicalDerivative` collision was a false-AGREE channel in the atom every verdict runs on.**
  Found while dormant: 2072 constructions on S10 across 162 names, zero carrying two identities.

## ⛔⛔ Two process findings that cost this session, both now governing

- **A test that "covers" an invariant demonstrates it on one example.** Three times a half-fix passed the
  new test, the whole suite, the full battery **and** produced byte-identical comparator output.
  ⇒ require the builder to build the weaker implementations and show the test fails on each
  ⇒ `feedback_test_must_fail_on_weaker_fixes`.
- **An unchanged regression baseline can be worth nothing.** Measured: `unfixed`, `fixed` and a
  knowingly-wrong variant all produced identical stdout on both configs, because no separator character
  occurs in any identity field. ⛔ Never report "the baseline did not move" as evidence of correctness.

---

## ⛔ OWED ON S9's SIGNED-OFF RECORD — ⚠ found 2026-08-09, ⛔ do NOT fold into S10

⛔⛔ **`steps/S9_light_requires_shear.md:244` calls the full-gradient control *"an ordinary elastic
solid"*.** ⚠ **That imports a substance claim the computation does not make** — and specifically grants
the medium a **shear modulus**.
⛔ **This model's bulk is a SUPERFLUID: zero shear modulus by definition.** ⭐ S9's own record says so
**twelve lines earlier** (`:41-43`): *"A GP/NLS superfluid is a fluid — zero shear modulus. The brane's
shear rigidity therefore requires the substructure."* ⇒ the record contradicts itself across 200 lines.
⚠ **And the brane's shear rigidity is `R-S1-02`, status OPEN**, resting on **one** external document with
⛔ **no executing script**.
⭐ **The control is a STIFFNESS FUNCTIONAL** — `Σ(∂_i u_j)²` against curl-only — ⛔ not a medium; and it is
a **single-modulus** functional, so it is not general isotropic elasticity either (⚠ which has two).
⇒ ⭐ **STANDING: name a control by its FUNCTIONAL, ⛔ never by a substance.** ⚠ The same phrase was in
S10's record and its removal there was correct.

## ⚠ What is still open on S10, deliberately

⇒ ⭐⭐ **SUPERSEDED — the full owed-items ledger, with a PHYSICS verdict on each, is in the `CURRENT STATE`
block above** (*"WHAT S10 STILL OWES"*). ⭐ Exactly one item survives the filter: the **three record
sentences corrected after the last leg reported** (the eight-of-twelve verdict count and its two
companions) — ⛔ they have not been through a leg.
⭐ `DEFECT_REGISTER.md#F7` is **DISCHARGED**: `reduction/` is deleted, and the new comparator was measured
and does **not** reproduce it.

⛔ **Two entries that stood here are DISCHARGED BY DELETION, ⛔ not by fixing:**
`declaration_load_ablation.py`'s 19-minute run, and `D_brane`'s mis-pointed `source_locus` in
`quantities.yaml`. ⭐ Both files are gone; ⛔ neither is owed.

## ⚠ Known limits of the OLD instrument — ⛔ provenance only, the instrument is deleted

⭐ **Kept for one reason: the first two are properties of the PHYSICS OBJECTS, not of the tooling**, so a
rebuild can reintroduce them.

- ⭐ The canonical derivative key used `sorted(set(arguments))`, so a repeated coordinate in a dependence
  list was not distinguished. Multiplicity is not physical; recorded rather than open.
- ⭐ A stiffness form the gradient substitution cannot cover had **no completeness guard**. ⚠ Still true of
  anything that substitutes on a form rather than deriving it.
- ⚠ The rest — exit codes carrying no signal, `_free_symbol_names` dead code, the comparator's
  `--config`/`--output` contract — described `reduction/` and ⛔ died with it.

## Pinned

- **`s11-as-built`** (`e7658d3a`) — both S11 engines and the two directives they were built from, with the
  survey findings in the tag annotation.
- **`s10-as-built`** (`9309da70`) — the S10 spec and both engines as they ran before that rebuild.
- **`wip-2026-08-05-unreviewed`** (`92461853`) — two builds committed without review and reverted. Prior art
  only; nothing in it is known-good.
