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
| S10 | ▶ **CHAIN + COMPARATOR BOTH DONE** — export `e644876c`+`c84263ed`, comparator `82443c95`. **20 legs, ⛔ no physics moved.** ⭐⭐ **NEXT: the STEP RECORD** (stale: 664 lines, 23 dead citations), then card → D12 → registers. `.wl` untouched |
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
| paper card `paper/steps/S10_two_transverse_photons.tex` | ⭐ **targeted edit, ⛔ NOT a rewrite** — 2 stale refs in 299 lines; physics verified unchanged |
| requirements + defect registers | not started |

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

## S11 — where it actually is

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

- **Three sentences in the S10 record were corrected after the last leg reported** — the eight-of-twelve
  verdict count and its two companions. ⛔ They have not been through a leg.
- `DEFECT_REGISTER.md#f7`'s owed measurement. Re-scope it or discharge it.

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
