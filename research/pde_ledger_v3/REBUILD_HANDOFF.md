# Script rebuild — state and what's next

Read `CLAUDE.md` first (how we work), then this (where we are).
Process detail: `docs/development_pipeline.md`. Defects: `DEFECT_REGISTER.md`.

Last updated 2026-08-08. ⛔ The registry design is ABANDONED — read `S9_REWRITE_PLAN.md` before anything else.

---

## The scope — three steps need rebuilding, and the debt is exactly these

The measured defect was in **three independently-built steps**: engines emitting physics conclusions as
typed sentences with no CAS object behind them, missed by eight review legs.

| step | state |
|---|---|
| S9 | rebuilt |
| S10 | **CLOSED** `e5a2c695` — harness, record, register pass 1, spec repair, card |
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

## S11 — where it actually is

### ⭐ Done

| | |
|---|---|
| **Shared physics spec** | `directives/S11_SHARED_PHYSICS.md`, committed **`f49a1684`**. Three review rounds, four legs. |
| **As-built baseline** | tag **`s11-as-built`** (annotated, at `e7658d3a`) — both old engines and the two directives they were built from. |

### ▶ In flight — 2026-08-06, uncommitted

| | |
|---|---|
| **Comparator derivative atoms** | ⭐ **CLOSED `681925bb`** — five fix rounds, eight legs. |
| **Engine directives** | `S11_wl_rebuild_directive.md`, `S11_py_rebuild_directive.md`. Two review rounds, four legs, **seven** blocking defects. |
| **Registry repoint** | ⭐ Done, two legs ⇒ `_asbuilt/README.md`. |
| **Engine 1 (WL)** | ⭐ Rebuilt and **settled** · **3750 tags, 351 `_LOCAL_` (351 listed), 84 s, exit 0, all 18 cells, 0 skipped**. **3 review rounds / 6 legs · 3 fix rounds · 9 defects.** ⚠ Fix round 3 (`F1`) is **UNREVIEWED** — fold it into the spec-repair review round. |
| **Engine 2 (PY)** | ▶ Building (~3 h; output growing, cells landing in sweep order). |

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

⭐ **S11 is ready to resume and does NOT depend on the instrument track below:**
- ✅ Both engines have run. `mathematica/out/S11_stray_longitudinal_mathematica_audit.out` (972 KB) and
  `scripts/out/S11_stray_longitudinal_sympy_audit.out` (1.15 MB) are on disk.
- ⛔ **`reduction/checks_S11.yaml` was never written.** That is the gap.
- ⭐ Write it **S9-shaped** — 194 lines, 12 declared pairs — ⛔ **not S10-shaped at 3,121 lines**, of which
  91% is a hand-written name→name pair list. S10 is where the apparatus outgrew the physics.
- ⭐ Apply **pattern item 5** (above): S11 already emits `rho_br`, `mu_R` and `B_comp` dimension **laws**
  symbolically, so it is better positioned than S9 or S10 were.

⚠ **The fingerprint/semantic-comparison decision recorded below was CORRECTED by review.** ⛔ A fingerprint
is an equality oracle for **already-paired** objects — ⛔ it is **not** a join key and cannot discover
which object corresponds to which. Equal content is not identity: this tree emits many rows with
identical payloads. ⭐ **Declared pairing stays; what dies is the symbol-spelling negotiation.**
⇒ `docs/method_prior_art_findings.md` carries the retractions.

---

### ⛔⛔ THE INSTRUMENT TRACK IS ABANDONED — 2026-08-08, user decision

⛔⛔ **Everything in the table below is SUPERSEDED. ⛔ Do not build on it, ⛔ do not resume it, ⛔ do not
"finish" W0/W1/W2.** It is kept only to explain what the files in `reduction/` are.

⭐⭐ **The replacement is `S9_REWRITE_PLAN.md`.** Read that, ⛔ not this.

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

### ⛔⛔ FIRST TASK NEXT SESSION — CLEANUP, ⛔ before any building

⚠ **`reduction/` accumulated 20 files in one commit** (`b8ade918`) and it is scaffolding for the abandoned
track. ⚠ 131 files in `directives/`. ⚠ Six stray `_review_prompt*.md` at the ledger root.

⭐ **Delete what served only the registry-reconciliation design.** Candidates, ⭐ verify each before removing:

`dimension_laws.py` · `dimension_law_check.py` · `dimension_law_able_to_fail.py` ·
`engine_dimension_pin.py` · `registry_dimension_witness*.py|yaml` · `registry_import_fence.py` ·
`w3_acceptance_ablations.py` · `w3_duplicate_pin_ablation.py` · `w3_runner_check.py` ·
`w3_weaker_fix_mutations.py` · `w4_pin_completeness_runs.py` · `w4_shifted_registry_printer.py` ·
`w4_weaker_implementation_runs.py` · `test_dimension_laws.py` · `test_registry_dimension_witness.py` ·
`test_registry_import_fence.py` · `test_runner_contract.py` · `W3_REGISTRY_D_LAWS_REPORT.md`

⚠ **Keep** the two engines' own machinery and anything S9/S10/S11 import at runtime. ⭐ Check with
`grep -rl <name> scripts/ mathematica/` before deleting.
⭐ Root-level: `_review_prompt*.md` ×6, `HARNESS_S9_PILOT_PLAN.md` (four legs said do not build it),
`CROSS_STEP_DIMENSION_SCOPE.md` and `INTEGRATION_TODO.md` (both scope the abandoned track),
`EXPORT_CHAIN_ARCHITECTURE.md` (superseded by `S9_REWRITE_PLAN.md`).
⛔ **Commit the deletion on its own**, ⛔ not folded into the S9 build.

⚠ **The one result worth carrying forward** before `engine_dimension_pin.py` goes: five independently-built
engines derive the same symbolic dimensions for `rho_br`, `mu_R` and `B_comp`, and three of them never read
the registry. ⭐ That is a physics corroboration; ⭐ record it in the S9/S10 step records, ⛔ do not keep the
script to preserve it.

---

### ⚠ SUPERSEDED — the instrument track, kept for provenance only

⭐ Front door for the abandoned integration debt: `INTEGRATION_TODO.md`. Scope:
`CROSS_STEP_DIMENSION_SCOPE.md`.

| id | what | state | owed |
|---|---|---|---|
| **W3** | registry declares a dimension as a **law** bound to `Q.brane.D_brane` (`dimension_laws.py`, three brane quantities) | ⭐ **closed with W4:** the typed test duplicate is gone and committed engine operands police a single-sided registry change | production homogeneity still does not anchor the absolute coefficient; shared/common-mode re-derivation remains outside this pin |
| **W2** | registry-dimension witness — compares **emitted** vs **declared** dimension vectors | ⛔ **BROKEN for S10**: a planted engine defect leaves `STATUS_COUNTS` byte-identical; the branch check short-circuits residual adjudication so the 4 S10 rows have **zero discriminating power** | ⛔ **fix round 2** — see below |
| **W4** | rebuild the witness against W3's law-aware registry | ⭐ **closed:** five S9/S10/S11 engine transcripts are parsed; every engine operand, registry-law operand, and residual is printed before its guard | S11-Python's standing `3-D` local residual is a `D=3` specialisation artefact, not a physics disagreement |
| **W1** | dimension type-check on registry substitutions | not started | — |
| **W0** | engines emit `c_γ` itself (**route A**, user-decided) — ⚠ neither engine emits a non-squared speed, so R4 is red at **both** S9 and S10 | not started | one engine at a time, legs each |

The production dimensional-homogeneity gate resolves the declared binding and
checks relation homogeneity; it does not independently determine the absolute
`D` coefficient.  Reduction tests now catch a registry-only coefficient change
by comparing it with fixed committed symbolic engine output.  A synchronized
registry/output change still requires independent re-derivation, and the
medium-sector law/grammar work remains deferred.  The rerun fence is computed
by `python reduction/registry_import_fence.py --list`, not by step number.

#### ⛔ W2 fix round 2 — the owed list, from two legs with executed counterexamples

1. ⛔⛔ **Row grain.** One row per `(artifact, quantity, branch_dimension)`. ⛔ The branch check must **not**
   precede residual adjudication. ⚠ Today all four D-cells fold into one row and any off-registry branch
   forces the whole row, hiding a real D=3 disagreement.
2. ⚠ `BRANCH_DIMENSION_MISMATCH` is decided by **manifest selection alone** — a tautology that fires
   unconditionally. ⛔ It must not be guarded while that is true. ⚠ The committed manifest ships the
   calibration's **perturbed** condition (`dimensions: [2,3,4,5]`) as its baseline.
3. ⛔ **Coverage invariant does not bind** — shrink the manifest to `S9-py` alone and everything is green,
   exit 0, `8 passed` unchanged. Pin the required artifact scope **outside** the manifest.
4. ⛔ `NOT_EMITTED` means *not selected*, ⛔ not *not emitted*. **`Q.brane.B_comp` has zero rows anywhere**
   though both S11 engines emit it, and **S11-py is absent from the manifest** ⇒ **S11 has no dual-engine
   comparison at all.**
5. ⚠ **Multiplier** is correct everywhere but a *coordinated* wrong declaration greens the row where no
   axis cross-check exists (exhibited on S11-wl `c_L`). ⭐ Require a second independent source or an axis
   measurement on any multiplier-bearing source.
6. ⚠ `python -m unittest test_registry_dimension_witness` → **`Ran 0 tests`** (free functions, not
   `TestCase`). ⭐ Also: **there is no runner** — `INTEGRATION_TODO.md` item 1.

⭐⭐ **The engine emits ONE law per quantity, identical in all four D-cells.** The four-value spread is the
**witness's own arithmetic**. ⛔ Do not read it as engine output.

#### ⭐ What is measured about the physics, at its real width

⛔ **No physics defect was found.** Three independent parties derived the same verdict on the D question:
engines correct per branch, registry correct at its declared `D_brane: 3`, and the gap was that nothing
said so — now closed by W3.
⛔⛔ **Two blindnesses are MEASURED, not hypothetical, and W3 does not touch either:**
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

### ⭐ The remaining order — do not reorder

1. Land the comparator fix (two legs, then commit).
2. **Build both engines**, independently, each from the spec plus a short engine-specific directive.
   Two legs per engine.
3. Run both into committed outputs under `mathematica/out/` and `scripts/out/`.
4. Write **`reduction/checks_S11.yaml`**. Two legs.
5. **Extend `harness_ablation.py` to cover it.** ⛔ A checks config with no battery is an unpoliced
   declaration — S10's battery is `ACCEPTANCE 1–19` and it is what let a leg find two gate holes.
6. Run the harness. Only now can it say anything about S11.
7. Step record · paper card · requirements-register pass 2 · defect-register entries.
8. **Repoint the registry provenance** for `B_comp`, `c_L` and `R5`.

### ⚠ What the next session must NOT be told, and must run itself

`W_XFORM_CURLONLY` is **literally S10's `MAIN` action** — same density, ansatz, phase average and `D`
sweep. ⇒ S11's curl-only package should reproduce S10's baseline exactly: the **same physical system
computed by four independently built engines**, which the ledger has never had.

⛔ The cross-step language was **deliberately removed from the spec** during the loosening round, because
pointing a builder at committed S10 rows is a target it can converge on. ⭐ **The comparison is the
orchestrator's to run after both engines exist.**

### ⚠ One trap for `checks_S11.yaml`

⛔ **Do not declare `N6_BASIS` as a cross-engine row.** A nullspace basis is not canonical; S10 has 11 rows
reporting `DISAGREE` on representation alone. The spec already says `N6` is display-only.

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

## ⚠ What is still open on S10, deliberately

- **`declaration_load_ablation.py` runs ~19 minutes on the S10 config.** ⚠ Under the clarified rule this is
  **acceptable if it prints incrementally** — the bar is observable progress, not elapsed time
  ⇒ `feedback_script_timeout_policy`. Check that it flushes; ⛔ do not narrow it.
- **Three sentences corrected after the last leg reported** — the eight-of-twelve verdict count and its two
  companions. They have not been through a leg.
- `DEFECT_REGISTER.md#f7`'s owed measurement. Re-scope it or discharge it.
- ⚠ **A second registry provenance dent**: `D_brane`'s `source_locus` (`reduction/quantities.yaml:177`)
  points at S9 for a quantity S10 introduces. Fix with the `B_comp`/`c_L` pass.

## Known limits of the instrument

- The canonical derivative key uses `sorted(set(arguments))`, so a repeated coordinate in a dependence list
  is not distinguished. Multiplicity is not physical; recorded rather than open.
- A stiffness form the gradient substitution cannot cover has no completeness guard.
- Exit code carries no signal on either config; both exit 2 in the healthy state. **Judge by counters.**
- `_free_symbol_names` has **no call site anywhere in the repository** — pre-existing dead code.
- ⚠ The comparator requires `--config` **and** one `--output ENGINE=PATH` per engine. It resolves the
  registry root from the config path, so a bare copy of `reduction/` fails with `RegistryValidationError`.

## Pinned

- **`s11-as-built`** (`e7658d3a`) — both S11 engines and the two directives they were built from, with the
  survey findings in the tag annotation.
- **`s10-as-built`** (`9309da70`) — the S10 spec and both engines as they ran before that rebuild.
- **`wip-2026-08-05-unreviewed`** (`92461853`) — two builds committed without review and reverted. Prior art
  only; nothing in it is known-good.
