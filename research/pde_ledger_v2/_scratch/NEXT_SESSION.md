# NEXT SESSION — handoff (authored 2026-07-30, after the DIRECTION CHANGE)

## ⛔ READ BUDGET — the largest single cost of the last two sessions

Read **this file** and **`STATUS.md`**. That is your whole orientation budget. Delegate everything else:
send an agent with a specific question and require **"≤20 lines + file:line loci"** back. ⛔ Do not open
`docs/development_pipeline.md`, `manifests/DIMENSION_REWRITE.md` or the stage notes "to get oriented" —
they are 300, 1500 and 1000+ lines. Re-read a *section* only when about to act on that exact section.

⭐ **ONE NAMED EXCEPTION — read it before resuming, not "to get oriented":** `manifests/DIMENSION_REWRITE.md`
**§1b (the D1–D5 decisions)** and **§3b (what those decisions REOPENED)**. Several recorded conclusions —
three waivers, four "impossible" stages, a coverage estimate — held only under constraints since lifted,
and will otherwise read as settled.

---

## STATE

Branch `ledger-v2-rebuild`. ⛔ Do not trust a hash written here — run `git log --oneline -5` and
`git status` first. **Dimension rewrite: 7 of 30 converted** (004, 011, 012, 013, 016, 018, 023).

---

## ⭐⭐ THE DIRECTION CHANGED — read this before doing anything

The verification apparatus had grown four layers above the physics: a review instrument, an acceptance
gate for it, a byte-authority over that gate, a documented procedure over the authority. A full session
was spent at those layers verifying **no physics**. The user stopped it.

**The governing test — now the first thing `docs/development_pipeline.md` states:**
- Catches a way the **PHYSICS** could be wrong? → **keep.**
- Only a way the **TOOLING** could be wrong, or a motivated adversary? → **cut.**

⛔ **Immutability is not a discipline here.** Files are freely editable; nothing is frozen, byte-perfect
or under custody. ⛔ Role separation is **one builder, one fresh reviewer**, plus a **second independent
review leg on physics-bearing artifacts only**. No three-session sequences, no round-robin, no four-model
bookend, no per-chunk user gates.

⚠ **You will find committed artifacts that predate this** — a conformance fixture suite, a freeze
authority, a re-acceptance procedure. They are inert and paused. ⛔ Do not build on them; do not spend a
session deleting them either unless asked.

---

## ⭐⭐ THE FINDING THAT MATTERS MOST — do not re-derive, do not soften

**Dual-engine agreement is VACUOUS where both sides are hand-declared literals.** ⚠ State it at the width
the measurement supports: **no dimensional input enters from outside the stage's own typed declarations.**
stage023 emits **29 records — 22 typed `SOURCED_DIMS` plus 7 live `dim_of` walks**, and those walks run
over exactly those 22 literals and nothing else. stage016 emits **21 — 12 typed rule-table entries plus 9
computed** by `dim_of` over real expression inputs, again sourced only from its own declarations.
⛔ It is **not** "29/29 literals" and **not** "zero computed" — both overstate; ✅ **`STATUS.md` was
corrected 2026-07-30 to these same counts**, so the two files now agree. So the comparator on these
stages catches a *transcription split between two typed copies of the same numbers*, and nothing else.

⚠⚠ **THE COUNTER-MEASUREMENT — read it WITH the finding, not after it.** Per `STATUS.md`, relabelling
stage016's basis leaves **that stage's own 82 assertions blind at exit 0** (82 PASS, printing
`measure: 'M^3'`) while the comparator catches **18 of 21**. ⇒ The comparator is the **sole instrument**
between a converted stage and a relabelled basis. ⛔ **Do not read the finding above as licence to cut it
under the governing test** — a relabelled basis is a way the PHYSICS goes wrong, not the tooling.

⇒ **Nothing that survives independently RE-DERIVES the physics** outside one fresh agent. The cross-engine
gate shows two implementations agree, never that either is right — and git shows all 43 pairs were first
added together, so even *independence of route* is discipline, not evidence.

**Hence the front (user decision, 2026-07-30): a DERIVED-vs-DECLARED census** — how much of this ledger
computes something rather than asserting it. Nobody can currently say.

⚠ **It is NOT a binary, cannot be measured as one, and must not be attempted as one.** A value can be
physically derived and still be *stored* as a literal. The real states are ~six (derived-in-form ·
executed derivation · branch-determined · input · gap · benchmark), and the method needs one stable ID
**per occurrence**, joined across stages by **adjudicated** identity, not by candidate name. ⚠ **Name
grouping is precisely what already fails:** `research/pde_ledger_v2/CANONICAL_DIMENSIONS.md` §"GROUPING
LIMITATIONS" lists **eight known-same cross-stage pairs it does not group** (011 CamelCase vs 012
snake_case: `CsSquaredDim`, `CorruptKDim`, `EnergyDim`, `FourVolumeDim`, `MassDim`, `OmegaDim`,
`PressureDim`, `RhoDim`), and 016's `T_w` → candidate `TW` never meets 013's `Tw` → `Tw`, so it cannot
see that pair at all.

**Substrates that already exist — the census EXTENDS them, it does not start clean. ⚠ Each is stated
WITH its limit, because extending them does NOT by itself enumerate documentary imports, verdicts,
controls, assertions, or non-dimensional computed objects:**
- `software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md` — the semantic classes
  (DERIVED / INPUT / gap / benchmark) already assigned per value. ⚠ **Limit: it covers a limited Path-A
  constant set** (its own scope line: *"the whole Path-A constant set"*), not the ledger's objects.
- `research/pde_ledger_v2/notes/parameter_register.md` — dependency and reduction edges per knob.
  ⚠ **Limit: it aggregates PARAMETERS and REDUCTION EDGES, not occurrences** — one row per knob, and its
  provenance classes are a **different seven-class** taxonomy (see P0-1).
- ⭐ **`research/pde_ledger_v2/notes/stages/ledger_stage043_irreducible_count_range.md` §10 — the
  `REGISTER_TO_COUNT_MANIFEST`, a BOUNDED ~152-ID manifest with its OWN identity and category system**
  (mutually-disjoint categories, engine-qualified row/knob IDs, exact total stated and asserted in
  `scripts/ledger_stage043_irreducible_count_range_sympy_audit.py`). ⛔ The census must **reconcile with
  it, not reinvent it** — an ID scheme and a category partition at corpus scale already exist here.
- the **§4-c1 per-quantity route tables** in the stage notes — per-quantity CORRECT/WRONG/UNDETERMINED
  *with its route*, plus the count of records that are declared literals in both engines. ⭐ **stage023
  already demonstrates the method** — **34** routed objects and two scoped tallies (**24/0/10** as
  identifications against the corpus, **27/0/7** inside the stage-local closure); **stage016** carries a
  tracked **21/0** verdict, **12** typed in both engines and **9** computed. ⚠ **Limit, twice over: these
  cover DIMENSION-VALUED quantities only, and tracked tables exist for exactly TWO stages — 016 and
  023.** ⛔ And their tallies are **dimensional-correctness verdicts, not provenance states** (see P0-2).
- ⚠ **Part VII's stage046 firewall requires a calibration map of this kind** — *"every constant
  DERIVED/INPUT/gap/benchmark"* (`research/pde_ledger_v2/notes/part7_integration_atomic_split.md`, the
  046 row). ⇒ The census is not new work bolted on; it **overlaps** a Part VII obligation. ⛔ **State the
  overlap, do not claim equivalence:** the same row also consumes **stage043's range** and the corrected
  R1 inventory, 046 assembles from the **044→045 spine** (same file, the **"Ordering (ratified)"** block:
  *"046/047/048 assemble from the spine"*), and its requirement is over
  **constants**, whereas the proposed **occurrence** census is broader than what 046 needs.

---

## ▶ THE WORK

⭐⭐ **DECIDED (user, 2026-07-30): the CENSUS is the front; the conversions continue BEHIND it.** ⛔ The
question is closed — do not re-ask it. ⛔ Do not delete the reasoning that made it a choice, because it
is why the choice was available: the conversions are the original minimal ask (one shared import so every
script's dimensions come from one place), the census is what says whether the physics is derived at all,
and **they are independent** — the census does **not** depend on conversion, so it can run against any
stage's existing declarations.

⭐⭐ **The finding that motivates the census, WITH its counter-measurement — read them together; either
one alone produces the wrong conclusion.** ⇒ **Nothing that survives independently RE-DERIVES the
physics:** the cross-engine gate shows two implementations agree, never that either is right, and git
shows all 43 pairs were first added together. ⚠ **AND the comparator remains the sole instrument between
a converted stage and a relabelled basis** — relabelling stage016's basis leaves that stage's own **82
assertions blind at exit 0**, while the comparator catches **18 of 21**. The first half is why the census
is the front; it is ⛔ **not** licence to cut the comparator.

⛔⛔ **THE WAY THIS CENSUS FINISHES AND RETURNS A WRONG NUMBER — read it before writing the schema, not
after the fan-out. Counting SYNTACTIC EXECUTION as PHYSICAL DERIVATION inflates the "computed"
fraction.** Measured, not feared: of stage016's **nine computed records, three are self-referential** —
`actual_M2`, `actual_K2`, `actual_K2_over_M2` walk a declaration back to the constant that defines it —
and **none of the nine is computed from any physical input**
(`research/pde_ledger_v2/notes/stages/ledger_stage016_l2_so3_covariance.md` §1.6(3); all nine are pure
functions of the twelve `dim_rules.*` declarations). stage023's **seven `dim_of` walks run over
`SOURCED_DIMS` and nothing else** — the stage's own 22 typed literals
(`scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py`, `run_dimension_check`). In both
cases a code path really ran and **no physics entered**.
⇒ **The EXTERNAL-INPUT CRITERION IS A REQUIRED FIELD OF THE SCHEMA, not a caveat in the prose:** every
row must answer *does this value trace to a field equation, or to an input from outside the stage's own
declarations?* ⛔ A row that cannot answer it is **not** counted as computed. Without that field the
census returns a flattering number — and the number is what gets quoted.

**P0 — THE CENSUS: spec → pilot → fan-out, in that order.**
1. **Specify the schema.** ⛔ One builder, one fresh reviewer; the census is physics-bearing, so it gets
   the **second independent review leg**. The spec must settle, at minimum:
   - the **states**, the **per-occurrence ID**, the **join rule**, and **what evidence each row cites**;
   - ⭐⭐ the **external-input criterion above, as a required field** — see the block immediately above;
   - the **census universe** — which objects are in scope at all;
   - what counts as **one occurrence**, and the **inclusion rules**;
   - the **denominator** the fraction is reported against, and the **aggregation rule** across stages;
   - an **ID registry** — where the IDs live and who mints them.
   ⚠ **The states OVERLAP, and the spec must resolve that rather than inherit it:** BRANCH-DETERMINED
   may *already* be computed, and `research/pde_ledger_v2/notes/parameter_register.md` classifies with a
   **different seven-class** taxonomy (`ACTION` / `DERIVED` / `CONV` / `FREE-UNREDUCED` /
   `CALIB`(`-ANCHOR`) / `CANDIDATE` / `GAP`). ⇒ give an explicit **precedence rule**, or split
   **provenance** and **execution** into **two axes**.
   ⚠ **The universe question is already consequential, not hypothetical:** stage023 carries **42** source
   objects under its declared membership rule (§1.6, `ENUM_COUNT|row_count|42`), **29** emitted records,
   and **34** physics-routed objects (§1.7(1)) — three different candidate universes inside one stage.
2. **Pilot on THREE stages — 023, 016 and 043.** 023 and 016 are the two that already carry §4-c1 route
   tables. ⭐ **043 is the third** because it exercises **DERIVED-in-form versus unexecuted debt**
   directly (its §11 keeps `R35` labelled `DERIVED-in-form` rather than flipping it to `PENDING-debt`)
   and it computes over a **corpus-scale bounded manifest** rather than one stage's records. Check the
   schema survives contact with all three.
   ⛔ **Do NOT require the pilot to reproduce 023's `24/0/10` and `27/0/7` or 016's `21/0`.** Those are
   **dimensional CORRECT/WRONG/UNDETERMINED verdicts — a different axis from provenance**, so a
   provenance schema has no reason to land on those distributions and reproducing them is not evidence
   about the schema either way. What the pilot **may** reuse is those tables' **row universes** and their
   **correctness fields**; the tallies are context, ⛔ not an oracle.
3. **Then fan out.** ⚠ Do not fan out before the pilot; a schema defect multiplied across 43 stages is
   the expensive failure here.

⏸ **P1–P4 are the CONVERSION work. It continues BEHIND the census** — recorded order and hazards intact,
⛔ but none of it is the next step any more.

**P1 — Pass-1b, small.** `composite_build.py` recovers dimension order by finding an AST class named
exactly `Dim`; the shared module exports `Dimension`, and the checker has **zero** references to
`ledger_dimensions`. Every conversion therefore *removes* that stage's recovery. **Verified: recovery is
10 of 43** — exactly 7 scripts carry an exact `class Dim` (005, 006, 007, 008, 009, 030, 031) plus 3
registered bare-tuple digests (032, 038, 042) — and **all 7 converted stages carry none**. Teach the
checker the module before converting more. ⏸ **Still worth doing whenever conversions resume, for exactly
that reason — but no longer urgent, and ⛔ not to be presented as the next step.**

**P2 — two conversion accelerators, both small, no ceremony. ⭐ Order per `STATUS.md`: (1) the ablation
driver, (2) the `DIM|` emitter — and both land BEFORE 027 begins.**
- **(1) the ablation driver** — mutate a declaration, confirm the declared assert fires, record it;
  reviewed by one fresh agent. Current trimmed spec:
  `research/pde_ledger_v2/notes/ablation_driver/REQUIREMENTS.md`. ⛔ `CONTRACT.md` and the `fixtures_v4/`
  suite in that same directory are **superseded history**, not authoritative.
- **(2) the shared Mathematica `DIM|` emitter** —
  `research/pde_ledger_v2/notes/wl_emitter/REQUIREMENTS.md` (v2, reviewed CLEAN). **One builder writes it,
  one fresh agent reviews it.** No contract, no fixtures. It addresses five measured hazards each stage
  currently re-derives: dead code after the terminal `Exit[]` (all 43 `.wl` files end in one), a derived
  label beside raw storage order, duplicate/corrupt records from a holder invoked many times, hardcoded
  exponent literals, fractions needing exact serialisation.

**P3 — the conversions.** Order: the stage027-shape decision, then 027, then 021 (heaviest). ⚠ 027 is
MIXED **in one respect only** — its single computed vector never reaches top level (it dies in `runAll`'s
`Module`). ⛔ Its **16 declared per-symbol `baseDims` vectors ARE already top-level** (`027 wl:183-191`),
so it is **not** forced to a 1-row `DIM|` stage and the `.wl` route does **not** fail to produce
per-symbol vectors. ⚠ 021's renderer reorders **display** to `L,T,M` (`wl:125`/`:139`) with each label
carrying its own exponent — **there is no M↔T value swap**, and the storage order `(L,M,T)` is stated
in-file (`wl:342`/`:384`/`:528`). The only hazard is scraping that display and labelling the sequence
`axes=L,M,T`. ⚠ 035/036 are unprototyped; **fix the canonical-table generator before 035/036/037** — it
*raises* on cross-engine axis-order disagreement, the pattern the 037 **spike** exhibited (the committed
037 `.wl` declares no axis order at all — see `STATUS.md`). ⭐ **Measure whether 027 comes in dramatically cheaper than 023.** If not, say so
plainly rather than grinding through 22 more.

**P4 — the manifests, trimmed.** 4 of 44 exist. Semantic core = quantity identity, dimensional relations,
dependency graph and cycles, consumption completeness, mutation — **plus C2/C3/C5/C8**, which were
wrongly cut once and restored: C2 detects changed consumed equations, C3 consumption of retired physics,
C5 wrong irreducible-count ranges, C8 calibrated/target-matched genesis.

**P5 — standing, none of it blocking.** The seven `WORK-023-*` adjudications, sharpest being the **`D0`
seam** — real, but ⚠ **stage017 has no dimension machinery in either engine**, so it "implies" nothing on
its own: the `M T⁻²` reading comes from `research/pde_ledger_v2/notes/parameter_register.md:185` (whose
ℓ=2 radial scalars are computed in **stage016's** dimension rules, `EXPECTED_K`), while 021/023/027
declare `M L⁻¹T⁻²`, and the propagation contradicts an asserted dimensionless target.
`DIMENSION_REWRITE.md` §12 states this correctly — cite it rather than the compressed form.
⛔ **Dimension libraries: CLOSED, NO-GO — do not reopen.**

---

## ⭐ ALREADY MEASURED — do not re-derive

- **Recovery is 10 of 43**, not the ~16 both front-door docs claimed until 2026-07-30.
- **The module digest cannot establish "you edited the module and did not re-run the stages"** — it
  compares hashes, and `--accept` resets without checking any producer.
- **"Consistent by construction" is wrong**: one shared import gives *representation unity*, not correct
  dimensions. Two stages can declare the same wrong exponents from one module.
- **U13** — a mis-binding *inside* an exponent-degeneracy class is invisible to the stage gate, the
  comparator and a re-run: 25 of 29 stage023 records share their triple, and five simultaneous same-class
  rebindings gave exit 0, 111 PASS and a green comparator. ⚠ A **source read does** reach it.
  ✅ **Evidenced — the artifact is retained** (rescued 2026-07-30 out of gitignored `_scratch/stage023_h/`
  to `notes/stage023_step_h_evidence/u13_u14_rescued/`): `u13_expE_same_class_rebinding.diff` carries
  exactly the five rebindings (`M0←R0`, `K0c←Z1ret`, `g_U←g_W`, `eta_null←q_free`, `T0←P0_physical`) as a
  five-line diff against the committed script, with `u13_expE.stage.out` (exit 0, 111 PASS) and the two
  `dimensions.txt` payloads, which reduce to the same md5 — only the header's `source_sha256` differs.
- **U14** — jointly, **20 of 22 declarations are unpinned by anything in the Python**; only `M0` and `D1`
  are constrained, and only against declared `EXPECTED_DIMS`. A one-at-a-time ablation shows 6, because
  moving one breaks a relation that moving them together preserves. ⚠ **Scope — do not over-read this as
  workstream-wide blindness.** The cross-engine comparator **does** catch the joint mutation:
  `RESULT|stage=stage023|status=FAIL|mismatches=16` (`u13_u14_rescued/u14_i7_joint.cmp.out`), because the
  `.wl` carries an independent literal `baseDims` table. U14 bounds the **Python** gate, not detection as a
  whole. (The genuinely undetected case is **U13**, whose payload is byte-identical.)
- **The dimensional gate is not what catches a corrupted declaration** — `base_verdict` ranks
  `dimensional` third, so the selector control or the baseline earned-verdict check fires first, in all
  16 caught rows.
- ⛔ **Open, module-level:** `ledger_dimensions._exact` **coerces** an `sp.Float` instead of rejecting it,
  so a float typo becomes an exact-but-wrong rational silently (`Dim(0.1,0,0)` ≠ `1/10`). Touches all
  seven converted stages.
- ⛔ **`run_all_audits.sh` prints `Fail: N` but EXITS 0.** It gates on the module pin only and never
  invokes the comparator or the generator (zero references to either). ⇒ Trusting its exit code is a
  **silent pass**; read the `Fail:` tally, not `$?`.
- **An agent CAN run Mathematica here** — a review agent ran `timeout 600 math -script …` to exit 0 and
  byte-diffed the transcript itself. ⚠ One agent, one occasion; the ≤2 concurrent seat cap must count
  agents too.

## ⛔ TRAPS — measured 2026-07-29/30

- ⭐ **A reviewer will read a prior verdict left in the working tree and echo it instead of re-deriving.**
  Measured this session: a confirm-pass reproduced an earlier verdict verbatim, session footer included.
  ⇒ Move prior verdicts **out of the project tree** before launching a review. Absence beats instruction.
- ⭐ **Verify a reviewer's RATIONALE separately from its FINDING.** Three times this session a true finding
  arrived with an unsupported reason attached, and twice that reason was relayed onward as fact.
- **Ask an applier for a COUNT, and change the METHOD.** A keyword sweep reported 1 stale instance; the
  same agent reading by tense found 3; an independent party found a 4th; widening to the *class* found 9.
- **An unanchored `grep` for a completion marker** matches that marker quoted inside another log. Anchor
  it: `grep -q '^___CODEX_DONE___'`.
- **A reaped waiter is not a finished job, and a killed waiter is not a dead job.** Check the artifact.
  Never `pgrep` by name — the waiter matches itself.
- **A commit cannot cite its own hash.** Cite the path; hashes are for *other* commits.
- **Line-number citations decay whenever a file grows** — the seven `WORK-023-*` loci are stale in **both**
  citing files, uniformly *within* each file but by **different** amounts: `STATUS.md` cited **+37 too
  high**, `parameter_register.md` **−62 too low**. ⛔ 99 was the *spread between the two files*, never one
  offset. ⚠ Both sets were stale at HEAD `968921dd` and are corrected in the working tree (uncommitted).
  Prefer a named anchor.

## ▶ OPEN DECISIONS for the user

⛔ **"Census or conversions first" is CLOSED** (user, 2026-07-30 — the census, see ▶ THE WORK). Do not
re-open it as a question.

1. **The numerical solver's freeze hash** — deliberately untouched by the scale-back. There the freeze is
   the *mechanism* of target-blindness, and post-hoc refitting is a way the **physics** goes wrong, so it
   plausibly passes the governing test. ⛔ Undecided; do not strip it by citing the scale-back.
2. **Whether to commit the `_scratch` review records** — gitignored, will not survive.

## OPERATING MODEL

One builder, one fresh reviewer; a **second** independent leg on physics-bearing artifacts only; the
**physics leg stays blocking**. Decision list → applier agent → **you review the diff**. Edit, never
rewrite. Agent reports ≤40 lines + a file path. Re-run every named acceptance command yourself and read
its literal exit code. ⛔ Commit only when the user asks.

⭐ **Three PER-STAGE steps the scale-back KEPT — `STATUS.md` still carries them and this handoff had
dropped them. They are physics-catching, not ceremony:**
- **The PHYSICS leg runs FIRST and is BLOCKING** — a NAMING decision is a physics decision, and it is
  derived from the model, never checked against a claim.
- **Regenerate the `.out` and BYTE-COMPARE it.** An uncaught `Throw` exits 0 with an empty or truncated
  transcript, and **only** the byte-compare catches it. ⛔ No longer an orchestrator-only duty — whoever
  runs the stage does it.
- **The co-authorship guard** — the party that wrote the `.py` must **not** be the party that adjusts the
  `.wl` until the comparator agrees. Tuning whichever side disagrees is the shortcut that resembles a
  pass, not a fix.
