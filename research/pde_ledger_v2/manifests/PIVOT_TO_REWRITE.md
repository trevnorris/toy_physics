# PIVOT — from "survey and accommodate" to "rewrite and gate" (2026-07-26)

⭐ **This supersedes the survey-first approach in `LEDGER_WIDE_PLAN.md` §2.1–§2.2.** Read this before
acting on that plan. Decision made by the USER after two measurements (below); the reasoning is
recorded here so it is not re-litigated.

---

## 1. THE DECISION

**Old plan.** Survey all 43 SymPy audit scripts to discover every dimension idiom they use, so a shared
module could be designed to *accommodate* them, then refactor onto it with "identical PASS counts" as
the behaviour-preservation gate.

**New plan.** Treat the scripts as **output to be regenerated**, not as a specification to be satisfied.
Design the shared dimension module as a **choice**, rewrite each script's dimension handling to conform,
and gate the rewrite on a **cross-engine comparison harness** plus the standing multi-AI red-team pass.

**The user's argument, which the measurements support:** accommodating 43 incompatible idioms consumed a
large budget producing machinery to *discover* constraints we do not have to honour. The same pipeline
that produced correct scripts the first time can produce correct scripts again, against a cleaner spec.

---

## 2. WHAT THE MEASUREMENTS FOUND

Two read-only measurements, 2026-07-26. Full reports copied to
`notes/measure_register_sufficiency.md` and `notes/measure_rewrite_feasibility.md` (the originals live
under gitignored `_scratch/pass1_dim_survey/measure/`).

### 2.1 The parameter register CANNOT seed the module
The proposal was to seed the shared module from `notes/parameter_register.md`. It does not work:

| Measure | Value |
|---|---|
| Coverage of scripts' dimension-bearing quantities | **105 / 226 = 46.5 %** |
| `pending` rows | 2 (`W_slab`, `κ_4`) |
| Distinct bases in the corpus | **6 conventions over 4 axis sets** |
| Register-vs-locus agreement (21 stages sampled) | 67 agree / 4 disagree |

- **It is a *parameter* register (knobs), not a *quantity* register.** The 121 misses are a missing
  *category*, not missing stages: 74 intermediates/derived composites, 23 fields/coordinates, 12 knobs
  under a different key, 8 convention-variants, 4 benchmarks.
- ⛔ **The identity key is not machine-readable** — prose cells with strikethrough and escaped pipes, no
  `quantity_id`. `κ`, `C`, `A`, `K_η`, `T_Ω` each name several distinct quantities; `c_E`/`C_E` differ
  only by case. A module keyed on this cannot be built.
- ⛔ **Basis handling fails on both extreme stages.** stage042's stiffness-basis quantities are **absent
  entirely**, and **no normalisation convention is stated anywhere** — including
  `manifests/DIM_ORDER_DECISION.md`. stage038's four-axis quantities are present but normalised into
  `{L,T,M}` **with different values**: its `mu_R` contradicts the register row, its `A_E` contradicts
  stage037's.
- ⚠ **Some provenance claims are false.** Rows for `T_Ω`, `β₂`, `M̃/K̃/T̃_Ω` claim "stage 017,
  dual-engine verified", but **stage017 has no dimension machinery in either engine** — it only cites a
  flag. The values live at stage016. Same shape for `Vp0/ℓ_c`, which **no engine computes**.

**Verdict: use the register as a correctness ORACLE for the 105 quantities it covers. Not as the seed.**

### 2.2 There is no single idiom to consolidate — and the `.wl` control has a hole
| Measure | Value |
|---|---|
| Dimension machinery | **3,752 / 39,511 lines = 9.5 %**, in 30 of 43 scripts |
| Per script | min 0 / median 81 / max 351 (stage021, 45.5 % of file) |
| `.wl` coverage of dimensional content | COVERED 30 · PARTIAL 5 · UNCOVERED 8 |
| Distinct orderings | `(L,T,M)`×14, `(L,M,T)`×8, `(M,L,T)`×5, `(L,T)` 2-axis ×1 (008), 4-axis ×1 (038), stiffness ×1 (042) |
| Scripts with fractional exponents | **8** — four of them in `{L,T,M}`, so fractions are not 042-only |

⭐ **THE LOAD-BEARING FINDING.** The `.wl` is **not** a sufficient gate for *this particular* rewrite:

> The `.wl` only reveals a `.py` error when the **declared dimension** and its **assertion target**
> disagree. A shared module supplies **both from one place**, so a module bug moves them together —
> `.py` stays green and `.out` stays byte-identical.

Compounding it: **there is no cross-comparison harness at all** (the two runners write to separate trees
and never diff), and only **9 of 30** `.out` files print a computed exponent triple — the other 21
collapse the whole check to `PASS <TOOTH>`. And **22 of 30 stages need their literals transposed**
between orderings, where a transposition error is symmetric by construction and therefore invisible.

**The PASS gate is also weaker than assumed:** 16 scripts emit no `PASS tally:` at all (001–003,
016–028); the *ordered* marker list is **not** stable (stage007 swaps two adjacent lines between
baseline runs, so an ordered gate false-alarms while a multiset gate holds); stage013 emits 5 duplicate
markers and 042 one, so a set comparison silently loses 6; and only **2 of 30** dimension-bearing
scripts guard their own check arity, so a check made vacuous keeps its `PASS`.

---

## 3. ⭐ RED-TEAM FINDINGS ABOUT THE EXISTING LEDGER (independent of the refactor)

These are defects in committed work, surfaced by the measurements. **Catalogue, do not silently fix.**

1. **stage031: 20 of 21 `PASS_UNITS_*` rows compare an expression to a copy of itself.** They cannot
   fail.
2. **stage017 asserts a hardcoded `True`** in both engines.
3. **`K_eta` carries three different dimensions** across stages 013 / 016 / 023. The register documents
   the 013-vs-023 pair but not 016.
4. **Two stages compute the same ratio `r_BA` in incompatible unit systems.** Both land dimensionless,
   so nothing catches it.
5. **`ε0/ε1` vs `Z0_ret/Z1_ret`** — the register calls these the same two dofs; they carry different
   dimensions (`1` vs `M T⁻²`).
6. **Register locus mis-attribution** (§2.1) — dual-engine verification claimed for stages that compute
   nothing dimensional.
7. **`A_E`**: stage038 says `L·E`, stage037 says `M L³T⁻²`.
8. **`μ_η`/`T_w`** line-vs-volume convention clash; **`M0`** is `T⁻¹` at 009/010 but `M T⁻¹` at 023.

⚠ **Nine same-name/different-dimension pairs total.** These are the real work of the rewrite: they must
be **adjudicated**, not transcribed. Several are likely genuine physics/bookkeeping errors.

---

## 4. THE NEW PLAN, IN ORDER

1. ⭐ **Build the cross-engine dimension comparison harness.** Every dimension-bearing stage prints its
   computed exponent triple; the harness diffs `.py` against `.wl` **on values**, not on pass/fail. This
   is the gate the rewrite needs, it does not exist today, and it would have caught every §3 conflict on
   its own. **Independently valuable regardless of the rewrite.**
2. **Fix the register's false provenance claims** (§2.1) — a maintained document asserting
   dual-engine verification for stages that verify nothing is a defect that misled this measurement.
3. **Design the shared dimension module.** Source of truth = each stage's physics, cross-checked against
   the register where it speaks (105 quantities). Must support: ≥4 axes, fractional/symbolic exponents,
   and a stated basis convention. **Design against stage042 (stiffness basis, fractions) and stage038
   (four axes) first** — they are the extremes and they broke the register.
4. **Rewrite stage-by-stage**, adjudicating §3 conflicts as they surface. Escalate anything that changes
   a physical claim rather than a convention.
5. **Standing multi-AI red-team pass** over the result, per repo-root `docs/development_pipeline.md`.

**Gate for the rewrite** (replaces "identical PASS counts", which §2.2 shows is weak): the cross-engine
harness agreeing on values, plus a **multiset** (not ordered, not set) PASS-marker comparison, plus
`mathematica/out/` staying byte-identical after `$NNNNN` normalisation.

---

## 5. WHAT IS PARKED, NOT DELETED

The survey apparatus is **not** deleted — it is parked, and parts stay useful:
- `schemas/` (committed, `35dc6aa0`) — record schema, verification schema, 3,106-line validator, 173
  green fixtures. Could still validate rewrite records if we want structured per-stage output.
- `_scratch/pass1_dim_survey/directive/SURVEY_DIRECTIVE.md` at **r18** — gitignored; backed up at
  `_scratch_backups/`.
- **Repair round 8 fixed real defects** and is committed: `Dim()` was unrecordable in 13 scripts
  (fabrication-forcing), the `basis_locus` applicability contradiction, and two false claims about
  Python scoping. Its exploit-derived fixtures remain valid tests of that validator.

⛔ **Dead: the premise that the module must accommodate existing idioms.** The 4-script pilot
(`PILOT_DIRECTIVE.md`), the 44-way agent dispatch, and tasks #12/#16/#17/#18/#19 were all in service of
that premise. #17/#18/#19 only matter if the survey is revived.

**Do NOT redo:** the `.out` determinism probe (44/44 reproduce), the SymPy baseline (43/43 exit 0), or
the two measurements in §2 — all recorded.

---

## 6. WHAT THE MEASUREMENTS ALREADY GIVE THE MODULE DESIGN
The survey's actual deliverable was "what must the module support". Much of it now exists:
- **226 distinct dimension-bearing quantities** enumerated across the 43 scripts.
- **6 basis conventions over 4 axis sets**, with the per-stage assignment.
- **8 scripts using fractional exponents.**
- **The 30 stages carrying dimension machinery**, with line counts.
- **9 cross-stage conflicts** to adjudicate.

Method caveat from the report: the 226/401 counts are exact under its stated criterion; bucketing of the
121 register misses is a judgement call at the margins (±5 between "intermediate" and "renamed knob").
