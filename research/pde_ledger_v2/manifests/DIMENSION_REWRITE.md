# DIMENSION REWRITE — the single canonical doc for this workstream

**This file replaces `REWRITE_HANDOFF.md`, `PIVOT_TO_REWRITE.md`, `RESTART_PROMPT.md` and
`notes/dimension_inventory_and_provenance.md`.** Fold new findings in HERE. Do not write a new doc to
explain something this one already covers — that habit produced ten overlapping files for one
workstream and is the reason this consolidation exists.

Branch `ledger-v2-rebuild`. Last updated at `f5ff1843`.

---

## 0. READ ORDER — and the failure this prevents

1. `docs/model_map.md` — the model in one place. **Read it before touching a script.**
2. `STATUS.md` — the program's front door.
3. **This file.**
4. `notes/rewrite_reference_table.md` — the per-stage data (basis, loci, hazards).
5. Evidence, only if you need to check a number: `notes/measure_register_sufficiency.md`,
   `notes/measure_rewrite_feasibility.md` (⚠ has known errors, see §8), `notes/harness_infeasibility_review.md`.

⚠ A session resumed on a detailed handoff prompt, skipped `model_map.md`, and began re-deriving
conclusions the corpus already recorded. **A specific handoff prompt is exactly when re-reading feels
skippable, and exactly when it bites.**

## 1. ⭐ THE CHARTER — a USER decision, do not re-litigate

> *"The shared dimension module is **SymPy-ONLY**. Mathematica stays an **independent verification
> engine**. Consequence: the `.py`-vs-`.wl` dimension comparison is a **permanent standing
> cross-check**, not a one-time pre-unification sweep."*

Everything below follows from this. The `.wl` files never import the module, never receive values from
the `.py`, and are never "standardised" to agree with it. The comparison is infrastructure we keep,
not scaffolding we discard once the corpus is unified.

## 1b. ⭐ USER DECISIONS, 2026-07-26 — "correctness is king"

> *"I don't care if the scripts are byte identical on the output or if we need to do more derivations.
> Correctness is king. If we need to rename or whatever, let's make it happen. Just because something
> will be hard to do doesn't mean we shouldn't do it."*

**D1 — stdout byte-identity is a DIAGNOSTIC, not a blocking gate.** Report it when it holds (it is good
evidence a refactor changed representation and not results), but **never let it constrain a correct
change**. It only ever existed for stages 001–028 anyway (§4).

**D2 — print-only on the `.wl` is RELAXED.** A `.wl` may gain genuine new computation in order to
expose a value. ⇒ **The four "impossible" stages reopen** (037, 036, 035, 044 — they were blocked *only*
by print-only), as does any `UNREACHABLE` waived because reaching it needed a data-flow change
(stage011 `MassDim`/`OmegaDim`, stage012 `mass_dim`). Re-derive the "~26 well-gated + 4 that can't be"
assessment; it was made under a constraint that no longer applies.

**D3 — ⛔ ENGINE INDEPENDENCE IS *NOT* RELAXED.** This is not a difficulty constraint, so "correctness
is king" does not loosen it — it *is* the correctness. A derivation added to a `.wl` must be **its own
route**; deriving it by transliterating the `.py` makes agreement meaningless. See
`research/pde_ledger/notes/MATHEMATICA_MIRROR_POLICY.md` and the §1 charter. Adding computation to a
`.wl` therefore requires stating **how the route differs** from the `.py`'s.

**D4 — RENAMING IS AUTHORIZED, in one direction only.**
- ✅ Rename so **the same quantity shares one name** across stages (`EnergyDim` → `energy_dim`). This
  *creates* cross-stage comparisons that currently go undetected — see §8's GROUPING LIMITATIONS.
- ⛔ **Never** rename so that **different quantities share a name**, and never rename to make a
  conflict disappear. That destroys the check and hides reduction debt (§7). These are opposite
  operations; only the first is sanctioned.
- Fold the standardisation into each stage's rewrite as it comes up.

**D5 — cost and difficulty are not reasons to choose a weaker approach.** Where a correct route needs
more Wolfram runs, more derivations, or more review legs, take it.

## 2. WHY THIS EXISTS (not refactoring hygiene)

The 44-stage manifest fanout is **blocked**: dimension recovery covers only ~16 of 43 scripts while the
schema requires all. The decision (`b5527062`, `aae5d389`) was to **fix the corpus, not weaken the
check** — thirteen dimension idioms across 43 scripts *is* drift and the checker was reporting it
correctly. One shared module lets the checker delete `BARE_TUPLE_DIM_ORDER_BY_SHA256`, the AST
bare-tuple recovery, live module execution, and the per-script order registry: **one import, one
recovery path.**

Three consumers: the per-stage cross-check · the composite checker's recovery path · **Part VII's
whole-system dimensional firewall**, a named ledger deliverable.

Standing rule from the plan of record: **never adjust the process because the corpus is inconvenient.**

## 3. STATE

**4 of 30 done** — stage004, stage011, stage012, stage013. (43 audit scripts; 13 have no dimension
machinery: 001, 014, 015, 017, 019, 020, 022, 024, 025, 026, 028, 033, 043.)

| stage | compared | waived | detectors (L↔M / M↔T / L↔T) | note |
|---|---|---|---|---|
| **004** | 20 / 20 | none | 17 of 20 (L↔M) | `(L,T,M)`, `render="symbolic"` |
| **011** | 10 / 12 | `MassDim`, `OmegaDim` | 10 of 10 | ⚠ **waivers REOPENED by D2** |
| **012** | 18 / 19 | `mass_dim` | 16 / 8 / 16 of 18 | ⚠ **waiver REOPENED by D2** |
| **013** | **15 / 15** | **none** | 12 / 12 / 10 of 15 | first stage with zero waivers |

Canonical table (`CANONICAL_DIMENSIONS.md`, regenerated): 66 quantity rows, 2 candidate groups
(`KDim`, `Tw` — both AGREE), 0 `NEEDS_ADJUDICATION`.

| commit | what |
|---|---|
| `adcfbdfd` | module v0 + stage004 + stage011 |
| `60e7032c` | sidecar + comparator **with a non-empty floor** |
| `d29055b6` | stage011 `.wl` emission — cross-check 2 → 10 |
| `57a9bf34` / `f5ff1843` | stage012 `.wl` emission / `.py` rewrite |
| `6ad649cc` / `cc7f785e` | doc consolidation — 9 dead docs deleted, STATUS.md 491 → 107 |
| `ae7e4a7d` / `f01de097` | the generated canonical table (output-only) |
| `50972622` | ⭐ **the D1–D5 user decisions** (§1b) |
| `2a55a2bb` / `4391c69c` | stage013 `.wl` (incl. `K_eta`, first D2 case) / `.py` — 15/15, no waivers |
| `535e8e4b` | ⭐ anti-substitution rule → `docs/development_pipeline.md` checklist |
| `005c8f46` | ⭐ orchestrator `.out`-regeneration step (e2) + measured gate limits |

⚠ **Read a stage's `COVERAGE|` line before believing it is done.** "2 of 30" was once claimed while
stage011's gate passed with its own able-to-fail teeth outside the comparison.

## 3b. ⚠⚠ REOPENED BY D1/D2/D4 — conclusions that are NO LONGER TRUE

These were correct under constraints the user has since lifted (§1b). **They will read as settled
unless you check here first.** Do not inherit them.

1. **stage011's `MassDim`/`OmegaDim` and stage012's `mass_dim` waivers.** Waived because reaching them
   needed a data-flow change in the `.wl` — which **D2 now permits**. `CANONICAL_DIMENSIONS.md` still
   shows them `ONE_SIDED_PY (WAIVED)`. stage013 proves the pattern: `K_eta` was waived-equivalent under
   print-only and is now derived by both engines. **Go back and close these.**
2. **"037, 036, 035, 044 are `.wl`-emission-impossible."** They were blocked *only* by print-only.
   They need genuine, independently-routed computation — which is now allowed.
3. **"~26 well-gated + 4 that can't be."** That estimate assumed (2). Re-derive it.
4. **The 5 `GROUPING LIMITATIONS` in `CANONICAL_DIMENSIONS.md`** (`EnergyDim`/`energy_dim`,
   `CsSquaredDim`/`cs_squared_dim`, `CorruptKDim`/`corrupt_K_dim`, `FourVolumeDim`/`four_volume_dim`,
   `MassDim`/`mass_dim`) are closable by **D4** renaming — same quantity, one name. Currently only
   `KDim` and `Tw` group; standardising would group all six stage011/012 shared quantities.
5. **stdout byte-identity as a gate.** D1 makes it a diagnostic. Do not contort code to preserve it.

## 4. ⭐ THE PER-STAGE LOOP

**(a)** Add **print-only** `DIM|` output to the stage's `.wl`. **(b)** Re-run; confirm exit 0, PASS
tally unchanged, `.out` reproduces byte-identically after `sed -E 's/\$[0-9]+/$N/g'`; re-baseline.
**Commit this before touching the `.py`** — freezing the reference first makes independence structural,
not disciplinary. **(c)** Write the prediction down in `notes/stage0NN_rewrite_prediction.md`.
**(d)** Rewrite the `.py` onto `scripts/ledger_dimensions.py`. **(e)** Compare, axis-labelled.
**(f)** Review: adversarial fresh agent → Grok physical-correctness pass. **(g)** Commit.

⭐ **(e2) THE ORCHESTRATOR MUST REGENERATE THE `.out` ITSELF, once per stage.** Verification agents are
barred from Mathematica, so they **cannot** confirm the reference side is genuine — as one put it, *"a
hand-edit and a real re-run are byte-identical, so its provenance rests on trust alone."* The `.out` is
the reference half of the only universal gate; if it could be hand-written, the cross-check proves
nothing. Run `math -script <the .wl>`, normalise with `sed -E 's/\$[0-9]+/$N/g'`, and confirm it
reproduces the committed `.out` byte-for-byte. Done for stage013 (sha `42ee1ad7fbf8283a`, exit 0).

⭐ **(f2) The physics leg is NOT optional, and here is the measured reason.** On stage013, **5 of the 15
compared records** (`symbol_dims.*`) are declared as literals in *both* engines. Comparing them catches
a transcription divergence between the two — but **a SHARED wrong declaration passes the comparator**.
Only an independent read against the model catches that. stage013 scored 9 CORRECT / 0 WRONG on that
leg; do not skip it because the cross-engine numbers look clean.

**Artifacts.** `.py` → `scripts/<basename>.dimensions.txt`, same flat format the `.wl` prints:
`DIM|axes=L,M,T|name=EnergyDim|exponents={2, 1, -2}`. One parser, both engines.
⛔ Not stdout · ⛔ not JSON (these are read by humans and review agents).

**Three gates, and each proves something different:**
1. **PASS multiset identical** — the audit still reports what it reported.
2. **stdout byte-identity** vs `tail -n +7 scripts/output/<basename>.txt` (diff exactly 6 wrapper
   lines) — behaviour preservation. ⚠ **A DIAGNOSTIC, NOT A BLOCKING GATE** (D1, §1b): report it, never
   let it prevent a correct change. Also **NOT a transposition detector** (§6).
   ⛔ **This gate only exists for stages 001–028.** `scripts/output/*.txt` covers 001–028 ONLY —
   **stages 030–044 have no committed Python transcript**, so gate 2 is unavailable there and gates 1
   and 3 carry the whole load. Plan for that before starting group B/C/D; do not discover it mid-stage.
3. **axis-labelled cross-engine comparison** — the only universal transposition detector.

**⛔ The comparator must have a non-empty floor, and does.** The first version reported `PASS` with
`compared=0` (header-only sidecar / no name matches / `.out` stripped of `DIM` lines). It was live:
stage011 passed at `py=12 wl=2` with both `Corrupt*` teeth uncompared. Now `compared==0` fails, and
unwaived `py_only`/`wl_only` fails. Waivers are per-stage, name every quantity, and are echoed in a
`WAIVERS` line every run.

## 5. ⛔ THE FABRICATION GUARD (step a) — where this can silently destroy itself

⚠ **Amended by D2/D3 (§1b).** New computation in a `.wl` is now ALLOWED where it is needed to expose a
value — but it must be an **independent route**, and the directive must require Codex to state *how
that route differs from the `.py`'s*. Everything below still stands; only the blanket "no new
computation" clause is lifted.

- **Prefer print-only** where a live binding already exists — it is cheaper and carries no
  independence risk. Add computation only when the value is genuinely unreachable otherwise.
- ⛔ **Never hardcode an exponent literal in a `Print`.** A constant compared to the `.py`'s constant is
  vacuous. This defect already exists at **stage013 `.wl:446-447`** and **stage018 `.wl:386`** — do not
  add a third. Use `ToString[InputForm[dim["KDim"]]]`.
- ⛔ **Never copy a value or the axis order from the `.py`.** Read the `.wl`'s axis order from the `.wl`
  and cite the locus. *Names* are join keys and may be aligned; *values* may not.
- ⭐ **`UNREACHABLE` is a correct answer** — with a reason and locus. But see §9: it must be earned per
  case, never granted by a rationale carried over from another stage.
- ⛔ **Prove the emission is live**: perturb the expression, watch the printed value move, restore.
  Source-grepping for literals is dodgeable; a runtime move is not.

## 6. ⭐ WHAT THE PILOT PROVED (why step (e) is mandatory)

Transposing stage011's basis left it **passing all markers, exit 0, stdout byte-identical**. A
`render="tuple"` transcript emits exponents in declared-axis order, so relabelling is **unobservable**;
`render="symbolic"` detects it but has blind spots (4 of 9 lines at stage004). Once the module supplies
both the declared value and its assertion target, a script's own checks are blind.
Re-confirmed independently on stage012: under **all three** transpositions its 82 assertions still
pass; only the comparator sees it (L↔M 16/18, M↔T 8/18, L↔T 16/18).
Corollary: a sum-of-squares residual is permutation-invariant and can never detect this.

## 7. THE INVENTORY — ALREADY MEASURED, DO NOT RE-RUN

`notes/measure_register_sufficiency.md`: **401 `(stage,quantity)` pairs = 226 distinct quantities**
over the 30 scripts with dimension machinery; the register covers **105/226 = 46.5 %**; **6 bases over
4 axis sets**; no symbolic exponents anywhere; fractional in 8 scripts. Counts exact; bucket labels ±5.

**⭐ Most "conflicts" are REDUCTION LEVELS, not drift.** The model spans a **4D bulk and a 3D brane**,
so integrating out a coordinate shifts L-powers legitimately: `K_eta` is volume-density (L³, 016) vs
line-density (L¹, 013) vs reduced scalar (L⁰, 023), all integrating to `M T⁻²`. Verified across four
quantities (`K_eta`, `T_w`, `μ_η`, `T_Omega`), each off by exactly its measure's L-power.
⛔ **Never resolve one by renaming the variants apart** — it destroys the check *and* hides reduction
debt, shrinking stage043's irreducible count `[40,49]`. The dimension is the symptom; the unperformed
reduction is the thing (register `:170` counts it as debt).

**⛔ One genuine error, not a convention:** stage038 computes `r_BA=q_T²c_γ²/(μ_R A_E)` with
`μ_R=M²L T⁻⁴E⁻¹`, `A_E=L·E`; 036/037/040/044/003/007 use `μ_R=M L⁻¹T⁻²`, `A_E=M L³T⁻²`. Both land
dimensionless so nothing catches it. **stage038's 4th axis is not independent**: `E = M L²T⁻²`,
uniquely determined (fixed by `A_E`, independently confirmed by `μ_R`), stated by no source.
`A_E` has **no dimension row**; the register's "already-registered (R63)" points at the bc_selection
R1, not a dimension declaration — *covered-by-edge ≠ own row*.

**Provenance is archaeology, not opinion.** "Is `lambda_T` new or renamed?" cannot be answered from
memory — Codex wrote the math. Ask it **quote-backed**, with `NONE_FOUND` a first-class answer, then
verify. Sources: `notes/stages/`, `*_source_map.md`, `software/*/reports/`, `decisions/`, 517 commit
messages, the `_scratch/` dirs, `research/` prior papers.

## 8. COVERAGE CENSUS — the "group A is free" premise was FALSE

| stage | py dims | **real** compared | labelled? | blind |
|---|---|---|---|---|
| 012 | 19 | **18** ✅ done | DIM records | L↔M 16, M↔T 8, L↔T 16 detect |
| 013 | 14 | **9**/14 (5 literal-vs-literal) | mixed | 9 |
| 018 | 5 | **3**/5 (2 literal-vs-literal) | labelled | 4 |
| 016 | 21 | 9/21 | labelled | 4 |
| 023 | 29 | 7/29 | BARE | 6 of 7, **5 dimensionless** |
| 027 | 17 | **1**/17 | bare + hardcoded gloss | 0 |
| 021 | 18 | 7/18 | labelled | **7 of 7 — verdict quantity dimensionless** |

⛔ **021's gate cannot detect a transposition at all.** 027 is 1-of-17. 013/018 print hardcoded literal
strings in both engines. ⇒ every group-A stage needs `.wl` emission; with the floor in place they now
**fail** rather than pass quietly. The old ordering was a cost heuristic, never a correctness one.
⚠ 021 also has cross-engine name collisions (`[P₀_raw]`/`[P0_raw]`, unicode `N₀` vs ASCII `N0`).

**Order:** 013, 018, 016, 023, 027, then 021 (heaviest). Then `(L,T,M)`, `(M,L,T)`, then 008 (2-axis),
038 (4-axis), 042 (stiffness). 013 shares stage012's scaffolding.

✅ **038 and 042 are already CLEARED as non-blocking — do not re-litigate this.** Proven by execution
(`scripts/probe_ledger_dimensions_extremes.py`, committed): the module handles 038's 4 axes
`(M,L,T,E)` with a heterogeneous `None`-sentinel slot, and 042's `(stiffness,L,T)` with
`fractions.Fraction` exponents plus the `len(set(...))` homogeneity test that requires hashability.
Each capability is ablation-verified (`_exact`→`sp.Float` trips the exactness check; `__hash__`→`id`
trips the set check). They are ordered LAST for effort, not for risk.

## 9. ⛔ LANDMINES
- **stage003's `.wl` is `(M,L,T)` and neither file says so.** Emitting `axes=L,T,M` corrupts every triple.
- **stage042's `.wl:816` comment says "MLT" — a mislabel;** its axes are `(stiffness, L, T)`.
  ⚠ Its guard runs **once**, not twice (an earlier claim was wrong) — 042 is likely recoverable.
- ⚠ **037, 036, 035, 044 were listed as `.wl`-emission-impossible — REOPENED by D2 (§1b).** They were
  blocked *only* by the print-only rule, which is lifted. Re-assess each; they need genuine (and
  independently-routed) computation, not a workaround.
- `lru_cache` — **5 stages, 11 sites**: 018 (×1), 022 (×4), 023 (×1), 027 (×1), **040 (×4)**. 043 has
  **zero**; 040's four are verified argument-pure.
- `grep -c '^PASS'` **over-counts by exactly 1** (the tally line self-matches).
- Codex needs `--sandbox danger-full-access` for Mathematica · **≤2 concurrent kernels**.
- ⛔ **Never `git checkout`/`restore`/`stash` to undo an ablation** on uncommitted work — restores from
  HEAD and destroys it. Use `cp` backups; verify with `git hash-object`.
- ⛔ **Never state the reason you expect in a directive** — it gets adopted instead of checked. The
  stage012 waiver wrongly covered `omega_dim` because the directive supplied stage011's rationale.
  Ask for the determination **plus its evidence**, and point the review leg at every escape hatch.
- ⚠ Process probes lie: `pgrep -f 'math -script'` matched Codex's own argv (the directive quoted the
  string); `[ -d /proc/$PID ]` with an empty PID tests `/proc/` and always says "alive". **Check the
  artifact, not the process.**

## 10. ⛔ WHAT NOT TO DO
- Do **not** generalise the module beyond the stage in front of you. stage011's expression walker is
  deliberately stage-local.
- Do **not** build a corpus-wide inventory, oracle, or completeness proof. Three such gates were
  specified and rejected; the per-stage loop replaces them.
- Do **not** repair tautologies during a refactor (stage031's 20-of-21 self-comparing rows, stage017's
  hardcoded `True`, stage004 `:302`/`:308`, stage013/018's literal prints, stage012's stale
  `source dict L522-L529` label). Catalogue them — task #22.
  ⭐ **stage013's hardcoded prints are now DEMONSTRATED stale-capable, not merely suspected.** Under a
  `Tw` ablation, `.py:588-589` / `.wl:446-447` printed `Tw=(1,1,-2)` and `K_eta=…=(-1,1,-2)` while the
  live `k_shared` in the *same run* had moved to `(0,1,-5)` — the strings contradicted the run that
  produced them and **nothing flagged it**. That is the strongest evidence yet for task #22, and it
  generalises: any hardcoded rendering of a computed value can silently go stale.
- Do **not** centralise the *values*. A single canonical table is wanted, but as **generated output,
  never an import** — if scripts read their values from it, agreement becomes tautological. The layers
  are: shared **machinery** (imported) · per-stage **declarations** (not shared) · a **generated**
  table (output only) · `parameter_register.md` (what values *should* be, checkable against the table).
- ⚠ Do **not** call cross-stage agreement independent verification. stage011 and stage012 agree on six
  quantities but both apply the *same recipe* — real transcription fidelity, not two derivations
  converging (Grok, `f5ff1843`).

## 11. WHAT ACTUALLY VERIFIES A STAGE
Engine-vs-engine agreement is necessary and **not sufficient** — both engines can be wrong together.
The check that closes the gap is deriving the emitted dimensions from the *model* (`model_map.md` §2:
4D bulk, `[ρ]=L⁻⁴`, `P=Kρ⁵`, `c_s²=5Kρ⁴/m`). On stage012 that returned **14 CORRECT / 0 WRONG**,
including `[K]=[P]/[ρ]⁵=ML¹⁸T⁻²` reproducing the declared primitive, and confirming `pressure`
`{-2,1,-2}` is right *as a 4D pressure* (`ML⁻²T⁻²`), not the 3D `ML⁻¹T⁻²`. **Make this a standing leg.**

## 12. OPEN / DEFERRED
- **Task #13 — semantic dimension adjudication, DEFERRED NOT DROPPED.** Descoped after four review
  rounds: the verdict was always computed from fields the artifact itself supplies, i.e. a
  document-shape validator was being asked to do symbolic mathematics. The generated table (§10) is
  the right home for it.
- The `r_BA` unit-system adjudication (§7) — a model question for the user.
- 12 "registered under a different key, **or new**" quantities (§7) — quote-backed archaeology, then gauntlet.
- `schemas/` + `validate_dimension_survey.py` are **parked** (survey-era, still committed).
