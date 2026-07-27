# DIMENSION REWRITE — HANDOFF (2026-07-26)

⭐ **READ THIS FIRST, then `notes/rewrite_reference_table.md` (463 lines, the per-stage data).**
Background/why: `manifests/PIVOT_TO_REWRITE.md`. ⛔ Do NOT read the survey-era docs — that approach is dead.

⭐⭐ **BEFORE touching any script, read `docs/model_map.md` + `STATUS.md` + the plan-of-record commits
`aae5d389`/`b5527062`, then `notes/dimension_inventory_and_provenance.md`.** That last file is the
inventory + provenance map: **the 226-quantity / 401-pair inventory ALREADY EXISTS** in
`notes/measure_register_sufficiency.md` — do not re-run it — and §3 there explains why most
"conflicts" are legitimate **reduction levels** (4D bulk → 3D brane → line → volume → reduced scalar),
not drift. A session that skipped this began re-deriving recorded conclusions.

**State: 2 of 30 stages done, and now genuinely so.** `adcfbdfd` (module v0 + 004 + 011) →
`60e7032c` (sidecar + comparator with a non-empty floor) → `d29055b6` (stage011 `.wl` emission).
**28 remain.** Branch `ledger-v2-rebuild`.

⚠ **Earlier "2 of 30 done" was optimistic.** stage004 is complete (20/20 compared). stage011 was NOT —
its `.wl` exposed 2 of 12 quantities and the comparator passed anyway, with both `Corrupt*` able-to-fail
teeth uncompared. Fixed: 10/12 compared, 2 honestly waived, and a basis transposition now trips **10**
mismatches instead of 1. **Read a stage's `COVERAGE|` line before believing it is done.**

---

## 1. WHAT EXISTS

`scripts/ledger_dimensions.py` (179 lines) — axis-**labelled** dimensions, per-stage basis binding,
exact `sympy.Rational` exponents, multiply/divide/power. `DimensionBasis` is called to construct
positionally; zero-arg call = dimensionless. Equality and residual pair by **axis name**; `components()`
and rendering use the stage's declared order.

Done: **stage004** `(L,T,M)`, **stage011** `(L,M,T)`.

## 2. ⭐ THE PER-STAGE LOOP

**(a)** Add **print-only** output to the stage's `.wl` for dimension values it computes but doesn't
expose. Format used by the pilot: `DIMENSIONS|axes=...` then one `DIM|axes=...|name=...|exponents={...}`
per value. ⛔ No logic, no new computation, no control-flow change.
**(b)** Re-run the `.wl`; confirm it exits 0, the PASS tally is unchanged, and the output reproduces
byte-identically after `sed -E 's/\$[0-9]+/$N/g'`. Re-baseline the `.out`.
**(c)** Record pre-rewrite Python values and **write the prediction down** — what should match, what
shouldn't, and why. A mismatch predicted beforehand is evidence; one explained afterward is a
rationalisation.
**(d)** Rewrite the `.py` onto the module. Acceptance: exit 0, **identical PASS multiset**.
**(e)** Compare `.py` values against the `.wl` printed values, **axis-labelled, never positional**.
Any *unpredicted* mismatch stops the stage.
**(f)** Commit.

✅ **SETTLED + PROVEN END-TO-END on stage011 (`60e7032c`, `d29055b6`). Group A follows this exactly.**

**The artifacts.** The `.py` emits `scripts/<basename>.dimensions.txt`, a **sidecar** in the *same flat
line format the `.wl` prints into its `.out`*: `DIM|axes=L,M,T|name=EnergyDim|exponents={2, 1, -2}`.
One parser (`scripts/compare_dimension_artifacts.py`) serves both engines. ⛔ **NOT stdout** — stdout
must stay byte-identical to the committed transcript (see below). ⛔ **NOT JSON** — these are read by
humans and review agents.

**Two independent gates, and you need both:**
1. **stdout byte-identity** — live stdout vs `tail -n +7 scripts/output/<basename>.txt` must diff to
   **exactly 6 lines** (wrapper only: `#`+blank at top, blank+`EXIT_CODE: 0` at bottom). This is a
   *behaviour-preservation* gate: it proves the rewrite changed representation, not results. Holds on
   004 and 011 across a 114-line rewrite. ⚠ It is **NOT** a transposition detector — see §3.
2. **the axis-labelled cross-engine comparison** — the only universal transposition detector.

**⛔ The comparator MUST have a non-empty floor, and it now does.** The first version reported
`status=PASS` with `compared=0` — reachable via a header-only sidecar, names that all miss, or an
`.out` stripped of `DIM` lines. It was **live**: stage011 passed at `py=12 wl=2` with 10 quantities
uncompared *including both `Corrupt*` able-to-fail teeth*. Now `compared==0` fails, and any unwaived
`py_only`/`wl_only` fails. Waivers are per-stage, name every waived quantity, and are **echoed in a
`WAIVERS` line every run** so a partially-compared stage can never read as complete.
⭐ **Always read the `COVERAGE|` line.** Cross-engine coverage varies wildly per stage — see §6.

**⛔ THE FABRICATION GUARD for step (a) — this is where the process can silently destroy itself:**
- **Print-only.** Print values the `.wl` *already computes*. No new computation, symbol, assignment or
  control flow.
- ⛔ **Never hardcode an exponent literal in the `Print`.** A constant compared to the `.py`'s constant
  is vacuous. This defect already exists at **stage013 `.wl:446-447`** and **stage018 `.wl:386`** —
  do not add a third. Every printed exponent must read a live computed expression
  (stage011's pattern: `ToString[InputForm[dim["KDim"]]]`).
- ⛔ **Never copy the value, name, or axis order from the `.py`.** Read the `.wl`'s axis order from the
  `.wl` and cite the locus. Copying it turns the independent engine into a mirror and deletes the gate.
- ⭐ **`UNREACHABLE` is a correct, expected answer.** If the `.wl` genuinely does not compute a
  quantity, report it with a reason + locus and waive it. A short honest list beats fabricated records.
  (stage011: `MassDim`/`OmegaDim` stay local to `buildDimensionalBlock[]`; exposing them would need a
  non-print data-flow change, so they are waived.)

**⛔ RESTORE RULE for the step-(e) ablation.** The files are often **uncommitted** mid-loop. Never use
`git checkout` / `git restore` / `git stash` — they restore from HEAD and **destroy the uncommitted
work** (this happened, 2026-07-26). Use `cp` backups and verify with `git hash-object` before/after.

## 3. ⭐ THE RULE THE PILOT ESTABLISHED
Transposing stage011's basis `(L,M,T)`→`(L,T,M)` left the script **passing all 60 markers, exit 0**.
Its own assertions cannot see a transposition once both ends come from the module. The labelled
cross-engine comparison caught 7 of 10 values. **Therefore step (e) is mandatory wherever basis order
can change — it is the only thing between us and a silent corpus-wide transposition.**
Corollary: stage004's `dim_residual` is a **sum of squares → permutation-invariant**. Never use it as
the axis check.

## 4. GOTCHAS THAT WOULD BURN A SESSION
- **`grep -c '^PASS'` over-counts by exactly 1** (the `PASS tally:` line matches its own pattern).
  Verified: `real == tally` for all 27 tally-emitting stages, no exceptions. **16 stages emit no tally
  line at all**: 001–003, 016–028.
- **Codex needs `--sandbox danger-full-access` for Mathematica.** `workspace-write` fails with a
  spurious licence error.
- ⚠ **Licence cap: never more than 2 concurrent `math -script`.**
- **Never wait on `pgrep` for a pattern your own waiter contains** — it self-matches and the loop never
  exits. Cost 6h43m on 2026-07-26. Wait on the captured PID from a pidfile, or artifact + log quiescence.
- `scripts/output/*.txt` covers **only 001–028**; stages 030–044 have no committed Python transcript.

## 5. ⛔ LANDMINES — read before touching these stages
- **stage003's `.wl` is `(M,L,T)`, not `(L,T,M)`, and neither file states it.** Emitting `axes=L,T,M`
  there **silently corrupts every triple**.
- **stage042's `.wl` comment (`:816`) says "MLT"** — a mislabel; the basis carries `{1/2, 3/2, -1}`.
- **stage021 is worse than previously recorded**: one `dim_record` emits **three renderings under two
  conventions** (`:213`, `:227`, `:238`, `:257`). Not merely two orders in one file.
- **Five stages where step (a) is genuinely impossible, not just costly:**
  - **037** — `.wl` contains zero integer tuples; a `DIM|` line would be a **new derivation**.
  - **036** — dimensions exist only as literal arguments; printing re-transcribes constants.
  - **035** — monomials needing `Exponent[]` conversion.
  - **044** — only **1 of ~96** vectors reachable.
  - **042** — ⚠ **CORRECTED 2026-07-26, half of this claim was wrong.** The base map IS a
    `Module`/function local (`.py:830`, `.wl:833`) — but the guard does **NOT** run twice. There is
    exactly **one** call site per engine (`.py:1868`, `.wl:1845`), invoked once inside
    `run_assertions`; the mutation flag only toggles `cE` *within* that single invocation. An in-body
    print emits **one clean copy per process**, not corrupted duplicates. 042 is therefore likely
    RECOVERABLE, not impossible. Re-verify before relying on it.
- **`lru_cache` — the earlier note was half wrong.** 043 has **zero**. 040's four are verified
  argument-pure. **Previously missed: 018, 022 (×4), 023, 027** — 5 stages, 11 sites total.
- **`ACTIVE_MUTATION` at import time in 15 stages**; 030 (`:107`) and 031 (`:94`) read it **late**.

## 6. REWRITE ORDER — ⚠ THE "GROUP A IS FREE" PREMISE IS FALSE (measured 2026-07-26)
The `(L,M,T)` stages *do* render some computed dimension values, so the ordering claimed they need no
`.wl` edit. **They render far too few, and some render constants.** Measured real coverage:

| stage | py dims | **real** compared | labelled? | blind to a transposition |
|---|---|---|---|---|
| 012 | 14 | 12/14 | **BARE** tuples | 8 of 12 |
| 013 | 14 | **9**/14 (5 are literal-vs-literal) | mixed | 9 |
| 018 | 5 | **3**/5 (2 are literal-vs-literal) | labelled | 4 |
| 016 | 21 | 9/21 | labelled | 4 |
| 023 | 29 | 7/29 | **BARE** (order only in a section heading) | 6 of 7, **5 dimensionless** |
| 027 | 17 | **1**/17 | bare + hardcoded gloss | 0 |
| 021 | 18 | 7/18 | labelled | **7 of 7 — verdict quantity dimensionless** |

⛔ **stage021's cross-engine gate cannot detect a transposition at all** (every printed value is blind;
the verdict-carrying `[P₀^phys]=1` is dimensionless). **027 is 1-of-17; 023 is effectively 1** real
detector. **013/018 print hardcoded literal dimension strings in BOTH engines** (`.wl:446-447`,
`.wl:386`) matched by identical `.py` literals — vacuous, and they inflate apparent coverage to 14/14
and 5/5. See §10 and task #22.

⇒ **Group A needs `.wl` `DIM|` emission too** (step (a)), a Mathematica seat and an `.out` re-baseline
per stage. With the comparator's non-empty floor these stages now **fail** rather than pass quietly, so
this is not optional. The old ordering was a *cost* heuristic, never a correctness one.
⚠ **stage021 also has real cross-engine name collisions** — `[P₀_raw]`/`[P0_raw]`,
`[P₀^phys]`/`[P0_physical]`/`P0Physical`/`p0Dim`, `mu_hat0`/`muHat0`, unicode `N₀` vs ASCII `N0` —
which surface as one-sided quantities and need a name mapping.

Still true: 012/013 share stage011's `build_dimensional_block` scaffolding — near-mechanical repeats.

**Group A** `(L,M,T)`: 012, 013, 018, 016, 023, 027, then **021 at the end of A** (heaviest, three
renderings). **Group B** `(L,T,M)`: 005 and the rest. **Group C** `(M,L,T)`: 003 (⛔ see landmine), 010,
039, 040, 041. **Group D** the awkward bases: 008 (2-axis `(L,T)`, no mass axis), 038, 042.
See `notes/rewrite_reference_table.md` for the exact ordering and per-stage reasoning.

## 7. ⭐ OPEN DECISION — settle on paper BEFORE group A
`PIVOT_TO_REWRITE.md` §4 said design the module against **042 and 038 first** (the extremes). The pilot
did the opposite and the module **explicitly declined four-axis and stiffness semantics**.
**Assessment:** the structural risk looks low — the module stores an arbitrary-length labelled mapping
and its contract probe already passed a four-axis construction. What was declined is *semantics*, not
structure. **So: write the 038 and 042 basis declarations on paper and confirm nothing in the module
blocks them.** That is a short check, not a redesign. Doing it costs little; discovering a block at
group D costs a re-run of A–C.
**Good news on 038:** its `unitState` returns all 8 four-axis vectors as a live top-level `List`;
`out:31` omits `mu_R` only because a **hardcoded string** is printed instead of the computed value. It
is easier than previously assessed.

## 8. THE CONFLICTS — ⚠ RE-FRAMED 2026-07-26, most are NOT conflicts
`PIVOT_TO_REWRITE.md` §3 lists nine; the register measurement adds **`μ_R`, `A_V`, `T_Omega`, `Q_E`**.
**But the flagship turned out not to be a conflict at all.**

⭐ **`K_eta` across 013/016/023 is a REDUCTION LEVEL, not drift** — line-density (`dw`~L¹) vs
volume-density (`a²dw dΩ`~L³) vs already-reduced scalar (L⁰); all three integrate to `M T⁻²`.
Confirmed independently across **four** quantities (`K_eta`, `T_w`, `μ_η`, `T_Omega`), each differing
by exactly the L-power its measure predicts, and already recorded at
`measure_register_sufficiency.md` §6.5 and in-corpus at `…stage016…py:837`. The model spans a **4D
bulk and a 3D brane**, so integrating out a coordinate legitimately shifts L-powers.
⇒ **Expect most of the remaining "conflicts" to be the same shape. Resolve them per-stage as they
arise; do NOT re-audit all 13 up front.** Full detail: `notes/dimension_inventory_and_provenance.md` §3.

⛔ **Never resolve one by renaming the variants apart.** Three names for three reduction levels makes
the conflict vanish *and destroys the consistency check* — worse, it hides reduction debt and silently
shrinks stage043's irreducible count `[40,49]`. The differing dimension is the SYMPTOM; the
unperformed reduction is the THING, and the register already counts it as debt (`:170`).

⛔ **ONE genuine error, NOT a convention — adjudicate it separately:** stage038 computes
`r_BA = q_T²c_γ²/(μ_R A_E)` with `μ_R = M²L T⁻⁴E⁻¹`, `A_E = L·E`, while 036/037/040/044/003/007 use
`μ_R = M L⁻¹T⁻²`, `A_E = M L³T⁻²`. **Two incompatible unit systems, and both land dimensionless so
nothing catches it.** See `dimension_inventory_and_provenance.md` §4.1.
⚠ When one surfaces: **stop and adjudicate**. It is a question about the model, not the code — bring the
values and both loci to the user. **Do not "standardize" the `.wl` to agree with the `.py`** — that
converts an independent engine into a mirror and destroys the only cross-check we have. The `.wl`
changes **only** when adjudication rules its value wrong, and then it is a physics fix.

## 9. CORRECTIONS TO CLAIMS IN OLDER DOCS
- "30 stages remain" → **28** (004, 011 done). Corpus dim lines 3,752 → **3,675**.
- "four fractional stages within `{L,T,M}`" → **wrong**; the four fractional-exponent stages in question
  are `(L,M,T)`: 012, 021, 023, 027. (Eight stages have fractional exponents: 004, 005, 006, 012, 021,
  023, 027, 042.)
- The `lru_cache` note (§5) and the conflict count (§8) as above.

## 10. WHAT NOT TO DO
- ⛔ Do not generalise the module beyond what a stage in front of you needs. stage011's SymPy expression
  walker was deliberately left **stage-local** — two stages is not evidence for a universal walker.
- ⛔ Do not repair the catalogued tautologies while refactoring (stage004 `:302`, `:308`; stage031's
  20-of-21 self-comparing rows; stage017's hardcoded `True`). Task #22. Fixing them inside a rewrite
  hides pre-existing defects.
- ⛔ Do not build a corpus-wide inventory, oracle, or completeness proof. Three such gates were
  specified and rejected; the per-stage loop replaces all of them.
