# DIMENSION REWRITE — HANDOFF (2026-07-26)

⭐ **READ THIS FIRST, then `notes/rewrite_reference_table.md` (463 lines, the per-stage data).**
Background/why: `manifests/PIVOT_TO_REWRITE.md`. ⛔ Do NOT read the survey-era docs — that approach is dead.

⭐⭐ **BEFORE touching any script, read `docs/model_map.md` + `STATUS.md` + the plan-of-record commits
`aae5d389`/`b5527062`, then `notes/dimension_inventory_and_provenance.md`.** That last file is the
inventory + provenance map: **the 226-quantity / 401-pair inventory ALREADY EXISTS** in
`notes/measure_register_sufficiency.md` — do not re-run it — and §3 there explains why most
"conflicts" are legitimate **reduction levels** (4D bulk → 3D brane → line → volume → reduced scalar),
not drift. A session that skipped this began re-deriving recorded conclusions.

**State: 2 of 30 stages done.** `adcfbdfd` = module v0 + stage004 + stage011, verified.
**28 remain.** Branch `ledger-v2-rebuild`.

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

⚠ **Step (e) has no `.py`-transcript side — including for the pilot.** stage004's `.py` emits **zero**
exponent triples; the pilot's "29 matches" came from an **external in-process comparison script**, not a
transcript diff. Either keep using an external comparator or decide to have the `.py` emit triples too —
**this is unresolved and must be settled early**, because it decides whether the cross-check is
reproducible from committed artifacts.

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

## 6. REWRITE ORDER (from the reference table, and the reason matters)
⭐ **`(L,M,T)` group first — because the 9 stages whose `.out` already renders computed dimension values
are exactly that group plus 004**: 004, 011, 012, 013, 016, 018, 021, 023, 027. Every `(L,T,M)` and
`(M,L,T)` stage needs a `.wl` edit + `.out` re-baseline + a Mathematica seat; these do not.
Also 012/013 share stage011's exact `build_dimensional_block` scaffolding — near-mechanical repeats.

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
