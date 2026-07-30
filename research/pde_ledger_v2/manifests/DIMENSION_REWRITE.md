# DIMENSION REWRITE — the single canonical doc for this workstream

**This file replaced four earlier docs, all now DELETED (recoverable via git):
`REWRITE_HANDOFF.md`, `PIVOT_TO_REWRITE.md`, `RESTART_PROMPT.md`,
`notes/dimension_inventory_and_provenance.md`. Those names are history, not pointers.** Fold new findings in HERE. Do not write a new doc to
explain something this one already covers — that habit produced ten overlapping files for one
workstream and is the reason this consolidation exists.

Branch `ledger-v2-rebuild`. Documentation state committed 2026-07-27 (`35710cee`); the
stage-conversion baseline before that fold was `1b645ed9`.

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

⛔⛔ **WHAT "INDEPENDENT ENGINE" RESTS ON — verified 2026-07-27, and it is NOT provenance.** For
stage016 both engine files were added in the **same commit** (`bfac580f`), and a history audit at
`da75c58f` found **all 43 current dual-engine pairs were first added pairwise in one commit**
(`SAME_FIRST_ADD=43`, `DIFFERENT_FIRST_ADD=0`). **Git therefore establishes no derivation order for
any of them: it cannot show that one engine was not written from the other.** Independence rests entirely on
**authoring discipline** — D3 below and `research/pde_ledger/notes/MATHEMATICA_MIRROR_POLICY.md` — and
on review legs that screen for transliteration. ⇒ State the claim precisely: a green comparator shows
that **two implementations agree**, not that they were **reached by two independent routes**. Neither
the digests, the comparator, nor any per-stage gate can upgrade that. Say "agreement", not
"independent verification", unless a specific leg established the route. ⚠ Precisely: neither the
current digests, the comparator, nor any current per-stage value gate establishes independent
authorship — a *future* gate is not ruled out, it just does not exist.

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

**7 of 30 converted** — stage004, stage011, stage012, stage013, stage016, stage018, stage023. ⚠ *Converted* is not *finished*: **all seven are waiver-free** (011/012's three reopened items CLOSED `8b006055`), but see §11 for what a green comparator does and does not establish. (43 audit scripts; 13 have no
dimension machinery: 001, 014, 015, 017, 019, 020, 022, 024, 025, 026, 028, 033, 043.)

| stage | compared | waived | detectors (L↔M / M↔T / L↔T) | note |
|---|---|---|---|---|
| **004** | 20 / 20 | none | 17 of 20 (L↔M) | `(L,T,M)`, `render="symbolic"` |
| **011** | **12 / 12** | **none** | 11 / 9 / 10 of 12 (measured) | ✅ waivers CLOSED `8b006055`; first `.wl` emitting every named binding in its block |
| **012** | **19 / 19** | **none** | 17 / 9 / 16 of 19 (measured) | ✅ waiver CLOSED `8b006055` |
| **013** | **15 / 15** | **none** | 12 / 12 / 10 of 15 | first stage with zero waivers |
| **016** | **21 / 21** | **none** | 18 / 15 / 16 of 21 (measured, adversarial leg) | ⭐ largest record set yet; **12 of 21 are declared literals in BOTH engines** and 0 are computed from a physical input |
| **018** | **6 / 6** | **none** | 3 / 5 / 6 of 6 (measured) | emits **6 of its 10** objects — the 4 omitted are enumerated in the stage note (§4-a) |
| **023** | **29 / 29** | **none** | ⛔ **not measured** — the 023 ablations measure decoy-declaration and record-repoint detection, a *different* measurement from this column's transposition counts; the cell is left empty rather than filled with it | ⭐ **29 of 29 are hand-typed literals in BOTH engines** — 22 declarations + 7 targets typed on both sides; the 7 computed records *are* live `dimOf` walks, but over exactly those literals and no other dimensional input (stage note §1.7(3)) — the least independent dimension block of the four for which that count exists |

⚠ **018's L↔M rate of 3/6 is set by PHYSICS, not by an omission** — *every* stage018 exponent has
`M = 0` (`[a]=L`, `[c_s]=L T⁻¹`, coefficients `T^n`). It declares a 3-axis basis and populates 2; the
M row is `0 == 0` six times and carries no information. Honest, and worth expecting again wherever a
slice has no mass content.

Canonical table (`research/pde_ledger_v2/CANONICAL_DIMENSIONS.md`, regenerated): **122** quantity rows,
**6** candidate groups, **3 `NEEDS_ADJUDICATION`** (`KEta`, `TOmega`, `muEta`), **0 `ONE_SIDED_PY`**
(every one of the 122 rows is `AGREE`).
⚠ **The adjudication count rose because stage023 landed, not because drift did.** `KEta` now carries
three reduction levels (013 line `M L⁻¹T⁻²`, 016 volume `M L⁻³T⁻²`, 023 reduced scalar `M T⁻²`) and
`TOmega` is the new group on exactly that pattern (016 volume `M L⁻³T⁻²` vs 023 reduced scalar
`M T⁻²`) — the **correct** output of surfacing reduction levels that were never renamed apart (§7),
adjudicated as such in the stage023 note §1.7(2b) (*"the names stand"*). `muEta` remains the
013-vs-016 line-vs-volume debt.

⭐⭐ **`ARTIFACT_NAME_WAIVERS` IS NOW EMPTY (`{}`).** Every converted stage compares every name it
emits, with no exemptions. ⚠ **That is a coverage statement, not a strength statement** — see §3b.1.

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
| `005c8f46` | ⭐ orchestrator `.out`-regeneration step (now §4-g2) + measured gate limits |
| `63dee5e4` | stage018 `.wl` — ⭐ the **false cross-stage merge** caught (`c_s0` vs `c_S`) |
| `1b645ed9` | stage018 `.py` — **5 of 30**; the R3 file-I/O prose fix; the **stale-sidecar** hole found |
| `9dd55c8a` | group-A reachability surveys + the corpus-wide terminal-`Exit[]` correction |
| `5b29f400` | stage016 `.wl` — 21 records; the **axis-label/vector decoupling** caught and fixed |
| `4dfadaa4` | ⭐ **MEASURED** fix to §4 gate 2's stdout recipe (`tail -n +9 \| head -n -2`) |
| `8c4b25b0` | stage016 `.py` — **6 of 30**; first `NEEDS_ADJUDICATION`; sidecar-forgery + module self-attestation found |

⚠ **`ARTIFACT_NAME_SET|` reports emitted-artifact name symmetry, not source coverage.** A quantity
omitted from both artifacts is invisible to it; stage011's gate once passed with able-to-fail teeth
outside the comparison.

## 3b. ⚠⚠ REOPENED BY D1/D2/D4 — conclusions that are NO LONGER TRUE

These were correct under constraints the user has since lifted (§1b). **They will read as settled
unless you check here first.** Do not inherit them.

1. ✅ **CLOSED `8b006055` — stage011's `MassDim`/`OmegaDim` and stage012's `mass_dim`.** All three are
   now emitted by both engines; the registry is empty. ⚠ **Read the closure honestly:** stage011's
   `OmegaDim` is consumed by *nothing* in stage011 (`dimOf` runs only on `5 K rhoStar^4/m`, which has
   no `ω`) and groups with nothing cross-stage until the D4 rename — the waiver is closed, the
   **quantity is not thereby checked**. All three are primitive declarations in both engines, so the
   comparison catches transcription divergence, not wrong upstream physics. Historical detail below.
   ~~**stage011's `MassDim`/`OmegaDim` and stage012's `mass_dim` waivers.**~~ Waived because reaching them
   needed a data-flow change in the `.wl` — which **D2 now permits**. `CANONICAL_DIMENSIONS.md` still
   shows them `ONE_SIDED_PY (WAIVED)`. stage013 proves the pattern: `K_eta` was waived-equivalent under
   print-only and is now derived by both engines. **Go back and close these.**
2. ~~**"037, 036, 035, 044 are `.wl`-emission-impossible."**~~ ⚠ **FALSIFIED for 037 by a working
   prototype; 035/036 have identified routes that are NOT yet prototyped** (2026-07-27) — see the
   spike result at the end of this section, and mind the distinction. 044 untested (frozen pending
   044-v2).
3. **"~26 well-gated + 4 that can't be."** That estimate assumed (2), which is now falsified.
   ⭐ **stage018 is DONE (`1b645ed9`); the three reopened waivers are CLOSED (`8b006055`).**
   ⭐ **016 and 023 are DONE as well** (both engines, comparator green). Next conversion — which is not
   the build queue (§12b): **027**, then 021 (heaviest) — per §8's recorded conversion-order line. The
   **build queue** is a different sequence: **ablation driver first, `DIM|` emitter second, both before
   027 begins.**

⭐⭐ **THE 037 SPIKE, 2026-07-27 — `ROUTE_EXISTS` for stage037, demonstrated with a working
prototype.** Run out of order, deliberately: the seven exceptional tail stages are **035, 036, 037,
044, 008, 038, 042**, which under the prior remaining order sat at positions **14–16 and 22–25** — so
an infeasibility would have surfaced only well into the run, after most of the cheap stages were spent.

⛔ **READ THE SCOPE PRECISELY — this is where an overstatement was caught and corrected.**
**Only 037 was prototyped.** For **035 and 036** the evidence is *source inspection and transfer
analysis*, not execution. Do not cite them as measured.

- **The route:** symbolically rescale the *live* expressions under the file's own native `unitRules`,
  cancel to a scale monomial — `Factor[Together[(expr /. unitRules[False])/expr]]` — then `Exponent`
  against the axis map. **21 of 21 quantities reachable, zero waivers.**
- **Independence (D3) satisfied, by three signs:** a different **algorithm** (Wolfram
  rescales-and-cancels; the `.py` walks the expression tree with an atom-dimension dictionary); a
  different **axis order chosen independently** (Wolfram `M,L,T` vs Python `(L,T,M)` — a live proof
  that axis-*labelled* pairing is doing real work); and a **seeded `.py`-only error caught by the real
  comparator** (`MISMATCH U_B` → `status=FAIL`), the probe loading the actual
  `compare_dimension_artifacts.py` via `importlib`, not a reimplementation.
- **MEASURED on 037:** ~27–33 `.wl` lines · ~5 s runtime · 16/16 PASS preserved · `.out` re-run
  **required**.
- **ESTIMATED, not measured:** **0.5–1 engineer-day per stage including review.** ⚠ Review effort was
  not measured for any stage, and not at all for 035/036. Treat it as a planning figure.
- **Review classification (a judgement, and a firm one):** **new mathematics, not print-only** — a
  reviewer must check scale-factor cancellation, Laurent exponent extraction, axis labelling, and the
  nonzero/homogeneity assumptions.
- **Transfer — ANALYSED, NOT PROTOTYPED.** **035 looks directly transferable and cheaper**
  (`unitScalingObject[]` already computes a nine-key association ≈ 20 component quantities; the work
  is flattening + name alignment). **036 looks transferable but dearer** — `dimensionResiduals[]`
  discards identities and returns only zero residuals, so D2 restructuring must build a named
  association over 13 groups / 29 component quantities. ⛔ The old *"printing them means
  re-transcribing constants"* verdict is **avoidable** — 036 has native primitive scaling rules
  (`.wl:232`) and live kernel/interaction/force expressions (`.wl:158`) to derive from — but that is
  an inspection finding. **Expect surprises on first contact with either.**
- ⚠ **The honest limit, from the spike itself:** 8 of 037's 21 records are effectively primitive
  declarations — comparing them detects **transcription divergence, not wrong upstream physics**; the
  other 13 are algebraically derived. **Independent corroboration that the physics leg (§4-c) is
  mandatory.**
- ⭐ **BONUS — the spike's emission pattern is better than the one in use.** It derives the `axes=`
  label from `Keys[dimensionAxes]` instead of typing it as a string literal. **By inspection, all five
  then-converted `.wl` files hardcoded their axes strings** — as written at the spike, 2026-07-27
  (`35710cee`), when the converted set was 004/011/012/013/018. ⚠ **That count is SUPERSEDED:**
  ⭐ **011 and 012 NOW USE the derived machine label** (`dimensionAxesLabel[]`, adopted during the
  waiver closure `8b006055` — their old hardcoded loci no longer exist, so cite the commit, not a line);
  **004, 013 and 018 still hardcode theirs** (`004 .wl:220`, `013 .wl:447`, `018 .wl:387-393`). 011/012
  also retain *typed
  human prose* labels (`.wl:412`, `.wl:571`). The stage018 adversarial leg demonstrated the
  stale-label *risk* **for stage018 only** — it did not audit the others. If a `.wl`'s internal order ever changed, a typed label
  would not follow. **Adopt the derived-label form from 016 onward**; the debt is **three** stages, not
  five — **004, 013 and 018** — bounded and recorded.

⭐ **Evidence: `research/pde_ledger_v2/notes/stage037_dimension_emission_spike.md`** (tracked) —
the route, the measured numbers, the seeded-error result, and the measured-vs-estimated split.
The working prototype itself remains in gitignored `_scratch/spike037/`, so the note reproduces
the emission block in full.
4. **The EIGHT `GROUPING LIMITATIONS` in `CANONICAL_DIMENSIONS.md`** are closable by **D4**
   renaming — same quantity, one name. The pairs (stage011 CamelCase vs stage012 snake_case):
   `CsSquaredDim`/`clean_walk.cs_squared_dim` · `CorruptKDim`/`corrupt_K_dim` ·
   `EnergyDim`/`energy_dim` · `FourVolumeDim`/`four_volume_dim` · `MassDim`/`mass_dim` ·
   `OmegaDim`/`omega_dim` · `PressureDim`/`pressure_dim` · `RhoDim`/`rho_dim`.
   Only `KDim` and `Tw` group today; standardising takes stage011/012 to **9** grouped shared
   quantities. ⭐ **Target spelling = snake_case** (stage012/013's convention, matching the `.py`
   variable names) ⇒ **stage011 is the file to rename.**
5. **stdout byte-identity as a gate.** D1 makes it a diagnostic. Do not contort code to preserve it.

## 4. ⭐ THE PER-STAGE LOOP

⛔ **THE STEPS ARE IN EXECUTION ORDER. Run them in this order; the numbering is not decorative.**

⛔⛔ **BEFORE YOU START: NINE PATHS THIS LOOP TOUCHES ARE NOW UNDER A SECOND FREEZE AUTHORITY, AND IT IS
NOT THIS WORKSTREAM'S.** As of 2026-07-29 the ablation driver's committed conformance suite
(`notes/ablation_driver/fixtures_v4/`) is frozen against
`notes/ablation_driver/fixtures_v4.accepted.sha256`, and that authority's 36 governed paths include nine
live dimension-rewrite paths — `scripts/ledger_dimensions.py`, `scripts/ledger_dimensions.accepted.sha256`,
`scripts/check_ledger_dimensions_pin.py`, `scripts/compare_dimension_artifacts.py`, stage023's audit `.py`
with its `.dimensions.txt` sidecar and its Mathematica `.out`, and
`notes/stage023_step_h_evidence/include_list.tsv` + `results.tsv` — because acceptance item A7 re-runs
stage023's producer. ⇒ **Any conversion that changes one of them invalidates the fixture freeze**, and
`fixtures_v4/run_conformance.py` then refuses to run at all until it is re-greened.
- **Where this loop hits it, concretely:** step (g)'s `check_ledger_dimensions_pin.py --accept` rewrites
  `ledger_dimensions.accepted.sha256`, which is itself governed; any edit to the shared module, the pin
  script or the comparator changes a governed path; and stage023's sidecar is governed, so a re-run that
  changes its bytes invalidates the freeze — a byte-identical regeneration does not, since the authority
  compares digests. The ordinary case — converting a *different* stage's own `.py`/`.wl` — touches nothing
  governed.
- **What to do about it — `verify_freeze.py` has NO `--accept` mode, so the authority is regenerated by
  hand and then committed.** The file is a plain digest manifest; the verifier rebuilds the whole
  inventory each run, so re-greening means making the manifest match the new bytes exactly:
  1. **Inventory** (`fixtures_v4/verify_freeze.py:54-61`) = every **file or symlink** under
     `fixtures_v4/` (`SUITE.rglob("*")`; directories are excluded) **plus** the 11 hard-coded `GOVERNING`
     paths (`:14-36`) — 36 entries today. A path added to or removed from either set must be added to or
     removed from the manifest: the two name sets are compared for **equality** (`:124`), not containment.
  2. **One line per entry:** the SHA-256 of the file's raw bytes (`:129`), then **exactly two spaces**
     (`:111`, `partition("  ")`), then the path **relative to the project root, POSIX form** (`:39-40`).
     The digest must be **64 lowercase hex** characters (`:114-115`), names must be unique (`:117`), and
     the whole file must be **ASCII** (`:107-109`). `sha256sum` run from the project root already emits
     exactly this two-space form.
  3. ⛔ **SORT THE PATHNAME LIST *BEFORE* HASHING — byte-wise, `LC_ALL=C`** (`:61`,
     `key=lambda path: relative(path).encode()`), and end the file with a single newline.
     ⚠ **Sorting `sha256sum`'s output instead sorts on the leading DIGEST and is wrong**, and the mistake
     is silent: the verifier parses the lines into a dict and **never checks their order** (`:110-120`), so
     a digest-sorted manifest **passes verification** while differing bytewise from what a correct
     regeneration produces. From the project root, this reproduced the committed authority byte-identically
     (verified 2026-07-29):
     ```bash
     F=research/pde_ledger_v2/notes/ablation_driver/fixtures_v4
     { find $F \( -type f -o -type l \) -print
       printf '%s\n' <the 11 GOVERNING paths from verify_freeze.py:14-36>
     } | LC_ALL=C sort | xargs -d '\n' sha256sum > $F.accepted.sha256
     ```
     Unless that `GOVERNING` tuple itself changed, those same 11 names are already the non-`fixtures_v4/`
     lines of the current authority — read them off it rather than retyping them.
  4. **Commit the governed byte change and the new authority as ONE self-consistent committed state.**
     The verifier compares the authority's worktree bytes against `HEAD:` (`:86-91`) and requires a clean
     worktree — but ⛔ **that cleanliness check is scoped to the suite and the authority only** (`:92-103`
     passes exactly `relative(SUITE)` and `relative(AUTHORITY)` as pathspecs). The 11 governing paths are
     **digest-compared against the worktree and never checked against git at all**, so the freeze can
     verify green while the governed source bytes it attests exist in **no commit** — a green freeze over
     an unreproducible tree. ⇒ Land the governed change and its re-accepted authority together (or the
     governed change first, then the authority), never the authority alone. An uncommitted authority is
     caught for you: it reports `PENDING-COMMIT` or *"freeze authority differs from HEAD"* however correct
     its digests are. ⛔ Suite-local `__pycache__`/`*.pyc` fails before any of this (`:64-72`, `:81-83`),
     and a **symlink** in the inventory is rejected outright (`:127-128`).
  5. **Confirm, from the project root** — the path is the one in `fixtures_v4/README.md`; a suite-relative
     or `notes/`-relative spelling exits **2** with no output:
     ```bash
     PYTHONDONTWRITEBYTECODE=1 python3 research/pde_ledger_v2/notes/ablation_driver/fixtures_v4/verify_freeze.py
     ```
     → `freeze verified: 36 files`, exit 0. That run is the acceptance test for your hand edit.
- ⛔ **Re-acceptance is a deliberate REVIEW STEP, not a digest refresh.** The reviewer is expected to check
  the **semantic** change before the baseline moves, and the reason is recorded in the **commit message**,
  because the authority file carries no reason field. The two existing re-acceptances (`41b66dd5`,
  `96a1a61b`) are written that way and are the template. Then consider whether **A7's expected projection
  must be re-derived**: A7 joins the driver's output against those committed legacy tables, so a change
  reaching them can move the expected answer rather than only its digest.
- ⭐ This coupling is **accepted and deliberate** (user decision 2026-07-29), not a defect, and it is not
  a licence to weaken the suite: per §12b(b) the building session may neither author nor weaken it, and
  what the freeze does and does **not** protect is stated in `fixtures_v4/README.md` and
  `fixtures_v4/FREEZE_LIMITS.md` — read it there rather than re-deriving it here.

⛔ **(a) ENUMERATE EVERY DIMENSION-VALUED OBJECT IN THE `.wl` FIRST — before deciding what to emit.**
The comparator checks **artifact name-set symmetry between two files**. A **symmetric** omission — a
quantity absent from *both* artifacts — produces neither `py_only` nor `wl_only`, so the waiver
mechanism *structurally cannot see it*, and the run prints `py_only=0|wl_only=0` with no waivers,
which reads as complete. stage018 emits **6 of its 10** objects and the line is identical either way.
⇒ List **every** dimension-valued object with its definition locus, and for each one you do not emit,
the reason plus its read locus. Home = the stage note (`notes/stages/ledger_stage0NN_*.md`);
**stage018's table is the template.** Doing this first is what makes (b) and (c) decidable.
⛔ A written record a reviewer checks — **NOT** a corpus-wide inventory, oracle or completeness proof.
Three of those were specified and rejected (§10); they fail the same way, by letting the enumeration
become artifact-supplied.

**(b)** Add `DIM|` output to the stage's `.wl` — **prefer print-only**, but under **D2** new computation is allowed where a value is otherwise unreachable (state the independent route, D3).
⭐⭐ **Derive BOTH the `axes=` label AND the exponent vector from ONE label→slot structure.** Deriving
only the label is **insufficient and was measured wrong**: stage016's first emitter used the mapped
label but printed the raw storage order, so permuting the map would have relabelled **all 21 records
while every standing gate stayed green**. Routing both through the same structure fixed it
(`5b29f400`). ⛔ A second, parallel axis list used only by the records is a defect, not a fix.

⭐⭐ **(c) THE PHYSICS LEG RUNS HERE — BEFORE THE FREEZE — AND IT IS BLOCKING.** On a fresh agent,
derive every proposed quantity's dimension from the **model** (`docs/model_map.md` §2, the stage's own
physics, `notes/parameter_register.md`), **and adjudicate every proposed NAME against D4.** This leg
is a **gate on step (d)** — it authorizes nothing by itself; (d) alone owns the re-baseline and the
reference commit.

⭐⭐ **(c1) THE PHYSICS LEG MUST LEAVE A TRACKED VERDICT — it is the only blocking leg that used to
leave none.** §4-a requires a tracked enumeration; this leg produced its verdicts into commit
messages and nothing else, which is why *"how many quantities are actually physics-verified"* was
not a number anyone could look up. Record in the stage note, alongside the (a) enumeration:
1. **per quantity** — `CORRECT` / `WRONG` / `UNDETERMINED`, with the derivation route, not just the
   verdict. A verdict without its route is an assertion.
2. **the D4 name determination** — same-quantity or different-quantity, with the evidence, for every
   name that could group with an existing one.
3. ⭐ **what this leg does NOT cover** — count the records that are *declared literals in both
   engines*, because for those the comparator catches transcription divergence and **nothing**
   catches wrong upstream physics. stage037: 8 of 21. stage018: 6 of 10. Without this line a green
   physics leg reads as full coverage.
4. ⭐⭐ **anything STRUCTURALLY UNCHECKABLE, which is the most valuable thing this leg produces.**
   Examples already found: stage018's 65 checks admit a **six-integer-parameter family** of
   declarations, so they cannot pin `[a]` or `[c_s]` at all; stage011's `OmegaDim` is consumed by
   nothing in stage011. **No amount of conversion fixes these** — they are a map of where the
   ledger's dimensional checks are hollow, and Part VII's whole-system firewall needs exactly that
   map. Treat it as a deliverable, not a complaint.
5. ⭐ **carry the scope inside every sharp scope claim.** Within a §4-(c1) tracked verdict, every
   `only` / `no other` / `never` / `none` claim must say in the same sentence which corpus,
   artifact, stage set, or engine set it covers. Four such claims in the stage023 note were true in
   their authors' intended scopes and false as written, and the stage023-session checks found none
   of the four; readers found them by opening the cited sources
   (`_scratch/stage023_orch/ORCH_FINDINGS.md:48-69`, **F5 — unqualified scope claims**).
   ⛔ Do not weaken a sharp claim into vagueness: its precision is the value; state the boundary
   (`_scratch/stage023_orch/ORCH_FINDINGS.md:60-63`, **F5's mitigation and wrong fix**).
6. ⭐ **make a same-note line reference recoverable by naming its target.** Within the growing stage
   note written by §4-(c1), a `:NNN` reference to a stable named section, numbered item, or named row
   must carry that name beside the number. Three same-note references decayed during the stage023
   session as insertions moved their targets, and validators in this dimension-rewrite workstream
   do not test whether a line number reaches the claimed sentence
   (`_scratch/stage023_orch/ORCH_FINDINGS.md:25-43`, **F6 — intra-file line references decay**).
   ⭐ **The same rule holds ACROSS files.** A `path:NNN` citation — into another note, into an engine,
   into the canonical table, into this manifest, or into the register — must name what it points at
   beside the number: a heading, a row key, a defined term, or a quoted fragment. **The measurement
   that motivates this:** four consecutive prose passes on stage023 were spent recomputing cross-file
   references that the next edit then re-invalidated — the note's references into its two engines,
   into the canonical table and into this manifest, and this manifest's and the register's references
   back into the note. A named target is what makes the reference survive that: a reader who finds the
   number stale can still reach the target, and the next editor can re-verify by content instead of by
   offset.
   ⛔ **This does not license dropping the line numbers.** They are still worth carrying — the number
   is what makes a fresh reference cheap to check, and the name is the recoverable anchor added beside
   it, never a replacement for it.
   ⛔⛔ **THE LIMIT, MEASURED — an anchor makes a stale reference RECOVERABLE; it does not make it
   CORRECT.** This clause can be satisfied while the pointer is false, so satisfying it is not evidence
   the pointer resolves. **This changeset is the counterexample.** The stage023 `.py` rewrite moved
   `SOURCED_DIMS`, and four §12 WORK citations that *already carried named targets* went on pointing at
   the wrong lines afterwards: `WORK-023-D0-SEAM`'s window for the `D0` declaration had come to land on
   `M0` instead, its two windows for the port formula and the dimensionless `P0_physical` target both
   stopped short of the lines holding them, the `[M0]`/`[D1]` window shared by
   `WORK-023-MOMENT-CONVENTION` and `WORK-023-STAGE009-MOMENT0` no longer reached `D1`, and
   `WORK-023-SOURCED-PROVENANCE`'s "in `SOURCED_DIMS`" range ended before the block's last
   declarations. ⇒ The number must be **re-verified by content whenever either file changes** — an edit
   to the *cited* file obligates the *citing* file — and **a named anchor beside a wrong number is not
   compliance**. The anchor buys the reader recovery; it buys the writer nothing.
   ⚠ A line range into a **regenerated** artifact cannot be made to resolve at all; cite the revision
   it was read at, plus what it showed, instead of live line numbers.

⛔ **Why this matters more than the conversion.** Cross-engine agreement has now been shown twice
to be necessary and **not sufficient** — it is blind to a same-dimension different-quantity merge
(`c_s0` vs `c_S`, stage018) and to a shared wrong declaration. **"30 of 30 converted" and "the
ledger's dimensions are demonstrably right" are different finish lines.** The conversion serves the
checker's recovery path — an engineering need. This leg serves the correctness claim, and it does
**not** depend on a stage being converted: it can be run against any stage's existing declarations.
⛔ **A naming decision is a PHYSICS decision. Never freeze one on dimensional evidence alone.**
**The measured reason (stage018, 2026-07-27):** the `.wl` proposed emitting the sound speed under
stage012's name, merging `c_s0` (bulk asymptotic `ρ0`) with `c_S` (wall `ρ*`). Both are `L T⁻¹`, so
**the comparator is blind to it by construction**, and the generator groups on the *scope-stripped*
name — it would have published a false `AGREE` into the table Part VII's dimensional firewall
consumes. **Two independent parties made this identical error** (the build, from the dimension table;
the orchestrator, in its pre-registration), neither having read the other. *"Same dimension + same
words sound speed"* is an **attractive** wrong merge. Running this leg at (f) instead cost a full
fix round to unfreeze the name; the leg is the same, only earlier.

**(d)** Re-run; confirm exit 0, PASS
tally unchanged, `.out` reproduces byte-identically after `sed -E 's/\$[0-9]+/$N/g'`; re-baseline.
**Commit this — plus the (a) enumeration — before touching the `.py`.** Freezing the reference first prevents both conversion sides from
changing together — it strengthens reference custody; it does **not** establish that the original
engine pair was independently authored (§1).
**(e)** Seal the prediction — ⛔ **OUTSIDE the repo** (scratchpad) until the stage's build *and*
its review legs have landed; copy it in and commit it at (i) as the record.
⚠ **Custody caveat, stated honestly:** committing it afterwards means **git cannot prove it predated
the build**. Its value is as a *working* pre-registration — it disciplines your own thinking and
records falsified predictions — not as cryptographic evidence. If that ever needs to be stronger,
record a hash of it outside the workspace before the build starts.
⛔ **A pre-registration left in the working tree is a SUPPLIED ANSWER.** At 018 a review leg read the
untracked note and cited it as authority for a claim that was wrong. Withholding it from the
*directive* achieves nothing — agents `ls`, `grep` and `git status` their way in. ⭐ What saved it: the
physics leg was told to **derive from the model**, not to *check a claim*. A leg pointed at first
principles survives contamination; one pointed at "verify this" does not.
**(f)** Rewrite the `.py` onto `scripts/ledger_dimensions.py`.
**(g)** **Re-run the `.py`, then** compare axis-labelled:
`python3 scripts/compare_dimension_artifacts.py <NNN>`. The sidecar is now **source-hash bound** — the
comparator recomputes sha256 for both the stage `.py` and `scripts/ledger_dimensions.py`, rejecting
either missing or mismatched assertion — but run the stage first anyway and say that you did.

⛔ **A bare stage run is a PRODUCER, not a validator.** It writes the dimension sidecar; its exit code
and PASS tally do **not** certify that artifact. Validation evidence comes from the module-pin control,
the axis-labelled comparator, and the canonical-table generator. After a legitimate, independently
reviewed edit to `scripts/ledger_dimensions.py`, an expected red pin is re-baselined explicitly with
the full, unabbreviated command:
```bash
python3 scripts/check_ledger_dimensions_pin.py --accept
```
Explicit acceptance intentionally hashes the current module bytes and replaces
`scripts/ledger_dimensions.accepted.sha256` (`scripts/check_ledger_dimensions_pin.py:101-115`); it
does not validate the change's correctness. By contrast, validation/check mode reads the expected
digest from the authority and hashes the module only for the actual digest
(`scripts/check_ledger_dimensions_pin.py:58-72`, `:75-98`): it never derives the expected digest
from the module. Inspect that one-line authority change, then re-run the standalone control, the stage
producer, the comparator, and `python3 scripts/generate_canonical_dimension_table.py`. The parser sets
`allow_abbrev=False`: spell `--accept` in full; `--acce` is rejected.
**(h)** Review: transliteration-fidelity fresh agent + adversarial-with-ablation fresh agent.
**(h2)** ⭐ **ONLY NOW, open the sealed prediction from (e) and adjudicate it** — record every
prediction confirmed and every one **falsified**, with evidence. Doing this before (h) would leak
the expected answers into the reviews; doing it never is how a pre-registration becomes decoration.
At 018 this is what caught P9 (the false `c_s0`/`c_S` merge) as a *wrong prediction of mine*, not
just a build defect.
**(i)** Commit **only after (g), (g2), (h) and (h2) are clean**, including the adjudicated
prediction note. ⚠ Any blocking finding at (h) stops the loop: remediate, then repeat every affected
execution and review step until both legs are clean.

*The two blocks below expand steps (g) and (c). **(g2) runs inside step (g)** — after the comparator,
before review; **(c2) runs inside step (c)**. Neither is an extra step after (i).*

⭐ **(g2) THE ORCHESTRATOR MUST REGENERATE THE `.out` ITSELF, once per stage.** The `.out` is the
reference half of the only universal gate; if it could be hand-written, the cross-check proves nothing —
as one agent put it, *"a hand-edit and a real re-run are byte-identical, so its provenance rests on trust
alone."* Run `math -script <the .wl>`, normalise with `sed -E 's/\$[0-9]+/$N/g'`, and confirm it
reproduces the committed `.out` byte-for-byte. Done for stage013 (sha `42ee1ad7fbf8283a`, exit 0).
⚠ **NEEDS RE-ADJUDICATION (2026-07-29) — the step's stated RATIONALE is refuted, not the step.** This
step used to read *"verification agents are barred from Mathematica, so they **cannot** confirm the
reference side is genuine."* That premise was measured false: a fresh review agent ran
`timeout 600 math -script …` to exit 0 in this environment and byte-diffed the transcript itself
(`research/pde_ledger_v2/_scratch/stage023_h/REPORT_ADVERSARIAL.md`, *"Positive control on the
Mathematica side"*). ⚠ **Do not overstate the measurement:** one agent, one occasion. It refutes
"cannot"; it establishes nothing about reliability or about licence contention — both remain
undetermined from a single observed run. ⛔ **The control stands meanwhile**, but that is the *practice*
continuing on a reason not yet re-established. ⚠ **A replacement rationale is PROPOSED, not
established:** that the reference half of the only universal gate must be regenerated by a party
independent of the build, since a build that produces both halves proves nothing by comparing them. It
is plausible; it is **not** what the step was originally justified by; the user must **adjudicate it
before it is cited as the reason**. ⚠ **Seat cap, now
relevant:** if agents run Mathematica, they consume licence seats — **never more than 2 concurrent
`math -script`**, counting every session running `math`, orchestrator included.

⭐ **(c2) MORE EVIDENCE THAT STEP (c) IS NOT OPTIONAL.** On stage013, **5 of the 15
compared records** (`symbol_dims.*`) are declared as literals in *both* engines. Comparing them catches
a transcription divergence between the two — but **a SHARED wrong declaration passes the comparator**.
Only an independent read against the model catches that. stage013 scored 9 CORRECT / 0 WRONG on that
leg; do not skip it because the cross-engine numbers look clean.

**Artifacts.** `.py` → `scripts/<basename>.dimensions.txt`, same flat format the `.wl` prints:
`DIM|axes=L,M,T|name=EnergyDim|exponents={2, 1, -2}`. One parser, both engines.
⛔ Not stdout · ⛔ not JSON (these are read by humans and review agents).

**Three gates, and each proves something different:**
1. **PASS multiset identical** — the audit still reports what it reported.
2. **stdout byte-identity** — behaviour preservation. ⚠ **A DIAGNOSTIC, NOT A BLOCKING GATE** (D1,
   §1b): report it, never let it prevent a correct change. Also **NOT a transposition detector** (§6).
   ⛔ **The recipe here was WRONG and is corrected — MEASURED 2026-07-27.** It read
   ~~`tail -n +7 scripts/output/<basename>.txt` (diff exactly 6 wrapper lines)~~. The committed
   transcripts actually carry an **8-line leading wrapper** (7 `#` lines **plus a blank**) and a
   **2-line trailing wrapper** (blank + `EXIT_CODE:`). The correct form is:
   ```bash
   tail -n +9 scripts/output/<basename>.txt | head -n -2
   ```
   ⭐ Verified by execution against **all five** converted stages — 004, 011, 012, 013, 018 — where it
   reproduces each live run **exactly**; the old form leaves a stray `#` and a blank line and so always
   reports a spurious difference. 26 of the 28 committed transcripts share this wrapper shape.
   ⛔ **This gate only exists for stages 001–028.** `scripts/output/*.txt` covers 001–028 ONLY —
   **stages 030–044 have no committed Python transcript**, so gate 2 is unavailable there and gates 1
   and 3 carry the whole load. Plan for that before starting group B/C/D; do not discover it mid-stage.
3. **axis-labelled cross-engine comparison** — the only universal transposition detector.

**⛔ The comparator must have a non-empty floor, and does.** The first version reported `PASS` with
`compared=0` (header-only sidecar / no name matches / `.out` stripped of `DIM` lines). It was live:
stage011 passed at `py=12 wl=2` with both `Corrupt*` teeth uncompared. Now a zero shared-name count fails, and
unwaived `py_only`/`wl_only` fails. Waivers are per-stage, name every quantity, and are echoed every
run. **Current output lines (renamed 2026-07-27 so they cannot be misread as source coverage):**
```
ARTIFACT_NAME_SET|stage=…|py=N|wl=N|shared=N|py_only=N|wl_only=N|source_coverage=not_checked
ARTIFACT_NAME_WAIVERS|stage=…|py_only=…|wl_only=…
```
⭐ **Sidecars are also source-hash bound**: `emit_dimension_sidecar` stamps separate `sha256`
assertions for the stage and `scripts/ledger_dimensions.py` into the header, and the comparator
**recomputes both**, rejecting either missing or mismatched digest. This closes both measured holes:
a transposed-but-not-re-run stage and an edited-shared-module-without-re-runs scored `PASS`.

## 5. ⛔ THE FABRICATION GUARD (step b) — where this can silently destroy itself

⚠ **Amended by D2/D3 (§1b).** New computation in a `.wl` is now ALLOWED where it is needed to expose a
value — but it must be an **independent route**, and the directive must require Codex to state *how
that route differs from the `.py`'s*. Everything below still stands; only the blanket "no new
computation" clause is lifted.

- **Prefer print-only** where a live binding already exists — it is cheaper and carries no
  independence risk. Add computation only when the value is genuinely unreachable otherwise.
- ⛔ **Never hardcode an exponent literal in a `Print`.** A constant compared to the `.py`'s constant is
  vacuous. This defect already exists at **stage013 `.wl:448-449`** and **stage018 `.wl:386`** — do not
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
verify. Sources: `notes/stages/`, `*_source_map.md`, `software/*/reports/`, `decisions/`, the commit
messages, the `_scratch/` dirs, `research/` prior papers.

## 8. COVERAGE CENSUS — the "group A is free" premise was FALSE

| stage | py dims | **real** compared | labelled? | blind |
|---|---|---|---|---|
| 012 | 19 | **18** ✅ done | DIM records | L↔M 16, M↔T 8, L↔T 16 detect |
| ~~013~~ | — | ✅ **DONE — 15/15, zero waivers** (`4391c69c`) | — | — |
| ~~018~~ | — | ✅ **DONE — 6/6, zero waivers** (`63dee5e4` + `1b645ed9`) | — | — |
| 016 | 21 | 9/21 | labelled | 4 |
| 023 | 29 | 7/29 | BARE | 6 of 7, **5 dimensionless** |
| 027 | 17 | **1**/17 | bare + hardcoded gloss | 0 |
| 021 | 18 | 7/18 | labelled | **7 of 7 — verdict quantity dimensionless** |

⛔ **021's gate cannot detect a transposition at all.** 027 is 1-of-17. 013/018 print hardcoded
**exponent-value** literals in both engines — dimension *values* typed into prose prints
(`013 .wl:448-449` / `.py:588-589`, `018 .wl:386` / `.py:543`), the §5 defect. ⚠ **This is NOT the
`axes=`-label claim and must not be merged with it:** on the label the measured set is **004/013/018**
(§3b), and 004 is absent here because its prose prints call `dimString[d[…]]` rather than typing values.
⇒ every group-A stage needs `.wl` emission; with the floor in place they now
**fail** rather than pass quietly. The old ordering was a cost heuristic, never a correctness one.
⚠ 021 also has cross-engine name collisions (`[P₀_raw]`/`[P0_raw]`, unicode `N₀` vs ASCII `N0`).

✅ **GAP CLOSED 2026-07-27 — the four remaining group-A stages are now surveyed.** They were the
fossil of the falsified *"group A is free"* premise: `notes/rewrite_reference_table.md` §5.3 surveyed
only the 21 stages then believed to need `.wl` work, so **016, 023, 027 and 021 had no reachability
verdict at all**. Three independent read-only agents filled it; loci spot-checked by the orchestrator.

| stage | verdict | values reachable | axis order — and HOW it is established | re-invocation | the hazard that decides the build |
|---|---|---|---|---|---|
| **016** | `REACHABLE` | **21 of 21** — 12 declared rule-table entries (`wl:305`) + 9 computed (`wl:306`) | `(L,M,T)`, ⭐ **code-evidenced**: the slot→label binding `{{"L",d[[1]]},{"M",d[[2]]},{"T",d[[3]]}}` at `wl:123`. Every slot carries a non-zero value somewhere | `evalDimensional` **6×** — 1 clean, 4 corrupted, 1 arity-self-check re-run | the `.out` today renders **zero** exponent vectors — only `dimText` monomial strings and prose |
| **023** | `REACHABLE` | **7 of 7** computed, via the top-level global `baselineDimAudit` (`wl:286-288`), + 22 declared `baseDims` + 7 declared `expectedDims` | `(L,M,T)`, ⭐ now **code-bound** in the committed `.wl`: the step-(d) emitter's `dimensionAxisSlots = {{"L",1},{"M",2},{"T",3}}` (`wl:212`) feeds both the `axes=` label and every exponent vector. ⚠ **The survey's observation was accurate for its referent** — the 769-line **pre-emission** `.wl`, where the order was *prose only* (then `wl:608`, now `wl:654`) and evidenced by non-zero slots in `baseDims`. ⚠ of the **7 computed outputs** only `A1` has a non-zero L — the order is evidenced by the inputs, not the results | `dimensionAudit` **17×** — 3 top-level + 14 memoized `caseFor` modes | ⛔ a **literal substitutes for a computed value** when the expression is 0 (`wl:272`): inert at baseline, **load-bearing in `"perfect"` mode**, where 4 of 7 rows then self-match against `expectedDims` |
| **027** | ⚠ **MIXED** — declared `REACHABLE`, computed **`LOCAL_ONLY`** | **0 of 1** computed vectors reach top level (it dies in `runAll`'s `Module`, `wl:742-751`); the 16 declared `baseDims` (`wl:183`) do | `(L,M,T)`, ⭐ **mechanically bound** through `uL/uM/uT` (`wl:200`, `:205`) — the **strongest axis evidence in the corpus**, code not prose | `evaluatePort` **19×**, 2 of them under a corrupted basis | ⛔ the `.wl`'s rescaling-ratio route **cannot produce per-symbol vectors at all** — 027 stays a **1-row** `DIM\|` stage unless new call sites are added (which D2 permits) |
| **021** | `REACHABLE` | **27 of 27** top-level bindings — ⭐ **no top-level `Module` anywhere**; ~21 clean named + ~64 mutation-scoped ≈ 85 | storage `(L,M,T)`, **prose only** (`wl:342`, `:384`, `:528`) — there is **no machine-readable axis binding**, so a D2 emitter must hardcode the header | `gateData` **13×** (8 under a corrupted map), `backSolveMutant` **5×** (all corrupted) | ⛔⛔ an **undocumented index-permuting renderer** — `{{"L",d[[1]]},{"T",d[[3]]},{"M",d[[2]]}}` at `wl:125`/`:139` — so scraping its printed output and labelling it `axes=L,M,T` **silently swaps M and T**. The `.py` documents its permutation (`py:224`); the `.wl` does not |

⛔⛔ **SURVEY COUNTS ARE PROVISIONAL, NOT COMPLETENESS PROOFS — measured on stage016, 2026-07-27.**
The stage016 survey correctly identified the **21 clean emitted quantities**; the per-stage §4-a work
then found two additional **non-emitted** source/control-flow cases — `lambdaRef`, a live numeric
factor whose dimension comes from the `NumericQ` fall-through, and a duplicate clean `evalDimensional`
re-run whose `Dims` are discarded after its `Ok` is read. ⚠ **The emitted count did NOT go 21 → 23**,
and ⚠ **neither case is unconsumed** (`lambdaRef` feeds `k2Ref = buildK2[lambdaRef]` at `.wl:229`, and
the re-run's `Ok` is checked). What the miss shows is that the **survey method did not close the
broader audit inventory**. ⇒ Treat the 023/027/021 counts as **provisional** pending each stage's
tracked §4-a enumeration plus adversarial review — which is reviewed evidence, not a completeness
proof (§4-a says so itself).

⛔⛔ **CORPUS-WIDE CORRECTION — "append the print at end-of-file" is DEAD CODE in all 43 `.wl`
files.** §5.3 defines `REACHABLE` as *"a print appended at end-of-file **or at the call site**
works"*. **Measured 2026-07-27: 43 of 43 `.wl` files terminate in `Exit[0]`/`Exit[1]`**, so the
end-of-file half of that definition is false everywhere — a print after the terminal `Exit[]` never
executes and emits nothing.
⇒ **The verdicts are unaffected** (the *values* really are reachable, and no restructuring is
needed); only the **insertion locus** changes: emit *before* the terminal `Exit[]` block, at or after
the top-level binding. The working precedent is stage018, whose `DIM|` prints sit at `.wl:387-393`
against its `Exit[]` at `.wl:499`/`:502`.
⚠ Three of the four surveys flagged this independently, which is why it is recorded here rather than
in a stage note: it applies to every stage still to be converted, including the 21 already surveyed.
⭐ The acceptance criterion that catches it without anyone having to remember it: **the count of
emitted `DIM|` records must equal the count of objects marked *emitted* in the §4-a enumeration** — a
dead print yields 0 and fails.

**Follow-on findings from the same surveys** (recorded, not acted on):
- ⛔ **task #22 material, 021** — `py:557` compares `data["raw_symbol_dims"][N0]` against the constant
  `SOURCED_N0_DIM`, which is **the same constant** (`py:145`, `py:381`): it cannot fail for any value.
  The Mathematica twin `wl:342` retypes the literal and **is** a real transcription guard. The two
  engines are not equally strong here, and the weaker side is the `.py`.
- **Fractional stored values confirmed in all three of 023/027/021**, so the emitter must serialise
  exact rationals, never floats. ⚠ **Correction — stage023's `gU`/`gW` content does not cancel against
  itself.** In `(L,M,T)`, the declarations are `[Ω_U]=[Ω_W]={0,0,-1}`,
  `[R_mix]={0,0,-2}`, `[g_U]=[g_W]={-1/2,1/2,-2}`, `[D0]={-1,1,-2}`,
  `[c_s]={1,0,-1}`, `[a]={1,0,0}` (`023 wl:238-244`), and the port is built at
  `023 wl:105-117`. The full arithmetic is:
  `[Ω_U²g_W]=2{0,0,-1}+{-1/2,1/2,-2}={-1/2,1/2,-4}` and
  `[R_mix g_U]={0,0,-2}+{-1/2,1/2,-2}={-1/2,1/2,-4}`, so
  `[pport]={-1/2,1/2,-4}`; `[dport]=[Ω_U²Ω_W²]=[R_mix²]={0,0,-4}`;
  `[nport]=2[pport]-2[dport]={-1,1,-8}-{0,0,-8}={-1,1,0}`;
  `[p0]=[nport]-[D0]={-1,1,0}-{-1,1,-2}={0,0,2}`; and
  `[P0Physical]=2([c_s]-[a])+[p0]={0,0,-2}+{0,0,2}={0,0,0}`.
  What actually cancels is the `{-1,1}` contribution from the squared `g` content in `nport`
  against `D0`. If `gU` and `gW` are wrongly made dimensionless, then
  `[pport]={0,0,-2}`, `[nport]={0,0,-4}-{0,0,-8}={0,0,4}`,
  `[p0]={0,0,4}-{-1,1,-2}={1,-1,6}`, and
  `[P0Physical]={1,-1,4}=L M⁻¹ T⁴`, not dimensionless. The other fractional cases remain
  `I25` `{5/2,0,0}` and `CouplingAPower` `-7/2` (`027 wl:186`, `:154`), and `muDim`
  `-1/2` (`021 wl:303`, rendered `L⁻¹ T⁻¹ M⁻¹ᐟ²`).
- **Name-pairing debt, measured.** 023: `.wl` `"P0Physical"` vs `.py` `"P0_physical"` — the other six
  join. 027: the `.wl` names objects by variable only and has **no counterpart** to the `.py`'s nine
  expression-embedded record names (`py:526-536`). 021: 3 of 10 symbol keys differ (`cs`/`c_s`,
  `chiQsym`/`chi_Q`, `muHat0`/`mu_hat0`), every derived-object name differs, and the `.wl` is split
  **against itself** — `[P₀_raw]` in a `Print` (`wl:335`) vs `[P0_raw]` in the assertion name
  (`wl:344`). These are D4 join-key work, and ⛔ aligning a *name* never licenses aligning a *value*.
- **021 dead weight:** `lhsRawDim` (`wl:308`) and `gamma5Dim` (`wl:310`) are computed and consumed by
  **nothing** — yet the `.py` *does* export `Gamma5` (`py:534`). A real cross-engine asymmetry.
  `dimText` (`wl:124`) is defined and never called.

**Conversion order — which stage converts next:** ~~013, 018,~~ ✅ done → ~~the three reopened waivers as one batch~~ ✅ **CLOSED `8b006055`**
→ ~~016~~ ✅ done → ~~023~~ ✅ done (both engines, comparator green) → **027 next**, then 021 (heaviest). Then `(L,T,M)`, `(M,L,T)`, then 008 (2-axis), 038 (4-axis), 042 (stiffness).
⚠ **This is the conversion order only, and it is not the build queue.** Which *tooling* gets built next
is a separate sequence, set by §12b: **the ablation driver first, then the shared Mathematica `DIM|`
emitter** — both built **before** the 027 conversion begins. Neither line overrides the other — they
order different things.
⚠ **The 037 spike ran out of order on purpose** (§3b) — that was a feasibility measurement, not a
conversion, and 037 stays in its group-B slot for the actual work.

✅ **038 and 042 are already CLEARED as non-blocking — do not re-litigate this.** Proven by execution
(`scripts/probe_ledger_dimensions_extremes.py`, committed): the module handles 038's 4 axes
`(M,L,T,E)` with a heterogeneous `None`-sentinel slot, and 042's `(stiffness,L,T)` with
`fractions.Fraction` exponents plus the `len(set(...))` homogeneity test that requires hashability.
Each capability is ablation-verified (`_exact`→`sp.Float` trips the exactness check; `__hash__`→`id`
trips the set check). They are ordered LAST for effort, not for risk.

## 9. ⛔ LANDMINES
- ⛔ **MEASURED — exact-rational SYNTAX is not validated.** The two-line input
  `DIMENSIONS|axes=a,b,c,d,e` /
  `DIM|axes=a,b,c,d,e|name=probe|exponents={0.5,1/2,-2,1.0e0,2.5}` was passed through
  `compare_dimension_artifacts.load_dimensions`. Literal output:
  `'0.5' -> Fraction(1, 2)`, `'1/2' -> Fraction(1, 2)`,
  `'-2' -> Fraction(-2, 1)`, `'1.0e0' -> Fraction(1, 1)`,
  `'2.5' -> Fraction(5, 2)`; exit **0**. The loader sends every stripped token straight to
  `fractions.Fraction` (`scripts/compare_dimension_artifacts.py:87-95`, called from
  `load_dimensions` at `:99-145`), and the canonical-table generator delegates to that same loader
  (`scripts/generate_canonical_dimension_table.py:233-234`). **Status today:** *serialise exact
  rationals, never floats* is an authoring rule only; no validator rejects decimal or scientific
  syntax, and normalization can make a non-conforming artifact look exact downstream. §12 tracks
  the validator fix; this is not checked today.
  ⛔ **MEASURED, end-to-end — float exponents now survive the whole comparator, not just the loader.**
  The earlier observation stopped at `load_dimensions`; a stage023 run with float exponents emitted into
  **both** artifacts, and the ledger-dimensions module pin re-accepted afterwards, reaches
  `compare_dimension_artifacts.py 023` reporting `status=PASS`, exit **0**. Same hole, one layer further
  out: `Fraction` normalization equates the float and rational spellings before anything compares them, so
  a green comparator is not evidence that either artifact serialised exact rationals. It does **not** imply
  the values are wrong — normalization preserved them here; what it implies is that the authoring rule is
  the only thing enforcing the syntax, end to end.
- ⛔ **MEASURED — the comparator has no record-coverage floor above 1.** With **28 of the 29** stage023
  records dropped from **both** artifacts, `python3 scripts/compare_dimension_artifacts.py 023` printed
  `shared=1` and `RESULT|stage=023|status=PASS|mismatches=0`, exit **0**. The earlier zero floor is in
  place and has fired — `compare_dimension_artifacts.py:313-314` appends `compared=0; no shared quantities
  were compared` — but it is a floor at **zero**, so a `compared=1 of 29` state passes unremarked. What
  this implies: a green comparator bounds only the records that reached it, never how many should have.
  What it does **not** imply: that the record still compared was wrong — it matched. The drop was applied
  to **both** artifacts at once, which is exactly the state no cross-engine comparison can see; the
  observation is about the missing floor, not about either emitter.
  ⚠ **A third validator-layer item, orthogonal to both §12 tracked
  specs** — the exact-rational spec governs exponent-token lexical exactness, the axis-order spec governs
  cross-engine axis-order metadata, and neither counts records. ⛔ Do not fold it into either acceptance
  list; a per-stage expected-record-count floor is a separate change.
- ⛔⛔ **MEASURED — an uncaught `Throw` can make `math -script` green with no trustworthy
  transcript.** The direct probe run was `math -script throwtest.wl`, where the complete program was
  `Print["before"]; Throw["boom", "sometag"]; Print["after"]; Exit[0];`. Stdout was literally
  `before` plus a newline — no diagnostic or `after` — and the exit status was **0**. In the
  stage023 perturbation, changing only `baseDims[K0c]` from `{0,1,-2}` to `{1,1,-2}` made the
  `K0c + Z0ret` sum heterogeneous: `math -script probe_K0c.wl` produced **0 stdout lines**, **0
  stderr lines**, exit **0**; the unmodified control, run the same way from the same directory,
  produced **214 stdout lines**, exit **0**, ending
  `TALLY mathematica: 117 pass + 0 fail = 117 checks` / `OVERALL PASS`.

  **Two routes, two observables — do not collapse them.** **Route A (load time):** an immediate
  top-level `Set` or bare expression reaches a throwing primitive before the top-level `Catch`
  exists; in these files nothing has printed yet, so the transcript is **empty**, exit 0.
  **Route B (tag mismatch):** a `Throw` tag is not named by the file's `Catch`; when it fires from
  inside the guarded run, the already-started transcript is **truncated**, exit 0. The static survey
  found **6 EXPOSED / 38 NOT EXPOSED / 0 UNDETERMINED = 44 `.wl` files**. Route A reaches 5 files
  (018, 021, 023, 038, 039); Route B reaches 3 (018, 021, 033), with 018/021 in both:

  | file | route and verified source loci | conversion position |
  |---|---|---|
  | `ledger_stage018_dtn_hankel_fingerprint_mathematica_audit.wl` | **A+B** — uncaught `"stage018DimError"` at `:100-120`; load-time `dimOf` calls at `:229-233`; top-level `Catch` names only `"ledgerStage018Failure"` at `:489-492` | **already converted** (§3) |
  | `ledger_stage021_dimensional_closure_mathematica_audit.wl` | **A+B** — `dimOf` throws `"dimOfFailure"` at `:87-105`; immediate calls span `:265-326`; the unshielded in-run call is `:483`; top-level `Catch` names only `"ledgerStage021Failure"` at `:539-542` | **not yet** — after 023/027 (§8) |
  | `ledger_stage023_nullspace_underdetermination_mathematica_audit.wl` | **A** — `fail[msg_] := Throw[msg, "ledgerStage023Failure"]` is `:38`; `dimOf` reaches it at `:222-245` (its four `fail` calls: `:226`, `:231`, `:239`, `:243`); the load-time `dimensionAudit` assignments `baselineDimAudit`/`corruptSourcedDimAudit`/`corruptFreeDimAudit` are `:286-296`; top-level `Catch[runAll[], "ledgerStage023Failure", …]` is `:802-806` | **already converted** (§3) — the `.wl` gained `emitDimensionRecords[]` (defined `:298-332`, called from `runAll` at `:789`), which does **not** change the route: the emitter is inside `runAll`, while the exposure is the pre-`Catch` load-time block above it |
  | `ledger_stage033_native_p_no_emergent_gauss_mathematica_audit.wl` | **B** — four `"pipelineFailure"` throws at `:180`, `:192`, `:440`, `:450`; the sole top-level `Catch` at `:727-1095` names `"ledgerStage033Failure"` (`:1093`) | **not a conversion target** — §3 records no dimension machinery |
  | `ledger_stage038_sealed_landing_electric_bc_r1_mathematica_audit.wl` | **A** — `manifestDisposition` fallthrough raises at `:819-833`; load-time `sourceManifest` is `:836-838`; top-level `Catch` is `:927-1185` | **not yet** — tail position (§8) |
  | `ledger_stage039_b_t_time_reversal_even_departure_mathematica_audit.wl` | **A** — `sourceDisposition` fallthrough raises at `:515-531`; load-time `sourceManifest` is `:534-540`; top-level `Catch` is `:578-871` | **not yet** — later axis-group work (§8) |

  ⭐⭐ **BOUND: this is a hole in exit-code gating, not in the calibrated conversion process.**
  §4-(d)/(g2) regenerates the Mathematica transcript and byte-compares it with the committed
  non-empty `.out`; an empty or truncated result differs at the first missing byte. Exit 0 alone is
  not evidence, but the committed evidence is not thereby untrustworthy. The six source hazards
  remain catalogued rather than corrected, and §12 tracks an executable transcript-liveness gate.

  **Survey method and limits.** This was static source analysis: position-preserving comment/string
  blanking and top-level statement splitting; under- and over-approximating name-based call graphs;
  pre-`Catch` prefix inspection; and literal tag matching. Both call-graph passes returned the same
  six files, bracket depth returned to zero in all 44, and none was left undetermined. ⚠
  **Reachable does not mean the guard fires on committed inputs** — firing is runtime-contingent
  (stage023 required a perturbation). Dynamic dispatch through constructed names could evade the
  name-based graph (`Symbol`/`ToExpression` were not found here); computed `Catch` tags would evade
  literal matching (none occurs here); and the survey does not analyze other termination mechanisms
  (`CheckAbort`/`Abort`/`Quit` do not occur in these 44). `NOT EXPOSED` is only the stated structural
  result: no immediate top-level path reached a throwing primitive under either pass.
- ⛔ **A print appended at end-of-file NEVER RUNS — measured, 43 of 43 `.wl` files end in
  `Exit[0]`/`Exit[1]`.** Emit *before* the terminal `Exit[]` block; stage018 (`.wl:387-393` vs its
  `Exit[]` at `:499`) is the working precedent. This kills the "or at end-of-file" half of §5.3's
  `REACHABLE` definition without changing any verdict — the values are still reachable. §8.
- ⛔ **`REACHABLE` describes the VALUE, not the print site.** Three separate things must each be
  checked before choosing an emission locus: (i) does the value survive to top level, (ii) does the
  chosen line actually execute, (iii) how many times does the holder run, and under what corruption
  flags. Measured re-invocation counts on the group-A tail: 016 **6×**, 021 **13×**, 023 **17×**,
  027 **19×** — an in-body print emits duplicate *and* deliberately-corrupted records in every one.
- **stage003's basis is `(M,L,T)`, NOT the more common `(L,T,M)`.** ⚠ The earlier claim that "neither
  file says so" was FALSE: `scripts/ledger_stage003_*.py:87` declares
  `def dim(m_power, l_power, t_power)`, stating the order explicitly. Read it; never assume.
  Emitting `axes=L,T,M` there would corrupt every triple.
- **stage042's `.wl:816` comment says "MLT" — a mislabel;** its axes are `(stiffness, L, T)`.
  ⚠ Its guard runs **once**, not twice (an earlier claim was wrong) — 042 is likely recoverable.
- ✅ **037, 036, 035 are NOT `.wl`-emission-impossible — that verdict is MEASURED FALSE** (spike,
  §3b). `ROUTE_EXISTS`, prototyped, 21/21 on 037, ~0.5–1 engineer-day per stage. 044 untested and
  frozen pending 044-v2. ⛔ `notes/rewrite_reference_table.md` §5.6 still carries the old "genuinely
  impossible" framing for these three — **it is stale; this section wins.**
- ⛔⛔ **THE SHARED MODULE WAS SELF-ATTESTING — THE ORIGINAL SINGLE-POINT-OF-FAILURE RESULT REMAINS
  LOAD-BEARING.** Measured 2026-07-27 before the external control existed: stubbing
  `ledger_dimensions.dim_residual` to `return sp.Integer(0)` made **nine** of stage016's dimensional
  gates vacuous while the stage still reported 82 PASS / exit 0 and the comparator exited 0. The
  sidecar's `ledger_dimensions_sha256` binding could not detect it because both asserted and current
  digests came from the same mutated file. ⇒ **That binding closes STALENESS, not TAMPERING.**

  **The missing validation/check control now has an expected side outside the module it polices.**
  In validation/check mode, `scripts/ledger_dimensions.accepted.sha256:1` supplies the expected
  digest and the checker hashes the current module source bytes only for the actual digest
  (`scripts/check_ledger_dimensions_pin.py:58-72`, `:75-98`); explicit acceptance intentionally
  derives the digest from the current module and writes that authority
  (`scripts/check_ledger_dimensions_pin.py:101-115`). The comparator invokes the check before
  comparison (`scripts/compare_dimension_artifacts.py:249-258`), the generator invokes it before
  loading stages and before writing the table (`scripts/generate_canonical_dimension_table.py:229-233`,
  `:459-468`), and `scripts/run_all_audits.sh:20` invokes it under `set -euo pipefail`
  (`scripts/run_all_audits.sh:12`). With the same `dim_residual` stub, the stage still produces 82 PASS
  / exit 0, but an un-rebaselined control and comparator now fail `MODULE_PIN_MISMATCH`. This
  establishes one narrow fact: the SHA-256 of the current `scripts/ledger_dimensions.py` source bytes
  equals the deliberately accepted SHA-256 (`scripts/ledger_dimensions.accepted.sha256:1`;
  `scripts/check_ledger_dimensions_pin.py:88-98`).

  **Residual risk — do not upgrade that narrow fact into execution or correctness evidence:**
  - **Bytecode is outside the pin** (`scripts/check_ledger_dimensions_pin.py:88-98);
    `scripts/__pycache__/` is gitignored and unverified (`.gitignore:1`). CPython's
    import path truncates source mtime to an integer
    (`/usr/lib/python3.10/importlib/_bootstrap_external.py:973`) and supplies source size to timestamp
    validation (`/usr/lib/python3.10/importlib/_bootstrap_external.py:1000-1005`), which rejects a
    stored-mtime mismatch (`/usr/lib/python3.10/importlib/_bootstrap_external.py:637-640`) or a stored
    source-size mismatch (`/usr/lib/python3.10/importlib/_bootstrap_external.py:641-643`). **Same
    integer mtime and same source size are both required** for stale timestamp bytecode to remain
    valid. The measured equal-size `sum`→`min` rewrite within one second still reused bytecode, so
    deliberate header construction is not necessary
    (`_scratch/modpin/REMEDIATION_REPORT.md:34-43`); equal-length audit edits such as sign flips,
    `sum`→`min`, and `+`→`-` make the size condition routine and the hazard realistic here
    (`_scratch/modpin/REMEDIATION_REPORT.md:41-43`).
  - **The check does not cover its own execution environment.** In an isolated copy, a targeted
    `sitecustomize.py` on `PYTHONPATH` spoofing `hashlib.sha256` only for the mutated module bytes made
    the standalone control, all six converted-stage comparator CLIs, and the generator CLI exit 0.
    This is inherent to an in-process check.
  - **For the common trust root, the pin moves the single point of failure; it does not abolish it.**
    The authority has one entry naming only `ledger_dimensions.py`
    (`scripts/ledger_dimensions.accepted.sha256:1`), and every consumer delegates to the shared
    decision function (`scripts/check_ledger_dimensions_pin.py:75-98`). Editing that shared decision
    to return success compromises all consumers. There are **four** individual invocation sites:
    `scripts/compare_dimension_artifacts.py:251`,
    `scripts/generate_canonical_dimension_table.py:231` and `:461`, and
    `scripts/run_all_audits.sh:20`. Deleting one generator invocation leaves the other active; deleting
    the comparator or runner invocation compromises only that path, while the standalone control and
    the other validators still reject through the shared decision
    (`scripts/check_ledger_dimensions_pin.py:145-146`). This source-integrity hole is separate from the
    execution-environment interposition above.
  - **The module pin itself covers no stage `.py`, no `.wl`, no `mathematica/out/*.out`, and no
    sidecar record content or production provenance.** The separate source/sidecar/comparator controls
    still apply, and the forgeable-sidecar hole immediately below remains open: no validator executes
    the stage.
  - **The accepted digest is a bare trust root, not a correctness witness.** `--accept` moves the
    baseline without a reason field, signature, or second witness. In the measured stub replay,
    accepting the mutant made the comparator green again.

  ⛔ **`run_all_audits.sh` is therefore not a general validation gate.** Its pin invocation does
  propagate failure because line 20 is a bare command under `set -euo pipefail`, but the runner only
  tallies stage failures at line 134, then reaches EOF after the line-140 summary write with no
  `exit` or `return` derived from `$fail` (`scripts/run_all_audits.sh:134-140`), so it can report a
  non-zero `Fail:` count and still exit 0. It invokes the comparator and generator **zero times**, so the
  cross-engine validators still have no aggregate runner.
- ⛔⛔ **THE SIDECAR IS FORGEABLE — the `.py` side has no (g2).** `emit_dimension_sidecar` digests the
  **source bytes**, never the record content, and the comparator never executes the stage.
  **Demonstrated by execution:** mutate the `.py` so it declares wrong values (runs green, 82 PASS),
  then hand-write a `.wl`-matching sidecar carrying the *mutated* `.py`'s sha256 → the honest comparator
  says `FAIL mismatches=2`, the forged one says `PASS mismatches=0`, exit 0, while the `.py` on disk
  still declares `a_dim: Dim(2,0,0)`. The digest proves *"this sidecar names the `.py` on disk"*, never
  *"this sidecar was produced by running it."* ⇒ §4-g's *"run the stage first anyway and say that you
  did"* rests on an **unverifiable assertion** — exactly the trust problem (g2) fixes for the `.out`,
  still open on the Python side. ⭐ **Interim control: the orchestrator regenerates the sidecar itself
  before the (i) commit**, as done for stage016.
- ✅ **A `.py`, the shared dimension module, and a sidecar are source-hash bound.**
  `emit_dimension_sidecar` asserts separate SHA-256 digests of the stage source and
  `scripts/ledger_dimensions.py`; the comparator computes both current hashes, and the
  canonical-table generator delegates to that same check. They reject absence or disagreement.
  This closes the measured transposed-but-not-re-run stage018 PASS and detects a sidecar made stale by
  a shared-module edit. It does not detect an accepted shared-module mutation, execute the stage,
  measure source coverage, bind transitive dependencies, or replace the separate orchestrator control
  that regenerates each Mathematica `.out`.
- ⛔ **stage013's `m_shared` / `m_dims.*` is the shared dimension of the `M_AB` matrix entries, NOT
  the physical mass `m_GNLS`** (`.wl:452,465,468`, `.py:592-593`). Both render as `M`. ⛔ **Never
  D4-rename it to a generic mass name** — that merges different quantities, and no dimensional
  check could catch it. Same shape as the `c_s0`/`c_S` merge caught at stage018.
- `lru_cache` — **5 stages, 11 sites**: 018 (×1), 022 (×4), 023 (×1), 027 (×1), **040 (×4)**. 043 has
  **zero**; 040's four are verified argument-pure.
- `grep -c '^PASS'` **over-counts by exactly 1 only for the 27 scripts that emit a `PASS tally:`
  line**, because that tally self-matches. It does not over-count a `TALLY sympy:` stage: a fresh
  stage016 run produced 82 `^PASS` lines and `TALLY sympy: 82 pass + 0 fail = 82 checks`.
- Codex needs `--sandbox danger-full-access` for Mathematica · **≤2 concurrent kernels**.
- ⛔ **Never `git checkout`/`restore`/`stash` to undo an ablation** on uncommitted work — restores from
  HEAD and destroys it. Use `cp` backups; verify with `git hash-object`.
- ⛔ **After every ablation edit and its restore on a shared Python module, clear `scripts/__pycache__/`
  before the next execution.** Equal-size edit/run/restore loops can otherwise reuse a
  timestamp-valid stale `.pyc` from either side of the loop.
- ⛔ **Never state the reason you expect in a directive** — it gets adopted instead of checked. The
  stage012 waiver wrongly covered `omega_dim` because the directive supplied stage011's rationale.
  Ask for the determination **plus its evidence**, and point the review leg at every escape hatch.
- ⚠ Process probes lie: `pgrep -f 'math -script'` matched Codex's own argv (the directive quoted the
  string); `[ -d /proc/$PID ]` with an empty PID tests `/proc/` and always says "alive". **Check the
  artifact, not the process.**
- ✅ **NOT A HAZARD — `*.out` tracking is already covered, and the alarm was FALSE.** The 2026-07-28
  provenance survey reported, as an incidental find, that a bare `*.out` rule would silently untrack a
  *new* stage's Mathematica transcript. ⚠ **The mechanism has since changed and is re-established here
  against the current file:** there is no bare `*.out` rule left to negate — `.gitignore:19-22` says so
  in terms and `.gitignore:23` narrows the LaTeX ignore to `**/paper/*.out`, while `.gitignore:130-131`
  records that the former `!research/pde_ledger_v2/mathematica/out/*.out` negation was deleted as
  redundant. **Verified by execution, not by reading:** `git check-ignore -v` on the real path
  `research/pde_ledger_v2/mathematica/out/ledger_stage023_…_mathematica_audit.out` prints nothing and
  exits 1 — no rule matches it — and all 44 transcripts in that directory are tracked (44/44).
  ⭐ **The instructive part is the failure mode, which is why this bullet stays:** a `.gitignore` verdict
  is worthless without running `git check-ignore -v` on a real path; the last matching rule wins.

## 10. ⛔ WHAT NOT TO DO
- Do **not** generalise the module beyond the stage in front of you. stage011's expression walker is
  deliberately stage-local.
- Do **not** build a corpus-wide inventory, oracle, or completeness proof. Three such gates were
  specified and rejected; the per-stage loop replaces them.
- Do **not** repair tautologies during a refactor (stage031's 20-of-21 self-comparing rows, stage017's
  hardcoded `True`, stage004 `:302`/`:308`, stage013/018's literal prints, stage012's stale
  `source dict L522-L529` label). Catalogue them — task #22.
  ⭐ **stage013's hardcoded prints are now DEMONSTRATED stale-capable, not merely suspected.** Under a
  `Tw` ablation, `.py:588-589` / `.wl:448-449` printed `Tw=(1,1,-2)` and `K_eta=…=(-1,1,-2)` while the
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

⛔ **PROCESS RULE: a bare `python3 scripts/ledger_stageNNN_..._sympy_audit.py` run is a PRODUCER, not a
validator.** It writes a dimension sidecar; its exit code and PASS tally are **not validation
evidence**. The validators are the comparator, the canonical-table generator, and the module-pin
control. None of them upgrades the producer's own green tally into proof that the sidecar came from a
real run; the forgeable-sidecar hole in §9 remains.

⭐⭐ **THE CROSS-CHECK EARNS ITS KEEP — measured on stage016, 2026-07-27.** Relabelling the stage's
basis leaves its **own** assertions completely blind: under `("M","L","T")` it prints
`{'measure': 'M^3', 'M2_integral': 'L', …}` and still reports **82 PASS / 0 FAIL / exit 0**. The
comparator catches it — **18 of 21** mismatches under L↔M, 15/21 under L↔T, 16/21 under M↔T, 21/21 on
an axis rename. ⇒ **The comparator is the SOLE instrument standing between a converted stage and a
relabelled basis.** Waiving it for a stage removes the only detector, which is why
`ARTIFACT_NAME_WAIVERS` staying `{}` matters. The records it cannot see are exactly those whose
exponents are equal on the swapped axes.

⛔⛔ **AND HERE IS ITS CEILING — a RECONSTRUCTION-TO-LIVE-EXPRESSION linkage gap.** Six of the nine
computed records (`measure`, `m2_integral`, the three `k_*` terms, `k2_integral`) walk a **stage-local
dimensional reconstruction** built at `sympy:335-340`; the other three (`actual_M2`, `actual_K2`, and
their ratio) walk the **live** stage expressions `m2_expr` / `k2_expr` (`sympy:347-349`).
⚠ **No gate binds the reconstruction to the live physics it mirrors.** Executed on stage016, two
dimension-preserving edits *confined to that reconstruction* — dropping `lambda_m` from the angular
term, and deleting the angular term from the reconstructed `k2_integral` — each left **82 PASS /
exit 0, byte-identical `DIM|` records and a green comparator**. ⇒ The conversion and the comparator
cannot detect a reconstruction-vs-source mismatch; the comparator says so itself with
`source_coverage=not_checked`. **No current automated gate in this workstream closes that linkage
gap.** ⛔ Do not overstate this as "records never walk the stage's own expressions" — three of them do.

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
- `schemas/` + `schemas/validate_dimension_survey.py` are **parked** (survey-era, still committed).
- ⛔ **WORK-023-MOMENT-CONVENTION — adjudicate the pathA_28 source-moment identity and propagate it.**
  ⚠⚠ **INSTRUMENT INSUFFICIENT — the dimensional instrument below does not decide the remaining identity.**
  - **Dispute.** Stage006 records `CONVENTION P10` for projected `ρ_B=L⁻⁴` and
    `S_leak=L⁻⁴T⁻¹`
    (`scripts/ledger_stage006_two_phase_chiB_ontology_sympy_audit.py:391-397`), and stage008 records
    the spatial measure by defining `M0=∫_brane S_leak d³x` and the corresponding `D1`
    (`scripts/ledger_stage008_monopole_dipole_return_spec_sympy_audit.py:311-317`).
    Together those loci establish the number-density/`d³x` convention for the 006→008 route:
    `[M0]=L⁻¹T⁻¹` follows before any downstream declaration. Stage008 separately says that its
    dimensional block invents no `M0/D1/Q2` triples
    (`scripts/ledger_stage008_monopole_dipole_return_spec_sympy_audit.py:523-546`).
    Stage009 nevertheless declares `MOMENT0=T⁻¹`
    (`scripts/ledger_stage009_flat_slab_return_residual_sympy_audit.py:464-468`), stage010 declares
    its own `dim_M0=T⁻¹`
    (`scripts/ledger_stage010_slab_localization_p2_nogo_sympy_audit.py:553-562`), and stage023
    declares `[M0]=MT⁻¹`, `[D1]=MLT⁻¹`
    (`scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:460-465`,
    **the `SOURCED_DIMS` opening plus its `M0` and `D1` entries**).
    The stage023 note asserts “one intended object flowing 008 → 009/010 → 022/023”
    (`notes/stages/ledger_stage023_nullspace_underdetermination.md:521-527`, **§1.7(2), `M0`/`D1`**),
    but no locus cited in this entry derives that alias chain. The carrier/measure convention is
    therefore recorded upstream; the open item is whether those later names denote that same object.
    ⭐ **This is more decidable than the former convention dispute:** the recorded 006→008 route
    removes the carrier/measure choice, and a cross-stage provenance/dataflow derivation can test
    the remaining identity. The hazard has the `c_s0`/`c_S` shape
    (`notes/stages/ledger_stage016_l2_so3_covariance.md:187-191`, **(2) D4 name determinations,
    sharpest silent-merge hazard**), which a dimensional or cross-engine equality cannot detect.
  - **Loci on both sides.** Recorded number-density/`d³x` route:
    `scripts/ledger_stage006_two_phase_chiB_ontology_sympy_audit.py:391-397` and
    `scripts/ledger_stage008_monopole_dipole_return_spec_sympy_audit.py:311-317,523-546`;
    later declaration route:
    `scripts/ledger_stage009_flat_slab_return_residual_sympy_audit.py:464-468`,
    `scripts/ledger_stage010_slab_localization_p2_nogo_sympy_audit.py:553-562`, and
    `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:460-465`
    (**the `SOURCED_DIMS` opening plus its `M0` and `D1` entries**); asserted
    cross-stage identity: `notes/stages/ledger_stage023_nullspace_underdetermination.md:521-527`
    (**§1.7(2), `M0`/`D1`**).
  - **What would settle it.** Trace the definitions and actual handoffs across the asserted
    008→009/010→022/023 chain. If they are aliases, propagate the recorded 006→008 normalization
    through `M0/D1`, `R0/R1`, and `A0/A1` and correct every incompatible declaration; if they are
    distinct, record the separate definitions and remove the false identity. This replaces the
    former request to choose a carrier or measure: those inputs are already recorded upstream.
    The trace must include stage010's own declaration and the 022 handoff because the asserted chain
    names both; a trace restricted to stages 006, 008, 009, and 023 would not cover the identity as
    stated in the stage023 note.
  - ⚠⚠ **Instrument — INSUFFICIENT for the remaining same-object identity.** The dimensional
    instrument derives the consequences of stage006's recorded projected `S_leak` and stage008's
    recorded `d³x` moment
    (`scripts/ledger_stage006_two_phase_chiB_ontology_sympy_audit.py:391-397`;
    `scripts/ledger_stage008_monopole_dipole_return_spec_sympy_audit.py:311-317`).
    **What it does compute:** the upstream `[M0]`/`[D1]` route and its propagated `R0/R1/A0/A1`
    consequences, exposing which downstream declarations conflict if the asserted identity holds.
    **Why it does NOT decide:** dimensional propagation cannot prove that equal-looking names in
    stages 008, 009, 010, 022, and 023 are aliases; §4-(c1)'s naming rule requires physical
    provenance rather than dimensional equality. **What would be required to decide:** a
    cross-stage provenance/dataflow derivation showing the actual handoff or separate definitions.
    Because the carrier and measure are already recorded, `INSTRUMENT INSUFFICIENT` remains correct
    only for this dimensional instrument and the remaining identity; it is not a corpus-absence claim.
    A separate provenance/dataflow instrument is therefore named and potentially decisive; the
    insufficiency label does not classify the whole work item as undecidable from the repository.
- ⛔ **WORK-023-STAGE009-MOMENT0 — replace stage009's same-literal assertion with an upstream-derived check.**
  - **Dispute.** Stage009 declares `MOMENT0=T⁻¹`
    (`scripts/ledger_stage009_flat_slab_return_residual_sympy_audit.py:464-468`) for the `M0`
    supplied by stage008's moment definition
    (`scripts/ledger_stage008_monopole_dipole_return_spec_sympy_audit.py:311-317`), while stage023
    declares `[M0]=MT⁻¹`
    (`scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:460-464`,
    **the `SOURCED_DIMS` opening plus its `M0` entry**). In each
    stage009 engine, the only `M0`-bearing dimensional assertion constructs `[Z]` from the same
    `dimM0` literal and compares it back with that literal, so that assertion cannot reject any
    replacement triple in the scope of those two implementations
    (`scripts/ledger_stage009_flat_slab_return_residual_sympy_audit.py:467-496`;
    `mathematica/ledger_stage009_flat_slab_return_residual_mathematica_audit.wl:440-458`).
  - **Loci on both sides.** Stage009 declaration/check side:
    `scripts/ledger_stage009_flat_slab_return_residual_sympy_audit.py:467-496` and
    `mathematica/ledger_stage009_flat_slab_return_residual_mathematica_audit.wl:440-458`;
    upstream/downstream same-object side:
    `scripts/ledger_stage008_monopole_dipole_return_spec_sympy_audit.py:311-317` and
    `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:460-464`
    (**the `SOURCED_DIMS` opening plus its `M0` entry**).
  - **What would settle it.** Derive stage009's expected `[M0]` from the pathA_28 source definition
    without reading `MOMENT0`, then compare the declaration with that independently built result;
    §1.7 records the present defect as W2 at
    `notes/stages/ledger_stage023_nullspace_underdetermination.md:716`
    (**§1.7(5), the `W2` fragment, *"stage009 `MOMENT0 = T⁻¹` plus an assertion that cannot fail"***).
  - **Instrument.** A dual-engine stage009 computation whose expected branch consumes stage008's
    named `S_leak·d³x` route and whose ablation changes the source carrier or measure without changing
    the declaration (`scripts/ledger_stage008_monopole_dipole_return_spec_sympy_audit.py:311-317`;
    `scripts/ledger_stage009_flat_slab_return_residual_sympy_audit.py:474-496`).
- ⛔ **WORK-023-D0-SEAM — derive the 017→019 denominator normalization or separate the objects.**
  - **Dispute.** Stage017 defines and exports `D0=K₂−(B̃0+Z̃0)`
    (`notes/stages/ledger_stage017_grouped_p2_lane_isotropy.md:27-34`,
    **§1.1 grouped-lane assembly**;
    `notes/stages/ledger_stage017_grouped_p2_lane_isotropy.md:113-115`,
    **Exported — the ℓ=2 PORT KERNEL**), and the
    register records that D-lane at `MT⁻²` (`notes/parameter_register.md:185`). Stages 021, 023,
    and 027 declare the identified denominator `ML⁻¹T⁻²`
    (`scripts/ledger_stage021_dimensional_closure_sympy_audit.py:138-146`;
    `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:468`
    (**`SOURCED_DIMS[D0]`**);
    `scripts/ledger_stage027_port_checks_closure_sympy_audit.py:199-208`); stage023's port formula
    then makes the register reading incompatible with its dimensionless `P0_physical` target
    (`scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:220-232,486-493`,
    **`build_port_kernel` and `EXPECTED_DIMS["P0_physical"]`**).
  - **Loci on both sides.** Stage017/register side:
    `notes/stages/ledger_stage017_grouped_p2_lane_isotropy.md:27-34,113-115`
    (**§1.1 grouped-lane assembly; Exported — the ℓ=2 PORT KERNEL**) and
    `notes/parameter_register.md:185`; 021/023/027 closure side:
    `scripts/ledger_stage021_dimensional_closure_sympy_audit.py:138-146`,
    `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:468`
    (**`SOURCED_DIMS[D0]`**), and
    `scripts/ledger_stage027_port_checks_closure_sympy_audit.py:199-208`.
  - **What would settle it.** Derive an explicit 017→019 normalization carrying the missing `L⁻¹`,
    derive that the two `D0` names denote distinct objects, or identify and correct the wrong member
    of `[D0]`, `[N0_from_port]`, and `[P0_physical]`; these three settlement routes are recorded at
    `notes/stages/ledger_stage023_nullspace_underdetermination.md:474-476`
    (**§1.7(1), the `D0` two-readings note, *"Settled by:"* bullet**).
  - **Instrument.** A cross-stage, dual-engine derivation that carries the denominator from named
    upstream stages 017 and 019 into the 021/023/027 closure and checks the independently derived
    `P0_physical` target (`notes/stages/ledger_stage017_grouped_p2_lane_isotropy.md:113-115`,
    **Exported — the ℓ=2 PORT KERNEL**;
    `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:220-232,486-493`,
    **`build_port_kernel` and `EXPECTED_DIMS["P0_physical"]`**).
- ⛔ **WORK-023-STIFFNESS-REDUCTION — execute the scalar reduction that W4's two loci describe differently.**
  ⚠⚠ **INSTRUMENT INSUFFICIENT — the instrument named below reaches the ℓ=0 half only, not `K1`.**
  - **Dispute.** The stage023 source map states the stage013/017→023 `K0c/K1` identification as
    performed (`notes/stage023_pathA34_nullspace_underdetermination_source_map.md:250-253`), while
    the register records the same reduction as `FREE-UNREDUCED`/`PENDING`
    (`notes/parameter_register.md:170`); §1.7 classifies the source-map statement as stale pre-build
    rather than a live result (`notes/stages/ledger_stage023_nullspace_underdetermination.md:719-723`,
    **§1.7(5), W4**).
  - **Loci on both sides.** Performed-identification side:
    `notes/stage023_pathA34_nullspace_underdetermination_source_map.md:250-253`; pending-reduction
    side: `notes/parameter_register.md:170`.
  - **What would settle it.** Execute the profile-and-measure scalar reduction of stage013's
    `K_AB` collective and stage017's harmonic sector to stage023's `K0c` and `K1`; the register
    states that discharge route at `notes/parameter_register.md:170`.
  - ⚠⚠ **Instrument — INSUFFICIENT for the `K1` half of this dispute.** The obvious instrument is a
    dual-engine reduction cross-checked against stage013's projected `K_AB` definition
    (`notes/stages/ledger_stage013_breathing_harmonic_mk_projection.md:75-82`,
    **M_AB / K_AB by real ∫dw operator projection**), stage016's measure-bearing stiffness assembly
    (`notes/stages/ledger_stage016_l2_so3_covariance.md:74-81`, **1.4 Angular dimensional
    consistency**), and stage017's exported grouped lane
    (`notes/stages/ledger_stage017_grouped_p2_lane_isotropy.md:27-34`, **§1.1 grouped-lane
    assembly**). **What it does compute** — worth having: stage013's genuine `∫dw` operator
    projection of the collective `K_AB` at `M T⁻²`, which is the ℓ=0 route to `K0c` and IS in reach
    of these named inputs; and stage016/017's reduced pair `K̃ + λ_m·T̃_Ω` at `λ_m = 6`, which is the
    **ℓ=2** assembly. ⛔ **Why it does NOT decide `K1`:** stage023's `K1 = K_eta + 2·T_Omega` is that
    same assembly at `λ_m = ℓ(ℓ+1) = 2`, i.e. **ℓ=1**, while the stage016/017 loci named here reduce
    against the frozen **ℓ=2** profile `β₂`
    (`notes/stages/ledger_stage023_nullspace_underdetermination.md:496-512`, **§1.7(2), merge to
    refuse**) — the ℓ=1 profile the `K1` reduction consumes is not among these inputs.
    ⭐ **What would be required to decide it:** a derived ℓ=1 radial profile `β₁`, from the same wall
    operator and boundary data that produced `β₂`. ⛔ **No derived or defined ℓ=1 radial profile
    `β₁` is recorded in the four sources stage023's leg searched** —
    `research/pde_ledger_v2/notes/` (incl. `stages/`), `research/pde_ledger_v2/paper/stages/`,
    `research/pde_ledger_v2/manifests/`, and `docs/model_map.md`; within that scope, `β₁` occurrences
    are meta-level open-question discussions or stage029's unrelated `β₁PN`, not a profile derivation
    or definition (`notes/stages/ledger_stage023_nullspace_underdetermination.md:496-512`,
    **§1.7(2), merge to refuse**). **So the `K1` half cannot be settled from those four sources as
    they currently stand** — it is gated on **WORK-023-L1-L2-PROFILE-IDENTITY** first producing
    `β₁`, and of the two halves only the ℓ=0 `K0c` route is executable from the loci named here.
- ⛔ **WORK-023-L1-L2-PROFILE-IDENTITY — derive whether the reduced ℓ=1 and ℓ=2 stiffnesses coincide.**
  - **Dispute.** Stage016's reduced scalars use the frozen `β₂` profile
    (`notes/stages/ledger_stage016_l2_so3_covariance.md:74-81,130-131`), while stage023's `K1`
    is the ℓ=1 assembly (`notes/stages/ledger_stage023_nullspace_underdetermination.md:496-512`,
    **§1.7(2), merge to refuse**). §1.7 leaves the same-number question undetermined pending a
    derivation of `β₁≡β₂` (`notes/stages/ledger_stage023_nullspace_underdetermination.md:496-512`,
    **§1.7(2), merge to refuse**).
  - **Loci on both sides.** ℓ=2 side:
    `notes/stages/ledger_stage016_l2_so3_covariance.md:74-81,130-131`; ℓ=1 side:
    `notes/stages/ledger_stage023_nullspace_underdetermination.md:496-512`
    (**§1.7(2), merge to refuse**).
  - **What would settle it.** Derive the ℓ=1 and ℓ=2 radial profiles from the same wall operator
    and boundary data, then prove `β₁≡β₂` or exhibit the derived distinction
    (`notes/stages/ledger_stage023_nullspace_underdetermination.md:496-512`,
    **§1.7(2), merge to refuse**).
  - **Instrument.** A dual-engine radial eigenproblem with separately represented `β₁` and `β₂`,
    solved from one specified wall operator, one set of boundary data and one normalization, which
    then **evaluates and compares the reduced scalar integrals** `K̃ = ∫[T_w β'² + K_η β²]dV` and
    `T̃_Ω = ∫T_Ω β² dV` at ℓ=1 against ℓ=2 — those reduced scalars, not the profiles, are what `K1`
    and stage016's `K₂` consume, and comparing the profiles alone would not decide the dispute
    because `β₁ ≠ β₂` does not by itself imply unequal reduced scalars
    (`notes/stages/ledger_stage016_l2_so3_covariance.md:74-81`, **1.4 Angular dimensional
    consistency**; `notes/stages/ledger_stage023_nullspace_underdetermination.md:496-512`,
    **§1.7(2), merge to refuse**).
- ⛔ **WORK-023-CS-EVALUATION — derive the state point of stage023's sound-speed carrier.**
  ⚠⚠ **INSTRUMENT INSUFFICIENT — the instrument named below does not decide this dispute.**
  - **Dispute.** Stage005 derives the state-dependent `c_s(ρ)` and names `c_s0` as its asymptote
    (`notes/stages/ledger_stage005_sound_speed_light_ratio.md:74-104`,
    **§1, EOS derivation; §2, velocity scales**); the register records that asymptotic carrier
    (`notes/parameter_register.md:129`). Stage023 consumes `c_s` in `z=aω/c_s` while §1.7 finds
    no evaluation point in the four sources it searched
    (`notes/stages/ledger_stage023_nullspace_underdetermination.md:513-520`,
    **§1.7(2), c_s forward hazard**).
  - **Loci on both sides.** Asymptotic `c_s0` side:
    `notes/stages/ledger_stage005_sound_speed_light_ratio.md:74-104`
    (**§1, EOS derivation; §2, velocity scales**) and `notes/parameter_register.md:129`;
    stage023 carrier side:
    `notes/stages/ledger_stage023_nullspace_underdetermination.md:513-520`
    (**§1.7(2), c_s forward hazard**) and
    `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:211-222,333-338`.
  - **What would settle it.** Derive the density state sampled by the stage023 outgoing-wave
    carrier; evaluation at `ρ0` would establish the `c_s0` identification, while a different
    derived state would establish distinct carriers
    (`notes/stages/ledger_stage005_sound_speed_light_ratio.md:74-104`,
    **§1, EOS derivation; §2, velocity scales**).
  - ⚠⚠ **Instrument — INSUFFICIENT for this dispute.** The obvious instrument is a cross-check
    against named upstream stage005 that carries an explicit state argument through a dual-engine
    reconstruction of stage023's `z` and residual amplitudes
    (`notes/stages/ledger_stage005_sound_speed_light_ratio.md:74-104`, **§1, EOS derivation; §2,
    velocity scales**; `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:211-222,333-338`).
    **What it does compute** — worth having: it re-expresses `z = aω/c_s` and the amplitudes
    `A0 = i·v₀·z·M0·(1−T0)`, `A1 = i·v₁·z³·D1·(1−T1)` as explicit functions of `ρ` through stage005's
    derived `c_s²(ρ) = 5Kρ⁴/m_GNLS` and its log-slope `d ln c_s/d ln ρ = 2`, making every stage023
    consumer's sensitivity to a shift of state point computable rather than assumed.
    ⛔ **Why it does NOT decide:** making the ρ-dependence explicit is not a selection of ρ —
    stage005 derives the *function* `c_s(ρ)` and names `c_s0` only as its asymptote, so neither that
    derivation nor the reconstruction states which `ρ` the stage023 outgoing carrier samples.
    ⭐ **What would be required to decide it:** a background/boundary condition for the **outgoing**
    wave — a stated far-field or matching condition fixing the density the radiated carrier
    propagates through (`ρ0` would establish the `c_s0` identification; a derived wall-side state
    would establish a distinct carrier). ⛔ **No such condition appears in the four sources §1.7
    searched** — the stage023 note's §4, the two stage023 engines, `docs/model_map.md:62`, and the
    register — and that set of four is the scope of this negative
    (`notes/stages/ledger_stage023_nullspace_underdetermination.md:513-520`, **§1.7(2), c_s forward
    hazard**). **So this cannot be settled from those four sources as they currently stand:** it
    needs the radiation background itself written down, as a derivation or as a recorded convention.
- ⛔ **WORK-023-SOURCED-PROVENANCE — trace the declaration lineage before using “sourced.”**
  - **Dispute.** The stage023 Mathematica artifact emits `sourced_dims.M0/D1`
    (`mathematica/ledger_stage023_nullspace_underdetermination_mathematica_audit.wl:298-310`)
    and the Python engine stores them in `SOURCED_DIMS`
    (`scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:460-483`,
    **the whole `SOURCED_DIMS` block, all 22 entries**), while named
    upstream stage008 explicitly declines to assign dimensions to `M0/D1/Q2`
    (`scripts/ledger_stage008_monopole_dipole_return_spec_sympy_audit.py:523-546`). §1.7 records
    that mismatch as W5 (`notes/stages/ledger_stage023_nullspace_underdetermination.md:723-725`,
    **§1.7(5), W5**).
  - **Loci on both sides.** “Sourced” side:
    `mathematica/ledger_stage023_nullspace_underdetermination_mathematica_audit.wl:298-310` and
    `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:460-483`
    (**the whole `SOURCED_DIMS` block, all 22 entries**); no-upstream-
    dimension side:
    `scripts/ledger_stage008_monopole_dipole_return_spec_sympy_audit.py:523-546`.
  - **What would settle it.** Trace each stage023 declaration to a named upstream dimension
    derivation and classify entries without such a derivation as stage-local declarations; the
    present term means only “declared in this stage's own dict”
    (`notes/stages/ledger_stage023_nullspace_underdetermination.md:723-725`,
    **§1.7(5), W5**).
  - **Instrument.** A per-record provenance cross-check against named upstream stage008 plus both
    stage023 declaration tables, with an able-to-fail fixture that removes or redirects an upstream
    locus (`scripts/ledger_stage008_monopole_dipole_return_spec_sympy_audit.py:523-546`;
    `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:460-483`
    (**the whole `SOURCED_DIMS` block, all 22 entries**);
    `mathematica/ledger_stage023_nullspace_underdetermination_mathematica_audit.wl:298-310`).
- ⛔ **TRACKED SPEC — exact-rational artifact syntax in
  `scripts/compare_dimension_artifacts.py`.** This is a validator change, not per-stage work, so it
  belongs here. Acceptance: (1) before `Fraction` normalization, accept only stripped signed integer
  or signed `numerator/denominator` tokens and keep zero-denominator rejection; (2) reject decimal
  points and scientific notation, including `0.5`, `1.0e0`, and `2.5`; (3) prove `-2`, `1/2`, and
  `-7/2` still load; (4) make both the comparator CLI and canonical-table generation fail nonzero on
  the bad-token fixture; (5) keep explicit non-empty record floors in both fixtures; (6) re-run every
  converted-stage comparator and regenerate the canonical table. **Relation to the existing tracked
  spec below:** both change the validator layer, but this one governs exponent-token lexical
  exactness while the existing spec governs cross-engine axis-order metadata. They are orthogonal;
  neither acceptance list subsumes or restates the other.
- ⛔ **TRACKED SPEC — Mathematica transcript liveness in a new
  `scripts/check_mathematica_transcript.py`.** This is likewise a harness/consumer change, not a
  per-stage conversion. Acceptance: (1) take the generated transcript and committed reference as
  explicit paths; (2) print byte and line counts for both and reject either artifact at zero bytes or
  zero lines; (3) apply the existing kernel-ID normalization, then require byte identity; (4) return
  nonzero on a missing, empty, truncated, or unequal generated transcript regardless of the producer's
  exit status; (5) fixtures prove exact non-empty PASS, empty-generated FAIL at the first byte,
  empty-reference FAIL, and truncated-prefix FAIL; (6) wire the named check into §4-(d)/(g2) without
  weakening the existing manual control. This instruments the control that already catches the
  hazard; it does not repair the six catalogued `Throw`/`Catch` source structures.
- ⛔ **TRACKED SPEC — the canonical-table generator's axis-order invariant must be REPLACED, not
  deleted, before 035/036/037.** It currently *raises* on cross-engine axis-order divergence
  (`generate_canonical_dimension_table.py:255-265`), which is exactly the 037 route (Wolfram `M,L,T`
  vs Python `(L,T,M)` — `notes/stage037_dimension_emission_spike.md:79-85`), so it will crash the
  toolchain. ⚠ Deleting the raise alone is **refuted**: it leaves a non-deterministic
  `next(iter(axis_orders))` selection, and the raise doubles as the **zero-record floor**.
  Acceptance: (1) keep per-file validation in `load_dimensions`; (2) store explicit `python_axes` +
  `wolfram_axes`, engine-tagged; (3) keep an EXPLICIT non-empty check rather than relying on
  `len(axis_orders)!=1`; (4) render stage coverage by the same rule as rows — one tuple when equal,
  `py …; wl …` when different; (5) when orders differ, render **both** positional exponent tuples
  (`row_exponents_text` currently collapses to Python's with no indication); (6) prose "each stage's
  declared order" → "each **engine's** declared order"; (7) tests for divergent-but-label-equivalent,
  divergent mismatch, empty artifact, deterministic coverage; (8) keep semantic validation at the
  emitters + review layer — deriving `axes=` from the live basis is a materially better guard than
  cross-engine equality. ⚠ The comparator pairs by axis **label** and never compares the two
  artifacts' declared axis **order**, so a reordered-but-relabelled sidecar passes today.

### 12b. ⭐ PRIOR-ART SURVEY, 2026-07-28 — three decisions, all "build the small thing"

Three read-only surveys asked whether this workstream is reinventing solved problems: a Python
dimensional-analysis library instead of `scripts/ledger_dimensions.py`; an off-the-shelf mutation-testing
tool instead of the hand-rolled per-tooth ablation; and an attestation scheme for §9's forgeable-sidecar
hole. ⭐ **All three returned "do not adopt" — each for a different and specific reason, and each leaves
a concrete project-owned deliverable in its place.** ⚠ The survey called all three *small*; (b) below
records that for the ablation driver that sizing proved low by a large factor.
Nothing was installed and no corpus file was modified. The
working notes are in `_scratch/prior_art/` (`SURVEY_dimension_libs.md`, `SURVEY_mutation_testing.md`,
`SURVEY_provenance.md`) — ⚠ they **still exist**, but they are gitignored, so they are not in git and are
not guaranteed to survive; what is load-bearing is recorded here.

**(a) Dimension libraries — NO-GO, keep `ledger_dimensions.py`.**
- ⭐ **`sympy.physics.units` discards the declared axis order.** `DimensionSystem` sorts its bases:
  declared `('M','L','T')`, `('L','T','M')` and `('L','M','T')` all collapse to `(L, M, T)`, and
  `dim_vector` is likewise sorted. stage004's committed sidecar declares `axes=L,T,M`, which the library
  **cannot produce** — the order must be re-imposed by a wrapper the library does not supply. That is
  decisive here and nowhere else: axis-labelled comparison is this workstream's **sole** transposition
  detector (§11).
- Four further measured regressions against the incumbent: **zero-exponent axes are dropped**
  (`get_dimensional_dependencies(L)` → `{L: 1}`, not a full-length vector, while
  `compare_dimension_artifacts.build_dimension` raises on a length mismatch); **no membership
  validation** — `get_dimensional_dependencies(Dimension("typo_axis"))` → `{typo_axis: 1}` with no
  error, where `basis.from_mapping` raises `missing=/extra=`, so a misspelled axis would silently pass
  and silently change the answer; **Float exponents silently accepted** — `L**0.5` →
  `{Dimension(L): 0.500000000000000}`, where `_exact` coerces to `sp.Rational` at construction so a
  Float can never enter; and **the walker keys off `Quantity`, not `Symbol`** —
  `SI._collect_factor_and_dimension(sp.Symbol('x'))` → `(x, Dimension(1))`, i.e. an **unmapped symbol is
  silently DIMENSIONLESS**, where each stage's `dim_of` raises
  `DimError("missing dimension for symbol …")`. (Also: `Dimension("L")` is a global singleton shared
  across bases, and `Dimension("mass") == ` SymPy's SI mass is `True`.)
- **`astropy.units` manufactures the exact defect this project treats as fatal, by design.**
  `astropy/units/utils.py sanitize_power` converts any power whose denominator is a power of two to a
  float: `elif (denom & (denom - 1)) == 0: p = float(p)`. Denominator 2 is a power of two, so `1/2` and
  `3/2` become `0.5` and `1.5` — hitting `probe_stage042`'s
  `charge = basis(Fraction(1,2), Fraction(3,2), Fraction(-1))` (which then asserts that no `sp.Float` is
  present) and stage004's `sp.Rational(1, 2)` sound-speed/healing-length powers. The library's
  normalisation *is* the defect, silently, below any wrapper.
- **`unyt` cannot register a new base dimension.** `unyt/dimensions.py`'s `base_dimensions` is a
  hardcoded module-level list of 9 SymPy Symbols; `define_unit`/`UnitRegistry.add` register *units* over
  existing bases. Adding `stiffness` means patching a third-party global — strictly worse than 222 pinned
  lines for a control that is supposed to be byte-reproducible.
- **`numericalunits` is non-deterministic and incompatible with SHA-pinned artifacts.** There is no
  dimension object at all (*"you cannot directly see what the units are. You are supposed to already know
  what the units are"*); base units are randomly-chosen floats, so detection requires running the
  calculation **twice in separate Python sessions** and comparing.
- **`pint` is plausible but was NOT EXECUTED.** Non-physical and extra base dimensions are documented as
  supported (*"primary dimensions don't need to be declared; they can be defined for the first time as
  part of a unit definition"*), but the default `UnitRegistry(non_int_type=float)` turns an exponent of
  `1/2` into float `0.5`, and the fix (`non_int_type=Fraction`) is **registry-wide** and drags quantity
  magnitudes with it; no declared axis order, no zero-padding, no SymPy interop, and a new unpinned
  dependency.
- ⭐ **The reframing worth keeping: the incumbent's edge is axis-order preservation and STRICTNESS, not
  expressiveness.** Both exotic bases *pass* on sympy — the 4-axis `(M,L,T,E)` with `E` genuinely
  independent (`mu_r = M²·L·T⁻⁴·E⁻¹`) and the non-physical `(stiffness,L,T)` with
  `charge = S^(1/2)·L^(3/2)·T⁻¹` both build, with exponent types `Rational`/`Half`/`NegativeOne` and zero
  `Float` atoms. The gap is never "it cannot express our basis"; it is that it will not keep the declared
  order and will not refuse bad input.
- ⚠ **Scope of the verdicts.** Only `sympy` (1.14.0, on CPython 3.10.12) is installed and installation
  was forbidden, so the `pint`/`unyt`/`astropy`/`numericalunits` verdicts rest on upstream source and
  documentation (`pint/util.py`, `pint/registry.py`, `unyt/dimensions.py`, `astropy/units/utils.py`),
  **not on execution**. Migration cost if this is ever reopened: 7 basis declarations, ≈214 dimension
  call sites, 7 sidecar regenerations (each header embeds `ledger_dimensions_sha256`, so all invalidate
  at once), a pin re-baseline, matching edits to the 7 hand-rolled `.wl` emitters — and a trust-root
  change, because the SHA would then cover only the wrapper while the library version becomes an
  **unpinned input the sidecar does not attest**.

**(b) Mutation testing — do not adopt a tool; build the ablation driver.**
⭐⭐ **NO LONGER A CANDIDATE — USER DECISION, 2026-07-29; the decision and its constraints are the last
four bullets of this subsection.**
⚠ **The "~100 lines" this survey recorded was an estimate of the hand-rolled loops being *replaced*, not
a size budget for the deliverable.** The contract since authored against the accepted requirements is
substantially larger, because exact resume, crash repair, an outcome truth table and a restore proof are
all required behaviour. The estimate was low by a large factor — which is itself worth a reader knowing,
and is why it is kept here rather than deleted.
- **`mutmut` writes SCHEMATA — the decisive blocker, because it makes a control pass vacuously rather
  than error.** Its `ARCHITECTURE.rst`: *"The mutated files contains the original code and the mutants.
  With the `MUTANT_UNDER_TEST` environment variable, we can specify (among other things) which mutant
  should be enabled."* Every mutant is therefore **byte-identical on disk**, so `emit_dimension_sidecar`'s
  `source_sha256` is a constant across all mutants and the freshness pin `require_fresh_stage_source`
  *"degenerates into 'the schemata file equals itself'"* — a control that silently stops being a control.
  ⚠ Its sibling fails the other way, and the two must not be collapsed: if `ledger_dimensions.py` falls
  inside `source_paths`, the schemata-rewritten bytes never equal `ledger_dimensions.accepted.sha256`, so
  every mutant becomes a **uniform false kill** (`COMPARISON_SKIPPED: ledger-dimensions module pin
  failed`) that tells you nothing. mutmut is additionally pytest-only since 3.0 with no runner key, and
  copies the source into `mutants/`, breaking every `__file__`-derived path.
- **`mutatest` mutates `__pycache__`, and CPython never caches the module run as `__main__` — so it is a
  silent no-op here.** Verified locally: running `python3 pcachetest/s.py` produced **no** `__pycache__`;
  a subsequent `import s` produced one. Every mutant would be reported "survived" while nothing was ever
  mutated — a uniform false negative. (`MutPy` fails similarly for a different reason: in-memory AST
  injection through the import machinery cannot reach a script run as a separate process, and it is
  unmaintained since 2019-11.)
- ⭐ **NO tool has a per-mutant hook to reset an emitted SIDE-EFFECT artifact.** cosmic-ray's
  `use_mutation` puts *"the unmutated **code** back in place"* and nothing else; mutmut's `mutants/` tree
  is per-run, not per-mutant, so mutant *n*'s artifact is visible to mutant *n+1*; only universalmutator
  handles it, and by **convention in the test command** (its own README example already `rm`s a generated
  artifact before each trial), not by a hook. This is the failure that actually bit: of the 22
  `cmp_*.txt` banked from the real stage023 ablation, **16 are the 673-byte freshness-failure text**
  (`COMPARISON_SKIPPED: Python dimension sidecar freshness failed`) and only 6 are a real `MISMATCH` —
  i.e. for 16 of 22 mutants the comparator's verdict was produced by the residual `.dimensions.txt`, not
  by a dimension comparison at all.
- **Targeting is subtractive or random everywhere, with no include-list.** cosmic-ray's three shipped
  filters (`cr-filter-operators`, `cr-filter-pragma`, `cr-filter-git`) all *skip* work items and there is
  **no "mutate only these lines" input**; mutmut offers file globs, pragmas and function-name globs (and
  `SOURCED_DIMS` is a module-level dict literal, so the function glob cannot reach it); MutPy's
  `--percentage` and mutatest's `--nlocations` select **random** sites. ⚠ universalmutator's line
  restriction is **unverified**. That collides head-on with `docs/development_pipeline.md:330-332` — the
  ablation target list is the **orchestrator's**, and any tool whose `init` enumerates the sites has
  taken that ownership. Scale: `cosmic-ray init` on the stage023 file would emit order-10³ work items
  against the wanted 22 + 29 = 51, so a custom filter would have to delete ~95%.
- **The bind ablation is not expressible by any shipped operator.** `bind.tsv`'s
  `SOURCED_DIMS[Z1ret] → SOURCED_DIMS[a]` is an identifier-for-identifier swap; cosmic-ray's
  `VariableReplacer` substitutes `Number(value=str(randint(-100, 100)))` — a **random integer**, and
  non-deterministic across runs, so a banked TSV row could not be reproduced.
- **cosmic-ray is the only serious candidate** (verdict PARTIAL; universalmutator is also PARTIAL, on the
  weakest evidence base of the five). Adopting it still requires a custom include-filter over
  `session.sqlite`, a *packaged* custom operator for the rebinding, a bespoke SQLite reader to rebuild the
  five TSV columns, a `bash -c` wrapper (`test-command` is `shlex.split`, not shelled), and the sidecar
  reset inside that wrapper — *"Items 2–5 are more code than the thing being replaced."*
- ⭐ **The one non-negotiable design constraint it establishes: MUTATE ON DISK, AT THE REAL PATH.**
  cosmic-ray does (*"applies a mutation to a file on disk, and after the with-block it put the unmutated
  code back in place"*) and the survey calls that *"the single biggest compatibility fact"* — because
  this project's controls hash the file's own bytes, `__file__`, the sidecar's destination path and the
  `sha256` freshness pin only stay honest when the real file at the real path is what changes.
- ⭐ **The three features the survey found missing from the hand-rolled loops:** (1) per-mutant reset
  of the source **and** of every emitted artifact, as a first-class step rather than something smuggled
  into a wrapper; (2) per-mutant capture retained, never overwritten; (3) resume after interrupt, by
  skipping rows already banked. Its input is an explicit **orchestrator-supplied include-list** of
  `(name, line, old_text, new_text)` rows, which satisfies the process clause by construction. Restore
  still proves itself by `cp` + `git hash-object` (§9), never `git checkout`. ⚠ **These three are a
  floor, not the requirement set** — the requirements since accepted against this decision add exact
  resume, crash repair, an outcome truth table and a restore proof (see the standing-position bullet
  below), which is where the size estimate went.
- ⭐⭐ **USER DECISION, 2026-07-29 — the ablation driver LEADS the BUILD QUEUE (which tooling gets built
  next), ahead of the shared Mathematica `DIM|` emitter, and it is not per-stage throwaway.** ⚠ **Two
  different sequences run here and this decision sets only the second.** The **conversion order** (which
  stage converts next) is §8's, and it is unchanged: **027** is still the next conversion. The **build
  queue** is this one, and it reads **driver first, emitter second** — so in time both are built
  **before** the 027 conversion begins.
  The existing hand-rolled ablation harnesses are to be switched over to it, and it is to be **reused
  for every subsequent stage rather than re-written per stage**. This is a commitment to build; the two
  remaining survey outcomes are not (see the closing paragraph of §12b).
- ⛔ **It is acceptance tooling, so it cannot grade itself.** Per `docs/development_pipeline.md:319-322`
  (**checklist item 1b — "Every checker named in the directive already existed before the build"**),
  when the deliverable *is* the checker, the positive and negative conformance fixtures and their
  expected outcomes must be authored and **frozen by a different session** before the build starts; and
  per this workstream's own floor rule (§4, `60e7032c` — *sidecar + comparator with a non-empty floor*),
  their **non-empty floors** are frozen with them. The building session may neither author nor weaken
  any of the three. ✅ **This rule has been SATISFIED for the driver, not merely stated** — see the
  step-2-done bullet below; it remains a standing rule for the `DIM|` emitter and anything after it.
- ⭐ **The retrofit's able-to-fail acceptance test has its ANSWER, and — since 2026-07-29 — its CHECKER
  too.** ⚠ Say which of the two, because "already available" was read both ways: the **answer** is the
  **committed oracle** (`notes/stage023_step_h_evidence/results.tsv`, in git) **and the legacy→new field
  mapping**, frozen in the driver's contract (`notes/ablation_driver/CONTRACT.md` §C-9). The **executable
  checker** no longer has to be written: it is item A7 of the frozen suite
  (`notes/ablation_driver/fixtures_v4/run_conformance.py` with `fixture_oracle.py`), which performs the
  live join against those committed legacy tables rather than trusting a precomputed digest — and, per the
  bullet above, it was authored by a session that is not the building one. ⚠ The phrase's original point
  still holds and is the durable one: **nobody gets to decide the right answer after the fact.** The
  comparison being automated came later and does not replace that.
  Re-running stage023's two ablation axes through the new driver must reproduce the
  committed `notes/stage023_step_h_evidence/results.tsv` — **22 `A1_DECLARATION` rows and 29
  `A2_BINDING` rows, 51 data rows under one 13-column header** — **exactly on every mapped field**.
  ⚠ **Not byte equality, which was never achievable:** that table carries the old hand-rolled harness's
  schema (`stage_exit`, `pass_count`/`fail_count`, `first_fail`, `sidecar_written`, and no outcome
  column), so a successor with its own result schema cannot be byte-equal to it. The **mapping** from
  legacy column to new field is part of the driver's **contract**, frozen with it rather than chosen at
  comparison time. ⛔ This is not a relaxation to "broadly matches" — exactness on the mapped fields is
  the point, and a disagreement on one is a **finding to report**, not a difference to reconcile. A
  driver that cannot reproduce those tables under that mapping has not replaced the harness it claims to
  replace, and the failure is visible without anyone deciding what the right answer was.
- ⭐ **A fourth required feature the list above does not state, learned this stage: whatever the driver
  emits must be written to be read from the COMMIT, not from the scratch directory it ran in.** Every
  path, count and row reference in its output must resolve against the committed tree; run-local scratch
  paths, and counts true only of the run, are not acceptable output. Four separate claims in this
  stage's hand-written evidence summary (`notes/stage023_step_h_evidence/ABLATION_SUMMARY.md`) were true
  of the run and false of the commit, and each had to be corrected under review — so this is a
  requirement on the driver's output format, where that class of error is designed out, not a reviewer's
  checklist item, where it is merely caught.
- ⭐ **Where the driver stands as of 2026-07-29, so the next session need not re-derive it.** The
  specifications are **committed** under `research/pde_ledger_v2/notes/ablation_driver/` — ⚠ **cite those
  paths, not the gitignored `_scratch/` originals** (`.gitignore:96`, the `research/**/_scratch/` rule),
  which will not survive the session.
  The requirements are at **v4** in `notes/ablation_driver/REQUIREMENTS.md:1` (heading **"Ablation
  driver — requirements and acceptance intent (v4)"**); **R8's original property** — that
  committed evidence alone lets a reader re-run a row — was **withdrawn by user decision on scope**
  (`notes/ablation_driver/REQUIREMENTS.md:8`, **"v4 — USER DECISION 2026-07-29: scope down, R8 is
  dropped"**), and R8 now states the narrower prose rule **"the evidence claims exactly what it
  supports"** (`notes/ablation_driver/REQUIREMENTS.md:108`, **R8**); and a **contract has been
  authored against v4** (`notes/ablation_driver/CONTRACT.md:5`,
  **"Applies to: `REQUIREMENTS.md` v4, 2026-07-29"**), alongside `CONTRACT_NOTES.md` /
  `CONTRACT_NOTES_V4.md` and `FIXTURES_REPORT.md` (⛔ the last describes a **superseded v1 suite** — see
  its status banner).
- ⭐⭐ **STEP 2 OF THE THREE-SESSION SHAPE IS DONE, 2026-07-29 — the fixtures are FROZEN AND COMMITTED.**
  The suite is `notes/ablation_driver/fixtures_v4/` and its byte authority is **outside** it, at
  `notes/ablation_driver/fixtures_v4.accepted.sha256`; `fixtures_v4/verify_freeze.py` reports **36 governed
  paths** and exits 0. **What is earned:** all nine acceptance items **A1–A9** pass a driver written from
  `CONTRACT.md` alone, A7 over the complete real stage023 row set; every rejection was probed to fire on its
  **named** property rather than on an incidental difference; and non-interference survives **SIGTERM and
  SIGKILL** at mutation boundaries.
  ⚠⚠ **THREE INDEPENDENT FULL-LADDER PASSES EXIST, BUT NOT ON ONE SET OF BYTES — say it this way, because
  the aggregate count reads as stronger than the provenance is.** **Two** passes are recorded at the freeze
  commit `abcb7f2b`: two independent reviewers each wrote a driver from `CONTRACT.md` alone and ran the real
  ladder, all nine items passing with A7 over the complete real row set (~25 min each). `4c22872a` then
  changed **two governed files** — `REQUIREMENTS.md` (+69 lines) and `run_conformance.py` (126 lines
  touched) — so those two passes do **not** cover the suite as it now stands. **One** pass is recorded on
  the bytes `41b66dd5` accepted: one independent party, a driver from `CONTRACT.md` alone, the full ladder
  *on those exact bytes*, all nine items, **51/51 A7 rows**, zero spurious wait trips. `0803392f` then
  changed `run_conformance.py` again (3 insertions, 2 deletions) and was verified only as far as **reaching
  driver invocation from a clean clone**, on the argument that the change is inert when the scratch
  directory already exists; `96a1a61b` re-accepted on that argument. ⇒ **Established: one full-ladder pass
  at the `41b66dd5` byte-state, two at `abcb7f2b`. OUTSTANDING: any full-ladder pass on current `HEAD`
  bytes.** That argument is recorded here, not adjudicated; a reader who wants the pass can run one.
  ⛔ **Do not restate here what
  the freeze does NOT cover** — the trust boundary and the deliberately accepted coverage gaps have a
  canonical home in `fixtures_v4/FREEZE_LIMITS.md` and `fixtures_v4/README.md`; cite those.
  ▶ **The remaining step is the BUILD, and it is HARD-GATED:** `fixtures_v4/run_conformance.py` runs the
  freeze verifier first and refuses conformance at all while verification fails, and the building session
  may neither author nor weaken the fixtures. ⚠ **It also couples to the per-stage loop:** nine of the 36
  governed paths are live dimension-rewrite paths — see the standing block at the head of **§4**, which is
  where a stage-converter meets it.

**(c) Provenance — the missing control is `(g2)` for Python, not a tool.**
- ⭐ **The structural finding, which decides the whole survey: no attestation scheme binds an artifact to
  its producer *by attesting alone*.** ⚠ Scope this precisely — schemes in this family do define and
  carry producer information (SLSA's provenance is verifiable information tracking an artifact to where
  and how it was produced), so the claim is **not** that they say nothing about the producer. It is that
  a digest, a signature, or a link-metadata record is computed **over bytes that already exist**, and
  nothing in the act of hashing or signing can distinguish bytes a process wrote from bytes a hand
  wrote — so the producer information is *asserted*, not established, by the attestation step itself. Every tool that closes the gap closes it by **re-executing the producer and comparing
  bytes** — which is exactly what §4-(g2) already does for the `.wl`. ⇒ **the missing control is (g2) for
  Python**, and unlike Mathematica it needs no licence seat, so it can be **mechanised** instead of left
  to "the orchestrator regenerates the sidecar itself before the (i) commit" (§9).
- **The citable negatives.** *in-toto*: the spec makes an `expected_command` mismatch a **warning** only,
  because the field *"can easily be forged (e.g. by changing the PATH environment variable in a host) and
  thus it should not be trusted for security checks"* — it hashes declared paths before and after the
  wrapped command, so `in-toto-run --products sidecar.txt -- python3 stage.py` yields a link that
  verifies even in a tree where the sidecar was hand-written and the stage never touches it. *DVC*:
  `dvc commit` *"stores the current contents of files … in the cache"* explicitly **without re-executing**,
  documented for *"after executing stage commands by hand"*. *Snakemake*: `--touch` exists *"to pretend
  that the rules were executed, in order to fool future invocations of snakemake."* *SLSA*: L3 makes the
  signature unforgeable **by the build steps** — a genuine by-construction property of the *signature*,
  not of the artifact — while a workflow whose only step copies a committed file still gets perfectly
  valid provenance, and producer-submitted bad code is explicitly out of scope. *sigstore/cosign* signs
  whatever bytes it is handed, and `commit.gpgsign=true` is already configured here.
- **Bazel and Nix DO close it for Python — and their closing act is still regenerate-and-diff.** Bazel's
  sandbox plus `diff_test`/`write_source_files` *"ensures the source tree file … to be written to is up
  to date"*; Nix's equivalent is `nix build && cmp`. The sandbox's contribution is only that the *freshly
  built* side is guaranteed genuine. ⛔ Both are disproportionate for a directory of standalone scripts
  run ad hoc, and **Mathematica's licence blocks hermeticity either way** (a floating licence server needs
  exactly the network the sandbox is designed to remove). Observed-execution provenance (ReproZip /
  Sciunit / noWorkflow) is the only family that binds without re-execution, but its trace is itself a
  file, so it relocates the key-custody question instead of escaping it.
- ⭐ **The cheap control, MEASURED — not a hypothetical.** Copy `scripts/` to a temp dir outside the repo;
  **delete the target `.dimensions.txt` there**; clear `__pycache__`; run the stage; `cmp` the emission
  against the committed bytes; exit 1 on any difference or on a missing emission, with a
  `SIDECAR_REGEN_OK|…` marker in the `MODULE_PIN_OK` style. Run that way — different directory, different
  time, committed sidecar absent — **6 of 6 then-converted stages regenerated BYTE-IDENTICAL, ≈123 s
  total** (004 0.5 s · 011 2.0 s · 012 15.8 s · 013 **82.3 s** · 016 11.3 s · 018 10.6 s). ⚠ Those six
  are the whole measured set: stage023 was not yet converted when this ran and is **not** covered by it
  (§3 now counts seven). The emission is
  deterministic by construction as well as by measurement: no set iteration, no float formatting, no
  timestamps, no absolute paths. ⚠ stage013 alone is ≈82 s and 23 stages remain (§3), so this belongs as
  a **per-stage** gate inside the §4 loop, never as an always-on whole-corpus gate.
- ⭐ **Why deleting the committed copy in the run tree is load-bearing:** it defeats **self-echo** (a
  script that reads its own committed sidecar and reprints it), because the committed bytes are never fed
  to the run. That is precisely why the naive in-place `regenerate && git diff --exit-code` form is
  **weaker** — in place, the committed file is present during the run.
- ⛔ **This closes PRODUCTION PROVENANCE ONLY** — "the artifact is what running that source emits".
  Correctness is a separate property and stays with the cross-engine comparator and the model read (§11).
  Do not let anyone upgrade a green regeneration into a correctness claim; that is the same category
  error §11 already warns about for a bare `python3 …` run.
- ⚠ **What it still does not close** (unchanged from §9): environment interposition — a `sitecustomize.py`
  on `PYTHONPATH` — is inherent to any in-process check, cheaply mitigated by running the regeneration
  inside a pinned `python:3.10-slim` container (the docker daemon is reachable here); and **deleting the
  check's own invocation** is fixed only by an executor outside the repo, the cheapest being a GitHub
  Actions workflow (the repo is public, so minutes are free) whose workflow file is itself at the commit
  under review. ⛔ Python-only — Mathematica cannot run on a hosted runner — but the Python side is
  exactly the side that lacks (g2).

⭐ **The common shape: all three point at a project-owned deliverable, not an adoption** — a 222-line
module kept, an ablation driver the survey sized at ~100 lines, a ~60-line regeneration check. ⚠ Those
are the survey's sizings of the work each would displace, not budgets, and the driver's has not held:
its accepted requirements and authored contract are substantially larger (see (b)). ⛔ **What this does NOT
license:** none of it changes the charter (§1), the per-stage loop (§4), or any existing verdict.

⚠ **The three are NO LONGER in the same state — do not read this section as three candidates.**
**(a) is CLOSED as NO-GO:** `ledger_dimensions.py` stands, no library is adopted, and the migration cost
is recorded only against a future reopening. **(b) is a USER DECISION, 2026-07-29, and it LEADS
the BUILD QUEUE — which is not the conversion order (§8), where 027 is still next, and which now holds
a second item behind it, the shared `DIM|` emitter (below)** — the ablation
driver is to be built, the hand-rolled harnesses switched over to it, and it
reused per stage thereafter, subject to the acceptance-tooling constraints in §12b(b). ⚠ **(b) is no
longer only a decision:** as of 2026-07-29 its requirements are at v4, a contract has been authored
against them, and **its conformance fixtures are frozen and committed** — only the **build** remains, and
it is gated on that freeze verifying (see (b)'s standing-position bullet). **(c) alone
remains a candidate deliverable** — the Python sidecar-regeneration control is recorded so the survey
need not be re-run, and is not a commitment to build.

---

⭐⭐ **SECOND IN THE BUILD QUEUE — the shared Mathematica `DIM|` emitter. USER DECISION, 2026-07-29.**
⚠ **Not a survey outcome** — it is recorded here because this is where the build queue lives, and the
queue now has **two** items: the **ablation driver first (b), the emitter second, both before the 027
conversion**. ⚠ The build queue is not the conversion order (§8): **027 remains the next *conversion***
regardless of where either build lands. The standing "correctness is king, regardless of cost" rule
(§1b) applies here with the same force as it does to the driver.
- **Why.** 23 stages remain (§3) and every one of them currently hand-rolls its own `.wl` emission
  block, re-deriving the same measured hazards each time: a print appended at end-of-file is **dead
  code**, because 43 of 43 `.wl` files terminate in `Exit[]` (§8); a derived `axes=` label printed
  **beside raw storage order** (stage016's axis-label/vector decoupling, `5b29f400`, §3); **duplicate and
  corrupted records from an emission holder invoked many times** (measured 6× on 016, 13× on 021, 17× on
  023, 19× on 027, §9); **hardcoded exponent literals in a `Print`**, which are vacuous against the
  `.py`'s constant (§5); and **fractional values needing exact serialization** (023's `g_U`/`g_W`, plus
  027 and 021, §8) — a float representation must not be reachable.
- **Where it lives.** `research/pde_ledger_v2/notes/wl_emitter/REQUIREMENTS.md:1`, heading
  **"Mathematica `DIM|` emitter — requirements and acceptance intent (v2)"** — at **v2**, design-reviewed
  **three times and CLEAN on the third**, which returned no blocking findings and recorded that no fourth
  specification round was warranted. ⚠ The review logs are in gitignored `_scratch/`
  (`codex_wl_emitter_review{,2,3}.log`) and will not survive; the requirements file is the durable record.
- **The same three-session shape as the driver**, and for the same reason — acceptance tooling cannot
  grade itself (`docs/development_pipeline.md` checklist 1b): **contract → fixtures frozen by a session
  that will not build → build by a third session that may neither author nor weaken them.**
- **Its retrofit oracle (A7).** Applying the emitter to **stage023's `.wl`** must reproduce that stage's
  committed `.out` `DIM|` block — **29 records, `axes=L,M,T`** — **byte-identical** under the existing
  kernel-ID normalization. ⛔ A disagreement is a **finding to report**, never a difference to adjust
  either side into. ⚠ A7 is not decidable until the contract names that extraction-and-normalization
  operation exactly.
- ⚠ **The open decision it carries — and the first thing its contract must settle: a loaded shared file,
  or a per-stage block?** (`REQUIREMENTS.md` contract question 1.) If it is **loaded**, it becomes an
  **unpinned trust root for every stage's reference half** — precisely the position
  `ledger_dimensions.py` occupied before its pin existed (§9), and one edit then compromises every stage
  that loads it. If instead it is a **copied block**, the exposure moves to **copy drift** between
  stages. The contract settles it; what is recorded here is that the question is **live**, and why it
  matters.
- ⭐ **The measured fact R1 rests on.** Of the seven converted stages, **004/013/018** render a **literal**
  `axes=` label while **011/012/016/023** **compute** one; and under the property *one structure feeding
  both the label and every exponent vector*, the compliant set is **at least 011/012/016/023** — 011/012
  build the label from `Keys[dimensionAxes]` and every emitted vector from that same association, while
  016/023 re-project through an explicit `dimensionAxisSlots` table. ⚠ The narrower reading — *"only 016
  and 023"* — is **false** under that property as written. Whether the narrower serialization-time
  re-projection is the property actually wanted is a contract question, not a settled one; and bringing
  004/013/018 across is a **separate per-stage decision** with its own review (§3b records their literal
  labels as bounded debt).
