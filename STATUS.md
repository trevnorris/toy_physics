# STATUS — where the Path-A program is (single front door)

**A thin pointer, not a copy.** Per-stage detail belongs in the per-Part split docs and per-stage
notes. History is in git. If this file starts growing narrative, cut it.

---

## ▶ YOU ARE HERE (2026-07-27)

**Current front: the DIMENSION REWRITE — 6 of 30 scripts converted** (stage004, 011, 012, 013,
**016**, 018) — and ⭐ **all six are WAIVER-FREE**: `ARTIFACT_NAME_WAIVERS` is empty, so every
converted stage compares every name it emits, with no exemptions.
⚠ That is a **coverage** statement, not a strength one.

⭐⭐ **THE CROSS-CHECK IS EARNING ITS KEEP — measured on stage016.** Relabelling that stage's basis
leaves its **own** 82 assertions completely blind (exit 0, 82 PASS, printing `measure: 'M^3'`), while
the comparator catches **18 of 21**. The comparator is the *sole* instrument between a converted stage
and a relabelled basis — which is why the empty waiver registry matters.

⚠ **First-ever `NEEDS_ADJUDICATION` in the canonical table — 2 groups, and both are correct.**
`K_eta` (013 line-density `M L⁻¹T⁻²` vs 016 volume-density `M L⁻³T⁻²`) and `mu_eta` (`M L⁻¹` vs
`M L⁻³`). These are **REDUCTION LEVELS, not drift**, surfaced rather than hidden because the variants
were not renamed apart (§7). ⭐ `a` now groups AGREE across 016/018 — the same throat radius, reached
independently by two blind legs. ⚠ `T_w` does **not** group (016 `T_w` → `TW` vs 013 `Tw` → `Tw`), so
that line-vs-volume debt stays invisible — a measured consequence of naming debt, not a theory.

⭐ **TRACK TWO COUNTERS, NOT ONE — they are different finish lines.**
- **converted** — on the shared module, cross-engine gate green: **6 of 30**.
- **physics-verification evidence** — ⚠ **quantity-level; there is NO defensible stage count yet.**
  Recorded: stage012 **14 CORRECT / 0 WRONG**, stage013 **9 / 0**, all six emitted stage018
  records, the three formerly-waived 011/012 records, and ⭐ **stage016's tracked 21 / 0 verdict**
  (12 declared literals in both engines, **0** computed from any physical input). Those pre-§4-c1 results were never
  normalised into per-stage tracked verdict tables, so **do not call whole stages verified and do
  not infer a "25 remaining" complement.** §4-c1 exists so this becomes countable going forward.

Cross-engine agreement is necessary and **not sufficient** — proven twice: it is blind to a
same-dimension different-quantity merge, and to a shared wrong declaration. The physics leg is
what establishes correctness, and **it does not depend on conversion** — it can run against any
stage's existing declarations. If the goal that matters is Part VII's firewall being trustworthy
rather than the manifests unblocking, that is the leg to sequence around. See `DIMENSION_REWRITE.md`
§4-c1.

Every SymPy audit script's dimension handling moves onto one shared module
`research/pde_ledger_v2/scripts/ledger_dimensions.py`, **one stage at a time**, each verified by an
independent `.py`-vs-`.wl` cross-check.

> ⭐ **READ `research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md`** — the single canonical doc for
> this workstream. And read `docs/model_map.md` **before touching any script**.

✅ **The module pin is DONE** (2026-07-27) — the self-attestation that blocked conversion is closed, and
conversion is unblocked. ⚠ It earns **one narrow fact**, not a clean bill of health; see the pin block
below for exactly what it does and does not cover.
**▶ NEXT: (1) stage023, (2) the stage027-shape decision, (3) 027, (4) 021.**
✅ **stage016 COMPLETE** — both engines, comparator `py=21|wl=21|shared=21|mismatches=0`.
✅ The four group-A stages have `.wl` reachability verdicts — that recorded GAP is **closed**;
⚠ but those counts are **PROVISIONAL, not completeness proofs**: 016's survey got its 21 emitted
quantities right, yet missed two **non-emitted** source/control-flow cases — so the method does not
close the broader inventory. Each stage's tracked §4-a enumeration + adversarial review is what does.
027 is the awkward one (**MIXED**: its computed vector never reaches top level, and its `.wl` route
cannot produce per-symbol vectors, so it stays a 1-row stage unless new call sites are added).
⛔ **Measured corpus-wide: all 43 `.wl` files end in `Exit[]`,** so a print appended at end-of-file is
dead code. Emit before the terminal `Exit[]`. See `DIMENSION_REWRITE.md` §8/§9.
⛔ **Before converting 035/036/037, fix the canonical-table generator** — it *raises* on cross-engine
axis-order disagreement, which is exactly the 037 pattern (`wl M,L,T` vs `py (L,T,M)`). Measured; it
will crash the toolchain.

Per stage: emit `DIM|` records into the `.wl` → ⭐ **run the PHYSICS leg BEFORE the commit** (it is
blocking, and a NAMING decision is a physics decision) → **commit the `.wl` before touching the `.py`**
(freezing the reference first is reference custody — NOT proof of independent authorship) → rewrite the `.py` onto the
module → **re-run the `.py`, then** compare axis-labelled → **orchestrator regenerates the `.out`**
(verify agents cannot run Mathematica, so only you can confirm the reference is genuine) → fidelity +
adversarial fresh agents → commit. Keep the prediction note **outside the repo** until the reviews land.

⭐ **The hard tail is bounded.** A spike **prototyped** an independent `.wl` route for **stage037**
(`ROUTE_EXISTS`, 21/21 quantities, a real comparator failure on a seeded error), so the old
"genuinely impossible" verdict is **false for 037**. **035 and 036 have identified routes that are
NOT yet prototyped** — source inspection only. ~0.5–1 engineer-day per stage is an **estimate**, not
a measurement. See `DIMENSION_REWRITE.md` §3b.

⚠ **Python sidecars are source-hash bound — but that closes STALENESS, not TAMPERING.** The digest
proves *"this sidecar names the `.py` on disk"*, never *"it was produced by running it"*, and the
comparator never executes the stage. **Demonstrated 2026-07-27:** a hand-written sidecar carrying a
mutated `.py`'s digest reaches `PASS`/exit 0 while that `.py` declares wrong values. ⛔ **That hole is
still OPEN.**
⭐ **Interim control: the orchestrator regenerates the sidecar itself before committing.**
See `DIMENSION_REWRITE.md` §9.

✅ **THE SHARED MODULE IS NO LONGER SELF-ATTESTING (2026-07-27).** `scripts/ledger_dimensions.py` is now
pinned by `scripts/check_ledger_dimensions_pin.py` against `scripts/ledger_dimensions.accepted.sha256`
— an authority **no producer writes**. ⚠ **State the independence with its mode qualifier:** in
*validation* mode the expected digest never comes from the module; the explicit `--accept` operation
*intentionally does* hash the module and rewrite the authority
(`check_ledger_dimensions_pin.py:101-115`). ⭐ **The module itself is UNCHANGED**,
so the pin lives entirely *outside* `ledger_dimensions.py` and cannot be removed by editing **that**
file. Stubbing `dim_residual` now fails the control (1), the comparator (1) and the generator (2) with a
`MODULE_PIN_MISMATCH` class distinct from sidecar staleness; ⭐ **and re-running every producer does not
launder it** — the property the old arrangement failed. Re-baselining is the explicit `--accept`
operation, documented in `DIMENSION_REWRITE.md` §4.

⚠ **Scope it precisely: this earns ONE narrow fact** — the **SHA-256 of** the current source bytes of
`scripts/ledger_dimensions.py` equals the deliberately accepted SHA-256. ⛔ No accepted *byte sequence*
is retained anywhere; the authority holds a digest only, so "the bytes are equal" is not what is
checked. It is **not** a guarantee that the dimensional gates ran honestly.

⚠ **What the pin does NOT cover — do not overstate it (the second time this workstream has had to say
so).** It pins the **source bytes of one file**. Not **bytecode** (`__pycache__` is gitignored;
⭐ measured: CPython rejects timestamp bytecode if the stored mtime **or** the stored source size
differs (`_bootstrap_external.py:637-643`), and it truncates mtime to a whole second (`:973`) — so
**BOTH conditions are required**, and two ordinary *equal-size* writes inside one second satisfy them
with no header construction. Equal-length edits — sign flips, `sum`→`min`, `+`→`-` — are exactly the
dangerous ones in a dimensional audit, which is what makes the size condition routine here). Not its own **execution
environment** (a `sitecustomize.py` spoofing `hashlib` for the module's bytes takes every validator
green — inherent to any in-process check). Not stage sources, the `.wl`, the `.out`, or the forgeable
sidecar. And the authority is a **bare trust root**: `--accept` moves the baseline with no reason field,
signature, or second witness. ⛔ **And not the pin's own code**, which is unpinned source — but state the
two cases separately, they are not equivalent:
- **the shared trust root** — the decision function `check_ledger_dimensions_pin.py:75-98` plus the
  authority file. Editing it compromises **every** consumer at once. ⇒ **that** is where the pin moves
  the single point of failure rather than abolishing it.
- **the four individual call sites** — `compare_dimension_artifacts.py:251`,
  `generate_canonical_dimension_table.py:231` **and** `:461`, `run_all_audits.sh:20`. Deleting one
  compromises **that path only**; the standalone control and the remaining validators still reject (the
  generator holds two, so removing one leaves the other live).

⭐ **A bare stage run is a PRODUCER, not a validator** (user decision, 2026-07-27) — its exit code and
`PASS` tally are **not** validation evidence. The validators are the comparator, the generator, and the
pin control. §11 had already measured why (a relabelled basis leaves stage016's own 82 assertions blind
at exit 0); this states it as a rule, and it is what lets the pin sit outside the module.
⚠ `run_all_audits.sh` tallies `Fail: N` but **exits 0** — it gates on the pin, not on audit failures,
and it never invokes the comparator or generator.

⚠ **Before resuming, read `DIMENSION_REWRITE.md` §1b (the D1–D5 decisions) and §3b (what those
decisions REOPENED).** Several recorded conclusions — three waivers, four "impossible" stages, a
coverage estimate — were correct only under constraints since lifted, and will read as settled.

**Why this is the front:** the 44-stage manifest fanout is blocked — dimension recovery covers only
~16 of 43 scripts. The decision (`b5527062`, `aae5d389`) was to **fix the corpus, not weaken the
check**. One shared module gives the checker one recovery path and feeds Part VII's whole-system
dimensional firewall.

## ⏸ PAUSED (resume after the rewrite; user confirms sequencing)

- **stage 044-v2** — redo stage044 with a DYNAMICAL-Σ sleeve (un-freeze `S_hold`, commit the
  `κ_bend / κ_anchor / collar-tension` bending knobs). User-decided, 044-LOCAL.
  Anchor: `research/pde_ledger_v2/notes/stage044_v2_unfreeze_prep.md`.
- **stage 045** (VII-2b) — the non-variational drain/return block + BCs + force partition, where the
  drain-placement crux and the USER mini-gate land. Drain = the dynamical `Γ_B`; frozen-wall ruled out.
  Anchor: `research/pde_ledger_v2/notes/stage045_nonvariational_block_prep.md`.
- **Manifest / integration-test system** — built + committed (`e849e303`), 4 of 44 manifests extracted.
  This is what the dimension rewrite unblocks. Docs: `research/pde_ledger_v2/manifests/MANIFEST_README.md` +
  `research/pde_ledger_v2/manifests/EXTRACTION_PROTOCOL.md`.

## LEDGER BUILD STATUS

| Part | Sector | Status |
|---|---|---|
| 0 | Conceptual | scaffolding only |
| I | Medium | ✅ built (004–007) |
| II | Gravity | ✅ COMPLETE (001–002, 008–029) — sector closed |
| III | Light | ✅ DONE = stage003 (surviving-solution rule) |
| IV | Charge | ✅ COMPLETE (030–033) |
| V | Magnetism | ✅ COMPLETE (034–039) |
| VI | Knit | ✅ COMPLETE (040–042) |
| VII | Integration | 🔄 **2 of 7** — 043 ✅, 044 ✅; 045–049 remain |

## STANDING RULES

- **Findings are the product; green is not the goal.** A result that breaks the concept is welcome and
  first-class. A clean "it all works" is suspicious.
- **The ledger shows the SURVIVING solution only.** Discarded approaches →
  `research/pde_ledger_v2/notes/ledger_exclusions_failures_paper_backlog.md`.
- **Never adjust the process because the corpus is inconvenient.**
- **Only held-out DIMENSIONLESS ratios test the model.** `G`, `c`, `ℏ`, `ℓ_P` are calibration.
- **AI prose never establishes a math fact.** Codex codes and runs dual-engine; verification happens on
  fresh agents; user gate per gate. See `docs/development_pipeline.md`.
- **One canonical doc per workstream.** Fold new findings in — do not write a new doc re-explaining
  ground an existing one covers. Delete dead docs (git preserves them), but extract unique content first.

## MAP — what you want → where it lives

| You want… | Look here |
|---|---|
| ⭐⭐ **The model** — throughline, per-sector derivation atlas, honest earned/calibrated/R1/departure ledger, glossary | `docs/model_map.md` |
| ⭐ **The current front** — the dimension rewrite | `research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md` |
| ⭐ **How we work** — pipeline, roles, the review gauntlet | `docs/development_pipeline.md` |
| The ledger-build resume detail | `research/pde_ledger_v2/notes/RESUME_ROADMAP.md` |
| What's left across ALL sectors → "simulation-ready" | `docs/development_plan.md` |
| Per-Part build history + user-gate records | `research/pde_ledger_v2/notes/part{1..7}_*_atomic_split.md` |
| Per-stage notes | `research/pde_ledger_v2/notes/stages/` |
| Every knob: dimension, class, provenance, reduction debt | `research/pde_ledger_v2/notes/parameter_register.md` |
| The irreducible-count audit | `research/pde_ledger_v2/notes/midway_knob_audit.md` |
| Every value classified DERIVED / INPUT / gap / benchmark | `software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md` |
| The simulation hand-off spec (equations + BC packet + R1→R4 ladder) | `docs/two_throat_simulation_handoff_spec.md` |
| Numbered decisions (esp. **16**, the `P`-retirement) | `software/stage1_solver/decisions/` |
| Retired approaches / the failures-paper backlog | `research/pde_ledger_v2/notes/ledger_exclusions_failures_paper_backlog.md` |
| The calibrate-predict methodology | `software/stage1_solver/decisions/09_calibrate_predict_methodology.md` |
| ⭐ **The EM-track record** — U1/U2 + Phase B/C, the `𝔅` boundary-operator verdict (144/144 UNRESOLVED). **Sole home of that narrative** | `docs/em_analog_next_phase_handoff.md` |
| The EM physical picture + MacCullagh template + leak findings | `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md` |
| Full current state + resume-after-`/compact` pointer (**sync this file with it**) | `software/stage1_solver/decisions/13_emergent_constants_derivation.md` §0 |
| The defect-regime + held-out-surplus roadmap | `docs/defect_interaction_map.md` |
| The gravity moving-throat PDE gate checklist | `research/pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md` |
| The pre-registration (what was committed to in advance) | `docs/pathA_preregistration.md` |

⛔ **Do NOT read `docs/conceptual_foundation.md`** — vision/history, superseded by `docs/model_map.md`. It
predates the EM reconsideration and the retired `P` field, and re-confuses.

## Reference — the `pathA_22b` verdict equation (the *earlier* verdict-count framing)

```
P0 · χ_Q · g_mhat² · λγ⁵ / g_G  =  54/5
 ✓     ✓     gap1     gap2  cal-on-G     (✓ = derived; gap1 g_mhat absorbs 54/5; gap2 λγ ← EM anchor)
G = (a·c_s²/m_GNLS)·g_G ,  m̂0 = (c_s/(a²·√m_GNLS))·g_mhat ,  c = λγ·c_s
```

Here `χ_Q ≈ 0.712` (pathA_22b Gate-3). The moving-throat ladder's Gate 4 (`pathA_33`) later derived
**`χ_Q = 1`** in the outgoing-DtN Hankel context — a different computation of the same-named factor.
**Reconciling the two is a live Part-VII debt**; the ladder reaches the same `54/5` via the
earned/calibrated split `54/5 = 2·27/5`, not via this product.
