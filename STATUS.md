# STATUS — where the Path-A program is (single front door)

**A thin pointer, not a copy.** Per-stage detail belongs in the per-Part split docs and per-stage
notes. History is in git. If this file starts growing narrative, cut it.

---

## ▶ YOU ARE HERE (2026-07-29)

**Current front: the DIMENSION REWRITE — 7 of 30 scripts converted** (stage004, 011, 012, 013,
016, 018, **023**) — and ⭐ **all seven are WAIVER-FREE**: `ARTIFACT_NAME_WAIVERS` is empty, so every
converted stage compares every name it emits, with no exemptions.
⚠ That is a **coverage** statement, not a strength one.

⭐⭐ **THE CROSS-CHECK IS EARNING ITS KEEP — measured on stage016.** Relabelling that stage's basis
leaves its **own** 82 assertions completely blind (exit 0, 82 PASS, printing `measure: 'M^3'`), while
the comparator catches **18 of 21**. The comparator is the *sole* instrument between a converted stage
and a relabelled basis — which is why the empty waiver registry matters.

⚠ **`NEEDS_ADJUDICATION` in the canonical table — 3 groups, and all three are correct.** `K_eta`
(three levels: 013 line `M L⁻¹T⁻²` / 016 volume `M L⁻³T⁻²` / 023 reduced scalar `M T⁻²`), `T_Omega`
(016 volume `M L⁻³T⁻²` vs 023 reduced scalar `M T⁻²`) and `mu_eta` (`M L⁻¹` vs `M L⁻³`). These are
**REDUCTION LEVELS, not drift**, surfaced rather than hidden because the variants
were not renamed apart (§7). ⭐ `a` now groups AGREE across 016/018/023 — three stages' `a`
rows now meet in one candidate group and carry the same throat-radius `L`. ⚠ That is **agreement
between implementations, not independence of route**: per `research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md`
§1 (**THE CHARTER**, *"a green comparator shows that two implementations agree, not that they were
reached by two independent routes"*), and git establishes no derivation order for any engine pair.
⚠ `T_w` does **not** group (016 `T_w` → `TW` vs 013 `Tw` → `Tw`), so
that line-vs-volume debt stays invisible — a measured consequence of naming debt, not a theory.

⭐ **TRACK TWO COUNTERS, NOT ONE — they are different finish lines.**
- **converted** — on the shared module, cross-engine gate green: **7 of 30**.
- **physics-verification evidence** — ⚠ **quantity-level; there is NO defensible stage count yet.**
  Recorded: stage012 **14 CORRECT / 0 WRONG**, stage013 **9 / 0**, all six emitted stage018
  records, the three formerly-waived 011/012 records, and ⭐ **stage016's tracked 21 / 0 verdict**
  (12 declared literals in both engines, **0** computed from any physical input). Those pre-§4-c1 results were never
  normalised into per-stage tracked verdict tables, so **do not call whole stages verified and do
  not infer a "25 remaining" complement.** §4-c1 exists so this becomes countable going forward.

⛔⛔ **DUAL-ENGINE AGREEMENT IS VACUOUS WHERE BOTH SIDES ARE HAND-DECLARED LITERALS — recorded
2026-07-30, and it is the most important thing on this page.** **stage023 is 29/29 literals on both
sides** (22 declarations + 7 targets typed into *both* engines; the 7 "computed" records are live
`dimOf` walks, but over exactly those literals and no other dimensional input). Beside the stage016
measurement above — **12 declared literals in both engines, zero computed from a physical input** —
that means the comparator on these stages catches a **transcription split between two typed copies of
the same numbers, and nothing else**. ⇒ **Nothing remaining independently RE-DERIVES the physics
outside one fresh agent.** ⛔ Do not soften this and ⛔ do not read a fix into it: it is stated here
because it is the reason the **derived-vs-declared census is the next front**.

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

✅ **The shared module is no longer self-attesting** (2026-07-27) — conversion is unblocked. ⚠ What that
digest earns is a **staleness** signal only (*the module changed and the stages were not re-run*); see the
digest block below.
✅ **stage023 IS CONVERTED — both halves.** Its `.py` is on the shared module, the comparator is green
(`py=29|wl=29|shared=29|py_only=0|wl_only=0|mismatches=0`, no waivers), the orchestrator regenerated
both engines' artifacts itself, and the sealed prediction is adjudicated (**5 fully confirmed · 1
falsified (P2) · 1 split** — P3's mechanism confirmed, its exclusivity falsified by U13).
Evidence: stage note §1.6 (§4-a enumeration), §1.7 (§4-c1 physics verdict), §5.1 (steps g/g2/h/h2).
⛔ The `.wl` never joins the module: the charter is **SymPy-only** and the Mathematica side is authored as
its own route (`DIMENSION_REWRITE.md` §1), which is why the comparison is a permanent standing cross-check.
⚠ Read that as the **authoring rule** it is — §1's own honest statement governs: a green comparator shows
**two implementations agree, not that they were reached independently**, and git establishes no order.
**▶ NEXT CONVERSION, per the recorded conversion order (`DIMENSION_REWRITE.md` §8):** (1) the
stage027-shape decision, (2) 027, (3) 021 (heaviest). Detail and the measured validator/harness hazards
are in `DIMENSION_REWRITE.md` §8/§9.
✅ **The ablation-fixture FREEZE AUTHORITY IS RETIRED** (user decision, 2026-07-29/30) — with it, the
coupling that made nine live dimension-rewrite paths untouchable. Convert freely; nothing here is frozen,
byte-perfect or under a custody rule. See `DIMENSION_REWRITE.md` §4.
**▶ NEXT BUILD — a different queue, not a competing answer. It has TWO items, in this order:**
**(1) the ablation driver**, **RE-SCOPED SMALL** (user decision, 2026-07-29/30, `DIMENSION_REWRITE.md`
§12b(b)): mutate a declaration, confirm the declared assert fires, record it — reviewed by one fresh
agent. ⛔ No contract, no frozen fixtures, no three-session shape. Requirements (trimmed) at
`research/pde_ledger_v2/notes/ablation_driver/REQUIREMENTS.md`; `CONTRACT.md` and the `fixtures_v4/` suite
are superseded and stay only as history — ⭐ §C-9's legacy mapping table is worth reading as a **reference**
for the retrofit, but ⛔ it is **not** authoritative and A7 does not require agreement with it.
**(2) the shared Mathematica `DIM|` emitter** (`DIMENSION_REWRITE.md` §12b, closing block) — **just write
it**, against `research/pde_ledger_v2/notes/wl_emitter/REQUIREMENTS.md` (v2, reviewed CLEAN). ⛔ No
contract, no frozen fixtures, no three-session shape. Which stage converts next and which
tooling gets built next are separate sequences; in time **both** builds land **before** 027 begins.

**▶ WHAT STAGE023 LEFT OPEN.** Stage023's tracked physics leg records 34 quantity routes and the
two scoped tallies `24/0/10` on corpus identifications and `27/0/7` inside the stage-local closure
(`research/pde_ledger_v2/notes/stages/ledger_stage023_nullspace_underdetermination.md:379-440`,
**§1.7(1), per-quantity verdict, tally, and the 34-row route table**). Its seven unresolved derivation questions now have named
work routes in `research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md`: `:938`
(**WORK-023-MOMENT-CONVENTION**), `:998` (**WORK-023-STAGE009-MOMENT0**), `:1027`
(**WORK-023-D0-SEAM**), `:1061` (**WORK-023-STIFFNESS-REDUCTION**), `:1100`
(**WORK-023-L1-L2-PROFILE-IDENTITY**), `:1124` (**WORK-023-CS-EVALUATION**), and `:1166`
(**WORK-023-SOURCED-PROVENANCE**), all in **§12**. W3 is outside that work list because its confirmed
arithmetic correction is already folded at
`research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md:648-663`
(**§8, stage023 `gU/gW` correction**); `q_free` is outside it because §1.7 classifies that record as
an unread control rather than a competing identification
(`research/pde_ledger_v2/notes/stages/ledger_stage023_nullspace_underdetermination.md:428,442-452`,
**§1.7(1), `q_free` verdict and tally explanation**). ⚠ These §12 items are derivation work, not user
gates, and they do **not** hold back the counter above: stage023 is converted on **both** halves and
is counted in the **7 of 30**.

✅ **stage016 COMPLETE** — both engines, comparator `py=21|wl=21|shared=21|mismatches=0`.
✅ The four group-A stages have `.wl` reachability verdicts — that recorded GAP is **closed**;
⚠ but those counts are **PROVISIONAL, not completeness proofs**: 016's survey got its 21 emitted
quantities right, yet missed two **non-emitted** source/control-flow cases — so the method does not
close the broader inventory. Each stage's tracked §4-a enumeration + adversarial review is what does.
027 is the awkward one (**MIXED**: its **computed** vector never reaches top level — it dies in `runAll`'s
`Module`. ⚠ **Corrected 2026-07-30:** its **16 declared per-symbol `baseDims` vectors DO reach top level**,
so 027 is **not** confined to a 1-row stage; only reaching the computed vector needs new call sites).
⛔ **Measured corpus-wide: all 43 `.wl` files end in `Exit[]`,** so a print appended at end-of-file is
dead code. Emit before the terminal `Exit[]`. See `DIMENSION_REWRITE.md` §8/§9.
⛔ **Before converting 035/036/037, fix the canonical-table generator** — it *raises* on cross-engine
axis-order disagreement. ⚠ **What was actually measured, stated precisely:** the committed stage037 `.wl`
declares **no axis order at all** — it is a formal `massScale`/`lengthScale`/`timeScale` rescaling with no
ordered vector and no `axes=` header — so the crash was **not** measured on the committed pair. The
`M,L,T` order exists only in the **gitignored** `_scratch/spike037/` prototype (`dimensionAxes` keyed
`M`,`L`,`T`), against the committed `.py`'s `(L,T,M)` (`ledger_stage037_route_b_boost_structural_relation_sympy_audit.py:604`).
⇒ The disagreement is real and will appear the moment 037 emits, but it is a **spike-vs-`.py`**
measurement, not a committed-pair one.

Per stage: emit `DIM|` records into the `.wl` → ⭐ **run the PHYSICS leg FIRST** (it is
blocking, and a NAMING decision is a physics decision) → rewrite the `.py` onto the
module → **re-run the `.py`, then** compare axis-labelled → **regenerate the `.out` and byte-compare it**
(⭐ load-bearing: an uncaught `Throw` exits 0 with an empty or truncated transcript and only the byte-compare
catches it — `DIMENSION_REWRITE.md` §9. ⛔ **No longer an orchestrator-only duty** — whoever runs the stage
does it; the old *"agents cannot run Mathematica"* reason was measured false and the
independent-party replacement was never adjudicated, so it is cut as a custody step, `DIMENSION_REWRITE.md`
§4-(g2)) → **two mutually independent fresh review legs** (each fidelity + per-tooth ablation; a stage's
math is physics-bearing) → commit.
⛔ **Retired 2026-07-29/30:** commit-the-`.wl`-before-the-`.py` reference custody. ⭐ **What replaces it and
is NOT ceremony: the co-authorship guard — the party that wrote the `.py` must not be the party that adjusts
the `.wl` until the comparator agrees** (`DIMENSION_REWRITE.md` §4-(d)); tuning whichever side disagrees is
the LLM-shortcut-that-resembles-a-pass, not a fix.
⭐ **The prediction goes in the SESSION SCRATCHPAD and is folded into the stage note after the reviews land**
(restored 2026-07-30) — because *before launching any agent, ask what is reachable from the working tree that
states your expected answer; absence beats instruction* [`never-supply-the-expected-reason`]. ⛔ That is the
whole reason: no custody, no sealing, no ordering ritual. The physics leg is still told to **derive from the
model**, never to *check a claim*.

⭐ **The hard tail is bounded.** A spike **prototyped** an independent `.wl` route for **stage037**
(`ROUTE_EXISTS`, 21/21 quantities, a real comparator failure on a seeded error), so the old
"genuinely impossible" verdict is **false for 037**. **035 and 036 have identified routes that are
NOT yet prototyped** — source inspection only. ~0.5–1 engineer-day per stage is an **estimate**, not
a measurement. See `DIMENSION_REWRITE.md` §3b.

⭐ **Python sidecars are source-hash bound, and that is exactly a STALENESS check — which is all it is
asked to be.** It catches the real, common error: *the `.py` changed and the sidecar was not
regenerated.* ⚠ It does **not** prove the sidecar came from a run — a hand-written one carrying the
right digest reaches `PASS` (demonstrated 2026-07-27). ⛔ **Retired 2026-07-29/30: the orchestrator-
regenerates-the-sidecar control, and calling this an open hole.** Forging a sidecar is a *motivated-
adversary* move, and this project hardens against **drift and honest error**, not that
(`docs/development_pipeline.md`, *THE POSTURE*). What catches a wrong sidecar is the comparator plus the
blocking physics leg. See `DIMENSION_REWRITE.md` §9.

⭐ **THE SHARED-MODULE DIGEST — A STALENESS PING (downgraded from an "authority" by user decision,
2026-07-29/30).** `scripts/check_ledger_dimensions_pin.py` compares `scripts/ledger_dimensions.py`'s
current source bytes against `scripts/ledger_dimensions.accepted.sha256`, and fails on any difference.
⛔⛔ **State only that.** It **cannot** establish *"you edited the module and did not re-run the stages"*:
it inspects no producer and no run — it compares two hashes, and `--accept` rewrites the recorded hash from
the current module bytes without checking that anything was re-run. A red digest means **the module differs
from the last accepted baseline**, which is usually a module edit whose downstream has not been refreshed —
useful, cheap, and *not* evidence about the stages either way. ⭐ **The check that actually detects a stage
not re-run after a module edit is the SIDECAR binding** (the module digest is stamped into each stage's
sidecar header and the comparator recomputes it). Stubbing `dim_residual` does trip the pin in the standalone
control, the comparator and the generator (class `MODULE_PIN_MISMATCH`, distinct from sidecar staleness).
⛔ **RETIRED FRAMING — do not restore it:** "an authority **no producer writes**", the digest as a
**validator** whose green says anything about correctness, and **re-acceptance as a review event with a
recorded reason and a second witness**. A red digest after a legitimate module edit means **refresh and
re-run the producers** (`--accept`, then stage → comparator → generator); that is a reset, not a trust
decision. Procedure: `DIMENSION_REWRITE.md` §4.
⚠ **Its bound:** it covers no stage source, no `.wl`, no `.out` and no sidecar content, and it never executes
a stage. ⇒ **Read a green digest as "not stale", nothing more.** ⛔ **Clear `scripts/__pycache__/` after any
ablation edit/restore loop** — equal-size edits (sign flips, `sum`→`min`) let CPython reuse timestamp-valid
stale bytecode by accident; that one is live and practical. *(The bytecode / `sitecustomize` / trust-root
analysis that used to sit here is cut 2026-07-30: those are motivated-adversary routes around a staleness
ping, which the governing test does not buy.)* Detail: `DIMENSION_REWRITE.md` §9.

⭐ **A bare stage run is a PRODUCER, not a validator** (user decision, 2026-07-27) — its exit code and
`PASS` tally are **not** validation evidence. The validators are the comparator and the generator (the
module digest is a staleness check, not a validator). §11 had already measured why (a relabelled basis leaves stage016's own 82 assertions blind
at exit 0); this states it as a rule, and it is what lets the pin sit outside the module.
⚠ `run_all_audits.sh` tallies `Fail: N` but **exits 0** — it gates on the pin, not on audit failures,
and it never invokes the comparator or generator.

⚠ **Before resuming, read `DIMENSION_REWRITE.md` §1b (the D1–D5 decisions) and §3b (what those
decisions REOPENED).** Several recorded conclusions — three waivers, four "impossible" stages, a
coverage estimate — were correct only under constraints since lifted, and will read as settled.

**Why this is the front:** ⭐ **one shared import, so every script's dimensions are written in ONE
REPRESENTATION and come from one place.** ⛔ **Not "consistent by construction" — that overstates.** One
shared module buys **representation unity** (one basis type, one exponent type, one recovery path); it does
**not** buy correct dimensions, because **two stages can declare the same wrong exponents through one module**
just as easily as through thirteen idioms. Correctness comes from the blocking physics leg.
Thirteen dimension idioms across 43 scripts *is* drift; the decision
(`b5527062`, `aae5d389`) was to **fix the corpus, not weaken the check**. Part VII's whole-system
dimensional firewall is **claimed** to consume the module directly rather than the manifests — ⚠ **an
assertion, not an established fact:** stage046 is unbuilt, and `notes/part7_integration_atomic_split.md`
(the 046 row) names the firewall without naming its input source. Do not cite it as settled. The manifests' semantic
core continues in parallel, trimmed (§ "PAUSED" below).
⛔ **The old justification — "the 44-stage manifest fanout is blocked, dimension recovery covers only
~16 of 43 scripts" — is RETIRED as false.** Verified 2026-07-29/30: the composite checker recovers
dimensions from **10 of 43** scripts — exactly **7** carry a `class Dim` the recovery walks (005, 006,
007, 008, 009, 030, 031) plus **3** registered bare-tuple digests (032, 038, 042) — and **all seven
converted stages carry none**, because the shared module exports `Dimension`, not `Dim`. Conversion
therefore does not raise that count, and the fanout was never what the rewrite unblocks.
⚠ **Precisely:** conversion *lowers* recovery only for the **seven `class Dim` stages**; for every other
script it merely fails to raise it, because there was nothing recoverable there to lose.

## ⏸ PAUSED (resume after the rewrite; user confirms sequencing)

- **stage 044-v2** — redo stage044 with a DYNAMICAL-Σ sleeve (un-freeze `S_hold`, commit the
  `κ_bend / κ_anchor / collar-tension` bending knobs). User-decided, 044-LOCAL.
  Anchor: `research/pde_ledger_v2/notes/stage044_v2_unfreeze_prep.md`.
- **stage 045** (VII-2b) — the non-variational drain/return block + BCs + force partition, where the
  drain-placement crux and the USER mini-gate land. Drain = the dynamical `Γ_B`; frozen-wall ruled out.
  ⭐ **That mini-gate STANDS** — it is a modelling DECISION for the user, which the reduced process still stops
  for; it is not a per-chunk gate.
  Anchor: `research/pde_ledger_v2/notes/stage045_nonvariational_block_prep.md`.
- **Manifest / integration-test system** — built + committed (`e849e303`), 4 of 44 manifests extracted.
  ⭐ **CONTINUES, ON ITS SEMANTIC CORE** (user decision, 2026-07-29/30; corrected 2026-07-30): quantity
  identity, **citation integrity**, **export/lifecycle enforcement**, dimensional relations, **the
  lifecycle census**, the dependency graph and its cycles, mutation, **genesis**, consumption completeness.
  ⛔ **The trim that dropped citation-integrity / lifecycle / genesis is WITHDRAWN** — each catches a way the
  *physics* could be wrong (a changed consumed equation · consumption of retired physics · a wrong
  irreducible-count range · calibrated-or-target-matched genesis, i.e. fit-vs-derive), not bookkeeping.
  ⚠ It is **not** what the dimension rewrite unblocks (see the front's justification above) — the two are
  independent. Docs:
  `research/pde_ledger_v2/manifests/MANIFEST_README.md` +
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

- ⭐⭐ **PHYSICS, NOT CEREMONY** (user decision, 2026-07-29/30). Two-person toy-physics self-consistency
  project; checks exist to catch a **wrong derivation**. ⛔ **Immutability is not a discipline here** —
  files are freely editable, and nothing is frozen, byte-perfect, or under a custody rule. **The test for
  any check or process rule: does it catch a way the PHYSICS could be wrong? → keep. Only a way the
  TOOLING could be wrong, or a motivated adversary? → cut.** Roles collapse to **one builder; one fresh
  reviewer for prose and process, TWO mutually independent fresh review legs for physics-bearing artifacts**
  (amended 2026-07-30), with the physics leg **blocking**. Owned by `docs/development_pipeline.md`.
- **Findings are the product; green is not the goal.** A result that breaks the concept is welcome and
  first-class. A clean "it all works" is suspicious.
- **The ledger shows the SURVIVING solution only.** Discarded approaches →
  `research/pde_ledger_v2/notes/ledger_exclusions_failures_paper_backlog.md`.
- **Never adjust the process because the corpus is inconvenient.**
- **Only held-out DIMENSIONLESS ratios test the model.** `G`, `c`, `ℏ`, `ℓ_P` are calibration.
- **AI prose never establishes a math fact.** The builder codes and runs dual-engine; verification happens
  on **fresh agents** — **two independent ones where the artifact is physics-bearing**, one for prose and
  process. ⚠ Stop for the user at a decision, a blocking finding or a no-go — **not** at
  every chunk boundary. See `docs/development_pipeline.md`.
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
