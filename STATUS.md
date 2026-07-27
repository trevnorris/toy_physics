# STATUS — where the Path-A program is (single front door)

**A thin pointer, not a copy.** Per-stage detail belongs in the per-Part split docs and per-stage
notes. History is in git. If this file starts growing narrative, cut it.

---

## ▶ YOU ARE HERE (2026-07-26)

**Current front: the DIMENSION REWRITE — 3 of 30 scripts done** (stage004, stage011, stage012).

Every SymPy audit script's dimension handling moves onto one shared module
`research/pde_ledger_v2/scripts/ledger_dimensions.py`, **one stage at a time**, each verified by an
independent `.py`-vs-`.wl` cross-check.

> ⭐ **READ `research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md`** — the single canonical doc for
> this workstream. And read `docs/model_map.md` **before touching any script**.

**▶ NEXT: stage013**, then 018, 016, 023, 027, 021. Per stage: emit print-only `DIM|` records into the
`.wl` and **commit that before touching the `.py`** (freezing the reference first makes engine
independence structural) → rewrite the `.py` onto the module → compare axis-labelled → verify
(adversarial fresh agent + a Grok pass checking the dimensions against the *model*, not just against
the other engine) → commit.

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
