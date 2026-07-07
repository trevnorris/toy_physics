# Ledger v2 — Build Blueprint (Phase B, step B1)

> **Status: APPROVED by user 2026-07-07** (design decisions A–D + completeness standards settled below).
> Codex design-reviewed (NEEDS_REVISION → folded → confirm CLEAN); revised with the user's design steer.
> This is the foundational directive for building the rebuilt PDE ledger.
> **Reads-first:** `notes/ledger_v2_rebuild_plan.md` (multi-phase plan + why-rebuild evidence), `STATUS.md`
> ▶ RESUME HERE, `software/stage1_solver/decisions/13` §0. Memory: [[project-pde-ledger-fold-in-plan]],
> [[project-brane-existence-defect-structure]].

---

## 0. Purpose, governing goals & guardrail

**What this delivers:** a concrete, reviewed plan for assembling the four force-sector packages (each an earned core +
its calibrated / sim-deferred / characterized-departure scopes) + the consistency knit into a rebuilt,
sector-organized ledger, reusing the old 253-stage ledger's *machinery* and none of its search-sprawl content.

**⭐ Governing goals (user, 2026-07-07) — these drive every design decision:**
1. **Zero context loss.** Every derivation step, calibration input, characterized departure, and no-go is preserved.
2. **One self-contained ledger** a future physicist can fully work the model from — the ledger *is* the content, not
   a set of pointers to `pathA_*` reports.
3. **A source corpus** from which focused papers get extracted per sub-topic (so stages must be atomic and citable).

**Guardrail (unchanged, load-bearing):** building the ledger completes the **SPEC, not the proof.** Gate-6 nonlinear
branch realization stays **SIM-DEFERRED**; every sim-dependent quantity is preserved as a precisely-posed open item
(Part VII); a no-go is still possible. The ledger is a faithful, dual-engine-verified *assembly* of what is already
earned + a bounded amount of genuinely-new integration synthesis (Part VII + the χ_Q reconciliation) — it does
**not** re-open or re-derive the earned physics.

---

## 1. Where it's built + the path swap (Phase B2) + machinery copy

**Build-alongside-then-swap, on a dedicated branch (user-approved 2026-07-07).** The entire rebuild happens on a
dedicated branch **`ledger-v2-rebuild`** (off `master`), so git history shows the full incremental development of the
ledger. Build the new ledger in a **transient staging dir `research/pde_ledger_v2/`** while the old
`research/pde_ledger/` stays intact as the **quarry + machinery source** through all of B2/B3 + Phase C. At the very
end (after Phase C hardening), the **culminating step on the branch** is: **replace `research/pde_ledger/` with the
new build in ONE atomic commit** (and delete the transient `_v2` dir). Then **merge the branch into `master` with a
merge commit (`--no-ff`)** so master's history records the whole development leading up to the overwrite.

- The `_v2` string lives **only** in the transient build-folder name — it never appears in stage stems, LaTeX, or the
  final path. **Final canonical path = `research/pde_ledger/`.** Rationale: git history + Zenodo **concept-DOI**
  versioning both preserve v1; a single canonical path keeps README/repo references valid; a permanent `_v2` suffix
  would be confusing cruft. (Zenodo publishes as a new *version* of the existing ledger DOI — transparent lineage.)
- **⚠ Pre-swap gate (Phase-C exit):** confirm **`origin/master` contains the full old ledger** (it is the "nothing
  lost" guarantee) *before* the overwrite commit. Currently `master` is committed but **not pushed** — push is a
  hard prerequisite of the swap. (No pushing without user ask; flagged here as the gate.)

**Machinery to COPY into `research/pde_ledger_v2/` (≈45 infra files; from the machinery inventory):**
- **Stage triplet template** — one representative note + `*_sympy_audit.py` + `*_mathematica_audit.wl` + paper card
  `stage_NNN.tex`, kept as an empty/example template.
- **Runners** — `scripts/run_all_audits.sh`, `mathematica/run_all_audits.sh`, `mathematica/run_one_audit.sh`
  (glob → invoke → write `output/{stem}.txt` → stamp `EXIT_CODE:`), plus the `output/` dirs. **⚠ B2 retargets the
  runner glob patterns + usage/help examples** from `moving_throat_pde_*` → the new stem (`scripts/run_all_audits.sh:44`,
  `mathematica/run_all_audits.sh:41`, `mathematica/run_one_audit.sh:30`).
- **`notes/MATHEMATICA_MIRROR_POLICY.md`** — the independent-route rule (§5).
- **LaTeX build chain** — `paper/pde_ledger.tex` (→ retitle, keep filename generic), `document_setup.tex`,
  `macros.tex` (`\stagefield`/`\claimstatus`/`\resultanchor`/`\StatusExact|Reduced|Open`/`\StageFile`),
  `main_parts.tex` skeleton, `frontmatter/00–05`, infra appendices (`reproducibility_map`, `source_file_index`,
  `parameter_provenance_audit`, `fill_workflow`), and empty `paper/parts/`, `paper/appendices/stage_appendix_partNN.tex`,
  `paper/stages/`.
- **Provenance/coverage scaffolds, RESET to zero** — `notes/STAGE_PROVENANCE_INDEX.md` (empty table + headers),
  `notes/STAGE_VERIFICATION_COVERAGE.md` (headers + `canonical_stage_count: 0`),
  `notes/LINEAR_STAGE_RENUMBERING_MANIFEST.json` (`canonical_stage_count: 0`), redteam `MANIFEST.yaml` stubs.
- **Do NOT copy (content):** the 253 stage notes/scripts/cards, the 8 physics `parts/part0N_*.tex`, redteam findings,
  `notes/em_projected/`, `materials.zip`.

**Stem/naming (decision D, resolved):** clean, reader-facing **`ledger_stageNNN_{slug}`**, globally numbered,
grouped into Parts in the LaTeX/appendix layer (as the old ledger grouped `part01 = stages 001–023`). No `v2`, no
`moving_throat_pde_`. Sector is carried by the Part grouping + provenance index; paper-extraction pulls by
Part/stage-range.

---

## 2. The v2 structure — Part 0 + Parts I–VII (sector map)

| Part | Sector | Nature | Audited stages (est.) |
|---|---|---|---|
| **0** | Conceptual foundation (the physical picture) | narrative frontmatter (no scripts) | — (prose chapters) |
| **I** | The medium | assemble | ~4–6 |
| **II** | Gravity | assemble (bulk of the algebra) | ~20–28 |
| **III** | Light | assemble (earned-photon headline + departures) | ~4–6 |
| **IV** | Charge | assemble | ~2–3 |
| **V** | Magnetism | assemble | ~6–8 |
| **VI** | The knit | assemble | ~4–6 |
| **VII** | Integration | **genuinely new synthesis** (+ χ_Q reconciliation + registers) | ~5–8 |

**Total ≈ 45–65 audited stages** (vs 253) — the ~225-stage search tail + ~52 EM/charge-drift stages do NOT carry over.

**⭐ Granularity rule (decision A, resolved) — physics-driven atomic, NOT a target number.** Each stage = **one
atomic, independently-verifiable derivation step.** A gate that bundles separable sub-derivations decomposes into
several stages (e.g. Gate 4 `pathA_33` → {DtN Hankel-series fingerprint → prefactor algebra → χ_Q closure →
`μ̂₀`-free dim-check}); a self-contained gate stays one stage. This preserves all intermediate algebra (goal 1) and
makes each step citable for extracted papers (goal 3). **The exact atomic split per gate is finalized per-Part during
B3, shown to the user before that Part is built.** Guardrail against old-ledger sprawl: granularity tracks the
*derivation's* natural steps, never a search.

**Honest-scope convention (used in §3):** many earned gates carry a `FAIL_*` *top-line* whose value is the **earned
content inside** + a **characterized departure**. Each such stage is represented by its **earned headline** with the
FAIL as the honest landing (e.g. light = "the medium carries 2 transverse photons at `c_γ²=μ_R/ρ_br`" with
`FAIL_CAUCHY_STRAY_LONGITUDINAL` the characterized departure). No softening; the departure is first-class.

---

## 3. Per-Part carry-over list (gate → stage; atomic split finalized in B3)

Verdict tokens are read from each report's top-line. Scope class ∈ {EARNED, CALIBRATED, SIM-DEFERRED,
CHARACTERIZED-DEPARTURE}. "Reshape" = the §5 work each fold needs.

### Part 0 — Conceptual foundation (narrative frontmatter)
A faithful, self-contained statement of the physical picture (what the medium / brane / four sectors / defect *are*),
folded from `docs/conceptual_foundation.md` so the ledger stands alone (goal 2). Native-terms mnemonic: gravity=flow
(one-way drain), light=brane shear, charge & magnetism = the throat (static vs moving). Prose + figures; no scripts.

### Part I — The medium (~4–6 stages)
| Gate row | Source gate(s) | Verdict token | Scope | Notes |
|---|---|---|---|---|
| GNLS action + EOS (`P=Kρ⁵`, `c_s²=5Kρ⁴/m`) | pathA_19 (dim foundation) + pathA_20 (c_s) | RETAIN_L_T_M / C_GAMMA_RATIO_UNDERDETERMINED | EARNED (form) + CALIBRATED (λγ) | both `.wl`-only → **owe `_sympy.py`** when folded |
| Two-phase material-state ontology (order field χ_B) | `notes/brane_bulk_handoff.md` + `docs/conceptual_foundation.md` | ACTION_SPECIFIED_CLASSIFIED | EARNED (structure) | postulated microstructure, labeled. **RECONCILED 2026-07-07:** `pathA_23_stage0` (EM-native brane-elastic fork) is DROPPED to history per §4 — its MacCullagh menu was superseded by the pathA_35 G0 freeze; the two-phase ontology here is the physically-distinct χ_B material-state closure. Fresh-authored audit (dim + recovery-reduction + θ-as-φ no-go). See `research/pde_ledger_v2/notes/part1_medium_atomic_split.md`. |
| Order field + shear-surface freeze (LIVE) | pathA_35 (T0_SHEAR_FROZEN) | T0_SHEAR_FROZEN + SECOND_MEDIUM_DRIFT_AT_FREEZE(11) | CALIBRATED (freeze inputs) | LIVE freeze; `pathA_25_G0` density-smectic dimcheck = HISTORY (route closed), not a stage |

### Part II — Gravity (~20–28 stages)
| Gate row | Source gate | Verdict token | Scope | Reshape / notes |
|---|---|---|---|---|
| Drain → `1/r²` inter-defect force | pathA_21c | FORCE_ATTRACTIVE_DERIVED | EARNED (form/sign) + CALIBRATED (magnitude) | dual-engine (A2 added `_sympy.py`) |
| Brane↔bulk return, `1/r²` survives slab + return-residual radiation | pathA_29 | RETURN_RESIDUAL_PREDICTION | EARNED + falsifiable residual | dual-engine |
| Monopole/dipole radiation constraint spec | pathA_28 | MONOPOLE_DIPOLE_RETURN_CONDITIONAL | EARNED (constraint) | dual-engine (both `pathA_28_monopole_sympy.py` + `.wl` exit 0) |
| Gate 1 — frozen-wall D/N unit test | pathA_30 | DN_UNITTEST_BC_DEPENDENT | CALIBRATED (BC input) | dual-engine |
| Gate 2 — scalar breathing `η₀₀` | pathA_31 | BREATHING_CALIBRATED | CALIBRATED | dual-engine |
| Gate 3 — grouped-`P2` / ℓ=2 isotropy | pathA_32 | ISOTROPY_CALIBRATED | CALIBRATED + EARNED (SO(3) theorem) | dual-engine |
| Gate 4 — quadrupole `54/5` normalization | pathA_33 | QUAD_CALIBRATED | EARNED (`27` fingerprint) + CALIBRATED (`2/5·G`) | dual-engine; `G=GENUINE_BLOCKED`; decomposes into several atomic stages |
| Gate 5 — cross-ℓ unification | pathA_34 | FAIL_UNDERDETERMINED_NOT_PREDICTIVE / CROSS_L_RESIDUAL_PREDICTION | CHARACTERIZED-DEPARTURE | dual-engine; Gate-6 selector = sim-deferred input |
| Density-mode ℓ=2 radiative port (A1) | pathA_43 | DENSITY_PORT_HOSTED | EARNED (structure, vector-free) + CALIBRATED (magnitude) | dual-engine; retires the `A_w` scaffold |
| 2.5PN match-back (A3) | pathA_2_5pn_matchback | (consistency over calibrated moments) | CALIBRATED-consistency | dual-engine artifact exists |
| PN ladder match-**CITE** (1PN→4PN + 2.5PN) | `research/4d_*pn*` (DOI'd papers) | (audited external corpus) | CITE-only | **no re-derivation**; cite by DOI (README) |
| Gate 6 — branch realization | (the sim) | — | **SIM-DEFERRED** | CITED reconnaissance; the wall |

### Part III — Light (~4–6 stages)
| Gate row | Source gate | Verdict token | Scope | Notes |
|---|---|---|---|---|
| Light on the shear-surface brane (couple-stress structure) | pathA_35 gateL | FAIL_COUPLE_STRESS_NOGO | CHARACTERIZED-DEPARTURE | earned: brane carries light |
| `c_γ²=μ_R/ρ_br`, 2 transverse photons | pathA_36 | FAIL_CAUCHY_STRAY_LONGITUDINAL | EARNED (2 photons) + CHARACTERIZED-DEPARTURE (stray longitudinal) | dual-engine |

### Part IV — Charge (~2–3 stages)
| Gate row | Source gate | Verdict token | Scope | Notes |
|---|---|---|---|---|
| Throat-body electric `1/r²` localized Coulomb, `±w` sign, Goldstone `h` | pathA_38 | THROAT_ELECTRIC_LOCALIZED_COULOMB | EARNED (`p=2` localization, sign) + CALIBRATED (`Q_E`) + SIM-DEFERRED (parity/compactness) | dual-engine |

### Part V — Magnetism (~6–8 stages)
| Gate row | Source gate | Verdict token | Scope | Notes |
|---|---|---|---|---|
| Scalar-admixture screen (stage 0/1) | pathA_39_stage01 | FAIL_OBSERVABLE_SCALAR_ADMIXTURE | CHARACTERIZED-DEPARTURE + SIM-DEFERRED | robust import-forced `h` floor (exact isolated Maxwell excluded) + density-admixture upgrade rests on sim-deferred `a_L≠0`/`q_L` |
| Magnetic force (stage 2) | pathA_39_magnetic | MAGNETIC_FORCE_DERIVED | EARNED (sign) + SIM-DEFERRED (`aT/aL`) | dual-engine |
| Operator parity under motion (stage 3) | pathA_39_stage3 | FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING | CHARACTERIZED-DEPARTURE | dual-engine |
| Field-coupling classification (stage 4) | pathA_39_stage4 | primary=FIELD_SCALAR_VECTOR_DEPARTURE (flags: scalar_sector_stable, density_charge_coupled, operator_parity_contamination) | EARNED (DOF count) + CHARACTERIZED-DEPARTURE | dual-engine; sector close |

### Part VI — The knit (~4–6 stages)
| Gate row | Source gate | Verdict token | Scope | Notes |
|---|---|---|---|---|
| Cone-lock `λγ=1`, `c_E=c_γ` | pathA_40 | CONE_LOCK_CALIBRATED | CALIBRATED (both locks are inputs, `Δr=2`) | dual-engine |
| NG5 second-medium drift | pathA_41 | SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) | CHARACTERIZED-DEPARTURE (incomplete 4D→3D reduction) | dual-engine |
| Charge-coupled scalar map | pathA_42 | SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED | CHARACTERIZED-DEPARTURE + SIM-DEFERRED (magnitude) | dual-engine |

### Part VII — Integration (~5–8 stages) — **genuinely new synthesis**
| Gate row | Content | Scope | Notes |
|---|---|---|---|
| Unified frozen parent action + BCs | the one action all sectors reduce from | new synthesis | full gauntlet |
| The calibration map | every postulate/constant labeled DERIVED/INPUT/gap/benchmark | assemble (cite `decisions/14`) | audit all constants accounted |
| Whole-system dimensional-firewall check | one dim-check across the assembled action | new | able-to-fail, units-restored |
| **The registers** (open-items / no-go / calibration) | every sim-dependent quantity + exact eq + confirming measurement; every characterized departure; every no-go — a first-class, permanent register (completeness standard #4) | new | **the falsification-first record; nothing lost** |
| The simulation hand-off spec | the Gate-6 / sim-deferred deliverable derived from the registers | new | *the literal deliverable* |
| χ_Q reconciliation | `pathA_22b ≈0.712` vs `pathA_33 =1` — same name, different computations | **new derivation work** | tracked open item; may need real algebra |

---

## 4. What is DROPPED (quarry for stray lemmas only)

- **EM import `004–020`** (17 stages) — the fundamental projected-Maxwell block; superseded (our EM is emergent, NOT
  exact Maxwell; the `S_EM[A]` scaffold is replaceable bookkeeping — diagnostic-confirmed).
- **The ~225-stage branch-selection/optimization SEARCH tail** — non-terminating Gate-6 reconnaissance; its converged
  result is the sim-deferred citation (Part II), not the search.
- **Charge-probe drift `219–253`** (35 stages).
- **Superseded pathA foundations** — only `pathA_19`/`pathA_20` are carried (Part-I dim/velocity foundation);
  `pathA_20b/21/21b/22a/22b` are reference-only. `pathA_23` (EM-native), `pathA_24` (little-arrows), `pathA_25`
  density-smectic route-closure (incl. its `G0_DIMCHECK_PASS` — the LIVE order-field freeze is `pathA_35`'s
  `T0_SHEAR_FROZEN`, **not** `pathA_25`), `pathA_26` (drain Derrick), `pathA_37` (retired flow-gate) all stay as
  **history** (`docs/conceptual_history.md`), not stages.

---

## 5. The reshape spec (per gate → stage) — the load-bearing per-fold work

Our `stage1_solver/tools/pathA_NN_*` pairs violate the ledger contract two ways: (1) `argparse --compare` + they
write JSON/YAML payloads; (2) the `.wl` is a **payload mirror** of the SymPy (exactly what `MATHEMATICA_MIRROR_POLICY`
rejects). Each fold requires:

1. **Strip** the `argparse`/`--compare` harness and **all file-writing** (no JSON/YAML payloads).
2. Make **each engine standalone**: print-only, **assert-zero-residual**, `raise`/`Exit[1]` on failure — ledger
   idioms (`banner`/`subbanner`/`expect_zero`/`pass`/`fail`).
3. **Re-author each `.wl` as a genuinely independent route** — different decomposition + native Wolfram primitives
   (`Solve`, `Reduce[ForAll]`, `Series`, `Residue`, `Inverse`, `NIntegrate`), **not** a transliteration of the `.py`.
   A genuine disagreement escalates to Claude+Codex ([[feedback-claude-codex-resolve-math]]).
4. **Dual-engine both exit 0** via the runners.
5. **Register** the stage: row in `STAGE_PROVENANCE_INDEX.md`, `pending`→`verified` in `STAGE_VERIFICATION_COVERAGE.md`
   (bump denominator), bump `canonical_stage_count`.
6. **Author the note + card as SELF-CONTAINED (completeness standard #1, HARD STANDARD).** The stage note
   (`notes/stages/{stem}.md`) carries the **full derivation inline** — a reader never has to open a `pathA_*` report
   to follow it. The `pathA_*` report is cited as *provenance*, not as the content. Decimals quarantined to a labeled
   benchmark slice.
7. **Update the parameter register (`research/pde_ledger_v2/notes/parameter_register.md`) — MANDATORY every stage, never
   skip.** Add the stage's new parameters (dimension `{L,T,M}` + provenance class) and every new relation edge
   (`DERIVED`/`IMPOSED`/`PENDING`/`CLOSED-NEG`/`CODIM-PROVEN`); record any reduction route the stage *discharges*
   (debt → DERIVED) or *closes* (→ CLOSED-NEG). This is the running seed for Part VII's calibration map + the irreducible
   codimension count, and it doubles as the cross-stage dimensional ledger (the whole-system dim check at Part VII then
   confirms rather than discovers). The register tracks *nominal vs irreducible* knobs relationally — a knob can be a
   manifestation of others.

**Dual-engine gaps to close during fold:** only the folded Part-I foundations **`pathA_19` and `pathA_20`** are
`.wl`-only and owe a `_sympy.py`. Everything else carried is already dual-engine (incl. `pathA_28`). `pathA_20b/21/
21b/22a/22b` are reference-only (not folded → owe nothing); `pathA_24`/`pathA_37` are history.

**Division of labor** ([[feedback-claude-reviews-codex-codes]], [[feedback-codex-is-fix-applier]]): **Codex** reshapes
the scripts + re-authors the independent `.wl` + runs to exit 0; **the orchestrator** authors the self-contained
note/card + registration + reviews.

---

## 6. Per-stage verification protocol (decision B, resolved)

- **Parts I–VI (repackaging of already-earned+tri-reviewed physics):**
  - *Pre-exec (lighter):* Claude drafts the reshape directive → **Codex design-review** (xhigh) → fold → Codex
    confirm. **No GLM tertiary** for a repackaging (no new physics verdict).
  - *Exec:* dual-engine, **both exit 0** via independent routes.
  - *Post-exec — the MANDATORY tri-review, all three, every reshaped script, never optional*
    ([[feedback-review-agents]]): (1) orchestrator **arbiter re-run** (unchanged scripts), (2) **fidelity audit**
    (fresh agent — the reshaped step verifies the *same* claims as the committed report; no dropped check/
    transliteration error), (3) **adversarial-with-ablation** (fresh agent). **The adversarial leg is scoped to
    reshape-integrity** — it hunts a reshape-introduced rig (a `.wl` that is secretly a mirror, dropped able-to-fail
    teeth, sneaked-in hardcoding, pass-by-construction), NOT re-litigating the already-earned physics. Uniform across
    every stage (goal 2). An in-script mutation probe does **not** substitute for leg (3).
- **Part VII (genuinely-new synthesis) + the χ_Q reconciliation:** full new-derivation gauntlet — directive → **Codex
  design-review → GLM-5.2 tertiary → Codex re-green** → dual-engine → full tri-review → user gate.

**Pilot-then-batch (decision C, resolved) — a TWO-STAGE pilot, user-gated after each:**
1. **II.1 `pathA_21c`** — the clean happy-path: locks the mechanical template (reshaped scripts + independent `.wl` +
   self-contained note + card + registration + the LaTeX/Part scaffolding). Already dual-engine.
2. **III.2 `pathA_36`** — a **FAIL-headline** stage: proves the earned-content-with-characterized-departure
   representation before batching (Parts III/V/VI are full of these).

Then batch Part-by-Part; per-batch user gate (sequential audit chunks, [[feedback-sequential-audit-chunks]]).

---

## 7. Build order

1. **B2** — stand up `research/pde_ledger_v2/` alongside the old ledger: copy the machinery (§1), reset counts to 0,
   retarget the runner globs/help-text to the new stem, verify the empty LaTeX build compiles + the runners run clean
   on the template.
2. **B3** — assemble in dependency order: **Part 0 (conceptual) → I medium → II gravity → III light → IV charge →
   V magnetism → VI knit → VII integration.**
   - **Pilot first:** II.1 `pathA_21c` → III.2 `pathA_36` (user-gate each), THEN batch.
   - Finalize each Part's atomic-stage split (the §2 granularity rule) and show the user before building that Part.
   - Part VII last (it needs I–VI assembled to synthesize the parent action + calibration map + registers).
3. **Phase C** — redteam the built v2 (hardening). Then, at the **pre-swap gate** (§1: confirm `origin/master` has the
   old ledger), **overwrite `research/pde_ledger/` with the build in one atomic commit** on the branch, delete the
   transient `_v2` dir, and **merge `ledger-v2-rebuild` → `master` with `--no-ff`.**

---

## 8. Tracked integration items (surfaced now so they aren't lost)

- **χ_Q reconciliation** (Part VII): `pathA_22b ≈0.712` (older minimal-combination context) vs `pathA_33 =1`
  (outgoing-DtN Hankel context) — same name, different computations. Reconcile explicitly; do **not** silently merge.
- **The registers** (Part VII, completeness standard #4): every sim-dependent quantity
  (`{L/a, ε per ℓ, stability, Z0_ret/Z1_ret selector, a_L, C_hu, M_h/c_E, aT/aL, …}`) + its exact equation + the
  confirming measurement; every characterized departure; every no-go — a permanent falsification-first record.
- **The calibration map** (Part VII): every postulate/constant labeled; cross-checked against
  `decisions/14_value_provenance_and_calibration_map.md`. **Running seed = `research/pde_ledger_v2/notes/parameter_register.md`**
  (the living knobs/dimensions/reductions register, updated every stage per §5 step 7). Part VII's job is the **irreducible
  codimension count** over it (the `Δr=2` Krull-dimension technique from pathA_40 scaled to the full constraint variety) —
  nominal knobs minus derived/convention/discharged-reduction = the predictive-surplus denominator.

---

## 9. Design decisions — RESOLVED (user, 2026-07-07)

- **A — Granularity:** physics-driven atomic (one independently-verifiable step per stage), ≈45–65 stages; per-Part
  atomic split finalized in B3 and shown to the user before building that Part.
- **B — Verification depth:** uniform full tri-review every reshaped stage, with the **adversarial leg scoped to
  reshape-integrity**; Part VII gets the full new-derivation gauntlet.
- **C — Pilot:** two-stage — II.1 `pathA_21c` (happy path) → III.2 `pathA_36` (FAIL-headline), user-gate each.
- **D — Path/naming:** build in transient `research/pde_ledger_v2/` alongside the old ledger; **overwrite
  `research/pde_ledger/` in one atomic commit at the end** (git history + Zenodo concept-DOI preserve v1; pre-swap
  push gate). Clean reader-facing stem `ledger_stageNNN_{slug}`; no `v2`, no `moving_throat_pde_`.
- **Completeness #1 — self-contained notes:** HARD STANDARD (§5.6).
- **Completeness #2 — conceptual Part 0:** included (§3 Part 0).
- **Completeness #3 — PN corpus:** cited by DOI (README), not inlined (the `4d_*pn*` papers are DOI'd + in `research/`).
- **Completeness #4 — open-items/no-go/calibration registers:** first-class in Part VII.

---

## 10. Size estimate + guardrails (restated)

- **≈ 45–65 audited stages vs 253** — the ~225-stage search tail + ~52 EM/charge-drift stages do **not** carry over.
- Completes the **SPEC, not the proof**; **Gate 6 stays sim-deferred**; a no-go is still possible (esp. the
  characterized departures in Parts III/V/VI and the χ_Q reconciliation).
- **v1 preserved by git history + Zenodo concept-DOI versioning;** the overwrite happens once, atomically, only after
  Phase C and the pre-swap push gate.

---

## 11. References
- Plan: `notes/ledger_v2_rebuild_plan.md`. Front door: `STATUS.md`. Live state: `decisions/13` §0.
- Machinery source (quarry): `research/pde_ledger/` — stage triplet + `notes/MATHEMATICA_MIRROR_POLICY.md` +
  LaTeX build (`paper/pde_ledger.tex` → `document_setup`/`macros`/`main_parts`/`parts`/`appendices`/`stages`) +
  `notes/STAGE_PROVENANCE_INDEX.md` / `STAGE_VERIFICATION_COVERAGE.md` / `LINEAR_STAGE_RENUMBERING_MANIFEST.json` +
  runners.
- Earned gates: `software/stage1_solver/{tools,reports,directives}/pathA_{21c,28,29,30–36,38,39,40,41,42,43}_*` +
  `tools/pathA_2_5pn_matchback_*` + `notes/pathA_43_2_5pn_matchback.md`.
- PN corpus (CITE by DOI, do not re-derive): `research/4d_{1pn_full,1pn_bridge,2pn,2_5pn,3pn,4pn}`,
  `research/1pn_orbital_dynamics` (DOIs in `README.md`).
- Conceptual: `docs/conceptual_foundation.md` (→ Part 0). Full scope: `docs/development_plan.md`. Process:
  `docs/development_pipeline.md`.
