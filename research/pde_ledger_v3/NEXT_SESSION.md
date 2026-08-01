# NEXT SESSION — v3 handoff

**Branch `ledger-v2-rebuild`.** ⚠ Run `git log --oneline -5` and `git status` first. ⛔ Do not trust a
hash written in any doc, including this one.

---

## ⛔ ORIENTATION BUDGET — read these five, in order, and nothing else

1. ⭐ `docs/model_map.md` — **READ-FIRST master map**: the throughline, the per-sector atlas, and the
   earned/calibrated/R1/departure ledger
2. `research/pde_ledger_v3/CHARTER.md` — what v3 is, and what its scope **excludes**
3. `research/pde_ledger_v3/V3_STEP_PLAN.md` — ⭐ **§0 first**, then S0.5 → S2
4. `research/pde_ledger_v3/DEFECT_REGISTER.md` — the known error surface
5. `research/pde_ledger_v3/SESSION_REASONING.md` — how we got here

⭐ **Why (1) was added (user decision, 2026-07-31):** this budget's omission of the master map is the
measured reason the model's recorded throat-support mechanism (trapped standing-wave pressure) was
invisible to v3 work, and the scope-widening argument itself came from `docs/model_map.md#shared-r1-throat-solve`.

⚠ Still a **bounded** list — five documents, not licence to read the corpus.
⛔ **Do not open** the v2 census, manifests, `archive/`, or the 43 stage notes "to get oriented".
Delegate any other read as a **specific question**, and require **≤20 lines + `file:line` loci** back.

⭐ `TECHNIQUES_THAT_WORKED.md` is worth ten minutes — fifteen moves, each with what it actually caught.

---

## ⭐⭐ THE ONE RULE THAT CHANGED — read this before doing anything

**Every step is walked SIDE BY SIDE with the user.** Not a stop-for-approval at the end — the user
participates in the derivation.

⛔ **Do not pre-derive a step and present conclusions.** Derive in the open, one move at a time, with the
reasoning stated **before** the result. Stop at every substantive move, not just step boundaries.

⭐⭐ **Flag every identification BEFORE making it:** *"I am about to treat X and Y as the same quantity —
here is why, and here is what would make that wrong."* **Thirteen** pin-shaped identifications are on record — one of which I created *while writing the register that warns about them*;
that move is what creates them.

**Why:** v2 was delegated heavily to go fast. That is the measured reason three structural findings —
the `a`-pin, the falsified lepton tower, the Spin Problem — went unnoticed for months by the only
physics reviewer this project has. ⚠ Slower is the point. ⛔ Do not optimise it back, and ⛔ do not read
it as licence to build apparatus instead.

Sub-agents: **lookup only, not judgement.** Codex and the external legs come **after** we have walked a
step — to check our work, not to do it.

---

## ⛔⛔ WHAT A FRESH SESSION GETS WRONG — measured, 2026-07-31

⛔ These are **not** hypothetical: each was actually believed and acted on during one session, by an
orchestrator who had read the orientation budget.

1. ⛔ **Wrong: *"the throat's support mechanism is an open question."***
   ⭐ Right: the throat is held open by **outward trapped standing-wave pressure** against inward brane
   tension plus superfluid backpressure, the trapped mode being a **brane-shear standing wave** with
   Wheeler-**geon** ancestry. ⚠ Contested by a later drain-flow reading, and the one computed result used
   a scalar surrogate. → `docs/conceptual_foundation.md` §4 (the 2026-06-24 sharpening);
   `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md` §16(B); and **PHASE 4b step
   Q2** in `V3_STEP_PLAN.md`, which now carries all three readings.

2. ⛔ **Wrong: *"the test is whether this matches real superfluid physics."***
   ⭐ Right: this is an **ANALOG**. The medium's properties are postulated as needed; departing from
   real-superfluid behaviour is expected and is ⛔ **not** a defect. Only an internal impossibility
   falsifies. → `CHARTER.md#falsification-standard`.

3. ⛔ **Wrong: *"`FAIL_CAUCHY_STRAY_LONGITUDINAL` is a defect"*** — and ⛔ equally wrong, *"the stray
   longitudinal becomes the `±w` displacement, i.e. charge."*
   ⭐ Right: the mode is **real and irremovable**, and its role was **CATALYTIC, ⛔ not constitutive**. Its
   irremovability **falsified the exact-Maxwell target** — exact Maxwell takes charge as an **external
   source put in by hand**, so matching it exactly would have inherited that gap — and thereby **forced
   the rethink** that produced the `±w` brane-tug mechanism. ⛔ **Different objects:** `u_L` is in-plane
   longitudinal brane displacement, `±w` is the throat's normal/orientation direction, and the charge the
   rethink produced is the **`h`-branon**, recorded *"(distinct from `u_L`)"*. ⚠ The token's name is
   misleading. → **S11** in `V3_STEP_PLAN.md` (the full history and the three distinctness loci).

4. ⛔ **Wrong: *"`r_cone = 9/2` conflicts with `λγ ≈ 1`, so the model is falsified."***
   ⭐ Right: **different objects.** `r_cone = c_E²/c_γ²` is **Lock B**; `λγ = c_γ/c_s` is **Lock A**. And
   `9/2` is a **witness constructed to demonstrate non-entailment**, ⛔ not a result.
   → row `A-CAUGHT` in `DEFECT_REGISTER.md`; `research/pde_ledger_v2/paper/stages/stage_040.tex:88-90`.

5. ⛔ **Wrong: *"whether one medium can carry both a longitudinal and a transverse mode is an open risk."***
   ⭐ Right: **computed, and it can** — the `μ_br > 0` branch carries two transverse modes and a
   longitudinal mode simultaneously and consistently. The obstruction is to **suppressing** the
   longitudinal one. → `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md` §15
   (the trilemma / pincer).

6. ⛔ **Wrong: *"B2's `FAIL_COUPLE_STRESS_NOGO` means light rests on a falsified foundation."***
   ⭐ Right: it rules out only **deriving** `μ_R` from the polar field `P`; light stands on the bare
   postulated modulus. → `DEFECT_REGISTER.md#B2`;
   `research/pde_ledger_v2/notes/stage030_pathA35_gateL_source_map.md:240-244`.

7. ⛔ **Wrong: *"`λγ` is a held-out test the model can fail."***
   ⭐ Right: it is currently a **calibration** — `CONE_LOCK_CALIBRATED`, Lock A untied. A calibrated value
   cannot falsify; it was fitted. Only a **derived** value ≠ 1 would. ⚠ It is still a **knob spent**, and
   `λγ³` enters the gravity normalization. → **S20a** in `V3_STEP_PLAN.md`;
   `docs/medium_requirements_and_prior_art.md:181` (Gate K).

8. ⛔ **Wrong: *"`R1` means one thing."***
   ⭐ Right: **three.** An **edge** in the parameter register (`:268`, the sound-speed derivation); **rung
   1 of the simulation ladder** (the one nonlinear throat solve) in `R1_REQUIRED(...)` (`:329`); and "the
   shared R1" in the model map. ⚠ Two of the three are in the same file.
   → `research/pde_ledger_v2/notes/parameter_register.md:268` and `:329`;
   `docs/model_map.md#shared-r1-throat-solve`.

9. ⛔ **Wrong: *"the orientation budget is sufficient."***
   ⭐ Right: it omitted `docs/model_map.md`, the READ-FIRST master map — which is why item 1 was missed,
   and where the scope-widening argument came from. ⚠ Now item 1 of the budget. → the budget section
   above.

10. ⛔ **Wrong: *"the launch prompt supersedes this handoff."***
    ⭐ Right: this session's launch prompt was a **lossy hand-compression** written after the handoff but
    from an earlier state — it reported round 2 when nine rounds had run, and omitted the carried
    residuals entirely. ⇒ ⭐ **Where a prompt and this file disagree, reconstruct the reasoning; ⛔ do not
    assume the later text is the more current one.**

11. ⛔ **Wrong: *"a review round's findings can be folded without re-reviewing."***
    ⭐ Right: measured twice this session — round 1's fold introduced **three** new physics errors that
    round 2 caught, and round 2's fold introduced more. ⇒ **A fold is not self-certifying.**

⭐ The common shape across items **2, 3, 6 and 7** is **importing a standard the project does not use**,
and across **3 and 4** it is **a surface match between different objects** (item **3** is both: an
imported exact-Maxwell target *and* a weld of `u_L` onto `±w`) ⇒ the questions that catch them are
*"what test does this project actually apply?"* and *"is this the same KIND of thing?"*

---

## ▶ WHERE WE ARE

**✅ Step ① is DONE** (`407eed94`). The `a`-pin is retired from everything that computes. Four parties
derived the new acceptance numbers independently, by four different methods: baseline `9/5/4`, C-M1
`9/6/3`, C-M2 replaced with a real entailment, C-M3 `9/4/5`. ⭐ `dim_after` is 5 on both sides — the
**Δ moving 5→4** is the signal.

**✅ v3 is OPEN** — charter, defect register, step plan, reasoning, techniques.

⭐⭐ **S0.5 was NOT reached.** The last session bought **corrections** instead — six commits, none of them
a walkthrough step:

- `c13f9329` — ⭐ **scope WIDENED to include charge and magnetism** (user decision). The charter's stated
  reason for deferring them did not survive checking; the stronger argument is arithmetic — **one**
  nonlinear throat solve is the shared `R1` for gravity `{μ_R, ρ_br}`, electric `bc_selection` **and**
  magnetism `q_T`, so a deliverable naming some siblings and omitting others is **incomplete, not
  scoped** (`docs/model_map.md#shared-r1-throat-solve`). **PHASE 4b (Q1–Q7)** is the new material.
  ⭐ Same commit corrected the **phase-6 knit rule** to **input-vs-consequence**
  (`docs/derivation_walkthrough_plan.md#knit-rule`).
- `b78dba6e` — **round-1 review folded: 8 blocking.** ⛔ The governing failure is this project's recurring
  one — the widening had been **APPENDED, not APPLIED**. Acceptance criterion ④ — the *definition of
  done* — still read *"the interior debts gravity depends on"*, so an executor could satisfy it with a
  gravity-only deliverable; and `STATUS.md`, the repository front door, still declared the old scope.
- `7e4dc71d` — **round-2 review folded: 6 blocking.** ⚠⚠ **Three of the six were errors the round-1 fold
  itself introduced** (an invented drain dependency, an over-sharpened A5, a false "four classes" in Q4).
  ⇒ This is the measurement behind item **11** above: a fold is not self-certifying.
- `3cfdc3c3` — the six stale `stage005` loci repaired (**residual #1**), the photon-sim track recorded,
  and `λγ`'s observational constraint recorded — ⚠ **which `18d92331` then CORRECTED**: `λγ` is a
  *calibration*, ⛔ not a test the model can fail (item **7** above).
- `931943d5` — ⭐ **name-based citation tags** for the files we own, plus `scripts/check_citation_tags.py`,
  **ablated per tooth** (renamed tag → `DANGLING CITATION`, duplicated tag → `DUPLICATE TAG`, bad path →
  `UNRESOLVED PATH`), each restored by exact inverse `sed` rather than `git checkout`.
- `18d92331` — the onboarding record above (the eleven-item list), ⭐ **CHARTER §1.3**
  (`CHARTER.md#falsification-standard`) stating for the first time what would count as **FAILURE**, and
  the correction of **a pin the orchestrator made**: the stray longitudinal is **CATALYTIC, ⛔ not
  constitutive** (item **3** above).

**▶ NEXT: S0.5 → S1 → S1.5 → S2** — side by side. ⛔ **None of it is started.**

⛔ **S14a is new and BLOCKING for phase 4.** Correcting S12 to the committed dynamical `Γ_B` law severed
the chain to the imported gravity results: `Γ_B` conserves total material and converts only *order*,
while the imported force derivation assumes a Gauss drain removing *number flux*. Until a bridge derives
the projected order-loss source, `J`, the `J → Q` map and the `ω→0` return law, **S14, S15 and S20 are
CONDITIONAL, not carried.**

### S0.5 is first and it is not optional

⛔ **Registry surgery before any step is banked.** The live seed already carries `Q.medium.c_gamma` and
`Q.medium.lambda_gamma` as *medium* quantities with `R3` `DERIVED-EXECUTED`. So S1 executed as-is would
**re-import the O-01 universe hole** — a quantity the walkthrough never introduced, counted as a medium
input.

⛔⛔ **DELETED, not retired** (user decision, 2026-07-31): `Q.medium.c_gamma`, `Q.medium.lambda_gamma`
**and `R3`** come out of the registry **outright**; `c_γ` is then introduced at **S9** and `λγ`/`R3` at
**S20a**, each with its own provenance. ⭐ The right word is **PREMATURE, ⛔ not "retired"** — they were
placed in the medium block *before the walkthrough reached the step that introduces them*. Recompute
acceptance; ⛔ never preserve the old numbers. ⚠ **This breaks code, in one pass:** `active_variables`,
`C-M1`, `able_to_fail.py` and `test_registry.py` all ride `R3` → **S0.5** in `V3_STEP_PLAN.md`.

⚠ **S1.5, the substrate-action step, runs AFTER S1 and before S2** (it consumes S1's primitives) — S1 gives the EOS but not the GNLS action,
quantum-gradient term or Madelung balances, yet S4 invokes a "core balance" and S12 needs the
momentum/energy partners.

---

## ⛔⛔ THE SIX RESIDUALS — re-checked against HEAD (3 ✅ · 1 deferral · 1 partial · 1 open)

Nine review rounds were run (Codex + Grok, iterating). **Grok returned CLEAN five times**, empty blocker
list. **Codex is still finding items** — these six. They were recorded rather than half-folded because the
session ran out of context, and a half-fold is what caused three of the previous rounds' findings.
⚠ **Status below is verified against HEAD**, ⛔ not carried from the previous handoff.

1. ✅ **DONE — `3cfdc3c3`.** The `WHOLLY HISTORICAL` banner on `stage005` shifted the file by **+5 lines**,
   so six `reduction/quantities.yaml` loci pointed at the wrong text while **every gate stayed green** —
   `_validate_loci` only checks that a range *fits inside the file* (**F4**). Corrected: `c_s0` `73-88 →
   78-93`, `c_gamma` `68-90 → 73-95`, `lambda_gamma` `189-203 → 194-208`, each verified by opening the new
   range. ⭐ **And the other 23 loci in `quantities.yaml` were swept** — nobody had checked them; all 16
   unique citations resolve, so the six were the only stale ones. **F4 itself is still OPEN**
   (`DEFECT_REGISTER.md` row F4) — see the pending queue below.
2. ⏸ **SPLITS — part moot, part DEFERRED, ⛔ nothing to fix today.** `c_s0`, `c_gamma` and `lambda_gamma`
   do all still cite `stage005`, which is now banner-marked do-not-consume. But `c_gamma`/`lambda_gamma`
   are **deleted outright at S0.5**, so their provenance is moot for them. ⛔ **`c_s0` stays live**, and its
   honest provenance is **S2's own record — which does not exist yet.** ⇒ Record this as an explicit
   **deferral to S2**, ⛔ not as a fix: re-pointing `c_s0` at anything before S2 derives it forward would
   invent a provenance.
3. ✅ **DONE — verified at HEAD.** **S16** in `V3_STEP_PLAN.md` block-quotes `4d_2_5pn.tex:605-614` **as
   written** (the remainder `O((a²/r²)·Φ_ext/r)` carries **no** mass factor, so force sits opposite
   acceleration), records that as a **typo**, and cites the dimensionally correct form **separately** at
   `4d_1pn_full.tex:886` — ⛔ never as a quotation. ⚠ The repaired form appears nowhere as a quote.
4. ✅ **DONE — and the residue was in TWO files, not one.** `docs/derivation_walkthrough_plan.md`'s
   phase-6 table cell and its replacement rule carry **only** the corrected input-vs-consequence form
   (`c13f9329`, `docs/derivation_walkthrough_plan.md#knit-rule`); **S21** in `V3_STEP_PLAN.md`, which still
   stated the **superseded first** replacement, now quotes the corrected rule. ⭐ Fixing the governing doc
   alone would have left the step that *executes* the knit still forbidding exactly what the corrected rule
   expects it to produce.
5. ⏸ **PARTLY DONE — the surplus list is marked, the ATLAS is not.** ✅ `docs/model_map.md`'s §4
   earned-surplus list now carries both `⛔ CONDITIONAL pending the S14a drain bridge, NOT earned` bullets,
   covering all **three** charter-conditional items (the `1/r²` exponent, the attractive sign, and the
   stage009 `RETURN_RESIDUAL_PREDICTION`) — `7e4dc71d`. ⛔ **Still open:** the pre-correction assertion
   survives **below the front-pointer** in the same file's derivation atlas, where `pathA_29` is still
   headed *"**held-out FALSIFIABLE PREDICTION** (`RETURN_RESIDUAL_PREDICTION`)"* with no S14a marking, and
   §2's sector table still states gravity's `1/r²` + attraction flat. ⇒ ⭐ **Sweep the whole file, not the
   ledger section** — that split is exactly the append-vs-apply failure. ⚠ The **count** in that list is a
   separate item and is ⛔ **not** an edit — see the pending queue.
6. ⛔ **OPEN — unchanged.** **A13's branch choice does not propagate through S6** (the kink step): it is
   gated at **S1**, **S5** and **S12** only, and **S6**'s `Register:` line names **A2** alone — while S6
   solves the EL equation for `δ` and `σ_wall` out of the very constants the A13 branch selects.

⇒ ⭐ **What is left here is #2's deferral, #5's atlas sweep and #6** — small, and both live ones are
residue from a correction that was *appended rather than applied*.

---

## ▶ PENDING WORK — queued, with the reason and the blocker

⭐ A **list**, ⛔ not machinery. Nothing here is started.

1. **Registry `anchor:` extension to `_validate_loci`** (the structural fix for **F4**). ⏸ **Deferred
   until AFTER S0.5, deliberately** — S0.5 deletes `Q.medium.c_gamma`, `Q.medium.lambda_gamma` and `R3`,
   so anchoring loci that are about to change is wasted work. **Blocker: S0.5.**
2. **The carried-results survey** — a **bounded** sweep for computed/earned results that `CHARTER.md` §3's
   *"what carries forward from v2"* list does not name. ⚠ It was held pending a clean fold; **the fold is
   now clean (`18d92331`), so it is UNBLOCKED.** ⭐ **Discipline: a list with loci, one pass.** ⛔ No
   schema, ⛔ no checker, ⛔ no manifest — that is exactly how the retired v2 census started.
3. ⭐ **The `FAIL_`-token audit** — ⚠ **a hypothesis from ONE confirmed instance, ⛔ NOT a finding.** Some
   `FAIL_*` / `*_NOGO` tokens may name **characterized DEPARTURES** rather than internal inconsistencies.
   ⭐ **The distinction that decides it** (`CHARTER.md#falsification-standard`): an **internal
   impossibility** — ghost mode, energy unbounded below, lost hyperbolicity, constraint-count failure —
   **versus producing MORE than an incomplete reference theory has a slot for**, which §1.3 says outright
   is ⛔ not failure. **`FAIL_UNSPECIFIED_SUBSTRUCTURE` is a candidate:** its own report says the quantity
   *"is therefore **not derived as zero and not derived as nonzero**"*
   (`software/stage1_solver/reports/pathA_23_stage2_constitutive_form.md:11`) — that is an
   **UNDETERMINED**, filed under `FAIL_`. ⛔ One instance is not a pattern; the audit is what would decide.
4. **`scripts/check_citation_tags.py` owes a review leg.** ⛔ The orchestrator **authored** it
   (`931943d5`), and **the builder is not the reviewer**. ⚠ It is per-tooth ablated but not
   independently reviewed. **Blocker: a fresh reviewer.**
5. ⚠ **`docs/model_map.md`'s earned-surplus COUNT.** The list is headed *"the earned surplus, **~6–7
   structural + 1 departure**"* (`docs/model_map.md#earned-surplus`), and the same figure recurs in the
   **Midway Knob Audit headline** further down the same section. **Three of its items are now
   CONDITIONAL** pending S14a. ⛔ **The count was DELIBERATELY not adjusted, and that is the correct
   state:** the bullets are ⛔ not one-to-one with counted items, and the method doc requires the **§7a
   closing certification** before any count is quoted (`docs/derivation_walkthrough_plan.md` §7a) —
   recomputing it here would be the exact move that rule forbids. ⇒ **It needs the certification, ⛔ not
   an edit.**
6. **The `R1` namespace collision** — three different objects share the name (item **8** above): the
   parameter-register **edge** `R1` (`parameter_register.md:268`, the sound-speed derivation), **rung 1 of
   the simulation ladder** in `R1_REQUIRED(...)` (`:329`), and *"the shared R1"* in the model map. ⚠ Two
   of the three are in the same file, and the two namespaces are **nested** — `R10`/`R30`/`R33` are
   *edges* that are themselves blocked on the *rung*. ⇒ **Its durable home is row `F5` in
   `DEFECT_REGISTER.md`**; the note inside **S0.5** (`V3_STEP_PLAN.md#s05-r-namespace-hazard`) is
   in-flight guidance only. ⛔ **Nothing has been fixed** — bare `R1` occurs **1295×** across **236**
   tracked files, so this is corpus-wide vocabulary, ⛔ not a small edit. ⭐ **Cheap mitigation, available
   now: never write bare `R1`** — write `R1_REQUIRED(...)` or "rung `R1`" in full at every use.

---

## ⛔ THE FOUR THINGS MOST LIKELY TO GO WRONG

1. **Treating a clean close of the gravity subledger as validation of the model.** It is not. Gravity is
   clean *because it is insulated from the interior*. See `CHARTER.md` §1.1 — and note the boundary is now
   **amended**: the worldtube result is **response-side** (how a body *moves in* a field), conditional
   on compactness, leading-order only, with a correction `O(a²/r²)` governed by the undetermined `a`.
2. **Rebuilding apparatus.** Two efforts died this way. The register is a **list**; `reduction/` already
   exists and must be **reused**, not rebuilt.
3. **Asserting absence from a partial read.** ⛔ This happened *this session*, on the headline
   conclusion: I read ~135 lines of a **3239-line** file and wrote "there is no surviving mass–radius
   slope". It was false. ⇒ `wc -l` before any universal negative.
4. **Re-deriving the PN ladder.** S17 is **cite-only**. ⚠ And `4d_2_5pn`/`4d_4pn` are *conditional*
   derivations by their own titles — do not carry them forward as more than that.

---

## ▶ DECISIONS OWED BY THE USER

1. **O-02** — is `K` + the EOS exponent one entry or two? (v2 left it open; the sim-input reading says
   two, the algebraic reading says one.)
2. **`ħ`'s class** — `postulated` (a property of the medium) or `calibrated` (an external constant we
   import)? In the sim-input set either way; the label moves it between tier 1 and tier 2.
3. **Scope confirmation** — the charter's boundary may not survive S16/S19 (the compactness premise vs
   the Spin Problem). If it doesn't, the honest move is to amend the charter, not preserve the boundary.
4. ⭐ **The `C-M1` replacement mutation.** `C-M1` **dies with `R3`**: its entire content is *"drop `R3`"*
   (`acceptance_check.py:54` builds it as `without_r3 + (lambda_gamma - lambda_gamma,)`, and
   `constraint_dimension` discards zero-valued constraints, so the appended term is **inert**). S0.5
   deletes `R3` ⇒ C-M1 has **nothing left to mutate**. Its replacement must be a **plausible physics
   error on a surviving relation** — `R1` (the sound speed), or `R2.xi_h` / `R2.h0` — ⛔ not a
   relabelling and ⛔ not another inert term. ⭐ **Which one is a USER call**, not the builder's.
   → **S0.5** in `V3_STEP_PLAN.md`.

---

## ▶ PHOTON-SIMULATION TRACK (user proposal, 2026-07-31) — ⛔ NOT a v3 walkthrough step

⛔ **This is a parallel research track, not a step in the v3 walkthrough.** It sits in no phase and gates
nothing.

**Why it is attractive.** Light is the one sector that does **not** need a defect. The corpus's deferral
language is specifically about the **throat interior / full 4+1D**, ⛔ **not** about the brane-shear
sector (`docs/two_throat_simulation_handoff_spec.md:671`, `:696`):

> *"The full object — a time-dependent, fully nonlinear, 4+1-dimensional solve of two moving throats
> radiating into an open medium — is **presently blocked** …"*; *"The **full** nonlinear two-throat sim
> stays deferred."*

**What already exists.** A working **linear** transverse-wave integrator —
`software/force_visualizer/physics/light.py:75-131` (explicit 2nd-order finite difference, 1D periodic
domain, 300 points). ⭐ It **ran**: measured `ω` at `k = 1–4` matched `c_γ·k` to **0.04%**, and both
transverse polarizations were recovered (`software/force_visualizer/output/verification_report.txt:42-48`).

⛔ **What is missing, and it is the whole job.** There is **no nonlinear shear equation anywhere in the
corpus**. Every form is quadratic, and
`research/pde_ledger_v2/notes/stages/ledger_stage034_transverse_move_action_row.md:147` records what
happens if you change that:

> *"(Breaking well-formedness — making the coupling nonlinear in `u_T`, non-local, or non-variational —
> fires the tooth.)"*

⚠ **You cannot reach it by "un-linearizing":** the quadratic Lagrangian was written **directly**, not
expanded from a nonlinear parent — which is **C6** (`docs/two_throat_simulation_handoff_spec.md:324`):

> *"**Closure status — no closed parent action exists.**"*

⭐ **The target is already a named open gate** (`atlas/nodes/open_gates/OPEN_NONLINEAR_S_SIGMA.md:63`):

> *"Promote or derive nonlinear throat action whose quadratic limit is S_eta^(2)."*

**Why it matters beyond the sim.** `docs/em_gravity_mined_verdicts.md:38` records that the trapped-shear
geon is *"intrinsically NONLINEAR — a linear/quadratic action can't establish it"*. ⇒ The geon (**C4**,
no equation) and the throat's holding-open mechanism are both downstream of this **same** missing
nonlinear action.

⛔ **Two traps, both already recorded:**

1. **Scale.** The linear equation needs only the **ratio** `c_γ² = μ_R/ρ_br`. Nonlinearity introduces an
   amplitude scale, so it needs `μ_R` and `ρ_br` **absolutely** — and both are calibrated/postulated
   (`[POSTULATE: stiffness value]`, `docs/two_throat_simulation_handoff_spec.md:70`).
2. ⛔ **Tautology.** If `μ_R` is supplied as an input, measuring `c_γ = √(μ_R/ρ_br)` returns the number
   that was put in. It becomes a real test **only if the brane's stiffness EMERGES** from the bulk medium
   rather than being imposed — and that is exactly where `pathA_35` returned `FAIL_COUPLE_STRESS_NOGO`
   (register **B2**, `DEFECT_REGISTER.md#B2`).

**Status:** recorded, ⛔ **not started**; ⚠ deliberately not begun before **S0.5**. ⇒ If it becomes real
work it earns its own document then.

---

## ⭐ THE QUESTION THAT IS NOT IN SCOPE, AND WILL STILL BE OPEN

**What makes a muon a muon?** Mass differs 207× while charge is *exactly* identical. The falsified
support-only tower was the only family label the model had, and killing it left none. ⛔ v3 will
not touch this, by construction — recorded so a debt list stays honest about what it is *not* paying
down.
