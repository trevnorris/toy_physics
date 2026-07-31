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

## ▶ WHERE WE ARE

**✅ Step ① is DONE** (`407eed94`). The `a`-pin is retired from everything that computes. Four parties
derived the new acceptance numbers independently, by four different methods: baseline `9/5/4`, C-M1
`9/6/3`, C-M2 replaced with a real entailment, C-M3 `9/4/5`. ⭐ `dim_after` is 5 on both sides — the
**Δ moving 5→4** is the signal.

**✅ v3 is OPEN** — charter, defect register, step plan, reasoning, techniques. Three review rounds folded (see below).

**▶ NEXT: S0.5 → S1 → S1.5 → S2** — side by side.

⚠ **Round 3 status: Grok CLEAN (empty blockers); Codex REGRESSED with five.** All five are folded but
**not re-reviewed**. ⛔ Treat S0.5–S2 as unverified and check as we go.

⛔ **S14a is new and BLOCKING for phase 4.** Correcting S12 to the committed dynamical `Γ_B` law severed
the chain to the imported gravity results: `Γ_B` conserves total material and converts only *order*,
while the imported force derivation assumes a Gauss drain removing *number flux*. Until a bridge derives
the projected order-loss source, `J`, the `J → Q` map and the `ω→0` return law, **S14, S15 and S20 are
CONDITIONAL, not carried.**

### S0.5 is first and it is not optional

⛔ **Registry surgery before any step is banked.** The live seed already carries `Q.medium.c_gamma` and
`Q.medium.lambda_gamma` as *medium* quantities with `R3` `DERIVED-EXECUTED`. So S1 executed as-is would
**re-import the O-01 universe hole** — a quantity the walkthrough never introduced, counted as a medium
input. `c_γ`/`λγ` leave the medium counting contract until S9 / S20a; recompute acceptance; ⛔ never
preserve the old numbers.

⚠ **S1.5, the substrate-action step, runs AFTER S1 and before S2** (it consumes S1's primitives) — S1 gives the EOS but not the GNLS action,
quantum-gradient term or Madelung balances, yet S4 invokes a "core balance" and S12 needs the
momentum/energy partners.

---

## ⛔⛔ FIRST TASK — six residuals, carried deliberately, NOT folded

Nine review rounds were run (Codex + Grok, iterating). **Grok returned CLEAN five times**, empty blocker
list. **Codex is still finding items** — these six. They are recorded rather than half-folded because the
session ran out of context, and a half-fold is what caused three of the previous rounds' findings.

1. ⭐⭐ **`stage005`'s registry loci are now SILENTLY WRONG — and this is a live demonstration of F4.**
   The `WHOLLY HISTORICAL` banner I added shifted `stage005`'s line numbers, so
   `reduction/quantities.yaml`'s loci (`:100, :104, :157, :161, :176 …`, ranges `73-88`, `68-90`,
   `189-203`) now point at the wrong text. ⛔ **Every gate still passes** — `_validate_loci` only checks
   the range *fits inside the file*. ⇒ Fix the loci, and treat this as the argument for fixing **F4**
   itself.
2. **S2's live registry provenance points into an artifact now marked do-not-consume.** `c_s0`,
   `c_gamma`, `lambda_gamma` all cite `stage005`. Re-point them as part of **S0.5**.
3. **S16's quotation** — the source block is now restored to the paper's actual (dimensionally loose)
   text, with the typo recorded. ⚠ Verify the repaired form is cited *separately* and never as a quote.
4. ✅ **DONE — and the residue was in TWO files, not one.** `docs/derivation_walkthrough_plan.md`'s
   phase-6 table cell and its replacement rule carry **only** the corrected input-vs-consequence form
   (`c13f9329`); **S21** in `V3_STEP_PLAN.md`, which still stated the **superseded first** replacement,
   now quotes the corrected rule. ⭐ Fixing the governing doc alone would have left the step that
   *executes* the knit still forbidding exactly what the corrected rule expects it to produce.
5. **`model_map.md` / control plane** — pre-correction assertions may remain below the new front-pointer.
6. **A13's branch choice does not propagate through S6** (the kink step) — it is gated at S1/S5/S12 only.

⇒ ⭐ **Do these first, side by side, before S0.5.** They are small, and every one of them is residue from
a correction that was *appended rather than applied*.

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
