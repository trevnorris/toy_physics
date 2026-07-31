# NEXT SESSION — v3 handoff

**Branch `ledger-v2-rebuild`.** ⚠ Run `git log --oneline -5` and `git status` first. ⛔ Do not trust a
hash written in any doc, including this one.

---

## ⛔ ORIENTATION BUDGET — read these four, in order, and nothing else

1. `research/pde_ledger_v3/CHARTER.md` — what v3 is, and what its scope **excludes**
2. `research/pde_ledger_v3/V3_STEP_PLAN.md` — ⭐ **§0 first**, then S0.5 → S2
3. `research/pde_ledger_v3/DEFECT_REGISTER.md` — the known error surface
4. `research/pde_ledger_v3/SESSION_REASONING.md` — how we got here

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
here is why, and here is what would make that wrong."* Ten pin-shaped identifications are on record;
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

**✅ v3 is OPEN** — charter, defect register, step plan, reasoning, techniques. Two review rounds folded.

**▶ NEXT: S0.5**, then S1, S2 — side by side.

### S0.5 is first and it is not optional

⛔ **Registry surgery before any step is banked.** The live seed already carries `Q.medium.c_gamma` and
`Q.medium.lambda_gamma` as *medium* quantities with `R3` `DERIVED-EXECUTED`. So S1 executed as-is would
**re-import the O-01 universe hole** — a quantity the walkthrough never introduced, counted as a medium
input. `c_γ`/`λγ` leave the medium counting contract until S9 / S20a; recompute acceptance; ⛔ never
preserve the old numbers.

⚠ A **substrate-action step** is also owed before S2 — S1 gives the EOS but not the GNLS action,
quantum-gradient term or Madelung balances, yet S4 invokes a "core balance" and S12 needs the
momentum/energy partners.

---

## ⛔ THE FOUR THINGS MOST LIKELY TO GO WRONG

1. **Treating a clean v3-gravity close as validation of the model.** It is not. Gravity is clean
   *because it is insulated from the interior*. See `CHARTER.md` §1.1 — and note the boundary is now
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

## ⭐ THE QUESTION THAT IS NOT IN SCOPE, AND WILL STILL BE OPEN

**What makes a muon a muon?** Mass differs 207× while charge is *exactly* identical. The falsified
support-only tower was the only family label the model had, and killing it left none. ⛔ v3-gravity will
not touch this, by construction — recorded so a debt list stays honest about what it is *not* paying
down.
