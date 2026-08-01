# Walkthrough decisions — user calls that shape what gets counted

Physics and modelling decisions made during the walkthrough. ⛔ Not derivations — those are the numbered
step records. This file exists so a decision is not silently re-litigated by a later session.

---

## D-01 — ⛔ The throat radius is NOT determined by the units pin (2026-07-30)

**Decision.** The `a`-pin does not determine the throat size, and treating it as though it does is wrong.
Two distinct things have been sharing one symbol:

| | what it is | status |
|---|---|---|
| **the pin `a`** | the nullity-1 residue of imposing four unit pins `{a, c_s0, ħ, m_GNLS}` on three base dimensions (`stage004:124-134`) | ⛔ a **units artifact**. Register: `CONV`, "never free", `A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT` |
| **the throat radius** | the physical size of a defect's throat | ⭐ a real physical quantity, ⛔ **not** obtainable this way |

⚠ **The collision is in the corpus and was recorded as a success.** `stage016:184` calls `a` *"the same
throat radius as stage018's — a physical-radius"*; `stage023:481` writes *"throat/pin radius"* as if they
are one thing; and `STATUS.md` logged *"`a` now groups AGREE across 016/018/023 … the same throat-radius
`L`"* as evidence of consistency.

⛔ **That agreement is DIMENSIONAL ONLY.** Two different lengths agree on being lengths, trivially. It is
a same-name-different-quantity collision at the foundation of the model, filed as confirmation. ⇒ This is
exactly what the plan's **check 4 (step-to-step identity, adjudicated not name-matched)** exists to catch.

### How the throat radius is actually fixed

**Now — calibration (the cheap way).** Calibrate the throat size against a **given flow rate of the
superfluid** and a **lepton mass**. ⇒ Class **`calibrated`**: chosen because the calibration is necessary
for the model to work. **TIER 2.**

**Eventually — derivation (the whole program).** ⭐ *"The throat radius is the end-all of this entire
program. Calculating the throat radius involves everything in our bag."* ⛔ Not available at this stage.

⭐ **Why it matters more than one row:** if the throat could be properly calculated, the payoff is
something like **the lepton ladder** — the mass sequence falling out rather than being fitted. ⚠ Recorded
as the aspiration it is, ⛔ **not** as a route, and ⛔ not as anything currently in reach.

### What this fixes in the count

| field | value |
|---|---|
| `is_tier` | **TIER 2 — calibrated** (flow rate + lepton mass) |
| `should_be_tier` | **TIER 3 — emergent** |
| `should_be_basis` | `physical-picture-expectation` — ⛔ **not** `named-route`; there is no executable route, only the program's goal |
| `delta` | **FLAGGED** |

⭐ **This is the single largest entry on the revisit list**, and the `is`/`should-be` pair is carrying
exactly the load it was designed for: the count is honest today (a calibrated input), and the gap says
where the whole program is pointed.

### ✅ Consequences — APPLIED 2026-07-31 (`407eed94`)

⛔ **This section is HISTORY, not instruction.** The pin was retired by removal: `R2.a_pin` and
`Q.medium.a_pin` deleted, `stage004` §4 rewritten to the complete three-pin basis, `parameter_register`
`:132`/`:269` corrected. ⛔ **Item 3 below — "the registry class is OPEN, do not resolve it" — is
SUPERSEDED**: there is no relation left to classify. Live workstream: `research/pde_ledger_v3/`.

<details><summary>original text, retained for provenance</summary>


1. ⛔ **The pin `a` and the throat radius must be separated** — different quantities, different symbols.
   ⚠ Touches `stage004`, `parameter_register.md:132`, and the `a` rows in 016/018/023.
2. **The pin is probably unnecessary.** Nothing consumes `a` in the medium block — it is the output of its
   own pin relation and appears in no `input_qids`. It is also `ξ_h/√2` (`stage004:133`), i.e. the healing
   scale up to a constant, and `ξ_h` is the one with a physical justification ("core balance") while `a` is
   explicitly a pin choice.
3. ⛔⛔ **`R2.a_pin`'s registry class is OPEN PENDING A USER DECISION, and the ambient count depends on
   it.** It was changed `CONVENTIONAL → DERIVED-EXECUTED` on a too-literal reading of the classification
   rule; a units-pin identity is ⛔ not "a defining equation in terms of other model quantities". ⚠ Note
   what this decision does and does not settle: it settles that **the throat radius** is calibrated, and
   it adds clause 4 below — it does ⛔ **not** by itself decide the pin relation's registry class, which
   is why that is listed here as a consequence *not yet applied*. ⛔ **Do not resolve it in passing.**
   ⚠ The registry data documents still carry `DERIVED-EXECUTED`; that is their present state, not a
   verdict, and they are deliberately left untouched until the call is made. If it reverts, the ambient
   count of 10 and the acceptance MATCH both need re-establishing.
   ⇒ **Aligned 2026-07-30 across every document that states a class:**
   `docs/derivation_walkthrough_plan.md` §1.1 + §8 (REOPENED) · `reduction/README.md` (UNRESOLVED, data
   as-is) · `_scratch/NEXT_SESSION.md` (open, and item 3 of OPEN FOR THE USER) ·
   `walkthrough/01_sound_speed.md` §A (its "both entries now say derived" flagged reopened) · here.
4. ⭐ **The classification rule gains a clause:** *a relation arising from imposing unit pins is not a
   defining equation*, however much it looks like one.

</details>

⚠ **`ħ` stays open — and ⛔ the retired pin is no longer the reason (updated 2026-08-01).** ⛔ **State it
precisely:** `ħ` is **not** unrelated to everything — `R2.xi_h` (`ξ_h = √2·ħ/(m_GNLS·c_s0)`) takes it as
an input. What is true is narrower: **no relation designates `ħ` as its OUTPUT**, so nothing derives it.
⚠⚠ **And note the shape:** `R2.xi_h` differs from the deleted pin only by a `√2` — **A7** records that
this is precisely what made the pin plausible for eleven months, so ⛔ do not read `R2.xi_h` as pinning
`ħ`; it designates `ξ_h`, and that direction is a **choice recorded in the registry
(`designated_output`), ⛔ not a derivation**. The relation is silent about which of the two is
fundamental. ⇒ The hypothesis that **`ħ` is a particle size in the medium** needs something *outside*
that relation to pin one of them, and nothing currently does. ⇒ Still not a bookkeeping gap; a real gap
in the model. What `ħ` does and does not touch: `research/pde_ledger_v3/DEFECT_REGISTER.md` row **C11**.

---

### D-01a — sequencing and blast radius (user, 2026-07-30)

**Order is fixed:** ✅ ~~archive the old apparatus~~ (**DONE** 2026-07-30 — `archive/`, see
`archive/README.md`) → ① **repair the `a`-pin damage** → ② resume the ledger walkthrough.
⛔ Do not resume phase 0 before ① is done. ⚠ Repair may require **revisiting
already-derived steps**, since anything that consumed the pin as a physical radius is suspect.

⛔⛔ **TRIAGE FIRST — a mention is not damage.** Many files reference `a` legitimately as a symbol. The
damage is only where one of these is true:
- **(a)** the pin is presented as *determining a physical radius*;
- **(b)** the pin and the throat radius are **conflated** as one quantity;
- **(c)** their dimensional coincidence was recorded as *evidence of consistency*.

⇒ Sort every hit into damaged / benign-mention **before** editing anything. Editing on a grep hit is how a
symbol-level change becomes a corpus-wide mess.

**Tier 1 — live v2 ledger, must be fixed**
| file | why |
|---|---|
| `notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md` + its `.py`/`.wl` | ⭐ the pin's **definition site** (`:124-134`) — start here |
| `notes/parameter_register.md:132` | the `a` (pin) row, `CONV`, "never free" |
| `notes/stages/ledger_stage016_l2_so3_covariance.md:184` | ⛔ **(b)** asserts `a` *"is the same throat radius as stage018's — a physical-radius"* |
| `notes/stages/ledger_stage023_nullspace_underdetermination.md:481` | ⛔ **(b)** *"throat/pin radius"* as one thing |
| `STATUS.md#position` | ⛔ **(c)** logs the three-stage grouping as agreement |
| `notes/stages/ledger_stage005_…md` + `notes/stage005_pathA20_source_map.md` | carry the anti-tautology caveat — ⭐ these are *correct*, keep and cite |
| `scripts/ledger_stage043_…_sympy_audit.py` + `.wl` | reference the pin in the count manifest |
| `reduction/` (`quantities.yaml`, `relations.yaml`) | our own `a_pin` quantity + `R2.a_pin`; ambient count depends on its class |

**Tier 2 — Path-A / solver cluster; scope decision needed before touching**
`software/stage1_solver/`: `decisions/13`, `decisions/14`, `reports/pathA_19`, `pathA_21`, `pathA_21b`,
`directives/pathA_18/19/20/20b/21`, `src/stage1_solver/dimensional_check.py`; plus
`docs/pathA_preregistration.md` (which already warns ⛔ *"do NOT rely on the `a=c_s=ħ=m=1` natural-unit
pins"* — ⭐ correct, keep).

**Tier 3 — EM-charge cluster, ~15 files** — `software/em_charge_attribute/u1_*`, `u2_*`. ⚠ Uses `a_pin`
heavily. ⛔ Scope unknown; do not touch without deciding whether that workstream is live.

**Tier 4 — out of scope** — `research/pde_ledger/redteam_adversarial/*` (the *old* ledger, not v2), and
`archive/census/pilot/*` (✅ archived 2026-07-30).

⭐ **Start at the definition site (`stage004`) and work outward.** Fixing consumers first leaves the
source of the error in place to re-infect them.

---

## OPEN ITEMS — recorded 2026-07-30, ⛔ not decided

⚠ These are **findings awaiting a call**, not decisions. They live here because they change what gets
counted, and a later session would otherwise re-litigate them from scratch.

### O-01 — ⭐ a universe hole in the seed: `c_γ` was never introduced by a step

**The finding.** `c_γ` is one of the six quantities the medium block says must be supplied
(`_scratch/NEXT_SESSION.md`, THE TWO SCRIPTS) — but ⛔ **no walkthrough step introduces it.** Step 1
(`00_medium_and_brane.md`) lists 4 scalars (`ħ`, `m_GNLS`, `K`, `ρ0`), 3 discrete choices, 3 fields and
1 BC, and `c_γ` is in none of them. Step 2 (`01_sound_speed.md`) records *"what's new: nothing"* and
merely *uses* `c_γ` in passing, in the ratio `λ_γ = c_γ/c_s`. The registry also already carries the
derived `λ_γ` and `h0`, which likewise no step reached.

⇒ **The seed carries quantities the walkthrough has not walked to.** ⭐ That is precisely the
**"universe hole"** the plan's §7a check 4 (top-down reconciliation) exists to catch — *a parameter that
no step ever named* — and it has surfaced at step 2, on a registry of eleven quantities.

⛔ **Do not fix it by back-filling `c_γ` into step 1.** `c_γ² = μ_R/ρ_br` is the light-sector cone
(`docs/model_map.md:63`) and belongs to the excitations phase; importing it into the medium's defining
properties would put a phase-2 input into the tier-1 core for bookkeeping convenience — the same error
step 1 avoided by refusing to record `c_s` there.

✅ **DECIDED 2026-07-31 — S0.5 PENDING EXECUTION.** The seed is **trimmed**: `c_γ` and `λγ` are retired
from the medium contract and re-introduced with provenance at v3 **S9** / **S20a**. ⛔ The "trim or keep"
question below is closed.

<details><summary>original wording</summary>
**What is open.** Whether the seed should be **trimmed back** to what the walkthrough has actually
reached, or **kept ahead** of it as a scaffold with the gap declared. ⚠ Either way the gap must be
visible, because `show_reduced.py`'s "5 simulation inputs" currently counts `c_γ` among them on no
step's authority.
</details>

### O-02 — steps 1 and 2 disagree on a count: is `K` + the exponent one entry or two?

**The finding.** Step 1 counts them as **two** separate "what's new" entries: `K` as a scalar primitive
(`00_medium_and_brane.md`, Scalars table) and the **EOS polytropic exponent** as a discrete/structural
choice (same file, Discrete table). Step 2 derives `[K] = M L^(4n−2) T⁻²` and concludes they are ⛔ **one
structural choice, not two independent inputs** (`01_sound_speed.md`, "the finding") — changing `n`
changes the *dimension* of `K`.

**Did step 1's headline change?** ⇒ **No — the headline stands at 4 scalars, 3 discrete**, and the entry
in `00_medium_and_brane.md` now says why. ⚠ **This is a judgement, made on §1.0's definition of the
count** (*"how many variables must I define for the simulation to run"*): a simulation must be given a
**numeric value for `K`** *and* a **selection of `n`** — supplying one does not supply the other, so both
remain members of the sim-input set, in different partitions. What step 2 establishes is that they are
⛔ **not independent** — the pair is constrained, not merged.

⛔ **This is the case §1.0 warns about**, and the reason it forbids reporting one integer: under the
*algebraic* residual-dimension reading the coupling plausibly removes a degree of freedom, while under
the *sim-input* reading it removes nothing. ⭐ Per §1.0, **that disagreement is itself the finding** —
⛔ do not reconcile it silently. **Open for the user:** whether the count should follow the sim-input
reading (headline unchanged, as recorded) or mark the pair as one coupled entry.
