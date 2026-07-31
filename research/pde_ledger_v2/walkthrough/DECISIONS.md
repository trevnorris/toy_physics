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

### Consequences, not yet applied

1. ⛔ **The pin `a` and the throat radius must be separated** — different quantities, different symbols.
   ⚠ Touches `stage004`, `parameter_register.md:132`, and the `a` rows in 016/018/023.
2. **The pin is probably unnecessary.** Nothing consumes `a` in the medium block — it is the output of its
   own pin relation and appears in no `input_qids`. It is also `ξ_h/√2` (`stage004:133`), i.e. the healing
   scale up to a constant, and `ξ_h` is the one with a physical justification ("core balance") while `a` is
   explicitly a pin choice.
3. ⚠ **`R2.a_pin`'s registry class is unresolved and the ambient count depends on it.** It was changed
   `CONVENTIONAL → DERIVED-EXECUTED` on a too-literal reading of the classification rule; a units-pin
   identity is ⛔ not "a defining equation in terms of other model quantities". If it reverts, the ambient
   count of 10 and the acceptance MATCH both need re-establishing.
4. ⭐ **The classification rule gains a clause:** *a relation arising from imposing unit pins is not a
   defining equation*, however much it looks like one.

⚠ **`ħ` stays open, and this is why.** `ħ` and `a` are joined only by that pin identity, so **neither is
derived from the other** — the relation is silent about which is fundamental. The hypothesis that `ħ` is a
particle size in the medium needs something *outside* that relation to pin one of them, and nothing
currently does. ⇒ Not a bookkeeping gap; a real gap in the model.
