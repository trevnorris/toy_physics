# Decision 16 — Retire the brane polar field `P` (the "little arrows") from the frozen action

**Decided by the user, 2026-07-13**, on the U1 Phase-A finding. Status: OPERATIVE for all EM-analog/U1+ work; the formal ledger-v2 amendment (stage007 T0 polar terms) is banked as PENDING until the paused ledger resumes.

## The decision
The brane polar orientation field `P` (introduced as medium Requirement A3 / pathA_24 T0 freeze; assigned payoffs "light + charge") is **removed from the frozen brane action as an independent dynamical field**, together with its couplings (`λ_Pu` twist–shear, `αAniso` easy-plane wall term) and T0 polar coefficients. Rationale (user's framing): the brane's degree-of-freedom budget is already spoken for — inflow/compression = gravity, shear = light, mass-flow swirl = gravitomagnetism — and charge/magnetism live on the **bulk side of the throat** (`±w` body beyond the mouth; its motion = magnetism). Piling brane-side EM machinery on top was over-stuffing, and every computed payoff for `P` failed independently.

## Evidence backing it (all computed + verified)
1. **Charge:** native-`P` constraint gate → `NATIVE_P_NO_EMERGENT_GAUSS` (arrows cannot host a first-class U(1); charge = the throat's `±w` orientation instead).
2. **Light:** `pathA_35 FAIL_COUPLE_STRESS_NOGO` (light's stiffness cannot be derived from `P` substructure); `pathA_36` derives both photons with no `P` sector present.
3. **Wall:** the polar-vector wall self-localization was falsified (`pathA_24` T1); the wall is `χ_B`/smectic.
4. **Active harm:** U1 Phase A (verified, contested, `INSTABILITY_CONFIRMED_STRUCTURAL`) — the frozen massless-`P` + `λ_Pu` baseline has a long-wavelength helical (Lifshitz-type) instability for ANY `λ_Pu ≠ 0`: `det = k²(2A_P μ_R k − λ_Pu²) < 0` for `k < λ_Pu²/(2A_P μ_R)`; the pre-registered `pathA_35` G0.3 `FAIL_NOT_BOUNDED_BELOW` exposure firing. The stage006 anisotropy is easy-plane as written and cannot gap it; no gyroscopic rescue exists.

## Scope + door left open
- Retirement is of the **brane** field's *independent dynamics*. The still-speculative **bulk**-order hypothesis (supersolid lattice as a possible U2-constraint host) is a different object (bulk positional order, other side of the mouth) and is NOT foreclosed. If orientational order is ever *absolutely required*, it re-enters through a new T0-level freeze with its own gauntlet — not by un-retiring this.
- The A3 concept that survives is the head≠tail (`+w ≠ −w`) polarity itself — migrated to the throat sleeve orientation, where charge now lives.

## Consequences / follow-ups
- **U1 Phase A is re-run on the `P`-retired action** (computed, not assumed): expected 6 channels; verdict recomputed honestly. Report/results updated in `software/em_charge_attribute/`.
- **Parameter register:** `λ_Pu`, `αAniso`, and the T0 polar coefficients drop from the live knob set (a knob REDUCTION; register update + read-only Codex verify with the U1 rows).
- **Ledger v2 (paused):** stage007's T0 polar terms get a formal amendment banking this decision when the ledger resumes; no built v2 stage is touched now.
- Instability finding is retained as a first-class result ABOUT the retired baseline (commit `a5c079eb`): it is the reason the retirement is not arbitrary.
