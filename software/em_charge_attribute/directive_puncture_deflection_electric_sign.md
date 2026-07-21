# Directive — the puncture-deflection mechanism for the electric SIGN

**Status:** foundational derivation. **Target-blind.** Full gauntlet (Codex→Grok→Codex design-review BEFORE any build). Method = audit §F: **define the piece the model's OWN way, compute, ACCEPT whatever falls out — do NOT import Maxwell, do NOT chase Coulomb, do NOT assume the sign or that the datum is pinned.**

**Why this exists.** The verified leftover-scalar build (`e2d84a88`) landed `NO_NATIVE_CLAMP` + a *calibrated-knob* sign, with two soft spots: (a) nothing in G0 *holds* a signed `u_L`/deflection datum, and (b) the sign rode a *free* clamp-to-`g_χh` ratio. This directive tests a physically-motivated candidate for BOTH: **a `±w` puncture geometrically bends the brane into `±w`, and that bend is held by the puncture's existence (topological `±w`), with a magnitude set by the throat geometry (`~r_e`) — not a free postulate.** If real, it supplies the "holder" (resolving `NO_NATIVE_CLAMP`) and fixes the ratio (turning the sign from calibrated to *earned*). It is the physical argument for the boundary operator `𝔅` that U2 proved the bare math does not pick.

**Process (banked).** `codex exec -c model_reasoning_effort=xhigh … < /dev/null`; Mathematica-running builds add `--sandbox danger-full-access`; background/`setsid` for long runs (survives the harness kill; [[project-em-sector-reconsideration]]); never kill mid-`WolframKernel` run (seat leak). Dual-engine (SymPy + independent Mathematica) wherever analytic; per-tooth ablation; dimensional check units-restored; NO fake/commentary scripts. **The BUILD session (when it comes) must produce deliverables and EXIT — it does NOT run its own verification legs ([[feedback-build-session-stops-at-deliverables]]).**

---

## 0. Read-first (build from the committed structure — do not paraphrase this directive)

- The verified build + its verdict: `software/em_charge_attribute/leftover_scalar_electric_sign_{result,VERIFICATION}.md` (the pin→repel / source→attract fork; the free-knob finding this directive tries to resolve; the coupled `(u_L,h)` action with `B_eff=ρ_B0²/χ_c`, `K_h`, `C_hu`, `D=7/4`).
- The committed action + the zeroed coupling this mechanism turns ON: `software/em_charge_attribute/g0_closure_card_v0.md` — the declared-zero `r_B·u_L=0` (and any sleeve/puncture↔deflection coupling) is `[POSTULATE]`; the mechanism claims that zero is physically wrong.
- The physical picture (for the MECHANISM only — the audit/handoff remain authoritative for program state): `docs/conceptual_foundation.md` §4 (the throat = a puncture into `±w`; the brane *bends into `w`* at a throat; `r_e` = throat-body size from `ke²/a=m_ec²`; "like resist 4D overlap → repel, opposite merge → attract/annihilate"; charge = universal/quantized) + `software/stage1_solver/reports/pathA_42_charge_coupled_scalar.md` (the charge scalar `h` = **the electric field's flex INTO `w`** — the deflection field is already in the model).
- The maintained map: `docs/model_definition_audit.md` §B/§F (I2 sign IMPORTED; `r_e`), `docs/em_analog_next_phase_handoff.md`, `docs/two_throat_simulation_handoff_spec.md` (R1).

## 1. The object under test (state it target-blind — do NOT assume pinned or the sign)

**Mechanism:** a throat is a puncture through the brane into `±w`. Forming/holding that puncture **deflects the surrounding brane into the same `±w`** (an extrinsic bend / "flex into `w`", the field pathA_42 calls `h`, coupled to the density mode `u_L` via `C_hu`). The claim to TEST, not assume: this signed deflection is (i) **held** (nonrelaxable while the puncture exists — tied to the topological `±w`), (ii) of a **magnitude set by the throat geometry** (`~r_e`), not a free postulate, and (iii) sourced as a far-field **monopole**. Two same-`w` deflections add → stored-energy rises on approach; opposite cancel → falls (and can heal/annihilate).

⚠ The deflection field identity is itself part of the question: it is most naturally the `h` "flex-into-`w`" (currently modeled as a *source* → attract) and/or the coupled `(u_L,h)` mode. The mechanism REINTERPRETS the throat's deflection BC from a relaxable SOURCE to a geometric PIN — that pin-vs-source distinction is the sign, and it must be DERIVED from the puncture geometry, not assumed.

## 2. Requirement — three target-blind questions

**Q-PIN (the crux — derive, don't assume).** From the puncture/throat structure (the extrinsic bend required to connect brane→bulk through `w`; the `±w` topology; the throat's finite geometry), derive whether the mouth deflection is a **held Dirichlet-type datum** (pinned: it cannot relax while the puncture exists) or a **relaxable source** (it can be screened/relaxed with the puncture intact). State what — geometrically/topologically — would hold it, and whether that holder is native to the puncture (as opposed to the G0 `[POSTULATE]` zero). Emit `PIN_FORCED(mechanism)` / `SOURCE_RELAXES` / `UNDETERMINED_ANALYTICALLY(named gap → R1)`.

**Q-FORCE.** Given the derived BC, compute the **far-field** (`R ≫ r_e`) two-body interaction for like (`s₁=s₂`) and opposite `w`: sign, falloff (`1/R^n`), and the coefficient — from the deflection's stored gradient energy in the coupled `(u_L,h)` sector (with `C_hu`, per the verified build's machinery). Target-blind: compare to Coulomb (like-repel, `1/R²`, opposite-attract) ONLY in §4.

**Q-MAG.** Is the deflection magnitude **set by the throat geometry** (the `r_e` / self-energy balance) so that the coupling strength is a *prediction*, not a free knob? Characterize the dependence; if it needs the nonlinear throat solve, say so and defer to R1 (do not fake a number). (This is what would upgrade the sign from calibrated to earned.)

## 3. Method freedom + honest tiering (Codex designs the route)

Codex designs the computation. Honest tiering (state which tier each result is in):
- **Tier-A (analytic, do now):** the far-field force given each candidate BC (pin/source) — largely in hand from the build; the *geometric* pin-vs-source argument as far as a reduced/model throat embedding allows (e.g. minimal-surface / catenoid-like throat far-field bend); the monopole/quantization structure; the density-dependence functional form (Q from §5); consistency (Tier-A of the build: `D>0`, decoupling, zero-ledger-with-`r_B·u_L`-now-ON, charge⊥mass).
- **Tier-B (R1-deferred, name it — do NOT fake):** the full nonlinear throat solve that would DECISIVELY fix pin-vs-source and the `~r_e` magnitude; the close-range regime (§5).
No clamping to Maxwell; derive the sector's own equations. The pin-vs-source and magnitude may honestly be `R1_REQUIRED`.

## 4. Target-blind verdict block (EM target appears ONLY here) — TOTAL precedence

Compare to real-EM (like charges repel, `1/R²`) only here. First match wins:
1. Tier-A inconsistency (turning on `r_B·u_L`/the bend breaks a committed sector) → **`NO_GO(sector)`**.
2. **`PIN_FORCED`** (the puncture geometrically holds the signed deflection) — judge Q-FORCE:
   - far field = like-repel `1/R²` (opposite-attract) AND Tier-A consistent AND (Q-MAG: magnitude set by geometry) → **`SIGN_EARNED`** (the win — sign predicted, not calibrated; triple-check; note the `r_e` link).
   - like-repel `1/R²` but magnitude/full-consistency Tier-B-deferred → **`SIGN_FORCED_MAGNITUDE_R1_REQUIRED`** (sign forced by the pin; magnitude/consistency pending R1).
   - wrong-sign / wrong-range → **`NO_GO(pin_wrong_sign|range)`** (mechanism falsified).
3. **`SOURCE_RELAXES`** → the deflection is a relaxable source → far field attracts → **`MECHANISM_FALSIFIED_SOURCE_ATTRACT`** (honest no; the puncture does not hold the bend).
4. **`UNDETERMINED_ANALYTICALLY(named gap)`** → `R1_REQUIRED`: state exactly what the throat solve must decide (pin-vs-source and/or magnitude), first-class.
No code path branches on the desired answer before this block.

## 5. Scope + prediction-hooks (preserve — these are the payoff)

- **FAR-FIELD ONLY.** `R ≫ r_e`. **Close-range (`R ~ r_e`) is OUT OF SCOPE here but is the high-value follow-on:** direct 4D-body overlap → a finite-size repulsive core, deviations from Coulomb/QED at short range (held-out vs `e⁻e⁻`/`e⁻e⁺` scattering), and opposite-`w` annihilation dynamics. Flag it; do not compute it here.
- **DENSITY-DEPENDENCE (characterize).** `B_eff=ρ_B0²/χ_c` is density-dependent → the effective electric coupling varies with local medium density; since mass depletes density (gravity=drain), this is a **charge↔density↔gravity coupling** (Maxwell-forbidden). Derive its functional form and flag it as a held-out prediction (EM near strong gravity / black holes / neutron stars). Magnitude sim-gated.
- **MONOPOLE + QUANTIZATION.** Report whether the mechanism naturally yields a single-signed radial **monopole** and a **universal/quantized** charge (one puncture = one throat-set bend). A byproduct win if it does.

## 6. Acceptance (able-to-fail, target-blind, per-tooth)

- Faithfulness: the coupled `(u_L,h)` action + the newly-ON `r_B·u_L`/deflection coupling + `r_e`/throat picture transcribed from the committed sources term-by-term.
- Q-PIN is a DERIVED result (geometric/topological), able-to-fail: a control where the deflection is manifestly relaxable must return `SOURCE_RELAXES`.
- Q-FORCE sign/falloff extracted symbolically, `s₁,s₂` free, target-blind; the pin-BC and source-BC both computed (both signs reachable — but which the geometry *forces* is Q-PIN's job, not assumed).
- Q-MAG: the `r_e`/geometry dependence is computed or honestly `R1_REQUIRED` (no fabricated magnitude).
- Density-dependence functional form asserted + able-to-fail; monopole/quantization reported.
- Verdict precedence TOTAL; truth-table over all §4 landings; target-blindness guard (no hardcoded sign); mutation self-test (each tooth dies at its own assert).
- Dual-engine agree term-by-term where analytic; both exit 0.

## 7. Deliverables
- `puncture_deflection_electric_sign_check.py` (SymPy) + `.wl` (independent Mathematica), both exit 0, emitting Q-PIN / Q-FORCE / Q-MAG / §5 hooks + the §4 landing + all teeth.
- `puncture_deflection_electric_sign_result.md` — ≤2 pages: the derived pin-vs-source verdict + reasoning, the far-field force, the magnitude status (earned vs R1), the density/monopole hooks, the §4 landing, and an explicit DECIDED-vs-R1 split. No softening.
- **After build (orchestrator arranges, NOT the build session):** fresh-agent 4-leg + Grok compute-verify.
