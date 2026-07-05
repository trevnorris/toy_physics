# FOLLOW-UP (post-pathA_39 Stage 4): investigate the charge-coupled scalar (the "extra scalar" departure)

> **✅ DONE (2026-07-05) — `pathA_42` = `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED`.** The scalar does NOT clearly break the model, nor
> is it trivially hidden: it's a MAPPED departure with the decisive magnitude SIM_GATED (walled on the deferred throat solve —
> reachability read = `HARD_WALL`; `M_h`/`c_E` un-pinnable now). 5 channels: `h_EP`=EP-safe-on-decoupled-floor against Gate-L's
> mass-coupled fifth force BUT only CONDITIONAL on `q_h/Q_E` throat universality (GLM F1 — the "absorb into Q_E" step needs universality);
> `radiation`=`SIM_GATED` exponent-3 on the current ledger (Codex adjudicated against pathA_38: the static Coulomb calibrates via `Q_E`+`N₀`
> NOT the dynamic stiffness `K_h`, so `M_h` stays free → radiation gated even at `c_E=c_γ`, NO order-one Lorentz-vs-radiation conflict — the
> GLM-F6 exponent-1 tension is CONDITIONAL on a future `K_h`-pinning); `universality`/`u_L_EP`(+mixing)/`preferred_frame`=`SIM_GATED` variants;
> vacuum-Cherenkov (`c_E<c_γ`) + stellar-cooling/cosmological bounds = throat-solve re-open obligations. Full gauntlet (Codex design-review
> 4 BLOCKERs + GLM-5.2 tertiary 8 findings + Codex normalization-adjudication, all folded); tri-review CAUGHT 3 stamped/theatrical
> mechanisms (guard-A theater, fabricated deletion-robustness, stamped exponent) → remediated (guard = genuinely WIRED firewall reading
> substrate provenance + scanning the results path; genuine deletion recompute; mechanized dipole/flux exponent) → re-tri-review CLEAN +
> final tuple-scan hardening. `ENGINE_AGREE 22`, deterministic. Gate pins nothing → NO NG5/pathA_40 reopen. Artifacts:
> `{directives,tools,reports}/pathA_42_charge_coupled_scalar*`. The investigation below is the ORIGINAL framing (kept as history).

**Trigger (user, 2026-07-04):** the magnetism sector lands as a *characterized departure* — magnetism exists as a real force
(pathA_39 Stage 2, `MAGNETIC_FORCE_DERIVED`, tri-reviewed) but is inseparable from an extra **charge-coupled scalar**. The user
asked: what is it, why don't we see it in everyday life, could it exist but only show up very specifically? Investigate AFTER
Stage 4 (do not derail the sector close).

## What the extra scalar IS (computed — pathA_39 Stage 0+1 + Stage 2, both tri-reviewed CLEAN)
A propagating scalar mode of the medium that the electric charge couples to. Two pieces:
- **`h`-branon** — the propagating version of pathA_38's STATIC Coulomb mediator `h` (brane flexing into `w`), speed `c_E`, charge
  residue `∝ q_h²/M_h`. Import-forced by pathA_38's `M_h>0` (the static Coulomb mediator propagates dynamically).
- **Density/`c_s` mode** — the compression wave = the gravity/density channel, speed `c_s`, charge residue `∝ q_L²` (rests on the
  sim-deferred moving-source amplitude `a_L`).
In Stage 2 it reappears as the **unavoidable attractive scalar-current admixture** (the longitudinal `u_L` channel; both channels
attractive for like currents, no cancellation).

## What we DON'T know (the investigation)
Its observable MAGNITUDE — which hinges on quantities currently sim-deferred/uncomputed: `M_h`/`c_E` (`h`-branon dynamism), `a_L`
(density current strength). The structure is computed; the size is not.

## Why it may be hidden (hypotheses to test)
1. **Nearly non-radiating (the clean one):** as `M_h→0`, `c_E→∞`, `h` becomes NON-dynamical = an instantaneous Coulomb constraint,
   exactly like real EM's scalar potential (which does not radiate). **The departure VANISHES in the stiff limit → the model reduces
   to ordinary EM.** Departure size = how dynamical `h` is (`M_h`/`c_E`, unpinned). If `c_E` large ⇒ scalar barely radiates ⇒
   invisible except at extreme accelerations/frequencies (the user's "shows up very specifically").
2. **Gravity-weak:** the density/`c_s` piece is the gravity/density channel ⇒ charge coupling to it radiates at gravitational
   strength (`~10⁻³⁶` of EM) ⇒ unobservable by construction.
3. **Fifth-force constrained:** a long-range charge-coupled scalar is tightly bounded by Eötvös/fifth-force experiments. Cuts both
   ways: the model MUST predict it weak to be consistent; an appreciable predicted magnitude = a genuine tension/falsification.

## The investigation to run (post-Stage-4)
- Pin down `M_h`/`c_E` (is the `h`-branon nearly non-radiating? what sets `c_E`?) and `a_L` (density-current strength) — some may
  need the sim-deferred moving-throat solve; identify what is algebra-reachable now.
- Compute the scalar's observable magnitude + its coupling to real observables (scalar radiation from accelerating charges; an
  extra static scalar force between charges beyond Coulomb).
- Map to experimental constraints (fifth-force / Eötvös / scalar-radiation bounds; equivalence-principle tests).
- Decide the verdict: does the model predict it's NATURALLY hidden (consistency + a possible new prediction at extreme
  accelerations/high frequencies) — or an appreciable, already-excluded departure (falsification)?
- Falsification framing ([[feedback-falsification-is-the-goal]]): this is first-class predictive surplus — a specific prediction
  that either consistently hides, newly predicts, or falsifiably conflicts. Report straight either way.

## Pointers
`reports/pathA_39_{scalar_admixture_screen,magnetic_force}.md`; `directives/pathA_39_magnetism_close_maxwell.md`;
`docs/conceptual_foundation.md` ⭐ v8 block; memories `[[project-brane-existence-defect-structure]]`,
`[[feedback-native-em-mechanisms]]`, `[[project-calibrated-pde-goal]]`.

---

## SHARPENED (2026-07-05, post-orientation) — corrected framing + the user-chosen scope

**Correction that reframes everything: `M_h` is the branon INERTIA (coeff of `(∂_t h)²`), NOT a mass gap.** There is no
`μ²h²` term ⇒ dispersion `ω²=c_E²k²` ⇒ **`h` is massless/gapless, long-range** (both scalar channels are; not screened).
So this is a genuine EXTRA physical DOF (Stage-4 counted 4 EM DOF vs Maxwell's 2; deleting the `h` block → exact Maxwell).

**The clean picture of the departure:**
- STATICALLY the model is SAFE — the static Coulomb 1/r² IS `h` (pathA_38, calibrated via `Q_E`); no spurious static fifth
  force. `h` couples to CHARGE only (mass residue exactly 0).
- The departure is DYNAMICAL: real EM's scalar potential `A₀` is a non-radiating gauge CONSTRAINT, but here the Coulomb
  mediator `h` is a genuine propagating DOF ⇒ **an accelerating charge radiates scalar `h`-waves** (real EM cannot). Size
  set by how dynamical `h` is = by `c_E` (and `M_h`).
- ⭐ SHARP TRADE-OFF: `c_E→∞` ⇒ `h` instantaneous/non-radiating (like `A₀`) BUT preferred-frame (Lorentz-violating);
  `c_E=c_γ` ⇒ Lorentz-invariant BUT `h` is a real radiating extra scalar. Real EM gets BOTH only because `A₀` is gauge, not
  physical. Here `h` is physical ⇒ you get Lorentz invariance OR a hidden scalar, not trivially both.

**Two departure channels, cleanly separated:** (1) `h`-branon — charge-only; observable = scalar radiation + preferred-frame
unless `c_E=c_γ`; residue floor `4Q_E²tanh²(b/ℓ)/(M_h b²)` is `a_L`- AND `C_hu`-INDEPENDENT ⇒ robust, not tunable away.
(2) `u_L` density mode — couples to charge AND MASS ⇒ the fifth-force/EP (Eötvös) channel, at gravitational/density strength;
magnitude rides on the sim-deferred `a_L`.

**Computable now:** static=Coulomb; the radiation/Lorentz trade-off; the analytic stiff-limit control (never run); the
`q_h/Q_E=2tanh(b/ℓ)/b` universality test (species-dependent `b/ℓ` ⇒ charge/α non-universality). **Sim-deferred:** density
magnitude (`a_L`); numeric `M_h,c_E` — BUT these are ANALYTIC branon-action targets (project the bulk action onto the known
`h` zero-mode: `M_h`=∫inertia·profile², `K_h`=∫stiffness·(∂profile)²), NOT the 4D sim.

**USER SCOPE DECISION (2026-07-05):** "algebra-reachable characterization PLUS attempt to pin `M_h`/`c_E` via the branon
action." ⭐ The payoff: computing `c_E²=K_h/M_h` from the SAME brane that gives `c_γ²=μ_R/ρ_br` may FORCE `c_E=c_γ` (or
confirm a genuine gap) — resolving the pathA_40 cone-lock FROM BELOW and settling the trade-off. NOTE: this EXECUTES part of
the NG5 pending Route-A ⇒ success triggers the NG5 forward re-open (flips `c_E`, maybe `ρ_br,μ_R`, SIM_DEFERRED→DERIVED).
Scoping with Codex first (`_scratch/pathA_42_scalar_scope_codex.md`), then directive → gauntlet. Prospective gate = `pathA_42`.
