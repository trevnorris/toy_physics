# Directive pathA_20 — Emergent constants, Step 2: derive `c_s` and `c` (the velocity structure)

**Status:** ✅ FINALIZED against the `pathA_19` foundation (landed + reviewed FIDELITY-CLEAN 2026-06-19; decision-13 §8).
Consumes the pathA_19 outputs verbatim: base dimensions = **`{L,T,M}`** (mass NOT shown emergent — the `m_defect`↔inflow
bridge is DEFERRED to `pathA_21`, not refuted); flux invariant **`J` with `[J]=T⁻¹`** (frame-independent, bulk+brane);
healing relations **`a=ħ/(m c_s0)`**, **`ξ_h=√2 ħ/(m c_s0)`**, **`h0=(5K/4)ρ0⁴=m c_s0²/4`**; the **`a→J` re-pin**
(`a`=branch-dependent mouth-moment, `J`=the invariant); the three R_norm-not-dimensionless conversion factors stay
DEFERRED to `pathA_22`. Ready for the Codex DIRECTIVE design-review (standing process for foundational directives), then
execution — gated by the user.
**Date:** 2026-06-19
**Owner:** Codex (DERIVES + consolidates + codes; iterates until scripts exit 0). Claude reviews afterward.
**Trigger:** decision-13 §4 items 1–2, Step 2 of 4. Chain: `pathA_19` (foundation) → **this (`c_s`+`c`)** →
`pathA_21` (`G`) → `pathA_22` (scale-map → `m̂0²·S_port` → re-run B2c).

## Why this step (context)
With the dimensional foundation fixed by `pathA_19`, derive the two propagation/limit speeds. Honest priors from the
pre-work: (a) `c_s` is stated everywhere as `c_s²=5Kρ⁴/m` but only relative to the IMPOSED polytropic EOS; (b) the
existing material (`research/em_fields/paper/em_fields.tex` ~160/168/178/699) derives `c` as the acoustic-metric
light-cone speed and identifies `c=c_s` in the weak-field limit — so the honest prior is `c=c_s`, and the
terminal-velocity / standing-wave reading is a hypothesis to TEST, with the math deciding.

## Scope & stance
DERIVATION + CONSOLIDATION. No model-formula changes, no freeze touch, no `m̂0²·S_port` un-pin. Extend
`dimensional_check.py` (new group; leave pathA_18 + the pathA_19 foundation group intact). Same infra constraints as
`pathA_19` (read-only sim dir; no `$RT exec`; `timeout 600`; ≤2 MMA seats; YAML/md). **Dual-engine required** where MMA
can verify the `c_s`/`c` algebra and the dimensional reductions. Use the base system + flux invariant from `pathA_19`.

## Work items

### S1 — `c_s` (sound speed): consolidate + make rigorous + dimension-check
- Consolidate `c_s²(ρ)=(1/m)(dP/dρ)=5Kρ⁴/m` from the EOS, in the `pathA_19` base dimensions; verify `[c_s]=L·T⁻¹`
  symbolically (SymPy + `.wl`).
- State-dependence: confirm `c_s(ρ)∝ρ²` for `n=5` and record `∂ln c_s/∂ln ρ`. `ρ` is NOT a free dial — it is set by
  the coupled stationary flow (continuity + Bernoulli, see S2), so `c_s` drifts because the flow makes `ρ` drift.
  Treat `c_s` as a profile `c_s(x)` with asymptotic constant `c_{s,0}`.
- **Honesty requirement:** the EOS `P=Kρ⁵` is an IMPOSED polytropic closure, so `c_s` is "derived" only relative to it.
  State this provenance plainly; deriving the EOS from the GNLS action is NOT required here but the gap must be flagged.

### S2 — the three velocity scales; derive `c` as the photon/standing-wave ceiling
This model has (at least) THREE distinct velocities, all `L·T⁻¹` (so dimensional analysis cannot separate them — only
the dynamics can). Name, dimension-pin, and RELATE them — do not collapse them prematurely:
  - `v_b` — **background condensate flow velocity** `v_b=(ħ/m)∇θ` (drift of the medium). GRAVITATIONAL-sector
    variable: gradients of `v_b` and `ρ` are the analog field; `v_b=c_s` is the acoustic horizon. Varies in space.
    **`ρ` is ALSO coupled to `v_b`:** stationary continuity `∇·(ρv_b)=0` + quantum-Bernoulli `½mv_b²+μ(ρ)+V_conf+Q=const`
    (`μ=h=(5K/4)ρ⁴`) lock `ρ, v_b, c_s(ρ)` into ONE self-consistent profile (flow up → `ρ` down → `c_s` down). LOCATE
    the existing `ρ(v_b)`/Bernoulli equations in the ledger and incorporate (do not reinvent). Name + dimension-pin
    `v_b`; its full gravitational role is `pathA_21`.
  - `c_s` — **bulk sound speed** (density/phonon waves), from S1.
  - `c_γ` — **photon/gauge-wave speed** (massless gauge excitation on the brane); the brane light cone.
- **`c` as the terminal-velocity ceiling, from the standing-wave ontology.** A massive particle (throat) is a
  **standing wave of the gauge/photon field** — two counter-propagating `c_γ`-waves. Derive: (a) rest internal
  oscillation = Compton frequency `ω₀=m*c_γ²/ħ`, so `E_rest=m*c_γ²` is trapped-wave energy; (b) driving the envelope at
  `v` Doppler-shifts the components and slows the internal clock as `ω₀/γ` (time dilation), freezing (`→0`) as
  `v→c_γ`; (c) the envelope cannot outrun its constituent `c_γ`-waves, so the terminal ceiling is `c=c_γ`,
  `lim_{E→∞}∂E/∂p=c_γ`. This is the model's physical origin of relativistic kinematics and of `E=m*c²`. The
  drag/enforcement view (Cherenkov/Mach, Landau `v_c=min_p ε(p)/p`) agrees in the infinite-drive limit.
- **The mass-bridge that `c` unlocks (record the FORM here; the derivation + M-collapse is `pathA_21`).** The
  `pathA_19` honest negative kept `{L,T,M}` ONLY because eliminating `M` requires `c`, which did not exist yet
  (circular). The natural bridge is **`m_defect = α_J · ħ J / c²`** — i.e. `E=m c²` with the INFLOW setting the energy:
  the drain rate `J` (a `T⁻¹` frequency) gives a rest energy `ħ·(α_J J)` and `/c²` returns a mass. This is the SAME
  statement as the standing-wave rest energy `E_rest=ħω₀` (so mass-as-inflow ≡ mass-as-standing-light: `ω₀=α_J J`).
  **de Broglie consistency:** in `ħ=c=1` units mass, rest frequency, and `J` are all `T⁻¹` (dimensionally identical);
  in full `{L,T,M}` they differ by the `ħ/c²` bridge. So the `m↔J` dimensional equivalence is a DERIVED destination
  (granted once `c` here + the `ħ` verdict in S3 are in hand), NOT a startable axiom. pathA_20's job is to put the
  derived `c` and `c_γ` on the table so `pathA_21` can run the bridge + M-collapse test non-circularly; do NOT collapse
  `M` in this step. Tie `m*` to the `pathA_19` flux/inflow mass (`J`, F1/F2) symbolically only.
- **The surviving open number: is `c_γ=c_s`?** The standing-wave argument forces ceiling = photon speed `c_γ`; whether
  THAT equals the bulk sound speed depends on whether gauge + density share the acoustic metric, or the localization
  profile `Z(w)`/width `a` rescales the brane gauge cone (the EM sector already shows `μ₀^eff=μ₀/Z_int`). Derive
  `c_γ/c_s` in closed form, AND determine whether it is `ρ`-INDEPENDENT (a pure number) even though each speed drifts.
- **`(c/c_s)³`:** identify the tail factor `R_tail=Θ_tail(c/c_s)³−1` with the derived `c_γ/c_s` (=1 if they coincide);
  reconcile with the `c=c_s` weak-field statement in `em_fields.tex`.
- **Constants vs profiles:** separate genuine constants (`K,m,ħ` + the flux `J`/asymptotic `ρ₀,c_{s,0}`, per
  `pathA_19`) from profiles (`ρ(x),v_b(x),c_s(x)`, possibly `c_γ(x)`); state which reference `c_s=1` denotes.

### S2b — the throat as a transonic drain: the flux law `J(ρ₀, geometry)` (replaces the "constant `J`" assumption)
`pathA_19` flagged `NO_NET_ACCRETION_BC_UNDERIVED` and did NOT prove the flux is a universal constant. Resolve the
THREE distinct "is it constant?" questions, from the coupled continuity+Bernoulli+EOS profile (S1/S2), and TEST the
transonic hypothesis rather than assuming it:
  - **Flux conservation ≠ velocity constancy.** In steady state `∇·(ρv_b)=0` makes `J=∮ρv_b·dΣ` the SAME through any
    surface enclosing the drain (Gauss, no sources/leakage) — but `v_b` is NOT constant: as the throat converges
    (`ρ` drops, cross-section narrows) `v_b` ACCELERATES to keep `J` fixed. The throat is a **nozzle**. State both
    explicitly; do not conflate the conserved integral `J` with the local speed `v_b`.
  - **Sonic point = acoustic horizon.** Locate where `v_b` crosses `c_s(ρ)` along the profile. A subsonic→supersonic
    (transonic) drain is the canonical analog black hole (Unruh draining-bathtub); `v_b=c_s` is the horizon (already
    named in S2). Determine whether the sonic point sits at/near the mouth.
  - **Choked flux ⇒ `J` is set by the background, not universal.** If the throat is transonic the flux is CHOKED —
    pinned by the critical (sonic) condition at the throat from the upstream stagnation state (the background
    condensate). Derive the critical-flux scaling and DIMENSION-CHECK it to `T⁻¹`, e.g. `J_crit ~ ρ₀ c_{s,0} A_throat`
    (`ρ₀ c_{s,0} a³` bulk → `T⁻¹`; `ρ_{3,0} c_{s,0} a²` brane → `T⁻¹`), and express its dependence on the background
    pressure `P₀=Kρ₀⁵`. This makes `J` (hence any inflow-mass) **environment-dependent on `ρ₀`** — a derived feature,
    not a frozen input. Confirm/deny `J_in=J_out` (true no-net-accretion) is still gated by the throat-bottom topology
    (open/closed/connected, undetermined) — carry it forward if unresolved.
  - **Honesty requirement:** "transonic choked drain" is a CANDIDATE mechanism. Derive the stationary profile and TEST
    whether the throat is actually transonic; a non-transonic (everywhere-subsonic) result is a valid outcome that
    changes the flux law. Do not assert the horizon/choking; let the profile decide.

### S3 — `ħ` provenance: fundamental input or emergent? (named fork; feeds the `pathA_21` M-collapse test)
`pathA_19` F3 derived `a=ħ/(m c_s0)`, i.e. **`ħ = (dimensionless factor)·m_GNLS·c_s0·a`** — so among `{ħ, m, c_s0, a}`
only THREE are dimensionally independent. Calling `ħ` fundamental and `a` derived (as F3 did) is a CONVENTION, not a
physics result. Investigate what `ħ` actually IS in this model and whether it can be expressed from the other genuine
constants (`K, ρ₀, m`) or a substrate length:
  - Catalog `ħ`'s structural roles and dimension-check each: quantum of circulation (`κ=2πħ/m`); the
    phase↔momentum exchange rate (`p=ħ∇θ`, with `θ` dimensionless — "momentum carried per radian of phase twist");
    the quantum-pressure / core-size scale (`Q=-ħ²/2m ∇²√ρ/√ρ`, sets the healing length `ξ_h=√2 ħ/(m c_s0)`).
  - **The fork (state which, or carry as a named residual):** if `ħ` is EMERGENT (a derived combination, or fixed by a
    microscopic substrate length so `a` is the fundamental scale and `ħ=m c_s0 a/√2`), then `ħ=1` is a derived relation
    and `M` can collapse cleanly in `pathA_21`. If `ħ` is a genuine FUNDAMENTAL input (irreducible action quantum of
    the medium), then `ħ=1` stays a unit choice and the eventual `M`-collapse rests on that fundamental. A "`ħ` stays
    fundamental, here is why" verdict with its blocking reason is a VALID result — do not force emergence.
  - Record the verdict + reasoning; it is a direct input to the `pathA_21` mass-bridge / M-collapse test. No base-set
    change is MADE in `pathA_20` (that is a `pathA_19`-class base decision); S3 only DETERMINES the status and hands it
    forward.

## Acceptance criteria (to finalize against the pathA_19 base)
1. `c_s` consolidated + machine-verified `[c_s]=L/T` in the established base; state-dependence + Bernoulli coupling
   recorded; EOS-closure provenance stated honestly.
2. Three velocities named + dimension-pinned; `ρ`–`v_b`–`c_s` coupling explicit (ledger equations located); constants
   vs profiles separated.
3. `c` derived as the standing-wave photon ceiling `c=c_γ` (Compton freeze, `m*`, `E=m*c²`, `m*` tied to the
   `pathA_19` inflow-mass); `c_γ/c_s` in closed form + its `ρ`-(in)dependence; `(c/c_s)³` identified and reconciled
   with `em_fields.tex`. The mass-bridge FORM `m_defect=α_J ħ J/c²` + the de Broglie `ħ=c=1` equivalence note are
   RECORDED (not executed); `M` is NOT collapsed here.
4. (S2b) the transonic hypothesis is TESTED on the derived profile; flux-conservation-vs-velocity-acceleration stated;
   if transonic, the choked flux law `J_crit(ρ₀, geometry)` is derived + machine-checked to `[J]=T⁻¹` + its `P₀=Kρ₀⁵`
   dependence recorded; the sonic-point/horizon location reported; no-net-accretion (`J_in=J_out`) status carried.
5. (S3) the `ħ`-provenance verdict (fundamental vs emergent) is stated with reasoning + dimension-checked role catalog;
   no base-set change is MADE; the verdict is handed forward as a `pathA_21` input.
6. Dual-engine `.wl` agreement; new harness group passes; scripts exit 0 within `timeout 600`; pathA_18 +
   pathA_19 groups untouched.

**Acceptance is PASS/FAIL with NAMED RESIDUALS (per `pathA_19`): script exit-0 is NECESSARY, NOT SUFFICIENT.** Every
rejected hypothesis (e.g. non-transonic throat, `c_γ≠c_s`, `ħ` stays fundamental) must leave a named residual + source
+ downstream consequence. Negative/`UNDETERMINED` results are valid outcomes, not execution failures.

## Out of scope
The mass-bridge DERIVATION + M-collapse test (`pathA_21`, using the `c` + `ħ`-verdict from here); emergent `G` and the
`m↔G` unification (`pathA_21`); scale-map → `m̂0²·S_port` (`pathA_22`); B2c rerun; freeze amendment. Do NOT collapse `M`
or change the base set in this step.

## Review (orchestrator, after Codex)
Transliteration-fidelity (one clean agent per new script); independent dimensional re-derivation of `[c_s]`,`[c]`,
the three velocities; adversarial pass (distrust-all-clean) targeting the `c=c_s` vs `c≠c_s` verdict and the
EOS-closure honesty. Claude reads only residuals. Then gate to `pathA_21`.
