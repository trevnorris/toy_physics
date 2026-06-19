# Directive pathA_20 — Emergent constants, Step 2: derive `c_s` and `c` (the velocity structure)

**Status:** ✅ FINALIZED against the `pathA_19` foundation (landed + reviewed FIDELITY-CLEAN 2026-06-19; decision-13 §8)
AND directive-design-reviewed (Codex `gpt-5.5` → SOUND-WITH-FIXES; all 10 fixes applied 2026-06-19; pending a confirm-
pass before execution). Consumes the pathA_19 outputs verbatim: base dimensions = **`{L,T,M}`** (mass NOT shown emergent
— the `m_defect`↔inflow bridge is DEFERRED to `pathA_21`, not refuted); flux invariant **`J` with `[J]=T⁻¹`** (frame-
independent, bulk+brane); healing relations **`a=ħ/(m_GNLS c_s0)` ⟹ `ħ=m_GNLS·c_s0·a`** and **`ξ_h=√2 ħ/(m_GNLS c_s0)`
⟹ `ħ=m_GNLS·c_s0·ξ_h/√2`** (NOTE `a` and `ξ_h` differ by `√2`: `a=ξ_h/√2` — do NOT mix the factor), **`h0=(5K/4)ρ0⁴=
m_GNLS c_s0²/4`**; the **`a→J` re-pin** (`a`=branch-dependent mouth-moment, `J`=the invariant); the three R_norm-not-
dimensionless conversion factors stay DEFERRED to `pathA_22`.
**Date:** 2026-06-19
**Owner:** Codex (DERIVES + consolidates + codes; iterates until scripts exit 0). Claude reviews afterward.
**Trigger:** decision-13 §4 items 1–2, Step 2 of 4. Chain: `pathA_19` (foundation) → **this (`c_s`+`c`)** →
`pathA_21` (`G` + mass-bridge) → `pathA_22` (scale-map → `m̂0²·S_port` → re-run B2c).

## Why this step (context)
With the dimensional foundation fixed by `pathA_19`, derive the two propagation/limit speeds. Honest priors from the
pre-work: (a) `c_s` is stated everywhere as `c_s²=5Kρ⁴/m_GNLS` but only relative to the IMPOSED polytropic EOS; (b) the
existing material (`research/em_fields/paper/em_fields.tex` ~160/168/178/699) derives `c` as the acoustic-metric
light-cone speed and identifies `c=c_s` in the weak-field limit — a USEFUL PRIOR but NOT sufficient proof (it lives in
the legacy 3D/SI acoustic-reuse frame flagged in `pathA_19`). So `c_γ=c_s` must be DERIVED from the wave-sector operator
or recorded UNDERDETERMINED; the terminal-velocity / standing-wave reading is a hypothesis to TEST, with the math
deciding.

## Scope & stance
DERIVATION + CONSOLIDATION. No model-formula changes, no freeze touch, no `m̂0²·S_port` un-pin, no `M`-collapse. Extend
`dimensional_check.py` (new group; leave pathA_18 + the pathA_19 foundation group intact). Same infra constraints as
`pathA_19` (read-only sim dir; no `$RT exec`; `timeout 600`; ≤2 MMA seats; YAML/md). **Mass-symbol discipline:** use
`m_GNLS` for ALL EOS/healing/Madelung/action/sound-speed formulas (the constituent mass) and `m_defect` for ALL
throat rest-energy/standing-wave/inflow-mass formulas; any equality between them is OUT OF SCOPE (`pathA_21`). **Dual-
engine** is REQUIRED for the ALGEBRAIC/DIMENSIONAL claims only (see acceptance §7); non-algebraic physics judgments get
human-readable residuals, NOT `.wl` certification. Use the base system + flux invariant from `pathA_19`.

## Work items

### S1 — `c_s` (sound speed): consolidate + make rigorous + dimension-check
- Consolidate `c_s²(ρ)=(1/m_GNLS)(dP/dρ)=5Kρ⁴/m_GNLS` from the EOS, in the `pathA_19` base dimensions; verify
  `[c_s]=L·T⁻¹` symbolically (SymPy + `.wl`).
- State-dependence: confirm `c_s(ρ)∝ρ²` for `n=5` and record `∂ln c_s/∂ln ρ`. `ρ` is NOT a free dial — it is set by
  the coupled stationary flow (continuity + Bernoulli, see S2), so `c_s` drifts because the flow makes `ρ` drift.
  Treat `c_s` as a profile `c_s(x)` with asymptotic constant `c_{s,0}`.
- **Honesty requirement:** the EOS `P=Kρ⁵` is an IMPOSED polytropic closure, so `c_s` is "derived" only relative to it.
  State this provenance plainly; deriving the EOS from the GNLS action is NOT required here but the gap must be flagged.

### S2 — the three velocity scales; derive `c` as the photon/standing-wave ceiling
This model has (at least) THREE distinct velocities, all `L·T⁻¹` (so dimensional analysis cannot separate them — only
the dynamics can). Name, dimension-pin, and RELATE them — do not collapse them prematurely:
  - `v_b` — **background condensate flow velocity** `v_b=(ħ/m_GNLS)∇θ` (drift of the medium). GRAVITATIONAL-sector
    variable: gradients of `v_b` and `ρ` are the analog field; `v_b=c_s` is the acoustic horizon. Varies in space.
    **`ρ` is ALSO coupled to `v_b`:** stationary continuity `∇·(ρv_b)=0` + quantum-Bernoulli
    `½m_GNLS v_b²+μ(ρ)+V_conf+Q=const` (`μ=h=(5K/4)ρ⁴`) lock `ρ, v_b, c_s(ρ)` into ONE self-consistent profile (flow up
    → `ρ` down → `c_s` down). LOCATE the existing `ρ(v_b)`/Bernoulli equations in the ledger and incorporate (do not
    reinvent). Name + dimension-pin `v_b`; its full gravitational role is `pathA_21`.
  - `c_s` — **bulk sound speed** (density/phonon waves), from S1.
  - `c_γ` — **photon/gauge-wave speed** (massless gauge excitation on the brane); the brane light cone.
- **`c` as the terminal-velocity ceiling — derive it WITHOUT assuming `E=m c²`.** A massive particle (throat) is a
  standing wave of the gauge/photon field (two counter-propagating `c_γ`-waves). Derive the ceiling from FIRST
  PRINCIPLES of the wave sector — the gauge-wave dispersion relation / group velocity, or a trapped-mode Hamiltonian
  for the bound standing mode — NOT from the Compton relation: (a) the terminal/group-velocity ceiling is the
  constituent wave speed `c_γ` (`lim_{E→∞}∂E/∂p=c_γ`); the envelope cannot outrun its constituent `c_γ`-waves;
  (b) for the bound mode, obtain the rest oscillation `ω₀` from the trapped-mode condition and DERIVE the dilation
  `ω₀→ω₀/γ` under drive. **The Compton relation `ω₀=m_defect c_γ²/ħ` and `E_rest=m_defect c²` are OUTPUTS to be checked
  for consistency, NOT premises — a derivation that STARTS from `E=m_defect c_γ²` is REJECTED.** The drag/enforcement
  view (Cherenkov/Mach, Landau `v_c=min_p ε(p)/p`) is a cross-check in the infinite-drive limit.
- **The mass-bridge `c` unlocks (RECORD the candidate form only; derivation + `α_J` + M-collapse = `pathA_21`).**
  pathA_19 kept `{L,T,M}` because eliminating `M` needs `c` (circular without it). The candidate bridge is
  **`m_defect = α_J · ħ J / c²`** — `E=mc²` with the INFLOW as the energy (drain rate `J` → rest energy → mass), the
  SAME statement as the standing-wave `E_rest=ħω₀` (`ω₀=α_J J`). **`h` vs `ħ` (the `2π`):** if `J` is a CYCLE-count rate
  (windings/drainage events per time, `ν`-like — see S3), it pairs with `h=2πħ`, so the bridge carries an explicit
  factor `m_defect = α_J h J_ν/c² = 2π α_J ħ J_ν/c²`; whether `α_J` absorbs the `2π` or it survives is a `pathA_21`
  derivation — FLAG it, do not resolve it here. **de Broglie note:** in `ħ=c=1` units mass, rest frequency, and `J` are
  all `T⁻¹`; in full `{L,T,M}` they differ by the `ħ/c²` bridge — so the `m↔J` equivalence is a **candidate conversion
  ONLY, physically granted iff `pathA_21` derives the source/boundary/Hamiltonian relation and `α_J`** (NOT granted by
  merely having `c`+`ħ` in hand). Do NOT collapse `M`. Tie `m_defect` (=`m*`) to the pathA_19 inflow mass (`J`, F1/F2)
  symbolically only; keep `m_GNLS` and `m_defect` distinct.
- **The open number: is `c_γ=c_s`?** The standing-wave argument forces ceiling = photon speed `c_γ`; whether THAT equals
  the bulk sound speed depends on whether gauge + density share the acoustic metric, or the localization profile
  `Z(w)`/width `a` rescales the brane gauge cone (the EM sector already shows `μ₀^eff=μ₀/Z_int`). **DERIVE `c_γ/c_s`
  from the gauge/density KINETIC OPERATOR (the wave-sector action)**, and determine whether it is `ρ`-INDEPENDENT (a
  pure number) even though each speed drifts. **Do NOT accept `c_γ=c_s` from shared dimensions or the weak-field
  `em_fields.tex` prose alone.** If the needed gauge/density kinetic-operator data is absent, record
  **`C_GAMMA_RATIO_UNDERDETERMINED`** with the missing operator/profile data + downstream consequence.
- **`(c/c_s)³`:** if `c_γ/c_s` is derived, identify the tail factor `R_tail=Θ_tail(c/c_s)³−1` with it (=1 if they
  coincide) and reconcile with the (legacy-frame) `c=c_s` weak-field statement in `em_fields.tex`; if underdetermined,
  carry the tail factor as conditional on the `c_γ/c_s` residual.
- **Constants vs profiles:** separate genuine constants (`K, m_GNLS, ħ` + the flux `J`/asymptotic `ρ₀, c_{s,0}`, per
  `pathA_19`) from profiles (`ρ(x), v_b(x), c_s(x)`, possibly `c_γ(x)`); state which reference `c_s=1` denotes.

### S2b — the throat as a transonic drain: the `flux_law_verdict` (replaces the "constant `J`" assumption)
`pathA_19` flagged `NO_NET_ACCRETION_BC_UNDERIVED` and did NOT prove the flux is a universal constant. From the coupled
continuity+Bernoulli+EOS profile (S1/S2), resolve the THREE distinct "is it constant?" questions and TEST the transonic
hypothesis rather than assuming it. Produce ONE **`flux_law_verdict`** (one of the three cases below):
  - **Flux conservation ≠ velocity constancy.** In steady state `∇·(ρv_b)=0` makes `J=∮ρv_b·dΣ` the SAME through any
    surface enclosing the drain (Gauss, no sources/leakage) — but `v_b` is NOT constant: as the throat converges
    (`ρ` drops, cross-section narrows) `v_b` ACCELERATES to keep `J` fixed. The throat is a **nozzle**. State both
    explicitly; do not conflate the conserved integral `J` with the local speed `v_b`.
  - **Sonic point = acoustic horizon.** Locate where `v_b` crosses `c_s(ρ)` along the profile. A subsonic→supersonic
    (transonic) drain is the canonical analog black hole (Unruh draining-bathtub); `v_b=c_s` is the horizon (named in
    S2). Determine whether the sonic point sits at/near the mouth.
  - **`flux_law_verdict` ∈ { `TRANSONIC_CHOKED` | `NONTRANSONIC_NO_CHOKED_FLUX` | `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA` }:**
      - `TRANSONIC_CHOKED`: the flux is choked at the sonic throat. Give it as a **critical-value law**
        `J_crit = C_{EOS,geom} · ρ_* c_{s,*} A_*` (local sonic values `ρ_*, c_{s,*}, A_*`, related to the upstream
        `ρ₀, c_{s,0}` by the SOLVED profile), with `C_{EOS,geom}` **DERIVED from the Bernoulli/EOS/geometry** (not just
        upstream dimensional scaling). Machine-check `[J_crit]=T⁻¹` (e.g. `ρ_* c_{s,*} a³` bulk / `ρ_{3,*} c_{s,*} a²`
        brane); record its dependence on the background pressure `P₀=Kρ₀⁵` (so `J`, hence any inflow-mass, is
        environment-dependent — a derived feature, not a frozen input).
      - `NONTRANSONIC_NO_CHOKED_FLUX`: everywhere-subsonic; give the alternate flux law / profile relation that DOES set
        `J`, with its `ρ₀` dependence.
      - `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`: the parent sources give continuity/Euler identities + a
        stationary surface definition but NOT a solved nozzle profile (no branch geometry/confinement/BC). NAME the
        missing data and its effect on `pathA_21`. This is a VALID outcome.
  - `J_in=J_out` (true no-net-accretion) stays gated by the throat-bottom topology (open/closed/connected, undetermined)
    — carry it forward if unresolved.
  - **Honesty requirement:** "transonic choked drain" is a CANDIDATE mechanism — derive the stationary profile and let
    it decide the verdict; do not assert the horizon/choking.

### S3 — `ħ` provenance + the `h`/`2π` decomposition (named fork; feeds the `pathA_21` mass-bridge + M-collapse)
pathA_19 F3 gave `a=ħ/(m_GNLS c_s0)`, i.e. **`ħ = m_GNLS·c_s0·a`** (for the pathA_19 mouth-moment `a`; the `√2` appears
ONLY for the healing length, `ħ = m_GNLS·c_s0·ξ_h/√2`, since `a=ξ_h/√2` — keep `a` and `ξ_h` distinct). So among
`{ħ, m_GNLS, c_s0, a}` only THREE are dimensionally independent; calling `ħ` fundamental and `a` derived (as F3 did) is
a CONVENTION, not physics. Investigate `ħ`'s structure + provenance — TEST, do not declare:
  - **Role catalog (dimension-check each):** quantum of circulation `κ=∮v_b·dl=(ħ/m_GNLS)∮∇θ·dl=(h/m_GNLS)·n`
    (topological, `n∈ℤ`); the phase↔momentum exchange `p=ħ∇θ` (`θ` dimensionless — momentum per radian of phase twist);
    the quantum-pressure / core scale `Q=-ħ²/(2 m_GNLS) ∇²√ρ/√ρ` (sets `ξ_h`).
  - **The `h`/`2π` decomposition (winding-quantization hypothesis — TEST it).** Single-valuedness of `ψ=√ρ e^{iθ}`
    forces phase to wind by `2πn`, so circulation is DISCRETE in steps `h/m_GNLS` (it cannot circulate continuously).
    Determine whether decomposing `ħ=h/(2π)` is physically meaningful here: (a) does `h` (action per complete winding)
    emerge as the natural quantum of the CHARGE/vorticity sector (charge = winding `n`) while mass stays CONTINUOUS via
    the inflow `J` (so same-`n` leptons can differ in mass)? (b) is the inflow/drainage rate `J` a CYCLE-count rate
    (`ν`-like → pairs with `h`, `2π` enters the mass-bridge) or an ANGULAR rate (`ω`-like → `ħ`)? (c) does defect energy
    scale as `n²` (since energy `∝κ²`), making the energy ladder non-uniform (a prediction: higher windings disfavored)?
    Keep `ħ` intact as the GNLS PDE coefficient; expand to `h/(2π)` ONLY inside winding/circulation integrals and
    cycle-vs-radian rate relations where the `2π` is physical. A "no meaningful split" verdict is valid.
  - **The provenance fork — emergence requires an INDEPENDENT `ħ`-free input (ANTI-TAUTOLOGY GATE).** `ħ` is EMERGENT
    ONLY if an independent microscopic/substrate/action relation fixes a length or the action coefficient WITHOUT
    already containing `ħ` (e.g. a substrate spacing sets `a`, making `ħ=m_GNLS c_s0 a` a derived value; or `h` is fixed
    by a substrate winding action). **Rearranging the pathA_19 pin relation `a=ħ/(m_GNLS c_s0)` is NOT a derivation of
    emergence (tautology) — dimension-checking the role catalog is NOT sufficient.** Background-dependence is
    CORROBORATING evidence only: with `c_s0∝ρ₀²`, an `a`-fixed branch makes `ħ` vary with `ρ₀` (⇒ not a true constant ⇒
    emergent); test it. If no independent `ħ`-free relation exists, the valid verdict is **`HBAR_FUNDAMENTAL`** or
    **`HBAR_PROVENANCE_UNDETERMINED`**, each with source + consequence.
  - Record the verdict + reasoning; it is a direct input to `pathA_21`. No base-set change is MADE here (that is a
    `pathA_19`-class decision); S3 only DETERMINES the status and hands it forward.

## Acceptance criteria
1. `c_s` consolidated + machine-verified `[c_s]=L/T` (`c_s²=5Kρ⁴/m_GNLS`) in the pathA_19 base; state-dependence
   (`∝ρ²`) + the `ρ`–`v_b`–`c_s` Bernoulli coupling recorded; EOS-closure provenance stated honestly.
2. Three velocities (`v_b, c_s, c_γ`) named + dimension-pinned; coupling explicit (ledger equations located); genuine
   constants vs profiles separated; `m_GNLS` (constituent) and `m_defect` (throat) kept distinct throughout.
3. `c=c_γ` derived from WAVE-sector first principles (dispersion/group-velocity or trapped-mode Hamiltonian), NOT from
   an assumed `E=m_defect c²`/Compton (a proof starting from `E=m_defect c_γ²` is REJECTED; `E=mc²`/Compton are
   consistency OUTPUTS). The mass-bridge CANDIDATE `m_defect=α_J ħ J/c²` (+ the `h`/`2π` cycle-rate caveat + de Broglie
   note) is RECORDED as a candidate conversion only; `M` NOT collapsed; `α_J` + the physical bridge deferred to
   `pathA_21`.
4. `c_γ/c_s` either DERIVED from the gauge/density kinetic operator (+ `ρ`-(in)dependence + the `(c/c_s)³` tail-factor
   identification) OR recorded `C_GAMMA_RATIO_UNDERDETERMINED` with the missing data; `c_γ=c_s` from shared dimensions /
   weak-field prose alone is REJECTED.
5. (S2b) the transonic hypothesis is TESTED on the derived profile; flux-conservation-vs-velocity-acceleration stated;
   the result is a single `flux_law_verdict` ∈ {`TRANSONIC_CHOKED` with the DERIVED critical-value law
   `J_crit=C_{EOS,geom}·ρ_* c_{s,*} A_*` (`[J]=T⁻¹` machine-checked, `P₀=Kρ₀⁵` dependence) ; `NONTRANSONIC_NO_CHOKED_FLUX`
   with the alternate law ; `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA` naming the missing data}; no-net-accretion
   (`J_in=J_out`) topology status carried.
6. (S3) the `ħ`-provenance verdict ∈ {`HBAR_EMERGENT` (ONLY with a named independent `ħ`-free relation) ;
   `HBAR_FUNDAMENTAL` ; `HBAR_PROVENANCE_UNDETERMINED`} stated with reasoning; the `h`/`2π` decomposition assessed
   (charge↔`h` vs mass↔`J`; `J` cycle-vs-angular; energy`∝n²`); rearranging the pin relation is NOT accepted as
   emergence; no base-set change MADE; verdict handed to `pathA_21`.
7. Dual-engine `.wl` agreement for the ALGEBRAIC/DIMENSIONAL claims ONLY (`[c_s]`, `[c_γ]`, flux dimensions, role-catalog
   dimensions, ratio algebra); non-algebraic physics judgments (transonic existence, `ħ` provenance, profile existence,
   no-net-accretion topology, standing-wave ontology) get human-readable residuals, NOT `.wl` certification. New harness
   group passes; scripts exit 0 within `timeout 600`; pathA_18 + pathA_19 groups untouched.

**Acceptance is PASS/FAIL with NAMED RESIDUALS: script exit-0 is NECESSARY, NOT SUFFICIENT.** Every rejected hypothesis
(non-transonic throat; `c_γ≠c_s` / underdetermined; `ħ` fundamental/undetermined; no meaningful `h`/`2π` split) must
leave a named residual + source + downstream consequence. Negative/`UNDETERMINED` results are VALID outcomes, not
execution failures.

## Out of scope
The mass-bridge DERIVATION + `α_J` + the `2π` resolution + M-collapse test (`pathA_21`, using the `c` + `ħ`-verdict +
`flux_law_verdict` from here); emergent `G` and the `m↔G` unification (`pathA_21`); scale-map → `m̂0²·S_port`
(`pathA_22`); B2c rerun; freeze amendment. Do NOT collapse `M` or change the base set in this step.

## Review (orchestrator, after Codex)
Transliteration-fidelity (one clean agent per new script); independent dimensional re-derivation of `[c_s]`, `[c_γ]`,
the three velocities, the flux dimensions; adversarial pass (distrust-all-clean) targeting: the `c_γ=c_s` vs
`C_GAMMA_RATIO_UNDERDETERMINED` verdict, the standing-wave derivation's NON-CIRCULARITY (no `E=mc²` premise), the
`flux_law_verdict` (transonic not assumed), the `ħ`-emergence ANTI-TAUTOLOGY gate, and the EOS-closure honesty. Claude
reads only residuals. Then gate to `pathA_21`.
