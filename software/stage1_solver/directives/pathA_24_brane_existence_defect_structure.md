# Directive pathA_24 — Testing the "little-arrows" (polar-orientation) brane mechanism: stable wall → light → magnetism → charge/throat → dark energy (T1–T5)

**Status:** DRAFT v2.2 (little-arrows T1–T5 rework, design-review-hardened + GLM-tertiary-folded, 2026-06-23). **Supersedes the generic-domain-wall v1.1** (which was Codex
`SOUND-AS-IS` but mechanism-agnostic — Phase A/B/C over an unspecified "double-well / second OP / self-trapping" candidate set).
v2 replaces that with the **concrete polar-orientation ("little-arrows") mechanism** agreed with the user (2026-06-23,
`docs/conceptual_foundation.md` §2), reorganizes into the **T1–T5 ladder** (cheapest-concept-fatal first), **corrects the charge
framing** (charge = puncture *direction* `±w`, NOT a topological winding / π₂ — the old B2 still imported winding-charge), and
adds **T5, the bulk↔brane / dark-energy cycle**. Because the restructure is substantial, prior Codex `SOUND` did NOT carry.
**Status:** Codex design-review (`-c model_reasoning_effort=xhigh`) ✅ round 1 → `NOT-SOUND` (14 fixes: T0 immutable freeze
artifact + branch rule; T1 no-overclaim-on-imposed-`w`; T1b core-texture-as-output; T2 mechanical shear-vs-Frank discriminator +
fixture battery; T3b alternative-BC battery; T4 stale-source firewall; T4a derive-the-symmetry; T4c derive-operator-before-Coulomb;
T4b/T4d preregistered thresholds; T5 conservation/cosmology scaffold; exhaustive label routing; tier renamed
`CONDITIONAL_SURPLUS_ELIGIBLE` + both-counts surplus; reports-only + integration directive; G10 wall-cosmology gate + G1 hardened)
→ applied → v2.1. Codex confirm-pass → `NOT-SOUND` (3 routing/hygiene residuals) → applied; Codex reconfirm → `NOT-SOUND`
(1 residual: `FAIL_CHARGE_FIREWALL` → CONCEPT-FATAL) → applied; Codex closing reconfirm ✅ `SOUND-AS-IS`. **GLM tertiary ✅
`SOUND-WITH-CONCERNS`** (structurally sound; 10 concerns = calibrate-expectations + tighten 2 gates) → folded → **v2.2**: honest-priors
block (three-way no-win → expect ≤2 of {light, stable wall, emergent-`w`}; emergent-`w` likely lost); **T2 reframed MacCullagh-aware**
(light = rotational-elastic curl-only photon, NOT Cauchy displacement shear; structural pre-check coupled-vs-free-director;
torque-vs-traction discriminator — Claude+GLM resolved: GLM's "T2 pre-decided FAIL" applied generic-LC default over the model's
MacCullagh template); T0 operational single-medium test + framework-circularity caveat; T1c preregistered long-lived threshold; T4a
universality≠quantization + symmetry-triviality flag + label sharpening; P2 operationalized independence; G10 ZKO σ-magnitude
contradiction. **NEXT: Codex confirm-pass on the GLM-driven v2.2 edits** → execute rung-by-rung, each **tri-reviewed**
(orchestrator re-run + transliteration-fidelity audit + adversarial review on clean agents) before its user gate.
**Owner:** orchestrator (Claude) — scaffolding only; **Codex codes/designs all script routes + writes/runs all scripts; Claude
reviews only.** **Engine:** Mathematica leads (wall profiles, fluctuation spectra, stability/Hessian, surface stress, puncture/
throat energetics, Green's-function/Coulomb correspondence, flow balances), SymPy cross-checks; dual-engine where Mathematica can
independently verify; ≤2 concurrent `math -script`; units restored throughout.
**Resume/context (read first):** `docs/conceptual_foundation.md` (canonical conceptual vision — §1 medium, §2 little-arrows wall,
§3 four sectors, §4 defect/throat, §5 dark energy, §7 open questions, §8 this ladder); `decisions/15` §17 (brane = domain wall)
+ §18 (defect = puncture); program state `decisions/13` §0; front door `STATUS.md`. Memory:
`[[project-brane-existence-defect-structure]]`, `[[feedback-native-em-mechanisms]]`, `[[project-em-brane-shear-picture]]`,
`[[project-single-medium-concept-vs-math-drift]]`, `[[reference-conceptual-foundation-doc]]`.

---

## The mechanism under test (concrete — this is what makes v2 testable, not generic)
Give the medium's constituents the **minimal non-trivial property beyond mass: a polar orientation** — each is a *little arrow*
with a head and a tail (**genuinely polar, `+w ≠ −w`; NOT a headless axis/director `n ≡ −n`**). Model it as a **polar
order-parameter field** carried by the one medium (an extra DOF of the same substance, ³He-style — NOT a second medium). The
single mechanism is claimed to deliver, in one shot:
- **The brane** = a wall between an arrows-point-`+w` domain and an arrows-point-`−w` domain (two mirror states of the *same*
  medium); surface tension from the orientation gradient + the medium's cohesion.
- **Why space is 3D** = the arrows align along ONE axis (which *defines* `w`); the surface ⊥ that axis is 3D (codim-1).
- **Light** = at the wall centre the arrows rotate through and **lie flat, in-plane** → the brane carries genuine in-plane
  orientational structure → in-plane **shear** waves (the photon).
- **Magnetism** = in the bulk arrows point `±w` (⊥ our space) → no in-plane structure → **bulk stays shear-free** → Magnus/
  magnetism untouched.
- **Charge** = a throat puncturing toward the `+w` domain vs the `−w` domain → the **two charge signs** (pure puncture
  *direction*; **two domains to puncture into ⇒ exactly two signs**; **no winding, no swirl**).
- **(Ambitious) dark energy** = a bulk↔brane medium cycle throttled by the wall tension (§5 / T5).

This is a **working hypothesis to falsify**, not a result. Each T-rung is an able-to-fail gate; the cheapest concept-fatal one
(T1, stable wall) runs FIRST and gates the rest.

## Falsification stance (load-bearing — [[feedback-falsification-is-the-goal]])
A genuine test; **FAIL is the default expectation until disproven**, and a clean "it all derives" is the *suspicious* outcome to
re-check HARDER. The flagship welcome FAILs, each first-class and never to be rescued:
- **T1:** a polar-arrow field **cannot make a stable codim-1 wall** without ad-hoc structure; OR stability is bought only by
  **imposing** a preferred `w`-axis (so `w` is *not* emergent — the "why 3D" claim collapses); OR the vacuum manifold is a sphere
  so `+w` unwinds to `−w` (no topological protection → no stable wall).
- **T2:** the wall is a **fluid membrane / wrong-mode** wall — tension + bending only (`FAIL_NO_SHEAR_FLUID_MEMBRANE`), or its in-plane modes are
  **Frank/orientational (spin-wave) waves with the wrong dispersion/count**, not the 2-polarization linear-dispersion transverse
  photon (`FAIL_WRONG_MODE_DISPERSION`); the Stage-2/A4 problem then persists.
- **T3:** the in-plane wall stiffness **leaks into the bulk** (bulk acquires shear rigidity) → **breaks Magnus/magnetism**
  (`FAIL_MAGNUS_CONTAMINATED`) — concept-fatal.
- **T4:** the two `±w` punctures are **not mirror-symmetric** (asymmetric charges) or do not exist; the throat balance gives **no
  stable finite radius**; charge comes out **mass-dependent** (breaks universality); brane-tension energy **does not reproduce
  Coulomb**.
- **T5:** the net bulk↔brane flow is **outward** (no expansion) / influx means **densification not area-growth** (drifts `c_s`) /
  no sane decel→accel crossover.
- **Concept-breakers (any rung):** the mechanism is really a **second medium** with no parent-action provenance
  (`FAIL_SECOND_MEDIUM_DRIFT`, near-fatal for single-medium); or a cross-consistency gate (Lorentz, gravity double-count, charge
  firewall, Maxwell double-count) fails.

## Honest priors / expected terminus (calibrate expectations BEFORE executing — GLM tertiary, 2026-06-23)
These do NOT change any gate's pass/fail logic or routing; they record what we *expect* so a partial result is read honestly and a
suspiciously-complete result is re-checked harder. The gates remain able-to-fail in both directions.
- **The three-way no-win (likely near-theorem for a classical polar-vector OP).** Three desiderata pull against each other:
  **(i) light** needs an *orientational* OP (a scalar has no in-plane structure — agreed); **(ii) a stable `±w` wall** needs
  *disconnected* vacua (π₀=ℤ₂), which for a vector OP needs an **easy-axis anisotropy** (O(4)→O(3)); **(iii) emergent `w`** needs
  *full rotational symmetry* (no anisotropy → sphere of vacua → wall not topologically protected). **Any two are reachable; all
  three are likely excluded at once.** The anisotropy that buys stability (ii) is exactly what imposes `w` and kills emergence
  (iii). ⇒ **Expected terminus: at most two of {light, stable wall, emergent-`w`}.** `W_EMERGENT_AND_STABLE` together with
  `WALL_DERIVES_MAXWELL_GRADE_SHEAR` would require *escaping* this near-theorem — treat any such combined result with **maximum
  scrutiny** (likely an error). The realistic best case is `W_IMPOSED_FOR_STABILITY` (stable wall, `w` imposed by an easy axis →
  the "why space is 3D" payoff weakens to "given a preferred axis, 3D follows"). **Plan to lose the emergent-`w` leg.**
- **On light (T2): the question is NOT "does light fall out" — it never would.** A bare scalar superfluid cannot carry light; we
  always intended to *postulate a mode* (the arrows) and **test whether that mode behaves like light**. The real fork (see T2) is
  **clean MacCullagh rotational-elastic photon vs. (Frank orientation wave OR Cauchy-with-stray-longitudinal)** — NOT "displacement
  shear or bust." The little-arrows mechanism is *purpose-built* for the MacCullagh side (the arrows = the rotation-resisting
  Cosserat/gyrostatic elements MacCullagh needed), so T2 is a genuine ~even coin-flip, not a foregone FAIL. The honest hard part is
  getting **curl-only** (no stray longitudinal mode) + gauge/constraint closure — the same wall pathA_23 Stage 2 hit
  (`decisions/15` §15). **Photon packets / quantization (why light comes in localized quanta) are OUT of scope here** — a separate
  later frontier (soliton-in-a-nonlinear-medium is the lead); T2 tests only mode-existence + EM-structure, not quantization.

## Honest classification
A **`NEW_PARENT_ACTION` refinement**: the polar order parameter is an **added DOF of the one medium**, not derivable from the
current single-scalar GNLS (`U(ρ)=(K/4)ρ⁵` is single-well; the arrow field is new structure). So nothing here can be called
"derived from the current action." The correct provenance frame: the *brane/light/charge/dark-energy* are tested as
**derived-from-the-little-arrows-postulate** — i.e. CONDITIONAL on the polar-OP postulate — and the postulate itself is scored for
**independent motivation** (³He precedent, minimality, single-medium fidelity) vs **ad-hoc-ness** (tuned to hit the payoff). The
key non-circular win available is **leverage**: ONE postulated structure (the arrows) yielding **many** otherwise-independent
phenomena (wall + 3D + light's 2-polarization + shear-free bulk + two charge signs) is a genuine surplus *even though the
postulate is not derived* — provided the structure is fixed **target-blind** (T0) and not retro-tuned per payoff. Carry the §14
conceptual costs (classical/mean-field scope; substructure beneath the mean-field — the `[[reference-conceptual-foundation-doc]]`
§1 "infer-substructure-from-effects" lens). Conditional-verdict rule: any postulated/ad-hoc ingredient ⇒ CONDITIONAL; **no
`decisions/14`/paper/verdict updates without explicit user acceptance.**

---

## T0 — ANTI-CIRCULARITY PREREGISTRATION GATE (complete before any T1 solve; binding on all of T1–T5)
The dominant risk is **operational circularity**: tuning the arrow Lagrangian (its potential, anisotropy, stiffness, medium
coupling) *because* it later yields a stable wall / 2 shear modes / two charge signs / `c_γ=c_s`. T0 freezes the rules
**target-blind**, before any payoff is known.
- **Freeze the minimal polar-OP Lagrangian, target-blind.** Write the *simplest* polar (`+w≠−w`) order parameter coupled to the
  one medium — its kinetic/Frank-stiffness term, its self-potential, its coupling to `ρ`/flow — using ONLY: current
  parent-action facts + the single independently-motivated substructure fact ("the constituents carry a polar orientation,
  ³He-style"). Enumerate the **few** discrete modelling choices that must be made (e.g. order-parameter target space: unit vector
  on `S³` vs. a polar field with easy-axis anisotropy vs. a two-component construction; whether `w` is a spontaneously chosen
  alignment direction or an explicit anisotropy axis). **List each choice and its allowed values BEFORE solving.**
- **Emit an IMMUTABLE T0 freeze artifact before any T1 solve** (`reports/pathA_24_T0_freeze.md`, treated as append-only): field
  content; action grammar; **explicitly allowed AND explicitly excluded invariants**; source citations for each term; the full
  discrete-choice list with allowed values; a prior ranking of choices; a **branch budget**; and a content stamp/hash. The freeze
  is the contract every later rung is audited against.
- **Branch rule (no silent down-selection).** Either **run every preregistered branch** through the rungs, OR prune branches
  **only** by target-blind criteria that were **recorded in the freeze artifact before solving**. Pruning a branch *after* seeing
  its (or a sibling's) payoff is forbidden. "We tried the minimal one" is not acceptable unless the freeze artifact shows it was
  the pre-ranked first branch.
- **Forbidden-information rule.** A modelling choice's *admissibility* may NOT cite: wall existence, `σ_brane`, `μ_brane>0`, shear
  mode count, `c_γ`, `λγ≈1`, two charge signs, Coulomb, `e`, `α`, throat size, or any T2–T5 payoff. Admissibility is on intrinsic
  merit only.
- **Framework-circularity caveat (GLM concern 6 — T0 freezes PARAMETERS, not the FRAMEWORK CHOICE).** The choice "polar **vector**
  OP" (vs scalar, tensor, headless director) is *itself* target-motivated (we want in-plane structure → light; two mirror domains →
  two signs). T0's freeze stops parameter-tuning circularity, NOT framework-bias. So the `INDEPENDENTLY_MOTIVATED` provenance class
  requires the **³He precedent (orientational order in a real superfluid) to be cited as motivation that does NOT reference light,
  charge, or any T2–T5 payoff.** If the *only* reason for "polar vector" is "it gives shear and two signs," the provenance degrades
  to `AD_HOC`. (A fully target-blind ideal would enumerate all OP types and rank by intrinsic merit before any payoff; record how
  far the actual freeze falls short of that ideal.)
- **OPERATIONAL single-medium test (GLM concern 3 — makes `FAIL_SECOND_MEDIUM_DRIFT` TESTED, not asserted).** "Carried by the one
  medium" vs "a second field" has a concrete criterion: the OP's parameters (orientational stiffness, self-potential depth,
  medium-coupling) must be **related to the medium's existing parameters** (`ρ0`, `c_s`, `K`), not all independent free constants.
  **If ≥2 OP parameters are independent of the medium's parameters ⇒ flag `FAIL_SECOND_MEDIUM_DRIFT`** (it is functionally a second
  field regardless of "carried-by" rhetoric). The ³He benchmark: ³He's orientational stiffness is set by the gap `Δ` and Fermi
  velocity `v_F` — both *medium* quantities — not by independent OP-specific constants. Record each OP parameter as
  medium-derived / medium-related / independent.
- **Scoring axes (target-blind):** field/DOF count; symmetry & provenance; locality; minimality vs current action; single-medium
  fidelity ([[project-single-medium-concept-vs-math-drift]]) per the operational test above; whether `w` ends up emergent or
  imposed (record, don't optimize).
- **Provenance class for the polar-OP postulate + each choice:** `INDEPENDENTLY_MOTIVATED_NEW_PARENT_ACTION` (³He-grade,
  minimal, target-blind) / `POSTULATED` / `AD_HOC` (introduced/tuned to hit a payoff).
- **Binding rule (hardened).** Any **post-T0 structural change that affects a payoff** — a new/retuned invariant, an added
  anisotropy, a different target space introduced after solving begins — does **not** merely "cap the affected rung." It forces
  **either** a brand-new directive with its own fresh T0 freeze, **or** a hard `AD_HOC_RESCUE` label on the whole run (which bars
  every predictive/surplus claim, not just one rung). The "puzzle = hint about richer substructure" lens
  (`conceptual_foundation.md` §1) may motivate a **future, separately-preregistered** directive — it is **NOT** a licence to
  rescue an in-flight FAIL by adding structure.
Outcomes: `T0_FROZEN(artifact_hash)` / `FAIL_NO_MINIMAL_POLAR_LAGRANGIAN` / `POLAR_OP_AD_HOC` / `AD_HOC_RESCUE` (terminal for all
predictive/surplus claims).

---

## T1 — Does a polar field form a STABLE wall? (the brane make-or-break; gates everything)
If the polar-arrow field cannot self-localize a stable wall, the rest of the ladder is moot — and that is itself the headline
result. Run this rung to completion and gate before T2.
- **T1a — baseline (confirm the gap).** Confirm the current single scalar `U(ρ)=(K/4)ρ⁵` is single-well ⇒ **no wall from the
  scalar alone** (the brane is presently imposed via `V_conf`/`Z`/`W`/`B_ℓ`/`k_w`). Pin precisely what is imposed vs derivable.
  Outcomes: `BRANE_IMPOSED_NOT_DERIVED` (expected) / `BRANE_SELF_LOCALIZES_AS_IS` (surprising — re-check hard).
- **T1b — wall profile + finite surface tension (core texture is an OUTPUT, not an ansatz).** With the T0-frozen polar-OP
  Lagrangian, solve the **unconstrained energy minimum** subject ONLY to the asymptotic boundary data (`+w` as `w→+∞`, `−w` as
  `w→−∞`) — do **not** prescribe a flat/in-plane centre. Derive the profile + a **finite surface tension `σ_brane`** (units
  restored — a 3-surface in 4 spatial dims; see T4-Units), then **classify the core**: does the minimizer actually rotate through
  an **in-plane (flat) core texture** (the precondition T2's light claim needs), or does it pass through `±w` discontinuously / via
  a non-flat path / by amplitude-collapse (`|OP|→0`)? Outcomes: `WALL_PROFILE_WITH_FLAT_CORE(σ_brane)` /
  `WALL_PROFILE_BUT_NO_FLAT_CORE_TEXTURE` (→ T2 little-arrows light route is DEAD) / `FAIL_FLAT_CORE_ASSUMED` (if the solve
  smuggled the flat core in as a constraint — re-do unconstrained) / `FAIL_NO_WALL_PROFILE`.
- **T1c — STABILITY (make-or-break sub-gate; this is GLM's correct a-priori point).** A polar *vector* whose vacua fill a sphere
  gives **no topologically stable wall** — `+w` can continuously rotate ("unwind") to `−w` through the in-plane directions, so the
  wall is not protected and can unwrap. Determine the **minimal** structure that makes the `±w` wall genuinely stable (or
  adequately long-lived): disconnected vacua / **π₀ = ℤ₂** via a **polar double-well or easy-axis anisotropy** along `w`. Compute
  the fluctuation spectrum of the wall: **no negative (unstable) mode**; characterize any near-zero unwinding mode and its
  lifetime. **PREREGISTER the "long-lived" threshold BEFORE the scan (GLM IMPORTANT 4 — mirror T4b; T1c gates T2–T5 so this is the
  highest-leverage spurious-pass hole):** the unwinding lifetime `τ_unwind` must exceed a *now-fixed* threshold (e.g. `τ_unwind ≫
  τ_Hubble` by a stated factor) — chosen and written down before solving, never post-hoc. `WALL_METASTABLE_LONGLIVED` requires
  `τ_unwind` above that threshold; otherwise `FAIL_WALL_UNWINDS_SPHERE_VACUA`. Outcomes: `WALL_TOPOLOGICALLY_STABLE(π₀=ℤ₂)` /
  `WALL_METASTABLE_LONGLIVED(τ_unwind)` / `FAIL_WALL_UNWINDS_SPHERE_VACUA` / `FAIL_WALL_UNSTABLE_MODE`.
- **T1d — is `w` EMERGENT or IMPOSED? (honesty sub-gate — do not paper over).** The "why space is 3D" payoff requires the
  alignment axis to be **spontaneously chosen** (the arrows pick a direction; `w` is *defined* by it). But T1c's stability often
  wants an **explicit easy-axis anisotropy** — which **imposes** a preferred `w`. These pull against each other; **report which
  actually happened**. Isotropic-stiffness + spontaneous pick ⇒ `w` emergent **but** sphere-of-vacua instability (T1c risk);
  easy-axis anisotropy ⇒ stable **but** `w` imposed (the "why 3D" claim weakens to "given a preferred axis, 3D follows"). Outcomes:
  `W_EMERGENT_AND_STABLE` (reserved ONLY for a demonstrated spontaneous-axis mechanism with genuinely disconnected `+w/−w`
  sectors and **no hidden fixed background vector** — the strong win, re-check hard) / `W_IMPOSED_FOR_STABILITY(anisotropy)` /
  `W_EMERGENT_BUT_UNSTABLE`.
- **T1e — confinement (why we don't flow into the bulk).** Show the wall is a potential well **binding zero modes** at `w≈0` (us
  and our matter); the energy to leave = the confinement gap. Connect to the existing `V_conf`/`Z`/`W`/`B_ℓ`/`k_w` as the
  *effective* description this wall would derive. Outcomes: `CONFINEMENT_FROM_WALL` / `FAIL_NO_BOUND_ZERO_MODES`.
**T1 verdict roll-up (no overclaim on imposed `w`):** `T1_STABLE_WALL_DERIVED` requires T1b=`WALL_PROFILE_WITH_FLAT_CORE` + T1c
stable + T1e bound modes + **T1d=`W_EMERGENT_AND_STABLE`**. Any of {`W_IMPOSED_FOR_STABILITY`, explicit easy-axis anisotropy, any
anisotropy beyond the independently-motivated T0 postulate} ⇒ `T1_STABLE_WALL_CONDITIONAL` or
`FAIL_STABILITY_REQUIRES_EXPLICIT_AXIS` — **never** `T1_STABLE_WALL_DERIVED` (an imposed easy axis loses the "why space is 3D"
payoff). Metastable ⇒ `T1_STABLE_WALL_CONDITIONAL`. No wall ⇒ `T1_FAIL_NO_STABLE_WALL`.

## T2 — Does the wall carry in-plane SHEAR → light? (the light make-or-break)
Only meaningful if T1 produced a wall (else MOOT). Compute the spectrum of fluctuations **localized on the wall** and ask whether
the in-plane arrow texture yields the photon. **The correct fork is MacCullagh-clean-light vs (Frank OR Cauchy) — NOT "Cauchy
displacement shear or bust."** Light in this model is the **MacCullagh rotational-elastic** mode: energy in the *local rotation
(curl)* of a genuine medium displacement, `∝(∇×u)²` → **exactly two transverse polarizations, NO longitudinal mode**, Maxwell
structure. That is *why* the postulated arrows (rotation-resisting Cosserat/gyrostatic elements) are the candidate: ordinary matter
resists compression and Cauchy shear but **not** local rotation, so a medium that resists rotation needs orientable substructure —
the arrows. (Note: a **Cauchy** displacement-shear wall is NOT the goal — it leaves a stray longitudinal "second photon"
`FAIL_CAUCHY_STRAY_LONGITUDINAL`; and a **pure Frank** director wave is not light either — it exerts torque, not traction.)
Confronts pathA_23's Stage-2 lesson (`decisions/15` §15) + GLM's a-priori objection head-on.
- **Structural pre-check (do this FIRST — GLM Q1, the decisive a-priori fork).** Inspect the frozen action: is the arrow's in-plane
  rotation **rigidly coupled to a genuine displacement of the medium** (tilting an arrow *is* a local rotation of the material →
  rotational elasticity is possible → MacCullagh route live), or is the OP a **free, decoupled director** whose only energy is
  orientation-gradient `(∇n)²` (→ a pure **Frank/spin wave** that exerts torque not traction → **pre-decided
  `FAIL_WRONG_MODE_DISPERSION`**, regardless of mode count or dispersion)? A decoupled director short-circuits T2 to FAIL. This is
  the same wall as pathA_23 Stage-2 `FAIL_UNSPECIFIED_SUBSTRUCTURE`: the medium must supply a genuine displacement for the rotation
  to act on. ("Carried by the one medium" — the T0 single-medium test — is precisely the coupling that *could* turn Frank into
  MacCullagh; this is why it is mechanically load-bearing, not just philosophical.)
- **Require ALL of (the photon signature):** (i) **exactly two** transverse in-plane modes; (ii) **linear, non-dispersive**
  dispersion `ω = c_γ k` at a single speed `c_γ`; (iii) the stiffness is **rotational/curl-type (MacCullagh)** — energy in
  `(∇×u)²` giving **no physical longitudinal photon** — and explicitly **not** Cauchy symmetric-strain elasticity (which leaves a
  stray longitudinal mode, `FAIL_CAUCHY_STRAY_LONGITUDINAL`); (iv) **no massless flexural/bending fifth-force mode** masquerading
  as light (the drumhead `u_w` bending mode must be a separate, identified `ω ∝ k²` scalar, not a second photon); (v) full-action
  **gauge/constraint closure** — a curl-only *potential energy* alone is NOT enough; the pathA_23 C5 obstruction must actually
  close; (vi) **source-coupling compatibility** with the charge firewall.
- **Three-way discriminator (Frank vs Cauchy vs MacCullagh — decide from the quadratic action, not by eyeballing dispersion).**
  Build: (a) the **DOF classifier** from the quadratic action's field content; (b) the **Noether momentum** / surface stress;
  (c) the decisive **torque-vs-traction test** — does an in-plane wiggle exert a genuine mechanical **traction** (force/area across
  a cut surface, conjugate to *displacement* → elastic: MacCullagh or Cauchy) or only an orientational **torque** (conjugate to
  *orientation* → Frank, not light)?; (d) the **principal symbol** of the mode operator (**curl-type → MacCullagh, transverse-only**;
  **symmetric-strain → Cauchy, with longitudinal**). PASS = **MacCullagh rotational-elastic** (traction + curl-only + two
  transverse + gauge-closed). Torque-only ⇒ Frank ⇒ `FAIL_WRONG_MODE_DISPERSION`; symmetric-strain ⇒ Cauchy ⇒
  `FAIL_CAUCHY_STRAY_LONGITUDINAL` — even if count/dispersion look right.
- **Negative-control fixture battery (each must return the stated result or the discriminator is tautological,
  [[feedback-decisive-test-not-tautological]]):** a **fluid membrane** → `FAIL_NO_SHEAR_FLUID_MEMBRANE`; a **free Frank/director
  wall** (decoupled OP) → `FAIL_WRONG_MODE_DISPERSION`; an **ordinary Cauchy elastic wall** → `FAIL_CAUCHY_STRAY_LONGITUDINAL`; a
  **curl-only-potential-energy wall WITHOUT full gauge closure** (the pathA_23 C5 obstruction) → `FAIL_CURL_ONLY_NOT_GAUGE`; and —
  so the test can SUCCEED, not only fail — a **genuine gauge-closed MacCullagh wall** → `WALL_DERIVES_MAXWELL_GRADE_SHEAR`. Require
  the explicit constraint/scalar-potential analog that **removes the longitudinal zero mode** (condition iii).
Outcomes: `WALL_DERIVES_MAXWELL_GRADE_SHEAR(μ_brane, c_γ, constraints)` (all six hold + the discriminator confirms **MacCullagh
rotational elasticity** — would **retire** the pathA_23 MacCullagh postulate by *deriving* it from the arrow substructure) /
`SHEAR_POSTULATED_STILL(→CONDITIONAL)` / `FAIL_NO_SHEAR_FLUID_MEMBRANE` / `FAIL_WRONG_MODE_DISPERSION` (Frank/spin-wave, torque not
traction) / `FAIL_CAUCHY_STRAY_LONGITUDINAL` / `FAIL_CURL_ONLY_NOT_GAUGE` / `FAIL_BENDING_MASSLESS_FIFTH_FORCE`.

## T3 — Bulk stays shear-free → magnetism preserved + separate-sector legitimacy
Only meaningful if T2 produced in-plane shear (else MOOT/post-mortem). Confirm the in-plane stiffness is **confined to the wall**
and does **not** leak into the bulk.
- **T3a — bulk shear-free.** With arrows `±w` in the bulk, confirm the bulk acquires **no in-plane shear rigidity**; the
  Magnus/vorticity (magnetism) mechanism is untouched. Outcomes: `BULK_SHEAR_FREE_MAGNUS_PRESERVED` / `FAIL_MAGNUS_CONTAMINATED`
  (concept-fatal).
- **T3b — separate-sector legitimacy (the pathA_23 Stage-3b answer, now DERIVED here not inherited — and it must be EARNED against
  the alternatives, not declared).** Stage-3b was explicitly *model-contingent*: if the brane `u^a` and bulk `v^a` are separated
  by fiat, nothing is earned. So run a **mandatory alternative-boundary-condition battery** on the *same* derived wall and report
  which survive: (1) **free-slip membrane-in-fluid**; (2) **strict single-medium no-slip matching `u̇^a = v^a`**; (3) **mixed
  throat-curvature** cases. From the derived wall: derive the surface collective coordinates; test tangential free-slip / no-locking
  between `u̇^a` and `v^a`; determine whether the wall's internal stress sources **only** the wall EOM or **also** bulk vorticity;
  carry the throat-curvature leak as an open/fail channel (relocated-not-retired from Stage-3b). Only if a single-medium-legitimate
  matching (not a declared separation) gives no-leak ⇒ `SEPARATE_SECTOR_DERIVED`; if only the *declared* separate-sector action
  passes ⇒ `SEPARATE_SECTOR_CONDITIONAL`; if strict matching reinstates the leak ⇒ `FAIL_STRICT_FLUID_LOCKING_LEAK_RETURNS`.
  Outcomes: `SEPARATE_SECTOR_DERIVED` / `SEPARATE_SECTOR_CONDITIONAL` / `FAIL_STRICT_FLUID_LOCKING_LEAK_RETURNS`. **(`FAIL_MAGNUS_CONTAMINATED` from T3a is concept-fatal — see routing table.)**

## T4 — Puncture into `±w` domains → two charge signs + finite throat
Only meaningful if T1 produced a wall (needs `σ_brane`; the charge/throat sub-rungs need `σ_brane`, NOT `μ_brane`, so they may run
even if T2 fails — flag light/`λγ` as DEAD in that branch). **This rung corrects the old B2 framing: charge is the puncture
DIRECTION (which `±w` domain the throat drains into), NOT a topological winding / π₂ classification. Do NOT compute a winding
number or invoke a vacuum-manifold homotopy charge — that conflates electric (puncture direction) with magnetic (swirl).**
**STALE-SOURCE FIREWALL (binding on executors):** `decisions/15` §18 and `[[project-brane-existence-defect-structure]]` still
contain the *old* "charge = winding/orientation / charge from topology" phrasing in places. Those passages are **SUPERSEDED**
here: `η_Q` maps ONLY to puncture direction into the `+w` vs `−w` domain. Any winding / circulation / swirl belongs to
**magnetism/spin** and is a **FAIL for electric charge** (`FAIL_DIRECTION_NO_FIELD_COUPLING` if a winding is what couples to the
charge field). The canonical conceptual source is `docs/conceptual_foundation.md` §3–§4, not the stale lines.
- **T4-Units — Units/Jacobian gate (terminal for T4 charge/Coulomb arithmetic).** Before any Coulomb/self-energy/calibration
  number: the brane is a **3-surface in 4 spatial dims**, so `σ_brane` is a 3D energy density and the `σ→charge` map is dominated
  by measure/Jacobian conventions (cf. the `m̂0²·S_port` natural-units saga — one missing volume factor manufactured a false
  result). Pin dimensions of `σ_brane`, surface energy, deformation amplitude, `ε₀/μ₀`, `e`, `Z_int`; brane- vs bulk-volume
  measures; Jacobians near the throat; consistency with the existing localized Maxwell term. Outcomes: `UNITS_CONSISTENT` /
  `UNITS_OR_JACOBIAN_FAIL` (terminal).
- **T4a — two charge signs from the two domains (DERIVE the symmetry; do NOT bake it in from the words "mirror states").** The
  two symmetric signs must be an **output of the same action**, not labels imposed by the setup (else `FAIL_NO_BOTH_SIGNS` /
  `FAIL_ASYMMETRIC_SIGNS` are unreachable). Concretely: construct **two finite-energy throat solutions of the one frozen action** —
  one draining into the `+w` domain, one into the `−w` domain — with **boundary data related by the domain-exchange symmetry**, and
  **derive** that the charge-coupling coefficient is equal in magnitude for the two (do not assume it from "mirror"). **If T0
  admitted any parity-breaking polar coupling, it must be run (or have been target-blind-pruned in the freeze artifact)** — a
  parity-breaking coupling is exactly what would yield `FAIL_ASYMMETRIC_SIGNS`, so it cannot be quietly dropped. Then check: the
  sign is **independent of throat radius / standing-wave energy / circulation**; it maps to the existing `q_eff = η_Q e_*/√Z_int`
  ontology + charge firewall (`pde.tex:279-312`) + mixed-core.
  - **Honest scope of "the symmetry derivation" (GLM concern 9).** If the frozen action is `w→−w` symmetric (natural for a polar OP
    with an even potential), the two-sign symmetry follows by that symmetry **by construction** — a one-line kinematic consequence,
    NOT a deep derived result; **flag it as such.** The gate's real bite is therefore (a) running any parity-breaking coupling T0
    admitted (→ `FAIL_ASYMMETRIC_SIGNS`), and (b) the field-coupling check below.
  - **Universality ≠ topological quantization (GLM concern 8 — state plainly).** Discarding the winding/π₂ picture **relocates**, it
    does not *resolve*, the quantization question. The puncture-direction mechanism gives **sign quantization** (binary `±w`,
    genuine) and **magnitude universality** (one `σ_brane` → one `|e|` for all punctures) — but **NOT** a topologically fixed
    magnitude: `σ_brane` is a *continuous* parameter, so `|e|` is a continuous prediction fixed by **calibration** (P2), not by a
    topological invariant. And if the wall is only metastable (T1c), **charge conservation is only metastable** too (a `+w` puncture
    could flip if the wall locally unwinds). Do not let "charge = direction" be read as supplying quantization of the magnitude.
  - **Sharpen the two field-coupling failure labels (GLM concern 10).** `TWO_LABELS_ONLY_NO_FIELD_COUPLING` (CONDITIONAL) = the
    puncture direction is a bare label that **couples to no field** (the charge sector is *empty*). `FAIL_DIRECTION_NO_FIELD_COUPLING`
    (SUBCLAIM-TERMINAL) = the puncture exists but the **wrong mechanism couples** — a winding/circulation/swirl is what sources the
    charge field instead of the direction (the charge sector is *misidentified*, i.e. it collapsed back onto magnetism). Both fail
    the charge claim; route by empty-vs-misidentified.
  Outcomes: `TWO_SIGNS_FROM_DIRECTION_DERIVED` / `FAIL_NO_BOTH_SIGNS` / `FAIL_ASYMMETRIC_SIGNS` / `TWO_LABELS_ONLY_NO_FIELD_COUPLING`
  (signs exist kinematically but do not source the charge ontology) / `FAIL_DIRECTION_NO_FIELD_COUPLING` / `FAIL_CHARGE_FIREWALL`.
- **T4b — stable finite throat size (PREREGISTER the thresholds before solving — no rhetorical softening of a marginal throat).**
  Surface-tension **closing** stress/energy vs the **holding-open** agent (geon standing-wave pressure and/or drain flow); solve
  the balance. A stationary point is NOT enough. **Before the scan, record in the report:** (1) the **Hessian** must be
  positive-definite over radius *and* profile modes; (2) energy **bounded** as `a→0` and `a→∞`; (3) a finite **breathing-mode
  frequency** ω_breath>0 with no runaway; (4) a **geon leakage time τ_leak** with the preregistered pass band — `τ_leak` compared
  to the relevant dynamical/cosmological scale, with the "long-lived" threshold fixed *now* (e.g. `τ_leak ≫ τ_dyn` by a stated
  factor), not chosen post-hoc (Wheeler geons are at best quasi-stable — confront it). Outcomes: `THROAT_STABLY_DERIVED(a)` /
  `METASTABLE_THROAT_CONDITIONAL` / `FAIL_UNSTABLE_THROAT` / `FAIL_GEON_LEAKS_TOO_FAST` / `FAIL_NO_FINITE_BALANCE`.
- **T4c — Coulomb must be EARNED (derive the OPERATOR before assuming radial symmetry/Laplacian).** First derive the **static
  field operator** from the wall/throat action + dimension-restored geometry — do NOT assume a 3D Laplacian, radial symmetry, or
  `φ~1/r` (any such assumption before solving makes Coulomb nearly automatic = circular; flag it). Classify the derived operator:
  a **massive, anisotropic, higher-gradient, nonlocal, or mixed bulk-brane** operator routes to an explicit non-Coulomb outcome.
  Only a genuine massless brane-3-space Laplacian gives Coulomb. Then test BOTH: single-puncture field energy + **finite
  self-energy `~q²/a`** (using `a` from T4b); two-puncture interaction incl. **sign for like/unlike** punctures + the **`1/r²`
  force**. **Negative control:** feed a known **non-Laplacian** operator and confirm the classifier returns `FAIL_NOT_COULOMB`.
  **Dependency split:** if T2/T3 failed (no Maxwell-grade light sector), T4c may at most be a *static conditional Coulomb test*,
  NOT a Maxwell-grade EM claim. Outcomes: `TENSION_REPRODUCES_COULOMB` (operator-derived + falloff + sign + finite self-energy) /
  `FAIL_NOT_COULOMB` / `FAIL_WRONG_FORCE_SIGN` / `SELF_ENERGY_STILL_DIVERGENT`.
- **T4d — charge ⊥ mass decoupling (PREREGISTER tolerance/order; bounded, not assumed).** T4b/T4c create *indirect* coupling (geon
  pressure sets `a`; `a` sets self-energy), so define allowed vs forbidden first: allowed = mass-dependence of throat
  size/self-energy; **forbidden = mass-dependence of charge** (`q_eff`, `η_Q`, the Coulomb coefficient). Before the scan, **declare
  the tolerance and order**, then build a **sensitivity matrix** of charge observables vs geon energy/radius (hold puncture
  direction fixed, vary geon energy). Outcomes: `DECOUPLING_HOLDS_TO_ORDER(order)` / `DECOUPLING_CONDITIONAL_SMALL_COUPLING` /
  `FAIL_CHARGE_MASS_DEPENDENT`.
- **T4e — defect-structure ↔ sign preference (open, exploratory; NOT a pass criterion).** Probe the §4 hypothesis: does a simple
  drain (electron) vs a knot (proton) **energetically prefer** one `±w` puncture direction via bulk interaction — with the
  opposite orientation the antiparticle (so the preference is energetic/statistical, NOT a rigid lock; antiparticles must remain
  allowed)? Record as `SIGN_PREFERENCE_OBSERVED(ΔE)` / `NO_SIGN_PREFERENCE` / `FAIL_FORBIDS_ANTIPARTICLES` (the matter-abundance
  asymmetry is explicitly OUT of scope — a separate question).

## T5 — The bulk↔brane cycle (dark energy; ambitious — only after T1–T4)
Only meaningful if T1 produced a wall with tension (else MOOT). **Also MOOT/post-mortem if T4b failed** (`FAIL_UNSTABLE_THROAT` /
`FAIL_GEON_LEAKS_TOO_FAST` / `FAIL_NO_FINITE_BALANCE`) — the matter-driven drainage presupposes stable throats. Set up the two
opposing medium flows and test the §5 claims. Signs and rates are everything; **failure here does NOT retract T1–T4.** This is the
easiest rung to produce plausible *words* with no checkable accounting — so it is gated on a real conservation/cosmology model, not
narrative.
- **T5-Model — REQUIRED accounting scaffold (before T5a–T5c claims).** Write down, units restored: (1) a **continuity equation**
  for medium on/off the brane (the conserved medium splits between bulk and brane reservoirs); (2) an explicit **map from areal
  density to the scale factor `a(t)`** (how "wall area grows" becomes expansion); (3) **`c_s(ρ)` drift bounds** (the densification
  alternative is a real, observable failure — bound it); (4) the **wall-tension stress-energy / gravitational backreaction** (the
  wall has its own energy budget — see G10); (5) declared **numerical sanity bands** for `Ω_DE` and the crossover redshift. Outcomes:
  `T5_MODEL_CLOSED` / `FAIL_NO_CONSERVATION_MODEL`.
- **T5a — the two flows.** (1) Matter-driven **brane→bulk drainage** (gravity = medium draining down throats — the **same** drain
  as G5/G9, not a second book-keeping), scaling with **matter content**. (2) Tension-throttled **bulk→brane areal leak** (bulk
  drawn toward matter; arrows must rotate `±w`→flat to enter, which the wall tension throttles into a slow steady leak), scaling
  with **brane area**. Derive both rates with units restored. Outcomes: `BOTH_FLOWS_DERIVED` / `FAIL_NO_INWARD_LEAK_MECHANISM`.
- **T5b — net sign + expansion.** Is the **net** flow **inward** (medium accumulating on the brane)? Our space *is* the wall;
  does inflow at ~fixed areal density force **area growth** (expansion) rather than **densification** (which would instead drift
  `c_s` over cosmic time — a distinct, checkable, and potentially falsifying alternative)? Outcomes: `NET_INWARD_AREA_GROWTH` /
  `FAIL_NET_OUTWARD` / `DENSIFICATION_NOT_EXPANSION(c_s drift)`.
- **T5c — acceleration + decel→accel crossover.** Influx ∝ **area** ⇒ exponential/de Sitter-like **acceleration** (dark-energy
  signature); drainage ∝ **matter** ⇒ early (matter-dense) **deceleration**, late (matter-dilute) **acceleration** — does the
  crossover timing come out even roughly sane? Outcomes: `ACCEL_CROSSOVER_PLAUSIBLE` / `FAIL_NO_ACCELERATION` /
  `FAIL_CROSSOVER_TIMING_INSANE`. (CONDITIONAL/exploratory throughout — T5 is the most speculative rung; never bank it.)

---

## PAYOFF BOOKKEEPING — verdict feed-forward + calibrate-predict (gated; predictive ONLY per the routing table below)
- **P1 — feed-forward to the verdict.** ONLY if T2 derived `μ_brane` (Maxwell-grade), T3 derived the separate sector, and T4 set a
  stable throat: feed `c_γ²=μ_brane/ρ_brane` + `σ_brane` into **`λγ`** and the GR-quadrupole verdict bookkeeping; update the value
  map ONLY after a full pass + explicit user acceptance of any CONDITIONAL status.
- **P2 — calibrate-predict: brane tension ↔ electric charge (surplus/rank accounting; the Gate-4 lesson).** The EM coupling is
  currently un-derived (Gate 4: gravity `g_G` absorbs `G` with zero surplus; `λγ` a gap). IF the puncture's tension-energy ties to
  charge, calibrate the single universal `σ_brane` on `e`/`α` and predict held-out values — BUT only if the accounting is
  predictive (else this is just another absorbed anchor).
  - **FIRST GATE — rank/surplus accounting (before ANY calibration number).** Freeze calibration constants, nuisance constants,
    and discrete branch choices *before* fitting; compute the **rank of independent held-out observables after constraints**.
    **Report BOTH counts** (the Gate-4 lesson — a NEW_PARENT_ACTION can manufacture apparent surplus): (i) **surplus within the
    frozen arrow action** (held-out observables beyond the arrow-postulate's own knobs); (ii) **surplus relative to the current
    scalar parent action** (what the arrow postulate buys *over* the pre-existing model — this is the honest "is the new structure
    earning its keep" number). Classify: `EM_COUPLING_ABSORPTIVE_NO_SURPLUS` / `EM_COUPLING_CONSTRAINED_SURPLUS(rank=N)` /
    `EM_COUPLING_PREDICTIVE_SURPLUS(rank_within=N, rank_vs_scalar=M)`. If ≥3 independent moduli knobs remain ⇒ `NO_CASCADE_ABSORPTIVE`.
    Settle BEFORE quoting any number. **Leverage prerequisite:** the cascade requires T1/T2's wall structure to **relate** the
    moduli (`σ_brane`, `μ_brane`→`c_γ`, `ρ_brane`); three independent knobs ⇒ no cascade. **No `decisions/14`/paper update may follow
    from a conditional surplus without explicit user acceptance.**
    - **OPERATIONALIZE "independent" (GLM concern 7 — "rank" is gameable by eyeballing).** Independence must be checked by
      **dimensional analysis + symmetry relations, not visual inspection**: (a) list **ALL** parameters of the frozen action (OP +
      medium), marking each medium-derived / medium-related / **new-independent**; (b) list **ALL quantitative observables**
      (NUMBERS — exclude qualitative structural outcomes like "3D" / "two signs" / "wall exists"); (c) check every observable pair
      for dimensional/symmetry relations that make them dependent (e.g. `c_γ` and `λγ=c_γ/c_s` are **not** independent if `c_s` is
      known); (d) **honest vs-scalar surplus = (independent quantitative observables) − (new independent parameters) − (calibration
      anchors used).** Adding the polar OP introduces ~3–4 new parameters, so even 5–6 explained *phenomena* can net **0–1**
      quantitative surplus after calibrating on `e` — that is absorptive, report it as such. **Qualitative structural outcomes
      (wall, 3D, two signs, shear-free bulk) count as *leverage*, NOT quantitative surplus.**
  - **Held-out checklist — TIERED (target-blind; AFTER calibration):**
    - **In-scope for pathA_24:** charge universality (mass-independent `e`); Coulomb coefficient + `1/r²` + force sign; finite
      self-energy scale; stable throat radius `a` (*which* scale — `r_e≈2.8 fm` vs Compton `λ_C`, differing by `α` — is itself
      diagnostic); `λγ=c_γ/c_s≈1` (GW170817 cone-lock) **only if** T1/T2 derive the moduli relation.
    - **Downstream (aspirational, NOT pathA_24 pass criteria):** lepton mass ratios (need a quantized geon-spectrum directive);
      `g−2` (spin/anomaly machinery); `α≈1/137` (a separate brane/bulk moduli-ratio derivation — moonshot, do NOT bank).

## DEPENDENCY ROUTING — explicit T-outcome → downstream eligibility (evaluate before each rung & before payoff)
Route by the *actual* emitted labels (no stale names). The ladder is strictly gated: a rung's FAIL does not auto-kill lower rungs
that need only its upstream products, but it does cap what can be *claimed*. **The authoritative per-label router is the
"EXHAUSTIVE LABEL ROUTING" table below; the surplus-eligibility gate is stated first.**
- **`CONDITIONAL_SURPLUS_ELIGIBLE`** (formerly "PREDICTIVE-ELIGIBLE"; renamed so the tier never reads as an ordinary, unconditioned
  prediction — the polar-OP is a NEW_PARENT_ACTION, so any surplus is **conditional on the single arrow postulate**, never "derived
  from the current action"). P1/P2 may be called *surplus* ONLY if ALL hold: `T0_FROZEN` (artifact hash, no `AD_HOC_RESCUE`, polar-OP
  `INDEPENDENTLY_MOTIVATED_NEW_PARENT_ACTION`, no post-payoff retuning); `T1_STABLE_WALL_DERIVED`; `CONFINEMENT_FROM_WALL`;
  `WALL_DERIVES_MAXWELL_GRADE_SHEAR`; `SEPARATE_SECTOR_DERIVED`; `THROAT_STABLY_DERIVED`; `TWO_SIGNS_FROM_DIRECTION_DERIVED`;
  `TENSION_REPRODUCES_COULOMB`; `DECOUPLING_HOLDS_TO_ORDER`; PLUS the moduli-relation leverage prerequisite (P2) and the
  **both-counts** surplus report (within-arrow-action AND vs-scalar-parent). Missing any ⇒ at best CONDITIONAL_EXPLORATORY.

**EXHAUSTIVE LABEL ROUTING.** Routing conventions (so the table is exhaustive by rule, not by listing every parametrized variant):
**(C1)** `LABEL(args)` routes identically to `LABEL` (the parenthetical payload is informational). **(C2)** Every rung's stated
**success/PASS outcome** (the `*_DERIVED` / `*_OK` / `*_PRESERVED` / `*_CONSISTENT` / `*_FOUND` / `T0_FROZEN` / `UNITS_CONSISTENT` /
`WALL_PROFILE_WITH_FLAT_CORE` / `*_CLOSED` / `BOTH_FLOWS_DERIVED` / `NET_INWARD_AREA_GROWTH` / `ACCEL_CROSSOVER_PLAUSIBLE` /
`SIGN_PREFERENCE_OBSERVED` / `NO_SIGN_PREFERENCE` outcomes) → **PROCEED** (carry the product forward; this is *not* a surplus claim
— surplus eligibility is the separate gate above). **(C3)** Any FAIL/conditional label **not explicitly classified below defaults to
CONCEPT-FATAL** (never optimistic). The non-PASS labels are classified into exactly four buckets:
- **CONCEPT-FATAL — stop ALL model claims; record the welcome FAIL plainly (off-concept, never "the model"):**
  `FAIL_NO_MINIMAL_POLAR_LAGRANGIAN`, `AD_HOC_RESCUE`, `FAIL_SECOND_MEDIUM_DRIFT`, `T1_FAIL_NO_STABLE_WALL`,
  `FAIL_WALL_UNWINDS_SPHERE_VACUA`, `FAIL_WALL_UNSTABLE_MODE`, `FAIL_FLAT_CORE_ASSUMED`, `FAIL_NO_WALL_PROFILE`,
  `FAIL_NO_BOUND_ZERO_MODES`; `FAIL_MAGNUS_CONTAMINATED`, `FAIL_PREFERRED_FRAME`, `FAIL_THROAT_ENERGY_DOUBLE_COUNT`,
  `FAIL_DRAIN_PICTURE_BROKEN`, `FAIL_CORE_ONTOLOGY_LOST`, `FAIL_EM_DOUBLE_COUNT`, `FAIL_DRAIN_DOUBLE_COUNT`, `FAIL_WALL_OVERCLOSES`,
  `FAIL_SIGMA_BACKREACTION_FATAL`, `FAIL_CHARGE_FIREWALL` (the charge firewall is a declared concept-breaker — a violation means the
  puncture picture is inconsistent with the existing charge ontology, not merely that one subclaim failed).
- **SUBCLAIM-TERMINAL — the *named payoff* dies, but the ladder CONTINUES (independent subclaims that need only upstream products
  may still run, conditionally); never a surplus claim:**
  - **kills LIGHT/`λγ`** (T4a–T4d charge/throat may still run — they need `σ_brane`, not `μ_brane`): all T2 FAILs
    (`FAIL_NO_SHEAR_FLUID_MEMBRANE`, `FAIL_WRONG_MODE_DISPERSION`, `FAIL_CAUCHY_STRAY_LONGITUDINAL`, `FAIL_CURL_ONLY_NOT_GAUGE`,
    `FAIL_BENDING_MASSLESS_FIFTH_FORCE`); `WALL_PROFILE_BUT_NO_FLAT_CORE_TEXTURE`.
  - **kills the separate-sector/no-leak claim** (near-fatal — flag hard): `FAIL_STRICT_FLUID_LOCKING_LEAK_RETURNS`.
  - **terminates all T4 charge/Coulomb arithmetic AND P2:** `UNITS_OR_JACOBIAN_FAIL`.
  - **kills the CHARGE claim:** `FAIL_NO_BOTH_SIGNS`, `FAIL_ASYMMETRIC_SIGNS`, `FAIL_DIRECTION_NO_FIELD_COUPLING`,
    `FAIL_CHARGE_MASS_DEPENDENT`. (`FAIL_CHARGE_FIREWALL` is **CONCEPT-FATAL** — see above, not here.)
  - **kills the THROAT/MASS claim (⇒ T5 also MOOT):** `FAIL_UNSTABLE_THROAT`, `FAIL_GEON_LEAKS_TOO_FAST`, `FAIL_NO_FINITE_BALANCE`,
    `FAIL_GEON_TRANSIENT`.
  - **kills the COULOMB claim:** `FAIL_NOT_COULOMB`, `FAIL_WRONG_FORCE_SIGN`, `SELF_ENERGY_STILL_DIVERGENT`.
  - **kills the T5/dark-energy claim only (T1–T4 untouched):** `FAIL_NO_CONSERVATION_MODEL`, `FAIL_NO_INWARD_LEAK_MECHANISM`,
    `FAIL_NET_OUTWARD`, `DENSIFICATION_NOT_EXPANSION`, `FAIL_NO_ACCELERATION`, `FAIL_CROSSOVER_TIMING_INSANE`.
  - **kills only the T4e sign-preference *hypothesis* (exploratory; ladder untouched):** `FAIL_FORBIDS_ANTIPARTICLES`.
- **CONDITIONAL_EXPLORATORY — rung may run; NO `decisions/14`/paper/verdict updates:** `POLAR_OP_AD_HOC`,
  `T1_STABLE_WALL_CONDITIONAL`, `WALL_METASTABLE_LONGLIVED`, `W_IMPOSED_FOR_STABILITY`, `W_EMERGENT_BUT_UNSTABLE`,
  `FAIL_STABILITY_REQUIRES_EXPLICIT_AXIS` (a FAIL of the *derived* claim — run T2+ only as exploratory), `SHEAR_POSTULATED_STILL`,
  `SEPARATE_SECTOR_CONDITIONAL`, `METASTABLE_THROAT_CONDITIONAL`, `DECOUPLING_CONDITIONAL_SMALL_COUPLING`,
  `TWO_LABELS_ONLY_NO_FIELD_COUPLING`, `EFFECTIVE_LORENTZ_SUPPRESSED_CONDITIONAL`, `WALL_BACKREACTION_CONDITIONAL`,
  `MAXWELL_DERIVED_AS_EFFECTIVE` (T4c static-conditional-Coulomb when T2/T3 failed falls here).
- **NO-SURPLUS-BUT-PROCEED (P2 accounting):** `EM_COUPLING_ABSORPTIVE_NO_SURPLUS`, `NO_CASCADE_ABSORPTIVE` (the calibration merely
  absorbs `e`/`α`, like Gate-4's `g_G`-absorbs-`G` — report plainly, no held-out surplus claim). `EM_COUPLING_CONSTRAINED_SURPLUS`
  / `EM_COUPLING_PREDICTIVE_SURPLUS` proceed to a *conditional* surplus claim only if the `CONDITIONAL_SURPLUS_ELIGIBLE` gate holds.

## CROSS-CONSISTENCY GAUNTLET (death-mode gates spanning T1–T5; any FAIL is first-class)
- **G1. Lorentz / preferred-`w` frame (tied hard to T1d).** Does a wall in a bulk with a (possibly imposed, per T1d) `w`-axis pick
  a frame? Can brane observers recover effective Lorentz symmetry? **If T1d returned `W_IMPOSED_FOR_STABILITY` (explicit easy
  axis), G1 DEFAULTS to `FAIL_PREFERRED_FRAME` unless an effective-Lorentz-suppression mechanism is *derived* (not assumed)** — an
  imposed easy-axis breaks isotropy at the level of the action and the burden is on showing the breaking is unobservable at brane
  energies. `EFFECTIVE_LORENTZ_OK` / `EFFECTIVE_LORENTZ_SUPPRESSED_CONDITIONAL(scale)` / `FAIL_PREFERRED_FRAME`.
- **G2. Geon stability/lifetime** (MASS = trapped standing wave; gated in T4b; cross-recorded). `GEON_STABLE_OR_LONGLIVED` /
  `FAIL_GEON_TRANSIENT`.
- **G3. Charge sign / antimatter** (both `±w` punctures exist, map to `η_Q=±1`; gated in T4a; cross-recorded).
- **G4. Two-puncture interaction** (correct attraction/repulsion sign + `1/r²`; gated in T4c; cross-recorded).
- **G5. Gravity consistency.** "Our 3D space = the wall" must stay compatible with **gravity = bulk flow into defects**, with **no
  double-counting** of throat energy and no breaking of the drain picture. `GRAVITY_CONSISTENT` / `FAIL_THROAT_ENERGY_DOUBLE_COUNT`
  / `FAIL_DRAIN_PICTURE_BROKEN`.
- **G6. Magnetism / Magnus** (wall shear must not contaminate the shear-free bulk; gated in T3a; cross-recorded).
- **G7. Mixed-core preservation.** The puncture picture must respect that `A_w`, `J^w`, `F_{μw}`, `E_w`, `C_a` are suppressed only
  in the strict far-field reduction, **not erased** from the microscopic ontology. `MIXED_CORE_PRESERVED` / `FAIL_CORE_ONTOLOGY_LOST`.
- **G8. Parent-Maxwell bookkeeping.** State whether brane-tension-as-electric-field-energy **replaces / derives / coexists with**
  the existing localized `A_M` Maxwell sector (`pde.tex:357-365`) — duplicating both **double-counts EM**.
  `MAXWELL_REPLACED_CLEANLY` / `MAXWELL_DERIVED_AS_EFFECTIVE` / `FAIL_EM_DOUBLE_COUNT`.
- **G9. Dark-energy ↔ gravity consistency.** The T5 brane→bulk drainage must be the **same** flow as the G5 gravity drain
  (one medium, one drain — not two book-keepings of the same throat). `DE_DRAIN_IS_GRAVITY_DRAIN` / `FAIL_DRAIN_DOUBLE_COUNT`.
- **G10. (NEW) Wall-energy cosmology / gravitational backreaction.** The wall itself has an energy budget — a stable domain-wall
  *network* or a large `σ_brane` is a classic cosmological hazard (**Zel'dovich–Kobzarev–Okun bound**), independent of (and prior
  to) the T5 dark-energy story. Check: the wall-tension **energy density** and whether it **overcloses** the universe;
  **domain-network formation/coarsening/nucleation** (do a network of `±w` walls survive — and is a *single* observed brane
  compatible?); the wall's **extrinsic-curvature / gravitational backreaction**; and whether `σ_brane` **contaminates the
  dark-energy/dark-matter budget** (i.e. is the T5 "dark energy" actually just the wall's own tension energy, or genuinely the
  bulk↔brane *flow*?).
  - **A-PRIORI CONCERN — the T4↔G10 tension-magnitude contradiction (GLM concern 5; this is a likely CONTRADICTION, not a routine
    check).** Charge calibration (T4c/P2) wants `σ_brane` **large** (the tension-energy around a puncture must reproduce the Coulomb
    self-energy `~q²/a` and the measured `e`/`α` — typically a high scale). Cosmology (G10) wants `σ_brane` **tiny** (a wall that
    dominates a Hubble volume needs `σ ≲ (meV)⁴`-ish). These differ by *many* orders of magnitude. The T5 "dark energy = the *flow*,
    not the static tension" separation must be **DERIVED, not assumed**, to dodge this — and **even then** the wall's own static
    tension energy must be shown not to overclose. **Expected G10 outcome: `FAIL_WALL_OVERCLOSES`** unless a mechanism suppresses
    the wall's gravitational backreaction. Treat a clean pass here with maximum scrutiny.
  Outcomes: `WALL_COSMOLOGY_OK` / `WALL_BACKREACTION_CONDITIONAL` / `FAIL_WALL_OVERCLOSES` / `FAIL_SIGMA_BACKREACTION_FATAL`.

---

## Discipline (binding, same as pathA_23)
Target-blind / able-to-fail at every step; failure is the default expectation; never soften/rescue a FAIL; a clean derivation
re-checked HARDER. **No tautology / no reverse-engineered structure** (T0 is the operational guard — an arrow Lagrangian tuned
*only* to give shear/charge/`c_γ` is circular and caps the rung at CONDITIONAL; the "puzzle = richer-substructure" lens is a
future-directive motivator, never an in-flight rescue). **Native-mechanism discipline** ([[feedback-native-em-mechanisms]]):
charge = puncture **direction** `±w` (NOT winding/π₂/swirl); magnetism = swirl; gravity = drain; light = brane shear; **no point
particles** (defects are extended finite throats — so classical point-charge/point-defect machinery, including vacuum-manifold
homotopy charges, is the wrong tool). Honest provenance + conditional-verdict rule (the polar-OP is a NEW_PARENT_ACTION ⇒ at best
CONDITIONAL-on-the-postulate; the win is *leverage*, not "derived"). Dimensional consistency with units restored (T4-Units and the
T5 rate units are gates, not footnotes). Dual-engine symbolic checks + transliteration-fidelity audit after every computational
step. Codex designs routes + writes/runs all scripts; Claude reviews only. No paper / `decisions/*` edits during execution
(reports → `reports/`; orchestrator owns canonical docs). Scripts ≤10 min with `timeout 600` **on the scripts** (never the Codex
session). **Gate language only** in this directive — vision/payoff framing lives in `docs/conceptual_foundation.md` +
`decisions/15`; nothing is called "derived"/"retired"/"closed" until its gates pass.

## Deliverables
**Execution is REPORTS-ONLY.** Per rung: `reports/pathA_24_T<n>_*.md` (verdict-first) + Mathematica `.wl` (+ SymPy cross-check) +
the T0 freeze artifact (`reports/pathA_24_T0_freeze.md`). The final synthesis report states, with honest provenance, **which of**
{stable wall, emergent-`w`, light/2-polarization shear, shear-free bulk, two charge signs, finite throat, dark-energy cycle} the
single polar-OP postulate actually delivers, and at what provenance grade.
**Canonical-doc edits are OUT of scope for execution and gated behind a SEPARATE integration directive** (`pathA_24_integration`),
to be written and run **only after explicit user acceptance** of the result and its provenance cost: edits to `decisions/15`,
`decisions/14`, `decisions/13` §0, `STATUS.md`, `docs/conceptual_foundation.md` (§7 open-questions → resolved/refined),
`pde.tex`, and any statement that a current gap (`λγ`, the EM anchor) has been retired must wait for that integration step. A
foundational CONDITIONAL result must NOT leak into canonical docs before the user accepts the provenance.

## Review plan
1. **Codex design-review** of the v2 DRAFT (`-c model_reasoning_effort=xhigh`) ✅ round 1 → `NOT-SOUND` (14 fixes) → applied → v2.1.
2. **Codex confirm-pass** on v2.1 ✅ → `NOT-SOUND` (3 routing/hygiene residuals: exhaustive-router by-rule + PASS labels, T2-fail
   severity split, `FAIL_NO_SHEAR` alias + version hygiene) → applied. **Codex reconfirm** ✅ → `NOT-SOUND` (1 residual:
   `FAIL_CHARGE_FIREWALL` mis-bucketed SUBCLAIM-TERMINAL — it is a declared concept-breaker → moved to CONCEPT-FATAL) → applied.
   **Codex closing reconfirm ✅ `SOUND-AS-IS`** (full concept-breaker sweep passed). Directive is Codex-solid.
3. **GLM tertiary** (foundational, user-mediated) ✅ → `SOUND-WITH-CONCERNS` (10 concerns; structurally sound, gates reachable-FAIL).
   Folded into v2.2: the three-way no-win honest-prior block; the MacCullagh-aware T2 reframe (Claude+GLM resolved — GLM's "T2
   pre-decided FAIL" applied a generic liquid-crystal default over the model's MacCullagh rotational-elasticity template, in which
   the photon IS a rotational mode; the honest fork is clean-MacCullagh vs Frank/Cauchy); T0 operational single-medium test +
   framework-circularity caveat; T1c preregistered long-lived threshold; T4a universality≠quantization + symmetry-triviality flag +
   sharpened field-coupling labels; P2 operationalized independence; G10 ZKO σ-magnitude contradiction.
4. **Codex confirm-pass on the GLM-driven v2.2 edits** (verify they landed soundly + the MacCullagh-aware T2 framing is able-to-fail
   in both directions, not a rescue) — NEXT.
5. **Execute rung-by-rung** (T1 FIRST) with per-rung tri-review (orchestrator re-run + transliteration-fidelity audit + adversarial
   review on clean agents) + user gate. **T1 gates everything** — if a polar field can't make a stable wall, T2–T5 are moot and
   that is the result. Honest expected harvest (per the priors block): a CONDITIONAL stable wall (axis likely imposed) + a decisive
   clean-MacCullagh-vs-Frank/Cauchy test of light.
