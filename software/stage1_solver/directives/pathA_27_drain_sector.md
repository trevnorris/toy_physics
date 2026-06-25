# pathA_27 — Drain-Sector Derivation (the "2a" throat-soliton gravity mechanism)

**STATUS: SCOPING / DESIGN — review-complete (full gauntlet SOUND, 2026-06-24); awaiting USER GATE → promotion to execution
directive.** This is the design artifact for track 2a. It has passed the normal gauntlet (Codex design-review r1–r3 → GLM
tertiary → fold 11 → Codex re-confirm r4–r5, all SOUND — see Review log) and the §8 decisions are closed. Next: user gate, then
promotion to an execution directive (concrete acceptance criteria + computation plan). Read `docs/conceptual_foundation.md`
§3 (gravity row), §4 (throat), §5 (bulk↔brane cycle) and `STATUS.md` ⭐⭐ first.

---

## 0. One-line goal

Derive the throat's **steady drain** — the background superfluid flow brane→bulk through a held-open throat — **conditional on
`Δμ ≠ 0`**, and from that single conditional derivation **establish the drain's far-field FORM and a stability UPPER BOUND**:
(a) its far-field **FORM** (structurally braneward `1/r²` — see Payoff A) and (b) a stability **UPPER BOUND** `g_max` (via an
independently-derived throttle). **Out of scope (the deeper question — §8):** proving `Δμ ≠ 0` absolutely, and its exact
magnitude, is the metastability/cosmology sub-question this pass does *not* settle.
- **This pass does NOT** derive the gravitational **coupling strength** (a calibration knob), the **equivalence principle**
  (pathA_21c carries `EP_NOT_DERIVED`), or **GR structure** (geodesics / lensing / metric) — all downstream. It is the drain's
  form + stability bound, not "the gravity sector closed."

The native picture this serves: **gravity = the flow between draining defects.**

## 1. Corrected physical picture (LOAD-BEARING — supersedes the "AC→DC rectification" framing in the notes)

The existing notes (`notes/inner_throat/4d_next_steps.md`) and the pathA_26 interpretation carry a **drift** — that a steady
drain must be obtained by *rectifying* an oscillating (AC) trans-brane current into a DC one. **That framing is retired here.**
The correct picture (user, 2026-06-24):

- **No sloshing.** The trapped shear wave (= mass) is **statically bound**: once the defect forms, the mode is larger than the
  throat neck → below waveguide cutoff → **evanescent → genuinely trapped**. It cannot exit the neck (sub-cutoff) and cannot
  enter the bulk (bulk is shear-free, light *is* shear → no medium for it off the brane). It carries **no net trans-brane
  current** — this is exactly what `J_w = 0 on the exact mode` means, and it is *expected and fine*, not a problem to rectify.
- **The photon is not the gravity.** The trapped wave's only job is **structural** — it holds the throat **open** (the pathA_26
  Phase-A force balance: outward wave pressure vs inward tension + backpressure → equilibrium radius `R*`).
- **Gravity is a separate, steady background flow.** With the throat held open, **ground-state superfluid de-structures**
  (ordered brane → unstructured bulk) and **escapes through the puncture into the bulk**. The **steady inflow field around a
  throat IS its gravitational field**; two throats each draining → the flow between them → **mutual pull = gravity** (§3 table:
  gravity = the medium's inflow/drain toward defects; test bodies are carried inward by the flow).
- **What drives it / what throttles it — two DIFFERENT physics, kept separate (load-bearing; see §3.1/§3.2).**
  - **`Δμ` = the THERMODYNAMIC drive** — the brane↔bulk **ground-state free-energy difference** (= the pathA_26 `μ_drive`). It
    is the *potential* that makes fluid want to convert brane→bulk.
  - **`κ_c` = the KINETIC rate / permeability** — flux per unit drive, set by the **de-structuring conversion BARRIER HEIGHT at
    the throat** (the energy cost of converting ordered brane → unstructured bulk *at the puncture*; an OUTWARD brane→bulk
    barrier). The barrier lives in `κ_c` (rate), **not** in `Δμ` (drive). A higher barrier → smaller `κ_c` → slower leak → the
    natural candidate for *why gravity is weak*.
  - Whether the resulting cap falls below `gcrit` is the open, **able-to-fail** question of §3 (it may not).
  - *Note on `conceptual_foundation.md` §5:* its "admits it as a slow, steady leak" describes the **INWARD bulk→brane
    cosmological return**, not the outward throat drainage throttled here. It is an analogy only — it does **not** supply the
    outward barrier (which sits in `κ_c`, §3.1) nor the drive bound (`Δμ_max` from brane metastability, §3.2).
- **Cross-ref / superseded text.** The pathA_26 report's *Interpretation* section (`reports/pathA_26_derrick.md:122–123`) still
  carries the retired **"AC→DC rectification"** language for the drain. **This §1 supersedes it** — there is no rectification;
  the drain is the steady background de-structuring flow described above.

## 2. Where pathA_26 left it (precise connection points)

- **Conservative soliton EXISTS (Phase A, the solid result):** `a* ≈ 1.8924`, `L* ≈ 1.8121`, Hessian eigs `[0.311, 1.616] > 0`,
  virial residuals ≈ 0. The tested pathA_26 object (dimensionless, `K_eos = 0.002`, `G=c=c_s=1`):
  ```
  E*(a,L) = Ω_fluid,excess + I·ω(a,L) + P_vac·V(a,L) + σ·A(a,L)
  ω(a,L)  = √( c_w²(π²/a² + (π/2)²/L²) + μ_in² )            [fixed-action standing wave]
  ```
  with the fixed-μ sharp-wall depletion contributing `+P(ρ0)·V(a,L)` (μ = U'(ρ0); in the execution slice `P(ρ0)=0.002`,
  `P_vac=0.001`, `σ=0.003`). (`reports/pathA_26_derrick.md:19–45`; `_results.yaml:85–191`.)
  - **`ρ0` = the background medium density** (= `1.0` in the frozen ρ0=1 slice). The `/ρ0` in the drain/`g` formulas below is a
    **no-op in this slice** and only matters once units are restored.
- **The drain entered only as a PARAMETER (Phase C):**
  ```
  slaved flux:        Φ_w   = κ_c·(a³/L)·μ_drive
  nonconservative F:  F_nc  = s·(κ_c·μ_drive)²/ρ0·(2a⁷/L , a⁸/L²)
  open stiffness:     K_open = −dF_nc/dq    (enters the Jacobian as H + K_open)
  dimensionless drain: g     = (κ_c·μ_drive)²/ρ0
  ```
  (`reports/pathA_26_derrick.md:60–66`.)
- **Stability threshold (ESTABLISHED — 2a does NOT re-derive it):** `gcrit_+ = 0.00131` (s=+1), `gcrit_− = 0.00606` (s=−1,
  conservative max). Above `max(gcrit) ≈ 0.006`, `det(H + K_open) < 0` → **static divergence (fold/saddle, not flutter)**; no
  passive damping recovers it (`reports/pathA_26_derrick.md:66,68`).
- **THE GAP 2a closes:** `g` was synthetic (the `g≈89` worst-case is not a particle). **Nobody derived the physical `κ_c` and
  `μ_drive` (=Δμ) for a real defect.** That derivation is the whole of 2a.

## 3. The derivation — ONE task, two payoffs (recast 2026-06-24 per user)

There is no "existence vs stability" split. There is a **single derivation — the throttle-limited drain strength** — and both
questions are read off its output. Why they are not separable: the stability check is `g_max < gcrit`, but the conservative
`min(gcrit₊,gcrit₋)=0.00131` is known and **`g_max = (κ_c·Δμ_max)²/ρ0` is not** — the barrier-set `κ_c` (§3.1) and the
metastability-bounded `Δμ_max` (§3.2) are exactly what this derivation produces. You cannot run the stability check without first deriving the drain
permeability (`κ_c`, §3.1) + drive bound (`Δμ_max`, §3.2). A stability "pass" would be **vacuous** if the drain were actually ~0
(a stable particle that does not gravitate is the concept failing quietly) — but in this pass drain-existence is **weak/out-of-
scope** (`κ_c>0` generically; absolute `Δμ≠0` is §8); existence is addressed *incidentally* by the §3.2 metastability derivation.
The live content is the **stability bound** + the **Mach check**.

**The single task — derive the throttle-limited drain:**
1. **Throat permeability `κ_c` — the KINETIC rate coefficient** (NOT the drive). It is the **non-geometric** part of the
   flux-per-unit-drive, set by the **de-structuring conversion BARRIER HEIGHT** at the throat. The *geometric* factor
   `G_c = a³/L` is already known from the held-open throat `(a*, L*)` (pathA_26); the quantity to DERIVE is `κ_c`. Then
   `Φ_w = κ_c·G_c(a*,L*)·Δμ`. (Drive `Δμ` is a separate, thermodynamic quantity — §3.2.)
   - **Logical-completeness branch (NOT a live falsifier this pass):** `κ_c = 0` would mean an *infinite* barrier (impermeable
     throat). But de-structuring is a phase conversion with **finite** energy cost → `κ_c > 0` **generically**; `κ_c = 0` is a
     **measure-zero limiting case** that fires *only* if the barrier is *provably* infinite. Keep it for completeness, but flag
     it physically likely-inaccessible — do **not** count it as a genuine able-to-fail gate (see Payoff A / §4).
2. **Bound on the drive `Δμ` — from BRANE METASTABILITY (a THERMODYNAMIC argument, NOT the kinetic barrier).** Derive the
   **bound** `Δμ ≤ Δμ_max` from the requirement that **the brane exists / does not spontaneously dissolve**: if the drive
   exceeded the spontaneous-bulk-conversion threshold, the ordered brane would convert wholesale into bulk and there would be no
   brane (no throat, no light). So a brane that exists at all implies `Δμ` sits below that threshold. This is **independent of
   `gcrit`** and **independent of the exact cosmology-coupled magnitude** (we derive the *bound*, not the value — see §8).
   - **The kinetic(`κ_c`)-vs-thermodynamic(`Δμ_max`) split is load-bearing.** The barrier (rate) lives in `κ_c` (§3.1); the
     metastability threshold (drive) lives in `Δμ_max`. They are **different physics**, not two names for one parameter — which
     is *what makes #3's `g_max/gcrit` ordering genuinely contingent* (it is not shared-parameter algebra). The metastability
     derivation also incidentally reveals whether `Δμ ≠ 0` at all (existence), consistent with §8 (bound, not exact magnitude).
   - **⚠ ANTI-TAUTOLOGY FIREWALL (this is what makes the gate non-tautological and able-to-fail) — TWO requirements:**
     - **(a) No `gcrit` smuggling.** `Δμ_max` (and `κ_c`) MUST be derived from the independent metastability/barrier physics
       above. Their derivation MUST NOT use `gcrit`, the Mach horizon, or any pathA_26 stability quantity — otherwise
       `g_max < gcrit` is true by construction. Comparison to `gcrit` happens **ONLY afterward.**
     - **(b) Parametric-contingency check (necessary AND sufficient).** Exhibit BOTH `g_max` and `gcrit` as explicit functions of
       the frozen set `{σ, K_eos, ρ0, P_vac, a*, L*}`, then CHECK that the ratio `g_max/gcrit` is **NOT a parameter-independent
       constant** (no structural near-identity). The ordering `g_max < gcrit` must be **contingent** on the parameter values. If
       `g_max/gcrit` comes out parameter-independent → the stability gate is **tautological** (fails the firewall) even if (a)
       holds.
   - **Able-to-fail:** the bound must be derived, not assumed, and the post-hoc comparison can come out *above* `gcrit`.
3. **Bound on the drain, NOT a pinned value.** A bound on the drive gives an **upper bound** on the drain:
   `g_max = (κ_c·Δμ_max)²/ρ0`. The actual drain is `g_phys = (κ_c·Δμ)²/ρ0` with the actual (unknown) `Δμ ≤ Δμ_max`, hence
   `g_phys ≤ g_max`. (Correspondingly `Φ_w,max = κ_c·G_c(a*,L*)·Δμ_max`.)

**Payoff A — gravity's far-field FORM, conditional on `Δμ ≠ 0` (honest scope: mostly structural this pass):**
- **`Φ_w ≠ 0` (drain-existence).** Since `Φ_w = κ_c·G_c·Δμ`, nonzero drain requires the two **VARIABLE** factors `κ_c` and `Δμ`
  nonzero (the third factor `G_c = a³/L ≈ 3.74` is known-nonzero from `(a*,L*)`). Both zero-routes are **WEAK falsifiers in this
  pass**, so do not over-claim them:
  - `Δμ = 0` (degenerate ground states) → no drain. But the absolute `Δμ ≠ 0` question is **out-of-scope (§8)** and this pass is
    *conditional* on `Δμ ≠ 0` — so this is not a live gate here (the §3.2 metastability derivation incidentally bears on it).
  - `κ_c = 0` → no drain. But `κ_c = 0` is a **measure-zero** infinite-barrier limit (§3.1) → physically likely-inaccessible.
  - **Net (honest):** the drain-EXISTENCE falsifier (`Φ_w = 0`) is WEAK this pass. The genuinely live falsifiers are the
    **stability bound** (`g_max` vs `gcrit`, Payoff B) and the **Mach check**, plus **existence-via-metastability** (§3.2).
- **One-sided-throat assumption (state + check).** `Φ_w ≠ 0` giving a brane-side inflow `v_r ≠ 0` requires the throat to drain
  the **brane one-sidedly** (brane → bulk in one `±w` direction), **not** to conduct **bulk-to-bulk through-flow** (bulk-1 →
  bulk-2 with `v_r = 0` on the brane). The **`±w` one-directional puncture** (charge = direction, §3/§4 of conceptual_foundation)
  supports one-sidedness. **Required check:** confirm the drain actually pulls from the brane (`v_r ≠ 0`), not merely conducts
  bulk-to-bulk.
- **`1/r²` far-field — STRUCTURAL-IN-SCOPE (consistency check, NOT a contingent falsifier).** In *this* model the brane is
  **impermeable except at throats**, so the far-field brane-side inflow is **structurally 3D** → reduced-3D continuity
  `ρ·v_r·4πr² = Φ_w` → `v_r ∝ 1/r²`. The `1/r²` is *just continuity given impermeability*, not a triumph and not a contingent
  gate. **The 4D-leak (`1/r³`) lane is physical ONLY with DISTRIBUTED brane permeability** (the §5 inward leak), which is
  **OUT OF SCOPE** here — so it cannot fire in this pass. (Distributed permeability is the parameter that *would* turn on the 4D
  lane.) This **supersedes pathA_21c's `FORCE_POWER_PROFILE_UNDERDETERMINED` for THIS specific model**, which fixes
  impermeability.
- **Inter-defect attraction — a structural FEATURE + consistency check, NOT a falsifier.** Every throat is a **drain** (a sink),
  and the one-directional `±w` puncture means like-signed throats have `Q1·Q2 > 0` **by construction** → the inter-defect
  far-field force is **universally attractive** (a *feature*: gravity is always attractive), not a contingent able-to-fail
  result. Reusing the pathA_21c far-field extraction here is a **formal consistency check**, not a gate. pathA_21c's sign carries
  `SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE`; **this pass INHERITS that residual and does NOT rely on it** (attraction is
  structural, not read off the residual-laden full sign).

**Payoff B — stability upper bound (the likely-easy half, now non-vacuous):**
- **The test:** `g_max < min(gcrit₊, gcrit₋) = 0.00131` (the conservative threshold) → since `g_phys ≤ g_max`, the soliton sits
  **deep in the stable corner** → the pathA_26 "destabilization" is the **black-hole/transonic regime, not the particle.** The
  looser threshold `0.00606` may be used **only if** the drain's back-reaction sign is shown to be `s=−1`. **Determining `s` is
  part of this derivation.** **Able-to-fail:** if `g_max ≳ gcrit` (at the applicable sign), the physical particle is destabilized
  by its own gravity → concept in serious trouble (and we would then need the exact `Δμ` magnitude to decide).
- **Mach cross-check** (Lane 5): define `v_b` = **the drain INFLOW SPEED AT THE THROAT MOUTH (`r ≈ a*`)** — the fastest point of
  the flow, where transonic onset `v_b → c_s` first occurs (this makes the `v_b(gcrit) ≈ c_s` test well-posed). **DERIVE** the
  map `g → Φ_w → v_b`. Then **TEST** whether `v_b` evaluated at `gcrit` is actually near `c_s` (the analog horizon) — do **not**
  assert the equivalence `g ≪ gcrit ⟺ v_b ≪ c_s`. **Able-to-fail:** `v_b(gcrit)` may NOT be near `c_s`; if the physical drain
  needs `v_b ≥ c_s`, there is no subsonic gravity regime.

## 4. Able-to-fail summary (this design must be able to break the concept — [[feedback-falsification-is-the-goal]])

**The genuinely LIVE (able-to-fail) gates this pass:**
- **Stability bound:** `g_max ≳ gcrit` (at the applicable sign `s`) → upper bound does not clear the threshold → particle
  plausibly destabilized by its own gravity. **The firewall makes this real (two prongs):** (a) `Δμ_max` from brane metastability
  + `κ_c` from the barrier are derived from INDEPENDENT physics (NOT `gcrit` / Mach / any pathA_26 stability quantity), so
  `g_max < gcrit` is not true by construction; AND (b) the parametric-contingency check — if `g_max/gcrit` over the frozen set
  `{σ,K_eos,ρ0,P_vac,a*,L*}` is parameter-independent, the gate is tautological and FAILS the firewall. The kinetic-vs-
  thermodynamic split (§3.1/§3.2) is what keeps this contingent.
- **Mach check:** derived `v_b(gcrit)` (inflow speed at the throat mouth `r≈a*`) not near `c_s`, or the drain requires
  `v_b ≥ c_s` at physical scale → no subsonic gravity regime.
- **Existence-via-metastability (§3.2):** the metastability derivation reveals whether `Δμ ≠ 0` at all; a degenerate result
  (`Δμ = 0`) would mean no drain.

**NOT live falsifiers this pass (honest — do not over-claim):**
- **Drain-existence `Φ_w = 0` is WEAK here:** `Δμ = 0` is out-of-scope (§8; this pass is conditional on `Δμ ≠ 0`) and `κ_c = 0`
  is a measure-zero infinite-barrier limit (§3.1). Kept for logical completeness only.
- **`1/r²` far-field is STRUCTURAL, not a gate:** the brane is impermeable except at throats → 3D inflow → `1/r²` by continuity.
  The 4D-leak (`1/r³`) lane needs DISTRIBUTED permeability (§5) = out-of-scope, so it cannot fire here.
- **Inter-defect attraction is a FEATURE, not a falsifier:** one-directional `±w` drains → `Q1Q2>0` → universal attraction by
  construction. pathA_21c reuse is a formal consistency check; its `SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE` is inherited,
  not relied on.
- **One-sidedness is a required CHECK** (drain pulls from the brane, `v_r≠0`; not bulk-to-bulk through-flow), supported by the
  `±w` one-directional puncture — a consistency check, not an able-to-fail gate.

A clean "it all works" is the suspicious outcome; the live gates above must each be able to fire.

## 5. Structural decisions (carry, don't force — honest scope)

- **Baseline = infinite bulk** (the longstanding model — the neighbor-brane was only ever an *observation* that one could form
  at some distance `w`, not a competing foundation). Conservation is **automatic** via single-medium de-structuring (§5(i)); we
  do **not** need a neighbor brane for the local drain + far-field. 
- **The neighbor-brane-at-distance-`w` is a deferred refinement,** load-bearing only for: (a) the **return/circulation**
  (one-way de-structuring vs steady-state — a *cosmology* question), (b) gravity's **range**, (c) **±w = charge sign**. Flag
  where it would change the answer — notably, recovering `1/r²` at **very long range from a finite slab** needs
  warping/graviton-zero-mode localization, which is **out of scope** for the near/mid-field existence question here.
- **`Δμ`'s sign/magnitude is where existence couples to metastability/cosmology** (is the bulk the true ground state → one-way
  drain?). For this pass we derive the far-field **form** with `Δμ` parameterized; whether `Δμ ≠ 0` and its size is the deep
  sub-question (§8).

## 6. Methodology (reuse + discipline)

- **Analytic far-field + magnitude-as-calibration** ([[feedback-analytic-far-field-beats-noisy-4d]]): derive the flow
  form/sign/power cleanly; leave magnitude a calibration knob; use numerics only for profile/magnitude.
- **Reuse:** pathA_26 functional + Jacobian/Hessian machinery (`tools/pathA_26_derrick_sympy.py`); the `Φ_w = ∫ρ v_w`
  diagnostic; the 4D GNLS solver + equilibrium-scan harness (`notes/inner_throat/`, `scripts/inner_throat/run_equilibrium_scan.py`);
  steady-outflow methods (transonic critical-point / Mach; pseudo-arclength continuation if a fold appears;
  shifted-Newton / pseudo-transient regularization — never `JᵀJ`; open/reservoir BCs) from
  `reports/pathA_throat_solver_literature_synthesis.md`.
- **Frozen spec (non-negotiable):** `n=5`, `P=Kρ⁵`, `c_s²=5Kρ⁴/m`, `h=(5K/4)ρ⁴`, `U=(K/4)ρ⁵`, `K_eos=0.002`, `G=c=c_s=1`.
  Conservative baseline `(a*,L*)=(1.892,1.812)`; virial identities exact; `E_total = m_eff c²` coherent.
- **Dual-engine** (SymPy + Mathematica) on every shared symbolic expression ([[feedback-dual-engine-required]]); **units
  restored** for any dimensional claim ([[feedback-dimensional-consistency-check]]). Codex codes, Claude reviews. Execution
  directive gets the full gauntlet, then tri-review (arbiter rerun + transliteration-fidelity + adversarial), then user gate.

## 7. Out of scope (for this first pass)

- The cosmology **return/circulation** (§5 inward leak), dark-energy crossover.
- **Lepton mass spectrum** (tower falsified; one soliton mass only).
- Full **4D coupled time-dependent** solve (use far-field + reduced (a,L) model first).
- **Charge-sign / ±w** mechanism.
- **Re-deriving `gcrit`** (established by pathA_26).

## 8. Resolved scoping decisions

1. **Existence and stability are ONE derivation** (the throttle-limited drain strength), not a sequence — both are payoffs of §3
   (user, 2026-06-24). The earlier D1/D2 split is retired. The derivation is **conditional on `Δμ ≠ 0`**: it yields the far-field
   FORM and a stability UPPER BOUND; proving `Δμ ≠ 0` absolutely + its magnitude is the deeper out-of-scope question (§8.2).
2. **`Δμ` handling = metastability-bound, not exact magnitude** (user, 2026-06-24). Derive the **brane-metastability** bound
   `Δμ_max` (the brane exists ⇒ `Δμ` below the spontaneous-bulk-conversion threshold — a THERMODYNAMIC drive bound, NOT the
   kinetic/de-structuring barrier, which lives in `κ_c`, §3.1) sufficient for the stability statement and the far-field form;
   leave the exact cosmology-coupled `Δμ` magnitude **out of scope** for this pass. Revisit only if `g_max ≳ gcrit`. **The bound
   carries the §3.2 anti-tautology firewall:** `Δμ_max` (and `κ_c`) are derived from INDEPENDENT physics (metastability for the
   drive, barrier for the rate) and MUST NOT use `gcrit` / the Mach horizon / any pathA_26 stability quantity — `gcrit` enters
   only in the post-hoc comparison, which also carries the §3.2(b) parametric-contingency check.
3. **Minimal model first:** reduce to the analytic sink/far-field model + the existing reduced `(a,L)` energy closure before any
   new 4D PDE solve (§6). Confirmed.
4. **Models are POSTULATED minimal analog models, not derived from the frozen spec (user, 2026-06-24 — "Option A").** The
   brane↔bulk free-energy (thermodynamic) model behind `Δμ`/`Δμ_max` and the kinetic conversion model behind `κ_c` are NOT present
   in the frozen pathA_26 spec (and the brane structure is itself unresolved post-B4/R/C). Per the analog methodology
   ([[feedback-analog-find-consistent-structure]]) we **POSTULATE minimal such models as DECLARED inputs** (calibration is fine —
   we are *finding a setup that works*, not proving) and **TEST whether a consistent working regime EXISTS.** The goal is to
   demonstrate a brane+soliton+drain scenario in which particles can exist; the falsification is a **NO-GO between the joint
   requirements** {nonzero drain, brane-metastable `Δμ<Δμ_max`, stable `g_max<gcrit`, subsonic `v_b<c_s`, structural `1/r²`}.
   `Δμ_max` is derived from the POSTULATED free-energy functional's metastability (NOT from `σ`/tension — that is the kinetic
   `κ_c`); `κ_c` from the POSTULATED kinetic conversion law. Derive the postulated models more thoroughly LATER if the working
   model reveals how. Honest cost: the two postulated models are new declared inputs (NG5 input-count pressure, booked).

These close the open decisions; the full design-review gauntlet is complete (SOUND — see Review log).

---

## Review log
- v0 (2026-06-24): drafted from the corrected physical picture (no-sloshing, drain = background de-structuring flow, photon ≠
  gravity) + the pathA_26 connection-point extraction.
- v0.1 (2026-06-24): user closed §8 — recast §3 to ONE derivation (the physical drain strength) with two payoffs
  (existence + stability); `Δμ` = throttle-bound not exact magnitude; D1/D2 split retired. Ready for Codex design-review (xhigh)
  → GLM → fold.
- Codex design-review r1 (2026-06-24, gpt-5.5, xhigh): VERDICT **NOT_SOUND**. 11-item fix-list: 1 MECHANICAL folded (stale
  "For D1 we derive…" → "For this pass we derive…" in §5); 10 PHYSICS items pending orchestrator sign-off (nonzero-vs-out-of-scope
  tension; `g_max`≠`g_phys`; throttle-bound tautology guard; leading "far below gcrit" wording; gcrit sign-split min vs s=−1;
  `κ_c` vs geometric `G_c=a³/L`; `1/r²` lane-selection check; `v_b`↔`g` map for the Mach cross-check; Phase-A functional shorthand;
  §5 "slow steady leak" is bulk→brane not brane→bulk). Not yet folded — orchestrator to verify before fold.
- Fold r1 (2026-06-24): all 10 PHYSICS items orchestrator-signed-off and FOLDED. §0 goal → conditional-on-`Δμ≠0` (FORM + UPPER
  BOUND; absolute `Δμ` out-of-scope); §1 leading-wording neutralized + §5 "slow steady leak" demoted to inward-return analogy only;
  §2 Phase-A object replaced with tested `E*(a,L)=Ω_fluid,excess+I·ω+P_vac·V+σ·A` (+`P(ρ0)·V` depletion); §3.1 `κ_c`=non-geometric
  permeability with `Φ_w=κ_c·G_c(a*,L*)·Δμ`; §3.2/§4 anti-tautology FIREWALL (`Δμ_max` from INDEPENDENT barrier, no `gcrit`/Mach/
  pathA_26-stability); §3.3 `g_max=(κ_c·Δμ_max)²/ρ0` ≥ `g_phys`; Payoff-A REQUIRED lane-selection (3D vs 4D) gate; Payoff-B test vs
  `min(gcrit)=0.00131` unless `s=−1` (deriving `s` is in scope) + Mach `v_b` defined & `g→Φ_w→v_b` map derived-then-tested; §8.1/8.2
  updated. → Codex confirm round next.
- Codex confirm r2 (2026-06-24, gpt-5.5, xhigh): all 10 prior PHYSICS items confirmed **RESOLVED**. VERDICT **NOT_SOUND** on ONE
  NEW physics item only: §3 Payoff A's "`Φ_w ≠ 0` given `Δμ ≠ 0`" ignores that `Φ_w = κ_c·G_c·Δμ` also vanishes if derived `κ_c=0`
  (or conversion forbidden) → a derived `κ_c=0`/conversion-forbidden must be an EXPLICIT drain-fail (gravity=drain falsified)
  outcome, not assumed away. Pending orchestrator sign-off — not folded.
- Fold r2 (2026-06-24): the new item was orchestrator-signed-off and FOLDED. Drain-nonzero falsification made SYMMETRIC over both
  factors of `Φ_w = κ_c·G_c·Δμ`: §3.1 states `κ_c` may derive to 0 (impermeable / conversion-forbidden) as a first-class verdict;
  §3 Payoff A and §4 add `κ_c = 0` → `Φ_w = 0` → gravity=drain FALSIFIED (regardless of `Δμ`) alongside the existing `Δμ = 0`
  branch. → final Codex confirm round next.
- Codex confirm r3 (2026-06-24, gpt-5.5, xhigh): **VERDICT SOUND; NEW_ITEM RESOLVED; no further issues.** Closes the
  Codex-to-green leg. Next per gauntlet ordering = ONE GLM tertiary pass → fold any corrections → Codex back-to-green if needed.
- GLM tertiary r1 (2026-06-24, glm-5.2, session `ses_1045da67affelb4DPNru2q7L3v`): the orchestrator's initial `opencode export`
  salvage hit the known emission hang (0-token final message); a subsequent salvage recovered the verdict: **SOUND_WITH_CONCERNS,
  11 findings** (7 PHYSICS + 4 MECHANICAL). All 11 orchestrator-signed-off.
- Fold r3 (2026-06-24): all 11 GLM findings FOLDED per orchestrator resolutions. PHYSICS: (1) inter-defect attraction reframed
  from falsifier → structural FEATURE + consistency check (one-directional ±w drains ⇒ Q1Q2>0 ⇒ always attractive; pathA_21c sign
  residual inherited not relied on); (2) `1/r²` lane downgraded to STRUCTURAL-in-scope (brane impermeable except at throats ⇒ 3D;
  4D-leak needs distributed permeability = §5 out-of-scope); (3) firewall STRENGTHENED — require g_max,gcrit as functions of the
  frozen set + check g_max/gcrit is NOT parameter-independent; (4) κ_c=KINETIC(barrier-height)/Δμ=THERMODYNAMIC(drive) split, with
  Δμ_max from BRANE METASTABILITY (not the kinetic barrier); (5) κ_c=0 falsifier downgraded to measure-zero logical-completeness +
  acknowledge drain-existence falsifier is WEAK this pass; (6) one-sided-throat (drains brane, v_r≠0; not bulk-to-bulk through-flow)
  stated + checked; (7) §0 softened "closes gravity sector"→"establishes far-field form + stability bound" (no coupling/EP/GR).
  MECHANICAL: (8) "both VARIABLE factors" + G_c≈3.74 nonzero; (9) ρ0 defined; (10) v_b = inflow speed at throat mouth r≈a*;
  (11) cross-ref pathA_26 Interpretation:122–123 still carries retired AC→DC, §1 supersedes. → Codex re-confirm round next.
- Codex re-confirm r4 (2026-06-24, gpt-5.5, xhigh): all 11 GLM folds confirmed **OK except #4 PROBLEM** — two STALE spots
  (§8.2 + §3-opening) still coupled `Δμ_max` to the "tension/barrier," contradicting the item-4 metastability reassignment.
  Scrutiny (a) split physically right, (b) downgrades honestly stated, (c) connection points consistent. VERDICT NOT_SOUND on
  that lone stale-wording contradiction. **Folded** as completion of the already-signed-off item-4 resolution: §8.2 → `Δμ_max`
  from brane-metastability (thermodynamic drive bound, NOT the kinetic barrier, which is `κ_c`); §3-opening "independently-
  barrier-bounded Δμ_max" → "metastability-bounded `Δμ_max` / barrier-set `κ_c`." → final Codex re-confirm round next.
- Codex re-confirm r5 (2026-06-24, gpt-5.5, xhigh, narrow scope): **VERDICT SOUND; R4_ITEM CLOSED; no further issues.** The
  item-4 completion-fix closes the r4 NOT_SOUND; `κ_c`=kinetic-barrier / `Δμ_max`=thermodynamic-metastability is now consistent
  across §3-opening, §3.1, §3.2, §4, §8.2. **This CLOSES the full review gauntlet** (Codex-to-green r1–r3 → GLM tertiary
  SOUND_WITH_CONCERNS [user-run] → fold 11 → Codex r4 NOT_SOUND [stale item-4 spots] → fold → Codex r5 SOUND). The scoping
  directive is review-complete and ready for the user gate / promotion to an execution directive.
