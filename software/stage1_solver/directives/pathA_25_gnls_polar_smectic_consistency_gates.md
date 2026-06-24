# Directive pathA_25 — The GNLS polar-smectic superfluid: a consistency-gate test for ONE structure that satisfies ALL the requirements

**Status:** DRAFT v4 (Codex design-review xhigh + GLM tertiary + Codex confirm-pass round-1 all folded — see §19 changelog; awaiting
Codex confirm-pass round-2 → user gate before execution).
**Supersedes the *method* of:** `directives/pathA_24_brane_existence_defect_structure.md` (the little-arrows T1–T5 *derivation ladder*).
**Keeps the *content* of:** `reports/pathA_24_T0_freeze.md` (the frozen GNLS + polar-OP action, hash `8fa41ac51e88`).
**Conceptual source (read first):** `docs/conceptual_foundation.md` (v3) and `docs/medium_requirements_and_prior_art.md` (requirements
A/B/C, prior-art survey, the candidate structure, the gates).

---

## 0. Provenance & what changed (why this is a NEW directive, not a pathA_24 rung)

pathA_24 **T1 falsified** the *static* little-arrows domain-wall brane (`T1_FAIL_NO_STABLE_WALL`, tri-reviewed genuine, commit
`2fa91886`): the frozen O(4)-isotropic soft-spin polar action has a **connected** vacuum manifold `S³` (π₀=0), so a `+w`/`−w` wall
**spreads to infinite width** (`σ_L = 5π²Ka²ρ0⁵/(4Lhalf) → 0`) and **unwinds with zero barrier** — no localized wall, no flat core, no
confinement. The φ⁴ negative control returned a *stable* kink (the test was able-to-fail both ways), so the failure is real.

Two user reframes (2026-06-23) drove the pivot, both now binding (memory: `feedback-analog-find-consistent-structure`):
1. **Dynamic, not static.** The medium is intrinsically self-interacting/time-dependent; "lowest-energy *static* configuration" was
   the wrong question.
2. **We are building a mathematical ANALOG, not deriving the universe.** We do **not** have to derive the medium/brane/constituents
   from first principles (they could be arbitrarily deep). The finite, sharper question is: **is there a single self-consistent
   superfluid structure that satisfies ALL the requirements (A/B/C) at once?** We may **postulate the structure freely**; the only
   test is internal **consistency**.

A 5-agent prior-art survey (`docs/medium_requirements_and_prior_art.md`) identified the concrete escape: the **smectic / lamellar**
mechanism makes a codim-1 surface the *spontaneous ground state* (broken translation), qualitatively immune to T1's spread-and-unwind.
That is the structural change this directive tests.

---

## 1. The methodology (BINDING — this replaces the T-ladder)

**Specify the FULL candidate structure (postulated freely — it's an analog), then test internal CONSISTENCY across all requirements
by HUNTING FOR A NO-GO.** A no-go = **two or more requirements are mutually incompatible in this structure** (e.g. "any layering term
strong enough to make a stable codim-1 stripe also destroys the long-wavelength `c_s`/gravity"). This is the new adversarial target.
It is **not** "the minimal postulate wasn't enough" (which just invites endless structure-adding) — that failure mode is retired.

Consequences that are binding for every gate below:
- **The FULL structure — including the light-sector constitutive package — is frozen UP FRONT in G0 (§7), before any gate is
  computed.** Postulating an *ingredient* is allowed; postulating an *outcome* is not, and **adding an ingredient after seeing a gate
  result is `AD_HOC_RESCUE`** (a structural change then requires a fresh G0). We may freely choose the action's terms (the *family* of
  layering drivers, the MacCullagh stiffness + anchoring package). We may **not** hand-draw the layer profile, assert a mode is
  bounded-below, or assume a leak is small — those are the *consistency questions the gates decide*. (Direct T1 lesson: never
  postulate the brane before deriving-or-failing it.)
- **Conditional-verdict rule (inherited from pathA_23/24).** Any gate that *relies on* a postulated ingredient (the driver-family
  choice, the MacCullagh package) returns a `…_CONDITIONAL` pass at best, with the postulate named in the verdict. A clean
  unconditional pass is only available where the result follows from the *kept* GNLS+OP action alone.
- **Every postulated knob is counted (§7 G0.3).** Because we postulate freely, the *only* guard against drifting into a second medium
  is a complete parameter ledger: every independent constant AND function (dimensionful or dimensionless — kernel shape, MacCullagh
  moduli, anchoring inertia/gap, impedance/leak, couplings) is counted. ≥2 independent new inputs ⇒ `SECOND_MEDIUM_DRIFT` pressure
  (NG5), reported plainly.
- **A no-go is a first-class success of this program** (`feedback-falsification-is-the-goal`). Never rescue or soften it. A clean
  "every gate passes" is the *suspicious* outcome and triggers extra adversarial scrutiny.

---

## 2. The candidate structure under test (concrete — this is what makes the gates testable, not generic)

One compressible superfluid `ψ = √ρ e^{iθ}` on `X^i=(x,y,z,w)`, `i=1..4`, whose ground state is postulated to develop **at once**
(i) a 1D density layering (smectic → the emergent-axis brane), (ii) a polar in-plane orientation within the layer (the arrows →
charge), and (iii) an in-plane rotational-elastic (MacCullagh) stiffness on the layer (→ light). Working name: **the GNLS
polar-smectic superfluid**. All three are frozen in G0; the gates test their mutual consistency.

### 2.1 KEPT, byte-for-byte (the existing program rides on this)
- **GNLS medium:** `ψ=√ρe^{iθ}`; quantum pressure `(ℏ²/8mρ)(∇ρ)²`; **single-well** EOS `U(ρ)=(K/4)ρ⁵` → `P(ρ)=Kρ⁵`,
  `c_s²(ρ)=5Kρ⁴/m` (B3); flow `v_i=(ℏ/m)∂_iθ−(q_*/m)A_i`, circulation/vortices → gravity-drain (B1) + magnetism (B2); gauge coupling
  through `v_i`.
- **The T0 polar OP** `L_pol` from `reports/pathA_24_T0_freeze.md` §2.2, **unchanged** (frozen-action SHA-256 `8fa41ac51e88…`):
  `L_pol = ½ mρ a²(D_t^v P)² − ½ mρ c_s²(ρ) a²(∂_j P^i)² − ¼ mρ c_s²(ρ)(|P|²−1)²`, `D_t^v P^i=∂_t P^i+v^j∂_j P^i`, `P^i∈R^4`, polar
  (`P≠−P`), carried by the medium (ρ-weighted, advected), one-constant Frank, O(4)-isotropic, no easy axis.

### 2.2 ADDED ingredient #1 — the layering / smectic driver (finite, target-blind family; exact set frozen in G0)
The local single-well `U(ρ)` **cannot layer on its own** (T1a: its only stationary endpoint is ρ=0). The smectic must come from a
**new finite-range / non-local term** that makes the density (or coupled density-OP) dispersion develop a **finite-`k` (roton)
minimum that softens** → a stripe ground state. Real-matter pedigree: dipolar-BEC supersolids; soft-core / Rydberg-dressed
condensates; spin–orbit-coupled stripe condensates.

To stay target-blind, G0 (§7) freezes a **finite-dimensional admissible set** with a **single pre-committed baseline branch** and a
small number of **named sensitivity branches** (the T0 §4 branch-budget discipline). The candidate families G0 must finitely
parameterize and choose among:
- **Family R (non-local density kernel):** `½∫∫ρ(X)V(|X−X′|)ρ(X′)` with `V` drawn from a **finite-parameter** form (e.g. a fixed
  soft-core/Gaussian-shell functional with ≤2 shape parameters), `Ṽ(k)` having a roton minimum at finite `k*`. Arbitrary free-form
  `V` is **forbidden** (it is an uncounted functional knob).
- **Family L (Lifshitz / Brazovskii gradient):** `−c₁|∇ρ|² + c₂|∇²ρ|²` (`c₁,c₂>0`).
- **Family C (polar↔density coupling):** a **specified covariant form** (the exact admissible operators frozen in G0, e.g. `λ(∇·P)δρ`)
  — not "forms to be fixed later."

**Target-blind selection rule (binding on G0):** admit/parameterize/rank a driver on locality / symmetry / minimality /
single-medium / real-matter-pedigree grounds **only** — never because it helps a downstream gate. A sensitivity branch may **not**
rescue the baseline unless the roll-up *explicitly states the baseline FAILED and a named frozen alternate passed* (no silent
cherry-picking).

### 2.3 ADDED ingredient #2 — the light-sector constitutive package (postulated; frozen in G0)
pathA_23 Stage 2 established that the medium's record does **not** *derive* a clean MacCullagh shear law (`FAIL_UNSPECIFIED_SUBSTRUCTURE`).
Under the analog reframe we therefore **postulate** the light-sector package and freeze it up front (so it cannot be smuggled in after
seeing B4). The package comprises:
- an **in-plane rotational-elastic (MacCullagh) stiffness** on the layer (energy in the curl of the in-plane material displacement,
  `∝μ_br(∇×u)²`), with modulus/moduli `μ_br` — and **either** a scalar-potential analog `φ` / explicit constraint that removes the
  longitudinal sector, **or no such device** (in which case the **C5 gauge obstruction** below is an expected able-to-fail outcome of
  Gate L);
- a **spin / couple-stress sector** with **specified rotational inertia and stiffness** (exact form + provenance frozen in G0).
  *Ingredient, not outcome:* it is specified as a physical sector, **not** as "a reservoir that achieves closure" — whether it
  *achieves* angular-momentum closure is Gate L's to decide. The **preferred single-medium form reuses the polar field `Pⁱ`'s own
  rotational inertia** as the Cosserat micro-rotation (**zero new DOF**, matching `conceptual_foundation §3` where the arrows ARE the
  gyrostat elements); any independent micro-rotation variables are **new DOF counted in G0.4**;
- a **specified coupling operator** tying `Pⁱ` to the in-plane material displacement `u` (exact admissible operators frozen in G0, as
  for Family C in §2.2). *Ingredient, not outcome:* specified as an operator, **not** as "a tie that makes the modes traction-carrying"
  — whether the modes *are* traction-carrying (vs Frank torque-only) is Gate L's to decide.

**The C5 MacCullagh gauge obstruction (named — this program's own verified failure mode, `decisions/15` §11).** A curl-only potential
`½μ_br(∇×u)²` is invariant under `u→u+∇χ`, but the kinetic term `½ρ(∂_t u)²` is **not** invariant for time-dependent `χ`; the EOM then
forces `∂_t²(∇·u)=0` — the longitudinal mode is a **constrained physical zero mode, not a removable gauge artifact** (Maxwell escapes
because the scalar potential `φ` compensates; MacCullagh has no `φ`). G0 must state whether the frozen package carries a
`φ`-analog/constraint that removes this zero mode; if not, C5 is an **expected able-to-fail outcome** of Gate L(a-iii).

All package constants/functions are **counted in the G0.4 ledger** (a package needing ≥2 independent new inputs is
`SECOND_MEDIUM_DRIFT` pressure, NG5). Gate L (§9) tests whether this *frozen* package is internally consistent (traction-not-torque,
no constrained longitudinal zero mode, bounded-below with angular-momentum closure, no leak) on the *actual* B4 layer — it does not
get to assume any of those.

### 2.4 Honest change from the prior picture (state it plainly, do not bury it)
**The density NOW MODULATES.** The earlier little-arrows claim ("the two states live in orientation, not density; `U(ρ)` stays
single-well; we never modulate density") is **refined**: `U(ρ)` stays single-well *locally*, but the new driver gives a finite-`k`
density modulation. "Density doesn't modulate" is no longer true. This is a deliberate, pedigreed structural addition, recorded as
such.

---

## 3. Falsification stance (load-bearing)

The program **wins information either way**:
- **A no-go between gates** (e.g. Gate B4 needs a driver strong enough that Gate B/gravity breaks; or bounded-below light requires a
  spin reservoir that is a second medium; or the smectic layer is in-plane-liquid so light cannot live on it) is the **headline
  result we are hunting** — a genuine structural impossibility, far stronger than "minimal postulate insufficient."
- **A full conditional pass** (all gates pass, MacCullagh/driver postulates named) means the analog *exists* — the bridge is
  demonstrable — and the payoff moves to held-out surplus under calibrate-predict.

Every gate below is constructed to be **able to FAIL**, with explicit FAIL labels and **negative controls that must return the
opposite verdict** so the test is demonstrably not rigged (`feedback-decisive-test-not-tautological`).

---

## 4. Honest priors / expected terminus (calibrate expectations BEFORE executing)

Recorded now so a later result cannot be back-rationalised:
- **Gate B4 (smectic exists):** prior = **plausibly PASS** (the survey's purpose-built escape; dipolar BECs really do stripe) — but a
  real risk of `FAIL_PHASE_SEPARATION` / `FAIL_NOT_CODIM1` / `SMECTIC_CONDITIONAL` (fine-tuned window). The emergent-axis-for-free leg
  remains the historically hardest (Davies–George–Volkas) — a stable *stack* with an emergent normal is the realistic best case;
  "why exactly 3D" is a hoped-for bonus, not a goal.
- **Gate L (light — THE CRUX):** prior = **genuine coin-flip, the most likely no-go.** A plain smectic is liquid in-plane (no shear);
  the rigidity must come entirely from the postulated package, and it must clear **three** independent hurdles at once
  (traction-not-Frank + 2-transverse/no-longitudinal; bounded-below WITH angular-momentum closure; no inter-layer leak). Any one
  failing is the no-go. This is where pathA_23 already saw `FAIL_UNSPECIFIED_SUBSTRUCTURE` and `LEAK_BOUNDED_CONDITIONAL`.
- **Gate B (gravity bundle on a layered background):** prior = **non-trivial risk** — a smectic generically makes `c_s` *anisotropic*
  and adds a layer-displacement phonon; whether the long-wavelength isotropic GR bundle survives is a real question, not a formality.
- **Gate K (cone-lock `c_γ≈c_s`):** prior = **CALIBRATION_GAP almost surely** (`λγ` stays a free input). A *derived* `c_γ=c_s` would
  be a too-clean surprise → extra scrutiny, not celebration; a calibration gap **cannot** count as a derived B11.
- **Inherited walls (NOT gates):** dynamics/`G`/`α` (calibrate-predict), emergent-axis/why-3D, Lorentz/preferred-frame — concede
  honestly; they are universal walls, not our shortfall.
- **Requirement coverage / out-of-scope (explicit, not silent):** B9 (spin) and B12 (light packets) are open frontiers — no gate in
  this directive. **B10 (dark energy / bulk↔brane cycle)** has a specific signature (area growth, decel→accel crossover) but is
  **deferred to a future cosmology-scope directive** — out of scope here. Recorded as a coverage gap, not an omission.

We do **not** bank any pass. Reaching a clean no-go early is a good outcome, not a failure of the directive.

---

## 5. Honest classification (verdict grammar)

Each gate emits a verdict from its own labelled set (below) **plus** a provenance tag: `DERIVED_FROM_KEPT_ACTION` (follows from
GNLS+OP alone) | `CONDITIONAL_ON_<postulate>` (relies on a postulated ingredient: a driver-family choice, the MacCullagh package…) |
`CALIBRATION_GAP` (a quantity the structure does not fix — declared, not hidden; cannot be reported as derived).

Roll-up verdict for the directive is one of:
- `STRUCTURE_CONSISTENT_CONDITIONAL` — all gates pass (incl. a **minimal B8 compatibility** pass, §12 Gate T), postulates named;
- `STRUCTURE_CONSISTENT_EXCEPT_T_DEFERRED` — all gates pass **except** the throat, where **even the minimal T-compat check was not
  run** (the whole throat analysis is deferred off the critical path). *(If T-compat passes and only the finite-radius deep-solve is
  deferred, the roll-up is `STRUCTURE_CONSISTENT_CONDITIONAL`, not this label — see §12 Gate T.)*;
- `NO_GO(<gate/requirement pair>)` — a named mutual incompatibility (the hunted result);
- `INCONCLUSIVE_<gate>` — a gate could not be decided; must be re-scoped, never silently passed.

---

## 6. Gate ordering & rationale (read before the gates)

Logical dependency forces most of the order; risk-priority sets the rest.

1. **G0 — Structure freeze** (precondition; reports-only + dual-engine dim-check). Freezes the FULL structure: kept action + driver
   family (finite) + the light-sector package + the complete parameter ledger.
2. **Gate B4 — does the driver yield a stable codim-1 emergent-axis smectic ground state?** The **direct T1-replacement** and the
   precondition for everything (no layer ⇒ no light, no charge, no throat). Run **first among computational gates** because (a) it is
   where we were just falsified, so it is the highest-value re-test, and (b) it **produces the actual layer profile** (period,
   amplitude, in-plane texture) that Gate L consumes. **We refuse to hand-draw the layer** (T1 lesson) → Gate L cannot precede B4.
   *(Codex design-review confirmed this ordering: running full Gate L on a generic slab would recreate the T1 hand-drawn-layer trap.)*
3. **Gate L — light on the established layer (THE CRUX / most-likely no-go).** Run **immediately after** Gate B4 on Gate B4's own
   output profile. This is the operational reading of "Gate L first": the first *crux*, run as early as its precondition allows. The
   bulk of adversarial budget goes here.
4. **Gate S, Gate B, Gate Q, Gate K, Gate T** — after the crux. Each is able-to-fail and several feed the §13 gauntlet.

**Optional early L precheck (non-verdict, user-approved):** a *profile-agnostic algebraic constitutive* check of the frozen
light-sector package — mode count, energy positivity, and angular-momentum closure for an *arbitrary* positive surface density /
support — MAY be run in parallel with B4. **The angular-momentum-closure item is included only if it is first shown able-to-fail for a
uniform slab** (if closure is an algebraic identity for uniform density it is tautological → drop it from the precheck and keep only
mode count + energy positivity, which remain profile-dependent through the effective moduli). It is explicitly **not a verdict**: it
may **not** test leak, confinement, or charge, may
**not** declare `LIGHT_OK`, and its result is provisional until re-run on B4's real profile in Gate L. It requires explicit user
approval before running (reconciles §16's one-gate-at-a-time rule).

---

## 7. G0 — STRUCTURE FREEZE (anti-circularity preregistration; binding on all gates)

**Goal:** freeze the FULL postulated action — kept GNLS + kept T0 `L_pol` + the layering driver + the light-sector package —
target-blind, *before* any gate is computed, so no gate can be tuned to its answer and no ingredient can be added later.

**Requirements (Codex authors the freeze artifact `reports/pathA_25_G0_freeze.md`):**
- **G0.1** Re-state the kept GNLS action and the T0 `L_pol` **byte-identically**, and **re-verify the T0 SHA-256 `8fa41ac51e88…`**
  (fidelity guard — both engines recompute it before every later gate, exactly as T1 did).
- **G0.2** Freeze the **layering driver**: a **finite-dimensional admissible set** (from §2.2 Families R/L/C), the **single
  pre-committed baseline branch**, the **named sensitivity branches**, the **target-blind admission/ranking criteria**, and the
  **priors** — in the T0 §4 style. No driver may be admitted/ranked/tuned by any downstream payoff. (Reproduce T0 §7's
  forbidden-information statement, extended to: smectic stability, light mode count/bounded-below/traction-vs-torque, angular-momentum
  closure, leak, magnetism, charge signs, `c_γ=c_s`/`λγ`, throat.) **Family-C note (NG2):** the `Pⁱ↔∇ρ` coupling form has a *direct
  NG2 consequence* — a `λ(P·∇ρ)²`-type form tends to **pin `Pⁱ` along the layer normal (out-of-plane)**, depleting the in-plane polar
  order Gate L needs, whereas a `λ(∇·P)δρ`-type form does not pin `Pⁱ`. G0 must **record which admitted forms carry this
  out-of-plane-pinning risk** (target-blind — recorded, not used to choose).
- **G0.3** Freeze the **light-sector package** (§2.3) **target-blind** (same rule as G0.2 — chosen on physical/symmetry/minimality
  grounds, never because it helps Gate L pass): the MacCullagh stiffness form + modulus/moduli; the spin/couple-stress sector's
  **specified inertia + stiffness** (state explicitly whether it **reuses `Pⁱ`** = zero new DOF, or adds **independent micro-rotation**
  = new DOF); the **specified `Pⁱ`↔`u` coupling operator**; and **whether a `φ`-analog/constraint is included** to remove the C5
  longitudinal zero mode (if not, record that C5 is an expected able-to-fail outcome of Gate L). Each with provenance (`POSTULATED`)
  and branch budget.
- **G0.4 COMPLETE PARAMETER LEDGER.** Tabulate **every independent parameter AND function** in the full frozen structure —
  **dimensionful and dimensionless** — including: driver strength + range + roton scale `k*` + any kernel shape parameters; MacCullagh
  moduli; anchoring inertia/gap; any impedance/leak parameter; Family-C couplings. Classify each as `medium-related` (tied to existing
  `ρ,K,m,a,c_s`) vs `independent new input`, and classify **derived** scales (e.g. `k*` emerging from the dispersion) **separately**
  from independent inputs. Carry the T0 §6 rule: **≥2 independent new inputs ⇒ raise `SECOND_MEDIUM_DRIFT` pressure** explicitly (NG5)
  — do not hide it under "one medium" language.
- **G0.5** Produce a **new combined freeze block** (kept action + driver + light package) with its own SHA-256, hashed by the same
  `awk | sha256sum` recipe T0 used. Every later gate re-verifies BOTH hashes (T0's and the new one) before computing.
- **G0.6** **Dual-engine dimensional check** (SymPy + Mathematica, units RESTORED per `feedback-dimensional-consistency-check`): every
  new term reduces to the action density `M L⁻² T⁻²`; `k*` has `L⁻¹`; the MacCullagh modulus and anchoring inertia carry the right
  dimensions. Report the **long-wavelength (`k→0`) limit explicitly**: the driver must be a *finite-`k`* feature that leaves the `k→0`
  `c_s`/EOS available for Gate B to check. Orchestrator re-runs both engines as arbiter + transliteration-fidelity audit.

**G0 outcome labels:** `G0_FROZEN(<newhash>)` | `FAIL_NO_ADMISSIBLE_DRIVER` (no finite-range family meets locality/symmetry/single-
medium without already importing a payoff) | `SECOND_MEDIUM_DRIFT_AT_FREEZE` (the full structure needs ≥2 independent new inputs → the
"one medium" claim is already strained — first-class finding, report it; routes per §14).

**Scope:** reports-only + the dim-check computation. No gate solve. (T0 scoping lesson: a "freeze" rung that makes arithmetic claims —
the dim-check and the `k→0` limit — still gets the full dual-engine + arbiter + fidelity treatment.)

---

## 8. Gate B4 — does the driver yield a STABLE, codim-1, EMERGENT-AXIS smectic ground state? (the T1-replacement; first computational gate)

**The question (able-to-fail):** does the **pre-committed baseline** driver make the GNLS polar-smectic ground state a **stable,
periodic, one-dimensional (codim-1) density modulation** — a *stack of 3-surfaces* with an **emergent** layer normal — rather than
uniform, collapsed, phase-separated, or higher-dimensional (columnar/crystal)?

**What must be established (Codex designs the route; dual-engine):**
- The dispersion of small fluctuations about the uniform state — the **full coupled (ρ, θ, Pⁱ) spectrum** (plus any light-sector
  field contributing at quadratic order in the ground state; **not** the density sector alone, since `L_pol`'s Frank stiffness /
  potential fluctuations can stabilize or destabilize the stripe) — locating the roton minimum `ω²(k*)` and the **softening threshold**
  `ω²(k*)→0`. Show the instability is at **finite `k*`** (smectic), not `k=0` (phase separation).
- Beyond threshold, the **nonlinear stationary modulated state** (amplitude/period selection) and its **stability** (bounded-below; no
  runaway → no collapse). The resulting layers' elasticity must be the **rotationally-invariant de Gennes smectic form** (normal
  *compression* + layer *curvature*, not a frame-dependent `|∂u|²+|∇²u|²`), with both moduli **`B,K>0` verified from the derived
  profile** (not assumed).
- **Codim-1 check:** the selected modulation is a single-`k*` stripe (one emergent axis), not a multi-`k` lattice/crystal. A
  1D-stripe *ansatz* alone is **not** an acceptable pass — the single-`k` state must be shown energetically preferred over the
  competing multi-`k` states. Record whether the axis is genuinely emergent (spontaneous translation+rotation breaking) vs pinned by
  the box/BCs (honesty per Davies–George–Volkas).

**Negative / positive controls (the test must be demonstrably able-to-fail):**
- **NC-B4a [negative]:** the pure kept GNLS (driver coefficient → 0) MUST return `FAIL_NO_MODULATION` (confirms the driver, not an
  artifact, makes the layers).
- **NC-B4b [positive method control]:** a textbook roton-softening kernel in its known supersolid window MUST return a stable
  finite-`k*` stripe (confirms the method can detect a real smectic). *(Labeled a positive control, per design-review — not a negative
  control.)*
- **NC-B4c [negative]:** a purely attractive `k=0` instability MUST return `FAIL_PHASE_SEPARATION` (finite-`k` vs `k=0`
  discrimination).
- **NC-B4d [negative]:** a driver favoring a 2D/3D modulation MUST return `FAIL_NOT_CODIM1` (the codim-1 vs lattice discrimination
  actually fires).

**Outcome labels:** `SMECTIC_GROUND_STATE_STABLE` (codim-1, emergent normal, `B,K>0`, bounded-below, single-`k` preferred over
multi-`k`) | `SMECTIC_CONDITIONAL` (stable only in a recorded fine-tuned window — note the window for the §13 gauntlet) |
`FAIL_NO_MODULATION` | `FAIL_COLLAPSE` (unbounded-below) | `FAIL_PHASE_SEPARATION` (`k=0`) | `FAIL_NOT_CODIM1`
(columnar/crystalline preferred).

**Dependency:** PASS (incl. CONDITIONAL) is required for Gate L, Q, T to be meaningful. A FAIL here is a candidate **NO-GO** for the
whole structure (the smectic escape didn't escape) — first-class; report and pause downstream gates pending user decision (enrich the
driver family only via a **fresh G0** — never in place).

---

## 9. Gate L — LIGHT on the layer (THE CRUX; run on Gate B4's actual profile)

**The question (able-to-fail, three independent hurdles — ALL required):** on the Gate-B4 layer, does the frozen light-sector package
(§2.3) actually carry *light*? Because the MacCullagh curl-only form is *postulated*, "2-transverse/no-longitudinal" is **not** by
itself the test — the genuine able-to-fail content is below.

- **L(a) traction-not-torque, mode structure, AND no constrained longitudinal zero mode.** Three sub-hurdles:
  - **L(a-i)** Are the in-plane modes **traction-carrying material-displacement** waves (MacCullagh light) — the polar order **rigidly
    tied to displacement** so the stiffness does mechanical work as traction — rather than **Frank torque-only orientation** waves that
    *look* transverse but carry no traction (not light; `conceptual_foundation §3`)?
  - **L(a-ii)** Are there exactly **2 transverse, non-dispersive (`ω=c_γ k`) modes and no *massive* longitudinal (Cauchy-type) third
    mode**?
  - **L(a-iii) C5 channel** Is the **longitudinal sector gauge** (no *physical* zero mode) — i.e. the **C5 obstruction**
    (`∂_t²(∇·u)=0` from the curl-only kinetic term, §2.3) is **either resolved by the frozen `φ`-analog/constraint or shown not to fire
    on the B4 profile**? A *constrained* longitudinal zero mode is a FAIL, distinct from L(a-ii)'s massive stray mode.
- **L(b) bounded-below WITH angular-momentum closure.** Is the rotational-elastic energy **positive-definite**, AND does the
  antisymmetric MacCullagh stress achieve **angular-momentum / couple-stress closure** through the frozen anchoring reservoir? (The
  pathA_23 lesson: positive curl energy is necessary but **not sufficient** — without a spin/couple reservoir the gyrostat
  instability returns even with positive stiffness.)
- **L(c) no inter-layer leak.** Does the in-plane shear stay **confined to the layer** — no leading-order traction leak into the
  shear-free inter-layer bulk? Quantify the leak (zero / curvature-localized / impedance-bounded) per the pathA_23 Stage-3/3b
  machinery — do **not** assume it is small.

**Controls (must return the opposite verdict):**
- **NC-La [Frank false-positive]:** the frozen `L_pol` pure-orientation (Frank torque) waves MUST return
  `FAIL_FRANK_TORQUE_NOT_MACCULLAGH_TRACTION` (the traction-vs-torque discriminator actually fires).
- **NC-La2 [C5 zero-mode]:** a curl-only package with **no `φ`-analog and no constraint** MUST return
  `FAIL_C5_LONGITUDINAL_ZERO_MODE` (the gate distinguishes a *constrained physical* longitudinal zero mode from a *gauge-removable* one
  — NC-Lb's massive Cauchy mode is a different channel).
- **NC-Lb [Cauchy]:** a Cauchy-elastic in-plane stiffness MUST return `FAIL_LONGITUDINAL_STRAY` (the stray *massive* third mode is
  seen).
- **NC-Lc [no stiffness]:** a genuinely in-plane-liquid layer (zero rotational stiffness) MUST return `FAIL_NO_INPLANE_STIFFNESS` (the
  "plain smectic is liquid in-plane" failure is seen).
- **NC-Ld [no reservoir]:** a package with **positive curl energy but no spin/couple reservoir** MUST return
  `FAIL_ANGULAR_MOMENTUM_CLOSURE` (the gyrostat detector fires on the sharp obstruction, not just on a reverse-sign ghost).

**Outcome labels:** `LIGHT_OK_CONDITIONAL` (L(a-i,ii,iii)+L(b)+L(c) all pass; CONDITIONAL on the frozen MacCullagh package, named) |
`FAIL_FRANK_TORQUE_NOT_MACCULLAGH_TRACTION` | `FAIL_C5_LONGITUDINAL_ZERO_MODE` | `FAIL_LONGITUDINAL_STRAY` | `FAIL_NO_INPLANE_STIFFNESS`
| `FAIL_NO_INERTIAL_ANCHOR` /
`FAIL_ANGULAR_MOMENTUM_CLOSURE` | and for L(c), one of: `LEAK_NONE` | `LIGHT_FREE_SLIPS_CURVATURE_LOCALIZED` |
`LEAK_BOUNDED_CONDITIONAL(ε_leak)` | `FAIL_LEAK_BREAKS_MAGNUS` (only when the leak is leading-order/unbounded — do **not** force a
curvature-localized or impedance-bounded leak into this label; align with pathA_23 Stage-3b). Any genuine FAIL is a candidate
**NO-GO(B6 ↔ {B4 smectic / B2 magnetism / C1 single-medium})** — report which hurdle failed and why; if L(b) needs an anchoring
reservoir that adds ≥2 independent inputs, that is the **NG5 / C1** no-go (record it).

**Provenance note:** even a full pass is `CONDITIONAL_ON_MacCullagh_package` — pathA_23 Stage-2 established the medium does **not**
*derive* this stiffness; here we test whether the *frozen* package can coexist with B4/B2 without the torque-impostor / negative-energy
/ leak / second-medium pathologies. The win is "consistent coexistence," not "derived light."

---

## 10. Gate S — magnetism preserved (distinct from L(c): the bulk reservoir, not the layer leak)

**The question (able-to-fail):** on the **B4 layered background**, does the **inter-layer bulk remain a clean shear-free superfluid
reservoir** that carries quantized vortices + Magnus (B2)? This is the dual of, but **distinct from**, L(c): L(c) asks "does shear
leak *out* of the layer"; Gate S asks "does the **density modulation itself** or the in-plane light stiffness degrade the bulk
vorticity / Magnus / superfluid reservoir between layers."

**What must be established:** the vortex/Magnus response on the modulated background; whether the layering reduces the inter-layer
superfluid phase stiffness enough to disrupt circulation; whether any L-stiffness bleeds into the bulk.

**Negative control:** a model with (wrongly) bulk-filling stiffness MUST return `FAIL_BULK_SHEAR` (the gate sees Magnus destroyed).

**Outcome labels:** `MAGNETISM_PRESERVED` | `FAIL_BULK_SHEAR` | `FAIL_INTERLAYER_CONNECTIVITY` (modulation cuts the bulk superfluid
reservoir) | `MAGNETISM_CONDITIONAL` (preserved only under the same impedance/`ε_leak` price as L(c) — tie the verdict to Gate L's leak
number).

---

## 11. Gate B — brane ↔ gravity compatibility (the layered background must not break the existing GR bundle)

**The question (able-to-fail):** the driver is a **finite-`k`** feature, but the **layered ground state breaks translation + rotation**
— so the long-wavelength physics is no longer guaranteed isotropic. Does the structure leave **`c_s`, flow, drain (B1/B3)** and the
**already-calibrated GR-quadrupole bundle** (`χ_Q≈0.712`, `P0`, `D0`) intact to stated tolerance?

**What must be established (homogenization on the B4 profile, dual-engine):**
- **scale separation (precondition):** the smectic period `2π/k*` must be **well below the GR-bundle extraction scales** (state the
  ratio); if it is not, homogenization is invalid *regardless* of the anisotropy magnitude → `FAIL_NO_SCALE_SEPARATION`.
- the **long-wavelength acoustic metric / sound *tensor*** on the layered background — is `c_s` still effectively isotropic at the
  scales the GR bundle uses, or does it split into in-plane vs cross-layer speeds beyond tolerance?
- the **new low-energy modes** the smectic introduces (the layer-displacement phonon `u`, and any OP Goldstone) — do they leave the
  `χ_Q`/`P0` extraction inputs unchanged, or does an extra propagating mode contaminate the bundle?
- the **inter-layer superfluid phase-stiffness connectivity** (does the drain/flow still connect across layers as the bundle assumes)?

**Negative control:** a driver with a spurious `k→0` contribution MUST return `FAIL_CS_DISTURBED` (the gate sees a long-wavelength
shift).

**Outcome labels:** `GRAVITY_PRESERVED` | `FAIL_NO_SCALE_SEPARATION` (layer period not ≪ bundle scale → homogenization invalid) |
`FAIL_CS_DISTURBED` (long-wavelength `c_s`/EOS moved beyond tolerance) | `FAIL_ANISOTROPIC_CS`
(`c_s` splits in-plane vs cross-layer beyond tolerance) | `FAIL_EXTRA_SMECTIC_MODE_BREAKS_BUNDLE` (the layer phonon / extra Goldstone
contaminates `χ_Q`/`P0`) | `FAIL_SUPERFLUID_CONNECTIVITY_LOST` (drain/flow no longer connects across layers) | `FAIL_BUNDLE_BROKEN`
(`χ_Q`/`P0` changed → would invalidate Gate-3) | `GRAVITY_CONDITIONAL` (preserved only if the driver strength is bounded below a
recorded threshold — feed the threshold to §13 NG1).

---

## 12. Gate Q, Gate K, Gate T (after the crux)

### Gate Q — two charge signs (B7) — on the B4 profile
**Question:** does the **actual B4 layered profile** supply **two mirror, mass-independent puncture directions** (charge = puncture
direction, not winding; `feedback-native-em-mechanisms`)? Polarity (`P≠−P`) is *necessary but not sufficient* — the gate must use the
B4 profile and confirm a genuine ±-puncture invariant.
**Controls / FAILs:** the headless-director control (`P~−P`) MUST return `FAIL_NO_TWO_SIGNS` (proves `P≠−P` matters); and the gate MUST
`FAIL_POLARITY_ONLY_NOT_CHARGE` if the B4 smectic has only an unoriented layer normal, a uniform in-plane `P` (no flip), no bulk on
both sides of the puncture, or no mass-independent puncture-direction invariant.
**Labels:** `TWO_SIGNS_OK` | `TWO_SIGNS_CONDITIONAL` | `FAIL_NO_TWO_SIGNS` | `FAIL_POLARITY_ONLY_NOT_CHARGE`.

### Gate K — cone-lock `c_γ≈c_s` (B11 / `λγ`)
**Question:** compute the in-plane light speed `c_γ` (from Gate L's rotational modulus + the layer inertia) and compare to the bulk
`c_s` at a **stated tolerance + provenance**. **Prior: CALIBRATION_GAP** — they are set by different moduli and will not be equal
automatically; **a calibration gap CANNOT be reported as a derived B11.**
**Controls / FAILs:** a wrong-dispersion / no-well-defined-`c_γ` package MUST return `FAIL_CONE_INCONSISTENT` (the gate can fail, not
just bookkeep).
**Labels:** `CONE_LOCK_CALIBRATION_GAP` (expected; `λγ` stays a declared free input) | `CONE_LOCK_DERIVED` (`c_γ=c_s` by structure
within tolerance — a too-clean surprise → extra adversarial scrutiny before believing it) | `FAIL_CONE_INCONSISTENT`. Do **not** bank
a derived cone-lock.

### Gate T — throat / mass (B8) — minimal compatibility required for the roll-up
**Question (two tiers):**
- **T-compat (REQUIRED for any structure-consistent roll-up):** is a defect **puncturing the layers** *compatible* with the B4 layered
  ground state and the charge⊥mass decoupling — i.e. does a finite-throat puncture exist as a consistent configuration (tension-closing
  vs holding-open balance is *well-posed*, bulk on both sides, charge=direction / mass=trapped-wave separable)?
- **T-deepsolve (DEFERRABLE off the critical path):** the full finite-radius profile + self-energy magnitude (the throat deep-solve
  was ruled off the critical path before — `feedback-rescope-blocker-vs-downstream-need`); run only if a downstream observable needs
  it.
**Labels:** `THROAT_COMPAT_OK` (T-compat passes; the finite-radius deep-solve may be deferred — this label already carries that
deferral) | `THROAT_CONSISTENT_CONDITIONAL` (compat + a partial deep-solve) | `FAIL_THROAT_INCOMPATIBLE` (no consistent puncture on the
layered ground state — a candidate NO-GO).
**Roll-up rule (resolves the §5 contradiction):** `STRUCTURE_CONSISTENT_CONDITIONAL` requires **at least `THROAT_COMPAT_OK`**. If even
T-compat is not run, the roll-up is `STRUCTURE_CONSISTENT_EXCEPT_T_DEFERRED`, never the unqualified consistent label.

---

## 13. CROSS-CONSISTENCY GAUNTLET (no-gos spanning gates — the real prize)

After each gate, re-evaluate these **mutual-incompatibility** checks. Any firing is a first-class `NO_GO(<pair>)` result:
- **NG1 (B4 ↔ B1):** does the driver strength Gate B4 *needs* for a stable stripe exceed the threshold Gate B allows before
  `c_s`/anisotropy/the extra-mode/connectivity/the GR bundle break (any of §11's FAILs)? → NO-GO (can't have the brane and keep
  gravity).
- **NG2 (B6 ↔ B4):** does the in-plane rotational stiffness Gate L needs require an in-plane order incompatible with the smectic layer
  Gate B4 actually selects? Examples: the layer is in-plane-liquid by construction; the polar order can't tie to displacement on it;
  **or — the most physically likely channel — the smectic density gradient (a Family-C `λ(P·∇ρ)²`-type coupling) pins `Pⁱ` along the
  layer normal (out-of-plane, the standard smectic-A tendency), depleting the in-plane polar order** (a Gate L
  `FAIL_NO_INPLANE_STIFFNESS` driven by this pinning *is* an NG2 firing, not an isolated gate failure). → NO-GO.
- **NG3 (B6 ↔ B2):** is the only way to make light bounded-below + traction-coupled one that **also** leaks shear into the bulk
  (leading-order) and kills Magnus (L(c)/Gate S)? → NO-GO.
- **NG4 (B6-internal trilemma — NOT a mutual no-go by itself):** does defeating the angular-momentum/gyrostat obstruction (L(b)) force
  a stiffness that reintroduces a longitudinal mode or a torque-impostor (L(a))? Label `FAIL_LIGHT_TRILEMMA`. It becomes a true
  **NO-GO** only when its resolution requires an extra spin reservoir / leak / second medium — i.e. when it **escalates to NG3 or
  NG5**; record that bridge explicitly.
- **NG5 (C1 single-medium ↔ everything):** does satisfying the gates require **≥2 independent new inputs** (the complete G0.4 ledger —
  driver shape + MacCullagh moduli + anchoring inertia/gap + impedance/leak + couplings, dimensionful OR dimensionless)? → the
  structure is effectively a second medium (`SECOND_MEDIUM_DRIFT`), undermining C1 — a first-class no-go.

The gauntlet is the **point** of the directive: we are more interested in finding one of NG1–NG5 than in a clean pass.

---

## 14. DEPENDENCY ROUTING (gate-outcome → downstream eligibility)

- G0 `FAIL_NO_ADMISSIBLE_DRIVER` → stop; re-author the freeze (the structure as posed is inadmissible — a finding).
- G0 `SECOND_MEDIUM_DRIFT_AT_FREEZE` → **do not stop**; record the NG5 pressure, proceed (the gates may still find a sharper no-go),
  but the best achievable roll-up is flagged with the drift. The user decides whether the drift alone is disqualifying.
- Gate B4 `FAIL_*` → candidate NO-GO; **do not** run Gate L/Q/T on a non-existent layer (T1 lesson). Report; user decides whether to
  enrich the driver family via a **fresh G0** (never in-place — that is `AD_HOC_RESCUE`).
- Gate B4 PASS/CONDITIONAL → Gate L eligible (on B4's profile); **Gate B, Gate S, Gate Q, and Gate T are also eligible and remain
  meaningful even if Gate L later fails** — they test B1/B2/B7/B8 on the B4 profile *independently of light*, and one of them could be
  the sharper no-go the directive is hunting. (Only **Gate K** genuinely needs Gate L, since it consumes `c_γ`.)
- Gate L `FAIL_*` → candidate crux NO-GO; record **magnetism (S), gravity (B), charge (Q), and throat (T)** as independent structural
  facts; Gate K is moot for *light*.
- Gate L PASS → Gate K eligible (Q and T already eligible after B4).
- Gate S `FAIL_*` → candidate **NO_GO(B2 ↔ {B4 smectic | B6 light package})** — choose the pair from the cause:
  `FAIL_INTERLAYER_CONNECTIVITY` → B2 ↔ B4 (layering cuts the bulk reservoir); `FAIL_BULK_SHEAR` → B2 ↔ B6 (light stiffness bleeds into
  the bulk); report; user decides.
- Gate Q `FAIL_*` → candidate **NO-GO(B7 ↔ B4 structure)**; report; user decides.
- Gate K `FAIL_CONE_INCONSISTENT` → candidate **NO-GO(B11 ↔ B6)**; `CONE_LOCK_CALIBRATION_GAP` is **not** a failure (proceed to §15
  payoff bookkeeping with `λγ` declared free).
- Gate T routing (explicit — a real B8 failure must NEVER be misfiled as a deferral):
  `FAIL_THROAT_INCOMPATIBLE` → candidate **NO_GO(B8 ↔ B4 structure)** (a real no-go, not a deferred roll-up);
  `THROAT_COMPAT_OK` **or** `THROAT_CONSISTENT_CONDITIONAL` → eligible for `STRUCTURE_CONSISTENT_CONDITIONAL` (the finite-radius
  deep-solve may still be deferred — that deferral is already carried by `THROAT_COMPAT_OK`);
  **T-compat not run** → `STRUCTURE_CONSISTENT_EXCEPT_T_DEFERRED`.
- Any §13 NG firing → `NO_GO(<pair>)` is the directive verdict; remaining gates run only if they add independent information.

---

## 15. PAYOFF BOOKKEEPING (gated; calibrate-predict — predictive ONLY after a conditional pass)

No payoff is computed until the structure reaches `STRUCTURE_CONSISTENT_CONDITIONAL`. Then, under calibrate-predict
(`feedback-calibrate-predict-methodology`): `c_γ` (Gate K) feeds `λγ`; the charge/throat structure (Gate Q/T) feeds the EM anchor
(`λγ`-pinning) the verdict equation needs (STATUS §"verdict equation"); held-out surplus (g−2, 5PN, ringdown, multi-defect) rides the
shared *derived* `χ_Q`+`P0/D0` bundle. **Non-negotiable (Gate-4 lesson):** any tension/charge or modulus/speed map must derive
**constant-free** (fewer constants than predictions) or it merely ABSORBS the calibration target (zero surplus) — report plainly,
never dress an absorption as a prediction.

---

## 16. Discipline (BINDING — same spine as pathA_23/24)

- **Codex codes + designs the routes + writes/runs all scripts; iterates until exit 0. Claude REVIEWS only** (audit + verify); the
  orchestrator owns scaffolding (directives/decisions/trackers) and never hand-applies *code* fixes
  (`feedback-claude-reviews-codex-codes`, `feedback-codex-is-fix-applier`).
- **AI prose NEVER establishes a math fact.** Every computed claim (dim-check, dispersion, mode count, traction-vs-torque,
  sign-definiteness, angular-momentum closure, leak number, anisotropy, `c_γ`) is **dual-engine** (SymPy + Mathematica wherever
  Mathematica *can* verify — `feedback-dual-engine-required`); the orchestrator independently **re-runs both engines as arbiter**, runs
  a **transliteration-fidelity audit** (code-vs-equations, `feedback-transliteration-fidelity-audit`), and a **clean adversarial
  agent** per gate (`feedback-review-agents`).
- **A FAIL/no-go verdict gets HARDER scrutiny than a pass** (`feedback-negative-verdict-short-circuit`, `feedback-contest-wall-verdict`):
  truth is in the OUTPUT (JSON/scripts), not the prose; contest a wall with full budget before banking it.
- **Units restored** for every dimensional claim (`feedback-dimensional-consistency-check`). **Freeze-fidelity guard:** both engines
  recompute the T0 hash `8fa41ac51e88` **and** the G0 hash before any gate, and assert the frozen blocks are present unchanged.
- **Reports-only.** Gates write to `reports/pathA_25_*`; scripts to `tools/pathA_25_*`; scratch to `_scratch/`. **No canonical-doc
  edits** (`pde.tex`, `decisions/*`, `docs/*`) without a separate integration directive + explicit user acceptance. **Stage explicit
  paths — never `git add -A`** in the solver dir. Commit only when the user asks; linear to master (no branch).
- **Mathematica:** ≤2 concurrent `math -script` seats. **Codex:** launch backgrounded at `-c model_reasoning_effort=xhigh` (verify the
  log shows `reasoning effort: xhigh`; if lower, TaskStop + relaunch); **never** wrap the codex session in shell `timeout`; the 600s
  cap is only for the audit/derivation scripts Codex runs. **Never `pkill`/`killall` by pattern** (shared box).
- **Sequential, user-gated:** one gate at a time; explicit user gate between gates; never roll forward autonomously
  (`feedback-sequential-audit-chunks`). The only parallel-eligible item is the §6 optional early-L *non-verdict* precheck, and only
  with explicit user approval.
- **Per-gate execution-prompt design-review (the T1 lesson):** before each (expensive) gate run, the gate's execution prompt is itself
  design-reviewed (Codex; + GLM on the crux gates B4/L) — a sound directive does not prevent a spurious run from a loose execution
  prompt.
- **YAML/markdown for all LLM-read/written I/O** (`feedback-no-json-for-llm-io`); JSON only for machine-to-machine engine outputs.

---

## 17. Deliverables

- `reports/pathA_25_G0_freeze.md` — the frozen full structure (kept GNLS + T0 `L_pol` + finite driver set + light-sector package) +
  branch budget + **complete parameter ledger (G0.4)** + new SHA-256; `reports/pathA_25_G0_dimcheck.md` — dual-engine dim-check.
- `reports/pathA_25_gateB4_smectic.md`, `…_gateL_light.md`, `…_gateS_magnetism.md`, `…_gateB_gravity.md`, `…_gateQ_charge.md`,
  `…_gateK_conelock.md`, `…_gateT_throat.md` — one verdict-first report per gate (label on line 1).
- `tools/pathA_25_*_sympy.py` + `tools/pathA_25_*.wl` — dual-engine scripts per computational gate (freeze-fidelity-guarded; both exit
  0; engine-agreement asserted).
- A short rolling `reports/pathA_25_STATUS.md` (verdict ledger across gates + the §13 gauntlet state).

---

## 18. Review plan (the per-gate gauntlet — T1 lesson baked in)

1. **This directive** → Codex design-review (xhigh) ✅ → fold fixes (v2) ✅ → **GLM tertiary** ✅ (SOUND_WITH_CONCERNS; folded → v3) →
   **Codex confirm-pass round-1** ✅ (NOT_SOUND — §14 Gate-T routing blocker; folded → v4) → **Codex confirm-pass round-2** → **user
   gate**.
2. **For EACH gate, before its (expensive) run:** Claude drafts the gate's **execution prompt**; **design-review the execution prompt**
   (Codex, + GLM on crux gates B4/L) — *not just the directive*. Fold fixes → run.
3. **After EACH gate run:** orchestrator arbiter re-run (both engines, identical) + transliteration-fidelity audit + clean adversarial
   agent → verdict → **user gate** before the next gate.
4. **Never alter the calibrated process unilaterally** (`feedback-never-alter-calibrated-process`): if Codex/a step is failing, HALT
   and bring it to the user.

---

## 19. Changelog

- **v1 → v2 (Codex design-review xhigh, folded by orchestrator):** (BLOCKERS) freeze the **light-sector MacCullagh + anchoring
  package up front in G0** (§2.3, §7 G0.3) so it cannot be added after B4; **finite-dimensional driver family** + baseline/sensitivity
  budget, no free-form kernel (§2.2, §7 G0.2); Gate L **traction-vs-Frank-torque** discriminator + control resolves the curl-only
  tautology (§9 L(a), NC-La); Gate L(b) **angular-momentum/spin-couple closure** (not reverse-sign ghost) + no-reservoir control (§9
  L(b), NC-Ld); Gate B **anisotropic-`c_s` / extra-smectic-mode / superfluid-connectivity** checks for the layered background (§11);
  **B8 minimal compatibility required for the roll-up** + new roll-up labels (§5, §12 Gate T). (FIXES) NC-B4b relabeled positive +
  added multi-`k` `FAIL_NOT_CODIM1` control + no 1D-ansatz-only pass (§8); rotationally-invariant de Gennes elasticity with `B,K>0`
  verified (§8); leak labels aligned to pathA_23 (`LEAK_BOUNDED_CONDITIONAL` / `LIGHT_FREE_SLIPS_CURVATURE_LOCALIZED`, §9 L(c)); Gate S
  made distinct (bulk reservoir, not layer leak) (§10); Gate Q must use the B4 profile + `FAIL_POLARITY_ONLY_NOT_CHARGE` (§12); Gate K
  tolerance + controls + "gap ≠ derived B11" (§12); NG4 relabeled intra-B6 trilemma with NG3/NG5 bridge (§13); **complete parameter
  ledger** incl. dimensionless knobs (§7 G0.4, §13 NG5); routing for `SECOND_MEDIUM_DRIFT_AT_FREEZE` + Gate B/S eligible after L fails
  (§14); §6/§16 early-L probe reconciled as a non-verdict user-approved precheck. (CONFIRMED) B4-before-L ordering; kept `L_pol` matches
  T0 + hash; T1 stated faithfully; native terms preserved.
- **v2 → v3 (GLM tertiary review, folded by orchestrator):** (BLOCKER) the **C5 MacCullagh gauge obstruction** — this program's own
  verified failure mode (`decisions/15` §11: curl-only kinetic term → `∂_t²(∇·u)=0`, a *constrained* longitudinal zero mode, not a
  removable gauge artifact) — is now named in §2.3, addressed in G0.3 (frozen package must state whether it carries a `φ`-analog/
  constraint), and tested by new sub-hurdle **L(a-iii)** + control **NC-La2** + label `FAIL_C5_LONGITUDINAL_ZERO_MODE` (distinct from
  NC-Lb's massive Cauchy mode). (FIXES) §2.3/G0.3 reservoir + `Pⁱ↔u` tie reworded **ingredient-not-outcome** + DOF accounting
  (preferred form reuses `Pⁱ`, zero new DOF) (F1); §5/§12 `EXCEPT_T_DEFERRED` contradiction resolved (F2); Gate Q + Gate T eligible
  **after B4** (not gated on L) so a Gate-L fail doesn't lose B7/B8 info (F3); routing added for Gate Q + Gate K failures (F4); NG2's
  most-likely channel (**smectic density-gradient pins `Pⁱ` out-of-plane**, smectic-A tendency) named + Family-C risk recorded in G0.2
  (F5). (CONSIDERS) early-L precheck AM-closure item kept only if able-to-fail on a uniform slab (C1); Gate B4 spectrum is the **full
  coupled (ρ,θ,Pⁱ)** system (C2); Gate B **scale-separation** precondition `2π/k* ≪ bundle scale` + `FAIL_NO_SCALE_SEPARATION` (C3);
  B9/B10/B12 coverage/deferral made explicit in §4 (C4). (CONFIRMED by GLM, native-fit) no standard-physics over-import; finite-`k` on
  a single superfluid coherent; de Gennes form, AM-closure make-or-break, NC-La/NC-Ld discriminators, NG1–NG5 all sound; NG5 expected
  to fire at G0 and §14's don't-stop routing is the right design.
- **v3 → v4 (Codex confirm-pass round-1, folded by orchestrator):** (BLOCKER) §14 Gate T routing split explicitly so a real
  `FAIL_THROAT_INCOMPATIBLE` routes to **`NO_GO(B8 ↔ B4 structure)`** and can **never** be misfiled as `EXCEPT_T_DEFERRED` (the v3
  Gate-T line had lumped a real B8 failure into the deferred roll-up — a contradiction with §5/§12 introduced by the v3 F2/F3 edits);
  the redundant standalone `T_DEEPSOLVE_DEFERRED` label removed (deep-solve deferral is carried by `THROAT_COMPAT_OK`). (FIX) §14 now
  routes **Gate S** failures → candidate `NO_GO(B2 ↔ {B4 | B6})` (cause-dependent pair). Codex confirm-pass otherwise verified all
  v1→v2 and v2→v3 folds correctly applied (incl. the C5 obstruction) and found no payoff smuggling / native-terms violation.
