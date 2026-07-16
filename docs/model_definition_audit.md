# Model definition audit — what is actually DEFINED vs IMPORTED vs OPEN

**Purpose.** A deliberately minimal, honest inventory. For every phenomenon: the model's native mechanism, and whether that mechanism is actually *defined* — or whether we've been *importing* standard-physics structure because the native piece is under-specified. The aspiration is a complete brane+bulk equation set we could simulate and read values off of; we are reverse-engineering it from observation, so this doc states exactly which pieces of those equations exist and which don't. This **supersedes** `conceptual_foundation.md` for the question "what do we actually have." **v2 — fact-checked by Codex + Grok (both `AUDIT_NEEDS_CORRECTIONS`, convergent); their label corrections folded.**

**Status vocabulary:**
- **NATIVE-DEFINED** — derived from the medium's own dynamics; an actual equation produces it.
- **CONDITIONALLY-DERIVED** — derived, but *given* a postulated stiffness / ansatz / source law (the conditioning piece is not itself derived).
- **CALIBRATED** — a free constant fit to observation (not first-principles).
- **POSTULATED** — imposed as a model input; internally consistent but not derived from the one medium.
- **⚠️ IMPORTED** — a value/sign/source-law borrowed from standard physics because the native piece is under-defined. (The habit we are hunting.)
- **OPEN** — not defined at all.

---

## A. Foundational objects

| Object | Native mechanism | Status | Defined / missing |
|---|---|---|---|
| **The medium** | compressible superfluid in 4+1D | **NATIVE-DEFINED (compression) + CALIBRATED** | Density/compression EOS specified (quintic, `n=5`, `c_s²=5Kρ⁴/m`). **Missing:** full nonlinear 4D solve (sim-deferred); the sub-mean-field constituents are inferred, not defined. |
| **The brane** (our 3D space) | ordered phase / domain wall at `w≈0` | **POSTULATED** | Wall put in by hand (postulated `χ_B` double-well). **Missing:** derivation from the one medium (route-c) OPEN; **slab width/stability not self-selected** (an unearned input). |
| **The throat** (a particle) | a puncture where brane drains into bulk | **PARTIALLY DEFINED / OPEN** | "Finite puncture" is the picture. **Missing:** size/aspect `L/a` *not* self-selected; the throat's shared **drain-vs-charge energy budget** OPEN; trapped-wave mass *spectrum* **falsified** (1:9:25 vs 1:207:3477). |

## B. The four forces

| Force | Native mechanism | Status | Defined / missing |
|---|---|---|---|
| **Gravity** | throat *drain* (inflow `v_r`); `1/r²` via a finite-slab zero mode; changes ride `c_s` | **NATIVE-DEFINED (form) + CALIBRATED (`G`)** | Best-defined sector. Drain-flow field *is* the gravity field; PN ladder 1PN→4PN matches GR (calibrated). `p=2` derived **within the postulated slab family** (`pathA_29`). **Missing:** nonlinear bulk↔brane **return closure**; `G` fit, not derived. |
| **Light** | brane **MacCullagh** (rotational-elastic) shear; 2 transverse, `c_γ` | **MODE CONDITIONALLY-DERIVED (on postulated stiffness) / stiffness POSTULATED** | The 2-photon mode `ω²=(μ_R/ρ_br)k²` is derived *given* the postulated MacCullagh stiffness `μ_R` (`pathA_36`). **Missing:** origin of that stiffness (structured wall / substructure) OPEN (route-c; `pathA_35 FAIL_COUPLE_STRESS_NOGO`); **light packet localization + quantization** OPEN. |
| **Electric charge** (static force) | two throats' 4D bodies beyond the mouth; `±w` = sign | **FALLOFF NATIVE-DEFINED (cond.) / mediator = a SCALAR / SIGN ⚠️ IMPORTED / MAG CALIBRATED** | Falloff `1/r²` (`p=2`) is derived (conditional on the scalar-`h` mediator ansatz; `pathA_38`). **The mediator is a scalar `h`, not Maxwell.** **Sign is not native:** "like-repel" holds only if a coupling `G₀>0`; the no-lock scalar limit *attracts*. `Q_E` + the nonzero-monopole source CALIBRATED. |
| **Magnetism** (moving-charge force) | the moving 4D throat-body under motion; `±w` | **SOURCE LAW ⚠️ IMPORTED / kernel CONDITIONALLY-DERIVED / native dynamics OPEN** | Given the **declared** source `j = sV`, the `R^-2` current–current force + like-current attraction *are* derived, and the kernel sign is healthy (reversing it needs a ghost sign). **What's imported:** the moving-source law `j=sV` itself. **What's OPEN:** deriving that source from the throat-body's *actual* motion; and under motion **EM ↔ gravitomagnetic operators mix** (`pathA_39` Stage 3). Object is *committed* (moving 4D throat-body, **not** bulk vorticity = gravitomagnetism); its native dynamics were never defined. |

## C. Cross-cutting — defined vs open

| Item | Status | Detail |
|---|---|---|
| **No native first-class `U(1)` / Gauss law** | **NATIVE-DEFINED result: `NATIVE_P_NO_EMERGENT_GAUSS`** | Computed + verified: the continuum little-arrows sector is *second-class* → **no emergent Gauss constraint natively.** Exact-Maxwell / ice-rule EM is a **postulated bolt-on**, not native. |
| **Charge quantization + two signs** | **`Z₂` orientation NATIVE / additive-integer charge + Gauss + conservation OPEN / `Q_E` CALIBRATED** | Bare `±w` gives only *two orientations*. It does **not** supply additive integer charge, a Gauss law, or a conserved current — those are OPEN (per the native-`U(1)` result above). |
| **The E/M dual sign (Ward tie)** | **⚠️ NOT native — TWO independent kernels** | Electric rides the scalar `h`; magnetism rides light's `u_T` — two independently-signed channels, not one Ward-tied Maxwell kernel. Maxwell's dual sign (like charges repel, like currents attract) is not natively produced. |
| **The charge ↔ medium coupling** (what the `±w` body *is* dynamically; how it moves through the medium) | **OPEN** | One of the several EM gaps (not "the" single gap — see §D). Sets the imported source law `j=sV`. |
| **Characterized EM scalar departure** | **NATIVE-DEFINED (as a departure) / magnitude SIM-GATED** | The EM field carries a propagating **charge-coupled `h`-branon scalar** (a 5th-force-like admixture; `pathA_39/42`); its magnitude, universality, mixing, preferred-frame effects are sim-gated. |
| **Mass** (mechanism + value + spectrum) | **mechanism POSTULATED / absolute value OPEN / spectrum FALSIFIED** | Trapped-geon mechanism postulated; the absolute mass scale is *undetermined* (a candidate relation, not a demonstrated value); the lepton mass tower is falsified. |
| **Throat motion / inertia** (core inertia, Berry) | **OPEN** | How a throat responds to a force / moves is undefined — part of why the magnetism source law had to be imported. |
| **Spin** | **OPEN** | Not placed in the picture. |
| **Gravitomagnetism** (mass-flow swirl) | **POSTULATED form + CALIBRATED strength** | The far-zone swirl was imposed by matching Kerr and its strength fixed by the same match — not a native derivation. Source-level distinct from EM-magnetism, but operator mixing under motion is open. |
| **`λγ` (light speed vs gravity-ripple speed)** | **CALIBRATED (`pathA_40 CONE_LOCK`) / native derivation OPEN** | Both `c_γ=c_s` and `c_E=c_γ` are treated as independent *calibrated* locks; the values are set, the derivation is open. |
| **Cosmology** (bulk↔brane dark-energy signs/rates; wall-tension) | **OPEN** | Hypothesis-level; signs and rates unconfirmed. |

---

## D. The finding this audit crystallizes (corrected)

Two sectors are genuinely built:
- **Compression → gravity** (native-defined form, calibrated `G`).
- **Brane shear → light** (mode derived *on* a postulated stiffness).

**All EM under-definition lives in ONE SECTOR — the charge / `±w`-body sector — but it is SEVERAL INDEPENDENT missing pieces, not one "lock":**
1. **No native first-class `U(1)`/Gauss** (`NATIVE_P_NO_EMERGENT_GAUSS`, computed) → no native additive charge / quantization / conservation.
2. **The moving-source dynamics of the `±w` body are undefined** → the magnetism source law `j=sV` was imported; the electric sign rests on an undefined lock.
3. **Two independently-signed kernels** (electric on `h`, magnetism on `u_T`) → Maxwell's dual sign is not natively tied.
4. **Throat inertia, the scalar-departure parameters, and the dynamical throat↔charge binding** — each independently open.

So the honest statement is **not** "one missing definition." It is: *the charge sector is where the model is under-defined, and it is a **cluster of independent gaps**; each imported EM specific traces to one of them.* Whenever an EM value was needed, standard EM got imported at whichever of these gaps was load-bearing.

**Consequence for how we proceed:** there is no single silver-bullet definition. The tractable order is probably: **(1)** decide whether the target is *exact Maxwell* (needs a bolt-on constrained DOF — computed non-native) or the *characterized-departure* analog (largely already computed); then **(2)** for whichever target, define the `±w`-body's own dynamics (motion/inertia) so the source law and signs stop being imported.

---

## E. Explicitly OPEN (the honest list)
1. Native first-class `U(1)` / additive charge / Gauss / conservation (computed *absent* on the continuum; a bolt-on if wanted).
2. The `±w`-body's native dynamics (motion, inertia, the source law behind `j=sV`).
3. The electric sign (rests on the undefined lock; no-lock limit attracts).
4. A single Ward-tied E/M kernel (currently two).
5. The brane's derivation (route-c) / light's stiffness origin.
6. Throat size/structure `L/a`; slab-width selection/stability; bulk↔brane return closure; the drain-vs-charge energy budget.
7. Throat motion/inertia; spin.
8. The lepton mass spectrum (falsified) and absolute mass scale (undetermined).
9. Native derivation of `λγ` (values calibrated).
10. Light packet localization + quantization.
11. Cosmology signs/rates (dark energy, wall tension).

---

## F. The EM-analog "to-define" list — the next research phase

**Target (decided):** the *characterized-departure* Maxwell analog — NOT exact Maxwell (computed non-native via the native-`P` gate). Method: **define each piece the model's own way, compute, and accept whatever departure falls out — do not import Maxwell, do not chase an exact match.**

The Maxwell analog rests on **three IMPORTED pieces, downstream of two definable foundations + one structural fact:**

**⚠️ IMPORTED (to be REPLACED by computed outputs):**
- **I1 — the moving-charge current law `j = sV`** (borrowed convection current; `pathA_39`).
- **I2 — the electric sign "like charges repel"** (holds only via the assumed lock `G₀>0`; the model's no-lock limit *attracts*; `pathA_38`).
- **I3 — the Maxwell exchange machinery + dual sign** (positive propagators, `U=−jGj`; two independently-signed kernels — electric on `h`, magnetism on `u_T` — imposed to a Ward tie).

**FOUNDATIONS to DEFINE natively (the imports are downstream of these):**
- **U1 ⭐ (keystone) — the `±w` throat-body as a dynamical object:** its structure (a puncture extended into `w`), its **effective inertia** (mass from the medium it displaces as it moves), and its **force-response**.
- **U2 — its medium-coupling:** the **roll-vs-slip boundary condition** of a moving puncture in the superfluid — does the medium stick (no-slip) or slip? — derived from the puncture's actual (drain/vortex-like) structure, not chosen.
- **U3 (structural fact, not upgradable) — charge is natively `Z₂`** (two signs); additive-integer charge + Gauss + conservation are computed *non-native* (native-`P`). If additivity matters, it's a separate, non-native question.

**Dependency + order:** define **U1** first (keystone) → **U2** → then the imports become **computed**: the **static** two-body force → the electric sign (replaces I2); the **moving** two-body force → the magnetism sign *and* whether the current is really `∝V` (replaces I1, I3). Whatever the signs and current-law come out to — Maxwell-like or a departure — that is the answer.

---

## G. The gravity↔EM force-strength hierarchy (the `~10⁴²`) — refines E.6

**The fact.** `F_electric/F_gravity ≈ 4×10⁴²` for two electrons (mass-dependent: `~10³⁶` for protons). This is the largest coupling gap in physics and is **unexplained in ALL of mainstream physics** (the "hierarchy problem" family) — not a defect unique to this model. In standard terms it is `α_G = (m/m_Planck)²`, i.e. *"the electron mass is absurdly tiny vs the Planck mass."* In native terms it maps to: **a throat's mass-flux coupling (gravity) is `~10⁴²` weaker than its topological-charge coupling (EM).**

**Conceptual resolution of the "conservation should tie them" worry (do NOT relapse into this).** A tempting error: *"if energy is conserved, the force into the throat = the force out, so gravity and EM can't differ by `10⁴²`."* This is wrong on two counts: (i) force ≠ energy, and (ii) **energy conservation runs WITHIN each force channel, not ACROSS channels** — there is no ledger forcing the gravity-channel energy to equal the EM-channel energy. The deeper resolution is that **gravity and EM are NOT two outputs of one drain**: gravity is sourced by the throat's **mass-flux** (drain throughput); EM by the throat's **`±w` topological charge** (puncture structure). *Different source properties of the same throat* → independent couplings → a `10⁴²` gap violates nothing. (A passive drain has one throughput; EM is not that throughput, it is the topology.)

**⭐ Structural parallel worth reproducing (`(v/c)²`).** Gravitomagnetism is weaker than gravity by `~(v/c)²`, the **identical** structure to magnetism being `~(v/c)²` weaker than electricity (measured: Gravity Probe B / LAGEOS frame-dragging, `~10⁻¹²` for Earth). This "moving source → velocity-suppressed second force" pattern appears to be **generic to how a source couples to a medium**, and the model should reproduce it *cleanly in both sectors* — a consistency check (mild positive test). **BUT** `(v/c)²` is a MODEST suppression (`~10⁻¹²` to `~10⁻²⁰`); it is emphatically **NOT** the mechanism for the `10⁴²` primary-channel gap.

**What `10⁴²` actually requires (and the trap).** The primary-channel gap needs the inflow/drain coupling to be intrinsically `~10⁴²` softer than the charge coupling — a *dramatic* suppression (exponential `e^{−S}`, or a high power of a length-scale ratio), not a velocity factor. ⚠️ **NUMEROLOGY TRAP:** matching a big number by eye (e.g. `(ℓ_P/r_e)² ~ 10⁻⁴⁰`) is worthless (the Dirac large-number graveyard); the model must *derive the exponent and the relationship*, or it does not count.

**Status + why it matters.** Because **one medium sources both forces from the same throat solution**, the ratio is **NOT two free knobs** — it is a *constraint the throat structure must satisfy*. That makes `F_e/F_g` one of the sharpest possible **held-out, both-sectors-at-once dimensionless tests**: derive the mass-flux coupling and the charge coupling from the *same* U1/U2 throat without tuning either, read off the ratio → it lands near `10⁴²` (spectacular win) or it doesn't (clean falsification). **Currently it is FIT, not tested** (`G` calibrated + EM couplings calibrated/imported → the ratio falls out of two independent tunings → no teeth yet). It becomes a real test only once both couplings come from one un-tuned throat solution. This is the operational face of open item **E.6** (the drain-vs-charge energy budget), judged under the dimensionless-held-out-residual rule.
