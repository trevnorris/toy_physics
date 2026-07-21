# Electric sector: the exclusion cascade and the leftover-scalar-as-charge hypothesis

> ⚠ **CORRECTION (2026-07-20, later the same day — after a reevaluation against the maintained EM-track docs). This writeup OVER-PROMISED the `u_L` clamp as "the leading electric direction."** Accurate status: the clamp is a **live CANDIDATE for the still-OPEN electric SIGN** (the audit's I2 "lock `G₀>0`") — ONE of the 144 UNRESOLVED boundary branches that U1 + U2 (both A9-verified) proved the committed model does *not* determine — **not** a native determination. The "field-exchange exclusions" below were EXPECTED, not a discovery: the committed electric mechanism is the `±w` 4D-body + scalar-`h` mediator (the `h`-branon is the committed charge-coupled scalar — the electric field's flex into `w` — **distinct** from `u_L`), not field-exchange. **READ-FIRST the maintained map — `docs/em_analog_next_phase_handoff.md` + `docs/model_definition_audit.md` (§F) + memory [[project-em-sector-reconsideration]] — NOT `conceptual_foundation.md`.** The exclusion cascade and the clamp derivation recorded below remain valid; only the "leading direction" framing is corrected. **User decision (2026-07-20): build the mechanism properly via the full gauntlet.**

**Date:** 2026-07-20. **Status:** synthesis of a full session's EM-track work. **Scope:** how the electric (Coulomb like-repel) sign was hunted, every mechanism that was *excluded* and why, the physics lessons, and the surviving leading hypothesis (electric charge = a stored unit of the leftover longitudinal scalar). Private-chat framing (engages the unification directly); the honest toy-analog scope applies to any public write-up.

Companion reading: `docs/conceptual_foundation.md` (§3–§4, the native picture), `docs/two_throat_simulation_handoff_spec.md` (the `𝔅`/mouth-ensemble fork, force decomposition), `docs/em_phaseC_force_decomposition.md`, `software/em_charge_attribute/g0_closure_card_v0*` (the excluded scalar-mouth mechanism + its committed static-prerequisite check suite `bdd0104d`).

---

## 0. TL;DR

- **Every FIELD-EXCHANGE mechanism for the electric force is excluded** — each gives like-*attract* (wrong Coulomb sign) or fails to close, for reasons that are spin-statistics / structural, not tuning accidents. Triple-checked (Codex + Grok compute + GLM + orchestrator analysis), and several were already documented in the ledger as settled fails.
- **The electric SIGN is not a prediction of the committed model.** It is controlled by the throat's boundary condition (`𝔅` / mouth ensemble), which U2 already proved the committed model does not determine, and which defaults to the wrong sign wherever the model computes anything.
- **The surviving leading hypothesis (user-originated):** electric charge = a **stored fixed unit of the leftover longitudinal scalar** `u_L` (the `FAIL_CAUCHY_STRAY_LONGITUDINAL` mode that refuses to close into Maxwell). If the throat **clamps** (stores) a signed unit of `u_L` rather than **sourcing** it, the stored gradient energy gives **like-repel at `1/R²` = Coulomb**, with **no gauge/`A⁰` smuggling** — both Codex and Grok derived this explicitly. The clamp is physically motivated because a **geon is a standing, trapped wave = stored energy** (a store, not an emitter).
- **The one obstruction that looked fatal (`u_L` is the "mass channel" → equivalence-principle conflict) was debunked:** `q_M → u_L` is a *retrofitted label of convenience* (a hand-placed free symbol `J_mass=[0,0,q_M,0]`, no formula/provenance), not a derived identification. Real mass lives in the throat/geon → **drain** (`pathA_29`), a bulk flow — a different object from `u_L`. The model itself states **"charge ⊥ mass."** So repurposing `u_L` for charge does *not* conflict with mass.

---

## 1. The exclusion cascade (settled)

The electric force is the **charge-odd** (`s₁s₂`-dependent, `s=±1` = which `±w` side the throat threads into) part of the two-throat force. Coulomb requires **like-repel / opposite-attract**. Every candidate:

| candidate mechanism | far-field | sign | why excluded |
|---|---|---|---|
| **compression scalar `h`** (G0 mouth source) | `1/R²` | like-**attract** | spin-0 exchange attracts like sources; `K_eff=K_h−C_hu²/B_eff=7/8>0`; reversing `J_m` can't help (pair energy ∝ source²). Triple-confirmed. |
| **bulk drain** (active flow) | `1/R³` even; `a²/R⁵` odd | attract | 4-D bulk object → subleading to `1/R²`; charge-odd part reinforces attract; responsive gate `Γ(h,P)` → `1/R⁴`, still subleading. |
| **transverse light-vector `u_T`** | — | no monopole; magnetic like-**attract** | divergence-free (`∮u_T·dA=0`) → carries no monopole → no central `1/R²`; where it mediates statically, pathA_39 gives magnetostatic like-attract. |
| **longitudinal/Gauss Maxwell completion** | would give Coulomb | — | doesn't close with provenance: `FAIL_CAUCHY_STRAY_LONGITUDINAL` (`B_eff>0`, `K_θ<0`); Maxwell only `BY_TUNING`; `conceptual_foundation §3` calls "curl-of-light closes Maxwell" a *predetermined fail*. |
| **naive wall-tension / topological defect** | short-range | right *relative* sign, but | no protected `±w` charge: `pathA_24` found connected `S³` vacua (`π₀=0`), charge-sign wall unwinds with zero barrier (`T1_FAIL_NO_STABLE_WALL`); `Z₂` gives walls not point-charges; the repulsive kink channel is short-range Yukawa. |

**Also important — pathA_38 was a false positive:** its encouraging `U₊₊>0 → "like_repel"` was a *positive Green-bilinear* (`q₊q₊G`), NOT the physical on-shell energy `Ω=−½⟨J,A⁻¹J⟩<0` (attract). The scalar never had Coulomb; the "like-repel" label was a sign-convention artifact.

---

## 2. Physics lessons (durable)

1. **Scalars attract, vectors repel (spin-statistics).** Coulomb like-repel is a spin-1 fact — but in *real* EM it comes specifically from the **`A⁰`/Gauss (charge-density) sector**, not from the transverse (light) vector. Naively putting charge on the transverse light-vector gives the *magnetic* (attractive) half.
2. **In THIS model, light and the Coulomb sector are NOT gauge-unified.** Light = brane transverse shear `u_T`; the electric sector is a *separate* leftover scalar. There is no single `A_μ` (never exact Maxwell). **Therefore the real-EM `A⁰`-gauge explanation of like-repel is unavailable and must not be imported** — any repulsion must come from the plain field's own stored-energy/boundary structure.
3. **The electric sign lives in the throat boundary condition (`𝔅`), which the committed model does not determine (U2).** Wherever the model *does* compute a boundary response, it computes **fixed-source → attract** (e.g. Codex's discriminant `∂A/∂K = −64/49 < 0`). So the like-repel sign requires either a *postulate* of the fixed-defect (clamp) BC, or its *derivation* from the (currently frozen) throat core.
4. **"Static freezing kills the physics" — again.** Coulomb repulsion needs the nonlinear/unfrozen throat core (the geon + sleeve microstructure the G0 MVP zeroed via `S_hold`, frozen `Σ`, zeroed geon couplings). Every static simplification deleted the sign-carrying structure.
5. **Charge ⊥ mass (model-native).** `conceptual_foundation §4`: mass = trapped geon standing wave → sources gravity via the **drain**; charge = the `±w` puncture direction; the two "do not directly interact at this level." Mass is *not* a property of the medium and is *not* the brane `u_L` mode.

---

## 3. The leading hypothesis: charge = a stored unit of the leftover longitudinal scalar

**Claim.** The electric force is carried by the **leftover longitudinal scalar `u_L`** — the real, gapless, second-class "stray" mode from pathA_36 that refuses to fold into Maxwell (`FAIL_CAUCHY_STRAY_LONGITUDINAL`). A throat **clamps** (stores) a signed fixed unit of `u_L` (oriented by `±w`); the electric force is the scalar's **stored gradient energy** `½T_L∫(∇u_L)²`.

**Derived sign (both Codex and Grok, target-blind, gauge-barred):**
- **Clamp / fixed-defect** (throat holds a nonrelaxable signed monopole; energy `~K·q²`): `U(R)=+4πT_L q²·s₁s₂/R` → **like-repel, opposite-attract, `1/R²` = Coulomb.** Positive cross-term `½T_L∫|∇(φ₁+φ₂)|²`, **no gauge structure used**.
- **Source / fixed-source** (throat emits; energy `~J²/K`): `U(R)=−(J₀²/4πT_L)·s₁s₂/R` → like-attract (the excluded branon).
- Discriminant (mouth-ensemble fork): `d(coupling)/d(wall stiffness K)` — *strengthens* ⇒ clamp ⇒ repel; *weakens* ⇒ source ⇒ attract.

**Why the clamp is physically motivated (not arbitrary):** a **geon is a standing, trapped wave = stored energy**. A trapped object holding a fixed oriented unit of the scalar *is* a clamp; an emitter would radiate. So "energy stored in the scalar in the direction of the throat" is what a geon does by its nature. Turning this from plausibility into derivation needs the **unfrozen geon core**.

**Why it survives where the others died:**
- Gapless `u_L` → `1/R²` range available (no Yukawa problem) — *if* a nonrelaxable monopole exists.
- No gauge/`A⁰` needed → consistent with "light and charge not unified here."
- Mass-consistent (see §4) and aligned with the model's own "charge ⊥ mass."
- It is the **first electric route that is not structurally excluded** — it is merely *underived* (the clamp BC), which is a derivation/consistency problem, not a wrong-sign wall.

---

## 4. The "mass channel" debunking (why the EP objection evaporates)

A prior review objected that `u_L` "is the mass channel (`q_M → u_L`)" and so cannot also carry charge (equivalence-principle conflict). **An objective code-level audit found this label is a retrofit, not a derivation:**

- `q_M` is a **hand-placed free symbol**: `J_mass = sp.Matrix([0,0,q_M,0])` (pathA_39 Stage 4 tool ~line 621; pathA_42 tool asserts `source_mass_vector=[q_M,0]`). It has **no formula, no provenance, undetermined `SIM_GATED` magnitude** — unlike the *computed* couplings `q_h, q_L, q_A` in the same derivation.
- It was inserted only to give the leftover scalar "a mass leg to run the equivalence-principle test against," motivated by the loose "mass = energy density → density mode."
- **Real mass** = the throat/geon standing wave sourcing gravity via the **drain** (`pathA_29`; `conceptual_foundation §4` lines 436, 478–479), a bulk `w`-transport object — a *different DOF* from the brane `u_L`.
- The model already had `u_L` carrying both a charge coupling `q_L` and the mass-leg `q_M` as a **`SIM_GATED` magnitude question, never a no-go**.

**Conclusion:** using `u_L` for charge does **not** genuinely conflict with mass; the "EP conflict" was an artifact of the mislabel. **Residual (honest, non-blocking):** a mass, being an energy density, will perturb `u_L` slightly → a possible small charge↔mass admixture (`MIXED_SCALAR_EP_RISK`), sim-gated and at most gravitational-strength — a *magnitude to characterize* (possibly a small predicted charge-mass coupling), not a blocker.

---

## 5. Forward plan (the new direction)

Adopt the clamp as the throat's charge-ansatz (well-motivated by geon = stored energy) and **work out what it forces** — predictive-surplus test, per the calibrate-and-build methodology (falsification = a NO-GO between requirements):

1. **Electric force** — derived: `1/R²` like-repel from the stored `u_L` gradient energy. (In hand, conditional on the clamp.)
2. **Magnetism sign** — a *moving* clamp = a current; does the same clamp postulate give the right magnetostatic sign? (Cross-check vs pathA_39.)
3. **Hierarchy** — clamp (electric) energy vs drain (gravity) energy: does the ratio give the `~10⁴²` electric/gravity hierarchy for free? (The "same wall stiffness sets both" thread.)
4. **EP admixture** — characterize the magnitude of the charge↔mass `u_L` admixture; confirm it is ≤ gravitational-strength (or a testable prediction).
5. **Upgrade postulate → prediction** — try to *derive* the clamp (the nonrelaxable signed `u_L` unit) from the **unfrozen geon/sleeve core** the MVP zeroed. This is the "unfreeze the core" move; it is where the sign stops being a postulate.

**Reusable assets:** the committed static-prerequisite check suite (`bdd0104d`, both engines) remains sound machinery. The `𝔅`/mouth-ensemble apparatus is exactly the clamp-vs-source fork.

---

## 6. Provenance (raw derivations)

Working logs (under `software/em_charge_attribute/_scratch/`, gitignored — this doc is the durable distillation): `g0_powercounting_{grok,glm}.log` (field-exchange exclusions); `g0_card_gauntlet_r1.log` (NO-GO); `electric_vector_{hypothesis_v0.md,codex_review.log}` + the fresh-Claude review (vector route excluded); `wall_tension_{sign_consult.md,codex.log}` + fresh-Claude (wall-tension NEEDS-X); `leftover_scalar_charge_derivation.md` + `leftover_scalar_{codex,grok}.log` (the clamp derivation, target-blind); and the mass-channel audit (fresh-Claude agent). Memory anchor: `[[reference-two-throat-sim-spec]]`, `[[project-em-sector-reconsideration]]`.
