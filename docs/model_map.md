# The Model Map — one-medium analog, derivation atlas + conceptual throughline

**What this is.** The single high-level map of the whole toy model: the conceptual picture, every earned derivation with a one-line note + pointer to the full doc/scripts, the honest ledger of what is *predicted* vs *calibrated* vs *unresolved (R1)* vs *departure*, and a glossary. Read this to hold the model in your head; read the cited sources to act on any specific number.

> **This is a synthesis map, not the source of truth.** Assembled 2026-07-21 from a six-agent sector fan-out over the committed repo. For any discrepancy the **cited source files are authoritative** — the sector reports under `software/stage1_solver/reports/` and `software/em_charge_attribute/`, the v2-ledger stages under `research/pde_ledger_v2/notes/stages/`, the blueprint `notes/ledger_v2_blueprint.md`, the trackers `research/pde_ledger_v2/notes/{parameter_register,midway_knob_audit}.md`, and the resume doc `research/pde_ledger_v2/notes/RESUME_ROADMAP.md`. Re-read those before trusting any claim here.
>
> **Framing.** This is a *toy analog*: the goal is ONE self-consistent structure (a calibrated PDE) that reproduces GR-like and EM-like far-field behavior — a working math **bridge**, not an ontology claim. A result that **breaks** the concept is welcome and first-class. Magnitudes are calibrated; only held-out **dimensionless** structure tests the model. A clean "it all works" is suspicious.

---

## 1. The model in one picture (the throughline)

Everything bottoms out in **ONE compressible superfluid** — a Gross–Pitaevskii / nonlinear-Schrödinger (GNLS) condensate in a 4D bulk, order parameter `ψ`, number density `ρ=|ψ|²`, closed by a stiff polytropic equation of state `P=Kρ⁵` (postulated closure).

That one medium exists in **two material states**, distinguished by an order field `χ_B∈[0,1]`:
- **the brane** (`χ_B=1`) — the ordered, *shear-supporting* phase; our 3D space is a domain-wall brane at `w=0` in the 4D bulk;
- **the bulk** (`χ_B=0`) — the SAME medium, de-structured and *shear-free*.

⚠ This two-phase split is **POSTULATED, not derived** — the parent potential `U(ρ)` is single-well, so the brane does not *emerge* from the one medium; it is put in by hand (see §3.1). Read the throughline below as the *intended architecture*, with the honest earned-vs-postulated status in §3–§4.

A **particle is a defect** — a topological **throat** that punctures the brane into `±w` (a phase-conversion site / open channel, not a suction force). "Defect", "throat", "puncture" — never "dent" (a dent doesn't puncture).

The four forces are **modes of the one medium**, not four theories bolted together:

| Force | Native mechanism | Carrier / cone |
|---|---|---|
| **Gravity** | the **drain/flow** between defects — each throat removes number-flux, two drains sit in each other's inflow → attraction (`1/r²` on the brane, `1/R³` in the bulk). NOT a ripple. | density mode, `c_s` |
| **Light** | in-plane **MacCullagh shear** of the brane (energy in `curl u`, on the brane so it leaves the bulk/throat undisturbed) → two transverse photons | brane shear `u_T`, `c_γ` |
| **Charge** | the **static ±w throat**: the puncture bends the brane into `±w`; the charge **sign = the ±w orientation** (Z₂-topological, NOT additive) | electric scalar `h`, `c_E` ⚠ |
| **Magnetism** | the **moving ±w throat** — the electric twin; a moving signed throat shears the brane through the *same* transverse sector as light | brane shear `u_T`, `c_γ` |

⚠ `c_E` has **no committed lock** to `c_γ` (`c_s/c_E=2` was chosen to *avoid* forcing one) — do NOT read the table as asserting `c_E=c_γ`; that equality is the open `r_cone` question (§2, §3.6). And gravity's drain/return is **non-variational** (a source/BC, not a Lagrangian term) — which is why the unified structure is a *synthesis*, not a sum (§2).

The unifying test is **boost-consistency / emergent Lorentz**: does *one* calibration make the electric and magnetic sectors two faces of one interaction? The magnetism sector is built precisely to probe this.

**Honest status in one breath:** all four far-field sectors now exist on one shared field set. Gravity is the most complete (form + fingerprints earned, magnitudes calibrated, one held-out falsifiable departure). Light earns two transverse photons but leaves a characterized stray-longitudinal departure. Charge and magnetism earn the *structure* (mechanism, falloff, tensor form) target-blind, but the electric **sign** — and with it the magnetic sign — is genuinely **unresolved (R1)**, waiting on a sim-deferred nonlinear throat solve. The model is a coherent, calibratable single-medium analog with a real, un-papered-over honesty ledger.

---

## 2. The shared substrate (fields, speeds, and why it is NOT a sum)

**Primitive constants of the medium (~4 declared universal):** `ħ` `[ML²T⁻¹]`, `m_GNLS` `[M]`, `K` `[ML¹⁸T⁻²]`, `ρ0` `[L⁻⁴]`.

**Shared field set:** `{ρ, θ, r_B, H/h, u_L, u_T}`
- `ρ=|ψ|²` — 4D number density `[L⁻⁴]`; `ρ0` = asymptotic datum. (In the two-phase writeup `n`, with `n_B=χ_B n` the ordered fraction.)
- `θ` — condensate phase (dimensionless); flow `v_b=(ħ/m)∇θ`.
- `χ_B` — two-phase order field `∈[0,1]`, brane=1/bulk=0; **independent scalar, explicitly NOT `|P_∥|²`**.
- `r_B / ρ_B` — wall (brane) amplitude / projected on-brane density `[L⁻⁴]` (`ρ_B=∫W χ_B n dw`); `ρ_br` = brane mass density `[ML⁻³]`.
- `H / h` — electric scalar; `H` parent `[L⁻¹]` localized by a Pöschl–Teller operator, `h=P₀H` its one gapless zero mode; embedding identity `ξ_w=ℓh`.
- `u_L` — longitudinal brane displacement (the light sector's stray mode).
- `u_T` — transverse brane displacement `∇·u_T=0` (the two photons; reused by magnetism).
- `P` — polar orientation field ("little arrows") — **RETIRED as an independent field by Decision 16** (see §6); survives only as the head≠tail `+w≠−w` polarity migrated to the throat sleeve.

**The three (+1) speeds — state carefully:**
- `c_s² = 5Kρ⁴/m` — phonon / density / **gravity-change** speed; DERIVED from the EOS (factor 5 from `dP/dρ`); `d ln c_s/d ln ρ=2`.
- `c_γ² = μ_R/ρ_br` — **light cone**, MacCullagh shear modulus over brane inertia; **SHARED by magnetism**.
- `c_E` — electric speed, `K_h=M_h c_E²`; **NO committed cone lock** (`c_s/c_E=2` was *chosen* to avoid forcing one).
- `v_b=(ħ/m)∇θ` — condensate flow. And `c=c_γ` is the wave-sector ceiling (earned from wave kinematics, NOT from `E=mc²`).
- ⚠ `λ_γ=c_γ/c_s` is a **DERIVED ratio, NOT a free knob** — the free-unreduced content is `c_γ`, ultimately `{μ_R, ρ_br}`. The `c_±²=(3±√2)/2` witness values are a G0 witness point in `c_E=1` units, **not** universal cones; there is no `c_L²=B_eff/ρ_br` "lock pair" (the scalar-block inertia is `A_eff=ρ_br+C_J²/κ_phase`).

**Why the unified structure is NOT an additive sum.** The four sectors share ONE field set, but Part VII must **synthesize**, not write `(medium)+(gravity)+(light)+(charge)+(magnetism)`:
- **gravity's drain/return is NON-variational** — carried as source / boundary-condition bookkeeping, not a Lagrangian term;
- the G0 card's `r_B` wall action must be **reconciled/de-duplicated** with Part I's `χ_B` wall action (compounded by the pending Decision-16 retirement).

---

## 3. The derivation atlas

Status tokens: **EARNED** (target-blind, within stated postulates) · **CALIBRATED** (magnitude tuned to an external anchor) · **R1_REQUIRED(x)** (unresolved, needs the named closure — usually the throat solve) · **DEPARTURE** (a characterized, first-class break from the reference theory) · **SIM-DEFERRED** · **cite-only** · **SUPERSEDED**.

### 3.1 Medium — Part I ✅ built (stages 004–007), tri-reviewed

*Lays the shared substrate; earns no force sector by itself. The brane/two-phase split is POSTULATED, not derived from the one medium.*

- **stage 004** — GNLS parent action + dimensional foundation. **EARNED (dim foundation).** Fixes `{L,T,M}` base + field dictionary (17 checks); mass is an explicit constant (not emergent). `[ψ]=L⁻²`, `[K]=ML¹⁸T⁻²`; `ξ_h=√2ħ/(mc_s0)`; `h0=mc_s0²/4`. → `research/pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md` (+ `scripts/…_sympy_audit.py`, `mathematica/…_mathematica_audit.wl`). *Open:* parent action & `m_GNLS` POSTULATED; `m_defect` NOT emergent (`INFLOW_MASS_SOURCE_MISSING`).
- **stage 005** — EOS sound speed + light/sound ratio. **EARNED (within imposed EOS), CALIBRATED landing.** `c_s²=5Kρ⁴/m` (factor 5 derived); `c=c_γ` ceiling. The **free-unreduced quantity is `c_γ`** (ultimately `{μ_R,ρ_br}`); the *ratio* `λ_γ=c_γ/c_s` is DERIVED, and only the *constraint* `λ_γ=1` is a calibrated cone lock (`pathA_40`). A reversible negative control proves `c_γ` is unpinned by the parent action (`C_GAMMA_RATIO_UNDERDETERMINED`). → `…/ledger_stage005_sound_speed_light_ratio.md`. *Open:* EOS `P=Kρ⁵` IMPOSED; `c_γ` numeric undetermined; `ħ` provenance a gap.
- **stage 006** — Two-phase material-state ontology (`χ_B`). **`ACTION_SPECIFIED_CLASSIFIED` (postulated functional).** Two-phase free energy + order-balance PDE (13 pins); single-kink wall admission; recovery reduction to the old `S_leak` law (assert-zero). → `…/ledger_stage006_two_phase_chiB_ontology.md`. *Open/departure:* **whole two-phase split POSTULATED** (parent `U(ρ)` is single-well); **θ-as-Maxwell-φ = FATAL_FLAW** (wrong-sign stiffness; Maxwell only `BY_TUNING`); slab width `W_slab` FREE-UNREDUCED. `α_aniso` **retired by Decision 16** (amendment landed 2026-07-21; stage006 operative `DRIFT(6)→DRIFT(5)`).
- **stage 007** — Shear-surface G0 freeze (frozen brane/light action). **EARNED (freeze fidelity + historical DOF=8 / operative DOF=4 + dim firewall), POSTULATED/CALIBRATED landing.** Hash-guarded frozen MacCullagh shear + P–u coupling + `u_w` gap; honest "second-medium drift = 11". `c_γ²=μ_R/ρ_br`. → `…/ledger_stage007_shear_surface_g0_freeze.md`. *Open:* the "11" = 4 constants {ρ_br, μ_R, λ_Pu, Ω_w} + `g_ℓ` + 6 postulates; `μ_R/ρ_br` Route-A reduction PENDING (R10); the `P`-machinery here was **retired by Decision 16** (amendment landed 2026-07-21): the freeze-as-run "11"/DOF=8 STANDS as the immutable historical tier, but the operative post-D16 layer is DOF=4, drift 7 (`POST_D16_DRIFT(7)`), `λ_Pu` retired — a −5 route-less reduction.

### 3.2 Gravity — Part II ✅ built (stages 001–002, 008–029), tri-reviewed, sector CLOSED

*Gravity = the drain/flow between defects. Calibrate-predict: form/sign/fingerprint earned target-blind; magnitudes (`G`, `54/5`, `2/5`) calibrated to a separately-audited PN corpus.*

- **pathA_21c** — inter-defect force from the Noether stress tensor. **EARNED (form+sign), CALIBRATED magnitude.** Like drains **attract**; `F_12 = −(m N_∞,3 Q1Q2/4πr²) r̂` (reduced-3D `r⁻²`; bulk `R⁻³`). → `software/stage1_solver/reports/pathA_21c_force_from_noether_stress_tensor.md`; folded into ledger stages 001/002. *Open:* full sign residual (quantum/`V_conf`/Maxwell pieces not all evaluated); `G` form not derived.
- **pathA_28** — monopole/dipole radiation constraint spec. **EARNED (DtN amplitudes/orders); constraint-spec not suppression proof.** `A_raw(ℓ0)∝ω`, `(ℓ1)∝ω³`, `(ℓ2)∝ω⁵` (the `i Q2 a⁵ω⁵/27c_s⁵`). → `…/pathA_28_monopole.md`; ledger stage008. *Open:* `cancellation_possible` is a literal flag, not computed.
- **pathA_29** — flat-slab brane↔bulk return residual. **held-out FALSIFIABLE PREDICTION** (`RETURN_RESIDUAL_PREDICTION`). A finite return leaves a bounded, lower-order ℓ=0/1 `c_s`-radiation residual: **`1−T_ℓ(0)=ε_ℓ/(1+ε_ℓ)`**, orders `p_res=1,3` — the departure from GR (which forbids ℓ=0/1). Localization: `p=2` (1/r² survives the slab) vs anti-localizing warp `RETURN_NOGO` (able-to-fail). → `…/pathA_29_brane_bulk_return.md`; ledger stages 009/010. *Open:* slab geometry POSTULATED; `ε_ℓ` FREE-UNREDUCED; full closure = Gate-6.
- **pathA_30–34** — moving-throat reduced-closure Gates 1–5. Gate1 DtN unit test (`DN_UNITTEST_BC_DEPENDENT`); Gate2 breathing mode (`BREATHING_CALIBRATED`); Gate3 SO(3)/ℓ=2 isotropy — EARNED slice (`−Δ_S²Y_2m=6Y_2m`, lane defects=0) + **`ISOTROPY_CALIBRATED`** landing (`{T_Ω,β₂}` calibrated pending R36); **Gate4** ℓ=2 quadrupole fingerprint — **EARNED:** `Ŷ₂ᵒᵘᵗ=1+z²/9+4z⁴/81+i z⁵/27+…`, **χ_Q=+1** outgoing; Gate5 cross-ℓ **honest NEGATIVE** (`FAIL_UNDERDETERMINED_NOT_PREDICTIVE` → `CROSS_L_RESIDUAL_PREDICTION`: the linear model does NOT pin the ℓ=0/1 return magnitude; needs a Gate-6 selector). The earned **cross-ℓ radiation fingerprint** is `{1, 1/2, 1/27}` at `{ω¹, ω³, ω⁵}` (the ℓ=0/1/2 leading amplitudes). → reports `pathA_3{0,1,2,3,4}_*.md`; ledger stages 011–023. 
- **pathA_43** — density-mode ℓ=2 quadrupole radiative port. **`DENSITY_PORT_HOSTED`.** Builds the ℓ=2 port natively on the density/`c_s` mode as a two-port over (`q2` wall, `Φ2` bulk) — **retires the old `A_w` vector scaffold**; vector-freedom + `a⁻⁵` + outgoing sign `+i z⁵/27` pass. → `…/pathA_43_density_quadrupole_port.md`; ledger stages 024–027. *Open:* moment/coefficient values SIM-DEFERRED; `G`, `54/5`, `2/5` CALIBRATED.
- **pathA_2.5PN matchback** — Burke–Thorne consistency cap. **`MATCHBACK_CONSISTENT` (CALIBRATED).** Density-port moments reproduce `Γ̄₅=2G/(5c⁵)` and `K̄4K̄0=4K̄2²`; the `/27` is the earned fingerprint. → ledger stage028. 
- **PN corpus (1PN→4PN + 2.5PN)** — **cite-only.** Seven Zenodo-DOI'd, separately-audited, GR-matched papers; `G=GENUINE_BLOCKED`; cited not re-derived. → ledger stage029; dirs `research/1pn_orbital_dynamics`, `research/4d_{1pn_bridge,1pn_full,2pn,2_5pn,3pn,4pn}`.

**Gravity's held-out surplus:** the bounded ℓ=0/1 `ε_ℓ/(1+ε_ℓ)` residual (THE falsifier — GR/Birkhoff forbid it), the ℓ=2 rational fingerprint `1/9,4/81,1/27` + `χ_Q=+1`, the cross-ℓ fingerprint `{1, 1/2, 1/27}` at `{ω¹,ω³,ω⁵}`. **NOT held-out (calibrated/external):** `G`, `54/5`, `2/5`, force-norm, and the Gate-2/3 wall packet `{μ_η,T_w,β,Vp0/ℓ_c,T_Ω,β₂}`.

### 3.3 Light — Part III 🌱 seeded only (pilot stage003)

*Light = in-plane MacCullagh shear of the brane; two transverse photons. Least affected by the EM reconsideration — it IS the shared transverse sector, and it establishes the `u_T`/`c_γ` foundation magnetism reuses.*

- **pathA_35 gate L** — couple-stress no-go. **DEPARTURE (`FAIL_COUPLE_STRESS_NOGO`).** Light's shear stiffness `μ_R` **cannot** be sourced from a polar-field `P` substructure (live `P`→hidden spin-waves/unbounded; gapped→nonzero residual; slaved→degraded provenance). Feeds `c_γ²=μ_R/ρ_br` forward. → `software/stage1_solver/reports/pathA_35_gateL_light.md` (+ freeze `pathA_35_G0_freeze.md`). *This is the direct evidence behind Decision 16.* `μ_R` is honestly a POSTULATED modulus.
- **pathA_36** — C5 phase-potential derivation. **EARNED (2 transverse photons) + DEPARTURE.** Two massless photons at `c_γ²=μ_R/ρ_br` (`PASS_TRANSVERSE_UNDISTURBED`, able-to-fail); longitudinal half → **`FAIL_CAUCHY_STRAY_LONGITUDINAL`**: a Dirac–Bergmann **second-class** pair (one propagating stray longitudinal DOF), NOT Maxwell's first-class Gauss. Maxwell locus reachable `BY_TUNING` only. → `…/pathA_36_c5_phase_potential.md`.
- **stage 003 (seed)** — transverse photons + stray longitudinal (the in-ledger FAIL-headline pilot). **seeded (CHARACTERIZED-DEPARTURE).** Re-derives pathA_36 in-ledger, closing its fidelity gap (the Josephson sign `C_J=−Jρ_B0` derived in-script, not hardcoded). → `research/pde_ledger_v2/notes/stages/ledger_stage003_transverse_photons_stray_longitudinal.md`. ⚠ its "Next step" cites a stale `c_L²=B_eff/ρ_br` cone-pair — superseded (see §2).

**Still to build:** `part3_light_atomic_split.md` (absent) → per-Part gate → ~4–6 stages from pathA_35+36+seed.

### 3.4 Charge — Part IV ⬜ not started (built + verified off-ledger)

*Charge = the static ±w throat (puncture-deflection): the puncture bends the brane into ±w; sign = ±w orientation (Z₂). Coexists with gravity on the same brane (`internal_inconsistency=none`, computed).*

- **puncture-deflection electric-sign** — **R1_REQUIRED(bc_selection).** A ±w puncture bends the brane (`ξ_w=ℓh`), sourcing `h` via `g_χh`; the far-field two-body force `F_X=s₁s₂A_X/4πR²` is **target-blind** but its **sign is not invariant across admissible boundary-condition ensembles**: `A_V` neutral-**positive** (repel), `A_J=−m_gg(j+g)²` **negative** (attract), `A_M=m_gg(q²−g²)` **indefinite** (rides free `q/g`). So §4 fires `bc_selection`. → `software/em_charge_attribute/puncture_deflection_electric_sign_result.md` (+ `directive_…md`, `…_check.py/.wl`, independent recompute). Verified: dual-engine + 4-leg + Grok. **The sign is neither earned NOR calibrated — the formal result is R1** (the informal "calibratable sign" is a user framing, not the build's verdict).
- **G0 shared minimal closure card v0** — **DRAFT v0 POSTULATED closure.** Fixes parent physics once (bulk Madelung, wall `r_B`, localized `H`+PT zero mode, `χ↔h` coupling, `M_h/K_h/C_hu` with `C_hu²<B_eff K_h`); only the sleeve-contact operator `𝔅` varies (E2/E3). → `software/em_charge_attribute/g0_closure_card_v0.md` (+ scope). *Postdates Part I; duplicates its bulk/wall/scalar structure → Part VII must de-dup, NOT additive-merge.*
- **Superseded predecessors:** leftover-scalar `u_L`-clamp (sign = a calibrated knob; `NO_NATIVE_CLAMP`) → `…/leftover_scalar_electric_sign_result.md`; defect/antidefect wall-healing → `docs/em_phaseC_force_decomposition.md`; **pathA_38** throat-body Coulomb (`THROAT_ELECTRIC_LOCALIZED_COULOMB`, the *old* Part-IV anchor) → `software/stage1_solver/reports/pathA_38_throat_body_electric_localization.md`.
- **DEPARTURE — `NATIVE_P_NO_EMERGENT_GAUSS`:** exact U(1)/Maxwell Gauss proven **non-native** (Dirac–Bergmann second-class at quadratic order; no first-class chain). → `software/em_charge_attribute/reports/native_p_constraint_gate.md`.

**Resolving the R1:** only the sim-deferred **nonlinear throat/core solve** (`s↦h_A` boundary functional/barrier) can pin the class, force a holder, and (upside) deliver `q>g`.

### 3.5 Magnetism — Part V ⬜ not started (built + verified off-ledger, commit 53cf049f)

*Magnetism = the moving ±w throat (the electric twin). Two independent far-field derivations; their COMPARISON is the boost-consistency / emergent-Lorentz result.*

- **Q-CURRENT (native source)** — **EARNED (tensor form).** From signed-throat continuity, the source is derived **natively** — `J_T=q_T s η V` — NOT the old barred `j∝sV`. Supplies the ledger-ready action row (no pre-existing G0 row changed, `internal_inconsistency=none`): `S_{T+move}=∫[½ρ_br|u̇_T|² − ½μ_R|∇×u_T|² + q_T Σ s_i η(x−X_i)V_i·u_T]`, `∇·u_T=0`, `c_γ²=μ_R/ρ_br`, `q_T=λ_T τ_d`.
- **Q-BOOST (Route A)** — **EARNED (structural), tier-A conditional.** Boosts the electric anchor to the transverse Darwin interaction: `U_A=(s₁s₂A_E/4πR)[1−(D_V+A_V)/2c_γ²]`, kernel `(δ_ij+n_in_j)/8πR`. Inherits the electric `A_E` R1.
- **Q-DIRECT (Route B, BLIND to A)** — **R1_REQUIRED(direct_moving_throat).** Solves the moving-throat field directly (`foreign_payload=None`, independence enforced); reproduces Route A's kernel: `U_B=−s₁s₂q_T²(D_V+A_V)/8πμ_R R`. Magnitude carries the unresolved `q_T`.
- **Q-COMPARE (the crux)** — **R1_REQUIRED(electric_bc_selection).** Tensor structure, `R⁻¹/R⁻²` falloff, and `O(V₁V₂)` order **AGREE target-blind**; full equality reduces to a single ratio `r_BA=q_T²/(ρ_br·A_E)` that is an **unresolved output** (needs `r_BA=1` AND cone `r_cone=1`, plus `q_T`, higher `v/c`, active-flux integrability — all open). → `software/em_charge_attribute/magnetism_moving_throat_result.md` (+ `…_check.py/.wl`). Verified: dual-engine + fresh-agent legs + GLM. The formal landing is the terminal **`R1_REQUIRED(electric_bc_selection)`** plus three explicitly-emitted co-blockers **`R1_REQUIRED(direct_moving_throat)`**, **`R1_REQUIRED(magnitude)`**, **`R1_REQUIRED(consistency)`** (informally "doubly-R1": it needs the throat solve AND the electric bc selection); the magnetic sign inherits the electric R1.
- **DEPARTURE — `B_TIME_REVERSAL_EVEN`:** `b_T=∇×u_T` is time-reversal **EVEN** (Maxwell `B` is T-odd), correctly axial — a concrete self-consistent prediction; magnetism requires the active-drain time-arrow `τ_d`.

**Still to build (IV+V):** re-scope onto these builds (the blueprint's pathA_38/39 scope is superseded — `pathA_39` rested on the barred `j∝sV`); author `part{4,5}_*_atomic_split.md` → per-Part gate → build.

### 3.6 Knit — Part VI ⬜ not started (gates pathA_40/41/42 exist, need reshape + re-adjudication)

*The knit VERIFIES "one medium, all emergent" instead of asserting it: shared cones? one irreducible parameter set with no smuggled second substance? does the propagating scalar stay consistent?*

- **pathA_40 — cone-lock.** **CALIBRATED-LOCK (`Δr=2`), to-RE-ADJUDICATE.** Neither Lorentz/Maxwell lock is derived on the earned ledger; both are calibrated: **(A) `λ_γ=1`** ⟺ `μ_R/ρ_br=5Kρ⁴/m`, **(B) `c_E=c_γ`** ⟺ `c_E²ρ_br=μ_R`. ⚠ The newer G0 card + magnetism build do **NOT** establish these — no committed cone lock — so Part VI must re-adjudicate. → `software/stage1_solver/reports/pathA_40_cone_lock.md`.
- **pathA_41 — NG5 second-medium drift.** **INCOMPLETE (`SECOND_MEDIUM_DRIFT`).** The one-medium claim doesn't fully close: `{ρ_B0, χ_c, C_hu}` remain irreducibly independent (unreduced brane-surface/embedding parameters — NOT a separate substance; `no_fourth_arena=True`). `{ρ_br, μ_R, c_E}` reducible-in-principle (throat solve). *This is where the knit's real falsification power lives — the cone-lock is near-vacuous; NG5 is the decisive reducibility test.* → `…/pathA_41_ng5_second_medium_drift.md`. (Note: pathA_41's drift is this NG5 trio `{ρ_B0,χ_c,C_hu}` — distinct from the stage007 `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)`, whose "11" was contested by an earlier GLM catch but is now RESOLVED — the count STANDS, no `ρ_br` overcount, dual-engine + anti-absorption-guarded per `parameter_register.md` + stage007.)
- **pathA_42 — charge-coupled scalar map.** **INCOMPLETE / SIM-GATED (`SCALAR_DEPARTURE_MAPPED`).** The propagating `h`-scalar doesn't clearly break the model, but its break-risk magnitude is sim-gated: `h_EP` earned-safe on the decoupled floor; radiation/universality/`u_L_EP`/preferred-frame all `SIM_GATED`. Sharp trade-off: `c_E→∞` non-radiating but preferred-frame; `c_E=c_γ` Lorentz-invariant but a real radiating extra scalar. → `…/pathA_42_charge_coupled_scalar.md`.

### 3.7 Integration — Part VII ⬜ not started (the literal deliverable)

The genuinely-new synthesis, needs I–VI assembled:
- the **unified equation set** — conservative parent action PLUS explicit **non-variational** drain/return sources/controllers/BCs (NOT one functional) — the completeness floor;
- the **calibration map** (every constant tagged DERIVED / INPUT / gap / benchmark);
- a whole-system **dimensional-firewall** check;
- the permanent **registers** (open-items / no-go / calibration — the falsification record);
- the **simulation hand-off spec** (each sim-dependent quantity + its exact equation + confirming measurement);
- the **χ_Q reconciliation** (`pathA_22b≈0.712` numeric vs `pathA_33=1` exact — same name, different computations; reconcile, don't merge);
- the full irreducible count via the pathA_40 `Δr=2` codimension technique.
Gets the full new-derivation gauntlet + GLM tertiary + user gate.

---

## 4. The honest ledger (surplus vs debt)

**Held-out / falsifiable predictions (the earned surplus, ~6–7 structural + 1 departure):**
- the gravity `1/r²` exponent + attractive sign;
- the ℓ=2 DtN rational fingerprint `{1/9, 4/81, 1/27}`, `χ_Q=+1`, the cross-ℓ fingerprint `{1, 1/2, 1/27}`, SO(3) `λ_m=6`, and the squared-denominator signature (the `−2D₂N₀` term). (The `27` is earned separately as part of the DtN fingerprint; in `54/5=2·27/5` the `27` is earned but the `2/5` is external/calibrated — these are two distinct results, not one.)
- **the bounded monopole/dipole `c_s`-radiation residual `ε_ℓ/(1+ε_ℓ)`** — the headline falsifier (GR/Birkhoff forbid ℓ=0/1).
- ⚠ `F_e/F_g` (the hierarchy capstone) is **NOT yet held-out** — it becomes the sharpest dimensionless test only once BOTH couplings come from one throat action; presently FIT/not-tested.

**Characterized departures (first-class, never softened):**
- EM is **NOT exact Maxwell** — `NATIVE_P_NO_EMERGENT_GAUSS` (exact U(1) proven non-native);
- the **stray longitudinal DOF** — `FAIL_CAUCHY_STRAY_LONGITUDINAL` (second-class, not first-class Gauss);
- `b_T=∇×u_T` is **time-reversal EVEN** vs Maxwell's T-odd `B`;
- a real **radiating extra scalar** `h` (accelerating charge radiates scalar waves);
- light's stiffness `μ_R` **cannot** come from `P` substructure — `FAIL_COUPLE_STRESS_NOGO`.

**Unresolved (R1) — the shared throat solve is the crux:**
- electric **sign** `R1_REQUIRED(bc_selection)`; magnetic sign inherits it;
- magnetism `r_BA=1?`, `r_cone=1?`, `q_T`, higher `v/c`, active-flux integrability;
- cross-ℓ return selector (Gate-6): the linear model does not pin `ε0/ε1`.
- **One nonlinear throat solve** is the shared R1 for gravity `{μ_R,ρ_br}` (audit R10 + R30 + R33 — *not* all six), electric `bc_selection`, and magnetism `q_T` — one interior solve collapses several knobs at once. SIM-deferred (tractability), not a knob the PDE sets.

**Calibrated / external (not held-out):** `G` (`GENUINE_BLOCKED`), `54/5`, `2/5`, the force-magnitude norm, the six Part-II wall-packet CALIB inputs `{μ_η, T_w, β, Vp0/ℓ_c, T_Ω, β₂}`, and the calibrated closure moments / `Γ̄₅`. (Not `λ_γ` — the *ratio* is DERIVED; what's calibrated is the *lock* `λ_γ=1`, and the free-unreduced quantity is `c_γ`.)

**Midway Knob Audit headline (Parts I–II dry-run, `midway_knob_audit.md`; post-Decision-16 −5):** irreducible total **~34–43** = universal 4 + reduction-debt 15–18 (route-ful) + postulated structure 14–20 (route-less) + force-mag 1. **Route-less liability ~15–21 ≫ held-out surplus ~6–7.** Sobering but **NOT a no-go** — a substantial portion of the debt is route-ful (reduction debt ~15–18: the throat solve + Gate-6 return closure convert CALIB→DERIVED in one stroke); post-Decision-16 the route-less liability (~15–21) is now **comparable to, no longer clearly exceeding,** the route-ful debt (~15–18), though it still robustly dominates the earned held-out surplus (~6–7). **Two parked user decisions** (C1/C2 counting convention; R35 label) pick the point in the range and must resolve before Part VII.

---

## 5. Ledger-v2 build status (thin — see RESUME_ROADMAP for detail)

Central equation-set rebuild at `research/pde_ledger_v2/` (branch `ledger-v2-rebuild`). Stages numbered by BUILD order (001–029), not Part order.

| Part | Sector | Status |
|---|---|---|
| 0 | Conceptual | scaffolding-only (placeholder `.tex`) |
| I | Medium | ✅ built (004–007), tri-reviewed |
| II | Gravity | ✅ built (001–002, 008–029), CLOSED |
| III | Light | 🌱 seeded (pilot 003); atomic-split unwritten |
| IV | Charge | ⬜ not started (build verified off-ledger) |
| V | Magnetism | ⬜ not started (build verified off-ledger) |
| VI | Knit | ⬜ not started (gates 40/41/42 need reshape) |
| VII | Integration | ⬜ not started |

**Build order:** III → IV → V → VI → VII. **Immediate next:** ✅ the Decision-16 Part-I amendment LANDED (2026-07-21) → author `part3_light_atomic_split.md` → per-Part user gate → build Part III. Each executable stage = the 6-artifact unit (note · TeX card · SymPy audit · independent Mathematica audit · source map · register entry), dual-engine both exit 0, per-stage Codex→Grok→Codex + tri-review, two user gates (per-Part + per-stage). Full detail: `research/pde_ledger_v2/notes/RESUME_ROADMAP.md`.

### 5.1 Ledger completion manifest — exactly what math folds into each remaining Part, and where it lives

*Everything below is verified to resolve on disk (2026-07-21). **Already folded:** Parts I+II = ledger stages 001–002 + 004–029, source = the §3.1/§3.2 pathA reports + PN corpus (stage003 is the Part III seed, folded into III — not I/II). This manifest covers the REMAINING Parts III–VII.*

**Reshape tiers (the per-gate cost of bringing a script into the ledger's print-only/assert-zero/independent-`.wl` contract):**
- **FULL** — the `pathA_*` tools carry an `argparse --compare` + JSON/YAML payload harness and a payload-mirror `.wl`; folding = strip the harness → print-only/assert-zero/exit-nonzero, AND re-author the `.wl` as a genuinely independent route.
- **LIGHT** — the `software/em_charge_attribute/` EM builds are already assert-heavy (only minor harness); folding = wrap into the 6-artifact stage form + confirm the `.wl` is an independent route.

**Part III — Light** → author `research/pde_ledger_v2/notes/part3_light_atomic_split.md` → ~4–6 stages (built from/after the in-ledger seed `ledger_stage003`).
- Source math: `software/stage1_solver/reports/pathA_35_gateL_light.md`, `…/pathA_35_G0_freeze.md`, `…/pathA_36_c5_phase_potential.md`
- Scripts (FULL reshape): `software/stage1_solver/tools/pathA_35_gateL_sympy.py` · `…/pathA_35_gateL.wl`; `…/pathA_35_G0_sympy.py` · `…/pathA_35_G0.wl`; `…/pathA_36_c5_sympy.py` · `…/pathA_36_c5.wl`
- Already in-ledger: `research/pde_ledger_v2/scripts/ledger_stage003_transverse_photons_stray_longitudinal_sympy_audit.py` · `research/pde_ledger_v2/mathematica/ledger_stage003_transverse_photons_stray_longitudinal_mathematica_audit.wl`
- Owed: fold both departures (`FAIL_COUPLE_STRESS_NOGO`, `FAIL_CAUCHY_STRAY_LONGITUDINAL`) first-class; carry `c_γ²=μ_R/ρ_br` forward as the `u_T` foundation Part V reuses.

**Part IV — Charge** → author `…/part4_charge_atomic_split.md` → ~3–4 stages.
- Source math: `software/em_charge_attribute/puncture_deflection_electric_sign_result.md`, `…/g0_closure_card_v0.md`, `…/reports/native_p_constraint_gate.md`
- Scripts (LIGHT reshape): `software/em_charge_attribute/puncture_deflection_electric_sign_check.py` · `…_check.wl` (+ `…_independent_recompute.py`); `…/g0_closure_card_v0_checks.py` · `…_checks.wl`; `…/native_p_gate_sympy.py` · `…/native_p_gate_dual.wl`
- Owed: fold the `(u_L,h)` block + `ξ_w=ℓh` + h-mediated `1/R²`; carry `R1_REQUIRED(bc_selection)` (sign unresolved) + the `NATIVE_P_NO_EMERGENT_GAUSS` departure; **G0-vs-Part-I de-dup is deferred to Part VII** (reconcile, not additive-merge).

**Part V — Magnetism** → author `…/part5_magnetism_atomic_split.md` → ~4–6 stages.
- Source math: `software/em_charge_attribute/magnetism_moving_throat_result.md`
- Scripts (LIGHT reshape): `software/em_charge_attribute/magnetism_moving_throat_check.py` · `…_check.wl`
- Owed: fold the `S_{T+move}` transverse-vector row + native `J_T=q_T sηV`; carry `R1_REQUIRED(electric_bc_selection)` + its 3 co-blockers + the `b_T` T-even departure; reuse Part III's `u_T`/`c_γ`. (Blueprint's `pathA_39`/`j∝sV` scope is superseded — do NOT fold it.)

**Part VI — Knit** → author `…/part6_knit_atomic_split.md` → ~4–6 stages.
- Source math: `software/stage1_solver/reports/pathA_40_cone_lock.md`, `…/pathA_41_ng5_second_medium_drift.md`, `…/pathA_42_charge_coupled_scalar.md`
- Scripts (FULL reshape): `software/stage1_solver/tools/pathA_40_cone_lock_sympy.py` · `…/pathA_40_cone_lock.wl`; `…/pathA_41_ng5_second_medium_drift_sympy.py` · `…/pathA_41_ng5_second_medium_drift.wl`; `…/pathA_42_charge_coupled_scalar_sympy.py` · `…/pathA_42_charge_coupled_scalar.wl`
- Owed: **RE-ADJUDICATE** the cone-lock — `pathA_40`'s calibrated `λ_γ=1`/`c_E=c_γ` are no longer established by the G0 card / magnetism build (no committed lock); carry the NG5 `{ρ_B0,χ_c,C_hu}` reduction-incompleteness + the `pathA_42` `SCALAR_DEPARTURE_MAPPED` sim-gated risk.

**Part VII — Integration** → author `…/part7_integration_atomic_split.md` → ~5–8 GENUINELY-NEW stages (synthesis; full new-derivation gauntlet + GLM tertiary).
- Source math: assembles Parts I–VI (no single prior script). χ_Q reconciliation sources: `software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md` (numeric ≈0.712) vs `software/stage1_solver/tools/pathA_33_quadrupole_normalization_sympy.py` · `…/pathA_33_quadrupole_normalization.wl` (exact =1).
- Owed: author the unified equation set (conservative action + explicit non-variational drain/return sources/BCs), the calibration map, the whole-system dimensional firewall, the open-items/no-go/calibration registers, the simulation hand-off spec, the χ_Q reconciliation, and the full irreducible count via the `pathA_40` `Δr=2` codimension technique.

---

## 6. Known staleness / to-reconcile (don't be misled)

- **Decision 16 — its Part-I amendment has LANDED (2026-07-21).** `software/stage1_solver/decisions/16_retire_brane_polar_field.md` is **decided and operative** (it retires the brane polar field `P` — the "little arrows" — plus `λ_Pu` + `α_aniso`); the formal ledger amendment is now applied — stage006/007 audit scripts + `parameter_register.md` are amended (operative DOF=4, drift 7 `POST_D16_DRIFT(7)`, stage006 `DRIFT(5)`; `λ_Pu`+`α_aniso` retired; audit scripts pass — stage006 SymPy 121/Math 119, stage007 SymPy 142/Math 140). The historical freeze-as-run "11"/DOF=8 record STANDS immutable. The amendment (a definite **−5** route-less knob REDUCTION) marks the blueprint's `χ_B=|P_∥|²` route OBSOLETE-as-carried **unless a new T0 freeze is authorized** (it remains a live, high-risk route-(c) future gate in the knob audit — one whose collapse would clear ~half of the postulated-structure liability — not a dead end). Evidence: every `P` payoff failed independently (charge no-Gauss, light couple-stress no-go, wall falsified, and an active `INSTABILITY_CONFIRMED_STRUCTURAL` for any `λ_Pu≠0`).
- **Blueprint Parts IV/V scope is stale by design** — scoped on `pathA_38`/`pathA_39` before the EM reconsideration; `pathA_39` rests on the barred `j∝sV`. Re-scope onto the puncture-deflection + moving-throat builds.
- **Blueprint Part-VI cone-lock framing is stale** — it records pathA_40 as a settled `CONE_LOCK_CALIBRATED`; the newer builds don't establish the locks → it's a *to-re-adjudicate*.
- **χ_Q reconciliation** (`0.712` vs `1`) — live Part-VII debt.
- **G0 card is DRAFT v0** — un-reconciled with Part I; de-dup at Part VII, not additive-merge.
- **stage003 "Next step"** cites a superseded `c_L²=B_eff/ρ_br` cone-pair — don't carry verbatim.
- **The old 253-stage `research/pde_ledger/`** is a SUPERSEDED quarry — cite the v2 Part-II stages + pathA reports as the earned source, not it.

*(Resolved, no longer a to-reconcile: the stage007 `drift=11` count was contested but STANDS — no `ρ_br` overcount; see §3.6.)*

---

## 7. Glossary

- **Primitive constants** — the medium's ~4 declared-universal inputs: `ħ`, `m_GNLS`, `K`, `ρ0`. **`m_GNLS`** is the condensate boson mass — a *declared universal constant*, NOT emergent; do **not** confuse it with **`m_defect`** (a particle's mass), which is a GAP (`INFLOW_MASS_SOURCE_MISSING`, not derived — only a dimensional bridge `m∼α_J ħJ/c_γ²` exists).
- **Defect / throat / puncture** — a particle: a topological puncture of the brane into `±w`. Never "dent".
- **Brane** — the ordered (`χ_B=1`), shear-supporting phase; our 3D space, a domain wall at `w=0`. **Bulk** — the same medium de-structured (`χ_B=0`), shear-free.
- **±w** — the two sides the throat punctures into; the charge **sign** is this orientation (Z₂, not additive).
- **Drain** — a throat removing number-flux; gravity is the flow between drains.
- **Cone / speed** — `c_s` (density/gravity-change), `c_γ` (light, shared with magnetism), `c_E` (electric, no committed lock). `λ_γ=c_γ/c_s`.
- **DtN** — Dirichlet-to-Neumann map; the exterior-wave operator whose rational coefficients are the earned fingerprint. **χ_Q** — the outgoing-wave sign (`+1`).
- **R1_REQUIRED(x)** — an honest *unresolved* landing: the result needs the named closure `x` (usually the sim-deferred nonlinear throat solve). Not a failure, not a fudge — a scoped open item.
- **Departure** — a characterized, first-class break from the reference theory (GR or Maxwell); welcome, not a bug.
- **Atomic split** — the plan that decomposes a Part into individual buildable **stages**, approved at the per-Part gate before building.
- **Stage / 6-artifact unit** — the smallest verifiable derivation: note + TeX card + SymPy audit + independent Mathematica audit + source map + register entry.
- **Knit (Part VI)** — the consistency check that the four sectors are ONE medium (shared cones `pathA_40`, no smuggled second substance `pathA_41`, consistent scalar `pathA_42`).
- **Integration (Part VII)** — the unified equation set + calibration map + registers + sim hand-off spec — the literal deliverable.
- **Cone-lock** — a calibrated equality between two cones (`λ_γ=1`, `c_E=c_γ`); not derived, must be re-adjudicated.
- **GENUINE_BLOCKED** — a magnitude (e.g. `G`) that is a genuine external anchor, not derivable at this level.
- **Tri-review** — arbiter re-run + fidelity audit + adversarial-with-per-tooth-ablation, each on a fresh agent.
- **The G0 card** — the shared minimal closure card (`g0_closure_card_v0.md`); a DRAFT v0 postulated closure the EM builds reduce from.
