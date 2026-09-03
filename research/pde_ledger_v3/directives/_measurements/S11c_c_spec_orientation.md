# S11c-c spec — orientation record (what the substrate fixes before authoring)

Rule 2 binds the orchestrator: this records the source reads behind the S11c-c SHARED_PHYSICS spec, with the
files/lines they rest on. It is not a result; it is the provenance of the spec's SUPPLIED setup.

## The object S11c-c computes (from `directives/S11c_decisions.md:43`, N1/N2 table row S11c-c)
> **S11c-c curved-interface bulk closure**: the perturbed two-face outgoing bulk **DtN/impedance** (flux,
> traction, permeability/memory, face parity); the `v_bulk_normal_0` carry-or-restrict decision (`N11`) ⇒
> exports **coupled nonlocal self-energy/operator**.

## The seam — where S11c-a and S11c-b explicitly stopped and handed to S11c-c
- **S11c-a T-i** (`directives/S11c_a_SHARED_PHYSICS.md:448-453`, `S11CA_CLOSURE_SHAPE_DERIV`): shape-
  differentiates the face closure but **"This object is NOT B0c's bulk-response assembly `δp=Z·v_bulk`; no
  bulk DtN, impedance, or pressure-response solve belongs in T-i."** ⇒ the bulk solve is S11c-c's.
- **S11c-b §0** (`directives/S11c_b_SHARED_PHYSICS.md:25-27,35-39`): excludes "solution of the curved two-face
  outgoing bulk problem and its nonlocal DtN/impedance/self-energy; no bulk-response solve (`δp=Z·v_bulk`), no
  permeability/memory kernel, no face-parity impedance object" — the flat-face kernels `Λ_I(ω)=Λ_I⁰/(1−iωτ_I)`
  stay **symbolic** in every §3b face/flux contribution. ⇒ S11c-b's operator carries a **symbolic `δp_s`**; the
  θ-row = `evolution_mass_balance − Σ closure_shape_deriv` (`steps/S11c_b_variable_coefficient_operator.md:39`)
  where `closure_shape_deriv` carries `Λ_A𝒜_s+Λ_V V_s` and `Λ_X𝒜_s` with `δp_s` un-eliminated. S11c-c supplies
  the bulk relation that closes `δp_s`.

## The flat-face reference S11c-c generalizes — S11b B0/B0c (`directives/S11b_SHARED_PHYSICS.md`)
- **Bulk acoustics (§2, :167,170):** `v=∇₄φ`, `δp=−ρ_m∂_tφ`, `∂_t²φ=c_s0²∇₄²φ`; both half-spaces `|w|>W₀/2`
  are bulk (two independent half-spaces — the slab interior is not bulk).
- **Impedance def (§1, :92):** `Z ≡ (bulk pressure perturbation at a face)/(bulk material's OUTWARD normal
  velocity at that face)`.
- **Radiation condition (§1, :93):** in each half-space retain only outgoing (real normal wavenumber) or
  decaying (imaginary) waves; NO incoming from infinity; boundedness alone does not select a branch when the
  normal wavenumber is real. State the selected sign per half-space.
- **Branch object `q_out` (§1b, :116-152):** the root of §2's bulk dispersion selected by §1; branch points on
  the real axis at the **sound cone `ω=±c_s0|k|`**; complex-ω continuation preserves the branch; at `q_out=0`
  the two bulk solutions coalesce (degenerate/grazing).
- **B0b (§9, :605-618):** flat impedance `Z`; THREE regimes — bulk normal wavenumber² **positive / negative /
  zero(grazing)**; report all three incl. behaviour AT the zero; per-FACE inertial loading `m_add` in the
  reactive (purely-imaginary-Z) regime; ⚠ the per-face `m_add` acts with the SAME sign against each face's
  OUTWARD acceleration — the global-`w` signed pair `(X,−X)` is a convention artifact. ⇒ `S11B_Z_IMPERMEABLE`,
  `S11B_Z_BY_REGIME`, `S11B_Z_BY_PARITY`, `S11B_ADDED_MASS`, `S11B_GRAZING_BEHAVIOUR`.
- **B0c (§9, :620-636):** permeable face response — close `δp=Z·v_bulk,±` (bulk) + `v_bulk,±=V_±+J_±/ρ_m`
  (interfacial mass balance, kinematics) + the §4 closure `J_s=Λ_A𝒜_s+Λ_V V_s`, `𝒜_s=μ_s−δp_s/ρ_m`,
  `μ_s=μ_θ/ρ_br⁰`, `t_s=−(δp_s+Λ_X𝒜_s)n̂_s`; solve together → `δp` as a function of `V_±,μ_θ`; degenerate loci
  via the §8 locus protocol; dissipation audit per regime+parity. ⇒ `S11B_FACE_RESPONSE`,
  `S11B_FACE_RESPONSE_COEFFS`, `S11B_DEGENERATE_LOCI_*`, `S11B_PERMEABLE_DISSIPATIVE_BY_REGIME_AND_PARITY`.
- ⚠ **Two measured traps (§9, :760-816):** (1) a **real** part of `Z` in the propagating regime is **bulk
  radiation** to infinity, present even with impermeable faces — ⛔ do NOT read `Re Z` as transfer THROUGH the
  interface; (2) **branch re-selection** can turn a sink into a source / damped into apparent growth — report
  the SIGN of every imaginary part; a growing AND a decaying root are both reportable; ⛔ never suppress or
  re-branch. These are rule-17 live-quantity hazards for S11c-c.

## S11c-a substrate S11c-c consumes (the tilted-face shape derivatives, `S11CA_*`, per the S11c-a step record)
`S11CA_FACE_NORMAL` (T-a; `n̂_s^α=(−½∇W_bg,s)+O(|∇W_bg|²)`), `S11CA_CONORMAL_DERIV` (T-a′; `n̂_s·∇₄`, the
operator a DtN produces), `S11CA_FACE_MEASURE_SHAPE_DERIV` (T-a″; `a_s^α=√(1+|∇h_s|²)`), `S11CA_FACE_VELOCITY`
(T-b; `V_s^α`), `S11CA_RELATIVE_FLUX` (T-c; `J_s^α`), `S11CA_KINEMATIC_BALANCE` (T-c′; `n̂_s·v_bulk,s=V_s+J_s/ρ_m`),
`S11CA_TRACTION` (T-d), `S11CA_FACE_SHIFT` (T-e; `δ[f(x,h_s)]=δf(x,h_s⁰)+δh_s∂_wf⁰`), `S11CA_CLOSURE_SHAPE_DERIV`
(T-i, the seam). Faces tilt OPPOSITELY: `∂h_s/∂x = (s/2)∂W_bg` ⇒ the curvature correction has a definite parity
under `s→−s` — this is the "face parity" object, COMPUTED not stated.

## S11c-b substrate S11c-c consumes (`S11CB_*`)
`S11CB_SLAB_OPERATOR` (+`_TERM_ORIGINS`; the divergence-form variable-coefficient operator with symbolic
`δp_s`/`Λ` in the face slots), `S11CB_MU_THETA_OPERATOR` (`μ_θ=δU/δθ`, held-fixed constitutive operand feeding
`μ_s=μ_θ/ρ_br,bg⁰`), `S11CB_COUPLING_KERNEL` (the off-diagonal transverse→{θ,e_W,u_L} block). ⚠ IMPORTED
per-engine-verified, cross-engine residual DEFERRED (`steps/S11c_b_variable_coefficient_operator.md:15-31`) —
so S11c-c cross-checks its OWN objects; it does not re-open S11c-b.

## N-series constraints binding the S11c-c spec
- **N1** (`:24`): own spec, two blind engines, frozen T7 comparator, `S11c_c_exports.py`; SymPy chains from
  `S11c_b_exports.py` (transitively the S11b LEDGER + S11c-a), WL blind re-derives. `BUILD_INPUT_DIGESTS` pins
  {S11c-c audit, imported exports, S11c-c spec}.
- **N5** (`:75-84`): ⛔ NO global dispersion relation `ω(k)`; emit the DtN **operator/kernel**, not a spectrum;
  profile-conditioning is S11c-d. Expansion: `O(ε)` wave × first shape order in `η`/`σ_W`.
- **N8** (`:146-151`): chain wiring + T7 contract + the S11c-a reconciliation schema
  (`steps/S11c_a_interface_shape_derivatives.md:233-253`) inherited verbatim; ⛔ no pre-registered fold.
- **N11a** (`:171-178`): `v_bulk_normal_0` INHERITED as a standing rest-frame limit; every S11c-c object
  conditional on the derived smallness domain `|q·v_bulk_normal_0/ω|≪1` (`steps/S11b_interface_coupling_law.md:159-161`;
  large `k c_s0/|ω|` necessary NOT sufficient). ⛔ Do NOT make the convective bulk operator an S11c task — the
  rest-frame PDE `∂_t²φ=c_s0²∇₄²φ` stands; `v_bulk_normal_0` never aliased to `v_0`.
- **N12** (`:115-120`): every object multigraded `(ε,η,σ_W)`; bulk constants `ρ_m,c_s0` do NOT vary in-plane
  (`directives/S11c_b_SHARED_PHYSICS.md:202`) — the nonlocality comes from the CURVED BOUNDARY, not a varying bulk.
- **N13/N7**: confinement (`N13`) + falsification/bench-optics (`N7`) are **S11c-e**, NOT here.
- **N14** (`:127-133`): fresh `S11CC_*` names; ⛔ never reuse imported keys `W_0,mu_R,rho_br,e_W,v_0` for a
  varying object.

## Process (N1/N2; matches S11c-a/b — no per-sub-step decision list)
Family decision list (`S11c_decisions.md`) is done. S11c-c's FIRST artifact = `directives/S11c_c_SHARED_PHYSICS.md`,
orchestrator-written ⇒ **2 decision legs (Codex + Grok), rule 7 TRIGGER**, folded ONCE, then the dual-engine
build. There was NO `S11c_a_decisions.md`/`S11c_b_decisions.md`; the SHARED_PHYSICS spec is the physics-bearing
artifact that carries the two legs.
