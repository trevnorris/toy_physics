# Citation Map

Reusable citation catalog for the moving-throat PDE archive. Organized by the
stable `MTDC-T*` public anchors. Future papers (framework, PN, retarded,
quotient, realization, material) should cite the anchor name first and use
stage numbers only for detailed audit routing.

See also:
- `CITATION_SUPPORT_SET.md` — rationale for the 25-stage support set
- `CHECKPOINT_TRUST_AUDIT.md` — current trust tier per checkpoint
- `STAGE_PROVENANCE_INDEX.md` — raw note and audit artifact paths
- `STAGE_VERIFICATION_COVERAGE.md` — repo-wide coverage control sheet
- `MATHEMATICA_MIRROR_POLICY.md` — which Mma audits count as independent vs replay

Snapshot date: `2026-04-21`

## How To Use This File

1. Find the anchor family your paper needs (below).
2. Read the one-line claim to confirm it's the thing you want to cite.
3. Use the recommended citation form in your bibtex / `\cite{}` stem.
4. If a referee asks "where is this verified?", point them at the audit paths.
5. If your paper needs a *quantitative* statement rather than the symbolic
   theorem, check the `Numerical` row — not every anchor has one.
6. Respect scope caveats. Do not cite an anchor outside its declared status.

## Quick-Scan: Safe To Cite

All 28 entries below have trustworthy audit support, explicit constant
provenance, and a human review note. Most are citation-grade symbolic anchors;
entries marked `Open` or `Numerical` must be cited only within the narrower
scope stated in their caveats.

| Anchor | Stage(s) | Trust | Numerical | Short claim |
|---|---|---|---|---|
| MTDC-T1 | 001 | strong | — | Geometry lift of breathing + modal wall + linearized Maxwell. |
| MTDC-T2 | 002 | strong | — | Breathing reduction to `(a, L, P_2)` conservative modules. |
| MTDC-T3 | 003 | strong | yes | First microscopic BdG support kernel; exact pole-shift packet. |
| MTDC-T4 | 005 | strong | — | Grouped-`P_2` outgoing-normalization bridge; invariant response product. |
| MTDC-T5 | 006, 007, 019 | strong | — | Full grouped bundle, overlap isotropy law, support-feasibility frontier. |
| MTDC-T6.4 | 032, 034 | strong | — | D/N overlap primitive and lowest-twin sufficiency criterion. |
| MTDC-T7.2 | 052 | strong | — | Reduced support/source three-zone verdict. |
| MTDC-T7.7 | 072, 073 | strong | — | Family-1 minimal isotropic verdict and theorem-status boundary. |
| MTDC-T8 (class) | 100 | strong | — | Compensated Robin-mixed classification endpoint (Stages 90-100). |
| MTDC-T8 (renorm) | 140 | strong | yes | Renormalized Family-1 co-evolving branch capstone (Numerical/Open; Stages 123-140). |
| MTDC-T8.1 | 079 | strong | — | Geometry-lane verdict (conservative 3/4 + 1/4 module). |
| MTDC-T8.2 | 088 | strong | — | Exact fixing of `χ_Q` via outgoing-DtN match. |
| MTDC-T8.3 | 095 | strong | — | Hybrid/Robin compensation law (`ρ_R=4σ_W`, `κ_W=1/3`, `γ_W=1/9`). |
| MTDC-T8.5 | 146 | strong | — | First explicit off-family drift scalar `δ_⊥` with microscopic formula. |
| MTDC-T9.3 | 168 | strong | yes | Microscopic monomial coordinates `(C_tr,*, C_nt,*, ε_η)`. |
| MTDC-T9.6 | 183 | strong | — | Four-scalar final verdict packet `Δ_full = (Δ_Q, q_tr, q_nt, q_η)`. |
| MTDC-T10.1 | 186 | strong | — | Scalar graph-slice theorem; IVT crossing for `χ_Q(y)=1`. |
| MTDC-T10.4 | 201 | strong | — | Final local mixed-ray closure; support-≤5 splice interval. |
| MTDC-T11.1 | 204 | strong | — | Resonance/linewidth dispersive no-free-lunch theorem. |
| MTDC-T11.6 | 222 | strong | — | Rigid-mouth physical `(U,V)` chart; Cartesian orbit-lock theorem. |
| MTDC-T11.7 | 225 | strong | — | Twin-support placement and coherent orbit-lock compiler. |
| MTDC-T12.1 | 226 | strong | — | Relaxed recovery map; short-range open-system firewall. |
| MTDC-T12.3 | 231 | strong | yes | Dynamic event-chain compiler; turning-point / threshold-speed / WKB. |
| MTDC-T12.5 | 236 | strong | yes | Physical calibration and material-threshold export packet. |

---

## Foundational Anchors (framework papers lean here first)

### MTDC-T1 — Geometry lift

**Claim**: The confinement potential perturbation `δV_conf` has a linear
functional form in `η(Ω, w, t)`; modal wall PDE separates harmonically with
`ℓ(ℓ+1)` angular eigenvalues; localized Maxwell equation carries a densitized
`Z(w)` weight.

**Backing stage**: `001`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl`

**Scope caveats**: The matter-side source `S_{lm}^{(ψ,A)}` is treated as an
abstract symbolic bilinear coupling; its microscopic origin is deferred to
Stages 028-031.

**Supports paper families**: framework, PN, retarded, quotient, realization, material.

**Recommended citation stem**: `MTDC-T1 (stage 001)`

---

### MTDC-T2 — Breathing reduction

**Claim**: Matrix Euler-Lagrange equations `M_{AB}Q̈^B + K_{AB}Q^B = 0` with
explicit overlap matrices `M_{AB} = 4π ∫ μ_η α_A α_B dw`,
`K_{AB} = 4π ∫ [T_w α'_A α'_B + K_0 α_A α_B] dw`. Real `P_2` lane carries the
grouped module with stiffness shift `K_η + 6 T_Ω`.

**Backing stage**: `002`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage002_breathing_reduction_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl`

**Scope caveats**: Mathematica basis restricted to `{Y_20, Y_21c, Y_22c}`; the
other two real `P_2` components are closed by an explicit phase-shift identity
rather than a 5×5 orthogonality integral.

**Supports paper families**: framework, PN, realization.

**Recommended citation stem**: `MTDC-T2 (stage 002)`

---

### MTDC-T3 — BdG pole-shift packet

**Claim**: Exact Schur-complement kernel `D_0^eff(ω) = K_0 - ω²M_0 - C(Ω_0² - ω²I)^{-1}C^T`;
closed-form two-pole solution `ω_±² = [Ω_η² + ϖ² ± √((Ω_η² - ϖ²)² + 4g²/M)]/2`;
perturbative shifts `δΩ_η² = -g²/[M(ϖ² - Ω_η²)] + O(g⁴)` and analogue for the
matter side.

**Backing stage**: `003`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage003_bdg_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl`
- Numerical: `scripts/numerical/stage003_004_foundational_stress.py`
- Numerical (Mma): `mathematica/numerical/stage003_004_foundational_stress.wl`

**Scope caveats**: The grouped-`P_2` kernel is audited with a single matter
mode per lane rather than an explicit `Σ_α` sum. Linear structure makes
extension to multi-mode trivial, but scripts do not spell it out.

**Supports paper families**: framework, PN.

**Recommended citation stem**: `MTDC-T3 (stage 003)`

---

### MTDC-T4 — Grouped-`P_2` normalization bridge

**Claim**: Outgoing fingerprint `A = a²/(9c_s²)`, `B = 4a⁴/(81c_s⁴)`,
`Γ_5 = a⁵/(27 c_s⁵)` (rebuilt from spherical Hankel `h_2^(1)` via DtN series).
Invariant normalization product
`m̂_0² N_0/(K - B_0 - Z_0) = 54 G c_s⁵ / (5 a⁵ c⁵)`.

**Backing stage**: `005`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage005_grouped_p2_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage005_grouped_p2_normalization_bridge_mathematica_audit.wl`

**Scope caveats**: None open. All six boxed identities verified exactly.

**Supports paper families**: framework, PN.

**Recommended citation stem**: `MTDC-T4 (stage 005)`

---

### MTDC-T5 — Full grouped bundle, overlap isotropy, support-feasibility frontier

**Claim (bundle)**: Exact projectors `P_x̄, P_a, P_b` with completeness and
orthogonality; coefficient decomposition `D_{An} = K_A - B_{An} - Z_{An}` etc.
Isotropic normalization product reconstructed.

**Claim (isotropy)**: Real-STF orthonormality `∫Y_A Y_B dΩ = δ_{AB}`;
axisymmetric overlap matrix `M^(20) = (√5/(7√π)) diag(1, 1/2, 1/2, -1, -1)`;
grouped splitting signature `(λ_{20}, λ_{21}, λ_{22}) = (1, 1/2, -1)`; defect
ratio `b = 3a`.

**Claim (frontier)**: Support feasibility condition `M_mix ≤ G(ξ_req, δ)` with
closed-form `G = 9ξ(ξ+δ)/(9δ+11ξ)`. Three-conjunct admissibility test
`R_target ≥ 1 ∧ F(ξ_req, δ) = R_target ∧ M_mix ≤ G(ξ_req, δ)`.

**Backing stages**: `006, 007, 019`

**Audit paths**:
- SymPy (006): `scripts/moving_throat_pde_stage006_full_grouped_bundle_sympy_audit.py`
- Mathematica (006): `mathematica/moving_throat_pde_stage006_full_grouped_bundle_mathematica_audit.wl`
- SymPy (007): `scripts/moving_throat_pde_stage007_overlap_isotropy_sympy_audit.py`
- Mathematica (007): `mathematica/moving_throat_pde_stage007_overlap_isotropy_mathematica_audit.wl`
- SymPy (019): `scripts/moving_throat_pde_stage019_support_feasibility_frontier_sympy_audit.py`
- Mathematica (019): `mathematica/moving_throat_pde_stage019_support_feasibility_frontier_mathematica_audit.wl`

**Scope caveats**: Stage 006 carries `B_{An}` as symbols anchored to Stage 003's
Schur derivation rather than re-deriving the sum at Stage 006. Not a defect
if cited together with `MTDC-T3`.

**Supports paper families**: framework, realization, quotient.

**Recommended citation stems**:
- `MTDC-T5 (bundle, stage 006)`
- `MTDC-T5 (isotropy, stage 007)`
- `MTDC-T5 (frontier, stage 019)`

---

## MTDC-T6 — Lowest-twin sufficiency and D/N overlap primitive

### MTDC-T6.4

**Claim (overlap)**: From the microscopic D/N mode
`χ_n(s) = √(2/L) sin((n+½)πs/L)`, overlap ratio `I_n/I_0 = 1/(2n+1)`; twin
formula `ζ_n^twin = 1/((2n+1)²(1 + x·n(n+1)))` where
`x = π² T_X/(L² K_W^eff)`; lowest-twin value `ζ_0^twin = 1`.

**Claim (sufficiency)**: `Π_tr ≤ 2 C_mix = 16Λ(1-ε)/π²` iff the lowest twin
supplies demand (`S_0 = 2`).

**Backing stages**: `032, 034`

**Audit paths**:
- SymPy (032): `scripts/moving_throat_pde_stage032_dn_overlap_zeta_sympy_audit.py`
- Mathematica (032): `mathematica/moving_throat_pde_stage032_dn_overlap_zeta_mathematica_audit.wl`
- SymPy (034): `scripts/moving_throat_pde_stage034_lowest_twin_criterion_sympy_audit.py`
- Mathematica (034): `mathematica/moving_throat_pde_stage034_lowest_twin_criterion_mathematica_audit.wl`

**Scope caveats**: Stage 032 is a Wave 5 hardening addition (previously lacked
executable audit). Now carries the microscopic overlap derivation; Stage 033
and downstream consumers import `twin_support_ratio` rather than re-declaring.

**Supports paper families**: framework, realization.

**Recommended citation stems**:
- `MTDC-T6.4 (D/N overlap primitive, stage 032)`
- `MTDC-T6.4 (lowest-twin sufficiency, stage 034)`

---

## MTDC-T7 — Reduced support/source verdict family

### MTDC-T7.2 — Reduced support/source three-zone verdict

**Claim**: Universal fail `W_wall < Pe_req/Δ_∞`; universal success
`W_wall > Pe_req/Δ_0`; profile-sensitive band
`Pe_req/Δ_∞ ≤ W_wall ≤ Pe_req/Δ_0`.

**Backing stage**: `052`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage052_final_reduced_verdict_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage052_final_reduced_verdict_mathematica_audit.wl`

**Scope caveats**: Script verifies `Δ_0 < Δ_∞` and `C_res² < 1` positivity from
carried regime assumptions. Prior tautology concern is closed.

**Supports paper families**: framework, realization.

**Recommended citation stem**: `MTDC-T7.2 (stage 052)`

---

### MTDC-T7.7 — Family-1 verdict and theorem-status boundary

**Claim (verdict)**: Minimal isotropic demand `ζ_req^min = 1/3 < A_F1 ≈ 1.00005`
where `A_F1 = (κ_F1 + π²/4)/(κ_F1 + y_F1²)` with `κ_F1 = 12321/5` and
`y_F1 tan(y_F1) = 37`. Consequence: `Pe_req = 0`.

**Claim (boundary)**: Natural contact-plus-pole reading gives `ρ_α = 4/3`,
`ζ_req = 1/3`; this lies strictly inside the Family-1 success region; the
explicit branch succeeds at zero transport bias.

**Backing stages**: `072, 073`

**Audit paths**:
- SymPy (072): `scripts/moving_throat_pde_stage072_family1_minimal_isotropic_verdict_sympy_audit.py`
- Mathematica (072): `mathematica/moving_throat_pde_stage072_family1_minimal_isotropic_verdict_mathematica_audit.wl`
- SymPy (073): `scripts/moving_throat_pde_stage073_updated_reduced_status_sympy_audit.py`
- Mathematica (073): `mathematica/moving_throat_pde_stage073_updated_reduced_status_mathematica_audit.wl`

**Scope caveats**: Status boundary is narrow — claim is valid under the
minimal-module hypothesis. Do not cite as a proof of the minimal module itself.

**Supports paper families**: framework, realization.

**Recommended citation stems**:
- `MTDC-T7.7 (verdict, stage 072)`
- `MTDC-T7.7 (status boundary, stage 073)`

---

## MTDC-T8 — Geometry lane and realization class

### MTDC-T8 (classification capstone)

**Claim**: Admissible low-frequency outlet deformations reduce to the
compensated Robin-mixed class. Five candidate classes (pure-scale,
pure-argument, Robin-only, mixed-only, Robin+mixed); only Robin+mixed
preserves both even-coefficient and odd-normalization constraints on the
core-balance surface with the D/N tube normalization
`L_W = πa√((1+r_c)/3)/2`.

**Backing stage**: `100` (capstone of Stages 90-100)

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage100_outlet_core_status_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage100_outlet_core_status_mathematica_audit.wl`

**Scope caveats**: Status `ExactClosure / Open`. The microscopic question —
whether the real core actually lands on the compensation surface — remains
open; the audit closes the classification theorem under that conditional. Do
not cite as a proof of realization.

**Supports paper families**: retarded, quotient, realization.

**Recommended citation stem**: `MTDC-T8 classification (stage 100)`

---

### MTDC-T8 (selection + renormalization capstone)

**Claim**: Regular Family-1 co-evolving branch selected; renormalized
canonical tuple `(Σ_0^can, T̂^can, S^can, Π^can)` with exact compensation
`R_can = 1/4`. Frozen-traction contrast `R(σ_fp) > 1/4` strictly — proof
that renormalization is *needed*, not just consistent. Bisection root
uniqueness via monotonicity scan.

**Backing stage**: `140` (capstone of Stages 123-140)

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage140_core_mouth_coevolution_status_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage140_core_mouth_coevolution_status_mathematica_audit.wl`
- Numerical: `scripts/numerical/stage140_core_mouth_coevolution_status_stress.py`

**Scope caveats**: Status `Numerical / Open`. The canonical tuple is computed
numerically via bisection on `[3, 6]`; CAS audits import the result. Frozen
traction contrast lives in the numerical harness. Do not cite as a full
symbolic closure.

**Supports paper families**: retarded, realization.

**Recommended citation stem**: `MTDC-T8 selection/renormalization (stage 140)`

---

### MTDC-T8.1 — Geometry-lane verdict

**Claim**: Isotropic ℓ=0 ↔ ℓ=2 decoupling gives the conservative module
`Ŷ_Q^cons = 3/4 + (1/4)/(1 - ω²/Ω_Q²)`; derived `ρ_α = 4/3`, `ζ_req = 1/3`.

**Backing stage**: `079`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage079_geometry_lane_check_verdict_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage079_geometry_lane_check_verdict_mathematica_audit.wl`

**Scope caveats**: Parity between CAS engines restored post-Wave-4; both
engines now verify ℓ=0/ℓ=2 orthogonality and the frequency-dependent carrier.

**Supports paper families**: framework, realization.

**Recommended citation stem**: `MTDC-T8.1 (stage 079)`

---

### MTDC-T8.2 — Fixing `χ_Q` via outgoing-DtN

**Claim**: Matching the canonical grouped module to the exact outgoing DtN
branch gives `χ_Q = 1` (exact, not perturbative).

**Backing stage**: `088`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage088_chiQ_fix_from_outgoing_dtn_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage088_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl`

**Scope caveats**: None open. Normalized z-expansion verified against
`(z²/9, 4z⁴/81, i z⁵ χ_Q/27)`.

**Supports paper families**: framework, PN.

**Recommended citation stem**: `MTDC-T8.2 (stage 088)`

---

### MTDC-T8.3 — Hybrid/Robin compensation law

**Claim**: Nontrivial compensated branch has `ρ_R = 4σ_W`, `κ_W = 1/3`,
and preserves odd normalization iff `γ_W = 1/9`.

**Backing stage**: `095`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage095_hybrid_robin_mixed_compensation_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage095_hybrid_robin_mixed_compensation_mathematica_audit.wl`

**Scope caveats**: None open. Both CAS engines run `sp.solve` / `Solve[..., Reals]`
independently and return the same two branches.

**Supports paper families**: realization, retarded.

**Recommended citation stem**: `MTDC-T8.3 (stage 095)`

---

### MTDC-T8.5 — First explicit off-family drift

**Claim**: Scalar `δ_⊥ = δg - g'(r)·δr` is the first true off-family drift.
Microscopic parent-action formula available; even-preservation and
tangent-motion conditions close to zero.

**Backing stage**: `146`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage146_off_family_normal_coordinate_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage146_off_family_normal_coordinate_mathematica_audit.wl`

**Scope caveats**: Stage cards scripts also verify outlet-defect transport
`(δC, δE_2, δE_4, Δ_Q)` and mouth-bias split `δΠ`; strictly stronger than the
paper's minimal prose.

**Supports paper families**: PN, realization.

**Recommended citation stem**: `MTDC-T8.5 (stage 146)`

---

## MTDC-T9 — Microscopic monomials and final verdict packet

### MTDC-T9.3 — Microscopic monomial coordinates

**Claim**: Branch composites push to direct monomials
`(C_tr,*, C_nt,*, ε_η)` with drift laws `δ ln C_tr,* = Σ_tr`,
`δ ln C_nt,* = Σ_nt`, `δ ln ε_η = Σ_η`. Exponent matrix `M_*` with
dependent-minor determinant `1 + χ_0,*`.

**Backing stage**: `168`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage168_microscopic_monomials_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage168_microscopic_monomials_mathematica_audit.wl`
- Numerical: `scripts/numerical/stage168_170_orbit_stress.py`
- Numerical (Mma): `mathematica/numerical/stage168_170_orbit_stress.wl`

**Scope caveats**: Paper's `δ_U := π² T_U/(L² K_U)` vs scripts' `T_U/K_U` is a
reparametrization; absorbed constants have zero logarithmic drift, so drift
laws match exactly.

**Supports paper families**: framework, PN, realization.

**Recommended citation stem**: `MTDC-T9.3 (stage 168)`

---

### MTDC-T9.6 — Four-scalar final verdict packet

**Claim**: Packet `Δ_full = (Δ_Q, q_tr, q_nt, q_η)` combines Packet A and
Packet B in additive, multiplicative, and mismatch charts. Cocycle law and
finite similarity-orbit lock hold; final zero packet is a *reduced compiler
condition*, not a proof of PDE realization.

**Backing stage**: `183`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage183_reference_free_home_stretch_theorem_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage183_reference_free_home_stretch_theorem_mathematica_audit.wl`

**Scope caveats**: Compiler condition, not realization. Do not overclaim.

**Supports paper families**: PN, realization.

**Recommended citation stem**: `MTDC-T9.6 (stage 183)`

---

## MTDC-T10 — Scalar graph-slice and mixed-ray closure

### MTDC-T10.1 — Scalar graph-slice theorem / IVT crossing

**Claim**: Graph slice `Z_* = { x_*^graph(y) : Δ̂_Q(y) = 0 }`; graph-tangent
kernel `M_* δx_graph = 0` (Packet B vanishes on graph); one-parameter IVT
crossing theorem with sign-change witness and unique interior root.

**Backing stage**: `186`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage186_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage186_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.wl`

**Scope caveats**: IVT is exhibited as a worked witness on an explicit path;
general real-analysis IVT theorem is not re-proven (it is a standard result).

**Supports paper families**: realization.

**Recommended citation stem**: `MTDC-T10.1 (stage 186)`

---

### MTDC-T10.4 — Final local mixed-ray closure

**Claim**: Boundary-identification theorem `∂Δ_5^+ = S_{≤4}^loc`; support
cardinality ceiling = 5; support-≤5 splice interval with explicit
three-family classification; preferred total budget 1464, fallback 2640.

**Backing stage**: `201`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage201_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage201_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl`

**Scope caveats**: None open.

**Supports paper families**: realization.

**Recommended citation stem**: `MTDC-T10.4 (stage 201)`

---

## MTDC-T11 — Resonance, physical chart, and placement

### MTDC-T11.1 — Resonance/linewidth dispersive no-free-lunch theorem

**Claim**: `|Re χ_*|/|Im χ_*| = |δ|/γ_*` — resonant conservative gain is
inseparable from absorptive loading; maximal balanced leverage at equal
conservative and absorptive magnitudes. Breit-Wigner reduction
`χ(ω) = A_*/(δ - iγ_*)`; low-loss survival window.

**Backing stage**: `204`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage204_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage204_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl`

**Scope caveats**: Numerical slice in script is explicitly probe-only.

**Supports paper families**: realization, material.

**Recommended citation stem**: `MTDC-T11.1 (stage 204)`

---

### MTDC-T11.6 — Rigid-mouth physical `(U, V)` chart

**Claim**: Physical log-variables `U = ln(T²/T_ref²)`, `V = ln(ε_η/ε_{η,ref})`
give a diagonal rigid-mouth packet. Dependent-plane compiler
`(Δ_T, Δ_{K_η}, Δ_μ) = (0, -V, U - V)`. Cartesian orbit-lock theorem: lock is
exactly `U = V = 0`.

**Backing stage**: `222`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage222_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage222_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.wl`

**Scope caveats**: Support-blindness driven through the Stage-221 carried
branch formula (not by symbol-free differentiation).

**Supports paper families**: realization, retarded, material.

**Recommended citation stem**: `MTDC-T11.6 (stage 222)`

---

### MTDC-T11.7 — Twin-support placement and coherent orbit-lock

**Claim**: Placement coordinate `ϱ_phys = 2(1-ε)/3`; support-lane
classification (selected branch above mixed-only, below non-twin). Finite
orbit packet with support/orbit split; two-packet compiler (direct + orbit)
invertible on the declared branch.

**Backing stage**: `225`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage225_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage225_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl`

**Scope caveats**: Rational sample point is probe-only. `N_Q-1` finish-line
datum is a labeling note, not an algebraic claim.

**Supports paper families**: realization, material.

**Recommended citation stem**: `MTDC-T11.7 (stage 225)`

---

## MTDC-T12 — Realization, dynamics, material

### MTDC-T12.1 — Relaxed recovery map and short-range firewall

**Claim**: Gaussian leakage/work lane values `S_leak = -√2 ℓ_w j_0/4`,
`W_w = √(2π) ℓ_w j_0 E_0/8`. Non-rigid stationarity
`U = k_V f_U/(k_U k_V - χ²)`, etc. Short-range one-port correction kernel
`δV_stat(r) = -½[C_6/r⁶ + 2 C_4 e^(-2κr)/r⁴ + C_2 e^(-4κr)/r²]`; firewall
`lim_{r→∞} r·δV_stat = 0` — relaxed branch is a short-range/open-system
bypass, NOT a new long-range same-charge force.

**Backing stage**: `226`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage226_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage226_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.wl`

**Scope caveats**: Short-range kernel rebuilt from declared primitives
`S_Q(r) = r^{-3}`, `S_Y(r) = e^{-2κr}/r`. Do not cite as a
new long-range force.

**Supports paper families**: realization, material.

**Recommended citation stem**: `MTDC-T12.1 (stage 226)`

---

### MTDC-T12.3 — Dynamic event-chain compiler

**Claim**: Exact energy conservation for reduced radial EOM; finite-radius
threshold law `v_crit,new = √(2(V_peak - V_0)/m_s)` and Coulomb contact
threshold `v_contact,Coul`; turning-point transport `dr_±/dE = 1/V'(r_±(E))`;
near-top parabolic action `I_top = π(V_peak - E)√(m_s/K_peak)/ℏ_eff + O(ΔE^{3/2})`.

**Backing stage**: `231`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage231_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage231_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`
- Numerical: `scripts/numerical/stage231_event_chain_stress.py`
- Numerical (Mma): `mathematica/numerical/stage231_event_chain_stress.wl`

**Scope caveats**: Benchmark numerics in CAS scripts are labeled; real
quantitative cross-checks live in the numerical harness.

**Supports paper families**: realization, material.

**Recommended citation stem**: `MTDC-T12.3 (stage 231)`

---

### MTDC-T12.5 — Physical calibration and material-threshold export packet

**Claim**: Physical dictionary `t^phys = t_* t`, `r^phys = (λ_phys/λ_ref) r`,
`E^phys = E_* E`. Safe-edge lattice event-equivalent rate
`γ_lat,safe^eq = f_lat(s_c) μ_η (s_0² - s_c²)/s_c`. Four screening ratios
`Π_ep, Π_χ, Π_k, Π_T` with `Π_* ≥ 1` admissibility stack; Korringa ceiling
`T_max = K_corr/t_cross^phys`.

**Backing stage**: `236`

**Audit paths**:
- SymPy: `scripts/moving_throat_pde_stage236_physical_calibration_and_material_threshold_companion_from_the_stage235_export_and_cold_survival_compiler_sympy_audit.py`
- Mathematica: `mathematica/moving_throat_pde_stage236_physical_calibration_and_material_threshold_companion_from_the_stage235_export_and_cold_survival_compiler_mathematica_audit.wl`
- Numerical: `scripts/numerical/stage236_material_threshold_stress.py`
- Numerical (Mma): `mathematica/numerical/stage236_material_threshold_stress.wl`

**Scope caveats**: Numerical stress includes genuine per-ratio fail cases
(`ep_fail_host`, `thermal_fail_host`, `geometry_fail_host`,
`stiffness_fail_host`). Screening ratios are the appropriate quantitative
lever.

**Supports paper families**: material, realization.

**Recommended citation stem**: `MTDC-T12.5 (stage 236)`

---

## Paper-Family Citation Menus

Quick reference: which anchors each downstream paper type is most likely to cite.

### Framework paper (parent ontology, geometry, linearization)
- MTDC-T1, T2 (geometry lift, breathing reduction) — foundational
- MTDC-T3 (BdG pole-shift packet) — first microscopic kernel
- MTDC-T4, T5 (normalization bridge, bundle, isotropy) — response machinery
- MTDC-T8.2 (`χ_Q = 1` fix) — outgoing-DtN anchor

### PN paper (wave response, far field, leading-order radiation)
- MTDC-T3, T4 (response packet and invariant normalization)
- MTDC-T8.2 (`χ_Q = 1`)
- MTDC-T8.5 (off-family drift)
- MTDC-T9.3, T9.6 (monomial coordinates, four-scalar packet)

### Retarded / mouth-branch paper
- MTDC-T8 classification (stage 100) — outlet deformation class
- MTDC-T8.3 (compensation law, stage 095)
- MTDC-T8 renormalization (stage 140) — co-evolving branch capstone
- MTDC-T11.6 (physical `(U, V)` chart, stage 222)

### Quotient paper (invariance, orbit structure)
- MTDC-T5 (isotropy and splitting signature)
- MTDC-T9.6 (four-scalar packet with orbit lock)
- MTDC-T11.6, T11.7 (physical chart and orbit-lock compilers)

### Realization / same-charge paper
- MTDC-T10.1, T10.4 (graph-slice, mixed-ray closure)
- MTDC-T11.* (resonance gate, chart, placement)
- MTDC-T12.1 (relaxed recovery + short-range firewall — **do not let this read as a new long-range force**)

### Material / experimental-contact paper
- MTDC-T12.3 (dynamic event chain, WKB, thresholds)
- MTDC-T12.5 (screening ratios, Korringa ceiling, material calibration)
- MTDC-T11.1 (resonance/linewidth tradeoff)

---

## Do Not Cite Yet

Stages that are in the archive but should **not** be cited as load-bearing
support in a paper unless they are hardened further.

### Stages with no executable audit (11 canonical stages)
`086, 096, 103, 107, 111, 115, 119, 124, 128, 132, 136`

**Policy**: These are paper-unreferenced (confirmed by a paper-grep
triage). They exist as appendix self-contained ledger entries and do not
appear in the citation graph of the canonical paper. Do not cite in new
papers unless hardening is added first. This count follows the canonical
`001--236` stage baseline in `STAGE_VERIFICATION_COVERAGE.md`; any
non-canonical status-only cards sit outside that numbered-stage count.

### Mathematica-only outliers (2 stages)
`067, 076` — have Mma audit but no SymPy mirror.

**Policy**: Not cited anywhere load-bearing. Coverage gap is cosmetic. Do
not introduce them as citation targets; they are internal consistency checks
only.

### SymPy-only late-stage region
`171--182, 184--185, 187--200, 202--203, 205--221, 223--224, 227--230, 232--235`

**Policy**: These are SymPy-verified derivation stages between the named
checkpoint anchors. They support the chain but are not intended as direct
citation handles. Future papers should cite the adjacent `MTDC-T*` anchor,
not the intermediate stage number.

### Status tags to respect
- `\StatusExactClosure` — safe to cite within declared scope.
- `\StatusOpen` — do not cite as a closed result; the stage is declaring an
  open branch. Affected anchors: **MTDC-T8 classification (stage 100)** and
  **MTDC-T8 renormalization (stage 140)** both carry `/ Open` — scope to
  the conditional they certify.
- `\StatusNumerical` — citation-ready only for quantitative bounds, not for
  symbolic closure statements. Affected anchors: MTDC-T8 renormalization
  (stage 140), some Part VII stages (215 etc.).

---

## Update Policy

When adding a new paper or new hardening round:

1. If a new checkpoint is added to the support set, add a new anchor entry
   here (in anchor-number order) with all six rows (claim, paths, caveats,
   paper families, citation stem).
2. If a status tag changes (`Open` → `ExactClosure`), update the **Scope
   caveats** line for the affected entry.
3. If a Wave-5-style hardening adds a previously-unaudited stage, add it as
   a new entry or append it to the existing anchor entry.
4. Keep the **Quick-Scan Safe To Cite** table in sync with the detailed
   entries — if it disagrees with the body, fix the table.
5. Keep the paper-family menus as the public interface for co-authors
   choosing what to cite. Most people will read only those menus plus the
   quick-scan table.
