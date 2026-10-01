# S11c_d_zinvariant_operator_blocks_directive.md (v4) — grounding commands (rule 2 / E1)

Mechanical lookups, run from the repo root on 2026-09-30, HEAD `c1e96e76`. Regenerated from the commands below;
nothing transcribed. Fences are longer than fences in quoted source text.

````
$ sed -n 8,12p AGENTS.md
For future S11c Python constructors, validators and export jobs, use
`scripts/s11c_guarded_run.py` around the existing supervisor. It requires a
host systemd user manager, verifies a whole-job 2 GiB memory cap, disables
job swap, limits process count and CPU affinity, lowers CPU/I/O priority,
records resource samples and fails closed. Keep native memory limits too. Never fall back to an unguarded launch if containment fails.
````

````
$ ls -la scripts/s11c_guarded_run.py
-rw-rw-r-- 1 trevnorris trevnorris 10714 Sep 30 09:51 scripts/s11c_guarded_run.py
````

````
$ sed -n 5p research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md
On the user's request to assess before continuing and stop for go/no-go, the existing central-balance-v2 scientific child was suspended in memory, not terminated or restarted. The pause receipt records PID 4097233, its process start identity, cgroup and exact command. It had completed two matrix groups after twenty restorations and saved 12,288 batches of the next group; no finite solve had been reached. The guard/supervisor remain active, no deadline was added, and pinned sources are unchanged. Resume requires the user's decision and rechecking that process identity. Do not launch another job alongside it.
````

````
$ sed -n 59,71p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
The slab degrees of freedom are S11b's (`S11b_SHARED_PHYSICS.md:69–80`), with the internal slab fields
`{u,δW,θ}` distinguished from the two independent face variables `{ζ_+,ζ_-}` (S11c-a §3a):

```text
u(x,t)     in-plane displacement, three in-plane components, no w-component ;
θ(x,t)     Eulerian densification, ρ_4D = ρ_4D⁰(1+θ) ;
ζ_+, ζ_-   the two independent face variables, combined as
           δW ≡ ζ_+ − ζ_- (thickness) ,   ζ_c ≡ (ζ_+ + ζ_-)/2 (centre shift) ,
           ζ_s = ζ_c + s δW/2 ,   e_W ≡ δW/W₀ .
```

⛔ `ζ_c` is an independent face DOF; no S11c-b computation may set `ζ_c=0` or replace the two face variables
by a thickness-only ansatz (S11c-a §3a), except in an explicitly named centre-fixed uniform regression.
````

````
$ sed -n 90,91p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
`v_bulk_normal_0` (the bulk normal drain, `S11b_SHARED_PHYSICS.md:104`) is a scope-limit parameter, not an
active DOF, and appears in no derived operator (§0).
````

````
$ sed -n 95,97p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
Inherited from S11c-a §1b unchanged: the rest-frame bulk fields `v_bulk=∇₄φ`, `δp=−ρ_m∂_tφ`,
`∂_t²φ=c_s0²∇₄²φ`; the current and conservation law `j=ρ_4D v_bulk`, `∂_tρ_4D+∇₄·j=0`; and the
dynamic, anchored slab window `Ω` supplied in S11c-a §3. S11c-b performs no curved-bulk response solve (§0).
````

````
$ sed -n 145,154p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
J_s = Λ_A(ω)𝒜_s + Λ_V(ω)V_s ,   Λ_I(ω)=Λ_I⁰/(1−iωτ_I) ,  I∈{A,V,X} ,
𝒜_s = μ_s − δp_s/ρ_m ,   μ_s = μ_θ/ρ_br⁰ ,
t_s = −(δp_s + Λ_X(ω)𝒜_s)n̂_s ,   n̂_s·v_bulk,s = V_s + J_s/ρ_m ,
∂_tΣ + ∇_x·(Σ v) = −(J₊+J₋) ,   Σ ≡ Σ_E ≡ ρ_4D W ,   v ≡ ∂_t u ,
δ_vΣ_mat = 0 ,   (uniform linearisation)  δ_vθ + δ_ve_W + ∇_x·δ_vu = 0 .
```

The three `τ_I` are independent. Equations of motion are obtained by S11b's method — balance laws, the
binding virtual-displacement rule, variational derivatives with held-fixed fields named, and prescribed
external virtual work — **not** by putting an irreversible response kernel in an ordinary action.
````

````
$ sed -n 180,187p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
### 2a · The background ansatz — inherited

The background ansatz is S11c-a §2 imported unchanged: the constant bindings `W̄₀≡W_0`, `μ̄_R≡mu_R`; the
fresh varying profiles on the anchor coordinate `y`,

```text
ξ ≡ y/L_W ,   W_bg(y) ≡ W̄₀[1+η w₁(ξ)] ,   μ_R,bg(y) ≡ μ̄_R[1+η m₁(ξ)] ,   σ_W ≡ η W̄₀/L_W ,
∂_{yᵢ}W_bg = σ_W ∂_{ξᵢ}w₁ ,   ∂_{yᵢ}μ_R,bg = (μ̄_R/W̄₀) σ_W ∂_{ξᵢ}m₁ ;
````

````
$ sed -n 195,196p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
𝔅⁰ ≡ {W_bg, μ_R,bg, ρ_4D,bg⁰, ρ_br,bg⁰, θ⁰, V_s⁰, J_s⁰, 𝒜_s⁰, boundary loads} ,   θ⁰ ≡ 0 ,
𝒮_hold⁰ ≡ {f_hold⁰(x), t_hold,s⁰(x)} ,   V_s⁰ = J_s⁰ = 𝒜_s⁰ = 0 .
````

````
$ grep -n '^### 3c' research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
312:### 3c · The off-diagonal coupling kernel
````

````
$ sed -n 312,346p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
### 3c · The off-diagonal coupling kernel

Extract the block of the §3b operator that couples the transverse structure to the `{θ,e_W,u_L}` structure
— the object whose uniform limit is S11b's decoupled zero (§1d) — as a **weak variational restriction** under
the supplied stored-energy/kinetic variational pairing of §1c. Use independent transverse trial and test
displacements `u_T,v_T` with `∇·u_T=∇·v_T=0`, independent longitudinal trial and test displacements `u_L,v_L`
with `∇×u_L=∇×v_L=0`, and independent trial and test fields for `θ` and `e_W`. A transverse trial paired with
a `{θ,e_W,u_L}` test defines the transverse→thickness block; a `{θ,e_W,u_L}` trial paired with a transverse
test defines the thickness→transverse block. These are restrictions to local differential-operator-labelled
trial/test spaces, not a global spectral or Helmholtz projection. They attribute both the gradient-structured
terms and the undifferentiated-`u` spurion couplings such as `g·u` (§1a/`N15`) to a sector without introducing
a global projector (`N5`).

⛔ Do **not** implement the split by setting only the **undifferentiated** field occurrences to zero: that
projection is inert on gradient content and leaves diagonal thickness/longitudinal dynamics inside the
block. ⛔ Extract both weak blocks from the §3b operator itself, not from a parallel direct-variation route.
For the in-plane domain, take all trial and test fields to have compact support in its interior, so the
in-plane integration-by-parts boundary term is fixed to zero; the inherited face boundary conditions still
apply.

Emit both off-diagonal blocks and their pairing-based adjointness residual. Define that residual by applying
the two weak operator restrictions to arbitrary independent cross-sector trial/test fields and comparing the
transverse→thickness pairing with the formal-adjoint thickness→transverse pairing. ⛔ It is **not** the mixed
second derivative of the scalar energy (`∂²U/∂u_T∂e_W − ∂²U/∂e_W∂u_T`), which is zero for any energy by the
commutation of partial derivatives and tests nothing (rule 2 corollary 3). If the two blocks are adjoint by
construction, emit the two blocks and state that there is no independent second route rather than dressing
a structural zero as a check. Emit the two blocks, the pairing-based adjointness residual when it is an
independent route, and the kernel's `(ε,η,σ_W)` multigrade and dimension, per anchoring and density
representative. ⛔ Do not filter the kernel to a single channel and do not state its coefficient, sign,
parity, or grade — the operator's own weak block extraction is the computation.

```text
⇒ S11CB_COUPLING_KERNEL , S11CB_COUPLING_KERNEL_TERM_ORIGINS .
```

````

````
$ sed -n 15,31p research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md
> ⚠⚠ **STATUS OF THE CLOSE (2026-09-03).** S11c-b closes on **per-engine leg-verification** + a coarse single-case
> cross-engine consistency check. **The full CROSS-ENGINE RESIDUAL is DEFERRED to a ≥64 GB box by USER CHOICE** — a
> lighter CORE-only residual (SLAB rows + `FINAL_KERNEL` + μ_θ + energy + admissibility, ~8 GB, which DOES fit this
> 30 GB box) was OFFERED and NOT taken (`66e8d021`; `DEFERRED_HEAVY_RUNS.md`). The blocker for the FULL residual: the
> fresh WL primaries `.out` also builds `COUPLING_KERNEL_TERM_ORIGINS` (a full `extractCoupling` per origin), which
> exceeded 30 GB (one case > 15.6 GB and growing → OOM; STEP 0's on-box claim was overturned —
> `_measurements/S11c_b_step0_residual_scope.md`). ⇒ the reviewed comparator decision lists P2a (slab-row join +
> `row_residual`) and P2b (§3a scale bridge) are committed as the SPEC for that run, but their BUILDS + the run + the
> #88 re-adjudication + the 2 owed control-hardenings are ≥64 GB work. **The committed in-tree `.out` are STALE for
> this record** (WL at `d4adbd99` = pre-#89b, frozen operator; PY = pre-fold PRIMARIES_ONLY) — a reader of them will
> not see the objects below; the fresh post-fold/#90 `.out` live only in `~/.s11_build/` scratch. `S11c_b_exports.py`
> was regenerated by the fresh PY run (faithful — `BUILD_INPUT_DIGESTS` match the committed folded+#90 engine + spec)
> and COMMITTED `af560257` (per user direction, satisfying N1's per-sub-step exports obligation); ⚠ its coupling
> content is per-engine leg-verified (fold pin B + #90) but NOT yet cross-engine-validated, so S11c-c imports
> per-engine-verified coupling with the cross-engine residual still deferred. The two
> whole-row SIGN conventions (kinetic −K PY vs +K WL; face generalized-force) and #90's two flags (closure-fold sign;
> uniform-limit Λ survivor) remain **cross-engine-UNVALIDATED**.
````

````
$ sed -n 35,37p research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md
density representatives ρ4D/ρbr, two anchorings LAB_HELD/MATERIAL_ADVECTED), S11c-b computes: (1) the **§3a energy
basis** — the O(3)-Kronecker field-bilinear invariant family, corrected to **40 = 10 uniform + 15 ∂W_bg-spurion + 15
∂μ_R,bg-spurion**; (2) the **variable-coefficient slab OPERATOR ROWS** — the equations of motion for U/θ/e_W with the
````

````
$ sed -n 112,114p research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md
  in-band-but-verified-out-of-band: the #89 PY skipped controls + the #89b WL heavy controls.) **Two whole-row SIGN
  CONVENTIONS** to adjudicate there (kinetic −K PY vs +K WL; face generalized-force PY `+diff` vs WL
  `−linearVirtualVariation`) — the comparator SURFACES them, does not normalize them (rule 1/6). **#90's two flags**
````

````
$ sed -n 343,354p research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
```text
δ_vx_s^α ≡ δ_vR_s^α|_X ,              v_face,s^α ≡ ∂_tR_s^α|_X ,
V_s^α ≡ V_{n,s}^α ≡ n̂_s^α·v_face,s^α .
```

Use the same `V_s^α` in every object below. With all bulk quantities traced at `R_s^α`, define once

```text
J_s^α ≡ ρ_m (v_bulk,s − v_face,s^α)·n̂_s^α ,
n̂_s^α·v_bulk,s = V_s^α + J_s^α/ρ_m ,
𝒜_s^α ≡ μ_s^α − δp_s^α/ρ_m ,
t_s^α ≡ −(δp_s^α + Λ_X(ω) 𝒜_s^α)n̂_s^α ,
````

````
$ sed -n 365,366p research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
δ_v𝒲_bulk^α ≡ Σ_s a_s^α t_s^α·δ_vx_s^α ,
∂_tΣ^α + ∇_x·(Σ^α v) = −Σ_s a_s^α J_s^α ,       v = ∂_t u .
````

````
$ grep -n '^DIRECTIONS =' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
57:DIRECTIONS = range(3)
````

````
$ sed -n 186,195p research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
grad_mu = tuple(
    inherited_symbol(
        f"mu_R_bg_d{i}", "DERIVED", f"background modulus first jet {i}", dim_div(DIM_ENERGY, DIM_L)
    )
    for i in range(1, 4)
)
w1_grad = tuple(
    inherited_symbol(f"w1_profile_d{i}", "KNOB", f"thickness-profile first jet {i}", DIM_ZERO)
    for i in range(1, 4)
)
````

````
$ sed -n 418,422p research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
    for j in range(i, 4):
        bind_additional_inherited(
            f"w1_profile_d{i}d{j}",
            "KNOB",
            f"dimensionless thickness-profile second jet {i},{j}",
````

````
$ grep -n 'CENTER_FACE_GENERALIZED_ROW' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
2237:        "CENTER_FACE_GENERALIZED_ROW": center_face,
3112:            Str("CENTER_FACE_GENERALIZED_ROW"),
3113:            sp.sympify(face_rows["CENTER_FACE_GENERALIZED_ROW"]),
````

````
$ grep -n 'spatialCoordinates = ' research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl
241:spatialCoordinates = {xOne, xTwo, xThree};
````

````
$ ls -la research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl research/pde_ledger_v3/scripts/S11c_b_exports.py
-rw-rw-r-- 1 trevnorris trevnorris   140124 Sep  3 02:10 research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl
-rw-rw-r-- 1 trevnorris trevnorris   203199 Sep 14 04:04 research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
-rw-rw-r-- 1 trevnorris trevnorris 61624077 Sep 14 05:49 research/pde_ledger_v3/scripts/S11c_b_exports.py
````

````
$ git log -1 --format='%h %s' af560257
af560257 S11c-b exports.py regen (folded + #90 LEDGER) — faithful, digests match committed inputs
````

````
$ ls research/pde_ledger_v3/steps/ | grep -i c1
S11c_c1_curved_bulk_closure.md
````

## v4 added lookups

````
$ sed -n '24,33p' AGENTS.md
Add explicit diff exclusions for newly generated data files exceeding 1 MiB;
keep handwritten code, physical input files and concise reports reviewable.

# Standing execution-time policy (user instruction, September 30, 2026)

All runs have no wall-clock, native alarm, CPU-time or inactivity deadline unless
the user explicitly requests a limit for a particular run. Do not add automatic
time caps to save model tokens or split a computation into timed continuations.
Keep memory containment, zero swap, host-memory protection, process/thread/CPU
controls, overlap refusal, scientific failure checks and durable checkpoints.
````

````
$ sed -n '143,169p;214,233p;569,577p' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
# Exact inherited identities required by the chain directive.
c_s0 = inherited_symbol("c_s0", "KNOB", "bulk sound speed", DIM_VELOCITY)
mu_R = inherited_symbol("mu_R", "KNOB", "uniform curl modulus", DIM_ENERGY)
rho_br = inherited_symbol("rho_br", "KNOB", "uniform integrated brane density", DIM_RHOBR)
W0 = inherited_symbol("W_0", "KNOB", "uniform reference thickness", DIM_L)
e_W = inherited_symbol("e_W", "COORDINATE", "uniform-reference thickness fraction", DIM_ZERO)
rho_m = inherited_symbol("rho_m", "KNOB", "bulk mass density", DIM_RHO4)
v_dr = inherited_symbol(
    "v_bulk_normal_0", "KNOB", "inert bulk-normal drain scope-limit speed", DIM_VELOCITY
)

# Other supplied constants needed by the independently constructed basis and balance laws.
B_rho_3 = inherited_symbol("B_rho_3", "KNOB", "uniform integrated compression modulus", DIM_ENERGY)
C = inherited_symbol("C", "KNOB", "densification-thickness coefficient", dim_div(DIM_ENERGY, DIM_L))
k_W = inherited_symbol("k_W", "KNOB", "thickness restoring coefficient", dim_div(DIM_ENERGY, 2 * DIM_L))
kappa_W = inherited_symbol(
    "kappa_W", "KNOB", "thickness-gradient coefficient", dim_div(DIM_ENERGY, 2 * DIM_L)
)
mu_W = inherited_symbol("mu_W", "KNOB", "thickness inertia", dim_mul(DIM_M, -3 * DIM_L))
mu_S = inherited_symbol("mu_S", "KNOB", "symmetric-gradient modulus", DIM_ENERGY)
G_theta_u = inherited_symbol("G_theta_u", "KNOB", "densification-divergence coefficient", DIM_ENERGY)
G_W_u = inherited_symbol("G_W_u", "KNOB", "thickness-divergence coefficient", DIM_ENERGY)
kappa_theta = inherited_symbol(
    "kappa_theta", "KNOB", "densification-gradient coefficient", dim_mul(DIM_ENERGY, 2 * DIM_L)
)
kappa_theta_W = inherited_symbol(
    "kappa_theta_W", "KNOB", "mixed scalar-gradient coefficient", dim_mul(DIM_ENERGY, 2 * DIM_L)
# Wave coordinates use the S11c-a identities where available.
theta = inherited_symbol("theta", "COORDINATE", "Eulerian densification", DIM_ZERO)
u = tuple(
    inherited_symbol(f"u_{a}", "COORDINATE", f"displacement component {a}", DIM_L)
    for a in range(1, 4)
)
grad_u = tuple(
    tuple(
        inherited_symbol(f"u_{a}_d{i}", "COORDINATE", f"displacement first jet {a},{i}", DIM_ZERO)
        for i in range(1, 4)
    )
    for a in range(1, 4)
)
grad_theta = tuple(
    inherited_symbol(f"grad_theta_{i}", "COORDINATE", f"densification first jet {i}", -DIM_L)
    for i in range(1, 4)
)
grad_e = tuple(
    inherited_symbol(f"e_W_d{i}", "COORDINATE", f"thickness-fraction first jet {i}", -DIM_L)
    for i in range(1, 4)
bind_additional_inherited("zeta_c", "COORDINATE", "face-centre displacement", DIM_L)
bind_additional_inherited("zeta_c_t", "COORDINATE", "face-centre velocity", DIM_VELOCITY)
bind_additional_inherited(
    "delta_v_zeta_c", "COORDINATE", "virtual face-centre displacement", DIM_L
)
for i in range(1, 4):
    bind_additional_inherited(
        f"zeta_c_d{i}", "COORDINATE", f"face-centre displacement jet {i}", DIM_ZERO
    )
````

````
$ sed -n '1760,1837p;4119,4141p' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
    selected_terms: list[tuple[str, sp.Expr]] = []
    basis_rows = []

    for index in UNIFORM_SELECTED:
        label, abstract_invariant = UNIFORM_CANDIDATES[index]
        invariant = sp.expand(abstract_invariant.subs(uniform_substitution, simultaneous=True))
        coefficient = uniform_coefficient(label, abstract_invariant)
        if ablation_source == "W_BG" and ablation_direction is not None:
            coefficient = coefficient.subs(grad_W[ablation_direction], 0)
        if ablation_source == "MU_R_BG" and ablation_direction is not None:
            coefficient = coefficient.subs(grad_mu[ablation_direction], 0)
        term = sp.expand(coefficient * invariant)
        selected_terms.append((label, term))
        coefficient_jet = sp.Tuple(
            *(
                total_derivative(
                    coefficient,
                    i,
                    background_depth=background_depth,
                    background_first=background_first,
                    zero_sources=zero_sources,
                )
                for i in DIRECTIONS
            )
        )
        basis_rows.append(
            sp.Tuple(
                Str(label),
                invariant,
                coefficient,
                epsilon**2 * term,
                coefficient_jet,
                dimension_of(epsilon**2 * term),
            )
        )

    new_rows = []
    new_omissions = []
    for source, g_vector, actual_vector in (
        ("W_BG", tuple(bg[i] if g_w[i] != 0 else sp.Integer(0) for i in DIRECTIONS), g_w),
        ("MU_R_BG", tuple(bg[i] if g_mu[i] != 0 else sp.Integer(0) for i in DIRECTIONS), g_mu),
    ):
        candidates = enumerate_new_candidates(g_vector)
        expressions = tuple(expression for _, expression in candidates)
        signatures = basis_euler_signatures(
            expressions,
            basis_fields,
            background_first_jets=bg,
            background_second_jets=basis_background_second,
            background_depth=background_depth,
        )
        selected, omitted = quotient_independent_indices(expressions, signatures)
        substitution = live_basis_substitution(
            actual_vector,
            source=source,
            background_depth=background_depth,
            background_first=background_first,
            zero_sources=zero_sources,
        )
        for index in selected:
            label, abstract_invariant = candidates[index]
            invariant = sp.expand(abstract_invariant.subs(substitution, simultaneous=True))
            coefficient = NEW_COEFFICIENTS[source][index]
            invariant_dimension = dimension_of(invariant)
            SYMBOL_DIMENSIONS[coefficient] = dim_div(DIM_ENERGY, invariant_dimension)
            term = sp.expand(coefficient * invariant)
            full_label = f"{source}_{label}"
            selected_terms.append((full_label, term))
            new_rows.append(
                sp.Tuple(
                    Str(source),
                    Str(label),
                    invariant,
                    coefficient,
                    epsilon**2 * term,
                    dimension_of(epsilon**2 * term),
                )
            )
def task_energy_basis() -> None:
    variable = {}
    counts = {}
    new_counts = {}
    new = {}
    omissions = {}
    for branch in BRANCHES:
        build = construct_energy(branch)
        ENERGY_PRIMARY_CASES[branch] = build
        variable[branch] = retained_energy_basis(build.basis)
        counts[branch] = build.count
        new_counts[branch] = build.new_counts
        new[branch] = retained_grade(build.new_invariants)
        omissions[branch] = retained_grade(build.omissions)
    emit_primary("ENERGY_BASIS_VARIABLE", variable, "energy_basis_variable")
    emit_primary("ENERGY_BASIS_COUNT", counts, "energy_basis_count")
    emit_primary(
        "ENERGY_BASIS_NEW_INVARIANT_COUNT",
        new_counts,
        "energy_basis_new_invariant_count",
    )
    emit_primary("ENERGY_BASIS_NEW_INVARIANTS", new, "energy_basis_new_invariants")
    emit_primary("ENERGY_BASIS_OMISSIONS", omissions, "energy_basis_omissions")
````

````
$ sed -n '239,242p' research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl

braneDimension = 3;
spatialCoordinates = {xOne, xTwo, xThree};
materialCoordinates = {capitalXOne, capitalXTwo, capitalXThree};
````

````
$ sed -n '95,100p;150,165p;578,592p;629,646p' research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py
theta = inherited("theta", "COORDINATE", "Eulerian fractional densification")
zeta_c = inherited("zeta_c", "COORDINATE", "face-centre displacement")
delta_v_theta = inherited("delta_v_theta", "COORDINATE", "virtual densification variation")
delta_v_e_W = inherited("delta_v_e_W", "COORDINATE", "virtual fractional-thickness variation")
delta_p_plus = inherited("delta_p_plus", "COORDINATE", "upper perturbation face pressure")
delta_p_minus = inherited("delta_p_minus", "COORDINATE", "lower perturbation face pressure")
# Traced bulk perturbations and their normal jets at the flat reference faces.
# These are wave-field coordinates, not premises about the physical
# background.  The background normal jets used by the supplied shifted-trace
# law are computed below from the members of the supplied background state.
delta_v_bulk = {
    s: tuple(symbol(f"delta_v_bulk_{'plus' if s == 1 else 'minus'}_{i}", "COORDINATE", f"bulk velocity perturbation at flat reference face {s}, component {i}") for i in range(1, 5))
    for s in FACES
}
dw_delta_v_bulk = {
    s: tuple(symbol(f"d_w_delta_v_bulk_{'plus' if s == 1 else 'minus'}_{i}", "COORDINATE", f"bulk velocity perturbation normal jet at flat reference face {s}, component {i}") for i in range(1, 5))
    for s in FACES
}
dw_delta_p = {
    1: symbol("d_w_delta_p_plus", "COORDINATE", "upper pressure-perturbation normal jet at the flat reference face"),
    -1: symbol("d_w_delta_p_minus", "COORDINATE", "lower pressure-perturbation normal jet at the flat reference face"),
}
def affine_bulk_perturbation(
    reference_value: sp.Expr,
    reference_normal_jet: sp.Expr,
    face: int,
) -> sp.Expr:
    """Construct the supplied wave ansatz around the flat reference face."""
    reference_height = sp.Rational(face, 2) * W0
    return reference_value + (w - reference_height) * reference_normal_jet


def compose_physical_trace(
    background_field: sp.Expr,
    perturbation_field: sp.Expr,
    h0: sp.Expr,
    dh: sp.Expr,
def physical_trace_fields(
    face: int,
    h0: sp.Expr,
    dh: sp.Expr,
    parameter: sp.Symbol,
    rho4_bg_exact: sp.Expr,
) -> tuple[sp.Expr, tuple[sp.Expr, ...], sp.Expr, tuple[sp.Expr, ...]]:
    """Build all physical traces from the common §3c field composition."""
    pressure_reference = delta_p_plus if face == 1 else delta_p_minus
    pressure_perturbation = affine_bulk_perturbation(
        pressure_reference, dw_delta_p[face], face,
    )
    velocity_perturbation = tuple(
        affine_bulk_perturbation(
            delta_v_bulk[face][i], dw_delta_v_bulk[face][i], face,
        )
        for i in range(4)
    )
````

````
$ sed -n '484,524p' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
for face_name in ("plus", "minus"):
    bind_additional_inherited(
        f"delta_p_{face_name}", "COORDINATE", f"{face_name} face pressure perturbation", DIM_PRESSURE
    )
    bind_additional_inherited(
        f"d_w_delta_p_{face_name}",
        "COORDINATE",
        f"{face_name} face pressure normal jet",
        dim_div(DIM_PRESSURE, DIM_L),
    )
    bind_additional_inherited(
        f"delta_rho_4D_face_{face_name}",
        "COORDINATE",
        f"{face_name} face density perturbation",
        DIM_RHO4,
    )
    bind_additional_inherited(
        f"d_w_delta_rho_4D_face_{face_name}",
        "COORDINATE",
        f"{face_name} face density normal jet",
        dim_div(DIM_RHO4, DIM_L),
    )
    for component in range(1, 5):
        bind_additional_inherited(
            f"delta_v_bulk_{face_name}_{component}",
            "COORDINATE",
            f"{face_name} bulk velocity perturbation {component}",
            DIM_VELOCITY,
        )
        bind_additional_inherited(
            f"d_w_delta_v_bulk_{face_name}_{component}",
            "COORDINATE",
            f"{face_name} bulk velocity normal jet {component}",
            dim_div(DIM_VELOCITY, DIM_L),
        )
        bind_additional_inherited(
            f"d_w_delta_j_bulk_{face_name}_{component}",
            "COORDINATE",
            f"{face_name} bulk-current normal jet {component}",
            dim_div(DIM_FLUX, DIM_L),
        )
````

````
$ sed -n '2135,2239p;2967,3021p;3095,3127p' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
def face_generalized_force_rows(
    bundle: sp.Tuple,
    branch: str,
    representative: str,
    evolution_origins: sp.Tuple,
    mu_theta_amplitude: sp.Expr,
) -> dict[str, object]:
    """Extract supplied face virtual-work coefficients as operator rows.

    The S11c-a virtual-work object is already a weak density.  Differentiating it
    with respect to its independent virtual displacements produces generalized
    force rows; the later §3c restriction acts on those rows exactly once.
    """

    thickness_work_source = selected_substrate_axes(
        bundle,
        "virtual_work_shape_deriv",
        (branch, "DELTA_W", "DELTA_W", representative),
    )
    center_work_source = selected_substrate_axes(
        bundle,
        "virtual_work_shape_deriv",
        (branch, "DELTA_W", "ZETA_C", representative),
    )
    thickness_work = bind_mu_theta_operand(
        thickness_work_source,
        branch,
        mu_theta_amplitude,
    )
    center_work = bind_mu_theta_operand(
        center_work_source,
        branch,
        mu_theta_amplitude,
    )
    thickness_density = sp.sympify(thickness_work[1])
    center_density = sp.sympify(center_work[1])
    u_face = tuple(
        sp.expand(sp.diff(thickness_density, delta_v_u[a]))
        for a in DIRECTIONS
    )
    e_face = sp.expand(sp.diff(thickness_density, delta_v_e_W))
    center_face = sp.expand(
        sp.diff(center_density, INCOMING_LEDGER["delta_v_zeta_c"]["value"])
    )
    face_flux = sp.sympify(
        named_tuple_row(evolution_origins, "TRUE_AREA_FACE_FLUX")
    )
    source_operands = sp.Tuple(
        sp.Tuple(
            Str("TRACTION"),
            sp.Tuple(
                *(
                    casify(
                        bind_mu_theta_operand(
                            selected_substrate_axes(
                                bundle,
                                "traction",
                                (branch, face, "DELTA_W", representative),
                            ),
                            branch,
                            mu_theta_amplitude,
                        )
                    )
                    for face in FACES
                )
            ),
        ),
        sp.Tuple(
            Str("VIRTUAL_WORK_SHAPE_DERIV"),
            casify(thickness_work),
        ),
        sp.Tuple(
            Str("CLOSURE_SHAPE_DERIV"),
            sp.Tuple(
                *(
                    casify(
                        bind_mu_theta_operand(
                            selected_substrate_axes(
                                bundle,
                                "closure_shape_deriv",
                                (branch, face, "DELTA_W", representative),
                            ),
                            branch,
                            mu_theta_amplitude,
                        )
                    )
                    for face in FACES
                )
            ),
        ),
        sp.Tuple(
            Str("MU_THETA_FACE_BINDING"),
            sp.Tuple(
                Str("LIVE_MU_THETA_AMPLITUDE"),
                mu_theta_amplitude,
            ),
        ),
    )
    return {
        "U": u_face,
        "E_W": e_face,
        "THETA_FACE_FLUX": face_flux,
        "CENTER_FACE_GENERALIZED_ROW": center_face,
        "SOURCE_OPERANDS": source_operands,
    }
    operator["U_BODY_BALANCE"] = sp.Tuple(
        sp.Tuple(
            Str("LOCAL"),
            sp.Tuple(
                *(
                    sp.expand(raw_u_local[a] + reaction_u[a]["LOCAL"])
                    for a in DIRECTIONS
                )
            ),
        ),
        sp.Tuple(
            Str("DIVERGENCE_FLUX"),
            sp.Tuple(
                *(
                    sp.Tuple(
                        *(
                            sp.expand(raw_u_flux[a][i] + first_jet_flux(reaction_u[a])[i])
                            for i in DIRECTIONS
                        )
                    )
                    for a in DIRECTIONS
                )
            ),
        ),
        sp.Tuple(
            Str("EXPANDED"),
            sp.Tuple(
                *(
                    sp.expand(reduced_u[a]["EXPANDED"] + u_kinetic[a])
                    for a in DIRECTIONS
                )
            ),
        ),
    )
    raw_e_flux = tuple(named_tuple_row(raw_e_balance, "DIVERGENCE_FLUX"))
    reaction_e_flux = first_jet_flux(reaction_e)
    operator["E_W_BALANCE"] = sp.Tuple(
        sp.Tuple(
            Str("LOCAL"),
            sp.expand(named_tuple_row(raw_e_balance, "LOCAL") + reaction_e["LOCAL"]),
        ),
        sp.Tuple(
            Str("DIVERGENCE_FLUX"),
            sp.Tuple(
                *(
                    sp.expand(raw_e_flux[i] + reaction_e_flux[i])
                    for i in DIRECTIONS
                )
            ),
        ),
        sp.Tuple(
            Str("EXPANDED"),
            sp.expand(reduced_e["EXPANDED"] + e_kinetic),
        ),
    )
    operator["THETA_BALANCE"] = sp.Tuple(
        sp.Tuple(Str("SOURCE_OPERAND"), mass_balance),
        sp.Tuple(Str("EVOLUTION_TERM_ORIGINS"), evolution_origins),
        sp.Tuple(Str("EXPANDED"), mass_balance),
    )
    operator["ADVECTIVE_MASS_OPERAND"] = named_tuple_row(
        evolution_origins, "BACKGROUND_ADVECTION"
    )
    operator["FACE_FLUX_BOUNDARY_OPERANDS"] = faces
    operator["FACE_GENERALIZED_FORCE_ROWS"] = sp.Tuple(
        sp.Tuple(Str("U"), sp.Tuple(*physical_face_u)),
        sp.Tuple(Str("E_W"), physical_face_e),
        sp.Tuple(
            Str("THETA_FACE_FLUX"),
            sp.sympify(face_rows["THETA_FACE_FLUX"]),
        ),
        sp.Tuple(
            Str("CENTER_FACE_GENERALIZED_ROW"),
            sp.sympify(face_rows["CENTER_FACE_GENERALIZED_ROW"]),
        ),
        sp.Tuple(Str("SOURCE_OPERANDS"), face_rows["SOURCE_OPERANDS"]),
    )

    if branch == "MATERIAL_ADVECTED":
        mu_theta_value = sp.Tuple(Str("mu_theta^alpha"), epsilon * mu_theta_amplitude)
    else:
        mu_theta_value = sp.Tuple(Str("mu_theta"), epsilon * mu_theta_amplitude)
    reserved_mu_theta = INCOMING_LEDGER[
        "mu_theta_M" if branch == "MATERIAL_ADVECTED" else "mu_theta_L"
    ]["value"]
    operator["MU_THETA_FACE_BINDING"] = sp.Tuple(
        reserved_mu_theta, mu_theta_amplitude
    )
````

````
$ sed -n '600,626p;878,944p' research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py
    face_height = h0 + parameter * dh
    return bulk_field.subs(w, face_height)


@dataclass(frozen=True)
class FaceSource:
    branch: str
    face: int
    dof: str
    representative: str
    route: str
    parameter: sp.Symbol
    h0: sp.Expr
    dh: sp.Expr
    normal_exact: tuple[sp.Expr, ...]
    measure_exact: sp.Expr
    virtual_displacement: tuple[sp.Expr, ...]
    face_velocity_exact: tuple[sp.Expr, ...]
    pressure_trace_exact: sp.Expr
    bulk_velocity_trace_exact: tuple[sp.Expr, ...]
    density_trace_exact: sp.Expr
    current_trace_exact: tuple[sp.Expr, ...]
    rho4_bg_exact: sp.Expr
    rhobr_bg_exact: sp.Expr


_FACE_CACHE: dict[tuple[object, ...], FaceSource] = {}
def face_normal_raw(source: FaceSource) -> sp.Tuple:
    background = tuple(component.subs(source.parameter, 0) for component in source.normal_exact)
    derivative = shape(source.normal_exact, source.parameter)
    return sp.Tuple(sp.ImmutableMatrix(background), epsilon * sp.ImmutableMatrix(derivative), sp.ImmutableMatrix(background) + epsilon * sp.ImmutableMatrix(derivative))


def face_measure_raw(source: FaceSource) -> sp.Tuple:
    background = source.measure_exact.subs(source.parameter, 0)
    derivative = shape(source.measure_exact, source.parameter)
    return sp.Tuple(background, epsilon * derivative, background + epsilon * derivative)


def face_velocity_raw(source: FaceSource) -> sp.Expr:
    return epsilon * shape(dot(source.normal_exact, source.face_velocity_exact), source.parameter)


def relative_flux_raw(source: FaceSource) -> sp.Expr:
    relative = tuple(source.bulk_velocity_trace_exact[i] - source.face_velocity_exact[i] for i in range(4))
    exact = rho_m * dot(relative, source.normal_exact)
    return epsilon * shape(exact, source.parameter)


def true_area_flux_raw(source: FaceSource) -> sp.Expr:
    relative = tuple(source.bulk_velocity_trace_exact[i] - source.face_velocity_exact[i] for i in range(4))
    exact_flux = rho_m * dot(relative, source.normal_exact)
    return epsilon * shape(source.measure_exact * exact_flux, source.parameter)


def pressure_trace_raw(source: FaceSource) -> sp.Expr:
    return epsilon * shape(source.pressure_trace_exact, source.parameter)


def mu_specific_raw(source: FaceSource) -> sp.Expr:
    exact = source.parameter * mu_theta_branch[source.branch] / source.rhobr_bg_exact
    return epsilon * shape(exact, source.parameter)


def affinity_raw(source: FaceSource) -> sp.Expr:
    return mu_specific_raw(source) - pressure_trace_raw(source) / rho_m


def traction_raw(source: FaceSource) -> sp.ImmutableMatrix:
    exact_pressure = source.pressure_trace_exact
    exact_mu = source.parameter * mu_theta_branch[source.branch] / source.rhobr_bg_exact
    exact_affinity = exact_mu - exact_pressure / rho_m
    exact = -(exact_pressure + Lambda_X * exact_affinity) * sp.ImmutableMatrix(source.normal_exact)
    return epsilon * sp.ImmutableMatrix(shape(exact, source.parameter))


def closure_raw(source: FaceSource) -> sp.Expr:
    return sp.factor_terms(relative_flux_raw(source) - Lambda_A * affinity_raw(source) - Lambda_V * face_velocity_raw(source))


def kinematic_raw(source: FaceSource) -> sp.Tuple:
    operand_a = epsilon * shape(
        dot(source.normal_exact, source.bulk_velocity_trace_exact),
        source.parameter,
    )
    flux_from_relative_law = relative_flux_raw(source)
    operand_b = sp.factor_terms(face_velocity_raw(source) + flux_from_relative_law / rho_m)
    residual = sp.factor_terms(operand_a - operand_b)
    return sp.Tuple(
        sp.Tuple(Str("OPERAND_A"), operand_a),
        sp.Tuple(Str("OPERAND_B"), operand_b),
        sp.Tuple(Str("RESIDUAL"), residual),
    )

````

````
$ sed -n '133,148p;566,575p;780,804p;840,861p' research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py
# Wave and virtual jets.  They remain distinct operands; in particular, no
# physical u is substituted for delta_v_u in T-g.
u = tuple(symbol(f"u_{i}", "COORDINATE", f"in-plane displacement component {i}") for i in range(1, 4))
u_t = tuple(symbol(f"u_{i}_t", "COORDINATE", f"time derivative of displacement component {i}") for i in range(1, 4))
grad_u = tuple(tuple(symbol(f"u_{i}_d{j}", "COORDINATE", f"displacement jet {i},{j}") for j in range(1, 4)) for i in range(1, 4))
grad_u_t = tuple(tuple(symbol(f"u_{i}_t_d{j}", "COORDINATE", f"velocity jet {i},{j}") for j in range(1, 4)) for i in range(1, 4))
theta_t = symbol("theta_t", "COORDINATE", "time derivative of densification")
grad_theta = tuple(symbol(f"theta_d{i}", "COORDINATE", f"densification first jet {i}") for i in range(1, 4))
e_W_t = symbol("e_W_t", "COORDINATE", "time derivative of uniform-reference fractional thickness")
grad_e_W = tuple(symbol(f"e_W_d{i}", "COORDINATE", f"fractional-thickness first jet {i}") for i in range(1, 4))
zeta_c_t = symbol("zeta_c_t", "COORDINATE", "time derivative of centre displacement")
grad_zeta_c = tuple(symbol(f"zeta_c_d{i}", "COORDINATE", f"centre-displacement first jet {i}") for i in range(1, 4))

delta_v_u = tuple(symbol(f"delta_v_u_{i}", "COORDINATE", f"virtual displacement component {i}") for i in range(1, 4))
grad_delta_v_u = tuple(tuple(symbol(f"delta_v_u_{i}_d{j}", "COORDINATE", f"virtual-displacement jet {i},{j}") for j in range(1, 4)) for i in range(1, 4))
delta_v_zeta_c = symbol("delta_v_zeta_c", "COORDINATE", "virtual centre-face displacement")
def dof_fields(dof: str, face: int) -> tuple[sp.Expr, sp.Expr, tuple[sp.Expr, ...], sp.Expr, sp.Expr]:
    if dof == "DELTA_W":
        return (
            face * W0 * e_W / 2,
            face * W0 * e_W_t / 2,
            tuple(face * W0 * grad_e_W[i] / 2 for i in range(3)),
            W0 * e_W,
            W0 * delta_v_e_W,
        )
    return zeta_c, zeta_c_t, grad_zeta_c, sp.Integer(0), sp.Integer(0)
    virtual_parameter = sp.Dummy("material_virtual_parameter", real=True)
    virtual_centre = delta_v_zeta_c if dof == "ZETA_C" else sp.Integer(0)
    virtual_thickness = virtual_delta_W if dof == "DELTA_W" else sp.Integer(0)
    if branch == "LAB_HELD":
        virtual_anchor = W_bg + virtual_parameter * dot(delta_v_u, tuple(profile_gradient))
    else:
        virtual_anchor = W_bg
    virtual_face_height = (
        virtual_parameter * virtual_centre
        + face * (virtual_anchor + virtual_parameter * virtual_thickness) / 2
    )
    virtual_vertical = shape(virtual_face_height, virtual_parameter)
    virtual_displacement = tuple(delta_v_u) + (virtual_vertical,)

    if branch == "LAB_HELD":
        anchor_time_derivative = dot(u_t, tuple(profile_gradient))
    else:
        anchor_time_derivative = sp.Integer(0)
    centre_time_derivative = zeta_t if dof == "ZETA_C" else sp.Integer(0)
    thickness_time_derivative = W0 * e_W_t if dof == "DELTA_W" else sp.Integer(0)
    face_vertical_velocity = (
        centre_time_derivative
        + face * (anchor_time_derivative + thickness_time_derivative) / 2
    )
    face_velocity_exact = tuple(parameter * item for item in tuple(u_t) + (face_vertical_velocity,))
    parameter = sp.Dummy("shape_parameter", real=True)
    scales = source_jet_scales(ablate_direction, reverse_upper_x1, face)
    zeta, zeta_t, _, _, virtual_delta_W = dof_fields(dof, face)
    advective_height = sp.Add(*(u[i] * scales[i] * grad_W[i] for i in range(3)))
    h0 = face * W_bg / 2
    dh = zeta if branch == "LAB_HELD" else zeta - face * advective_height / 2
    h_exact = h0 + parameter * dh

    grad_h = tuple(dx(h_exact, i, scales) for i in range(3))
    denominator = sp.sqrt(1 + dot(grad_h, grad_h))
    normal_exact = tuple([-face * component / denominator for component in grad_h] + [face / denominator])
    measure_exact = denominator

    virtual_zeta = delta_v_zeta_c if dof == "ZETA_C" else face * virtual_delta_W / 2
    if branch == "LAB_HELD":
        virtual_vertical = virtual_zeta + face * dot(delta_v_u, tuple(scales[i] * grad_W[i] for i in range(3))) / 2
        velocity_vertical = zeta_t + face * dot(u_t, tuple(scales[i] * grad_W[i] for i in range(3))) / 2
    else:
        virtual_vertical = virtual_zeta
        velocity_vertical = zeta_t
    virtual_displacement = tuple(delta_v_u) + (virtual_vertical,)
    face_velocity_exact = tuple(parameter * item for item in tuple(u_t) + (velocity_vertical,))
````

````
$ sed -n '182,199p' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
grad_W = tuple(
    inherited_symbol(f"W_bg_d{i}", "DERIVED", f"background thickness first jet {i}", DIM_ZERO)
    for i in range(1, 4)
)
grad_mu = tuple(
    inherited_symbol(
        f"mu_R_bg_d{i}", "DERIVED", f"background modulus first jet {i}", dim_div(DIM_ENERGY, DIM_L)
    )
    for i in range(1, 4)
)
w1_grad = tuple(
    inherited_symbol(f"w1_profile_d{i}", "KNOB", f"thickness-profile first jet {i}", DIM_ZERO)
    for i in range(1, 4)
)
m1_grad = tuple(
    inherited_symbol(f"m1_profile_d{i}", "KNOB", f"modulus-profile first jet {i}", DIM_ZERO)
    for i in range(1, 4)
)
````
