# S11c_d_zinvariant_operator_blocks_directive.md (v3) — grounding commands (rule 2 / E1)

Mechanical lookups, run from the repo root on 2026-09-30, HEAD `e5f0dfda`. Regenerated from the commands below; nothing transcribed.

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

