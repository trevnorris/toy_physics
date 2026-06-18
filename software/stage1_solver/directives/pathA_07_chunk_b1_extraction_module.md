# Directive pathA_07 — Chunk B1: §5 field→D0 extraction module + independent validation

**Owner:** Codex (codes, designs the route). Claude reviews afterward (transliteration-fidelity per math layer +
adversarial). **Status:** post-freeze (GATE-A `frozen: YES`, commit `1703f4c`, `candidate_freeze_hash ed3585…`).
**Scope of THIS chunk:** build + INDEPENDENTLY VALIDATE the extraction *machinery* that maps a converged closed
background `{ψ0,A0,R0}` (and the l=2 wall eigenmode `χ`) to `{K,M,B_n,Z_n,N_n,D_n,P_n,R_norm,R_pole}` per
decision-11 §5. **Do NOT run the frozen physical closed solve or the `R_norm(τ)=0` root-find here — that is chunk
B2.** B1 proves the calculator is correct against independent oracles; B2 runs the real frozen background through it.

## Why split this out
The recurring project sin is a tautological / can't-fail gate and a faithful-but-wrong operator that MMS can't
catch. The extraction algebra must be proven correct against an INDEPENDENT oracle *before* it is ever trusted to
produce a physical `R_norm`. So B1 validates the machinery on manufactured / fixture inputs with known answers; it
must not depend on the physical answer being right.

## Build (per decision-11 §5)
A new module (e.g. `src/stage1_solver/patha_extraction.py`) with three layers. Reuse the existing algebra
verbatim where it already exists; the v2 fixtures are READ-ONLY references (do not modify them, do not import from
`research/pde_audit/simulation/`, do not touch `physical_export_permitted`).

1. **l=2 wall-stiffness eigenproblem (NEW).** On a converged `R0(w)` + a resolved `SSigmaProvider`, assemble the
   symmetric stiffness `A₂` for
   `L₂ η = −∂_w(T_{w,Σ} ∂_w η) + [K_η(w) + 6 T_{Ω,Σ}] η`, with
   `K_η = U_{Σ,RR} − ∂_w(T_{w,Σ,R} R0') + ½ T_{w,Σ,RR} (R0')²`, modal BCs `η(0)=0` (Dirichlet) and
   `T_{w,Σ}η'(L)=0` (natural). Assemble the norm matrix `W` and mass matrix `M_μ` (weight `μ_Σ`). Solve the
   generalized eigenproblem `A₂ χ = K W χ` for the **lowest positive** eigenpair; normalize `χᵀ W χ = 1` and fix
   orientation `∫χ dw > 0`. Return `K = χᵀ A₂ χ` and `M = χᵀ M_μ χ`. (For the frozen `homogeneous_isotropic_hooke_v1`
   family this reduces to `L₂ = −τ η'' + 7τ/a² η`, so `K = τ·κ̂` with `κ̂` a fixed geometric eigenvalue — useful as
   the analytic oracle below; but keep the assembly GENERAL in the provider's constitutive derivatives so it is not
   hard-wired to the harmonic case.)
2. **Field→coupling overlaps (reuse v2_22a algebra).** Using the SAME `χ`, form the overlap integrals
   (`I_{η,φ_j}`, `I_{η,u}`, `I_{η,w}`, `I_{u,w}`) against the converged `{ψ0,A0}` mode profiles and convert to BdG +
   mixed-port reduced couplings (`c_j, g_U, g_W, R_mix, …`), exactly as
   `research/pde_audit/scripts/stage_v2_22a_profile_to_coefficient_adapter.py` does. Then the moments
   `B_n = Σ_j c_j²/ϖ_j^{2(n+1)}`; per-port `Δ=Ω_U²Ω_W²−R_mix²`, `Q=…`, `P=…`; `Z_n`, `N_n` per the v2_21/v2_22a
   expansions.
3. **Observable algebra (reuse v2_21 algebra).** `D0=K−B0−Z0`, `D2=−(M+B2+Z2)`, `D4=−(B4+Z4)`; `P0=N0/D0`,
   `P2=(D0 N2−2 D2 N0)/D0²`, `P4=(D0² N4−2 D0(D2 N2+D4 N0)+3 D2² N0)/D0³`;
   `R_norm=m̂0² S_port P0 − 54 G c_s⁵/(5 a⁵ c⁵)`; `R_pole=D0(B4+Z4)−3(M+B2+Z2)²`. These MUST be byte-identical in
   form to `research/pde_audit/scripts/stage_v2_21_branch_extraction_fixture.py` (already audited faithful).

## Independent validation (the heart of B1 — every gate must be falsifiable)
1. **Eigenproblem vs analytic oracle.** For the harmonic family on `[0,L]` with `a=1`, `L=37/20`, the BVP
   `−τη''+7(τ/a²)η = K η` (Dirichlet at 0, natural/zero-slope at L) has an ANALYTIC lowest eigenvalue/eigenfunction
   (derive `κ̂` in closed form). Check: discrete `K/τ → κ̂` and `χ → χ_analytic` under ≥3-level grid refinement at
   **second order** (report observed order). This is a genuine, non-circular oracle — the discrete solver is not
   allowed to define its own target. Confirm `K` is independent of `R0(w)` for this family (the §0/§6 reportable
   finding) by re-running with a perturbed `R0` and showing `K` unchanged to tolerance.
2. **Coupling + observable layers vs the shipped fixtures.** Reproduce the `stage_v2_22a` profile→coefficient
   outputs and the `stage_v2_21` `{K,B,Z,N}→{D,P,R_norm,R_pole}` packet to their published tolerances on the
   fixtures' OWN manufactured inputs (analytic isotropic branch + at least one other fixture branch). Byte/`1e-9`
   agreement. (This exercises the `R_norm`/`R_pole` formulas on manufactured inputs — NOT on the frozen physical
   background — so it stays target-blind for the physical question.)
3. **Manufactured field→coefficient round-trip (MMS for the NEW overlap layer).** Build a manufactured
   `{χ, ψ0, A0}` with analytically known overlaps → known `{B0,Z0,N0}`; confirm the field→coupling layer recovers
   them at second order. This guards the genuinely new wiring (v2_22a validates the algebra; this validates the
   field-sampling/quadrature into it).
4. **No can't-fail gates.** Every validation must be able to FAIL: no verbatim-copy "independent" recompute (the
   1c lesson), no comparing a number to itself, no tolerance loose enough to pass a wrong operator. State for each
   gate what wrong answer it would catch.

## Acceptance criteria (Codex iterates until ALL pass, exit 0; any script wrapped `timeout 600`)
1. The module resolves the frozen `homogeneous_isotropic_hooke_v1` provider and assembles `A₂,W,M_μ`; the lowest
   positive eigenpair is found, normalized (`χᵀWχ=1`, `∫χ>0`).
2. Eigenproblem grid-convergence to the analytic `κ̂` at order ≈2 (≥3 levels; report the order + the residual to
   `κ̂`); `K`-independence-of-`R0` demonstrated for the harmonic family.
3. v2_22a and v2_21 fixture reproductions agree to ≤1e-9 (report the max abs/rel diff per packet field).
4. Manufactured field→coefficient round-trip recovers `{B0,Z0,N0}` at order ≈2.
5. A focused test module under `tests/` covers the above; full `pytest` for the patha_* suite stays green (no
   regression). Firewall untouched; no `git add`/commit (orchestrator commits).
6. NO physical frozen closed solve and NO `R_norm(τ)=0` root-find in this chunk (those are B2). The module MAY
   expose the entry points B2 will call, but B1 does not execute them on the physical background.

## Report back
Per-layer fidelity result (with the analytic `κ̂` value + observed orders + fixture max-diffs), the files
created/modified, the test results, and an explicit statement of what each validation gate would catch if the
operator were wrong. Flag any place the v2_21/v2_22a algebra had to be adapted (and why) rather than reused verbatim.
