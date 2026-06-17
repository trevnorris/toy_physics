# Directive — M1b: DERIVE the BdG spectrum + wall coefficients; POSIT mixed ports transparently

You are the **CODER** (Codex codes / Claude reviews). Build the genuine-derivation core of the Mathematica
branch producer. M0 (env check) and M1a (geometry + V2-22B handshake) are done. A Claude+Codex consult fixed
the **derivable scope** (recorded in `software/stage1_solver/decisions/03_m1b_derivable_scope.md`): the BdG
matter sector and the wall coefficients ARE canonically derivable; the **mixed-Maxwell ports are NOT yet**
(no canonical 1-D eigenproblem; `g_U/g_W` need the open `S_η^(A)` = Path A). So M1b derives what is derivable
and carries the mixed ports as **explicitly-labeled posited data** — it does NOT invent the mixed kernel.

**This is still pre-freeze and target-blind** (placeholder constitutive params; GATE A and the frozen run
come later in M1c). The deliverable is the genuine *derivation machinery* (a real eigensolve + real weak-form
integrals), not a physics result. **Review axis #1 = DERIVED-not-posited** (the BdG modes must come from a
genuine eigensolve, demonstrable), and a **transliteration-fidelity audit** will check the BdG operator
term-by-term against the canonical source — so transliterate faithfully.

## Scope IN — derive these

### A. Target-blind stationary background
Define a well-defined, target-blind open-throat stationary background to linearize around (reuse the M1a
geometry `R0(w)=a+(R_exit−a)(3x²−2x³)`). A full 2-D coupled WP1 Newton solve is **NOT required** for M1b — a
tractable, clearly-labeled stationary matter profile `ρ0`/`ψ0` consistent with the geometry/confinement is
acceptable (state exactly what you use and that it is engineering-smoke, target-blind). The fidelity that
matters at M1b is the *operator*, not the background's full self-consistency (that tightens at M1c).

### B. DERIVE the BdG (Bogoliubov) eigenproblem  — the genuine eigensolve
Transliterate the canonical linearized matter operator around the background and solve it for positive-energy
stable normal modes `{varpi_α > 0, φ_α(·)}` with a genuine eigensolver (`Eigensystem` / generalized
`Eigenvalues[{A,B}]` / sparse `Arnoldi` — your choice). Terms (transliterate faithfully):
- Kinetic `−(ℏ²/2m)[∇² − ℓ(ℓ+1)/r²]` with the ℓ=2 centrifugal `ℓ(ℓ+1)=6`.
- Confinement `V_conf` (the `(r/R0)⁴` radial wall + axial trap).
- Quintic EOS enthalpy `h=(5K/4)ρ⁴` and its linearization `δh` (P=Kρ⁵ / U=(K/4)ρ⁵ — do NOT simplify).
- Chemical potential `μ`, and the covariant gauge coupling (minimal `∂→∂−i(q/ℏ)A`).
**Canonical source-of-truth (transliterate against, quote in your report):** compact `L_BdG`
`notes/moving_throat_pde_program_compact.md` §4.7 ~L1406–1428; EOS L578–582. **Term-by-term cross-reference
(must match):** torch `coupled_branch.py:coupled_pde_residual` L231–313 (matter L263–270),
`p2_tangent.py` L186–281 (ℓ=2 centrifugal `matter_angular` L240–245), `physics.py:quintic_enthalpy` L46–47.
Report the spectrum + per-mode eigen-residual (proof the modes are genuine), and demonstrate stability
(`varpi_α>0`) and the Stieltjes moment positivity `B0·B4 − B2² ≥ 0` (V2-20 L233–243).

### C. DERIVE the wall coefficients K, M
Compute the weak-form wall matrices on a chosen frozen basis `{χ_i(w)}` (your choice — finite elements,
splines, or computed eigenfunctions):
`M_ij = ∫₀ᴸ μ_η χ_i χ_j dw`, `K_ij = ∫₀ᴸ [T_w χ'_i χ'_j + V_l χ_i χ_j] dw + Y_L(0) χ_i(L) χ_j(L)`,
with `V_l = K_η + l(l+1) T_Ω` (l=2). **Source:** V2-20 `stage_v2_20_weak_form_branch_extraction_derivation.md`
L118–136. Constitutive `μ_η,T_w,T_Ω,K_η` = placeholder target-blind (frozen at GATE A later). Report K, M and
confirm M>0.

### D. DERIVE the BdG–wall coupling c_α
Compute the physical wall→matter coupling overlap and `c_α` for each mode. The wall drive is
`δV_conf = −(V'_wall(Σ0/ℓc)/ℓc)·η` (compact L1080–1085, L1424–1428; torch
`p2_tangent.py:wall_to_matter_coefficient` L165–183 = `4·V_radial·r⁴/R0⁵`). Emit `lambda_B` ONLY as a
normalization-convention factor in the `c_α = λ_B·I_ηφ` split (state in the report that `λ_B` is convention,
`c_α` is the physical object). Report the overlaps + `c_α` + the derived `B0/B2/B4`.

### E. POSIT the mixed-Maxwell ports — transparently labeled, NOT derived, NOT tuned
Emit `Ω_U, Ω_W, R, g_U, g_W` (+ `u/w` profiles) as **explicitly-labeled posited placeholders** (sane,
arbitrary, NOT reverse-engineered from any target). Add packet metadata
`mixed_ports_status: "posited_not_derived"` and a loud note in the report that this sector is the open
frontier (no canonical 1-D mixed eigenproblem; `g_U/g_W` Path-A-blocked per decision 03). Keep `Δ_eff>1e-12`.

### F. Emit the V2-22B packet, run the V2 chain, label the residual honestly
Emit a schema-valid, target-blind packet (no forbidden target fields; `parent_action_status=
effective_wall_closure`; `boundary_class=open_impedance`; `R_exit>0`). Run the real
`research/pde_audit/scripts/stage_v2_22c_end_to_end_smoke_pipeline.py` on it (`timeout 600`). **LABEL the
residual loudly:** "BdG + wall DERIVED; mixed ports POSITED → `R_norm` is a target-blind PARTIALLY-POSITED
diagnostic, NOT a clean falsification test" (per decision 03). Do not present R_norm as a result; do not tune
anything to any target.

## Scope OUT
Deriving the mixed-Maxwell sector (a separate frontier consult is running); GATE A / the frozen run / M1c;
the full 2-D coupled WP1 (optional, not required); inventing `S_η`/any kernel; committing.

## Working rules
≤2 `math`/`wolframscript` kernels (1 here); `timeout 600` around every script RUN; iterate to exit 0; never
wrap your own session in a shell `timeout`. Files under `software/stage1_solver/mathematica/`; packet + run
artifacts → gitignored `runs/`. Keep substantive outputs deterministic (timestamp-only nondeterminism OK at
this pre-freeze stage). No `git add`/commit. No network/GPU. Do NOT read/import the firewalled
`research/pde_audit/simulation/`.

## Report back to the reviewer
The background you used (+ that it's target-blind smoke); the DERIVED BdG spectrum (varpi_α + eigen-residuals
proving genuineness) + stability + Stieltjes check; K, M (+ M>0); the overlaps + c_α + derived B0/B2/B4; the
POSITED mixed-port values (clearly labeled); the packet `validation_pass`; the V2 residual with the loud
partial-diagnostic label; and a term-by-term map of your BdG operator to the canonical source for the fidelity
audit. Then Claude reviews (transliteration-fidelity + adversarial derived-not-posited/honesty) + arbiter.
