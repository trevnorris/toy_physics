# Directive — Spike-1: ℓ=2 vector-harmonic localized-Maxwell operator + gauge deflation + MMS

You are the **CODER** (Codex codes / Claude reviews). Build the FOUNDATION of the mixed-Maxwell port spike:
the genuine ℓ=2 vector-spherical-harmonic (VSH) localized-Maxwell operator — the real replacement for the
torch *scalarized engineering-smoke* wrapper. This is Spike-1 of 3 (decision 05). **Operator + gauge handling
+ MMS + self-adjointness ONLY** — NO transfer, NO N0/Z0, NO V2 packet yet (those are Spike-2/Spike-3). This
is the highest fidelity-risk chunk; a transliteration-fidelity audit will check every lane against the parent
weak form, so transliterate faithfully.

The DESIGN is fixed and accepted — read it and follow it:
`software/stage1_solver/decisions/05_mixed_maxwell_spike_design.md` (D1, D2, D6) and the full derivation in
`software/stage1_solver/_scratch/mixed_maxwell_spike_design.log` (sections D1, D2, D6). Do not redesign the
physics; implement D1/D2 and DESIGN the numerical route.

## What to build (per D1/D2)

### A. The ℓ=2 VSH gauge-field decomposition (L = l(l+1) = 6)
`δA_0 = φ(r,w) Y`, `δA_w = a_w(r,w) Y`, brane 3-vector `= a_r Y e_r + a_E Ψ_lm + a_B Φ_lm`
(`Ψ_lm = ∇_Ω Y` poloidal/electric, `Φ_lm = e_r×∇_Ω Y` toroidal/magnetic). Five lanes `(φ, a_r, a_E, a_B, a_w)`
on the `(r,w)` grid. Reuse the M1a/M1b target-blind geometry `R0(w)` and a target-blind localization weight
`Z(w)` consistent with the torch localization (floor + Gaussian; state your placeholder).

### B. Assemble from the PARENT quadratic weak form (NOT scalarized components)
`S_A^(2) = (1/2μ0) ∫ r² dr dw Z(w) [ |E|² − |B|² − |C|² − ξ⁻¹|G|² ] − ∫ A·δJ` with the angular weights
`|E|² = |E_r|² + L|E_E|² + L|E_B|² + |E_w|²` (similarly |C|, and `|B|² = L|B_B|² + |B_r|² + L|B_E|²`), and the
gauge-invariant components exactly as in D1 (`G, E_r/E_E/E_B/E_w, C_r/C_E/C_B, B_B/B_r/B_E`). The resulting
operator MUST contain the vector-Laplacian cross-terms the smoke dropped: `−2a_r/r²`, `+2L a_E/r²`,
`+2a_r/r²` (D1 vector-Laplacian cross-check: `(Δ₃A)_r = D_L a_r − 2a_r/r² + 2L a_E/r²`,
`(Δ₃A)_E = D_L a_E + 2a_r/r²`, `(Δ₃A)_B = D_L a_B`, `D_L f = f_rr + 2f_r/r − Lf/r²`). For THIS chunk use the
conservative operator with **closed (Dirichlet/manufactured) test BCs** (the open-impedance DtN exit is
Spike-2).

### C. Gauge handling
Implement the Lorenz/Gauss residual `G`; deflate the pure-gauge null space (pure gradients
`(φ,a_r,a_E,a_w) = (iωχ, ∂_rχ, χ/r, ∂_wχ)` have `F=0`); apply one scalar gauge anchor. Report the physical
quotient (`G=0`) dimension.

## Acceptance criteria (these are the real tests — D6)
1. **MMS on ALL FIVE lanes, with discriminating power for the cross-terms.** Construct a manufactured solution
   with NONZERO coupled `a_r` AND `a_E` so the vector-Laplacian cross-terms (`−2a_r/r²`, `2L a_E/r²`,
   `2a_r/r²`) are exercised. The genuine operator must converge at the expected order; **demonstrate the
   teeth** — show that DROPPING the cross-terms (the scalarized-smoke form) makes this same MMS FAIL (so the
   test actually distinguishes the genuine reduction from the smoke). Report both.
2. **Pure-gauge deflation:** a pure-gradient input yields field strength `F ≈ 0` (to discretization order),
   zero physical norm/flux after quotienting, and is deflated from the physical operator. A non-deflated
   pure-gauge mode is a FAILURE.
3. **Self-adjointness:** the conservative operator is symmetric/self-adjoint under the `Z(w) r²` inner product
   with closed test BCs (report the asymmetry norm → machine-ε).
4. **Term map (for the fidelity audit):** map each operator term to its parent-action source (the VSH-reduced
   `S_A^(2)`, the field-strength components, the vector-Laplacian cross-terms) with file:line into decision 05 /
   the design log + compact §2.5.

## Scope OUT
The outgoing DtN/open-impedance BC, the forward `δJ_ψ` transfer, `N0/Z0` extraction, the V2-22B packet/
extension, any residual — all Spike-2/Spike-3. Do NOT invent `S_η`. Do NOT touch the firewalled
`research/pde_audit/simulation/`. No GATE-A freeze.

## Working rules
≤2 `math`/`wolframscript` kernels (1 here); `timeout 600` around every script RUN; iterate to exit 0; never
wrap your own session in a shell `timeout`. Files under `software/stage1_solver/mathematica/` (e.g.
`mt15_03_spike1_vsh_maxwell_operator.wls` + a report); run artifacts → gitignored `runs/`. Keep substantive
outputs deterministic (timestamp-only nondeterminism OK pre-freeze). No `git add`/commit. No network/GPU.
Target-blind throughout.

## Report back to the reviewer
The 5-lane operator structure + which terms carry the cross-terms; the gauge-deflation result (physical
quotient dim, pure-gauge F≈0); the self-adjointness asymmetry norm; the MMS convergence orders for all 5
lanes AND the teeth demonstration (smoke form fails the same MMS); and the parent-source term map. Then Claude
runs the transliteration-fidelity audit (VSH operator vs parent weak form) + adversarial (genuine MMS teeth?
gauge deflation real? self-adjoint? any can't-fail check?) + arbiter.
