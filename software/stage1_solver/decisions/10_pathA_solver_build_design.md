# Decision 10 — self-consistent CLOSED Path-A solver: build design

**Date:** 2026-06-17
**Mechanism:** Claude+Codex design consult (read-only). Prompt + log:
`software/stage1_solver/_scratch/pathA_solver_design_consult.{md,log}`. Codex codes the build; Claude reviews +
fidelity-audits. User gave the conceptual GO (promote `S_Σ`); this is the target-blind Stage-1 build (no
calibration, no `R_norm`, no physical export; firewall + export guard untouched).

## Design (agreed — sound)

**Engine split:** the existing torch WP1 solver is the solve engine (it owns the conservative `(r,w)` residual,
matter current, Maxwell residual, Newton/JFNK, the target-blind harnesses — `coupled_branch.py`, `newton.py`).
**Mathematica = dual-engine CROSS-CHECK** of the new continuum/static-balance algebra + the MMS forcing (NOT a
second nonlinear solver).

**New unknown + operator:** add `R0(w)` as a 1D wall-grid unknown on `grid.w_centers`; state becomes
`[ψ_R, ψ_I, A0, Ar, Aw, R0, μ]`. Keep the frozen-geometry helper for the effective-closure paths (additive — do
NOT globally replace `coupled_branch.py:124`); add a closed-path confinement call that receives the solved `R0`.
Discrete static balance, conservative wall form (faithful to derivation D3):
`−(F[j+1]−F[j])/dw + ½ T_{w,Σ,R}(Rj,wj)(R0'_j)² + U_{Σ,R}(Rj,wj) − S_j = 0`, `F[j+½] = T_{w,Σ}(R_face,w_face)
(R[j+1]−R[j])/dw`. Reuse the wall-FV pattern (`operators.py:485`) but as a NEW nonlinear `R0` operator (not the
linear `η` operator). BC: mouth FV Dirichlet `R0(0)=a`; exit default natural open Neumann / zero traction (see
pre-1b item 3).

**Return source:** the matter return enters the wall residual RHS only, using the SAME confinement derivative as
the forward wall→matter coupling (reciprocity literal; the coefficient already appears in step-8a as
`4 V_radial r⁴/R0⁵`, `p2_tangent.py:165`). The gauge return is NOT a separate local kernel — the monolithic
Newton captures the matter-mediated route through the coupled matter/Maxwell/current equations (no double-count;
derivation D2/D4, `compact:627`).

**Parameterized `S_Σ`:** a **serializable constitutive spec + registry** (hashable for freeze — NOT raw callables
in config) providing `mu, T_w, T_w_R, T_w_RR, T_Omega, U, U_R, U_RR` of `(R,w)`. Chunk-1 defaults = smooth
positive placeholder families + MMS manufactured families only. No GATE-A forms, no calibration, no `R_norm`, no
export.

**Chunk decomposition (each: clean review + transliteration-fidelity audit):**
- **1a** — static-balance operator + MMS + Mathematica cross-check. Manufacture nonconstant `R0(w)`,
  `T_w(R,w)`, nonzero `T_w_R`, `U_R`, nonzero gradient-square term; derive the continuum forcing INDEPENDENTLY in
  SymPy and Mathematica; require order ~2 + term-by-term nonzero diagnostics + target-token scan. **Does NOT need
  the 3 pre-1b items** (MMS uses a manufactured `S_j`, Dirichlet both ends).
- **1b** — closed residual + Newton: add `R0` to the torch residual + the radial-reduced return source; adapt
  pack/unpack, initial guess, resampling, JVP check, colored sparse preconditioner (layout `5*cells+nw+1`,
  `preconditioners.py:195`). **Needs the 3 pre-1b items resolved first.**
- **1c** — self-consistent-balance validation: recompute wall LHS/source from the converged `{ψ,A,R0}`, report
  component L∞ norms; reuse step-4 convergence + step-6 conservation style, R0-aware raw-field observables, no
  extraction map.

**Risks / stop conditions:** stop if `R0` goes nonpositive; if the exit BC must be changed to force convergence;
if the source sign / radial reduction is unresolved; if Newton converges only after hidden clamps. The
near-`D0→0` regime is OUT OF SCOPE for this placeholder build (stay away — generic placeholder `S_Σ` away from
the pole). CPU sparse LU caps ~30k DOF class (maybe lower with the denser `R0` column coupling).

## Pre-1b items to resolve (Claude+Codex math determination — BEFORE the closed-solve chunk)

1. **`δρ` background subtraction** in the static wall RHS: `ρ`, `ρ−ρ_ref`, or another declared subtraction?
   (Note: the static self-consistent balance is the FULL nonlinear `δS_total/δR=0`, so the RHS is the full static
   confinement force, not a linearized perturbation — confirm the exact form + reference.)
2. **Radial reduction convention** of `−k1 δρ` to the 1D wall grid: `Σ_i 4π r_i² dr (−k1 δρ)` in the flat-`dw`
   wall convention?
3. **Open exit BC**: natural Neumann / zero traction, or a specific Robin impedance (with what parameters)?

## Pre-1b items — RESOLVED (Claude+Codex consult, agreed; `_scratch/pathA_pre1b_items_consult.{md,log}`)

1. **Static RHS source = `S = −V_conf,Σ(r−R0(w))·ρ0 = −k1·ρ0`** — the FULL `ρ0` (the static balance is the full
   nonlinear `δS_total/δR=0`; `−k1 δρ` of D1 is its linearization). Sign: `L_int=−V_conf·ρ`, `∂_R L_int =
   +V_conf,Σ ρ`, wall source `= −δL/δR = −V_conf,Σ ρ` (matches D1's `−k1`). **Absolute `ρ0`, NO reference
   subtraction** (`ρ_ref=0`; a counterpressure/`ρ−ρ_ref` would be a new declared counterterm → GATE-A, not
   default). Same derivative convention as `p2_tangent.py:165` (`dV/dR0 = −4 V_radial r⁴/R0⁵`).
2. **Radial reduction = full spherical `4πr²` measure → per-`dw` wall source:**
   `S_j = Σ_i ΔV_i^r·[−V_conf,Σ(r_i−R0(w_j)) ρ0_ij]`, `ΔV_i^r = (4π/3)(r_{i+½}³−r_{i−½}³)` (the solver's shell
   volumes, `grid.py:80`), landed on `grid.w_centers` in the flat-`dw` wall convention (`operators.py:485`).
   The localized `V_wall` derivative concentrates the contribution near `r≈R0(w)`.
3. **Exit BC = natural open-end zero-traction** `T_{w,Σ}(R0(L),L) R0'(L)=0` (`Y_L=0` Robin limit), `R0(L)>0`,
   mouth `R0(0)=a`. Consistent with the open finite conduit (`compact:310-325`, no hard cap; Neumann = `Y_L→0`).
   A nonzero `Y_L` impedance would be a GATE-A `free_choice` DOF (deferred).

**Still GATE-A (freeze-sheet) choices, NOT defaults:** the actual `S_Σ` constitutive family; any nonzero exit
impedance `Y_L`; any background counterpressure/reference subtraction.

## Status
1a building (directive `directives/pathA_02_chunk1a_*`). Pre-1b items RESOLVED → 1b unblocked spec-wise (still
needs 1a done+reviewed first, since 1b extends the 1a operator into the closed Newton). Build is target-blind;
calibration run remains USER-freeze-gated (pre-reg §M Stage 2).
